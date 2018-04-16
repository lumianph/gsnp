#ifndef __COMPRESS_KERNEL_CU__
#define __COMPRESS_KERNEL_CU__

#include "compress.h"
#include <gsnp_util.h>
//#include <thrust/device_ptr>
#include <thrust/unique.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>
#include "cuda_header.h"

#define PRINT_COMPRESS_INFO

#define _1MB (1024.0*1024.0)
#define _1KB (1024.0)


/* data type definition */
#define NUM_BASE_PER_BYTE (4) //2 bits per base

//typedef uint32 code_t;
#define CODE_MAX (UINT32_MAX)


inline bool isGoodBase(const char base) {

        if((base == 'A') || (base == 'a') ||
           (base == 'C') || (base == 'c') ||
           (base == 'G') || (base == 'g') ||
           (base == 'T') || (base == 't')) {
                return true;
        } else {
                return false;
        }
}


inline uint32 base2int(const char base) {
        if(base == 'A' || base == 'a') {
                return 0;
        }

        if(base == 'C' || base == 'c') {
                return 1;
        }

        if(base == 'G' || base == 'g') {
                return 2;
        }

        if(base == 'T' || base == 't') {
                return 3;
        }

        ERROR_EXIT("unknown base when converting to an integer");
        return 0xffffffff;
}


static char baseTable[4] = {'A', 'C', 'G', 'T'};

inline char int2base(const int code) {

	/*
	if(code == 0) {
		return 'A';
	}
	if(code == 1) {
		return 'C';
	}
	if(code == 2) {
		return 'G';
	}
	if(code == 3) {
		return 'T';
	}

	ERROR_EXIT("unknown code");
	return 0xff;
	*/

	return baseTable[code];
}

uint64 getNumBinWord(const uint64 numElement) {
        const uint32 numBasePerWord = NUM_BASE_PER_BYTE*sizeof(code_t);
        const uint64 numWord = ceil(numElement/((double)numBasePerWord));

        return numWord;
}


/**
 * encode a sequence (A, C, G, T a, c, g, t only) to a binary code
 * A, a: 00
 * C, c: 01
 * G, g: 10
 * T, t: 11
 * @return the number of words allocated for the code array
 * @param code the encoded words
 * @param seq the original base array, size: numElement
 * @param numElement #base in the seq
 * */
uint64 encode(code_t* code, const char* seq, const uint64 numElement) {
        const uint32 numBasePerWord = NUM_BASE_PER_BYTE*sizeof(code_t);
        const uint64 numWord = getNumBinWord(numElement);
        memset(code, 0, sizeof(code_t)*numWord);

        for(uint64 wordId = 0; wordId < numWord; wordId++) {
                const uint32 numBaseThisWord = (wordId == (numWord - 1)) ?
                                (numElement - wordId*numBasePerWord) :
                                (numBasePerWord);
                assert(numBaseThisWord <= numBasePerWord);

                uint32 word = 0;
                for(uint32 baseId = 0; baseId < numBaseThisWord; baseId++) {
                        const char base = seq[wordId*numBasePerWord + baseId];

                        if(isGoodBase(base)) {
                                word = word<<2;
                                word |= base2int(base);
                        } else {
                                ERROR_EXIT("unknown base when encoding.");
                        }
                }
                code[wordId] = word;
        }

        return numWord;
}


void decode(char* seq, const uint64 numElement, 
	    const code_t* code, const uint64 numCodeWord) {

	const uint32 numBasePerWord = NUM_BASE_PER_BYTE*sizeof(code_t);
	const code_t mask = 0x3;

	//all complete code words
	for(uint64 wordId = 0; wordId < numCodeWord - 1; wordId++) {
		for(uint32 baseId = 0; baseId < numBasePerWord; baseId++) {
			code_t c = (code[wordId]>>(2*numBasePerWord - (baseId + 1)*2))&mask;
			seq[wordId*numBasePerWord + baseId] = int2base(c);
		}
	}

	const uint32 numBaseLeft = (numElement - (numCodeWord - 1)*numBasePerWord);
	assert(numBaseLeft <= numBasePerWord);
	for(uint32 baseId = 0; baseId < numBaseLeft; baseId++) {
		code_t c = (code[numCodeWord - 1]>>(2*(numBaseLeft - baseId - 1)))&mask;
		seq[(numCodeWord - 1)*numBasePerWord + baseId] = int2base(c);
	}
}

uint64 encode(code_t** code, const char* seq, const uint64 numElement) {
        const uint64 numWord = getNumBinWord(numElement);
        *code = (code_t*)malloc(sizeof(code_t)*numWord);

        return encode(*code, seq, numElement);
}



static byte* d_inData = NULL;
static word* d_outData = NULL;
static int* d_tmpBuf = NULL;
static byte* d_dict = NULL;

void dictComprStart(const uint64 numElement, const uint64 dataSize) {
	GPUMALLOC((void**)&d_inData, dataSize*numElement);
	GPUMALLOC((void**)&d_outData, sizeof(word)*numElement);
	GPUMALLOC((void**)&d_dict, dataSize*numElement);
	GPUMALLOC((void**)&d_tmpBuf, sizeof(int)*numElement);
}


void dictComprEnd() {
	GPUFREE(d_inData);
	GPUFREE(d_outData);
	GPUFREE(d_dict);
	GPUFREE(d_tmpBuf);
}


__global__
void dictEncoding_kernel(word* outData, const int* idx, const int numWord,
                        const int numElementPerWord, const int numBitPerElement) {

        const int numTotalThread = NUM_TOTAL_THREAD;
        const int globalThreadOffset = GLOBAL_THREAD_OFFSET;
        int elementId = 0;

        for(int wordId = globalThreadOffset; wordId < numWord; wordId += numTotalThread) {

                word out = 0;
                for(int i = 0; i < numElementPerWord; i++) {
                        elementId = wordId*numElementPerWord + i;

                        out <<= numBitPerElement;
                        out |= idx[elementId];
                }
                outData[wordId] = out;
        }
}

//TODO: use shared memory for the better performance
void cudaDictEncoding(word* outData, int* idx, const int numElementInData,
                        const int numBitPerElement, const int numElementPerWord, const int numWordCompressed) {
        int numBlock = 2400;
        int numThread = 256;

        //except the last word
        dictEncoding_kernel<<<numBlock, numThread>>>(outData, idx,
                                        numWordCompressed - 1, numElementPerWord, numBitPerElement);
        cutilCheckMsg("dictEncoding_kernel 1");

        //the last one
        dictEncoding_kernel<<<1, 1>>>(outData + numWordCompressed - 1, idx + numElementPerWord*(numWordCompressed - 1),
                                1, numElementInData - numElementPerWord*(numWordCompressed - 1), numBitPerElement);
        cutilCheckMsg("dictEncoding_kernel 2");
}


template<class DataType>
void dictCompressionKernel(word* outData, uint64* numWordOutData /*in and out*/,
                        thrust::device_ptr<DataType> dict, uint64* numElementDict /*in and out*/,
                        thrust::device_ptr<DataType> inData, const uint64 numElementInData,
                        thrust::device_ptr<int> buf /*with the same #element as input*/) {

	//PRINT_FUNC_NAME;
	//unsigned int timer = 0;


        /*build the dictionary*/
	//startTimer(&timer);
        if(*numElementDict < numElementInData) {
		printf("numElementDict = %d, numElementInData = %d\n", 
			*numElementDict, numElementInData);

                ERROR_EXIT("*numElementDict < numElementInData");
        }
        thrust::copy(inData, inData + numElementInData, dict);
        thrust::sort(dict, dict + numElementInData);
        thrust::device_ptr<DataType> dictEnd = thrust::unique(dict, dict + numElementInData);
        *numElementDict = dictEnd - dict;
	//endTimer(&timer, "building dictionary");

        /* search the code */
	//timer = 0;
	//startTimer(&timer);
        const int numBitPerElement = getNumBitPerElement(*numElementDict);
        const int numElementPerWord = getNumElementPerWord(*numElementDict);
        assert(numElementPerWord > 0);
        const int numWordCompressed = (int)ceil(numElementInData/(double)numElementPerWord);
        if(*numWordOutData < numWordCompressed) {
                ERROR_EXIT("numWordOutData < numWordCompressed\n");
        }
        *numWordOutData = numWordCompressed;
        thrust::lower_bound(dict, dictEnd, inData, inData + numElementInData, buf);
	//endTimer(&timer, "search the index");

        /*encode to the binary representation*/
	//timer = 0;
	//startTimer(&timer);
        cudaDictEncoding(outData, thrust::raw_pointer_cast(buf),
                        numElementInData, numBitPerElement, numElementPerWord, numWordCompressed);
	//endTimer(&timer, "index compaction");
}

void cuda_dictCompression(word* outData, uint64* numWordOutData, 
                        uint32* dict, uint64* numElementDict,
                        uint32* inData, const uint64 numElementInData) {

	thrust::device_ptr<uint32> dictPtr((uint32*)d_dict);
	thrust::device_ptr<uint32> inDataPtr((uint32*)d_inData);
	thrust::copy(inData, inData + numElementInData, inDataPtr);
	thrust::device_ptr<int> bufPtr(d_tmpBuf);

	dictCompressionKernel(d_outData, numWordOutData, dictPtr, numElementDict, inDataPtr, numElementInData, bufPtr);

	FROMGPU(outData, d_outData, sizeof(word)*(*numWordOutData));
	FROMGPU(dict, (uint32*)d_dict, sizeof(uint32)*(*numElementDict));

	cutilSafeCall(cudaThreadSynchronize());
}

void cuda_dictCompression(word* outData, uint64* numWordOutData,
                        double* dict, uint64* numElementDict,
                        double* inData, const uint64 numElementInData) {
        thrust::device_ptr<double> dictPtr((double*)d_dict);
        thrust::device_ptr<double> inDataPtr((double*)d_inData);
        thrust::copy(inData, inData + numElementInData, inDataPtr);
        thrust::device_ptr<int> bufPtr(d_tmpBuf);

        dictCompressionKernel(d_outData, numWordOutData, dictPtr, numElementDict, inDataPtr, numElementInData, bufPtr);

        FROMGPU(outData, d_outData, sizeof(word)*(*numWordOutData));
        FROMGPU(dict, (double*)d_dict, sizeof(double)*(*numElementDict));

	cutilSafeCall(cudaThreadSynchronize());
}


void cuda_dictCompression(word* outData, uint64* numWordOutData,
                        char* dict, uint64* numElementDict,
                        char* inData, const uint64 numElementInData) {
        thrust::device_ptr<char> dictPtr((char*)d_dict);
        thrust::device_ptr<char> inDataPtr((char*)d_inData);
        thrust::copy(inData, inData + numElementInData, inDataPtr);
        thrust::device_ptr<int> bufPtr(d_tmpBuf);

        dictCompressionKernel(d_outData, numWordOutData, dictPtr, numElementDict, inDataPtr, numElementInData, bufPtr);

        FROMGPU(outData, d_outData, sizeof(word)*(*numWordOutData));
        FROMGPU(dict, (char*)d_dict, sizeof(char)*(*numElementDict));

	cutilSafeCall(cudaThreadSynchronize());
}

void cuda_dictCompression(word* outData, uint64* numWordOutData,
                        int* dict, uint64* numElementDict,
                        int* inData, const uint64 numElementInData) {
        thrust::device_ptr<int> dictPtr((int*)d_dict);
        thrust::device_ptr<int> inDataPtr((int*)d_inData);
        thrust::copy(inData, inData + numElementInData, inDataPtr);
        thrust::device_ptr<int> bufPtr(d_tmpBuf);

        dictCompressionKernel(d_outData, numWordOutData, dictPtr, numElementDict, inDataPtr, numElementInData, bufPtr);
                        
        FROMGPU(outData, d_outData, sizeof(word)*(*numWordOutData));
        FROMGPU(dict, (int*)d_dict, sizeof(int)*(*numElementDict));

	cutilSafeCall(cudaThreadSynchronize());
}

uint32 RLECompression(int* data, const uint64 numElement, uint32* runLength) {

        /*0. check parameters*/
        assert(numElement <= (uint64)0xffffffff); //only support 32 bit in this version

        /*1.apply run-length encoding*/
        for(uint32 i = 0; i < numElement; i++) {
                runLength[i] = i;
        }
        thrust::pair<int*, uint32*> end = thrust::unique_by_key(data, data + numElement, runLength);
        const uint32 numRun = end.first - data;
        for(uint32 runId = 0; runId < numRun - 1; runId++) {
                runLength[runId] = runLength[runId + 1] - runLength[runId];
        }
        runLength[numRun - 1] = numElement - runLength[numRun - 1]; //the last one

        return numRun;
}


uint32 RLECompression(double* data, const uint64 numElement, uint32* runLength) {

        /*0. check parameters*/
        assert(numElement <= (uint64)0xffffffff); //only support 32 bit in this version

        /*1.apply run-length encoding*/
        for(uint32 i = 0; i < numElement; i++) {
                runLength[i] = i;
        }
        thrust::pair<double*, uint32*> end = thrust::unique_by_key(data, data + numElement, runLength);
        const uint32 numRun = end.first - data;
        for(uint32 runId = 0; runId < numRun - 1; runId++) {
                runLength[runId] = runLength[runId + 1] - runLength[runId];
        }
        runLength[numRun - 1] = numElement - runLength[numRun - 1]; //the last one

        return numRun;
}


uint32 RLECompression(char* data, const uint64 numElement, int* runLength) {

        /*0. check parameters*/
        assert(numElement <= (uint64)0xffffffff); //only support 32 bit in this version

        /*1.apply run-length encoding*/
        for(uint32 i = 0; i < numElement; i++) {
                runLength[i] = i;
        }
        thrust::pair<char*, int*> end = thrust::unique_by_key(data, data + numElement, runLength);
        const uint32 numRun = end.first - data;
        for(uint32 runId = 0; runId < numRun - 1; runId++) {
                runLength[runId] = runLength[runId + 1] - runLength[runId];
        }
        runLength[numRun - 1] = numElement - runLength[numRun - 1]; //the last one

        return numRun;
}

#endif /* __COMPRESS_KERNEL_CU__ */


