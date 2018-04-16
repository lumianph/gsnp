#ifndef __DICT_COMPRESSION_H__
#define __DICT_COMPRESSION_H__

#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <gsnp_util.h>

//#include <thrust/unique.h>

using namespace std;

typedef uint32 code_t;

uint64 getNumBinWord(const uint64 numElement);
uint64 encode(code_t* code, const char* seq, const uint64 numElement);
void decode(char* seq, const uint64 numElement, const code_t* code, const uint64 numCodeWord);
void getNRegion(vector< pair<int, int> >* nRegion, const bool* nflag,
                const char* base1, const char* base2, const uint64 call_length);


void dictComprStart(const uint64 numElement, const uint64 dataSize);
void dictComprEnd();

void cuda_dictCompression(word* outData, uint64* numWordOutData,
                        char* dict, uint64* numElementDict,
                        char* inData, const uint64 numElementInData);
void cuda_dictCompression(word* outData, uint64* numWordOutData,
                        int* dict, uint64* numElementDict,
                        int* inData, const uint64 numElementInData);
void cuda_dictCompression(word* outData, uint64* numWordOutData,
                        uint32* dict, uint64* numElementDict,
                        uint32* inData, const uint64 numElementInData);
void cuda_dictCompression(word* outData, uint64* numWordOutData,
                        double* dict, uint64* numElementDict,
                        double* inData, const uint64 numElementInData);


#define MAX_DICT_SIZE (32*8192)

const static word maskTable[] = {0x0, 0x1, 0x3, 0x7, 0xf, 0x1f, 0x3f, 0x7f, 0xff, 
								 0x1ff, 0x3ff, 0x7ff, 0xfff, 0x1fff, 0x3fff, 0x7fff, 0xffff, 
								 0x1ffff, 0x3ffff, 0x7ffff, 0xfffff, 0x1fffff, 0x3fffff, 0x7fffff, 0xffffff, 
								 0x1ffffff, 0x3ffffff, 0x7ffffff, 0xfffffff, 0x1fffffff, 0x3fffffff, 0x7fffffff, 0xffffffff};

inline int getNumBitPerElement(const uint64 numElementDict) {
	int numBit = (int)ceil(log((double)(numElementDict))/log((double)2.f));
	if(numBit == 0) {	//cannot be zero
		numBit = 1;
	}

	return numBit;
}

inline int getNumElementPerWord(const uint64 numElementDict) {
	return (int)floor((sizeof(word)*8)/(double)getNumBitPerElement(numElementDict));
}


template<class DataType>
void dictDecompression(DataType* outData, const uint64 numElementOutData, 
					   const word* inData, const uint64 numWordInData, 
					   const DataType* dict, const uint64 numElementDict)
{
	const int numBitPerElement = getNumBitPerElement(numElementDict);
	const int numElementPerWord = getNumElementPerWord(numElementDict);

	uint64 elementId = 0;
	uint64 dictIdx = 0;
	uint64 numElementInWord = 0;
	const word mask = maskTable[numBitPerElement];

	for(uint64 wordId = 0; wordId < numWordInData; wordId++) {
		numElementInWord = (wordId == (numWordInData - 1)) ? 
					(numElementOutData - numElementPerWord*wordId) : (numElementPerWord);

		for(uint64 elementIdInWord = 0; elementIdInWord < numElementInWord; elementIdInWord++) {
			elementId = wordId*numElementPerWord + elementIdInWord;
			dictIdx = (inData[wordId]>>(numBitPerElement*(numElementInWord - elementIdInWord - 1)))&mask;
			assert(dictIdx < numElementDict);

			outData[elementId] = dict[dictIdx];
		}
	}
}



/**
 * RLE + dictionary compression
 * @param data input and output, the data is rewritten as unique elements
 * @param runLength output, the length is AT LEAST numElement
 * */
uint32 RLECompression(int* data, const uint64 numElement, uint32* runLength);
uint32 RLECompression(double* data, const uint64 numElement, uint32* runLength);
uint32 RLECompression(char* data, const uint64 numElement, int* runLength);


template<class DataType>
void RLEDecompression(DataType* data, const DataType* runValue, 
		      const uint32* runLength, const uint64 numRun) {
	uint64 idx = 0;
	for(uint64 runId = 0; runId < numRun; runId++) {
		for(uint32 offset = 0; offset < runLength[runId]; offset++) {
			data[idx] = runValue[runId];
			idx++;
		}
	}
}


/**
 * dictionary compression
 * */
template<class DataType>
void dictCompression(word* outData, uint64* numWordOutData /*in and out*/, 
			DataType* dict, uint64* numElementDict /*in and out*/, 
			const DataType* inData, const uint64 numElementInData) 
{

	/*1. find the size of dictionary*/
        if(*numElementDict < numElementInData) {
                printf("numElementDict = %d, numElementInData = %d\n",
                        *numElementDict, numElementInData);
                ERROR_EXIT("*numElementDict < numElementInData");
        }
	memcpy(dict, inData, sizeof(DataType)*numElementInData);
	sort(dict, dict + numElementInData);
	DataType* dictEnd = unique(dict, dict + numElementInData);
	const uint64 dictSize = dictEnd - dict;

	if(*numElementDict < dictSize) {
		ERROR_EXIT("numElementDict < dictSize\n");
	}
	*numElementDict = dictSize;

	/*2. encode the data using bits according to the dictionary*/
	const int numBitPerElement = getNumBitPerElement(*numElementDict);
	assert(numBitPerElement <= 8*sizeof(word));
	assert(pow(2.0, (int)numBitPerElement) >= (*numElementDict));
	const int numElementPerWord = getNumElementPerWord(*numElementDict);
	assert(numElementPerWord > 0);
	const int numWordCompressed = (int)ceil(numElementInData/(double)numElementPerWord);
	if(*numWordOutData < numWordCompressed) {
		ERROR_EXIT("numWordOutData < numWordCompressed\n");
	}
	*numWordOutData = numWordCompressed;

	
	uint64 wordIdx = 0;
	DataType* p = NULL;
	long idx = 0;
	memset(outData, 0, sizeof(word)*(*numWordOutData));
	for(uint64 dataIdx = 0; dataIdx < numElementInData; dataIdx++) {
		//find index in the dictionary
		wordIdx = dataIdx/numElementPerWord;
		p = find(dict, dictEnd, inData[dataIdx]);
		idx = p - dict;
		assert(idx >= 0 && idx < *numElementDict);

		//encoding
		outData[wordIdx] <<= (numBitPerElement);
		outData[wordIdx] |= idx;
	}
}



template<class DataType>
struct RLEDictData {
	word* cRunValue;		//compressed run value
	word* cRunLength;		//compressed run length
	DataType* runValueDict;		//the dictionary of cRunValue
	uint32* runLengthDict;		//the dictionary of cRunLength
	uint32* runLengthBuf;

	//for RLE compressed data
	uint64 numRun;

	//for dict compressed data
	uint64 numWordRunValue;		
	uint64 numWordRunLength;	
	uint64 numEleRunValueDict;
	uint64 numEleRunLengthDict;
	
	uint64 maxNumElement;
	uint64 numElement;

	RLEDictData(const uint64 in_maxNumElement, bool alloc = true) {
		maxNumElement = in_maxNumElement;
		numElement = 0;

		if(alloc) {
			cRunValue = (word*)malloc(sizeof(word)*maxNumElement);
			cRunLength = (word*)malloc(sizeof(word)*maxNumElement);
			runValueDict = (DataType*)malloc(sizeof(DataType)*maxNumElement);
			runLengthDict = (uint32*)malloc(sizeof(DataType)*maxNumElement);
		}

		runLengthBuf = (uint32*)malloc(sizeof(uint32)*maxNumElement);
		numWordRunValue = 0;
		numWordRunLength = 0;
		numEleRunValueDict = 0;
		numEleRunLengthDict = 0;
	}

	

	void close() {
		free(cRunValue);
		free(cRunLength);
		free(runValueDict);
		free(runLengthDict);
		free(runLengthBuf);

		numRun = 0;
                numWordRunValue = 0;
                numWordRunLength = 0;
                numEleRunValueDict = 0;
                numEleRunLengthDict = 0;
	}

	uint64 getDataSize() {
		return sizeof(word)*numWordRunValue + sizeof(word)*numWordRunLength + 
			sizeof(DataType)*numEleRunValueDict + sizeof(uint32)*numEleRunLengthDict;
	}

	/**
	 * RLE compression + dictionary compression
	 * @param inData the content will be rewritten after return
	 * @param numElement
	 * @param buf at least size sizeof(uint32)*numElement
	 * */
	void compress(DataType* data, const uint64 in_numElement) {


		//INIT_TIMER;
	

		/*0. parameters checking*/
		assert(in_numElement <= maxNumElement);
		numElement = in_numElement;

		/*1. run length encoding, the output is in data and runLength*/
		//cRunLength can hold all since the type size of word and uint32 are the same
		//START_TIMER;
		numRun = RLECompression(data, numElement, runLengthBuf);
		//END_TIMER;
		//PRINT_TIMER_SEC("RLE compression");
		
		/*2. compress run value(stored in the data) using dict*/
		numWordRunValue = numRun;
		numEleRunValueDict = MAX_DICT_SIZE;
		//START_TIMER;
		cuda_dictCompression(cRunValue, &numWordRunValue, runValueDict, &numEleRunValueDict,
                                data, numRun);
		//END_TIMER;
		//PRINT_TIMER_SEC("dict. compression 1");
		

		/*3. compress run length(stored in the buf) using dict*/
		numWordRunLength = numRun;
		numEleRunLengthDict = MAX_DICT_SIZE;
		//START_TIMER;
		cuda_dictCompression(cRunLength, &numWordRunLength, runLengthDict, &numEleRunLengthDict,
				runLengthBuf, numRun);
		//END_TIMER;
		//PRINT_TIMER_SEC("dict. compression 2");
	}

	/**
 	* @param data output, should be allocated appropriately
 	* @param numElement
 	* */
	void decompress(DataType* data, byte* buf, const uint64 bufSize) {
		/*0. parameters checking*/
		if(bufSize < (sizeof(DataType)*numRun + sizeof(uint32)*numRun)) {
			ERROR_EXIT("bufSize is too small to hold decompressed run value and run length.");
		}
		DataType* runValue = (DataType*)buf;
		uint32* runLength = (uint32*)(buf + sizeof(DataType)*numRun);

		/* 1. decompress run value using dict*/
		dictDecompression(runValue, numRun, cRunValue, numWordRunValue,
                                  runValueDict, numEleRunValueDict);

		/* 2. decompress run length using dict*/
		dictDecompression(runLength, numRun, cRunLength, numWordRunLength, 
				  runLengthDict, numEleRunLengthDict);

	
		/*3. decompress data using RLE*/
		RLEDecompression(data, runValue, runLength, numRun);
	}
};

#endif /* __DICT_COMPRESSION_H__ */
