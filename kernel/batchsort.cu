#ifndef __BATCH_SORT_CU__
#define __BATCH_SORT_CU__


#include <stdio.h>
#include <stdlib.h>
#include <algorithm>


//#include <cutil.h>
#include <gsnp_util.h>
#include "cuda_header.h"
#include "bitonic_kernel.cu"

#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <thrust/fill.h>
#include <thrust/reduce.h>
#include <thrust/unique.h>

using namespace std;


#define POW2_TABLE_SIZE (14)
uint32 pow2Table[POW2_TABLE_SIZE] = {
        2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384
};

uint32 getReduceSize(const uint32 n) {
        bool isFind = false;
        uint32 m = 0;

        for(uint32 i = 0; i < POW2_TABLE_SIZE; i++) {
                if(pow2Table[i] >= n) {
                        isFind = true;
                        m = pow2Table[i];
                        break;
                }
        }

        if(!isFind) {
                ERROR_EXIT("input number is too large or small!");
        }

        return m;
}


uint32 factorRadix2(uint32 *log2L, uint32 L){
    if(!L){
        *log2L = 0;
        return 0;
    }else{
        for(*log2L = 0; (L & 1) == 0; L >>= 1, *log2L++);
        return L;
    }
}



void bitonicSortSP(uint32* d_data, const int numArray, const int arrayLength) {

	/*parameter checking*/	
        uint32 log2L;
        uint32 factorizationRemainder = factorRadix2(&log2L, arrayLength);
        if( factorizationRemainder != 1 ) {
		ERROR_EXIT("only support the array length is power of 2!");
	}
	const int numBlock = 2400;//multiple pass in the kernel when numArray > numBlock
	const int numThread = arrayLength;
	const int sharedMemSize = sizeof(uint32)*arrayLength;
	assert(numBlock <= MAX_NUM_BLOCK);
	assert(numThread <= MAX_NUM_THREAD);
	assert(sharedMemSize <= MAX_SHARED_MEM_SIZE);

	/*call the kernel*/
	bitonicSort_kernel<<<numBlock, numThread, sharedMemSize>>>(d_data, numArray, arrayLength);
	cutilCheckMsg("bitonicSort_kernel");

	cutilSafeCall(cudaThreadSynchronize());
}

void batchSortSP(uint32* d_data, const int arrayPitch, 
		const int numArray, const uint32* d_arrayLength, const int maxArrayLength) {

	assert(maxArrayLength <= arrayPitch);
	assert(maxArrayLength <= MAX_NUM_THREAD);
        uint32 log2L;
        uint32 factorizationRemainder = factorRadix2(&log2L, maxArrayLength);
        if( factorizationRemainder != 1 ) {
                ERROR_EXIT("only support the max. array length is power of 2!");
        }

        const int numBlock = 2400;

        //length (1, 32]
        const int numThread = maxArrayLength;
        const int sharedMemSize = sizeof(uint32)*numThread;
	assert(sharedMemSize <= MAX_SHARED_MEM_SIZE);
        batchSort_kernel<<<numBlock, numThread, sharedMemSize>>>
                (d_data, arrayPitch, numArray, d_arrayLength, 1, maxArrayLength);
        cutilCheckMsg("batchSort_kernel");

	cutilSafeCall(cudaThreadSynchronize());
}


void batchSortMP(uint32* d_data, const int arrayPitch, const int numArray, const uint32* d_arrayLength) {

	int numBlock, numThread, sharedMemSize;

	/*parameter checking*/
        uint32 log2L;
        uint32 factorizationRemainder = factorRadix2(&log2L, arrayPitch);
        if( factorizationRemainder != 1 ) {
                ERROR_EXIT("only support the array length is power of 2!");
        }

	numBlock = 2400;


	//length (1, 32]
	numThread = 32;
	sharedMemSize = sizeof(uint32)*numThread;
	batchSort_kernel<<<numBlock, numThread, sharedMemSize>>>
		(d_data, arrayPitch, numArray, d_arrayLength, 1, 32);
	cutilCheckMsg("batchSort_kernel 1");


	//length (32, 64]
        numThread = 64;
        sharedMemSize = sizeof(uint32)*numThread;
        batchSort_kernel<<<numBlock, numThread, sharedMemSize>>>
                (d_data, arrayPitch, numArray, d_arrayLength, 32, 64);
        cutilCheckMsg("batchSort_kernel 2");


	//length > 64
        numThread = arrayPitch;
        sharedMemSize = sizeof(uint32)*numThread;
	assert(numThread <= MAX_NUM_THREAD);
	assert(sharedMemSize <= MAX_SHARED_MEM_SIZE);
        batchSort_kernel<<<numBlock, numThread, sharedMemSize>>>
                (d_data, arrayPitch, numArray, d_arrayLength, 64, arrayPitch);
        cutilCheckMsg("batchSort_kernel 3");

	cutilSafeCall(cudaThreadSynchronize());
}



/**
* multiple passes
* pass 0: [0, 1]
* pass 1: (1, 8]
* pass 2: (8, 16]
* pass 3: (16, 32]
* pass 4: (32, 64]
* pass 5: > 64
*/

static int maxNumArray = 0;

int* d_arrayId = NULL;
int* d_passId = NULL;
int* d_offset = NULL;
int* h_passId = NULL;
int* h_offset = NULL;
const static int maxNumPass = 6;
static int upTable[6] = {1, 8, 16, 32, 64, 0};

void batchSortBegin(const int in_maxNumArray) {
	maxNumArray = in_maxNumArray;

	GPUMALLOC((void**)&d_arrayId, sizeof(int)*maxNumArray);
	GPUMALLOC((void**)&d_passId, sizeof(int)*maxNumArray);
	GPUMALLOC((void**)&d_offset, sizeof(int)*maxNumArray);

	h_passId = (int*)malloc(sizeof(int)*maxNumPass);
	h_offset = (int*)malloc(sizeof(int)*maxNumPass);	
}

void batchSortFinalize() {
	GPUFREE(d_arrayId);
	GPUFREE(d_passId);
	GPUFREE(d_offset);
	maxNumArray = 0;

	free(h_passId);
	free(h_offset);
}


__global__
void assignPassId_kernel(int* d_offset, int* d_arrayId, int* d_passId, const uint32* d_arrayLength, const int numArray) {
	const int globalThreadOffset = GLOBAL_THREAD_OFFSET;
	const int numTotalThread = NUM_TOTAL_THREAD;
	int len, passId;	

	for(int arrayId = globalThreadOffset; arrayId < numArray; arrayId += numTotalThread) {
		len = d_arrayLength[arrayId];

		if(len <= 1) {
			passId = 0;
		} else if (len <= 8) {
			passId = 1;
		} else if (len <= 16) {
			passId = 2;
		} else if (len <= 32) {
			passId = 3;
		} else if (len <= 64) {
			passId = 4;
		} else {
			passId = 5;
		}

		d_offset[arrayId] = arrayId;
		d_arrayId[arrayId] = arrayId;
		d_passId[arrayId] = passId;
	}
}


int preparePass(int* d_passId, int* d_arrayId, int* d_offset, const uint32* d_arrayLength, const int numArray) {
	

	//assigne pass id for each array
	assignPassId_kernel<<<1200, 256>>>(d_offset, d_arrayId, d_passId, d_arrayLength, numArray);
	cutilCheckMsg("assignPassId_kernel");

	//sort by key according to passId, the value is stored in d_arrayId
	thrust::device_ptr<int> dptr_arrayId(d_arrayId);
	thrust::device_ptr<int> dptr_passId(d_passId);
	thrust::sort_by_key(dptr_passId, dptr_passId + numArray, dptr_arrayId);

	//unique_by_key, calculate the offset of each pass on the d_arrayId 
	thrust::device_ptr<int> dptr_offset(d_offset);
	thrust::pair< thrust::device_ptr<int>, thrust::device_ptr<int> > new_end = 
		thrust::unique_by_key(dptr_passId, dptr_passId + numArray, dptr_offset);


	return new_end.first - dptr_passId;
}

void batchSortMP2(uint32* d_data, const int arrayPitch, const int numArray, const uint32* d_arrayLength) {

	/*0. parameters checking*/
        uint32 log2L;
        uint32 factorizationRemainder = factorRadix2(&log2L, arrayPitch);
        if( factorizationRemainder != 1 ) {
                ERROR_EXIT("only support the array pitch is power of 2!");
        }
	assert(numArray <= maxNumArray);
	
	upTable[5] = arrayPitch;	//assigne the max. upper bound

	/*1. get the ids for different passes*/
	int numPass = preparePass(d_passId, d_arrayId, d_offset, d_arrayLength, numArray);
	assert(numPass <= maxNumPass);
	FROMGPU(h_passId, d_passId, sizeof(int)*numPass);
	FROMGPU(h_offset, d_offset, sizeof(int)*numPass);
	cutilSafeCall(cudaThreadSynchronize());

	/*2. sorting for each pass*/
	const int numBlock = 1200;
	const int numThread = 256;
	const int sharedMemSize = sizeof(uint32)*numThread;

	//printf("numPass = %d\n", numPass);
	for(int pass = 0; pass < numPass; pass++) {
		const int begin = h_offset[pass];
		const int end = (pass == numPass - 1) ? (numArray) : (h_offset[pass + 1]);
		const int length = end - begin;
		const int passId = h_passId[pass];
		assert(passId < 6);
		if(passId == 0) {
			continue; //0 or 1, do not need to be sorted
		}
		const int up = upTable[passId];


		//printf("pass = %d, arrayPitch = %d, begin = %d, length = %d, up = %d, numArray = %d\n", 
		//	pass, arrayPitch, begin, length, up, numArray);

		batchSortMP2_kernel<<<numBlock, numThread, sharedMemSize>>>
			(d_data, arrayPitch, d_arrayId + begin, length, up, numArray);
		cutilCheckMsg("batchSortMP2_kernel");
	}

	cutilSafeCall(cudaThreadSynchronize());
}


void batchSortSP2(uint32* d_data, const int numArray, const int arrayLength, const int arrayPitch) {

        const int numThread = 256;
	const int numBlock = 1200;
        const int sharedMemSize = sizeof(uint32)*numThread;
        batchSort_kernel<<<numBlock, numThread, sharedMemSize>>>(d_data, numArray, arrayLength, arrayLength);
        cutilCheckMsg("batchSort_kernel");
        cutilSafeCall(cudaThreadSynchronize());
}

void test_batchSort(const int numArray, const int arrayPitch) {
	
	PRINT_FUNC_NAME;

	srand(time(NULL));

	batchSortBegin(numArray);


	const int numElement = numArray*arrayPitch;
	uint32* h_data = (uint32*)malloc(sizeof(uint32)*numElement);
	uint32* h_arrayLength = (uint32*)malloc(sizeof(uint32)*numArray);
	int table[6] = {0};
	for(int i = 0; i < numArray; i++) {
		h_arrayLength[i] = rand()%arrayPitch;

                const int len = h_arrayLength[i];

                if(len <= 1) {
                        table[0]++;
                } else if (len <= 8) {
			table[1]++;
                } else if (len <= 16) {
                        table[2]++;
                } else if (len <= 32) {
                        table[3]++;
                } else if (len <= 64) {
                        table[4]++;
                } else {
                        table[5]++;
                }
	}

	printf("array length distribution: \n");
	for(int i = 0; i < 6; i++) {
		printf("%d, %d\n", i, table[i]);
	}
	printLine();

	memset(h_data, 0xffff, sizeof(uint32)*numElement);
	for(int arrayId = 0; arrayId < numArray; arrayId++) {
		for(int i = 0; i < h_arrayLength[arrayId]; i++) {
			h_data[arrayId*arrayPitch + i] = rand();
		}
	}

	unsigned int timer = 0;

	uint32* d_data = NULL;
	GPUMALLOC((void**)&d_data, sizeof(uint32)*numElement);
	TOGPU(d_data, h_data, sizeof(uint32)*numElement);
	uint32* d_arrayLength = NULL;
	GPUMALLOC((void**)&d_arrayLength, sizeof(uint32)*numArray);
	TOGPU(d_arrayLength, h_arrayLength, sizeof(uint32)*numArray);

	//GPU batchSortMP	
	timer = 0;
	startTimer(&timer);
	batchSortMP2(d_data, arrayPitch, numArray, d_arrayLength);
	double gpu_sec = endTimer(&timer, "GPU batchSortMP")/1000.0f;
	uint32* gpu_out = (uint32*)malloc(sizeof(uint32)*numElement);
	FROMGPU(gpu_out, d_data, sizeof(uint32)*numElement);

	//GPU batchSortSP
/*	timer = 0;
	startTimer(&timer);
	bitonicSortSP(d_data, numArray, arrayPitch);
	endTimer(&timer, "GPU batchSortSP");
*/

	//CPU
	timer = 0;
	startTimer(&timer);
	for(int arrayId = 0; arrayId < numArray; arrayId++) {
		sort(h_data + arrayId*arrayPitch, h_data + arrayId*arrayPitch + (h_arrayLength[arrayId]));
	}
	endTimer(&timer, "CPU sort");

/*

	int offset = 0;
	for(int arrayId = 0; arrayId < numArray; arrayId++) {
		printf("arrayId = %d, arrayLength = %d\n", arrayId, h_arrayLength[arrayId]);
		for(int i = 0; i < arrayPitch; i++) {
			printf("%d: %u\t%u\n", offset, h_data[arrayId*arrayPitch + i], gpu_out[arrayId*arrayPitch + i]);
			offset++;
		}		

		printLine();
	}

*/


	//CHECK RESULTS	
	simpleCheck(h_data, gpu_out, numElement);

	//output
        printf("numArray = %d, arrayPitch = %d, numElement = %d\n", numArray, arrayPitch, numElement);
        printf("bandwidth = %.3f M elements per sec\n", (numElement/(1000000.0f))/gpu_sec);



	batchSortFinalize();
}


void test_bandwidth(const int numArray, const int arrayLength) {
	PRINT_FUNC_NAME;
	printf("numArray = %d, arrayLength = %d\n", numArray, arrayLength);

	const uint64 numElement = numArray*arrayLength;
	uint32* h_data = (uint32*)malloc(sizeof(uint32)*numElement);
	for(uint64 i = 0; i < numElement; i++) {
		h_data[i] = rand();
	}

	unsigned int timer = 0;

	//GPU bitonic sort
	uint32* d_data = NULL;
	GPUMALLOC((void**)&d_data, sizeof(uint32)*numElement);
	TOGPU(d_data, h_data, sizeof(uint32)*numElement);
	startTimer(&timer);
		bitonicSortSP(d_data, numArray, arrayLength);
		cutilSafeCall(cudaThreadSynchronize());	
	double gpu_sec = endTimer(&timer, "GPU bitonic parallel sort")/1000.0f;
	uint32* gpu_out = (uint32*)malloc(sizeof(uint32)*numElement);
	FROMGPU(gpu_out, d_data, sizeof(uint32)*numElement);

	//GPU sequential sort for each array
	uint32* d_data2 = NULL;
	GPUMALLOC((void**)&d_data2, sizeof(uint32)*numElement);
	TOGPU(d_data2, h_data, sizeof(uint32)*numElement);
	thrust::device_ptr<uint32> d_data2Ptr(d_data2);
        timer = 0;
        startTimer(&timer);
                for(int arrayId = 0; arrayId < numArray; arrayId++) {
                        thrust::sort(d_data2Ptr + arrayId*arrayLength, d_data2Ptr + (arrayId + 1)*arrayLength);
                }
                cutilSafeCall(cudaThreadSynchronize());
        double gpu_sec2 = endTimer(&timer, "GPU	thrust seqeuntial sort")/1000.0f;
	uint32* gpu_out2 = (uint32*)malloc(sizeof(uint32)*numElement);
	FROMGPU(gpu_out2, d_data2, sizeof(uint32)*numElement);

	//CPU
	uint32* gold_out = h_data;
	timer = 0;
	startTimer(&timer);
		for(int arrayId = 0; arrayId < numArray; arrayId++) {
			sort(gold_out + arrayId*arrayLength, gold_out + (arrayId + 1)*arrayLength);
		}
		cutilSafeCall(cudaThreadSynchronize());
	double cpu_sec = endTimer(&timer, "CPU")/1000.0f;

	//check results
	simpleCheck(gold_out, gpu_out, numElement, "check sorting result 1");
	simpleCheck(gold_out, gpu_out2, numElement, "check sorting result 2");
	
	//output
	printf("GPU bitonic sort bandwidht = %.3f M elements per sec\n", numArray*arrayLength/double(1024.0*1024.0)/gpu_sec);
	printf("GPU thrust sequential bandwidth: %.3f M elements per sec\n", numArray*arrayLength/double(1024.0*1024.0)/gpu_sec2);
	printf("CPU qsort bandwidth = %.3f M elements per sec\n", numArray*arrayLength/double(1024.0*1024.0)/cpu_sec);
}


void test_batchSort2(const int numArray, const int arrayLength) {
	PRINT_FUNC_NAME;

	const int numElement = numArray*arrayLength;
	uint32* h_data = (uint32*)malloc(sizeof(uint32)*numElement);
	for(int i = 0; i < numElement; i++) {
		h_data[i] = rand();	
	}

	printf("numArray = %d, arrayLength = %d, numElement = %d\n", 
		numArray, arrayLength, numElement);

	unsigned int timer = 0;

	//GPU
	uint32* d_data = NULL;
	GPUMALLOC((void**)&d_data, sizeof(uint32)*numElement);
	TOGPU(d_data, h_data, sizeof(uint32)*numElement);
	const int numThread = 256;
	//const int numArrayPerBlock = numThread/arrayLength;
	const int numBlock = 1200;
	const int sharedMemSize = sizeof(uint32)*numThread;
	printf("numBlock = %d, numThread = %d\n", numBlock, numThread);
	startTimer(&timer);
		batchSort_kernel<<<numBlock, numThread, sharedMemSize>>>(d_data, numArray, arrayLength, arrayLength);
		cutilCheckMsg("batchSort_kernel");
		cutilSafeCall(cudaThreadSynchronize());
	double sec = endTimer(&timer, "batchSort_kernel")/1000.0f;
	uint32* gpu_data = (uint32*)malloc(sizeof(uint32)*numElement);
	FROMGPU(gpu_data, d_data, sizeof(uint32)*numElement);


	//CPU
	uint32* gold_data = h_data;
	for(int arrayId = 0; arrayId < numArray; arrayId++) {
		sort(gold_data + arrayId*arrayLength, gold_data + (arrayId + 1)*arrayLength);
	}

	//check results
	simpleCheck(gold_data, gpu_data, numElement);

	//output 
	printf("sortintg throughput: %.3f M elements per sec \n", numArray*arrayLength/(1024.0f*1024.0f)/sec);
	
	GPUFREE(d_data);
	free(h_data);
}

/*

int main(int argc, char** argv) {

	cudaSetDevice(3);

	const int numElement = 16*1024*1024;
	for(int arrayLength = 8; arrayLength <= 256; arrayLength *= 2) {
		const int numArray = numElement/arrayLength;
		//test_bandwidth(numArray, arrayLength);
		test_batchSort(numArray, arrayLength);
		//test_batchSort2(numArray, arrayLength);
		printLine();
	}

//	test_batchSort(20, 128);

	return EXIT_SUCCESS;
}

*/

#endif /*__BATCH_SORT_CU__*/

