#ifndef __SEG_SORT_CU__
#define __SEG_SORT_CU__


#include <algorithm>
#include <gsnp_util.h>
#include "cuda_header.cu"

using namespace std;

#define BIT_SORT_START (8) 	//sorting from the bit 8
#define NUM_BIT_SORT (17)	//sorting bits from 8 to 24, 17 bits in total
#define DATA_IDX(arrayId, elementId) (elementId*pitch + arrayId)



__device__
void swap(uint32* a, uint32* b) {
	uint32 tmp = *a;
	*a = *b;
	*b = tmp;
}

__device__
void bubbleSort(uint32* data, const uint32 numElement, const uint32 pitch) {
	uint32 data_j_1 = 0;
	uint32 data_j = 0;

	for(uint32 i = 0; i < numElement; i++) {
		for(uint32 j = numElement - 1; j > i; j--) {

			data_j_1 = data[(j - 1)*pitch];
			data_j = data[j*pitch];

			if(data_j_1 > data_j) {
				//swap(&data_j_1, &data_j);	
	                        data[(j - 1)*pitch] = data_j;
        	                data[j*pitch] = data_j_1;	
			}
			__syncthreads();

			//data[(j - 1)*pitch] = data_j_1;
			//data[j*pitch] = data_j;
		}
	}
}

__global__
void segSortBubble_kernel(uint32* d_data, const uint32 numArray, const uint32 arraySize, 
			const uint32 pitch, const uint32* d_arraySizePerArray) {

	const uint32 numTotalThread = NUM_TOTAL_THREAD;
	const uint32 globalThreadOffset = GLOBAL_THREAD_OFFSET;

	for(uint32 arrayId = globalThreadOffset; arrayId < numArray; arrayId += numTotalThread) {
		bubbleSort(d_data + arrayId, d_arraySizePerArray[arrayId], pitch);
	}
}



//TODO: fix the bank conflict
//TODO: use multiple pass scatter
//TODO: calculate a block size first
__global__
void radixSort_kernel1(uint32* d_out, const uint32* d_in, 
			const uint32 numArray, 
			const uint32 arraySize, const uint32* d_arraySizePerArray,
			const uint32 pitch, 
			const uint32 mask, const uint32 numBitRightShift, const uint32 counterArraySize) {
	const uint32 numTotalThread = NUM_TOTAL_THREAD;
	const uint32 globalThreadOffset = GLOBAL_THREAD_OFFSET;
	extern __shared__ uint8 s_data[];
	uint8* l_sdata = s_data + threadIdx.x*(sizeof(uint16)*counterArraySize + sizeof(uint32)*counterArraySize);
	uint16* s_count = (uint16*)(l_sdata);
	uint32* s_offset = (uint32*)(l_sdata + sizeof(uint16)*counterArraySize);
	uint32 size = 0;
	const uint32 numPass = (uint32)ceil(numArray/(double)numTotalThread);
	uint32 element = 0;


	//for(uint32 arrayId = globalThreadOffset; arrayId < numArray; arrayId += numTotalThread) {
	for(uint32 passId = 0; passId < numPass; passId++) {

		const uint32 arrayId = passId*numTotalThread + globalThreadOffset;


		if(arrayId < numArray) {
 
			size = d_arraySizePerArray[arrayId];
			//0. set the counter as 0
			for(uint32 i = 0; i < counterArraySize; i++) {
			        s_count[i] = 0;
	       	 	}
		}

		__syncthreads();

		if(arrayId < numArray) {
			//1. count the number of elements in each partition, store in the s_count
			for(uint32 i = 0; i < arraySize; i++) {
				__syncthreads();
				element = d_in[DATA_IDX(arrayId, i)];
							
				if(i < size) {
					const uint16 pid = (element&mask)>>numBitRightShift;
					s_count[pid]++;
				}
			}		
		}


		if(arrayId < numArray) {
			//2. do the prefix sum on s_count, store the result in s_offset;
			s_offset[0] = 0;
			for(uint32 i = 1; i < counterArraySize; i++) {
				s_offset[i] = s_count[i - 1] + s_offset[i - 1];
			}		
			//3. write to the correct possition
			//3.1 clean the count to 0 again
			for(uint32 i = 0; i < counterArraySize; i++) {
		                s_count[i] = 0;
		        }
		}

		__syncthreads();

		if(arrayId < numArray) {
			//3.2 write to the correct positions
		        for(uint32 i = 0; i < arraySize; i++) {
				__syncthreads();
				element = d_in[DATA_IDX(arrayId, i)];

				if(i < size) { 
		                	const uint16 pid = (element&mask)>>numBitRightShift;
			
					d_out[pitch*(s_offset[pid] + s_count[pid]) + arrayId] = element;
					//d_out[i] = element;
			
		                	s_count[pid]++;
				}
		        }
	
		}
	}	
}


__global__
void radixSortPos_kernel(uint32* d_buf, const uint32* d_in,
                        const uint32 numArray,
                        const uint32 arraySize, const uint32* d_arraySizePerArray,
                        const uint32 pitch,
                        const uint32 mask, const uint32 numBitRightShift, const uint32 counterArraySize) {
        const uint32 numTotalThread = NUM_TOTAL_THREAD;
        const uint32 globalThreadOffset = GLOBAL_THREAD_OFFSET;
        extern __shared__ uint8 s_data[];
        uint8* l_sdata = s_data + threadIdx.x*(sizeof(uint16)*counterArraySize + sizeof(uint32)*counterArraySize);
        uint16* s_count = (uint16*)(l_sdata);
        uint32* s_offset = (uint32*)(l_sdata + sizeof(uint16)*counterArraySize);
        uint32 size = 0;
        const uint32 numPass = (uint32)ceil(numArray/(double)numTotalThread);
        uint32 element = 0;


        //for(uint32 arrayId = globalThreadOffset; arrayId < numArray; arrayId += numTotalThread) {
        for(uint32 passId = 0; passId < numPass; passId++) {

                const uint32 arrayId = passId*numTotalThread + globalThreadOffset;


                if(arrayId < numArray) {

                        size = d_arraySizePerArray[arrayId];
                        //0. set the counter as 0
                        for(uint32 i = 0; i < counterArraySize; i++) {
                                s_count[i] = 0;
                        }
                }

                __syncthreads();

                if(arrayId < numArray) {
                        //1. count the number of elements in each partition, store in the s_count
                        for(uint32 i = 0; i < arraySize; i++) {
                                __syncthreads();
                                element = d_in[DATA_IDX(arrayId, i)];

                                if(i < size) {
                                        const uint16 pid = (element&mask)>>numBitRightShift;
                                        s_count[pid]++;
                                }
                        }
                }


                if(arrayId < numArray) {
                        //2. do the prefix sum on s_count, store the result in s_offset;
                        s_offset[0] = 0;
                        for(uint32 i = 1; i < counterArraySize; i++) {
                                s_offset[i] = s_count[i - 1] + s_offset[i - 1];
                        }
                        //3. write to the correct possition
                        //3.1 clean the count to 0 again
                        for(uint32 i = 0; i < counterArraySize; i++) {
                                s_count[i] = 0;
                        }
                }

                __syncthreads();

                if(arrayId < numArray) {
                        //3.2 write to the correct positions
                        for(uint32 i = 0; i < arraySize; i++) {
                                __syncthreads();
                                element = d_in[DATA_IDX(arrayId, i)];

                                if(i < size) {
                                        const uint16 pid = (element&mask)>>numBitRightShift;

                                        //d_out[pitch*(s_offset[pid] + s_count[pid]) + arrayId] = element;
					d_buf[pitch*i + arrayId] = s_offset[pid] + s_count[pid];

                                        s_count[pid]++;
                                }
                        }

                }
        }
}



//naive scatter
__global__
void radixSortScatter_kernel1(uint32* d_out, const uint32* d_in, const uint32* d_offset, 
			const uint32 numArray, 
			const uint32 arraySize, const uint32* d_arraySizePerArray, const uint32 pitch) {
	const uint32 numTotalThread = NUM_TOTAL_THREAD;
	const uint32 globalThreadOffset = GLOBAL_THREAD_OFFSET;
	uint32 size = 0;

	for(uint32 arrayId = globalThreadOffset; arrayId < numArray; arrayId += numTotalThread) {
		size = d_arraySizePerArray[arrayId];

		//naive scatter
		for(uint32 i = 0; i < size; i++) {
			uint32 pos = pitch*d_offset[i*pitch + arrayId] + arrayId;
			d_out[pos] = d_in[i*pitch + arrayId];
		}
	}	
}


//multipass scatter
__global__
void radixSortScatter_kernel2(uint32* d_out, const uint32* d_in, const uint32* d_offset,
                        const uint32 numArray, const uint32 arraySize, const uint32* d_arraySizePerArray, 
			const uint32 pitch, 
			const uint32 winSize) {
        const uint32 numTotalThread = NUM_TOTAL_THREAD;
        const uint32 globalThreadOffset = GLOBAL_THREAD_OFFSET;
        uint32 size = 0;

        for(uint32 arrayId = globalThreadOffset; arrayId < numArray; arrayId += numTotalThread) {
                size = d_arraySizePerArray[arrayId];
		const uint32 numWin = (uint32)ceil(size/(double)winSize);

		for(uint32 winId = 0; winId < numWin; winId++) {
			const uint32 start = winId*winSize;
			const uint32 end = start + winSize;

			for(uint32 i = 0; i < arraySize; i++) {
				const uint32 offset = d_offset[i*pitch + arrayId];
				const uint32 element = d_in[i*pitch + arrayId];
			
				if(i < size) {
					if((offset >= start) && (offset < end)) {
						d_out[pitch*offset + arrayId] = element;
					}
				}
			}
		}
     	}
}




//shared mem. + multipass scatter
__global__
void radixSortScatter_kernel3(uint32* d_out, const uint32* d_in, const uint32* d_offset,
                        const uint32 numArray, const uint32 arraySize, const uint32* d_arraySizePerArray,
                        const uint32 pitch,
                        const uint32 winSize) {
        const uint32 numTotalThread = NUM_TOTAL_THREAD;
        const uint32 globalThreadOffset = GLOBAL_THREAD_OFFSET;
	const uint32 numWin = (uint32)ceil(arraySize/(double)winSize);
        uint32 size, arrayId, start, end, offset, element;
	extern __shared__ uint32 s_data2[];
	uint32* l_data = s_data2 + threadIdx.x*winSize;
	//uint32 l_data[32];

	for(uint32 arrayId = globalThreadOffset; arrayId < numArray; arrayId += numTotalThread) {

		size = d_arraySizePerArray[arrayId];

		for(uint32 winId = 0; winId < numWin; winId++) {

			start = winId*winSize;
			end = start + winSize;

			for(uint32 i = 0; i < arraySize; i++) {

				offset = d_offset[i*pitch + arrayId];
				element = d_in[i*pitch + arrayId];
				
				//if(i < size) {
				//	d_out[i*pitch + arrayId] = element + winId;
				//}

				if((i < size) && (offset >= start) && (offset < end)) {
					l_data[offset - start] = element;
				}
			}

			__syncthreads();
			
			for(uint32 i = 0; i < winSize; i++) {
				d_out[(start + i)*pitch + arrayId] = l_data[i];
			}
		}
	}	
}



void radixSort_v2(uint32* d_out, uint32* d_data, uint32* d_buf,
		const uint32 numArray, 
		const uint32 arraySize, const uint32* d_arraySizePerArray,
		const uint32 pitch, const uint32 unitNumBit,
		uint32 numBlock, uint32 numThread) {

        PRINT_FUNC_NAME;

        assert(arraySize <= 65536 /*UINT16_MAX*/);      //since we use uint16 to store the count
        assert(unitNumBit <= 16);                       //since we use the count array to store the partition id

        const uint32 numPass = ceil(NUM_BIT_SORT/(double)unitNumBit);
        const uint32 counterSizePerThread = pow(2, unitNumBit);
        const uint32 sharedMemSize = numThread*(sizeof(uint16)*counterSizePerThread + sizeof(uint32)*counterSizePerThread);
        printf("\tnumBlock = %d, numThread = %d\n", numBlock, numThread);
        printf("\tunitNumBit = %d, numPass = %d\n", unitNumBit, numPass);
        printf("\tshared memory size allocated: %.3f KB\n", sharedMemSize/1024.0f);
        assert(sharedMemSize <= MAX_SHARED_MEM_SIZE);


        //get the mask
        uint32 mask = 1;
        for(uint32 i = 1; i < unitNumBit; i++) {
                mask = mask<<1;
                mask |= 1;
        }
        mask = mask<<BIT_SORT_START;
        uint32 numBitRightShift = BIT_SORT_START;

        //ping-pong buffer between d_out and d_data
        uint32* d_read = NULL;
        uint32* d_write = NULL;
	unsigned int timer = 0;

        for(uint32 passId = 0; passId < numPass; passId++) {
                printf("pass %d: mask = %0x, numBitShiftRight = %d\n", passId, mask, numBitRightShift);

                //we use buffer ping-pong
                if(passId%2 == 0) {
                        d_read = d_data;
                        d_write = d_out;
                } else {
                        d_read = d_out;
                        d_write = d_data;
                }


                //radixSort_kernel1<<<numBlock, numThread, sharedMemSize>>>
                //        (d_write, d_read, numArray, arraySize, d_arraySizePerArray,
                //        pitch, mask, numBitRightShift, counterSizePerThread);
                //cutilCheckMsg("radixSort_kernel1");

		//1. get the output positions
		startTimer(&timer);
                radixSortPos_kernel<<<numBlock, numThread, sharedMemSize>>>
                        (d_buf, d_read, numArray, arraySize, d_arraySizePerArray,
                        pitch, mask, numBitRightShift, counterSizePerThread);
		cutilCheckMsg("radixSortPos_kernel");
		cutilSafeCall(cudaThreadSynchronize());
		endTimer(&timer, "\tradixSortPos_kernel");
		

		//2. scatter according to the output positions, scattering into d_write
		
		//navie scatter kernel
		/*timer = 0;
		startTimer(&timer);
		radixSortScatter_kernel1<<<numBlock, numThread>>>
			(d_write, d_read, d_buf, numArray, arraySize, d_arraySizePerArray, pitch);
		cutilCheckMsg("radixSortScatter_kernel");
		cutilSafeCall(cudaThreadSynchronize());
		endTimer(&timer, "\tradixSortScatter_kernel1");
		*/


		//multi-pass scatter kernel
		/*timer = 0;
		const uint32 winSize = 64;
		startTimer(&timer);
		radixSortScatter_kernel2<<<numBlock, numThread>>>
			(d_write, d_read, d_buf,
                        numArray, arraySize, d_arraySizePerArray,
                        pitch, winSize);
		cutilCheckMsg("radixSortScatter_kernel2");
		cutilSafeCall(cudaThreadSynchronize());
		endTimer(&timer, "\tradixSortScatter_kernel2");
		*/


		
		//multipass + shared memory scatter kernel
		timer = 0;
		const uint32 winSize = 32;
		numBlock = 512;
		numThread = 64;
		const uint32 scatterSharedMemSize = sizeof(uint32)*winSize*numThread;
		const uint32 numWin = (uint32)ceil(arraySize/(double)winSize);
		printf("winSize = %d, numWin = %d, numBlock = %d, numThread = %d, scatterSharedMemSize = %.3f KB\n", 
			winSize, numWin, numBlock, numThread, scatterSharedMemSize/1024.0f); 
		assert(scatterSharedMemSize <= MAX_SHARED_MEM_SIZE);
		assert(numArray%(numBlock*numThread) == 0); //complete pass
		startTimer(&timer);
		radixSortScatter_kernel3<<<numBlock, numThread, scatterSharedMemSize>>>
			(d_write, d_read, d_buf, 
			numArray, arraySize, d_arraySizePerArray, 
			pitch, winSize);
		cutilCheckMsg("radixSortScatter_kernel3");
		cutilSafeCall(cudaThreadSynchronize());
		endTimer(&timer, "\tradixSortScatter_kernel3");
		



                mask = mask<<unitNumBit;
                numBitRightShift += unitNumBit;
        }


        if(numPass%2 == 0) {
                GPUTOGPU(d_out, d_data, sizeof(uint32)*pitch*arraySize);
        }

}




/**
* the radix sort, v1, we just write to the globle memory
* note, we only need to sort from bit 9 to 25, 17 bits in total;
*/

void radixSort_v1(uint32* d_out, uint32* d_data, 
		  const uint32 numArray, 
		  const uint32 arraySize, const uint32* d_arraySizePerArray,
		  const uint32 pitch, const uint32 unitNumBit,
		  const uint32 numBlock, const uint32 numThread) {

	PRINT_FUNC_NAME;

	assert(arraySize <= 65536 /*UINT16_MAX*/); 	//since we use uint16 to store the count
	assert(unitNumBit <= 16);			//since we use the count array to store the partition id

	const uint32 numPass = ceil(NUM_BIT_SORT/(double)unitNumBit);
	const uint32 counterSizePerThread = pow(2, unitNumBit);
	const uint32 sharedMemSize = numThread*(sizeof(uint16)*counterSizePerThread + sizeof(uint32)*counterSizePerThread);
	printf("\tnumBlock = %d, numThread = %d\n", numBlock, numThread);
	printf("\tunitNumBit = %d, numPass = %d\n", unitNumBit, numPass);
	printf("\tshared memory size allocated: %.3f KB\n", sharedMemSize/1024.0f);
	assert(sharedMemSize <= MAX_SHARED_MEM_SIZE);


	//get the mask
	uint32 mask = 1;
	for(uint32 i = 1; i < unitNumBit; i++) {
		mask = mask<<1;
		mask |= 1;
	}
	mask = mask<<BIT_SORT_START;
	uint32 numBitRightShift = BIT_SORT_START;
	
	//ping-pong buffer between d_out and d_data
	uint32* d_read = NULL;
	uint32* d_write = NULL;

	for(uint32 passId = 0; passId < numPass; passId++) {
		printf("pass %d: mask = %0x, numBitShiftRight = %d\n", passId, mask, numBitRightShift);

		//we use buffer ping-pong
		if(passId%2 == 0) {
			d_read = d_data;
			d_write = d_out;
		} else {
			d_read = d_out;
			d_write = d_data;
		}

		radixSort_kernel1<<<numBlock, numThread, sharedMemSize>>>
			(d_write, d_read, numArray, arraySize, d_arraySizePerArray,
			pitch, mask, numBitRightShift, counterSizePerThread);
		cutilCheckMsg("radixSort_kernel1");

		mask = mask<<unitNumBit;
		numBitRightShift += unitNumBit;
	}


	if(numPass%2 == 0) {
		GPUTOGPU(d_out, d_data, sizeof(uint32)*pitch*arraySize);		
	}
}


__device__
void insertionSort(uint32* data, const uint32 n, const uint32 pitch = 1) {
	int i, j;
	uint32 temp;	

	for(i = 1; i < n; i++) {
		temp = data[i*pitch];
		j = i - 1;
		while(temp < data[j*pitch] && j >= 0) {
			data[(j + 1)*pitch] = data[j*pitch];
			j = j - 1;
		}
		data[pitch*(j + 1)] = temp;
	}
}

__global__
void cacheSort_kernel(uint32* d_data, const uint32 numArray, 
			const uint32 arraySize, const uint32* d_arraySize,
			const uint32 pitch) {
	const uint32 globalThreadOffset = GLOBAL_THREAD_OFFSET;
	const uint32 numTotalThread = NUM_TOTAL_THREAD;
	extern __shared__ uint32 s_data3[];
	//uint32* l_data = s_data3 + threadIdx.x*arraySize;
	uint32* l_data = s_data3 + threadIdx.x;
	const uint32 l_pitch = blockDim.x;
	
	for(uint32 i = 0; i < arraySize; i++) {
		l_data[i*l_pitch] = d_data[i*pitch + globalThreadOffset];
	}
	__syncthreads();


	//bubbleSort(l_data, d_arraySize[globalThreadOffset], l_pitch);
	insertionSort(l_data, d_arraySize[globalThreadOffset], l_pitch);

	__syncthreads();

        for(uint32 i = 0; i < arraySize; i++) {
                d_data[i*pitch + globalThreadOffset] = l_data[i*l_pitch];
        }
}


/**
* to address the jth element in array i
* the address is pitch*j + i
*/
void gpu_segSort(uint32* d_data, const uint32 numArray, const uint32 arraySize, 
		const uint32 pitch, const uint32* d_arraySizePerArray,
		uint32 numBlock = 2400, uint32 numThread = 256) {

	PRINT_FUNC_NAME;
	printf("numArray = %u, arraySize = %u, pitch = %u\n", numArray, arraySize, pitch);

	unsigned int timer = 0;

	/* the simple bubble sort */
	/*timer = 0;
	cutilSafeCall(cudaThreadSynchronize());
	startTimer(&timer);
	segSortBubble_kernel<<<numBlock, numThread>>>(d_data, numArray, arraySize, pitch, d_arraySizePerArray);
	cutilCheckMsg("segSortBubble_kernel");
	cutilSafeCall(cudaThreadSynchronize());
	endTimer(&timer, "segSortBubble_kernel");
	*/
	


	/* the radix sort, v1 */
	uint32 unitNumBit = 4;
	numBlock = 1200;
	numThread = 128;
	uint32* d_out = NULL;
	GPUMALLOC((void**)&d_out, sizeof(uint32)*pitch*arraySize);
	cutilSafeCall(cudaMemset(d_out, 0, sizeof(uint32)*pitch*arraySize));
	timer = 0;
	startTimer(&timer);
		radixSort_v1(d_out, d_data, numArray, arraySize, d_arraySizePerArray, pitch, unitNumBit, numBlock, numThread);
		cutilSafeCall(cudaThreadSynchronize());
	endTimer(&timer, "radixSort_v1");
	GPUTOGPU(d_data, d_out, sizeof(uint32)*pitch*arraySize);
	GPUFREE(d_out);
	


	/* the radix sort, v2 */
        /*uint32 unitNumBit = 4;
        numBlock = 1200;
        numThread = 128;
        uint32* d_out = NULL;
        GPUMALLOC((void**)&d_out, sizeof(uint32)*pitch*arraySize);
	uint32* d_buf = NULL;
	GPUMALLOC((void**)&d_buf, sizeof(uint32)*pitch*arraySize);
        cutilSafeCall(cudaMemset(d_out, 0, sizeof(uint32)*pitch*arraySize));
        startTimer(&timer);
        radixSort_v2(d_out, d_data, d_buf, numArray, arraySize, d_arraySizePerArray, pitch, unitNumBit, numBlock, numThread);
        cutilSafeCall(cudaThreadSynchronize());
        endTimer(&timer, "radixSort_v2");
        GPUTOGPU(d_data, d_out, sizeof(uint32)*pitch*arraySize);
        GPUFREE(d_out);
	GPUFREE(d_buf);*/


	/*the pure cache (registers or shared memory) based sort*/
	/*printTitle("cache-based sorting");
	numThread = 64;
	assert(numArray%numThread == 0);
	numBlock = numArray/numThread;
	const uint32 sharedMemSize = sizeof(uint32)*arraySize*numThread;
	printf("numBlock = %d, numThread = %d, arraySize = %d, sharedMemSize = %.3f KB\n",
		numBlock, numThread, arraySize, sharedMemSize/(1024.0f));
	assert(sharedMemSize <= MAX_SHARED_MEM_SIZE);
	uint32* d_out = NULL;
	GPUMALLOC((void**)&d_out, sizeof(uint32)*pitch*arraySize);
	timer = 0;
	startTimer(&timer);
		cacheSort_kernel<<<numBlock, numThread, sharedMemSize>>>
				(d_data, numArray, arraySize, d_arraySizePerArray, pitch);
		cutilCheckMsg("cacheSort_kernel");
		cutilSafeCall(cudaThreadSynchronize());
	endTimer(&timer, "cacheSort_kernel");*/
	
}


void printArrayGroup(uint32* data, const uint32 numArray, const uint32 arraySize, const uint32* arraySizePerArray, 
			const uint32 pitch = 1) {
	for(uint32 arrayId = 0; arrayId < numArray; arrayId++) {

		printf("Array %u (%u): \n", arrayId, arraySizePerArray[arrayId]);
		for(uint32 i = 0; i < arraySizePerArray[arrayId]; i++) {
			printf("%u\n", data[DATA_IDX(arrayId, i)]);
		}
		printLine();
	}
}


int uint32_cmp(const void * a, const void * b) {
	return ( *(uint32*)a - *(uint32*)b );
}

int main() {

	cudaSetDevice(3);
	
	uint32 numArray = 1024*64;
	uint32 arraySize = 120;
	uint32 pitch = numArray + 35;
	uint64 numTotalElement = pitch*arraySize;

	printf("numArray = %d, arraySize = %d, pitch = %d\n", 
		numArray, arraySize, pitch);
	printf("numTotalElement = %lu, memory allocated: %.3f MB\n", 
		numTotalElement, sizeof(uint32)*numTotalElement/(1024.0f*1024.0f));

	uint32* h_data = (uint32*)malloc(sizeof(uint32)*numTotalElement);
	memset(h_data, 0, sizeof(uint32)*numTotalElement);
	const uint32 mask = 0x01ffff00;
	uint32* arraySizePerArray = (uint32*)malloc(sizeof(uint32)*numArray);
	for(uint32 i = 0; i < numArray; i++) {
		arraySizePerArray[i] = arraySize/2 + rand()%(arraySize/2) - 1;
		assert(arraySizePerArray[i] <= arraySize);
		
		for(uint64 j = 0; j < arraySizePerArray[i]; j++) {
                	h_data[j*pitch + i] = rand_uint32()&mask; //we only need from bit 9 to bit 25     
        	}
	}
	uint32* d_data = NULL;
	GPUMALLOC((void**)&d_data, sizeof(uint32)*numTotalElement);
	TOGPU(d_data, h_data, sizeof(uint32)*numTotalElement);
	uint32* d_arraySizePerArray = NULL;
	GPUMALLOC((void**)&d_arraySizePerArray, sizeof(uint32)*numArray);
	TOGPU(d_arraySizePerArray, arraySizePerArray, sizeof(uint32)*numArray);

	//printTitle("before sorting");
	//printArrayGroup(h_data, numArray, arraySize, arraySizePerArray, pitch);


	
	gpu_segSort(d_data, numArray, arraySize, pitch, d_arraySizePerArray);
	uint32* gpu_result = (uint32*)malloc(sizeof(uint32)*numTotalElement);	
	FROMGPU(gpu_result, d_data, sizeof(uint32)*numTotalElement);
	//printTitle("after GPU sorting");
	//printArrayGroup(gpu_result, numArray, arraySize, arraySizePerArray, pitch);


	
	//gold sorting
	//convert to linear memory
	uint32* gold_result = (uint32*)malloc(sizeof(uint32)*numArray*arraySize);
	for(uint32 arrayId = 0; arrayId < numArray; arrayId++) {
		for(uint32 i = 0; i < arraySize; i++) {
			gold_result[arrayId*arraySize + i] = h_data[DATA_IDX(arrayId, i)];
		}
	}
	//sort for each array
	printTitle("Gold sort");
	unsigned int timer = 0;
	startTimer(&timer);
	for(uint32 arrayId = 0; arrayId < numArray; arrayId++) {
		uint32* start = gold_result + arrayId*arraySize;
		sort(start, start + arraySizePerArray[arrayId]);

		
		//for(uint32 i = 0; i < arraySizePerArray[arrayId]; i++) {
		//	printf("%u\n", gold_result[arrayId*arraySize + i]);
		//}
		//printLine();
		
	}
	endTimer(&timer, "Gold sort");

	
	//check result
	printf("check result...\n");
	for(uint32 arrayId = 0; arrayId < numArray; arrayId++) {
		for(uint32 i = 0; i < arraySizePerArray[arrayId]; i++) {
			if(gpu_result[DATA_IDX(arrayId, i)] != gold_result[arrayId*arraySize + i]) {
				ERROR_EXIT("FAILED !!!");
			}

		}
	}
	printf("\tPASS :)\n");

	
	free(h_data);
	free(gold_result);
	free(arraySizePerArray);
	GPUFREE(d_data);
	GPUFREE(d_arraySizePerArray);

	

	return EXIT_SUCCESS;
}


#endif /* __SEG_SORT_CU__ */


