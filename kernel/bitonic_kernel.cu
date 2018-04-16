
#ifndef _BITONIC_KERNEL_H_
#define _BITONIC_KERNEL_H_

#define NUM    256


__device__ inline void swap(uint32 & a, uint32 & b)
{
	// Alternative swap doesn't use a temporary register:
//	 a ^= b;
//	 b ^= a;
//	 a ^= b;
	
    uint32 tmp = a;
    a = b;
    b = tmp;
}


__global__
void batchSort_kernel(uint32* d_data, const int arrayPitch, const int numArray, 
		const uint32* d_arrayLength, const int low, const int up) {

        extern __shared__ uint32 shared[];
        const unsigned int tid = threadIdx.x;
        const int numChunk = ceil(numArray/(double)gridDim.x);



	for(int chunkId = 0; chunkId < numChunk; chunkId++) {
		const int arrayId = chunkId*gridDim.x + blockIdx.x;

		if(arrayId < numArray) {
			const int arrayLength = d_arrayLength[arrayId];
	                if(!((arrayLength > low) && (arrayLength <= up))) {
	                        continue; //skip this block
	                }	
		}

		if(arrayId < numArray) {
			shared[tid] = d_data[arrayId*arrayPitch + tid];
		}
		__syncthreads();

		for (unsigned int k = 2; k <= up; k <<= 1) {
        
			for (unsigned int j = k / 2; j > 0; j >>= 1) {
				unsigned int ixj = tid ^ j;

				if (ixj > tid) {
					if ((tid & k) == 0) {
						if (shared[tid] > shared[ixj]) {
							 swap(shared[tid], shared[ixj]);
						}
					} else {
						if (shared[tid] < shared[ixj]) {
							swap(shared[tid], shared[ixj]);
						}	
					}
				}
				__syncthreads();
			}
		}
		
		__syncthreads();
		if(arrayId < numArray) {
			d_data[arrayId*arrayPitch + tid] = shared[tid];
		}
	}

}


__global__
void bitonicSort_kernel(uint32* d_data, const int numArray, const int arrayLength) {

	extern __shared__ uint32 shared[];
	const unsigned int tid = threadIdx.x;
	const int numChunk = ceil(numArray/(double)gridDim.x);

	for(int chunkId = 0; chunkId < numChunk; chunkId++) {
		const int arrayId = chunkId*gridDim.x + blockIdx.x;

		if(arrayId < numArray) {
			shared[tid] = d_data[arrayId*arrayLength + tid];
		}
		__syncthreads();

		for (unsigned int k = 2; k <= arrayLength; k <<= 1) {
        
			for (unsigned int j = k / 2; j > 0; j >>= 1) {
				unsigned int ixj = tid ^ j;

				if (ixj > tid) {
					if ((tid & k) == 0) {
						if (shared[tid] > shared[ixj]) {
							 swap(shared[tid], shared[ixj]);
						}
					} else {
						if (shared[tid] < shared[ixj]) {
							swap(shared[tid], shared[ixj]);
						}	
					}
				}
				__syncthreads();
			}
		}
		
		__syncthreads();
		if(arrayId < numArray) {
			d_data[arrayId*arrayLength + tid] = shared[tid];
		}
	}
}



__global__
void batchSortMP2_kernel(uint32* d_data, const int arrayPitch, 
			const int* d_arrayId, const int numArray, 
			const int arrayLength, const int maxArrayId) {
        const int tid = threadIdx.x;
        const int arrayOffset = tid%arrayLength;
        const int numArrayPerBlock = blockDim.x/arrayLength;
        const int numArrayPerChunk = numArrayPerBlock*gridDim.x;
        const int numChunk = ceil(numArray/(double)numArrayPerChunk);
        int globalArrayId, blockArrayId, arrayBase;

        extern __shared__ uint32 shared[];

        for(int chunkId = 0; chunkId < numChunk; chunkId++) {
                blockArrayId = tid/arrayLength;
                globalArrayId = d_arrayId[chunkId*numArrayPerChunk + GLOBAL_THREAD_OFFSET/arrayLength];
                arrayBase = blockArrayId*arrayLength;



                //load the data to the shared memory
                if(globalArrayId < maxArrayId) {
                        shared[arrayBase + arrayOffset] = d_data[globalArrayId*arrayPitch + arrayOffset];
                }
                __syncthreads();


                //do the bitonic sort
                for(unsigned int k = 2; k <= arrayLength; k <<= 1) {
                        for(unsigned int j = k/2; j > 0; j >>= 1) {
                                unsigned int ixj = arrayOffset^j;
                                if(ixj > arrayOffset) {
                                        if((arrayOffset & k) == 0) {
                                                if(shared[arrayBase + arrayOffset] >
                                                  shared[arrayBase + ixj]) {
                                                        swap(shared[arrayBase + arrayOffset],
                                                        shared[arrayBase + ixj]);
                                                }
                                        } else {
                                                if(shared[arrayBase + arrayOffset] <
                                                  shared[arrayBase + ixj]) {
                                                        swap(shared[arrayBase + arrayOffset],
                                                        shared[arrayBase + ixj]);
                                                }
                                        }
                                }
                                __syncthreads();
                        }
                }


                __syncthreads();
                //copy the data back
                if(globalArrayId < maxArrayId) {
                        d_data[globalArrayId*arrayPitch + arrayOffset] = shared[arrayBase + arrayOffset];
                }
        }
}


__global__
void batchSort_kernel(uint32* d_data, const int numArray, 
			const int arrayLength, const int arrayPitch) {
	const int tid = threadIdx.x;
	const int arrayOffset = tid%arrayLength;
	const int numArrayPerBlock = blockDim.x/arrayLength;
	const int numArrayPerChunk = numArrayPerBlock*gridDim.x;
	const int numChunk = ceil(numArray/(double)numArrayPerChunk);
	int globalArrayId, blockArrayId, arrayBase;
	
	extern __shared__ uint32 shared[];


	for(int chunkId = 0; chunkId < numChunk; chunkId++) {
       		blockArrayId = tid/arrayLength;
		globalArrayId = chunkId*numArrayPerChunk + GLOBAL_THREAD_OFFSET/arrayLength;
		arrayBase = blockArrayId*arrayLength;



		//load the data to the shared memory
		if(globalArrayId < numArray) {
			shared[arrayBase + arrayOffset] = d_data[globalArrayId*arrayPitch + arrayOffset];
		}
		__syncthreads();



		//do the bitonic sort
		for(unsigned int k = 2; k <= arrayLength; k <<= 1) {
	        	for(unsigned int j = k/2; j > 0; j >>= 1) {
	                	unsigned int ixj = arrayOffset^j;
       	                	if(ixj > arrayOffset) {
        	                       	if((arrayOffset & k) == 0) {
                	                       	if(shared[arrayBase + arrayOffset] >
                        	                  shared[arrayBase + ixj]) {
                                                	swap(shared[arrayBase + arrayOffset],
                                       	             	shared[arrayBase + ixj]);
                                       		}
                                	} else {
                                       		if(shared[arrayBase + arrayOffset] <
                                       	   	  shared[arrayBase + ixj]) {
                                       	        	swap(shared[arrayBase + arrayOffset],
                                       	            	shared[arrayBase + ixj]);
                                       		}
                                	}
                        	}
                        	__syncthreads();
                	}
        	}


		__syncthreads();
		//copy the data back
		if(globalArrayId < numArray) {
			d_data[globalArrayId*arrayPitch + arrayOffset] = shared[arrayBase + arrayOffset];
		}
	}
}







#endif // _BITONIC_KERNEL_H_
