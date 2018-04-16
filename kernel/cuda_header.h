
#ifndef _CUDA_HEADER_CU_
#define _CUDA_HEADER_CU_

#include <stdio.h>
#include <stdlib.h>
#include <cutil_inline.h>


/** for CUDA 2.0 and below versions*/
#ifdef CUDA_2
    #define cutilCheckMsg CUT_CHECK_ERROR
    #define cutilSafeCall CUDA_SAFE_CALL
#endif


/* kernel macros */
#define NUM_TOTAL_THREAD (gridDim.x*blockDim.x)
#define GLOBAL_THREAD_OFFSET (blockDim.x*blockIdx.x + threadIdx.x)

/** macro utility */
static double gMemSize = 0.0f;
#define GPUMALLOC(D_DATA, MEM_SIZE) cutilSafeCall(cudaMalloc(D_DATA, MEM_SIZE)); gMemSize += MEM_SIZE
#define TOGPU(D_DATA, H_DATA, MEM_SIZE) cutilSafeCall(cudaMemcpy(D_DATA, H_DATA, MEM_SIZE, cudaMemcpyHostToDevice))
#define FROMGPU( H_DATA, D_DATA, MEM_SIZE ) cutilSafeCall(cudaMemcpy( H_DATA, D_DATA, MEM_SIZE, cudaMemcpyDeviceToHost))
#define GPUTOGPU( DEST, SRC, MEM_SIZE ) cutilSafeCall(cudaMemcpy( DEST, SRC, MEM_SIZE, cudaMemcpyDeviceToDevice ))
#define GPUFREE( MEM ) cutilSafeCall(cudaFree(MEM));
#define PRINT_G_MEM_SIZE printf("gMemSize = %f\n", gMemSize);


/* hardware specification */
#define MAX_SHARED_MEM_SIZE (16000) // in bytes
#define MAX_NUM_THREAD (512)
#define MAX_NUM_BLOCK (65536)

/** timer utility */
inline void startTimer( unsigned int*  timer ) {
    	*timer = 0;
	CUT_SAFE_CALL( cutCreateTimer( timer));
	CUT_SAFE_CALL( cutStartTimer( *timer));
}

inline float endTimer( unsigned int*  timer, char* title ) {
    	CUT_SAFE_CALL( cutStopTimer( *timer));
	float ms = cutGetTimerValue(*timer);	
	printf( "*** %s processing time: %.6f sec ***\n", title, ms/1000.0);
	CUT_SAFE_CALL( cutDeleteTimer( *timer));

	return ms;
}


#endif

