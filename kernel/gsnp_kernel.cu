#ifndef __GSNP_KERNEL_CU__
#define __GSNP_KERNEL_CU__

#include "gsnp.h"
#include "cuda_header.h"
#include <thrust/device_ptr.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>
#include "likelihood_kernel.cu"
#include "update_kernel.cu"
#include "posterior_kernel.cu"
#include "recycle_kernel.cu"

using namespace thrust;



void GPUInit() {
	cutilSafeCall(cudaSetDevice(DEVICE_ID));
}

void GSNPKernelInit(const uint64 winSize, const int readLength, 
		    const rate_t pcr_dependency, const rate_t global_dependency,
		    const rate_t* h_pmatrix,
		    const rate_t* p_prior, const rate_t* p_rank
){

	likelihoodKernelInit(winSize, readLength, pcr_dependency, global_dependency, h_pmatrix);
	updateKernelInit(winSize, readLength);
	posteriorKernelInit(p_prior, p_rank, winSize, readLength);
}


void GSNPKernelFinalize() {
	likelihoodKernelFinalize();
	updateKernelFinalize();
	posteriorKernelFinalize();

	cutilSafeCall(cudaThreadExit());
}


#endif /* __GSNP_KERNEL_CU__ */
