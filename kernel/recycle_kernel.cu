#ifndef RECYCLE_KERNEL_
#define RECYCLE_KERNEL_
#include "gsnp.h"
#include "common_kernel.h"
#include "cuda_header.h"
#include "likelihood_kernel.cu"
#include "update_kernel.cu"
//#include "bitonic_sort.cu"
#include "common.h"
#include <string>
void cuda_init(
        ubit64_t win_size,
        ubit64_t read_len,
	ubit64_t start,
	int* posArr,
	const uint32 numBlock, const uint32 numThread)
{
	for(int i = 0; i < read_len + win_size; i++){
		posArr[i] = i + start;
	}

	cutilSafeCall(cudaMemset(d_sparseBaseInfo, 0xffff, sizeof(uint32) * (read_len + win_size) * MAX_NUM_NON_ZERO));
        cutilSafeCall(cudaMemset(d_oriArr, 0xFF, sizeof(char) * (read_len + win_size)));
        cutilSafeCall(cudaMemset(d_depthArr, 0, sizeof(int) * (read_len + win_size)));
        cutilSafeCall(cudaMemset(d_repeat_timeArr, 0, sizeof(int) * (read_len + win_size)));
        cutilSafeCall(cudaMemset(d_dep_uniArr, 0, sizeof(int) * (read_len + win_size)));
        cutilSafeCall(cudaMemset(d_count_uniArr, 0, sizeof(int) * 4 * (read_len + win_size)));
        cutilSafeCall(cudaMemset(d_q_sumArr, 0, sizeof(int) * 4 * (read_len + win_size)));
        cutilSafeCall(cudaMemset(d_count_allArr, 0, sizeof(int) * 4 * (read_len + win_size)));
        cutilSafeCall(cudaMemset(d_numNonZeroPerSite, 0, sizeof(uint32) * (read_len + win_size)));

}

/*
__global__
void recycle_kernel(
	int* d_depthArr,
	int* d_repeat_timeArr,
	int* d_dep_uniArr,
	int* d_count_uniArr,
	int* d_q_sumArr,
	int* d_count_allArr,
	uint32* d_numNonZeroPerSite,
	uint32* d_sparseBaseInfo,
	ubit64_t win_size,
	ubit64_t read_len,
	int global_win_size,
	uint32 MAX_NUM_NON_ZERO		
)
{
        const uint32 numTotalThread = NUM_TOTAL_THREAD;
        const uint32 globalThreadOffset = GLOBAL_THREAD_OFFSET;
	for(int i = globalThreadOffset; i < read_len; i += numTotalThread){
		d_depthArr[i] = d_depthArr[i + win_size];
		d_repeat_timeArr[i] = d_repeat_timeArr[i + win_size];
		d_dep_uniArr[i] = d_dep_uniArr[i + win_size];
		d_numNonZeroPerSite[i] = d_numNonZeroPerSite[i + win_size];

		for(int j = 0; j < 4; j++){
			d_count_uniArr[i * 4 + j] = d_count_uniArr[(i + win_size) * 4 + j];
			d_q_sumArr[i * 4 + j] = d_q_sumArr[(i + win_size) * 4 + j];
			d_count_allArr[i * 4 + j] = d_count_allArr[(i + win_size) * 4 + j];
		}

		for(uint32 offset = 0; offset < MAX_NUM_NON_ZERO; offset++){
			d_sparseBaseInfo[BASEINFO_IDX(i, offset, win_size + read_len, MAX_NUM_NON_ZERO)] = 
				d_sparseBaseInfo[BASEINFO_IDX(i + global_win_size, offset, win_size + read_len, MAX_NUM_NON_ZERO)];
		}
	}

	for(int i = globalThreadOffset + read_len; i < read_len + win_size; i += numTotalThread){
		d_depthArr[i] = 0;
		d_repeat_timeArr[i] = 0;
		d_dep_uniArr[i] = 0;
		d_numNonZeroPerSite[i] = 0;

                for(int j = 0; j < 4; j++){
                        d_count_uniArr[i * 4 + j] = 0;
                        d_q_sumArr[i * 4 + j] = 0;
                        d_count_allArr[i * 4 + j] = 0;
                }

		for(uint32 offset = 0; offset < MAX_NUM_NON_ZERO; offset++) {
                        d_sparseBaseInfo[BASEINFO_IDX(i, offset, numSite, MAX_NUM_NON_ZERO)] = BASEINFO_MAX;
                }
	}
}
*/


void cuda_recycle(
	ubit64_t win_size,
        ubit64_t read_len,
	int* posArr,
        const uint32 numBlock, const uint32 numThread)
{
        memcpy(posArr, posArr + win_size, read_len * sizeof(int));
        for(int i = read_len; i < read_len + win_size; i++){
                posArr[i] = posArr[i - 1] + 1;
        }


	GPUTOGPU(d_depthArr, d_depthArr + win_size, read_len * sizeof(int));
	cutilSafeCall(cudaMemset(d_depthArr + read_len, 0, win_size * sizeof(int)));

	GPUTOGPU(d_repeat_timeArr, d_repeat_timeArr + win_size, read_len * sizeof(int));
	cutilSafeCall(cudaMemset(d_repeat_timeArr + read_len, 0, win_size * sizeof(int)));

	GPUTOGPU(d_dep_uniArr, d_dep_uniArr + win_size, read_len * sizeof(int));
	cutilSafeCall(cudaMemset(d_dep_uniArr + read_len, 0, win_size * sizeof(int)));

	GPUTOGPU(d_numNonZeroPerSite, d_numNonZeroPerSite + win_size, read_len * sizeof(uint32));
	cutilSafeCall(cudaMemset(d_numNonZeroPerSite + read_len, 0, win_size * sizeof(uint32)));

	GPUTOGPU(d_count_uniArr, d_count_uniArr + win_size * 4, read_len * 4 * sizeof(int));
	cutilSafeCall(cudaMemset(d_count_uniArr + read_len * 4, 0, win_size * 4 * sizeof(int)));

	GPUTOGPU(d_q_sumArr, d_q_sumArr + win_size * 4, read_len * 4 * sizeof(int));
	cutilSafeCall(cudaMemset(d_q_sumArr + read_len * 4, 0, win_size * 4 * sizeof(int)));

	GPUTOGPU(d_count_allArr, d_count_allArr + win_size * 4, read_len * 4 * sizeof(int));
	cutilSafeCall(cudaMemset(d_count_allArr + read_len * 4, 0, win_size * 4 * sizeof(int)));

	GPUTOGPU(d_sparseBaseInfo, d_sparseBaseInfo + win_size * MAX_NUM_NON_ZERO, read_len * MAX_NUM_NON_ZERO * sizeof(uint32));
	cutilSafeCall(cudaMemset(d_sparseBaseInfo + read_len * MAX_NUM_NON_ZERO, BASEINFO_MAX, win_size * MAX_NUM_NON_ZERO * sizeof(uint32)));

//	recycle_kernel<<<numBlock, numThread>>>(d_depthArr, d_repeat_timeArr, d_dep_uniArr, d_count_uniArr, d_q_sumArr, d_count_allArr, d_numNonZeroPerSite, d_sparseBaseInfo, win_size, read_len, global_win_size, MAX_NUM_NON_ZERO);
	cudaThreadSynchronize();
	cutilCheckMsg("cuda_recycle!!!!");
/*
        memcpy(posArr, posArr + win_size, read_len * sizeof(int));
	for(int i = read_len; i < read_len + win_size; i++){
		posArr[i] = posArr[i - 1] + 1;
	}
*/	
}

/*
__global__
void quick_recycle(
        int* d_depthArr,
        int* d_repeat_timeArr,
        int* d_dep_uniArr,
        int* d_count_uniArr,
        int* d_q_sumArr,
        int* d_count_allArr,
        uint32* d_numNonZeroPerSite,
        uint32* d_sparseBaseInfo,
        ubit64_t win_size,
        ubit64_t read_len,
        int global_win_size,
        uint32 MAX_NUM_NON_ZERO
){
        const uint32 numTotalThread = NUM_TOTAL_THREAD;
        const uint32 globalThreadOffset = GLOBAL_THREAD_OFFSET;

        for(int i = globalThreadOffset; i < read_len; i += numTotalThread){
                d_depthArr[i] = d_depthArr[i + win_size];
                d_repeat_timeArr[i] = d_repeat_timeArr[i + win_size];
                d_dep_uniArr[i] = d_dep_uniArr[i + win_size];
                d_numNonZeroPerSite[i] = d_numNonZeroPerSite[i + win_size];

                for(int j = 0; j < 4; j++){
                        d_count_uniArr[i * 4 + j] = d_count_uniArr[(i + win_size) * 4 + j];
                        d_q_sumArr[i * 4 + j] = d_q_sumArr[(i + win_size) * 4 + j];
                        d_count_allArr[i * 4 + j] = d_count_allArr[(i + win_size) * 4 + j];
                }


                for(uint32 offset = 0; offset < MAX_NUM_NON_ZERO; offset++){
                        d_sparseBaseInfo[BASEINFO_IDX(i, offset, win_size + read_len, MAX_NUM_NON_ZERO)] =
                                d_sparseBaseInfo[BASEINFO_IDX(i + global_win_size, offset, win_size + read_len, MAX_NUM_NON_ZERO)];
                }
        }

        for(int i = globalThreadOffset + read_len; i < read_len + win_size; i += numTotalThread){
                d_depthArr[i] = 0;
                d_repeat_timeArr[i] = 0;
                d_dep_uniArr[i] = 0;
                d_numNonZeroPerSite[i] = 0;

                for(int j = 0; j < 4; j++){
                        d_count_uniArr[i * 4 + j] = 0;
                        d_q_sumArr[i * 4 + j] = 0;
                        d_count_allArr[i * 4 + j] = 0;
                }
                for(uint32 offset = 0; offset < MAX_NUM_NON_ZERO; offset++) {
                        d_sparseBaseInfo[BASEINFO_IDX(i, offset, numSite, MAX_NUM_NON_ZERO)] = BASEINFO_MAX;
                }

	
        }

}
*/

void cuda_quickRecycle(
        ubit64_t win_size,
        ubit64_t read_len,
	int* posArr,

        int* count_uniArr,
        int* q_sumArr,
        int* count_allArr,
        int* depthArr,
        int* repeat_timeArr,

	const uint32 numBlock, const uint32 numThread
)
{
        GPUTOGPU(d_depthArr, d_depthArr + win_size, read_len * sizeof(int));
        cutilSafeCall(cudaMemset(d_depthArr + read_len, 0, win_size * sizeof(int)));

        GPUTOGPU(d_repeat_timeArr, d_repeat_timeArr + win_size, read_len * sizeof(int));
        cutilSafeCall(cudaMemset(d_repeat_timeArr + read_len, 0, win_size * sizeof(int)));

        GPUTOGPU(d_dep_uniArr, d_dep_uniArr + win_size, read_len * sizeof(int));
        cutilSafeCall(cudaMemset(d_dep_uniArr + read_len, 0, win_size * sizeof(int)));

        GPUTOGPU(d_numNonZeroPerSite, d_numNonZeroPerSite + win_size, read_len * sizeof(uint32));
        cutilSafeCall(cudaMemset(d_numNonZeroPerSite + read_len, 0, win_size * sizeof(uint32)));

        GPUTOGPU(d_count_uniArr, d_count_uniArr + win_size * 4, read_len * 4 * sizeof(int));
        cutilSafeCall(cudaMemset(d_count_uniArr + read_len * 4, 0, win_size * 4 * sizeof(int)));

        GPUTOGPU(d_q_sumArr, d_q_sumArr + win_size * 4, read_len * 4 * sizeof(int));
        cutilSafeCall(cudaMemset(d_q_sumArr + read_len * 4, 0, win_size * 4 * sizeof(int)));

        GPUTOGPU(d_count_allArr, d_count_allArr + win_size * 4, read_len * 4 * sizeof(int));
        cutilSafeCall(cudaMemset(d_count_allArr + read_len * 4, 0, win_size * 4 * sizeof(int)));

        GPUTOGPU(d_sparseBaseInfo, d_sparseBaseInfo + win_size * MAX_NUM_NON_ZERO, read_len * MAX_NUM_NON_ZERO * sizeof(uint32));
        cutilSafeCall(cudaMemset(d_sparseBaseInfo + read_len * MAX_NUM_NON_ZERO, BASEINFO_MAX, win_size * MAX_NUM_NON_ZERO * sizeof(uint32)));

	cudaThreadSynchronize();
	cutilCheckMsg("cuda_quickRecycle!!!!!!!!!!!!!!");

	memcpy(posArr, posArr + win_size, read_len * sizeof(int));
        for(int i = read_len; i < read_len + win_size; i++){
                posArr[i] = posArr[i - 1] + 1;
        }

	const uint32 _counting_length = global_win_size + read_len;
        FROMGPU(count_uniArr, d_count_uniArr, sizeof(int) * COUNTUNI_SIZE * _counting_length);
        FROMGPU(q_sumArr, d_q_sumArr, sizeof(int) * QSUM_SIZE * _counting_length);
        FROMGPU(count_allArr, d_count_allArr, sizeof(int) * COUNTALL_SIZE * _counting_length);
        FROMGPU(depthArr, d_depthArr, sizeof(int) * _counting_length);
        FROMGPU(repeat_timeArr, d_repeat_timeArr, sizeof(int) * _counting_length);



}
#endif
