#ifndef UPDATE_CU_
#define UPDATE_CU_

#include "common_kernel.h"
#include <thrust/extrema.h>



char d_qmin = 0;
int inputSizeChar = 0;
int inputSizeInt = 0;


double memcpyTimeUpdate = 0.0;
double kernelTimeUpdate = 0.0;
double kernelTimeMaxNum = 0.0;
double matchTime = 0.0;
double totalFunTime = 0.0;
void gold_pos_match(int* positionArr, int* siteMatchArr, int readCount)
{

	int siteID = 0;
	memset(siteMatchArr, 0, sizeof(int)*global_win_size * 2);
	for(int readId = 0; readId < readCount; readId++){
		siteID = positionArr[readId]%global_win_size;
		if(siteMatchArr[siteID*2 + 1] == 0){
			siteMatchArr[siteID*2] = readId;
			siteMatchArr[siteID*2 + 1] = 1;
		}else{
			siteMatchArr[siteID*2 + 1] += 1;
		}
	}

}


void updateKernelInit(const uint32 winSize, const uint32 _read_length) {

	const uint32 _counting_length = winSize + _read_length;
        inputSizeChar = sizeof(char)* _read_length * _read_num_per_window
                + sizeof(char) * _read_length * _read_num_per_window
                + sizeof(char) * _read_num_per_window;
        inputSizeInt = sizeof(int) * _read_num_per_window
                + sizeof(int) * global_win_size * 2;

        GPUMALLOC((void**)&d_inputUpdateChar, inputSizeChar);
        GPUMALLOC((void**)&d_inputUpdateInt, inputSizeInt);

        d_readArr = d_inputUpdateChar;
        d_qualArr = d_inputUpdateChar + _read_length * _read_num_per_window;
        d_strandArr = d_inputUpdateChar + _read_length * _read_num_per_window
                                        + _read_length * _read_num_per_window;

	d_readLen_hitArr = d_inputUpdateInt;
        d_siteMatchArr = d_inputUpdateInt + _read_num_per_window;

	GPUMALLOC((void**)&d_sparseBaseInfo, sizeof(uint32) * _counting_length  * _ini_max_num_non_zero_per_site);

	GPUMALLOC((void**)&d_numNonZeroPerSite, sizeof(uint32) * _counting_length);

	GPUMALLOC((void**)&d_count_uniArr, sizeof(int) * COUNTUNI_SIZE * _counting_length);

	GPUMALLOC((void**)&d_q_sumArr, sizeof(int) * QSUM_SIZE * _counting_length);

	GPUMALLOC((void**)&d_count_allArr, sizeof(int) * COUNTALL_SIZE * _counting_length);

	GPUMALLOC((void**)&d_depthArr, sizeof(int) * _counting_length);
	
	GPUMALLOC((void**)&d_positionArr, sizeof(int) * winSize);

	GPUMALLOC((void**)&d_repeat_timeArr, sizeof(int) * _counting_length);

	GPUMALLOC((void**)&d_dep_uniArr, sizeof(int) * _counting_length);

	GPUMALLOC((void**)&d_maxNumNonZero, sizeof(uint32));
	cutilSafeCall(cudaMemset(d_maxNumNonZero, 0, sizeof(uint32)));
}


void updateKernelFinalize() {
	GPUFREE(d_inputUpdateChar);
	GPUFREE(d_inputUpdateInt);
	GPUFREE(d_sparseBaseInfo);
	GPUFREE(d_numNonZeroPerSite);
	GPUFREE(d_count_uniArr);
	GPUFREE(d_q_sumArr);
	GPUFREE(d_count_allArr);
	GPUFREE(d_depthArr);
	GPUFREE(d_repeat_timeArr);
	GPUFREE(d_dep_uniArr);
	GPUFREE(d_maxNumNonZero);
	GPUFREE(d_positionArr);
}


int gold_updateSparse(
	struct st_update* updateData,
	uint32** sparse_base_infoArr,
	uint32** numNonZeroPerSite,
	int maxNumNonZero,
	int read_len,

	int* count_uniArr,
	int* q_sumArr,
	int* count_allArr,

	int* depthArr,
	int* repeat_timeArr,
	int* dep_uniArr,
	int readCount
)
{


   const uint32 _read_length = read_len;
    int readNum = readCount;
    int supposeMaxNumNonZero = _ini_max_num_non_zero_per_site;	
    maxNumNonZero = 0;




    uint64_t coord, sub;
    for(int id=0; id < readNum; id++){
        for (coord = 0; coord < updateData->read_lenArr[id]; coord++) {
                if ((updateData->positionArr[id] + coord) / global_win_size == updateData->positionArr[id] / global_win_size) {
                        sub = (updateData->positionArr[id] + coord) % global_win_size;
                } else {
                        sub = (updateData->positionArr[id] + coord) % global_win_size + global_win_size; // Use the tail to store the info so that it won't intervene the uncalled bases
                }
                depthArr[sub] += 1;
                repeat_timeArr[sub] += updateData->hitArr[id];
                if ((updateData->readArr[id * _read_length + coord] == 'N') || updateData->qualArr[id * _read_length + coord] < updateData->qmin
                                || dep_uniArr[sub] >= 0xFF) {// An N, low quality or meaningless huge depth
                        continue;
                }
                if (updateData->hitArr[id] == 1) {
                        dep_uniArr[sub] += 1;// Update the covering info: 4x2x64x64 matrix, base x strand x q_score x read_pos, 2-1-6-6 bits for each
                        if (updateData->strandArr[id] == '+') {// Binary strand: 0 for plus and 1 for minus
		
                                const uint32 prefix_word = (uint32)(
                                                            (((updateData->readArr[id * _read_length + coord] & 0x6)>>1)) << 15
                                                            |(MAX_QSCORE - (updateData->qualArr[id * _read_length + coord] - updateData->qmin))<<9
                                                            |coord << 1
                                                            |0
                                                           );
	
				(*numNonZeroPerSite)[sub] += 1;
      
                                if((*numNonZeroPerSite)[sub] > maxNumNonZero) {
                                	maxNumNonZero = (*numNonZeroPerSite)[sub];
                                        assert(maxNumNonZero <= _ini_max_num_non_zero_per_site);
                                        if(maxNumNonZero > supposeMaxNumNonZero){
                                		ERROR_EXIT("exceeds the pre-allcoated memory, fix it later.");       
					}				
                                }

				small_int base = 1;

				(*sparse_base_infoArr)[BASEINFO_IDX(sub, ((*numNonZeroPerSite)[sub] - 1), global_num_site, MAX_NUM_NON_ZERO)] = prefix_word<<8|base;
			
                        } else {

                                const uint32 prefix_word = (uint32)(
                                                            (((updateData->readArr[id * _read_length + coord] & 0x6)>>1)) << 15
                                                            |(MAX_QSCORE - (updateData->qualArr[id * _read_length + coord] - updateData->qmin))<<9
                                                            |(updateData->read_lenArr[id] - 1 - coord)<<1
                                                            | 1
                                                           );

                        	(*numNonZeroPerSite)[sub] += 1;
				
                                if((*numNonZeroPerSite)[sub] > maxNumNonZero) {
                                        maxNumNonZero = (*numNonZeroPerSite)[sub];
                                        assert(maxNumNonZero <= _ini_max_num_non_zero_per_site);
                                        if(maxNumNonZero > supposeMaxNumNonZero){
                                                ERROR_EXIT("exceeds the pre-allcoated memory, fix it later.");                                        }                               
                                }

                                small_int base = 1;
                                (*sparse_base_infoArr)[BASEINFO_IDX(sub, ((*numNonZeroPerSite)[sub] - 1), global_num_site, MAX_NUM_NON_ZERO)] = prefix_word<<8|base;
			
                        }
                        count_uniArr[sub*COUNTUNI_SIZE +( (updateData->readArr[id * _read_length + coord] >> 1) & 3)] += 1;
                      	q_sumArr[sub*QSUM_SIZE +( (updateData->readArr[id * _read_length + coord] >> 1) & 3)]
                                        += (updateData->qualArr[id * _read_length + coord] - updateData->qmin);
                } else {
                        ;// Repeats
                }

                count_allArr[sub*COUNTALL_SIZE +( (updateData->readArr[id * _read_length + coord] >> 1) & 3)] += 1;
        }
    }
    return maxNumNonZero;
}




bool sparse_base_info_check(
                        uint32* sparse_new,
                        uint32* numPerSite_new,
                        uint32* sparse_base,
                        uint32* numPerSite_base,
                        uint32  numSite
){
	bool is_the_same = true;
	uint64 sum_new = 0, sum_base = 0;
	for(uint32 siteID = 0; siteID < numSite; siteID++){
		if(numPerSite_new[siteID] != numPerSite_base[siteID]){
			printf("numPerSite diff at site: %d [%d/%d]\n ",siteID,numPerSite_new[siteID], numPerSite_base[siteID]);
			is_the_same = false;
		}
	}
	if(is_the_same){
		for(uint32 siteID = 0; siteID < numSite; siteID++){
			for(uint32 offset =0; offset < numPerSite_base[siteID]; offset++)
			{
				sum_new += sparse_new[offset*numSite + siteID];
				sum_base += sparse_base[offset*numSite +siteID];
			}
			if(sum_new != sum_base){
				printf("sparse diff at site: %d [%d/%d]\n",siteID, sum_new, sum_base);
				is_the_same = false;
			}
		}
	}
	return is_the_same;	
}



__global__
void updateSparse_kernel(
		int* d_positionArr,
                char* d_readArr,
                char* d_qualArr,
		int* d_readLen_hitArr,
		int* d_siteMatchArr,
               	char* d_strandArr,
	        uint32* d_sparseBaseInfo,
        	uint32* d_numNonZeroPerSite,
              	int* d_count_uniArr,
		int* d_q_sumArr,
              	int* d_count_allArr,
                int* d_depthArr,
                int* d_repeat_timeArr,
                int* d_dep_uniArr,
		uint32* d_maxNumNonZero,

		int  _counting_length,
		int  global_win_size,
		int  _read_length,
		int readCount,
		char qmin, 
		const uint32 offsetPitch)
{
    	const uint32 numTotalThread = NUM_TOTAL_THREAD;
    	const uint32 globalThreadOffset = GLOBAL_THREAD_OFFSET;
    	for(int siteID = globalThreadOffset; siteID < global_win_size + _read_length; siteID += numTotalThread){
		int coord = 0, strand_val = 0, coord_val = 0;
		int depth = 0, repeat_time = 0, dep_uni = 0, num = 0;
		

		if(siteID < _read_length){
			dep_uni = d_dep_uniArr[siteID];
			num = d_numNonZeroPerSite[siteID];
		}
		int siteIdCov = siteID - _read_length + 1;
		if(siteIdCov >= global_win_size){break;}
		if(siteIdCov < 0){siteIdCov = 0;}

		for(; siteIdCov <= siteID; siteIdCov++ ){
			if(siteIdCov >= global_win_size){break;}
			int siteMatch = d_siteMatchArr[siteIdCov * 2];
			int siteMatchMax = siteMatch + d_siteMatchArr[siteIdCov * 2 + 1];
			for(int readId = siteMatch; readId < siteMatchMax; readId++){
				coord = siteID - siteIdCov;	
				int readLen_hit = d_readLen_hitArr[readId];
				int read_len = (readLen_hit&0xffff0000)>>16;
				if((coord >= read_len) || (coord < 0)){continue;}
				int hit = (readLen_hit&0x0000ffff); 
				char read = d_readArr[readId * _read_length + coord];
                                char qual = d_qualArr[readId * _read_length + coord];

                		depth += 1;
				repeat_time += hit;
                		if((read == 'N') || (qual < qmin)|| (dep_uni >= 0xFF)) {// An N, low quality or meaningless huge depth
                        		continue;
                		}
                		if (hit == 1) {
                        		dep_uni += 1;
                                        if (d_strandArr[readId] == '+') {// Binary strand: 0 for plus and 1 for minus
                                                strand_val = 0;
						coord_val = coord;
                                        }
                                        else{ 
						strand_val = 1;
						coord_val = read_len - 1 - coord;
       					}

					const uint32 prefix_word = (uint32)( 
							    	(((read & 0x6)>>1)) << 15
							    	|(MAX_QSCORE - (qual - qmin))<<9
							    	|coord_val << 1
							    	|strand_val
							   	);
					num += 1;
					small_int base = 1;
					d_sparseBaseInfo[siteID * offsetPitch + num - 1] = prefix_word<<8|base;
					d_count_uniArr[siteID*COUNTUNI_SIZE + ((read >> 1) & 3)] += 1;
                     			d_q_sumArr[siteID * QSUM_SIZE + ((read >> 1) & 3)] += (qual - qmin);
                		} else {
                        		;// Repeats
                		}
                		d_count_allArr[siteID*COUNTALL_SIZE + ((read >> 1) & 3)] += 1;
			}
        	}
//		if(d_positionArr[siteID] == 203685889){printf("BEGIN: depth = %d, d_depth = %d", depth, d_depthArr[siteID]);}
		d_depthArr[siteID] += depth;
//		if(d_positionArr[siteID] == 203685889){printf("END: depth = %d, d_depth = %d", depth, d_depthArr[siteID]);}
		d_repeat_timeArr[siteID] += repeat_time;
	 	d_dep_uniArr[siteID] = dep_uni;
		d_numNonZeroPerSite[siteID] = num;
	}  
}

uint32 cuda_updateSparse(
        struct st_update* updateData,
        uint32* sparseBaseInfo,
        uint32* numNonZeroPerSite,
        int read_len,

        int* count_uniArr,
        int* q_sumArr,
        int* count_allArr,
        int* depthArr,
        int* repeat_timeArr,
        int* dep_uniArr,
	int readCount, 
	const uint32 numBlock, 
	const uint32 numThread) {


	const uint32 _counting_length = global_win_size + read_len;
        TOGPU(d_inputUpdateChar, updateData->inputUpdateChar, inputSizeChar);
        TOGPU(d_inputUpdateInt, updateData->inputUpdateInt, inputSizeInt);
	TOGPU(d_positionArr, updateData->positionArr, global_win_size * sizeof(int));
	d_qmin = updateData->qmin;

	updateSparse_kernel<<<numBlock, numThread>>>(d_positionArr, d_readArr, d_qualArr, d_readLen_hitArr, d_siteMatchArr, d_strandArr, d_sparseBaseInfo, d_numNonZeroPerSite, d_count_uniArr, d_q_sumArr, d_count_allArr, d_depthArr, d_repeat_timeArr, d_dep_uniArr, d_maxNumNonZero,  _counting_length, global_win_size, read_len, readCount, d_qmin, MAX_NUM_NON_ZERO);

	cudaThreadSynchronize();

        thrust::device_ptr<uint32> numPtr(d_numNonZeroPerSite);
        thrust::device_ptr<uint32> maxPtr(d_maxNumNonZero);
        maxPtr = thrust::max_element(numPtr, numPtr + global_win_size);
        d_maxNumNonZero = thrust::raw_pointer_cast(maxPtr);


	cudaThreadSynchronize();
	cutilCheckMsg("updateSparse_kernel");
	
	uint32* tmpMaxNum = (uint32*)malloc(sizeof(uint32));
	FROMGPU(tmpMaxNum, d_maxNumNonZero, sizeof(uint32));
	uint32 maxNumNonZero = tmpMaxNum[0];
	free(tmpMaxNum);
	FROMGPU(count_uniArr, d_count_uniArr, sizeof(int) * COUNTUNI_SIZE * _counting_length);
	FROMGPU(q_sumArr, d_q_sumArr, sizeof(int) * QSUM_SIZE * _counting_length);
	FROMGPU(count_allArr, d_count_allArr, sizeof(int) * COUNTALL_SIZE * _counting_length);
	FROMGPU(depthArr, d_depthArr, sizeof(int) * _counting_length);
	FROMGPU(repeat_timeArr, d_repeat_timeArr, sizeof(int) * _counting_length);

	return maxNumNonZero;
}


#endif
