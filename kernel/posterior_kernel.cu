#ifndef POSTERIOR_KERNEL_
#define POSTERIOR_KERNEL_
#include "gsnp.h"
#include "common_kernel.h"
#include "cuda_header.h"
#include "likelihood_kernel.cu"
#include "update_kernel.cu"
//#include "bitonic_sort.cu"
#include "common.h"
#include <string>
double memcpyTimePost = 0.0;
double memcpySizePost = 0.0;
double kernelTimePost = 0.0;
int* d_debug = NULL;
rate_t* d_type_prob = NULL;
char* d_type2 = NULL;
int* d_qual1 = NULL;
int* d_qual2 = NULL;
int* d_qual3 = NULL;
void posteriorKernelInit(
        const rate_t* p_prior,
        const rate_t* p_rank,
	int win_size,
	int read_len)
{
	int numSite = win_size + read_len;
	GPUMALLOC((void**)&d_p_prior, sizeof(rate_t) * 8 * 4 * 4);
	GPUMALLOC((void**)&d_p_rank, sizeof(rate_t) * 64 * 64 * 2048);
	GPUMALLOC((void**)&d_oriArr, sizeof(char) * numSite);

	GPUMALLOC((void**)&d_base1, sizeof(char) * win_size);
	cutilSafeCall(cudaMemset(d_base1, 0, sizeof(char) * win_size));
	GPUMALLOC((void**)&d_base2, sizeof(char) * win_size);
	cutilSafeCall(cudaMemset(d_base2, 0, sizeof(char) * win_size));
	GPUMALLOC((void**)&d_type1, sizeof(char) * win_size);
	cutilSafeCall(cudaMemset(d_type1, 0, sizeof(char) * win_size));
	GPUMALLOC((void**)&d_q_cns, sizeof(int) * win_size);
	cutilSafeCall(cudaMemset(d_q_cns, 0, sizeof(int) * win_size));
	GPUMALLOC((void**)&d_rank_sum_test_value, sizeof(double) * win_size);
	cutilSafeCall(cudaMemset(d_rank_sum_test_value, 0, sizeof(double) * win_size));
	GPUMALLOC((void**)&d_validatedArr, sizeof(bool) * win_size);
	GPUMALLOC((void**)&d_hapmap_siteArr, sizeof(bool) * win_size);
	GPUMALLOC((void**)&d_indel_siteArr, sizeof(bool) * win_size);
	GPUMALLOC((void**)&d_freqArr, sizeof(rate_t) * 4 * win_size);
	GPUMALLOC((void**)&d_nflag, sizeof(bool) * win_size);


        GPUMALLOC((void**)&d_debug, sizeof(int) * 20);
        cutilSafeCall(cudaMemset(d_debug, 0, sizeof(int) * 20));

	GPUMALLOC((void**)&d_type_prob, sizeof(rate_t) * 17 * win_size);
	cutilSafeCall(cudaMemset(d_type_prob, 0.0, sizeof(rate_t) * 17 * win_size));
	GPUMALLOC((void**)&d_type2, sizeof(char) * win_size);
	cutilSafeCall(cudaMemset(d_type2, 0, sizeof(char) * win_size));
	GPUMALLOC((void**)&d_qual1, sizeof(int) * win_size);
	cutilSafeCall(cudaMemset(d_qual1, 0, sizeof(int) * win_size));
	GPUMALLOC((void**)&d_qual2, sizeof(int) * win_size);
	cutilSafeCall(cudaMemset(d_qual2, 0, sizeof(int) * win_size));
	GPUMALLOC((void**)&d_qual3, sizeof(int) * win_size);
	cutilSafeCall(cudaMemset(d_qual3, 0, sizeof(int) * win_size));

        TOGPU(d_p_rank, p_rank, sizeof(double) * 64 * 64 * 2048);
        TOGPU(d_p_prior, p_prior, sizeof(rate_t) * 8 * 4 * 4);

}
void posteriorKernelFinalize()
{
	GPUFREE(d_p_prior);
	GPUFREE(d_p_rank);
	GPUFREE(d_oriArr);
	GPUFREE(d_base1);
	GPUFREE(d_base2);
	GPUFREE(d_type1);
	GPUFREE(d_q_cns);
	GPUFREE(d_rank_sum_test_value);	
	GPUFREE(d_validatedArr);
	GPUFREE(d_hapmap_siteArr);
	GPUFREE(d_indel_siteArr);
	GPUFREE(d_freqArr);
	GPUFREE(d_nflag);

	GPUFREE(d_debug);
	GPUFREE(d_type_prob);
	GPUFREE(d_type2);
	GPUFREE(d_qual1);
	GPUFREE(d_qual2);
	GPUFREE(d_qual3);
}

double gold_table_test(double *p_rank, int n1, int n2, double T1, double T2) {
        if(n1<=n2) {
                return p_rank[(n1+n2)<<17|n1<<11|(int)(T1)]+(T1-(int)T1)*(p_rank[(n1+n2)<<16|n1<<11|(int)(T1+1)]-p_rank[(n1+n2)<<17|n1<<11|(int)(T1)]);
        }
        else {
                return p_rank[(n1+n2)<<17|n2<<11|(int)(T2)]+(T2-(int)T2)*(p_rank[(n1+n2)<<16|n2<<11|(int)(T2+1)]-p_rank[(n1+n2)<<17|n2<<11|(int)(T2)]);
        }
}
double gold_normal_value(double z) {
        if (z>6.0 || z<-6.0) {
                return 0.0;
        }
        else {
                double b1 = 0.31938153;
                double b2 = -0.356563782;
                double b3 = 1.781477937;
                double b4 = -1.821255978;
                double b5 = 1.330274429;
                double p = 0.2316419;
                double c2 = 0.39894228;

                double a = fabs(z);
                double t = 1.0/(1.0+a*p);
                double b = c2*exp((-z)*(z/2.0));
                double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
                n = 1.0 - b*n;
                if (z < 0.0) n = 1.0 - n;
                return n>0.5?1-n:n;
        }
}
double gold_normal_test(int n1, int n2, double T1, double T2) {
        double u1, u2;
        u1 = (T1 - n1*(n1+n2+1)/2) / sqrt(n1*n2*(n1+n2+1)/(double)12);
        u2 = (T2 - n2*(n1+n2+1)/2) / sqrt(n1*n2*(n1+n2+1)/(double)12);
        return gold_normal_value(fabs(u1)>fabs(u2)?u1:u2);
}

double gold_rank_test(
        const uint32* sparse_base_infoArr,
        const uint64 siteID,
        uint32 numNonZero,
        uint64 read_len,
        int* count_uni_buf,
        int q_max_min,
        char best_type,
        double* p_rank)
{
        if( (best_type&3) == ((best_type>>2)&3) ) {// HOM
                return 1.0;
        }
        if( count_uni_buf[siteID*COUNTUNI_SIZE + (best_type&3)]==0 || count_uni_buf[siteID*COUNTUNI_SIZE + ((best_type>>2)&3)]==0) {// HET with one allele...
                return 0.0;
        }
        int *same_qual_count = new int [64];
        double *rank_array = new double [64];
        memset(same_qual_count,0,sizeof(int)*64);
        memset(rank_array,0,sizeof(double)*64);

        int rank(0);
        double T[4]={0.0, 0.0, 0.0, 0.0};
        bool is_need[4] ={false,false,false,false};
        is_need[(best_type&3)]=true; is_need[((best_type>>2)&3)]=true;
        int  q_score;

       // uint64 numSite = global_win_size + read_len;
        for(uint32 offset = 0; offset < numNonZero; offset++)
        {

                //const uint32 word = sparse_base_infoArr[offset*numSite + siteID];

                const uint32 word = sparse_base_infoArr[BASEINFO_IDX(siteID, offset, numSite, MAX_NUM_NON_ZERO)];
                const int o_base = (word&0x01800000)>>23;
                const int q_score = MAX_QSCORE - ((word&0x007e0000)>>17);
                const int base = (word&0x000000ff);
                if(count_uni_buf[siteID*COUNTUNI_SIZE+o_base] == 0 || !is_need[o_base])continue;
                same_qual_count[q_score] += base;

        }
        rank = 0;
        for(q_score=0;q_score<=(ubit64_t)(q_max_min + 1);q_score++) {
                rank_array[q_score]= rank+(1+same_qual_count[q_score])/2.0;
                rank += same_qual_count[q_score];
        }

        for(uint32 offset =0; offset < numNonZero; offset++)
        {
                //uint32 word = sparse_base_infoArr[offset*numSite + siteID];
                uint32 word = sparse_base_infoArr[BASEINFO_IDX(siteID, offset, numSite, MAX_NUM_NON_ZERO)];
                int o_base = (word&0x01800000)>>23;
                int q_score = MAX_QSCORE - ((word&0x007e0000)>>17);
                int base = (word&0x000000ff);
                if(count_uni_buf[siteID*COUNTUNI_SIZE+o_base] == 0 || !is_need[o_base])continue;
                T[o_base] += rank_array[q_score] * base;
		if(siteID == 14429 && offset == 45){cerr<<word<<"#"<<BASEINFO_IDX(siteID, offset, numSite, MAX_NUM_NON_ZERO)<<"#"<<base<<"#"<<T[o_base]<<"#"<<rank_array[q_score]<<"#"<<same_qual_count[q_score]<<endl;}

        }

        delete [] same_qual_count;
        delete [] rank_array;

        if (count_uni_buf[siteID*COUNTUNI_SIZE + (best_type&3)]+count_uni_buf[siteID*COUNTUNI_SIZE + ((best_type>>2)&3)]<64) {
                return gold_table_test(p_rank, count_uni_buf[siteID*COUNTUNI_SIZE + (best_type&3)], count_uni_buf[siteID*COUNTUNI_SIZE + ((best_type>>2)&3)], T[best_type&3], T[(best_type>>2)&3]);
        }
        else {
                return gold_normal_test(count_uni_buf[siteID*COUNTUNI_SIZE + (best_type&3)], count_uni_buf[siteID*COUNTUNI_SIZE + ((best_type>>2)&3)],T[best_type&3], T[(best_type>>2)&3]);
        }

}
int gold_snp_p_prior_gen(double * real_p_prior, char ref, 
				bool is_indel, bool is_validated, bool is_hapmap, rate_t* freq,
				rate_t para_het_val_r, rate_t para_het_unval_r, rate_t para_althom_val_r, rate_t para_althom_unval_r 
				) {
        if (is_indel) {
                return 0;
        }
        char base, allele1, allele2;
        int allele_count;
        allele_count = 0;
        for (base=0; base != 4; base ++) {
                if(freq[base]>0) {
                        // The base is found in dbSNP
                        allele_count += 1;
                }
        }
        if(allele_count <= 1) {
                // Should never occur
                cerr<<"Previous Extract SNP error."<<endl;
                exit(255);
                return -1;
        }
        char t_base = (ref&0x3);
        for(allele1=0;allele1!=4;allele1++) {
                for(allele2=allele1;allele2!=4;allele2++) {
                        if (! is_hapmap) {
                                if(freq[allele1]>0 && freq[allele2]>0) {
                                        // Here the frequency is just a tag to indicate SNP alleles in non-HapMap sites
                                        if(allele1 == allele2 && allele1 == t_base) {
                                                // refHOM
                                                real_p_prior[allele1<<2|allele2] = 1;
                                        }
                                        else if (allele1 == t_base || allele2 == t_base) {
                                                // refHET: 1 ref 1 alt
                                                real_p_prior[allele1<<2|allele2] = is_validated ? para_het_val_r : para_het_unval_r;
                                        }
                                        else if (allele1 == allele2) {
                                                real_p_prior[allele1<<2|allele2] =  is_validated ? para_althom_val_r : para_althom_unval_r;
                                        }
                                        else {
                                                // altHET: 2 diff alt base
                                                real_p_prior[allele1<<2|allele2] = is_validated ? para_het_val_r : para_het_unval_r;
                                        }
                                }
                        }
                        else {
                                // Real HapMap Sites
                                if(freq[allele1]>0 && freq[allele2]>0) {
                                        real_p_prior[allele1<<2|allele2] = (allele1==allele2 ? 1 : (2*para_het_val_r)) * freq[allele1] * freq[allele2];
                                }
                        }
                }
        }
        return 1;
}

void gold_posterior(
	st_likelihood* likelihoodData,
	rate_t* p_prior,
	rate_t* type_prob,//(rate_t*)malloc(sizeof(rate_t) * LIKELY_SIZE)
	rate_t* p_rank,//(rate_t*)malloc(sizeof(rate_t) * PRIOR_SIZE)
	uint32* sparse_base_infoArr,
	uint32* numNonZeroPerSite,
        int* count_uni_buf,
	int* q_sumArr,
	int* dep_uniArr,
	int* count_allArr,
	char* oriArr,
	int* posArr,

	char* base1,
	char* base2,
	char* type1,
	int* q_cns,
	double* rank_sum_test_value,
	
	bool* validatedArr, bool* hapmap_siteArr, bool* indel_siteArr, rate_t* freqArr,
		
	rate_t para_het_val_r, rate_t para_het_unval_r, rate_t para_althom_val_r, rate_t para_althom_unval_r, bool para_rank_sum_mode, bool para_is_monoploid, bool para_refine_mode, char para_q_max_min,

	ubit64_t read_len,
	uint64_t call_length	
)
{
        char allele1, allele2, genotype, 
                        type2 /*, base1, base2, base3*/;
        int i, qual1, qual2, qual3, /*q_cns, */all_count1,
                        all_count2, all_count3;
	//char base3;
	double real_p_prior[16];	
        for (int siteID= 0; siteID != call_length; siteID++) {
                if(likelihoodData->nflag[siteID]) {
                        continue;
                }
	
                //adjust base1, base2, and base3
                base1[siteID] = 0, base2[siteID] = 0; //, base3 = 0;
                qual1 = -1, qual2 = -2, qual3 = -3;
                all_count1 = 0, all_count2 = 0, all_count3 = 0;
                if (dep_uniArr[siteID]) {
                        for (i = 0; i != 4; i++) {
                                if (q_sumArr[siteID * QSUM_SIZE + i] >= qual1) {
                                        //base3 = base2[siteID];
                                        qual3 = qual2;
                                        base2[siteID] = base1[siteID];
                                        qual2 = qual1;
                                        base1[siteID] = i;
                                        qual1 = q_sumArr[siteID * QSUM_SIZE + i];
                                } else if (q_sumArr[siteID * QSUM_SIZE + i] >= qual2) {
                                        //base3 = base2[siteID];
                                        qual3 = qual2;
                                        base2[siteID] = i;
                                        qual2 = q_sumArr[siteID * QSUM_SIZE + i];
                                } else if (q_sumArr[siteID * QSUM_SIZE + i] >= qual3) {
                                        //base3 = i;
                                        qual3 = q_sumArr[siteID * QSUM_SIZE + i];
                                } else {
                                        ;
                                }
                        }
                        if (qual1 == 0) {
                                base1[siteID] = (oriArr[siteID] & 7);
                        } else if (qual2 == 0 && base1[siteID] != (oriArr[siteID] & 7)) {
                                base2[siteID] = (oriArr[siteID] & 7);
                        } else {
                                ;
                        }
                } else {
                        for (i = 0; i != 4; i++) {
                                if (count_allArr[siteID * COUNTALL_SIZE + i] >= all_count1) {
                                        //base3 = base2[siteID];
                                        all_count3 = all_count2;
                                        base2[siteID] = base1[siteID];
                                        all_count2 = all_count1;
                                        base1[siteID] = i;
                                        all_count1 = count_allArr[siteID * COUNTALL_SIZE + i];
                                } else if (count_allArr[siteID * COUNTALL_SIZE + i] >= all_count2) {
                                        //base3 = base2[siteID];
                                        all_count3 = all_count2;
                                        base2[siteID] = i;
                                        all_count2 = count_allArr[siteID * COUNTALL_SIZE + i];
                                } else if (count_allArr[siteID * COUNTALL_SIZE + i] >= all_count3) {
                                        //base3 = i;
                                        all_count3 = count_allArr[siteID * COUNTALL_SIZE + i];
                                } else {
                                        ;
                                }
                        }
                        if (all_count1 == 0) {
                                base1[siteID] = (oriArr[siteID] & 7);
                        } else if (all_count2 == 0 && base1[siteID] != (oriArr[siteID] & 7)) {
                                base2[siteID] = (oriArr[siteID] & 7);
                        } else {
                                ;
                        }
                }

                //calculate posterier probability
                memcpy(real_p_prior, &p_prior[((ubit64_t) oriArr[siteID] & 0x7) << 4], sizeof(double) * 16);
                if ((oriArr[siteID] & 0x8) && para_refine_mode) {
			gold_snp_p_prior_gen(real_p_prior, oriArr[siteID], 
						indel_siteArr[siteID], validatedArr[siteID], hapmap_siteArr[siteID], &freqArr[siteID * 4], 
						para_het_val_r, para_het_unval_r, para_althom_val_r, para_althom_unval_r);
		}
                memset(type_prob, 0, sizeof(rate_t) * 17);
                type2 = type1[siteID] = 16;
               	for (allele1 = 0; allele1 != 4; allele1++) {
                       	for (allele2 = allele1; allele2 != 4; allele2++) {
                               	genotype = allele1 << 2 | allele2;
                               	if (para_is_monoploid && allele1 != allele2) {
                                       	continue;
                               	}
				type_prob[genotype] = likelihoodData->typeLikelihood[siteID * LIKELY_SIZE + genotype]+ log10(
                                               real_p_prior[genotype]);

                               	if (type_prob[genotype] >= type_prob[type1[siteID]] || type1[siteID]
                                                	== 16) {
                                       	type2 = type1[siteID];
                                       	type1[siteID] = genotype;
                               	} else if (type_prob[genotype] >= type_prob[type2]
                                               	|| type2 == 16) {
                                       	type2 = genotype;
                               	} else {
                                       	;
                               	}
                       	}
               	}
                //rank_sum test and something else
                if (para_rank_sum_mode) {
                       	rank_sum_test_value[siteID]
                       		= gold_rank_test(sparse_base_infoArr, siteID, numNonZeroPerSite[siteID], read_len,
                                              count_uni_buf, para_q_max_min, type1[siteID],  p_rank);
                } else {
                       	rank_sum_test_value[siteID] = 1.0;
                }
                if (rank_sum_test_value[siteID] == 0.0) {
                       	q_cns[siteID] = 0;
                } else {
                       	q_cns[siteID] = (int) (10 * (type_prob[type1[siteID]] - type_prob[type2])
                                       + 10 * log10(rank_sum_test_value[siteID]));
                }
                if ((type1[siteID] & 3) == ((type1[siteID] >> 2) & 3)) {
                       	if (qual1 > 0 && base1[siteID] != (type1[siteID] & 3)) {
                               	q_cns[siteID] = 0;
                       	} else if (q_cns[siteID] > qual1 - qual2) {
                               	q_cns[siteID] = qual1 - qual2;
                	} else {
                               	;
                       	}
                } else {
                       	if (q_sumArr[siteID * QSUM_SIZE + base1[siteID]] > 0 
				&& q_sumArr[siteID * QSUM_SIZE + base2[siteID]] > 0 
				&& type1[siteID]  == (base1[siteID] < base2[siteID] ? (base1[siteID] << 2 | base2[siteID]) : (base2[siteID] << 2
                                                       | base1[siteID]))) {
                        	if (q_cns[siteID] > qual2 - qual3) {
                                	q_cns[siteID] = qual2 - qual3;
                        	}
                       	} else {
                                	q_cns[siteID] = 0;
                        }
                }
                if (q_cns[siteID] > 99) {
                       	q_cns[siteID] = 99;
                }
                if (q_cns[siteID] < 0) {
                       	q_cns[siteID] = 0;
                }
	}
}


__inline__ __device__
double gpu_table_test(double *d_p_rank, int n1, int n2, double T1, double T2, int* d_debug) {
	int idx1 = (n1+n2)<<17|n1<<11|(int)(T1);
	int idx2 = (n1+n2)<<16|n1<<11|(int)(T1+1);
	int idx3 = (n1+n2)<<17|n2<<11|(int)(T2);
	int idx4 = (n1+n2)<<16|n2<<11|(int)(T2+1);
        if(n1<=n2) {
//		return d_p_rank[(n1+n2)<<17|n1<<11|(int)(T1)]+(T1-(int)T1)*(d_p_rank[(n1+n2)<<16|n1<<11|(int)(T1+1)]-d_p_rank[(n1+n2)<<17|n1<<11|(int)(T1)]);
		return d_p_rank[idx1] + (T1-(int)T1)*(d_p_rank[idx2] - d_p_rank[idx1]);
        }
        else {
//	        return d_p_rank[(n1+n2)<<17|n2<<11|(int)(T2)]+(T2-(int)T2)*(d_p_rank[(n1+n2)<<16|n2<<11|(int)(T2+1)]-d_p_rank[(n1+n2)<<17|n2<<11|(int)(T2)]);
		return d_p_rank[idx3] + (T2-(int)T2)*(d_p_rank[idx4] - d_p_rank[idx3]);
        }
}

__inline__ __device__
double gpu_normal_value(double z) {
        if (z>6.0 || z<-6.0) {
                return 0.0;
        }
        else {
                double b1 = 0.31938153;
                double b2 = -0.356563782;
                double b3 = 1.781477937;
                double b4 = -1.821255978;
                double b5 = 1.330274429;
                double p = 0.2316419;
                double c2 = 0.39894228;

                double a = fabs(z);
                double t = 1.0/(1.0+a*p);
                double b = c2*exp((-z)*(z/2.0));
                double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
                n = 1.0 - b*n;
                if (z < 0.0) n = 1.0 - n;
                return n>0.5?1-n:n;
        }
}
__inline__ __device__
double gpu_normal_test(int n1, int n2, double T1, double T2) {
        double u1, u2;
        u1 = (T1 - n1*(n1+n2+1)/2) / sqrt(n1*n2*(n1+n2+1)/(double)12);
        u2 = (T2 - n2*(n1+n2+1)/2) / sqrt(n1*n2*(n1+n2+1)/(double)12);
        return gpu_normal_value(fabs(u1)>fabs(u2)?u1:u2);
}
__inline__ __device__
uint32 gpuGetReduceSize(const uint32 n){
    

	return n;
   
/*
	 bool isFind = false;
        uint32 m = 0;
	uint32 pow2Table[POW2_TABLE_SIZE] = {
        	2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384
	};

        for(uint32 i = 0; i < POW2_TABLE_SIZE; i++) {
                if(pow2Table[i] >= n) {
                        isFind = true;
                        m = pow2Table[i];
                        break;
                }
        }
*/
/*
        if(!isFind) {
        //      ERROR_EXIT("input number is too large or small!");
        }
*/
  //      return m;

}


__inline__ __device__
double gpu_rank_test(
                const uint32* d_sparseBaseInfo,
                const uint64 siteID,
                uint32 numNonZero,
                int* d_count_uniArr,
                double* d_p_rank,
                char best_type,
                uint64 read_len,
                int q_max_min,
		uint32 maxNumNonZero,
		int* d_debug
)
{
        int same_qual_count[64];
        double rank_array[64];
        double T[4];
        bool is_need[4];
        int rank = 0;
        //uint64 numSite = global_win_size + read_len;
        for(int i = 0; i < 64; i++){same_qual_count[i] = 0; rank_array[i] = 0.0;}
        for(int i = 0; i < 4; i++){T[i] = 0.0; is_need[i] = false;}
        is_need[(best_type&3)]=true; is_need[((best_type>>2)&3)]=true;

        if( (best_type&3) == ((best_type>>2)&3) ) {// HOM
                return 1.0;
        }

        if( d_count_uniArr[siteID*COUNTUNI_SIZE + (best_type&3)]==0 || d_count_uniArr[siteID*COUNTUNI_SIZE + ((best_type>>2)&3)]==0) {
                return 0.0;
        }

        for(uint32 offset = 0; offset < numNonZero; offset++)
        {
		const uint32 word = d_sparseBaseInfo[BASEINFO_IDX(siteID, offset, numSite, gpuGetReduceSize(MAX_NUM_NON_ZERO))];
//                const uint32 word = d_sparseBaseInfo[offset*numSite + siteID];
                const int o_base = (word&0x01800000)>>23;
                const int q_score = MAX_QSCORE - ((word&0x007e0000)>>17);
                const int base = (word&0x000000ff);
                if(d_count_uniArr[siteID*COUNTUNI_SIZE+o_base] == 0 || !is_need[o_base])continue;
                same_qual_count[q_score] += base;
        }

       rank = 0;
       for(int q_score=0;q_score<=(ubit64_t)(q_max_min + 1);q_score++) {
                rank_array[q_score]= rank+(1+same_qual_count[q_score])/2.0;
                rank += same_qual_count[q_score];
       }

       for(uint32 offset =0; offset < numNonZero; offset++)
       {
                //uint32 word = d_sparseBaseInfo[offset*numSite + siteID];
		uint32 word = d_sparseBaseInfo[BASEINFO_IDX(siteID, offset, numSite, gpuGetReduceSize(MAX_NUM_NON_ZERO))];
                int o_base = (word&0x01800000)>>23;
                int q_score = MAX_QSCORE - ((word&0x007e0000)>>17);
                int base = (word&0x000000ff);
                if(d_count_uniArr[siteID*COUNTUNI_SIZE+o_base] == 0 || !is_need[o_base])continue;
                T[o_base] += rank_array[q_score] * base;
       }

        if (d_count_uniArr[siteID*COUNTUNI_SIZE + (best_type&3)]+d_count_uniArr[siteID*COUNTUNI_SIZE + ((best_type>>2)&3)]<64) {
                return gpu_table_test(d_p_rank, d_count_uniArr[siteID*COUNTUNI_SIZE + (best_type&3)], d_count_uniArr[siteID*COUNTUNI_SIZE + ((best_type>>2)&3)], T[best_type&3], T[(best_type>>2)&3], d_debug);
        }
        else {
                return gpu_normal_test(d_count_uniArr[siteID*COUNTUNI_SIZE + (best_type&3)], d_count_uniArr[siteID*COUNTUNI_SIZE + ((best_type>>2)&3)],T[best_type&3], T[(best_type>>2)&3]);
        }
}


__inline__ __device__
int gpu_snp_p_prior_gen(double * real_p_prior, char ref,
                                bool is_indel, bool is_validated, bool is_hapmap, rate_t* freq,
                                rate_t para_het_val_r, rate_t para_het_unval_r, rate_t para_althom_val_r, rate_t para_althom_unval_r
                                ) {
        if (is_indel) {
                return 0;
        }
        char base, allele1, allele2;
        int allele_count;
        allele_count = 0;
        for (base=0; base != 4; base ++) {
                if(freq[base]>0) {
                        // The base is found in dbSNP
                        allele_count += 1;
                }
        }
//
//        if(allele_count <= 1) {
  //              // Should never occur
    //            cerr<<"Previous Extract SNP error."<<endl;
      //          exit(255);
        //        return -1;
//        }
//
        char t_base = (ref&0x3);
        for(allele1=0;allele1!=4;allele1++) {
                for(allele2=allele1;allele2!=4;allele2++) {
                        if (! is_hapmap) {
                                if(freq[allele1]>0 && freq[allele2]>0) {
                                        // Here the frequency is just a tag to indicate SNP alleles in non-HapMap sites
                                        if(allele1 == allele2 && allele1 == t_base) {
                                                // refHOM
                                                real_p_prior[allele1<<2|allele2] = 1;
                                        }
                                        else if (allele1 == t_base || allele2 == t_base) {
                                                // refHET: 1 ref 1 alt
                                                real_p_prior[allele1<<2|allele2] = is_validated ? para_het_val_r : para_het_unval_r;
                                        }
                                        else if (allele1 == allele2) {
                                                real_p_prior[allele1<<2|allele2] =  is_validated ? para_althom_val_r : para_althom_unval_r;
                                        }
                                        else {
                                                // altHET: 2 diff alt base
                                                real_p_prior[allele1<<2|allele2] = is_validated ? para_het_val_r : para_het_unval_r;
                                        }
                                }
                        }
                        else {
                                // Real HapMap Sites
                                if(freq[allele1]>0 && freq[allele2]>0) {
                                        real_p_prior[allele1<<2|allele2] = (allele1==allele2 ? 1 : (2*para_het_val_r)) * freq[allele1] * freq[allele2];
                                }
                        }
                }
        }

        return 1;
}

__global__
void posterior_kernel(
                rate_t* d_type_prob,
                char* d_type2,
                int* d_qual1,
                int* d_qual2,
                int* d_qual3,

                rate_t* d_typeLikelihood,
                rate_t* d_p_prior,
                rate_t* d_p_rank,
                uint32* d_sparseBaseInfo,
                uint32* d_numNonZeroPerSite,
                int* d_count_uniArr,
                int* d_q_sumArr,
                int* d_dep_uniArr,
                int* d_count_allArr,
                char* d_oriArr,

                char* d_base1,
                char* d_base2,
                char* d_type1,

                bool* d_validatedArr, bool* d_hapmap_siteArr, bool* d_indel_siteArr, rate_t* d_freqArr,
		
		bool* d_nflag,
                rate_t para_het_val_r, rate_t para_het_unval_r, rate_t para_althom_val_r, rate_t para_althom_unval_r, bool para_is_monoploid, bool para_refine_mode, 

                uint64_t call_length,
		uint32 maxNumNonZero

)
{
        const uint32 numTotalThread = NUM_TOTAL_THREAD;
        const uint32 globalThreadOffset = GLOBAL_THREAD_OFFSET;

        for (int siteID = globalThreadOffset; siteID < call_length; siteID += numTotalThread) {

	        char allele1, allele2, genotype, type2;
        	int qual1, qual2, qual3, all_count1, all_count2, all_count3;
        	//char base3;
	        double real_p_prior[16];
		rate_t type_prob[16 + 1];
		for(int i = 0; i < 17; i++){type_prob[i] = 0.0;}

                if(d_nflag[siteID]) {
                       continue;
                }
                //adjust base1, base2, and base3
                d_base1[siteID] = 0, d_base2[siteID] = 0; //, base3 = 0;
                qual1 = -1, qual2 = -2, qual3 = -3;
                all_count1 = 0, all_count2 = 0, all_count3 = 0;
                if (d_dep_uniArr[siteID]) {
                        for (int i = 0; i != 4; i++) {
                                if (d_q_sumArr[siteID * QSUM_SIZE + i] >= qual1) {
                                        //base3 = d_base2[siteID];
                                        qual3 = qual2;
                                        d_base2[siteID] = d_base1[siteID];
                                        qual2 = qual1;
                                        d_base1[siteID] = i;
                                        qual1 = d_q_sumArr[siteID * QSUM_SIZE + i];
                                } else if (d_q_sumArr[siteID * QSUM_SIZE + i] >= qual2) {
                                        //base3 = d_base2[siteID];
                                        qual3 = qual2;
                                        d_base2[siteID] = i;
                                        qual2 = d_q_sumArr[siteID * QSUM_SIZE + i];
                                } else if (d_q_sumArr[siteID * QSUM_SIZE + i] >= qual3) {
                                        //base3 = i;
                                        qual3 = d_q_sumArr[siteID * QSUM_SIZE + i];
                                } else {
                                        ;
                                }
                        }
                        if (qual1 == 0) {
                                d_base1[siteID] = (d_oriArr[siteID] & 7);
                        } else if (qual2 == 0 && d_base1[siteID] != (d_oriArr[siteID] & 7)) {
                                d_base2[siteID] = (d_oriArr[siteID] & 7);
                        } else {
                                ;
                        }
                } else {
                        for (int i = 0; i != 4; i++) {
                                if (d_count_allArr[siteID * COUNTALL_SIZE + i] >= all_count1) {
                                        //base3 = d_base2[siteID];
                                        all_count3 = all_count2;
                                        d_base2[siteID] = d_base1[siteID];
                                        all_count2 = all_count1;
                                        d_base1[siteID] = i;
                                        all_count1 = d_count_allArr[siteID * COUNTALL_SIZE + i];
                                } else if (d_count_allArr[siteID * COUNTALL_SIZE + i] >= all_count2) {
                                        //base3 = d_base2[siteID];
                                        all_count3 = all_count2;
                                        d_base2[siteID] = i;
                                        all_count2 = d_count_allArr[siteID * COUNTALL_SIZE + i];
                                } else if (d_count_allArr[siteID * COUNTALL_SIZE + i] >= all_count3) {
                                        //base3 = i;
                                        all_count3 = d_count_allArr[siteID * COUNTALL_SIZE + i];
                                } else {
                                        ;
                                }
                        }
                        if (all_count1 == 0) {
                                d_base1[siteID] = (d_oriArr[siteID] & 7);
                        } else if (all_count2 == 0 && d_base1[siteID] != (d_oriArr[siteID] & 7)) {
                                d_base2[siteID] = (d_oriArr[siteID] & 7);
                        } else {
                                ;
                        }
                }

                //calculate posterier probability
//                memcpy(real_p_prior, &d_p_prior[((ubit64_t) d_oriArr[siteID] & 0x7) << 4], sizeof(double) * 16);
		for(int i = 0; i < 16; i++){real_p_prior[i] = d_p_prior[(((ubit64_t) d_oriArr[siteID] & 0x7) << 4) + i];}
                if ((d_oriArr[siteID] & 0x8) && para_refine_mode) {
			gpu_snp_p_prior_gen(real_p_prior, d_oriArr[siteID], 
						d_indel_siteArr[siteID], d_validatedArr[siteID], d_hapmap_siteArr[siteID], &d_freqArr[siteID * 4], 
						para_het_val_r, para_het_unval_r, para_althom_val_r, para_althom_unval_r);
		}

                type2 = 16;
		d_type1[siteID] = 16;
               	for (allele1 = 0; allele1 != 4; allele1++) {
                       	for (allele2 = allele1; allele2 != 4; allele2++) {
                               	genotype = allele1 << 2 | allele2;
                               	if (para_is_monoploid && allele1 != allele2) {
                                       	continue;
                               	}
				
				type_prob[genotype] = d_typeLikelihood[siteID * LIKELY_SIZE + genotype]+ log10(
                                               real_p_prior[genotype]);
                               	if (type_prob[genotype] >= type_prob[d_type1[siteID]] || d_type1[siteID] == 16) {
                                       	type2 = d_type1[siteID];
                                       	d_type1[siteID] = genotype;
                               	} else if (type_prob[genotype] >= type_prob[type2]
                                               	|| type2 == 16) {
                                       	type2 = genotype;
                               	} else {
                                       	;
                               	}
                       	}
               	}
		for(int i = 0; i < 17; i++){
			d_type_prob[siteID * 17 + i] = type_prob[i];
		}
		d_type2[siteID] = type2;
		d_qual1[siteID] = qual1;
		d_qual2[siteID] = qual2;
		d_qual3[siteID] = qual3;
                //rank_sum test and something else
        	
	}
}
__global__ 
void posterior_kernel2(
		rate_t* d_type_prob,
		char* d_type2,
		int* d_qual1,
		int* d_qual2,
		int* d_qual3,

                rate_t* d_p_rank,
                uint32* d_sparseBaseInfo,
                uint32* d_numNonZeroPerSite,
                int* d_count_uniArr,
                int* d_q_sumArr,

                char* d_base1,
                char* d_base2,
                char* d_type1,
                int* d_q_cns,
                double* d_rank_sum_test_value,

                bool* d_nflag,
		bool para_rank_sum_mode,
		char para_q_max_min,
                ubit64_t read_len,
                uint64_t call_length,
                uint32 maxNumNonZero,
                int* d_debug
	
)
{
        const uint32 numTotalThread = NUM_TOTAL_THREAD;
        const uint32 globalThreadOffset = GLOBAL_THREAD_OFFSET;

        for (int siteID = globalThreadOffset; siteID < call_length; siteID += numTotalThread) {

                if(d_nflag[siteID]) {
                       continue;
                }

                if (para_rank_sum_mode) {
                        d_rank_sum_test_value[siteID]
                                = gpu_rank_test(d_sparseBaseInfo, siteID, d_numNonZeroPerSite[siteID], d_count_uniArr, d_p_rank, d_type1[siteID], read_len, para_q_max_min, maxNumNonZero, d_debug);
                } else {
                        d_rank_sum_test_value[siteID] = 1.0;
                }

                if (d_rank_sum_test_value[siteID] == 0.0) {
                        d_q_cns[siteID] = 0;
                } else {
                        d_q_cns[siteID] = (int) (10 * (d_type_prob[siteID * 17 + d_type1[siteID]] - d_type_prob[siteID * 17 + d_type2[siteID]])
                                       + 10 * log10(d_rank_sum_test_value[siteID]));
                }
                if ((d_type1[siteID] & 3) == ((d_type1[siteID] >> 2) & 3)) {
                        if (d_qual1[siteID] > 0 && d_base1[siteID] != (d_type1[siteID] & 3)) {
                                d_q_cns[siteID] = 0;
                        } else if (d_q_cns[siteID] > d_qual1[siteID] - d_qual2[siteID]) {
                                d_q_cns[siteID] = d_qual1[siteID] - d_qual2[siteID];
                        } else {
                                ;
                        }
                } else {
                        if (d_q_sumArr[siteID * QSUM_SIZE + d_base1[siteID]] > 0
                                && d_q_sumArr[siteID * QSUM_SIZE + d_base2[siteID]] > 0
                                && d_type1[siteID]  == (d_base1[siteID] < d_base2[siteID] ? (d_base1[siteID] << 2 | d_base2[siteID]) : (d_base2[siteID] << 2
                                                       | d_base1[siteID]))) {
                                if (d_q_cns[siteID] > d_qual2[siteID] - d_qual3[siteID]) {
                                        d_q_cns[siteID] = d_qual2[siteID] - d_qual3[siteID];
                                }
                        } else {
                                        d_q_cns[siteID] = 0;
                        }
                }
                if (d_q_cns[siteID] > 99) {
                        d_q_cns[siteID] = 99;
                }
                if (d_q_cns[siteID] < 0) {
                        d_q_cns[siteID] = 0;
                }
	}
}
void cuda_posterior(
	st_likelihood* likelihoodData,
	char* oriArr,

	char* base1,
	char* base2,
	char* type1,
	int* q_cns,
	double* rank_sum_test_value,
	
	bool* validatedArr, bool* hapmap_siteArr, bool* indel_siteArr, rate_t* freqArr,
		
	rate_t para_het_val_r, rate_t para_het_unval_r, rate_t para_althom_val_r, rate_t para_althom_unval_r, bool para_rank_sum_mode, bool para_is_monoploid, bool para_refine_mode, char para_q_max_min,

	ubit64_t read_len,
	int win_size,
	uint64_t call_length,
	uint32 maxNumNonZero,
	uint32 numBlock, uint32 numThread
)
{
//unsigned int timer = 0;
//startTimer(&timer);

	TOGPU(d_oriArr, oriArr, sizeof(char) * win_size);
	TOGPU(d_validatedArr, validatedArr, sizeof(bool) * win_size);
	TOGPU(d_hapmap_siteArr, hapmap_siteArr, sizeof(bool) * win_size);
	TOGPU(d_indel_siteArr, indel_siteArr, sizeof(bool) * win_size);
	TOGPU(d_freqArr, freqArr, sizeof(rate_t) * 4 * win_size);
	TOGPU(d_nflag, likelihoodData->nflag, sizeof(bool) * win_size);

cudaThreadSynchronize();
cutilCheckMsg("post copy 1!!!");
//memcpyTimePost += endTimer(&timer, "copy 1");
//memcpySizePost += sizeof(double) * 64*64*2048 + sizeof(rate_t) * 8 * 4 * 4 + sizeof(char) * win_size + sizeof(int) * win_size + sizeof(bool) * win_size + sizeof(bool) * win_size * 3 + sizeof(rate_t) * 4 * win_size;

//timer = 0;
//startTimer(&timer);
	posterior_kernel<<<numBlock, numThread>>>(
                d_type_prob,
                d_type2,
                d_qual1,
                d_qual2,
                d_qual3,

		d_typeLikelihood,
		d_p_prior,
		d_p_rank,
		d_sparseBaseInfo,
		d_numNonZeroPerSite,
		d_count_uniArr,
		d_q_sumArr,
		d_dep_uniArr,
		d_count_allArr,
		d_oriArr,
		d_base1,
		d_base2,
		d_type1,
				
	        d_validatedArr, d_hapmap_siteArr, d_indel_siteArr, d_freqArr,
		d_nflag,

        	para_het_val_r, para_het_unval_r, para_althom_val_r, para_althom_unval_r, para_is_monoploid, para_refine_mode,

        	call_length,
		maxNumNonZero
	);
        posterior_kernel2<<<numBlock, numThread>>>(
		d_type_prob,
		d_type2,
		d_qual1,
		d_qual2,
		d_qual3,
                d_p_rank,
                d_sparseBaseInfo,
                d_numNonZeroPerSite,
                d_count_uniArr,
                d_q_sumArr,
                d_base1,
                d_base2,
                d_type1,
                d_q_cns,
                d_rank_sum_test_value,

                d_nflag,
		para_rank_sum_mode,
		para_q_max_min,
                read_len,
                call_length,
                maxNumNonZero,
                d_debug
        );

cudaThreadSynchronize();
//cutilCheckMsg("cuda_posterior!!!!!!!!!!");
//kernelTimePost += endTimer(&timer, "posterior kernel");


//timer = 0;
//startTimer(&timer);
	
	FROMGPU(base1, d_base1, sizeof(char) * win_size);
	FROMGPU(base2, d_base2, sizeof(char) * win_size);
	FROMGPU(type1, d_type1, sizeof(char) * win_size);
	FROMGPU(q_cns, d_q_cns, sizeof(int) * win_size);
	FROMGPU(rank_sum_test_value, d_rank_sum_test_value, sizeof(double) * win_size);

/*
        int* debug = NULL;
        debug = (int*)malloc(sizeof(int)*20);
        memset(debug, 0, sizeof(int)*20);
        FROMGPU(debug, d_debug, sizeof(int) * 20);
        for(int i = 0; i < 5; i++){cerr<<debug[i]<<"$"<<endl;}
*/
cudaThreadSynchronize();
cutilCheckMsg("post copy 2!!!1");
//memcpyTimePost += endTimer(&timer, "copy 2");
//memcpySizePost += sizeof(char) * win_size * 3 + sizeof(int) * win_size + sizeof(double) * win_size;

//cerr<<"memcpyTimePost:"<<memcpyTimePost/1000.0<<endl;
//cerr<<"kernelTimePost:"<<kernelTimePost/1000.0<<endl;
//cerr<<"memcpySizePost:"<<memcpySizePost/(1024 * 1024 * 1024)<<endl;
	

}

#endif
