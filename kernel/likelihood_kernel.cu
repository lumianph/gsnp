#ifndef __LIKELIHOOD_KERNEL_CU__
#define __LIKELIHOOD_KERNEL_CU__


#include "common_kernel.h"
#include "batchsort.cu"

static unsigned long g_siteId = 0;
#define NEWPMATRIX_SIZE ((PMATRIX_SIZE/4)*10)



/* for better performance */
void materializePMatrix(const rate_t* pmatrix, double* new_pmatrix);

/* using the dense representation for the base_info */
void gold_likelihoodDense(st_likelihood* likelihoodData, rate_t* typeLikelihood);

bool has_nflag(const bool* nflag, const uint64 n);

bool has_allnflag(const bool* nflag, const uint64 n);


__inline__ __host__ __device__
uint32 encodeBaseWord(const uint32 coord, const uint32 strand,
                       const uint32 o_base, const uint32 q_score, const uint32 occ) {
       uint32 word = 0;
       word = (o_base<<23)|((MAX_QSCORE - q_score)<<17)|(coord<<9)|(strand<<8)|(occ);
       return word;
}

__inline__ __host__ __device__
void decodeBaseWord(const uint32 word,
               uint32* coord, uint32* strand, uint32* o_base,
               uint32* q_score, uint32* occ) {

       *o_base = (word&0x01800000)>>23;
       *q_score = MAX_QSCORE - ((word&0x007e0000)>>17);
       *coord = (word&0x0001fe00)>>9;
       *strand = (word&0x00000100)>>8;
       *occ = (word&0x000000ff);
}




//////////////////////////////////////////////////////////




void likelihoodKernelInit(const uint64 winSize, const int readLength, 
			  const rate_t pcr_dependency, const rate_t global_dependency,
			  const rate_t* h_pmatrix) {


	/*variables*/
	g_winSize = winSize;
	g_readLength = readLength;
	g_pitch = g_winSize + g_readLength;
	g_maxNumNonZeroPerSite = MAX_NUM_NON_ZERO;
	g_pcrDependency = pcr_dependency;
	g_globalDependency = global_dependency;
	g_count = 0;

	/*generate new p_matrix on the CPU and GPU*/
	rate_t* h_newpmatrix = (rate_t*)malloc(sizeof(rate_t)*NEWPMATRIX_SIZE);
	materializePMatrix(h_pmatrix, h_newpmatrix);
	GPUMALLOC((void**)&d_newpmatrix, sizeof(rate_t)*NEWPMATRIX_SIZE);
	TOGPU(d_newpmatrix, h_newpmatrix, sizeof(rate_t)*NEWPMATRIX_SIZE);
	free(h_newpmatrix);	

	/* GPU memory prepare */
	GPUMALLOC((void**)&d_typeLikelihood, sizeof(rate_t)*g_winSize*LIKELY_SIZE);	
	GPUMALLOC((void**)&d_sparseBaseInfo, sizeof(uint32)*g_pitch*g_maxNumNonZeroPerSite);
	GPUMALLOC((void**)&d_sparseBaseInfo2, sizeof(uint32)*g_pitch*g_maxNumNonZeroPerSite);
	GPUMALLOC((void**)&d_numNonZeroPerSite, sizeof(uint32)*g_winSize);
	GPUMALLOC((void**)&d_pcrDepCount, sizeof(int)*g_winSize*2*g_readLength);
	batchSortBegin(g_winSize);


	/* the lookup tables */
        double h_log10Table[MAX_QSCORE] = {0};
        for(uint32 i = 0; i < MAX_QSCORE; i++) {
                h_log10Table[i] = log10((double)i);
        }
        cudaMemcpyToSymbol(d_log10Table, h_log10Table, sizeof(double)*MAX_QSCORE);
}


void likelihoodKernelFinalize() {
	g_winSize = 0;
	g_readLength = 0;
	g_pcrDependency = 0;
	g_globalDependency = 0;
	g_pitch = 0;
	g_maxNumNonZeroPerSite = 0;
	g_count = 0;

	GPUFREE(d_newpmatrix);
	GPUFREE(d_typeLikelihood);
	GPUFREE(d_sparseBaseInfo);
	GPUFREE(d_numNonZeroPerSite);
	GPUFREE(d_pcrDepCount);

	batchSortFinalize();
}


__device__ __host__
__inline__
int device_get_countuni(const int* d_countuni, const uint64 siteId, const int i) {
	return d_countuni[siteId*COUNTUNI_SIZE + i];
}

__device__ __host__
__inline__
small_int device_get_baseInfo(const small_int* d_baseInfo, const uint64 siteId, const int i) {
	return d_baseInfo[siteId*BASEINFO_SIZE + i];
}



bool has_nflag(const bool* nflag, const uint64 n) {
	bool is = false;
	for(uint64 i = 0; i < n; i++) {
		if(nflag[i]) {
			is = true;
			break;
		}
	}
	return is;
}



void materializePMatrix(const rate_t* pmatrix, double* new_pmatrix) {


        for(int q_adj = 0; q_adj < 256; q_adj++) {
                for(int coord = 0; coord < 256; coord++) {
                        for(int o_base = 0; o_base < 4; o_base++) {
                                const uint32 new_idxBase = (((uint32)(q_adj))<< 10)|((coord)<< 2)|(o_base);

                                int tmp = 0;
                                for(int allele1 = 0; allele1 < 4; allele1++) {
                                        for(int allele2 = allele1; allele2 < 4; allele2++) {
                                                new_pmatrix[new_idxBase*10 + tmp] =
                                                                        log10(
                                                                        0.5*pmatrix[P_MAT_IDX(q_adj, coord, allele1, o_base)]
                                                                                        +
                                                                        0.5*pmatrix[P_MAT_IDX(q_adj, coord, allele2, o_base)]);

                                                tmp++;
                                        }
                                }
                                assert(tmp == 10);
                        }
                }
        }
}






__global__
void likelihoodSparse_kernel(rate_t* d_typeLikelihood, 
			     const uint32* d_newBaseInfo, const uint32* d_depth, const uint32 maxDepth,
			     int* d_pcrDepCount, const double* d_newpmatrix,
                             const uint32 callLength, const uint32 pitch,
			     const uint32 readLength, 
                             const rate_t pcrDependency, const rate_t globalDependency ) { 
	const uint32 numTotalThread = NUM_TOTAL_THREAD;
	const uint32 globalThreadOffset = GLOBAL_THREAD_OFFSET;
	int global_dep_count = 0;
	int* pcr_dep_count = NULL;
	int last_o_base, q_adjusted;
	uint32 o_base, strand, q_score, coord, occ;
	double qscore_log10;
	int tmp = 0;
	uint64 idx = 0;
	rate_t s_typeLikelihood[16];


	for(uint64 siteId = globalThreadOffset; siteId < callLength; siteId += numTotalThread) {
	
#pragma unroll 16
                for(int i = 0; i < 16; i++) {
                        s_typeLikelihood[i] = 0.0;
                }


		global_dep_count = -1;
		pcr_dep_count = d_pcrDepCount + siteId;
		const uint32 depth = d_depth[siteId];


		for(uint64 offset = 0; offset < maxDepth; offset++) {
			//const uint32 word = d_newBaseInfo[offset*callLength + siteId];
		
			//do the real computation
			if(offset < depth) {
				const uint32 word = d_newBaseInfo[BASEINFO_IDX(siteId, offset, pitch, maxDepth)];

				decodeBaseWord(word, &coord, &strand, &o_base, &q_score, &occ);
				qscore_log10 = d_log10Table[q_score];

				if(last_o_base != o_base) {
					global_dep_count = -1;
					for(uint32 i = 0; i < 2*readLength; i++) {
			                        pcr_dep_count[i*callLength] = 0;
                			}
					last_o_base = o_base;
				}

				idx = callLength*(strand*readLength + coord);
				tmp = pcr_dep_count[idx];
				for(small_int k = 0; k < occ; k++) {
					if(tmp == 0) {
						global_dep_count += 1;
					}

                                	tmp++;
                                	q_adjusted = (int)(pow(10.0, (qscore_log10
                                             + (tmp - 1)*pcrDependency
                                             + global_dep_count*globalDependency)) + 0.5);
	                                if(q_adjusted < 1) {
         	                               q_adjusted = 1;
         	                      	}

			
					int tmp = 0;
                                	for (int allele1 = 0; allele1 != 4; allele1++) {
                                        	for (int allele2 = allele1; allele2 != 4; allele2++) {
                                                	s_typeLikelihood[(allele1<<2)|allele2]
                                                        +=
							d_newpmatrix[((((uint32)(q_adjusted))<< 10)|((coord)<< 2)|(o_base))*10+ tmp]; 
							tmp++;
                                        	}
                               		}
				}
				pcr_dep_count[idx] = tmp;
			}
		
			__syncthreads();
		}
		

#pragma unroll 16
                for(int i = 0; i < 16; i++) {
                        d_typeLikelihood[siteId*LIKELY_SIZE + i] = s_typeLikelihood[i];
                }
	}
}


uint64 numTotalSite = 0;


void cuda_likelihoodSparse(rate_t* h_typeLikelihood,
			   const uint32 numSite, const uint32 maxNumNonZero, 
			   const uint32 sitePitch, const uint32 offsetPitch, 
			   const uint32 numBlock, const uint32 numThread) {

	g_siteId += numSite;

        if(maxNumNonZero == 0) {
                memset(h_typeLikelihood, 0, sizeof(rate_t)*LIKELY_SIZE*numSite);
                g_count++;

                return;
        }

	/* 0. parameter checking */
	assert(numSite <= g_winSize);
	assert(maxNumNonZero <= offsetPitch);


	/*single pass*/
	//const uint32 maxArrayLength = getReduceSize(maxNumNonZero);
	//batchSortSP(d_sparseBaseInfo, offsetPitch, numSite, d_numNonZeroPerSite, maxArrayLength); 
	/*multipass*/
	batchSortMP2(d_sparseBaseInfo, offsetPitch, numSite, d_numNonZeroPerSite);


	/* 3. likelihood computation kernel */
	//timer = 0;
	//startTimer(&timer);
	        cutilSafeCall(cudaMemset(d_pcrDepCount, 0, sizeof(int)*2*g_readLength*numSite));
	        cutilSafeCall(cudaMemset(d_typeLikelihood, 0, sizeof(rate_t)*LIKELY_SIZE*numSite));
	        likelihoodSparse_kernel<<<numBlock, numThread>>>(d_typeLikelihood,
                                                d_sparseBaseInfo, d_numNonZeroPerSite, offsetPitch,
                                                d_pcrDepCount, d_newpmatrix,
                                                numSite, sitePitch,
                                                g_readLength,
                                                g_pcrDependency, g_globalDependency);
	        cutilCheckMsg("likelihoodSparse_kernel");
	        cutilSafeCall(cudaThreadSynchronize());
	//endTimer(&timer, "likelihood computation");

	g_count++;
}


//////////////////////////////
//temp functions 

void copy_likelihoodSparse(const uint32* h_sparseBaseInfo, const uint64 size,
                           const uint32* h_numNonZeroPerSite, const uint32 callLength) {

	assert(size <= sizeof(uint32)*g_pitch*g_maxNumNonZeroPerSite);

        TOGPU(d_sparseBaseInfo, h_sparseBaseInfo, size);
        TOGPU(d_numNonZeroPerSite, h_numNonZeroPerSite, sizeof(uint32)*callLength);
}


#endif /* __LIKELIHOOD_KERNEL_CU__ */
 

