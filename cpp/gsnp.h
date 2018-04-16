/*
 * gsnp_kernel.h
 * this is a warp header for the GPU kenrel functions
 *
 * Created on:	Dec 23, 2010
 * Author:	lumian
 */

#ifndef GSNP_KERNEL_H_
#define GSNP_KERNEL_H_

/* the original header file */
//#include "soap_snp.h"

/* the improved code */
#define P_MAT_IDX(q_adj, coord, allele, o_base) (((ubit64_t)(q_adj))<< 12)|((coord)<< 4)|((allele)<< 2)|(o_base)


/* used for the GPU kernel files */
#include <gsnp_util.h>

#include "common.h"

#include "soap_snp.h"

void GPUInit();


struct st_likelihood {

	/* global parameters */
	char qmax;
	char qmin;
	int readLength;
        rate_t globalDependency;
        rate_t pcrDependency;
	uint64 maxCallLength;
	uint64 winSize;
	uint64 numSite;	//max., for a window (window_size + read_length)
	rate_t* p_matrix;

	/* window-based parameters */
	uint64 callLength;
	bool* nflag;
	//two of four information
	int* countuni;
//	small_int* baseInfo;

	/* output */
	rate_t* typeLikelihood;

	st_likelihood(const Parameter* param, const uint64 in_winSize, const rate_t* mat, 
			int* countuni_buf) {

		this->qmax = param->q_max;
		this->qmin = param->q_min;
		this->readLength = param->read_length;
		this->globalDependency = param->global_dependency;
		this->pcrDependency = param->pcr_dependency;
		this->winSize = in_winSize;
		this->maxCallLength = in_winSize;
		this->numSite = winSize + readLength;
		this->p_matrix = (rate_t*)malloc(sizeof(rate_t)*PMATRIX_SIZE);
		memcpy(this->p_matrix, mat, sizeof(rate_t)*PMATRIX_SIZE);

		this->callLength = 0;
		this->nflag = (bool*)malloc(sizeof(bool)*maxCallLength);
		//this->countuni = (int*)malloc(sizeof(int)*COUNTUNI_SIZE*numSite);
		//this->baseInfo = (small_int*)malloc(sizeof(small_int)*BASEINFO_SIZE*numSite);
		this->countuni = countuni_buf;
//		this->baseInfo = baseInfo_buf;


		this->typeLikelihood = (rate_t*)malloc(sizeof(rate_t)*LIKELY_SIZE*maxCallLength);
	};


	//rearrange the data layout
	//TODO: the rearrange can be avoided, just use the same memory layout on the GPU and CPU
	//in fact, we only need a memory copy rather than data rearrange
	void prepare(const uint64 in_callLength) {
		assert(in_callLength <= maxCallLength);
		this->callLength = in_callLength;
	}

	void release() {
                this->qmax = 0;
                this->qmin = 0;
                this->readLength = 0;
                this->globalDependency = 0;
                this->pcrDependency = 0;
                this->winSize = 0;
                this->maxCallLength = 0;
                this->numSite = 0;

		if(nflag) {
			free(nflag);
                	this->nflag = NULL;
		}

		if(p_matrix) {
			free(p_matrix);
			this->p_matrix = NULL;
		}
	
		/*
		if(countuni) {
			free(countuni);
			countuni = NULL;
		}

		if(baseInfo) {
			free(baseInfo);
			baseInfo = NULL;
		}
		*/

		if(typeLikelihood) {
                	free(typeLikelihood);
			typeLikelihood = NULL;
		}

	};
};


void GSNPKernelInit(const uint64 winSize, const int readLength,
                    const rate_t pcr_dependency, const rate_t global_dependency,
                    const rate_t* h_pmatrix, const rate_t* p_prior, const rate_t* p_rank);


void cuda_likelihood(st_likelihood* likelihood, const int numBlock = 1200, const int numThread = 64);

struct st_update{
    char *inputUpdateChar;
    int *inputUpdateInt;
    
    char *readArr;
    char *qualArr;
    string *char_nameArr;
    int *hitArr;
    int *read_lenArr;
    int *positionArr;
    int *siteMatchArr;
    int *mismatchArr;
    char *abArr;
    char *strandArr;
    int *readLen_hitArr;
    char qmin;
    uint32 _read_length;

    st_update(const Parameter* param){
	this->_read_length = param->read_length;

        int inputSizeChar = sizeof(char)* _read_length * _read_num_per_window
                + sizeof(char) * _read_length * _read_num_per_window
                + sizeof(char) * _read_num_per_window;

        int inputSizeInt = sizeof(int) * _read_num_per_window
                + sizeof(int) * global_win_size * 2;


        this->qmin = param->q_min;
        this->inputUpdateChar = (char*)malloc(inputSizeChar);
        this->inputUpdateInt = (int*)malloc(inputSizeInt);

        this->readArr = this->inputUpdateChar;
        this->qualArr = this->inputUpdateChar + _read_length * _read_num_per_window;
        this->strandArr = this->inputUpdateChar +  _read_length * _read_num_per_window
                                        +  _read_length * _read_num_per_window;

	this->readLen_hitArr = this->inputUpdateInt;
        this->siteMatchArr = this->inputUpdateInt +  _read_num_per_window;
        this->mismatchArr = (int *)malloc(sizeof(int) * _read_num_per_window);
        memset(mismatchArr,0,sizeof(int) * _read_num_per_window);

        this->char_nameArr = new string[sizeof(string) * _read_num_per_window];
        this->abArr = (char *)malloc(sizeof(char) * _read_num_per_window);
        memset(abArr,0,sizeof(char) * _read_num_per_window);
        this->positionArr = (int *)malloc(sizeof(int) * _read_num_per_window);
        memset(positionArr,0,sizeof(int) * _read_num_per_window);

        this->hitArr = (int *)malloc(sizeof(int) * _read_num_per_window);
        memset(hitArr,0,sizeof(int) * _read_num_per_window);

        this->read_lenArr = (int *)malloc(sizeof(int) * _read_num_per_window);
        memset(read_lenArr,0,sizeof(int) * _read_num_per_window);


/*
	this->readArr = (char*)malloc(sizeof(char) * _read_length * _read_num_per_window);
	this->qualArr = (char*)malloc(sizeof(char) * _read_length * _read_num_per_window);

        this->char_nameArr = new string[sizeof(string) * _read_num_per_window];

        this->hitArr = (int *)malloc(sizeof(int) * _read_num_per_window);
        memset(hitArr,0,sizeof(int) * _read_num_per_window);

        this->read_lenArr = (int *)malloc(sizeof(int) * _read_num_per_window);
        memset(read_lenArr,0,sizeof(int) * _read_num_per_window);

        this->positionArr = (int *)malloc(sizeof(int) * _read_num_per_window);
        memset(positionArr,0,sizeof(int) * _read_num_per_window);

	this->siteMatchArr = (int*)malloc(sizeof(int) * global_win_size);
	memset(siteMatchArr, 0, sizeof(int) * global_win_size);

        this->mismatchArr = (int *)malloc(sizeof(int) * _read_num_per_window);
        memset(mismatchArr,0,sizeof(int) * _read_num_per_window);

        this->abArr = (char *)malloc(sizeof(char) * _read_num_per_window);
        memset(abArr,0,sizeof(char) * _read_num_per_window);

        this->strandArr = (char *)malloc(sizeof(char) * _read_num_per_window);
        memset(strandArr,0,sizeof(char) * _read_num_per_window);
*/
    };
    void updateSoap(Soap_format soap, int readCount){

	assert(readCount <= _read_num_per_window);
	assert(soap.qual.size() <= _read_length);

        memcpy(this->readArr + readCount * _read_length, soap.read.c_str(), _read_length);
	memcpy(this->qualArr + readCount * _read_length, soap.qual.c_str(), soap.qual.size());
        this->char_nameArr[readCount] = soap.chr_name;
	this->readLen_hitArr[readCount] = (int)(
						((soap.read_len) & 0xffff)<<16
						|(soap.hit)
						);
        this->positionArr[readCount] = soap.position;
        this->mismatchArr[readCount] = soap.mismatch;
        this->abArr[readCount] = soap.ab;
        this->strandArr[readCount] = soap.strand;
    };

    void updateSoap2(int readCount, int tmpVal){

        assert(readCount <= _read_num_per_window);

        memcpy(this->readArr + readCount * _read_length, this->readArr + tmpVal * _read_length, _read_length);
        memcpy(this->qualArr + readCount * _read_length, this->qualArr + tmpVal * _read_length, _read_length);
        this->char_nameArr[readCount] = this->char_nameArr[tmpVal];
        this->readLen_hitArr[readCount] = (int)(
                                                ((this->read_lenArr[tmpVal]) & 0xffff)<<16
                                                |(this->hitArr[tmpVal])
                                                );
        this->positionArr[readCount] = this->positionArr[tmpVal];
        this->mismatchArr[readCount] = this->mismatchArr[tmpVal];
        this->abArr[readCount] = this->abArr[tmpVal];
        this->strandArr[readCount] = this->strandArr[tmpVal];
    };


    void release(){
                this->qmin = 0;
		if(inputUpdateChar){
			free(inputUpdateChar);
			this->inputUpdateChar = NULL;
		}
		if(inputUpdateInt){
			free(inputUpdateInt);
			this->inputUpdateInt = NULL;
		}
		if(readLen_hitArr){
			free(readLen_hitArr);
			this->readLen_hitArr = NULL;
		}
		if(siteMatchArr){
			free(siteMatchArr);
			this->siteMatchArr = NULL;
		}
                if(readArr) {
			free(readArr);
                        this->readArr = NULL;
                }
                if(qualArr) {
			free(qualArr);
                        this->qualArr = NULL;
                }
                if(char_nameArr) {
                        delete[] char_nameArr;
                        this->char_nameArr = NULL;
                }
                if(hitArr) {
                        free(hitArr);
                        this->hitArr = NULL;
                }
                if(read_lenArr) {
                        free(read_lenArr);
                        this->read_lenArr = NULL;
                }
                if(positionArr) {
                        free(positionArr);
                        this->positionArr = NULL;
                }
                if(mismatchArr) {
                        free(mismatchArr);
                        this->mismatchArr = NULL;
                }
                if(abArr) {
                        free(abArr);
                        this->abArr = NULL;
                }
                if(strandArr) {
                        free(strandArr);
                        this->strandArr = NULL;
                }
    };

};
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
);


void gold_update(
	struct st_update* updateData,
	small_int* base_infoArr,

        int* count_uniArr,
        int* q_sumArr,
        int* count_allArr,
        int* depthArr,
        int* repeat_timeArr,
        int* dep_uniArr,
        int readCount
);
double gold_table_test(double *p_rank, int n1, int n2, double T1, double T2);
double gold_normal_value(double z);
double gold_normal_test(int n1, int n2, double T1, double T2);
double gold_rank_test(
	const uint32* sparse_base_infoArr, 
	const uint64 siteID, 
	uint32 numNonZero,
	uint64 read_len, 
	int* count_uni_buf, 
	int q_max_min,
	char best_type, 
	double* p_rank
);  
void gold_likelihoodSparse(rate_t* typeLikelihood,
                           uint32* sparseBaseInfo, 
			   const uint32 numSite, const uint32 pitch, 
			   const uint32* numNonZeroPerSite,
                           const uint32 readLength, const int q_maxmin,
                           const rate_t pcrDependency, const rate_t globalDependency,
                           const rate_t* p_matrix
);
uint32 gold_generateSparseBaseinfo(uint32** sparseBaseInfo, uint32* numNonZeroPerSite, 
			      small_int* baseInfo, 
			      const uint32 readLength,
			      const uint32 numSite, const int q_maxmin 
);
bool sparse_base_info_check(
			uint32* sparse_new,
			uint32* numPerSite_new,
			uint32* sparse_base,
			uint32* numPerSite_base,
			uint32  numSite
);



/* temp APIs, for test and debug */
//void copy_likelihoodSparse(const uint32* h_sparseBaseInfo, const uint32 pitch, const uint32 maxNumNonZeroPerSite,
//                           const uint32* h_numNonZeroPerSite, const uint32 callLength);

void copy_likelihoodSparse(const uint32* h_sparseBaseInfo, const uint64 size,
                           const uint32* h_numNonZeroPerSite, const uint32 callLength);

/* the real GPU APIs, cleanup the code above later */

void gold_pos_match(int* positionArr, int* siteMatchArr, int readCount);

void gpu_likelihoodSparse(rate_t* h_typeLikelihood,
                          const uint32 callLength, const uint32 pitch,
                          const uint32 maxNumNonZero, 
                          const uint32 numBlock = 1200, const uint32 numThread = 128);

void gpu_segSort(const uint32 numSite, const uint32 offsetPitch );

void gold_segSort(uint32* baseInfo, const uint32 numSite, const uint32* siteSize,
                 const uint32 numSitePitch, const uint32 offsetPitch);


void cuda_likelihoodSparse(rate_t* h_typeLikelihood,
                           const uint32 numSite, const uint32 maxNumNonZero,
                           const uint32 sitePitch, const uint32 offsetPitch,
                           const uint32 numBlock = 1200, const uint32 numThread = 128);

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
        const uint32 numBlock = 1200,
        const uint32 numThread = 256);

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
);

struct st_snp{
        bool* validatedArr;
        bool* hapmap_siteArr;
        bool* indel_siteArr;
        rate_t* freqArr;
        st_snp(int win_size){
                validatedArr = (bool*) malloc(sizeof(bool) * win_size);
                memset(validatedArr, false, sizeof(bool) * win_size);
                hapmap_siteArr = (bool*) malloc(sizeof(bool) * win_size);
                memset(hapmap_siteArr, false, sizeof(bool) * win_size);
                indel_siteArr = (bool*) malloc(sizeof(bool) * win_size);
                memset(indel_siteArr, false, sizeof(bool) * win_size);
                freqArr = (rate_t*) malloc(sizeof(rate_t) * 4 * win_size);
                memset(freqArr, 0, sizeof(rate_t) * 4 * win_size);
        };
        void release(){
                free(validatedArr);
                free(hapmap_siteArr);
                free(indel_siteArr);
                free(freqArr);
        };
};

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
        const uint32 numBlock = 1200, const uint32 numThread = 256
);


void cuda_init(
        ubit64_t win_size,
        ubit64_t read_len,
        ubit64_t start,
        int* posArr,
        const uint32 numBlock = 1200, const uint32 numThread = 256
);
void cuda_recycle(
        ubit64_t win_size,
        ubit64_t read_len,
        int* posArr,
        const uint32 numBlock = 1200, const uint32 numThread = 256
);
void cuda_quickRecycle(
        ubit64_t win_size,
        ubit64_t read_len,
        int* posArr,
        int* count_uniArr,
        int* q_sumArr,
        int* count_allArr,
        int* depthArr,
        int* repeat_timeArr,

        const uint32 numBlock = 1200, const uint32 numThread = 256
);


void printTable();

#endif /* GSNP_KERNEL_H_ */
