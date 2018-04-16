#ifndef __KERNEL_COMMON_H__
#define __KERNEL_COMMON_H__


/* global variables */
static uint64 g_winSize = 0;
static uint64 g_readLength = 0;
static uint64 g_pitch = 0;
static uint64 g_maxNumNonZeroPerSite = 0;
static uint64 g_count = 0;      //a helper counter
static rate_t g_pcrDependency = 0;
static rate_t g_globalDependency = 0;


/* global GPU memory buffer */
rate_t* d_typeLikelihood = NULL;
uint32* d_sparseBaseInfo = NULL;
uint32* d_sparseBaseInfo2 = NULL;
uint32* d_numNonZeroPerSite = NULL;
int* d_pcrDepCount = NULL;
rate_t* d_newpmatrix = NULL;
__constant__ double d_log10Table[MAX_QSCORE];

/*buffer used in kernel/posterior_kernel.cu*/
rate_t* d_p_prior = NULL;
rate_t* d_p_rank = NULL;
char* d_oriArr = NULL;
char* d_base1 = NULL;
char* d_base2 = NULL;
char* d_type1 = NULL;
int* d_q_cns = NULL;
double* d_rank_sum_test_value = NULL;
bool* d_validatedArr = NULL;
bool* d_hapmap_siteArr = NULL;
bool* d_indel_siteArr = NULL;
rate_t* d_freqArr = NULL;
bool* d_nflag = NULL;

//int* d_posArr = NULL;

/*buffer used in kernel/update_kernel.cu*/
uint32* d_maxNumNonZero = NULL;
int* d_count_uniArr = NULL;
int* d_q_sumArr = NULL;
int* d_count_allArr = NULL;
int* d_depthArr = NULL;
int* d_repeat_timeArr = NULL;
int* d_dep_uniArr = NULL;
char* d_readArr = NULL;
char* d_qualArr = NULL;
int* d_hitArr = NULL;
int* d_read_lenArr = NULL;
int* d_readLen_hitArr = NULL;
int* d_siteMatchArr = NULL;
//int* d_mismatchArr = NULL;
char* d_strandArr = NULL;

char* d_inputUpdateChar = NULL;
int* d_inputUpdateInt = NULL;

int* d_positionArr = NULL;
#endif 


