#ifndef __COMMON_H__
#define __COMMON_H__

#include <gsnp_util.h>

#include <string>

using namespace std;


/**
 * configureation parameters
 * */
#define         DEVICE_ID 		(0)
#define         MAX_NUM_NON_ZERO 	(256)			//adjust this parameter larger, if any
const int 	global_win_size = 	32*8192;		//set this variable smaller if you have a small GPU memory
const int 	_qual_length 	= 	200;

/**
 * temp file names
 * */
//in.ab        in.hit       in.mismatch  in.position  in.qual      in.read      in.readlen   in.strand
#define AB_FILE			"in.ab"
#define HIT_FILE		"in.hit"
#define MISMATCH_FILE		"in.mismatch"
#define POSITION_FILE		"in.position"
#define QUAL_FILE		"in.qual"
#define READ_FILE 		"in.read"
#define READLEN_FILE 		"in.readlen"
#define STRAND_FILE		"in.strand"



/***
 * do not modify these global variables
 * */
const int  _ini_max_num_non_zero_per_site = MAX_NUM_NON_ZERO; 	//this value should be updated after reallocation
const int _read_num_per_window = 4*global_win_size+100;
static string global_chr_name;
static uint64 g_numTotalSite = 0;
///////////////


typedef unsigned long long ubit64_t;
typedef unsigned int ubit32_t;
typedef double rate_t;
typedef unsigned char small_int;
using namespace std;
const size_t capacity = sizeof(ubit64_t)*8/4;
const char abbv[17]={'A','M','W','R','M','C','Y','S','W','Y','T','K','R','S','K','G','N'};
const ubit64_t glf_base_code[8]={1,2,8,4,15,15,15,15}; // A C T G
const ubit64_t glf_type_code[10]={0,5,15,10,1,3,2,7,6,11};// AA,CC,GG,TT,AC,AG,AT,CG,CT,GT

//typedef struct {int x; int y;} pair;





static uint64 n_cns = 0;
static uint64 n_recycle = 0;
static uint64 n_qrecycle = 0;
static double recycle_timer = 0.0f;
static double qrecycle_timer = 0.0f;
static double cns1_timer = 0.0f;
static double cns2_timer = 0.0f;
static double cns3_timer = 0.0f;
static double cns4_timer = 0.0f;
static double update_timer = 0.0f;

#define MAX_QSCORE (64)
#define BASEINFO_MAX (0xffffffff) //memset as the max. value, for sorting
#define MB_1 (1024.0*1024.0) //1 MB
#define LIKELY_SIZE (17)
#define MAX_READ_LENGTH (256)
#define BASEINFO_SIZE (4*2*64*MAX_READ_LENGTH)
#define COUNTUNI_SIZE (4)
#define QSUM_SIZE (4)
#define COUNTALL_SIZE (4)

#define PMATRIX_SIZE (uint64)(256*256*4*4)
//#define max_readNum (3*global_win_size+100)

/* 2D-array access */
#define SITE_MAJOR

#define IDX_2D(row, column, pitch) ((row)*(pitch) + (column))
#ifdef SITE_MAJOR
	#define BASEINFO_IDX(siteId, offset, numsite_pitch, numoffset_pitch) ((siteId)*(numoffset_pitch) + (offset))        // site-major
#else
	#define BASEINFO_IDX(siteId, offset, numsite_pitch, numoffset_pitch) ((offset)*(numsite_pitch) + (siteId))          // offset-major
#endif

inline void set_countuni(int* countuni, const uint64 siteId, const int i, const int value) {
        countuni[IDX_2D(siteId, i, COUNTUNI_SIZE)] = value;
}

inline int get_countuni(int* countuni, const uint64 siteId, const int i) {
        return countuni[IDX_2D(siteId, i, COUNTUNI_SIZE)];
}

inline void set_baseInfo(small_int* baseInfo, const uint64 siteId, const int i, const small_int value) {
        baseInfo[IDX_2D(siteId, i, BASEINFO_SIZE)] = value;
}

inline small_int get_baseInfo(small_int* baseInfo, const uint64 siteId, const int i) {
        return baseInfo[IDX_2D(siteId, i, BASEINFO_SIZE)];
}

inline bool isAllN(const bool* nflag, const uint64 numElement) {
	bool allN = true;
	
	for(uint64 i = 0; i < numElement; i++) {
		if(!nflag[i]) {
			allN = false;
			break;
		}
	}

	return allN;
}

class Parameter {
public:
	char q_min; // The char stands for 0 in fastq
	char q_max; // max quality score
	small_int read_length; // max read length
	bool is_monoploid; // Is it an monoploid? chrX,Y,M in man.
	bool is_snp_only;  // Only output possible SNP sites?
	bool refine_mode; // Refine prior probability using dbSNP
	bool rank_sum_mode; // Use rank sum test to refine HET quality
	bool binom_mode; // Use binomial test to refine HET quality
	bool transition_dominant; // Consider transition/transversion ratio?
	int glf_format; // Generate Output in GLF format: New File Definitions! Since May 1, 2009
	bool region_only; // Only report consensus in specified region
	string glf_header; // Header of GLF format
	rate_t althom_novel_r, het_novel_r; // Expected novel prior
	rate_t althom_val_r, het_val_r; // Expected Validated dbSNP prior
	rate_t althom_unval_r, het_unval_r; // Expected Unvalidated dbSNP prior
	rate_t global_dependency, pcr_dependency; // Error dependencies, 1 is NO dependency
// Default onstruction
	Parameter(){
		q_min = 64;
		q_max = 64+40;
		read_length = 45;
		is_monoploid = is_snp_only = refine_mode = rank_sum_mode = binom_mode = transition_dominant = region_only =false;
		glf_format = 0;
		glf_header = "";
		althom_novel_r=0.0005, het_novel_r=0.0010;
		althom_val_r=0.05, het_val_r=0.10;
		althom_unval_r=0.01, het_unval_r=0.02;
		global_dependency= log10(0.9), pcr_dependency= log10(0.5); // In Log10 Scale
	};
};


#endif /* __COMMON_H__ */


