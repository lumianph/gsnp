#include <kernel/cuda_header.h>
#include <gsnp_util.h>
#include <assert.h>
#include <vector>
//#include "soap_snp.h"
#include "gsnp.h"
#include "compress.h"
#include "io.h"


using namespace std;


/*for compression*/
Buffer* buffer = NULL;
Column_1_2* c_1_2 = NULL;
Column_3_4_6* c_3_4_6 = NULL;
Column_RLEDict* c_5_7_8_9_14_16 = NULL;
Column_10* c_10 = NULL;
Column_11_12_13* c_11_12_13 = NULL;
Column_15_17* c_15_17 = NULL;

static struct st_likelihood* likelihoodData = NULL;	
static struct st_snp* snpData = NULL;

double real_p_prior[16];
double likelihoods[10];
int* pcr_dep_count = NULL;
char* base1 = NULL;
char* base2 = NULL;
char* type1 = NULL;
int* q_cns = NULL;
double* rank_sum_test_value = NULL; 


bool useCompression = false;
void setCompression() {
	useCompression = true;
}


int Call_win::initialize(ubit64_t start) {
	std::string::size_type i;
	for (i = 0; i != read_len + win_size; i++) {
		sites[i].pos = i + start;
	}

        //important!!!
        memset(sparse_base_infoArr, 0xffff, sizeof(uint32)*(read_len + win_size)*MAX_NUM_NON_ZERO);

	return 1;
}

int Call_win::deep_init(ubit64_t start) {
	int i;
	for (i = 0; i != read_len + win_size; i++) {
		sites[i].pos = i + start;
		sites[i].ori = 0xFF;
		sites[i].depth = 0;
		sites[i].repeat_time = 0;
		sites[i].dep_uni = 0;
		memset(sites[i].count_uni, 0, sizeof(int) * 4);
		memset(sites[i].q_sum, 0, sizeof(int) * 4);
//		memset(sites[i].base_info, 0, sizeof(small_int) * 4 * 2 * 64 * 256);
		memset(sites[i].count_all, 0, sizeof(int) * 4);
	}
        for(i=0; i < read_len + win_size; i++){
             depthArr[i] = sites[i].depth;
             repeat_timeArr[i] = sites[i].repeat_time;
             dep_uniArr[i] = sites[i].dep_uni;
        }



	memset(numNonZeroPerSite, 0, sizeof(uint32)*(win_size+read_len));

	//important!!!
	memset(sparse_base_infoArr, 0xffff, sizeof(uint32)*(read_len + win_size)*MAX_NUM_NON_ZERO);

	return 1;
}

int Call_win::recycle() {
	n_recycle++;
	INIT_TIMER;
	START_TIMER;

	std::string::size_type i;
	for (i = 0; i != read_len; i++) {
		sites[i].pos = sites[i + win_size].pos;
		sites[i].ori = sites[i + win_size].ori;
		sites[i].depth = sites[i + win_size].depth;
		sites[i].repeat_time = sites[i + win_size].repeat_time;
		sites[i].dep_uni = sites[i + win_size].dep_uni;
//		memcpy(sites[i].base_info, sites[i + win_size].base_info,
//				sizeof(small_int) * 4 * 2 * 64 * 256); // 4 types of bases, 2 strands, max quality score is 64, and max read length 256
		memcpy(sites[i].count_uni, sites[i + win_size].count_uni, sizeof(int)*4);
		memcpy(sites[i].q_sum, sites[i + win_size].q_sum, sizeof(int) * 4);
		memcpy(sites[i].count_all, sites[i + win_size].count_all, sizeof(int)*4);
	}

	memcpy(numNonZeroPerSite,numNonZeroPerSite + win_size,sizeof(uint32)*read_len);
	memset(numNonZeroPerSite + read_len, 0, sizeof(uint32)*win_size);
// uint64 numSite = win_size + read_len;

//	memset(sparse_base_info, 0, sizeof(uint32)*numSite*_ini_max_num_non_zero_per_site);
        uint64 numSite = win_size + read_len;
        for(uint32 siteID = 0; siteID < read_len; siteID++){
             //for(uint32 offset = 0; offset < numNonZeroPerSite[siteID]; offset++){
             for(uint32 offset = 0; offset < MAX_NUM_NON_ZERO; offset++){
                     //sparse_base_infoArr[offset*numSite + siteID] = sparse_base_infoArr[offset*numSite + siteID + global_win_size];
                     sparse_base_infoArr[BASEINFO_IDX(siteID, offset, numSite, MAX_NUM_NON_ZERO)] = 
				sparse_base_infoArr[BASEINFO_IDX(siteID + global_win_size, offset, numSite, MAX_NUM_NON_ZERO)];

             }
        }

	for (i = read_len; i != read_len + win_size; i++) {
		sites[i].ori = 0xFF;
		sites[i].pos = sites[i - 1].pos + 1;
		sites[i].depth = 0;
		sites[i].repeat_time = 0;
		sites[i].dep_uni = 0;
		memset(sites[i].count_uni, 0, sizeof(int) * 4);
		memset(sites[i].q_sum, 0, sizeof(int) * 4);
//		memset(sites[i].base_info, 0, sizeof(small_int) * 4 * 2 * 64 * 256);
		memset(sites[i].count_all, 0, sizeof(int) * 4);


		for(uint32 offset = 0; offset < MAX_NUM_NON_ZERO; offset++) {
			sparse_base_infoArr[BASEINFO_IDX(i, offset, numSite, MAX_NUM_NON_ZERO)] = BASEINFO_MAX;
		}
	}

	

        for(i=0; i < read_len + win_size; i++){
	        depthArr[i] = sites[i].depth;
        	repeat_timeArr[i] = sites[i].repeat_time;
             	dep_uniArr[i] = sites[i].dep_uni;
        }

	END_TIMER;


	return 1;
}


int Call_win::quick_recycle() {

	n_qrecycle++;
	INIT_TIMER;
	START_TIMER;

        std::string::size_type i;
        for(i=0; i != read_len ; i++) {
                sites[i].pos = sites[i+win_size].pos;
                sites[i].ori = sites[i+win_size].ori;
                sites[i].depth = sites[i+win_size].depth;
                sites[i].repeat_time = sites[i+win_size].repeat_time;
                sites[i].dep_uni = sites[i+win_size].dep_uni;
        }
        for(i=read_len; i != read_len+win_size; i++) {
                sites[i].ori = 0xFF;
                sites[i].pos = sites[i-1].pos+1;
                sites[i].depth = 0;
                sites[i].repeat_time = 0;
                sites[i].dep_uni = 0;
        }
        for(i=0; i < read_len + win_size; i++){
             	depthArr[i] = sites[i].depth;
             	repeat_timeArr[i] = sites[i].repeat_time;
             	dep_uniArr[i] = sites[i].dep_uni;
        }

        memcpy(numNonZeroPerSite,numNonZeroPerSite + win_size,sizeof(uint32)*read_len);
        memset(numNonZeroPerSite + read_len, 0, sizeof(uint32)*win_size);

	uint64 numSite = win_size + read_len;
        for(uint32 siteID = 0; siteID < read_len; siteID++){
             //for(uint32 offset = 0; offset < numNonZeroPerSite[siteID]; offset++){
             for(uint32 offset = 0; offset < MAX_NUM_NON_ZERO; offset++){
                     //sparse_base_infoArr[offset*numSite + siteID] = sparse_base_infoArr[offset*numSite + siteID + global_win_size];
                     sparse_base_infoArr[BASEINFO_IDX(siteID, offset, numSite, MAX_NUM_NON_ZERO)] =
                                sparse_base_infoArr[BASEINFO_IDX(siteID + global_win_size, offset, numSite, MAX_NUM_NON_ZERO)];
             }
        }


	END_TIMER;

        return 1;
}


int Call_win::new_call_cns(Chr_name call_name, Chr_info* call_chr,
                ubit64_t call_length, Prob_matrix * mat, Parameter * para,
                std::ofstream & consensus, std::ofstream & baseinfo) {
	
	return call_cns2(call_name, call_chr, call_length, mat, para, consensus, baseinfo);
}


/**
 * this is a modified version
 */

static unsigned long winId = 0;

int Call_win::call_cns2(Chr_name call_name, Chr_info* call_chr,
		ubit64_t call_length, Prob_matrix * mat, Parameter * para,
		std::ofstream & consensus, std::ofstream & baseinfo) {


	double /*rank_sum_test_value, */binomial_test_value;
	rate_t* type_likely_buf = likelihoodData->typeLikelihood;
	INIT_TIMER;
	n_cns++;

	/*1. memory prepare and initialization*/
	for (std::string::size_type j = 0; j != call_length; j++) {
		oriArr[j] = (call_chr->get_bin_base(posArr[j])) & 0xF;
		if (((oriArr[j] & 4) != 0)&&depthArr[j] == 0) { // an N
			likelihoodData->nflag[j] = true;
			continue;
		} else {
			likelihoodData->nflag[j] = false;
		}
	}

	/* 2. likelihood calculation  */
	likelihoodData->prepare(call_length);	
	
	if(maxNumNonZero > 0) {
#ifdef CPU_COUNTING
                copy_likelihoodSparse(sparse_base_infoArr, sizeof(uint32)*(win_size + read_len)*MAX_NUM_NON_ZERO,
                                numNonZeroPerSite, call_length);
#endif
        }
	cuda_likelihoodSparse(likelihoodData->typeLikelihood,
                              call_length, maxNumNonZero, win_size + read_len, MAX_NUM_NON_ZERO); 

	




	/* 3. calculate postier probability */
        for(int siteID = 0; siteID < call_length; siteID++)
        {
                if ((oriArr[siteID] & 0x8) && para->refine_mode) {
                        Snp_info* snp = call_chr->find_snp(posArr[siteID]);

                        snpData->validatedArr[siteID] = call_chr->find_snp(posArr[siteID])->is_validated();
                        snpData->hapmap_siteArr[siteID] = call_chr->find_snp(posArr[siteID])->is_hapmap();
                        snpData->indel_siteArr[siteID] = call_chr->find_snp(posArr[siteID])->is_indel();
                        for(char ii = 0; ii != 4; ii++){
                                snpData->freqArr[siteID * 4 + ii] = call_chr->find_snp(posArr[siteID])->get_freq(ii);
                        }
                }
        }
#ifdef CPU_POSTERIOR
        gold_posterior(
                        likelihoodData,
                        mat->p_prior,
                        mat->type_prob,//(rate_t*)malloc(sizeof(rate_t) * LIKELY_SIZE)
                        mat->p_rank,//(rate_t*)malloc(sizeof(rate_t) * PRIOR_SIZE)
                        sparse_base_infoArr,
                        numNonZeroPerSite,
                        count_uni_buf,
                        q_sumArr,
                        dep_uniArr,
                        count_allArr,
                        oriArr,
                        posArr,
                        base1,
                        base2,
                        type1,
                        q_cns,
                        rank_sum_test_value,

                        snpData->validatedArr, snpData->hapmap_siteArr, snpData->indel_siteArr, snpData->freqArr,

                        para->het_val_r, para->het_unval_r, para->althom_val_r, para->althom_unval_r, para->rank_sum_mode, para->is_monoploid, para->refine_mode, para->q_max - para->q_min,

                        read_len,
                        call_length
        );
#else
        cuda_posterior(
                likelihoodData,
                oriArr,

                base1,
                base2,
                type1,
                q_cns,
                rank_sum_test_value,

                snpData->validatedArr, snpData->hapmap_siteArr, snpData->indel_siteArr, snpData->freqArr,

                para->het_val_r, para->het_unval_r, para->althom_val_r, para->althom_unval_r, para->rank_sum_mode, para->is_monoploid, para->refine_mode, para->q_max - para->q_min,

                read_len,
                win_size,
                call_length,
                maxNumNonZero
        );
#endif

if(useCompression) {
	/*4.1 compression*/
        //column 1 and 2
        c_1_2->push_back(call_name, posArr[0] + 1, call_length);

        //column 3, 4, and 6
        char* col_3 = (char*)buffer->getBuf(sizeof(char)*call_length);
        char* col_4 = (char*)buffer->getBuf(sizeof(char)*call_length);
        char* col_6 = (char*)buffer->getBuf(sizeof(char)*call_length);
        extract_3_4_6(col_3, col_4, col_6, likelihoodData->nflag, oriArr, type1, base1, call_length);
        c_3_4_6->compress(col_3, col_4, col_6, likelihoodData->nflag, call_length);
        c_3_4_6->write();
        buffer->clear();


        //column 5, 7, 8, 9, 14, 16 
        c_5_7_8_9_14_16->compress_write(q_cns, q_sumArr, count_uni_buf, count_allArr,
                                        depthArr, repeat_timeArr,
                                        base1, likelihoodData->nflag, call_length);


        //column 10, encoding
        char* col_10 = (char*)buffer->getBuf(sizeof(char)*call_length);
        extract_10(col_10, base1, base2, likelihoodData->nflag, call_length);
        c_10->compress_write(col_10, likelihoodData->nflag, call_length);
        buffer->clear();


        //column 11, 12, 13
        int* col11 = (int*)buffer->getBuf(sizeof(int)*call_length);
        int* col12 = (int*)buffer->getBuf(sizeof(int)*call_length);
        int* col13 = (int*)buffer->getBuf(sizeof(int)*call_length);
        extract_11_12_13(col11, col12, col13, q_sumArr, count_uni_buf,
                        count_allArr, base1, base2, likelihoodData->nflag, call_length);
        c_11_12_13->compress_write(col11, col12, col13, likelihoodData->nflag, call_length);
        buffer->clear();


        //column 15, 17
        c_15_17->compress_write(base1, base2, likelihoodData->nflag, rank_sum_test_value, oriArr, call_length);

} else {
	/* 4.2 the original output*/
	for (std::string::size_type j = 0; j != call_length; j++) {

	
		if (likelihoodData->nflag[j]) { 
                                consensus << call_name << '\t' << (posArr[j] + 1)
                                                << "\tN\tN\t0\tN\t0\t0\t0\tN\t0\t0\t0\t0\t1.000\t255.000\t0"
                                                << endl;
                
				continue;
		} else if (base1[j] < 4 && base2[j] < 4) {
                                        consensus << call_name << '\t' << (posArr[j] + 1)
                                                        << '\t' << ("ACTGNNNN"[(oriArr[j] & 0x7)])
                                                        << '\t' << abbv[type1[j]] << '\t' << q_cns[j] << '\t'
                                                        << ("ACTGNNNN"[base1[j]]) << '\t'
                                                        << (q_sumArr[j * 4 + base1[j]] == 0 ? 0
                                                                        : q_sumArr[j * 4 + base1[j]]
                                                                                        / count_uni_buf[j * 4 + base1[j]])
                                                        << '\t' << count_uni_buf[j * 4 + base1[j]] << '\t'
                                                        << count_allArr[j * 4 + base1[j]] << '\t'
                                                        << ("ACTGNNNN"[base2[j]]) << '\t'
                                                        << (q_sumArr[j * 4 + base2[j]] == 0 ? 0
                                                                        : q_sumArr[j * 4 + base2[j]]
                                                                                        / count_uni_buf[j * 4 + base2[j]])
                                                        << '\t' << count_uni_buf[j * 4 + base2[j]] << '\t'
                                                        << count_allArr[j * 4 + base2[j]] << '\t'
                                                        << depthArr[j] << '\t' << showpoint
                                                        << rank_sum_test_value[j] << '\t' << (depthArr[j]
                                                        == 0 ? 255 : (double) (repeat_timeArr[j])
                                                        / depthArr[j]) << '\t'
                                                        << ((oriArr[j] & 8) ? 1 : 0) << endl;

		} else if (base1[j] < 4) {
                                        consensus << call_name << '\t' << (posArr[j] + 1)
                                                        << '\t' << ("ACTGNNNN"[(oriArr[j] & 0x7)])
                                                        << '\t' << abbv[type1[j]] << '\t' << q_cns[j] << '\t'
                                                        << ("ACTGNNNN"[base1[j]]) << '\t'
                                                        << (q_sumArr[j * 4 + base1[j]] == 0 ? 0
                                                                        : q_sumArr[j * 4 + base1[j]]
                                                                                        / count_uni_buf[j * 4 + base1[j]])
                                                        << '\t' << count_uni_buf[j * 4 + base1[j]] << '\t'
                                                        << count_allArr[j * 4 + base1[j]] << '\t'
                                                        << "N\t0\t0\t0\t" << depthArr[j] << '\t'
                                                        << showpoint << rank_sum_test_value[j] << '\t'
                                                        << (depthArr[j] == 0 ? 255
                                                                        : (double) (repeat_timeArr[j])
                                                                                        / depthArr[j]) << '\t'
                                                        << ((oriArr[j] & 8) ? 1 : 0) << endl;

		} else {
                                        consensus << call_name << '\t' << (posArr[j] + 1)
                                                        << "\tN\tN\t0\tN\t0\t0\t0\tN\t0\t0\t0\t0\t1.000\t255.000\t0"
                                                        << endl;

		}

	}
}


	return 1;
}






int Call_win::new_soap2cns(std::ifstream & alignment, std::ofstream & consensus,
		std::ofstream & baseinfo, Genome * genome, Prob_matrix * mat,
		Parameter * para, const string cnsFileName) {
	int g_win = 0;

	/* data and memory buffer initialization */
        GSNPKernelInit(win_size, para->read_length, para->pcr_dependency, para->global_dependency, mat->p_matrix, mat->p_prior, mat->p_rank);
	likelihoodData = new st_likelihood(para, win_size, mat->p_matrix, count_uni_buf);	
	snpData = new st_snp(win_size);
	pcr_dep_count = new int[para->read_length * 2];
	base1 = (char*)malloc(sizeof(char)*win_size);
	base2 = (char*)malloc(sizeof(char)*win_size);
	type1 = (char*)malloc(sizeof(char)*win_size);
	q_cns = (int*)malloc(sizeof(int)*win_size);
	rank_sum_test_value = (double*)malloc(sizeof(double)*win_size);

		
	/* for compression */
	buffer = new Buffer(1024*1024*1024);	//initialized as 1 GB memory
	c_1_2 = new Column_1_2();
	c_3_4_6 = new Column_3_4_6(global_win_size);
	c_3_4_6->startWrite(cnsFileName);
	c_5_7_8_9_14_16 = new Column_RLEDict(global_win_size);
	c_5_7_8_9_14_16->startWrite(cnsFileName);
	c_10 = new Column_10(global_win_size);
	c_10->startWrite(cnsFileName);
	c_11_12_13 = new Column_11_12_13(global_win_size);
	c_11_12_13->startWrite(cnsFileName);
	c_15_17 = new Column_15_17(global_win_size);	
	c_15_17->startWrite(cnsFileName);
	dictComprStart(global_win_size, sizeof(double));


	Soap_format soap;
	map<Chr_name, Chr_info*>::iterator current_chr, prev_chr;
	current_chr = prev_chr = genome->chromosomes.end();
	int coord, sub;
	int last_start(0);
	bool recycled(false);
	unsigned long count = 0;
	uint64 readCount = 0;
        struct st_update* updateData = new st_update(para);
        bool mark;


	SoapData soapData(para->read_length, READ);
	const int readLen = para->read_length;


	printf("Data preprocessing done, start window-by-window processing ...\n");

	//while(!feof(hitFile)){	
	while(true) {
 
		if(!soapData.fetch(updateData->readArr + readCount*readLen,
				 updateData->qualArr + readCount*readLen,
				 (updateData->hitArr + readCount), 
				 (updateData->read_lenArr + readCount),
				 (updateData->positionArr + readCount),
				 (updateData->mismatchArr + readCount), 
				 (updateData->strandArr + readCount))) {
			break;
		}


        	updateData->readLen_hitArr[readCount] = (int)(
 	        	                                ((updateData->read_lenArr[readCount]) & 0xffff)<<16
                        	                        |(updateData->hitArr[readCount])
       	                         	                );

                count++;

	
                        if (updateData->positionArr[readCount] < 0) {
                                continue;
                        }
                        if (current_chr == genome->chromosomes.end()) {
			
                                current_chr = genome->chromosomes.find(mat->chr_name);
#ifdef CPU_RECYCLE
                                initialize(0);
#else
                                cuda_init(win_size, read_len, 0, posArr);
#endif
                                last_start = 0;
                        } else {
                                ;
                        }


                        if (updateData->positionArr[readCount] + updateData->read_lenArr[readCount] >= current_chr->second->length()) {
                                continue;
                        }

                        if (updateData->positionArr[readCount] < last_start) {
                                cerr << "Errors in sorting:" << updateData->positionArr[readCount] << "<"  << last_start << endl;
                                exit(255);
                        }


                        recycled = false;
                        unsigned long nrecycle = 0;
                        mark = false;
                        if(( readCount!=0 ) && (updateData->positionArr[readCount] /win_size >last_start /win_size))
                        {
				g_win++;
                               	assert(readCount <= _read_num_per_window);
	                        gold_pos_match(updateData->positionArr, updateData->siteMatchArr, readCount);
        	                maxNumNonZero = cuda_updateSparse(updateData, sparse_base_infoArr, numNonZeroPerSite, read_len, count_uni_buf, q_sumArr, count_allArr, depthArr, repeat_timeArr, dep_uniArr, readCount);
                                mark = true;
                        }

                        while (updateData->positionArr[readCount] / win_size > last_start / win_size) {
                                if (!para->region_only  || current_chr->second->is_in_region_win(last_start)) {

                                        if (last_start > posArr[win_size - 1]) {
                                                cuda_init(win_size, read_len, (last_start / win_size) * win_size, posArr);
                                        }
                                        new_call_cns(current_chr->first, current_chr->second,
                                                        win_size, mat, para, consensus, baseinfo);

                                        last_start = posArr[win_size - 1];
                                        recycled = false;
                                }

                                if (!para->region_only) {
                                        if(nrecycle < 2) {
                                                cuda_recycle(win_size, read_len, posArr);
                                                n_recycle++;
                                        }
                                        else {
						cuda_recycle(win_size, read_len, posArr);
                                                n_qrecycle++;
                                         }

                                        recycled = true;
                                        last_start = posArr[win_size - 1];
                                } else {
                                        ERROR_EXIT("Never here 2!");
                                }

                                nrecycle++;
                        }//end while

                        last_start = updateData->positionArr[readCount];
                        if(mark){
				int tmpVal = readCount;
                                readCount = 0;
                                updateData -> updateSoap2(readCount, tmpVal);
                        }
                        readCount +=1;
//                }//end if

                                        
	}//end while

        if(!mark){
		assert(readCount <= _read_num_per_window);

                gold_pos_match(updateData->positionArr, updateData->siteMatchArr, readCount);
		maxNumNonZero = cuda_updateSparse(updateData, sparse_base_infoArr, numNonZeroPerSite, 
			read_len, count_uni_buf, q_sumArr, count_allArr, depthArr, repeat_timeArr, dep_uniArr, readCount);

        }


	//The unprocessed tail of chromosome
	printf("Processing the tail ...\n");
	recycled = false;
	unsigned long nrecycle = 0;
	while (current_chr->second->length() > posArr[win_size - 1]) {
                if (!para->region_only || current_chr->second->is_in_region_win(last_start)) {
                        if (last_start > posArr[win_size - 1]) {
#ifdef CPU_RECYCLE
                                initialize((last_start / win_size) * win_size);
#else
                                cuda_init(win_size, read_len, (last_start / win_size) * win_size, posArr);
#endif
                        }

                       	new_call_cns(current_chr->first, current_chr->second, win_size, mat,
                       	        para, consensus, baseinfo);
			
                        last_start = posArr[win_size - 1];
                        recycled = false;
                }
                if (!para->region_only) {
                        if(nrecycle < 2) {
#ifdef CPU_RECYCLE
                                recycle();
#else
                                cuda_quickRecycle(win_size, read_len, posArr, count_uni_buf, q_sumArr, count_allArr, depthArr, repeat_timeArr);
                                n_recycle++;
#endif
                        }
                        else {
#ifdef CPU_RECYCLE
                                quick_recycle();
#else
				cuda_quickRecycle(win_size, read_len, posArr, count_uni_buf, q_sumArr, count_allArr, depthArr, repeat_timeArr);          
                                n_qrecycle++;
#endif
                        }

                        recycled = true;
                        last_start = posArr[win_size - 1];
                } else {
                        ERROR_EXIT("Never here 3!");
                }

                nrecycle++;
		
	}
	// The last window
	if (last_start > posArr[win_size - 1]) {
#ifdef CPU_RECYCLE
                initialize((last_start / win_size) * win_size);
#else
                cuda_init(win_size, read_len, (last_start / win_size) * win_size, posArr);
#endif

	}
	new_call_cns(current_chr->first, current_chr->second, current_chr->second->length() % win_size, mat, para, consensus, baseinfo);


	/* compressed data output and close */
	c_1_2->write(cnsFileName);
	c_1_2->close();
	c_3_4_6->endWrite();
	c_3_4_6->close();
	c_5_7_8_9_14_16->endWrite();
	c_5_7_8_9_14_16->close();
	c_10->endWrite();
	c_10->close();
	c_11_12_13->endWrite();
	c_11_12_13->close();
	c_15_17->endWrite();
	c_15_17->close();
	dictComprEnd();
	soapData.close();


	alignment.close();
	consensus.close();
	baseinfo.close();

	return 1;
}
