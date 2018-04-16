#ifndef __COMPRESS_H__
#define __COMPRESS_H__


#include "soap_snp.h"
#include <gsnp_util.h>
#include <vector>
#include <string>

#include "dict_compression.h"

using namespace std;

/*
typedef uint32 code_t;

uint64 getNumBinWord(const uint64 numElement);
uint64 encode(code_t* code, const char* seq, const uint64 numElement);
void decode(char* seq, const uint64 numElement, const code_t* code, const uint64 numCodeWord);
void getNRegion(vector< pair<int, int> >* nRegion, const bool* nflag,
                const char* base1, const char* base2, const uint64 call_length);
*/


/*GPU kernels*/
void dictComprStart(const uint64 numElement, const uint64 dataSize);
void dictComprEnd();





struct Buffer {
	byte* buf;
	uint64 maxBufSize;
	uint64 usedBufSize;

	Buffer(const uint64 in_maxBufSize) {
		maxBufSize = in_maxBufSize;
		buf = (byte*)malloc(maxBufSize);
		usedBufSize = 0;
	}

        byte* getBuf(const uint64 sizeInByte) {
       
		if((usedBufSize + sizeInByte) > maxBufSize) {
			byte* tempBuf = (byte*)malloc(usedBufSize);
			memcpy(tempBuf, buf, usedBufSize);
			free(buf);

			maxBufSize = (usedBufSize + sizeInByte)*2;
			buf = (byte*)malloc(maxBufSize);
			memcpy(buf, tempBuf, usedBufSize);
			free(tempBuf);

                        byte* p = buf + usedBufSize;
                        usedBufSize += sizeInByte;
                        return p;			
		} else {

                	byte* p = buf + usedBufSize;
                	usedBufSize += sizeInByte;
                	return p;
		}
        };

	//dump the useful data into the disk
	void dumpToDisk(FILE* file) {
		safe_fwrite(buf, sizeof(byte), usedBufSize, file);
	}

	void clear() {
		usedBufSize = 0;
	}

	void close() {
		if(buf) {
			free(buf);
			buf = NULL;
		}
		maxBufSize = 0;
		usedBufSize = 0;
	}
};

inline
void extract_3_4_6(char* col3, char* col4, char* col6,
		   bool* nflag, const char* ori, char* type1, char* base1, const uint64 winSize) {
	for (int j = 0; j < winSize; j++) {
		if(nflag[j]) { 
			col3[j] = 'N';
			col4[j] = 'N';
			col6[j] = 'N';               
			continue;
		} else if (base1[j] < 4) {
                        col3[j] = ("ACTGNNNN"[(ori[j] & 0x7)]);
                        col4[j] = abbv[type1[j]];
                        col6[j] = ("ACTGNNNN"[base1[j]]);
		} else {
                        col3[j] = 'N';
                        col4[j] = 'N';
                        col6[j] = 'N';
		}
	}
}

struct Column_1_2 {
	uint64 maxNumWindow;
	int maxCallNameSize;

	uint64 numWindow;
	char** callName;
	char* callNameBuf;
	uint64* posStart;
	uint64* winSize;

	Column_1_2(const uint64 in_maxNumWindow = 1024*1024, const int in_maxCallNameSize = 128) {
		maxNumWindow = in_maxNumWindow;
		maxCallNameSize = in_maxCallNameSize;
		numWindow = 0;

		callName = (char**)malloc(sizeof(char*)*maxNumWindow);
		callNameBuf = (char*)malloc(sizeof(char)*maxCallNameSize*maxNumWindow);
		for(uint64 i = 0; i < maxNumWindow; i++) {
			callName[i] = callNameBuf + maxCallNameSize*i;
		}
		posStart = (uint64*)malloc(sizeof(uint64)*maxNumWindow);
		winSize = (uint64*)malloc(sizeof(uint64)*maxNumWindow);	
	}

	void close() {
		free(callName);
		free(callNameBuf);
		free(posStart);
		free(winSize);
		maxNumWindow = 0;
		winSize = 0;
	}

	void push_back(string in_callName, uint64 in_posStart, uint64 in_winSize) {

		strcpy(callName[numWindow], in_callName.c_str());
		posStart[numWindow] = in_posStart;
		winSize[numWindow] = in_winSize;

		numWindow++;
		assert(numWindow <= maxNumWindow);
	}


	void write(string cnsFileName) {
		FILE* file = NULL;
		safe_fopen(&file, (cnsFileName + ".1_2").c_str(), "w");

		/* 1. write some meta information */
		fprintf(file, "%lu\n", numWindow);
			
		/*2. each line three items per window*/
		for(uint64 i = 0; i < numWindow; i++) {
			fprintf(file, "%s\t%lu\t%lu\n", 
				callName[i], posStart[i], winSize[i]);	
		}

		safe_fclose(file);
	}

	void read(const string cnsFileName) {
                FILE* file = NULL;
                safe_fopen(&file, (cnsFileName + ".1_2").c_str(), "r");

		fscanf(file, "%lu", &numWindow);
	
		for(uint64 i = 0; i < numWindow; i++) {
			fscanf(file, "%s\t%lu\t%lu\n", (callName[i]), &(posStart[i]), &(winSize[i]));
		}

		safe_fclose(file);
	}


	//decompress one window of data for the column 1 and 2
	void decompress() {

	}
};

struct Column_3_4_6 {

	//structure parameters
	uint64 maxNumElement;
	bool bufferInitialized;
	FILE* file;

	/**
 	* helper index to access param
 	* */
	static const int numParam = 7;
	static const int __isAllN = 0;		//1 yes, all N
	static const int __dataSize = 1;	//the total memory size allocated
	static const int __winSize = 2;
	static  const int __col_3_codeLen = 3;	//number of codes
	static const int __col_3_N_len = 4;	//number of N
	static const int __col_4_len = 5;	//number of col 4 bases differs from col 3
	static const int __col_6_len = 6;	//number of col 6 bases differs from col 3

	/**
 	* output parameters
 	* */
	uint64 param[7]; 
	Buffer* buf;

	//temp pointer
	int* col_4_id;
	char* col_4_val;
	int* col_6_id;
	char* col_6_val;
	int* col_3_N_id;
	code_t* col_3_code;

	Column_3_4_6(const uint64 in_maxNumElement) {
		maxNumElement = in_maxNumElement;
		uint64 bufSize = sizeof(int)*maxNumElement*3 + 	//idx for column 4 and 6, idx for N in column 3
			  sizeof(char)*maxNumElement*2 + 	//base for column 4 and 6
			  sizeof(code_t)*maxNumElement;		//code for column 3
		buf = new Buffer(bufSize);
	
	};

	void close() {
		buf->close();
	}

	void startWrite(const string cnsFileName) {
		safe_fopen(&file, string(cnsFileName + ".3_4_6").c_str(), "wb");
	}

	void endWrite() {
		safe_fclose(file);
	}


	void startRead(const string cnsFileName) {
		safe_fopen(&file, string(cnsFileName + ".3_4_6").c_str(), "rb");	
	}

	void endRead() {
		safe_fclose(file);
	}

	//read a block of window data
	void read() {
		clear();

		/*1. read the array of parameters*/
		safe_fread(param, sizeof(uint64), numParam, file);

		/*2. assign the memory pointer*/	
		assignMemory();

		/*3. read the real data*/
		safe_fread(buf->buf, sizeof(byte), param[__dataSize], file);
	}

	//write one block data, and clean the data stored in the structure
	void write() {
		/*1. write the parameter array*/
		safe_fwrite(param, sizeof(uint64), numParam, file);

		/* 2. write the data content*/
		buf->dumpToDisk(file);		

		/* 3. clean the data content for the next window */
		clear();
	}

	void clear() {
       		param[__isAllN] = 1;
        	param[__dataSize] = 0;
        	param[__winSize] = 0;
        	param[__col_3_codeLen] = 0;
        	param[__col_3_N_len] = 0;
        	param[__col_4_len] = 0;
        	param[__col_6_len] = 0;
	
		buf->clear();
	}

	void decompress(char* col_3, char* col_4, char* col_6) {

		if(param[__isAllN] == 1) {
			for(uint64 i = 0; i < param[__winSize]; i++) {
				col_3[i] = 'N';
				col_4[i] = 'N';
				col_6[i] = 'N';
			}
		} else {
			const uint64 winSize = param[__winSize];	

			/* 1. decode for column 3 */	
			decode(col_3, winSize, col_3_code, param[__col_3_codeLen]);

			/* 2. set N values */
			for(uint64 i = 0; i < param[__col_3_N_len]; i++) {
				col_3[col_3_N_id[i]] = 'N';
			}

			/*3. set differences of column 4*/
			memcpy(col_4, col_3, sizeof(char)*param[__winSize]);
			for(uint64 i = 0; i < param[__col_4_len]; i++) {
				col_4[col_4_id[i]] = col_4_val[i];
			}	
	
			/*4. set differences of column 6*/
			memcpy(col_6, col_3, sizeof(char)*param[__winSize]);
			for(uint64 i = 0; i < param[__col_6_len]; i++) {
				col_6[col_6_id[i]] = col_6_val[i];
			}
		}
	}

	void assignMemory() {
                col_4_id = (int*)buf->getBuf(sizeof(int)*param[__col_4_len]);
                col_4_val = (char*)buf->getBuf(sizeof(char)*param[__col_4_len]);
                col_6_id = (int*)buf->getBuf(sizeof(int)*param[__col_6_len]);
                col_6_val = (char*)buf->getBuf(sizeof(char)*param[__col_6_len]);
                col_3_N_id = (int*)buf->getBuf(sizeof(int)*param[__col_3_N_len]);
                col_3_code = (code_t*)buf->getBuf(sizeof(code_t)*param[__col_3_codeLen]);
	}

	void compress(char* col_3, char* col_4, char* col_6, const bool* nflag,
			const uint64 numElement) {
		
				
		/*0. check parameters */
		if(isAllN(nflag, numElement)) {
			param[__isAllN] = 1;
			param[__dataSize] = 0;
			param[__winSize] = numElement;
			param[__col_3_codeLen] = 0;
			param[__col_3_N_len] = 0;
			param[__col_4_len] = 0;
			param[__col_6_len] = 0;
			//printf("All N, skip the compression ...\n");
		
			return;
		} else {
			param[__isAllN] = 0;
		}
		assert(numElement <= maxNumElement);
		param[__winSize] = numElement;

		/*1. find the size of col 4 and col 6*/
		uint64 col_3_N_len = 0;
		uint64 col_4_len = 0;
		uint64 col_6_len = 0;
		for(uint64 i = 0; i < numElement; i++) {
			if(col_3[i] == 'N') {
				col_3_N_len++;
			}
			if(col_4[i] != col_3[i]) {
				col_4_len++;
			}
			if(col_6[i] != col_3[i]) {
				col_6_len++;
			}
		}
		param[__col_4_len] = col_4_len;
		param[__col_6_len] = col_6_len;
		param[__col_3_N_len] = col_3_N_len;
		param[__col_3_codeLen] = getNumBinWord(numElement);

		/* 2. assigne memory */
		assignMemory();
		param[__dataSize] = buf->usedBufSize;

		/* 3. replace N with A for col 3, store differences for col 4 and 6 */
		uint64 col_3_N_idx = 0;
		uint64 col_4_idx = 0;
		uint64 col_6_idx = 0;
		for(uint64 i = 0; i < numElement; i++) {
			if(col_4[i] != col_3[i]) {
				col_4_id[col_4_idx] = i;
				col_4_val[col_4_idx] = col_4[i];
				col_4_idx++;
			}
			if(col_6[i] != col_3[i]) {
				col_6_id[col_6_idx] = i;
				col_6_val[col_6_idx] = col_6[i];
				col_6_idx++;
			}

                        if(col_3[i] == 'N') {
                                col_3_N_id[col_3_N_idx] = i;
				col_3[i] = 'A';
                                col_3_N_idx++;
                        }
		}
		assert(col_3_N_idx == param[__col_3_N_len]);
		assert(col_4_idx == param[__col_4_len]);
		assert(col_6_idx == param[__col_6_len]);

		/*4. encode*/
		encode(col_3_code, col_3, numElement);
	}

	uint64 getDataSize() {
		return buf->usedBufSize;
	}
};



inline void extract_5(int* col5, const int* q_cns, const char* base1, const bool* nflag, const int call_length) {
	for (int j = 0; j != call_length; j++) {
		if (nflag[j]) { 
			col5[j] = 0;               
			continue;
		} else if (base1[j] < 4) {
			col5[j] = q_cns[j];
		} else {
			col5[j] = 0;
		}
	}
}

inline void extract_7(int* col7, const int* q_sumArr, const int* count_uni_buf,
			 const char* base1, const bool* nflag, const int call_length) {
	for (int j = 0; j != call_length; j++) {
		if (nflag[j]) { 
			col7[j] = 0;
      			continue;
		} else if (base1[j] < 4) {
			//col7[j] = (sites[j].q_sum[base1[j]] == 0 ? 0 : sites[j].q_sum[base1[j]]/ sites[j].count_uni[base1[j]]);
                        col7[j] = (q_sumArr[j * 4 + base1[j]] == 0 ? 0 : q_sumArr[j * 4 + base1[j]]/count_uni_buf[j * 4 + base1[j]]);
		} else {
			col7[j] = 0;
		}
	}
}

inline void extract_8(int* col8, const int* count_uni_buf, const char* base1, const bool* nflag, const int call_length) {
        for (int j = 0; j != call_length; j++) {
                if (nflag[j]) {                                 
                        col8[j] = 0;
                        continue;
                } else if (base1[j] < 4) {
			//col8[j] = sites[j].count_uni[base1[j]];
			col8[j] = count_uni_buf[j * 4 + base1[j]];
                } else {
                        col8[j] = 0;
                }
        }
}

inline void extract_9(int* col9, const int* count_allArr, const char* base1, const bool* nflag, const int call_length) {
        for (int j = 0; j != call_length; j++) {
                if (nflag[j]) {
                        col9[j] = 0;
                        continue;
                } else if (base1[j] < 4) {
			//col9[j] = sites[j].count_all[base1[j]];
			col9[j] = count_allArr[j * 4 + base1[j]];
                } else {
                        col9[j] = 0;
                }
        }
}

inline void extract_14(int* col14, const int* depthArr, const char* base1, const bool* nflag, const int call_length) {
        for (int j = 0; j != call_length; j++) {
                if (nflag[j]) {
                        col14[j] = 0;
                        continue;
                } else if (base1[j] < 4) {
                        //col14[j] = sites[j].depth;
                        col14[j] = depthArr[j];
                } else {
                        col14[j] = 0;
                }
        }
}

inline void extract_16(double* col16, const int* depthArr, const int* repeat_timeArr,
			const char* base1, const bool* nflag, const int call_length) {
        for (int j = 0; j != call_length; j++) {
                if (nflag[j]) {
                        col16[j] = 255.0;
                        continue;
                } else if (base1[j] < 4) {
                        //col16[j] = (sites[j].depth == 0 ? 255 : (double) (sites[j].repeat_time) / sites[j].depth);
                	col16[j] = (depthArr[j] == 0 ? 255 : (double) (repeat_timeArr[j]) / depthArr[j]);
                } else {
                        col16[j] = 255.0;
                }
        }
}



//5, 7, 8, 9, 14, 16
struct Column_RLEDict {
	
	uint64 maxWinSize;

	RLEDictData<int>* RLEDictInt;
	RLEDictData<double>* RLEDictDbl;

	//for extracted data
	byte* tmpBuf;
	uint64 tmpBufSize; 
	int* intBuf;
	double* dblBuf;

	//for one column of the compressed data before writing
	Buffer* cData;

	FILE* file;	

	static const int numParam = 6;
	uint64 param[6];
	static const int __numRun = 0;
	static const int __numWordValue = 1;
	static const int __numWordLength = 2;
	static const int __numEleValueDict = 3;
	static const int __numEleLengthDict = 4;
	static const int __dataSize = 5;	//total data size allocated

	Column_RLEDict(const uint64 in_maxWinSize) {
		maxWinSize = in_maxWinSize;

		//RLEDictInt = new RLEDictData<int>(maxWinSize);
		//RLEDictDbl = new RLEDictData<double>(maxWinSize);
		tmpBufSize = sizeof(int)*maxWinSize + sizeof(double)*maxWinSize;
		tmpBuf = (byte*)malloc(tmpBufSize);
		intBuf = (int*)tmpBuf;
		dblBuf = (double*)(tmpBuf + sizeof(int)*maxWinSize);

		cData = new Buffer(sizeof(double)*maxWinSize);
	};

	void clear() {
		cData->clear();
	}

	void close() {
		//RLEDictInt->close();
		//RLEDictDbl->close();
		//free(intBuf);
		//free(dblBuf);
	}

        void startWrite(const string cnsFileName) {
                RLEDictInt = new RLEDictData<int>(maxWinSize, true);
                RLEDictDbl = new RLEDictData<double>(maxWinSize, true);

                safe_fopen(&file, string(cnsFileName + ".5_7_8_9_14_16").c_str(), "wb");
        }

        void endWrite() {
                safe_fclose(file);
        }


        void startRead(const string cnsFileName) {
                RLEDictInt = new RLEDictData<int>(maxWinSize, false);
                RLEDictDbl = new RLEDictData<double>(maxWinSize, false);

                safe_fopen(&file, string(cnsFileName + ".5_7_8_9_14_16").c_str(), "rb");
        }

        void endRead() {
                safe_fclose(file);
        }	

	//compress and write the compressed data column by column
	void compress_write(const int* q_cns, const int* q_sumArr, const int* count_uni_buf, 
			const int* count_allArr, const int* depthArr, const int* repeat_timeArr, 
			const char* base1, const bool* nflag, const int call_length) {


		//INIT_TIMER;

		//first write a global N indicator
		int NIndicator = 0;
		if(isAllN(nflag, call_length)) {
			NIndicator = 1;
			fwrite(&NIndicator, sizeof(int), 1, file);
			return ;
		} else {
                        NIndicator = 0;
                        fwrite(&NIndicator, sizeof(int), 1, file);
		}

		assert(call_length <= maxWinSize);

		////5, 7, 8, 9, 14, 16

		//column 5
		int* col5 = intBuf;
		extract_5(col5, q_cns, base1, nflag, call_length);
		//START_TIMER;
			RLEDictInt->compress(col5, call_length);
		//END_TIMER;
		//PRINT_TIMER_SEC("column 5");
		copyFromData(RLEDictInt);
		safe_fwrite(param, sizeof(uint64), numParam, file);
		cData->dumpToDisk(file);
		clear();

		//column 7
                int* col7 = intBuf;
		extract_7(col7, q_sumArr, count_uni_buf, base1, nflag, call_length);
		//START_TIMER;
                	RLEDictInt->compress(col7, call_length);
		//END_TIMER;
		//PRINT_TIMER_SEC("column 7");
                copyFromData(RLEDictInt);
                safe_fwrite(param, sizeof(uint64), numParam, file);
                cData->dumpToDisk(file);
                clear();


		//column 8
                int* col8 = intBuf;
                extract_8(col8, count_uni_buf, base1, nflag, call_length);
		//START_TIMER
                	RLEDictInt->compress(col8, call_length);
		//END_TIMER;
		//PRINT_TIMER_SEC("column 8");
                copyFromData(RLEDictInt);
                safe_fwrite(param, sizeof(uint64), numParam, file);
                cData->dumpToDisk(file);
                clear();

		//column 9
                int* col9 = intBuf;
                extract_9(col9, count_allArr, base1, nflag, call_length);
		//START_TIMER;
                	RLEDictInt->compress(col9, call_length);
		//END_TIMER;
		//PRINT_TIMER_SEC("column 9");
                copyFromData(RLEDictInt);
                safe_fwrite(param, sizeof(uint64), numParam, file);
                cData->dumpToDisk(file);
                clear();

		//column 14
                int* col14 = intBuf;
                extract_14(col14, depthArr, base1, nflag, call_length);
		//START_TIMER;
                	RLEDictInt->compress(col14, call_length);
		//END_TIMER;
		//PRINT_TIMER_SEC("column 14");
                copyFromData(RLEDictInt);
                safe_fwrite(param, sizeof(uint64), numParam, file);
                cData->dumpToDisk(file);
                clear();

		//column 16
                double* col16 = dblBuf;
                extract_16(col16, depthArr, repeat_timeArr, base1, nflag, call_length);
		//START_TIMER
                	RLEDictDbl->compress(col16, call_length);
		//END_TIMER;
		//PRINT_TIMER_SEC("column 16");
                copyFromData(RLEDictDbl);
                safe_fwrite(param, sizeof(uint64), numParam, file);
                cData->dumpToDisk(file);
                clear();
	}

	//5_7_8_9_14_16
	void read_decompress(int* col5, int * col7, int* col8, int* col9, int* col14, double* col16, const uint64 winSize) {
		
		//first read the N indicator
		int NIndicator = 0;
		safe_fread(&NIndicator, sizeof(int), 1, file);
		if(NIndicator == 1) {
			for(uint64 i = 0; i < winSize; i++) {
				col5[i] = 0;
				col7[i] = 0;
				col8[i] = 0;
				col9[i] = 0;
				col14[i] = 0;
				col16[i] = 255.00;
			}			

			return;
		}

		byte* dataP = NULL;

		//column 5
		safe_fread(param, sizeof(uint64), numParam, file);
		dataP = cData->getBuf(param[__dataSize]);
		safe_fread(dataP, sizeof(byte), param[__dataSize], file);
		copyToData(RLEDictInt, dataP);
		RLEDictInt->decompress(col5, tmpBuf, tmpBufSize);
		clear();	

		//column 7
                safe_fread(param, sizeof(uint64), numParam, file);
                dataP = cData->getBuf(param[__dataSize]);
                safe_fread(dataP, sizeof(byte), param[__dataSize], file);
                copyToData(RLEDictInt, dataP);
                RLEDictInt->decompress(col7, tmpBuf, tmpBufSize);
                clear();

		//column 8
                safe_fread(param, sizeof(uint64), numParam, file);
                dataP = cData->getBuf(param[__dataSize]);
                safe_fread(dataP, sizeof(byte), param[__dataSize], file);
                copyToData(RLEDictInt, dataP);
                RLEDictInt->decompress(col8, tmpBuf, tmpBufSize);
                clear();

		//column 9
                safe_fread(param, sizeof(uint64), numParam, file);
                dataP = cData->getBuf(param[__dataSize]);
                safe_fread(dataP, sizeof(byte), param[__dataSize], file);
                copyToData(RLEDictInt, dataP);
                RLEDictInt->decompress(col9, tmpBuf, tmpBufSize);
                clear();

		//column 14
                safe_fread(param, sizeof(uint64), numParam, file);
                dataP = cData->getBuf(param[__dataSize]);
                safe_fread(dataP, sizeof(byte), param[__dataSize], file);
                copyToData(RLEDictInt, dataP);
                RLEDictInt->decompress(col14, tmpBuf, tmpBufSize);
                clear();

		//column 16
                safe_fread(param, sizeof(uint64), numParam, file);
                dataP = cData->getBuf(param[__dataSize]);
                safe_fread(dataP, sizeof(byte), param[__dataSize], file);
                copyToData(RLEDictDbl, dataP);
                RLEDictDbl->decompress(col16, tmpBuf, tmpBufSize);
                clear();
	}


	template<class DataType>
	void copyToData(RLEDictData<DataType>* outData, byte* inData) {
	
		//parameters
                outData->numRun = param[__numRun];
                outData->numWordRunValue = param[__numWordValue];
                outData->numWordRunLength = param[__numWordLength];
                outData->numEleRunValueDict = param[__numEleValueDict];
                outData->numEleRunLengthDict = param[__numEleLengthDict];

		//assign the memory
		outData->cRunValue = (word*)inData;
		outData->cRunLength = (word*)(inData + sizeof(word)*param[__numWordValue]);
		outData->runValueDict = (DataType*)(inData + sizeof(word)*param[__numWordValue] + sizeof(word)*param[__numWordLength]);
		outData->runLengthDict = (uint32*)(inData + sizeof(word)*param[__numWordValue] + sizeof(word)*param[__numWordLength] + 
						sizeof(DataType)*param[__numEleValueDict]);
	}

	//copy and arrange parameters and data
	template<class DataType>
	void copyFromData(RLEDictData<DataType>* inData) {
		//parameters
		param[__numRun] = inData->numRun;
		param[__numWordValue] = inData->numWordRunValue;
		param[__numWordLength] = inData->numWordRunLength;
		param[__numEleValueDict] = inData->numEleRunValueDict;
		param[__numEleLengthDict] = inData->numEleRunLengthDict;
		param[__dataSize] = inData->getDataSize();

		//data
		byte* buf = cData->getBuf(param[__dataSize]);
		memcpy(buf, inData->cRunValue, sizeof(word)*param[__numWordValue]);
		memcpy(buf + sizeof(word)*param[__numWordValue], 
			inData->cRunLength, sizeof(word)*param[__numWordLength]);	
		memcpy(buf + sizeof(word)*param[__numWordValue] + sizeof(word)*param[__numWordLength], 
			inData->runValueDict, sizeof(DataType)*param[__numEleValueDict]);
		memcpy(buf + sizeof(word)*param[__numWordValue] + sizeof(word)*param[__numWordLength] 
			+ sizeof(DataType)*param[__numEleValueDict], 
			inData->runLengthDict, sizeof(uint32)*param[__numEleLengthDict]);
	}
};


struct Column_10 {
	uint64 maxWinSize;

	FILE* file;
	Buffer* buf;

	static const int numParam = 3;
	uint64 param[3];
	static const int __numN = 0;
	static const int __numWord = 1;
	static const int __numElement = 2;

	Column_10(const uint64 in_maxWinSize) {
		maxWinSize = in_maxWinSize;
		buf = new Buffer(sizeof(char)*maxWinSize + sizeof(int)*maxWinSize);
	};

	void close() {
		buf->close();
	}

	void startWrite(string fileName) {
		safe_fopen(&file, string(fileName + ".10").c_str(), "wb");
	}

	void endWrite() {
		safe_fclose(file);
	}

        void startRead(string fileName) {
                safe_fopen(&file, string(fileName + ".10").c_str(), "rb");
        }

        void endRead() {
                safe_fclose(file);
        }


	void read_decompress(char* col10, const int call_length) {
		assert(call_length <= maxWinSize);

		int NIndicator = 0;
		safe_fread(&NIndicator, sizeof(int), 1, file);
		if(NIndicator == 1) {
			for(int i = 0; i < call_length; i++) {
				col10[i] = 'N';
			}
			return;
		} 

		safe_fread(param, sizeof(uint64), numParam, file);
		assert(param[__numElement] == call_length);
		code_t* code = (code_t*)buf->getBuf(sizeof(code_t)*param[__numWord]);
		uint32* nindex = (uint32*)buf->getBuf(sizeof(uint32)*param[__numN]);
		safe_fread(code, sizeof(code_t), param[__numWord], file);
		safe_fread(nindex, sizeof(uint32), param[__numN], file);

		decode(col10, call_length, code, param[__numWord]);
		
		for(uint64 i = 0; i < param[__numN]; i++) {
			col10[nindex[i]] = 'N';
		}
		
		buf->clear();
	}

	void compress_write(char* col10, const bool* nflag, const int call_length) {
		assert(call_length <= maxWinSize);

		//write a N indicator
		int NIndicator = 0;
		if(isAllN(nflag, call_length)) {
			NIndicator = 1;
			safe_fwrite(&NIndicator, sizeof(int), 1, file);
			return;
		} else {
                        NIndicator = 0;
                        safe_fwrite(&NIndicator, sizeof(int), 1, file);
		}

		//the memory pointer
		const uint64 numWordCode = getNumBinWord(call_length);
		code_t* code = (code_t*)buf->getBuf(sizeof(code_t)*numWordCode);
		uint32* nindex = (uint32*)buf->getBuf(sizeof(uint32)*call_length);

		//count the N
		uint64 numN = 0;
		for(uint64 i = 0; i < call_length; i++) {
			if(col10[i] == 'N') {
				nindex[numN] = i;
				col10[i] = 'A';		//replace 'N' as 'A' for encoding
				numN++;
			}
		}	

		//set all parameters
		param[__numN] = numN;
		param[__numWord] = numWordCode;
		param[__numElement] = call_length;

		//encoding
		encode(code, col10, call_length);
		
		//write
		safe_fwrite(param, sizeof(uint64), numParam, file);
		safe_fwrite(code, sizeof(code_t), numWordCode, file);
		safe_fwrite(nindex, sizeof(uint32), numN, file);
	
		//clear the buffer for the next window
		buf->clear();
	}	
};

inline void extract_10(char* col10, const char* base1, const char* base2, const bool* nflag, const int call_length) {
	for (int j = 0; j != call_length; j++) {
		if (nflag[j]) { 
               		col10[j] = 'N';
			continue;
		} else if (base1[j] < 4 && base2[j] < 4) {
			col10[j] = ("ACTGNNNN"[base2[j]]);
		} else if (base1[j] < 4) {
			col10[j] = 'N';
		} else {
			col10[j] = 'N';
		}

	}
}

inline void extract_11_12_13(int* col11, int* col12, int* col13,
				const int* q_sumArr, const int* count_uni_buf, const int* count_allArr, 
			    const char* base1, const char* base2, const bool* nflag, const int call_length) {
                for (int j = 0; j != call_length; j++) {
                        if (nflag[j]) {
                                col11[j] = 0;
                                col12[j] = 0;
                                col13[j] = 0;
                                continue;
                        } else if (base1[j] < 4 && base2[j] < 4) {
                                        //col11[j] = (sites[j].q_sum[base2[j]] == 0?0:sites[j].q_sum[base2[j]]/sites[j].count_uni[base2[j]]);
                                        //col12[j] = sites[j].count_uni[base2[j]];
                                        //col13[j] = sites[j].count_all[base2[j]];
                                          
				col11[j] = (q_sumArr[j*4 + base2[j]] == 0 ? 0 : q_sumArr[j*4 + base2[j]]/count_uni_buf[j*4 + base2[j]]);
                                col12[j] = count_uni_buf[j * 4 + base2[j]];
                                col13[j] = count_allArr[j * 4 + base2[j]];


                        } else if (base1[j] < 4) {
                                col11[j] = 0;
                                col12[j] = 0;
                                col13[j] = 0;
                        } else {
                                col11[j] = 0;
                                col12[j] = 0;
                                col13[j] = 0;
                        }
                }
}

/**
 * strategy:
 * we first record the differences id of column 11 and 12
 *
 * */
struct Column_11_12_13 {

	FILE* file;

	uint64 maxWinSize;
	Buffer* buf;
	word* wordBuf;
	int* dictBuf;

	//parameters
	const static int numParam = 4;
	const static int __winSize = 0;
	const static int __numCommon = 1;
	const static int __numExtraCol12 = 2;
	const static int __numExtraCol13 = 3; 
	uint64 param[4];

	//memory pointers
	int* commonId;
	int* commonCol11Val;
	int* commonCol12Val;
	int* commonCol13Val;
	int* extraCol12Id;
	int* extraCol12Val;
	int* extraCol13Id;
	int* extraCol13Val;

	Column_11_12_13(const uint64 in_maxWinSize) {
		maxWinSize = in_maxWinSize;

		const uint64 bufSize = sizeof(uint32)*maxWinSize + 	//the common id
					sizeof(int)*maxWinSize*3 +   	//the three columns values with common id
					sizeof(int)*maxWinSize*4;	//the additional column ids and values in column 12 and 13
		buf = new Buffer(bufSize);
		wordBuf = (word*)malloc(sizeof(word)*maxWinSize);
		dictBuf = (int*)malloc(sizeof(int)*maxWinSize);
	};

	void close() {
		buf->close();
		free(wordBuf);
		free(dictBuf);
	}

	void clear() {
		buf->clear();
	}

        void startWrite(string fileName) {
                safe_fopen(&file, string(fileName + ".11_12_13").c_str(), "wb");
        }

        void endWrite() {
                safe_fclose(file);
        }

        void startRead(string fileName) {
                safe_fopen(&file, string(fileName + ".11_12_13").c_str(), "rb");
        }

        void endRead() {
                safe_fclose(file);
        }


	void assignMemory(const int call_length) {
                commonId = (int*)buf->getBuf(sizeof(int)*call_length);
                commonCol11Val = (int*)buf->getBuf(sizeof(int)*call_length);
                commonCol12Val = (int*)buf->getBuf(sizeof(int)*call_length);
                commonCol13Val = (int*)buf->getBuf(sizeof(int)*call_length);
                extraCol12Id = (int*)buf->getBuf(sizeof(int)*call_length);
                extraCol12Val = (int*)buf->getBuf(sizeof(int)*call_length);
                extraCol13Id = (int*)buf->getBuf(sizeof(int)*call_length);
                extraCol13Val = (int*)buf->getBuf(sizeof(int)*call_length);
	}

	void read_decompress(int* col11, int* col12, int* col13, const int winSize) {
		assert(winSize <= maxWinSize);

		int NIndicator = 0;
		safe_fread(&NIndicator, sizeof(int), 1, file);
		if(NIndicator == 1) {
			for(int i = 0; i < winSize; i++) {
				col11[i] = 0;
				col12[i] = 0;
				col13[i] = 0;
			}	
			return;
		}

		assignMemory(winSize);

		//read the uncompressed data
		safe_fread(param, sizeof(uint64), numParam, file);
		assert(winSize == param[__winSize]);
                safe_fread(extraCol12Id, sizeof(int), param[__numExtraCol12], file);
                safe_fread(extraCol12Val, sizeof(int), param[__numExtraCol12], file);
                safe_fread(extraCol13Id, sizeof(int), param[__numExtraCol13], file);
                safe_fread(extraCol13Val, sizeof(int), param[__numExtraCol13], file);
                safe_fread(commonId, sizeof(int), param[__numCommon], file);	

		//read the compressed data

		uint64 dictParam[2];
                safe_fread(dictParam, sizeof(uint64), 2, file);
                safe_fread(wordBuf, sizeof(word), dictParam[0], file);
                safe_fread(dictBuf, sizeof(int), dictParam[1], file);
		dictDecompression(commonCol11Val, param[__numCommon], wordBuf, dictParam[0], dictBuf, dictParam[1]);

                safe_fread(dictParam, sizeof(uint64), 2, file);
                safe_fread(wordBuf, sizeof(word), dictParam[0], file);
                safe_fread(dictBuf, sizeof(int), dictParam[1], file);                
		dictDecompression(commonCol12Val, param[__numCommon], wordBuf, dictParam[0], dictBuf, dictParam[1]);

                safe_fread(dictParam, sizeof(uint64), 2, file);
                safe_fread(wordBuf, sizeof(word), dictParam[0], file);
                safe_fread(dictBuf, sizeof(int), dictParam[1], file);
		dictDecompression(commonCol13Val, param[__numCommon], wordBuf, dictParam[0], dictBuf, dictParam[1]);

	
		//extraction
		memset(col11, 0, sizeof(int)*winSize);
		memset(col12, 0, sizeof(int)*winSize);	
		memset(col13, 0, sizeof(int)*winSize);
		for(int i = 0; i < param[__numCommon]; i++) {
			col11[commonId[i]] = commonCol11Val[i];
			col12[commonId[i]] = commonCol12Val[i];	
			col13[commonId[i]] = commonCol13Val[i];
		}
		for(int i = 0; i < param[__numExtraCol12]; i++) {
			col12[extraCol12Id[i]] = extraCol12Val[i];
		}
                for(int i = 0; i < param[__numExtraCol13]; i++) {
                        col13[extraCol13Id[i]] = extraCol13Val[i];
                }

		clear();		
	}

	void compress_write(int* col11, int* col12, int* col13, const bool* nflag, const uint64 call_length) {
		
		assert(call_length <= maxWinSize);

		int NIndicator = 0;
		if(isAllN(nflag, call_length)) {
			NIndicator = 1;
			safe_fwrite(&NIndicator, sizeof(int), 1, file);

			return;
		} else {
                        NIndicator = 0;
                        safe_fwrite(&NIndicator, sizeof(int), 1, file);
		}

		//assign memory pointer
		assignMemory(call_length);		

		uint64 numCommon = 0;
		uint64 numExtraCol12 = 0;
		uint64 numExtraCol13 = 0;
	
		for(int i = 0; i < call_length; i++) {

			//record common id and values
			if(col11[i] != 0) {
				commonId[numCommon] = i;
				commonCol11Val[numCommon] = col11[i];
				commonCol12Val[numCommon] = col12[i];
				commonCol13Val[numCommon] = col13[i];	

				numCommon++;
			} else if(col11[i] != col12[i]) {
				extraCol12Id[numExtraCol12] = i;
				extraCol12Val[numExtraCol12] = col12[i];
				numExtraCol12++;
			} else if(col11[i] != col13[i]) {
				extraCol13Id[numExtraCol13] = i;
				extraCol13Val[numExtraCol13] = col13[i];
				numExtraCol13++;
			} 
		}

		param[__winSize] = call_length;
		param[__numCommon] = numCommon;
		param[__numExtraCol12] = numExtraCol12;
		param[__numExtraCol13] = numExtraCol13;
 
		//write parameters and extranCol12Id, extraCol12Val, extraCol13Id, extraCol13Val
		safe_fwrite(param, sizeof(uint64), numParam, file);
		safe_fwrite(extraCol12Id, sizeof(int), numExtraCol12, file);	
		safe_fwrite(extraCol12Val, sizeof(int), numExtraCol12, file);
		safe_fwrite(extraCol13Id, sizeof(int), numExtraCol13, file);
		safe_fwrite(extraCol13Val, sizeof(int), numExtraCol13, file);
		safe_fwrite(commonId, sizeof(int), numCommon, file);
	
		//compress the commonCol1Val, commonCol2Val, commonCol3Val using dictionary encoding
		//numWordOutData, numElementDict

		uint64 dictParam[2];

                dictParam[0] = maxWinSize;
                dictParam[1] = maxWinSize;
                dictCompression(wordBuf, &(dictParam[0]), dictBuf, &(dictParam[1]), commonCol11Val, numCommon);
                safe_fwrite(dictParam, sizeof(uint64), 2, file);
                safe_fwrite(wordBuf, sizeof(word), dictParam[0], file);
                safe_fwrite(dictBuf, sizeof(int), dictParam[1], file);

                dictParam[0] = maxWinSize;
                dictParam[1] = maxWinSize;
                dictCompression(wordBuf, &(dictParam[0]), dictBuf, &(dictParam[1]), commonCol12Val, numCommon);
                safe_fwrite(dictParam, sizeof(uint64), 2, file);
                safe_fwrite(wordBuf, sizeof(word), dictParam[0], file);
                safe_fwrite(dictBuf, sizeof(int), dictParam[1], file);


                dictParam[0] = maxWinSize;
                dictParam[1] = maxWinSize;
                dictCompression(wordBuf, &(dictParam[0]), dictBuf, &(dictParam[1]), commonCol13Val, numCommon);
                safe_fwrite(dictParam, sizeof(uint64), 2, file);
                safe_fwrite(wordBuf, sizeof(word), dictParam[0], file);
                safe_fwrite(dictBuf, sizeof(int), dictParam[1], file);
		
		clear();
	}
};


struct Column_15_17 {

	uint64 maxWinSize;
	Buffer* buf;
	FILE* file;

	
	Column_15_17(const uint64 in_maxWinSize) {
		maxWinSize = in_maxWinSize;
		buf = new Buffer(sizeof(int)*maxWinSize + sizeof(double)*maxWinSize);
	}

	void close() {
		buf->close();
	}

	void clear() {
		buf->clear();
	}

        void startWrite(string fileName) {
                safe_fopen(&file, string(fileName + ".15_17").c_str(), "wb");
        }

        void endWrite() {
                safe_fclose(file);
        }

        void startRead(string fileName) {
                safe_fopen(&file, string(fileName + ".15_17").c_str(), "rb");
        }

        void endRead() {
                safe_fclose(file);
        }

	void read_decompress(double* col15, int* col17, const uint64 winSize) {
		assert(winSize <= maxWinSize);

		//if all N
		int NIndicator = 0;
		safe_fread(&NIndicator, sizeof(int), 1, file);
		if(NIndicator == 1) {
			for(uint64 i = 0; i < winSize; i++) {
				col15[i] = 1.0f;
				col17[i] = 0;
			}
			return;
		}

		//set the base value
                for(uint64 i = 0; i < winSize; i++) {
               		col15[i] = 1.0f;
                	col17[i] = 0;
                }

		//read and decompress column 15
		int n = 0;
		safe_fread(&n, sizeof(int), 1, file);
                int* id = (int*)buf->getBuf(sizeof(int)*n);
                double* val = (double*)buf->getBuf(sizeof(double)*n);
		safe_fread(id, sizeof(int), n, file);
		safe_fread(val, sizeof(double), n, file);
		for(uint64 i = 0; i < n; i++) {
			col15[id[i]] = val[i];
		}
		buf->clear();

		//read and decompress column 17
                safe_fread(&n, sizeof(int), 1, file);
                id = (int*)buf->getBuf(sizeof(int)*n);
		safe_fread(id, sizeof(int), n, file);
		for(uint64 i = 0; i < n; i++) {
			col17[id[i]] = 1;
		}
		buf->clear();
	}


	void compress_write(const char* base1, const char* base2, const bool* nflag, 
			    double* rank_sum_test_value, const char* oriArr, const uint64 call_length) {
		assert(call_length <= maxWinSize);

		//N indicator
		int NIndicator = 0;
		if(isAllN(nflag, call_length)) {
			NIndicator = 1;
			safe_fwrite(&NIndicator, sizeof(int), 1, file);
			return;
		} else {
                        NIndicator = 0;
                        safe_fwrite(&NIndicator, sizeof(int), 1, file);
		}

		//compress and write column 15
		int* id = (int*)buf->getBuf(sizeof(int)*call_length);
		double* val = (double*)buf->getBuf(sizeof(double)*call_length);
		int n = 0;

		for (uint64 j = 0; j != call_length; j++) {
			if (nflag[j]) { 
				continue;
			} else if (base1[j] < 4) {
        	                if(rank_sum_test_value[j] != 1.0) {
        	                        id[n] = j;
               	                	val[n] = rank_sum_test_value[j];
                                	n++;
                        	}
			} else {}
		}
		//printf("======> column 15, n = %d\n", n);
		safe_fwrite(&n, sizeof(int), 1, file);
		safe_fwrite(id, sizeof(int), n, file);
		safe_fwrite(val, sizeof(double), n, file);
		buf->clear();

		//compress and write column 17
		id = (int*)buf->getBuf(sizeof(int)*call_length);
		n = 0;
                for (uint64 j = 0; j != call_length; j++) {
                        if (nflag[j]) {        
                                continue;
                        } else if (base1[j] < 4) {
                                if(((oriArr[j] & 8) ? 1 : 0) != 0) {
                                        id[n] = j;
                                        n++;
                                }
                        } else {}
                }
		//printf("======> column 17, n = %d\n", n);
		safe_fwrite(&n, sizeof(int), 1, file);
		safe_fwrite(id, sizeof(int), n, file);
		buf->clear();
	}	
};


#endif /* __COMPRESS_H__ */
