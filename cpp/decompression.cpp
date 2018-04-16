

#include "compress.h"
#include "common.h"
#include <gsnp_util.h>



struct CNSDataWindow {
	static const int maxChrName = 128;

	uint64 maxNumRows;
	uint64 numRows;

	static const int numColumns = 17;

	//17 columns
	char* chrIdBuf;
	char** chrId;
	int* tPos;
	char* tRef;
	char* tCns;
	int* tQual;
	char* tBase1;
	int* tAvgQ1;
	int* tCountUni1;
	int* tCountAll1;
	char* tBase2;
	int* tAvgQ2;
	int* tCountUni2;
	int* tCountAll2;
	int* tDepth;
	double* rankSum;
	double* copySum;
	int* SNPStatus;

	void setNumRows(const uint64 in_numRows) {
		if(in_numRows > maxNumRows) {
			ERROR_EXIT("in_numRows > maxNumRows");
		} else {
			numRows = in_numRows;
		}
	}	

	CNSDataWindow(const uint64 in_maxNumRows) {
		maxNumRows = in_maxNumRows;
		numRows = 0;

	        chrId = (char**)malloc(sizeof(char*)*maxNumRows);//1
		chrIdBuf = (char*)malloc(sizeof(char)*maxChrName*maxNumRows);
		for(uint64 i = 0; i < maxNumRows; i++) {
			chrId[i] = chrIdBuf + maxChrName*i;
		}
	        tPos = (int*)malloc(sizeof(int)*maxNumRows);//2
	        tRef = (char*)malloc(sizeof(char)*maxNumRows);//3
	        tCns = (char*)malloc(sizeof(char)*maxNumRows);//4
	        tQual = (int*)malloc(sizeof(int)*maxNumRows);//5
	        tBase1 = (char*)malloc(sizeof(char)*maxNumRows);//6
	        tAvgQ1 = (int*)malloc(sizeof(int)*maxNumRows);//7
	        tCountUni1 = (int*)malloc(sizeof(int)*maxNumRows);//8
	        tCountAll1 = (int*)malloc(sizeof(int)*maxNumRows);//9
	        tBase2 = (char*)malloc(sizeof(char)*maxNumRows);//10
	        tAvgQ2 = (int*)malloc(sizeof(int)*maxNumRows);//11
	        tCountUni2 = (int*)malloc(sizeof(int)*maxNumRows);//12
	        tCountAll2 = (int*)malloc(sizeof(int)*maxNumRows);//13
	        tDepth = (int*)malloc(sizeof(int)*maxNumRows);//14
       		rankSum = (double*)malloc(sizeof(double)*maxNumRows);//15
        	copySum = (double*)malloc(sizeof(double)*maxNumRows);//16
        	SNPStatus = (int*)malloc(sizeof(int)*maxNumRows);//17
	}


        void print() {

		printf("numRows = %lu\n", numRows);
		printLine();

		for(uint64 i = 0; i < numRows; i++) {
                	printf("%s\t%d\t%c\t%c\t%d\t%c\t%d\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%.5f\t%.3f\t%d\t\n",
                                chrId[i], (tPos[i]), (tRef[i]), (tCns[i]), (tQual[i]),
                                (tBase1[i]), (tAvgQ1[i]), (tCountUni1[i]), (tCountAll1[i]),
                                (tBase2[i]), (tAvgQ2[i]), (tCountUni2[i]), (tCountAll2[i]),
                                (tDepth[i]), (rankSum[i]), (copySum[i]), (SNPStatus[i]));

		}
        }

	void close() {
		//TODO
	}
};



struct CNSFile {
	FILE* file;

	CNSFile(const char* fileName) {
		safe_fopen(&file, fileName, "r");
	}


	//coninuous read a window of data
	void read(CNSDataWindow* data, const uint64 numRows) {
		data->setNumRows(numRows);

		int numItemRead = 0;
		for(uint64 i = 0; i < numRows; i++) {
			numItemRead = fscanf(file, "%s\t%d\t%c\t%c\t%d\t%c\t%d\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%lf\t%lf\t%d\n", 
				data->chrId[i], &(data->tPos[i]), &(data->tRef[i]), &(data->tCns[i]), &(data->tQual[i]), 
				&(data->tBase1[i]), &(data->tAvgQ1[i]), &(data->tCountUni1[i]), &(data->tCountAll1[i]),
				&(data->tBase2[i]), &(data->tAvgQ2[i]), &(data->tCountUni2[i]), &(data->tCountAll2[i]), 
				&(data->tDepth[i]), &(data->rankSum[i]), &(data->copySum[i]), &(data->SNPStatus[i]));

			if(numItemRead != data->numColumns) {
				printf("i = %lu, numItemRead = %d, numColumns = %d\n", 
					i, numItemRead, data->numColumns);

				ERROR_EXIT("numItemRead != data->numColumns");
			}
		}
	}


	void close() {
		safe_fclose(file);
	}
};

int decompressAll(int argc, char** argv) {

	if(argc != 3) {
		printf("Two parameters:\n");
		printf("1. the input compressed file.\n");
		printf("2. the output file name.\n");
		return EXIT_FAILURE;		
	}

	const char* cnsFileName = argv[1]; 
	const char* outFileName = argv[2];
	ofstream outFile(outFileName);

	
	/*read the meta information, column 1 and column 2*/
	Column_1_2* c_1_2 = new Column_1_2();
	c_1_2->read(cnsFileName);
	const uint64 maxWinSize = max(c_1_2->winSize, c_1_2->numWindow);
	
	/*prepare the memory*/
	char* ref_col3 = (char*)malloc(sizeof(char)*maxWinSize);
	char* ref_col4 = (char*)malloc(sizeof(char)*maxWinSize);
	char* ref_col6 = (char*)malloc(sizeof(char)*maxWinSize);
	int* ref_col5 = (int*)malloc(sizeof(int)*maxWinSize);
	int* ref_col7 = (int*)malloc(sizeof(int)*maxWinSize);
	int* ref_col8 = (int*)malloc(sizeof(int)*maxWinSize);
	int* ref_col9 = (int*)malloc(sizeof(int)*maxWinSize);
	int* ref_col14 = (int*)malloc(sizeof(int)*maxWinSize);
	double* ref_col16 = (double*)malloc(sizeof(double)*maxWinSize);
	char* ref_col10 = (char*)malloc(sizeof(char)*maxWinSize);
	int* ref_col11 = (int*)malloc(sizeof(int)*maxWinSize);
	int* ref_col12 = (int*)malloc(sizeof(int)*maxWinSize);
	int* ref_col13 = (int*)malloc(sizeof(int)*maxWinSize);
	double* ref_col15 = (double*)malloc(sizeof(double)*maxWinSize);
	int* ref_col17 = (int*)malloc(sizeof(int)*maxWinSize);
	Column_3_4_6* c_3_4_6 = new Column_3_4_6(maxWinSize);
	c_3_4_6->startRead(cnsFileName);
	Column_RLEDict* c_5_7_8_9_14_16 = new Column_RLEDict(maxWinSize);	
	c_5_7_8_9_14_16->startRead(cnsFileName);
	Column_10* c_10 = new Column_10(maxWinSize);
	c_10->startRead(cnsFileName);
	Column_11_12_13* c_11_12_13 = new Column_11_12_13(maxWinSize);
	c_11_12_13->startRead(cnsFileName);
	Column_15_17* c_15_17 = new Column_15_17(maxWinSize);
	c_15_17->startRead(cnsFileName);

	/*decompress window by window*/
	for(uint64 winId = 0; winId < c_1_2->numWindow; winId++) {

		/*column 1 and 2*/
		//already read
		
		/*column 3, 4, 6*/
		c_3_4_6->read();
		c_3_4_6->decompress(ref_col3, ref_col4, ref_col6);


		/*column 5, 7, 8, 9, 14, 16*/
		c_5_7_8_9_14_16->read_decompress(ref_col5, ref_col7, ref_col8, 
				ref_col9, ref_col14, ref_col16, c_1_2->winSize[winId]);

		/*column 10*/
		c_10->read_decompress(ref_col10, c_1_2->winSize[winId]);

		/*column 11, 12, 13*/
		c_11_12_13->read_decompress(ref_col11, ref_col12, ref_col13, c_1_2->winSize[winId]);
 
		/*column 15, 17*/
		c_15_17->read_decompress(ref_col15, ref_col17, c_1_2->winSize[winId]);
		
	
		/*now we have all data in the window, output to a file*/	
		for(uint64 i = 0; i < c_1_2->winSize[winId]; i++) {


			if(ref_col3[i] == 'N' && ref_col4[i] == 'N' && ref_col6[i] == 'N') {
				outFile << c_1_2->callName[winId] << '\t' <<
					c_1_2->posStart[winId] + i << '\t' << 
					"N\tN\t0\tN\t0\t0\t0\tN\t0\t0\t0\t0\t1.000\t255.000\t0" << endl;
			} else {
				outFile << c_1_2->callName[winId] << '\t' << c_1_2->posStart[winId] + i
                                                        << '\t' << ref_col3[i] << '\t' << ref_col4[i] << '\t' 
							<< ref_col5[i] << '\t' << ref_col6[i] << '\t'
                                                        << ref_col7[i] << '\t' << ref_col8[i] << '\t'
                                                        << ref_col9[i] << '\t' << ref_col10[i] << '\t'
                                                        << ref_col11[i] << '\t' << ref_col12[i] << '\t'
                                                        << ref_col13[i] << '\t' << ref_col14[i] << '\t' << showpoint
                                                        << ref_col15[i] << '\t' << ref_col16[i] << '\t'
                                                        << ref_col17[i] << endl;

			}	
		}
	}

	
	/*cleanup*/
	//safe_fclose(outFile);
	outFile.close();
	c_1_2->close();
	c_3_4_6->endRead();
	c_3_4_6->close();
	c_5_7_8_9_14_16->close();	
	c_10->endRead();
	c_10->close();
	c_11_12_13->endRead();
	c_11_12_13->close();
	c_15_17->endRead();
	c_15_17->close();
	free(ref_col3);
	free(ref_col4);
	free(ref_col5);
	free(ref_col6);
	free(ref_col7);
	free(ref_col8);
	free(ref_col9);
	free(ref_col10);
	free(ref_col11);
	free(ref_col12);
	free(ref_col13);
	free(ref_col14);
	free(ref_col15);
	free(ref_col16);
	free(ref_col17);
	
	return EXIT_SUCCESS;
}

int check(int argc, char** argv) {

	if(argc == 1) {
		printf("1. the gold cns file name.\n");
		printf("2. the reference cns file name.\n");
		return EXIT_FAILURE;
	}

	//const uint64 winSize = global_win_size;
	const char* cnsFileName = argv[2];
	const char* goldFileName = argv[1];
	INIT_TIMER;


	double woCompSec = 0.0;
	double wCompSec = 0.0;


	/*read the meta information*/
	Column_1_2* c_1_2 = new Column_1_2();
	c_1_2->read(cnsFileName);
	const uint64 maxWinSize = max(c_1_2->winSize, c_1_2->numWindow);


	printf("gold cns file: %s; cns file: %s\n", goldFileName, cnsFileName);
	printf("maxWinSize = %d, numWindow = %d\n", maxWinSize, c_1_2->numWindow);

	/*prepare the memory*/
	CNSFile* cnsFile = new CNSFile(goldFileName);
	CNSDataWindow* cnsData = new CNSDataWindow(maxWinSize);

	char* ref_col3 = (char*)malloc(sizeof(char)*maxWinSize);
	char* ref_col4 = (char*)malloc(sizeof(char)*maxWinSize);
	char* ref_col6 = (char*)malloc(sizeof(char)*maxWinSize);
	int* ref_col5 = (int*)malloc(sizeof(int)*maxWinSize);
	int* ref_col7 = (int*)malloc(sizeof(int)*maxWinSize);
	int* ref_col8 = (int*)malloc(sizeof(int)*maxWinSize);
	int* ref_col9 = (int*)malloc(sizeof(int)*maxWinSize);
	int* ref_col14 = (int*)malloc(sizeof(int)*maxWinSize);
	double* ref_col16 = (double*)malloc(sizeof(double)*maxWinSize);
	char* ref_col10 = (char*)malloc(sizeof(char)*maxWinSize);
	int* ref_col11 = (int*)malloc(sizeof(int)*maxWinSize);
	int* ref_col12 = (int*)malloc(sizeof(int)*maxWinSize);
	int* ref_col13 = (int*)malloc(sizeof(int)*maxWinSize);
	double* ref_col15 = (double*)malloc(sizeof(double)*maxWinSize);
	int* ref_col17 = (int*)malloc(sizeof(int)*maxWinSize);

	Column_3_4_6* c_3_4_6 = new Column_3_4_6(maxWinSize);
	c_3_4_6->startRead(cnsFileName);
	Column_RLEDict* c_5_7_8_9_14_16 = new Column_RLEDict(maxWinSize);	
	c_5_7_8_9_14_16->startRead(cnsFileName);
	Column_10* c_10 = new Column_10(maxWinSize);
	c_10->startRead(cnsFileName);
	Column_11_12_13* c_11_12_13 = new Column_11_12_13(maxWinSize);
	c_11_12_13->startRead(cnsFileName);
	Column_15_17* c_15_17 = new Column_15_17(maxWinSize);
	c_15_17->startRead(cnsFileName);


	for(uint64 winId = 0; winId < c_1_2->numWindow; winId++) {

		printf("winId = %d\n", winId);

		/*read the gold data*/
		START_TIMER;
		cnsFile->read(cnsData, c_1_2->winSize[winId]);
		END_TIMER;
		PRINT_TIMER_SEC("gold data read");
		woCompSec += GET_TIMER_SEC;

		/*column 1 and 2*/
		for(uint64 i = 0; i < c_1_2->winSize[winId]; i++) {

			assert(strcmp(cnsData->chrId[i], c_1_2->callName[winId]) == 0);
			assert(cnsData->tPos[i] == c_1_2->posStart[winId] + i);
		}
	

		/*column 3, 4, 6*/
		START_TIMER;
		c_3_4_6->read();
		c_3_4_6->decompress(ref_col3, ref_col4, ref_col6);
		END_TIMER;
		wCompSec += GET_TIMER_SEC;
		simpleCheck(cnsData->tRef, ref_col3, c_1_2->winSize[winId], "column 3");
		simpleCheck(cnsData->tCns, ref_col4, c_1_2->winSize[winId], "column 4");
		simpleCheck(cnsData->tBase1, ref_col6, c_1_2->winSize[winId], "column 6");

		/*column 5, 7, 8, 9, 14, 16*/
		START_TIMER
		c_5_7_8_9_14_16->read_decompress(ref_col5, ref_col7, ref_col8, 
				ref_col9, ref_col14, ref_col16, c_1_2->winSize[winId]);
		END_TIMER;
		wCompSec += GET_TIMER_SEC;
		simpleCheck(cnsData->tQual, ref_col5, c_1_2->winSize[winId], "column 5");
		simpleCheck(cnsData->tAvgQ1, ref_col7, c_1_2->winSize[winId], "column 7");
		simpleCheck(cnsData->tCountUni1, ref_col8, c_1_2->winSize[winId], "column 8");
		simpleCheck(cnsData->tCountAll1, ref_col9, c_1_2->winSize[winId], "column 9");
		simpleCheck(cnsData->tDepth, ref_col14, c_1_2->winSize[winId], "column 14");
		simpleCheck(cnsData->copySum, ref_col16, c_1_2->winSize[winId], 0.00001, "column 16");

		/*column 10*/
		START_TIMER;
		c_10->read_decompress(ref_col10, c_1_2->winSize[winId]);
		END_TIMER;
		wCompSec += GET_TIMER_SEC;
		simpleCheck(cnsData->tBase2, ref_col10, c_1_2->winSize[winId], "column 10");

		/*column 11, 12, 13*/
		START_TIMER;
		c_11_12_13->read_decompress(ref_col11, ref_col12, ref_col13, c_1_2->winSize[winId]);
		END_TIMER;
		wCompSec += GET_TIMER_SEC;
		simpleCheck(cnsData->tAvgQ2, ref_col11, c_1_2->winSize[winId], "column 11");
		simpleCheck(cnsData->tCountUni2, ref_col12, c_1_2->winSize[winId], "column 12");
		simpleCheck(cnsData->tCountAll2, ref_col13, c_1_2->winSize[winId], "column 13");
 
		/*column 15, 17*/
		START_TIMER;
		c_15_17->read_decompress(ref_col15, ref_col17, c_1_2->winSize[winId]);
		END_TIMER;
		wCompSec += GET_TIMER_SEC;
		simpleCheck(cnsData->rankSum, ref_col15, c_1_2->winSize[winId], 0.00001, "column 15");
		simpleCheck(cnsData->SNPStatus, ref_col17, c_1_2->winSize[winId], "column 17");
		printLine();

	}

	
	c_1_2->close();
	c_3_4_6->endRead();
	c_3_4_6->close();
	c_5_7_8_9_14_16->close();	
	c_10->endRead();
	c_10->close();
	c_11_12_13->endRead();
	c_11_12_13->close();
	c_15_17->endRead();
	c_15_17->close();

	printf("=====> all done.\n");

	printf("w/o compression, read time: %.3f sec\n", woCompSec);
	printf("w/ compression, read time: %.3f sec\n", wCompSec);

	return EXIT_SUCCESS;
}


int main(int argc, char** argv) {

	return decompressAll(argc, argv);
	//return check(argc, argv);
}

