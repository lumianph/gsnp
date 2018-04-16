#ifndef __GSNP_UTIL_H__
#define __GSNP_UTIL_H__

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sys/time.h>
#include <assert.h>
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <safe_driver.h>

using namespace std;

/* data types */
typedef unsigned char byte;
typedef unsigned int word32;
typedef unsigned long word64;
typedef unsigned char uint8;
typedef unsigned short uint16;
typedef unsigned int uint32;
typedef unsigned long uint64;

typedef word32 word;

#define NUM_BYTE_PER_WORD32 (sizeof(word32))
#define NUM_BYTE_PER_WORD64 (sizeof(word64))


/* timing functions */
inline double getSec( struct timeval tvStart, struct timeval tvEnd ) {
        double tStart = (double)tvStart.tv_sec + 1e-6*tvStart.tv_usec;
        double tEnd = (double)tvEnd.tv_sec + 1e-6*tvEnd.tv_usec;
        return (tEnd - tStart);
}

#define INIT_TIMER struct timeval start, end;
#define START_TIMER gettimeofday(&start, NULL);
#define END_TIMER   gettimeofday(&end, NULL);
#define PRINT_TIMER_SEC(msg) printf("*** %s: %.3f sec ***\n", msg, getSec(start, end));
#define GET_TIMER_SEC getSec(start, end);


/* result checking functions */
template<class DataType>
void simpleCheck(DataType* data1, DataType* data2, 
		const uint64 numElement, const char* title = NULL) {
	printf("%s checking result .............................\n", title);
	for(uint64 i = 0; i < numElement; i++) {
		if(data1[i] != data2[i]) {
			cout << "\ti = " << i << ", data1 = " << data1[i] << ", data2 = " << data2[i] << endl;
			printf("\tFAILED!!!!!!!!!\n");
			exit(1);
		}
	}
	
	printf("\tPassed :)\n");
}

template<class DataType>
void simpleCheck(DataType* data1, DataType* data2,
                const uint64 numElement, const double eps, 
		const char* title = NULL) {
        printf("%s checking result .............................\n", title);
        for(uint64 i = 0; i < numElement; i++) {
                if(fabs((data1[i] - data2[i])/(double)data1[i]) > eps) {
                        cout << "\ti = " << i << ", data1 = " << data1[i] << ", data2 = " << data2[i] << endl;
                        printf("\tFAILED!!!!!!!!!\n");
                        exit(1);
                }
        }

        printf("\tPassed :)\n");
}

/* print helper functions */
#define PRINT_FUNC_NAME printf("*** FUNCTION: %s ***\n", __func__);

inline void printLine() {
	printf("---------------------------------------\n");
}

inline void printTitle(const char* msg) {
	printf("------------- %s ---------------\n", msg);
}

inline void printString(const char* str, const uint64 numChar) {
	for(uint i = 0; i < numChar; i++) {
		printf("%c", str[i]);
	}
	printf("\n");
}


/* error handling functions */
#define ERROR_EXIT(msg) printf("!!!ERROR@%s: %s\n", __func__, msg); exit(EXIT_FAILURE);


/* random number generators */
inline byte randByte() {
	return rand()%256;
}


inline uint32 rand_uint32() {
	return (((uint32)randByte())<<24)|(((uint32)randByte())<<16)|(((uint32)randByte())<<8)|(randByte());
}


/* bit set and get functions, for a byte array, position is from the left to right, 
 * the position index is starting from 0
 * */
//look-up table
const static byte bit1Table[] = {0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01};
const static uint8 bit1CountTable[256] = {
0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
};

//!!!NOTE: these macros are all based on byte-base (8-bit) memory representation
#define BIT_SET1(v, i) ((v)[(i)/8] |= bit1Table[(i)%8])                 
#define BIT_GET(v, i) (((v)[(i)/8])&bit1Table[(i)%8])                  //return 0 if the bit is 0, othereise return non-0
#define BIT_GET_FLAG(v, i) (BIT_GET(v, i) == 0) ? (0) : (1)
//return the number of bit "1" in the vector
inline uint64 countBit1(const byte* buf, const uint64 numByte) {
	uint64 count1 = 0;
	for(uint64 i = 0; i < numByte; i++) {
		count1 += bit1CountTable[buf[i]];
	}

	return count1;
}


/*comparison functions */
template<class DataType>
DataType max(const DataType* data, const uint64 numElement) {
	DataType maxValue = data[0];

	for(uint64 i = 1; i < numElement; i++) {
		if(data[i] > maxValue) {
			maxValue = data[i];
		}
	}
	
	return maxValue;
}

#endif /* __GSNP_UTIL_H__ */

