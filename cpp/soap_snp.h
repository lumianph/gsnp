#ifndef SOAP_SNP_HH_
#define SOAP_SNP_HH_
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <cmath>
#include <iomanip>
#include <cassert>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <gsnp_util.h>
#include "common.h"
#include "dict_compression.h"
//#include "compress.h"

/*
typedef unsigned long long ubit64_t;
typedef unsigned int ubit32_t;
typedef double rate_t;
typedef unsigned char small_int;
using namespace std;
const size_t capacity = sizeof(ubit64_t)*8/4;
const char abbv[17]={'A','M','W','R','M','C','Y','S','W','Y','T','K','R','S','K','G','N'};
const ubit64_t glf_base_code[8]={1,2,8,4,15,15,15,15}; // A C T G
const ubit64_t glf_type_code[10]={0,5,15,10,1,3,2,7,6,11};// AA,CC,GG,TT,AC,AG,AT,CG,CT,GT
const int global_win_size = 4000;
*/


void setCompression();

// Some global variables
//string global_chr_name;
class Files {
public:
	ifstream soap_result, ref_seq, dbsnp, region;
	ofstream consensus, baseinfo, o_region;
	fstream matrix_file;
	Files(){
		soap_result.close();
		ref_seq.close();
		dbsnp.close();
		consensus.close();
		baseinfo.close();
		matrix_file.close();
		region.close();
		o_region.close();
	};
};



#define MAX_NUM_LINE_BUF (1000000) 
struct LineBuffer {
	int bufSize;
	int maxNumLine;
	int numLineBuf;
	string lineBuf[MAX_NUM_LINE_BUF];
	std::ifstream* ifs;
	int numRead;
	double read_sec;
	
	LineBuffer(std::ifstream & in_ifs) {
		bufSize = MAX_NUM_LINE_BUF;
		numLineBuf = 0;
		maxNumLine = 0;
		ifs = &in_ifs;

		numRead = 0;
		read_sec = 0.0f;
	}

	bool getLine(string** line) {
		if(numLineBuf == 0) {	//if the buffer is empyt, then read from the disk
			read();
		}

		if(numLineBuf == 0) {	//if it is still empty, it indicates that there is no data in the file
			return false;
		} else {
			*line = &(lineBuf[maxNumLine - numLineBuf]);
			numLineBuf--;

			return true;
		}
	}

	void read() {

		printf("==============> read the data from disk.\n");

		assert(numLineBuf == 0);	//only read when the buffer is empty
		maxNumLine = 0;

		INIT_TIMER;
		START_TIMER;
		
		while(getline(*ifs, lineBuf[numLineBuf])) {
			numLineBuf++;
			maxNumLine++;

			if(numLineBuf == bufSize) {
				break;
			}
		}

		numRead++;
		END_TIMER;
		read_sec += GET_TIMER_SEC;
	}
};



class Soap_format {
	// Soap alignment result
public:
	std::string read_id, read, qual, chr_name;
	int hit, read_len, position, mismatch;
	char ab, strand;
	double read_sec;

//public:
	Soap_format(){
		read_sec = 0.0f;
	};

	friend std::istringstream & operator>>(std::istringstream & alignment, Soap_format & soap) {
		//INIT_TIMER;
		//START_TIMER;

		alignment>>soap.read_id>>soap.read>>soap.qual>>soap.hit>>soap.ab>>soap.read_len>>soap.strand>>soap.chr_name>>soap.position>>soap.mismatch;
		//cerr<<soap<<endl;
		//exit(1);
		if(soap.mismatch>200) {
			int indel_pos,indel_len;
			string temp("");
			alignment>>indel_pos;
			indel_len = soap.mismatch-200;
			for(int i=0; i!=indel_len; i++) {
				temp = temp+'N';
			}
			soap.read = soap.read.substr(0,indel_pos)+temp+soap.read.substr(indel_pos,soap.read_len-indel_pos);
			soap.qual = soap.qual.substr(0,indel_pos)+temp+soap.qual.substr(indel_pos,soap.read_len-indel_pos);
			//cerr<<soap<<endl;
		}
		else if (soap.mismatch>100) {
			int indel_pos,indel_len;
			alignment>>indel_pos;
			indel_len = soap.mismatch-100;
			soap.read = soap.read.substr(0,indel_pos) + soap.read.substr(indel_pos+indel_len, soap.read_len-indel_pos-indel_len);
			soap.qual = soap.qual.substr(0,indel_pos) + soap.qual.substr(indel_pos+indel_len, soap.read_len-indel_pos-indel_len);
			//cerr<<soap<<endl;
		}
		//cerr<<soap.position<<'\t';
		soap.position -= 1;
		//cerr<<soap.position<<endl;
		//
		//

		//END_TIMER;
		//soap.read_sec += GET_TIMER_SEC;

		return alignment;
	}
	friend std::ostream & operator<<(std::ostream & o, Soap_format & soap) {
		o<<soap.read_id<<'\t'<<soap.read<<'\t'<<soap.qual<<'\t'<<soap.hit<<'\t'<<soap.ab<<'\t'<<soap.read_len<<'\t'<<soap.strand<<'\t'<<soap.chr_name<<'\t'<<soap.position<<'\t'<<soap.mismatch;
		return o;
	}
	char get_base(std::string::size_type coord) {
		return read[coord];
	}
	char get_qual(std::string::size_type coord) {
		return qual[coord];
	}
	bool is_fwd(){
		return (strand=='+');
	}
	int get_read_len(){
		return read_len;
	}
	inline int get_pos(){
		return position;
	}
	std::string get_chr_name(){
		return chr_name;
	}
	int get_hit(){
		return hit;
	}
	bool is_unique(){
		return (hit==1);
	}
	bool is_N(int coord) {
		return (read[coord] == 'N');
	}
};



enum DataMode {READ, WRITE};

struct SoapData {

	DataMode mode;
	static const int maxNumReadBuf = 1000000;
	double write_time;
	double comp_time;
	double compwrite_time;
	double total_time;
	struct timeval totalStart, totalEnd;


	int numReadInBuf;
	int readId;
	int readLength;

	char* read;
	char* qual;
	int* hit;
	int* read_len;
	int* position;
	int* mismatch;
	//char* ab;
	char* strand;
	string chr_name;
	code_t* code;

	FILE* readFile;
	FILE* qualFile;
	FILE* hitFile;
	FILE* readlenFile;
	FILE* positionFile;
	FILE* mismatchFile;
	FILE* abFile;
	FILE* strandFile;	
	FILE* metaFile;

	word* cData;
	int* intDict;
	char* charDict;


	int n_write, n_read;
	
	//initilization to compress and write or read
        SoapData(const int in_readLength, DataMode in_mode) {
		numReadInBuf = 0;	
		readId = 0;
		readLength = in_readLength;
		mode = in_mode;
		n_write = 0;
		n_read = 0;
		write_time = 0.0f;
		comp_time = 0.0f;
		compwrite_time = 0.0f;

		read = (char*)malloc(sizeof(char)*readLength*maxNumReadBuf);
		qual = (char*)malloc(sizeof(char)*readLength*maxNumReadBuf);
		hit = (int*)malloc(sizeof(int)*maxNumReadBuf);
		read_len = (int*)malloc(sizeof(int)*maxNumReadBuf);
		position = (int*)malloc(sizeof(int)*maxNumReadBuf);
		mismatch = (int*)malloc(sizeof(int)*maxNumReadBuf);
		//ab = (char*)malloc(sizeof(char)*maxNumReadBuf);
		strand = (char*)malloc(sizeof(char)*maxNumReadBuf);

                const uint64 numWordCode = getNumBinWord(readLength*maxNumReadBuf);
                code = (code_t*)malloc(sizeof(code_t)*numWordCode);

		if(mode == WRITE) {
			safe_fopen(&readFile, READ_FILE, "wb");
			safe_fopen(&qualFile, QUAL_FILE, "wb");
			safe_fopen(&hitFile, HIT_FILE, "wb");
			safe_fopen(&readlenFile, READLEN_FILE, "wb");
			safe_fopen(&positionFile, POSITION_FILE, "wb");
			safe_fopen(&mismatchFile, MISMATCH_FILE, "wb");
			//safe_fopen(&abFile, AB_FILE, "wb");
			safe_fopen(&strandFile, STRAND_FILE, "wb");
			safe_fopen(&metaFile, "in.meta", "wb");
		} else if(mode == READ) {
                        safe_fopen(&readFile, READ_FILE, "rb");
                        safe_fopen(&qualFile, QUAL_FILE, "rb");
                        safe_fopen(&hitFile, HIT_FILE, "rb");
                        safe_fopen(&readlenFile, READLEN_FILE, "rb");
                        safe_fopen(&positionFile, POSITION_FILE, "rb");
                        safe_fopen(&mismatchFile, MISMATCH_FILE, "rb");
                        //safe_fopen(&abFile, AB_FILE, "wb");
                        safe_fopen(&strandFile, STRAND_FILE, "rb");
			safe_fopen(&metaFile, "in.meta", "rb");

			safe_fread(&n_write, sizeof(int), 1, metaFile);
		} else {
			ERROR_EXIT("Never here!!!!");
		}

		//memset(read, 0, sizeof(char)*readLength*maxNumReadBuf);
		for(uint64 i = 0; i < readLength*maxNumReadBuf; i++) {
			read[i] = 'A';
		}
		memset(qual, 0, sizeof(char)*readLength*maxNumReadBuf);
		memset(hit, 0, sizeof(int)*maxNumReadBuf);
		memset(read_len, 0, sizeof(int)*maxNumReadBuf);
		memset(position, 0, sizeof(int)*maxNumReadBuf);
		memset(mismatch, 0, sizeof(int)*maxNumReadBuf);
		memset(strand, 0, sizeof(char)*maxNumReadBuf);


                cData = (word*)malloc(sizeof(word)*readLength*maxNumReadBuf);
                intDict = (int*)malloc(sizeof(int)*maxNumReadBuf);
                charDict = (char*)malloc(sizeof(char)*readLength*maxNumReadBuf);

	
		if(mode == WRITE) {	
			dictComprStart(readLength*maxNumReadBuf, sizeof(int));	
        	}
		
	}

	

	void close() {


		if(mode == WRITE) {
			compress_write();
			safe_fwrite(&n_write, sizeof(int), 1, metaFile);
			//printf("writing n_write = %d done.\n", n_write);
			//printf("comp_time = %.3f sec\n", comp_time);
			//printf("compwrite_time = %.3f sec\n", compwrite_time);
		}


        	free(read);
	        free(qual);
	        free(hit);
	        free(read_len);
	        free(position);
	        free(mismatch);
	        //free(ab);
	        free(strand);
		free(code);


	        safe_fclose(readFile);
	        safe_fclose(qualFile);
	        safe_fclose(hitFile);
	      	safe_fclose(readlenFile);
       		safe_fclose(positionFile);
	        safe_fclose(mismatchFile);
	        //safe_fclose(abFile);
	        safe_fclose(strandFile);	
		safe_fclose(metaFile);

                free(cData);
                free(intDict);
                free(charDict);

		if(mode == WRITE) {	
			dictComprEnd();
		}
	}




	void insert(Soap_format* soap) {
		//insert a read
		memcpy(read + numReadInBuf*readLength, soap->read.c_str(), sizeof(char)*soap->read_len);
		memcpy(qual + numReadInBuf*readLength, soap->qual.c_str(), sizeof(char)*soap->read_len);
		hit[numReadInBuf] = soap->hit;
		read_len[numReadInBuf] = soap->read_len;
		position[numReadInBuf] = soap->position;
		mismatch[numReadInBuf] = soap->mismatch;
		//ab[numReadInBuf] = soap->ab;
		strand[numReadInBuf] = soap->strand;
		chr_name = soap->chr_name;
		numReadInBuf++;

		//if the buffer is full, then flush it to the disk
		if(numReadInBuf == maxNumReadBuf) {
			//write();	
			compress_write();
			
			assert(numReadInBuf == 0);
		}
	}
 
	//copy a record for updateData
	bool fetch(char* out_read, char* out_qual, int* out_hit, int* out_read_len, 
		   int* out_position, int* out_mismatch, char* out_strand) {

		if(numReadInBuf == 0) {
			decompress_read();
			if(numReadInBuf == 0) {
				return false;
			}
		}

		memcpy(out_read, read + readId*readLength, sizeof(char)*readLength);
		memcpy(out_qual, qual + readId*readLength, sizeof(char)*readLength);
		*out_hit = hit[readId];
		*out_read_len = read_len[readId];
		*out_position = position[readId];
		*out_mismatch = mismatch[readId];
		*out_strand = strand[readId];

		readId++;
		numReadInBuf--;

		return true;
	}

	void write() {

		safe_fwrite(read, sizeof(char), readLength*numReadInBuf, readFile);
		safe_fwrite(qual, sizeof(char), readLength*numReadInBuf, qualFile);
		safe_fwrite(hit, sizeof(int), numReadInBuf, hitFile);
		safe_fwrite(read_len, sizeof(int), numReadInBuf, readlenFile);
		safe_fwrite(position, sizeof(int), numReadInBuf, positionFile);
		safe_fwrite(mismatch, sizeof(int), numReadInBuf, mismatchFile);
		//safe_fwrite(ab, sizeof(char), numReadInBuf, abFile);
		safe_fwrite(strand, sizeof(char), numReadInBuf, strandFile);

		numReadInBuf = 0;
	}

	
	//decompress a block of data
	void decompress_read() {


                uint64 numWordCData = 0;
                uint64 numElementDict;
                uint64 numElementDData = 0;

		if(n_read == n_write) {

			numReadInBuf = 0;
			readId = 0;
			return;
		}

		
		/*
                safe_fread(&numElementDData, sizeof(uint64), 1, readFile);
                safe_fread(&numWordCData, sizeof(uint64), 1, readFile);
                safe_fread(&numElementDict, sizeof(uint64), 1, readFile);	
                safe_fread(cData, sizeof(word), numWordCData, readFile);
                safe_fread(charDict, sizeof(char), numElementDict, readFile);
		dictDecompression(read, numElementDData, cData, numWordCData, charDict, numElementDict);
		*/
		decompress_base();

		/*
                safe_fread(&numElementDData, sizeof(uint64), 1, qualFile);
                safe_fread(&numWordCData, sizeof(uint64), 1, qualFile);
                safe_fread(&numElementDict, sizeof(uint64), 1, qualFile);
                safe_fread(cData, sizeof(word), numWordCData, qualFile);
                safe_fread(charDict, sizeof(char), numElementDict, qualFile);
                dictDecompression(qual, numElementDData, cData, numWordCData, charDict, numElementDict);
		*/
		decompress_qual();

                safe_fread(&numElementDData, sizeof(uint64), 1, hitFile);
                safe_fread(&numWordCData, sizeof(uint64), 1, hitFile);
                safe_fread(&numElementDict, sizeof(uint64), 1, hitFile);
                safe_fread(cData, sizeof(word), numWordCData, hitFile);
                safe_fread(intDict, sizeof(int), numElementDict, hitFile);
                dictDecompression(hit, numElementDData, cData, numWordCData, intDict, numElementDict);

                safe_fread(&numElementDData, sizeof(uint64), 1, readlenFile);
                safe_fread(&numWordCData, sizeof(uint64), 1, readlenFile);
                safe_fread(&numElementDict, sizeof(uint64), 1, readlenFile);
                safe_fread(cData, sizeof(word), numWordCData, readlenFile);
                safe_fread(intDict, sizeof(int), numElementDict, readlenFile);
                dictDecompression(read_len, numElementDData, cData, numWordCData, intDict, numElementDict);

		
		/*
                safe_fread(&numElementDData, sizeof(uint64), 1, positionFile);
                safe_fread(&numWordCData, sizeof(uint64), 1, positionFile);
                safe_fread(&numElementDict, sizeof(uint64), 1, positionFile);
                safe_fread(cData, sizeof(word), numWordCData, positionFile);
                safe_fread(intDict, sizeof(int), numElementDict, positionFile);
                dictDecompression(position, numElementDData, cData, numWordCData, intDict, numElementDict);
		*/
		int numElement = 0;
                safe_fread(&numElement, sizeof(int), 1, positionFile);
                safe_fread(position, sizeof(int), numElement, positionFile);


                safe_fread(&numElementDData, sizeof(uint64), 1, mismatchFile);
                safe_fread(&numWordCData, sizeof(uint64), 1, mismatchFile);
                safe_fread(&numElementDict, sizeof(uint64), 1, mismatchFile);
                safe_fread(cData, sizeof(word), numWordCData, mismatchFile);
                safe_fread(intDict, sizeof(int), numElementDict, mismatchFile);
                dictDecompression(mismatch, numElementDData, cData, numWordCData, intDict, numElementDict);

                safe_fread(&numElementDData, sizeof(uint64), 1, strandFile);
                safe_fread(&numWordCData, sizeof(uint64), 1, strandFile);
                safe_fread(&numElementDict, sizeof(uint64), 1, strandFile);
                safe_fread(cData, sizeof(word), numWordCData, strandFile);
                safe_fread(charDict, sizeof(char), numElementDict, strandFile);
                dictDecompression(strand, numElementDData, cData, numWordCData, charDict, numElementDict);


                n_read++;

		numReadInBuf = numElementDData; //but do not equal to the read and qual elements
		readId = 0;

		assert(numReadInBuf > 0 && numReadInBuf <= maxNumReadBuf);
	}


	void compress_base() {
		
		/*dict encoding*/
                /*uint64 numWordCData = readLength*numReadInBuf;
                uint64 numElementDict = readLength*numReadInBuf;
                uint64 numElementDData = readLength*numReadInBuf;
                cuda_dictCompression(cData, &numWordCData, charDict, &numElementDict, read, readLength*numReadInBuf);
                safe_fwrite(&numElementDData, sizeof(uint64), 1, readFile);
                safe_fwrite(&numWordCData, sizeof(uint64), 1, readFile);
                safe_fwrite(&numElementDict, sizeof(uint64), 1, readFile);
                safe_fwrite(cData, sizeof(word), numWordCData, readFile);
                safe_fwrite(charDict, sizeof(char), numElementDict, readFile);
                safe_fflush(readFile);*/

		/*diirect encoding*/
		const uint64 numWordCode = getNumBinWord(readLength*numReadInBuf);
		uint64 numElementDData = readLength*numReadInBuf;
		encode(code, read, readLength*numReadInBuf); 
		safe_fwrite(&numElementDData, sizeof(uint64), 1, readFile);
		safe_fwrite(&numWordCode, sizeof(uint64), 1, readFile);
		safe_fwrite(code, sizeof(code_t), numWordCode, readFile);
		safe_fflush(readFile);

	}

	void decompress_base() {

		/*
		uint64 numElementDData, numWordCData, numElementDict;

                safe_fread(&numElementDData, sizeof(uint64), 1, readFile);
                safe_fread(&numWordCData, sizeof(uint64), 1, readFile);
                safe_fread(&numElementDict, sizeof(uint64), 1, readFile);       
                safe_fread(cData, sizeof(word), numWordCData, readFile);
                safe_fread(charDict, sizeof(char), numElementDict, readFile);
                dictDecompression(read, numElementDData, cData, numWordCData, charDict, numElementDict);
		*/

		uint64 numWordCode, numElementDData;
		safe_fread(&numElementDData, sizeof(uint64), 1, readFile);
		safe_fread(&numWordCode, sizeof(uint64), 1, readFile);
		safe_fread(code, sizeof(code_t), numWordCode, readFile);
		decode(read, numElementDData, code, numWordCode);
	}


	void compress_qual() {
                /*uint64 numWordCData = readLength*numReadInBuf;
                uint64 numElementDict = readLength*numReadInBuf;
                uint64 numElementDData = readLength*numReadInBuf;
                cuda_dictCompression(cData, &numWordCData, charDict, &numElementDict, qual, readLength*numReadInBuf);
                safe_fwrite(&numElementDData, sizeof(uint64), 1, qualFile);
                safe_fwrite(&numWordCData, sizeof(uint64), 1, qualFile);
                safe_fwrite(&numElementDict, sizeof(uint64), 1, qualFile);
                safe_fwrite(cData, sizeof(word), numWordCData, qualFile);
                safe_fwrite(charDict, sizeof(char), numElementDict, qualFile);
                safe_fflush(qualFile);
		*/

		uint64 numElement = readLength*numReadInBuf;
		safe_fwrite(&numElement, sizeof(uint64), 1, qualFile);
		safe_fwrite(qual, sizeof(char), numElement, qualFile);
		safe_fflush(qualFile);
	}

	void decompress_qual() {
		uint64 numElement = 0;
		safe_fread(&numElement, sizeof(uint64), 1, qualFile);
		safe_fread(qual, sizeof(char), numElement, qualFile);
	}

	//compress a block of data
	void compress_write() {

		struct timeval tstart, tend;
		gettimeofday(&tstart, NULL);

		n_write++;
		//INIT_TIMER;

		uint64 numWordCData = 0;
		uint64 numElementDict;
		uint64 numElementDData = 0;
	
			
		//START_TIMER;	
		compress_base();
		//END_TIMER;
		//PRINT_TIMER_SEC("read");


		//START_TIMER;
		compress_qual();
		//END_TIMER;
		//PRINT_TIMER_SEC("qual");
		


		//START_TIMER;
                numWordCData = numReadInBuf;
                numElementDict = numReadInBuf;
		numElementDData = numReadInBuf;
                cuda_dictCompression(cData, &numWordCData, intDict, &numElementDict, hit, numReadInBuf);
		safe_fwrite(&numElementDData, sizeof(uint64), 1, hitFile);
                safe_fwrite(&numWordCData, sizeof(uint64), 1, hitFile);
                safe_fwrite(&numElementDict, sizeof(uint64), 1, hitFile);
                safe_fwrite(cData, sizeof(word), numWordCData, hitFile);
                safe_fwrite(intDict, sizeof(int), numElementDict, hitFile);
		safe_fflush(hitFile);
		//END_TIMER;
		//PRINT_TIMER_SEC("hit");


		//START_TIMER;
                numWordCData = numReadInBuf;
                numElementDict = numReadInBuf;
		numElementDData = numReadInBuf;
                cuda_dictCompression(cData, &numWordCData, intDict, &numElementDict, read_len, numReadInBuf);
		safe_fwrite(&numElementDData, sizeof(uint64), 1, readlenFile);
                safe_fwrite(&numWordCData, sizeof(uint64), 1, readlenFile);
                safe_fwrite(&numElementDict, sizeof(uint64), 1, readlenFile);
                safe_fwrite(cData, sizeof(word), numWordCData, readlenFile);
                safe_fwrite(intDict, sizeof(int), numElementDict, readlenFile);
		safe_fflush(readlenFile);
		//END_TIMER;
		//PRINT_TIMER_SEC("read_len");


		//START_TIMER;
                /*numWordCData = numReadInBuf;
                numElementDict = numReadInBuf;
		numElementDData = numReadInBuf;
                cuda_dictCompression(cData, &numWordCData, intDict, &numElementDict, position, numReadInBuf);
		safe_fwrite(&numElementDData, sizeof(uint64), 1, positionFile);
                safe_fwrite(&numWordCData, sizeof(uint64), 1, positionFile);
                safe_fwrite(&numElementDict, sizeof(uint64), 1, positionFile);
                safe_fwrite(cData, sizeof(word), numWordCData, positionFile);
                safe_fwrite(intDict, sizeof(int), numElementDict, positionFile);
		safe_fflush(positionFile);*/
		safe_fwrite(&numReadInBuf, sizeof(int), 1, positionFile);
		safe_fwrite(position, sizeof(int), numReadInBuf, positionFile);
		safe_fflush(positionFile);
		//END_TIMER;
		//PRINT_TIMER_SEC("position");


		//START_TIMER;
                numWordCData = numReadInBuf;
                numElementDict = numReadInBuf;
		numElementDData = numReadInBuf;
                cuda_dictCompression(cData, &numWordCData, intDict, &numElementDict, mismatch, numReadInBuf);
		safe_fwrite(&numElementDData, sizeof(uint64), 1, mismatchFile);
                safe_fwrite(&numWordCData, sizeof(uint64), 1, mismatchFile);
                safe_fwrite(&numElementDict, sizeof(uint64), 1, mismatchFile);
                safe_fwrite(cData, sizeof(word), numWordCData, mismatchFile);
                safe_fwrite(intDict, sizeof(int), numElementDict, mismatchFile);
		safe_fflush(mismatchFile);
		//END_TIMER;
		//PRINT_TIMER_SEC("mismatch");


		//START_TIMER;
                numWordCData = numReadInBuf;
                numElementDict = numReadInBuf;
		numElementDData = numReadInBuf;
                cuda_dictCompression(cData, &numWordCData, charDict, &numElementDict, strand, numReadInBuf);
		safe_fwrite(&numElementDData, sizeof(uint64), 1, strandFile);
                safe_fwrite(&numWordCData, sizeof(uint64), 1, strandFile);
                safe_fwrite(&numElementDict, sizeof(uint64), 1, strandFile);
                safe_fwrite(cData, sizeof(word), numWordCData, strandFile);
                safe_fwrite(charDict, sizeof(char), numElementDict, strandFile);
		safe_fflush(strandFile);
		//END_TIMER;
		//PRINT_TIMER_SEC("strand");


		numReadInBuf = 0;


		gettimeofday(&tend, NULL);
		compwrite_time += getSec(tstart, tend);


		//memset for the next pass
                for(uint64 i = 0; i < readLength*maxNumReadBuf; i++) {
                        read[i] = 'A';
                }
	}

};



// dbSNP information
class Snp_info {
	bool validated;
	bool hapmap_site;
	bool indel_site;
	rate_t * freq; // Arrary of 4 elements recording frequency of ACTG
public:
	Snp_info(){
		validated=hapmap_site=indel_site=false;
		freq = new rate_t [4];
		memset(freq,0,sizeof(rate_t)*4);
	}
	Snp_info(const Snp_info & other) {
		validated = other.validated;
		hapmap_site = other.hapmap_site;
		indel_site = other.indel_site;
		freq = new rate_t [4];
		memcpy(freq, other.freq, sizeof(rate_t)*4);
	}
	~Snp_info(){
		delete [] freq;
	}
	friend std::istringstream& operator>>(std::istringstream & s, Snp_info & snp_form) {
		s>>snp_form.hapmap_site>>snp_form.validated>>snp_form.indel_site>>snp_form.freq[0]>>snp_form.freq[1]>>snp_form.freq[2]>>snp_form.freq[3];
		return s;
	}
	Snp_info & operator=(Snp_info& other) {
		this->validated = other.validated;
		this->hapmap_site = other.hapmap_site;
		this->indel_site = other.indel_site;
		this->freq = new rate_t [4];
		memcpy(this->freq, other.freq, sizeof(rate_t)*4);
		return *this;

	}
	bool is_validated(){
		return validated;
	}
	bool is_hapmap(){
		return hapmap_site;
	}
	bool is_indel(){
		return indel_site;
	}
	rate_t get_freq(char bin_base_2bit) {
		return freq[bin_base_2bit];
	}
};

// Chromosome(Reference) information
class Chr_info {
	ubit32_t len;
	ubit64_t* bin_seq; // Sequence in binary format
	ubit64_t* region_mask;
	ubit64_t* region_win_mask;
	// 4bits for one base: 1 bit dbSNPstatus, 1bit for N, followed two bit of base A: 00, C: 01, T: 10, G:11,
	// Every ubit64_t could store 16 bases
	map<ubit64_t, Snp_info*> dbsnp;
public:
	Chr_info(){
		len = 0;
		bin_seq = NULL;
		region_mask = NULL;
		region_win_mask = NULL;
	};
	Chr_info(const Chr_info & other);
	~Chr_info(){
		delete [] bin_seq;
		delete [] region_mask;
		delete [] region_win_mask;
	}
	ubit32_t length() {
		return len;
	}
	ubit64_t get_bin_base(std::string::size_type pos) {
		return (bin_seq[pos/capacity]>>(pos%capacity*4))&0xF; // All 4 bits
	}
	int binarize(std::string & seq);
	int insert_snp(std::string::size_type pos, Snp_info & new_snp);
	int region_mask_ini();
	bool is_in_region(std::string::size_type pos) {
		return ((region_mask[pos/64]>>(63-pos%64))&1);
	}
	bool is_in_region_win(std::string::size_type pos) {
		pos /= global_win_size; // Calculate in which windows the site is
		//cerr<<pos<<endl;
		//exit(1);
		return ((region_win_mask[pos/64]>>(63-pos%64))&1);
	}
	int set_region(int start, int end);
	Snp_info * find_snp(ubit64_t pos) {
		return dbsnp.find(pos)->second;
	}
	ubit64_t * get_region() {
		return region_mask;
	}
};

typedef std::string Chr_name;
class Genome {
public:
	map<Chr_name, Chr_info*> chromosomes;

	Genome(ifstream & fasta, ifstream & known_snp);
	~Genome();

	bool add_chr(Chr_name &);
	int read_region(std::ifstream & region, Parameter * para);
};

class Prob_matrix {
public:
	rate_t *p_matrix, *p_prior; // Calibration matrix and prior probabilities
	rate_t *base_freq, *type_likely, *type_prob; // Estimate base frequency, conditional probability, and posterior probablity
	rate_t *p_rank, *p_binom; // Ranksum test and binomial test on HETs
	string chr_name;
	Prob_matrix();
	~Prob_matrix();
	int matrix_gen(std::ifstream & alignment, Parameter * para, Genome * genome);
	int matrix_read(std::fstream & mat_in, Parameter * para);
	int matrix_write(std::fstream & mat_out, Parameter * para);
	int prior_gen(Parameter * para);
	int rank_table_gen();

};

class Pos_info {
public:
	unsigned char ori;
	small_int *base_info;
	int pos, *count_uni, *q_sum, depth, dep_uni, repeat_time, *count_all;

	Pos_info(){
		ori = 0xFF;
		//base_info = new small_int [4*2*64*256]; // base info : 4x2x64x64 matrix, base x strand x qual x read_pos
		//memset(base_info,0,sizeof(small_int)*4*2*64*256);
		pos = -1;
		//count_uni = new int [4]; // Count of unique bases
		//memset(count_uni,0,sizeof(int)*4);
		q_sum = new int [4]; // Sum of quality of unique bases
		memset(q_sum,0,sizeof(int)*4);
		depth = 0;
		dep_uni = 0;
		repeat_time = 0;
		count_all = new int [4]; // Count of all bases
		memset(count_all,0,sizeof(int)*4);
	}
	~Pos_info(){
		//delete [] base_info;
		//delete [] count_uni;
		delete [] q_sum;
		delete [] count_all;
	}
};

class Call_win {
public:

	///////////////////
	
	int* count_uni_buf;
	///////////////////
	
        char* oriArr;
        int* posArr;
	int* q_sumArr;
        int* count_allArr;
        int* depthArr;
        int* repeat_timeArr;
        int* dep_uniArr;
	uint32* sparse_base_infoArr;
	uint32* numNonZeroPerSite;
	uint32 maxNumNonZero;
	/////////////////////	
	ubit64_t win_size;
	ubit64_t read_len;
	Pos_info * sites;
	Call_win(ubit64_t read_length, ubit64_t window_size=global_win_size) {
		sites = new Pos_info [window_size+read_length];
		win_size = window_size;
		read_len = read_length;

		//////////////////////////
		
		const uint64 numSite = window_size + read_length;
		count_uni_buf = (int*)malloc(sizeof(int)*COUNTUNI_SIZE*numSite);
		memset(count_uni_buf, 0, sizeof(int)*COUNTUNI_SIZE*numSite);
		for(uint64 i = 0; i < numSite; i++) {
			sites[i].count_uni = count_uni_buf + i*COUNTUNI_SIZE;
		}
		/////////////////////////
		
		maxNumNonZero = 0;
                sparse_base_infoArr = (uint32*)malloc(sizeof(uint32) * numSite * _ini_max_num_non_zero_per_site);
                memset(sparse_base_infoArr, 0xffff, sizeof(uint32) * numSite * _ini_max_num_non_zero_per_site);
                oriArr = (char*)malloc(sizeof(char) * numSite);
                memset(oriArr, 0xFF, sizeof(char) * numSite);

                posArr = (int*)malloc(sizeof(int) * numSite);
                memset(posArr, -1, sizeof(int) * numSite);

                numNonZeroPerSite = (uint32*)malloc(numSite * sizeof(uint32));
                memset(numNonZeroPerSite, 0, sizeof(uint32) * numSite);

                q_sumArr = (int*)malloc(numSite * QSUM_SIZE *sizeof(int));
                memset(q_sumArr,0,sizeof(int) * QSUM_SIZE * numSite);

                count_allArr = (int*)malloc(numSite * COUNTALL_SIZE * sizeof(int));
                memset(count_allArr, 0, sizeof(int) * COUNTALL_SIZE * numSite);

                depthArr = (int*)malloc(numSite * sizeof(int));
                memset(depthArr, 0, numSite * sizeof(int));

                repeat_timeArr = (int*)malloc(numSite * sizeof(int));
                memset(repeat_timeArr, 0, numSite * sizeof(int));

                dep_uniArr = (int*)malloc(numSite *sizeof(int));
                memset(dep_uniArr, 0, numSite * sizeof(int));

                for(uint64 i=0; i < numSite; i++){
                        sites[i].q_sum = q_sumArr + i*4;
                        sites[i].count_all = count_allArr + i*4;
                }

 		////////////////////////////

	}
	~Call_win(){
                free(oriArr);
                free(posArr);
                free(count_uni_buf);
                free(sparse_base_infoArr);
                free(numNonZeroPerSite);
                free(q_sumArr);
                free(count_allArr);
                free(depthArr);
                free(repeat_timeArr);
                free(dep_uniArr);

	}


	int initialize(ubit64_t start);
	int deep_init(ubit64_t start);
	int recycle();
	int quick_recycle();
        int new_call_cns(Chr_name call_name, Chr_info* call_chr, ubit64_t call_length, Prob_matrix * mat, Parameter * para, std::ofstream & consensus, std::ofstream & baseinfo);

	int call_cns(Chr_name call_name, Chr_info* call_chr, ubit64_t call_length, Prob_matrix * mat, Parameter * para, std::ofstream & consensus, std::ofstream & baseinfo);
	int call_cns2(Chr_name call_name, Chr_info* call_chr, ubit64_t call_length, Prob_matrix * mat, Parameter * para, std::ofstream & consensus, std::ofstream & baseinfo);
	int soap2cns(std::ifstream & alignment, std::ofstream & consensus, std::ofstream & baseinfo, Genome * genome, Prob_matrix * mat, Parameter * para);
        int new_soap2cns(std::ifstream & alignment, std::ofstream & consensus, std::ofstream & baseinfo, Genome * genome, Prob_matrix * mat, Parameter * para, const string cnsFileName);

	int snp_p_prior_gen(double * real_p_prior, Snp_info* snp, Parameter * para, char ref);
	double rank_test(Pos_info & info, char best_type, double * p_rank, Parameter * para);
	double normal_value(double z);
	double normal_test(int n1, int n2, double T1, double T2);
	double table_test(double *p_rank, int n1, int n2, double T1, double T2);
};

#endif /*SOAP_SNP_HH_*/
