#ifndef __GSNP_IO_H__
#define __GSNP_IO_H__

#include <gsnp_util.h>



struct ReadBuf{
	
	FILE* file;

	ReadBuf(const char* fileName, const uint64 bufSize) {
		safe_fopen(&file, fileName, "rb");
	}

	void close()  {

	}
};


#endif /*__GSNP_IO_H__*/



