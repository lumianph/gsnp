#ifndef __SAFE_DRIVER_H__
#define __SAFE_DRIVER_H__

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>


inline void* safe_calloc(size_t num, size_t size) {
	void* p = calloc(num, size);
	if(p == NULL) {
		printf("!!!safe_calloc failed, num = %ld, size = %d.\n", num, size);
		exit(1);
	}
	return p;
}

inline void* safe_malloc(size_t size) {
	void* p = malloc(size);
	if(p == NULL) {
		printf("!!!safe_malloc failed, size = %ld.\n", size);
		exit(1);
	}
	return p;
}

inline void safe_free(void* ptr) {
	if(ptr != NULL) {
		free(ptr);
		ptr = NULL;
	}
}

inline void* safe_memset(void* ptr, int value, size_t num) {
	return memset(ptr, value, num);
}

inline void safe_rename( const char* oldName, const char* newName )
{
	int rc = rename(oldName, newName);
	if( rc != 0 )
	{
		printf( "file re-name error: %s -> %s\n", oldName, newName );
		exit(1);
	}
}

inline void safe_fopen(FILE** file, const char* fileName, const char* mode)
{
	*file = fopen( fileName, mode );
	if( (*file) == NULL )
	{
		printf( "cannot open the file: %s\n", fileName );
		exit(1);
	}
};


//delete a file
inline void safe_remove(const char* fileName) {
	if(remove(fileName) != 0) {
		printf("Error deleting file: %s\n", fileName);
		exit(1);
	}
}

inline int safe_fclose(FILE* file)
{
	int rc = fclose( file );
	if( rc != 0 )
	{
		printf( "fclose error\n" );
		exit(1);
	}
	
	return rc;
};


inline int safe_fseek( FILE* file, long offset, int origin )
{
	int rc = fseek( file, offset, origin );
	if( rc != 0 )
	{
		printf( "fseek error\n" );
		exit(1);
	}
	
	return rc;
}

inline size_t safe_fwrite( const void* str, size_t size, size_t count, FILE* file )
{
	size_t rc = fwrite( str, size, count, file );
	if( rc != count )
	{
		printf( "fwrite error: to write, %d; written count %d\n", count, rc );
		exit(1);
	}

	return rc;
};

inline size_t safe_fread( void* dstBuf, size_t size, size_t count, FILE* file )
{
	size_t rc = fread( dstBuf, size, count, file );
	if( rc != count )
	{
		printf( "fread error: to read %d; read count %d\n", count, rc );
		exit(1);
	}

	return rc;
};

inline int safe_fflush( FILE* file )
{
	int rc = fflush( file );
	if( rc != 0 )
	{
		printf( "fflush error.\n" );
		exit(0);
	}
	return rc;
};


#endif //__SAFE_DRIVER_H__



