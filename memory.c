/**
 * safe calls to alloc function : abort in case of error
 * Pierre Lindenbaum PhD.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define TEST_MEMORY(PTR,SIZEOF) if(PTR==NULL) {fprintf(stderr,"[%s:%d]Out of memory(%ld bytes).\n",file,line,SIZEOF); exit(EXIT_FAILURE);}

void* _safeMalloc(size_t n,const char* file,int line)
	{
	void* p=malloc(n);
	TEST_MEMORY(p,n);
	return p;
	}
	
void* _safeCalloc(size_t n,size_t m,const char* file,int line)
	{
	void* p=calloc(n,m);
	TEST_MEMORY(p,n*m);
	return p;
	}
	
void* _safeRealloc(void* ptr,size_t N,const char* file,int line)
	{
	void* p=realloc(ptr,N);
	TEST_MEMORY(p,N);
	return p;
	}
	
char* _safeStrdup(const char* s,const char* file,int line)	
	{
	char* p=strdup(s);
	TEST_MEMORY(p,strlen(s));
	return p;
	}
