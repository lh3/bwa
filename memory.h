/**
 * safe calls to alloc function : abort in case of error
 * Pierre Lindenbaum PhD.
 */
#ifndef BWA_MEMORY_H
#define BWA_MEMORY_H
#include <stddef.h>
#include <stdlib.h>

extern void* _safeMalloc(size_t,const char* file,int line);
extern void* _safeCalloc(size_t,size_t,const char* file,int line);
extern void* _safeRealloc(void*,size_t,const char* file,int line);
extern char* _safeStrdup(const char*,const char* file,int line);

#define SAFE_MALLOC(SIZEOF) _safeMalloc(SIZEOF,__FILE__,__LINE__)
#define SAFE_CALLOC(N,SIZEOF) _safeCalloc(N,SIZEOF,__FILE__,__LINE__)
#define SAFE_REALLOC(PTR,SIZEOF) _safeRealloc(PTR,SIZEOF,__FILE__,__LINE__)
#define SAFE_STRDUP(s) _safeStrdup(s,__FILE__,__LINE__)

#ifdef AC_KSEQ_H
#error "insert memory.h before kseq.h"
#endif

#ifdef __AC_KHASH_H
#error "insert memory.h before kash.h"
#endif

#ifdef __AC_KBTREE_H
#error "insert memory.h before kbtree.h"
#endif

/** for khash.h */

#ifndef kcalloc
#define kcalloc(N,Z) SAFE_CALLOC(N,Z)
#endif
#ifndef kmalloc
#define kmalloc(Z) SAFE_MALLOC(Z)
#endif
#ifndef krealloc
#define krealloc(P,Z) SAFE_REALLOC(P,Z)
#endif

#endif

