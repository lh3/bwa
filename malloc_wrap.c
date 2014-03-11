#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#ifdef USE_MALLOC_WRAPPERS
/* Don't wrap ourselves */
#  undef USE_MALLOC_WRAPPERS
#endif
#include "malloc_wrap.h"

void *align_calloc(size_t nmemb, size_t size) {
  size_t n = nmemb*size;
  void *p;
  if (posix_memalign(&p,64,n)) p=NULL;
  else memset(p,0,n);
  return p;
}

void *align_malloc(size_t size) {
  void *p;
  if (posix_memalign(&p,64,size)) p=NULL;
  return p;
}

void *wrap_calloc(size_t nmemb, size_t size,
				  const char *file, unsigned int line, const char *func) {
#ifdef USE_ALIGNED_MALLOC
  void *p = align_calloc(nmemb, size);
#else
	void *p = calloc(nmemb, size);
#endif
	if (NULL == p) {
		fprintf(stderr,
				"[%s] Failed to allocate %zd bytes at %s line %u: %s\n",
				func, nmemb * size, file, line, strerror(errno));
		exit(EXIT_FAILURE);
	}
	return p;
}

void *wrap_malloc(size_t size,
				  const char *file, unsigned int line, const char *func) {
#ifdef USE_ALIGNED_MALLOC
  void *p = align_malloc(size);
#else
	void *p = malloc(size);
#endif
	if (NULL == p) {
		fprintf(stderr,
				"[%s] Failed to allocate %zd bytes at %s line %u: %s\n",
				func, size, file, line, strerror(errno));
		exit(EXIT_FAILURE);
	}
	return p;
}

void *wrap_realloc(void *ptr, size_t size,
				   const char *file, unsigned int line, const char *func) {
	void *p = realloc(ptr, size);
	if (NULL == p) {
		fprintf(stderr,
				"[%s] Failed to allocate %zd bytes at %s line %u: %s\n",
				func, size, file, line, strerror(errno));
		exit(EXIT_FAILURE);
	}
	return p;
}

char *wrap_strdup(const char *s,
				  const char *file, unsigned int line, const char *func) {
	char *p = strdup(s);
	if (NULL == p) {
		fprintf(stderr,
				"[%s] Failed to allocate %zd bytes at %s line %u: %s\n",
				func, strlen(s), file, line, strerror(errno));
		exit(EXIT_FAILURE);
	}
	return p;
}
