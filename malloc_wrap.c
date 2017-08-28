#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#ifdef USE_MALLOC_WRAPPERS
/* Don't wrap ourselves */
#  undef USE_MALLOC_WRAPPERS
#endif
#include "malloc_wrap.h"

void *wrap_calloc(size_t nmemb, size_t size,
				  const char *file, unsigned int line, const char *func) {
	void *p = calloc(nmemb, size);
	if (NULL == p) {
		fprintf(stderr,
				"[%s] Failed to allocate %zu bytes at %s line %u: %s\n",
				func, nmemb * size, file, line, strerror(errno));
		exit(EXIT_FAILURE);
	}
	return p;
}

void *wrap_malloc(size_t size,
				  const char *file, unsigned int line, const char *func) {
	void *p = malloc(size);
	if (NULL == p) {
		fprintf(stderr,
				"[%s] Failed to allocate %zu bytes at %s line %u: %s\n",
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
				"[%s] Failed to allocate %zu bytes at %s line %u: %s\n",
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
				"[%s] Failed to allocate %zu bytes at %s line %u: %s\n",
				func, strlen(s), file, line, strerror(errno));
		exit(EXIT_FAILURE);
	}
	return p;
}
