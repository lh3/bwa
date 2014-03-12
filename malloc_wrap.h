#ifndef MALLOC_WRAP_H
#define MALLOC_WRAP_H

#include <stdlib.h>  /* Avoid breaking the usual definitions */
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

	void *wrap_calloc(size_t nmemb, size_t size,
					  const char *file, unsigned int line, const char *func);
	void *wrap_malloc(size_t size,
					  const char *file, unsigned int line, const char *func);
	void *wrap_realloc(void *ptr, size_t size,
					   const char *file, unsigned int line, const char *func);
	char *wrap_strdup(const char *s,
					  const char *file, unsigned int line, const char *func);

  void *align_calloc(size_t nmemb, size_t size);
  void *align_malloc(size_t size);

#ifdef __cplusplus
}
#endif

#ifdef USE_MALLOC_WRAPPERS
#  ifdef calloc
#    undef calloc
#  endif
#  define calloc(n, s)  wrap_calloc( (n), (s), __FILE__, __LINE__, __func__)

#  ifdef malloc
#    undef malloc
#  endif
#  define malloc(s)     wrap_malloc( (s),      __FILE__, __LINE__, __func__)

#  ifdef realloc
#    undef realloc
#  endif
#  define realloc(p, s) wrap_realloc((p), (s), __FILE__, __LINE__, __func__)

#  ifdef strdup
#    undef strdup
#  endif
#  define strdup(s)     wrap_strdup( (s),      __FILE__, __LINE__, __func__)

#else if defined(USE_ALIGNED_MALLOC)
#  ifdef calloc
#    undef calloc
#  endif
#  define calloc(n, s)  align_calloc( (n), (s) )

#  ifdef malloc
#    undef malloc
#  endif
#  define malloc(s)     align_malloc( (s) )

#endif /* USE_MALLOC_WRAPPERS */

#endif /* MALLOC_WRAP_H */
