#include <stdarg.h>
#include <stdio.h>
#include "kstring.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

int bwa_kvsprintf(kstring_t *s, const char *fmt, va_list ap)
{
	va_list ap2;
	int l;
	va_copy(ap2, ap);
	l = vsnprintf(s->s + s->l, s->m - s->l, fmt, ap);
	if (l + 1 > s->m - s->l) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
		l = vsnprintf(s->s + s->l, s->m - s->l, fmt, ap2);
	}
	va_end(ap2);
	s->l += l;
	return l;
}

#ifdef KSTRING_MAIN
#include <stdio.h>
int main()
{
	kstring_t *s;
	s = (kstring_t*)calloc(1, sizeof(kstring_t));
	ksprintf(s, "abcdefg: %d", 100);
	printf("%s\n", s->s);
	free(s);
	return 0;
}
#endif
