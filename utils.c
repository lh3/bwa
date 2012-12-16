/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <errno.h>
#include <sys/resource.h>
#include <sys/time.h>
#include "utils.h"

FILE *err_xopen_core(const char *func, const char *fn, const char *mode)
{
	FILE *fp = 0;
	if (strcmp(fn, "-") == 0)
		return (strstr(mode, "r"))? stdin : stdout;
	if ((fp = fopen(fn, mode)) == 0) {
		err_fatal(func, "fail to open file '%s' : %s", fn, strerror(errno));
	}
	return fp;
}
FILE *err_xreopen_core(const char *func, const char *fn, const char *mode, FILE *fp)
{
	if (freopen(fn, mode, fp) == 0) {
		err_fatal(func, "fail to open file '%s' : %s", fn, strerror(errno));
	}
	return fp;
}
gzFile err_xzopen_core(const char *func, const char *fn, const char *mode)
{
	gzFile fp;
	if (strcmp(fn, "-") == 0) {
		fp = gzdopen(fileno((strstr(mode, "r"))? stdin : stdout), mode);
		/* According to zlib.h, this is the only reason gzdopen can fail */
		if (!fp) err_fatal(func, "Out of memory");
		return fp;
	}
	if ((fp = gzopen(fn, mode)) == 0) {
		err_fatal(func, "fail to open file '%s' : %s", fn, errno ? strerror(errno) : "Out of memory");
	}
	return fp;
}

void err_fatal(const char *header, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	fprintf(stderr, "[%s] ", header);
	vfprintf(stderr, fmt, args);
	fprintf(stderr, "\n");
	va_end(args);
	exit(EXIT_FAILURE);
}

void err_fatal_core(const char *header, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	fprintf(stderr, "[%s] ", header);
	vfprintf(stderr, fmt, args);
	fprintf(stderr, " Abort!\n");
	va_end(args);
	abort();
}

void _err_fatal_simple(const char *func, const char *msg)
{
	fprintf(stderr, "[%s] %s\n", func, msg);
	exit(EXIT_FAILURE);
}

void _err_fatal_simple_core(const char *func, const char *msg)
{
	fprintf(stderr, "[%s] %s Abort!\n", func, msg);
	abort();
}

size_t err_fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream)
{
    size_t ret = fwrite(ptr, size, nmemb, stream);
    if (ret != nmemb) 
    {
        _err_fatal_simple("fwrite", strerror(errno));
    }
    return ret;
}

size_t err_fread_noeof(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
	size_t ret = fread(ptr, size, nmemb, stream);
	if (ret != nmemb)
	{
		_err_fatal_simple("fread", ferror(stream) ? strerror(errno) : "Unexpected end of file");
	}
	return ret;
}

int err_gzread(gzFile file, void *ptr, unsigned int len)
{
	int ret = gzread(file, ptr, len);

	if (ret < 0)
	{
		int errnum = 0;
		const char *msg = gzerror(file, &errnum);
		_err_fatal_simple("gzread", Z_ERRNO == errnum ? strerror(errno) : msg);
	}

	return ret;
}

int err_fseek(FILE *stream, long offset, int whence)
{
	int ret = fseek(stream, offset, whence);
	if (0 != ret)
	{
		_err_fatal_simple("fseek", strerror(errno));
	}
	return ret;
}

long err_ftell(FILE *stream)
{
	long ret = ftell(stream);
	if (-1 == ret)
	{
		_err_fatal_simple("ftell", strerror(errno));
	}
	return ret;
}

int err_printf(const char *format, ...) 
{
    va_list arg;
    int done;

    va_start(arg, format);
    done = vfprintf(stdout, format, arg);
    int saveErrno = errno;
    va_end(arg);

    if (done < 0) 
    {
        _err_fatal_simple("vfprintf(stdout)", strerror(saveErrno));
    }
    return done;
}

int err_fprintf(FILE *stream, const char *format, ...) 
{
    va_list arg;
    int done;

    va_start(arg, format);
    done = vfprintf(stream, format, arg);
    int saveErrno = errno;
    va_end(arg);

    if (done < 0) 
    {
        _err_fatal_simple("vfprintf", strerror(saveErrno));
    }
    return done;
}

int err_fflush(FILE *stream) 
{
    int ret = fflush(stream);
    if (ret != 0) 
    {
        _err_fatal_simple("fflush", strerror(errno));
    }
    return ret;
}

int err_fclose(FILE *stream) 
{
    int ret = fclose(stream);
    if (ret != 0) 
    {
        _err_fatal_simple("fclose", strerror(errno));
    }
    return ret;
}

void *err_calloc(size_t nmemb, size_t size, const char *file, unsigned int line, const char *func)
{
	void *p = calloc(nmemb, size);
	if (NULL == p)
	{
		err_fatal(func, "Failed to allocate %zd bytes at %s line %u: %s\n", nmemb * size, file, line, strerror(errno));
	}
	return p;
}

void *err_malloc(size_t size, const char *file, unsigned int line, const char *func)
{
	void *p = malloc(size);
	if (NULL == p)
	{
		err_fatal(func, "Failed to allocate %zd bytes at %s line %u: %s\n", size, file, line, strerror(errno));
	}
	return p;
}

void *err_realloc(void *ptr, size_t size, const char *file, unsigned int line, const char *func)
{
	void *p = realloc(ptr, size);
	if (NULL == p)
	{
		err_fatal(func, "Failed to allocate %zd bytes at %s line %u: %s\n", size, file, line, strerror(errno));
	}
	return p;
}

char *err_strdup(const char *s, const char *file, unsigned int line, const char *func)
{
	char *p = strdup(s);

	if (NULL == p)
	{
		err_fatal(func, "Failed to allocate %zd bytes at %s line %u: %s\n", strlen(s), file, line, strerror(errno));
	}
	return p;
}

double cputime()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}
