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
#include <sys/time.h>
#include <time.h>
#include <zlib.h>
#include "utils.h"

extern time_t _prog_start;

FILE *err_xopen_core(const char *func, const char *fn, const char *mode)
{
	FILE *fp = 0;
	if (strcmp(fn, "-") == 0)
		return (strstr(mode, "r"))? stdin : stdout;
	if ((fp = fopen(fn, mode)) == 0) {
		fprintf(stderr, "[%s] fail to open file '%s'. Abort!\n", func, fn);
		abort();
	}
	// setvbuf(fp, NULL, _IOFBF, 1048576);
	return fp;
}
FILE *err_xreopen_core(const char *func, const char *fn, const char *mode, FILE *fp)
{
	if (freopen(fn, mode, fp) == 0) {
		fprintf(stderr, "[%s] fail to open file '%s': ", func, fn);
		perror(NULL);
		fprintf(stderr, "Abort!\n");
		abort();
	}
	// setvbuf(fp, NULL, _IOFBF, 1048576);
	return fp;
}
gzFile err_xzopen_core(const char *func, const char *fn, const char *mode)
{
	gzFile fp;
	if (strcmp(fn, "-") == 0)
		return gzdopen(fileno((strstr(mode, "r"))? stdin : stdout), mode);
	if ((fp = gzopen(fn, mode)) == 0) {
		fprintf(stderr, "[%s] fail to open file '%s'. Abort!\n", func, fn);
		abort();
	}
	// gzbuffer(fp, 524288);
	return fp;
}
void err_fatal(const char *header, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	fprintf(stderr, "[%s] ", header);
	vfprintf(stderr, fmt, args);
	fprintf(stderr, " Abort!\n");
	va_end(args);
	abort();
}

void err_fatal_simple_core(const char *func, const char *msg)
{
	fprintf(stderr, "[%s] %s Abort!\n", func, msg);
	abort();
}

clock_t clock(void)
{
	clock_t clks = 0;
	struct timeval st;
	time_t time_now;

	// mck - use wall time ...

	gettimeofday(&st, NULL);
	time_now = st.tv_sec * 1000000L + (time_t)st.tv_usec;
	clks = (clock_t)(((double)(time_now - _prog_start) / 1000000.0) * (double)CLOCKS_PER_SEC);

	return clks;
}

int getmaxrss(long *maxrsskb)
{
  int len = 0;
  int srtn = 0;
  char procf[257] = { "" };
  FILE *fp = NULL;
  char line[2001] = { "" };
  char crap[2001] = { "" };
  char units[2001] = { "" };
  long maxrss = 0L;

  if(maxrsskb == NULL){
    return -1;
  }

  sprintf(procf,"/proc/%d/status",getpid());

  fp = fopen(procf, "r");
  if(fp == NULL){
    return -1;
  }

  while(fgets(line, 2000, fp) != NULL){
    if(strncasecmp(line,"VmPeak:",7) == 0){
      len = (int)strlen(line);
      line[len-1] = '\0';
      srtn = sscanf(line,"%s%ld%s",crap,&maxrss,units);
      if(srtn == 2){
        *maxrsskb = maxrss / 1024L;
      }else if(srtn == 3){
        if( (strcasecmp(units,"B") == 0) || (strcasecmp(units,"BYTES") == 0) ){
          *maxrsskb = maxrss / 1024L;
        }else if( (strcasecmp(units,"k") == 0) || (strcasecmp(units,"kB") == 0) ){
          *maxrsskb = maxrss * 1L;
        }else if( (strcasecmp(units,"m") == 0) || (strcasecmp(units,"mB") == 0) ){
          *maxrsskb = maxrss * 1024L;
        }else if( (strcasecmp(units,"g") == 0) || (strcasecmp(units,"gB") == 0) ){
          *maxrsskb = maxrss * 1024L * 1024L;
        }else{
          *maxrsskb = maxrss * 1L;
        }
      }
      break;
    }
  }

  fclose(fp);
 
  return 0;
}
