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

#ifndef LH3_UTILS_H
#define LH3_UTILS_H

#include <stdio.h>
#include <zlib.h>

#define err_fatal_simple(msg) err_fatal_simple_core(__func__, msg)
#define xopen(fn, mode) err_xopen_core(__func__, fn, mode)
#define xzopen(fn, mode) err_xzopen_core(__func__, fn, mode)
#define xassert(cond, msg) if ((cond) == 0) err_fatal_simple_core(__func__, msg)

#ifdef __cplusplus
extern "C" {
#endif

	void err_fatal(const char *header, const char *fmt, ...);
	void err_fatal_simple_core(const char *func, const char *msg);
	FILE *err_xopen_core(const char *func, const char *fn, const char *mode);
	gzFile err_xzopen_core(const char *func, const char *fn, const char *mode);

#ifdef __cplusplus
}
#endif

#endif
