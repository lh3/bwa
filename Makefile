CC=			gcc
# CC=			clang --analyze
CFLAGS=		-g -Wall -Wno-unused-function -O2
WRAP_MALLOC=-DUSE_MALLOC_WRAPPERS
AR=			ar
DFLAGS=		-DHAVE_PTHREAD $(WRAP_MALLOC)
LOBJS=		src/utils.o \
			src/kthread.o \
			src/kstring.o \
			src/ksw.o \
			src/bwt.o \
			src/bntseq.o \
			src/bwa.o \
			src/bwamem.o \
			src/bwamem_pair.o \
			src/bwamem_extra.o \
			src/malloc_wrap.o \
			src/QSufSort.o \
			src/bwt_gen.o \
			src/rope.o \
			src/rle.o \
			src/is.o \
			src/bwtindex.o

AOBJS=		src/bwashm.o \
			src/bwase.o \
			src/bwaseqio.o \
			src/bwtgap.o \
			src/bwtaln.o \
			src/bamlite.o \
			src/bwape.o \
			src/kopen.o \
			src/pemerge.o \
			src/maxk.o \
			src/bwtsw2_core.o \
			src/bwtsw2_main.o \
			src/bwtsw2_aux.o \
			src/bwt_lite.o \
			src/bwtsw2_chain.o \
			src/fastmap.o \
			src/bwtsw2_pair.o

PROG=		bwa
INCLUDES=	
LIBS=		-lm -lz -lpthread
SUBDIRS=	.

ifeq ($(shell uname -s),Linux)
	LIBS += -lrt
endif

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

bwa:libbwa.a $(AOBJS) src/main.o
		$(CC) $(CFLAGS) $(DFLAGS) $(AOBJS) src/main.o -o $@ -L. -lbwa $(LIBS)

bwamem-lite:libbwa.a src/example.o
		$(CC) $(CFLAGS) $(DFLAGS) src/example.o -o $@ -L. -lbwa $(LIBS)

libbwa.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

clean:
		rm -f gmon.out src/*.o a.out $(PROG) *~ *.a

depend:
	( LC_ALL=C ; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c )

# DO NOT DELETE THIS LINE -- make depend depends on it.

src/QSufSort.o: src/QSufSort.h
src/bamlite.o: src/bamlite.h src/malloc_wrap.h
src/bntseq.o: src/bntseq.h src/utils.h src/kseq.h src/malloc_wrap.h src/khash.h
src/bwa.o: src/bntseq.h src/bwa.h src/bwt.h src/ksw.h src/utils.h src/kstring.h src/malloc_wrap.h src/kvec.h
src/bwa.o: src/kseq.h
src/bwamem.o: src/kstring.h src/malloc_wrap.h src/bwamem.h src/bwt.h src/bntseq.h src/bwa.h src/ksw.h src/kvec.h
src/bwamem.o: src/ksort.h src/utils.h src/kbtree.h
src/bwamem_extra.o: src/bwa.h src/bntseq.h src/bwt.h src/bwamem.h src/kstring.h src/malloc_wrap.h
src/bwamem_pair.o: src/kstring.h src/malloc_wrap.h src/bwamem.h src/bwt.h src/bntseq.h src/bwa.h src/kvec.h
src/bwamem_pair.o: src/utils.h src/ksw.h
src/bwape.o: src/bwtaln.h src/bwt.h src/kvec.h src/malloc_wrap.h src/bntseq.h src/utils.h src/bwase.h src/bwa.h
src/bwape.o: src/ksw.h src/khash.h
src/bwase.o: src/bwase.h src/bntseq.h src/bwt.h src/bwtaln.h src/utils.h src/kstring.h src/malloc_wrap.h
src/bwase.o: src/bwa.h src/ksw.h
src/bwaseqio.o: src/bwtaln.h src/bwt.h src/utils.h src/bamlite.h src/malloc_wrap.h src/kseq.h
src/bwashm.o: src/bwa.h src/bntseq.h src/bwt.h
src/bwt.o: src/utils.h src/bwt.h src/kvec.h src/malloc_wrap.h
src/bwt_gen.o: src/QSufSort.h src/malloc_wrap.h
src/bwt_lite.o: src/bwt_lite.h src/malloc_wrap.h
src/bwtaln.o: src/bwtaln.h src/bwt.h src/bwtgap.h src/utils.h src/bwa.h src/bntseq.h src/malloc_wrap.h
src/bwtgap.o: src/bwtgap.h src/bwt.h src/bwtaln.h src/malloc_wrap.h
src/bwtindex.o: src/bntseq.h src/bwa.h src/bwt.h src/utils.h src/rle.h src/rope.h src/malloc_wrap.h
src/bwtsw2_aux.o: src/bntseq.h src/bwt_lite.h src/utils.h src/bwtsw2.h src/bwt.h src/kstring.h
src/bwtsw2_aux.o: src/malloc_wrap.h src/bwa.h src/ksw.h src/kseq.h src/ksort.h
src/bwtsw2_chain.o: src/bwtsw2.h src/bntseq.h src/bwt_lite.h src/bwt.h src/malloc_wrap.h src/ksort.h
src/bwtsw2_core.o: src/bwt_lite.h src/bwtsw2.h src/bntseq.h src/bwt.h src/kvec.h src/malloc_wrap.h
src/bwtsw2_core.o: src/khash.h src/ksort.h
src/bwtsw2_main.o: src/bwt.h src/bwtsw2.h src/bntseq.h src/bwt_lite.h src/utils.h src/bwa.h
src/bwtsw2_pair.o: src/utils.h src/bwt.h src/bntseq.h src/bwtsw2.h src/bwt_lite.h src/kstring.h
src/bwtsw2_pair.o: src/malloc_wrap.h src/ksw.h
src/example.o: src/bwamem.h src/bwt.h src/bntseq.h src/bwa.h src/kseq.h src/malloc_wrap.h
src/fastmap.o: src/bwa.h src/bntseq.h src/bwt.h src/bwamem.h src/kvec.h src/malloc_wrap.h src/utils.h src/kseq.h
src/is.o: src/malloc_wrap.h
src/kopen.o: src/malloc_wrap.h
src/kstring.o: src/kstring.h src/malloc_wrap.h
src/ksw.o: src/ksw.h src/malloc_wrap.h
src/main.o: src/kstring.h src/malloc_wrap.h src/utils.h
src/malloc_wrap.o: src/malloc_wrap.h
src/maxk.o: src/bwa.h src/bntseq.h src/bwt.h src/bwamem.h src/kseq.h src/malloc_wrap.h
src/pemerge.o: src/ksw.h src/kseq.h src/malloc_wrap.h src/kstring.h src/bwa.h src/bntseq.h src/bwt.h src/utils.h
src/rle.o: src/rle.h
src/rope.o: src/rle.h src/rope.h
src/utils.o: src/utils.h src/ksort.h src/malloc_wrap.h src/kseq.h
