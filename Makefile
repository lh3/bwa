CC=			gcc
CFLAGS=		-g -Wall -O2
AR=			ar
DFLAGS=		-DHAVE_PTHREAD
LOBJS=		utils.o kstring.o ksw.o bwt.o bntseq.o bwa.o bwamem.o bwamem_pair.o
AOBJS=		QSufSort.o bwt_gen.o bwase.o bwaseqio.o bwtgap.o bwtaln.o bamlite.o \
			is.o bwtindex.o bwape.o kopen.o pemerge.o \
			bwtsw2_core.o bwtsw2_main.o bwtsw2_aux.o bwt_lite.o \
			bwtsw2_chain.o fastmap.o bwtsw2_pair.o
PROG=		bwa
INCLUDES=	
LIBS=		-lm -lz -lpthread
SUBDIRS=	.

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

bwa:libbwa.a $(AOBJS) main.o
		$(CC) $(CFLAGS) $(DFLAGS) $(AOBJS) main.o -o $@ -L. -lbwa $(LIBS)

bwamem-lite:libbwa.a example.o
		$(CC) $(CFLAGS) $(DFLAGS) example.o -o $@ -L. -lbwa $(LIBS)

libbwa.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

ksw.o:ksw.h
kstring.o:kstring.h
utils.o:utils.h ksort.h kseq.h
bntseq.o:bntseq.h
bwt.o:bwt.h utils.h
bwa.o:bwa.h bwt.h bntseq.h
bwamem.o:ksw.h kbtree.h ksort.h kvec.h kstring.h utils.h bwamem.h
bwamem_pair.o:ksw.h kvec.h kstring.h utils.h bwamem.h

QSufSort.o:QSufSort.h
bwt_gen.o:QSufSort.h

fastmap.o:bwt.h bwamem.h

bwtaln.o:bwt.h bwtaln.h kseq.h
bwtgap.o:bwtgap.h bwtaln.h bwt.h

bwtsw2_core.o:bwtsw2.h bwt.h bwt_lite.h
bwtsw2_aux.o:bwtsw2.h bwt.h bwt_lite.h
bwtsw2_main.o:bwtsw2.h

clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a
