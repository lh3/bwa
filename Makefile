CC=			gcc
CXX=		g++
CFLAGS=		-g -Wall -O2
CXXFLAGS=	$(CFLAGS)
AR=			ar
DFLAGS=		-DHAVE_PTHREAD #-D_NO_SSE2 #-D_FILE_OFFSET_BITS=64
LOBJS=		bwa.o bamlite.o utils.o bwt.o bwtio.o bwtaln.o bwtgap.o bntseq.o stdaln.o \
			bwaseqio.o bwase.o kstring.o
AOBJS=		QSufSort.o bwt_gen.o \
			is.o bwtmisc.o bwtindex.o ksw.o simple_dp.o \
			bwape.o cs2nt.o \
			bwtsw2_core.o bwtsw2_main.o bwtsw2_aux.o bwt_lite.o \
			bwtsw2_chain.o fastmap.o bwtsw2_pair.o
PROG=		bwa
INCLUDES=	
LIBS=		-lm -lz -lpthread
SUBDIRS=	.

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
.cc.o:
		$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

bwa:libbwa.a $(AOBJS) main.o
		$(CC) $(CFLAGS) $(DFLAGS) $(AOBJS) main.o -o $@ -L. -lbwa $(LIBS)

libbwa.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a

QSufSort.o: QSufSort.h
bamlite.o: bamlite.h utils.h
bntseq.o: bntseq.h kseq.h main.h utils.h
bwa.o: bntseq.h bwa.h bwt.h bwtaln.h bwtgap.h stdaln.h utils.h
bwape.o: bntseq.h bwase.h bwt.h bwtaln.h khash.h ksort.h kvec.h stdaln.h
bwape.o: utils.h
bwase.o: bntseq.h bwase.h bwt.h bwtaln.h kstring.h stdaln.h utils.h
bwaseqio.o: bamlite.h bwt.h bwtaln.h kseq.h stdaln.h utils.h
bwt.o: bwt.h kvec.h utils.h
bwt_gen.o: QSufSort.h utils.h
bwt_lite.o: bwt_lite.h utils.h
bwtaln.o: bwt.h bwtaln.h bwtgap.h stdaln.h utils.h
bwtgap.o: bwt.h bwtaln.h bwtgap.h stdaln.h utils.h
bwtindex.o: bntseq.h bwt.h main.h utils.h
bwtio.o: bwt.h utils.h
bwtmisc.o: bntseq.h bwt.h main.h utils.h
bwtsw2_aux.o: bntseq.h bwt.h bwt_lite.h bwtsw2.h kseq.h ksort.h kstring.h
bwtsw2_aux.o: stdaln.h utils.h
bwtsw2_chain.o: bntseq.h bwt.h bwt_lite.h bwtsw2.h ksort.h utils.h
bwtsw2_core.o: bntseq.h bwt.h bwt_lite.h bwtsw2.h khash.h ksort.h kvec.h
bwtsw2_core.o: utils.h
bwtsw2_main.o: bntseq.h bwt.h bwt_lite.h bwtsw2.h utils.h
bwtsw2_pair.o: bntseq.h bwt.h bwt_lite.h bwtsw2.h kstring.h ksw.h utils.h
cs2nt.o: bwt.h bwtaln.h stdaln.h utils.h
fastmap.o: bntseq.h bwt.h kseq.h kvec.h utils.h
is.o: utils.h
kstring.o: kstring.h utils.h
ksw.o: ksw.h utils.h
main.o: main.h utils.h
simple_dp.o: kseq.h stdaln.h utils.h
stdaln.o: stdaln.h utils.h
utils.o: utils.h
