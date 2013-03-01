CC=			gcc
CFLAGS=		-g -Wall -O2 -msse2
CXXFLAGS=	$(CFLAGS)
AR=			ar
DFLAGS=		-DHAVE_PTHREAD #-D_NO_SSE2 #-D_FILE_OFFSET_BITS=64
LOBJS=		utils.o kstring.o ksw.o kopen.o bwt.o bntseq.o bwa.o bwamem.o bwamem_pair.o
AOBJS=		QSufSort.o bwt_gen.o stdaln.o bwase.o bwaseqio.o bwtgap.o bwtaln.o bamlite.o \
			is.o bwtindex.o bwape.o \
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
		$(CC) $(CFLAGS) $(DFLAGS) $(AOBJS) main.o -o $@ $(LIBS) -L. -lbwa

libbwa.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a

depend:
	( LC_ALL=C ; export LC_ALL; makedepend -Y -- $(CFLAGS) -- *.c )

# DO NOT DELETE THIS LINE -- make depend depends on it.

QSufSort.o: QSufSort.h
bamlite.o: utils.h bamlite.h
bntseq.o: bntseq.h main.h utils.h kseq.h
bwa.o: bntseq.h bwa.h bwt.h ksw.h utils.h kseq.h
bwamem.o: kstring.h utils.h bwamem.h bwt.h bntseq.h bwa.h ksw.h kvec.h
bwamem.o: ksort.h kbtree.h
bwamem_pair.o: kstring.h utils.h bwamem.h bwt.h bntseq.h bwa.h kvec.h ksw.h
bwape.o: bwtaln.h bwt.h stdaln.h kvec.h bntseq.h utils.h bwase.h bwa.h
bwape.o: khash.h
bwase.o: stdaln.h bwase.h bntseq.h bwt.h bwtaln.h utils.h kstring.h bwa.h
bwaseqio.o: bwtaln.h bwt.h stdaln.h utils.h bamlite.h kseq.h
bwt.o: utils.h bwt.h kvec.h
bwt_gen.o: QSufSort.h utils.h
bwt_lite.o: bwt_lite.h utils.h
bwtaln.o: bwtaln.h bwt.h stdaln.h bwtgap.h utils.h bwa.h bntseq.h
bwtgap.o: bwtgap.h bwt.h bwtaln.h stdaln.h utils.h
bwtindex.o: bntseq.h bwt.h main.h utils.h
bwtsw2_aux.o: bntseq.h bwt_lite.h utils.h bwtsw2.h bwt.h stdaln.h kstring.h
bwtsw2_aux.o: bwa.h kseq.h ksort.h
bwtsw2_chain.o: bwtsw2.h bntseq.h bwt_lite.h bwt.h utils.h ksort.h
bwtsw2_core.o: bwt_lite.h bwtsw2.h bntseq.h bwt.h kvec.h utils.h khash.h
bwtsw2_core.o: ksort.h
bwtsw2_main.o: bwt.h bwtsw2.h bntseq.h bwt_lite.h utils.h bwa.h
bwtsw2_pair.o: utils.h bwt.h bntseq.h bwtsw2.h bwt_lite.h kstring.h ksw.h
fastmap.o: bwa.h bntseq.h bwt.h bwamem.h kvec.h utils.h kseq.h
is.o: utils.h
kopen.o: utils.h
kstring.o: kstring.h utils.h
ksw.o: ksw.h utils.h
main.o: main.h utils.h
stdaln.o: stdaln.h utils.h
utils.o: utils.h ksort.h kseq.h
