CC=		gcc
CXX=		g++
CFLAGS=		-g -Wall -O2 #-Wno-unused-but-set-variable
CXXFLAGS=	$(CFLAGS)
DFLAGS=		-DHAVE_PTHREAD #-D_NO_SSE2 #-D_FILE_OFFSET_BITS=64
OBJS=		QSufSort.o bwt_gen.o utils.o bwt.o bwtio.o bwtaln.o bwtgap.o \
			is.o bntseq.o bwtmisc.o bwtindex.o ksw.o stdaln.o simple_dp.o \
			bwaseqio.o bwase.o bwape.o kstring.o cs2nt.o \
			bwase1.o bwase4.o bwaseio1.o bwapese1.o \
			bwape1.o bwape2.o bwape3.o bwape4.o bwapeio1.o \
			bwtsw2_core.o bwtsw2_main.o bwtsw2_aux.o bwt_lite.o \
			bwtsw2_chain.o bamlite.o fastmap.o bwtsw2_pair.o
PROG=		bwa
INCLUDES=	
LIBS=		-lm -lz -lpthread
SUBDIRS=	. bwt_gen

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
.cc.o:
		$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

bwa:$(OBJS) main.c main.h utils.h
		d=`date` ;\
		$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDES) -DBLDDATE="$$d" -c main.c -o main.o ;\
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS) main.o -o $@ $(LIBS)

depend:
		makedepend $(DFLAGS) -Y *.c

clean:
		$(RM) gmon.out core.* *.o a.out $(PROG) *~ *.a

cleanall:
		$(MAKE) clean

# DO NOT DELETE

QSufSort.o: QSufSort.h
bamlite.o: bamlite.h
bntseq.o: bntseq.h main.h utils.h kseq.h
bwape.o: bwtaln.h bwt.h stdaln.h kvec.h bntseq.h utils.h bwase.h bwatpx.h
bwape.o: kstring.h khash.h ksort.h
bwape1.o: bwatpx.h bwtaln.h bwt.h stdaln.h bntseq.h kvec.h kstring.h khash.h
bwape1.o: ksort.h
bwape2.o: bwatpx.h bwtaln.h bwt.h stdaln.h bntseq.h kvec.h kstring.h
bwape3.o: bwatpx.h bwtaln.h bwt.h stdaln.h bntseq.h kvec.h kstring.h
bwape4.o: bwatpx.h bwtaln.h bwt.h stdaln.h bntseq.h kvec.h kstring.h
bwapeio1.o: bwatpx.h bwtaln.h bwt.h stdaln.h bntseq.h kvec.h kstring.h
bwapese1.o: bwatpx.h bwtaln.h bwt.h stdaln.h bntseq.h kvec.h kstring.h
bwase.o: stdaln.h bwase.h bntseq.h bwt.h bwtaln.h utils.h kstring.h bwatpx.h
bwase.o: kvec.h
bwase1.o: bwatpx.h bwtaln.h bwt.h stdaln.h bntseq.h kvec.h kstring.h
bwase4.o: bwatpx.h bwtaln.h bwt.h stdaln.h bntseq.h kvec.h kstring.h
bwaseio1.o: bwatpx.h bwtaln.h bwt.h stdaln.h bntseq.h kvec.h kstring.h
bwaseqio.o: bwtaln.h bwt.h stdaln.h utils.h bamlite.h kseq.h
bwt.o: utils.h bwt.h kvec.h
bwt_gen.o: QSufSort.h
bwt_lite.o: bwt_lite.h
bwtaln.o: bwtaln.h bwt.h stdaln.h bwtgap.h utils.h
bwtgap.o: bwtgap.h bwt.h bwtaln.h stdaln.h
bwtindex.o: bntseq.h bwt.h main.h utils.h
bwtio.o: bwt.h utils.h
bwtmisc.o: bntseq.h utils.h main.h bwt.h
bwtsw2_aux.o: bntseq.h bwt_lite.h utils.h bwtsw2.h bwt.h stdaln.h kstring.h
bwtsw2_aux.o: kseq.h ksort.h
bwtsw2_chain.o: bwtsw2.h bntseq.h bwt_lite.h bwt.h ksort.h
bwtsw2_core.o: bwt_lite.h bwtsw2.h bntseq.h bwt.h kvec.h khash.h ksort.h
bwtsw2_main.o: bwt.h bwtsw2.h bntseq.h bwt_lite.h utils.h
bwtsw2_pair.o: bwt.h bntseq.h bwtsw2.h bwt_lite.h kstring.h ksw.h ksort.h
cs2nt.o: bwtaln.h bwt.h stdaln.h
fastmap.o: bntseq.h bwt.h kvec.h kseq.h
kstring.o: kstring.h
ksw.o: ksw.h
main.o: main.h utils.h
simple_dp.o: stdaln.h utils.h kseq.h
stdaln.o: stdaln.h
utils.o: utils.h
