AM_CFLAGS = -Wall
LIBBWASRC = utils.c kstring.c ksw.c bwt.c bntseq.c bwa.c bwamem.c bwamem_pair.c malloc_wrap.c
A1 = QSufSort.c bwt_gen.c bwase.c bwaseqio.c bwtgap.c bwtaln.c bamlite.c \
                        is.c bwtindex.c bwape.c kopen.c pemerge.c \
                        bwtsw2_core.c bwtsw2_main.c bwtsw2_aux.c bwt_lite.c \
                        bwtsw2_chain.c fastmap.c bwtsw2_pair.c

noinst_LIBRARIES = libbwa.a
libbwa_a_SOURCES = $(LIBBWASRC)

bin_PROGRAMS = bwa bwamem-lite
bwa_SOURCES = $(A1) main.c
bwa_LDADD = libbwa.a

bwamem_lite_SOURCES = example.c
bwamem_lite_LDADD = libbwa.a
