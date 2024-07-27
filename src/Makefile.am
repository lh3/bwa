# for building both the bwa library and the binary

libheaders = utils.h kstring.h ksw.h bwt.h bntseq.h bwa.h bwamem.h malloc_wrap.h QSufSort.h rope.h rle.h ksort.h kbtree.h khash.h kvec.h bwtsw2.h
libsources= utils.c kthread.c kstring.c ksw.c bwt.c bntseq.c bwa.c bwamem.c bwamem_pair.c bwamem_extra.c malloc_wrap.c QSufSort.c bwt_gen.c rope.c rle.c is.c bwtindex.c

appsources = bwashm.c bwase.c bwaseqio.c bwtgap.c bwtaln.c bamlite.c bwape.c kopen.c pemerge.c maxk.c bwtsw2_core.c bwtsw2_main.c bwtsw2_aux.c bwt_lite.c bwtsw2_chain.c fastmap.c bwtsw2_pair.c
appheaders = bwase.h bwtgap.h bwtaln.h bamlite.h bwt_lite.h 

pkginclude_HEADERS = $(libheaders)
lib_LTLIBRARIES = libbwa.la
libbwa_la_SOURCES = $(libheaders) $(libsources)
libbwa_la_LDFLAGS = -version-info 0:0:0

bin_PROGRAMS = bwa bwamem-light
LDADD = libbwa.la

bwa_SOURCES = main.c $(appsources) $(libheaders) $(appheaders)
bwamem_light_SOURCES = example.c bwamem.h kseq.h

