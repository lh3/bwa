cmake_minimum_required(VERSION 2.8)
project(bwa)

set(CFLAGS
    -g
    -Wall
    -Wno-unused-function
    -O2)

set(WRAP_MALLOC
    -DUSE_MALLOC_WRAPPERS)

set(LIBS
    m
    z
    pthread)

exec_program(uname
    ARGS -s
    OUTPUT_VARIABLE target)

if(target STREQUAL "Linux")
    set(LIBS
        rt
        ${LIBS})
endif()

add_definitions(-DHAVE_PTHREAD ${WRAP_MALLOC})
string(REPLACE ";" " " CMAKE_C_FLAGS "${CFLAGS}")

set(LOBJS_FILES
    utils.c
    kthread.c
    kstring.c
    ksw.c
    bwt.c
    bntseq.c
    bwa.c
    bwamem.c
    bwamem_pair.c
    bwamem_extra.c
    malloc_wrap.c
    QSufSort.c
    bwt_gen.c
    rope.c
    rle.c
    is.c
    bwtindex.c)

set(AOBJS_FILES
    bwashm.c
    bwase.c
    bwaseqio.c
    bwtgap.c
    bwtaln.c
    bamlite.c
    bwape.c
    kopen.c
    pemerge.c
    maxk.c
    bwtsw2_core.c
    bwtsw2_main.c
    bwtsw2_aux.c
    bwt_lite.c
    bwtsw2_chain.c
    fastmap.c
    bwtsw2_pair.c)


add_library(bwa_static ${LOBJS_FILES})

add_executable(bwa ${AOBJS_FILES} main.c)
target_link_libraries(bwa bwa_static) 

add_executable(bwamem-lite example.c)
target_link_libraries(bwamem-lite bwa_static) 

target_link_libraries(bwa ${LIBS})
target_link_libraries(bwamem-lite ${LIBS})
