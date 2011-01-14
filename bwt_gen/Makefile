CC=			gcc
CFLAGS=		-g -Wall -O2 -m64 # comment out `-m64' for 32-bit compilation
DFLAGS=		-D_FILE_OFFSET_BITS=64
OBJS=		bwt_gen.o QSufSort.o
INCLUDES=	
VERSION=	0.1.0
LIBS=		
SUBDIRS=

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

lib:libbwtgen.a

libbwtgen.a:$(OBJS)
		$(AR) -cru $@ $(OBJS)

cleanlocal:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a

clean:cleanlocal
