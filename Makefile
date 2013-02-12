CC=		gcc
CXX=		g++
# CFLAGS=		-ggdb -O0 -Wall -DUSE_MMAP
CFLAGS=		-ggdb -O2 -Wall -DUSE_MMAP
# CFLAGS=		-ggdb -O2 -Wall
CFLAGS+=	`pkg-config --cflags libzmq`
# CFLAGS+=    `pkg-config --cflags libczmq`
CXXFLAGS=	$(CFLAGS)
DFLAGS=		-DHAVE_PTHREAD #-D_FILE_OFFSET_BITS=64
OBJS=		utils.o bwt.o bwtio.o bwtaln.o bwtgap.o is.o \
			bntseq.o bwtmisc.o bwtindex.o stdaln.o simple_dp.o \
			bwaseqio.o bwase.o bwape.o kstring.o cs2nt.o \
			bwtsw2_core.o bwtsw2_main.o bwtsw2_aux.o bwt_lite.o \
			bwtsw2_chain.o bamlite.o bam2bam.o bgzf.o
PROG=		bwa
INCLUDES=	
LIBS=		-lm -lz -lpthread -Lbwt_gen -lbwtgen
LIBS+=		`pkg-config --libs libzmq`
# LIBS+=     	`pkg-config --libs libczmq`
SUBDIRS=	. bwt_gen
prefix=		/usr/local

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
.cc.o:
		$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

install: all
	install -d         $(prefix)/bin
	install -m 775 bwa $(prefix)/bin/
	install -d           $(prefix)/share/man/man1
	install -m 644 bwa.1 $(prefix)/share/man/man1/

lib-recur all-recur clean-recur cleanlocal-recur install-recur:
		@target=`echo $@ | sed s/-recur//`; \
		wdir=`pwd`; \
		list='$(SUBDIRS)'; for subdir in $$list; do \
			cd $$subdir; \
			$(MAKE) CC="$(CC)" CXX="$(CXX)" DFLAGS="$(DFLAGS)" CFLAGS="$(CFLAGS)" \
				INCLUDES="$(INCLUDES)" $$target || exit 1; \
			cd $$wdir; \
		done;

lib:

bwa:lib-recur $(OBJS) main.o
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS) main.o -o $@ $(LIBS)


depend.mk: *.c
		gcc -MM $^ > $@

-include depend.mk

cleanlocal:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a

clean:cleanlocal-recur
