include $(FSLCONFDIR)/default.mk

PROJNAME = misc_c

USRINCFLAGS = -I${INC_PROB} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_PROB} -L${LIB_ZLIB}

IOLIBS = -lfslio -lniftiio -lznz
#LIBS = -lss_8UI $(IOLIBS) -lprob -lm -lz

XFILES = whirlgif hist2prob morph_kernel ptoz tmpnam ztop 

all: $(XFILES)

hist2prob: hist2prob.c
	$(CC) $(CFLAGS) -o hist2prob hist2prob.c -L${LIB_ZLIB} -lm -lz

morph_kernel: morph_kernel.c
	$(CC) -o morph_kernel morph_kernel.c

ptoz: ptoz.cc
	$(CXX) $(CFLAGS) -o ptoz ptoz.cc $(LDFLAGS) -lprob -lm -lz

ztop: ztop.cc
	$(CXX) $(CFLAGS) -o ztop ztop.cc $(LDFLAGS) -lprob -lm -lz

tmpnam: tmpnam.c
	$(CC) $(CFLAGS) -o tmpnam tmpnam.c

whirlgif: gifdecod.c gifencod.c whirlgif.c
	$(CC) $(CFLAGS) -o whirlgif whirlgif.c gifdecod.c gifencod.c
