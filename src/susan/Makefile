include $(FSLCONFDIR)/default.mk

PROJNAME = susan

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_ZLIB}

IOLIBS = -lfslio -lniftiio -lznz -lm -lz
SSLIBS = -lss_32R $(IOLIBS)


LIBS = -lss_16SI -lfslio -lniftiio -lznz -lm -lz
LIBS++ = -lnewimage -lmiscmaths -lfslio -lnewmat -lutils -lniftiio -lznz -lm -lz


RUNTCLS = Susan

XFILES = susan 

all: susan


susan: susan.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ susan.o ${LIBS++}

susan_smooth_16SI: susan_smooth.c
	$(CC) $(CFLAGS) -DFDT="signed short" -o susan_smooth_16SI susan_smooth.c $(LDFLAGS) $(LIBS)

