include $(FSLCONFDIR)/default.mk

PROJNAME = siena

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_ZLIB}

#LIBS = -lss_32R -lfslio -lniftiio -lznz -lm -lz
LIBCC = -lnewimage -lmiscmaths -lfslio -lniftiio -lznz -lnewmat -lutils -lm -lz

SN_OBJS = surface_norm.o
CS_OBJS = cylsamp.o
GT_OBJS = groupttest.o

XFILES = siena_diff
TESTXFILES = groupttest surface_norm cylsamp 
SCRIPTS = siena siena_flirt siena_cal sienax siena_flow2std viena_quant viena_createpng

all:	${XFILES}

# DON'T REMOVE THE -O0!!!! siena_diff doesn't work with optimisation
siena_diff: siena_diff.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -O0 -o siena_diff siena_diff.cc ${LIBCC}

surface_norm: ${SN_OBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${SN_OBJS} ${LIBCC}

cylsamp: ${CS_OBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${CS_OBJS} ${LIBCC}

groupttest: ${GT_OBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${GT_OBJS} ${LIBCC}

