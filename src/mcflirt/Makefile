include ${FSLCONFDIR}/default.mk

PROJNAME = mcflirt

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L${LIB_ZLIB}

MOBJS = mcflirt.o Globaloptions.o Log.o

LIBS = -lnewimage -lmiscmaths -lutils -lfslio -lniftiio -lznz -lnewmat -lprob -lm -lz

XFILES = mcflirt 

all:	${XFILES}

mcflirt:${MOBJS}
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${MOBJS} ${LIBS}