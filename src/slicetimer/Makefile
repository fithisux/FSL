include $(FSLCONFDIR)/default.mk

PROJNAME = slicetimer

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L${LIB_ZLIB}

OBJS = slicetimer.o

LIBS = -lnewimage -lmiscmaths -lfslio -lniftiio -lznz -lnewmat -lprob -lm -lutils -lz

XFILES = slicetimer

all:	${XFILES}

slicetimer:${OBJS}
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${OBJS} ${LIBS}

