# Makefile for LESIONS

include ${FSLCONFDIR}/default.mk

PROJNAME = lesions

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB} -I${INC_BOOST}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_ZLIB}

LIBS = -lwarpfns -lbasisfield -lnewimage -lmiscmaths -lfslio -lniftiio -lznz -lnewmat -lutils -lm -lz

RF_OBJS = rf_lda.o
PP_OBJS = preproc.o
LF_OBJS = lesion_filling.o

RUNTCLS = 
XFILES = lesion_filling
TESTXFILES = 
HFILES =
SCRIPTS = 

all:	${XFILES} 

rf_lda:	${RF_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${RF_OBJS} ${LIBS}

preproc:	${PP_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${PP_OBJS} ${LIBS}

lesion_filling:	${LF_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${LF_OBJS} ${LIBS}

