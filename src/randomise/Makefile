include ${FSLCONFDIR}/default.mk

PROJNAME = randomise

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L${LIB_ZLIB}

LIBS = -lnewimage -lmiscmaths -lutils -lm -lnewmat -lfslio -lniftiio -lznz -lprob -lz

XFILES = randomise fdr unconfound

SCRIPTS = design_ttest2 randomise_parallel tfce_support randomise_combine setup_masks 

R_OBJS = randomise.o ranopts.o

F_OBJS = fdr.o

#OPTFLAGS =

all:	${XFILES}

randomise:  ${R_OBJS}
	    ${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@ ${R_OBJS} ${LIBS}

unconfound: unconfound.o
	     ${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@ unconfound.o ${LIBS}

fdr: ${F_OBJS}
	    ${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@ ${F_OBJS} ${LIBS}
