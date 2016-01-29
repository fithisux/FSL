# A Makefile for miscplot

include ${FSLCONFDIR}/default.mk

PROJNAME = libvis

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_GD} -I${INC_GDC} -I${INC_PNG} -I${INC_ZLIB} -I${INC_BOOST}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L${LIB_GD} -L${LIB_GDC} -L${LIB_PNG} -L${LIB_ZLIB}

LIBS = -lnewimage -lmiscplot -lmiscpic -lmiscmaths -lfslio -lniftiio -lznz -lnewmat -lprob -lm  -lgdc -lgd -lpng -lz
CXXFLAGS += -std=c++98

PLOTOBJS =  miscplot.o  

PICOBJS  =  miscpic.o

TOBJS = test.o  

TESTXFILES = test

all: miscplot miscpic

miscplot: ${PLOTOBJS}
	${AR} -r libmiscplot.a ${PLOTOBJS}

miscpic: ${PICOBJS}
	${AR} -r libmiscpic.a ${PICOBJS}

test: ${TOBJS}
	$(CXX) ${CXXFLAGS} ${LDFLAGS} -o $@ ${TOBJS} ${LIBS}
