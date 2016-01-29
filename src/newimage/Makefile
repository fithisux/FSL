# A Makefile for the image class

include ${FSLCONFDIR}/default.mk

PROJNAME = newimage

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_BOOST} -I${INC_ZLIB} -DEXPOSE_TREACHEROUS
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_ZLIB}

LIBS =-lmiscmaths -lfslio -lniftiio -lutils -lnewmat -lznz -lm -lz

OBJS = lazy.o newimage.o generalio.o newimagefns.o complexvolume.o imfft.o
NOBJS = ${OBJS} costfns.o

TOBJS = ${NOBJS} testprog.o  
SOBJS = ${NOBJS} skeletonapp.o  

TESTXFILES = testprog skeletonapp
XFILES = 

all: libnewimage.a

libnewimage.a: ${NOBJS}
	${AR} -r libnewimage.a ${NOBJS}

skeletonapp: ${SOBJS}
	$(CXX)  ${CXXFLAGS} ${LDFLAGS} -o $@ ${SOBJS} ${LIBS}

testprog: ${TOBJS}
	$(CXX)  ${CXXFLAGS} ${LDFLAGS} -o $@ ${TOBJS} ${LIBS}

