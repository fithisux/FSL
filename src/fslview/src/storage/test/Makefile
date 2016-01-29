include ${FSLCONFDIR}/default.mk

PROJNAME = testStorage

CPPUNIT = ${HOME}/usr/local
INC_CPPUNIT = ${CPPUNIT}/include
LIB_CPPUNIT = ${CPPUNIT}/lib

USRINCFLAGS = -I${HOME}/${FSLMACHTYPE}/include -I${INC_CPPUNIT} -I${INC_BOOST} -I${INC_ZLIB}
USRLDFLAGS  = -L${LIB_CPPUNIT}

GOBJS = genimage.o
HOBJS = testhistfns.o
LIBS = -lstorage -lfslio -lniftiio -lznz -lz

all: testhistfns 

genimage: ${GOBJS}
	${CXX} ${LDFLAGS} -o $@ ${GOBJS} -lcppunit -ldl ${LIBS}

testhistfns: ${HOBJS}
	${CXX} ${LDFLAGS} -o $@ ${HOBJS} ${LIBS}

check: genimage
	@ ./genimage

