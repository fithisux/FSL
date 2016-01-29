# A Makefile for warpfns library

include ${FSLCONFDIR}/default.mk

PROJNAME = warpfns

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_BOOST} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_ZLIB}

LIBS = -lbasisfield -lmeshclass -lnewimage -lmiscmaths -lfslio -lniftiio -lutils -lnewmat -lznz -lm -lz

LIB_OBJS = warpfns.o fnirt_file_reader.o fnirt_file_writer.o point_list.o

TESTXFILES = test_NewimageCoord2NewimageCoord

all: libwarpfns.a 

libwarpfns.a: ${LIB_OBJS}
	${AR} -r libwarpfns.a ${LIB_OBJS}
	${RANLIB} libwarpfns.a

