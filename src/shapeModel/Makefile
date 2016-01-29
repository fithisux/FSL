# $Id: Makefile,v 1.6 2009/03/20 10:12:28 brian Exp $
include ${FSLCONFDIR}/default.mk

PROJNAME = shapeModel

SHAPEOBJS=shapeModel.o
USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT}  -L${LIB_ZLIB}
 

all: shapeModel.a 

shapeModel.a: ${SHAPEOBJS}
	${AR} -r libshapeModel.a ${SHAPEOBJS}




