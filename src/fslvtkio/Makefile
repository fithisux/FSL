# $Id: Makefile,v 1.7 2009/05/15 18:46:21 brian Exp $
include ${FSLCONFDIR}/default.mk

PROJNAME = fslvtkio

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB} -I${HOME}/fsldev/lib
USRLDFLAGS = -L${LIB_NEWMAT}  -L${LIB_ZLIB} 

LIBS= -lmeshclass -lnewimage -lmiscmaths -lfslio -lniftiio -lznz -lnewmat -lutils -lz

SHAPEOBJS=fslvtkio.o
 
MODELOBJS=fslvtkio.o

all: fslvtkio.a 

fslvtkio.a: ${SHAPEOBJS}
	${AR} -r libfslvtkio.a ${SHAPEOBJS}




