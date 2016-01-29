# $Id: Makefile,v 1.5 2009/07/09 10:49:22 mwebster Exp $
include ${FSLCONFDIR}/default.mk

PROJNAME = MVdisc

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT}  -L${LIB_ZLIB} 

LIBS= -lfslvtkio -lmeshclass -lnewimage -lmiscmaths -lfslio -lniftiio -lznz -lnewmat -lutils -lz

SHAPEOBJS=MVdisc.o
 
MODELOBJS=MVdisc.o

all: MVdisc.a 

MVdisc.a: ${SHAPEOBJS}
	${AR} -r libMVdisc.a ${SHAPEOBJS}




