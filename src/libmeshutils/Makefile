# $Id: Makefile,v 1.19 2012/07/23 13:25:20 mwebster Exp $
include ${FSLCONFDIR}/default.mk

PROJNAME = meshUtils

LD_LIBRARY_PATH=${FSLDIR}/lib

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB} -I${INC_PROB} -I${INC_BOOST} 
USRLDFLAGS =  -L${LIB_PROB} -L${LIB_NEWMAT} -L${LIB_ZLIB}

BUILT_VTKCONV = $(shell [ -f ../../libfslvtkconv* ]; echo $$?)

ifeq ($(BUILT_VTKCONV), 0)
	USRINCFLAGS+=-I${VTKDIR_INC}
endif

ifeq ($(BUILT_VTKCONV), 0)
        USRLDFLAGS+=-L${VTKDIR_LIB}
endif

LIBS+=-lwarpfns -lbasisfield -lfslvtkio -lMVdisc

ifeq ($(BUILT_VTKCONV), 0)
        LIBS+=-lfslvtkconv 
endif

LIBS+=-lfirst_lib -lshapeModel -lmeshclass -lnewimage -lmiscmaths -lprob -lfslio -lniftiio -lznz -lutils -lnewmat  -lz

SWOBJS=meshUtils.o
RMUOBJS=meshUtils.o run_mesh_utils.o
XFILES=run_mesh_utils

all: meshUtils.a ${XFILES}

meshUtils.a: ${SWOBJS} 
	${AR} -r libmeshUtils.a  ${SWOBJS}
	@${MAKE} libinstall

run_mesh_utils: ${RMUOBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${RMUOBJS} ${LIBS}


