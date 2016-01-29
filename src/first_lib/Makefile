# $Id: Makefile,v 1.3 2009/03/23 10:05:21 cowboy Exp $
include ${FSLCONFDIR}/default.mk

PROJNAME = first_lib

USRINCFLAGS = -I${INC_NEWMAT} 

USRLDFLAGS = -L${LIB_NEWMAT} 


LIBS= -lMVdisc -lmeshUtils -lfslvtkIO -lfirst_lib -lshapeModel -lfslvtkio -lmeshclass -lnewimage -lprob -lmiscmaths -lfslio -lniftiio -lznz -lnewmat  -lutils -lz 

XFILES=

SCRIPTS = 

all: ${XFILES} first_lib.a

first_lib.a: first_mesh.o first_newmat_vec.o
	${AR} -r libfirst_lib.a  first_mesh.o  first_newmat_vec.o

