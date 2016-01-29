#
# Makfile for objects in the BASISFIELD namespace
#

include ${FSLCONFDIR}/default.mk

PROJNAME = basisfield

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_BOOST} 
USRLDFLAGS = -L${LIB_NEWMAT} 

OBJS = dctfield.o splinefield.o basisfield.o splines.o

LIBS = -lnewimage -lmiscmaths -lutils -lfslio -lniftiio -lznz -lz -lnewmat -lprob -lm

all: libbasisfield.a

libbasisfield.a: ${OBJS}
	${AR} -r $@ ${OBJS}
	${RANLIB} libbasisfield.a

