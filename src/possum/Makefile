# $Id: Makefile,v 1.63 2012/12/26 20:38:37 ivana Exp $
include ${FSLCONFDIR}/default.mk

PROJNAME = possum

#USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_CEPHES}
#USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L${LIB_CEPHES}

#with the next line it will not do the table version automatically
#USRCXXFLAGS = -DNOTABLE

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L${LIB_ZLIB}

LIBS=-lnewimage -lmiscmaths -lfslio -lnewmat -lutils -lprob -lniftiio -lznz -lm -lz

BOBJS=b0calc.o
SOBJS=spharm_rm.o
IOBJS=possum.o possumfns.o 
SIOBJS=signal2image.o
PSOBJS=test_possum.o
AOBJS=pulse.o possumfns.o
SNOBJS=systemnoise.o
PAROBJS=possum_sum.o
RIXOBJS=possum_matrix.o possumfns.o
THEORYOBJ=tcalc.o

RUNTCLS=Possum
XFILES=possum spharm_rm signal2image pulse systemnoise possum_sum b0calc possum_matrix tcalc
TESTXFILES=test_possum 
SCRIPTS=possumX possumX_postproc.sh generate_b0 generate_brain generate_b0calc
MFILES=read_pulse.m write_pulse.m

DBGFLAGS=-DNDEBUG

all: ${XFILES} matlabfiles

b0calc: ${BOBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${BOBJS} ${LIBS} 

spharm_rm: ${SOBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${SOBJS} ${LIBS} 

possum:   ${IOBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${IOBJS} ${LIBS} 

signal2image:   ${SIOBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${SIOBJS} ${LIBS} 

test_possum:   ${PSOBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${PSOBJS} ${LIBS} 

pulse:   ${AOBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${AOBJS} ${LIBS} 

systemnoise:	${SNOBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${SNOBJS} ${LIBS} 

possum_sum: ${PAROBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${PAROBJS} ${LIBS} 

possum_matrix: ${RIXOBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${RIXOBJS} ${LIBS}

tcalc: ${THEORYOBJ}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${THEORYOBJ} ${LIBS}

matlabfiles:
	@if [ ! -d ${DESTDIR}/etc/matlab ] ; then ${MKDIR} -p ${DESTDIR}/etc/matlab ; ${CHMOD} -R g+w ${DESTDIR}/etc ; fi
	${CP} ${MFILES} ${DESTDIR}/etc/matlab
