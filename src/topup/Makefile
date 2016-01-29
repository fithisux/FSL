include ${FSLCONFDIR}/default.mk

PROJNAME=topup

USRINCFLAGS=-I${INC_NEWMAT} -I${INC_BOOST} 
USRLDFLAGS=-L${LIB_NEWMAT} 

OBJS=topup_file_io.o topup_costfunctions.o topupfns.o topup.o
APPLYOBJS=topup_file_io.o displacement_vector.o applytopup.o
TESTOBJS=topup_file_io.o displacement_vector.o test_displacement_vector.o
LIBS=-lwarpfns -lmeshclass -lbasisfield -lnewimage -lmiscmaths -lprob -lfslio -lniftiio -lznz -lutils -lnewmat -lm -lz

XFILES=topup applytopup 
FXFILES=test_displacement_vector

all: ${XFILES} schedule

schedule:
	@if [ ! -d ${DESTDIR}/etc ] ; then ${MKDIR} ${DESTDIR}/etc ; ${CHMOD} g+w ${DESTDIR}/etc ; fi
	@if [ ! -d ${DESTDIR}/etc/flirtsch ] ; then ${MKDIR} ${DESTDIR}/etc/flirtsch ; ${CHMOD} g+w ${DESTDIR}/etc/flirtsch ; fi
	${CP} -rf flirtsch/* ${DESTDIR}/etc/flirtsch/.

topup: ${OBJS}
	${CXX} ${LDFLAGS} -o $@ ${OBJS} ${LIBS}

applytopup: ${APPLYOBJS}
	${CXX} ${LDFLAGS} -o $@ ${APPLYOBJS} ${LIBS}

test_displacement_vector: ${TESTOBJS}
	${CXX} ${LDFLAGS} -o $@ ${TESTOBJS} ${LIBS}
