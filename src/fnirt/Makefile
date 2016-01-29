include ${FSLCONFDIR}/default.mk

PROJNAME=fnirt

USRINCFLAGS=-I${INC_NEWMAT} -I${INC_BOOST} 
USRLDFLAGS=-L${LIB_NEWMAT} 

FNIRT_OBJS=fnirt_costfunctions.o intensity_mappers.o matching_points.o fnirtfns.o fnirt.o
INTMAP_OBJS=fnirt_costfunctions.o fnirtfns.o intensity_mappers.o matching_points.o mapintensities.o
INVWARP_OBJS=new_invwarp.o
LIBS=-lwarpfns -lmeshclass -lbasisfield -lnewimage -lmiscmaths -lfslio -lniftiio -lznz -lutils -lnewmat -lm -lz

XFILES=fnirt new_invwarp

all: ${XFILES} config

config:
	@if [ ! -d ${DESTDIR}/etc ] ; then ${MKDIR} ${DESTDIR}/etc ; ${CHMOD} g+w ${DESTDIR}/etc ; fi
	@if [ ! -d ${DESTDIR}/etc/flirtsch ] ; then ${MKDIR} ${DESTDIR}/etc/flirtsch ; ${CHMOD} g+w ${DESTDIR}/etc/flirtsch ; fi
	${CP} -rf fnirtcnf/* ${DESTDIR}/etc/flirtsch/.

fnirt: ${FNIRT_OBJS}
	${CXX} ${LDFLAGS} -o fnirt ${FNIRT_OBJS} ${LIBS}

new_invwarp: ${INVWARP_OBJS}
	${CXX} ${LDFLAGS} -o new_invwarp ${INVWARP_OBJS} ${LIBS}

mapintensities: ${INTMAP_OBJS}
	${CXX} ${LDFLAGS} -o mapintensities ${INTMAP_OBJS} ${LIBS}

