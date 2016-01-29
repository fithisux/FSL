include ${FSLCONFDIR}/default.mk

PROJNAME = cluster

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_ZLIB} -I${INC_BOOST} 
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L${LIB_ZLIB}

LIBS = -lutils -lwarpfns -lbasisfield -lnewimage -lmiscmaths -lprob -lfslio -lniftiio -lznz -lnewmat -lz

IOBJS = infertest.o infer.o
ILIBS = -lprob

SMOOTHEST = smoothest
SOBJS = smoothest.o


CLUSTER = cluster
COBJS = cluster.o infer.o
CLIBS = -lutils -lwarpfns -lbasisfield -lnewimage -lmiscmaths -lfslio -lniftiio -lznz -lnewmat -lprob -lz

CCOBJS = connectedcomp.o

XFILES = smoothest cluster connectedcomp

SCRIPTS = cluster2html

all: ${XFILES}

smoothest: ${SOBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${SOBJS} ${LIBS}

infertest: ${IOBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${IOBJS} ${ILIBS}

cluster: ${COBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${COBJS} ${CLIBS} 

connectedcomp: ${CCOBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${CCOBJS} ${LIBS} 
