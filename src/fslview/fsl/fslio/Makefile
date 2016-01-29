include ${FSLCONFDIR}/default.mk

PROJNAME = fslio

USRINCFLAGS = -I${INC_ZLIB}

OBJS=fslio.o

SCRIPTS = remove_ext fsloutputtype imtest imglob imcp imln imrm immv

all: libfslio.a

libfslio.a: ${OBJS}
	${AR} -r libfslio.a ${OBJS}

TESTOBJS = exttest.o
TESTLIBS = -lfslio -lniftiio -lznz -lz

test: all exttest

exttest: ${TESTOBJS}
	${CC} ${LDFLAGS} -o exttest ${TESTOBJS} ${TESTLIBS}

