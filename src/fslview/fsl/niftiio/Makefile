include ${FSLCONFDIR}/default.mk

PROJNAME = niftiio

USRINCFLAGS = -I${INC_ZLIB} -DHAVE_ZLIB

LIBS = -lm -lznz -lz

OBJS=nifti1_io.o
NOBJS=nifti1_test.o nifti1_io.o
MOBJS=mjtest.o

SCRIPTS = 
TESTXFILES = nifti1_test mjtest

all: libniftiio.a

test: ${TESTXFILES}

libniftiio.a: ${OBJS}
	${AR} -r libniftiio.a ${OBJS}

nifti1_test: ${NOBJS}
	${CC} ${CFLAGS} ${LDFLAGS} -o $@ ${NOBJS} ${LIBS}

mjtest: ${MOBJS}
	${CC} ${CFLAGS} ${LDFLAGS} -o $@ ${MOBJS} ${LIBS}

