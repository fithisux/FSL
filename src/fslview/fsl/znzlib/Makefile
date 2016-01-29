include ${FSLCONFDIR}/default.mk

PROJNAME = znzlib

USRINCFLAGS = -I${INC_ZLIB} -DHAVE_ZLIB

LIBS = -lznz -lm -lz 

OBJS=znzlib.o

SCRIPTS = 
TESTXFILES = testprog

all: libznz.a

test: ${TESTXFILES}

libznz.a: ${OBJS}
	${AR} -r libznz.a ${OBJS}

testprog: libznz.a testprog.c
	${CC} ${CFLAGS} ${LDFLAGS} -o testprog testprog.c ${LIBS}
