include ${FSLCONFDIR}/default.mk

PROJNAME = utils

USRINCFLAGS = -I${INC_NEWMAT}
USRLIBFLAGS = -L${LIB_NEWMAT}

USRCXXFLAGS = 

OptionsTestTarget = opttst
OptionsTestObjects = opttst.o
OptionsTestLibrarys = -lutils

LibraryTarget = libutils.a
LibraryObjects = matches.o functions.o usage.o check.o parse.o log.o time_tracer.o

all: ${LibraryTarget}

${OptionsTestTarget}: ${LibraryTarget} ${OptionsTestObjects}
	${CXX} ${LDFLAGS} -o $@  ${OptionsTestObjects} \
		${OptionsTestLibrarys}

${LibraryTarget}: ${LibraryObjects}
	${RM} -f ${LibraryTarget}
	${AR} -r $@ ${LibraryObjects}
	${RANLIB} ${LibraryTarget}

XFILES = 
AFILES = ${LibraryTarget}

