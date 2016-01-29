include ${FSLCONFDIR}/default.mk

PROJNAME = storage

USRINCFLAGS = -I${INC_BOOST} -I${INC_ZLIB}

HFILES = *.h *.inc histogramfns.hpp

OBJS = image.o volume.o timeseries.o imageinfo.o histogram.o error.o

all: libstorage.a

libstorage.a: ${OBJS}
	${AR} -r libstorage.a ${OBJS}



