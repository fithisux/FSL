
include $(FSLCONFDIR)/default.mk

PROJNAME = miscvis

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_GD} -I${INC_GDC} -I${INC_PNG} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L${LIB_GD} -L${LIB_GDC} -L${LIB_PNG} -L${LIB_ZLIB}

LIBS = -lutils -lnewimage -lmiscplot -lmiscpic -lmiscmaths -lfslio -lniftiio -lznz -lnewmat -lprob -lm  -lgdc -lgd -lpng -lz

WOBJS = writepng.o wpng.o

DOBJS = readpng.o 

XFILES = slicer overlay fsl_tsplot fsl_boxplot fsl_histogram create_lut wpng pngappend 

SCRIPTS =  slices slicesdir slicesmask slices_summary

RUNTCLS = Renderstats

all:	slicer overlay fsl_tsplot fsl_boxplot fsl_histogram create_lut wpng pngappend 

slicer:	slicer.o
	$(CXX) ${CXXFLAGS} ${LDFLAGS} -o slicer slicer.o  ${LIBS}

overlay: overlay.o
	$(CXX) ${CXXFLAGS} ${LDFLAGS} -o overlay overlay.o ${LIBS}

fsl_tsplot:	fsl_tsplot.o
	$(CXX) ${CXXFLAGS} ${LDFLAGS} -o fsl_tsplot fsl_tsplot.o ${LIBS}

fsl_boxplot:	fsl_boxplot.o
	$(CXX) ${CXXFLAGS} ${LDFLAGS} -o fsl_boxplot fsl_boxplot.o ${LIBS}

fsl_histogram:	fsl_histogram.o
	$(CXX) ${CXXFLAGS} ${LDFLAGS} -o fsl_histogram fsl_histogram.o ${LIBS}

create_lut: create_lut.c
	$(CC) $(CFLAGS) -o create_lut create_lut.c
	@if [ ! -d ${DESTDIR}/etc ] ; then \
	  ${MKDIR} ${DESTDIR}/etc ; \
	  ${CHMOD} g+w ${DESTDIR}/etc ; \
	fi
	@if [ ! -d ${DESTDIR}/etc/luts ] ; then \
	  ${MKDIR} ${DESTDIR}/etc/luts ; \
	  ${CHMOD} g+w ${DESTDIR}/etc/luts ; \
	fi
	${CP} luts+pics/*.gif luts+pics/*.lut $(DESTDIR)/etc/luts
	./create_lut $(DESTDIR)/etc/luts/render

pngappend: pngappend.cc
	$(CXX)  ${CXXFLAGS} ${LDFLAGS} -o pngappend pngappend.cc ${LIBS}

wpng: ${WOBJS}
	$(CC)  ${CXXFLAGS} ${LDFLAGS} -o $@ ${WOBJS} ${LIBS} 