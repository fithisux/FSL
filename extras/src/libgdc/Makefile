include ${FSLCONFDIR}/default.mk

PROJNAME = libgdc
DBGFLAGS = 

USRINCFLAGS = -I${INC_GD} -I${INC_PNG} -I${INC_ZLIB}
ifeq "$(CC)" "gcc"
   GCCVERSIONGTEQ4 := $(shell expr `gcc -dumpversion | awk -F. '{ printf "%i.%i", $$1, $$2 }'` = 4.0)

   ifeq "$(GCCVERSIONGTEQ4)" "1"
      # MUST compile without optimisation as gcc-4.0.1 has a loop unrolling bug (or something)
      USRINCFLAGS += O0
   endif
endif

TEST_PROGRAMS= gdc_samp1 gdc_samp2 gdc_pie_samp

OBJS = gifencode.o price_conv.o gdc.o gdc_pie.o gdchart.o array_alloc.o

LIBS = -lgd -lpng -lz -lm

TESTXFILES = gdc_samp1 gdc_samp2 gdc_pie_samp

all: libgdc.a 

test: $(TEST_PROGRAMS)

libgdc.a: ${OBJS}
	${AR} -r libgdc.a ${OBJS}

gdc_samp1: gdc_samp1.o
	$(CC)  ${CXXFLAGS} ${LDFLAGS} -o $@ gdc_samp1.o ${LIBS} -lgdc

gdc_samp2: gdc_samp2.o
	$(CC)  ${CXXFLAGS} ${LDFLAGS} -o $@ gdc_samp2.o ${LIBS} -lgdc

gdc_pie_samp: gdc_pie_samp.o
	$(CC)  ${CXXFLAGS} ${LDFLAGS} -o $@ gdc_pie_samp.o ${LIBS} -lgdc
