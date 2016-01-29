include ${FSLCONFDIR}/default.mk

PROJNAME = mm

USRINCFLAGS = -I${INC_NEWRAN} -I${INC_NEWMAT} -I${INC_PROB} -I${INC_GD} -I${INC_GDC} -I${INC_PNG} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWRAN} -L${LIB_NEWMAT} -L${LIB_PROB} -L${LIB_GD} -L${LIB_GDC} -L${LIB_PNG} -L${LIB_ZLIB}

LIBS = -lnewimage -lmiscmaths -lutils -lmiscplot -lmiscpic -lfslio -lniftiio -lznz -lm -lnewmat -lprob -lgdc -lgd -lpng -lz

XFILES = mm

OBJS = mmoptions.o mixture_model.o

#OPTFLAGS =

all:	${XFILES} libmm.a

libmm.a: ${OBJS}
	${AR} -r libmm.a ${OBJS}

mm:${OBJS} mm.o
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} mm.o ${LIBS}

mm2:${OBJS} mm.o
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} mm.o ${LIBS}
