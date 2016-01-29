include ${FSLCONFDIR}/default.mk

PROJNAME = filmbabe

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_GD} -I${INC_GDC} -I${INC_PNG} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L${LIB_GD} -L${LIB_GDC} -L${LIB_PNG} -L${LIB_ZLIB}

LIBS = -lutils -lmiscplot -lmiscpic -lnewimage -lmiscmaths -lnewmat -lfslio -lniftiio -lprob -lgdc -lgd -lpng -lznz -lz -lm

SCRIPTS = filmbabescript

XFILES = filmbabe halfcosbasis

TCLFILES = *.tcl *.ppm

RUNTCLS = Make_flobs

OBJS = halfcosbasisoptions.o filmbabe_manager.o filmbabeoptions.o filmbabe_vb_flobs.o

#OPTFLAGS =

all:	${XFILES}

filmbabe:${OBJS} filmbabe.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} filmbabe.o ${LIBS}

halfcosbasis:halfcosbasisoptions.o halfcosbasis.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ halfcosbasisoptions.o halfcosbasis.o ${LIBS}

