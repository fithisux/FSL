include $(FSLCONFDIR)/default.mk

PROJNAME = qboot

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_NEWRAN} -I${INC_CPROB} -I${INC_PROB} -I${INC_BOOST} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_NEWRAN} -L${LIB_CPROB} -L${LIB_PROB} -L${LIB_ZLIB}

DLIBS = -lnewimage -lutils -lmiscmaths -lnewmat -lnewran -lfslio -lniftiio -lznz -lcprob -lprob -lm -lz

QBOOT=qboot

QBOOTOBJS=qboot.o qbootOptions.o

SGEQBOOT = qboot_parallel qboot_postproc.sh qboot_preproc.sh
SCRIPTS = ${SGEQBOOT} 

XFILES = qboot 


all: qboot
${QBOOT}:	   ${QBOOTOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${QBOOTOBJS} ${DLIBS}