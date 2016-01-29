# A Makefile for MELODIC

include ${FSLCONFDIR}/default.mk

OPTFLAGS = -O3 -Wno-deprecated -ggdb
OPTFLAGS_alphaev6-dec-osf5.0-gcc2.95.2 = -O3 -mieee -mfp-trap-mode=sui

PROJNAME = melodic

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_GD} -I${INC_GDC} -I${INC_PNG} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L${LIB_GD} -L${LIB_GDC} -L${LIB_PNG} -L${LIB_ZLIB}

LIBS = -lutils -lnewimage -lmiscplot -lmiscpic -lmiscmaths -lfslio -lniftiio -lznz -lnewmat -lprob -lm  -lgdc -lgd -lpng -lz

TEST_OBJS = test.o 

GGMIX_OBJS = ggmix.o

FSL_GLM_OBJS = melhlprfns.o fsl_glm.o

FSL_SBCA_OBJS = melhlprfns.o fsl_sbca.o

FSL_MVLM_OBJS = melhlprfns.o fsl_mvlm.o

FSL_BWTF_OBJS = melhlprfns.o fsl_bwtf.o

FSL_REGFILT_OBJS = melhlprfns.o fsl_regfilt.o

FSL_SCHURPROD_OBJS = melhlprfns.o fsl_schurprod.o

MELODIC_OBJS =  meloptions.o melhlprfns.o melgmix.o meldata.o melpca.o melica.o melreport.o melodic.o 

TESTXFILES = test

XFILES = fsl_glm fsl_sbca fsl_mvlm fsl_regfilt fsl_schurprod fsl_bwtf melodic

OTHERFILES =

BUILD_SBCA = $(shell [ -f fsl_sbca.cc ]; echo $$?)
ifeq ($(BUILD_SBCA), 0)
	OTHERFILES+=fsl_sbca
endif

BUILD_MVLM = $(shell [ -f fsl_mvlm.cc ]; echo $$?)
ifeq ($(BUILD_MVLM), 0)
        OTHERFILES+=fsl_mvlm
endif

BUILD_BWTF = $(shell [ -f fsl_bwtf.cc ]; echo $$?)
ifeq ($(BUILD_BPTF), 0)
        OTHERFILES+=fsl_bwtf
endif

BUILD_SCHUR = $(shell [ -f fsl_schurprod.cc ]; echo $$?)
ifeq ($(BUILD_SCHUR), 0)
	OTHERFILES+=fsl_schurprod
endif

XFILES+=${OTHERFILES}


RUNTCLS = Melodic

SCRIPTS = melodicreport dual_regression

all: ggmix fsl_regfilt fsl_glm melodic ${OTHERFILES}

ggmix: ${GGMIX_OBJS}
	${AR} -r libggmix.a ${GGMIX_OBJS} 

melodic: ${MELODIC_OBJS}
	$(CXX) ${CXXFLAGS} ${LDFLAGS} -o $@ ${MELODIC_OBJS} ${LIBS}

test: ${TEST_OBJS}
	$(CXX) ${CXXFLAGS} ${LDFLAGS} -o $@ ${TEST_OBJS} ${LIBS}

fsl_glm: ${FSL_GLM_OBJS}
	$(CXX) ${CXXFLAGS} ${LDFLAGS} -o $@ ${FSL_GLM_OBJS} ${LIBS}

fsl_sbca: ${FSL_SBCA_OBJS}
		$(CXX) ${CXXFLAGS} ${LDFLAGS} -o $@ ${FSL_SBCA_OBJS} ${LIBS}

fsl_schurprod: ${FSL_SCHURPROD_OBJS}
				$(CXX) ${CXXFLAGS} ${LDFLAGS} -o $@ ${FSL_SCHURPROD_OBJS} ${LIBS}

fsl_mvlm: ${FSL_MVLM_OBJS}
		$(CXX) ${CXXFLAGS} ${LDFLAGS} -o $@ ${FSL_MVLM_OBJS} ${LIBS}

fsl_bwtf: ${FSL_BWTF_OBJS}
		$(CXX) ${CXXFLAGS} ${LDFLAGS} -o $@ ${FSL_BWTF_OBJS} ${LIBS}

fsl_regfilt: ${FSL_REGFILT_OBJS}
	$(CXX) ${CXXFLAGS} ${LDFLAGS} -o $@ ${FSL_REGFILT_OBJS} ${LIBS}

