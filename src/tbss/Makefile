include ${FSLCONFDIR}/default.mk

PROJNAME = tbss

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_CPROB} -I${INC_PROB} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_CPROB} -L${LIB_PROB} -L${LIB_ZLIB}

LIBS = -lnewimage -lmiscmaths -lprob -lfslio -lniftiio -lznz -lnewmat -lutils -lm -lz

XFILES = distancemap tbss_skeleton swap_voxelwise swap_subjectwise

SCRIPTS = tbss_1_preproc tbss_2_reg tbss_3_postreg tbss_4_prestats fsl_reg tbss_non_FA tbss_fill tbss_3_postreg_new tbss_non_FA_new tbss_deproject tbss_sym tbss_x

all: ${XFILES} 

distancemap: distancemap.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ distancemap.cc ${LIBS} 

tbss_skeleton: tbss_skeleton.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ tbss_skeleton.cc ${LIBS} 

swap_voxelwise: swap_voxelwise.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ swap_voxelwise.cc ${LIBS} 

swap_subjectwise: swap_subjectwise.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ swap_subjectwise.cc ${LIBS} 

