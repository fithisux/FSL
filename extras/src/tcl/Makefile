include ${FSLCONFDIR}/common/vars.mk
include ${FSLCONFDIR}/${FSLMACHTYPE}/systemvars.mk
include ${FSLCONFDIR}/${FSLMACHTYPE}/externallibs.mk

PROJNAME = tcl

SCRIPTS = fsltclsh

all:

install: $(SCRIPTS)
	${MKDIR} -p ${FSLDIR}/bin
	${INSTALL} -m 0775 ${SCRIPTS} ${FSLDIR}/bin
	${MV} ${FSLEXTBIN}/tclsh8.4 ${FSLEXTBIN}/fsltclsh8.4
