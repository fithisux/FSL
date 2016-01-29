include ${FSLCONFDIR}/common/vars.mk
include ${FSLCONFDIR}/${FSLMACHTYPE}/systemvars.mk
include ${FSLCONFDIR}/${FSLMACHTYPE}/externallibs.mk

PROJNAME = tk

SCRIPTS = fslwish

all:

install: $(SCRIPTS)
	${MKDIR} -p ${FSLDIR}/bin
	${INSTALL} -m 0775 ${SCRIPTS} ${FSLDIR}/bin
	${MV} ${FSLEXTBIN}/wish8.4 ${FSLEXTBIN}/fslwish8.4
