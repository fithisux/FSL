
FSLCONFDIR=.

include default.mk

PROJNAME = config

install:
	@${MAKE} "DESTDIR=${FSLDEVDIR}" install-general

install-fmrib:
	@${MAKE} "DESTDIR=${FSLDIR}" install-general

install-general:
	${CP} -R -f . ${DESTDIR}/config
	${RM} -rf ${DESTDIR}/config/CVS ${DESTDIR}/config/*/CVS

