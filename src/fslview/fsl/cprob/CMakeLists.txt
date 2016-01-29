SET(LIBPROB_SOURCES
  bdtr.c btdtr.c chdtr.c drand.c fdtr.c gamma.c gdtr.c igam.c
  igami.c incbet.c incbi.c mtherr.c nbdtr.c ndtr.c ndtri.c pdtr.c
  stdtr.c unity.c polevl.c const.c)

ADD_LIBRARY(prob STATIC ${LIBPROB_SOURCES})

IF(FSL_BUILD)
  INSTALL_FILES(/extras/include/libprob .*\\.h$)
  INSTALL_TARGETS(/extras/lib prob)
ENDIF(FSL_BUILD)
