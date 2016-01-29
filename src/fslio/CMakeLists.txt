ADD_LIBRARY(fslio STATIC fslio.c)

IF(FSL_BUILD)
  INSTALL_FILES(/include/fslio .*\\.h$)
  INSTALL_TARGETS(/lib fslio)
ENDIF(FSL_BUILD)
