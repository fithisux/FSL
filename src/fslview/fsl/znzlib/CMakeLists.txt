ADD_LIBRARY(znz STATIC znzlib.c)

IF(FSL_BUILD)
  INSTALL_FILES(/include/znzlib .*\\.h$)
  INSTALL_TARGETS(/lib znz)
ENDIF(FSL_BUILD)
