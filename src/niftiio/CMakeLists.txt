ADD_LIBRARY(niftiio nifti1_io.c)

IF(FSL_BUILD)
  INSTALL_FILES(/include/niftiio .*\\.h$)
  INSTALL_TARGETS(/lib niftiio)
ENDIF(FSL_BUILD)
