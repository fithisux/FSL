IF(FSL_BUILD)
  INSTALL_FILES(/extras/include/libprob .*\\.h$)
  INSTALL_TARGETS(/extras/lib prob)
ENDIF(FSL_BUILD)
