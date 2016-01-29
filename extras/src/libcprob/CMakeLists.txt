IF(FSL_BUILD)
  INSTALL_FILES(/extras/include/libcprob .*\\.h$)
  INSTALL_TARGETS(/extras/lib cprob)
ENDIF(FSL_BUILD)
