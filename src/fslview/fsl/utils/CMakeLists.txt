PROJECT(utils)

IF(NOT LIBRARY_OUTPUT_PATH)
  SET(LIBRARY_OUTPUT_PATH ${utils_BINARY_DIR}/lib)
ENDIF(NOT LIBRARY_OUTPUT_PATH)

#SET(utils_SOURCES
#  matches.cc functions.cc usage.cc check.cc parse.cc log.cc time_tracer.cc)
SET(utils_SOURCES
  matches.cc functions.cc usage.cc check.cc parse.cc)

ADD_LIBRARY(utils ${utils_SOURCES})
ADD_EXECUTABLE(opttst opttst.cc)
ADD_DEPENDENCIES(opttst utils)

IF(UNIX)
  TARGET_LINK_LIBRARIES(opttst ${LIBRARY_OUTPUT_PATH}/libutils.a)
ELSE(UNIX)
  TARGET_LINK_LIBRARIES(opttst utils)
ENDIF(UNIX)

IF(FSL_BUILD)
  INSTALL_FILES(/include/utils .*\\.h$)
  INSTALL_TARGETS(/lib utils)
ENDIF(FSL_BUILD)

ENABLE_TESTING()

ADD_TEST(opttst1     ${EXECUTABLE_OUTPUT_PATH}/opttst -d 7 -m mask -r res)
ADD_TEST(opttst2     ${EXECUTABLE_OUTPUT_PATH}/opttst -V)
ADD_TEST(opttst3     ${EXECUTABLE_OUTPUT_PATH}/opttst -verb)
