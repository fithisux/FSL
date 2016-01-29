SET(Storage_SOURCES
	image.cpp
	volume.cpp
	timeseries.cpp
	imageinfo.cpp
	histogram.cpp
	error.cpp)

ADD_LIBRARY(storage ${Storage_SOURCES})

FIND_PATH(FSL_DIR etc/fslconf/fsl.sh /usr/local/fsl)
INCLUDE_DIRECTORIES(${FSL_DIR}/include)

FIND_PATH(FSL_EXTRAS include/newmat/newmatap.h ${FSL_DIR}/extras)
INCLUDE_DIRECTORIES(${FSL_EXTRAS}/include/newmat)

FIND_PATH(BOOST_DIR boost/shared_ptr.hpp
	/Developer/boost
	/usr/local/boost
	/usr/local/boost_1_33_0
        /usr/include)
INCLUDE_DIRECTORIES(${BOOST_DIR})

IF(FSL_BUILD)
  INSTALL_FILES(/include/storage .*\\.h$)
  INSTALL_FILES(/include/storage .*\\.hpp$)
  INSTALL_FILES(/include/storage .*\\.inc$)
  INSTALL_TARGETS(/lib storage)
ENDIF(FSL_BUILD)

SUBDIRS(test)

