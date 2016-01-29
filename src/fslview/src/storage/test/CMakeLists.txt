PROJECT(FSLVIEW)
SET(testhist_SOURCES testhistfns.cpp)
SET(testimage_SOURCES testimage.cpp)

FIND_PATH(BOOST_DIR boost/shared_ptr.hpp
	/Developer/boost
	/usr/local/boost
	/usr/local/boost_1_33_0
        /usr/include)
INCLUDE_DIRECTORIES(${BOOST_DIR})

FIND_PATH(FSL_DIR etc/fslconf/fsl.sh /usr/local/fsl)
INCLUDE_DIRECTORIES(${FSL_DIR}/include)

INCLUDE_DIRECTORIES(${FSLVIEW_SOURCE_DIR}/fsl)
LINK_DIRECTORIES( ${FSL_DIR}/lib ${FSL_DIR}/extras/lib )

ADD_EXECUTABLE (testhistfns ${testhist_SOURCES})
TARGET_LINK_LIBRARIES (testhistfns storage miscmaths newmat niftiio fslio znz z)

ADD_EXECUTABLE(testimage ${testimage_SOURCES})
TARGET_LINK_LIBRARIES(testimage utils storage miscmaths newmat niftiio fslio znz z)

ENABLE_TESTING()

SET(data_DIR $ENV{HOME}/fslview-test-data)
ADD_TEST(Image6G1 ${EXECUTABLE_OUTPUT_PATH}/testimage -i ${data_DIR}/all_FA_skeletonised -c ${data_DIR}/list1.txt)
ADD_TEST(Image6G2 ${EXECUTABLE_OUTPUT_PATH}/testimage -i ${data_DIR}/all_FA_skeletonised -c ${data_DIR}/list2.txt)

