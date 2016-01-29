PROJECT(FSLVIEW)

SET(utility_srcs
  preferences.cpp
  atlas.cpp
  filemanager.cpp
  cursor.cpp
  lookuptable.cpp
  clusterdata.cpp)

SET(srcs
  main.cpp
  preferences.cpp
  atlas.cpp
  filemanager.cpp
  logger.cpp
  application.cpp
  assistantclient.cpp
  splashscreen.cpp 
  version.cpp
  tracker.cpp
  cursor.cpp
  bricon.cpp
  imagegroup.cpp 
  overlaylist.cpp 
  metaimage.cpp 
  imagedata.cpp 
  drawsettings.cpp 
  imagedisplaysetting.cpp 
  lookuptable.cpp
  curvedatalist.cpp 
  imagedatastore.cpp 
  imagebuffer.cpp 
  graphmanager.cpp 
  briconwidget.cpp 
  cursorwidget.cpp 
  drawwidget.cpp 
  imagewidget.cpp 
  viewwidget.cpp 
  orthowidget.cpp 
  singlewidget.cpp 
  lightboxwidget.cpp 
  slicewidget.cpp 
  sliceview.cpp 
  timeserieswidget.cpp 
  singleserieswidget.cpp 
  overlaywidget.cpp 
  histogramwidget.cpp 
  clusterbrowser.cpp
  histogramtoolbar.cpp 
  overlayinfodialog.cpp
  meshoptionsdialog.cpp
  atlasoptionsdialog.cpp
  atlasselectiondialog.cpp
  viewoptionsdialog.cpp
  histogramoptionsdialogimpl.cpp 
  properties.cpp 
  propertiesdialogimpl.cpp 
  vtkpropertydialog.cpp 
  vtktoolbar.cpp 
  timeseriestoolbar.cpp
  talairachwidget.cpp
  shape.cpp 
  rect.cpp 
  command.cpp 
  modelfit.cpp
  clusterdata.cpp
  tsplotcode.cpp 
  maintoolbar.cpp 
  modetoolbar.cpp 
  vtkmeshsurface.cpp
  vtkwidget.cpp)

QT_WRAP_CPP( ui moc_srcs
  application.h
  assistantclient.h
  overlayinfodialog.h
  meshoptionsdialog.h
  atlasoptionsdialog.h
  atlasselectiondialog.h
  viewoptionsdialog.h
  imagewidget.h
  viewwidget.h
  sliceview.h
  orthowidget.h
  singlewidget.h
  lightboxwidget.h
  briconwidget.h
  cursorwidget.h
  drawwidget.h
  slicewidget.h
  splashscreen.h
  timeserieswidget.h
  talairachwidget.h
  clusterbrowser.h
  singleserieswidget.h
  maintoolbar.h
  modetoolbar.h
  histogramwidget.h
  vtkwidget.h
  overlaywidget.h)

QT_WRAP_UI(ui ui_hdrs ui_srcs
  sliceviewbase.ui
  applicationwindowbase.ui
  clusterbrowserbase.ui
  histogramoptionsdialog.ui 
  histogramtoolbarbase.ui 
  meshoptionsdialogbase.ui 
  atlasoptionsdialogbase.ui 
  atlasselectiondialogbase.ui 
  viewoptionsdialogbase.ui 
  talairachwidgetbase.ui 
  propertiesdialog.ui 
  overlaywidgetbase.ui 
  createmaskdialog.ui 
  overlayinfodialogbase.ui 
  timeseriestoolbarbase.ui 
  smartpendialogbase.ui
  maintoolbarbase.ui 
  modetoolbarbase.ui 
  vtktoolbarbase.ui 
  vtkpropertydialogbase.ui 
  drawtoolbarbase.ui 
  cursorwidgetbase.ui
  briconwidgetbase.ui
  timeserieswindowbase.ui)

FIND_PACKAGE( Qt3 REQUIRED )
FIND_PACKAGE( VTK REQUIRED )
FIND_PATH( BOOST_DIR boost/shared_ptr.hpp
	   /Developer/boost
	   /usr/local/boost
	   /usr/local/boost_1_33_0 )
FIND_PATH( QWT_DIR include/qwt_wheel.h 
	   /Developer/qwt 
	   /usr/local/qwt
	   /usr/local/qwt-4.2.0)
IF(NOT WINDOWS)
  FIND_PATH( FSL_DIR etc/fslconf/fsl.sh
	   /usr/local/fsl )
ENDIF(NOT WINDOWS)

ADD_DEFINITIONS( ${QT_DEFINITIONS} )

IF(VTK_DIR)
  INCLUDE( ${VTK_USE_FILE} )
  ADD_DEFINITIONS( -DHAVE_VTK )
ENDIF(VTK_DIR)

SET(HAS_SYSLOG ON CACHE BOOL "Use syslog found on many Unix systems")

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h )

INCLUDE_DIRECTORIES( ${QT_INCLUDE_DIR} ${QT_INCLUDE_DIR}/qwt )
INCLUDE_DIRECTORIES( ${QWT_DIR}/include )
INCLUDE_DIRECTORIES( ${BOOST_DIR} )
INCLUDE_DIRECTORIES( ${FSL_DIR}/include )
INCLUDE_DIRECTORIES( ${FSL_DIR}/extras/include/newmat )
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} )

LINK_DIRECTORIES( ${X11_LIBRARY_DIR} )
LINK_DIRECTORIES( ${VTK_LIBRARY_DIRS} )
LINK_DIRECTORIES( ${QWT_DIR}/lib )
LINK_DIRECTORIES( ${FSL_DIR}/lib ${FSL_DIR}/extras/lib )
#LINK_DIRECTORIES( ${FSLVIEW_BINARY_DIR}/lib )

ADD_LIBRARY( ui ${ui_srcs} ${moc_srcs} )
ADD_LIBRARY( utilities ${utility_srcs} )

ADD_EXECUTABLE( atlasquery atlasquery.cc tracker.cpp )
TARGET_LINK_LIBRARIES( atlasquery ${QT_LIBRARIES} utilities storage miscmaths newmat utils fslio niftiio znz z )

IF(APPLE)
  SET( MACOSX_BUNDLE_ICON_FILE application.icns )
  INSTALL(FILES application.icns DESTINATION bin/fslview.app/Contents/Resources)
  SET( MACOSX_BUNDLE_COPYRIGHT "Copyright (c) 2005, University of Oxford" )
  ADD_EXECUTABLE( fslview MACOSX_BUNDLE ${srcs} )
  FILE( GLOB vtklibs ${VTK_DIR}/*.dylib* )
  FILE( GLOB qwtlibs ${QWT_DIR}/lib/*.dylib* )
  FILE( GLOB qtlibs ${QT_INCLUDE_DIR}/../lib/*.dylib* )
  INSTALL(FILES ${vtklibs} ${qtlibs} ${qwtlibs} DESTINATION bin/fslview.app/Contents/Frameworks)
#  INSTALL(DIRECTORY $ENV{QTDIR}/bin/assistant.app DESTINATION lib)
#   INSTALL(TARGETS fslview atlasquery 
# 	RUNTIME DESTINATION /bin
# 	LIBRARY DESTINATION /lib
# 	ARCHIVE DESTINATION /lib )
ELSE (APPLE)
  ADD_EXECUTABLE( fslview ${srcs} )
  FILE( GLOB vtklibs ${VTK_DIR}/../libvtk*.so* ${VTK_DIR}/../libQVTK.so* )
  FILE( GLOB osmesa /usr/lib*/libOSMesa*.so* )
  FILE( GLOB qwtlibs ${QWT_DIR}/lib/*.so* )
  FILE( GLOB qtlibs ${QT_QT_LIBRARY}* ${QT_QASSISTANTCLIENT_LIBRARY}* )
  FILE( GLOB assistant ${QT_INCLUDE_DIR}/../bin/assistant )
  INSTALL(FILES ${vtklibs} ${osmesa} ${qtlibs} ${qwtlibs} DESTINATION lib)
#  INSTALL(FILES ${assistant} DESTINATION lib)
#   INSTALL(TARGETS fslview atlasquery RUNTIME DESTINATION /bin )
ENDIF(APPLE)

TARGET_LINK_LIBRARIES( fslview utilities ui storage)
TARGET_LINK_LIBRARIES( fslview ${QT_LIBRARIES} fslio niftiio znz z )
TARGET_LINK_LIBRARIES( fslview miscmaths newmat utils )
TARGET_LINK_LIBRARIES( fslview qassistantclient qwt )
TARGET_LINK_LIBRARIES( fslview vtkHybrid vtkWidgets QVTK )

INSTALL(TARGETS fslview atlasquery 
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib )
IF(UNIX)
#  TARGET_LINK_LIBRARIES( fslview ${FSLVIEW_BINARY_DIR}/lib/libfslio.a )
#   INSTALL_FILES ( /lib ${VTK_LIBRARY_DIR}/.*\\.so$ )
ENDIF(UNIX)
IF(APPLE)
#   FILE( GLOB vtklibs ${VTK_DIR}/../*.dylib* )
#   INSTALL_FILES ( /bin/fslview.app/Contents/Frameworks FILES ${vtklibs} )
ENDIF(APPLE)
IF(WINDOWS)
    LINK_DIRECTORIES( ${FSLVIEW_BINARY_DIR}/lib/Debug )
ENDIF(WINDOWS)

# Now for the test harness

CREATE_TEST_SOURCELIST(Tests
  atlastests.cpp
  atlastst.cpp
  agtest.cpp)

REMOVE(srcs main.cpp)

#ADD_EXECUTABLE( atlastests ${srcs} ${Tests} )
#TARGET_LINK_LIBRARIES(atlastests ui storage
#  fslio niftiio znz z
#  miscmaths newmat utils qassistantclient qwt QVTK vtkHybrid vtkWidgets )

#SET(TestsToRun ${Tests})
#REMOVE(TestsToRun atlastests.cpp)

#FOREACH (test ${TestsToRun})
#  GET_FILENAME_COMPONENT( TName ${test} NAME_WE )
#  ADD_TEST(${TName} ${EXECUTABLE_OUTPUT_PATH}/atlastests ${TName})
#ENDFOREACH (test)
