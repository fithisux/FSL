
TEMPLATE = app
LANGUAGE = C++
CONFIG += qt warn_off release exceptions qwt fsl boost vtk
#CONFIG += qt warn_on debug exceptions qwt fsl boost vtk
INCLUDEPATH += .
TARGET = fslview
SOURCES += application.cpp \
           version.cpp \
           tracker.cpp \
           main.cpp \
           cursor.cpp \
           bricon.cpp \
           filemanager.cpp \
           imagegroup.cpp \
           overlaylist.cpp \
           overlayinfodialog.cpp \
           briconwidget.cpp \
           cursorwidget.cpp \
           drawwidget.cpp \
           drawsettings.cpp \
           imagedisplaysetting.cpp \
           orthowidget.cpp \
           viewwidget.cpp \
           imagewidget.cpp \
           singlewidget.cpp \
           lightboxwidget.cpp \
           lookuptable.cpp \
           metaimage.cpp \
           imagedata.cpp \
           curvedatalist.cpp \
           graphmanager.cpp \
           imagedatastore.cpp \
           imagebuffer.cpp \
           rect.cpp \
           slicewidget.cpp \
           splashscreen.cpp \
           timeserieswidget.cpp \
           gridserieswidget.cpp \
           singleserieswidget.cpp \
           cubeserieswidget.cpp \
           overlaywidget.cpp \
           histogramwidget.cpp \
           histogramtoolbar.cpp \
           histogramoptionsdialogimpl.cpp \
           properties.cpp \
           propertiesdialogimpl.cpp \
           vtkpropertydialog.cpp \
	   vtktoolbar.cpp \
           shape.cpp \
           command.cpp \
           modelfit.cpp \
	   tsplotcode.cpp \
           featmodel.cpp \
	   maintoolbar.cpp \
	   plotoptions.cpp \
           vtkwidget.cpp

HEADERS += application.h \
           version.h \
	   options.h \
           tracker.h \
           cursor.h \
           bricon.h \
           filemanager.h \
           imagegroup.h \
           overlaylist.h \
           overlayinfodialog.h \
           briconwidget.h \
           cursorwidget.h \
           drawwidget.h \
           drawsettings.h \
           imagedisplaysetting.h \
           orthowidget.h \
           viewwidget.h \
           imagewidget.h \
           singlewidget.h \
           lightboxwidget.h \
           lookuptable.h \
           metaimage.h \
           imagedata.h \
           curvedatalist.h \
           graphmanager.h \
           imagedatastore.h \
           imagebuffer.h \
           rect.h \
           slicewidget.h \
           splashscreen.h \
           timeserieswidget.h \
           gridserieswidget.h \
           singleserieswidget.h \
           cubeserieswidget.h \
           overlaywidget.h \
           histogramwidget.h \
	   histogramtoolbar.h \
           histogramoptionsdialogimpl.h \
           properties.h \
           propertiesdialogimpl.h \
           vtkpropertydialog.h \
	   vtktoolbar.h \
           shape.h \
           command.h \
           modelfit.h \
	   tsplotcode.h \
           featmodel.h \
	   maintoolbar.h \
	   plotoptions.h \
           vtkwidget.h

FORMS += histogramoptionsdialog.ui \
	 histogramtoolbarbase.ui \
         propertiesdialog.ui \
         overlaywidgetbase.ui \
         createmaskdialog.ui \
         overlayinfodialogbase.ui \
	 maintoolbarbase.ui \
	 vtktoolbarbase.ui \
	 vtkpropertydialogbase.ui \
	 drawtoolbarbase.ui \
	 briconwidgetbase.ui

macx{
  QMAKE_LFLAGS_SONAME += -Wl,-install_name,@executable_path/../Frameworks/
  RC_FILE = application.icns
}

mac{
  QMAKE_LFLAGS_SONAME += -Wl,-install_name,@executable_path/../Frameworks/
  RC_FILE = application.icns
  LIBS +=-framework Carbon
}

fsl{
  INCLUDEPATH += ${FSLDEVDIR}/include
  INCLUDEPATH += ${FSLDIR}/include
  INCLUDEPATH += ${FSLDIR}/extras/include/newmat
  LIBS += -L${FSLDEVDIR}/lib -lstorage -lfslio -lniftiio -lznz
  LIBS += -L${FSLDIR}/lib -lmiscmaths -lnewmat -lutils
  LIBS += -L${FSLDIR}/extras/lib
}

boost{
  INCLUDEPATH += ${BOOSTDIR}
}

qwt{
  QMAKE_CXXFLAGS += -DHAVE_QWTSTDVECTORDATA
  INCLUDEPATH += ${QWTDIR}/include
  LIBS += -L${QWTDIR}/lib
  LIBS += -lqwt
}

unix{
  LIBS += -lz
  QMAKE_CXXFLAGS_RELEASE -= -fno-exceptions
#  QMAKE_LFLAGS_RELEASE += -Wa,-rpath=${QWTDIR}/lib
  debug{
    QMAKE_CXXFLAGS -= -O2
  }
}

mem_debug{
  LIBS += -lefence
}

vtk{
  QMAKE_CXXFLAGS += -DHAVE_VTK
  INCLUDEPATH += ${VTKDIR}/include/vtk
  LIBS += -L${VTKDIR}/lib/vtk
  LIBS += -lQVTK
  LIBS += -lvtkCommon
  LIBS += -lvtkGraphics
  LIBS += -lvtkRendering
  LIBS += -lvtkFiltering
  LIBS += -lvtkHybrid
  LIBS += -lvtkWidgets
  LIBS += -lvtkImaging
  LIBS += -lvtkIO
}
