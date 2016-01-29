# $Id: externallibs.mk,v 1.3 2008/11/26 14:15:14 mwebster Exp $

# External Library and Include Paths

FSLEXTLIB=${FSLDIR}/extras/lib
FSLEXTINC=${FSLDIR}/extras/include
FSLEXTBIN=${FSLDIR}/extras/bin

# CEPHES library
LIB_CEPHES = ${FSLEXTLIB}
INC_CEPHES = ${FSLEXTINC}/cephes

# GD library
LIB_GD = ${FSLEXTLIB}
INC_GD = ${FSLEXTINC}

# GDC library
LIB_GDC = ${FSLEXTLIB}
INC_GDC = ${FSLEXTINC}/libgdc

# GSL library
LIB_GSL = ${FSLEXTLIB}
INC_GSL = ${FSLEXTINC}/gsl

# PNG library
LIB_PNG = ${FSLEXTLIB}
INC_PNG = ${FSLEXTINC}

# PROB library
LIB_PROB = ${FSLEXTLIB}
INC_PROB = ${FSLEXTINC}/libprob

# CPROB library
LIB_CPROB = ${FSLEXTLIB}
INC_CPROB = ${FSLEXTINC}/libcprob

# NEWMAT library
LIB_NEWMAT = ${FSLEXTLIB}
INC_NEWMAT = ${FSLEXTINC}/newmat

# NEWRAN library
LIB_NEWRAN = ${FSLEXTLIB}
INC_NEWRAN = ${FSLEXTINC}/newran

# ZLIB library
LIB_ZLIB = ${FSLEXTLIB}
INC_ZLIB = ${FSLEXTINC}

# BOOST library
BOOSTDIR = ${FSLEXTINC}/boost
LIB_BOOST = ${BOOSTDIR}
INC_BOOST = ${BOOSTDIR}

# QT library
QTDIR = /usr/lib/qt3
LIB_QT = ${QTDIR}/lib
INC_QT = ${QTDIR}/include

# QWT library
QWTDIR = /usr/local/qwt
LIB_QWT = ${QWTDIR}/lib
INC_QWT = ${QWTDIR}/include
 
# FFTW3 library
LIB_FFTW3 = ${FSLEXTLIB}
INC_FFTW3 = ${FSLEXTINC}/fftw3

# VTK library
VTKDIR_INC = /usr/local/vtk/include/vtk-5.0/
VTKDIR_LIB = /usr/local/vtk/lib/
