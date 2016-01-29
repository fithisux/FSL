# $Id: externallibs.mk,v 1.3 2014/02/10 17:05:45 mwebster Exp $

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
LIB_ZLIB = /usr/lib
INC_ZLIB = /usr/include

# QT library
QTDIR = /sw
LIB_QT = ${QTDIR}/lib
INC_QT = ${QTDIR}/include

# BOOST library
BOOSTDIR = ${FSLEXTINC}/boost
LIB_BOOST = ${BOOSTDIR}
INC_BOOST = ${BOOSTDIR}

# QWT library
QWTDIR = /usr/local/qwt
INC_QWT = ${QWTDIR}/include
LIB_QWT = ${QWTDIR}/lib
 
# FFTW3 library
LIB_FFTW3 = ${FSLEXTLIB}
INC_FFTW3 = ${FSLEXTINC}/fftw3

# VTK library
VTKDIR_INC = /usr/local/vtk-universal/include/vtk-5.3/
VTKDIR_LIB = /usr/local/vtk-universal/lib/vtk-5.3/
