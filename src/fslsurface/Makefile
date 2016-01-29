# $Id: Makefile,v 1.20 2014/02/24 11:09:05 cowboy Exp $
include ${FSLCONFDIR}/default.mk

PROJNAME = fslsurface

FMLIBS= -lshapeModel -lfslvtkio -lmeshclass -lnewimage -lmiscmaths  -lfslio -lniftiio -lznz -lnewmat -lutils -lprob -lz
isapple = $(findstring apple, $(FSLMACHTYPE))

ifneq ($(isapple),apple)
        FMLIBS += -lGL -lGLU
endif

ifeq ($(isapple),apple)
        applebitness = -arch x86_64
        SYSROOT = -I/usr/X11R6/include
        SYSROOTL =

        ifeq ($(FSLMACHTYPE),apple-darwin8-gcc4.0)
                SYSROOT = -isysroot /Developer/SDKs/MacOSX10.4u.sdk -I/usr/X11R6/include/
                SYSROOTL = -L/Developer/SDKs/MacOSX10.4u.sdk/usr/X11R6/lib/ 
                applebitness = -arch i386
        endif

        ifeq ($(FSLMACHTYPE),apple-darwin10-gcc4.2)
                applebitness = -arch i386 -arch x86_64
        endif

        ifeq ($(isapple),apple)
                ARCHFLAGS =  $(applebitness) $(SYSROOT) 
                ARCHLDFLAGS = -Wl,-search_paths_first ${applebitness} $(SYSROOT) $(SYSROOTL) -framework OpenGL -framework GLUT
        endif 

endif

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_ZLIB} -I${FSLDIR}/include/niftiio 

USRLDFLAGS = -L${LIB_NEWMAT}  -L${LIB_PROB} -L${LIB_ZLIB}


LIBS=-lgiftiio -lexpat -lfirst_lib -lmeshclass

APPLY_ZLIB = -DHAVE_ZLIB


SHAPEOBJS=fslsurface.o 
MODELOBJS=fslsurface.o

BRIANLIBFILES=fslsurfacegl.a
LIBFILES=fslsurface.a fslsurface_backcompat.a
XFILES=fslsurfacemaths

all:  ${LIBFILES} ${XFILES}

fslsurface.a: fslsurfacefns.o fslsurfaceio.o fslsurface.o 
	${AR} -r libfslsurface.a  fslsurfacefns.o fslsurfaceio.o fslsurface.o
fslsurfacegl.a: fslsurfacegl.o fslsurface.a
	${AR} -r libfslsurfacegl.a fslsurfacegl.o
fslsurface_backcompat.a: fslsurface_first.o fslsurface_dataconv.o fslsurface.a
	${AR} -r libfslsurface_backcompat.a fslsurface_first.o fslsurface_dataconv.o
fslsurfacemaths: fslsurface.a fslsurfacemaths.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslsurfacemaths.o fslsurfaceio.o fslsurfacefns.o fslsurface_first.o fslsurface.o ${LIBS} ${FMLIBS}
marching_cubes: fslsurface.a marching_cubes.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ marching_cubes.o fslsurfaceio.o fslsurfacefns.o ${LIBS} ${FMLIBS}
