# $Id: systemvars.mk,v 1.3 2007/07/13 11:00:21 duncan Exp $

# for SHELL, do not change the type of shell - only use Bourne or BASH
SHELL = /bin/sh

# Compiler dependent variables

CC = gcc
CXX = c++

ARCHFLAGS = -mv8 -ffast-math -fomit-frame-pointer
ARCHLDFLAGS = -static
DEPENDFLAGS = -MM

OPTFLAGS = -O6 -fexpensive-optimizations ${ARCHFLAGS}
MACHDBGFLAGS =
GNU_ANSI_FLAGS = -Wall -ansi -pedantic
SGI_ANSI_FLAGS = -ansi -fullwarn
ANSI_FLAGS = ${GNU_ANSI_FLAGS}

INSTALL = ginstall
RM = /bin/rm
CP = /bin/cp
MV = /bin/mv
CHMOD = /bin/chmod
MKDIR = /bin/mkdir
RANLIB = echo
TCLSH = ${FSLDIR}/bin/fsltclsh

