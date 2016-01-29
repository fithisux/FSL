
# for SHELL, do not change the type of shell - only use Bourne or BASH
SHELL = /bin/sh

# System dependent commands (NB: the first two are the most platform dependent)

INSTALL = ginstall -p
RANLIB = ranlib

RM = /bin/rm
CP = /bin/cp
MV = /bin/mv
CHMOD = /bin/chmod
MKDIR = /bin/mkdir
TCLSH = ${FSLDIR}/bin/fsltclsh

# Compiler dependent variables

CC = gcc
CXX = c++
CSTATICFLAGS = -static
CXXSTATICFLAGS = -static

ARCHFLAGS = 

DEPENDFLAGS = -MM

OPTFLAGS =  -O3 -fexpensive-optimizations ${ARCHFLAGS}
MACHDBGFLAGS =
GNU_ANSI_FLAGS = -Wall -ansi -pedantic
SGI_ANSI_FLAGS = -ansi -fullwarn
ANSI_FLAGS = ${GNU_ANSI_FLAGS}


