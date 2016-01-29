# for SHELL, do not change the type of shell - only use Bourne or BASH
SHELL = /bin/sh

# Compiler dependent variables

CC = cc
CXX = c++
CSTATICFLAGS = 
CXXSTATICFLAGS = 

ARCHFLAGS = -arch ppc -arch ppc64 -arch i386 -arch x86_64
ARCHLDFLAGS = -Wl,-search_paths_first -arch ppc -arch ppc64 -arch i386 -arch x86_64
PER_ARCH_CFLAGS_i386 =  -msse3
PER_ARCH_CFLAGS_ppc = 
PER_ARCH_CFLAGS_ppc64 = -mcpu=G5
PER_ARCH_CFLAGS_x86_64= -msse3

DEPENDFLAGS = -MM

OPTFLAGS =  -O3 -fexpensive-optimizations
MACHDBGFLAGS = -g
GNU_ANSI_FLAGS = -Wall -Wno-long-long -Wno-long-double -ansi -pedantic
SGI_ANSI_FLAGS = -ansi -fullwarn
ANSI_FLAGS = ${GNU_ANSI_FLAGS}

# Variables determined with AUTOCONFIG: 

INSTALL = install -p -c
RM = /bin/rm
CP = /bin/cp
MV = /bin/mv
CHMOD = /bin/chmod
MKDIR = /bin/mkdir
RANLIB = ranlib
TCLSH = ${FSLDIR}/bin/fsltclsh
