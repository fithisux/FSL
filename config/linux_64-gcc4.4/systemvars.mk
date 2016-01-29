# $Id: systemvars.mk,v 1.4 2015/09/07 10:55:55 cowboy Exp $

# System dependent paths

RM = /bin/rm
CHMOD = /bin/chmod
MKDIR = /bin/mkdir
CP = /bin/cp
MV = /bin/mv
INSTALL = install -p
TCLSH = ${FSLDIR}/bin/fsltclsh
RANLIB = echo

FSLML = ${FSLDIR}/bin/fslml

# for SHELL, do not change the type of shell - only use Bourne or BASH
SHELL = /bin/sh

# Compiler dependent variables

CC = gcc
CXX = c++
CSTATICFLAGS = -static
CXXSTATICFLAGS = -static

ARCHFLAGS = -m64 

PARALLELFLAGS = -fopenmp

DEPENDFLAGS = -MM

OPTFLAGS = -g -O3 -fexpensive-optimizations ${ARCHFLAGS}
MACHDBGFLAGS = -g
GNU_ANSI_FLAGS = -Wall -ansi -pedantic -Wno-long-long
SGI_ANSI_FLAGS = -ansi -fullwarn
ANSI_FLAGS = ${GNU_ANSI_FLAGS}

# CUDA development environment
CUDA_INSTALLATION = /opt/cuda-6.0
LIB_CUDA = ${CUDA_INSTALLATION}/lib64
INC_CUDA = ${CUDA_INSTALLATION}/inc
NVCC = ${CUDA_INSTALLATION}/bin/nvcc