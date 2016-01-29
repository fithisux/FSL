# for SHELL, do not change the type of shell - only use Bourne or BASH
SHELL = /bin/sh

# Compiler dependent variables

CC = gcc
CXX = g++
CSTATICFLAGS = 
CXXSTATICFLAGS = 

ARCHFLAGS = -arch x86_64
ARCHLDFLAGS = -Wl,-search_paths_first -arch x86_64
PER_ARCH_CFLAGS_x86_64 = -msse3
CFLAGS = -std=c99

DEPENDFLAGS = -MM

OPTFLAGS =  -O3
MACHDBGFLAGS = -g
GNU_ANSI_FLAGS = -Wall -Wno-long-long -ansi -pedantic
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
