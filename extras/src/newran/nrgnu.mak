CXX = g++
CXXFLAGS = -O -Dbool_LIB -Wall

everything:   tryrand

%.o:          %.cpp
	      $(CXX) $(CXXFLAGS) -c $*.cpp

OBJ_TR = tryrand.o tryrand1.o tryrand2.o tryrand3.o tryrand4.o hist.o \
      newran.o extreal.o myexcept.o

tryrand:      $(OBJ_TR)
	      $(CXX) -o $@ $(OBJ_TR) -L. -lm

newranxx = include.h newran.h boolean.h myexcept.h

tryrand.o:    $(newranxx) tryrand.cpp 

tryrand1.o:   $(newranxx) tryrand1.cpp 

tryrand2.o:   $(newranxx) tryrand2.cpp 

tryrand3.o:   $(newranxx) tryrand3.cpp 

tryrand4.o:   $(newranxx) tryrand4.cpp 

hist.o:       $(newranxx) hist.cpp 

newran.o:     $(newranxx) newran.cpp 

extreal.o:    include.h boolean.h extreal.h myexcept.h extreal.cpp 

myexcept.o:   include.h boolean.h myexcept.h myexcept.cpp 

