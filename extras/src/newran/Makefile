include ${FSLCONFDIR}/default.mk

PROJNAME = newran

OPTFLAGS = -O2
LIBS = -lnewran -lm

AFILES = libnewran.a
TESTXFILES = tryrand

all:   libnewran.a 
tests: all ${TESTXFILES}

%.o:          %.cpp
	      ${CP} $*.cpp $*.cxx  
	      $(CXX) $(CXXFLAGS) -c $*.cxx
	      ${RM} -f $*.cxx  


OBJ_LIB = hist.o newran.o extreal.o myexcept.o
libnewran.a:	$(OBJ_LIB)
	        $(AR) cr $@ $(OBJ_LIB)
	        ${RANLIB} $@

OBJ_T = tryrand.o tryrand1.o tryrand2.o tryrand3.o tryrand4.o

tryrand:    	$(OBJ_T) libnewran.a
	        $(CXX)  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJ_T} ${LIBS}




