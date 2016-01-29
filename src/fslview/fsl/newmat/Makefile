include ${FSLCONFDIR}/default.mk

PROJNAME = newmat

#OPTFLAGS = -O2
LIBS = -lnewmat -lm

AFILES = libnewmat.a
TESTXFILES = tmt example nl_ex sl_ex garch test_exc

all:   libnewmat.a 
tests: all ${TESTXFILES}

%.o:          %.cpp
	      ${CP} $*.cpp $*.cxx  
	      $(CXX) $(CXXFLAGS) -c $*.cxx
	      ${RM} -f $*.cxx  


OBJ_LIB = bandmat.o cholesky.o evalue.o fft.o hholder.o     \
  jacobi.o myexcept.o newmat1.o newmat2.o newmat3.o         \
  newmat4.o newmat5.o newmat6.o newmat7.o newmat8.o         \
  newmat9.o newmatex.o newmatnl.o newmatrm.o solution.o     \
  sort.o submat.o svd.o newfft.o

libnewmat.a:	$(OBJ_LIB)
	        $(AR) cr $@ $(OBJ_LIB)
	        ${RANLIB} $@

OBJ_T = tmt.o tmt1.o tmt2.o tmt3.o tmt4.o tmt5.o tmt6.o    \
  tmt7.o tmt8.o tmt9.o tmta.o tmtb.o tmtc.o tmtd.o tmte.o  \
  tmtf.o tmtg.o tmth.o tmti.o tmtj.o tmtk.o tmtl.o

tmt:    	$(OBJ_T) libnewmat.a
	        $(CXX)  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJ_T} ${LIBS}

OBJ_E = example.o

example:    	$(OBJ_E) libnewmat.a
	        $(CXX)  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJ_E} ${LIBS}

OBJ_N = nl_ex.o

nl_ex:    	$(OBJ_N) libnewmat.a
	        $(CXX)  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJ_N} ${LIBS}

OBJ_S = sl_ex.o

sl_ex:    	$(OBJ_S) libnewmat.a
	        $(CXX)  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJ_S} ${LIBS}

OBJ_G = garch.o

garch:          $(OBJ_G) libnewmat.a
	        $(CXX)  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJ_G} ${LIBS}

OBJ_X = test_exc.o

test_exc:	$(OBJ_X) libnewmat.a
	        $(CXX)  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJ_X} ${LIBS}






