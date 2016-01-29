
include ${FSLCONFDIR}/default.mk

PROJNAME = avwutils

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L${LIB_ZLIB}

LIBS = -lnewimage -lmiscmaths -lprob -lfslio -lnewmat -lutils -lniftiio -lznz -lm -lz

IOLIBS = -lfslio -lniftiio -lznz -lm -lz

XFILES = fslorient fslstats fslmerge \
         fslpspec fslroi fslnvols fsl2ascii fslascii2img  \
         fslsplit fslcc fslinterleave \
         fslhd fslcpgeom fslcreatehd fslcorrecthd fslmaths \
         fslcomplex fslfft fslmeants fslslice fslswapdim_exe fslchfiletype_exe calc_grad_perc_dev fslsmoothfill fslselectvols

SCRIPTS = fslval fslchpixdim fslanimate fslsize sliceanimate fslinfo fsledithd avw2fsl fslswapdim fslchfiletype fslreorient2std fslmodhd fsladd

FSCRIPTS = fmriborient 

all: ${XFILES}

fslstats: fslstats.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslstats.o ${LIBS} 

fslfft: fslfft.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o fslfft fslfft.cc ${LIBS}

fslmeants: fslmeants.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o fslmeants fslmeants.cc ${LIBS}

avwfixfloat: avwfixfloat.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o avwfixfloat avwfixfloat.cc ${LIBS}

fslcomplex: fslcomplex.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o fslcomplex fslcomplex.cc ${LIBS}

fslswapdim_exe: fslswapdim.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o fslswapdim_exe fslswapdim.o ${LIBS} 

fslorient: fslorient.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslorient.o ${LIBS} 

fslslice: fslslice.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o fslslice fslslice.cc ${LIBS}

fslpspec: fslpspec.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslpspec.o ${LIBS} 


#start of matts c++ rewrites
fslroi: fslroi.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslroi.o ${LIBS}

fslnvols: fslnvols.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslnvols.o ${LIBS}

fsl2ascii: fsl2ascii.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fsl2ascii.o ${LIBS}

fslascii2img: fslascii2img.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslascii2img.o ${LIBS}

fslsplit: fslsplit.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslsplit.o ${LIBS}

fslmerge: fslmerge.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslmerge.o ${LIBS}

fslcc: fslcc.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslcc.o ${LIBS}

fslinterleave: fslinterleave.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslinterleave.o ${LIBS}

fslhd: fslhd.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslhd.o ${LIBS}

fslcpgeom: fslcpgeom.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslcpgeom.o ${LIBS}

fslcorrecthd: fslcorrecthd.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslcorrecthd.o ${LIBS}

fslcreatehd: fslcreatehd.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslcreatehd.o ${LIBS}

fslmaths: fslmaths.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslmaths.o ${LIBS}

fsldecorr4d: fsldecorr4d.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fsldecorr4d.o ${LIBS}

fslchfiletype_exe: fslchfiletype.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o fslchfiletype_exe fslchfiletype.o ${LIBS}

calc_grad_perc_dev: calc_grad_perc_dev.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o calc_grad_perc_dev calc_grad_perc_dev.o ${LIBS}

fslsmoothfill: fslsmoothfill.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslsmoothfill.o ${LIBS}

fslselectvols: fslselectvols.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslselectvols.o ${LIBS}
