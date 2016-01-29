# Makefile for FLIRT

include ${FSLCONFDIR}/default.mk

PROJNAME = flirt

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB} -I${INC_BOOST}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_ZLIB}

LIBS = -lwarpfns -lbasisfield -lnewimage -lmiscmaths -lfslio -lniftiio -lznz -lnewmat -lutils -lm -lz

FL_OBJS = globaloptions.o flirt.o 
C_OBJS = convert_xfm.o
A_OBJS = avscale.o
R_OBJS = rmsdiff.o
STD_OBJS = std2imgcoord.o
EPI_OBJS = img2stdcoord.o
IMG_OBJS = img2imgcoord.o
X_OBJS = applyxfm4D.o
P_OBJS = pointflirt.o
M_OBJS = makerot.o
I_OBJS = imapper.o
MT_OBJS = midtrans.o

RUNTCLS = Flirt InvertXFM ApplyXFM ConcatXFM Nudge
XFILES = flirt convert_xfm avscale rmsdiff std2imgcoord img2stdcoord \
	img2imgcoord applyxfm4D pointflirt makerot midtrans
TESTXFILES = 
HFILES =
SCRIPTS = extracttxt pairreg standard_space_roi flirt_average epi_reg aff2rigid

all:	${XFILES} schedule

schedule:
	@if [ ! -d ${DESTDIR}/etc ] ; then ${MKDIR} ${DESTDIR}/etc ; ${CHMOD} g+w ${DESTDIR}/etc ; fi
	@if [ ! -d ${DESTDIR}/etc/flirtsch ] ; then ${MKDIR} ${DESTDIR}/etc/flirtsch ; ${CHMOD} g+w ${DESTDIR}/etc/flirtsch ; fi
	${CP} -rf flirtsch/* ${DESTDIR}/etc/flirtsch/.

flirt:    	${FL_OBJS}
	        $(CXX)  ${CXXFLAGS} ${LDFLAGS} -o $@ ${FL_OBJS} ${LIBS}


convert_xfm:    ${C_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${C_OBJS} ${LIBS}

avscale:        ${A_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${A_OBJS} ${LIBS}

rmsdiff:        ${R_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${R_OBJS} ${LIBS}

std2imgcoord:	${STD_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${STD_OBJS} ${LIBS}

img2stdcoord:	${EPI_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${EPI_OBJS} ${LIBS}

img2imgcoord:	${IMG_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${IMG_OBJS} ${LIBS}

applyxfm4D:	${X_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${X_OBJS} ${LIBS}

pointflirt:	${P_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${P_OBJS} ${LIBS}

makerot:	${M_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${M_OBJS} ${LIBS}

imapper:	${I_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${I_OBJS} ${LIBS}

midtrans:	${MT_OBJS}
		${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${MT_OBJS} ${LIBS}

