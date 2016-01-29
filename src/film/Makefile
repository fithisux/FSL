include ${FSLCONFDIR}/default.mk

PROJNAME = film

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L${LIB_ZLIB}

LIBS = -lfslsurface -lfslvtkio -lmeshclass  -lfirst_lib -lnewimage -lmiscmaths -lutils -lm -lnewmat -lgiftiio -lfslio -lniftiio -lznz -lprob -lz -lexpat

XFILES = film_gls contrast_mgr ttoz ftoz ttologp

OBJS =  ContrastMgrOptions.o ContrastMgr.o gaussComparer.o AutoCorrEstimator.o glm.o paradigm.o FilmOlsOptions.o FilmGlsOptionsRes.o FilmGlsOptions.o glimGls.o

all:	${XFILES}

contrast_mgr:${OBJS} contrast_mgr.o
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} contrast_mgr.o ${LIBS}

ttoz:${OBJS} ttoz.o
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} ttoz.o ${LIBS}

ttologp:${OBJS} ttologp.o
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} ttologp.o ${LIBS}

ftoz:${OBJS} ftoz.o
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} ftoz.o ${LIBS}

film_gls:${OBJS} film_gls.o
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} film_gls.o ${LIBS} -l giftiio

film_gls_res:${OBJS} film_gls_res.o
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} film_gls_res.o ${LIBS}

quick:${OBJS} quick.o
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} quick.o ${LIBS}

testSmooth:testSmooth.o
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@ testSmooth.o ${LIBS} -lgiftiio
