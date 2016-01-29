include $(FSLCONFDIR)/default.mk

ifeq ($(COMPILE_GPU), 1)
	COMPILE_WITH_GPU=libbedpostx_cuda.so merge_parts_gpu xfibres_gpu CUDA/split_parts_gpu
	SCRIPTS_GPU=CUDA/bedpostx_gpu CUDA/bedpostx_postproc_gpu.sh
endif

PROJNAME = fdt

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_NEWRAN} -I${INC_CPROB} -I${INC_PROB} -I${INC_BOOST} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_NEWRAN} -L${LIB_CPROB} -L${LIB_PROB} -L${LIB_ZLIB}

DLIBS = -lwarpfns -lbasisfield -lmeshclass -lnewimage -lutils -lmiscmaths -lnewmat -lnewran -lfslio -lniftiio -lznz -lcprob -lprob -lm -lz

DTIFIT=dtifit
CCOPS=ccops
PTX=probtrackx
MED=medianfilter
ROM=reord_OM
SAUS=sausages
XFIBRES=xfibres
XFIBRES2=xfibres_2
RV=replacevols
MDV=make_dyadic_vectors
FMO=fdt_matrix_ops
INDEXER=indexer
TEST=testfile
VECREG=vecreg
KURTOSIS=kurtosis
SWAPDYADS=swap_dyadic_vectors
PVMFIT=pvmfit
DTIGEN=dtigen
BASGEN=basgen
RARNG=rearrange
XPRED=xfibres_pred
RUBIX=rubix
RBXPRED=rubix_pred
EDDYCOMBINE=eddy_combine
LIBBEDPOSTX_CUDA=libbedpostx_cuda.so
MERGE_PARTS_GPU=merge_parts_gpu
SPLIT_PARTS_GPU=CUDA/split_parts_gpu
XFIBRES_GPU=xfibres_gpu

DTIFITOBJS=dtifit.o dtifitOptions.o diffmodels.o Bingham_Watson_approx.o
CCOPSOBJS=ccops.o ccopsOptions.o
PTXOBJS=probtrackx.o probtrackxOptions.o streamlines.o ptx_simple.o ptx_seedmask.o ptx_twomasks.o ptx_nmasks.o ptx_meshmask.o
MEDOBJS=medianfilter.o 
ROMOBJS=reord_OM.o
SAUSOBJS=sausages.o
XFIBOBJS=xfibres.o xfibresoptions.o diffmodels.o Bingham_Watson_approx.o 
XFIBOBJS2=xfibres_2.o xfibresoptions.o
RVOBJS=replacevols.o
MDVOBJS=make_dyadic_vectors.o
FMOOBJS=fdt_matrix_ops.o
INDEXEROBJS=indexer.o
TESTOBJS=testfile.o 
VECREGOBJS=vecreg.o
KURTOSISOBJS=kurtosis.o dtifitOptions.o
SWAPDYADSOBJS=swap_dyadic_vectors.o
PVMFITOBJS=pvmfit.o pvmfitOptions.o diffmodels.o Bingham_Watson_approx.o
DTIGENOBJS=dtigen.o
BASGENOBJS=basgen.o diffmodels.o Bingham_Watson_approx.o
RARNGOBJS=rearrange.o
XPREDOBJS=xfibres_pred.o
RUBIXOBJS=rubix.o diffmodels.o rubixvox.o rubixoptions.o Bingham_Watson_approx.o
RBXPREDOBJS=rubix_pred.o
EDDYCOMBINEOBJS=eddy_combine.o
MERGE_PARTS_GPUOBJS=merge_parts_gpu.o xfibresoptions.o
SPLIT_PARTS_GPUOBJS=CUDA/split_parts_gpu.o
XFIBRES_GPUOBJS=xfibres_gpu.o xfibresoptions.o diffmodels.o Bingham_Watson_approx.o

SGEBEDPOST = bedpost 
SGEBEDPOSTX = bedpostx bedpostx_postproc.sh bedpostx_preproc.sh bedpostx_single_slice.sh bedpostx_datacheck

SCRIPTS = eddy_correct zeropad maskdyads probtrack fdt_rotate_bvecs select_dwi_vols ${SGEBEDPOST} ${SGEBEDPOSTX} ${SCRIPTS_GPU}
FSCRIPTS = correct_and_average ocmr_preproc

XFILES = dtifit ccops medianfilter make_dyadic_vectors vecreg xfibres probtrackx pvmfit dtigen eddy_combine ${COMPILE_WITH_GPU}

FXFILES = reord_OM sausages replacevols fdt_matrix_ops indexer rearrange xfibres_pred


RUNTCLS = Fdt

all: ${XFILES} ${FXFILES} 

${PTX}:		   ${PTXOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${PTXOBJS} ${DLIBS}

${PT}:		   ${PTOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${PTOBJS} ${DLIBS} 

${FTB}:    	${FTBOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${FTBOBJS} ${DLIBS} 

${PJ}:    	${PJOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${PJOBJS} ${DLIBS} 

${MED}:    	${MEDOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${MEDOBJS} ${DLIBS} 

${DTIFIT}:    	${DTIFITOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${DTIFITOBJS} ${DLIBS}

${CCOPS}:    	${CCOPSOBJS}	
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${CCOPSOBJS} ${DLIBS}

${ROM}:    	${ROMOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${ROMOBJS} ${DLIBS}

${SAUS}:    	${SAUSOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${SAUSOBJS} ${DLIBS}

${XFIBRES}:    	${XFIBOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${XFIBOBJS} ${DLIBS}

${XFIBRES2}:    	${XFIBOBJS2}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${XFIBOBJS2} ${DLIBS}

${RV}:    	${RVOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${RVOBJS} ${DLIBS}

${MDV}:    	${MDVOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${MDVOBJS} ${DLIBS}

${FMO}:    	${FMOOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${FMOOBJS} ${DLIBS}

${INDEXER}:    	${INDEXEROBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${INDEXEROBJS} ${DLIBS}

${TEST}:    	${TESTOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${TESTOBJS} ${DLIBS}

${VECREG}:    	${VECREGOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${VECREGOBJS} ${DLIBS}


${KURTOSIS}:   ${KURTOSISOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${KURTOSISOBJS} ${DLIBS}

${SWAPDYADS}: ${SWAPDYADSOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${SWAPDYADSOBJS} ${DLIBS}

${PVMFIT}:    	${PVMFITOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${PVMFITOBJS} ${DLIBS}

${DTIGEN}:    	${DTIGENOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${DTIGENOBJS} ${DLIBS}

${BASGEN}:    	${BASGENOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${BASGENOBJS} ${DLIBS}

${RARNG}: 	${RARNGOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${RARNGOBJS} ${DLIBS}

${XPRED}: 	${XPREDOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${XPREDOBJS} ${DLIBS}

${RUBIX}: 	${RUBIXOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${RUBIXOBJS} ${DLIBS}

${RBXPRED}: 	${RBXPREDOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${RBXPREDOBJS} ${DLIBS}

${EDDYCOMBINE}: ${EDDYCOMBINEOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${EDDYCOMBINEOBJS} ${DLIBS}

${LIBBEDPOSTX_CUDA}: 
		${CUDA}/bin/nvcc --shared --compiler-options '-fPIC' -o CUDA/libbedpostx_cuda.so CUDA/init_gpu.cu CUDA/samples.cu CUDA/diffmodels.cu CUDA/runmcmc.cu  CUDA/xfibres_gpu.cu -O3  -gencode=arch=compute_20,code=\"sm_20,compute_20\" -gencode=arch=compute_30,code=\"sm_30,compute_30\" -gencode=arch=compute_35,code=\"sm_35,compute_35\" -gencode=arch=compute_50,code=\"sm_50,compute_50\" -lcudart -lcuda -lcurand -I. -L${CUDA}/lib64 -L${CUDA}/lib -ICUDA/options -I${CUDA}/include/thrust -I${FSLDIR}/extras/include/newmat -I${FSLDIR}/include -I${FSLDIR}/extras/include/boost -maxrregcount=64
		@if [ ! -d ${FSLDEVDIR}/lib/ ] ; then ${MKDIR} ${FSLDEVDIR}/lib ; fi
		${CP} -rf CUDA/libbedpostx_cuda.so ${FSLDEVDIR}/lib

${MERGE_PARTS_GPU}: ${MERGE_PARTS_GPUOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${MERGE_PARTS_GPUOBJS} ${DLIBS}

${SPLIT_PARTS_GPU}: ${SPLIT_PARTS_GPUOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${SPLIT_PARTS_GPUOBJS} ${DLIBS}

${XFIBRES_GPU}: ${XFIBRES_GPUOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${XFIBRES_GPUOBJS} ${DLIBS} -lcudart -lcuda -lcurand -lbedpostx_cuda -LCUDA -L${CUDA}/lib64 -L${CUDA}/lib
