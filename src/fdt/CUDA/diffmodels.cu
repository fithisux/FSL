/*  diffmodels.cu

    Tim Behrens, Saad Jbabdi, Stam Sotiropoulos, Moises Hernandez  - FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk
    
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
    
    LICENCE
    
    FMRIB Software Library, Release 5.0 (c) 2012, The University of
    Oxford (the "Software")
    
    The Software remains the property of the University of Oxford ("the
    University").
    
    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.
    
    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.
    
    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.
    
    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. If you are
    interested in using the Software commercially, please contact Isis
    Innovation Limited ("Isis"), the technology transfer company of the
    University, to negotiate a licence. Contact details are:
    innovation@isis.ox.ac.uk quoting reference DE/9564. */

#include "newmat.h"
#include "newimage/newimageall.h"
#include "diffmodels.h"

#include "PVM_single.cu"
#include "PVM_single_c.cu"
#include "PVM_multi.cu"
#include "fit_gpu_kernels.h"
#include "sync_check.h"

#include <time.h>
#include <sys/time.h>
#include "init_gpu.h"

////////////////////////////////////////////////////// 
//   FIT IN GPU
////////////////////////////////////////////////////// 
void fit_PVM_single(	//INPUT
			const vector<ColumnVector> 	datam_vec, 
			const vector<Matrix> 		bvecs_vec,
			const vector<Matrix> 		bvals_vec,
			thrust::device_vector<float> 	datam_gpu, 
			thrust::device_vector<float>	bvecs_gpu, 
			thrust::device_vector<float>	bvals_gpu,
			int				ndirections,	
			int 				nfib,
			bool 				m_include_f0,
			bool				gradnonlin,
			string 				output_file,		
			//OUTPUT
			thrust::device_vector<float>&	params_gpu)
{
	std::ofstream myfile;
	myfile.open (output_file.data(), ios::out | ios::app );

	int nvox = datam_vec.size();
	int nparams;
	if (m_include_f0)
      		nparams = nfib*3 + 3; 
    	else
      		nparams = nfib*3 + 2;

	thrust::host_vector<float> params_host;
	params_host.resize(nvox*nparams);
	
	for(int vox=0;vox<nvox;vox++){
		// initialise with a tensor
		int pos_bv;
		if (gradnonlin) pos_bv=vox;
		else pos_bv=0;
  		DTI dti(datam_vec[vox],bvecs_vec[pos_bv],bvals_vec[pos_bv]);
  		dti.linfit();

  		// set starting parameters for nonlinear fitting
  		float _th,_ph;
  		cart2sph(dti.get_v1(),_th,_ph);

  		params_host[vox*nparams] = dti.get_s0();
		params_host[vox*nparams+1] = dti.get_md()>0?dti.get_md()*2:0.001; // empirically found that d~2*MD
  		//params_host[vox*nparams+1] = dti.get_l1()>0?dti.get_l1():0.002; // empirically found that d~L1
	  	params_host[vox*nparams+2] = dti.get_fa()<1?f2x(dti.get_fa()):f2x(0.95); // first pvf = FA 
	  	params_host[vox*nparams+3] = _th;
	  	params_host[vox*nparams+4] = _ph;
	  	float sumf=x2f(params_host[vox*nparams+2]);
	  	float tmpsumf=sumf;
	  	for(int ii=2,i=5;ii<=nfib;ii++,i+=3){
		    	float denom=2;
		    	do{
		      		params_host[vox*nparams+i] = f2x(x2f(params_host[vox*nparams+i-3])/denom);
		      		denom *= 2;
		      		tmpsumf = sumf + x2f(params_host[vox*nparams+i]);
		    	}while(tmpsumf>=1);
		    	sumf += x2f(params_host[vox*nparams+i]);
		    	cart2sph(dti.get_v(ii),_th,_ph);
		    	params_host[vox*nparams+i+1] = _th;
		    	params_host[vox*nparams+i+2] = _ph;
	  	}
	  	if (m_include_f0)
	    		params_host[vox*nparams+nparams-1]=f2x(FSMALL);
	}

	thrust::copy(params_host.begin(), params_host.end(), params_gpu.begin());	

	int blocks = nvox;
   	dim3 Dim_Grid(blocks,1);
  	dim3 Dim_Block(THREADS_BLOCK_FIT,1);

	int amount_shared = 6*sizeof(double)+(THREADS_BLOCK_FIT*nparams+THREADS_BLOCK_FIT+5*nparams+2*nparams*nparams+4*nfib+2)*sizeof(float)+(2+nparams)*sizeof(int);

	myfile << "Shared Memory Used in fit_PVM_single: " << amount_shared << "\n"; 

	fit_PVM_single_kernel<<<Dim_Grid, Dim_Block, amount_shared>>>(thrust::raw_pointer_cast(datam_gpu.data()), thrust::raw_pointer_cast(bvecs_gpu.data()), thrust::raw_pointer_cast(bvals_gpu.data()) ,nvox, ndirections, nfib, nparams, m_include_f0, gradnonlin,thrust::raw_pointer_cast(params_gpu.data()));
	sync_check("fit_PVM_single_kernel");

	myfile.close();
}

void fit_PVM_single_c(	//INPUT
			const vector<ColumnVector> 	datam_vec, 
			const vector<Matrix> 		bvecs_vec,
			const vector<Matrix> 		bvals_vec,
			thrust::device_vector<float> 	datam_gpu, 
			thrust::device_vector<float>	bvecs_gpu, 
			thrust::device_vector<float>	bvals_gpu,
			int				ndirections,
			int 				nfib,		
			bool 				m_include_f0,
			bool				gradnonlin,
			string 				output_file,	
			//OUTPUT
			thrust::device_vector<float>&	params_gpu)
{
	std::ofstream myfile;
	myfile.open (output_file.data(), ios::out | ios::app );

	int nvox = datam_vec.size(); 
	int nparams;
	if (m_include_f0)
      		nparams = nfib*3 + 3; 
    	else
      		nparams = nfib*3 + 2;

	thrust::host_vector<float> params_host;
	params_host.resize(nvox*nparams);

	for(int vox=0;vox<nvox;vox++){
		// initialise with a tensor
		int pos_bv;
		if (gradnonlin) pos_bv=vox;
		else pos_bv=0;
  		DTI dti(datam_vec[vox],bvecs_vec[pos_bv],bvals_vec[pos_bv]);
  		dti.linfit();

  		// set starting parameters for nonlinear fitting
  		float _th,_ph;
  		cart2sph(dti.get_v1(),_th,_ph);

  		ColumnVector start(nparams);
  		//Initialize the non-linear fitter. Use the DTI estimates for most parameters, apart from the volume fractions
  		start(1) = dti.get_s0();
  		//start(2) = d2lambda(dti.get_md()>0?dti.get_md()*2:0.001); // empirically found that d~2*MD
  		start(2) = d2lambda(dti.get_l1()>0?dti.get_l1():0.002); // empirically found that d~L1
  		start(4) = _th;
  		start(5) = _ph;
  		for(int ii=2,i=6;ii<=nfib;ii++,i+=3){
    			cart2sph(dti.get_v(ii),_th,_ph);
    			start(i+1) = _th;
    			start(i+2) = _ph;
  		}
  
  		// do a better job for initializing the volume fractions
		PVM_single_c pvm(datam_vec[vox],bvecs_vec[pos_bv],bvals_vec[pos_bv],nfib,false,m_include_f0,false);
  		pvm.fit_pvf(start);

		for(int i=0;i<nparams;i++){ 
			params_host[vox*nparams+i]=start(i+1);
		}
	}

	thrust::copy(params_host.begin(), params_host.end(), params_gpu.begin());	

	int blocks = nvox;
   	dim3 Dim_Grid(blocks,1);
  	dim3 Dim_Block(THREADS_BLOCK_FIT,1);

	int amount_shared = 6*sizeof(double)+(THREADS_BLOCK_FIT*nparams+THREADS_BLOCK_FIT+5*nparams+2*nparams*nparams+4*nfib+nfib*nfib+6)*sizeof(float)+(2+nparams)*sizeof(int);

	myfile << "Shared Memory Used in fit_PVM_single_c: " << amount_shared << "\n"; 

	fit_PVM_single_c_kernel<<<Dim_Grid, Dim_Block, amount_shared>>>(thrust::raw_pointer_cast(datam_gpu.data()), thrust::raw_pointer_cast(bvecs_gpu.data()), thrust::raw_pointer_cast(bvals_gpu.data()) ,nvox, ndirections, nfib, nparams, false, m_include_f0, false, gradnonlin,thrust::raw_pointer_cast(params_gpu.data()));
	sync_check("fit_PVM_single_c_kernel");

	myfile.close();
}

void fit_PVM_multi(	//INPUT
			thrust::device_vector<float> 	datam_gpu, 
			thrust::device_vector<float>	bvecs_gpu, 
			thrust::device_vector<float>	bvals_gpu,	
			int 				nvox,	
			int				ndirections,
			int				nfib,		
			bool 				m_include_f0,
			bool				gradnonlin,
			float				R_prior_mean,
			int				Gamma_ball_only,
			string 				output_file,
			//OUTPUT
			thrust::device_vector<float>&	params_gpu)
{
	std::ofstream myfile;
	myfile.open (output_file.data(), ios::out | ios::app );

	int blocks = nvox;
   	dim3 Dim_Grid(blocks,1);
  	dim3 Dim_Block(THREADS_BLOCK_FIT,1);

	int nparams;
	if (m_include_f0)
      		nparams = nfib*3 + 4; 
    	else
      		nparams = nfib*3 + 3;

	thrust::device_vector<float> params_PVM_single_c_gpu; 	//copy params to an auxiliar structure because there are different number of nparams
	params_PVM_single_c_gpu.resize(nvox*nparams);		//between single_c and multi. We must read and write in different structures, 
	thrust::copy(params_gpu.begin(), params_gpu.end(), params_PVM_single_c_gpu.begin());	
								//maybe 1 block finish before other one read their params.

	int amount_shared = 6*sizeof(double)+(THREADS_BLOCK_FIT*nparams+THREADS_BLOCK_FIT+5*nparams+2*nparams*nparams+4*nfib+3)*sizeof(float)+(2+nparams)*sizeof(int);

	myfile << "Shared Memory Used in fit_PVM_multi: " << amount_shared << "\n"; 

	float invR =  1.0/(2.0*R_prior_mean+1);

	fit_PVM_multi_kernel<<<Dim_Grid, Dim_Block, amount_shared>>>(thrust::raw_pointer_cast(datam_gpu.data()), thrust::raw_pointer_cast(params_PVM_single_c_gpu.data()), thrust::raw_pointer_cast(bvecs_gpu.data()), thrust::raw_pointer_cast(bvals_gpu.data()), R_prior_mean, invR, nvox, ndirections, nfib, nparams, Gamma_ball_only, m_include_f0, gradnonlin,thrust::raw_pointer_cast(params_gpu.data()));
	sync_check("fit_PVM_multi_kernel");

	myfile.close();
}

void calculate_tau(	//INPUT
			thrust::device_vector<float> 	datam_gpu, 
			thrust::device_vector<float>	params_gpu,
			thrust::device_vector<float>	bvecs_gpu, 
			thrust::device_vector<float>	bvals_gpu,
			thrust::host_vector<int>	vox_repeat,
			int				nrepeat,
			int				ndirections,
			int				nfib,
			int 				model,
			bool 				m_include_f0,
			bool 				nonlin,
			bool				gradnonlin,
			float				R_prior_mean,
			int				Gamma_ball_only,
			string 				output_file,	
			//OUTPUT
			thrust::host_vector<float>&	tau)
{
	std::ofstream myfile;
	myfile.open (output_file.data(), ios::out | ios::app );
   	myfile << "--------- CALCULATE TAU/RESIDULAS IN GPU ------------ " << "\n"; 

	struct timeval t1,t2;
   	double time;
   	gettimeofday(&t1,NULL);

	int nvox = vox_repeat.size(); 

	int nparams;
	if (m_include_f0)
      		nparams = nfib*3 + 3; 
    	else
      		nparams = nfib*3 + 2;
	if(model>=2) nparams++;
	
	thrust::device_vector<bool> includes_f0_gpu;
	includes_f0_gpu.resize(nvox);
	thrust::fill(includes_f0_gpu.begin(), includes_f0_gpu.end(), m_include_f0);

	if(m_include_f0){
		for(int i=0;i<nrepeat;i++){
			includes_f0_gpu[vox_repeat[i]]=false;	//if has been reprocessed f0 will be 0.
		}
	}

	int blocks = nvox;
   	dim3 Dim_Grid(blocks,1);
  	dim3 Dim_Block(THREADS_BLOCK_FIT,1);

	int amount_shared = (nparams+4*nfib+3)*sizeof(float) + sizeof(int);

	thrust::device_vector<float> residuals_gpu;
	residuals_gpu.resize(nvox*ndirections);

	if(model==1){
		if(nonlin){ 
			get_residuals_PVM_single_kernel<<<Dim_Grid, Dim_Block,amount_shared>>>(thrust::raw_pointer_cast(datam_gpu.data()), thrust::raw_pointer_cast(params_gpu.data()), thrust::raw_pointer_cast(bvecs_gpu.data()), thrust::raw_pointer_cast(bvals_gpu.data()), nvox, ndirections, nfib, nparams, m_include_f0,gradnonlin, thrust::raw_pointer_cast(includes_f0_gpu.data()), thrust::raw_pointer_cast(residuals_gpu.data()));
			sync_check("get_residuals_PVM_single_kernel");

		}else{
			get_residuals_PVM_single_c_kernel<<<Dim_Grid, Dim_Block,amount_shared>>>(thrust::raw_pointer_cast(datam_gpu.data()), thrust::raw_pointer_cast(params_gpu.data()), thrust::raw_pointer_cast(bvecs_gpu.data()), thrust::raw_pointer_cast(bvals_gpu.data()), nvox, ndirections, nfib, nparams, m_include_f0, gradnonlin,thrust::raw_pointer_cast(includes_f0_gpu.data()), thrust::raw_pointer_cast(residuals_gpu.data()));
			sync_check("get_residuals_PVM_single_c_kernel");
		}
	}else{
		float invR =  1.0/(2.0*R_prior_mean+1);
      		//model 2 or 3 
		get_residuals_PVM_multi_kernel<<<Dim_Grid, Dim_Block,amount_shared>>>(thrust::raw_pointer_cast(datam_gpu.data()), thrust::raw_pointer_cast(params_gpu.data()), thrust::raw_pointer_cast(bvecs_gpu.data()), thrust::raw_pointer_cast(bvals_gpu.data()), R_prior_mean, invR,nvox, ndirections, nfib, nparams, Gamma_ball_only, m_include_f0, gradnonlin,thrust::raw_pointer_cast(includes_f0_gpu.data()), thrust::raw_pointer_cast(residuals_gpu.data()));
		sync_check("get_residuals_PVM_multi_kernel");
	}

	thrust::host_vector<float> residuals_host;
	residuals_host.resize(nvox*ndirections);
	thrust::copy(residuals_gpu.begin(), residuals_gpu.end(), residuals_host.begin());	

	ColumnVector res(ndirections);
	for(int vox=0;vox<nvox;vox++){	
		for(int i=0;i<ndirections;i++) res(i+1)= residuals_host[vox*ndirections+i];  

		float variance=var(res).AsScalar();
		tau[vox]=1.0/variance;
	}

	gettimeofday(&t2,NULL);
    	time=timeval_diff(&t2,&t1);
   	myfile << "TIME TOTAL: " << time << " seconds\n"; 
	myfile << "-----------------------------------------------------" << "\n\n" ; 	
	myfile.close();			
}

