/*  levenberg_marquardt.cu

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

#ifndef __LEVENBERG
#define __LEVENBERG

#include "solver_mult_inverse.cu"
#include "diffmodels.cuh"
#include "options.h"

//CPU version in nonlin.h
__device__ const double EPS_gpu = 2.0e-16;       	//Losely based on NRinC 20.1

//CPU version in nonlin.cpp
__device__ inline bool zero_cf_diff_conv(double* cfo,double* cfn,double* cftol){
  	return(2.0*abs(*cfo-*cfn) <= *cftol*(abs(*cfo)+abs(*cfn)+EPS_gpu));
}

__device__ void levenberg_marquardt_PVM_single_gpu(	//INPUT
							const float*		mydata, 
							const float*		bvecs, 
							const float*		bvals, 
							const int		ndirections,
							const int		nfib,
							const int 		nparams,
							const bool 		m_include_f0,
							const int		idSubVOX,
							float* 			step,		//shared memory
							float*			grad,           //shared memory     	          
						   	float* 			hess,		//shared memory
							float* 			inverse,	//shared memory
							double* 		pcf,		//shared memory
							double* 		ncf,		//shared memory
							double* 		lambda,		//shared memory
							double* 		cftol,		//shared memory
							double* 		ltol,		//shared memory
							double* 		olambda,	//shared memory
							int* 			success,    	//shared memory
							int* 			end,    	//shared memory
							float*			J,		//shared memory
							float*			reduction,	//shared memory
							float* 			fs,		//shared memory
						  	float*			x,		//shared memory
							float* 			_d,		//shared memory
						  	float* 			sumf,		//shared memory
							float*			C,		//shared memory
							float*			el,		//shared memory
							int*			indx,		//shared memory
							//INPUT-OUTPUT
							float*			myparams)	//shared memory
{
	int niter=0; 
	int maxiter=200;

   	if(idSubVOX==0){
		*end=false;
   		*lambda=0.1;
   		*cftol=1.0e-8;
   		*ltol=1.0e20;                  
   		*success = true;               
   		*olambda = 0.0;              
   		*ncf=0;  
	}

   	cf_PVM_single(myparams,mydata,bvecs,bvals,ndirections,nfib,nparams,m_include_f0,idSubVOX,reduction,fs,x,_d,sumf,pcf);  
	__syncthreads();

   	while (!(*success&&niter++>=maxiter)){ 	//if success we don't increase niter (first condition is true)
						//function cost has been decreased, we have advanced.
   		if(*success){
    			grad_PVM_single(myparams,mydata,bvecs,bvals,ndirections,nfib,nparams,m_include_f0,idSubVOX,J,reduction,fs,x,_d,sumf,grad); 
			__syncthreads(); 
    			hess_PVM_single(myparams,bvecs,bvals,ndirections,nfib,nparams,m_include_f0,idSubVOX,J,reduction,fs,x,_d,sumf,hess);  
    		}

		if(idSubVOX==0){
    			for (int i=0; i<nparams; i++) {                         
				hess[(i*nparams)+i]+=*lambda-*olambda;	//Levenberg LM_L
    			}

    			solver(hess,grad,nparams,C,el,indx,inverse);

    			for (int i=0;i<nparams;i++){
				step[i]=-inverse[i];		
    			}

   			for(int i=0;i<nparams;i++){
				step[i]=myparams[i]+step[i];
   			}
		}
		
		__syncthreads();
   		cf_PVM_single(step,mydata,bvecs,bvals,ndirections,nfib,nparams,m_include_f0,idSubVOX,reduction,fs,x,_d,sumf,ncf); 

		if(idSubVOX==0){
   			if (*success = (*ncf < *pcf)){ 
				*olambda = 0.0;
        			for(int i=0;i<nparams;i++){
					myparams[i]=step[i];
   				}
        			*lambda=*lambda/10.0;

				if (zero_cf_diff_conv(pcf,ncf,cftol)){
					*end=true;
				}
				*pcf=*ncf;
    			}else{
				*olambda=*lambda;
				*lambda=*lambda*10.0;
				if(*lambda> *ltol){ 
					*end=true;
				}
			}
    		}	
		__syncthreads();
		if(*end) return;		
   	}
}

__device__ void levenberg_marquardt_PVM_single_c_gpu(	//INPUT
							const float*		mydata, 
							const float*		bvecs, 
							const float*		bvals,
							const int		ndirections, 
							const int		nfib,
							const int 		nparams,
							const bool 		m_include_f0,
							const int		idSubVOX,
							float* 			step,		//shared memory
							float*			grad,           //shared memory     	          
						   	float* 			hess,		//shared memory
							float* 			inverse,	//shared memory
							double* 		pcf,		//shared memory
							double* 		ncf,		//shared memory
							double* 		lambda,		//shared memory
							double* 		cftol,		//shared memory
							double* 		ltol,		//shared memory
							double* 		olambda,	//shared memory
							int* 			success,    	//shared memory
							int* 			end,    	//shared memory
							float*			J,		//shared memory
							float*			reduction,	//shared memory
							float* 			fs,		//shared memory
							float*			f_deriv,	//shared memory
						  	float*			x,		//shared memory
							float* 			_d,		//shared memory
						  	float* 			sumf,		//shared memory
							float*			C,		//shared memory
							float*			el,		//shared memory
							int*			indx,		//shared memory
							//INPUT-OUTPUT
							float*			myparams)	//shared memory
{
	int niter=0; 
	int maxiter=200;

   	if(idSubVOX==0){
		*end=false;
   		*lambda=0.1;
   		*cftol=1.0e-8;
   		*ltol=1.0e20;                  
   		*success = true;               
   		*olambda = 0.0;              
   		*ncf=0;  
	}
			
	cf_PVM_single_c(myparams,mydata,bvecs,bvals,ndirections,nfib,nparams,m_include_f0,idSubVOX,reduction,fs,x,_d,sumf,pcf);  
	__syncthreads();
	
   	while (!(*success&&niter++ >= maxiter)){ 	//if success we don't increase niter (first condition is true)
							//function cost has been decreased, we have advanced.
   		if(*success){
			grad_PVM_single_c(myparams,mydata,bvecs,bvals,ndirections,nfib,nparams,m_include_f0,idSubVOX,J,reduction,fs,f_deriv,x,_d,sumf,grad);  
			__syncthreads();
    			hess_PVM_single_c(myparams,bvecs,bvals,ndirections,nfib,nparams,m_include_f0,idSubVOX,J,reduction,fs,f_deriv,x,_d,sumf,hess);  
    		}

		if(idSubVOX==0){
    			for (int i=0; i<nparams; i++) {                         
				hess[(i*nparams)+i]+=*lambda-*olambda;	//Levenberg LM_L
    			}

    			solver(hess,grad,nparams,C,el,indx,inverse);

    			for (int i=0;i<nparams;i++){
				step[i]=-inverse[i];		
    			}

   			for(int i=0;i<nparams;i++){
				step[i]=myparams[i]+step[i];
   			}
		}

		__syncthreads();
   		cf_PVM_single_c(step,mydata,bvecs,bvals,ndirections,nfib,nparams,m_include_f0,idSubVOX,reduction,fs,x,_d,sumf,ncf); 

		if(idSubVOX==0){
   			if (*success = (*ncf < *pcf)) {
				*olambda = 0.0;
        			for(int i=0;i<nparams;i++){
					myparams[i]=step[i];
   				}
        			*lambda=*lambda/10.0;

				if (zero_cf_diff_conv(pcf,ncf,cftol)){
					*end=true;
				}
				*pcf=*ncf;
    			}else{
				*olambda=*lambda;
				*lambda=*lambda*10.0;
				if(*lambda> *ltol){ 
					*end=true;
				}
    			}
		}
		__syncthreads();
		if(*end) return;		
   	}
}


__device__ void levenberg_marquardt_PVM_multi_gpu(	//INPUT
							const float*		mydata, 
							const float*		bvecs, 
							const float*		bvals, 
							const float		R,
							const float		invR,
							const int		ndirections,
							const int		nfib,
							const int 		nparams,
							const bool 		m_include_f0,
							const int		idSubVOX,
							const int		Gamma_for_ball_only,
							float* 			step,		//shared memory
							float*			grad,           //shared memory     	          
						   	float* 			hess,		//shared memory
							float* 			inverse,	//shared memory
							double* 		pcf,		//shared memory
							double* 		ncf,		//shared memory
							double* 		lambda,		//shared memory
							double* 		cftol,		//shared memory
							double* 		ltol,		//shared memory
							double* 		olambda,	//shared memory
							int* 			success,    	//shared memory
							int* 			end,    	//shared memory
							float*			J,		//shared memory
							float*			reduction,	//shared memory
							float* 			fs,		//shared memory
						  	float*			x,		//shared memory
							float* 			_a,		//shared memory
							float* 			_b,		//shared memory
						  	float* 			sumf,		//shared memory
							float*			C,		//shared memory
							float*			el,		//shared memory
							int*			indx,		//shared memory
							//INPUT-OUTPUT
							float*			myparams)	//shared memory
{
	int niter=0; 
	int maxiter=200;

   	if(idSubVOX==0){
		*end=false;
   		*lambda=0.1;
   		*cftol=1.0e-8;
   		*ltol=1.0e20;                  
   		*success = true;               
   		*olambda = 0.0;              
   		*ncf=0;  
	}

	cf_PVM_multi(myparams,mydata,bvecs,bvals,R,invR,ndirections,nfib,nparams,m_include_f0,idSubVOX,Gamma_for_ball_only,reduction,fs,x,_a,_b,sumf,pcf);  
	__syncthreads();
	
   	while (!(*success&&niter++ >= maxiter)){ 	//if success we don't increase niter (first condition is true)
							//function cost has been decreased, we have advanced.
   		if(*success){
			grad_PVM_multi(myparams,mydata,bvecs,bvals,R,invR,ndirections,nfib,nparams,m_include_f0,
			idSubVOX,Gamma_for_ball_only,J,reduction,fs,x,_a,_b,sumf,grad);  

			__syncthreads(); 
    			hess_PVM_multi(myparams,bvecs,bvals,R,invR,ndirections,nfib,nparams,m_include_f0,idSubVOX,Gamma_for_ball_only,J,reduction,fs,x,_a,_b,sumf,hess);  
    		}

		if(idSubVOX==0){
    			for (int i=0; i<nparams; i++) {                         
				hess[(i*nparams)+i]+=*lambda-*olambda;	//Levenberg LM_L
    			}

    			solver(hess,grad,nparams,C,el,indx,inverse);

    			for (int i=0;i<nparams;i++){
				step[i]=-inverse[i];		
    			}

   			for(int i=0;i<nparams;i++){
				step[i]=myparams[i]+step[i];
   			}
		}

		__syncthreads();
   		cf_PVM_multi(step,mydata,bvecs,bvals,R,invR,ndirections,nfib,nparams,m_include_f0,idSubVOX,Gamma_for_ball_only,reduction,fs,x,_a,_b,sumf,ncf); 

		if(idSubVOX==0){
   			if (*success = (*ncf < *pcf)) {
				*olambda = 0.0;
        			for(int i=0;i<nparams;i++){
					myparams[i]=step[i];
   				}
        			*lambda=*lambda/10.0;

				if (zero_cf_diff_conv(pcf,ncf,cftol)){
					*end=true;
				}
				*pcf=*ncf;
    			}else{
				*olambda=*lambda;
				*lambda=*lambda*10.0;
				if(*lambda> *ltol){ 
					*end=true;
				}
    			}
		}
		__syncthreads();
		if(*end) return;				
   	}
}
#endif
