/*  PVM_single_c.cu

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

#include "diffmodels_utils.h"
#include "levenberg_marquardt.cu"
#include "options.h"

#include <fstream>

/////////////////////////////////////
/////////////////////////////////////
/// 	    PVM_single_c	  /// 
/////////////////////////////////////
/////////////////////////////////////

__device__ 
inline float isoterm_PVM_single_c(const int pt,const float* _d,const float *bvals){
  	return exp(-bvals[pt]**_d);
}

__device__ 
inline float isoterm_lambda_PVM_single_c(const int pt,const float lambda,const float *bvals){
  	return(-2*bvals[pt]*lambda*exp(-bvals[pt]*lambda*lambda));
}

__device__ 
inline float anisoterm_PVM_single_c(const int pt,const float* _d,const float3 x, const float *bvecs, const float *bvals, const int ndirections){
	float dp = bvecs[pt]*x.x+bvecs[ndirections+pt]*x.y+bvecs[(2*ndirections)+pt]*x.z;
	return exp(-bvals[pt]**_d*dp*dp);
}

__device__ 
inline float anisoterm_lambda_PVM_single_c(const int pt,const float lambda,const float3 x, const float *bvecs, const float *bvals, const int ndirections){
	float dp = bvecs[pt]*x.x+bvecs[ndirections+pt]*x.y+bvecs[(2*ndirections)+pt]*x.z;
	return(-2*bvals[pt]*lambda*dp*dp*exp(-bvals[pt]*lambda*lambda*dp*dp));
}

__device__ 
inline float anisoterm_th_PVM_single_c(const int pt,const float* _d,const float3 x, const float _th,const float _ph,const float *bvecs, const float *bvals, const int ndirections){
	float sinth,costh,sinph,cosph;
	sincos(_th,&sinth,&costh);
	sincos(_ph,&sinph,&cosph);
	float dp = bvecs[pt]*x.x+bvecs[ndirections+pt]*x.y+bvecs[(2*ndirections)+pt]*x.z;
	float dp1 = costh*(bvecs[pt]*cosph+bvecs[ndirections+pt]*sinph)-bvecs[(2*ndirections)+pt]*sinth;
  	return(-2*bvals[pt]**_d*dp*dp1*exp(-bvals[pt]**_d*dp*dp));
}

__device__ 
inline float anisoterm_ph_PVM_single_c(const int pt,const float* _d,const float3 x, const float _th,const float _ph,const float *bvecs, const float *bvals, const int ndirections){
	float sinth,sinph,cosph;
	sinth=sin(_th);
	sincos(_ph,&sinph,&cosph);
  	float dp = bvecs[pt]*x.x+bvecs[ndirections+pt]*x.y+bvecs[(2*ndirections)+pt]*x.z;
	float dp1 = sinth*(-bvecs[pt]*sinph+bvecs[ndirections+pt]*cosph);
  	return(-2*bvals[pt]**_d*dp*dp1*exp(-bvals[pt]**_d*dp*dp));
}

//If the sum of the fractions is >1, then zero as many fractions
//as necessary, so that the sum becomes smaller than 1.
//in diffmodel.cc
__device__ void fix_fsum_PVM_single_c(		//INPUT 
						int nfib,
						//INPUT - OUTPUT){
						float *fs)
{
  	float sumf=0.0;
  	for(int i=0;i<nfib;i++){
    		sumf+=fs[i];
    		if(sumf>=1){
      			for(int j=i;j<nfib;j++) 
				fs[j]=FSMALL_gpu;  //make the fraction almost zero
      			break;
    		}
  	}
}

//in diffmodel.cc
__device__ void sort_PVM_single_c(int nfib,float* params)
{
	float temp_f, temp_th, temp_ph;
	// Order vector descending using f parameters as index
  	for(int i=1; i<(nfib); i++){ 
    		for(int j=0; j<(nfib-i); j++){ 
      			if (params[2+j*3] < params[2+(j+1)*3]){ 
        			temp_f = params[2+j*3];
				temp_th = params[2+j*3+1];
				temp_ph = params[2+j*3+2];
        			params[2+j*3] = params[2+(j+1)*3]; 
				params[2+j*3+1] = params[2+(j+1)*3+1]; 
				params[2+j*3+2] = params[2+(j+1)*3+2]; 
        			params[2+(j+1)*3] = temp_f; 
				params[2+(j+1)*3+1] = temp_th; 
				params[2+(j+1)*3+2] = temp_ph; 
      			} 
    		} 
  	} 
}

__device__  void fractions_deriv_PVM_single_c(	//INPUT
						const float*	params,
						const float* 	fs, 
						const int	nfib,
						const int	idSubVOX,
						//OUTPUT
						float* 		Deriv) 
{
	int nparams_per_fibre=3;
  	float fsum;
	int k=idSubVOX%nfib;
	for (int j=0; j<nfib; j++){
		Deriv[j*nfib+k]=0;
    	}

  	int kk = 2+(k*nparams_per_fibre);
	float sinparamkk = sin(2*params[kk]);

	for (int j=0; j<nfib; j++){
		int jj = 2+(j*nparams_per_fibre);
      		if (j==k){
			fsum=1; 
			for (int n=0; n<=(j-1); n++){
	  			fsum-=fs[n];
			}
			Deriv[j*nfib+k]=sinparamkk*fsum;
      		}else if (j>k){
			float sinparam = sin(params[jj]);
			fsum=0;
			for (int n=0; n<=(j-1); n++){
	  			fsum+=Deriv[n*nfib+k];
			}
			Deriv[j*nfib+k]=  -(sinparam*sinparam)*fsum;
      		}
    	}
}

//cost function PVM_single_c
__device__ void cf_PVM_single_c(	//INPUT
					const float*		params,
					const float*		mdata,
					const float*		bvecs, 
					const float*		bvals,
					const int 		ndirections,
					const int		nfib,
					const int 		nparams,
					const bool 		m_include_f0,
					const int		idSubVOX,
					float*			reduction,	//shared memory
					float* 			fs,		//shared memory
					float*			x,		//shared memory	
					float* 			_d,		//shared memory
					float* 			sumf,		//shared memory
					//OUTPUT
					double*			cfv)
{
	if(idSubVOX<nfib){
		int kk = 2+3*(idSubVOX);
		float sinth,costh,sinph,cosph;
		sincos(params[kk+1],&sinth,&costh);
		sincos(params[kk+2],&sinph,&cosph);
		x[idSubVOX*3] = sinth*cosph;
    		x[idSubVOX*3+1] = sinth*sinph;
    		x[idSubVOX*3+2] = costh;
  	}
	if(idSubVOX==0){
		*_d = lambda2d_gpu(params[1]);
		*cfv = 0.0;
		*sumf=0;
		float partial_fsum;
		for(int k=0;k<nfib;k++){
			int kk = 2+3*(k);
    			//partial_fsum ///////////
			partial_fsum=1.0;
			for(int j=0;j<k;j++)
				partial_fsum-=fs[j];
    			//////////////////////////
			fs[k] = beta2f_gpu(params[kk])*partial_fsum;
    			*sumf += fs[k];
		}
	}

	int ndir = ndirections/THREADS_BLOCK_FIT;
	if(idSubVOX<(ndirections%THREADS_BLOCK_FIT)) ndir++;
	
	float err;
	float3 x2;
	int dir_iter=idSubVOX;

	__syncthreads();
	
	reduction[idSubVOX]=0;
	for(int dir=0;dir<ndir;dir++){
		err = 0.0;
    		for(int k=0;k<nfib;k++){
			x2.x=x[k*3];
			x2.y=x[k*3+1];
			x2.z=x[k*3+2];	
			err += fs[k]*anisoterm_PVM_single_c(dir_iter,_d,x2,bvecs,bvals,ndirections); 
    		}
		if(m_include_f0){
			//partial_fsum ///////////
			float partial_fsum=1.0;
			for(int j=0;j<nfib;j++)
				partial_fsum-=fs[j];
	     		//////////////////////////
			float temp_f0=beta2f_gpu(params[nparams-1])*partial_fsum;
			err= (params[0]*((temp_f0+(1-*sumf-temp_f0)*isoterm_PVM_single_c(dir_iter,_d,bvals))+err))-mdata[dir_iter];
		}else{
			err = params[0]*((1-*sumf)*isoterm_PVM_single_c(dir_iter,_d,bvals)+err)-mdata[dir_iter];
		}
		reduction[idSubVOX]+= err*err;  
		dir_iter+=THREADS_BLOCK_FIT;
  	}  

	__syncthreads();

	if(idSubVOX==0){
		for(int i=0;i<THREADS_BLOCK_FIT;i++){
			*cfv+=reduction[i];
		}	
	}	
}


//gradient function PVM_single_c
__device__ void grad_PVM_single_c(	//INPUT
					const float*		params,
					const float*		mdata,
					const float*		bvecs, 
					const float*		bvals,
					const int		ndirections,
					const int		nfib,
					const int 		nparams,
					const bool 		m_include_f0,
					const int		idSubVOX,	
					float*			J,		//shared memory	
					float*			reduction,	//shared memory
					float* 			fs,		//shared memory
					float*			f_deriv,	//shared memory
					float*			x,		//shared memory
					float* 			_d,		//shared memory
					float* 			sumf,		//shared memory
					//OUTPUT
					float*			grad)
{
	if(idSubVOX<nfib){
		int kk = 2+3*(idSubVOX);
		float sinth,costh,sinph,cosph;
		sincos(params[kk+1],&sinth,&costh);
		sincos(params[kk+2],&sinph,&cosph);
    		x[idSubVOX*3] = sinth*cosph;
    		x[idSubVOX*3+1] = sinth*sinph;
    		x[idSubVOX*3+2] = costh;
  	}
	if(idSubVOX==0){
		*_d = lambda2d_gpu(params[1]);
		*sumf=0;
		float partial_fsum;
		for(int k=0;k<nfib;k++){
			int kk = 2+3*(k);
    			//partial_fsum ///////////
			partial_fsum=1.0;
			for(int j=0;j<k;j++)
				partial_fsum-=fs[j];
    			//////////////////////////
			fs[k] = beta2f_gpu(params[kk])*partial_fsum;
    			*sumf += fs[k];
		}
		for (int p=0;p<nparams;p++) grad[p]=0;
	}

	__syncthreads();

  	if(idSubVOX<nfib){ 
		fractions_deriv_PVM_single_c(params,fs,nfib,idSubVOX,f_deriv); 
	} 

	int ndir = ndirections/THREADS_BLOCK_FIT;
	if(idSubVOX<(ndirections%THREADS_BLOCK_FIT)) ndir++;
	int max_dir = ndirections/THREADS_BLOCK_FIT;
	if(ndirections%THREADS_BLOCK_FIT) max_dir++;

	float* myJ = &J[idSubVOX*nparams];
	float diff;
  	float sig;
	float Iso_term;
	float3 xx;
	int dir_iter=idSubVOX;
  	//float Aniso_terms[MAXNFIBRES];  //reuse Shared J --- myJ[kk+1]

	__syncthreads();

  	for(int dir=0;dir<max_dir;dir++){
		for (int p=0; p<nparams; p++) myJ[p]=0;
		if(dir<ndir){
    			for(int k=0;k<nfib;k++){
				int kk = 2+3*(k) +1;
      				xx.x=x[k*3];
      				xx.y=x[k*3+1];
     				xx.z=x[k*3+2];	
      				//Aniso_terms[k]=anisoterm_PVM_single_c(dir_iter,_d,xx,bvecs,bvals,ndirections);
				myJ[kk] = anisoterm_PVM_single_c(dir_iter,_d,xx,bvecs,bvals,ndirections);
    			}
			Iso_term=isoterm_PVM_single_c(dir_iter,_d,bvals);  //Precompute some terms for this datapoint
    			sig = 0;
    			for(int k=0;k<nfib;k++){
     				int kk = 2+3*(k);
      				xx.x=x[k*3];
      				xx.y=x[k*3+1];
      				xx.z=x[k*3+2];		
      				sig += fs[k]*myJ[kk+1];//Aniso_terms[k];
     				myJ[1] += params[0]*fs[k]*anisoterm_lambda_PVM_single_c(dir_iter,params[1],xx,bvecs,bvals,ndirections);
     				myJ[kk] = 0;
      				for (int j=0;j<nfib;j++){
					if(f_deriv[j*nfib+k]!=0){
	  					//myJ[kk] += params[0]*(Aniso_terms[j]-Iso_term)*f_deriv[j*nfib+k]; 
						myJ[kk] += params[0]*(myJ[2+j*3+1]-Iso_term)*f_deriv[j*nfib+k]; 
					}
      				}
			}
			for(int k=0;k<nfib;k++){
     				int kk = 2+3*(k);
      				xx.x=x[k*3];
      				xx.y=x[k*3+1];
      				xx.z=x[k*3+2];		
      				myJ[kk+1] = params[0]*fs[k]*anisoterm_th_PVM_single_c(dir_iter,_d,xx,params[kk+1],params[kk+2],bvecs,bvals,ndirections);
      				myJ[kk+2] = params[0]*fs[k]*anisoterm_ph_PVM_single_c(dir_iter,_d,xx,params[kk+1],params[kk+2],bvecs,bvals,ndirections);
    			}
    			if(m_include_f0){
				//partial_fsum ///////////
    				float partial_fsum=1.0;
    				for(int j=0;j<(nfib);j++)
					partial_fsum-=fs[j];
				//////////////////////////
				float temp_f0=beta2f_gpu(params[nparams-1])*partial_fsum;

    				//derivative with respect to f0
    				myJ[nparams-1]= params[0]*(1-Iso_term)*sin(float(2*params[nparams-1]))*partial_fsum; 
				sig=params[0]*((temp_f0+(1-*sumf-temp_f0)*Iso_term)+sig);
    				myJ[1] += params[0]*(1-*sumf-temp_f0)*isoterm_lambda_PVM_single_c(dir_iter,params[1],bvals);
    			}else{
				sig = params[0]*((1-*sumf)*Iso_term+sig);
	    			myJ[1] += params[0]*(1-*sumf)*isoterm_lambda_PVM_single_c(dir_iter,params[1],bvals);
    			}
    			diff = sig - mdata[dir_iter];
    			myJ[0] = sig/params[0]; 
		}

		for (int p=0;p<nparams;p++){ 
			reduction[idSubVOX]=2*myJ[p]*diff;

			__syncthreads();
			if(idSubVOX==0){
				for(int i=0;i<THREADS_BLOCK_FIT;i++){
					grad[p] += reduction[i];
				}
			}
			__syncthreads(); 
		} 
		dir_iter+=THREADS_BLOCK_FIT;
  	}
}


//hessian function PVM_single_c
__device__ void hess_PVM_single_c(	//INPUT
					const float*		params,
					const float*		bvecs, 
					const float*		bvals,
					const int 		ndirections,
					const int		nfib,
					const int 		nparams,
					const bool 		m_include_f0,
					const int		idSubVOX,		
					float*			J,		//shared memory
					float*			reduction,	//shared memory
					float* 			fs,		//shared memory
					float*			f_deriv,	//shared memory
					float*			x,		//shared memory
					float* 			_d,		//shared memory
					float* 			sumf,		//shared memory
					//OUTPUT
					float*			hess)
{
	if(idSubVOX<nfib){
		int kk = 2+3*(idSubVOX);
		float sinth,costh,sinph,cosph;
		sincos(params[kk+1],&sinth,&costh);
		sincos(params[kk+2],&sinph,&cosph);
    		x[idSubVOX*3] = sinth*cosph;
    		x[idSubVOX*3+1] = sinth*sinph;
    		x[idSubVOX*3+2] = costh;
  	}
	if(idSubVOX==0){
		*_d = lambda2d_gpu(params[1]);
		*sumf=0;
		float partial_fsum;
		for(int k=0;k<nfib;k++){
			int kk = 2+3*(k);
    			//partial_fsum ///////////
			partial_fsum=1.0;
			for(int j=0;j<k;j++)
				partial_fsum-=fs[j];
    			//////////////////////////
			fs[k] = beta2f_gpu(params[kk])*partial_fsum;
    			*sumf += fs[k];
		}
		for (int p=0;p<nparams;p++){
			for (int p2=0;p2<nparams;p2++){ 
				hess[p*nparams+p2] = 0;
			}
		}
	}

	__syncthreads();

  	if(idSubVOX<nfib){ 
		fractions_deriv_PVM_single_c(params,fs,nfib,idSubVOX,f_deriv); 
	} 

  	int ndir = ndirections/THREADS_BLOCK_FIT;
	if(idSubVOX<(ndirections%THREADS_BLOCK_FIT)) ndir++;
	int max_dir = ndirections/THREADS_BLOCK_FIT;
	if(ndirections%THREADS_BLOCK_FIT) max_dir++;

	float* myJ = &J[idSubVOX*nparams];
  	float sig;
	float Iso_term;
	float3 xx;
	int dir_iter=idSubVOX;
  	//float Aniso_terms[MAXNFIBRES]; //reuse Shared J --- myJ[kk+1]

	__syncthreads();

  	for(int dir=0;dir<max_dir;dir++){
		for (int p=0; p<nparams; p++) myJ[p]=0;
		if(dir<ndir){	
    			for(int k=0;k<nfib;k++){
				int kk = 2+3*(k) +1;
      				xx.x=x[k*3];
      				xx.y=x[k*3+1];
      				xx.z=x[k*3+2];	
      				//Aniso_terms[k]=anisoterm_PVM_single_c(dir_iter,_d,xx,bvecs,bvals,ndirections);
				myJ[kk] = anisoterm_PVM_single_c(dir_iter,_d,xx,bvecs,bvals,ndirections);
    			}
			Iso_term=isoterm_PVM_single_c(dir_iter,_d,bvals);  //Precompute some terms for this datapoint
    			sig = 0;
    			for(int k=0;k<nfib;k++){
      				int kk = 2+3*(k);
      				xx.x=x[k*3];
      				xx.y=x[k*3+1];
      				xx.z=x[k*3+2];		 
      				sig += fs[k]*myJ[kk+1];//Aniso_terms[k];
      				myJ[1] += params[0]*fs[k]*anisoterm_lambda_PVM_single_c(dir_iter,params[1],xx,bvecs,bvals,ndirections);	 
      				for (int j=0; j<nfib; j++){
					if (f_deriv[j*nfib+k]!=0)
	  				//myJ[kk] += params[0]*(Aniso_terms[j]-Iso_term)*f_deriv[j*nfib+k]; 
					myJ[kk] += params[0]*(myJ[2+3*j+1]-Iso_term)*f_deriv[j*nfib+k]; 
      				}
			}
			for(int k=0;k<nfib;k++){
				int kk = 2+3*(k);
      				xx.x=x[k*3];
      				xx.y=x[k*3+1];
      				xx.z=x[k*3+2];
      				myJ[kk+1] = params[0]*fs[k]*anisoterm_th_PVM_single_c(dir_iter,_d,xx,params[kk+1],params[kk+2],bvecs,bvals,ndirections);
      				myJ[kk+2] = params[0]*fs[k]*anisoterm_ph_PVM_single_c(dir_iter,_d,xx,params[kk+1],params[kk+2],bvecs,bvals,ndirections);
    			}
    			if(m_include_f0){
				//partial_fsum ///////////
	    			float partial_fsum=1.0;
	    			for(int j=0;j<(nfib);j++)
					partial_fsum-=fs[j];
	    			//////////////////////////
    				float temp_f0=beta2f_gpu(params[nparams-1])*partial_fsum;
    				//derivative with respect to f0
    				myJ[nparams-1]= params[0]*(1-Iso_term)*sin(float(2*params[nparams-1]))*partial_fsum; 
				sig= params[0]*((temp_f0+(1-*sumf-temp_f0)*Iso_term)+sig);
    				myJ[1] += params[0]*(1-*sumf-temp_f0)*isoterm_lambda_PVM_single_c(dir_iter,params[1],bvals);
    			}else{
	    			sig = params[0]*((1-*sumf)*Iso_term+sig);
	    			myJ[1] += params[0]*(1-*sumf)*isoterm_lambda_PVM_single_c(dir_iter,params[1],bvals);
    			}
    			myJ[0] = sig/params[0]; 
		}

		for (int p=0;p<nparams;p++){
			for (int p2=p;p2<nparams;p2++){ 

				reduction[idSubVOX]=2*(myJ[p]*myJ[p2]);
				__syncthreads();
				if(idSubVOX==0){
					for(int i=0;i<THREADS_BLOCK_FIT;i++){
						hess[p*nparams+p2] += reduction[i];
					}
				}
				__syncthreads(); 
			}
		}
		dir_iter+=THREADS_BLOCK_FIT;
  	}

	if(idSubVOX==0){
  		for (int j=0; j<nparams; j++) {
    			for (int i=j+1; i<nparams; i++) {
     				hess[i*nparams+j]=hess[j*nparams+i];	
    			}
  		}
	}
}
//in diffmodel.cc
extern "C" __global__ void fit_PVM_single_c_kernel(	//INPUT
							const float* 		data, 
							const float* 		bvecs, 
							const float* 		bvals, 
							const int 		nvox,
							const int		ndirections, 
							const int 		nfib,
							const int		nparams,
							const bool		m_eval_BIC,
							const bool 		m_include_f0,
							const bool	 	m_return_fanning,
							const bool		gradnonlin,
							//INPUT - OUTPUT
							float* 			params)
{
	int idSubVOX = threadIdx.x;
	int idVOX = blockIdx.x;
	int threadsBlock = blockDim.x;

	////////// DYNAMIC SHARED MEMORY ///////////
	extern __shared__ double shared[];
	double* pcf = (double*) shared;					//1   
	double* ncf = (double*) &pcf[1];				//1   
	double* lambda = (double*) &ncf[1];				//1  
	double* cftol = (double*) &lambda[1];				//1  
	double* ltol = (double*) &cftol[1];				//1  
	double* olambda = (double*) &ltol[1];				//1  

	float* J = (float*)&olambda[1];					//threadsBlock*nparams
	float* reduction = (float*)&J[threadsBlock*nparams];		//threadsBlock
	float* myparams = (float*) &reduction[threadsBlock];		//nparams
	float* grad = (float*) &myparams[nparams];			//nparams      
   	float* hess = (float*) &grad[nparams];				//nparams*nparams   
	float* step = (float*) &hess[nparams*nparams];			//nparams      
	float* inverse = (float*) &step[nparams];			//nparams   

	float* fs = (float*) &inverse[nparams];				//nfib
	float* f_deriv = (float*) &fs[nfib];				//nfib*nfib
  	float* x = (float*) &f_deriv[nfib*nfib];			//nfib*3
	float* _d = (float*) &x[nfib*3];				//1
  	float* sumf = (float*) &_d[1];					//1

	float* C = (float*)&sumf[1];					//nparams*nparams;
	float* el =  (float*)&C[nparams*nparams];			//nparams

	int* indx = (int*)&el[nparams];					//nparams
	int* success = (int*) &indx[nparams];				//1
	int* end = (int*) &success[1];					//1   
	////////// DYNAMIC SHARED MEMORY ///////////

	if(idSubVOX<nparams){
		myparams[idSubVOX]=params[(idVOX*nparams)+idSubVOX];
	}

	__syncthreads();

	int pos_bvals, pos_bvecs;
	if(gradnonlin){ 
		pos_bvals=idVOX*ndirections;
		pos_bvecs=idVOX*3*ndirections;
	}else{ 
		pos_bvals=0;
		pos_bvecs=0;
	}
	//do the fit
	levenberg_marquardt_PVM_single_c_gpu(&data[idVOX*ndirections],&bvecs[pos_bvecs],&bvals[pos_bvals],ndirections,nfib,nparams,m_include_f0,idSubVOX,step,grad,hess,inverse, pcf,ncf,lambda,cftol,ltol,olambda,success,end,J,reduction,fs,f_deriv,x,_d,sumf,C,el,indx,myparams);

	__syncthreads();

	// finalise parameters
	// m_s0-myparams[0] 	m_d-myparams[1] 	m_f-m_th-m_ph-myparams[2,3,4,5, etc..]   	m_f0-myparams[nparams-1]
	
	if(idSubVOX==0){
  		myparams[1] = lambda2d_gpu(myparams[1]); 
  		for(int k=0;k<nfib;k++){
    			int kk = 2 + 3*(k);
    			//partial_fsum ///////////
			float partial_fsum=1.0;
			for(int j=0;j<k;j++)
				partial_fsum-=myparams[2 + 3*j];
    			//////////////////////////
    			myparams[kk]  = beta2f_gpu(myparams[kk])*partial_fsum;
  		}
  
  		if (m_include_f0){
			//partial_fsum ///////////
	    		float partial_fsum=1.0;
	    		for(int j=0;j<(nfib);j++){
				partial_fsum-=myparams[2 + 3*j];
			}
	    		//////////////////////////
    			myparams[nparams-1]= beta2f_gpu(myparams[nparams-1])*partial_fsum;
		}
		sort_PVM_single_c(nfib,myparams);
	}
	__syncthreads();

	if(idSubVOX<nparams){
		params[(idVOX*nparams)+idSubVOX] = myparams[idSubVOX];
	}
}

//in diffmodel.cc
extern "C" __global__ void get_residuals_PVM_single_c_kernel(	//INPUT
								const float* 		data, 
								const float* 		params,
								const float* 		bvecs, 
								const float* 		bvals, 
								const int 		nvox, 
								const int		ndirections,
								const int 		nfib, 
								const int		nparams,
								const bool 		m_include_f0,
								const bool		gradnonlin,
								const bool* 		includes_f0,								
								//OUTPUT
								float*			residuals)
{
	int idSubVOX = threadIdx.x;
	int idVOX = blockIdx.x;
	int threadsBlock = blockDim.x;

	////////// DYNAMIC SHARED MEMORY ///////////
	extern __shared__ double shared[];
	float* myparams = (float*) shared;			//nparams
	float* fs = (float*) &myparams[nparams];		//nfib
  	float* x = (float*) &fs[nfib];				//nfib*3
	float* _d = (float*) &x[nfib*3];			//1
  	float* sumf = (float*) &_d[1];				//1
	int* my_include_f0 = (int*) &sumf[1];			//1		
	////////// DYNAMIC SHARED MEMORY ///////////
	
	float val; 
	float predicted_signal;
	float mydata;
	

	if(idSubVOX==0){
		*my_include_f0 = includes_f0[idVOX];

		//m_s0-myparams[0]  m_d-myparams[1] m_f-m_th-m_ph-myparams[2,3,4,5 etc..]  m_f0-myparams[nparams-1]
  		
  		myparams[0]=params[(idVOX*nparams)+0];
		if(myparams[1]<0)  myparams[1] = 0;	//This can be due to numerical errors..sqrt
  		else myparams[1] = d2lambda_gpu(params[(idVOX*nparams)+1]);

		float partial_fsum;	
  		for(int k=0;k<nfib;k++){
    			int kk = 2+3*k;
			//partial_fsum ///////////
			partial_fsum=1.0;
			for(int j=0;j<k;j++)
				partial_fsum-=fs[j];
	     		//////////////////////////
			fs[k] = params[(idVOX*nparams)+kk];
			float tmpr=fs[k]/partial_fsum;
    			if (tmpr>1.0) tmpr=1; //This can be due to numerical errors
			if (tmpr<0.0) tmpr=0; //This can be due to numerical errors..sqrt
    			myparams[kk]   = f2beta_gpu(tmpr);
    			myparams[kk+1] = params[(idVOX*nparams)+kk+1];
    			myparams[kk+2] = params[(idVOX*nparams)+kk+2];
  		}
  		if (*my_include_f0){
			//partial_fsum ///////////
			partial_fsum=1.0;
			for(int j=0;j<nfib;j++)
				partial_fsum-=fs[j];
	     		//////////////////////////	
			float tmpr=params[(idVOX*nparams)+nparams-1]/partial_fsum;
    			if (tmpr>1.0) tmpr=1; //This can be due to numerical errors..asin
			if (tmpr<0.0) tmpr=0; //This can be due to numerical errors..sqrt
    			myparams[nparams-1]= f2beta_gpu(tmpr);	
		}
	}

	__syncthreads();

	if(idSubVOX<nfib){
		int kk = 2+3*idSubVOX;
		float sinth,costh,sinph,cosph;
		sincos(myparams[kk+1],&sinth,&costh);
		sincos(myparams[kk+2],&sinph,&cosph);
    		x[idSubVOX*3] = sinth*cosph;
    		x[idSubVOX*3+1] = sinth*sinph;
    		x[idSubVOX*3+2] = costh;
  	}

	if(idSubVOX==0){
		float partial_fsum;	
		*sumf=0;
		for(int k=0;k<nfib;k++){
    			int kk = 2+3*k;
			////// partial_fsum //////
			partial_fsum=1.0;
			for(int j=0;j<k;j++)
				partial_fsum-=fs[j];
    			//////////////////////////
	    		fs[k] = beta2f_gpu(myparams[kk])*partial_fsum;
	    		*sumf += fs[k];
		}
		*_d = lambda2d_gpu(myparams[1]);
	}

	int ndir = ndirections/threadsBlock;
	if(idSubVOX<(ndirections%threadsBlock)) ndir++;
	
	float3 x2;
	int dir_iter=idSubVOX; 

	__syncthreads();

	int pos_bvals, pos_bvecs;
	if(gradnonlin){ 
		pos_bvals=idVOX*ndirections;
		pos_bvecs=idVOX*3*ndirections;
	}else{ 
		pos_bvals=0;
		pos_bvecs=0;
	}

	for(int dir=0;dir<ndir;dir++){
		mydata = data[(idVOX*ndirections)+dir_iter];
		predicted_signal=0;	//pred = 0;
		val = 0.0;
    		for(int k=0;k<nfib;k++){
			x2.x=x[k*3];
			x2.y=x[k*3+1];
			x2.z=x[k*3+2];	 
      			val += fs[k]*anisoterm_PVM_single_c(dir_iter,_d,x2,&bvecs[pos_bvecs],&bvals[pos_bvals],ndirections);
    		}	
    		if (*my_include_f0){
			//partial_fsum ///////////
			float partial_fsum=1.0;
			for(int j=0;j<nfib;j++)
				partial_fsum-=fs[j];
	     		//////////////////////////
      			float temp_f0= beta2f_gpu(myparams[nparams-1])*partial_fsum;
      			predicted_signal = myparams[0]*(temp_f0+(1-*sumf-temp_f0)*isoterm_PVM_single_c(dir_iter,_d,&bvals[pos_bvals])+val);
    		}else{
      			predicted_signal = myparams[0]*((1-*sumf)*isoterm_PVM_single_c(dir_iter,_d,&bvals[pos_bvals])+val); 
		}

		//residuals=m_data-predicted_signal;
		residuals[idVOX*ndirections+dir_iter]= mydata - predicted_signal;

		dir_iter+=threadsBlock;
	}
}
