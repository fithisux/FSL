/*  PVM_multi.cu

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

/////////////////////////////////////
/////////////////////////////////////
/// 	    PVM_multi	 	  /// 
/////////////////////////////////////
/////////////////////////////////////

__device__ inline float isoterm_PVM_multi(const int pt,const float* _a,const float* _b, const float *bvals){
	return exp(-*_a*log(1+bvals[pt]**_b));
}

__device__ inline float isoterm_a_PVM_multi(const int pt,const float* _a,const float* _b, const float *bvals){
    	return  -log(1+bvals[pt]**_b)*exp(-*_a*log(1+bvals[pt]**_b));
}

__device__ inline float isoterm_b_PVM_multi(const int pt,const float* _a,const float* _b, const float *bvals){
      	return -*_a*bvals[pt]/(1+bvals[pt]**_b)*exp(-*_a*log(1+bvals[pt]**_b));
}

__device__ inline float anisoterm_PVM_multi(const int pt,const float* _a,const float* _b,const float3 x,const float *bvecs, const float *bvals, const float R, const float invR, const int ndirections,const int Gamma_for_ball_only){
	float dp = bvecs[pt]*x.x+bvecs[ndirections+pt]*x.y+bvecs[(2*ndirections)+pt]*x.z;
	if(Gamma_for_ball_only==1){
		return exp(-bvals[pt]**_a**_b*dp*dp);
	}else if(Gamma_for_ball_only==2){
		return exp(-bvals[pt]*3**_a**_b*invR*((1-R)*dp*dp+R));		
	}else{
  		return exp(-*_a*log(1+bvals[pt]**_b*(dp*dp)));
	}
}
 
__device__ inline float anisoterm_a_PVM_multi(const int pt,const float* _a,const float* _b,const float3 x,const float *bvecs, const float *bvals, const float R, const float invR, const int ndirections,const int Gamma_for_ball_only){
	float dp = bvecs[pt]*x.x+bvecs[ndirections+pt]*x.y+bvecs[(2*ndirections)+pt]*x.z;
	if(Gamma_for_ball_only==1){
		return (-bvals[pt]**_b*dp*dp* exp(-bvals[pt]**_a**_b*dp*dp));
  	}else if(Gamma_for_ball_only==2){
		float dp2=bvals[pt]*3**_b*invR*((1-R)*dp*dp+R);
		return(-dp2*exp(-dp2**_a));
	}else{
		return -log(1+bvals[pt]*(dp*dp)**_b)* exp(-*_a*log(1+bvals[pt]*(dp*dp)**_b));
  	}
  	
}

__device__ inline float anisoterm_b_PVM_multi(const int pt,const float* _a,const float* _b,const float3 x,const float *bvecs, const float *bvals, const float R, const float invR, const int ndirections,const int Gamma_for_ball_only){
  	float dp = bvecs[pt]*x.x+bvecs[ndirections+pt]*x.y+bvecs[(2*ndirections)+pt]*x.z;
	if(Gamma_for_ball_only==1){
		return(-bvals[pt]**_a*dp*dp*exp(-bvals[pt]**_a**_b*dp*dp));
  	}else if(Gamma_for_ball_only==2){
		float dp2=bvals[pt]*3**_a*invR*((1-R)*dp*dp+R);
		return(-dp2*exp(-dp2**_b));
  	}else{
		return (-*_a*bvals[pt]*(dp*dp)/ (1+bvals[pt]*(dp*dp)**_b)*exp(-*_a*log(1+bvals[pt]*(dp*dp)**_b)));
  	}
}

__device__ inline float anisoterm_th_PVM_multi(const int pt,const float* _a,const float* _b,const float3 x,const float _th,const float _ph,const float *bvecs, const float *bvals, const float R, const float invR, const int ndirections,const int Gamma_for_ball_only){
	float sinth,costh,sinph,cosph;
	sincos(_th,&sinth,&costh);
	sincos(_ph,&sinph,&cosph);
  	float dp = bvecs[pt]*x.x+bvecs[ndirections+pt]*x.y+bvecs[(2*ndirections)+pt]*x.z;
  	float dp1 = costh* (bvecs[pt]*cosph + bvecs[ndirections+pt]*sinph) - bvecs[(2*ndirections)+pt]*sinth;
	if(Gamma_for_ball_only==1){
  		return(-2*bvals[pt]**_a**_b*dp*dp1*exp(-bvals[pt]**_a**_b*dp*dp));
  	}else if(Gamma_for_ball_only==2){
		float dp2=2*bvals[pt]*3**_a**_b*invR*(1-R)*dp1;
		return(-dp2*exp(-bvals[pt]*3**_a**_b*invR*((1-R)*dp*dp+R)));
  	}else{
		return  (-*_a**_b*bvals[pt]/(1+bvals[pt]*(dp*dp)**_b)*exp(-*_a*log(1+bvals[pt]*(dp*dp)**_b))*2*dp*dp1);	
	}
}

__device__ inline float anisoterm_ph_PVM_multi(const int pt,const float* _a,const float* _b,const float3 x,const float _th,const float _ph,const float *bvecs, const float *bvals, const float R, const float invR, const int ndirections,const int Gamma_for_ball_only){
	float sinth,sinph,cosph;
	sinth=sin(_th);
	sincos(_ph,&sinph,&cosph);
  	float dp = bvecs[pt]*x.x+bvecs[ndirections+pt]*x.y+bvecs[(2*ndirections)+pt]*x.z;
  	float dp1 = sinth* (-bvecs[pt]*sinph + bvecs[ndirections+pt]*cosph);
	if(Gamma_for_ball_only==1){
  		return(-2*bvals[pt]**_a**_b*dp*dp1*exp(-bvals[pt]**_a**_b*dp*dp));
  	}else if(Gamma_for_ball_only==2){
		float dp2=2*bvals[pt]*3**_a**_b*invR*(1-R)*dp1;
		return(-dp2*exp(-bvals[pt]*3**_a**_b*invR*((1-R)*dp*dp+R)));
 	}else{
		return  (-*_a**_b*bvals[pt]/(1+bvals[pt]*(dp*dp)**_b)*exp(-*_a*log(1+bvals[pt]*(dp*dp)**_b))*2*dp*dp1);
  	}
}

//in diffmodel.cc
__device__ void fix_fsum_PVM_multi(	//INPUT 
					bool m_include_f0, 
					int nfib,
					int nparams,
					//INPUT - OUTPUT){
					float *params)
{
  	float sumf=0;
  	if (m_include_f0) 
    		sumf=params[nparams-1];
  	for(int i=0;i<nfib;i++){
    		if (params[3+(i*3)]==0) 
			params[3+(i*3)]=FSMALL_gpu;
    		sumf+=params[3+(i*3)];
    		if(sumf>=1){
			for(int j=i;j<nfib;j++)
				params[3+(j*3)]=FSMALL_gpu;
			break;
		}
	}
}

//in diffmodel.cc
__device__ void sort_PVM_multi(int nfib,float* params)
{
	float temp_f, temp_th, temp_ph;
	// Order vector descending using f parameters as index
	for(int i=1; i<(nfib); i++){ 
    		for(int j=0; j<(nfib-i); j++){ 
      			if (params[3+j*3] < params[3+(j+1)*3]){ 
        			temp_f = params[3+j*3];
				temp_th = params[3+j*3+1];
				temp_ph = params[3+j*3+2];
        			params[3+j*3] = params[3+(j+1)*3]; 
				params[3+j*3+1] = params[3+(j+1)*3+1]; 
				params[3+j*3+2] = params[3+(j+1)*3+2]; 
        			params[3+(j+1)*3] = temp_f; 
				params[3+(j+1)*3+1] = temp_th; 
				params[3+(j+1)*3+2] = temp_ph; 
      			} 
    		} 
  	} 
}

//cost function PVM_multi
__device__ void cf_PVM_multi(		//INPUT
					const float*		params,
					const float*		mdata,
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
					float*			reduction,	//shared memory
					float* 			fs,		//shared memory
					float*			x,		//shared memory	
					float* 			_a,		//shared memory
					float* 			_b,		//shared memory
					float* 			sumf,		//shared memory
					//OUTPUT
					double*			cfv)
{
	if(idSubVOX<nfib){
		int kk = 3+3*(idSubVOX);
		float sinth,costh,sinph,cosph;
		sincos(params[kk+1],&sinth,&costh);
		sincos(params[kk+2],&sinph,&cosph);
		fs[idSubVOX] = x2f_gpu(params[kk]);
		x[idSubVOX*3] = sinth*cosph;
    		x[idSubVOX*3+1] = sinth*sinph;
    		x[idSubVOX*3+2] = costh;
  	}
	if(idSubVOX==0){
		*_a= abs(params[1]);
		*_b= abs(params[2]); 
		*cfv = 0.0;
		*sumf=0;
		for(int k=0;k<nfib;k++) *sumf+= fs[k];
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
			err += fs[k]*anisoterm_PVM_multi(dir_iter,_a,_b,x2,bvecs,bvals,R,invR,ndirections,Gamma_for_ball_only); 
    		}
		if(m_include_f0){
			float temp_f0=x2f_gpu(params[nparams-1]);
			err = (abs(params[0])*(temp_f0+((1-*sumf-temp_f0)*isoterm_PVM_multi(dir_iter,_a,_b,bvals)+err)))-mdata[dir_iter];
		}else{
			err = abs(params[0])*((1-*sumf)*isoterm_PVM_multi(dir_iter,_a,_b,bvals)+err)-mdata[dir_iter];
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

//gradient function PVM_multi
__device__ void grad_PVM_multi(		//INPUT
					const float*		params,
					const float*		mdata,
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
					float*			J,		//shared memory
					float*			reduction,	//shared memory
					float* 			fs,		//shared memory
					float*			x,		//shared memory	
					float* 			_a,		//shared memory
					float* 			_b,		//shared memory
					float* 			sumf,		//shared memory
					//OUTPUT
					float*			grad)
{	
	if(idSubVOX<nfib){
		int kk = 3+3*(idSubVOX);
		float sinth,costh,sinph,cosph;
		sincos(params[kk+1],&sinth,&costh);
		sincos(params[kk+2],&sinph,&cosph);
		fs[idSubVOX] = x2f_gpu(params[kk]);
		x[idSubVOX*3] = sinth*cosph;
    		x[idSubVOX*3+1] = sinth*sinph;
    		x[idSubVOX*3+2] = costh;
  	}
	if(idSubVOX==0){
		*_a= abs(params[1]);
		*_b= abs(params[2]); 
		*sumf=0;
		for(int k=0;k<nfib;k++) *sumf+= fs[k];
		for (int p=0;p<nparams;p++) grad[p]=0;
	}

  	int ndir = ndirections/THREADS_BLOCK_FIT;
	if(idSubVOX<(ndirections%THREADS_BLOCK_FIT)) ndir++;
	int max_dir = ndirections/THREADS_BLOCK_FIT;
	if(ndirections%THREADS_BLOCK_FIT) max_dir++;

	float* myJ = &J[idSubVOX*nparams];
	float diff;
  	float sig;
	float3 xx;
	int dir_iter=idSubVOX;

	__syncthreads();

  	for(int dir=0;dir<max_dir;dir++){
		for (int p=0; p<nparams; p++) myJ[p]=0;
		if(dir<ndir){
    			sig = 0;
    			for(int k=0;k<nfib;k++){
      				int kk = 3+3*(k);
      				xx.x=x[k*3];
      				xx.y=x[k*3+1];
      				xx.z=x[k*3+2];		
      				sig += fs[k]*anisoterm_PVM_multi(dir_iter,_a,_b,xx,bvecs,bvals,R,invR,ndirections,Gamma_for_ball_only);

      				myJ[1] += (params[1]>0?1.0:-1.0)*abs(params[0])*fs[k]*
				anisoterm_a_PVM_multi(dir_iter,_a,_b,xx,bvecs,bvals,R,invR,ndirections,Gamma_for_ball_only); 

				myJ[2] += (params[2]>0?1.0:-1.0)*abs(params[0])*fs[k]*
				anisoterm_b_PVM_multi(dir_iter,_a,_b,xx,bvecs,bvals,R,invR,ndirections,Gamma_for_ball_only);

				myJ[kk] = abs(params[0])*(anisoterm_PVM_multi(dir_iter,_a,_b,xx,bvecs,bvals,R,invR,ndirections,Gamma_for_ball_only)
				-isoterm_PVM_multi(dir_iter,_a,_b,bvals))*two_pi_gpu*sign_gpu(params[kk])*1/(1+params[kk]*params[kk]); 

      				myJ[kk+1] = abs(params[0])*fs[k]*
				anisoterm_th_PVM_multi(dir_iter,_a,_b,xx,params[kk+1],params[kk+2],bvecs,bvals,R,invR,ndirections,Gamma_for_ball_only);  

      				myJ[kk+2] = abs(params[0])*fs[k]*
				anisoterm_ph_PVM_multi(dir_iter,_a,_b,xx,params[kk+1],params[kk+2],bvecs,bvals,R,invR,ndirections,Gamma_for_ball_only);
    			}
    			if(m_include_f0){
				float temp_f0=x2f_gpu(params[nparams-1]);
				myJ[nparams-1]= abs(params[0])*(1-isoterm_PVM_multi(dir_iter,_a,_b,bvals))*
				two_pi_gpu*sign_gpu(params[nparams-1])*1/(1+params[nparams-1]*params[nparams-1]);

				sig=abs(params[0])*((temp_f0+(1-*sumf-temp_f0)*isoterm_PVM_multi(dir_iter,_a,_b,bvals))+sig);
    				myJ[1] += (params[1]>0?1.0:-1.0)*abs(params[0])*(1-*sumf-temp_f0)*isoterm_a_PVM_multi(dir_iter,_a,_b,bvals);
				myJ[2] += (params[2]>0?1.0:-1.0)*abs(params[0])*(1-*sumf-temp_f0)*isoterm_b_PVM_multi(dir_iter,_a,_b,bvals);
    			}else{
	    			sig = abs(params[0]) * ((1-*sumf)*isoterm_PVM_multi(dir_iter,_a,_b,bvals)+sig);
	    			myJ[1] += (params[1]>0?1.0:-1.0)*abs(params[0])*(1-*sumf)*isoterm_a_PVM_multi(dir_iter,_a,_b,bvals);
	    			myJ[2] += (params[2]>0?1.0:-1.0)*abs(params[0])*(1-*sumf)*isoterm_b_PVM_multi(dir_iter,_a,_b,bvals);	
    			}
    
    			diff = sig - mdata[dir_iter];
    			myJ[0] = (params[0]>0?1.0:-1.0)*sig/params[0]; 
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

//hessian function PVM_multi 
__device__ void hess_PVM_multi(		//INPUT
					const float*		params,
					const float*		bvecs, 
					const float*		bvals,
					const float		R,
					const float		invR,
					const int 		ndirections,
					const int		nfib,
					const int 		nparams,
					const bool 		m_include_f0,
					const int		idSubVOX,
					const int		Gamma_for_ball_only,
					float*			J,		//shared memory
					float*			reduction,	//shared memory
					float* 			fs,		//shared memory
					float*			x,		//shared memory	
					float* 			_a,		//shared memory
					float* 			_b,		//shared memory
					float* 			sumf,		//shared memory
					//OUTPUT
					float*			hess)
{
	if(idSubVOX<nfib){
		int kk = 3+3*(idSubVOX);
		float sinth,costh,sinph,cosph;
		sincos(params[kk+1],&sinth,&costh);
		sincos(params[kk+2],&sinph,&cosph);
		fs[idSubVOX] = x2f_gpu(params[kk]);
		x[idSubVOX*3] = sinth*cosph;
    		x[idSubVOX*3+1] = sinth*sinph;
    		x[idSubVOX*3+2] = costh;
  	}
	if(idSubVOX==0){
		*_a= abs(params[1]);
		*_b= abs(params[2]); 
		*sumf=0;
		for(int k=0;k<nfib;k++) *sumf+= fs[k];
		for (int p=0;p<nparams;p++){
			for (int p2=0;p2<nparams;p2++){ 
				hess[p*nparams+p2] = 0;
			}
		}
	}

  	int ndir = ndirections/THREADS_BLOCK_FIT;
	if(idSubVOX<(ndirections%THREADS_BLOCK_FIT)) ndir++;
	int max_dir = ndirections/THREADS_BLOCK_FIT;
	if(ndirections%THREADS_BLOCK_FIT) max_dir++;

	float* myJ = &J[idSubVOX*nparams];
  	float sig;
	float3 xx;
	int dir_iter=idSubVOX; 

	__syncthreads(); 
	
  	for(int dir=0;dir<max_dir;dir++){
		for (int p=0; p<nparams; p++) myJ[p]=0;
		if(dir<ndir){
    			sig = 0;
    			for(int k=0;k<nfib;k++){
      				int kk = 3+3*(k);
      				xx.x=x[k*3];
      				xx.y=x[k*3+1];
      				xx.z=x[k*3+2];		
      				sig += fs[k]*anisoterm_PVM_multi(dir_iter,_a,_b,xx,bvecs,bvals,R,invR,ndirections,Gamma_for_ball_only);

      				float cov = two_pi_gpu*sign_gpu(params[kk])*1/(1+params[kk]*params[kk]);	
      				myJ[1] += (params[1]>0?1.0:-1.0)*abs(params[0])*fs[k]*
				anisoterm_a_PVM_multi(dir_iter,_a,_b,xx,bvecs,bvals,R,invR,ndirections,Gamma_for_ball_only);

				myJ[2] += (params[2]>0?1.0:-1.0)*abs(params[0])*fs[k]*
				anisoterm_b_PVM_multi(dir_iter,_a,_b,xx,bvecs,bvals,R,invR,ndirections,Gamma_for_ball_only);

				myJ[kk] = abs(params[0])*
				(anisoterm_PVM_multi(dir_iter,_a,_b,xx,bvecs,bvals,R,invR,ndirections,Gamma_for_ball_only)-
				isoterm_PVM_multi(dir_iter,_a,_b,bvals))*cov;

      				myJ[kk+1] = abs(params[0])*fs[k]*
				anisoterm_th_PVM_multi(dir_iter,_a,_b,xx,params[kk+1],params[kk+2],bvecs,bvals,R,invR,ndirections,Gamma_for_ball_only);

      				myJ[kk+2] = abs(params[0])*fs[k]*
				anisoterm_ph_PVM_multi(dir_iter,_a,_b,xx,params[kk+1],params[kk+2],bvecs,bvals,R,invR,ndirections,Gamma_for_ball_only);
    			}
    			if(m_include_f0){
				float temp_f0=x2f_gpu(params[nparams-1]);
				myJ[nparams-1]= abs(params[0])*(1-isoterm_PVM_multi(dir_iter,_a,_b,bvals))*two_pi_gpu*sign_gpu(params[nparams-1])*1/(1+params[nparams-1]*params[nparams-1]);
	    			sig = abs(params[0])* (temp_f0+(1-*sumf-temp_f0)*isoterm_PVM_multi(dir_iter,_a,_b,bvals)+sig);
    				myJ[1] += (params[1]>0?1.0:-1.0)*abs(params[0])*(1-*sumf-temp_f0)*isoterm_a_PVM_multi(dir_iter,_a,_b,bvals);
				myJ[2] += (params[2]>0?1.0:-1.0)*abs(params[0])*(1-*sumf-temp_f0)*isoterm_b_PVM_multi(dir_iter,_a,_b,bvals);
    			}else{
				sig = abs(params[0])*((1-*sumf)*isoterm_PVM_multi(dir_iter,_a,_b,bvals)+sig);
	    			myJ[1] += (params[1]>0?1.0:-1.0)*abs(params[0])*(1-*sumf)*isoterm_a_PVM_multi(dir_iter,_a,_b,bvals);
	    			myJ[2] += (params[2]>0?1.0:-1.0)*abs(params[0])*(1-*sumf)*isoterm_b_PVM_multi(dir_iter,_a,_b,bvals);	
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
extern "C" __global__ void fit_PVM_multi_kernel(	//INPUT
							const float* 		data, 
							const float* 		params_PVM_single_c,
							const float* 		bvecs, 
							const float* 		bvals, 
							const float		R,
							const float		invR,
							const int 		nvox, 
							const int		ndirections,
							const int 		nfib, 	
							const int		nparams,
							const int		Gamma_for_ball_only,			
							const bool 		m_include_f0,
							const bool		gradnonlin,
							//OUTPUT
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
  	float* x = (float*) &fs[nfib];					//nfib*3
	float* _a = (float*) &x[nfib*3];				//1
	float* _b = (float*) &_a[1];					//1
  	float* sumf = (float*) &_b[1];					//1

	float* C = (float*)&sumf[1];					//nparams*nparams;
	float* el =  (float*)&C[nparams*nparams];			//nparams

	int* indx = (int*)&el[nparams];					//nparams
	int* success = (int*) &indx[nparams];				//1
	int* end = (int*) &success[1];					//1   
	////////// DYNAMIC SHARED MEMORY ///////////

	if(idSubVOX==0){
		
		int nparams_single_c = nparams-1;

		myparams[0] = params_PVM_single_c[(idVOX*nparams_single_c)+0];			//pvm1.get_s0();
  		myparams[1] = 1.0;								//start with d=d_std
  		for(int i=0,ii=3;i<nfib;i++,ii+=3){
    			myparams[ii] = f2x_gpu(params_PVM_single_c[(idVOX*nparams_single_c)+ii-1]);
    			myparams[ii+1] = params_PVM_single_c[(idVOX*nparams_single_c)+ii];
    			myparams[ii+2] = params_PVM_single_c[(idVOX*nparams_single_c)+ii+1];
  		}
		myparams[2] = params_PVM_single_c[(idVOX*nparams_single_c)+1] ; 		//pvm1.get_d();
  		if (m_include_f0)
			myparams[nparams-1]=f2x_gpu(params_PVM_single_c[(idVOX*nparams_single_c)+nparams_single_c-1]);
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
	levenberg_marquardt_PVM_multi_gpu(&data[idVOX*ndirections],&bvecs[pos_bvecs],&bvals[pos_bvals],R,invR, 
	ndirections,nfib,nparams,m_include_f0,idSubVOX,Gamma_for_ball_only,
	step,grad,hess,inverse, pcf,ncf,lambda,cftol,ltol,olambda,success,end,J,
	reduction,fs,x,_a,_b,sumf,C,el,indx,myparams);

	__syncthreads();

  	// finalise parameters
	//m_s0-myparams[0] 	m_d-myparams[1] 	m_d_std-myparams[2]		m_f-m_th-m_ph-myparams[3,4,5,6 etc..]   	m_f0-myparams[nparams-1]

	if(idSubVOX==0){  	
		float aux = myparams[1];

  		myparams[1] = abs(aux*myparams[2]);
		myparams[2] = sqrt(float(abs(aux*myparams[2]*myparams[2])));
  		for(int i=3,k=0;k<nfib;i+=3,k++){
    			myparams[i]  = x2f_gpu(myparams[i]);
  		}
  		if (m_include_f0)
    			myparams[nparams-1]=x2f_gpu(myparams[nparams-1]);

		sort_PVM_multi(nfib,myparams);
  		fix_fsum_PVM_multi(m_include_f0,nfib,nparams,myparams);
	}
	__syncthreads();

	if(idSubVOX<nparams){
		params[(idVOX*nparams)+idSubVOX] = myparams[idSubVOX];
	}
}

//in diffmodel.cc
extern "C" __global__ void get_residuals_PVM_multi_kernel(	//INPUT
								const float* 		data, 
								const float* 		params,
								const float* 		bvecs, 
								const float* 		bvals, 
								const float		R,
								const float		invR,
								const int 		nvox, 
								const int		ndirections,
								const int 		nfib, 
								const int		nparams,
								const int		Gamma_for_ball_only,
								const bool 		m_include_f0,
								const bool		gradnonlin,
								const bool* 		includes_f0,								
								//OUTPUT
								float*			residuals)
{
	int idSubVOX = threadIdx.x;
	int idVOX = blockIdx.x;

	////////// DYNAMIC SHARED MEMORY ///////////
	extern __shared__ double shared[];
	float* myparams = (float*) shared;			//nparams
	float* fs = (float*) &myparams[nparams];		//nfib
  	float* x = (float*) &fs[nfib];				//nfib*3
	float* _a = (float*) &x[nfib*3];			//1
	float* _b = (float*) &_a[1];				//1
  	float* sumf = (float*) &_b[1];				//1
	int* my_include_f0 = (int*) &sumf[1];			//1	
	////////// DYNAMIC SHARED MEMORY ///////////

	float val;
	float predicted_signal;
	float mydata;

	if(idSubVOX==0){
		*my_include_f0 = includes_f0[idVOX];

  		//m_s0-myparams[0]  m_d-myparams[1]  m_d_std-myparams[2]  m_f-m_th-m_ph-myparams[3,4,5,6 etc..]  m_f0-myparams[nparams-1]

  		myparams[0] = params[(idVOX*nparams)+0];
		float aux1 = params[(idVOX*nparams)+1];
		float aux2 = params[(idVOX*nparams)+2];
		
  		myparams[1] = aux1*aux1/aux2/aux2;		//m_d*m_d/m_d_std/m_d_std;
  		myparams[2] = aux2*aux2/aux1;			//m_d_std*m_d_std/m_d; // =1/beta
  		
  		if (*my_include_f0)
    			myparams[nparams-1]=f2x_gpu(params[(idVOX*nparams)+nparams-1]);
	}

	if(idSubVOX<nfib){
		int kk = 3+3*idSubVOX;
		float sinth,costh,sinph,cosph;
	
		myparams[kk]   = f2x_gpu(params[(idVOX*nparams)+kk]);
    		myparams[kk+1] = params[(idVOX*nparams)+kk+1];
    		myparams[kk+2] = params[(idVOX*nparams)+kk+2];

		sincos(myparams[kk+1],&sinth,&costh);
		sincos(myparams[kk+2],&sinph,&cosph);		
    		fs[idSubVOX] = x2f_gpu(myparams[kk]);
    		x[idSubVOX*3] = sinth*cosph;
    		x[idSubVOX*3+1] = sinth*sinph;
    		x[idSubVOX*3+2] = costh;
  	}

	__syncthreads(); 

	if(idSubVOX==0){
  		*_a = abs(myparams[1]);
  		*_b = abs(myparams[2]);
  		*sumf=0;
  		for(int k=0;k<nfib;k++){
	    		*sumf += fs[k];
		}
  	}
  	
	int ndir = ndirections/THREADS_BLOCK_FIT;
	if(idSubVOX<(ndirections%THREADS_BLOCK_FIT)) ndir++;
	
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
      			val += fs[k]*anisoterm_PVM_multi(dir_iter,_a,_b,x2,&bvecs[pos_bvecs],&bvals[pos_bvals],R,invR,ndirections,Gamma_for_ball_only);
    		}	
    		if (*my_include_f0){
      			float temp_f0=x2f_gpu(myparams[nparams-1]);
      			predicted_signal = abs(myparams[0])*(temp_f0+(1-*sumf-temp_f0)*isoterm_PVM_multi(dir_iter,_a,_b,&bvals[pos_bvals])+val);
    		}else{
      			predicted_signal = abs(myparams[0])*((1-*sumf)*isoterm_PVM_multi(dir_iter,_a,_b,&bvals[pos_bvals])+val); 
  		}   

		//residuals=m_data-predicted_signal;
		residuals[idVOX*ndirections+dir_iter]= mydata - predicted_signal;

		dir_iter+=THREADS_BLOCK_FIT;
	}
}
