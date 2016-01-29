/*  runmcmc_kernels.cu

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

#include <iostream>
#include <fstream>
#include <stdio.h>
#include "fibre_gpu.h"
#include <math.h>
#include <string.h>
#include <string>
#include <cuda_runtime.h>
#include <cuda.h>
#include <curand.h>

#include <options.h>

#define maxfloat 1e10
#define UPPERDIFF 0.005

__device__ inline void propose(float* param, float* old, float prop, float random){
	*old=*param;
	*param = *param + random*prop;
}
__device__ inline void reject(float* param, float* prior, float* old, float* m_prior_en, float* m_old_prior_en, float* m_energy, float* m_old_energy, int* rej){
	*param=old[0];
	*prior=old[1];
	*m_prior_en=*m_old_prior_en;
	*rej=*rej+1;
	//restore_energy()
      	*m_energy=*m_old_energy;
}
__device__ inline void rejectF(float* param, float* prior, float* old, float* m_prior_en, float* m_old_prior_en, float* fm_prior_en, float* fm_old_prior_en, float* m_energy, float* m_old_energy, int* rej){
	*param=old[0];
	*prior=old[1];
	*fm_prior_en=*fm_old_prior_en;
	*m_prior_en=*m_old_prior_en;
	*rej=*rej+1;
	//restore_energy()
      	*m_energy=*m_old_energy;
}				
__device__ inline void getfsum(float* fsum, float* m_f, float m_f0, int nfib){
	*fsum=m_f0;
	for(int f=0;f<nfib;f++){  
		*fsum = *fsum + m_f[f];
	}
}
__device__ inline bool compute_test_energy(float *m_energy, float* m_old_energy, float m_prior_en, float m_likelihood_en, float random){
	*m_old_energy=*m_energy;
      	*m_energy=m_prior_en+m_likelihood_en;

	double tmp=exp(double(*m_old_energy-*m_energy));
	return (tmp>random);
}
__device__ inline void compute_signal(double *signals,double *oldsignals,float mbvals,float* m_d, float* m_dstd, float* m_R, float angtmp, int model){
	*oldsignals=*signals;
	if(model==1 || (*m_dstd<1e-5 && model==2)){	
		*signals=exp(double(-*m_d*mbvals*angtmp));
	}else if(model==2){
		//double dbeta= *m_d/(*m_dstd**m_dstd);
	 	//double dalpha= *m_d*dbeta;   
		//*signals=exp(double(log(double(dbeta/(dbeta+mbvals*angtmp)))*dalpha));   
		float sig2=*m_dstd**m_dstd;
	 	float dalpha=*m_d**m_d/sig2;      
		*signals=exp(log(double(*m_d/(*m_d + mbvals*angtmp*sig2)))*dalpha); // more stable
	}else if(model==3){
		float invR=1.0f/(2.0f**m_R+1.0f);
	   	*signals=exp(-mbvals*3**m_d*invR*((1-*m_R)*angtmp+*m_R));
       	}
}
__device__ inline void compute_iso_signal(double *isosignals,double *oldisosignals, float mbvals,float* m_d, float* m_dstd, int model){
	*oldisosignals=*isosignals;
	if(model==1 || *m_dstd<1e-5){
	 	*isosignals=exp(double(-m_d[0]*mbvals));	
	}else if(model>=2){
		//double dbeta= *m_d/(*m_dstd**m_dstd);
	  	//double dalpha= *m_d*dbeta;
		//*isosignals=exp(double(log(double(dbeta/(dbeta+mbvals)))*dalpha));
		float sig2=*m_dstd**m_dstd;
		float dalpha=*m_d**m_d/sig2;	
		*isosignals=exp(log(double(*m_d/(*m_d+mbvals*sig2)))*dalpha); // more numerically stable
	}
}
__device__ inline void restore_signals(double* signals, double* oldsignals, int idVOX, int idSubVOX, int mydirs, int nfib, int ndirections, int threadsBlock){
	for(int f=0;f<nfib;f++){
		for(int i=0; i<mydirs; i++){
			int pos = idVOX*ndirections*nfib + f*ndirections + idSubVOX + i*threadsBlock;
			signals[pos] = oldsignals[pos];
		}	
	}
}
__device__ inline void restore_isosignals(double* isosignals, double* oldisosignals, int idVOX, int idSubVOX, int mydirs, int ndirections, int threadsBlock){
	for(int i=0; i<mydirs; i++){
		int pos = idVOX*ndirections + idSubVOX + i*threadsBlock;
		isosignals[pos]=oldisosignals[pos];
	}
}
__device__ inline void restore_angtmp_signals(double* signals, double* oldsignals,double* angtmp, double* oldangtmp, int idVOX, int idSubVOX, int mydirs, int nfib, int fibre, int ndirections, int threadsBlock){
	for(int i=0; i<mydirs; i++){
		int pos = idVOX*ndirections*nfib + fibre*ndirections + idSubVOX + i*threadsBlock;
		int pos2 = idVOX*ndirections + idSubVOX + i*threadsBlock;
		angtmp[pos]=oldangtmp[pos2];
		signals[pos] = oldsignals[pos];	
	}
}
__device__  inline void compute_prior(float *m_prior_en, float *m_prior_en_old,float* m_d_prior,float* m_S0_prior,float *m_prior_enf, float* m_f0_prior, float* m_tau_prior, float* m_dstd_prior, float* m_R_prior, int nfib){			
        *m_prior_en_old=*m_prior_en;
	*m_prior_en=*m_d_prior+*m_S0_prior+*m_dstd_prior+*m_R_prior+*m_tau_prior+*m_f0_prior;
	for(int f=0;f<nfib;f++){
		*m_prior_en=*m_prior_en+m_prior_enf[f];
	}	
}

__device__ inline float logIo(const float& x){
    	float y,b;
    	b= fabs(x);
    	if (b<3.75){
      		float a=x/3.75;
      		a*=a;
      		//Bessel function evaluation
		y=1.0+a*(3.5156229+a*(3.0899424+a*(1.2067492+a*(0.2659732+a*(0.0360768+a*0.0045813)))));
      		y=log(double(y));
    	}else{
      		float a=3.75/b; 
      		//Bessel function evaluation
      		//y=(exp(b)/sqrt(b))*(0.39894228+a*(0.01328592+a*(0.00225319+a*(-0.00157565+a*(0.00916281+a*(-0.02057706+a*(0.02635537+a*(-0.01647633+a*0.00392377))))))));
      		//Logarithm of Bessel function

		y=b+log(double((0.39894228+a*(0.01328592+a*(0.00225319+a*(-0.00157565+a*(0.00916281+a*(-0.02057706+a*(0.02635537+a*(-0.01647633+a*0.00392377))))))))/sqrt(b)));
    	}

    	return y;
}
__device__ inline void compute_likelihood(int idSubVOX,float* m_S0,float *m_likelihood_en,float *m_f,double *signals,double *isosignals,const float *mdata,float* fsum,double *reduction, float* m_f0, const bool rician, float* m_tau,int mydirs, int threadsBlock, int ndirections, int nfib){
	
	double pred;
	int pos;

	reduction[idSubVOX]=0;
	for(int i=0; i<mydirs; i++){
		pred=0;
	      	for(int f=0;f<nfib;f++){
			pos = f*ndirections + idSubVOX + i*threadsBlock;
			pred= pred+m_f[f]*signals[pos];
	     	}
		pos = idSubVOX + i*threadsBlock;
		pred= *m_S0*(pred+(1-*fsum)*isosignals[pos]+*m_f0); //F0
	
		if(!rician){
			double diff = mdata[pos]-pred;
			reduction[idSubVOX] = reduction[idSubVOX]+(diff*diff);
		}else{
			pred= log(mdata[pos])+(-0.5**m_tau*(mdata[pos]*mdata[pos]+pred*pred)+logIo(*m_tau*pred*mdata[pos]));  
			reduction[idSubVOX] = reduction[idSubVOX]+pred;
		}
	}

	__syncthreads();

	unsigned int s2=threadsBlock;
	for(unsigned int s=threadsBlock>>1; s>0; s>>=1) {
		if((s2%2)&&(idSubVOX==(s-1))) reduction[idSubVOX]= reduction[idSubVOX] + reduction[idSubVOX + s +1]; 
        	if (idSubVOX < s){
            		reduction[idSubVOX] = reduction[idSubVOX] + reduction[idSubVOX + s];
       	 	}
		s2=s;
        	__syncthreads();
    	}
	if(idSubVOX==0){
		double sumsquares=0;
		sumsquares+=reduction[0];
		if(!rician){ 
		 	*m_likelihood_en=(ndirections/2.0)*log(sumsquares/2.0); 
		}else{
			*m_likelihood_en= -ndirections*log(*m_tau)-sumsquares;
		}
	}
}
			  
extern "C" __global__ void init_Fibres_Multifibres_kernel(	//INPUT
								const float*			datam,
								const float*			params,
								const float*			tau,
								const float*			bvals,
								const double*			alpha,
								const double*			beta,
								const float			R_priormean,
								const float			R_priorstd,
								const float			R_priorfudge,		
								const int			ndirections,
								const int 			nfib,
								const int 			nparams_fit,
								const int 			model,
								const float 			fudgevalue,
								const bool			m_includef0,
								const bool			rician,
								const bool 			m_ardf0,	// opts.ardf0.value()
								const bool 			ard_value,	// opts.all_ard.value()
								const bool 			no_ard_value,	// opts.no_ard.value()
								const bool			gradnonlin,
								//TO USE
								double*				angtmp,
								//OUTPUT
								FibreGPU*			fibres,
								MultifibreGPU*			multifibres,
								double*				signals,
								double*				isosignals)
{
	int idSubVOX= threadIdx.x;
	int threadsBlock = blockDim.x;
	int idVOX= blockIdx.x;
	bool leader = (idSubVOX==0);

	////////// DYNAMIC SHARED MEMORY ///////////
	extern __shared__ double shared[];
	double* reduction = (double*)shared;				//threadsBlock

	float* m_S0 = (float*) &reduction[threadsBlock];		//1
	float* m_d = (float*) &m_S0[1];					//1
	float* m_dstd =(float*) &m_d[1];				//1
	float* m_R =(float*) &m_dstd[1];				//1	
	float* m_f0 = (float*) &m_R[1];					//1
	float* m_tau = (float*) &m_f0[1];				//1
	float* m_th = (float*) &m_tau[1];				//nfib
	float* m_ph = (float*) &m_th[nfib];				//nfib
	float* m_f = (float*) &m_ph[nfib];				//nfib

	float* fsum = (float*) &m_f[nfib];				//1
	float* m_likelihood_en = (float*) &fsum[1];			//1
	float* m_prior_en = (float*) &m_likelihood_en[1];		//1

	int* posBV = (int*) &m_prior_en[1];				//1
	////////// DYNAMIC SHARED MEMORY ///////////
	
	// m_s0-params[0]	m_d-params[1]	m_f-m_th-m_ph-params[add+2,3,4,5, etc..]	m_f0-params[nparams-1]
	if(leader){
		if(gradnonlin) *posBV = (idVOX*ndirections);
		else *posBV = 0;

		*m_S0 = params[idVOX*nparams_fit];
		multifibres[idVOX].m_S0 = *m_S0;
		multifibres[idVOX].m_S0_prior = 0;
		multifibres[idVOX].m_S0_acc = 0;
		multifibres[idVOX].m_S0_rej = 0;
	
		*m_d=params[idVOX*nparams_fit+1];
		if(*m_d<0 || *m_d> UPPERDIFF) *m_d=2e-3;			//this is in xfibres...after fit
		multifibres[idVOX].m_d = *m_d;
		multifibres[idVOX].m_d_prior = 0;
		multifibres[idVOX].m_d_acc = 0;
		multifibres[idVOX].m_d_rej = 0;

		if(model>=2){ 
			*m_dstd=params[idVOX*nparams_fit+2];
			float upper_d_std=0.01;
			if (model==3) upper_d_std=0.004;
      			if(*m_dstd<0 || *m_dstd>upper_d_std) *m_dstd=*m_d/10; 	//this is in xfibres...after fit
			if (model==3){ 
				*m_R=R_priormean;	
			}else{ 
				*m_R=0;
			}
		}
		else *m_dstd = 0;
		multifibres[idVOX].m_dstd = *m_dstd;
		multifibres[idVOX].m_dstd_prior = 0;
		multifibres[idVOX].m_dstd_acc = 0;
		multifibres[idVOX].m_dstd_rej = 0;

		multifibres[idVOX].m_R = *m_R;
		multifibres[idVOX].m_R_prior = 0;
		multifibres[idVOX].m_R_acc = 0;
		multifibres[idVOX].m_R_rej = 0;

		if (m_includef0) *m_f0=params[idVOX*nparams_fit+nparams_fit-1];
		else *m_f0=0;
		multifibres[idVOX].m_f0 = *m_f0;
		multifibres[idVOX].m_f0_prior = 0;
		multifibres[idVOX].m_f0_acc = 0;
		multifibres[idVOX].m_f0_rej = 0;

		*m_tau = tau[idVOX];
		multifibres[idVOX].m_tau = *m_tau;
		multifibres[idVOX].m_tau_prior = 0;
		multifibres[idVOX].m_tau_acc = 0;
		multifibres[idVOX].m_tau_rej = 0;
	}
	__syncthreads();

	int mydirs = ndirections/threadsBlock;
	int mod = ndirections%threadsBlock;
	if(mod&&(idSubVOX<mod)) mydirs++;

	//------ Fibre constructor ------
	if(idSubVOX<nfib){
		int add=0;
		if(model>=2) add=1;		// if model 2 we have d_std and then 1 more parameter in position 2
		int pos = (idVOX*nfib)+idSubVOX;

		m_th[idSubVOX]=params[idVOX*nparams_fit+2+3*idSubVOX+1+add];
		fibres[pos].m_th = m_th[idSubVOX];
		fibres[pos].m_th_prop = 0.2;
		float m_th_prior = 0;
		fibres[pos].m_th_acc = 0;
		fibres[pos].m_th_rej = 0;
		
		//compute_th_prior();
	      	if(m_th==0){
			m_th_prior=0;
		}else{
			m_th_prior=-log(double(fabs(sin(double(m_th[idSubVOX]))/2)));
	      	}
		fibres[pos].m_th_prior = m_th_prior;
		
		float m_ph_prior=0;	//compute_ph_prior();
		m_ph[idSubVOX]=params[idVOX*nparams_fit+2+3*idSubVOX+2+add];
		fibres[pos].m_ph = m_ph[idSubVOX];
		fibres[pos].m_ph_prop = 0.2;
		fibres[pos].m_ph_prior = 0;	//compute_ph_prior();
		fibres[pos].m_ph_acc = 0;
		fibres[pos].m_ph_rej = 0;

		m_f[idSubVOX] = params[idVOX*nparams_fit+2+3*idSubVOX+add]; 
		fibres[pos].m_f=m_f[idSubVOX];
		fibres[pos].m_f_prop = 0.2;
		float m_f_prior = 0;
		fibres[pos].m_f_acc = 0;
		fibres[pos].m_f_rej = 0;
			
		if(idSubVOX==0){
			fibres[pos].m_lam_jump = ard_value;
		}else{
			fibres[pos].m_lam_jump = !no_ard_value;
		}

		//compute_f_prior();
      		if (m_f[idSubVOX]<=0 | m_f[idSubVOX]>=1 ){
      		}else{
	  		if(fibres[pos].m_lam_jump){              
	    			m_f_prior=log(double(m_f[idSubVOX]));
	  		}else{
	    			m_f_prior=0;
			}
			m_f_prior= fudgevalue* m_f_prior;
      		}
		fibres[pos].m_f_prior = m_f_prior;

		//fibres[vox].m_lam = m_lam; ??
		//fibres[vox].m_lam_prop = 1;
		//fibres[vox].m_lam_prior = 0;
		//compute_lam_prior();

		//compute_prior();
		fibres[pos].m_prior_en= m_th_prior + m_ph_prior + m_f_prior;
	}

	__syncthreads();

	//compute_signal_pre
	for(int f=0;f<nfib;f++){	
		for(int i=0; i<mydirs; i++){
			double myalpha = alpha[*posBV+idSubVOX+i*threadsBlock];
			double cos_alpha_minus_theta=cos(double(myalpha-m_th[f]));   
			double cos_alpha_plus_theta=cos(double(myalpha+m_th[f]));
			int pos = idVOX*ndirections*nfib + f*ndirections + idSubVOX + i*threadsBlock;
			double aux = (cos(double(m_ph[f]-beta[*posBV+idSubVOX+i*threadsBlock]))*(cos_alpha_minus_theta-cos_alpha_plus_theta)/2)+(cos_alpha_minus_theta+cos_alpha_plus_theta)/2;
		     	aux =  aux*aux;
		 	angtmp[pos]= aux;
		}
	}
	//------ Fibre constructor ------
	//compute_signal()
	double old;
	for(int f=0;f<nfib;f++){
		for(int i=0; i<mydirs; i++){
			int pos = idVOX*ndirections*nfib + f*ndirections + idSubVOX + i*threadsBlock;
			compute_signal(&signals[pos],&old,bvals[*posBV+idSubVOX+i*threadsBlock],m_d,m_dstd,m_R,angtmp[pos],model);
		}
	}

	//------ initialise_energies ------
	if(leader){
		getfsum(fsum,m_f,*m_f0,nfib);
		
	      	//compute_d_prior(); 
		if(*m_d>=0 && *m_d<=UPPERDIFF){
			if (model==3){
 	          		//float alpha=3.0; float beta=4000;  //Gamma_prior around 0.5-1E-3
 	          		multifibres[idVOX].m_d_prior =(1.0f-3.0f)*log(*m_d)+4000.0f**m_d;
 	        	}
       		}

	      	if(model>=2){
			//compute_d_std_prior();
			float upper_d_std=0.01;
			if (model==3) upper_d_std=0.004;
			if(*m_dstd>0 && *m_dstd<=upper_d_std){
				multifibres[idVOX].m_dstd_prior=log(*m_dstd);
			}
			if (model==3){
	  			//compute_R_prior();
				float upper_R=2.0f*R_priormean;
				float lower_R=R_priormean-2.0f*R_priorstd;
				if (R_priormean>0.5f)
					upper_R=1.0f;
				if (lower_R<0.0f)
					lower_R=1E-8f;
				if (R_priorfudge>0.0f && *m_d>UPPERDIFF/2.0f){
					//then use an ARD prior to avoid competition with the isotropic compartments
 	        			if (*m_R>=1E-8f && *m_R<=upper_R){
 	          				multifibres[idVOX].m_R_prior=R_priorfudge*log(*m_R);
 	        			}
 	      			}else{	
      					if(*m_R>lower_R && *m_R<=upper_R){
						float Rstd2=R_priorstd*R_priorstd; 
						multifibres[idVOX].m_R_prior=(*m_R-R_priormean)*(*m_R-R_priormean)/Rstd2;  //Gaussian prior
      					}
				}
			}
		}
	      	//compute_tau_prior(); m_tau_prior=0; so it doesn't do nothing, it is already 0
	      	if (m_includef0){
			//compute_f0_prior();
			if (*m_f0<=0 || *m_f0>=1){
	      		}else{
				if(!m_ardf0){}     	//Without ARD
				else              	//With ARD
		  			multifibres[idVOX].m_f0_prior= log(double(*m_f0));
	      		}
		}
	      	//compute_S0_prior(); m_S0_prior=0; so i don't do nothing, it is already 0
		//*m_prior_en = 0;
	      	//compute_prior();
	      	*m_prior_en=multifibres[idVOX].m_d_prior+multifibres[idVOX].m_S0_prior;
	      	if(model>=2)
			*m_prior_en= *m_prior_en+multifibres[idVOX].m_dstd_prior;
		if(model==3)
			*m_prior_en= *m_prior_en+multifibres[idVOX].m_R_prior;
	      	//if(m_rician) m_prior_en=m_prior_en+m_tau_prior; is 0
	      	if (m_includef0)
			*m_prior_en=*m_prior_en+multifibres[idVOX].m_f0_prior;
	      	for(int fib=0;fib<nfib;fib++){
			*m_prior_en=*m_prior_en+ fibres[idVOX*nfib+fib].m_prior_en;
	      	} 
		multifibres[idVOX].m_prior_en = *m_prior_en;
	}
	//------ initialise_energies ------
	//compute_iso_signal()
	for(int i=0; i<mydirs; i++){
		int pos = idVOX*ndirections + idSubVOX + i*threadsBlock;	
		compute_iso_signal(&isosignals[pos],&old,bvals[*posBV+idSubVOX+i*threadsBlock],m_d,m_dstd,model);
	}		
 
	__syncthreads();

	//------ initialise_energies ------
	//compute_likelihood()
	compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,&signals[idVOX*nfib*ndirections],&isosignals[idVOX*ndirections],&datam[idVOX*ndirections],fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);

	__syncthreads();

	if(leader){
		multifibres[idVOX].m_likelihood_en = *m_likelihood_en;
	      	//------ initialise_energies ------
		//compute_energy();	
		multifibres[idVOX].m_energy = *m_prior_en+*m_likelihood_en;

	    	//initialise_props();
	      	multifibres[idVOX].m_S0_prop=multifibres[idVOX].m_S0/10.0; 
	      	multifibres[idVOX].m_d_prop=*m_d/10.0;
	      	multifibres[idVOX].m_dstd_prop=*m_dstd/10.0;
	      	multifibres[idVOX].m_tau_prop=*m_tau/2.0;
	      	multifibres[idVOX].m_f0_prop=0.2;
		multifibres[idVOX].m_R_prop=*m_R/10.0;
	}
}

extern "C" __global__ void runmcmc_kernel(	//INPUT 
						const float*			datam,
						const float*			bvals,
						const double*			alpha,
						const double*			beta,
						float*				randomsN,
						float*				randomsU,
						const float			R_priormean,
						const float			R_priorstd,	
						const float			R_priorfudge,			
						const int			ndirections,
						const int			nfib,
						const int			nparams,
						const int 			model,
						const float 			fudgevalue,
						const bool 			m_include_f0,
						const bool 			m_ardf0,
						const bool 			can_use_ard, 
						const bool 			rician,
						const bool			gradnonlin,
						const int 			updateproposalevery, 	//update every this number of iterations	
						const int 			iterations,		//num of iterations to do this time (maybe is a part of the total)
						const int 			current_iter,		//the number of the current iteration over the total iterations
						const int 			iters_burnin,		//iters in burin, we need it to continue the updates at the correct time. 
						const int 			record_every, 		//record every this number
						const int 			totalrecords,		//total number of records to do
						//TO USE
						double*				oldsignals,
						double*				oldisosignals,
						double*				angtmp,
						double*				oldangtmp,
						//INPUT-OUTPUT
						FibreGPU*			fibres,
						MultifibreGPU*			multifibres,
						double*				signals,
						double*				isosignals,
						//OUTPUT
						float*				rf0,			//record of parameters
						float*				rtau,
						float*				rs0,
						float*				rd,
						float*				rdstd,
						float*				rR,
						float*				rth,
						float*				rph, 
						float*				rf)
{	
	int idSubVOX= threadIdx.x;
	int threadsBlock = blockDim.x;
	bool leader = (idSubVOX==0);

	////////// DYNAMIC SHARED MEMORY ///////////
	extern __shared__ double shared[];
	double* reduction = (double*)shared;				//threadsBlock

	float* m_S0 = (float*) &reduction[threadsBlock];		//1
	float* m_d = (float*) &m_S0[1];					//1
	float* m_dstd =(float*) &m_d[1];				//1
	float* m_R =(float*) &m_dstd[1];				//1	
	float* m_f0 = (float*) &m_R[1];					//1
	float* m_tau = (float*) &m_f0[1];				//1
	float* m_th = (float*) &m_tau[1];				//nfib
	float* m_ph = (float*) &m_th[nfib];				//nfib
	float* m_f = (float*) &m_ph[nfib];				//nfib

	float* m_S0_prior = (float*) &m_f[nfib];			//1
	float* m_d_prior = (float*) &m_S0_prior[1];			//1
	float* m_dstd_prior = (float*) &m_d_prior[1];			//1
	float* m_R_prior = (float*) &m_dstd_prior[1];			//1	
	float* m_f0_prior = (float*) &m_R_prior[1];			//1
	float* m_tau_prior = (float*) &m_f0_prior[1];			//1
	float* m_th_prior = (float*) &m_tau_prior[1];			//nfib
	float* m_ph_prior = (float*) &m_th_prior[nfib];			//nfib
	float* m_f_prior = (float*) &m_ph_prior[nfib];			//nfib

	float* m_S0_prop = (float*) &m_f_prior[nfib];			//1
	float* m_d_prop = (float*) &m_S0_prop[1];			//1
	float* m_dstd_prop = (float*) &m_d_prop[1];			//1
	float* m_R_prop = (float*) &m_dstd_prop[1];			//1
	float* m_f0_prop = (float*) &m_R_prop[1];			//1
	float* m_tau_prop = (float*) &m_f0_prop[1];			//1
	float* m_th_prop = (float*) &m_tau_prop[1];			//nfib
	float* m_ph_prop = (float*) &m_th_prop[nfib];			//nfib
	float* m_f_prop = (float*) &m_ph_prop[nfib];			//nfib

	float* fsum = (float*) &m_f_prop[nfib];				//1
	float* m_likelihood_en = (float*) &fsum[1];			//1
	float* m_prior_en = (float*) &m_likelihood_en[1];		//1
	float* m_old_prior_en = (float*) &m_prior_en[1];		//1
	float* fm_prior_en = (float*) &m_old_prior_en[1];		//nfib		
	float* fm_old_prior_en = (float*) &fm_prior_en[nfib];		//1		
	float* m_energy = (float*) &fm_old_prior_en[1];			//1
	float* m_old_energy = (float*) &m_energy[1];			//1
	float* old = (float*) &m_old_energy[1];				//2

	float* prerandN = (float*) &old[2];				//nparams
	float* prerandU = (float*) &prerandN[nparams];			//nparams

	int* m_S0_acc = (int*) &prerandU[nparams];			//1
	int* m_d_acc = (int*) &m_S0_acc[1];				//1
	int* m_dstd_acc = (int*) &m_d_acc[1];				//1
	int* m_R_acc = (int*) &m_dstd_acc[1];				//1
	int* m_f0_acc = (int*) &m_R_acc[1];				//1
	int* m_tau_acc = (int*) &m_f0_acc[1];				//1
	int* m_th_acc = (int*) &m_tau_acc[1];				//nfib
	int* m_ph_acc = (int*) &m_th_acc[nfib];				//nfib
	int* m_f_acc = (int*) &m_ph_acc[nfib];				//nfib

	int* m_S0_rej = (int*) &m_f_acc[nfib];				//1
	int* m_d_rej = (int*) &m_S0_rej[1];				//1
	int* m_dstd_rej = (int*) &m_d_rej[1];				//1	
	int* m_R_rej = (int*) &m_dstd_rej[1];				//1	
	int* m_f0_rej = (int*) &m_R_rej[1];				//1
	int* m_tau_rej = (int*) &m_f0_rej[1];				//1	
	int* m_th_rej = (int*) &m_tau_rej[1];				//nfib
	int* m_ph_rej = (int*) &m_th_rej[nfib];				//nfib
	int* m_f_rej = (int*) &m_ph_rej[nfib];				//nfib
	
	int* rejflag = (int*) &m_f_rej[nfib];				//3					
	int* m_lam_jump = (int*) &rejflag[3];				//nfib	
	int* idVOX = (int*) &m_lam_jump[nfib];				//1		
	int* count_update = (int*) &idVOX[1];				//1	
	int* recordcount = (int*) &count_update[1];			//1
	int* sample = (int*) &recordcount[1];				//1		
	int* localrand = (int*) &sample[1];				//1
	int* posBV = (int*) &localrand[1];				//1
	////////// DYNAMIC SHARED MEMORY ///////////
	
	if (leader){
		*idVOX= blockIdx.x;
		*count_update = current_iter+iters_burnin;	//count for updates
		*recordcount = current_iter;	
		if(record_every) *sample=1+(current_iter/record_every);		//the next number of sample.....the index start in 0

		if(gradnonlin)*posBV = (*idVOX*ndirections);
		else *posBV = 0;

		*m_prior_en=multifibres[*idVOX].m_prior_en;

		if(model>=2){
			*m_dstd_acc=multifibres[*idVOX].m_dstd_acc;
			*m_dstd_rej=multifibres[*idVOX].m_dstd_rej;
			*m_dstd_prior=multifibres[*idVOX].m_dstd_prior;
			*m_dstd_prop=multifibres[*idVOX].m_dstd_prop;
			*m_dstd=multifibres[*idVOX].m_dstd;
		}else{
			*m_dstd_acc=0;
			*m_dstd_rej=0;
			*m_dstd_prior=0;
			*m_dstd_prop=0;
			*m_dstd=0;
		}
		if(model==3){
			*m_R_acc=multifibres[*idVOX].m_R_acc;
			*m_R_rej=multifibres[*idVOX].m_R_rej;
			*m_R_prior=multifibres[*idVOX].m_R_prior;
			*m_R_prop=multifibres[*idVOX].m_R_prop;
			*m_R=multifibres[*idVOX].m_R;
		}else{
			*m_R_acc=0;
			*m_R_rej=0;
			*m_R_prior=0;
			*m_R_prop=0;
			*m_R=0;
		}
	
		*m_d=multifibres[*idVOX].m_d;
		*m_energy=multifibres[*idVOX].m_energy;
		*m_d_prop=multifibres[*idVOX].m_d_prop;
		*m_d_prior=multifibres[*idVOX].m_d_prior;
		*m_S0_prior=multifibres[*idVOX].m_S0_prior;
		*m_S0=multifibres[*idVOX].m_S0;
		*m_likelihood_en=multifibres[*idVOX].m_likelihood_en;
		*m_d_acc=multifibres[*idVOX].m_d_acc;
		*m_d_rej=multifibres[*idVOX].m_d_rej;
		*m_S0_acc=multifibres[*idVOX].m_S0_acc;
		*m_S0_rej=multifibres[*idVOX].m_S0_rej;
		*m_S0_prop=multifibres[*idVOX].m_S0_prop;

		if(m_include_f0){
			*m_f0_acc=multifibres[*idVOX].m_f0_acc;
			*m_f0_rej=multifibres[*idVOX].m_f0_rej;
			*m_f0_prop=multifibres[*idVOX].m_f0_prop;
			*m_f0_prior=multifibres[*idVOX].m_f0_prior;
			*m_f0=multifibres[*idVOX].m_f0;
		}else{ 
			*m_f0_acc=0;
			*m_f0_rej=0;
			*m_f0_prop=0;
			*m_f0_prior=0;
			*m_f0=0;
		}
				
		if(rician){
			*m_tau_acc=multifibres[*idVOX].m_tau_acc;
			*m_tau_rej=multifibres[*idVOX].m_tau_rej;
			*m_tau_prop=multifibres[*idVOX].m_tau_prop;
			*m_tau_prior=multifibres[*idVOX].m_tau_prior;
			*m_tau=multifibres[*idVOX].m_tau;	
		}else{ 
			*m_tau_acc=0;
			*m_tau_rej=0;
			*m_tau_prop=0;
			*m_tau_prior=0;
			*m_tau=0;
		}
	}

	__syncthreads();

	int mydirs = ndirections/threadsBlock;
	int mod = ndirections%threadsBlock;
	if(mod&&(idSubVOX<mod)) mydirs++;

	if(idSubVOX<nfib){
		int pos = (*idVOX*nfib)+idSubVOX;
		m_th[idSubVOX]=fibres[pos].m_th;
		m_ph[idSubVOX]=fibres[pos].m_ph;
		m_f[idSubVOX]=fibres[pos].m_f;
	
		m_th_acc[idSubVOX]=fibres[pos].m_th_acc;
		m_th_rej[idSubVOX]=fibres[pos].m_th_rej;
		m_ph_acc[idSubVOX]=fibres[pos].m_ph_acc;
		m_ph_rej[idSubVOX]=fibres[pos].m_ph_rej;
		m_f_acc[idSubVOX]=fibres[pos].m_f_acc;
		m_f_rej[idSubVOX]=fibres[pos].m_f_rej;

		fm_prior_en[idSubVOX]=fibres[pos].m_prior_en;
		m_th_prior[idSubVOX]=fibres[pos].m_th_prior;
		m_ph_prior[idSubVOX]=fibres[pos].m_ph_prior;
		m_f_prior[idSubVOX]=fibres[pos].m_f_prior;

		m_th_prop[idSubVOX]=fibres[pos].m_th_prop;
		m_ph_prop[idSubVOX]=fibres[pos].m_ph_prop;
		m_f_prop[idSubVOX]=fibres[pos].m_f_prop;

		m_lam_jump[idSubVOX]=fibres[pos].m_lam_jump;		
	}
	__syncthreads();
		
	//compute_signal_pre
	for(int f=0;f<nfib;f++){
		for(int i=0; i<mydirs; i++){	
			double myalpha = alpha[*posBV+idSubVOX+i*threadsBlock];		
			double cos_alpha_minus_theta=cos(double(myalpha-m_th[f]));   
			double cos_alpha_plus_theta=cos(double(myalpha+m_th[f]));
			int pos = *idVOX*ndirections*nfib + f*ndirections + idSubVOX + i*threadsBlock;
			double aux = (cos(double(m_ph[f]-beta[*posBV+idSubVOX+i*threadsBlock]))*(cos_alpha_minus_theta-cos_alpha_plus_theta)/2)+(cos_alpha_minus_theta+cos_alpha_plus_theta)/2;
		     	aux =  aux*aux;
		 	angtmp[pos]= aux;
		}
	}

	if (leader) getfsum(fsum,m_f,*m_f0,nfib);
	for (int niter=0; niter<iterations; niter++){
		//code jump()

		//prefetching randoms
		if (leader){
			*count_update=*count_update+1;
			*recordcount=*recordcount+1;
			int idrand = *idVOX*iterations*nparams+(niter*nparams);
			*localrand = 0;
			if(m_include_f0){
				prerandN[*localrand]=randomsN[idrand];
				prerandU[*localrand]=randomsU[idrand];
				idrand++;
				*localrand=*localrand+1;
			}
			if(rician){
				prerandN[*localrand]=randomsN[idrand];
				prerandU[*localrand]=randomsU[idrand];
				idrand++;
				*localrand=*localrand+1;
			}
			//d
			prerandN[*localrand]=randomsN[idrand];
			prerandU[*localrand]=randomsU[idrand];
			idrand++;
			*localrand=*localrand+1;
			//d_std
			if(model>=2){
				prerandN[*localrand]=randomsN[idrand];
				prerandU[*localrand]=randomsU[idrand];
				idrand++;
				*localrand=*localrand+1;
			}
			//R
			if(model==3){
				prerandN[*localrand]=randomsN[idrand];
				prerandU[*localrand]=randomsU[idrand];
				idrand++;
				*localrand=*localrand+1;
			}
			//S0
			prerandN[*localrand]=randomsN[idrand];
			prerandU[*localrand]=randomsU[idrand];
			idrand++;
			*localrand=*localrand+1;

			for(int f=0;f<nfib;f++){
				prerandN[*localrand]=randomsN[idrand];
				prerandU[*localrand]=randomsU[idrand];
				idrand++;
				*localrand=*localrand+1;

				prerandN[*localrand]=randomsN[idrand];
				prerandU[*localrand]=randomsU[idrand];
				idrand++;
				*localrand=*localrand+1;

				prerandN[*localrand]=randomsN[idrand];					
				prerandU[*localrand]=randomsU[idrand];	
				idrand++;	
				*localrand=*localrand+1;	
			}
			*localrand = 0;
		}
////////////////////////////////////////////////////////////////// F0
		if(m_include_f0){
			if (leader){
				propose(m_f0,old,*m_f0_prop,prerandN[*localrand]);
				//compute_f0_prior()     
				old[1]=*m_f0_prior;
	      			if(*m_f0<=0 || *m_f0 >=1){ 
					rejflag[0]=true;
				}else{ 	
					rejflag[0]=false;
					if(!m_ardf0){
						*m_f0_prior=0;
	      				}else{
						*m_f0_prior=log(double(*m_f0));
					}
				}
				getfsum(fsum,m_f,*m_f0,nfib);
				//compute_prior()
				compute_prior(m_prior_en,m_old_prior_en,m_d_prior,m_S0_prior,fm_prior_en,m_f0_prior,m_tau_prior,m_dstd_prior,m_R_prior,nfib);
				//reject_f_sum()
				rejflag[1]=(*fsum>1);
			}
			__syncthreads();
			//compute_likelihood()
			compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,&signals[*idVOX*nfib*ndirections],&isosignals[*idVOX*ndirections],&datam[*idVOX*ndirections],fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);
			__syncthreads();

			if (leader){
				rejflag[2]=compute_test_energy(m_energy,m_old_energy,*m_prior_en,*m_likelihood_en,prerandU[*localrand]);
				*localrand=*localrand+1;	
				if(!rejflag[0]){
					if(!rejflag[1]){
						if(rejflag[2]){
		  					*m_f0_acc=*m_f0_acc+1;   
						}else{
							reject(m_f0,m_f0_prior,old,m_prior_en,m_old_prior_en,m_energy,m_old_energy,m_f0_rej);
							getfsum(fsum,m_f,*m_f0,nfib);
						}
					}else{
						reject(m_f0,m_f0_prior,old,m_prior_en,m_old_prior_en,m_energy,m_old_energy,m_f0_rej);
						getfsum(fsum,m_f,*m_f0,nfib);
					}
				}else{ 
					reject(m_f0,m_f0_prior,old,m_prior_en,m_old_prior_en,m_energy,m_old_energy,m_f0_rej);
					getfsum(fsum,m_f,*m_f0,nfib);
				}
			}
		}
////////////////////////////////////////////////////////////////// TAU
		if(rician){
			if (leader){
				propose(m_tau,old,*m_tau_prop,prerandN[*localrand]);
				//compute_tau_prior()     
				old[1]=*m_tau_prior;
	      			if(*m_tau<=0){ 
					rejflag[0]=true;
				}else{ 	
					rejflag[0]=false;
					*m_tau_prior=0;
				}
				//compute_prior()
				compute_prior(m_prior_en,m_old_prior_en,m_d_prior,m_S0_prior,fm_prior_en,m_f0_prior,m_tau_prior,m_dstd_prior,m_R_prior,nfib);
			}
			__syncthreads();
			//compute_likelihood()
			compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,&signals[*idVOX*nfib*ndirections],&isosignals[*idVOX*ndirections],&datam[*idVOX*ndirections],fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);
			__syncthreads();

			if (leader){
				rejflag[1]=compute_test_energy(m_energy,m_old_energy,*m_prior_en,*m_likelihood_en,prerandU[*localrand]);
				*localrand=*localrand+1;
				if(!rejflag[0]){
					if(rejflag[1]){
		  				*m_tau_acc=*m_tau_acc+1;   
					}else{ 
						reject(m_tau,m_tau_prior,old,m_prior_en,m_old_prior_en,m_energy,m_old_energy,m_tau_rej);
					}
				}else{ 
					reject(m_tau,m_tau_prior,old,m_prior_en,m_old_prior_en,m_energy,m_old_energy,m_tau_rej);
				}
			}	
		}
////////////////////////////////////////////////////////////////// D
		if (leader){
			propose(m_d,old,*m_d_prop,prerandN[*localrand]);
			//compute_d_prior()      
			old[1]=*m_d_prior;	
			if(*m_d<0 || *m_d>UPPERDIFF){
				rejflag[0]=true;
			}else{
				if (model==3){
					//float alpha=3.0; float beta=4000;  //Gamma_prior around 0.5-1E-3
					*m_d_prior=(1.0f-3.0f)*log(*m_d)+4000.0f**m_d;
				}else{
					*m_d_prior=0;
				}
				rejflag[0]=false;
			}
		}
		__syncthreads();	
		//compute_signal()
		for(int f=0;f<nfib;f++){
			for(int i=0; i<mydirs; i++){
				int pos = *idVOX*ndirections*nfib + f*ndirections + idSubVOX + i*threadsBlock;
				compute_signal(&signals[pos],&oldsignals[pos],bvals[*posBV+idSubVOX+i*threadsBlock],m_d,m_dstd,m_R,angtmp[pos],model);
			}
		}
		//compute_iso_signal()
		for(int i=0; i<mydirs; i++){
			int pos = *idVOX*ndirections + idSubVOX + i*threadsBlock;
			compute_iso_signal(&isosignals[pos],&oldisosignals[pos],bvals[*posBV+idSubVOX+i*threadsBlock],m_d,m_dstd,model);
		}				

		if (leader){
			//compute_prior()
			compute_prior(m_prior_en,m_old_prior_en,m_d_prior,m_S0_prior,fm_prior_en,m_f0_prior,m_tau_prior,m_dstd_prior,m_R_prior,nfib);
		}
		__syncthreads();	
		//compute_likelihood()
		compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,&signals[*idVOX*nfib*ndirections],&isosignals[*idVOX*ndirections],&datam[*idVOX*ndirections],fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);		
		__syncthreads();
				
		if (leader){
			rejflag[1]=compute_test_energy(m_energy,m_old_energy,*m_prior_en,*m_likelihood_en,prerandU[*localrand]);
			*localrand=*localrand+1;
		}	
		__syncthreads();

       		if(!rejflag[0]){
			if(rejflag[1]){
	  			if (leader) *m_d_acc=*m_d_acc+1;   
			}else{
				if (leader){
					reject(m_d,m_d_prior,old,m_prior_en,m_old_prior_en,m_energy,m_old_energy,m_d_rej);
				}
				restore_signals(signals,oldsignals,*idVOX,idSubVOX,mydirs,nfib,ndirections,threadsBlock);
				restore_isosignals(isosignals,oldisosignals,*idVOX,idSubVOX,mydirs,ndirections,threadsBlock);
			}
        	}else{ 
			if (leader){
				reject(m_d,m_d_prior,old,m_prior_en,m_old_prior_en,m_energy,m_old_energy,m_d_rej);
			}
      			restore_signals(signals,oldsignals,*idVOX,idSubVOX,mydirs,nfib,ndirections,threadsBlock);
			restore_isosignals(isosignals,oldisosignals,*idVOX,idSubVOX,mydirs,ndirections,threadsBlock);
        	}
////////////////////////////////////////////////////////////////// D_STD
		if(model>=2){
			if (leader){	
				propose(m_dstd,old,*m_dstd_prop,prerandN[*localrand]);
				//compute_d_std_prior()     
				old[1]=*m_dstd_prior;
				float upper_d_std=0.01;
				if (model==3) upper_d_std=0.004;
				if(*m_dstd<=0 || *m_dstd>upper_d_std){
					rejflag[0]=true;
				}else{
					*m_dstd_prior=log(*m_dstd);
					rejflag[0]=false;	
				}
			}
			__syncthreads();
			//compute_signal()
			if(model==2){
				for(int f=0;f<nfib;f++){
					for(int i=0; i<mydirs; i++){				
						int pos = *idVOX*ndirections*nfib + f*ndirections + idSubVOX + i*threadsBlock;
						compute_signal(&signals[pos],&oldsignals[pos],bvals[*posBV+idSubVOX+i*threadsBlock],m_d,m_dstd,m_R,angtmp[pos],model);
					}
				}
			}
			//compute_iso_signal()
			for(int i=0; i<mydirs; i++){
				int pos = *idVOX*ndirections + idSubVOX + i*threadsBlock;
				compute_iso_signal(&isosignals[pos],&oldisosignals[pos],bvals[*posBV+idSubVOX+i*threadsBlock],m_d,m_dstd,model);
			}
			if (leader){
				//compute_prior()
				compute_prior(m_prior_en,m_old_prior_en,m_d_prior,m_S0_prior,fm_prior_en,m_f0_prior,m_tau_prior,m_dstd_prior,m_R_prior,nfib);
			}
			__syncthreads();
			//compute_likelihood()
			compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,&signals[*idVOX*nfib*ndirections],&isosignals[*idVOX*ndirections],&datam[*idVOX*ndirections],fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);
			__syncthreads();

			if (leader){
				rejflag[1]=compute_test_energy(m_energy,m_old_energy,*m_prior_en,*m_likelihood_en,prerandU[*localrand]);
				*localrand=*localrand+1;					
			}
			__syncthreads();
				
			if(!rejflag[0]){
				if(rejflag[1]){
		  			if (leader) *m_dstd_acc=*m_dstd_acc+1;   
				}else{ 
					if (leader){
						reject(m_dstd,m_dstd_prior,old,m_prior_en,m_old_prior_en,m_energy,m_old_energy,m_dstd_rej);
					}
					if(model==2){
						restore_signals(signals,oldsignals,*idVOX,idSubVOX,mydirs,nfib,ndirections,threadsBlock);
					}
					restore_isosignals(isosignals,oldisosignals,*idVOX,idSubVOX,mydirs,ndirections,threadsBlock);
				}
			}else{ 
				if (leader){
					reject(m_dstd,m_dstd_prior,old,m_prior_en,m_old_prior_en,m_energy,m_old_energy,m_dstd_rej);
				}
				if(model==2){
					restore_signals(signals,oldsignals,*idVOX,idSubVOX,mydirs,nfib,ndirections,threadsBlock);
				}
				restore_isosignals(isosignals,oldisosignals,*idVOX,idSubVOX,mydirs,ndirections,threadsBlock);
			}
////////////////////////////////////////////////////////////////// R
			if(model==3){
				if (leader){	
					propose(m_R,old,*m_R_prop,prerandN[*localrand]);
					//compute_R_prior()     
					old[1]=*m_R_prior;
					float upper_R=2.0f*R_priormean;
					float lower_R=R_priormean-2.0f*R_priorstd;
      					if (R_priormean>0.5f){
						upper_R=1.0f;
					}
					if (lower_R<0.0f)
						lower_R=1e-8f;
					if (R_priorfudge>0.0f && *m_d>UPPERDIFF/2.0f){
					//then use an ARD prior to avoid competition with the isotropic compartments
 	        				if (*m_R<1E-8f || *m_R>upper_R)
 	          					rejflag[0]=true;
 	        				else{
 	          					*m_R_prior=R_priorfudge*log(*m_R);
 	          					rejflag[0]=false;
 	        				}
 	      				}else{
	      					if(*m_R<=lower_R || *m_R>upper_R){  
							//Truncate prior to avoid too spherical (high m_R) or too anisitropic (small m_R) profiles 
							rejflag[0]=true;
						}else{
							float Rstd2=R_priorstd*R_priorstd; 
							*m_R_prior=(*m_R-R_priormean)*(*m_R-R_priormean)/Rstd2;  //Gaussian prior
							rejflag[0]=false;
	      					}
					}
				}
				__syncthreads();
				//compute_signal()
				for(int f=0;f<nfib;f++){
					for(int i=0; i<mydirs; i++){				
						int pos = *idVOX*ndirections*nfib + f*ndirections + idSubVOX + i*threadsBlock;
						compute_signal(&signals[pos],&oldsignals[pos],bvals[*posBV+idSubVOX+i*threadsBlock],m_d,m_dstd,m_R,angtmp[pos],model);
					}
				}
				if (leader){
					//compute_prior()
					compute_prior(m_prior_en,m_old_prior_en,m_d_prior,m_S0_prior,fm_prior_en,m_f0_prior,m_tau_prior,m_dstd_prior,m_R_prior,nfib);
				}
				__syncthreads();
				//compute_likelihood()
				compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,&signals[*idVOX*nfib*ndirections],&isosignals[*idVOX*ndirections],&datam[*idVOX*ndirections],fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);
				__syncthreads();

				if (leader){
					rejflag[1]=compute_test_energy(m_energy,m_old_energy,*m_prior_en,*m_likelihood_en,prerandU[*localrand]);
					*localrand=*localrand+1;					
				}
				__syncthreads();
				
				if(!rejflag[0]){
					if(rejflag[1]){
		  				if (leader) *m_R_acc=*m_R_acc+1;   
					}else{ 
						if (leader){
							reject(m_R,m_R_prior,old,m_prior_en,m_old_prior_en,m_energy,m_old_energy,m_R_rej);
						}
						restore_signals(signals,oldsignals,*idVOX,idSubVOX,mydirs,nfib,ndirections,threadsBlock);
					}
				}else{ 
					if (leader){
						reject(m_R,m_R_prior,old,m_prior_en,m_old_prior_en,m_energy,m_old_energy,m_R_rej);
					}
					restore_signals(signals,oldsignals,*idVOX,idSubVOX,mydirs,nfib,ndirections,threadsBlock);
				}
			}
		}
////////////////////////////////////////////////////////////////// S0
		if (leader){
			propose(m_S0,old,*m_S0_prop,prerandN[*localrand]);
			//compute_S0_prior()
			old[1]=*m_S0_prior;
        		if(*m_S0<0) rejflag[0]=true;
        		else{    
				*m_S0_prior=0;
	  			rejflag[0]=false;
        		}
			//compute_prior()
			compute_prior(m_prior_en,m_old_prior_en,m_d_prior,m_S0_prior,fm_prior_en,m_f0_prior,m_tau_prior,m_dstd_prior,m_R_prior,nfib);
		}
		__syncthreads();
		//compute_likelihood()
		compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,&signals[*idVOX*nfib*ndirections],&isosignals[*idVOX*ndirections],&datam[*idVOX*ndirections],fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);
		__syncthreads();

		if (leader){
			rejflag[1]=compute_test_energy(m_energy,m_old_energy,*m_prior_en,*m_likelihood_en,prerandU[*localrand]);
			*localrand=*localrand+1;

        		if(!rejflag[0]){
				if(rejflag[1]){
	  				*m_S0_acc=*m_S0_acc+1;   
				}else{
					reject(m_S0,m_S0_prior,old,m_prior_en,m_old_prior_en,m_energy,m_old_energy,m_S0_rej);
				}
        		}else{ 
				reject(m_S0,m_S0_prior,old,m_prior_en,m_old_prior_en,m_energy,m_old_energy,m_S0_rej);
			}
        	}
////////////////////////////////////////////////////////////////////////////     TH
     		for(int fibre=0;fibre<nfib;fibre++){  
			if (leader){ 
				propose(&m_th[fibre],old,m_th_prop[fibre],prerandN[*localrand]);
				//compute_th_prior()
				old[1]=m_th_prior[fibre];
      	   			if(m_th[fibre]==0){
					m_th_prior[fibre]=0;
		   		}else{
					m_th_prior[fibre]=-log(double(fabs(sin(double(m_th[fibre]))/2)));
	      	   		}
		  		//rejflag[0]=false; /////////////////always false
				//compute_prior()
				*fm_old_prior_en=fm_prior_en[fibre];
	      	   		fm_prior_en[fibre]=m_th_prior[fibre]+m_ph_prior[fibre]+m_f_prior[fibre];	
			}
			__syncthreads();
			//compute_signal()
			//compute_signal_pre	
			for(int i=0; i<mydirs; i++){
				double myalpha = alpha[*posBV+idSubVOX+i*threadsBlock];
				double cos_alpha_minus_theta=cos(double(myalpha-m_th[fibre]));   
				double cos_alpha_plus_theta=cos(double(myalpha+m_th[fibre]));
				int pos = *idVOX*ndirections*nfib + fibre*ndirections + idSubVOX + i*threadsBlock;
				int pos2 = *idVOX*ndirections + idSubVOX + i*threadsBlock;
				oldangtmp[pos2]=angtmp[pos];
				double aux = (cos(double(m_ph[fibre]-beta[*posBV+idSubVOX+i*threadsBlock]))*(cos_alpha_minus_theta-cos_alpha_plus_theta)/2)+(cos_alpha_minus_theta+cos_alpha_plus_theta)/2;
		     		aux =  aux*aux;
		 		angtmp[pos]= aux;

				compute_signal(&signals[pos],&oldsignals[pos],bvals[*posBV+idSubVOX+i*threadsBlock],m_d,m_dstd,m_R,angtmp[pos],model);
			}
			if (leader){
				//compute_prior()
				compute_prior(m_prior_en,m_old_prior_en,m_d_prior,m_S0_prior,fm_prior_en,m_f0_prior,m_tau_prior,m_dstd_prior,m_R_prior,nfib);	
			}
			__syncthreads();
			//compute_likelihood()
			compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,&signals[*idVOX*nfib*ndirections],&isosignals[*idVOX*ndirections],&datam[*idVOX*ndirections],fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);
			__syncthreads();

			if (leader){ 
				rejflag[1]=compute_test_energy(m_energy,m_old_energy,*m_prior_en,*m_likelihood_en,prerandU[*localrand]);
				*localrand=*localrand+1;	
			}
			__syncthreads();
			
			if(rejflag[1]){
		  		if (leader) m_th_acc[fibre]++;   
			}else{
				if (leader){
					rejectF(&m_th[fibre],&m_th_prior[fibre],old,m_prior_en,m_old_prior_en,&fm_prior_en[fibre],fm_old_prior_en,m_energy,m_old_energy,&m_th_rej[fibre]);
				}
				//compute_signal_pre undo
				restore_angtmp_signals(signals,oldsignals,angtmp,oldangtmp,*idVOX,idSubVOX,mydirs,nfib,fibre,ndirections,threadsBlock);
			}
			__syncthreads();
///////////////////////////////////////     PH
			if (leader){
				propose(&m_ph[fibre],old,m_ph_prop[fibre],prerandN[*localrand]);
				//compute_ph_prior()
				old[1]=m_ph_prior[fibre];
      				m_ph_prior[fibre]=0;
      				//rejflag[0]=false;
				//compute_prior()
				*fm_old_prior_en=fm_prior_en[fibre];
      	   			fm_prior_en[fibre]=m_th_prior[fibre]+m_ph_prior[fibre]+m_f_prior[fibre];
			}
			__syncthreads();
			//compute_signal()
			//compute_signal_pre
			for(int i=0; i<mydirs; i++){
				double myalpha = alpha[*posBV+idSubVOX+i*threadsBlock];
				double cos_alpha_minus_theta=cos(double(myalpha-m_th[fibre]));   
			  	double cos_alpha_plus_theta=cos(double(myalpha+m_th[fibre]));
				int pos = *idVOX*ndirections*nfib + fibre*ndirections + idSubVOX + i*threadsBlock;
				int pos2 = *idVOX*ndirections + idSubVOX + i*threadsBlock;
				oldangtmp[pos2]=angtmp[pos];
				double aux = (cos(double(m_ph[fibre]-beta[*posBV+idSubVOX+i*threadsBlock]))*(cos_alpha_minus_theta-cos_alpha_plus_theta)/2)+(cos_alpha_minus_theta+cos_alpha_plus_theta)/2;
		     		aux =  aux*aux;
		 		angtmp[pos]= aux;
				
				compute_signal(&signals[pos],&oldsignals[pos],bvals[*posBV+idSubVOX+i*threadsBlock],m_d,m_dstd,m_R,angtmp[pos],model);
			}

			if (leader){
				//compute_prior()
				compute_prior(m_prior_en,m_old_prior_en,m_d_prior,m_S0_prior,fm_prior_en,m_f0_prior,m_tau_prior,m_dstd_prior,m_R_prior,nfib);
			}
			__syncthreads();
			//compute_likelihood()
			compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,&signals[*idVOX*nfib*ndirections],&isosignals[*idVOX*ndirections],&datam[*idVOX*ndirections],fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);
			__syncthreads();

			if (leader){
				rejflag[1]=compute_test_energy(m_energy,m_old_energy,*m_prior_en,*m_likelihood_en,prerandU[*localrand]);
				*localrand=*localrand+1;
			}
			__syncthreads();

			//if(!rejflag[0]){
			if(rejflag[1]){
		  		if (leader) m_ph_acc[fibre]++;   
			}else{
				if (leader){
					rejectF(&m_ph[fibre],&m_ph_prior[fibre],old,m_prior_en,m_old_prior_en,&fm_prior_en[fibre],fm_old_prior_en,m_energy,m_old_energy,&m_ph_rej[fibre]);
				}
				//compute_signal_pre undo
				restore_angtmp_signals(signals,oldsignals,angtmp,oldangtmp,*idVOX,idSubVOX,mydirs,nfib,fibre,ndirections,threadsBlock);
			}

			__syncthreads();
////////////////////////////////////////////             F
			if (leader){
				propose(&m_f[fibre],old,m_f_prop[fibre],prerandN[*localrand]);

	     			//compute_f_prior()
	        		old[1]=m_f_prior[fibre];
				if (m_f[fibre]<=0 || m_f[fibre]>=1) rejflag[0]=true;
	        		else{
		      			if(!can_use_ard ){
		  				m_f_prior[fibre]=0;
					}else{
		  				if(m_lam_jump[fibre]){
							m_f_prior[fibre]=log(double(m_f[fibre]));
						}else{
		    					m_f_prior[fibre]=0;
		  				}
					}
					m_f_prior[fibre]=fudgevalue*m_f_prior[fibre];
					rejflag[0]=false;
	      			}
				//compute_prior()
				*fm_old_prior_en=fm_prior_en[fibre];
      	   			fm_prior_en[fibre]=m_th_prior[fibre]+m_ph_prior[fibre]+m_f_prior[fibre];
						
				getfsum(fsum,m_f,*m_f0,nfib);
				//reject_f_sum()
				rejflag[1]=(*fsum>1);
				//compute_prior()
				compute_prior(m_prior_en,m_old_prior_en,m_d_prior,m_S0_prior,fm_prior_en,m_f0_prior,m_tau_prior,m_dstd_prior,m_R_prior,nfib);	
			}

			__syncthreads();
			//compute_likelihood()
			compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,&signals[*idVOX*nfib*ndirections],&isosignals[*idVOX*ndirections],&datam[*idVOX*ndirections],fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);	
			__syncthreads();

			if (leader){
				rejflag[2]=compute_test_energy(m_energy,m_old_energy,*m_prior_en,*m_likelihood_en,prerandU[*localrand]);
				*localrand=*localrand+1;

		      		if(!rejflag[0]){
					if(!rejflag[1]){
						if(rejflag[2]){
			  				m_f_acc[fibre]++;   
						}else{
							rejectF(&m_f[fibre],&m_f_prior[fibre],old,m_prior_en,m_old_prior_en,&fm_prior_en[fibre],fm_old_prior_en,m_energy,m_old_energy,&m_f_rej[fibre]);
							getfsum(fsum,m_f,*m_f0,nfib);
						}
					}else{ 
						rejectF(&m_f[fibre],&m_f_prior[fibre],old,m_prior_en,m_old_prior_en,&fm_prior_en[fibre],fm_old_prior_en,m_energy,m_old_energy,&m_f_rej[fibre]);
						getfsum(fsum,m_f,*m_f0,nfib);
					}
				}else{
					rejectF(&m_f[fibre],&m_f_prior[fibre],old,m_prior_en,m_old_prior_en,&fm_prior_en[fibre],fm_old_prior_en,m_energy,m_old_energy,&m_f_rej[fibre]);
					getfsum(fsum,m_f,*m_f0,nfib);
				}
			}
			__syncthreads();	

        	}//end while nfib

		if((record_every)&&((*recordcount%record_every)==0)&&leader){
			rd[(*idVOX*totalrecords)+*sample-1]= *m_d;
			if(m_include_f0) rf0[(*idVOX*totalrecords)+*sample-1]= *m_f0;
			if(rician) rtau[(*idVOX*totalrecords)+*sample-1]= *m_tau;
			if(model>=2) rdstd[(*idVOX*totalrecords)+*sample-1]= *m_dstd;
			if(model==3) rR[(*idVOX*totalrecords)+*sample-1]= *m_R;	
			rs0[(*idVOX*totalrecords)+*sample-1]= *m_S0;
			for(int j=0;j<nfib;j++){
				rth[(*idVOX*totalrecords*nfib)+(j*totalrecords)+*sample-1]=m_th[j];
				rph[(*idVOX*totalrecords*nfib)+(j*totalrecords)+*sample-1]=m_ph[j];
				rf[(*idVOX*totalrecords*nfib)+(j*totalrecords)+*sample-1]=m_f[j];
			}
			*sample=*sample+1;
        	}

        	if(((*count_update%updateproposalevery)==0)&&leader){
			//m_multifibre.update_proposals();
			*m_d_prop*=sqrt(float(*m_d_acc+1)/float(*m_d_rej+1));
			*m_d_prop=min(*m_d_prop,maxfloat);

			if(rician){
				*m_tau_prop*=sqrt(float(*m_tau_acc+1)/float(*m_tau_rej+1));
				*m_tau_prop=min(*m_tau_prop,maxfloat);
				*m_tau_acc=0; 
				*m_tau_rej=0;	
			}

			if(m_include_f0){
				*m_f0_prop*=sqrt(float(*m_f0_acc+1)/float(*m_f0_rej+1));
				*m_f0_prop=min(*m_f0_prop,maxfloat);
				*m_f0_acc=0; 
				*m_f0_rej=0;	
			}	

			if(model>=2){
				*m_dstd_prop*=sqrt(float(*m_dstd_acc+1)/float(*m_dstd_rej+1));
				*m_dstd_prop=min(*m_dstd_prop,maxfloat);
				*m_dstd_acc=0; 
				*m_dstd_rej=0;	
				if(model==3){
					*m_R_prop*=sqrt(float(*m_R_acc+1)/float(*m_R_rej+1));
					*m_R_prop=min(*m_R_prop,maxfloat);
					*m_R_acc=0; 
					*m_R_rej=0;
				}
			}

			*m_S0_prop*=sqrt(float(*m_S0_acc+1)/float(*m_S0_rej+1));
			*m_S0_prop=min(*m_S0_prop,maxfloat);
			*m_d_acc=0; 
			*m_d_rej=0;
			*m_S0_acc=0; 
			*m_S0_rej=0;
			for(int f=0; f<nfib;f++){
				//m_fibres[f].update_proposals();
				m_th_prop[f]*=sqrt(float(m_th_acc[f]+1)/float(m_th_rej[f]+1));
				m_th_prop[f]=min(m_th_prop[f],maxfloat);
		      		m_ph_prop[f]*=sqrt(float(m_ph_acc[f]+1)/float(m_ph_rej[f]+1));
		      		m_ph_prop[f]=min(m_ph_prop[f],maxfloat);
		      		m_f_prop[f]*=sqrt(float(m_f_acc[f]+1)/float(m_f_rej[f]+1));
		      		m_f_prop[f]=min(m_f_prop[f],maxfloat);
			      
		      		m_th_acc[f]=0; 
		      		m_th_rej[f]=0;
		      		m_ph_acc[f]=0; 
		      		m_ph_rej[f]=0;
		      		m_f_acc[f]=0; 
		      		m_f_rej[f]=0;
			}
		}

		__syncthreads();	

        } //end while iterations

	if(leader){
		multifibres[*idVOX].m_S0=*m_S0;
		multifibres[*idVOX].m_S0_prior=*m_S0_prior;
		multifibres[*idVOX].m_S0_prop=*m_S0_prop;
		multifibres[*idVOX].m_S0_acc=*m_S0_acc;
		multifibres[*idVOX].m_S0_rej=*m_S0_rej;

		multifibres[*idVOX].m_d=*m_d;
		multifibres[*idVOX].m_d_prior=*m_d_prior;
		multifibres[*idVOX].m_d_prop=*m_d_prop;
		multifibres[*idVOX].m_d_acc=*m_d_acc;
		multifibres[*idVOX].m_d_rej=*m_d_rej;
	
		multifibres[*idVOX].m_prior_en=*m_prior_en;
		multifibres[*idVOX].m_energy=*m_energy;
		multifibres[*idVOX].m_likelihood_en=*m_likelihood_en;

		if(m_include_f0){
			multifibres[*idVOX].m_f0_prior=*m_f0_prior;
			multifibres[*idVOX].m_f0=*m_f0;
			multifibres[*idVOX].m_f0_acc=*m_f0_acc;
			multifibres[*idVOX].m_f0_rej=*m_f0_rej;
			multifibres[*idVOX].m_f0_prop=*m_f0_prop;
		}
		if(rician){
			multifibres[*idVOX].m_tau_prior=*m_tau_prior;
			multifibres[*idVOX].m_tau=*m_tau;
			multifibres[*idVOX].m_tau_acc=*m_tau_acc;
			multifibres[*idVOX].m_tau_rej=*m_tau_rej;
			multifibres[*idVOX].m_tau_prop=*m_tau_prop;
		}
		if(model>=2){
			multifibres[*idVOX].m_dstd_prior=*m_dstd_prior;
			multifibres[*idVOX].m_dstd=*m_dstd;
			multifibres[*idVOX].m_dstd_acc=*m_dstd_acc;
			multifibres[*idVOX].m_dstd_rej=*m_dstd_rej;
			multifibres[*idVOX].m_dstd_prop=*m_dstd_prop;
			if(model==3){
				multifibres[*idVOX].m_R_prior=*m_R_prior;
				multifibres[*idVOX].m_R=*m_R;
				multifibres[*idVOX].m_R_acc=*m_R_acc;
				multifibres[*idVOX].m_R_rej=*m_R_rej;
				multifibres[*idVOX].m_R_prop=*m_R_prop;
			}
		}
	}
	
	if(idSubVOX<nfib){
		int pos = (*idVOX*nfib)+idSubVOX;
	
		fibres[pos].m_th=m_th[idSubVOX];
		fibres[pos].m_ph=m_ph[idSubVOX];
		fibres[pos].m_f=m_f[idSubVOX];

		fibres[pos].m_th_acc=m_th_acc[idSubVOX];
		fibres[pos].m_th_rej=m_th_rej[idSubVOX];
		fibres[pos].m_ph_acc=m_ph_acc[idSubVOX];
		fibres[pos].m_ph_rej=m_ph_rej[idSubVOX];
		fibres[pos].m_f_acc=m_f_acc[idSubVOX];
		fibres[pos].m_f_rej=m_f_rej[idSubVOX];

		fibres[pos].m_prior_en=fm_prior_en[idSubVOX];
		fibres[pos].m_th_prior=m_th_prior[idSubVOX];
		fibres[pos].m_ph_prior=m_ph_prior[idSubVOX];
		fibres[pos].m_f_prior=m_f_prior[idSubVOX];

		fibres[pos].m_th_prop=m_th_prop[idSubVOX];
		fibres[pos].m_ph_prop=m_ph_prop[idSubVOX];
		fibres[pos].m_f_prop=m_f_prop[idSubVOX];

		fibres[pos].m_lam_jump=m_lam_jump[idSubVOX];		
	}
}

