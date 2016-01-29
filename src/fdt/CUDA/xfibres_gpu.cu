/*  xfibres_gpu.cu

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
#include "xfibresoptions.h"

#include "xfibres_gpu.cuh"
#include "diffmodels.cuh"
#include "runmcmc.h"
#include "samples.h"
#include "options.h"

#include <host_vector.h>
#include <device_vector.h> 

#include <time.h>
#include <sys/time.h>
#include "init_gpu.h"
#include <fstream>

using namespace Xfibres;

void xfibres_gpu(	//INPUT
			const Matrix			datam,
			const Matrix			bvecs,
			const Matrix			bvals,
			const Matrix	 		gradm, 
			int				idpart,
			int				idSubpart,
			string				subjdir)
{
	//write num of part in a string for log file
	char part_str[8];
	char subpart_str[8];
	char aux[8];
	sprintf(part_str,"%d",idpart);
	while(strlen(part_str)<4){
		strcpy(aux,"0");
		strcat(aux,part_str);
		strcpy(part_str,aux);
	}
	sprintf(subpart_str,"%d",idSubpart);
	while(strlen(subpart_str)<4){
		strcpy(aux,"0");
		strcat(aux,subpart_str);
		strcpy(subpart_str,aux);
	}
	string gpu_log(subjdir);		//logfile
	gpu_log += ".bedpostX/logs/logs_gpu/part_";
	gpu_log += part_str;
	gpu_log += "-subpart_";
	gpu_log += subpart_str;
	std::ofstream myfile;
	myfile.open (gpu_log.data(), ios::out | ios::app );
	myfile << "----------------------------------------------------- " << "\n"; 
   	myfile << "---------------- PART " << idpart  << " SUBPART "<< idSubpart << " ------------------- " << "\n"; 
   	myfile << "----------------------------------------------------- " << "\n"; 
	myfile.close();

	xfibresOptions& opts = xfibresOptions::getInstance();

	int nvox = datam.Ncols();
	int ndirections = datam.Nrows();
	int nfib= opts.nfibres.value(); 
	bool gradnonlin=opts.grad_file.set();

	if(nvox>0){
		thrust::host_vector<float> datam_host, bvecs_host, bvals_host, params_host;
		thrust::host_vector<double> alpha_host, beta_host;
		thrust::host_vector<float> tau_host;
		vector<ColumnVector> datam_vec;
		vector<Matrix> bvecs_vec, bvals_vec;

		///// FIT /////
		prepare_data_gpu_FIT(datam,bvecs,bvals,gradm,datam_vec, bvecs_vec, bvals_vec, datam_host, bvecs_host,  bvals_host, alpha_host, beta_host, params_host, tau_host);	

		thrust::device_vector<float> datam_gpu=datam_host;
		thrust::device_vector<float> bvecs_gpu=bvecs_host;
		thrust::device_vector<float> bvals_gpu=bvals_host;	
		thrust::device_vector<float> params_gpu=params_host;
		thrust::host_vector<int> vox_repeat;	//contains the id's of voxels repeated
		vox_repeat.resize(nvox);
		int nrepeat=0;

		fit(datam_vec,bvecs_vec,bvals_vec,datam_host,bvecs_host,bvals_host,datam_gpu,bvecs_gpu,bvals_gpu,ndirections,gpu_log,params_gpu,vox_repeat,nrepeat);

		if(opts.rician.value()){
			float R = opts.R_prior_mean.value();
			int Gamma_ball_only=0;  //That flag for diffmodels means default model2
      			if (opts.modelnum.value()==3) Gamma_ball_only=2;  //That flag for diffmodels means default model3 (with constant R)	
			calculate_tau(datam_gpu,params_gpu,bvecs_gpu,bvals_gpu,vox_repeat,nrepeat,ndirections,nfib,opts.modelnum.value(),opts.f0.value(),opts.nonlin.value(),gradnonlin,R,Gamma_ball_only,gpu_log,tau_host);
		}

		bvecs_gpu.clear();		//free bvecs_gpu
		bvecs_gpu.shrink_to_fit();
	
		//////   RUN MCMC  //////
		thrust::host_vector<double> signals_host,isosignals_host;
		thrust::host_vector<FibreGPU> fibres_host;
		thrust::host_vector<MultifibreGPU> multifibres_host;
	
		prepare_data_gpu_MCMC(nvox, ndirections, nfib, signals_host, isosignals_host, fibres_host, multifibres_host);

		thrust::device_vector<double> signals_gpu=signals_host;
		thrust::device_vector<double> isosignals_gpu=isosignals_host;
		thrust::device_vector<FibreGPU> fibres_gpu=fibres_host;
		thrust::device_vector<MultifibreGPU> multifibres_gpu=multifibres_host;
		thrust::device_vector<float> tau_gpu = tau_host;
		thrust::device_vector<double> alpha_gpu=alpha_host;
		thrust::device_vector<double> beta_gpu=beta_host;

		init_Fibres_Multifibres(datam_gpu, params_gpu, tau_gpu, bvals_gpu, alpha_gpu, beta_gpu, ndirections, gpu_log, fibres_gpu, multifibres_gpu, signals_gpu, isosignals_gpu);

		srand(opts.seed.value());  //randoms seed

		runmcmc_burnin(datam_gpu, bvals_gpu, alpha_gpu, beta_gpu, ndirections, rand(), gpu_log, fibres_gpu,multifibres_gpu, signals_gpu, isosignals_gpu);

		thrust::device_vector<float> rf0_gpu,rtau_gpu,rs0_gpu,rd_gpu,rdstd_gpu,rR_gpu,rth_gpu,rph_gpu,rf_gpu;

		prepare_data_gpu_MCMC_record(nvox,rf0_gpu,rtau_gpu,rs0_gpu,rd_gpu,rdstd_gpu,rR_gpu,rth_gpu,rph_gpu,rf_gpu);

		runmcmc_record(datam_gpu, bvals_gpu, alpha_gpu,beta_gpu, fibres_gpu, multifibres_gpu, signals_gpu, isosignals_gpu, ndirections, rand(), gpu_log, rf0_gpu, rtau_gpu, rs0_gpu, rd_gpu, rdstd_gpu, rR_gpu, rth_gpu, rph_gpu, rf_gpu);

		/////// FINISH ALL VOXELS  ///////
		record_finish_voxels(rf0_gpu,rtau_gpu,rs0_gpu,rd_gpu,rdstd_gpu,rR_gpu,rth_gpu,rph_gpu,rf_gpu,nvox,idSubpart);
	}else{
		/////// FINISH EMPTY SLICE  ///////	
		Samples samples(nvox,ndirections);
		samples.save(idSubpart);
	}
}


// Correct bvals/bvecs accounting for Gradient Nonlinearities
// ColumnVector grad_nonlin has 9 entries, corresponding to the 3 components of each of the x,y and z gradient deviation
void correct_bvals_bvecs(const Matrix& bvals,const Matrix& bvecs, const ColumnVector& grad_nonlin, Matrix& bvals_c, Matrix& bvecs_c){
  	bvals_c=bvals; bvecs_c=bvecs;
  	Matrix L(3,3);  //gradient coil tensor
  	float mag;
  	L(1,1)=grad_nonlin(1);  L(1,2)=grad_nonlin(4);  L(1,3)=grad_nonlin(7);
  	L(2,1)=grad_nonlin(2);  L(2,2)=grad_nonlin(5);  L(2,3)=grad_nonlin(8);
  	L(3,1)=grad_nonlin(3);  L(3,2)=grad_nonlin(6);  L(3,3)=grad_nonlin(9);

  	IdentityMatrix Id(3);
  
  	for (int l=1; l<=bvals.Ncols(); l++){
    		if (bvals(1,l)>0){ //do not correct b0s
     		 	bvecs_c.Column(l)=(Id+L)*bvecs.Column(l);
      			mag=sqrt(bvecs_c(1,l)*bvecs_c(1,l)+bvecs_c(2,l)*bvecs_c(2,l)+bvecs_c(3,l)*bvecs_c(3,l));
      			if (mag!=0)
				bvecs_c.Column(l)=bvecs_c.Column(l)/mag;
      			bvals_c(1,l)=mag*mag*bvals(1,l);//mag^2 as b propto |G|^2
    		}
  	}
}

//////   FIT  //////
void fit(	//INPUT
		const vector<ColumnVector> 	datam_vec, 
		const vector<Matrix> 		bvecs_vec,
		const vector<Matrix> 		bvals_vec,
		thrust::host_vector<float> 	datam_host,
		thrust::host_vector<float>	bvecs_host, 
		thrust::host_vector<float>	bvals_host,
		thrust::device_vector<float> 	datam_gpu, 
		thrust::device_vector<float>	bvecs_gpu, 
		thrust::device_vector<float>	bvals_gpu,
		int 				ndirections,
		string 				output_file,
		//OUTPUT
		thrust::device_vector<float>&	params_gpu,
		thrust::host_vector<int>&	vox_repeat,	//for get residuals with or withot f0
		int&				nrepeat)
{
	std::ofstream myfile;
	myfile.open (output_file.data(), ios::out | ios::app );
   	myfile << "------------------- FIT IN GPU ---------------------- " << "\n"; 

	struct timeval t1,t2;
   	double time;
   	gettimeofday(&t1,NULL);

	xfibresOptions& opts = xfibresOptions::getInstance();
	int nvox = datam_vec.size();
	int nfib= opts.nfibres.value();
	int nparams_fit = 2+3*opts.nfibres.value();
	if(opts.modelnum.value()>=2) nparams_fit++;
	if(opts.f0.value()) nparams_fit++;
	bool gradnonlin=opts.grad_file.set();

	if(opts.modelnum.value()==1){
		if(opts.nonlin.value()){ 
			fit_PVM_single(datam_vec,bvecs_vec,bvals_vec,datam_gpu,bvecs_gpu,bvals_gpu,ndirections,nfib,opts.f0.value(),gradnonlin,output_file,params_gpu);

			if (opts.f0.value()){
				float md,mf,f0;	
				thrust::host_vector<float> params_host;
				params_host.resize(nvox*nparams_fit);
				thrust::copy(params_gpu.begin(), params_gpu.end(), params_host.begin());	
				for(int vox=0;vox<nvox;vox++){			
					md = params_host[vox*nparams_fit+(1)];
					mf = params_host[vox*nparams_fit+(2)];
					f0 = params_host[vox*nparams_fit+(nparams_fit-1)];
					if ((opts.nfibres.value()>0 && mf<0.05) || md>0.007 || f0>0.4){		//if true we need to repeat this voxel
						vox_repeat[nrepeat]=vox;
						nrepeat++;
					}
				}
				if(nrepeat>0){
					//prepare structures for the voxels that need to be reprocessed
					vector<ColumnVector> 	datam_repeat_vec; 
					vector<Matrix> 		bvecs_repeat_vec;
					vector<Matrix> 		bvals_repeat_vec;
					thrust::host_vector<float> 	datam_repeat_host;
					thrust::host_vector<float> 	bvecs_repeat_host;	
					thrust::host_vector<float> 	bvals_repeat_host;	
					thrust::host_vector<float> 	params_repeat_host;	
								
					prepare_data_gpu_FIT_repeat(datam_host, bvecs_host, bvals_host, vox_repeat, nrepeat, ndirections, datam_repeat_vec, bvecs_repeat_vec, bvals_repeat_vec, datam_repeat_host, bvecs_repeat_host, bvals_repeat_host, params_repeat_host);

					thrust::device_vector<float> datam_repeat_gpu=datam_repeat_host;
					thrust::device_vector<float> bvecs_repeat_gpu=bvecs_repeat_host;
					thrust::device_vector<float> bvals_repeat_gpu=bvals_repeat_host;	
					thrust::device_vector<float> params_repeat_gpu=params_repeat_host;
				
		 			fit_PVM_single(datam_repeat_vec,bvecs_repeat_vec,bvals_repeat_vec,datam_repeat_gpu,bvecs_repeat_gpu,bvals_repeat_gpu,ndirections,nfib,false,gradnonlin,output_file,params_repeat_gpu);
					thrust::copy(params_repeat_gpu.begin(), params_repeat_gpu.end(), params_repeat_host.begin());	
					//mix all the parameteres: repeated and not repeated
					mix_params(params_repeat_host,vox_repeat, nrepeat, nvox, params_gpu);
				}
	  		}
		}else{
			fit_PVM_single_c(datam_vec,bvecs_vec,bvals_vec,datam_gpu,bvecs_gpu,bvals_gpu,ndirections,nfib,opts.f0.value(),gradnonlin,output_file,params_gpu);

			if (opts.f0.value()){
				float md,mf,f0;	
				thrust::host_vector<float> params_host;
				params_host.resize(nvox*nparams_fit);
				thrust::copy(params_gpu.begin(), params_gpu.end(), params_host.begin());	
				for(int vox=0;vox<nvox;vox++){		
					md = params_host[vox*nparams_fit+(1)];
					mf = params_host[vox*nparams_fit+(2)];
					f0 = params_host[vox*nparams_fit+(nparams_fit-1)];
					if ((opts.nfibres.value()>0 && mf<0.05) || md>0.007 || f0>0.4){		//if true we need to repeat this voxel
						vox_repeat[nrepeat]=vox;
						nrepeat++;
					}
				}
				if(nrepeat>0){
					//prepare structures for the voxels that need to be reprocessed
					vector<ColumnVector> 	datam_repeat_vec; 
					vector<Matrix> 		bvecs_repeat_vec;
					vector<Matrix> 		bvals_repeat_vec;
					thrust::host_vector<float> 	datam_repeat_host;
					thrust::host_vector<float> 	bvecs_repeat_host;	
					thrust::host_vector<float> 	bvals_repeat_host;	
					thrust::host_vector<float> 	params_repeat_host;	
								
					prepare_data_gpu_FIT_repeat(datam_host, bvecs_host, bvals_host, vox_repeat, nrepeat, ndirections, datam_repeat_vec, bvecs_repeat_vec, bvals_repeat_vec, datam_repeat_host, bvecs_repeat_host, bvals_repeat_host, params_repeat_host);

					thrust::device_vector<float> datam_repeat_gpu=datam_repeat_host;
					thrust::device_vector<float> bvecs_repeat_gpu=bvecs_repeat_host;
					thrust::device_vector<float> bvals_repeat_gpu=bvals_repeat_host;	
					thrust::device_vector<float> params_repeat_gpu=params_repeat_host;
				
		 			fit_PVM_single_c(datam_repeat_vec,bvecs_repeat_vec,bvals_repeat_vec,datam_repeat_gpu,bvecs_repeat_gpu,bvals_repeat_gpu,ndirections,nfib,false,gradnonlin,output_file,params_repeat_gpu);
					thrust::copy(params_repeat_gpu.begin(), params_repeat_gpu.end(), params_repeat_host.begin());	

					//mix all the parameteres: repeated and not repeated
					mix_params(params_repeat_host ,vox_repeat, nrepeat, nvox, params_gpu);
				}
	  		}
		}
	}else{
      		//model 2 or 3
		fit_PVM_single(datam_vec,bvecs_vec,bvals_vec,datam_gpu,bvecs_gpu,bvals_gpu,ndirections,nfib,opts.f0.value(),gradnonlin,output_file,params_gpu);
	
		float R = opts.R_prior_mean.value();
		int Gamma_ball_only=0;  //That flag for diffmodels means default model2
      		if (opts.modelnum.value()==3) Gamma_ball_only=2;  //That flag for diffmodels means default model3 (with constant R)

		fit_PVM_multi(datam_gpu,bvecs_gpu,bvals_gpu,nvox,ndirections,nfib,opts.f0.value(),gradnonlin,R,Gamma_ball_only,output_file,params_gpu);	

		if (opts.f0.value()){
				float md,mf,f0;	
				thrust::host_vector<float> params_host;
				params_host.resize(nvox*nparams_fit);
				thrust::copy(params_gpu.begin(), params_gpu.end(), params_host.begin());	
				for(int vox=0;vox<nvox;vox++){			
					md = params_host[vox*nparams_fit+(1)];
					mf = params_host[vox*nparams_fit+(3)];
					f0 = params_host[vox*nparams_fit+(nparams_fit-1)];
					if ((opts.nfibres.value()>0 && mf<0.05) || md>0.007 || f0>0.4){		//if true we need to repeat this voxel
						vox_repeat[nrepeat]=vox;
						nrepeat++;
					}
				}
				if(nrepeat>0){
					//prepare structures for the voxels that need to be reprocessed
					vector<ColumnVector> 	datam_repeat_vec; 
					vector<Matrix> 		bvecs_repeat_vec;
					vector<Matrix> 		bvals_repeat_vec;
					thrust::host_vector<float> 	datam_repeat_host;
					thrust::host_vector<float> 	bvecs_repeat_host;	
					thrust::host_vector<float> 	bvals_repeat_host;	
					thrust::host_vector<float> 	params_repeat_host;		
								
					prepare_data_gpu_FIT_repeat(datam_host, bvecs_host, bvals_host, vox_repeat, nrepeat, ndirections, datam_repeat_vec, bvecs_repeat_vec, bvals_repeat_vec, datam_repeat_host, bvecs_repeat_host,  bvals_repeat_host, params_repeat_host);

					thrust::device_vector<float> datam_repeat_gpu=datam_repeat_host;
					thrust::device_vector<float> bvecs_repeat_gpu=bvecs_repeat_host;
					thrust::device_vector<float> bvals_repeat_gpu=bvals_repeat_host;	
					thrust::device_vector<float> params_repeat_gpu=params_repeat_host;
				
		 			fit_PVM_single_c(datam_repeat_vec,bvecs_repeat_vec,bvals_repeat_vec,datam_repeat_gpu,bvecs_repeat_gpu,bvals_repeat_gpu,ndirections,nfib,false,gradnonlin,output_file,params_repeat_gpu);

					fit_PVM_multi(datam_repeat_gpu,bvecs_repeat_gpu,bvals_repeat_gpu,nrepeat,ndirections,nfib,false,gradnonlin,R,Gamma_ball_only,output_file,params_repeat_gpu);	
					thrust::copy(params_repeat_gpu.begin(), params_repeat_gpu.end(), params_repeat_host.begin());	
		
					//mix all the parameteres: repeated and not repeated
					mix_params(params_repeat_host ,vox_repeat, nrepeat,  nvox, params_gpu);
				}
	  		}	
	}

	gettimeofday(&t2,NULL);
    	time=timeval_diff(&t2,&t1);
   	myfile << "TIME TOTAL: " << time << " seconds\n"; 
	myfile << "-----------------------------------------------------" << "\n\n" ; 
	myfile.close();
}

void remove_NonPositive_entries(ColumnVector& Voxdata){  
//Zero, Negative Entries can be obtained from spline interpolation 
  	int pos; 
  	float MinS=Voxdata.Minimum1(pos); 
  	float MaxS=Voxdata.Maximum();
	if (MinS<=0 && MaxS>0){  
	//when there are some non-positive entries, but not all are zero
    		vector<int> minpositions;
	    	while (MinS<=0){
	      		minpositions.push_back(pos);
	      		Voxdata(pos)=MaxS;    //temporarilly make the non-positive values Max
	      		MinS=Voxdata.Minimum1(pos);
	    	}
	    	MinS=Voxdata.Minimum(); //Now find the Minimum of positive entries
	    	for (unsigned int i=0; i<minpositions.size(); i++)
	      		Voxdata(minpositions[i])=MinS; //Replace non-positive entries with that minimum
  	}
}

//prepare the structures for copy all neccesary data to FIT in GPU
void prepare_data_gpu_FIT(	//INPUT
				const Matrix				datam,
				const Matrix				bvecs,
				const Matrix				bvals,
				const Matrix	 			gradm, 
				//OUTPUT
				vector<ColumnVector>&			datam_vec,
				vector<Matrix>&				bvecs_vec,
				vector<Matrix>&				bvals_vec,
				thrust::host_vector<float>&   		datam_host,	//data prepared for copy to GPU
				thrust::host_vector<float>&		bvecs_host,				
				thrust::host_vector<float>&		bvals_host,
				thrust::host_vector<double>&		alpha_host,
				thrust::host_vector<double>&		beta_host,
				thrust::host_vector<float>&		params_host,
				thrust::host_vector<float>&		tau_host)
{
	xfibresOptions& opts = xfibresOptions::getInstance();
	int nvox = datam.Ncols(); 
	int ndirections = datam.Nrows(); 

	datam_vec.resize(nvox);
	datam_host.resize(nvox*ndirections); 
	for(int vox=0;vox<nvox;vox++){
		ColumnVector voxdata;
		voxdata=datam.Column(vox+1);
     		if(opts.rician.value()) remove_NonPositive_entries(voxdata); //So that log(data) does not give infinity in the likelihood
		datam_vec[vox]=voxdata;
		for(int j=0;j<ndirections;j++){
			datam_host[vox*ndirections+j]=voxdata(j+1);
		}
	}

	if (opts.grad_file.set()){
		bvecs_vec.resize(nvox);
		bvals_vec.resize(nvox);
		bvecs_host.resize(nvox*bvecs.Nrows()*bvecs.Ncols());
		bvals_host.resize(nvox*bvals.Ncols());
		alpha_host.resize(nvox*bvecs.Ncols());
		beta_host.resize(nvox*bvecs.Ncols());
	}else{
		bvecs_vec.resize(1);
		bvals_vec.resize(1);
		bvecs_host.resize(1*bvecs.Nrows()*bvecs.Ncols());
		bvals_host.resize(1*bvals.Ncols());
		alpha_host.resize(1*bvecs.Ncols());
		beta_host.resize(1*bvecs.Ncols());
	}

	ColumnVector alpha,beta;

	if (opts.grad_file.set()){
		for(int vox=0;vox<nvox;vox++){
			correct_bvals_bvecs(bvals,bvecs, gradm.Column(vox+1),bvals_vec[vox],bvecs_vec[vox]); //correct for gradient nonlinearities
 			MISCMATHS::cart2sph(bvecs_vec[vox],alpha,beta);
			
			for(int dir=0;dir<ndirections;dir++){
				bvecs_host[vox*ndirections*3+dir] = bvecs_vec[vox](1,dir+1);
				bvecs_host[vox*ndirections*3+ndirections+dir] = bvecs_vec[vox](2,dir+1);
				bvecs_host[vox*ndirections*3+ndirections*2+dir] = bvecs_vec[vox](3,dir+1);
				bvals_host[vox*ndirections+dir] = bvals_vec[vox](1,dir+1);

				alpha_host[vox*ndirections+dir] = alpha(dir+1);
        			beta_host[vox*ndirections+dir] = beta(dir+1);
			}
		}
		
	}else{
 		MISCMATHS::cart2sph(bvecs,alpha,beta);
		bvecs_vec[0]=bvecs;
		bvals_vec[0]=bvals;
		for(int dir=0;dir<ndirections;dir++){
			bvecs_host[dir] = bvecs(1,dir+1);
			bvecs_host[ndirections+dir] = bvecs(2,dir+1);
			bvecs_host[ndirections*2+dir] = bvecs(3,dir+1);
        		bvals_host[dir] = bvals(1,dir+1);
			
			alpha_host[dir] = alpha(dir+1);
        		beta_host[dir] = beta(dir+1);
		}
	}
	
	int nfib= opts.nfibres.value();
	int nparams;

	if(opts.f0.value()) nparams=3+nfib*3;
	else nparams=2+nfib*3;	
	if(opts.modelnum.value()>=2) nparams++;

	params_host.resize(nvox*nparams);
	tau_host.resize(nvox);
}

//prepare the structures for copy all neccesary data to FIT in GPU when is repeated because f0. Only some voxels
void prepare_data_gpu_FIT_repeat(	//INPUT
					thrust::host_vector<float>   		datam_host,	
					thrust::host_vector<float>		bvecs_host,				
					thrust::host_vector<float>		bvals_host,
					thrust::host_vector<int>		vox_repeat,
					int					nrepeat,
					int					ndirections,
					//OUTPUT
					vector<ColumnVector>&			datam_repeat_vec,
					vector<Matrix>&				bvecs_repeat_vec,
					vector<Matrix>&				bvals_repeat_vec,
					thrust::host_vector<float>&   		datam_repeat_host,	//data prepared for copy to GPU
					thrust::host_vector<float>&		bvecs_repeat_host,				
					thrust::host_vector<float>&		bvals_repeat_host,
					thrust::host_vector<float>&		params_repeat_host)
{
	xfibresOptions& opts = xfibresOptions::getInstance();

	ColumnVector datam(ndirections);
	Matrix	bvecs(3,ndirections);
	Matrix	bvals(1,ndirections);

	datam_repeat_vec.resize(nrepeat);
	datam_repeat_host.resize(nrepeat*ndirections); 
	
	if (opts.grad_file.set()){
		bvecs_repeat_vec.resize(nrepeat);
		bvals_repeat_vec.resize(nrepeat);
		bvecs_repeat_host.resize(nrepeat*3*ndirections);
		bvals_repeat_host.resize(nrepeat*ndirections);
	}else{
		bvecs_repeat_vec.resize(1);
		bvals_repeat_vec.resize(1);
		bvecs_repeat_host.resize(1*3*ndirections);
		bvals_repeat_host.resize(1*ndirections);
	}

	
	for(int vox=0;vox<nrepeat;vox++){
		for(int dir=0;dir<ndirections;dir++){
			datam(dir+1)= datam_host[vox_repeat[vox]*ndirections+dir]; 
			datam_repeat_host[vox*ndirections+dir]=datam_host[vox_repeat[vox]*ndirections+dir];
		}
		datam_repeat_vec[vox]=datam;
	}

	if (opts.grad_file.set()){
		for(int vox=0;vox<nrepeat;vox++){
			for(int dir=0;dir<ndirections;dir++){
				bvecs_repeat_host[vox*ndirections*3+dir] = bvecs_host[vox_repeat[vox]*ndirections*3+dir];
				bvecs_repeat_host[vox*ndirections*3+ndirections+dir] = bvecs_host[vox_repeat[vox]*ndirections*3+ndirections+dir];
				bvecs_repeat_host[vox*ndirections*3+ndirections*2+dir] = bvecs_host[vox_repeat[vox]*ndirections*3+ndirections*2+dir];
				bvals_repeat_host[vox*ndirections+dir] = bvals_host[vox_repeat[vox]*ndirections+dir];
			
				bvecs(1,dir+1)= bvecs_host[vox_repeat[vox]*ndirections*3+dir];
				bvecs(2,dir+1)= bvecs_host[vox_repeat[vox]*ndirections*3+ndirections+dir];
				bvecs(3,dir+1)= bvecs_host[vox_repeat[vox]*ndirections*3+ndirections*2+dir];
				bvals(1,dir+1)= bvals_host[vox_repeat[vox]*ndirections+dir];
			}
			bvecs_repeat_vec[vox]=bvecs;
			bvals_repeat_vec[vox]=bvals;
		}
	}else{
		for(int dir=0;dir<ndirections;dir++){
			bvecs_repeat_host[dir] = bvecs_host[dir];
			bvecs_repeat_host[ndirections+dir] = bvecs_host[ndirections+dir];
			bvecs_repeat_host[ndirections*2+dir] = bvecs_host[ndirections*2+dir];
			bvals_repeat_host[dir] = bvals_host[dir];
			
			bvecs(1,dir+1)= bvecs_host[dir];
			bvecs(2,dir+1)= bvecs_host[ndirections+dir];
			bvecs(3,dir+1)= bvecs_host[ndirections*2+dir];
			bvals(1,dir+1)= bvals_host[dir];
		}
		bvecs_repeat_vec[0]=bvecs;
		bvals_repeat_vec[0]=bvals;
	}
	
	int nfib= opts.nfibres.value();
	int nparams;

	nparams=2+nfib*3;	
	if(opts.modelnum.value()>=2) nparams++;

	params_repeat_host.resize(nrepeat*nparams);
}


void mix_params(	//INPUT
			thrust::host_vector<float>   		params_repeat_host,
			thrust::host_vector<int>		vox_repeat,
			int					nrepeat,
			int					nvox,
			//INPUT-OUTPUT
			thrust::device_vector<float>&   	params_gpu)
{
	xfibresOptions& opts = xfibresOptions::getInstance();
	int nfib= opts.nfibres.value();
	int nparams = 2+3*opts.nfibres.value();
	if(opts.modelnum.value()>=2) nparams++;

	thrust::host_vector<float> params_host;
	params_host.resize(nvox*(nparams+1));
	thrust::copy(params_gpu.begin(), params_gpu.end(), params_host.begin());	

	for(int vox=0;vox<nrepeat;vox++){
		for(int par=0;par<nparams;par++){
			params_host[vox_repeat[vox]*(nparams+1)+par] = params_repeat_host[vox*nparams+par]; //(nparams+1) to count f0
		}
		params_host[vox_repeat[vox]*(nparams+1)+nparams] = 0.001;	//pvmf0=0.001
	}
	thrust::copy(params_host.begin(), params_host.end(), params_gpu.begin());	
}

void prepare_data_gpu_MCMC(	//INPUT
				int 					nvox,
				int					ndirections,
				int 					nfib,
				//OUTPUT
				thrust::host_vector<double>&		signals_host,
				thrust::host_vector<double>&		isosignals_host,
				thrust::host_vector<FibreGPU>& 		fibres_host,
				thrust::host_vector<MultifibreGPU>& 	multifibres_host)
{ 	
	signals_host.resize(nvox*nfib*ndirections);
	isosignals_host.resize(nvox*ndirections);	
	fibres_host.resize(nvox*nfib);	
	multifibres_host.resize(nvox);
}

void prepare_data_gpu_MCMC_record(	//INPUT
					int 						nvox,
					//OUTPUT
					thrust::device_vector<float>&			rf0_gpu,
					thrust::device_vector<float>&			rtau_gpu,
					thrust::device_vector<float>&			rs0_gpu,
					thrust::device_vector<float>&			rd_gpu,
					thrust::device_vector<float>&			rdstd_gpu,
					thrust::device_vector<float>&			rR_gpu,
					thrust::device_vector<float>&			rth_gpu,
					thrust::device_vector<float>&			rph_gpu,
					thrust::device_vector<float>&			rf_gpu)
{ 	
	xfibresOptions& opts = xfibresOptions::getInstance();

	int nfib = opts.nfibres.value();	
	int nsamples = (opts.njumps.value()/opts.sampleevery.value());   
	
	if(opts.f0.value()) rf0_gpu.resize(nvox*nsamples); 
	if(opts.rician.value()) rtau_gpu.resize(nvox*nsamples);  
	rs0_gpu.resize(nvox*nsamples);  
	rd_gpu.resize(nvox*nsamples);
	if(opts.modelnum.value()>=2) rdstd_gpu.resize(nvox*nsamples);  
	if(opts.modelnum.value()==3) rR_gpu.resize(nvox*nsamples); 
	rth_gpu.resize(nvox*nsamples*nfib);  
	rph_gpu.resize(nvox*nsamples*nfib);  
	rf_gpu.resize(nvox*nsamples*nfib);  
}

void record_finish_voxels(	//INPUT
				thrust::device_vector<float>&			rf0_gpu,
				thrust::device_vector<float>&			rtau_gpu,
				thrust::device_vector<float>&			rs0_gpu,
				thrust::device_vector<float>&			rd_gpu,
				thrust::device_vector<float>&			rdstd_gpu,
				thrust::device_vector<float>&			rR_gpu,
				thrust::device_vector<float>&			rth_gpu,
				thrust::device_vector<float>&			rph_gpu,
				thrust::device_vector<float>&			rf_gpu,
				int 						nvox,
				int						idpart)
{
	xfibresOptions& opts = xfibresOptions::getInstance();

	int nfib = opts.nfibres.value();	
	int nsamples = (opts.njumps.value()/opts.sampleevery.value());   

	thrust::host_vector<float> rf0_host,rtau_host,rs0_host,rd_host,rdstd_host,rR_host,rth_host,rph_host,rf_host;

	rf0_host.resize(nvox*nsamples);
	rtau_host.resize(nvox*nsamples);
	rs0_host.resize(nvox*nsamples);
	rd_host.resize(nvox*nsamples);
	if(opts.modelnum.value()>=2) rdstd_host.resize(nvox*nsamples);
	if(opts.modelnum.value()==3) rR_host.resize(nvox*nsamples);
	rth_host.resize(nvox*nfib*nsamples);
	rph_host.resize(nvox*nfib*nsamples);
	rf_host.resize(nvox*nfib*nsamples);

	if(opts.f0.value()) thrust::copy(rf0_gpu.begin(), rf0_gpu.end(), rf0_host.begin());
	if(opts.rician.value()) thrust::copy(rtau_gpu.begin(), rtau_gpu.end(), rtau_host.begin());
	thrust::copy(rs0_gpu.begin(), rs0_gpu.end(), rs0_host.begin());
	thrust::copy(rd_gpu.begin(), rd_gpu.end(), rd_host.begin());
	if(opts.modelnum.value()>=2) thrust::copy(rdstd_gpu.begin(), rdstd_gpu.end(), rdstd_host.begin());
	if(opts.modelnum.value()==3) thrust::copy(rR_gpu.begin(), rR_gpu.end(), rR_host.begin());
	thrust::copy(rth_gpu.begin(), rth_gpu.end(), rth_host.begin());
	thrust::copy(rph_gpu.begin(), rph_gpu.end(), rph_host.begin());
	thrust::copy(rf_gpu.begin(), rf_gpu.end(), rf_host.begin());	

	Samples samples(nvox,nsamples);

	float ard,arf0,artau,ardstd,arR,ars0;	
	float *arth = new float[nfib];
    	float *arph = new float[nfib]; 
    	float *arf = new float[nfib];

	for(int vox=0;vox<nvox;vox++){
		for(int rec=0;rec<nsamples;rec++){	
			ard=rd_host[(vox*nsamples)+rec];
			if(opts.f0.value()){	
				arf0=rf0_host[(vox*nsamples)+rec];
			}

			if(opts.rician.value()){	
				artau=rtau_host[(vox*nsamples)+rec];
			}

			if(opts.modelnum.value()>=2){	
				ardstd=rdstd_host[(vox*nsamples)+rec];
			}
			if(opts.modelnum.value()==3){	
				arR=rR_host[(vox*nsamples)+rec];
			}
		
			ars0=rs0_host[(vox*nsamples)+rec];

			for(int j=0;j<nfib;j++){
				arth[j]=rth_host[(vox*nfib*nsamples)+(j*nsamples)+rec];
				arph[j]=rph_host[(vox*nfib*nsamples)+(j*nsamples)+rec];
				arf[j]=rf_host[(vox*nfib*nsamples)+(j*nsamples)+rec];

			}
			samples.record(ard,arf0,artau,ardstd,arR,ars0,arth,arph,arf,vox+1,rec+1);
		}	
		samples.finish_voxel(vox+1);
   	}
	samples.save(idpart);
}

