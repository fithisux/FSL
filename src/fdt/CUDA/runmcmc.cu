/*  runmcmc.cu

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

#include "xfibresoptions.h"
#include <curand.h>
#include "runmcmc_kernels.cu"
#include "sync_check.h"

#include <host_vector.h>
#include <device_vector.h> 

#include <time.h>
#include <sys/time.h>
#include "init_gpu.h"

using namespace Xfibres;

////////////////////////////////////////////////////// 
//   MCMC IN GPU
////////////////////////////////////////////////////// 

void init_Fibres_Multifibres(	//INPUT
				thrust::device_vector<float> 			datam_gpu,
				thrust::device_vector<float> 			params_gpu,
				thrust::device_vector<float> 			tau_gpu,
				thrust::device_vector<float> 			bvals_gpu,
				thrust::device_vector<double> 			alpha_gpu,
				thrust::device_vector<double> 			beta_gpu,
				const int 					ndirections,
				string 						output_file, 
				//OUTPUT
				thrust::device_vector<FibreGPU>& 		fibres_gpu,
				thrust::device_vector<MultifibreGPU>& 		multifibres_gpu,
				thrust::device_vector<double>&			signals_gpu,
				thrust::device_vector<double>&			isosignals_gpu)
{
	std::ofstream myfile;
	myfile.open (output_file.data(), ios::out | ios::app );
   	myfile << "----- MCMC ALGORITHM PART INITIALITATION IN GPU ----- " << "\n";  	

   	struct timeval t1,t2;
   	double time;
   	gettimeofday(&t1,NULL);

	int nvox = multifibres_gpu.size();

	xfibresOptions& opts = xfibresOptions::getInstance();
	int nfib= opts.nfibres.value();
	int nparams_fit = 2+3*opts.nfibres.value();
	if(opts.modelnum.value()>=2) nparams_fit++;
	if(opts.f0.value()) nparams_fit++;

	thrust::device_vector<double> angtmp_gpu;
	angtmp_gpu.resize(nvox*ndirections*nfib);
	

	bool gradnonlin = opts.grad_file.set();

	int blocks = nvox; 
  	dim3 Dim_Grid_MCMC(blocks, 1);
  	dim3 Dim_Block_MCMC(THREADS_BLOCK_MCMC ,1);	///dimensions for MCMC

	float *datam_ptr = thrust::raw_pointer_cast(datam_gpu.data());
	float *params_ptr = thrust::raw_pointer_cast(params_gpu.data());	
	float *tau_ptr = thrust::raw_pointer_cast(tau_gpu.data());	
	float *bvals_ptr = thrust::raw_pointer_cast(bvals_gpu.data());
	double *alpha_ptr = thrust::raw_pointer_cast(alpha_gpu.data());
	double *beta_ptr = thrust::raw_pointer_cast(beta_gpu.data());
	FibreGPU *fibres_ptr =  thrust::raw_pointer_cast(fibres_gpu.data());
	MultifibreGPU *multifibres_ptr = thrust::raw_pointer_cast(multifibres_gpu.data());
	double *signals_ptr = thrust::raw_pointer_cast(signals_gpu.data());
	double *isosignals_ptr = thrust::raw_pointer_cast(isosignals_gpu.data());
	double *angtmp_ptr = thrust::raw_pointer_cast(angtmp_gpu.data());

	int amount_shared = (THREADS_BLOCK_MCMC)*sizeof(double) + (3*nfib + 9)*sizeof(float) + sizeof(int);

	myfile << "Shared Memory Used in init_Fibres_Multifibres: " << amount_shared << "\n";

	init_Fibres_Multifibres_kernel<<< Dim_Grid_MCMC, Dim_Block_MCMC, amount_shared>>>(datam_ptr, params_ptr, tau_ptr, bvals_ptr, alpha_ptr, beta_ptr, opts.R_prior_mean.value(), opts.R_prior_std.value(),opts.R_prior_fudge.value(), ndirections, nfib, nparams_fit, opts.modelnum.value(), opts.fudge.value(), opts.f0.value(), opts.rician.value(), opts.ardf0.value(), opts.all_ard.value(), opts.no_ard.value(), gradnonlin, angtmp_ptr, fibres_ptr, multifibres_ptr, signals_ptr, isosignals_ptr);
	sync_check("init_Fibres_Multifibres_kernel");

	gettimeofday(&t2,NULL);
    	time=timeval_diff(&t2,&t1);
   	myfile << "TIME TOTAL: " << time << " seconds\n"; 
	myfile << "-----------------------------------------------------" << "\n\n" ; 
	myfile.close();
}

void runmcmc_burnin(	//INPUT
			thrust::device_vector<float> 			datam_gpu,
			thrust::device_vector<float> 			bvals_gpu,
			thrust::device_vector<double> 			alpha_gpu,
			thrust::device_vector<double> 			beta_gpu,
			const int 					ndirections,
			double 						seed,
			string 						output_file, 
			//INPUT-OUTPUT
			thrust::device_vector<FibreGPU>& 		fibres_gpu,
			thrust::device_vector<MultifibreGPU>& 		multifibres_gpu,
			thrust::device_vector<double>&			signals_gpu,
			thrust::device_vector<double>&			isosignals_gpu)
{
	xfibresOptions& opts = xfibresOptions::getInstance();
	
	std::ofstream myfile;
	myfile.open (output_file.data(), ios::out | ios::app ); 
   	myfile << "--------- MCMC ALGORITHM PART BURNIN IN GPU --------- " << "\n";  	

   	struct timeval t1,t2,t_tot1,t_tot2;
   	double time,timecurand,timemcmc;
   	time=0;
   	timecurand=0;
   	timemcmc=0;

   	gettimeofday(&t_tot1,NULL);

   	size_t free,total;
	
	int nvox = multifibres_gpu.size();
   	int nfib= opts.nfibres.value();
	int nparams;

	bool gradnonlin=opts.grad_file.set();

	if(opts.f0.value()) nparams=3+nfib*3;
	else nparams=2+nfib*3;	
	if(opts.modelnum.value()>=2) nparams++;
	if(opts.modelnum.value()==3) nparams++;	
	if(opts.rician.value()) nparams++;

	thrust::device_vector<float> recors_null_gpu;
	recors_null_gpu.resize(1);

	thrust::device_vector<double> angtmp_gpu;
	thrust::device_vector<double> oldangtmp_gpu;
	thrust::device_vector<double> oldsignals_gpu;
	thrust::device_vector<double> oldisosignals_gpu;
	
	angtmp_gpu.resize(nvox*ndirections*nfib);
	oldangtmp_gpu.resize(nvox*ndirections);
	oldsignals_gpu.resize(nvox*ndirections*nfib);
	oldisosignals_gpu.resize(nvox*ndirections);
   
   	unsigned int totalrandoms=(opts.nburn.value() * nvox * nparams);

   	cuMemGetInfo(&free,&total);
   	myfile << "Free memory Before Randoms: "<< free <<  " ---- Total memory: " << total << "\n";
   	//4 bytes each float, 2 random arrays, and 80% of total memory at this moment 
   	unsigned int maxrandoms=((free*0.8)/(4*2)); 

   	myfile << "Total randoms: " << totalrandoms << "\n"; 
   	myfile << "Max randoms: " << maxrandoms << "\n"; 
   
   	int steps; //num iter if not enough memory
   	int minrandoms; //min num of randoms ensamble
   	minrandoms= nvox * nparams;

   	int iters_step=0;
	int nrandoms=0;	

   	if(totalrandoms>maxrandoms){ 
		iters_step = maxrandoms / minrandoms; 			//iterations in each step
		nrandoms = iters_step*minrandoms;			//nrandoms for each step
		steps =  (opts.nburn.value()/iters_step);  		//repeat process steps times, no enough memory for all randoms 			
   	}else{ 
		nrandoms = totalrandoms;
		iters_step= opts.nburn.value();
		steps = 0;
  	}
	if(nrandoms%2){							//CURAND must generates multiples of 2 randoms
		nrandoms++;
	}
	
   	myfile << "Process " << opts.nburn.value() << " iterations divided in "<< steps << " steps with "<< iters_step << " iterations in each one" << "\n";    

   	int last_step = opts.nburn.value() - (iters_step*steps);
   	int last_randoms = (last_step*minrandoms);
	if(last_randoms%2){						//CURAND must generates multiples of 2 randoms
		last_randoms++;
	}

   	myfile << "Last step with " << last_step << " iterations" << "\n"; 
	
	thrust::device_vector<float> randomsN_gpu;
	thrust::device_vector<float> randomsU_gpu;	
	randomsN_gpu.resize(nrandoms);
	randomsU_gpu.resize(nrandoms);

   	cuMemGetInfo(&free,&total);
   	myfile << "Free memory after Malloc Randoms: "<< free <<  " ---- Total memory: " << total << "\n";
   
  	int blocks = nvox;        
  	dim3 Dim_Grid(blocks, 1);
  	dim3 Dim_Block(THREADS_BLOCK_MCMC,1);	//dimensions for MCMC   

   	myfile << "\n" << "NUM BLOCKS: " << blocks << "\n"; 
   	myfile << "THREADS PER BLOCK : " << THREADS_BLOCK_MCMC << "\n\n"; 	

   	curandGenerator_t gen;
   	curandCreateGenerator(&gen,CURAND_RNG_PSEUDO_DEFAULT);
   	curandSetPseudoRandomGeneratorSeed(gen,seed);

	//get pointers
	float *datam_ptr = thrust::raw_pointer_cast(datam_gpu.data());
	float *bvals_ptr = thrust::raw_pointer_cast(bvals_gpu.data());
	double *alpha_ptr = thrust::raw_pointer_cast(alpha_gpu.data());
	double *beta_ptr = thrust::raw_pointer_cast(beta_gpu.data());
	float *randomsN_ptr = thrust::raw_pointer_cast(randomsN_gpu.data());
	float *randomsU_ptr = thrust::raw_pointer_cast(randomsU_gpu.data());
	FibreGPU *fibres_ptr =  thrust::raw_pointer_cast(fibres_gpu.data());
	MultifibreGPU *multifibres_ptr = thrust::raw_pointer_cast(multifibres_gpu.data());
	double *signals_ptr = thrust::raw_pointer_cast(signals_gpu.data());
	double *isosignals_ptr = thrust::raw_pointer_cast(isosignals_gpu.data());

	double *angtmp_ptr = thrust::raw_pointer_cast(angtmp_gpu.data());
	double *oldangtmp_ptr = thrust::raw_pointer_cast(oldangtmp_gpu.data());
	double *oldsignals_ptr = thrust::raw_pointer_cast(oldsignals_gpu.data());
	double *oldisosignals_ptr = thrust::raw_pointer_cast(oldisosignals_gpu.data());

	float *records_null = thrust::raw_pointer_cast(recors_null_gpu.data());

	int amount_shared = (THREADS_BLOCK_MCMC)*sizeof(double) + (10*nfib + 2*nparams + 27)*sizeof(float) + (7*nfib + 21)*sizeof(int);

	myfile << "Shared Memory Used in runmcmc_burnin: " << amount_shared << "\n";
	
   	for(int i=0;i<steps;i++){

   		gettimeofday(&t1,NULL);

	   	curandStatus_t status = curandGenerateNormal(gen,randomsN_ptr,nrandoms,0,1); 
		if (status != CURAND_STATUS_SUCCESS)
		{
			printf("Failure generating cuda random numbers: %d\n",status);
			exit(1);
		}
	   	status = curandGenerateUniform(gen,randomsU_ptr,nrandoms);	//generate randoms
		if (status != CURAND_STATUS_SUCCESS)
		{
			printf("Failure generating cuda random numbers: %d\n",status);
			exit(1);
		}

 	   	gettimeofday(&t2,NULL);
    	   	timecurand+=timeval_diff(&t2,&t1);

	   	gettimeofday(&t1,NULL);

	   	runmcmc_kernel<<< Dim_Grid, Dim_Block, amount_shared >>>(datam_ptr, bvals_ptr, alpha_ptr, beta_ptr, randomsN_ptr, randomsU_ptr, opts.R_prior_mean.value(), opts.R_prior_std.value(),opts.R_prior_fudge.value(), ndirections, nfib, nparams, opts.modelnum.value(), opts.fudge.value(), opts.f0.value(), opts.ardf0.value(), !opts.no_ard.value(), opts.rician.value(), gradnonlin, opts.updateproposalevery.value(), iters_step, (i*iters_step), 0, 0, 0, oldsignals_ptr, oldisosignals_ptr, angtmp_ptr, oldangtmp_ptr, fibres_ptr, multifibres_ptr, signals_ptr, isosignals_ptr,records_null,records_null,records_null,records_null,records_null,records_null,records_null,records_null, records_null);   
	   	sync_check("runmcmc_burnin_kernel");

 	   	gettimeofday(&t2,NULL);
    	   	timemcmc+=timeval_diff(&t2,&t1);
   	}

   	gettimeofday(&t1,NULL); 

   	if(nvox!=0){
   		curandStatus_t status = curandGenerateNormal(gen,randomsN_ptr,last_randoms,0,1);
		if (status != CURAND_STATUS_SUCCESS)
		{
			printf("Failure generating cuda random numbers: %d\n",status);
			exit(1);
		}
   		status = curandGenerateUniform(gen,randomsU_ptr,last_randoms); 	//generate randoms
		if (status != CURAND_STATUS_SUCCESS)
		{
			printf("Failure generating cuda random numbers: %d\n",status);
			exit(1);
		}
   	}
	
   	gettimeofday(&t2,NULL);
   	timecurand+=timeval_diff(&t2,&t1);

   	gettimeofday(&t1,NULL);

   	if(nvox!=0){
		runmcmc_kernel<<< Dim_Grid, Dim_Block, amount_shared >>>(datam_ptr, bvals_ptr, alpha_ptr, beta_ptr, randomsN_ptr, randomsU_ptr, opts.R_prior_mean.value(), opts.R_prior_std.value(),opts.R_prior_fudge.value(), ndirections, nfib, nparams, opts.modelnum.value(), opts.fudge.value(), opts.f0.value(), opts.ardf0.value(), !opts.no_ard.value(), opts.rician.value(), gradnonlin, opts.updateproposalevery.value(), last_step, (steps*iters_step), 0, 0, 0, oldsignals_ptr, oldisosignals_ptr, angtmp_ptr, oldangtmp_ptr, fibres_ptr, multifibres_ptr, signals_ptr, isosignals_ptr,records_null,records_null,records_null,records_null,records_null,records_null,records_null, records_null,records_null); 
   		sync_check("runmcmc_burnin_kernel");
   	}

   	gettimeofday(&t2,NULL);
   	timemcmc+=timeval_diff(&t2,&t1);

    	myfile << "TIME CURAND: " << timecurand << " seconds\n"; 
    	myfile << "TIME RUNMCMC: " << timemcmc << " seconds\n"; 
   
   	curandDestroyGenerator(gen);

	gettimeofday(&t_tot2,NULL);
    	time=timeval_diff(&t_tot2,&t_tot1);
   	myfile << "TIME TOTAL: " << time << " seconds\n"; 
	myfile << "-----------------------------------------------------" << "\n\n" ; 
	myfile.close();

   	sync_check("runmcmc_burnin");
}


void runmcmc_record(	//INPUT
			thrust::device_vector<float> 			datam_gpu,
			thrust::device_vector<float> 			bvals_gpu,
			thrust::device_vector<double> 			alpha_gpu,
			thrust::device_vector<double> 			beta_gpu,
			thrust::device_vector<FibreGPU> 		fibres_gpu,
			thrust::device_vector<MultifibreGPU> 		multifibres_gpu,
			thrust::device_vector<double>			signals_gpu,
			thrust::device_vector<double>			isosignals_gpu,
			const int 					ndirections,
			double 						seed,
			string 						output_file, 
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
	
	std::ofstream myfile;
	myfile.open (output_file.data(), ios::out | ios::app );
   	myfile << "--------- MCMC ALGORITHM PART RECORD IN GPU --------- " << "\n"; 	

   	struct timeval t1,t2,t_tot1,t_tot2;
   	double time,timecurand,timemcmc;
   	time=0;
   	timecurand=0;
   	timemcmc=0;

   	gettimeofday(&t_tot1,NULL);

   	size_t free,total;

	int totalrecords = (opts.njumps.value()/opts.sampleevery.value()); 
	
	int nvox = multifibres_gpu.size();
   	int nfib= opts.nfibres.value();
	int nparams;

	bool gradnonlin=opts.grad_file.set();

	if(opts.f0.value()) nparams=3+nfib*3;
	else nparams=2+nfib*3;	
	if(opts.modelnum.value()>=2) nparams++;
	if(opts.modelnum.value()==3) nparams++;	
	if(opts.rician.value()) nparams++;

	thrust::device_vector<double> angtmp_gpu;
	thrust::device_vector<double> oldangtmp_gpu;
	thrust::device_vector<double> oldsignals_gpu;
	thrust::device_vector<double> oldisosignals_gpu;
	
	angtmp_gpu.resize(nvox*ndirections*nfib);
	oldangtmp_gpu.resize(nvox*ndirections);
	oldsignals_gpu.resize(nvox*ndirections*nfib);
	oldisosignals_gpu.resize(nvox*ndirections);
   
   	unsigned int totalrandoms=(opts.njumps.value() * nvox * nparams);

   	cuMemGetInfo(&free,&total);
   	myfile << "Free memory Before Randoms: "<< free <<  " ---- Total memory: " << total << "\n";
   	//4 bytes each float, 2 random arrays, and 80% of total memory at this moment 
   	unsigned int maxrandoms=((free*0.8)/(4*2)); 

   	myfile << "Total randoms: " << totalrandoms << "\n"; 
   	myfile << "Max randoms: " << maxrandoms << "\n"; 
   
   	int steps; //num iter if not enough memory
   	int minrandoms; //min num of randoms ensamble
   	minrandoms= nvox * nparams;

   	int iters_step=0;
	int nrandoms=0;	

   	if(totalrandoms>maxrandoms){ 
		iters_step = maxrandoms / minrandoms; 		//iterations in each step
		nrandoms = iters_step*minrandoms;		//nrandoms for each step
		steps =  (opts.njumps.value()/iters_step);  	//repeat process steps times, no enough memory for all randoms 			
   	}else{ 
		nrandoms = totalrandoms;
		iters_step= opts.njumps.value();
		steps = 0;
  	}   
	if(nrandoms%2){						//CURAND must generates multiples of 2 randoms
		nrandoms++;
	}

   	myfile << "Process " << opts.njumps.value() << " iterations divided in "<< steps << " steps with "<< iters_step << " iterations in each one" << "\n";    

   	int last_step = opts.njumps.value() - (iters_step*steps);
   	int last_randoms = (last_step*minrandoms); 
	if(last_randoms%2){					//CURAND must generates multiples of 2 randoms
		last_randoms++;
	}

   	myfile << "Last step with " << last_step << " iterations" << "\n"; 
	
	thrust::device_vector<float> randomsN_gpu;
	thrust::device_vector<float> randomsU_gpu;	
	randomsN_gpu.resize(nrandoms);
	randomsU_gpu.resize(nrandoms);

   	cuMemGetInfo(&free,&total);
   	myfile << "Free memory after Malloc Randoms: "<< free <<  " ---- Total memory: " << total << "\n";
   
  	int blocks = nvox;        
  	dim3 Dim_Grid(blocks, 1);
  	dim3 Dim_Block(THREADS_BLOCK_MCMC,1);	//dimensions for MCMC   

   	myfile << "\n" << "NUM BLOCKS: " << blocks << "\n"; 
   	myfile << "THREADS PER BLOCK : " << THREADS_BLOCK_MCMC << "\n\n"; 	

   	curandGenerator_t gen;
   	curandCreateGenerator(&gen,CURAND_RNG_PSEUDO_DEFAULT);
   	curandSetPseudoRandomGeneratorSeed(gen,seed);

	//get pointers
	float *datam_ptr = thrust::raw_pointer_cast(datam_gpu.data());
	float *bvals_ptr = thrust::raw_pointer_cast(bvals_gpu.data());
	double *alpha_ptr = thrust::raw_pointer_cast(alpha_gpu.data());
	double *beta_ptr = thrust::raw_pointer_cast(beta_gpu.data());
	float *randomsN_ptr = thrust::raw_pointer_cast(randomsN_gpu.data());
	float *randomsU_ptr = thrust::raw_pointer_cast(randomsU_gpu.data());
	FibreGPU *fibres_ptr =  thrust::raw_pointer_cast(fibres_gpu.data());
	MultifibreGPU *multifibres_ptr = thrust::raw_pointer_cast(multifibres_gpu.data());
	double *signals_ptr = thrust::raw_pointer_cast(signals_gpu.data());
	double *isosignals_ptr = thrust::raw_pointer_cast(isosignals_gpu.data());

	double *angtmp_ptr = thrust::raw_pointer_cast(angtmp_gpu.data());
	double *oldangtmp_ptr = thrust::raw_pointer_cast(oldangtmp_gpu.data());
	double *oldsignals_ptr = thrust::raw_pointer_cast(oldsignals_gpu.data());
	double *oldisosignals_ptr = thrust::raw_pointer_cast(oldisosignals_gpu.data());
	
	float *rf0_ptr = thrust::raw_pointer_cast(rf0_gpu.data());
	float *rtau_ptr = thrust::raw_pointer_cast(rtau_gpu.data());
	float *rs0_ptr = thrust::raw_pointer_cast(rs0_gpu.data());
	float *rd_ptr = thrust::raw_pointer_cast(rd_gpu.data());
	float *rdstd_ptr = thrust::raw_pointer_cast(rdstd_gpu.data());
	float *rR_ptr = thrust::raw_pointer_cast(rR_gpu.data());	
	float *rth_ptr = thrust::raw_pointer_cast(rth_gpu.data());
	float *rph_ptr = thrust::raw_pointer_cast(rph_gpu.data());
	float *rf_ptr = thrust::raw_pointer_cast(rf_gpu.data());

	int amount_shared = (THREADS_BLOCK_MCMC)*sizeof(double) + (10*nfib + 2*nparams + 27)*sizeof(float) + (7*nfib + 21)*sizeof(int);

	myfile << "Shared Memory Used in runmcmc_record: " << amount_shared << "\n";

   	for(int i=0;i<steps;i++){

   		gettimeofday(&t1,NULL);

	   	curandStatus_t status = curandGenerateNormal(gen,randomsN_ptr,nrandoms,0,1);
		if (status != CURAND_STATUS_SUCCESS)
		{
			printf("Failure generating cuda random numbers: %d\n",status);
			exit(1);
		}
	   	status = curandGenerateUniform(gen,randomsU_ptr,nrandoms);	//generate randoms
		if (status != CURAND_STATUS_SUCCESS)
		{
			printf("Failure generating cuda random numbers: %d\n",status);
			exit(1);
		}

 	   	gettimeofday(&t2,NULL);
    	   	timecurand+=timeval_diff(&t2,&t1);

	   	gettimeofday(&t1,NULL);

	   	runmcmc_kernel<<< Dim_Grid, Dim_Block, amount_shared >>>(datam_ptr, bvals_ptr, alpha_ptr, beta_ptr, randomsN_ptr, randomsU_ptr, opts.R_prior_mean.value(), opts.R_prior_std.value(),opts.R_prior_fudge.value(), ndirections, nfib, nparams, opts.modelnum.value(), opts.fudge.value(), opts.f0.value(), opts.ardf0.value(), !opts.no_ard.value(), opts.rician.value(), gradnonlin, opts.updateproposalevery.value(), iters_step, (i*iters_step), opts.nburn.value(), opts.sampleevery.value(), totalrecords, oldsignals_ptr, oldisosignals_ptr, angtmp_ptr, oldangtmp_ptr, fibres_ptr, multifibres_ptr, signals_ptr, isosignals_ptr, rf0_ptr, rtau_ptr, rs0_ptr, rd_ptr, rdstd_ptr, rR_ptr, rth_ptr, rph_ptr, rf_ptr);
	   	sync_check("runmcmc_record_kernel");

 	   	gettimeofday(&t2,NULL);
    	   	timemcmc+=timeval_diff(&t2,&t1);
   	}

   	gettimeofday(&t1,NULL);

   	if(nvox!=0){
   		curandStatus_t status = curandGenerateNormal(gen,randomsN_ptr,last_randoms,0,1);
		if (status != CURAND_STATUS_SUCCESS)
		{
			printf("Failure generating cuda random numbers: %d\n",status);
			exit(1);
		}
   		status = curandGenerateUniform(gen,randomsU_ptr,last_randoms); 	//generate randoms
		if (status != CURAND_STATUS_SUCCESS)
		{
			printf("Failure generating cuda random numbers: %d\n",status);
			exit(1);
		}
   	}
	
   	gettimeofday(&t2,NULL);
   	timecurand+=timeval_diff(&t2,&t1);

   	gettimeofday(&t1,NULL);

   	if(nvox!=0){
		runmcmc_kernel<<< Dim_Grid, Dim_Block, amount_shared >>>(datam_ptr, bvals_ptr, alpha_ptr, beta_ptr,randomsN_ptr, randomsU_ptr, opts.R_prior_mean.value(), opts.R_prior_std.value(),opts.R_prior_fudge.value(), ndirections, nfib, nparams, opts.modelnum.value(), opts.fudge.value(), opts.f0.value(), opts.ardf0.value(), !opts.no_ard.value(), opts.rician.value(), gradnonlin, opts.updateproposalevery.value(), last_step, (steps*iters_step), opts.nburn.value(), opts.sampleevery.value(), totalrecords, oldsignals_ptr, oldisosignals_ptr, angtmp_ptr, oldangtmp_ptr, fibres_ptr, multifibres_ptr, signals_ptr, isosignals_ptr, rf0_ptr, rtau_ptr, rs0_ptr, rd_ptr, rdstd_ptr, rR_ptr, rth_ptr, rph_ptr, rf_ptr);   
   		sync_check("runmcmc_record_kernel");
   	}

   	gettimeofday(&t2,NULL);
   	timemcmc+=timeval_diff(&t2,&t1);


    	myfile << "TIME CURAND: " << timecurand << " seconds\n"; 
    	myfile << "TIME RUNMCMC: " << timemcmc << " seconds\n"; 
   
   	curandDestroyGenerator(gen);

	gettimeofday(&t_tot2,NULL);
    	time=timeval_diff(&t_tot2,&t_tot1);
   	myfile << "TIME TOTAL: " << time << " seconds\n"; 
	myfile << "-----------------------------------------------------" << "\n" ;
	myfile.close(); 
	
   	sync_check("runmcmc_record");
}
