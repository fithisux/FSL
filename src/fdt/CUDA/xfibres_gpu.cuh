/*  xfibres_gpu.cuh

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

#include "newimage/newimageall.h"
#include <host_vector.h>
#include <device_vector.h> 

#include "fibre_gpu.h"

//implemented and used in xfibres_gpu.cu
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
		int				ndirections,
		string 				output_file,
		//OUTPUT
		thrust::device_vector<float>&	params_gpu,
		thrust::host_vector<int>&	vox_repeat,
		int&				nrepeat);

//implemented and used in xfibres_gpu.cu
void prepare_data_gpu_FIT(	//INPUT
				const Matrix				datam,
				const Matrix				bvecs,
				const Matrix				bvals,
				const Matrix	 			gradm, 
				//OUTPUT
				vector<ColumnVector>&			datam_vec,
				vector<Matrix>&				bvecs_vec,
				vector<Matrix>&				bvals_vec,
				thrust::host_vector<float>&   		datam_host,	
				thrust::host_vector<float>&		bvecs_host,				
				thrust::host_vector<float>&		bvals_host,
				thrust::host_vector<double>&		alpha_host,
				thrust::host_vector<double>&		beta_host,
				thrust::host_vector<float>&		params_host,
				thrust::host_vector<float>&		tau_host);


//implemented and used in xfibres_gpu.cu
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
					thrust::host_vector<float>&   		datam_repeat_host,
					thrust::host_vector<float>&		bvecs_repeat_host,				
					thrust::host_vector<float>&		bvals_repeat_host,
					thrust::host_vector<float>&		params_repeat_host);

//implemented and used in xfibres_gpu.cu
void mix_params(	//INPUT
			thrust::host_vector<float>   		params_repeat_gpu,
			thrust::host_vector<int>			vox_repeat,
			int						nrepeat,
			int						nvox,
			//INPUT-OUTPUT
			thrust::device_vector<float>&   		params_gpu);


//implemented and used in xfibres_gpu.cu
void prepare_data_gpu_MCMC(	//INPUT
				int 					nvox,
				int					ndirections,
				int 					nfib,
				//OUTPUT
				thrust::host_vector<double>&		signals_host,
				thrust::host_vector<double>&		isosignals_host,
				thrust::host_vector<FibreGPU>& 		fibres_host,
				thrust::host_vector<MultifibreGPU>& 	multifibres_host);

//implemented and used in xfibres_gpu.cu
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
					thrust::device_vector<float>&			rf_gpu);

//implemented and used in xfibres_gpu.cu
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
				int						idpart);

