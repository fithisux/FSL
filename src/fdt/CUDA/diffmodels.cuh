/*  diffmodels.cuh

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

#include <device_vector.h>

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
			thrust::device_vector<float>&	params_gpu);

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
			thrust::device_vector<float>&	params_gpu);

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
			thrust::device_vector<float>&	params_gpu);

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
			thrust::host_vector<float>&	tau);


__device__ void cf_PVM_single(		//INPUT
					const float*			params,
					const float*			data,
					const float*			bvecs, 
					const float*			bvals,	
					const int			ndirections,
					const int			nfib,
					const int 			nparams, 
					const bool 			m_include_f0,
					const int			idSubVOX,
					float*				reduction,					
					float* 				fs,
					float*				x,
					float*				_d,
					float*				sumf,
					//OUTPUT
					double* 			cfv);

__device__ void grad_PVM_single(	//INPUT
					const float*			params,
					const float*			data,
					const float*			bvecs, 
					const float*			bvals,
					const int			ndirections,
					const int			nfib,
					const int 			nparams,
					const bool 			m_include_f0,
					const int			idSubVOX,
					float*				J,
					float*				reduction,					
					float* 				fs,
					float*				x,
					float* 				_d,
					float* 				sumf,
					//OUTPUT
					float*				grad);

__device__ void hess_PVM_single(	//INPUT
					const float*			params,
					const float*			bvecs, 
					const float*			bvals,
					const int			ndirections,
					const int			nfib,
					const int 			nparams,
					const bool 			m_include_f0,
					const int			idSubVOX,
					float*				J,
					float*				reduction,
					float* 				fs,
					float*				x,
					float* 				_d,
					float* 				sumf,
					//OUTPUT
					float*				hess);

__device__ void cf_PVM_single_c(	//INPUT
					const float*			params,
					const float*			data,
					const float*			bvecs, 
					const float*			bvals,
					const int			ndirections,
					const int			nfib,
					const int 			nparams, 
					const bool 			m_include_f0,
					const int			idSubVOX,
					float*				reduction,
					float* 				fs,
					float*				x,
					float* 				_d,
					float* 				sumf,
					//OUTPUT
					double* 			cfv);


__device__ void grad_PVM_single_c(	//INPUT
					const float*			params,
					const float*			data,
					const float*			bvecs, 
					const float*			bvals,
					const int			ndirections,
					const int			nfib,
					const int 			nparams,
					const bool 			m_include_f0,
					const int			idSubVOX,
					float*				J,
					float*				reduction,					
					float* 				fs,
					float* 				f_deriv,
					float*				x,
					float* 				_d,
					float* 				sumf,
					//OUTPUT
					float*				grad);

__device__ void hess_PVM_single_c(	//INPUT
					const float*			params,
					const float*			bvecs, 
					const float*			bvals,
					const int			ndirections,
					const int			nfib,
					const int 			nparams,
					const bool 			m_include_f0,
					const int			idSubVOX,
					float*				J,
					float*				reduction,					
					float* 				fs,
					float* 				f_deriv,
					float*				x,
					float* 				_d,
					float* 				sumf,
					//OUTPUT
					float*				hess);

__device__ void cf_PVM_multi(		//INPUT
					const float*			params,
					const float*			data,
					const float*			bvecs, 
					const float*			bvals,
					const float			R,
					const float			invR,	
					const int			ndirections,
					const int			nfib,
					const int 			nparams, 
					const bool 			m_include_f0,
					const int			idSubVOX,
					const int			Gamma_for_ball_only,
					float*				reduction,					
					float* 				fs,
					float*				x,
					float* 				_a,
					float* 				_b,
					float* 				sumf,
					//OUTPUT
					double*				cfv);

__device__ void grad_PVM_multi(		//INPUT
					const float*			params,
					const float*			data,
					const float*			bvecs, 
					const float*			bvals,
					const float			R,
					const float			invR,
					const int			ndirections,
					const int			nfib,
					const int 			nparams,
					const bool 			m_include_f0,
					const int			idSubVOX,
					const int			Gamma_for_ball_only,
					float*				J,
					float*				reduction,					
					float* 				fs,
					float*				x,
					float* 				_a,
					float* 				_b,
					float* 				sumf,
					//OUTPUT
					float*				grad);

__device__ void hess_PVM_multi(		//INPUT
					const float*			params,
					const float*			bvecs, 
					const float*			bvals,
					const float			R,
					const float			invR,
					const int			ndirections,
					const int			nfib,
					const int 			nparams,
					const bool 			m_include_f0,
					const int			idSubVOX,
					const int			Gamma_for_ball_only,
					float*				J,
					float*				reduction,					
					float* 				fs,
					float*				x,
					float* 				_a,
					float*				_b,
					float* 				sumf,
					//OUTPUT
					float*				hess);
