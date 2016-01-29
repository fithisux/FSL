/*  MELODIC - Multivariate exploratory linear optimized decomposition into 
              independent components
    
    melhlprfns.cc - misc functions

    Christian F. Beckmann, FMRIB Analysis Group
    
    Copyright (C) 1999-2013 University of Oxford */

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

#ifndef __MELODICHLPR_h
#define __MELODICHLPR_h

#include "newimage/newimageall.h"

	#ifdef __APPLE__
	#include <mach/mach.h>
	#define mmsg(msg) { \
	  struct task_basic_info t_info; \
	  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT; \
	  if (KERN_SUCCESS == task_info(mach_task_self(), TASK_BASIC_INFO, (task_info_t) &t_info, &t_info_count)) \
		{ \
			cout << " MEM: " << msg << " res: " << t_info.resident_size/1000000 << " virt: " << t_info.virtual_size/1000000 << "\n"; \
			} \
	}
	#else
	#define mmsg(msg) { \
	   cout << msg; \
	}
	#endif

using namespace NEWIMAGE;

namespace Melodic{

  void update_mask(volume<float>& mask, Matrix& Data);
  void del_vols(volume4D<float>& in, int howmany);

  Matrix smoothColumns(const Matrix& inp);
  Matrix calc_FFT(const Matrix& Mat, const bool logpwr = 0);

  Matrix convert_to_pbsc(Matrix& Mat);

  RowVector varnorm(Matrix& in, int dim = 30, float level = 1.6, int econ = 20000);
       void varnorm(Matrix& in, const RowVector& vars);
  RowVector varnorm(Matrix& in, SymmetricMatrix& Corr, int dim = 30, float level = 1.6, int econ = 20000);

  Matrix SP2(const Matrix& in, const Matrix& weights, int econ = 20000);
  void SP3(Matrix& in, const Matrix& weights);

  RowVector Feta(int n1,int n2);
  RowVector cumsum(const RowVector& Inp);

  Matrix corrcoef(const Matrix& in1, const Matrix& in2);
  Matrix corrcoef(const Matrix& in1, const Matrix& in2, const Matrix& part);
  float calc_white(const Matrix& tmpE, const RowVector& tmpD, const RowVector& PercEV, int dim, Matrix& param, Matrix& paramS, Matrix& white, Matrix& dewhite);
  float calc_white(const Matrix& tmpE, const RowVector& tmpD, const RowVector& PercEV, int dim, Matrix& white, Matrix& dewhite);
  void calc_white(const Matrix& tmpE, const RowVector& tmpD, int dim, Matrix& param, Matrix& paramS, Matrix& white, Matrix& dewhite);
  void calc_white(const Matrix& tmpE, const RowVector& tmpD, int dim, Matrix& white, Matrix& dewhite);
  void calc_white(const SymmetricMatrix& Corr, int dim, Matrix& white, Matrix& dewhite);
  
  void std_pca(const Matrix& Mat, SymmetricMatrix& Corr, Matrix& evecs, RowVector& evals, int econ = 20000);
  void std_pca(const Matrix& Mat, const Matrix& weights, SymmetricMatrix& Corr, Matrix& evecs, RowVector& evals, int econ = 20000);
  void em_pca(const Matrix& Mat, Matrix& evecs, RowVector& evals, int num_pc = 1, int iter = 20);
  void em_pca(const Matrix& Mat, Matrix& guess, Matrix& evecs, RowVector& evals, int num_pc = 1, int iter = 20);

  float rankapprox(const Matrix& Mat, Matrix& cols, Matrix& rows, int dim = 1);
  RowVector krfact(const Matrix& Mat, Matrix& cols, Matrix& rows);
  RowVector krfact(const Matrix& Mat, int colnum, Matrix& cols, Matrix& rows);
  Matrix krprod(const Matrix& cols, const Matrix& rows);
  Matrix krapprox(const Matrix& Mat, int size_col, int dim = 1);

  void adj_eigspec(const RowVector& in, RowVector& out1, RowVector& out2, RowVector& out3, int& out4, int num_vox, float resels);
  void adj_eigspec(const RowVector& in, RowVector& out1, RowVector& out2);

  int ppca_dim(const Matrix& in, const Matrix& weights, Matrix& PPCA, RowVector& AdjEV, RowVector& PercEV, SymmetricMatrix& Corr, Matrix& tmpE, RowVector &tmpD, float resels, string which);
  int ppca_dim(const Matrix& in, const Matrix& weights, Matrix& PPCA, RowVector& AdjEV, RowVector& PercEV, float resels, string which);
  int ppca_dim(const Matrix& in, const Matrix& weights, float resels, string which);
  ColumnVector ppca_select(Matrix& PPCAest, int& dim, int maxEV, string which);
  Matrix ppca_est(const RowVector& eigenvalues, const int N1, const float N2);
  Matrix ppca_est(const RowVector& eigenvalues, const int N);

  ColumnVector acf(const ColumnVector& in, int order);
  ColumnVector pacf(const ColumnVector& in, int maxorder = 1);
  Matrix est_ar(const Matrix& Mat, int maxorder);
  ColumnVector gen_ar(const ColumnVector& in, int maxorder = 1);
  Matrix gen_ar(const Matrix& in, int maxorder);
  Matrix gen_arCorr(const Matrix& in, int maxorder);
 
	class basicGLM{
		public:
		
			//constructor
			basicGLM(){}
		
			//destructor
			~basicGLM(){}
		
			void olsfit(const Matrix& data, const Matrix& design, 
				const Matrix& contrasts, int DOFadjust = -1);

			inline Matrix& get_t(){return t;}
			inline Matrix& get_z(){return z;}
			inline Matrix& get_p(){return p;}
			inline Matrix& get_f_fmf(){return f_fmf;}
			inline Matrix& get_pf_fmf(){return pf_fmf;}
			inline Matrix& get_cbeta(){return cbeta;}
			inline Matrix& get_beta(){return beta;}
			inline Matrix& get_varcb(){return varcb;}
			inline Matrix& get_sigsq(){return sigsq;}
			inline Matrix& get_residu(){return residu;}
			inline int get_dof(){return dof;}
			
		private:
			Matrix beta;
			Matrix residu;
			Matrix sigsq;
			Matrix varcb;
			Matrix cbeta;
			Matrix f_fmf, pf_fmf;
			int dof;
			Matrix t;
			Matrix z;
			Matrix p;
  };
//	Matrix glm_ols(const Matrix& dat, const Matrix& design);
}

#endif
