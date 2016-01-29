/*  samples.h

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


using namespace Xfibres;

////////////////////////////////////////////
//       MCMC SAMPLE STORAGE
////////////////////////////////////////////

class Samples{
  	xfibresOptions& opts;
  	Matrix m_dsamples;
  	Matrix m_d_stdsamples;
	Matrix m_Rsamples;
  	Matrix m_S0samples;
  	Matrix m_f0samples;

	//   // storing signal
	//   Matrix m_mean_sig;
	//   Matrix m_std_sig;
	//   Matrix m_sig2;

  	vector<Matrix> m_thsamples;
  	vector<Matrix> m_phsamples;
  	vector<Matrix> m_fsamples;
  	vector<Matrix> m_lamsamples;

  	//for storing means
  	RowVector m_mean_dsamples;
  	RowVector m_mean_d_stdsamples;
	RowVector m_mean_Rsamples;
  	RowVector m_mean_S0samples;
  	RowVector m_mean_f0samples;
  	RowVector m_mean_tausamples;
  	vector<Matrix> m_dyadic_vectors;
  	vector<RowVector> m_mean_fsamples;
  	vector<RowVector> m_mean_lamsamples;

  	//float m_sum_d;  changed GPU version
  	//float m_sum_d_std;  changed GPU version
  	//float m_sum_S0;  changed GPU version
  	//float m_sum_f0;  changed GPU version
  	//float m_sum_tau;  changed GPU version
  	//vector<SymmetricMatrix> m_dyad;  changed GPU version
  	//vector<float> m_sum_f;  changed GPU version
  	//vector<float> m_sum_lam;  changed GPU version
  	//ColumnVector m_vec;  changed GPU version

  	/////////////// GPU version /////////////////////
  	float *m_sum_d;
  	float *m_sum_S0;
  	float *m_sum_d_std;
	float *m_sum_R;
  	float *m_sum_f0;
  	float *m_sum_tau;

  	vector<SymmetricMatrix> *m_dyad;
  	vector<float>  *m_sum_f;
  	vector<float> *m_sum_lam;
  	ColumnVector *m_vec;
  	////////////////////////////////////////////////
  
  	int m_nsamps;

	public:

  	Samples(int nvoxels,int nmeasures);
    
	void record(float rd,float rf0,float rtau,float rdstd,float rR,float rs0,float *rth,float *rph, float *rf, int vox, int samp);
  
  	void finish_voxel(int vox);
  
  	void save(int idpart);
};
