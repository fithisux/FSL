/*  gsmanager.h

    Mark Woolrich, Tim Behrens - FMRIB Image Analysis Group

    Copyright (C) 2002 University of Oxford  */

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

#if !defined(gsmanager_h)
#define gsmanager_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "gsoptions.h"
#include "newimage/newimageall.h"
#include "design.h"

using namespace NEWIMAGE;
using namespace MISCMATHS;

namespace Gs {
    
  // Give this class a file containing
  class Gsmanager
    {
    public:

      // constructor
      Gsmanager() : 
	opts(GsOptions::getInstance()),
	nmaskvoxels(0)
	{
	}

      // load data from file in from file and set up starting values
      void setup();

      // initialise
      void initialise();

      void run();

      // saves results in logging directory
      void save();

      // Destructor
      virtual ~Gsmanager() {}
 
    private:

      const Gsmanager& operator=(Gsmanager& par);     
      Gsmanager(Gsmanager& des);      

      void multitfit(const Matrix& x, ColumnVector& m, SymmetricMatrix& covar, float& v, bool fixmean=false) const;

      float log_likelihood(float beta, const ColumnVector& gam, const ColumnVector& y, const Matrix& z, const ColumnVector& S);
      float log_likelihood_outlier(float beta, float beta_outlier, const ColumnVector& gam, const ColumnVector& y, const Matrix& z, const ColumnVector& S, float global_prob_outlier, const ColumnVector& prob_outlier);

      float marg_posterior_energy(float x, const ColumnVector& y, const Matrix& z, const ColumnVector& S);
      float marg_posterior_energy_outlier(float logbeta, float logbeta_outlier, const ColumnVector& y, const Matrix& z, const ColumnVector& S, const ColumnVector& prob_outlier);

      float solveforbeta(const ColumnVector& y, const Matrix& z, const ColumnVector& S);

      bool pass_through_to_mcmc(float zlowerthresh, float zupperthresh, int px, int py, int pz);

      // functions to calc contrasts
      void ols_contrasts(const ColumnVector& gammean, const SymmetricMatrix& gamS, int px, int py, int pz);
      void fe_contrasts(const ColumnVector& gammean, const SymmetricMatrix& gamS, int px, int py, int pz);
      void flame1_contrasts(const ColumnVector& gammean, const SymmetricMatrix& gamS, int px, int py, int pz);
      void flame1_contrasts_with_outliers(const ColumnVector& mn, const SymmetricMatrix& covariance, int px, int py, int pz);
      void flame2_contrasts(const Matrix& gamsamples, int px, int py, int pz);

      void t_ols_contrast(const ColumnVector& gammean, const SymmetricMatrix& gamS, const RowVector& tcontrast, float& cope, float& varcope, float& t, float& dof, float& z, int px, int py, int pz);	
      void f_ols_contrast(const ColumnVector& gammean, const SymmetricMatrix& gamS, const Matrix& fcontrast, float& f, float& dof1, float& dof2, float& z, int px, int py, int pz);
	
	
      void t_mcmc_contrast(const Matrix& gamsamples, const RowVector& tcontrast, float& cope, float& varcope, float& t, float& dof, float& z, int px, int py, int pz);
      void f_mcmc_contrast(const Matrix& gamsamples, const Matrix& fcontrast, float& f, float& dof1, float& dof2, float& z, int px, int py, int pz);


      // voxelwise functions to perform the different inference approaches
      void fixed_effects_onvoxel(const ColumnVector& Y, const Matrix& z, const ColumnVector& S, ColumnVector& gam, SymmetricMatrix& gamcovariance);
      void flame_stage1_onvoxel(const vector<ColumnVector>& Yg, const ColumnVector& Y, const vector<Matrix>& zg, const Matrix& z, const vector<ColumnVector>& Sg, const ColumnVector& S, ColumnVector& beta, ColumnVector& gam, SymmetricMatrix& gamcovariance, vector<float>& marg, vector<ColumnVector>& weights_g, int& nparams, int px, int py, int pz);
      void flame_stage1_onvoxel_inferoutliers(const vector<ColumnVector>& Yg, const ColumnVector& Y, const vector<Matrix>& zg, const Matrix& z, const vector<ColumnVector>& Sg, const ColumnVector& S, ColumnVector& beta, ColumnVector& gam, SymmetricMatrix& gamcovariance, ColumnVector& global_prob_outlier, vector<ColumnVector>& prob_outlier_g,  ColumnVector& prob_outlier, ColumnVector& beta_outlier, vector<float>& marg, vector<ColumnVector>& weights_g, int& nparams, vector<bool>& no_outliers, int px, int py, int pz);
      void flame_stage1_inferoutliers();
      void init_flame_stage1_inferoutliers();

      // functions to perform the different inference approaches
      void fixed_effects(); 
      void ols(); 
      void flame_stage1();
      void flame_stage2();
      void do_kmeans(const Matrix& data,vector<int>& z,const int k,Matrix& means);
      void randomise(vector< pair<float,int> >& r);
      vector< pair<float,int> > randomise(const int n);

      void regularise_flame2_contrasts();
     
      volume<float> mcmc_mask;
     
      // intermediates
      Design design;

      vector<volume<float> > beta_b;
      vector<volume<float> > beta_c;
      vector<volume<float> > beta_mean;
      vector<volume<float> > beta_outlier_mean;
      vector<volume<float> > global_prob_outlier_mean;
      vector<volume4D<float> > prob_outlier_mean;
      vector<volume4D<float> > weights;

      volume4D<float> cov_pes;

      // outputs

      vector<volume<float> > pes;
      vector<volume<float> > ts;
      vector<volume<float> > tdofs;
      vector<volume<float> > zts;
      vector<volume<float> > zflame1upperts;
      vector<volume<float> > zflame1lowerts;
      vector<volume<float> > tcopes;
      vector<volume<float> > tvarcopes;      
      vector<volume<float> > fs;
      vector<volume<float> > fdof1s;
      vector<volume<float> > fdof2s;
      vector<volume<float> > zfs;
      vector<volume<float> > zflame1upperfs;
      vector<volume<float> > zflame1lowerfs;

      // intermediates
      int ngs;
      int nevs;
      int ntpts;
      int xsize;
      int ysize;
      int zsize;

      GsOptions& opts;

      int nmaskvoxels;

      bool dofpassedin;
    };

  bool compare(const pair<float,int> &r1,const pair<float,int> &r2);

}   
#endif







