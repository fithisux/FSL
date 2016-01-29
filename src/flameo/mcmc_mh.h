/*  mcmc_mh.h

    Mark Woolrich, Tim Behrens, FMRIB Image Analysis Group

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

#if !defined(mcmc_mh_h)
#define mcmc_mh_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "gsoptions.h"
#include "newimage/newimageall.h"
#include "newmat.h"
#include "design.h"

using namespace NEWIMAGE;
using namespace NEWMAT;

namespace Gs {
    
  class Mcmc_Mh
    {
    public:

      // constructor
      Mcmc_Mh(const ColumnVector pcopedata, 
	      const ColumnVector pvarcopedata, 
	      const ColumnVector pdofvarcopedata, 
	      const Design& pdesign, 
	      ColumnVector pgamma_mean, 
	      SymmetricMatrix pgamma_S, 
	      ColumnVector pbeta_b, 
	      ColumnVector pbeta_c, 
	      Matrix& pgamma_samples,
	      Matrix& pbeta_samples,
	      Matrix& pphi_samples,
	      ColumnVector& plikelihood_samples,
	      vector<ColumnVector>& pss_samples,
	      int pnsamples,
	      int px, int py, int pz,
	      const ColumnVector& pprob_outlier,
	      const vector<float>& pglobal_prob_outlier,
	      const vector<float>& pbeta_outlier,
	      bool pinfer_outliers) :
      copedata(pcopedata),
	varcopedata(pvarcopedata),
	dofvarcopedata(pdofvarcopedata),
	design(pdesign),
	design_matrix(pdesign.getdm(px,py,pz)),
	nevs(design_matrix.Ncols()),
	gamma_mean(pgamma_mean),
	gamma_S(pgamma_S),
	gamma_latest(nevs),
	gamma_samples(pgamma_samples),
	gamma_nrejected(nevs),
	gamma_naccepted(nevs),
	gamma_proposal_std(nevs),
	beta_b(pbeta_b),
	beta_c(pbeta_c),
	beta_latest(design.getngs()),
	beta_samples(pbeta_samples),
	beta_nrejected(design.getngs()),
	beta_naccepted(design.getngs()),
	beta_proposal_std(design.getngs()),
	beta_prior_energy_old(design.getngs()),
	phi_latest(design.getntpts()),
	phi_samples(pphi_samples),
	phi_nrejected(design.getntpts()),
	phi_naccepted(design.getntpts()),
	phi_proposal_std(design.getntpts()),
	phi_prior_energy_old(design.getntpts()),
	likelihood_energy_old(0.0),
	likelihood_samples(plikelihood_samples),
	//	ss_samples(pss_samples),
	ngs(design.getngs()),
	ntpts(design.getntpts()),
	nsamples(pnsamples),
	opts(GsOptions::getInstance()),
	delta(design.getntpts()),
	sampcount(0),
	subsampcount(0),
	sumovere(design.getntpts()),
	prec_ontwo(design.getntpts()),
	logprec_ontwo(design.getntpts()),
	uncertainty_in_varcopes(opts.dofvarcopefile.value() != string("")),
	voxx(px),
	voxy(py),
	voxz(pz),
	prob_outlier(pprob_outlier),
	global_prob_outlier(pglobal_prob_outlier),
	beta_outlier(pbeta_outlier),
	infer_outliers(pinfer_outliers)
	{ 
	}

      // load data from file in from file and set up starting values
      void setup();

      // runs the chain
      void run();

      // jumps
      void jump();

      // sample chain
      void sample(int samp);
      
      // DIC
/*       void dic(float& DIC, float& pd); */

      // getters
      const int getnsamples() const { return nsamples; }   

      const ColumnVector& getgamma_naccepted() const { return gamma_naccepted; }
      const ColumnVector& getgamma_nrejected() const { return gamma_nrejected; } 
      const ColumnVector& getbeta_naccepted() const { return beta_naccepted; }
      const ColumnVector& getbeta_nrejected() const { return beta_nrejected; } 
      const ColumnVector& getphi_naccepted() const { return phi_naccepted; }
      const ColumnVector& getphi_nrejected() const { return phi_nrejected; }  
      // Destructor
      virtual ~Mcmc_Mh() {}

      ColumnVector c_samples; 
    private:
    
      void beta_jump();
      void phi_jump();
      void gamma_jump();
      void all_jump();

      float likelihood_energy(const int echanged, const float gamma_old, const bool betachanged);

      float likelihood_energy_phichanged(const int t);

      float beta_prior_energy(int g);
      float phi_prior_energy(int g);

/*       void sample_sumsquares(int samp); */
/*       float sumsquare_residuals(const Matrix& pdm, const ColumnVector& pdata, const ColumnVector& ppes); */

      Mcmc_Mh();
      const Mcmc_Mh& operator=(Mcmc_Mh& mcmc_mh);     
      Mcmc_Mh(Mcmc_Mh& mcmc_mh);

      const ColumnVector copedata;
      const ColumnVector varcopedata;
      const ColumnVector dofvarcopedata;
      const Design& design;
      Matrix design_matrix;
      int nevs;

      // mean
      ColumnVector gamma_mean;

      // Covariance:
      SymmetricMatrix gamma_S;

      ColumnVector gamma_latest;
      Matrix& gamma_samples;
      ColumnVector gamma_nrejected;
      ColumnVector gamma_naccepted;
      ColumnVector gamma_proposal_std;

      ColumnVector beta_b;
      ColumnVector beta_c;
      ColumnVector beta_latest;
      Matrix& beta_samples;
      ColumnVector beta_nrejected;
      ColumnVector beta_naccepted;
      ColumnVector beta_proposal_std;
      ColumnVector beta_prior_energy_old;

      ColumnVector phi_latest;
      Matrix& phi_samples;
      ColumnVector phi_nrejected;
      ColumnVector phi_naccepted;
      ColumnVector phi_proposal_std;
      ColumnVector phi_prior_energy_old;

      float likelihood_energy_old;

      ColumnVector& likelihood_samples;

/*       vector<ColumnVector>& ss_samples; */

      int ngs;
      int ntpts;

      int nsamples;

      GsOptions& opts;

      float c_latest;

      ColumnVector delta;
      
      int sampcount;
      int subsampcount;

      ColumnVector sumovere;
      ColumnVector prec_ontwo;
      ColumnVector logprec_ontwo;

      bool uncertainty_in_varcopes;

      int voxx; int voxy; int voxz;

      const ColumnVector& prob_outlier;
      const vector<float>& global_prob_outlier;
      const vector<float>& beta_outlier;

      bool infer_outliers;
 
    };
}   
#endif

