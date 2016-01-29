/*  mcmc.h

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

#if !defined(mcmc_h)
#define mcmc_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "gsoptions.h"
#include "newimage/newimageall.h"
#include "newmat.h"
#include "design.h"
#include "gam.h"

using namespace NEWIMAGE;
using namespace NEWMAT;

namespace Gs {
    
  class Mcmc
    {
    public:

      // constructor
      Mcmc(const ColumnVector pcopedata, 
	   const ColumnVector pvarcopedata, 
	   const Design& pdesign, 
	   ColumnVector pgamma_mean, 
	   ColumnVector pgamma_S, 
	   ColumnVector pbeta_b, 
	   ColumnVector pbeta_c, 
	   ColumnVector pm_mean, 
	   ColumnVector pm_var, 
	   Matrix& pgamma_samples, 
	   int pnsamples, bool pswitched = false) :
	copedata(pcopedata),
	varcopedata(pvarcopedata),
	design(pdesign),
	gamma_mean(pgamma_mean),
	gamma_S(),
	gamma_samples(pgamma_samples),
	beta_b(pbeta_b),
	beta_c(pbeta_c),
	m_mean(pm_mean),
	m_var(pm_var),
	ngs(design.getngs()),
	nevs(design.getnevs()),
	ntpts(design.getntpts()),
	nsamples(pnsamples),
	opts(GsOptions::getInstance()),
	gam(Gam::getInstance()),	
	switched(pswitched),
	prior_dominating(false)
	{ 
	  reshape(gamma_S,pgamma_S,nevs,nevs);
	}

      // load data from file in from file and set up starting values
      void setup();

      // runs the chain
      void run();

      // jumps
      void jump(bool relax);

      // sample chain
      void sample(int samp);

      // getters
      const int getnsamples() const { return nsamples; }      
      const bool is_prior_dominating() const { return prior_dominating; }

      // Destructor
      virtual ~Mcmc() {}

      ColumnVector c_samples; 
    private:
    
      void beta_jump(bool relax);
      void m_jump(bool relax);
      void gamma_jump(bool relax);

      void beta_jump_switched(bool relax);
      void m_jump_switched(bool relax);
      void gamma_jump_switched(bool relax);
      
      Mcmc();
      const Mcmc& operator=(Mcmc& mcmc);     
      Mcmc(Mcmc& mcmc);

      const ColumnVector copedata;
      const ColumnVector varcopedata;
      const Design& design;

      ColumnVector gamma_mean;
      Matrix gamma_S;
      ColumnVector gamma_latest;
      Matrix& gamma_samples;

      ColumnVector beta_b;
      ColumnVector beta_c;
      ColumnVector beta_latest;

      ColumnVector m_mean;
      ColumnVector m_var;
      ColumnVector m_latest;

      int ngs;
      int nevs;
      int ntpts;

      int nsamples;

      GsOptions& opts;

      Gam& gam;

      bool switched;

      float c_latest;

      bool prior_dominating;
      
    };
}   
#endif

