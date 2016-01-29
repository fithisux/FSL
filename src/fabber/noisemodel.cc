/*  noisemodel.cc - Class implementation for generic noise models

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

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

#include "noisemodel.h"

// Handy mathematical function, used in some free energy calculations
#include <math.h>

// Calculate log-gamma from a Taylor expansion; good to one part in 2e-10.
double gammaln(double x)
{
  ColumnVector series(7);
  series << 2.5066282746310005 << 76.18009172947146
	 << -86.50532032941677 << 24.01409824083091
	 << -1.231739572450155 << 0.1208650973866179e-2
	 << -0.5395239384953e-5;

  double total = 1.000000000190015;
  for (int i = 2; i <= series.Nrows(); i++)
    total += series(i) / (x+i-1);

  return log( series(1) * total / x ) + (x + 0.5)*log(x + 5.5) - x - 5.5;
}

#include "noisemodel_ar.h"
#include "noisemodel_white.h"

NoiseModel* NoiseModel::NewFromName(const string& name, ArgsType& args)
{
    // Update this to add your own models to the code
    
  if (name == "ar")
    {
      string nPhis = args.ReadWithDefault("num-echoes","(default)");
      if (nPhis == "(default)")
	{
	  nPhis = "1";
	  Warning::IssueOnce("Defaulting to --num-echoes=1");
	}
      string ar1CrossTerms = args.ReadWithDefault("ar1-cross-terms","none");
      
      return new Ar1cNoiseModel(ar1CrossTerms, convertTo<int>(nPhis));
    }
  else if (name == "ar1" || name == "ar1c")
    {
      Warning::IssueOnce("--noise="+name+" is depreciated; use --noise=ar --num-echoes=2 for dual-echo ASL.");
      string ar1CrossTerms = args.ReadWithDefault("ar1-cross-terms","none");
      string nPhis = args.ReadWithDefault("num-echoes","2");
      
      return new Ar1cNoiseModel(ar1CrossTerms, convertTo<int>(nPhis));
    }
    else if (name == "white")
    {
      //      string pattern = args.ReadWithDefault("noise-pattern","1");
      //      return new WhiteNoiseModel(pattern);
      return new WhiteNoiseModel(args);
    }
    // Your models go here!
    else
    {
      throw Invalid_option("Unrecognized noise model '" + name + "'");
    }
}

// ARD stuff
double NoiseModel::SetupARD(vector<int> ardindices,
			  const MVNDist& theta,
			  MVNDist& thetaPrior) const {
  Tracer_Plus tr("Noisemodel::SetupARD");
  double Fard=0;

  if (~ardindices.empty()) {
    SymmetricMatrix PriorPrec;
    PriorPrec = thetaPrior.GetPrecisions();
    SymmetricMatrix PostCov = theta.GetCovariance();

    for (int i=1; i<= ardindices.size(); i++)
      {
	PriorPrec(ardindices[i],ardindices[i]) = 1e-12; //set prior to be initally non-informative
	thetaPrior.means(ardindices[i]) = 0;
	
	//set the Free energy contribution from ARD term
	double b = 2/(theta.means(ardindices[i])*theta.means(ardindices[i]) + PostCov(ardindices[i],ardindices[i]));
	Fard += -1.5*(log(b) + digamma(0.5)) - 0.5 - gammaln(0.5) - 0.5*log(b); //taking c as 0.5 - which it will be!
      }

	thetaPrior.SetPrecisions(PriorPrec);
  }

  return Fard;
}

double NoiseModel::UpdateARD(vector<int> ardindices,
			  const MVNDist& theta,
			  MVNDist& thetaPrior) const {
  Tracer_Plus tr("Noisemodel::UpdateARD");
  double Fard=0;

  if (~ardindices.empty()) {
    SymmetricMatrix PriorCov;
      SymmetricMatrix PostCov;
      PriorCov = thetaPrior.GetCovariance();
      PostCov = theta.GetCovariance();

    for (int i=1; i<= ardindices.size(); i++)
      {
	PriorCov(ardindices[i],ardindices[i]) = theta.means(ardindices[i])*theta.means(ardindices[i]) + PostCov(ardindices[i],ardindices[i]);
	
	//set the Free energy contribution from ARD term
	double b = 2/(theta.means(ardindices[i])*theta.means(ardindices[i]) + PostCov(ardindices[i],ardindices[i]));
	Fard += -1.5*(log(b) + digamma(0.5)) - 0.5 - gammaln(0.5) - 0.5*log(b); //taking c as 0.5 - which it will be!
      }

	thetaPrior.SetCovariance(PriorCov);
  }

  return Fard;
}
