/*  noisemodel.h - Class declaration for generic noise models

    Adrian Groves and Michael Chappell, FMRIB Image Analysis Group

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

#pragma once
#include "dist_mvn.h"
#include "fwdmodel_linear.h"

class NoiseParams {
    /* Base class -- derived classes will store the data for NoiseModel derivatives */
public:
//    virtual void Load( const string& filename ) = 0; //??
//    virtual void Save( const string& filename ) const = 0; //??
    
    virtual NoiseParams* Clone() const = 0;
    virtual const NoiseParams& operator=(const NoiseParams& in) = 0;
    
    virtual const MVNDist OutputAsMVN() const = 0;
    virtual void InputFromMVN(const MVNDist& mvn) = 0;
       
    // Human-readable debug output (dump internal state to LOG)
    virtual void Dump(const string indent = "") const = 0;

    virtual ~NoiseParams() { return; }   
};

class NoiseModel {
    /* This class & derived classes should be essentially data-free, instead storing
     * the relevant noise parameters in a NoiseParams-derived subclass. */
     
 public:

  // Create a new identical copy of this object (e.g. for spatial vb)
//  virtual NoiseModel* Clone() const = 0;
    virtual NoiseParams* NewParams() const = 0;

  // Load priors from file, and also initialize posteriors
//  virtual void LoadPrior( const string& filename ) = 0;

    // Suggest some nice default values for noise parameters:
    virtual void HardcodedInitialDists(NoiseParams& prior, NoiseParams& posterior) const = 0;

  // Output your internal posterior distribution as an MVN.
  // (bit of a hack -- some noise models don't fit into the MVN framework well)
//  virtual const MVNDist GetResultsAsMVN() const = 0;

  // Some noise models might want to precalculate things (for efficiency 
  //  reasons), based on the length of the data... if you don't know what
  //  this is for then just ignore it.
  virtual void Precalculate( NoiseParams& noise, const NoiseParams& noisePrior, 
    const ColumnVector& sampleData ) const { return; }

  // The obligatory virtual destructor
  virtual ~NoiseModel() { return; }

  // VB Updates
  
  // The following could potentially be split into substeps; but since 
  // these would necessarily be model-specific, it's nice to have a 
  // general catch-all update step.  Presumably this function
  // would call all the other functions in some order.
  
  virtual void UpdateNoise(
    NoiseParams& noise, 
    const NoiseParams& noisePrior, 
  	const MVNDist& theta,
  	const LinearFwdModel& model,
  	const ColumnVector& data) const = 0;

  virtual void UpdateTheta(
    const NoiseParams& noise, 
//    const NoiseParams& noisePrior, 
  	MVNDist& theta,
  	const MVNDist& thetaPrior,
  	const LinearFwdModel& model,
        const ColumnVector& data,
    MVNDist* thetaWithoutPrior = NULL ,
                 // for --spatial-prior-output-correction
    float LMalpha = 0
    ) const = 0;

  virtual double CalcFreeEnergy(
    const NoiseParams& noise, 
    const NoiseParams& noisePrior, 
	const MVNDist& theta,
  	const MVNDist& thetaPrior,
  	const LinearFwdModel& model,
  	const ColumnVector& data) const = 0;

  // ARD things
  double SetupARD(vector<int> ardindices,
		  const MVNDist& theta,
		  MVNDist& thetaPrior) const;

 double UpdateARD(vector<int> ardindices,
		  const MVNDist& theta,
		  MVNDist& thetaPrior) const;

  // Potentially other functions could go here, 
  // e.g. likelihood at a point (for MCMC) or sampling function (for Gibbs)

//  virtual void SaveParams(const MVNDist& theta) { /* do nothing */ }
//  virtual void RevertParams(MVNDist& theta)
//        { throw Invalid_option("This noise model does not support reverting (don't use the trial-mode convergence detector with it)\n"); }

  // Static member function
  // If you're given a noise model name, this returns a new NoiseModel 
  // of the appropriate subclass.
  static NoiseModel* NewFromName(const string& name, ArgsType& args);

 private:
  // Prevent copying using anything other than the Clone() function.
  // Could implement it, but not particularly useful and the default
  // shallow copy is not right.
  const NoiseModel& operator=(const NoiseModel&) const
    { assert(false); return *this; } // = operator not allowed
  // don't need a private copy constructor -- abstract class.
};

// Handy mathematical function, used by some free energy calculations
double gammaln(double xx);

