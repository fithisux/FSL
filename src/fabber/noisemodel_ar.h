/*  noisemodel_ar.cc - Class implementation for the AR(1) noise model

    Adrian Groves and Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2007-2008 University of Oxford  */

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
#include "dist_gamma.h"
#include <vector>

using namespace std;

class Ar1cParams;

// Helper class -- caches some of the AR matrices
class Ar1cMatrixCache {
public:
  const SymmetricBandMatrix& GetMatrix(unsigned n, unsigned a12pow, 
                       unsigned a3pow) const;
  const SymmetricBandMatrix& GetMarginal(unsigned n) const;

  void Update(const Ar1cParams& dist, int nTimes);

  Ar1cMatrixCache(int numPhis) : nPhis(numPhis) { return; }
  Ar1cMatrixCache(const Ar1cMatrixCache& from)
    : alphaMarginals(from.alphaMarginals), alphaMatrices(from.alphaMatrices),
      nPhis(from.nPhis)
    { return; }

private:
  vector<SymmetricBandMatrix> alphaMarginals; 
       // recalculated whenever alpha changes
  unsigned FlattenIndex(unsigned n, unsigned a12pow, unsigned a34pow) const
    { assert(n==1 || n==2 && a12pow<=2 && a34pow<=2);
      return n-1 + 2*( a12pow + 3*(a34pow) ); } 

  vector<SymmetricBandMatrix> alphaMatrices; 
       // should only be calculated once
       // Note that if more than one model is being inferred upon at a time,
       // this will be unnecessarily duplicated in every one of them --
       // might speed things up considerably by sharing.
       
  int nPhis;
};

 
// Parameter-storage class -- it's really just an enhanced structure
class Ar1cParams : public NoiseParams {
public:
    virtual Ar1cParams* Clone() const
        { return new Ar1cParams(*this); } 

    virtual const Ar1cParams& operator=(const NoiseParams& in)
      { const Ar1cParams& from = dynamic_cast<const Ar1cParams&>(in);
	alpha = from.alpha; phis = from.phis; alphaMat = from.alphaMat; return *this; }

    virtual const MVNDist OutputAsMVN() const;
    virtual void InputFromMVN(const MVNDist& mvn);
       
    // Human-readable debug output (dump internal state to LOG)
    virtual void Dump(const string indent = "") const;

    // Constructor/destructor
    Ar1cParams(int nAlpha, int nPhi) : 
        alpha(nAlpha), phis(nPhi), alphaMat(nPhi) { return; }
    Ar1cParams(const Ar1cParams& from) : 
        alpha(from.alpha), phis(from.phis), alphaMat(from.alphaMat) { return; }
    virtual ~Ar1cParams() { return; }

private:
    friend class Ar1cNoiseModel; // Needs to use this class like it's a structure
    friend class Ar1cMatrixCache;
    MVNDist alpha;
    vector<GammaDist> phis;
    
    Ar1cMatrixCache alphaMat;
};


class Ar1cNoiseModel : public NoiseModel {
 public:

//  virtual Ar1cNoiseModel* Clone() const;
  // makes a new identical copy of this object
  
    virtual Ar1cParams* NewParams() const
        { return new Ar1cParams( NumAlphas(), nPhis ); }

    virtual void HardcodedInitialDists(NoiseParams& prior, NoiseParams& posterior) const;


//  virtual void LoadPrior( const string& filename );
  // loads priors from file, and also initializes posteriors
  
    virtual void Precalculate( NoiseParams& noise, const NoiseParams& noisePrior,
        const ColumnVector& sampleData ) const;
    // Used to pre-evaluate the alpha matrices in the cache

  // virtual void AdjustPrior(...) might be needed for multi-voxel methods...
  // probably best for that to go in a derived class. 
  
//  virtual void Dump(const string indent = "") const;
//  virtual void DumpPrior(const string indent = "") const;  
//  virtual void DumpPosterior(const string indent = "") const; 
  // human-readable debug output

//  virtual const MVNDist GetResultsAsMVN() const;

  // Constructor/destructor

    Ar1cNoiseModel(const string& ar1CrossTerms, int numPhis );
    // ar1CrossTerms must be either "none", "dual", or "same". 
    
    virtual ~Ar1cNoiseModel() { return; }
  
  // VB Updates
    
  virtual void UpdateNoise(
    NoiseParams& noise, 
    const NoiseParams& noisePrior, 
  	const MVNDist& theta,
  	const LinearFwdModel& linear,
  	const ColumnVector& data) const
  { UpdateAlpha(noise, noisePrior, theta, linear, data);
    UpdatePhi(noise, noisePrior, theta, linear, data); } 

  virtual void UpdateAlpha(
    NoiseParams& noise, 
    const NoiseParams& noisePrior, 
    const MVNDist& theta,
    const LinearFwdModel& model,
    const ColumnVector& data) const;
    
  virtual void UpdatePhi(
    NoiseParams& noise, 
    const NoiseParams& noisePrior,   
    const MVNDist& theta,
    const LinearFwdModel& model,
    const ColumnVector& data) const;

  virtual void UpdateTheta(
    const NoiseParams& noise, 
//    const NoiseParams& noisePrior,   
  	MVNDist& theta,
  	const MVNDist& thetaPrior,
  	const LinearFwdModel& model,
        const ColumnVector& data,
    MVNDist* thetaWithoutPrior = NULL,
    float LMalpha = 0
    ) const;

  virtual double CalcFreeEnergy(
    const NoiseParams& noise, 
    const NoiseParams& noisePrior,   
	const MVNDist& theta,
  	const MVNDist& thetaPrior,
  	const LinearFwdModel& model,
  	const ColumnVector& data) const;

//  void SaveParams(const MVNDist& theta) {};		 
//  void RevertParams(MVNDist& theta) {};

 protected: 
//  Ar1cParameters* prior;
//  Ar1cParameters* posterior;  
    // Whenever this changes, call alphaMat.Update!

//  Ar1cMatrixCache alphaMat;
  const string ar1Type;
  int NumAlphas() const; // converts the above string into a number
  const int nPhis;
};

