/*  inference_spatialvb.h - implementation of VB with spatial priors

    Adrian Groves, FMRIB Image Analysis Group

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

#include "inference_vb.h"
#ifndef __FABBER_LIBRARYONLY
#include "newimage/newimageall.h"
#endif //__FABBER_LIBRARYONLY

class CovarianceCache {
 public:
#ifndef __FABBER_LIBRARYONLY
  void CalcDistances(const NEWIMAGE::volume<float>& mask, const string& distanceMeasure);
#endif //__FABBER_LIBRARYONLY
  void CalcDistances(const NEWMAT::Matrix& voxelCoords, const string& distanceMeasure);
  const SymmetricMatrix& GetDistances() const { return distances; }

  const ReturnMatrix GetC(double delta) const; // quick to calculate
  const SymmetricMatrix& GetCinv(double delta) const;

  //  const Matrix& GetCiCodist(double delta) const;
  const SymmetricMatrix& GetCiCodistCi(double delta, double* CiCodistTrace = NULL) const;

  bool GetCachedInRange(double *guess, double lower, double upper, bool allowEndpoints = false) const;
  // If there's a cached value in (lower, upper), set *guess = value and 
  // return true; otherwise return false and don't change *guess.

 private:
  SymmetricMatrix distances;
  typedef map<double, SymmetricMatrix> Cinv_cache_type;
  mutable Cinv_cache_type Cinv_cache; 
  
  typedef map<double, pair<SymmetricMatrix,double> > CiCodistCi_cache_type;
  //  mutable CiCodist_cache_type CiCodist_cache; // only really use the Trace
  mutable CiCodistCi_cache_type CiCodistCi_cache;
};




class SpatialVariationalBayes : public VariationalBayesInferenceTechnique {
public:
    SpatialVariationalBayes() : 
        VariationalBayesInferenceTechnique(), 
        spatialDims(-1) { return; }
    virtual void Setup(ArgsType& args); // no changes needed
    virtual void DoCalculations(const DataSet& data);
//    virtual ~SpatialVariationalBayes();

 protected:

    int spatialDims; // 0 = no spatial norm; 2 = slice only; 3 = volume
    bool continuingFromFile;

    //    bool useDataDrivenSmoothness;
    //    bool useShrinkageMethod;
    //    bool useDirichletBC;
    //    bool useMRF;
    //    bool useMRF2; // without the dirichlet bcs

    double maxPrecisionIncreasePerIteration; // Should be >1, or -1 = unlimited

    vector<vector<int> > neighbours; // Sparse matrix would be easier
    vector<vector<int> > neighbours2; // Sparse matrix would be easier
#ifndef __FABBER_LIBRARYONLY
    void CalcNeighbours(const NEWIMAGE::volume<float>& mask);
#endif //__FABBER_LIBRARYONLY
    void CalcNeighbours(const Matrix& voxelCoords);

    //vector<string> imagepriorstr; now inherited from spatialvb
    
    // For the new (Sahani-based) smoothing method:    
    CovarianceCache covar;
    string distanceMeasure;

    double fixedDelta;
    double fixedRho;
    bool updateSpatialPriorOnFirstIteration;  
    bool useEvidenceOptimization;
    double alwaysInitialDeltaGuess;
    
    bool useFullEvidenceOptimization;
    bool useSimultaneousEvidenceOptimization;
    int firstParameterForFullEO;
    bool useCovarianceMarginalsRatherThanPrecisions;
    bool keepInterparameterCovariances;

    int newDeltaEvaluations;

    string spatialPriorsTypes; // one character per parameter
    //    bool spatialPriorOutputCorrection;

    bool bruteForceDeltaSearch;

    double OptimizeSmoothingScale(
      const DiagonalMatrix& covRatio,
      //const SymmetricMatrix& covRatioSupplemented,
      const ColumnVector& meanDiffRatio, 
      double guess, double* optimizedRho = NULL, 
      bool allowRhoToVary = true,
      bool allowDeltaToVary = true) const;

    double OptimizeEvidence(
      // const vector<MVNDist>& fwdPriorVox, // used for parameters other than k
      const vector<MVNDist*>& fwdPosteriorWithoutPrior, // used for parameter k
      // const vector<SymmetricMatrix>& Si,
      int k, const MVNDist* initialFwdPrior, double guess,
      bool allowRhoToVary = false,
      double* rhoOut = NULL) const;
};


#ifdef __FABBER_LIBRARYONLY_TESTWITHNEWIMAGE
#include "newimage/newimage.h"
using namespace NEWIMAGE;
// Helper function, useful elsewhere:
void ConvertMaskToVoxelCoordinates(const volume<float>& mask, Matrix& voxelCoords);
#endif //__FABBER_LIBRARYONLY_TESTWITHNEWIMAGE
