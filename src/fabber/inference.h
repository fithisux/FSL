/*  inference.h - General inference technique base class

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


#if !defined(__inference_h)
#define __inference_h


#pragma once
#include <map>
#include <string>
#include <vector>
#include "fwdmodel.h"
#include "noisemodel.h"
#include "easylog.h"
#include "easyoptions.h"
#include "dataset.h"
#ifdef __FABBER_MOTION
 #include "Update_deformation.h"
#endif //__FABBER_MOTION

class InferenceTechnique {
    
 public:
  static InferenceTechnique* NewFromName(const string& method);
    // returns a *new* instance of an Inference-Technique-derived class,
    // as determined by the name given in "method".
    
 public:
  InferenceTechnique() : model(NULL), noise(NULL) { return; }
  virtual void Setup(ArgsType& args);
  virtual void SetOutputFilenames(const string& output)
    { outputDir = output; }
  virtual void DoCalculations(const DataSet& data) = 0;
  virtual void SaveResults(const DataSet& data) const;
  virtual ~InferenceTechnique();

 protected:
  FwdModel* model;
  NoiseModel* noise;
  string outputDir;
  bool saveModelFit;
  bool saveResiduals;
  
  vector<MVNDist*> resultMVNs;
  vector<MVNDist*> resultMVNsWithoutPrior; // optional; used by Adrian's spatial priors research
  vector<double> resultFs;

  void InitMVNFromFile(vector<MVNDist*>& continueFromDists,string continueFromFile, const DataSet& allData, string paramFilename);
  
  // Motion related stuff
  int Nmcstep; // number of motion correction steps to run

private:
    const InferenceTechnique& operator=(const InferenceTechnique& from)
        { assert(false); return from; } // just not allowed. 

};


// Motion Correction class
#ifdef __FABBER_MOTION
//   NB: for now the mask should cover the *entire* image as we zero everything
//       outside of the mask, which is not good for registration
//       In future we'd need allData to be able to provide the original image (or something to)

class MCobj {
public:
  MCobj(const DataSet& allData);
  void run_mc(const Matrix& modelpred_mat, Matrix& finalimage_mat);
  void set_num_iter(int nit) { num_iter=nit; }
private:
  int num_iter;  // default 10
  volume<float> mask;
  volume4D<float> defx;
  volume4D<float> defy;
  volume4D<float> defz;
  // things below are kept for efficiency (?) in order to avoid repeated allocation/destruction
  volume4D<float> tmpx;
  volume4D<float> tmpy;
  volume4D<float> tmpz;
  volume4D<float> modelpred;
  volume4D<float> finalimage;
  volume4D<float> wholeimage;
};

#endif // __FABBER_MOTION

#endif

