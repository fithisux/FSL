// Topup - FMRIB's Tool for correction of susceptibility induced distortions
//
// topup.cpp
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2012 University of Oxford 

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

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/shared_ptr.hpp>
#include "newmat.h"
#include "newmatio.h"
#ifndef EXPOSE_TREACHEROUS
#define EXPOSE_TREACHEROUS           // To allow us to use .set_sform etc
#endif
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/nonlin.h"
#include "utils/stack_dump.h"
#include "warpfns/warpfns.h"
#include "basisfield/basisfield.h"
#include "basisfield/splinefield.h"
#include "basisfield/dctfield.h"
#include "topup_costfunctions.h"
#include "topupfns.h"

using namespace std;
using namespace NEWMAT;
using namespace NEWIMAGE;
using namespace BASISFIELD;
using namespace TOPUP;

int main(int   argc,
         char  *argv[])
{
  StackDump::Install(); // Gives us informative stack dump if/when program crashes

  // Read command line input
  boost::shared_ptr<topup_clp> clp;
  try {
    clp = parse_topup_command_line(argc,argv);
  }
  catch (const std::exception& error) {
    cerr << "Topup: Error occured when parsing the command line" << endl;
    cerr << "Exception thrown with message: " << error.what() << endl; 
    exit(EXIT_FAILURE);
  }
  
  // Read input images
  NEWIMAGE::volume4D<float>  in;
  try {
    read_volume4D(in,clp->ImaFname());
  }
  catch (const std::exception& error) {
    cerr << "Error occurred when reading --imain file " << clp->ImaFname() << endl;
    cerr << "Exception thrown with message: " << error.what() << endl; 
    exit(EXIT_FAILURE);
  }

  // Prepare for running  

  // SCale images to obtain a consistent weighting against regularisation,
  // and possibly between scans if the scale has not been preserved.
  
  std::vector<double> means(in.tsize());
  double gmean = 0.0;
  for (int i=0; i<in.tsize(); i++) {
    means[i] = in[i].mean();
    gmean += means[i];
  }
  gmean /= in.tsize();
  if (clp->IndividualScaling()) {
    for (int i=0; i<in.tsize(); i++) in[i] *= 100.0/means[i];
  }
  else in *= 100.0/gmean;

  boost::shared_ptr<TopupCF>      cf;
  boost::shared_ptr<NonlinParam>  nlpar;
  try {
    // Create cost-function object and
    // set properties for first level

    cf = boost::shared_ptr<TopupCF>(new TopupCF(in,clp->PhaseEncodeVectors(),clp->ReadoutTimes(),clp->WarpRes(1),clp->SplineOrder()));
    cf->SetTracePrint(clp->Trace());
    cf->SetVerbose(clp->Verbose());
    cf->SetDebug(clp->DebugLevel());
    cf->SetLevel(1);
    cf->SetInterpolationModel(clp->InterpolationModel());
    cf->SetRegridding(clp->Regridding(in));
    cf->SubSample(clp->SubSampling(1));
    cf->Smooth(clp->FWHM(1));
    cf->SetMovementsFixed(!clp->EstimateMovements(1));
    cf->SetRegularisation(clp->Lambda(1),clp->RegularisationModel());
    cf->SetSSQLambda(clp->SSQLambda());
    cf->SetHessianPrecision(clp->HessianPrecision());
    
    // Create non-linear parameters object
    ColumnVector spar(cf->NPar());
    spar = 0;
    nlpar = boost::shared_ptr<NonlinParam>(new NonlinParam(cf->NPar(),clp->OptimisationMethod(1),spar));
    if (nlpar->Method() == MISCMATHS::NL_LM) {
      nlpar->SetEquationSolverMaxIter(500);
      nlpar->SetEquationSolverTol(1.0e-3);
    }
    nlpar->SetMaxIter(clp->MaxIter(1));
  }
  catch (const std::exception& error) {
    cerr << "Error occurred when preparing to run topup" << endl;
    cerr << "Exception thrown with message: " << error.what() << endl; 
    exit(EXIT_FAILURE);
  }

  // Run minimisation at first level

  try {
    nonlin(*nlpar,*cf);
  }
  catch (const std::exception& error) {
    cerr << "Error occurred when running first level of topup" << endl;
    cerr << "Exception thrown with message: " << error.what() << endl; 
    exit(EXIT_FAILURE);
  }

  // Run remaining levels to refine solution.

  unsigned int l=2;
  try {
    for ( ; l<=clp->NoOfLevels(); l++) {
      if (clp->Verbose()) cout << "***Going to next resolution level***" << endl;
      // Change settings for cost-function object
      cf->SetLevel(l);
      cf->SubSample(clp->SubSampling(l));
      cf->Smooth(clp->FWHM(l));
      cf->SetWarpResolution(clp->WarpRes(l));
      cf->SetMovementsFixed(!clp->EstimateMovements(l));
      cf->SetRegularisation(clp->Lambda(l),clp->RegularisationModel());
      // Make a new nonlinear object
      nlpar = boost::shared_ptr<NonlinParam>(new NonlinParam(cf->NPar(),clp->OptimisationMethod(l),cf->Par()));
      if (nlpar->Method() == MISCMATHS::NL_LM) {
        nlpar->SetEquationSolverMaxIter(500);
        nlpar->SetEquationSolverTol(1.0e-3);
      }
      nlpar->SetMaxIter(clp->MaxIter(l));
      try {
        nonlin(*nlpar,*cf);
      }
      catch (const std::exception& error) {
        cerr << "Error occurred when running level " << l << " of topup" << endl;
        cerr << "Exception thrown with message: " << error.what() << endl; 
        exit(EXIT_FAILURE);
      }
    }
  }
  catch (const std::exception& error) {
    cerr << "Error occurred when preparing to run level " << l << " of topup" << endl;
    cerr << "Exception thrown with message: " << error.what() << endl; 
    exit(EXIT_FAILURE);
  }

  // Save Everything we have been asked so save
  try {
    if (clp->SubSampling(clp->NoOfLevels()) > 1) {
      cf->Smooth(0.0);
      cf->SubSample(1);
    }
    cf->WriteCoefficients(clp->CoefFname());
    cf->WriteMovementParameters(clp->MovParFname());
    if (clp->FieldFname().size()) cf->WriteField(clp->FieldFname(),in);
    if (clp->ImaOutFname().size()) cf->WriteUnwarped(clp->ImaOutFname(),in,gmean/100.0);
    if (clp->DisplacementFieldBaseFname().size()) {
      cf->WriteDisplacementFields(clp->DisplacementFieldBaseFname());
    }
    if (clp->RigidBodyBaseFname().size()) {
      cf->WriteRigidBodyMatrices(clp->RigidBodyBaseFname());
    }
    if (clp->JacobianBaseFname().size()) {
      cf->WriteJacobians(clp->JacobianBaseFname());
    }
  }
  catch(const std::exception& error) {
    cerr << "Error occured while writing output from topup" << endl;
    cerr << "Exception thrown with message: " << error.what() << endl; 
    exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}
