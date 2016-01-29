// FNIRT - FMRIB's Non-linear Image Registration Tool
//
// fnirt.cpp
//
// Jesper Andersson, FMRIB Image Analysis Group
//
/*    Copyright (C) 2012 University of Oxford  */

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
#include "warpfns/warpfns.h"
#include "basisfield/basisfield.h"
#include "basisfield/splinefield.h"
#include "basisfield/dctfield.h"
#include "matching_points.h"
#include "fnirt_costfunctions.h"
#include "intensity_mappers.h"
#include "fnirtfns.h"

using namespace std;
using namespace NEWMAT;
using namespace NEWIMAGE;
using namespace BASISFIELD;
using namespace FNIRT;


int main(int   argc,
         char  *argv[])
{
  // Parse command line and return object that we can 
  // later ask what the user wants us to do.
  boost::shared_ptr<fnirt_clp> clp;
  try {
    clp = parse_fnirt_command_line(argc,argv);
  }
  catch (const std::exception& error) {
    cerr << "Error occured when parsing the command line" << endl;
    cerr << "Exception thrown with message: " << error.what() << endl; 
    exit(EXIT_FAILURE);
  }

  // Read --ref and --in volumes
  boost::shared_ptr<volume<float> >     ref(new volume<float>);
  boost::shared_ptr<volume<float> >     obj(new volume<float>);
  try {
    read_volume(*ref,clp->Ref());
    read_volume(*obj,clp->Obj());
  }
  catch (const std::exception& error) {
    cerr << "Error occurred when reading --ref or --obj file" << endl;
    cerr << "Exception thrown with message: " << error.what() << endl; 
    exit(EXIT_FAILURE);
  }

  // Check for registration to self (sometimes done in TBSS)
  if (trying_to_register_to_self(clp->Ref(),*ref,clp->Obj(),*obj,clp->Affine())) {
    try {
      write_self_results(*clp,*ref);
    } 
    catch(const std::exception& error) {
      cerr << "Error occurred when writing results of registration to self" << endl;
      cerr << "Exception thrown with message: " << error.what() << endl; 
    }
    if (clp->Verbose()) cout << "Input to fnirt was two identical images. Result is unity transform." << endl;
    exit(EXIT_SUCCESS);  // Sort of
  }

  // Prepare for fnirting
  boost::shared_ptr<SSD_fnirt_CF>       cf;
  boost::shared_ptr<NonlinParam>        nlpar;
  double                                objmean = 0.0;
  try {
    // Create refmask as combination of explicit and implicit masks.
    boost::shared_ptr<volume<char> > refmask;
    if (clp->UseRefMask(1)) refmask = make_mask(clp->RefMask(),InclusiveMask,*ref,clp->UseImplicitRefMask(),clp->ImplicitRefValue());
    else refmask = make_mask(clp->RefMask(),IgnoreMask,*ref,clp->UseImplicitRefMask(),clp->ImplicitRefValue());

    // Normalise ref and obj global mean to 100
    double refmean = spmlike_mean(*ref);
    objmean = spmlike_mean(*obj);
    (*ref) *= (100.0/refmean);     
    (*obj) *= (100.0/objmean);          

    // Create objmask as a combination of explicit and implicit masks.
    boost::shared_ptr<volume<char> > objmask;
    if (clp->UseObjMask(1)) objmask = make_mask(clp->ObjMask(),InclusiveMask,*obj,clp->UseImplicitObjMask(),clp->ImplicitObjValue());
    else objmask = make_mask(clp->ObjMask(),IgnoreMask,*obj,clp->UseImplicitObjMask(),clp->ImplicitObjValue());

    // Initialise the field that describes the warps.
    std::vector<boost::shared_ptr<basisfield> >    field = init_warpfield(*clp);

    // Set up model for intensity mapping
    boost::shared_ptr<IntensityMapper>   mymapper = init_intensity_mapper(*clp);

    // Create cost-function object
    cf = boost::shared_ptr<SSD_fnirt_CF>(new SSD_fnirt_CF(*ref,*obj,clp->Affine(),field,mymapper));
    cf->SetVerbose(clp->Verbose());
    cf->SetRegularisationModel(clp->RegularisationModel());
    cf->SetLambda(clp->Lambda(1));
    cf->SetLevel(1);
    cf->SetIntensityMappingFixed(!clp->EstimateIntensity(1));
    cf->SetHessianPrecision(clp->HessianPrecision());
    cf->SetInterpolationModel(clp->InterpolationModel());
    if (clp->WeightLambdaBySSD()) cf->WeightLambdaBySSD(); 
    if (clp->UseRefDeriv()) cf->UseRefDerivs();
    if (clp->Debug()) cf->SetDebug(clp->Debug());

    // Smooth and sub-sample
    cf->SmoothRef(clp->RefFWHM(1));
    cf->SmoothObj(clp->ObjFWHM(1));
    cf->SubsampleRef(clp->SubSampling(1));

    // Set masks (if any)
    if (refmask) cf->SetRefMask(*refmask);
    if (objmask) cf->SetObjMask(*objmask);

    // Check if we should include point-lists in the warping
    if (clp->RefPointList().length()) {
      PointList         rpl(clp->RefPointList(),clp->Ref());
      rpl.SetAffine(clp->Affine().i());
      PointList         opl(clp->ObjPointList(),clp->Obj());
      MatchingPoints    mpl(rpl,opl);
      cf->SetMatchingPoints(mpl);
      cf->SetMatchingPointsLambda(clp->MatchingPointsLambda(1));
    }

    // Initialise nonlin parameter object
    NEWMAT::ColumnVector             spar = cf->Par();   // Start guesses for parameters
    nlpar = boost::shared_ptr<NonlinParam>(new NonlinParam(cf->NPar(),clp->MinimisationMethod(),spar));
    set_nlpars(*nlpar);
    nlpar->SetMaxIter(clp->MaxIter(1));
  }
  catch (const std::exception& error) {
    cerr << "Error occurred when preparing to fnirt" << endl;
    cerr << "Exception thrown with message: " << error.what() << endl; 
    exit(EXIT_FAILURE);
  }

  // And now find parameters for the first level of subsampling
  try {
    // MISCMATHS::NonlinOut status = nonlin(*nlpar,*cf);
    nonlin(*nlpar,*cf);
    if (!constrain_warpfield(*cf,*clp,5)) { // N.B. arbitrary number 5
      std::pair<double,double> range = cf->JacobianRange();
      cout << "Warning, Jacobian not within prescribed range. Prescription is " << clp->JacLowerBound() << " -- " << clp->JacUpperBound();
      cout << " and obtained range is " << range.first << " -- " << range.second << endl;
    }
  }
  catch (const std::exception& error) {
    cerr << "Error occured during estimation at first level of subsampling" << endl;
    cerr << "Exception thrown with message: " << error.what() << endl; 
    exit(EXIT_FAILURE);
  }

  // Loop over remaining levels of sub-sampling
  try {
    for (int ssl=2; ssl<=clp->NoLevels(); ssl++) {
      if (clp->Verbose()) cout << "***Going to next resolution level***" << endl;
      cf->SetLevel(ssl);

      // New mask use
      if (clp->UseRefMask(ssl) != clp->UseRefMask(ssl-1)) {
        boost::shared_ptr<volume<char> > refmask;
        if (clp->UseRefMask(ssl)) refmask = make_mask(clp->RefMask(),InclusiveMask,*ref,clp->UseImplicitRefMask(),clp->ImplicitRefValue());
        else refmask = make_mask(clp->RefMask(),IgnoreMask,*ref,clp->UseImplicitRefMask(),clp->ImplicitRefValue());
        if (refmask) cf->SetRefMask(*refmask);      
      }
      if (clp->UseObjMask(ssl) != clp->UseObjMask(ssl-1)) {
        boost::shared_ptr<volume<char> > objmask;
        if (clp->UseObjMask(1)) objmask = make_mask(clp->ObjMask(),ExclusiveMask,*obj,clp->UseImplicitObjMask(),clp->ImplicitObjValue());
        else objmask = make_mask(clp->ObjMask(),IgnoreMask,*obj,clp->UseImplicitObjMask(),clp->ImplicitObjValue());
        if (objmask) cf->SetObjMask(*objmask);
      }

      // New smoothness
      if (clp->RefFWHM(ssl) != clp->RefFWHM(ssl-1)) cf->SmoothRef(clp->RefFWHM(ssl));
      if (clp->ObjFWHM(ssl) != clp->ObjFWHM(ssl-1)) cf->SmoothObj(clp->ObjFWHM(ssl));

      // New subsampling
      if (clp->Verbose()) cout << "Setting subsampling" << endl;
      if (clp->SubSampling(ssl) != clp->SubSampling(ssl-1)) cf->SubsampleRef(clp->SubSampling(ssl));
      if (clp->Verbose()) cout << "Setting reg mode" << endl;
      cf->SetRegularisationModel(clp->RegularisationModel());
      if (clp->Verbose()) cout << "Setting lambda" << endl;
      cf->SetLambda(clp->Lambda(ssl));  // New (possibly) lambda
      cf->SetMatchingPointsLambda(clp->MatchingPointsLambda(ssl));
      cf->SetIntensityMappingFixed(!clp->EstimateIntensity(ssl));

      // Splash out on brand new nonlin object
      nlpar = boost::shared_ptr<NonlinParam>(new NonlinParam(cf->NPar(),clp->MinimisationMethod(),cf->Par()));
      set_nlpars(*nlpar);
      nlpar->SetMaxIter(clp->MaxIter(ssl));

      // MISCMATHS::NonlinOut status = nonlin(*nlpar,*cf);
      nonlin(*nlpar,*cf);
      if (!constrain_warpfield(*cf,*clp,10)) { // N.B. arbitrary number 10
        std::pair<double,double> range = cf->JacobianRange();
        cout << "Warning, Jacobian not within prescribed range. Prescription is " << clp->JacLowerBound() << " -- " << clp->JacUpperBound();
        cout << " and obtained range is " << range.first << " -- " << range.second << endl;
      }
    }

    // If we stopped short of full resolution, 
    // up-sample field to full resolution.
    if (clp->SubSampling(clp->NoLevels()) > 1) cf->SubsampleRef(1);
  }
  catch (const std::exception& error) {
    cerr << "Error occured during estimation at subsampling level > 1" << endl;
    cerr << "Exception thrown with message: " << error.what() << endl; 
    exit(EXIT_FAILURE);
  }

  // Save everything we have been asked to save
  try {
    cf->SaveDefCoefs(clp->CoefFname());                                      // Coefficients
    if (clp->FieldFname().length()) cf->SaveDefFields(clp->FieldFname());    // Field
    if (clp->JacFname().length()) cf->SaveJacobian(clp->JacFname());         // Jacobian
    if (clp->RefOutFname().length()) cf->SaveScaledRef(clp->RefOutFname());  // Intensity modulated ref scan
    // Intensity-mapping
    if (clp->IntensityMappingFname().length()) cf->SaveIntensityMapping(clp->IntensityMappingFname());
    // Warped object image
    if (clp->ObjOutFname().length()) {
      if (clp->ObjFWHM(clp->NoLevels())!=0.0) cf->SmoothObj(0.0);          // "Un-smooth" it
      cf->IntensityScaleObj(objmean/100.0);                                // "Un-scale" it
      cf->SaveRobj(clp->ObjOutFname());                                  
    }
  }
  catch(const std::exception& error) {
    cerr << "Error occured while writing output from fnirt" << endl;
    cerr << "Exception thrown with message: " << error.what() << endl; 
    exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);  // R.I.P.
}
    
