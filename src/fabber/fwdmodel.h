/*  fwdmodel.h - The base class for generic forward models

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


/* fwdmodel.h
 * Class declaration for generic forward models and related classes.
 * Written by Adrian Groves, 2007
 * FMRIB Centre, University of Oxford
 *
 * Last modified: $Date: 2013/04/29 12:38:19 $ $Author: chappell $ $Revision: 1.20 $
 */

//#pragma once // prevent multiple includes
#ifndef __FABBER_FWDMODEL_H
#define __FABBER_FWDMODEL_H 1

#include "assert.h"
#include "newmatap.h"
#include <string>
#include <vector>
#include "dist_mvn.h"
#include "easyoptions.h"

using namespace NEWMAT;
using namespace std;

/*
class FwdModelIdStruct {
public:
  virtual void DumpVector(const ColumnVector& vec, const string& indent = "") = 0;
  virtual ~FwdModelIdStruct() { return; }
};*/

class FwdModel {
public:
  // Virtual functions: common to all FwdModels
  
  virtual void Evaluate(const ColumnVector& params, 
			      ColumnVector& result) const = 0;
  // Evaluate the forward model

  virtual int Gradient(const ColumnVector& params, Matrix& grad) const;
  // evaluate the gradient, the int return is to indicate whether a valid gradient is returned by the model
                  
  virtual string ModelVersion() const; 
  // Return a CVS version info string
  // See fwdmodel.cc for an example of how to implement this.

  virtual int NumParams() const = 0;
  // How long should the parameter vector be?
  
  virtual int NumOutputs() const;
  // How long is output vector?  Default implementation uses Evaluate.

  // Various other useful functions:
  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const = 0;
  // Load up some sensible suggestions for initial prior & posterior values

  virtual void Initialise(MVNDist& posterior) const {};
  // voxelwise initialization of the posterior
 
  virtual void NameParams(vector<string>& names) const = 0;
  // Name each of the parameters -- see fwdmodel_linear.h for a generic implementation
  
  virtual void DumpParameters(const ColumnVector& params, 
                              const string& indent="") const;
  // Describe what a given parameter vector means (to LOG)
  // Default implementation uses NameParams to give reasonably meaningful output 
  
  
  // Static member function, to pick a forward model from a name
  static FwdModel* NewFromName(const string& name, ArgsType& args);
  
  // Usage information for this model
  static void ModelUsageFromName(const string& name, ArgsType& args);

  // An ARD update step can be specified in the model
  virtual void UpdateARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const { return; };

  // Steup function for the ARD process (forces the prior on the parameter that is subject to ARD to be correct) - really a worst case scenario if people are loading in their own priors
  virtual void SetupARD( const MVNDist& posterior, MVNDist& prior, double& Fard ) const { return; };

  //vector of indicies of parameters to which ARD should be applied;
  vector<int> ardindices;
  
  // For models that need the data values in the voxel to calculate
  virtual void pass_in_data( const ColumnVector& voxdata ) { data=voxdata; return; };
  virtual void pass_in_data( const ColumnVector& voxdata, const ColumnVector& voxsuppdata )
  { data = voxdata; suppdata = voxsuppdata; return; };

  // For models that need to know the voxel co-ordinates of the data
  virtual void pass_in_coords( const ColumnVector& coords);

  virtual ~FwdModel() { return; };
  // Virtual destructor
  
  // Your derived classes should have storage for all constants that are
  // implicitly part of g() -- e.g. pulse sequence parameters, any parameters
  // that are assumed to take known values, and basis functions.  Given these
  // constants, NumParams() should have a fixed value.

 protected:
  // storage for voxel co-ordinates
  int coord_x;
  int coord_y;
  int coord_z;
  //storage for data
  ColumnVector data;
  ColumnVector suppdata;
};

#endif /* __FABBER_FWDMODEL_H */

