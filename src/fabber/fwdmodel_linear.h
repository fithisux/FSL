/*  fwdmodel_linear.h - Linear forward model and related classes

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

#pragma once

#include "fwdmodel.h"

class LinearFwdModel : public FwdModel {
 public:
  // Virtual function overrides
  virtual void Evaluate(const ColumnVector& params, 
                      ColumnVector& result) const;
  virtual int NumParams() const { return centre.Nrows(); }
  virtual void DumpParameters(const ColumnVector& vec,
                              const string& indent = "") const;                            
  virtual void NameParams(vector<string>& names) const;

  ReturnMatrix Jacobian() const { return jacobian; }
  ReturnMatrix Centre() const { return centre; }
  ReturnMatrix Offset() const { return offset; }

  LinearFwdModel(const Matrix& jac, 
		 const ColumnVector& ctr, 
		 const ColumnVector& off) 
    : jacobian(jac), centre(ctr), offset(off) 
    { assert(jac.Nrows() == ctr.Ncols()); assert(jac.Ncols() == off.Ncols()); }
    
  // Upgrading to a full externally-accessible model type
  LinearFwdModel(ArgsType& args);
  virtual string ModelVersion() const;
  static void ModelUsage();
  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

 protected:
  LinearFwdModel() { return; } // Leave uninitialized; derived classes only

  Matrix jacobian;     // J (tranposed?)
  ColumnVector centre; // m
  ColumnVector offset; // g(m)
    // The amount to effectively subtract from Y is g(m)-J*m
};

class LinearizedFwdModel : public LinearFwdModel {
public:
  // Virtual function overrides
  virtual void DumpParameters(const ColumnVector& vec,
                              const string& indent = "") const;
  virtual void NameParams(vector<string>& names) const
    { assert(fcn); fcn->NameParams(names); }
  using FwdModel::ModelVersion; // Tell the compiler we want both ours and the base version
  string ModelVersion() { assert(fcn != NULL); return fcn->ModelVersion(); }

  // Constructor (leaves centre, offset and jacobian empty)
  LinearizedFwdModel(const FwdModel* model) : fcn(model) { return; }
  
  // Copy constructor (needed for using vector<LinearizedFwdModel>)
  // NOTE: This is a reference, not a pointer... and it *copies* the
  // given LinearizedFwdModel, rather than using it as its nonlinear model!
  LinearizedFwdModel(const LinearizedFwdModel& from) 
    : LinearFwdModel(from), fcn(from.fcn) { return; }

  void ReCentre(const ColumnVector& about);
  // centre=about; offset=fcn(about); 
  // jacobian = numerical differentiation about centre

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const
    { assert(fcn); fcn->HardcodedInitialDists(prior, posterior); }

  
private:
  const FwdModel* fcn;  
};

