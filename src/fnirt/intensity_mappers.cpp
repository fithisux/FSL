// Definitions of intensity-mapper classes used by FNIRT
//
// intensity_mappers.cpp
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2007 University of Oxford 
//
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
#include <string>
#include <vector>
#include <time.h>
#include <boost/shared_ptr.hpp>
#include "newmat.h"
#include "newmatio.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "basisfield/basisfield.h"
#include "basisfield/splinefield.h"
#include "basisfield/dctfield.h"
#include "fnirtfns.h"
#include "intensity_mappers.h"

using namespace std;
using namespace boost;
using namespace NEWMAT;
using namespace NEWIMAGE;
using namespace BASISFIELD;
using namespace FNIRT;

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Member functions for class IntensityMapper
//
///////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Returns the parameters that determine the intensity mapping
//
///////////////////////////////////////////////////////////////////////////////////////////////

NEWMAT::ReturnMatrix IntensityMapper::GetPar()
const
{
  NEWMAT::ColumnVector    par(NPar());
  if (_fixed) return(par);             // If fixed, no parameters (just constants ;-)).
  unsigned int csz = 0;

  if (_sfld.size()) {
    csz = _sfld[0]->CoefSz();
    for (unsigned int i=0; i<_sfld.size(); i++) {
      par.Rows(i*csz+1,(i+1)*csz) = *(_sfld[i]->GetCoef());
    }
  }
  for (unsigned int i=0; i<_sfac.size(); i++) {
    par(_sfld.size()*csz+i+1) = _sfac[i] / _presc[i];
  }

  par.Release();
  return(par);
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Set the parameters that determine the intensity mapping
//
///////////////////////////////////////////////////////////////////////////////////////////////

void IntensityMapper::SetPar(const NEWMAT::ColumnVector&  par)
{
  if (static_cast<unsigned int>(par.Nrows()) != NPar()) throw IntensityMapperException("IntensityMapper::SetPar: Input vector par of wrong size");
  if (_fixed) throw IntensityMapperException("IntensityMapper::SetPar: Trying to set parameters of fixed intensity-mapper");

  unsigned int  csz = 0;
  if (_sfld.size()) csz = _sfld[0]->CoefSz();

  for (unsigned int i=0; i<_sfld.size(); i++) {
    _sfld[i]->SetCoef(par.Rows(i*csz+1,(i+1)*csz));
  }
  for (unsigned int i=0; i<_sfac.size(); i++) {
    _sfac[i] = _presc[i] * par(_sfld.size()*csz+i+1);
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Saves the global intensity mapping vector in NEWMAT ColumnVector style.
//
///////////////////////////////////////////////////////////////////////////////////////////////

void IntensityMapper::SaveGlobalMapping(const std::string& fname) 
const
{
  if (_sfac.size()) {
    NEWMAT::ColumnVector  tmpvec(_sfac.size());
    for (unsigned int i=0; i<_sfac.size(); i++) tmpvec(i+1) = _sfac[i];
    write_ascii_matrix(tmpvec,fname);
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Saves the local intensity mapping in image format.
//
///////////////////////////////////////////////////////////////////////////////////////////////

void IntensityMapper::SaveLocalMapping(const std::string& fname) 
const
{
  if (_sfld.size()) {
    NEWIMAGE::volume4D<float>  volfld(_sfld[0]->FieldSz_x(),_sfld[0]->FieldSz_y(),_sfld[0]->FieldSz_z(),_sfld.size());
    volfld.setxdim(_sfld[0]->Vxs_x()); volfld.setydim(_sfld[0]->Vxs_y()); volfld.setzdim(_sfld[0]->Vxs_z());
    for (unsigned int i=0; i<_sfld.size(); i++) {
      _sfld[i]->AsVolume(volfld[i]);
    }
    save_volume4D(volfld,fname); 
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Print out all global parameters (for debugging)
//
///////////////////////////////////////////////////////////////////////////////////////////////

void IntensityMapper::DebugPrint()
const
{
  if (_sfac.size()) {
    cout << "There are " << _sfac.size() << " global scaling parameters" << endl << "Global parameters = ";
    for (unsigned int i=0; i<_sfac.size(); i++) {
      cout << _sfac[i] << "  ";
    }
    cout << endl;
  }
  if (_sfld.size()) {
    cout << "There are " << _sfld.size() << " intensity fields being modelled" << endl;
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Do an intensity scaling of invol
//
///////////////////////////////////////////////////////////////////////////////////////////////

void IntensityMapper::ScaleMe(const NEWIMAGE::volume<float>&  ivol,
                              NEWIMAGE::volume<float>&        ovol)
const
{
  if (ivol.xsize() != ovol.xsize() ||
      ivol.ysize() != ovol.ysize() ||
      ivol.ysize() != ovol.ysize()) {
    throw IntensityMapperException("IntensityMapper::ScaleMe: Size mismatch between input- and output-volume.");
  }
  if (_sfld.size()) {
    if (static_cast<unsigned int>(ivol.xsize()) != _sfld[0]->FieldSz_x() ||
        static_cast<unsigned int>(ivol.ysize()) != _sfld[0]->FieldSz_y() ||
        static_cast<unsigned int>(ivol.zsize()) != _sfld[0]->FieldSz_z()) {
      throw IntensityMapperException("IntensityMapper::ScaleMe: Size mismatch between scaling-field and input volume.");
    }
  }

  // This is really effectively a switch-list
  // with all the different mapping options.

  if (_mt == NO_SCALING) {
    ovol = ivol;
  }
  else if (_mt == GLOBAL_LINEAR) {
    for (int k=0; k<ivol.zsize(); k++) {
      for (int j=0; j<ivol.ysize(); j++) {
	for (int i=0; i<ivol.xsize(); i++) {
          ovol(i,j,k) = _sfac[0] * ivol(i,j,k);
	}
      }
    }
  }
  else if (_mt == GLOBAL_NON_LINEAR || _mt == LOCAL_BIAS_WITH_GLOBAL_NON_LINEAR) {
    for (int k=0; k<ivol.zsize(); k++) {
      for (int j=0; j<ivol.ysize(); j++) {
        unsigned int bi = k*ivol.xsize()*ivol.ysize() + j*ivol.xsize();
	for (int i=0; i<ivol.xsize(); i++) {
          ovol(i,j,k) = _sfac[0];
          double prod = 1.0;
          for (unsigned int order=1; order<_sfac.size(); order++) {
            prod *= ivol(i,j,k);
	    ovol(i,j,k) += _sfac[order] * prod;
	  }
          if (_mt == LOCAL_BIAS_WITH_GLOBAL_NON_LINEAR) {
	    ovol(i,j,k) *= _sfld[0]->Peek(bi+i);
	  }
	}
      }
    }
  }
  else if (_mt == LOCAL_LINEAR) {
    for (int k=0; k<ivol.zsize(); k++) {
      for (int j=0; j<ivol.ysize(); j++) {
        unsigned int bi = k*ivol.xsize()*ivol.ysize() + j*ivol.xsize();
	for (int i=0; i<ivol.xsize(); i++) {
          ovol(i,j,k) = ivol(i,j,k) * _sfld[0]->Peek(bi+i);
	}
      }
    }
  }
  else if (_mt == LOCAL_NON_LINEAR) {
    for (int k=0; k<ivol.zsize(); k++) {
      for (int j=0; j<ivol.ysize(); j++) {
        unsigned int bi = k*ivol.xsize()*ivol.ysize() + j*ivol.xsize();
	for (int i=0; i<ivol.xsize(); i++) {
	  ovol(i,j,k) = _sfld[0]->Peek(bi+i);
          double ival = ivol(i,j,k);
          double prod = 1.0;
          for (unsigned int order=1; order<_sfld.size(); order++) {
	    prod *= ival;
            ovol(i,j,k) += _sfld[order]->Peek(bi+i) * prod;
	  }
	}
      }
    }
  }
  else {
    throw IntensityMapperException("Now, that is an error you don't see every day");
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Returns the membrane energy of the current configuration of local
// mappings. It is usually a good idea to add this to ones cost-function.
//
///////////////////////////////////////////////////////////////////////////////////////////////

double IntensityMapper::CFContribution() 
const
{
  double memen = 0.0;
  if (!Fixed()) {
    for (unsigned int i=0; i<_sfld.size(); i++) {
      memen += _lambda * _sfld[i]->BendEnergy();
      /*
      if (_sfld.size()==1) memen = _lambda * _sfld[i]->BendEnergy();
      else {
        memen += _lambda * pow(10.0,i-1.0) * _sfld[i]->BendEnergy();
      }
      */
    }
  }
  return(memen);
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Updates the fields (if any) in accordance with the new parameters.
// At the moment it assumes a splinefield, and I am not very happy
// with how I have solved it. Definitely a candidate for a clean-up.
//
///////////////////////////////////////////////////////////////////////////////////////////////

void IntensityMapper::NewSubSampling(const std::vector<unsigned int>&    ms,       // New matrix size
                                     const std::vector<double>&          vxs,      // New voxel-size
                                     const std::vector<unsigned int>&    new_ss,   // New level of sub-sampling
                                     const std::vector<unsigned int>&    old_ss)   // Old level of sub-sampling
{
  for (unsigned int i=0; i<_sfld.size(); i++) {
    BASISFIELD::splinefield& tmp_fld = dynamic_cast<BASISFIELD::splinefield &>(*(_sfld[i]));
    std::vector<unsigned int>  ksp(3,0);
    ksp[0] = tmp_fld.Ksp_x(); ksp[1] = tmp_fld.Ksp_y(); ksp[2] = tmp_fld.Ksp_z();
    for (unsigned int j=0; j<3; j++) ksp[j] = static_cast<unsigned int>((double(old_ss[j]) / double(new_ss[j])) * double(ksp[j]) + 0.5);
    _sfld[i] = _sfld[i]->ZoomField(ms,vxs);       // First step changes voxel-size
    _sfld[i] = _sfld[i]->ZoomField(ms,vxs,ksp);   // Second step changes knot-spacing 
  }
}
///////////////////////////////////////////////////////////////////////////////////////////////
//
// Returns the (empty) gradient for the case where there is no scaling.
//
///////////////////////////////////////////////////////////////////////////////////////////////

NEWMAT::ReturnMatrix IntensityMapper::Gradient(const NEWIMAGE::volume<float>&   ref,
                                               const NEWIMAGE::volume<float>&   diff,
                                               const NEWIMAGE::volume<char>     *mask)
{
  if (_mt != NO_SCALING) throw IntensityMapperException("IntensityMapper::Gradient: Must use derived class for intensity mapping");

  NEWMAT::ColumnVector   empty;
  empty.Release();
  return(empty);
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Returns the (empty) Hessian for the case where there is no scaling.
//
///////////////////////////////////////////////////////////////////////////////////////////////

boost::shared_ptr<MISCMATHS::BFMatrix> IntensityMapper::Hessian(const NEWIMAGE::volume<float>&    ref,
                                                                const NEWIMAGE::volume<char>      *mask,
                                                                MISCMATHS::BFMatrixPrecisionType  prec) 
const
{
  if (_mt != NO_SCALING) throw IntensityMapperException("IntensityMapper::Hessian: Must use derived class for intensity mapping");

  boost::shared_ptr<MISCMATHS::BFMatrix>   hess(new MISCMATHS::FullBFMatrix(0,0));  // Empty matrix
  return(hess);
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Returns the (empty) "Cross-Hessian" for the case where there is no scaling.
//
///////////////////////////////////////////////////////////////////////////////////////////////

boost::shared_ptr<MISCMATHS::BFMatrix> IntensityMapper::CrossHessian(const BASISFIELD::basisfield&        dfield,
                                                                     const NEWIMAGE::volume<float>&       dima,
                                                                     const NEWIMAGE::volume<float>&       ref,
                                                                     const NEWIMAGE::volume<char>         *mask,
                                                                     MISCMATHS::BFMatrixPrecisionType     prec) 
const
{
  if (_mt != NO_SCALING) throw IntensityMapperException("IntensityMapper::CrossHessian: Must use derived class for intensity mapping");

  boost::shared_ptr<MISCMATHS::BFMatrix>   hess(new MISCMATHS::FullBFMatrix(0,0));  // Empty matrix
  return(hess);
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Member functions for class SSDIntensityMapper
//
///////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Calculate the gradient of a sum-of-squared-differences cost-function with respect
// to the paramaters of the intensity mapping.
//
///////////////////////////////////////////////////////////////////////////////////////////////

NEWMAT::ReturnMatrix SSDIntensityMapper::Gradient(const NEWIMAGE::volume<float>&      ref,
                                                  const NEWIMAGE::volume<float>&      diff,
                                                  const NEWIMAGE::volume<char>        *mask)
{
  if (!samesize(ref,diff) || (mask && !samesize(ref,*mask))) throw IntensityMapperException("SSDIntensityMapper::Gradient: Image size mismatch");
  if (_sfld.size() && (_sfld[0]->FieldSz_x() != static_cast<unsigned int>(ref.xsize()) || 
                       _sfld[0]->FieldSz_y() != static_cast<unsigned int>(ref.ysize()) ||
                       _sfld[0]->FieldSz_z() != static_cast<unsigned int>(ref.zsize()))) {
    throw IntensityMapperException("SSDIntensityMapper::Gradient: Size mismatch between intensity field and input image");
  }

  NEWMAT::ColumnVector    grad;  // Zero-size vector
  if (Fixed()) return(grad);     // Return empty gradient if we are not estimating parameters

  // This is really effectively a switch-list
  // with all the different mapping options.

  if (_mt == NO_SCALING) {
  }                         // Intentional empty statement
  else if (_mt == GLOBAL_LINEAR) {
    grad.ReSize(1);
    grad(1) = - dotproduct(ref,diff,mask);
  }
  else if (_mt == GLOBAL_NON_LINEAR) {
    grad.ReSize(_sfac.size());
    double mv = mean(ref,mask);
    _presc[0] = mv;
    grad(1) =  - _presc[0] * sum(diff,mask);
    for (unsigned int i=1; i<_sfac.size(); i++) {
      _presc[i] = pow(mv,double(1.0-i));
      grad(i+1) =  - _presc[i] * powerdotproduct(ref,i,diff,1,mask);
    }
  }
  else if (_mt == LOCAL_LINEAR) {
    grad = - _sfld[0]->Jte(ref,diff,mask);
    grad += 0.5 * _lambda * _sfld[0]->BendEnergyGrad();
  }
  else if (_mt == LOCAL_BIAS_WITH_GLOBAL_NON_LINEAR) {
    NEWIMAGE::volume<float>                  ref_power(ref);
    NEWIMAGE::volume<float>                  wgtrefpwr(ref.xsize(),ref.ysize(),ref.zsize());
    wgtrefpwr.copyproperties(ref);
    wgtrefpwr = static_cast<float>(_sfac[0]);
    NEWIMAGE::volume<float>                  biasfield(ref.xsize(),ref.ysize(),ref.zsize());
    biasfield.copyproperties(ref);
    biasfield.insert_vec(*_sfld[0]->Get());
    grad.ReSize(_sfld[0]->CoefSz() + _sfac.size());
    grad(_sfld[0]->CoefSz()+1) = - _presc[0]*dotproduct(biasfield,diff,mask);
    for (unsigned int i=1; i<_sfac.size(); i++) {
      grad(_sfld[0]->CoefSz()+i+1) = - _presc[i]*funny_dotproduct(ref_power,1,biasfield,diff,mask);
      wgtrefpwr += ref_power * static_cast<float>(_sfac[i]);
      if (i < (_sfac.size()-1)) ref_power *= ref;
    }
    grad.Rows(1,_sfld[0]->CoefSz()) = - _sfld[0]->Jte(wgtrefpwr,diff,mask);
    grad.Rows(1,_sfld[0]->CoefSz()) += 0.5 * _lambda * _sfld[0]->BendEnergyGrad();
  }
  else if (_mt == LOCAL_NON_LINEAR) {
    grad.ReSize(_sfld.size() * _sfld[0]->CoefSz());
    NEWIMAGE::volume<float>  ref_power(ref.xsize(),ref.ysize(),ref.zsize());
    ref_power.copyproperties(ref);
    ref_power = 1.0;
    grad.Rows(1,_sfld[0]->CoefSz()) = - _sfld[0]->Jte(ref_power,diff,mask);
    grad.Rows(1,_sfld[0]->CoefSz()) += 0.5 * _lambda * _sfld[0]->BendEnergyGrad();
    //grad.Rows(1,_sfld[0]->CoefSz()) += 0.5 * _lambda * 0.1 * _sfld[0]->BendEnergyGrad();
    for (unsigned int i=1; i<_sfld.size(); i++) {
      ref_power *= ref;
      grad.Rows(i*_sfld[0]->CoefSz()+1,(i+1)*_sfld[0]->CoefSz()) = - _sfld[i]->Jte(ref_power,diff,mask);
      grad.Rows(i*_sfld[0]->CoefSz()+1,(i+1)*_sfld[0]->CoefSz()) += 0.5 * _lambda * _sfld[i]->BendEnergyGrad();
      //grad.Rows(i*_sfld[0]->CoefSz()+1,(i+1)*_sfld[0]->CoefSz()) += 0.5 * _lambda * pow(10.0,i-1.0) * _sfld[i]->BendEnergyGrad();
    }
  }

  grad.Release();
  return(grad);
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Calculate the Hessian of a sum-of-squared-differences cost-function 
// with respect to the paramaters of the intensity mapping.
//
///////////////////////////////////////////////////////////////////////////////////////////////

boost::shared_ptr<MISCMATHS::BFMatrix> SSDIntensityMapper::Hessian(const NEWIMAGE::volume<float>&    ref,
                                                                   const NEWIMAGE::volume<char>      *mask,
                                                                   MISCMATHS::BFMatrixPrecisionType  prec)
const
{
  if (mask && !samesize(ref,*mask)) throw IntensityMapperException("SSDIntensityMapper::Gradient: Image-Mask size mismatch");
  if (_sfld.size() && (_sfld[0]->FieldSz_x() != static_cast<unsigned int>(ref.xsize()) || 
                       _sfld[0]->FieldSz_y() != static_cast<unsigned int>(ref.ysize()) ||
                       _sfld[0]->FieldSz_z() != static_cast<unsigned int>(ref.zsize()))) {
    throw IntensityMapperException("SSDIntensityMapper::Gradient: Size mismatch between intensity field and input image");
  }
  
  boost::shared_ptr<MISCMATHS::BFMatrix>   hess; // Null-pointer
  if (Fixed()) {
    hess = boost::shared_ptr<MISCMATHS::BFMatrix>(new MISCMATHS::FullBFMatrix(0,0));  // Empty matrix
    return(hess);
  }

  // This is really effectively a switch-list
  // with all the different mapping options.

  if (_mt == NO_SCALING) {
    hess = boost::shared_ptr<MISCMATHS::BFMatrix>(new MISCMATHS::FullBFMatrix(0,0));  // Empty matrix
  }
  else if (_mt == GLOBAL_LINEAR) {
    double tmp = dotproduct(ref,ref,mask);
    hess = boost::shared_ptr<MISCMATHS::BFMatrix>(new MISCMATHS::FullBFMatrix(1,1));
    hess->Set(1,1,tmp);
  }
  else if (_mt == GLOBAL_NON_LINEAR) {
    NEWMAT::Matrix         FtF(_sfac.size(),_sfac.size());
    for (unsigned int i=0; i<_sfac.size(); i++) {
      for (unsigned int j=i; j<_sfac.size(); j++) {
        FtF(i+1,j+1) = _presc[i] * _presc[j] * powerdotproduct(ref,i,ref,j,mask);
        if (i != j) FtF(j+1,i+1) = FtF(i+1,j+1);
      }
    }
    hess = boost::shared_ptr<MISCMATHS::BFMatrix>(new MISCMATHS::FullBFMatrix(FtF));
  }
  else if (_mt == LOCAL_LINEAR) {
    boost::shared_ptr<MISCMATHS::BFMatrix>   reg = _sfld[0]->BendEnergyHess(prec);
    reg->MulMeByScalar(0.5 * _lambda);
    hess = _sfld[0]->JtJ(ref,mask,prec);
    hess->AddToMe(*reg);
  }
  else if (_mt == LOCAL_BIAS_WITH_GLOBAL_NON_LINEAR) {
    // FtF pertaining to global scaling factors
    NEWIMAGE::volume<float>                  biasfield(ref.xsize(),ref.ysize(),ref.zsize());
    biasfield.copyproperties(ref);
    biasfield.insert_vec(*(_sfld[0]->Get()));
    NEWMAT::Matrix         FtF(_sfac.size(),_sfac.size());
    for (unsigned int i=0; i<_sfac.size(); i++) {
      for (unsigned int j=i; j<_sfac.size(); j++) {
        FtF(i+1,j+1) = _presc[i] * _presc[j] * powerdotproduct(biasfield,2,ref,i+j,mask);
        if (i != j) FtF(j+1,i+1) = FtF(i+1,j+1);
      }
    }
    // Calculating entities needed for local and cross-part
    NEWIMAGE::volume<float>                  ref_power(ref);
    NEWIMAGE::volume<float>                  wgtrefpwr(ref.xsize(),ref.ysize(),ref.zsize());
    wgtrefpwr.copyproperties(ref);
    wgtrefpwr = static_cast<float>(_sfac[0]);
    for (unsigned int i=1; i<_sfac.size(); i++) {
      wgtrefpwr += ref_power * static_cast<float>(_sfac[i]);
      if (i < (_sfac.size()-1)) ref_power *= ref;
    }
    // Local part
    boost::shared_ptr<MISCMATHS::BFMatrix>   reg = _sfld[0]->BendEnergyHess(prec);
    reg->MulMeByScalar(0.5 * _lambda);
    hess = _sfld[0]->JtJ(wgtrefpwr,mask,prec);
    hess->AddToMe(*reg);    
    // Cross (between local and global) part
    NEWMAT::Matrix                           JtF(_sfld[0]->CoefSz(),_sfac.size());
    wgtrefpwr *= biasfield;
    JtF.Column(1) = _presc[0] * _sfld[0]->Jte(wgtrefpwr,mask);
    ref_power = ref;
    for (unsigned int i=1; i<_sfac.size(); i++) {
      JtF.Column(i+1) = _presc[i] * _sfld[0]->Jte(wgtrefpwr,ref_power,mask);
      if (i < (_sfac.size()-1)) ref_power *= ref;
    }
    // Put all the partitions together
    hess->HorConcat2MyRight(JtF);
    JtF = JtF.t() | FtF;
    hess->VertConcatBelowMe(JtF);
  }
  else if (_mt == LOCAL_NON_LINEAR) {
    std::vector<std::vector<boost::shared_ptr<MISCMATHS::BFMatrix> > >  JtJ(_sfld.size());
    NEWIMAGE::volume<float>     ref_power_i(ref.xsize(),ref.ysize(),ref.zsize());
    ref_power_i.copyproperties(ref);
    ref_power_i = 1.0;
    NEWIMAGE::volume<float>     ref_power_j(ref_power_i);
    // Fill in upper right part of vector-vector 
    // "matrix" of matrices.
    for (unsigned int i=0; i<_sfld.size(); i++) {
      if (i) ref_power_i *= ref;
      JtJ[i].resize(_sfld.size());
      ref_power_j = ref_power_i;
      for (unsigned int j=i; j<_sfld.size(); j++) {
        if (j==i) {
          boost::shared_ptr<MISCMATHS::BFMatrix>   reg = _sfld[i]->BendEnergyHess(prec);
          reg->MulMeByScalar(0.5 * _lambda);
          // reg->MulMeByScalar(0.5 * _lambda * pow(10.0,i-1.0));
          JtJ[i][j] = _sfld[i]->JtJ(ref_power_i,mask,prec);
          JtJ[i][j]->AddToMe(*reg);
	}
        else {
          ref_power_j *= ref;
          JtJ[i][j] = _sfld[0]->JtJ(ref_power_i,ref_power_j,mask,prec);
	}
      }
    }
    // Concatenate top->bottom so that the top element
    // of each column has the concatenation of all
    // matrices of that column.
    for (unsigned int col=0; col<_sfld.size(); col++) {
      for (unsigned int row=1; row<=col; row++) {
        JtJ[0][col]->VertConcatBelowMe(*JtJ[row][col]);
      }
      for (unsigned int col2=col+1; col2<_sfld.size(); col2++) {
        JtJ[0][col]->VertConcatBelowMe(*JtJ[col][col2]);
      } 
    }
    // Concatenate top row left->right to yield full Hessian
    hess = JtJ[0][0];
    for (unsigned int col=1; col<_sfld.size(); col++) {
      hess->HorConcat2MyRight(*JtJ[0][col]);
    }
  }

  return(hess);    
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Calculates the cross-terms (off-diagonal) of the Hessian
// of the cost-function. I.e. the mixed partial derivatives
// where the differentiation is w.r.t. one displacement parameter
// and one scaling parameter.
//
///////////////////////////////////////////////////////////////////////////////////////////////

boost::shared_ptr<MISCMATHS::BFMatrix> SSDIntensityMapper::CrossHessian(const BASISFIELD::basisfield&        dfield,
                                                                        const NEWIMAGE::volume<float>&       dima,
                                                                        const NEWIMAGE::volume<float>&       ref,
                                                                        const NEWIMAGE::volume<char>         *mask,
                                                                        MISCMATHS::BFMatrixPrecisionType     prec) 
const
{
  if (!samesize(dima,ref)) throw IntensityMapperException("SSDIntensityMapper::Gradient: Image-Image size mismatch");
  if (mask && !samesize(ref,*mask)) throw IntensityMapperException("SSDIntensityMapper::Gradient: Image-Mask size mismatch");
  if (_sfld.size() && (_sfld[0]->FieldSz_x() != static_cast<unsigned int>(ref.xsize()) || 
                       _sfld[0]->FieldSz_y() != static_cast<unsigned int>(ref.ysize()) ||
                       _sfld[0]->FieldSz_z() != static_cast<unsigned int>(ref.zsize()))) {
    throw IntensityMapperException("SSDIntensityMapper::Gradient: Size mismatch between intensity field and input image");
  }
  if (_sfld.size() && (_sfld[0]->FieldSz_x() != dfield.FieldSz_x() || 
                       _sfld[0]->FieldSz_y() != dfield.FieldSz_y() || 
                       _sfld[0]->FieldSz_z() != dfield.FieldSz_z())) {
    throw IntensityMapperException("SSDIntensityMapper::Hessian: Size mismatch between intensity field and displacement field");
  }

  boost::shared_ptr<MISCMATHS::BFMatrix>   hess;  // Null-pointer
  if (Fixed()) {
    hess = boost::shared_ptr<MISCMATHS::BFMatrix>(new MISCMATHS::FullBFMatrix(0,0));  // Empty matrix
    return(hess);
  }

  // This is really effectively a switch-list
  // with all the different mapping options.

  if (_mt == NO_SCALING) {
    hess = boost::shared_ptr<MISCMATHS::BFMatrix>(new MISCMATHS::FullBFMatrix(0,0));  // Empty matrix
  }
  else if (_mt == GLOBAL_LINEAR) {
    hess = boost::shared_ptr<MISCMATHS::BFMatrix>(new MISCMATHS::FullBFMatrix(-dfield.Jte(dima,ref,mask)));  // N.B. the minus sign
  }
  else if (_mt == GLOBAL_NON_LINEAR) {
    NEWMAT::Matrix           JtF(dfield.CoefSz(),_sfac.size()); 
    NEWIMAGE::volume<float>  ref_power(ref);
    JtF.Column(1) = - _presc[0] * dfield.Jte(dima,mask);                  // N.B. the minus sign
    for (unsigned int i=1; i<_sfac.size(); i++) {
      JtF.Column(i+1) = - _presc[i] * dfield.Jte(dima,ref_power,mask);    // N.B. the minus sign
      if (i < (_sfac.size()-1)) ref_power *= ref;
    }
    hess = boost::shared_ptr<MISCMATHS::BFMatrix>(new MISCMATHS::FullBFMatrix(JtF));
  }
  else if (_mt == LOCAL_LINEAR) {
    hess = dfield.JtJ(dima,*(_sfld[0]),-ref,mask,prec);    // N.B the minus sign
  }
  else if (_mt == LOCAL_BIAS_WITH_GLOBAL_NON_LINEAR) {
    NEWIMAGE::volume<float>                  biasfield(ref.xsize(),ref.ysize(),ref.zsize());
    biasfield.copyproperties(ref);
    biasfield.insert_vec(*_sfld[0]->Get());
    NEWIMAGE::volume<float>                  ref_power(ref);
    NEWIMAGE::volume<float>                  wgtrefpwr(ref.xsize(),ref.ysize(),ref.zsize());
    wgtrefpwr.copyproperties(ref);
    wgtrefpwr = _sfac[0];
    NEWMAT::Matrix                           JtF(dfield.CoefSz(),_sfac.size());
    JtF.Column(1) = - _presc[0] * dfield.Jte(dima,biasfield,mask);               // Sign
    for (unsigned int i=1; i<_sfac.size(); i++) {
      wgtrefpwr += ref_power * _sfac[i];
      JtF.Column(i+1) = - _presc[i] * dfield.Jte(dima,biasfield*ref_power,mask); // Sign
      if (i < (_sfac.size()-1)) ref_power *= ref;
    }
    hess = dfield.JtJ(dima,*(_sfld[0]),-wgtrefpwr,mask,prec);                    // Sign
    hess->HorConcat2MyRight(JtF);
  }
  else if (_mt == LOCAL_NON_LINEAR) {
    NEWIMAGE::volume<float>   ref_power(ref.xsize(),ref.ysize(),ref.zsize());
    ref_power.copyproperties(ref);
    ref_power = 1.0;
    hess = dfield.JtJ(dima,*(_sfld[0]),-ref_power,mask,prec);                        // Sign
    for (unsigned int i=1; i<_sfld.size(); i++) {
      ref_power *= ref;
      hess->HorConcat2MyRight(*(dfield.JtJ(dima,*(_sfld[i]),-ref_power,mask,prec))); // Sign
    }
  }

  return(hess);
}

// Utility function that calculates some sums of elementwise products
// of images. If you consider vol1, vol2 and vol3 below as column vectors
// then it returns ((v1.^n).*v2)'*v3

double SSDIntensityMapper::funny_dotproduct(const NEWIMAGE::volume<float>&   vol1,
                                            unsigned int                     n,
                                            const NEWIMAGE::volume<float>&   vol2,
                                            const NEWIMAGE::volume<float>&   vol3,
                                            const NEWIMAGE::volume<char>&    mask) const
{
  if (!samesize(vol1,vol2) || !samesize(vol1,vol3) || !samesize(vol1,mask)) {
    throw IntensityMapperException("SSDIntensityMapper::funny_dotproduct: Image dimension mismatch");
  }

  double rval = 0.0;
  NEWIMAGE::volume<char>::fast_const_iterator  itm = mask.fbegin();
  for (NEWIMAGE::volume<float>::fast_const_iterator it1=vol1.fbegin(), it2=vol2.fbegin(), it3=vol3.fbegin(), it_end=vol1.fend(); it1 != it_end; ++it1, ++it2, ++it3, ++itm) {
    if (*itm > 0.5) {
      double val = 1.0;
      for (unsigned int i=0; i<n; i++) val *= *it1;
      rval += static_cast<double>(val * (*it2) * (*it3));
    }
  }
  return(rval); 
}

double SSDIntensityMapper::funny_dotproduct(const NEWIMAGE::volume<float>&   vol1,
                                            unsigned int                     n,
                                            const NEWIMAGE::volume<float>&   vol2,
                                            const NEWIMAGE::volume<float>&   vol3) const
{
  if (!samesize(vol1,vol2) || !samesize(vol1,vol3)) {
    throw IntensityMapperException("SSDIntensityMapper::funny_dotproduct: Image dimension mismatch");
  }

  double rval = 0.0;
  for (NEWIMAGE::volume<float>::fast_const_iterator it1=vol1.fbegin(), it2=vol2.fbegin(), it3=vol3.fbegin(), it_end=vol1.fend(); it1 != it_end; ++it1, ++it2, ++it3) {
    double val = 1.0;
    for (unsigned int i=0; i<n; i++) val *= *it1;
    rval += static_cast<double>(val * (*it2) * (*it3));
  }
  return(rval); 
}

double SSDIntensityMapper::funny_dotproduct(const NEWIMAGE::volume<float>&   vol1,
                                            unsigned int                     n,
                                            const NEWIMAGE::volume<float>&   vol2,
                                            const NEWIMAGE::volume<float>&   vol3,
                                            const NEWIMAGE::volume<char>     *mask) const
{
  if (mask) return(funny_dotproduct(vol1,n,vol2,vol3,*mask));
  else return(funny_dotproduct(vol1,n,vol2,vol3));
}

// Utility function that returns mean of a vector, possibly ignoring zeros.

double SSDIntensityMapper::vector_mean(const NEWMAT::ColumnVector&  vec,
                                       bool                         exclude_zeros) const
{
  double       sum = double(vec.Sum());
  unsigned int n = 0;

  if (exclude_zeros) {
    for (int i=0; i<vec.Nrows(); i++) {
      if (vec(i+1)) n++;
    }
  }
  else n = static_cast<unsigned int>(vec.Nrows());

  return(sum/double(n));
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Member functions for class IntensityMapperReader
//
///////////////////////////////////////////////////////////////////////////////////////////////

vector<unsigned int> IntensityMapperReader::LocalFieldSize() const
{
  vector<unsigned int> rvec(3,0);
  if (!_has_local) throw IntensityMapperReaderException("LocalFieldSize: No local intensity mapping info available.");
  rvec[0] = static_cast<unsigned int>(_local.xsize()); 
  rvec[1] = static_cast<unsigned int>(_local.ysize()); 
  rvec[2] = static_cast<unsigned int>(_local.zsize()); 
  return(rvec);
}

vector<double> IntensityMapperReader::LocalFieldVoxelSize() const
{
  if (!_has_local) throw IntensityMapperReaderException("LocalFieldVoxelSize: No local intensity mapping info available.");
  vector<double> rvec(3,0.0);
  rvec[0] = _local.xdim(); rvec[1] = _local.ydim(); rvec[2] = _local.zdim();
  return(rvec);
}

vector<boost::shared_ptr<basisfield> > IntensityMapperReader::GetLocalAsSplinefieldVector(const vector<unsigned int>& ksp) const
{
  if (!_has_local) throw IntensityMapperReaderException("GetLocalAsSplinefieldVector: No local intensity mapping info available.");
  vector<boost::shared_ptr<basisfield> >  rvec(_local.tsize());
  vector<unsigned int>             size(3,0);
  vector<double>                   vxs(3,0.0);
  size[0] = static_cast<unsigned int>(_local.xsize());
  size[1] = static_cast<unsigned int>(_local.ysize());
  size[2] = static_cast<unsigned int>(_local.zsize());
  vxs[0] = _local.xdim(); vxs[1] = _local.ydim(); vxs[2] = _local.zdim(); 
  for (int i=0; i<_local.tsize(); i++) {
    rvec[i] = boost::shared_ptr<basisfield>(new splinefield(size,vxs,ksp));
    rvec[i]->Set(_local[i]);
  }
  return(rvec);
}

splinefield IntensityMapperReader::GetLocalAsSingleSplinefield(const vector<unsigned int>& ksp,
                                                               unsigned int                indx) const
{
  if (!_has_local) throw IntensityMapperReaderException("GetLocalAsSingleSplinefield: No local intensity mapping info available.");
  vector<unsigned int>             size(3,0);
  vector<double>                   vxs(3,0.0);
  size[0] = static_cast<unsigned int>(_local[indx].xsize());
  size[1] = static_cast<unsigned int>(_local[indx].ysize());
  size[2] = static_cast<unsigned int>(_local[indx].zsize());
  vxs[0] = _local[indx].xdim(); vxs[1] = _local[indx].ydim(); vxs[2] = _local[indx].zdim(); 
  splinefield                      field(size,vxs,ksp);
  field.Set(_local[indx]);
  return(field);
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
// File reader for internal use
//
///////////////////////////////////////////////////////////////////////////////////////////////

void IntensityMapperReader::common_read(const string& fname)
{
  // If file has been fully specified (i.e. including extension)
  // we should make sure that the file actually exists.
  if (FNIRT::extension(fname).size()) {
    if (!FNIRT::check_exist(fname)) {
      throw IntensityMapperReaderException(string("common_read: File ")+fname+string(" does not exist"));
    }
  }
  // First read .txt file, if we are asked to
  if ((FNIRT::extension(fname).size() && FNIRT::extension(fname)==string(".txt")) || 
      !FNIRT::extension(fname).size()) {
    string tmpname;
    if (!FNIRT::extension(fname).size()) tmpname = fname + string(".txt");
    else tmpname = fname;
    if (FNIRT::check_exist(tmpname)) {
      ColumnVector  tmp_vec = read_ascii_matrix(tmpname);
      _global.resize(tmp_vec.Nrows());
      for (int i=0; i<tmp_vec.Nrows(); i++) _global[i] = tmp_vec(i+1);
      _has_global = true;
      _fname = fname;
    }
  }

  // Then read image-file, if we are asked to
  if ((FNIRT::extension(fname).size() && (FNIRT::extension(fname)==string(".hdr") ||
                                          FNIRT::extension(fname)==string(".nii") ||
                                          FNIRT::extension(fname)==string(".nii.gz"))) ||
      !FNIRT::extension(fname).size()) {
    // Use MJ's read routines to check for existence
    try {
      read_volume4D(_local,fname);
      _has_local = true;
      if (!_has_global) _fname = fname;
    }
    catch(...) {
      if (!_has_global) { // If there is no image file, and there were no .txt, we consider that an error.
        throw IntensityMapperReaderException(string("No valid file matching ")+fname+string(" found"));
      }
    }
  }
}
