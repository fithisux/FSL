// Definitions of cost-function classes used by topup
//
// topup_costfunctions.cpp
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2009 University of Oxford 
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
#include <ctime>
#include <cmath>
#include <boost/shared_ptr.hpp>
#include "newmat.h"
#include "newmatio.h"
#ifndef EXPOSE_TREACHEROUS
#define EXPOSE_TREACHEROUS           // To allow us to use .set_sform etc
#endif
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "warpfns/warpfns.h"
#include "basisfield/basisfield.h"
#include "basisfield/splinefield.h"
#include "basisfield/dctfield.h"
#include "warpfns/fnirt_file_reader.h"
#include "topup_file_io.h"
#include "topup_costfunctions.h"


#ifndef SQR
#define SQR(A)  (A)*(A)
#endif

using namespace TOPUP;

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class TopupScan
//
// Has the data, acquisition info and movement parameters for
// one scan. Responsible for managing the scan, resampling the
// scan (and serve it up) and keeping track of update status.
//
// {{{ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

TopupScan::TopupScan(const NEWIMAGE::volume<float>& scan, 
                     const NEWMAT::ColumnVector&    pevec, 
                     double                         rotime)
: _pevec(pevec), _rotime(rotime), _uptodate(false), _tp(false)
{
  if (_pevec.Nrows() != 3) throw TopupException("TopupScan::TopupScan: pevec must have three elements");
  if (_pevec(3) != 0.0) throw TopupException("TopupScan::TopupScan: third element of pevec must be zero");
  _orig = boost::shared_ptr<NEWIMAGE::volume<float> >(new NEWIMAGE::volume<float>(scan));
  _orig->setextrapolationmethod(NEWIMAGE::periodic);
  if (_pevec(1) && !_pevec(2)) _orig->setextrapolationvalidity(true,false,false);
  else if (!_pevec(1) && _pevec(2)) _orig->setextrapolationvalidity(false,true,false);
  else _orig->setextrapolationvalidity(false,false,false);
  _regrid = _orig;
  _subsamp = _orig;
  _smooth = _orig;
  _mp.ReSize(6);
  _mp = 0.0;
}

std::vector<unsigned int> TopupScan::ImageSize(SizeType type) const
{
  if (TracePrint()) cout << "Entering TopupScan::ImageSize" << endl;

  std::vector<unsigned int> size(3,0);
  if (type == Original) {
    size[0] = _orig->xsize(); size[1] = _orig->ysize(); size[2] = _orig->zsize();
  }
  else if (type == Regridded) {
    size[0] = _regrid->xsize(); size[1] = _regrid->ysize(); size[2] = _regrid->zsize();
  }
  else if (type == Subsampled) {
    size[0] = _subsamp->xsize(); size[1] = _subsamp->ysize(); size[2] = _subsamp->zsize();
  }
  else { // Implies Target
    unsigned int ss = _regrid->xsize() / _subsamp->xsize();
    if (_regrid->ysize()/_subsamp->ysize() != int(ss) || _regrid->zsize() / _subsamp->zsize() != int(ss)) {
      throw TopupException("TopupScan::ImageSize(Target): Inconsistent sub-sampling");
    }
    if (_orig->xsize()%ss || _orig->ysize()%ss || _orig->zsize()%ss) {
      throw TopupException("TopupScan::ImageSize(Target): Sub-sampling incompatible with original image size");
    }
    size[0] = _orig->xsize() / ss; size[1] = _orig->ysize() / ss; size[2] = _orig->zsize() / ss;
  }
  if (TracePrint()) cout << "Leaving TopupScan::ImageSize" << endl;

  return(size);
}

std::vector<double> TopupScan::ImageVxs(SizeType type) const
{
  if (TracePrint()) cout << "Entering TopupScan::ImageVxs" << endl;

  std::vector<double> vxs(3,0.0);
  if (type == Original) {
    vxs[0] = _orig->xdim(); vxs[1] = _orig->ydim(); vxs[2] = _orig->zdim();
  }
  else if (type == Regridded) {
    vxs[0] = _regrid->xdim(); vxs[1] = _regrid->ydim(); vxs[2] = _regrid->zdim();
  }
  else if (type == Subsampled) {
    vxs[0] = _subsamp->xdim(); vxs[1] = _subsamp->ydim(); vxs[2] = _subsamp->zdim();
  }
  else { // Implies Target
    unsigned int ss = _regrid->xsize() / _subsamp->xsize();
    if (_regrid->ysize()/_subsamp->ysize() != int(ss) || _regrid->zsize() / _subsamp->zsize() != int(ss)) {
      throw TopupException("TopupScan::ImageVxs(Target): Inconsistent sub-sampling");
    }
    if (_orig->xsize()%ss || _orig->ysize()%ss || _orig->zsize()%ss) {
      throw TopupException("TopupScan::ImageVxs(Target): Sub-sampling incompatible with original image size");
    }
    vxs[0] = ss * _orig->xdim(); vxs[1] = ss * _orig->ydim(); vxs[2] = ss * _orig->zdim();
  }
  if (TracePrint()) cout << "Leaving TopupScan::ImageVxs" << endl;

  return(vxs);
}

NEWIMAGE::volume<float> TopupScan::GetResampled(const BASISFIELD::splinefield& field) const
{
  if (TracePrint()) cout << "Entering TopupScan::GetResampled" << endl;

  update(field);
  if (TracePrint()) cout << "Leaving TopupScan::GetResampled" << endl;

  return(_jac * _resampled);
}

NEWIMAGE::volume<char> TopupScan::GetMask(const BASISFIELD::splinefield& field) const
{
  if (TracePrint()) cout << "Entering TopupScan::GetMask" << endl;

  update(field);
  if (TracePrint()) cout << "Leaving TopupScan::GetMask" << endl;

  return(_mask);
}

NEWIMAGE::volume<float> TopupScan::GetAlpha(const BASISFIELD::splinefield& field) const
{
  if (TracePrint()) cout << "Entering TopupScan::GetAlpha" << endl;

  update(field);
  // We will compensate for sub-sampling and change of voxel size here.
  // The derivative of the images are in units of /voxel, where voxel is the
  // re-gridded and sub-sampled space.
  std::vector<double> ovxs = ImageVxs(Original);
  std::vector<double> svxs = ImageVxs(Subsampled);
  double sf[3]; for (int i=0; i<3; i++) sf[i] = svxs[i] / ovxs[i];

  NEWIMAGE::volume<float>  rval;

  if (_pevec(1) && !_pevec(2)) rval = float(_rotime / sf[0] * _pevec(1)) * _derivs[0] * _jac;
  else if (!_pevec(1) && _pevec(2)) rval = float(_rotime / sf[1] * _pevec(2)) * _derivs[1] * _jac;
  else rval = float(_rotime) * ((float(_pevec(1)) * _derivs[0])/sf[0] + (float(_pevec(2)) * _derivs[1])/sf[1]) * _jac;

  if (TracePrint()) cout << "Leaving TopupScan::GetAlpha" << endl;

  return(rval);
}

NEWIMAGE::volume<float> TopupScan::GetBeta(const BASISFIELD::splinefield& field) const
{
  if (TracePrint()) cout << "Entering TopupScan::GetBeta" << endl;

  update(field);
  // We will compensate for sub-sampling here, though strictly speaking it really
  // belongs with the spatial derivative of the splines. Putting it here means
  // that there are fewer places where we need to do it (and maybe forget).
  // The derivative of the splines is in the voxel-space of the target.
  std::vector<double> ovxs = ImageVxs(Original);
  std::vector<double> tvxs = ImageVxs(Target);
  float ss = tvxs[0] / ovxs[0];

  NEWIMAGE::volume<float> rval;

  if (_pevec(1)) rval = float(_rotime / ss * _pevec(1)) * _resampled;
  else {
    rval = _resampled;
    rval = 0.0;
  }

  if (TracePrint()) cout << "Leaving TopupScan::GetBeta" << endl;

  return(rval);
}

NEWIMAGE::volume<float> TopupScan::GetGamma(const BASISFIELD::splinefield& field) const
{
  if (TracePrint()) cout << "Entering TopupScan::GetGamma" << endl;

  update(field);

  NEWIMAGE::volume<float> rval;
  // We will compensate for sub-sampling here, though strictly speaking it really
  // belongs with the spatial derivative of the splines. Putting it here means
  // that there are fewer places where we need to do it (and maybe forget).
  // The derivative of the splines is in the voxel-space of the target.
  std::vector<double> ovxs = ImageVxs(Original);
  std::vector<double> tvxs = ImageVxs(Target);
  float ss = tvxs[1] / ovxs[1];

  if (_pevec(2)) rval = float(_rotime / ss * _pevec(2)) * _resampled;
  else {
    rval = _resampled;
    rval = 0.0;
  }

  if (TracePrint()) cout << "Leaving TopupScan::GetGamma" << endl;

  return(rval);
}

// Returns the derivative w.r.t. the i'th movement parameter. 0...5 corresponds to
// dx, dy, dz, rx, ry, rz. It is in /mm for the translations and /radian for the rotations.
NEWIMAGE::volume<float> TopupScan::GetMovementDerivative(unsigned int i,
                                                         const BASISFIELD::splinefield& field) const
{
  if (TracePrint()) cout << "Entering TopupScan::GetMovementDerivative" << endl;

  if (i >= 6) throw TopupException("TopupScan::GetMovementDerivative: i must be in range 0...5");

  update(field);

  std::vector<unsigned int> isz = this->ImageSize(Target);
  std::vector<double> idim = this->ImageVxs(Target); 
  NEWIMAGE::volume<float>  rval; rval.reinitialize(int(isz[0]),int(isz[1]),int(isz[2]));
  rval.setdims(float(idim[0]),float(idim[1]),float(idim[2]));  

  // We need a rescaling since the derivaties are in units of /voxel for the
  // regridded/subsampled space whereas the voxel deviations below are in
  // voxels in the target space.
  std::vector<double> sdim = this->ImageVxs(Subsampled);
  double sf[3]; for (int ii=0; ii<3; ii++) sf[ii] = idim[ii]/sdim[ii];

  double tiny = 1e-4;      // If translation
  if (i > 2) tiny = 1e-5;  // If rotation
  NEWMAT::ColumnVector p(6);
  p = _mp; p(i+1) += tiny;
  NEWMAT::Matrix T1 = rval.sampling_mat().i() * mp_to_matrix(_mp).i() * rval.sampling_mat();
  NEWMAT::Matrix T2 = rval.sampling_mat().i() * mp_to_matrix(p).i() * rval.sampling_mat();

  NEWMAT::ColumnVector x(4);
  x(4) = 1.0;  
  for (int k=0; k<rval.zsize(); k++) {
    x(3) = k;
    for (int j=0; j<rval.ysize(); j++) {
      x(2) = j;
      for (int i=0; i<rval.xsize(); i++) {
        x(1) = i;
	NEWMAT::ColumnVector y1 = T1*x;
	NEWMAT::ColumnVector y2 = T2*x;
        rval(i,j,k) = float(sf[0]*(y2(1)-y1(1))*_derivs[0](i,j,k) + sf[1]*(y2(2)-y1(2))*_derivs[1](i,j,k) + sf[2]*(y2(3)-y1(3))*_derivs[2](i,j,k));
      }
    }
  }
  rval /= tiny;
  rval *= GetJacobian(field);

  if (TracePrint()) cout << "Leaving TopupScan::GetMovementDerivative" << endl;

  return(rval);
}

// Mimics the functionality of the routine above, but calculates the derivative vector
// numerically instead. Used for validation and debugging
NEWIMAGE::volume<float> TopupScan::GetNumericalMovementDerivative(unsigned int i,
                                                                  const BASISFIELD::splinefield& field) const
{
  NEWIMAGE::volume<float> resampled = GetResampled(field);
  NEWMAT::ColumnVector tmp_mp = _mp;
  double tiny = 1e-4;
  if (i > 2) tiny = 1e-4;
  tmp_mp(i+1) += tiny;
  // cout << "i = " << i << ", mp(i+1) = " << _mp(i+1) << ", tmp_mp(i+1) = " << tmp_mp(i+1) << endl;
   
  NEWMAT::Matrix rb = mp_to_matrix(tmp_mp);

  // Map field(Hz)->displacement_fields
  // Note that general_transform expects displacement fields in mm, hence the multiplication with voxel-size.
  NEWIMAGE::volume4D<float> df(_smooth->xsize(),_smooth->ysize(),_smooth->zsize(),3);
  copybasicproperties(*_smooth,df[0]); copybasicproperties(df[0],df[1]); copybasicproperties(df[0],df[2]);
  // I am forced to cast away constness on field here. When I have more
  // time I should address this in the splinefield class instead.
  BASISFIELD::splinefield& tmpfield = const_cast<BASISFIELD::splinefield& >(field);
  tmpfield.AsVolume(df[0]);
  df[1] = float(_rotime * _pevec(2) * _orig->ydim()) * df[0];
  df[0] *= _rotime * _pevec(1) * _orig->xdim(); 
  df[2] = 0; // Don't allow any component in the z-direction

  NEWIMAGE::volume<float> resampled2(_smooth->xsize(),_smooth->ysize(),_smooth->zsize());
  copybasicproperties(*_smooth,resampled2);
  resampled2 = 0.0;
  general_transform(*_smooth,rb,df,resampled2);
  resampled2 *= GetJacobian(field);
  resampled2 -= resampled;
  resampled2 /= tiny;

  return(resampled2);
}

NEWIMAGE::volume<float> TopupScan::GetJacobian(const BASISFIELD::splinefield& field) const
{
  if (TracePrint()) cout << "Entering TopupScan::GetJacobian" << endl;

  update(field);

  if (TracePrint()) cout << "Leaving TopupScan::GetJacobian" << endl;

  return(_jac);
}

void TopupScan::SetMovementParameters(const NEWMAT::ColumnVector& mp) const
{
  if (TracePrint()) cout << "Entering TopupScan::SetMovementParameters" << endl;

  if (mp.Nrows() != 6) throw TopupException("TopupScan::SetMovementParameters: mp must have six elements");
  for (unsigned int i=0; i<6; i++) {
    if (fabs(mp(i+1)-_mp(i+1))>1e-9) { _uptodate = false; _mp = mp; break; }
  } 

  if (TracePrint()) cout << "Leaving TopupScan::SetMovementParameters" << endl;
}

void TopupScan::SetInterpolationModel(TopupInterpolationType it) const
{
  if (TracePrint()) cout << "Entering TopupScan::SetInterpolationModel" << endl;

  _orig->setinterpolationmethod(translate_interp_type(it));
  _regrid->setinterpolationmethod(translate_interp_type(it));
  _subsamp->setinterpolationmethod(translate_interp_type(it));
  _smooth->setinterpolationmethod(translate_interp_type(it));
  _uptodate = false;

  if (TracePrint()) cout << "Leaving TopupScan::SetInterpolationModel" << endl;
}

void TopupScan::Smooth(double fwhm)
{
  if (TracePrint()) cout << "Entering TopupScan::Smooth" << endl;

  if (fwhm == 0.0) {
    _smooth = _subsamp;
  }
  else {
    if (_smooth == _subsamp) { // If no copy made yet
      _smooth = boost::shared_ptr<NEWIMAGE::volume<float> >(new NEWIMAGE::volume<float>(_subsamp->xsize(),_subsamp->ysize(),_subsamp->zsize()));
    }
    copybasicproperties(*_subsamp,*_smooth);
    *_smooth = NEWIMAGE::smooth(*_subsamp,fwhm/sqrt(8.0*log(2.0)));
  }
  _uptodate = false;

  if (TracePrint()) cout << "Leaving TopupScan::Smooth" << endl;
}

void TopupScan::SubSample(unsigned int ss)
{
  if (TracePrint()) cout << "Entering TopupScan::SubSample" << endl;

  if (_regrid->xsize()%ss || _regrid->ysize()%ss || _regrid->zsize()%ss) throw TopupException("TopupScan::SubSample: invalid subsampling factor");
  if (ss == 1) {
    _subsamp = _regrid;
  }
  else {
    _subsamp = boost::shared_ptr<NEWIMAGE::volume<float> >(new NEWIMAGE::volume<float>(_regrid->xsize()/ss,_regrid->ysize()/ss,_regrid->zsize()/ss));
    _subsamp->setdims(ss*_regrid->xdim(),ss*_regrid->ydim(),ss*_regrid->zdim());
    _subsamp->setinterpolationmethod(_regrid->getinterpolationmethod());
    _subsamp->setextrapolationmethod(_regrid->getextrapolationmethod());
    std::vector<bool> epval = _regrid->getextrapolationvalidity();
    _subsamp->setextrapolationvalidity(epval[0],epval[1],epval[2]);
    *(_subsamp) = 0.0;
    for (int k=0; k<_subsamp->zsize(); k++) {
      for (int j=0; j<_subsamp->ysize(); j++) {
        for (int i=0; i<_subsamp->xsize(); i++) {
          for (unsigned int kk=0; kk<ss; kk++) {
            for (unsigned int jj=0; jj<ss; jj++) {
              for (unsigned int ii=0; ii<ss; ii++) {
                (*_subsamp)(i,j,k) += (*_regrid)(i*ss+ii,j*ss+jj,k*ss+kk);
              }
            }
          }
          (*_subsamp)(i,j,k) /= double(ss*ss*ss);
        }
      }
    }
  }
  _uptodate = false;
  _smooth = _subsamp;

  if (TracePrint()) cout << "Leaving TopupScan::SubSample" << endl;
}

void TopupScan::ReGrid(int xsz, int ysz, int zsz)
{
  if (TracePrint()) cout << "Entering TopupScan::Regrid" << endl;

  if (xsz == _orig->xsize() && ysz == _orig->ysize() && zsz == _orig->zsize()) {
    _regrid = _orig;
  }
  else {
    _regrid = boost::shared_ptr<NEWIMAGE::volume<float> >(new NEWIMAGE::volume<float>(xsz,ysz,zsz));
    double xfov = (_orig->xsize()-1)*_orig->xdim() - 1e-6;
    double yfov = (_orig->ysize()-1)*_orig->ydim() - 1e-6;
    double zfov = (_orig->zsize()-1)*_orig->zdim() - 1e-6;
    _regrid->setdims(xfov/double(xsz),yfov/double(ysz),zfov/double(zsz));
    _regrid->setinterpolationmethod(_orig->getinterpolationmethod());
    _regrid->setextrapolationmethod(_orig->getextrapolationmethod());
    std::vector<bool> epval = _orig->getextrapolationvalidity();
    _regrid->setextrapolationvalidity(epval[0],epval[1],epval[2]);
    double xfac = _regrid->xdim()/_orig->xdim();
    for (int k=0; k<_regrid->zsize(); k++) {
      double z = k*_regrid->zdim()/_orig->zdim();
      for (int j=0; j<_regrid->ysize(); j++) {
	double y = j*_regrid->ydim()/_orig->ydim();
	for (int i=0; i<_regrid->xsize(); i++) {
	  (*_regrid)(i,j,k) = _orig->interpolate(i*xfac,y,z);
	}
      }
    }
    _uptodate = false;
    _subsamp = _regrid;
    _smooth = _subsamp;
  }

  if (TracePrint()) cout << "Leaving TopupScan::ReGrid" << endl;
}

NEWIMAGE::volume4D<float> TopupScan::GetDisplacementField(const BASISFIELD::splinefield&  field) const
{
  std::vector<unsigned int> isz = this->ImageSize(Target);
  std::vector<double> idim = this->ImageVxs(Target); 
  NEWIMAGE::volume4D<float> df(isz[0],isz[1],isz[2],3);
  df.setdims(float(idim[0]),float(idim[1]),float(idim[2]),1.0);  
  // I am forced to cast away constness on field here. When I have more
  // time I should address this in the splinefield class instead.
  BASISFIELD::splinefield& tmpfield = const_cast<BASISFIELD::splinefield& >(field);
  tmpfield.AsVolume(df[0]);
  df[1] = float(_rotime * _pevec(2) * _orig->ydim()) * df[0];
  df[0] *= _rotime * _pevec(1) * _orig->xdim(); 
  df[2] = 0; // Don't allow any component in the z-direction

  return(df);
}

NEWMAT::Matrix TopupScan::mp_to_matrix(const NEWMAT::ColumnVector& mp) const
{
  if (TracePrint()) cout << "Entering TopupScan::mp_to_matrix" << endl;

  NEWMAT::Matrix mat = MovePar2Matrix(mp,*_subsamp); // Defined in topup_file_io

  if (TracePrint()) cout << "Leaving TopupScan::mp_to_matrix" << endl;

  return(mat);
}

//
// This routine will update all the derived images
//
void TopupScan::update(const BASISFIELD::splinefield&  field) const
{
  if (TracePrint()) cout << "Entering TopupScan::update" << endl;

  if (!_uptodate) {
    if (TracePrint()) cout << "TopupScan::update: Things not up to date, procceeds with updating." << endl;

    // cout << "_mp = " << _mp(1) << ", " << _mp(2) << ", " << _mp(3) << ", " << _mp(4) << ", " << _mp(5) << ", " << _mp(6) << endl;
    NEWMAT::Matrix rb = mp_to_matrix(_mp);

    // Map field(Hz)->displacement_fields _for_ _the_ _non-subsampled_ _data_
    // Note that general_transform expects displacement fields in mm, hence the multiplication with voxel-size.
    std::vector<unsigned int> isz = this->ImageSize(Target);
    std::vector<double> idim = this->ImageVxs(Target); 
    NEWIMAGE::volume4D<float> df(isz[0],isz[1],isz[2],3);
    df.setdims(float(idim[0]),float(idim[1]),float(idim[2]),1.0);  
    // I am forced to cast away constness on field here. When I have more
    // time I should address this in the splinefield class instead.
    BASISFIELD::splinefield& tmpfield = const_cast<BASISFIELD::splinefield& >(field);
    tmpfield.AsVolume(df[0]);
    df[1] = float(_rotime * _pevec(2) * _orig->ydim()) * df[0];
    df[0] *= _rotime * _pevec(1) * _orig->xdim(); 
    df[2] = 0; // Don't allow any component in the z-direction

    // Get resampled data, resampled derivatives and mask
    _resampled.reinitialize(int(isz[0]),int(isz[1]),int(isz[2])); _resampled.setdims(float(idim[0]),float(idim[1]),float(idim[2])); 
    _mask.reinitialize(int(isz[0]),int(isz[1]),int(isz[2])); _mask.setdims(float(idim[0]),float(idim[1]),float(idim[2])); 
    _derivs.reinitialize(int(isz[0]),int(isz[1]),int(isz[2]),3); _derivs.setdims(float(idim[0]),float(idim[1]),float(idim[2]),1.0);  
    _resampled = 0.0; _derivs = 0.0; _mask = 0;
    general_transform_3partial(*_smooth,rb,df,_resampled,_derivs,_mask);

    // Set a "frame" in the non-phase-encode directions of
    // the mask to zero. This is to prevent small movements
    // from incurring a great cost because of a new set of 
    // voxels being excluded.
    for (int j=0; j<_mask.ysize(); j++) {
      for (int i=0; i<_mask.xsize(); i++) {
        _mask(i,j,0) = 0;
        _mask(i,j,_mask.zsize()-1) = 0;
      }
    }
    if (_pevec(1) && !_pevec(2)) {
      for (int k=0; k<_mask.zsize(); k++) {
	for (int i=0; i<_mask.xsize(); i++) {
	  _mask(i,0,k) = 0;
          _mask(i,_mask.ysize()-1,k) = 0;
	}
      }
    } 
    else if (!_pevec(1) && _pevec(2)) {
      for (int k=0; k<_mask.zsize(); k++) {
	for (int j=0; j<_mask.ysize(); j++) {
	  _mask(0,j,k) = 0;
          _mask(_mask.xsize()-1,j,k) = 0;
	}
      }
    } 

    // Get Jacobian
    BASISFIELD::splinefield xcomp = field;
    BASISFIELD::splinefield ycomp = field;
    BASISFIELD::splinefield zcomp = field;
    xcomp.ScaleField(_rotime * _pevec(1) * _orig->xdim());
    ycomp.ScaleField(_rotime * _pevec(2) * _orig->ydim());
    zcomp.ScaleField(0.0);
    _jac.reinitialize(int(isz[0]),int(isz[1]),int(isz[2])); _jac.setdims(float(idim[0]),float(idim[1]),float(idim[2])); 
    NEWIMAGE::deffield2jacobian(xcomp,ycomp,zcomp,_jac);
  
    _uptodate = true;
  }

  if (TracePrint()) cout << "Leaving TopupScan::update" << endl;
}

// }}} End of fold.

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class TopupScanManager
//
// Responsible for managing the scans, updating the resampled ones
// as neccessary, keep track of update-status and to serve up
// resampled scans and a mask that is the intersection of the
// individual masks.
//
// {{{ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

TopupScanManager::TopupScanManager(const NEWIMAGE::volume4D<float>&   scans,
                                   const NEWMAT::Matrix&              pevecs,
                                   const NEWMAT::ColumnVector&        rotimes)
  : _ss(1), _fwhm(0.0), _mpi(scans.tsize()), _tp(false)
{
  if (scans.tsize() != pevecs.Ncols() || scans.tsize() != rotimes.Nrows()) throw TopupException("TopupScanManager::TopupScanManage: Mismatched input parameters");

  _scans.resize(scans.tsize());
  for (int i=0; i<scans.tsize(); i++) {
    _scans[i] = new TopupScan(scans[i],pevecs.Column(i+1),rotimes(i+1));
  }
  copybasicproperties(scans[0],_mean);
  _mean = 0.0;
  copybasicproperties(scans[0],_mask);
  _mask = 0;
  _up_to_date = false;
  _regrid_sz = _scans[0]->ImageSize(Original);
  // Work out which movement parameters can be estimated
  _hasx=false; _hasy=false;
  double first_xvec = 0; double first_yvec = 0;
  unsigned int index_five = 0; // Index of file for which we can only estimate 5 movement parameters
  for (unsigned int i=0; i<_scans.size(); i++) {
    NEWMAT::ColumnVector pevec = _scans[i]->GetPhaseEncodeDirection();
    if (pevec(1)) { _hasx = true; if (!first_xvec) first_xvec = pevec(1); else if (have_different_signs(first_xvec,pevec(1)) && !index_five) index_five=i; }
    if (pevec(2)) { _hasy = true; if (!first_yvec) first_yvec = pevec(2); else if (have_different_signs(first_yvec,pevec(2)) && !index_five) index_five=i; }
  }
  // cout << "first_xvec = " << first_xvec << ", first_yvec = " << first_yvec << ", index_five = " << index_five << endl;
  _mpi[0].resize(0);      // First scan is always reference
  if (_hasx && _hasy) {   // If so we can estimate all paramaters
    for (unsigned int i=1; i<NoOfScans(); i++) { 
      _mpi[i].resize(6); 
      for (unsigned int j=0; j<6; j++) _mpi[i][j] = j;
    } 
  }
  else {
    for (unsigned int i=1; i<NoOfScans(); i++) { 
      if (i==index_five) {
	_mpi[i].resize(5);
	if (_hasx) { for (unsigned int j=0; j<5; j++) _mpi[i][j] = j+1; }
	else if (_hasy) { _mpi[i][0] = 0; for (unsigned int j=1; j<5; j++) _mpi[i][j] = j+1; }
      }
      else {
	_mpi[i].resize(6); 
	for (unsigned int j=0; j<6; j++) _mpi[i][j] = j;
      }
    } 
  }
}

TopupScanManager::~TopupScanManager()
{
  for (unsigned int i=0; i<_scans.size(); i++) delete _scans[i];
}

unsigned int TopupScanManager::NoOfMovementParametersForScan(unsigned int scan) const
{
  if (scan >= _scans.size()) throw TopupException("TopupScanManager::NoOfMovementParametersForScan: scan index out of range");
  return(_mpi[scan].size());
}


NEWMAT::ReturnMatrix TopupScanManager::GetMovementParameters() const
{
  if (TracePrint()) cout << "Entering TopupScanManager::GetMovementParameters" << endl;

  NEWMAT::ColumnVector ovec(NoOfMovementParameters());

  unsigned int indx = 1;
  for (unsigned int i=0; i<NoOfScans(); i++) {
    NEWMAT::ColumnVector tmp = _scans[i]->GetMovementParameters();
    for (unsigned int j=0; j<_mpi[i].size(); j++) {
      ovec(indx++) = tmp(_mpi[i][j]+1);
    }
  }

  if (TracePrint()) cout << "Leaving TopupScanManager::GetMovementParameters" << endl;

  ovec.Release();
  return(ovec);
}

NEWMAT::ReturnMatrix TopupScanManager::GetAllMovementParameters() const
{
  if (TracePrint()) cout << "Entering TopupScanManager::GetAllMovementParameters" << endl;

  NEWMAT::ColumnVector ovec(6*NoOfScans());

  for (unsigned int i=0; i<NoOfScans(); i++) {
    ovec.Rows(i*6+1,(i+1)*6) = _scans[i]->GetMovementParameters();
  }

  if (TracePrint()) cout << "Leaving TopupScanManager::GetAllMovementParameters" << endl;

  ovec.Release();
  return(ovec);
}

NEWMAT::ReturnMatrix TopupScanManager::GetRigidBodyMatrix(unsigned int i) const
{
  if (TracePrint()) cout << "Entering TopupScanManager::GetRigidBodyMatrix" << endl;

  if (i >= _scans.size()) throw TopupException("TopupScanManager::GetRigidBodyMatrix: index i out of range");

  NEWMAT::Matrix omat = _scans[i]->GetRigidBodyMatrix();

  if (TracePrint()) cout << "Leaving TopupScanManager::GetRigidBodyMatrix" << endl;

  omat.Release();
  return(omat);
}


bool TopupScanManager::HasBeta(unsigned int i) const
{
  if (TracePrint()) cout << "Entering TopupScanManager::HasBeta" << endl;

  if (i >= _scans.size()) throw TopupException("TopupScanManager::HasBeta: index i out of range");

  if (TracePrint()) cout << "Leaving TopupScanManager::HasBeta" << endl;

  return(_scans[i]->HasBeta());
}

bool TopupScanManager::HasGamma(unsigned int i) const
{
  if (TracePrint()) cout << "Entering TopupScanManager::HasGamma" << endl;

  if (i >= _scans.size()) throw TopupException("TopupScanManager::HasGamma: index i out of range");

  if (TracePrint()) cout << "Leaving TopupScanManager::HasGamma" << endl;

  return(_scans[i]->HasGamma());
}

NEWIMAGE::volume<float> TopupScanManager::GetScan(unsigned int i,
                                                  const BASISFIELD::splinefield& field) const
{
  if (TracePrint()) cout << "Entering TopupScanManager::GetScan" << endl;

  if (i >= _scans.size()) throw TopupException("TopupScanManager::GetScan: index i out of range");

  if (TracePrint()) cout << "Leaving TopupScanManager::GetScan" << endl;

  return(_scans[i]->GetResampled(field));
}

NEWIMAGE::volume<float> TopupScanManager::GetAlpha(unsigned int i,
                                                   const BASISFIELD::splinefield& field) const
{
  if (TracePrint()) cout << "Entering TopupScanManager::GetAlpha" << endl;

  if (i >= _scans.size()) throw TopupException("TopupScanManager::GetAlpha: index i out of range");

  if (TracePrint()) cout << "Leaving TopupScanManager::GetAlpha" << endl;

  return(_scans[i]->GetAlpha(field));
}

NEWIMAGE::volume<float> TopupScanManager::GetBeta(unsigned int i,
                                                  const BASISFIELD::splinefield& field) const
{
  if (TracePrint()) cout << "Entering TopupScanManager::GetBeta" << endl;

  if (i >= _scans.size()) throw TopupException("TopupScanManager::GetBeta: index i out of range");

  if (TracePrint()) cout << "Leaving TopupScanManager::GetBeta" << endl;

  return(_scans[i]->GetBeta(field));
}

NEWIMAGE::volume<float> TopupScanManager::GetGamma(unsigned int i,
                                                   const BASISFIELD::splinefield& field) const
{
  if (TracePrint()) cout << "Entering TopupScanManager::GetGamma" << endl;

  if (i >= _scans.size()) throw TopupException("TopupScanManager::GetGamma: index i out of range");

  if (TracePrint()) cout << "Leaving TopupScanManager::GetGamma" << endl;

  return(_scans[i]->GetGamma(field));
}

NEWIMAGE::volume<float> TopupScanManager::GetMovementDerivative(unsigned int scan,
                                                                unsigned int deriv,
                                                                const BASISFIELD::splinefield& field) const
{
  if (TracePrint()) cout << "Entering TopupScanManager::GetMovementDerivative" << endl;

  if (scan >= _scans.size()) throw TopupException("TopupScanManager::GetMovementDerivative: index scan out of range");
  if (deriv >= _mpi[scan].size()) throw TopupException("TopupScanManager::GetMovementDerivative: index deriv out of range");

  if (TracePrint()) cout << "Leaving TopupScanManager::GetMovementDerivative" << endl;

  return(_scans[scan]->GetMovementDerivative(_mpi[scan][deriv],field));
}

NEWIMAGE::volume<float> TopupScanManager::GetNumericalMovementDerivative(unsigned int scan,
                                                                         unsigned int deriv,
                                                                         const BASISFIELD::splinefield& field) const
{
  if (TracePrint()) cout << "Entering TopupScanManager::GetNumericalMovementDerivative" << endl;

  if (scan >= _scans.size()) throw TopupException("TopupScanManager::GetNumericalMovementDerivative: index scan out of range");
  if (deriv >= _mpi[scan].size()) throw TopupException("TopupScanManager::GetNumericalMovementDerivative: index deriv out of range");

  if (TracePrint()) cout << "Leaving TopupScanManager::GetNumericalMovementDerivative" << endl;

  // cout << "scan = " << scan << ", deriv = " << deriv << ", _mpi[scan][deriv] = " << _mpi[scan][deriv] << endl;
  return(_scans[scan]->GetNumericalMovementDerivative(_mpi[scan][deriv],field));
}

NEWIMAGE::volume<float> TopupScanManager::GetJacobian(unsigned int i,
                                                      const BASISFIELD::splinefield& field) const
{
  if (TracePrint()) cout << "Entering TopupScanManager::GetJacobian" << endl;

  if (i >= _scans.size()) throw TopupException("TopupScanManager::GetJacobian: index i out of range");

  if (TracePrint()) cout << "Leaving TopupScanManager::GetJacobian" << endl;

  return(_scans[i]->GetJacobian(field));
}

void TopupScanManager::FieldUpdated() const
{
  if (TracePrint()) cout << "Entering TopupScanManager::FieldUpdated" << endl;

  _up_to_date = false;
  for (unsigned int i=0; i<_scans.size(); i++) _scans[i]->SetUpToDate(false);

  if (TracePrint()) cout << "Leaving TopupScanManager::FieldUpdated" << endl;
}

void TopupScanManager::ReGrid(unsigned int xsz, unsigned int ysz, unsigned int zsz)
{
  if (TracePrint()) cout << "Entering TopupScanManager::ReGrid" << endl;

  if (xsz != _regrid_sz[0] || ysz != _regrid_sz[1] || zsz != _regrid_sz[2]) {
    _regrid_sz[0] = xsz; _regrid_sz[1] = ysz; _regrid_sz[2] = zsz; 
    _up_to_date = false;
    for (unsigned int i=0; i<NoOfScans(); i++) {
      _scans[i]->ReGrid(int(xsz),int(ysz),int(zsz));
      _scans[i]->SubSample(_ss);
      _scans[i]->Smooth(_fwhm);
    }    
  }
  

  if (TracePrint()) cout << "Leaving TopupScanManager::ReGrid" << endl;
}

void TopupScanManager::SubSample(unsigned int ss)
{
  if (TracePrint()) cout << "Entering TopupScanManager::SubSample" << endl;

  if (ss != _ss) {
    _ss = ss;
    _up_to_date = false;
    for (unsigned int i=0; i<NoOfScans(); i++) {
      _scans[i]->SubSample(ss);
      _scans[i]->Smooth(_fwhm);
    }
  }

  if (TracePrint()) cout << "Leaving TopupScanManager::SubSample" << endl;
}

void TopupScanManager::Smooth(double fwhm)
{
  if (TracePrint()) cout << "Entering TopupScanManager::Smooth" << endl;

  if (fwhm != _fwhm) {
    _fwhm = fwhm;
    _up_to_date = false;
    for (unsigned int i=0; i<NoOfScans(); i++) _scans[i]->Smooth(fwhm);
  }

  if (TracePrint()) cout << "Leaving TopupScanManager::Smooth" << endl;
}

void TopupScanManager::SetMovementParameters(const NEWMAT::ColumnVector& p) const
{
  if (TracePrint()) cout << "Entering TopupScanManager::SetMovementParameters" << endl;

  if (p.Nrows() != int(NoOfMovementParameters())) throw TopupException("TopupScanManager::SetMovementParameters: wrong size parameter vector p");
  unsigned int indx=1;
  for (unsigned int i=0; i<NoOfScans(); i++) {
    NEWMAT::ColumnVector tmp=p.Rows(indx,indx+_mpi[i].size()-1);
    indx += _mpi[i].size();
    if (_mpi[i].size()==6) _scans[i]->SetMovementParameters(tmp);
    else if (_mpi[i].size()) {
      NEWMAT::ColumnVector mp(6); mp = 0.0;
      for (unsigned int j=0; j<_mpi[i].size(); j++) {
	mp(_mpi[i][j]+1) = tmp(j+1);
      }
      _scans[i]->SetMovementParameters(mp);
    }
    if (!_scans[i]->UpToDate()) _up_to_date = false;
  }

  if (TracePrint()) cout << "Leaving TopupScanManager::SetMovementParameters" << endl;
}

void TopupScanManager::SetInterpolationModel(TopupInterpolationType it) const
{
  if (TracePrint()) cout << "Entering TopupScanManager::SetInterpolationModel" << endl;

  if (it != _it) {
    _up_to_date = false;
    _it = it;
    for (unsigned int i=0; i<NoOfScans(); i++) {
      _scans[i]->SetInterpolationModel(it);
    }
  }

  if (TracePrint()) cout << "Leaving TopupScanManager::SetInterpolationModel" << endl;
}

void TopupScanManager::update(const BASISFIELD::splinefield& field) const
{
  if (TracePrint()) cout << "Entering TopupScanManager::update" << endl;

  if (!_up_to_date) {
    if (TracePrint()) cout << "TopupScanManager::update: Things not up to date, proceeds with updating." << endl;

    _mean = _scans[0]->GetResampled(field);
    _mask = _scans[0]->GetMask(field);
    _mean_alpha = _scans[0]->GetAlpha(field);
    if (HasBeta()) _mean_beta = _scans[0]->GetBeta(field);
    if (HasGamma()) _mean_gamma = _scans[0]->GetGamma(field);
    for (unsigned int i=1; i<_scans.size(); i++) {
      _mean += _scans[i]->GetResampled(field);
      _mask *= _scans[i]->GetMask(field);
      _mean_alpha += _scans[i]->GetAlpha(field);
      if (HasBeta(i)) _mean_beta += _scans[i]->GetBeta(field);
      if (HasGamma(i)) _mean_gamma += _scans[i]->GetGamma(field);
    }
    _mean /= _scans.size();
    _mean_alpha /= _scans.size();
    _mean_beta /= _scans.size();
    _mean_gamma /= _scans.size();

    _up_to_date = true;
  }

  if (TracePrint()) cout << "Leaving TopupScanManager::update" << endl;
}

// }}} End of fold.

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class TopupCF
//
// The "main" class of the Topup application. Passed to the nonlin
// library.
//
// {{{ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

TopupCF::TopupCF(const NEWIMAGE::volume4D<float>&     scans,
                 const NEWMAT::Matrix&                pevecs,
                 const NEWMAT::ColumnVector&          rotimes,
                 double                               warpres,
                 unsigned int                         sporder)
  : _sm(scans,pevecs,rotimes), _field(field_factory(scans,warpres,sporder)), _wr(warpres), _lambda(10), _ssql(true), _rt(BendingEnergy), _mf(false), _hp(MISCMATHS::BFMatrixDoublePrecision), _dl(0), _level(0), _iter(0), _attempt(0)
{
  _sm.SetInterpolationModel(LinearInterp);
}

NEWMAT::ReturnMatrix TopupCF::Par() const
{
  NEWMAT::ColumnVector  ovec(NPar());

  ovec.Rows(1,NDefPar()) = *(_field.GetCoef());
  if (!_mf) ovec.Rows(NDefPar()+1,NDefPar()+NMovPar()) = _sm.GetMovementParameters();

  ovec.Release();
  return(ovec);
}

void TopupCF::SetRegridding(const std::vector<unsigned int>& rims)
{
  std::vector<unsigned int> curims = _sm.ImageSize(Regridded);
  if (rims[0] != curims[0] || rims[1] != curims[1] || rims[2] != curims[2]) {
    _sm.ReGrid(rims);
  }
}

void TopupCF::SubSample(unsigned int ss)
{
  if (TracePrint()) cout << "Entering TopupCF::SubSample" << endl;

  if (ss != _sm.SubSampling()) {
    double step = double(ss)/double(_sm.SubSampling());
    _sm.SubSample(ss); // Will throw if ss is incompatible with image size
    // Make new field
    std::vector<unsigned int> size = _sm.ImageSize(Target);
    std::vector<double> vxs = _sm.ImageVxs(Target);
    std::vector<unsigned int> ksp(3,0);
    for (unsigned int i=0; i<3; i++) ksp[i] = (static_cast<unsigned int>(MISCMATHS::round(double(_wr/vxs[i])))) ? static_cast<unsigned int>(MISCMATHS::round(double(_wr/vxs[i]))) : 1;
    BASISFIELD::splinefield new_field(size,vxs,ksp,_field.Order());
    // Insert old field into new
    NEWIMAGE::volume<float> vol(size[0],size[1],size[2]);
    vol.setdims(vxs[0],vxs[1],vxs[2]);
    double start = (step-1.0)/2.0;
    double kk=start;
    for (unsigned int k=0; k<size[2]; k++, kk+=step) {
      double jj=start;
      for (unsigned int j=0; j<size[1]; j++, jj+=step) {
        double ii=start;
	for (unsigned int i=0; i<size[0]; i++, ii+=step) {
          vol(i,j,k) = _field(ii,jj,kk);
	}
      }
    }
    new_field.Set(vol);
    _field = new_field;
  }

  if (TracePrint()) cout << "Leaving TopupCF::SubSample" << endl;
}

void TopupCF::SetWarpResolution(double wr)
{
  if (TracePrint()) cout << "Entering TopupCF::SetWarpResolution" << endl;

  if (wr != _wr) {
    _wr = wr;
    _sm.FieldUpdated();  // To notify scans they are not up to date
    // Make new field
    std::vector<unsigned int> size = _sm.ImageSize(Target);
    std::vector<double> vxs = _sm.ImageVxs(Target);
    std::vector<unsigned int> ksp(3,0);
    for (unsigned int i=0; i<3; i++) ksp[i] = (static_cast<unsigned int>(MISCMATHS::round(double(_wr/vxs[i])))) ? static_cast<unsigned int>(MISCMATHS::round(double(_wr/vxs[i]))) : 1;
    BASISFIELD::splinefield new_field(size,vxs,ksp,_field.Order());
    // Get old field in image format
    NEWIMAGE::volume<float>  vol(size[0],size[1],size[2]);
    vol.setdims(vxs[0],vxs[1],vxs[2]);
    _field.AsVolume(vol);
    new_field.Set(vol);
    _field = new_field;    
  }

  if (TracePrint()) cout << "Leaving TopupCF::SetWarpResolution" << endl;
}

void TopupCF::WriteCoefficients(const std::string& fname) const
{
  if (TracePrint()) cout << "Entering TopupCF::WriteCoefficients" << endl;

  TopupFileWriter   write_it(fname,_field);

  if (TracePrint()) cout << "Leaving TopupCF::WriteCoefficients" << endl;
}

void TopupCF::WriteUnwarped(const std::string& fname) const
{
  if (TracePrint()) cout << "Entering TopupCF::WriteUnwarped" << endl;

  std::vector<unsigned int> imsz = _sm.ImageSize(Subsampled);
  NEWIMAGE::volume4D<float> out(imsz[0],imsz[1],imsz[2],_sm.NoOfScans());
  for (unsigned int i=0; i<_sm.NoOfScans(); i++) {
    out[i] = _sm.GetScan(i,_field);
  }
  write_volume4D(out,fname);

  if (TracePrint()) cout << "Leaving TopupCF::WriteUnwarped" << endl;
}

void TopupCF::WriteUnwarped(const std::string&               fname,
			    const NEWIMAGE::volume4D<float>& hdr,
			    double                           sf) const
{
  if (TracePrint()) cout << "Entering TopupCF::WriteUnwarped" << endl;

  NEWIMAGE::volume4D<float> out = hdr;
  
  for (unsigned int i=0; i<_sm.NoOfScans(); i++) {
    out[i] = _sm.GetScan(i,_field);
    out[i] *= sf;
  }
  NEWIMAGE::volume<char> mask = _sm.GetMask(_field);
  for (int t=0; t<out.tsize(); t++) {
    for (int k=0; k<out.zsize(); k++) {
      for (int j=0; j<out.ysize(); j++) {
        for (int i=0; i<out.xsize(); i++) {
          out(i,j,k,t) = (mask(i,j,k)) ? out(i,j,k,t) : 0.0;
        }
      }
    }
  }
  
  NEWIMAGE::copybasicproperties(hdr,out);
  write_volume4D(out,fname);

  if (TracePrint()) cout << "Leaving TopupCF::WriteUnwarped" << endl;
}

void TopupCF::WriteJacobiansForDebug(const std::string& fname) const
{
  if (TracePrint()) cout << "Entering TopupCF::WriteJacobiansForDebug" << endl;

  std::vector<unsigned int> imsz = _sm.ImageSize(Subsampled);
  NEWIMAGE::volume4D<float> out(imsz[0],imsz[1],imsz[2],_sm.NoOfScans());
  for (unsigned int i=0; i<_sm.NoOfScans(); i++) {
    out[i] = _sm.GetJacobian(i,_field);
  }
  write_volume4D(out,fname);

  if (TracePrint()) cout << "Leaving TopupCF::WriteJacobiansForDebug" << endl;
}

void TopupCF::WriteJacobians(const std::string& fname) const
{
  if (TracePrint()) cout << "Entering TopupCF::WriteJacobians" << endl;

  for (unsigned int i=0; i<_sm.NoOfScans(); i++) {
    NEWIMAGE::volume<float> out = _sm.GetJacobian(i,_field);
    std::stringstream tmp(fname+std::string("_"),ios_base::out|ios_base::app);
    tmp.width(2); tmp.fill('0'); tmp << i+1;
    write_volume(out,tmp.str());
  }

  if (TracePrint()) cout << "Leaving TopupCF::WriteJacobians" << endl;
}

void TopupCF::WriteField(const std::string& fname) const
{
  if (TracePrint()) cout << "Entering TopupCF::WriteField" << endl;

  TopupFileWriter  write_it(fname,_sm.GetScan(0,_field),_field);

  if (TracePrint()) cout << "Leaving TopupCF::WriteField" << endl;
}

void TopupCF::WriteField(const std::string&               fname,
			 const NEWIMAGE::volume4D<float>& hdr) const
{
  if (TracePrint()) cout << "Entering TopupCF::WriteField" << endl;

  TopupFileWriter  write_it(fname,hdr[0],_field);

  if (TracePrint()) cout << "Leaving TopupCF::WriteField" << endl;
}

void TopupCF::WriteDisplacementFields(const std::string& fname) const
{
  if (TracePrint()) cout << "Entering TopupCF::WriteDisplacementFields" << endl;

  for (unsigned int i=0; i<_sm.NoOfScans(); i++) {
    NEWIMAGE::volume4D<float> out = _sm.GetDisplacementField(i,_field);
    std::stringstream tmp(fname+std::string("_"),ios_base::out|ios_base::app);
    tmp.width(2); tmp.fill('0'); tmp << i+1;
    write_volume4D(out,tmp.str());
  }

  if (TracePrint()) cout << "Leaving TopupCF::WriteDisplacementFields" << endl;
}

void TopupCF::WriteMask(const std::string& fname) const
{
  if (TracePrint()) cout << "Entering TopupCF::WriteMask" << endl;

  const NEWIMAGE::volume<char>&    mask = _sm.GetMask(_field);
  write_volume(mask,fname);

  if (TracePrint()) cout << "Leaving TopupCF::WriteMask" << endl;
}

void TopupCF::WriteMaskedDiff(const std::string& fname) const
{
  if (TracePrint()) cout << "Entering TopupCF::WriteMaskedDiff" << endl;

  std::vector<unsigned int>        imsz = _sm.ImageSize(Subsampled);
  NEWIMAGE::volume4D<float>        out(imsz[0],imsz[1],imsz[2],_sm.NoOfScans());
  const NEWIMAGE::volume<float>&   mean = _sm.GetMean(_field);
  const NEWIMAGE::volume<char>&    mask = _sm.GetMask(_field);
  for (unsigned int s=0; s<_sm.NoOfScans(); s++) {
    out[s] = 0.0;
    const NEWIMAGE::volume<float>& scan = _sm.GetScan(s,_field);
    for (int k=0; k<mean.zsize(); k++) {
      for (int j=0; j<mean.ysize(); j++) {
	for (int i=0; i<mean.xsize(); i++) {
	  if (mask(i,j,k)) {
	    out[s](i,j,k) = mean(i,j,k)-scan(i,j,k);
	  }
	}
      }
    }
  }
  write_volume4D(out,fname);

  if (TracePrint()) cout << "Leaving TopupCF::WriteMaskedDiff" << endl;
}

void TopupCF::WriteMovementParameters(const std::string& fname) const
{
  if (TracePrint()) cout << "Entering TopupCF::WriteMovementParameters" << endl;

  NEWMAT::ColumnVector   mp = _sm.GetAllMovementParameters();
  NEWMAT::Matrix         out(_sm.NoOfScans(),6);
  for (unsigned int s=0; s<_sm.NoOfScans(); s++) {
    for (unsigned int p=0; p<6; p++) {
      out(s+1,p+1) = mp(s*6+p+1);
    }
  }
  
  TopupFileWriter   write_it(fname,out);

  if (TracePrint()) cout << "Leaving TopupCF::WriteMovementParameters" << endl;
}

void TopupCF::WriteRigidBodyMatrices(const std::string& fname) const
{
  if (TracePrint()) cout << "Entering TopupCF::WriteRigidBodyMatrices" << endl;

  for (unsigned s=0; s<_sm.NoOfScans(); s++) {
    NEWMAT::Matrix M = _sm.GetRigidBodyMatrix(s);
    std::stringstream tmp(fname+std::string("_"),ios_base::out|ios_base::app);
    tmp.width(2); tmp.fill('0'); tmp << s+1 << ".mat";
    MISCMATHS::write_ascii_matrix(tmp.str(),M);
  }

  if (TracePrint()) cout << "Leaving TopupCF::WriteRigidBodyMatrices" << endl;
}

double TopupCF::FieldEnergy() const
{
  if (TracePrint()) cout << "Entering TopupCF::FieldEnergy" << endl;

  if (_rt==BendingEnergy) return(_field.BendEnergy());
  else return(_field.MemEnergy());

  if (TracePrint()) cout << "Leaving TopupCF::FieldEnergy" << endl;
}

NEWMAT::ReturnMatrix TopupCF::FieldEnergyGrad() const
{
  if (TracePrint()) cout << "Entering TopupCF::FieldEnergyGrad" << endl;

  NEWMAT::ColumnVector rval;

  if (_rt==BendingEnergy) rval = _field.BendEnergyGrad();
  else rval = _field.MemEnergyGrad();

  if (TracePrint()) cout << "Leaving TopupCF::FieldEnergyGrad" << endl;

  rval.Release();
  return(rval);  
}

boost::shared_ptr<MISCMATHS::BFMatrix> TopupCF::FieldEnergyHess() const
{
  if (TracePrint()) cout << "Entering TopupCF::FieldEnergyHess" << endl;

  boost::shared_ptr<MISCMATHS::BFMatrix> hptr;
  if (_rt==BendingEnergy) hptr = _field.BendEnergyHess(_hp);
  else hptr = _field.MemEnergyHess(_hp);
    
  if (TracePrint()) cout << "Leaving TopupCF::FieldEnergyHess" << endl;

  return(hptr);
}

double TopupCF::cf(const NEWMAT::ColumnVector& p) const
{
  if (TracePrint()) cout << "Entering TopupCF::cf" << endl;
 
  if (p.Nrows() != int(NPar())) throw TopupException("Topup_CF::cf: wrong size parameter vector p");

  set_field_params(p);
  set_movement_params(p);

  const NEWIMAGE::volume<float>&   mean = _sm.GetMean(_field);
  const NEWIMAGE::volume<char>&    mask = _sm.GetMask(_field);

  double ssd = 0.0;
  unsigned int n = static_cast<unsigned int>(mask.sum());
  for (unsigned int s=0; s<_sm.NoOfScans(); s++) {
    const NEWIMAGE::volume<float>&  scan = _sm.GetScan(s,_field);
    for (int k=0; k<mean.zsize(); k++) {
      for (int j=0; j<mean.ysize(); j++) {
	for (int i=0; i<mean.xsize(); i++) {
	  if (mask(i,j,k)) {
	    ssd += SQR(mean(i,j,k)-scan(i,j,k));
	  }
	}
      }
    }
  }
  ssd /= double(n*(_sm.NoOfScans()-1));
  set_latest_ssd(ssd);

  if (debug_level()) {
    set_attempt(attempt()+1);
    WriteUnwarped(string("TopupDebugUnwarped")+debug_string());    
    WriteJacobiansForDebug(string("TopupDebugJacobians")+debug_string());   
    WriteField(string("TopupDebugField")+debug_string());
    WriteMovementParameters(string("TopupDebugMovementParameters")+debug_string()+string(".txt"));
    if (debug_level() > 1) {
      WriteMask(string("TopupDebugMask")+debug_string());
      WriteMaskedDiff(string("TopupDebugMaskedDiff")+debug_string());   
    }
  }

  if (Verbose()) cout << "SSD = " << ssd << "\tn = " << n;

  double reg = Lambda() * FieldEnergy();  // May be membrane- or bending-energy.
  if (Verbose()) cout << "\tReg = " << reg;

  double cost = ssd + reg;
  if (Verbose()) cout << "\tCost = " << cost << endl;

  if (TracePrint()) cout << "Leaving TopupCF::cf" << endl;

  return(cost);
}

NEWMAT::ReturnMatrix TopupCF::grad(const NEWMAT::ColumnVector& p) const
{
  if (TracePrint()) cout << "Entering TopupCF::grad" << endl;

  if (p.Nrows() != int(NPar())) throw TopupException("Topup_CF::grad: wrong size parameter vector p");

  set_field_params(p);
  if (!MovementsFixed()) set_movement_params(p);
  NEWMAT::ColumnVector gradient(NPar());
  gradient = 0.0;

  // The top part of the gradient pertains to the non-linear displacements

  const NEWIMAGE::volume<float>&   mean = _sm.GetMean(_field);
  const NEWIMAGE::volume<float>&   mean_alpha = _sm.GetMeanAlpha(_field);
  const NEWIMAGE::volume<float>&   mean_beta = _sm.GetMeanBeta(_field);
  const NEWIMAGE::volume<float>&   mean_gamma = _sm.GetMeanGamma(_field);
  const NEWIMAGE::volume<char>&    mask = _sm.GetMask(_field);
  unsigned int n = static_cast<unsigned int>(mask.sum());
  unsigned int m = _sm.NoOfScans();

  // Here we calculate the entities used in equation 18 in the tech-report
  //
  // First Allocate memory for the ones we will use
  std::vector<NEWIMAGE::volume<float> *> abf(3,static_cast<NEWIMAGE::volume<float> *>(0));
  abf[0] = new NEWIMAGE::volume<float>(mean.xsize(),mean.ysize(),mean.zsize());
  copybasicproperties(mean,*(abf[0])); *(abf[0]) = 0.0;
  if (_sm.HasBeta()) {
    abf[1] = new NEWIMAGE::volume<float>(mean.xsize(),mean.ysize(),mean.zsize());
    copybasicproperties(mean,*(abf[1])); *(abf[1]) = 0.0;
  }
  if (_sm.HasGamma()) {
    abf[2] = new NEWIMAGE::volume<float>(mean.xsize(),mean.ysize(),mean.zsize());
    copybasicproperties(mean,*(abf[2])); *(abf[2]) = 0.0;
  }

  // Then calculate them

  for (unsigned int i=0; i<_sm.NoOfScans(); i++) {
    NEWIMAGE::volume<float>  diff = mean - _sm.GetScan(i,_field);    
    *(abf[0]) += (mean_alpha - _sm.GetAlpha(i,_field)) * diff;
    if (_sm.HasBeta(i)) *(abf[1]) += (mean_beta - _sm.GetBeta(i,_field)) * diff;
    if (_sm.HasGamma(i)) *(abf[2]) += (mean_gamma - _sm.GetGamma(i,_field)) * diff;
  }
  gradient.Rows(1,NDefPar()) += 2.0 * _field.Jte(*(abf[0]),&mask) / (n*(m - 1));
  std::vector<unsigned int> deriv(3,0);
  if (_sm.HasBeta()) { deriv[0] = 1; gradient.Rows(1,NDefPar()) += 2.0 * _field.Jte(deriv,*(abf[1]),&mask) / (n*(m - 1)); deriv[0] = 0; }
  if (_sm.HasGamma()) { deriv[1] = 1; gradient.Rows(1,NDefPar()) += 2.0 * _field.Jte(deriv,*(abf[2]),&mask) / (n*(m - 1)); deriv[1] = 0; }

  delete abf[0];
  if (_sm.HasBeta()) delete abf[1];
  if (_sm.HasGamma()) delete abf[2];

  // Now do the movement bit

  if (!MovementsFixed()) {
    unsigned int gi = NDefPar()+1;
    for (unsigned int s=0; s<_sm.NoOfScans(); s++) {
      NEWIMAGE::volume<float>  diff = mean - _sm.GetScan(s,_field);
      // cout << "Scan s = " << s << ", No d = " << _sm.NoOfMovementParametersForScan(s) << endl;
      for (unsigned int d=0; d<_sm.NoOfMovementParametersForScan(s); d++) {
	NEWIMAGE::volume<float> deriv = _sm.GetMovementDerivative(s,d,_field);
        gradient(gi) = - (2.0 / (n*(m - 1))) * sum_of_prod(diff,deriv,mask);
        gi++;
      }
    }
    // gradient.Rows(NDefPar()+1,NPar()) = numerical_gradient(p,NDefPar()+1,NPar(),1e-4,false); // Uncomment for numerical derivatives
  }
 
  // Look at numerical and analytical derivatives

  /*
  static int counter=0;
  cout << "counter = " << counter << endl;
  if (counter == 14) {
    
    NEWMAT::ColumnVector comp_grad(NPar()-NDefPar());
    unsigned int gi = 1;
    for (unsigned int s=0; s<_sm.NoOfScans(); s++) {
      NEWIMAGE::volume<float>  diff = mean - _sm.GetScan(s,_field);
      cout << "Scan s = " << s << ", No d = " << _sm.NoOfMovementParametersForScan(s) << endl;
      for (unsigned int d=0; d<_sm.NoOfMovementParametersForScan(s); d++) {
	NEWIMAGE::volume<float> aderiv = _sm.GetMovementDerivative(s,d,_field);
	NEWIMAGE::volume<float> nderiv = _sm.GetNumericalMovementDerivative(s,d,_field);
        char fname[256]; sprintf(fname,"Analytical_derivative_image_s%d_d%d",s+1,d+1);
        write_volume(aderiv,string(fname));
        sprintf(fname,"Numerical_derivative_image_s%d_d%d",s+1,d+1);
        write_volume(nderiv,string(fname));
        comp_grad(gi) = - (2.0 / (n*(m - 1))) * sum_of_prod(diff,nderiv,mask);
        gi++;
      }
    }
    
    // Derivatives w.r.t. the warps
    MISCMATHS::write_ascii_matrix(string("numerical_warp_derivatives.txt"),numerical_gradient(p,NDefPar()/2,NDefPar()/2+99,1e-1,false));
    MISCMATHS::write_ascii_matrix(string("analytical_warp_derivatives.txt"),gradient.Rows(NDefPar()/2,NDefPar()/2+99));
    // Movement derivatives
    MISCMATHS::write_ascii_matrix(string("numerical_movement_derivatives.txt"),numerical_gradient(p,NDefPar()+1,NPar(),1e-4,false));
    MISCMATHS::write_ascii_matrix(string("analytical_movement_derivatives.txt"),gradient.Rows(NDefPar()+1,NPar()));
    // MISCMATHS::write_ascii_matrix(string("composite_movement_derivatives.txt"),comp_grad);
    exit(EXIT_SUCCESS);
  }
  else counter++;
  */

  if (debug_level()) {
    set_iter(iter()+1);
    if (debug_level() > 0) {
      MISCMATHS::write_ascii_matrix(string("TopupDebugGradient_")+debug_string()+string(".txt"),gradient);
    }
  }

  // Add regularisation to the non-linear bit

  gradient.Rows(1,NDefPar()) += Lambda()*FieldEnergyGrad();

  if (debug_level() > 1) {
    MISCMATHS::write_ascii_matrix(string("TopupDebugGradientWithReg_")+debug_string()+string(".txt"),gradient);
  }

  if (TracePrint()) cout << "Leaving TopupCF::grad" << endl;
  
  gradient.Release();
  return(gradient);
}

// The resulting hessian is organised as
// [d2f/dwdw d2f/dwdp]
// [d2f/dpdw d2f/dpdp]

boost::shared_ptr<MISCMATHS::BFMatrix> TopupCF::hess(const NEWMAT::ColumnVector&             p,
                                                     boost::shared_ptr<MISCMATHS::BFMatrix>  iptr) const
{
  if (TracePrint()) cout << "Entering TopupCF::hess" << endl;

  if (p.Nrows() != int(NPar())) throw TopupException("Topup_CF::grad: wrong size parameter vector p");

  set_field_params(p);
  if (!MovementsFixed()) set_movement_params(p);

  const NEWIMAGE::volume<float>&   mean = _sm.GetMean(_field);
  const NEWIMAGE::volume<float>&   mean_alpha = _sm.GetMeanAlpha(_field);
  const NEWIMAGE::volume<float>&   mean_beta = _sm.GetMeanBeta(_field);
  const NEWIMAGE::volume<float>&   mean_gamma = _sm.GetMeanGamma(_field);
  const NEWIMAGE::volume<char>&    mask = _sm.GetMask(_field);
  unsigned int n = static_cast<unsigned int>(mask.sum());
  unsigned int m = _sm.NoOfScans();
  
  // Here we calculate the entities in equation 20 in the tech-report
  //
  // The order of the 6 elements in abc will be aa ab ac bb bc cc
  // Space will be allocated only for those that actually have non-zero values
  std::vector<NEWIMAGE::volume<float> *> abc(6,static_cast<NEWIMAGE::volume<float> *>(0));
  abc[0] = new NEWIMAGE::volume<float>(mean.xsize(),mean.ysize(),mean.zsize());
  copybasicproperties(mean,*(abc[0])); *(abc[0]) = 0.0;
  if (_sm.HasBeta()) {
    abc[1] = new NEWIMAGE::volume<float>(mean.xsize(),mean.ysize(),mean.zsize());
    copybasicproperties(mean,*(abc[1])); *(abc[1]) = 0.0;
    abc[3] = new NEWIMAGE::volume<float>(mean.xsize(),mean.ysize(),mean.zsize());
    copybasicproperties(mean,*(abc[3])); *(abc[3]) = 0.0;
  }   
  if (_sm.HasGamma()) {
    abc[2] = new NEWIMAGE::volume<float>(mean.xsize(),mean.ysize(),mean.zsize());
    copybasicproperties(mean,*(abc[2])); *(abc[2]) = 0.0;
    abc[5] = new NEWIMAGE::volume<float>(mean.xsize(),mean.ysize(),mean.zsize());
    copybasicproperties(mean,*(abc[5])); *(abc[5]) = 0.0;
  }
  if (_sm.HasBeta() && _sm.HasGamma()) {
    abc[4] = new NEWIMAGE::volume<float>(mean.xsize(),mean.ysize(),mean.zsize());
    copybasicproperties(mean,*(abc[4])); *(abc[4]) = 0.0;    
  }
  
  // Now calculate aa ab ac bb bc cc as needed
  for (unsigned int i=0; i<_sm.NoOfScans(); i++) {
    NEWIMAGE::volume<float>  alpha_diff = mean_alpha - _sm.GetAlpha(i,_field);
    *(abc[0]) += (alpha_diff * alpha_diff);
    if (_sm.HasBeta()) {
      NEWIMAGE::volume<float> beta_diff = mean_beta - _sm.GetBeta(i,_field);
      *(abc[1]) += (alpha_diff * beta_diff);
      *(abc[3]) += (beta_diff * beta_diff);
      if (_sm.HasGamma()) {
	NEWIMAGE::volume<float> gamma_diff = mean_gamma - _sm.GetGamma(i,_field);
        *(abc[2]) += (alpha_diff * gamma_diff);
        *(abc[4]) += (beta_diff * gamma_diff);
        *(abc[5]) += (gamma_diff * gamma_diff);
      }
    }
    else if (_sm.HasGamma()) {
      NEWIMAGE::volume<float> gamma_diff = mean_gamma - _sm.GetGamma(i,_field);
      *(abc[2]) += (alpha_diff * gamma_diff);
      *(abc[5]) += (gamma_diff * gamma_diff);
    }
  }
        
  // And use these to calculate the non-linear part of the Hessian
  NEWIMAGE::volume<float> ones(mean.xsize(),mean.ysize(),mean.zsize());
  copybasicproperties(mean,ones); ones = 1.0;
  boost::shared_ptr<MISCMATHS::BFMatrix> nonlin_bit = _field.JtJ(*(abc[0]),ones,&mask,HessianPrecision());
  std::vector<unsigned int> deriv1(3,0);
  std::vector<unsigned int> deriv2(3,0);
  if (_sm.HasBeta()) {
    deriv2[0] = 1;      
    boost::shared_ptr<MISCMATHS::BFMatrix> tmp = _field.JtJ(deriv1,*(abc[1]),deriv2,ones,&mask,HessianPrecision());
    nonlin_bit->AddToMe(*tmp);
    nonlin_bit->AddToMe(*(tmp->Transpose()));
    tmp = _field.JtJ(deriv2,*(abc[3]),ones,&mask,HessianPrecision());
    nonlin_bit->AddToMe(*tmp);
    deriv2[0] = 0;
  }
  if (_sm.HasGamma()) {
    deriv2[1] = 1;
    boost::shared_ptr<MISCMATHS::BFMatrix> tmp = _field.JtJ(deriv1,*(abc[2]),deriv2,ones,&mask,HessianPrecision());
    nonlin_bit->AddToMe(*tmp);
    nonlin_bit->AddToMe(*(tmp->Transpose()));
    tmp = _field.JtJ(deriv2,*(abc[5]),ones,&mask,HessianPrecision());
    nonlin_bit->AddToMe(*tmp);
    deriv2[1] = 0;
  }
  if (_sm.HasBeta() && _sm.HasGamma()) {
    deriv1[0] = 1;
    deriv2[1] = 1;
    boost::shared_ptr<MISCMATHS::BFMatrix> tmp = _field.JtJ(deriv1,*(abc[4]),deriv2,ones,&mask,HessianPrecision());
    nonlin_bit->AddToMe(*tmp);
    nonlin_bit->AddToMe(*(tmp->Transpose()));
  }
  nonlin_bit->MulMeByScalar(2.0/(n*(m-1)));

  // Add regularisation
  if (Lambda()) nonlin_bit->AddToMe(*FieldEnergyHess(),Lambda());

  // Debug printouts
  if (debug_level()>1) {
    NEWIMAGE::volume4D<float> out(mean.xsize(),mean.ysize(),mean.zsize(),6);
    copybasicproperties(mean,out);
    out = 0.0;
    out[0] = *(abc[0]);
    if (_sm.HasBeta()) { out[1] = *(abc[1]); out[3] = *(abc[3]); }
    if (_sm.HasGamma()) { out[2] = *(abc[2]); out[5] = *(abc[5]); }
    if (_sm.HasBeta() && _sm.HasGamma()) out[4] = *(abc[4]);
    write_volume4D(out,string("TopupDebugABC")+debug_string());
  }

  // Clean up a little
  delete abc[0];
  if (_sm.HasBeta()) {
    delete abc[1];
    delete abc[3];
  }   
  if (_sm.HasGamma()) {
    delete abc[2];
    delete abc[5];
  }
  if (_sm.HasBeta() && _sm.HasGamma()) {
    delete abc[4];
  }

  if (!MovementsFixed()) {
    // Now calculate the movement bit (lower right corner).
    std::vector<unsigned int> tile_sizes(_sm.NoOfScans());
    for (unsigned int s=1; s<_sm.NoOfScans(); s++) tile_sizes[s] = _sm.NoOfMovementParametersForScan(s);
    TiledMatrix tiled_move_bit(tile_sizes);
    for (unsigned int s1=1; s1<_sm.NoOfScans(); s1++) {
      for (unsigned int s2=s1; s2<_sm.NoOfScans(); s2++) {
        if (s1 == s2) tiled_move_bit.SetTile(s1,s2,((2.0/(n*(m-1))) * (1.0 - 1.0/m)) * movement_hessian(s1,s2,mask,_field));
        else {
          tiled_move_bit.SetTile(s1,s2,- ((1.0/(n*(m-1))) * (2.0/m)) * movement_hessian(s1,s2,mask,_field));
          tiled_move_bit.SetTile(s2,s1,tiled_move_bit.GetTile(s1,s2).t());
        }
      }
    }
    boost::shared_ptr<MISCMATHS::BFMatrix> move_bit = boost::shared_ptr<MISCMATHS::BFMatrix>(new FullBFMatrix(tiled_move_bit.Untile()));

    // And finally the interaction (the left and bottom "stripes");

    NEWMAT::Matrix interaction(NDefPar(),NMovPar());
    NEWIMAGE::volume<float> beta_diff;
    NEWIMAGE::volume<float> gamma_diff;
    std::vector<unsigned int> deriv(3,0);
    unsigned int offset = 0;
    for (unsigned int s=1; s<_sm.NoOfScans(); s++) {
      NEWIMAGE::volume<float> alpha_diff = mean_alpha - _sm.GetAlpha(s,_field);
      if (_sm.HasBeta()) beta_diff = mean_beta - _sm.GetBeta(s,_field);
      if (_sm.HasGamma()) gamma_diff = mean_gamma - _sm.GetGamma(s,_field);
      for (unsigned d=0; d<_sm.NoOfMovementParametersForScan(s); d++) {
	NEWIMAGE::volume<float> deriv_wrt_move = _sm.GetMovementDerivative(s,d,_field);
        interaction.Column(offset+d+1) = - _field.Jte(deriv_wrt_move,alpha_diff,&mask);
        if (_sm.HasBeta()) { deriv[0]=1; interaction.Column(offset+d+1) -=  _field.Jte(deriv,deriv_wrt_move,beta_diff,&mask); deriv[0]=0; }
        if (_sm.HasGamma()) { deriv[1]=1; interaction.Column(offset+d+1) -=  _field.Jte(deriv,deriv_wrt_move,gamma_diff,&mask); deriv[1]=0; }
      }
      offset += _sm.NoOfMovementParametersForScan(s);
    }
    interaction *= 2.0/(n*(m-1));

    boost::shared_ptr<MISCMATHS::BFMatrix> bf_interaction = boost::shared_ptr<MISCMATHS::BFMatrix>(new FullBFMatrix(interaction));

    // And then put them together

    nonlin_bit->VertConcatBelowMe(*(bf_interaction->Transpose()));
    bf_interaction->VertConcatBelowMe(*move_bit);
    nonlin_bit->HorConcat2MyRight(*bf_interaction);
  }

  /*
  static int counter=0;
  cout << "counter = " << counter << endl;
  if (counter == 24) {
    // 2nd derivatives w.r.t. the warps
    MISCMATHS::write_ascii_matrix(string("numerical_warp_hessian.txt"),numerical_hessian(p,NDefPar()/2,NDefPar()/2+10,1e-1,false));
    MISCMATHS::write_ascii_matrix(string("analytical_warp_hessian.txt"),nonlin_bit->SubMatrix(NDefPar()/2,NDefPar()/2+10,NDefPar()/2,
											      NDefPar()/2+10));
    // Movement 2nd derivatives
    MISCMATHS::write_ascii_matrix(string("numerical_movement_hessian.txt"),numerical_hessian(p,NDefPar()+1,NPar(),1e-3,false));
    MISCMATHS::write_ascii_matrix(string("analytical_movement_hessian.txt"),nonlin_bit->SubMatrix(NDefPar()+1,NPar(),
												  NDefPar()+1,NPar()));
    // Interaction 2nd derivatives
    MISCMATHS::write_ascii_matrix(string("numerical_interaction_hessian.txt"),numerical_hessian(p,NDefPar()/2+1,NDefPar()/2+10,
												NDefPar()+1,NPar(),1e-1,1e-3,false));
    MISCMATHS::write_ascii_matrix(string("analytical_interaction_hessian.txt"),nonlin_bit->SubMatrix(NDefPar()/2+1,NDefPar()/2+10,
												     NDefPar()+1,NPar()));
    exit(EXIT_SUCCESS);
  }
  else counter++;
  */

  if (debug_level() > 2) nonlin_bit->Print(string("TopupDebugHess")+debug_string()+string(".txt"))
;
  if (TracePrint()) cout << "Leaving TopupCF::hess" << endl;
  
  return(nonlin_bit);
}

BASISFIELD::splinefield TopupCF::field_factory(const NEWIMAGE::volume4D<float>& scans,
                                               double                           warpres,
                                               unsigned int                     sporder) const
{
  std::vector<unsigned int>  size(3,0);
  size[0]=scans.xsize(); size[1]=scans.ysize(); size[2]=scans.zsize();
  std::vector<double>        vxs(3,0.0); 
  vxs[0]=scans.xdim(); vxs[1]=scans.ydim(); vxs[2]=scans.zdim();
  std::vector<unsigned int>  ksp(3,0);
  for (unsigned int i=0; i<3; i++) ksp[i] = (static_cast<unsigned int>(MISCMATHS::round(double(warpres/vxs[i])))) ? static_cast<unsigned int>(MISCMATHS::round(double(warpres/vxs[i]))) : 1;
  BASISFIELD::splinefield field(size,vxs,ksp,sporder);

  return(field);
}

void TopupCF::set_field_params(const NEWMAT::ColumnVector& p) const
{
  if (TracePrint()) cout << "Entering TopupCF::set_field_params" << endl;

  if (p.Nrows() != int(NPar())) throw TopupException("Topup_CF::set_field_params: wrong size parameter vector p");
  for (unsigned int i=0; i<_field.CoefSz(); i++) {
    if (fabs(p(i+1)-_field.GetCoef(i)) > 1e-6) {
      _field.SetCoef(p.Rows(1,_field.CoefSz()));
      _sm.FieldUpdated();
      break;
    }
  }

  if (TracePrint()) cout << "Leaving TopupCF::set_field_params" << endl;
}
 
void TopupCF::set_movement_params(const NEWMAT::ColumnVector& p) const
{
  if (TracePrint()) cout << "Entering TopupCF::set_movement_params" << endl;

  if (!_mf) {
    if (p.Nrows() != int(NPar())) throw TopupException("Topup_CF::set_movement_params: wrong size parameter vector p");
    NEWMAT::ColumnVector mp = p.Rows(_field.CoefSz()+1,p.Nrows());
    _sm.SetMovementParameters(mp);
  }

  if (TracePrint()) cout << "Leaving TopupCF::set_movement_params" << endl;
}

double TopupCF::sum_of_prod(const NEWIMAGE::volume<float>& i1,
                            const NEWIMAGE::volume<float>& i2,
                            const NEWIMAGE::volume<char>&  mask) const
{
  if (TracePrint()) cout << "Entering TopupCF::sum_of_prod" << endl;

  double sum = 0.0;

  for (int k=0; k<i1.zsize(); k++) {
    for (int j=0; j<i1.ysize(); j++) {
      for (int i=0; i<i1.xsize(); i++) {
        if (mask(i,j,k)) sum += i1(i,j,k)*i2(i,j,k);
      }
    }
  }

  if (TracePrint()) cout << "Leaving TopupCF::sum_of_prod" << endl;

  return(sum);
}

NEWMAT::ReturnMatrix TopupCF::movement_hessian(unsigned int s1, // Scan 1, row
                                               unsigned int s2, // Scan 2, col
                                               const NEWIMAGE::volume<char>&   mask,
                                               const BASISFIELD::splinefield&  field) const
{
  if (TracePrint()) cout << "Entering TopupCF::movement_hessian" << endl;

  NEWMAT::Matrix rval(_sm.NoOfMovementParametersForScan(s1),_sm.NoOfMovementParametersForScan(s2));
  std::vector<unsigned int> isz = _sm.ImageSize(Subsampled);
  NEWIMAGE::volume4D<float> deriv1(isz[0],isz[1],isz[2],_sm.NoOfMovementParametersForScan(s1));
  for (unsigned int d=0; d<_sm.NoOfMovementParametersForScan(s1); d++) {
    deriv1[d] = _sm.GetMovementDerivative(s1,d,field);
  }
  if (s1 == s2) { // If it is on the diagonal
    for (unsigned int d1=0; d1<_sm.NoOfMovementParametersForScan(s1); d1++) {
      for (unsigned int d2=d1; d2<_sm.NoOfMovementParametersForScan(s1); d2++) {
        rval(d1+1,d2+1) = sum_of_prod(deriv1[d1],deriv1[d2],mask);
        if (d1 != d2) rval(d2+1,d1+1) = rval(d1+1,d2+1);
      }
    }
  }
  else {
    NEWIMAGE::volume4D<float> deriv2(isz[0],isz[1],isz[2],_sm.NoOfMovementParametersForScan(s2));
    for (unsigned int d=0; d<_sm.NoOfMovementParametersForScan(s2); d++) {
      deriv2[d] = _sm.GetMovementDerivative(s2,d,field);
    }
    for (unsigned int d1=0; d1<_sm.NoOfMovementParametersForScan(s1); d1++) {
      for (unsigned int d2=0; d2<_sm.NoOfMovementParametersForScan(s2); d2++) {
        rval(d1+1,d2+1) = sum_of_prod(deriv1[d1],deriv2[d2],mask);
      }
    }
  }

  if (TracePrint()) cout << "Leaving TopupCF::movement_hessian" << endl;

  rval.Release();
  return(rval);
}

NEWMAT::ColumnVector TopupCF::numerical_gradient(const NEWMAT::ColumnVector& p, unsigned int fi, 
						 unsigned int li, double delta, bool reg) const
{
  NEWMAT::ColumnVector  tmp_p = p;
  double old_lambda = _lambda;
  if (!reg) _lambda = 0;

  if (delta < 0) {
    if (fi < NDefPar()) delta = 1e-1;
    else delta = 1e-4;
  }

  double lcf = cf(p);
  NEWMAT::ColumnVector   numd(li-fi+1);

  for (unsigned int gi=fi, i=1; gi<=li; gi++, i++) {
    tmp_p(gi) += delta;
    numd(i) = (cf(tmp_p) - lcf) / delta;
    tmp_p(gi) -= delta;
  }
  if (!reg) _lambda = old_lambda;

  return(numd);
}

NEWMAT::Matrix TopupCF::numerical_hessian(const NEWMAT::ColumnVector& p, unsigned int fi1, unsigned int li1, unsigned int fi2, 
					  unsigned int li2, double delta1, double delta2, bool reg) const
{

  NEWMAT::Matrix hess(li1-fi1+1,li2-fi2+1);
  NEWMAT::ColumnVector  tmp_p = p;
  double old_lambda = _lambda;
  if (reg) _lambda = 0;
  if (delta1 < 0) {
    if (fi1 < NDefPar()) delta1 = 1e-1;
    else delta1 = 1e-4;
  }
  if (delta2 < 0) {
    if (fi2 < NDefPar()) delta2 = 1e-1;
    else delta2 = 1e-4;
  }

  double lcf = cf(p);
  for (unsigned int row=fi1, i=1; row<=li1; row++, i++) {
    tmp_p(row) += delta1;
    double first = (cf(tmp_p) - lcf) / delta1;
    tmp_p(row) -= delta1;
    for (unsigned int col=fi2, j=1; col<=li2; col++, j++) {
      tmp_p(col) += delta2;
      double tmp_cf = cf(tmp_p);
      tmp_p(row) += delta1;
      double second = (cf(tmp_p) - tmp_cf) / delta1;
      tmp_p(col) -= delta2;
      tmp_p(row) -= delta1;
      hess(i,j) = (second-first) / delta2;
    }  
  }
  if (!reg) _lambda = old_lambda;

  return(hess);
}

// }}} End of fold
