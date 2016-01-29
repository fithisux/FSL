// Definitions for class basisfield
//
// basisfield.cpp
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
//
#include <string>
#include <iostream>
#include "newmat.h"
#include "newimage/newimage.h"
#include "miscmaths/bfmatrix.h"
#include "basisfield.h"

using namespace std;
using namespace NEWMAT;
using namespace NEWIMAGE;

namespace BASISFIELD {

// Constructor, assignement and destructor

basisfield::basisfield(const std::vector<unsigned int>& psz, const std::vector<double>& pvxs)
: ndim(psz.size()), sz(3,1), vxs(3,0.0), coef(), futd(4,false), field(4)
{
  if (psz.size()<1 || psz.size()>3) {throw BasisfieldException("basisfield::basisfield::Invalid dimensionality of field");}
  if (psz.size() != pvxs.size()) {throw BasisfieldException("basisfield::basisfield:: Dimensionality mismatch between psz and pvxs");}
  for (int i=0; i<3; i++) {
    if (i<int(psz.size()) && (psz[i] < 1 || psz[i] > MAX_SIZE)) {throw BasisfieldException("basisfield::basisfield::Invalid size of field");}
    sz[i] = (i<int(psz.size())) ? psz[i] : 1;
    vxs[i] = (i<int(pvxs.size())) ? pvxs[i] : 0.0;
  }
}

basisfield::basisfield(const basisfield& inf)
: ndim(inf.ndim), sz(3,1), vxs(3,0.0), coef(), futd(4,false), field(4)
{
  assign_basisfield(inf);
}

basisfield& basisfield::operator=(const basisfield& inf)
{
  if (&inf == this) {return(*this);} // Detect self
  assign_basisfield(inf);
  return(*this);
}


basisfield::~basisfield() {}

// General utility functions

ReturnMatrix basisfield::mm2vox(unsigned int sz) const
{
  Matrix rmat = IdentityMatrix(sz);
  rmat(1,1) = 1.0/Vxs_x();
  rmat(2,2) = 1.0/Vxs_y();
  rmat(3,3) = 1.0/Vxs_z();
  rmat.Release();
  return(rmat);
}

ReturnMatrix basisfield::vox2mm(unsigned int sz) const
{
  Matrix rmat = IdentityMatrix(sz);
  rmat(1,1) = Vxs_x();
  rmat(2,2) = Vxs_y();
  rmat(3,3) = Vxs_z();
  rmat.Release();
  return(rmat);
}

double basisfield::Peek(unsigned int x, unsigned int y, unsigned int z, FieldIndex fi)
{
  if ( x>=FieldSz_x() || y>=FieldSz_y() || z>=FieldSz_z()) { 
    throw BasisfieldException("basisfield::PeekField:: Co-ordinates out of bounds");
  }
  if (!coef) {return(0.0);} // Consider field as zero if no coefficients set
  if (!UpToDate(fi)) {Update(fi);}

  return(UnsafePeek(z*FieldSz_x()*FieldSz_y()+y*FieldSz_x()+x,fi));
}    

double basisfield::Peek(unsigned int vi, FieldIndex fi)
{
  if (vi>=FieldSz()) {throw BasisfieldException("basisfield::PeekField:: Voxel index out of bounds");}
  if (!coef) {return(0.0);} // Consider field as zero if no coefficients set
  if (!UpToDate(fi)) {Update(fi);}

  return(UnsafePeek(vi,fi));
}

double basisfield::PeekWide(int i, int j, int k, FieldIndex fi)
{
  if (!(i<0 || j<0 || k<0 || static_cast<unsigned int>(i)>=FieldSz_x() || static_cast<unsigned int>(j)>=FieldSz_y() || static_cast<unsigned int>(k)>=FieldSz_z())) {  // Inside "valid" FOV
    return(Peek(static_cast<unsigned int>(i),static_cast<unsigned int>(j),static_cast<unsigned int>(k),fi));
  }
  else {
    return(peek_outside_fov(i,j,k,fi));
  }
}
  
void basisfield::SetCoef(const ColumnVector& pcoef) 
{
  if (pcoef.Nrows() != int(CoefSz())) {throw BasisfieldException("basisfield::SetCoef::Mismatch between input vector and # of coefficients");}
  if (!coef) {coef = boost::shared_ptr<NEWMAT::ColumnVector>(new NEWMAT::ColumnVector(pcoef));}
  else {*coef = pcoef;}
  futd.assign(4,false);
}

/////////////////////////////////////////////////////////////////////
//
// Calulates and sets coefficients such that the field is the best
// possible approximation to the supplied field. The current
// implementation is not as efficient as it could be.
//
/////////////////////////////////////////////////////////////////////

void basisfield::Set(const volume<float>& pfield)
{
  if (int(FieldSz_x()) != pfield.xsize() || int(FieldSz_y()) != pfield.ysize() || int(FieldSz_z()) != pfield.zsize()) {
    throw BasisfieldException("basisfield::Set:: Matrix size mismatch beween basisfield class and supplied field");
  }
  if (Vxs_x() != pfield.xdim() || Vxs_y() != pfield.ydim() || Vxs_z() != pfield.zdim()) {
    throw BasisfieldException("basisfield::Set:: Voxel size mismatch beween basisfield class and supplied field");
  }

  volume<float>   volume_of_ones(pfield.xsize(),pfield.ysize(),pfield.zsize());
  volume_of_ones.copyproperties(pfield);
  volume_of_ones = 1.0;

  double lambda = 0.001;
  ColumnVector y = Jte(pfield,0);
  boost::shared_ptr<MISCMATHS::BFMatrix>  XtX = JtJ(volume_of_ones);
  boost::shared_ptr<MISCMATHS::BFMatrix>  BeEn = BendEnergyHess();
  XtX->AddToMe(*BeEn,lambda);
  ColumnVector coef_roof = XtX->SolveForx(y,SYM_POSDEF,1e-6,500);
  SetCoef(coef_roof);
}

/////////////////////////////////////////////////////////////////////
//
// Optional gateway to the Set( volume<float>& ) method. 
//
/////////////////////////////////////////////////////////////////////

void basisfield::Set(const ColumnVector& pfield) 
{
  if (pfield.Nrows() != int(FieldSz())) {throw BasisfieldException("basisfield::Set::Mismatch between input vector and size of field");}

  volume<float>  vol_pfield(FieldSz_x(),FieldSz_y(),FieldSz_z());
  vol_pfield.setdims(Vxs_x(),Vxs_y(),Vxs_z());
  vol_pfield.insert_vec(pfield);
  
  Set(vol_pfield);
}


void basisfield::AsVolume(volume<float>& vol, FieldIndex fi)
{
  if (int(FieldSz_x()) != vol.xsize() || int(FieldSz_y()) != vol.ysize() || int(FieldSz_z()) != vol.zsize()) {
    throw BasisfieldException("basisfield::AsVolume:: Matrix size mismatch beween field and volume");
  }
  if (Vxs_x() != vol.xdim() || Vxs_y() != vol.ydim() || Vxs_z() != vol.zdim()) {
    throw BasisfieldException("basisfield::AsVolume:: Voxel size mismatch beween field and volume");
  }
  if (!coef) {throw BasisfieldException("basisfield::AsVolume: Coefficients undefined");}

  if (!UpToDate(fi)) {Update(fi);}

  const boost::shared_ptr<NEWMAT::ColumnVector> tmptr = Get(fi);
  int vindx=0;
  for (unsigned int k=0; k<FieldSz_z(); k++) {
    for (unsigned int j=0; j<FieldSz_y(); j++) {
      for (unsigned int i=0; i<FieldSz_x(); i++) {
        vol(i,j,k) = tmptr->element(vindx++);
      }
    }
  }
}

// Functions that are declared private or protected

void basisfield::assign_basisfield(const basisfield& inf) // Helper function for copy constructor and assignment
{
  futd = inf.futd;
  ndim = inf.ndim;
  vxs = inf.vxs;
  sz = inf.sz;
  coef = boost::shared_ptr<NEWMAT::ColumnVector>(new NEWMAT::ColumnVector(*(inf.coef)));
  for (int i=0; i<int(inf.field.size()); i++) {
    if (inf.field[i]) {field[i] = boost::shared_ptr<NEWMAT::ColumnVector>(new NEWMAT::ColumnVector(*(inf.field[i])));}
    else {field[i] = inf.field[i];}
  }
}

boost::shared_ptr<NEWMAT::ColumnVector> basisfield::get(FieldIndex fi)
{
  if (!coef) {throw BasisfieldException("basisfield::Get: Coefficients undefined");}

  if (!UpToDate(fi)) {
    Update(fi);
  }
  return(field[fi]);
}

boost::shared_ptr<NEWMAT::ColumnVector> basisfield::get_ptr(FieldIndex fi)
{
  if (!field[fi]) {
    field[fi] = boost::shared_ptr<NEWMAT::ColumnVector>(new NEWMAT::ColumnVector(FieldSz()));
  }
  return(field[fi]);
}

boost::shared_ptr<NEWMAT::ColumnVector> basisfield::get_coef() const
{
  if (!coef) {throw BasisfieldException("basisfield::get_coef: Coefficients undefined");}
  return(coef);  
}

double basisfield::get_coef(unsigned int i, unsigned int j, unsigned int k) const
{
  if (i >= CoefSz_x() || j >= CoefSz_y() || k >= CoefSz_z()) throw BasisfieldException("basisfield::get_coef: i, j or k out of range");
  return(coef->element(k*CoefSz_x()*CoefSz_y()+j*CoefSz_x()+i));
}

double basisfield::get_coef(unsigned int i) const
{
  if (i >= CoefSz()) throw BasisfieldException("basisfield::get_coef: i out of range");
  return(coef->element(i));
}


} // End namespace BASISFIELD
