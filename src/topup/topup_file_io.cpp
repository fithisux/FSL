// Definitions of classes used to write and
// read files written by topup, and potentially
// by other pieces of software as long as they
// are valid magnetic-field files.
//
// topup_file_io.cpp
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2009 University of Oxford 
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


#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
#include "newmat.h"

#ifndef EXPOSE_TREACHEROUS
#define EXPOSE_TREACHEROUS           // To allow us to use .sampling_mat()
#endif

#include "newimage/newimageall.h"
#include "basisfield/basisfield.h"
#include "basisfield/splinefield.h"
#include "basisfield/dctfield.h"
#include "topup_file_io.h"

namespace TOPUP {

// Global functions

NEWMAT::Matrix MovePar2Matrix(const NEWMAT::ColumnVector&     mp,
                              const NEWIMAGE::volume<float>&  vol) 
{
  if (mp.Nrows() != 6) throw TopupFileIOException("MovePar2Matrix: mp must have 6 elements");

  NEWMAT::ColumnVector tmp(6);
  tmp(1) = mp(4); tmp(2) = mp(5); tmp(3) = mp(6);
  tmp(4) = mp(1); tmp(5) = mp(2); tmp(6) = mp(3);
  
  NEWMAT::ColumnVector cntr(3);
  cntr(1) = ((vol.xsize()-1)*vol.xdim())/2.0;
  cntr(2) = ((vol.ysize()-1)*vol.ydim())/2.0;
  cntr(3) = ((vol.zsize()-1)*vol.zdim())/2.0;

  NEWMAT::Matrix mat(4,4);
  MISCMATHS::construct_rotmat_euler(tmp,6,mat,cntr);

  return(mat);
}

NEWMAT::ColumnVector Matrix2MovePar(const NEWMAT::Matrix&           M,
				    const NEWIMAGE::volume<float>&  vol)
{
  NEWMAT::ColumnVector mp(6); mp = 0.0;
  NEWMAT::ColumnVector rot(3);

  MISCMATHS::rotmat2euler(rot,M);
  mp.Rows(4,6) = rot;

  NEWMAT::Matrix MM = MovePar2Matrix(mp,vol);
  mp.Rows(1,3) = M.SubMatrix(1,3,4,4) - MM.SubMatrix(1,3,4,4);

  return(mp);  
}

/////////////////////////////////////////////////////////////////////
//
//  Definitions for class TopupFileWriter
//
/////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
//
// Constructor for coefficient file
//
// The coefficients are saved in a slightly dodgy format where the fields of the niftii
// header are used to store information that is neccessary for us to reconstruct the 
// displacement fields from the coefficients. E.g. the knot-spacings are stored in the 
// "pixdims" and the matrix size of the field/template in the offsets of the qform.
// For this reason the coefficient files need to be read/written usinge the _orig_ forms
// of the read/write_volume functions. Otherwise the i/o-functions would interpret these
// fields and potentially left-right swap on read/write.
//
////////////////////////////////////////////////////////////////////////////////////////////////
TopupFileWriter::TopupFileWriter(const std::string&                 fname, 
                                 const BASISFIELD::splinefield&     field)
{
  NEWIMAGE::volume<float>         coefs(int(field.CoefSz_x()),int(field.CoefSz_y()),int(field.CoefSz_z()));
  std::vector<float>              ksp(3,1.0);

  ksp[0] = float(field.Ksp_x()); ksp[1] = float(field.Ksp_y()); ksp[2] = float(field.Ksp_z());
  if (field.Order() == 2) coefs.set_intent(FSL_TOPUP_QUADRATIC_SPLINE_COEFFICIENTS,field.Vxs_x(),field.Vxs_y(),field.Vxs_z());
  else if (field.Order() == 3) coefs.set_intent(FSL_TOPUP_CUBIC_SPLINE_COEFFICIENTS,field.Vxs_x(),field.Vxs_y(),field.Vxs_z());
  coefs.setxdim(ksp[0]); coefs.setydim(ksp[1]); coefs.setzdim(ksp[2]);
  NEWMAT::Matrix  qform(4,4); 
  qform = IdentityMatrix(4);
  qform(1,4) = field.FieldSz_x(); 
  qform(2,4) = field.FieldSz_y(); 
  qform(3,4) = field.FieldSz_z();
  coefs.set_qform(NIFTI_XFORM_SCANNER_ANAT,qform);
  for (unsigned int k=0; k<field.CoefSz_z(); k++) {
    for (unsigned int j=0; j<field.CoefSz_y(); j++) {
      for (unsigned int i=0; i<field.CoefSz_x(); i++) {
        coefs(i,j,k) = field.GetCoef(i,j,k);
      }
    }
  }

  save_orig_volume(coefs,fname);
}
// Constructor for field file
TopupFileWriter::TopupFileWriter(const std::string&                   fname,
                                 const NEWIMAGE::volume<float>&       ref,
                                 const BASISFIELD::splinefield&       field)
{
  // Make sure matrix- and voxel-sizes are the same in reference and field files
  if (ref.xsize() != int(field.FieldSz_x()) || ref.ysize() != int(field.FieldSz_y()) || ref.zsize() != int(field.FieldSz_z())) {
    throw TopupFileIOException("TopupFileWriter::TopupFileWriter: Ref scan and field has different matrix size");
  }
  if (ref.xdim() != field.Vxs_x() || ref.ydim() != field.Vxs_y() || ref.zdim() != field.Vxs_z()) {
    throw TopupFileIOException("TopupFileWriter::TopupFileWriter: Ref scan and field has different voxel size");
  }
  NEWIMAGE::volume<float>   out = ref;
  BASISFIELD::splinefield&  non_const_field = const_cast<BASISFIELD::splinefield &>(field);  // Must rewrite field classes
  non_const_field.AsVolume(out);
  out.set_intent(FSL_TOPUP_FIELD,out.intent_param(0),out.intent_param(1),out.intent_param(2));
  out.setDisplayMaximum(0.0);
  out.setDisplayMinimum(0.0);
  save_volume(out,fname);
}
// Constructor for movement parameter file
TopupFileWriter::TopupFileWriter(const std::string&       fname,
                                 const NEWMAT::Matrix&    mp)
{
  if (write_ascii_matrix(fname,mp) < 0) throw TopupFileIOException(string("TopupFileWriter::TopupFileWriter: Failed to write movement parameter file ")+fname);
}
// Constructor for rigid body matrix files
TopupFileWriter::TopupFileWriter(const std::string&                    fname,
                                 const std::vector<NEWMAT::Matrix>&    M)
{
  NEWMAT::Matrix  omat(M.size()*M[0].Nrows(),M[0].Ncols());
  for (unsigned int i=0; i<M.size(); i++) {
    if (M[i].Nrows() != M[0].Nrows() || M[i].Ncols() != M[0].Ncols())
      omat.Rows((i-1)*M[0].Nrows()+1,i*M[0].Nrows()) = M[i];
  }
  if (write_ascii_matrix(fname,omat) < 0) throw TopupFileIOException(string("TopupFileWriter::TopupFileWriter: Failed to write matrix file ")+fname);
}

/////////////////////////////////////////////////////////////////////
//
//  Definitions for class TopupFileReader
//
/////////////////////////////////////////////////////////////////////

TopupFileReader::TopupFileReader(const std::string& fname)
{
  // Make sure no extension was given
  string::size_type dotidx = fname.find_last_of(".");
  if (dotidx != string::npos) { // If there is a dot
    string::size_type eopidx = fname.find_last_of("/");
    if (eopidx == string::npos || (eopidx != string::npos && dotidx > eopidx)) {
      throw TopupFileIOException(string("TopupFileReader::TopupFileReader: Filename must be given without extension")+fname);
    }
  } 
  common_read(fname+string("_fieldcoef"));
  ReadMovements(fname+string("_movpar.txt"));
}

void TopupFileReader::ReadMovements(const std::string& fname)
{
  if (_vol_rep) {
    read_movement(fname,*_vol_rep);
  }
  else if (_coef_rep) {
    NEWIMAGE::volume<float>  tmp(_coef_rep->FieldSz_x(),_coef_rep->FieldSz_y(),_coef_rep->FieldSz_z());
    tmp.setdims(_coef_rep->Vxs_x(),_coef_rep->Vxs_y(),_coef_rep->Vxs_z());
    read_movement(fname,tmp);
  }
  else throw TopupFileIOException("TopupFileReader::ReadMovements: Attempting to read movement parameters before knowing FOV");
}

void TopupFileReader::common_read(const std::string& fname)
{
  // Read volume given by fname
  NEWIMAGE::volume<float>  vol;
  read_volume_hdr_only(vol,fname);

  // Decode it according to intent code of volume
  switch (vol.intent_code()) {
  case FSL_TOPUP_CUBIC_SPLINE_COEFFICIENTS:
  case FSL_TOPUP_QUADRATIC_SPLINE_COEFFICIENTS:
    read_orig_volume(vol,fname);
    _coef_rep = read_coef_file(vol);
    if (_vol_rep) _vol_rep = boost::shared_ptr<NEWIMAGE::volume<float> >();
    break;
  case FSL_TOPUP_FIELD:
    read_volume(vol,fname);
    _vol_rep = boost::shared_ptr<NEWIMAGE::volume<float> >(new NEWIMAGE::volume<float>(vol));
    if (_coef_rep) _coef_rep = boost::shared_ptr<BASISFIELD::splinefield>();
    break;
  default: // Assume that it is a field (e.g. from fugue)
    read_volume(vol,fname);
    _vol_rep = boost::shared_ptr<NEWIMAGE::volume<float> >(new NEWIMAGE::volume<float>(vol));
    if (_coef_rep) _coef_rep = boost::shared_ptr<BASISFIELD::splinefield>();
    break;
  }  
}

void TopupFileReader::ensure_volume() const
{
  if (!_vol_rep) {
    if (!_coef_rep) TopupFileIOException("TopupFileReader::ensure_volume: Attempt to access uninitialised field");
    _vol_rep = boost::shared_ptr<NEWIMAGE::volume<float> >(new NEWIMAGE::volume<float>(_coef_rep->FieldSz_x(),_coef_rep->FieldSz_y(),_coef_rep->FieldSz_z()));
    _vol_rep->setdims(_coef_rep->Vxs_x(),_coef_rep->Vxs_y(),_coef_rep->Vxs_z());
    _coef_rep->AsVolume(*_vol_rep);
  }
}

void TopupFileReader::ensure_field() const
{
  if (!_coef_rep) {
    TopupFileIOException("TopupFileReader::ensure_field: Attempt to access uninitialised field");
  }
}
/////////////////////////////////////////////////////////////////////
//
// Read coefficient file generated by Topup
//
/////////////////////////////////////////////////////////////////////

boost::shared_ptr<BASISFIELD::splinefield> TopupFileReader::read_coef_file(const NEWIMAGE::volume<float>& vcoef)
{
  // Collect info needed to re-create the field
  NEWMAT::Matrix  qform = vcoef.qform_mat();
  std::vector<unsigned int>  sz(3,0);
  std::vector<double>        vxs(3,0.0);
  for (int i=0; i<3; i++) {
    sz[i] = static_cast<unsigned int>(qform(i+1,4));
    vxs[i] = static_cast<double>(vcoef.intent_param(i+1));
  }
  std::vector<unsigned int>  ksp(3,0);
  unsigned int          order = 3;
  if (vcoef.intent_code() == FSL_TOPUP_QUADRATIC_SPLINE_COEFFICIENTS) order = 2;
  ksp[0] = static_cast<unsigned int>(vcoef.xdim() + 0.5);
  ksp[1] = static_cast<unsigned int>(vcoef.ydim() + 0.5);
  ksp[2] = static_cast<unsigned int>(vcoef.zdim() + 0.5);
  // Create the field
  boost::shared_ptr<BASISFIELD::splinefield>  field(new BASISFIELD::splinefield(sz,vxs,ksp,order));
  // Check for internal consistency
  if (field->CoefSz_x() != static_cast<unsigned int>(vcoef.xsize()) ||
      field->CoefSz_y() != static_cast<unsigned int>(vcoef.ysize()) ||
      field->CoefSz_z() != static_cast<unsigned int>(vcoef.zsize())) {
    throw TopupFileIOException("read_coef_file: Internally inconsistent coef-file");
  }
  // Set the coefficients from the file
  field->SetCoef(vcoef.vec());

  return(field);
}

/////////////////////////////////////////////////////////////////////
//
// Read movement parameters generated by Topup. This can be either
// an nx6 matrix with 6 movement parameters for n scans, or an 
// n*4x4 matrix with a transformation matrix for each of n scans.
//
/////////////////////////////////////////////////////////////////////

void TopupFileReader::read_movement(const std::string&              fname,
                                    const NEWIMAGE::volume<float>&  vol)
{
  NEWMAT::Matrix tmp = read_ascii_matrix(fname);
  if (tmp.Ncols() == 6) { // Assume 6 movement parameters per row
    _move.resize(tmp.Nrows());
    _mp.resize(tmp.Nrows());
    for (int i=0; i<tmp.Nrows(); i++) {
      _mp[i] = tmp.Row(i+1).t();
      _move[i] = MovePar2Matrix(_mp[i],vol);
      _mp_valid = true;
    }
  }
  else if (tmp.Ncols() == 4 && !(tmp.Nrows()%4)) {
    _move.resize(tmp.Nrows()/4);
    for (int i=0; i<(tmp.Nrows()/4); i++) {
      _move[i] = tmp.SubMatrix(i*4+1,(i+1)*4,1,4);
      _mp_valid = false;
    }
  }
  return;
}

} // End namespace TOPUP
