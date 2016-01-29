// Definitions of PointList classes used by FNIRT
//
// point_list.cpp
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2008-2012 University of Oxford 
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
#include <string>
#include <fstream>
#include <iostream>
#include "newmat.h"
#include "newmatio.h"

#ifndef EXPOSE_TREACHEROUS
#define I_DEFINED_ET
#define EXPOSE_TREACHEROUS
#endif

#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "meshclass/meshclass.h"
#include "point_list.h"

using namespace std;
using namespace NEWMAT;
using namespace NEWIMAGE;
using namespace mesh;

PointList::PointList(const std::string& pfname, const std::string& ifname) : _pfname(pfname), _ifname(ifname) 
{
  _affine = IdentityMatrix(4);

  // First let us read image file
  volume<float>  vol;
  try { read_volume_hdr_only(vol,ifname); }
  catch (...) { throw PointListException(string("PointList: Can not read file: ")+ifname); }
  _vox2mm = vol.sampling_mat();
  if (vol.sform_code() != NIFTI_XFORM_UNKNOWN) _SQform = vol.sform_mat();
  else if (vol.qform_code() != NIFTI_XFORM_UNKNOWN) _SQform = vol.qform_mat();
  else _SQform = IdentityMatrix(4);

  if (is_maybe_vtk(pfname)) {
    if (PointList::read_first_file(pfname)<0) {
      throw PointListException(string("PointList: Can not read vtk-file: ")+pfname);
    }
  }
  else { // Assume text file
    if (PointList::read_ascii_file(pfname)<0) {
      throw PointListException(string("PointList: Can not read text-file: ")+pfname);
    }
  }  
}

NEWMAT::ReturnMatrix PointList::Point(unsigned int indx, bool one_ext) const
{
  if (indx >= NPoints()) throw PointListException("Point: Index out of range");
  ColumnVector vec(4);
  vec.Rows(1,3) = _points.Column(indx+1); vec(4) = 1.0;
  vec = _affine * vec;
  if (one_ext) {
    vec.Release();
    return(vec);
  }
  ColumnVector rvec = vec.Rows(1,3);
  rvec.Release();
  return(rvec);
}

NEWMAT::ReturnMatrix PointList::RawPoint(unsigned int indx, bool one_ext) const
{
  if (indx >= NPoints()) throw PointListException("Point: Index out of range");
  if (one_ext) {
    ColumnVector vec(4);
    vec.Rows(1,3) = _points.Column(indx+1); vec(4) = 1.0;
    vec.Release();
    return(vec);
  }
  ColumnVector vec(3);
  vec = _points.Column(indx+1);
  vec.Release();
  return(vec);
}

NEWMAT::ReturnMatrix PointList::Vox2mm() const
{
  Matrix rmat = _vox2mm;
  rmat.Release();
  return(rmat);
}

NEWMAT::ReturnMatrix PointList::mm2Vox() const
{
  Matrix rmat = _vox2mm.i();
  rmat.Release();
  return(rmat);
}

void PointList::SetAffine(const NEWMAT::Matrix& aff)
{
  if (aff.Nrows() != 4 || aff.Ncols() != 4) throw PointListException("SetAffine: Matrix must be 4x4");
  _affine = aff;
}

// This routine reads a list of points from a text-file that 
// the user has (presumably) entered manually. It will be assumed
// that the user has used fslview to deduce these points, and that
// he/she has entered the coordinates from the second column of the
// coordinate box. In this box will be "FSL mm-coordinates" if no
// sform/qform has been set or "world coordinates" if there is
// a valid sform/qform.
// The PointList is expected to return "FSL mm-coordinates" so
// if it is in "world coordinates" we will need to undo that.
int PointList::read_ascii_file(const std::string& fname)
{
  Matrix pl = MISCMATHS::read_ascii_matrix(fname);
  if (pl.Nrows() == 0 || pl.Ncols() == 0 || (pl.Nrows() != 3 && pl.Ncols() != 3)) return(-1);
  if (pl.Ncols() == 3) pl = pl.t();

  if ((_SQform-IdentityMatrix(4)).MaximumAbsoluteValue() > 1e-6) { // If q/s form is set.
    _points.ReSize(pl.Nrows(),pl.Ncols());
    Matrix M2mm = _vox2mm * _SQform.i();
    for (int i=1; i<=pl.Ncols(); i++) {
      ColumnVector pt(4);
      pt.Rows(1,3) = pl.Column(i); pt(4) = 1.0;
      pt = M2mm * pt;
      _points.Column(i) = pt.Rows(1,3);
    }
  }
  else _points = pl;  
  return(1);
}

// This routines reads a list of points generated by first
// or other "first-related" routines that uses the same
// output format. They are expected to already be in
// "FSL mm-coordinates".
int PointList::read_first_file(const std::string& fname)
{
  // Wrapping this in a try since the meshclass reader sometimes
  // returns -1 and sometimes crashes when encountering invalid
  // .vtk files. And sometimes it crashes without there being a
  // proper exception, so it is not possible to catch. Phew!
  try {
    Mesh pl;
    if (pl.load(fname)<0) return(-1);
    _points.ReSize(3,pl.nvertices());
    for (vector<Mpoint *>::iterator i=pl._points.begin(); i!=pl._points.end(); i++) {
      int indx = (*i)->get_no();
      Pt p = pl.get_point(indx)->get_coord();
      _points(1,indx+1) = p.X; _points(2,indx+1) = p.Y; _points(3,indx+1) = p.Z; 
    }
    return(1);
  }
  catch (const std::exception& error) {
    cerr << "PointList::read_first_file: Exception thrown with message: " << error.what() << endl; 
    return(-1);
  }
  return(-1);
}

bool PointList::is_maybe_vtk(const std::string& fname) const
{
  ifstream fs;
  fs.open(fname.c_str());
  if (fs.is_open()) {
    string tag;
    getline(fs,tag);
    fs.clear(); fs.close();
    if (tag.find("# vtk DataFile") == string::npos) return(false);
    else return(true);
  }
  return(false);
}

#ifdef I_DEFINED_ET
#undef I_DEFINED_ET
#undef EXPOSE_TREACHEROUS   // Avoid exporting dodgy routines
#endif
