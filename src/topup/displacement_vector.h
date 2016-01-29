//  
//  Declarations for displacement vector class DispVec
//
//  displacement_vector.h
//
//  Implements a displacement vector class that can be
//  used to obtain inverses, K-matrices etc.
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2010 University of Oxford 
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
#ifndef displacement_vector_h
#define displacement_vector_h

#include <string>
#include <vector>
#include <cmath>
#include "newmat.h"
#include "newimage/newimage.h"
#include "miscmaths/SpMat.h"

namespace TOPUP {

class DispVecException: public std::exception
{
private:
  std::string m_msg;
public:
  DispVecException(const std::string& msg) throw(): m_msg(msg) {}

  virtual const char * what() const throw() {
    return string("DispVec::" + m_msg).c_str();
  }

  ~DispVecException() throw() {}
};

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class CircularArray:
// Helper class that facilitates implementation of wrap-around.
// Interface includes:
// operator[i]: Allows you to access it as a circular array, so e.g.
//              ca[-1] will return the v[n-1] and ca[n] would return
//              v[0].
// IndexInRange(i): Translates the index i into the index of the 
//                  actual element that it would access. So e.g.
//                  ca.IndexInRange(-1) would return n-1 and
//                  ca.IndexInRange(n) returns 0.
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class CircularArray
{
public:
  CircularArray() : _n(0), _v(0), _sf(1.0) {}
  CircularArray(unsigned int n) : _n(n), _v(new double[_n]), _sf(1.0) { for (unsigned int i=0; i<_n; i++) _v[i]=0.0; }
  CircularArray(const std::vector<float>& v) : _n(v.size()), _v(new double[_n]), _sf(1.0) { for (unsigned int i=0; i<_n; i++) _v[i]=static_cast<double>(v[i]); }
  CircularArray(const std::vector<double>& v) : _n(v.size()), _v(new double[_n]), _sf(1.0) { for (unsigned int i=0; i<_n; i++) _v[i]=v[i]; }
  CircularArray(const NEWMAT::ColumnVector& v) : _n(v.Nrows()), _v(new double[_n]), _sf(1.0) { for (unsigned int i=0; i<_n; i++) _v[i]=static_cast<double>(v(i+1)); }
  ~CircularArray() { delete [] _v; }
  void Set(const std::vector<float>& v) { if (_n!=v.size()) { if (_v) delete [] _v; _n=v.size(); _v=new double[_n]; } for (unsigned int i=0; i<_n; i++) _v[i]=double(v[i]); }
  void Set(const std::vector<double>& v) { if (_n!=v.size()) { if (_v) delete [] _v; _n=v.size(); _v=new double[_n]; } for (unsigned int i=0; i<_n; i++) _v[i]=v[i]; }
  void Set(const NEWMAT::ColumnVector& v) { if (int(_n)!=v.Nrows()) { if (_v) delete [] _v; _n=v.Nrows(); _v=new double[_n]; } for (unsigned int i=0; i<_n; i++) _v[i]=double(v(i+1)); }
  void SetFromRow(const NEWIMAGE::volume<float>& ima, unsigned int k, unsigned int j) { set_from_row_or_col(ima,k,j,true); }
  void SetFromColumn(const NEWIMAGE::volume<float>& ima, unsigned int k, unsigned int i) { set_from_row_or_col(ima,k,i,false); }
  void Print(const std::string& fname) const;
  unsigned int N() const { return(_n); }
  void SetScaleFactor(double sf) const { _sf=sf; }
  double GetScaleFactor() const { return(_sf); }
  int Find(double x) const;
  double operator[](int i) const { int j=i%_n; if (j<0) return(i+_sf*_v[_n+j]); else return(i+_sf*_v[j]); }
  double Inv(double x) const;
  unsigned int IndexInRange(int i) const { int j=i%static_cast<int>(_n); if (j<0) return(static_cast<unsigned int>(_n+j)); else return(static_cast<unsigned int>(j)); }
private:
  unsigned int    _n;
  double          *_v;
  mutable double  _sf;

  void set_from_row_or_col(const NEWIMAGE::volume<float>& ima, int k, int ij, bool row);
};

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class DispVec:
// Holds one row/column of a displacement field. The interface
// includes:
// SetFromRow: Picks values from a row in a volume.
// SetFromColumn: Picks values from a column in a volume.
// GetK_Matrix: Returns a corresponding K-matrix such that
//              y=K*x gives a distorted vector y from a "true"
//              vector x. For this to work the values in it
//              should be scaled to units of voxels. If this is
//              not the case (let's say it is in units of Hz)
//              one can use the form K=v.GetK_Matrix(scale_fac)
//              where in this case scale_fac would be readout
//              time in seconds. It can also be negative to
//              facilitate its use for top-down-bottom-up 
//              calculations.
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class DispVec
{
public:
  DispVec() {}
  DispVec(unsigned int n) : _ca(n) {}
  DispVec(const std::vector<float>& v) : _ca(v) {}
  DispVec(const std::vector<double>& v) : _ca(v) {}
  DispVec(const NEWMAT::ColumnVector& v) : _ca(v) {}
  ~DispVec() {}
  void Set(const std::vector<float>& v) { _ca.Set(v); }
  void Set(const std::vector<double>& v) { _ca.Set(v); }
  void Set(const NEWMAT::ColumnVector& v) { _ca.Set(v); }
  void SetScaleFactor(double sf) { _ca.SetScaleFactor(sf); }
  double GetScaleFactor() const { return(_ca.GetScaleFactor()); }
  void SetFromRow(const NEWIMAGE::volume<float>& ima, int k, int j) { _ca.SetFromRow(ima,k,j); } 
  void SetFromColumn(const NEWIMAGE::volume<float>& ima, int k, int i) { _ca.SetFromColumn(ima,k,i); }
  bool RowIsAlright(const NEWIMAGE::volume<float>& mask, int slice, int row) const;
  bool ColumnIsAlright(const NEWIMAGE::volume<float>& mask, int slice, int col) const;
  void Print(const std::string fname=string("")) const { _ca.Print(fname); }
  double operator[](int i) const { return(_ca[i]); }
  double Inv(int i) const { return(_ca.Inv(static_cast<double>(i))); }
  NEWMAT::ReturnMatrix GetDisplacements() const;
  NEWMAT::ReturnMatrix GetInverseDisplacements(double sf=1.0) const;
  NEWMAT::ReturnMatrix GetK_Matrix(double sf=1.0) const;
  NEWMAT::ReturnMatrix GetS_Matrix(bool wrap=true) const;
  MISCMATHS::SpMat<double> GetSparseK_Matrix(double sf=1.0) const { NEWMAT::Matrix K=GetK_Matrix(sf); return(MISCMATHS::SpMat<double>(K)); } // Should be re-written for efficiency

private:
  CircularArray   _ca;

  unsigned int get_non_zero_entries_of_row(unsigned int i, unsigned int *indx, double *val) const;
};


} // End of namespace TOPUP
#endif
