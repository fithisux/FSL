//  
//  fsl_splines.h
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

#ifndef fsl_splines_h
#define fsl_splines_h

#include <vector>
#include <string>
#include <cmath>
#include "newmat.h"
#include "miscmaths/miscmaths.h"

namespace BASISFIELD {

class FslSplinesException: public std::exception
{
private:
  std::string m_msg;
public:
  FslSplinesException(const std::string& msg) throw(): m_msg(msg) {}

  virtual const char * what() const throw() {
    return string("FslSplines::" + m_msg).c_str();
  }

  ~FslSplinesException() throw() {}
};

// Declare ahead
template<class T>
class Spline1D;

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class Spline3D:
// This class is used by class splinefield to calculate field,
// derivative of field w.r.t. coefficients etc. It is implemented
// through Spline1D, containing a vector of three of those.  
//
// Its interface includes
// Operator (i,j,k) :
// Zero offset read access with range check.
// Peek(i,j,k):
// Zero offset read access without range check
// Operator [i] :
// Zero offset read access into vectorised spline without range check.
// RangeInField:
// Given a coef/spline index returns range of non-zero indicies in a 3D spline field
// OffsetIntoKernel:
// Given a coef/spline index returns first indicies in spline that falls within field
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

template<class T>
class Spline3D
{
public:
  // Constructors 
  Spline3D(unsigned int order, unsigned int ksp, std::vector<unsigned int> deriv = std::vector<unsigned int>(3,0)) : _vsp1D(3), _premul(false)
  {
    std::vector<unsigned int> vksp(3,ksp);
    common_construction(order,vksp,deriv);
  }
  Spline3D(unsigned int order, const std::vector<unsigned int>& ksp, std::vector<unsigned int> deriv = std::vector<unsigned int>(3,0)) : _vsp1D(3), _premul(false)
  {
    common_construction(order,ksp,deriv);
  }
  Spline3D(const Spline3D& sp) : _vsp1D(3)
  {
    common_construction(sp.Order(),sp.KnotSpacing(),sp.DerivativeOrder());
    if (sp.IsPremul()) {
      memcpy(_sp3D,sp._sp3D,sp.TotalKernelSize()*sizeof(T));
      _premul = true;
    }
  }
  // Destructor
  ~Spline3D() {delete [] _osp3D; delete [] _sp3D;}

  Spline3D<T>& operator=(const Spline3D<T>& insp)
  {
    if (this != &insp) {
      delete[] _osp3D;
      delete[] _sp3D;
      common_construction(insp.Order(),insp.KnotSpacing(),insp.DerivativeOrder());
      if (insp.IsPremul()) {
        memcpy(_sp3D,insp._sp3D,insp.TotalKernelSize()*sizeof(T));
        _premul = true;
      }
    }
    return(*this);
  }  

  // Functions for enquiring about properties of spline
  unsigned int Order() const {return(_vsp1D[0].Order());}
  unsigned int KnotSpacing(unsigned int dir) const {return(_vsp1D[dir].KnotSpacing());}
  std::vector<unsigned int> KnotSpacing() const {std::vector<unsigned int> vec(3,0); for (int i=0; i<3; i++) vec[i]=_vsp1D[i].KnotSpacing(); return(vec);}
  unsigned int DerivativeOrder(unsigned int dir) const {return(_vsp1D[dir].DerivativeOrder());}
  std::vector<unsigned int> DerivativeOrder() const {std::vector<unsigned int> vec(3,0); for (int i=0; i<3; i++) vec[i]=_vsp1D[i].DerivativeOrder(); return(vec);}
  unsigned int KernelSize(unsigned int dir) const {return(_vsp1D[dir].KernelSize());}
  unsigned int TotalKernelSize() const {return(_vsp1D[0].KernelSize()*_vsp1D[1].KernelSize()*_vsp1D[2].KernelSize());}
  unsigned int NCoef(unsigned int dir, unsigned int isz) const {return(_vsp1D[dir].NCoef(isz));}
  unsigned int FullFOV(unsigned int dir, unsigned int csz) const {return(_vsp1D[dir].FullFOV(csz));}
  unsigned int TotalFullFOV(const std::vector<unsigned int>& csz) const {return(_vsp1D[0].FullFOV(csz[0])*_vsp1D[1].FullFOV(csz[1])*_vsp1D[2].FullFOV(csz[2]));}
  bool IsPremul() const {return(_premul);}
 
  // Functions for read access to individual element when sampled on
  // a regular grid corresponding to voxels where knot-spacing is an
  // integer number of voxels
  T operator()(unsigned int i, unsigned int j, unsigned int k) const;  // Zero offset with range check
  T Peek(unsigned int i, unsigned int j, unsigned int k) const  // Zero offset without range check
  {
    return(_sp3D[k*_vsp1D[1].KernelSize()*_vsp1D[0].KernelSize()+j*_vsp1D[0].KernelSize()+i]);
  }
  T operator[](unsigned int indx) const {return(_sp3D[indx]);} // Zero offset vectorised access vithout range check

  // Functions for read access to individual elements when sampled
  // on any (non-regular) grid. Index should be in units of knot-spacings.
  T operator()(double x, double y, double z) const {return(_vsp1D[0](x) * _vsp1D[1](y) * _vsp1D[2](z));}
  T operator()(float x, float y, float z) const {return(this->operator()(static_cast<double>(x),static_cast<double>(y),static_cast<double>(z)));}

  // Function for read access to value of spline given by cindx
  // at voxel location given by vox.
  T SplineValueAtVoxel(const std::vector<double>&        vox,          // Voxel index
                       const std::vector<unsigned int>&  cindx) const  // Spline index
  {
    return(_vsp1D[2].SplineValueAtVoxel(vox[2],cindx[2]) * _vsp1D[1].SplineValueAtVoxel(vox[1],cindx[1]) * _vsp1D[0].SplineValueAtVoxel(vox[0],cindx[0]));
  }
  T SplineValueAtVoxel(const std::vector<float>&        vox,           // Voxel index
                       const std::vector<unsigned int>&  cindx) const  // Spline index
  {
    std::vector<double>  dvox(3,0);
    for (int i=0; i<3; i++) dvox[i] = static_cast<double>(vox[i]);
    return(SplineValueAtVoxel(dvox,cindx));
  }
  T SplineValueAtVoxel(const std::vector<unsigned int>&  vox,           // Voxel index
                       const std::vector<unsigned int>&  cindx) const  // Spline index
  {
    std::vector<double>  dvox(3,0);
    for (int i=0; i<3; i++) dvox[i] = static_cast<double>(vox[i]);
    return(SplineValueAtVoxel(dvox,cindx));
  }

  // Functions for access to entire spline kernel
  NEWMAT::ReturnMatrix AsNewmat() const;

  // Arithmetic operations on splines.

  Spline3D<T>& operator*=(double s);                             // Multiplication of self by scalar
  Spline3D<T>& operator/=(double s) {return((*this)*=(1.0/s));}  // Division of self by scalar

  // Functions returning the range of spline-kernels that has an overlap
  // with that given by cindx. The first form assumes that "this" spline
  // has the same knot-spacing as those it overlaps with. The second form
  // allows for passing a different spline by which "this" overlaps.
  // Note that I am using the C indexing convention here so that last
  // is one past. I.e. the range is first -- (last-1).

  bool RangeOfOverlappingSplines(const std::vector<unsigned int>&  cindx,
                                 const std::vector<unsigned int>&  isz,
                                 std::vector<unsigned int>&        first,
                                 std::vector<unsigned int>&        last) const;  

  bool RangeOfOverlappingSplines(const std::vector<unsigned int>&  cindx,
                                 const std::vector<unsigned int>&  isz,
                                 const Spline3D<T>&                sp2,
                                 std::vector<unsigned int>&        first,
                                 std::vector<unsigned int>&        last) const;  

  // The range functions will return the first and last index (zero-offset) 
  // into a field in which the spline is the cindx'th (zero-offset) spline.
  // Note that the range is really first -- (last-1).

  void RangeInField(const std::vector<unsigned int>& cindx,        // Index of spline/coefficient
                    const std::vector<unsigned int>& isz,          // Size of image/field
                    std::vector<unsigned int>&       first,        // First index in x-, y- and direction
                    std::vector<unsigned int>&       last) const   // Last index in x-, y- and direction
  {
    for (unsigned int i=0; i<3; i++) _vsp1D[i].Range(cindx[i],isz[i],first[i],last[i]);
  }         
  void XRangeInField(unsigned int cindx, unsigned int isz, unsigned int& first, unsigned int& last) const { _vsp1D[0].Range(cindx,isz,first,last); }
  void YRangeInField(unsigned int cindx, unsigned int isz, unsigned int& first, unsigned int& last) const { _vsp1D[1].Range(cindx,isz,first,last); }
  void ZRangeInField(unsigned int cindx, unsigned int isz, unsigned int& first, unsigned int& last) const { _vsp1D[2].Range(cindx,isz,first,last); }

  // The RangeOfSplines functions will return the first and (one past) last indicies 
  // of spline kernels that have a support (non-zero) value at a voxel location.
  
  void RangeOfSplines(const std::vector<double>&        vox,         // Voxel index (does not need to be integer)
                      const std::vector<unsigned int>&  csz,         // Coefficient matrix size
                      std::vector<unsigned int>&        first,       // Index of first spline with support for vox
                      std::vector<unsigned int>&        last) const  // Index of last spline with support for vox
  {
    for (unsigned int i=0; i<3; i++) _vsp1D[i].RangeOfSplines(vox[i],csz[i],first[i],last[i]);
  }
  void RangeOfSplines(const std::vector<float>&         vox,         // Voxel index (does not need to be integer)
                      const std::vector<unsigned int>&  csz,         // Coefficient matrix size
                      std::vector<unsigned int>&        first,       // Index of first spline with support for vox
                      std::vector<unsigned int>&        last) const  // Index of last spline with support for vox
  {
    std::vector<double>  dvox(3,0);
    for (int i=0; i<3; i++) dvox[i] = static_cast<double>(vox[i]);
    RangeOfSplines(dvox,csz,first,last);
  }
  void RangeOfSplines(const std::vector<unsigned int>&  vox,         // Voxel index
                      const std::vector<unsigned int>&  csz,         // Coefficient matrix size
                      std::vector<unsigned int>&        first,       // Index of first spline with support for vox
                      std::vector<unsigned int>&        last) const  // Index of last spline with support for vox
  {
    std::vector<double>  dvox(3,0);
    for (int i=0; i<3; i++) dvox[i] = static_cast<double>(vox[i]);
    RangeOfSplines(dvox,csz,first,last);
  }

  // The offset functions returns the offset into the spline kernel to the first
  // index in the spline that falls inside the field.

  void OffsetIntoKernel(const std::vector<unsigned int>&  indx,        // Index of spline/coefficient
                        const std::vector<unsigned int>&  isz,         // Size of image/field
                        std::vector<unsigned int>&        offs) const  // Offset in x-, y- and z-direction
  {
    for (unsigned int i=0; i<3; i++) offs[i] = _vsp1D[i].Offset(indx[i],isz[i]);
  }
  unsigned int XOffsetIntoKernel(unsigned int indx, unsigned int isz) const { return(_vsp1D[0].Offset(indx,isz)); } 
  unsigned int YOffsetIntoKernel(unsigned int indx, unsigned int isz) const { return(_vsp1D[1].Offset(indx,isz)); } 
  unsigned int ZOffsetIntoKernel(unsigned int indx, unsigned int isz) const { return(_vsp1D[2].Offset(indx,isz)); } 

  // The premul function is useful for multiplying a spline with a field once and for all
  // before multiplying it with other splines. Note that all the read-access methods 
  // returns the pre-multiplied values after a call to this function.

  template <class S>
  void Premul(const std::vector<unsigned int>&  indx,       // Index of spline/coefficient into field of same size as ima
              const std::vector<unsigned int>&  isz,        // Size of ima
              const S                           *ima);      // Well.

  // The NzMax routines calculates the maximum # of non-zero elements of A'*B where
  // A and B are matrices with one column per spline and each column containing a 
  // spline kernel. The first form assume identical splines in A and B.

  unsigned int NzMax(const std::vector<unsigned int>&  isz) const;
  unsigned int NzMax(const std::vector<unsigned int>&  isz,
                     const Spline3D<T>&                sp2) const;

  T MulByOther(const std::vector<unsigned int>&  cindx1,
               const std::vector<unsigned int>&  cindx2,
               const std::vector<unsigned int>&  isz,
               const Spline3D<T>&                sp2) const;

private:
  std::vector<Spline1D<T> >   _vsp1D;       // Vector of 1D splines
  T                           *_osp3D;      // 3D spline kernel
  T                           *_sp3D;       // 3D spline kernel, possibly pre-multiplied
  bool                        _premul;      // True if _sp3D has been pre-multiplied by image (or scalar)

  void common_construction(unsigned int order, const std::vector<unsigned int>& ksp, const std::vector<unsigned int>& deriv);

};

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class Spline1D:
// This is for all purposes a helper class for Spline3D. 
// It implements some basic characteristics of 1D splines.
// Its interface includes
//
// Operator () :
// One offset read access with range check.
// Peek:
// One offset read access without range check
// Operator [] :
// Zero offset read access without range check.
// Range:
// Given a coef/spline index returns range of non-zero indicies in a 1D spline field
// Offset:
// Given a coef/spline index returns first index in spline that falls within field
// GetAMatrix:
// Returns A such that A*coef gives 1D field, where coef ColumnVector of coefficients
// GetMMatrix:
// Returns M that can be used to zoom up and down between knot-spacings. 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

template<class T>
class Spline1D
{
public:
  Spline1D(unsigned int order, unsigned int ksp, unsigned int deriv=0) : _order(order), _ksp(ksp), _deriv(deriv), _n(spline_size()), _sp(spline_kernel()) 
  {
    if (_order < 2 || _order > 3) throw FslSplinesException("Spline1D::Spline1D: only orders 2 and 3 implemented");
    if (_ksp < 1) throw FslSplinesException("Spline1D::Spline1D: knot-spacing has to be > 0");
    if (_deriv > (_order)) throw FslSplinesException("Spline1D::Spline1D: only differentiable order times");
  }
  Spline1D(const Spline1D& insp) : _order(insp._order), _ksp(insp._ksp), _deriv(insp._deriv), _n(insp._n), _sp(spline_kernel()) {}    
  Spline1D() : _order(3), _ksp(1), _deriv(0), _n(spline_size()), _sp(spline_kernel()) {}
  ~Spline1D() {delete[] _sp;}
  
  Spline1D<T>& operator=(const Spline1D<T>& insp)
  {
    if (this != &insp) { 
      _order = insp._order; _ksp = insp._ksp; _deriv = insp._deriv; _n = insp._n;
      delete[] _sp; _sp = spline_kernel();
    }
    return(*this);
  }

  // Functions for enquring about properties of spline

  unsigned int Order() const {return(_order);}
  unsigned int KnotSpacing() const {return(_ksp);}
  unsigned int DerivativeOrder() const {return(_deriv);}
  unsigned int KernelSize() const {return(_n);}
  unsigned int FullFOV(unsigned int csz) const;
  unsigned int NCoef(unsigned int isz) const;

  // Returns the offset (in the field, zero-offset) for the center
  // of the specified spline (default is first spline).
  int OffsetIntoField(unsigned int cindx=0) const;

  // Returns the range of voxels for which cindx has 
  // support (i.e. non-zero value).
  // Note that range is really first -- (last-1).
  void Range(// Input
             unsigned int   cindx,       // Index of coefficient
             unsigned int   isz,         // Size of image/field
             // Output
             unsigned int&  first,       // First index (into field) of spline 
             unsigned int&  last) const; // Last index (into field) of spline

  // Returns the range of splines that has support
  // (i.e. non-zero value) for vox.
  // Note that the range is really first -- (last-1).
  void RangeOfSplines(double         vox,           // Voxel index (can be non-integer)
                      unsigned int   csz,           // Coefficient matrix size
                      unsigned int&  first,         // Index of first spline with support for vox
                      unsigned int&  last) const;   // Index of last spline with support for vox

  // Returns first index (into spline-kernel that falls within field/image
  unsigned int Offset(unsigned int   cindx,       // Index of coefficient 
		      unsigned int   isz) const;  // Image/field-size 

  // Returns the rance of spline/coefficient indicies that overlap with *this
  // Note that range is really first -- (last-1)
  bool RangeOfOverlappingSplines(// Input
                                 unsigned int   cindx,       // Index of coefficient/spline
                                 unsigned int   isz,         // Size of image/field
                                 // Output
                                 unsigned int&  first,       // First coefficient/spline with overlap 
                                 unsigned int&  last) const; // Last coefficient/spline with overlap

  bool RangeOfOverlappingSplines(// Input
				 unsigned int          cindx,       // Index of coefficient/spline
                                 unsigned int          isz,         // Size of image/field
                                 const Spline1D<T>&    sp2,         // Type of spline that we want to know overlap with
                                 // Output
                                 unsigned int&         first,       // First coefficient/spline with overlap 
                                 unsigned int&         last) const; // Last coefficient/spline with overlap

  // Returns the value of spline given by cindx at voxel vox
  T SplineValueAtVoxel(double         vox,         // Voxel index
                       unsigned int   cindx) const // Spline coefficient index
  {
    return(kernel_value(vox2splindx(vox,cindx),_order,_ksp,_deriv));
  }

  // Functions for read access to individual elements when sampled
  // on a regular grid corresponding to voxels where knot-spacing
  // is an integer number of voxels.
  T operator()(unsigned int indx) const // One offset with range check
  {
    if (indx>0 && indx<(_n+1)) return(_sp[indx-1]);
    else throw FslSplinesException("Spline1D::operator(): index out of range");
  }
  T Peek(unsigned int indx) const {return(_sp[indx-1]);}      // One offset without range check
  T operator[](unsigned int indx) const {return(_sp[indx]);}  // Zero offset without range check

  // Functions for read access to individual elements when sampled
  // on any (non-regular) grid. Index should be in units of knot-spacings.
  T operator()(double x) const {return(this->kernel_value(x,_order,_ksp,_deriv));}
  T operator()(float x) const {return(this->operator()(static_cast<double>(x)));}

  // Functions for obtaining reference to/copy of entire kernel
  const T* Kernel() const {return(_sp);}
  NEWMAT::ReturnMatrix AsNewmat() const;

  // Functions for generating/zooming 1D fields 
  NEWMAT::ReturnMatrix GetAMatrix(unsigned int isz,                    // Image/field size
                                  unsigned int csz=0) const;           // # of splines/coefficients

  NEWMAT::ReturnMatrix GetMMatrix(const Spline1D<T>&  s,               // Spline for which we want coefficients
                                  unsigned int        isz,             // Size of image/field
                                  unsigned int        csz1=0,          // # of coefficients we know for this spline
                                  unsigned int        csz2=0) const;   // # of coefficients we want to know for s

private:
  unsigned int                 _order;
  unsigned int                 _ksp;
  unsigned int                 _deriv;
  unsigned int                 _n;
  T                            *_sp;  

  unsigned int spline_size(unsigned int order=0, unsigned int ksp=0) 
  {
    if (!order) order = _order;
    if (!ksp) ksp = _ksp;
    if (!(order%2) && ksp%2) {   // If even order (e.g. quadratic splines) and odd knot-spacing
      return((order+1)*ksp);
    }
    else { 
      return((order+1)*ksp-1); 
    }
  }

  double vox2splindx(double vox, unsigned int cindx) const
  {
    return((vox-OffsetIntoField(cindx))/static_cast<double>(_ksp));
  }

  T* spline_kernel(unsigned int order=0, unsigned int ksp=0, unsigned int deriv=0);

  T kernel_value(double x, unsigned int order, unsigned int ksp, unsigned int deriv) const;
};


/////////////////////////////////////////////////////////////////////
//
// Here starts memeber functions for Spline3D
//
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
//
// Multiplication by scalar
//
/////////////////////////////////////////////////////////////////////

template<class T>
Spline3D<T>& Spline3D<T>::operator*=(double s)
{
  if (_premul) throw FslSplinesException("Spline3D::operator*=: You should not combine multiplication by scalar and image");

  for (unsigned int i=0; i<TotalKernelSize(); i++) _sp3D[i] *= s;

  return(*this);
}

// N.B. Next three functions non-members

template<class T>
Spline3D<T> operator*(const Spline3D<T>& sp, double s)
{
  Spline3D<T>   spo(sp);
  spo *= s;
  return(spo);
}

template<class T>
Spline3D<T> operator*(double s, const Spline3D<T>& sp)
{
  return(sp*s);
}

template<class T>
Spline3D<T> operator/(const Spline3D<T>& sp, double s)
{
  Spline3D<T>  spo(sp);
  spo /= s;
  return(spo);
}

/////////////////////////////////////////////////////////////////////
//
// Zero-offset read-only access with range check
//
/////////////////////////////////////////////////////////////////////

template<class T>
inline T Spline3D<T>::operator()(unsigned int i, unsigned int j, unsigned int k) const
{
  if (i<_vsp1D[0].KernelSize() && j<_vsp1D[1].KernelSize() && k<_vsp1D[2].KernelSize()) {
    return(_sp3D[k*_vsp1D[1].KernelSize()*_vsp1D[0].KernelSize()+j*_vsp1D[0].KernelSize()+i]);
  }
  else {
    throw FslSplinesException("Spline3D::operator(): index out of range");
  }
}

/////////////////////////////////////////////////////////////////////
//
// Return entire 3D kernel in NEWMAT matrix format
//
/////////////////////////////////////////////////////////////////////

template<class T>
NEWMAT::ReturnMatrix Spline3D<T>::AsNewmat() const
{
  NEWMAT::Matrix   M(KernelSize(1)*KernelSize(2),KernelSize(0));

  for (unsigned int k=0; k<KernelSize(2); k++) {
    for (unsigned int j=0; j<KernelSize(1); j++) {
      for (unsigned int i=0; i<KernelSize(0); i++) {
        M(k*KernelSize(1)+j+1,i+1) = _sp3D[k*KernelSize(1)*KernelSize(0)+j*KernelSize(0)+i];
      }
    }
  }
  M.Release();
  return(M);
}
/////////////////////////////////////////////////////////////////////
//
// Sum of pair-wise products of all overlapping voxels of
// *this and sp2.
//
/////////////////////////////////////////////////////////////////////

template<class T>
T Spline3D<T>::MulByOther(const std::vector<unsigned int>&  cindx1,
                          const std::vector<unsigned int>&  cindx2,
                          const std::vector<unsigned int>&  isz,
                          const Spline3D<T>&                sp2) const
{
  std::vector<unsigned int>   f1(3,0);      // First index into *this of overlapping part
  std::vector<unsigned int>   f2(3,0);      // First index into sp2 of overlapping part
  std::vector<unsigned int>   l1(3,0);      // Last index into *this of overlapping part

  unsigned int if1 = 0, if2 = 0;            // First index given cindx1[i] and cindx2[i] in (some) ima of size isz[i]
  unsigned int il1 = 0, il2 = 0;            // Last index given cindx1[i] and cindx2[i] in (some) ima of size isz[i]
  unsigned int o1 = 0, o2 = 0;              // Offset into spline kernel for cindx1[i] and cindx2[i] respectively
  T            prod  = static_cast<T>(0.0); // Product
  
  for (unsigned int d=0; d<3; d++) {
    _vsp1D[d].Range(cindx1[d],isz[d],if1,il1);
    o1 = _vsp1D[d].Offset(cindx1[d],isz[d]);
    sp2._vsp1D[d].Range(cindx2[d],isz[d],if2,il2);
    o2 = sp2._vsp1D[d].Offset(cindx2[d],isz[d]);
    f1[d] = max(static_cast<int>(if2)-static_cast<int>(if1),0) + o1;
    f2[d] = max(static_cast<int>(if1)-static_cast<int>(if2),0) + o2;
    l1[d] = min(il1,il2) - if1 + o1;
  }

  for (unsigned int k1=f1[2], k2=f2[2]; k1<l1[2]; k1++, k2++) {
    unsigned int bi1 = k1*KernelSize(1)*KernelSize(0); 
    unsigned int bi2 = k2*sp2.KernelSize(1)*sp2.KernelSize(0); 
    for (unsigned int j1=f1[1], j2=f2[1]; j1<l1[1]; j1++, j2++) {
      unsigned int bi12 = bi1 + j1*KernelSize(0);
      unsigned int bi22 = bi2 + j2*sp2.KernelSize(0);
      for (unsigned int i1=f1[0], i2=f2[0]; i1<l1[0]; i1++, i2++) {
        prod += this->operator[](bi12+i1) * sp2[bi22+i2];
      }
    }
  }
  return(prod);
}

/////////////////////////////////////////////////////////////////////
//
// Get maximum # of non-zero elements in A'*B.
//
/////////////////////////////////////////////////////////////////////

template<class T>
unsigned int Spline3D<T>::NzMax(const std::vector<unsigned int>&  isz) const
{
  unsigned int nzmax=0;

  for (unsigned int k=0; k<_vsp1D[2].NCoef(isz[2]); k++) {
    unsigned int first=0, last=0;
    _vsp1D[2].RangeOfOverlappingSplines(k,isz[2],first,last);
    unsigned int zn = last-first;
    for (unsigned int j=0; j<_vsp1D[1].NCoef(isz[1]); j++) {
      _vsp1D[1].RangeOfOverlappingSplines(j,isz[1],first,last);
      unsigned int yn = last-first;
      for (unsigned int i=0; i<_vsp1D[0].NCoef(isz[0]); i++) {
        _vsp1D[0].RangeOfOverlappingSplines(i,isz[0],first,last);
        nzmax += zn*yn*(last-first);
      }
    }
  }
  return(nzmax);    
}

template<class T>
unsigned int Spline3D<T>::NzMax(const std::vector<unsigned int>&  isz,
                                const Spline3D&                   sp2) const
{
  unsigned int nzmax=0;

  for (unsigned int k=0; k<_vsp1D[2].NCoef(isz[2]); k++) {
    unsigned int first=0, last=0;
    _vsp1D[2].RangeOfOverlappingSplines(k,isz[2],sp2._vsp1D[2],first,last);
    unsigned int zn = last-first;
    for (unsigned int j=0; j<_vsp1D[1].NCoef(isz[1]); j++) {
      _vsp1D[1].RangeOfOverlappingSplines(j,isz[1],sp2._vsp1D[1],first,last);
      unsigned int yn = last-first;
      for (unsigned int i=0; i<_vsp1D[0].NCoef(isz[0]); i++) {
        _vsp1D[0].RangeOfOverlappingSplines(i,isz[0],sp2._vsp1D[0],first,last);
        nzmax += zn*yn*(last-first);
      }
    }
  }
  return(nzmax);    
}

/////////////////////////////////////////////////////////////////////
//
// Get range of indicies of splines that has an overlap with *this
//
/////////////////////////////////////////////////////////////////////

template<class T>
bool Spline3D<T>::RangeOfOverlappingSplines(const std::vector<unsigned int>&  cindx,
                                            const std::vector<unsigned int>&  isz,
                                            std::vector<unsigned int>&        first,
                                            std::vector<unsigned int>&        last) const
{
  for (unsigned int i=0; i<3; i++) { 
    if (!_vsp1D[i].RangeOfOverlappingSplines(cindx[i],isz[i],first[i],last[i])) return(false);
  }
  return(true);
}

template<class T>
bool Spline3D<T>::RangeOfOverlappingSplines(const std::vector<unsigned int>&  cindx,
                                            const std::vector<unsigned int>&  isz,
                                            const Spline3D<T>&                sp2,
                                            std::vector<unsigned int>&        first,
                                            std::vector<unsigned int>&        last) const  
{
  for (unsigned int i=0; i<3; i++) { 
    if (!_vsp1D[i].RangeOfOverlappingSplines(cindx[i],isz[i],sp2._vsp1D[i],first[i],last[i])) return(false);
  }
  return(true);
}

/////////////////////////////////////////////////////////////////////
//
// Pre-multiply spline with image.
//
/////////////////////////////////////////////////////////////////////

template<class T>
template<class S>
void Spline3D<T>::Premul(const std::vector<unsigned int>&  cindx,      // Index of spline/coefficient into field of same size as ima
                         const std::vector<unsigned int>&  isz,        // Size of ima
                         const S                           *ima)       // Well.
{
  std::vector<unsigned int>    fi(3,0);    // First indicies into ima
  std::vector<unsigned int>    li(3,0);    // Last indicies into ima
  std::vector<unsigned int>    so(3,0);    // Offset into spline

  RangeInField(cindx,isz,fi,li);
  OffsetIntoKernel(cindx,isz,so);
  memset(_sp3D,0,TotalKernelSize()*sizeof(T));

  for (unsigned int ik=fi[2], sk=so[2]; ik<li[2]; ik++, sk++) {
    for (unsigned int ij=fi[1], sj=so[1]; ij<li[1]; ij++, sj++) {
      unsigned int ibi = ik*isz[1]*isz[0] + ij*isz[0];
      unsigned int sbi = sk*KernelSize(1)*KernelSize(0) + sj*KernelSize(0);
      for (unsigned int ii=fi[0], si=so[0]; ii<li[0]; ii++, si++) {
        _sp3D[sbi+si] = _osp3D[sbi+si] * static_cast<T>(ima[ibi+ii]);
      }
    }
  }
  _premul = true;         
}

/////////////////////////////////////////////////////////////////////
//
// Private helper function for constructors of Spline3D
//
/////////////////////////////////////////////////////////////////////

template<class T>
void Spline3D<T>::common_construction(unsigned int order, const std::vector<unsigned int>& ksp, const std::vector<unsigned int>& deriv)
{
  unsigned int size = 1;
  for (unsigned int i=0; i<3; i++) {
    _vsp1D[i] = Spline1D<T>(order,ksp[i],deriv[i]);
    size *= _vsp1D[i].KernelSize();
  }
  _osp3D = new T[size];
  for (unsigned int k=0;k<_vsp1D[2].KernelSize(); k++) {
    for (unsigned int j=0; j<_vsp1D[1].KernelSize(); j++) {
      unsigned int bi = k*_vsp1D[1].KernelSize()*_vsp1D[0].KernelSize()+j*_vsp1D[0].KernelSize();
      for (unsigned int i=0; i<_vsp1D[0].KernelSize(); i++) {
        _osp3D[bi+i] = _vsp1D[2][k]*_vsp1D[1][j]*_vsp1D[0][i];
      }
    }
  }
  _sp3D = new T[size];
  memcpy(_sp3D,_osp3D,TotalKernelSize()*sizeof(T)); 
}

/////////////////////////////////////////////////////////////////////
//
// Here starts memeber functions for Spline1D
//
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
//
// Returns the range of indicies of splines that overlap with self,
// when the other splines are of same type as self.
//
/////////////////////////////////////////////////////////////////////

template<class T>
inline bool Spline1D<T>::RangeOfOverlappingSplines(// Input
                                                   unsigned int   cindx,       // Index of coefficient/spline
                                                   unsigned int   isz,         // Size of image/field
                                                   // Output
                                                   unsigned int&  first,       // First coefficient/spline with overlap 
                                                   unsigned int&  last) const  // One past last coefficient/spline with overlap
{
  if (_ksp==1) {
    first = max(0,static_cast<int>(cindx) - static_cast<int>(2*(_order/2)));
    last = min(NCoef(isz),cindx+2*(_order/2)+1);
  }
  else {
    first = max(0,static_cast<int>(cindx) - static_cast<int>(_order));
    last = min(NCoef(isz),cindx+_order+1);
  }
  if (last > first) return(true);
  else return(false);
}

template<class T>
bool Spline1D<T>::RangeOfOverlappingSplines(// Input
				            unsigned int          cindx,       // Index of coefficient/spline
                                            unsigned int          isz,         // Size of image/field
                                            const Spline1D<T>&    sp2,         // Type of spline that we want to know overlap with
                                            // Output
                                            unsigned int&         first,       // First coefficient/spline with overlap 
                                            unsigned int&         last) const  // Last coefficient/spline with overlap
{
  unsigned int fme, lme;       // First and last index of *this into ima of size isz
  unsigned int fsp2, lsp2;     // First and lasp index of sp2 into ima of size isz
  int          ifirst;         // Signed version of first

  Range(cindx,isz,fme,lme);
  // cout << "Range is " << fme << " -- " << lme << endl;
  ifirst = -1; 
  for (unsigned int i=0; i<sp2.NCoef(isz); i++) {
    sp2.Range(i,isz,fsp2,lsp2);
    // cout << "sp2.Range is " << fsp2 << " -- " << lsp2 << endl;
    int overlap = min(lme-1,lsp2-1) - max(fme,fsp2);
    if (overlap >= 0) {  // If there is at least one common voxel
      if (ifirst == -1) ifirst = i;
      else last = i+1;
    }
    else {               // If there is no overlap
      if (ifirst != -1) break;
    }
  }

  if (ifirst == -1) {    // If there were no overlap at all
    first = 0; last = 0; // Return zero range
    return(false);
  }
  else {
    first = static_cast<unsigned int>(ifirst);   
    return(true);
  }
}

/////////////////////////////////////////////////////////////////////
//
// Returns the first and (one past) last index into a zero-offset 
// field for the spline with index cindx.
//
/////////////////////////////////////////////////////////////////////

template<class T>
inline void Spline1D<T>::Range(// Input
                               unsigned int   cindx,       // Index of coefficient
                               unsigned int   isz,         // Size of image/field
                               // Output
                               unsigned int&  first,       // First index (into field) of spline 
                               unsigned int&  last) const  // One past last index (into field) of spline
{
  switch(_order) {
  case 2: case 3:
    first = std::max(0,OffsetIntoField(cindx)-static_cast<int>(KernelSize()/2));
    last = std::min(static_cast<int>(isz),OffsetIntoField(cindx)+static_cast<int>(KernelSize()/2)+1);
    //  cout << "New range = " << first << " -- " << last << endl;
    //  first = std::max(0,(static_cast<int>(cindx)-3)*static_cast<int>(_ksp) + 1);
    //  last = min(isz,(cindx+1)*_ksp - 1);  // Scary change made 21/2-08. J.A.
    //  last = std::min(isz,(cindx+1)*_ksp);
    //  cout << "Old range = " << first << " -- " << last << endl;
    break;
  default:
    throw FslSplinesException("Spline1D::Range: Only orders 2 and 3 implemented");
  }
}

/////////////////////////////////////////////////////////////////////
//
// Returns the first and (one past) last coefficient/spline index 
// that has support for voxel vox.
//
/////////////////////////////////////////////////////////////////////

template<class T>
inline void Spline1D<T>::RangeOfSplines(// Input 
                                        double         vox,           // Voxel index (can be non-integer)
                                        unsigned int   csz,           // Coefficient matrix size
                                        // Output
                                        unsigned int&  first,         // Index of first spline with support for vox
                                        unsigned int&  last) const   // Index of last spline with support for vox
{
  first = std::max(0,static_cast<int>(ceil((vox - static_cast<double>(OffsetIntoField(0)) - (static_cast<double>(KernelSize())+1.0)/2.0) / _ksp)));
  int tmp = std::min(static_cast<int>(csz),static_cast<int>(ceil((vox - static_cast<double>(OffsetIntoField(0)) + (static_cast<double>(KernelSize())+1.0)/2.0) / _ksp)));
  if (tmp < 0) last = 0;
  else last = static_cast<unsigned int>(tmp);
}

/////////////////////////////////////////////////////////////////////
//
// Returns the offset (in the field, zero-offset) for the center 
// (middle) element of the specified spline. 
//
/////////////////////////////////////////////////////////////////////

template<class T>
inline int Spline1D<T>::OffsetIntoField(unsigned int cindx) const   
{
  int oif = 0;
  if (_ksp == 1) oif = cindx*_ksp;
  else oif = (cindx - _order/2)*_ksp;
  return(oif);
}

/////////////////////////////////////////////////////////////////////
//
// Returns the offset (in the spline) for the first 
// element of the spline that falls within the field.
//
/////////////////////////////////////////////////////////////////////

template<class T>
inline unsigned int Spline1D<T>::Offset(unsigned int   cindx,       // Index of coefficient 
					unsigned int   isz) const   // Image/field-size 
{
  unsigned int os=0;
  switch(_order) {
  case 2: case 3:
    os = -std::min(0,OffsetIntoField(cindx)-static_cast<int>(KernelSize())/2);
    // cout << "New Offset = " << os << endl;
    // os = -std::min(0,(static_cast<int>(cindx)-3)*static_cast<int>(_ksp) + 1);
    // cout << "Old Offset = " << os << endl;
    break;
  default:
    throw FslSplinesException("Spline1D::Offset: Only orders 2 and 3 implemented");
  }
  return(os);
}

/////////////////////////////////////////////////////////////////////
//
// Returns a copy of the spline kernel as a NEWMAT ColumnVector.
//
/////////////////////////////////////////////////////////////////////

template<class T>
NEWMAT::ReturnMatrix Spline1D<T>::AsNewmat() const
{
  NEWMAT::ColumnVector s(_n);
  for (unsigned int i=0; i<_n; i++) s.element(i) = static_cast<double>(_sp[i]);
  s.Release();
  return(s);
}


/////////////////////////////////////////////////////////////////////
//
// Returns the "full" FOV (int # of voxels) for csz # of splines.
// With "full" I mean all those voxels for which there is a non-zero
// representation of any spline. 
//
/////////////////////////////////////////////////////////////////////

template<class T>
unsigned int Spline1D<T>::FullFOV(unsigned int csz) const
{
  if (KernelSize() % 2) { // If kernel size odd # of voxels
    return((csz-1) * KnotSpacing() + 1 + KernelSize() - 1);
  }
  else {  // If kernel size even # of voxels
    return((csz-1) * KnotSpacing() + 1 + KernelSize());
  } 
}

/////////////////////////////////////////////////////////////////////
//
// Returns the smallest number of coefficients/splines neccesary
// to completely specify a 1D spline function over the range 1-isz
// (or 0-(isz-1) if you prefer).
// Quite infuriatingly I am using a value that is sometimes 
// non-optimal for _order=3. This is in order to ensure backwards
// compatibilty of some applications. The "correct" value is that
// which is given in the else branch.
//
/////////////////////////////////////////////////////////////////////

template<class T>
unsigned int Spline1D<T>::NCoef(unsigned int isz) const
{
  unsigned int nc=0;

  if (_ksp==1) nc = isz;
  else if (_order==3) {
    nc = static_cast<unsigned int>(ceil(double(isz+1)/double(_ksp))) + 2;
  }
  else {
    nc = isz / _ksp + _order/2;
    while (OffsetIntoField(nc)-KernelSize()/2 < isz-1) nc++;
  }

  return(nc);
}

/////////////////////////////////////////////////////////////////////
//
// Returns matrix M that allows you to map the coefficients for a 
// 1D spline function from one knot-spacing to another. It is intended
// for use in zooming (in and out) splinefields. Let us say we have
// a spline sp1 and a set of coefficients c1 giving a 1D function
// over a range 1-64 through f = sp1.A(64)*c1; Let us further say
// that we have another spline (with a different knot-spacing) sp2
// and that we want to find the coefficients c2 such that 
// f = sp2.A(64)*c2; We can then find c2 from c2 = sp2.M(sp1,64)*c1;
// If sp2.KnotSpacing()<=sp1.KnotSpacing() the result will be "exact",
// i.e. f is identical. If sp2.KnotSpacing()>sp1.KnotSpacing() the 
// result is optimal in a least-squares sense.
//
/////////////////////////////////////////////////////////////////////

template<class T>
NEWMAT::ReturnMatrix Spline1D<T>::GetMMatrix(const Spline1D<T>&  s, 
                                             unsigned int        isz,
                                             unsigned int        csz1,
                                             unsigned int        csz2) const
{
  if (!csz1) csz1 = NCoef(isz);
  if (!csz2) csz2 = s.NCoef(isz);
  NEWMAT::Matrix A1 = GetAMatrix(isz,csz1);
  NEWMAT::Matrix A2 = s.GetAMatrix(isz,csz2);
  // NEWMAT::Matrix M = (A1.t()*A1).i() * (A1.t()*A2);
  NEWMAT::Matrix M = MISCMATHS::pinv(A1.t()*A1) * (A1.t()*A2);
  M.Release();
  return(M);
}

/////////////////////////////////////////////////////////////////////
//
// Returns matrix A that allows you to calculate a 1D spline function
// f over a range 1-isz given a set of coefficients c from f = A*c,
// where f and c are NEWMAT::ColumnVector's.
//
/////////////////////////////////////////////////////////////////////


template<class T>
NEWMAT::ReturnMatrix Spline1D<T>::GetAMatrix(unsigned int isz,
					     unsigned int csz) const
{
  if (!csz) csz = NCoef(isz);
  NEWMAT::Matrix A(isz,csz);
  A = 0.0;

  unsigned int fr=0, lr=0;
  unsigned int so=0;
  for (unsigned int i=0; i<csz; i++) {   // Loop over splines
    Range(i,isz,fr,lr);                  // Get range of row indicies
    so = Offset(i,isz);
    for (unsigned int j=fr, si=so; j<lr; j++, si++) {
      A(j+1,i+1) = (*this)[si];
    }
  }
  A.Release();
  return(A); 
}


/*
template<class T>
NEWMAT::ReturnMatrix Spline1D<T>::GetAMatrix(unsigned int isz, 
                                             unsigned int csz) const
{
  if (!csz) csz = NCoef(isz);
  NEWMAT::Matrix A(isz,csz);
  A = 0.0;
  for (unsigned int r=1; r<=isz; r++) {
    if ((r % _ksp) == 1) {
      for (int i=1; i<4; i++) {
        unsigned int c = (r / _ksp) + i;          // Intentional integer division
	if (c>=1 && c<=csz) {
          A(r,c) = _sp[_n - i*_ksp];
	}
      }
    }
    else {
      for (int i=0; i<4; i++) {
	unsigned int c = ((r-1) / _ksp) + i + 1;  // Intentional integer division
	if (c>=1 && c<=csz) {
	  A(r,c) = _sp[_n - (i+1)*_ksp + ((r-1)%_ksp)];
	}
      }
    }
  }
  A.Release();
  return(A);    
}
*/

/////////////////////////////////////////////////////////////////////
//
// Internal (private) function that returns the value of a spline
// for the abscissa x. x should be in units of knot-spacings from
// the centre of the spline. 
//
/////////////////////////////////////////////////////////////////////

template<class T>
T Spline1D<T>::kernel_value(double x, unsigned int order, unsigned int ksp, unsigned int deriv) const
{
  double ax = fabs(x);

  switch (order) {
  case 0:
  case 1:
    switch (deriv) {
    case 0:
      if (ax < 1) return(static_cast<T>(1 - ax));
      else return(static_cast<T>(0.0));
      // break; // Commented away to stop warnings from nvcc compiler
    case 1:
      if (ax < 1e-6) return(static_cast<T>(0.0));                                        // Smack on center discontinuity
      else if (fabs(ax-1.0) < 1e-6) return(static_cast<T>(-(x/ax) * 0.5 / double(ksp))); // Smack on edge discontinuity
      else if (ax < 1.0) return(static_cast<T>(-(x/ax) / double(ksp)));
      break;
    }
    break;
  case 2:
    switch (deriv) {
    case 0:
      if (ax <= 0.5) return(static_cast<T>(0.75 - ax*ax));
      else if (ax < 1.5) return(static_cast<T>(0.5 * (1.5-ax) * (1.5-ax)));
      else return(static_cast<T>(0.0));
      // break; // Commented away to stop warnings from nvcc compiler
    case 1:
      if (ax < 1e-6) return(static_cast<T>(0.0));
      else if (ax <= 0.5) return(static_cast<T>(- (x/ax) * 2.0 * (ax / double(ksp))));
      else if (ax < 1.5) return(static_cast<T>((x/ax) * 0.5 * ((2.0*ax - 3)/double(ksp))));
      else return(static_cast<T>(0.0));
      // break; // Commented away to stop warnings from nvcc compiler
    case 2:
      if (fabs(ax-1.5) < 1e-6) return(static_cast<T>(0.5 / double(ksp*ksp)));          // Smack on edge discontinuity
      else if (fabs(ax-0.5) < 1e-6) return(static_cast<T>(-0.5 / double(ksp*ksp)));    // Smack on inner discontinuity
      else if (ax < 0.5) return(static_cast<T>(-2.0 / double(ksp*ksp)));
      else if (ax < 1.5) return(static_cast<T>(1.0 / double(ksp*ksp)));
      else return(static_cast<T>(0.0));
      // break; // Commented away to stop warnings from nvcc compiler
    default:
      break;
    }
    break;
  case 3:
    switch (deriv) {
    case 0:
      if (ax <= 1.0) return(static_cast<T>((2.0/3.0) + (ax*ax)*(ax/2.0 - 1)));
      else if (ax < 2.0) return(static_cast<T>((1.0/6.0) * (2.0-ax)*(2.0-ax)*(2.0-ax)));
      else return(static_cast<T>(0.0));
      // break; // Commented away to stop warnings from nvcc compiler
    case 1:
      if (ax < 1e-6) return(static_cast<T>(0.0));
      else if (ax <= 1.0) return(static_cast<T>((x/ax) * (1.5*ax*ax - 2.0*ax) / double(ksp)));
      else if (ax < 2.0) return(static_cast<T>((x/ax) * (-0.5*(2.0-ax)*(2.0-ax)) / double(ksp)));
      else return(static_cast<T>(0.0));
      // break; // Commented away to stop warnings from nvcc compiler
    case 2:
      if (ax <= 1.0) return(static_cast<T>((3.0*ax - 2.0)/double(ksp*ksp)));
      else if (ax < 2.0) return(static_cast<T>((2.0-ax)/double(ksp*ksp)));
      else return(static_cast<T>(0.0));
      // break; // Commented away to stop warnings from nvcc compiler
    default:
      break;
    }
    break;
  default:
    break; 
  }
  return(0.0);
}

/////////////////////////////////////////////////////////////////////
//
// Internal (private) function that is used to calculate and 
// store the spline function on a regular grid.
//
/////////////////////////////////////////////////////////////////////

template<class T>
T* Spline1D<T>::spline_kernel(unsigned int order, unsigned int ksp, unsigned int deriv)
{
  if (!order) order = _order;
  if (!ksp) ksp = _ksp;
  if (!deriv) deriv = _deriv;
  unsigned int spsz = spline_size(order,ksp);
  T *tsp = new T[spsz];

  double coa = static_cast<double>(spsz-1) / 2.0; // Centre of array in voxels
  for (unsigned int i=0; i<spsz; i++) {
    double x = static_cast<double>(i-coa)/static_cast<double>(ksp);
    tsp[i] = kernel_value(x,order,ksp,deriv);
  }
  return(tsp);    
}

/////////////////////////////////////////////////////////////////////
//
// Obsolete Internal (private) function that used to calculate 
// and store the spline function on a regular grid.
//
/////////////////////////////////////////////////////////////////////

/*
template<class T>
T* Spline1D<T>::spline_kernel(unsigned int order, unsigned int ksp, unsigned int deriv)
{
  if (!order) order = _order;
  if (!ksp) ksp = _ksp;
  if (!deriv) deriv = _deriv;
  // cout << "order = " << order << ", ksp = " << ksp << ", spline_size(order,ksp) = " << spline_size(order,ksp) << endl;
  T  *tsp = new T[spline_size(order,ksp)];

  switch (order) {
  case 1:
    break;
  case 2:
    // cout << "I'm constructing 2nd order spline" << endl;
    {
      unsigned int coffs = spline_size(order,ksp) / 2;  // Deliberate truncation
      for (unsigned int i=0; i<coffs+1; i++) {
	double tmp = double(i) / double(ksp);
	switch (deriv) { 
	case 0:
	  if (tmp <= 0.5) tsp[coffs+i] = 0.75 - tmp*tmp; 
	  else if (tmp < 1.5) tsp[coffs+i] = 0.5 * (1.5-tmp) * (1.5-tmp);
	  if (i) tsp[coffs-i] = tsp[coffs+i];
	  break;
	case 1:
	  if (tmp <= 0.5) tsp[coffs+i] = - 2.0 * tmp / double(ksp);
	  else if (tmp < 1.5) tsp[coffs+i] = 0.5 * (2.0*tmp - 3) / double(ksp);
	  if (i) tsp[coffs-i] = -tsp[coffs+i];
	  break;
	case 2:
          if (fabs(tmp-1.5) < 1e-6) tsp[coffs+1] = 0.5 / double(ksp*ksp);    // If smack on "edge" discontinuity 
	  else if (fabs(tmp-0.5) < 1e-6) tsp[coffs+i] = -0.5 / double(ksp*ksp);   // If smack on "inner" discontinuity
	  else if (tmp < 0.5) tsp[coffs+i] = - 2.0 / double(ksp*ksp);
	  else if (tmp < 1.5) tsp[coffs+i] = 1.0 / double(ksp*ksp);
	  if (i) tsp[coffs-i] = tsp[coffs+i];
	  break;
	}
      }
    }
    break;
  case 3:
    // cout << "I'm constructing 3rd order spline" << endl;
    for (unsigned int i=1; i<ksp; i++) {
      double tmp = double(i) / double(ksp);
      switch (deriv) {
      case 0:
        tsp[i-1] = (1.0/6.0) * tmp*tmp*tmp;
        break;
      case 1:
        tsp[i-1] = (0.5/double(ksp)) * tmp*tmp;
        break;
      case 2:
        tsp[i-1] = (1.0/(double(ksp)*double(ksp))) * tmp;
        break;
      }
    }
    for (unsigned int i=ksp; i<(2*ksp+1); i++) {
      double tmp = double(i) / double(ksp);
      switch (deriv) {
      case 0:
        tsp[i-1] = (2.0/3.0) - 0.5*(2.0-tmp)*(2.0-tmp)*tmp;
        break;
      case 1:
        tsp[i-1] = (1.0/double(ksp)) *(-2.0*(tmp-2.0) - 1.5*(tmp-2.0)*(tmp-2.0));
        break;
      case 2:
        tsp[i-1] = (1.0/(double(ksp)*double(ksp))) * (4.0 - 3.0*tmp);
        break;
      }
    }
    for (unsigned int i=(2*ksp+1); i<(4*ksp); i++) {
      switch (deriv) {
      case 0: case 2:
        tsp[i-1] = tsp[spline_size(order,ksp)-i];
        break;
      case 1:
        tsp[i-1] = -tsp[spline_size(order,ksp)-i];
        break;
      }
    }
    break;
  default:
    break;
  }
  return(tsp);
}
*/

} // End namespace BASISFIELD

#endif // End #ifndef fsl_splines_h
