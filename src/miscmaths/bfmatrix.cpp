//
//    Definitions for class BFMatrix
//
//    Jesper Andersson, FMRIB Image Analysis Group
//
//    Copyright (C) 2007 University of Oxford 
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
#include <iostream>
#include <iomanip>
#include <fstream>
#include <boost/shared_ptr.hpp>
#include "newmat.h"
#include "newmatio.h"
#include "miscmaths.h"
#include "bfmatrix.h"

namespace MISCMATHS {

//
// Member functions for BFMatrix
//

NEWMAT::Matrix BFMatrix::SubMatrix(unsigned int fr, unsigned int lr, unsigned int fc, unsigned int lc) const
{
  if (fr<1 || fc<1 || lr>Nrows() || lc>Ncols() || fr>lr || fc>lc) throw BFMatrixException("BFMatrix::SubMatrix: index out of range");
  NEWMAT::Matrix omat(lr-fr+1,lc-fc+1);
  for (unsigned int r=fr, ri=1; r<=lr; r++, ri++) {
    for (unsigned int c=fc, ci=1; c<=lc; c++, ci++) {
      omat(ri,ci) = this->Peek(r,c);
    }
  }
  return(omat);
}

//
// Member functions for FullBFMatrix
//

void FullBFMatrix::Print(const std::string fname) const
{
  if (!fname.length()) cout << endl << *mp << endl;
  else write_ascii_matrix(fname,*mp);
}

boost::shared_ptr<BFMatrix> FullBFMatrix::Transpose() 
const
{
  boost::shared_ptr<FullBFMatrix>  tm(new FullBFMatrix(mp->t()));
  return(tm);
}

//
// Concatenate two matrices yielding a third
//

void FullBFMatrix::HorConcat(const BFMatrix& B, BFMatrix& AB) const
{
  if (B.Nrows() && Nrows() != B.Nrows()) {throw BFMatrixException("FullBFMatrix::HorConcat: Matrices must have same # of rows");}

  FullBFMatrix *pAB = dynamic_cast<FullBFMatrix *>(&AB);  
  if (pAB) { // This means output is a full matrix
    *pAB = *this;
    pAB->HorConcat2MyRight(B);
  }
  else {
    SparseBFMatrix<double> *psdAB = dynamic_cast<SparseBFMatrix<double> *>(&AB);
    if (psdAB) { 
      *psdAB = SparseBFMatrix<double>(this->AsMatrix());
      psdAB->HorConcat2MyRight(B);
    }
    else {
      SparseBFMatrix<float> *psfAB = dynamic_cast<SparseBFMatrix<float> *>(&AB);
      if (psfAB) { 
        *psfAB = SparseBFMatrix<float>(this->AsMatrix());
        psfAB->HorConcat2MyRight(B);
      }
      else throw BFMatrixException("FullBFMatrix::HorConcat: dynamic cast error"); 
    }
  }
}

void FullBFMatrix::HorConcat(const NEWMAT::Matrix& B, BFMatrix& AB) const
{
  if (B.Nrows() && int(Nrows()) != B.Nrows()) {throw BFMatrixException("FullBFMatrix::HorConcat: Matrices must have same # of rows");}

  FullBFMatrix *pAB = dynamic_cast<FullBFMatrix *>(&AB);  
  if (pAB) { // This means output is a full matrix
    *pAB = *this;
    pAB->HorConcat2MyRight(B);
  }
  else {
    SparseBFMatrix<double> *psdAB = dynamic_cast<SparseBFMatrix<double> *>(&AB);
    if (psdAB) {
      *psdAB = SparseBFMatrix<double>(this->AsMatrix());
      psdAB->HorConcat2MyRight(B);
    }
    else {
      SparseBFMatrix<float> *psfAB = dynamic_cast<SparseBFMatrix<float> *>(&AB);
      if (psfAB) {
        *psfAB = SparseBFMatrix<float>(this->AsMatrix());
        psfAB->HorConcat2MyRight(B);
      }
      else throw BFMatrixException("FullBFMatrix::HorConcat: dynamic cast error"); 
    }
  }
}

void FullBFMatrix::VertConcat(const BFMatrix& B, BFMatrix& AB) const
{
  if (B.Ncols() && Ncols() != B.Ncols()) {throw BFMatrixException("FullBFMatrix::VertConcat: Matrices must have same # of columns");}

  FullBFMatrix *pAB = dynamic_cast<FullBFMatrix *>(&AB);  
  if (pAB) { // This means output is a full matrix
    *pAB = *this;
    pAB->VertConcatBelowMe(B);
  }
  else {
    SparseBFMatrix<double> *psdAB = dynamic_cast<SparseBFMatrix<double> *>(&AB);
    if (psdAB) {
      *psdAB = SparseBFMatrix<double>(this->AsMatrix());
      psdAB->VertConcatBelowMe(B);
    }
    else {
      SparseBFMatrix<float> *psfAB = dynamic_cast<SparseBFMatrix<float> *>(&AB);
      if (psfAB) {
        *psfAB = SparseBFMatrix<float>(this->AsMatrix());
        psfAB->VertConcatBelowMe(B);
      }
      else throw BFMatrixException("FullBFMatrix::VertConcat: dynamic cast error"); 
    }
  }
}

void FullBFMatrix::VertConcat(const NEWMAT::Matrix& B, BFMatrix& AB) const
{
  if (B.Ncols() && int(Ncols()) != B.Ncols()) {throw BFMatrixException("FullBFMatrix::VertConcat: Matrices must have same # of columns");}

  FullBFMatrix *pAB = dynamic_cast<FullBFMatrix *>(&AB);  
  if (pAB) { // This means output is a full matrix
    *pAB = *this;
    pAB->VertConcatBelowMe(B);
  }
  else {
    SparseBFMatrix<double> *psdAB = dynamic_cast<SparseBFMatrix<double> *>(&AB);
    if (psdAB) {
      *psdAB = SparseBFMatrix<double>(this->AsMatrix());
      psdAB->VertConcatBelowMe(B);
    }
    else {
      SparseBFMatrix<float> *psfAB = dynamic_cast<SparseBFMatrix<float> *>(&AB);
      if (psfAB) {
        *psfAB = SparseBFMatrix<float>(this->AsMatrix());
        psfAB->VertConcatBelowMe(B);
      }
      else throw BFMatrixException("FullBFMatrix::VertConcat: dynamic cast error"); 
    }
  }
}

//  
// Concatenation of another matrix to *this
//
void FullBFMatrix::HorConcat2MyRight(const BFMatrix& B)
{
  if (!B.Nrows()) return;

  if (Nrows() != B.Nrows()) {throw BFMatrixException("FullBFMatrix::HorConcat2MyRight: Matrices must have same # of rows");}

  const FullBFMatrix *pB = dynamic_cast<const FullBFMatrix *>(&B);
  if (pB) { // If B was full
    *mp |= *(pB->mp);
  }
  else {
    const SparseBFMatrix<double> *psdB = dynamic_cast<const SparseBFMatrix<double> *>(&B);
    if (psdB) {
      this->HorConcat2MyRight(psdB->AsMatrix());
    }
    else { 
      const SparseBFMatrix<float> *psfB = dynamic_cast<const SparseBFMatrix<float> *>(&B);
      if (psfB) {
        this->HorConcat2MyRight(psfB->AsMatrix());
      }
      else throw BFMatrixException("FullBFMatrix::HorConcat2MyRight: dynamic cast error"); 
    }
  }
}

void FullBFMatrix::HorConcat2MyRight(const NEWMAT::Matrix& B)
{
  if (!B.Nrows()) return;

  if (int(Nrows()) != B.Nrows()) {throw BFMatrixException("FullBFMatrix::HorConcat2MyRight: Matrices must have same # of rows");}
  *mp |= B;
}

void FullBFMatrix::VertConcatBelowMe(const BFMatrix& B)
{
  if (!B.Ncols()) return;

  if (Ncols() != B.Ncols()) {throw BFMatrixException("FullBFMatrix::VertConcatBelowMe: Matrices must have same # of columns");}

  const FullBFMatrix *pB = dynamic_cast<const FullBFMatrix *>(&B);
  if (pB) { // Means B is full
    *mp &= *(pB->mp);
  }
  else {
    const SparseBFMatrix<double> *psdB = dynamic_cast<const SparseBFMatrix<double> *>(&B);
    if (psdB) {
      this->VertConcatBelowMe(psdB->AsMatrix());
    }
    else {
      const SparseBFMatrix<float> *psfB = dynamic_cast<const SparseBFMatrix<float> *>(&B);
      if (psfB) {
        this->VertConcatBelowMe(psfB->AsMatrix());
      }
      else throw BFMatrixException("FullBFMatrix::HorConcatBelowMe: dynamic cast error"); 
    }
  }
}

void FullBFMatrix::VertConcatBelowMe(const NEWMAT::Matrix& B)
{
  if (!B.Ncols()) return;

  if (int(Ncols()) != B.Ncols()) {throw BFMatrixException("FullBFMatrix::VertConcatBelowMe: Matrices must have same # of columns");}
  *mp &= B;
}

// Multiply this matrix with scalar

void FullBFMatrix::MulMeByScalar(double s)
{
  *mp = s*(*mp);
}

// Multiply by vector
NEWMAT::ReturnMatrix FullBFMatrix::MulByVec(const NEWMAT::ColumnVector& invec) const
{
  if (invec.Nrows() != int(Ncols())) {throw BFMatrixException("FullBFMatrix::MulByVec: Matrix-vector size mismatch");}
  NEWMAT::ColumnVector  ret;
  ret = (*mp)*invec;
  ret.Release();
  return(ret);
}

// Add another matrix to this one
void FullBFMatrix::AddToMe(const BFMatrix& m, double s)
{
  if (Ncols() != m.Ncols() || Nrows() != m.Nrows()) {
    throw BFMatrixException("FullBFMatrix::AddToMe: Matrix size mismatch");
  }

  const FullBFMatrix *pm = dynamic_cast<const FullBFMatrix *>(&m);
  if (pm) { // If m is full matrix
    *mp += s*(*(pm->mp));
  }
  else {
    const SparseBFMatrix<double> *psdm = dynamic_cast<const SparseBFMatrix<double> *>(&m);
    if (psdm) *mp += s*psdm->AsMatrix();
    else {
      const SparseBFMatrix<float> *psfm = dynamic_cast<const SparseBFMatrix<float> *>(&m);
      if (psfm) *mp += s*psfm->AsMatrix();
      else throw BFMatrixException("FullBFMatrix::AddToMe: dynamic cast error");
    }
  }
}

// Given A*x=b, solve for x
NEWMAT::ReturnMatrix FullBFMatrix::SolveForx(const NEWMAT::ColumnVector& b,       // Ignoring all parameters except b
					     MISCMATHS::MatrixType       type,
					     double                      tol,
                                             int                         miter) const
{
  if (int(Nrows()) != b.Nrows()) {throw BFMatrixException("FullBFMatrix::SolveForx: Matrix-vector size mismatch");}
  NEWMAT::ColumnVector  ret;
  ret = mp->i()*b;
  ret.Release();
  return(ret);
}

BFMatrixColumnIterator BFMatrix::begin(unsigned int col) const 
{
  if (col > Ncols()) throw BFMatrixException("BFMatrix:begin col out of range"); 
  return(BFMatrixColumnIterator(*this,col));
}

BFMatrixColumnIterator BFMatrix::end(unsigned int col) const 
{
  if (col > Ncols()) throw BFMatrixException("BFMatrix:begin col out of range"); 
  return(BFMatrixColumnIterator(*this,col,true));
}
    
} // End namespace MISCMATHS
