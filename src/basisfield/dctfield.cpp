// Definitions for class dctfield
//
// dctfield.cpp
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
#include <iomanip>
#include <fstream>
#include <vector>
#include <boost/shared_ptr.hpp>
#include "newmat.h"
#include "miscmaths/bfmatrix.h"
#include "dctfield.h"

#ifndef SQR
#define SQR(A) ((A)*(A))
#endif

using namespace std;
using namespace NEWMAT;

namespace BASISFIELD {

// Constructor, assignement and destructor

dctfield::dctfield(const std::vector<unsigned int>& psz, const std::vector<double>& pvxs, const std::vector<unsigned int>& porder)
  : basisfield(psz,pvxs), dctbas(3)
{
  if (porder.size() != NDim()) {throw BasisfieldException("dctfield::dctfield: Dimensionality mismatch");}
  if (porder[0]<1 || porder[0]>FieldSz_x() || (NDim()>1 && (porder[1]>FieldSz_y())) || (NDim()>2 && (porder[2]>FieldSz_z()))) {
    throw BasisfieldException("dctfield::dctfield: Invalid order of DCT transform");
  }

  unsigned int lsz[] = {FieldSz_x(), FieldSz_y(), FieldSz_z()};
  for (unsigned int i=0; i<3; i++) {
    dctbas[i] = std::vector<boost::shared_ptr<NEWMAT::Matrix> >(3);
    int lorder = (i<NDim()) ? porder[i] : 1;
    for (int d=0; d<3; d++) {
      dctbas[i][d] = make_dctbas(lsz[i],lorder,d);
    }
  }

  boost::shared_ptr<NEWMAT::ColumnVector>  lcoef(new NEWMAT::ColumnVector(CoefSz()));
  *lcoef = 0.0;
  set_coef_ptr(lcoef);  
}

dctfield::dctfield(const dctfield& inf) 
: basisfield(inf), dctbas(3)
{
  basisfield::assign_basisfield(inf);
  dctfield::assign_dctfield(inf);
}

void dctfield::assign_dctfield(const dctfield& inf)
{
  for (int i=0; i<3; i++) {
    dctbas[i] = std::vector<boost::shared_ptr<NEWMAT::Matrix> >(3);
    for (int j=0; j<3; j++) {
      dctbas[i][j] = boost::shared_ptr<NEWMAT::Matrix>(new NEWMAT::Matrix(*(inf.dctbas[i][j])));
    }
  }
}

dctfield& dctfield::operator=(const dctfield& inf)
{
  if (&inf == this) {return(*this);} // Detect self

  dctbas = std::vector<std::vector<boost::shared_ptr<NEWMAT::Matrix> > >(3);
  basisfield::assign_basisfield(inf);   // Assign common part
  dctfield::assign_dctfield(inf);     // Assign dctfield specific bits

  return(*this);
}

// General utility functions

// Functions that actually do some work

void dctfield::SetToConstant(double fv)
{
  double cfac = dctbas[0][0]->element(0,0) * dctbas[1][0]->element(0,0) * dctbas[2][0]->element(0,0);
  double coefval = fv / cfac;
  NEWMAT::ColumnVector  lcoef(CoefSz());
  lcoef = 0.0;
  lcoef(1) = coefval;
  SetCoef(lcoef);
}

std::vector<double> dctfield::SubsampledVoxelSize(const std::vector<unsigned int>&  subsampling,
					          std::vector<double>               ovxs,
					          std::vector<unsigned int>         isize) const
{
  std::vector<unsigned int>     osize(3), nsize(3);
  if (isize.size()) nsize = SubsampledMatrixSize(subsampling,osize=isize);
  else { 
    osize[0] = FieldSz_x(); osize[1] = FieldSz_y(); osize[2] = FieldSz_z();
    nsize = SubsampledMatrixSize(subsampling,osize);
  }
    
  std::vector<double>  nvxs(ovxs);
  if (!nvxs.size()) {nvxs.resize(3); nvxs[0]=Vxs_x(); nvxs[1]=Vxs_y(); nvxs[2]=Vxs_z();}
  for (int i=0; i<int(nvxs.size()); i++) {
    nvxs[i] = nvxs[i]*double(osize[i]+1)/double(nsize[i]+1);
  }
  return(nvxs); 
}

std::vector<unsigned int> dctfield::SubsampledMatrixSize(const std::vector<unsigned int>&   subsampling,
                                                         std::vector<unsigned int>          isize) const
{
  std::vector<unsigned int>  osize(3,0);
  if (isize.size()) osize = isize;
  else {osize[0] = FieldSz_x(); osize[1] = FieldSz_y(); osize[2] = FieldSz_z();}
  for (unsigned int i=0; i<osize.size(); i++) {
    while (osize[i]%subsampling[i]) osize[i] += 1;
    osize[i] /= subsampling[i];
  }
  return(osize);
}

void dctfield::Update(FieldIndex fi)
{
  if (fi>int(NDim())) {throw BasisfieldException("dctfield::Update: Cannot take derivative in singleton direction");}

  if (UpToDate(fi)) {return;} // Field already fine.

  const boost::shared_ptr<NEWMAT::ColumnVector> lcoef = GetCoef();
  if (!lcoef) {throw BasisfieldException("dctfield::Update: No coefficients set yet");}

  boost::shared_ptr<NEWMAT::ColumnVector> fptr = get_ptr(fi);
  int tmp[] = {0,0,0};
  if (fi) {tmp[fi-1] = 1;}
  AkBkCxb(*(dctbas[0][tmp[0]]),*(dctbas[1][tmp[1]]),*(dctbas[2][tmp[2]]),*lcoef,*fptr);
  set_update_flag(true,fi);
}

//
// All the Jte and JtJ routines are presently just stubs. However, all the 
// functionality has been properly implemented. Have a look at the 
// deprecated ColumnVector versions of these routines in dead_code.cpp
// to see how to get these routines functional again.
//

NEWMAT::ReturnMatrix dctfield::Jte(const NEWIMAGE::volume<float>&  ima1,
                                   const NEWIMAGE::volume<float>&  ima2,
                                   const NEWIMAGE::volume<char>    *mask) const
{
  NEWMAT::ColumnVector tmp(1);
  tmp.Release();
  return(tmp);
}

NEWMAT::ReturnMatrix dctfield::Jte(const std::vector<unsigned int>&  deriv,
                                   const NEWIMAGE::volume<float>&    ima1,
                                   const NEWIMAGE::volume<float>&    ima2,
                                   const NEWIMAGE::volume<char>      *mask) const
{
  NEWMAT::ColumnVector tmp(1);
  tmp.Release();
  return(tmp);
}

NEWMAT::ReturnMatrix dctfield::Jte(const NEWIMAGE::volume<float>&    ima,
                                   const NEWIMAGE::volume<char>      *mask) const
{
  NEWMAT::ColumnVector tmp(1);
  tmp.Release();
  return(tmp);
}


NEWMAT::ReturnMatrix dctfield::Jte(const std::vector<unsigned int>&  deriv,
                                   const NEWIMAGE::volume<float>&    ima,
                                   const NEWIMAGE::volume<char>      *mask) const
{
  NEWMAT::ColumnVector tmp(1);
  tmp.Release();
  return(tmp);
}

boost::shared_ptr<MISCMATHS::BFMatrix> dctfield::JtJ(const NEWIMAGE::volume<float>&    ima,
                                                     const NEWIMAGE::volume<char>      *mask,
                                                     MISCMATHS::BFMatrixPrecisionType  prec) const
{
  boost::shared_ptr<BFMatrix> H = boost::shared_ptr<BFMatrix>(new FullBFMatrix);
  return(H);
}

boost::shared_ptr<MISCMATHS::BFMatrix> dctfield::JtJ(const NEWIMAGE::volume<float>&       ima1,
                                                     const NEWIMAGE::volume<float>&       ima2,
                                                     const NEWIMAGE::volume<char>         *mask,
                                                     MISCMATHS::BFMatrixPrecisionType     prec) const
{
  boost::shared_ptr<BFMatrix> H = boost::shared_ptr<BFMatrix>(new FullBFMatrix);
  return(H);
}

boost::shared_ptr<MISCMATHS::BFMatrix> dctfield::JtJ(const std::vector<unsigned int>&       deriv, 
                                                     const NEWIMAGE::volume<float>&         ima,
                                                     const NEWIMAGE::volume<char>           *mask,
                                                     MISCMATHS::BFMatrixPrecisionType       prec) const
{
  boost::shared_ptr<BFMatrix> H = boost::shared_ptr<BFMatrix>(new FullBFMatrix);
  return(H);
}

boost::shared_ptr<MISCMATHS::BFMatrix> dctfield::JtJ(const std::vector<unsigned int>&       deriv, 
                                                     const NEWIMAGE::volume<float>&         ima1,
                                                     const NEWIMAGE::volume<float>&         ima2,
                                                     const NEWIMAGE::volume<char>           *mask,
                                                     MISCMATHS::BFMatrixPrecisionType       prec) const
{
  boost::shared_ptr<BFMatrix> H = boost::shared_ptr<BFMatrix>(new FullBFMatrix);
  return(H);
}

boost::shared_ptr<MISCMATHS::BFMatrix> dctfield::JtJ(const std::vector<unsigned int>&         deriv1,
                                                     const NEWIMAGE::volume<float>&           ima1,
                                                     const std::vector<unsigned int>&         deriv2,
                                                     const NEWIMAGE::volume<float>&           ima2,
                                                     const NEWIMAGE::volume<char>             *mask,
                                                     MISCMATHS::BFMatrixPrecisionType         prec) const
{
  boost::shared_ptr<BFMatrix> H = boost::shared_ptr<BFMatrix>(new FullBFMatrix);
  return(H);
}

boost::shared_ptr<MISCMATHS::BFMatrix> dctfield::JtJ(const NEWIMAGE::volume<float>&        ima1,
                                                     const basisfield&                     bf2,
                                                     const NEWIMAGE::volume<float>&        ima2,
                                                     const NEWIMAGE::volume<char>          *mask,
                                                     MISCMATHS::BFMatrixPrecisionType      prec) const
{
  boost::shared_ptr<BFMatrix> H = boost::shared_ptr<BFMatrix>(new FullBFMatrix);
  return(H);
}


double dctfield::MemEnergy() const // Membrane energy of the field
{
  const boost::shared_ptr<NEWMAT::ColumnVector> lcoef = GetCoef();
  if (!lcoef) {throw BasisfieldException("dctfield::MemEnergy: No coefficients set yet");}

  NEWMAT::DiagonalMatrix H(CoefSz());
  memen_H(H);
  return(static_cast<double>(DotProduct(*(lcoef),H*(*lcoef))));
}

NEWMAT::ReturnMatrix dctfield::MemEnergyGrad() const // Gradient of the membrane energy
{
  const boost::shared_ptr<NEWMAT::ColumnVector> lcoef = GetCoef();
  if (!lcoef) {throw BasisfieldException("dctfield::MemEnergy: No coefficients set yet");}

  ColumnVector  grad(CoefSz());
  NEWMAT::DiagonalMatrix H(CoefSz());
  memen_H(H);
  grad = H*(*lcoef);
  grad.Release();
  return(grad);
}
 
boost::shared_ptr<BFMatrix> dctfield::MemEnergyHess(MISCMATHS::BFMatrixPrecisionType prec) const // Hessian of membrane energy
{
  if (prec!=BFMatrixDoublePrecision) throw BasisfieldException("dctfield::MemEnergyHess: Hessian must be double precision"); 

  boost::shared_ptr<FullBFMatrix>     H(new FullBFMatrix);
  NEWMAT::DiagonalMatrix              tmp(CoefSz());
  memen_H(tmp);
  boost::shared_ptr<NEWMAT::Matrix>   mp(new Matrix(tmp));   // Wasteful, but can't make it work any other way. 
  H->SetMatrixPtr(mp);

  return(H);
}

boost::shared_ptr<BASISFIELD::basisfield> dctfield::ZoomField(const std::vector<unsigned int>&    psz,
                                                              const std::vector<double>&          pvxs,
                                                              std::vector<unsigned int>           porder) const
{
  // Make sure requested matrix-size is valid
  for (int i=0; i<3; i++) {
    if (!valid_size(psz[i],pvxs[i],porder[i],i)) throw BasisfieldException("dctfield::ZoomField: Invalid zoom requsted");
  }  
 
  boost::shared_ptr<BASISFIELD::dctfield> tptr(new BASISFIELD::dctfield(psz,pvxs,porder));

  const boost::shared_ptr<NEWMAT::ColumnVector> ocptr = GetCoef();
  boost::shared_ptr<NEWMAT::ColumnVector> ncptr = tptr->get_coef();

  // Repack old coefficients into relevant part of new coefficient vector

  int oorder[3] = {CoefSz_x(), CoefSz_y(), CoefSz_z()};
  *ncptr = 0.0;
  for (int z=0; z<oorder[2]; z++) {
    for (int y=0; y<oorder[1]; y++) {
      for (int x=0; x<oorder[0]; x++) {
        ncptr->element(z*porder[0]*porder[1]+y*porder[0]+x) = ocptr->element(z*oorder[0]*oorder[1]+y*oorder[0]+x);
      }
    }
  }

  return(tptr);
}

// Functions that are declared private

boost::shared_ptr<NEWMAT::Matrix> dctfield::make_dctbas(unsigned int size, unsigned int order, unsigned int deriv) const
{
  boost::shared_ptr<NEWMAT::Matrix>  tmptr(new NEWMAT::Matrix(size,order));
  for (unsigned int i=0; i<size; i++) {
    tmptr->element(i,0) = (deriv) ? 0 : 1.0/sqrt(double(size));
  }
  for (unsigned int i=0; i<size; i++) {
    for (unsigned int j=1; j<order; j++) {
      switch (deriv) {
      case 0:
        tmptr->element(i,j) = sqrt(2.0/size) * cos((PI*(2.0*i+1.0)*j)/(2.0*size));
	break;
      case 1:
	tmptr->element(i,j) = -sqrt(2.0/size) * sin((PI*(2.0*i+1.0)*j)/(2.0*size)) * PI*j/double(size);
        break;
      case 2:
	tmptr->element(i,j) = -sqrt(2.0/size) * cos((PI*(2.0*i+1.0)*j)/(2.0*size)) * SQR((PI*j)/double(size));
	break;
      }
    }
  }
  return(tmptr);
}

bool dctfield::valid_size(unsigned int psz, double pvxs, unsigned int porder, unsigned int pdim) const
{
  // Make sure we are upsampling or leaving the same
  unsigned int osz[3] = {FieldSz_x(), FieldSz_y(), FieldSz_z()};
  if (psz < osz[pdim]) return(false);
  // Make sure we have same or higher order
  unsigned int oorder[3] = {CoefSz_x(), CoefSz_y(), CoefSz_z()};
  if (porder < oorder[pdim]) return(false);
  // Make sure upsampling step is valid
  std::vector<unsigned int>   nsz(1,0);
  bool valid = false;
  unsigned int ss = 1;
  for (; ss<64; ss*=2) {
    nsz[0] = psz;
    nsz = SubsampledMatrixSize(ss,nsz);
    if (nsz[0] == osz[pdim]) {
      valid = true;
      break;
    }
  }
  // Make sure voxel-size matches
  if (valid) {
    double eps = 1e-6;
    double ovxs[] = {Vxs_x(), Vxs_y(), Vxs_z()};
    std::vector<double> nvxs(1,pvxs);
    nsz[0] = psz;
    nvxs = SubsampledVoxelSize(ss,nvxs,nsz);
    valid = (fabs(nvxs[0]-ovxs[pdim]) < eps); 
  }
  return(valid);
}

//
// Calculates (AkBkC)*b where A, B and C are matrices, k denotes
// Kronecker tensor product and b is a colum vector with 
// A.Ncols()*B.Ncols()*C.Ncols() rows.
//
// Useful for calulating the field, and also for calculating the
// gradient.
// 
void dctfield::AkBkCxb(const NEWMAT::Matrix&                   A, 
                       const NEWMAT::Matrix&                   B, 
                       const NEWMAT::Matrix&                   C, 
                       const NEWMAT::ColumnVector&             b, 
                       NEWMAT::ColumnVector&                   ret) const
{
  ret = 0.0;
  for (int mk=0; mk<C.Ncols(); mk++) {
    NEWMAT::Matrix tmpM(A.Nrows(),B.Ncols());
    for (int mj=0; mj<B.Ncols(); mj++) {
      tmpM.Column(mj+1) = A * b.Rows(mk*A.Ncols()*B.Ncols()+mj*A.Ncols()+1,mk*A.Ncols()*B.Ncols()+(mj+1)*A.Ncols());
    }
    NEWMAT::ColumnVector tmpV(A.Nrows()*B.Nrows());
    for (int nj=0; nj<B.Nrows(); nj++) {
      tmpV.Rows(nj*A.Nrows()+1,(nj+1)*A.Nrows()) = tmpM*(B.Row(nj+1)).t();
    }
    for (int nk=0; nk<C.Nrows(); nk++) {
      ret.Rows(nk*A.Nrows()*B.Nrows()+1,(nk+1)*A.Nrows()*B.Nrows()) += C.element(nk,mk) * tmpV;
    }
  }    
}

//
// The following routines could be written much more clearly if I 
// used the capabilities of the NEWMAT classes. But it would also
// mean that a LOT of constructors would be invoked, and a lot of
// bounds would be checked. So I have chosen to (partly) revert to 
// good ole C programming.
//

void dctfield::AtA(const NEWMAT::Matrix&        Bx,
                   const NEWMAT::Matrix&        By,
                   const NEWMAT::Matrix&        Bz,
                   const NEWMAT::ColumnVector&  ima,
                   NEWMAT::Matrix&              AtA) const
{
  double  *istore = ima.Store();
  double  *AtAp = AtA.Store();
  int     znc = Bz.Ncols();
  int     xysz = Bx.Ncols()*By.Ncols();
  int     AtAsz = Bx.Ncols()*By.Ncols()*Bz.Ncols();
  double  *AtA1p = new double[SQR(Bx.Ncols()*By.Ncols())];

  AtA = 0.0;
  for (int z=0; z<Bz.Nrows(); z++) {
    // Get contribution for this slice
    double *iptr = &(istore[z*Bx.Nrows()*By.Nrows()]);
    one_slice_AtA(Bx,By,iptr,AtA1p);
    // Tile scaled versions of AtA1p
    for (int zr=0; zr<znc; zr++) {
      for (int zc=zr; zc<znc; zc++) {
        double scl = Bz.element(z,zr)*Bz.element(z,zc);
        for (int xyr=0; xyr<xysz; xyr++) {
          for (int xyc=xyr; xyc<xysz; xyc++) {
	    AtAp[(zr*xysz+xyr)*AtAsz + zc*xysz+xyc] += scl*AtA1p[xyr*xysz+xyc];
	  }
	}
      }
    }
  }
  //
  // We have now filled in the upper diagonal part of the
  // upper diagonal set of tiles. Now we need to do the rest.
  //
  // Fill in lower diagonals of upper diagonal tiles
  //
  for (int zr=0; zr<znc; zr++) {
    for (int zc=zr; zc<znc; zc++) {
      for (int xyr=0; xyr<xysz; xyr++) {
	for (int xyc=0; xyc<xyr; xyc++) {
	  AtAp[(zr*xysz+xyr)*AtAsz + zc*xysz+xyc] = AtAp[(zr*xysz+xyc)*AtAsz + zc*xysz+xyr];
	}
      }
    }
  }
  // Fill in lower diagonal tiles
  for (int zr=0; zr<znc; zr++) {
    for (int zc=0; zc<zr; zc++) {
      for (int xyr=0; xyr<xysz; xyr++) {
	for (int xyc=0; xyc<xysz; xyc++) {
	  AtAp[(zr*xysz+xyr)*AtAsz + zc*xysz+xyc] = AtAp[(zc*xysz+xyr)*AtAsz + zr*xysz+xyc];
	}
      }
    }
  }  
  
  delete[] AtA1p;
}

void dctfield::one_slice_AtA(const NEWMAT::Matrix&   Bx,
                             const NEWMAT::Matrix&   By,
                             const double            *ima,
                             double                  *AtA) const
{
  NEWMAT::Matrix           XtX(Bx.Ncols(),Bx.Ncols());
  NEWMAT::Matrix           Bxt = Bx.t();
  NEWMAT::DiagonalMatrix   pima(Bx.Nrows());
  int                      xnc = Bx.Ncols();
  int                      ync = By.Ncols();
  int                      AtAsz = xnc*ync;

  memset(AtA,0,SQR(AtAsz)*sizeof(double));
  for (int my=0; my<By.Nrows(); my++) {
    const double *iptr = &(ima[my*Bx.Nrows()]);
    pima << iptr;
    XtX = Bxt*pima*Bx;     // This line can be made more efficient
    double *XtXp = XtX.Store();
    // Now we tile scaled versions of XtX
    for (int yr=0; yr<ync; yr++) {
      for (int yc=yr; yc<ync; yc++) {
        double scl = By.element(my,yr)*By.element(my,yc);
	for (int xr=0; xr<xnc; xr++) {
	  for (int xc=xr; xc<xnc; xc++) {
            AtA[(yr*xnc+xr)*AtAsz + yc*xnc+xc] += scl*XtXp[xr*xnc+xc];
	  }
	}
      }
    }
  }
  //
  // We have now filled in the upper diagonal part of the
  // upper diagonal set of tiles. Now we need to do the rest.
  //
  // Fill in lower diagonals of upper diagonal tiles
  //
  for (int yr=0; yr<ync; yr++) {
    for (int yc=yr; yc<ync; yc++) {
      for (int xr=0; xr<xnc; xr++) {
	for (int xc=0; xc<xr; xc++) {
	  AtA[(yr*xnc+xr)*AtAsz + yc*xnc+xc] = AtA[(yr*xnc+xc)*AtAsz + yc*xnc+xr];
	}
      }
    }
  }
  //
  // Fill in lower diagonal tiles
  //
  for (int yr=0; yr<ync; yr++) {
    for (int yc=0; yc<yr; yc++) {
      for (int xr=0; xr<xnc; xr++) {
	for (int xc=0; xc<xnc; xc++) {
	  AtA[(yr*xnc+xr)*AtAsz + yc*xnc+xc] = AtA[(yc*xnc+xc)*AtAsz + yr*xnc+xr];
	}
      }
    }
  }
  /*
  NEWMAT::Matrix  AA(AtAsz,AtAsz);
  AA << AtA;
  cout << std::setprecision(20) << AA << endl;      
  */
}

void dctfield::AtB(const NEWMAT::Matrix&        Ax,
		   const NEWMAT::Matrix&        Ay,
		   const NEWMAT::Matrix&        Az,
                   const NEWMAT::Matrix&        Bx,
                   const NEWMAT::Matrix&        By,
                   const NEWMAT::Matrix&        Bz,
                   const NEWMAT::ColumnVector&  ima,
                   NEWMAT::Matrix&              AtB) const
{
  const double  *istore = ima.Store();
  double        *AtBp = AtB.Store();
  int           znc = Az.Ncols();
  int           xysz = Ax.Ncols()*Ay.Ncols();
  int           AtBsz = Ax.Ncols()*Ay.Ncols()*Az.Ncols();
  double        *AtB1p = new double[SQR(Ax.Ncols()*Ay.Ncols())];

  for (int z=0; z<Az.Nrows(); z++) {
    // Get contribution for this slice
    const double *iptr = &(istore[z*Ax.Nrows()*Ay.Nrows()]);
    one_slice_AtB(Ax,Ay,Bx,By,iptr,AtB1p);
    // Tile scaled versions of AtB1p
    for (int zr=0; zr<znc; zr++) {
      for (int zc=0; zc<znc; zc++) {
        double scl = Az.element(z,zr)*Bz.element(z,zc);
        for (int xyr=0; xyr<xysz; xyr++) {
          for (int xyc=0; xyc<xysz; xyc++) {
	    AtBp[(zr*xysz+xyr)*AtBsz + zc*xysz+xyc] += scl*AtB1p[xyr*xysz+xyc];
	  }
	}
      }
    }
  }

  delete[] AtB1p;
}

void dctfield::one_slice_AtB(const NEWMAT::Matrix&   Ax,
                             const NEWMAT::Matrix&   Ay,
                             const NEWMAT::Matrix&   Bx,
                             const NEWMAT::Matrix&   By,
                             const double            *ima,
                             double                  *AtB) const
{
  NEWMAT::Matrix           XtX(Ax.Ncols(),Bx.Ncols());
  NEWMAT::Matrix           Axt = Ax.t();
  NEWMAT::DiagonalMatrix   pima(Ax.Nrows());
  int                      xnc = Ax.Ncols();
  int                      ync = Ay.Ncols();
  int                      AtBsz = xnc*ync;

  memset(AtB,0,SQR(AtBsz)*sizeof(double));
  for (int my=0; my<By.Nrows(); my++) {
    const double *iptr = &(ima[my*Ax.Nrows()]);
    pima << iptr;
    XtX = Axt*pima*Bx;     
    double *XtXp = XtX.Store();
    // Now we tile scaled versions of XtX
    for (int yr=0; yr<ync; yr++) {
      for (int yc=0; yc<ync; yc++) {
        double scl = Ay.element(my,yr)*By.element(my,yc);
	for (int xr=0; xr<xnc; xr++) {
	  for (int xc=0; xc<xnc; xc++) {
            AtB[(yr*xnc+xr)*AtBsz + yc*xnc+xc] += scl*XtXp[xr*xnc+xc];
	  }
	}
      }
    }
  }
}

void dctfield::memen_H(NEWMAT::DiagonalMatrix&  mH) const
{
  
  mH = 0.0;
  for (int d=0; d<3; d++) {
    std::vector<int>  dv(3,0);
    dv[d] = 1;
    std::vector<DiagonalMatrix> AtAv(3);
    int csz[] = {CoefSz_x(), CoefSz_y(), CoefSz_z()};
    for (int dim=0; dim<3; dim++) {
      AtAv[dim] = NEWMAT::DiagonalMatrix(CoefSz_x());	
      for (int i=1; i<=csz[dim]; i++) {
        AtAv[dim](i) = DotProduct(dctbas[dim][dv[0]]->Column(i),dctbas[dim][dv[0]]->Column(i));
      }
    }
    mH += KP(AtAv[2],KP(AtAv[1],AtAv[0]));
  }      
}

} // End namespace BASISFIELD
