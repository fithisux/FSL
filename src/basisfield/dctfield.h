// Declarations of class dctfield
//
// dctfield.h
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

#ifndef dctfield_h
#define dctfield_h

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
#include "newmat.h"
#include "miscmaths/bfmatrix.h"
#include "splines.h"
#include "basisfield.h"

namespace BASISFIELD {

static const double PI = 3.141592653589793;

class dctfield: public basisfield 
{
private:
  // Think of dctbas as a vector of vectors of matrices that define
  // the DCT basis set. The first index indicate direction (i.e. x-
  // y- or z-direction). The second index indicates derivative, i.e.
  // 0th, 1st or 2nd order derivative. E.g. *(dctbas[1][0]) is of
  // type Matrix and is the 0th order derivative of the DCT basis
  // set for the y-direction.
  std::vector<std::vector<boost::shared_ptr<NEWMAT::Matrix> > >     dctbas;    // DCT basis in the x-, y- and z-directions

  // Functions for internal use

  boost::shared_ptr<NEWMAT::Matrix>  make_dctbas(unsigned int size, unsigned int order, unsigned int deriv=0) const;
  void AkBkCxb(const NEWMAT::Matrix&                   A, 
               const NEWMAT::Matrix&                   B, 
               const NEWMAT::Matrix&                   C, 
               const NEWMAT::ColumnVector&             b, 
               NEWMAT::ColumnVector&                   ret) const;
  void AtA(const NEWMAT::Matrix&        Bx,
           const NEWMAT::Matrix&        By,
           const NEWMAT::Matrix&        Bz,
           const NEWMAT::ColumnVector&  ima,
           NEWMAT::Matrix&              AtA) const;
  void one_slice_AtA(const NEWMAT::Matrix&   Bx,
                     const NEWMAT::Matrix&   By,
                     const double            *ima,
                     double                  *AtA) const;
  void AtB(const NEWMAT::Matrix&        Ax,
	   const NEWMAT::Matrix&        Ay,
	   const NEWMAT::Matrix&        Az,
           const NEWMAT::Matrix&        Bx,
           const NEWMAT::Matrix&        By,
           const NEWMAT::Matrix&        Bz,
           const NEWMAT::ColumnVector&  ima,
           NEWMAT::Matrix&              AtB) const;
  void one_slice_AtB(const NEWMAT::Matrix&   Ax,
                     const NEWMAT::Matrix&   Ay,
                     const NEWMAT::Matrix&   Bx,
                     const NEWMAT::Matrix&   By,
                     const double            *ima,
                     double                  *AtB) const;
  void memen_H(NEWMAT::DiagonalMatrix&  mH) const;



protected:

  // Functions for use in this and derived classes
  virtual void assign_dctfield(const dctfield& inf);

public:

  // Constructors and destructors, including assignment

  dctfield(const std::vector<unsigned int>& psz, const std::vector<double>& pvxs, const std::vector<unsigned int>& order);
  dctfield(const dctfield& inf);
  dctfield& operator=(const dctfield& inf);
  virtual ~dctfield() {} // This should drop straight through to base-class

  bool valid_size(unsigned int psz, double pvxs, unsigned int porder, unsigned int pdim) const;

  // Explicit instruction to compiler that we intend not to refine some Peek functions
  using basisfield::Peek;
  // Getting the value for a non-integer voxel location
  virtual double Peek(double x, double y, double z, FieldIndex fi=FIELD) const {return(0.0);}

  // General utility functions

  virtual unsigned int CoefSz_x() const {return(static_cast<unsigned int>(dctbas[0][0]->Ncols()));}
  virtual unsigned int CoefSz_y() const {return(static_cast<unsigned int>(dctbas[1][0]->Ncols()));}
  virtual unsigned int CoefSz_z() const {return(static_cast<unsigned int>(dctbas[2][0]->Ncols()));}

  virtual bool HasGlobalSupport() const {return(true);}

  // Functions that actually do some work

  virtual void SetToConstant(double fv);

  virtual void RangeOfBasesWithSupportAtXyz(const NEWMAT::ColumnVector&       xyz,
                                            std::vector<unsigned int>&        first,
                                            std::vector<unsigned int>&        last) const
  {
    first[0] = first[1] = first[2] = 0;
    last[0] = CoefSz_x(); last[1] = CoefSz_y(); last[2] = CoefSz_z();
  }

  // Get the value of basis lmn at point xyz
  virtual double ValueOfBasisLmnAtXyz(const std::vector<unsigned int>&  lmn,
                                      const NEWMAT::ColumnVector&       xyz) const
  {
    return(0);
  }

  virtual std::vector<double> SubsampledVoxelSize(unsigned int               ss,
			                          std::vector<double>        vxs = std::vector<double>(),
				                  std::vector<unsigned int>  ms = std::vector<unsigned int>()) const
  {
    std::vector<unsigned int>  ssv(NDim(),ss);
    return(SubsampledVoxelSize(ssv,vxs,ms));
  }
  virtual std::vector<double> SubsampledVoxelSize(const std::vector<unsigned int>&  ss,
			                          std::vector<double>               vxs = std::vector<double>(),
				                  std::vector<unsigned int>         ms = std::vector<unsigned int>()) const;

  virtual std::vector<unsigned int> SubsampledMatrixSize(unsigned int               ss,
                                                         std::vector<unsigned int>  ms = std::vector<unsigned int>()) const
  {
    std::vector<unsigned int>  ssv(NDim(),ss);
    return(SubsampledMatrixSize(ssv,ms));
  }
  virtual std::vector<unsigned int> SubsampledMatrixSize(const std::vector<unsigned int>&  ss,
                                                         std::vector<unsigned int>         ms = std::vector<unsigned int>()) const;

  virtual void Update(FieldIndex fi);

  virtual NEWMAT::ReturnMatrix Jte(const NEWIMAGE::volume<float>&  ima1,
                                   const NEWIMAGE::volume<float>&  ima2,
                                   const NEWIMAGE::volume<char>    *mask) const;

  virtual NEWMAT::ReturnMatrix Jte(const std::vector<unsigned int>&  deriv,
                                   const NEWIMAGE::volume<float>&    ima1,
                                   const NEWIMAGE::volume<float>&    ima2,
                                   const NEWIMAGE::volume<char>      *mask) const;

  virtual NEWMAT::ReturnMatrix Jte(const NEWIMAGE::volume<float>&    ima,
                                   const NEWIMAGE::volume<char>      *mask) const;

  virtual NEWMAT::ReturnMatrix Jte(const std::vector<unsigned int>&  deriv,
                                   const NEWIMAGE::volume<float>&    ima,
                                   const NEWIMAGE::volume<char>      *mask) const;

  virtual boost::shared_ptr<MISCMATHS::BFMatrix> JtJ(const NEWIMAGE::volume<float>&    ima,
                                                     const NEWIMAGE::volume<char>      *mask,
                                                     MISCMATHS::BFMatrixPrecisionType  prec) const;

  virtual boost::shared_ptr<MISCMATHS::BFMatrix> JtJ(const NEWIMAGE::volume<float>&       ima1,
                                                     const NEWIMAGE::volume<float>&       ima2,
                                                     const NEWIMAGE::volume<char>         *mask,
                                                     MISCMATHS::BFMatrixPrecisionType     prec) const;

  virtual boost::shared_ptr<MISCMATHS::BFMatrix> JtJ(const std::vector<unsigned int>&       deriv, 
                                                     const NEWIMAGE::volume<float>&         ima,
                                                     const NEWIMAGE::volume<char>           *mask,
                                                     MISCMATHS::BFMatrixPrecisionType       prec) const;

  virtual boost::shared_ptr<MISCMATHS::BFMatrix> JtJ(const std::vector<unsigned int>&       deriv, 
                                                     const NEWIMAGE::volume<float>&         ima1,
                                                     const NEWIMAGE::volume<float>&         ima2,
                                                     const NEWIMAGE::volume<char>           *mask,
                                                     MISCMATHS::BFMatrixPrecisionType       prec) const;

  virtual boost::shared_ptr<MISCMATHS::BFMatrix> JtJ(const std::vector<unsigned int>&         deriv1,
                                                     const NEWIMAGE::volume<float>&           ima1,
                                                     const std::vector<unsigned int>&         deriv2,
                                                     const NEWIMAGE::volume<float>&           ima2,
                                                     const NEWIMAGE::volume<char>             *mask,
                                                     MISCMATHS::BFMatrixPrecisionType         prec) const;

  virtual boost::shared_ptr<MISCMATHS::BFMatrix> JtJ(const NEWIMAGE::volume<float>&        ima1,
                                                     const basisfield&                     bf2,
                                                     const NEWIMAGE::volume<float>&        ima2,
                                                     const NEWIMAGE::volume<char>          *mask,
                                                     MISCMATHS::BFMatrixPrecisionType      prec) const;


  virtual double MemEnergy() const;
  virtual double BendEnergy() const {throw BasisfieldException("dctfield::BendEnergy not yet implemented"); } // return(0.0);} // nvcc complained

  virtual NEWMAT::ReturnMatrix MemEnergyGrad() const;
  virtual NEWMAT::ReturnMatrix BendEnergyGrad() const
  {
    // Matrix  skrutt(1,1); // nvcc complained
    throw BasisfieldException("dctfield::BendEnergyGrad not yet implemented");
    // return(skrutt); // nvcc complained
  }

  virtual boost::shared_ptr<MISCMATHS::BFMatrix> MemEnergyHess(MISCMATHS::BFMatrixPrecisionType   prec) const;
  virtual boost::shared_ptr<MISCMATHS::BFMatrix> BendEnergyHess(MISCMATHS::BFMatrixPrecisionType   prec) const
  {
    // boost::shared_ptr<MISCMATHS::BFMatrix>    skrutt(new MISCMATHS::FullBFMatrix()); // nvcc complained
    throw BasisfieldException("dctfield::BendEnergyHess not yet implemented");
    // return(skrutt); // nvcc complained
  }

  virtual boost::shared_ptr<BASISFIELD::basisfield> ZoomField(const std::vector<unsigned int>&     psz,
                                                              const std::vector<double>&           pvxs,
                                                              std::vector<unsigned int>            order=std::vector<unsigned int>()) const;

};

} // End namespace BASISFIELD

#endif
