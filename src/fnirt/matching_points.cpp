// Definitions of class used for point-matching in FNIRT
//
// matching_points.cpp
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2008 University of Oxford 
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
#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
#include "newmat.h"
#include "newimage/newimageall.h"
#include "miscmaths/bfmatrix.h"
#include "basisfield/basisfield.h"
#include "basisfield/splinefield.h"
#include "basisfield/dctfield.h"
#include "warpfns/point_list.h"
#include "matching_points.h"

using namespace std;
using namespace boost;
using namespace NEWMAT;
using namespace NEWIMAGE;
using namespace BASISFIELD;
using namespace FNIRT;

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Member functions for class MatchingPoints
//
///////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Returns the difference in position in dimension dim for point-pair given by index
//
///////////////////////////////////////////////////////////////////////////////////////////////

double MatchingPoints::Diff(unsigned int indx, unsigned int dim, const BASISFIELD::basisfield& bf) const
{
  ColumnVector ref_pnt = _in_ref->Point(indx);
  ColumnVector rrp = _in_ref->RawPoint(indx);  // Raw Ref Point
  rrp = bf.mm2vox(3) * rrp;
  ColumnVector in_pnt = _in_in->Point(indx);
  return(in_pnt(dim+1) - (ref_pnt(dim+1) + bf(rrp(1),rrp(2),rrp(3))));
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Returns the sum-of-squared differences for the dimension given by dim
//
///////////////////////////////////////////////////////////////////////////////////////////////

double MatchingPoints::SSD(unsigned int dim, const basisfield& bf) const
{
  if (dim > 2) throw MatchingPointsException("SSD: dim must be 0--2");

  double ssd = 0.0;
  for (unsigned int i=0; i<_in_ref->NPoints(); i++) {
    ssd += MISCMATHS::Sqr(Diff(i,dim,bf));
  }
  return(ssd);
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Returns the gradient of the sum-of-squared differences for the dimension given by dim
//
///////////////////////////////////////////////////////////////////////////////////////////////

ReturnMatrix MatchingPoints::SSD_Gradient(unsigned int dim, const basisfield& bf) const
{
  if (dim > 2) throw MatchingPointsException("SSD: dim must be 0--2");

  ColumnVector grad(bf.CoefSz());
  grad = 0.0;
  for (unsigned int pti=0; pti<_in_ref->NPoints(); pti++) {
    ColumnVector rrp = _in_ref->RawPoint(pti);  // Raw Ref Point
    rrp = bf.mm2vox(3) * rrp;
    vector<unsigned int>  first(3,0);
    vector<unsigned int>  last(3,0);
    bf.RangeOfBasesWithSupportAtXyz(rrp,first,last);
    // cout << "rrp = " << rrp(1) << "  " << rrp(2) << "  " << rrp(3) << endl;
    // cout << "first = " << first[0] << "  " << first[1] << "  " << first[2] << endl;
    // cout << "last = " << last[0] << "  " << last[1] << "  " << last[2] << endl;
    vector<unsigned int>  ci(3,0);
    double diff = Diff(pti,dim,bf);
    for (ci[2]=first[2]; ci[2]<last[2]; ci[2]++) {
      for (ci[1]=first[1]; ci[1]<last[1]; ci[1]++) {
	for (ci[0]=first[0]; ci[0]<last[0]; ci[0]++) {
          unsigned int lindx = ci[2]*bf.CoefSz_x()*bf.CoefSz_y()+ci[1]*bf.CoefSz_x()+ci[0];
          grad(lindx+1) += -2.0 * diff * bf.ValueOfBasisLmnAtXyz(ci,rrp);
	}
      }
    }
  }
  grad.Release();
  return(grad);  
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Returns the Hessian of the sum-of-squared differences for the dimension given by dim
//
///////////////////////////////////////////////////////////////////////////////////////////////

boost::shared_ptr<BFMatrix> MatchingPoints::SSD_Hessian(unsigned int dim, const basisfield& bf, BFMatrixPrecisionType prec) const
{
  if (dim > 2) throw MatchingPointsException("SSD: dim must be 0--2");

  // Make sure matrix is of correct type
  boost::shared_ptr<BFMatrix>  hess;
  if (bf.HasGlobalSupport()) hess = boost::shared_ptr<BFMatrix>(new FullBFMatrix(bf.CoefSz(),bf.CoefSz()));
  else {
    if (prec == BFMatrixDoublePrecision) hess = boost::shared_ptr<BFMatrix>(new SparseBFMatrix<double>(bf.CoefSz(),bf.CoefSz()));
    else hess = boost::shared_ptr<BFMatrix>(new SparseBFMatrix<float>(bf.CoefSz(),bf.CoefSz()));
  }
  
  for (unsigned int pti=0; pti<_in_ref->NPoints(); pti++) {
    ColumnVector rrp = _in_ref->RawPoint(pti);  // Raw Ref Point
    rrp = bf.mm2vox(3) * rrp;
    vector<unsigned int>  first(3,0);
    vector<unsigned int>  last(3,0);
    bf.RangeOfBasesWithSupportAtXyz(rrp,first,last);
    vector<unsigned int>  rci(3,0);
    for (rci[2]=first[2]; rci[2]<last[2]; rci[2]++) {
      for (rci[1]=first[1]; rci[1]<last[1]; rci[1]++) {
	for (rci[0]=first[0]; rci[0]<last[0]; rci[0]++) {
          unsigned int row_lindx = rci[2]*bf.CoefSz_x()*bf.CoefSz_y()+rci[1]*bf.CoefSz_x()+rci[0];
          double row_val = bf.ValueOfBasisLmnAtXyz(rci,rrp);
	  vector<unsigned int>  cci(3,0);
          for (cci[2]=first[2]; cci[2]<last[2]; cci[2]++) {
            for (cci[1]=first[1]; cci[1]<last[1]; cci[1]++) {
	      for (cci[0]=first[0]; cci[0]<last[0]; cci[0]++) {
                unsigned int col_lindx = cci[2]*bf.CoefSz_x()*bf.CoefSz_y()+cci[1]*bf.CoefSz_x()+cci[0];
                double col_val = bf.ValueOfBasisLmnAtXyz(cci,rrp);
                hess->AddTo(row_lindx+1,col_lindx+1,2.0*row_val*col_val);
	      }
	    }
	  }
	}
      }
    }
  }
  return(hess);        
}
