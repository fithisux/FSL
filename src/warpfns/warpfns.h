/*  warpfns.h

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2001 University of Oxford  */

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

#if !defined(__warpfns_h)
#define __warpfns_h

#include <cstdlib>
#include <ctime>
#include <string>
#include <iostream>
#include <vector>
#define WANT_STREAM
#define WANT_MATH

#include "newmatap.h"
#include "newmatio.h"

#ifndef EXPOSE_TREACHEROUS
#define I_DEFINED_ET
#define EXPOSE_TREACHEROUS           // To allow us to use .sampling_mat()
#endif

#include "newimage/newimageall.h"
#include "basisfield/basisfield.h"

namespace NEWIMAGE {

class WarpFnsException: public std::exception
{
private:
  std::string m_msg;
public:
  WarpFnsException(const std::string& msg) throw(): m_msg(msg) {}

  virtual const char * what() const throw() {
    return string("warpfns::" + m_msg).c_str();
  }

  ~WarpFnsException() throw() {}
};


  //
  // Here stars declarations of functions inherited from initial
  // version of warpfns.h.
  //
  int affine2warp(const NEWMAT::Matrix& affmat, volume4D<float>& warpvol,
	  	  const volume<float>& outvol);

  int shift2warp(const volume<float>& shiftmap, 
	         volume4D<float>& warp, const string& shiftdir);

  int convertwarp_rel2abs(volume4D<float>& warpvol);
  int convertwarp_abs2rel(volume4D<float>& warpvol);

  int concat_warps(const volume4D<float>& prewarp, 
                   const volume4D<float>& postwarp,
		   volume4D<float>&       totalwarp);

  // default value for gammabar is for rad/s units; te in seconds; 
  //   lrgrad in rad/s/voxel
  volume<float> calc_sigloss(volume4D<float>& lrgrad, float te, 
			     float gammabar=0.5/M_PI);


  // try to determine if the warp is stored in absolute (vs relative) convention
  bool is_abs_convention(const volume4D<float>& warpvol);

  void jacobian_check(volume4D<float>& jvol,
		      ColumnVector& jacobian_stats, 
	  	      const volume4D<float>& warp,
		      float minJ, float maxJ, bool use_vol=true);

  volume4D<float> jacobian_check(ColumnVector& jacobian_stats, 
				 const volume4D<float>& warp,
				 float minJ, float maxJ);

  ColumnVector jacobian_quick_check(const volume4D<float>& warp,
				    float minJ, float maxJ); 


  void constrain_topology(volume4D<float>& warp, float minJ, float maxJ);

  void constrain_topology(volume4D<float>& warp);

//////////////////////////////////////////////////////////////////////////
//
// Here starts declarations of coordinate-transform functions
// that will be defined below.
//
// The general format of the NewimageCoord2NewimageCoord (here 
// abbreviated to N2N) is
// 
// trgt_coord = N2N(some_transforms,src_vol,trgt_vol,src_coord)
//
// The purpose of the routines is to supply a voxel-coordinate in one
// space (e.g. example_func) and obtain the corresponding voxel-
// coordinate in another space (e.g. standard space).
//
// The set of transforms that are indicated in "some_transforms" above
// Are as follows. Imagine we have four spaces a, b, c and d and that we
// have a linear transform M1 mapping a onto b, a non-linear transform
// w mapping c onto b (N.B. the order here) and a linear transform M2
// mapping c onto d.
//
// -------       -------       -------       -------
// |     |   M1  |     |   w   |     |   M2  |     |
// |  a  |------>|  b  |<------|  c  |------>|  d  |
// |     |  lin  |     | non-  |     |  lin  |     |
// -------       ------- lin   -------       -------
//
// A common example would be that a=example_func, b=highres
// c=standard, M1=example_func2highres.mat and w=highres2standard_warp.nii.gz
// In this example there is nothing corresponding to D or M2.
// Note also that the affine part of the mapping between b and c (i.e.
// highres2standard.mat) is incorporated into w.
//
// Let us now say we have a coordinate xf in the space of example_func (a), 
// and we want to map that to xs in standard space (c). We would then use a call
// that can be schematically described as
//
// xs = N2N(Matrix&             M1 = example_func2highres.mat,
//          volume4D<float>&    warps = highres2standard_warp.nii.gz,
//          bool                invert_warps = true
//          volume<D>&          src = example_func.nii.gz,
//          volume<S>&          dest = standard.nii.gz
//          ColumnVector&       xf)
//
// The important points to realise here is that first we use M1 to get
// from a->b, hence we pass M1 in first. Secondly we want to go from b->c,
// so we pass in w. *BUT* w maps c->b, which is why we have set the 
// invert_warps flag.
//
// Let us now assume we have a coordinate xs in standard space (c) and that
// we want to map it to a coordinate xf in functional space (a). We would
// then use a call like
//
// xf = N2N(volume4D<float>&    warps = highres2standard_warp.nii.gz,
//          bool                invert_warps = false,
//          Matrix&             M1 = inverse(example_func2highres.mat,
//          volume<D>&          src = standard.nii.gz,
//          volume<S>&          dest = example_func.nii.gz,
//          ColumnVector&       xs)
//
// The points to note here are that we first map from c->b, so we
// pass in w first. And this time w goes in the right direction
// so we set the invert_warps flag to false. After that we go from
// b->a, so we pass in THE INVERSE OF M1 to take us that step.
//
// These examples should hopefully explain the thoughts behind the
// various overloaded versions of N2N below. But before I stop I will
// also give you the same example in "proper code".
//
// volume<float>    funcvol; read_volume(funcvol,"example_func");
// volume<float>    stdvol; read_volume(stdvol,"standard");
// Matrix           M1 = read_ascii_matrix("example_func2highres.mat");
// FnirtFileReader  reader("highres2standard_warp");
//
// ColumnVector     xf(3);
// ColumnVector     xs(3);
//
// xf << 36 << 42 << 23;    // Coordinate 36,42,23 in example_func
// xs = NewimageCoord2NewimageCoord(M1,reader.FieldAsNewimageVolume4D(true),
//                                  true,funcvol,stdvol,xf);
//
//
// xs << 63 << 69 << 41;    // Coordinate 63,69,41 in standard
// xf = NewimageCoord2NewimageCoord(reader.FieldAsNewimageVolume4D(true),
//                                  false,M1.i(),stdvol,funcvol,xs);
//
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
//
//  Declarations
//
//////////////////////////////////////////////////////////////////////////

//
// For use e.g. when using affine mapping highres->standard
//
template <class D, class S>
NEWMAT::ReturnMatrix NewimageCoord2NewimageCoord(const NEWMAT::Matrix&        M,
                                                 const volume<D>&             srcvol,
                                                 const volume<S>&             destvol,
                                                 const NEWMAT::ColumnVector&  srccoord);

//
// For use e.g. when transforming (non-linearly) between highres and standard
//
template <class D, class S>
NEWMAT::ReturnMatrix NewimageCoord2NewimageCoord(const volume4D<float>&       warps,
                                                 bool                         inv_flag,
                                                 const volume<D>&             srcvol,
                                                 const volume<S>&             destvol,
                                                 const NEWMAT::ColumnVector&  srccoord);

//
// For use e.g. when transforming from standard space to functional space
//
template <class D, class S>
NEWMAT::ReturnMatrix NewimageCoord2NewimageCoord(const volume4D<float>&       warps,
						 bool                         inv_flag,
						 const NEWMAT::Matrix&        M,
						 const volume<D>&             srcvol,
						 const volume<S>&             destvol,
						 const NEWMAT::ColumnVector&  srccoord);

//
// For use e.g. when transforming from functional space to standard space
//
template <class D, class S>
NEWMAT::ReturnMatrix NewimageCoord2NewimageCoord(const NEWMAT::Matrix&        M,
						 const volume4D<float>&       warps,
						 bool                         inv_flag,
						 const volume<D>&             srcvol,
						 const volume<S>&             destvol,
						 const NEWMAT::ColumnVector&  srccoord);

//
// For use when doing the whole shabang. Whenever that might be.
//
template <class D, class S>
NEWMAT::ReturnMatrix NewimageCoord2NewimageCoord(const NEWMAT::Matrix&        M1,
						 const volume4D<float>&       warps,
						 bool                         inv_flag,
						 const NEWMAT::Matrix&        M2,
						 const volume<D>&             srcvol,
						 const volume<S>&             destvol,
						 const NEWMAT::ColumnVector&  srccoord);

//
// Internal function providing functionality for 
// all the overloaded functions above.
//
template <class D, class S>
int raw_newimagecoord2newimagecoord(const NEWMAT::Matrix         *M1,
				    const volume4D<float>        *warps,
				    bool                         inv_flag,
				    const NEWMAT::Matrix         *M2,
				    const volume<D>&             src,
				    const volume<S>&             trgt,
				    NEWMAT::ColumnVector&        coord);

//
// Function to calculate the inverse lookup for a single coordinate
// Uses newimage voxel coordinates everywhere and takes relative, mm warps
//
template <class T>
NEWMAT::ColumnVector inv_coord(const volume4D<float>&      warp, 
			       const volume<T>&            srcvol,
			       const NEWMAT::ColumnVector& coord);


//////////////////////////////////////////////////////////////////////////
//
//  Definitions
//
//////////////////////////////////////////////////////////////////////////

//
// For use e.g. when using affine mapping highres->standard
//
template <class D, class S>
NEWMAT::ReturnMatrix NewimageCoord2NewimageCoord(const NEWMAT::Matrix&        M,
						 const volume<D>&             srcvol,
						 const volume<S>&             destvol,
						 const NEWMAT::ColumnVector&  srccoord)
{
  // Check that matrix is kosher
  if (M.Nrows()!=4 || M.Ncols()!=4) imthrow("NewimageCoord2NewimageCoord: M must be a 4x4 matrix",11);
  // Allocate output coordinates
  NEWMAT::ColumnVector   tmp(4);
  tmp.Rows(1,3) = srccoord.Rows(1,3);
  tmp(4) = 1.0;
  // Do the conversion
  raw_newimagecoord2newimagecoord(&M,0,false,0,srcvol,destvol,tmp);
  // Return coordinate-vector of same size as incord
  NEWMAT::ColumnVector  destcord(srccoord);
  destcord.Rows(1,3) = tmp.Rows(1,3);
  destcord.Release();
  return(destcord);
}

//
// For use e.g. when transforming (non-linearly) between highres and standard
//
template <class D, class S>
NEWMAT::ReturnMatrix NewimageCoord2NewimageCoord(const volume4D<float>&       warps,
						 bool                         inv_flag,
						 const volume<D>&             srcvol,
						 const volume<S>&             destvol,
						 const NEWMAT::ColumnVector&  srccoord)
{
  // Check and prep warp-field
  if (warps.tsize() != 3) imthrow("NewimageCoord2NewimageCoord: warps must be a 4D volume with three points in fourth dimension",11);
  if (warps.nvoxels() <= 0) imthrow("NewimageCoord2NewimageCoord: warps must have a non-zero size",11);
  // Allocate output coordinates
  NEWMAT::ColumnVector   tmp(4);
  tmp.Rows(1,3) = srccoord.Rows(1,3);
  tmp(4) = 1.0;
  // Do the conversion
  raw_newimagecoord2newimagecoord(0,&warps,inv_flag,0,srcvol,destvol,tmp);
  // Return coordinate-vector of same size as incord
  NEWMAT::ColumnVector  destcord(srccoord);
  destcord.Rows(1,3) = tmp.Rows(1,3);
  destcord.Release();
  return(destcord);
}

//
// For use e.g. when transforming from standard space to functional space
//
template <class D, class S>
NEWMAT::ReturnMatrix NewimageCoord2NewimageCoord(const volume4D<float>&       warps,
						 bool                         inv_flag,
						 const NEWMAT::Matrix&        M,
						 const volume<D>&             srcvol,
						 const volume<S>&             destvol,
						 const NEWMAT::ColumnVector&  srccoord)
{
  // Check that matrix is kosher
  if (M.Nrows()!=4 || M.Ncols()!=4) imthrow("NewimageCoord2NewimageCoord: M must be a 4x4 matrix",11);
  // Check and prep warp-field
  if (warps.tsize() != 3) imthrow("NewimageCoord2NewimageCoord: warps must be a 4D volume with three points in fourth dimension",11);
  if (warps.nvoxels() <= 0) imthrow("NewimageCoord2NewimageCoord: warps must have a non-zero size",11);
  // Allocate output coordinates
  NEWMAT::ColumnVector   tmp(4);
  tmp.Rows(1,3) = srccoord.Rows(1,3);
  tmp(4) = 1.0;
  // Do the conversion
  raw_newimagecoord2newimagecoord(0,&warps,inv_flag,&M,srcvol,destvol,tmp);
  // Return coordinate-vector of same size as incord
  NEWMAT::ColumnVector  destcord(srccoord);
  destcord.Rows(1,3) = tmp.Rows(1,3);
  destcord.Release();
  return(destcord);
}

//
// For use e.g. when transforming from functional space to standard space
//
template <class D, class S>
NEWMAT::ReturnMatrix NewimageCoord2NewimageCoord(const NEWMAT::Matrix&        M,
						 const volume4D<float>&       warps,
						 bool                         inv_flag,
						 const volume<D>&             srcvol,
						 const volume<S>&             destvol,
						 const NEWMAT::ColumnVector&  srccoord)
{
  // Check that matrix is kosher
  if (M.Nrows()!=4 || M.Ncols()!=4) imthrow("NewimageCoord2NewimageCoord: M must be a 4x4 matrix",11);
  // Check and prep warp-field
  if (warps.tsize() != 3) imthrow("NewimageCoord2NewimageCoord: warps must be a 4D volume with three points in fourth dimension",11);
  if (warps.nvoxels() <= 0) imthrow("NewimageCoord2NewimageCoord: warps must have a non-zero size",11);
  // Allocate output coordinates
  NEWMAT::ColumnVector   tmp(4);
  tmp.Rows(1,3) = srccoord.Rows(1,3);
  tmp(4) = 1.0;
  // Do the conversion
  raw_newimagecoord2newimagecoord(&M,&warps,inv_flag,0,srcvol,destvol,tmp);
  // Return coordinate-vector of same size as incord
  NEWMAT::ColumnVector  destcord(srccoord);
  destcord.Rows(1,3) = tmp.Rows(1,3);
  destcord.Release();
  return(destcord);
}

//
// For use when doing the whole shabang. Whenever that might be.
//
template <class D, class S>
NEWMAT::ReturnMatrix NewimageCoord2NewimageCoord(const NEWMAT::Matrix&        M1,
						 const volume4D<float>&       warps,
						 bool                         inv_flag,
						 const NEWMAT::Matrix&        M2,
						 const volume<D>&             srcvol,
						 const volume<S>&             destvol,
						 const NEWMAT::ColumnVector&  srccoord)
{
  // Check that matrices are kosher
  if (M1.Nrows()!=4 || M1.Ncols()!=4) imthrow("NewimageCoord2NewimageCoord: M1 must be a 4x4 matrix",11);
  if (M2.Nrows()!=4 || M2.Ncols()!=4) imthrow("NewimageCoord2NewimageCoord: M2 must be a 4x4 matrix",11);
  // Check and prep warp-field
  if (warps.tsize() != 3) imthrow("NewimageCoord2NewimageCoord: warps must be a 4D volume with three points in fourth dimension",11);
  if (warps.nvoxels() <= 0) imthrow("NewimageCoord2NewimageCoord: warps must have a non-zero size",11);
  // Allocate output coordinates
  NEWMAT::ColumnVector   tmp(4);
  tmp.Rows(1,3) = srccoord.Rows(1,3);
  tmp(4) = 1.0;
  // Do the conversion
  raw_newimagecoord2newimagecoord(&M1,&warps,inv_flag,&M2,srcvol,destvol,tmp);
  // Return coordinate-vector of same size as incord
  NEWMAT::ColumnVector  destcord(srccoord);
  destcord.Rows(1,3) = tmp.Rows(1,3);
  destcord.Release();
  return(destcord);
}

template <class D, class S>
int raw_newimagecoord2newimagecoord(const NEWMAT::Matrix         *M1,
				    const volume4D<float>        *warps,
				    bool                         inv_flag,
				    const NEWMAT::Matrix         *M2,
				    const volume<D>&             src,
				    const volume<S>&             trgt,
				    NEWMAT::ColumnVector&        coord)
{
  int  rval = 1;
  //
  // First go from voxel-coordinate in src to mm-space of warps
  //
  coord = src.sampling_mat()*coord;
  if (M1) coord = (*M1)*coord;
  //
  // Then find displacement-vector at that point
  //
  if (warps) {
    if (inv_flag) {
      // Here we will use a fake volume to take us back and 
      // forth between mm and newimage coordinates. This is
      // just so that we shall be able to pass voxel-coordinates
      // into inv_coord even though here we are really in mm-space.
      volume<float>  fake_vol = (*warps)[0];
      coord = fake_vol.sampling_mat().i() * coord;  // mm->vox in fake-space
      coord = inv_coord(*warps,fake_vol,coord);     // vox_in_fake->vox_in_ref_of_warps
      coord = warps->sampling_mat() * coord;        // vox->mm in ref_of_warps-space
    }
    else {
      NEWMAT::ColumnVector  vd_coord = warps->sampling_mat().i() * coord;
      extrapolation oldex = warps->getextrapolationmethod();
      if (oldex!=periodic) warps->setextrapolationmethod(extraslice);
      NEWMAT::ColumnVector  dvec(4);
      dvec = 0.0;
      dvec(1) = (*warps)[0].interpolate(vd_coord(1),vd_coord(2),vd_coord(3));
      dvec(2) = (*warps)[1].interpolate(vd_coord(1),vd_coord(2),vd_coord(3));
      dvec(3) = (*warps)[2].interpolate(vd_coord(1),vd_coord(2),vd_coord(3));
      coord += dvec;
      rval = (warps->in_bounds(float(vd_coord(1)),float(vd_coord(2)),float(vd_coord(3)))) ? 1 : -1;   // Indicate invalid warp-value
      warps->setextrapolationmethod(oldex);
    }
  }
  //
  // Finally go to voxel-space of trgt
  //
  if (M2) coord = (*M2)*coord;
  coord = trgt.sampling_mat().i()*coord;

  return(rval);
}


template <class T>
NEWMAT::ColumnVector inv_coord(const volume4D<float>&      warp, 
                               const volume<T>&            srcvol,
		               const NEWMAT::ColumnVector& coord)
{
  // coord is in source (x1) space - same as srcvol
  // need the srcvol to work out the voxel<->mm conversion for the x1 space
  int N=5;
  int N_2=(int)(N/2);

  NEWMAT::ColumnVector coord4(4);
  coord4 << coord(1) << coord(2) << coord(3) << 1.0;

  volume<float> idx(N,N,N);  // NB: different coordinate conventions to matlab
  idx=0.0f;
  int c=1;
  for (int z=0; z<N; z++) {
    for (int y=0; y<N; y++) {
      for (int x=0; x<N; x++) {
	idx(x,y,z)=c++;
      }
    }
  }

  // Form matrix M which stores trilinear coefficients needed to interpolate inverse warpfield
  // The coordlist is the list of x2 coordinates mapping into the neighbourhood (in voxel coords)
  NEWMAT::Matrix M;
  NEWMAT::Matrix coordlist;
  float dx,dy,dz;
  float dist, mindist=Max(srcvol.maxx(),Max(srcvol.maxy(),srcvol.maxz()));
  mindist *= mindist;  // conservative estimate
  float minptx=-1, minpty=-1, minptz=-1;
  // loop around all voxels in ref space: x2
  for (int z2=0; z2<=warp.maxz(); z2++) {
    for (int y2=0; y2<=warp.maxy(); y2++) {
      for (int x2=0; x2<=warp.maxx(); x2++) {
	// if warp points towards voxel in the neighbour of the target coord
	dx=(warp[0](x2,y2,z2)+x2*warp.xdim())/srcvol.xdim()-coord4(1);
	dy=(warp[1](x2,y2,z2)+y2*warp.ydim())/srcvol.ydim()-coord4(2);
	dz=(warp[2](x2,y2,z2)+z2*warp.zdim())/srcvol.zdim()-coord4(3);
	dist=dx*dx+dy*dy+dz*dz;  // in voxels - might be better in mm?
	if (dist<mindist) { mindist=dist; minptx=x2; minpty=y2; minptz=z2; }
	if ((fabs(dx)<N_2) && (fabs(dy)<N_2) && (fabs(dz)<N_2)) {
	  addrow(M,N*N*N);
	  // add current voxel coord (x2) to coordlist
	  addrow(coordlist,3);
	  coordlist(coordlist.Nrows(),1)=x2;
	  coordlist(coordlist.Nrows(),2)=y2;
	  coordlist(coordlist.Nrows(),3)=z2;
	  // loop around all voxels in target coords' neighbourhood: x1
	  for (int z1=-N_2; z1<=N_2; z1++) {
	    for (int y1=-N_2; y1<=N_2; y1++) {
	      for (int x1=-N_2; x1<=N_2; x1++) {
		if ((fabs(dx-x1)<1.0) && (fabs(dy-y1)<1.0) && (fabs(dz-z1)<1.0)) {
		  //   M(newrow,idx(x1)) = (1-fabs(dx))*(1-fabs(dy))*(1-fabs(dz))
		  M(M.Nrows(),MISCMATHS::round(idx(x1+N_2,y1+N_2,z1+N_2))) = (1.0-fabs(dx-x1))*(1.0-fabs(dy-y1))*(1.0-fabs(dz-z1));
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  int ncols=N*N*N;  
  NEWMAT::Matrix L;
  // form regularisation matrix (L) - isotropic assumption for now - probably unimportant
  for (int z=0; z<N; z++) {
    for (int y=0; y<N; y++) {
      for (int x=0; x<N; x++) {
	if ((x>0) && (x<N-1)) {
	  addrow(L,ncols);
	  L(L.Nrows(),MISCMATHS::round(idx(x,y,z)))=2;
	  L(L.Nrows(),MISCMATHS::round(idx(x-1,y,z)))=-1; 
	  L(L.Nrows(),MISCMATHS::round(idx(x+1,y,z)))=-1;
	}
	if ((y>0) && (y<N-1)) {
	  addrow(L,ncols);
	  L(L.Nrows(),MISCMATHS::round(idx(x,y,z)))=2;
	  L(L.Nrows(),MISCMATHS::round(idx(x,y-1,z)))=-1; 
	  L(L.Nrows(),MISCMATHS::round(idx(x,y+1,z)))=-1;
	}
	if ((z>0) && (z<N-1)) {
	  addrow(L,ncols);
	  L(L.Nrows(),MISCMATHS::round(idx(x,y,z)))=2;
	  L(L.Nrows(),MISCMATHS::round(idx(x,y,z-1)))=-1; 
	  L(L.Nrows(),MISCMATHS::round(idx(x,y,z+1)))=-1;
	}
      }
    }
  }

  // now calculate new coordinate (in voxel coords)
  NEWMAT::ColumnVector newcoord(4);
  if (M.Nrows()>8) {
    // normalise both M and L (on a per-voxel basis)
    M *= 1;
    //L *= 1*M.Nrows()/L.Nrows();
    float lambda=0.5;   // this is the default in invwarp.cc
    NEWMAT::ColumnVector coordsx, coordsy, coordsz;
    NEWMAT::CroutMatrix X = M.t()*M + lambda*L.t()*L;  // carries out LU decomposition - for efficiency
    coordsx = X.i()*M.t()*coordlist.SubMatrix(1,coordlist.Nrows(),1,1);
    coordsy = X.i()*M.t()*coordlist.SubMatrix(1,coordlist.Nrows(),2,2);
    coordsz = X.i()*M.t()*coordlist.SubMatrix(1,coordlist.Nrows(),3,3);
    newcoord << coordsx(MISCMATHS::round(idx(N_2,N_2,N_2))) 
	     << coordsy(MISCMATHS::round(idx(N_2,N_2,N_2))) 
	     << coordsz(MISCMATHS::round(idx(N_2,N_2,N_2))) << 1.0;
  } else {
    // just copy the nearest value of relative warp (wherever it was) and add that to the existing voxel
    //  this should deal with any constant translation and maybe even rotation...
    newcoord << (warp[0].interpolate(minptx,minpty,minptz) + coord4(1)*srcvol.xdim())/warp.xdim()
	     << (warp[1].interpolate(minptx,minpty,minptz) + coord4(2)*srcvol.ydim())/warp.ydim()
	     << (warp[2].interpolate(minptx,minpty,minptz) + coord4(3)*srcvol.zdim())/warp.zdim() << 1.0;
  }

  if (coord.Nrows()==3) {
    NEWMAT::ColumnVector nc3(3);  nc3 << newcoord(1) << newcoord(2) << newcoord(3);  newcoord=nc3; 
  }
  return newcoord;
}


/////////////////////////////////////////////////////////////////////
//
// Here starts declarations of the various overloaded versions
// of raw_general_transform.
//
/////////////////////////////////////////////////////////////////////


  //
  // Here starts declarations of the templated functions that will
  // be defined below.
  //

  template <class T>    
  void raw_general_transform(// Input
                             const volume<T>&         s,          // Input volume
                             const NEWMAT::Matrix&    A,          // Mapping of t onto in
                             const volume4D<float>&   d,          // Displacement fields
                             const vector<int>&       defdir,     // Directions of displacements. 
                             const vector<int>&       derivdir,   // Directions of derivatives
                             const NEWMAT::Matrix     *TT,        // Mapping of out onto t
                             const NEWMAT::Matrix     *M,         // Mapping of in onto s
                             // Output
                             volume<T>&               out,        // Output volume
                             volume4D<T>&             deriv,      // Partial derivative directions
                             volume<char>             *invol);    // Mask indicating what voxels fell inside original volume

/////////////////////////////////////////////////////////////////////
//
// The following two routines are slightly simplified interfaces
// such that users should not have to pass in zero-pointers when
// they do not want to have an "infov" mask outpuy, or when there
// are only a template and an inout volume such that we have no
// need for M or TT.
//
/////////////////////////////////////////////////////////////////////

  template<class T>
  void raw_general_transform(// Input
                             const volume<T>&         vin,        // Input volume
                             const NEWMAT::Matrix&    aff,        // 4x4 affine transformation matrix
                             const volume4D<float>&   df,         // Displacement fields
                             const vector<int>&       defdir,     // Directions of displacements. 
                             const vector<int>&       derivdir,   // Directions of derivatives
                             // Output
                             volume<T>&               vout,       // Output volume
                             volume4D<T>&             deriv);     // Partial derivative directions

  template <class T>
  void raw_general_transform(// Input
                             const volume<T>&         vin,        // Input volume
                             const NEWMAT::Matrix&    aff,        // 4x4 affine transformation matrix
                             const volume4D<float>&   df,         // Displacement fields
                             const vector<int>&       defdir,     // Directions of displacements. 
                             const vector<int>&       derivdir,   // Directions of derivatives
                             // Output
                             volume<T>&               vout,       // Output volume
                             volume4D<T>&             deriv,      // Partial derivative directions
                             volume<char>&            invol);     // Mask indicating what voxels fell inside original volume (or valid extrapolation)

  template<class T>
  void apply_warp(// Input
                  const volume<T>&        vin,         // Input volume
                  const NEWMAT::Matrix&   A,           // 4x4 affine transform
                  const volume4D<float>   d,           // Displacement fields
                  const NEWMAT::Matrix&   TT,
                  const NEWMAT::Matrix&   M,
                  // Output
                  volume<T>&              vout);       // Resampled output volume

  template<class T>
  void apply_warp(// Input
                  const volume<T>&        vin,         // Input volume
                  const NEWMAT::Matrix&   A,           // 4x4 affine transform
                  const volume4D<float>   d,           // Displacement fields
                  const NEWMAT::Matrix&   TT,
                  const NEWMAT::Matrix&   M,
                  // Output
                  volume<T>&              vout,        // Resampled output volume
		  volume<char>&           mask);       // Set when inside original volume


/////////////////////////////////////////////////////////////////////
//
// The following three routines are interafaces to mimick the old
// functions apply_warp and raw_apply_warp. These are used manily 
// for resampling of images that are typically related to the
// input to e.g. fnirt or fugue through some rigid matrix M. 
// Examples of M would be a rigid mapping between a subjects 
// functional volumes and his/her structural or between a subjects
// functional volumes and his/her field map. 
//
/////////////////////////////////////////////////////////////////////

  template <class T>
  int apply_warp(const volume<T>&        invol,
                 volume<T>&              outvol,
                 const volume4D<float>&  warpvol);

  template <class T>
  int apply_warp(const volume<T>&        invol, 
                 volume<T>&              outvol,
	         const volume4D<float>&  warpvol, 
	         const NEWMAT::Matrix&   premat, 
                 const NEWMAT::Matrix&   postmat);

  template <class T>
  int raw_apply_warp(const volume<T>&        invol, 
                     volume<T>&              outvol,
		     const volume4D<float>&  warpvol, 
		     const NEWMAT::Matrix&   premat, 
                     const NEWMAT::Matrix&   postmat);

  template <class T>
  void affine_transform(// Input
                        const volume<T>&         vin,
                        const NEWMAT::Matrix&    aff,
                        // Output
                        volume<T>&               vout);

  template <class T>
  void affine_transform(// Input
                        const volume<T>&         vin,
                        const NEWMAT::Matrix&    aff,
                        // Output
                        volume<T>&               vout,
                        volume<char>&            invol);

  template <class T>
  void affine_transform_3partial(// Input
                                 const volume<T>&       vin,
                                 const NEWMAT::Matrix&  aff,
                                 // Output
                                 volume<T>&             vout,
                                 volume4D<T>&           deriv);

  template <class T>
  void affine_transform_3partial(// Input
                                 const volume<T>&       vin,
                                 const NEWMAT::Matrix&  aff,
                                 // Output
                                 volume<T>&             vout,
                                 volume4D<T>&           deriv,
                                 volume<char>&          invol);
  template <class T>
  void displacement_transform_1D(// Input
                                 const volume<T>&       vin,
                                 const NEWMAT::Matrix&  aff,
                                 const volume<float>&   df,
                                 int                    defdir,
                                 // Output
                                 volume<T>&             vout);

  template <class T>
  void displacement_transform_1D(// Input
                                 const volume<T>&       vin,
                                 const NEWMAT::Matrix&  aff,
                                 const volume<float>&   df,
                                 int                    defdir,
                                 // Output
                                 volume<T>&             vout,
                                 volume<char>&          invol);

  template <class T>
  void displacement_transform_1D_3partial(// Input
                                          const volume<T>&       vin,
                                          const NEWMAT::Matrix&  aff,
                                          const volume<float>&   df,
                                          int                    dir,
                                          // Output
                                          volume<T>&             vout,
                                          volume4D<T>&           deriv);

  template <class T>
  void displacement_transform_1D_3partial(// Input
                                          const volume<T>&       vin,
                                          const NEWMAT::Matrix&  aff,
                                          const volume<float>&   df,
                                          int                    dir,
                                          // Output
                                          volume<T>&             vout,
                                          volume4D<T>&           deriv,
                                          volume<char>&          invol);

  template <class T>
  void general_transform(// Input
                         const volume<T>&         vin,
                         const NEWMAT::Matrix&    aff,
                         const volume4D<float>&   df,
                         // Output
                         volume<T>&               vout);

  template <class T>	       
  void general_transform(// Input
                         const volume<T>&         vin,
                         const NEWMAT::Matrix&    aff,
                         const volume4D<float>&   df,
                         // Output
                         volume<T>&               vout,
                         volume<char>&            invol);

  template <class T>
  void general_transform_3partial(// Input
                                  const volume<T>&         vin,
                                  const NEWMAT::Matrix&    aff,
                                  const volume4D<float>&   df,
                                  // Output
                                  volume<T>&               vout,
                                  volume4D<T>&             deriv);

  template <class T>
  void general_transform_3partial(// Input
                                  const volume<T>&         vin,
                                  const NEWMAT::Matrix&    aff,
                                  const volume4D<float>&   df,
                                  // Output
                                  volume<T>&               vout,
                                  volume4D<T>&             deriv,
                                  volume<char>&            invol);

  template <class T> 
  void get_displacement_fields(// Input
                                  const volume<T>&         vin,
                                  const NEWMAT::Matrix&    aff,
                                  const volume4D<float>&   dfin,
                                  // Output
                                  volume4D<float>&         dfout);

  ///////////////////////////////////////////////////////////////////////////
  // IMAGE PROCESSING ROUTINES
  ///////////////////////////////////////////////////////////////////////////

  // General Transform
  //
  // The routine "raw_general_transform" is the heart of the "warping" 
  // functions. It provides functionality for calculating warped images
  // and partial derivatives in warped space for use by routines for 
  // non-linear registration as well as distortion correction. In addition
  // it is also used for final resampling of images given a pre-determined
  // displacement field.
  // 
  // In the most general case we might have registered some volume i to a
  // template s that already had an affine transformation matrix A mapping
  // i onto s. The non-linear mapping of i onto s is given by a displacement
  // field d. We may in addition also have a volume f that maps linearly onto
  // the volume in through a linear transform M. A typical example would be
  // that s is e.g. the avg152 (implementing the MNI space), i is a structural
  // from some subject and f is a functional volume from that same subject.
  // A is a matrix generated by flirt mapping i onto s, M is another matrix
  // generated by flirt mapping the functional onto the structural and d is 
  // a displacement field calculated by fnirt.
  // Finally there is another volume out, which defines the space to which we
  // ultimately want to resample f (or i, if there is no f). There is an
  // affine transform T that maps s onto out. An example of out might be the
  // Talairach space and T might be a "known" matrix that effects an approximate
  // mapping MNI->Talairach. Another example of out might simply be a volume in
  // the space of s, but with a different voxel/matrix size. For example if one
  // wants to use a template with a 1mm isotropic resolution (for high resolution
  // nonlinear registration) but wants to resample the functional data to a
  // volume with 2mm resolution (because 1mm might be a little over the top
  // for functional data).
  //
  // In the code I have retained the notation I have sketched above. To recap
  // 
  // volume<T>       f   // "Final" volume in chain. The volume that we want to map some 
  //                     // cordinate x_out (in space of out) into so that we can can interpolate
  //                     // intensity values from s and write into out
  // Matrix          A;  // Affine mapping of i onto s
  // volume4D<float> d;  // Non-linear mapping of i onto s
  // volume<T>       s;  // Volume used as template. A coordinate x_i in volume i is related
  //                     // to a coordinate x_s in s as x_i = inv(B_i)*inv(A)*B_s*x_s + inv(B_i)*d(x_s)
  //                     // where B_i and B_s are voxel->mm transforms for i and s respectively.
  // volume<T>       out // Volume into which we ultimately want to map f (or i, if no f)
  // Matrix          T   // Affine (or subset of) mapping of t onto s
  // Matrix          M   // Affine (typically rigid sub-set of) mapping of in onto s
  //
  // Together with some other options (such as derivatives or no derivatives, 3D or 1D non-linear
  // transform etc) this all makes the calling interface a little messy. There is therefore a set
  // of alterantive calls with a reduced interface that will be more convenient for many
  // applications.
  //

/////////////////////////////////////////////////////////////////////
//
// This is as raw as it gets, with pointers and all
//
/////////////////////////////////////////////////////////////////////

  template <class T>    
  void raw_general_transform(// Input
                             const volume<T>&         f,          // Input volume
                             const NEWMAT::Matrix&    A,          // 4x4 affine transformation matrix
                             const volume4D<float>&   d,          // Displacement fields (also defines space of t). Note that
                                                                  // these are "relative" fields in units of mm.
                             const vector<int>&       defdir,     // Directions of displacements. 
                             const vector<int>&       derivdir,   // Directions of derivatives
                             const NEWMAT::Matrix     *TT,        // Mapping of out onto t
                             const NEWMAT::Matrix     *M,         // Mapping of in onto s
                             // Output
                             volume<T>&               out,        // Output volume
                             volume4D<T>&             deriv,      // Partial derivatives. Note that the derivatives are in units
                                                                  // "per voxel".
                             volume<char>             *valid)     // Mask indicating what voxels fell inside original fov OR that extrapolation is valid
  {
    if (int(defdir.size()) != d.tsize()) {imthrow("NEWIMAGE::raw_general_transform: Mismatch in input. defdir.size() must equal d.tsize()",11);}
    for (int i=0; i<int(defdir.size()); i++) {
      if (defdir[i] < 0 || defdir[i] > 2) {imthrow("NEWIMAGE::raw_general_transform: Mismatch in input. Displacements can only be specified for directions 0,1 or 2.",11);}
    }
    if (int(derivdir.size()) != deriv.tsize()) {imthrow("NEWIMAGE::raw_general_transform: Mismatch in input. derivdir.size() must equal deriv.tsize()",11);}
    for (int i=0; i<int(derivdir.size()); i++) {
      if (derivdir[i] < 0 || derivdir[i] > 2) {imthrow("NEWIMAGE::raw_general_transform: Mismatch in input. Displacements can only be specified for directions 0,1 or 2.",11);}
    }
    if (out.nvoxels() <= 0) {imthrow("NEWIMAGE::raw_general_transform: Size of vout must be set",8);}
    if (derivdir.size() && !samesize(out,deriv[0]))  {imthrow("NEWIMAGE::raw_general_transform: vout and deriv must have same dimensions",11);}
    if (valid && !samesize(out,*valid)) {imthrow("NEWIMAGE::raw_general_transform: vout and valid must have same dimensions",11);}
    if (A.Nrows() != 4 || A.Ncols() != 4) {imthrow("NEWIMAGE::raw_general_transform: A must be 4x4 matrix",11);}
    if (TT) {
      if (TT->Nrows() != 4 || TT->Ncols() != 4) {imthrow("NEWIMAGE::raw_general_transform: T must be 4x4 matrix",11);}
    }
    if (M) {
      if (M->Nrows() != 4 || M->Ncols() != 4) {imthrow("NEWIMAGE::raw_general_transform: M must be 4x4 matrix",11);}
    }
  
    extrapolation oldex = f.getextrapolationmethod();
    extrapolation d_oldex = extraslice;  // Assign arbitrary value to silence compiler
    vector<bool>  d_old_epvalidity;
    if ((oldex==boundsassert) || (oldex==boundsexception)) {f.setextrapolationmethod(constpad);}
    if (d.tsize()) {
      d_oldex = d.getextrapolationmethod();
      if (oldex==periodic) d.setextrapolationmethod(periodic);
      else d.setextrapolationmethod(extraslice);
      d_old_epvalidity = d.getextrapolationvalidity();
      vector<bool> epvalidity = f.getextrapolationvalidity();
      d.setextrapolationvalidity(epvalidity[0],epvalidity[1],epvalidity[2]);
    }

    // Repackage info about displacement directions in more convenient form
    int xd, yd, zd;
    xd=yd=zd= -1;
    for (int i=0; i<int(defdir.size()); i++) {
      if (defdir[i]==0) {xd=i;} else if (defdir[i]==1) {yd=i;} else {zd=i;}
    }
    // Repackage info about derivative directions in more convenient form
    int xp, yp, zp;
    xp=yp=zp= -1;
    for (int i=0; i<int(derivdir.size()); i++) {
      if (derivdir[i]==0) {xp=i;} else if (derivdir[i]==1) {yp=i;} else {zp=i;}
    }
    
    // Create a matrix iM mapping from voxel coordinates in volume out
    // to voxel-coordinates in volume s (same space as d).

    float T11=0.0, T12=0.0, T13=0.0, T14=0.0, T21=0.0, T22=0.0;
    float T23=0.0, T24=0.0, T31=0.0, T32=0.0, T33=0.0, T34=0.0;
    NEWMAT::Matrix iT(4,4);
    if (TT) {
      if (d.tsize()) iT = d.sampling_mat().i() * TT->i() * out.sampling_mat();
      else iT = out.sampling_mat().i() * TT->i() * out.sampling_mat();
    }
    else {
      if (d.tsize()) iT = d.sampling_mat().i() * out.sampling_mat();
      else iT = IdentityMatrix(4);
    }
    //
    // If iT is different from the identity matrix we should indicate this,
    // and also specify values for T11 , T12 etc.
    //
    bool useiT = false;
    if ((iT-IdentityMatrix(4)).MaximumAbsoluteValue() > 1e-6) {
      useiT = true;
      T11=iT(1,1), T12=iT(1,2), T13=iT(1,3), T14=iT(1,4);
      T21=iT(2,1), T22=iT(2,2), T23=iT(2,3), T24=iT(2,4);
      T31=iT(3,1), T32=iT(3,2), T33=iT(3,3), T34=iT(3,4);
    }

    // Create a matrix iA mapping from voxel coordinates in volume t to
    // mm-coordinates in volume i. Affine part only of course 

    NEWMAT::Matrix iA = A.i();

    // cout << "It is still I" << endl;
    // if (f.left_right_order()==FSL_NEUROLOGICAL) {iA = f.swapmat(-1,2,3) * iA;}      // Swap if input neurological
    // if (out.left_right_order()==FSL_NEUROLOGICAL) {iA = iA * out.swapmat(-1,2,3);}  // Swap if output neurological

    if (defdir.size()) iA = iA * d[0].sampling_mat();  // If we have a displacement field
    else iA = iA * out.sampling_mat();                 // Else
    
    float A11=iA(1,1), A12=iA(1,2), A13=iA(1,3), A14=iA(1,4);
    float A21=iA(2,1), A22=iA(2,2), A23=iA(2,3), A24=iA(2,4);
    float A31=iA(3,1), A32=iA(3,2), A33=iA(3,3), A34=iA(3,4); 

    // Create a matrix mapping from mm-coordinates in volume i
    // to voxel coordinates in volume f. If the matrix M is empty
    // this is simply a mm->voxel mapping for volume f

    NEWMAT::Matrix iM(4,4);
    if (M) {
      iM = f.sampling_mat().i() * M->i();
    }
    else {
      iM = f.sampling_mat().i();
    }
    float M11=iM(1,1), M12=iM(1,2), M13=iM(1,3), M14=iM(1,4);
    float M21=iM(2,1), M22=iM(2,2), M23=iM(2,3), M24=iM(2,4);
    float M31=iM(3,1), M32=iM(3,2), M33=iM(3,3), M34=iM(3,4); 

    float o1,o2,o3;
  
    // I have put some "outer if's" leading to code multiplication here.
    // I have done so to ensure that we don't pay a performance penalty
    // For example for the cases where we want to use the routine to get
    // derivatives for the affine only case, or when we want to perform
    // a general transformation without calculating any derivatives.

    if (!defdir.size()) { // If we have an affine only transform
      // Affine only means we can combine all three matrices into one
      NEWMAT::Matrix iB = iM*iA*iT;
      A11=iB(1,1), A12=iB(1,2), A13=iB(1,3), A14=iB(1,4);
      A21=iB(2,1), A22=iB(2,2), A23=iB(2,3), A24=iB(2,4);
      A31=iB(3,1), A32=iB(3,2), A33=iB(3,3), A34=iB(3,4); 
      if (!derivdir.size()) { // If we don't need to calculate derivatives
        for (int z=0; z<out.zsize(); z++) { 
          for (int x=0; x<out.xsize(); x++) {
	    o1=x*A11 + z*A13 + A14;  // y=0
	    o2=x*A21 + z*A23 + A24;  // y=0
	    o3=x*A31 + z*A33 + A34;  // y=0
 	    for (int y=0; y<out.ysize(); y++) {
	      out(x,y,z) = ((T) f.interpolate(o1,o2,o3));
              if (valid) {valid->operator()(x,y,z) = (f.valid(o1,o2,o3)) ? 1 : 0;}
	      o1 += A12;
	      o2 += A22;
	      o3 += A32;
	    }
          }
        }
      }
      else { // If we need derivatives in at least one direction
        for (int z=0; z<out.zsize(); z++) { 
          for (int x=0; x<out.xsize(); x++) { 
	    o1=x*A11 + z*A13 + A14;  // y=0
	    o2=x*A21 + z*A23 + A24;  // y=0
	    o3=x*A31 + z*A33 + A34;  // y=0
 	    for (int y=0; y<out.ysize(); y++) {
              if (derivdir.size() == 1) { // If we want a single partial
                float tmp;
	        out(x,y,z) = ((T) f.interp1partial(o1,o2,o3,derivdir[0],&tmp));
                deriv(x,y,z,0) = ((T) tmp);
	      }
              else { // More than one derivative
		float tmp1,tmp2,tmp3;
		out(x,y,z) = ((T) f.interp3partial(o1,o2,o3,&tmp1,&tmp2,&tmp3));
                if (!(xp<0)) {deriv(x,y,z,xp)=((T) tmp1);}
                if (!(yp<0)) {deriv(x,y,z,yp)=((T) tmp2);}
                if (!(zp<0)) {deriv(x,y,z,zp)=((T) tmp3);}
	      }
              if (valid) {valid->operator()(x,y,z) = (f.valid(o1,o2,o3)) ? 1 : 0;}
	      o1 += A12;
	      o2 += A22;
	      o3 += A32;
	    }
          }
        }
      }
    }
    else { // We have displacements in at least one direction
      float oo1,oo2,oo3;
      if (derivdir.size()) { // If we need to calculate derivatives in at least one direction
        for (int z=0; z<out.zsize(); z++) { 
          for (int y=0; y<out.ysize(); y++) {
            for (int x=0; x<out.xsize(); x++) {
              if (useiT) {
                o1 = T11*x + T12*y + T13*z + T14; 
                o2 = T21*x + T22*y + T23*z + T24; 
                o3 = T31*x + T32*y + T33*z + T34;
                if (xd<0) oo1 = A11*o1 + A12*o2 + A13*o3 + A14;
                else oo1 = A11*o1 + A12*o2 + A13*o3 + A14 + d[xd].interpolate(o1,o2,o3);
                if (yd<0) oo2 = A21*o1 + A22*o2 + A23*o3 + A24;
                else oo2 = A21*o1 + A22*o2 + A23*o3 + A24 + d[yd].interpolate(o1,o2,o3); 
                if (zd<0) oo3 = A31*o1 + A32*o2 + A33*o3 + A34;
                else oo3 = A31*o1 + A32*o2 + A33*o3 + A34 + d[zd].interpolate(o1,o2,o3);
                if (valid) valid->operator()(x,y,z) = (d.valid(o1,o2,o3)) ? 1 : 0;   // Label as outside FOV if no info on warp
	      }
	      else {
		o1 = A11*x + A12*y + A13*z + A14;
		o2 = A21*x + A22*y + A23*z + A24;
		o3 = A31*x + A32*y + A33*z + A34;
                if (xd<0) {oo1=o1;} else {oo1=o1+d(x,y,z,xd);}
                if (yd<0) {oo2=o2;} else {oo2=o2+d(x,y,z,yd);}
                if (zd<0) {oo3=o3;} else {oo3=o3+d(x,y,z,zd);}
                if (valid) valid->operator()(x,y,z) = 1;  // So far, so good
	      }
              o1 = M11*oo1 + M12*oo2 + M13*oo3 + M14;
              o2 = M21*oo1 + M22*oo2 + M23*oo3 + M24;
              o3 = M31*oo1 + M32*oo2 + M33*oo3 + M34;
              if (derivdir.size() == 1) { // If we want a single partial
                float tmp;
	        out(x,y,z) = ((T) f.interp1partial(o1,o2,o3,derivdir[0],&tmp));
                deriv(x,y,z,0) = ((T) tmp);
	      }
              else { // More than one derivative
		float tmp1,tmp2,tmp3;
		out(x,y,z) = ((T) f.interp3partial(o1,o2,o3,&tmp1,&tmp2,&tmp3));
                if (!(xp<0)) {deriv(x,y,z,xp)=((T) tmp1);}
                if (!(yp<0)) {deriv(x,y,z,yp)=((T) tmp2);}
                if (!(zp<0)) {deriv(x,y,z,zp)=((T) tmp3);}
	      }
              if (valid) valid->operator()(x,y,z) &= (f.valid(o1,o2,o3)) ? 1 : 0; // Kosher only if valid in both d and s
	    }
          }
        }
      }
      else { // If we don't need derivatives
        for (int z=0; z<out.zsize(); z++) { 
          for (int y=0; y<out.ysize(); y++) {
            for (int x=0; x<out.xsize(); x++) {
              if (useiT) {
                o1 = T11*x + T12*y + T13*z + T14; 
                o2 = T21*x + T22*y + T23*z + T24; 
                o3 = T31*x + T32*y + T33*z + T34;
                if (xd<0) oo1 = A11*o1 + A12*o2 + A13*o3 + A14;
                else oo1 = A11*o1 + A12*o2 + A13*o3 + A14 + d[xd].interpolate(o1,o2,o3);
                if (yd<0) oo2 = A21*o1 + A22*o2 + A23*o3 + A24;
                else oo2 = A21*o1 + A22*o2 + A23*o3 + A24 + d[yd].interpolate(o1,o2,o3); 
                if (zd<0) oo3 = A31*o1 + A32*o2 + A33*o3 + A34;
                else oo3 = A31*o1 + A32*o2 + A33*o3 + A34 + d[zd].interpolate(o1,o2,o3);
                if (valid) valid->operator()(x,y,z) = (d.valid(o1,o2,o3)) ? 1 : 0;   // Label as outside FOV if no info on warp
	      }
	      else {
		o1 = A11*x + A12*y + A13*z + A14;
		o2 = A21*x + A22*y + A23*z + A24;
		o3 = A31*x + A32*y + A33*z + A34;
                if (xd<0) {oo1=o1;} else {oo1=o1+d(x,y,z,xd);}
                if (yd<0) {oo2=o2;} else {oo2=o2+d(x,y,z,yd);}
                if (zd<0) {oo3=o3;} else {oo3=o3+d(x,y,z,zd);}
                if (valid) valid->operator()(x,y,z) = 1;  // So far, so good
	      }
              o1 = M11*oo1 + M12*oo2 + M13*oo3 + M14;
              o2 = M21*oo1 + M22*oo2 + M23*oo3 + M24;
              o3 = M31*oo1 + M32*oo2 + M33*oo3 + M34;
	      out(x,y,z) = ((T) f.interpolate(o1,o2,o3));
              if (valid) valid->operator()(x,y,z) &= (f.valid(o1,o2,o3)) ? 1 : 0; // Kosher only if valid in both d and s
	    }
          }
        }
      }
    }

    //
    // Set the sform and qform appropriately.
    // 1. If the outvol has it's codes set, then leave it.
    // 2. If outvol doesn't have the codes set AND there
    //    is a warpfield AND the warpfield has the codes
    //    set then copy codes warpfield->outvol.
    // 3. If the outvol doesn't have its codes set AND 
    //    there is NO warpfield AND invol has its codes
    //    set, then set the q and sform to the transformed
    //    version of those in invol. Set codes to ALIGNED_ANAT.
    //
    NEWMAT::Matrix nmat;
    if ( (out.sform_code()==NIFTI_XFORM_UNKNOWN) &&       // qform is set in outvol
         (out.qform_code()!=NIFTI_XFORM_UNKNOWN) ) {
      out.set_sform(out.qform_code(), out.qform_mat());   // Copy qform->sform
    }
    else if ( (out.qform_code()==NIFTI_XFORM_UNKNOWN) &&  // sform is set in outvol
	      (out.sform_code()!=NIFTI_XFORM_UNKNOWN) ) {
      out.set_qform(out.sform_code(), out.sform_mat());   // Copy sform->qform
    }
    else if ( (out.qform_code()==NIFTI_XFORM_UNKNOWN) &&  // Neither is set in outvol
	      (out.sform_code()==NIFTI_XFORM_UNKNOWN) ) {
      if (defdir.size()) {                                // If there is a warp-field
        if (d.sform_code()!=NIFTI_XFORM_UNKNOWN) {        // If sform of warp-field is known
          out.set_sform(d.sform_code(),d.sform_mat());
          out.set_qform(d.sform_code(),d.sform_mat());
	}
        else if (d.qform_code()!=NIFTI_XFORM_UNKNOWN) {
          out.set_sform(d.qform_code(),d.qform_mat());
          out.set_qform(d.qform_code(),d.qform_mat());
        }
      }
      else { // I there is no warp-field, i.e. an affine transform.
        if (f.sform_code()!=NIFTI_XFORM_UNKNOWN) {
          if (TT) iA = iA*TT->i();
          if (M) iA = M->i()*iA;
          iA = f.sform_mat()*iA;
          out.set_sform(NIFTI_XFORM_ALIGNED_ANAT,iA);
          out.set_qform(NIFTI_XFORM_ALIGNED_ANAT,iA);
        }
        else if (f.sform_code()!=NIFTI_XFORM_UNKNOWN) {
          if (TT) iA = iA*TT->i();
          if (M) iA = M->i()*iA;
          iA = f.qform_mat()*iA;
          out.set_sform(NIFTI_XFORM_ALIGNED_ANAT,iA);
          out.set_qform(NIFTI_XFORM_ALIGNED_ANAT,iA);
	}
      }
    }
    
    // restore settings and return
    f.setextrapolationmethod(oldex);
    if (d.tsize()) {
      d.setextrapolationmethod(d_oldex);
      d.setextrapolationvalidity(d_old_epvalidity[0],d_old_epvalidity[1],d_old_epvalidity[2]); 
    }
   // All done!  
  }

/////////////////////////////////////////////////////////////////////
//
// The following two routines are slightly simplified interfaces
// such that users should not have to pass in zero-pointers when
// they do not want to have an "infov" mask outpuy, or when there
// are only a template and an inout volume such that we have no
// need for M or TT.
//
/////////////////////////////////////////////////////////////////////


  template <class T>
  void raw_general_transform(// Input
                             const volume<T>&         vin,        // Input volume
                             const NEWMAT::Matrix&    A,          // 4x4 affine transformation matrix
                             const volume4D<float>&   d,          // Displacement fields
                             const vector<int>&       defdir,     // Directions of displacements. 
                             const vector<int>&       derivdir,   // Directions of derivatives
                             // Output
                             volume<T>&               vout,       // Output volume
                             volume4D<T>&             deriv)      // Partial derivatives
  {
    raw_general_transform(vin,A,d,defdir,derivdir,0,0,vout,deriv,0);
  }

  template <class T>
  void raw_general_transform(// Input
                             const volume<T>&         vin,        // Input volume
                             const NEWMAT::Matrix&    A,          // 4x4 affine transformation matrix
                             const volume4D<float>&   d,          // Displacement fields
                             const vector<int>&       defdir,     // Directions of displacements. 
                             const vector<int>&       derivdir,   // Directions of derivatives
                             // Output
                             volume<T>&               vout,       // Output volume
                             volume4D<T>&             deriv,      // Partial derivative directions
                             volume<char>&            invol)      // Mask indicating what voxels fell inside original volume
  {
    raw_general_transform(vin,A,d,defdir,derivdir,0,0,vout,deriv,&invol);    
  }

  // These routines supply a convenient interface for applywarp.

  template<class T>
  void apply_warp(// Input
                  const volume<T>&        vin,         // Input volume
                  const NEWMAT::Matrix&   A,           // 4x4 affine transform
                  const volume4D<float>   d,           // Displacement fields
                  const NEWMAT::Matrix&   TT,
                  const NEWMAT::Matrix&   M,
                  // Output
                  volume<T>&              vout)        // Resampled output volume
  {
    std::vector<int>  defdir(3);
    for (int i=0; i<3; i++) defdir[i] = i;
    std::vector<int>          derivdir;
    volume4D<T>               deriv;
    const NEWMAT::Matrix      *Tptr = NULL;
    const NEWMAT::Matrix      *Mptr = NULL;

    if ((TT-IdentityMatrix(4)).MaximumAbsoluteValue() > 1e-6) Tptr = &TT;
    if ((M-IdentityMatrix(4)).MaximumAbsoluteValue() > 1e-6) Mptr = &M;

    raw_general_transform(vin,A,d,defdir,derivdir,Tptr,Mptr,vout,deriv,NULL);
  }

  template<class T>
  void apply_warp(// Input
                  const volume<T>&        vin,         // Input volume
                  const NEWMAT::Matrix&   A,           // 4x4 affine transform
                  const volume4D<float>   d,           // Displacement fields
                  const NEWMAT::Matrix&   TT,
                  const NEWMAT::Matrix&   M,
                  // Output
                  volume<T>&              vout,        // Resampled output volume
		  volume<char>&           mask)        // Set when inside original volume
  {
    std::vector<int>  defdir(3);
    for (int i=0; i<3; i++) defdir[i] = i;
    std::vector<int>          derivdir;
    volume4D<T>               deriv;
    const NEWMAT::Matrix      *Tptr = NULL;
    const NEWMAT::Matrix      *Mptr = NULL;
    mask.reinitialize(vout.xsize(),vout.ysize(),vout.zsize());  // Just to be certain
    copybasicproperties(vout,mask);

    if ((TT-IdentityMatrix(4)).MaximumAbsoluteValue() > 1e-6) Tptr = &TT;
    if ((M-IdentityMatrix(4)).MaximumAbsoluteValue() > 1e-6) Mptr = &M;

    raw_general_transform(vin,A,d,defdir,derivdir,Tptr,Mptr,vout,deriv,&mask);
  }

/////////////////////////////////////////////////////////////////////
//
// The following three routines are interfaces to mimick the old
// functions apply_warp and raw_apply_warp. These are used mainly 
// for resampling of images that are typically related to the
// input to e.g. fnirt or fugue through some rigid matrix M. 
// Examples of M would be a rigid mapping between a subjects 
// functional volumes and his/her structural or between a subjects
// functional volumes and his/her field map. 
//
/////////////////////////////////////////////////////////////////////

  template <class T>
  int apply_warp(const volume<T>&        invol,
                 volume<T>&              outvol,
                 const volume4D<float>&  warpvol)
  {
    NEWMAT::IdentityMatrix eye(4);
    return(apply_warp(invol,outvol,warpvol,eye,eye));
  }

  template <class T>
  int apply_warp(const volume<T>&                invol, 
                 volume<T>&                      outvol,
	         const volume4D<float>&          warpvol, 
	         const NEWMAT::Matrix&           premat, 
                 const NEWMAT::Matrix&           postmat)
  {
  // set the desired extrapolation settings
  extrapolation oldin = invol.getextrapolationmethod();
  extrapolation oldwarp = warpvol.getextrapolationmethod();
  warpvol.setextrapolationmethod(extraslice);
  invol.setextrapolationmethod(extraslice);
  float oldpad = invol.getpadvalue();
  invol.setpadvalue(invol.backgroundval());

  int retval = raw_apply_warp(invol,outvol,warpvol,premat,postmat);

  // restore extrapolation settings
  warpvol.setextrapolationmethod(oldwarp);
  invol.setextrapolationmethod(oldin);
  invol.setpadvalue(oldpad);
  
  return retval;
  }

  template <class T>
  int raw_apply_warp(const volume<T>&                invol, 
                     volume<T>&                      outvol,
		     const volume4D<float>&          warpvol, 
		     const NEWMAT::Matrix&           premat, 
                     const NEWMAT::Matrix&           postmat)
  {
    NEWMAT::IdentityMatrix   A(4);
    vector<int>              defdir(3);
    vector<int>              derivdir;
    volume4D<T>              deriv;

    defdir[0] = 0; defdir[1] = 1; defdir[2] = 2;

    raw_general_transform(invol,A,warpvol,defdir,derivdir,&postmat,&premat,outvol,deriv,NULL);

    return(0);
  }

  
  // The following handful of routines are simplified interfaces to 
  // raw_general_transform that may be convenient for certain
  // specific applications.

  template <class T>
  void affine_transform(// Input
                        const volume<T>&         vin,
                        const NEWMAT::Matrix&    aff,
                        // Output
                        volume<T>&               vout)
  {
    volume4D<float>  pdf;
    vector<int>      pdefdir;
    volume4D<float>  deriv;
    vector<int>      pderivdir;

    raw_general_transform(vin,aff,pdf,pdefdir,pderivdir,vout,deriv);
  }
  template <class T>
  void affine_transform(// Input
                        const volume<T>&         vin,
                        const NEWMAT::Matrix&    aff,
                        // Output
                        volume<T>&               vout,
                        volume<char>&            invol)
  {
    volume4D<float>  pdf;
    vector<int>      pdefdir;
    volume4D<float>  deriv;
    vector<int>      pderivdir;

    raw_general_transform(vin,aff,pdf,pdefdir,pderivdir,vout,deriv,invol);
  }
  template <class T>
  void affine_transform_3partial(// Input
                                 const volume<T>&             vin,
                                 const NEWMAT::Matrix&        aff,
                                 // Output
                                 volume<T>&                   vout,
                                 volume4D<T>&                 deriv)
  {
    volume4D<float>  pdf;
    vector<int>      pdefdir;
    vector<int>      pderivdir(3);
    for (int i=0; i<3; i++) {pderivdir[i] = i;}
    raw_general_transform(vin,aff,pdf,pdefdir,pderivdir,vout,deriv);
  }
  
  template <class T>
  void affine_transform_3partial(// Input
                                 const volume<T>&       vin,
                                 const NEWMAT::Matrix&  aff,
                                 // Output
                                 volume<T>&             vout,
                                 volume4D<T>&           deriv,
                                 volume<char>&          invol)
  {
    volume4D<float>  pdf;
    vector<int>      pdefdir;
    vector<int>      pderivdir(3);
    for (int i=0; i<3; i++) {pderivdir[i] = i;}
    raw_general_transform(vin,aff,pdf,pdefdir,pderivdir,vout,deriv,invol);
  }
  
  template <class T>
  void displacement_transform_1D(// Input
                                 const volume<T>&       vin,
                                 const NEWMAT::Matrix&  aff,
                                 const volume<float>&   df,
                                 int                    defdir,
                                 // Output
                                 volume<T>&             vout)
  {
    volume4D<float>  pdf;
    vector<int>      pdefdir(1,defdir);
    vector<int>      pderivdir;
    volume4D<T>      pderiv;

    pdf.addvolume(df);
    raw_general_transform(vin,aff,pdf,pdefdir,pderivdir,vout,pderiv);
  }

  template <class T>
  void displacement_transform_1D(// Input
                                 const volume<T>&       vin,
                                 const NEWMAT::Matrix&  aff,
                                 const volume<float>&   df,
                                 int                    defdir,
                                 // Output
                                 volume<T>&             vout,
                                 volume<char>&          invol)
  {
    volume4D<float>  pdf;
    vector<int>      pdefdir(1,defdir);
    vector<int>      pderivdir;
    volume4D<T>      pderiv;

    pdf.addvolume(df);
    raw_general_transform(vin,aff,pdf,pdefdir,pderivdir,vout,pderiv,invol);
  }
  
  template <class T>
  void displacement_transform_1D_3partial(// Input
                                          const volume<T>&       vin,
                                          const NEWMAT::Matrix&  aff,
                                          const volume<float>&   df,
                                          int                    dir,
                                          // Output
                                          volume<T>&             vout,
                                          volume4D<T>&           deriv)
  {
    volume4D<float>      pdf;
    vector<int>          pdefdir(1,dir);
    vector<int>          pderivdir(3);
    
    for (int i=0; i<3; i++) {pderivdir[i] = i;}
    pdf.addvolume(df);
    raw_general_transform(vin,aff,pdf,pdefdir,pderivdir,vout,deriv);    
  }
			
  template <class T>
  void displacement_transform_1D_3partial(// Input
                                          const volume<T>&       vin,
                                          const NEWMAT::Matrix&  aff,
                                          const volume<float>&   df,
                                          int                    dir,
                                          // Output
                                          volume<T>&             vout,
                                          volume4D<T>&           deriv,
                                          volume<char>&          invol)
  {
    volume4D<float>      pdf;
    vector<int>          pdefdir(1,dir);
    vector<int>          pderivdir(3);
    
    for (int i=0; i<3; i++) {pderivdir[i] = i;}
    pdf.addvolume(df);
    raw_general_transform(vin,aff,pdf,pdefdir,pderivdir,vout,deriv,invol);    
  }
			
  template <class T>	       
  void general_transform(// Input
                         const volume<T>&         vin,
                         const NEWMAT::Matrix&    aff,
                         const volume4D<float>&   df,
                         // Output
                         volume<T>&               vout)
  {
    vector<int>       pdefdir(3);
    vector<int>       pderivdir;
    volume4D<T>       pderiv;

    for (int i=0; i<3; i++) {pdefdir[i] = i;}
    raw_general_transform(vin,aff,df,pdefdir,pderivdir,vout,pderiv);
  }

  template <class T>	       
  void general_transform(// Input
                         const volume<T>&         vin,
                         const NEWMAT::Matrix&    aff,
                         const volume4D<float>&   df,
                         // Output
                         volume<T>&               vout,
                         volume<char>&            invol)
  {
    vector<int>       pdefdir(3);
    vector<int>       pderivdir;
    volume4D<T>       pderiv;

    for (int i=0; i<3; i++) {pdefdir[i] = i;}
    raw_general_transform(vin,aff,df,pdefdir,pderivdir,vout,pderiv,invol);
  }

  template <class T>
  void general_transform_3partial(// Input
                                  const volume<T>&         vin,
                                  const NEWMAT::Matrix&    aff,
                                  const volume4D<float>&   df,
                                  // Output
                                  volume<T>&               vout,
                                  volume4D<T>&             deriv)
  {
    vector<int>   dir(3);
    for (int i=0; i<3; i++) {dir[i] = i;}

    raw_general_transform(vin,aff,df,dir,dir,vout,deriv);
  }
                                
  template <class T>
  void general_transform_3partial(// Input
                                  const volume<T>&         vin,
                                  const NEWMAT::Matrix&    aff,
                                  const volume4D<float>&   df,
                                  // Output
                                  volume<T>&               vout,
                                  volume4D<T>&             deriv,
                                  volume<char>&            invol)
  {
    vector<int>   dir(3);
    for (int i=0; i<3; i++) {dir[i] = i;}

    raw_general_transform(vin,aff,df,dir,dir,vout,deriv,invol);
  }
 
  // 
  // This routine mimics some of the code in raw_general_transform.
  // Ideally they should be re-written to allow for code re-use
  // but this solution is the least "invasive" in the sense of
  // risking to introduce any new bugs.
  // It will return a 3D displacement field in mm.
  //
  template <class T> 
  void get_displacement_fields(// Input
                                  const volume<T>&         vin,
                                  const NEWMAT::Matrix&    aff,
                                  const volume4D<float>&   dfin,
                                  // Output
                                  volume4D<float>&         dfout)
  {
    if (dfin.tsize() != 3) imthrow("NEWIMAGE::get_displacement_field: dfin.tsize() must be 3.",11);
    if (!NEWIMAGE::samesize(vin,dfin[0])) imthrow("NEWIMAGE::get_displacement_field: mismatch between vin and dfin.",11);
    if (!NEWIMAGE::samesize(dfin,dfout)) imthrow("NEWIMAGE::get_displacement_field: mismatch between dfin and dfout.",11);
    NEWMAT::Matrix v2w = vin.sampling_mat();  // voxel->world      
    NEWMAT::Matrix iA = aff.i();              // world->transformed_world
    NEWMAT::ColumnVector vox(4); vox(4) = 1.0;
    NEWMAT::ColumnVector mm(4);
    NEWMAT::ColumnVector tmm(4);
    for (int k=0; k<vin.zsize(); k++) {
      vox(3) = k;
      for (int j=0; j<vin.ysize(); j++) {
	vox(2) = j;
	for (int i=0; i<vin.xsize(); i++) {
	  vox(1) = i;
	  mm = v2w*vox;
	  tmm = iA*mm;
          dfout(i,j,k,0) = tmm(1) + dfin(i,j,k,0) - mm(1);
          dfout(i,j,k,1) = tmm(2) + dfin(i,j,k,1) - mm(2);
          dfout(i,j,k,2) = tmm(3) + dfin(i,j,k,2) - mm(3);
	}
      }
    }
  }                             
}

#ifdef I_DEFINED_ET
#undef I_DEFINED_ET
#undef EXPOSE_TREACHEROUS   // Avoid exporting dodgy routines
#endif

#endif
