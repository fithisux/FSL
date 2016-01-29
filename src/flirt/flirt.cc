/*  FLIRT - FMRIB's Linear Image Registration Tool

    flirt.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

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
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <time.h>
#include <vector>
#include <algorithm>
#define WANT_STREAM
#define WANT_MATH

#ifndef EXPOSE_TREACHEROUS
#define EXPOSE_TREACHEROUS
#endif

#include "newmatap.h"
#include "newmatio.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/optimise.h"
#include "newimage/costfns.h"
#include "newimage/newimageall.h"
#include "defaultschedule.h"
#include "globaloptions.h"

#ifndef NO_NAMESPACE
 using namespace MISCMATHS;
 using namespace NEWMAT;
 using namespace NEWIMAGE;
#endif

// Put current version number here:
const string version = "6.0";

////////////////////////////////////////////////////////////////////////////

// GLOBAL VOLUMES FOR WEIGHTING INFO

volume<float> global_refweight, global_testweight, global_refweight1;
volume<float> global_refweight2, global_refweight4, global_refweight8;
volume<float> global_seg, global_init_testvol, global_init_testweight;
volume<float> global_fmap, global_fmap_mask;
Matrix global_coords, global_norms;
bool global_scale1OK=true, read_testvol=false;
float global_sampling=1.0f;

////////////////////////////////////////////////////////////////////////////

void print_vector(float x, float y, float z)
{
  cerr << "(" << x << "," << y << "," << z << ")"; 
}


//------------------------------------------------------------------------//

void setupsinc(const volume<float>& invol)
{
  // the following full-width is in voxels
  int w = MISCMATHS::round(globaloptions::get().sincwidth);
  if (globaloptions::get().sincwindow==Hanning) {
    invol.definesincinterpolation("hanning",w);
  } else if (globaloptions::get().sincwindow==Blackman) {
    invol.definesincinterpolation("blackman",w);
  } else if (globaloptions::get().sincwindow==Rect) {
    invol.definesincinterpolation("rectangular",w);
  }
}


void set_basescale(const string& filenameA, const string& filenameB)
{
  if (!(globaloptions::get().force_basescale)) { 
    // only try to automatically determine if it was not requested by the user
    volume<float> volA, volB;
    read_volume_hdr_only(volA,filenameA);
    read_volume_hdr_only(volB,filenameB);
    float maxdimA = Max(Max(volA.xdim(),volA.ydim()),volA.zdim());
    float maxdimB = Max(Max(volB.xdim(),volB.ydim()),volB.zdim());
    if (Max(maxdimA,maxdimB)>12) {  // over 150% of largest (8mm) scale
      globaloptions::get().basescale = Max(maxdimA,maxdimB)/8;
    }
    float mindimA = Min(Min(volA.xdim(),volA.ydim()),volA.zdim());
    float mindimB = Min(Min(volB.xdim(),volB.ydim()),volB.zdim());
    if (Min(mindimA,mindimB)<0.75) { // between 0.5 and 1 mm 
      globaloptions::get().basescale = Min(mindimA,mindimB);
    }
  }
}


int FLIRT_read_volume4D(volume4D<float>& target, const string& filename)
{
  // make voxels bigger if basescale is smaller than 1.0 (and vice versa)
  int retval = read_volume4D(target,filename);  // as radiological
  // if basescale != 1.0
  if (fabs(globaloptions::get().basescale - 1.0)>1e-5) {
    target.setxdim(target.xdim() / globaloptions::get().basescale);
    target.setydim(target.ydim() / globaloptions::get().basescale);
    target.setzdim(target.zdim() / globaloptions::get().basescale);
  }
  return retval;
}

int FLIRT_read_volume(volume<float>& target, const string& filename)
{
  // make voxels bigger if basescale is smaller than 1.0 (and vice versa)
  int retval = read_volume(target,filename);  // as radiological
  // if basescale != 1.0
  if (fabs(globaloptions::get().basescale - 1.0)>1e-5) {
    target.setxdim(target.xdim() / globaloptions::get().basescale);
    target.setydim(target.ydim() / globaloptions::get().basescale);
    target.setzdim(target.zdim() / globaloptions::get().basescale);
  }
  return retval;
}


//------------------------------------------------------------------------//


int vector2affine(const ColumnVector& params, int n, const ColumnVector& centre,
		  Matrix& aff)
{
  if (n<=0) return 0;
  // order of parameters is 3 rotation + 3 translation + 3 scales + 3 skews
  // angles are in radians

  switch (globaloptions::get().anglerep) 
    {
    case Euler:
      compose_aff(params,n,centre,aff,construct_rotmat_euler);
      break;
    case Quaternion:
      compose_aff(params,n,centre,aff,construct_rotmat_quat);
      break;
    default:
      cerr << "Invalid Rotation Representation" << endl;
      return -1;
    }
  return 0;
}  


int vector2affine(const ColumnVector& params, int n, Matrix& aff)
{
  return vector2affine(params,n,globaloptions::get().impair->testvol.cog("scaled_mm"),aff);
}


int vector2affine(const float params[], int n, Matrix& aff)
{
  ColumnVector vparams(12);
  for (int i=1; i<=n; i++) {
    vparams(i) = params[i];
  }
  return vector2affine(vparams,n,aff);
}


int affmat2vector(const Matrix& aff, int n, const ColumnVector& centre,
		  ColumnVector& params)
{
  switch (globaloptions::get().anglerep) 
    {
    case Euler:
      decompose_aff(params,aff,centre,rotmat2euler);
      break;
    case Quaternion:
      decompose_aff(params,aff,centre,rotmat2quat);
      break;
    default:
      cerr << "Invalid Rotation Representation" << endl;
    }
  return 0;
}


int affmat2vector(const Matrix& aff, int n, ColumnVector& params)
{
  return affmat2vector(aff,n,globaloptions::get().impair->testvol.cog("scaled_mm"),params);
}


void set_param_basis(Matrix &parambasis, int no_params)
{
  parambasis = 0.0;
  for (int i=1; i<=no_params; i++) {
    parambasis(i,i)=1.0;
  }
}


ColumnVector default_nonlin_params(void)
{
  int pe_dir = globaloptions::get().pe_dir;
  int N=0;
  if (abs(pe_dir)==1) N=globaloptions::get().impair->testvol.xsize();
  if (abs(pe_dir)==2) N=globaloptions::get().impair->testvol.ysize();
  if (abs(pe_dir)==3) N=globaloptions::get().impair->testvol.zsize();
  float fmapscaling = globaloptions::get().echo_spacing * N / (2.0*M_PI);
  if (pe_dir<0) fmapscaling *= -1;
  ColumnVector nonlin_params(1);
  nonlin_params=fmapscaling;
  return nonlin_params;
}


// cost function interfaces

int setcostfntype(Costfn* imagepair, costfns ctype) {
  imagepair->set_costfn(ctype);
  if ((imagepair->get_costfn()==BBR) && (!imagepair->is_bbr_set())) {
    if (globaloptions::get().debug) {
      cerr << "Setting bbr seg 1" << endl;
      cerr << "Cost is set to " << imagepair->get_costfn() << endl;
      cerr << "Cost == BBR is " << (imagepair->get_costfn() == BBR) << endl;
      cerr << "Result (pre) of is_bbr_set is " << imagepair->is_bbr_set() << endl;
    }
    if (globaloptions::get().usecoords) {
      imagepair->set_bbr_coords(global_coords,global_norms);
    } else {
      imagepair->set_bbr_seg(global_seg);  // only run BBR at one scale so not wasteful
    }
    imagepair->set_bbr_fmap(global_fmap,global_fmap_mask,globaloptions::get().pe_dir);
    imagepair->set_bbr_type(globaloptions::get().bbr_type);
    imagepair->set_bbr_slope(globaloptions::get().bbr_slope);
    if (globaloptions::get().debug) {
      cerr << "Result (post) of is_bbr_set is " << imagepair->is_bbr_set() << endl;
    }
  }
  return 0;
}

int setcostfntype(costfns ctype) {
  if (globaloptions::get().impair) {  // only do this if impair is set
    return setcostfntype(globaloptions::get().impair, ctype);
  }
  return -1;
}

int setcostfntype(const string& cname) {
  return setcostfntype(costfn_type(cname));
}



int setup_costfn(Costfn* imagepair, costfns curcostfn, int no_bins, float smoothsize, float fuzzyfrac)
{
  setcostfntype(imagepair,curcostfn);
  imagepair->set_no_bins(no_bins);
  imagepair->smoothsize = smoothsize;
  imagepair->fuzzyfrac = fuzzyfrac;
  return 0;
}


float costfn(const Matrix& uninitaffmat, const ColumnVector& nonlin_params)
{
  Tracer tr("costfn");
  Matrix affmat = uninitaffmat * globaloptions::get().initmat;  // apply initial matrix
  setcostfntype(globaloptions::get().currentcostfn);
  float retval = 0.0;
  retval = globaloptions::get().impair->cost(affmat,nonlin_params);
  return retval;
}


float costfn(const Matrix& uninitaffmat)
{
  Tracer tr("costfn");
  float retval = 0.0;
  if ((globaloptions::get().currentcostfn==BBR) && (globaloptions::get().pe_dir!=0)) {
    // call the non-linear version of costfn, which will apply the initmat there
    retval = costfn(uninitaffmat,default_nonlin_params());
  } else {
    Matrix affmat = uninitaffmat * globaloptions::get().initmat;  // apply initial matrix
    setcostfntype(globaloptions::get().currentcostfn);
    retval = globaloptions::get().impair->cost(affmat);
  }
  return retval;
}


float costfn(const ColumnVector& params)
{
  Tracer tr("costfn");
  Matrix affmat(4,4);
  vector2affine(params,globaloptions::get().no_params,affmat);
  float retval;
  int pe_dir=globaloptions::get().pe_dir;
  if ((globaloptions::get().impair->get_costfn()==BBR) && (pe_dir!=0)) {
    retval = costfn(affmat,default_nonlin_params());
  } else {
    retval = costfn(affmat);
  }
  if (globaloptions::get().verbose>=5) {
    cout << globaloptions::get().impair->count() << " : ";
    cout << retval << " :: ";
    for (int i=1; i<=globaloptions::get().no_params; i++) 
      { cout << params(i) << " "; }
    cout << endl;
  }
  return retval;
}
  

//----------------------------------------------------------------------//

void affine_and_fmap_transform(const volume<float>& testvol, const volume<float>& refvol,
			       volume<float>& outputvol, const Matrix& finalmat, const ColumnVector& nonlin_params)
{
  if (globaloptions::get().debug) { cerr << "Pre bbr_resamp" << endl; }
  globaloptions::get().impair->bbr_resamp(finalmat,nonlin_params,outputvol);
  if (globaloptions::get().debug) { cerr << "Post bbr_resamp" << endl; }
}


void final_transform(const volume<float>& testvol, const volume<float>& refvol,
		     const Matrix& finalmat, volume<float>& outputvol) 
{
  if (globaloptions::get().interpmethod == NearestNeighbour) {
    testvol.setinterpolationmethod(nearestneighbour);
  } else if (globaloptions::get().interpmethod == NEWIMAGE::Sinc) {
    setupsinc(testvol);
    testvol.setinterpolationmethod(sinc);
  } else if (globaloptions::get().interpmethod == NEWIMAGE::Spline) {
    testvol.setinterpolationmethod(spline);
  } else {
    testvol.setinterpolationmethod(trilinear);
  }

  if (globaloptions::get().forcebackgnd) {
    testvol.setpadvalue(globaloptions::get().backgndval);
    testvol.setextrapolationmethod(constpad);
  } else {
    testvol.setpadvalue(testvol.backgroundval());
    testvol.setextrapolationmethod(extraslice);
  }


  float paddingsize = globaloptions::get().paddingsize;
  if (globaloptions::get().mode2D) {
    paddingsize = Max(1.0,paddingsize);
  }
  if (globaloptions::get().pe_dir==0) {  // test to see if fieldmap is being used
    affine_transform(testvol,outputvol,finalmat,paddingsize,false);
  } else {
    // Only setup costfn if it isn't already done (normally first time around in a 4D)
    if (globaloptions::get().debug) { cerr << "Start affine_and_fmap_transform" << endl; }
    if (globaloptions::get().impair==NULL) {
      if (globaloptions::get().debug) { cerr << "Pre new costfn" << endl; }
      globaloptions::get().impair = new Costfn(refvol,testvol);
      if (globaloptions::get().debug) { cerr << "Pre setup_costfn" << endl; }
      // the following sets up all BBR related stuff like fieldmaps, etc.  (will also recalculate wm edge points)
      setup_costfn(globaloptions::get().impair, BBR,
		   globaloptions::get().no_bins,
		   globaloptions::get().smoothsize,globaloptions::get().fuzzyfrac);  
    }
    affine_and_fmap_transform(testvol,refvol,outputvol,finalmat,default_nonlin_params());
  }
}

template <class T>
int safe_save_volume(const volume<T>& source, const string& filename)
{
  if (!globaloptions::get().nosave) {
    const_cast< volume<T>& >(source).setDisplayMaximumMinimum(0,0);
    save_volume_dtype(source,filename,globaloptions::get().datatype);
  }
  return 0;
}

void save_matrix_data(const Matrix& matresult)
{
  Matrix outfmat = matresult;
  write_ascii_matrix(outfmat,globaloptions::get().outputmatascii); 
}

void save_matrix_data(const Matrix& matresult, const volume<float>& initvol, 
		      const volume<float>& finalvol)
{
  save_matrix_data(matresult);
}

//----------------------------------------------------------------------------//


//////////////////////////////////////////////////////////////////////////

// OPTIMISATION SUPPORT


float estimate_scaling(const volume<float>& vol) {
  Tracer tr("estimate_scaling");
  return Min(Min(vol.xdim(),vol.ydim()),vol.zdim());
}

float estimate_scaling() {
  Tracer tr("estimate_scaling");
  return estimate_scaling(globaloptions::get().impair->refvol);
}

void set_param_tols(ColumnVector &param_tol, int no_params)
{
      // Tolerances are: 0.57 degrees (0.005 radians), 0.2mm translation
      //    0.005 scale and 0.001 skew
//    float diagonal[12]={0.005, 0.005, 0.005, 0.2, 0.2, 0.2, 0.002, 0.002, 0.002,
//    		      0.001, 0.001, 0.001};
  if (param_tol.Nrows()<no_params) {
    param_tol.ReSize(no_params);
  }
  for (int i=1; i<=no_params; i++) {
    //    param_tol(i)=diagonal[i-1];
    param_tol(i)=globaloptions::get().tolerance(i);
  }
  param_tol *= globaloptions::get().requestedscale; 
    // scale it up by the current scaling
}



void initialise_params(ColumnVector& params)
{
  Real paramsf[12] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0};
  params << paramsf;
}


void optimise(ColumnVector& params, int no_params, ColumnVector& param_tol, 
	      int &no_its, float *fans, 
	      float (*costfunc)(const ColumnVector &), int itmax=4)
{
  // sets up the initial parameters and calls the optimisation routine
  if ((params.MaximumAbsoluteValue() < 0.001) && (params.Nrows()>=12) )  
    initialise_params(params);
  { 
    Matrix affmattst(4,4);
    vector2affine(params,no_params,affmattst);
    if (globaloptions::get().verbose>=5) {
      cout << "Starting with : " << endl << affmattst;
    }
  }
  Matrix parambasis(no_params,no_params);
  set_param_basis(parambasis,no_params);
  float ptol[13];
  for (int i=1; i<=no_params; i++) { ptol[i] = param_tol(i); }
  
  *fans = MISCMATHS::optimise(params,no_params,param_tol,costfunc,no_its,itmax,
			      globaloptions::get().boundguess,
			      globaloptions::get().optimisationtype);
}



////////////////////////////////////////////////////////////////////////////


void params12toN(ColumnVector& params)
{
  // Convert the full 12 dof param vector to a small param vector
  Tracer tr("params12toN");
  ColumnVector nparams;
  nparams = pinv(globaloptions::get().parammask)*(params - globaloptions::get().refparams);
  params = nparams;
}


void paramsNto12(ColumnVector& params)
{
  // Convert small param vector to full 12 dof param vector
  Tracer tr("paramsNto12");
  ColumnVector param12;
  param12 = globaloptions::get().parammask*params + globaloptions::get().refparams;
  params = param12;
}



float subset_costfn(const ColumnVector& params)
{
  Tracer tr("subset_costfn");
  ColumnVector param12;
  param12 = params;
  paramsNto12(param12);
  float retval = costfn(param12);
  if (globaloptions::get().verbose>=7) {
    cout << globaloptions::get().impair->count() << " : ";
    cout << retval << " :: " << param12.t() << endl;
  }
  return retval;
}

//------------------------------------------------------------------------//




void find_cost_minima(Matrix& bestpts, const volume<float>& cost) {
  Tracer tr("find_cost_minima");
  volume<float> minv(cost);
  ColumnVector bestpt(3);
  minv = 0.0;
  bestpt = 0.0;
  int xb=cost.xsize(), yb=cost.ysize(), zb=cost.zsize(), minima=0;
  // calculate the number of neighbouring voxels *less* than the centre one
  for (int z=0; z<zb; z++) {
    for (int y=0; y<yb; y++) {
      for (int x=0; x<xb; x++) {
	for (int zo=-1; zo<=1; zo++) {
	  for (int yo=-1; yo<=1; yo++) {
	    for (int xo=-1; xo<=1; xo++) {
	      if (cost.in_bounds(x+xo,y+yo,z+zo)) {
		if (cost(x+xo,y+yo,z+zo) < cost(x,y,z)) {
		  minv(x,y,z) += 1.0;
		}
	      }
	    }
	  }
	}
	if (cost.in_bounds(x,y,z)) {
	  if (cost(x,y,z) < 
	      cost(MISCMATHS::round(bestpt(1)),
		   MISCMATHS::round(bestpt(2)),
		   MISCMATHS::round(bestpt(3))))
	    {
	      bestpt(1) = (float) x;
	      bestpt(2) = (float) y;
	      bestpt(3) = (float) z;
	    }
	  if ( cost.in_bounds(x+1,y+1,z+1) && (minv(x,y,z) < 0.5) )
	    minima++;
	}
      }
    }
  }
  int idx=1;
  if (minima<=0) {
    bestpts.ReSize(1,16);
    bestpts(1,1) = bestpt(1);
    bestpts(1,2) = bestpt(2);
    bestpts(1,3) = bestpt(3);
    return;
  }
  bestpts.ReSize(minima,16);
  for (int z=0; z<zb; z++) {
    for (int y=0; y<yb; y++) {
      for (int x=0; x<xb; x++) {
	if ( cost.in_bounds(x,y,z) && cost.in_bounds(x+1,y+1,z+1) ) {
	  if (minv(x,y,z) < 0.5) {
	    bestpts(idx,1) = (float) x; 
	    bestpts(idx,2) = (float) y; 
	    bestpts(idx,3) = (float) z;
	    idx++;
	    if (globaloptions::get().verbose>=3)
	      cout << "COST minima at : " << x << "," << y << "," << z << endl;
	  }
	}
      }
    }  
  }
}
  

void set_rot_sampling(ColumnVector& rots, float lowerbound, float upperbound) {
  Tracer tr("set_rot_sampling");
  // this function RELIES on the number of rows in the vector rots being set
  int nmax = rots.Nrows();
  if (nmax==1) { rots(1) = (upperbound + lowerbound) / 2.0;  return; }
  for (int n=1; n<=nmax; n++) {
    rots(n) = lowerbound + ((float) (n-1)) * (upperbound - lowerbound) / 
                                                        ((float) (nmax-1));
  }
}
    

void set_rot_samplings(ColumnVector& rxcoarse, ColumnVector& rycoarse,
		       ColumnVector& rzcoarse, ColumnVector& rxfine,
		       ColumnVector& ryfine, ColumnVector& rzfine) {
  Tracer tr("set_rot_samplings");
  //int coarsesize = 4, finesize = 11;
  // sets the number of rows (angle samples) for each axis
  rxcoarse.ReSize(MISCMATHS::round((globaloptions::get().searchrx(2) 
			 - globaloptions::get().searchrx(1))
			/globaloptions::get().coarsedelta)+1);
  rycoarse.ReSize(MISCMATHS::round((globaloptions::get().searchry(2) 
			 - globaloptions::get().searchry(1))
			/globaloptions::get().coarsedelta)+1);
  rzcoarse.ReSize(MISCMATHS::round((globaloptions::get().searchrz(2) 
			 - globaloptions::get().searchrz(1))
			/globaloptions::get().coarsedelta)+1);
  rxfine.ReSize(MISCMATHS::round((globaloptions::get().searchrx(2) 
		       - globaloptions::get().searchrx(1))
			/globaloptions::get().finedelta)+1);
  ryfine.ReSize(MISCMATHS::round((globaloptions::get().searchry(2) 
		       - globaloptions::get().searchry(1))
			/globaloptions::get().finedelta)+1);
  rzfine.ReSize(MISCMATHS::round((globaloptions::get().searchrz(2) 
		       - globaloptions::get().searchrz(1))
			/globaloptions::get().finedelta)+1);
  // now get the appropriate angle sample values
  set_rot_sampling(rxcoarse,globaloptions::get().searchrx(1),
		   globaloptions::get().searchrx(2));
  set_rot_sampling(rycoarse,globaloptions::get().searchry(1),
		   globaloptions::get().searchry(2));
  set_rot_sampling(rzcoarse,globaloptions::get().searchrz(1),
		   globaloptions::get().searchrz(2));
  set_rot_sampling(rxfine,globaloptions::get().searchrx(1),
		   globaloptions::get().searchrx(2));
  set_rot_sampling(ryfine,globaloptions::get().searchry(1),
		   globaloptions::get().searchry(2));
  set_rot_sampling(rzfine,globaloptions::get().searchrz(1),
		   globaloptions::get().searchrz(2));
  if (globaloptions::get().verbose>=4) {
    cout << "Coarse rotation samplings are:\n" << rxcoarse.t() << rycoarse.t() 
	 << rzcoarse.t() << " and fine rotation samplings are:\n"
	 << rxfine.t() << ryfine.t() << rzfine.t() << endl;
  }
}



void search_cost(Matrix& paramlist, volume<float>& costs, volume<float>& tx, 
		 volume<float>& ty, volume<float>& tz, volume<float>& scale) {
  Tracer tr("search_cost");
  int storedverbose = globaloptions::get().verbose;
  int storeddof = globaloptions::get().dof;
  anglereps useranglerep = globaloptions::get().anglerep;
  globaloptions::get().verbose -= 2;
  globaloptions::get().anglerep = Euler;  // a workaround hack
  globaloptions::get().currentcostfn = globaloptions::get().searchcostfn;
  globaloptions::get().dof = globaloptions::get().searchdof;

  ColumnVector coarserx, coarsery, coarserz, finerx, finery, finerz;
  set_rot_samplings(coarserx,coarsery,coarserz,finerx,finery,finerz);

  // set up the type of parameter subset (via the global mask)
  // here 3 translations and 1 (common) scaling are used
  if (globaloptions::get().dof>6) {
    globaloptions::get().parammask.ReSize(12,4);  // was 3
    globaloptions::get().parammask = 0.0;
    globaloptions::get().parammask(7,1) = 1.0;  // didn't used to exist
    globaloptions::get().parammask(8,1) = 1.0;  // didn't used to exist
    globaloptions::get().parammask(9,1) = 1.0;  // didn't used to exist
    globaloptions::get().parammask(4,2) = 1.0;
    globaloptions::get().parammask(5,3) = 1.0;
    globaloptions::get().parammask(6,4) = 1.0;
  } else {
    globaloptions::get().parammask.ReSize(12,3);
    globaloptions::get().parammask = 0.0;
    globaloptions::get().parammask(4,1) = 1.0;
    globaloptions::get().parammask(5,2) = 1.0;
    globaloptions::get().parammask(6,3) = 1.0;
  }
	
  ColumnVector param_tol, param_tol0(12), param_tol1(12), params_8(12);
  globaloptions::get().no_params = 12; // necessary for any subset_costfn call
  param_tol0 = globaloptions::get().refparams;
  params12toN(param_tol0);
  set_param_tols(param_tol1,12);
  param_tol1 = param_tol1 + globaloptions::get().refparams;
  params12toN(param_tol1);
  param_tol = param_tol1 - param_tol0;

  // search coarsely, optimising each point and storing the results
  int no_its=0;
  float fans=0.0, rx,ry,rz;
  Matrix affmat(4,4);
  ColumnVector trans(3), testv(4), testv2(4);
  tx.reinitialize(coarserx.Nrows(),coarsery.Nrows(),coarserz.Nrows());
  ty.reinitialize(coarserx.Nrows(),coarsery.Nrows(),coarserz.Nrows());
  tz.reinitialize(coarserx.Nrows(),coarsery.Nrows(),coarserz.Nrows());
  scale.reinitialize(coarserx.Nrows(),coarsery.Nrows(),coarserz.Nrows());
  // fix the reference parameter (starting estimates)
  globaloptions::get().refparams = 0.0;
  globaloptions::get().refparams(7) = 1.0;
  globaloptions::get().refparams(8) = 1.0;
  globaloptions::get().refparams(9) = 1.0;
  // set the initial translation (to align cog's)
  //  trans = refvol.cog("scaled_mm") - initmat * testvol.cog("scaled_mm")
  ColumnVector testcog(4), tcog(3);
  tcog = globaloptions::get().impair->testvol.cog("scaled_mm");
  testcog(1)=tcog(1); testcog(2)=tcog(2); testcog(3)=tcog(3); testcog(4)=1.0;
  testcog = globaloptions::get().initmat * testcog;
  trans = globaloptions::get().impair->refvol.cog("scaled_mm");
  trans(1) -= testcog(1);
  trans(2) -= testcog(2);
  trans(3) -= testcog(3);
  globaloptions::get().refparams(4) = trans(1);
  globaloptions::get().refparams(5) = trans(2);
  globaloptions::get().refparams(6) = trans(3);
  for (int ix=0; ix<coarserx.Nrows(); ix++) {
    for (int iy=0; iy<coarsery.Nrows(); iy++) {
      for (int iz=0; iz<coarserz.Nrows(); iz++) {
	rx = coarserx(ix+1);
	ry = coarsery(iy+1);
	rz = coarserz(iz+1);
	globaloptions::get().refparams(1) = rx;
	globaloptions::get().refparams(2) = ry;
	globaloptions::get().refparams(3) = rz;
	params_8 = globaloptions::get().refparams;
	if (globaloptions::get().verbose>=4) {
	  cout << "Starting with " << params_8.t();
	  cout << "  and tolerance " << param_tol.t();
	}
	params12toN(params_8);
	optimise(params_8,globaloptions::get().parammask.Ncols(),
		 param_tol,no_its,&fans,subset_costfn);
	//		 param_tol,&no_its,&fans,subset_costfn);
	paramsNto12(params_8);
	tx(ix,iy,iz) = params_8(4);
	ty(ix,iy,iz) = params_8(5);
	tz(ix,iy,iz) = params_8(6);
	scale(ix,iy,iz) = params_8(7);
	
	if (globaloptions::get().verbose>=4) {
	  cout << " dearranged: " << params_8.t();
	}
      }
      if (globaloptions::get().verbose>=2) cout << "*";
    }
  }
  if (globaloptions::get().verbose>=2) cout << endl;

  // scale = 1.0;  // for now disallow non-unity scalings
  float medianscale = scale.percentile(0.50);
  scale = medianscale;  // try a constant, median scale for all
  if (globaloptions::get().verbose>=2) 
    { cout << "Median scale = " << medianscale << endl; }

  if (globaloptions::get().verbose>=4) {
    safe_save_volume(tx,"tx");
    safe_save_volume(ty,"ty");
    safe_save_volume(tz,"tz");
    safe_save_volume(scale,"sc");
  }

  // use the optimised translation (and scale) to interpolate parameters
  //  and find the costs at many more sites (without excessive time)
  float xf, yf, zf, txv=0, tyv=0, tzv=0, scv=0;
  float factorx = ((float) coarserx.Nrows()-1) 
                      / Max((float) 1.0,((float) finerx.Nrows()-1));
  float factory = ((float) coarsery.Nrows()-1)
                      / Max((float) 1.0,((float) finery.Nrows()-1));
  float factorz = ((float) coarserz.Nrows()-1)
                      / Max((float) 1.0,((float) finerz.Nrows()-1));
  costs.reinitialize(finerx.Nrows(),finery.Nrows(),finerz.Nrows());
  for (int ix=0; ix<finerx.Nrows(); ix++) {
    for (int iy=0; iy<finery.Nrows(); iy++) {
      for (int iz=0; iz<finerz.Nrows(); iz++) {
	rx = finerx(ix+1);
	ry = finery(iy+1);
	rz = finerz(iz+1);
	globaloptions::get().refparams(1) = rx;
	globaloptions::get().refparams(2) = ry;
	globaloptions::get().refparams(3) = rz;
	xf = ((float) ix)*factorx;
	yf = ((float) iy)*factory;
	zf = ((float) iz)*factorz;
	txv = tx.interpolate(xf,yf,zf);
	tyv = ty.interpolate(xf,yf,zf);
	tzv = tz.interpolate(xf,yf,zf);
	scv = scale.interpolate(xf,yf,zf);
	if ((scv<0.5) || (scv>2.0))  scv = 1.0;
	if (globaloptions::get().dof<=6) scv = 1.0;
	globaloptions::get().refparams(4) = txv;
	globaloptions::get().refparams(5) = tyv;
	globaloptions::get().refparams(6) = tzv;
	globaloptions::get().refparams(7) = scv;
	globaloptions::get().refparams(8) = scv;
	globaloptions::get().refparams(9) = scv;
	params_8 = globaloptions::get().refparams;
	costs(ix,iy,iz) = costfn(params_8);
      }
      if (globaloptions::get().verbose>=2) cout << "*";
    }
  }
  if (globaloptions::get().verbose>=2) cout << endl;

  if (globaloptions::get().verbose>=4) {
    safe_save_volume(costs,"costs");
  }


  float costmin=costs.min(), costmax=costs.max();
  // the following assumes that costmin > 0
  float factor = 0.2;
  float costthresh = Min(costmin + factor*(costmax-costmin),
			 costs.percentile(0.20));
  // avoid the percentile giving the costmin (or less)
  if (costthresh <= costmin)  costthresh = Max(costmin*1.0001,costmin*0.9999);
  if (globaloptions::get().verbose>=4) {
    cout << "Cost threshold = " << costthresh << " and minimum is "
	 << costmin << endl;
  }

  // for all costs less than 150% of the best cost so far, optimise...
  int numsubcost = 0;
  for (int ix=0; ix<costs.xsize(); ix++) {
    for (int iy=0; iy<costs.ysize(); iy++) {
      for (int iz=0; iz<costs.zsize(); iz++) {
	if (costs(ix,iy,iz) <= costthresh) {
	  numsubcost++;
	}
      }
    }
  }
  if (numsubcost<=0) {
    cout << "WARNING: Found 0 or less sub-threshold costs" << endl;
    numsubcost = 1;
  }
  Matrix bestparams(numsubcost,13);
  int n=1;
  for (int ix=0; ix<finerx.Nrows(); ix++) {
    for (int iy=0; iy<finery.Nrows(); iy++) {
      for (int iz=0; iz<finerz.Nrows(); iz++) {
	if (costs(ix,iy,iz) < costthresh) {
	  rx = finerx(ix+1);
	  ry = finery(iy+1);
	  rz = finerz(iz+1); 
	  xf = ((float) ix)*factorx;
	  yf = ((float) iy)*factory;
	  zf = ((float) iz)*factorz;
	  txv = tx.interpolate(xf,yf,zf);
	  tyv = ty.interpolate(xf,yf,zf);
	  tzv = tz.interpolate(xf,yf,zf);
	  scv = scale.interpolate(xf,yf,zf);
	  if ((scv<0.5) || (scv>2.0))  scv = 1.0;
	  if (globaloptions::get().dof<=6) scv = 1.0;
	  params_8 = 0.0;
	  params_8(1) = rx;  params_8(2) = ry;  params_8(3) = rz;
	  params_8(4) = txv; params_8(5) = tyv; params_8(6) = tzv;
	  params_8(7) = scv; params_8(8) = scv; params_8(9) = scv;
	  globaloptions::get().refparams = params_8;
	  params12toN(params_8);
	  optimise(params_8,globaloptions::get().parammask.Ncols(),
		   param_tol,no_its,&fans,subset_costfn);
	           //param_tol,&no_its,&fans,subset_costfn);
	  paramsNto12(params_8);
	  costs(ix,iy,iz) = fans;
	  bestparams(n,1) = fans;
	  bestparams.SubMatrix(n,n,2,13) = params_8.t();
	  n++;
	  if (globaloptions::get().verbose>=3) {
	    cout << "(" << ix << "," << iy << "," << iz << ") => " << fans
		 << " with " << params_8.t();
	  }
	}
      }
    }
  }

  if (globaloptions::get().verbose>=3) {
    safe_save_volume(costs,"costs");
    cout << "Costs (1st column) are:\n" << bestparams << endl;
  }

  // find the cost minima and return these
  Matrix bestpts;
  find_cost_minima(bestpts,costs);
  paramlist.ReSize(bestpts.Nrows(),12);
  int ix,iy,iz;
  for (int n=1; n<=paramlist.Nrows(); n++) {
    ix = MISCMATHS::round(bestpts(n,1));
    iy = MISCMATHS::round(bestpts(n,2));
    iz = MISCMATHS::round(bestpts(n,3));
    if (globaloptions::get().verbose>=3) 
      cout << "Cost minima at : " << ix << "," << iy << "," << iz << endl;
    rx = finerx(ix+1);
    ry = finery(iy+1);
    rz = finerz(iz+1);
    xf = ((float) ix)*factorx;
    yf = ((float) iy)*factory;
    zf = ((float) iz)*factorz;
    txv = tx.interpolate(xf,yf,zf);
    tyv = ty.interpolate(xf,yf,zf);
    tzv = tz.interpolate(xf,yf,zf);
    scv = scale.interpolate(xf,yf,zf);
    if ((scv<0.5) || (scv>2.0))  scv = 1.0;
    if (globaloptions::get().dof<=6) scv = 1.0;
    params_8 = 0.0;
    params_8(1) = rx;  params_8(2) = ry;  params_8(3) = rz;
    params_8(4) = txv; params_8(5) = tyv; params_8(6) = tzv;
    params_8(7) = scv; params_8(8) = scv; params_8(9) = scv;
    paramlist.SubMatrix(n,n,1,12) = params_8.t();
  }

  if (globaloptions::get().verbose>=3) {
    cout << "Chosen parameters:\n" << paramlist << endl;
  }

  globaloptions::get().anglerep = useranglerep;
  globaloptions::get().verbose = storedverbose;
  globaloptions::get().dof = storeddof;
}


////////////////////////////////////////////////////////////////////////////

float measure_cost(Matrix& affmat, int input_dof)
{
  Tracer tr("measure_cost");
  // the most basic strategy - just do a single optimisation run at the
  //  specified dof
  int dof=input_dof;
  if (dof<6) { 
    cerr << "Erroneous dof " << dof << " : using 6 instead\n"; 
    dof=6; 
  }
  if (dof>12) {
    cerr << "Erroneous dof " << dof << " : using 12 instead\n"; 
    dof=12;
  }

  return costfn(affmat);
}  


void aligncog(Matrix& affmat)
{
  Tracer tr("aligncog");
  Matrix resmat = affmat;
  resmat(1,4) = resmat(2,4) = resmat(3,4) = 0.0;
  ColumnVector transoff;
  transoff = resmat.SubMatrix(1,3,1,3)
    * globaloptions::get().impair->testvol.cog("scaled_mm")
    - globaloptions::get().impair->refvol.cog("scaled_mm");
  resmat(1,4) = -transoff(1);
  resmat(2,4) = -transoff(2);
  resmat(3,4) = -transoff(3);
  affmat = resmat;
}


void alignpaxes(Matrix& affmat)
{
  Tracer tr("alignpaxes");
  Matrix paxref, paxtest;
  paxref = globaloptions::get().impair->refvol.principleaxes_mat();
  paxtest = globaloptions::get().impair->testvol.principleaxes_mat();
  ColumnVector rv, tv;
  for (int n=1; n<=2; n++) {
    rv = paxref.SubMatrix(1,3,n,n);
    tv = paxtest.SubMatrix(1,3,n,n);
    if (dot(rv,tv)<0) {
      rv=-rv;
      paxref.SubMatrix(1,3,n,n)=rv;
    }
  }
  if (paxref.Determinant()<0) 
    paxref.SubMatrix(1,3,3,3) = -paxref.SubMatrix(1,3,3,3);
  if (paxtest.Determinant()<0) 
    paxtest.SubMatrix(1,3,3,3) = -paxtest.SubMatrix(1,3,3,3);
  affmat = IdentityMatrix(4);
  affmat.SubMatrix(1,3,1,3) = paxref * paxtest.t();
  aligncog(affmat);
}


////////////////////////////////////////////////////////////////////////////

int optimise_strategy0(Matrix& matresult, float& fans, int max_iterations=4)
{
  Tracer tr("optimise_strategy0");
  // the most basic strategy - just do a single optimisation run at the
  //  specified dof

  if (globaloptions::get().verbose>3) {
    cout << "Using subset cost function" << endl;
  }

  ColumnVector params(12), params_N(12), param_tol(12), param_tol0(12), 
    param_tol1(12);
  int no_its=0;

  globaloptions::get().no_params = 12; // necessary for any subset_costfn call

  affmat2vector(matresult,12,params);
  globaloptions::get().refparams = params;

  // calculate the parameter tolerances
  param_tol0 = globaloptions::get().refparams;
  params12toN(param_tol0);
  set_param_tols(param_tol1,12);
  param_tol1 = param_tol1 + globaloptions::get().refparams;
  params12toN(param_tol1);
  param_tol = param_tol1 - param_tol0;


  params_N = globaloptions::get().refparams;
  params12toN(params_N);
  if (globaloptions::get().verbose>6) {
    cout << "Starting with " << params_N.t();
    cout << "  and tolerance " << param_tol.t();
  }
  optimise(params_N,globaloptions::get().parammask.Ncols(),
	   param_tol,no_its,&fans,subset_costfn,max_iterations);
  paramsNto12(params_N);
  params = params_N;

  vector2affine(params,12,matresult);
  return no_its;
}  



int optimise_strategy1(Matrix& matresult, float& fans, int input_dof, 
		       int max_iterations=4)
{
  Tracer tr("optimise_strategy1");
  // the most basic strategy - just do a single optimisation run at the
  //  specified dof
  int dof=input_dof;
  if (dof<6) { 
    cerr << "Erroneous dof " << dof << " : using 6 instead\n"; 
    dof=6; 
  }
  if (dof>12) {
    cerr << "Erroneous dof " << dof << " : using 12 instead\n"; 
    dof=12;
  }

  ColumnVector params(12), param_tol(12);
  int no_its=0;
  globaloptions::get().no_params = dof;
  set_param_tols(param_tol,12);  // 12 used to be dof
  affmat2vector(matresult,dof,params);
  //optimise(params,dof,param_tol,&no_its,&fans,costfn,max_iterations);
  if (globaloptions::get().verbose>6) {
    cout << "Starting with " << params.t();
    cout << "  and tolerance " << param_tol.t();
  }
  optimise(params,dof,param_tol,no_its,&fans,costfn,max_iterations);
  vector2affine(params,dof,matresult);
  return no_its;
}  

//-------------------------------------------------------------------------//


int sorted_posn(const float costval, const Matrix& opt_matrixlist)
{
  Tracer tr("sorted_posn");
  int n=1;
  if (opt_matrixlist.Nrows()<1) { return 1; }
  for (n=1; n<=opt_matrixlist.Nrows(); n++) {
    if (costval < opt_matrixlist(n,1)) {
      return n;
    }
  }
  return opt_matrixlist.Nrows()+1;
}
 

void delete_row(Matrix& mat, const int row) 
{
  Tracer tr("delete_row");
  if ((row<1) || (row>mat.Nrows()))  return;
  int cols = mat.Ncols();
  if (cols<1) return;
  Matrix temp(mat.Nrows()-1,cols);
  if (row>=2) {
    temp.SubMatrix(1,row-1,1,cols) = mat.SubMatrix(1,row-1,1,cols);
  }
  if (row<mat.Nrows()) {
    temp.SubMatrix(row,mat.Nrows()-1,1,cols) 
      = mat.SubMatrix(row+1,mat.Nrows(),1,cols);
  }
  mat = temp;
}
 

void insert_row(Matrix& mat, const int row, const Matrix& rowmat)
{
  Tracer tr("insert_row");
  Matrix temp;
  if ((row<1) || (row>(mat.Nrows()+1)))  return;
  int cols = mat.Ncols(), nmax = mat.Nrows();
  if (cols<1) return;
  if ((rowmat.Ncols()!=cols) || (rowmat.Nrows()!=1)) { return; }
  temp.ReSize(nmax + 1,cols);
  if ((nmax>0) && (row>=2)) {
    temp.SubMatrix(1,row-1,1,cols) = mat.SubMatrix(1,row-1,1,cols);
  }
  if ((nmax>0) && (row<=nmax)) {
    temp.SubMatrix(row+1,nmax+1,1,cols) = mat.SubMatrix(row,nmax,1,cols);
  }
  temp.SubMatrix(row,row,1,cols) = rowmat;
  mat = temp;
}


int add2list(const Matrix& rowmat, Matrix& opt_matrixlist, float rms_min=1.0) 
{
  Tracer tr("add2list");
  if (rowmat.Nrows()!=1) {
    cerr << "WARNING (add2list): cannot add matrix with " << rowmat.Nrows() 
	 << " rows to matrix\n";
    return -1;
  }
  Matrix mat1(4,4), mat2(4,4);
  reshape(mat1,rowmat.SubMatrix(1,1,2,17),4,4);
  float costval = rowmat(1,1);
  int nmax = opt_matrixlist.Nrows();
  for (int n=1; n<=nmax; n++) {
    reshape(mat2,opt_matrixlist.SubMatrix(n,n,2,17),4,4);
    if (rms_deviation(mat1,mat2)<rms_min) {
      if (costval < opt_matrixlist(n,1)) {
	// if it is close but better, then delete, and add in appropriate posn
	delete_row(opt_matrixlist,n);
	break;
      } else {
	// it is close but worse, so don't add it
	return -1;
      }
    }
  }
  // otherwise add this one to maintain sorted order (ascending) of cost
  int pos = sorted_posn(costval,opt_matrixlist);
  insert_row(opt_matrixlist,pos,rowmat);
  return pos;
}


void optimise_strategy3(Matrix& opt_matrixlist)
{
  Tracer tr("optimise_strategy3");
  // the basic low-level search strategy - it searches the cost function
  //  space (across rotations)  then selects the best candidates and optimises
  //  these too   : best used for the 8x subsampled (far too slow otherwise)
  // the returned matrixlist contains in each row the cost and reshaped matrix 
  //  for both optimised, and pre-optimised positions, ranked in ascending
  //  order of the OPTIMISED cost
  
  
  Matrix paramlist;
  volume<float> costs,tx,ty,tz,scale;
  globaloptions::get().currentcostfn = globaloptions::get().searchcostfn;
  search_cost(paramlist,costs,tx,ty,tz,scale);
  
  int dof = Min(globaloptions::get().dof,7);
  ColumnVector params_8(12);
  float costval=0.0, rms_min = estimate_scaling();
  Matrix reshapedmat(1,16), matresult(4,4), premat(4,4), matrow(1,34);
  opt_matrixlist.ReSize(0,34);
  // freely optimise each member of the parameter list (allow rotns to vary)
  int verbose = globaloptions::get().verbose;
  globaloptions::get().verbose-=2;
  if (verbose>=3) { 
    cout << "After free optimisation, parameters are:" << endl;
  }
  for (int n=1; n<=paramlist.Nrows(); n++) {
    params_8 = paramlist.SubMatrix(n,n,1,12).t();
    vector2affine(params_8,dof,matresult);
    premat = matresult;
    optimise_strategy1(matresult,costval,dof);
    // form the optimised and pre-optimised costs and matrices into a row
    reshape(reshapedmat,matresult,1,16);
    matrow(1,1) = costval;
    matrow.SubMatrix(1,1,2,17) = reshapedmat;
    matrow(1,18) = costfn(premat);
    reshape(reshapedmat,premat,1,16);
    matrow.SubMatrix(1,1,19,34) = reshapedmat;
    // add the row to the optimised list (ordered by optimised cost)
    add2list(matrow,opt_matrixlist,rms_min);
    if (verbose>=3) { 
      affmat2vector(matresult,12,params_8);
      cout << costval << " ::: " << params_8.t(); 
    }
  }
  globaloptions::get().verbose = verbose;
  globaloptions::get().currentcostfn = globaloptions::get().maincostfn;

  if (globaloptions::get().verbose>=3) {
    cout << "Parameters (1st column costs) are:\n" << opt_matrixlist << endl;
  }
}

//-------------------------------------------------------------------------//

void set_perturbations(Matrix& delta, Matrix& perturbmask)
{
  Tracer tr("set_perturbations");
  // set the perturbations to be applied to the pre-optimised
  //  solutions
  // the perturbations are given by the elementwise product of
  //  delta (a single column) and the individual columns of perturbmask
  // the number of columns of perturbmask set the number of perturbations
  //  to be tried - it should always include a zero column = unperturbed case
  ColumnVector coarserx, coarsery, coarserz, finerx, finery, finerz;
  set_rot_samplings(coarserx,coarsery,coarserz,finerx,finery,finerz);
  ColumnVector param_tol(12);
  set_param_tols(param_tol,12);
  // set the magnitude of the variations
  float delscale, delrx, delry, delrz;
  delrx = 3.0*param_tol(1); 
  delry = 3.0*param_tol(2);
  delrz = 3.0*param_tol(3);
  if (finerx.Nrows()>1) { delrx = 0.5*(finerx(2) - finerx(1)); }
  if (finery.Nrows()>1) { delry = 0.5*(finery(2) - finery(1)); }
  if (finerz.Nrows()>1) { delrz = 0.5*(finerz(2) - finerz(1)); }
  delscale=0.0;  
  if (globaloptions::get().dof>=7) delscale = 0.1;  // scale only set if allowed
  delta.ReSize(12,1);
  delta = 0.05;  // translation default
  delta(1,1) = delrx; delta(2,1) = delry;  delta(3,1) = delrz;
  delta(7,1) = delscale;  delta(8,1) = delscale;  delta(9,1) = delscale;
  // set the directions of the variations 
  // (+/- for each rotation and  4 scale settings: +/- 1 and 2 times delscale)
  perturbmask.ReSize(12,11);
  perturbmask = 0.0;
  perturbmask(1,1) = 1.0;
  perturbmask(1,2) = -1.0;
  perturbmask(2,3) = 1.0;
  perturbmask(2,4) = -1.0;
  perturbmask(3,5) = 1.0;
  perturbmask(3,6) = -1.0;
  perturbmask(7,7) = 1.0; 
  perturbmask(8,7) = 1.0; 
  perturbmask(9,7) = 1.0;
  perturbmask(7,8) = 2.0; 
  perturbmask(8,8) = 2.0; 
  perturbmask(9,8) = 2.0;
  perturbmask(7,9) = -1.0; 
  perturbmask(8,9) = -1.0; 
  perturbmask(9,9) = -1.0;
  perturbmask(7,10)=-2.0; 
  perturbmask(8,10)=-2.0; 
  perturbmask(9,10)=-2.0;
}


////////////////////////////////////////////////////////////////////////////
// SUPPORT FUNCTIONS FOR READING AND INITIALISING IMAGES AND MATRICES 
//  - VERY FLIRT SPECIFIC
////////////////////////////////////////////////////////////////////////////


Matrix scalemat(const Matrix& mat)
{
  // Adjust for scaled coordinates
  //   e.g. basescale=0.1mm => orig coord = 0.2 --> scaled coord = 2
  //   all external matrices in orig coord (mm) but internally work in
  //   scaled coords as the input volumes have their voxel size scaled
  Matrix orig2scaled_coord(4,4), returnmat;
  orig2scaled_coord=IdentityMatrix(4);
  orig2scaled_coord(1,1)=1/globaloptions::get().basescale;
  orig2scaled_coord(2,2)=1/globaloptions::get().basescale;
  orig2scaled_coord(3,3)=1/globaloptions::get().basescale;
  //  need internal version of initmat to go to/from scaled coords
  returnmat = orig2scaled_coord * mat * orig2scaled_coord.i();
  return returnmat;
}


Matrix qsform_init_mat(const volume<float>& refvol, const volume<float>& testvol) 
{
  Matrix returnmat;
  returnmat = IdentityMatrix(4);
  if ( ((refvol.qform_code() + refvol.sform_code())>0) && 
       ((testvol.qform_code() + testvol.sform_code())>0) ) {
    returnmat = refvol.sampling_mat()
      * refvol.newimagevox2mm_mat().i() * testvol.newimagevox2mm_mat()
      * testvol.sampling_mat().i();
    returnmat = scalemat(returnmat);
  }
  return returnmat;
}


void set_initmat(const volume<float>& refvol, const volume<float>& testvol)
{
  // Initialise with user-supplied matrix
  if (globaloptions::get().initmatfname.size()>0) {
    globaloptions::get().initmat =
      read_ascii_matrix(globaloptions::get().initmatfname);
  } else {
    // If not matrix then use s/q form info (unless told to ignore it)
    if (globaloptions::get().initmatsqform) {
      globaloptions::get().initmat = qsform_init_mat(refvol,testvol);
    }
  }

  // Adjust for scaled coordinates
  globaloptions::get().initmat = scalemat(globaloptions::get().initmat);

  if (globaloptions::get().verbose>=2) {
    cout << "Init Matrix = \n" << globaloptions::get().initmat << endl;
  }
}


void double_end_slices(volume<float>& testvol)
{
  // this is necessary for single slice volumes so that interpolation can
  //  be done (in general it is good to do for small number of slices so
  //  that the end ones get counted and not de-weighted by the cost fns)
  volume<float> newtestvol(testvol.xsize(),testvol.ysize(),testvol.zsize()+2);
  newtestvol.setdims(testvol.xdim(),testvol.ydim(),8.0f);
  for (int z=0; z<= testvol.zsize()+1; z++) {
    for (int y=0; y<testvol.ysize(); y++) {
      for (int x=0; x<testvol.xsize(); x++) {
	int ez = z-1;
	if (ez<0) ez=0;
	if (ez>=testvol.zsize())  ez=testvol.zsize()-1;
	newtestvol(x,y,z) = testvol(x,y,ez);
      }
    }
  }
  testvol = newtestvol;
}


int get_testvol(volume<float>& testvol)
{
  Tracer tr("get_testvol");
  short dtype=0;
  float minval=0.0, maxval=0.0;
  if (!read_testvol) {
    FLIRT_read_volume(testvol,globaloptions::get().inputfname);
    dtype = NEWIMAGE::dtype(globaloptions::get().inputfname);
    if (!globaloptions::get().forcedatatype)
      globaloptions::get().datatype = dtype;
    if (testvol.zsize()==1) {
      double_end_slices(testvol);
    }

    minval = testvol.robustmin();
    maxval = testvol.robustmax();
    if (globaloptions::get().clamping) clamp(testvol,minval,maxval);
    
    if (globaloptions::get().useweights) {
      if (globaloptions::get().testweightfname.length()>0) {
	FLIRT_read_volume(global_testweight,globaloptions::get().testweightfname);
	if (global_testweight.zsize()==1) {
	  double_end_slices(global_testweight);
	}
      } else {
	global_testweight = testvol;  // Fix size and dimensions
	global_testweight = 1.0;      // Set all elements to unity weighting
      }
    }
    global_init_testvol = testvol;
    global_init_testweight = global_testweight;
    read_testvol = true;
  } else {
    testvol = global_init_testvol;
    global_testweight = global_init_testweight;
  }

  if (globaloptions::get().verbose>=2) {
    cout << "Testvol sampling matrix =\n" << testvol.sampling_mat() << endl;
    cout << "Testvol Data Type = " << dtype << endl;
    cout << "Testvol intensity ";
    if (globaloptions::get().clamping) {
      cout << "clamped between " << minval << " and " << maxval << endl;
    } else {
      cout << "between " << testvol.min() << " and " << testvol.max() << endl;
    }
  }
  return 0;
}  


int get_refvol(volume<float>& refvol)
{
  Tracer tr("get_refvol");
  FLIRT_read_volume(refvol,globaloptions::get().reffname);
  if ((refvol.zsize()==1) && (globaloptions::get().do_optimise)) {
    double_end_slices(refvol);
  }

  float minval=0.0, maxval=0.0;
  minval = refvol.robustmin();
  maxval = refvol.robustmax();
  if (globaloptions::get().clamping) clamp(refvol,minval,maxval);

  if (globaloptions::get().useweights) {
    if (globaloptions::get().refweightfname.length()>0) {
      FLIRT_read_volume(global_refweight,globaloptions::get().refweightfname);
      if (global_refweight.zsize()==1) { 
	double_end_slices(global_refweight);
      }
    } else {
      global_refweight = refvol;   // Fix size and dimensions
      global_refweight = 1.0;      // Set all elements to unity weighting
    }
  }

  if (globaloptions::get().fmapfname.length()>0) {
      FLIRT_read_volume(global_fmap,globaloptions::get().fmapfname);
      if (globaloptions::get().fmapmaskfname.length()>0) {
	FLIRT_read_volume(global_fmap_mask,globaloptions::get().fmapmaskfname);
      } else {
	global_fmap_mask = global_fmap*0.0f + 1.0f;
      }
  }

  if (globaloptions::get().useseg) {
    if (!globaloptions::get().usecoords) {
      FLIRT_read_volume(global_seg,globaloptions::get().wmsegfname);
      if (global_seg.zsize()==1) { 
	double_end_slices(global_seg);
      }
    } else {
      global_coords = read_ascii_matrix(globaloptions::get().wmcoordsfname);
      global_norms = read_ascii_matrix(globaloptions::get().wmnormsfname);
      if ( (global_coords.Nrows()==0) || (global_norms.Nrows()==0) || 
	   (global_coords.Nrows() != global_norms.Nrows()) ) 
	{
	  cerr << "Coordinate and Normal matrices " << globaloptions::get().wmcoordsfname << " and " << globaloptions::get().wmnormsfname << " are either zero size or different sizes" << endl; 
	  exit(1);
	}
    }
  }
  
  if (globaloptions::get().verbose>=2) {
    cout << "Refvol intensity ";
    if (globaloptions::get().clamping) {
      cout << "clamped between " << minval << " and " << maxval << endl;
    } else {
      cout << "between " << refvol.min() << " and " << refvol.max() << endl;
    }
  }
  return 0;
}

volume<float> filter_subsample_by_2(const volume<float>& vin)
{
  return subsample_by_2(vin);
}

volume<float> filter_blur(const volume<float>& vin)
{
  volume<float> tmpvol;
  tmpvol = blur(vin,global_sampling);
  return tmpvol;
}

volume<float> filter_resamp_blur(const volume<float>& vin)
{
  volume<float> tmpvol;
  tmpvol = blur(vin,global_sampling);
  tmpvol = isotropic_resample(tmpvol,global_sampling);
  return tmpvol;
}

int filter_weight(volume<float>& blur_w, const volume<float>& weight,
		  volume<float> (*filter_func)(const volume<float>&))
{
  // implements filter(W).*Thresh(filter(B)) 
  //                  where W is the weight and B is a binarised weight
  // filter(...) = isotropic_resample(blur(...))
  float thresh1=0.01, thresh2=0.9;
  volume<float> tmpvol;
  // form B
  tmpvol = binarise(weight,thresh1);
  // form filter(B)
  tmpvol = (*filter_func)(tmpvol);
  // form Thresh(filter(B))
  tmpvol.binarise(thresh2);
  // now form filter(W)
  // NB: if blur_w = weight on input, then the following changes weight!
  blur_w = (*filter_func)(weight);
  // now form filter(W).*Thresh(filter(B))
  blur_w *= tmpvol;
  return 0;
}


int filter_weight(volume<float>& blur_w, const volume<float>& weight, float sampling,
		  volume<float> (*filter_func)(const volume<float>&))
{
  global_sampling = sampling;
  return filter_weight(blur_w,weight,filter_func);
}


int filter_image(volume<float>& blur_im, const volume<float>& orig_im,
		 const volume<float>& weight, 
		 volume<float> (*filter_func)(const volume<float>&))
{
  // implements filter(I.*B)./filter(B) where I in the image and B is a binarised weight
  // filter(...) = isotropic_resample(blur(...))
  // form B then filter(I.*B)
  float thresh=0.01;
  volume<float> tmpvol;
  tmpvol = binarise(weight,thresh);
  // NB: if blur_im = orig_im on input, then the following line changes orig_im!
  blur_im = (*filter_func)(orig_im*tmpvol);  
  // form filter(B)
  tmpvol = binarise(weight,thresh);
  tmpvol = (*filter_func)(tmpvol);
  // now form filter(I.*B)./filter(B) - with safe (non-zero) division
  blur_im = divide(blur_im,tmpvol,tmpvol);
  return 0;
}


int filter_image(volume<float>& blur_im, const volume<float>& orig_im,
		 const volume<float>& weight, bool useweight,
		 volume<float> (*filter_func)(const volume<float>& ))
{
  if (useweight) {
    filter_image(blur_im,orig_im,weight,filter_func);
  } else {
    blur_im = (*filter_func)(orig_im);
  }
  return 0;
}


int filter_image(volume<float>& blur_im, const volume<float>& orig_im,
		 const volume<float>& weight, float sampling, bool useweight,
		 volume<float> (*filter_func)(const volume<float>& ))
{
  global_sampling = sampling;
  return filter_image(blur_im, orig_im, weight, useweight, filter_func);
}



int resample_refvol(volume<float>& refvol, float sampling=1.0)
{
  Tracer tr("resample_refvol");
  float sampl=sampling;
  if (sampl<1e-5) { 
    cerr << "WARNING: sampling " << sampl << " less than 0.00001 mm" << endl
	 << "         Setting to 1.0 sampling instead" << endl;
    sampl=1.0;
  }
  if (sampl<1e-3) { 
    cerr << "WARNING: under minimum sampling of 0.001 mm" << endl;
  }

  if (globaloptions::get().verbose >= 2) 
    cout << "Resampling refvol isotropically: " << sampling << endl;

  // isotropically resample the volume
  filter_image(refvol,refvol,global_refweight,sampl,
	       globaloptions::get().useweights,filter_resamp_blur);

  if (globaloptions::get().verbose>=2) print_volume_info(refvol,"Refvol");      
  
  return 0;
}


////////////////////////////////////////////////////////////////////////////


// template for either volume or volume4D (and no other)
template <class V>
int output_dtype(const V& outvol) 
{
  // used to determine whether a float should be forced for a mask
  //  output (only if not overriden by user forcedatatype)
  int dtype;
  dtype = globaloptions::get().datatype;
  // do not change this if the user is forcing a datatype
  if (!globaloptions::get().forcedatatype) {
    // only need to change it if it is not currently a floating type
    if ( (dtype!=DT_FLOAT) && (dtype!=DT_DOUBLE) ) {
      // only force floating output if range is currently < 1.5
      if ((outvol.max() - outvol.min())<1.5) {
	dtype=DT_FLOAT;
      }
    }
  }
  return dtype;
}

////////////////////////////////////////////////////////////////////////////

// this does the applyxfm!
void do_applyxfm()
{
  Tracer tr("do_applyxfm");
  volume<float> refvol;
  volume4D<float> testvol;

  // want unity basescale for transformed output
  globaloptions::get().basescale = 1.0;

  // set up image pair and global pointer
  
  FLIRT_read_volume(refvol,globaloptions::get().reffname);
  FLIRT_read_volume4D(testvol,globaloptions::get().inputfname);

  if ( (refvol.sform_code()!=NIFTI_XFORM_UNKNOWN) && 
       (testvol[0].sform_code()!=NIFTI_XFORM_UNKNOWN) ) {
    if (globaloptions::get().verbose>0) {
      cerr << "WARNING: Both reference and input images have an sform matrix set" << endl;
    }
  }

  short dtype;
  dtype = NEWIMAGE::dtype(globaloptions::get().inputfname);
  if (!globaloptions::get().forcedatatype)
    globaloptions::get().datatype = dtype;
  

  set_initmat(refvol,testvol[0]);

  if (globaloptions::get().verbose>0) {
    if (refvol.sform_code()!=NIFTI_XFORM_UNKNOWN) {
      cout << "The output image will use the sform from the reference image" << endl;    
    } else if (testvol[0].sform_code()!=NIFTI_XFORM_UNKNOWN) {
      cout << "The output image will use the transformed sform from the input image" << endl;    
    }
  }

  if ( (globaloptions::get().verbose>0) || (globaloptions::get().printinit)) {
    cout << "Init Matrix = \n" << globaloptions::get().initmat << endl;
  }
  
  if (globaloptions::get().iso) {
    resample_refvol(refvol,globaloptions::get().isoscale);
  }

  float min_sampling_ref=1.0;
  min_sampling_ref = Min(refvol.xdim(),Min(refvol.ydim(),refvol.zdim()));

  // // need the following to set up BBR fieldmap parameters (MJ TODO - UNTESTED!!!)
  // setup_costfn(globaloptions::get().impair, globaloptions::get().maincostfn,
  // 	       globaloptions::get().no_bins,
  // 	       globaloptions::get().smoothsize,globaloptions::get().fuzzyfrac);  

  if (globaloptions::get().outputfname.size()>0) {
    volume4D<float> outputvol;
    for (int t0=testvol.mint(); t0<=testvol.maxt(); t0++) {
      int tref=t0-testvol.mint();
      outputvol.addvolume(refvol);
      if ((globaloptions::get().interpmethod != NearestNeighbour) &&
	  (globaloptions::get().interpblur)) {
	filter_image(testvol[t0],testvol[t0],testvol[t0],min_sampling_ref,
		     false,filter_blur);
      }
      
      if (globaloptions::get().verbose>=2) { 
	print_volume_info(refvol,"refvol"); 
	print_volume_info(testvol,"inputvol"); 
      }
      
      final_transform(testvol[t0],refvol,globaloptions::get().initmat,outputvol[tref]);
    }
    int outputdtype = output_dtype(outputvol);
    outputvol.setDisplayMaximumMinimum(0,0);
    outputvol.settdim(testvol.tdim());
    save_volume4D_dtype(outputvol,globaloptions::get().outputfname.c_str(),
			outputdtype);
    if (globaloptions::get().verbose>=2) {
      print_volume_info(outputvol,"Resampled volume");
    }
  }

  if (globaloptions::get().outputmatascii.size()>0) {
    save_matrix_data(globaloptions::get().initmat);
  }

  exit(0);
}


////////////////////////////////////////////////////////////////////////////

// ROUTINES FOR DEALING WITH THE SCHEDULE FILE

int firstelementless(const RowVector& r1, const RowVector& r2)
{
  Tracer tr("firstelementless");
  return (r1(1) < r2(1));
}


void stripleadingspace(string& line)
{
  Tracer tr("stripleadingspace");
  if (line.length()<1) return;
  string::size_type i2=line.find_first_not_of("	 ");
  if (i2==string::npos) {
    return;
  }
  if (i2!=0) {
    line = line.substr(i2,string::npos);
  }
}


void striptrailingspace(string& line)
{
  Tracer tr("striptrailingspace");
  if (line.length()<1) return;
  string::size_type i2=line.find_last_not_of("	 ");
  if (i2==string::npos) {
    line="";
    return;
  }
  if (i2!=0) {
    line = line.substr(0,i2+1);
  }
}



int parseline(const string& inname, std::vector<string>& words) 
{
  Tracer tr("parseline");
  string name = inname;
  string::size_type i1;
  words.clear();
  while (name.length()>=1) {
    stripleadingspace(name);
    striptrailingspace(name);
    i1=name.find_first_of("	 ");
    words.push_back(name.substr(0,i1));
    name.erase(0,i1);
  }
  return 0;
}



int setmatvariable(const string& name, MatVecPtr& namedmat)
{
  Tracer tr("setmatvariable");
  if (name.size()<1) return -2;
  if (name == "S") {
    namedmat = &globaloptions::get().searchoptmat;
    return 0;
  } else if (name == "P") {
    namedmat = &globaloptions::get().preoptsearchmat;
    return 0;
  } else if (name == "U") {
    namedmat = &(globaloptions::get().usrmat[0]);
    return 0;
  } else if (name[0] == 'U') {
    int idx = name[1] - 'A' + 1;
    if ((idx<1) || (idx>27)) {
      if ((name[1]<'0') || (name[1]>'9')) {
	cerr << "Can only accept UA to UZ, not " << name << endl;
	return -3;
      } else {
	idx = 0;
      }
    }
    namedmat = &(globaloptions::get().usrmat[idx]);
    return 0;
  } else {
    cerr << "Cannot interpret " << name << " as a Matrix name" << endl;
    return -1;
  }
}


int setscalarvariable(const string& name, float& namedscalar)
{
  Tracer tr("setscalarvariable");
  if (name.size()<1) return -2;
  if (name == "MAXDOF") {
    namedscalar = (float) globaloptions::get().dof;
    return 0;
  } else if (name == "MINSAMPLING") {
    namedscalar = (float) globaloptions::get().min_sampling;
    return 0;
  } else if (isalpha(name[0])) {
    if (globaloptions::get().verbose>20) {
      cerr << "Cannot interpret " << name << " as a scalar variable name" << endl;
    }
    namedscalar = 0.0;
    return -1;
  } else {
    namedscalar = atof(name.c_str());
  }
  return 0;
}


int setscalarvariable(const string& name, int& namedscalar)
{
  Tracer tr("setscalarvariable");
  float temp=0.0;
  int retval = setscalarvariable(name,temp);
  namedscalar = MISCMATHS::round(temp);
  return retval;
}


int parsematname(const string& inname, MatVecPtr& usrdefmat, int& r1, int& r2) 
{
  Tracer tr("parsematname");
  string name = inname;
  stripleadingspace(name);
  striptrailingspace(name);
  string::size_type colon, hyphen;
  if (name.length()<1) {
    return 0;
  }
  colon=name.find_first_of(":");
  hyphen=name.find_first_of("-");
  string basename=name, row1="", row2="";
  if (colon!=string::npos) {
    basename = name.substr(0,colon);
  }
  if ((colon!=string::npos) && (hyphen!=string::npos)) {
    // there is a colon and a hyphen
    row1 = name.substr(colon+1,hyphen-colon-1);
    row2 = name.substr(hyphen+1);
  } else {
    if ((colon!=string::npos) && (hyphen==string::npos)) {
      // there is a colon but no hyphen
      row1 = name.substr(colon+1);
      row2 = row1;
    }
  }
  if (row1.length()<1) row1="1";
  if (row2.length()<1)  row2="999999";

  setmatvariable(basename, usrdefmat);
  setscalarvariable(row1,r1);
  setscalarvariable(row2,r2);
  if (r1<1) r1=1;
  if (r2<r1) r2=r1;
  return 0;
}


//-----------------------------------------------------------------------//

void usrcopy(MatVecPtr usrsrcmat, MatVecPtr usrdestmat, 
	     unsigned int usrsrcrow1, unsigned int usrsrcrow2) 
{
  Tracer tr("usrcopy");
  // COPY src -> dest
  for (unsigned int crow=usrsrcrow1; crow<=Min(usrsrcrow2,usrsrcmat->size()); 
       crow++)
    {
      usrdestmat->push_back((*usrsrcmat)[crow-1]);
    }
}


void usrclear(MatVecPtr usrsrcmat) 
{
  Tracer tr("usrclear");
  // CLEAR src
  if (usrsrcmat->size() > 0) {
    usrsrcmat->erase(usrsrcmat->begin(),usrsrcmat->end());
  }
}



void usrprint(MatVecPtr usrsrcmat, unsigned int row1, unsigned int row2) 
{
  Tracer tr("usrprint");
  // PRINT src
  for (unsigned int r=row1; r<=Min(usrsrcmat->size(),row2); r++) {
    cout << (*usrsrcmat)[r-1];
  }
}


int usrsave(string filename, MatVecPtr usrsrcmat, 
	     unsigned int row1, unsigned int row2) 
{
  Tracer tr("usrsave");
  // SAVE src

  ofstream fptr(filename.c_str());
  if (!fptr) { 
    cerr << "Could not open file " << filename << " for writing" << endl;
    return -1;
  }
  for (unsigned int r=row1; r<=Min(usrsrcmat->size(),row2); r++) {
    fptr << (*usrsrcmat)[r-1];
  }
  fptr << endl;
  fptr.close();
  return 0;
}


int usrread(string filename, MatVecPtr usrsrcmat)
{
  Tracer tr("usrread");
  // READ src

  ifstream fptr(filename.c_str());
  if (!fptr) { 
    cerr << "Could not open file " << filename << " for reading" << endl;
    return -1;
  }
  usrsrcmat->clear();
  RowVector currow(17);
  while (!fptr.eof()) {
    currow = 0.0;
    for (unsigned int c=1; c<=17; c++) {
      fptr >> currow(c);
    }
    if (!fptr.eof()) {
      usrsrcmat->push_back(currow);
    }
  }
  fptr.close();
  return 0;
}


int usrsetrow(MatVecPtr usrsrcmat,const std::vector<string> &words)
{
  Tracer tr("usrsetrow");
  // SETROW mat p1 p2 .. p16

  RowVector currow(17);
  currow = 0.0;
  for (unsigned int c=2; c<=17; c++) {
    currow(c) = atof(words[c].c_str());
  }
  usrsrcmat->push_back(currow);
  return 0;
}

int usrsetrowqsform(MatVecPtr usrsrcmat, const volume<float>& refvol, 
		    const volume<float>& testvol)
{
  Tracer tr("usrsetrowqsform");
  // SETROWQSFORM mat

  Matrix qsinitmat(4,4);
  // set the qsform initialisation matrix (if not already the -init)
  if (globaloptions::get().initmatsqform) {
    qsinitmat = IdentityMatrix(4);
  } else {
    // account for user-defined initmat (when it is not the sqform init)
    qsinitmat = qsform_init_mat(refvol,testvol) * globaloptions::get().initmat.i();
  }
  RowVector currow(17);
  currow = 0.0;
  unsigned int c=2;
  for (unsigned int col=1; col<=4; col++) {
    for (unsigned int row=1; row<=4; row++) {
      currow(c++) = qsinitmat(row,col);
    }
  }
  usrsrcmat->push_back(currow);
  return 0;
}

int usrsetrowparams(MatVecPtr usrsrcmat,const std::vector<string> &words)
{
  Tracer tr("usrsetrowparams");
  // SETROWPARAMS mat p1 p2 ...

  ColumnVector params(12);
  Matrix affmat(4,4), reshaped(1,16);
  RowVector rowresult(17);
  float costval=0.0;

  params = 0.0;
  for (unsigned int c=2; c<=13; c++) {
    params(c-1) = atof(words[c].c_str());
  }

  vector2affine(params,12,affmat);
  reshape(reshaped,affmat,1,16);
  rowresult(1) = costval;
  rowresult.SubMatrix(1,1,2,17) = reshaped;
  usrsrcmat->push_back(rowresult);

  return 0;
}



int usrsetoption(const std::vector<string> &words)
{
  Tracer tr("usrsetoption");
  // SETOPTION name values

  int len = words.size();
  if (len<3) {
    cerr << "Must specify an option and a value" << endl;
    return -2;
  }
    
  string option = words[1];

  ColumnVector fvalues(len-2);
  float fval=0.0;

  for (int i=2; i<len; i++) {
    setscalarvariable(words[i],fval);
    fvalues(i-1) = fval;
  }

  if (globaloptions::get().verbose>2) {
    cout << "setoptions: " << option << " " << fvalues.t() << endl;
  }

  if (option=="smoothing") {
    globaloptions::get().smoothsize = fvalues(1);
    if (globaloptions::get().impair)
      globaloptions::get().impair->smoothsize = globaloptions::get().smoothsize;
    return 0;
  } else if (option=="fuzzyfraction") {
    globaloptions::get().fuzzyfrac = fvalues(1);
    if (globaloptions::get().impair)
      globaloptions::get().impair->fuzzyfrac = globaloptions::get().fuzzyfrac;
    return 0;
  } else if (option=="optimisationtype") {
    globaloptions::get().optimisationtype = words[2];
    return 0;
  } else if (option=="costfunction") {
    if (costfn_type(words[2])==Unknown) {
      cerr << "Unrecognised cost function: " << words[2] << endl;
      exit(-1);
    }
    globaloptions::get().currentcostfn = costfn_type(words[2]);
    globaloptions::get().maincostfn = costfn_type(words[2]);
    globaloptions::get().searchcostfn = costfn_type(words[2]);
    if (globaloptions::get().impair)
      setcostfntype(globaloptions::get().currentcostfn);
    return 0;
  } else if (option=="tolerance") {
    globaloptions::get().tolerance 
       = fvalues/globaloptions::get().requestedscale;
    // Note: division by requestedscale used so that the absolute
    //       tolerance used at this scale is that specified by the user
    return 0;
  } else if (option=="rescaletolerance") {
    globaloptions::get().tolerance *= fvalues;
    // Note: this will apply for all scales hereon
    return 0;
  } else if (option=="boundguess") {
    globaloptions::get().boundguess = fvalues;
    return 0;
  } else if (option=="nosubset") {
    globaloptions::get().usrsubset = false;
    return 0;
  } else if (option=="paramsubset") {
    globaloptions::get().usrsubset = true;
    int noparams = (int) fvalues(1);
    if (fvalues.Nrows()<(12*noparams+1)) {
      cerr << "Must specify at least " << 12*noparams << " for a " 
	   << noparams << " parameter subset mask" << endl;
      return -3;
    }
    globaloptions::get().parammask.ReSize(12,noparams);
    globaloptions::get().parammask = 0;
    for (int n=1; n<=noparams; n++) {
      globaloptions::get().parammask.SubMatrix(1,12,n,n) = 
	fvalues.SubMatrix((n-1)*12+2,n*12+1,1,1);
    }
    return 0;
  } else if (option=="minsampling") {
    globaloptions::get().min_sampling = fvalues(1);
    return 0;
  } else if (option=="bbrstep") {
    globaloptions::get().impair->set_bbr_step(fvalues(1));
    return 0;
  } else {
    cerr << "Option " << option << " is unrecognised - ignoring" << endl;
  }
  return -1;
}


void usraligncog(MatVecPtr stdresultmat, 
		 MatVecPtr usrmatptr, 
		 unsigned int usrrow1, unsigned int usrrow2)
{
  Tracer tr("usraligncog");
  Matrix delta, perturbmask, matresult;
  ColumnVector params(12);
  RowVector rowresult(17);
  for (unsigned int crow=usrrow1; crow<=Min(usrrow2,usrmatptr->size()); crow++)
    {
      Matrix reshaped = (*usrmatptr)[crow-1].SubMatrix(1,1,2,17);
      reshape(matresult,reshaped,4,4);
      aligncog(matresult);
      reshape(reshaped,matresult,1,16);
      rowresult(1) = (*usrmatptr)[crow-1].SubMatrix(1,1,1,1).AsScalar();
      rowresult.SubMatrix(1,1,2,17) = reshaped;
      // store result
      stdresultmat->push_back(rowresult);
    }
}


void usralignpaxes(MatVecPtr stdresultmat, 
		   MatVecPtr usrmatptr, 
		   unsigned int usrrow1, unsigned int usrrow2)
{
  Tracer tr("usralignpaxes");
  Matrix delta, perturbmask, matresult;
  ColumnVector params(12);
  RowVector rowresult(17);
  for (unsigned int crow=usrrow1; crow<=Min(usrrow2,usrmatptr->size()); crow++)
    {
      Matrix reshaped = (*usrmatptr)[crow-1].SubMatrix(1,1,2,17);
      reshape(matresult,reshaped,4,4);
      alignpaxes(matresult);
      reshape(reshaped,matresult,1,16);
      rowresult(1) = (*usrmatptr)[crow-1].SubMatrix(1,1,1,1).AsScalar();
      rowresult.SubMatrix(1,1,2,17) = reshaped;
      // store result
      stdresultmat->push_back(rowresult);
    }
}


void usrsort(MatVecPtr usrsrcmat) 
{
  Tracer tr("usrsort");
  // SORT src
  sort(usrsrcmat->begin(),usrsrcmat->end(),firstelementless);
}


void usrdualsort(MatVecPtr usrsrcmat1, MatVecPtr usrsrcmat2) 
{
  Tracer tr("usrdualsort");
  // DUALSORT
  if (usrsrcmat1->size() != usrsrcmat2->size()) {
    cerr << "Cannot dual sort matrices of unequal size";
    return;
  }
  RowVector tmprow(34);
  MatVec optmatsorted(0);
  for (unsigned int i=0; i<usrsrcmat1->size(); i++) {
    tmprow.SubMatrix(1,1,1,17) = (*usrsrcmat1)[i];
    tmprow.SubMatrix(1,1,18,34) = (*usrsrcmat2)[i];
    optmatsorted.push_back(tmprow);
  }
  sort(optmatsorted.begin(),optmatsorted.end(),firstelementless);
  usrclear(usrsrcmat1);
  usrclear(usrsrcmat2);
  tmprow.ReSize(17);
  for (unsigned int i=0; i<optmatsorted.size(); i++) {
    tmprow = optmatsorted[i].SubMatrix(1,1,1,17);
    usrsrcmat1->push_back(tmprow);
    tmprow = optmatsorted[i].SubMatrix(1,1,18,34);
    usrsrcmat2->push_back(tmprow);
  }
}


void usrsearch(MatVecPtr searchoptmat, MatVecPtr preoptsearchmat, int sdof=12) 
{
  Tracer tr("usrsearch");
  // SEARCH
  Matrix opt_matrixlist;
  globaloptions::get().searchdof = sdof;
  optimise_strategy3(opt_matrixlist);
  MatVec optmatsorted(0);
  RowVector tmprow;
  for (int i=1; i<=opt_matrixlist.Nrows(); i++) {
    tmprow = opt_matrixlist.SubMatrix(i,i,1,34);
    optmatsorted.push_back(tmprow);
  }
  sort(optmatsorted.begin(),optmatsorted.end(),firstelementless);
  for (unsigned int i=0; i<optmatsorted.size(); i++) {
    tmprow = optmatsorted[i].SubMatrix(1,1,1,17);
    searchoptmat->push_back(tmprow);
    tmprow = optmatsorted[i].SubMatrix(1,1,18,34);
    preoptsearchmat->push_back(tmprow);
  }
}


int usrreadparams(string filename, MatVecPtr usrsrcmat)
{
  Tracer tr("usrreadparams");
  // READ src

  ifstream fptr(filename.c_str());
  if (!fptr) { 
    cerr << "Could not open file " << filename << " for reading" << endl;
    return -1;
  }
  usrsrcmat->clear();
  ColumnVector params(12);
  Matrix affmat(4,4), reshaped(1,16);
  RowVector rowresult(17);
  float costval;
  while (!fptr.eof()) {
    params = 0.0;
    fptr >> costval;
    for (unsigned int c=1; c<=12; c++) {
      fptr >> params(c);
    }
    if (!fptr.eof()) {
      vector2affine(params,12,affmat);
      reshape(reshaped,affmat,1,16);
      rowresult(1) = costval;
      rowresult.SubMatrix(1,1,2,17) = reshaped;
      usrsrcmat->push_back(rowresult);
    }
  }
  fptr.close();
  return 0;
}


int usrsaveparams(string filename, MatVecPtr usrmatptr, 
		    unsigned int usrrow1, unsigned int usrrow2)
{
  Tracer tr("usrsaveparams");
  // PRINTPARAMS
  ofstream fptr(filename.c_str());
  if (!fptr) { 
    cerr << "Could not open file " << filename << " for writing" << endl;
    return -1;
  }
  Matrix matresult;
  ColumnVector params(12);
  for (unsigned int crow=usrrow1; crow<=Min(usrrow2,usrmatptr->size()); crow++)
    {
      // the pre-optimised case with perturbations
      Matrix reshaped = (*usrmatptr)[crow-1].SubMatrix(1,1,2,17);
      reshape(matresult,reshaped,4,4);
      affmat2vector(matresult,12,params);
      fptr << ((*usrmatptr)[crow-1])(1) << " " << params.t() << endl;
    }
  fptr << endl;
  fptr.close();
  return 0;
}


void usrprintparams(MatVecPtr usrmatptr, 
		    unsigned int usrrow1, unsigned int usrrow2)
{
  Tracer tr("usrprintparams");
  // PRINTPARAMS
  Matrix matresult;
  ColumnVector params(12);
  for (unsigned int crow=usrrow1; crow<=Min(usrrow2,usrmatptr->size()); crow++)
    {
      // the pre-optimised case with perturbations
      Matrix reshaped = (*usrmatptr)[crow-1].SubMatrix(1,1,2,17);
      reshape(matresult,reshaped,4,4);
      affmat2vector(matresult,12,params);
      cout << ((*usrmatptr)[crow-1])(1) << " " << params.t() << endl;
    }
}


void usrmeasurecost(MatVecPtr stdresultmat, 
		    MatVecPtr usrmatptr, 
		    unsigned int usrrow1, unsigned int usrrow2, int usrdof, 
		    ColumnVector& usrperturbation, bool usrperturbrelative)
{
  Tracer tr("usrmeasurecost");
  // MEASURE COST
  Matrix delta, perturbmask, matresult;
  set_perturbations(delta,perturbmask);
  int dof = Min(globaloptions::get().dof,usrdof);
  ColumnVector params(12);
  RowVector rowresult(17);
  for (unsigned int crow=usrrow1; crow<=Min(usrrow2,usrmatptr->size()); crow++)
    {
      // the pre-optimised case with perturbations
      Matrix reshaped = (*usrmatptr)[crow-1].SubMatrix(1,1,2,17);
      reshape(matresult,reshaped,4,4);
      affmat2vector(matresult,12,params);
      // use the elementwise product to produce the perturbation
      if (usrperturbrelative) {
	params += SP(usrperturbation,delta); // rel
      } else {
	params += usrperturbation; // abs
      }
      vector2affine(params,12,matresult);
      
      float costval=0.0;
      costval = measure_cost(matresult,dof);
      reshape(reshaped,matresult,1,16);
      rowresult(1) = costval;
      rowresult.SubMatrix(1,1,2,17) = reshaped;
      // store result
      stdresultmat->push_back(rowresult);
    }
}


void usrgridmeasurecost(MatVecPtr stdresultmat, 
			MatVecPtr usrmatptr, 
			unsigned int usrrow1, unsigned int usrrow2, int usrdof, 
			ColumnVector& usrperturbation1, ColumnVector& usrpertstep,
			ColumnVector& usrperturbation2, bool usrperturbrelative)
{
  Tracer tr("usrgridmeasurecost");
  // OPTIMISE
  Matrix delta, perturbmask, matresult;
  set_perturbations(delta,perturbmask);
  int dof = Min(globaloptions::get().dof,usrdof);
  ColumnVector params0(12), params(12), usrperturbation;
  RowVector rowresult(17);
  for (unsigned int crow=usrrow1; crow<=Min(usrrow2,usrmatptr->size()); crow++)
    {
      Matrix reshaped = (*usrmatptr)[crow-1].SubMatrix(1,1,2,17);
      reshape(matresult,reshaped,4,4);
      affmat2vector(matresult,12,params0);
      ColumnVector nsteps(usrperturbation1.Nrows());
      for (int n=1; n<=nsteps.Nrows(); n++) {
	nsteps(n) = ceil(Max((usrperturbation2(n)-usrperturbation1(n))/usrpertstep(n),1e-6)+0.000001);
      }
      int nit = 1;
      for (int n=1; n<=nsteps.Nrows(); n++) { nit *= MISCMATHS::round(nsteps(n)); }
      ColumnVector nvec;
      nvec = nsteps * 0.0f;
      for (int n=1; n<=nit; n++) {
	params = params0;  // start with the unperturbed params
	// the following increments the step pointer (to simulate n nested loops)
	int step=1;
	for (int m=1; m<=nsteps.Nrows(); m++) {
	  nvec(m)+=step;
	  if (MISCMATHS::round(nvec(m))<MISCMATHS::round(nsteps(m))) { step=0; } 
	  else { nvec(m)=0.0f; }
	}
	usrperturbation = usrperturbation1 + SP(nvec,usrpertstep);
	// use the elementwise product to produce the perturbation
	if (usrperturbrelative) {
	  params += SP(usrperturbation,delta); // rel
	} else {
	  params += usrperturbation; // abs
	}
	vector2affine(params,12,matresult);
	
	float costval=0.0;
	costval = measure_cost(matresult,dof);
	reshape(reshaped,matresult,1,16);
	rowresult(1) = costval;
	rowresult.SubMatrix(1,1,2,17) = reshaped;
	// store result
	stdresultmat->push_back(rowresult);
      }
    }
}


void usroptimise(MatVecPtr stdresultmat, 
		 MatVecPtr usrmatptr, 
		 unsigned int usrrow1, unsigned int usrrow2, int usrdof, 
		 ColumnVector& usrperturbation, bool usrperturbrelative,
		 int usrmaxitn)
{
  Tracer tr("usroptimise");
  // OPTIMISE
  Matrix delta, perturbmask, matresult;
  set_perturbations(delta,perturbmask);
  int dof = Min(globaloptions::get().dof,usrdof);
  ColumnVector params(12);
  RowVector rowresult(17);
  for (unsigned int crow=usrrow1; crow<=Min(usrrow2,usrmatptr->size()); crow++)
    {
      // the pre-optimised case with perturbations
      Matrix reshaped = (*usrmatptr)[crow-1].SubMatrix(1,1,2,17);
      reshape(matresult,reshaped,4,4);
      affmat2vector(matresult,12,params);
      // use the elementwise product to produce the perturbation
      if (usrperturbrelative) {
	params += SP(usrperturbation,delta); // rel
      } else {
	params += usrperturbation; // abs
      }
      vector2affine(params,12,matresult);
      
      float costval=0.0;
      if (globaloptions::get().usrsubset) {
	optimise_strategy0(matresult,costval,usrmaxitn);
      } else {
	optimise_strategy1(matresult,costval,dof,usrmaxitn);
      }
      reshape(reshaped,matresult,1,16);
      rowresult(1) = costval;
      rowresult.SubMatrix(1,1,2,17) = reshaped;
      // store result
      stdresultmat->push_back(rowresult);
    }
}



void usrsetscale(float usrscale, bool usrforce,
		 volume<float>& testvol, volume<float>& refvol, 
		 volume<float>& refvol_2, volume<float>& refvol_4, 
		 volume<float>& refvol_8) {
  // SETSCALE (int usrscale = 8,4,2,1)
  // testvol must be passed in as a static storage is needed so that
  //  the object pointed to in globaloptions::get().impair does not go 
  //  out of scope
  Tracer tr("usrsetscale");
  float scale = usrscale;
  bool forcescale=usrforce;
  if (globaloptions::get().force_scaling) forcescale=true;
  Costfn *globalpair=0;
  if (usrscale > 0.0) globaloptions::get().requestedscale = usrscale;

  if (globaloptions::get().debug) {
    cerr << "MJ DEBUG OUTPUT: scale = " << scale << endl;
    cerr << "MJ DEBUG OUTPUT: min_sampling = " << globaloptions::get().min_sampling 
	 << endl;
    print_volume_info(testvol,"testvol DEBUG");
  }

  // refresh the testvol (get rid of previous blurred version)
  get_testvol(testvol);

  if ( (forcescale) || (globaloptions::get().min_sampling<=1.25 * scale) ) {  // MJ NOTE: SHOULD THE SCALE BE FORCED TO BE THE NEXT HIGHEST SENSIBLE ONE TO STOP IT STARTING AT 8MM (THE DEFAULT) AND THEN NOT GOING TO 1MM (BUT INSTEAD MAKING IT GO TO 2MM?)
    globaloptions::get().lastsampling = scale;
    // blur test volume to correct scale
    volume<float> testvolnew;
    filter_image(testvolnew,testvol,global_testweight,scale,
		 globaloptions::get().useweights,filter_blur);
    if (globaloptions::get().useweights) {
      filter_weight(global_testweight,global_testweight,scale,filter_blur);
    }
    testvol = testvolnew;
    
    // select correct refvol
    volume<float> *refvolnew=0;
    if (fabs(scale-8.0)<0.01) {
      refvolnew = &refvol_8;
      if (globaloptions::get().useweights) {
	global_refweight = global_refweight8;
      }
    } else if (fabs(scale-4.0)<0.01) {
      refvolnew = &refvol_4;
      if (globaloptions::get().useweights) {
	global_refweight = global_refweight4;
      }
    } else if (fabs(scale-2.0)<0.01) {
      refvolnew = &refvol_2;
      if (globaloptions::get().useweights) {
	global_refweight = global_refweight2;
      }
    } else if (global_scale1OK && (fabs(scale-1.0)<0.01) ) {
      refvolnew = &refvol;
      if (globaloptions::get().useweights) {
	global_refweight = global_refweight1;
      }
    } else {
      // make a new refvol at this requested scale
      global_scale1OK = false;
      volume<float> tmpvol;
      get_refvol(tmpvol);  // gets raw refvol and global_refweight
      resample_refvol(tmpvol,scale);
      refvol = tmpvol;  // destroy base refvol!
      refvolnew = &refvol;
      if (globaloptions::get().useweights) {
	resample_refvol(global_refweight,scale);
	global_refweight1 = global_refweight;  // destroy global_refweight
      }
    }
    if (globaloptions::get().useweights) {
      globalpair = new Costfn(*refvolnew,testvol,
			      global_refweight,global_testweight);
    } else {
      globalpair = new Costfn(*refvolnew,testvol);
    }
    
    setup_costfn(globalpair,globaloptions::get().currentcostfn,
		 int(globaloptions::get().no_bins/scale),
		 globaloptions::get().smoothsize,globaloptions::get().fuzzyfrac);
    if (globaloptions::get().verbose>=3) {
      if (globaloptions::get().impair) {
	cout << "Previous scale used " << globaloptions::get().impair->count()
	     << " cost function evaluations" << endl;
      }
    }
    if (globaloptions::get().impair)  delete globaloptions::get().impair;
    globaloptions::get().impair = globalpair;
  }
}

void interpretcommand(const string& comline, bool& skip,
		      volume<float>& testvol, volume<float>& refvol, 
		      volume<float>& refvol_2, volume<float>& refvol_4, 
		      volume<float>& refvol_8)
{
  Tracer tr("interpretcommand");
  std::vector<string> words(0);
  parseline(comline,words);
  if (words.size()<1) return;
  if ((words[0])[0] == '#') return;  // comment line

  if (skip) {
    skip=false;
    return;
  }

  if (words[0]=="copy") {
    // COPY
    if (words.size()<3) {
      cerr << "Wrong number of args to COPY" << endl;
      exit(-1);
    }
    MatVecPtr src, dest;
    int d1, d2, d3, d4;
    parsematname(words[1],src,d1,d2);
    parsematname(words[2],dest,d3,d4);
    usrcopy(src,dest,d1,d2);
  } else if (words[0]=="clear") {
    // CLEAR
    if (words.size()<2) {
      cerr << "Wrong number of args to CLEAR" << endl;
      exit(-1);
    }
    MatVecPtr src;
    int d1, d2;
    parsematname(words[1],src,d1,d2);
    usrclear(src);
  } else if (words[0]=="exit") {
    // EXIT
    exit(0);
  } else if (words[0]=="print") {
    // PRINT
    if (words.size()<2) {
      cerr << "Wrong number of args to PRINT" << endl;
      exit(-1);
    }
    MatVecPtr src;
    int d1, d2;
    parsematname(words[1],src,d1,d2);
    usrprint(src,d1,d2);
  } else if (words[0]=="save") {
    // SAVE
    if (words.size()<3) {
      cerr << "Wrong number of args to SAVE" << endl;
      exit(-1);
    }
    MatVecPtr src;
    int d1, d2;
    parsematname(words[1],src,d1,d2);
    usrsave(words[2],src,d1,d2);
  } else if (words[0]=="read") {
    // READ
    if (words.size()<3) {
      cerr << "Wrong number of args to READ" << endl;
      exit(-1);
    }
    MatVecPtr src;
    int d1, d2;
    parsematname(words[1],src,d1,d2);
    usrread(words[2],src);
  } else if (words[0]=="printparams") {
    // PRINTPARAMS
    if (words.size()<2) {
      cerr << "Wrong number of args to PRINTPARAMS" << endl;
      exit(-1);
    }
    MatVecPtr src;
    int d1, d2;
    parsematname(words[1],src,d1,d2);
    usrprintparams(src,d1,d2);
  } else if (words[0]=="saveparams") {
    // SAVEPARAMS
    if (words.size()<3) {
      cerr << "Wrong number of args to SAVEPARAMS" << endl;
      exit(-1);
    }
    MatVecPtr src;
    int d1, d2;
    parsematname(words[1],src,d1,d2);
    usrsaveparams(words[2],src,d1,d2);
  } else if (words[0]=="readparams") {
    // READPARAMS
    if (words.size()<3) {
      cerr << "Wrong number of args to READPARAMS" << endl;
      exit(-1);
    }
    MatVecPtr src;
    int d1, d2;
    parsematname(words[1],src,d1,d2);
    usrreadparams(words[2],src);
  } else if (words[0]=="dualsort") {
    // DUALSORT
    if (words.size()<3) {
      cerr << "Wrong number of args to DUALSORT" << endl;
      exit(-1);
    }
    MatVecPtr mat1, mat2;
    int d1, d2, d3, d4;
    parsematname(words[1],mat1,d1,d2);
    parsematname(words[2],mat2,d3,d4);
    usrdualsort(mat1,mat2);
  } else if (words[0]=="sort") {
    // SORT
    if (words.size()<2) {
      cerr << "Wrong number of args to SORT" << endl;
      exit(-1);
    }
    MatVecPtr src;
    int d1, d2;
    parsematname(words[1],src,d1,d2);
    usrsort(src);
  } else if (words[0]=="search") {
    // SEARCH
    int usrdof=12;
    if (words.size()>2) {
      cerr << "Wrong number of args to SEARCH" << endl;
      exit(-1);
    }
    if (words.size()==2) {
      setscalarvariable(words[1],usrdof);
    }
    usrsearch(&globaloptions::get().searchoptmat,
	      &globaloptions::get().preoptsearchmat,usrdof);
  } else if (words[0]=="optimise") {
    // OPTIMISE
    if (words.size()<5) {
      cerr << "Wrong number of args to OPTIMISE" << endl;
      exit(-1);
    }
    int usrdof=12, usrmaxitn=4, usrrow1=1, usrrow2=999999;
    ColumnVector usrperturbation(12);
    usrperturbation = 0.0;
    MatVecPtr usrdefmat;
    bool usrperturbrelative = true;
    setscalarvariable(words[1],usrdof);
    parsematname(words[2],usrdefmat,usrrow1,usrrow2);
    unsigned int wordno=3;
    float tempparam=0.0;
    while ((wordno<words.size()) 
	   && (words[wordno]!="rel") && (words[wordno]!="abs")) {
      setscalarvariable(words[wordno],tempparam);
      usrperturbation(wordno-2) = tempparam;
      wordno++;
    }
    if (wordno<words.size()) {
      if (words[wordno]=="abs")  usrperturbrelative = false;
      wordno++;
    }
    if (wordno<words.size()) {
      setscalarvariable(words[wordno],usrmaxitn);
    }
    usroptimise(&(globaloptions::get().usrmat[0]),
		usrdefmat,usrrow1,usrrow2,usrdof, 
		usrperturbation,usrperturbrelative,usrmaxitn);
  } else if (words[0]=="measurecost") {
    // MEASURECOST
    if (words.size()<5) {
      cerr << "Wrong number of args to MEASURECOST" << endl;
      exit(-1);
    }
    int usrdof=12, usrrow1=1, usrrow2=999999;
    ColumnVector usrperturbation(12);
    usrperturbation = 0.0;
    MatVecPtr usrdefmat;
    bool usrperturbrelative = true;
    setscalarvariable(words[1],usrdof);
    parsematname(words[2],usrdefmat,usrrow1,usrrow2);
    unsigned int wordno=3;
    float tempparam=0.0;
    while ((wordno<words.size()) 
	   && (words[wordno]!="rel") && (words[wordno]!="abs")) {
      setscalarvariable(words[wordno],tempparam);
      usrperturbation(wordno-2) = tempparam;
      wordno++;
    }
    if (wordno<words.size()) {
      if (words[wordno]=="abs")  usrperturbrelative = false;
      wordno++;
    }
    usrmeasurecost(&(globaloptions::get().usrmat[0]),
		   usrdefmat,usrrow1,usrrow2,usrdof, 
		   usrperturbation,usrperturbrelative);
  } else if (words[0]=="gridmeasurecost") {
    // GRIDMEASURECOST
    if (words.size()<7) {
      cerr << "Wrong number of args to GRIDMEASURECOST" << endl;
      exit(-1);
    }
    int usrdof=12, usrrow1=1, usrrow2=999999;
    ColumnVector usrperturbation1(12), usrperturbation2(12), usrpertstep(12);
    usrperturbation1 = 0.0;
    usrperturbation2 = 0.0;
    usrpertstep = 0.0;
    MatVecPtr usrdefmat;
    bool usrperturbrelative = true;
    setscalarvariable(words[1],usrdof);
    parsematname(words[2],usrdefmat,usrrow1,usrrow2);
    unsigned int wordno=3;
    float tempparam=0.0;
    int compidx=1;
    while ((wordno+2<words.size()) 
	   && (words[wordno]!="rel") && (words[wordno]!="abs")) {
      setscalarvariable(words[wordno],tempparam);
      usrperturbation1(compidx) = tempparam;
      wordno++;
      setscalarvariable(words[wordno],tempparam);
      usrpertstep(compidx) = tempparam;
      wordno++;
      setscalarvariable(words[wordno],tempparam);
      usrperturbation2(compidx) = tempparam;
      wordno++;
      compidx++;
    }
    while (wordno<words.size()) {
      if (words[wordno]=="abs")  usrperturbrelative = false;
      wordno++;
    }
    usrgridmeasurecost(&(globaloptions::get().usrmat[0]),
		       usrdefmat,usrrow1,usrrow2,usrdof, 
		       usrperturbation1,usrpertstep,
		       usrperturbation2,usrperturbrelative);
  } else if (words[0]=="aligncog") {
    // ALIGNCOG
    if (words.size()<2) {
      cerr << "Wrong number of args to ALIGNCOG" << endl;
      exit(-1);
    }
    MatVecPtr src;
    int d1, d2;
    parsematname(words[1],src,d1,d2);
    usraligncog(&(globaloptions::get().usrmat[0]),src,d1,d2);
  } else if (words[0]=="alignpaxes") {
    // ALIGNPAXES
    if (words.size()<2) {
      cerr << "Wrong number of args to ALIGNPAXES" << endl;
      exit(-1);
    }
    MatVecPtr src;
    int d1, d2;
    parsematname(words[1],src,d1,d2);
    usralignpaxes(&(globaloptions::get().usrmat[0]),src,d1,d2);
  } else if (words[0]=="setscale") {
    // SETSCALE
    if (words.size()<2) {
      cerr << "Wrong number of args to SETSCALE" << endl;
      exit(-1);
    }
    bool usrforce=false;
    float usrscale;
    setscalarvariable(words[1],usrscale);
    if ((words.size()>=3) && (words[2]=="force")) usrforce=true;
    usrsetscale(usrscale,usrforce,testvol,refvol,refvol_2,refvol_4,refvol_8);
  } else if (words[0]=="setrow") {
    // SETROW
    if (words.size()<18) {
      cerr << "Wrong number of args to SETROW" << endl;
      exit(-1);
    }
    MatVecPtr src;
    int d1, d2;
    parsematname(words[1],src,d1,d2);
    usrsetrow(src,words);
  } else if (words[0]=="setrowqsform") {
    // SETROWQSFORM
    if (words.size()<2) {
      cerr << "Wrong number of args to SETROWQSFORM" << endl;
      exit(-1);
    }
    MatVecPtr src;
    int d1, d2;
    parsematname(words[1],src,d1,d2);
    usrsetrowqsform(src,refvol,testvol);
  } else if (words[0]=="setrowparams") {
    // SETROWPARAMS
    if (words.size()<14) {
      cerr << "Wrong number of args to SETROWPARAMS" << endl;
      exit(-1);
    }
    MatVecPtr src;
    int d1, d2;
    parsematname(words[1],src,d1,d2);
    usrsetrowparams(src,words);
  } else if (words[0]=="setoption") {
    // SETOPTION
    if (words.size()<3) {
      cerr << "Wrong number of args to SETOPTION" << endl;
      exit(-1);
    }
    usrsetoption(words);
  } else if (words[0]=="if") {
    // IF
    if (words.size()<4) {
      cerr << "Wrong number of args to IF" << endl;
      exit(-1);
    }
    float arg1, arg2;
    setscalarvariable(words[1],arg1);
    setscalarvariable(words[3],arg2);
    if (words[2]==">") {
      skip = !(arg1 > arg2);
    } else if (words[2]=="<") {
      skip = !(arg1 < arg2);
    } else if (words[2]=="==") {
      skip = !(arg1 == arg2);
    } else if (words[2]=="!=") {
      skip = !(arg1 != arg2);
    } else if (words[2]=="<=") {
      skip = !(arg1 <= arg2);
    } else if (words[2]==">=") {
      skip = !(arg1 >= arg2);
    } else {
      cerr << "Cannot recognise operator " << words[2] << " in IF statement\n";
      exit(-1);
    }
  } else {
    cerr << "Unrecognised command " << words[0] << endl;
    cerr << " ... ignoring" << endl;
    //cerr << "Quitting." << endl;
    //exit(-1);
  }
}


////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{
  Tracer tr("main");

  try {

    globaloptions::get().parse_command_line(argc, argv,version);

    if (!globaloptions::get().do_optimise) {   
      do_applyxfm();
    }
    // reset the basescale for images where voxels are quite different from the
    //   usual human brain size (this must be done before any volumes are read, 
    //   since part of the reading process uses the basescale for re-scaling)
    set_basescale(globaloptions::get().reffname,globaloptions::get().inputfname);

    // READ IN THE VOLUMES

    volume<float> refvol, testvol;
    get_refvol(refvol);
    get_testvol(testvol);
    set_initmat(refvol,testvol);

    if ( (refvol.sform_code()!=NIFTI_XFORM_UNKNOWN) && 
	 (testvol.sform_code()!=NIFTI_XFORM_UNKNOWN) ) {
      if (globaloptions::get().verbose>0) {
	cerr << "WARNING: Both reference and input images have an sform matrix set" << endl;
      }
    }

    if (globaloptions::get().verbose>0) {
      if (refvol.sform_code()!=NIFTI_XFORM_UNKNOWN) {
	cout << "The output image will use the sform from the reference image" << endl;    
      }
      if (testvol.sform_code()!=NIFTI_XFORM_UNKNOWN) {
	cout << "The output image will use the transformed sform from the input image" << endl;    
      }
    }

    if ( (globaloptions::get().verbose>0) || (globaloptions::get().printinit)) {
      cout << "Init Matrix = \n" << globaloptions::get().initmat << endl;
    }

    // Calculate quantities used to work out the correct resolution/sampling
    //  - especially important for very high-res volumes where a 1mm scale is too big

    float min_sampling_ref=1.0, min_sampling_test=1.0, min_sampling=1.0;
    min_sampling_ref = Min(refvol.xdim(),Min(refvol.ydim(),refvol.zdim()));
    min_sampling_test = Min(testvol.xdim(),Min(testvol.ydim(),testvol.zdim()));
    min_sampling = (float) ceil(Max(min_sampling_ref,min_sampling_test));
    if (!globaloptions::get().force_scaling) {
      // take the MAXIMUM of the user specified minimum and the estimated min
      if (globaloptions::get().min_sampling < min_sampling) 
	globaloptions::get().min_sampling = min_sampling;
    }
    
    if (globaloptions::get().verbose>=3) {
      cout << "CoG for refvol is:  " << refvol.cog("scaled_mm").t();
      cout << "CoG for testvol is:  " << testvol.cog("scaled_mm").t();
    }
  
    // CREATE THE VARIOUS SUB-SAMPLED REFERENCE VOLUMES FOR THE MULTI-SCALE

    // REFVOL RESAMPLING
    //    QUESTION: IS IT OK TO RECURSIVELY DEFINE WEIGHTS AND SUBSAMPLE LIKE THIS?
    //              OR WOULD DIRECT IMPLEMENTATION OF 4 AND 8 TIMES SUBSAMPLING BE 
    //              BETTER?  (SEP2010)
    volume<float> refvol_2, refvol_4, refvol_8;
    if (globaloptions::get().resample) {
      // set up subsampled volumes by factors of 2, 4 and 8
      if (globaloptions::get().verbose >= 2) 
	cout << "Subsampling the volumes" << endl;
      
      resample_refvol(refvol,globaloptions::get().min_sampling);
      filter_weight(global_refweight,global_refweight,
		    globaloptions::get().min_sampling,filter_resamp_blur);
      // the following tests enforce a maximum subsampling (i.e. for large voxel sizes, do not subsample as much as if the voxels are small)
      global_refweight1 = global_refweight;
      // SCALE 2
      if (globaloptions::get().min_sampling < 1.9) {
	filter_image(refvol_2,refvol,global_refweight,
		     globaloptions::get().useweights,filter_subsample_by_2);
	if (globaloptions::get().useweights) {
	  filter_weight(global_refweight2,global_refweight1,filter_subsample_by_2);
	}
      } else {
	refvol_2 = refvol;
	if (globaloptions::get().useweights) { global_refweight2 = global_refweight1; }
      }
      // SCALE 4
      if (globaloptions::get().min_sampling < 3.9) {
	filter_image(refvol_4,refvol_2,global_refweight2,
		     globaloptions::get().useweights,filter_subsample_by_2);
	if (globaloptions::get().useweights) {
	  filter_weight(global_refweight4,global_refweight2,filter_subsample_by_2);
	}
      } else {
	refvol_4 = refvol_2;
	if (globaloptions::get().useweights) { global_refweight4 = global_refweight2; }
      }
      // SCALE 8
      if (globaloptions::get().min_sampling < 7.9) {
	filter_image(refvol_8,refvol_4,global_refweight4,
		     globaloptions::get().useweights,filter_subsample_by_2);
	if (globaloptions::get().useweights) {
	  filter_weight(global_refweight8,global_refweight4,filter_subsample_by_2);
	}
      } else {
	refvol_8 = refvol_4;
	if (globaloptions::get().useweights) { global_refweight8 = global_refweight4; }
      }

    } else {
      // if no resampling chosen, then refvol is simply copied and nothing done to testvol
      refvol_8 = refvol;
      refvol_4 = refvol;
      refvol_2 = refvol;
      if (globaloptions::get().useweights) {
	global_refweight1 = global_refweight;
	global_refweight2 = global_refweight;
	global_refweight4 = global_refweight;
	global_refweight8 = global_refweight;
      }
    }


    // TESTVOL RESAMPLING
    if (globaloptions::get().resample) {
      volume<float> testvol_8;
      filter_image(testvol_8,testvol,global_testweight,8.0,
		   globaloptions::get().useweights,filter_blur);
      filter_weight(global_testweight,global_testweight,8.0,filter_blur);
    
      testvol = testvol_8;
    }

    if (globaloptions::get().debug) {
      save_volume(refvol_8,"refvol_8");
      save_volume(refvol_4,"refvol_4");
      save_volume(refvol_2,"refvol_2");
      save_volume(refvol,"refvol");
      save_volume(global_refweight1,"global_refweight1");
      save_volume(global_refweight2,"global_refweight2");
      save_volume(global_refweight4,"global_refweight4");
      save_volume(global_refweight8,"global_refweight8");
      save_volume(testvol,"testvol");
      save_volume(global_testweight,"testweight");
    }


    // set up image pair and global pointer, plus setup cost function params
    if (globaloptions::get().impair)  delete globaloptions::get().impair;
    global_refweight = global_refweight8;
    globaloptions::get().lastsampling = 8;
    if (globaloptions::get().useweights) {
      globaloptions::get().impair  = new Costfn(refvol_8,testvol,
						global_refweight,
						global_testweight);
    } else {
      globaloptions::get().impair = new Costfn(refvol_8,testvol);
    }

    globaloptions::get().currentcostfn = globaloptions::get().maincostfn;   
    setup_costfn(globaloptions::get().impair, globaloptions::get().currentcostfn,
		 globaloptions::get().no_bins/8,
		 globaloptions::get().smoothsize,globaloptions::get().fuzzyfrac);  
    if (globaloptions::get().verbose>=2) print_volume_info(testvol,"TESTVOL");
  

    Matrix matresult(4,4);
    ColumnVector params_8(12), param_tol(12);


    // PERFORM THE OPTIMISATION

    std::vector<string> schedulecoms(0);
    string comline;
    if (globaloptions::get().schedulefname.length()<1) {
      if (globaloptions::get().mode2D) {
	set2Ddefaultschedule(schedulecoms);
      } else {
	setdefaultschedule(schedulecoms);
      }
    } else {
      // open the schedule file
      ifstream schedulefile(globaloptions::get().schedulefname.c_str());
      if (!schedulefile) {
	cerr << "Could not open file" << globaloptions::get().schedulefname << endl;
	return -1;
      }
      while (!schedulefile.eof()) {
	getline(schedulefile,comline);
	schedulecoms.push_back(comline);
      }
      schedulefile.close();
    }

    // interpret each line in the schedule command vector
    bool skip=false;
    for (unsigned int i=0; i<schedulecoms.size(); i++) {
      comline = schedulecoms[i];
      if (globaloptions::get().verbose>=1) {
	cout << " >> " << comline << endl;
      }
      interpretcommand(comline,skip,testvol,refvol,refvol_2,refvol_4,refvol_8);
    }

    if (globaloptions::get().debug) {  // run this to save out any cost function debug info
      cerr << "Final DEBUG call in FLIRT" << endl;
      cerr << "Pointer to impair is " << globaloptions::get().impair << endl;
      globaloptions::get().impair->set_debug_mode(true);
      cerr << "Final DEBUG call in FLIRT 2" << endl;
      //cerr << "Reshaped matrix is:" << endl << reshaped << endl;
      cerr << "Scale is:" << globaloptions::get().requestedscale << endl;
      interpretcommand("measurecost 12 U:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 8",skip,testvol,refvol,refvol_2,refvol_4,refvol_8);
      // costfn(reshaped); // this currently causes an Abort!  Why?!?  No idea!!  :(
      cerr << "Final DEBUG call in FLIRT 3" << endl;
      globaloptions::get().impair->set_debug_mode(false);
      cerr << "Final DEBUG call in FLIRT 4" << endl;
    }
    
    // FINISHED OPTIMISATION - NOW GENERATE OUTPUTS

    // re-read the initial volume, and transform it by the optimised result

    Matrix reshaped;
    if (globaloptions::get().usrmat[0].size()>0) {
      reshaped = (globaloptions::get().usrmat[0])[0].SubMatrix(1,1,2,17);
      reshape(matresult,reshaped,4,4);

      // want unity basescale for transformed output
      float oldbasescale = globaloptions::get().basescale;
      globaloptions::get().basescale = 1.0;
      // make sure the old images don't get used
      if (globaloptions::get().impair) {
	delete globaloptions::get().impair;
	globaloptions::get().impair = NULL;
      }
      FLIRT_read_volume(testvol,globaloptions::get().inputfname);
      FLIRT_read_volume(refvol,globaloptions::get().reffname);
      if (globaloptions::get().verbose>=2) {
	print_volume_info(testvol,"testvol");
	print_volume_info(refvol,"refvol");
      }

      Matrix finalmat = matresult * globaloptions::get().initmat;
      finalmat(1,4) *= oldbasescale;
      finalmat(2,4) *= oldbasescale;
      finalmat(3,4) *= oldbasescale;
      if (globaloptions::get().verbose>=2) {
	cout << "Final transform matrix is:" << endl << finalmat << endl;
      }
      save_matrix_data(finalmat,testvol,refvol);
  
      // generate the outputvolume (not safe_save st -out overrides -nosave)
      if (globaloptions::get().outputfname.size()>0) {
	volume<float> newtestvol = refvol;
	float min_sampling_ref=1.0;
	min_sampling_ref = Min(refvol.xdim(),Min(refvol.ydim(),refvol.zdim()));
	if ((globaloptions::get().interpmethod != NearestNeighbour) &&
	    (globaloptions::get().interpblur)) {
	  filter_image(testvol,testvol,testvol,min_sampling_ref,
		       false,filter_blur);
	}
	final_transform(testvol,refvol,finalmat,newtestvol);      
	if (globaloptions::get().verbose>=2) {
	  print_volume_info(newtestvol,"Transformed testvol");
	}
	int outputdtype = output_dtype(newtestvol);
	newtestvol.setDisplayMaximumMinimum(0,0);
	save_volume_dtype(newtestvol,globaloptions::get().outputfname.c_str(),
			  outputdtype);
      }
      if ( (globaloptions::get().outputmatascii.size()<=0) ) {
	cout << endl << "Final result: " << endl << finalmat << endl;
      }
    } else {
      cerr << "Failed to calculate any transformation matrix" << endl;
    }

  }  
  catch(std::exception &e) {
    cerr << e.what() << endl;
  } 

  return(0);
}










