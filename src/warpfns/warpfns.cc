/*  warpfns.cc

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


// *** NOTE *** //
// MUST DEAL WITH THE PROBLEM OF EXTRAPOLATING VALUES FOR WARPS AS IT IS
// NOT GUARANTEED THAT THE FOV IS SUFFICIENT


#ifndef EXPOSE_TREACHEROUS
#define EXPOSE_TREACHEROUS           // To allow us to use .sampling_mat()
#endif

#include "newmatio.h"
#include "newmat.h"
#include "newimage/newimage.h"
#include "basisfield/basisfield.h"
#include "basisfield/splinefield.h"
#include "basisfield/dctfield.h"
#include "warpfns.h"

using namespace NEWMAT;

namespace NEWIMAGE {

////////////////////////////////////////////////////////////////////////////

int affine2warp(const Matrix& affmat, volume4D<float>& warpvol,
		const volume<float>& outvol)
{
  // all warps in absolute format
  if (outvol.nvoxels() <= 0) {
    cerr << "Cannot do affine2warp as outvol has no size" << endl;
    return -1;
  }
  warpvol.reinitialize(outvol.xsize(),outvol.ysize(),outvol.zsize(),3);
  warpvol[0] = outvol;
  warpvol[1] = outvol;
  warpvol[2] = outvol;

  ColumnVector xin(4), xout(4);
  xin(4) = 1.0;  xout(4)=1.0;

  for (int z=outvol.minz(); z<=outvol.maxz(); z++) {
    for (int y=outvol.miny(); y<=outvol.maxy(); y++) {
      for (int x=outvol.minx(); x<=outvol.maxx(); x++) {
	//   convert x,y,z to mm coords (xout)
	xout(1) = x;  xout(2) = y;  xout(3) = z;
	xout = outvol.sampling_mat() * xout;
	xin = affmat.i() * xout;
	// use the mm coordinates to store the results
	warpvol[0](x,y,z) = xin(1);
	warpvol[1](x,y,z) = xin(2);
	warpvol[2](x,y,z) = xin(3);
      }
    }
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////

int calc_dir(const string& shiftdir, int& dir, int& sign)
{
  // dir is 0,1,2 for x,y,z  :  sign is +/- 1 for x/x-  etc..
  if (shiftdir=="x") {
    dir=0;  sign=1;
  } else if (shiftdir=="y") {
    dir=1;  sign=1;
  } else if (shiftdir=="z") {
    dir=2;  sign=1;
  } else if (shiftdir=="x-") {
    dir=0;  sign=-1;
  } else if (shiftdir=="y-") {
    dir=1;  sign=-1;
  } else if (shiftdir=="z-") {
    dir=2;  sign=-1;
  } else {
    cerr << "Cannot interpret shift direction = " << shiftdir << endl;
    return -1;
  }
  return 0;
}


int shift2warp(const volume<float>& shiftmap, 
	       volume4D<float>& warp, const string& shiftdir)
{
  affine2warp(IdentityMatrix(4),warp,shiftmap);  // use shiftmap as refvol (set size)
  int dir, sign;
  calc_dir(shiftdir,dir,sign);
  float voxdim = shiftmap.sampling_mat()(dir+1,dir+1);

  for (int z=shiftmap.minz(); z<=shiftmap.maxz(); z++) {
    for (int y=shiftmap.miny(); y<=shiftmap.maxy(); y++) {
      for (int x=shiftmap.minx(); x<=shiftmap.maxx(); x++) {
	// get amount of shift in mm
	float shift = shiftmap(x,y,z) * voxdim * sign;
	warp[dir](x,y,z) += shift;
      }
    }
  }
  return 0;
}


////////////////////////////////////////////////////////////////////////////

int convertwarp_rel2abs(volume4D<float>& warpvol)
{
  // conversion is: w(x) = x + u(x)  (all in mm)
  for (int z=0; z<warpvol.zsize(); z++) {
    for (int y=0; y<warpvol.ysize(); y++) {
      for (int x=0; x<warpvol.xsize(); x++) {
	warpvol(x,y,z,0) += x*warpvol.xdim();
	warpvol(x,y,z,1) += y*warpvol.ydim();
	warpvol(x,y,z,2) += z*warpvol.zdim();
      }
    }
  }
  return 0;
}

int convertwarp_abs2rel(volume4D<float>& warpvol)
{
  // conversion is: w(x) = x + u(x)  (all in mm)
  for (int z=0; z<warpvol.zsize(); z++) {
    for (int y=0; y<warpvol.ysize(); y++) {
      for (int x=0; x<warpvol.xsize(); x++) {
	warpvol(x,y,z,0) -= x*warpvol.xdim();
	warpvol(x,y,z,1) -= y*warpvol.ydim();
	warpvol(x,y,z,2) -= z*warpvol.zdim();
      }
    }
  }
  return 0;
}


bool is_abs_convention(const volume4D<float>& inwarpvol)
{
  volume4D<float>&  warpvol = const_cast<volume4D<float>& >(inwarpvol);

  // try to determine if the warp is stored in absolute (vs relative) convention
  bool abs_warp=false;
  float stddev0 = warpvol[0].stddev()+warpvol[1].stddev()+warpvol[2].stddev();
  convertwarp_abs2rel(warpvol);
  float stddev1 = warpvol[0].stddev()+warpvol[1].stddev()+warpvol[2].stddev();
  // restore to the original form
  convertwarp_rel2abs(warpvol);
  // assume that relative warp always has less stddev
  if (stddev0>stddev1) {
    // the initial one (greater stddev) was absolute
    abs_warp=true;
  } else {
    // the initial one was relative
    abs_warp=false;
  }
  return abs_warp;
}


////////////////////////////////////////////////////////////////////////////


Matrix best_fit_aff(const volume4D<float>& warp)
{
  Matrix XX(4,4), XY(4,3);
  double xx=0, yy=0, zz=0, xy=0, xz=0, yz=0, xt=0, yt=0, zt=0, x0, y0, z0;
  double xxtot=0, yytot=0, zztot=0, xytot=0, xztot=0, yztot=0, xtot=0, ytot=0, ztot=0;
  long int n=0, nlim, ntot=0;
  nlim = (long int) sqrt((double) warp[0].nvoxels());
  if (nlim<1000) nlim=1000;
  for (int z=warp.minz(); z<=warp.maxz(); z++) {
    for (int y=warp.miny(); y<=warp.maxy(); y++) {
      for (int x=warp.minx(); x<=warp.maxx(); x++) {
	x0=x*warp.xdim(); y0=y*warp.ydim(); z0=z*warp.zdim();
	xx += x0*x0;
	yy += y0*y0;
	zz += z0*z0;
	xy += x0*y0;
	xz += x0*z0;
	yz += y0*z0;
	xt += x0;
	yt += y0;
	zt += z0;
	n++;
	if (n>nlim) {
	  ntot+=n; xxtot+=xx; yytot+=yy; zztot+=zz; xytot+=xy; xztot+=xz; yztot+=yz; 
	  xtot+=xt; ytot+=yt; ztot+=zt;
	  n=0; xx=0; yy=0; zz=0; xy=0; xz=0; yz=0; xt=0; yt=0; zt=0;
	}
      }
    }
  }
  ntot+=n; xxtot+=xx; yytot+=yy; zztot+=zz; xytot+=xy; xztot+=xz; yztot+=yz;
  xtot+=xt; ytot+=yt; ztot+=zt;
  XX(1,1)=xxtot/ntot;
  XX(2,2)=yytot/ntot;
  XX(3,3)=zztot/ntot;
  XX(4,4)=1;
  XX(1,2)=xytot/ntot;
  XX(2,1)=XX(1,2);
  XX(1,3)=xztot/ntot;
  XX(3,1)=XX(1,3);
  XX(2,3)=yztot/ntot;
  XX(3,2)=XX(2,3);
  XX(1,4)=xtot/ntot;
  XX(4,1)=XX(1,4);
  XX(2,4)=ytot/ntot;
  XX(4,2)=XX(2,4);
  XX(3,4)=ztot/ntot;
  XX(4,3)=XX(3,4);

  for (int m=0; m<3; m++) {
    double vxtot=0, vytot=0, vztot=0, vtot=0, vx=0, vy=0, vz=0, vt=0, val;
    for (int z=warp.minz(); z<=warp.maxz(); z++) {
      for (int y=warp.miny(); y<=warp.maxy(); y++) {
	for (int x=warp.minx(); x<=warp.maxx(); x++) {
	  val = (double) warp(x,y,z,m);
	  x0=x*warp.xdim(); y0=y*warp.ydim(); z0=z*warp.zdim();
	  vx += val*x0;
	  vy += val*y0;
	  vz += val*z0;
	  vt += val;
	  if (n>nlim) {
	    vxtot+=vx; vytot+=vy; vztot+=vz; vtot+=vt;
	    n=0; vx=0; vy=0; vz=0; vt=0;
	  }
	}
      }
    }
    vxtot+=vx; vytot+=vy; vztot+=vz; vtot+=vt;
    XY(1,m+1)=vxtot/ntot;
    XY(2,m+1)=vytot/ntot;
    XY(3,m+1)=vztot/ntot;
    XY(4,m+1)=vtot/ntot;
  }
  
  Matrix Beta;
  Beta = pinv(XX) * XY;
  Matrix aff(4,4);
  for (int r=1; r<=3; r++) {
    for (int n=1; n<=4; n++) {
      aff(r,n) = Beta(n,r);
    }
  }
  aff(4,1)=0;  aff(4,2)=0;  aff(4,3)=0;  aff(4,4)=1;

  return aff;
}


int concat_warps(const volume4D<float>& prewarp, 
                 const volume4D<float>& postwarp,
		 volume4D<float>& totalwarp)
{
  // all warps are in absolute convention
  totalwarp = postwarp;  // set size
  totalwarp = 0.0;
  ColumnVector xmid(4), xmid_vox(4), xpre(4);
  Matrix extrap_aff;
  extrap_aff = best_fit_aff(prewarp);
  xmid(4) = 1.0; xmid_vox(4)=1.0; xpre(4)=1.0;
  for (int z=postwarp.minz(); z<=postwarp.maxz(); z++) {
    for (int y=postwarp.miny(); y<=postwarp.maxy(); y++) {
      for (int x=postwarp.minx(); x<=postwarp.maxx(); x++) {
	xmid(1) = postwarp[0](x,y,z);
	xmid(2) = postwarp[1](x,y,z);
	xmid(3) = postwarp[2](x,y,z);
	// convert xmid from mm to voxels (of prewarp image)
	xmid_vox = prewarp[0].sampling_mat().i() * xmid;
	if (prewarp.in_bounds((float) xmid_vox(1),(float) xmid_vox(2),(float) xmid_vox(3))) {
	  // look up the coordinates in prewarp
	  xpre(1) = prewarp[0].interpolate(xmid_vox(1),xmid_vox(2),xmid_vox(3));
	  xpre(2) = prewarp[1].interpolate(xmid_vox(1),xmid_vox(2),xmid_vox(3));
	  xpre(3) = prewarp[2].interpolate(xmid_vox(1),xmid_vox(2),xmid_vox(3));
	} else {
	  // alternative interpolation for warp fields
	  xpre = extrap_aff * xmid;
	}
	// set these mm coordinates as the result
	totalwarp[0](x,y,z) = xpre(1);
	totalwarp[1](x,y,z) = xpre(2);
	totalwarp[2](x,y,z) = xpre(3);
      }
    }
  }
  return 0;
}

volume<float> calc_sigloss(volume4D<float>& lrgrad, float te, float gammabar)
{
  // input gradient is in units of rad/s/voxel
  float gbarte_2 = gammabar * te / 2.0;
  volume<float> sigloss = lrgrad[0] * 0.0f;
  for (int z=lrgrad.minz(); z<=lrgrad.maxz(); z++) {
    for (int y=lrgrad.miny(); y<=lrgrad.maxy(); y++) {
      for (int x=lrgrad.minx(); x<=lrgrad.maxx(); x++) {
	float sincl = Sinc(gbarte_2*lrgrad[0](x,y,z));
	float sincr = Sinc(gbarte_2*lrgrad[1](x,y,z));
	float thetal = M_PI*gbarte_2*lrgrad[0](x,y,z);
	float thetar = M_PI*gbarte_2*lrgrad[1](x,y,z);
	float sigloss_re = 0.5 * ( sincl*cos(thetal) + sincr*cos(thetar) );
	float sigloss_im = 0.5 * ( sincl*sin(thetal) + sincr*sin(thetar) );
	sigloss(x,y,z) = sqrt( Sqr(sigloss_re) + Sqr(sigloss_im) );
      }
    }
  }
  return sigloss;
}


// topology preservation code

void jacobian_check(volume4D<float>& jvol,
		    ColumnVector& jacobian_stats, 
		    const volume4D<float>& warp,
		    float minJ, float maxJ, bool use_vol)
{
  // set up jacobian stats to contain: min, max, num < minJ, num > maxJ
  if (jacobian_stats.Nrows()!=4) { jacobian_stats.ReSize(4); }
  jacobian_stats = 0.0;
  jacobian_stats(1)=1.0; jacobian_stats(2)=1.0; 
  if (use_vol) {
    if ((jvol.tsize()!=8) || !samesize(jvol[0],warp[0])) {
      jvol = warp;  // set up all the right properties
      jvol = 0.0f; 
      for (int n=1; n<=5; n++) { jvol.addvolume(jvol[0]); }
    }
  }
  float Jfff, Jbff, Jfbf, Jffb, Jbbf, Jbfb, Jfbb, Jbbb;
  float wx000=0,wx001=0,wx010=0,wx011=0,wx100=0,wx101=0,wx110=0,wx111=0;
  float wy000=0,wy001=0,wy010=0,wy011=0,wy100=0,wy101=0,wy110=0,wy111=0;
  float wz000=0,wz001=0,wz010=0,wz011=0,wz100=0,wz101=0,wz110=0,wz111=0;
  float volscale=1.0/(warp.xdim() * warp.ydim() * warp.zdim());
  for (int z=warp.minz(); z<=warp.maxz()-1; z++) {
    for (int y=warp.miny(); y<=warp.maxy()-1; y++) {
      for (int x=warp.minx(); x<=warp.maxx()-1; x++) {
	warp[0].getneighbours(x,y,z,wx000,wx001,wx010,wx011,
			      wx100,wx101,wx110,wx111);
	warp[1].getneighbours(x,y,z,wy000,wy001,wy010,wy011,
			      wy100,wy101,wy110,wy111);
	warp[2].getneighbours(x,y,z,wz000,wz001,wz010,wz011,
			      wz100,wz101,wz110,wz111);
	Jfff = (wx100-wx000) *
	  ((wy010-wy000) * (wz001-wz000) - (wy001-wy000) * (wz010-wz000)) 
	  - (wx010-wx000) *
	  ((wy100-wy000) * (wz001-wz000) - (wy001-wy000) * (wz100-wz000)) 
	  + (wx001-wx000) *
	  ((wy100-wy000) * (wz010-wz000) - (wy010-wy000) * (wz100-wz000));
	Jfff *= volscale;
	if (Jfff<jacobian_stats(1)) { jacobian_stats(1)=Jfff; }
	if (Jfff>jacobian_stats(2)) { jacobian_stats(2)=Jfff; }
	if (Jfff<minJ) { jacobian_stats(3)+=1.0; }
	if (Jfff>maxJ) { jacobian_stats(4)+=1.0; }
	Jbff = (wx100-wx000) *
	  ((wy110-wy100) * (wz101-wz100) - (wy101-wy100) * (wz110-wz100)) 
	  - (wx110-wx100) *
	  ((wy100-wy000) * (wz101-wz100) - (wy101-wy100) * (wz100-wz000)) 
	  + (wx101-wx100) *
	  ((wy100-wy000) * (wz110-wz100) - (wy110-wy100) * (wz100-wz000));
	Jbff *= volscale;
	if (Jbff<jacobian_stats(1)) { jacobian_stats(1)=Jbff; }
	if (Jbff>jacobian_stats(2)) { jacobian_stats(2)=Jbff; }
	if (Jbff<minJ) { jacobian_stats(3)+=1.0; }
	if (Jbff>maxJ) { jacobian_stats(4)+=1.0; }
	Jfbf = (wx110-wx010) *
	  ((wy010-wy000) * (wz011-wz010) - (wy011-wy010) * (wz010-wz000)) 
	  - (wx010-wx000) *
	  ((wy110-wy010) * (wz011-wz010) - (wy011-wy010) * (wz110-wz010)) 
	  + (wx011-wx010) *
	  ((wy110-wy010) * (wz010-wz000) - (wy010-wy000) * (wz110-wz010));
	Jfbf *= volscale;
	if (Jfbf<jacobian_stats(1)) { jacobian_stats(1)=Jfbf; }
	if (Jfbf>jacobian_stats(2)) { jacobian_stats(2)=Jfbf; }
	if (Jfbf<minJ) { jacobian_stats(3)+=1.0; }
	if (Jfbf>maxJ) { jacobian_stats(4)+=1.0; }
	Jffb = (wx101-wx001) *
	  ((wy011-wy001) * (wz001-wz000) - (wy001-wy000) * (wz011-wz001)) 
	  - (wx011-wx001) *
	  ((wy101-wy001) * (wz001-wz000) - (wy001-wy000) * (wz101-wz001)) 
	  + (wx001-wx000) *
	  ((wy101-wy001) * (wz011-wz001) - (wy011-wy001) * (wz101-wz001));
	Jffb *= volscale;
	if (Jffb<jacobian_stats(1)) { jacobian_stats(1)=Jffb; }
	if (Jffb>jacobian_stats(2)) { jacobian_stats(2)=Jffb; }
	if (Jffb<minJ) { jacobian_stats(3)+=1.0; }
	if (Jffb>maxJ) { jacobian_stats(4)+=1.0; }
	Jfbb = (wx111-wx011) *
	  ((wy011-wy001) * (wz011-wz010) - (wy011-wy010) * (wz011-wz001)) 
	  - (wx011-wx001) *
	  ((wy111-wy011) * (wz011-wz010) - (wy011-wy010) * (wz111-wz011)) 
	  + (wx011-wx010) *
	  ((wy111-wy011) * (wz011-wz001) - (wy011-wy001) * (wz111-wz011));
	Jfbb *= volscale;
	if (Jfbb<jacobian_stats(1)) { jacobian_stats(1)=Jfbb; }
	if (Jfbb>jacobian_stats(2)) { jacobian_stats(2)=Jfbb; }
	if (Jfbb<minJ) { jacobian_stats(3)+=1.0; }
	if (Jfbb>maxJ) { jacobian_stats(4)+=1.0; }
	Jbfb = (wx101-wx001) *
	  ((wy111-wy101) * (wz101-wz100) - (wy101-wy100) * (wz111-wz101)) 
	  - (wx111-wx101) *
	  ((wy101-wy001) * (wz101-wz100) - (wy101-wy100) * (wz101-wz001)) 
	  + (wx101-wx100) *
	  ((wy101-wy001) * (wz111-wz101) - (wy111-wy101) * (wz101-wz001));
	Jbfb *= volscale;
	if (Jbfb<jacobian_stats(1)) { jacobian_stats(1)=Jbfb; }
	if (Jbfb>jacobian_stats(2)) { jacobian_stats(2)=Jbfb; }
	if (Jbfb<minJ) { jacobian_stats(3)+=1.0; }
	if (Jbfb>maxJ) { jacobian_stats(4)+=1.0; }
	Jbbf = (wx110-wx010) *
	  ((wy110-wy100) * (wz111-wz110) - (wy111-wy110) * (wz110-wz100)) 
	  - (wx110-wx100) *
	  ((wy110-wy010) * (wz111-wz110) - (wy111-wy110) * (wz110-wz010)) 
	  + (wx111-wx110) *
	  ((wy110-wy010) * (wz110-wz100) - (wy110-wy100) * (wz110-wz010));
	Jbbf *= volscale;
	if (Jbbf<jacobian_stats(1)) { jacobian_stats(1)=Jbbf; }
	if (Jbbf>jacobian_stats(2)) { jacobian_stats(2)=Jbbf; }
	if (Jbbf<minJ) { jacobian_stats(3)+=1.0; }
	if (Jbbf>maxJ) { jacobian_stats(4)+=1.0; }
	Jbbb = (wx111-wx011) *
	  ((wy111-wy101) * (wz111-wz110) - (wy111-wy110) * (wz111-wz101)) 
	  - (wx111-wx101) *
	  ((wy111-wy011) * (wz111-wz110) - (wy111-wy110) * (wz111-wz011)) 
	  + (wx111-wx110) *
	  ((wy111-wy011) * (wz111-wz101) - (wy111-wy101) * (wz111-wz011));
	Jbbb *= volscale;
	if (Jbbb<jacobian_stats(1)) { jacobian_stats(1)=Jbbb; }
	if (Jbbb>jacobian_stats(2)) { jacobian_stats(2)=Jbbb; }
	if (Jbbb<minJ) { jacobian_stats(3)+=1.0; }
	if (Jbbb>maxJ) { jacobian_stats(4)+=1.0; }
	if (use_vol) {
	  // the following must be consistent with get_jac_offset()
	  jvol(x,y,z,0) = Jfff;
	  jvol(x,y,z,1) = Jbff;
	  jvol(x,y,z,2) = Jfbf;
	  jvol(x,y,z,3) = Jffb;
	  jvol(x,y,z,4) = Jfbb;
	  jvol(x,y,z,5) = Jbfb;
	  jvol(x,y,z,6) = Jbbf;
	  jvol(x,y,z,7) = Jbbb;
	}
      }
    }
  }
}


volume4D<float> jacobian_check(ColumnVector& jacobian_stats, 
			       const volume4D<float>& warp,
			       float minJ, float maxJ)
{
  volume4D<float> jvol;
  jacobian_check(jvol,jacobian_stats,warp,minJ,maxJ,true);
  return jvol;
}


ColumnVector jacobian_quick_check(const volume4D<float>& warp,
				  float minJ, float maxJ)
{
  volume4D<float> dummy;
  ColumnVector jacobian_stats;
  jacobian_check(dummy,jacobian_stats,warp,minJ,maxJ,false);
  return jacobian_stats;
}

void grad_calc(volume4D<float>& gradvols, const volume4D<float>& warp)
{
  // returns gradients in the order: dx'/dx, dx'/dy, dx'/dz, dy'/dx, etc
  if ((gradvols.tsize()!=9) || !samesize(gradvols[0],warp[0])) {
    gradvols = warp;  // set up all the right properties
    gradvols=0.0f; 
    for (int n=1; n<=6; n++) { gradvols.addvolume(gradvols[0]); }
  }
  float dx=warp.xdim(), dy=warp.ydim(), dz=warp.zdim();
  float wx000=0,wx001=0,wx010=0,wx011=0,wx100=0,wx101=0,wx110=0,wx111=0;
  float wy000=0,wy001=0,wy010=0,wy011=0,wy100=0,wy101=0,wy110=0,wy111=0;
  float wz000=0,wz001=0,wz010=0,wz011=0,wz100=0,wz101=0,wz110=0,wz111=0;
  for (int z=warp.minz(); z<=warp.maxz(); z++) {
    for (int y=warp.miny(); y<=warp.maxy(); y++) {
      for (int x=warp.minx(); x<=warp.maxx(); x++) {
	if ((z<warp.maxz()) && (y<warp.maxy()) && (x<warp.maxx())) {
	  warp[0].getneighbours(x,y,z,wx000,wx001,wx010,wx011,
				wx100,wx101,wx110,wx111);
	  warp[1].getneighbours(x,y,z,wy000,wy001,wy010,wy011,
				wy100,wy101,wy110,wy111);
	  warp[2].getneighbours(x,y,z,wz000,wz001,wz010,wz011,
				wz100,wz101,wz110,wz111);
	  /*
	  gradvols[0](x,y,z) = (wx100-wx000)/dx;
	  gradvols[1](x,y,z) = (wx010-wx000)/dy;
	  gradvols[2](x,y,z) = (wx001-wx000)/dz;
	  gradvols[3](x,y,z) = (wy100-wy000)/dx;
	  gradvols[4](x,y,z) = (wy010-wy000)/dy;
	  gradvols[5](x,y,z) = (wy001-wy000)/dz;
	  gradvols[6](x,y,z) = (wz100-wz000)/dx;
	  gradvols[7](x,y,z) = (wz010-wz000)/dy;
	  gradvols[8](x,y,z) = (wz001-wz000)/dz;
	  */
	  gradvols[0](x,y,z) = (wx100-wx000)/dx;
	  gradvols[1](x,y,z) = (wx010-wx000)/dx;
	  gradvols[2](x,y,z) = (wx001-wx000)/dx;
	  gradvols[3](x,y,z) = (wy100-wy000)/dy;
	  gradvols[4](x,y,z) = (wy010-wy000)/dy;
	  gradvols[5](x,y,z) = (wy001-wy000)/dy;
	  gradvols[6](x,y,z) = (wz100-wz000)/dz;
	  gradvols[7](x,y,z) = (wz010-wz000)/dz;
	  gradvols[8](x,y,z) = (wz001-wz000)/dz;
	} else {
	  int x2=x+1,y2=y+1,z2=z+1;
	  if (x2>warp.maxx()) { x2 -= warp.maxx() + 1 - warp.minx(); }
	  if (y2>warp.maxy()) { y2 -= warp.maxy() + 1 - warp.miny(); }
	  if (z2>warp.maxz()) { z2 -= warp.maxz() + 1 - warp.minz(); }
	  /*
 	  gradvols[0](x,y,z) = (warp[0](x2,y,z) - warp[0](x,y,z))/dx;
	  gradvols[1](x,y,z) = (warp[0](x,y2,z) - warp[0](x,y,z))/dy;
	  gradvols[2](x,y,z) = (warp[0](x,y,z2) - warp[0](x,y,z))/dz;
	  gradvols[3](x,y,z) = (warp[1](x2,y,z) - warp[1](x,y,z))/dx;
	  gradvols[4](x,y,z) = (warp[1](x,y2,z) - warp[1](x,y,z))/dy;
	  gradvols[5](x,y,z) = (warp[1](x,y,z2) - warp[1](x,y,z))/dz;
	  gradvols[6](x,y,z) = (warp[2](x2,y,z) - warp[2](x,y,z))/dx;
	  gradvols[7](x,y,z) = (warp[2](x,y2,z) - warp[2](x,y,z))/dy;
	  gradvols[8](x,y,z) = (warp[2](x,y,z2) - warp[2](x,y,z))/dz;
	  */
	  gradvols[0](x,y,z) = (warp[0](x2,y,z) - warp[0](x,y,z))/dx;
	  gradvols[1](x,y,z) = (warp[0](x,y2,z) - warp[0](x,y,z))/dx;
	  gradvols[2](x,y,z) = (warp[0](x,y,z2) - warp[0](x,y,z))/dx;
	  gradvols[3](x,y,z) = (warp[1](x2,y,z) - warp[1](x,y,z))/dy;
	  gradvols[4](x,y,z) = (warp[1](x,y2,z) - warp[1](x,y,z))/dy;
	  gradvols[5](x,y,z) = (warp[1](x,y,z2) - warp[1](x,y,z))/dy;
	  gradvols[6](x,y,z) = (warp[2](x2,y,z) - warp[2](x,y,z))/dz;
	  gradvols[7](x,y,z) = (warp[2](x,y2,z) - warp[2](x,y,z))/dz;
	  gradvols[8](x,y,z) = (warp[2](x,y,z2) - warp[2](x,y,z))/dz;
	}
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
// I am using warpfns_ifft3 and warpfns_fft3 in lieu of the more general ifft3 and fft3 that
// are a part of the NEWIMAGE namespace. This is because by sacrificing "generality" I can
// make them ~30% faster, which I considered worth it since they are potentially called many
// times as part of projecting a deformation field onto a "Jacobian-limited surface".
//
// Note also that I have already tried _lots_ of ways to speed them up further by e.g. 
// using pointers into NEWMAT and NEWIMAGE data and by using NEWMAT::RealFFT and
// NEWMAT::RealFFTI for the first and last dimension transforms respectively. In neither
// case was the speed up sufficient to motivate the increased complexity of the code.
//
///////////////////////////////////////////////////////////////////////////////////////////////

void warpfns_ifft3(const NEWIMAGE::volume<float>&    rein,
                   const NEWIMAGE::volume<float>&    imin,
                   NEWIMAGE::volume<float>&          out)
{
  if (!samesize(rein,out)) out = rein;

  NEWIMAGE::volume<float>  retmp(rein.xsize(),rein.ysize(),rein.zsize());
  NEWIMAGE::volume<float>  imtmp(rein.xsize(),rein.ysize(),rein.zsize());
  out = 0.0; retmp = 0.0; imtmp = 0.0;
  // In x-direction
  NEWMAT::ColumnVector reinv(rein.xsize());
  NEWMAT::ColumnVector iminv(rein.xsize());
  NEWMAT::ColumnVector reoutv(rein.xsize());
  NEWMAT::ColumnVector imoutv(rein.xsize());
  for (int k=0; k<rein.zsize(); k++) {
    for (int j=0; j<rein.ysize(); j++) {
      for (int i=0; i<rein.xsize(); i++) {
        reinv(i+1) = rein(i,j,k);
        iminv(i+1) = imin(i,j,k);
      }
      FFTI(reinv,iminv,reoutv,imoutv);
      for (int i=0; i<rein.xsize(); i++) {
        retmp(i,j,k) = reoutv(i+1);
        imtmp(i,j,k) = imoutv(i+1);
      }
    }
  }
  // Y-direction
  reinv.ReSize(rein.ysize());
  iminv.ReSize(rein.ysize());
  reoutv.ReSize(rein.ysize());
  imoutv.ReSize(rein.ysize());
  for (int k=0; k<rein.zsize(); k++) {
    for (int i=0; i<rein.xsize(); i++) {
      for (int j=0; j<rein.ysize(); j++) {
        reinv(j+1) = retmp(i,j,k);
        iminv(j+1) = imtmp(i,j,k);
      }
      FFTI(reinv,iminv,reoutv,imoutv);
      for (int j=0; j<rein.ysize(); j++) {
        retmp(i,j,k) = reoutv(j+1);
        imtmp(i,j,k) = imoutv(j+1);
      }
    }
  }
  // And z-direction
  reinv.ReSize(rein.zsize());
  iminv.ReSize(rein.zsize());
  reoutv.ReSize(rein.zsize());
  imoutv.ReSize(rein.zsize());
  for (int j=0; j<rein.ysize(); j++) {
    for (int i=0; i<rein.xsize(); i++) {
      for (int k=0; k<rein.zsize(); k++) {
        reinv(k+1) = retmp(i,j,k);
        iminv(k+1) = imtmp(i,j,k);
      }
      FFTI(reinv,iminv,reoutv,imoutv);
      for (int k=0; k<rein.zsize(); k++) {
        out(i,j,k) = reoutv(k+1);
      }
    }
  }
  return;
}

void warpfns_fft3(const NEWIMAGE::volume<float>&    in,
                  NEWIMAGE::volume<float>&          reout,
                  NEWIMAGE::volume<float>&          imout)
{
  if (!samesize(in,reout)) reout = in;
  if (!samesize(in,imout)) imout = in;
  reout = 0.0;
  imout = 0.0;
  // cout << "x-direction" << endl;
  // In x-direction
  NEWMAT::ColumnVector reinv(in.xsize());
  NEWMAT::ColumnVector iminv(in.xsize());
  NEWMAT::ColumnVector reoutv(in.xsize());
  NEWMAT::ColumnVector imoutv(in.xsize());
  reinv = 0.0; iminv = 0.0;
  for (int k=0; k<in.zsize(); k++) {
    for (int j=0; j<in.ysize(); j++) {
      for (int i=0; i<in.xsize(); i++) {
        reinv(i+1) = in(i,j,k);
      }
      FFT(reinv,iminv,reoutv,imoutv);
      for (int i=0; i<in.xsize(); i++) {
        reout(i,j,k) = reoutv(i+1);
        imout(i,j,k) = imoutv(i+1);
      }
    }
  }
  // In y-direction
  // cout << "y-direction" << endl;
  reinv.ReSize(in.ysize());
  iminv.ReSize(in.ysize());
  reoutv.ReSize(in.ysize());
  imoutv.ReSize(in.ysize());
  reinv = 0.0; iminv = 0.0;
  for (int k=0; k<in.zsize(); k++) {
    for (int i=0; i<in.xsize(); i++) {
      for (int j=0; j<in.ysize(); j++) {
        reinv(j+1) = reout(i,j,k);
        iminv(j+1) = imout(i,j,k);
      }
      FFT(reinv,iminv,reoutv,imoutv);
      for (int j=0; j<in.ysize(); j++) {
        reout(i,j,k) = reoutv(j+1);
        imout(i,j,k) = imoutv(j+1);
      }
    }
  }
  // And z-direction
  // cout << "z-direction" << endl;
  reinv.ReSize(in.zsize());
  iminv.ReSize(in.zsize());
  reoutv.ReSize(in.zsize());
  imoutv.ReSize(in.zsize());
  reinv = 0.0; iminv = 0.0;
  for (int j=0; j<in.ysize(); j++) {
    for (int i=0; i<in.xsize(); i++) {
      for (int k=0; k<in.zsize(); k++) {
        reinv(k+1) = reout(i,j,k);
        iminv(k+1) = imout(i,j,k);
      }
      FFT(reinv,iminv,reoutv,imoutv);
      for (int k=0; k<in.zsize(); k++) {
        reout(i,j,k) = reoutv(k+1);
        imout(i,j,k) = imoutv(k+1);
      }
    }
  }
  return;
}

void integrate_gradient_field(volume4D<float>& newwarp, 
			      const volume4D<float>& grad,
			      float warpmeanx, float warpmeany, float warpmeanz)
{
  // enforces integrability constraints and returns the integrated grad field
  // Note that the mean of the newwarp will be equal to warpmean{x,y,z}
  //  pass in: oldwarp[0].mean(), oldwarp[1].mean(), oldwarp[2].mean()
  
  int Nx, Ny, Nz;
  Nx = grad.xsize();
  Ny = grad.ysize();
  Nz = grad.zsize();
  if ((newwarp.tsize()!=3) || !samesize(newwarp[0],grad[0])) {
    newwarp = grad;
    for (int n=8; n>2; n--) { newwarp.deletevolume(n); }
    newwarp = 0.0f;
  }
  volume4D<float> gradkre(newwarp), gradkim(newwarp);
  float gradkrealx, gradkimagx, gradkrealy, gradkimagy, gradkrealz, gradkimagz;
  // enforce things separately for gradients of warp[0], warp[1] and warp[2]
  for (int n=0; n<3; n++) {
    // take FFT of the x,y,z gradient fields (of warp[n])
    warpfns_fft3(grad[n*3+0],gradkre[0],gradkim[0]);
    warpfns_fft3(grad[n*3+1],gradkre[1],gradkim[1]);
    warpfns_fft3(grad[n*3+2],gradkre[2],gradkim[2]);
    for (int z=grad.minz(); z<=grad.maxz(); z++) {
      float argz=2.0*M_PI*z/Nz;
      float cosz=cos(argz), sinz=sin(argz);
      for (int y=grad.miny(); y<=grad.maxy(); y++) {
        float argy=2.0*M_PI*y/Ny;
        float cosy=cos(argy), siny=sin(argy);
	for (int x=grad.minx(); x<=grad.maxx(); x++) {
	  float argx=2.0*M_PI*x/Nx;
          float cosx=cos(argx), sinx=sin(argx);
	  float norm = 6.0 - 2.0*cosz - 2.0*cosy - 2.0*cosx;
	  gradkrealx = gradkre[0](x,y,z);
	  gradkimagx = gradkim[0](x,y,z);
	  gradkrealy = gradkre[1](x,y,z);
	  gradkimagy = gradkim[1](x,y,z);
	  gradkrealz = gradkre[2](x,y,z);
	  gradkimagz = gradkim[2](x,y,z);
	  float dotprodre = gradkrealx * (cosx-1) + gradkimagx * sinx;
	  float dotprodim = gradkimagx * (cosx-1) - gradkrealx * sinx;
	  dotprodre += gradkrealy * (cosy-1) + gradkimagy * siny;
	  dotprodim += gradkimagy * (cosy-1) - gradkrealy * siny;
	  dotprodre += gradkrealz * (cosz-1) + gradkimagz * sinz;
	  dotprodim += gradkimagz * (cosz-1) - gradkrealz * sinz;
	  // write back values into gradkre[0] and gradkim[0]
	  if (fabs(norm)>1e-12) {
	    gradkre[0](x,y,z) = dotprodre / norm;
	    gradkim[0](x,y,z) = dotprodim / norm;
	  } else {
	    gradkre[0](x,y,z) = 0;
	    gradkim[0](x,y,z) = 0;
	  }
	}
      }
    }
    // take IFFT to get the integrated gradient field
    volume<float> dummy1, dummy2;
    warpfns_ifft3(gradkre[0],gradkim[0],newwarp[n]);
  }
  // rescale newwarp by the voxel dimensions
  newwarp[0] *= grad.xdim();
  newwarp[1] *= grad.ydim();
  newwarp[2] *= grad.zdim();
  // adjust the mean values
  newwarp[0] += warpmeanx;
  newwarp[1] += warpmeany;
  newwarp[2] += warpmeanz;
}


void get_jac_offset(int jacnum, int* xoff, int* yoff, int* zoff)
{
  xoff[1]=0; yoff[1]=0; zoff[1]=0; 
  xoff[2]=0; yoff[2]=0; zoff[2]=0; 
  xoff[3]=0; yoff[3]=0; zoff[3]=0; 

  if (jacnum==0) {  // Jfff
    xoff[1]=0; yoff[1]=0; zoff[1]=0; 
    xoff[2]=0; yoff[2]=0; zoff[2]=0; 
    xoff[3]=0; yoff[3]=0; zoff[3]=0; 
  }
  if (jacnum==1) { // Jbff
    xoff[1]=0; yoff[1]=0; zoff[1]=0; 
    xoff[2]=1; yoff[2]=0; zoff[2]=0; 
    xoff[3]=1; yoff[3]=0; zoff[3]=0; 
  }
  if (jacnum==2) { // Jfbf
    xoff[1]=0; yoff[1]=1; zoff[1]=0; 
    xoff[2]=0; yoff[2]=0; zoff[2]=0; 
    xoff[3]=0; yoff[3]=1; zoff[3]=0; 
  }
  if (jacnum==3) { // Jffb
    xoff[1]=0; yoff[1]=0; zoff[1]=1; 
    xoff[2]=0; yoff[2]=0; zoff[2]=1; 
    xoff[3]=0; yoff[3]=0; zoff[3]=0; 
  }
  if (jacnum==4) { // Jfbb
    xoff[1]=0; yoff[1]=1; zoff[1]=1; 
    xoff[2]=0; yoff[2]=0; zoff[2]=1; 
    xoff[3]=0; yoff[3]=1; zoff[3]=0; 
  }
  if (jacnum==5) { // Jbfb
    xoff[1]=0; yoff[1]=0; zoff[1]=1; 
    xoff[2]=1; yoff[2]=0; zoff[2]=1; 
    xoff[3]=1; yoff[3]=0; zoff[3]=0; 
  }
  if (jacnum==6) { // Jbbf
    xoff[1]=0; yoff[1]=1; zoff[1]=0; 
    xoff[2]=1; yoff[2]=0; zoff[2]=0; 
    xoff[3]=1; yoff[3]=1; zoff[3]=0; 
  }
  if (jacnum==7) { // Jbbb
    xoff[1]=0; yoff[1]=1; zoff[1]=1; 
    xoff[2]=1; yoff[2]=0; zoff[2]=1; 
    xoff[3]=1; yoff[3]=1; zoff[3]=0; 
  }
}

void limit_grad(volume4D<float>& grad, const volume4D<float>& jvol, 
		float minJ, float maxJ, const Matrix& p_initaffmat)
{
  Matrix J(3,3), Jnew(3,3), J0(3,3), initaffmat;
  initaffmat=p_initaffmat;
  if ((initaffmat.Nrows()!=4) || (initaffmat.Ncols()!=4)) {
    initaffmat = IdentityMatrix(4);
  }
  J0 = initaffmat.SubMatrix(1,3,1,3);
  float alpha, detJ;
  int xoff[4], yoff[4], zoff[4];
  for (int z=grad.minz(); z<=grad.maxz()-1; z++) {
    for (int y=grad.miny(); y<=grad.maxy()-1; y++) {
      for (int x=grad.minx(); x<=grad.maxx()-1; x++) {
	for (int jacnum=0; jacnum<8; jacnum++) {
	  if ((jvol[jacnum](x,y,z)<minJ) || (jvol[jacnum](x,y,z)>maxJ)) {
	    get_jac_offset(jacnum,xoff,yoff,zoff);
	    // interpolate between current J matrix and initaffmat
	    for (int n1=1; n1<=3; n1++) { for (int n2=1; n2<=3; n2++) {
		J(n1,n2) = grad[3*n1+n2-4](x+xoff[n2],y+yoff[n2],z+zoff[n2]);
	    } }
	    alpha = 0.0;
	    Jnew = J;
	    detJ = Jnew.Determinant();
// 	      if (fabs(detJ - jvol[jacnum](x,y,z))>0.01) {
// 		cout << "ERROR: miscmatched jacobians " << detJ << " and "
// 		     << jvol[jacnum](x,y,z) << " for jacnum = "<<jacnum<<endl;
// 	      }
	    while ( (detJ > maxJ) || (detJ < minJ) ) {
	      alpha += 0.1;
	      if (alpha>1.0) alpha=1.0;
	      Jnew = (1 - alpha ) * J + alpha * J0;
	      detJ = Jnew.Determinant();
	    }
	    // and another one for luck!
	    alpha += 0.1;
	    if (alpha>1.0) alpha=1.0;
	    // rescale gradients as required
	    for (int n1=1; n1<=3; n1++) { for (int n2=1; n2<=3; n2++) {
		grad[3*n1+n2-4](x+xoff[n2],y+yoff[n2],z+zoff[n2]) = 
		(1 - alpha) * J(n1,n2) + alpha * J0(n1,n2);
	    } }
// 	      cout << "Rescaling with alpha = " << alpha << endl; 
// 	      if (alpha>0.05) {
// 		cout << "New J = " << Jnew << endl;
// 	      }
	  }
	}
      }
    }
  }
}


void limit_grad(volume4D<float>& grad, const volume4D<float>& jvol, 
		float minJ, float maxJ)
{
  Matrix initaffmat;
  initaffmat = IdentityMatrix(4);
  limit_grad(grad,jvol,minJ,maxJ,initaffmat);

}

void constrain_topology(volume4D<float>& warp, float minJ, float maxJ)
{
  ColumnVector jstats(4);
  jstats=jacobian_quick_check(warp,minJ,maxJ);
  volume4D<float> grad, jvol;
  int n=1, maxit=10;
  while ( (n++<maxit) && ( (jstats(3)>0.5) || (jstats(4)>0.5) ) ) {
    grad_calc(grad,warp);
    jacobian_check(jvol,jstats,warp,minJ,maxJ);
    // cout << "Jacobian stats of (min,max,#<min,#>max): "<<jstats.t()<<endl;
    limit_grad(grad,jvol,minJ,maxJ);
    integrate_gradient_field(warp,grad,warp[0].mean(),warp[1].mean(),warp[2].mean());
    jstats=jacobian_quick_check(warp,minJ,maxJ);
    // cout << "Jacobian quick stats of (min,max,#<min,#>max): "<<jstats.t()<<endl;
  }
}

void constrain_topology(volume4D<float>& warp)
{
  constrain_topology(warp,0.01,100.0);  // mainly just enforcing positivity
}

}

