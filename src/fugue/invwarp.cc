/*  invwarp.cc

    Mark Jenkinson and Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 2001-2006 University of Oxford  */

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

#include "utils/options.h"
#include "miscmaths/miscmaths.h"
#include "warpfns/warpfns.h"
#include "warpfns/fnirt_file_reader.h"

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

using namespace Utilities;
using namespace NEWMAT;
using namespace MISCMATHS;
using namespace NEWIMAGE;

volume4D<float> global_costim;
volume4D<float> global_warp;

bool abs_warp=true;

////////////////////////////////////////////////////////////////////////////

// COMMAND LINE OPTIONS

string title="invwarp \nCopyright(c) 2007, University of Oxford (Mark Jenkinson)";
string examples="invwarp -w warpvol -o invwarpvol -r refvol";

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<bool> debug(string("--debug"), false,
		  string("turn on debugging output"),
		  false, no_argument);
Option<bool> abswarp(string("--abs"), false,
		  string("use absolute warp convention (default): x' = w(x)"),
		  false, no_argument);
Option<bool> relwarp(string("--rel"), false,
		  string("use relative warp convention: x' = x + w(x)"),
		  false, no_argument);
Option<bool> nojaccon(string("--noconstraint"), false,
		  string("do not apply the Jacobian constraint"),
		  false, no_argument);
Option<int> niter(string("-n,--niter"), 10,
			string("number of gradient-descent iterations (default=10)"),
			false, requires_argument);
Option<float> lambdaval(string("--regularise"), 1.0,
			string("regularisation strength (default=1.0)"),
			false, requires_argument);
Option<float> jmin(string("--jmin"), 0.01,
			string("minimum acceptable Jacobian value for constraint (default 0.01)"),
			false, requires_argument);
Option<float> jmax(string("--jmax"), 100.0,
			string("maximum acceptable Jacobian value for constraint (default 100.0)"),
			false, requires_argument);
Option<string> initwarpname(string("--initwarp"), string(""),
			string("filename for initial warp transform (volume)"),
			false, requires_argument);
Option<string> refvolname(string("-r,--ref"), string(""),
		       string("filename for new reference image, i.e., what was originally the input image (determines inverse warpvol's FOV and pixdims)"),
		       true, requires_argument);
Option<string> outname(string("-o,--out"), string(""),
		       string("filename for output (inverse warped) image"),
		       true, requires_argument);
Option<string> warpname(string("-w,--warp"), string(""),
			string("filename for warp/shiftmap transform (volume)"),
			true, requires_argument);


////////////////////////////////////////////////////////////////////////////


// Optimisation functions

  // A temporary fix of including the std:: in front of all abs() etc
  //  has been done for now
  using std::abs;

  bool estquadmin(float &xnew, float x1, float xmid, float x2, 
		   float y1, float ymid, float y2)
  {
    // Finds the estimated quadratic minimum's position
    float ad=0.0, bd=0.0, det=0.0;
    ad = (xmid - x2)*(ymid - y1) - (xmid - x1)*(ymid - y2);
    bd = -(xmid*xmid - x2*x2)*(ymid - y1) + (xmid*xmid - x1*x1)*(ymid - y2);
    det = (xmid - x2)*(x2 -x1)*(x1 - xmid);
    if ((fabs(det)>1e-15) && (ad/det < 0)) {  // quadratic only has a maxima
      xnew = 0.0;
      return false;
    }
    if (fabs(ad)>1e-15) {
      xnew = -bd/(2*ad);
      return true;
    } else {  // near linear condition -> get closer to an end point
      xnew = 0.0;
      return false;
    }
    return false;
  }


  float extrapolatept(float x1, float xmid, float x2)
  {
    // xmid must be between x1 and x2
    // use the golden ratio (scale similar result)
    const float extensionratio = 0.3819660;
    float xnew;
    if (fabs(x2-xmid)>fabs(x1-xmid)) {
      xnew = extensionratio * x2 + (1 - extensionratio) * xmid;
    } else {
      xnew = extensionratio * x1 + (1 - extensionratio) * xmid;
    }
    return xnew;
  }
  


  float nextpt(float x1, float xmid, float x2, float y1, float ymid, float y2)
  {
    // x1 and x2 are the bounds, xmid is between them

    float xnew;
    bool quadok=false;
    quadok = estquadmin(xnew,x1,xmid,x2,y1,ymid,y2);

    // check to see that the quadratic result is in the range
    if ((!quadok) || (xnew < Min(x1,x2)) || (xnew > Max(x1,x2))) {
      xnew = extrapolatept(x1,xmid,x2);
    }
    return xnew;
  }

      

  void findinitialbound(float &x1, float &xmid, float &x2, 
			float &y1, float &ymid, float &y2, 
			float (*func)(const volume4D<float> &),
			const volume4D<float> &unitdir, 
			const volume4D<float> &pt)
  {
    const float extrapolationfactor = 1.6;
    const float maxextrap = extrapolationfactor*2;
    if (y1==0)  y1 = (*func)(x1*unitdir + pt);
    if (ymid==0)  ymid = (*func)(xmid*unitdir + pt);
    if (y1<ymid) {   // swap a and b if this is the case
      float tempx = x1, tempy = y1;
      x1 = xmid;     y1 = ymid;
      xmid = tempx;  ymid = tempy;
    }

    float newx2 = 0.0, newy2=0.0, maxx2=0.0;
    float dir=1.0;
    if (xmid<x1) dir=-1.0;

    bool quadok;

    x2 = xmid + extrapolationfactor*(xmid - x1);
    y2 = (*func)(x2*unitdir + pt);

    while (ymid > y2) {  // note: must maintain y1 >= ymid
	
      // cout << "    <" << Min(x1,x2) << "," << xmid 
      //   << "," << Max(x1,x2) << ">" << endl;
      maxx2 = xmid + maxextrap*(x2 - xmid);
      quadok = estquadmin(newx2,x1,xmid,x2,y1,ymid,y2);
      if ((!quadok) || ((newx2 - x1)*dir<0) || ((newx2 - maxx2)*dir>0)) {
	newx2 = xmid + extrapolationfactor*(x2-x1);
      }
      
      newy2 = (*func)(newx2*unitdir + pt);

      if ((newx2 - xmid)*(newx2 - x1)<0) {  // newx2 is between x1 and xmid
	if (newy2 < ymid) {  // found a bracket!
	  x2 = xmid;  y2 = ymid;
	  xmid = newx2;  ymid = newy2;
	  break;
	} else {  // can use newx2 as a new value for x1 (as newy2 >= ymid)
	  x1 = newx2;  y1 = newy2;
	}
      } else {  // newx2 is between xmid and maxx2
	if (newy2 > ymid) { // found a bracket!
	  x2 = newx2;  y2 = newy2;
	  break;
	} else if ((newx2 - x2)*dir<0) {  // newx2 closer to xmid than old x2
	  x1 = xmid;  y1 = ymid;
	  xmid = newx2;  ymid = newy2;
	} else {
	  x1 = xmid;  y1 = ymid;
	  xmid = x2;  ymid = y2;
	  x2 = newx2;  y2 = newy2;
	}
      }
	
    }

    if ( (y2<ymid) || (y1<ymid) ) {
      cerr << "findinitialbound failed to bracket: current triplet is" << endl;
    }
  }
  

  float optimise1d(volume4D<float> &pt, const volume4D<float>& unitdir, 
		  float unittol, int &iterations_done, 
		  float (*func)(const volume4D<float>&), int max_iter,
		  float init_value, float boundguess) 
  {
    // Golden Search Routine
    // Must pass in the direction vector in N-space (dir), the initial
    //  N-dim point (pt), the acceptable tolerance (tol) and other
    //  stuff
    // Note that the length of the direction vector is unimportant
    // Pass in previous costfn value as init_value, if known, otherwise
    //  pass in 0.0 and it will force the calculation
    // Unlike the version in optimise.cc the boundguess is in absolute
    //  units, not in units of unittol

    float y1,y2,ymid;
    float x1,x2,xmid;

    // set up initial points
    xmid = 0.0;
    x1 = boundguess;  // initial guess (bound)
    if (init_value==0.0) ymid = (*func)(xmid*unitdir + pt);
    else ymid = init_value;
    y1 = (*func)(x1*unitdir + pt);
    findinitialbound(x1,xmid,x2,y1,ymid,y2,func,unitdir,pt);

    if (verbose.value()) {
      cout << "BOUND = (" << x1 << "," << y1 << ")  ";
      cout << "(" << xmid << "," << ymid << ")  ";
      cout << "(" << x2 << "," << y2 << ")" << endl;
    }

    float min_dist = 0.1 * unittol;
    float xnew, ynew;
    int it=0;
    while ( ((++it)<=max_iter) && (fabs((x2-x1)/unittol)>1.0) )
      {
	// cout << "  [" << Min(x1,x2) << "," << Max(x1,x2) << "]" << endl;

	if (it>0) {
	  xnew = nextpt(x1,xmid,x2,y1,ymid,y2);
	} else {
	  xnew = extrapolatept(x1,xmid,x2);
	}

	float dirn=1.0;
	if (x2<x1) dirn=-1.0;

	if (fabs(xnew - x1)<min_dist) {
	  xnew = x1 + dirn*min_dist;
	}

	if (fabs(xnew - x2)<min_dist) {
	  xnew = x2 - dirn*min_dist;
	}

	if (fabs(xnew - xmid)<min_dist) {
	  xnew = extrapolatept(x1,xmid,x2);
	}

	if (fabs(xmid - x1)<0.4*unittol) {
	  xnew = xmid + dirn*0.5*unittol;
	}

	if (fabs(xmid - x2)<0.4*unittol) {
	  xnew = xmid - dirn*0.5*unittol;
	}

	if (verbose.value()) { cout << "xnew = " << xnew << endl; }
	ynew = (*func)(xnew*unitdir + pt);

	if ((xnew - xmid)*(x2 - xmid) > 0) {  // is xnew between x2 and xmid ?
	  // swap x1 and x2 so that xnew is between x1 and xmid
	  float xtemp = x1;  x1 = x2;  x2 = xtemp;
	  float ytemp = y1;  y1 = y2;  y2 = ytemp;
	}
	if (ynew < ymid) {
	  // new interval is [xmid,x1] with xnew as best point in the middle
	  x2 = xmid;  y2 = ymid;
	  xmid = xnew;  ymid = ynew;
	} else {
	  // new interval is  [x2,xnew] with xmid as best point still
	  x1 = xnew;  y1 = ynew;
	}
      }
    iterations_done = it;
    pt = xmid*unitdir + pt;
    return ymid;
  }


///////////////////////////////////////////////////////////////////////////

int fill_image(volume4D<float>& invol, const volume<float>& mask)
{
  invol = sparseinterpolate(invol,mask);
  return 0;
}


int fill_warp(volume4D<float>& warpvol, const volume<float>& mask)
{
  int retval=0;
  convertwarp_abs2rel(warpvol);
  retval = fill_image(warpvol,mask);
  convertwarp_rel2abs(warpvol);
  return retval;
}


///////////////////////////////////////////////////////////////////////////

volume<float> initialise_invwarp(volume4D<float>& invwarp, const volume4D<float>& warpvol)
{
  if (invwarp.nvoxels()<=1) {
    invwarp = warpvol*0.0f;
  }
  float maxdist = (float) (2.0);
  volume<float> pixdist;
  pixdist=invwarp[0]*0.0f + maxdist;

  for (int z=warpvol.minz(); z<=warpvol.maxz(); z++) {
    for (int y=warpvol.miny(); y<=warpvol.maxy(); y++) {
      for (int x=warpvol.minx(); x<=warpvol.maxx(); x++) {
	float xd, yd, zd;
	// convert to voxel coordinates (in invwarp)
	xd = warpvol(x,y,z,0)/invwarp.xdim();
	yd = warpvol(x,y,z,1)/invwarp.ydim();
	zd = warpvol(x,y,z,2)/invwarp.zdim();
	int xn, yn, zn;
	xn = MISCMATHS::round(xd);
	yn = MISCMATHS::round(yd);
	zn = MISCMATHS::round(zd);
	if (invwarp.in_bounds(xn,yn,zn)) {
	  float dist=fabs((xd-xn)*(yd-yn)*(zd-zn));
	  if (dist<pixdist(xn,yn,zn)) {
	    pixdist(xn,yn,zn)=dist;
	    invwarp(xn,yn,zn,0)=x*warpvol.xdim();
	    invwarp(xn,yn,zn,1)=y*warpvol.ydim();
	    invwarp(xn,yn,zn,2)=z*warpvol.zdim();
	  }
	}
      }
    }
  }

  // make pixdist the valid mask
  pixdist = 1.0f - binarise(pixdist,maxdist -0.5f);
  if (debug.value() && outname.set()) {
    save_volume(pixdist,fslbasename(outname.value())+"_initmask");
  }

  fill_warp(invwarp,pixdist);

  return pixdist;
}

float calc_cost(volume4D<float>& costim, 
		const volume4D<float>& warp, 
		const volume4D<float>& invwarp, float lambda=0.0)
{
  // C = \sum (x_d - x_u)^2 / N = \sum c_x
  float cost=0.0;
  costim=0.0;
  for (int z=warp.minz(); z<=warp.maxz(); z++) {
    for (int y=warp.miny(); y<=warp.maxy(); y++) {
      for (int x=warp.minx(); x<=warp.maxx(); x++) {
	// convert to voxel coordinates (in invwarp)
	float xu=warp(x,y,z,0)/invwarp.xdim();
	float yu=warp(x,y,z,1)/invwarp.ydim();
	float zu=warp(x,y,z,2)/invwarp.zdim();
	if ((xu>=0) && (yu>=0) && (zu>=0) && (xu<=invwarp.maxx()) && 
	    (yu<=invwarp.maxy()) && (zu<=invwarp.maxz())) {
	  float xw=invwarp[0].interpolate(xu,yu,zu);
	  float yw=invwarp[1].interpolate(xu,yu,zu);
	  float zw=invwarp[2].interpolate(xu,yu,zu);
	  costim(x,y,z,0) = xw - x*warp.xdim();
	  costim(x,y,z,1) = yw - y*warp.ydim();
	  costim(x,y,z,2) = zw - z*warp.zdim();
	}
      }
    }
  }
  cost = costim.sumsquares()/costim.nvoxels();

  float TotLap=0.0, Lap=0.0, d2wdx2=0.0;
  if (lambda>0.0) {
    // now add in square Laplacian regularisation
    for (int c=0; c<=2; c++) {
      for (int z=invwarp.minz(); z<=invwarp.maxz(); z++) {
	for (int y=invwarp.miny(); y<=invwarp.maxy(); y++) {
	  for (int x=invwarp.minx(); x<=invwarp.maxx(); x++) {
	    if ((x>invwarp.minx()) && (x<invwarp.maxx())) {
	      d2wdx2 = 2.0*invwarp(x,y,z,c)-invwarp(x-1,y,z,c)-invwarp(x+1,y,z,c);
	      Lap+= d2wdx2*d2wdx2;
	    }
	    if ((y>invwarp.miny()) && (y<invwarp.maxy())) {
	      d2wdx2 = 2.0*invwarp(x,y,z,c)-invwarp(x,y-1,z,c)-invwarp(x,y+1,z,c);
	      Lap+= d2wdx2*d2wdx2;
	    }
	    if ((z>invwarp.minz()) && (z<invwarp.maxz())) {
	      d2wdx2 = 2.0*invwarp(x,y,z,c)-invwarp(x,y,z-1,c)-invwarp(x,y,z+1,c);
	      Lap+= d2wdx2*d2wdx2;
	    }
	  }
	}
	TotLap += Lap;
	Lap=0.0;
      }
    }    
    TotLap+=Lap;  // should be unnecessary
    if (debug.value()) { 
      cout << "Lap cost = " << ( lambda / costim.nvoxels() ) * TotLap << endl; 
      cout << "TotLap = " << TotLap << endl;
    }
    cost += ( lambda / costim.nvoxels() ) * TotLap;
  }

  if (verbose.value()) cout << "Cost = " << cost << endl;
  return cost;
}


float calc_cost_stub(const volume4D<float>& invwarp)
{
  return calc_cost(global_costim, global_warp, invwarp, lambdaval.value());
}

volume4D<float> cost_deriv(const volume4D<float>& warp, 
			   const volume4D<float>& invwarp, 
			   float& cost, float lambda=0.0)
{
  // C = \sum (x_d - x_u)^2 / N = \sum c_x
  // dC/dx_p = \sum 2*(x_d - x_u)*dx_d/dx_p / N
  //         = \sum 2 * c_x * dx_d/dx_p 
  if (debug.value()) { cout << "Lambda = " << lambda << endl; }
  volume4D<float> costim(warp);
  cost = calc_cost(costim,warp,invwarp,lambda);
  volume4D<float> deriv(invwarp);
  deriv = 0.0f;
  for (int z=warp.minz(); z<=warp.maxz(); z++) {
    for (int y=warp.miny(); y<=warp.maxy(); y++) {
      for (int x=warp.minx(); x<=warp.maxx(); x++) {
	// convert to voxel coordinates (in invwarp)
	float xu=warp(x,y,z,0)/invwarp.xdim();
	float yu=warp(x,y,z,1)/invwarp.ydim();
	float zu=warp(x,y,z,2)/invwarp.zdim();
	int xu0 = MISCMATHS::round(floor(xu));
	int yu0 = MISCMATHS::round(floor(yu));
	int zu0 = MISCMATHS::round(floor(zu));
	if (deriv.in_bounds(xu0,yu0,zu0) && 
	    deriv.in_bounds(xu0+1,yu0+1,zu0+1)) {
	  for (int zd=0; zd<=1; zd++) {
	    for (int yd=0; yd<=1; yd++) {
	      for (int xd=0; xd<=1; xd++) {
		float k=(1-fabs(xu0+xd-xu))*(1-fabs(yu0+yd-yu))*(1-fabs(zu0+zd-zu));
		deriv.value(xu0+xd,yu0+yd,zu0+zd,0) += 2*costim(x,y,z,0)*k;
		deriv.value(xu0+xd,yu0+yd,zu0+zd,1) += 2*costim(x,y,z,1)*k;
		deriv.value(xu0+xd,yu0+yd,zu0+zd,2) += 2*costim(x,y,z,2)*k;
	      }
	    }
	  }
	}
      }
    }
  }

  if (lambda>0.0) {
    // L = \sum (2*x1 - x0 - x2)^2 / N
    // dL/dx_p = \sum f*(2*x1-x0-x2) / N  with f=4 for x_p=x1, f=-2 otherwise 
    volume4D<float> dLap(deriv);
    dLap*=0.0;
    float d2wdx2=0.0;
    // now add in derivative of square Laplacian regularisation
    for (int c=0; c<=2; c++) {
      for (int z=deriv.minz(); z<=deriv.maxz(); z++) {
	for (int y=deriv.miny(); y<=deriv.maxy(); y++) {
	  for (int x=deriv.minx(); x<=deriv.maxx(); x++) {
	    if ((x>invwarp.minx()) && (x<invwarp.maxx())) {
	      d2wdx2 = 2.0*invwarp(x,y,z,c)-invwarp(x-1,y,z,c)-invwarp(x+1,y,z,c);
	      dLap(x,y,z,c) += 4.0*d2wdx2;
	      dLap(x-1,y,z,c) += -2.0*d2wdx2;
	      dLap(x+1,y,z,c) += -2.0*d2wdx2;
	    }
	    if ((y>invwarp.miny()) && (y<invwarp.maxy())) {
	      d2wdx2 = 2.0*invwarp(x,y,z,c)-invwarp(x,y-1,z,c)-invwarp(x,y+1,z,c);
	      dLap(x,y,z,c) += 4.0*d2wdx2;
	      dLap(x,y-1,z,c) += -2.0*d2wdx2;
	      dLap(x,y+1,z,c) += -2.0*d2wdx2;
	    }
	    if ((z>invwarp.minz()) && (z<invwarp.maxz())) {
	      d2wdx2 = 2.0*invwarp(x,y,z,c)-invwarp(x,y,z-1,c)-invwarp(x,y,z+1,c);
	      dLap(x,y,z,c) += 4.0*d2wdx2;
	      dLap(x,y,z-1,c) += -2.0*d2wdx2;
	      dLap(x,y,z+1,c) += -2.0*d2wdx2;
	    }
	  }
	}
      }
    }    
    deriv += (lambda / costim.nvoxels()) * dLap;
  }

  return deriv;
}


int gradient_descent(volume4D<float>& invwarp, const volume4D<float>& warpvol, int n_iter,
		     volume4D<float>& deriv, volume4D<float>& costim)
{
  // gradient descent: du = a*dC/du ; dC = dC/du . du = a*|dC/du|^2
  //                   dC = -C = a*|dC/du|^2 =>  a = -C / |dC/du|^2
  int iterations=0, totiter=0;
  float cost;

  for (int iter=1; iter<=n_iter; iter++) {
    deriv = cost_deriv(warpvol,invwarp,cost,lambdaval.value());
    // normalise derivative and calculate a1 as a guess of the 
    //  amount needed along the derivative to cause a useful single voxel shift
    if (deriv.sumsquares()<1e-12) {
      // no derivative probably means the warp is trivial and so nothing to do
      return 0;
    }
    deriv /= sqrt(deriv.sumsquares()/deriv.nvoxels());
    float a1=1.0/Max(fabs(deriv.max()),fabs(deriv.min()));
    global_costim = costim;
    global_warp = warpvol;
    cost = optimise1d(invwarp, deriv, 0.1*a1, iterations, 
		      calc_cost_stub, 10, 0.0, 10.0*a1); 
    totiter += iterations;
    if (verbose.value()) { cout << "Iteration #"<<iter<<" : cost = " << cost << endl; }
  }
  return totiter;
}


///////////////////////////////////////////////////////////////////////////


int invwarp()
{
 
 // read in images
  volume4D<float>    invwarp;
  volume4D<float>    warpvol;
  AbsOrRelWarps      spec_wt=UnknownWarps;        // Specified warp convention
  Matrix             skrutt=IdentityMatrix(4);    // Not used

  read_volume4D(warpvol,warpname.value());

  if (abswarp.value()) spec_wt = AbsoluteWarps;
  else if (relwarp.value()) spec_wt = RelativeWarps;
  FnirtFileReader    fnirtfile(warpname.value(),spec_wt,verbose.set());
  warpvol = fnirtfile.FieldAsNewimageVolume4D(true);
  convertwarp_rel2abs(warpvol); 

  // read reference volume and set size of invwarp
  read_volume4D(invwarp,refvolname.value());
  while (invwarp.tsize()<3) { invwarp.addvolume(invwarp[0]); }
  while (invwarp.tsize()>3) { invwarp.deletevolume(invwarp.maxt()); }
  // inwarp.tsize() == 3 here
  invwarp=0.0f;
  invwarp.setDisplayMaximumMinimum(0,0);
  if (debug.value()) { print_volume_info(invwarp,"invwarp"); }

  // set up the initial warp field
  volume<float> mask;
  if (initwarpname.set()) {
    read_volume4D(invwarp,initwarpname.value());
    mask=invwarp[0]*0.0f + 1.0f;
  } else {
    mask=initialise_invwarp(invwarp,warpvol);
  }
  if (debug.value()) { 
    print_volume_info(invwarp,"invwarp");
    if (outname.set()) {   // save results
      if (fnirtfile.AbsOrRel()==RelativeWarps) { convertwarp_abs2rel(invwarp); }
      save_volume4D(invwarp,fslbasename(outname.value())+"_init");
      if (fnirtfile.AbsOrRel()==RelativeWarps) { convertwarp_rel2abs(invwarp); }
    }
  }

  // do work
  volume4D<float> deriv, costim(warpvol);
  if (verbose.value()) { cout << "Perform gradient descent" << endl; }
  gradient_descent(invwarp,warpvol,niter.value(),deriv,costim);
  if (verbose.value()) { cout << "After gradient descent" << endl; }

  if (debug.value()) { 
      if (fnirtfile.AbsOrRel()==RelativeWarps) { convertwarp_abs2rel(invwarp); }
      save_volume4D(invwarp,fslbasename(outname.value())+"_prefill");
      if (fnirtfile.AbsOrRel()==RelativeWarps) { convertwarp_rel2abs(invwarp); }
  }

  if (!nojaccon.value()) {
    if (verbose.value()) { cout << "Re-fill image using valid mask" << endl; }
    fill_warp(invwarp,mask);
    
    if (debug.value()) { 
      if (fnirtfile.AbsOrRel()==RelativeWarps) { convertwarp_abs2rel(invwarp); }
      save_volume4D(invwarp,fslbasename(outname.value())+"_postfill");
      if (fnirtfile.AbsOrRel()==RelativeWarps) { convertwarp_rel2abs(invwarp); }
    }
    if (verbose.value()) { cout << "Constrain Jacobian" << endl; }
    constrain_topology(invwarp,jmin.value(),jmax.value());
    
//     // apply another round of gradient descent
//     gradient_descent(invwarp,warpvol,1,deriv,costim);
    
//     // reconstrain Jacobian (hopefully unnecessary)
//     constrain_topology(invwarp,jmin.value(),jmax.value());
  }


  if (verbose.value()) { cout <<"Constrained Jacobian"<< endl; }

  if (fnirtfile.AbsOrRel()==RelativeWarps) { convertwarp_abs2rel(invwarp); }
  save_volume4D(invwarp,outname.value());
  if (debug.value()) {
    copybasicproperties(invwarp,deriv);
    copybasicproperties(invwarp,costim);
    save_volume4D(deriv,fslbasename(outname.value())+"_deriv");
    save_volume4D(costim,fslbasename(outname.value())+"_costim");
  }

  return(EXIT_SUCCESS);
}




int main(int argc, char *argv[])
{

  Tracer tr("main");

  OptionParser options(title, examples);

  try {
    options.add(warpname);
    options.add(outname);
    options.add(refvolname);
    options.add(relwarp);
    options.add(abswarp);
    options.add(niter);
    options.add(lambdaval);
    options.add(initwarpname);
    options.add(nojaccon);
    options.add(jmin);
    options.add(jmax);
    options.add(debug);
    options.add(verbose);
    options.add(help);
    
    int nparsed = options.parse_command_line(argc, argv);
    if (nparsed < argc) {
      for (; nparsed<argc; nparsed++) {
        cerr << "Unknown input: " << argv[nparsed] << endl;
      }
      exit(EXIT_FAILURE);
    }

    if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
      {
	options.usage();
	exit(EXIT_FAILURE);
      }
  } catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } 

  if (abswarp.value() && relwarp.value()) {
    cerr << "--abs and --rel flags cannot both be set" << endl;
    exit(EXIT_FAILURE);
  }

  volume4D<float>   warpvol;
  try {
    read_volume4D_hdr_only(warpvol,warpname.value());
  }
  catch(...) {
    cerr << "invwarp: Problem reading warp-file " << warpname.value() << endl;
    exit(EXIT_FAILURE);
  }
  if ((warpvol.intent_code()==FSL_CUBIC_SPLINE_COEFFICIENTS || 
       warpvol.intent_code()==FSL_QUADRATIC_SPLINE_COEFFICIENTS || 
       warpvol.intent_code()==FSL_DCT_COEFFICIENTS ||
       warpvol.intent_code()==FSL_FNIRT_DISPLACEMENT_FIELD) &&
      (abswarp.value() || relwarp.value())) {
    cout << "--abs and --rel flags ignored when reading fnirt coefficient files" << endl;
  }

  return invwarp();
}

