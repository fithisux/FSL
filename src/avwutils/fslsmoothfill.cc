/*  smoothfill.cc

    Mark Jenkinson, FMRIB Image Analysis Group

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
#include "newimage/newimageall.h"

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

using namespace Utilities;
using namespace NEWMAT;
using namespace MISCMATHS;
using namespace NEWIMAGE;

volume<float> global_mask;

////////////////////////////////////////////////////////////////////////////

// COMMAND LINE OPTIONS

string title="smoothfill (Version 1.0)\nCopyright(c) 2012, University of Oxford (Mark Jenkinson)";
string examples="smoothfill -i inputimage -m mask -o smoothedresult";

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<bool> debug(string("--debug"), false,
		  string("turn on debugging output"),
		  false, no_argument);
Option<int> n_iter(string("-n,--niter"), 10,
		       string("number of iterations"),
		       false, requires_argument);
Option<string> maskname(string("-m,--mask"), string(""),
		       string("filename for mask image"),
		       true, requires_argument);
Option<string> outname(string("-o,--out"), string(""),
		       string("filename for output (inverse warped) image"),
		       true, requires_argument);
Option<string> inputname(string("-i,--in"), string(""),
			string("filename for input image"),
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
			float (*func)(const volume<float> &),
			const volume<float> &unitdir, 
			const volume<float> &pt)
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
  

  float optimise1d(volume<float> &pt, const volume<float>& unitdir, 
		  float unittol, int &iterations_done, 
		  float (*func)(const volume<float>&), int max_iter,
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


struct vec3p { int x; int y; int z; int n; };


double dilateval(const volume<float>& im, const volume<float>& mask, int x, int y, int z) 
{
  double sum=0;
  int n=0;
  for (int zz=Max(0,z-1); zz<=Min(z+1,im.maxz()); zz++) {
    for (int yy=Max(0,y-1); yy<=Min(y+1,im.maxy()); yy++) {
      for (int xx=Max(0,x-1); xx<=Min(x+1,im.maxx()); xx++) {
	if ((mask(xx,yy,zz)>(float) 0.5) && !((xx==x) && (yy==y) && (zz==z))) {
	  sum += im(xx,yy,zz);
	  n++;
	}
      }
    }
  }
  return sum/(Max(n,1));
}


bool ext_edge(const volume<float>& mask, int x, int y, int z) 
{
  if (mask(x,y,z)>(float) 0.5) { return false; }
  for (int zz=Max(0,z-1); zz<=Min(z+1,mask.maxz()); zz++) {
    for (int yy=Max(0,y-1); yy<=Min(y+1,mask.maxy()); yy++) {
      for (int xx=Max(0,x-1); xx<=Min(x+1,mask.maxx()); xx++) {
	if ((mask(xx,yy,zz)>(float) 0.5) && !((xx==x) && (yy==y) && (zz==z))) {
	  return true;
	}
      }
    }
  }
  return false;
}

int dilall_extra(volume<float>& im, volume<float>& mask) 
{
  if (!samesize(im,mask)) { cerr << "ERROR::dilall::image are not the same size" << endl; return 1; }
  deque<vec3p> ptlist;
  vec3p v, newv;
  // initial pass
  for (int z=0; z<=im.maxz(); z++) {
    for (int y=0; y<=im.maxy(); y++) {
      for (int x=0; x<=im.maxx(); x++) {
	if (ext_edge(mask,x,y,z)) {
	  v.x=x; v.y=y; v.z=z; v.n=2;
	  ptlist.push_front(v);
	}
      }
    }
  }
  while (!ptlist.empty()) {
    v = ptlist.back();
    ptlist.pop_back();
    if (mask(v.x,v.y,v.z)<=(float) 0.5) {  // only do it if the voxel is still unset
      im(v.x,v.y,v.z)=dilateval(im,mask,v.x,v.y,v.z);
      mask(v.x,v.y,v.z)=(float)v.n;
      // check neighbours and add them to the list if necessary
      for (int zz=Max(0,v.z-1); zz<=Min(v.z+1,im.maxz()); zz++) {
	for (int yy=Max(0,v.y-1); yy<=Min(v.y+1,im.maxy()); yy++) {
	  for (int xx=Max(0,v.x-1); xx<=Min(v.x+1,im.maxx()); xx++) {
	    if (ext_edge(mask,xx,yy,zz)) { newv.x=xx; newv.y=yy; newv.z=zz; newv.n=v.n+1; ptlist.push_front(newv); }
	  }
	}
      }
    }
  }
  return 0;
} 




///////////////////////////////////////////////////////////////////////////

float calc_cost(const volume<float>& invol)
{
  // Cost = Laplacian of image (only counting areas _outside_ of the mask - the portion to be filled)
  float Lap=0.0, d2wdx2=0.0;
  // now add in square Laplacian regularisation
  for (int z=invol.minz(); z<=invol.maxz(); z++) {
    for (int y=invol.miny(); y<=invol.maxy(); y++) {
      for (int x=invol.minx(); x<=invol.maxx(); x++) {
	if (global_mask(x,y,z)<0.5) {
	  if ((x>invol.minx()) && (x<invol.maxx())) {
	    d2wdx2 = 2.0*invol(x,y,z)-invol(x-1,y,z)-invol(x+1,y,z);
	    Lap+= d2wdx2*d2wdx2;
	  }
	  if ((y>invol.miny()) && (y<invol.maxy())) {
	    d2wdx2 = 2.0*invol(x,y,z)-invol(x,y-1,z)-invol(x,y+1,z);
	    Lap+= d2wdx2*d2wdx2;
	  }
	  if ((z>invol.minz()) && (z<invol.maxz())) {
	    d2wdx2 = 2.0*invol(x,y,z)-invol(x,y,z-1)-invol(x,y,z+1);
	    Lap+= d2wdx2*d2wdx2;
	  }
	}
      }
    }
  }

  Lap/=invol.nvoxels();

  if (verbose.value()) cout << "Cost = " << Lap << endl;
  return Lap;
}


volume<float> cost_deriv(const volume<float>& invol)
{
  // L = \sum (2*x1 - x0 - x2)^2 / N
  // dL/dx_p = \sum f*(2*x1-x0-x2) / N  with f=4 for x_p=x1, f=-2 otherwise 
  volume<float> dLap(invol);
  dLap*=0.0;
  float d2wdx2=0.0;
  bool gm0=false;
  // now add in derivative of square Laplacian regularisation
  for (int z=dLap.minz(); z<=dLap.maxz(); z++) {
    for (int y=dLap.miny(); y<=dLap.maxy(); y++) {
      for (int x=dLap.minx(); x<=dLap.maxx(); x++) {
	gm0=global_mask(x,y,z)<0.5;
	if ((x>invol.minx()) && (x<invol.maxx())) {
	  d2wdx2 = 2.0*invol(x,y,z)-invol(x-1,y,z)-invol(x+1,y,z);
	  if (gm0) { dLap(x,y,z) += 4.0*d2wdx2; }
	  if (global_mask(x-1,y,z)<0.5) { dLap(x-1,y,z) += -2.0*d2wdx2; }
	  if (global_mask(x+1,y,z)<0.5) { dLap(x+1,y,z) += -2.0*d2wdx2; }
	}
	if ((y>invol.miny()) && (y<invol.maxy())) {
	  d2wdx2 = 2.0*invol(x,y,z)-invol(x,y-1,z)-invol(x,y+1,z);
	  if (gm0) { dLap(x,y,z) += 4.0*d2wdx2; }
	  if (global_mask(x,y-1,z)<0.5) { dLap(x,y-1,z) += -2.0*d2wdx2; }
	  if (global_mask(x,y+1,z)<0.5) { dLap(x,y+1,z) += -2.0*d2wdx2; }
	}
	if ((z>invol.minz()) && (z<invol.maxz())) {
	  d2wdx2 = 2.0*invol(x,y,z)-invol(x,y,z-1)-invol(x,y,z+1);
	  if (gm0) { dLap(x,y,z) += 4.0*d2wdx2; }
	  if (global_mask(x,y,z-1)<0.5) { dLap(x,y,z-1) += -2.0*d2wdx2; }
	  if (global_mask(x,y,z+1)<0.5) { dLap(x,y,z+1) += -2.0*d2wdx2; }
	}
      }
    }
  }
  dLap /= invol.nvoxels();

  return dLap;
}


int gradient_descent(volume<float>& invol, volume<float>& deriv)
{
  // gradient descent: du = a*dC/du ; dC = dC/du . du = a*|dC/du|^2
  //                   dC = -C = a*|dC/du|^2 =>  a = -C / |dC/du|^2
  int iterations=0, totiter=0;
  float cost=0;

  for (int iter=1; iter<=n_iter.value(); iter++) {
    deriv = cost_deriv(invol);
    // normalise derivative and calculate a1 as a guess of the 
    //  amount needed along the derivative to cause a useful single voxel shift
    if (deriv.sumsquares()<1e-12) {
      // no derivative probably means the input is trivial and so nothing to do
      return 0;
    }
    // // either call the optimise1d function or just do some simple scaled gradient updates
    // deriv /= sqrt(deriv.sumsquares()/deriv.nvoxels());
    // float a1=1.0/Max(fabs(deriv.max()),fabs(deriv.min()));
    // cost = optimise1d(invol, deriv, 0.1*a1, iterations, 
    // 		      calc_cost, 10, 0.0, 10.0*a1); 
    if (iter==1) cost=calc_cost(invol);
    invol += ((float) (-0.1*cost/deriv.sumsquares()))*deriv;
    cost = calc_cost(invol);
    totiter += iterations;
    if (verbose.value()) { cout << "Iteration #"<<iter<<" : cost = " << cost << endl; }
  }
  return totiter;
}


///////////////////////////////////////////////////////////////////////////


int do_work()
{
  // read in images
  volume<float> invol;
  read_volume(invol,inputname.value());
  
  // read reference volume and set size of invol
  read_volume(global_mask,maskname.value());

  // do work
  volume<float> deriv;
  if (verbose.value()) { cout << "Perform gradient descent" << endl; }
  //gradient_descent(invol,deriv);
  // dilate the given values, but this has "streaks", so blur afterwards
  volume<float> mask;
  mask=global_mask;
  dilall_extra(invol,mask);
  save_volume(invol,fslbasename(outname.value())+"_init");
  float cost = calc_cost(invol);
  cout << "Cost is " << cost << endl;
  // generate subsampled versions of the dilated volume for spatial blurring (distance dependent)
  volume<float> vol2, vol4, vol8, vol16, vol32;
  invol.setextrapolationmethod(extraslice);
  vol2 = subsample_by_2(invol,true);
  vol4 = subsample_by_2(vol2,true);
  vol8 = subsample_by_2(vol4,true);
  vol16 = subsample_by_2(vol8,true);
  vol32 = subsample_by_2(vol16,true);
  save_volume(vol2,fslbasename(outname.value())+"_vol2");
  save_volume(vol32,fslbasename(outname.value())+"_vol32");
  for (int z=0; z<=invol.maxz(); z++) {
    for (int y=0; y<=invol.maxy(); y++) {
      for (int x=0; x<=invol.maxx(); x++) {
	// interpolate between subsampled versions of the input image (gives smooth, but spatially non-uniform, blurring)
	float idx=log(mask(x,y,z))/log(2.0);   // values in mask encode "distance"
	idx-=1;  // so that it only starts interpolating at higher values
	if (idx>4) idx=4;
	if (idx>0) {
	  int idx1=(int) floor(idx);
	  int idx2=(int) ceil(idx);
	  float val1=0, val2=0; 
	  if (idx1==0) { val1=vol2.interpolate(x/2.0,y/2.0,z/2.0); val2=vol4.interpolate(x/4.0,y/4.0,z/4.0); }
	  if (idx2==0) val2=val1;  // in case idx==0
	  if (idx1==1) { val1=vol4.interpolate(x/4.0,y/4.0,z/4.0); val2=vol8.interpolate(x/8.0,y/8.0,z/8.0); }
	  if (idx1==2) { val1=vol8.interpolate(x/8.0,y/8.0,z/8.0); val2=vol16.interpolate(x/16.0,y/16.0,z/16.0); }
	  if (idx1==3) { val1=vol16.interpolate(x/16.0,y/16.0,z/16.0); val2=vol32.interpolate(x/32.0,y/32.0,z/32.0); }
	  if (idx1>=4) { val1=vol32.interpolate(x/32.0,y/32.0,z/32.0); val2=val1; }
	  invol(x,y,z) = val1 + (idx-idx1)*(val2-val1);
	}
      }
    }
  }
  cost = calc_cost(invol);
  cout << "Cost is " << cost << endl;
  if (verbose.value()) { cout << "After gradient descent" << endl; }
  
  save_volume(invol,outname.value());
  save_volume(mask,fslbasename(outname.value())+"_idxmask");
  
  return(EXIT_SUCCESS);
}




int main(int argc, char *argv[])
{

  Tracer tr("main");

  OptionParser options(title, examples);

  try {
    options.add(inputname);
    options.add(outname);
    options.add(maskname);
    options.add(n_iter);
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

  volume<float>   invol;
  try {
    read_volume_hdr_only(invol,inputname.value());
  }
  catch(...) {
    cerr << "smoothfill: Problem reading input image " << inputname.value() << endl;
    exit(EXIT_FAILURE);
  }

  return do_work();
}

