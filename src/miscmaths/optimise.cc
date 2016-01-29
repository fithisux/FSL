/*  optimise.cc

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

// Mathematical optimisation functions


#include <iostream>
#include <cmath>
#include "optimise.h"
#include "miscmaths.h"

namespace MISCMATHS {

  // The following lines are ignored by the current SGI compiler
  //  (version egcs-2.91.57)
  // A temporary fix of including the std:: in front of all abs() etc
  //  has been done for now
  using std::abs;
  using std::sqrt;
  using std::exp;
  using std::log;

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
			float (*func)(const ColumnVector &),
			const ColumnVector &unitdir, const ColumnVector &pt)
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
  

  float optimise1d(ColumnVector &pt, const ColumnVector dir, 
		  const ColumnVector tol, int &iterations_done, 
		  float (*func)(const ColumnVector &), int max_iter,
		  float &init_value, float boundguess) 
  {
    // Golden Search Routine
    // Must pass in the direction vector in N-space (dir), the initial
    //  N-dim point (pt), the acceptable tolerance (tol) and other
    //  stuff
    // Note that the length of the direction vector is unimportant

    float y1,y2,ymid;
    float x1,x2,xmid;

    // Calculate dot product of dir by tol
    //  st (x1-x2)*dir_tol = average number of tolerances between x1 and x2
    float dir_tol = 0.0;
    ColumnVector unitdir;
    unitdir = dir/std::sqrt(dir.SumSquare());
    for (int n=1; n<=tol.Nrows(); n++) {
      if (fabs(tol(n))>1e-15) {
	dir_tol += fabs(unitdir(n)/tol(n));
      }
    }
    float unittol = fabs(1/dir_tol);

    // set up initial points
    xmid = 0.0;
    x1 = boundguess * unittol;  // initial guess (bound)
    if (init_value==0.0) { init_value = (*func)(xmid*unitdir + pt); }
    ymid = init_value;
    y1 = (*func)(x1*unitdir + pt);
    findinitialbound(x1,xmid,x2,y1,ymid,y2,func,unitdir,pt);

    // cout << "(" << x1 << "," << y1 << ")  ";
    // cout << "(" << xmid << "," << ymid << ")  ";
    // cout << "(" << x2 << "," << y2 << ")" << endl;

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



  
  float optimise(ColumnVector &pt, int numopt, const ColumnVector &tol, 
		 float (*func)(const ColumnVector &), int &iterations_done, 
		 int max_iter, const ColumnVector& boundguess, 
		 const string type)
  {
    // Note that numopt can be less than pt.Nrows() - e.g. 6 dof optimisation
    //  but with a 12 dimensional vector

    // Calculate dot product of dir by tol
    //  st (x1-x2)*dir_tol = average number of tolerances between x1 and x2
    ColumnVector inv_tol(tol.Nrows());
    inv_tol = 0.0;
    for (int n=1; n<=tol.Nrows(); n++) {
      if (fabs(tol(n))>1e-15) {
	inv_tol(n) = fabs(1.0/tol(n));
      }
    }
    inv_tol /= (float) tol.Nrows();

    Matrix dirs(pt.Nrows(),pt.Nrows());
    dirs = IdentityMatrix(pt.Nrows());
    ColumnVector dir(pt.Nrows()), initpt, deltaf(pt.Nrows());
    deltaf=0.0f;
    int lit=0, littot=0, it=0;
    float fval=0.0, fval2=0.0, bndguess, finit=0.0, fend=0.0, fextrap=0.0;
    while ((++it)<=max_iter)
      {
	initpt = pt;
	bndguess = boundguess(Min(it,boundguess.Nrows()));  // ceiling of nrows
	for (int n=1; n<=numopt; n++) {
	  for (int m=1; m<=pt.Nrows(); m++) { dir(m) = dirs(m,n); }
	  fval2 = optimise1d(pt,dir,tol,lit,func,100,fval,bndguess);
	  deltaf(n)=fval2-fval;
	  if (n==1) { finit = fval; }
	  fval=fval2;
	  littot += lit;
	}

	// check to see if the point has moved more than one average tolerance
	float avtol = SP((initpt - pt),inv_tol).SumAbsoluteValue();
	if (avtol < 1.0) break;

	// if continuing then change the directions if using Powell's method
	if (type=="powell") 
	  {
	    // find direction of maximal change
	    int bestm=1;
	    for (int m=1; m<=numopt; m++) { 
	      if (deltaf(m)<deltaf(bestm)) bestm=m; 
	    }
	    fend=fval;
	    fextrap=(*func)(initpt + 2*(pt-initpt));
	    float df=fabs(deltaf(bestm));
	    if ( (2 * (finit-2*fend+fextrap) * (finit-fend-df)*(finit-fend-df)) < ( (finit-fextrap)*(finit-fextrap)*df ) ) {
	      if (fextrap<finit) {
		cout << "Applying POWELL correction" << endl;
		cout << "finit, fend, fextrap = " << finit << " , " << fend << " , " << fextrap << endl;
		// do another minimisation
		fval2 = optimise1d(pt,pt-initpt,tol,lit,func,100,fval,bndguess);
		fval=fval2;
		cout << "fval = " << fval << endl;
		littot += lit;
		// replace direction of maximum change with pt-initpt
		for (int m=1; m<=pt.Nrows(); m++) {
		  dirs(m,bestm)=pt(m)-initpt(m);
		}
	      }
	    }
	  }
	
      }
    // cout << endl << "Major iterations = " << it << endl;
    iterations_done = littot;
    return fval;
  }


}








