/*  rungekutta.cc

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 2002 University of Oxford  */

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

#include "rungekutta.h"

using namespace std;

namespace MISCMATHS {

void rk(ColumnVector& ret, const ColumnVector& y, const ColumnVector& dy, float x, float h, const Derivative& deriv,const ColumnVector& paramvalues)
{ 
  Tracer tr("rk"); 

  float hh=h*0.5;
  float xh=x+hh;
  
  //first step
  ColumnVector yt=y+hh*dy;
  
  //second step
  ColumnVector dyt = deriv.evaluate(xh,yt,paramvalues);
  yt=y+hh*dyt;
  
  //third step
  ColumnVector dym = deriv.evaluate(xh,yt,paramvalues);
  yt=y+h*dym;
  dym=dym+dyt;
  
  //fourth step
  dyt = deriv.evaluate(x+h,yt,paramvalues);
  
  //addup
  ret = y+h*(dy+dyt+2*dym)/6;
}

void rkqc(ColumnVector& y, float& x, float& hnext, ColumnVector& dy, float htry, float eps, const Derivative& deriv,const ColumnVector& paramvalues)
{
  Tracer tr("rkqc"); 

  float xsav = x;
  ColumnVector dysav = dy;
  ColumnVector ysav = y;
  float h = htry;
  float hdid;
  ColumnVector ytemp;

  while(true)
    {
      // take 2 1/2 step sizes
  
      // first 1/2 step
      float hh=h*0.5;

      rk(ytemp,ysav,dysav,xsav,hh,deriv,paramvalues);
  
      // second 1/2 step
      x=xsav+hh;
      dy = deriv.evaluate(x,ytemp,paramvalues);
      rk(y,ytemp,dysav,xsav,hh,deriv,paramvalues);
      x=xsav+h;
      if(x==xsav) cerr << "step size too small" << endl;

      // take large step size
      rk(ytemp,ysav,dysav,xsav,h,deriv,paramvalues);
   
      // eval accuracy
      float errmax = 0.0;
      for(int i=1; i<=y.Nrows(); i++)
	{
	  //errmax=max(abs((y-ytemp)./y));
	  
	  float tmp = fabs((y(i)-ytemp(i))/y(i));
	  if(tmp > errmax) errmax = tmp;
	}

      errmax=errmax/eps;
      
      if(errmax <=1.0) 
	{
	  // step OK, compute step size for next step
	  hdid=h;
	  
	  if(errmax>6e-4)
	    hnext=h*std::exp(-0.2*std::log(errmax));
	  else
	    hnext=4*h;
	  
	  break;
      }
      else 
	{
	  // step too large,
	  h=h*std::exp(-0.25*std::log(errmax));
	}
    }

  y = y+(y-ytemp)/15;
}

void runge_kutta(Matrix& yp, ColumnVector& xp, ColumnVector& hp, const ColumnVector& ystart, float x1, float x2, float eps, float hmin, const Derivative& deriv,const ColumnVector& paramvalues)
{
  Tracer tr("runge_kutta"); 

  int MAXSTEP=1000;

  ColumnVector y = ystart;

  float x=x1;
  xp.ReSize(MAXSTEP,1);
  xp = 0;
  xp(1) =x1;

  float h=hp(1);
  hp.ReSize(MAXSTEP,1);
  hp = 0;

  yp.ReSize(MAXSTEP,y.Nrows());
  yp = 0;

  int kout=1;
  
  ColumnVector dy;

  for(int k=1; k <= MAXSTEP; k++)
    { 
      dy = deriv.evaluate(x,y,paramvalues);

      // store results:
      xp(kout)=x;
      yp.Row(kout)=y;
      hp(kout)=h;
  
      kout=kout+1;
    
      // stop overshoot of step past x2:
      if((x+h-x2)*(x+h-x1)>0) h=x2-x;

      float hnext = 0.0;
      rkqc(y,x,hnext,dy,h,eps,deriv,paramvalues);

      if((x-x2)*(x2-x1) >= 0.0)
      {
	xp(kout)=x;
	yp.Row(kout)=y;
	hp(kout)=h;
	//kout=kout+1;
	
        xp = xp.Rows(1,kout);
	yp = yp.Rows(1,kout);

	return;
      }
      else
      {
	if(hnext<=hmin) cerr << "step size too small" << endl;
	h=hnext;
      } 
      
    }
  cerr << "too many steps" << endl;
}

}
