/*  f2z.cc

    Mark Woolrich & Mark Jenkinson, FMRIB Image Analysis Group

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

#include <cmath>
#include "f2z.h"
#include "utils/log.h"
#include "utils/tracer_plus.h"
#include <stdexcept>
#include "libprob.h"

using namespace NEWMAT;
using namespace Utilities;

namespace MISCMATHS {

  F2z* F2z::f2z = NULL;
 
  float F2z::largef2logp(float f, int d1, int d2)
  {
    Tracer_Plus ts("F2z::largef2logp");

    // no of iterations:
    int N = 20;

//     cout << f<< endl;
//     cout << d1<< endl;
//     cout << d2<< endl;

    if (f<=0.0) {
      cerr << "f cannot be zero or negative!" << endl;
      return 0.0;
    }
  
    if (d1<=0 || d2<=0) { 
      cerr << "DOFs cannot be zero or negative!" << endl;
      return 0.0; 
    }
 
    double alpha=d1/(double)d2;
    double m=(d1+d2)/2.0;
    double n=(1-d1/2.0);
    double loggam = (d1/2.0)*(::log(d1/(double)d2)-logbeta(d2/2.0,d1/2.0));

    //iter=f^(-n)/(alpha*(n+m-1)) + n*f^(-(n+1))/(alpha^2*(n+m-1)*(n+m)) + n*(n+1)*f^(-(n+2))/(alpha^3*(n+m-1)*(n+m)*(n+m+1));

    double top = 1.0;
    double bot = n+m-1;
    double iter = 0.0;
 
//     cerr << "logbeta(d2/2.0,d1/2.0)=" << logbeta(d2/2.0,d1/2.0) << endl;
//     cerr << "loggam = " << loggam << endl;
//     cerr << "n = " << n << endl;
//     cerr << "m = " << m << endl;

    for(int i = 1; i <= N; i++)
      {
	// cerr << "i=" << i;
		  iter = iter + top* ( std::pow( f,float(-(n+i-1)) ) / ( std::pow(alpha,double(i))*bot ) );	
	top = top*(n-1+i)*(-1);
	bot = bot*(n+m-1+i);
// 	cerr << "iter=" << iter;
      }


    if(iter <= 0) throw Exception("iter negative");

    float logp = loggam-(m-1)*(::log(1+alpha*f))+::log(iter);

//     cerr << "iter = " << iter << endl;
//     cerr << "logp = " << logp << endl;

    return logp;
  }

  bool F2z::islargef(float f, int d1, int d2, float &logp) {
   
    if(f > 2.0 && d1>1)
      {

	try
	  {
	    logp=largef2logp(f,d1,d2);	
	  }
	catch(Exception& p_excp) 
	  {
	    cerr << "Negative iter in F2z::largef2logp" << endl;
	    return false;
	  }

	return issmalllogp(logp);
      }
    else
      return false;
  }

  bool F2z::issmalllogp(float logp) {
    return (logp < -14.5);
  }

  float F2z::convert(float f, int d1, int d2) 
  {
    Tracer_Plus ts("F2z::convert");

    float z = 0.0, logp=0.0;

    if(!islargef(f,d1,d2,logp)) {

      double p = MISCMATHS::fdtr(d1, d2, f);

      z = MISCMATHS::ndtri(p);
    }
      else {

	z = logp2largez(logp);
      }

      return z;
    }

  void F2z::ComputeFStats(const ColumnVector& p_fs, int p_dof1, int p_dof2, ColumnVector& p_zs)
  {
    ColumnVector dof2 = p_fs;
    dof2 = p_dof2;
    ComputeFStats(p_fs,p_dof1,dof2,p_zs);
  }
  
  void F2z::ComputeFStats(const ColumnVector& p_fs, int p_dof1, const ColumnVector& p_dof2, ColumnVector& p_zs)
  {
    Tracer_Plus ts("F2z::ComputeFStats");
    
    int numTS = p_fs.Nrows();

    p_zs.ReSize(numTS);
    F2z& f2z = F2z::getInstance();
    
    for(int i = 1; i <= numTS; i++)
      {		  	
	if (p_fs(i) > 0.0)
	  {

// 	    cerr << "i=" << i;
// 	    cerr << ",p_fs(i)=" << p_fs(i);
// 	    cerr << ",p_dof1=" << p_dof1;
// 	    cerr << ",p_dof2=" << p_dof2(i) << endl;

	    p_zs(i) = f2z.convert(p_fs(i),int(p_dof1),int(p_dof2(i)));  
	  }
	else
	  {
	    p_zs(i) = 0.0;
	  }     
      }
  }
   void F2z::ComputeFStats(const ColumnVector& p_fs, const ColumnVector& p_dof1, const ColumnVector& p_dof2, ColumnVector& p_zs)
  {
    Tracer_Plus ts("F2z::ComputeFStats");
    
    int numTS = p_fs.Nrows();

    p_zs.ReSize(numTS);
    F2z& f2z = F2z::getInstance();
    
    for(int i = 1; i <= numTS; i++)
      {		  	
	if (p_fs(i) > 0.0)
	  {

// 	    cerr << "i=" << i;
// 	    cerr << ",p_fs(i)=" << p_fs(i);
// 	    cerr << ",p_dof1=" << p_dof1;
// 	    cerr << ",p_dof2=" << p_dof2(i) << endl;

	    p_zs(i) = f2z.convert(p_fs(i),int(p_dof1(i)),int(p_dof2(i)));  
	  }
	else
	  {
	    p_zs(i) = 0.0;
	  }     
      }
  }
  
}






























