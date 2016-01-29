/*  testgam.cc

    Mark Woolrich, Tim Behrens - FMRIB Image Analysis Group

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

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#define WANT_STREAM
#define WANT_MATH
#include "newmatap.h"
#include "newmatio.h"
#include <string>
#include <math.h>
#include "utils/log.h"
#include "gam.h"
#include "miscmaths/miscprob.h"
#include "miscmaths/miscmaths.h"
#include "newimage/newimageall.h"


#ifndef NO_NAMESPACE
using namespace Utilities;
using namespace Gs;
using namespace NEWMAT;
using namespace MISCMATHS;
using namespace NEWIMAGE;
#endif


 

  
float csevl(const float x, const ColumnVector& cs, const int n) 
{
 
  float b0 = 0;
  float b1 = 0;
  float b2 = 0;
  const float twox=2*x;
  
  for(int i=1; i<=n; i++)
    {
      b2=b1;
      b1=b0;
      b0=twox*b1-b2+cs(n+1-i);
    }
  
  return 0.5*(b0-b2);
}

float digamma(const float x) 
{ 
   ColumnVector psics(23);
  ColumnVector apsics(16);
  int ntapsi=16,ntpsi=23;
   psics << -.038057080835217922E0<<
	    .49141539302938713E0<<
	    -.056815747821244730E0<<
	    .008357821225914313E0<<
	    -.001333232857994342E0<<
	    .000220313287069308E0<<
	    -.000037040238178456E0<<
	    .000006283793654854E0<<
	    -.000001071263908506E0<<
	    .000000183128394654E0<<
	    -.000000031353509361E0<<
	    .000000005372808776E0<<
	    -.000000000921168141E0<<
	    .000000000157981265E0<<
	    -.000000000027098646E0<<
	    .000000000004648722E0<<
	    -.000000000000797527E0<<
	    .000000000000136827E0<<
	    -.000000000000023475E0<<
	    .000000000000004027E0<<
	    -.000000000000000691E0<<
	    .000000000000000118E0<<
	    -.000000000000000020E0;
	  
	  apsics <<-.0204749044678185E0<<
	    -.0101801271534859E0<<
	    .0000559718725387E0<<
	    -.0000012917176570E0<<
	    .0000000572858606E0<<
	    -.0000000038213539E0<<
	    .0000000003397434E0<<
	    -.0000000000374838E0<<
	    .0000000000048990E0<<
	    -.0000000000007344E0<<
	    .0000000000001233E0<<
	    -.0000000000000228E0<<
	    .0000000000000045E0<<
	    -.0000000000000009E0<<
	    .0000000000000002E0<<
	    -.0000000000000000E0;


  float y = fabs(x);
  float psi;

  if(y<2.0)
    {
      // do we need to deal with the following case?
      // c psi(x) for -2. .lt. x .lt. 2.

      int n = int(floor(x));
      y = x - n;
      n = n - 1;
      psi = csevl(2*y-1, psics, ntpsi);
      if(n==-1)
	{
	  psi = psi - 1.0/x;
	}
    }
  else
    {
      const float aux = csevl(8/(Sqr(y))-1, apsics, ntapsi);
      psi = log(x) - 0.5/x + aux;
    }
    
  return psi;
}

float mgradpt(const float v, const ColumnVector& xsq, const int P) 
{
    

  int n = xsq.Nrows();

  float sumoveri = 0.0;

  for(int i=1; i<=n; i++)
    {
      sumoveri += (v+P)*(xsq(i)/Sqr(v))/(2*(1+xsq(i)/v)) - 0.5*log(1+xsq(i)/v);
    }
    
  return -(n*(0.5*digamma((v+P)/2.0) - 0.5*P/v - 0.5*digamma(v/2.0)) + sumoveri);	     

}

float mdofls(const ColumnVector& xsq, const float phi, const int P) 
{    
    

  ColumnVector xsqscaled = xsq/phi;

  float vbig=1000;
  float vmid=50;
  float vsm=0.1;
  float fbig=mgradpt(vbig,xsqscaled,P);
  float fmid=mgradpt(vmid,xsqscaled,P);
  float fsm=mgradpt(vsm,xsqscaled,P);

  while((vbig-vmid)>0.5 && (vmid-vsm)>0.5)
    {  
      if(sign(fsm)!=sign(fmid))
	{
	  vbig=vmid;
	  fbig=fmid;
	  vmid=(vmid+vsm)/2.0;
	  fmid=mgradpt(vmid,xsqscaled,P);
	}
      else if(sign(fmid)!=sign(fbig))
	{
	  vsm=vmid;
	  fsm=fmid;
	  vmid=(vbig+vmid)/2.0;
	  fmid=mgradpt(vmid,xsqscaled,P);
	}
      else
	{
	  vsm=100;
	  vmid=100;
	  vbig=100;
	}
    }

  return round(vmid);
}
void multitfit(const Matrix& x, ColumnVector& m, SymmetricMatrix& covar, float& v) 
{

  int nevs=x.Nrows();
  int n = x.Ncols();
  v = 10;
  OUT("FOO1");
  Matrix mx;
  remmean(x,mx,m,2);
  
  covar = cov(mx.t());
  float tmp = pow(covar.Determinant(),(1.0/nevs));    
  covar=covar/tmp;

  // xsq(i) = x(i)'*inv(c)*x(i)
  ColumnVector xsq(n);
  SymmetricMatrix invcovar = covar.i();
  OUT("Smelly Germans");  
  for(int i = 1; i <= n; i++)
    {
      xsq(i) = (mx.Column(i).t()*invcovar*mx.Column(i)).AsScalar();
    }

  OUT("BLAHDEBLAH");
  float phi = tmp;

  for(int i = 1; i <= 20; i++)
    {
      float newphi = 0.0;
      for(int j = 1; j <= n; j++)
	{
	  float tau = phi*(v+nevs)/(v*phi+xsq(j));
	  newphi += tau*xsq(j);
	}
      phi = newphi/float(n*nevs);

      // 	write_ascii_matrix(xsq,'xsq');
      // 	OUT(nevs);
      // 	OUT(phi);
      v = mdofls(xsq,phi,nevs);
      //	OUT(v);
    } 
  OUT("BUGGER");
  covar = covar*phi;
}


int main(int argc, char *argv[])
{
  try{
    ColumnVector m;
    SymmetricMatrix covar;
    float v;
    volume4D<float> xavw;
    read_volume4D(xavw,"x");
    OUT("BUGGER1");
    ColumnVector xtmp = xavw.voxelts(0,0,0);
    OUT("BUGGER2");
    RowVector x=xtmp.t();
    OUT("BUGGER2");
      
    //    cout << x << endl;

    OUT("here");
    multitfit(x,m,covar,v);
    OUT("here2");
    OUT(m);
    OUT(covar);
    OUT(v);
    

  }
  catch(Exception p_excp) 
    {
      cerr << p_excp.what() << endl;
    }
  return 0;
}
