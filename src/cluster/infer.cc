/*  infer.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2000-2004 University of Oxford  */

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
#include <cstdlib>
#include <cmath>

#include "infer.h"
#include "libprob.h"
#include "miscmaths/miscmaths.h"

#define POSIX_SOURCE 1

#if !defined(M_PI)
#define M_PI (4 * atan(1.0))
#endif

using namespace std;

Infer::Infer(float udLh, float ut, unsigned int uV) {
  // the following bounds are checked to ensure that the exponent
  //  does not underflow, which is assumed to occur for results
  //  of less than 1e-37  => abs(t)<13.0

  // assign to copies
  dLh = udLh;
  t = ut;
  V = uV;
  if (V<=0.0) V=1.0;
  // dimensionality
  D = 3.0;
  // if (zsize <= 1)  D = 2.0;  // to be done by calling program

  // NB: the (sqr(t) -1) is previous D=3 version (from where??)
  if (fabs(t)<13.0) {
    Em_ = V * pow(double(2*M_PI),double(-(D+1)/2)) * dLh * pow((MISCMATHS::Sqr(t) - 1), (D-1)/2) *
      exp(-MISCMATHS::Sqr(t)/2.0); 
  } else {
    Em_ = 0.0;  // underflowed exp()
  }

  if (fabs(t)<8.0) {
    B_ = pow((MISCMATHS::gamma(1.0+D/2.0)*Em_)/(V*(0.5 + 0.5*MISCMATHS::erf(-t/sqrt(2.0)))),(2.0/D));
  } else {
    // the large t approximation  (see appendix below)
    float a1 = V * dLh * pow(double(2*M_PI),double(-(D+1)/2));
    float a3 = pow((MISCMATHS::gamma(1+D/2.0)  / V ),(2.0/D));
    float tsq = t*t;
    float c = pow(2*M_PI,-1.0/2.0) * t / ( 1.0 - 1.0/tsq + 3.0/(tsq*tsq)) ;
    float Em_q = a1 * pow(double(tsq - 1.0),double(D-1)/2) * c;
    B_ = a3 * pow(double(Em_q),double(2.0/D));
  }
  

//      cout << "E{m} " << Em_ << endl;
//      cout << "Beta = " << B_ << endl;
}
  
//////////////////////////////////////////////////////////////////////////////

// Calculate and return log(p)

float Infer::operator() (unsigned int k) {
  // ideally returns the following:
  //    return 1 - exp(-Em_ * exp(-B_ * pow( k , 2.0 / D)));
  // but in practice must be careful about ranges
  // Assumes that exp(+/-87) => 1e+/38 is OK for floats
  float exponent_thresh = 80.0;
  float arg1 = -B_ * pow(k , 2.0 / D);
  if (fabs(arg1)>exponent_thresh) {
    // approximation for logp
    float logp = arg1 + log(Em_);
    return logp;
  } else {
    float exp1 = exp(arg1);
    float arg2 = -Em_ * exp1;
    if (fabs(arg2)>exponent_thresh) {
      // approximation of  1 - exp(arg2)
      float p = -arg2; 
      if (p>0) return log(p);
    } else {
      float exp2 = exp(arg2);
      float p = 1.0 - exp2;
      if ( (p==0.0) && (arg2<0.0) ) { p = -arg2; } // approx for 1-exp2
      return log(p);
    }
  }
  cerr << "Warning: could not compute p-value accurately." << endl;
  return -500;
}



// MATHEMATICAL APPENDIX

/*

The formulas that need to be calculated are:
(1) E_m = V * dLh * (2*pi)^(-(D+1)/2) * (t^2 -1)^((D-1)/2) * exp(-t^2 /2)
(2) Beta = (Gamma(D/2+1)/V * E_m / Phi(-t) )^(2/D)
(3) p = 1 - exp( - E_m * exp(-Beta*k^(2/D)))

where Phi(-t) = Gaussian cumulant = (1/2 + 1/2*MISCMATHS::erf(-t/sqrt(2)))

These are approximated by:

(2a) Beta = (Gamma(D/2+1)/V)^(2/D) * (Em1)^(2/D) * Ct^(2/D)
where Em1 = V * dLh * (2*pi)^(-(D+1)/2) * (t^2 -1)^((D-1)/2)
      Ct = (2*pi)^(-1/2) * t / ( 1.0 - 1.0/t^2 + 3.0/t^4 )
      which approximates ( exp(-t^2 /2) / Phi(-t) )^(2/D)
      using 1/2 - 1/2*MISCMATHS::erf(t/sqrt(2)) = (2*pi)^(1/2) * exp(-t^2 /2) * 
                                        (1-1/t^2+3/t^4) / t
					(this is derived in TR00MJ1)

and

(3a) log(p) = (- Beta * k^(2/D)) + log(Em)

 */

