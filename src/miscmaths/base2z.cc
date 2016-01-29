/*  base2z.cc

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
#include "base2z.h"
#include "libprob.h"

namespace MISCMATHS {

  Base2z* Base2z::base2z = NULL;

  float Base2z::logbeta(float v, float w)
  {
    return MISCMATHS::lgam(v)+MISCMATHS::lgam(w)-MISCMATHS::lgam(v+w);
  }

  float Base2z::logp2largez(float logp) 
  {
    // Large Z extrapolation routine for converting log(p) to Z values
    //  written by Mark Jenkinson, March 2000
    //
    // Equations were derived by using integration by parts and give the
    //  following formulae:
    //   Z to log(p)
    //       log(p) = -1/2*z*z - 1/2*log(2*pi) - log(z) 
    //                + log(1 - 1/(z*z) + 3/(z*z*z*z))
    // this equation is then solved by the recursion:
    //   z_0 = sqrt(2*(-log(p) - 1/2*log(2*pi)))
    //   z_{n+1} = sqrt(2*(-log(p) - 1/2*log(2*pi) - log(z_n)  
    //             + log(1 - 1/(zn*zn) + 3/(zn*zn*zn*zn)) ))
    // In practice this recursion is quite accurate in 3 to 5 iterations
    // The equation is accurate to 1 part in 10^3 for Z>3.12  (3 iterations)


    static const float pi = 3.141592653590;
    static const float log2pi = log(2*pi);
    
    float z0, zn;
    // iteratively solve for z given log p
    float b = -2*logp - log2pi; 
    z0 = sqrt(b);
    zn = z0;
    for (int m=1; m<=3; m++) {
      // zn = sqrt(b + 2*log(1/zn - 1/(zn*zn*zn) + 3/(zn*zn*zn*zn*zn)));
      zn = sqrt(b + 2*log(((3/(zn*zn) - 1)/(zn*zn) + 1)/zn) );
    }
  
    return zn;
  }
  
  float Base2z::convertlogp2z(float logp) 
    {
      // logp must be the *natural* logarithm of p, not base 10
      float z = 0.0;
      
      if(!issmalllogp(logp)) {
	z = MISCMATHS::ndtri(exp(logp));
      }
      else {
	z = logp2largez(logp);
      }

      return z;
    }

}






























