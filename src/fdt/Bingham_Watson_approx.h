/*  Directional Statistics Functions

    Bingham and Watson Distributions and functions to approximate their normalizing constant
    
    Stam Sotiropoulos  - FMRIB Image Analysis Group
 
    Copyright (C) 2011 University of Oxford  */

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

#if !defined (Bingham_Watson_approx_h)
#define Bingham_Watson_approx_h

#include <iostream>
#include <fstream>
#include <iomanip>
#define WANT_STREAM
#define WANT_MATH
#include <string>
#include "miscmaths/miscmaths.h"
#include "miscmaths/nonlin.h"
#include "stdlib.h"


using namespace NEWMAT;
using namespace MISCMATHS;


  #ifndef M_PI
  #define M_PI 3.14159265358979323846
  #endif

  #define INV3 0.333333333333333   // 1/3
  #define SQRT3 1.732050807568877 //sqrt(3)
  #define INV54 0.018518518518519 // 1/54
 
  #define min3(a,b,c)  (a < b ? min(a,c) : min(b,c) )
 

  //Saddlepoint approximation of confluent hypergeometric function of a matrix argument B
  //Vector x has the eigenvalues of B.
  float hyp_Sapprox(ColumnVector& x);

 
  //Saddlepoint approximation of confluent hypergeometric function of a matrix argument, with its eigenvalues being l1,l1,l2 or l1,l2,l2 with l1!=l2.
  //Vector x has the three eigenvalues. This function can be also used to approximate a confluent hypergeometric function of a scalar argument k
  //by providing x=[k 0 0].
  float hyp_Sapprox_twoequal(ColumnVector& x);

 
  //Saddlepoint approximation of the ratio of two hypergeometric functions, with matrix arguments L and B (3x3). Vectors xL & xB contain the eigenvalues of L and B.
  //Used for the ball & Binghams model. 
  float hyp_SratioB(ColumnVector& xL,ColumnVector& xB);


  //Saddlepoint aproximation of the ratio ot two hypergeometric functions with matrix arguments L and B in two steps: First denominator, then numerator. 
  //This allows them to be updated independently, used for the ball & Binghams model to compute the likelihood faster.
  //This function returns values used in the denominator approximation. xB containes the two non-zero eigenvalues of matrix B.
  ReturnMatrix approx_denominatorB(ColumnVector& xB);


  //Second step for saddlepoint approximation of the ratio of two hypergeometric functions with matrix arguments L and B (xL has the eigenvalues of L). 
  //Assume that the denominator has already been approximated by the function above and the parameters are stored in denomvals.
  //Here approximate the numerator and return the total ratio approximation.
  float hyp_SratioB_knowndenom(ColumnVector &xL,ColumnVector& denomvals);


  //Saddlepoint approximation of the ratio of two hypergeometric functions, one with matrix argument L and another with scalar argument k. Vector xL contains the eigenvalues of L.
  //Used for the ball & Watsons model. 
  float hyp_SratioW(ColumnVector& xL,const double k);


  //Saddlepoint aproximation of the ratio ot two hypergeometric functions, one with matrix arguments L and the other with scalar argument k in two steps: 
  //First denominator, then numerator. This allows them to be updated independently, used for the ball & Watsons model to compute the likelihood faster.
  //This function returns values used in the denominator approximation.
  ReturnMatrix approx_denominatorW(const double k);


  //Second step for saddlepoint approximation of the ratio of two hypergeometric functions, with matrix argument L and scalar argument k (xL has the eigenvalues of L). 
  //Assume that the denominator has already been approximated by the function above and the parameters are stored in denomvals.
  //Here approximate the numerator and return the total ratio approximation.
  float hyp_SratioW_knowndenom(ColumnVector &xL,ColumnVector& denomvals);


  //Using the values of vector x, construct a qubic equation and solve it analytically.
  //Solution used for the saddlepoint approximation of the confluent hypergeometric function with matrix argument B (3x3) (See Kume & Wood, 2005)
  //Vector x contains the eigenvalues of B. 
  float find_t(const ColumnVector& x);

  //cubic root
  float croot(const float& x);


#endif
