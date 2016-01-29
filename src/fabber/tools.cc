/*  tools.cc - miscellaneous useful functions & classes

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

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

#include "tools.h"
#include "easylog.h"
#include <limits>

double DescendingZeroFinder::FindZero() const
{
    Tracer_Plus tr("DescendingZeroFinder::FindZero");
    
    double lower = searchMin;
    double upper = searchMax;
    double atLower, atUpper;
    
    double atSearchGuess = fcn(searchGuess);

    if (verbosity >= 2)
      LOG_ERR("IG: f(" << (searchGuess) 
	      << ") == " << atSearchGuess << endl);

    if (atSearchGuess < 0) 
    {
        upper = searchGuess;
        atUpper = atSearchGuess;
        atLower = fcn(lower);
	if (verbosity >= 2)
	  LOG_ERR("LG: f(" << (lower) << ") == " << atLower << endl);
        if (atLower <= 0)
            return lower;  // hit the limit
    }
    else
    {
        lower = searchGuess;
        atLower = atSearchGuess;
        atUpper = fcn(upper);
	if (verbosity >= 2)
	  LOG_ERR("UG: f(" << (upper) << ") == " << atUpper << endl);
        if (atUpper >= 0)
            return upper; // hit the limit
    }
    
    int evals = maxEvaluations - 2;
    double maxJump = searchScale;
    double guess, prevGuess = searchGuess;
    //    double nonlinearity = 10; // force a bisection first time
    
    // Interpolation only
    while ( (evals > 0) && 
	    (upper-lower > tolX || 
	     atLower-atUpper > tolY || 
	     upper/lower > ratioTolX ||
	     atUpper/atLower > ratioTolY) )
    {
        guess = guesstimator->GetGuess(lower, upper, atLower, atUpper);
	
	if (lower == guess || guess == upper)
	  {
	    // This should only happen if we're near the limits of
	    // double precision (constant factor probably depends on
	    // the guesstimator)
	    assert(upper-lower <= 2*fabs(lower)*numeric_limits<double>::epsilon());
	    LOG_ERR("DescendingZeroFinder: giving up without reaching tolerances because we're at the limits of double precision!\n");
	    LOG_ERR("Lower: f(" << lower << ") == " << atLower << endl <<
		    "Upper: f(" << upper << ") == " << atUpper << endl);
	    break;
	  }

	assert( lower < guess && guess < upper );
	//cout << lower << "<" << guess << "<" << upper << endl;
	if (!fcn.PickFasterGuess(&guess, lower, upper)) 
	  evals--; // only count the non-cached evaluations.
	//cout << lower << "<" << guess << "<" << upper << endl;
	assert(lower < guess && guess < upper);

        if ( guess-prevGuess > maxJump)
            guess = prevGuess + maxJump;
        else if ( guess-prevGuess < -maxJump )
            guess = prevGuess - maxJump;
            
        maxJump *= searchScaleGrowth;
        prevGuess = guess;
         
        // never mind -- already out of bounds
        // from checking the limits above.

	//        cout << "upper" << "\t" << "guess" << "\t" << "lower" << endl;
	//        cout << upper << "\t" << guess << "\t" << lower << endl;
        assert(guess > lower); assert(guess < upper);
        double atGuess = fcn(guess);
	//        cout << atUpper << "\t" << atGuess << "\t" << atLower << endl << endl;

	if (verbosity >= 2)
	  LOG_ERR("NG: f(" << (guess) << ") == " << atGuess << endl);


	//        double atGuessIfLinear = 
	  //            ( atLower+(atUpper-atLower)*(guess-lower)/(upper-lower) - atGuess );
	//        if (atGuess > atGuessIfLinear)
	  //            nonlinearity = (atGuess-atGuessIfLinear)/(atGuess-atUpper);
	//        else
	  //            nonlinearity = (atGuessIfLinear-atGuess)/(atLower-atGuess);
            
   
	//        cout << "Nonlinearity = " << nonlinearity << endl;
	//        cout << "atGuess = " << atGuess
	//             << "\natGuessIfLinear = " << atGuessIfLinear
	//             << "\natLower = " << atLower
	//             << "\natUpper = " << atUpper << endl;
       
        if (atGuess < 0) 
        {
            upper = guess;
            atUpper = atGuess;
        }
        else
        {
            lower = guess;
            atLower = atGuess;
        }
    } 

    /*    
    // One final interpolation -- not necessary, we could pick anything
    // between lower and upper really.
    guess = guesstimator->GetGuess(lower, upper, atLower, atUpper);
    assert( lower <= guess && guess <= upper );
    fcn.PickFasterGuess(&guess, lower, upper, true);
    assert( lower <= guess && guess <= upper );
    */

    // Pick either lower or upper bound, depending on which is closer to zero
    assert(atLower >= 0 && -atUpper >= 0);
    if (atLower < -atUpper)
      guess = lower;
    else
      guess = upper;

    if (verbosity >= 1)
      LOG_ERR("Final upper/lower ratio: " << (upper/lower) << endl);

    return guess;
}

double RiddlersGuesstimator::GetGuess(double lower, double upper, double atLower, double atUpper)
{
  Tracer_Plus tr("RiddlersGuesstimator::GetGuess");
  // equations below: from NRIC, section 9.2.  Simpler than Brent, slightly less reliable.

  if (halfDone)
    {
      // Phase two: fancy estimation.
      halfDone = false;
      double x3, fx3;

      if (x1 == lower)
	{
	  x3 = upper;
	  fx3 = atUpper;
	}
      else
	{
	  assert(x2 == upper);
	  x3 = lower;
	  fx3 = atLower;
	}
      
      assert(x1 < x3 && x3 < x2);
      assert(fx2 < fx1);

      Warning::IssueOnce("Riddler's Method; No special cases!");

      if (false) //x3 != (x1 + x2)/2)
	{
	  LOG_ERR("x3 == " << x3 << ", x1 == " << x1 << ", x2 = " << x2 << endl);
	  LOG_ERR("x3 - (x1+x2)/2 == " << x3 - (x1+x2)/2 << endl);
	  Warning::IssueAlways("Riddler's Method: x3 != (x1+x2)/2");
	}
      else if (true) //(fx2 < fx3 && fx3 < fx1)
	{
	  double s = (fx1-fx2 > 0) ? +1.0 : -1.0; // s = sign(fx1-fx2)
	  double x4 = x3 + (x3-x1)*s*fx3/sqrt(fx3*fx3-fx1*fx2);
      
	  Warning::IssueAlways("Riddler's Method: phase two");

	  assert(lower < x4 && x4 < upper);

	  return x4;
	}
      else
	{
	  Warning::IssueAlways("Riddler's Method cheat: dropping back to the bisection method!");
	}
    }
  
  halfDone = true;
  // Phase one: just pick the midpoint, but save the values for phase two
  x1 = lower;
  fx1 = atLower;
  x2 = upper;
  fx2 = atUpper;
  
  Warning::IssueAlways("Riddler's Method: phase one");
  
  double x3 = (x1 + x2)/2;
  return x3;

}
