/*  histogram.cc

    Mark Woolrich, FMRIB Image Analysis Group

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

#include "miscmaths.h"
#include "histogram.h"

using namespace std;

#ifndef NO_NAMESPACE
namespace MISCMATHS {
#endif

  void Histogram::generate()
    {
      Tracer ts("Histogram::generate");

      int size = sourceData.Nrows();
      
      if(calcRange)
	{
	  // calculate range automatically
	  histMin=histMax=sourceData(1);
	  for(int i=1; i<=size; i++)
	    {
	      if (sourceData(i)>histMax)
		histMax=sourceData(i);
	      if (sourceData(i)<histMin)
		histMin=sourceData(i);
	    }
	}
      
      // zero histogram
      histogram.ReSize(bins);
      histogram=0;
      
      // create histogram; the MIN is so that the maximum value falls in the
      // last valid bin, not the (last+1) bin
      for(int i=1; i<=size; i++)
	{
	   histogram(getBin(sourceData(i)))++;
	}
    }

  void Histogram::smooth()
    {
      Tracer ts("Histogram::smooth");

      ColumnVector newhist=histogram;

      // smooth in i direction
      newhist=0;
      ColumnVector kernel(3); 
      // corresponds to Gaussian with sigma=0.8 voxels
      //       kernel(1)=0.5;
      //       kernel(2)=0.2283;      
      //       kernel(3)=0.0219;
      // corresponds to Gaussian with sigma=0.6 voxels
      //       kernel(1)=0.6638;
      //       kernel(2)=0.1655;      
      //       kernel(3)=0.0026;

      //gauss(0.5,5,1)
      kernel(1)=0.7866;
      kernel(2)=0.1065;      
      kernel(3)=0.0003;

      for(int i=1; i<=bins; i++)
	  {
	    float val=0.5*histogram(i);
	    float norm=kernel(1);

	    if(i>1)
	      {
		val+=kernel(2)*(histogram(i-1));
		norm+=kernel(2);
	      }
	    if(i>2)
	      {
		val+=kernel(3)*(histogram(i-2));
		norm+=kernel(3);		
	      }
	    if(i<bins)
	      {
		val+=kernel(2)*(histogram(i+1));
		norm+=kernel(2);
	      }
	    if(i<bins-1)
	      {
		val+=kernel(3)*(histogram(i+2));
		norm+=kernel(3);		
	      }
	    val/=norm;

	    newhist(i)=val;
	  }

      histogram=newhist;

    }

  int Histogram::integrate(float value1, float value2) const
    {
      int upperLimit = getBin(value2);
      int sum = 0;

      for(int i = getBin(value1)+1; i< upperLimit; i++)
	{
	  sum += (int)histogram(i);
	}
      return sum;
    }

  float Histogram::mode() const
    {
      int maxbin = 0;
      int maxnum = 0;

      for(int i = 1; i< bins; i++)
	{
	  if((int)histogram(i) > maxnum) {
	    maxnum = (int)histogram(i);
	    maxbin = i;
	  }
	}

      return getValue(maxbin);
    }

#ifndef NO_NAMESPACE
}
#endif




































