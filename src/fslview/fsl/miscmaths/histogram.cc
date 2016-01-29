/*  histogram.cc

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

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




































