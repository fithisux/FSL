/*  histogram.h

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

#if !defined(__histogram_h)
#define __histogram_h

#include <iostream>
#include <fstream>
#define WANT_STREAM
#define WANT_MATH

#include "newmatap.h"
#include "newmatio.h"
#include "miscmaths.h"

using namespace NEWMAT;

namespace MISCMATHS {
 
  class Histogram
    {
    public:
      Histogram(){};
      const Histogram& operator=(const Histogram& in){
	sourceData=in.sourceData; calcRange=in.calcRange; histMin=in.histMin; histMax=in.histMax; bins=in.bins;
	return *this;
      }

      Histogram(const Histogram& in){*this=in;}

      Histogram(const ColumnVector& psourceData, int numBins)
	: sourceData(psourceData), calcRange(true), bins(numBins){}

      Histogram(const ColumnVector& psourceData, float phistMin, float phistMax, int numBins) 
	: sourceData(psourceData), calcRange(false), histMin(phistMin), histMax(phistMax), bins(numBins){}
      
      void set(const ColumnVector& psourceData, int numBins) {	
	sourceData=psourceData; calcRange=true; bins=numBins;
      }

      void set(const ColumnVector& psourceData, float phistMin, float phistMax, int numBins) {	
	sourceData=psourceData; calcRange=false; histMin=phistMin; histMax=phistMax; bins=numBins;
      }

      void generate();
      
      float getHistMin() const {return histMin;}
      float getHistMax() const {return histMax;}
      void setHistMax(float phistMax) {histMax = phistMax;}
      void setHistMin(float phistMin) {histMin = phistMin;}
      void smooth();

      int integrateAll() {return sourceData.Nrows();}

      const ColumnVector& getData() {return histogram;}
      void setData(const ColumnVector& phist) { histogram=phist;}

      int integrateToInf(float value) const { return integrate(value, histMax); }
      int integrateFromInf(float value) const { return integrate(histMin, value); }
      int integrate(float value1, float value2) const;

      float mode() const;

      int getBin(float value) const;
      float getValue(int bin) const;

    protected:

    private:

      ColumnVector sourceData;
      ColumnVector histogram;

      bool calcRange;

      float histMin;
      float histMax;

      int bins; // number of bins in histogram
    };

  inline int Histogram::getBin(float value) const
    {
      float binwidth=(histMax-histMin)/bins;
      return Max(1, Min((int)((((float)bins)*((float)(value-(histMin-binwidth))))/((float)(histMax-histMin))),bins));
    }
  
  inline float Histogram::getValue(int bin) const
    {
      return (bin*(histMax-histMin))/(float)bins + histMin;
    }

}

#endif





