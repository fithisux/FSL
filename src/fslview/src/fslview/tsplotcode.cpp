#include "tsplotcode.h"
#include <math.h>



void TsPlotCode::preWhitenModel(const ColumnVector& ac,
                                const Matrix& designMatrix,
                                Matrix& preWhitenedMatrix)
{          
  int nevs = designMatrix.Ncols();
  int npts = designMatrix.Nrows();

  preWhitenedMatrix = designMatrix;
  
  int zeropad = (int)pow(2,ceil(log(npts)/log(2)));

  // get prewhitening filter from ac
  ColumnVector pwfilter;

  establishPwFilter(ac, pwfilter, zeropad, npts);

  // prewhiten each of the evs
  for(int ev = 1; ev <= nevs; ev++)
    {
      ColumnVector pwdata;
      ColumnVector test = designMatrix.Column(ev);
      preWhitenData(designMatrix.Column(ev), 
                    pwdata, pwfilter, zeropad, npts);
      preWhitenedMatrix.Column(ev) = pwdata;
    }
}

void TsPlotCode::establishPwFilter(const ColumnVector& ac, 
                                   ColumnVector& pwfilter, 
                                   int zeropad, int npts)
{
 // FFT auto corr estimate
  ColumnVector dummy(zeropad);
  dummy = 0;
  ColumnVector vrow(zeropad);
  vrow = 0;

  // ac maybe cutoff to be smaller than npts
  int sizeTS = ac.Nrows();
  if (sizeTS > npts/2)
    sizeTS = npts/2;
  
  vrow.Rows(1,sizeTS) = ac.Rows(1,sizeTS);
  vrow.Rows(zeropad - sizeTS + 2, zeropad) = ac.Rows(2, sizeTS).Reverse();

  ColumnVector ac_fft_imag;

  FFT(vrow, dummy, pwfilter, ac_fft_imag);      

  // inverse auto corr to give prewhitening filter
  // no DC component so set first value to 0
  pwfilter(1) = 0.0;
  
  for(int j = 2; j <= zeropad; j++)
    {
      if (pwfilter(j)<0)
	{
	  cout << "Warning: possible high autocorrelation in time series" << endl;
	  pwfilter(j)=0;
	}
      else
		  pwfilter(j) = 1.0/::sqrt(pwfilter(j));	      
    }  

  // normalise pwfilter such that sum(j)((pwfilter)^2/zeropad)) = 1
  pwfilter /= ::sqrt(pwfilter.SumSquare()/zeropad);
}

void TsPlotCode::preWhitenData(const ColumnVector& data, ColumnVector& pwdata, ColumnVector& pwfilter, int zeropad, int npts)
{ 

 ColumnVector data_fft_real, data_fft_imag, realifft, dummy;
  dummy.ReSize(zeropad);
  dummy = 0;

  // Remove and store mean
  float mn = MISCMATHS::mean(data).AsScalar();
  pwdata.ReSize(zeropad);
  pwdata = 0;


  pwdata.Rows(1,npts) = data - mn;

  // FFT data
  FFT(pwdata, dummy, data_fft_real, data_fft_imag);
  FFTI(SP(pwfilter, data_fft_real), SP(pwfilter, data_fft_imag), realifft, dummy); 

  // take first npts and restore mean
  pwdata = realifft.Rows(1,npts) + mn;  
  
}






TimeSeries::Handle TsPlotCode::preWhitenTimeseries(const ColumnVector& ac, 
                                   TimeSeries::Handle& ts)
{
  int npts = ts->inqVolCount();
  
  ColumnVector  tsCv(convertTimeSeries(ts));

  int zeropad = (int)pow(2,ceil(log(npts)/log(2)));
  
  // get prewhitening filter from ac
  ColumnVector pwfilter;
  establishPwFilter(ac, pwfilter, zeropad, npts);

  // prewhiten
  ColumnVector pwts;
  preWhitenData(tsCv, pwts, pwfilter, zeropad, npts);

  return convertColumnVector(pwts,ts->inqX(),ts->inqY(),ts->inqZ());
}
  

TimeSeries::Handle TsPlotCode::convertColumnVector(const ColumnVector& cv,
                                           short x, short y, short z)
{
  TimeSeries::Handle modelTs = TimeSeriesD::create(x,y,z,cv.Nrows());
  for(int n = 0;n < cv.Nrows();n++)
    modelTs->setValue(n,cv(n+1));

  return modelTs;
}

TimeSeries::Handle TsPlotCode::convertMatrix(const Matrix& mat, int col,
                                           short x, short y, short z)
{
  TimeSeries::Handle modelTs = TimeSeriesD::create(x,y,z,
                                                  mat.Nrows());
  for(int n = 0;n < mat.Nrows();n++)
    modelTs->setValue(n,mat(n+1,col));

  return modelTs;
}

ColumnVector TsPlotCode::convertTimeSeries(const TimeSeries::Handle& ts)
{
  ColumnVector colVec(ts->inqVolCount());
  
  for(int j = 0;j < ts->inqVolCount();++j)
    colVec(j+1) = ts->value(j);

  return colVec;
}
