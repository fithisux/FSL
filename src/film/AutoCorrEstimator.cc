/*  AutoCorrEstimator.cc

    Mark Woolrich and Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 1999-2008 University of Oxford  */

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
#include <sstream>
#define WANT_STREAM

#include "AutoCorrEstimator.h"
#include "utils/log.h"
#include "miscmaths/histogram.h"
#include "newimage/newimageall.h"
#include "glm.h"

using namespace Utilities;
using namespace NEWIMAGE;

namespace FILM {

  void AutoCorrEstimator::setDesignMatrix(const Matrix& dm) {
    Tracer tr("AutoCorrEstimator::setDesignMatrix");

    int numPars = dm.Ncols();
    
    dminFFTReal.ReSize(zeropad, numPars);
    dminFFTImag.ReSize(zeropad, numPars);
      
    ColumnVector dmrow;
    dmrow.ReSize(zeropad);

    ColumnVector dm_fft_real, dm_fft_imag;
    ColumnVector dummy(zeropad);
    ColumnVector realifft(zeropad);
    dm_mn.ReSize(numPars);
    dm_mn=0;

    // FFT design matrix
    for(int k = 1; k <= numPars; k++)
      {
	dummy = 0;
	dmrow = 0;
	dm_mn(k) = MISCMATHS::mean(ColumnVector(dm.Column(k))).AsScalar();
	dmrow.Rows(1,sizeTS) = dm.Column(k) - dm_mn(k);

	FFT(dmrow, dummy, dm_fft_real, dm_fft_imag);
	
	dminFFTImag.Column(k) = dm_fft_imag;
	dminFFTReal.Column(k) = dm_fft_real;
      } 
  }

  void AutoCorrEstimator::preWhiten(const ColumnVector& in, ColumnVector& ret, int i, Matrix& dmret, bool highfreqremovalonly) {
      Tracer tr("AutoCorrEstimator::preWhiten");

      int numPars = dminFFTReal.Ncols();

      ret.ReSize(sizeTS);
      dmret.ReSize(sizeTS, numPars);
      // FFT auto corr estimate
      dummy = 0;
      vrow = 0;
      vrow.Rows(1,sizeTS/2) = acEst.Column(i).Rows(1,sizeTS/2);
      vrow.Rows(zeropad - sizeTS/2 + 2, zeropad) = acEst.Column(i).Rows(2, sizeTS/2).Reverse();

      FFT(vrow, dummy, ac_fft_real, ac_fft_im);      

      float norm = ac_fft_real.SumSquare();
      
      // Compare with raw FFT to detect high frequency artefacts:
      bool violate = false;
      ColumnVector violators(zeropad);
      violators = 1;
      for(int j = 1; j <= zeropad; j++)
	{
	  if(highfreqremovalonly)
	    {
	      E(j,i) = sqrt(E(j,i)/((ac_fft_real(j)*ac_fft_real(j))/norm));
	      
	      // look for high frequency artefacts
	      if(E(j,i) > 4 && j > zeropad/4 && j < 3*zeropad/4)
		{	  
		  violate = true;
		  violators(j) = 0;
		  countLargeE(j) = countLargeE(j) + 1;
		}
	    }
	}

      // FFT x data
      dummy = 0;
      xrow = 0;
      xrow.Rows(1,sizeTS) = in;
     
      FFT(xrow, dummy, x_fft_real, x_fft_im);
      
      if(highfreqremovalonly)
	{
	  ac_fft_real = violators;
	}
      else
	{
	  // inverse auto corr to give prewhitening filter
	  // no DC component so set first value to 0
	  ac_fft_real(1) = 0.0;
      
	  for(int j = 2; j <= zeropad; j++)
	    {
	      ac_fft_real(j) = 1.0/sqrt(fabs(ac_fft_real(j)));	      
	    }
	}

      // normalise ac_fft such that sum(j)(ac_fft_real)^2 = 1
      ac_fft_real /= sqrt(ac_fft_real.SumSquare()/zeropad);

      // filter design matrix
      for(int k = 1; k <= numPars; k++)
	{
	  dm_fft_real = dminFFTReal.Column(k);
	  dm_fft_imag = dminFFTImag.Column(k);
	  FFTI(SP(ac_fft_real, dm_fft_real), SP(ac_fft_real, dm_fft_imag), realifft, dummy);

	  // place result into ret:
	  dmret.Column(k) = realifft.Rows(1,sizeTS) + dm_mn(k);
	  //float std = pow(MISCMATHS::var(ColumnVector(dmret.Column(k))),0.5);
	  //dmret.Column(k) = (dmret.Column(k)/std) + mn(k);
	}

      // Do filtering of data:
      FFTI(SP(ac_fft_real, x_fft_real), SP(ac_fft_real, x_fft_im), realifft, dummy);
      
      // place result into ret:
      ret = realifft.Rows(1,sizeTS);     

    }

  Matrix AutoCorrEstimator::fitAutoRegressiveModel()
    {
      Tracer trace("AutoCorrEstimator::fitAutoRegressiveModel");
      
      cout << "Fitting autoregressive model..." << endl;
      
      const int maxorder = 15;
      const int minorder = 1;
      
      // setup temp variables
      ColumnVector x(sizeTS);
      ColumnVector order(numTS);
      Matrix betas(maxorder, numTS);
      betas = 0;
      acEst.ReSize(sizeTS, numTS);
      acEst = 0;
      int co = 1;
       
      for(int i = 1; i <= numTS; i++)
	{
	  x = xdata.Column(i);
	  ColumnVector betastmp;

	  order(i) = pacf(x, minorder, maxorder, betastmp);

	  if(order(i) != -1)
	    {
	      // Calculate auto corr:
	      ColumnVector Krow(sizeTS);
	      Krow = 0;
	      Krow(sizeTS) = 1;
	      if(order(i) > 0)
		{
		  Krow.Rows(sizeTS-int(order(i)), sizeTS-1) = -betastmp.Rows(1,int(order(i))).Reverse();
		  betas.SubMatrix(1,int(order(i)),i,i) = betastmp.Rows(1,int(order(i)));
		}
	      
	      Matrix Kinv(sizeTS, sizeTS);
	      Kinv = 0;
  
	      for(int j = 1; j <= sizeTS; j++)
		{
		  if(order(i)==1)
		    {
		      float arone = betastmp(1);
		      for(int k = 1; k <= sizeTS; k++)
			{
			  Kinv(j,k) = MISCMATHS::pow(float(arone),int(abs(k-j)));
			}
		    }
		  else
		    Kinv.SubMatrix(j,j,1,j) = Krow.Rows(sizeTS-j+1,sizeTS).t();		  
		}

	      // Kinv now becomes V:
	      if(order(i)!=1)
		Kinv = (Kinv.t()*Kinv).i();

	      acEst.SubMatrix(1,sizeTS/2+1,i,i) = (Kinv.SubMatrix(sizeTS/2, sizeTS/2, sizeTS/2, sizeTS)/Kinv.MaximumAbsoluteValue()).AsColumn();
	      
	      if(co > 200)
		{
		  co = 1;
		  cout<< (float)i/(float)numTS << ",";
		}
	      else
		co++;
	    }
	}
      
      write_ascii_matrix(LogSingleton::getInstance().appendDir("order"), order);
      write_ascii_matrix(LogSingleton::getInstance().appendDir("betas"), betas);
      countLargeE = 0;
      cout<< " Completed" << endl; 
      return(betas);
    }

  int AutoCorrEstimator::pacf(const ColumnVector& x, int minorder, int maxorder, ColumnVector& betas)
    { 
      Tracer ts("pacf");
      int order = -1;

      // Set c
      Matrix c(1,1);
      c(1,1) = 1;
      
      Glm glm;

      for(int i = minorder+1; i <= maxorder+1; i++)
	{
	  ColumnVector y = x.Rows(i+1,sizeTS);

	  // Setup design matrix
	  Matrix X(sizeTS-i, i);
	  X = 0;
	  
	  for(int j = 1; j <= i; j++)
	    {
	      X.Column(j) = x.Rows(i+1-j,sizeTS-j).AsColumn();	    
	    }
	  
	  glm.setData(y, X, c);
	  
	  glm.ComputeResids();
	  betas = glm.Getb();
	  
	  if((abs(betas(i)) < (1/sizeTS) + (2/pow(sizeTS,0.5)) && order == -1) 
	     || i == maxorder+1)
          {      
	    order = i-1;
	    break;
	  }
	}
      return order; 
    }
 
  double AutoCorrEstimator::establishUsanThresh(const ColumnVector& epivol)
    {
      int num = epivol.Nrows();
      Histogram hist(epivol, max(num/200,1));
      hist.generate();
      float mode = hist.mode();
      cout<< "mode = " << mode << endl;

      double sum = 0.0;
      int count = 0;

      // Work out standard deviation from mode for values greater than mode:
      for(int i = 1; i <= num; i++) {
	if(epivol(i) > mode) {
	  sum += (epivol(i) - mode)*(epivol(i) - mode);
	  count++;
	}
      }

      double sig = pow(sum/num, 0.5);
      cout<< "sig = " << static_cast<int>(sig) << endl;
      return sig/3;
    } 

  void AutoCorrEstimator::spatiallySmooth(const ColumnVector& epivol, int masksize, double usan_thresh, const volume<float>& usan_vol, int lag) {
    Tracer trace("AutoCorrEstimator::spatiallySmooth");
    if(numTS<=1)
	cerr << "Warning: Number of voxels = " << numTS << ". Spatial smoothing of autocorrelation estimates is not carried out" << endl;
    else
      {
	
	if(lag==0)
	  lag = MISCMATHS::Min(40,int(sizeTS/4));

	if(usan_thresh == 0) usan_thresh = establishUsanThresh(epivol); // Establish epi thresh to use:
	volume4D<float> susan_vol(mask.xsize(),mask.ysize(),mask.zsize(),1);
	volume<float> usan_area(mask.xsize(),mask.ysize(),mask.zsize());
	volume<float> kernel;
	kernel = gaussian_kernel3D(masksize,mask.xdim(),mask.ydim(),mask.zdim(),2.0);	
	int factor = 10000;
	cout<< "Spatially smoothing auto corr estimates" << endl;
	
	for(int i=2 ; i <= lag; i++)
	  {
	    // setup susan input
	    susan_vol.setmatrix(acEst.Row(i),mask); 
	    susan_vol*=factor;
	    susan_vol[0]=susan_convolve(susan_vol[0],kernel,1,0,1,&usan_area,usan_vol,usan_thresh*(float)usan_thresh);
	    // insert output back into acEst
            susan_vol/=factor;
	    acEst.Row(i)=susan_vol.matrix(mask);
	    cout<< ".";
	  }
	
	cout<< endl << "Completed" << endl;
      }
  }

  void AutoCorrEstimator::spatiallySmooth(fslSurface<float, unsigned int> surfaceData, const float sigma, const float extent,int lag) {
   Tracer trace("AutoCorrEstimator::spatiallySmooth_GIFTI");
   if(numTS<=1)
     cerr << "Warning: Number of vertices = " << numTS << ". Spatial smoothing of autocorrelation estimates is not carried out" << endl;
   else {	
     if(lag==0)
       lag = MISCMATHS::Min(40,int(sizeTS/4));
     cout<< "Spatially smoothing auto corr estimates for surface" << endl;
     for(int i=2 ; i <= lag; i++) {
       vector<float> scalars;
       for ( unsigned int point = 1; point<=surfaceData.getNumberOfVertices(); point++ )
	 scalars.push_back( acEst(i,point) );
       surfaceData.clearScalars();
       surfaceData.insertScalars(scalars,0,"input field" );  //Add data to surface
       sc_smooth_gaussian_geodesic(  surfaceData , 0, sigma, extent , false);
       int point(1);
       for (std::vector< float >::const_iterator i_sc = surfaceData.const_scbegin(0); i_sc != surfaceData.const_scend(0); i_sc++, point++)
	 acEst(i,point) = *i_sc; // insert output back into acEst
       cout<< ".";
     }
     cout<< endl << "Completed" << endl;
   }
}

  
  void AutoCorrEstimator::calcRaw(int lag) { 
    
    cout<< "Calculating raw AutoCorrs...";      

    MISCMATHS::xcorr(xdata, acEst, lag, zeropad);

    cout<< " Completed" << endl;  
  }
  
  void AutoCorrEstimator::filter(const ColumnVector& filterFFT) {

    Tracer tr("AutoCorrEstimator::filter");

    cout<< "Combining temporal filtering effects with AutoCorr estimates... ";

    // This function adjusts the autocorrelations as if the
    // xdata has been filtered by the passed in filterFFT
    // DOES NOT filter the xdata itself
    ColumnVector vrow;

    // make sure p_vrow is cyclic (even function)
    vrow.ReSize(zeropad);

    ColumnVector fft_real;
    ColumnVector fft_im;
    ColumnVector dummy(zeropad);    
    ColumnVector realifft(zeropad);

    for(int i = 1; i <= numTS; i++)
      {
	dummy = 0;
	vrow = 0;
	vrow.Rows(1,sizeTS/2) = acEst.Column(i).Rows(1,sizeTS/2);
	vrow.Rows(zeropad - sizeTS/2 + 2, zeropad) = acEst.Column(i).Rows(2, sizeTS/2).Reverse();
      
	FFT(vrow, dummy, fft_real, fft_im);

	FFTI(SP(fft_real, filterFFT), dummy, realifft, dummy);
	
	// place result into acEst:
	acEst.Column(i) = realifft.Rows(1,sizeTS)/realifft(1);
      }

    cout<< " Completed" << endl;
    
  }

  void AutoCorrEstimator::multitaper(int M) {
    Tracer tr("AutoCorrEstimator::multitaper");
    
    cout<< "Multitapering... ";

    Matrix slepians;
    getSlepians(M, sizeTS, slepians);

    //LogSingleton::getInstance().out("slepians", slepians, false);
    
    ColumnVector x(zeropad);
    x = 0;
    ColumnVector fft_real;
    ColumnVector fft_im;
    ColumnVector dummy(zeropad);
    ColumnVector dummy2;
    ColumnVector realifft(zeropad);
    dummy = 0;
    
    Matrix Sk(zeropad, slepians.Ncols());
    acEst.ReSize(sizeTS, numTS);
    acEst = 0;

    for(int i = 1; i <= numTS; i++) 
      {
	// Compute FFT for each slepian taper
	for(int k = 1; k <= slepians.Ncols(); k++) 
	  {
	    x.Rows(1,sizeTS) = SP(slepians.Column(k), xdata.Column(i));
	   
	    FFT(x, dummy, fft_real, fft_im);
	    for(int j = 1; j <= zeropad; j++)
	    {
	      // (x+iy)(x-iy) = x^2 + y^2
	      fft_real(j) = fft_real(j)*fft_real(j) + fft_im(j)*fft_im(j);
	      Sk(j,k) = fft_real(j);
	    }
	  }

	// Pool multitaper FFTs
	fft_im = 0;
	for(int j = 1; j <= zeropad; j++)
	  {
	    fft_real(j) = MISCMATHS::mean(ColumnVector(Sk.Row(j).t())).AsScalar();
	    
	  }

	// IFFT to get autocorr
	FFTI(fft_real, fft_im, realifft, dummy2);
	//LogSingleton::getInstance().out("Sk", Sk, false);
	//LogSingleton::getInstance().out("realifft", realifft);
	//LogSingleton::getInstance().out("fftreal", fft_real);
	
	float varx = MISCMATHS::var(ColumnVector(x.Rows(1,sizeTS))).AsScalar();
	acEst.Column(i)=realifft.Rows(1,sizeTS)/varx;
      }
    countLargeE = 0;
    cout<< "Completed" << endl;
  }

  void AutoCorrEstimator::getSlepians(int M, int sizeTS, Matrix& slepians) {
    
    Tracer tr("AutoCorrEstimator::getSlepians");
    slepians.ReSize(sizeTS, 2*M);
    
    ifstream in;
    
    ostringstream osc;
    osc << sizeTS << "_" << M;
	
    string fname("/usr/people/woolrich/parads/mt_" + osc.str());
    in.open(fname.c_str(), ios::in);
    if(!in)
      throw Exception("Multitapering: Slepians file not found");
   
    for(int j = 1; j <= sizeTS; j++) 
      {
	for(int i = 1; i <= 2*M; i++) 
	  {
	    in >> slepians(j,i);
	  }
      }	

    in.close();
  }

  void AutoCorrEstimator::tukey(int M) {
    
    Tracer tr("AutoCorrEstimator::tukey");

    cout<< "Tukey M = " << M << endl;

    cout<< "Tukey estimates... ";
	
    ColumnVector window(M);

    for(int j = 1; j <= M; j++)
      {
	window(j) = 0.5*(1+cos(M_PI*j/(float(M))));
      }
    
    for(int i = 1; i <= xdata.Ncols(); i++) {
	
	acEst.SubMatrix(1,M,i,i) = SP(acEst.SubMatrix(1,M,i,i),window);
	acEst.SubMatrix(M+1,sizeTS,i,i) = 0;
	
    }
    countLargeE = 0;
    cout<< "Completed" << endl;
  }

  void AutoCorrEstimator::pava() {
    
    Tracer tr("AutoCorrEstimator::pava");
    
    cout<< "Using New PAVA on AutoCorr estimates... ";

    for(int i = 1; i <= numTS; i++) {
	int stopat = (int)sizeTS/2;

	// 5% point of distribution of autocorr about zero
	const float th = (-1/sizeTS)+(2/sqrt(sizeTS));

	ColumnVector values = acEst.Column(i);
	ColumnVector zero(1);
	zero = 0;
	values = values.Rows(1,stopat) & zero;

	ColumnVector gm(stopat + 1);
	for(int j = 1; j <= stopat + 1; ++j)
	  gm(j) = j;
	  
	ColumnVector weights(stopat+1);
	weights = 1;

	bool anyviolators = true;

	while(anyviolators) {
	  anyviolators = false;
	  
	  for(int k = 2; k <= values.Nrows(); k++) {
	    if(values(k) > values(k-1)) {
	      anyviolators = true;
	      values(k-1) = (values(k-1)*weights(k-1) + values(k)*weights(k))/(weights(k-1) + weights(k));
	      values = values.Rows(1,k-1) & values.Rows(k+1,values.Nrows());
	      weights(k-1) = weights(k) + weights(k-1);
	      weights = weights.Rows(1,k-1) & weights.Rows(k+1,weights.Nrows());

	      for(int j = 1; j <= stopat + 1; j++) {
		if(gm(j) >= k)  
		  gm(j) = gm(j)-1;
	      }

	      break;
	    }
	  }
	}
	
	acEst.Column(i) = 0.0;
	int j=1;

	for(; j <= stopat; j++) {
	  
	  acEst(j,i) = values(int(gm(j)));
	  if(acEst(j,i) <= 0.0)
	    {
	      acEst(j,i) = 0.0;
	      break;
 	    }
	}
	
	if(acEst(2,i) < th/2)
	{
	acEst.SubMatrix(2,stopat,i,i) = 0;
	}

	else if(j > 2)
	  //if(j > 2)  
	  {
	    int endst = j;
	    int stst = j-(int)(1+(j/8.0));

	    const int expwidth = MISCMATHS::Max((endst - stst)/2,1);
	    const int exppow = 2;
	
	    for(j = stst; j <= endst; j++)
	      {
		acEst(j,i) = acEst(j,i)*exp(-MISCMATHS::pow((j-stst)/float(expwidth),int(exppow)));
	      }	
	  }
	
    }
    countLargeE = 0;

    cout<< " Completed" << endl;
  }
  
  void AutoCorrEstimator::applyConstraints() {

    Tracer tr("AutoCorrEstimator::applyConstraints");

    cout<< "Applying constraints to AutoCorr estimates... ";

    for(int i = 1; i <= numTS; i++)
      {
	int j = 3;
	int stopat = (int)sizeTS/4;

	// found1 is last valid value above threshold
	int found1 = stopat;

	// 5% point of distribution of autocorr about zero
	const float thresh = (-1/sizeTS)+(2/sqrt(sizeTS));

	acEst(2,i) = (acEst(2,i)+ acEst(3,i))/2;
	if(acEst(2,i) < 0)
	  {
	    acEst(2,i) = 0;
	  }
	else
	  {
	    float grad = 0.0;

	    while(j <= stopat && j < found1 + 2)
	      {    
		grad = ((acEst(j,i) + acEst(j+1,i))/2 - acEst(j-1,i))/1.5;
		if(grad < 0)
		  acEst(j,i) = grad + acEst(j-1,i);
		else
		  acEst(j,i) = acEst(j-1,i);
		
		// look for threshold
		if(acEst(j,i) < thresh/3.0 && found1 == stopat)
		  {
		    found1 = j;
		  }
	   
		if(acEst(j,i) < 0)
		  {
		    acEst(j,i) = 0;
		  }
		
		j++;
	      }
	  }

	// set rest to zero:
	for(; j <= sizeTS; j++)
	  {
	      acEst(j,i) = 0;
	  }
      }
    cout<< "Completed" << endl;
  }

  void AutoCorrEstimator::getMeanEstimate(ColumnVector& ret)
    {
      Tracer tr("AutoCorrEstimator::getMeanEstimate");

      ret.ReSize(acEst.Nrows());
      // Calc global Vrow:
      for(int i = 1; i <= acEst.Nrows(); i++)
	{
	  ret(i) = MISCMATHS::mean(ColumnVector(acEst.Row(i).AsColumn())).AsScalar();
	}
    }
}












