/*  AutoCorrEstimator.h

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

#if !defined(__AutoCorrEstimator_h)
#define __AutoCorrEstimator_h

#include <iostream>
#include <fstream>
#define WANT_STREAM
#define WANT_MATH

#include "miscmaths/miscmaths.h"
#include "newimage/newimageall.h"
#include "fslsurface/fslsurface.h"

using namespace NEWMAT;
using namespace MISCMATHS;
using namespace fslsurface_name;

namespace FILM {
     
  class AutoCorrEstimator
    {
    public:
      AutoCorrEstimator(const Matrix& pxdata) : 
        sizeTS(pxdata.Nrows()),
	numTS(pxdata.Ncols()),
	xdata(pxdata),
	acEst(pxdata.Nrows(), pxdata.Ncols()),
	vrow(),
	xrow(),
	dummy(),
	realifft(), 
	dm_mn(),
	zeropad(0)
	{ 
	  zeropad = MISCMATHS::nextpow2(pxdata.Nrows());
	  vrow.ReSize(zeropad);
	  xrow.ReSize(zeropad);
	  dummy.ReSize(zeropad);
	  realifft.ReSize(zeropad);
	  countLargeE.ReSize(zeropad);
	  countLargeE = 0;
	  acEst = 0;
	  acEst.Row(1) = 1;
	}

      void calcRaw(int lag = 0);
      void spatiallySmooth(const ColumnVector& epivol, int masksize, double usanthresh, const NEWIMAGE::volume<float>& usan_vol, int lag=0);
      void spatiallySmooth(const fslSurface<float, unsigned int> surfaceData, const float sigma, const float extent, int lag=0); 
      void applyConstraints();
      void filter(const ColumnVector& filterFFT);
      Matrix fitAutoRegressiveModel();
      void pava();
      void preWhiten(const ColumnVector& in, ColumnVector& ret, int i, Matrix& dmret, bool highfreqremovalonly=false);
      void setDesignMatrix(const Matrix& dm);
      double establishUsanThresh(const ColumnVector& epivol);

      void getMeanEstimate(ColumnVector& ret);

      Matrix& getEstimates() { return acEst; }
      Matrix& getE() { return E; }
      ColumnVector& getCountLargeE(){ return countLargeE; }

      int getZeroPad() { return zeropad; }
      void tukey(int M);
      void multitaper(int M);
      int pacf(const ColumnVector& x, int minorder, int maxorder, ColumnVector& betas);
      
      NEWIMAGE::volume<float> mask;

    private:
      const int sizeTS;
      const int numTS;
      AutoCorrEstimator();
      const AutoCorrEstimator& operator=(AutoCorrEstimator&);
      AutoCorrEstimator(AutoCorrEstimator&);
      void getSlepians(int M, int sizeTS, Matrix& slepians);

      const Matrix& xdata;
      Matrix acEst;
      Matrix E;
      ColumnVector countLargeE;

      Matrix dminFFTReal;
      Matrix dminFFTImag;

      ColumnVector vrow;
      ColumnVector xrow;     

      ColumnVector dm_fft_real, dm_fft_imag;
      ColumnVector x_fft_real, ac_fft_real;
      ColumnVector x_fft_im, ac_fft_im;
     
      ColumnVector dummy;
      ColumnVector realifft;
      ColumnVector dm_mn;

      int zeropad;
    };
 
}

#endif


