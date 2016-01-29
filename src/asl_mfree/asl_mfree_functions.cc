/*   asl_mfree_functions.cc various functions for the model-free ASL

      Michael Chappell - IBME and FMIRB Image Analysis Group

      Copyright (C) 2010 University of Oxford */

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

#include "asl_mfree_functions.h"

namespace OXASL {

  ReturnMatrix SVDdeconv(const Matrix& data, const Matrix& aif, float dt) {
    // do a singular value deconvolution of the data to get residue function

    int nti = data.Nrows();
    int nvox = data.Ncols();
    float truncfac = 0.2;

    // voxelwise SVD deconvolution
    Matrix aifconv; Matrix residue(nti,nvox);
    DiagonalMatrix S; DiagonalMatrix D; Matrix U; Matrix V;
    for (int v=1; v<=nvox; v++) {
      //make convolution matrix
      aifconv = dt * convmtx(aif.Column(v));
      //SVD
      SVD(aifconv,S,U,V);
      // invert the singular values
      D = S.i();
      // truncate (zero all singular values below threshold)
      for (int i=2; i<=D.Nrows(); i++) {
	if (S(i,i) < truncfac*S(1,1)) D(i,i)=0;
      }
      // calculate resdiue
      residue.Column(v) = V*D*U.t()*data.Column(v);
    }
      
    return residue; 
  }

  ReturnMatrix convmtx(const ColumnVector& invec){
    // create a (simple) convolution matrix

    int nentry = invec.Nrows();
    Matrix cmat(nentry,nentry);
    cmat=0.0;
    for (int i=1; i<=nentry; i++) {
      cmat.SubMatrix(i,i,1,i) = ((invec.Rows(1,i)).Reverse()).AsRow();
    }

    return cmat;
  }

  ReturnMatrix SVDdeconv_circular(const Matrix& data, const Matrix& aif, float dt) {
    // do a singular value deconvolution of the data to get residue function

    int nti = data.Nrows();
    int nvox = data.Ncols();
    float truncfac = 0.2;

    int nextra = floor(nti*1.2);
    ColumnVector padding(nextra);
    padding = 0.0;

    // voxelwise SVD deconvolution
    Matrix aifconv; Matrix residue(nti+nextra,nvox);
    DiagonalMatrix S; DiagonalMatrix D; Matrix U; Matrix V;
    for (int v=1; v<=nvox; v++) {
      //make convolution matrix
      aifconv = dt * convmtx_circular(aif.Column(v) & padding);
      //SVD
      SVD(aifconv,S,U,V);
      // invert the singular values
      D = S.i();
      // truncate (zero all singular values below threshold)
      for (int i=2; i<=D.Nrows(); i++) {
	if (S(i,i) < truncfac*S(1,1)) D(i,i)=0;
      }
      // calculate resdiue
      residue.Column(v) = V*D*U.t()*(data.Column(v) & padding);
    }
      
    residue = residue.Rows(1,nti); // only keep the number of timepoints that we have in the input.

    return residue; 
  }

 ReturnMatrix SVDdeconv_wu(const Matrix& data, const Matrix& aif, float dt) {
    // do a singular value deconvolution of the data to get residue function
   // circular convolution matrix
   // OI index
   // after Wu MRM 2003.

    int nti = data.Nrows();
    int nvox = data.Ncols();
    float oi_thresh = 0.1;

    int nextra = floor(nti*1.2);
    ColumnVector padding(nextra);
    padding = 0.0;

    // voxelwise SVD deconvolution
    Matrix aifconv; Matrix residue(nti+nextra,nvox);
    ColumnVector resid(nti+nextra);
    DiagonalMatrix S; DiagonalMatrix D; Matrix U; Matrix V;
    for (int v=1; v<=nvox; v++) {
      //make convolution matrix
      aifconv = dt * convmtx_circular(aif.Column(v) & padding);
      //SVD
      SVD(aifconv,S,U,V);
      // invert the singular values
      D = S.i();

      //first try with all singular values
      resid = V*D*U.t()*(data.Column(v) & padding);
      float oi;
      oi = 1/(nti*resid.Maximum())* (resid.Rows(3,nti+nextra)-2*resid.Rows(2,nti+nextra-1)+resid.Rows(1,nti+nextra-2)).SumAbsoluteValue();

      // start removing singular values one by one until we get the OI we want
      int i=nti+nextra;
      while (oi>oi_thresh & i>1) {
	D(i,i) = 0;
	resid = V*D*U.t()*(data.Column(v) & padding);
	oi = 1/((nti+nextra)*resid.Maximum())* (resid.Rows(3,nti+nextra)-2*resid.Rows(2,nti+nextra-1)+resid.Rows(1,nti+nextra-2)).SumAbsoluteValue();
	i--;
      }

      residue.Column(v) = resid;
    }


    return residue; 
  }

  ReturnMatrix convmtx_circular(const ColumnVector& invec){
    // create a (simple) convolution matrix

    int nentry = invec.Nrows();
    Matrix cmat(nentry,nentry);
    cmat=0.0;
    for (int i=1; i<=nentry; i++) {
      cmat.SubMatrix(i,i,1,i) = ((invec.Rows(1,i)).Reverse()).AsRow();
      cmat.SubMatrix(i,i,i+1,nentry) = ((invec.Rows(i+1,nentry)).Reverse()).AsRow();
    }

    return cmat;
  }

  void Deconv(const Matrix& data, const Matrix& aif, float dt, ColumnVector& mag, Matrix& resid) {
    //perform deconvolution and output the magnitude and residue function

    //do the deconvolution
    resid = SVDdeconv_wu(data,aif,dt);
    // extract magntiude and residue separately
    int nvox = data.Ncols();
    mag.ReSize(nvox);
    for (int v=1; v<=nvox; v++) {
      mag(v) = (resid.Column(v)).Maximum();
      resid.Column(v) /= mag(v);
    }

  }

  void BootStrap(const Matrix& aif, const Matrix& data, float dt, const ColumnVector& mag, const Matrix& resid, int Nwb, ColumnVector& magsd) {
    // do wild boot strapping estimation of CBF precision - returns the s.d. estimate
    int nvox = data.Ncols();
    int ntpts = data.Nrows();

    // calculate the residuals
    Matrix modelfit(data);
    Matrix aifconv; ColumnVector mfittemp;
    int nextra = resid.Nrows() - aif.Nrows();
    ColumnVector padding(nextra);
    padding = 0.0;
    for (int v=1; v<=nvox; v++) {
      //make convolution matrix
      aifconv = dt * convmtx_circular(aif.Column(v) & padding);
      // calculate resdiue
      mfittemp = aifconv*resid.Column(v);
      modelfit.Column(v) = mfittemp.Rows(1,ntpts);
    }
    Matrix residuals(data);
    residuals = data - modelfit;

    // wild bootstrapping
    ColumnVector Radevec(ntpts);
    Matrix WBdata(data);
    Matrix estresid(resid);
    Matrix WBresiduals(residuals);
    Matrix magdist(Nwb,nvox);
    cout << "WB step (of " << Nwb << "): ";
    for (int b=1; b<=Nwb; b++) {
      cout << b << " " << flush;
      // sample from the Rademacher distrbution
      Radevec=1.0;
      for (int i=1; i<=ntpts; i++) {
	if (rand() > RAND_MAX/2) { Radevec(i)=-1.0; }
      }

      // apply this to the residuals to get WB residuals
      for (int v=1; v<=nvox; v++) {
	WBresiduals.Column(v) = SP(residuals.Column(v),Radevec);
      }
      // make the WB data
      WBdata = modelfit + WBresiduals;

      // do the deconvolution to get magntiude estimates
      estresid = SVDdeconv_wu(WBdata,aif,dt);
      // extract magntiude and residue separately

      for (int v=1; v<=nvox; v++) {
	magdist(b,v) = (estresid.Column(v)).Maximum();
      }
    }
    

    magsd.ReSize(nvox);
    magsd = (stdev(magdist)).t();
    cout << endl << "WB done" << endl;
  }

  void Prepare_AIF(volume4D<float>& aif, const volume<float>& metric, const volume<float>& mask, const float mthresh) {
    // Prepeare the AIF by retaining only the AIF that meet our metric and filling in with the nearest good AIF if necessary
    int nx = aif.xsize();
    int ny = aif.ysize();
    int nz = aif.zsize();

    float bestdist; float dist;
    for (int x=0; x<=nx; x++) { for (int y=0; y<=ny; y++) { for (int z=0; z<=nz; z++) {
      // if metric is not above threshold then search for nearest voxel where it is satisfied
	  //cout << "Voxel: " << x << " " << y << " " << z << ": " << metric(x,y,z) << endl;
	  if ( mask(x,y,z)>0 ) {
	  if (metric(x,y,z) < mthresh) {
	bestdist = 1e12;
	Matrix aifcand;
	for (int sx=0; sx<=nx; sx++) { for (int sy=0;sy<=ny; sy++) { for (int sz=0; sz<=nz; sz++) {
	      // cout << sx << " " << sy << " " << sz << endl;
	      if ( (mask(x,y,z)>0) && (metric(sx,sy,sz) >= mthresh) ) {
		dist = pow(sx-x,2.0) + pow(sy-y,2.0) + pow(sz-z,2.0); //strictly this is distance squared
		//cout << dist << endl;
	    if (dist < bestdist) {
	      aifcand = aif.voxelts(sx,sy,sz); 
	      //cout << "New cand: " << aifcand.t() << endl;
	      bestdist=dist;
	    }
	    else if (dist == bestdist) {
	      aifcand |= aif.voxelts(sx,sy,sz);
	      //cout << aifcand.t() << endl;
	    }
	  }
	    } } }
	//cout << "AIF: " << mean(aifcand,2).t() << endl;
	    aif.setvoxelts(mean(aifcand,2),x,y,z);
	  }
	  }
	} } }
    }

  void Correct_magnitude(ColumnVector& mag, const ColumnVector& batd, const float T1, const float dt=0.0, const float fa=0.0) {
    //magnitude correction to take acocutn fop differences in BAT between aif and tissue

    mag = SP(mag,exp( batd/T1 ));
    if (fa>0) {
      for (int v=1; v<=mag.Nrows(); v++) {
	mag(v) *= 1/pow( cos(fa/180*M_PI),floor((batd(v)-1e-3)/dt) ); //the 1e-3 deals with the case where batd is a integer mul;tiple of dt
      }
    }
    cout << endl;
  }

  void Estimate_BAT_difference(const Matrix& resid, ColumnVector& batd, const float dt) {
    //Estimate the BAT difference between AIF and tissue curve from residue function

    batd.ReSize(resid.Ncols());

    float magtemp;
    int battemp;
    for (int v=1; v<=resid.Ncols(); v++) {
      magtemp = (resid.Column(v)).Maximum1(battemp);
      if (battemp > resid.Nrows()/2) {
	battemp -= resid.Nrows();
      }
      batd(v) = dt*(battemp-1);
    }

  }

  void Estimate_onset(const Matrix& curves, ColumnVector& bate, const float dt) {
    //Estimate the time of onset of the curve using edge detection

    int ntpts=curves.Nrows();
    bate.ReSize(curves.Ncols());

    ColumnVector kernel(ntpts);
    kernel = 0.0;
    ColumnVector kern(7);
    kern << 0.006 << 0.061 << 0.242 << 0.383 << 0.242 << 0.061 << 0.006;
    kernel.Rows(1,kern.Nrows()) = kern;
    
    Matrix kernmtx;
    kernmtx = convmtx(kernel);
    
    ColumnVector smoothdata(curves.Column(1));
    ColumnVector dgrad(ntpts-1);
    for (int v=1; v<=curves.Ncols(); v++) {
      //apply smoothing kernel
      smoothdata = kernmtx*curves.Column(v);
      // take gradient (forward difference)
      dgrad = smoothdata.Rows(2,ntpts) - smoothdata.Rows(1,ntpts-1);
      float gthresh = 0.2*dgrad.Maximum();
      int i=1; bool cont=true;
      while ((i<ntpts) & cont) {
	if (dgrad(i)>gthresh) cont=false;
	else i++;
	       }
      //cout << i << " " ;
      bate.Row(v) = i*dt;
    }
    cout << endl;
  }
}
