/*   asl_functions.cc various functions for the manipulation of ASL data

      Michael Chappell - FMIRB Image Analysis Group

      Copyright (C) 2009 University of Oxford */

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

#include "asl_functions.h"

namespace OXASL {

  void data2stdform(Matrix& datamtx, vector<Matrix>& asldata, int ntis, bool isblocked, bool ispairs) {
    int nvox=datamtx.Ncols();
    int nmeas=datamtx.Nrows()/ntis;
    int nrpts;
    if (ispairs) nrpts=nmeas/2;
    else nrpts=nmeas;

    //cout << nmeas << " " << nrpts << " " << endl;

    if (isblocked)
      {
	Matrix thisti(nmeas,nvox);
	for (int ti=1; ti<=ntis; ti++) 
	  {
	    thisti=0;
	    //asldata[ti-1].ReSize(nvox,nmeas);
	    //extract the measurements for this TI
	    for (int i=1; i<=nrpts; i++)
	      {
		if (ispairs) {
		  thisti.Row(2*i-1) = datamtx.Row(2*(i-1)*ntis+2*ti-1);
		  thisti.Row(2*i) = datamtx.Row(2*(i-1)*ntis+2*ti);
		}
		else {
		  thisti.Row(i) = datamtx.Row((i-1)*ntis+ti);
		}
	      }
	  
	    asldata.push_back(thisti);
	  }

	/*
	//deal with orphan data at end (if present)
	if (datamtx.Nrows() > nmeas*ntis) {
	  int i=nrpts*2*ntis+1;
	  int ti=0;
	  while (i<= datamtx.Nrows()) {
	    if (ispairs) {
	      asldata[ti] = asldata[ti] & datamtx.Row(i);
	      i++;
	      if (i > datamtx.Nrows()) throw Exception("Cannot process this dataset as pairs: odd number of spare measurements at end");
	      asldata[ti] = asldata[ti] & datamtx.Row(i);
	    }
	    else {
	      asldata[ti] = asldata[ti] & datamtx.Row(i);
	    }
	    i++;
	    ti++;
	  }
	}
	*/
	
      }
    else 
      {

	for (int ti=1; ti<=ntis; ti++) 
	  //asldata[ti-1].ReSize(nvox,nmeas);
	  asldata.push_back(datamtx.Rows((ti-1)*nmeas+1,ti*nmeas));

	if (datamtx.Nrows() > nmeas*ntis) throw Exception("Orphaned data found at end of file - this is not logical when data is in TI blocks");
      }
  }

  void stdform2data(vector<Matrix>& asldata, Matrix& datareturn, bool outblocked, bool outpairs) {
    int ntis = asldata.size();
    int nvox = asldata[0].Ncols();
    int nmeas = asldata.back().Nrows(); //safer to determine this from the very last TI (in case of orphan measurements when nodiscard is turned on)
    int ninc=1;
    if (outpairs) ninc=2;
    
    if (outblocked) {
      datareturn.ReSize(ntis*nmeas,nvox);
      int idx=1;
      for (int m=1; m<=nmeas; m+=ninc)
	for (int ti=1; ti<=ntis; ti++) {
	  datareturn.Row(idx) = asldata[ti-1].Row(m);
	  idx++;
	  if (outpairs) {
	    datareturn.Row(idx) = asldata[ti-1].Row(m+1);
	    idx++;
	  }
	}
      assert(idx-1==ntis*nmeas);

      /*
      //deal with orphans
      if (asldata.front().Nrows() > nmeas) {
	int ti=0;
	int finalm = nmeas+ninc;
	while (asldata[ti].Nrows() > nmeas) {
	  datareturn.Row(idx) = asldata[ti].Row(finalm);
	  idx++;
	  if (outpairs) {
	    datareturn.Row(idx) = asldata[ti].Row(finalm+1);
	    idx++;
	  }
	  ti++;
	}
      }
      */
    }
    else {
      datareturn = asldata[0];
      for (int ti=1; ti<ntis; ti++) {
	datareturn &= asldata[ti];
      }
    }
  }

  void separatepairs(vector<Matrix>& asldata, vector<Matrix>& asldataodd, vector<Matrix>& asldataeven) {
    int ntis = asldata.size();
    int nmeas = asldata[0].Nrows();
    int nrpts=nmeas/2; //if we are using this function then the data must contain pairs

    // just in case the vectors are not empty to start with
    asldataodd.clear();
    asldataeven.clear();

    int idx;

    for (int ti=0; ti<ntis; ti++) {

	Matrix oddmtx;
	oddmtx = asldata[ti].Row(1);
	Matrix evenmtx;
	evenmtx = asldata[ti].Row(2);

	for (int r=2; r<=nrpts; r++) {
	  idx=(r-1)*2+1;
	  oddmtx &= asldata[ti].Row(idx);
	  evenmtx &= asldata[ti].Row(idx+1);
	}
	asldataodd.push_back(oddmtx);
	asldataeven.push_back(evenmtx);
      }
  }

  void mergepairs(vector<Matrix>& asldata, vector<Matrix>& asldataodd, vector<Matrix>&asldataeven) {
    int ntis = asldataodd.size();
    int nmeas = asldataodd[0].Nrows();
    int nrpts=nmeas; //asldataodd does not contain pairs

    asldata.clear(); //make sure this is clear

    for (int ti=0; ti<ntis; ti++) {
      Matrix aslmtx;
      aslmtx = asldataodd[ti].Row(1);
      aslmtx &= asldataeven[ti].Row(1);
      for (int r=2; r<=nrpts; r++) {
	aslmtx &= asldataodd[ti].Row(r);
	aslmtx &= asldataeven[ti].Row(r);
      }
      asldata.push_back(aslmtx);
    }
  }

  void timeans(vector<Matrix>& asldata, vector<Matrix>& meanreturn) {
    int ntis = asldata.size();
    int nvox = asldata[0].Ncols();
    Matrix meanti(ntis,nvox);
      for (int n=0; n<ntis; n++) 
	{
	  meanreturn.push_back(mean(asldata[n],1));
	}
  }

  void timeanout(vector<Matrix>& asldata, volume<float>& mask, string fname, bool outpairs) {
      cout << "Outputting ASL data mean at each TI" << endl;
      
	vector<Matrix> meanti;
      if (!outpairs) {
	timeans(asldata,meanti);
      }
      else {
	//need to preserve pairs in the data - separate
	vector<Matrix> asldataodd;
	vector<Matrix> asldataeven;
	separatepairs(asldata,asldataodd,asldataeven);

	//take mean of odds and evens
	vector<Matrix> meanodd;
	timeans(asldataodd,meanodd);
	vector<Matrix> meaneven;
	timeans(asldataeven,meaneven);

	//recombine
	mergepairs(meanti,meanodd,meaneven);
      }

      Matrix meantimat;
      stdform2data(meanti, meantimat, false, false);

      volume4D<float> meanout;
      meanout.setmatrix(meantimat,mask);
      save_volume4D(meanout,fname);
      
    }

  void splitout(vector<Matrix>& asldata, volume<float>& mask, string froot) {
      cout << "Splitting ASL data into files for each TI" << endl;
      int ntis = asldata.size();
      volume4D<float> blockout;
      for (int n=0; n<ntis; n++) 
	{
	  char cstr [5];
	  if (n<10) sprintf(cstr,"00%d",n);
	  else if (n<100) sprintf(cstr,"0%d",n);
	  else if (n<1000) sprintf(cstr,"%d",n);
	  else throw Exception("More than 1000 measurements in this ASL data file, sorry cannot handle this operation");
	  string tino(cstr);
	  
	  blockout.setmatrix(asldata[n],mask);
	  save_volume4D(blockout,froot+tino);
	}
    }

  void genepochs(vector<Matrix>& asldata, vector<Matrix>& epochreturn, int eplen, int epol) {
      
      int ntis = asldata.size();
      int nvox = asldata[0].Ncols();
      int nmeas = asldata[0].Nrows();

      epochreturn.clear();

      if (epol>=eplen) throw Exception("The epoch overlap may not exceed or equal the length of the epoch");
      int epadv = eplen-epol;

      int e=1;
      Matrix epoch_temp(ntis,nvox);
      while ((e-1)*epadv+eplen<=nmeas)
	{
	  //Matrix ti_temp(eplen,nvox);
	  //go through the TIs
	  for (int ti=1; ti<=ntis; ti++) {
	    epoch_temp.Row(ti) = mean(asldata[ti-1].Rows((e-1)*epadv+1,(e-1)*epadv+eplen),1);
	  }

	  epochreturn.push_back(epoch_temp);

	  e++;
	}

      e--; //correct for final e++
      int unused=nmeas-( (e-1)*epadv+eplen );
      if (unused>0) cout << "Number of measurements from end of data discarded: " << unused << endl;
    }

  void gentiepochs(Matrix& asldata, vector<Matrix>& epochreturn, int eplen, int epol) {
      // expects data in matrix form (blocks of repeats)

      int nvox = asldata.Ncols();
      int nmeas = asldata.Nrows();

      epochreturn.clear();

      if (epol>=eplen) throw Exception("The epoch overlap may not exceed or equal the length of the epoch");
      int epadv = eplen-epol;

      int e=1;
      Matrix epoch_temp(eplen,nvox);
      while ((e-1)*epadv+eplen<=nmeas)
	{
	  epoch_temp = asldata.Rows( (e-1)*epadv+1,(e-1)*epadv+eplen );
	  epochreturn.push_back(epoch_temp);
	  e++;
	}

      e--;
      
      int unused=nmeas-( (e-1)*epadv+eplen );
      if (unused>0) cout << "Number of measurements from end of data discarded: " << unused << endl;
    }


    void epochout(vector<Matrix>& asldata, volume<float>& mask, string froot, int eplen, int epol, bool outpairs, bool tiunit) {
      //Save out epochs of data.
      // tiunit indicates that we want epochs to be done over the Tis rather than over the repeats
    cout << "Ouput ASL data epochs" << endl;
 
    //generate the epochs
    vector<Matrix> theepochs;
    Matrix alldata;
    if (!outpairs) {
      if (!tiunit) {
	genepochs(asldata, theepochs, eplen, epol);
      }
      else {
	// epochs over the TIs need to convert the data to matrix form
	stdform2data(asldata,alldata,true,false);
	gentiepochs(alldata,theepochs,eplen,epol);
      }
    }
    else {
      //need to preserve pairs in the data - separate
      vector<Matrix> asldataodd;
      vector<Matrix> asldataeven;
      separatepairs(asldata,asldataodd,asldataeven);

      vector<Matrix> oddepochs;
      if (!tiunit) {
	genepochs(asldataodd, oddepochs, eplen, epol);
      }
      else {
	stdform2data(asldataodd,alldata,true,false);
	gentiepochs(alldata, oddepochs, eplen, epol);
      }

      vector<Matrix> evenepochs;
      if (!tiunit) {
	genepochs(asldataeven, evenepochs, eplen, epol);
      }
      else {
	stdform2data(asldataodd,alldata,true,false);
	gentiepochs(alldata, evenepochs, eplen, epol);
	
      }

      mergepairs(theepochs,oddepochs,evenepochs);
    }

    int nepoch = theepochs.size();
    volume4D<float> epochout;
    for (int e=0; e<nepoch; e++) {
      char cstr [5];
	  if (e<10) sprintf(cstr,"00%d",e);
	  else if (e<100) sprintf(cstr,"0%d",e);
	  else if (e<1000) sprintf(cstr,"%d",e);
	  else throw Exception("More than 1000 epochs in this ASL data file, sorry cannot handle this operation");
	  string epno(cstr);

	  epochout.setmatrix(theepochs[e],mask);
	  save_volume4D(epochout,froot+epno);
    }
  }

  ReturnMatrix SVDdeconv(const Matrix& data, const Matrix& aif) {
    // do a singular value deconvolution of the data to get residue function

    int nti = data.Nrows();
    int nvox = data.Ncols();
    float truncfac = 0.2;

    // voxelwise SVD deconvolution
    Matrix aifconv; Matrix residue(nti,nvox);
    DiagonalMatrix S; DiagonalMatrix D; Matrix U; Matrix V;
    for (int v=1; v<=nvox; v++) {
      //make convolution matrix
      aifconv = convmtx(aif.Column(v));
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

  void deconvout(vector<Matrix>& asldata, volume<float>& mask, Matrix& aif, string fname) {
    //perform deconvolution and output the magnitude and residue function

    //take the mean in each TI (we dont want mulitple repeats here)
    vector<Matrix> meandata;
    timeans(asldata, meandata);
    Matrix data;
    stdform2data(meandata,data,false,false); //way to get mean result from standard form into a matrix that can now be processed

    //do the deconvolution
    Matrix resid = SVDdeconv(data,aif);
    // extract magntiude and residue separately
    int nvox = data.Ncols();
    ColumnVector mag(nvox);
    for (int v=1; v<=nvox; v++) {
      mag(v) = (resid.Column(v)).Maximum();
      resid.Column(v) /= mag(v);
    }

    //output 
    volume4D<float> residout;
    residout.setmatrix(resid,mask);
    save_volume4D(residout,fname+"_residuals");

    volume4D<float> magout;
    magout.setmatrix(mag.AsMatrix(1,nvox),mask);
    save_volume4D(magout,fname+"_magntiude");
  }



}
