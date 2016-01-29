/*  MELODIC - Multivariate exploratory linear optimized decomposition into 
              independent components
    
    melhlprfns.cc - misc functions

    Christian F. Beckmann, FMRIB Analysis Group
    
    Copyright (C) 1999-2013 University of Oxford */

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

#include "melhlprfns.h"
#include "libprob.h"
#include "miscmaths/miscprob.h"
#include "miscmaths/t2z.h"
#include "miscmaths/f2z.h"

namespace Melodic{

  void update_mask(volume<float>& mask, Matrix& Data)
  {
    Matrix DStDev=stdev(Data);
    volume4D<float> tmpMask, RawData;
    tmpMask.setmatrix(DStDev,mask);

    float tMmax;
    volume<float> tmpMask2;
    tmpMask2 = tmpMask[0];

    tMmax = tmpMask2.max();
    double st_mean = DStDev.Sum()/DStDev.Ncols();
    double st_std  = stdev(DStDev.t()).AsScalar();
      
	volume4D<float> newmask;
    newmask = binarise(tmpMask2,(float) max((float) st_mean-3*st_std,(float) 0.01*st_mean),tMmax);

	Matrix newmaskM,newData;
	newmaskM = newmask.matrix(mask);
	int N = Data.Nrows();
	
	if(int(newmaskM.Row(1).SumAbsoluteValue() + 0.3) < newmaskM.Ncols()){
		RawData.setmatrix(Data.Row(1),mask);
		newData = RawData.matrix(newmask[0]);
		for(int r=2; r <= N; r++){
			RawData.setmatrix(Data.Row(r),mask);
			newData &= RawData.matrix(newmask[0]);
		}
		Data = newData;
		mask = newmask[0];  
	}
  }

  void del_vols(volume4D<float>& in, int howmany)
  {
    for(int ctr=1; ctr<=howmany; ctr++){
	in.deletevolume(ctr);
    }    
  }

  Matrix calc_FFT(const Matrix& Mat, const bool logpwr)
  {
    Matrix res;
      for(int ctr=1; ctr <= Mat.Ncols(); ctr++){
	      ColumnVector tmpCol;  
	      tmpCol=Mat.Column(ctr);
	      ColumnVector FtmpCol_real;
	      ColumnVector FtmpCol_imag;
	      ColumnVector tmpPow;
	      if(tmpCol.Nrows()%2 != 0){
	        Matrix empty(1,1); empty=0;
	        tmpCol &= empty;
	      }
	      RealFFT(tmpCol,FtmpCol_real,FtmpCol_imag);
	      tmpPow = pow(FtmpCol_real,2)+pow(FtmpCol_imag,2);
	      tmpPow = tmpPow.Rows(2,tmpPow.Nrows());
	      if(logpwr) tmpPow = log(tmpPow);
	      if(res.Storage()==0){res= tmpPow;}else{res|=tmpPow;}
      }
      return res;
  }  //Matrix calc_FFT()

  Matrix smoothColumns(const Matrix& inp)
  {
    Matrix temp(inp);
    int ctr1 = temp.Nrows();
    Matrix temp2(temp);
    temp2=0;
    
    temp = temp.Row(4) & temp.Row(3) & temp.Row(2) & temp & temp.Row(ctr1-1) 
      & temp.Row(ctr1-2) &temp.Row(ctr1-3);
    
    double kern[] ={0.0045 , 0.055, 0.25, 0.4, 0.25, 0.055, 0.0045};
    double fac = 0.9090909;

    
    for(int cc=1;cc<=temp2.Ncols();cc++){
      for(int cr=1;cr<=temp2.Nrows();cr++){
	temp2(cr,cc) = fac*( kern[0] * temp(cr,cc) + kern[1] * temp(cr+1,cc) + 
			     kern[2] * temp(cr+2,cc) + kern[3] * temp(cr+3,cc) + 
			     kern[4] * temp(cr+4,cc) + kern[5] * temp(cr+5,cc) + 
			     kern[6] * temp(cr+6,cc));
      }
    }
    return temp2;
  }  //Matrix smoothColumns()

  Matrix convert_to_pbsc(Matrix& inp)
  {
    Matrix meanimg;
    meanimg = mean(inp);
    float eps = 0.00001;

    for(int ctr=1; ctr <= inp.Ncols(); ctr++){
      if(meanimg(1,ctr) < eps) 
	meanimg(1,ctr) = eps;
    }

    for(int ctr=1; ctr <= inp.Nrows(); ctr++){
      Matrix tmp;
      tmp << inp.Row(ctr);
      inp.Row(ctr) << 100 * SP((tmp - meanimg),pow(meanimg,-1));   
    }

    inp = remmean(inp);
    return meanimg;
  }  //void convert_to_pbsc   

  RowVector varnorm(Matrix& in, int dim, float level, int econ)
  {
    SymmetricMatrix Corr(cov_r(in,false,econ));
    RowVector out;
    out = varnorm(in,Corr,dim,level, econ);
    return out;
  }  //RowVector varnorm

  void varnorm(Matrix& in, const RowVector& vars)
  {
    for(int ctr=1; ctr <=in.Nrows();ctr++)
      in.Row(ctr) = SD(in.Row(ctr),vars);
  }
	
  RowVector varnorm(Matrix& in, SymmetricMatrix& Corr, int dim, float level, int econ)
  { 
	
    Matrix tmpE, white, dewhite;
    RowVector tmpD, tmpD2;

    std_pca(remmean(in,2), Corr, tmpE, tmpD, econ);
    calc_white(tmpE,tmpD, dim, white, dewhite);
    
    Matrix ws = white * in;
    for(int ctr1 = 1; ctr1<=ws.Ncols(); ctr1++)
      for(int ctr2 = 1; ctr2<=ws.Nrows(); ctr2++)
				if(std::abs(ws(ctr2,ctr1)) < level)
	  			ws(ctr2,ctr1)=0;
    tmpD = stdev(in - dewhite*ws);
    for(int ctr1 = 1; ctr1<=tmpD.Ncols(); ctr1++)
      if(tmpD(ctr1) < 0.01){
				tmpD(ctr1) = 1.0;
				in.Column(ctr1) = 0.0*in.Column(ctr1);
      }
	  varnorm(in,tmpD);

    return tmpD;
  }  //RowVector varnorm



  Matrix SP2(const Matrix& in, const Matrix& weights, int econ)
  {
    Matrix Res;
    Res = in;
    if(econ>0){
      ColumnVector tmp;
      for(int ctr=1; ctr <= in.Ncols(); ctr++){
	tmp = in.Column(ctr);
	tmp = tmp * weights(1,ctr);
	Res.Column(ctr) = tmp;
      }
    }
    else{
      Res = ones(in.Nrows(),1)*weights.Row(1);
      Res = SP(in,Res);
    }
    return Res;
  }  //Matrix SP2

  void SP3(Matrix& in, const Matrix& weights)
  {
	for(int ctr=1; ctr <= in.Nrows(); ctr++){
		in.Row(ctr) << SP(in.Row(ctr),weights.AsRow());
	} 
  }

  Matrix corrcoef(const Matrix& in1, const Matrix& in2){
		Matrix tmp = in1;
		tmp |= in2;
		Matrix out;
		out=MISCMATHS::corrcoef(tmp,0);
		return out.SubMatrix(1,in1.Ncols(),in1.Ncols()+1,out.Ncols());
  }

  Matrix corrcoef(const Matrix& in1, const Matrix& in2, const Matrix& part){
	Matrix tmp1 = in1, tmp2 = in2, out;	
	if(part.Storage()){
		tmp1 = tmp1 - part * pinv(part) * tmp1;
		tmp2 = tmp2 - part * pinv(part) * tmp2;
	}
	
	out = corrcoef(tmp1,tmp2);
	return out;
  }

  float calc_white(const Matrix& tmpE, const RowVector& tmpD, const RowVector& PercEV,  int dim, Matrix& param, Matrix& paramS, Matrix& white, Matrix& dewhite)
  {

//	tmpD2= tmpD | tmpPD.AsRow().Columns(tmpPE.Ncols()-param.Ncols()+1,tmpPE.Ncols());
//  cerr << tmpPD.AsRow().Columns(tmpPE.Ncols()-param.Ncols()+1,tmpPE.Ncols()) << endl;

//    

//	Matrix tmpPE;
//	tmpPE = SP(param,ones(param.Nrows(),1)*pow(stdev(param,1)*std::sqrt((float)param.Nrows()),-1));

//	RE |= tmpPE;
//	RowVector tmpD2;
//	tmpD2 = tmpD | stdev(param,1).AsRow()*std::sqrt((float)param.Nrows());
//	RD << abs(diag(tmpD2.t()));
//    RD << RD.SymSubMatrix(N-dim+1,N);

	Matrix RE;
    DiagonalMatrix RD;
    int N = tmpE.Ncols();
    dim = std::min(dim,N);

//	cerr << stdev(param) << endl;
    RE = tmpE.Columns(std::min(N-dim+1+param.Ncols(),N-2),N);
	RE |= param;

//	cerr << paramS.Nrows() << " x " << paramS.Ncols() << endl;
//	cerr << paramS << endl;
	RowVector tmpD2;
	tmpD2 = tmpD | pow(paramS,2).AsRow();
    RD << abs(diag(tmpD2.t()));

//	cerr << " " <<tmpD2.Ncols() << " " << N << " " << dim << endl;
    RD << RD.SymSubMatrix(N-dim+1+param.Ncols(),N+param.Ncols());    

    float res = 1.0;    
    white = sqrt(abs(RD.i()))*RE.t();
    dewhite = RE*sqrt(RD);

    if(dim < N)
      res = PercEV(dim);
    return res;
  }  //Matrix calc_white

  float calc_white(const Matrix& tmpE, const RowVector& tmpD, const RowVector& PercEV, int dim, Matrix& white, Matrix& dewhite)
  {
    Matrix RE;
    DiagonalMatrix RD;
    int N = tmpE.Ncols();
    dim = std::min(dim,N);
    RE = tmpE.Columns(N-dim+1,N);
    RD << abs(diag(tmpD.t()));
    RD << RD.SymSubMatrix(N-dim+1,N);    

    float res = 1.0;    
    white = sqrt(abs(RD.i()))*RE.t();
    dewhite = RE*sqrt(RD);

    if(dim < N)
      res = PercEV(dim);
    return res;
  }  //Matrix calc_white

  void calc_white(const Matrix& tmpE, const RowVector& tmpD, int dim, Matrix& white, Matrix& dewhite)
  {
    RowVector tmp(tmpE.Ncols());
    float tmp2;
    tmp2 = calc_white(tmpE,tmpD, tmp, dim, white, dewhite); 
  }  //Matrix calc_white

  void calc_white(const Matrix& tmpE, const RowVector& tmpD, int dim, Matrix& param, Matrix& paramS, Matrix& white, Matrix& dewhite)
  {
    RowVector tmp(tmpE.Ncols());
    float tmp2;
    tmp2 = calc_white(tmpE,tmpD, tmp, dim, param, paramS, white, dewhite); 
  }  //Matrix calc_white

  void calc_white(const SymmetricMatrix& Corr, int dim, Matrix& white, Matrix& dewhite)
  {
    Matrix RE;
    DiagonalMatrix RD;
    RowVector tmp2;
    EigenValues(Corr,RD,RE);
    tmp2 = diag(RD).t();
    calc_white(RE,tmp2, dim, white, dewhite); 
  }  //Matrix calc_white
  
 
  void std_pca(const Matrix& Mat, const Matrix& weights, SymmetricMatrix& Corr, Matrix& evecs, RowVector& evals, int econ)
  {
    if(weights.Storage()>0)
      Corr = cov_r(Mat, weights, econ);
    else
      Corr = cov_r(Mat,false,econ);

    DiagonalMatrix tmpD;
    EigenValues(Corr,tmpD,evecs);
    evals = tmpD.AsRow();
  }  //void std_pca

  void std_pca(const Matrix& Mat, SymmetricMatrix& Corr, Matrix& evecs, RowVector& evals, int econ)
  {
    Matrix weights;
    std_pca(Mat,weights,Corr,evecs,evals, econ);
  }  //void std_pca

  void em_pca(const Matrix& Mat, Matrix& evecs, RowVector& evals, int num_pc, int iter)
  {
    Matrix guess;
    guess = normrnd(Mat.Nrows(),num_pc);
    em_pca(Mat,guess,evecs,evals,num_pc,iter);
  }  //void em_pca

  void em_pca(const Matrix& Mat, Matrix& guess, Matrix& evecs, RowVector& evals, int num_pc, int iter)
  {
    Matrix C;
    if(guess.Ncols() < num_pc){
      C=normrnd(Mat.Nrows(),num_pc);
      C.Columns(1,guess.Ncols()) = guess;
    }
    else
      C = guess;

    Matrix tmp, tmp2;
    for(int ctr=1; ctr <= iter; ctr++){
      // E-Step
      tmp = C.t()*C;
      tmp = tmp.i();
      tmp = tmp * C.t();
      tmp = tmp * Mat;
      // M-Step
      tmp2 = tmp * tmp.t();
      tmp2 = tmp2.i();
      tmp2 = Mat*tmp.t()*tmp2;
      C = tmp2;
    }
    
    symm_orth(C);
    Matrix Evc;
    SymmetricMatrix tmpC;
    RowVector Evl;
    tmp = C.t() * Mat;
    std_pca(tmp,tmpC,Evc,Evl);
    evals = Evl;
    evecs = C * Evc;
  }  //void em_pca

  float rankapprox(const Matrix& Mat, Matrix& cols, Matrix& rows, int dim)
  { 
    SymmetricMatrix Corr;
    Matrix Evecs, tmpWM, tmpDWM, tmp;
    RowVector Evals;
    std_pca(Mat.t(), Corr, Evecs, Evals);
    calc_white(Corr, dim, tmpWM, tmpDWM);
    tmp = tmpWM * Mat.t();
    cols = tmp.t();
    rows << tmpDWM;
		float res;
		Evals=fliplr(Evals);
		res = sum(Evals.Columns(1,dim),2).AsScalar()/sum(Evals,2).AsScalar()*100;
		return res;
  } // rankapprox

  RowVector krfact(const Matrix& Mat, Matrix& cols, Matrix& rows)
  {
		Matrix out; RowVector res(Mat.Ncols());
    for(int ctr1 = 1; ctr1 <= Mat.Ncols(); ctr1++){
			Matrix tmpVals(cols.Nrows(),rows.Nrows());
			for(int ctr2 = 1; ctr2 <= rows.Nrows(); ctr2++)
	  		tmpVals.Column(ctr2) << Mat.SubMatrix(cols.Nrows() * 
				(ctr2 - 1) + 1,cols.Nrows()*ctr2 ,ctr1,ctr1);
	
			Matrix tmpcols, tmprows;
	 		res(ctr1) =rankapprox(tmpVals, tmpcols, tmprows);
			cols.Column(ctr1) = tmpcols;
			rows.Column(ctr1) = tmprows;
    }
		return res;
  } // krfact

  RowVector krfact(const Matrix& Mat, int colnum, Matrix& cols, Matrix& rows)
  {
		RowVector res;
    cols = zeros(colnum,Mat.Ncols());
    rows = zeros(int(Mat.Nrows() / colnum),Mat.Ncols());
    res = krfact(Mat,cols,rows);
		return res;
  } // krfact

  Matrix krprod(const Matrix& cols, const Matrix& rows)
  {
    Matrix out;
    out = zeros(cols.Nrows()*rows.Nrows(),cols.Ncols());
    for(int ctr1 = 1; ctr1 <= cols.Ncols(); ctr1++)
      for(int ctr2 = 1; ctr2 <= rows.Nrows(); ctr2++)
	{
	  out.SubMatrix(cols.Nrows()*(ctr2-1)+1,cols.Nrows()*ctr2,ctr1,ctr1) << cols.Column(ctr1) * rows(ctr2,ctr1);
	}
    return out;
  } // krprod

  Matrix krapprox(const Matrix& Mat, int size_cols, int dim)
  {
    Matrix out, cols, rows;
    out = zeros(Mat.Nrows(), Mat.Ncols());
    cols = zeros(size_cols,Mat.Ncols());
    rows = zeros(int(Mat.Nrows() / size_cols), Mat.Ncols());
    krfact(Mat,cols,rows);
    out = krprod(cols, rows);
    return out;
  } // krapprox

  void adj_eigspec(const RowVector& in, RowVector& out1, RowVector& out2, RowVector& out3, int& out4, int num_vox, float resels)
  {
    RowVector AdjEV;
    AdjEV << in.AsRow();
    AdjEV = AdjEV.Columns(3,AdjEV.Ncols());
    AdjEV = AdjEV.Reverse();

    RowVector CircleLaw;
    int NumVox = (int) floor(num_vox/(2.5*resels));

    CircleLaw = Feta(int(AdjEV.Ncols()), NumVox);

    for(int ctr=1;ctr<=CircleLaw.Ncols(); ctr++){
      if(CircleLaw(ctr)<5*10e-10){CircleLaw(ctr) = 5*10e-10;}
    } 

    /*    float slope;
    slope = CircleLaw.Columns(int(AdjEV.Ncols()/4),AdjEV.Ncols() - 
			      int(AdjEV.Ncols()/4)).Sum() /  
      AdjEV.Columns(int(AdjEV.Ncols()/4),AdjEV.Ncols() - 
      int(AdjEV.Ncols()/4)).Sum();*/

    RowVector PercEV(AdjEV);
    PercEV = cumsum(AdjEV / sum(AdjEV,2).AsScalar());

    AdjEV << SP(AdjEV,pow(CircleLaw.Columns(1,AdjEV.Ncols()),-1));

    SortDescending(AdjEV);
    int maxEV = 1;
    float threshold = 0.98;
    for(int ctr_i = 1; ctr_i < PercEV.Ncols(); ctr_i++){ 
      if((PercEV(ctr_i)<threshold)&&(PercEV(ctr_i+1)>=threshold)){maxEV=ctr_i;}
    }

    if(maxEV<3){maxEV=PercEV.Ncols()/2;}
    RowVector NewEV;
    Matrix temp1;
    temp1 = abs(AdjEV);
    NewEV << temp1.SubMatrix(1,1,1,maxEV);

    AdjEV = (AdjEV - min(AdjEV).AsScalar())/(max(AdjEV).AsScalar() - min(AdjEV).AsScalar());

    out1 = AdjEV;
    out2 = PercEV;
    out3 = NewEV;
    out4 = maxEV;
  }  //adj_eigspec

 void adj_eigspec(const RowVector& in, RowVector& out1, RowVector& out2)
  {
    RowVector AdjEV, PercEV;
    AdjEV = in.Reverse();
    SortDescending(AdjEV);
  
    PercEV = cumsum(AdjEV / sum(AdjEV,2).AsScalar());
    AdjEV = (AdjEV - min(AdjEV).AsScalar())/(max(AdjEV).AsScalar() - min(AdjEV).AsScalar());
    out1 = AdjEV;
    out2 = PercEV;
  }  //adj_eigspec

  RowVector Feta(int n1, int n2)
  {
    float nu = (float) n1/n2; 
    float bm = pow((1-sqrt(nu)),2.0);
    float bp = pow((1+sqrt(nu)),2.0);

    float lrange = 0.9*bm;
    float urange = 1.1*bp;

    RowVector eta(30*n1);
    float rangestepsize = (urange - lrange) / eta.Ncols(); 
    for(int ctr_i = 1; ctr_i <= eta.Ncols(); ctr_i++){ 
      eta(ctr_i) = lrange + rangestepsize * (ctr_i);
    }

    RowVector teta(10*n1);
    teta = 0;
    float stepsize = (bp - bm) / teta.Ncols();
    for(int ctr_i = 1; ctr_i <= teta.Ncols(); ctr_i++){ 
      teta(ctr_i) = stepsize*(ctr_i);
    }  
    
    RowVector feta(teta);
    feta = SP(pow(2*M_PI*nu*(teta + bm),-1), pow(SP(teta, bp-bm-teta),0.5));
   
    teta = teta + bm;

    RowVector claw(eta);
    claw = 0;
    for(int ctr_i = 1; ctr_i <= eta.Ncols(); ctr_i++){
      double tmpval = 0.0;
      for(int ctr_j = 1; ctr_j <= teta.Ncols(); ctr_j++){
	if(( double(teta(ctr_j))/double(eta(ctr_i)) )<1)
	  tmpval += feta(ctr_j);
      }
      claw(ctr_i) = n1*(1-stepsize*tmpval);
    }
    
    RowVector Res(n1); //invert the CDF
    Res = 0;
    for(int ctr_i = 1; ctr_i < eta.Ncols(); ctr_i++){ //Should this be <= instead of <?
      if(floor(claw(ctr_i))>floor(claw(ctr_i+1))){
	Res(int(floor(claw(ctr_i)))) = eta(ctr_i);
      }
    }
 
    return Res;
  }  //RowVector Feta

  RowVector cumsum(const RowVector& Inp)
  {
    UpperTriangularMatrix UT(Inp.Ncols());
    UT=1.0;
    RowVector Res;
    Res = Inp * UT;
    return Res;
  }  //RowVector cumsum

  int ppca_dim(const Matrix& in, const Matrix& weights, Matrix& PPCA, RowVector& AdjEV, RowVector& PercEV,  SymmetricMatrix& Corr, Matrix& tmpE, RowVector &tmpD, float resels, string which)
  {   
    std_pca(in,weights,Corr,tmpE,tmpD);

    int maxEV = 1;
    RowVector NewEV;
    adj_eigspec(tmpD.AsRow(),AdjEV,PercEV,NewEV,maxEV,in.Ncols(),resels);
    
    int res;
		PPCA = ppca_est(NewEV, in.Ncols(),resels);
    ColumnVector tmp = ppca_select(PPCA, res, maxEV, which);
		
		PPCA = tmp | PPCA;
    return res;
  }  //int ppca_dim

  int ppca_dim(const Matrix& in, const Matrix& weights, Matrix& PPCA, RowVector& AdjEV, RowVector& PercEV, float resels, string which)
  {   
    RowVector tmpD;
    Matrix tmpE;
    SymmetricMatrix Corr;

    int res = ppca_dim(in, weights, PPCA, AdjEV, PercEV, Corr, tmpE, tmpD, resels, which);
    return res;
  }  //int ppca_dim

  int ppca_dim(const Matrix& in, const Matrix& weights, float resels, string which)
  {
    ColumnVector PPCA;
    RowVector AdjEV, PercEV;
    int res = ppca_dim(in,weights,PPCA,AdjEV,PercEV,resels,which);
    return res;
  }  //int ppca_dim

  ColumnVector ppca_select(Matrix& PPCAest, int& dim, int maxEV, string which)
  {
    RowVector estimators(5);
    estimators = 1.0;
    
    for(int ctr=1; ctr<=PPCAest.Ncols(); ctr++){
      PPCAest.Column(ctr) = (PPCAest.Column(ctr) - 
			   min(PPCAest.Column(ctr)).AsScalar()) / 
	( max(PPCAest.Column(ctr)).AsScalar() - 
	  min(PPCAest.Column(ctr)).AsScalar());
    }
    
    int ctr_i = 1;
    while((ctr_i< PPCAest.Nrows()-1)&&
	  (PPCAest(ctr_i,2) < PPCAest(ctr_i+1,2))&&(ctr_i<maxEV))
      {estimators(1)=ctr_i+1;ctr_i++;}
    ctr_i = 1;
    while((ctr_i< PPCAest.Nrows()-1)&&
	  (PPCAest(ctr_i,3) < PPCAest(ctr_i+1,3))&&(ctr_i<maxEV))
      {estimators(2)=ctr_i+1;ctr_i++;}
    ctr_i = 1;
    while((ctr_i< PPCAest.Nrows()-1)&&
	  (PPCAest(ctr_i,4) < PPCAest(ctr_i+1,4))&&(ctr_i<maxEV))
      {estimators(3)=ctr_i+1;ctr_i++;}
    ctr_i = 1;
    while((ctr_i< PPCAest.Nrows()-1)&&
	  (PPCAest(ctr_i,5) < PPCAest(ctr_i+1,5))&&(ctr_i<maxEV))
      {estimators(4)=ctr_i+1;ctr_i++;}
    ctr_i = 1;
    while((ctr_i< PPCAest.Nrows()-1)&&
	  (PPCAest(ctr_i,6) < PPCAest(ctr_i+1,6))&&(ctr_i<maxEV))
      {estimators(5)=ctr_i+1;ctr_i++;}

    int res = 0;
    ColumnVector PPCA;
 		RowVector PercEV(PPCAest.Column(1).t());
	  PercEV = cumsum(PercEV / sum(PercEV,2).AsScalar());

	  if(which == string("aut")) {
			if(int(estimators(2)) < int(estimators(1)) && 
				float(PercEV(int(estimators(2))))>0.8){
				res=int(estimators(2));
	      PPCA << PPCAest.Column(3);
			}else{
				res = int(estimators(1));
	      PPCA << PPCAest.Column(2);
			}
	  }			
    if(which == string("lap")){
      res = int(estimators(1));
      PPCA << PPCAest.Column(2);
    }
    if(which == string("bic")){
      res = int(estimators(2));
      PPCA << PPCAest.Column(3);
    }
    if(which == string("mdl")){
      res = int(estimators(3));
      PPCA << PPCAest.Column(4);
    }
    if(which == string("rrn")){
      res = int(estimators(4));
      PPCA << PPCAest.Column(5);
    }
    if(which == string("aic")){
      res = int(estimators(5));
      PPCA << PPCAest.Column(6);
    }
		if(which == string("median")){
			RowVector unsorted(estimators);	
			SortAscending(unsorted);
			ctr_i=1;
			res=int(unsorted(3));
			while(res != int(estimators(ctr_i)))			
				ctr_i++;
			PPCA << PPCAest.Column(ctr_i);
		}
    if(res==0 || which == string("mean")){//median estimator
      PPCA = mean(PPCAest.Columns(2,6),2);
			res=int(mean(estimators,2).AsScalar());
    }

    dim = res;
    return PPCA;
  }  //RowVector ppca_select

  Matrix ppca_est(const RowVector& eigenvalues, const int N1, const float N2)
  { 
    Matrix Res;
    Res = ppca_est(eigenvalues, (int) floor(N1/(2.5*N2)));
    return Res;
  }  //Matrix ppca_est

  Matrix ppca_est(const RowVector& eigenvalues, const int N)
  {
    RowVector logLambda(eigenvalues);
    logLambda = log(eigenvalues);

    int d = logLambda.Ncols();

    RowVector k(d);
    for(int ctr_i = 1; ctr_i <=d; ctr_i++){
      k(ctr_i)=ctr_i;
    }
   
    RowVector m(d);
    m=d*k-0.5*SP(k,k+1); 

    RowVector loggam(d);
    loggam = 0.5*k.Reverse();
    for(int ctr_i = 1; ctr_i <=d; ctr_i++){
      loggam(ctr_i)=lgam(loggam(ctr_i));
    }
    loggam = cumsum(loggam); 

    RowVector l_probU(d);
    l_probU = -log(2)*k + loggam - cumsum(0.5*log(M_PI)*k.Reverse());

    RowVector tmp1;
    tmp1 = -cumsum(eigenvalues).Reverse()+sum(eigenvalues,2).AsScalar();
    tmp1(1) = 0.95*tmp1(2);
    tmp1=tmp1.Reverse();

    RowVector tmp2;
    tmp2 = -cumsum(logLambda).Reverse()+sum(logLambda,2).AsScalar();
    tmp2(1)=tmp2(2);
    tmp2=tmp2.Reverse();

    RowVector tmp3;
    tmp3 = d-k;
    tmp3(d) = 1.0;

    RowVector tmp4;
    tmp4 = SP(tmp1,pow(tmp3,-1));    
    for(int ctr_i = 1; ctr_i <=d; ctr_i++){
      if(tmp4(ctr_i)<0.01){tmp4(ctr_i)=0.01;}
      if(tmp3(ctr_i)<0.01){tmp3(ctr_i)=0.01;}
      if(tmp1(ctr_i)<0.01){tmp1(ctr_i)=0.01;}
    }

    RowVector l_nu;
    l_nu = SP(-N/2*(d-k),log(tmp4));
    l_nu(d) = 0;

    RowVector l_lam;
    l_lam = -(N/2)*cumsum(logLambda);

    RowVector l_lhood;
    l_lhood = SP(pow(tmp3,-1),tmp2)-log(SP(pow(tmp3,-1),tmp1));

    Matrix t1,t2, t3;
    UpperTriangularMatrix triu(d);
    triu = 1.0;
    for(int ctr_i = 1; ctr_i <= triu.Ncols(); ctr_i++){
      triu(ctr_i,ctr_i)=0.0;
    }
    t1 = (ones(d,1) * eigenvalues);
    t1 = SP(triu,t1.t() - t1);
    t2 = pow(tmp4,-1).t()*ones(1,d);
    t3 = ones(d,1)*pow(eigenvalues,-1);
    t2 = SP(triu, t2.t()-t3.t());
    for(int ctr_i = 1; ctr_i <= t1.Ncols(); ctr_i++){
      for(int ctr_j = 1; ctr_j <= t1.Nrows(); ctr_j++){
	if(t1(ctr_j,ctr_i)<=0){t1(ctr_j,ctr_i)=1;}
      } 
    }
    for(int ctr_i = 1; ctr_i <= t2.Ncols(); ctr_i++){
      for(int ctr_j = 1; ctr_j <= t2.Nrows(); ctr_j++){
	if(t2(ctr_j,ctr_i)<=0){t2(ctr_j,ctr_i)=1;}
      }
    } 
    t1 = cumsum(sum(log(t1),2).AsRow());
    t2 = cumsum(sum(log(t2),2).AsRow());

    RowVector l_Az(d);
    l_Az << (t1+t2);

    RowVector l_lap;
    l_lap = l_probU + l_nu +l_Az + l_lam + 0.5*log(2*M_PI)*(m+k)-0.5*log(N)*k;
 
    RowVector l_BIC;
    l_BIC = l_lam + l_nu - 0.5*log(N)*(m+k);

    RowVector l_RRN;
    l_RRN = -0.5*N*SP(k,log(SP(cumsum(eigenvalues),pow(k,-1))))+l_nu;

    RowVector l_AIC;
    l_AIC = -2*N*SP(tmp3,l_lhood)+ 2*(1+d*k+0.5*(k-1));
    l_AIC = -l_AIC;

    RowVector l_MDL;
    l_MDL = -N*SP(tmp3,l_lhood)+ 0.5*(1+d*k+0.5*(k-1))*log(N);
    l_MDL = -l_MDL;

    Matrix Res;

    Res = eigenvalues.t();
    Res |= l_lap.t();
    Res |= l_BIC.t();
    Res |= l_MDL.t();
    Res |= l_RRN.t();
    Res |= l_AIC.t();
    
   
    return Res;
  }  //Matrix ppca_est

  ColumnVector acf(const ColumnVector& in, int order)
  {
    double meanval;
    meanval = mean(in).AsScalar();
    int tpoints = in.Nrows();
    
    ColumnVector y, res;
    Matrix X, tmp;

    y = in.Rows(order+1,tpoints) - meanval;
    X = zeros(order + 1, order);
    for(int ctr1 = 1; ctr1 <= order; ctr1++)
      X.Column(ctr1) = in.Rows(order + 1 - ctr1, tpoints - ctr1) - meanval;
    tmp = X.t()*X;
    tmp = tmp.i();
    tmp = tmp * X.t();
    res << tmp * y;
    return res;
  }  //ColumnVector acf

  ColumnVector pacf(const ColumnVector& in, int maxorder)
  {
    int tpoint = in.Nrows();
    ColumnVector res;
    res = acf(in, maxorder);
    for(int ctr1 = 1; ctr1 <= maxorder; ctr1++)
      if ( res.Column(ctr1).AsScalar() <  (1/tpoint) + 2/(float)std::pow(tpoint,0.5)) 
	res.Column(ctr1) = 0;
    return res;
  }  //ColumnVector pacf
  
  Matrix est_ar(const Matrix& Mat, int maxorder)
  {
    Matrix res;
    res = zeros(maxorder, Mat.Ncols());
    ColumnVector tmp;
    for (int ctr = 1; ctr <= Mat.Ncols(); ctr++){
      tmp = pacf(Mat.Column(ctr));
      res.Column(ctr) = tmp;
    }
    return res;
  }  //Matrix est_ar

  ColumnVector gen_ar(const ColumnVector& in, int maxorder)
  {
    float sdev;
    sdev = stdev(in).AsScalar();
    ColumnVector x, arcoeff, scaled;
    scaled = in / sdev;
    arcoeff = pacf( scaled, maxorder);
    x = normrnd(in.Nrows(),1).AsColumn() * sdev;
    for(int ctr1=2; ctr1 <= in.Nrows(); ctr1++)
      for(int ctr2 = 1; ctr2 <= maxorder; ctr2++)
	x(ctr1) = arcoeff(ctr2) * x(std::max(1,int(ctr1-ctr2))) + x(ctr1);
    return x;
  }  //ColumnVector gen_ar

  Matrix gen_ar(const Matrix& in, int maxorder)
  {
    Matrix res;
    res = in;
    ColumnVector tmp;
    for(int ctr=1; ctr <= in.Ncols(); ctr++){
      tmp = in.Column(ctr);
      res.Column(ctr) = gen_ar(tmp, maxorder);
    } 
    return res;
  }  //Matrix gen_ar

  Matrix gen_arCorr(const Matrix& in, int maxorder)
  {
    Matrix res;
    res = zeros(in.Nrows(), in.Nrows());
    ColumnVector tmp;
    for(int ctr=1; ctr<= in.Ncols(); ctr++){
      tmp = in.Column(ctr);
      tmp = gen_ar(tmp, maxorder);
      res += tmp * tmp.t();
    }
    return res;
  }  //Matrix gen_arCorr

	void basicGLM::olsfit(const Matrix& data, const Matrix& design, 
		const Matrix& contrasts, int requestedDOF)
	{
		beta = zeros(design.Ncols(),1); 
		residu = zeros(1); sigsq = -1.0*ones(1); varcb = -1.0*ones(1); 
		t = zeros(1); z = zeros(1); p=-1.0*ones(1);
		dof = (int)-1; cbeta = -1.0*ones(1); 

		if(data.Nrows()==design.Nrows()){
			Matrix dat = data;
			Matrix tmp = design.t()*design;
			Matrix pinvdes = tmp.i()*design.t();
			
			beta = pinvdes * dat;
			residu = dat - design*beta;
			dof = ols_dof(design);
			if ( requestedDOF>0)
			  dof = requestedDOF;
			sigsq = sum(SP(residu,residu))/dof;
			
			float fact = float(dof) / design.Ncols();
			f_fmf =  SP(sum(SP(design*beta,design*beta)),pow(sum(SP(residu,residu)),-1)) * fact;
		
			pf_fmf = f_fmf.Row(1); 
			for(int ctr1=1;ctr1<=f_fmf.Ncols();ctr1++)
				pf_fmf(1,ctr1) = 1.0-MISCMATHS::fdtr(design.Ncols(),dof,f_fmf.Column(ctr1).AsScalar());
				
			if(contrasts.Storage()>0 && contrasts.Ncols()==beta.Nrows()){
				cbeta = contrasts*beta;
				Matrix tmp = contrasts*pinvdes*pinvdes.t()*contrasts.t();
				varcb = diag(tmp)*sigsq;
				t = SP(cbeta,pow(varcb,-0.5));
				z = t; p=t; 
				for(int ctr1=1;ctr1<=t.Ncols();ctr1++){
					ColumnVector tmp = t.Column(ctr1);
					T2z::ComputeZStats(varcb.Column(ctr1),cbeta.Column(ctr1),dof, tmp);
					z.Column(ctr1) << tmp;
					T2z::ComputePs(varcb.Column(ctr1),cbeta.Column(ctr1),dof, tmp);
					p.Column(ctr1) << exp(tmp);
				}
			}
		}	
	
	}


}
