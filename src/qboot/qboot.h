/*  Copyright (C) 2010 University of Oxford  

    Stamatios Sotiropoulos - FMRIB Image Analysis Group  */

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

#if !defined(qboot_h)
#define qboot_h


#include <iostream>
#include "stdlib.h"
#include "libprob.h"
#include <cmath>
#include "miscmaths/miscprob.h"
#include "miscmaths/miscmaths.h"
#include "newimage/newimageall.h"
#include "qbootOptions.h"
#include "peakfinder.h"
#include "sphere_tessellation.h"


using namespace std; 
using namespace NEWMAT;
using namespace NEWIMAGE;
using namespace MISCMATHS;

namespace ODFs{
  #define b0max 50   //b values up to b0max will be considered b0s
  #define shellrange 50  //how much +/- b value variability in each q shell is allowed. E.g. a shell at b=1000, can range from b=950 to b=1050. 

  
  /////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////
  //Class to estimate the ODF SH coefficients in a single voxel using single-shell data
  class ODF_voxel{
    Matrix m_coeff_samples;       //Stores the ODF SH coefficients (N_coeff x N_samples)
    
    ColumnVector m_coeff;         //Vector that holds the current ODF SH coefficients
    ColumnVector m_coeff_signal;  //Vector that holds the signal SH coefficients
    ColumnVector m_current_data;  //Vectors used during bootstrapping to keep the current data state
    ColumnVector m_current_residuals; //and the current residuals
    ColumnVector m_residuals;     //Vector that holds the original residuals that will be bootstrapped
    
    const ColumnVector& m_data;   //Vector with signal attenuations (for CSA-ODFs, transform to log of the -log of the signal attenuations)
    const int& m_modelnum;        //1 for Descoteaux's ODFs, 2 for CSA-ODFs
    const Matrix& m_SHT;          //Pseudoinverse of spherical harmonic transform matrix, used to get the signal SH coefficients 
    const Matrix& m_HAT;          //Hat matrix, used during bootstrap 
    const DiagonalMatrix& m_scaling_factors;//Holds the scaling factors to obtain the ODF SH coefficients from the signal SH ones
     __attribute__((unused)) const int& m_SH_order;        //maximum SH order (must be even)  
    const int& Num_of_coeff;      //Number of coefficients that will be estimated.
    const int& n_samples;         //Number of samples per coefficient
    const float& m_delta;         //Regularization parameter for signal intensities

  public:
    //constructor
  ODF_voxel(const ColumnVector& data, const Matrix& SHT, const Matrix& HAT, const DiagonalMatrix& scaling_factors, const int& modelnum, const int& SH_order, const int& Num_coeff, const int& nsamples, const float& delta): m_data(data), m_modelnum(modelnum), m_SHT(SHT), m_HAT(HAT), m_scaling_factors(scaling_factors), m_SH_order(SH_order),Num_of_coeff(Num_coeff),n_samples(nsamples), m_delta(delta) {
      m_current_data=m_data;
      m_coeff_samples.ReSize(Num_coeff,nsamples);
      m_coeff_samples=0;
    }
    
    //destructor
    ~ODF_voxel(){}
    
    Matrix& get_SHcoeff_ref() { return m_coeff_samples;}
    Matrix get_SHcoeff() { return m_coeff_samples;}


    void compute_signal_SH_coeff(){ 
      m_coeff_signal=m_SHT*m_current_data; 
    }
    

    //Regularizes signal attenuations (provided by Iman Aganj)
    void regularize_intensities (ColumnVector& attenuation, const float& delta) const{
      if (delta==0){   //Clip attenuation values. If att<0 => att=0, if att>1 => att=1 
	for (int i=1; i<=attenuation.Nrows(); i++)
	  attenuation(i)=(attenuation(i)>=0 && attenuation(i)<=1)*attenuation(i)+(attenuation(i)>1);
      } 
      else{            //Use function from Aganj et al, MRM, 2010
	for (int i=1; i<=attenuation.Nrows(); i++)
	  attenuation(i)=(attenuation(i)<0)*(0.5*delta) + (attenuation(i)>=0 && attenuation(i)<delta)*(0.5*delta+0.5*(attenuation(i)*attenuation(i))/delta) 
	                                        + (attenuation(i)>=delta && attenuation(i)<1-delta)*attenuation(i)  
	                                        + (attenuation(i)>=1-delta && attenuation(i)<1)*(1-0.5*delta-0.5*((1-attenuation(i))*(1-attenuation(i)))/delta) 
	                                        + (attenuation(i)>=1)*(1-0.5*delta);
	}
    }


    //Take the log(-log()) of the signal attenuations
    void log_log(ColumnVector& signal) const{
      for (int i=1; i<=signal.Nrows(); i++)
	signal(i)=log(-log(signal(i)));
    }

    //Transform log intensities back to signal attenuations by taking the exp(-exp())
    void exp_exp(ColumnVector& signal) const{
      for (int i=1; i<=signal.Nrows(); i++)
	signal(i)=exp(-exp(signal(i)));
    }


    //Performs Linear Regression and estimates SH coefficients
    void compute_ODF_SH_coeff(){ 
      m_coeff_signal = m_SHT*m_current_data;
      m_coeff=m_scaling_factors*m_coeff_signal;
      if (m_modelnum==1){                    
	double norm_const=m_coeff(1)*2*sqrt(M_PI);  //Normalize coefficients
	for (int m=1; m<=Num_of_coeff; m++)
	  m_coeff(m)=m_coeff(m)/norm_const;  
      }
      else if (m_modelnum==2)
	m_coeff(1)=1/(2*sqrt(M_PI));
    }

    
    void bootstrapping() {
      if (m_modelnum==1)
	bootstrapping_modelnum1();
      else if (m_modelnum==2)
	bootstrapping_modelnum2();
    }

    //Perform wild residual bootstrap, according to (Whitcher et al, 2008), for Tuch's ODFs
    void bootstrapping_modelnum1(){
      ColumnVector predicted;
      
      regularize_intensities(m_current_data,0);//Regularize attenuations
      compute_ODF_SH_coeff();               //Fit the model to the original data
      m_coeff_samples.Column(1)=m_coeff;    //Save the estimates for the SH coefficients (deterministic estimates)
      predicted=m_HAT*m_current_data;       //Get the predicted data
      m_residuals=m_current_data-predicted; //And then the residuals
      modify_residuals();                   //Modify the residuals
      m_current_residuals=m_residuals;
   
      for (int n=2; n<=n_samples; n++){
	bootstrap_residuals();                //Get boostrapped residuals
	m_current_data=predicted+m_current_residuals; //Get the new data state
	regularize_intensities(m_current_data,0);//Regularize attenuations
	compute_ODF_SH_coeff();               //Fit the model to the current data state
	m_coeff_samples.Column(n)=m_coeff;    //Save the estimates for the SH coefficients
      }
    }


    //Perform wild residual bootstrap for CSA-ODFs
    void bootstrapping_modelnum2(){
      ColumnVector predicted;
      
      regularize_intensities(m_current_data,m_delta);//Regularize attenuations
      log_log(m_current_data);
      compute_ODF_SH_coeff();               //Fit the model to the original data
      m_coeff_samples.Column(1)=m_coeff;    //Save the estimates for the SH coefficients (deterministic estimates)
      predicted=m_HAT*m_current_data;       //Get the predicted data
      m_residuals=m_current_data-predicted; //And then the residuals
      modify_residuals();                   //Modify the residuals
      m_current_residuals=m_residuals;
  
      for (int n=2; n<=n_samples; n++){
	bootstrap_residuals();                //Get boostrapped residuals
	m_current_data=predicted+m_current_residuals; //Get the new data state
	exp_exp(m_current_data);
	regularize_intensities(m_current_data,m_delta);//Regularize attenuations
	log_log(m_current_data);
	compute_ODF_SH_coeff();               //Fit the model to the current data state
	m_coeff_samples.Column(n)=m_coeff;    //Save the estimates for the SH coefficients
      }
    }

    void modify_residuals() {
      for (int m=1; m<=m_residuals.Nrows(); m++){
	m_residuals(m)/=sqrt(1-m_HAT(m,m));   //Correct for heteroscedasticity
      }
    }
  

    void bootstrap_residuals() {
      for (int m=1; m<=m_residuals.Nrows(); m++){
	if (unifrnd().AsScalar()<0.5)         //Change the sign with probability 0.5
	  m_current_residuals(m)=-m_residuals(m); 
	else
	  m_current_residuals(m)=m_residuals(m); 
      }
    }

  };


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Class to estimate the ODF SH coefficients in a single voxel using the biexponential model and multi-shell data
  class ODF_voxel_biexp{
    Matrix m_coeff_samples;       //Stores the ODF SH coefficients (N_coeff x N_samples)
    
    ColumnVector m_coeff;         //Vector that holds the current ODF SH coefficients
    ColumnVector m_coeff_signal;  //Vector that holds the signal SH coefficients
    ColumnVector m_current_data;  //Vectors used during bootstrapping to keep the current data state
    ColumnVector m_current_residuals; //and the current residuals
    ColumnVector m_residuals;     //Vector that holds the original residuals that will be bootstrapped
    
    const Matrix& m_data;         //NxM Matrix with signal attenuations. Each of the N rows corresponds to measurements along the same direction at M(=3) different b-values 
    const Matrix& m_SHT;          //Pseudoinverse of spherical harmonic transform matrix, used to get the signal SH coefficients 
    const Matrix& m_HAT;          //Hat matrix, used during bootstrap 
    const DiagonalMatrix& m_scaling_factors;//Holds the scaling factors to obtain the ODF SH coefficients from the signal SH ones
    __attribute__((unused)) const int& m_SH_order;        //maximum SH order (must be even)  
    __attribute__((unused)) const int& Num_of_coeff;      //Number of coefficients that will be estimated.
    const int& n_samples;         //Number of samples per coefficient
    const float& m_delta;         //Regularization parameter for signal intensities

  public:
    //constructor
  ODF_voxel_biexp(const Matrix& data, const Matrix& SHT, const Matrix& HAT, const DiagonalMatrix& scaling_factors, const int& SH_order, const int& Num_coeff, const int& nsamples, const float& delta): m_data(data), m_SHT(SHT), m_HAT(HAT), m_scaling_factors(scaling_factors), m_SH_order(SH_order),Num_of_coeff(Num_coeff),n_samples(nsamples), m_delta(delta) {
      m_current_data.ReSize(data.Nrows()); m_current_data=0;
      m_coeff_samples.ReSize(Num_coeff,nsamples);
      m_coeff_samples=0;
    }
    
    //destructor
    ~ODF_voxel_biexp(){}
    
    Matrix& get_SHcoeff_ref() { return m_coeff_samples;}
    Matrix get_SHcoeff() { return m_coeff_samples;}

    void compute_signal_SH_coeff(){ 
      m_coeff_signal=m_SHT*m_current_data; 
    }
    

    //Clip attenuation values. If att<0 => att=0, if att>1 => att=1 
    void clip_intensities (ColumnVector& attenuation) const{
      for (int i=1; i<=attenuation.Nrows(); i++)
	attenuation(i)=(attenuation(i)>=0 && attenuation(i)<=1)*attenuation(i)+(attenuation(i)>1);
    }


    //Regularizes signal attenuations (provided by Iman Aganj)
    void proj1(Matrix& attenuation, const float& delta) const{
      if (delta==0){   //Clip attenuation values. If att<0 => att=0, if att>1 => att=1 
	for (int i=1; i<=attenuation.Nrows(); i++)
	  for (int j=1; j<=attenuation.Ncols(); j++)
	    attenuation(i,j)=(attenuation(i,j)>=0 && attenuation(i,j)<=1)*attenuation(i,j)+(attenuation(i,j)>1);
      } 
      else{            //Use function from Aganj et al, MRM, 2010
	for (int i=1; i<=attenuation.Nrows(); i++)
	  for (int j=1; j<=attenuation.Ncols(); j++)
	    attenuation(i,j)=(attenuation(i,j)<0)*(0.5*delta) + (attenuation(i,j)>=0 && attenuation(i,j)<delta)*(0.5*delta+0.5*(attenuation(i,j)*attenuation(i,j))/delta) 
	      + (attenuation(i,j)>=delta && attenuation(i,j)<1-delta)*attenuation(i,j)+(attenuation(i,j)>=1-delta && attenuation(i,j)<1)*(1-0.5*delta-0.5*((1-attenuation(i,j))*(1-attenuation(i,j)))/delta) 
	      + (attenuation(i,j)>=1)*(1-0.5*delta);
	}
    }

    //Implements the first part of the projection onto the subspace defined by
    //inequalities [31] (Aganj, 2010). Taking the logarithm of the first three inequalities
    //converts the quadratic inequalities to linear ones.
    void proj2(Matrix& attenuation, const float& delta) const{
      float sF=sqrt(5);
      ColumnVector s0, a0, b0,ta,tb,e,m,a,b; 
      Matrix T0, E, C(attenuation.Nrows(),7), A(attenuation.Nrows(),7), B(attenuation.Nrows(),7);
      
      T0=-(MISCMATHS::log(attenuation));
      s0=MISCMATHS::sum(T0,2);
      a0=MISCMATHS::SD(T0.Column(1),s0);
      b0=MISCMATHS::SD(T0.Column(2),s0);
      ta=3*a0;   tb=3*b0;
      e=tb-2*ta; m=2*tb+ta;
       
      for (int j=1; j<=attenuation.Nrows(); j++){
	C(j,1)=(tb(j)<1+3*delta) && (0.5+1.5*(sF+1)*delta<ta(j)) && (ta(j)<1-3*(sF+2)*delta);
	C(j,2)=(e(j)<=-1+3*(2*sF+5)*delta) && (ta(j)>=1-3*(sF+2)*delta);
	C(j,3)=(m(j)>3-3*sF*delta) && (-1+3*(2*sF+5)*delta<e(j)) && (e(j)<-3*sF*delta);
	C(j,4)=(m(j)>=3-3*sF*delta) && (e(j)>=-3*sF*delta);
	C(j,5)=(2.5+1.5*(5+sF)*delta<m(j)) && (m(j)<3-3*sF*delta) && (e(j)>-3*sF*delta);
	C(j,6)=(ta(j)<=0.5+1.5*(sF+1)*delta) && (m(j)<=2.5+1.5*(5+sF)*delta);
	C(j,7)=!(C(j,1) || C(j,2) || C(j,3) || C(j,4) || C(j,5) || C(j,6));

	A(j,1)=C(j,1)*a0(j);
	A(j,2)=C(j,2)*(1.0/3.0-(sF+2)*delta);
	A(j,3)=C(j,3)*(0.2+0.8*a0(j)-0.4*b0(j)-delta/sF);
	A(j,4)=C(j,4)*(0.2+delta/sF);
	A(j,5)=C(j,5)*(0.2*a0(j)+0.4*b0(j)+2*delta/sF);
	A(j,6)=C(j,6)*(1.0/6.0+0.5*(sF+1)*delta);
        A(j,7)=C(j,7)*a0(j);
   
	B(j,1)=C(j,1)*(1.0/3.0+delta);
	B(j,2)=C(j,2)*(1.0/3.0+delta);
	B(j,3)=C(j,3)*(0.4-0.4*a0(j)+0.2*b0(j)-2*delta/sF);
	B(j,4)=C(j,4)*(0.4-3*delta/sF); 
	B(j,5)=C(j,5)*(0.4*a0(j)+0.8*b0(j)-delta/sF);
	B(j,6)=C(j,6)*(1.0/3.0+delta);
	B(j,7)=C(j,7)*b0(j);
      }
      a=MISCMATHS::sum(A,2); b=MISCMATHS::sum(B,2);
      E=(SP(a,s0) | SP(b,s0) | SP(1-a-b,s0));
      attenuation=exp(-E); 
    }

    
    // Implements the second part of the projection onto the subspace defined by 
    //inequalities [31]. (Aganj, 2010). 
    void proj3(ColumnVector& A, ColumnVector& a, ColumnVector& b, const float& delta){
      float del, s6 = sqrt(6), s15 = s6/2.0; 
      Matrix AM, aM, bM;
      ColumnVector d(A.Nrows()), I(A.Nrows());
      d=delta;
      
      AM = (A | A | A | d | ((A+a-b-s6*d)/3.0) | d | d | d | A | (0.2*(2*a+A-2*(s6+1)*d)) | (0.2*(-2*b+A+2-2*(s6+1)*d)) | d | d | d | (0.5-(1+s15)*d) );
      aM = (a | a | (1-d) | a | ((2*A+5*a+b+s6*d)/6.0) | a | (1-d) | (0.5*(a+b)+(1+s15)*d) | (1-d) | (0.2*(4*a+2*A+(s6+1)*d)) | (1-d) | ((s6+3)*d) | (1-d) | (1-d) | (1-d) );
      bM = (b | d | b | b | ((-2*A+a+5*b-s6*d)/6.0) | d | b | (0.5*(a+b)-(1+s15)*d) | d | d | (0.2*(4*b-2*A+1-(s6+1)*d)) | d | d | (1-(s6+3)*d) | d);
      
      Matrix R2(AM.Nrows(),AM.Ncols());
      
      del = delta*.99;
      for (int i=1; i<=AM.Nrows(); i++){
	for (int j=1; j<=AM.Ncols(); j++){
	  if (del<AM(i,j) && 2*(AM(i,j)+del*s15)<aM(i,j)-bM(i,j) && bM(i,j)>del && aM(i,j)<1-del)
	    R2(i,j) = (AM(i,j)-A(i))*(AM(i,j)-A(i))+ (aM(i,j)-a(i))*(aM(i,j)-a(i))+(bM(i,j)-b(i))*(bM(i,j)-b(i));
	  else
	    R2(i,j) = 1e20;
	}
	int ind;
	R2.Row(i).Minimum1(ind);
	I(i)=ind;
      }
      for (int i=1; i<=A.Nrows(); i++){
	A(i) = AM(i,(int)I(i));
	a(i) = aM(i,(int)I(i));
	b(i) = bM(i,(int)I(i));
      }
    } 


    void filter_and_param_estimation(){
      float f;
      Matrix filt_data; ColumnVector P2(m_data.Nrows()), A(m_data.Nrows()), B2(m_data.Nrows()), P(m_data.Nrows()), B(m_data.Nrows()),alpha, beta;
      filt_data=m_data;  //m_data assumed to be a Nx3 matrix, with N directions measured at three b-values

      proj1(filt_data,m_delta);
      proj2(filt_data,m_delta);
      for (int i=1; i<=m_data.Nrows(); i++){  //for each direction
	//Estimate alpha,beta,f
	P2(i)=filt_data(i,2)-filt_data(i,1)*filt_data(i,1);
	A(i)=(filt_data(i,3)-filt_data(i,1)*filt_data(i,2))/(2*P2(i));
	B(i)=A(i)*A(i)-(filt_data(i,1)*filt_data(i,3)-filt_data(i,2)*filt_data(i,2))/P2(i);
	if (B(i)<0)
	  B(i)=0;
	if (P2(i)<0)
	  P2(i)=0;
      }
      P2=sqrt(P2); B=sqrt(B);
      alpha=A+B; beta=A-B;
      proj3(P2, alpha, beta, m_delta);  //Change P2, alpha, beta

      for (int i=1; i<=m_data.Nrows(); i++){  //for each direction
	//fraction is computed in a bit different way from Eq. [30]. Since the next
	//line is correct with both plus and minus signs, we need to test both
	//answers to find out which one to use.
	float ER1, ER2, temp=(2*P2(i)/(alpha(i)-beta(i)));
	f=0.5+0.5*sqrt(1-temp*temp);
	ER1= fabs(f*(alpha(i)-beta(i))+beta(i)-filt_data(i,1))+fabs(f*(alpha(i)*alpha(i)-beta(i)*beta(i))+beta(i)*beta(i)-filt_data(i,2))+fabs(f*(alpha(i)*alpha(i)*alpha(i)-beta(i)*beta(i)*beta(i))+beta(i)*beta(i)*beta(i)-filt_data(i,3));
	ER2= fabs((1-f)*(alpha(i)-beta(i))+beta(i)-filt_data(i,1))+fabs((1-f)*(alpha(i)*alpha(i)-beta(i)*beta(i))+beta(i)*beta(i)-filt_data(i,2))+fabs((1-f)*(alpha(i)*alpha(i)*alpha(i)-beta(i)*beta(i)*beta(i))+beta(i)*beta(i)*beta(i)-filt_data(i,3));
	f = f*(ER1<ER2)+(1-f)*(ER1>=ER2);

	//Then obtain the expression for the bi-exponential model of the signal attenuation 
	m_current_data(i)=f*log(-log(alpha(i)))+(1-f)*log(-log(beta(i)));
      }
   } 


    //Performs Linear Regression and estimates SH coefficients
    void compute_ODF_SH_coeff(){ 
      m_coeff_signal = m_SHT*m_current_data;
      m_coeff=m_scaling_factors*m_coeff_signal;
      m_coeff(1)=1/(2*sqrt(M_PI));
    }

    
    //Perform wild residual bootstrap, according to (Whitcher et al, 2008), for multi-shell CSA ODFs.
    //The uncertainty is not given the data directly, but given a biexponential model of the data.
    //Bootstrapping is performed on Eq. [24] of (Aganj, MRM, 2010), i.e. residuals are added to the biexponential
    //model approximation, not the data!
    void bootstrapping(){
      ColumnVector predicted;
      filter_and_param_estimation();        //Filter attenuations by doing inequality projection and return 
                                            //the biexponential model approximation to m_current_data
      compute_ODF_SH_coeff();               //Fit the model to the original approximation
      m_coeff_samples.Column(1)=m_coeff;    //Save the estimates for the SH coefficients (deterministic estimates)
      predicted=m_HAT*m_current_data;       //Get the predicted "data"
      m_residuals=m_current_data-predicted; //And then the residuals
      modify_residuals();                   //Modify the residuals
      m_current_residuals=m_residuals;

      for (int n=2; n<=n_samples; n++){
	bootstrap_residuals();                //Get boostrapped residuals
	m_current_data=predicted+m_current_residuals; //Get the new "data" state
	compute_ODF_SH_coeff();               //Fit the model to the current "data" state
	m_coeff_samples.Column(n)=m_coeff;    //Save the estimates for the SH coefficients
      } 
    } 


    void modify_residuals() {
      for (int m=1; m<=m_residuals.Nrows(); m++){
	m_residuals(m)/=sqrt(1-m_HAT(m,m));   //Correct for heteroscedasticity
      }
    }
  

    void bootstrap_residuals() {
      for (int m=1; m<=m_residuals.Nrows(); m++){
	if (unifrnd().AsScalar()<0.5)         //Change the sign with probability 0.5
	  m_current_residuals(m)=-m_residuals(m); 
	else
	  m_current_residuals(m)=m_residuals(m); 
      }
    }

  };



  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Class that computes generalised FA value for the mean ODF obtained across N_samples of ODF coefficients.
  //A sphere tesselation is provided to indicate points on the sphere and ODF values are computed at these points.
  //GFA is then computed numerically as in (Tuch, 2004).
  class ODF_GFA{
    const Matrix& SHcoeff_samples;        //SH coefficients that describe the ODF (Num_coeff x N_samples)
    __attribute__((unused)) const int& lmax;                      //maximum SH order employed
    const int Num_of_coeff;               //number of SH coefficients 
    __attribute__((unused)) const int& nsamples;                  //number of bootstrap samples
    const Matrix& Eval_SH;                //(i,j) element: for each tesselation point i, contains the jth spherical harmonic evaluated at i (Matrix common to all voxels)

    float gfa;                            //holds the generalised FA value, corresponding to the mean ODF shape 
 
  public:
    //Constructor
    ODF_GFA(const Matrix& coeff, const int& SH_order, const int& nsamp,  const Matrix& EvalSH): SHcoeff_samples(coeff), lmax(SH_order), Num_of_coeff((SH_order+1)*(SH_order+2)/2), nsamples(nsamp), Eval_SH(EvalSH) {
      gfa=0;
    }

    //Destructor
    ~ODF_GFA(){}

    float get_gfa(){ return gfa; }
			  
    //Compute numerically the generalised FA (as in Tuch,2004), using the coefficients of the mean ODF across samples
    void compute_gfa(){
      int num_points=Eval_SH.Nrows();
      ColumnVector func_vals(num_points);
      ColumnVector m_avg_coeff;
      m_avg_coeff.ReSize(Num_of_coeff); //compute the mean SH coefficients across samples
      for (int i=1; i<=Num_of_coeff; i++)        
	m_avg_coeff(i)=SHcoeff_samples.Row(i).Sum()/SHcoeff_samples.Ncols(); 

      func_vals=0;
      for (int i=1; i<=num_points; i++)	      //compute the mean ODF value at each tesselation point
	for (int j=1; j<=Num_of_coeff; j++)
	  func_vals(i)+=m_avg_coeff(j)*Eval_SH(i,j);

      float meanf=func_vals.Sum()/num_points; float sum1=0, sum2=0;
      for (int i=1; i<=num_points; i++){	      //compute the gfa
      	sum1+=(func_vals(i)-meanf)*(func_vals(i)-meanf);
      	sum2+=func_vals(i)*func_vals(i);
      }
      gfa=sqrt(num_points*sum1/((num_points-1)*sum2));
    }

  };



  /////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////
  //Class to record and save samples
  class Samples{
    //for storing samples
    vector<Matrix> coeff_samples;            //The vector will be Num_of_coeff x nsamples x nvoxels 
  
    vector<Matrix> th_samples;               //The vectors will be max_num_of_peaks x nsamples x nvoxels
    vector<Matrix> ph_samples;               
    vector<Matrix> f_samples;               
    
    //for storing means
    vector<RowVector> mean_fsamples;         //max_num_of_peaks x nvoxels
    vector<Matrix> dyadic_vectors;           //max_num_of_peaks x 3 x nvoxels
    Matrix mean_SHcoeff;                     //Num_of_coeff x nvoxels, saves for each voxel the mean ODF coefficient across all ODF samples
    RowVector gfa;                           //If requested, 1 x nvoxels vectors that contains the GFA value for each shell 

    const volume< float >& m_mask;           //The brain mask that will be used to create Output images 
    const int& n_samples;                    //Number of samples
    const int& n_coeff;                      //Number of coefficients to store
    const int& n_voxels;                     //Number of valid non-zero voxels
    int zero_max_num_peaks;                  //A kludge
    const int& max_num_peaks;                //Maximum number of peaks allowed in each voxel
    const bool& m_meancoeff;                 //Flag to indicate whether the mean ODF coefficients will be saved 
    const bool& m_gfa;                       //Flag to indicate whether the GFA for the mean ODF will be saved 


  public:
  //Constructor (for saving coefficient samples)
  Samples(const int& ncoeff, const volume<float>& mask, const int& nvoxels, const int& nsamples, const bool& meancoeff=false, const bool& gfa_flag=false):
    m_mask(mask),n_samples(nsamples),n_coeff(ncoeff),n_voxels(nvoxels),zero_max_num_peaks(0),max_num_peaks(zero_max_num_peaks),m_meancoeff(meancoeff), m_gfa(gfa_flag){
      if (m_meancoeff){
	mean_SHcoeff.ReSize(ncoeff,nvoxels);  
	mean_SHcoeff=0;
      }
      else{                              //If m_meancoeff is set, do not save all coefficient samples 
	Matrix temp; 
	temp.ReSize(nsamples,nvoxels);   //Initialize structure that will hold the samples (n_coeff x n_samples x n_voxels)
	temp=0;
	for (int f=0; f<ncoeff; f++)
	  coeff_samples.push_back(temp);
      }
      if (m_gfa){
	gfa.ReSize(nvoxels);
	gfa=0;
      }
    }


  //Constructor (for saving directions samples)
  Samples(const volume<float>& mask, const int& nvoxels, const int& nsamples, const int& maxNumPeaks, const int& ncoeff, const bool& meancoeff=false, const bool& gfa_flag=false):
    m_mask(mask),n_samples(nsamples),n_coeff(ncoeff),n_voxels(nvoxels),max_num_peaks(maxNumPeaks), m_meancoeff(meancoeff), m_gfa(gfa_flag){  
      Matrix temp; RowVector temp2;
      temp.ReSize(nsamples,nvoxels);   //Initialize structure that will hold the samples (n_samples x n_voxels)
      temp=0;
      temp2.ReSize(nvoxels);
      temp2=0;
      for (int f=0; f<maxNumPeaks; f++){
	th_samples.push_back(temp);
	ph_samples.push_back(temp);
	f_samples.push_back(temp);
	mean_fsamples.push_back(temp2);
      }
      temp.ReSize(3,nvoxels); temp=0;
      for (int f=0; f<maxNumPeaks; f++)
      	dyadic_vectors.push_back(temp);
      if (m_meancoeff){
	mean_SHcoeff.ReSize(ncoeff,nvoxels);  
	mean_SHcoeff=0;
      }
      if (m_gfa){
	gfa.ReSize(nvoxels);
	gfa=0;
      }
    }

  //Destructor
  ~Samples(){}
  

  //Return the principle eigenvector of a dyadic tensor
  void Dyadic_tensor_e1(const SymmetricMatrix& dyad, ColumnVector& e1){
    DiagonalMatrix dyad_L; //eigenvalues
    Matrix dyad_V; //eigenvectors
      
    EigenValues(dyad,dyad_L,dyad_V);
    int maxeig;
    if(dyad_L(1)>dyad_L(2)){
      if(dyad_L(1)>dyad_L(3)) 
	maxeig=1;
      else 
	maxeig=3;
    }
    else{
      if(dyad_L(2)>dyad_L(3)) 
	maxeig=2;
      else 
	maxeig=3;
    }
    e1(1)=dyad_V(1,maxeig);  e1(2)=dyad_V(2,maxeig);  e1(3)=dyad_V(3,maxeig);
  }


  //Record the coefficient samples obtained for voxel "vox"
  //The input Matrix is Num_coeff x n_samples
  void record(const Matrix& Coeff, const int vox){
    if (m_meancoeff){
      for (int p=1; p<=n_coeff; p++)         //Get the average SH coefficients across samples 
	mean_SHcoeff(p,vox)=Coeff.Row(p).Sum()/Coeff.Ncols();
    }
    else{
      for (int p=0; p<n_coeff; p++)
	for (int n=1; n<=n_samples; n++)
	  coeff_samples[p](n,vox)=Coeff(p+1,n);
    }
  } 


  //Record the GFA obtained for voxel "vox"
  void record(const float gfaVal, const int vox){
    gfa(vox)=gfaVal;
  } 

  
  //Record the direction samples obtained for voxel "vox" and get the average dyads and fs
  //The input matrices are (max_num_peaks) x (n_samples+1) - the last entry is the mean ODF peaks  
  //Coeff is Num_coeff x n_samples, used to store the mean ODF coefficients if requested
  void record(const Matrix& Coeff, Matrix& theta, Matrix& phi, Matrix& f, const int vox){
    if (m_meancoeff){
      for (int p=1; p<=n_coeff; p++)         //Get the average SH coefficients across samples 
	mean_SHcoeff(p,vox)=Coeff.Row(p).Sum()/Coeff.Ncols();
    }
   
    //Sort the peaks across samples using their angular distance
    sort_peaks(theta, phi, f);

    //Record peaks
    for (int p=0; p<max_num_peaks; p++)
      for (int n=1; n<=n_samples; n++){
	th_samples[p](n,vox)=theta(p+1,n);
        ph_samples[p](n,vox)=phi(p+1,n);
        f_samples[p](n,vox)=f(p+1,n);
      }
 
    //Get the mean values using the mean ODF shape (which is last in the Matrices entries)
    for (int p=0; p<max_num_peaks; p++){
      mean_fsamples[p](vox)=f(p+1,n_samples+1);   
      if (f(p+1,n_samples+1)!=0){
	dyadic_vectors[p](1,vox)=sin(theta(p+1,n_samples+1))*cos(phi(p+1,n_samples+1));
	dyadic_vectors[p](2,vox)=sin(theta(p+1,n_samples+1))*sin(phi(p+1,n_samples+1));
	dyadic_vectors[p](3,vox)=cos(theta(p+1,n_samples+1));
      }
      else{
	dyadic_vectors[p](1,vox)=0; 
	dyadic_vectors[p](2,vox)=0; 
	dyadic_vectors[p](3,vox)=0;
      }
    }

    /*
    //After sorting peaks, compute the mean dyad and mean f for each peak
    for (int p=0; p<max_num_peaks; p++){
      dyad=0;
      int flag=0;
      for (int n=1; n<=n_samples; n++)
	if (f(p+1,n)!=0){
	  mean_fsamples[p](vox)+=f(p+1,n);   //Add each sample to the mean f
	  tempvec << sin(theta(p+1,n))*cos(phi(p+1,n)) << sin(theta(p+1,n))*sin(phi(p+1,n)) << cos(theta(p+1,n));
	  dyad << dyad + tempvec*tempvec.t();   //Create mean dyadic tensor of all direction samples
	  flag=1;
	}
      if (flag!=0){    //if there are samples for this population
	mean_fsamples[p](vox)/=n_samples;
	Dyadic_tensor_e1(dyad, tempvec);
	dyadic_vectors[p](1,vox)=tempvec(1);
	dyadic_vectors[p](2,vox)=tempvec(2);
	dyadic_vectors[p](3,vox)=tempvec(3);
	} 
	} */
  }  



  //Sort peaks across samples using their absolute dot product
  //The input matrices are (max_num_peaks) x (n_samples+1) - the last entry is the mean ODF peaks  
  void sort_peaks(Matrix& theta, Matrix& phi, Matrix& f){
    ColumnVector my_temp; 
    const int nsamp=theta.Ncols();
    Matrix vecx(max_num_peaks,nsamp);  Matrix vecy(max_num_peaks,nsamp);  Matrix vecz(max_num_peaks,nsamp);
    Matrix theta_new(max_num_peaks,nsamp);  Matrix phi_new(max_num_peaks,nsamp);  Matrix f_new(max_num_peaks,nsamp);
    int max_ref_count=0;   int ref_samp=1;

    //Choose a reference sample (exclude the mean sample, which is last) 
    for (int n=1; n<=nsamp-1; n++){
      int count=0;
      for (int p=1; p<=max_num_peaks; p++)  //For each sample, count how many peaks exist
       	if (f(p,n)!=0) count++;
      if (count>max_ref_count){             //Keep the sample index with the maximum non-zero peaks
	max_ref_count=count;                     
	ref_samp=n;
      }
    }

    //Swap the reference sample with the first  
    if (ref_samp!=1){
      my_temp=theta.Column(1); theta.Column(1)=theta.Column(ref_samp); theta.Column(ref_samp)=my_temp;
      my_temp=phi.Column(1); phi.Column(1)=phi.Column(ref_samp); phi.Column(ref_samp)=my_temp;
      my_temp=f.Column(1); f.Column(1)=f.Column(ref_samp); f.Column(ref_samp)=my_temp;
    }
    
    vecx=0; vecy=0; vecz=0;
    //Convert (theta,phi) to vectors with cartesian coordinates
    for (int n=1; n<=nsamp; n++)
      for (int p=1; p<=max_num_peaks; p++)
	if (f(p,n)!=0){
	  vecx(p,n)=sin(theta(p,n))*cos(phi(p,n)); 
	  vecy(p,n)=sin(theta(p,n))*sin(phi(p,n)); 
	  vecz(p,n)=cos(theta(p,n));
	}
    
    //Sort the remaining samples using the angular agreement with the reference peaks
   Matrix dots(max_num_peaks,max_num_peaks);
   Matrix valid(max_num_peaks,max_num_peaks);
   
   theta_new=0; phi_new=0; f_new=0;
   theta_new.Column(1)=theta.Column(1);   phi_new.Column(1)=phi.Column(1);   f_new.Column(1)=f.Column(1);  //save the reference peaks
  
   for (int n=2; n<=nsamp; n++){           //For each sample
     dots=0; 
     valid=1;
     for (int q=1; q<=max_num_peaks; q++){     //For each non-zero sample peak
       if (f(q,n)!=0)
	 for (int p=1; p<=max_num_peaks; p++){  //Get the dot product with each reference peak and store to a matrix
	   if (f(p,1)!=0)
	     dots(p,q)=fabs(vecx(p,1)*vecx(q,n)+vecy(p,1)*vecy(q,n)+vecz(p,1)*vecz(q,n)); 
	   else
	     valid.Row(p)=0;
	 }
       else
	 valid.Column(q)=0;                   //Indicate which entries of the matrix dots are not valid, because peaks have zero value 
     }

     for (int q=1; q<=max_num_peaks; q++){
       int sp,sq;
       dots.Maximum2(sp,sq);        //Return the location of the maximum dot product
       if (valid(sp,sq)==1){
	 theta_new(sp,n)=theta(sq,n);            //Swap peak location     
	 phi_new(sp,n)=phi(sq,n);     
	 f_new(sp,n)=f(sq,n);         
	 dots.Column(sq)=0;                      //Neglect the other dot products of this sample peak with the reference ones
	 dots.Row(sp)=0;                         //and this reference peak, as it has been paired
	 valid.Column(sq)=0; valid.Row(sp)=0;
       }
     }
   }
   
  /* for (int n=2; n<=nsamp; n++)             //For each sample
      for (int p=1; p<=max_num_peaks; p++){  //For each reference peak 
	maxdot=-1; pos=1;
	for (int q=1; q<=max_num_peaks; q++) //Get the dot product with all non-zero sample peaks
	  if (f(q,n)!=0){
	    dot=fabs(vecx(p,1)*vecx(q,n)+vecy(p,1)*vecy(q,n)+vecz(p,1)*vecz(q,n)); 
	    if (dot>=maxdot){                //and pick the closest one
	      maxdot=dot;
	      pos=q;                            
	    }
	  }
	  theta_new(p,n)=theta(pos,n);
	  phi_new(p,n)=phi(pos,n);
	  f_new(p,n)=f(pos,n);
	  f(pos,n)=0;                      //Indicate that this sample peak has been already assigned 
	  } */
   theta=theta_new; phi=phi_new; f=f_new;  
  }



  //Save all recorded coefficient samples in images
  void save_coeff(){
    volume4D<float> tmp;
    Log& logger=LogSingleton::getInstance();

    if (m_meancoeff){                  //Save mean coefficients in a single 4D image
      tmp.setmatrix(mean_SHcoeff,m_mask);
      string oname="meanSHcoeff";
      save_volume4D(tmp,logger.appendDir(oname)); 
    }
    else{
      for(int m=0; m<n_coeff; m++){   //For each coefficient create a new 4D image
	tmp.setmatrix(coeff_samples[m],m_mask);
	string oname="c"+num2str(m+1)+"samples";
	save_volume4D(tmp,logger.appendDir(oname));
      }
    }
    if (m_gfa){
      tmp.setmatrix(gfa,m_mask);
      string oname="GFA";
      save_volume4D(tmp,logger.appendDir(oname)); 
    }
  }  
 

  //Save all recorded direction samples in images
  void save_dir(){
    volume4D<float> tmp;
    Matrix thsamples_tmp(max_num_peaks,n_samples);     //temporary structures used when rearranging peaks
    Matrix phsamples_tmp(max_num_peaks,n_samples);    
    Matrix fsamples_tmp(max_num_peaks,n_samples);    
    Matrix dyadic_vectors_tmp(max_num_peaks,3);    
    RowVector mean_fsamples_tmp(max_num_peaks);

    Log& logger=LogSingleton::getInstance();
    
    //Sort the output according to the value of the mean ODF across the samples (mean_fsamples) (i.e. determine which population is 1,2,3...)
    for(int vox=1;vox<=n_voxels; vox++){  //for each voxel
      vector<pair<float,int> > sfs;
      pair<float,int> ftmp;
      
      for(int f=0;f<max_num_peaks;f++){  
	ftmp.first=mean_fsamples[f](vox);//this is the value according which sorting occurs
	ftmp.second=f;                   //this is the index of the peak
	sfs.push_back(ftmp);
      }
      sort(sfs.begin(),sfs.end());       //sort in ascending order
      
      //Rearrange peaks 
      for(int f=0;f<max_num_peaks;f++){
	for(int samp=1;samp<=n_samples;samp++){
     	  thsamples_tmp(f+1,samp)=th_samples[sfs[(sfs.size()-1)-f].second](samp,vox);
	  phsamples_tmp(f+1,samp)=ph_samples[sfs[(sfs.size()-1)-f].second](samp,vox);
	  fsamples_tmp(f+1,samp)=f_samples[sfs[(sfs.size()-1)-f].second](samp,vox);
	}
	mean_fsamples_tmp(f+1)=mean_fsamples[sfs[(sfs.size()-1)-f].second](vox);
	dyadic_vectors_tmp(f+1,1)=dyadic_vectors[sfs[(sfs.size()-1)-f].second](1,vox);
	dyadic_vectors_tmp(f+1,2)=dyadic_vectors[sfs[(sfs.size()-1)-f].second](2,vox);
	dyadic_vectors_tmp(f+1,3)=dyadic_vectors[sfs[(sfs.size()-1)-f].second](3,vox);
      }

      for(int f=0;f<max_num_peaks;f++){
	for(int samp=1;samp<=n_samples;samp++){
     	  th_samples[f](samp,vox)=thsamples_tmp(f+1,samp);
     	  ph_samples[f](samp,vox)=phsamples_tmp(f+1,samp);
     	  f_samples[f](samp,vox)=fsamples_tmp(f+1,samp);
	}
	mean_fsamples[f](vox)=mean_fsamples_tmp(f+1);
	dyadic_vectors[f](1,vox)=dyadic_vectors_tmp(f+1,1);
	dyadic_vectors[f](2,vox)=dyadic_vectors_tmp(f+1,2);
	dyadic_vectors[f](3,vox)=dyadic_vectors_tmp(f+1,3);
      }
    } 

    //Save sorted fibres 
    for(int m=0; m<max_num_peaks; m++){ //For each coefficient create a new 4D image
      tmp.setmatrix(th_samples[m],m_mask);
      string oname="merged_th"+num2str(m+1)+"samples";
      save_volume4D(tmp,logger.appendDir(oname));
        
      tmp.setmatrix(ph_samples[m],m_mask);
      oname="merged_ph"+num2str(m+1)+"samples";
      save_volume4D(tmp,logger.appendDir(oname));
      
      tmp.setmatrix(f_samples[m],m_mask);
      oname="merged_f"+num2str(m+1)+"samples";
      save_volume4D(tmp,logger.appendDir(oname)); 
    
      tmp.setmatrix(dyadic_vectors[m],m_mask);
      oname="dyads"+num2str(m+1);
      save_volume4D(tmp,logger.appendDir(oname));
    
      tmp.setmatrix(mean_fsamples[m],m_mask);
      oname="mean_f"+num2str(m+1)+"samples";
      save_volume4D(tmp,logger.appendDir(oname)); 
    }

    if (m_meancoeff){                  //Save mean coefficients in a single 4D image
      tmp.setmatrix(mean_SHcoeff,m_mask);
      string oname="mean_SHcoeff";
      save_volume4D(tmp,logger.appendDir(oname)); 
    } 

    if (m_gfa){
      tmp.setmatrix(gfa,m_mask);
      string oname="GFA";
      save_volume4D(tmp,logger.appendDir(oname)); 
    }
  }
  
  };

 
  /////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////
  //Class that perfors operations common for all voxels. These include preparing the data
  //computing the design and Hat matrix for spherical harmonic decomposition and scaling
  //factors for transforming signal SH to ODF SH coefficients  
  class ODF_Volume_Manager{
    qbootOptions& opts;          //Options set during function call
    Matrix m_SHT;                //Pseudoinverse of spherical harmonic transform matrix, used to get the signal SH coefficients (Num_coeff x Num_datapoints) 
    Matrix m_HAT;                //Hat matrix used for residual-bootstrap (Num_datapoints x Num_datapoints)
    DiagonalMatrix m_scaling_factors;//Holds the scaling factors to obtain the ODF SH coefficients from the signal SH ones
    ColumnVector m_theta;        //Euler Angles of bvecs
    ColumnVector m_phi;
    
    Matrix m_data;               //Matrix with signal attenuations
    Matrix m_bvecs;              //bvecs and bvals
    Matrix m_bvals;
    volume<float> m_mask;        //The brain mask that will be used to create Output images 

    const int npeaks;            //Maximum number of ODF peaks detected 
    const float peak_threshold;  //Mimimum ODF value for a point to be considered a peak
    const int m_modelnum;        //1 for Descoteaux's ODFs, 2 for CSA-ODFs, 3 for multi-shell CSA ODFs
    const int lmax;              //maximum SH order (must be even)
    const int Num_of_coeff;      //Number of SH coefficients to be estimated
    const int nsamples;          //Number of bootstrap samples
    const float m_lambda;        //Laplace-Beltrami regularization parameter
    const float m_delta;         //Regularization parameter for signal attenuations utilized in CSA-ODF estimation
    const float m_alpha;         //Laplacian sharpening parameter for model=1
    const bool m_savecoeff;      //Flag to indicate whether the ODF SH coeff or the ODF peaks (default) will be saved
    const bool m_savemeancoeff;  //Flag to indicate whether the mean ODF SH coeff (along with the peaks) will be saved
    const bool m_gfa_flag;       //Flag to indicate whether the GFA will be computed and saved for each voxel

  public:
    //constructor
  ODF_Volume_Manager():
    opts(qbootOptions::getInstance()),npeaks(opts.npeaks.value()), peak_threshold(opts.peak_threshold.value()), m_modelnum(opts.modelnum.value()), lmax(opts.lmax.value()), Num_of_coeff((lmax+1)*(lmax+2)/2), 
      nsamples(opts.nsamples.value()), m_lambda(opts.lambda.value()), m_delta(opts.delta.value()), m_alpha(opts.alpha.value()), m_savecoeff(opts.savecoeff.value()), m_savemeancoeff(opts.savemeancoeff.value()), m_gfa_flag(opts.gfa.value()) {
    }

    //destructor
    ~ODF_Volume_Manager(){}  


    //Read bvals,bvecs, data and nodif_brain_mask
    void Read_data(){
      if (opts.verbose.value())
	cout<<"reading data, bvals and bvecs"<<endl;
      m_bvals=read_ascii_matrix(opts.bvalsfile.value());    //Read bvals and bvecs
      m_bvecs=read_ascii_matrix(opts.bvecsfile.value());
      if(m_bvecs.Nrows()>3) m_bvecs=m_bvecs.t();
      if(m_bvals.Nrows()>1) m_bvals=m_bvals.t();
      for(int i=1;i<=m_bvecs.Ncols();i++){
	float tmpsum=sqrt(m_bvecs(1,i)*m_bvecs(1,i)+m_bvecs(2,i)*m_bvecs(2,i)+m_bvecs(3,i)*m_bvecs(3,i));
	if(tmpsum!=0){
	  m_bvecs(1,i)=m_bvecs(1,i)/tmpsum;
	  m_bvecs(2,i)=m_bvecs(2,i)/tmpsum;
	  m_bvecs(3,i)=m_bvecs(3,i)/tmpsum;
	}  
      }
      volume4D<float> data_vol;                            //Read data and brain mask
      read_volume4D(data_vol,opts.datafile.value());
      read_volume(m_mask,opts.maskfile.value());
      m_data=data_vol.matrix(m_mask);  
    }

    
    //Convert signals to signal attenuations. Also Reshape the data when multiple shells are available.
    void Prepare_data_for_ODFs(){
      int count;
      float bval3=0, bval1=0;
      bool multishell=false;

      for (int n=1; n<=m_bvals.Ncols(); n++)   //Check whether more than one non-zero b-values exist 
	if (m_bvals(1,n)>b0max)
	  bval3=m_bvals(1,n);                  //The last b-value will be stored
      
      for (int n=1; n<=m_bvals.Ncols(); n++)
	if (m_bvals(1,n)>b0max && fabs(m_bvals(1,n)-bval3)>2*shellrange){ //Then we have a b value from a different shell
	  bval1=m_bvals(1,n);
	  multishell=true;
	  break;
	} 
      
      
      if (!multishell){ //Single-shell analysis
	count=Sig2SigAttenuation(m_data, m_bvals,m_bvecs); //Convert signal to signal attenuation
	if (m_modelnum!=1 && m_modelnum!=2){
	  cerr<<"Wrong modelnum chosen! Choose 1 or 2 for single-shell data!"<<endl<<"Exiting now!"<<endl;
	  exit(1);
	}
	if (count<Num_of_coeff){
	  cerr<<"Wrong lmax chosen! "<<Num_of_coeff<<" coefficients need to be estimated, but only "
	      <<count<<" datapoints are available!"<<endl<<"Exiting now!"<<endl;
	  exit(1);
	}
      }
      else{  //Multi-shell analysis
	cout<<"Multiple b-values have been detected. Three shells are assumed, acquired one after the other!"<<endl;
      	float bval2=0;
	for (int n=1; n<=m_bvals.Ncols(); n++)  //Store the second b-value
	  if (m_bvals(1,n)>b0max && fabs(m_bvals(1,n)-bval1)>2*shellrange && fabs(m_bvals(1,n)-bval3)>2*shellrange){
	    bval2=m_bvals(1,n);
	    break;
	  } 

	//Find the mean b-value per shell
	float btemp=0; int bcnt=0;
	for (int n=1; n<=m_bvals.Ncols(); n++) 
	  if (m_bvals(1,n)>b0max && fabs(m_bvals(1,n)-bval1)<=2*shellrange){
	    btemp+=m_bvals(1,n); bcnt++;
	  } 
	bval1=btemp/bcnt;

	btemp=0; bcnt=0;
	for (int n=1; n<=m_bvals.Ncols(); n++) 
	  if (m_bvals(1,n)>b0max && fabs(m_bvals(1,n)-bval2)<=2*shellrange){
	    btemp+=m_bvals(1,n); bcnt++;
	  } 
	bval2=btemp/bcnt;

	btemp=0; bcnt=0;
	for (int n=1; n<=m_bvals.Ncols(); n++) 
	  if (m_bvals(1,n)>b0max && fabs(m_bvals(1,n)-bval3)<=2*shellrange){
	    btemp+=m_bvals(1,n); bcnt++;
	  } 
	bval3=btemp/bcnt;
	cout<<"Mean b-values for detected shells are "<<bval1<<", "<<bval2<<" and "<<bval3<<endl;
      

	//Assume the same number of samples in each shell. If not, need to provide a text file with the number of samples in each shell
	vector<int> shell_num(3);  //Keep the Number of samples in each shell (including the b=0)
	if (opts.qshellsfile.value().empty()){ //If no qshells file provided, assume same number of directions in each shell
	  shell_num[0]=(int) (m_bvals.Ncols()/3.0); shell_num[1]=(int) (m_bvals.Ncols()/3.0); shell_num[2]=(int) (m_bvals.Ncols()/3.0); 
	  if (opts.verbose.value())
	    cout<<"No qshells file provided. Assuming same number of directions in each shell!"<<endl;;
	}
	else{  //A qshell file has 3 entries in a row, indicating hoe many samples were taken at each shell.
	  Matrix qshells=read_ascii_matrix(opts.qshellsfile.value());
	  shell_num[0]=(int)qshells(1,1); shell_num[1]=(int)qshells(1,2); shell_num[2]=(int)qshells(1,3);
	  if (opts.verbose.value())
	    cout<<"Reading qshells file..."<<endl;;
	}
	
	//Get the signal attenuations for each shell. Do it individually, as the b=0 may have different intensity scalings for each shell!
	vector<Matrix> data_shell, bvecs_shell, bvals_shell;
	//Copying shell1
	data_shell.push_back(m_data.SubMatrix(1,shell_num[0],1,m_data.Ncols()));              
	bvecs_shell.push_back(m_bvecs.SubMatrix(1,3,1,shell_num[0]));    
	bvals_shell.push_back(m_bvals.SubMatrix(1,1,1,shell_num[0]));
	//Copying shell2
	data_shell.push_back(m_data.SubMatrix(1+shell_num[0],shell_num[0]+shell_num[1],1,m_data.Ncols())); 
	bvecs_shell.push_back(m_bvecs.SubMatrix(1,3,1+shell_num[0],shell_num[0]+shell_num[1])); 
	bvals_shell.push_back(m_bvals.SubMatrix(1,1,1+shell_num[0],shell_num[0]+shell_num[1]));
	//Copying shell3
	data_shell.push_back(m_data.SubMatrix(1+shell_num[0]+shell_num[1],shell_num[0]+shell_num[1]+shell_num[2],1,m_data.Ncols()));     
	bvecs_shell.push_back(m_bvecs.SubMatrix(1,3,1+shell_num[0]+shell_num[1],shell_num[0]+shell_num[1]+shell_num[2]));  
	bvals_shell.push_back(m_bvals.SubMatrix(1,1,1+shell_num[0]+shell_num[1],shell_num[0]+shell_num[1]+shell_num[2]));

	for (int n=0; n<3; n++)
	  Sig2SigAttenuation(data_shell[n],bvals_shell[n],bvecs_shell[n]); 
     
	//Check whether interpolation is needed (i.e. whether different number of directions or directions between shells, then need to interpolate)
	bool interp_flag=false; 	int max, max_shell=0, min, min_shell=0;
	if (bvecs_shell[0].Ncols()!=bvecs_shell[1].Ncols() || bvecs_shell[0].Ncols()!=bvecs_shell[2].Ncols() || bvecs_shell[1].Ncols()!=bvecs_shell[2].Ncols())
	  interp_flag=true;
	else{ //if same number of directions, check whether the directions between shells are different (consider them different only if they differ more than 1 degree)
	  for (int i=1; i<=bvecs_shell[0].Ncols(); i++)
	    if (fabs(dot(bvecs_shell[0].Column(i),bvecs_shell[1].Column(i)))<=0.9998) {interp_flag=true; break;}
	  for (int i=1; i<=bvecs_shell[0].Ncols(); i++)
	    if (fabs(dot(bvecs_shell[0].Column(i),bvecs_shell[2].Column(i)))<=0.9998) {interp_flag=true; break;}
	  for (int i=1; i<=bvecs_shell[1].Ncols(); i++)
	    if (fabs(dot(bvecs_shell[1].Column(i),bvecs_shell[2].Column(i)))<=0.9998) {interp_flag=true; break;}
	} 

	if (interp_flag){
	  //Find which shell has the fewer directions
	  if (bvecs_shell[0].Ncols()<=bvecs_shell[1].Ncols()){ min=bvecs_shell[0].Ncols(); min_shell=0;}
	  else{ min=bvecs_shell[1].Ncols(); min_shell=1;}
	  if (min>bvecs_shell[2].Ncols()) min_shell=2;

	  //Find which shell has the most directions
	  if (bvecs_shell[0].Ncols()>=bvecs_shell[1].Ncols()){ max=bvecs_shell[0].Ncols(); max_shell=0;}
	  else{ max=bvecs_shell[1].Ncols(); max_shell=1;}
	  if (max<bvecs_shell[2].Ncols()) max_shell=2;

	  //Determine which order will be used for interpolation (common to all shells). Use 6, unless the data are not enough 
	  int lmax_interp=6; //10 
	  while (((lmax_interp+1)*(lmax_interp+2)/2)>bvecs_shell[min_shell].Ncols()/3.0 && lmax_interp>lmax)
	    lmax_interp-=2;

	  if (opts.verbose.value())
	    cout<<"Interpolating data in each shell using spherical harmonics of order "<<lmax_interp<<endl;

	  //Interpolate all shells along these directions
	  Interpolate_Shell(data_shell[0], bvecs_shell[0], bvecs_shell[max_shell],lmax_interp);
	  Interpolate_Shell(data_shell[1], bvecs_shell[1], bvecs_shell[max_shell],lmax_interp);
	  Interpolate_Shell(data_shell[2], bvecs_shell[2], bvecs_shell[max_shell],lmax_interp);
	}

	m_data<< (data_shell[0] & data_shell[1] & data_shell[2]);
	m_bvecs<< bvecs_shell[max_shell]; //These are the directions common to all shells
	if (m_modelnum==2){  //Get the geometric mean of the data
	  cout<<"Monoexponential model chosen, data from multiple shells will be averaged."<<endl;
	  float exp1=1.0/bval1, exp2=1.0/bval2, exp3=1.0/bval3, exp4=1.0/(exp1+exp2+exp3);  
	  for(int vox=1;vox<=m_data.Ncols();vox++)      //For each voxel
	    for (int i=1; i<=m_bvecs.Ncols(); i++){     //For each direction, get the geometric mean of the measurements
	      float geom_mean=pow((pow(data_shell[0](i,vox),exp1)*pow(data_shell[1](i,vox),exp2)*pow(data_shell[2](i,vox),exp3)),exp4);
	      data_shell[0](i,vox)=geom_mean;
	  }
	  m_data<<data_shell[0];
	}
	else if (m_modelnum==3){}  //Do nothing for now
	else{
	  cerr<<"Wrong modelnum chosen! Choose 2 or 3 for multi-shell data!"<<endl<<"Exiting now!"<<endl;
	  exit(1);
	}
      }
    }



    //For each voxel of a data matrix Num_of_dir x Num_of_voxels, interpolate the data measured along measured_bvecs at points 
    //indicated by target_bvecs. New data are returned to the same Matrix, Interpolate using spherical harmonics up to order l_max;
    void Interpolate_Shell(Matrix& data, const Matrix& measured_bvecs, const Matrix& target_bvecs, const int& l_max) const{
      ColumnVector measured_theta, measured_phi, target_theta, target_phi;
      Matrix tempSHT, SHT, interpSHT,new_data(target_bvecs.Ncols(),data.Ncols());

      cart2sph(measured_bvecs,measured_theta,measured_phi);
      cart2sph(target_bvecs,target_theta,target_phi);
      
      tempSHT=fill_SH_matrix(measured_theta,measured_phi,l_max);   //Fill Matrix SHT, the design matrix for the spherical harmonic transform
      SHT=pinv(tempSHT);   
      interpSHT=fill_SH_matrix(target_theta,target_phi,l_max);     //Fill the interpolation matrix 
 
      for (int vox=1; vox<=data.Ncols(); vox++){
	ColumnVector vox_data=data.Column(vox);
	//ColumnVector coeff_signal = SHT*vox_data;                //Get the SH coefficients of the signal
	new_data.Column(vox)= interpSHT* SHT *vox_data;            //Get the new interpolated data
      }
      data << new_data;
    }



    //For a set of data, remove all b=0 entries from bvals and bvecs and the data. Divide each DW signal with the mean S0 to 
    //get the signal attenuation. Return the new data, bvals and bvecs (overwrite the old) and the number of non-b=0 DW directions.
    //b0s are detected as volumes with a low b value (smaller than b0max)
    int Sig2SigAttenuation(Matrix& data, Matrix& bvals, Matrix& bvecs) const{
      int count=0;
      RowVector S0;
      Matrix bvals_new, bvecs_new, data_new;
      int Nd=data.Nrows();   //Number of total datapoints 
      int Mv=data.Ncols();   //Total number of non-background voxels
      
      for (int n=1; n<=Nd; n++)
	if (bvals(1,n)>b0max)
	  count++;             //Count how many DWs have been acquired (excluding b=0's)
      
      if (count==Nd){
	cerr<<"At least one b=0 image is required! Exiting now!"<<endl;
	exit(1);
      }

      bvals_new.ReSize(1,count); 
      bvecs_new.ReSize(3,count); 
      data_new.ReSize(count,Mv);
      data_new=0;
      S0.ReSize(Mv);
      S0=0;
      
      int m=1;
      for (int n=1; n<=Nd; n++){
	if (bvals(1,n)>b0max){
	  bvals_new(1,m)=bvals(1,n);    //Remove all b=0 entries from bvals and bvecs
	  bvecs_new.Column(m)=bvecs.Column(n);
	  for(int vox=1;vox<=Mv;vox++)
	    data_new(m,vox)=data(n,vox);
	  m++;
	}
	else{   //If it is a b=0
	  for(int vox=1;vox<=Mv;vox++)
	    S0(vox)+=data(n,vox);
	}
      }
      
      S0=S0/(Nd-count);                  //Get the average S0 intensity in each voxel
      bvals<<bvals_new;
      bvecs<<bvecs_new;
      
      for(int vox=1;vox<=Mv;vox++)       //Get the signal attenuation
	for (int m=1; m<=count; m++)
	  data_new(m,vox)/=S0(vox);
      data<<data_new;

      return count;
    }
 

    //Computes scaling factors that convert signal SH coefficients to ODF SH coefficients
    //model_num defines the model used (1 for Descoteaux's ODFs, 2 for CSA-ODFs)
    void compute_scaling_factors(){
      double prod, prod2, P;
      int j,n;

      j=(lmax+1)*(lmax+2)/2;
      m_scaling_factors.ReSize(j);
      switch (m_modelnum){
      case 1:
	P=2*M_PI;
	m_scaling_factors(1)=P;
	for (int l=2; l<=lmax; l=l+2)
	  for (int m=-l; m<=l; m++){
	    j=l*(l+1)/2+m+1;
	    prod=P*pow(-1,(float)l/2); 
	    prod2=1;
	    for (n=3; n<=l-1; n=n+2)
	      prod=prod*n;
	    for (n=2; n<=l; n=n+2)
	      prod2=prod2*n;
	    m_scaling_factors(j)=prod/prod2;
	  }
	if (m_alpha!=0)
	  for (int l=2; l<=lmax; l=l+2)
	    for (int m=-l; m<=l; m++){
	      j=l*(l+1)/2+m+1;
	      m_scaling_factors(j)*=(1-m_alpha*(-l)*(l+1)); //Laplace-Beltrami sharpening 
	    }                                               //according to (Descoteaux, 2005) (Subtract from the ODF a portion m_alpha of its Laplacian) 
	break;
      default:  //for modelnum 2 and 3
	P=2/(16*M_PI);
	m_scaling_factors(1)=0;
	for (int l=2; l<=lmax; l=l+2)
	  for (int m=-l; m<=l; m++){
	    j=l*(l+1)/2+m+1;
	    prod=P*pow(-1,(float)l/2); 
	    prod2=1;
	    for (n=3; n<=l-1; n=n+2)
	      prod=prod*n;
	    for (n=2; n<=l; n=n+2)
	      prod2=prod2*n;
	    m_scaling_factors(j)=-l*(l+1)*prod/prod2;
	  }
	break;
      }
    }


    //Computes the pseudoinverse of the spherical harmonic transform matrix
    //Only even SH are considered
    //The spherical harmonics basis is modified as in (Descoteaux,MRM,2007) to be symmetric and real
    void invert_SH_matrix(){
      Matrix tempSHT;
      
      tempSHT=fill_SH_matrix(m_theta,m_phi,lmax);   //Fill Matrix SHT, the design matrix for the spherical harmonic transform
      if (m_lambda!=0)                              //Apply Laplace-Beltrami Regularization during estimation
	m_SHT=LB_reg_pinv(tempSHT);
      else
	m_SHT=pinv(tempSHT);                        //Return the pseudoinverse of the matrix
      m_HAT=tempSHT*m_SHT;                          //This is the Hat matrix
    } 


    //Returns the pseudoinverse of tmp, but also performs Laplace-Beltrami regularization
    ReturnMatrix LB_reg_pinv(const Matrix& tmp) const{
      Matrix tmp_t, res;
      DiagonalMatrix L(Num_of_coeff);
      
      tmp_t=tmp.t();
      int j=1;
      for (int l=0; l<=lmax; l=l+2)
	for (int m=-l; m<=l; m++){
	  L(j)=l*l*(l+1)*(l+1);
	  j++;
	}
      res=(tmp_t*tmp+m_lambda*L).i()*tmp_t;
      res.Release();
      return res;
    }


    //Get points evenly distributed on the sphere using icosahedral tessellation.  
    //Also fill in structures that are common to all Peak_finder_finite_diff objects  
    void Init_Peak_finder_finite_diff(Matrix& index, Matrix& Eval_SH, ColumnVector& p_theta, ColumnVector& p_phi, ColumnVector& p_z) const {
      const int tessel_degree=5; //Tessellation degree for obtaining points on the sphere
      const int num_neighbours=20; //25 for tessel_degree=5, 6 for tessel_degree=4 (angular distance of 4.5 degrees)
 
      Points m_points(tessel_degree,num_neighbours);            //Get points on the sphere
      index = m_points.get_index();
      cart2sph(m_points.get_TessPoints_Ref(),p_theta, p_phi);   //For each point of the tesselation get spherical angles
      p_z = m_points.get_TessPoints_coord(3);                   //Save the z coordinate of each point
      Eval_SH=fill_SH_matrix(p_theta,p_phi,lmax);               //Compute the value of spherical harmonics for each point and store to a matrix
      // cout<<m_points.get_ang_dist()<<" "<<m_points.get_Eucl_dist()<<endl;
    }


    //Get points evenly distributed on the sphere using icosahedral tessellation.  
    //Return a Matrix that contains SHs evaluated at each point  
    void GFA_tessel_points(Matrix& Eval_SH) const {
      const int tessel_degree=3; //Tessellation degree for obtaining points on the sphere
      const int num_neighbours=1; //25 for tessel_degree=5, 6 for tessel_degree=4 (angular distance of 4.5 degrees)
      ColumnVector p_theta, p_phi;

      Points m_points(tessel_degree,num_neighbours);            //Get points on the sphere
      cart2sph(m_points.get_TessPoints_Ref(),p_theta, p_phi);   //For each point of the tesselation get spherical angles
      Eval_SH=fill_SH_matrix(p_theta,p_phi,lmax);               //Compute the value of spherical harmonics for each point and store to a matrix
    }



    //Performs ODF estimation for all voxels
    void run_all(){
      Read_data();
      Prepare_data_for_ODFs();   //Prepare data and matrices for processing
      cart2sph(m_bvecs,m_theta,m_phi);
      if (opts.verbose.value())
	cout<<"computing design matrix and scaling factors for ODF estimation"<<endl;
      compute_scaling_factors();
      invert_SH_matrix();

      if (m_savecoeff)   //Save all coefficients or, if requested, the mean ODF coefficients
	run_all_coeff();
      else              //Save peaks and, if requested, the mean ODF coefficients
	run_all_peaks();
    }


    void run_all_coeff(){
      Matrix Eval_SH;              //Matrix that is used for GFA calculation and holds the value of spherical harmonics at the points on the sphere
      Samples samp(Num_of_coeff, m_mask, m_data.Ncols(), nsamples, m_savemeancoeff,m_gfa_flag);
      float gfa;

      if (opts.verbose.value())
	cout<<"running inference for each voxel"<<endl;
      if (m_savecoeff && m_savemeancoeff)
	cout<<"--savemeancoeff has been set, only mean ODF coefficients will be saved"<<endl;

      if (m_gfa_flag)
	GFA_tessel_points(Eval_SH);
      
      for(int vox=1;vox<=m_data.Ncols();vox++){
	cout <<vox<<"/"<<m_data.Ncols()<<endl;
	Matrix SHcoeff;
	if (m_modelnum==3){
	  Matrix temp_data(m_bvecs.Ncols(),3);  //3 is the number of q-shells
	  temp_data.Column(1)=(m_data.Column(vox)).Rows(1,m_bvecs.Ncols());  //Get multi-shell data for that voxel
	  temp_data.Column(2)=(m_data.Column(vox)).Rows(m_bvecs.Ncols()+1,2*m_bvecs.Ncols());
	  temp_data.Column(3)=(m_data.Column(vox)).Rows(2*m_bvecs.Ncols()+1,3*m_bvecs.Ncols());
	  ODF_voxel_biexp local_ODF(temp_data, m_SHT, m_HAT, m_scaling_factors, lmax, Num_of_coeff, nsamples, m_delta);
	  local_ODF.bootstrapping();
	  SHcoeff=local_ODF.get_SHcoeff_ref();
	  samp.record(SHcoeff, vox);
	}
	else{
	  ColumnVector temp_data=m_data.Column(vox);  //data for that voxel
	  ODF_voxel local_ODF(temp_data, m_SHT, m_HAT, m_scaling_factors, m_modelnum, lmax, Num_of_coeff, nsamples,m_delta);
	  local_ODF.bootstrapping();
	  SHcoeff=local_ODF.get_SHcoeff_ref();
	  samp.record(SHcoeff, vox);
	}

	if (m_gfa_flag){
	  ODF_GFA local_vox(SHcoeff, lmax, nsamples, Eval_SH);
	  local_vox.compute_gfa();
	  gfa=local_vox.get_gfa();
	  samp.record(gfa,vox);
	}

      }
      if (opts.verbose.value())
	cout<<"saving samples"<<endl;
      samp.save_coeff();
    }


    void run_all_peaks(){
      //Peak finder structures
      ColumnVector p_theta;        //Spherical angles of points on the sphere used during peak finder
      ColumnVector p_phi;
      ColumnVector p_z;            //z coordinates of points on the sphere (used to separate points in two hemispheres)
      Matrix Eval_SH;              //Matrix that is used during peak finder and holds the value of spherical harmonics at the points on the sphere
      Matrix index;                //Matrix used during peak finder and stores the closest neighbours of each point on the sphere
      int peak_finder=1;           //Flag that keeps which peak finder is used
      
      Samples samp(m_mask, m_data.Ncols(),nsamples, npeaks, Num_of_coeff, m_savemeancoeff,m_gfa_flag);
      
      if (opts.peak_finder.value()==2){
	peak_finder=2;
	if (lmax!=4){
	  cout<<"Cannot use the continuous peak finder with lmax!=4. Switching to the discrete one!"<<endl;
	  peak_finder=1;
	}
      }

      if (peak_finder==1){
	if (opts.verbose.value())
	  cout<<"preparing Peak Finder structures"<<endl;
	Init_Peak_finder_finite_diff(index, Eval_SH, p_theta, p_phi, p_z);
	}

      if (opts.verbose.value())
	cout<<"running inference for each voxel"<<endl;

      for(int vox=1;vox<=m_data.Ncols();vox++){
	cout <<vox<<"/"<<m_data.Ncols()<<endl;
	Matrix SHcoeff; float gfa;
	if (m_modelnum==3){
	  Matrix temp_data(m_bvecs.Ncols(),3);  //3 is the number of q-shells
	  temp_data.Column(1)=(m_data.Column(vox)).Rows(1,m_bvecs.Ncols());  //Get multi-shell data for that voxel
	  temp_data.Column(2)=(m_data.Column(vox)).Rows(m_bvecs.Ncols()+1,2*m_bvecs.Ncols());
	  temp_data.Column(3)=(m_data.Column(vox)).Rows(2*m_bvecs.Ncols()+1,3*m_bvecs.Ncols());
	  ODF_voxel_biexp local_ODF(temp_data, m_SHT, m_HAT, m_scaling_factors, lmax, Num_of_coeff, nsamples, m_delta);
	  local_ODF.bootstrapping();
	  SHcoeff=local_ODF.get_SHcoeff();
	}
	else{
	  ColumnVector temp_data=m_data.Column(vox);  //data for that voxel
	  ODF_voxel local_ODF(temp_data, m_SHT, m_HAT, m_scaling_factors, m_modelnum, lmax, Num_of_coeff, nsamples,m_delta);
	  local_ODF.bootstrapping();
	  SHcoeff=local_ODF.get_SHcoeff();
	} 

	if (peak_finder==1){
	  Peak_finder_finite_diff Peaks(SHcoeff, npeaks, peak_threshold, lmax, nsamples, p_theta, p_phi, p_z, Eval_SH, index); 
	  Peaks.run();
	  Matrix theta=Peaks.get_theta();
	  Matrix phi=Peaks.get_phi();
	  Matrix funcvals=Peaks.get_funcvals();
	  samp.record(SHcoeff, theta, phi,funcvals, vox);
	}
	else{
	  Peak_finder_analytic Peaks(SHcoeff, npeaks, peak_threshold, lmax, nsamples); 
	  Peaks.run();
	  Matrix theta=Peaks.get_theta();
	  Matrix phi=Peaks.get_phi();
	  Matrix funcvals=Peaks.get_funcvals();
	  samp.record(SHcoeff,theta, phi,funcvals, vox);
	}

	if (m_gfa_flag){
	  ODF_GFA local_vox(SHcoeff, lmax, nsamples, Eval_SH);
	  local_vox.compute_gfa();
	  gfa=local_vox.get_gfa();
	  samp.record(gfa,vox);
	}
      }
      if (opts.verbose.value())
	cout<<"saving samples"<<endl;
      samp.save_dir();
    }

  };
 
}
#endif


