/*  Diffusion model fitting

    Timothy Behrens, Saad Jbabdi, Stam Sotiropoulos  - FMRIB Image Analysis Group
 
    Copyright (C) 2005 University of Oxford  */

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

#if !defined (diffmodels_h)
#define diffmodels_h

#include <iostream>
#include <fstream>
#include <iomanip>
#define WANT_STREAM
#define WANT_MATH
#include <string>
#include "utils/log.h"
#include "utils/tracer_plus.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/nonlin.h"
#include "stdlib.h"



using namespace NEWMAT;
using namespace MISCMATHS;


#define two_pi 0.636619772
#define FSMALL 0.001

#define f2x(x) (std::tan((x)/two_pi))   //fraction transformation used in the old model 1
#define x2f(x) (std::abs(two_pi*std::atan((x))))

#define f2beta(f) (std::asin(std::sqrt(f))) //fraction transformation used in the new model 1
#define beta2f(beta) (std::pow(std::sin(beta),2.0))
#define d2lambda(d) (std::sqrt(d))     //diffusivity transformation used in the new model 1
#define lambda2d(lambda) (lambda*lambda)   //d->lambda^2>=0

#define lowlim 4.0                        //lowlim>0
//#define k12l1(k1) (std::sqrt(k1-lowlim))  //transformation used in the fanning model for the Bingham principal eigenvalue k1
//#define l12k1(l1) (l1*l1+lowlim)          //k1->l1^2+lowlim>=lowlim>0 
#define UL 10000                        //lowlim>0
#define k12l1(k1) (std::asin(std::sqrt((k1-lowlim)/UL)))    
#define l12k1(l1) (pow(sin(l1),2.0)*UL+lowlim)        
                     
#define upperlim 100.0                     //upperlim>1
//#define Invupperlim 0.034482758620690    //1/(upperlim-1)
#define w2gam(w) (std::asin(std::sqrt((w-1)/(upperlim-1))))       //transformation used in the fanning model for the Bingham eigenvalue ratio k1/k2
#define gam2w(gam) (1.0+std::pow(std::sin(gam),2.0)*(upperlim-1))  //w->1+(upperlim-1)*sin^2(gam), so that 1<=w<=upperlim


#define bigger(a,b) ((a)>(b)?(a):(b))
#define smaller(a,b) ((a)>(b)?(b):(a))

#define nonzerosign(a) ((a)<0?-1:1)
#define tiny 1.0e-5                             //Used for numerical diffenetiation
#define SQRTtiny (sqrt(tiny))



////////////////////////////////////////////////
//       DIFFUSION TENSOR MODEL
////////////////////////////////////////////////

class DTI : public NonlinCF{
public: 
  DTI(const ColumnVector& iY,
      const Matrix& ibvecs,const Matrix& ibvals){
    Y = iY;
    npts = Y.Nrows();
    m_v1.ReSize(3);
    m_v2.ReSize(3);
    m_v3.ReSize(3);
    bvecs=ibvecs;
    bvals=ibvals;
    form_Amat();
    nparams=7;
  }
  DTI(const ColumnVector& iY,
      const Matrix& inAmat):Amat(inAmat){
    Y = iY;
    npts = Y.Nrows();
    m_v1.ReSize(3);
    m_v2.ReSize(3);
    m_v3.ReSize(3);
    nparams=7;
    iAmat = pinv(Amat);
  }
  ~DTI(){}
  void linfit();
  void nonlinfit();
  void calc_tensor_parameters();
  void sort();
  void set_data(const ColumnVector& data){Y=data;}
  float get_fa()const{return m_fa;}
  float get_md()const{return m_md;}
  float get_s0()const{return m_s0;}
  float get_mo()const{return m_mo;}
  ColumnVector get_v1()const{return m_v1;}
  ColumnVector get_v2()const{return m_v2;}
  ColumnVector get_v3()const{return m_v3;}
  float get_l1()const{return m_l1;}
  float get_l2()const{return m_l2;}
  float get_l3()const{return m_l3;}
  ColumnVector get_eigen()const{ColumnVector x(3);x<<m_l1<<m_l2<<m_l3;return x;}
  ColumnVector get_tensor()const{
    ColumnVector x(6);
    x << m_tens(1,1)
      << m_tens(2,1)
      << m_tens(3,1)
      << m_tens(2,2)
      << m_tens(3,2)
      << m_tens(3,3);
    return x;
  }
  ColumnVector get_v(const int& i)const{if(i==1)return m_v1;else if(i==2)return m_v2;else return m_v3;}
  ReturnMatrix get_prediction()const;
  SymmetricMatrix get_covar()const{return m_covar;}
  ColumnVector get_data()const{return Y;}
  Matrix get_Amat()const{return Amat;}

  // derivatives of tensor functions w.r.t. tensor parameters
  ReturnMatrix calc_fa_grad(const ColumnVector& _tens)const;
  float calc_fa_var()const;
  ColumnVector calc_md_grad(const ColumnVector& _tens)const;
  ColumnVector calc_mo_grad(const ColumnVector& _tens)const;

  // conversion between rotation matrix and angles
  void rot2angles(const Matrix& rot,float& th1,float& th2,float& th3)const;
  void angles2rot(const float& th1,const float& th2,const float& th3,Matrix& rot)const;

  void print()const{
    cout << "DTI FIT RESULTS " << endl;
    cout << "S0   :" << m_s0 << endl;
    cout << "MD   :" << m_md << endl;
    cout << "FA   :" << m_fa << endl;
    cout << "MO   :" << m_mo << endl;
    ColumnVector x(3);
    x=m_v1;
    if(x(3)<0)x=-x;
    float _th,_ph;cart2sph(x,_th,_ph);
    cout << "TH   :" << _th*180.0/M_PI << " deg" << endl; 
    cout << "PH   :" << _ph*180.0/M_PI << " deg" << endl; 
    cout << "V1   : " << x(1) << " " << x(2) << " " << x(3) << endl;
  }
  void form_Amat(){
    Amat.ReSize(bvecs.Ncols(),7);
    Matrix tmpvec(3,1), tmpmat;
    for( int i = 1; i <= bvecs.Ncols(); i++){
      tmpvec << bvecs(1,i) << bvecs(2,i) << bvecs(3,i);
      tmpmat = tmpvec*tmpvec.t()*bvals(1,i);
      Amat(i,1) = tmpmat(1,1);
      Amat(i,2) = 2*tmpmat(1,2);
      Amat(i,3) = 2*tmpmat(1,3);
      Amat(i,4) = tmpmat(2,2);
      Amat(i,5) = 2*tmpmat(2,3);
      Amat(i,6) = tmpmat(3,3);
      Amat(i,7) = 1;
    }
    iAmat = pinv(Amat);
  }
  void vec2tens(const ColumnVector& Vec){
    m_tens.ReSize(3);
    m_tens(1,1)=Vec(1);
    m_tens(2,1)=Vec(2);
    m_tens(3,1)=Vec(3);
    m_tens(2,2)=Vec(4);
    m_tens(3,2)=Vec(5);
    m_tens(3,3)=Vec(6);
  }
  void vec2tens(const ColumnVector& Vec,SymmetricMatrix& Tens)const{
    Tens.ReSize(3);
    Tens(1,1)=Vec(1);
    Tens(2,1)=Vec(2);
    Tens(3,1)=Vec(3);
    Tens(2,2)=Vec(4);
    Tens(3,2)=Vec(5);
    Tens(3,3)=Vec(6);
  }
  void tens2vec(const SymmetricMatrix& Tens,ColumnVector& Vec)const{
    Vec.ReSize(6);
    Vec<<Tens(1,1)<<Tens(2,1)<<Tens(3,1)<<Tens(2,2)<<Tens(3,2)<<Tens(3,3);
  }

  // nonlinear fitting routines
  NEWMAT::ReturnMatrix grad(const NEWMAT::ColumnVector& p)const;
  boost::shared_ptr<BFMatrix> hess(const NEWMAT::ColumnVector&p,boost::shared_ptr<BFMatrix> iptr)const;
  double cf(const NEWMAT::ColumnVector& p)const;
  NEWMAT::ReturnMatrix forwardModel(const NEWMAT::ColumnVector& p)const;
  
  ColumnVector rotproduct(const ColumnVector& x,const Matrix& R)const;
  ColumnVector rotproduct(const ColumnVector& x,const Matrix& R1,const Matrix& R2)const;
  float anisoterm(const int& pt,const ColumnVector& ls,const Matrix& xx)const;
  
private:
  Matrix bvecs;
  Matrix bvals;
  ColumnVector Y;
  Matrix Amat,iAmat;
  int npts,nparams;
  ColumnVector m_v1,m_v2,m_v3;
  float m_l1,m_l2,m_l3;
  float m_fa,m_s0,m_md,m_mo;
  float m_sse;
  SymmetricMatrix m_tens;
  SymmetricMatrix m_covar;
};


////////////////////////////////////////////////
//       Partial Volume Models
////////////////////////////////////////////////

// Generic class
class PVM {
public:
  PVM(const ColumnVector& iY,
      const Matrix& ibvecs, const Matrix& ibvals,
      const int& nfibres):Y(iY),bvecs(ibvecs),bvals(ibvals){
    
    npts    = Y.Nrows();
    nfib    = nfibres;
    
    cart2sph(ibvecs,alpha,beta);
    
    cosalpha.ReSize(npts);
    sinalpha.ReSize(npts);
    for(int i=1;i<=npts;i++){
      sinalpha(i) = sin(alpha(i));
      cosalpha(i) = cos(alpha(i));
    }
    
  }
  virtual ~PVM(){}
  
  // PVM virtual routines
  virtual void fit()  = 0;
  virtual void sort() = 0;
  virtual void print()const = 0;
  virtual void print(const ColumnVector& p)const = 0;
  
  virtual ReturnMatrix get_prediction()const = 0;

protected:
  const ColumnVector& Y;
  const Matrix& bvecs;
  const Matrix& bvals;
  ColumnVector alpha;
  ColumnVector sinalpha;
  ColumnVector cosalpha;
  ColumnVector beta;
  
  int npts;
  int nfib;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Model 1 : mono-exponential (for single shell). Constrained optimization for the diffusivity, fractions and their sum<1
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class PVM_single_c : public PVM, public NonlinCF {
public:
   PVM_single_c(const ColumnVector& iY,
	     const Matrix& ibvecs, const Matrix& ibvals,
		const int& nfibres, bool m_BIC=false, bool incl_f0=false, bool m_fan_angle=false):PVM(iY,ibvecs,ibvals,nfibres),m_include_f0(incl_f0),m_eval_BIC(m_BIC),m_return_fanning(m_fan_angle){

    if (m_include_f0)
      nparams = nfib*3 + 3; 
    else
      nparams = nfib*3 + 2;
    
    m_f.ReSize(nfib);
    m_th.ReSize(nfib);
    m_ph.ReSize(nfib);
    if (m_return_fanning)
      m_fanning_angles.ReSize(nfib);
  }
  ~PVM_single_c(){}

  // routines from NonlinCF
  NEWMAT::ReturnMatrix grad(const NEWMAT::ColumnVector& p)const;
  boost::shared_ptr<BFMatrix> hess(const NEWMAT::ColumnVector&p,boost::shared_ptr<BFMatrix> iptr)const;
  double cf(const NEWMAT::ColumnVector& p)const;
  NEWMAT::ReturnMatrix forwardModel(const NEWMAT::ColumnVector& p)const;
 

  // other routines
  void fit();
  void sort();                                         //Sort compartments according to their volume fraction
  void fit_pvf(ColumnVector& x)const;                  //Estimate the volume fractions given all the other parameters using Linear Least Squares. Used to better initialize the Nonlinear fitter
  void fix_fsum(ColumnVector& fs) const;
  float partial_fsum(ColumnVector& fs, int ii) const;  //Returns 1-Sum(f_j), 1<=j<=ii. (ii<=nfib). Used for transforming beta to f and vice versa
  void print()const;                                   //Print the final estimates (after having them transformed)
  void print(const ColumnVector& p)const;              //Print the estimates using a vector with the untransformed parameter values 
  ReturnMatrix get_prediction()const;                  //Applies the forward model and gets the model predicted signal using the estimated parameter values (true,non-transformed space)  
  
  float get_s0()const{return m_s0;}
  float get_f0()const{return m_f0;}
  float get_d()const{return m_d;}
  ColumnVector get_f()const{return m_f;}
  ColumnVector get_th()const{return m_th;}
  ColumnVector get_ph()const{return m_ph;}
  float get_f(const int& i)const{return m_f(i);}
  float get_th(const int& i)const{return m_th(i);}
  float get_ph(const int& i)const{return m_ph(i);}
  float get_BIC() const{return m_BIC;}
  ColumnVector get_fanning_angles() const{return m_fanning_angles;}
  float get_fanning_angle(const int& i) const{return m_fanning_angles(i);}
  vector<ColumnVector> get_invHes_e1() const{return m_invprHes_e1;}
  ColumnVector get_invHes_e1(const int& i) const{return m_invprHes_e1[i-1];}
  vector<Matrix> get_Hessian() const{return m_Hessian;}
  Matrix get_Hessian(const int& i) const{return m_Hessian[i-1];}


  //Functions used to obtain a prediction of the fanning angle (if any), associated with each fibre compartment
  void eval_Hessian_at_peaks();        //For each fibre, compute a 3x3 Hessian of the cartesian (x,y,z) coordinates of the orientation
                                       //evaluated at the estimated parameters
  void Fanning_angles_from_Hessian();  //For each fibre, get the projection of the 2nd eigenvector of the Hessian to the fanning plane and get a fanning angle in [0,pi).
 

  // useful functions for calculating signal and its derivatives
  // functions
  float isoterm(const int& pt,const float& _d)const;
  float anisoterm(const int& pt,const float& _d,const ColumnVector& x)const;
  // 1st order derivatives
  float isoterm_lambda(const int& pt,const float& lambda)const;
  float anisoterm_lambda(const int& pt,const float& lambda,const ColumnVector& x)const;
  float anisoterm_th(const int& pt,const float& _d,const ColumnVector& x,const float& _th,const float& _ph)const;
  float anisoterm_ph(const int& pt,const float& _d,const ColumnVector& x,const float& _th,const float& _ph)const;
  ReturnMatrix fractions_deriv(const int& nfib, const ColumnVector& fs, const ColumnVector& bs) const;
  
private:  
  int   nparams;
  float m_s0;
  float m_d;
  float m_f0;
  ColumnVector m_f;
  ColumnVector m_th;
  ColumnVector m_ph;
  const bool m_include_f0;       //Indicate whether f0 will be used in the model (an unattenuated signal compartment). That will be added as the last parameter
  const bool m_eval_BIC;         //Indicate whether the Bayesian Information Criterion for the fitted model is computed
  const bool m_return_fanning;   //Indicate whether fanning angles predictions are made. For each fitted fibre compartment i, use the second eigenvector of the inverse Hessian 
                                 //evaluated at this fibre orientation to predict fanning angle for i.  
  float m_BIC;                   //Bayesian Information Criterion for the fitted model
  ColumnVector m_fanning_angles; //Use the second eigenvector of the inverse Hessian evaluated at each fibre orientation i to predict fanning angle for fibre compartment i.  
  vector<Matrix> m_Hessian;      //Vector that keeps the Hessian matrix for each fibre orientation w.r.t. the Cartesian coordinates x,y,z, evaluated at the estimated orientation
  vector<ColumnVector> m_invprHes_e1;  //Vector that keeps the first eigenvector of the projected inverse Hessian for each fibre orientation w.r.t. the Cartesian coordinates x,y,z, evaluated at the estimated orientation
};



///////////////////////////////////////////////////////////////////////////
//       Old Model 1 with no constraints for the sum of fractions
//////////////////////////////////////////////////////////////////////////
// Model 1 : mono-exponential (for single shell)
class PVM_single : public PVM, public NonlinCF {
public:
  PVM_single(const ColumnVector& iY,
	     const Matrix& ibvecs, const Matrix& ibvals,
	     const int& nfibres, bool incl_f0=false):PVM(iY,ibvecs,ibvals,nfibres), m_include_f0(incl_f0){

    if (m_include_f0)
      nparams = nfib*3 + 3; 
    else
      nparams = nfib*3 + 2;

    m_f.ReSize(nfib);
    m_th.ReSize(nfib);
    m_ph.ReSize(nfib);
  }
  ~PVM_single(){}

  // routines from NonlinCF
  NEWMAT::ReturnMatrix grad(const NEWMAT::ColumnVector& p)const;
  boost::shared_ptr<BFMatrix> hess(const NEWMAT::ColumnVector&p,boost::shared_ptr<BFMatrix> iptr)const;
  double cf(const NEWMAT::ColumnVector& p)const;
  NEWMAT::ReturnMatrix forwardModel(const NEWMAT::ColumnVector& p)const;

  // other routines
  void fit();
  void sort();
  void fix_fsum();
  void print()const{
    cout << "PVM (Single) FIT RESULTS " << endl;
    cout << "S0   :" << m_s0 << endl;
    cout << "D    :" << m_d << endl;
    for(int i=1;i<=nfib;i++){
      cout << "F" << i << "   :" << m_f(i) << endl;
      ColumnVector x(3);
      x << sin(m_th(i))*cos(m_ph(i)) << sin(m_th(i))*sin(m_ph(i)) << cos(m_th(i));
      if(x(3)<0)x=-x;
      float _th,_ph;cart2sph(x,_th,_ph);
      cout << "TH" << i << "  :" << _th*180.0/M_PI << " deg" << endl; 
      cout << "PH" << i << "  :" << _ph*180.0/M_PI << " deg" << endl; 
      cout << "DIR" << i << "   : " << x(1) << " " << x(2) << " " << x(3) << endl;
    }
  }

  void print(const ColumnVector& p)const{
    cout << "PARAMETER VALUES " << endl;
    cout << "S0   :" << p(1) << endl;
    cout << "D    :" << p(2) << endl;
    for(int i=3,ii=1;ii<=nfib;i+=3,ii++){
      cout << "F" << ii << "   :" << x2f(p(i)) << endl;
      cout << "TH" << ii << "  :" << p(i+1)*180.0/M_PI << " deg" << endl; 
      cout << "PH" << ii << "  :" << p(i+2)*180.0/M_PI << " deg" << endl; 
    }
    if (m_include_f0)
      cout << "f0    :" << x2f(p(nparams)) << endl;
  }

  // getters
  float get_s0()const{return m_s0;}
  float get_f0()const{return m_f0;}
  float get_d()const{return m_d;}
  ColumnVector get_f()const{return m_f;}
  ColumnVector get_th()const{return m_th;}
  ColumnVector get_ph()const{return m_ph;}
  float get_f(const int& i)const{return m_f(i);}
  float get_th(const int& i)const{return m_th(i);}
  float get_ph(const int& i)const{return m_ph(i);}
  ReturnMatrix get_prediction()const;

  // setters	 
  void set_s0(const float& s0){m_s0=s0;}	 
  void set_f0(const float& f0){m_f0=f0;}	 
  void set_d(const float& d){m_d=d;}	 
  void set_f(const ColumnVector& f){m_f=f;}	 
  void set_th_ph(const Matrix& dyads){	 
    MISCMATHS::cart2sph(dyads,m_th,m_ph);	 
  }

  // useful functions for calculating signal and its derivatives
  // functions
  float isoterm(const int& pt,const float& _d)const;
  float anisoterm(const int& pt,const float& _d,const ColumnVector& x)const;
  float bvecs_fibre_dp(const int& pt,const float& _th,const float& _ph)const;
  float bvecs_fibre_dp(const int& pt,const ColumnVector& x)const;
  // 1st order derivatives
  float isoterm_d(const int& pt,const float& _d)const;
  float anisoterm_d(const int& pt,const float& _d,const ColumnVector& x)const;
  float anisoterm_th(const int& pt,const float& _d,const ColumnVector& x,const float& _th,const float& _ph)const;
  float anisoterm_ph(const int& pt,const float& _d,const ColumnVector& x,const float& _th,const float& _ph)const;
  // 2nd order derivatives
  float isoterm_dd(const int& pt,const float& _d)const;
  float anisoterm_dd(const int& pt,const float& _d,const ColumnVector& x)const;
  float anisoterm_dth(const int& pt,const float& _d,const ColumnVector& x,const float& _th,const float& _ph)const;
  float anisoterm_dph(const int& pt,const float& _d,const ColumnVector& x,const float& _th,const float& _ph)const;
  float anisoterm_thth(const int& pt,const float& _d,const ColumnVector& x,const float& _th,const float& _ph)const;
  float anisoterm_phph(const int& pt,const float& _d,const ColumnVector& x,const float& _th,const float& _ph)const;
  float anisoterm_thph(const int& pt,const float& _d,const ColumnVector& x,const float& _th,const float& _ph)const;

private:  
  int   nparams;
  float m_s0;
  float m_d;
  float m_f0;
  ColumnVector m_f;
  ColumnVector m_th;
  ColumnVector m_ph;
  const bool m_include_f0;   //Indicate whether f0 will be used in the model (an unattenuated signal compartment)
};



////////////////////////////////////////////////
//       Partial Volume Models
////////////////////////////////////////////////

// Model 2 : non-mono-exponential (for multiple shells)
class PVM_multi : public PVM, public NonlinCF {
public:
  PVM_multi(const ColumnVector& iY,
	    const Matrix& ibvecs, const Matrix& ibvals,
	    const int& nfibres, int Gamma_for_ball_only=0, float R=0.13, bool incl_f0=false):PVM(iY,ibvecs,ibvals,nfibres),m_Gamma_for_ball_only(Gamma_for_ball_only),m_R(R),m_include_f0(incl_f0){

    if (m_include_f0)
      nparams = nfib*3 + 4; 
    else    
      nparams = nfib*3 + 3;

    m_invR=1.0/(2.0*m_R+1);
    m_f.ReSize(nfib);
    m_th.ReSize(nfib);
    m_ph.ReSize(nfib);
  }
  ~PVM_multi(){}

  // routines from NonlinCF
  NEWMAT::ReturnMatrix grad(const NEWMAT::ColumnVector& p)const;
  boost::shared_ptr<BFMatrix> hess(const NEWMAT::ColumnVector&p,boost::shared_ptr<BFMatrix> iptr)const;
  double cf(const NEWMAT::ColumnVector& p)const;
  NEWMAT::ReturnMatrix forwardModel(const NEWMAT::ColumnVector& p)const;

  // other routines
  void fit();
  void sort();
  void fix_fsum();
  void print()const{
    cout << "PVM (MULTI) FIT RESULTS " << endl;
    cout << "S0    :" << m_s0 << endl;
    cout << "D     :" << m_d << endl;
    cout << "D_STD :" << m_d_std << endl;
    for(int i=1;i<=nfib;i++){
      cout << "F" << i << "    :" << m_f(i) << endl;
      ColumnVector x(3);
      x << sin(m_th(i))*cos(m_ph(i)) << sin(m_th(i))*sin(m_ph(i)) << cos(m_th(i));
      if(x(3)<0)x=-x;
      cout << "TH" << i << "   :" << m_th(i) << endl; 
      cout << "PH" << i << "   :" << m_ph(i) << endl; 
      cout << "DIR" << i << "   : " << x(1) << " " << x(2) << " " << x(3) << endl;
    }
  }
  void print(const ColumnVector& p)const{
    cout << "PARAMETER VALUES " << endl;
    cout << "S0    :" << p(1) << endl;
    cout << "D     :" << p(2) << endl;
    cout << "D_STD :" << p(3) << endl;
    for(int i=3,ii=1;ii<=nfib;i+=3,ii++){
      cout << "F" << ii << "    :" << x2f(p(i)) << endl;
      cout << "TH" << ii << "   :" << p(i+1) << endl; 
      cout << "PH" << ii << "   :" << p(i+2) << endl; 
    }
    if (m_include_f0)
      cout << "f0    :" << x2f(p(nparams)) << endl;
  }

  float get_s0()const{return m_s0;}
  float get_d()const{return m_d;}
  float get_f0()const{return m_f0;}
  float get_d_std()const{return m_d_std;}
  ColumnVector get_f()const{return m_f;}
  ColumnVector get_th()const{return m_th;}
  ColumnVector get_ph()const{return m_ph;}
  float get_f(const int& i)const{return m_f(i);}
  float get_th(const int& i)const{return m_th(i);}
  float get_ph(const int& i)const{return m_ph(i);}

  // setters	 
  void set_s0(const float& s0){m_s0=s0;}	 
  void set_f0(const float& f0){m_f0=f0;}	 
  void set_d(const float& d){m_d=d;}	 
  void set_d_std(const float& d_std){m_d_std=d_std;}	 
  void set_f(const ColumnVector& f){m_f=f;}	 
  void set_th_ph(const Matrix& dyads){	 
    MISCMATHS::cart2sph(dyads,m_th,m_ph);	 
  }

  ReturnMatrix get_prediction()const;

  // useful functions for calculating signal and its derivatives
  // functions
  float isoterm(const int& pt,const float& _a,const float& _b)const;
  float anisoterm(const int& pt,const float& _a,const float& _b,const ColumnVector& x, const int)const;
  // 1st order derivatives
  float isoterm_a(const int& pt,const float& _a,const float& _b)const;
  float anisoterm_a(const int& pt,const float& _a,const float& _b,const ColumnVector& x,const int)const;
  float isoterm_b(const int& pt,const float& _a,const float& _b)const;
  float anisoterm_b(const int& pt,const float& _a,const float& _b,const ColumnVector& x,const int)const;
  float anisoterm_th(const int& pt,const float& _a,const float& _b,const ColumnVector& x,const float& _th,const float& _ph,const int)const;
  float anisoterm_ph(const int& pt,const float& _a,const float& _b,const ColumnVector& x,const float& _th,const float& _ph,const int)const;
  
private:
  int   nparams;
  float m_s0;
  float m_d;
  float m_d_std;
  float m_f0;
  ColumnVector m_f;
  ColumnVector m_th;
  ColumnVector m_ph;
  const int m_Gamma_for_ball_only; //Model2 Twists: 0 (default), 1 sets off the Gamma distr of diff. for the fibres, 2 as in (1) but also use fixed tensor kernels for fibres
  const float m_R;                 //If m_Gamma_for_ball_only=2, then m_R defines the anisotropy of the tensor kernels. This is kept constant, used to initialise xfibres model=3.
  float m_invR;                    //If m_Gamma_for_ball_only=2, 1/(2*m_R+1) is precomputed as used a lot in calculations
  const bool m_include_f0;         //Indicate whether f0 will be used in the model (an unattenuated signal compartment)
};




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Ball & Binghams Fanning Model : Pseudo-Constrained optimization for the diffusivity, fractions and their sum<1, the eigenvalues
//of the Bingham Matrices. Use Levenberg-Marquardt to fit and get the gradient numerically 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class PVM_Ball_Binghams : public PVM, public NonlinCF {
public:
   PVM_Ball_Binghams(const ColumnVector& iY,
	     const Matrix& ibvecs, const Matrix& ibvals,
		     const int& nfibres, bool eval_BIC=false, bool incl_f0=false, bool grid_search=false):PVM(iY,ibvecs,ibvals,nfibres),m_eval_BIC(eval_BIC),m_include_f0(incl_f0), m_gridsearch(grid_search){

    nparams_per_fibre=6;

    if (m_include_f0)
      nparams = nfib*nparams_per_fibre + 3; 
    else
      nparams = nfib*nparams_per_fibre + 2;
    
    Matrix temp; ColumnVector gvec(3);
    //For each DW direction contains the dyadic product Matrix scaled by the b value: bvals(i)*bvecs(i)*bvecs(i)^T 
    for(int i=1;i<=npts;i++){
      gvec<< ibvecs(1,i) << ibvecs(2,i) << ibvecs(3,i);
      temp<< gvec*gvec.t();
      temp=ibvals(1,i)*temp;
      bvecs_dyadic.push_back(temp); 
    }

    m_f.ReSize(nfib);
    m_th.ReSize(nfib);
    m_ph.ReSize(nfib);
    m_psi.ReSize(nfib);
    m_k1.ReSize(nfib);
    m_k2.ReSize(nfib);
    m_BIC=1e15;
  }

  ~PVM_Ball_Binghams(){}

  // routines from NonlinCF
  NEWMAT::ReturnMatrix grad(const NEWMAT::ColumnVector& p) const;
  boost::shared_ptr<BFMatrix> hess(const NEWMAT::ColumnVector&p,boost::shared_ptr<BFMatrix> iptr)const;
  double cf(const NEWMAT::ColumnVector& p)const;
  NEWMAT::ReturnMatrix forwardModel(const NEWMAT::ColumnVector& p)const; //Applies the forward model and gets a model predicted signal using the parameter values in p (transformed parameter space)                                                                               
  //Instead of returning the model predicted signal for each direction returns the individual signal contributions weighted by their fractions
  //i.e. isotropic, anisotropic1, anisotropic2, etc. Summing those gives the signal
  NEWMAT::ReturnMatrix forwardModel_compartments(const NEWMAT::ColumnVector& p)const;

  //Builds up the model predicted signal for each direction by using precomputed individual compartment signals, stored in Matrix Sig. 
  //Weights them with the fractions, scales with S0 and sums to get the signal.
  NEWMAT::ReturnMatrix pred_from_compartments(const NEWMAT::ColumnVector& p, const NEWMAT::Matrix& Sig) const;
  
  //Builds up the model predicted signal for each direction by using precomputed individual compartment signals, stored in Matrix Sig. 
  //Weights them with the fractions, scales with S0 and sums to get the signal.
  //The signal of the fibre compartment with index fib is recalculated.
  NEWMAT::ReturnMatrix pred_from_compartments(const NEWMAT::ColumnVector& p, const NEWMAT::Matrix& Sig,const int& fib) const;


  // other routines
  void fit();
  void sort();                                         //Sort compartments according to their volume fraction
  float partial_fsum(ColumnVector& fs, int ii) const;  //Returns 1-Sum(f_j), 1<=j<=ii. (ii<=nfib). Used for transforming beta to f and vice versa
  void print()const;                                   //Print the final estimates (after having them untransformed)
  void print(const ColumnVector& p)const;              //Print the estimates using a vector with the transformed parameter values (i.e. need to untransform to get d,fs etc)
  ReturnMatrix get_prediction()const;                  //Applies the forward model and gets the model predicted signal using the estimated parameter values  (true,non-transformed space)  
  
  float get_s0()const{return m_s0;}
  float get_f0()const{return m_f0;}
  float get_d()const{return m_d;}
  ColumnVector get_f()const{return m_f;}
  ColumnVector get_th()const{return m_th;}
  ColumnVector get_ph()const{return m_ph;}
  ColumnVector get_psi()const{return m_psi;}
  ColumnVector get_k1()const{return m_k1;}
  ColumnVector get_k2()const{return m_k2;}
  float get_f(const int& i)const{return m_f(i);}
  float get_th(const int& i)const{return m_th(i);}
  float get_ph(const int& i)const{return m_ph(i);}
  float get_psi(const int& i)const{return m_psi(i);}
  ReturnMatrix get_fanning_vector(const int& i) const;  //Returns a vector that indicates the fanning orientation
  float get_k1(const int& i)const{return m_k1(i);}
  float get_k2(const int& i)const{return m_k2(i);}
  float get_BIC() const{return m_BIC;}
  
 

  // useful functions for calculating signal and its derivatives
  float isoterm(const int& pt,const float& _d)const;
  ReturnMatrix fractions_deriv(const int& nfib, const ColumnVector& fs, const ColumnVector& bs) const;

private:  
  vector<Matrix> bvecs_dyadic;   //For each DW direction contains the dyadic product Matrix scaled by the b value: bvals(i)*bvecs(i)*bvecs(i)^T 
  int nparams_per_fibre;
  int   nparams;
  float m_s0;
  float m_d;
  float m_f0;
  ColumnVector m_f;
  ColumnVector m_th;
  ColumnVector m_ph;
  ColumnVector m_psi;
  ColumnVector m_k1;
  ColumnVector m_k2;

  float m_BIC;                   //Bayesian Information Criterion for the fitted model

  const bool m_eval_BIC;         //Indicate whether the Bayesian Information Criterion for the fitted model is computed
  const bool m_include_f0;       //Indicate whether f0 will be used in the model (an unattenuated signal compartment). That will be added as the last parameter
  const bool m_gridsearch;       //Indicate whether a grid search will be used to initialize the fanning eigenvalues. This can take much longer to run. 
};




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Ball & Watsons Fanning Model : Pseudo-Constrained optimization for the diffusivity, fractions and their sum<1, the spread of the 
//Watson distribution. Use Levenberg-Marquardt to fit and get the gradient numerically 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class PVM_Ball_Watsons : public PVM, public NonlinCF {
public:
   PVM_Ball_Watsons(const ColumnVector& iY,
	     const Matrix& ibvecs, const Matrix& ibvals,
		    const int& nfibres, bool eval_BIC=false, bool incl_f0=false, bool grid_search=false):PVM(iY,ibvecs,ibvals,nfibres),m_eval_BIC(eval_BIC),m_include_f0(incl_f0), m_gridsearch(grid_search){

    nparams_per_fibre=4;

    if (m_include_f0)
      nparams = nfib*nparams_per_fibre + 3; 
    else
      nparams = nfib*nparams_per_fibre + 2;
    
    m_f.ReSize(nfib);
    m_th.ReSize(nfib);
    m_ph.ReSize(nfib);
    m_k.ReSize(nfib);
    m_BIC=1e15;
  }

  ~PVM_Ball_Watsons(){}

  // routines from NonlinCF
  NEWMAT::ReturnMatrix grad(const NEWMAT::ColumnVector& p) const;
  boost::shared_ptr<BFMatrix> hess(const NEWMAT::ColumnVector&p,boost::shared_ptr<BFMatrix> iptr)const;
  double cf(const NEWMAT::ColumnVector& p)const;
  NEWMAT::ReturnMatrix forwardModel(const NEWMAT::ColumnVector& p)const; //Applies the forward model and gets a model predicted signal using the parameter values in p (transformed parameter space)                                                                               
  //Instead of returning the model predicted signal for each direction returns the individual signal contributions weighted by their fractions
  //i.e. isotropic, anisotropic1, anisotropic2, etc. Summing those gives the signal
  NEWMAT::ReturnMatrix forwardModel_compartments(const NEWMAT::ColumnVector& p)const;

  //Builds up the model predicted signal for each direction by using precomputed individual compartment signals, stored in Matrix Sig. 
  //Weights them with the fractions, scales with S0 and sums to get the signal.
  NEWMAT::ReturnMatrix pred_from_compartments(const NEWMAT::ColumnVector& p, const NEWMAT::Matrix& Sig) const;
  
  //Builds up the model predicted signal for each direction by using precomputed individual compartment signals, stored in Matrix Sig. 
  //Weights them with the fractions, scales with S0 and sums to get the signal.
  //The signal of the fibre compartment with index fib is recalculated.
  NEWMAT::ReturnMatrix pred_from_compartments(const NEWMAT::ColumnVector& p, const NEWMAT::Matrix& Sig,const int& fib) const;


  // other routines
  void fit();
  void sort();                                         //Sort compartments according to their volume fraction
  float partial_fsum(ColumnVector& fs, int ii) const;  //Returns 1-Sum(f_j), 1<=j<=ii. (ii<=nfib). Used for transforming beta to f and vice versa
  void print()const;                                   //Print the final estimates (after having them untransformed)
  void print(const ColumnVector& p)const;              //Print the estimates using a vector with the transformed parameter values (i.e. need to untransform to get d,fs etc)
  ReturnMatrix get_prediction()const;                  //Applies the forward model and gets the model predicted signal using the estimated parameter values  (true,non-transformed space)  
  
  float get_s0()const{return m_s0;}
  float get_f0()const{return m_f0;}
  float get_d()const{return m_d;}
  ColumnVector get_f()const{return m_f;}
  ColumnVector get_th()const{return m_th;}
  ColumnVector get_ph()const{return m_ph;}
  ColumnVector get_k()const{return m_k;}
  float get_f(const int& i)const{return m_f(i);}
  float get_th(const int& i)const{return m_th(i);}
  float get_ph(const int& i)const{return m_ph(i);}
  float get_k(const int& i)const{return m_k(i);}
  float get_BIC() const{return m_BIC;}
  
  // useful functions for calculating signal and its derivatives
  float isoterm(const int& pt,const float& _d)const;
  ReturnMatrix fractions_deriv(const int& nfib, const ColumnVector& fs, const ColumnVector& bs) const;

private:  
  int nparams_per_fibre;
  int   nparams;
  float m_s0;
  float m_d;
  float m_f0;
  ColumnVector m_f;
  ColumnVector m_th;
  ColumnVector m_ph;
  ColumnVector m_k;

  float m_BIC;                   //Bayesian Information Criterion for the fitted model

  const bool m_eval_BIC;         //Indicate whether the Bayesian Information Criterion for the fitted model is computed
  const bool m_include_f0;       //Indicate whether f0 will be used in the model (an unattenuated signal compartment). That will be added as the last parameter
  const bool m_gridsearch;       //Indicate whether a grid search will be used to initialize the fanning eigenvalues. This can take much longer to run. 

};



#endif
