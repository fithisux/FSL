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


#if !defined(peakfinder_h)
#define peakfinder_h


#include <iostream>
#include "stdlib.h"
#include "libprob.h"
#include <cmath>
#include "miscmaths/miscmaths.h"

namespace ODFs{


  /////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////Useful Functions///////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////

  //Computes and returns a! (a factorial, a>=0)
  double factorial(int a){
    double res=1.0;
    int j;

    if (a==0)			//0!=1
      return res;
    else{
      for (j=1; j<=a; j++)
	res=res*j;
      return res;
    }
  }


  //Computes the associated Legendre polynomial Plm(x). Here m and l should be integers satisfying 0<=m<=l, while x lies in the range -1<=x<=1.
  //Function obtained from Numerical Recipes.
  double plgndr(int l, int m, double x){
    double fact,pll=0,pmm,pmmp1,somx2;
    int i,ll;
    
    //if (m < 0 || m > l || fabs(x) > 1.0){  cerr<<"Bad arguments for Legendre Polynomial"; return 0.0; }
    pmm=1.0; //Compute Pmm (l=m)
    if (m > 0){
      somx2=sqrt((1.0-x)*(1.0+x));
      fact=1.0;
      for (i=1;i<=m;i++){
	pmm *= -fact*somx2;
	fact += 2.0;
      }
    }
    if (l == m)
      return pmm;
    else{ //Compute P(m+1)m  (l=m+1)
      pmmp1=x*(2*m+1)*pmm;
      if (l == (m+1))
	return pmmp1;
      else{ //Compute Plm (l>m+1)
	for (ll=m+2;ll<=l;ll++){
	  pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
	  pmm=pmmp1;
	  pmmp1=pll;
	}
	return pll;
      }
    }
  } 


  //For a set of M points on the sphere (described by the spherical angles theta and phi)
  //fill a matrix with even spherical harmonics up to order max_order. The matrix is MxR, where R is the number of SH coefficients
  ReturnMatrix fill_SH_matrix(const ColumnVector& theta, const ColumnVector& phi, const int& max_order) {
    int R=(max_order+1)*(max_order+2)/2;
    int M=theta.Nrows();
    int j,l,m, p; float Plm,mag;
    Matrix tempSHT;
    tempSHT.ReSize(M,R);
      
    for (p=1; p<=M; p++){		//for each point (i.e. row of the matrix)
      j=1;
      for (l=0; l<=max_order; l=l+2)
	for (m=-l; m<=l; m++){
	  Plm=plgndr(l,abs(m),cos(theta(p)));
	  mag=sqrt((2*l+1)/(4.0*M_PI)*factorial(l-abs(m))/factorial(l+abs(m)))*Plm;
	  if (m<0)
	    tempSHT(p,j)=sqrt(2)*mag*cos(fabs(m)*phi(p));
	  else if (m==0)
	    tempSHT(p,j)=mag;
	  else
	    tempSHT(p,j)=pow(-1,(float)m)*sqrt(2)*mag*sin(m*phi(p));   //pow(-1,m+1.0) according to Maxime
	  j++;
	}
    }
    tempSHT.Release();
    return tempSHT;
  }



  /////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////
  //Class that estimates the peaks of Nsamples of an ODF, given Nsamples of the ODF SH coefficients, 
  //using analytic calculations (according to Aganj et al, MICCAI, 2010) 
  class Peak_finder_analytic{
    const Matrix& SHcoeff_samples;        //SH coefficients that describe the ODF (Num_coeff x N_samples)
    
                                          //The last entry in the following matrices will correspond to the mean ODF peaks
    Matrix m_theta_samples;               //Stores the theta angles of the extracted ODF peaks (maxNum_peaks x N_samples+1)
    Matrix m_phi_samples;                 //Stores the phi angles of the extracted ODF peaks (maxNum_peaks x N_samples+1)
    Matrix m_f_samples;                   //Stores the ODF amplitude at the extracted ODF peaks (maxNum_peaks x N_samples+1)
  
    const int& max_num_peaks;             //maximum number of peaks to be stored (if more detected, assume isotopic)
    const float& m_peak_threshold;        //Threshold on the minimum ODF value that a peak should have
    const int& lmax;                      //maximum SH order employed
    const int Num_of_coeff;               //number of SH coefficients 
    const int& nsamples;                  //number of bootstrap samples

    //temporary structures
    vector<ColumnVector> peaks_mx;        //Keeps the candidate peaks (theta, phi)

  public:
    //Constructor
  Peak_finder_analytic(const Matrix& coeff, const int& npeaks, const float& threshold, const int& SH_order, const int& nsamp): SHcoeff_samples(coeff), max_num_peaks(npeaks), m_peak_threshold(threshold), lmax(SH_order), Num_of_coeff((SH_order+1)*(SH_order+2)/2), nsamples(nsamp){
      m_theta_samples.ReSize(max_num_peaks, nsamples+1);  //The last sample will be the peaks of the mean ODF shape
      m_phi_samples.ReSize(max_num_peaks, nsamples+1);
      m_f_samples.ReSize(max_num_peaks, nsamples+1);
      m_theta_samples=0; m_phi_samples=0; m_f_samples=0;
    }
    
    //Destructor
    ~Peak_finder_analytic(){}


    Matrix& get_theta_ref()  {return m_theta_samples;}
    Matrix& get_phi_ref()   {return m_phi_samples;}
    Matrix& get_funcvals_ref()  {return m_f_samples;}
    Matrix get_theta()  {return m_theta_samples;}
    Matrix get_phi()   {return m_phi_samples;}
    Matrix get_funcvals()  {return m_f_samples;}


    void Find_candidate_peaks(const ColumnVector& SHcoeff) {
      const double thr=0.03;             //Threshold on the derivative of the ODF with respect to theta 
      const double phi_step=0.005;        //step size for 1D exhaustive search on phi
      bool HighRes;                      //When close to maxima increase resolution
      double mag, Y, Yp, sn, cs;
      double phi, dPhi;
      double A,B,C,D,E,F,G,H, Bp, Cp, Ep, Fp, Gp, Hp, Bs, Cs, Es, Fs, Gs, Hs;
      ColumnVector a,ap;
      a=SHcoeff; ap=SHcoeff;
      
      peaks_mx.clear();
      for (int constRes=0; constRes<=1; constRes++){
	phi=0;
	while (phi<(2*M_PI)){        //Iterate over phi
	  //For each phi, compute certain expressions
	 for (int l=0; l<=4; l=l+2){  
	    for (int m=-l; m<=l; m++){
	      int j=l*(l+1)/2+m+1;
	      if (m<0){
		mag=sqrt(((2*l+1)/(2*M_PI))*factorial(l+m)/factorial(l-m));
		Y=mag*cos(m*phi);
		Yp=-m*mag*sin(m*phi);
	      }
	      else if (m==0){
		Y=sqrt((2*l+1)/(4*M_PI));
		Yp=0;
	      }
	      else{
		mag=pow(-1,(float)m)*sqrt(((2*l+1)/(2*M_PI))*factorial(l-m)/factorial(l+m));
		Y=mag*sin(m*phi);
		Yp=m*mag*cos(m*phi);
	      }
	      a(j) =SHcoeff(j)*Y;
	      ap(j)=SHcoeff(j)*Yp;
	    }
	  }

	  //The following are functions only of phi
	  A=0.5*a(4);          B=-3*(a(3)+a(5));    C=3*(a(2)+a(6));  D=0.125*a(11); E=-2.5*(a(10)+a(12));  
	  F=7.5*(a(9)+a(13)); G=-105*(a(8)+a(14)); H=105*(a(7)+a(15));
        
	  Bp=-3*(ap(3)+ap(5));      Cp=3*(ap(2)+ap(6));      Ep=-2.5*(ap(10)+ap(12)); 
	  Fp=7.5*(ap(9)+ap(13)); Gp=-105*(ap(8)+ap(14));  Hp=105*(ap(7)+ap(15));
        
	  Bs=-B;    Cs=-4*C;  Es=-E; 
	  Fs=-4*F;  Gs=-9*G;  Hs=-16*H;

	  //Solve cubic for tan(theta)
	  vector<double> tth=solveCubic(Hp+Cp-Fp, Gp+Bp-3*Ep, 6*Fp+Cp, Bp+4*Ep); 
	 
	  HighRes = false;
	  dPhi = phi_step;
	  
	  //for each real cubic solution for tan(theta)
	  for (int n=0; n<(int)tth.size(); n++){   
	    double tmp=atan(tth[n]);
	    double th = floor(tmp/(float)M_PI);    //Get the value for theta
	    th=tmp-th*M_PI;                       //As the modulo of the division atan(tth[n])/M_PI
	    sn=sin(2*th); cs=cos(2*th);
	    tmp=ODF_dtheta(sn, cs, A, B, C, D, E, F, G, H);
	    if (fabs(tmp) < thr){  //If Eq. 12 is true
	      Matrix Hes(2,2);
	      Hes(1,1) = ODF_dtheta2(sn, cs, A, B, C, D, E, F, G, H);  //Compute the Hessian
	      Hes(1,2) = ODF_dtheta(sn, cs, 0, Bp, Cp, 0, Ep, Fp, Gp, Hp); Hes(2,1)=Hes(1,2);
	      Hes(2,2) = ODF_dphi2(sn, cs, 0, Bs, Cs, 0, Es, Fs, Gs, Hs);
	      double det=Hes.Determinant(); 
	      double tr=Hes.Trace();
	      HighRes=true;
	      if (det>=0 && tr<=0){   //Then we have a local maximum
		ColumnVector peak(2);
		peak<< th << phi;  //Store Results
		peaks_mx.push_back(peak);
	      }
	    }
	    if (constRes==1){
	      double t2=tth[n]*tth[n];  double t3=t2*tth[n]; double t4=t3*tth[n];
	      double const_step=phi_step*(1+t2)/sqrt(t2+t4+pow((((Hs+Cs-Fs)*t3+(Gs+Bs-3*Es)*t2+(6*Fs+Cs)*tth[n]+(Bs+4*Es))/(3*(Hp+Cp-Fp)*t2+2*(Gp+Bp-3*Ep)*tth[n]+(6*Fp+Cp))),2.0));
	      if (const_step<dPhi)
		dPhi=const_step;
	    }
	  }
	  //Update phi
	  if (HighRes)   
	    phi=phi+dPhi*0.5;
	  else
	    phi=phi+dPhi;
	}
      }
    }



    //Cluster candidate peaks (vectors close to each other within dist are clustered and averaged)
    //Then the ODF value is checked at a candidate peak, as it must be above threshold
    //If remaining peaks are still more than max_num_peaks, ODF is assumed isotropic and 1 peak is returned. 
    //A vector is returned with entries (theta, phi, ODF_value) for each peak that survived filtering
    vector<ColumnVector> Filter_peaks(const ColumnVector& SHcoeff) const {
      const double dist=0.4; 
      int npeaks=0, nM=0; double p, dp,dn,d;
      ColumnVector u(3), v_theta, v_phi;
      vector<ColumnVector> v;
      Matrix v_mx, Eval_SH; 
      vector < vector < ColumnVector > > clust;
      
      clust.resize((int)peaks_mx.size());
      for (int i=0; i<(int)peaks_mx.size(); i++){  //For each candidate peak
	u << sin(peaks_mx[i](1))*cos(peaks_mx[i](2)) << sin(peaks_mx[i](1))*sin(peaks_mx[i](2)) << cos(peaks_mx[i](1)); //get the corresponding vector u
	p=1e20;
	for (int n=0; n<npeaks; n++){  //for each other maximum v already visited
	  dp=sqrt((v[n](1)-u(1))*(v[n](1)-u(1))+(v[n](2)-u(2))*(v[n](2)-u(2))+(v[n](3)-u(3))*(v[n](3)-u(3))); //Get the euclidean distance
	  dn=sqrt((v[n](1)+u(1))*(v[n](1)+u(1))+(v[n](2)+u(2))*(v[n](2)+u(2))+(v[n](3)+u(3))*(v[n](3)+u(3))); //with u and -u
	  d=min(dp,dn);
	  if (d<p){           //Choose the closest v
	    p=d;     nM=n;    //store its index nM
	    if (dn<dp)
	      u=-u;           //and flip u if neccesary
	  }
	}
	if (p<dist)           //If u is very close to any other maximum v
	  clust[nM].push_back(u);  //store it with all other vectors that are close to v vector (with index nM)
	else {                //Otherwise store u as it is for output
	  v.push_back(u);
	  npeaks++;
	}
      }
      
      for (int i=0; i<(int)peaks_mx.size(); i++)   //for each cluster, store the mean vector 
	if (clust[i].size()!=0){
	  v[i]=0;
	  for (int vc=0; vc<(int)clust[i].size(); vc++)
	    v[i]=v[i]+clust[i][vc];
	  v[i]/=clust[i].size();
	  double mag=sqrt(v[i](1)*v[i](1)+v[i](2)*v[i](2)+v[i](3)*v[i](3));
	  v[i]/=mag;
	} 

      
      if (npeaks!=0){
	//Check the ODF amplitudes at each candidate peak
	Matrix v_mx(3,npeaks), Eval_SH;
	for (int i=0; i<npeaks; i++)
	  v_mx.Column(i+1)=v[i];
	cart2sph(v_mx,v_theta,v_phi);           //Convert candidate peaks to spherical angles
	Eval_SH=fill_SH_matrix(v_theta, v_phi, lmax);  //Evaluate spherical harmonics at each peak
	ColumnVector func_vals(npeaks);
	func_vals=0;
	for (int i=1; i<=npeaks; i++)	      //compute the ODF value at each peak
	  for (int j=1; j<=Num_of_coeff; j++)
	    func_vals(i)+=SHcoeff(j)*Eval_SH(i,j);
     
      
	v.clear();
	int pos;
	double maxval=func_vals.Maximum1(pos);  //Get the maximum ODF peak value and the corresponding peak index
	for (int i=1; i<=npeaks; i++)	      //Keep only peaks with high enough amplitude
	  if (func_vals(i)>=m_peak_threshold*maxval){
	    u<<v_theta(i)<<v_phi(i)<<func_vals(i);
	    v.push_back(u);
	  }
	npeaks=v.size(); 

	if (npeaks>max_num_peaks){             //If still too many peaks, keep only the max_num_peaks with maximum value
	  v.clear();
	  for (int i=1; i<=max_num_peaks; i++){
	    maxval=func_vals.Maximum1(pos);  //Get the maximum ODF peak value and the corresponding peak index
	    u<<v_theta(pos)<<v_phi(pos)<<func_vals(pos); 
	    v.push_back(u);
	    func_vals(pos)=0;               //zero that entry in order to find the next maximum    
	  }
	}
      }
      else{
	 u<<0<<0<<0; 
	 v.push_back(u);
      }   
      return v;
    }
    

    //Routine that finds and filters peaks and returns a (theta,phi,f) set for each peak
    void run(){
      for (int n=1; n<=nsamples; n++){
	Find_candidate_peaks(SHcoeff_samples.Column(n));
	vector<ColumnVector> peaks=Filter_peaks(SHcoeff_samples.Column(n));
	for (int i=0; i<(int)peaks.size(); i++){   //Save each returned peak
	  m_theta_samples(i+1,n)=peaks[i](1);
	  m_phi_samples(i+1,n)=peaks[i](2);
	  m_f_samples(i+1,n)=peaks[i](3);
	}
      }

      //Find the peaks of the average ODF shape across samples
      ColumnVector m_avg_coeff;
      m_avg_coeff.ReSize(Num_of_coeff);
      for (int i=1; i<=Num_of_coeff; i++)         //Get the average SH coefficients across samples 
	m_avg_coeff(i)=SHcoeff_samples.Row(i).Sum()/SHcoeff_samples.Ncols();
      
      Find_candidate_peaks(m_avg_coeff);
      vector<ColumnVector> peaks=Filter_peaks(m_avg_coeff);
      for (int i=0; i<(int)peaks.size(); i++){   //Save each returned peak
	m_theta_samples(i+1,nsamples+1)=peaks[i](1);
	m_phi_samples(i+1,nsamples+1)=peaks[i](2);
	m_f_samples(i+1,nsamples+1)=peaks[i](3);
      }
    }


    //Compute analytically elements of the Hessian of the ODF (for lmax=4). Hessian elements are returned as dphi2, dtheta_dphi and dtheta2.
    //Hessian is computed at a given (theta,phi) point. For computational speed, expressions of theta and phi encountered in the Hessian calculation are provided directly.
    //sn and cs are expressions of theta, A,B,C,D,E,F,G,H are expressions of phi.
    double ODF_dtheta2(const double& sn, const double& cs, const double& A, const double& B, const double& C, const double& D, const double& E, const double& F, const double& G, const double& H){
      double dtheta2=4*(G-7*E)*sn*cs + 2*(7*F-35*D-H)*(2*cs*cs-1) + 2*(H+C-F-3*A-5*D)*cs -(E+2*B+G)*sn;
      return dtheta2;
    }

    double ODF_dphi2(const double& sn, const double& cs, const double& A, const double& B, const double& C, const double& D, const double& E, const double& F, const double& G, const double& H){
      double dphi2=35*D*((1+cs)*(1+cs)/4)+(3*A-30*D)*(1+cs)/2.0+3*D-A + 0.5*(7*E*(1+cs)/2.0-3*E+B)*sn + (7*F*(1+cs)/2+C-F)*(1-cs)/2.0 + G*sn*(1-cs)/4.0 + H*((1-cs)*(1-cs)/4);
      return dphi2;
    }

    double ODF_dtheta(const double& sn, const double& cs, const double& A, const double& B, const double& C, const double& D, const double& E, const double& F, const double& G, const double& H){
      double dtheta=(G-7*E)*sn*sn + (7*F-35*D-H)*sn*cs + (H+C-F-3*A-5*D)*sn + (0.5*E+B+0.5*G)*cs -0.5*G+3.5*E;
      return dtheta;
    }

    //Cubic root
    double croot(const double& x){
      double res;
      
      if (x>=0) res=pow(x,1.0/3.0);
      else res=-pow(-x,1.0/3.0);
    
      return res;
    }


    //Solves analytically a polynomial equation up to degree 3, i.e. a*x^3+b*x^2+c*x+d=0 and returns
    //to a vector any real solutions
    vector<double> solveCubic(const double& a, const double& b, const double& c, const double& d){
      double p,q,b_a,c_a,d_a, discrim, offset,ee,tmp, root;
      vector<double> roots;
      double inv3=1.0/3.0;

      if (a!=0){    //Solve Cubic
	b_a=b/a; c_a=c/a; d_a=d/a;
	p=c_a-b_a*b_a*inv3;
	q=(2.0*b_a*b_a*b_a-9.0*b_a*c_a+27.0*d_a)/27.0;
	p=p*inv3;
	q=q*0.5;
	discrim=q*q+p*p*p;
	offset=b_a*inv3;
	if (discrim>0.0){  //One real root
	  ee=sqrt(discrim);
	  tmp=-q+ee;  root =croot(tmp);
	  tmp=-q-ee;  root+=croot(tmp);
	  root-=offset;
	  roots.push_back(root);
	}
	else if (discrim<0.0){  //Three real roots
	  ee=sqrt(-discrim);
	  double tmp2=-q; 
	  double angle= 2.0*inv3*atan(ee/(sqrt(tmp2*tmp2+ee*ee)+tmp2));
	  double sqrt3=sqrt(3.0);
	  tmp=cos(angle);
	  tmp2=sin(angle);
	  ee=sqrt(-p);
	  root=2*ee*tmp-offset; roots.push_back(root);
	  root=-ee*(tmp+sqrt3*tmp2)-offset; roots.push_back(root);
	  root=-ee*(tmp-sqrt3*tmp2)-offset; roots.push_back(root);
	}
	else{      //One or Two roots
	  tmp=-q;
	  tmp=croot(tmp);
	  root=2*tmp-offset; roots.push_back(root);
	  if (p!=0 || q!=0){
	    root=-tmp-offset; roots.push_back(root);
	  }
	}
      }
    
      else if (b!=0){  //Solve Quadratic
	discrim=c*c-4*b*d;
	if (discrim>0){
	  tmp=sqrt(discrim);
	  root=(-c+tmp)/(2.0*b); roots.push_back(root);
	  root=(-c-tmp)/(2.0*b); roots.push_back(root);
	}
	else if (discrim==0){
	  root=-c/(2.0*b); roots.push_back(root);
	}
      }
    
      else if (c!=0){  //Solve Linear
	root=-d/c; roots.push_back(root);
      }
      return roots;
    }

  };




  /////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////
  //Class that estimates the peaks of Nsamples of an ODF, given Nsamples of the ODF SH coefficients, 
  //using finite differences on the sphere. A sphere tesselation is used to
  //get points on the sphere and points that have values larger than their 
  //neighbours are candidate peaks. These are then filtered to give the final peaks
  class Peak_finder_finite_diff{
    const Matrix& SHcoeff_samples;        //SH coefficients that describe the ODF (Num_coeff x N_samples)
    
                                          //The last entry in the following matrices will correspond to the mean ODF peaks
    Matrix m_theta_samples;               //Stores the theta angles of the extracted ODF peaks (maxNum_peaks x N_samples+1)
    Matrix m_phi_samples;                 //Stores the phi angles of the extracted ODF peaks (maxNum_peaks x N_samples+1)
    Matrix m_f_samples;                   //Stores the ODF amplitude at the extracted ODF peaks (maxNum_peaks x N_samples+1)
  
    const int& max_num_peaks;             //maximum number of peaks to be stored (if more detected, assume isotopic)
    const float& m_peak_threshold;        //Threshold on the minimum ODF value that a peak should have
    const int& lmax;                      //maximum SH order employed
    const int Num_of_coeff;               //number of SH coefficients 
    const int& nsamples;                  //number of bootstrap samples
    
    //Structures common to all voxels
    const ColumnVector& p_theta;          //Theta and phi angles for the tesselation points 
    const ColumnVector& p_phi;
    const ColumnVector& p_z;              //z coordinate of the tesselation points (used to constrain seacrh in one hemisphere)
    const Matrix& Eval_SH;                //(i,j) element: for each tesselation point i, contains the jth spherical harmonic evaluated at i
    const Matrix& index;                  //(i,j) element: for each tesselation point i, contains the indices for the 6 closest points to i      

    //temporary structures
    vector<int> peaks_index;              //Keep an index to the peaks found
    vector<double> peaks_vals;             //Keeps the ODF value at the respective peaks

  public:
    //Constructor
  Peak_finder_finite_diff(const Matrix& coeff, const int& npeaks, const float& threshold, const int& SH_order, const int& nsamp, const ColumnVector& ptheta, const ColumnVector& pphi, const ColumnVector& pz, const Matrix& EvalSH, const Matrix& ind): SHcoeff_samples(coeff), max_num_peaks(npeaks), m_peak_threshold(threshold), lmax(SH_order), Num_of_coeff((SH_order+1)*(SH_order+2)/2), nsamples(nsamp), p_theta(ptheta),p_phi(pphi), p_z(pz), Eval_SH(EvalSH), index(ind) {
      m_theta_samples.ReSize(max_num_peaks, nsamples+1);  //The last sample will be the peaks of the mean ODF shape
      m_phi_samples.ReSize(max_num_peaks, nsamples+1);
      m_f_samples.ReSize(max_num_peaks, nsamples+1);
      m_theta_samples=0; m_phi_samples=0; m_f_samples=0;
    }
    
    //Destructor
    ~Peak_finder_finite_diff(){}

    Matrix& get_theta_ref()  {return m_theta_samples;}
    Matrix& get_phi_ref()   {return m_phi_samples;}
    Matrix& get_funcvals_ref()  {return m_f_samples;}
    Matrix get_theta()  {return m_theta_samples;}
    Matrix get_phi()   {return m_phi_samples;}
    Matrix get_funcvals()  {return m_f_samples;}



    //Finds local ODF maxima that represent fibre orientations 
    void Find_peaks(const ColumnVector& SHcoeff) {
      int num_neighbours=index.Ncols();
      int num_points=Eval_SH.Nrows();
      ColumnVector func_vals(num_points);
      peaks_index.clear();
      peaks_vals.clear();
      
      func_vals=0;
      for (int i=1; i<=num_points; i++)	      //compute the function value at each tesselation point
	for (int j=1; j<=Num_of_coeff; j++)
	  func_vals(i)+=SHcoeff(j)*Eval_SH(i,j);
  
      //Find the peaks using a finite difference method
      double maxFunc,maxFunc_prop; int flag;
      maxFunc=0;			     //Find the maximum of the function values
      for (int i=1; i<=num_points; i++)
	if (func_vals(i)>maxFunc)
	  maxFunc=func_vals(i); 
      maxFunc_prop=maxFunc*m_peak_threshold;

      //int pcount=1;					//Find the peaks
      for (int i=1; i<=num_points; i++){	//for each point
	if (p_z(i)>=0 && func_vals(i)>maxFunc_prop){      //search only in one hemisphere and for those points that have a value larger than a proportion of the Max_Func
	  flag=0; 
	  for (int j=1; j<=num_neighbours; j++)	//If the value at the focal point is smaller or equal than at least one neighbour, discard the point
	    if (func_vals(i)<=func_vals((int)index(i,j)))
	      flag=1;
	  if (flag==0){	                        //else the point is a peak
	    peaks_index.push_back(i);           //The stored peak index starts from 1! 
	    peaks_vals.push_back(func_vals(i)); //keep the function value at the peak
	  } 
	}
      } 
    }



    //Reduce threshold to get more filtering and separate peaks more
    vector<ColumnVector> Filter_peaks(const ColumnVector& SHcoeff){
      const double dist=0.4; 
      int npeaks=0, nM=0, num; double p, dp,dn,d;
      ColumnVector u(3), v_theta, v_phi;
      vector < vector < ColumnVector > > clust;
      vector<ColumnVector> v;
      
      num=peaks_index.size(); //Number of tentative peaks
      clust.resize(num);
      for (int i=1; i<=num; i++){  //For each candidate peak
	u << sin(p_theta(peaks_index[i-1]))*cos(p_phi(peaks_index[i-1])) << sin(p_theta(peaks_index[i-1]))*sin(p_phi(peaks_index[i-1])) << cos(p_theta(peaks_index[i-1])); //get the corresponding vector u
	p=1e20;
	for (int n=0; n<npeaks; n++){  //for each other maximum v already visited
	  dp=sqrt((v[n](1)-u(1))*(v[n](1)-u(1))+(v[n](2)-u(2))*(v[n](2)-u(2))+(v[n](3)-u(3))*(v[n](3)-u(3))); //Get the euclidean distance
	  dn=sqrt((v[n](1)+u(1))*(v[n](1)+u(1))+(v[n](2)+u(2))*(v[n](2)+u(2))+(v[n](3)+u(3))*(v[n](3)+u(3))); //with u and -u
	  d=min(dp,dn);
	  if (d<p){           //Choose the closest v
	    p=d;     nM=n;    //store its index nM
	    if (dn<dp)
	      u=-u;           //and flip u if neccesary
	  }
	}
	if (p<dist)          //If u is very close to any other maximum v
	  clust[nM].push_back(u);  //store it with all other vectors that are close to v vector (with index nM)
	else {                 //Otherwise store u as it is for output
	  v.push_back(u);
	  npeaks++;
	}
      }
      
      for (int i=0; i<num; i++)   //for each cluster, store the mean vector 
	if (clust[i].size()!=0){
	  v[i]=0;
	  for (int vc=0; vc<(int)clust[i].size(); vc++)
	    v[i]=v[i]+clust[i][vc];
	  v[i]/=clust[i].size();
	  double mag=sqrt(v[i](1)*v[i](1)+v[i](2)*v[i](2)+v[i](3)*v[i](3));
	  v[i]/=mag;
	} 
      
      if (npeaks!=0){
	//Check the ODF amplitudes at each candidate peak
	Matrix v_mx(3,npeaks), Eval_SH;
	for (int i=0; i<npeaks; i++)
	  v_mx.Column(i+1)=v[i];
	cart2sph(v_mx,v_theta,v_phi);           //Convert candidate peaks to spherical angles
	Eval_SH=fill_SH_matrix(v_theta, v_phi, lmax);  //Evaluate spherical harmonics at each peak
	ColumnVector func_vals(npeaks);
	func_vals=0;
	for (int i=1; i<=npeaks; i++)	      //compute the ODF value at each peak
	  for (int j=1; j<=Num_of_coeff; j++)
	    func_vals(i)+=SHcoeff(j)*Eval_SH(i,j);
     
      
	v.clear();
	int pos;
	double maxval=func_vals.Maximum1(pos);  //Get the maximum ODF peak value and the corresponding peak index
	for (int i=1; i<=npeaks; i++)	      //Keep only peaks with high enough amplitude
	  if (func_vals(i)>=m_peak_threshold*maxval){
	    u<<v_theta(i)<<v_phi(i)<<func_vals(i);
	    v.push_back(u);
	  }
	npeaks=v.size(); 

	if (npeaks>max_num_peaks){             //If still too many peaks, keep only the max_num_peaks with maximum value
	  v.clear();
	  for (int i=1; i<=max_num_peaks; i++){
	    maxval=func_vals.Maximum1(pos);  //Get the maximum ODF peak value and the corresponding peak index
	    u<<v_theta(pos)<<v_phi(pos)<<func_vals(pos); 
	    v.push_back(u);
	    func_vals(pos)=0;               //zero that entry in order to find the next maximum    
	  }
	}
      }
      else{
	 u<<0<<0<<0; 
	 v.push_back(u);
      }   
      return v;
    }



    //Reduce threshold to get more filtering and separate peaks more
    vector<ColumnVector> Filter_peaks_old(const float min_angle_separation=0.8){
      int total_flag,maxpos,num;
      ColumnVector counter;
      Matrix temp_peaks;
      double dot,max;
      vector<ColumnVector> v;

      num=peaks_index.size(); //Number of tentative peaks
      
      temp_peaks.ReSize(num,3);
      counter.ReSize(num);
      
      for (int i=1; i<=num; i++){  //Get the Cartesian coordinates of the peaks
	temp_peaks(i,1)=sin(p_theta(peaks_index[i-1]))*cos(p_phi(peaks_index[i-1]));
	temp_peaks(i,2)=sin(p_theta(peaks_index[i-1]))*sin(p_phi(peaks_index[i-1]));
	temp_peaks(i,3)=cos(p_theta(peaks_index[i-1]));
      }
    	  
      total_flag=1;
      while (total_flag==1){
	total_flag=0;
	for (int i=1; i<=num; i++){  //for each peak 
	  counter(i)=0;
	  for (int j=1; j<=num; j++){   //compute the dot product with all the rest
	    dot=fabs(temp_peaks(i,1)*temp_peaks(j,1)+temp_peaks(i,2)*temp_peaks(j,2)+temp_peaks(i,3)*temp_peaks(j,3));
	    if (j!=i && dot>min_angle_separation){  //if the angle is smaller than the threshold, i.e. the dot product is larger than the threshold
	      counter(i)=counter(i)+1;
	      total_flag=1;
	    }
	  }
	}
	if (total_flag==1){   //if a vector to be filtered exists
	  max=0;
	  maxpos=0;
	  for (int i=1; i<=num; i++){	//Find the vector with the maximum counter value
	    if (counter(i)>max){
	      max=counter(i);
	      maxpos=i;
	    }
	    else{                //if there are many vectors with the same counter value, remove the one with the smallest function value
	      if (counter(i)==max && peaks_vals[i-1]<=peaks_vals[maxpos-1]){
		max=counter(i);
		maxpos=i;
	      }
	    }
	  }
	  temp_peaks(maxpos,1)=temp_peaks(num,1); //Erase it by taking the last entry of peaks, put it in its place and reduce
	  temp_peaks(maxpos,2)=temp_peaks(num,2); //the size of peaks by 1.
	  temp_peaks(maxpos,3)=temp_peaks(num,3);
	  peaks_vals[maxpos-1]=peaks_vals[num-1];
	  peaks_index[maxpos-1]=peaks_index[num-1];
	  num=num-1;
	}
      }
      
      max=0;
      maxpos=0;
      if (num>max_num_peaks){		  //If too many peaks still exist, then assume the ODF is isotropic and keep only the peak
	for (int i=0; i<(int)peaks_index.size(); i++)	  //with the maximum func_value
	  if (peaks_vals[i]>max){
	    max=peaks_vals[i];
	    maxpos=i;
	  }
	num=1;
	peaks_vals[0]=peaks_vals[maxpos];
	peaks_index[0]=peaks_index[maxpos];
      } 
     
      return v;
    }


    //Routine that finds and filters peaks and returns a (theta,phi,f) set for each peak
    void run(){
      for (int n=1; n<=nsamples; n++){
	Find_peaks(SHcoeff_samples.Column(n));
	vector<ColumnVector> v=Filter_peaks(SHcoeff_samples.Column(n));
	for (int i=0; i<(int)v.size(); i++){   //Save each returned peak
	  m_theta_samples(i+1,n)=v[i](1);
	  m_phi_samples(i+1,n)=v[i](2);
	  m_f_samples(i+1,n)=v[i](3);

	}
      }

      //Find the peaks of the average ODF shape across samples
      ColumnVector m_avg_coeff;
      m_avg_coeff.ReSize(Num_of_coeff);
      for (int i=1; i<=Num_of_coeff; i++)         //Get the average SH coefficients across samples 
	m_avg_coeff(i)=SHcoeff_samples.Row(i).Sum()/SHcoeff_samples.Ncols();
      
      Find_peaks(m_avg_coeff);
      vector<ColumnVector> v=Filter_peaks(m_avg_coeff);
      for (int i=0; i<(int)v.size(); i++){   //Save each returned peak
	m_theta_samples(i+1,nsamples+1)=v[i](1);
	m_phi_samples(i+1,nsamples+1)=v[i](2);
	m_f_samples(i+1,nsamples+1)=v[i](3);
      }
    }


   }; 


}
#endif
