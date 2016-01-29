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

#include "Bingham_Watson_approx.h"
#include "diffmodels.h"



////////////////////////////////////////////////
//       DIFFUSION TENSOR MODEL
////////////////////////////////////////////////
void DTI::linfit(){
  ColumnVector logS(npts);
  ColumnVector Dvec(7);
  for (int i=1;i<=npts; i++){
    if(Y(i)>0)
      logS(i)=log(Y(i));
    else
      logS(i)=0;
  }
  Dvec = -iAmat*logS;
  if(Dvec(7)>-23)
    m_s0=exp(-Dvec(7));
  else
    m_s0=Y.MaximumAbsoluteValue();
  for (int i=1;i<=Y.Nrows();i++){
    if(m_s0<Y.Sum()/Y.Nrows()){ m_s0=Y.MaximumAbsoluteValue();  }
    logS(i)=(Y(i)/m_s0)>0.01 ? log(Y(i)):log(0.01*m_s0);
  }
  Dvec = -iAmat*logS;
  m_sse = (Amat*Dvec+logS).SumSquare();
  m_s0=exp(-Dvec(7));
  if(m_s0<Y.Sum()/Y.Nrows()){ m_s0=Y.Sum()/Y.Nrows();  }
  vec2tens(Dvec);
  calc_tensor_parameters();

  m_covar.ReSize(7);
  float dof=logS.Nrows()-7;
  float sig2=m_sse/dof;
  m_covar << sig2*(Amat.t()*Amat).i();
}
ColumnVector DTI::calc_md_grad(const ColumnVector& _tens)const{
  ColumnVector g(6);
  g = 0;
  g(1) = 1/3.0;
  g(4) = 1/3.0;
  g(6) = 1/3.0;
  return g;
}
// this will only work if the determinant is strictly positive
ReturnMatrix DTI::calc_fa_grad(const ColumnVector& _intens)const{
  ColumnVector gradv(6),ik(6),k(6);
  float m = (_intens(1)+_intens(4)+_intens(6))/3.0;
  SymmetricMatrix K(3),iK(3),M(3);

  // rescale input matrix
  vec2tens(_intens,M);
  //M /=m;

  m = M.Trace()/3.0;

  K = M - m*IdentityMatrix(3);
  tens2vec(K,k);
  iK << K.i();
  tens2vec(iK,ik);

  float p   = K.SumSquare()/6.0;
  float q   = K.Determinant()/2.0;
  float h   = std::sqrt(p*p*p-q*q)/q;
  float phi = std::atan(h)/3.0;
  if(q<0)phi+=M_PI;


  float _l1 = m + 2.0*std::sqrt(p)*std::cos(phi);
  float _l2 = m - std::sqrt(p)*(std::cos(phi)+std::sqrt(3.0)*std::sin(phi));
  float _l3 = m - std::sqrt(p)*(std::cos(phi)-std::sqrt(3.0)*std::sin(phi));

  float t  = 6.0/9.0*(_l1*_l1+_l2*_l2+_l3*_l3 - _l1*_l2-_l1*_l3-_l2*_l3);
  float b  = _l1*_l1+_l2*_l2+_l3*_l3;

  float _fa  = std::sqrt(3.0/2.0)*std::sqrt(t/b);
  

  float dfadl1 = 3.0/4.0/_fa * ( 6.0/9.0*(2.0*_l1-_l2-_l3)/b - t/b/b*2.0*_l1 );
  float dfadl2 = 3.0/4.0/_fa * ( 6.0/9.0*(2.0*_l2-_l1-_l3)/b - t/b/b*2.0*_l2 );
  float dfadl3 = 3.0/4.0/_fa * ( 6.0/9.0*(2.0*_l3-_l1-_l2)/b - t/b/b*2.0*_l3 );

  

  // determine dkdx
  ColumnVector dkdx(6);
  dkdx << 2.0/3.0 << 1.0 << 1.0 << 2.0/3.0 << 1.0 << 2.0/3.0;

  for(int i=1;i<=6;i++){
    float dL1dx=0,dL2dx=0,dL3dx=0;
    if(i==1||i==4||i==6){
      dL1dx=1.0/3.0;dL2dx=1.0/3.0;dL3dx=1.0/3.0;
    }
    //
    float p_p = k(i)/3.0 * dkdx(i);
    float q_p = q*ik(i) * dkdx(i);
    float h_p = (3.0*p*p*p_p - 2.0*q_p*q*(1+h*h))/2.0/h/q/q;

    float phi_p = h_p/(1+h*h)/3.0;

    dL1dx += p_p/std::sqrt(p)*std::cos(phi) - 2.0*std::sqrt(p)*phi_p*std::sin(phi);
    dL2dx -= std::sqrt(p)*(.5*p_p*(m-_l2)+phi_p*(std::sin(phi)-std::sqrt(3.0)*std::cos(phi)));
    dL3dx -= std::sqrt(p)*(.5*p_p*(m-_l3)+phi_p*(std::sin(phi)+std::sqrt(3.0)*std::cos(phi)));

    //
    gradv(i) = dfadl1*dL1dx + dfadl2*dL2dx + dfadl3*dL3dx;
  }
  gradv.Release();
  return gradv;
}
float DTI::calc_fa_var()const{
  ColumnVector grd;
  ColumnVector vtens;
  tens2vec(m_tens,vtens);
  grd = calc_fa_grad(vtens);
  ColumnVector g(7);
  g.SubMatrix(1,6,1,1) = grd;
  g(7) = 0;
  
  return((g.t()*m_covar*g).AsScalar());
}

void DTI::rot2angles(const Matrix& rot,float& th1,float& th2,float& th3)const{
  if(rot(3,1)!=1 && rot(3,1)!=-1){
    th2 = -asin(rot(3,1));
    float c=std::cos(th2);
    th1 = atan2(rot(3,2)/c,rot(3,3)/c);
    th3 = atan2(rot(2,1)/c,rot(1,1)/c);
  }
  else{
    th1 = atan2(rot(1,2),rot(1,3));
    th2 = -sign(rot(3,1))*M_PI/2;
    th3 = 0;
  }
}

void DTI::angles2rot(const float& th1,const float& th2,const float& th3,Matrix& rot)const{
  float c1=std::cos(th1),s1=std::sin(th1);
  float c2=std::cos(th2),s2=std::sin(th2);
  float c3=std::cos(th3),s3=std::sin(th3);

  rot(1,1) = c2*c3;    rot(1,2) = s1*s2*c3-c1*s3;    rot(3,1) = c1*s2*c3+s1*s3;
  rot(2,1) = c2*s3;    rot(2,2) = s1*s2*s3+c1*c3;    rot(3,2) = c1*s2*s3-s1*c3;
  rot(3,1) = -s2;      rot(3,2) = s1*c2;             rot(3,3) = c1*c2;
}


// nonlinear tensor fitting 
void DTI::nonlinfit(){
  // initialise with linear fitting
  linfit();

  print();

  // set starting parameters
  // params = s0, log(l1),log(l2), log(l3), th1, th2, th3
  ColumnVector start(nparams);

  start(1) = m_s0;
  // eigenvalues
  start(2) = m_l1>0?std::log(m_l1):std::log(1e-5);
  start(3) = m_l2>0?std::log(m_l2):std::log(1e-5);
  start(4) = m_l3>0?std::log(m_l3):std::log(1e-5);
  // angles
  float th1,th2,th3;
  Matrix rot(3,3);
  rot.Row(1) = m_v1.t();
  rot.Row(2) = m_v2.t();
  rot.Row(3) = m_v3.t();
  rot2angles(rot,th1,th2,th3);
  start(5) = th1;
  start(6) = th2;
  start(7) = th3;


  // do the fit
  NonlinParam  lmpar(start.Nrows(),NL_LM); 
  lmpar.SetGaussNewtonType(LM_L);
  lmpar.SetStartingEstimate(start);


  NonlinOut status;
  status = nonlin(lmpar,(*this));
  ColumnVector final_par(nparams);
  final_par = lmpar.Par();


  // finalise parameters
  m_s0 = final_par(1);
  m_l1 = exp(final_par(2));
  m_l2 = exp(final_par(3));
  m_l3 = exp(final_par(4));

  angles2rot(final_par(5),final_par(6),final_par(7),rot);
  m_v1 = rot.Row(1).t();
  m_v2 = rot.Row(2).t();
  m_v3 = rot.Row(3).t();

  sort();

  m_tens << m_l1*m_v1*m_v1.t() + m_l2*m_v2*m_v2.t() + m_l3*m_v3*m_v3.t();
  calc_tensor_parameters();

  print();
  //exit(1);

}
void DTI::sort(){
  vector< pair<float,int> > ls(3);
  vector<ColumnVector> vs(3);
  ls[0].first=m_l1;
  ls[0].second=0;
  ls[1].first=m_l2;
  ls[1].second=1;
  ls[2].first=m_l3;
  ls[2].second=2;
  vs[0]=m_v1;vs[1]=m_v2;vs[2]=m_v3;
  
  std::sort(ls.begin(),ls.end());

  m_l1 = ls[2].first;
  m_v1 = vs[ ls[2].second ];
  m_l2 = ls[1].first;
  m_v2 = vs[ ls[1].second ];
  m_l3 = ls[0].first;
  m_v3 = vs[ ls[0].second ];
  
}
void DTI::calc_tensor_parameters(){
  Matrix Vd;
  DiagonalMatrix Dd(3);
  // mean, eigenvalues and eigenvectors
  EigenValues(m_tens,Dd,Vd);
  m_md = Dd.Sum()/Dd.Nrows();
  m_l1 = Dd(3,3);
  m_l2 = Dd(2,2);
  m_l3 = Dd(1,1);
  m_v1 = Vd.Column(3);
  m_v2 = Vd.Column(2);
  m_v3 = Vd.Column(1);
  // mode
  float e1=m_l1-m_md, e2=m_l2-m_md, e3=m_l3-m_md;
  float n = (e1 + e2 - 2*e3)*(2*e1 - e2 - e3)*(e1 - 2*e2 + e3);
  float d = (e1*e1 + e2*e2 + e3*e3 - e1*e2 - e2*e3 - e1*e3);
  d = sqrt(bigger(0, d));
  d = 2*d*d*d;
  m_mo = smaller(bigger(d ? n/d : 0.0, -1),1);
  // fa
  float numer=1.5*((m_l1-m_md)*(m_l1-m_md)+(m_l2-m_md)*(m_l2-m_md)+(m_l3-m_md)*(m_l3-m_md));
  float denom=(m_l1*m_l1+m_l2*m_l2+m_l3*m_l3);
  if(denom>0) m_fa=numer/denom;
  else m_fa=0;
  if(m_fa>0){m_fa=sqrt(m_fa);}
  else{m_fa=0;}
}
// now the nonlinear fitting routines
double DTI::cf(const NEWMAT::ColumnVector& p)const{
  //cout << "CF" << endl;
  //OUT(p.t());
  double cfv = 0.0;
  double err = 0.0;
  ////////////////////////////////////
  ColumnVector ls(3);
  Matrix rot(3,3);
  angles2rot(p(5),p(6),p(7),rot);
  for(int k=2;k<=4;k++){
    ls(k-1) = exp(p(k));
  }
  ////////////////////////////////////
  for(int i=1;i<=Y.Nrows();i++){
    err = p(1)*anisoterm(i,ls,rot) - Y(i); 
    cfv += err*err; 
  }  
  //OUT(cfv);
  //cout<<"--------"<<endl;
  return(cfv);
}

NEWMAT::ReturnMatrix DTI::grad(const NEWMAT::ColumnVector& p)const{
  NEWMAT::ColumnVector gradv(p.Nrows());

  cout<<"grad"<<endl;
  OUT(p.t());

  gradv = 0.0;
  ////////////////////////////////////
  ColumnVector ls(3);
  Matrix rot(3,3);
  Matrix rot1(3,3),rot2(3,3),rot3(3,3);
  angles2rot(p(5),p(6),p(7),rot);

  angles2rot(p(5)+M_PI/2.0,p(6),p(7),rot1);
  angles2rot(p(5),p(6)+M_PI/2.0,p(7),rot2);
  angles2rot(p(5),p(6),p(7)+M_PI/2.0,rot3);
  for(int k=2;k<=4;k++){
    ls(k-1) = exp(p(k));
  }
  ////////////////////////////////////
  Matrix J(npts,nparams);
  ColumnVector x(3);
  ColumnVector diff(npts);
  float sig;
  for(int i=1;i<=Y.Nrows();i++){
    sig = p(1)*anisoterm(i,ls,rot);
    
    J(i,1) = sig/p(1);

    x = rotproduct(bvecs.Column(i),rot);
    J(i,2) = -bvals(1,i)*x(1)*sig*ls(1);
    J(i,3) = -bvals(1,i)*x(2)*sig*ls(2);
    J(i,4) = -bvals(1,i)*x(3)*sig*ls(3);

    x = rotproduct(bvecs.Column(i),rot1,rot);
    J(i,5) = -2.0*bvals(1,i)*(ls(1)*x(1)+ls(2)*x(2)+ls(3)*x(3))*sig;
    x = rotproduct(bvecs.Column(i),rot2,rot);
    J(i,6) = -2.0*bvals(1,i)*(ls(1)*x(1)+ls(2)*x(2)+ls(3)*x(3))*sig;
    x = rotproduct(bvecs.Column(i),rot3,rot);
    J(i,7) = -2.0*bvals(1,i)*(ls(1)*x(1)+ls(2)*x(2)+ls(3)*x(3))*sig;

    diff(i) = sig - Y(i);
  }

  OUT(diff.t());
  OUT(J.t());
  
  gradv = 2.0*J.t()*diff;

  OUT(gradv.t());
  cout<<"------"<<endl;

  gradv.Release();
  return gradv;
}

//this uses Gauss-Newton approximation
boost::shared_ptr<BFMatrix> DTI::hess(const NEWMAT::ColumnVector& p,boost::shared_ptr<BFMatrix> iptr)const{
  boost::shared_ptr<BFMatrix>   hessm;
  if (iptr && iptr->Nrows()==(unsigned int)p.Nrows() && iptr->Ncols()==(unsigned int)p.Nrows()) hessm = iptr;
  else hessm = boost::shared_ptr<BFMatrix>(new FullBFMatrix(p.Nrows(),p.Nrows()));

  cout<<"hess"<<endl;
  OUT(p.t());
  
  ////////////////////////////////////
  ColumnVector ls(3);
  Matrix rot(3,3);
  Matrix rot1(3,3),rot2(3,3),rot3(3,3);
  angles2rot(p(5),p(6),p(7),rot);

  angles2rot(p(5)+M_PI/2,p(6),p(7),rot1);
  angles2rot(p(5),p(6)+M_PI/2,p(7),rot2);
  angles2rot(p(5),p(6),p(7)+M_PI/2,rot3);
  for(int k=2;k<=4;k++){
    ls(k-1) = exp(p(k));
  }
  ////////////////////////////////////
  Matrix J(npts,nparams);
  ColumnVector x(3);
  float sig;
  for(int i=1;i<=Y.Nrows();i++){
    sig = p(1)*anisoterm(i,ls,rot);
    
    J(i,1) = sig/p(1);

    x = rotproduct(bvecs.Column(i),rot);
    J(i,2) = -bvals(1,i)*x(1)*sig*ls(1);
    J(i,3) = -bvals(1,i)*x(2)*sig*ls(2);
    J(i,4) = -bvals(1,i)*x(3)*sig*ls(3);

    x = rotproduct(bvecs.Column(i),rot1,rot);
    J(i,5) = -2.0*bvals(1,i)*(ls(1)*x(1)+ls(2)*x(2)+ls(3)*x(3))*sig;
    x = rotproduct(bvecs.Column(i),rot2,rot);
    J(i,6) = -2.0*bvals(1,i)*(ls(1)*x(1)+ls(2)*x(2)+ls(3)*x(3))*sig;
    x = rotproduct(bvecs.Column(i),rot3,rot);
    J(i,7) = -2.0*bvals(1,i)*(ls(1)*x(1)+ls(2)*x(2)+ls(3)*x(3))*sig;
  }
  

  for (int i=1; i<=p.Nrows(); i++){
    for (int j=i; j<=p.Nrows(); j++){
      sig = 0.0;
      for(int k=1;k<=J.Nrows();k++)
	sig += 2.0*(J(k,i)*J(k,j));
      hessm->Set(i,j,sig);
    }
  }
  for (int j=1; j<=p.Nrows(); j++) {
    for (int i=j+1; i<=p.Nrows(); i++) {
      hessm->Set(i,j,hessm->Peek(j,i));
    }
  }

  hessm->Print();
  cout<<"------"<<endl;

  return(hessm);
}

ColumnVector DTI::rotproduct(const ColumnVector& x,const Matrix& R)const{
  ColumnVector ret(3);
  
  for(int i=1;i<=3;i++)
    ret(i) = x(1)*x(1)*R(1,i)*R(1,i)+x(2)*x(2)*R(2,i)*R(2,i)+x(3)*x(3)*R(3,i)*R(3,i)
      +2.0*( x(1)*R(1,i)*(x(2)*R(2,i)+x(3)*R(3,i)) +x(2)*x(3)*R(2,i)*R(3,i) );   
  
  return ret;
}
ColumnVector DTI::rotproduct(const ColumnVector& x,const Matrix& R1,const Matrix& R2)const{
  ColumnVector ret(3);
  
  for(int i=1;i<=3;i++)
    ret(i) = x(1)*x(1)*R1(1,i)*R2(1,i)+x(2)*x(2)*R1(2,i)*R2(2,i)+x(3)*x(3)*R1(3,i)*R2(3,i)
      +( x(1)*R1(1,i)*(x(2)*R2(2,i)+x(3)*R2(3,i)) +x(2)*x(3)*R1(2,i)*R2(3,i) )
      +( x(1)*R2(1,i)*(x(2)*R1(2,i)+x(3)*R1(3,i)) +x(2)*x(3)*R2(2,i)*R1(3,i) );   
  
  return ret;
}

float DTI::anisoterm(const int& pt,const ColumnVector& ls,const Matrix& rot)const{
  ColumnVector x(3);
  x = rotproduct(bvecs.Column(pt),rot);

  return exp(-bvals(1,pt)*(ls(1)*x(1)+ls(2)*x(2)+ls(3)*x(3)));
}




/////////////////////////////////////////////////////////////////////////
//       PARTIAL VOLUME MODEL - SINGLE SHELL 
// Constrained Optimization for the diffusivity, fractions and their sum<1
//////////////////////////////////////////////////////////////////////////

void PVM_single_c::fit(){

  // initialise with a tensor
  DTI dti(Y,bvecs,bvals);
  dti.linfit();

  // set starting parameters for nonlinear fitting
  float _th,_ph;
  cart2sph(dti.get_v1(),_th,_ph);

  ColumnVector start(nparams);
  //Initialize the non-linear fitter. Use the DTI estimates for most parameters, apart from the volume fractions
  start(1) = dti.get_s0();
  //start(2) = d2lambda(dti.get_md()>0?dti.get_md()*2:0.001); // empirically found that d~2*MD
  start(2) = d2lambda(dti.get_l1()>0?dti.get_l1():0.002); // empirically found that d~L1
  start(4) = _th;
  start(5) = _ph;
  for(int ii=2,i=6;ii<=nfib;ii++,i+=3){
    cart2sph(dti.get_v(ii),_th,_ph);
    start(i+1) = _th;
    start(i+2) = _ph;
  }
  
  // do a better job for initializing the volume fractions
  fit_pvf(start);

  // do the fit
  NonlinParam  lmpar(start.Nrows(),NL_LM); 
  lmpar.SetGaussNewtonType(LM_L);
  lmpar.SetStartingEstimate(start);

  NonlinOut status;
  status = nonlin(lmpar,(*this));
  ColumnVector final_par(nparams);
  final_par = lmpar.Par();
  
  if (m_eval_BIC){ 
    double RSS=cf(final_par); //get the sum of squared residuals
    m_BIC=npts*log(RSS/npts)+log(npts)*nparams; //evaluate BIC
  }

  // finalise parameters
  m_s0 = final_par(1);
  m_d  = lambda2d(final_par(2));
  for(int k=1;k<=nfib;k++){
    int kk = 3 + 3*(k-1);
    m_f(k)  = beta2f(final_par(kk))*partial_fsum(m_f,k-1);
    m_th(k) = final_par(kk+1);
    m_ph(k) = final_par(kk+2);
  }
  
  if (m_return_fanning)
    Fanning_angles_from_Hessian();
  
  if (m_include_f0)
    m_f0=beta2f(final_par(nparams))*partial_fsum(m_f,nfib);
  
  sort(); 
}

void PVM_single_c::sort(){
  vector< pair<float,int> > fvals(nfib);
  ColumnVector ftmp(nfib),thtmp(nfib),phtmp(nfib),fantmp;
  vector<ColumnVector> Hess_vec_tmp; vector<Matrix> Hess;
  ftmp=m_f;thtmp=m_th;phtmp=m_ph; 
  
  if (m_return_fanning){
      fantmp=m_fanning_angles;
      Hess_vec_tmp=m_invprHes_e1;
      Hess=m_Hessian;
  }
  
  for(int i=1;i<=nfib;i++){
    pair<float,int> p(m_f(i),i);
    fvals[i-1] = p;
  }
  std::sort(fvals.begin(),fvals.end());
  for(int i=1,ii=nfib-1;ii>=0;i++,ii--){
    m_f(i)  = ftmp(fvals[ii].second);
    m_th(i) = thtmp(fvals[ii].second);
    m_ph(i) = phtmp(fvals[ii].second);
    if (m_return_fanning){
      m_fanning_angles(i)=fantmp(fvals[ii].second);
      m_invprHes_e1[i-1]=Hess_vec_tmp[fvals[ii].second-1];
      m_Hessian[i-1]=Hess[fvals[ii].second-1];
    }
  }
}


//Returns 1-Sum(f_j), 1<=j<=ii. (ii<=nfib)
//Used for transforming beta to f and vice versa
float PVM_single_c::partial_fsum(ColumnVector& fs, int ii) const{
  float fsum=1.0;
  for(int j=1;j<=ii;j++)
    fsum-=fs(j);
  return fsum;
}


//If the sum of the fractions is >1, then zero as many fractions
//as necessary, so that the sum becomes smaller than 1.
void PVM_single_c::fix_fsum(ColumnVector& fs)const{
  float sumf=0;
  for(int i=1;i<=nfib;i++){
    sumf+=fs(i);
    if(sumf>=1){
      for(int j=i;j<=nfib;j++) 
	fs(j)=FSMALL;  //make the fraction almost zero
      break;
    }
  }
}


//Find the volume fractions given all the other model 
//parameters using Linear Least Squares
void PVM_single_c::fit_pvf(ColumnVector& x)const{
  ColumnVector fs(nfib);
  ColumnVector Y_I(npts);
  Matrix       M(npts,nfib),dir(3,nfib);
  float s0=x(1),d=lambda2d(x(2)), f0;
  for(int k=1;k<=nfib;k++){
    int kk = 3+3*(k-1);
    dir(1,k) = sin(x(kk+1))*cos(x(kk+2));
    dir(2,k) = sin(x(kk+1))*sin(x(kk+2));
    dir(3,k) = cos(x(kk+1));
  }
  ////////////////////////////////////
  for(int i=1;i<=npts;i++){
    float Iso_term=isoterm(i,d);
    Y_I(i) = Y(i)-s0*Iso_term;
    for(int k=1;k<=nfib;k++){
      M(i,k)=s0*(anisoterm(i,d,dir.Column(k))-Iso_term);
    }
    //if (m_include_f0)
    //M(i,f_num)=s0*(1-Iso_term);
  }
  fs = pinv(M)*Y_I;
  if (m_include_f0){
    f0=FSMALL; fs(1)-=f0;  //Initialize f0 with a very small value
  }

  for(int k=1;k<=nfib;k++)
    fs(k)=fabs(fs(k));    //make sure that the initial values for the fractions are positive

  fix_fsum(fs);
  
  for(int k=1;k<=nfib;k++){
    float tmpr=fs(k)/partial_fsum(fs,k-1);
    if (tmpr>1.0) tmpr=1; //This can be due to numerical errors
    x(3+3*(k-1))=f2beta(tmpr);
  }

  if (m_include_f0){
    float tmpr=f0/partial_fsum(fs,nfib);
    if (tmpr>1.0) tmpr=1; //This can be due to numerical errors
    x(nparams)=f2beta(tmpr);
  }
}



//Print the final estimates (after having them transformed)
void PVM_single_c::print()const{
  cout << "PVM (Single) FIT RESULTS " << endl;
  cout << "S0   :" << m_s0 << endl;
  cout << "D    :" << m_d << endl;
  for(int i=1;i<=nfib;i++){
    cout << "F" << i << "   :" << m_f(i) << endl;
    ColumnVector x(3);
    x << sin(m_th(i))*cos(m_ph(i)) << sin(m_th(i))*sin(m_ph(i)) << cos(m_th(i));
    if(x(3)<0)x=-x;
    float _th,_ph;cart2sph(x,_th,_ph);
    cout << "TH" << i << "  :" << _th << endl; 
    cout << "PH" << i << "  :" << _ph << endl; 
    cout << "DIR" << i << "   : " << x(1) << " " << x(2) << " " << x(3) << endl;
  }
  if (m_include_f0)
    cout << "F0    :" << m_f0 << endl;
  if (m_eval_BIC)
    cout<< "BIC  :"<<m_BIC<<endl;
}



//Print the estimates using a vector with the untransformed parameter values
void PVM_single_c::print(const ColumnVector& p)const{
  ColumnVector f(nfib);

  cout << "PARAMETER VALUES " << endl;
  cout << "S0   :" << p(1) << endl;
  cout << "D    :" << lambda2d(p(2)) << endl;
  for(int i=3,ii=1;ii<=nfib;i+=3,ii++){
    f(ii) = beta2f(p(i))*partial_fsum(f,ii-1);
    cout << "F" << ii << "   :" << f(ii) << endl;
    cout << "TH" << ii << "  :" << p(i+1)*180.0/M_PI << " deg" << endl; 
    cout << "PH" << ii << "  :" << p(i+2)*180.0/M_PI << " deg" << endl; 
  }
  if (m_include_f0)
    cout << "F0    :" << beta2f(p(nparams))*partial_fsum(f,nfib);
}




ReturnMatrix PVM_single_c::get_prediction()const{
  ColumnVector pred(npts);
  ColumnVector p(nparams);
  ColumnVector fs(nfib);
  
  fs=m_f;
  p(1) = m_s0;
  p(2) = d2lambda(m_d);
  for(int i=3,ii=1;ii<=nfib;i+=3,ii++){
    float tmpr=m_f(ii)/partial_fsum(fs,ii-1);
    if (tmpr>1.0) tmpr=1; //This can be due to numerical errors
    p(i)   = f2beta(tmpr);
    p(i+1) = m_th(ii);
    p(i+2) = m_ph(ii);
  }
  if (m_include_f0){
    float tmpr=m_f0/partial_fsum(fs,nfib);
    if (tmpr>1.0) tmpr=1; //This can be due to numerical errors
    p(nparams)=f2beta(tmpr);
  }

  pred = forwardModel(p);
  pred.Release();
  return pred;
}


NEWMAT::ReturnMatrix PVM_single_c::forwardModel(const NEWMAT::ColumnVector& p)const{
  ColumnVector pred(npts);
  pred = 0;
  float val;
  float _d = lambda2d(p(2));
  ////////////////////////////////////
  ColumnVector fs(nfib);
  Matrix x(nfib,3);
  float sumf=0;
  for(int k=1;k<=nfib;k++){
    int kk = 3+3*(k-1);
    fs(k) = beta2f(p(kk))*partial_fsum(fs,k-1);
    sumf += fs(k);
    x(k,1) = sin(p(kk+1))*cos(p(kk+2));
    x(k,2) = sin(p(kk+1))*sin(p(kk+2));
    x(k,3) = cos(p(kk+1));
  }
  ////////////////////////////////////
  for(int i=1;i<=Y.Nrows();i++){
    val = 0.0;
    for(int k=1;k<=nfib;k++){
      val += fs(k)*anisoterm(i,_d,x.Row(k).t());
    }
    if (m_include_f0){
      float temp_f0=beta2f(p(nparams))*partial_fsum(fs,nfib);
      pred(i) = p(1)*(temp_f0+(1-sumf-temp_f0)*isoterm(i,_d)+val);
    } 
    else
      pred(i) = p(1)*((1-sumf)*isoterm(i,_d)+val); 
  }  
  pred.Release();
  return pred;
}




//Cost Function, sum of squared residuals
//assume that parameter values p are untransformed (e.g. need to transform them to get d, f's)
double PVM_single_c::cf(const NEWMAT::ColumnVector& p)const{
  //cout<<"CF"<<endl;
  //OUT(p.t());
  double cfv = 0.0;
  double err;
  float _d = lambda2d(p(2));
  ////////////////////////////////////
  ColumnVector fs(nfib);
  Matrix x(nfib,3);
  float sumf=0;
  for(int k=1;k<=nfib;k++){
    int kk = 3+3*(k-1);
    fs(k) = beta2f(p(kk))*partial_fsum(fs,k-1);
    sumf += fs(k);
    x(k,1) = sin(p(kk+1))*cos(p(kk+2));
    x(k,2) = sin(p(kk+1))*sin(p(kk+2));
    x(k,3) = cos(p(kk+1));
  }
  ////////////////////////////////////
  for(int i=1;i<=Y.Nrows();i++){
    err = 0.0;
    for(int k=1;k<=nfib;k++){
      err += fs(k)*anisoterm(i,_d,x.Row(k).t());
    }
    if (m_include_f0){
      float temp_f0=beta2f(p(nparams))*partial_fsum(fs,nfib);
      err = (p(1)*(temp_f0+(1-sumf-temp_f0)*isoterm(i,_d)+err) - Y(i)); 
    }
    else
      err = (p(1)*((1-sumf)*isoterm(i,_d)+err) - Y(i)); 
    cfv += err*err; 
  }  
  return(cfv);
}

NEWMAT::ReturnMatrix PVM_single_c::grad(const NEWMAT::ColumnVector& p)const{
  //cout<<"GRAD"<<endl;
  //OUT(p.t());
  NEWMAT::ColumnVector gradv(p.Nrows());
  gradv = 0.0;
  float _d = lambda2d(p(2));
  ////////////////////////////////////
  ColumnVector fs(nfib);
  ColumnVector bs(nfib);
  Matrix x(nfib,3);ColumnVector xx(3); ColumnVector yy(3);
  float sumf=0;
  
  for(int k=1;k<=nfib;k++){
    int kk = 3+3*(k-1);
    bs(k)=p(kk);
    fs(k) = beta2f(p(kk))*partial_fsum(fs,k-1);
    sumf += fs(k);
    x(k,1) = sin(p(kk+1))*cos(p(kk+2));
    x(k,2) = sin(p(kk+1))*sin(p(kk+2));
    x(k,3) = cos(p(kk+1));
  }

  ////////////////////////////////////
  Matrix f_deriv;
  //Compute the derivatives with respect to betas, i.e the transformed volume fraction variables
  f_deriv=fractions_deriv(nfib, fs, bs);  

  Matrix J(npts,nparams);     //Get the Jacobian of the model equation. The derivatives of the cost function
                              //for parameter j will then be: Grad_j=Sum(2*(F(x_i)-Y_i)J(i,j)), Sum across data points i
  ColumnVector diff(npts);
  float sig, Iso_term; 
  ColumnVector Aniso_terms(nfib);

  for(int i=1;i<=Y.Nrows();i++){
    Iso_term=isoterm(i,_d);  //Precompute some terms for this datapoint
    for(int k=1;k<=nfib;k++){
      xx = x.Row(k).t();
      Aniso_terms(k)=anisoterm(i,_d,xx);
    }
    sig = 0;
    J.Row(i)=0;
    for(int k=1;k<=nfib;k++){
      int kk = 3+3*(k-1);
      xx = x.Row(k).t();
      sig += fs(k)*Aniso_terms(k); //Total signal
      // other stuff for derivatives
      // lambda (i.e. d)
      J(i,2) += p(1)*fs(k)*anisoterm_lambda(i,p(2),xx);
      
      // beta (i.e. f)
      J(i,kk)=0;
      for (int j=1; j<=nfib; j++){
	if (f_deriv(j,k)!=0)
	  J(i,kk) += p(1)*(Aniso_terms(j)-Iso_term)*f_deriv(j,k);
      }
      // th
      J(i,kk+1) = p(1)*fs(k)*anisoterm_th(i,_d,xx,p(kk+1),p(kk+2));
      // ph
      J(i,kk+2) = p(1)*fs(k)*anisoterm_ph(i,_d,xx,p(kk+1),p(kk+2));
    }
    if (m_include_f0){
      float temp_f0=beta2f(p(nparams))*partial_fsum(fs,nfib);
      //derivative with respect to f0
      J(i,nparams)= p(1)*(1-Iso_term)*sin(2*p(nparams))*partial_fsum(fs,nfib);
      sig=p(1)*(temp_f0+(1-sumf-temp_f0)*Iso_term+sig);
      J(i,2) += p(1)*(1-sumf-temp_f0)*isoterm_lambda(i,p(2));
    }
    else{
      sig = p(1)*((1-sumf)*Iso_term+sig);
      J(i,2) += p(1)*(1-sumf)*isoterm_lambda(i,p(2)); //lambda
    }
    diff(i) = sig - Y(i);
    J(i,1) = sig/p(1);  //S0
  }
  gradv = 2*J.t()*diff;
  gradv.Release();
  return gradv;
}


//this uses Gauss-Newton approximation, i.e Hij ~ Sum(Gj*Gk), with G the derivative of the cost function with respect to parameter j
//and the Sum over all data points.
boost::shared_ptr<BFMatrix> PVM_single_c::hess(const NEWMAT::ColumnVector& p,boost::shared_ptr<BFMatrix> iptr)const{
  //cout<<"HESS"<<endl;
  //OUT(p.t());
  boost::shared_ptr<BFMatrix>   hessm;
  if (iptr && iptr->Nrows()==(unsigned int)p.Nrows() && iptr->Ncols()==(unsigned int)p.Nrows()) hessm = iptr;
  else hessm = boost::shared_ptr<BFMatrix>(new FullBFMatrix(p.Nrows(),p.Nrows()));

  float _d = lambda2d(p(2));
  ////////////////////////////////////
  ColumnVector fs(nfib);
  ColumnVector bs(nfib);
  Matrix x(nfib,3);ColumnVector xx(3); ColumnVector yy(3);
  float sumf=0;
  for(int k=1;k<=nfib;k++){
    int kk = 3+3*(k-1);
    bs(k)=p(kk);
    fs(k) = beta2f(p(kk))*partial_fsum(fs,k-1);
    sumf += fs(k);
    x(k,1) = sin(p(kk+1))*cos(p(kk+2));
    x(k,2) = sin(p(kk+1))*sin(p(kk+2));
    x(k,3) = cos(p(kk+1));
  }
  ////////////////////////////////////
  Matrix f_deriv;
  f_deriv=fractions_deriv(nfib, fs, bs);

  Matrix J(npts,nparams);
  float sig, Iso_term; 
  ColumnVector Aniso_terms(nfib);

  for(int i=1;i<=Y.Nrows();i++){
    Iso_term=isoterm(i,_d);  //Precompute some terms for this datapoint
    for(int k=1;k<=nfib;k++){
      xx = x.Row(k).t();
      Aniso_terms(k)=anisoterm(i,_d,xx);
    }
    sig = 0;
    J.Row(i)=0;
    for(int k=1;k<=nfib;k++){
      int kk = 3+3*(k-1);
      xx = x.Row(k).t();
      sig += fs(k)*Aniso_terms(k); //Total signal
      // other stuff for derivatives
      // lambda (i.e. d)
      J(i,2) += p(1)*fs(k)*anisoterm_lambda(i,p(2),xx);
      
      // beta (i.e. f)
      J(i,kk)=0;
      for (int j=1; j<=nfib; j++){
	if (f_deriv(j,k)!=0)
	  J(i,kk) += p(1)*(Aniso_terms(j)-Iso_term)*f_deriv(j,k);
      }
      // th
      J(i,kk+1) = p(1)*fs(k)*anisoterm_th(i,_d,xx,p(kk+1),p(kk+2));
      // ph
      J(i,kk+2) = p(1)*fs(k)*anisoterm_ph(i,_d,xx,p(kk+1),p(kk+2));
    }
    if (m_include_f0){
      float temp_f0=beta2f(p(nparams))*partial_fsum(fs,nfib);
      //derivative with respect to f0
      J(i,nparams)= p(1)*(1-Iso_term)*sin(2*p(nparams))*partial_fsum(fs,nfib);
      sig=p(1)*(temp_f0+(1-sumf-temp_f0)*Iso_term+sig);
      J(i,2) += p(1)*(1-sumf-temp_f0)*isoterm_lambda(i,p(2));
    }
    else{
      sig = p(1)*((1-sumf)*Iso_term+sig);
      J(i,2) += p(1)*(1-sumf)*isoterm_lambda(i,p(2)); //lambda
    }
    J(i,1) = sig/p(1);  //S0
  }

  for (int i=1; i<=p.Nrows(); i++){
    for (int j=i; j<=p.Nrows(); j++){
      sig = 0.0;
      for(int k=1;k<=J.Nrows();k++)
	sig += 2*(J(k,i)*J(k,j));
      hessm->Set(i,j,sig);
    }
  }
  for (int j=1; j<=p.Nrows(); j++) {
    for (int i=j+1; i<=p.Nrows(); i++) {
      hessm->Set(i,j,hessm->Peek(j,i));
    }
  }
  return(hessm);
}


//Once the model is fitted: For each fibre, compute a 3x3 Hessian of the Cost function at the cartesian (x,y,z) coordinates of the orientation,
//evaluated at the estimated parameters. Return the second eigenvector of the (inverse) Hessian - No need to invert though! 
void PVM_single_c::eval_Hessian_at_peaks(){
  vector<ColumnVector> V;
  ColumnVector temp_vec(3);
  float fiso=1,dot, Saniso,dSaniso_dx,dSaniso_dy,dSaniso_dz;
  ColumnVector S(npts), F(npts);
  Matrix dS_dx(npts,nfib), dS_dy(npts,nfib), dS_dz(npts,nfib);
  Matrix dS_dxdx(npts,nfib), dS_dydy(npts,nfib), dS_dzdz(npts,nfib), dS_dxdy(npts,nfib),dS_dxdz(npts,nfib),dS_dydz(npts,nfib); 
  
  for (int n=1; n<=nfib; n++){
    temp_vec<< cos(m_ph(n))*sin(m_th(n)) << sin(m_ph(n))*sin(m_th(n)) << cos(m_th(n));
    V.push_back(temp_vec); 
    fiso-=m_f(n);
  }
  
  for (int i=1; i<=npts; i++){      
    float bd=bvals(1,i)*m_d;
    S(i)=fiso*exp(-bd);
    for (int n=1; n<=nfib; n++){
      dot=V[n-1](1)*bvecs(1,i)+V[n-1](2)*bvecs(2,i)+V[n-1](3)*bvecs(3,i);
      
      Saniso=exp(-bd*dot*dot);
      S(i)=S(i)+m_f(n)*Saniso;
  
      float dot1=-2*bd*dot*Saniso; float dot2=-2*bd*dot; float dot3=-2*bd*Saniso; float fS0=m_s0*m_f(n);
      //First Derivatives of the signal wrt x,y,z
      dSaniso_dx=dot1*bvecs(1,i);  dSaniso_dy=dot1*bvecs(2,i); dSaniso_dz=dot1*bvecs(3,i);
      dS_dx(i,n)=fS0*dSaniso_dx;   dS_dy(i,n)=fS0*dSaniso_dy;  dS_dz(i,n)=fS0*dSaniso_dz;
      //Second Derivatives of the signal wrt x,y,z             
      dS_dxdx(i,n)=fS0*bvecs(1,i)*(dSaniso_dx*dot2+dot3*bvecs(1,i));  dS_dydy(i,n)=fS0*bvecs(2,i)*(dSaniso_dy*dot2+dot3*bvecs(2,i));
      dS_dzdz(i,n)=fS0*bvecs(3,i)*(dSaniso_dz*dot2+dot3*bvecs(3,i));  dS_dxdy(i,n)=fS0*bvecs(1,i)*(dSaniso_dy*dot2+dot3*bvecs(2,i));
      dS_dxdz(i,n)=fS0*bvecs(1,i)*(dSaniso_dz*dot2+dot3*bvecs(3,i));  dS_dydz(i,n)=fS0*bvecs(2,i)*(dSaniso_dz*dot2+dot3*bvecs(3,i));
    }
    F(i)=Y(i)-m_s0*S(i);
  }

  for (int n=1; n<=nfib; n++){   //For each fibre, compute the Hessian matrix of the cost function at (x,y,z)
    SymmetricMatrix H(3); H=0;
    for (int i=1; i<=npts; i++){       
      H(1,1)-=dS_dx(i,n)*dS_dx(i,n)-F(i)*dS_dxdx(i,n);
      H(2,2)-=dS_dy(i,n)*dS_dy(i,n)-F(i)*dS_dydy(i,n);
      H(3,3)-=dS_dz(i,n)*dS_dz(i,n)-F(i)*dS_dzdz(i,n);
      H(1,2)-=dS_dx(i,n)*dS_dy(i,n)-F(i)*dS_dxdy(i,n);
      H(1,3)-=dS_dx(i,n)*dS_dz(i,n)-F(i)*dS_dxdz(i,n);
      H(2,3)-=dS_dy(i,n)*dS_dz(i,n)-F(i)*dS_dydz(i,n);
    }
    H=-2*H;
    m_Hessian.push_back(H); //store the Hessian
  } 
}      



//For each fibre, get the projection of the inverse Hessian to the fanning plane. 
//Its first eigenvector will be utilized to get a fanning angle in [0,pi).
void PVM_single_c::Fanning_angles_from_Hessian(){
  Matrix Rth(3,3), Rph(3,3), R(3,3), H, A(3,3), P(3,3), E;
  DiagonalMatrix L; SymmetricMatrix Q(3);
  ColumnVector e1(3),vfib(3),v2(3),v3(3);
 
  eval_Hessian_at_peaks();  //Compute the Hessian for each fibre orientation
  
  //Then project the inverse Hessian to the fanning plane (perpendicular to the orientation) and obtain its first eigenvector
  for (int n=1; n<=nfib; n++){      //For each fitted fibre
    P << 0 << 0 << 0
      << 0 << 1 << 0
      << 0 << 0 << 1;
    H=m_Hessian[n-1];
    float sinth=sin(m_th(n)), costh=cos(m_th(n));
    float sinph=sin(m_ph(n)), cosph=cos(m_ph(n));
    vfib<< sinth*cosph << sinth*sinph << costh; //Corresponding fibre orientation
      //Define two vectors that are orthogonal to vfib
    if (vfib(1)==0 && vfib(2)==0) //then, we have a [0 0 1] orientation
      v2<< 1 << 0 << 0;
    else
      v2 << vfib(2) << -vfib(1) << 0; //define v2, so that vfib*v2=0;
    float mag=sqrt(v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3));
    v2=v2/mag;
    v3 << vfib(2)*v2(3)-vfib(3)*v2(2) << vfib(3)*v2(1)-vfib(1)*v2(3) << vfib(1)*v2(2)-vfib(2)*v2(1); //Now get the cross product
 
    //Define a Projection Matrix to the plane perpendicular to the fibre orientation vfib   
    A.Row(1)<<vfib.t(); A.Row(2)<<v2.t(); A.Row(3)<<v3.t();
    P= P*A;
    
    try{                   //If the Hessian is invertible (might not be in some background voxels)
      Q << P*H.i()*P.t();  //Project the inverse Hessian to this plane
      EigenValues(Q,L,E);  //Eigendecompose the projected inverse Hessian
      e1 = E.Column(3);    //Projected Hessian 1st eigenvector
    }
    catch(...){            //Otherwise give a random value and proceed        
      e1<<1<<0<<0;
    }
    e1 << A.t()*e1;
    m_invprHes_e1.push_back(e1);
    //float dot=DotProduct(e1,vfib); ColumnVector p;
    //if (dot<0) {e1=-e1; dot=-dot; }
    //p=e1-dot*vfib; //Projection of e2 on a plane perpendicular to vfib
    //p=p/sqrt(p(1)*p(1)+p(2)*p(2)+p(3)*p(3));
    //e1=p; 
    
    Rth<<costh << 0 << -sinth
       <<0     << 1 << 0
       <<sinth << 0 << costh;
    Rph<<cosph << sinph <<0
       <<-sinph<< cosph <<0
       <<0     << 0     <<1;
    R=Rth*Rph;      //Rotation Matrix for vfib to become parallel to z
    
    e1<< R*e1;      //Rotate the fanning plane effectively to xy plane
    m_fanning_angles(n)=atan2(-e1(1),e1(2));  //this gives [-pi,pi] range
    if (m_fanning_angles(n)<0)  m_fanning_angles(n)+=M_PI;  //transform it to [0,pi]
    //cout<<m_fanning_angles(n)<<endl<<endl;  
  }
}



////////////////////////////////////////////////
//       PARTIAL VOLUME MODEL - SINGLE SHELL (OLD)
////////////////////////////////////////////////

void PVM_single::fit(){

  // initialise with a tensor
  DTI dti(Y,bvecs,bvals);
  dti.linfit();
  //dti.print();

  // set starting parameters for nonlinear fitting
  float _th,_ph;
  cart2sph(dti.get_v1(),_th,_ph);

  ColumnVector start(nparams);
  start(1) = dti.get_s0();
  start(2) = dti.get_md()>0?dti.get_md()*2:0.001; // empirically found that d~2*MD
  //start(2) = dti.get_l1()>0?dti.get_l1():0.002; // empirically found that d~L1
  start(3) = dti.get_fa()<1?f2x(dti.get_fa()):f2x(0.95); // first pvf = FA 
  start(4) = _th;
  start(5) = _ph;
  float sumf=x2f(start(3));
  float tmpsumf=sumf;
  for(int ii=2,i=6;ii<=nfib;ii++,i+=3){
    float denom=2;
    do{
      start(i) = f2x(x2f(start(i-3))/denom);
      denom *= 2;
      tmpsumf = sumf + x2f(start(i));
    }while(tmpsumf>=1);
    sumf += x2f(start(i));
    cart2sph(dti.get_v(ii),_th,_ph);
    start(i+1) = _th;
    start(i+2) = _ph;
  }
  if (m_include_f0)
    start(nparams)=f2x(FSMALL);
 
  // do the fit
  NonlinParam  lmpar(start.Nrows(),NL_LM); 
  lmpar.SetGaussNewtonType(LM_L);
  lmpar.SetStartingEstimate(start);

  NonlinOut status;
  status = nonlin(lmpar,(*this));
  ColumnVector final_par(nparams);
  final_par = lmpar.Par();


  // finalise parameters
  m_s0 = final_par(1);
  m_d  = std::abs(final_par(2));
  for(int k=1;k<=nfib;k++){
    int kk = 3 + 3*(k-1);
    m_f(k)  = x2f(final_par(kk));
    m_th(k) = final_par(kk+1);
    m_ph(k) = final_par(kk+2);
  }
  if (m_include_f0)
    m_f0=x2f(final_par(nparams));
  sort();
  fix_fsum();
}

void PVM_single::sort(){
  vector< pair<float,int> > fvals(nfib);
  ColumnVector ftmp(nfib),thtmp(nfib),phtmp(nfib);
  ftmp=m_f;thtmp=m_th;phtmp=m_ph;
  for(int i=1;i<=nfib;i++){
    pair<float,int> p(m_f(i),i);
    fvals[i-1] = p;
  }
  std::sort(fvals.begin(),fvals.end());
  for(int i=1,ii=nfib-1;ii>=0;i++,ii--){
    m_f(i)  = ftmp(fvals[ii].second);
    m_th(i) = thtmp(fvals[ii].second);
    m_ph(i) = phtmp(fvals[ii].second);
  }
}

void PVM_single::fix_fsum(){
  float sumf=0;
  if (m_include_f0) 
    sumf=m_f0;
  for(int i=1;i<=nfib;i++){
    sumf+=m_f(i);
    if(sumf>=1){for(int j=i;j<=nfib;j++)m_f(j)=FSMALL; break;}
  }
}

ReturnMatrix PVM_single::get_prediction()const{
  ColumnVector pred(npts);
  ColumnVector p(nparams);
  p(1) = m_s0;
  p(2) = m_d;
  for(int i=3,ii=1;ii<=nfib;i+=3,ii++){
    p(i)   = f2x(m_f(ii));
    p(i+1) = m_th(ii);
    p(i+2) = m_ph(ii);
  }
  if (m_include_f0)
    p(nparams)=f2x(m_f0);
  pred = forwardModel(p);

  pred.Release();
  return pred;
}

NEWMAT::ReturnMatrix PVM_single::forwardModel(const NEWMAT::ColumnVector& p)const{
  //cout<<"FORWARD"<<endl;
  //OUT(p.t());
  ColumnVector pred(npts);
  pred = 0;
  float val;
  float _d = std::abs(p(2));
  
  ////////////////////////////////////
  ColumnVector fs(nfib);
  Matrix x(nfib,3);
  float sumf=0;
  for(int k=1;k<=nfib;k++){
    int kk = 3+3*(k-1);
    fs(k) = x2f(p(kk));
    sumf += fs(k);
    x(k,1) = sin(p(kk+1))*cos(p(kk+2));
    x(k,2) = sin(p(kk+1))*sin(p(kk+2));
    x(k,3) = cos(p(kk+1));
  }
  ////////////////////////////////////
  for(int i=1;i<=Y.Nrows();i++){
    val = 0.0;
    for(int k=1;k<=nfib;k++){
      val += fs(k)*anisoterm(i,_d,x.Row(k).t());
    }
    if (m_include_f0){
      float temp_f0=x2f(p(nparams));
      pred(i) = p(1)*(temp_f0+(1-sumf-temp_f0)*isoterm(i,_d)+val);
    } 
    else
      pred(i) = p(1)*((1-sumf)*isoterm(i,_d)+val); 
  }  
  pred.Release();
  //cout<<"----"<<endl;
  return pred;
}


double PVM_single::cf(const NEWMAT::ColumnVector& p)const{
  //cout<<"CF"<<endl;
  //OUT(p.t());
  double cfv = 0.0;
  double err;
  float _d = std::abs(p(2));
  ////////////////////////////////////
  ColumnVector fs(nfib);
  Matrix x(nfib,3);
  float sumf=0;
  for(int k=1;k<=nfib;k++){
    int kk = 3+3*(k-1);
    fs(k) = x2f(p(kk));
    sumf += fs(k);
    x(k,1) = sin(p(kk+1))*cos(p(kk+2));
    x(k,2) = sin(p(kk+1))*sin(p(kk+2));
    x(k,3) = cos(p(kk+1));
  }
  ////////////////////////////////////
  for(int i=1;i<=Y.Nrows();i++){
    err = 0.0;
    for(int k=1;k<=nfib;k++){
      err += fs(k)*anisoterm(i,_d,x.Row(k).t());
    }
    if (m_include_f0){
      float temp_f0=x2f(p(nparams));
      err = (p(1)*(temp_f0+(1-sumf-temp_f0)*isoterm(i,_d)+err) - Y(i)); 
    }
    else
      err = (p(1)*((1-sumf)*isoterm(i,_d)+err) - Y(i)); 
    cfv += err*err; 
  }  
  //OUT(cfv);
  //cout<<"----"<<endl;
  return(cfv);
}


NEWMAT::ReturnMatrix PVM_single::grad(const NEWMAT::ColumnVector& p)const{
  //cout<<"GRAD"<<endl;
  //OUT(p.t());
  NEWMAT::ColumnVector gradv(p.Nrows());
  gradv = 0.0;
  float _d = std::abs(p(2));
  ////////////////////////////////////
  ColumnVector fs(nfib);
  Matrix x(nfib,3);ColumnVector xx(3);
  float sumf=0;
  for(int k=1;k<=nfib;k++){
    int kk = 3+3*(k-1);
    fs(k) = x2f(p(kk));
    sumf += fs(k);
    x(k,1) = sin(p(kk+1))*cos(p(kk+2));
    x(k,2) = sin(p(kk+1))*sin(p(kk+2));
    x(k,3) = cos(p(kk+1));
  }
  ////////////////////////////////////
  Matrix J(npts,nparams);
  ColumnVector diff(npts);
  float sig;
  for(int i=1;i<=Y.Nrows();i++){
    sig = 0;
    J.Row(i)=0;
    for(int k=1;k<=nfib;k++){
      int kk = 3+3*(k-1);
      xx = x.Row(k).t();
      sig += fs(k)*anisoterm(i,_d,xx);
      // other stuff for derivatives
      // d
      J(i,2) += (p(2)>0?1.0:-1.0)*p(1)*fs(k)*anisoterm_d(i,_d,xx);
      // f
      J(i,kk) = p(1)*(anisoterm(i,_d,xx)-isoterm(i,_d)) * two_pi*sign(p(kk))*1/(1+p(kk)*p(kk));
      // th
      J(i,kk+1) = p(1)*fs(k)*anisoterm_th(i,_d,xx,p(kk+1),p(kk+2));
      // ph
      J(i,kk+2) = p(1)*fs(k)*anisoterm_ph(i,_d,xx,p(kk+1),p(kk+2));
    }
    if (m_include_f0){
      float temp_f0=x2f(p(nparams));
      //derivative with respect to f0
      J(i,nparams)= p(1)*(1-isoterm(i,_d)) * two_pi*sign(p(nparams))*1/(1+p(nparams)*p(nparams));
      sig=p(1)*(temp_f0+(1-sumf-temp_f0)*isoterm(i,_d)+sig);
      J(i,2) += (p(2)>0?1.0:-1.0)*p(1)*(1-sumf-temp_f0)*isoterm_d(i,_d);
    }
    else{
      sig = p(1)*((1-sumf)*isoterm(i,_d)+sig);
      J(i,2) += (p(2)>0?1.0:-1.0)*p(1)*(1-sumf)*isoterm_d(i,_d);
    }
    diff(i) = sig - Y(i);
    J(i,1) = sig/p(1);
  }
  
  gradv = 2*J.t()*diff;
  gradv.Release();
  return gradv;


}

//this uses Gauss-Newton approximation
boost::shared_ptr<BFMatrix> PVM_single::hess(const NEWMAT::ColumnVector& p,boost::shared_ptr<BFMatrix> iptr)const{
  //cout<<"HESS"<<endl;
  //OUT(p.t());
  boost::shared_ptr<BFMatrix>   hessm;
  if (iptr && iptr->Nrows()==(unsigned int)p.Nrows() && iptr->Ncols()==(unsigned int)p.Nrows()) hessm = iptr;
  else hessm = boost::shared_ptr<BFMatrix>(new FullBFMatrix(p.Nrows(),p.Nrows()));

  float _d = std::abs(p(2));
  ////////////////////////////////////
  ColumnVector fs(nfib);
  Matrix x(nfib,3);ColumnVector xx(3);
  float sumf=0;
  for(int k=1;k<=nfib;k++){
    int kk = 3+3*(k-1);
    fs(k) = x2f(p(kk));
    sumf += fs(k);
    x(k,1) = sin(p(kk+1))*cos(p(kk+2));
    x(k,2) = sin(p(kk+1))*sin(p(kk+2));
    x(k,3) = cos(p(kk+1));
  }
   ////////////////////////////////////
  Matrix J(npts,nparams);
  float sig;
  for(int i=1;i<=Y.Nrows();i++){
    sig = 0;
    J.Row(i)=0;
    for(int k=1;k<=nfib;k++){
      int kk = 3+3*(k-1);
      xx = x.Row(k).t();
      sig += fs(k)*anisoterm(i,_d,xx);
      // other stuff for derivatives
      // d
      J(i,2) += (p(2)>0?1.0:-1.0)*p(1)*fs(k)*anisoterm_d(i,_d,xx);
      // f
      J(i,kk) = p(1)*(anisoterm(i,_d,xx)-isoterm(i,_d)) * two_pi*sign(p(kk))*1/(1+p(kk)*p(kk));
      // th
      J(i,kk+1) = p(1)*fs(k)*anisoterm_th(i,_d,xx,p(kk+1),p(kk+2));
      // ph
      J(i,kk+2) = p(1)*fs(k)*anisoterm_ph(i,_d,xx,p(kk+1),p(kk+2));
    }
    if (m_include_f0){
      float temp_f0=x2f(p(nparams));
      //derivative with respect to f0
      J(i,nparams)= p(1)*(1-isoterm(i,_d)) * two_pi*sign(p(nparams))*1/(1+p(nparams)*p(nparams));
      sig=p(1)*(temp_f0+(1-sumf-temp_f0)*isoterm(i,_d)+sig);
      J(i,2) += (p(2)>0?1.0:-1.0)*p(1)*(1-sumf-temp_f0)*isoterm_d(i,_d);
    }
    else{
      sig = p(1)*((1-sumf)*isoterm(i,_d)+sig);
      J(i,2) += (p(2)>0?1.0:-1.0)*p(1)*(1-sumf)*isoterm_d(i,_d);
    }
    J(i,1) = sig/p(1);
  }
  
  for (int i=1; i<=p.Nrows(); i++){
    for (int j=i; j<=p.Nrows(); j++){
      sig = 0.0;
      for(int k=1;k<=J.Nrows();k++)
	sig += 2*(J(k,i)*J(k,j));
      hessm->Set(i,j,sig);
    }
  }
  for (int j=1; j<=p.Nrows(); j++) {
    for (int i=j+1; i<=p.Nrows(); i++) {
      hessm->Set(i,j,hessm->Peek(j,i));
    }
  }
  //hessm->Print();
  //cout<<"----"<<endl;
  return(hessm);
}





////////////////////////////////////////////////
//       PARTIAL VOLUME MODEL - MULTIPLE SHELLS
////////////////////////////////////////////////

void PVM_multi::fit(){

  // initialise with simple pvm
  PVM_single pvm1(Y,bvecs,bvals,nfib,m_include_f0);
  pvm1.fit();
  //cout<<"Init with single"<<endl;
  //pvm1.print();

  float _a,_b;
  _a = 1.0; // start with d=d_std
  _b = pvm1.get_d();

  ColumnVector start(nparams);
  start(1) = pvm1.get_s0();
  start(2) = _a;
  start(3) = _b;
  for(int i=1,ii=4;i<=nfib;i++,ii+=3){
    start(ii) = f2x(pvm1.get_f(i));
    start(ii+1) = pvm1.get_th(i);
    start(ii+2) = pvm1.get_ph(i);
  }
  if (m_include_f0)
    start(nparams)=f2x(pvm1.get_f0());
  
  // do the fit
  NonlinParam  lmpar(start.Nrows(),NL_LM); 
  lmpar.SetGaussNewtonType(LM_L);
  lmpar.SetStartingEstimate(start);

  NonlinOut status;
  status = nonlin(lmpar,(*this));
  ColumnVector final_par(nparams);
  final_par = lmpar.Par();

  // finalise parameters
  m_s0     = final_par(1);
  m_d      = std::abs(final_par(2)*final_par(3));
  m_d_std  = std::sqrt(std::abs(final_par(2)*final_par(3)*final_par(3)));
  for(int i=4,k=1;k<=nfib;i+=3,k++){
    m_f(k)  = x2f(final_par(i));
    m_th(k) = final_par(i+1);
    m_ph(k) = final_par(i+2);
  }
  if (m_include_f0)
    m_f0=x2f(final_par(nparams));
  sort();
  fix_fsum();
}


void PVM_multi::sort(){
  vector< pair<float,int> > fvals(nfib);
  ColumnVector ftmp(nfib),thtmp(nfib),phtmp(nfib);
  ftmp=m_f;thtmp=m_th;phtmp=m_ph;
  for(int i=1;i<=nfib;i++){
    pair<float,int> p(m_f(i),i);
    fvals[i-1] = p;
  }
  std::sort(fvals.begin(),fvals.end());
  for(int i=1,ii=nfib-1;ii>=0;i++,ii--){
    m_f(i)  = ftmp(fvals[ii].second);
    m_th(i) = thtmp(fvals[ii].second);
    m_ph(i) = phtmp(fvals[ii].second);
  }
}

void PVM_multi::fix_fsum(){
  float sumf=0;
  if (m_include_f0) 
    sumf=m_f0;
  for(int i=1;i<=nfib;i++){
    if (m_f(i)==0) m_f(i)=FSMALL;
    sumf+=m_f(i);
    if(sumf>=1){for(int j=i;j<=nfib;j++)m_f(j)=FSMALL;break;}
  }
}


ReturnMatrix PVM_multi::get_prediction()const{
  ColumnVector pred(npts);
  ColumnVector p(nparams);
  p(1) = m_s0;
  p(2) = m_d*m_d/m_d_std/m_d_std;
  p(3) = m_d_std*m_d_std/m_d; // =1/beta
  for(int k=1;k<=nfib;k++){
    int kk = 4+3*(k-1);
    p(kk)   = f2x(m_f(k));
    p(kk+1) = m_th(k);
    p(kk+2) = m_ph(k);
  }
  if (m_include_f0)
    p(nparams)=f2x(m_f0);
  pred = forwardModel(p);
  pred.Release();
  return pred;
}


NEWMAT::ReturnMatrix PVM_multi::forwardModel(const NEWMAT::ColumnVector& p)const{
  ColumnVector pred(npts);
  pred = 0;
  float val;
  float _a = std::abs(p(2));
  float _b = std::abs(p(3));
  ////////////////////////////////////
  ColumnVector fs(nfib);
  Matrix x(nfib,3);
  float sumf=0;
  for(int k=1;k<=nfib;k++){
    int kk = 4+3*(k-1);
    fs(k) = x2f(p(kk));
    sumf += fs(k);
    x(k,1) = sin(p(kk+1))*cos(p(kk+2));
    x(k,2) = sin(p(kk+1))*sin(p(kk+2));
    x(k,3) = cos(p(kk+1));
  }
  ////////////////////////////////////
  for(int i=1;i<=Y.Nrows();i++){
    val = 0.0;
    for(int k=1;k<=nfib;k++){
      val += fs(k)*anisoterm(i,_a,_b,x.Row(k).t(),m_Gamma_for_ball_only);
    }
    if (m_include_f0){
      float temp_f0=x2f(p(nparams));
      pred(i) = std::abs(p(1))*(temp_f0+(1-sumf-temp_f0)*isoterm(i,_a,_b)+val);
    } 
    else
      pred(i) = std::abs(p(1))*((1-sumf)*isoterm(i,_a,_b)+val); 
  }   
  pred.Release();
  return pred;
}


double PVM_multi::cf(const NEWMAT::ColumnVector& p)const{
  //cout<<"CF"<<endl;
  //OUT(p.t());
  double cfv = 0.0;
  double err;
  float _a = std::abs(p(2));
  float _b = std::abs(p(3));
  ////////////////////////////////////
  ColumnVector fs(nfib);
  Matrix x(nfib,3);
  float sumf=0;
  for(int k=1;k<=nfib;k++){
    int kk = 4+3*(k-1);
    fs(k) = x2f(p(kk));
    sumf += fs(k);
    x(k,1) = sin(p(kk+1))*cos(p(kk+2));
    x(k,2) = sin(p(kk+1))*sin(p(kk+2));
    x(k,3) = cos(p(kk+1));
  }
  ////////////////////////////////////
  for(int i=1;i<=Y.Nrows();i++){
    err = 0.0;
    for(int k=1;k<=nfib;k++){
      err += fs(k)*anisoterm(i,_a,_b,x.Row(k).t(),m_Gamma_for_ball_only);
    } 
    if (m_include_f0){
      float temp_f0=x2f(p(nparams));
      err = (std::abs(p(1))*(temp_f0+(1-sumf-temp_f0)*isoterm(i,_a,_b)+err) - Y(i)); 
    }
    else
      err = (std::abs(p(1))*((1-sumf)*isoterm(i,_a,_b)+err) - Y(i)); 
    cfv += err*err; 
  }  
  //OUT(cfv);
  //cout<<"----"<<endl;
  return(cfv);
}


NEWMAT::ReturnMatrix PVM_multi::grad(const NEWMAT::ColumnVector& p)const{
  //cout<<"GRAD"<<endl;
  //OUT(p.t());
  NEWMAT::ColumnVector gradv(p.Nrows());
  gradv = 0.0;
  float _a = std::abs(p(2));
  float _b = std::abs(p(3));
  ////////////////////////////////////
  ColumnVector fs(nfib);
  Matrix x(nfib,3);ColumnVector xx(3);
  float sumf=0;
  for(int k=1;k<=nfib;k++){
    int kk = 4+3*(k-1);
    fs(k) = x2f(p(kk));
    sumf += fs(k);
    x(k,1) = sin(p(kk+1))*cos(p(kk+2));
    x(k,2) = sin(p(kk+1))*sin(p(kk+2));
    x(k,3) = cos(p(kk+1));
  }
  ////////////////////////////////////
  Matrix J(npts,nparams);   
  ColumnVector diff(npts);
  float sig;
  for(int i=1;i<=Y.Nrows();i++){
    sig = 0;
    J.Row(i)=0;
    for(int k=1;k<=nfib;k++){
      int kk = 4+3*(k-1);
      xx = x.Row(k).t();
      sig += fs(k)*anisoterm(i,_a,_b,xx,m_Gamma_for_ball_only);
      // other stuff for derivatives
      // alpha
      J(i,2) += (p(2)>0?1.0:-1.0)*std::abs(p(1))*fs(k)*anisoterm_a(i,_a,_b,xx,m_Gamma_for_ball_only);
      // beta
      J(i,3) += (p(3)>0?1.0:-1.0)*std::abs(p(1))*fs(k)*anisoterm_b(i,_a,_b,xx,m_Gamma_for_ball_only);
      // f
      J(i,kk) = std::abs(p(1))*(anisoterm(i,_a,_b,xx,m_Gamma_for_ball_only)-isoterm(i,_a,_b)) * two_pi*sign(p(kk))*1/(1+p(kk)*p(kk));
      // th
      J(i,kk+1) = std::abs(p(1))*fs(k)*anisoterm_th(i,_a,_b,xx,p(kk+1),p(kk+2),m_Gamma_for_ball_only);
      // ph
      J(i,kk+2) = std::abs(p(1))*fs(k)*anisoterm_ph(i,_a,_b,xx,p(kk+1),p(kk+2),m_Gamma_for_ball_only);
    }
    if (m_include_f0){
      float temp_f0=x2f(p(nparams));
      //derivative with respect to f0
      J(i,nparams)= std::abs(p(1))*(1-isoterm(i,_a,_b))*two_pi*sign(p(nparams))*1/(1+p(nparams)*p(nparams));
      sig=std::abs(p(1))*(temp_f0+(1-sumf-temp_f0)*isoterm(i,_a,_b)+sig);
      J(i,2) += (p(2)>0?1.0:-1.0)*std::abs(p(1))*(1-sumf-temp_f0)*isoterm_a(i,_a,_b);
      J(i,3) += (p(3)>0?1.0:-1.0)*std::abs(p(1))*(1-sumf-temp_f0)*isoterm_b(i,_a,_b);
    }
    else{
      sig = std::abs(p(1))*((1-sumf)*isoterm(i,_a,_b)+sig);
      J(i,2) += (p(2)>0?1.0:-1.0)*std::abs(p(1))*(1-sumf)*isoterm_a(i,_a,_b);
      J(i,3) += (p(3)>0?1.0:-1.0)*std::abs(p(1))*(1-sumf)*isoterm_b(i,_a,_b);
    }
    diff(i) = sig - Y(i);

    J(i,1) = (p(1)>0?1.0:-1.0)*sig/p(1);
  }
  
  gradv = 2*J.t()*diff;
  //OUT(gradv.t());
  //cout<<"----"<<endl;
  gradv.Release();
  return gradv;
}


//this uses Gauss-Newton approximation
boost::shared_ptr<BFMatrix> PVM_multi::hess(const NEWMAT::ColumnVector& p,boost::shared_ptr<BFMatrix> iptr)const{
  //cout<<"HESS"<<endl;
  //OUT(p.t());
  boost::shared_ptr<BFMatrix>   hessm;
  if (iptr && iptr->Nrows()==(unsigned int)p.Nrows() && iptr->Ncols()==(unsigned int)p.Nrows()) hessm = iptr;
  else hessm = boost::shared_ptr<BFMatrix>(new FullBFMatrix(p.Nrows(),p.Nrows()));

  float _a = std::abs(p(2));
  float _b = std::abs(p(3));
  ////////////////////////////////////
  ColumnVector fs(nfib);
  Matrix x(nfib,3);ColumnVector xx(3);
  float sumf=0;
  for(int k=1;k<=nfib;k++){
    int kk = 4+3*(k-1);
    fs(k) = x2f(p(kk));
    sumf += fs(k);
    x(k,1) = sin(p(kk+1))*cos(p(kk+2));
    x(k,2) = sin(p(kk+1))*sin(p(kk+2));
    x(k,3) = cos(p(kk+1));
  }
  ////////////////////////////////////
  Matrix J(npts,nparams);
  float sig;
  for(int i=1;i<=Y.Nrows();i++){
    sig = 0;
    J.Row(i)=0;
    for(int k=1;k<=nfib;k++){
      int kk = 4+3*(k-1);
      xx = x.Row(k).t();
      sig += fs(k)*anisoterm(i,_a,_b,xx,m_Gamma_for_ball_only);
      // other stuff for derivatives
      // change of variable
      float cov = two_pi*sign(p(kk))*1/(1+p(kk)*p(kk));
      // alpha
      J(i,2) += (p(2)>0?1.0:-1.0)*std::abs(p(1))*fs(k)*anisoterm_a(i,_a,_b,xx,m_Gamma_for_ball_only);
      // beta
      J(i,3) += (p(3)>0?1.0:-1.0)*std::abs(p(1))*fs(k)*anisoterm_b(i,_a,_b,xx,m_Gamma_for_ball_only);
      // f
      J(i,kk) = std::abs(p(1))*(anisoterm(i,_a,_b,xx,m_Gamma_for_ball_only)-isoterm(i,_a,_b)) * cov;
      // th
      J(i,kk+1) = std::abs(p(1))*fs(k)*anisoterm_th(i,_a,_b,xx,p(kk+1),p(kk+2),m_Gamma_for_ball_only);
      // ph
      J(i,kk+2) = std::abs(p(1))*fs(k)*anisoterm_ph(i,_a,_b,xx,p(kk+1),p(kk+2),m_Gamma_for_ball_only);
    }
    if (m_include_f0){
      float temp_f0=x2f(p(nparams));
      //derivative with respect to f0
      J(i,nparams)= std::abs(p(1))*(1-isoterm(i,_a,_b))*two_pi*sign(p(nparams))*1/(1+p(nparams)*p(nparams));
      sig=std::abs(p(1))*(temp_f0+(1-sumf-temp_f0)*isoterm(i,_a,_b)+sig);
      J(i,2) += (p(2)>0?1.0:-1.0)*std::abs(p(1))*(1-sumf-temp_f0)*isoterm_a(i,_a,_b);
      J(i,3) += (p(3)>0?1.0:-1.0)*std::abs(p(1))*(1-sumf-temp_f0)*isoterm_b(i,_a,_b);
    }
    else{
      sig = std::abs(p(1))*((1-sumf)*isoterm(i,_a,_b)+sig);
      J(i,2) += (p(2)>0?1.0:-1.0)*std::abs(p(1))*(1-sumf)*isoterm_a(i,_a,_b);
      J(i,3) += (p(3)>0?1.0:-1.0)*std::abs(p(1))*(1-sumf)*isoterm_b(i,_a,_b);
    }  
    J(i,1) = (p(1)>0?1.0:-1.0)*sig/p(1);
  }
  
  for (int i=1; i<=p.Nrows(); i++){
    for (int j=i; j<=p.Nrows(); j++){
      sig = 0.0;
      for(int k=1;k<=J.Nrows();k++)
	sig += 2*(J(k,i)*J(k,j));
      hessm->Set(i,j,sig);
    }
  }
  for (int j=1; j<=p.Nrows(); j++) {
    for (int i=j+1; i<=p.Nrows(); i++) {
      hessm->Set(i,j,hessm->Peek(j,i));
    }
  }
  //hessm->Print();
  //cout<<"----"<<endl;
  return(hessm);
}




///////////////////////////////////////////////////////////////////////////
// FANNING MODEL - BALL & BINGHAMS 
// Constrained Optimization for the diffusivity, fractions and their sum<1,
// and the Bingham eigenvalues
//////////////////////////////////////////////////////////////////////////

void PVM_Ball_Binghams::fit(){
  // Fit the ball & stick first to initialize some of the parameters
  PVM_single_c pvmbs(Y,bvecs,bvals,nfib,false,m_include_f0,true);  //Return a fanning angle estimate as well
  pvmbs.fit();
  // pvmbs.print();

  ColumnVector k1_init, w_init;
  ColumnVector final_par(nparams);
  double minRSS=1e20;
  if (!m_gridsearch){    //Initialize the fanning eigenvalues using a grid or a set of intermediate values  
    k1_init.ReSize(1); k1_init<<20;
    w_init.ReSize(1); w_init<<5.0;
  }
  else{
    k1_init.ReSize(6); k1_init<< 10 << 20 << 50 << 100 << 500 << 1000;
    w_init.ReSize(6); w_init<< 1 << 5 << 10 << 20 << 50 << 90;
  }

  for (int n1=1; n1<=k1_init.Nrows(); n1++)
    for (int n2=1; n2<=w_init.Nrows(); n2++){

      ColumnVector start(nparams);
      ColumnVector fs(nfib); fs=0;
      //Initialize the non-linear fitter. Transform all initial values to the unconstrained parameter space 
      start(1) = pvmbs.get_s0();
      start(2) = d2lambda(pvmbs.get_d());
      for(int n=1,i=3; n<=nfib; n++,i+=nparams_per_fibre){
	fs(n)=pvmbs.get_f(n);
	float tmpr=fs(n)/partial_fsum(fs,n-1);
	if (tmpr>1) tmpr=1; //This can be true due to numerical errors
	start(i) = f2beta(tmpr); 
	start(i+1) = pvmbs.get_th(n);
	start(i+2) = pvmbs.get_ph(n);
	start(i+3) = pvmbs.get_fanning_angle(n);
	start(i+4) = k12l1(k1_init(n1));
	start(i+5) = w2gam(w_init(n2)); 
      } 
      if (m_include_f0){
	float tmpr=pvmbs.get_f0()/partial_fsum(fs,nfib);
	if (tmpr>1) tmpr=1; //This can be true due to numerical errors
	start(nparams)=f2beta(tmpr);
      }

      // do the fit
      NonlinParam  lmpar(start.Nrows(),NL_LM); 
      lmpar.SetGaussNewtonType(LM_LM);
      lmpar.SetStartingEstimate(start);
  
      //lmpar.LogCF(true);
 
      NonlinOut status;
      status = nonlin(lmpar,(*this));
      ColumnVector tmp_par(nparams);
      tmp_par = lmpar.Par();
      
      //cout<<"Number of Iterations: "<<lmpar.NIter()<<endl;
      //vector<double> Cf=lmpar.CFHistory();
      //for (int n=0; n<(int)Cf.size(); n++)
      //  cout<<Cf[n]<<" ";
      //cout<<endl;
      double RSS=cf(tmp_par); //get the sum of squared residuals
      if (RSS<=minRSS){
	final_par=tmp_par;
	minRSS=RSS;
      }
    }
 
  if (m_eval_BIC){ 
    m_BIC=npts*log(minRSS/npts)+log(npts)*nparams; //evaluate BIC
    //cout<<"RSS="<<RSS<<". BIC="<<m_BIC<<endl;
  }

  // finalise parameters
  m_s0 = final_par(1);
  m_d  = lambda2d(final_par(2));
  for(int n=1; n<=nfib; n++){
    int kk=3+nparams_per_fibre*(n-1);

    m_f(n)  = beta2f(final_par(kk))*partial_fsum(m_f,n-1);
    m_th(n) = final_par(kk+1);
    m_ph(n) = final_par(kk+2);
    m_psi(n)= final_par(kk+3);
    m_k1(n) = l12k1(final_par(kk+4));
    m_k2(n) = m_k1(n)/gam2w(final_par(kk+5));
  }

  if (m_include_f0)
    m_f0=beta2f(final_par(nparams))*partial_fsum(m_f,nfib);
  
  sort(); 
}


void PVM_Ball_Binghams::sort(){
  vector< pair<float,int> > fvals(nfib);
  ColumnVector ftmp(nfib),thtmp(nfib),phtmp(nfib),psitmp(nfib),k1tmp(nfib),k2tmp(nfib);
  ftmp=m_f;thtmp=m_th;phtmp=m_ph; psitmp=m_psi; k1tmp=m_k1; k2tmp=m_k2; 
  
  for(int i=1;i<=nfib;i++){
    pair<float,int> p(m_f(i),i);
    fvals[i-1] = p;
  }
  std::sort(fvals.begin(),fvals.end());
  for(int i=1,ii=nfib-1;ii>=0;i++,ii--){
    m_f(i)  = ftmp(fvals[ii].second);
    m_th(i) = thtmp(fvals[ii].second);
    m_ph(i) = phtmp(fvals[ii].second);
    m_psi(i)= psitmp(fvals[ii].second);
    m_k1(i)=  k1tmp(fvals[ii].second);
    m_k2(i)=  k2tmp(fvals[ii].second);
  }
}


//Returns 1-Sum(f_j), 1<=j<=ii. (ii<=nfib)
//Used for transforming beta to f and vice versa
float PVM_Ball_Binghams::partial_fsum(ColumnVector& fs, int ii) const{
  float fsum=1.0;
  for(int j=1;j<=ii;j++)
    fsum-=fs(j);
    if (fsum==0) //Very rare cases
      fsum=tiny;
  return fsum;
}



//Print the final estimates (after having them transformed)
void PVM_Ball_Binghams::print()const{
  cout << endl<<"Ball & Bingham FIT RESULTS " << endl;
  cout << "S0   :" << m_s0 << endl;
  cout << "D    :" << m_d << endl;
  for(int i=1;i<=nfib;i++){
    cout << "F" << i << "   :" << m_f(i) << endl;
    ColumnVector x(3),fan_vec;
    x << sin(m_th(i))*cos(m_ph(i)) << sin(m_th(i))*sin(m_ph(i)) << cos(m_th(i));
    fan_vec=get_fanning_vector(i);
    float _th,_ph,_psi;cart2sph(x,_th,_ph); _psi=m_psi(i);
    if(x(3)<0) {x=-x; _psi=M_PI-m_psi(i); }
    cout << "TH"  << i << " : " << _th*180.0/M_PI << " deg" << endl; 
    cout << "PH"  << i << " : " << _ph*180.0/M_PI << " deg" << endl; 
    cout << "PSI" << i << " : " <<_psi<<endl;
    cout << "DIR" << i << " : " << x(1) << " " << x(2) << " " << x(3) << endl;
    cout << "FAN_DIR" << i << " : " << fan_vec(1) << " " << fan_vec(2) << " " << fan_vec(3) << endl;
    cout << "K1_" << i << " : " <<m_k1(i)<<endl;
    cout << "K2_" << i << " : " <<m_k2(i)<<endl;
  }
  if (m_include_f0)
    cout << "F0    :" << m_f0 << endl;
  if (m_eval_BIC)
    cout<< "BIC  :"<<m_BIC<<endl;
}



//Print the estimates using a vector that contains the transformed parameter values
//i.e. need to untransform them to get d,f's etc
void PVM_Ball_Binghams::print(const ColumnVector& p)const{
  ColumnVector f(nfib);
  
  cout << "PARAMETER VALUES " << endl;
  cout << "S0   :" << p(1) << endl;
  cout << "D    :" << lambda2d(p(2)) << endl;
  for(int i=3,ii=1;ii<=nfib;i+=3,ii++){
    f(ii) = beta2f(p(i))*partial_fsum(f,ii-1);
    float k1=l12k1(p(i+4));
    cout << "F" << ii << "   :" << f(ii) << endl;
    cout << "TH" << ii << "  :" << p(i+1)*180.0/M_PI << " deg" << endl; 
    cout << "PH" << ii << "  :" << p(i+2)*180.0/M_PI << " deg" << endl; 
    cout << "PSI" << ii << "  :"<< p(i+3)*180.0/M_PI << " deg" << endl; 
    cout << "k1_" << ii << "  :"<< k1 << endl; 
    cout << "k2_" << ii << "  :"<< k1/gam2w(p(i+4))<< endl; 
  }
  if (m_include_f0)
    cout << "F0    :" << beta2f(p(nparams))*partial_fsum(f,nfib);
}



//Applies the forward model and gets the model predicted signal using the estimated parameter values  (true,non-transformed space)  
ReturnMatrix PVM_Ball_Binghams::get_prediction()const{
  ColumnVector pred(npts);
  ColumnVector p(nparams);
  ColumnVector fs(nfib);
  
  fs=m_f;
  p(1) = m_s0;   //Transform parameters to the space where they are uncostrained
  p(2) = d2lambda(m_d);
  for(int i=3,ii=1;ii<=nfib;i+=nparams_per_fibre,ii++){
    float tmpr=m_f(ii)/partial_fsum(fs,ii-1);
    if (tmpr>1.0) tmpr=1; //This can be due to numerical errors
    p(i)   = f2beta(tmpr);
    p(i+1) = m_th(ii);
    p(i+2) = m_ph(ii);
    p(i+3) = m_psi(ii);
    p(i+4) = k12l1(m_k1(ii));
    p(i+5) = w2gam(m_k1(ii)/m_k2(ii));
  }
  if (m_include_f0){
    float tmpr=m_f0/partial_fsum(fs,nfib);
    if (tmpr>1.0) tmpr=1; //This can be due to numerical errors
    p(nparams)=f2beta(tmpr);
  }

  pred = forwardModel(p);
  pred.Release();
  return pred;
} 



//Applies the forward model and gets a model predicted signal using the parameter values in p (transformed parameter space)  
NEWMAT::ReturnMatrix PVM_Ball_Binghams::forwardModel(const NEWMAT::ColumnVector& p)const{
  ColumnVector pred(npts);
  pred = 0;
  float val;
  float _d = lambda2d(p(2));
  ////////////////////////////////////
  ColumnVector fs(nfib);  ColumnVector temp_vec(3), denom(3);
  Matrix Rpsi(3,3), Rth(3,3), Rph(3,3), R(3,3);
  vector<Matrix> B; vector<ColumnVector> approx_denomB;
  DiagonalMatrix L(3);       SymmetricMatrix Q(3);
  L=0; Rpsi=0; Rth=0; Rph=0; Rth(2,2)=1; Rph(3,3)=1; Rpsi(3,3)=1;
  float sumf=0; fs=0;
  for(int k=1;k<=nfib;k++){
    int kk = 3+nparams_per_fibre*(k-1);
    fs(k) = beta2f(p(kk))*partial_fsum(fs,k-1);
    sumf += fs(k);
    float cosph, sinph,cospsi,sinpsi,costh,sinth,k1,k2;
    costh=cos(p(kk+1)); sinth=sin(p(kk+1)); cosph=cos(p(kk+2)); sinph=sin(p(kk+2));
    cospsi=cos(p(kk+3)); sinpsi=sin(p(kk+3)); 
    k1=l12k1(p(kk+4)); k2=k1/gam2w(p(kk+5));
    L(1)=-k1; L(2)=-k2; denom<<L(1)<<L(2)<<0;
    Rth(1,1)=costh; Rth(1,3)=-sinth; Rth(3,1)=sinth; Rth(3,3)=costh;
    Rph(1,1)=cosph; Rph(1,2)=sinph; Rph(2,1)=-sinph; Rph(2,2)=cosph;
    Rpsi(1,1)=cospsi; Rpsi(1,2)=sinpsi; Rpsi(2,1)=-sinpsi; Rpsi(2,2)=cospsi;
    R=Rpsi*Rth*Rph;
    R<<R.t()*L*R;
    B.push_back(R);
    temp_vec=approx_denominatorB(denom);
    approx_denomB.push_back(temp_vec);
  }
  ////////////////////////////////////
  for(int i=1;i<=Y.Nrows();i++){
    val = 0.0;
    for(int k=1;k<=nfib;k++){
      Q<<B[k-1]-_d*bvecs_dyadic[i-1];
      EigenValues(Q,L);  temp_vec<<L(1)<<L(2)<<L(3);
      val += fs(k)*hyp_SratioB_knowndenom(temp_vec,approx_denomB[k-1]);
    }
    if (m_include_f0){
      float temp_f0=beta2f(p(nparams))*partial_fsum(fs,nfib);
      pred(i) = p(1)*(temp_f0+(1-sumf-temp_f0)*isoterm(i,_d)+val);
    } 
    else
      pred(i) = p(1)*((1-sumf)*isoterm(i,_d)+val); 
  }
  pred.Release();
  return pred;
}


//Instead of returning the model predicted signal for each direction
//returns the individual signal contributions i.e. isotropic, anisotropic1, anisotropic2,etc. 
//Weighting with the fractions, scaling with S0 and summing those gives the signal.
//A Matrix npts x (nfib+1) is returned 
NEWMAT::ReturnMatrix PVM_Ball_Binghams::forwardModel_compartments(const NEWMAT::ColumnVector& p) const{
  Matrix Sig(npts,nfib+1);

  float _d = lambda2d(p(2));
  ////////////////////////////////////
  ColumnVector fs(nfib);  ColumnVector temp_vec(3), denom(3);
  Matrix Rpsi(3,3), Rth(3,3), Rph(3,3), R(3,3);
  vector<Matrix> B; vector<ColumnVector> approx_denomB;
  DiagonalMatrix L(3);       SymmetricMatrix Q(3);
  L=0; Rpsi=0; Rth=0; Rph=0; Rth(2,2)=1; Rph(3,3)=1; Rpsi(3,3)=1;
  
  float sumf=0; fs=0;
  for(int k=1;k<=nfib;k++){
    int kk = 3+nparams_per_fibre*(k-1);
    fs(k) = beta2f(p(kk))*partial_fsum(fs,k-1);
    sumf += fs(k);
    float cosph, sinph,cospsi,sinpsi,costh,sinth,k1,k2;
    costh=cos(p(kk+1)); sinth=sin(p(kk+1)); cosph=cos(p(kk+2)); sinph=sin(p(kk+2));
    cospsi=cos(p(kk+3)); sinpsi=sin(p(kk+3)); 
    k1=l12k1(p(kk+4)); k2=k1/gam2w(p(kk+5));
    L(1)=-k1; L(2)=-k2; denom<<L(1)<<L(2)<<0;
    Rth(1,1)=costh; Rth(1,3)=-sinth; Rth(3,1)=sinth; Rth(3,3)=costh;
    Rph(1,1)=cosph; Rph(1,2)=sinph; Rph(2,1)=-sinph; Rph(2,2)=cosph;
    Rpsi(1,1)=cospsi; Rpsi(1,2)=sinpsi; Rpsi(2,1)=-sinpsi; Rpsi(2,2)=cospsi;
    R=Rpsi*Rth*Rph;
    R<<R.t()*L*R;
    B.push_back(R);
    temp_vec=approx_denominatorB(denom);
    approx_denomB.push_back(temp_vec);
  }
  ////////////////////////////////////
  for(int i=1;i<=Y.Nrows();i++){
      Sig(i,1) = isoterm(i,_d);
 
    for(int k=1;k<=nfib;k++){
      Q<<B[k-1]-_d*bvecs_dyadic[i-1];
      EigenValues(Q,L);  temp_vec<<L(1)<<L(2)<<L(3);
      Sig(i,k+1)= hyp_SratioB_knowndenom(temp_vec,approx_denomB[k-1]);
    }
  }  

  Sig.Release();
  return Sig;
}


//Builds up the model predicted signal for each direction by using precomputed individual compartment signals, stored in Matrix Sig. 
//Weights them with the fractions, scales with S0 and sums to get the signal.
NEWMAT::ReturnMatrix PVM_Ball_Binghams::pred_from_compartments(const NEWMAT::ColumnVector& p, const NEWMAT::Matrix& Sig) const{
  ColumnVector pred(npts);
  float val;
  ColumnVector fs(nfib);

  float sumf=0; fs=0;
  for(int k=1;k<=nfib;k++){
    int kk = 3+nparams_per_fibre*(k-1);
    fs(k) = beta2f(p(kk))*partial_fsum(fs,k-1);
    sumf += fs(k);
  }
  
  for(int i=1;i<=Y.Nrows();i++){
    val = 0.0;
    for(int k=1;k<=nfib;k++)
      val += fs(k)*Sig(i,k+1);
    
    if (m_include_f0){
      float temp_f0=beta2f(p(nparams))*partial_fsum(fs,nfib);
      pred(i) = p(1)*(temp_f0+(1-sumf-temp_f0)*Sig(i,1)+val);
    } 
    else
      pred(i) = p(1)*((1-sumf)*Sig(i,1)+val); 
  }
  
  pred.Release();
  return pred;
}



//Builds up the model predicted signal for each direction by using precomputed individual compartment signals, stored in Matrix Sig. 
//Weights them with the fractions, scales with S0 and sums to get the signal.
//The signal of the fibre compartment with index fib is recalculated.
NEWMAT::ReturnMatrix PVM_Ball_Binghams::pred_from_compartments(const NEWMAT::ColumnVector& p, const NEWMAT::Matrix& Sig,const int& fib) const{
  ColumnVector pred(npts); Matrix newSig;
  float val;

  float _d = lambda2d(p(2));
  ////////////////////////////////////
  ColumnVector fs(nfib);  ColumnVector temp_vec(3), denom(3);
  Matrix Rpsi(3,3), Rth(3,3), Rph(3,3), R(3,3);
  DiagonalMatrix L(3);       SymmetricMatrix Q(3);
  L=0; Rpsi=0; Rth=0; Rph=0; Rth(2,2)=1; Rph(3,3)=1; Rpsi(3,3)=1;
  
  float sumf=0; fs=0;
  for(int k=1;k<=nfib;k++){
    int kk = 3+nparams_per_fibre*(k-1);
    fs(k) = beta2f(p(kk))*partial_fsum(fs,k-1);
    sumf += fs(k);
  }
  ///////////////////////////////////////
  int kk = 3+nparams_per_fibre*(fib-1);  float cosph, sinph,cospsi,sinpsi,costh,sinth,k1,k2;
  costh=cos(p(kk+1)); sinth=sin(p(kk+1)); cosph=cos(p(kk+2)); sinph=sin(p(kk+2));
  cospsi=cos(p(kk+3)); sinpsi=sin(p(kk+3)); 
  k1=l12k1(p(kk+4)); k2=k1/gam2w(p(kk+5));
  L(1)=-k1; L(2)=-k2; denom<<L(1)<<L(2)<<0;
  Rth(1,1)=costh; Rth(1,3)=-sinth; Rth(3,1)=sinth; Rth(3,3)=costh;
  Rph(1,1)=cosph; Rph(1,2)=sinph; Rph(2,1)=-sinph; Rph(2,2)=cosph;
  Rpsi(1,1)=cospsi; Rpsi(1,2)=sinpsi; Rpsi(2,1)=-sinpsi; Rpsi(2,2)=cospsi;
  R=Rpsi*Rth*Rph;
  R<<R.t()*L*R;
  temp_vec=approx_denominatorB(denom);
  denom=temp_vec;
  
  newSig=Sig;  //Get the new Signal for compartment fib
  for(int i=1;i<=Y.Nrows();i++){
    Q<<R-_d*bvecs_dyadic[i-1];
    EigenValues(Q,L);  temp_vec<<L(1)<<L(2)<<L(3);
    newSig(i,fib+1)=hyp_SratioB_knowndenom(temp_vec,denom);
  }  
  ///////////////////////////////////////

  for(int i=1;i<=Y.Nrows();i++){
    val = 0.0;
    for(int k=1;k<=nfib;k++)
      val += fs(k)*newSig(i,k+1);
    
    if (m_include_f0){
      float temp_f0=beta2f(p(nparams))*partial_fsum(fs,nfib);
      pred(i) = p(1)*(temp_f0+(1-sumf-temp_f0)*newSig(i,1)+val);
    } 
    else
      pred(i) = p(1)*((1-sumf)*newSig(i,1)+val); 
  }
  
  pred.Release();
  return pred;
}



//Cost Function, sum of squared residuals
//assume that parameter values p are transformed (e.g. need to untransform them to get d, f's,etc)
double PVM_Ball_Binghams::cf(const NEWMAT::ColumnVector& p)const{
  double cfv = 0.0;
  double err;
  ColumnVector S;

  S=forwardModel(p);  //Model predictions
  for(int i=1;i<=npts;i++){
    err=S(i)-Y(i);    //Residual
    cfv+=err*err;     //Sum of squared residuals
  }  
  //cout<<"CF="<<cfv<<endl; OUT(p.t());
  return(cfv);
}

/*

//Slower implementation
NEWMAT::ReturnMatrix PVM_Ball_Binghams::grad(const ColumnVector& p)const
{
  ColumnVector gradv(nparams);//This is the gradient of the cost function
  Matrix J(npts,nparams);     //This is the Jacobian matrix of the model equation
                              //The derivative of the cost function w.r.t. parameter j 
                              //will then be: Grad_j=Sum(2*(F(x_i)-Y_i)*J(i,j)), with Sum across data points i
  ColumnVector diff(npts);    //Residuals
  ColumnVector p_plus_h, S_trial,S;

  //Compute the Jacobian first using finite differences for each element
  double step;
  ColumnVector typical_scale(nparams); 
  typical_scale=1;
  for (int k=1; k<=nfib; k++){
    int kk = 3+nparams_per_fibre*(k-1); 
    typical_scale(kk+4)=100;
    typical_scale(kk+5)=1;
  }
  
  S=forwardModel(p);
  diff=S-Y;
  for (int i=1; i<=npts; i++)
    J(i,1)=S(i)/p(1);  //derivatives with respect to S0 are analytic: S_i/S0
    
  for (int n=2; n<=nparams; n++) {
    p_plus_h = p;
    step = SQRTtiny*nonzerosign(p(n))*max(fabs(p(n))*typical_scale(n),1.0);
    step = nonzerosign(tiny)*min(max(fabs(step),1.0e-8),0.1);  //check that 1e-8<step<0.1

    p_plus_h(n)=p(n)+step;
    S_trial=forwardModel(p_plus_h);
    for (int i=1; i<=npts; i++)
      J(i,n)=(S_trial(i)-S(i))/step;
  } 
  
  gradv =2*J.t()*diff;
  gradv.Release();
  return(gradv);    
}


//Slower implementation
//this uses Gauss-Newton approximation
boost::shared_ptr<BFMatrix> PVM_Ball_Binghams::hess(const NEWMAT::ColumnVector& p,boost::shared_ptr<BFMatrix> iptr)const{
  boost::shared_ptr<BFMatrix>   hessm;
  if (iptr && iptr->Nrows()==(unsigned int)p.Nrows() && iptr->Ncols()==(unsigned int)p.Nrows()) hessm = iptr;
  else hessm = boost::shared_ptr<BFMatrix>(new FullBFMatrix(p.Nrows(),p.Nrows()));
  
  Matrix J(npts,nparams);     //This is the Jacobian matrix of the model equation
  ColumnVector p_plus_h, S_trial,S;
  
  //Compute the Jacobian first using finite differences for each element
  double step,sig;
  ColumnVector typical_scale(nparams); 
  typical_scale=1;
  for (int k=1; k<=nfib; k++){
    int kk = 3+nparams_per_fibre*(k-1); 
    typical_scale(kk+4)=100;
    typical_scale(kk+5)=1;
  }

  S=forwardModel(p);

  for (int i=1; i<=npts; i++)
    J(i,1)=S(i)/p(1);  //derivatives with respect to S0 are analytic: S_i/S0
    
  for (int n=2; n<=nparams; n++) {
    p_plus_h = p;
    step = SQRTtiny*nonzerosign(p(n))*max(fabs(p(n))*typical_scale(n),1.0);
    step = nonzerosign(tiny)*min(max(fabs(step),1.0e-8),0.1);  //check that 1e-8<step<0.1

    p_plus_h(n)=p(n)+step;
    S_trial=forwardModel(p_plus_h);
    for (int i=1; i<=npts; i++)
      J(i,n)=(S_trial(i)-S(i))/step;
  } 
  
  for (int i=1; i<=p.Nrows(); i++){
    for (int j=i; j<=p.Nrows(); j++){
      sig = 0.0;
      for(int k=1;k<=J.Nrows();k++)
	sig += 2*(J(k,i)*J(k,j));
      hessm->Set(i,j,sig);
    }
  }
  for (int j=1; j<=p.Nrows(); j++) {
    for (int i=j+1; i<=p.Nrows(); i++) {
      hessm->Set(i,j,hessm->Peek(j,i));
    }
  }
  return(hessm);
}

*/

NEWMAT::ReturnMatrix PVM_Ball_Binghams::grad(const ColumnVector& p)const
{
  ColumnVector gradv(nparams);//This is the gradient of the cost function
  Matrix J(npts,nparams);     //This is the Jacobian matrix of the model equation
                              //The derivative of the cost function w.r.t. parameter j 
                              //will then be: Grad_j=Sum(2*(F(x_i)-Y_i)*J(i,j)), with Sum across data points i
  ColumnVector diff(npts);    //Residuals
  ColumnVector p_plus_h, S_trial,S;
  Matrix Sig;
  ColumnVector fs(nfib), bs(nfib);

  //Compute the Jacobian first using finite differences for each element. Derivatives are analytic for S0 and the volume fractions
  double step;
  ColumnVector typical_scale(nparams); 
  typical_scale=1;
  for (int k=1; k<=nfib; k++){
    int kk = 3+nparams_per_fibre*(k-1); 
    bs(k)=p(kk);
    fs(k) = beta2f(p(kk))*partial_fsum(fs,k-1);
    typical_scale(kk+4)=100;
    typical_scale(kk+5)=1;
  }

  //Compute the derivatives with respect to betas, i.e the transformed volume fraction variables
  Matrix f_deriv;
  f_deriv=fractions_deriv(nfib, fs, bs);  

  Sig=forwardModel_compartments(p);
  S=pred_from_compartments(p, Sig);
  diff=S-Y;
  for (int i=1; i<=npts; i++)
    J(i,1)=S(i)/p(1);  //derivatives with respect to S0 are analytic: S_i/S0
    
  for (int n=2; n<=nparams; n++) {
    p_plus_h = p;
    step = SQRTtiny*nonzerosign(p(n))*max(fabs(p(n))*typical_scale(n),1.0);
    step = nonzerosign(tiny)*min(max(fabs(step),1.0e-8),0.1);  //check that 1e-8<step<0.1
    p_plus_h(n)=p(n)+step;
      
    if (n<3 || (n==nparams && m_include_f0)){ //for d all compartments will change. Also for f0, if included. Use for those params numerical differentiation
      S_trial=forwardModel(p_plus_h);
      for (int i=1; i<=npts; i++)
	J(i,n)=(S_trial(i)-S(i))/step; 
    }
    else{  //for the other params, update only the signal of the relevant fibre compartment
      int fib_indx=(int)ceil((float)(n-2)/(float)nparams_per_fibre); //index indicating in which fibre compartment parameter n belongs to
      if (n==3+nparams_per_fibre*(fib_indx-1)){ //Then we have a volume fraction, derivative is analytic.
	for (int i=1; i<=npts; i++){
	  J(i,n)=0;
	  for (int j=1; j<=nfib; j++){
	    if (f_deriv(j,fib_indx)!=0)
	      J(i,n) += p(1)*(Sig(i,j+1)-Sig(i,1))*f_deriv(j,fib_indx);
	  }
	}
      }
      else{       //for all other params, use numerical differentiation
	S_trial=pred_from_compartments(p_plus_h, Sig,fib_indx);
	for (int i=1; i<=npts; i++){
	  J(i,n)=(S_trial(i)-S(i))/step;
	  //if (J(i,n)==0) cout<<"Zero gradient!!"<<endl;//J(i,n)=1e-8;  //stabilize LM in case differentiation has failed due to step size
	}
      }
    }
  } 
  gradv =2*J.t()*diff;
  gradv.Release();
  return(gradv);    
}



//this uses Gauss-Newton approximation
boost::shared_ptr<BFMatrix> PVM_Ball_Binghams::hess(const NEWMAT::ColumnVector& p,boost::shared_ptr<BFMatrix> iptr)const{
  boost::shared_ptr<BFMatrix>   hessm;
  if (iptr && iptr->Nrows()==(unsigned int)p.Nrows() && iptr->Ncols()==(unsigned int)p.Nrows()) hessm = iptr;
  else hessm = boost::shared_ptr<BFMatrix>(new FullBFMatrix(p.Nrows(),p.Nrows()));
  
  Matrix J(npts,nparams);     //This is the Jacobian matrix of the model equation
  ColumnVector p_plus_h, S_trial,S, fs(nfib), bs(nfib);
  Matrix Sig;
  
  //Compute the Jacobian first using finite differences for each element. Derivatives are analytic for S0 and the volume fractions
  double step,sigt;
  ColumnVector typical_scale(nparams); 
  typical_scale=1;
  for (int k=1; k<=nfib; k++){
    int kk = 3+nparams_per_fibre*(k-1); 
    bs(k)=p(kk);
    fs(k) = beta2f(p(kk))*partial_fsum(fs,k-1);
    typical_scale(kk+4)=100;
    typical_scale(kk+5)=1;
  }

  //Compute the derivatives with respect to betas, i.e the transformed volume fraction variables
  Matrix f_deriv;
  f_deriv=fractions_deriv(nfib, fs, bs);  

  Sig=forwardModel_compartments(p);
  S=pred_from_compartments(p, Sig);

  for (int i=1; i<=npts; i++)
    J(i,1)=S(i)/p(1);  //derivatives with respect to S0 are analytic: S_i/S0
    
  for (int n=2; n<=nparams; n++) {
    p_plus_h = p;
    step = SQRTtiny*nonzerosign(p(n))*max(fabs(p(n))*typical_scale(n),1.0);
    step = nonzerosign(tiny)*min(max(fabs(step),1.0e-8),0.1);  //check that 1e-8<step<0.1
    p_plus_h(n)=p(n)+step;
      
    if (n<3 || (n==nparams && m_include_f0)){ //for d all compartments will change. Also for f0, if included. Use for those params numerical differentiation
      S_trial=forwardModel(p_plus_h);
      for (int i=1; i<=npts; i++)
	J(i,n)=(S_trial(i)-S(i))/step; 
    }
    else{  //for the other params, update only the signal of the relevant fibre compartment
      int fib_indx=(int)ceil((float)(n-2)/(float)nparams_per_fibre); //index indicating in which fibre compartment parameter n belongs to
      if (n==3+nparams_per_fibre*(fib_indx-1)){ //Then we have a volume fraction, derivative is analytic.
	for (int i=1; i<=npts; i++){
	  J(i,n)=0;
	  for (int j=1; j<=nfib; j++){
	    if (f_deriv(j,fib_indx)!=0)
	      J(i,n) += p(1)*(Sig(i,j+1)-Sig(i,1))*f_deriv(j,fib_indx); 
	  }
	}
      }
      else{       //for all other params, use numerical differentiation
	S_trial=pred_from_compartments(p_plus_h, Sig,fib_indx);
	for (int i=1; i<=npts; i++){
	  J(i,n)=(S_trial(i)-S(i))/step;
	  //if (J(i,n)==0) J(i,n)=1e-8; //stabilize LM in case differentiation has failed due to step size
	}
      }
    }
  } 
  
  for (int i=1; i<=p.Nrows(); i++){
    for (int j=i; j<=p.Nrows(); j++){
      sigt = 0.0;
      for(int k=1;k<=J.Nrows();k++)
	sigt += 2*(J(k,i)*J(k,j));
      hessm->Set(i,j,sigt);
    }
  }
  for (int j=1; j<=p.Nrows(); j++) {
    for (int i=j+1; i<=p.Nrows(); i++) {
      hessm->Set(i,j,hessm->Peek(j,i));
    }
  }
  return(hessm);
}



float PVM_Ball_Binghams::isoterm(const int& pt,const float& _d)const{
  return(std::exp(-bvals(1,pt)*_d));
}


NEWMAT::ReturnMatrix PVM_Ball_Binghams::fractions_deriv(const int& nfib, const ColumnVector& fs, const ColumnVector& bs) const{
  NEWMAT::Matrix Deriv(nfib,nfib);
  float fsum;
  Deriv=0;
  for (int j=1; j<=nfib; j++)
    for (int k=1; k<=nfib; k++){
      if (j==k){
	fsum=1; 
	for (int n=1; n<=j-1; n++)
	  fsum-=fs(n);
	Deriv(j,k)=sin(2*bs(k))*fsum;
      }
      else if (j>k){
	fsum=0;
	for (int n=1; n<=j-1; n++)
	  fsum+=Deriv(n,k);
	Deriv(j,k)=-pow(sin(bs(j)),2.0)*fsum;  
      }
    } 	   
  Deriv.Release();
  return Deriv;
}



//Returns a vector that indicates the fanning orientation
NEWMAT::ReturnMatrix PVM_Ball_Binghams:: get_fanning_vector(const int& i) const{ 
  ColumnVector fan_vec(3); 
  float t_th=m_th(i);  float t_ph=m_ph(i);  float t_psi=m_psi(i);

  float costh=cos(t_th); float sinth=sin(t_th); 
  float cosph=cos(t_ph); float sinph=sin(t_ph);
  float cospsi=cos(t_psi); float sinpsi=sin(t_psi);
  /*
  Matrix Rpsi(3,3), Rth(3,3), Rph(3,3), R(3,3);
  Rpsi=0; Rth=0; Rph=0; Rth(2,2)=1; Rph(3,3)=1; Rpsi(3,3)=1;
  Rth(1,1)=costh; Rth(1,3)=-sinth; Rth(3,1)=sinth; Rth(3,3)=costh;
  Rph(1,1)=cosph; Rph(1,2)=sinph; Rph(2,1)=-sinph; Rph(2,2)=cosph;
  Rpsi(1,1)=cospsi; Rpsi(1,2)=sinpsi; Rpsi(2,1)=-sinpsi; Rpsi(2,2)=cospsi;
  R=Rpsi*Rth*Rph; */
  
  //fan_vec(1)=R(2,1); fan_vec(2)=R(2,2); fan_vec(3)=R(2,3);
  fan_vec(1)=-sinpsi*costh*cosph-cospsi*sinph; fan_vec(2)=-sinpsi*costh*sinph+cospsi*cosph; fan_vec(3)=sinpsi*sinth;

  fan_vec.Release();
  return fan_vec;
}




///////////////////////////////////////////////////////////////////////////
// FANNING MODEL - BALL & WATSONS 
// Constrained Optimization for the diffusivity, fractions and their sum<1,
// and the Bingham eigenvalues
//////////////////////////////////////////////////////////////////////////

void PVM_Ball_Watsons::fit(){
  // Fit the ball & stick first to initialize some of the parameters
  PVM_single_c pvmbs(Y,bvecs,bvals,nfib,false,m_include_f0);  
  pvmbs.fit();
  //  pvmbs.print();
  
  ColumnVector k_init;
  ColumnVector final_par(nparams);
  double minRSS=1e20;
  if (!m_gridsearch){
    k_init.ReSize(1); k_init<<20;
  }
  else{
    k_init.ReSize(6); k_init<< 10 << 20 << 50 << 100 << 500 << 1000;
  }

  for (int n1=1; n1<=k_init.Nrows(); n1++){
    ColumnVector start(nparams);
    ColumnVector fs(nfib); fs=0;
    //Initialize the non-linear fitter. Transform all initial values to the uncostrained parameter space 
    start(1) = pvmbs.get_s0();
    start(2) = d2lambda(pvmbs.get_d());
    for(int n=1,i=3; n<=nfib; n++,i+=nparams_per_fibre){
      fs(n)=pvmbs.get_f(n);
      float tmpr=fs(n)/partial_fsum(fs,n-1);
      if (tmpr>1) tmpr=1; //This can be true due to numerical errors
      start(i) = f2beta(tmpr); 
      start(i+1) = pvmbs.get_th(n);
      start(i+2) = pvmbs.get_ph(n);
      start(i+3) = k12l1(k_init(n1));
    } 
   
    if (m_include_f0){
      float tmpr=pvmbs.get_f0()/partial_fsum(fs,nfib);
      if (tmpr>1) tmpr=1; //This can be true due to numerical errors
      start(nparams)=f2beta(tmpr);
    }

    // do the fit
    NonlinParam  lmpar(start.Nrows(),NL_LM); 
    lmpar.SetGaussNewtonType(LM_LM);
    lmpar.SetStartingEstimate(start);
  
    //lmpar.LogCF(true);
 
    NonlinOut status;
    status = nonlin(lmpar,(*this));
    ColumnVector tmp_par(nparams);
    tmp_par = lmpar.Par();

    /*cout<<"Number of Iterations: "<<lmpar.NIter()<<endl;
      vector<double> Cf=lmpar.CFHistory();
      for (int n=0; n<(int)Cf.size(); n++)
      cout<<Cf[n]<<" ";
      cout<<endl;
    */
    double RSS=cf(tmp_par); //get the sum of squared residuals
    if (RSS<=minRSS){
      final_par=tmp_par;
      minRSS=RSS;
    }
  }
  
  if (m_eval_BIC){ 
    m_BIC=npts*log(minRSS/npts)+log(npts)*nparams; //evaluate BIC
  }
 
  // finalise parameters
  m_s0 = final_par(1);
  m_d  = lambda2d(final_par(2));
  for(int n=1; n<=nfib; n++){
    int kk=3+nparams_per_fibre*(n-1);

    m_f(n)  = beta2f(final_par(kk))*partial_fsum(m_f,n-1);
    m_th(n) = final_par(kk+1);
    m_ph(n) = final_par(kk+2);
    m_k(n) = l12k1(final_par(kk+3));
  }

  if (m_include_f0)
    m_f0=beta2f(final_par(nparams))*partial_fsum(m_f,nfib);
  
  sort();
}


void PVM_Ball_Watsons::sort(){
  vector< pair<float,int> > fvals(nfib);
  ColumnVector ftmp(nfib),thtmp(nfib),phtmp(nfib),ktmp(nfib);
  ftmp=m_f;thtmp=m_th;phtmp=m_ph; ktmp=m_k; 
  
  for(int i=1;i<=nfib;i++){
    pair<float,int> p(m_f(i),i);
    fvals[i-1] = p;
  }
  std::sort(fvals.begin(),fvals.end());
  for(int i=1,ii=nfib-1;ii>=0;i++,ii--){
    m_f(i)  = ftmp(fvals[ii].second);
    m_th(i) = thtmp(fvals[ii].second);
    m_ph(i) = phtmp(fvals[ii].second);
    m_k(i)=  ktmp(fvals[ii].second);
  }
}


//Returns 1-Sum(f_j), 1<=j<=ii. (ii<=nfib)
//Used for transforming beta to f and vice versa
float PVM_Ball_Watsons::partial_fsum(ColumnVector& fs, int ii) const{
  float fsum=1.0;
  for(int j=1;j<=ii;j++)
    fsum-=fs(j);
  if (fsum==0) //Very rare cases
    fsum=tiny;
  return fsum;
}



//Print the final estimates (after having them transformed)
void PVM_Ball_Watsons::print()const{
  cout << endl<<"Ball & Watson FIT RESULTS " << endl;
  cout << "S0   :" << m_s0 << endl;
  cout << "D    :" << m_d << endl;
  for(int i=1;i<=nfib;i++){
    cout << "F" << i << "   :" << m_f(i) << endl;
    ColumnVector x(3);
    x << sin(m_th(i))*cos(m_ph(i)) << sin(m_th(i))*sin(m_ph(i)) << cos(m_th(i));
    float _th,_ph;cart2sph(x,_th,_ph);
    if(x(3)<0) x=-x; 
    cout << "TH"  << i << " : " << _th*180.0/M_PI << " deg" << endl; 
    cout << "PH"  << i << " : " << _ph*180.0/M_PI << " deg" << endl; 
    cout << "DIR" << i << " : " << x(1) << " " << x(2) << " " << x(3) << endl;
    cout << "K_" << i << " : " <<m_k(i)<<endl;
  }
  if (m_include_f0)
    cout << "F0    :" << m_f0 << endl;
  if (m_eval_BIC)
    cout<< "BIC  :"<<m_BIC<<endl;
}



//Print the estimates using a vector that contains the transformed parameter values
//i.e. need to untransform them to get d,f's etc
void PVM_Ball_Watsons::print(const ColumnVector& p)const{
  ColumnVector f(nfib);
  
  cout << "PARAMETER VALUES " << endl;
  cout << "S0   :" << p(1) << endl;
  cout << "D    :" << lambda2d(p(2)) << endl;
  for(int i=3,ii=1;ii<=nfib;i+=3,ii++){
    f(ii) = beta2f(p(i))*partial_fsum(f,ii-1);
    float _k=l12k1(p(i+3));
    cout << "F" << ii << "   :" << f(ii) << endl;
    cout << "TH" << ii << "  :" << p(i+1)*180.0/M_PI << " deg" << endl; 
    cout << "PH" << ii << "  :" << p(i+2)*180.0/M_PI << " deg" << endl; 
    cout << "K_" << ii << "  :"<< _k << endl; 
  }
  if (m_include_f0)
    cout << "F0    :" << beta2f(p(nparams))*partial_fsum(f,nfib);
}



//Applies the forward model and gets the model predicted signal using the estimated parameter values  (true,non-transformed space)  
ReturnMatrix PVM_Ball_Watsons::get_prediction()const{
  ColumnVector pred(npts);
  ColumnVector p(nparams);
  ColumnVector fs(nfib);
  
  fs=m_f;
  p(1) = m_s0;   //Transform parameters to the space where they are uncostrained
  p(2) = d2lambda(m_d);
  for(int i=3,ii=1;ii<=nfib;i+=nparams_per_fibre,ii++){
    float tmpr=m_f(ii)/partial_fsum(fs,ii-1);
    if (tmpr>1.0) tmpr=1; //This can be due to numerical errors
    p(i)   = f2beta(tmpr);
    p(i+1) = m_th(ii);
    p(i+2) = m_ph(ii);
    p(i+3) = k12l1(m_k(ii));
  }
  if (m_include_f0){
    float tmpr=m_f0/partial_fsum(fs,nfib);
    if (tmpr>1.0) tmpr=1; //This can be due to numerical errors
    p(nparams)=f2beta(tmpr);
  }

  pred = forwardModel(p);
  pred.Release();
  return pred;
} 



//Applies the forward model and gets a model predicted signal using the parameter values in p (transformed parameter space)  
NEWMAT::ReturnMatrix PVM_Ball_Watsons::forwardModel(const NEWMAT::ColumnVector& p)const{
  ColumnVector pred(npts);
  pred = 0;
  float val;
  float _d = lambda2d(p(2));
  ////////////////////////////////////
  ColumnVector fs(nfib), ks(nfib);  ColumnVector temp_vec(3), denom(3);
  Matrix v(nfib,3); vector<ColumnVector> approx_denomW; Matrix A(3,3);
  float sumf=0; fs=0;
  
  for(int n=1;n<=nfib;n++){
    int nn = 3+nparams_per_fibre*(n-1);
    float cosph, sinph,costh,sinth;
    fs(n) = beta2f(p(nn))*partial_fsum(fs,n-1);
    sumf += fs(n);
    costh=cos(p(nn+1)); sinth=sin(p(nn+1)); cosph=cos(p(nn+2)); sinph=sin(p(nn+2));
    v(n,1) = cosph*sinth; v(n,2) = sinph*sinth; v(n,3) = costh;
    ks(n)=l12k1(p(nn+3)); 
    temp_vec=approx_denominatorW(ks(n));
    approx_denomW.push_back(temp_vec);
  }

  ////////////////////////////////////
  for(int i=1;i<=Y.Nrows();i++){
    val = 0.0;
    float bd=bvals(1,i)*_d;
    for(int n=1;n<=nfib;n++){
      A=(v.Row(n)).t()*(bvecs.Column(i)).t();
      float Q=-2*bd*ks(n)*(pow(A(1,1)+A(2,2)+A(3,3),2.0) - pow(A(2,1)-A(1,2),2.0) - pow(A(1,3)-A(3,1),2.0)-pow(A(2,3)-A(3,2),2.0)); 
      Q=sqrt(ks(n)*ks(n)+bd*bd+Q);
      temp_vec<<0.5*(ks(n)-bd+Q)<<0.5*(ks(n)-bd-Q)<<0;
      val+= fs(n)*hyp_SratioW_knowndenom(temp_vec,approx_denomW[n-1]);
    }
    if (m_include_f0){
      float temp_f0=beta2f(p(nparams))*partial_fsum(fs,nfib);
      pred(i) = p(1)*(temp_f0+(1-sumf-temp_f0)*isoterm(i,_d)+val);
    } 
    else
      pred(i) = p(1)*((1-sumf)*isoterm(i,_d)+val); 
  }  
  pred.Release();
  return pred;
}


//Instead of returning the model predicted signal for each direction
//returns the individual signal contributions i.e. isotropic, anisotropic1, anisotropic2,etc. 
//Weighting with the fractions, scaling with S0 and summing those gives the signal.
//A Matrix npts x (nfib+1) is returned 
NEWMAT::ReturnMatrix PVM_Ball_Watsons::forwardModel_compartments(const NEWMAT::ColumnVector& p) const{
  Matrix Sig(npts,nfib+1);

  float _d = lambda2d(p(2));
  ////////////////////////////////////
  ColumnVector ks(nfib);  ColumnVector temp_vec(3), denom(3);
  Matrix v(nfib,3); vector<ColumnVector> approx_denomW; Matrix A(3,3);
  
  for(int n=1;n<=nfib;n++){
    int nn = 3+nparams_per_fibre*(n-1);
    float cosph, sinph,costh,sinth;
    costh=cos(p(nn+1)); sinth=sin(p(nn+1)); cosph=cos(p(nn+2)); sinph=sin(p(nn+2));
    v(n,1) = cosph*sinth; v(n,2) = sinph*sinth; v(n,3) = costh;
    ks(n)=l12k1(p(nn+3)); 
    temp_vec=approx_denominatorW(ks(n));
    approx_denomW.push_back(temp_vec);
  }

  ////////////////////////////////////
  for(int i=1;i<=Y.Nrows();i++){
      Sig(i,1) = isoterm(i,_d);
      
      float bd=bvals(1,i)*_d;
      for(int n=1;n<=nfib;n++){
	A=(v.Row(n)).t()*(bvecs.Column(i)).t();
	float Q=-2*bd*ks(n)*(pow(A(1,1)+A(2,2)+A(3,3),2.0) - pow(A(2,1)-A(1,2),2.0) - pow(A(1,3)-A(3,1),2.0)-pow(A(2,3)-A(3,2),2.0)); 
	Q=sqrt(ks(n)*ks(n)+bd*bd+Q);
	temp_vec<<0.5*(ks(n)-bd+Q)<<0.5*(ks(n)-bd-Q)<<0;
	Sig(i,n+1)= hyp_SratioW_knowndenom(temp_vec,approx_denomW[n-1]);
      } 
  }  

  Sig.Release();
  return Sig;
}


//Builds up the model predicted signal for each direction by using precomputed individual compartment signals, stored in Matrix Sig. 
//Weights them with the fractions, scales with S0 and sums to get the signal.
NEWMAT::ReturnMatrix PVM_Ball_Watsons::pred_from_compartments(const NEWMAT::ColumnVector& p, const NEWMAT::Matrix& Sig) const{
  ColumnVector pred(npts);
  float val;
  ColumnVector fs(nfib);

  float sumf=0; fs=0;
  for(int n=1;n<=nfib;n++){
    int nn = 3+nparams_per_fibre*(n-1);
    fs(n) = beta2f(p(nn))*partial_fsum(fs,n-1);
    sumf += fs(n);
  }
  
  for(int i=1;i<=Y.Nrows();i++){
    val = 0.0;
    for(int n=1;n<=nfib;n++)
      val += fs(n)*Sig(i,n+1);
    
    if (m_include_f0){
      float temp_f0=beta2f(p(nparams))*partial_fsum(fs,nfib);
      pred(i) = p(1)*(temp_f0+(1-sumf-temp_f0)*Sig(i,1)+val);
    } 
    else
      pred(i) = p(1)*((1-sumf)*Sig(i,1)+val); 
  }
  
  pred.Release();
  return pred;
}



//Builds up the model predicted signal for each direction by using precomputed individual compartment signals, stored in Matrix Sig. 
//Weights them with the fractions, scales with S0 and sums to get the signal.
//The signal of the fibre compartment with index fib is recalculated.
NEWMAT::ReturnMatrix PVM_Ball_Watsons::pred_from_compartments(const NEWMAT::ColumnVector& p, const NEWMAT::Matrix& Sig,const int& fib) const{
  ColumnVector pred(npts); Matrix newSig;
  float val;

  float _d = lambda2d(p(2));
  ////////////////////////////////////
  ColumnVector fs(nfib);  ColumnVector temp_vec(3), denom(3), v(3); Matrix A(3,3);
  
  float sumf=0; fs=0;
  for(int n=1;n<=nfib;n++){
    int nn = 3+nparams_per_fibre*(n-1);
    fs(n) = beta2f(p(nn))*partial_fsum(fs,n-1);
    sumf += fs(n);
  }
  ///////////////////////////////////////
  int nn = 3+nparams_per_fibre*(fib-1);  float cosph, sinph,costh,sinth,_k;
  costh=cos(p(nn+1)); sinth=sin(p(nn+1)); cosph=cos(p(nn+2)); sinph=sin(p(nn+2));
  v(1) = cosph*sinth; v(2) = sinph*sinth; v(3) = costh;
  _k=l12k1(p(nn+3)); 
  temp_vec=approx_denominatorW(_k);
  denom=temp_vec;
  
  newSig=Sig;  //Get the new Signal for compartment fib
  for(int i=1;i<=Y.Nrows();i++){
    float bd=bvals(1,i)*_d;
    A=v*(bvecs.Column(i)).t();
    float Q=-2*bd*_k*(pow(A(1,1)+A(2,2)+A(3,3),2.0) - pow(A(2,1)-A(1,2),2.0) - pow(A(1,3)-A(3,1),2.0)-pow(A(2,3)-A(3,2),2.0)); 
    Q=sqrt(_k*_k+bd*bd+Q);
    temp_vec<<0.5*(_k-bd+Q)<<0.5*(_k-bd-Q)<<0;
    newSig(i,fib+1)= hyp_SratioW_knowndenom(temp_vec,denom);
  }  
  ///////////////////////////////////////

  for(int i=1;i<=Y.Nrows();i++){
    val = 0.0;
    for(int n=1;n<=nfib;n++)
      val += fs(n)*newSig(i,n+1);
    
    if (m_include_f0){
      float temp_f0=beta2f(p(nparams))*partial_fsum(fs,nfib);
      pred(i) = p(1)*(temp_f0+(1-sumf-temp_f0)*newSig(i,1)+val);
    } 
    else
      pred(i) = p(1)*((1-sumf)*newSig(i,1)+val); 
  }
  
  pred.Release();
  return pred;
}



//Cost Function, sum of squared residuals
//assume that parameter values p are transformed (e.g. need to untransform them to get d, f's,etc)
double PVM_Ball_Watsons::cf(const NEWMAT::ColumnVector& p)const{
  double cfv = 0.0;
  double err;
  ColumnVector S;

  S=forwardModel(p);  //Model predictions
  for(int i=1;i<=npts;i++){
    err=S(i)-Y(i);    //Residual
    cfv+=err*err;     //Sum of squared residuals
  }  
  //cout<<"CF="<<cfv<<endl; OUT(p.t());
  return(cfv);
}


NEWMAT::ReturnMatrix PVM_Ball_Watsons::grad(const ColumnVector& p)const
{
  ColumnVector gradv(nparams);//This is the gradient of the cost function
  Matrix J(npts,nparams);     //This is the Jacobian matrix of the model equation
                              //The derivative of the cost function w.r.t. parameter j 
                              //will then be: Grad_j=Sum(2*(F(x_i)-Y_i)*J(i,j)), with Sum across data points i
  ColumnVector diff(npts);    //Residuals
  ColumnVector p_plus_h, S_trial,S;
  Matrix Sig;
  ColumnVector fs(nfib), bs(nfib);

  //Compute the Jacobian first using finite differences for each element. Derivatives are analytic for S0 and the volume fractions
  double step;
  ColumnVector typical_scale(nparams); 
  typical_scale=1;
  for (int k=1; k<=nfib; k++){
    int kk = 3+nparams_per_fibre*(k-1); 
    bs(k)=p(kk);
    fs(k) = beta2f(p(kk))*partial_fsum(fs,k-1);
    typical_scale(kk+3)=100;
  }

  //Compute the derivatives with respect to betas, i.e the transformed volume fraction variables
  Matrix f_deriv;
  f_deriv=fractions_deriv(nfib, fs, bs);  

  Sig=forwardModel_compartments(p);
  S=pred_from_compartments(p, Sig);
  diff=S-Y;
  for (int i=1; i<=npts; i++)
    J(i,1)=S(i)/p(1);  //derivatives with respect to S0 are analytic: S_i/S0
    
  for (int n=2; n<=nparams; n++) {
    p_plus_h = p;
    step = SQRTtiny*nonzerosign(p(n))*max(fabs(p(n))*typical_scale(n),1.0);
    step = nonzerosign(tiny)*min(max(fabs(step),1.0e-8),0.1);  //check that 1e-8<step<0.1
    p_plus_h(n)=p(n)+step;
      
    if (n<3 || (n==nparams && m_include_f0)){ //for d all compartments will change. Also for f0, if included. Use for those params numerical differentiation
      S_trial=forwardModel(p_plus_h);
      for (int i=1; i<=npts; i++)
	J(i,n)=(S_trial(i)-S(i))/step; 
    }
    else{  //for the other params, update only the signal of the relevant fibre compartment
      int fib_indx=(int)ceil((float)(n-2)/(float)nparams_per_fibre); //index indicating in which fibre compartment parameter n belongs to
      if (n==3+nparams_per_fibre*(fib_indx-1)){ //Then we have a volume fraction, derivative is analytic.
	for (int i=1; i<=npts; i++){
	  J(i,n)=0;
	  for (int j=1; j<=nfib; j++){
	    if (f_deriv(j,fib_indx)!=0)
	      J(i,n) += p(1)*(Sig(i,j+1)-Sig(i,1))*f_deriv(j,fib_indx);
	  }
	}
      }
      else{       //for all other params, use numerical differentiation
	S_trial=pred_from_compartments(p_plus_h, Sig,fib_indx);
	for (int i=1; i<=npts; i++)
	  J(i,n)=(S_trial(i)-S(i))/step;
      }
    }
  } 
  gradv =2*J.t()*diff;
  gradv.Release();
  return(gradv);    
}



//this uses Gauss-Newton approximation
boost::shared_ptr<BFMatrix> PVM_Ball_Watsons::hess(const NEWMAT::ColumnVector& p,boost::shared_ptr<BFMatrix> iptr)const{
  boost::shared_ptr<BFMatrix>   hessm;
  if (iptr && iptr->Nrows()==(unsigned int)p.Nrows() && iptr->Ncols()==(unsigned int)p.Nrows()) hessm = iptr;
  else hessm = boost::shared_ptr<BFMatrix>(new FullBFMatrix(p.Nrows(),p.Nrows()));
  
  Matrix J(npts,nparams);     //This is the Jacobian matrix of the model equation
  ColumnVector p_plus_h, S_trial,S, fs(nfib), bs(nfib);
  Matrix Sig;
  
  //Compute the Jacobian first using finite differences for each element. Derivatives are analytic for S0 and the volume fractions
  double step,sigt;
  ColumnVector typical_scale(nparams); 
  typical_scale=1;
  for (int k=1; k<=nfib; k++){
    int kk = 3+nparams_per_fibre*(k-1); 
    bs(k)=p(kk);
    fs(k) = beta2f(p(kk))*partial_fsum(fs,k-1);
    typical_scale(kk+3)=100;
  }

  //Compute the derivatives with respect to betas, i.e the transformed volume fraction variables
  Matrix f_deriv;
  f_deriv=fractions_deriv(nfib, fs, bs);  

  Sig=forwardModel_compartments(p);
  S=pred_from_compartments(p, Sig);

  for (int i=1; i<=npts; i++)
    J(i,1)=S(i)/p(1);  //derivatives with respect to S0 are analytic: S_i/S0
    
  for (int n=2; n<=nparams; n++) {
    p_plus_h = p;
    step = SQRTtiny*nonzerosign(p(n))*max(fabs(p(n))*typical_scale(n),1.0);
    step = nonzerosign(tiny)*min(max(fabs(step),1.0e-8),0.1);  //check that 1e-8<step<0.1
    p_plus_h(n)=p(n)+step;
      
    if (n<3 || (n==nparams && m_include_f0)){ //for d all compartments will change. Also for f0, if included. Use for those params numerical differentiation
      S_trial=forwardModel(p_plus_h);
      for (int i=1; i<=npts; i++)
	J(i,n)=(S_trial(i)-S(i))/step; 
    }
    else{  //for the other params, update only the signal of the relevant fibre compartment
      int fib_indx=(int)ceil((float)(n-2)/(float)nparams_per_fibre); //index indicating in which fibre compartment parameter n belongs to
      if (n==3+nparams_per_fibre*(fib_indx-1)){ //Then we have a volume fraction, derivative is analytic.
	for (int i=1; i<=npts; i++){
	  J(i,n)=0;
	  for (int j=1; j<=nfib; j++){
	    if (f_deriv(j,fib_indx)!=0)
	      J(i,n) += p(1)*(Sig(i,j+1)-Sig(i,1))*f_deriv(j,fib_indx); 
	  }
	}
      }
      else{       //for all other params, use numerical differentiation
	S_trial=pred_from_compartments(p_plus_h, Sig,fib_indx);
	for (int i=1; i<=npts; i++)
	  J(i,n)=(S_trial(i)-S(i))/step;
      }
    }
  } 
  
  for (int i=1; i<=p.Nrows(); i++){
    for (int j=i; j<=p.Nrows(); j++){
      sigt = 0.0;
      for(int k=1;k<=J.Nrows();k++)
	sigt += 2*(J(k,i)*J(k,j));
      hessm->Set(i,j,sigt);
    }
  }
  for (int j=1; j<=p.Nrows(); j++) {
    for (int i=j+1; i<=p.Nrows(); i++) {
      hessm->Set(i,j,hessm->Peek(j,i));
    }
  }
  return(hessm);
}



float PVM_Ball_Watsons::isoterm(const int& pt,const float& _d)const{
  return(std::exp(-bvals(1,pt)*_d));
}


NEWMAT::ReturnMatrix PVM_Ball_Watsons::fractions_deriv(const int& nfib, const ColumnVector& fs, const ColumnVector& bs) const{
  NEWMAT::Matrix Deriv(nfib,nfib);
  float fsum;
  Deriv=0;
  for (int j=1; j<=nfib; j++)
    for (int k=1; k<=nfib; k++){
      if (j==k){
	fsum=1; 
	for (int n=1; n<=j-1; n++)
	  fsum-=fs(n);
	Deriv(j,k)=sin(2*bs(k))*fsum;
      }
      else if (j>k){
	fsum=0;
	for (int n=1; n<=j-1; n++)
	  fsum+=Deriv(n,k);
	Deriv(j,k)=-pow(sin(bs(j)),2.0)*fsum;  
      }
    } 	   
  Deriv.Release();
  return Deriv;
}






///////////////////////////////////////////////////////////////////////////////////////////////
//               USEFUL FUNCTIONS TO CALCULATE DERIVATIVES
///////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////
///////Model 1 (Constrained)
/////////////////////////////////////
// functions
float PVM_single_c::isoterm(const int& pt,const float& _d)const{
  return(std::exp(-bvals(1,pt)*_d));
}
float PVM_single_c::anisoterm(const int& pt,const float& _d,const ColumnVector& x)const{
  float dp = bvecs(1,pt)*x(1)+bvecs(2,pt)*x(2)+bvecs(3,pt)*x(3);
  return(std::exp(-bvals(1,pt)*_d*dp*dp));
}
// 1st order derivatives
float PVM_single_c::isoterm_lambda(const int& pt,const float& lambda)const{
  return(-2*bvals(1,pt)*lambda*std::exp(-bvals(1,pt)*lambda*lambda));
}
float PVM_single_c::anisoterm_lambda(const int& pt,const float& lambda,const ColumnVector& x)const{
  float dp = bvecs(1,pt)*x(1)+bvecs(2,pt)*x(2)+bvecs(3,pt)*x(3);
  return(-2*bvals(1,pt)*lambda*dp*dp*std::exp(-bvals(1,pt)*lambda*lambda*dp*dp));
}
float PVM_single_c::anisoterm_th(const int& pt,const float& _d,const ColumnVector& x,const float& _th,const float& _ph)const{
  float dp  = bvecs(1,pt)*x(1)+bvecs(2,pt)*x(2)+bvecs(3,pt)*x(3);
  float dp1 = cos(_th)*(bvecs(1,pt)*cos(_ph) + bvecs(2,pt)*sin(_ph)) - bvecs(3,pt)*sin(_th);
  return(-2*bvals(1,pt)*_d*dp*dp1*std::exp(-bvals(1,pt)*_d*dp*dp));
}
float PVM_single_c::anisoterm_ph(const int& pt,const float& _d,const ColumnVector& x,const float& _th,const float& _ph)const{
  float dp  = bvecs(1,pt)*x(1)+bvecs(2,pt)*x(2)+bvecs(3,pt)*x(3);
  float dp1 = sin(_th)*(-bvecs(1,pt)*sin(_ph) + bvecs(2,pt)*cos(_ph));
  return(-2*bvals(1,pt)*_d*dp*dp1*std::exp(-bvals(1,pt)*_d*dp*dp));
}

NEWMAT::ReturnMatrix PVM_single_c::fractions_deriv(const int& nfib, const ColumnVector& fs, const ColumnVector& bs) const{
  NEWMAT::Matrix Deriv(nfib,nfib);
  float fsum;
  Deriv=0;
  for (int j=1; j<=nfib; j++)
    for (int k=1; k<=nfib; k++){
      if (j==k){
	fsum=1; 
	for (int n=1; n<=j-1; n++)
	  fsum-=fs(n);
	Deriv(j,k)=sin(2*bs(k))*fsum;
      }
      else if (j>k){
	fsum=0;
	for (int n=1; n<=j-1; n++)
	  fsum+=Deriv(n,k);
	Deriv(j,k)=-pow(sin(bs(j)),2.0)*fsum;  
      }
    } 	   
  Deriv.Release();
  return Deriv;
}


/////////////////////////////////////
////////Model 1 (Old)
/////////////////////////////////////
//functions
float PVM_single::isoterm(const int& pt,const float& _d)const{
  return(std::exp(-bvals(1,pt)*_d));
}
float PVM_single::anisoterm(const int& pt,const float& _d,const ColumnVector& x)const{
  float dp = bvecs(1,pt)*x(1)+bvecs(2,pt)*x(2)+bvecs(3,pt)*x(3);
  return(std::exp(-bvals(1,pt)*_d*dp*dp));
}
float PVM_single::bvecs_fibre_dp(const int& pt,const float& _th,const float& _ph)const{
  float angtmp = cos(_ph-beta(pt))*sinalpha(pt)*sin(_th) + cosalpha(pt)*cos(_th);
  return(angtmp*angtmp);
}
float PVM_single::bvecs_fibre_dp(const int& pt,const ColumnVector& x)const{
  float dp = bvecs(1,pt)*x(1)+bvecs(2,pt)*x(2)+bvecs(3,pt)*x(3);
  return(dp*dp);
}
// 1st order derivatives
float PVM_single::isoterm_d(const int& pt,const float& _d)const{
  return(-bvals(1,pt)*std::exp(-bvals(1,pt)*_d));
}
float PVM_single::anisoterm_d(const int& pt,const float& _d,const ColumnVector& x)const{
  float dp = bvecs(1,pt)*x(1)+bvecs(2,pt)*x(2)+bvecs(3,pt)*x(3);
  return(-bvals(1,pt)*dp*dp*std::exp(-bvals(1,pt)*_d*dp*dp));
}
float PVM_single::anisoterm_th(const int& pt,const float& _d,const ColumnVector& x,const float& _th,const float& _ph)const{
  float dp  = bvecs(1,pt)*x(1)+bvecs(2,pt)*x(2)+bvecs(3,pt)*x(3);
  float dp1 = cos(_th)*(bvecs(1,pt)*cos(_ph) + bvecs(2,pt)*sin(_ph)) - bvecs(3,pt)*sin(_th);
  return(-2*bvals(1,pt)*_d*dp*dp1*std::exp(-bvals(1,pt)*_d*dp*dp));
}
float PVM_single::anisoterm_ph(const int& pt,const float& _d,const ColumnVector& x,const float& _th,const float& _ph)const{
  float dp  = bvecs(1,pt)*x(1)+bvecs(2,pt)*x(2)+bvecs(3,pt)*x(3);
  float dp1 = sin(_th)*(-bvecs(1,pt)*sin(_ph) + bvecs(2,pt)*cos(_ph));
  return(-2*bvals(1,pt)*_d*dp*dp1*std::exp(-bvals(1,pt)*_d*dp*dp));
}
// 2nd order derivatives
float PVM_single::isoterm_dd(const int& pt,const float& _d)const{
  return(bvals(1,pt)*bvals(1,pt)*std::exp(-bvals(1,pt)*_d));
}
float PVM_single::anisoterm_dd(const int& pt,const float& _d,const ColumnVector& x)const{
  float dp = bvecs(1,pt)*x(1)+bvecs(2,pt)*x(2)+bvecs(3,pt)*x(3);
  dp *= dp;
  return(bvals(1,pt)*dp*bvals(1,pt)*dp*std::exp(-bvals(1,pt)*_d*dp));
}
float PVM_single::anisoterm_dth(const int& pt,const float& _d,const ColumnVector& x,const float& _th,const float& _ph)const{
  float dp  = bvecs(1,pt)*x(1)+bvecs(2,pt)*x(2)+bvecs(3,pt)*x(3);
  float dp1 = cos(_th)*(bvecs(1,pt)*cos(_ph) + bvecs(2,pt)*sin(_ph)) - bvecs(3,pt)*sin(_th);
  return( -2*bvals(1,pt)*dp*dp1*(1-bvals(1,pt)*_d*dp*dp)*std::exp(-bvals(1,pt)*_d*dp*dp) );
}
float PVM_single::anisoterm_dph(const int& pt,const float& _d,const ColumnVector& x,const float& _th,const float& _ph)const{
  float dp  = bvecs(1,pt)*x(1)+bvecs(2,pt)*x(2)+bvecs(3,pt)*x(3);
  float dp1 = sin(_th)*(-bvecs(1,pt)*sin(_ph) + bvecs(2,pt)*cos(_ph));
  return( -2*bvals(1,pt)*dp*dp1*(1-bvals(1,pt)*_d*dp*dp)*std::exp(-bvals(1,pt)*_d*dp*dp) );
}
float PVM_single::anisoterm_thth(const int& pt,const float& _d,const ColumnVector& x,const float& _th,const float& _ph)const{
  float dp  = bvecs(1,pt)*x(1)+bvecs(2,pt)*x(2)+bvecs(3,pt)*x(3);
  return( -2*bvals(1,pt)*_d*std::exp(-bvals(1,pt)*_d*dp*dp)* ( (1-2*bvals(1,pt)*dp*dp) -dp*dp ) );
}
float PVM_single::anisoterm_phph(const int& pt,const float& _d,const ColumnVector& x,const float& _th,const float& _ph)const{
  float dp  = bvecs(1,pt)*x(1)+bvecs(2,pt)*x(2)+bvecs(3,pt)*x(3);
  float dp1 = (1-cos(2*_th))/2.0;
  float dp2 = -bvecs(1,pt)*x(1) - bvecs(2,pt)*x(2);
  return( -2*bvals(1,pt)*_d*std::exp(-bvals(1,pt)*_d*dp*dp)* ( (1-2*bvals(1,pt)*dp*dp)*dp1 +dp*dp2 ) );
}
float PVM_single::anisoterm_thph(const int& pt,const float& _d,const ColumnVector& x,const float& _th,const float& _ph)const{
  float dp  = bvecs(1,pt)*x(1)+bvecs(2,pt)*x(2)+bvecs(3,pt)*x(3);
  float dp2 = cos(_th)*(-bvecs(1,pt)*sin(_ph) + bvecs(2,pt)*cos(_ph));
  return( -2*bvals(1,pt)*_d*std::exp(-bvals(1,pt)*_d*dp*dp)* ( dp*dp2 ) );
}



////// NOW FOR MULTISHELL
// functions
float PVM_multi::isoterm(const int& pt,const float& _a,const float& _b)const{
  return(std::exp(-_a*std::log(1+bvals(1,pt)*_b)));
}
float PVM_multi::anisoterm(const int& pt,const float& _a,const float& _b,const ColumnVector& x,const int Gamma_for_ball_only)const{
  float dp = bvecs(1,pt)*x(1)+bvecs(2,pt)*x(2)+bvecs(3,pt)*x(3);
  switch (Gamma_for_ball_only){
  case 1:
    return(std::exp(-bvals(1,pt)*_a*_b*dp*dp));
  case 2:
    return(std::exp(-bvals(1,pt)*3*_a*_b*m_invR*((1-m_R)*dp*dp+m_R)));
  default:
    return(std::exp(-_a*std::log(1+bvals(1,pt)*_b*(dp*dp))));
  }
}

// 1st order derivatives
float PVM_multi::isoterm_a(const int& pt,const float& _a,const float& _b)const{
    return(-std::log(1+bvals(1,pt)*_b)*std::exp(-_a*std::log(1+bvals(1,pt)*_b)));
}
float PVM_multi::isoterm_b(const int& pt,const float& _a,const float& _b)const{
      return(-_a*bvals(1,pt)/(1+bvals(1,pt)*_b)*std::exp(-_a*std::log(1+bvals(1,pt)*_b)));
}

float PVM_multi::anisoterm_a(const int& pt,const float& _a,const float& _b,const ColumnVector& x,const int Gamma_for_ball_only)const{
  float dp = bvecs(1,pt)*x(1)+bvecs(2,pt)*x(2)+bvecs(3,pt)*x(3); float dp2;
  switch (Gamma_for_ball_only){
  case 1:
    return(-bvals(1,pt)*_b*dp*dp*std::exp(-bvals(1,pt)*_a*_b*dp*dp));
  case 2:
    dp2=bvals(1,pt)*3*_b*m_invR*((1-m_R)*dp*dp+m_R);
    return(-dp2*std::exp(-dp2*_a));
  default:
    return(-std::log(1+bvals(1,pt)*(dp*dp)*_b)*std::exp(-_a*std::log(1+bvals(1,pt)*(dp*dp)*_b)));
  }
}

float PVM_multi::anisoterm_b(const int& pt,const float& _a,const float& _b,const ColumnVector& x,const int Gamma_for_ball_only)const{
  float dp = bvecs(1,pt)*x(1)+bvecs(2,pt)*x(2)+bvecs(3,pt)*x(3); float dp2;
  switch (Gamma_for_ball_only){
  case 1:
    return(-bvals(1,pt)*_a*dp*dp*std::exp(-bvals(1,pt)*_a*_b*dp*dp));
  case 2:
    dp2=bvals(1,pt)*3*_a*m_invR*((1-m_R)*dp*dp+m_R);
    return(-dp2*std::exp(-dp2*_b));
  default:
    return(-_a*bvals(1,pt)*(dp*dp)/(1+bvals(1,pt)*(dp*dp)*_b)*std::exp(-_a*std::log(1+bvals(1,pt)*(dp*dp)*_b)));
  }
}
float PVM_multi::anisoterm_th(const int& pt,const float& _a,const float& _b,const ColumnVector& x,const float& _th,const float& _ph,const int Gamma_for_ball_only)const{
  float dp = bvecs(1,pt)*x(1)+bvecs(2,pt)*x(2)+bvecs(3,pt)*x(3);
  float dp1 = cos(_th)*(bvecs(1,pt)*cos(_ph) + bvecs(2,pt)*sin(_ph)) - bvecs(3,pt)*sin(_th); float dp2;
  switch (Gamma_for_ball_only){
  case 1:
    return(-2*bvals(1,pt)*_a*_b*dp*dp1*std::exp(-bvals(1,pt)*_a*_b*dp*dp));
  case 2:
    dp2=2*bvals(1,pt)*3*_a*_b*m_invR*(1-m_R)*dp1;
    return(-dp2*std::exp(-bvals(1,pt)*3*_a*_b*m_invR*((1-m_R)*dp*dp+m_R)));
  default:
    return(-_a*_b*bvals(1,pt)/(1+bvals(1,pt)*(dp*dp)*_b)*std::exp(-_a*std::log(1+bvals(1,pt)*(dp*dp)*_b))*2*dp*dp1);
  }
}
float PVM_multi::anisoterm_ph(const int& pt,const float& _a,const float& _b,const ColumnVector& x,const float& _th,const float& _ph,const int Gamma_for_ball_only)const{
  float dp = bvecs(1,pt)*x(1)+bvecs(2,pt)*x(2)+bvecs(3,pt)*x(3);
  float dp1 = sin(_th)*(-bvecs(1,pt)*sin(_ph) + bvecs(2,pt)*cos(_ph)); float dp2;
  switch (Gamma_for_ball_only){
  case 1:
    return(-2*bvals(1,pt)*_a*_b*dp*dp1*std::exp(-bvals(1,pt)*_a*_b*dp*dp));
  case 2:
    dp2=2*bvals(1,pt)*3*_a*_b*m_invR*(1-m_R)*dp1;
    return(-dp2*std::exp(-bvals(1,pt)*3*_a*_b*m_invR*((1-m_R)*dp*dp+m_R)));
  default:
    return(-_a*_b*bvals(1,pt)/(1+bvals(1,pt)*(dp*dp)*_b)*std::exp(-_a*std::log(1+bvals(1,pt)*(dp*dp)*_b))*2*dp*dp1);
  }
}

