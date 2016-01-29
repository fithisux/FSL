
/*  POSSUM
    Ivana Drobnjak & Mark Jenkinson
    Copyright (C) 2005-2007 University of Oxford  */

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

//POSSUM-FUNCTIONS

#include <assert.h>

#include <iostream>//standard c++ library
#include <string>  //include string class
#include <fstream> //to read and write to the file
#include <unistd.h>//what is this?
#include <vector> //for precalcuation of some variables

#include "newmatap.h"//including NEWMAT library
#include "newmatio.h"//
#include "newimage/newimageall.h"//including NEW IMAGE library
#include "possumfns.h"//including possum functions I made 
#include "libprob.h"
#include "miscmaths/miscprob.h"
#include "newimage/costfns.h" 
#include <vector>

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

using namespace std;
using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace NEWMAT;

//string title="possumfns (Version 2.0)\nCopyright(c) 2003, University of Oxford (Ivana Drobnjak)";

vector<RowVector> gstatic;
double* g1static;
double* g2static;
double* g3static;

vector<double> g1motion;
vector<double> g2motion;
vector<double> g3motion;
vector<double> g4motion;

vector<double> rotmotion1;
vector<double> rotmotion2;
vector<double> rotmotion3;
vector<double> transmotion;

vector<double> b0tmp111(1,0.0);
vector<double> b0tmp112(1,0.0);
vector<double> b0tmp113(1,0.0);
vector<double> b0tmp221(1,0.0);
vector<double> b0tmp222(1,0.0);
vector<double> b0tmp223(1,0.0);
vector<double> b0tmp331(1,0.0);
vector<double> b0tmp332(1,0.0);
vector<double> b0tmp333(1,0.0);

vector<double> b0tmp111freq(1,0.0);
vector<double> b0tmp112freq(1,0.0);
vector<double> b0tmp113freq(1,0.0);
vector<double> b0tmp221freq(1,0.0);
vector<double> b0tmp222freq(1,0.0);
vector<double> b0tmp223freq(1,0.0);
vector<double> b0tmp331freq(1,0.0);
vector<double> b0tmp332freq(1,0.0);
vector<double> b0tmp333freq(1,0.0);

double glo_zero=1e-20;
//TABLE SINC
double* table_sinc=0;//table for SINC  in range [0,Dsinc]
int glo_Nsinc=5000000;//no of elements
double glo_Dsinc=500.00;//domain
double glo_dsinc=glo_Dsinc/glo_Nsinc;//step size for SINC
double glo_idsinc=1/glo_dsinc;//inverse of step size for fster calc

//TABLE SIN AND COS
double* table_sin=0;//table for SIN in range [0,2pi]
double* table_cos=0;//table for COS in range [0,2pi]
int glo_Nsin=6000000;//no of elements
double glo_Dsin=2.0*M_PI;//domain
double glo_dsin=glo_Dsin/glo_Nsin;//step size for SIN and COS
double glo_idsin=1/glo_dsin;//inverse of the step size
double glo_twopi=2*M_PI;
double glo_itwopi=1/glo_twopi;

double glo_cx;
double glo_cy;
double glo_cz;
    
//SOME CONSTANTS
const double gammabar=42.58*1e06;//(in Hz/T)
const double gama=2*M_PI*gammabar;


/////////
//GENERAL
///////////////////////////////
double round_ivana(const double x, const int n){
  //rounds a number up to n digits of precision
  double xx,xxx,xxxx,nn;
  nn=MISCMATHS::pow(10.0f,(double) n);
  xx=nn*x;
  if (x>0) xxx=floor((float)xx+0.5);
  else xxx=ceil(xx-0.5);
  xxxx=xxx/nn;
  return xxxx;
}

/////////////////////////////
double norm(const RowVector q){
  int x=q.Ncols();
  double sqsum=0.0;
  for (int i=1;i<=x;i++){
    sqsum+=q(i)*q(i);
  }
  double norm_q=sqrt(sqsum);
  return norm_q;
}
////////////////////////////////////////
void coeff(const double xold, const double xnew, const double told, const double tnew, double& a, double& b){//needs improvement
  // return slope and intercept of a line
  //cout<<"xold="<<xold<<" xnew="<<xnew<<" told="<<told<<" tnew="<<tnew<<" a="<<a<<" b="<<b<<endl;
  if ((tnew-told)>1e-12){
    a=xold-told*(xnew-xold)/(tnew-told);	//intercept
    b=(xnew-xold)/(tnew-told);			//slope
  }
  //cout<<"xold="<<xold<<" xnew="<<xnew<<" told="<<told<<" tnew="<<tnew<<" a="<<a<<" b="<<b<<endl;
}
//////////////////////////////////////////from mj code costfns.cc in NEWIMAGE
double mj_sinc(const double x){
  if (fabs(x)<1e-7) { return 1.0-fabs(x); }
  double y=M_PI*x;
  return sin(y)/y;
}
/////////////////////////////////////////////////
Matrix rot(const double alpha,const string axes){//WE USE RIGHT HAND RULE FOR THE ROTATIONS AND ORDER IN ROTATION SEQUENCE Rz*Ry*Rx
  Matrix R(3,3);
  if (axes=="z"){
    //rotation right hand rule
    R <<cos(alpha)<<-sin(alpha)<<0
      <<sin(alpha)<<cos(alpha)<<0
      <<0<<0<<1;
  }
  if (axes=="x"){
    //rotation right hand rule
   R  <<1<<0<<0
      <<0<<cos(alpha)<<-sin(alpha)
      <<0<<sin(alpha)<<cos(alpha);
  }
  if (axes=="y"){
    //rotation right hand rule
    R  <<cos(alpha)<<0<<-sin(alpha)
       <<0<<1<<0
       <<sin(alpha)<<0<<cos(alpha);
  }
 return R;
 }
//////////////////////////////////////////////////////////////
//CONVERSIONS BETWEEN EULER ANGLES, QUATERNIONS and MATRICES
//////////////////////////////////////////////////////////////////////////
RowVector mult_quaternions(const RowVector a,const RowVector b){//the order of the multiplication is important as it is not comutative!!
  //quaternion is of the form c=(cos(angle/2),axisx*sin(angle/2),axisy*sin(angle/2),axisz*sin(angle/2))=(w,x,y,z)
  RowVector c(4);
  c(1)=a(1)*b(1)-a(2)*b(2)-a(3)*b(3)-a(4)*b(4);
  c(2)=a(1)*b(2)+a(2)*b(1)+a(3)*b(4)-a(4)*b(3);
  c(3)=a(1)*b(3)-a(2)*b(4)+a(3)*b(1)+a(4)*b(2);
  c(4)=a(1)*b(4)+a(2)*b(3)-a(3)*b(2)+a(4)*b(1);
  return c;
}
////////////////////////////////////////////////////////////////////////////////
RowVector euler_to_quaternion(const double a,const double b,const double c){//ok
  //constructed by having one basic quaternion for each rotation Rz(c), Ry(b) and Rx(a) and then multiplying them in order qz*qy*qx 
  //RowVector qx(4),qy(4),qz(4),q(4);
  //qx<<cos(a/2)<<sin(a/2)<<0.0<<0.0; //the norm is one
  //qy<<cos(b/2)<<0.0<<sin(b/2)<<0.0; //the norm is one
  //qz<<cos(c/2)<<0.0<<0.0<<sin(c/2); //the norm is one 
  //RowVector tmp(4);
  //tmp=mult_quaternions(qy,qx);
  //q=mult_quaternions(qz,tmp);// the norm of qz*qy*qx has to be also one as the norm(q1*q2)=norm(q1)*norm(q2)
  double cz=cos(c/2);
  double cy=cos(b/2);
  double cx=cos(a/2);
  double sz=sin(c/2);
  double sy=sin(b/2);
  double sx=sin(a/2);
  RowVector q(4);
  q<<cz*cy*cx+sz*sy*sx<<cz*cy*sx-sz*sy*cx<<cz*sy*cx+sz*cy*sx<<sz*cy*cx-cz*sy*sx;
  double norm_q=norm(q);//has to be one from the way we constructed q
  if (fabs(norm_q-1)>1e-12) {
      cout<<"Warning, euler_to_quat is producing quaternions with the norm different from 1!"<<"Norm difference from one is "<<norm_q-1<<"Quat is "<<q<<endl; 
  }
  return q;
}
/////////////////////////////////////////////////////
RowVector quaternion_to_angleaxis(const RowVector q){//ok, but be carefull when applying it after angleaxis_to_quaternion to same vectors
   //angleaxis is of the form (angle,axisx,axisy,axisz)
  double norm_q=norm(q);
  RowVector q_n(4);
  q_n=q;
  //cout<<norm_q-1<<endl;
  //cout<<"quaternion is "<<q_n<<endl;
  if (fabs(norm_q-1)>1e-12) {
    cout<<"Warning, quaternion is not normalised !"<<"Norm_q-1 is "<<norm_q-1<<"q is "<<q<<endl;
    q_n(1)=q(1)/norm_q;
    q_n(2)=q(2)/norm_q;
    q_n(3)=q(3)/norm_q;
    q_n(4)=q(4)/norm_q;
  }
  double angle;
  RowVector axis(3);
  if ((q_n(1)-1)>0) q_n(1)=1; //the difference can be really small i.e. 1e-16, but acos gives nan for it (for 1+1e-16) 
  if ((q_n(1)+1)<0) q_n(1)=-1;
  angle=acos(q_n(1))*2;
  //cout<<q_n(1)-1<<"  angle  "<<angle<<endl;
  //cout<<"angle  "<<angle<<endl;//0<=acos()<=pi! Although we always get positive angle out, it is the axis that determines the rotation, so if angle<0 then negative axis
  //double s=sqrt(1-q_n(1)*q_n(1));is not good as it can create 0 0 0 axis
  double s= sqrt(q_n(2)*q_n(2)+q_n(3)*q_n(3)+q_n(4)*q_n(4));
  //s=fabs(sin(angle/2));
  //cout<<"s  "<<s<<endl;
  if (fabs(s)<1e-12){
    angle=0;
    axis(1)=0; 
    axis(2)=0; 
    axis(3)=1; //taken so that the axis is always normalised but has no real meaning as the angle iz zero
  } else {
      axis(1)=q_n(2)/s;
      axis(2)=q_n(3)/s;
      axis(3)=q_n(4)/s;    
  }
  double norm_axis=norm(axis);//has to be one if the quaternion is normalised, easily proved that normalised quaternion iff normalised axis
  if (fabs(norm_axis)<1e-12){
    cout<<"Warning in quaternion_to_angleaxis: Norm of the axis is ZERO!"<<endl;
  }
  RowVector angleaxis(4);
  angleaxis(1)=angle;
  angleaxis(2)=axis(1);
  angleaxis(3)=axis(2);
  angleaxis(4)=axis(3);
  return angleaxis; 
} 
////////////////////////////////////////////////////////////////////////
RowVector angleaxis_to_quaternion(const RowVector angleaxis){
  //angleaxis is of the form (angle,x,y,z)
  double angle=angleaxis(1);
  RowVector axis(3);
  axis<<angleaxis(2)<<angleaxis(3)<<angleaxis(4);
  double norm_axis=norm(axis);
  if (fabs(norm_axis-1)>1e-12) {
    axis(1)=axis(1)/norm_axis;
    axis(2)=axis(2)/norm_axis;
    axis(3)=axis(3)/norm_axis;
  }
  double x=axis(1)*sin(angle/2);
  double y=axis(2)*sin(angle/2);
  double z=axis(3)*sin(angle/2);
  double w=cos(angle/2);
  RowVector q(4);
  q<<w<<x<<y<<z;
  return q;
}
////////////////////////////////////////////////////////
Matrix quaternion_to_matrix(const RowVector q){//not using
 double norm_q=norm(q);
 RowVector q_n(4);
 q_n=q;
  if (fabs(norm_q-1)>1e-12) {
    q_n(1)=q(1)/norm_q;
    q_n(2)=q(2)/norm_q;
    q_n(3)=q(3)/norm_q;
    q_n(4)=q(4)/norm_q;
  }
  double x=q_n(2);
  double y=q_n(3);
  double z=q_n(4);
  double w=q_n(1);
  Matrix R(3,3);
  R<<1-2*y*y-2*z*z<<2*x*y-2*z*w<<2*x*z+2*y*w
   <<2*x*y+2*z*w<<1-2*x*x-2*z*z<<2*y*z-2*x*w
   <<2*x*z-2*y*w<<2*y*z+2*x*w<<1-2*x*x-2*y*y;
  return R;
}
////////////////////////////////////////////////////////////////////////////////
RowVector matrix_to_quaternion(const Matrix M){//not using
  double trace = M(1,1) + M(2,2) + M(3,3) + 1;
  double s,x,y,z,w;
  if( trace > 0 ) {
    s = 0.5 / sqrt(trace);
    w = 0.25f/s;
    x = (M(3,2) - M(2,3))*s;
    y = (M(1,3) - M(3,1))*s;
    z = (M(2,1) - M(1,2))*s;
  }  
  else {
    if ( M(1,1) > M(2,2) && M(1,1) > M(3,3) ) {
      double s = 2 * sqrt( 1 + M(1,1) - M(2,2) - M(3,3));
      x = 0.25 * s;
      y = (M(1,2)+M(2,1))/s;
      z = (M(1,3)+M(3,1))/s;
      w = (M(2,3)-M(3,2))/s;
    } 
    else if (M(2,2) > M(3,3)) {
      double s = 2 * sqrt( 1 + M(2,2) - M(1,1) - M(3,3));
      x = (M(1,2) + M(2,1) ) / s;
      y = 0.25f * s;
      z = (M(2,3) + M(3,2) ) / s;
      w = (M(1,3) - M(3,1) ) / s;
    } 
    else {
      double s = 2 * sqrt( 1 + M(3,3) - M(1,1) - M(2,2) );
      x = (M(1,3) + M(3,1) ) / s;
      y = (M(2,3) + M(3,2) ) / s;
      z = 0.25 * s;
      w = (M(1,2) - M(2,1) ) / s;
    }
  }
  RowVector q(4);
  q<<w<<x<<y<<z;
  return q;
}
/////////////////////////////////////////////////////////////////////////////////////////////
Matrix euler_to_matrix(const double a,const double b, const double c){//not using
  Matrix R(3,3);
  R(1,1)=cos(b)*cos(c);
  R(1,2)=-cos(b)*sin(c);
  R(1,3)=-sin(b);
  R(2,1)=-sin(a)*sin(b)*cos(c)+cos(a)*sin(c);
  R(2,2)=sin(a)*sin(b)*sin(c)+cos(a)*cos(c);
  R(2,3)=-sin(a)*cos(b);
  R(3,1)=cos(a)*sin(b)*cos(c)+sin(a)*sin(c);
  R(3,2)=-cos(a)*sin(b)*sin(c)+sin(a)*cos(c);
  R(3,3)=cos(a)*cos(b);
  return R;
}
////////////////////////////////
//STUFF FOR THE ROTATIONS
////////////////////////////////////////////////////////////////////////////////////////////
ReturnMatrix axismat(const RowVector angleaxis){//look more into it!
  // returns A where R = I + sin(angle) A + (1 - cos(angle)) A^2
  RowVector axis(3);
  axis<<angleaxis(2)<<angleaxis(3)<<angleaxis(4);
  double norm_axis=norm(axis);
  if (fabs(norm_axis-1)>1e-12) {
    if (fabs(norm_axis)>1e-12){
      axis(1)=axis(1)/norm_axis;
      axis(2)=axis(2)/norm_axis;
      axis(3)=axis(3)/norm_axis;
    }
    else {
      //cout<<"Warning in axismatm: Norm of the axis is ZERO!"<<endl;
    }
  }
  Matrix m(3,3);
  m <<0.0<<-axis(3)<<axis(2)
    <<axis(3)<<0.0<<-axis(1)
    <<-axis(2)<<axis(1)<<0.0;
  m.Release(); 
  return m;
}
///////////////////////////////////////////////////////////////////////////////////////////////
ReturnMatrix rotmat(const RowVector angleaxis){
  // returns R where R = I + sin(angle) A + (1 - cos(angle) A^2
  Matrix m(3,3);
  Matrix d(3,3);
  d=0;
  d(1,1)=1;d(2,2)=1;d(3,3)=1;
  Matrix A(3,3);
  A=axismat(angleaxis);
  double angle=angleaxis(1);
  m=d+sin(angle)*A+(1-cos(angle))*A*A;
  m.Release();
  return m;
}
////////////////////////////////////////////////////////////////////////////////
ColumnVector free(const ColumnVector m,const double time_i,const RowVector tissue,const double phase_i,const double actinttt){
 // m is the magnetization vector just after the last rf pulse. this modul calculates the magnetization vector after the free precession, just before the next rf pulse
 // tissue:T1,T2,rho
 // time_i=time-rftime where time is time from the bigining of the acquisition and rftime is the time of the last rf pulse
 // phase_i=phase-rfphase the same logic as the time 
  Matrix W(3,3);
  ColumnVector Q(3);
  double e2=exp(-time_i/tissue(2)+actinttt);
  double e1=exp(-time_i/tissue(1));
  W <<e2<<0<<0
    <<0<<e2<<0
    <<0<<0<<e1;
  Q=0.0;
  Q(3)=tissue(3)*(1-e1);
  ColumnVector m1(3);
  m1=W*rot(phase_i,"z")*m+Q;
  return m1;
}
////////////////////////////////////////////////////////////////////////////////
ColumnVector free_cout(const ColumnVector m,const double time_i,const RowVector tissue,const double phase_i,const double actinttt){
 // m is the magnetization vector just after the last rf pulse. this modul calculates the magnetization vector after the free precession, just before the next rf pulse
 // tissue:T1,T2,rho
 // time_i=time-rftime where time is time from the bigining of the acquisition and rftime is the time of the last rf pulse
 // phase_i=phase-rfphase the same logic as the time 
  Matrix W(3,3);
  ColumnVector Q(3);
  double e2=exp(-time_i/tissue(2)+actinttt);
  double e1=exp(-time_i/tissue(1));
  W <<e2<<0<<0
    <<0<<e2<<0
    <<0<<0<<e1;
  Q=0.0;
  Q(3)=tissue(3)*(1-e1);
  ColumnVector m1(3);
  m1=W*rot(phase_i,"z")*m+Q;
  cout<<"Input: m="<<m<<", time="<<time_i<<", tissue="<<tissue<<", phase="<<phase_i<<", actint="<<actinttt<<". Output: e1="<<e1<<", e2="<<e2<<", rot output="<<rot(phase_i,"z")<<", W="<<W<<", Q="<<Q<<", m1="<<m1<<endl;
  return m1;
}
///////////////////////////
//INTERGRALS
//////////////////////////////////////////////////////////////////////////////////////////////////////
double i1(const double gold,const double gnew,const double told,const double tnew){
  // integral: i1 = \int G(t) dt
  double i11;
  double g1,g2;
  g1=0.0;g2=0.0;
  coeff(gold,gnew,told,tnew,g1,g2);
  //i11=tnew*(g1+g2*tnew/2)-told*(g1+g2*told/2);

//tejas-changed 26.10.12
  i11=(tnew-told)*g1+(tnew-told)*(tnew+told)*g2/2;
//tejas-end
  return i11;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
double i1new(const double gold,const double gnew,const double told,const double tnew){
  // integral: i1 = \int G(t) dt
  double i11;
  double g1,g2;
  g1=0.0;g2=0.0;
  coeff(gold,gnew,told,tnew,g1,g2);
  cout<<"Gradient is g1+g2*(tnew-told): g1="<<g1<<"; g2="<<g2<<"; tnew-told="<<tnew-told<<endl;
  //i11=tnew*(g1+g2*tnew/2)-told*(g1+g2*told/2);
  i11=(tnew-told)*g1+(tnew-told)*(tnew+told)*g2/2;

  return i11;
}
//////////////////////////////////////////////////////////////////////////////////////////////
double i2(const double gold,const double gnew,const double aold, const double anew, const double told,const double tnew){
  // integral: i2 = \int sin(\alpha(t)) G(t) dt
  double i22;
  double g1,g2;
  g1=0.0;g2=0.0;
  coeff(gold,gnew,told,tnew,g1,g2);
  double a1,a2;
  a1=0.0;a2=0.0;
  coeff(aold,anew,told,tnew,a1,a2);
  //if ((tnew-told)>0.1)cout<<a2<<endl;
  if (fabs(a2)*(tnew-told)<0.01){
    i22=(2*g1*cos(a2*(tnew-told)/4)*sin(a1+a2*(told+tnew)/2)+g2*tnew*sin(a1+a2*(3*tnew+told)/4)+g2*told*sin(a1+a2*(tnew+3*told)/4))*(Sinc(a2*(tnew-told)/(4*M_PI))*(tnew-told)/2);//teylor done for sin(o(1e-04)) error expected to be less than 0.0016%
      //i22=sin(a1)*i1(gold,gnew,told,tnew);
  }
  else {
    i22=g1*2*(sin(a2*(tnew-told)/2)/a2)*sin(a1+a2*(tnew+told)/2)+g2*(sin(a1)*(tnew*sin(a2*tnew)/a2-told*sin(a2*told)/a2-2*(sin(a2*(tnew-told)/2)/a2)*(sin(a2*(tnew+told)/2)/a2))+cos(a1)*(-tnew*cos(a2*tnew)/a2+told*cos(a2*told)/a2+2*(sin(a2*(tnew-told)/2)/a2)*(cos(a2*(tnew+told)/2)/a2)));
    if (fabs(a2)<0.00001) cout<<"WARNING check: a2="<<a2<<"i22="<<i22<<endl;
    //i22=(-(g1+g2*tnew)*cos(a1+a2*tnew)+(g1+g2*told)*cos(a1+a2*told)+(g2*sin(a1+a2*tnew)-g2*sin(a1+a2*told))/a2)/a2;
  }
  return i22;
}
////////////////////////////////////////////////////////////////////////////////////////
double i2new(const double gold,const double gnew,const double aold, const double anew, const double told,const double tnew){
  //USED ONLY DURING TESTING FOR OUTPUTING VALUES WHEN -V VERSION OF THE CODE
  // integral: i2 = \int sin(\alpha(t)) G(t) dt
  double i22;
  double g1,g2;
  g1=0.0;g2=0.0;
  coeff(gold,gnew,told,tnew,g1,g2);
  double a1,a2;
  a1=0.0;a2=0.0;
  coeff(aold,anew,told,tnew,a1,a2);
  cout<<"Angle is a1+a2*(tnew-told): a1="<<a1<<"; a2="<<a2<<"; tnew-told="<<tnew-told<<endl;
   cout<<"aold= "<<aold<<" anew= "<<anew<<endl;
  //if ((tnew-told)>0.1)cout<<a2<<endl;
  if (fabs(a2)*(tnew-told)<0.01){
    cout<<"fabs(a2)*(tnew-told)="<<fabs(a2)*(tnew-told)<<" < 0.01"<<endl;
    i22=(2*g1*cos(a2*(tnew-told)/4)*sin(a1+a2*(told+tnew)/2)+g2*tnew*sin(a1+a2*(3*tnew+told)/4)+g2*told*sin(a1+a2*(tnew+3*told)/4))*(Sinc(a2*(tnew-told)/(4*M_PI))*(tnew-told)/2);//teylor done for sin(o(1e-04)) error expected to be less than 0.0016%
      //i22=sin(a1)*i1(gold,gnew,told,tnew);
  }
  else {
    cout<<"fabs(a2)*(tnew-told)="<<fabs(a2)*(tnew-told)<<" > 0.01"<<endl;
    i22=g1*2*(sin(a2*(tnew-told)/2)/a2)*sin(a1+a2*(tnew+told)/2)+g2*(sin(a1)*(tnew*sin(a2*tnew)/a2-told*sin(a2*told)/a2-2*(sin(a2*(tnew-told)/2)/a2)*(sin(a2*(tnew+told)/2)/a2))+cos(a1)*(-tnew*cos(a2*tnew)/a2+told*cos(a2*told)/a2+2*(sin(a2*(tnew-told)/2)/a2)*(cos(a2*(tnew+told)/2)/a2)));
    if (fabs(a2)<0.00001) cout<<"WARNING check: a2="<<a2<<"i22="<<i22<<endl;
    //i22=(-(g1+g2*tnew)*cos(a1+a2*tnew)+(g1+g2*told)*cos(a1+a2*told)+(g2*sin(a1+a2*tnew)-g2*sin(a1+a2*told))/a2)/a2;
  }
  return i22;
}
///////////////////////////////////////////////////////////////////////////////////////////
double i3(const double gold,const double gnew,const double aold,const double anew,const double told,const double tnew){
  // integral: i3 = \int cos(\alpha(t)) G(t) dt
  double i33;
  double g1,g2;
  g1=0.0;g2=0.0;
  coeff(gold,gnew,told,tnew,g1,g2);
  double a1,a2;
  a1=0.0;a2=0.0;
  coeff(aold,anew,told,tnew,a1,a2);
  if (fabs(a2)*(tnew-told)<0.01){
    i33=(2*g1*cos(a2*(tnew-told)/4)*cos(a1+a2*(told+tnew)/2)+g2*tnew*cos(a1+a2*(3*tnew+told)/4)+g2*told*cos(a1+a2*(tnew+3*told)/4))*(Sinc(a2*(tnew-told)/(4*M_PI))*(tnew-told)/2);//teylor done for sin(o(1e-04))
    //i33=cos(a1)*i1(gold,gnew,told,tnew);
  }
  else {
    i33=g1*2*(sin(a2*(tnew-told)/2)/a2)*cos(a1+a2*(tnew+told)/2)+g2*(sin(a1)*(tnew*cos(a2*tnew)/a2-told*cos(a2*told)/a2-2*(sin(a2*(tnew-told)/2)/a2)*(cos(a2*(tnew+told)/2)/a2))+cos(a1)*(tnew*sin(a2*tnew)/a2-told*sin(a2*told)/a2-2*(sin(a2*(tnew-told)/2)/a2)*(sin(a2*(tnew+told)/2)/a2)));
    if (fabs(a2)<0.00001) cout<<"WARNING check: a2="<<a2<<"i33="<<i33<<endl;
    //i33=((g1+g2*tnew)*sin(a1+a2*tnew)-(g1+g2*told)*sin(a1+a2*told)+(g2*cos(a1+a2*tnew)-g2*cos(a1+a2*told))/a2)/a2;
  }
  return i33;
}
///////////////////////////////////////////////////////////////////////////////////////////
double i3new(const double gold,const double gnew,const double aold,const double anew,const double told,const double tnew){
  // integral: i3 = \int cos(\alpha(t)) G(t) dt
  double i33;
  double g1,g2;
  g1=0.0;g2=0.0;
  coeff(gold,gnew,told,tnew,g1,g2);
  double a1,a2;
  a1=0.0;a2=0.0;
  coeff(aold,anew,told,tnew,a1,a2);
  if (fabs(a2)*(tnew-told)<0.01){
    i33=(2*g1*cos(a2*(tnew-told)/4)*cos(a1+a2*(told+tnew)/2)+g2*tnew*cos(a1+a2*(3*tnew+told)/4)+g2*told*cos(a1+a2*(tnew+3*told)/4))*(Sinc(a2*(tnew-told)/(4*M_PI))*(tnew-told)/2);//teylor done for sin(o(1e-04))
    //i33=cos(a1)*i1(gold,gnew,told,tnew);
  }
  else {
    i33=g1*2*(sin(a2*(tnew-told)/2)/a2)*cos(a1+a2*(tnew+told)/2)+g2*(sin(a1)*(tnew*cos(a2*tnew)/a2-told*cos(a2*told)/a2-2*(sin(a2*(tnew-told)/2)/a2)*(cos(a2*(tnew+told)/2)/a2))+cos(a1)*(tnew*sin(a2*tnew)/a2-told*sin(a2*told)/a2-2*(sin(a2*(tnew-told)/2)/a2)*(sin(a2*(tnew+told)/2)/a2)));
    if (fabs(a2)<0.00001) cout<<"WARNING check: a2="<<a2<<"i33="<<i33<<endl;
    //i33=((g1+g2*tnew)*sin(a1+a2*tnew)-(g1+g2*told)*sin(a1+a2*told)+(g2*cos(a1+a2*tnew)-g2*cos(a1+a2*told))/a2)/a2;
  }
  return i33;
}
/////////////////////////////////////////////////////////////////////////////////////////////
double i4(const double gold,const double gnew,const double trold,const double trnew, const double told,const double tnew){
  // integral: i4 = \int T(t) G(t) dt    : T(t) is translation
  double i44;
  double g1,g2;
  g1=0.0;g2=0.0;
  coeff(gold,gnew,told,tnew,g1,g2);
  double tr1,tr2;
  tr1=0.0;tr2=0.0;
  coeff(trold,trnew,told,tnew,tr1,tr2);
  double e=tnew-told;
  double a=g1*tr1;
  double b=(g1*tr2+g2*tr1)/2;
  double c=g2*tr2/3;
  i44=e*(a+e*(b+e*c))+told*e*(2*b+3*c*(e+told));
  //i44=tnew*(g1*tr1+tnew*((g1*tr2+g2*tr1)/2+tnew*g2*tr2/3))-told*(g1*tr1+told*((g1*tr2+g2*tr1)/2+told*g2*tr2/3));
  return i44;
}
/////////////////////////////////////////////////////////////////////////////////////////////
double i4new(const double gold,const double gnew,const double trold,const double trnew, const double told,const double tnew){
  //Used ONLY when -v option for testing is on to see soome of the variables. 
  //Matrix M=rotmat(r2);
  
// integral: i4 = \int T(t) G(t) dt    : T(t) is translation
  double i44;
  double g1,g2;
  g1=0.0;g2=0.0;
  coeff(gold,gnew,told,tnew,g1,g2);
  double tr1,tr2;
  tr1=0.0;tr2=0.0;
  coeff(trold,trnew,told,tnew,tr1,tr2);
  cout<<"coeff tr:"<<tr1<<" "<<tr2<<endl;
  double e=tnew-told;
  double a=g1*tr1;
  double b=(g1*tr2+g2*tr1)/2;
  double c=g2*tr2/3;
  i44=e*(a+e*(b+e*c))+told*e*(2*b+3*c*(e+told));
  //i44=tnew*(g1*tr1+tnew*((g1*tr2+g2*tr1)/2+tnew*g2*tr2/3))-told*(g1*tr1+told*((g1*tr2+g2*tr1)/2+told*g2*tr2/3));
  return i44;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double I(const int h,const RowVector& r1,const RowVector& r2, const Matrix& g){
  //Matrix M=rotmat(r2);
  Matrix M=rotmat(r2);
  Matrix A=axismat(r1);
  Matrix AM=A*M;
  Matrix A2M=A*AM;
  Matrix v=g.Row(1)*M.Column(h)+g.Row(2)*AM.Column(h)+g.Row(3)*A2M.Column(h);
  return v(1,1);
}
///////////////////////////////
double Inew(const int h,const RowVector& r1,const RowVector& r2, const Matrix& g){
  //Used ONLY when -v option for testing is on to see soome of the variables. 
  //Matrix M=rotmat(r2);
  Matrix M=rotmat(r2);
  cout<<"M matrix "<<M<<endl;
  Matrix A=axismat(r1);
  cout<<"A matrix "<<A<<endl;
  Matrix AM=A*M;
  Matrix A2M=A*AM;
  Matrix v=g.Row(1)*M.Column(h)+g.Row(2)*AM.Column(h)+g.Row(3)*A2M.Column(h);
  return v(1,1);
}
////////////////////////////////
//STUFF FOR THE SORTER
////////////////////////////////////////////////////////////////////////////////////////////////
Matrix intertranslation(const double tx1,const double ty1,const double tz1,
			const double tx2, const double ty2, const double tz2,
			const ColumnVector vectortime){
  // interpolate between two positions (translations)
  // vectortime is a list of times to interpolate at *except* first point
  //  specifies when tx1,ty1,tz1 occurs and last point for tx2,ty2,tz2
  // output is (tx,ty,tz) per row, with first row for *second* input time
  //  but still including last input time => (tx2,ty2,tz2)
  int d=vectortime.Nrows();
  double dtt=0;
  Matrix inttra(d,3);
  double dt=vectortime(d)-vectortime(1);  
  for (int k=1;k<=d-1;k++){
    inttra(k,1)=dtt*(tx2-tx1)/dt+tx1;
    inttra(k,2)=dtt*(ty2-ty1)/dt+ty1;
    inttra(k,3)=dtt*(tz2-tz1)/dt+tz1;
  dtt=dtt+vectortime(k+1)-vectortime(k); 
  }
  inttra(d,1)=tx2;
  inttra(d,2)=ty2;
  inttra(d,3)=tz2;
  return inttra;
}
///////////////////////////
Matrix intertranslation_old(const double tx1,const double ty1,const double tz1,
			const double tx2, const double ty2, const double tz2,
			const ColumnVector vectortime){
  // interpolate between two positions (translations)
  // vectortime is a list of times to interpolate at *except* first point
  //  specifies when tx1,ty1,tz1 occurs and last point for tx2,ty2,tz2
  // output is (tx,ty,tz) per row, with first row for *second* input time
  //  but still including last input time => (tx2,ty2,tz2)
  int d=vectortime.Nrows();
  double dtt=0;
  Matrix inttra(d,3);
  double dt=vectortime(d)-vectortime(1);  
  for (int k=1;k<=d-1;k++){
    dtt=dtt+vectortime(k+1)-vectortime(k); 
    inttra(k,1)=dtt*(tx2-tx1)/dt+tx1;
    inttra(k,2)=dtt*(ty2-ty1)/dt+ty1;
    inttra(k,3)=dtt*(tz2-tz1)/dt+tz1;
  }
  inttra(d,1)=inttra(d-1,1);
  inttra(d,2)=inttra(d-1,2);
  inttra(d,3)=inttra(d-1,3);
  return inttra;
}
////////////////////
Matrix interrotation(const double a1,const double b1,const double c1,
		     const double a2,const double b2,const double c2,
		     const ColumnVector vectortime){
   // as above but for rotations (input in euler angles, output in angle/axis)
  int d=vectortime.Nrows();
  Matrix introt(d,8);
  RowVector aa(8);
  RowVector q1=euler_to_quaternion(a1,b1,c1); 
  RowVector q2=euler_to_quaternion(a2,b2,c2);  
  RowVector q1c(4); //conjugate of q1
  q1c <<q1(1)<<-q1(2)<<-q1(3)<<-q1(4);
  RowVector q=mult_quaternions(q2,q1c); 
  RowVector qa=quaternion_to_angleaxis(q);  
  RowVector q1a=quaternion_to_angleaxis(q1);
  RowVector q2a=quaternion_to_angleaxis(q2); 
  aa<<0.0<<qa(2)<<qa(3)<<qa(4)<<q1a(1)<<q1a(2)<<q1a(3)<<q1a(4);
  introt.Row(1)=aa;
  if (d>2){
    double dtt=0;
    double dt=vectortime(d)-vectortime(1);
    RowVector angle(d-2);
    for (int k=1;k<=d-2;k++){ 
      dtt=dtt+vectortime(k+1)-vectortime(k); 
      angle(k)=dtt*qa(1)/dt;
      aa<<angle(k)<<qa(2)<<qa(3)<<qa(4)<<q1a(1)<<q1a(2)<<q1a(3)<<q1a(4);
      introt.Row(k+1)=aa;
    }
  }
  aa<<0.0<<qa(2)<<qa(3)<<qa(4)<<q2a(1)<<q2a(2)<<q2a(3)<<q2a(4); //first angle just for help but it means nothing, with all zeros axis rot matrix witll be zero
  introt.Row(d)=aa;//begining of new interval
  return introt;
}
/////////////////////////////////////////////////////////////////////////////////////////////
Matrix interrotation_old(const double a1,const double b1,const double c1,
		     const double a2,const double b2,const double c2,
		     const ColumnVector vectortime){
   // as above but for rotations (input in euler angles, output in angle/axis)
  int d=vectortime.Nrows();
  Matrix introt(d,8);
  RowVector aa(8);
  RowVector q1=euler_to_quaternion(a1,b1,c1); 
  RowVector q2=euler_to_quaternion(a2,b2,c2);  
  RowVector q1c(4); //conjugate of q1
  q1c <<q1(1)<<-q1(2)<<-q1(3)<<-q1(4);
  RowVector q=mult_quaternions(q2,q1c); 
  RowVector qa=quaternion_to_angleaxis(q);  
  RowVector q1a=quaternion_to_angleaxis(q1);
  RowVector q2a=quaternion_to_angleaxis(q2); 
  if (d>2){
    double dtt=0;
    double dt=vectortime(d)-vectortime(1);
    RowVector angle(d-2);
    for (int k=1;k<=d-2;k++){ 
      dtt=dtt+vectortime(k+1)-vectortime(k); 
      angle(k)=dtt*qa(1)/dt;
      aa<<angle(k)<<qa(2)<<qa(3)<<qa(4)<<q1a(1)<<q1a(2)<<q1a(3)<<q1a(4);
      introt.Row(k)=aa;
    }
  }
  aa<<qa(1)<<qa(2)<<qa(3)<<qa(4)<<q1a(1)<<q1a(2)<<q1a(3)<<q1a(4);
  introt.Row(d-1)=aa;
  aa<<0.0<<qa(2)<<qa(3)<<qa(4)<<q2a(1)<<q2a(2)<<q2a(3)<<q2a(4); //first angle just for help but it means nothing, with all zeros axis rot matrix witll be zero
  introt.Row(d)=aa;//begining of new interval
  return introt;
}
////////////////////////////////////////////////////////////////////////////////
RowVector interpolation_gradients(const RowVector a,const RowVector b,const double c){
  //interpolate between a (at t=ta) and b (t=tb), at time t=c
  // NOTE: ta and tb are stored in the first elements of a and b
  // rowvectors are from EPI sequence
  RowVector p(8);
  p=0;
  p(1)=c;
  if (fabs(b(1)-a(1))<1e-12) cout<<"Warning: gradients in the EPI sequence are spaced too close to each other"<<endl;
  else{
    p(6)=(c-a(1))*(b(6)-a(6))/(b(1)-a(1))+a(6);
    p(7)=(c-a(1))*(b(7)-a(7))/(b(1)-a(1))+a(7);
    p(8)=(c-a(1))*(b(8)-a(8))/(b(1)-a(1))+a(8);
  }
  return p;
}
/////////////////////////////////////////////////////////////////////////////////////////////////
void sorter(PMatrix& pulse, const Matrix& motion, const string pulsefile,const int opt_test){//ok
  //matrix mainmatrix M  is the output of this function
  //(1)=time (s),
  //(2)=rf angle(rad),(3)=rf frequency bandwidth df(Hz),(4)=rf center frequecy fc(Hz),
  //(5)=readout (1/0),
  //(6)=x gradient (T/m),(7)=y gradient (T/m),(8)=z gradient (T/m),
  //(9)=Tx translation (m),(10)=Ty translation (m), (11)=Tz translation (m), 
  //(12)=b angle of rotation (rad),(13)=Bx,(14)=By,(15)=Bz rotation axis (m)---rotation of the interpolated motion point between A(k) and A(k+1)--relative to A(m)  
  //(16)=a angle of rotation (rad),(17)=Ax,(18)=Ay,(19)=Az rotation axis (m)---rotation at the control motion point A(k)
  
  cout<<"Started the sorter"<<endl;
  int dimp=pulse.Nrows();//size in time for pulse
  int tp=1;// counter for the pulse seq
  int dimm=motion.Nrows();//size in time for motion
  //Counting the size of the new matrix
  int ts=0;
  for (int tm=1;tm<=dimm;tm++){
    if (opt_test==1) cout<<"And here? "<<dimm<<"  tm= "<<tm<<endl;
    if (tp<=dimp){
       if (opt_test==1) cout<<"And here? "<<dimp<<";tp= "<<tp<<" value at tp "<<pulse.time(tp)<<endl;
      while (tp<=dimp && motion(tm,1)-pulse.time(tp)>1e-6){
        ts++;tp++;
      }
      ts++;
      if (opt_test==1) cout<<"tp before equal"<<tp<<"   "<<pulse.time(tp)<<endl;
      if (tp<=dimp && fabs(motion(tm,1)-pulse.time(tp))<1e-06){
	 if (opt_test==1) cout<<"Equal timepoint for pulse and motion: "<<motion(tm,1)<<endl;
	tp++;//cout<<"tp after equal"<<tp<<endl;
      }
    }
  }
  for (int k=tp;k<=dimp;k++){
    ts++;
  }
  int dimnew=ts;
  if (opt_test==1) cout<<"ts= "<<ts<<endl;
  //Storing the order of the rows from pulse and motion into one column. Positive int tp shows where the pulse rows go and negative int tm shows where the new rows go
  if (opt_test==1) cout<<"Creating the tsort vector"<<endl;
  vector<signed int> tsort(dimnew);
  if (opt_test==1) cout<<"tsort was just assigned "<<dimnew<<" pplaces in memory"<<endl;
  ts=0;tp=1;
  for (int tm=1;tm<=dimm;tm++){
    if (tp<=dimp){
      while (tp<=dimp && motion(tm,1)-pulse.time(tp)>1e-6){
	ts++;tsort[ts-1]=tp;tp++;
      }
      ts++;tsort[ts-1]=-tm; 
      if (opt_test==1) cout<<"ts= "<<ts<<"  tp= "<<tp<<endl;
      if (opt_test==1) cout<<"tm elements of tsort= "<<tsort[ts-1]<<endl;
      if (opt_test==1) cout<<"motion time is "<<motion(tm,1)<<"; pulse time is "<<pulse.time(tp)<<endl;
      if (tp<=dimp && fabs(motion(tm,1)-pulse.time(tp))<1e-06){
	tp++;  if (opt_test==1) cout<<"And now tp is "<<tp<<endl;
      }
    }
  }
  if (opt_test==1) cout<<"And tp outside the loop is "<<tp<<endl;
  for (int k=tp;k<=dimp;k++){
    ts++;tsort[ts-1]=k;
  }
  cout<<"Pulse dim="<<dimp<<endl;
  cout<<"Motion dim="<<dimm<<endl; 
  cout<<"Matrix dim="<<dimnew<<endl;
  //for (int k=0;k<=dimnew-1;k++){
  // cout<<tsort[k]<<endl;
  //}
  pulse.destroy();
  pulse.ReSize(dimnew,19);
  pulse=0.0;
  read_binary_matrix(pulse,pulsefile);
  //if the motion points are not matching the pulse seq points
  if (dimnew!=dimp){
     if (opt_test==1) cout<<"Still going ... 1 "<<endl;
    for (int ii=dimnew-1;ii>=0;ii--){
      if (opt_test==1) cout<<"And ii index is "<<ii<<endl;
      if (tsort[ii]>0){
         if (opt_test==1) cout<<"No idea fuck "<<endl;
         if (opt_test==1) cout<<pulse.time(1)<<endl;
         if (opt_test==1) cout<<tsort[0]<<endl;
         if (opt_test==1) cout<<pulse.time(tsort[0])<<endl;
        pulse.time(ii+1)=pulse.time(tsort[ii]);
         if (opt_test==1) cout<<"No idea fuck 2 "<<endl;
        for (int kk=2;kk<=8;kk++){
          if (opt_test==1) cout<<"No idea fuck 3 "<<endl;
	  pulse(ii+1,kk)=pulse(tsort[ii],kk);
         if (opt_test==1) cout<<"No idea fuck  4"<<endl;
        }
      }
      if (opt_test==1){
	cout<<"Still going 2"<<endl;
	cout<<"Accessing element number "<<ii<<endl;
      }
      if (tsort[ii]<0 && ii>=1){
        if (opt_test==1) cout<<"ii is"<<ii<<endl;
        pulse.time(ii+1)=motion(-tsort[ii],1);
	RowVector a(8),b(8),G(8);
	b<<pulse.time(ii+1)<<pulse(ii+1,2)<<pulse(ii+1,3)<<pulse(ii+1,4)<<pulse(ii+1,5)<<pulse(ii+1,6)<<pulse(ii+1,7)<<pulse(ii+1,8);
	a<<pulse.time(ii)<<pulse(ii,2)<<pulse(ii,3)<<pulse(ii,4)<<pulse(ii,5)<<pulse(ii,6)<<pulse(ii,7)<<pulse(ii,8);
	G=interpolation_gradients(a,b,motion(-tsort[ii],1)); 
        for (int kk=2;kk<=8;kk++){
	  pulse(ii+1,kk)=G(kk);
        }
       }
     if (opt_test==1) cout<<"Still going ... 3"<<endl;
    }
   if (opt_test==1) cout<<"Still going ... 4"<<endl;
  }
   if (opt_test==1) cout<<"Finding dimensions of the first tmpvec"<<endl;
  int ii=1; //counter for tsort vector.first element tsort[0] always -1 
  int cc=1; //counter for the new matrix pulse
  while (ii<=dimnew-1 && tsort[ii]>0){
    cc++;ii++;
  }
  int ccdim=cc;
  cc=1;
  if (opt_test==1) cout<<"ccdim="<<ccdim<<endl;
  ColumnVector tmpvec(ccdim+1);
  tmpvec(1)=0;
  ii=1;
  int iii=0;
  while (ii<=dimnew-1){
    if (opt_test==1) cout<<"In the loop: cc="<<cc<<"ccdim= "<<ccdim<<endl;
    if (tsort[ii]<0){
      int m=-tsort[ii];
      cc++;
      tmpvec(cc)=motion(m,1);
      if (opt_test==1) cout<<"motion file= "<<motion(m,1)<<"  "<<motion(m,2)<<"  "<<motion(m,3)<<"  "<<motion(m,4)<<"  "<<motion(m,5)<<"  "<<motion(m,6)<<"  "<<motion(m,7)<<endl;
      //string bla="/Users/ivana/scratch/karolina/testsorter/simdir_new/tmpvec";
      //if (m==4) write_ascii_matrix(tmpvec,bla);
	//cout<<"ii= "<<ii<<"  TmpVec="<<tmpvec<<endl;
      Matrix R=interrotation(motion(m-1,5),motion(m-1,6),motion(m-1,7),motion(m,5),motion(m,6),motion(m,7),tmpvec); 
      for (int kk=1;kk<=8;kk++){
	for (int ll=1;ll<=cc-1;ll++){
	  pulse(iii+ll,kk+11)=R(ll,kk);
	}
      }
      if (opt_test==1) cout<<"After rotation"<<endl;
      Matrix T=intertranslation(motion(m-1,2),motion(m-1,3),motion(m-1,4),motion(m,2),motion(m,3),motion(m,4),tmpvec);
      //if (m==2) write_ascii_matrix(T,bla);
      for (int kk=1;kk<=3;kk++){
	for (int ll=1;ll<=cc-1;ll++){
	  pulse(iii+ll,kk+8)=T(ll,kk);
	}
      }
      if (opt_test==1) cout<<"After translation"<<endl;
      if (opt_test==1) cout<<"m= "<<m<<"; ii= "<<ii<<endl;
      if (m<dimm && ii<dimnew-1){
	iii=ii;
	cc=1;
	if (opt_test==1) cout<<"Values of ii= "<<ii<<" and "<<tsort[ii]<<endl;
	while (ii<=dimnew-1 && tsort[ii+1]>0){
	  if (opt_test==1) cout<<"Values of ii= "<<ii<<endl;
	  cc++;ii++;
	}
	ii=iii;
	ccdim=cc;
	cc=1;
	tmpvec.ReSize(ccdim+1);
	tmpvec(1)=motion(m,1);
	ii++;
      }else{
	if (opt_test==1) cout<<"Are yoiu here?"<<endl;
        for (int kk=1;kk<=3;kk++){
	  int ll=cc;
	  while (iii+ll<=dimnew){
	    pulse(iii+ll,kk+8)=T(cc,kk);
	    ll++;
	  }
	}
	for (int kk=1;kk<=8;kk++){
	  int ll=cc;	 
	  while (iii+ll<=dimnew){
	    pulse(iii+ll,kk+11)=R(cc,kk);
	    ll++;
	  }
	}
	ii=dimnew;
      }
    }else{
      cc++;
      tmpvec(cc)=pulse.time(tsort[ii]);
      //cout<<"Got here? ii= "<<ii<<endl; 
      ii++;
    }
  }
 cout<<"new dimensions of pulse:"<<pulse.Nrows()<<"    "<<pulse.Ncols()<<endl;
//cout<<"new elements:"<<pulse.time(dimp)<<pulse(100,2)<<pulse(100,3)<<pulse(100,4)<<pulse(100,5)<<pulse(100,6)<<pulse(100,7)<<endl;
}

//////////
Matrix sorter_old(const Matrix& epi,const Matrix& motion){//ok
  //matrix mainmatrix M  is the output of this function
  //(1)=time (s),
  //(2)=rf angle(rad),(3)=rf frequency bandwidth df(Hz),(4)=rf center frequecy fc(Hz),
  //(5)=readout (1/0),
  //(6)=x gradient (T/m),(7)=y gradient (T/m),(8)=z gradient (T/m),
  //(9)=Tx translation (m),(10)=Ty translation (m), (11)=Tz translation (m), 
  //(12)=b angle of rotation (rad),(13)=Bx,(14)=By,(15)=Bz rotation axis (m)---rotation of the interpolated motion point between A(k) and A(k+1)--relative to A(m)  
  //(16)=a angle of rotation (rad),(17)=Ax,(18)=Ay,(19)=Az rotation axis (m)---rotation at the control motion point A(k)
  
  int dim1=epi.Nrows();//EPI
  int t1=2;
  int dim2=motion.Nrows();//MOTION
  Matrix mainmatrix(dim1+dim2*2,19);//MAINMATRIX
  mainmatrix=0;
  int t=1;
  
  for (int t2=1;t2<=dim2-1;t2++){
    cout<<"Motion counter= "<<t2<<" till "<<dim2-1<<endl;  
    ColumnVector timevector(dim1+dim2);
    timevector=0;
    int l=1;
    timevector(l)=motion(t2,1);
    while (motion(t2+1,1)-epi(t1,1)>1e-10){
      l=l+1;
      timevector(l)=epi(t1,1);
      t1=t1+1;
    }
    l=l+1;
    timevector(l)=motion(t2+1,1);
    int tmp6=0;
    RowVector G(8);
    G=0;
    if (fabs(motion(t2+1,1)-epi(t1,1))<=1e-10){
      t1=t1+1;
      tmp6=0;
    }
    else{
      G=interpolation_gradients(epi.Row(t1-1),epi.Row(t1),motion(t2+1,1)); 
      tmp6=1;
    }
    Matrix R=interrotation_old(motion(t2,5),motion(t2,6),motion(t2,7),motion(t2+1,5),motion(t2+1,6),motion(t2+1,7),timevector.Rows(1,l)); 
    Matrix T=intertranslation_old(motion(t2,2),motion(t2,3),motion(t2,4),motion(t2+1,2),motion(t2+1,3),motion(t2+1,4),timevector.Rows(1,l));
    for (int tt=1;tt<=l-2;tt++){
      t=t+1;
      RowVector tmp1(19);
      tmp1.Columns(1,8)=epi.Row(t1-l+tt+tmp6);
      tmp1.Columns(9,11)=T.Row(tt);
      tmp1.Columns(12,19)=R.Row(tt);
      mainmatrix.Row(t)=tmp1;
    }
    t=t+1;
    if (tmp6==1){
      RowVector Gpre(8);
      Gpre<<G(1)-1e-10<<0<<0<<0<<0<<G(6)<<G(7)<<G(8);
      mainmatrix.SubMatrix(t,t,1,8)=Gpre;
      mainmatrix.SubMatrix(t+1,t+1,1,8)=G;
    }
    else {
      RowVector Epre(8);
      RowVector E(8);
      Epre<<epi(t1-1,1)-1e-10<<0<<0<<0<<0<<epi(t1-1,6)<<epi(t1-1,7)<<epi(t1-1,8);
      E<<epi(t1-1,1)<<0<<0<<0<<0<<epi(t1-1,6)<<epi(t1-1,7)<<epi(t1-1,8);
      mainmatrix.SubMatrix(t,t,1,8)=Epre;
      mainmatrix.SubMatrix(t+1,t+1,1,8)=E;
    }
    mainmatrix.SubMatrix(t,t+1,9,11)=T.Rows(l-1,l);
    mainmatrix.SubMatrix(t,t+1,12,19)=R.Rows(l-1,l);
    t=t+1;
  }
  Matrix tmp2(1,11);
  tmp2=mainmatrix.SubMatrix(t,t,9,19); 
  for (int a=1;a<=dim1-t1+1;a++){
    mainmatrix.SubMatrix(t+a,t+a,9,19)=tmp2;
  }
  mainmatrix.SubMatrix(t+1,dim1-t1+t+1,1,8)=epi.SubMatrix(t1,dim1,1,8);
  mainmatrix=mainmatrix.Rows(1,dim1-t1+t+1); 
  cout<<"main matrix dim"<<dim1-t1+t+1<<"  .  "<<"epi matrix dim"<<dim1<<"....."<<endl;
  return mainmatrix;
}
//////////////////
//B0 FIELD STUFF
///////////////////////////////////////////////////////////////////////////////
int calc_gradients(const volume<double>& b, volume<double>& b0gx, volume<double>& b0gy, volume<double>& b0gz){
  // b is in mT, and dimensions of voxels are in mm 
  b0gx = b*0;
  b0gy = b0gx;
  b0gz = b0gx;
  for (int z=1;z<b.zsize()-1; z++) {
    for (int y=1; y<b.ysize()-1; y++) {
      for (int x=1; x<b.xsize()-1; x++) {
	b0gx(x,y,z) = 1*(b(x+1,y+1,z+1) + b(x+1,y-1,z+1) + b(x+1,y-1,z-1) + b(x+1,y+1,z-1) - b(x-1,y+1,z+1) - b(x-1,y-1,z+1) - b(x-1,y-1,z-1) - b(x-1,y+1,z-1)) + 6*(b(x+1,y,z+1) + b(x+1,y,z-1) + b(x+1,y+1,z) + b(x+1,y-1,z)- b(x-1,y,z+1) - b(x-1,y,z-1) - b(x-1,y+1,z) - b(x-1,y-1,z))+ 36*(b(x+1,y,z) - b(x-1,y,z));
	b0gy(x,y,z) = 1*(b(x+1,y+1,z+1) + b(x-1,y+1,z+1) + b(x-1,y+1,z-1) + b(x+1,y+1,z-1) - b(x+1,y-1,z+1) - b(x-1,y-1,z+1) - b(x-1,y-1,z-1) - b(x+1,y-1,z-1))+ 6*(b(x,y+1,z+1) + b(x,y+1,z-1) + b(x+1,y+1,z) + b(x-1,y+1,z)- b(x,y-1,z+1) - b(x,y-1,z-1) - b(x+1,y-1,z) - b(x-1,y-1,z)) + 36*(b(x,y+1,z) - b(x,y-1,z));
	b0gz(x,y,z) = 1*(b(x+1,y+1,z+1) + b(x-1,y+1,z+1) + b(x-1,y-1,z+1)+ b(x+1,y-1,z+1) - b(x+1,y+1,z-1)- b(x-1,y+1,z-1) - b(x-1,y-1,z-1) - b(x+1,y-1,z-1))+ 6*(b(x,y+1,z+1) + b(x,y-1,z+1) + b(x+1,y,z+1) + b(x-1,y,z+1)- b(x,y+1,z-1) - b(x,y-1,z-1) - b(x+1,y,z-1) - b(x-1,y,z-1))+ 36*(b(x,y,z+1) - b(x,y,z-1));
      }
    }
  }
  b0gx /= (128.0*b.xdim()*1e-3);  // in mT/m
  b0gy /= (128.0*b.ydim()*1e-3);
  b0gz /= (128.0*b.zdim()*1e-3);
  
  return 0;
}
///////////////////////////////////////////////////////////////
int calc_gradientsROI(volume<double>& b0, volume<double>& b0x, 
                      volume<double>& b0y, volume<double>& b0z,
                      const int myid, const int Nxx,const int numprocs){
  calc_gradients(b0,b0x,b0y,b0z);
  int Nx=b0.xsize();
  int Ny=b0.ysize();
  int Nz=b0.zsize();
  int Ntmp=(int)((Nxx-myid-1)/numprocs)+1;
  if (Ntmp<1){
   cout<<"WARNING:Number of processors bigger than the number of voxels in x-direction."<<endl;
   exit(EXIT_FAILURE);
 }
  volume<double> tmpvol1(Ntmp,Ny,Nz);
  volume<double> tmpvol2(Ntmp,Ny,Nz);
  volume<double> tmpvol3(Ntmp,Ny,Nz);
  volume<double> tmpvol4(Ntmp,Ny,Nz);
  int x0=0, x1=0, y0=0, y1=0, z0=0, z1=0;
  for (z0=0, z1=0; z0<Nz; z0++, z1++) {
    for (y0=0, y1=0; y0<Ny; y0++, y1++) {
      for (x0=myid, x1=0; x0<Nx; x0+=numprocs, x1++) {
	tmpvol1(x1,y1,z1) = b0(x0,y0,z0);
        tmpvol2(x1,y1,z1) = b0x(x0,y0,z0);
        tmpvol3(x1,y1,z1) = b0y(x0,y0,z0);
        tmpvol4(x1,y1,z1) = b0z(x0,y0,z0);
      }
    }
   }
   b0= tmpvol1;
   b0x= tmpvol2;
   b0y= tmpvol3;
   b0z= tmpvol4;

   return 0;
}
///////////////////////////////////////////////////////////////////////////////
int calc_gradients4D(const volume4D<double>& b, volume4D<double>& b0gx, volume4D<double>& b0gy, volume4D<double>& b0gz){
  // b is in mT, and dimensions of voxels are in mm 
  b0gx = b*0;
  b0gy = b0gx;
  b0gz = b0gx;
 for (int t=0;t<b.tsize(); t++) {
  for (int z=1;z<b.zsize()-1; z++) {
    for (int y=1; y<b.ysize()-1; y++) {
      for (int x=1; x<b.xsize()-1; x++) {
	b0gx(x,y,z,t) = 1*(b(x+1,y+1,z+1,t) + b(x+1,y-1,z+1,t) + b(x+1,y-1,z-1,t) + b(x+1,y+1,z-1,t) - b(x-1,y+1,z+1,t) - b(x-1,y-1,z+1,t) - b(x-1,y-1,z-1,t) - b(x-1,y+1,z-1,t)) + 6*(b(x+1,y,z+1,t) + b(x+1,y,z-1,t) + b(x+1,y+1,z,t) + b(x+1,y-1,z,t)- b(x-1,y,z+1,t) - b(x-1,y,z-1,t) - b(x-1,y+1,z,t) - b(x-1,y-1,z,t))+ 36*(b(x+1,y,z,t) - b(x-1,y,z,t));
	b0gy(x,y,z,t) = 1*(b(x+1,y+1,z+1,t) + b(x-1,y+1,z+1,t) + b(x-1,y+1,z-1,t) + b(x+1,y+1,z-1,t) - b(x+1,y-1,z+1,t) - b(x-1,y-1,z+1,t) - b(x-1,y-1,z-1,t) - b(x+1,y-1,z-1,t))+ 6*(b(x,y+1,z+1,t) + b(x,y+1,z-1,t) + b(x+1,y+1,z,t) + b(x-1,y+1,z,t)- b(x,y-1,z+1,t) - b(x,y-1,z-1,t) - b(x+1,y-1,z,t) - b(x-1,y-1,z,t)) + 36*(b(x,y+1,z,t) - b(x,y-1,z,t));
	b0gz(x,y,z,t) = 1*(b(x+1,y+1,z+1,t) + b(x-1,y+1,z+1,t) + b(x-1,y-1,z+1,t)+ b(x+1,y-1,z+1,t) - b(x+1,y+1,z-1,t)- b(x-1,y+1,z-1,t) - b(x-1,y-1,z-1,t) - b(x+1,y-1,z-1,t))+ 6*(b(x,y+1,z+1,t) + b(x,y-1,z+1,t) + b(x+1,y,z+1,t) + b(x-1,y,z+1,t)- b(x,y+1,z-1,t) - b(x,y-1,z-1,t) - b(x+1,y,z-1,t) - b(x-1,y,z-1,t))+ 36*(b(x,y,z+1,t) - b(x,y,z-1,t));
      }
    }
  }
 }
  b0gx /= (128.0*b.xdim()*1e-3);  // in mT/m
  b0gy /= (128.0*b.ydim()*1e-3);
  b0gz /= (128.0*b.zdim()*1e-3);
  
  return 0;
}
///////////////////////////////////////////////////////////////
int calc_gradients4DROI(volume4D<double>& b0, volume4D<double>& b0x, 
                      volume4D<double>& b0y, volume4D<double>& b0z,
			const int myid, const int Nxx, const int numprocs){
  calc_gradients4D(b0,b0x,b0y,b0z);
  int Nx=b0.xsize();
  int Ny=b0.ysize();
  int Nz=b0.zsize();
  int Nt=b0.tsize();
  volume4D<double> tmpvol1(Nx,Ny,Nz,Nt);
  volume4D<double> tmpvol2(Nx,Ny,Nz,Nt);
  volume4D<double> tmpvol3(Nx,Ny,Nz,Nt);
  volume4D<double> tmpvol4(Nx,Ny,Nz,Nt);
  int x0=0, x1=0, y0=0, y1=0, z0=0, z1=0, t0=0, t1=0;
  for (t0=0, t1=0; t0<Nt; t0++, t1++) {
    for (z0=0, z1=0; z0<Nz; z0++, z1++) {
      for (y0=0, y1=0; y0<Ny; y0++, y1++) {
        for (x0=myid, x1=0; x0<Nxx; x0+=numprocs, x1++) {
	  tmpvol1(x1,y1,z1,t1) = b0(x0,y0,z0,t0);
          tmpvol2(x1,y1,z1,t1) = b0x(x0,y0,z0,t0);
          tmpvol3(x1,y1,z1,t1) = b0y(x0,y0,z0,t0);
          tmpvol4(x1,y1,z1,t1) = b0z(x0,y0,z0,t0);
	}
      }
    }
  }
   tmpvol1.copyproperties(b0);
   tmpvol2.copyproperties(b0x);
   tmpvol3.copyproperties(b0y);
   tmpvol4.copyproperties(b0z);
   b0= tmpvol1;
   b0x= tmpvol2;
   b0y= tmpvol3;
   b0z= tmpvol4;

   return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double b0int(const double v1i,const double v1j,const double v2i, const double v2j,const double v3i,const double v3j,
             const double aold ,const double anew,const double told,const double tnew){
  //integral is  p1*\int(1)+p2*\int(sin(a1+a2*t))-p3*\int(cos(a1+a2*t))+p4*\int(sin^2(a1+a2*t))-p5*\int(sin(a1+a2*t)*cos(a1+a2*t))
  double val;
  double p1=v1i*v1j+v1i*v3j+v3i*v1j+2*v3i*v3j;
  double p2=v1i*v2j+v2i*v1j+v3i*v2j+v2i*v3j;
  double p3=v1i*v3j+2*v3i*v3j+v3i*v1j;
  double p4=v2i*v2j-v3i*v3j;
  double p5=v2i*v3j+v3i*v2j;
  double a1,a2;
  a1=0.0;a2=0.0;
  coeff(aold,anew,told,tnew,a1,a2);
  if (fabs(a2)*(tnew-told)<0.01){
    val=(p1+p4/2+(p2*sin(a1)-p3*cos(a1))*cos(a2*(tnew+told)/2)+(p2*cos(a1)+p3*sin(a1))*sin(a2*(tnew+told)/2)-(p4*cos(2*a1)+p5*sin(2*a1))*cos(a2*(tnew+told))/2+(p4*sin(2*a1)-p5*cos(2*a1))*sin(a2*(tnew+told))/2)*(tnew-told);
    //val=(p1+p2*sin(a1)-p3*cos(a1)-p4*sin(a1)*sin(a1)-p5*cos(a1)*sin(a1))*(tnew-told);
  }
  else {
    val=(p1+p4/2)*(tnew-told)+((p2*sin(a1)-p3*cos(a1))*cos(a2*(tnew+told)/2)*sin(a2*(tnew-told)/2)*2+(p2*cos(a1)-p3*sin(a1))*sin(a2*(tnew+told)/2)*sin(a2*(tnew-told)/2)*2-(p4*cos(2*a1)+p5*sin(2*a1))*cos(a2*(tnew+told))*sin(a2*(tnew-told))/2+(p4*sin(2*a1)-p5*cos(2*a1))*sin(a2*(tnew+told))/2*sin(a2*(tnew-told)))/a2;
  if (fabs(a2)<0.00001) cout<<"WARNING check: a2="<<a2<<"val"<<val<<endl;
    // val=p1*(tnew-told)-(p2*(cos(a1+a2*tnew)-cos(a1+a2*told))-p3*(sin(a1+a2*tnew)-sin(a1+a2*told))-p4*(sin(2*(a1+a2*tnew))-sin(2*(a1+a2*told)))/4+p5*(cos(2*(a1+a2*tnew))-cos(2*(a1+a2*told)))/4)/a2;
  }
    return val;
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double b0freq(const double v1i,const double v1j,const double v2i, const double v2j,const double v3i,const double v3j,
              const double anew){
  double val;
  double p1=v1i*v1j+v1i*v3j+v3i*v1j+3/2*v3i*v3j+1/2*v2i*v2j;
  double p2=v1i*v2j+v2i*v1j+v3i*v2j+v2i*v3j;
  double p3=v1i*v3j+v3i*v3j+v3i*v1j;
  double p4=v2i*v2j-v3i*v3j;
  double p5=v2i*v3j+v3i*v2j;
  val=p1+p2*sin(anew)-p3*cos(anew)-p4*sin(anew)*sin(anew)-p5*cos(anew)*sin(anew);
  return val;
}
/////////////////
//MAIN FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void voxel1(const double x,const double y,double z, 
            const RowVector& tissue,const PMatrix& H,
            const int nreadp,const int v,
            const double xdim,const double ydim,const double zdim,
            const double b0, const double b0x,const double b0y,const double b0z,
            const double* timecourse, const double* activation,const int Nact,
	    const string outputname,  
	    const double* table_slcprof, const double dslcp, const double dslcp_first, const int Nslc,
            const double den,const double RFtrans,const int opt_test,
            const int nospeedup,
            const int save_kcoord,
            double* sreal, double* simag)
{
  //  Returns a signal (one rowvector for the real part and one for the 
  //  imaginary part) for one voxel at place x,y,z(m) of the phantom through time
  // - x,y,z are the corrdinates at the center of the voxel
  // - tissue = (T_1(s),T_2(s),\rho) ; H = event matrix (pulse sequence)
  // - nreadp = number of readoutpoints !!! pay attention as you are not using this at all here;  
  // - v = internal voxel index number (used to check if it is the first voxel)
  // - xdim, ydim and zdim are dimensions of the voxel
  // - b0 is inhomogeneity in the field that can be due to either b0sus or chemical shift, it is in T
  // - b0x,b0y,b0z are values of gradients of b0 inhomogeneity field for each voxel (T/m) 
  // - timecourse (time when happens), activation (amount of change)and Nact (number of changes) define change in T2*  
  // - outputname is the name of the result once the programm finishes 
  // - table_slcprof is the vector specifying the sclive profile for the RF excitation
  // - dslcp is the timestep in the the slice profile
  // - dslcp_first is the first element in the slice profile
  // - Nslc is the number of elements in the slice profile
  // - den=phantom(xx,yy,zz,tt)*RFrec(xx,yy,zz)*xdim*ydim*zdim;
  // - RFtrans is the inhomogeneity in the RF transmit coil. Influences the flip angle : FAnew=FAold*RFtrans (1 is perfc. hom)
  // - if opt_test=1 then write out all the help lines. For debugging.
  // - zz_slc is the physical slice to which this voxel should belong to in ideal conditions
  // - numslc is the total number of slices to be imaged
  // - extra_down_slc and extra_up_slc show how many more extra slices are involved due to motion and slice profile
  // - Nrfpslc is number of RF pulses per slice acquisition: 1 for epi and resY for ge
  // - save_kcoord is yes or no option if we want to save the coordinates of the k-space as an additonal output
  // - sreal and simag are two channels of the signal ouput

///////////////////////////////////////////////////////////////////////////
//INITIALIZATION OF THE MAGNETIZATION VECTOR 
///////////////////////////////////////////////////////////////////////////

	ColumnVector m(3);	//magnetization vector
	m(1)=0;
	m(2)=0;
	m(3)=tissue(3);
	double m00=0;		//the magnitude of the transverse magnetization vector

/////////////////////////////////////////////////////////////////////////// 
//PHASE,GRADIENTS,COUNTERS,CONSTANTS
///////////////////////////////////////////////////////////////////////////

	double chshift = tissue(4);	//chemical shift
	int numpoints = H.Nrows();	//number of timepoints in the eventsequencer matrix
	int readstep = 0;		//keeps track of readout points
	int excitation = 0;
	double trf = 0;			//rftime
	double T2 = tissue(2);
	double iT2 = 1/T2;
	int rftest = 0;
	double g1,g2,g3,grf1,grf2,grf3; 
	g1=0.0; g2=0.0; g3=0.0; grf1=0.0; grf2=0.0; grf3=0.0;
	Matrix coord(3,nreadp);

	if ( v==1 )
	{
		g1static = new double[numpoints];
		g2static = new double[numpoints];
		g3static = new double[numpoints];

		////////////////////////////////////////////////////////////////////////
		//LOOK UP TABLES
		////////////////////////////////////////////////////////////////////////
		table_sinc = new double[glo_Nsinc+1];
		table_sin = new double[glo_Nsin+1];
		table_cos = new double[glo_Nsin+1];
		if (opt_test==1)
		{
			cout<<"Stepsize for Table for SINC: dsinc= "<<glo_dsinc<<endl;
			cout<<"Stepsize for Table for SIN: dsin="<<glo_dsin<<endl;
		}

		for (int n=0;n<=glo_Nsinc;n++)
			table_sinc[n]=Sinc(n*glo_dsinc);

		for (int n=0;n<=glo_Nsin;n++)
		{
			table_sin[n]=sin(n*glo_dsin);
			table_cos[n]=cos(n*glo_dsin);
		}

		glo_cx=gammabar*xdim;	//(Hz*m/T)
		glo_cy=gammabar*ydim;
		glo_cz=gammabar*zdim;
	} 

///////////////////////////////////////////////////////////////////////////
//  ACTIVATION PARAMETERS
///////////////////////////////////////////////////////////////////////////

	double actint = 0.0;
	int actstep = 0;
	double dT2_1 = 0.0;
	double dT2_2 = 0.0;

///////////////////////////////////////////////////////////////////////////
//SIGNAL
/////////////////////////////////////////////////////////////////////////// 
	for (int step=2; step <= numpoints; step++)
	{
		double tnew = H.time(step);
		double told = H.time(step-1);
		double rfangle = H(step,2);
		double read = H(step,5);
		double gxnew = H(step,6);
		double gynew = H(step,7);
		double gznew = H(step,8);
		if (v==1)
		{
			double gxold = H(step-1,6);
			double gyold = H(step-1,7);
			double gzold = H(step-1,8);

			g1 += i1(gxold,gxnew,told,tnew); //i1 - integral function
			if (gynew!=0 || gyold!=0)	g2+=i1(gyold,gynew,told, tnew);
			if (gznew!=0 || gzold!=0)	g3+=i1(gzold,gznew,told, tnew);

			g1static[step-2] = g1;
			g2static[step-2] = g2;
			g3static[step-2] = g3;
		}
		else
		{
			g1 = g1static[step-2];
			g2 = g2static[step-2];
			g3 = g3static[step-2];
		}

		double gg1 = g1 - grf1;
		double gg2 = g2 - grf2;
		double gg3 = g3 - grf3;

		double tt = tnew - trf;//time since the last rf pulse
		cout.precision(20);

//		if (v==1 && readstep%4096==1 && opt_test==1)cout<<"tnew "<<tnew<<";told "<<told<<";trf "<<trf<<";tt=tnew-trf "<<tt<<";tnew-told "<<tnew-told<<endl; 
//		if (v==1 && readstep%4096==1 && opt_test==1)cout<<"g1 "<<g1<<";grf1 "<<grf1<<";gg1 "<<gg1<<endl;

		if( told >= timecourse[actstep] && actstep <= (Nact-2) )
		{
			coeff(activation[actstep],activation[actstep+1],timecourse[actstep],timecourse[actstep+1],dT2_1,dT2_2);
			dT2_1 = dT2_1*iT2*iT2;
			dT2_2 = dT2_2*iT2*iT2;
			actstep = actstep+1;
		}
		actint += (dT2_1+dT2_2*(tnew + told)/2)*(tnew - told);
//		if (dT2_2<1e-12) actint+=(tnew-told)/(dT2_1+T2);
//		else actint+=log(1+(tnew-told)/(told+dT2_1/dT2_2))/dT2_2
	
//		double phase = 0;
		double phase = gama*gg1*x + gama*gg2*y + gama*gg3*z+gama*b0*tt + gama*chshift*tt; 

		if (rfangle != 0)
		{
			excitation = 0;
			double df = H(step,3);
			double fc = H(step,4);
			double f = gammabar*(gxnew*x+gynew*y+gznew*z+b0+chshift);

			if (opt_test==1)	cout<<"gxnew="<<gxnew<<"; x="<<x<<"; gynew="<<gynew<<"; y="<<y<<"; gznew="<<gznew<<"; z="<<z<<"; fc="<<fc<<"df="<<df<<". For voxel v="<<v<<": f="<<f<<endl;
			rftest = 0;
			double fval = (f - fc)/df;
			double off = (fval - dslcp_first)/dslcp;
			int nf = (int) off;

//			if (opt_test==1)	cout<<"fc="<<fc<<"; f"<<f<<"; df"<<df<<"; fval"<<fval<<"; dslcp_first"<<dslcp_first<<"; dslcp"<<dslcp<<"; off="<<off<<"nf"<<nf<<endl;
			if (nf >= 0 && nf <= (Nslc-2))
			{
				off -= nf;
				double ts = table_slcprof[nf];
				double sx = (table_slcprof[nf+1]-ts)*off + ts;
				double rfangle_f = sx*rfangle*RFtrans;//RFtrans are values 0 to 1 to derscribe the inhomogeneity f the receive RF field. 1 is for perfectly homog;
				if (opt_test==1 && readstep%4096==0 && v==1)
				{
					cout<<"table_slcprof[nf]= "<<table_slcprof[nf]<<"; rfangle= "<<rfangle<<"; RFtrans= "<<RFtrans<<endl;
					cout<<"rfangle_f=table_slcprof[nf]*rfangle*RFtrans= "<<rfangle_f<<endl;
				}

				if (rfangle_f>0)	//new stuff mon dec 19
				{

					excitation = 1;
					m = free(m,tt,tissue,phase,actint);
//new stuff mon dec 19
//due to crushers or any gradient induced dephasing over the voxel a new initial magnetisation is introduced which is the average of the magnetisations over the voxel  

					double xvalrf = fabs(glo_cx*(gg1 + b0x*tt));
					double yvalrf = fabs(glo_cy*(gg2 + b0y*tt));
					double zvalrf = fabs(glo_cz*(gg3 + b0z*tt));
					double xyzrf = Sinc(xvalrf)*Sinc(yvalrf)*Sinc(zvalrf);

					m(1) = m(1)*xyzrf;
					m(2) = m(2)*xyzrf;
					m = rot(rfangle_f,"x")*m;
					//new stuff
					m00 = sqrt( m(1)*m(1)+m(2)*m(2) );

					if (opt_test==1 && readstep%4096==0 && v==1 )	cout<<"Projection of the magnetisation vector into the xy plane after flipping is "<<m00<<endl;

					trf = tnew;
					grf1 = g1;
					grf2 = g2;
					grf3 = g3;


					actint = 0.0;
					rftest = 1;

				} //new stuff
			}
		}

		if (read != 0)
		{
			readstep = readstep+1;
			if (excitation == 1 || nospeedup == 1)
			{
//				phase = gama*gg1*x + gama*gg2*y + gama*gg3*z+gama*b0*tt + gama*chshift*tt; 

				double xval = fabs(glo_cx*(gg1 + b0x*tt));
				double yval = fabs(glo_cy*(gg2 + b0y*tt));
				double zval = fabs(glo_cz*(gg3 + b0z*tt));
				if (v==1 && save_kcoord==1)
				{
					coord(1,readstep) = gammabar*(gg1 + b0x*tt);//to record the distorted coordinatres in the k-space
					coord(2,readstep) = gammabar*(gg2 + b0y*tt);
					coord(3,readstep) = gammabar*(gg3 + b0z*tt);
				}
				double tmp;

	#ifdef NOTABLE
				if (v==1 && readstep==1)	cout<<"No table sinc calc"<<endl;

				tmp=m00*exp(-tt*iT2+actint)*Sinc(xval)*Sinc(yval)*Sinc(zval);

				sreal[readstep-1] += den*tmp*cos(phase);
				simag[readstep-1] += den*tmp*sin(phase);

	#else
				if (v==1 && readstep==1)	cout<<"Table sinc calc"<<endl;
				if (xval<=glo_Dsinc && yval<=glo_Dsinc && zval<=glo_Dsinc)
				{
				//TABLES for SINC
				//x
					double off = xval*glo_idsinc;
					int nbin = (int) off;
					off-=nbin;
					double ts = table_sinc[nbin];
					double sx = (table_sinc[nbin+1]-ts)*off + ts;
				//y
					off = yval*glo_idsinc;
					nbin = (int) off;
					off-=nbin;
					ts = table_sinc[nbin];
					double sy = (table_sinc[nbin+1]-ts)*off + ts;
				//z        
					off = zval*glo_idsinc;
					nbin = (int) off;
					off-=nbin;
					ts = table_sinc[nbin];
					double sz = (table_sinc[nbin+1]-ts)*off + ts;
				//CALCULATING TMP
					tmp = m00*exp(-tt*iT2+actint)*sx*sy*sz;
				//TESTING 
					if (opt_test==1 && v==1 && readstep%4096==2081)
					{
						cout.precision(20);
						cout<<"table_sinc(xval)= "<<sx<<"; table_sinc(yval)= "<<sy<<"; table_sinc(zval)= "<<sz<<endl;
						cout<<"sinc-table(xval)= "<<Sinc(xval)-sx<<"; sinc-table(yval)= "<<Sinc(yval)-sy<<"; sinc-tabel(zval)= "<<Sinc(zval)-sz<<endl;
					}
				}
				else
				{
					tmp = m00*exp(-tt*iT2+actint)*Sinc(xval)*Sinc(yval)*Sinc(zval);
				}
				//TABLES SIN AND COS CALCULATION

				double phase_2pi;
				if (phase>0)
					phase_2pi = phase - ((int) (phase*glo_itwopi))*glo_twopi;//one solution when phase exceedes Dsin , this is faster
				else
			 		phase_2pi = phase - ((int) (phase*glo_itwopi))*glo_twopi+glo_twopi;

				double off = phase_2pi*glo_idsin;
				int nphase = (int) off;		off -= nphase;
				double ts1 = table_sin[nphase],	tc1 = table_cos[nphase];

				double wanted_cos = (table_cos[nphase+1] - tc1)*off+tc1;
				double wanted_sin = (table_sin[nphase+1] - ts1)*off+ts1;

				//SIGNAL CALCULATION
				sreal[readstep-1] += den*tmp*wanted_cos;
				simag[readstep-1] += den*tmp*wanted_sin;

				//TESTING
				if(opt_test==1 && v==1 && readstep%16384==8192){
				{
				          cout.precision(20);
				          cout<<"readstep= "<<readstep<<endl;
				          cout<<"gama= "<<gama<<endl;
				          cout<<"gg1= "<<gg1<<"; gg2= "<<gg2<<"; gg3= "<<gg3<<endl;
				          cout<<"x ="<<x<<"; y= "<<y<<"; z= "<<z<<endl;
				          cout<<"b0= "<<b0<<"; chshift= "<<chshift<<endl;
				          cout<<"tnew= "<<tnew<<"; trf= "<<trf<<endl;
				          cout<<"tt=tnew-trf= "<<tt<<endl;
				          cout<<"phase=gama*(gg1*x+gg2*y+gg3*z+(b0+chshift)*tt)= "<<phase<<endl; 
					  cout<<"phase_2pi= "<<phase_2pi<<endl;
					  cout<<"table_cos(phase)= "<<wanted_cos<<"; cos-table(phase)= "<<cos(phase)-wanted_cos<<endl;
				          cout<<"table_sin(phase)= "<<wanted_sin<<"; sin-table(phase)= "<<sin(phase)-wanted_sin<<endl;
					  cout<<"glo_cx= "<<glo_cx<<"; glo_cy= "<<glo_cy<<"; glo_cz= "<<glo_cz<<endl;
					  cout<<"b0x= "<<b0x<<"; b0y= "<<b0y<<"; b0z= "<<b0z<<endl;
				          cout<<"Ival=fabs(glo_cI*(ggI+b0I*tt)); I=x,y,z"<<endl;
				          cout<<"xval= "<<xval<<"; yval= "<<yval<<"; zval= "<<zval<<endl;
				          cout<<"m00= "<<m00<<endl;
				          cout<<"actint= "<<actint<<"; iT2= "<<iT2<<" so exp(-tt*iT2+actint)= "<<exp(-tt*iT2+actint)<<endl;
					  cout<<"tmp= m00*exp(-tt*iT2+actint)*Sinc(xval)*Sinc(yval)*Sinc(zval)= "<<tmp<<endl;
					  cout<<"den= "<<den<<endl;
				          cout<<"sreal_1voxel("<<readstep<<")= den*tmp*table_cos[nphase]= "<<sreal[readstep-1]<<endl;
				          cout<<"simag_1voxel("<<readstep<<")= den*tmp*table_sin[nphase]= "<<simag[readstep-1]<<endl;
					  cout<<"-------------------------------------------------------------------------------"<<endl;
					}
				}
#endif
			}	//end of nospeedup
		}	//end of if read
	}	//end of main loop

	if (v==1 && save_kcoord==1)
	{
		write_binary_matrix(coord,outputname+"_kcoord" );
	}
}

/////////////////////////////////////////////////////////////////////////////

void voxel2(const double x,const double y,const double z, 
		      const RowVector& tissue,const PMatrix& H,
               	      const int nrf,const int nreadp,const int v,                               
		      const double xdim,const double ydim,const double zdim,
                      const double b0, const double b0x,const double b0y,const double b0z,
	              const double* timecourse, const double* activation,const int Nact,
	              const string outputname, 
	              const double* table_slcprof,const double dslcp, const double dslcp_first, const int Nslc,
                      const double den,const double RFtrans,const int opt_test,
	              const int nospeedup,
                      const int save_kcoord,
                      double* sreal, double* simag) {
  // Returns a signal (one rowvector for the real part and one for the imaginary part) for one voxel at place x,y,z of the phantom through time
  // tissue = (T_1,T_2,\rho) ; H = event matrix (pulse sequence + motion)
  // nrf=number of rf pulses ; v = internal voxel index number (used to check if this is the first voxel)
  // xdim, ydim and zdim are dimensions of the voxel
  // timecourse is a 2-column matrix where time is in first and multiplicative factor for the activation const "activation" in the second
  // sreal and simag are the outputs 
  //////////////////////////////////////////////////////////////////////////
  //INITIALIZATION OF THE MAGNETIZATION VECTOR
  //////////////////////////////////////////////////////////////////////////
  ColumnVector m(3);//magnetization vector: M
  m(1)=0;
  m(2)=0;
  m(3)=tissue(3);
  double m00=0; // magnitude of M_{xy}
  /////////////////////////////////////////////////////////////////////////
  //PHASE,GRADIENTS,COUNTERS,CONSTATNTS
  ////////////////////////////////////////////////////////////////////////
  double chshift=tissue(4);//chemical shift
  int numpoints=H.Nrows();//number of timepoints in the eventsequencer matrix
  int readstep=0;//keeps track of readout points
  int rfstep=0;//keeps track of rf pulses
  int excitation=0;
  //int rftest1=0;
  //int rftest2=0;
  double trf=0;//rftime
  double T2=tissue(2);
  double iT2=1/T2;
  double g1,g2,g3,g4;
  double grf1=0.0;
  double grf2=0.0;
  double grf3=0.0;
  double grf4=0.0;
  double rr1,rr2,rr3,trr;
  Matrix g(4,3); //integrated gradients for the first voxel starting from t=0 until tnew
  g=0.0;
  RowVector rnew(4),rmnew(4);//angle, axis values
  double trnew1=0;//translation
  double trnew2=0;
  double trnew3=0;
  double gxnew=0;//gradient strength at tnew
  double gynew=0;
  double gznew=0;
  Matrix coord(3,nreadp);
  if (v==1) {
    g1motion.resize(numpoints,0.0);
    g2motion.resize(numpoints,0.0);
    g3motion.resize(numpoints,0.0);
    g4motion.resize(numpoints,0.0);
    rotmotion1.resize(nrf,0.0);
    rotmotion2.resize(nrf,0.0);
    rotmotion3.resize(nrf,0.0);
    transmotion.resize(nrf,0.0);
     ////////////////////////////////////////////////////////////////////////
     //LOOK UP TABLES
     ///////////////////////////////////////////////////////////////////////
     table_sinc = new double[glo_Nsinc+1];
     table_sin = new double[glo_Nsin+1];
     table_cos = new double[glo_Nsin+1];
     if (opt_test==1){
       cout<<"Stepsize for Table for SINC: dsinc= "<<glo_dsinc<<endl;
       cout<<"Stepsize for Table for SIN: dsin="<<glo_dsin<<endl;
     }
     for (int n=0;n<=glo_Nsinc;n++){
       table_sinc[n]=Sinc(n*glo_dsinc);
     }
     for (int n=0;n<=glo_Nsin;n++){
       table_sin[n]=sin(n*glo_dsin-glo_Dsin);
       table_cos[n]=cos(n*glo_dsin-glo_Dsin);
     }
     glo_cx=gammabar*xdim; //(Hz*m/T)
     glo_cy=gammabar*ydim;
     glo_cz=gammabar*zdim;
  } 
  //we assume that for the step=1 all the values in the mainmatrix H are zero
  //////////////////////////////////////////////////////////////////////////
  //TESTING
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //  ACTIVATION PARAMETERS
  //////////////////////////////////////////////////////////////////////////
  double actint=0.0;
  int actstep=0;
  double dT2_1=0.0;
  double dT2_2=0.0;
  //////////////////////////////////////////////////////////////////////////
  //SIGNAL
  //////////////////////////////////////////////////////////////////////////
  //MAIN LOOP
  for (int step=2;step<=numpoints;step++){
    double told=H.time(step-1);//time at the previous step
    double tnew=H.time(step);//time at this step
    double rfangle=H(step,2);//flip angle
    double aold=H(step-1,12);//angle from last given point to the interpolated one (they have the same axis of rotation) at told
    double anew=H(step,12);//angle from last given point to the interpolated one (they have the same axis of rotation) at tnew
    double read=H(step,5);
    if (v==1){
      rnew<<H(step,12)<<H(step,13)<<H(step,14)<<H(step,15);//rotation betwen given motion points at tnew (angle & axis)
      rmnew<<H(step,16)<<H(step,17)<<H(step,18)<<H(step,19);//rotation at control motion point at tnew (angle & axis)
      double gxold=H(step-1,6);//gradients at told 
      double gyold=H(step-1,7);//
      double gzold=H(step-1,8);//
      gxnew=H(step,6);// gradients at tnew
      gynew=H(step,7);//
      gznew=H(step,8);//
      double trold1=H(step-1,9);//translation at told
      double trold2=H(step-1,10);//
      double trold3=H(step-1,11);
      trnew1=H(step,9);//translation at tnew
      trnew2=H(step,10);//
      trnew3=H(step,11);//
        g(1,1)+=i1(gxold,gxnew,told, tnew);
        g(2,1)+=i2(gxold,gxnew,aold,anew,told,tnew);
        double a1=0.0;
        double a2=0.0;
        coeff(aold,anew,told,tnew,a1,a2);
        g(3,1)+=i1(gxold,gxnew,told, tnew)- i3(gxold,gxnew,aold,anew,told,tnew);
	
        g(4,1)+=i4(gxold,gxnew,trold1,trnew1,told,tnew);
        if(opt_test==1 && readstep%4096==2081 && v==1){
	  cout<<"i4(gxold,gxnew,trold1,trnew1,told,tnew)"<<i4new(gxold,gxnew,trold1,trnew1,told,tnew)<<endl;
	}

      if (gyold!=0 || gynew!=0){
        g(1,2)+=i1(gyold,gynew,told,tnew);
        g(2,2)+=i2(gyold,gynew,aold,anew,told,tnew);

        g(3,2)+=i1(gyold,gynew,told,tnew)- i3(gyold,gynew,aold,anew,told,tnew);

        g(4,2)+=i4(gyold,gynew,trold2,trnew2,told,tnew);
      }
      if (gzold!=0 || gznew!=0){
        g(1,3)+=i1(gzold,gznew,told,tnew);
        g(2,3)+=i2(gzold,gznew,aold,anew,told,tnew);

        g(3,3)+=i1(gzold,gznew,told,tnew)- i3(gzold,gznew,aold,anew,told,tnew);

        g(4,3)+=i4(gzold,gznew,trold3,trnew3,told,tnew);
      }
      //Matrix R(3,3);
      //Matrix A(3,3);
      //R=rotmat(rmnew);//control matrix, changes only when the axis changes 
      //A=axismat(rnew);//moving matrix, controls rotation between control matrices
      g1=I(1,rnew,rmnew,g);
      if(opt_test==1 && readstep%4096==2081 && v==1){
	cout<<"g1=I(1,rnew,rmnew,g)= "<<g1<<endl;
	cout<<"Inew(1,rnew,rmnew,g)= "<<Inew(1,rnew,rmnew,g)<<endl;
        cout<<"rnew= "<<rnew<<endl;
        cout<<"rmnew= "<<rmnew<<endl;
        cout<<"g matrix "<<g<<endl;   
      }
      g2=I(2,rnew,rmnew,g);
      g3=I(3,rnew,rmnew,g);
      g4=g(4,1)+g(4,2)+g(4,3);
      g1motion[step-2]=g1;
      g2motion[step-2]=g2;
      g3motion[step-2]=g3;
      g4motion[step-2]=g4;
    }
    else {
      g1=g1motion[step-2];
      g2=g2motion[step-2];
      g3=g3motion[step-2];
      g4=g4motion[step-2];
    }
    double gg1=g1-grf1;
    double gg2=g2-grf2;
    double gg3=g3-grf3;
    double gg4=g4-grf4;
    if(opt_test==1 && readstep%4096==2081 && v==1){
      cout<<"g4="<<g4<<"; grf4="<<grf4<<"; gg4="<<gg4<<endl;
    }
    double tt=tnew-trf;
    if (told>=timecourse[actstep] && actstep<=(Nact-2)){
      coeff(activation[actstep],activation[actstep+1],timecourse[actstep],timecourse[actstep+1],dT2_1,dT2_2);
      dT2_1=dT2_1*iT2*iT2;
      dT2_2=dT2_2*iT2*iT2;
      actstep=actstep+1;
    }
    actint+=(dT2_1+dT2_2*(tnew+told)/2)*(tnew-told);
    if(opt_test==1 && readstep%4096==2081 && v==1){
      cout.precision(20);
      cout<<"actstep= "<<actstep<<"; actval= "<<activation[actstep]<<"; acttime="<<timecourse[actstep]<<"; dT2_1= "<<dT2_1<<"; dT2_2= "<<dT2_2<<"; iT2= "<<iT2<<"; actint= "<<actint<<endl;
    }
    double phase=gama*(gg1*x+gg2*y+gg3*z+gg4+(b0+chshift)*tt);
    if (fabs(rfangle)>1e-05){
      excitation=0;
      double df=H(step,3);//frequency width
      double fc=H(step,4);//center frequency
      rfstep=rfstep+1;
      if (v==1){
        //trnew1=H(step,9);//translation at tnew
        //trnew2=H(step,10);//
        //trnew3=H(step,11);//
        RowVector gradnew(3);
	gradnew <<gxnew<<gynew<<gznew;
        RowVector rr=gammabar*gradnew*(rotmat(rnew)*rotmat(rmnew));
	rr1=rr(1);rr2=rr(2);rr3=rr(3);
        trr=gammabar*(gxnew*trnew1+gynew*trnew2+gznew*trnew3);
	if (opt_test==1 && readstep%4096==0 && v==1) {
          cout<<"gxnew= "<<gxnew<<"; trnew1= "<<trnew1<<"; gammabar= "<<gammabar<<endl;
          cout<<"gynew= "<<gynew<<"; trnew2= "<<trnew2<<"; gammabar= "<<gammabar<<endl;
          cout<<"gznew= "<<gznew<<"; trnew3= "<<trnew3<<"; gammabar= "<<gammabar<<endl;
	  cout<<"trr=gammabar*(gxnew*trnew1+gynew*trnew2+gznew*trnew3)= "<<trr<<endl;
        }
        rotmotion1[rfstep-1]=rr1;
        rotmotion2[rfstep-1]=rr2;
        rotmotion3[rfstep-1]=rr3;
        transmotion[rfstep-1]=trr;
      }
      else {
        rr1=rotmotion1[rfstep-1];
        rr2=rotmotion2[rfstep-1];
        rr3=rotmotion3[rfstep-1];
        trr=transmotion[rfstep-1];
      }
      double f=rr1*x+rr2*y+rr3*z+trr+(b0+chshift)*gammabar;//gammabar is already included and grad included in trr
      double fval=(f-fc)/df;
      double off=(fval-dslcp_first)/dslcp;
      int nf=(int) off;
      //int nf=(int) (round_ivana(((f-fc)/df-dslcp_first)/dslcp,0));
       if (opt_test==1 && readstep%4096==0 && v==1 ){
        cout<<"fc_slc= "<<fc<<"; df_slc= "<<df<<endl;
        cout<<"rr1*x+rr2*y+rr3*z= "<<rr1*x+rr2*y+rr3*z<<endl;
        cout<<"trr= "<<trr<<endl;
        cout<<"(b0+chshift)*gammabar= "<<(b0+chshift)*gammabar<<endl;
	cout<<"fc_vox=rr1*x+rr2*y+rr3*z+trr+(b0+chshift)*gammabar= "<<f<<endl;
	cout<<"f in the table  nf= (int)((f-fc)*321/df+400+0.5)= "<<nf<<endl;
      }
       if (nf>=0 && nf<=(Nslc-2)) { 
        off-=nf;
        double ts=table_slcprof[nf];
        double sx=(table_slcprof[nf+1]-ts)*off + ts;
        double rfangle_f=sx*rfangle*RFtrans;//RFtrans are values 0 to 1 to derscribe the inhomogeneity f the receive RF field. 1 is for perfectly homog;
	if (opt_test==1 && readstep%4096==0 && v==1 ){
	  cout<<"table_slcprof[nf]= "<<sx<<"; rfangle= "<<rfangle<<"; RFtrans= "<<RFtrans<<endl;
	  cout<<"rfangle_f=table_slcprof[nf]*rfangle*RFtrans= "<<rfangle_f<<endl;
	}
        // if (fc-df/2<=f && f<=fc+df/2) {
        if (fabs(rfangle_f)>1e-06){//new stuff mon dec 19 //Edit: 29.11.12
	  excitation=1;
         if (opt_test==1 && readstep%4096==0 && readstep>7*4096){ 
          cout<<free_cout(m,tt,tissue,phase,actint)<<endl;
	 }
          m=free(m,tt,tissue,phase,actint);
	  //new stuff mon dec 19
	  //due to crushers or any gradient induced dephasing over the voxel a new initial magnetisation is introduced which is the average of the magnetisations over the voxel  
          double xvalrf=fabs(glo_cx*(gg1+b0x*tt));
          double yvalrf=fabs(glo_cy*(gg2+b0y*tt));
          double zvalrf=fabs(glo_cz*(gg3+b0z*tt));
	  double xyzrf=Sinc(xvalrf)*Sinc(yvalrf)*Sinc(zvalrf);
	  m(1)=m(1)*xyzrf;
          m(2)=m(2)*xyzrf;
          m=rot(rfangle_f,"x")*m;
          //new stuff
          m00=sqrt(m(1)*m(1)+m(2)*m(2));
          if (opt_test==1 && readstep%4096==0 && v==1){ 
            cout<<"Projection of the magnetisation vector into the xy plane after flipping is "<<m00<<endl;
	  }
          trf=tnew;
          grf1=g1;
          grf2=g2;
          grf3=g3;
	  grf4=g4;
	  actint=0.0;
	  //}
	}//new stuff
      }
    }
    if (read!=0){
      readstep=readstep+1;
      if (excitation==1 || nospeedup==1){
        double xval=fabs(glo_cx*(gg1+b0x*tt));
        double yval=fabs(glo_cy*(gg2+b0y*tt));
        double zval=fabs(glo_cz*(gg3+b0z*tt));
        if (v==1 && save_kcoord==1){
          coord(1,readstep)=gammabar*(gg1+b0x*tt);//to record the distorted coordinatres in the k-space
          coord(2,readstep)=gammabar*(gg2+b0y*tt);
          coord(3,readstep)=gammabar*(gg3+b0z*tt);
        }
	double tmp;
#ifdef NOTABLE
        tmp=m00*exp(-tt*iT2+actint)*Sinc(xval)*Sinc(yval)*Sinc(zval);
        sreal[readstep-1]+=den*tmp*cos(phase);
        simag[readstep-1]+=den*tmp*sin(phase);
#else
        if (xval<=glo_Dsinc && yval<=glo_Dsinc && zval<=glo_Dsinc){
          //TABLES for SINC
	  //x
          double off=xval*glo_idsinc;
	  int nbin=(int) off;
	  off-=nbin;
	  double ts=table_sinc[nbin];
	  double sx=(table_sinc[nbin+1]-ts)*off + ts;
	  //y
          off=yval*glo_idsinc;
	  nbin=(int) off;
	  off-=nbin;
	  ts=table_sinc[nbin];
	  double sy=(table_sinc[nbin+1]-ts)*off + ts;
	  //z        
          off=zval*glo_idsinc;
	  nbin=(int) off;
	  off-=nbin;
	  ts=table_sinc[nbin];
	  double sz=(table_sinc[nbin+1]-ts)*off + ts;
	  //
          tmp=m00*exp(-tt*iT2+actint)*sx*sy*sz;
	  if (opt_test==1 && readstep%4096==2081 && v==1){
	    cout.precision(20);
	    cout<<"table_sinc(xval)= "<<sx<<"; table_sinc(yval)= "<<sy<<"; table_sinc(zval)= "<<sz<<endl;
	    cout<<"sinc-table(xval)= "<<Sinc(xval)-sx<<"; sinc-table(yval)= "<<Sinc(yval)-sy<<"; sinc-tabel(zval)= "<<Sinc(zval)-sz<<endl;
	  }
        }
        else {
          tmp=m00*exp(-tt*iT2+actint)*Sinc(xval)*Sinc(yval)*Sinc(zval);
        }
        //TABLES SIN AND COS CALCULATION
        double phase_2pi;
        if (phase>0) phase_2pi=phase- ((int) (phase*glo_itwopi))*glo_twopi;//one solution when phase exceedes Dsin , this is faster
        else phase_2pi=phase- ((int) (phase*glo_itwopi))*glo_twopi+glo_twopi;
        double off=phase_2pi*glo_idsin;
        int nphase=(int) off;
        off-=nphase;
        double ts1=table_sin[nphase], tc1=table_cos[nphase];
        double wanted_cos=(table_cos[nphase+1]-tc1)*off+tc1;
        double wanted_sin=(table_sin[nphase+1]-ts1)*off+ts1;
        //END TABLES SIN AND COS CALCULATION
        sreal[readstep-1]+=den*tmp*wanted_cos;
        simag[readstep-1]+=den*tmp*wanted_sin;
        if(opt_test==1 && readstep%4096==2081 && v==1){
          cout.precision(20);
          cout<<"readstep= "<<readstep<<endl;
          cout<<"gama= "<<gama<<endl;
          cout<<"gg1= "<<gg1<<"; gg2= "<<gg2<<"; gg3= "<<gg3<<"; gg4="<<gg4<<endl;
          cout<<"x ="<<x<<"; y= "<<y<<"; z= "<<z<<endl;
          cout<<"b0= "<<b0<<"; chshift= "<<chshift<<endl;
          cout<<"tnew= "<<tnew<<"; trf= "<<trf<<endl;
          cout<<"tt=tnew-trf= "<<tt<<endl;
          cout<<"phase=gama*(gg1*x+gg2*y+gg3*z+gg4+(b0+chshift)*tt)= "<<phase<<endl; 
	  cout<<"phase_2pi= "<<phase_2pi<<endl;
	  cout<<"table_cos(phase)= "<<wanted_cos<<"; cos-table(phase)= "<<cos(phase)-wanted_cos<<endl;
          cout<<"table_sin(phase)= "<<wanted_sin<<"; sin-table(phase)= "<<sin(phase)-wanted_sin<<endl;
	  cout<<"glo_cx= "<<glo_cx<<"; glo_cy= "<<glo_cy<<"; glo_cz= "<<glo_cz<<endl;
	  cout<<"b0x= "<<b0x<<"; b0y= "<<b0y<<"; b0z= "<<b0z<<endl;
          cout<<"Ival=fabs(glo_cI*(ggI+b0I*tt)); I=x,y,z"<<endl;
          cout<<"xval= "<<xval<<"; yval= "<<yval<<"; zval= "<<zval<<endl;
          cout<<"m00= "<<m00<<endl;
          cout<<"actint= "<<actint<<"; iT2= "<<iT2<<" so exp(-tt*iT2+actint)= "<<exp(-tt*iT2+actint)<<endl;
	  cout<<"m00*exp(-tt*iT2+actint)="<<m00*exp(-tt*iT2+actint)<<endl;
          cout<<"Sinc(xval)*Sinc(yval)*Sinc(zval)= "<<Sinc(xval)*Sinc(yval)*Sinc(zval)<<endl;
          cout<<"tmp= m00*exp(-tt*iT2+actint)*Sinc(xval)*Sinc(yval)*Sinc(zval)= "<<tmp<<endl;
	  cout<<"den= "<<den<<endl;
	  cout<<"rfstep= "<<rfstep<<endl;
          cout<<"sreal_1voxel("<<readstep<<")= den*tmp*table_cos[nphase]= "<<den*tmp*wanted_cos<<"; alltillnow= "<<sreal[readstep-1]<<endl;
          cout<<"simag_1voxel("<<readstep<<")= den*tmp*table_sin[nphase]= "<<den*tmp*wanted_sin<<"; alltillnow= "<<simag[readstep-1]<<endl;
	  cout<<"-------------------------------------------------------------------------------"<<endl;
	}
#endif
      }         
    }//end of if read    
  }//end of main loop
 if (v==1 && save_kcoord==1){
    write_binary_matrix(coord,outputname+"_kcoord" );
 }
  }

//5.12.12///////////////////////additional option for RF angle averaging //////////////////////
void voxel3(const double x,const double y,const double z, 
	    const RowVector& tissue, const PMatrix& H, int const segA,
	    const int nrf,const int nreadp, const int nonzero, const int v,
	    const double xdim,const double ydim,const double zdim,
            const double b11,const double b12,const double b13,
	    const double b21,const double b22,const double b23,
	    const double b31,const double b32,const double b33, 
	    const double bx11,const double bx12,const double bx13,
	    const double bx21,const double bx22,const double bx23,
	    const double bx31,const double bx32,const double bx33,
	    const double by11,const double by12,const double by13,
	    const double by21,const double by22,const double by23,
	    const double by31,const double by32,const double by33,
            const double bz11,const double bz12,const double bz13,
	    const double bz21,const double bz22,const double bz23,
	    const double bz31,const double bz32,const double bz33,
            const double* timecourse, const double* activation,const int Nact,
	    const string outputname,const double* table_slcprof, const double dslcp, 
	    const double dslcp_first, const int Nslc,
            const double den, const double RFtrans,const int opt_test,
            const int nospeedup,const int save_kcoord, const bool rfavg,
            double* sreal, double* simag) {
  // Returns a signal (one rowvector for the real part and one for the imaginary part) 
  //for one voxel at place x,y,z of the phantom through time
  // tissue = (T_1,T_2,\rho) ; H = event matrix (pulse sequence + motion)
  // nrf=number of rf pulses ; v = internal voxel index number (used to check if this is 
  //the first voxel)
  // xdim, ydim and zdim are dimensions of the voxel
  // b1-9 are values of base of the perturbed field in the center of the voxel
  //(calculated from Maxwell's equations)  
  // bx1-9, by1-9 and bz1-9 are spatial gradients of the perturbed field 
  //(calculated with trilinear interpolation)
  // timecourse is a 2-column matrix where time is in first and multiplicative 
  //factor for the activation const "activation" in the second
  // sreal and simag are the outputs 

  //////////////////////////////////////////////////////////////////////////
  //READING IN THE SAVED STATE
  ///////////////////////////////////////////////////////////////////////// 
  int npst=16;
  static vector< vector<double> > pstate(nonzero, vector<double>(npst,0.0));
  if (segA==1){
    pstate[v-1][2]=tissue(3);
  }
  ColumnVector m(3);//magnetization vector: M
  m(1)=pstate[v-1][0];
  m(2)=pstate[v-1][1];
  m(3)=pstate[v-1][2];
  int excitation=(int)pstate[v-1][3];
  double trf=pstate[v-1][4];//last RF time
  double grf1=pstate[v-1][5];
  double grf2=pstate[v-1][6];
  double grf3=pstate[v-1][7];
  double grf4=pstate[v-1][8];
  double b0susrf=pstate[v-1][9];
  double b0xsusrf=pstate[v-1][10];
  double b0ysusrf=pstate[v-1][11];
  double b0zsusrf=pstate[v-1][12];
  double actint=pstate[v-1][13];
  int actstep=(int)pstate[v-1][14];
  double m00=pstate[v-1][15];
  //for V=1 
  int npste=13;
  static vector<double> pstatev(npste,0.0);
  static Matrix g(4,3); 
  int readstep=(int) pstatev[0];
  static int ars=0;
  static int arsf=0;
  static double b0tmp11=0.0; 
  static double b0tmp12=0.0;
  static double b0tmp13=0.0;
  static double b0tmp21=0.0;
  static double b0tmp22=0.0;
  static double b0tmp23=0.0;
  static double b0tmp31=0.0;
  static double b0tmp32=0.0;
  static double b0tmp33=0.0;
  static double b0tmp11freq=0.0;
  static double b0tmp12freq=0.0;
  static double b0tmp13freq=0.0;
  static double b0tmp21freq=0.0;
  static double b0tmp22freq=0.0;
  static double b0tmp23freq=0.0;
  static double b0tmp31freq=0.0;
  static double b0tmp32freq=0.0;
  static double b0tmp33freq=0.0;
  /////////////////////////////////////////////////////////////////////////
  //PHASE,GRADIENTS,COUNTERS,CONSTATNTS
  ////////////////////////////////////////////////////////////////////////
  double chshift=tissue(4);//chemical shift
  int numpoints=H.Nrows();//number of timepoints in the eventsequencer matrix
  double T2=tissue(2);
  double iT2=1/T2;
  double g1,g2,g3,g4;
  double rr1,rr2,rr3,trr;
  Matrix coord(3,nreadp);
  RowVector rnew(4),rmnew(4);//angle, axis values
  double trnew1=0.0;//translation
  double trnew2=0.0;
  double trnew3=0.0;
  //slice-profile weight def for the sections
  if (v==1) {
     g<<pstatev[1]<<pstatev[2]<<pstatev[3]
     <<pstatev[4]<<pstatev[5]<<pstatev[6]
     <<pstatev[7]<<pstatev[8]<<pstatev[9]
     <<pstatev[10]<<pstatev[11]<<pstatev[12];
     ars=b0tmp111.size();
     arsf=b0tmp111freq.size();
     b0tmp11=b0tmp111[ars-1]; 
     b0tmp12=b0tmp112[ars-1]; 
     b0tmp13=b0tmp113[ars-1];
     b0tmp21=b0tmp221[ars-1]; 
     b0tmp22=b0tmp222[ars-1]; 
     b0tmp23=b0tmp223[ars-1]; 
     b0tmp31=b0tmp331[ars-1]; 
     b0tmp32=b0tmp332[ars-1]; 
     b0tmp33=b0tmp333[ars-1]; 
     b0tmp11freq=b0tmp111freq[arsf-1]; 
     b0tmp12freq=b0tmp112freq[arsf-1]; 
     b0tmp13freq=b0tmp113freq[arsf-1];
     b0tmp21freq=b0tmp221freq[arsf-1]; 
     b0tmp22freq=b0tmp222freq[arsf-1]; 
     b0tmp23freq=b0tmp223freq[arsf-1]; 
     b0tmp31freq=b0tmp331freq[arsf-1]; 
     b0tmp32freq=b0tmp332freq[arsf-1]; 
     b0tmp33freq=b0tmp333freq[arsf-1]; 
     g1motion.resize(numpoints-1,0.0);
     g2motion.resize(numpoints-1,0.0);
     g3motion.resize(numpoints-1,0.0);
     g4motion.resize(numpoints-1,0.0);
     rotmotion1.resize(nrf,0.0);
     rotmotion2.resize(nrf,0.0);
     rotmotion3.resize(nrf,0.0);
     transmotion.resize(nrf,0.0);
     b0tmp111.resize(numpoints-1,0.0);
     b0tmp112.resize(numpoints-1,0.0);
     b0tmp113.resize(numpoints-1,0.0);
     b0tmp221.resize(numpoints-1,0.0);
     b0tmp222.resize(numpoints-1,0.0);
     b0tmp223.resize(numpoints-1,0.0);
     b0tmp331.resize(numpoints-1,0.0);
     b0tmp332.resize(numpoints-1,0.0);
     b0tmp333.resize(numpoints-1,0.0);
     b0tmp111freq.resize(nrf,0.0);
     b0tmp112freq.resize(nrf,0.0);
     b0tmp113freq.resize(nrf,0.0);
     b0tmp221freq.resize(nrf,0.0);
     b0tmp222freq.resize(nrf,0.0);
     b0tmp223freq.resize(nrf,0.0);
     b0tmp331freq.resize(nrf,0.0);
     b0tmp332freq.resize(nrf,0.0);
     b0tmp333freq.resize(nrf,0.0);
     ////////////////////////////////////////////////////////////////////////
     //LOOK UP TABLES
     ///////////////////////////////////////////////////////////////////////
     if (table_sinc==0) {
	 table_sinc = new double[glo_Nsinc];
	 for (int n=0;n<=glo_Nsinc;n++){
	   table_sinc[n]=Sinc(n*glo_dsinc);
	 }
     }
     if ((table_sin==0) || (table_cos==0)) {
       table_sin = new double[glo_Nsin];
       table_cos = new double[glo_Nsin];
       for (int n=0;n<=glo_Nsin;n++){
	 table_sin[n]=sin(n*glo_dsin-glo_Dsin);
	 table_cos[n]=cos(n*glo_dsin-glo_Dsin);
       }
     }
     if (opt_test==1){
       cout<<"Stepsize for Table for SINC: dsinc= "<<glo_dsinc<<endl;
       cout<<"Stepsize for Table for SIN: dsin="<<glo_dsin<<endl;
     }
     glo_cx=gammabar*xdim; //(Hz*m/T)
     glo_cy=gammabar*ydim;
     glo_cz=gammabar*zdim;
  }
  //we assume that for the step=1 all the values in the mainmatrix H are zero
  /////////////////////////////////////////////////////////////////////////
  // MOTION PARAMETERS, MAGNETIC SUSCEPTIBILITY PARAMETERS
  /////////////////////////////////////////////////////////////////////////
  double tt=0;
  double b0sus=0.0;
  double b0xsus=0.0;
  double b0ysus=0.0;
  double b0zsus=0.0;
  double b00=0.0;
  double b0x=0.0;
  double b0y=0.0;
  double b0z=0.0;
  double fb0=0.0;
  int rfstep=0;
  RowVector b0tmp1(3);//[0 0 1]*R
  RowVector b0tmp2(3);//[0 0 1]*A*R
  RowVector b0tmp3(3);//[0 0 1]*A*A*R
  //////////////////////////////////////////////////////////////////////////
  //  ACTIVATION PARAMETERS
  //////////////////////////////////////////////////////////////////////////
  double dT2_1=0.0;
  double dT2_2=0.0;
  //////////////////////////////////////////////////////////////////////////
  //MAIN LOOP
  for (int step=2;step<=numpoints;step++){
    double told=H.time(step-1);//time at the previous step
    double tnew=H.time(step);//time at this step
    double rfangle=H(step,2);//flip angle
    double aold=H(step-1,12);//angle from last given point to the 
    //interpolated one (they have the same axis of rotation) at told
    double anew=H(step,12);//angle from last given point to the 
    //interpolated one (they have the same axis of rotation) at tnew
    double read=H(step,5);
    if (v==1){
      rnew<<H(step,12)<<H(step,13)<<H(step,14)<<H(step,15);//rotation
      // betwen given motion points at tnew (angle & axis)
      rmnew<<H(step,16)<<H(step,17)<<H(step,18)<<H(step,19);//rotation at 
      //control motion point at tnew (angle & axis)
      double gxold=H(step-1,6);//gradients at told 
      double gyold=H(step-1,7);//
      double gzold=H(step-1,8);//
      double gxnew=H(step,6);// gradients at tnew
      double gynew=H(step,7);//
      double gznew=H(step,8);//
      double trold1=H(step-1,9);//translation at told
      double trold2=H(step-1,10);//
      double trold3=H(step-1,11);//
      double trnew1=H(step,9);//translation at tnew
      double trnew2=H(step,10);//
      double trnew3=H(step,11);//
      if (fabs(gxold)>1e-16 || fabs(gxnew)>1e-16){
        g(1,1)+=i1(gxold,gxnew,told, tnew);
        g(2,1)+=i2(gxold,gxnew,aold,anew,told,tnew);
        g(3,1)+=i1(gxold,gxnew,told, tnew)- i3(gxold,gxnew,aold,anew,told,tnew);
        g(4,1)+=i4(gxold,gxnew,trold1,trnew1,told,tnew);
      }
      if (fabs(gyold)>1e-16 || fabs(gynew)>1e-16){
        g(1,2)+=i1(gyold,gynew,told,tnew);
        g(2,2)+=i2(gyold,gynew,aold,anew,told,tnew);
        g(3,2)+=i1(gyold,gynew,told,tnew)- i3(gyold,gynew,aold,anew,told,tnew);
        g(4,2)+=i4(gyold,gynew,trold2,trnew2,told,tnew);
      }
      if (fabs(gzold)>1e-16 || fabs(gznew)>1e-16){
        g(1,3)+=i1(gzold,gznew,told,tnew);
        g(2,3)+=i2(gzold,gznew,aold,anew,told,tnew);
        g(3,3)+=i1(gzold,gznew,told,tnew)- i3(gzold,gznew,aold,anew,told,tnew);
        g(4,3)+=i4(gzold,gznew,trold3,trnew3,told,tnew);
      }
      Matrix R(3,3);
      Matrix A(3,3);
      R=rotmat(rmnew);//control matrix, changes only when the axis changes 
      A=axismat(rnew);//moving matrix, controls rotation between control matrices
      RowVector b0(3);
      b0 <<(double)0<<(double)0<<(double)1;//static magnetic field Bo
      b0tmp1=b0*R;
      b0tmp2=b0*A*R;
      b0tmp3=b0*A*A*R;
      b0tmp11+=b0int(b0tmp1(1),b0tmp1(1),b0tmp2(1),b0tmp2(1),
		     b0tmp3(1),b0tmp3(1),aold,anew,told,tnew);
      b0tmp12+=b0int(b0tmp1(1),b0tmp1(2),b0tmp2(1),b0tmp2(2),
		     b0tmp3(1),b0tmp3(2),aold,anew,told,tnew);
      b0tmp13+=b0int(b0tmp1(1),b0tmp1(3),b0tmp2(1),b0tmp2(3),
		     b0tmp3(1),b0tmp3(3),aold,anew,told,tnew);
      b0tmp21+=b0int(b0tmp1(2),b0tmp1(1),b0tmp2(2),b0tmp2(1),
		     b0tmp3(2),b0tmp3(1),aold,anew,told,tnew);
      b0tmp22+=b0int(b0tmp1(2),b0tmp1(2),b0tmp2(2),b0tmp2(2),
		     b0tmp3(2),b0tmp3(2),aold,anew,told,tnew);
      b0tmp23+=b0int(b0tmp1(2),b0tmp1(3),b0tmp2(2),b0tmp2(3),
		     b0tmp3(2),b0tmp3(3),aold,anew,told,tnew);
      b0tmp31+=b0int(b0tmp1(3),b0tmp1(1),b0tmp2(3),b0tmp2(1),
		     b0tmp3(3),b0tmp3(1),aold,anew,told,tnew);
      b0tmp32+=b0int(b0tmp1(3),b0tmp1(2),b0tmp2(3),b0tmp2(2),
		     b0tmp3(3),b0tmp3(2),aold,anew,told,tnew);
      b0tmp33+=b0int(b0tmp1(3),b0tmp1(3),b0tmp2(3),b0tmp2(3),
		     b0tmp3(3),b0tmp3(3),aold,anew,told,tnew);
      assert ( (step-2) < numpoints);
      b0tmp111[step-2]=b0tmp11;
      b0tmp112[step-2]=b0tmp12;
      b0tmp113[step-2]=b0tmp13;
      b0tmp221[step-2]=b0tmp21;
      b0tmp222[step-2]=b0tmp22;
      b0tmp223[step-2]=b0tmp23;
      b0tmp331[step-2]=b0tmp31;
      b0tmp332[step-2]=b0tmp32;
      b0tmp333[step-2]=b0tmp33;
      //b0 ends
      g1=I(1,rnew,rmnew,g);
      if(opt_test==1 && readstep%4096==2080 && v==1){
	cout<<"Inew(1,rnew,rmnew,g)= "<<Inew(1,rnew,rmnew,g)<<endl;
        cout<<"rnew= "<<rnew<<endl;
        cout<<"rmnew= "<<rmnew<<endl;
        cout<<"g matrix "<<g<<endl;
        cout<<"g INTEGRALS:"<<endl;
	cout<<"i1(gxold,gxnew,told, tnew)="<<i1new(gxold,gxnew,told, tnew)<<endl;
        cout<<"i2(gxold,gxnew,aold,anew,told,tnew)="<<i2new(gxold,gxnew,aold,
							    anew,told,tnew)<<endl;
        cout<<"i3(gxold,gxnew,aold,anew,told,tnew)="<<i3new(gxold,gxnew,aold,
							    anew,told,tnew)<<endl;
        cout<<"i4(gxold,gxnew,aold,anew,told,tnew)="<<i4new(gxold,gxnew,aold,
							    anew,told,tnew)<<endl;
      }
      g2=I(2,rnew,rmnew,g);
      g3=I(3,rnew,rmnew,g);
      g4=g(4,1)+g(4,2)+g(4,3);
      g1motion[step-2]=g1;
      g2motion[step-2]=g2;
      g3motion[step-2]=g3;
      g4motion[step-2]=g4;
    }
    else {
      g1=g1motion[step-2];
      g2=g2motion[step-2];
      g3=g3motion[step-2];
      g4=g4motion[step-2];
      //b0 starts
      b0tmp11=b0tmp111[step-2];
      b0tmp12=b0tmp112[step-2];
      b0tmp13=b0tmp113[step-2];
      b0tmp21=b0tmp221[step-2];
      b0tmp22=b0tmp222[step-2];
      b0tmp23=b0tmp223[step-2];
      b0tmp31=b0tmp331[step-2];
      b0tmp32=b0tmp332[step-2];
      b0tmp33=b0tmp333[step-2];
    }
    double gg1=g1-grf1;
    double gg2=g2-grf2;
    double gg3=g3-grf3;
    double gg4=g4-grf4;
    if(opt_test==1 && readstep%4096==2080 && v==1){
      cout<<"g1="<<g1<<"; grf1="<<grf1<<"; gg1="<<gg1<<endl;
      cout<<"g2="<<g2<<"; grf2="<<grf2<<"; gg2="<<gg2<<endl;
      cout<<"g3="<<g3<<"; grf3="<<grf3<<"; gg3="<<gg3<<endl;
      cout<<"g4="<<g4<<"; grf4="<<grf4<<"; gg4="<<gg4<<endl;
    }
    tt=tnew-trf;
    if (told>=timecourse[actstep] && actstep<=(Nact-2)){
	coeff(activation[actstep],activation[actstep+1],timecourse[actstep],
	      timecourse[actstep+1],dT2_1,dT2_2);
	dT2_1=dT2_1*iT2*iT2;
	dT2_2=dT2_2*iT2*iT2;
        actstep=actstep+1;
    }
    actint+=(dT2_1+dT2_2*(tnew+told)/2)*(tnew-told);
    if(opt_test==1 && readstep%4096==2080 && v==1){
      cout<<"actstep= "<<actstep<<"; actval= "<<activation[actstep]<<endl;
      cout<<"; acttime="<<timecourse[actstep]<<"; dT2_1= "<<dT2_1<<endl;
      cout<<"; dT2_2= "<<dT2_2<<"; iT2= "<<iT2<<"; actint= "<<actint<<endl;
    }
    //b0 starts
    b0sus= b0tmp11*b11+b0tmp12*b12+b0tmp13*b13+b0tmp21*b21+b0tmp22*b22+
      b0tmp23*b23+b0tmp31*b31+b0tmp32*b32+b0tmp33*b33;
    b0xsus= b0tmp11*bx11+b0tmp12*bx12+b0tmp13*bx13+b0tmp21*bx21+b0tmp22*
      bx22+b0tmp23*bx23+b0tmp31*bx31+b0tmp32*bx32+b0tmp33*bx33;
    b0ysus= b0tmp11*by11+b0tmp12*by12+b0tmp13*by13+b0tmp21*by21+b0tmp22*
      by22+b0tmp23*by23+b0tmp31*by31+b0tmp32*by32+b0tmp33*by33;
    b0zsus= b0tmp11*bz11+b0tmp12*bz12+b0tmp13*bz13+b0tmp21*bz21+b0tmp22*
      bz22+b0tmp23*bz23+b0tmp31*bz31+b0tmp32*bz32+b0tmp33*bz33;
    b00=b0sus-b0susrf;
    b0x=b0xsus-b0xsusrf;
    b0y=b0ysus-b0ysusrf;
    b0z=b0zsus-b0zsusrf;
    //b0 ends
    double phase=gama*(gg1*x+gg2*y+gg3*z+gg4+b00+chshift*tt);
    if (fabs(rfangle)>1e-06){
      excitation=0;
      double df=H(step,3);//frequency width
      double fc=H(step,4);//center frequency
      rfstep=rfstep+1;
      if (v==1){
        RowVector gradnew(3);
	gradnew<<H(step,6)<<H(step,7)<<H(step,8);
        RowVector rr=gammabar*gradnew*(rotmat(rnew)*rotmat(rmnew));
	rr1=rr(1);rr2=rr(2);rr3=rr(3);
        trr=gammabar*(gradnew(1)*trnew1+gradnew(2)*trnew2+gradnew(3)*trnew3);
	assert( (rfstep - 1) < nrf );
        rotmotion1[rfstep-1]=rr1;
        rotmotion2[rfstep-1]=rr2;
        rotmotion3[rfstep-1]=rr3;
        transmotion[rfstep-1]=trr;
        b0tmp11freq=b0freq(b0tmp1(1),b0tmp1(1),b0tmp2(1),b0tmp2(1),b0tmp3(1),b0tmp3(1),anew);
        b0tmp12freq=b0freq(b0tmp1(1),b0tmp1(2),b0tmp2(1),b0tmp2(2),b0tmp3(1),b0tmp3(2),anew);
        b0tmp13freq=b0freq(b0tmp1(1),b0tmp1(3),b0tmp2(1),b0tmp2(3),b0tmp3(1),b0tmp3(3),anew);
        b0tmp21freq=b0freq(b0tmp1(2),b0tmp1(1),b0tmp2(2),b0tmp2(1),b0tmp3(2),b0tmp3(1),anew);
        b0tmp22freq=b0freq(b0tmp1(2),b0tmp1(2),b0tmp2(2),b0tmp2(2),b0tmp3(2),b0tmp3(2),anew);
        b0tmp23freq=b0freq(b0tmp1(2),b0tmp1(3),b0tmp2(2),b0tmp2(3),b0tmp3(2),b0tmp3(3),anew);
        b0tmp31freq=b0freq(b0tmp1(3),b0tmp1(1),b0tmp2(3),b0tmp2(1),b0tmp3(3),b0tmp3(1),anew);
        b0tmp32freq=b0freq(b0tmp1(3),b0tmp1(2),b0tmp2(3),b0tmp2(2),b0tmp3(3),b0tmp3(2),anew);
        b0tmp33freq=b0freq(b0tmp1(3),b0tmp1(3),b0tmp2(3),b0tmp2(3),b0tmp3(3),b0tmp3(3),anew);
        b0tmp111freq[rfstep-1]=b0tmp11freq;
        b0tmp112freq[rfstep-1]=b0tmp12freq;
        b0tmp113freq[rfstep-1]=b0tmp13freq;
        b0tmp221freq[rfstep-1]=b0tmp21freq;
        b0tmp222freq[rfstep-1]=b0tmp22freq;
        b0tmp223freq[rfstep-1]=b0tmp23freq;
        b0tmp331freq[rfstep-1]=b0tmp31freq;
        b0tmp332freq[rfstep-1]=b0tmp32freq;
        b0tmp333freq[rfstep-1]=b0tmp33freq;
      }
      else {
        rr1=rotmotion1[rfstep-1];
        rr2=rotmotion2[rfstep-1];
        rr3=rotmotion3[rfstep-1];
        trr=transmotion[rfstep-1];
        b0tmp11freq=b0tmp111freq[rfstep-1];
        b0tmp12freq=b0tmp112freq[rfstep-1];
        b0tmp13freq=b0tmp113freq[rfstep-1];
        b0tmp21freq=b0tmp221freq[rfstep-1];
        b0tmp22freq=b0tmp222freq[rfstep-1];
        b0tmp23freq=b0tmp223freq[rfstep-1];
        b0tmp31freq=b0tmp331freq[rfstep-1];
        b0tmp32freq=b0tmp332freq[rfstep-1];
        b0tmp33freq=b0tmp333freq[rfstep-1];
      }
      fb0=gammabar*(b0tmp11freq*b11+b0tmp12freq*b12+b0tmp13freq*b13+
		    b0tmp21freq*b21+b0tmp22freq*b22+b0tmp23freq*b23+
		    b0tmp31freq*b31+b0tmp32freq*b32+b0tmp33freq*b33+chshift);
      double f=rr1*x+rr2*y+rr3*z+trr+fb0;//frequency at the centre of the voxel
      //finding table index for the center, 
      double fval=(f-fc)/df;
      double off=(fval-dslcp_first)/dslcp;
      int nf=(int) off;
       if (nf>=0 && nf<=(Nslc-2)) { 
        off-=nf;
    //29.11.12: added option for either slcprof signal weighting or rfangle weighting
	double prfweight = 1;
	if(!rfavg){
		double ts=table_slcprof[nf];
        double sx=(table_slcprof[nf+1]-ts)*off + ts;
        prfweight=sx;
	}
	else{
        double dfvox = rr1*(xdim/2)+rr2*(ydim/2)+rr3*(zdim/2);//we are assumming that motion
        //parameters are the same everywhere in the voxel
		//table index for the start of the voxel
        double fstart = (f-dfvox-fc)/df;
        double offstart=(fstart-dslcp_first)/dslcp;
        int nfstart=(int) offstart;
        //table index for the end of the voxel
        double fend = (f+dfvox-fc)/df;
        double offend=(fend-dslcp_first)/dslcp;
        int nfend=(int) offend;
        offstart-=nfstart;
		offend-=nfend;
        int prfcount = 0;
        for(int i=nfstart;i<=nfend;i++){
			prfweight += table_slcprof[i];
			prfcount++;
		}
		prfweight = prfweight/prfcount;
	}
        //cout<<"prfweight = "<<prfweight<<endl;
        double rfangle_f=prfweight*rfangle*RFtrans;//RFtrans are values 0 to 1 to derscribe 
        if (fabs(rfangle_f)>1e-06){//new stuff mon dec 19
          excitation=1;
          m=free(m,tt,tissue,phase,actint);
	  //new stuff mon dec 19
	  //due to crushers or any gradient induced dephasing over the voxel a new 
	  //initial magnetisation is introduced which is the average of the 
	  //magnetisations over the voxel  
          double xvalrf=fabs(glo_cx*(gg1+b0x*tt));
          double yvalrf=fabs(glo_cy*(gg2+b0y*tt));
          double zvalrf=fabs(glo_cz*(gg3+b0z*tt));
	  double xyzrf=Sinc(xvalrf)*Sinc(yvalrf)*Sinc(zvalrf);
	  m(1)=m(1)*xyzrf;
          m(2)=m(2)*xyzrf;
          m=rot(rfangle_f,"x")*m;
          //new stuff
          m00=sqrt(m(1)*m(1)+m(2)*m(2));
          trf=tnew;
          grf1=g1;
          grf2=g2;
          grf3=g3;
	  grf4=g4;
          b0susrf=b0sus;
	  b0xsusrf=b0xsus;
	  b0ysusrf=b0ysus;
	  b0zsusrf=b0zsus;
	  actint=0.0;
	  //}
	}//new stuff
      }
    }
    if (read!=0){
      readstep=readstep+1;
      if (excitation==1 || nospeedup==1){
        double xval=fabs(glo_cx*(gg1+b0x));
        double yval=fabs(glo_cy*(gg2+b0y));
        double zval=fabs(glo_cz*(gg3+b0z));
        if (v==1 && save_kcoord==1){
          coord(1,readstep)=gammabar*(gg1+b0x);//to record the distorted coordinatres 
	  //in the k-space
          coord(2,readstep)=gammabar*(gg2+b0y);
          coord(3,readstep)=gammabar*(gg3+b0z);
        }
	double tmp;
#ifdef NOTABLE
        tmp=m00*exp(-tt*iT2+actint)*Sinc(xval)*Sinc(yval)*Sinc(zval);
        sreal[readstep-1]+=den*tmp*cos(phase);
        simag[readstep-1]+=den*tmp*sin(phase);
#else
        if (xval<glo_Dsinc && yval<glo_Dsinc && zval<glo_Dsinc){
        //TABLES for SINC
	//x
        double off=xval*glo_idsinc;
	int nbin=(int) off;
	off-=nbin;
	assert ((nbin)<glo_Nsinc);
	double ts=table_sinc[nbin];
	double sx=(table_sinc[nbin+1]-ts)*off + ts;
	//y
        off=yval*glo_idsinc;
	nbin=(int) off;
	off-=nbin;
	ts=table_sinc[nbin];
	double sy=(table_sinc[nbin+1]-ts)*off + ts;
	//z        
        off=zval*glo_idsinc;
	nbin=(int) off;
	off-=nbin;
	ts=table_sinc[nbin];
	double sz=(table_sinc[nbin+1]-ts)*off + ts;
	//
        tmp=m00*exp(-tt*iT2+actint)*sx*sy*sz;

      }
      else {
        tmp=m00*exp(-tt*iT2+actint)*Sinc(xval)*Sinc(yval)*Sinc(zval);
      }
      //TABLES SIN AND COS CALCULATION
      double phase_2pi;
      if (phase>0) phase_2pi=phase- ((int) (phase*glo_itwopi))*glo_twopi;//one solution 
      //when phase exceedes Dsin , this is faster
      else phase_2pi=phase- ((int) (phase*glo_itwopi))*glo_twopi+glo_twopi;
      double off=phase_2pi*glo_idsin;
      int nphase=(int) off;
      off-=nphase;
      if (opt_test==1 && readstep%4096==2081 && v==1) {
        cout<<"nphase= "<<nphase<<"; Phase= "<<phase<<endl;
      	cout<<"((int) (phase*glo_itwopi))*glo_twopi"<<endl;
        cout<<((int) (phase*glo_itwopi))*glo_twopi<<" glo_twopi="<<glo_twopi<<endl;
        cout<<"; glo_itwopi="<<glo_itwopi<<"table_sin[nphase]"<<table_sin[nphase]<<endl;
      }
      assert ((nphase+1)<glo_Nsin+1);
      double ts1=table_sin[nphase], tc1=table_cos[nphase];
      double wanted_cos=(table_cos[nphase+1]-tc1)*off+tc1;
      double wanted_sin=(table_sin[nphase+1]-ts1)*off+ts1;
      //CALCULATING SIGNAL
      assert( (readstep-1) < nreadp );
      sreal[readstep-1]+=den*tmp*wanted_cos;
      simag[readstep-1]+=den*tmp*wanted_sin;
      if(opt_test==1 && readstep%4096==2081 && v==1){
          cout.precision(20);
          cout<<"readstep= "<<readstep<<endl;
          cout<<"gama= "<<gama<<endl;
          cout<<"gg1= "<<gg1<<"; gg2= "<<gg2<<"; gg3= "<<gg3<<"; gg4="<<gg4<<endl;
          cout<<"x ="<<x<<"; y= "<<y<<"; z= "<<z<<endl;
          cout<<"b0sus="<<b0sus<<"; b0susrf="<<b0susrf<<"; b0tmp11="<<b0tmp11<<endl;
	  cout<<"b11"<<b11<<endl;
          cout<<"b0= "<<b00<<"; chshift= "<<chshift<<endl;
          cout<<"tnew= "<<tnew<<"; trf= "<<trf<<endl;
          cout<<"tt=tnew-trf= "<<tt<<endl;
          cout<<"phase=gama*(gg1*x+gg2*y+gg3*z+gg4+(b0+chshift)*tt)= "<<phase<<endl; 
	  cout<<"phase_2pi= "<<phase_2pi<<endl;
	  cout<<"table_cos(phase)= "<<wanted_cos<<endl;
	  cout<<"cos-table(phase)= "<<cos(phase)-wanted_cos<<endl;
          cout<<"table_sin(phase)= "<<wanted_sin<<endl;
	  cout<<"; sin-table(phase)= "<<sin(phase)-wanted_sin<<endl;
	  cout<<"glo_cx= "<<glo_cx<<"; glo_cy= "<<glo_cy<<endl;
	  cout<<"; glo_cz= "<<glo_cz<<endl;
	  cout<<"b0x= "<<b0x<<"; b0y= "<<b0y<<"; b0z= "<<b0z<<endl;
          cout<<"Ival=fabs(glo_cI*(ggI+b0I*tt)); I=x,y,z"<<endl;
          cout<<"xval= "<<xval<<"; yval= "<<yval<<"; zval= "<<zval<<endl;
          cout<<"m00= "<<m00<<endl;
          cout<<"actint= "<<actint<<"; iT2= "<<iT2<<endl;
	  cout<<" so exp(-tt*iT2+actint)= "<<exp(-tt*iT2+actint)<<endl;
	  cout<<"m00*exp(-tt*iT2+actint)="<<m00*exp(-tt*iT2+actint)<<endl;
          cout<<"Sinc(xval)*Sinc(yval)*Sinc(zval)= "<<Sinc(xval)*Sinc(yval)*Sinc(zval)<<endl;
          cout<<"tmp= m00*exp(-tt*iT2+actint)*Sinc(xval)*Sinc(yval)*Sinc(zval)= "<<tmp<<endl;
	  cout<<"den= "<<den<<endl;
	  cout<<"rfstep= "<<rfstep<<endl;
          cout<<"sreal1vox("<<readstep<<")= den*tmp*wanted_cos="<<den*tmp*wanted_cos<<endl;
	  cout<<"alltillnow= "<<sreal[readstep-1]<<endl;
          cout<<"simag1vox("<<readstep<<")= den*tmp*wanted_sin"<<den*tmp*wanted_sin<<endl;
	  cout<<"alltillnow= "<<simag[readstep-1]<<endl;
	  cout<<"------------------------------------------------------------"<<endl;
      }
#endif
      }//end of if slcnum
    }//end of if read
  }//end of main loop
 //////////////////////////////////////////////////////////////////////////
  //SAVING THE STATE
  ///////////////////////////////////////////////////////////////////////// 
  pstate[v-1][0]=m(1);
  pstate[v-1][1]=m(2);
  pstate[v-1][2]=m(3);
  pstate[v-1][3]=(double)excitation;
  pstate[v-1][4]=trf;//last RF time
  pstate[v-1][5]=grf1;
  pstate[v-1][6]=grf2;
  pstate[v-1][7]=grf3;
  pstate[v-1][8]=grf4;
  pstate[v-1][9]=b0susrf;
  pstate[v-1][10]=b0xsusrf;
  pstate[v-1][11]=b0ysusrf;
  pstate[v-1][12]=b0zsusrf;
  pstate[v-1][13]=actint;
  pstate[v-1][14]=(double)actstep;
  pstate[v-1][15]=m00;
  //FOR V=1
  if (v==nonzero){
    pstatev[0]=(int)readstep;
    pstatev[1]=g(1,1);
    pstatev[2]=g(1,2);
    pstatev[3]=g(1,3);
    pstatev[4]=g(2,1);
    pstatev[5]=g(2,2);
    pstatev[6]=g(2,3);
    pstatev[7]=g(3,1);
    pstatev[8]=g(3,2);
    pstatev[9]=g(3,3);
    pstatev[10]=g(4,1);
    pstatev[11]=g(4,2);
    pstatev[12]=g(4,3);
  }
 if (v==1 && save_kcoord==1){
    write_binary_matrix(coord,outputname+"_kcoord" );
 }
}

////////////////////////////////
void voxel4(const double x,const double y,double z, 
            const RowVector& tissue,const PMatrix& H,
            const int nreadp,const int v,
            const double xdim,const double ydim,const double zdim,
            const double* b0time, const double* b0xtime,const double* b0ytime,const double* b0ztime,
            const double* b0timecourse, const int Nb0,
            const double b0, const double b0x,const double b0y,const double b0z,
            const double* timecourse, const double* activation,const int Nact,
	    const string outputname,  
	    const double* table_slcprof, const double dslcp, const double dslcp_first, const int Nslc,
            const double den,const double RFtrans,const int opt_test,
            const int nospeedup,
            const int save_kcoord,
            double* sreal, double* simag){
  // Returns a signal (one rowvector for the real part and one for the 
  //  imaginary part) for one voxel at place x,y,z(m) of the phantom through time
  // tissue = (T_1(s),T_2(s),\rho) ; H = event matrix (pulse sequence)
  // nreadp = number of readoutpoints !!! pay attention as you are not using this at all here;  
  // v = internal voxel index number (used to check if it is the first voxel)
  // idsx, idsy, idsz are constants = cx/ds where cx=gammabar dim_x (Hz*m/T) and ds=1/100000 interval in sinc table
  //b0shift is inhomogeneity in the field that can be due to either b0sus or chemical shift, it is in T
  //b0x,b0y,b0z are values of gradients of b0 inhomogeneity field for each voxel (T/m) 
  //table_sinc,sin,cos are there for use instead of conventional sinc, sin, cos --speed up
  //idss=1/dss, dss=1/40000 is an interval in the sin and cos table ; idsshelp = idss*2*pi;

  ///////////////////////////////////////////////////////////////////////////
  //INITIALIZATION OF THE MAGNETIZATION VECTOR 
  ///////////////////////////////////////////////////////////////////////////
  ColumnVector m(3);//magnetization vector
  m(1)=0;
  m(2)=0;
  m(3)=tissue(3);
  double m00=0;//the magnitude of the transverse magnetization vector
  /////////////////////////////////////////////////////////////////////////// 
  //PHASE,GRADIENTS,COUNTERS,CONSTANTS
  ///////////////////////////////////////////////////////////////////////////
  double chshift=tissue(4);//chemical shift
  int numpoints=H.Nrows();//number of timepoints in the eventsequencer matrix
  int readstep=0;//keeps track of readout points
  int excitation=0;
  double trf=0;//rftime
  double T2=tissue(2);
  double iT2=1/T2;
  int rftest=0;
  double g1,g2,g3,grf1,grf2,grf3; 
  g1=0.0; g2=0.0; g3=0.0; grf1=0.0; grf2=0.0; grf3=0.0;
  Matrix coord(3,nreadp); 
  if (v==1) {
     g1static=new double[numpoints];
     g2static=new double[numpoints];
     g3static=new double[numpoints];
     ////////////////////////////////////////////////////////////////////////
     //LOOK UP TABLES
     ////////////////////////////////////////////////////////////////////////
     table_sinc = new double[glo_Nsinc+1];
     table_sin = new double[glo_Nsin+1];
     table_cos = new double[glo_Nsin+1];
     if (opt_test==1){
       cout<<"Stepsize for Table for SINC: dsinc= "<<glo_dsinc<<endl;
       cout<<"Stepsize for Table for SIN: dsin="<<glo_dsin<<endl;
     }
     for (int n=0;n<=glo_Nsinc;n++){
       table_sinc[n]=Sinc(n*glo_dsinc);
     }
     for (int n=0;n<=glo_Nsin;n++){
       table_sin[n]=sin(n*glo_dsin);
       table_cos[n]=cos(n*glo_dsin);
     }
     glo_cx=gammabar*xdim; //(Hz*m/T)
     glo_cy=gammabar*ydim;
     glo_cz=gammabar*zdim;
     ////////////////////////////////////////////////////////////////////////
  } 
  ///////////////////////////////////////////////////////////////////////////
  //  ACTIVATION PARAMETERS
  ///////////////////////////////////////////////////////////////////////////
  double actint=0.0;
  int actstep=0;
  double dT2_1=0.0;
  double dT2_2=0.0;
  ///////////////////////////////////////////////////////////////////////////
  // b0 PARAMETERS (in the case b0 changes in TIME )
  ///////////////////////////////////////////////////////////////////////////
  double b0int=0.0;
  int b0step=0;
  double db0_1=0.0;
  double db0_2=0.0;
  double b0xint=0.0;
  double db0x_1=0.0;
  double db0x_2=0.0;
  double b0yint=0.0;
  double db0y_1=0.0;
  double db0y_2=0.0;
  double b0zint=0.0;
  double db0z_1=0.0;
  double db0z_2=0.0;
  double b0val=0.0;
  ///////////////////////////////////////////////////////////////////////////
  //SIGNAL
  /////////////////////////////////////////////////////////////////////////// 
  for (int step=2;step<=numpoints;step++){
    double tnew=H.time(step);
    double told=H.time(step-1);
    double rfangle=H(step,2);
    double read=H(step,5);
    double gxnew=H(step,6);
    double gynew=H(step,7);
    double gznew=H(step,8);
    if (v==1){
      double gxold=H(step-1,6);
      double gyold=H(step-1,7);
      double gzold=H(step-1,8);
        g1+=i1(gxold,gxnew,told, tnew);
      if (gynew!=0 || gyold!=0){
        g2+=i1(gyold,gynew,told, tnew);
      }
      if (gznew!=0 || gzold!=0){ 
        g3+=i1(gzold,gznew,told, tnew);
      }
      g1static[step-2]=g1;
      g2static[step-2]=g2;
      g3static[step-2]=g3;
    }  
    else{
	g1=g1static[step-2];
	g2=g2static[step-2];
	g3=g3static[step-2];
    } 
    double gg1=g1-grf1;
    double gg2=g2-grf2;
    double gg3=g3-grf3;
    double tt=tnew-trf;//time since the last rf pulse
    cout.precision(20);
    //if (v==1 && readstep%4096==1 && opt_test==1)cout<<"tnew "<<tnew<<";told "<<told<<";trf "<<trf<<";tt=tnew-trf "<<tt<<";tnew-told "<<tnew-told<<endl; 
    //if (v==1 && readstep%4096==1 && opt_test==1)cout<<"g1 "<<g1<<";grf1 "<<grf1<<";gg1 "<<gg1<<endl;
     if (told>=timecourse[actstep] && actstep<=(Nact-2)){
      coeff(activation[actstep],activation[actstep+1],timecourse[actstep],timecourse[actstep+1],dT2_1,dT2_2);
      dT2_1=dT2_1*iT2*iT2;
      dT2_2=dT2_2*iT2*iT2;
      actstep=actstep+1;
    }
    actint+=(dT2_1+dT2_2*(tnew+told)/2)*(tnew-told);
     if (told>=b0timecourse[b0step] && b0step<=(Nb0-2)){
       coeff(b0time[b0step],b0time[b0step+1],b0timecourse[b0step],b0timecourse[b0step+1],db0_1,db0_2);
       coeff(b0xtime[b0step],b0xtime[b0step+1],b0timecourse[b0step],b0timecourse[b0step+1],db0x_1,db0x_2);
       coeff(b0ytime[b0step],b0ytime[b0step+1],b0timecourse[b0step],b0timecourse[b0step+1],db0y_1,db0y_2);
       coeff(b0ztime[b0step],b0ztime[b0step+1],b0timecourse[b0step],b0timecourse[b0step+1],db0z_1,db0z_2);
       if (opt_test==1){
       cout<<"Voxel number="<<v<<"; b0step="<<b0step<<endl;
       cout<<"b0time[b0step]="<<b0time[b0step]<<"b0time[b0step+1]="<<b0time[b0step+1]<<endl;
       cout<<"b0timecourse[b0step]"<<b0timecourse[b0step]<<"b0timecourse[b0step+1]"<<b0timecourse[b0step+1]<<endl;
       cout<<"Coefficients="<<db0_1<<" "<<db0_2<<endl;
       cout<<"b0val="<<db0_1+db0_2*tnew<<endl;
       }
       b0step=b0step+1;
    }
    b0val=db0_1+db0_2*tnew;
    b0int+=(db0_1+db0_2*(tnew+told)/2)*(tnew-told);
    b0xint+=(db0x_1+db0x_2*(tnew+told)/2)*(tnew-told);
    b0yint+=(db0y_1+db0y_2*(tnew+told)/2)*(tnew-told);
    b0zint+=(db0z_1+db0z_2*(tnew+told)/2)*(tnew-told);
    double phase=gama*(gg1*x+gg2*y+gg3*z+b0int+chshift*tt+b0*tt); 
    if (rfangle!=0){
      excitation=0;
      double df=H(step,3);
      double fc=H(step,4);
      double f=gammabar*(gxnew*x+gynew*y+gznew*z+b0val+chshift+b0);
      //if (v==1 && readstep==2081)cout<<"For the first voxel: f_orig=  "<<gammabar*gznew*z<<"  f_b0=  "<<gammabar*b0<<" & df/2= "<<df/2<<endl;
      rftest=0;
      double fval=(f-fc)/df;
      double off=(fval-dslcp_first)/dslcp;
      int nf=(int) off;
      if (nf>=0 && nf<=(Nslc-2)) { 
        off-=nf;
        double ts=table_slcprof[nf];
        double sx=(table_slcprof[nf+1]-ts)*off + ts;
        double rfangle_f=sx*rfangle*RFtrans;//RFtrans are values 0 to 1 to derscribe the inhomogeneity f the receive RF field. 1 is for perfectly homog;
	if (opt_test==1 && fabs(z)<0.0005){
	  cout<<"table_slcprof[nf]= "<<table_slcprof[nf]<<"; rfangle= "<<rfangle<<"; RFtrans= "<<RFtrans<<endl;
	  cout<<"rfangle_f=table_slcprof[nf]*rfangle*RFtrans= "<<rfangle_f<<endl;
	}
        // if (fc-df/2<=f && f<=fc+df/2) {
        if (rfangle_f>0){//new stuff mon dec 19
          excitation=1;
          m=free(m,tt,tissue,phase,actint);
	  //new stuff mon dec 19
	  //due to crushers or any gradient induced dephasing over the voxel a new initial magnetisation is introduced which is the average of the magnetisations over the voxel  
          double xvalrf=fabs(glo_cx*(gg1+b0xint+b0x*tt));
          double yvalrf=fabs(glo_cy*(gg2+b0yint+b0y*tt));
          double zvalrf=fabs(glo_cz*(gg3+b0zint+b0z*tt));
	  double xyzrf=Sinc(xvalrf)*Sinc(yvalrf)*Sinc(zvalrf);
	  m(1)=m(1)*xyzrf;
          m(2)=m(2)*xyzrf;
          m=rot(rfangle_f,"x")*m;
          //new stuff
          m00=sqrt(m(1)*m(1)+m(2)*m(2));
          if (opt_test==1 && fabs(z)<0.0005 ){ 
            cout<<"Projection of the magnetisation vector into the xy plane after flipping is "<<m00<<endl;
	  }
          trf=tnew;
          grf1=g1;
          grf2=g2;
          grf3=g3;
	  actint=0.0;
	  b0int=0.0;
          b0xint=0.0;
          b0yint=0.0;
          b0zint=0.0;
	  rftest=1;
	  //}
	}//new stuff
      }
    }
    if (read!=0){
      readstep=readstep+1;
      if (excitation==1 || nospeedup==1){
	double xval=fabs(glo_cx*(gg1+b0xint+b0x*tt));//b0 TIME
	double yval=fabs(glo_cy*(gg2+b0yint+b0y*tt));//b0 TIME
	double zval=fabs(glo_cz*(gg3+b0zint+b0z*tt));//b0 TIME add b0int like for actint. the same type of term.
      if (v==1 && save_kcoord==1){
        coord(1,readstep)=gammabar*(gg1+b0xint+b0x*tt);//to record the distorted coordinatres in the k-space//b0 TIME
        coord(2,readstep)=gammabar*(gg2+b0yint+b0y*tt);// b0 TIME
        coord(3,readstep)=gammabar*(gg3+b0zint+b0z*tt);
      }
      double tmp;
#ifdef NOTABLE
      tmp=m00*exp(-tt*iT2+actint)*Sinc(xval)*Sinc(yval)*Sinc(zval);
      sreal[readstep-1]+=den*tmp*cos(phase);
      simag[readstep-1]+=den*tmp*sin(phase);
#else
      if (xval<=glo_Dsinc && yval<=glo_Dsinc && zval<=glo_Dsinc){
        //TABLES for SINC
	//x
        double off=xval*glo_idsinc;
	int nbin=(int) off;
	off-=nbin;
	double ts=table_sinc[nbin];
	double sx=(table_sinc[nbin+1]-ts)*off + ts;
	//y
        off=yval*glo_idsinc;
	nbin=(int) off;
	off-=nbin;
	ts=table_sinc[nbin];
	double sy=(table_sinc[nbin+1]-ts)*off + ts;
	//z        
        off=zval*glo_idsinc;
	nbin=(int) off;
	off-=nbin;
	ts=table_sinc[nbin];
	double sz=(table_sinc[nbin+1]-ts)*off + ts;
	//CALCULATING TMP
        tmp=m00*exp(-tt*iT2+actint)*sx*sy*sz;
	//TESTING 
	if (opt_test==1 && v==1 && readstep%4096==2081){
	cout.precision(20);
	cout<<"table_sinc(xval)= "<<sx<<"; table_sinc(yval)= "<<sy<<"; table_sinc(zval)= "<<sz<<endl;
	cout<<"sinc-table(xval)= "<<Sinc(xval)-sx<<"; sinc-table(yval)= "<<Sinc(yval)-sy<<"; sinc-tabel(zval)= "<<Sinc(zval)-sz<<endl;
	}
      }
      else {
        tmp=m00*exp(-tt*iT2+actint)*Sinc(xval)*Sinc(yval)*Sinc(zval);
      }
      //TABLES SIN AND COS CALCULATION
      double phase_2pi;
      if (phase>0) phase_2pi=phase- ((int) (phase*glo_itwopi))*glo_twopi;//one solution when phase exceedes Dsin , this is faster
      else phase_2pi=phase- ((int) (phase*glo_itwopi))*glo_twopi+glo_twopi;
      double off=phase_2pi*glo_idsin;
      int nphase=(int) off;
      off-=nphase;
      double ts1=table_sin[nphase], tc1=table_cos[nphase];
      double wanted_cos=(table_cos[nphase+1]-tc1)*off+tc1;
      double wanted_sin=(table_sin[nphase+1]-ts1)*off+ts1;
      //SIGNAL CALCULATION
      sreal[readstep-1]+=den*tmp*wanted_cos;
      simag[readstep-1]+=den*tmp*wanted_sin;
      //TESTING  
      if(opt_test==1 && v==1 && readstep%4096==2081){
          cout.precision(20);
          cout<<"readstep= "<<readstep<<endl;
          cout<<"gama= "<<gama<<endl;
          cout<<"gg1= "<<gg1<<"; gg2= "<<gg2<<"; gg3= "<<gg3<<endl;
          cout<<"x ="<<x<<"; y= "<<y<<"; z= "<<z<<endl;
          cout<<"b0= "<<b0<<"; chshift= "<<chshift<<endl;
          cout<<"b0val= "<<b0val<<"; b0timecourse "<<b0timecourse<<endl;
          cout<<"b0int= "<<b0int<<endl;
          cout<<"tnew= "<<tnew<<"; trf= "<<trf<<endl;
          cout<<"tt=tnew-trf= "<<tt<<endl;
          cout<<"phase=gama*(gg1*x+gg2*y+gg3*z+(b0+chshift)*tt)= "<<phase<<endl; 
	  cout<<"phase_2pi= "<<phase_2pi<<endl;
	  cout<<"table_cos(phase)= "<<wanted_cos<<"; cos-table(phase)= "<<cos(phase)-wanted_cos<<endl;
          cout<<"table_sin(phase)= "<<wanted_sin<<"; sin-table(phase)= "<<sin(phase)-wanted_sin<<endl;
	  cout<<"glo_cx= "<<glo_cx<<"; glo_cy= "<<glo_cy<<"; glo_cz= "<<glo_cz<<endl;
          cout<<"b0x= "<<b0x<<"; b0y= "<<b0y<<"; b0z= "<<b0z<<endl;
	  cout<<"b0xint= "<<b0xint<<"; b0yint= "<<b0yint<<"; b0zint= "<<b0zint<<endl;
          cout<<"Ival=fabs(glo_cI*(ggI+b0I*tt)); I=x,y,z"<<endl;
          cout<<"xval= "<<xval<<"; yval= "<<yval<<"; zval= "<<zval<<endl;
          cout<<"m00= "<<m00<<endl;
          cout<<"actint= "<<actint<<"; iT2= "<<iT2<<" so exp(-tt*iT2+actint)= "<<exp(-tt*iT2+actint)<<endl;
	  cout<<"tmp= m00*exp(-tt*iT2+actint)*Sinc(xval)*Sinc(yval)*Sinc(zval)= "<<tmp<<endl;
	  cout<<"den= "<<den<<endl;
          cout<<"sreal_1voxel("<<readstep<<")= den*tmp*table_cos[nphase]= "<<sreal[readstep-1]<<endl;
          cout<<"simag_1voxel("<<readstep<<")= den*tmp*table_sin[nphase]= "<<simag[readstep-1]<<endl;
	  cout<<"-------------------------------------------------------------------------------"<<endl;
      }
#endif
      }//end of slcnum..
    }//end of if read
  }//end of main loop
  if (v==1 && save_kcoord==1){
    write_binary_matrix(coord,outputname+"_kcoord" );
  }
}

  PMatrix::PMatrix() { dmat=0; fmat=0; rows=0; cols=0; }

  void PMatrix::initialize(int nrows, int ncols) 
  { 
    rows=nrows; 
    cols=ncols; 
    if ( (cols>1) && (rows>0) ) {
      dmat=new double[rows];
      fmat=new float[(cols-1)*rows];
      for (int n=0; n<rows; n++) dmat[n]=0.0;
    }
  }

  PMatrix::PMatrix(int nrows, int ncols)
  { this->initialize(nrows,ncols); }

  void PMatrix::destroy() 
  {
    if ( (rows>0) && (cols>1) ) {
      delete [] dmat;
      delete [] fmat;
    }
    rows=0;
    cols=0;
  }

  PMatrix::~PMatrix() 
  { 
    this->destroy(); 
  }

  void PMatrix::ReSize(int nrows, int ncols)
  {
    this->destroy();
    this->initialize(nrows, ncols);
  }

  PMatrix& PMatrix::operator=(const PMatrix& src)
  {
    if ( (rows!=src.rows) || (cols!=src.cols) ) 
      { this->initialize(src.rows,src.cols); }
    if (rows>0) {
      if (cols>1) {
	for (int n=0; n<rows; n++) {
	  dmat[n]=src.dmat[n];
	}
	for (int n=0; n<rows*(cols-1); n++) {
	  fmat[n]=src.fmat[n];
	}
      }
    }
    return *this;
  }

  PMatrix& PMatrix::operator=(double val)
  {
     if (rows>0) {
      if (cols>1) {
	for (int n=0; n<rows*(cols-1); n++) {
	  fmat[n]=(float) val;
	}
      }
    }
    return *this;
  }

  double& PMatrix::time(int r)
  {
    return dmat[r-1];
  }

  double PMatrix::time(int r) const
  {
    return dmat[r-1];
  }

  float& PMatrix::operator()(int r, int c)
  {
    return fmat[rows*(c-2) + r-1];
  }

  float PMatrix::operator()(int r, int c) const
  {
    return fmat[rows*(c-2) + r-1];
  }

  float& PMatrix::at(int r, int c)
  {
    if (!this->bounds_check(r,c)) { 
      cerr << "INAVLID ACCESS TO ELEMENT " << r << " , " << c << " OF PMatrix MATRIX" << endl; 
      return fmat[0]; 
    }
    return this->operator()(r,c);
  }

  float PMatrix::at(int r, int c) const
  {
    if (!this->bounds_check(r,c)) { 
      cerr << "INAVLID ACCESS TO ELEMENT " << r << " , " << c << " OF PMatrix MATRIX" << endl; 
      return fmat[0]; 
    }
    return this->operator()(r,c);
  }

  bool PMatrix::bounds_check(int r, int c) const
  {
    if ( (r<=0) || (c<=1) || (r>rows) || (c>cols) ) return false;
    else return true;
  }



////////////////////////////////////////////////////////////////////////////
#define BINFLAGP 43

  int read_binary_matrix(PMatrix& mres, const string& filename)
  {
    if ( filename.size()<1 ) return 1;
    ifstream fs(filename.c_str(), ios::in | ios::binary);
    if (!fs) { 
      cerr << "Could not open matrix file " << filename << endl;
      return 2;
    }
    
    bool swapbytes = false;
    unsigned int testval;
    // test for byte swapping
    fs.read((char*)&testval,sizeof(testval));
    if (testval!=BINFLAGP) {
      swapbytes = true;
      Swap_Nbytes(1,sizeof(testval),&testval);
      if (testval!=BINFLAGP) { 
	cerr << "Unrecognised binary matrix file format" << endl;
	return 2;
      }
    }

    // read matrix dimensions (num rows x num cols)
    unsigned int ival,nx,ny;
    fs.read((char*)&ival,sizeof(ival));
    // ignore the padding (reserved for future use)
    fs.read((char*)&ival,sizeof(ival));
    if (swapbytes) Swap_Nbytes(1,sizeof(ival),&ival);
    nx = ival;
    fs.read((char*)&ival,sizeof(ival));
    if (swapbytes) Swap_Nbytes(1,sizeof(ival),&ival);
    ny = ival;

    // set up and read matrix (rows fast, cols slow)
    double dval;
    float fval;
    if ( (((unsigned int) mres.Ncols())<ny) || (((unsigned int) mres.Nrows())<nx) ) {
      mres.ReSize(nx,ny);
    }
    for (unsigned int y=1; y<=ny; y++) {
      for (unsigned int x=1; x<=nx; x++) {
	if (y==1) {
	  fs.read((char*)&dval,sizeof(dval));
	  if (swapbytes) Swap_Nbytes(1,sizeof(dval),&dval);
	  mres.time(x)=dval;
	} else {
	  fs.read((char*)&fval,sizeof(fval));
	  if (swapbytes) Swap_Nbytes(1,sizeof(fval),&fval);
	  mres(x,y)=fval;
	}
      }
    }
    
    fs.close();
    return 0;
  }




  int write_binary_matrix(const PMatrix& mat, const string& filename)
  {
    Tracer tr("write_binary_matrix");
    if ( (filename.size()<1) ) return -1;
    ofstream fs(filename.c_str(), ios::out | ios::binary);
    if (!fs) { 
      cerr << "Could not open file " << filename << " for writing" << endl;
      return -1;
    }

    unsigned int ival, nx, ny;

    ival = BINFLAGP;
    fs.write((char*)&ival,sizeof(ival));
    ival = 0;  // padding (reserved for future use)
    fs.write((char*)&ival,sizeof(ival));
    ival = mat.Nrows();
    fs.write((char*)&ival,sizeof(ival));
    ival = mat.Ncols();
    fs.write((char*)&ival,sizeof(ival));

    nx = mat.Nrows();
    ny = mat.Ncols();

    double dval;
    float fval;
#ifdef PPC64	
    int n=0;
#endif
    for (unsigned int y=1; y<=ny; y++) {
      for (unsigned int x=1; x<=nx; x++) {
	if (y==1) { 
	  dval = mat.time(x);
	  fs.write((char*)&dval,sizeof(dval));
#ifdef PPC64	
	  if ((n++ % 50) == 0) fs.flush();
#endif
	} else {
	  fval = mat(x,y);
	  fs.write((char*)&fval,sizeof(fval));
#ifdef PPC64	
	  if ((n++ % 50) == 0) fs.flush();
#endif
	}
      }
    }

    fs.close();
    return 0;
  }


int read_binary_matrix(PMatrix& mres, const string& filename,
              int row1=-1, int row2=-1, int col1=-1, int col2=-1)
 {
   if ( filename.size()<1 ) return 1;
   ifstream fs(filename.c_str(), ios::in | ios::binary);
   if (!fs) {
     cerr << "Could not open matrix file " << filename << endl;
     return 2;
   }
   bool swapbytes = false;
   int testval;
   // test for byte swapping
   fs.read((char*)&testval,sizeof(testval));
   if (testval!=BINFLAGP) {
     swapbytes = true;
     Swap_Nbytes(1,sizeof(testval),&testval);
     if (testval!=BINFLAGP) {
   cerr << "Unrecognised binary matrix file format" << endl;
   return 2;
     }
   }

   // read matrix dimensions (num rows x num cols)
   int ival,nx,ny;
   fs.read((char*)&ival,sizeof(ival));
   // ignore the padding (reserved for future use)
   fs.read((char*)&ival,sizeof(ival));
   if (swapbytes) Swap_Nbytes(1,sizeof(ival),&ival);
   nx = ival;
   fs.read((char*)&ival,sizeof(ival));
   if (swapbytes) Swap_Nbytes(1,sizeof(ival),&ival);
   ny = ival;

   int r1=row1, r2=row2, c1=col1, c2=col2;
   if ((r1<1) || (r1>nx)) r1=1;
   if ((r2<1) || (r2>nx)) r2=nx;
   if ((c1<1) || (c1>ny)) c1=1;
   if ((c2<1) || (c2>ny)) c2=ny;

   // set up and read matrix (rows fast, cols slow)
   double dval;
   float fval;
   if ( (((int) mres.Ncols())<(c2-c1+1)) || (((int) mres.Nrows())<(r2-r1+1)) ) {
     mres.ReSize(c2-c1+1,r2-r1+1);
   }
   for (int y=1; y<=ny; y++) {
     // cerr<<" y= "<<y<<endl;
     for (int x=1; x<=nx; x++) {
         if (y==1) {
           fs.read((char*)&dval,sizeof(dval));
	   if ((x>=r1) && (x<=r2) && (y>=c1) && (y<=c2)) {
	     if (swapbytes) Swap_Nbytes(1,sizeof(dval),&dval);
	     mres.time(x-r1+1)=dval;
	   }
	 } else {
	   fs.read((char*)&fval,sizeof(fval));
	   if ((x>=r1) && (x<=r2) && (y>=c1) && (y<=c2)) {
	     if (swapbytes) Swap_Nbytes(1,sizeof(fval),&fval);
	     mres(x-r1+1,y-c1+1)=fval;
	   }
	 }
     }
   }
   fs.close();
   return 0;
 }




