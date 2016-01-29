/*  Copyright (C) 2008 University of Oxford  */

/* S.Jbabdi */

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
#include <cmath>
#include "miscmaths/miscmaths.h"
#include "miscmaths/minimize.h"
#include "newmat.h"
#include "newimage/newimageall.h"
#include "dtifitOptions.h"

using namespace std;
using namespace NEWMAT;
using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace DTIFIT;

const float maxfloat=1e10;
const float minfloat=1e-10;
const float maxlogfloat=23;
const float minlogfloat=-23;
const int maxint=1000000000; 


ReturnMatrix form_Kmat(const Matrix& r){
  Matrix K(r.Ncols(),15);

  Matrix ind(15,4);
  ind << 1 << 1 << 1 << 1
      << 2 << 2 << 2 << 2
      << 3 << 3 << 3 << 3
      << 1 << 1 << 1 << 2
      << 1 << 1 << 1 << 3
      << 1 << 2 << 2 << 2
      << 2 << 2 << 2 << 3
      << 1 << 3 << 3 << 3
      << 2 << 3 << 3 << 3
      << 1 << 1 << 2 << 2
      << 1 << 1 << 3 << 3
      << 2 << 2 << 3 << 3
      << 1 << 1 << 2 << 3
      << 1 << 2 << 2 << 3
      << 1 << 2 << 3 << 3;

  for(int i=1;i<=15;i++){
    for(int j=1;j<=r.Ncols();j++){
      K(j,i) = r((int)ind(i,1),j) * r((int)ind(i,2),j) * r((int)ind(i,3),j) * r((int)ind(i,4),j);
    }
  }

 //  for(int j=1;j<=r.Ncols();j++){
//     float x=r(1,j),y=r(2,j),z=r(3,j);
//     K(j,1) = MISCMATHS::pow(x,4);
//     K(j,2) = MISCMATHS::pow(y,4);
//     K(j,3) = MISCMATHS::pow(z,4);
//     K(j,4) = 4*MISCMATHS::pow(x,3)*y;
//     K(j,5) = 4*MISCMATHS::pow(x,3)*z;
//     K(j,6) = 4*MISCMATHS::pow(y,3)*x;
//     K(j,7) = 4*MISCMATHS::pow(y,3)*z;
//     K(j,8) = 4*MISCMATHS::pow(z,3)*x;
//     K(j,9) = 4*MISCMATHS::pow(z,3)*y;
//     K(j,10) = 6*MISCMATHS::pow(x,2)*MISCMATHS::pow(y,2);
//     K(j,11) = 6*MISCMATHS::pow(x,2)*MISCMATHS::pow(z,2);
//     K(j,12) = 6*MISCMATHS::pow(y,2)*MISCMATHS::pow(z,2);
//     K(j,13) = 12*MISCMATHS::pow(x,2)*y*z;
//     K(j,14) = 12*MISCMATHS::pow(y,2)*x*z;
//     K(j,15) = 12*MISCMATHS::pow(z,2)*x*y;
//     j+=1;
//   }


  K.Release();
  return K;
}


// note the order of the variable parameters
// D11,D12,D13,D22,D23,D33,logS0
// W1111,W2222,W333,W1112,W1113,W1222,W2223,W1333,
// W2333,W1122,W1133,W2233,W1123,W1223,W1233
class KurtosisNonlinCF : public gEvalFunction
{
protected:
  ColumnVector m_A;
  ColumnVector m_B;
  Matrix       m_C;
  Matrix       m_D;
  int          m_n;

public:
  KurtosisNonlinCF(const ColumnVector& data,const Matrix& bvals,const Matrix& bvecs):gEvalFunction()
  {

    m_n = data.Nrows();

    m_A.ReSize(m_n);
    m_B.ReSize(m_n);
    m_C.ReSize(m_n,6);
    m_D.ReSize(m_n,15);

    Matrix K = form_Kmat(bvecs);

    for (int i=1;i<=m_n;i++){
      if(data(i)>0){
	m_A(i)=-log(data(i));
      }
      else{
	m_A(i)=0;
      }
      m_B(i) = 1.0;
      
      m_C(i,1) = -bvals(1,i)*bvecs(1,i)*bvecs(1,i);
      m_C(i,2) = -2*bvals(1,i)*bvecs(1,i)*bvecs(2,i);
      m_C(i,3) = -2*bvals(1,i)*bvecs(1,i)*bvecs(3,i);
      m_C(i,4) = -bvals(1,i)*bvecs(2,i)*bvecs(2,i);
      m_C(i,5) = -2*bvals(1,i)*bvecs(2,i)*bvecs(3,i);
      m_C(i,6) = -bvals(1,i)*bvecs(3,i)*bvecs(3,i);

      for(int j=1;j<=15;j++)
	m_D(i,j) = (bvals(1,i)*bvals(1,i)/6) * K(i,j) / 9;

    }

    //        m_D=0;
    
  }

  virtual ~KurtosisNonlinCF(){};

  float evaluate(const ColumnVector& x) const{
    float res=0;
    res = ( m_A + m_B*x(7) + m_C*x.SubMatrix(1,6,1,1)
	    + m_D*x.SubMatrix(8,22,1,1)*(x(1)+x(4)+x(6))*(x(1)+x(4)+x(6))).SumSquare();

    OUT(x.t());
    OUT(res);

    return res;
  }
  ReturnMatrix g_evaluate(const ColumnVector& x) const{
    ColumnVector sj_g(x.Nrows());
    //    ColumnVector sj_gg;
    
    //  sj_gg = MISCMATHS::gradient(x,*this,1e-4);
 
    ColumnVector sj_d(6);
    ColumnVector sj_w(15);
    sj_d = x.SubMatrix(1,6,1,1);
    sj_w = x.SubMatrix(8,22,1,1);
    double sj_t = x(1)+x(4)+x(6);
    double sj_t2=sj_t*sj_t;

    ColumnVector sj_func(m_n);
    sj_func = m_A + m_B*x(7) + m_C*sj_d + m_D*sj_w*sj_t2;

    sj_g(1) = 2*NEWMAT::SP(sj_func,m_C.Column(1)+2*sj_t*m_D*sj_w).Sum();
    sj_g(2) = 2*NEWMAT::SP(sj_func,m_C.Column(2)).Sum();
    sj_g(3) = 2*NEWMAT::SP(sj_func,m_C.Column(3)).Sum();
    sj_g(4) = 2*NEWMAT::SP(sj_func,m_C.Column(4)+2*sj_t*m_D*sj_w).Sum();
    sj_g(5) = 2*NEWMAT::SP(sj_func,m_C.Column(5)).Sum();
    sj_g(6) = 2*NEWMAT::SP(sj_func,m_C.Column(6)+2*sj_t*m_D*sj_w).Sum();

    sj_g(7) = 2*NEWMAT::SP(sj_func,m_B).Sum();
 
    for(int i=1,j=8;j<=x.Nrows();i++,j++)
      sj_g(j) = 2*NEWMAT::SP(sj_func,sj_t2*m_D.Column(i)).Sum();

    OUT(sj_g.t());

     sj_g.Release();
     return sj_g;
     //    sj_gg.Release();
     //return sj_gg;

  }

  const KurtosisNonlinCF& operator=(const KurtosisNonlinCF& par)
  {
    m_A   = par.m_A;
    m_B   = par.m_B;
    m_C   = par.m_C;
    m_D   = par.m_D;
    m_n = par.m_n;

    return *this;
    
  }
  KurtosisNonlinCF(const KurtosisNonlinCF& rhs):
    m_A(rhs.m_A),m_B(rhs.m_B),m_C(rhs.m_C),m_D(rhs.m_D),m_n(rhs.m_n){
    *this=rhs;
  }

};

inline Matrix Anis()
{ 
  Matrix A(3,3);
  A << 1 << 0 << 0
    << 0 << 0 << 0
    << 0 << 0 << 0;
  return A;
}

Matrix form_Amat(const Matrix& r,const Matrix& b)
{
  Matrix A(r.Ncols(),7);
  Matrix tmpvec(3,1), tmpmat;
  
  for( int i = 1; i <= r.Ncols(); i++){
    tmpvec << r(1,i) << r(2,i) << r(3,i);
    tmpmat = tmpvec*tmpvec.t()*b(1,i);
    A(i,1) = tmpmat(1,1);
    A(i,2) = 2*tmpmat(1,2);
    A(i,3) = 2*tmpmat(1,3);
    A(i,4) = tmpmat(2,2);
    A(i,5) = 2*tmpmat(2,3);
    A(i,6) = tmpmat(3,3);
    A(i,7) = 1;
  }
  return A;
}
inline SymmetricMatrix vec2tens(ColumnVector& Vec){
  SymmetricMatrix tens(3);
  tens(1,1)=Vec(1);
  tens(2,1)=Vec(2);
  tens(3,1)=Vec(3);
  tens(2,2)=Vec(4);
  tens(3,2)=Vec(5);
  tens(3,3)=Vec(6);
  return tens;
}

void kurtosisfit(DiagonalMatrix& Dd,ColumnVector& evec1,ColumnVector& evec2, ColumnVector& evec3,
		 float& f,float& s0,ColumnVector& Dvec, float& mk, ColumnVector& tens4, 
		 const Matrix& Amat,const Matrix& Kmat,const ColumnVector& S,const Matrix& bvals,const Matrix& bvecs){


  // calculate DT and KT using non-linear fitting
  KurtosisNonlinCF KNL(S,bvals,bvecs);
  ColumnVector xmin(22);
  xmin=0.0;
  xmin << .002 << 0 << 0 << .001 << 0 << .001 << 1000 << .5 << 0.5 << 0.1 << 0 << 0 << 0 << 0 << 0 << 0 <<0<<0<<0<<0<<0<<0;
  KNL.minimize(xmin);


  Dvec.SubMatrix(1,6,1,1) = xmin.SubMatrix(1,6,1,1);
  tens4 = xmin.SubMatrix(8,22,1,1);
  Dvec(7) = exp(xmin(7));  
  s0 = Dvec(7);

  // Tensor Stuff
  float mDd, fsquared;
  SymmetricMatrix tens;
  DiagonalMatrix Ddsorted(3);
  Matrix Vd;
  tens = vec2tens(Dvec);
  EigenValues(tens,Dd,Vd);
  mDd = Dd.Sum()/Dd.Nrows();
  int maxind = Dd(1) > Dd(2) ? 1:2;   //finding max,mid and min eigenvalues
  maxind = Dd(maxind) > Dd(3) ? maxind:3;
  int midind;
  if( (Dd(1)>=Dd(2) && Dd(2)>=Dd(3)) || (Dd(1)<=Dd(2) && Dd(2)<=Dd(3)) ){midind=2;}
  else if( (Dd(2)>=Dd(1) && Dd(1)>=Dd(3)) || (Dd(2)<=Dd(1) && Dd(1)<=Dd(3)) ){midind=1;}
  else {midind=3;}
  int minind = Dd(1) < Dd(2) ? 1:2;   //finding maximum eigenvalue
  minind = Dd(minind) < Dd(3) ? minind:3;
  Ddsorted << Dd(maxind) << Dd(midind) << Dd(minind);
  Dd=Ddsorted;
  evec1 << Vd(1,maxind) << Vd(2,maxind) << Vd(3,maxind);
  evec2 << Vd(1,midind) << Vd(2,midind) << Vd(3,midind);
  evec3 << Vd(1,minind) << Vd(2,minind) << Vd(3,minind);
  float numer=1.5*((Dd(1)-mDd)*(Dd(1)-mDd)+(Dd(2)-mDd)*(Dd(2)-mDd)+(Dd(3)-mDd)*(Dd(3)-mDd));
  float denom=(Dd(1)*Dd(1)+Dd(2)*Dd(2)+Dd(3)*Dd(3));
  if(denom>0) fsquared=numer/denom;
  else fsquared=0;
  if(fsquared>0){f=sqrt(fsquared);}
  else{f=0;}

  // Kurtosis Stuff
  mk = 0;
  ColumnVector vec(S.Nrows());
  vec = Kmat*tens4;
  for(int i=1;i<=S.Nrows();i++){
    if(bvals(1,i)>0)
      mk+=vec(i)/(bvecs.Column(i).t()*tens*bvecs.Column(i)).AsScalar();
  }
  mk *= mDd*mDd;
  mk /= float(S.Nrows());
  
}

int main(int argc, char** argv)
{
  //parse command line
  dtifitOptions& opts = dtifitOptions::getInstance();
  int success=opts.parse_command_line(argc,argv);
  if(!success) return 1;

  if(opts.verbose.value()){
    cout<<"data file "<<opts.dtidatafile.value()<<endl;
    cout<<"mask file "<<opts.maskfile.value()<<endl;
    cout<<"bvecs     "<<opts.bvecsfile.value()<<endl;
    cout<<"bvals     "<<opts.bvalsfile.value()<<endl;
    if(opts.littlebit.value()){
      cout<<"min z     "<<opts.z_min.value()<<endl;
      cout<<"max z     "<<opts.z_max.value()<<endl;
      cout<<"min y     "<<opts.y_min.value()<<endl;
      cout<<"max y     "<<opts.y_max.value()<<endl;
      cout<<"min x     "<<opts.x_min.value()<<endl;
      cout<<"max x     "<<opts.x_max.value()<<endl;
    }
  }
  /////////////////////////////////////////
  // read bvecs and bvals
  // correct transpose and normalise bvecs
  Matrix r = read_ascii_matrix(opts.bvecsfile.value());
  if(r.Nrows()>3) r=r.t();
  for(int i=1;i<=r.Ncols();i++){
    float tmpsum=sqrt(r(1,i)*r(1,i)+r(2,i)*r(2,i)+r(3,i)*r(3,i));
    if(tmpsum!=0){
      r(1,i)=r(1,i)/tmpsum;
      r(2,i)=r(2,i)/tmpsum;
      r(3,i)=r(3,i)/tmpsum;
    }  
  }
  Matrix b = read_ascii_matrix(opts.bvalsfile.value());
  if(b.Nrows()>1) b=b.t();
  //////////////////////////////////////////

  volume4D<float> data;
  volume<int> mask;
  if(opts.verbose.value()) cout<<"reading data"<<endl;
  read_volume4D(data,opts.dtidatafile.value());
  if(opts.verbose.value()) cout<<"reading mask"<<endl;
  read_volume(mask,opts.maskfile.value());
  if(opts.verbose.value()) cout<<"ok"<<endl;
  int minx=opts.littlebit.value() ? opts.x_min.value():0;
  int maxx=opts.littlebit.value() ? opts.x_max.value():mask.xsize();
  int miny=opts.littlebit.value() ? opts.y_min.value():0;
  int maxy=opts.littlebit.value() ? opts.y_max.value():mask.ysize();
  int minz=opts.littlebit.value() ? opts.z_min.value():0;
  int maxz=opts.littlebit.value() ? opts.z_max.value():mask.zsize();
  cout<<minx<<" "<<maxx<<" "<<miny<<" "<<maxy<<" "<<minz<<" "<<maxz<<endl;
  if(opts.verbose.value()) cout<<"setting up vols"<<endl;

  volume<float>   l1(maxx-minx,maxy-miny,maxz-minz);
  volume<float>   l2(maxx-minx,maxy-miny,maxz-minz);
  volume<float>   l3(maxx-minx,maxy-miny,maxz-minz);
  volume<float>   MD(maxx-minx,maxy-miny,maxz-minz);
  volume<float>   FA(maxx-minx,maxy-miny,maxz-minz);
  volume<float>   S0(maxx-minx,maxy-miny,maxz-minz);
  volume4D<float> V1(maxx-minx,maxy-miny,maxz-minz,3);
  volume4D<float> V2(maxx-minx,maxy-miny,maxz-minz,3);
  volume4D<float> V3(maxx-minx,maxy-miny,maxz-minz,3);
  volume4D<float> Delements(maxx-minx,maxy-miny,maxz-minz,6);
  volume<float>   MK(maxx-minx,maxy-miny,maxz-minz);
  volume4D<float> KurtTens(maxx-minx,maxy-miny,maxz-minz,15);

  if(opts.verbose.value()) cout<<"copying input properties to output volumes"<<endl;
  copybasicproperties(data[0],l1);
  copybasicproperties(data[0],l2);
  copybasicproperties(data[0],l3);
  copybasicproperties(data[0],MD);
  copybasicproperties(data[0],FA);
  copybasicproperties(data[0],S0);
  copybasicproperties(data[0],V1[0]);
  copybasicproperties(data[0],V2[0]);
  copybasicproperties(data[0],V3[0]);
  copybasicproperties(data[0],Delements[0]);
  copybasicproperties(data[0],MK);
  copybasicproperties(data[0],KurtTens[0]);

  if(opts.verbose.value()) cout<<"zeroing output volumes"<<endl;
  l1=0;l2=0;l3=0;MD=0;FA=0;S0=0;V1=0;V2=0;V3=0;Delements=0;MK=0;KurtTens=0;
  if(opts.verbose.value()) cout<<"ok"<<endl;


  DiagonalMatrix evals(3);
  ColumnVector evec1(3),evec2(3),evec3(3);
  ColumnVector tens4(15);
  ColumnVector S(data.tsize());
  float fa,s0,mk;
  if(opts.verbose.value()) cout<<"Forming A matrix"<<endl;
  Matrix Amat = form_Amat(r,b);
  Matrix Kmat = form_Kmat(r);
  if(opts.verbose.value()) cout<<"starting the fits"<<endl;
  ColumnVector Dvec(7); Dvec=0;
  for(int k = minz; k < maxz; k++){
    cout<<k<<" slices processed"<<endl;
      for(int j=miny; j < maxy; j++){
	for(int i =minx; i< maxx; i++){
	
	  if(mask(i,j,k)>0){
	    
	    for(int t=0;t < data.tsize();t++){
	      S(t+1)=data(i,j,k,t);
	    }
	    //tensorfit(evals,evec1,evec2,evec3,fa,s0,Dvec,Amat,S);
	    kurtosisfit(evals,evec1,evec2,evec3,fa,s0,Dvec,mk,tens4,Amat,Kmat,S,b,r);


	    l1(i-minx,j-miny,k-minz)=evals(1);
	    l2(i-minx,j-miny,k-minz)=evals(2);
	    l3(i-minx,j-miny,k-minz)=evals(3);
	    MD(i-minx,j-miny,k-minz)=(evals(1)+evals(2)+evals(3))/3;
	    FA(i-minx,j-miny,k-minz)=fa;
	    S0(i-minx,j-miny,k-minz)=s0;
	    V1(i-minx,j-miny,k-minz,0)=evec1(1);
	    V1(i-minx,j-miny,k-minz,1)=evec1(2);
	    V1(i-minx,j-miny,k-minz,2)=evec1(3);
	    V2(i-minx,j-miny,k-minz,0)=evec2(1);
	    V2(i-minx,j-miny,k-minz,1)=evec2(2);
	    V2(i-minx,j-miny,k-minz,2)=evec2(3);
	    V3(i-minx,j-miny,k-minz,0)=evec3(1);
	    V3(i-minx,j-miny,k-minz,1)=evec3(2);
	    V3(i-minx,j-miny,k-minz,2)=evec3(3);
	    Delements(i-minx,j-miny,k-minz,0)=Dvec(1);
	    Delements(i-minx,j-miny,k-minz,1)=Dvec(2);
	    Delements(i-minx,j-miny,k-minz,2)=Dvec(3);
	    Delements(i-minx,j-miny,k-minz,3)=Dvec(4);
	    Delements(i-minx,j-miny,k-minz,4)=Dvec(5);
	    Delements(i-minx,j-miny,k-minz,5)=Dvec(6);

	    MK(i-minx,j-miny,k-minz)=mk;
	    for(int iii=0;iii<15;iii++)
	      KurtTens(i-minx,j-miny,k-minz,iii) = tens4(iii+1);

	  }
	}
      }
  }
  
    string fafile=opts.ofile.value()+"_FA";
    string s0file=opts.ofile.value()+"_S0";
    string l1file=opts.ofile.value()+"_L1";
    string l2file=opts.ofile.value()+"_L2";
    string l3file=opts.ofile.value()+"_L3";
    string v1file=opts.ofile.value()+"_V1";
    string v2file=opts.ofile.value()+"_V2";
    string v3file=opts.ofile.value()+"_V3";
    string MDfile=opts.ofile.value()+"_MD";
    string tensfile=opts.ofile.value()+"_tensor";
    string MKfile=opts.ofile.value()+"_MK";
    string kurtosisfile=opts.ofile.value()+"_kurtosis";
    if(opts.littlebit.value()){
      fafile+="littlebit";
      s0file+="littlebit";
      l1file+="littlebit";
      l2file+="littlebit";
      l3file+="littlebit";
      v1file+="littlebit";
      v2file+="littlebit";
      v3file+="littlebit";
      MDfile+="littlebit";
      tensfile+="littlebit";
      MKfile+="littlebit";
      kurtosisfile+="littlebit";
    }
  
    save_volume(FA,fafile);
    save_volume(S0,s0file);
    save_volume(l1,l1file);
    save_volume(l2,l2file);
    save_volume(l3,l3file);
    save_volume(MD,MDfile);
    save_volume4D(V1,v1file);
    save_volume4D(V2,v2file);
    save_volume4D(V3,v3file);
    save_volume(MK,MKfile);

    if(opts.savetensor.value()){
      save_volume4D(Delements,tensfile);
      save_volume4D(KurtTens,kurtosisfile);
    }
      
    
  return 0;
}













