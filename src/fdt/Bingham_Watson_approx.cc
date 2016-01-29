/*  Directional Statistics Functions

    Bingham and Watson Distributions and functions to approximate their normalizing constant
    
    Stam Sotiropoulos  - FMRIB Image Analysis Group
 
    Copyright (C) 2011 University of Oxford  */

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


  //Cubic root
  float croot(const float& x){
    float res;
    if (x>=0) 
      res=pow(x,INV3);
    else
      res=-(pow(-x,INV3));
    return res;
  }


  //Saddlepoint approximation of confluent hypergeometric function of a matrix argument B (Kume & Wood, 2005)
  //Vector x has the eigenvalues of B.
  float hyp_Sapprox(ColumnVector &x){
    float c1;

    //SortDescending(x);         //Not needed??
    if (x(1)==0 && x(2)==0 && x(3)==0)
      c1=1;
    else{
      float t=find_t(x);
      float T,R=1, K2=0, K3=0, K4=0;
      for (int i=1; i<=3; i++){
        R=R/sqrt(-x(i)-t);
	K2=K2+1.0/(2*pow(x(i)+t,2.0));
	K3=K3-1.0/pow(x(i)+t,3.0);
	K4=K4+3.0/pow(x(i)+t,4.0);
      }
      T=K4/(8*K2*K2)-5*K3*K3/(24*K2*K2*K2);  
      //c1=sqrt(2.0/K2)*pi*R*exp(-t);
      //float c3=c1*exp(T);
      //c1=c3/(4*pi);
      c1=0.25*sqrt(2.0/K2)*R*exp(-t+T);
    }
    return c1;
  }


  //Saddlepoint approximation of confluent hypergeometric function of a matrix argument, with its eigenvalues being l1,l1,l2 or l1,l2,l2 with l1!=l2.
  //Vector x has the three eigenvalues. This function can be also used to approximate a confluent hypergeometric function of a scalar argument k
  //by providing x=[k 0 0].
  float hyp_Sapprox_twoequal(ColumnVector& x){
    float c1, R, t, Bm, Bp, K2, K3, K4,T,q;

    //SortDescending(x);  //Not needed??
    if (x(1)==x(2)){
      q=1;
      x(2)=x(3);
    }
    else 
      q=2;

    R=sqrt(4*(x(1)-x(2))*(x(1)-x(2))+9+4*(2*q-3)*(x(2)-x(1)));
    t=(-2*(x(1)+x(2))-3-R)/4.0;
    Bm=1/(R+3-2*(x(2)-x(1)));
    Bp=1/(R+3+2*(x(2)-x(1)));
    K2=1.0/(2*pow(x(1)+t,2.0))+1.0/pow(x(2)+t,2.0);
    K3=-1.0/pow(x(1)+t,3.0)-2.0/pow(x(2)+t,3.0);
    K4=3.0/pow(x(1)+t,4.0)+6.0/pow(x(2)+t,4.0);
    T=K4/(8.0*K2*K2)-5.0*K3*K3/(24.0*K2*K2*K2);
    c1=sqrt(pow(Bm,q-1)*pow(Bp,2-q)/R)*exp(-t+T);
    return c1;
  }



 //Saddlepoint approximation of the ratio of two hypergeometric functions, with matrix arguments L and B (3x3). Vectors xL & xB contain the eigenvalues of L and B. 
 //Used for the ball & Binghams model, B has two non-zero eigenvalues that represent fanning indices.
  float hyp_SratioB(ColumnVector& xL,ColumnVector& xB){
    float c, T1, t1, Norm1, T2,t2,Norm2;
    float R,K2,K3,K4,tmp,tmp2,tmp3;

    //Approximate Numerator 
    if (xL(1)==0 && xL(2)==0 && xL(3)==0){
      Norm1=1; t1=0; T1=0; 
    }
    else{
      //SortDescending(xL);  //Not needed, performed by eigen-decomposition?
      t1=find_t(xL);
      R=1; K2=0; K3=0; K4=0; 
      for (int i=1; i<=3; i++){
        tmp=xL(i)+t1;
        tmp2=tmp*tmp; tmp3=tmp2*tmp;
        R=-R*tmp;
        K2=K2+0.5/tmp2;
        K3=K3-1.0/tmp3;
        K4=K4+3.0/(tmp3*tmp);
      }
      T1=K4/(8*K2*K2)-5*K3*K3/(24.0*K2*K2*K2);  
      Norm1=K2*R;
    }
    
    //Approximate Denominator
    //SortDescending(xB);  //Not needed, performed by eigen-decomposition?
    t2=find_t(xB);
    R=1; K2=0; K3=0; K4=0; 
    for (int i=1; i<=3; i++){
      tmp=xB(i)+t2;
      tmp2=tmp*tmp; tmp3=tmp2*tmp;
      R=-R*tmp;
      K2=K2+0.5/tmp2;
      K3=K3-1.0/tmp3;
      K4=K4+3.0/(tmp3*tmp);
    }
    T2=K4/(8*K2*K2)-5*K3*K3/(24*K2*K2*K2);  
    Norm2=K2*R;

    //Final Ratio
    c=sqrt(Norm2/Norm1)*exp(-t1+T1+t2-T2);
    return c;
  }


  //Saddlepoint aproximation of the ratio ot two hypergeometric functions with matrix arguments L and B in two steps: First denominator, then numerator. 
  //This allows them to be updated independently, used for the ball & Binghams model to compute the likelihood faster.
  //This function returns values used in the denominator approximation. xB containes the two non-zero eigenvalues of matrix B.
  ReturnMatrix approx_denominatorB(ColumnVector& xB){
    float t2,R,K2,K3,K4,tmp,tmp2,tmp3, T2, Norm2;
    ColumnVector Res(3);

    SortDescending(xB); //Not needed?
    t2=find_t(xB);
    R=1; K2=0; K3=0; K4=0; 
    for (int i=1; i<=3; i++){
      tmp=xB(i)+t2;
      tmp2=tmp*tmp; tmp3=tmp2*tmp;
      R=-R*tmp;
      K2=K2+0.5/tmp2;
      K3=K3-1.0/tmp3;
      K4=K4+3.0/(tmp3*tmp);
    }
    T2=K4/(8.0*K2*K2)-5.0*K3*K3/(24.0*K2*K2*K2);  
    Norm2=K2*R;
    Res<< t2 << T2 << Norm2;
    
    Res.Release();
    return Res;
  }


  //Second step for saddlepoint approximation of the ratio of two hypergeometric functions with matrix arguments L and B (xL has the eigenvalues of L). 
  //Assume that the denominator has already been approximated by the function above and the parameters are stored in denomvals.
  //Here approximate the numerator and return the total ratio approximation.
  float hyp_SratioB_knowndenom(ColumnVector &xL,ColumnVector& denomvals){
    float c,Norm1,t1,T1,R,K2,K3,K4,tmp,tmp2,tmp3;

    //Approximate Numerator 
    SortDescending(xL); //Not needed?
    if (xL(1)==0 && xL(2)==0 && xL(3)==0){
      Norm1=1; t1=0; T1=0; 
    }
    else{
      t1=find_t(xL);
      R=1; K2=0; K3=0; K4=0; 
      for (int i=1; i<=3; i++){
        tmp=xL(i)+t1;
        tmp2=tmp*tmp; tmp3=tmp2*tmp;
        R=-R*tmp;
        K2=K2+0.5/tmp2;
        K3=K3-1.0/tmp3;
        K4=K4+3.0/(tmp3*tmp);
      }
      T1=K4/(8.0*K2*K2)-5.0*K3*K3/(24.0*K2*K2*K2);  
      Norm1=K2*R;
    }
    //Final Ratio
    if (denomvals(3)/Norm1>0)  
      c=sqrt(denomvals(3)/Norm1)*exp(-t1+T1+denomvals(1)-denomvals(2));
    else
      c=0; //Signify that the approximation has failed. c should be >=0. How to treat these voxels???
    if (isinf(c))  //In this case again the approximation has failed 
      c=1e10;
    return c;
  }



  //Saddlepoint approximation of the ratio of two hypergeometric functions, one with matrix argument L and another with scalar argument k. Vector xL contains the eigenvalues of L.
  //Used for the ball & Watsons model, L has up to two non-zero eigenvalues and k!=0 is the fanning index.  
  float hyp_SratioW(ColumnVector& xL,const double k){
    float Norm1,t1,T1,R,K2,K3,K4,tmp,tmp2,tmp3,Rf, t2,T2, Norm2,Bm,c;

    //Approximate Numerator 
    if (xL(1)==0 && xL(2)==0){
      Norm1=1; t1=0; T1=0; 
    }
    else{
      SortDescending(xL);  //Not needed, performed by eigen-decomposition?
      t1=find_t(xL);
      R=1; K2=0; K3=0; K4=0; 
      for (int i=1; i<=3; i++){
	tmp=xL(i)+t1;
	tmp2=tmp*tmp; tmp3=tmp2*tmp;
        R=-R*tmp;
        K2=K2+0.5/tmp2;
        K3=K3-1.0/tmp3;
        K4=K4+3.0/(tmp3*tmp);
	}
      T1=K4/(8.0*K2*K2)-5.0*K3*K3/(24.0*K2*K2*K2);  
      Norm1=0.125/(K2*R);
    }
    
    //Approximate Denominator
    R=sqrt(4*k*k+9-4*k);
    Rf=R+3+2*k;
    t2=-0.25*Rf;
    tmp=k+t2; tmp2=tmp*tmp; tmp3=tmp2*tmp;
    Bm=1.0/Rf;
    K2=0.5/tmp2+1.0/(t2*t2);
    K3=-1.0/tmp3-2.0/(t2*t2*t2);
    K4=3.0/(tmp3*tmp)+6.0/(t2*t2*t2*t2);
    T2=K4/(8.0*K2*K2)-5.0*K3*K3/(24.0*K2*K2*K2);
    Norm2=Bm/R;
    
    //Final Ratio
    c=sqrt(Norm1/Norm2)*exp(-t1+T1+t2-T2);
    return c;
  }


  //Saddlepoint aproximation of the ratio ot two hypergeometric functions, one with matrix arguments L and the other with scalar argument k>0 in two steps: 
  //First denominator, then numerator. This allows them to be updated independently, used for the ball & Watsons model to compute the likelihood faster.
  //This function returns values used in the denominator approximation.
  ReturnMatrix approx_denominatorW(const double k){
    float t2,R,Bm,Rf,K2,K3,K4,tmp,tmp2,tmp3, T2, Norm2;
    ColumnVector Res(3);

    R=sqrt(4*k*k+9-4*k);
    Rf=R+3+2*k;
    t2=-0.25*Rf;
    tmp=k+t2; tmp2=tmp*tmp; tmp3=tmp2*tmp;
    Bm=1.0/Rf;
    K2=0.5/tmp2+1.0/(t2*t2);
    K3=-1.0/tmp3-2.0/(t2*t2*t2);
    K4=3.0/(tmp3*tmp)+6.0/(t2*t2*t2*t2);
    T2=K4/(8.0*K2*K2)-5.0*K3*K3/(24.0*K2*K2*K2);
    Norm2=Bm/R;
    Res<< t2 << T2 << Norm2;

    Res.Release();
    return Res;
  }


  //Second step for saddlepoint approximation of the ratio of two hypergeometric functions, with matrix argument L and scalar argument k (xL has the eigenvalues of L). 
  //Assume that the denominator has already been approximated by the function above and the parameters are stored in denomvals.
  //Here approximate the numerator and return the total ratio approximation.
  float hyp_SratioW_knowndenom(ColumnVector &xL,ColumnVector& denomvals){
    float c,Norm1,t1,T1,R,K2,K3,K4,tmp,tmp2,tmp3;

    //Approximate Numerator 
    if (xL(1)==0 && xL(2)==0){
      Norm1=1; t1=0; T1=0; 
    }
    else{
      SortDescending(xL); 
      t1=find_t(xL);
      R=1; K2=0; K3=0; K4=0; 
      for (int i=1; i<=3; i++){
	tmp=xL(i)+t1;
	tmp2=tmp*tmp; tmp3=tmp2*tmp;
        R=-R*tmp;
        K2=K2+0.5/tmp2;
        K3=K3-1.0/tmp3;
        K4=K4+3.0/(tmp3*tmp);
	}
      T1=K4/(8.0*K2*K2)-5.0*K3*K3/(24.0*K2*K2*K2);  
      Norm1=0.125/(K2*R);
    }
    
    //Final Ratio
    c=sqrt(Norm1/denomvals(3))*exp(-t1+T1+denomvals(1)-denomvals(2));
    return c;
  }



  //Using the values of vector x, construct a qubic equation and solve it analytically.
  //Solution used for the saddlepoint approximation of the confluent hypergeometric function with matrix argument B (3x3) (See Kume & Wood, 2005)
  //Vector x contains the eigenvalues of B.
  float find_t(const ColumnVector& x){

    float l1=-x(1), l2=-x(2), l3=-x(3);
    float a0,a1,a2,a3,ee,tmp,z1,z2,z3,p,q,D,offset;

    a3=l1*l2+l2*l3+l1*l3;
    a2=1.5-l1-l2-l3;
    a1=a3-l1-l2-l3;
    a0=0.5*(a3-2*l1*l2*l3);

    p=(a1-a2*a2*INV3)*INV3;
    q=(-9*a2*a1+27*a0+2*a2*a2*a2)*INV54;
    D=q*q+p*p*p;
    offset=a2*INV3;
    if (D>0){
      ee=sqrt(D);
      tmp=-q+ee; z1=croot(tmp);
      tmp=-q-ee; z1=z1+croot(tmp);
      z1=z1-offset; z2=z1; z3=z1;
    }
    else if (D<0){
      ee=sqrt(-D);
      float tmp2=-q; 
      float angle=2*INV3*atan(ee/(sqrt(tmp2*tmp2+ee*ee)+tmp2));
      tmp=cos(angle);
      tmp2=sin(angle);
      ee=sqrt(-p);
      z1=2*ee*tmp-offset; 
      z2=-ee*(tmp+SQRT3*tmp2)-offset; 
      z3=-ee*(tmp-SQRT3*tmp2)-offset; 
    }
    else{
      tmp=-q;
      tmp=croot(tmp);
      z1=2*tmp-offset; z2=z1; z3=z1;
      if (p!=0 || q!=0)
	z2=-tmp-offset; z3=z2;
    }

    z1=min3(z1,z2,z3);  //Pick the smallest root

    //    if (z2<=l1)
    // z1=z2;
    //else if (z3<=l1)
    //  z1=z3;
    return z1;
  }






