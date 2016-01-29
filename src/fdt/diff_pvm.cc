/* Diffusion Partial Volume Model  

    Tim Behrens - FMRIB Image Analysis Group

    Copyright (C) 2002 University of Oxford  */

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
#include <fstream>
#include <iomanip>
#include <sstream>
#define WANT_STREAM
#define WANT_MATH
//  #include "newmatap.h"
//  #include "newmatio.h"
#include <string>
#include <math.h>
#include "utils/log.h"
#include "diff_pvmoptions.h"
#include "utils/tracer_plus.h"
#include "miscmaths/miscprob.h"
#include "stdlib.h"
#include "bint/model.h"
#include "bint/lsmcmcmanager.h"
#include "bint/lslaplacemanager.h"

using namespace Bint;
using namespace Utilities;
using namespace NEWMAT;
using namespace MISCMATHS;



const float maxfloat=1e10;
const float minfloat=1e-10;
const float maxlogfloat=23;
const float minlogfloat=-23;


inline float min(float a,float b){
  return a<b ? a:b;}
inline float max(float a,float b){
  return a>b ? a:b;}
inline Matrix Anis()
{ 
  Matrix A(3,3);
  A << 1 << 0 << 0
    << 0 << 0 << 0
    << 0 << 0 << 0;
  return A;
}

inline Matrix Is()
{ 
  Matrix I(3,3);
  I << 1 << 0 << 0
    << 0 << 1 << 0
    << 0 << 0 << 1;
  return I;
}

inline ColumnVector Cross(const ColumnVector& A,const ColumnVector& B)
{
  ColumnVector res(3);
  res << A(2)*B(3)-A(3)*B(2)
      << A(3)*B(1)-A(1)*B(3)
      << A(1)*B(2)-B(1)*A(2);
  return res;
}

inline Matrix Cross(const Matrix& A,const Matrix& B)
{
  Matrix res(3,1);
  res << A(2,1)*B(3,1)-A(3,1)*B(2,1)
      << A(3,1)*B(1,1)-A(1,1)*B(3,1)
      << A(1,1)*B(2,1)-B(1,1)*A(2,1);
  return res;
}

float mod(float a, float b){
  while(a>b){a=a-b;}
  while(a<0){a=a+b;} 
  return a;
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


class Diff_pvmModel : public ForwardModel  
  {
  public:
   
    Diff_pvmModel(const Matrix& pbvecs,const Matrix& pbvals,int pdebuglevel)
      : ForwardModel(pdebuglevel), r(pbvecs) , b(pbvals), alpha(pbvals.Ncols()), beta(pbvals.Ncols()), debuglevel(pdebuglevel) 
	
    {
      Amat=form_Amat(r,b);
      cart2sph(r,alpha,beta);
    }
    
    ~Diff_pvmModel(){}
  
    virtual void setparams();
    ReturnMatrix nonlinearfunc(const ColumnVector& paramvalues) const; 
    void initialise(const ColumnVector& S);
    
    
  protected:
    
    const Matrix& r;
    const Matrix& b;
    ColumnVector alpha;
    ColumnVector beta;
    Matrix Amat;
    int debuglevel;
};  

void Diff_pvmModel::setparams()
  {
    Tracer_Plus tr("Diff_pvmModel::setdata");
    if(debuglevel>2){
      cout << "Diff_pvmModel::setparams"<<endl;
    }
    clear_params();
  
    SinPrior thtmp(1,-1000*M_PI,1000*M_PI);
    add_param("th",0.2,0.02,thtmp,true,true); //Will unwrap th param before saving
    UnifPrior phtmp(-1000*M_PI,1000*M_PI);
    add_param("ph",0.2, 0.02,phtmp,true,true); //Will unwrap th param before saving
    UnifPrior ftmp(0,1);
    add_param("f",0.5,0.02,ftmp,true,true);
    GammaPrior dtmp(4,1.0/0.0003); //test this out,
    add_param("d",0.005,0.00005,dtmp,true,true);
    UnifPrior S0tmp(0,100000);
    add_param("S0",10000,100,S0tmp,true,true);//false);
    
  }

ReturnMatrix Diff_pvmModel::nonlinearfunc(const ColumnVector& paramvalues) const
  {
    Tracer_Plus trace("Diff_pvmModel::nonlinearfunc");    
    if(debuglevel>2){
      cout << "Diff_pvmModel::nonlinearfunc"<<endl;
      cout<<paramvalues<<endl;
    }

    float th=paramvalues(1);
    float ph=paramvalues(2);
    float f=paramvalues(3);
    float D=paramvalues(4);
    float S0=paramvalues(5);
    //    cout <<" nlf "<<S0<<endl;
    
    ColumnVector ret(b.Ncols());
    float angtmp;
    for (int i = 1; i <= ret.Nrows(); i++){
      angtmp=cos(ph-beta(i))*sin(alpha(i))*sin(th) + cos(alpha(i))*cos(th);
      angtmp=angtmp*angtmp;
      //      cout <<angtmp<<endl;
      //      cout <<i<<endl;
      ret(i)=S0*(f*exp(-D*b(1,i)*angtmp)+(1-f)*exp(-D*b(1,i)));
    }


    if(debuglevel>2){
      cout<<ret<<endl;
      cout <<"done"<<endl;
    }
    ret.Release();
    return ret; 
  }



void Diff_pvmModel::initialise(const ColumnVector& S){

  Tracer_Plus trace("Diff_pvmModel::initialise");    
  if(debuglevel>2){
    cout << "Diff_pvmModel::initialise"<<endl;
  }


  ColumnVector logS(S.Nrows()),tmp(S.Nrows()),Dvec(7),dir(3);
  SymmetricMatrix tens;   //Basser's Diffusion Tensor;
  DiagonalMatrix Dd;   //eigenvalues
  Matrix Vd;   //eigenvectors
  float mDd,fsquared;
  float th,ph,f,D,S0;
  
  for ( int i = 1; i <= S.Nrows(); i++)
    {
      if(S(i)>0){
	logS(i)=log(S(i));
      }
      else{
	logS(i)=0;
      }
    }
  Dvec = -pinv(Amat)*logS;
  if(  Dvec(7) >  -maxlogfloat ){
    S0=exp(-Dvec(7));
  }
  else{
    S0=S.MaximumAbsoluteValue();
  }

  for ( int i = 1; i <= S.Nrows(); i++)
    {
      if(S0<S.Sum()/S.Nrows()){ S0=S.MaximumAbsoluteValue();  }
      logS(i)=(S(i)/S0)>0.01 ? log(S(i)):log(0.01*S0);
    }
  Dvec = -pinv(Amat)*logS;
  S0=exp(-Dvec(7));
  if(S0<S.Sum()/S.Nrows()){ S0=S.Sum()/S.Nrows();  }
  tens = vec2tens(Dvec);
  EigenValues(tens,Dd,Vd);
  mDd = Dd.Sum()/Dd.Nrows();
  int maxind = Dd(1) > Dd(2) ? 1:2;   //finding maximum eigenvalue
  maxind = Dd(maxind) > Dd(3) ? maxind:3;
  dir << Vd(1,maxind) << Vd(2,maxind) << Vd(3,maxind);
  cart2sph(dir,th,ph);
  th= mod(th,M_PI);
  ph= mod(ph,2*M_PI);
  D = Dd(maxind);

  float numer=(1.5*(Dd(1)-mDd)*(Dd(1)-mDd)+(Dd(2)-mDd)*(Dd(2)-mDd)+(Dd(3)-mDd)*(Dd(3)-mDd));
  float denom=(Dd(1)*Dd(1)+Dd(2)*Dd(2)+Dd(3)*Dd(3));
  if(denom>0) fsquared=numer/denom;
  else fsquared=0;
  if(fsquared>0){f=sqrt(fsquared);}
  else{f=0;}
  if(f>=0.95) f=0.95;
  if(f<=0.001) f=0.001;
  //cout<<"S0 "<<S0<<endl;
  //cout<<"S1 "<<S(1)<<endl;
  getparam(0).setinitvalue(th);
  getparam(1).setinitvalue(ph);
  getparam(2).setinitvalue(f);
  getparam(3).setinitvalue(D);
  getparam(4).setinitvalue(S0);
}


int main(int argc, char *argv[])
{
  try{  

    // Setup logging:
    Log& logger = LogSingleton::getInstance();
    
    // parse command line - will output arguments to logfile
    Diff_pvmOptions& opts = Diff_pvmOptions::getInstance();
    opts.parse_command_line(argc, argv, logger);

    srand(Diff_pvmOptions::getInstance().seed.value());
    
    if(opts.debuglevel.value()==1)
      Tracer_Plus::setrunningstackon();
    
    if(opts.timingon.value())
      Tracer_Plus::settimingon();
    
    // read data

    VolumeSeries data;
    data.read(opts.datafile.value());   
    // data.writeAsFloat(LogSingleton::getInstance().appendDir("data"));
//     cout<<"done"<<endl;
//     return 0;
    int ntpts = data.tsize();
    Matrix bvecs = read_ascii_matrix(opts.bvecsfile.value());
    if(bvecs.Nrows()>3) bvecs=bvecs.t();
    for(int i=1;i<=bvecs.Ncols();i++){
      float tmpsum=sqrt(bvecs(1,i)*bvecs(1,i)+bvecs(2,i)*bvecs(2,i)+bvecs(3,i)*bvecs(3,i));
      if(tmpsum!=0){
	bvecs(1,i)=bvecs(1,i)/tmpsum;
	bvecs(2,i)=bvecs(2,i)/tmpsum;
	bvecs(3,i)=bvecs(3,i)/tmpsum;
      }  
    }

    Matrix bvals = read_ascii_matrix(opts.bvalsfile.value());
    // mask:
    Volume mask;
    mask.read(opts.maskfile.value());
    mask.threshold(1e-16);
    
    // threshold using mask:
    data.setPreThresholdPositions(mask.getPreThresholdPositions());
    data.thresholdSeries();
    
    cout << "ntpts=" << ntpts << endl;
    cout << "nvoxels=" << mask.getVolumeSize() << endl;
    
    Diff_pvmModel model(bvecs,bvals,Diff_pvmOptions::getInstance().debuglevel.value());
    
    LSMCMCManager lsmcmc(Diff_pvmOptions::getInstance(),model,data,mask);
    LSLaplaceManager lslaplace(Diff_pvmOptions::getInstance(),model,data,mask);
    
    
    if(Diff_pvmOptions::getInstance().inference.value()=="mcmc")
      {
	lsmcmc.setup();
	lsmcmc.run();
	element_mod_n(lsmcmc.getsamples(0),M_PI);
	element_mod_n(lsmcmc.getsamples(1),2*M_PI);
	lsmcmc.save();
	
      }
    else
      {
	lslaplace.setup();
	lslaplace.run();
	lslaplace.save();
      }
    
    if(opts.timingon.value())
      Tracer_Plus::dump_times(logger.getDir());

    cout << endl << "Log directory was: " << logger.getDir() << endl;
  }
  catch(Exception& e) 
    {
      cerr << endl << e.what() << endl;
    }
  catch(X_OptionError& e) 
    {
      cerr << endl << e.what() << endl;
    }

  return 0;
}
