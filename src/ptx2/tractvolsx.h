/*  tractvolsx.h

    Tim Behrens, Saad Jbabdi, FMRIB Image Analysis Group

    Copyright (C) 2004 University of Oxford  */

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

#ifndef __TRACTVOLSX_H_
#define __TRACTVOLSX_H_

/////////////////////////////////////////////////////////
//         Class TractVolsx                             //
/////////////////////////////////////////////////////////

#include "newimage/newimageall.h"
#include <iostream>
#include "stdlib.h"
#include "probtrackxOptions.h"
#include "utils/tracer_plus.h"
using namespace std;
using namespace NEWIMAGE;
using namespace TRACT;
using namespace Utilities;

namespace TRACTVOLSX{
  class Tractvolsx
    {
    private:
      probtrackxOptions& opts;
      Log&               logger;

      vector<Matrix> thsamples;
      vector<Matrix> phsamples;
      vector<Matrix> fsamples;

      volume<int>    lut_vol2mat;

      int            nfibres;
      int            nsamples;

      bool           init_sample;
      int            fibst;
      bool           usef;
      
      volume<int>    locfibchoice;

    public:
      //constructors::
      Tractvolsx(const bool& usefin=false):opts(probtrackxOptions::getInstance()),
					   logger(LogSingleton::getInstance()),
					   init_sample(true),fibst(0),usef(usefin){}
      ~Tractvolsx(){}
      int get_nfibres()const{return nfibres;}
      int get_nsamples()const{return nsamples;}
      
      void reset(const int& fibst_in){
	init_sample=true;
	fibst=fibst_in;
      }

      int sample_fibre(int col,int samp,const ColumnVector& dir){
	float th,ph;ColumnVector x(3);
	vector<int> fibvec;
	for(int fib=0;fib<nfibres;fib++){	    
	  float ft=fsamples[fib](samp,col);
	  if(ft>opts.fibthresh.value()){
	    th=thsamples[fib](samp,col);
	    ph=phsamples[fib](samp,col);
	    x<<sin(th)*cos(ph)<<sin(th)*sin(ph)<<cos(th);
	    if( fabs( x(1)*dir(1)+x(2)*dir(2)+x(3)*dir(3) ) > 0.766 ){ //hard-coded 40 deg threshold
	      fibvec.push_back(fib);
	    }
	  }
	}
	if(fibvec.size()==0){
	  return 0;
	}
	else{
	  double rtmp=(rand()/(double(RAND_MAX)+1)) * fibvec.size();
	  return (fibvec[ (int)floor(rtmp) ]);
	}
      }

      int sample_fibre(int col,int samp,const int& mode=2){
	if(mode==0){
	  return 0;
	}
	if(mode==3){//sample all
	  double rtmp=(rand()/(double(RAND_MAX)+1)) * nfibres;
	  return int(floor(rtmp));
	}
	else{
	  if(mode==1){//sample all>thresh
	    vector<int> fibvec;
	    for(int fib=0;fib<nfibres;fib++){	    
	      float ft=fsamples[fib](samp,col);
	      if(ft>opts.fibthresh.value()){
		fibvec.push_back(fib);
	      }
	    }
	    if(fibvec.size()==0){
	      return 0;
	    }
	    else{
	      double rtmp=(rand()/(double(RAND_MAX)+1)) * fibvec.size();
	      return (fibvec[ (int)floor(rtmp) ]);
	    }
	  }
	  else if(mode==2){//sample all>thresh in proportion of f (default)
	    float fsumtmp=0;
	    for(int fib=0;fib<nfibres;fib++){	    
	      float ft=fsamples[fib](samp,col);
	      if(ft>opts.fibthresh.value()){
		fsumtmp+=ft;  //count total weight of f in this voxel. 
	      }
	    } 
	    if(fsumtmp==0){
	      return(0);
	    }
	    else{
	      float ft,fsumtmp2=0;
	      float rtmp=fsumtmp * (float)rand()/float(RAND_MAX);	      
	      for(int fib=0;fib<nfibres;fib++){
		ft=fsamples[fib](samp,col);
		if(ft>opts.fibthresh.value())
		  fsumtmp2 += ft;
		if(rtmp<=fsumtmp2){
		  return(fib); 
		}
	      }
	    }
	  }
	  else{
	    cerr<<"TRACTVOLSX::sample_fibre:Error - unknown mode = "<<mode<<endl;
	    exit(1);
	  }
	}
	return 0;
      }

      int sample_ang_prob(const vector<float>& probs){
	float sum=0;ColumnVector cumsum(probs.size());cumsum=0;
	int ind=0;
	for (unsigned int i=0;i<probs.size();i++){
	  sum += probs[i];
	  cumsum(i+1)=sum;
	}
	float U=rand()/float(RAND_MAX);
	U *= sum;
	for(unsigned int k=1;k<=probs.size();k++){
	  if(U<cumsum(k)){
	    ind=k-1;
	    break;
	  }
	}
	return ind;
      }

      //Initialise
      void initialise(const string& basename,const volume<float>& mask){
	volume4D<float> tmpvol;
	Matrix          tmpmat;		

	cout<<"Load bedpostx samples"<<endl;
	if(fsl_imageexists(basename+"_thsamples")){
	  cout<<"1"<<endl;
	  read_volume4D(tmpvol,basename+"_thsamples");
	  tmpmat=tmpvol.matrix(mask);
	  thsamples.push_back(tmpmat);
	  cout<<"2"<<endl;
	  read_volume4D(tmpvol,basename+"_phsamples");
	  tmpmat=tmpvol.matrix(mask);
	  phsamples.push_back(tmpmat);
	  cout<<"3"<<endl;
	  read_volume4D(tmpvol,basename+"_fsamples");
	  tmpmat=tmpvol.matrix(mask);
	  fsamples.push_back(tmpmat);

	  lut_vol2mat = tmpvol.vol2matrixkey(mask);
	  nsamples    = tmpmat.Nrows();
	  nfibres     = 1;
	}
	else{
	  int fib=1;
	  bool fib_existed=true;
	  while(fib_existed){
	    if(fsl_imageexists(basename+"_th"+num2str(fib)+"samples")){
	      cout<<fib<<"_1"<<endl;
	      read_volume4D(tmpvol,basename+"_th"+num2str(fib)+"samples");
	      tmpmat=tmpvol.matrix(mask);
	      thsamples.push_back(tmpmat);
	      cout<<fib<<"_2"<<endl;
	      read_volume4D(tmpvol,basename+"_ph"+num2str(fib)+"samples");
	      tmpmat=tmpvol.matrix(mask);
	      phsamples.push_back(tmpmat);
	      cout<<fib<<"_3"<<endl;
	      read_volume4D(tmpvol,basename+"_f"+num2str(fib)+"samples");
	      tmpmat=tmpvol.matrix(mask);
	      fsamples.push_back(tmpmat);
	      fib++;
	    }
	    else{
	      fib_existed=false;
	    }
	  }
	  if(fib==1){
	      cerr<<"Could not find samples to load. Exit without doing anything"<<endl;
	      exit(1);
	  }
	  lut_vol2mat = tmpvol.vol2matrixkey(mask);
	  nsamples = thsamples[0].Nrows();
	  nfibres  = (int)thsamples.size();
	}
	copybasicproperties(mask,lut_vol2mat);

	cout<<endl;
	cout<<"nfibres  : "<<nfibres<<endl;
	cout<<"nsamples : "<<nsamples<<endl;
	cout<<endl;
	cout<<"Done loading samples."<<endl;

	if(opts.locfibchoice.value()!=""){
	  read_volume(locfibchoice,opts.locfibchoice.value());	  
	}
      }
      
      
      ColumnVector sample(const float& x,const float& y,const float&z,
			  const float& r_x,const float& r_y,const float& r_z,
			  float& prefer_x,float& prefer_y,float& prefer_z,
			  const int& sample_fib,int& sampled_fib,
			  int& newx,int& newy,int& newz){

	//Tracer_Plus tr("sample");
	////////Probabilistic interpolation
	int cx =(int) ceil(x),fx=(int) floor(x);
	int cy =(int) ceil(y),fy=(int) floor(y);
	int cz =(int) ceil(z),fz=(int) floor(z);
	
	float pcx = (cx==fx)?1:(x-fx)/(cx-fx);
	float pcy = (cy==fy)?1:(y-fy)/(cy-fy);
	float pcz = (cz==fz)?1:(z-fz)/(cz-fz);
	
	newx = ((float)rand()/(float)RAND_MAX)>pcx?fx:cx;
	newy = ((float)rand()/(float)RAND_MAX)>pcy?fy:cy;
	newz = ((float)rand()/(float)RAND_MAX)>pcz?fz:cz;
	////////////////////////////////////	

	ColumnVector th_ph_f(3);	

	int col = lut_vol2mat(newx,newy,newz);
	if(col==0){//outside brain mask
	  th_ph_f=0;
	  return th_ph_f;
	}

	int samp=(int)MISCMATHS::round((float)rand()/float(RAND_MAX)*(float)(nsamples-1))+1;

	float theta=0,phi=0;
	float dotmax=0,dottmp=0;
	int fibind=0;
	if(nfibres>1){//more than 1 fibre
	  if(init_sample){//go for the specified fibre on the first jump or generate at random
	    if(!opts.fibst.set())
	      fibst=sample_fibre(col,samp,opts.randfib.value());

	    theta=thsamples[fibst](samp,col);
	    phi=phsamples[fibst](samp,col);
	    init_sample=false;
	  }
	  else{
	    if(sample_fib>0){ // pick specified fibre
	      fibind=sample_fibre(col,samp,sample_fib);	      
	      theta=thsamples[fibind](samp,col);
	      phi=phsamples[fibind](samp,col);
	    }
	    else{ 
	      if((fabs(prefer_x)+fabs(prefer_y)+fabs(prefer_z))==0){
		prefer_x=r_x;prefer_y=r_y;prefer_z=r_z;
	      }
	      int locrule=0;
	      if(opts.locfibchoice.value()!=""){locrule=locfibchoice(newx,newy,newz);}
	      if(locrule==1){
		fibind=sample_fibre(col,samp,1);
		theta=thsamples[fibind](samp,col);
		phi=phsamples[fibind](samp,col);
	      }
	      else if(locrule==2){ // like locrule=1 but with angle threshold
		ColumnVector dir(3);dir<<r_x<<r_y<<r_z;
		fibind=sample_fibre(col,samp,dir);
		theta=thsamples[fibind](samp,col);
		phi=phsamples[fibind](samp,col);
	      }
	      else if (locrule==3) { 
		fibind=sample_fibre(col,samp,2);
		theta=thsamples[fibind](samp,col);
		phi=phsamples[fibind](samp,col);
	      }
	      else{ // pick closest direction
		for(int fib=0;fib<nfibres;fib++){
		  if(fsamples[fib](samp,col)>opts.fibthresh.value()){
		    float phtmp=phsamples[fib](samp,col);
		    float thtmp=thsamples[fib](samp,col);
		    dottmp=fabs(sin(thtmp)*(cos(phtmp)*prefer_x + sin(phtmp)*prefer_y) + cos(thtmp)*prefer_z);
		    if(dottmp>dotmax){
		      dotmax=dottmp;
		      theta=thtmp;
		      phi=phtmp;
		      fibind=fib;
		    }
		  }
		}
		if(dotmax==0){
		  theta=thsamples[0](samp,col);
		  phi=phsamples[0](samp,col);
		  fibind=0;
		}
	      }
	    }
	  }
	}
	else{
	  theta=thsamples[0](samp,col);
	  phi=phsamples[0](samp,col);
	}
	
	float f;	
	if(usef){
	  f = fsamples[fibind](samp,col);
	}
	else{
	  f=1;
	}

	sampled_fib = fibind+1;

	th_ph_f(1)=theta;
	th_ph_f(2)=phi;
	th_ph_f(3)=f;
	return th_ph_f;
      }

      ColumnVector dimensions() const{
	ColumnVector dims(3);
	dims << lut_vol2mat.xdim() <<lut_vol2mat.ydim() << lut_vol2mat.zdim();
	return dims;
      }
    };
}

#endif



