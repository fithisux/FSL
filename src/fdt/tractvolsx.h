/*  tractvolsx.h

    Tim Behrens, FMRIB Image Analysis Group

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
using namespace std;
using namespace NEWIMAGE;
using namespace TRACT;

namespace TRACTVOLSX{
  class Tractvolsx
    {
    private:
      probtrackxOptions& opts;
      vector<volume4D<float>* > thsamples;
      vector<volume4D<float>* > phsamples;
      vector<volume4D<float>* > fsamples;
      bool init_sample;
      int fibst;
      bool usef;
      
    public:
      //constructors::
      Tractvolsx(const bool& usefin=false):opts(probtrackxOptions::getInstance()),init_sample(true),fibst(0),usef(usefin){}
      Tractvolsx():opts(probtrackxOptions::getInstance()){}
      ~Tractvolsx(){
	for(unsigned int m=0;m<thsamples.size();m++)
	  delete thsamples[m]; //ask flitney, do you just delete the ptr??
	for(unsigned int m=0;m<phsamples.size();m++)
	  delete phsamples[m];
	for(unsigned int m=0;m<fsamples.size();m++)
	    delete fsamples[m];
      }
      inline int nfibres()const{return (int)thsamples.size();}
      
      void reset(const int& fibst_in){
	init_sample=true;
	fibst=fibst_in;
      }
      //Initialise
      void initialise(const string& basename){
	

	if(fsl_imageexists(basename+"_thsamples")){
	  volume4D<float> *tmpthptr= new volume4D<float>;
	  volume4D<float> *tmpphptr= new volume4D<float>;
	  volume4D<float> *tmpfptr= new volume4D<float>;
	  cout<<"1"<<endl;
	  read_volume4D(*tmpthptr,basename+"_thsamples");
	  cout<<"2"<<endl;
	  thsamples.push_back(tmpthptr);
	  cout<<"3"<<endl;
	  read_volume4D(*tmpphptr,basename+"_phsamples");
	  cout<<"4"<<endl;
	  phsamples.push_back(tmpphptr);
	  cout<<"5"<<endl;
	  read_volume4D(*tmpfptr,basename+"_fsamples");
	  fsamples.push_back(tmpfptr);
	  cout<<"6"<<endl;
	}
	else{
	  int fib=1;
	  bool fib_existed=true;
	  while(fib_existed){
	    if(fsl_imageexists(basename+"_th"+num2str(fib)+"samples")){
	      volume4D<float> *tmpthptr= new volume4D<float>;
	      volume4D<float> *tmpphptr= new volume4D<float>;
	      volume4D<float> *tmpfptr= new volume4D<float>;
	      cout<<fib<<"_1"<<endl;
	      read_volume4D(*tmpthptr,basename+"_th"+num2str(fib)+"samples");
	      thsamples.push_back(tmpthptr);
	      cout<<fib<<"_2"<<endl;
	      read_volume4D(*tmpphptr,basename+"_ph"+num2str(fib)+"samples");
	      phsamples.push_back(tmpphptr);
	      cout<<fib<<"_3"<<endl;
	      read_volume4D(*tmpfptr,basename+"_f"+num2str(fib)+"samples");
	      fsamples.push_back(tmpfptr);
	      fib++;
	    }
	    else{
	      fib_existed=false;
	    }
	  }
	  
	}
	cout<<"7"<<endl;
      }
      
      
      ColumnVector sample(const float& x,const float& y,const float&z,const float& r_x,const float& r_y,const float& r_z,
			  float& prefer_x,float& prefer_y,float& prefer_z){

	////////Probabilistic interpolation
	int cx =(int) ceil(x),fx=(int) floor(x);
	int cy =(int) ceil(y),fy=(int) floor(y);
	int cz =(int) ceil(z),fz=(int) floor(z);
	
	//cerr<<x<<" "<<y<<" "<<z<<" "<<cx<<" "<<cy<<" "<<cz<<" "<<fx<<" "<<fy<<" "<<fz<<endl;
	float pcx,pcy,pcz;
	if(cx==fx)
	  pcx=1;
	else
	  pcx=(x-fx)/(cx-fx);
	
	if(cy==fy)
	  pcy=1;
	else
	  pcy=(y-fy)/(cy-fy);
	
	if(cz==fz)
	  pcz=1;
	else
	  pcz=(z-fz)/(cz-fz);
	
	///////new xyz values from probabilistic interpolation
	int newx,newy,newz; 
	float tmp=(float)rand()/float(RAND_MAX);
	if(tmp>pcx)
	  newx=fx;
	else
	  newx=cx;
	
	tmp=(float)rand()/float(RAND_MAX);
	if(tmp>pcy)
	  newy=fy;
	else
	  newy=cy;
	
	tmp=(float)rand()/float(RAND_MAX);
	if(tmp>pcz)
	  newz=fz;
	else
	  newz=cz;
 
	ColumnVector th_ph_f(3);	
	float samp=(float)rand()/float(RAND_MAX);
	samp=MISCMATHS::round(samp*((*thsamples[0]).tsize()-1));
	float theta=0,phi=0;
	float dotmax=0,dottmp=0;
	int fibind=0;
	if(thsamples.size()>1){//more than 1 fibre
	  if(init_sample){//go for the specified fibre on the first jump or generate at random
	    if(opts.randfib.value()==1){//this generates startfib at random (except for fibres where f<fibthresh)
	      vector<int> fibvec;
	      for(unsigned int fib=0;fib<thsamples.size();fib++){	    
		float ft=(*fsamples[fib])(int(newx),int(newy),int(newz),int(samp));
		if(ft>opts.fibthresh.value()){
		  fibvec.push_back(fib);
		}
	      }
	      
	      if(fibvec.size()==0){
		fibst=0;
	      }
	      else{
		float rtmp=(float)rand()/float(RAND_MAX) * float(fibvec.size()-1);
		fibst = fibvec[ (int)MISCMATHS::round(rtmp) ];	      
	      }
	      
	    }
	    else if(opts.randfib.value()==2){ //this generates startfib with probability proportional to f (except for fibres where f<fibthresh). 
	      //this chooses at random but in proportion to fsamples. 
	      float fsumtmp=0;
	      for(unsigned int fib=0;fib<thsamples.size();fib++){	    
		float ft=(*fsamples[fib])(int(newx),int(newy),int(newz),int(samp));
		if(ft>opts.fibthresh.value()){
		  fsumtmp+=ft;  //count total weight of f in this voxel. 
		}
	      }
	      
	      if(fsumtmp==0){
		fibst=0;
	      }
	      else{
		float ft,fsumtmp2=0;
		float rtmp=fsumtmp * (float)rand()/float(RAND_MAX);
		
		for(unsigned int fib=0;fib<thsamples.size();fib++){
		  ft=(*fsamples[fib])(int(newx),int(newy),int(newz),int(samp));
		  if(ft>opts.fibthresh.value())
		    fsumtmp2 += ft;
		  if(rtmp<=fsumtmp2){
		    fibst=(int)fib;
		    break;
		  }
		}		
	      }
	    }  

	    theta=(*thsamples[fibst])(int(newx),int(newy),int(newz),int(samp));
	    phi=(*phsamples[fibst])(int(newx),int(newy),int(newz),int(samp));
	    init_sample=false;
	  }
	  else{
	    if((fabs(prefer_x)+fabs(prefer_y)+fabs(prefer_z))==0){
	      prefer_x=r_x;prefer_y=r_y;prefer_z=r_z;
	    }
	    for(unsigned int fib=0;fib<thsamples.size();fib++){
	      if((*fsamples[fib])(int(newx),int(newy),int(newz),int(samp))>opts.fibthresh.value()){
		float phtmp=(*phsamples[fib])(int(newx),int(newy),int(newz),int(samp));
		float thtmp=(*thsamples[fib])(int(newx),int(newy),int(newz),int(samp));
		dottmp=fabs(sin(thtmp)*cos(phtmp)*prefer_x + sin(thtmp)*sin(phtmp)*prefer_y + cos(thtmp)*prefer_z);
		if(dottmp>dotmax){
		  dotmax=dottmp;
		  theta=thtmp;
		  phi=phtmp;
		  fibind=fib;
		}
		
	      }
	      
	    }
	    if(dotmax==0){
	      theta=(*thsamples[0])(int(newx),int(newy),int(newz),int(samp));
	      phi=(*phsamples[0])(int(newx),int(newy),int(newz),int(samp));
	    }
	  }
	}
	else{
	  theta=(*thsamples[0])(int(newx),int(newy),int(newz),int(samp));
	  phi=(*phsamples[0])(int(newx),int(newy),int(newz),int(samp));
	}

	
	float f;
	
	if(usef){
	  f = (*fsamples[fibind])(int(newx),int(newy),int(newz),int(samp));
	}
	else
	  f=1;
	
	th_ph_f(1)=theta;
	th_ph_f(2)=phi;
	th_ph_f(3)=f;
	return th_ph_f;
      }

      ColumnVector dimensions() const{
	ColumnVector dims(3);
	dims << (*thsamples[0]).xdim() <<(*thsamples[0]).ydim() << (*thsamples[0]).zdim();
	return dims;
      }
    };
}

#endif



