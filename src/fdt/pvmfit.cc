/*  Copyright (C) 2009 University of Oxford  */



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
#include "miscmaths/nonlin.h"
#include "newmat.h"
#include "pvmfitOptions.h"
#include "newimage/newimageall.h"
#include "diffmodels.h"

using namespace std;
using namespace NEWMAT;
using namespace MISCMATHS;
using namespace PVMFIT;
using namespace NEWIMAGE;




int main(int argc, char** argv)
{
  //parse command line
  pvmfitOptions& opts = pvmfitOptions::getInstance();
  int success=opts.parse_command_line(argc,argv);
  if(!success) return 1;
   if(opts.verbose.value()){
    cout<<"data file "<<opts.datafile.value()<<endl;
    cout<<"mask file "<<opts.maskfile.value()<<endl;
    cout<<"bvecs     "<<opts.bvecsfile.value()<<endl;
    cout<<"bvals     "<<opts.bvalsfile.value()<<endl;
  }
  
  // Set random seed:
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
  if(bvals.Nrows()>1) bvals=bvals.t();


  volume4D<float> data;
  volume<int> mask;

  if(opts.verbose.value()) cout<<"reading data"<<endl;
  read_volume4D(data,opts.datafile.value());

  if(opts.verbose.value()) cout<<"reading mask"<<endl;
  read_volume(mask,opts.maskfile.value());

  if(opts.verbose.value()) cout<<"ok"<<endl;
  int minx=0;
  int maxx=mask.xsize();
  int miny=0;
  int maxy=mask.ysize();
  int minz=0;
  int maxz=mask.zsize();
  cout<<minx<<" "<<maxx<<" "<<miny<<" "<<maxy<<" "<<minz<<" "<<maxz<<endl;

  if(opts.verbose.value()) cout<<"setting up vols"<<endl;
  volume<float> S0(maxx-minx,maxy-miny,maxz-minz), BIC, f0;
  volume<float> dvol(maxx-minx,maxy-miny,maxz-minz);
  volume<float> tmpvol(maxx-minx,maxy-miny,maxz-minz);
  volume4D<float> tmpvol4D(maxx-minx,maxy-miny,maxz-minz,3);

  vector< volume<float> > fvol,thvol,phvol,k1vol,k2vol, psivol;
  vector< volume4D<float> > dyads;
  vector< volume4D<float> > fanning_vecs; 

  if(opts.verbose.value()) cout<<"copying input properties to output volumes"<<endl;
  copybasicproperties(data[0],S0);
  copybasicproperties(data[0],dvol);
  copybasicproperties(data[0],tmpvol);
  copybasicproperties(data[0],tmpvol4D);

  tmpvol = 0;
  tmpvol4D = 0;
  for(int i=0;i<opts.nfibres.value();i++){
    fvol.push_back(tmpvol);
    thvol.push_back(tmpvol);
    phvol.push_back(tmpvol);
    dyads.push_back(tmpvol4D);
    if (opts.cnonlinear_Fanning.value() || opts.modelnum.value()==4)
      fanning_vecs.push_back(tmpvol4D);
    if (opts.modelnum.value()==3)
      k1vol.push_back(tmpvol);
    if (opts.modelnum.value()==4){
      k1vol.push_back(tmpvol);
      k2vol.push_back(tmpvol);
      psivol.push_back(tmpvol);
    }
  }

  if(opts.verbose.value()) cout<<"zeroing output volumes"<<endl;
  S0=0;dvol=0;
  volume<float> dvol_std;
  if(opts.modelnum.value()==2){
    dvol_std.reinitialize(maxx-minx,maxy-miny,maxz-minz);
    dvol_std=0;
  }
  BIC.reinitialize(maxx-minx,maxy-miny,maxz-minz);
  BIC=0;
  f0.reinitialize(maxx-minx,maxy-miny,maxz-minz);
  f0=0;
  
  

  if(opts.verbose.value()) cout<<"ok"<<endl;

  ColumnVector S(bvals.Ncols());
  if(opts.verbose.value()) cout<<"starting the fits"<<endl;
  for(int k = minz; k < maxz; k++){
    cout<<k<<" slices processed"<<endl;
    for(int j=miny; j < maxy; j++){
      for(int i =minx; i< maxx; i++){
	if(mask(i,j,k)==0)continue;

	for(int t=0;t < data.tsize();t++)
	  S(t+1)=data(i,j,k,t);
	///////////////////////////////////////////////////
	//Fit Ball & sticks with "nfibres" compartments
	///////////////////////////////////////////////////
	if(opts.modelnum.value()==1){
	  if (opts.cnonlinear.value()){  //Use pseudo-constrained optimization
	    PVM_single_c pvm(S,bvecs,bvals,opts.nfibres.value(),opts.saveBIC.value(),opts.use_f0.value());
	    pvm.fit();
	  
	    S0(i-minx,j-miny,k-minz)   = pvm.get_s0();
	    dvol(i-minx,j-miny,k-minz) = pvm.get_d();
	    BIC(i-minx,j-miny,k-minz) = pvm.get_BIC();
	    f0(i-minx,j-miny,k-minz) = pvm.get_f0();
	    for(int f=0;f<opts.nfibres.value();f++){
	      fvol[f](i-minx,j-miny,k-minz)  = pvm.get_f(f+1);
	      thvol[f](i-minx,j-miny,k-minz) = pvm.get_th(f+1);
	      phvol[f](i-minx,j-miny,k-minz) = pvm.get_ph(f+1);
	    }
	  }
	  else if (opts.cnonlinear_Fanning.value()){  //Use pseudo-constrained optimization and return fanning angle estimates using the Hessian of the cost function
	    PVM_single_c pvm(S,bvecs,bvals,opts.nfibres.value(),opts.saveBIC.value(),opts.use_f0.value(),true);
	    pvm.fit();
	  
	    S0(i-minx,j-miny,k-minz)   = pvm.get_s0();
	    dvol(i-minx,j-miny,k-minz) = pvm.get_d();
	    BIC(i-minx,j-miny,k-minz) = pvm.get_BIC();
	    f0(i-minx,j-miny,k-minz) = pvm.get_f0();
	    for(int f=0;f<opts.nfibres.value();f++){
	      fvol[f](i-minx,j-miny,k-minz)  = pvm.get_f(f+1);
	      thvol[f](i-minx,j-miny,k-minz) = pvm.get_th(f+1);
	      phvol[f](i-minx,j-miny,k-minz) = pvm.get_ph(f+1);
	      ColumnVector tmp_vec= pvm.get_invHes_e1(f+1);
	      fanning_vecs[f](i-minx,j-miny,k-minz,0)=tmp_vec(1);
	      fanning_vecs[f](i-minx,j-miny,k-minz,1)=tmp_vec(2);
	      fanning_vecs[f](i-minx,j-miny,k-minz,2)=tmp_vec(3);
	    }
	  }
	  else{  //Use original optimization
	    PVM_single pvm(S,bvecs,bvals,opts.nfibres.value(), opts.use_f0.value());
	    pvm.fit();
	  
	    S0(i-minx,j-miny,k-minz)   = pvm.get_s0();
	    dvol(i-minx,j-miny,k-minz) = pvm.get_d();
	    f0(i-minx,j-miny,k-minz) = pvm.get_f0();
	    for(int f=0;f<opts.nfibres.value();f++){
	      fvol[f](i-minx,j-miny,k-minz)  = pvm.get_f(f+1);
	      thvol[f](i-minx,j-miny,k-minz) = pvm.get_th(f+1);
	      phvol[f](i-minx,j-miny,k-minz) = pvm.get_ph(f+1);
	    }
	  }
	}
	/////////////////////////////////////////////////////////
	//Fit Ball & sticks (model 2) with "nfibres" compartments
	/////////////////////////////////////////////////////////
	else if (opts.modelnum.value()==2){ 
	  PVM_multi pvm(S,bvecs,bvals,opts.nfibres.value());
	  pvm.fit();

	  S0(i-minx,j-miny,k-minz)   = pvm.get_s0();
	  dvol(i-minx,j-miny,k-minz) = pvm.get_d();
	  dvol_std(i-minx,j-miny,k-minz) = pvm.get_d_std();
	  for(int f=0;f<opts.nfibres.value();f++){
	    fvol[f](i-minx,j-miny,k-minz)  = pvm.get_f(f+1);
	    thvol[f](i-minx,j-miny,k-minz) = pvm.get_th(f+1);
	    phvol[f](i-minx,j-miny,k-minz) = pvm.get_ph(f+1);
	  }
	}

	///////////////////////////////////////////////////
	//Fit Ball & Watsons with "nfibres" compartments
	///////////////////////////////////////////////////
	else if (opts.modelnum.value()==3){ 
	  cout<<i<<" "<<j<<" "<<k<<endl;
	  if (!opts.all.value()){
	    PVM_Ball_Watsons pvm(S,bvecs,bvals,opts.nfibres.value(),opts.saveBIC.value(),opts.use_f0.value(),opts.gridsearch.value());
	    pvm.fit();
	    
	    S0(i-minx,j-miny,k-minz)   = pvm.get_s0();
	    dvol(i-minx,j-miny,k-minz) = pvm.get_d();
	    BIC(i-minx,j-miny,k-minz) = pvm.get_BIC();
	    for (int f=0;f<opts.nfibres.value();f++){
	      fvol[f](i-minx,j-miny,k-minz)  = pvm.get_f(f+1);
	      thvol[f](i-minx,j-miny,k-minz) = pvm.get_th(f+1);
	      phvol[f](i-minx,j-miny,k-minz) = pvm.get_ph(f+1);
	      k1vol[f](i-minx,j-miny,k-minz) = pvm.get_k(f+1);
	    } 
	  }
	  else{  //Fit all Ball & Watsons with up to "nfibres" compartments and choose the best using BIC
	    float bestBIC=1.0e20;
	    for (int n=1; n<=opts.nfibres.value(); n++){
	      PVM_Ball_Watsons pvmn(S,bvecs,bvals,n,true,opts.use_f0.value(),opts.gridsearch.value());
	      pvmn.fit();	    
	      if (pvmn.get_BIC()<bestBIC){ //Keep the model with the smallest BIC
		bestBIC=pvmn.get_BIC();
		S0(i-minx,j-miny,k-minz)   = pvmn.get_s0();
		dvol(i-minx,j-miny,k-minz) = pvmn.get_d();
		BIC(i-minx,j-miny,k-minz) = pvmn.get_BIC();
		for (int f=0;f<n;f++){   //compartments are not sorted here! 
		  fvol[f](i-minx,j-miny,k-minz)  = pvmn.get_f(f+1);
		  thvol[f](i-minx,j-miny,k-minz) = pvmn.get_th(f+1);
		  phvol[f](i-minx,j-miny,k-minz) = pvmn.get_ph(f+1);
		  k1vol[f](i-minx,j-miny,k-minz) = pvmn.get_k(f+1);
		}
	      }
	    }
	  }
	}

	///////////////////////////////////////////////////
	//Fit Ball & Binghams with "nfibres" compartments
	///////////////////////////////////////////////////
	else if (opts.modelnum.value()==4){ 
	  cout<<i<<" "<<j<<" "<<k<<endl;
	  
	  if (!opts.all.value()){
	    PVM_Ball_Binghams pvm(S,bvecs,bvals,opts.nfibres.value(),opts.saveBIC.value(),opts.use_f0.value(),opts.gridsearch.value());
	    pvm.fit();
	    S0(i-minx,j-miny,k-minz)   = pvm.get_s0();
	    dvol(i-minx,j-miny,k-minz) = pvm.get_d();
	    BIC(i-minx,j-miny,k-minz) = pvm.get_BIC();
	    f0(i-minx,j-miny,k-minz) = pvm.get_f0();
	    for (int f=0;f<opts.nfibres.value();f++){
	      fvol[f](i-minx,j-miny,k-minz)  = pvm.get_f(f+1);
	      thvol[f](i-minx,j-miny,k-minz) = pvm.get_th(f+1);
	      phvol[f](i-minx,j-miny,k-minz) = pvm.get_ph(f+1);
	      k1vol[f](i-minx,j-miny,k-minz) = pvm.get_k1(f+1);
	      k2vol[f](i-minx,j-miny,k-minz) = pvm.get_k2(f+1);
	      psivol[f](i-minx,j-miny,k-minz) = pvm.get_psi(f+1);
	      ColumnVector tmp_vec= pvm.get_fanning_vector(f+1);
	      fanning_vecs[f](i-minx,j-miny,k-minz,0)=tmp_vec(1);
	      fanning_vecs[f](i-minx,j-miny,k-minz,1)=tmp_vec(2);
	      fanning_vecs[f](i-minx,j-miny,k-minz,2)=tmp_vec(3);
	    } 
	  }
	  else{	//Fit all Ball & Binghams with up to "nfibres" compartments and choose the best using BIC
	    float bestBIC=1.0e20;
	    for (int n=1; n<=opts.nfibres.value(); n++){
	      PVM_Ball_Binghams pvmn(S,bvecs,bvals,n,true,opts.use_f0.value(),opts.gridsearch.value());
	      pvmn.fit();
	      if (pvmn.get_BIC()<bestBIC){ //Keep the model with the smallest BIC
		bestBIC=pvmn.get_BIC();
		S0(i-minx,j-miny,k-minz)   = pvmn.get_s0();
		dvol(i-minx,j-miny,k-minz) = pvmn.get_d();
		BIC(i-minx,j-miny,k-minz) = pvmn.get_BIC();
		f0(i-minx,j-miny,k-minz) = pvmn.get_f0();
		for (int f=0;f<n;f++){   //compartments are not sorted here! 
		  fvol[f](i-minx,j-miny,k-minz)  = pvmn.get_f(f+1);
		  thvol[f](i-minx,j-miny,k-minz) = pvmn.get_th(f+1);
		  phvol[f](i-minx,j-miny,k-minz) = pvmn.get_ph(f+1);
		  k1vol[f](i-minx,j-miny,k-minz) = pvmn.get_k1(f+1);
		  k2vol[f](i-minx,j-miny,k-minz) = pvmn.get_k2(f+1);
		  psivol[f](i-minx,j-miny,k-minz) = pvmn.get_psi(f+1);
		  ColumnVector tmp_vec= pvmn.get_fanning_vector(f+1);
		  fanning_vecs[f](i-minx,j-miny,k-minz,0)=tmp_vec(1);
		  fanning_vecs[f](i-minx,j-miny,k-minz,1)=tmp_vec(2);
		  fanning_vecs[f](i-minx,j-miny,k-minz,2)=tmp_vec(3);
		}
	      }
	    }
	  }
	}

	for(int f=0;f<opts.nfibres.value();f++)
	  if (fvol[f](i-minx,j-miny,k-minz)!=0){
	    dyads[f](i-minx,j-miny,k-minz,0) = sin(thvol[f](i-minx,j-miny,k-minz))*cos(phvol[f](i-minx,j-miny,k-minz));
	    dyads[f](i-minx,j-miny,k-minz,1) = sin(thvol[f](i-minx,j-miny,k-minz))*sin(phvol[f](i-minx,j-miny,k-minz));
	    dyads[f](i-minx,j-miny,k-minz,2) = cos(thvol[f](i-minx,j-miny,k-minz));
	  }
      }
    }
  }

  if(opts.verbose.value())
    cout << "saving results" << endl;

  S0.setDisplayMaximumMinimum(S0.max(),0);
  save_volume(S0,opts.ofile.value()+"_S0");

  dvol.setDisplayMaximumMinimum(dvol.max(),0);
  save_volume(dvol,opts.ofile.value()+"_d");

  if(opts.modelnum.value()==2){
    dvol_std.setDisplayMaximumMinimum(dvol_std.max(),0);
    save_volume(dvol_std,opts.ofile.value()+"_d_std");
  }

  if(opts.saveBIC.value()){
    BIC.setDisplayMaximumMinimum(BIC.max(),BIC.min());
    save_volume(BIC,opts.ofile.value()+"_BIC");
  }

  if(opts.use_f0.value()){
    f0.setDisplayMaximumMinimum(1,0);
    save_volume(f0,opts.ofile.value()+"_f0");
  }

  for(int f=1;f<=opts.nfibres.value();f++){
    fvol[f-1].setDisplayMaximumMinimum(1,0);
    save_volume(fvol[f-1],opts.ofile.value()+"_f"+num2str(f));
    thvol[f-1].setDisplayMaximumMinimum(thvol[f-1].max(),thvol[f-1].min());
    save_volume(thvol[f-1],opts.ofile.value()+"_th"+num2str(f));
    phvol[f-1].setDisplayMaximumMinimum(phvol[f-1].max(),phvol[f-1].min());
    save_volume(phvol[f-1],opts.ofile.value()+"_ph"+num2str(f));
    dyads[f-1].setDisplayMaximumMinimum(1,-1);
    save_volume4D(dyads[f-1],opts.ofile.value()+"_dyads"+num2str(f));
    if (opts.cnonlinear_Fanning.value()){
      fanning_vecs[f-1].setDisplayMaximumMinimum(1,-1);
      save_volume4D(fanning_vecs[f-1],opts.ofile.value()+"_fans"+num2str(f));
    }
    if (opts.modelnum.value()==3){
      k1vol[f-1].setDisplayMaximumMinimum(k1vol[f-1].max(),k1vol[f-1].min());
      save_volume(k1vol[f-1],opts.ofile.value()+"_k1_"+num2str(f));
    }
    if (opts.modelnum.value()==4){
      k1vol[f-1].setDisplayMaximumMinimum(k1vol[f-1].max(),k1vol[f-1].min());
      save_volume(k1vol[f-1],opts.ofile.value()+"_k1_"+num2str(f));
      k2vol[f-1].setDisplayMaximumMinimum(k2vol[f-1].max(),k2vol[f-1].min());
      save_volume(k2vol[f-1],opts.ofile.value()+"_k2_"+num2str(f));
      psivol[f-1].setDisplayMaximumMinimum(psivol[f-1].max(),psivol[f-1].min());
      save_volume(psivol[f-1],opts.ofile.value()+"_psi_"+num2str(f));
      fanning_vecs[f-1].setDisplayMaximumMinimum(1,-1);
      save_volume4D(fanning_vecs[f-1],opts.ofile.value()+"_fans"+num2str(f));
    }
  }
 return 0;
}













