/*  MELODIC - Multivariate exploratory linear optimized decomposition into 
              independent components
    
    melgmix.cc - Gaussian Mixture Model

    Christian F. Beckmann, FMRIB Analysis Group
    
    Copyright (C) 1999-2013 University of Oxford */

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

#include "newimage/newimageall.h"
#include "melgmix.h"
#include "utils/log.h"
#include "miscmaths/miscprob.h"
#include <time.h>
#include "libvis/miscplot.h"
#include "libvis/miscpic.h"

using namespace Utilities;
using namespace NEWIMAGE;
using namespace MISCPLOT;
using namespace MISCPIC;

string float2str(float f,int width, int prec, bool scientif){
 	ostringstream os;
  int redw = int(std::abs(std::log10(std::abs(f))))+1;
  if(width>0)
    os.width(width);
  if(scientif)
    os.setf(ios::scientific);
  os.precision(redw+std::abs(prec));
  os.setf(ios::internal, ios::adjustfield);
  os << f;
  return os.str();
} 
  
namespace Melodic{
  
  void MelGMix::setup(const RowVector& dat, 
		const string dirname,
		int cnum, volume<float> themask, 
		volume<float> themean, 
		int num_mix, float eps, bool fixit){
    	cnumber = cnum;
    	Mask = themask;
    	Mean = themean;
    	prefix = string("IC_")+num2str(cnum);
    	
    	fitted = false;
    	nummix = num_mix;
    	numdata = dat.Ncols();
    
    	//normalise data 
    	datamean = mean(dat,2).AsScalar();
    	datastdev= stdev(dat,2).AsScalar();
    	data=(dat - datamean)/datastdev;

			dbgmsg(" mapdat; mean: " << datamean << " std: " <<datastdev << endl);

    	props=zeros(1,nummix);
    	vars=zeros(1,nummix);
    	means=zeros(1,nummix);
    	Params=zeros(1,nummix);
    	logprobY = 1.0;

    	props = std::pow(float(nummix),float(-1.0));

    	Matrix tmp1;
    	tmp1 = data * data.t() / numdata;
    	vars = tmp1.AsScalar();

    	float Dmin, Dmax, IntSize;
    	Dmin =  min(data).AsScalar(); Dmax = max(data).AsScalar();
    	IntSize = Dmax / nummix;

    	means(1)=mean(data,2).AsScalar(); 
    	for (int ctr=2; ctr <= means.Ncols(); ctr++){
      	means(ctr) =  Dmax - (ctr - 1.5) * IntSize; 
    	}
    	means(2)=means(1)+2*sqrt(vars(1));
    	if(nummix>2)
        means(3)=means(1)-2*sqrt(vars(1));

    	epsilon = eps;
    	if(epsilon >=0 && epsilon < 0.0000001)
      	{epsilon = std::log(float(dat.Ncols()))/1000 ;}
    	fixdim = fixit;
			dbgmsg(" parameters; " << means << " : " << vars << " : " << props << endl);
		}

  Matrix MelGMix::threshold(const RowVector& dat, string levels){ 
    Matrix Res;
    Res = 1.0;
    string tmpstr;
    tmpstr = levels;
    //cerr << " Levels : " << levels << endl << endl;
    Matrix levls(1,4);
    levls = 0;
    Matrix fpr;
    Matrix fdr;
    Matrix nht;

    char *p;
    char t[1024];
    const char *discard = ", [];{(})abceghijklmoqstuvwxyzABCEGHIJKLMNOQSTUVWXYZ~!@#$%^&*_-=+|\':></?";

    strcpy(t, tmpstr.c_str());
    p=strtok(t,discard);
    while(p){
      Matrix New(1,1);
      New(1,1) = atof(p);
      if(strchr(p,'d')){
				levls(1,3)++;
				if(fdr.Storage()>0){
	  			fdr = fdr | New;
				}else{
	  			fdr = New;
				}
      }else{
				if(strchr(p,'p')){
	  			levls(1,2)++;
	  			if(fpr.Storage()>0){
	    			fpr = fpr | New;
	  			}else{
	    			fpr = New;
	  			}
				}else{
	  			if(strchr(p,'n')){
	    			levls(1,4)++;
	    			if(nht.Storage()>0){
	      			nht = nht | New;
	    			}else{
	      			nht = New;
	    			}
	  			}else{
	    			levls(1,1)++;
	    			levls = levls | New;
	  			}
				}
      }
      p=strtok(NULL,discard);
    }

    if(fpr.Storage()>0){levls = levls | fpr;}
    if(fdr.Storage()>0){levls = levls | fdr;}
    if(nht.Storage()>0){levls = levls | nht;}
    
  //  cerr << " levles : " << levls << endl << endl;
    Res = threshold(data, levls);
    set_threshmaps(Res);

    return Res;
  }

  Matrix MelGMix::threshold(const RowVector& dat, Matrix& levels){  
  	Matrix tests;
    tests=levels;
    Matrix Nprobs;

    //if only single Gaussian: resort to nht
    if(nummix==1||props(1)>0.999||probmap.Sum()<0.05){
    	if(levels(1,4)==0){
				Matrix New(1,6);
				New(1,5)=0.05;
				New(1,6)=0.01;
				New(1,4)=2;New(1,1)=0;New(1,2)=0;New(1,3)=0;tests=New;
      }else{
				Matrix New;
				New = levels.Columns(int(1+levels(1,1)+levels(1,2)
				 	+levels(1,3)),levels.Ncols());
				New(1,4) = levels(1,4);
				New(1,1)=0;New(1,1)=0;New(1,3)=0;
				tests=New;
      }
    }

    int numtests = int(tests(1,1)+tests(1,2)+tests(1,3)+tests(1,4));    
    Matrix Res(numtests,numdata);
    Res = 0.0;
    int next = 1;

    for(int ctr1=1;ctr1<=tests(1,1);ctr1++){
      if(4+next <= tests.Ncols()){
				message("   alternative hypothesis test at p > " << tests(1,4+next) << endl);
				add_infstr(" alternative hypothesis test at p > "+float2str(tests(1,4+next),0,2,false));
				Matrix tmpNull;
				tmpNull = dat;
/*	
				float cutoffpos, cutoffneg;
				cutoffpos = means(1)+100*std::sqrt(vars(1)+0.0000001);
				cutoffneg = means(1)-100*std::sqrt(vars(1)+0.0000001);
	
				for(int ctr=1; ctr<=tmpNull.Ncols(); ctr++)
	  			if( probmap(ctr) > tests(1,4+next) ){
	    			if( dat(ctr) > means(1) )
	      			cutoffpos = std::min(cutoffpos, float(dat(ctr)));
	    			else
	      			cutoffneg = std::max(cutoffneg, float(dat(ctr)));
	 				}
	
				for(int ctr=1; ctr<=tmpNull.Ncols(); ctr++)
	  			if( (dat(ctr) > cutoffneg) && (dat(ctr) < cutoffpos) )
	    			tmpNull(1,ctr)=0.0;
*/
				for(int ctr=1; ctr<=tmpNull.Ncols(); ctr++)
					if( probmap(ctr) < tests(1,4+next) ){	
						tmpNull(1,ctr)=0.0;
						}
							
				Res.Row(next) << tmpNull;
      }
      next++;
    }
    
    for(int ctr1=1;ctr1<=tests(1,2);ctr1++){
      if(4+next <=tests.Ncols()){
				cerr << " false positives control " << tests(1,4+next)<<endl;
				Matrix tmp1;
				tmp1 = normcdf(dat,means(1),vars(1));
				Matrix tmpNull;
				tmpNull = dat; 
				for(int ctr=1; ctr<=tmp1.Ncols(); ctr++)
	  			if(tmp1(1,ctr) < tests(1,4+next))
	    		tmpNull(1,ctr)=0.0;
				Res.Row(next) << tmpNull;
      }
      next++;
    }

    for(int ctr1=1;ctr1<=tests(1,3);ctr1++){
      if(4+next <=tests.Ncols()){
				message("   Local False Discovery Rate control at p < " << tests(1,4+next) << endl);
				add_infstr(" Local False Discovery Rate control at p < "+float2str(tests(1,4+next),0,2,false));
				RowVector tmp=dat;
				SortAscending(tmp);
				RowVector newcdf(tmp);
	  		newcdf << normcdf(tmp,means(1),vars(1));

				float thrp = tmp(tmp.Ncols())+0.01;
				float thrn = tmp(1)-0.01;
				int ctr=tmp.Ncols();
				do{
					thrp = tmp(ctr);
					ctr-=1;
				}while(ctr>0 && ( (1.0-newcdf(ctr))*tmp.Ncols() < (tests(1,4+next)*(tmp.Ncols()-ctr+1))   ));

				ctr=1;
				do{
					thrn = tmp(ctr);
					ctr+=1;
				}while(ctr<=tmp.Ncols() && ( (newcdf(ctr))*tmp.Ncols() < (tests(1,4+next)*ctr)));

				tmp = dat;
				for(ctr=1; ctr<=tmp.Ncols();ctr++)
					if((tmp(ctr) < thrp)&&(tmp(ctr) > thrn))
					tmp(ctr) = 0.0;
				Res.Row(next) << tmp;
      }
			next++;
    }

    for(int ctr1=1;ctr1<=tests(1,4);ctr1++){
      if(4+next <=tests.Ncols()){
				message("   2-sided null hypothesis test at " << tests(1,4+next)<<endl);
				add_infstr(" 2-sided null hypothesis test at "+float2str(tests(1,4+next),0,2,false));
				double mu, sig;
				mu  = dat.Sum()/numdata;
				sig = var(dat,2).AsScalar();
				Matrix tmp1;
				tmp1 = normcdf(dat,mu,std::abs(sig));
				Matrix tmpNull;
				tmpNull = dat; 
				for(int ctr=1; ctr<=tmp1.Ncols(); ctr++)
	  			if((tmp1(1,ctr) < 1-0.5*(tests(1,4+next))&&
	      		(tmp1(1,ctr) > 0.5*(tests(1,4+next)))))
	    				tmpNull(1,ctr)=0.0;
				Res.Row(next) << tmpNull;
      }
      next++;
    }
   
    return Res;
	}

  /* GMM fitting  */

  void MelGMix::gmmupdate(){
    int it_ctr = 1;
    bool exitloop = false;
    float oldll;

    Matrix tmp0;Matrix tmp1;Matrix prob_K__y_theta;
    Matrix kdata;
    RowVector prob_Y__theta;RowVector Nbar;
    RowVector mubar;RowVector sigmabar;RowVector pibar;
    
    do{
      oldll = logprobY;

      //make sure all variances are acceptable
      for(int ctr1=1; ctr1 <=vars.Ncols(); ctr1++)
      	if(vars(ctr1)<0.0001){
      	  vars(ctr1) = 0.0001;
      	}

      tmp0 = normpdf(data,means,vars);
      tmp1 = SP(props.t()*ones(1,numdata),tmp0);      
      prob_Y__theta << sum(tmp1,1);
      logprobY = log(prob_Y__theta).Sum();
      prob_K__y_theta = SP(tmp1,pow(ones(nummix,1)*prob_Y__theta,-1));
      Nbar << sum(prob_K__y_theta,2).t();
      pibar = Nbar / numdata;
      kdata = ones(nummix,1)*data;
      mubar <<SP(sum(SP(kdata,prob_K__y_theta),2).t(),pow(Nbar,-1));    
      kdata -= mubar.t()*ones(1,numdata);
      kdata = pow(kdata,2);
      sigmabar << SP(sum(SP(kdata,prob_K__y_theta),2).t(),pow(Nbar,-1));
      
      means = mubar;
      vars  = sigmabar;
      props = pibar;

      if(epsilon<0){exitloop = it_ctr >= -epsilon;}
      else{exitloop = (((logprobY-oldll < epsilon)&&(it_ctr>20))
		  	||(it_ctr>400));}      
      it_ctr++;
    }while(!exitloop);
  }

  void MelGMix::gmmfit(){
    int i,j;

    if(fixdim){
      if(nummix>1){
				gmmupdate();
				add_params(means,vars,props,logprobY,MDL,Evi,true);
      }else{
				means.ReSize(1);
				means = data.Sum()/numdata;
				vars.ReSize(1);
				vars = var(data,2);
				props.ReSize(1);
				props = 1.0;
				gmmevidence();
      }
    }else{
      RowVector Score(Params.Ncols());
      do{
				gmmupdate();
				Score(nummix) = gmmevidence();    
				int idx1,idx2;
				RowVector pitmp = props;
     
				pitmp.MaximumAbsoluteValue1(idx1);
				pitmp(idx1)=0.0;
				pitmp.MaximumAbsoluteValue1(idx2);
	
				if(props.MaximumAbsoluteValue1(i)<0.9){
	  			if((vars(idx2)>0.15)&&
	     			(std::abs(means(idx2)-means(idx1))<0.5*vars(idx1))){
	    				Score(nummix) = Score(nummix)/(2*(means(idx1)));
	  				}	     
				}
	
				add_params(means,vars,props,logprobY,MDL,Evi,true);
     
				gmmreducemm();
				means = means.Columns(1,nummix);
				vars  = vars.Columns(1,nummix);
				props = props.Columns(1,nummix);

      }while(nummix>1);
      means.ReSize(1);
      means = data.Sum()/numdata;
      vars.ReSize(1);
      vars = var(data,2);
      props.ReSize(1);
      props = 1.0;
      Score(nummix) = gmmevidence();
      add_params(means,vars,props,logprobY,MDL,Evi,true);
      //identify best MM
      Score.MinimumAbsoluteValue2(i,j);
      means.ReSize(j);
      vars.ReSize(j);
      props.ReSize(j);
      nummix = j;
      int index; index = 3 + (j-1)*5; 
      means = Params.SubMatrix(index,index,1,j);
      vars  = Params.SubMatrix(index+1,index+1,1,j);
      props = Params.SubMatrix(index+2,index+2,1,j);
    }

    props.MaximumAbsoluteValue2(i,j);
    if(j>1){
      float tmp;
      tmp = means(1);means(1) = means(j);means(j)=tmp;
      tmp = vars(1);vars(1) = vars(j);vars(j)=tmp;
      tmp = props(1);props(1) = props(j);props(j)=tmp;
    }
    
    add_params(means,vars,props,logprobY,MDL,Evi,true);

    if(nummix==1)
      probmap << normcdf(data,means(1),vars(1));
    else{
      Matrix Nprobs;
      Matrix tmp0;
      tmp0 = normpdf(data,means,vars);
      Nprobs = SP(props.t()*ones(1,numdata),tmp0);
      tmp0 = ones(Nprobs.Nrows(),1)*pow(sum(Nprobs,1),-1);
      Nprobs = SP(tmp0,Nprobs);
      probmap << SP(sum(Nprobs.Rows(2,Nprobs.Nrows()),1),
		    pow(sum(Nprobs,1),-1));
    } 
  }

  float MelGMix::gmmevidence(){
    Matrix tmp0;
    if(means.Ncols()>1){
      tmp0 = normpdf(data,means,vars); 
    }else{
      tmp0 = normpdf(data,means.AsScalar(),vars.AsScalar());
    }
    Matrix tmp1;
    tmp1 = SP(props.t()*ones(1,numdata),tmp0);
    tmp0 = SP(tmp0,pow(ones(nummix,1)*sum(tmp1,1),-1));
    tmp0 = pow(tmp0-ones(nummix,1)*tmp0.Row(nummix),2);
    float logH = 0;
    if(means.Ncols()>1){
      logH = sum(log(sum(tmp0.Rows(1,nummix-1),2)),1).AsScalar();
    }
    logH = logH + 2*sum(log(std::sqrt(2.0)*numdata*props),2).AsScalar();
    logH = logH - 2*sum(props,2).AsScalar();
    
    RowVector prob_Y__theta;
    prob_Y__theta << sum(tmp1,1);    
    logprobY = log(prob_Y__theta).Sum();     
    MDL = -logprobY + (1.5*nummix-0.5)* std::log(float(numdata));   
    Evi = -logprobY +nummix*std::log(2.0)-std::log(MISCMATHS::gamma(nummix))
      -3*nummix/2*std::log(M_PI)+0.5*logH;
  
    return Evi;
  }

  void MelGMix::gmmreducemm(){
    Matrix dlm(nummix,nummix);
    Matrix mus(nummix,nummix);
    Matrix rs(nummix,nummix);

    for(int ctri=1;ctri<=nummix; ctri++){
      for(int ctrj=1;ctrj<=nummix; ctrj++){
				mus(ctri,ctrj) = (props(ctri)*means(ctri)+props(ctrj)*means(ctrj))
	      	/(props(ctri)+props(ctrj));
        rs(ctri,ctrj)  = (props(ctri)*(vars(ctri)+
			  std::pow(means(ctri)-mus(ctri,ctrj),2) ) 
        	+ props(ctrj)*(vars(ctrj) 
         	+ std::pow(means(ctrj)-mus(ctri,ctrj),2))) 
	        / (props(ctri)+props(ctrj));
				dlm(ctri,ctrj) = 0.5*numdata*(
			 		props(ctri)*std::log(
          std::abs(rs(ctri,ctrj))/std::abs(vars(ctri))) 
			 		+ props(ctrj)*std::log(std::abs(rs(ctri,ctrj))
          / std::abs(vars(ctrj))));
      }
    }

    dlm += IdentityMatrix(nummix)*dlm.Maximum();

    int i,j;
    float val;
    val=dlm.MinimumAbsoluteValue2(i,j);

    nummix--;
    
    RowVector newmean(nummix);
    RowVector newvars(nummix);
    RowVector newprop(nummix);

    int ctr0 = 1;
    for(int ctr=1; ctr<=nummix+1; ctr++){
      if(ctr!=i&&ctr!=j){
				newmean(ctr0) = means(ctr);
				newvars(ctr0) = vars(ctr);
				newprop(ctr0) = props(ctr);
				ctr0++;
      }
    }
    //cerr << "ctr0   " << ctr0 << endl;
    if(ctr0<=nummix){
      newmean(ctr0) = mus(i,j);
      newvars(ctr0) = rs(i,j);
      newprop(ctr0) = props(i)+props(j);
      
      means = newmean;    
      vars=newvars;
      props=newprop;
    }
  }

  void MelGMix::ggmfit(){
  	// fit a mixture of a Gaussian and multiple Gamma functions to the histogram
  
    float log_p_y_theta = 1.0;
    float old_ll = 2.0;
    float g_eps = 0.000001;
    int it_ctr = 0;
    double Dmax, Dmin;
   
    Dmax = 2 * data.Maximum();
    Dmin = -2 * data.Minimum();

    //resize means, vars and props
    if(nummix > 3)
      nummix = 3;
    means = means.Columns(1,nummix);
    vars  = vars.Columns(1,nummix);
    props = props.Columns(1,nummix);

    means(1) = -2*mean(data,2).AsScalar();

    Matrix tmp1;Matrix prob_K__y_theta;
    Matrix kdata;
    RowVector prob_Y__theta;RowVector Nbar;
    RowVector mubar;RowVector sigmabar;RowVector pibar;
    Matrix p_ygx(nummix,numdata);
    offset = 0.0;
    float const2;
    Matrix negdata(data);
    negdata = -data;

    while((it_ctr<30) ||
	  ((std::abs(old_ll - log_p_y_theta)>g_eps) && (it_ctr<500))){ // fit GGM	
 			it_ctr++;
			//offset = (std::min(0.2,1-props(1)))*std::sqrt(vars(1));

			//make sure all variances are acceptable
 			for(int ctr1=1; ctr1 <=nummix; ctr1++)
 	  		if(vars(ctr1)<0.0001){
 	    		vars(ctr1) = 0.0001;
 	  		}

 				p_ygx = 0.0;
 				p_ygx.Row(1) << normpdf(data,means(1),vars(1));
       
 				const2 = (2.6-props(1))*sqrt(vars(1))+means(1); //min. nht level
 
				means(2) = (std::max(means(2), std::max(0.001,
	   			0.5 * ( const2 + std::sqrt( const2*const2 + 4*vars(2))))));
				vars(2)  = std::max(std::min(vars(2), 0.5*std::pow(means(2),2)),0.0001);
				p_ygx.Row(2) << gammapdf(data,means(2),vars(2));
   
 				if(nummix>2){
	  			const2 = (2.6-props(1))*sqrt(vars(1))-means(1);
	
	  			means(3) = -(std::max(-means(3), std::max(0.001,
	      		0.5 * ( const2 + std::sqrt( const2*const2 + 4*vars(3))))));
	  			vars(3)  = std::max(std::min(vars(3), 0.5*std::pow(means(3),2)),0.0001);
 	  			p_ygx.Row(3) << gammapdf(negdata,-means(3),vars(3));
				}

 				tmp1 = SP(props.t()*ones(1,numdata),p_ygx);
 				prob_Y__theta << sum(tmp1,1);
	
				//deal with non-modelled voxels
				for(int ctr=1; ctr<=tmp1.Ncols(); ctr++)
	  			if(prob_Y__theta(ctr) < 0.0001)
	    			prob_Y__theta(ctr) = 0.0001;

 				old_ll = log_p_y_theta;
 				log_p_y_theta = log(prob_Y__theta).Sum();
 				if((it_ctr<30) ||
	   			((std::abs(old_ll - log_p_y_theta)>g_eps) && (it_ctr<300))){//update
	  
 	  			prob_K__y_theta = SP(tmp1,pow(ones(nummix,1)*prob_Y__theta,-1));
 	  			Nbar << sum(prob_K__y_theta,2).t();
	  			for(int ctr=1; ctr<=Nbar.Ncols(); ctr++)
	    			if(Nbar(ctr) < 0.0001 * numdata)
	      			Nbar = Nbar + 0.0001;
 	  			pibar= Nbar / numdata;
	  			// 	  cerr << "pibar :" << pibar << endl;
 	  			kdata = ones(nummix,1)*data;
 	  			mubar <<SP(sum(SP(kdata,prob_K__y_theta),2).t(),pow(Nbar,-1)); 
	  			// 	  cerr << "mubar :" << mubar << endl;

 	  			kdata -= mubar.t()*ones(1,numdata);
 	  			kdata = pow(kdata,2);
 	  			sigmabar << SP(sum(SP(kdata,prob_K__y_theta),2).t(),pow(Nbar,-1));
      
 	  			means = mubar;
 	  			vars  = sigmabar;
 	  			props = pibar;
 					}//update
    } //while loop

    props = props / sum(props,2).AsScalar();
    add_params(means,vars,props,logprobY,MDL,Evi,true);
	
    probmap << SP(sum(tmp1.Rows(2,tmp1.Nrows()),1),
		  pow(sum(tmp1,1),-1));

	dbgmsg("   mu: " << means << "  sig: " << vars << " prop: " << props << endl);

   	if(props(1)<0.4 ){
      //set up GMM
      message("    try Gaussian Mixture Model " << endl);
      props=zeros(1,nummix);
      vars=zeros(1,nummix);
      means=zeros(1,nummix);
      Params=zeros(1,nummix);
      logprobY = 1.0;  
      props = std::pow(float(nummix),float(-1.0));

      tmp1 = data * data.t() / numdata;
      vars = tmp1.AsScalar();
      float Dmin, Dmax, IntSize;
      Dmin = min(data).AsScalar(); Dmax = max(data).AsScalar();
      IntSize = Dmax / nummix;
      means(1)=mean(data,2).AsScalar(); 
      for (int ctr=2; ctr <= means.Ncols(); ctr++){
		means(ctr) =  Dmax - (ctr - 1.5) * IntSize; 
      }
      means(2)=means(1)+sqrt(vars(1));
      if(nummix>2)
		means(3)=means(1)-sqrt(vars(1));
      
      fit(string("GMM"));
    }

  }

  /*  INPUT / OUTPUT  */

  void MelGMix::add_params(Matrix& mu, Matrix& sig, Matrix& pi, 
		float logLH, float MDL, float Evi, bool advance){ 
    int size = Params.Ncols();
    if(size<2){size=2;}
    Matrix New(5,size);
    New = 0;
    
    New.SubMatrix(3,3,1,mu.Ncols())=mu;
    New.SubMatrix(4,4,1,mu.Ncols())=sig;
    New.SubMatrix(5,5,1,mu.Ncols())=pi;
    New(1,1)=nummix;
  
    New(1,2)=-logLH;
    New(2,1)=Evi;
    New(2,2)=MDL;
    if(Params.Storage()>nummix){ 
      Params = New & Params;
    }else{ 
      Params =  New;
    }
    }

  void MelGMix::get_params(int index, Matrix& mu, Matrix& sig, Matrix& pi, 
		float logLH, float MDL, float Evi){ 
   
    }

  void MelGMix::save(){
    
  }

  void MelGMix::status(const string &txt){
    cerr << txt << "epsilon " << epsilon << endl;
    cerr << txt << "nummix  " << nummix << endl;
    cerr << txt << "numdata " << numdata << endl;
    cerr << txt << "means   " << means << endl;
    cerr << txt << "vars    " << vars << endl;
    cerr << txt << "props   " << props << endl;
    //write_ascii_matrix(logger.appendDir(string(txt + "means")),means);
    //write_ascii_matrix(logger.appendDir(string(txt + "vars")),vars);
    //write_ascii_matrix(logger.appendDir(string(txt + "props")),props);
  }

  void MelGMix::create_rep(){
 
  }
  

}



