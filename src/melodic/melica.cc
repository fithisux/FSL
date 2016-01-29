/*  MELODIC - Multivariate exploratory linear optimized decomposition into 
              independent components
    
    melica.cc - ICA estimation 

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

#include <stdlib.h>
#include "utils/log.h"
#include "meloptions.h"
#include "meldata.h"
#include "melodic.h"
#include "newmatap.h"
#include "newmatio.h"
#include "melica.h"
#include "melpca.h"
#include "melhlprfns.h"
#include "miscmaths/miscprob.h"

using namespace Utilities;
using namespace NEWIMAGE;

namespace Melodic {
    
  void MelodicICA::ica_fastica_symm(const Matrix &Data){
    // based on Aapo Hyvärinen's fastica method
    // see www.cis.hut.fi/projects/ica/fastica/
    
    //initialize matrices
    Matrix redUMM_old, rank1_old;
    Matrix tmpU;    
    //srand((unsigned int)timer(NULL));
    redUMM = melodat.get_white()*
       unifrnd(melodat.get_white().Ncols(),dim); // got to start somewhere
    
	if(opts.debug.value())
		cerr << "redUMM init submatrix : " << endl << redUMM.SubMatrix(1,2,1,2) << endl;
		
    if(opts.guessfname.value().size()>1){
      message("  Use columns in " << opts.guessfname.value() 
	      << " as initial values for the mixing matrix " <<endl);
      Matrix guess ;
      guess  = melodat.get_white()*read_ascii_matrix(opts.guessfname.value());
      redUMM.Columns(1,guess.Ncols()) = guess;
    }
    
    symm_orth(redUMM);

    int itt_ctr,itt_ctr2=1,cum_itt=0,newmaxitts = opts.maxNumItt.value(); 
    double minAbsSin = 1.0, minAbsSin2 = 1.0;
    if(opts.approach.value() == string("tica")) {
      opts.maxNumItt.set_T(opts.rank1interval.value());
      rank1_old = melodat.get_dewhite()*redUMM;
      rank1_old = melodat.expand_dimred(rank1_old);
      rank1_old = krapprox(rank1_old,int(rank1_old.Nrows()/melodat.get_numfiles()));
    }
    do{// TICA loop
      itt_ctr = 1;
      do{ // da loop!!!
				redUMM_old = redUMM;      
				//calculate IC estimates
				tmpU = Data.t() * redUMM;
					
				//update redUMM depending on nonlinearity
				if(opts.nonlinearity.value()=="pow4"){
	  			redUMM = (Data * pow(tmpU,3.0)) / samples - 3 * redUMM;
				}
				if(opts.nonlinearity.value()=="pow3"){
	  			tmpU /= opts.nlconst1.value();
	  			redUMM = 3 * (Data * pow(tmpU,2.0)) / samples  - 
	    			(SP(ones(dim,1)*sum(tmpU,1),redUMM))/ samples;
				}
				if(opts.nonlinearity.value()=="tanh"){
	  			Matrix hyptanh;
	  			hyptanh = tanh(opts.nlconst1.value()*tmpU);
	  			redUMM = (Data * hyptanh - opts.nlconst1.value()*SP(ones(dim,1)*
						sum(1-pow(hyptanh,2),1),redUMM))/samples;						
				}
				if(opts.nonlinearity.value()=="gauss"){
	  			Matrix tmpUsq;
	  			Matrix tmpU2;
	  			Matrix gauss;
	  			Matrix dgauss;
	  			tmpUsq = pow(tmpU,2);
	  			tmpU2 = exp(-(opts.nlconst2.value()/2) * tmpUsq);
	  			gauss = SP(tmpU,tmpU2);
	  			dgauss = SP(1-opts.nlconst2.value()*tmpUsq,tmpU2);
	  			redUMM = (Data * gauss - SP(ones(dim,1)*
				    sum(dgauss,1),redUMM))/samples;
				}
           
				// orthogonalize the unmixing-matrix 
				symm_orth(redUMM);
				if(opts.approach.value() == string("tica")){
					message("");
				}

				//termination condition : angle between old and new < epsilon
				minAbsSin = 1 - diag(abs(redUMM.t()*redUMM_old)).Minimum();

				message("  Step no. " << cum_itt + itt_ctr << " change : " << minAbsSin << endl);
		//		if((abs(minAbsSin) < opts.epsilon.value())&&
		//		  (opts.approach.value()!=string("tica"))){ break;}
				if((abs(minAbsSin) < opts.epsilon.value())){ break;}
	
				itt_ctr++;
      } while(itt_ctr < opts.maxNumItt.value());
      cum_itt += itt_ctr;
      itt_ctr2++;
      if(opts.approach.value() == string("tica")){
				message(" Rank-1 approximation of the time courses; ");
	  		Matrix temp(melodat.get_dewhite() * redUMM);
	  		temp = melodat.expand_dimred(temp);
			  temp = krapprox(temp,int(temp.Nrows()/melodat.get_numfiles())); 
				minAbsSin2 = 1 - diag(abs(corrcoef(temp,rank1_old))).Minimum();
				rank1_old = temp;
				temp = melodat.reduce_dimred(temp);
				redUMM = melodat.get_white() * temp;

				message(" change : " << minAbsSin2 << endl);	
				if(abs(minAbsSin2) < opts.epsilonS.value() && abs(minAbsSin) < opts.epsilon.value()){ break;}
			}
    } while(
      (itt_ctr2 < newmaxitts/opts.maxNumItt.value()) && 
			(opts.approach.value() ==  string("tica")) && 
			cum_itt < newmaxitts);

    if((itt_ctr>=opts.maxNumItt.value() && (opts.approach.value()!=string("tica")))
			|| (cum_itt >= newmaxitts && opts.approach.value()==string("tica"))){
      cerr << "  No convergence after " << cum_itt  <<" steps "<<endl;
    } else {
      message("  Convergence after " << cum_itt  <<" steps " << endl << endl);
      no_convergence = false;

      {Matrix temp(melodat.get_dewhite() * redUMM);
       melodat.set_mix(temp);}
      {Matrix temp(redUMM.t()*melodat.get_white());
      melodat.set_unmix(temp);}
    } 
  }
  
  void MelodicICA::ica_fastica_defl(const Matrix &Data){
    if(!opts.explicitnums || opts.numICs.value()>dim){
      opts.numICs.set_T(dim); 
      message("  Using numICs:" << opts.numICs.value() << endl);
    }
     
    //redUMM = zeros(dim); // got to start somewhere
    redUMM = melodat.get_white()*
      unifrnd(melodat.get_white().Ncols(),opts.numICs.value());
    redUMM = zeros(melodat.get_white().Nrows(),opts.numICs.value());
    Matrix guess;
    int guesses=0;
    if(opts.guessfname.value().size()>1){
       message("  Use columns in " << opts.guessfname.value() << " as initial values for the mixing matrix " <<endl);
       guess  = melodat.get_white()*read_ascii_matrix(opts.guessfname.value());
       guesses = guess.Ncols();
    }

    int ctrIC = 1;
    int numRestart = 0;
    while(ctrIC<=opts.numICs.value()){   
     	message("  Extracting IC " << ctrIC << "  ... ");
      ColumnVector w;
      ColumnVector w_old;   
      ColumnVector tmpU;
      if(ctrIC <= guesses){
      	w = w - redUMM * redUMM.t() * w;
      	w = w / norm2(w);  
      	w_old = zeros(w.Nrows(),1);
      	int itt_ctr = 1; 
      	do{
	 				w_old = w;
	 				tmpU = Data.t() * w; 
					if(opts.nonlinearity.value()=="pow4"){
	  				w =  (Data * pow(tmpU,3.0)) / samples - 3 * w;
					}
					if(opts.nonlinearity.value()=="tanh"){
	  				ColumnVector hyptanh;
          	hyptanh = tanh(opts.nlconst1.value()*tmpU);
 	  				w = (Data * hyptanh - opts.nlconst1.value()*SP(ones(dim,1)*
          		sum(1-pow(hyptanh,2),1),w))/samples;
					} 
					if(opts.nonlinearity.value()=="pow3"){
	  				tmpU /= opts.nlconst1.value();
 	  				w = 3*(Data * pow(tmpU,2.0)) / samples - 2*(SP(ones(dim,1)*
           		sum(tmpU,1),w))/samples;
					} 
					if(opts.nonlinearity.value()=="gauss"){
	  				ColumnVector tmpUsq;
	  				ColumnVector tmpU2;
	  				ColumnVector gauss;
	  				ColumnVector dgauss;
	  				tmpUsq = pow(tmpU,2);
	  				tmpU2 = exp(-(opts.nlconst2.value()/2) * tmpUsq);
	  				gauss = SP(tmpU,tmpU2);
	  				dgauss = SP(1-opts.nlconst2.value()*tmpUsq,tmpU2);
          	w = (Data * gauss - SP(ones(dim,1)*
				 			sum(dgauss,1),w))/samples;
					}

					// orthogonalize w
					w = w - redUMM * redUMM.t() * w;
					w = w / norm2(w);  

					//termination condition : angle between old and new < epsilon
					if(((norm2(w-w_old) < 0.001*opts.epsilon.value())&&(itt_ctr>10)) || 
					   ((norm2(w+w_old) < 0.001*opts.epsilon.value())&&(itt_ctr>10))){
	 					break;
				  	}
        	//cout << norm2(w-w_old) << "   " << norm2(w+w_old) << endl;
					itt_ctr++;
     	  } while(itt_ctr <= opts.maxNumItt.value());

      if(itt_ctr<opts.maxNumItt.value()){
				redUMM.Column(ctrIC) = w;
        message(" estimated using " << itt_ctr << " iterations " << endl);
        ctrIC++; 
        numRestart = 0;
      } else{
        if(numRestart > opts.maxRestart.value()){
	  			message(endl << "  Estimation failed after " 
		  			<< numRestart << " attempts " 
		  			<< " giving up " << endl);
	  			break;
				}else{
          numRestart++;
	  			message(endl <<"  Estimation failed - restart " 
		  			<< numRestart << endl);
				}
      }
     	}
     	if(numRestart < opts.maxRestart.value()){
       	no_convergence = false;
       	{Matrix temp(melodat.get_dewhite() * redUMM);
       		melodat.set_mix(temp);}
       	{Matrix temp(redUMM.t()*melodat.get_white());
       		melodat.set_unmix(temp);}
     	}
		}
  }

  void MelodicICA::ica_maxent(const Matrix &Data){
    // based on Aapo Hyvärinen's fastica method
    // see www.cis.hut.fi/projects/ica/fastica/ 
    message(" MAXENT " << endl);
    //initialize matrices
    Matrix redUMM_old;
    Matrix tmpU;    
    Matrix gtmpU;
    double lambda = 0.015/std::log((double)melodat.get_white().Ncols());
    
    //srand((unsigned int)timer(NULL));
    redUMM = melodat.get_white()*
     	unifrnd(melodat.get_white().Ncols(),dim); // got to start somewhere
    
    if(opts.guessfname.value().size()>1){
      message("  Use columns in " << opts.guessfname.value() 
	      << " as initial values for the mixing matrix " <<endl);
      Matrix guess ;
      guess  = melodat.get_white()*read_ascii_matrix(opts.guessfname.value());
      redUMM.Columns(1,guess.Ncols()) = guess;
    }
    
    //    symm_orth(redUMM);
    int itt_ctr=1; 
    double minAbsSin = 1.0;
    Matrix Id;
    Id = IdentityMatrix(redUMM.Ncols());
    //cerr << " nonlinearity : " <<    opts.nonlinearity.value() << endl;

    do{ // da loop!!!
      redUMM_old = redUMM;      
      //calculate IC estimates
      tmpU = Data.t() * redUMM;
      if(opts.nonlinearity.value()=="tanh"){
				//Matrix hyptanh;
				//hyptanh = tanh(opts.nlconst1.value()*tmpU);
				//redUMM = (Data * hyptanh - opts.nlconst1.value()*SP(ones(dim,1)*
				//sum(1-pow(hyptanh,2),1),redUMM))/samples;
				gtmpU = tanh(opts.nlconst1.value()*tmpU);
				redUMM = redUMM + lambda*(Id+(1-2*gtmpU.t()*tmpU))*redUMM;
      }
      if(opts.nonlinearity.value()=="gauss"){
				gtmpU = pow(1+exp(-(opts.nlconst2.value()/2) * tmpU),-1);
				redUMM = redUMM + lambda*(Id - (gtmpU.t()-tmpU.t())*tmpU)*redUMM;
      }
 
      //termination condition : angle between old and new < epsilon
      minAbsSin = abs(1 - diag(abs(redUMM.t()*redUMM_old)).Minimum());
      message("  Step no. " << itt_ctr << " change : " << minAbsSin << endl);
      if(abs(minAbsSin) < opts.epsilon.value()){ break;}
      
      itt_ctr++;
    } while(itt_ctr < opts.maxNumItt.value());
    
    if(itt_ctr>=opts.maxNumItt.value()){
      cerr << "  No convergence after " << itt_ctr <<" steps "<<endl;
    } else {
      message("  Convergence after " << itt_ctr <<" steps " << endl << endl);
      no_convergence = false;
      {Matrix temp(melodat.get_dewhite() * redUMM);
       melodat.set_mix(temp);}
      {Matrix temp(redUMM.t()*melodat.get_white());
      melodat.set_unmix(temp);}
    } 
  }
  
  void MelodicICA::ica_jade(const Matrix &Data){ 
    int dim_sym = (int) (dim*(dim+1))/2;  
    int num_CM = dim_sym;
    Matrix CM;
    Matrix R; R = IdentityMatrix(dim);
    Matrix Qij; Qij = zeros(dim);
    Matrix Xim;
    Matrix Xjm;
    Matrix scale; scale = ones(dim,1)/samples;

    for (int im =1; im <= dim; im++){
      Xim = Data.Row(im);
			write_ascii_matrix("Xim",Data.Row(1));
      //Qij = SP((scale * pow(Xim,2)),Data) * Data.t();//- R - 2*R.Column(im)*R.Column(im).t();
      Qij = (pow(Xim,2)) * Data.t();//- R - 2*R.Column(im)*R.Column(im).t();
      if(im==1){CM = Qij; write_ascii_matrix("CM",CM);exit(2);}else{CM |= Qij;}
      for (int jm = 1; jm < im; jm++){
				Xjm = Data.Row(jm);
				Qij = (SP((scale * SP(Xim,Xjm)),Data) * Data.t() - 
					R.Column(im)*R.Column(jm).t() - R.Column(jm)*R.Column(im).t());
				Qij *= sqrt(2);
				CM  |= Qij;
      }
    }

    write_ascii_matrix("CM",CM);
    Matrix redUMM; redUMM = IdentityMatrix(dim);
  
    bool exitloop = false;
    int ctr_itt = 0;
    int ctr_updates = 0;
    Matrix Givens; Givens = zeros(2,num_CM);
    Matrix Givens_ip; Givens_ip = zeros(2);
    Matrix Givens_ro; Givens_ro = zeros(2);
    double det_on, det_off;
    double cos_theta, sin_theta, theta;

    while(!exitloop && ctr_itt <= opts.maxNumItt.value()){
      ctr_itt++;
      cout << "loop" <<endl;
      for(int ctr_p = 1; ctr_p < dim; ctr_p++){
				for(int ctr_q = ctr_p+1; ctr_q <= dim; ctr_q++){

	  			for(int ctr_i = 0; ctr_i < num_CM; ctr_i++){
	    			int Ip = ctr_p + ctr_i * dim;
	    			int Iq = ctr_q + ctr_i * dim;
	    			Givens(1,ctr_i + 1) = CM(ctr_p,Ip) - CM(ctr_q,Iq);
	    			Givens(2,ctr_i + 1) = CM(ctr_p,Iq) - CM(ctr_q,Ip);
	  			}
	  
	  			Givens_ip = Givens * Givens.t();
	  			det_on = Givens_ip(1,1) - Givens_ip(2,2);
	  			det_off = Givens_ip(1,2) + Givens_ip(2,1);
	  			theta = 0.5 * atan2(det_off, det_on + sqrt(det_on*det_on + det_off*det_off));

	  			cout << theta << endl;

	  			if(abs(theta) > opts.epsilon.value()){
	    			ctr_updates++;
	    			message("  Step No. "<< ctr_itt << " change : " << theta << endl);
			
	    			//create Givens rotation matrix
	    			cos_theta = cos(theta);
	    			sin_theta = sin(theta);
	    			Givens_ro(1,1) = cos_theta;
	    			Givens_ro(1,2) = -sin_theta;
	    			Givens_ro(2,1) = sin_theta;
	    			Givens_ro(2,2) = cos_theta;

	    			//update 2 entries of redUMM
	    			{Matrix tmp_redUMM;
	    				tmp_redUMM = redUMM.Column(ctr_p);
	    				tmp_redUMM |= redUMM.Column(ctr_q);
	    				tmp_redUMM *= Givens_ro;
	    				redUMM.Column(ctr_p) = tmp_redUMM.Column(1);
	    				redUMM.Column(ctr_q) = tmp_redUMM.Column(2);}

	    			//update Cumulant matrix
	    			{Matrix tmp_CM;
	    				tmp_CM = CM.Row(ctr_p);
	    				tmp_CM &= CM.Row(ctr_q);
	    				tmp_CM = Givens_ro.t() * tmp_CM;
	    				CM.Row(ctr_p) = tmp_CM.Row(1);
	    				CM.Row(ctr_q) = tmp_CM.Row(2);}

	    			//update Cumulant matrices
	    			for(int ctr_i = 0; ctr_i < num_CM; ctr_i++){
	      			int Ip = ctr_p + ctr_i * dim;
	      			int Iq = ctr_q + ctr_i * dim;
	      			CM.Column(Ip) = cos_theta*CM.Column(Ip)+sin_theta*CM.Column(Iq);
	      			CM.Column(Iq) = cos_theta*CM.Column(Iq)-sin_theta*CM.Column(Ip);
	    			}
	  			}else{
	    			exitloop = true;
	  			}
				}
      }
    }//while loop
    if(ctr_itt > opts.maxNumItt.value()){
    	cerr << "  No convergence after " << ctr_itt <<" steps "<<endl;
    } else {
     	message("  Convergence after " << ctr_itt <<" steps " << endl << endl);
      no_convergence = false;
      {Matrix temp(melodat.get_dewhite() * redUMM);
       	melodat.set_mix(temp);}
      {Matrix temp(redUMM.t()*melodat.get_white());
       	melodat.set_unmix(temp);}
    }
  }

  Matrix MelodicICA::sign(const Matrix &Inp){
    Matrix Res = Inp;
    Res = 1;
    for(int ctr_i = 1; ctr_i <= Inp.Ncols(); ctr_i++){
      for(int ctr_j = 1; ctr_j <= Inp.Nrows(); ctr_j++){
				if(Inp(ctr_j,ctr_i)<0){Res(ctr_j,ctr_i)=-1;}
      }
    } 
    return Res;
  }

  void MelodicICA::perf_ica(const Matrix &Data){ 
    message("Starting ICA estimation using " << opts.approach.value() 
	    << endl << endl);
    dim = Data.Nrows();
    samples = Data.Ncols();
    no_convergence = true;
    //switch to the chosen method   
    if(opts.approach.value()==string("symm") ||
       opts.approach.value()==string("tica") ||
       opts.approach.value()==string("parafac") ||
       opts.approach.value()==string("concat"))
      ica_fastica_symm(Data);
    if(opts.approach.value()==string("defl"))
      ica_fastica_defl(Data);
    if(opts.approach.value()==string("jade"))
      ica_jade(Data);
    if(opts.approach.value()==string("maxent"))
      ica_maxent(Data);
    
    if(!no_convergence){//calculate the IC
      Matrix temp(melodat.get_unmix()*melodat.get_Data());
      //  Add the mean time course again  
      //      temp += melodat.get_unmix()*melodat.get_meanC()*ones(1,temp.Ncols());

      //re-normalise the decomposition to std(mix)=1
      Matrix scales;
      scales = stdev(melodat.get_mix());   

      //cerr << " SCALES 1 " << scales << endl;
      Matrix tmp, tmp2;
      tmp = SP(melodat.get_mix(),ones(melodat.get_mix().Nrows(),1)*pow(scales,-1));
      temp = SP(temp,scales.t()*ones(1,temp.Ncols()));

      scales = scales.t();
  
      melodat.set_mix(tmp);

      melodat.set_IC(temp);

      melodat.set_ICstats(scales);
      melodat.sort();

	  //message("Calculating T- and S-modes " << endl);
      melodat.set_TSmode();
		
    }
  }
}


