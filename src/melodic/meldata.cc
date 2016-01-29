/*  MELODIC - Multivariate exploratory linear optimized decomposition into 
              independent components
    
    meldata.cc - data handler / container class

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
#include "meloptions.h"
#include "meldata.h"
#include "melodic.h" 
#include "utils/log.h"
#include <time.h>
#include <algorithm>
#include "miscmaths/miscprob.h"
#include "melhlprfns.h" 

using namespace Utilities;
using namespace NEWIMAGE;
 
namespace Melodic{
// {{{ Setup

  ReturnMatrix MelodicData::process_file(string fname, int numfiles)
  {
	dbgmsg(string("START: process_file") << endl);	
	
	Matrix tmpData;
	{
    	volume4D<float> RawData;

		memmsg(" before reading file "<< fname);

    	//read data
    	message("Reading data file " << fname << "  ... ");
    	read_volume4D(RawData,fname);
    	message(" done" << endl);
		memmsg(" after reading file "<< fname);
		
		del_vols(RawData,opts.dummy.value());
    
    	Mean += meanvol(RawData)/numfiles;

		//estimate smoothness
		memmsg(" before est smoothness ");
    	if((Resels == 0)&&(!opts.filtermode))
      		Resels = est_resels(RawData,Mask);
		memmsg(" after smoothness ");
	
    	//reshape
		memmsg(" before reshape ");
    	tmpData = RawData.matrix(Mask);
		memmsg(" after reshape ");	  
	}    

  // If a time series model design was specified, check 
  // that the data dimensions match the model dimensions
  if (Tdes.Storage() && (tmpData.Nrows() != Tdes.Nrows())) {

    cerr << "ERROR: " << fname << " " << 
      "- data dimensions (" << tmpData.Nrows() << ") "  << 
      "do not match model dimensions (" << Tdes.Nrows() << ")" << endl;
    exit(2);
  } 
        
    //convert to percent BOLD signal change
    if(opts.pbsc.value()){
      message("  Converting data to percent BOLD signal change ...");
      Matrix meanimg = convert_to_pbsc(tmpData);
      meanR = meanimg.Row(1);
      message(" done" << endl);
    } 
	else{
		if(opts.remove_meanvol.value())
		{	      
			message(string("  Removing mean image ..."));
			memmsg(" before remmean ");
      		remmean(tmpData,meanR,1);
			memmsg(" after remmean ");
      		message(" done" << endl);
		}
		else meanR=ones(1,tmpData.Ncols());
    }

	if(opts.remove_meantc.value()){
		remmean(tmpData,meanC,2);
	}
	
    //convert to power spectra
    if(opts.pspec.value()){
      message("  Converting data to powerspectra ...");
      tmpData = calc_FFT(tmpData);
      message(" done" << endl);
    }
	
    //switch dimension in case temporal ICA is required
    if(opts.temporal.value()){
      message(string("  Switching dimensions for temporal ICA") << endl);
      tmpData = tmpData.t();
      Matrix tmp;
      tmp = meanC;
      meanC = meanR.t();
      meanR = tmp.t();
      message("  Data size : " << Data.Nrows() << " x " << Data.Ncols() <<endl);
    }
      
    //variance - normalisation

    if(opts.varnorm.value()){
	  memmsg(" before VN ");
      message("  Normalising by voxel-wise variance ..."); 
		outMsize("stdDev",stdDev);
//			if(stdDev.Storage()==0)
      	stdDev = varnorm(tmpData,std::min(30,tmpData.Nrows()-1),
					opts.vn_level.value(), opts.econ.value());
//			else 	
//				stdDev += varnorm(tmpData,std::min(30,tmpData.Nrows()-1),
//					opts.vn_level.value(), opts.econ.value())/numfiles;
      stdDevi = pow(stdDev,-1); 
	  memmsg(" in VN ");
      message(" done" << endl);
    }

	//convert to instacorrs
	if(opts.insta_fn.value()>""){
		Matrix vscales = pow(stdev(tmpData,1),-1);
		varnorm(tmpData,vscales);
		
		Matrix tmpTC = tmpData * insta_mask.t();
		varnorm(tmpTC,pow(stdev(tmpTC),-1));
		
		for(int ctr=1; ctr <=tmpData.Ncols();ctr++)
			tmpData.Column(ctr) = SP(tmpData.Column(ctr),tmpTC);

	}

	tmpData.Release();
	dbgmsg(string("END: process_file") << endl);	
    return tmpData;
  }

  Matrix MelodicData::expand_mix()
  {
    Matrix out;
    out = expand_dimred(mixMatrix);
    return out;
  }

  Matrix MelodicData::expand_dimred(const Matrix& Mat)
  {
    int first, last;
    first = 1;
    last = DWM.at(0).Ncols();
    Matrix tmp = DWM.at(0) * Mat.Rows(first,last);
    for(unsigned int ctr = 1; ctr < DWM.size(); ctr++){
      first = last + 1;
      last += DWM.at(ctr).Ncols();
      tmp &= DWM.at(ctr) * Mat.Rows(first, last);
    }
    return tmp;
  }

  Matrix MelodicData::reduce_dimred(const Matrix& Mat)
  {
    int first, last;
    first = 1;
    last = WM.at(0).Ncols();
    Matrix tmp = WM.at(0) * Mat.Rows(first,last);
    for(unsigned int ctr = 1; ctr < WM.size(); ctr++){
      first = last + 1;
      last += WM.at(ctr).Ncols();
      tmp &= WM.at(ctr) * Mat.Rows(first, last);
    }
    return tmp;
  }

  void MelodicData::set_TSmode_depr()
  {
		Matrix tmp, tmpT, tmpS, tmpT2, tmpS2, tmpT3;

	    tmp = expand_dimred(mixMatrix);
	    tmpT = zeros(tmp.Nrows()/numfiles, tmp.Ncols());
	    tmpS = ones(numfiles, tmp.Ncols());

		outMsize("tmp",tmp);
		outMsize("tmpT",tmpT);
		outMsize("tmpS",tmpS);

	  	dbgmsg(string("   approach ") << opts.approach.value() << endl);	

		if(opts.approach.value()==string("tica")){
	      message("Calculating T- and S-modes " << endl);
	      explained_var = krfact(tmp,tmpT,tmpS);
	      Tmodes.clear(); Smodes.clear();
	      for(int ctr = 1; ctr <= tmp.Ncols(); ctr++){
		    tmpT3 << reshape(tmp.Column(ctr),tmpT.Nrows(),numfiles);
			outMsize("tmpT3", tmpT3);
	      	tmpT2 << tmpT.Column(ctr);
	      	tmpS2 << tmpS.Column(ctr);
			tmpT3 << SP(tmpT3,pow(ones(tmpT3.Nrows(),1)*tmpS2.t(),-1));
			if(numfiles>1)
		      tmpT2 |= tmpT3;
			if(mean(tmpS2,1).AsScalar()<0){
			  tmpT2*=-1.0;
			  tmpS2*=-1.0;
			}
	      	add_Tmodes(tmpT2);
	      	add_Smodes(tmpS2);
		  }
		}
		else{
			Tmodes.clear();
			Smodes.clear();
			for(int ctr = 1; ctr <= tmp.Ncols(); ctr++){
				tmpT3 << tmp.Column(ctr);
				add_Tmodes(tmpT3);
			}
		}

	    if(opts.approach.value()!=string("concat")){
		  //add GLM OLS fit
		  dbgmsg(string(" GLM fitting ") << endl);

		  if(Tdes.Storage()){
		    Matrix alltcs = Tmodes.at(0).Column(1);
		    for(int ctr=1; ctr < (int)Tmodes.size();ctr++)
			  alltcs|=Tmodes.at(ctr).Column(1);
		    if((alltcs.Nrows()==Tdes.Nrows())&&(Tdes.Nrows()>Tdes.Ncols()))
			  glmT.olsfit(alltcs,Tdes,Tcon);
		  }
		  if(Sdes.Storage()){
		    Matrix alltcs = Smodes.at(0);
		    for(int ctr=1; ctr < (int)Smodes.size();ctr++)
		  	  alltcs|=Smodes.at(ctr);
		    if((alltcs.Nrows()==Sdes.Nrows())&&(Sdes.Nrows()>Sdes.Ncols()&&alltcs.Nrows()>2))
			  glmS.olsfit(alltcs,Sdes,Scon);
		  }

	    }
	  //    else{
	//		dbgmsg(string(" Bypassing krfac ") << endl);
	//        add_Tmodes(tmp);
	//		add_Smodes(tmpS);
	//      }
  }

  void MelodicData::dual_regression()
  {
	dbgmsg(string("START: dual_regression") << endl);	

	Tmodes.clear();
	Smodes.clear();

	bool tmpvarnorm = opts.varnorm.value();
	// Switch off variance normalisation
	opts.varnorm.set_T(false);

	Log drO;

	if(opts.dr_out.value())
		drO.makeDir(logger.appendDir("dr"),"dr.log");

	Matrix tmpcont = diag(ones(IC.Nrows(),1)), s1,s2, tmpData, alltcs;
	basicGLM tmpglm;
	for(int ctr = 0; ctr < numfiles; ctr++){
		tmpData = process_file(opts.inputfname.value().at(ctr), numfiles);
		//may want to remove the spatial means first
		tmpglm.olsfit(remmean(tmpData.t(),1),remmean(IC.t(),1),tmpcont);
		s1=tmpglm.get_beta().t();

		outMsize("s1",s1);
	    outMsize("alltcs",alltcs);
		if(alltcs.Storage()==0)
			alltcs=s1;
		else
			alltcs&=s1;
			
		// output DR
		if(opts.dr_out.value()){
			
			dbgmsg(string("START: dual_regression output") << endl);
			write_ascii_matrix(drO.appendDir("dr_stage1_subject"+num2str(ctr,4)+".txt"),s1);
			//des_norm
			s1 =  SP(s1,ones(s1.Nrows(),1)*pow(stdev(s1,1),-1));
			tmpglm.olsfit(remmean(tmpData),remmean(s1,1),tmpcont);
			s2=tmpglm.get_beta();
			save4D(s2,string("dr/dr_stage2_subject"+num2str(ctr,4)));
			s2=tmpglm.get_z();
			save4D(s2,string("dr/dr_stage2_subject"+num2str(ctr,4)+"_Z"));
		}
    }

	for(int ctr = 1; ctr <= alltcs.Ncols(); ctr++){
		tmpcont << alltcs.Column(ctr);
		add_Tmodes(tmpcont);
	}
	
	for(int ctrC = 1; ctrC <=IC.Nrows(); ctrC++){
		Matrix tmpall = zeros(numfiles,IC.Ncols());
		string fnout = string("dr/dr_stage2_ic"+num2str(ctrC-1,4));
		for(int ctrS = 0; ctrS < numfiles; ctrS++){
			string fnin = logger.appendDir(string("dr/dr_stage2_subject"+num2str(ctrS,4)));
			dbgmsg(fnout << endl << fnin << endl);
			volume4D<float> vol;
			read_volume4DROI(vol,fnin,0,0,0,ctrC-1,-1,-1,-1,ctrC-1);
			
			Matrix tmp2 = vol.matrix(Mask);
		    tmpall.Row(ctrS+1) << vol.matrix(Mask);
		}
		save4D(tmpall,fnout);
	}

    opts.varnorm.set_T(tmpvarnorm);	
	dbgmsg(string("END: dual_regression") << endl);	
  }

  void MelodicData::set_TSmode()
  {
   	dbgmsg(string("START: set_TSmode")<< endl);	
	if(opts.dr.value())
		dual_regression();
	else
		set_TSmode_depr();
	
	dbgmsg(string("END: set_TSmode")<< endl);	
  }

  void MelodicData::setup_classic()
  {
	    dbgmsg(string("START: setup_classic") << endl);	
    	Matrix alldat, tmpData;
		bool tmpvarnorm = opts.varnorm.value();

		if(numfiles > 1 && opts.joined_vn.value()){
			opts.varnorm.set_T(false);
		}
    	alldat = process_file(opts.inputfname.value().at(0), numfiles) / numfiles;
		memmsg(" after process_file ");

		if(opts.pca_dim.value() > alldat.Nrows()-2){
			cerr << "ERROR:: too many components selected \n\n";
			exit(2);
		}

		for(int ctr = 1; ctr < numfiles; ctr++){
    		tmpData = process_file(opts.inputfname.value().at(ctr), numfiles) / numfiles;
			if(tmpData.Ncols() == alldat.Ncols() && tmpData.Nrows() == alldat.Nrows())
      			alldat = alldat + tmpData;	
			else{
				    if(opts.approach.value()==string("tica")){
						cerr << "ERROR:: data dimensions do not match, TICA not possible \n\n";
						exit(2); 
					}

					if(tmpData.Ncols() == alldat.Ncols()){
						int mindim = min(alldat.Nrows(),tmpData.Nrows());
						alldat = alldat.Rows(1,mindim);
						tmpData = tmpData.Rows(1,mindim);
						alldat += tmpData;	
					}				
					else 	
					message("Data dimensions do not match - ignoring "+opts.inputfname.value().at(ctr) << endl);
			}
   		}	

    	//update mask
    	if(opts.update_mask.value()){
      		message("Excluding voxels with constant value ...");
      		update_mask(Mask, alldat); 
      		message(" done" << endl);
    	}

		if((numfiles > 1 ) && opts.joined_vn.value() && tmpvarnorm){	
			//variance - normalisation
	    	message(endl<<"Normalising jointly by voxel-wise variance ..."); 
	    	stdDev = varnorm(alldat,alldat.Nrows(),opts.vn_level.value(),opts.econ.value());
	    	stdDevi = pow(stdDev,-1); 
	    	message(" done" << endl);
		}

		if(numfiles>1)
    		message(endl << "Initial data size : "<<alldat.Nrows()<<" x "<<alldat.Ncols()<<endl<<endl);

		if(opts.debug.value())
			save4D(alldat,"alldat");
    	//estimate model order
    	Matrix tmpPPCA;
    	RowVector AdjEV, PercEV;
    	Matrix tmpE;
	SymmetricMatrix Corr;
    	int order;

    	order = ppca_dim(remmean(alldat,2), RXweight, tmpPPCA, AdjEV, PercEV, Corr, pcaE, pcaD, Resels, opts.pca_est.value());	  
		if (opts.paradigmfname.value().length()>0)
			order += param.Ncols();

	  	if(opts.pca_dim.value() == 0){
      		opts.pca_dim.set_T(order);
			PPCA=tmpPPCA;
  		}
	  	if(opts.pca_dim.value() < 0){
      		opts.pca_dim.set_T(min(order,-1*opts.pca_dim.value()));
			PPCA=tmpPPCA;
  		}
    	order = opts.pca_dim.value();
		dbgmsg(endl << "Model order : "<<order<<endl);

		if (opts.paradigmfname.value().length()>0){
			Matrix tmpPscales;
			tmpPscales = param.t() * alldat;
			paramS = stdev(tmpPscales.t());

    		calc_white(pcaE, pcaD, order, param, paramS, whiteMatrix, dewhiteMatrix);
		}else
    		calc_white(pcaE, pcaD, order, whiteMatrix, dewhiteMatrix);

		if(opts.debug.value()){
			outMsize("pcaE",pcaE); saveascii(pcaE,"pcaE");
			outMsize("pcaD",pcaD); saveascii(pcaD,"pcaD");
			outMsize("AdjEV",AdjEV); saveascii(AdjEV,"AdjEV");
			outMsize("PercEV",PercEV); saveascii(PercEV,"PercEV");
			outMsize("tmpPPCA",tmpPPCA); saveascii(tmpPPCA,"tmpPPCA");
			outMsize("whiteMatrix",whiteMatrix); saveascii(whiteMatrix,"whiteMatrix");
			outMsize("dewhiteMatrix",dewhiteMatrix); saveascii(dewhiteMatrix,"dewhiteMatrix");
		}

		EV = AdjEV;
		EVP = PercEV;

    	if(numfiles == 1){
      		Data = alldat;
      		Matrix tmp = IdentityMatrix(Data.Nrows());
      		DWM.push_back(tmp);
      		WM.push_back(tmp);
    	} 
		else {

			dbgmsg("Multi-Subject ICA");
	        	//stdDev.CleanUp();
  			for(int ctr = 0; ctr < numfiles; ctr++){
				tmpData = process_file(opts.inputfname.value().at(ctr), numfiles);

				if(opts.joined_vn.value() && tmpvarnorm){
					dbgmsg("tmpData normalisation"<< endl);
					dbgmsg("stdDev "  << stdDev(1,2)<< endl);
					dbgmsg("tmpData " << tmpData.SubMatrix(1,1,1,2)<< endl);
					SP3(tmpData,pow(stdDev,-1));
				}
				//  whiten (separate / joint)
				Matrix newWM,newDWM; 
				if(!opts.joined_whiten.value()){	  
					message("  Individual whitening in a " << order << " dimensional subspace " << endl);
      	  			std_pca(tmpData, RXweight, Corr, pcaE, pcaD, opts.econ.value());
	  				calc_white(pcaE, pcaD, order, newWM, newDWM);
				}else{
					if(!opts.dr_pca.value()){
						std_pca(whiteMatrix*tmpData, RXweight, Corr, pcaE, pcaD, opts.econ.value());
						calc_white(pcaE, pcaD, order, newWM, newDWM);		
						newDWM=(dewhiteMatrix*newDWM);
						newWM=(newWM*whiteMatrix);
					}
					else{
					  if(opts.debug.value())
					    dbgmsg(" --mod_pca ");
						Matrix tmp1, tmp2;
						tmp1 = whiteMatrix * alldat;
						remmean(tmp1,2);
						tmp1 *= tmpData.t();
						tmp2 = pinv(tmp1.t()).t();  
						std_pca(tmp1 * tmpData, RXweight, Corr, pcaE, pcaD, opts.econ.value());
						calc_white(pcaE, pcaD, order, newWM, newDWM);		
						newDWM=(tmp2*newDWM);
						newWM=(newWM * tmp1);
					}
				}
				DWM.push_back(newDWM);
				WM.push_back(newWM);
				tmpData = newWM * tmpData;

				//concatenate Data
				if(Data.Storage() == 0)
	  			Data = tmpData;
				else
	  			Data &= tmpData;
      	    }
        }
    	opts.varnorm.set_T(tmpvarnorm);
   	    dbgmsg(string("END: setup_classic") << endl);	
 
  }

  void MelodicData::setup_migp()
  {
    dbgmsg(string("START: setup_migp") << endl);	
	
	std::vector<int> myctr;
	for (int i=0; i< numfiles ; ++i) myctr.push_back(i); 

	if(opts.migp_shuffle.value()){
		message("Randomising input file order" << endl);
		std::random_shuffle ( myctr.begin(), myctr.end() );
	}
	
	Matrix tmpData;
	bool tmpvarnorm = opts.varnorm.value();

	if(numfiles > 1 && opts.joined_vn.value()){
		opts.varnorm.set_T(false);
	}
		
	for(int ctr = 0; ctr < numfiles; ctr++){
		tmpData = process_file(opts.inputfname.value().at(myctr.at(ctr)), numfiles) / numfiles;
		
		if (opts.migpN.value()==0){
			opts.migpN.set_T(2*tmpData.Nrows()-1);
		}
		if(opts.debug.value())
			save4D(tmpData,string("preproc_dat") + num2str(ctr+1));
			
		if(Data.Storage()==0)
			Data = tmpData;
		else
  			Data &= tmpData;

		outMsize("Data", Data);
		//reduce dim down to manageable level
		if(Data.Nrows() > opts.migp_factor.value()*opts.migpN.value() || ctr==numfiles-1){
			message("  Reducing data matrix to a  " << opt.migpN.value() << " dimensional subspace " << endl);
			Matrix pcaE;
			SymmetricMatrix Corr;
			RowVector pcaD;
			std_pca(Data, RXweight, Corr, pcaE, pcaD, opts.econ.value());
		    pcaE = pcaE.Columns(pcaE.Ncols()-opts.migpN.value()+1,pcaE.Ncols());
		    Data = pcaE.t() * Data;	
		}
		outMsize("Data", Data);
		
    }

  	//update mask
    if(opts.update_mask.value()){
      message(endl<< "Excluding voxels with constant value ...");
      update_mask(Mask, Data); 
      message(" done" << endl);
    }

	Matrix tmp = IdentityMatrix(Data.Nrows());
	DWM.push_back(tmp);
	WM.push_back(tmp);
   	opts.varnorm.set_T(tmpvarnorm);

	if(opts.varnorm2.value()){
	  message("  Normalising by voxel-wise variance ..."); 
      stdDev = varnorm(Data,std::min(30,Data.Nrows()-1),
					opts.vn_level.value(), opts.econ.value());
	  message(" done" << endl);
	}
	
    dbgmsg(string("END: setup_migp") << endl);	
  }

  void MelodicData::setup()
  { 
	dbgmsg(string("START: setup") << endl);	
		
	numfiles = (int)opts.inputfname.value().size();
    setup_misc();
	if(opts.debug.value())
		memmsg(" after setup_misc ");
		
	if(opts.filtermode){ // basic setup for filtering only
		Data = process_file(opts.inputfname.value().at(0));
	}
	else{
		if((numfiles > 1) && (opts.approach.value()==string("defl") || opts.approach.value()==string("symm")))
	       opts.approach.set_T("concat");
		if(opts.migp.value())
			setup_migp();
		else
			setup_classic();   
    }

    message(endl << "  Data size : "<<Data.Nrows()<<" x "<<Data.Ncols()<<endl<<endl);
	outMsize("stdDev",stdDev);
	//meanC=mean(Data,2);
	if(opts.debug.value())
		save4D(Data,"concat_data");    
    //save the mean & mask
    save_volume(Mask,logger.appendDir("mask"));
    save_volume(Mean,logger.appendDir("mean"));

	dbgmsg(string("END: setup") << endl);	    
  } // void setup()
	
  void MelodicData::setup_misc()
  {
	dbgmsg(string("START: setup_misc") << endl);	

    //initialize Mean
    read_volume(Mean,opts.inputfname.value().at(0));

    //create mask
    create_mask(Mask);

	//setup background image
	if(opts.bgimage.value()>""){
		read_volume(background,opts.bgimage.value());
    	if(!samesize(Mean,background)){
        	cerr << "ERROR:: background image and data have different dimensions  \n\n";
        	exit(2);
    	}
	}else{
		background = Mean;
	}
		
    if(!samesize(Mean,Mask)){
      cerr << "ERROR:: mask and data have different dimensions  \n\n";
      exit(2);
    }

    //reset mean
    Mean *= 0;
     
    //set up weighting
    if(opts.segment.value().length()>0){
      create_RXweight();
    }

	//set up instacorr mask image
	if(opts.insta_fn.value()>""){
		dbgmsg(string(" Setting up instacorr mask") << endl);
		volume4D<float> tmp_im;
		read_volume4D(tmp_im,opts.insta_fn.value());
	
		if(!samesize(Mean,tmp_im[0])){
	        	cerr << "ERROR:: instacorr mask and data have different voxel dimensions  \n\n";
	        	exit(2);
		}	
		insta_mask = tmp_im.matrix(Mask); 
	}

    //seed the random number generator
    double tmptime = time(NULL);
    if ( opts.seed.value() != -1 ) {
      tmptime = opts.seed.value(); 
    }
    srand((unsigned int) tmptime);

	if(opts.paradigmfname.value().length()>0){
		message("  Use columns in " << opts.paradigmfname.value() 
	      << " for PCA initialisation" <<endl);
		param = read_ascii_matrix(opts.paradigmfname.value());
		
	    Matrix tmpPU, tmpPV;
		DiagonalMatrix tmpPD;
		SVD(param, tmpPD, tmpPU, tmpPV);
		param = tmpPU;
		
		opts.pca_dim.set_T(std::max(opts.pca_dim.value(), param.Ncols()+3));		
		if(opts.debug.value()){
			outMsize("Paradigm",param); saveascii(param,"param");
		}
		//opts.guessfname.set_T(opts.paradigmfname.value());
	}

	//read in post-proc design matrices etc
	if(opts.fn_Tdesign.value().length()>0)
		Tdes = read_ascii_matrix(opts.fn_Tdesign.value());
	if(opts.fn_Sdesign.value().length()>0)
		Sdes = read_ascii_matrix(opts.fn_Sdesign.value());
	if(opts.fn_Tcon.value().length()>0)
		Tcon = read_ascii_matrix(opts.fn_Tcon.value());
	if(opts.fn_Scon.value().length()>0)
		Scon = read_ascii_matrix(opts.fn_Scon.value());
	if(opts.fn_TconF.value().length()>0)
		TconF = read_ascii_matrix(opts.fn_TconF.value());
	if(opts.fn_SconF.value().length()>0)
		SconF = read_ascii_matrix(opts.fn_SconF.value());

  // Check that the number of input 
  // files matches the session design
  if (Sdes.Storage()) {
    if (Sdes.Nrows() != numfiles) {
      cerr << "ERROR: Number of input files (" << numfiles << ") " <<
        "does not match subject/session design (" << Sdes.Nrows() << ")" << endl;
      exit(2);
    }
  }			
  
  // Or create a default session design 
  // if one was not specified
  else if(numfiles>1){		
 		Sdes = ones(numfiles,1);
		if(Scon.Storage() == 0){
			Scon = ones(1,1);
			Scon &= -1*Scon;
		}
	}
	remmean(Tdes);
	
	dbgmsg(string("END: setup_misc") << endl);	
    
  }

  void MelodicData::save()
  {   

    //check for temporal ICA
    if(opts.temporal.value()){
      message(string("temporal ICA: transform back the data ... "));
      Matrix tmpIC = mixMatrix.t();
      mixMatrix=IC.t();
      IC=tmpIC;

      unmixMatrix=pinv(mixMatrix);
      Data = Data.t();
      tmpIC = meanC;
      meanC = meanR.t();
      meanR = tmpIC.t();
      //  whiteMatrix = whiteMatrix.t;
      //  dewhiteMatrix = dewhiteMatrix.t();
      message(string("done") << endl);
      opts.temporal.set_T(false); // Do not switch again!
    }
 
    message(endl << "Writing results to : " << endl);

    //Output IC	
    if((IC.Storage()>0)&&(opts.output_origIC.value())&&(after_mm==false))
      save4D(IC,opts.outputfname.value() + "_oIC");
      
    //Output IC -- adjusted for noise	
      if(IC.Storage()>0){
				volume4D<float> tempVol;	
   
				//Matrix ICadjust;
				if(after_mm){
	  			save4D(IC,opts.outputfname.value() + "_IC");
	  			// ICadjust = IC;
				}	
				else{
					
					Matrix resids = stdev(Data - mixMatrix * IC);
					for(int ctr=1;ctr<=resids.Ncols();ctr++)
						if(resids(1,ctr) < 0.05)
							resids(1,ctr)=1;
	  //			stdNoisei = pow(stdev(Data - mixMatrix * IC)*
		//				std::sqrt((float)(Data.Nrows()-1))/
		//				std::sqrt((float)(Data.Nrows()-IC.Nrows())),-1);
	  			stdNoisei = pow(resids*
						std::sqrt((float)(Data.Nrows()-1))/
						std::sqrt((float)(Data.Nrows()-IC.Nrows())),-1);
	  	
	  			ColumnVector diagvals;
	  			diagvals=pow(diag(unmixMatrix*unmixMatrix.t()),-0.5);
	
	  			save4D(SP(IC,diagvals*stdNoisei),opts.outputfname.value() + "_IC");
				}

				if(opts.output_origIC.value())
	  			save4D(stdNoisei,string("Noise__inv"));
      }
     
    //Output T- & S-modes
 
    save_Tmodes();
    save_Smodes();

    //Output mixMatrix
    if(mixMatrix.Storage()>0){
      saveascii(expand_mix(), opts.outputfname.value() + "_mix");
      mixFFT=calc_FFT(expand_mix(), opts.logPower.value());
      saveascii(mixFFT,opts.outputfname.value() + "_FTmix");      
    }

    //Output PPCA
    if(PPCA.Storage()>0)
      saveascii(PPCA, opts.outputfname.value() + "_PPCA");
  
    //Output ICstats
    if(ICstats.Storage()>0)
      saveascii(ICstats,opts.outputfname.value() + "_ICstats"); 
      
    //Output unmixMatrix
    if(opts.output_unmix.value() && unmixMatrix.Storage()>0)
      saveascii(unmixMatrix,opts.outputfname.value() + "_unmix");

    //Output Mask
    message("  "<< logger.appendDir("mask") <<endl);

    //Output mean
    if(opts.output_mean.value() && meanC.Storage()>0 && meanR.Storage()>0){
      saveascii(meanR,opts.outputfname.value() + "_meanR");
      saveascii(meanC,opts.outputfname.value() + "_meanC");
    }

    //Output white
    if(opts.output_white.value() && whiteMatrix.Storage()>0&&
      dewhiteMatrix.Storage()>0){
      	saveascii(whiteMatrix,opts.outputfname.value() + "_white");
      	saveascii(dewhiteMatrix,opts.outputfname.value() + "_dewhite");
      	Matrix tmp;
      	tmp=calc_FFT(dewhiteMatrix, opts.logPower.value());
      	saveascii(tmp,opts.outputfname.value() + "_FTdewhite");
	  }

    //Output PCA
    if(opts.output_pca.value() && pcaD.Storage()>0&&pcaE.Storage()>0){
      saveascii(pcaE,opts.outputfname.value() + "_pcaE");
      saveascii((Matrix) diag(pcaD),opts.outputfname.value() + "_pcaD");
      
      Matrix PCAmaps;
      if(whiteMatrix.Ncols()==Data.Ncols())
				PCAmaps = dewhiteMatrix.t();
      else
				PCAmaps = whiteMatrix * Data;

      save4D(PCAmaps,opts.outputfname.value() + "_pca");
    }
 
		message("...done" << endl);
  } //void save()

  int MelodicData::remove_components()
  {  
    message("Reading " << opts.filtermix.value() << endl) 
    mixMatrix = read_ascii_matrix(opts.filtermix.value());
    if (mixMatrix.Storage()<=0) {
      cerr <<" Please specify the mixing matrix correctly" << endl;
      exit(2);
    }
    
    unmixMatrix = pinv(mixMatrix);
    IC = unmixMatrix * Data;

    string tmpstr;
    tmpstr = opts.filter.value();

    Matrix noiseMix;
    Matrix noiseIC;

    int ctr=0;    
    char *p;
    char t[1024];
    const char *discard = ", [];{(})abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ~!@#$%^&*_-=+|\':><./?";
    
    message("Filtering the data...");
    strcpy(t, tmpstr.c_str());
    p=strtok(t,discard);
    ctr = atoi(p);
    if(ctr>0 && ctr<=mixMatrix.Ncols()){
      message(" "<< ctr );
      noiseMix = mixMatrix.Column(ctr);
      noiseIC  = IC.Row(ctr).t();    
    }
    else{
      cerr << endl<< "component number "<<ctr<<" does not exist" << endl;
    }
    
    do{
      p=strtok(NULL,discard);
      if(p){
				ctr = atoi(p);
	
        if(ctr>0 && ctr<=mixMatrix.Ncols()){
	  			message(" "<<ctr);
	  			noiseMix |= mixMatrix.Column(ctr);
	  			noiseIC  |= IC.Row(ctr).t();
				}
				else{
	  			cerr << endl<< "component number "<<ctr<<" does not exist" << endl;
				}
      }
    }while(p);
    message(endl);
    Matrix newData;

		outMsize("DATA",Data);
		outMsize("IC",IC);
		outMsize("noiseIC",noiseIC);
		outMsize("noiseMix",noiseMix);
		outMsize("meanR",meanR);
		outMsize("meanC",meanC);

    newData = Data - noiseMix * noiseIC.t();

		if(meanR.Storage()>0)
    	newData = newData + ones(newData.Nrows(),1)*meanR;
    
    volume4D<float> tmp;
    read_volume4D(tmp,opts.inputfname.value().at(0)); 
    tmp.setmatrix(newData,Mask);
    save_volume4D(tmp,logger.appendDir(opts.outputfname.value() + "_ICAfiltered")); 
   
    return 0;
  } // int remove_components()

  void MelodicData::create_RXweight()
  {
    message("Reading the weights for the covariance R_X from file "<< opts.segment.value() << endl);
  
    volume4D<float> tmpRX;
    read_volume4D(tmpRX,opts.segment.value());
    RXweight = tmpRX.matrix(Mask);
  } 

  void MelodicData::est_smoothness()
  {
    if(Resels == 0){
      string SM_path = opts.binpath + "smoothest";
      string Mask_fname = logger.appendDir("mask");

      if(opts.segment.value().length()>0){
				Mask_fname =  opts.segment.value();
      } 

      // Setup external call to smoothest:
      char callSMOOTHESTstr[1000];
      ostrstream osc(callSMOOTHESTstr,1000);
      osc  << SM_path << " -d " << data_dim()
	   		<< " -r " << opts.inputfname.value().at(0) << " -m " 
	   		<< Mask_fname << " > " << logger.appendDir("smoothest") << '\0';
      
      message("  Calling Smoothest: " << callSMOOTHESTstr << endl);
      system(callSMOOTHESTstr);

      //read back the results
      ifstream in;
      string str;
      Resels = 1.0;
      
      in.open(logger.appendDir("smoothest").c_str(), ios::in);
      if(in>0){
				for(int ctr=1; ctr<7; ctr++)
					in >> str;
				in.close();
				if(str!="nan")
	  			Resels = atof(str.c_str());
      }
    }
  }

  unsigned long MelodicData::standardise(volume<float>& mask, volume4D<float>& R)
  {
    	unsigned long count = 0;
    	int M=R.tsize();
    
    	for (int z=mask.minz(); z<=mask.maxz(); z++) {
      	for (int y=mask.miny(); y<=mask.maxy(); y++) {
					for (int x=mask.minx(); x<=mask.maxx(); x++) {
	  				if( mask(x,y,z) > 0.5) {
	    				count ++;
	    				if( M > 2 ) {
	      				// For each voxel 
	      				//    calculate mean and standard deviation...
	      				double Sx = 0.0, SSx = 0.0;	      
	      				for ( int t = 0; t < M; t++ ) {
									float R_it = R(x,y,z,t);
									Sx += R_it;
									SSx += (R_it)*(R_it);
	      				}
	      
	      				float mean = Sx / M;
	      				float sdsq = (SSx - ((Sx)*(Sx) / M)) / (M - 1) ;
	      
	      				if (sdsq<=0) {
									// trap for differences between mask and invalid data
									mask(x,y,z)=0;
									count--;
	      				} else {
	      					//    ... and use them to standardise to N(0, 1).
									for ( unsigned short t = 0; t < M; t++ ) {
		  							R(x,y,z,t) = (R(x,y,z,t) - mean) / sqrt(sdsq);
									}
	      				} 
	    				}
	  				}
					}
      	}
    	}  
    	return count;
  }

  float MelodicData::est_resels(volume4D<float> R, volume<float> mask)
  {
    message("  Estimating data smoothness ... ");
    unsigned long mask_volume = standardise(mask, R);

    int dof = R.tsize();
    unsigned long N = mask_volume;

    // MJ additions to make it cope with 2D images
    bool usez = true;
    if (R.zsize() <= 1) { usez = false; }

    enum {X = 0, Y, Z};
    double SSminus[3] = {0, 0, 0}, S2[3] = {0, 0, 0};

    int zstart=1;
    if (!usez) zstart=0;
    for ( unsigned short z = zstart; z < R.zsize() ; z++ )
      for ( unsigned short y = 1; y < R.ysize() ; y++ )
				for ( unsigned short x = 1; x < R.xsize() ; x++ )
	  			// Sum over N
	  			if( (mask(x, y, z)>0.5) &&
	      		(mask(x-1, y, z)>0.5) && 
	      		(mask(x, y-1, z)>0.5) && 
	      		( (!usez) || (mask(x, y, z-1)>0.5) ) ) {
	    
	    				N++;
	  
	    				for ( unsigned short t = 0; t < R.tsize(); t++ ) {
	      				// Sum over M
	      				SSminus[X] += R(x, y, z, t) * R(x-1, y, z, t);
	      				SSminus[Y] += R(x, y, z, t) * R(x, y-1, z, t);
	      				if (usez) SSminus[Z] += R(x, y, z, t) * R(x, y, z-1, t);

	      				S2[X] += 0.5 * (R(x, y, z, t)*R(x, y, z, t) + 
									R(x-1, y, z, t)*R(x-1, y, z, t));
	      				S2[Y] += 0.5 * (R(x, y, z, t)*R(x, y, z, t) + 
									R(x, y-1, z, t)*R(x, y-1, z, t));
	      				if (usez) S2[Z] += 0.5 * (R(x, y, z, t)*R(x, y, z, t) + 
									R(x, y, z-1, t)*R(x, y, z-1, t));
	    				}
	  			}

    double norm = 1.0/(double) N;
    double v = dof;	// v - degrees of freedom (nu)  
    if(R.tsize() > 1) {
      norm = (v - 2) / ((v - 1) * N * R.tsize());
    }
 
    // for extreme smoothness 
    if (SSminus[X]>=0.99999999*S2[X]) 
      SSminus[X]=0.99999*S2[X];  
    if (SSminus[Y]>=0.99999999*S2[Y]) 
      SSminus[Y]=0.99999*S2[Y];
    if (usez) 
      if (SSminus[Z]>=0.99999999*S2[Z]) 
				SSminus[Z]=0.99999*S2[Z];
    // Convert to sigma squared
    double sigmasq[3] = {0,0,0};
    sigmasq[X] = -1.0 / (4 * log(fabs(SSminus[X]/S2[X])));
    sigmasq[Y] = -1.0 / (4 * log(fabs(SSminus[Y]/S2[Y])));
    if (usez) { sigmasq[Z] = -1.0 / (4 * log(fabs(SSminus[Z]/S2[Z]))); }
    
    // Convert to full width half maximum
    double FWHM[3] = {0,0,0};
    FWHM[X] = sqrt(8 * log(2) * sigmasq[X]);
    FWHM[Y] = sqrt(8 * log(2) * sigmasq[Y]);
    if (usez) { FWHM[Z] = sqrt(8 * log(2) * sigmasq[Z]); }
    double resels = FWHM[X] * FWHM[Y];
    if (usez) resels *= FWHM[Z];

    message(" done " <<endl);
    return resels;
  }

  void MelodicData::create_mask(volume<float>& theMask)
  {
    if(opts.use_mask.value() && opts.maskfname.value().size()>0){   // mask provided 
      read_volume(theMask,opts.maskfname.value());
      message("Mask provided : " << opts.maskfname.value()<<endl<<endl);
    }
    else{
      if(opts.perf_bet.value() && opts.use_mask.value()){ //use BET
				message("Create mask ... ");
    		//save first image
    		tmpnam(Mean_fname); // generate a tmp name
    		save_volume(Mean,Mean_fname);    

				// set up all strings
				string BET_outputfname = string(Mean_fname)+"_brain";

				string BET_path = opts.binpath + "bet";
				string BET_optarg = "-m -f 0.4"; // see man bet
				string Mask_fname = BET_outputfname+"_mask";

				// Setup external call to BET:

		//		char callBETstr[1000];
	//			ostrstream betosc(callBETstr,1000);
//				betosc  << BET_path << " " << Mean_fname << " " 
//	     		<< BET_outputfname << " " << BET_optarg << " > /dev/null " << '\0';
	
//        message("  Calling BET: " << callBETstr << endl);
//				system(callBETstr);

				string tmpstr = BET_path + string(" ") + 
				                Mean_fname + string(" ") + BET_outputfname + string(" ") + 
				                BET_optarg + string(" > /dev/null ");
				system(tmpstr.c_str());
								
				// read back the Mask file   
				read_volume(theMask,Mask_fname);
				
				// clean /tmp
		    char callRMstr[1000];
		    ostrstream osc(callRMstr,1000);
		    osc  << "rm " << string(Mean_fname) <<"*  " << '\0';
		    system(callRMstr);
		   
				message("done" << endl);
      }  
      else{
				if(opts.use_mask.value()){   //just threshold the Mean
	  			message("Create mask ... ");
	  			float Mmin, Mmax, Mtmp;
	  			Mmin = Mean.min(); Mmax = Mean.max();
	  			theMask = binarise(Mean,Mmin + opts.threshold.value()* 
						(Mmax-Mmin),Mmax);
          Mtmp = Mmin + opts.threshold.value()* (Mmax-Mmin);
	  			message("done" << endl);
				}
				else{ //well, don't threshold then
	  		  theMask = Mean;
	  		  theMask = 1.0;
				}
      }
    }
    if(opts.remove_endslices.value()){ 
      // just in case mc introduced something nasty
      message("  Deleting end slices" << endl);
      for(int ctr1=theMask.miny(); ctr1<=theMask.maxy(); ctr1++){
				for(int ctr2=theMask.minx(); ctr2<=theMask.maxx(); ctr2++){   
	  			theMask(ctr2,ctr1,Mask.minz()) = 0.0;
	  			theMask(ctr2,ctr1,Mask.maxz()) = 0.0;
				}
      }
    }
  } //void create_mask()

  void MelodicData::sort()
  {
    int numComp = mixMatrix.Ncols(), numVox = IC.Ncols(), 
        numTime = mixMatrix.Nrows(), i,j;

	//flip IC maps to be positive (on max)
	//flip Subject/Session modes to be positive (on average)
	//flip time courses accordingly

    for(int ctr_i = 1; ctr_i <= numComp; ctr_i++)
      if(IC.Row(ctr_i).MaximumAbsoluteValue()>IC.Row(ctr_i).Maximum()){
				flipres(ctr_i);	
	  }
    message("Sorting IC maps" << endl);  
    Matrix tmpscales, tmpICrow, tmpMIXcol;
	if(numfiles > 1 && opts.approach.value()==string("tica")){
		set_TSmode();
		Matrix allmodes = Smodes.at(0);
		for(int ctr = 1; ctr < (int)Smodes.size();++ctr)
			allmodes |= Smodes.at(ctr);
		tmpscales = median(allmodes).t();	
	} else {
    // re-order wrt standard deviation of IC maps 
    tmpscales = stdev(IC,2);
	}
		
    ICstats = tmpscales;

    double max_val, min_val = tmpscales.Minimum()-1;

    for(int ctr_i = 1; ctr_i <= numComp; ctr_i++){
      max_val = tmpscales.Maximum2(i,j);
      ICstats(ctr_i,1)=max_val;
  
      tmpICrow = IC.Row(ctr_i);
      tmpMIXcol = mixMatrix.Column(ctr_i);
      
      IC.SubMatrix(ctr_i,ctr_i,1,numVox) = IC.SubMatrix(i,i,1,numVox);
      mixMatrix.SubMatrix(1,numTime,ctr_i,ctr_i) = 
				mixMatrix.SubMatrix(1,numTime,i,i);

      IC.SubMatrix(i,i,1,numVox) = tmpICrow.SubMatrix(1,1,1,numVox);
      mixMatrix.SubMatrix(1,numTime,i,i) = tmpMIXcol.SubMatrix(1,numTime,1,1);
  
      tmpscales(i,1)=tmpscales(ctr_i,1);
      tmpscales(ctr_i,1)=min_val;
    }

    ICstats /= ICstats.Column(1).Sum();
    ICstats *= 100;
    
    if(EVP.Storage()>0){
      tmpscales = ICstats.Column(1).AsMatrix(ICstats.Nrows(),1) * EVP(1,numComp);
      ICstats |= tmpscales;
    }

    if(Data.Storage()>0&&stdDev.Storage()>0){
 
      Matrix copeP(tmpscales), copeN(tmpscales);
      Matrix max_ICs(tmpscales), min_ICs(tmpscales);
			
      for(int ctr_i = 1; ctr_i <= numComp; ctr_i++){
				int i,j;
				max_ICs(ctr_i,1) = IC.Row(ctr_i).Maximum2(i,j);
				copeP(ctr_i,1) = std::abs((pinv(mixMatrix)*
					Data.Column(j)).Row(ctr_i).AsScalar()*stdDev(1,j)*100*
					(mixMatrix.Column(ctr_i).Maximum()-
					mixMatrix.Column(ctr_i).Minimum())/meanR(1,j));

				min_ICs(ctr_i,1) = IC.Row(ctr_i).Minimum2(i,j);
				copeN(ctr_i,1) = -1.0*std::abs((pinv(mixMatrix)*
					Data.Column(j)).Row(ctr_i).AsScalar()*stdDev(1,j)*100*
					(mixMatrix.Column(ctr_i).Maximum()-
					mixMatrix.Column(ctr_i).Minimum())/meanR(1,j));
      }
      ICstats |= copeP;
      ICstats |= copeN;
    }
			
    mixFFT=calc_FFT(expand_mix(), opts.logPower.value());
    unmixMatrix = pinv(mixMatrix);
  }

  void MelodicData::status(const string &txt)
  {
    cout << "MelodicData Object " << txt << endl;
    if(Data.Storage()>0){cout << "Data: " << Data.Nrows() <<"x" << Data.Ncols() << endl;}else{cout << "Data empty " <<endl;}
    if(pcaE.Storage()>0){cout << "pcaE: " << pcaE.Nrows() <<"x" << pcaE.Ncols() << endl;}else{cout << "pcaE empty " <<endl;}
    if(pcaD.Storage()>0){cout << "pcaD: " << pcaD.Nrows() <<"x" << pcaD.Ncols() << endl;}else{cout << "pcaD empty " <<endl;}
    if(whiteMatrix.Storage()>0){cout << "white: " << whiteMatrix.Nrows() <<"x" << whiteMatrix.Ncols() << endl;}else{cout << "white empty " <<endl;}
    if(dewhiteMatrix.Storage()>0){cout << "dewhite: " << dewhiteMatrix.Nrows() <<"x" << dewhiteMatrix.Ncols() << endl;}else{cout << "dewhite empty " <<endl;}
    if(mixMatrix.Storage()>0){cout << "mix: " << mixMatrix.Nrows() <<"x" << mixMatrix.Ncols() << endl;}else{cout << "mix empty " <<endl;}
    if(unmixMatrix.Storage()>0){cout << "unmix: " << unmixMatrix.Nrows() <<"x" << unmixMatrix.Ncols() << endl;}else{cout << "unmix empty " <<endl;}
    if(IC.Storage()>0){cout << "IC: " << IC.Nrows() <<"x" << IC.Ncols() << endl;}else{cout << "IC empty " <<endl;}
    
  } //void status()

}




