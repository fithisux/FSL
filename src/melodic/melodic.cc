/*  MELODIC - Multivariate exploratory linear optimized decomposition into 
              independent components
    
    melodic.cc - main program file

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

#include <iostream>
#include <iomanip>
#include "newmatap.h"
#include "newmatio.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include <string>
#include <math.h>
#include "utils/options.h"
#include "utils/log.h"
#include "meloptions.h"
#include "meldata.h"
#include "melpca.h"
#include "melica.h"
#include "melodic.h"
#include "melreport.h"
#include "melhlprfns.h"
#include "melgmix.h"

using namespace Utilities;
using namespace NEWMAT;
using namespace NEWIMAGE;
using namespace Melodic;
using namespace MISCPLOT;

string myfloat2str(float f, int width, int prec, bool scientif){
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

Matrix mmall(Log& logger, MelodicOptions& opts,
	     MelodicData& melodat, MelodicReport& report, Matrix& probs);

void mmonly(Log& logger, MelodicOptions& opts,
	    MelodicData& melodat, MelodicReport& report);

int main(int argc, char *argv[]){
  try{
    // Setup logging:
    Log& logger  =   LogSingleton::getInstance();
 
    // parse command line - will output arguments to logfile
    MelodicOptions& opts = MelodicOptions::getInstance();
    opts.parse_command_line(argc, argv, logger, Melodic::version); 

    //set up data object
    MelodicData melodat(opts,logger);
    //set up report object     
    MelodicReport report(melodat,opts,logger);
   
    if (opts.filtermode || opts.filtermix.value().length()>0 || opts.ICsfname.value().length()>0){
      if(opts.filtermode){ // just filter out some noise from a previous run
    		melodat.setup();
				if(opts.debug.value())
					message("  Denoising data setup completed "<< endl);
				melodat.remove_components();
      } else
				mmonly(logger,opts,melodat,report);
    } else {  // standard PICA now
      int retry = 0;
      bool no_conv;
      bool leaveloop = false;

      melodat.setup();

	  if (opts.maxRestart.value()<0)
		opts.maxRestart.set_T(melodat.data_dim());
 
      do{
				//do PCA pre-processing
				MelodicPCA pcaobj(melodat,opts,logger,report);
				pcaobj.perf_pca();
				pcaobj.perf_white();

				//do ICA
				MelodicICA icaobj(melodat,opts,logger,report);
				icaobj.perf_ica(melodat.get_white()*melodat.get_Data());
    
				no_conv = icaobj.no_convergence;

				if(no_conv){
				    retry++;
					if((opts.approach.value()=="symm")&&(retry == opts.maxRestart.value())){
						// try final round with defl
						opts.approach.set_T("defl");
					    message(endl << "Restarting MELODIC using deflation approach" << endl << endl);	
					}
					else{
					    // try using different dim	
						if((int)opts.pca_dim.value()*opts.retryfactor.value() > (int)(0.05*melodat.data_dim()+1)){
		      				opts.pca_dim.set_T((int)opts.pca_dim.value()*opts.retryfactor.value());
		    			}
     					else{
							if((int)opts.pca_dim.value()/opts.retryfactor.value() > (int)(melodat.data_dim())){
			      				opts.pca_dim.set_T((int)opts.pca_dim.value()/opts.retryfactor.value());
			    			}
							else{
								leaveloop = TRUE;
							}
						}	
					}				
				}
      } while (no_conv && retry<opts.maxRestart.value() && !leaveloop);	
     
      if(!no_conv){
				//save raw IC results
				melodat.save();
				Matrix pmaps;//(melodat.get_IC());
				Matrix mmres;
			
				message("Creating report index page ...");
				if( bool(opts.genreport.value()) ){
		  		report.analysistxt();
					if(melodat.get_numfiles()>1)
						report.Smode_rep();
		  		report.PPCA_rep();
				}
					
				message("done"<< endl <<endl);
				if(opts.perf_mm.value())
	  			mmres = mmall(logger,opts,melodat,report,pmaps);
				else{
	  			if( bool(opts.genreport.value()) ){
	    			message(endl 
		    			<< "Creating web report in " << report.getDir() 
		    			<< " " << endl);
	    			for(int ctr=1; ctr<= melodat.get_IC().Nrows(); ctr++){
	      			string prefix = "IC_"+num2str(ctr);
	      			message("  " << ctr);
	      			report.IC_simplerep(prefix,ctr,melodat.get_IC().Nrows());
	    			}

	    			message(endl << endl <<
		    			" To view the output report point your web browser at " <<
		    			report.getDir() + "/00index.html" << endl<< endl); 
	  			}
				}		 

				message("finished!" << endl << endl);
      } 
	  else { 
				message(endl <<"No convergence -- giving up " << endl << endl);
      }	     
    }
  }
  catch(Exception e) {
      cerr << endl << e.what() << endl;
  }
  catch(X_OptionError& e) {
      cerr << endl << e.what() << endl;
  }

  return 0;
}

void mmonly(Log& logger, MelodicOptions& opts,
	MelodicData& melodat, MelodicReport& report){

  Matrix ICs;
  Matrix mixMatrix;
  Matrix fmixMatrix;
  volume<float> Mask;
  volume<float> Mean;
  
  {
    volume4D<float> RawData;
    message("Reading data file " << opts.inputfname.value().at(0) << "  ... ");
    read_volume4D(RawData,opts.inputfname.value().at(0));
    message(" done" << endl);
    Mean = meanvol(RawData);
  }

  {
    volume4D<float> RawIC;
    message("Reading components " << opts.ICsfname.value() << "  ... ");
    read_volume4D(RawIC,opts.ICsfname.value());
    message(" done" << endl);

    message("Creating mask   ... ");
    Mask = binarise(RawIC[0],float(RawIC[0].min()),float(RawIC[0].max()));

    ICs = RawIC.matrix(Mask);
    if(ICs.Nrows()>1){
      Matrix DStDev=stdev(ICs);
      
      volume4D<float> tmpMask;
      tmpMask.setmatrix(DStDev,Mask);
      
      float tMmax;
      volume<float> tmpMask2;
      tmpMask2 = tmpMask[0];
      tMmax = tmpMask2.max();
      double st_mean = DStDev.Sum()/DStDev.Ncols();
      double st_std  = stdev(DStDev.t()).AsScalar();
      
      Mask = binarise(tmpMask2,(float) max((float) st_mean-3*st_std,
					   (float) 0.01*st_mean),tMmax);  
      ICs = RawIC.matrix(Mask);
    }
    else{
      Mask = binarise(RawIC[0],float(0.001),float(RawIC[0].max())) 
	+ binarise(RawIC[0],float(RawIC[0].min()),float(-0.001));
      ICs = RawIC.matrix(Mask);
    }

    //cerr << "ICs : " << ICs.Ncols() << ICs.Nrows() << endl;
    message(" done" << endl);
  }

  if(opts.filtermix.value().length() > 0){
    message("Reading mixing matrix " << opts.filtermix.value() << " ... ");
    mixMatrix = read_ascii_matrix(opts.filtermix.value());
    if (mixMatrix.Storage()<=0) {
      cerr <<" Please specify the mixing matrix correctly" << endl;
      exit(2);
    }
    message(" done" << endl);
  }else{
	mixMatrix=unifrnd(ICs.Nrows()+1,ICs.Nrows());	
  }


  if(opts.smodename.value().length() > 0){
    message("Reading matrix of subject modes: " << opts.smodename.value());
    Matrix tmp;
    tmp = read_ascii_matrix(opts.smodename.value());
    if (tmp.Storage()<=0) {
      cerr <<" Please specify the mixing matrix correctly" << endl;
      exit(2);
    }
    message(" done" << endl);
    for (int ctr = 1; ctr <= tmp.Ncols(); ctr++){
      Matrix tmp2 = tmp.Column(ctr);
      melodat.add_Smodes(tmp2);
    }
  }

  melodat.set_mask(Mask);
  melodat.set_mean(Mean);
  melodat.set_IC(ICs);
  melodat.set_mix(mixMatrix);
  fmixMatrix = calc_FFT(mixMatrix, opts.logPower.value());
  melodat.set_fmix(fmixMatrix);
  fmixMatrix = pinv(mixMatrix);
  melodat.set_unmix(fmixMatrix);
  //  melodat.sort();
  //  write_ascii_matrix("ICs",ICs);
  
  Matrix mmres;
  Matrix pmaps;//(ICs);
 if(opts.perf_mm.value())
    mmres = mmall(logger,opts,melodat,report,pmaps);
	}

Matrix mmall(Log& logger, MelodicOptions& opts,MelodicData& melodat, MelodicReport& report, Matrix& pmaps){
  
  Matrix mmpars(5*melodat.get_IC().Nrows(),5);
  mmpars = 0;
  Log stats;
  
  if(opts.output_MMstats.value()){  
    stats.makeDir(logger.appendDir("stats"),"stats.log");
  }

  message(endl 
	  << "Running Mixture Modelling on Z-transformed IC maps ..." 
	  << endl);
  
  ColumnVector diagvals;
  diagvals=pow(diag(melodat.get_unmix()*melodat.get_unmix().t()),-0.5);
  
  for(int ctr=1; ctr <= melodat.get_IC().Nrows(); ctr++){
    MelGMix mixmod(opts, logger);
    
    message("  IC map " << ctr << " ... "<< endl;);
    
    Matrix ICmap;
    if(melodat.get_stdNoisei().Storage()>0)
 		dbgmsg(" stdNoisei max : "<< melodat.get_stdNoisei().Maximum() <<" "<< melodat.get_stdNoisei().Minimum() << endl);

    if(opts.varnorm.value()&&melodat.get_stdNoisei().Storage()>0){
      ICmap = SP(melodat.get_IC().Row(ctr),diagvals(ctr)*melodat.get_stdNoisei());
    }else{
      ICmap = melodat.get_IC().Row(ctr);
	}
    string wherelog;
    if(opts.genreport.value())
      wherelog = report.getDir();
    else
      wherelog = logger.getDir();

		dbgmsg(" ICmap max : "<< mean(ICmap,2).AsScalar() << endl);
    mixmod.setup( ICmap,
		  wherelog,ctr,
		  melodat.get_mask(), 
		  melodat.get_mean(),3);
    message("   calculating mixture-model fit "<<endl);
    mixmod.fit("GGM");

    if(opts.output_MMstats.value()){
      message("   saving probability map:  ");
      melodat.save4D(mixmod.get_probmap(),
		     string("stats/probmap_")+num2str(ctr));
	  message("   saving mixture model fit:");
	  melodat.saveascii(mixmod.get_params(),
			 string("stats/MMstats_")+num2str(ctr));  
    }

    //re-scale spatial maps to mean 0 for nht
    if(opts.rescale_nht.value()){
      message("   re-scaling spatial maps ... "<< endl); 
      RowVector tmp;
      tmp = mixmod.get_means();
      float dmean  = tmp(1);
      tmp = mixmod.get_vars();
      float dstdev = sqrt(tmp(1));

      tmp = (mixmod.get_means() - dmean)/dstdev;
      mixmod.set_means(tmp);
      tmp = (mixmod.get_vars() / (dstdev*dstdev));
      mixmod.set_vars(tmp);

      //tmp = (mixmod.get_data() - dmean)/dstdev;
      tmp = (ICmap - dmean)/dstdev;
      mixmod.set_data(tmp);
      //if(opts.varnorm.value()&&melodat.get_stdNoisei().Storage()>0)
      //	tmp = SP(tmp,pow(diagvals(ctr)*melodat.get_stdNoisei(),-1));
     
      melodat.set_IC(ctr,tmp);
    }

    if(opts.smooth_probmap.value()<0.0){
      message("   smoothing probability map ... "<< endl);
      mixmod.smooth_probs(0.5*(std::min(std::min(std::abs(melodat.get_mean().xdim()),std::abs(melodat.get_mean().ydim())),std::abs(melodat.get_mean().zdim()))));  
    }

    if(opts.smooth_probmap.value()>0.0){
      message("   smoothing probability map ... "<< endl);
      mixmod.smooth_probs(opts.smooth_probmap.value());
    }

    message("   thresholding ... "<< endl);
    mixmod.threshold(opts.mmthresh.value());  

    Matrix tmp;
    tmp=(mixmod.get_threshmaps().Row(1));
    float posint = SP(tmp,gt(tmp,zeros(1,tmp.Ncols()))).Sum();
    float negint = -SP(tmp,lt(tmp,zeros(1,tmp.Ncols()))).Sum();
    
    if((posint<0.01)&&(negint<0.01)){
      mixmod.clear_infstr();
      mixmod.threshold("0.05n "+opts.mmthresh.value());
      posint = SP(tmp,gt(tmp,zeros(1,tmp.Ncols()))).Sum();
      negint = -SP(tmp,lt(tmp,zeros(1,tmp.Ncols()))).Sum();
    }
    if(negint>posint){//flip map
      //  melodat.flipres(ctr);
      //  mixmod.flipres(ctr);
    }

    //save mixture model stats 
    if(opts.output_MMstats.value()){
      stats << " IC " << num2str(ctr) << " " << mixmod.get_type() << endl
	    << " Means :  " << mixmod.get_means() << endl
	    << " Vars. :  " << mixmod.get_vars()  << endl
	    << " Prop. :  " << mixmod.get_pi()    << endl << endl;
      message("   saving thresholded Z-stats image:");
      melodat.save4D(mixmod.get_threshmaps(),
		     string("stats/thresh_zstat")+num2str(ctr));
    }

    //save mmpars
    // mmpars((ctr-1)*5+1,1) = ctr;
		//     if(mixmod.get_type()=="GGM")
		//       mmpars((ctr-1)*5+1,2) = 1.0;
		//     else
		//       mmpars((ctr-1)*5+1,2) = 0.0;
		//     mmpars((ctr-1)*5+1,2) = mixmod.get_means().Ncols();
		//     tmp =  mixmod.get_means();
		//     for(int ctr2=1;ctr2<=mixmod.get_means().Ncols();ctr2++)
		//       mmpars((ctr-1)*5+2,ctr2) = tmp(1,ctr2);
		//     tmp =  mixmod.get_vars();
		//     for(int ctr2=1;ctr2<=mixmod.get_vars().Ncols();ctr2++)
		//       mmpars((ctr-1)*5+3,ctr2) = tmp(1,ctr2);
		//     tmp =  mixmod.get_pi(); 
		//     for(int ctr2=1;ctr2<=mixmod.get_pi().Ncols();ctr2++)
		//       mmpars((ctr-1)*5+4,ctr2) = tmp(1,ctr2);
		//     mmpars((ctr-1)*5+5,1) = mixmod.get_offset();

    if( bool(opts.genreport.value()) ){
      message("   creating report page ... ");
      report.IC_rep(mixmod,ctr,melodat.get_IC().Nrows(),melodat.get_ICstats());	    
      message("   done" << endl);
    }
  }
  
  if(!opts.filtermode&&opts.ICsfname.value().length()==0){
    //now safe new data
    //    bool what = opts.verbose.value();
    //opts.verbose.set_T(false);
    
    melodat.set_after_mm(TRUE);
    melodat.save();
    //opts.verbose.set_T(what);

    //if(melodat.get_IC().Storage()>0){
    //  volume4D<float> tempVol;	
    //  tempVol.setmatrix(melodat.get_IC(),melodat.get_mask());
    //  save_volume4D(tempVol,logger.appendDir(opts.outputfname.value() 
  	//					     + "_IC"),melodat.tempInfo);
	  //  message(endl<< endl << " Saving " << logger.appendDir(opts.outputfname.value() + "_IC") <<endl);
   	//}

   if( bool(opts.genreport.value()) ){
    message(endl << endl << 
	    " To view the output report point your web browser at " <<
	    report.getDir() + "/00index.html" << endl << endl); 
   }   
  }
  return mmpars;
  }
