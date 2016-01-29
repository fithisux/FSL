/*  MELODIC - Multivariate exploratory linear optimized decomposition into 
              independent components
    
    meloptions.cc - class for command line options

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
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/log.h"
#include "meloptions.h"
#include "newimage/newimageall.h"

using namespace Utilities;
using namespace NEWIMAGE;

namespace Melodic {

MelodicOptions* MelodicOptions::gopt = NULL;

  void MelodicOptions::parse_command_line(int argc, char** argv, Log& logger, 
		const string &p_version){
  		//set version number and some other stuff
  		version = p_version;
  		filtermode = false;
  		explicitnums = false;
  		logfname = string("log.txt");

  		// work out the path to the $FSLDIR/bin directory
  		if(getenv("FSLDIR")!=0){
    		binpath = (string) getenv("FSLDIR") + "/bin/";
  		} else{
    		binpath = argv[0];
    		binpath = binpath.substr(0,binpath.length()-7);
  		}

  		// parse once to establish log directory name
  		for(int a = options.parse_command_line(argc, argv); a < argc; a++);

  		// act on simple command line arguments
  		if(help.value()){
      	print_usage(argc, argv);
      	exit(0);
    	} 
  		if(vers.value()){
      	print_version();
      	cout << endl;
      	exit(0);
    	} 
  		if(copyright.value()){
      	print_copyright();
      	exit(0);
    	} 
  		if(! options.check_compulsory_arguments()){
				print_usage(argc, argv);
      	exit(2);
    	}     

  		// check for invalid values
  		if (inputfname.value().size()<1) {
    		cerr << "ERROR:: Input volume not found\n\n";
    		print_usage(argc,argv);
    		exit(2);
  		}
  		if (approach.value() != "symm" && approach.value() != "defl"  && 
      	approach.value() != "jade" && approach.value() != "maxent" &&
      	approach.value() != "tica" && approach.value() != "concat"){
    			cerr << "ERROR:: unknown approach \n\n";
    			print_usage(argc,argv);
    			exit(2);
  			}
  		if (nonlinearity.value() != "pow3" && nonlinearity.value() != "pow4" && 
				nonlinearity.value() != "tanh"  && nonlinearity.value() != "gauss" ){
    			cerr << "ERROR:: unknown nonlinearity \n\n";
    			print_usage(argc,argv);
    			exit(2);
  			}
  		if (maxNumItt.value() < 1){
    		cerr << "ERROR:: maxItt too small \n\n";
    		print_usage(argc,argv);
    		exit(2);
  		}
  		if (epsilon.value() < 0.000000001){
    		cerr << "ERROR:: epsilon too small  \n\n";
    		print_usage(argc,argv);
    		exit(2);    
  		}
  		if (epsilon.value() >= 0.01){
    		cerr << "ERROR:: epsilon too large  \n\n";
    		print_usage(argc,argv);
    		exit(2);
  		}
  		if (nlconst1.value() <= 0){
    		cerr << "ERROR:: nlconst1 negative  \n\n";
    		print_usage(argc,argv);
    		exit(2);
  		}
  		if (!remove_meanvol.value()){
    		varnorm.set_T(false);
  		}
  		if (filter.set()){
    		if (filtermix.set()){
      		cerr << "ERROR:: no mixing matrix for filtering (use --mix='filename') \n\n"; 
      		print_usage(argc,argv);
      		exit(2);
    		} else {   
					temporal.set_T(false);
      		filtermode = true;
      		varnorm.set_T(false);
					pbsc.set_T(false);
					cerr << "WARNING: melodic denoising is deprecated, please use fsl_regfilt instead!" <<endl;
    		} 
  		}
  		if (output_all.value()){
    		output_unmix.set_T(true);
    		output_MMstats.set_T(true);
    		output_pca.set_T(true);
    		output_white.set_T(true);
    		output_origIC.set_T(true);
    		output_mean.set_T(true);
  		}
  		if (output_pca.value()){
    		output_white.set_T(true);
  		}
  		if(threshold.value()>=1){
    		threshold.set_T(threshold.value()/100);
    		if(threshold.value()>=1){
      		cerr << "ERROR:: threshold level not a percentage value  \n\n";
      		print_usage(argc,argv);
      		exit(2);
    		}
  		}
  		if (nlconst2.value() <= 0){
    		cerr << "ERROR:: nlconst2 negative  \n\n";
    		print_usage(argc,argv);
    		exit(2);
  		}
  		if (dummy.value() < 0){
    		cerr << "ERROR:: negative dummy value  \n\n";
    		print_usage(argc,argv);
    		exit(2);
  		}
  		if (repeats.value() < 1){
    		cerr << "ERROR:: repeats < 1 \n\n";
    		print_usage(argc,argv);
    		exit(2);
  		}
  		if (numICs.value() > 0){
    		explicitnums = true;
  		}
		if (insta_fn.value() > ""){
			varnorm.set_T(false);
		}
  
  		//in the case of indirect inputs, create the vector of input names here
  		if(!fsl_imageexists(inputfname.value().at(0))){
    		std::vector< string > tmpfnames;
    		ifstream fs(inputfname.value().at(0).c_str());
    		string cline;
    		while (!fs.eof()) {
		  getline(fs,cline);
		  if(cline.length()>0) {
		    while ( cline.find(' ') != cline.npos )
		      cline.erase( cline.find(' '),1);
		    tmpfnames.push_back(cline);
		  }
    		}	
    		fs.close();
    		inputfname.set_T(tmpfnames);
  		}

  		//transform inputfnames to their basenames
  		std::vector< string > tmpfnames;
  		for(int ctr=0; ctr < (int)inputfname.value().size() ; ctr++){
    		string basename;
    		basename = inputfname.value().at(ctr);
    		make_basename(basename);
    		tmpfnames.push_back(basename);
  		}
  		inputfname.set_T(tmpfnames);

  		//create melodic directory name
  		if(logdir.value()==string("log.txt")){
    		logdir.set_T(string(inputfname.value().at(0)+".ica"));
    		logger.makeDir(logdir.value(),logfname);
 			} else{
    		// setup logger directory
    		system(("mkdir "+ logdir.value() + " 2>/dev/null").c_str());
    		logger.setDir(logdir.value(),logfname);
  		}
  		message(endl << "Melodic Version " << version << endl << endl);

  		// parse again so that options are logged
  		for(int a = 0; a < argc; a++)
    		logger.str() << argv[a] << " ";
  			logger.str() << endl << "---------------------------------------------" 
					<< endl << endl;
  			message("Melodic results will be in " << logger.getDir() << endl << endl);
    }

void MelodicOptions::print_usage(int argc, char *argv[]){
 // print_copyright();
  options.usage();
  /* cout << "Usage: " << argv[0] << " ";
    Have own usage output here
  */
}

void MelodicOptions::print_version(){
	cout << endl <<"MELODIC Version " << version << endl;
  cout.flush();
}

void MelodicOptions::print_copyright(){
  cout << endl << title << endl;
  cout.flush();
}

void MelodicOptions::status(){
  cout << " version = " << version << endl;
 
  cout << " logdir = "  << logdir.value() << endl;
  cout << " inputfname = "  << inputfname.value().at(0) << inputfname.value().at(inputfname.value().size()-1) << endl;
  cout << " outputfname = "  << outputfname.value() << endl;
  cout << " guessfname = "  << guessfname.value() << endl;
  cout << " maskfname = "  << maskfname.value() << endl;
  cout << " paradigmfname = "  <<  paradigmfname.value() << endl;
  cout << " nonlinearity = "  << nonlinearity.value() << endl;
  cout << " approach = "  << approach.value() << endl;
  
  cout << " pca_dim = "  << pca_dim.value() << endl;
  cout << " segment = "  << segment.value() << endl;
  cout << " numICs = "  << numICs.value() << endl;
  cout << " maxNumItt = "  << maxNumItt.value() << endl;
  cout << " maxRestart = "  << maxRestart.value() << endl;
  cout << " dummy = "  << dummy.value() << endl;
  cout << " repeats = "  << repeats.value() << endl;
  cout << " nlconst1 = "  << nlconst1.value() << endl;
  cout << " nlconst2 = "  << nlconst2.value() << endl;
  cout << " smooth_probmap = "  << smooth_probmap.value() << endl;
  cout << " epsilon = "  << epsilon.value() << endl;
  cout << " threshold = "  << threshold.value() << endl;
  
  cout << " help = "  << help.value() << endl;
  cout << " verbose = "  << verbose.value() << endl;
  cout << " vers = "  << vers.value() << endl;
  cout << " copyright = "  << copyright.value() << endl;
  cout << " perf_bet = "  << perf_bet.value() << endl;
  cout << " remove_meanvol = "  << remove_meanvol.value() << endl;
  cout << " remove_endslices  = " << remove_endslices.value() << endl;
  cout << " output_pca = "  << output_pca.value() << endl;
  cout << " output_white = "  << output_white.value() << endl;
  cout << " output_mean = "  << output_mean.value() << endl;
  cout << " use_mask = "  << use_mask.value() << endl;
  cout << " output_unmix = "  << output_unmix.value() << endl;
  cout << " guess_remderiv = "  << guess_remderiv.value() << endl;
  cout << " filter  = " << filter.value() << endl;
  cout << " logPower  = " << logPower.value() << endl;
}

}
