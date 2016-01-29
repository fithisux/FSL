/*  gsoptions.h

    Mark Woolrich, Tim Behrens - FMRIB Image Analysis Group

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

#if !defined(GsOptions_h)
#define GsOptions_h

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "utils/log.h"

using namespace Utilities;

namespace Gs {

  class GsOptions {
  public:
    static GsOptions& getInstance();
    ~GsOptions() { delete gopt; }
  
    Option<bool> verbose;
    Option<int> debuglevel;
    Option<bool> timingon;
    Option<bool> help;
    Option<string> copefile;
    Option<string> varcopefile;
    Option<string> dofvarcopefile;
    Option<string> maskfile;
    Option<string> designfile;
    Option<string> covsplitfile;
    Option<string> tcontrastsfile;
    Option<string> fcontrastsfile;
    Option<string> logdir;
    Option<int> njumps;
    Option<int> burnin;
    Option<int> sampleevery;
    Option<float> zlowerthreshold;
    Option<float> zupperthreshold;
    Option<bool> fixmeanfortfit;
    Option<bool> infer_outliers;
    Option<bool> no_pe_output;
    Option<string> runmode;
    Option<int> seed;
    Option<vector<int> > voxelwise_ev_numbers;
    Option<vector<string> > voxelwise_ev_filenames;
    Option<float> sigma_smooth_flame2_dofs;
    Option<float> sigma_smooth_globalproboutlier;
    Option<int> io_niters;
    Option<bool> outputDof;
  
    void parse_command_line(int argc, char** argv, Log& logger);
  
  private:
    GsOptions();  
    const GsOptions& operator=(GsOptions&);
    GsOptions(GsOptions&);

    OptionParser options; 
      
    static GsOptions* gopt;
  
  };

  inline GsOptions& GsOptions::getInstance(){
    if(gopt == NULL)
      gopt = new GsOptions();
   
    return *gopt;
  }

  inline GsOptions::GsOptions() :
    verbose(string("-v,-V,--verbose"), false, 
	    string("switch on diagnostic messages"), 
	    false, no_argument),
    debuglevel(string("--debug,--debuglevel"), 0,
	       string("set debug level"), 
	       false, requires_argument),
    timingon(string("--to,--timingon"), false, 
	     string("turn timing on"), 
	     false, no_argument),
    help(string("-h,--help"), false,
	 string("display this message"),
	 false, no_argument),
    copefile(string("--cope,--copefile"), string("cope"),
	     string("cope regressor data file"),
	     true, requires_argument),  
    varcopefile(string("--vc,--vcope,--varcope,--varcopefile"), string(""),
		string("varcope weightings data file"),
		false, requires_argument),
    dofvarcopefile(string("--dvc,--dvcope,--dofvarcope,--dofvarcopefile"), string(""),
		   string("DOF data file for the varcope data"),
		   false, requires_argument),
    maskfile(string("--mask,--maskfile"), string(""),
	     string("mask file"),
	     true, requires_argument),
    designfile(string("--dm,--designfile"), string("design.mat"),
	       string("design matrix file"),
	       true, requires_argument), 
    covsplitfile(string("--cs,--csf,--covsplitfile"), string("design.gs"),
		 string("file containing an ASCII matrix specifying the groups the covariance is split into"),
		 true, requires_argument),
    tcontrastsfile(string("--tc,--tcf,--tcontrastsfile"), string("design.con"),
		   string("file containing an ASCII matrix specifying the t contrasts"),
		   true, requires_argument),  
    fcontrastsfile(string("--fc,--fcf,--fcontrastsfile"), string(""),
		   string("file containing an ASCII matrix specifying the f contrasts"),
		   false, requires_argument),  
    logdir(string("--ld,--logdir"), string("logdir"),
	   string("log directory"),
	   false, requires_argument),  
    njumps(string("--nj,--njumps"), 5000,
	   string("Num of jumps to be made by MCMC"),
	   false, requires_argument),
    burnin(string("--bi,--burnin"), 500,
	   string("Num of jumps at start of MCMC to be discarded"),
	   false, requires_argument),
    sampleevery(string("--se,--sampleevery"), 1,
		string("Num of jumps for each sample"),
		false, requires_argument),
    zlowerthreshold(string("--zlt"), 0,
		    string("Absolute value of lower threshold on FLAME stage 1 z-stats, used to see which voxels are passed to FLAME stage 2"),
		    false, requires_argument),
    zupperthreshold(string("--zut"), 0,
		    string("Absolute value of upper threshold on FLAME stage 1 z-stats, used to see which voxels are passed to FLAME stage 2"),
		    false, requires_argument),
    fixmeanfortfit(string("--fm,--fixmean"), false, 
		   string("Fix mean for tfit"), 
		   false, no_argument),
    infer_outliers(string("--io,--inferoutliers"), false, 
		   string("Infer outliers -- note that does not apply to fixed effects"), 
		   false, no_argument),
    no_pe_output(string("--npo,--nopeoutput"), false, 
		   string("Do not output PE files"), 
		   false, no_argument),
    runmode(string("--runmode"), string(""), 
	    string("Inference to perform: runmode=fe (fixed effects), runmode=ols (mixed effects - OLS), runmode=flame1 (mixed effects - FLAME stage 1), runmode=flame12 (mixed effects - FLAME stage 1+2)"), 
	    true, requires_argument),
    seed(string("--seed"), 10, 
	 string("seed for pseudo random number generator"), 
	 false, requires_argument),
    voxelwise_ev_numbers(string("--voxelwise_ev_numbers,--ven"), vector<int>(), 
	 string("List of numbers indicating voxelwise EVs position in the design matrix (list corresponds in order to files in voxelwise_ev_filenames)"), 
	 false, requires_argument),
    voxelwise_ev_filenames(string("--voxelwise_ev_filenames,--vef"), vector<string>(), 
	 string("List of 4D niftii files containing voxelwise EVs (list corresponds in order to numbers in voxelwise_ev_numbers)"), 
	 false, requires_argument),
    sigma_smooth_flame2_dofs(string("--sdof,--sigma_dofs"), 1,
				   string("sigma (in mm) to use for Gaussian smoothing the DOFs in FLAME 2. Default is 1mm, -1 indicates no smoothing"),
	 false, requires_argument),
    sigma_smooth_globalproboutlier(string("--so"), -1,
				   string(""),
	 false, requires_argument),
    io_niters(string("--ioni"), 12,
				   string("Number of max iterations to use when inferring outliers. Default is 12."),
	 false, requires_argument),
    outputDof(string("--outputdof"), false, 
	    string("output dof for lower-level compatibility usage"), 
	    false, no_argument),
    options("flameo"," flameo --cope=filtered_func_data --mask=mask --dm=design.mat --tc=design.con --cs=design.grp --runmode=ols\n flameo --cope=filtered_func_data --varcope=var_filtered_func_data --mask=mask --dm=design.mat --tc=design.con --cs=design.grp --runmode=flame1")
  {
    try {
      options.add(verbose);
      options.add(debuglevel);
      options.add(timingon);
      options.add(help);
      options.add(copefile);
      options.add(varcopefile);
      options.add(dofvarcopefile);
      options.add(maskfile);
      options.add(designfile);
      options.add(tcontrastsfile);
      options.add(fcontrastsfile);
      options.add(covsplitfile);
      options.add(logdir);
      options.add(njumps);
      options.add(burnin);
      options.add(sampleevery);
      options.add(zlowerthreshold);
      options.add(zupperthreshold);
      options.add(fixmeanfortfit);
      options.add(infer_outliers);
      options.add(no_pe_output);
      options.add(runmode);
      options.add(seed);
      options.add(voxelwise_ev_numbers);
      options.add(voxelwise_ev_filenames);
      options.add(sigma_smooth_flame2_dofs);
      options.add(sigma_smooth_globalproboutlier);
      options.add(io_niters);
      options.add(outputDof);
    }
    catch(X_OptionError& e) {
      options.usage();
      cerr << endl << e.what() << endl;
    } 
    catch(std::exception &e) {
      cerr << e.what() << endl;
    }    
     
  }
}

#endif



