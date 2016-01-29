 /*  MELODIC - Multivariate exploratory linear optimized decomposition into 
              independent components
    
    meloptions.h - class for command line options

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

#ifndef __MELODICOPTIONS_h
#define __MELODICOPTIONS_h

#include <string>
#include <strstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "utils/log.h"
#include "melodic.h"

using namespace Utilities;

namespace Melodic {

class MelodicOptions {
	public:
  	static MelodicOptions& getInstance();
  	~MelodicOptions() { delete gopt; }
  
  	string version;
  	string binpath;
  	string logfname;
  	bool   filtermode;
  	bool   explicitnums;
  
  	Option<string> logdir;
  	Option< std::vector<string> > inputfname;

  	Option<string> outputfname;

  	Option<string> maskfname;
  	Option<bool>   use_mask;
  	Option<bool>   update_mask;
  	Option<bool>   perf_bet;
  	Option<float>  threshold;

  	Option<int>    pca_dim;
  	Option<string> pca_est;
  	Option<bool>   joined_whiten;
  	Option<bool>   joined_vn;
	Option<bool>   dr_pca;
	Option<bool>   migp;
	Option<int>    migpN;
	Option<bool>   migp_shuffle;
	Option<int>	   migp_factor;
	Option<bool>   dr;
	Option<bool>   dr_out;
	Option<float>  vn_level;
  	Option<int>    numICs;
  	Option<string> approach;
  	Option<string> nonlinearity;

  	Option<bool>   varnorm;
 	Option<bool>   varnorm2;
  	Option<bool>   pbsc;
  	Option<bool>   pspec;
  	Option<string> segment;
  	Option<bool>   tsmooth;
  	Option<float>  epsilon;
  	Option<float>  epsilonS;
  	Option<int>    maxNumItt;
  	Option<int>    maxRestart;
  	Option<int>    rank1interval;

  	Option<string> mmthresh;
  	Option<bool>   perf_mm;
  	Option<string> ICsfname;
  	Option<string> filtermix; 
  	Option<string> smodename; 
  	Option<string> filter; 

  	Option<bool>   genreport;
	Option<string> guireport;
	Option<string> bgimage;
  	Option<float>  tr;
  	Option<bool>   logPower;
	Option<bool>   addsigchng;
	Option<bool>   allPPCA;
	Option<bool>   varplots;
	Option<bool>   varvals;

	Option<string> fn_Tdesign;
	Option<string> fn_Tcon;
	Option<string> fn_TconF;
	Option<string> fn_Sdesign;
	Option<string> fn_Scon;	
	Option<string> fn_SconF;	
	
  	Option<bool>   output_all;
  	Option<bool>   output_unmix;
  	Option<bool>   output_MMstats;
  	Option<bool>   output_pca;
  	Option<bool>   output_white;
  	Option<bool>   output_origIC;
  	Option<bool>   output_mean;

  	Option<bool> verbose;  
  	Option<bool> vers;
  	Option<bool> copyright;
  	Option<bool> help;
  	Option<bool> debug;

  	Option<string> guessfname;
  	Option<string> paradigmfname;
  	Option<string> axials_str;

  	Option<int>   dummy;
  	Option<int>   repeats;
	Option<int>   seed;
  	Option<float> nlconst1;
  	Option<float> nlconst2;
  	Option<float> smooth_probmap;
	Option<string> insta_fn;

  	Option<bool> remove_meanvol;
  	Option<bool> remove_meantc;
 	Option<bool> remove_endslices;
  	Option<bool> rescale_nht;

  	Option<bool> guess_remderiv;
  	Option<bool> temporal;

  	Option<float> retryfactor;
  	Option<int> econ;

  	void parse_command_line(int argc, char** argv, Log& logger,  const string &p_version);

 	private:
  	MelodicOptions();  
  	const MelodicOptions& operator=(MelodicOptions&);
  	MelodicOptions(MelodicOptions&);

  	OptionParser options; 
      
  	static MelodicOptions* gopt;

  	void print_usage(int argc, char *argv[]);  
  	void print_version(); 
  	void print_copyright();
  	void status();

};

 inline MelodicOptions& MelodicOptions::getInstance(){
   if(gopt == NULL)
     gopt = new MelodicOptions();
   
   return *gopt;
 }

 inline MelodicOptions::MelodicOptions() :
   logdir(string("-o,--outdir"), string("log.txt"),
	   string("output directory name\n"), 
	   false, requires_argument),
   inputfname(string("-i,--in"), std::vector<string>(),
	   string("input file names (either single file name or comma-separated list or text file)"), 
	   true, requires_argument),
   outputfname(string("-O,--out"), string("melodic"),
	   string("output file name"), 
	   false, requires_argument,false),
   maskfname(string("-m,--mask"), string(""),
	   string("file name of mask for thresholding"), 
	   false, requires_argument),
   use_mask(string("--nomask"), true,
	   string("switch off masking"), 
	   false, no_argument),
   update_mask(string("--update_mask"), true,
	   string("switch off mask updating"), 
	   false, no_argument),
   perf_bet(string("--nobet"), true,
	   string("\tswitch off BET"), 
	   false, no_argument),
   threshold(string("--bgthreshold"),  0.01,
	   string("brain / non-brain threshold (only if --nobet selected)\n"), 
	   false, requires_argument),
   pca_dim(string("-d,--dim"), 0,
	   string("dimensionality reduction into #num dimensions (default: automatic estimation)"), 
	   false, requires_argument),
   pca_est(string("--dimest"), string("lap"),
	   string("use specific dim. estimation technique: lap, bic, mdl, aic, mean (default: lap)"), 
	   false, requires_argument),
   joined_whiten(string("--sep_whiten"), false,
	   string("switch on separate whitening"), 
	   false, no_argument, false),
   joined_vn(string("--sep_vn"), true,
   	   string("switch on separate variance nomalisation (as opposed to separate VN)"), 
       false, no_argument),
   dr_pca(string("--mod_pca"), true,
	   string("switch off modified PCA for concat ICA"),
	   false, no_argument, false),
   migp(string("--migp"), false,
	   string("switch on MIGP data reduction"),
	   false, no_argument),	
   migpN(string("--migpN"), 0,
	   string("Number of internal Eigenmaps"),
	   false, requires_argument),
   migp_shuffle(string("--migp_shuffle"), true,
	   string("Randomise MIGP file order (default: TRUE)"),
	   false, no_argument),
   migp_factor(string("--migp_factor"), 2,
	   string("Internal Factor of mem-threshold relative to number of Eigenmaps (default: 2)"),
       false, requires_argument),
   dr(string("--dr"), false,
	   string("Dual Regression (default: false)"),
	   false, no_argument, false),
   dr_out(string("--dr_out"), false,
	   string("Dual Regression output for MIGP/concat ICA"),
	   false, no_argument, false),	
   vn_level(string("--vn_level"), float(2.3),
	   string("variance nomalisation threshold level (Z> value is ignored)"), 
	   false, requires_argument, false),
   numICs(string("-n,--numICs"), -1,
	   string("numer of IC's to extract (for deflation approach)"), 
	   false, requires_argument),
   approach(string("-a,--approach"),  string("symm"),
	   string("approach for decomposition, 2D: defl, symm (default), 3D: tica, concat (default)"),
	   false, requires_argument),
   nonlinearity(string("--nl"), string("pow3"),
	   string("\tnonlinearity: gauss, tanh, pow3, pow4"), 
	   false, requires_argument),
   varnorm(string("--vn,--varnorm"), true,
	   string("switch off variance normalisation"), 
	   false, no_argument),
   varnorm2(string("--vn2"), true,
		string("switch off 2nd level variance normalisation"), 
		false, no_argument, false),
   pbsc(string("--pbsc"), false,
	   string("        switch on conversion to percent BOLD signal change"), 
	   false, no_argument, false),
   pspec(string("--pspec"), false,
	   string("        switch on conversion to powerspectra"), 
	   false, no_argument, false),
   segment(string("--covarweight"), string(""),
	   string("voxel-wise weights for the covariance matrix (e.g. segmentation information)"),
	   false, requires_argument),
   tsmooth(string("--spca"),  false,
	   string("smooth the eigenvectors prior to IC decomposition"), 
	    false, no_argument, false),
   epsilon(string("--eps"), 0.00005,
	   string("minimum error change"), 
	   false, requires_argument),
   epsilonS(string("--epsS"), 0.03,
	   string("minimum error change for rank-1 approximation in TICA"), 
	   false, requires_argument),
   maxNumItt(string("--maxit"),  500,
	   string("\tmaximum number of iterations before restart"), 
	   false, requires_argument),
   maxRestart(string("--maxrestart"),  -1,
	   string("maximum number of restarts\n"), 
	   false, requires_argument),
   rank1interval(string("--rank1interval"),  10,
	   string("number of iterations between rank-1 approximation (TICA)\n"), 
		 false, requires_argument,false),
   mmthresh(string("--mmthresh"), string("0.5"),
	   string("threshold for Mixture Model based inference"), 
	   false, requires_argument),
   perf_mm(string("--no_mm"), true,
	   string("\tswitch off mixture modelling on IC maps\n "), 
	   false, no_argument),
   ICsfname(string("--ICs"), string(""),
	   string("\tinput filename of the IC components file for mixture modelling"), 
	   false, requires_argument),
   filtermix(string("--mix"),  string(""),
	   string("\tinput filename of mixing matrix for mixture modelling / filtering"), 
	   false, requires_argument),
   smodename(string("--smode"),  string(""),
	   string("\tinput filename of matrix of session modes for report generation"), 
	   false, requires_argument),
   filter(string("-f,--filter"),  string(""),
	   string("list of component numbers to remove\n "), 
	   false, requires_argument),
   genreport(string("--report"), false,
	   string("generate Melodic web report"), 
	   false, no_argument),
   guireport(string("--guireport"), string(""),
	   string("modify report for GUI use"), 
	   false, requires_argument, false),
   bgimage(string("--bgimage"),  string(""),
	   string("specify background image for report (default: mean image)\n "), 
	   false, requires_argument),
   tr(string("--tr"),  0.0,
	   string("\tTR in seconds"), 
	   false, requires_argument),
   logPower(string("--logPower"),  false,
	   string("calculate log of power for frequency spectrum\n"), 
	   false, no_argument),
   addsigchng(string("--sigchng"),  false,
	   string("add signal change estimates to report pages\n"), 
       false, no_argument, false),
   allPPCA(string("--allPPCA"),  false,
	   string("add all PPCA plots\n"), 
	   false, no_argument, false),
   varplots(string("--varplots"),  false,
	   string("add std error envelopes to time course plots\n"), 
	   false, no_argument, false),
   varvals(string("--varvals"),  false,
	   string("add rank1 values after plots\n"), 
	   false, no_argument, false),
   fn_Tdesign(string("--Tdes"), string(""),
	   string("        design matrix across time-domain"),
	   false, requires_argument),
   fn_Tcon(string("--Tcon"), string(""),
       string("        t-contrast matrix across time-domain"),
	   false, requires_argument),
   fn_TconF(string("--Tconf"), string(""),
	   string("        F-contrast matrix across time-domain"),
	   false, requires_argument, false),
   fn_Sdesign(string("--Sdes"), string(""),
	   string("        design matrix across subject-domain"),
	   false, requires_argument),
   fn_Scon(string("--Scon"), string(""),
	   string("        t-contrast matrix across subject-domain"),
	   false, requires_argument),	
   fn_SconF(string("--Sconf"), string(""),
	   string("        F-contrast matrix across subject-domain"),
	   false, requires_argument,false),	
   output_all(string("--Oall"),  false,
	   string("        output everything"), 
	   false, no_argument),
   output_unmix(string("--Ounmix"),  false,
	   string("output unmixing matrix"), 
	   false, no_argument),
   output_MMstats(string("--Ostats"),  false,
	   string("output thresholded maps and probability maps"), 
	   false, no_argument),
   output_pca(string("--Opca"),  false,
	   string("\toutput PCA results"), 
	   false, no_argument),
   output_white(string("--Owhite"),  false,
	   string("output whitening/dewhitening matrices"), 
	   false, no_argument),
   output_origIC(string("--Oorig"),  false,
	   string("\toutput the original ICs"), 
	   false, no_argument),
   output_mean(string("--Omean"),  false,
	   string("\toutput mean volume\n"), 
	   false, no_argument),
   verbose(string("-v,--verbose"), false,
	   string("switch on diagnostic messages"), 
	   false, no_argument),
   vers(string("-V,--version"), false,
	   string("prints version information"), 
	   false, no_argument),
   copyright(string("--copyright"), false,
	   string("prints copyright information"), 
	   false, no_argument),
   help(string("-h,--help"),  false,
	   string("prints this help message"), 
	   false, no_argument),
   debug(string("--debug"),  false,
	   string("        switch on debug messages"), 
	   false, no_argument),
   guessfname(string("--init_ica"), string(""),
	   string("file name of FEAT paradigm file (design.mat) for ICA initialisation"), 
	   false, requires_argument, false),
   paradigmfname(string("--init_pca"),  string(""),
	   string("file name of FEAT paradigm file (design.mat) for PCA initialisation"), 
	   false, requires_argument, false),
   axials_str(string("--report_maps"),  string(" -s 2 -A 950 "),
		   string("control string for spatial map images (see slicer)"), 
		   false, requires_argument),
   dummy(string("--dummy"),  0,
	   string("number of dummy volumes"), 
	   false, requires_argument,false),
   repeats(string("--repeats"), 1,
	   string("number of repeats (multistart)"), 
	   false, requires_argument, false),
   seed(string("--seed"), -1,
	   string("integer seed for random number generator within melodic"), 
	   false, requires_argument, false),
   nlconst1(string("--nl1,--nlconst1"),  1.0,
	   string("nonlinear constant 1"), 
	   false, requires_argument, false),
   nlconst2(string("--nl2,--nlconst2"),  1.0,
	   string("nonlinear constant 2"), 
	   false, requires_argument, false),
   smooth_probmap(string("--smooth_pm"),  0.0,
	   string("width of smoothing kernel for probability maps"), 
	   false, requires_argument, false),
   insta_fn(string("--insta_fn"), string(""),
	   string(" mask file name for instacorr calculation"),
       false, requires_argument, false),
   remove_meanvol(string("--keep_meanvol"), true,
	   string("do not subtract mean volume"), 
	   false, no_argument, false),
   remove_meantc(string("--remove_meantc"), false,
	   string("remove mean time course"), 
	   false, no_argument, false),
   remove_endslices(string("--remEndslices"),  false,
	   string("delete end slices (motion correction artefacts)"), 
	   false, no_argument,false),
   rescale_nht(string("--rescale_nht"),  true,
	   string("switch off map rescaling after mixture-modelling"), 
	   false, no_argument,false),
   guess_remderiv(string("--remove_deriv"),  false,
	   string("removes every second entry in paradigm file (EV derivatives)"), 
	   false, no_argument, false),
   temporal(string("--temporal"),  false,
	   string("perform temporal ICA"), 
	   false, no_argument, false),
   retryfactor(string("--retryfactor"), float(0.95),
		string("multiplicative factor for determining new dim if estimated dim fails to converge"),
		false, requires_argument, false),
   econ(string("--econ"), 20000, 
	   string("set ctrl parameter for helperfns econ mode"),
       false, requires_argument, false),
   options(title, usageexmpl)
   {
     try {  
      options.add(logdir);
	    options.add(inputfname);
	    options.add(outputfname);
	    options.add(guessfname);
	    options.add(maskfname);
	    options.add(use_mask);
	    options.add(update_mask);
	    options.add(perf_bet);
	    options.add(threshold);
	    options.add(pca_dim);
	    options.add(pca_est);
	    options.add(joined_whiten);
	    options.add(joined_vn);
		options.add(dr_pca);
		options.add(migp);
		options.add(migpN);
		options.add(migp_shuffle);
		options.add(migp_factor);
		options.add(dr);
		options.add(dr_out);
	    options.add(vn_level);
	    options.add(numICs);
	    options.add(approach);
	    options.add(nonlinearity);
	    options.add(varnorm);
		options.add(varnorm2);
	    options.add(pbsc);
	    options.add(pspec);
	    options.add(segment);
	    options.add(tsmooth);
	    options.add(epsilon);
	    options.add(epsilonS);
	    options.add(maxNumItt);
	    options.add(maxRestart);
	    options.add(rank1interval);
	    options.add(mmthresh);
	    options.add(perf_mm);
	    options.add(ICsfname);
	    options.add(filtermix);
	    options.add(smodename);
	    options.add(filter);
	    options.add(genreport);
	    options.add(guireport);
		options.add(bgimage);
	    options.add(tr);
	    options.add(logPower);
	    options.add(addsigchng);
	    options.add(allPPCA);
	    options.add(varplots);
	    options.add(varvals);
		options.add(fn_Tdesign);
		options.add(fn_Tcon);
		options.add(fn_TconF);
		options.add(fn_Sdesign);
		options.add(fn_Scon);
		options.add(fn_SconF);
	    options.add(output_all);
	    options.add(output_unmix);
	    options.add(output_MMstats);
	    options.add(output_pca);
	    options.add(output_white);
	    options.add(output_origIC);
	    options.add(output_mean);
	    options.add(verbose);
	    options.add(vers);
	    options.add(copyright);
	    options.add(help);
	    options.add(debug);
	   
	    options.add(guessfname);
	    options.add(paradigmfname); 
	    options.add(axials_str); 
	    options.add(dummy);
	    options.add(repeats);
	    options.add(seed);
	    options.add(nlconst1);
	    options.add(nlconst2);
	    options.add(smooth_probmap);
		options.add(insta_fn);
	    options.add(remove_meanvol);
	    options.add(remove_meantc);
	    options.add(remove_endslices);
	    options.add(rescale_nht);
	    options.add(guess_remderiv);
	    options.add(temporal);
		options.add(retryfactor);
		options.add(econ);
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



