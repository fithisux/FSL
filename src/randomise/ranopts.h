/*  ranopts.h

    Matthew Webster, Tim Behrens & Steve Smith (FMRIB) & Tom Nichols (UMich)

    Copyright (C) 2008-2010 University of Oxford  */

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

#if !defined(ranopts_h)
#define ranopts_h

#include <string>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "utils/log.h"

using namespace Utilities;

namespace RANDOMISE {

class ranopts {
 public:
  static ranopts& getInstance();
  ~ranopts() { delete gopt; }
  
  Option<bool> demean_data;
  Option<bool> one_samp;
  Option<string> in_fileroot;
  Option<string> maskname;
  Option<string> out_fileroot;
  Option<string> dm_file;
  Option<string> tc_file;
  Option<string> fc_file;
  Option<string> gp_file;
  Option<string> effectiveDesignFile;
  Option<bool> how_many_perms;
  Option<bool> parallelData;
  Option<int> n_perm;
  Option<bool> voxelwiseOutput;
  Option<bool> doFOnly;
  Option<bool> tfce;
  Option<bool> tfce2D;
  Option<float> cluster_thresh;
  Option<float> clustermass_thresh;
  Option<float> f_thresh;
  Option<float> fmass_thresh;
  Option<float> tfce_height;
  Option<float> tfce_delta;
  Option<float> tfce_size;
  Option<int> tfce_connectivity;
  Option<float> var_sm_sig;
  Option<bool> help;
  Option<bool> verbose;
  Option<bool> cluster_norm;
  Option<bool> outputRaw;
  Option<bool> outputUncorr;
  Option<bool> outputTextPerm;
  Option<bool> outputTextNull;
  Option<bool> output_permstat;
  Option<bool> disableNonConstantMask;
  Option<int> randomSeed;
  Option<vector<int> > voxelwise_ev_numbers;
  Option<vector<string> > voxelwise_ev_filenames;
  Option<int> nMultiVariate;
  Option<bool> isDebugging;
  Option<int> confoundMethod;
  Option<bool> detectNullSubjects;
  Option<bool> permuteBlocks;
  Option<bool> verbose_old;
  Option<int>  skipTo;
  Option<bool> outputGlm;

  void parse_command_line(int argc, char** argv,Log& logger);
  
 private:
  ranopts();  
  const ranopts& operator=(ranopts&);
  ranopts(ranopts&);

  OptionParser options; 
      
  static ranopts* gopt;
  
};

 inline ranopts& ranopts::getInstance(){
   if(gopt == NULL)
     gopt = new ranopts();
   
   return *gopt;
 }

 inline ranopts::ranopts() :
   demean_data(string("-D"), false,
	string("\tdemean data temporally before model fitting ( demean model as well if required )"),
	false, no_argument),
   one_samp(string("-1"), false,
	    string("\tperform 1-sample group-mean test instead of generic permutation test"),
	    false, no_argument),
   in_fileroot(string("-i"), "",
       string("~<input>\t4D input image"),
       true, requires_argument),
   maskname(string("-m"), "",
       string("~<mask>\tmask image"),
       false, requires_argument),
   out_fileroot(string("-o"), string(""),
	    string("~<out_root>\toutput file-rootname"),
	    true, requires_argument),  
   dm_file(string("-d"), string(""),
	    string("~<design.mat>\tdesign matrix file"),
	    false, requires_argument),  
   tc_file(string("-t"), string(""),
	    string("~<design.con>\tt contrasts file"),
	    false, requires_argument),  
   fc_file(string("-f"), string(""),
            string("~<design.fts>\tf contrasts file"),
	   false, requires_argument),  
   gp_file(string("-e"), string(""),
            string("~<design.grp>\texchangeability block labels file"),
	   false, requires_argument),  
   effectiveDesignFile(string("--effective_design"), string(""),
            string("~<design2.mat>\talternative design for determining valid permutations"),
	   false, requires_argument),
   how_many_perms(string("-q"), false,
	    string("\tprint out how many unique permutations would be generated and exit"),
	    false, no_argument), 
   parallelData(string("-Q"), false,
	    string("\tprint out information required for parallel mode and exit"),
	    false, no_argument),
   n_perm(string("-n"), 5000,
	    string("~<n_perm>\tnumber of permutations (default 5000, set to 0 for exhaustive)"),
	    false, requires_argument),  
   voxelwiseOutput(string("-x"),false,
	    string("\toutput voxelwise corrected p-value images"),
		 false, no_argument),  
   doFOnly(string("--fonly"), false, 
	   string("\tcalculate f-statistics only"), 
	   false, no_argument),
   tfce(string("-T"), false, 
	   string("\tcarry out Threshold-Free Cluster Enhancement"), 
	   false, no_argument),
   tfce2D(string("--T2"), false, 
	   string("\tcarry out Threshold-Free Cluster Enhancement with 2D optimisation (e.g. for TBSS data); H=2, E=1, C=26"), 
	   false, no_argument),
   cluster_thresh(string("-c"), -1,
	  string("~<thresh>\tcarry out cluster-based thresholding"),
	  false, requires_argument),
   clustermass_thresh(string("-C"), -1,
	  string("~<thresh>\tcarry out cluster-mass-based thresholding"),
	  false, requires_argument),
   f_thresh(string("-F"), -1,
	  string("~<thresh>\tcarry out f cluster thresholding"),
	  false, requires_argument),
   fmass_thresh(string("-S"), -1,
	  string("~<thresh>\tcarry out f cluster-mass thresholding"),
	  false, requires_argument),
   tfce_height(string("--tfce_H"), 2, string("~<H>\tTFCE height parameter (default=2)"), false, requires_argument),
   tfce_delta(string("--tfce_D"), 1, string("~<H>\tTFCE delta parameter overide"), false, requires_argument),
   tfce_size(string("--tfce_E"), 0.5, string("~<E>\tTFCE extent parameter (default=0.5)"), false, requires_argument),
   tfce_connectivity(string("--tfce_C"), 6, string("~<C>\tTFCE connectivity (6 or 26; default=6)"), false, requires_argument),
   var_sm_sig(string("-v"), 0,
	    string("~<std>\tuse variance smoothing for t-stats (std is in mm)"),
	     false, requires_argument),
   help(string("-h,--help"), false,
	string("display this message"),
	false, no_argument),
   verbose(string("--quiet"), true, 
	   string("\tswitch off diagnostic messages"),
	   false, no_argument),
   cluster_norm(string("--twopass"), false, 
	   string("carry out cluster normalisation thresholding"), 
		false, no_argument),
   outputRaw(string("-R"), false, 
	   string("\toutput raw ( unpermuted ) statistic images"), 
		false, no_argument),
   outputUncorr(string("--uncorrp"), false, 
	   string("output uncorrected p-value images"), 
		false, no_argument),
   outputTextPerm(string("-P"), false, 
	   string("\toutput permutation vector text file"), 
		false, no_argument),
   outputTextNull(string("-N"), false, 
	   string("\toutput null distribution text files"), 
		false, no_argument),
   output_permstat(string("--permout"), false, 
	   string("\toutput permuted tstat"), 
		false, no_argument,false),
   disableNonConstantMask(string("--norcmask"), false, 
	   string("don't remove constant voxels from mask"), 
		false, no_argument),
   randomSeed(string("--seed"),0,
	    string("~<seed>\tspecific integer seed for random number generator"),
		 false, requires_argument),
   voxelwise_ev_numbers(string("--vxl"), vector<int>(), 
         string("\tlist of numbers indicating voxelwise EVs position in the design matrix (list order corresponds to files in vxf option). caution BETA option."), 
         false, requires_argument),
   voxelwise_ev_filenames(string("--vxf"), vector<string>(), 
         string("\tlist of 4D images containing voxelwise EVs (list order corresponds to numbers in vxl option). caution BETA option."), 
         false, requires_argument),
   nMultiVariate(string("--multi"),1,
	    string("~<dim>\tmultivariate dimension (default 1). caution BETA option."),
		 false, requires_argument, false),
   isDebugging(string("--debug"), false, 
	   string("\tOutput debug information"),
	       false, no_argument,false),
   confoundMethod(string("-U"),1,
	    string("~<mode>\tconfound mode. 0: Kennedy Y_a on X_a (old) 1: Freedman-Lane Y_a on X|Z (default) 2: Y on X|Z. 3: ter Braak Y_aFull on X|Z Caution BETA option."),
		  false, requires_argument, false),
   detectNullSubjects(string("--detectNull"), false, 
	   string("attempt to detect uninformative rows in the effective regressor and not permute them"), 
		      false, no_argument, false),
   permuteBlocks(string("--permuteBlocks"), false, 
	   string("permute exchangeability blocks. Caution BETA option"), 
		 false, no_argument),
   verbose_old(string("-V"), false, 
	   string("\tswitch on diagnostic messages (deprecated: now always on unless quiet)"),
	       false, no_argument, false),
   skipTo(string("--skipTo"),0,
	    string("~<number>\tdo only contrast<number>. A setting of 0 will do all contrasts, numbering starts with f- and then t-contrasts"),
		false, requires_argument, false),
   outputGlm(string("--glm_output"), false, 
	   string("output glm information for t-statistics ( unpermuted case only )"), 
		false, no_argument),

   options("randomise v2.9", "randomise -i <input> -o <output> -d <design.mat> -t <design.con> [options]")
     {
    
     try {
       options.add(demean_data);
       options.add(one_samp);
       options.add(in_fileroot);
       options.add(maskname);
       options.add(out_fileroot);
       options.add(dm_file);
       options.add(tc_file);
       options.add(fc_file);
       options.add(gp_file);
       options.add(effectiveDesignFile);
       options.add(how_many_perms);
       options.add(parallelData);
       options.add(n_perm);    
       options.add(voxelwiseOutput);   
       options.add(doFOnly);
       options.add(tfce);
       options.add(tfce2D);
       options.add(cluster_thresh);
       options.add(clustermass_thresh);
       options.add(f_thresh);     
       options.add(fmass_thresh);
       options.add(var_sm_sig);
       options.add(help);
       options.add(verbose);
       options.add(cluster_norm);
       options.add(outputRaw);
       options.add(outputUncorr);
       options.add(outputTextPerm);
       options.add(outputTextNull);
       options.add(output_permstat);
       options.add(disableNonConstantMask);
       options.add(randomSeed);
       options.add(tfce_height);   
       options.add(tfce_delta);  
       options.add(tfce_size);     
       options.add(tfce_connectivity); 
       options.add(voxelwise_ev_numbers);
       options.add(voxelwise_ev_filenames); 
       options.add(nMultiVariate); 
       options.add(isDebugging);
       options.add(confoundMethod);
       options.add(detectNullSubjects);
       options.add(permuteBlocks);
       options.add(verbose_old);
       options.add(skipTo);
       options.add(outputGlm);
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

