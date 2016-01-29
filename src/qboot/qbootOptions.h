/*  qbootOptions.h

    Stamatios Sotiropoulos - FMRIB Image Analysis Group

    Copyright (C) 2010 University of Oxford  */

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

#if !defined(qbootOptions_h)
#define qbootOptions_h

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "utils/log.h"
#include "utils/tracer_plus.h"

using namespace Utilities;

namespace ODFs {

class qbootOptions {
 public:
  static qbootOptions& getInstance();
  ~qbootOptions() { delete gopt; }
  
  Option<string> logdir;
  Option<bool> forcedir;
  Option<string> datafile;
  Option<string> maskfile;
  Option<string> bvecsfile;
  Option<string> bvalsfile;
  Option<string> qshellsfile;
  Option<int> modelnum;
  Option<int> lmax;
  Option<int> npeaks;
  Option<float> peak_threshold;
  FmribOption<int> peak_finder;
  Option<int> nsamples;
  Option<float> lambda;
  Option<float> delta;
  Option<float> alpha;
  Option<int> seed;
  Option<bool> gfa;
  Option<bool> savecoeff;
  Option<bool> savemeancoeff;
  Option<bool> verbose;
  Option<bool> help;
  
  //Explanation of possible Output files
  //a) If --savecoeff, then only the samples of the coefficients will be saved and no peaks. This is very fast, but the output files can be huge.
  //b) If --savecoeff --savemeancoeff, then only the mean ODF coefficients (across the samples) will be saved for each voxel and again no peaks. This is very fast and also requires less memory/storage.
  //c) If --savemeancoeff, then the mean ODF coefficients (across the samples) will be saved for each voxel, along with the ODF peaks.
  //d) If none of the above is set, only the ODF peaks are saved.


  void parse_command_line(int argc, char** argv,  Log& logger);
  
 private:
  qbootOptions();  
  const qbootOptions& operator=(qbootOptions&);
  qbootOptions(qbootOptions&);

  OptionParser options; 
      
  static qbootOptions* gopt;
};

 inline qbootOptions& qbootOptions::getInstance(){
   if(gopt == NULL)
     gopt = new qbootOptions();
   
   return *gopt;
 }

 inline qbootOptions::qbootOptions() :
   logdir(string("--ld,--logdir"), string("logdir"),
	 string("Output directory (default is logdir)"),
	 false, requires_argument),
   forcedir(string("--forcedir"),false,string("Use the actual directory name given - i.e. don't add + to make a new directory"),false,no_argument),
   datafile(string("-k,--data"), string("data"),
	      string("Data file"),
	      true, requires_argument),  
   maskfile(string("-m,--mask"), string("nodif_brain_mask"),
	    string("Mask file"),
	    true, requires_argument),
   bvecsfile(string("-r,--bvecs"), string("bvecs"),
	     string("b vectors file"),
	     true, requires_argument),  
   bvalsfile(string("-b,--bvals"), string("bvals"),
	     string("b values file"),
	     true, requires_argument), 
   qshellsfile(string("--q"), string(""),
   	     string("\tFile provided with multi-shell data. Indicates the number of directions for each shell"),
   	     false, requires_argument), 
   modelnum(string("--model"),2, 
   	    string("\tWhich model to use. 1=Tuch's ODFs, 2=CSA ODFs (default), 3=multi-shell CSA ODFs"), 
   	    false,requires_argument),
   lmax(string("--lmax"),4, 
   	    string("\tMaximum spherical harmonic oder employed (must be even, default=4)"), 
   	    false,requires_argument),
   npeaks(string("--npeaks"),2, 
   	    string("Maximum number of ODF peaks to be detected (default 2)"), 
   	    false,requires_argument),
   peak_threshold(string("--thr"),0.4, 
   	    string("\tMinimum threshold for a local maxima to be considered an ODF peak. Expressed as a fraction of the maximum ODF value (default 0.4)"), 
   	    false,requires_argument),
   peak_finder(string("--pf"),1, 
   	    string("\tWhich peak finder to use. 1=Discrete, 2=Semi-continuous (can be only used with lmax=4) (default=1)"), 
   	    false,requires_argument),
   nsamples(string("--ns,--nsamples"),50, 
	   string("Number of bootstrap samples (default is 50)"), 
	   false,requires_argument),
   lambda(string("--lambda"),0,  //use 0.006 for model1
	   string("Laplace-Beltrami regularization parameter (default is 0)"), 
	   false,requires_argument),
   delta(string("--delta"),0.01, 
	   string("\tSignal attenuation regularization parameter for models=2,3 (default is 0.01)"), 
	   false,requires_argument),
   alpha(string("--alpha"),0, 
	   string("\tLaplacian sharpening parameter for model=1 (default is 0, should be smaller than 1)"), 
	   false,requires_argument),
   seed(string("--seed"),8665904,
	   string("\tSeed for pseudo-random number generator"),
           false,requires_argument),
   gfa(string("--gfa"),false,
	   string("\tCompute a generalised FA, using the mean ODF in each voxel"),
           false,no_argument),
   savecoeff(string("--savecoeff"),false,
	   string("Save the ODF coefficients instead of the peaks. WARNING: These can be huge files, please use a few bootstrap samples and a low lmax!"),
           false,no_argument),
   savemeancoeff(string("--savemeancoeff"),false,
	   string("Save the mean ODF coefficients across all samples"),
           false,no_argument),
   verbose(string("-V,--verbose"), false, 
	  string("Switch on diagnostic messages"), 
	  false, no_argument),
   help(string("-h,--help"), false,
	string("Display this message"),
	false, no_argument),
  
   options("qboot", "qboot --help (for list of options)\n")
   {
    
     try {
       options.add(logdir);
       options.add(forcedir);
       options.add(datafile);
       options.add(maskfile);
       options.add(bvecsfile);
       options.add(bvalsfile);
       options.add(qshellsfile);
       options.add(modelnum);
       options.add(lmax);
       options.add(npeaks);
       options.add(peak_threshold);
       options.add(peak_finder);
       options.add(nsamples);
       options.add(lambda);
       options.add(delta);
       options.add(alpha);
       options.add(seed);
       options.add(gfa);
       options.add(savecoeff);
       options.add(savemeancoeff);
       options.add(verbose);
       options.add(help);
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
