/*  bintoptions.h

    Mark Woolrich - FMRIB Image Analysis Group

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

#if !defined(BintOptions_h)
#define BintOptions_h

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "utils/log.h"

using namespace Utilities;

namespace Bint {

class BintOptions {

 public:

  virtual ~BintOptions(){};

  Option<bool> verbose;
  Option<int> debuglevel;
  Option<bool> timingon;
  Option<bool> help;
  Option<string> datafile;
  Option<string> maskfile;
  Option<string> logdir;
  Option<bool> forcedir;
  Option<string> inference;
  Option<int> njumps;
  Option<int> burnin;
  Option<int> sampleevery;
  Option<int> updateproposalevery;
  Option<float> acceptancerate;
  Option<int> seed;
  Option<float> prec;
  Option<bool> analmargprec;

  void parse_command_line(int argc, char** argv, Log& logger);
  
 protected:

  BintOptions(const string& str1, const string& str2);  
  OptionParser options; 

 private:

  BintOptions();
  const BintOptions& operator=(BintOptions&);
  BintOptions(BintOptions&);
  
};

 inline BintOptions::BintOptions(const string& str1, const string& str2) :
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
   datafile(string("--data,--datafile"), string("data"),
			  string("data regressor data file"),
		     true, requires_argument),
   maskfile(string("--mask,--maskfile"), string(""),
			  string("mask file"),
		     true, requires_argument),
   logdir(string("--ld,--logdir"), string("logdir"),
			  string("log directory (default is logdir)"),
		     false, requires_argument), 
   forcedir(string("--forcedir"), false,
		    string("Use the actual directory name given - i.e. don't add + to make a new directory"),
		    false, no_argument),
   inference(string("--inf,--inference"), string("mcmc"),
			  string("inference technique: mcmc\n laplace\n (default is mcmc)"),
		     false, requires_argument),  
   njumps(string("--nj,--njumps"), 5000,
			  string("Num of jumps to be made by MCMC (default is 5000)"),
		     false, requires_argument),
   burnin(string("--bi,--burnin"), 500,
			  string("Num of jumps at start of MCMC to be discarded (default is 500)"),
		     false, requires_argument),
   sampleevery(string("--se,--sampleevery"), 1,
			  string("Num of jumps for each sample (MCMC) (default is 1)"),
		     false, requires_argument),
   updateproposalevery(string("--upe,--updateproposalevery"), 40,
		       string("Num of jumps for each update to the proposal density std (MCMC) (default is 40)"),
		     false, requires_argument),
   acceptancerate(string("--arate,--acceptancerate"), 0.6,
			  string("Acceptance rate to aim for (MCMC) (default is 0.6)"),
		     false, requires_argument),
   seed(string("--seed"), 10, 
      string("seed for pseudo random number generator"), 
      false, requires_argument),
   prec(string("--prec"), -1, 
      string("value to fix error precision to (default is -1, which means error precision is not fixed)"), 
	false, requires_argument),
   analmargprec(string("--noamp"), true, 
      string("turn off Analytical Marginalisation of error Precision"), 
	false, no_argument),
   options(str1,str2)
   {
     try {
       options.add(verbose);
       options.add(debuglevel);
       options.add(timingon);
       options.add(help);
       options.add(datafile);
       options.add(maskfile);
       options.add(logdir);
       options.add(forcedir);
       options.add(inference);
       options.add(njumps);
       options.add(burnin);
       options.add(sampleevery);
       options.add(updateproposalevery);
       options.add(acceptancerate);
       options.add(seed);
       options.add(prec);
       options.add(analmargprec);
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



