/*  mmoptions.h

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

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

#if !defined(mmoptions_h)
#define mmoptions_h

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "utils/log.h"

using namespace Utilities;

namespace Mm {

class MmOptions {
 public:
  static MmOptions& getInstance();
  ~MmOptions() { delete gopt; }
  
  Option<bool> verbose;
  Option<int> debuglevel;
  Option<bool> timingon;
  Option<bool> help;
  Option<string> spatialdatafile;
  Option<string> epiexampledatafile;
  Option<string> maskfile;
  Option<string> logdir;
  Option<bool> nonspatial;
  Option<bool> fixmrfprec;
  Option<float> mrfprecstart;
  Option<float> mrfprecmultiplier;
  Option<float> initmultiplier;
  Option<bool> updatetheta;
  Option<bool> zfstatmode;
  Option<float> phi;
  Option<int> niters;
  Option<float> threshold;

  void parse_command_line(int argc, char** argv, Log& logger);
  
 private:
  MmOptions();  
  const MmOptions& operator=(MmOptions&);
  MmOptions(MmOptions&);

  OptionParser options; 
      
  static MmOptions* gopt;
  
};

 inline MmOptions& MmOptions::getInstance(){
   if(gopt == NULL)
     gopt = new MmOptions();
   
   return *gopt;
 }

 inline MmOptions::MmOptions() :
   verbose(string("-V,--verbose"), false, 
	   string("switch on diagnostic messages"), 
	   false, no_argument),
   debuglevel(string("--debug,--debuglevel"), 0,
		       string("set debug level"), 
		       false, requires_argument),
   timingon(string("--timingon"), false, 
		       string("turn timing on"), 
		       false, no_argument),
   help(string("-h,--help"), false,
		    string("display this message"),
		    false, no_argument),
   spatialdatafile(string("--sdf,--spatialdatafile"), string(""),
			  string("spatial map data file"),
		     true, requires_argument),
   epiexampledatafile(string("--edf,--epiexampledatafile"), string(""),
			  string("example epi data file"),
		     false, requires_argument),  
   maskfile(string("-m,--mask"), string(""),
			  string("mask file"),
		     true, requires_argument),
   logdir(string("-l,--ld,--logdir"), string("logdir"),
			  string("log directory"),
		     false, requires_argument),    
   nonspatial(string("--ns,--nonspatial"), false, 
	   string("Nonspatial mixture model"),		       
	   false, no_argument),
   fixmrfprec(string("--fmp,--fixmrfprec"), false, 
	   string("Fix MRF precision to mrfprecstart throughout"),		       
	   false, no_argument),
   mrfprecstart(string("--mps,--mrfprecstart"), 10, 
	   string("MRF precision initial value (default is 10)"),		       
	   false, requires_argument),
   mrfprecmultiplier(string("--mpm,--mrfprecmultiplier"), -1, 
	   string("Update multiplier for MRF precision (default is -1, do not multipy)"),  
	   false, requires_argument),
   initmultiplier(string("--im,--initmultiplier"), 0.3, 
	   string("Init multiplier (default is 0.3)"),		       
	   false, requires_argument),
   updatetheta(string("--nut,--noupdatetheta"), true, 
	       string("Turn off updating of distribution parameters after non-spatial fit"),false, no_argument),  
   zfstatmode(string("--zfstatmode"), false, 
	      string("Turn on zfstat mode - this enforces no deactivation class"),false, no_argument),  
   phi(string("--phi"), 0.05,
       string("phi (default is 0.05)"), 
       false, requires_argument), 
   niters(string("--ni,--niters"), -1, 
	string("niters (default is -1: auto stop)"), 
	false, requires_argument), 
   threshold(string("--th,--thresh"), 0.5, 
	string("threshold for use when displaying classification maps in MM.html report (default is 0.5, -1 indicates no thresholding)"), 
	false, requires_argument),
   options("mm", "mm -d <filename>\nmm --verbose\nmm --mask=<mask> --sdf=<filename> --logdir=<logdir>")
   {
     try {
       options.add(verbose);
       options.add(debuglevel);
       options.add(timingon);
       options.add(help);
       options.add(spatialdatafile);
       options.add(epiexampledatafile);
       options.add(maskfile);
       options.add(logdir);
       options.add(nonspatial);
       options.add(fixmrfprec);
       options.add(mrfprecstart);
       options.add(mrfprecmultiplier);
       options.add(initmultiplier);
       options.add(updatetheta);
       options.add(zfstatmode);
       options.add(phi);
       options.add(niters);   
       options.add(threshold);       
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



