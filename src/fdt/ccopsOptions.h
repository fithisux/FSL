/*  BpmOptions.h

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

#if !defined(ccopsOptions_h)
#define ccopsOptions_h

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "commonopts.h"
using namespace Utilities;

namespace CCOPS {

class ccopsOptions {
 public:
  static ccopsOptions& getInstance();
  ~ccopsOptions() { delete gopt; }
  
  Option<bool> help;
  Option<string> inmatrix;
  Option<string> basename;
  Option<string> directory;
  Option<string> excl_mask;
  Option<bool>  reord1;
  Option<bool>  reord2;
  Option<bool>  reord3;
  Option<float> connexity;
  Option<int>   bin;
  Option<float> power;
  Option<string> mask;
  Option<string> scheme;
  Option<int>    nclusters;
  bool parse_command_line(int argc, char** argv);
  
 private:
  ccopsOptions();  
  const ccopsOptions& operator=(ccopsOptions&);
  ccopsOptions(ccopsOptions&);

  OptionParser options; 
      
  static ccopsOptions* gopt;
  
};

 inline ccopsOptions& ccopsOptions::getInstance(){
   if(gopt == NULL)
     gopt = new ccopsOptions();
   
   return *gopt;
 }

 inline ccopsOptions::ccopsOptions() :
   help(string("-h,--help"), false,
	string("display this message"),
	false, no_argument),
   inmatrix(string("-i,--in"), string("fdt_matrix2"),
	       string("input matrix"),
	       false, requires_argument),  
   basename(string("-b,--basename"), string(""),
	       string("Output basename"),
	       true, requires_argument),
   directory(string("-d,--dir"), string("."),
	       string("Tractography Results Directory"),
	       false, requires_argument),
   excl_mask(string("-x"), string(""),
	     string("exclusion mask (in tract space)"),
	     false, requires_argument),  
   reord1(string("--r1"), bool(false),
	     string("do seedspace reordering (default no)"),
	     false, no_argument), 
   reord2(string("--r2"), bool(false),
	     string("do tractspace reordering (default no)"),
	     false, no_argument), 
   reord3(string("--tractreord"), bool(false),
	     string("propagate seed reordering onto tract space"),
	     false, no_argument), 
   connexity(string("--con"), 0.0,
	     string("add connexity constraint - value between 0 and 1 (0 is no constraint). default=0"),
	     false, requires_argument), 
   bin(string("--bin"), 0, 
	 string("binarise at (default 0 - no binarisation)"), 
	 false, requires_argument),
   power(string("-p,--power"), 1, 
	 string("power to raise the correlation matrix to (default 1)"), 
	 false, requires_argument),
   mask(string("-m,--mask"), "", 
	 string("brain mask used to output the clustered roi mask (not necessary if --dir set)"), 
	 false, requires_argument),
   scheme(string("-s,--scheme"), "spectral", 
	 string("Reordering algorithm. Can be either spectral (default) or kmeans or fuzzy"), 
	 false, requires_argument),
   nclusters(string("-k,--nclusters"), 2, 
	  string("Number of clusters to be used in kmeans or fuzzy"), 
	  false, requires_argument),
   options("ccops","")
   {
     
    
     try {
       options.add(help);
       options.add(inmatrix);
       options.add(basename);
       options.add(directory);
       options.add(excl_mask);
       options.add(reord1);
       options.add(reord2);
       options.add(reord3);
       options.add(connexity);
       options.add(bin);
       options.add(power);
       options.add(mask);
       options.add(scheme);
       options.add(nclusters);
       
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





