/*  pvmfitOptions.h

    Saad Jbabdi, FMRIB Image Analysis Group

    Copyright (C) 1999-2009 University of Oxford  */

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

#if !defined(pvmfitOptions_h)
#define pvmfitOptions_h

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "commonopts.h"
//#include "newmatall.h"

using namespace Utilities;

namespace PVMFIT {

class pvmfitOptions {
 public:
  static pvmfitOptions& getInstance();
  ~pvmfitOptions() { delete gopt; }
  
  Option<bool>   verbose;
  Option<bool>   help;
  Option<string> datafile;
  Option<string> ofile;
  Option<string> maskfile;
  Option<string> bvecsfile;
  Option<string> bvalsfile;
  Option<int>    nfibres;
  Option<int>    modelnum;
  Option<bool>   all;
  Option<bool>   cnonlinear;
  Option<bool>   cnonlinear_Fanning;
  Option<bool>   gridsearch;
  Option<bool>   use_f0;
  Option<bool>   saveBIC;

  bool parse_command_line(int argc, char** argv);
  
 private:
  pvmfitOptions();  
  const pvmfitOptions& operator=(pvmfitOptions&);
  pvmfitOptions(pvmfitOptions&);

  OptionParser options; 
      
  static pvmfitOptions* gopt;
  
};

 inline pvmfitOptions& pvmfitOptions::getInstance(){
   if(gopt == NULL)
     gopt = new pvmfitOptions();
   
   return *gopt;
 }

 inline pvmfitOptions::pvmfitOptions() :
  verbose(string("-V,--verbose"), false, 
	  string("switch on diagnostic messages"), 
	  false, no_argument),
   help(string("-h,--help"), false,
	string("display this message"),
	false, no_argument),
   datafile(string("-k,--data"), "",
	       string("data file"),
	       true, requires_argument),  
   ofile(string("-o,--out"), string("pvm"),
	       string("Output basename - default='pvm'"),
	       false, requires_argument),
   maskfile(string("-m,--mask"), "",
	    string("Bet binary mask file"),
	    true, requires_argument),
   bvecsfile(string("-r,--bvecs"), "",
	     string("b vectors file"),
	     true, requires_argument),  
   bvalsfile(string("-b,--bvals"), "",
	     string("b values file"),
	     true, requires_argument), 
   nfibres(string("-n,--nfibres"), 1,
	     string("number of fibres to fit - default=1"),
	     false, requires_argument), 
   modelnum(string("--model"), 1,
	     string("\t1:Ball-Sticks (single-shell); 2:Ball-Sticks (multi-shells); 4:Ball-Binghams"),
	     false, requires_argument), 
   all(string("--all"),false, 
	      string("\tFanning models: Fit all models from 1 up to N fibres and choose the best using BIC"),
	      false,no_argument),
   cnonlinear(string("--cnonlinear"),false, 
	      string("Model1: Apply constrained nonlinear optimization on the diffusivity, volume fractions and their sum"),
	      false,no_argument),
   cnonlinear_Fanning(string("--cnonlinear_F"),false, 
	      string("Model1: Apply constrained nonlinear optimization on the diffusivity, volume fractions and their sum.\n\t\t\tReturn n fanning angle estimates, using the Hessian of the cost function"),
	      false,no_argument),
   gridsearch(string("--gridsearch"),false, 
	      string("Use grid search (on the fanning eigenvalues). Default=off"),
	      false,no_argument),
   use_f0(string("--f0"),false, 
	      string("\tInclude noise floor parameter in the model"),
	      false,no_argument),
   saveBIC(string("--BIC"),false, 
	      string("\tSave BIC for certain models"),
	      false,no_argument),
   options("pvmfit", "pvmfit -k <datafile> -m <maskfile> -r <bvecsfile> -b <bvalsfile> [-n 2]\n")
   {
     
    
     try {
       options.add(verbose);
       options.add(help);
       options.add(datafile);
       options.add(ofile);
       options.add(maskfile);
       options.add(bvecsfile);
       options.add(bvalsfile);
       options.add(nfibres);
       options.add(modelnum);
       options.add(all);
       options.add(cnonlinear);
       options.add(cnonlinear_Fanning);
       options.add(gridsearch);
       options.add(use_f0);
       options.add(saveBIC);
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





