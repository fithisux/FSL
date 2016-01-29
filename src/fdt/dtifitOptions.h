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

#if !defined(dtifitOptions_h)
#define dtifitOptions_h

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "commonopts.h"
//#include "newmatall.h"

using namespace Utilities;

namespace DTIFIT {

class dtifitOptions {
 public:
  static dtifitOptions& getInstance();
  ~dtifitOptions() { delete gopt; }
  
  Option<bool> verbose;
  Option<bool> help;
  Option<string> dtidatafile;
  Option<string> ofile;
  Option<string> maskfile;
  Option<string> bvecsfile;
  Option<string> bvalsfile;
  Option<string> cni; //confounds of no interest. 
  Option<bool> sse; // Sum of squared errors
  Option<bool> wls; //Perform Weighted Least squares for tensor fitting 
  Option<bool> littlebit;
  Option<bool> savetensor;
  Option<int> z_min;
  Option<int> z_max;
  Option<int> y_min;
  Option<int> y_max;
  Option<int> x_min;
  Option<int> x_max;
  Option<string> grad_file;
  FmribOption<bool> save_bvals;
  bool parse_command_line(int argc, char** argv);
  
 private:
  dtifitOptions();  
  const dtifitOptions& operator=(dtifitOptions&);
  dtifitOptions(dtifitOptions&);

  OptionParser options; 
      
  static dtifitOptions* gopt;
  
};

 inline dtifitOptions& dtifitOptions::getInstance(){
   if(gopt == NULL)
     gopt = new dtifitOptions();
   
   return *gopt;
 }

 inline dtifitOptions::dtifitOptions() :
  verbose(string("-V,--verbose"), false, 
	  string("switch on diagnostic messages"), 
	  false, no_argument),
   help(string("-h,--help"), false,
	string("display this message"),
	false, no_argument),
   dtidatafile(string("-k,--data"), string("data"),
	       string("dti data file"),
	       true, requires_argument),  
   ofile(string("-o,--out"), string("dti"),
	       string("Output basename"),
	       true, requires_argument),
   maskfile(string("-m,--mask"), string("mask"),
	    string("Bet binary mask file"),
	    true, requires_argument),
   bvecsfile(string("-r,--bvecs"), string("bvecs"),
	     string("b vectors file"),
	     true, requires_argument),  
   bvalsfile(string("-b,--bvals"), string("bvals"),
	     string("b values file"),
	     true, requires_argument), 
   cni(string("--cni"), string(""),
	     string("Input confound regressors"),
	     false, requires_argument), 
   sse(string("--sse"), false,
	     string("Output sum of squared errors"),
	     false, no_argument), 
   wls(string("-w,--wls"),false, 
             string("Fit the tensor with weighted least squares"), 
             false, no_argument),
   littlebit(string("--littlebit"), false, 
	     string("Only process small area of brain"), 
	     false, no_argument),
   savetensor(string("--save_tensor"), false, 
	     string("Save the elements of the tensor"), 
	     false, no_argument),
   z_min(string("-z,--zmin"), 0, 
	 string("min z"), 
	 false, requires_argument),
   z_max(string("-Z,--zmax"), 42, 
	 string("max z"), 
	 false, requires_argument),
   y_min(string("-y,--ymin"), 0, 
	 string("min y"), 
	 false, requires_argument),
   y_max(string("-Y,--ymax"), 128, 
	 string("max y"), 
	 false, requires_argument),
   x_min(string("-x,--xmin"), 0, 
	 string("min x"), 
	 false, requires_argument),
   x_max(string("-X,--xmax"), 128, 
	 string("max x"), 
	 false, requires_argument),
   grad_file(string("--gradnonlin"), string("gradnonlin"),
	     string("Gradient Nonlinearity Tensor file"),
	     false, requires_argument),
   save_bvals(string("--savebvals"), false,
	     string("Save 4D file with bvalues, corrected for gradient nonlinearities"),
	     false,  no_argument),
   options("dtifit", "dtifit -k <filename>\n dtifit --verbose\n")
   {
     
    
     try {
       options.add(verbose);
       options.add(help);
       options.add(dtidatafile);
       options.add(ofile);
       options.add(maskfile);
       options.add(bvecsfile);
       options.add(bvalsfile);
       options.add(cni);
       options.add(sse);
       options.add(wls);
       options.add(littlebit);
       options.add(savetensor);
       options.add(z_min);
       options.add(z_max);
       options.add(y_min);
       options.add(y_max);
       options.add(x_min);
       options.add(x_max);
       options.add(grad_file);
       options.add(save_bvals);
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





