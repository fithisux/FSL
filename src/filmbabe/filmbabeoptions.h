/*  FilmbabeOptions.h

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

#if !defined(FilmbabeOptions_h)
#define FilmbabeOptions_h

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "utils/log.h"

using namespace Utilities;

namespace Filmbabe {

class FilmbabeOptions {
 public:
  static FilmbabeOptions& getInstance();
  ~FilmbabeOptions() { delete gopt; }
  
  Option<bool> verbose;
  Option<int> debuglevel;
  Option<bool> timingon;
  Option<bool> help;
  Option<bool> flobsprioroff;
  Option<string> datafile;
  Option<string> maskfile;
  Option<string> designfile;
  Option<string> flobsregressorsfile;
  Option<string> flobsdir;
  Option<string> priorcovarfile;
  Option<string> priormeanfile;
  Option<string> logdir;
  Option<int> niters;
  Option<float> tarmrfprec;
  Option<bool> tarard;
  Option<int> ntracesamps;
  Option<int> ntar;

  void parse_command_line(int argc, char** argv, Log& logger);
  
 private:
  FilmbabeOptions();  
  const FilmbabeOptions& operator=(FilmbabeOptions&);
  FilmbabeOptions(FilmbabeOptions&);

  OptionParser options; 
      
  static FilmbabeOptions* gopt;
  
};

 inline FilmbabeOptions& FilmbabeOptions::getInstance(){
   if(gopt == NULL)
     gopt = new FilmbabeOptions();
   
   return *gopt;
 }

 inline FilmbabeOptions::FilmbabeOptions() :
   verbose(string("-V,--verbose"), false, 
	   string("switch on diagnostic messages"), 
	   false, no_argument),
   debuglevel(string("--db,--debug,--debuglevel"), 0,
		       string("set debug level"), 
		       false, requires_argument),
   timingon(string("--to,--timingon"), false, 
		       string("turn timing on"), 
		       false, no_argument),
   help(string("-h,--help"), false,
		    string("display this message"),
		    false, no_argument),   
   flobsprioroff(string("--fpo,--flobsprioroff"), false,
		    string("Turn FLOBS prior off"),
		    false, no_argument),
   datafile(string("--df,--datafile"), string(""),
			  string("data file"),
		     true, requires_argument),  
   maskfile(string("-m,--mask"), string(""),
	    string("mask file"),
	    true, requires_argument),
   designfile(string("-d,--dm,--designfile"), string(""),
	      string("design matrix file"),
	      true, requires_argument),
   flobsregressorsfile(string("--frf"), string(""),
	      string("file indicating which regressors belong to which original ev design matrix file (a -1 label indicates a non-flobs regerssor)"),
	      true, requires_argument),
   flobsdir(string("--fd"), string(""),
	      string("flobs directory, when using flobs constraints need this or to explicitly specify --pcf and--pmf"),
	    false, requires_argument),
   priorcovarfile(string("--pcf,--priorcovarfile"), string(""),
	      string("priorcovar matrix file"),
	      false, requires_argument),
   priormeanfile(string("--pmf,--priormeanfile"), string(""),
	      string("priormean matrix file"),
	      false, requires_argument),
   logdir(string("-l,--ld,--logdir"), string("logdir"),
	  string("log directory"),
	  false, requires_argument),  
   niters(string("--ni"), 5,
	  string("Num pf VB iterations (default is 5)"),
	  false, requires_argument),
   tarmrfprec(string("--tmp,--tarmrfprec"), -1, 
	      string("MRF precision to impose on temporal AR maps, default is -1 for a proper full Bayes approach"),
	      false, requires_argument),
   tarard(string("--tarard"), false, 
	  string("Impose ARD/MRF on temporal AR"),	       
	  false, no_argument),
   ntracesamps(string("--nts,--ntracesamps"), 0,
	       string("No of samples to take to estimate trace (default is 0, which uses only diagonal elements of the precision matrix to estimate trace)"),
	       false, requires_argument),
   ntar(string("--ntar"), 3,
	string("Order of temporal AR (default is 3)"),
	false, requires_argument),
   options("filmbabe", "filmbabe -df <filename>\n filmbabe --verbose\n")
 {
     try {
       options.add(verbose);
       options.add(debuglevel);
       options.add(timingon);
       options.add(help);
       options.add(flobsprioroff);
       options.add(datafile);
       options.add(maskfile);
       options.add(designfile);
       options.add(flobsregressorsfile);
       options.add(flobsdir);
       options.add(priorcovarfile);
       options.add(priormeanfile);
       options.add(logdir);
       options.add(niters);
       options.add(tarmrfprec);
       options.add(tarard);
       options.add(ntracesamps);
       options.add(ntar);

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



