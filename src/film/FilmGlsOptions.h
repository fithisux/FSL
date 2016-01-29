/*  FilmGlsOptions.h

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

#if !defined(__FilmGlsOptions_h)
#define __FilmGlsOptions_h
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "utils/log.h"

using namespace Utilities;

namespace FILM {

class FilmGlsOptions {
 public:
  static FilmGlsOptions& getInstance();
  ~FilmGlsOptions() { delete gopt; }
  
  Option<bool> ac_only;
  Option<float> thresh;  
  Option<bool> fitAutoRegressiveModel;
  Option<bool> help;
  Option<bool> noest;
  Option<bool> output_pwdata;
  Option<bool> pava;
  Option<bool> smoothACEst;
  Option<bool> verbose;
  Option<string> datadir;
  Option<string> analysisMode;
  Option<string> inputDataName;
  Option<string> inputDataName2;
  Option<string> meanInputFile;
  Option<string> minimumTimepointFile;
  Option<string> paradigmfname;
  Option<string> contrastFile;
  Option<string> fContrastFile;
  Option<int> epith;
  Option<int> ms;
  Option<int> tukeysize;
  Option<float> multitapersize;
  Option<vector<int> > voxelwise_ev_numbers;
  Option<vector<string> > voxelwiseEvFilenames;

  void parse_command_line(int argc, char** argv,Log& logger);
  
private:
  FilmGlsOptions();
  const FilmGlsOptions& operator=(FilmGlsOptions&);
  FilmGlsOptions(FilmGlsOptions&);
  OptionParser options; 
  static FilmGlsOptions* gopt;
};

inline FilmGlsOptions& FilmGlsOptions::getInstance() {
   if(gopt == NULL)
     gopt = new FilmGlsOptions();
   return *gopt;
}

 
inline FilmGlsOptions::FilmGlsOptions() :
  ac_only(string("--ac"), false,
        string("\tperform autocorrelation estimation only"),
        false, no_argument),
  thresh(string("--thr"),0,
	 string("~<num>\tinitial threshold to apply to input data"),
	 false, requires_argument),
  fitAutoRegressiveModel(string("--ar"), false,
	string("\tfits autoregressive model - default is to use tukey with M=sqrt(numvols)"),
        false, no_argument),
  help(string("--help"), false,
	string("\tprints this message"),
        false, no_argument),
  noest(string("--noest"), false,
	string("\tdo not estimate auto corrs"),
	false, no_argument),
  output_pwdata(string("--outputPWdata"), false,
	string("\toutput prewhitened data and average design matrix"),
	false, no_argument),
  pava(string("--pava"), false,
	string("\testimates autocorr using PAVA - default is to use tukey with M=sqrt(numvols)"),
	false, no_argument),
  smoothACEst(string("--sa"), false,
	string("\tsmooths auto corr estimates"),
	false, no_argument),
  verbose(string("-v"), false,
	string("\toutputs full data"),
	false, no_argument),
  datadir(string("--rn"), string("results"),
       string("~<file>\tdirectory name to store results in, default is results"),
       false, requires_argument),
  analysisMode(string("--mode"), string("volumetric"),
       string("~<mode>\tanalysis mode, options are volumetric ( default ) or surface. Caution: surface-based functionality is still BETA"), false, requires_argument),
  inputDataName(string("--in"), string(""),string("~<file>\tinput data file ( NIFTI for volumetric, GIFTI for surface )"),true, requires_argument),
  inputDataName2(string("--in2"), string(""),string("~<file>\tinput surface for autocorr smoothing in surface-based analyses"),false, requires_argument),
   meanInputFile(string("--mf"), string(""),
       string("~<file>\tre-estimate mean_func baseline - for use with perfusion subtraction"),
       false, requires_argument), 
   minimumTimepointFile(string("--mft"), string(""),
       string("~<file>\tminimum timepoint file"),
       false, requires_argument),
   paradigmfname(string("--pd"), string(""),
       string("~<file>\tparadigm file"),
        false, requires_argument),
   contrastFile(string("--con"), string(""),
       string("~<file>\tt-contrasts file"),
        false, requires_argument),
   fContrastFile(string("--fcon"), string(""),
       string("~<file>\tf-contrasts file"),
        false, requires_argument),

   epith(string("--epith"), 0,
        string("~<num>\tsusan brightness threshold for volumetric analysis/smoothing sigma for surface analysis"),
        false, requires_argument),
   ms(string("--ms"), 4,
        string("~<num>\tsusan mask size for volumetric analysis/smoothing extent for surface analysis"),
        false, requires_argument),
   tukeysize(string("--tukey"), 0,
        string("~<num>\tuses tukey window to estimate autocorr with window size num - default is to use tukey with M=sqrt(numvols)"),
        false, requires_argument),
   multitapersize(string("--mt"), -1,
        string("~<num>\tuses multitapering with slepian tapers and num is the time-bandwidth product - default is to use tukey with M=sqrt(numvols)"),
        false, requires_argument),
   voxelwise_ev_numbers(string("--ven"), vector<int>(), 
         string("\tlist of numbers indicating voxelwise EVs position in the design matrix (list order corresponds to files in vxf option). Caution BETA option, only use with volumetric analysis."), 
         false, requires_argument),
   voxelwiseEvFilenames(string("--vef"), vector<string>(), 
         string("\tlist of 4D images containing voxelwise EVs (list order corresponds to numbers in vxl option). Caution BETA option, only use with volumetric analysis."), 
			false, requires_argument),
  options("film_gls","film_gls")
   {
     try {
       options.add(ac_only);
       options.add(thresh);
       options.add(fitAutoRegressiveModel);
       options.add(help);
       options.add(noest);
       options.add(output_pwdata);
       options.add(pava);
       options.add(smoothACEst);
       options.add(verbose);
       options.add(datadir);      
       options.add(analysisMode);
       options.add(inputDataName);
       options.add(inputDataName2);
       options.add(meanInputFile);
       options.add(minimumTimepointFile);
       options.add(paradigmfname);
       options.add(contrastFile);
       options.add(fContrastFile); 
       options.add(epith);
       options.add(ms);
       options.add(tukeysize);
       options.add(multitapersize);
       options.add(voxelwise_ev_numbers);
       options.add(voxelwiseEvFilenames);
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










