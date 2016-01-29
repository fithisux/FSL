/*  probtrackxOptions.h

    Tim Behrens, FMRIB Image Analysis Group

    Copyright (C) 2004 University of Oxford  */

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

#if !defined(probtrackxOptions_h)
#define probtrackxOptions_h

#include <string> 
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "utils/log.h"
#include "commonopts.h"

//#include "newmatall.h"
using namespace Utilities;

namespace TRACT {

class probtrackxOptions {
 public:
  static probtrackxOptions& getInstance();
  ~probtrackxOptions() { delete gopt; }
  
  Option<int> verbose;
  Option<bool> help;
  Option<string> basename;
  Option<string> maskfile;
  Option<string> seedfile; 
  Option<string> mode;
  Option<string> targetfile;
  Option<bool> simpleout;
  Option<bool> pathdist;
  Option<bool> s2tout;
  FmribOption<bool> matrix1out;
  FmribOption<bool> matrix2out;
  FmribOption<bool> matrix3out;
  FmribOption<string> maskmatrix3;
  FmribOption<bool> maskmatrixout;
  Option<string> outfile;
  Option<string> rubbishfile;
  Option<string> stopfile;
  FmribOption<string> prefdirfile;
  Option<string> seeds_to_dti;
  Option<string> dti_to_seeds;
  FmribOption<string> skipmask;
  Option<string> seedref;
  Option<string> mask2;
  Option<string> waypoints;
  Option<bool> network;
  Option<string> meshfile;
  FmribOption<string> meshspace;
  FmribOption<string> lrmask;
  Option<string> logdir; 
  Option<bool> forcedir;
  Option<int> nparticles;
  Option<int> nsteps;
  Option<float> distthresh;
  FmribOption<float> distthresh3;
  Option<float> c_thr;
  Option<float> fibthresh;
  Option<bool> sampvox;
  Option<float> steplength;
  Option<bool> loopcheck;
  Option<bool> usef;
  Option<int> randfib;
  Option<int> fibst;
  Option<bool> modeuler;
  Option<int> rseed;
  Option<bool> seedcountastext;
  FmribOption<bool> splitmatrix2;


  void parse_command_line(int argc, char** argv,Log& logger);
  void modecheck();
  void modehelp();
  void matrixmodehelp();
  void status();
 private:
  probtrackxOptions();  
  const probtrackxOptions& operator=(probtrackxOptions&);
  probtrackxOptions(probtrackxOptions&);

  OptionParser options; 
      
  static probtrackxOptions* gopt;
  
};


 inline probtrackxOptions& probtrackxOptions::getInstance(){
   if(gopt == NULL)
     gopt = new probtrackxOptions();
   
   return *gopt;
 }

 inline probtrackxOptions::probtrackxOptions() :
  verbose(string("-V,--verbose"), 0, 
	  string("Verbose level, [0-2]"), 
	  false, requires_argument),
   help(string("-h,--help"), false,
	string("Display this message"),
	false, no_argument),
   basename(string("-s,--samples"), string("merged"),
	       string("Basename for samples files"),
	       true, requires_argument),  
   maskfile(string("-m,--mask"), string("mask"),
	    string("Bet binary mask file in diffusion space"),
	    true, requires_argument),
   seedfile(string("-x,--seed"), string("Seed"),
	    string("Seed volume, or voxel, or ascii file with multiple volumes, or freesurfer label file"),
	    true, requires_argument),
   mode(string("--mode"), string(""),
	string("\tUse --mode=simple for single seed voxel"),
	    false, requires_argument),
  targetfile(string("--targetmasks"), string("cmasks"),
	    string("File containing a list of target masks - required for seeds_to_targets classification"),
	    false, requires_argument),
  simpleout(string("--opd"), false,
	    string("\tOutput path distribution"),
	    false, no_argument), 
  pathdist(string("--pd"), false,
	   string("\tCorrect path distribution for the length of the pathways"),
	   false, no_argument), 
  s2tout(string("--os2t"), false,
	 string("\tOutput seeds to targets"),
	 false, no_argument),
  matrix1out(string("--omatrix1"), false,
	  string("Output matrix1"),
	  false, no_argument), 
  matrix2out(string("--omatrix2"), false,
	  string("Output matrix2"),
	  false, no_argument), 
  matrix3out(string("--omatrix3"), false,
	  string("Output matrix3 (NxN connectivity matrix)"),
	  false, no_argument), 
  maskmatrix3(string("--mask3"), "",
	  string("Mask used for NxN connectivity matrix"),
	  false, requires_argument), 
  maskmatrixout(string("--omaskmatrix"), false,
		string("Output maskmatrix"),
		false, no_argument), 
   outfile(string("-o,--out"), string("fdt_paths"),
	   string("Output file (default='fdt_paths')"),
	   false, requires_argument),
   rubbishfile(string("--avoid"), string(""),
	       string("\tReject pathways passing through locations given by this mask"),
	       false, requires_argument),
   stopfile(string("--stop"), string(""),
	       string("\tStop tracking at locations given by this mask file"),
	       false, requires_argument),
   prefdirfile(string("--prefdir"), string(""),
	       string("Prefered orientation preset in a 4D mask"),
	       false, requires_argument),
   seeds_to_dti(string("--xfm"), string(""),
		string("\tTransform taking seed space to DTI space (either FLIRT matrix or FNIRT warpfield) - default is identity"),
		false, requires_argument),
   dti_to_seeds(string("--invxfm"), string(""),
		string("Transform taking DTI space to seed space (compulsory when using a warpfield for seeds_to_dti)"),
		false, requires_argument),
   skipmask(string("--no_integrity"), string(""),
	   string("no explanation needed"),
	   false, requires_argument),
  seedref(string("--seedref"), string(""),
	 string("Reference vol to define seed space in simple mode - diffusion space assumed if absent"),
	 false, requires_argument),
  mask2(string("--mask2"), string(""),
	 string("\tSecond mask in twomask_symm mode."),
       false, requires_argument),
 waypoints(string("--waypoints"), string(""),
	 string("Waypoint mask or ascii list of waypoint masks - only keep paths going through ALL the masks"),
       false, requires_argument),
 network(string("--network"), false,
	 string("Activate network mode - only keep paths going through at least one seed mask (required if multiple seed masks)"),
       false, no_argument),
   meshfile(string("--mesh"), string(""),
	 string("\tFreesurfer-type surface descriptor (in ascii format)"),
       false, requires_argument),
   meshspace(string("--meshspace"), string("freesurfer"),
	 string("Mesh reference space - either 'freesurfer' (default) or 'caret' or 'first'"),
       false, requires_argument),
  lrmask(string("--lrmask"), string(""),
	 string("Low resolution binary brain mask for storing connectivity distribution in matrix2 mode"),
       false, requires_argument),
  logdir(string("--dir"), string("logdir"),
	    string("\tDirectory to put the final volumes in - code makes this directory - default='logdir'"),
	    false, requires_argument),
  forcedir(string("--forcedir"), false,
	 string("Use the actual directory name given - i.e. don't add + to make a new directory"),
	 false, no_argument),
  nparticles(string("-P,--nsamples"), 5000,
	 string("Number of samples - default=5000"),
	 false, requires_argument),
   nsteps(string("-S,--nsteps"), 2000,
	    string("Number of steps per sample - default=2000"),
	    false, requires_argument),
   distthresh(string("--distthresh"), 0,
	    string("Discards samples shorter than this threshold (in mm - default=0)"),
	    false, requires_argument),
   distthresh3(string("--distthresh3"), 0,
	    string("Discards samples (in matrix3) shorter than this threshold (in voxels - default=0)"),
	    false, requires_argument),
   c_thr(string("-c,--cthr"), 0.2, 
	 string("Curvature threshold - default=0.2"), 
	 false, requires_argument),
  fibthresh(string("--fibthresh"), 0.01, 
	    string("Volume fraction before subsidary fibre orientations are considered - default=0.01"), 
	 false, requires_argument),
  sampvox(string("--sampvox"), false, 
	    string("Sample random points within seed voxels"), 
	 false, no_argument),
   steplength(string("--steplength"), 0.5, 
	 string("Steplength in mm - default=0.5"), 
	 false, requires_argument),
   loopcheck(string("-l,--loopcheck"), false, 
	 string("Perform loopchecks on paths - slower, but allows lower curvature threshold"), 
	 false, no_argument),
   usef(string("-f,--usef"), false, 
	 string("Use anisotropy to constrain tracking"), 
	 false, no_argument),
  randfib(string("--randfib"), 0, 
	 string("Default 0. Set to 1 to randomly sample initial fibres (with f > fibthresh). \n                        Set to 2 to sample in proportion fibres (with f>fibthresh) to f. \n                        Set to 3 to sample ALL populations at random (even if f<fibthresh)"), 
	 false, requires_argument),
  fibst(string("--fibst"),1, 
	 string("\tForce a starting fibre for tracking - default=1, i.e. first fibre orientation. Only works if randfib==0"), 
	 false, requires_argument),
  modeuler(string("--modeuler"), false, 
	   string("Use modified euler streamlining"), 
	   false, no_argument),
  rseed(string("--rseed"), 12345,
	string("\tRandom seed"),
	false, requires_argument), 
  seedcountastext(string("--s2tastext"), false,
		  string("Output seed-to-target counts as a text file (useful when seeding from a mesh)"),
		  false, no_argument), 
  splitmatrix2(string("--splitmatrix2"), false,
		  string("Split matrix 2 along seed dimension (in case it is too large)"),
		  false, no_argument), 
   options("probtrackx","probtrackx -s <basename> -m <maskname> -x <seedfile> -o <output> --targetmasks=<textfile>\n probtrackx --help\n")
   {
     
    
     try {
       options.add(verbose);
       options.add(help);
       options.add(basename);
       options.add(maskfile);
       options.add(seedfile); 
       options.add(mode);
       options.add(targetfile);
       options.add(skipmask);
       options.add(mask2);
       options.add(waypoints);
       options.add(network);
       options.add(meshfile);
       options.add(meshspace);
       options.add(lrmask);
       options.add(seedref);
       options.add(logdir); 
       options.add(forcedir); 
       options.add(simpleout);
       options.add(pathdist);
       options.add(s2tout);
       options.add(matrix1out);
       options.add(matrix2out);
       options.add(matrix3out);
       options.add(maskmatrix3);
       options.add(maskmatrixout);
       options.add(outfile);
       options.add(rubbishfile);
       options.add(stopfile);
       options.add(prefdirfile);
       options.add(seeds_to_dti);
       options.add(dti_to_seeds);
       options.add(nparticles);
       options.add(nsteps);
       options.add(distthresh);
       options.add(distthresh3);
       options.add(c_thr);
       options.add(fibthresh);
       options.add(sampvox);
       options.add(steplength);
       options.add(loopcheck);
       options.add(usef);
       options.add(randfib);
       options.add(fibst);
       options.add(modeuler);
       options.add(rseed);
       options.add(seedcountastext);
       options.add(splitmatrix2);
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







