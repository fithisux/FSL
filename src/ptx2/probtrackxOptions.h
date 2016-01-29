/*  probtrackxOptions.h

    Tim Behrens, Saad Jbabdi, FMRIB Image Analysis Group

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

using namespace Utilities;

namespace TRACT {

class probtrackxOptions {
 public:
  static probtrackxOptions& getInstance();
  ~probtrackxOptions() { delete gopt; }
  
  Option<int>              verbose;
  Option<bool>             help;

  Option<string>           basename;
  Option<string>           outfile;
  Option<string>           logdir; 
  Option<bool>             forcedir;

  Option<string>           maskfile;
  Option<string>           seedfile; 

  Option<bool>             simple;
  Option<bool>             network;
  Option<bool>             simpleout;
  Option<bool>             pathdist;
  Option<bool>             omeanpathlength; 
  Option<string>           pathfile;
  Option<bool>             s2tout;
  Option<bool>             s2tastext;
  Option<bool>		   closestvertex;

  Option<string>           targetfile;
  Option<string>           waypoints;
  Option<string>           waycond;
  Option<bool>             wayorder;
  Option<bool>             onewaycondition;  // apply waypoint conditions to half the tract
  Option<string>           rubbishfile;
  Option<string>           stopfile;
  Option<string>           wtstopfiles;

  Option<bool>             matrix1out;
  Option<float>            distthresh1;  
  Option<bool>             matrix2out;
  Option<string>           lrmask;
  Option<bool>             matrix3out;
  Option<string>           mask3;
  Option<string>           lrmask3;
  Option<float>            distthresh3;  
  Option<bool>             matrix4out;
  Option<string>           mask4;
  Option<string>           dtimask;

  Option<string>           seeds_to_dti;
  Option<string>           dti_to_seeds;
  Option<string>           seedref;
  Option<string>           meshspace;

  Option<int>              nparticles;
  Option<int>              nsteps;
  Option<float>            steplength;


  Option<float>            distthresh;  
  Option<float>            c_thr;
  Option<float>            fibthresh;
  Option<bool>             loopcheck;
  Option<bool>             usef;
  Option<bool>             modeuler;

  Option<float>            sampvox;
  Option<int>              randfib;
  Option<int>              fibst;
  Option<int>              rseed;

  // hidden options
  FmribOption<string>      prefdirfile;      // inside this mask, pick orientation closest to whatever is in here
  FmribOption<string>      skipmask;         // inside this mask, ignore data (inertia)
  FmribOption<bool>        forcefirststep;   // always take at least one step 
  FmribOption<bool>        osampfib;         // not yet
  FmribOption<bool>        onewayonly;       // in surface mode, track towards the brain (assumes surface normal points towards the brain)
  FmribOption<bool>        opathdir;         // like fdt_paths but with average local tract orientation
  FmribOption<bool>        save_paths;       // save paths to ascii file
  FmribOption<string>      locfibchoice;     // inside this mask, define local rules for fibre picking
  FmribOption<string>      loccurvthresh;    // inside this mask, define local curvature threshold
  FmribOption<bool>        targetpaths;      // output separate fdt_paths for each target
  FmribOption<bool>        noprobinterpol;   // turn off probabilistic interpolation

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
	string("Display this message\n\n"),
	false, no_argument),

   basename(string("-s,--samples"),"",
	    string("Basename for samples files - e.g. 'merged'"),
	    true, requires_argument),  

   outfile(string("-o,--out"), string("fdt_paths"),
	   string("Output file (default='fdt_paths')"),
	   false, requires_argument),
   logdir(string("--dir"), string("logdir"),
	  string("\tDirectory to put the final volumes in - code makes this directory - default='logdir'"),
	  false, requires_argument),
   forcedir(string("--forcedir"), false,
	    string("Use the actual directory name given - i.e. don't add + to make a new directory\n\n"),
	    false, no_argument),

   maskfile(string("-m,--mask"),"",
	    string("Bet binary mask file in diffusion space"),
	    true, requires_argument),   
   seedfile(string("-x,--seed"),"",
	    string("Seed volume or list (ascii text file) of volumes and/or surfaces"),
	    true, requires_argument),

   simple(string("--simple"),false,
	string("Track from a list of voxels (seed must be a ASCII list of coordinates)"),
	false, no_argument),
   network(string("--network"), false,
	   string("Activate network mode - only keep paths going through at least one of the other seed masks"),
	   false, no_argument),
   simpleout(string("--opd"), false,
	     string("\tOutput path distribution"),
	     false, no_argument), 
   pathdist(string("--pd"), false,
	    string("\tCorrect path distribution for the length of the pathways"),
	    false, no_argument), 
   omeanpathlength(string("--ompl"), false,
	    string("\tOutput mean path length from seed"),
	    false, no_argument),
   pathfile(string("--fopd"), "",
	    string("\tOther mask for binning tract distribution"),
	    false, requires_argument), 
   s2tout(string("--os2t"), false,
	  string("\tOutput seeds to targets"),
	  false, no_argument),
   s2tastext(string("--s2tastext"), false,
	     string("Output seed-to-target counts as a text file (default in simple mode)\n\n"),
	     false, no_argument), 
   closestvertex(string("--closestvertex"), false,
	     string("Count only nearest neighbour vertex when the face of a surface is crossed."),
	     false, no_argument), 

   targetfile(string("--targetmasks"),"",
	      string("File containing a list of target masks - for seeds_to_targets classification"),
	      false, requires_argument),
   waypoints(string("--waypoints"), string(""),
	     string("Waypoint mask or ascii list of waypoint masks - only keep paths going through ALL the masks"),
	     false, requires_argument),
   waycond(string("--waycond"),"AND",
	   string("Waypoint condition. Either 'AND' (default) or 'OR'"),
	   false, requires_argument),
   wayorder(string("--wayorder"),false,
	    string("Reject streamlines that do not hit waypoints in given order. Only valid if waycond=AND"),
	    false,no_argument),
   onewaycondition(string("--onewaycondition"),false,
	    string("Apply waypoint conditions to each half tract separately"),
	    false, no_argument),
   rubbishfile(string("--avoid"), string(""),
	       string("\tReject pathways passing through locations given by this mask"),
	       false, requires_argument),
   stopfile(string("--stop"), string(""),
	       string("\tStop tracking at locations given by this mask file"),
	       false, requires_argument),
   wtstopfiles(string("--wtstop"), string(""),
	       string("One mask or text file with mask names. Allow propagation within mask but terminate on exit. If multiple masks, non-overlapping volumes expected\n\n"),
	       false, requires_argument),

   matrix1out(string("--omatrix1"), false,
	      string("Output matrix1 - SeedToSeed Connectivity"),
	      false, no_argument), 
   distthresh1(string("--distthresh1"), 0,
	       string("Discards samples (in matrix1) shorter than this threshold (in mm - default=0)"),
	       false, requires_argument),
   matrix2out(string("--omatrix2"), false,
	      string("Output matrix2 - SeedToLowResMask"),
	      false, no_argument), 
   lrmask(string("--target2"), string(""),
	  string("Low resolution binary brain mask for storing connectivity distribution in matrix2 mode"),
	  false, requires_argument),
   matrix3out(string("--omatrix3"), false,
	      string("Output matrix3 (NxN connectivity matrix)"),
	      false, no_argument), 
   mask3(string("--target3"), "",
	 string("Mask used for NxN connectivity matrix (or Nxn if lrtarget3 is set)"),
	 false, requires_argument), 
   lrmask3(string("--lrtarget3"), "",
	 string("Column-space mask used for Nxn connectivity matrix"),
	 false, requires_argument), 
   distthresh3(string("--distthresh3"), 0,
	       string("Discards samples (in matrix3) shorter than this threshold (in mm - default=0)"),
	       false, requires_argument),
   matrix4out(string("--omatrix4"), false,
	      string("Output matrix4 - DtiMaskToSeed (special Oxford Sparse Format)"),
	      false, no_argument), 
   mask4(string("--colmask4"), string(""),
	 string("Mask for columns of matrix4 (default=seed mask)"),
	 false, requires_argument),
   dtimask(string("--target4"), string(""),
	  string("Brain mask in DTI space\n\n"),
	  false, requires_argument),

   seeds_to_dti(string("--xfm"),"",
		string("\tTransform taking seed space to DTI space (either FLIRT matrix or FNIRT warpfield) - default is identity"),
		false, requires_argument),
   dti_to_seeds(string("--invxfm"), string(""),
		string("Transform taking DTI space to seed space (compulsory when using a warpfield for seeds_to_dti)"),
		false, requires_argument),
   seedref(string("--seedref"),"",
	   string("Reference vol to define seed space in simple mode - diffusion space assumed if absent"),
	   false, requires_argument),
   meshspace(string("--meshspace"), string("caret"),
	     string("Mesh reference space - either 'caret' (default) or 'freesurfer' or 'first' or 'vox' \n\n"),
	     false, requires_argument),

   nparticles(string("-P,--nsamples"), 5000,
	      string("Number of samples - default=5000"),
	      false, requires_argument),
   nsteps(string("-S,--nsteps"), 2000,
	  string("Number of steps per sample - default=2000"),
	  false, requires_argument),
   steplength(string("--steplength"), 0.5, 
	      string("Steplength in mm - default=0.5\n\n"), 
	      false, requires_argument),

   distthresh(string("--distthresh"), 0,
	      string("Discards samples shorter than this threshold (in mm - default=0)"),
	      false, requires_argument),
   c_thr(string("-c,--cthr"), 0.2, 
	 string("Curvature threshold - default=0.2"), 
	 false, requires_argument),
   fibthresh(string("--fibthresh"), 0.01, 
	     string("Volume fraction before subsidary fibre orientations are considered - default=0.01"), 
	     false, requires_argument),
   loopcheck(string("-l,--loopcheck"), false, 
	     string("Perform loopchecks on paths - slower, but allows lower curvature threshold"), 
	     false, no_argument),
   usef(string("-f,--usef"), false, 
	 string("Use anisotropy to constrain tracking"), 
	 false, no_argument),
   modeuler(string("--modeuler"), false, 
	    string("Use modified euler streamlining\n\n"), 
	    false, no_argument),


   sampvox(string("--sampvox"), 0, 
	   string("Sample random points within a sphere with radius x mm from the center of the seed voxels (e.g. --sampvox=0.5, 0.5 mm radius sphere). Default=0"), 
	   false, requires_argument),
   randfib(string("--randfib"), 0, 
	   string("Default 0. Set to 1 to randomly sample initial fibres (with f > fibthresh). \n                        Set to 2 to sample in proportion fibres (with f>fibthresh) to f. \n                        Set to 3 to sample ALL populations at random (even if f<fibthresh)"), 
	   false, requires_argument),
   fibst(string("--fibst"),1, 
	 string("\tForce a starting fibre for tracking - default=1, i.e. first fibre orientation. Only works if randfib==0"), 
	 false, requires_argument),
   rseed(string("--rseed"), 12345,
	 string("\tRandom seed"),
	 false, requires_argument), 


   prefdirfile(string("--prefdir"), string(""),
	       string("Prefered orientation preset in a 4D mask"),
	       false, requires_argument),
   skipmask(string("--no_integrity"), string(""),
	    string("No explanation needed"),
	    false, requires_argument),
   forcefirststep(string("--forcefirststep"),false,
		  string("In case seed and stop masks are the same. Also needed for surface wt_stop masks that are also seeds."),
		  false, no_argument),
   osampfib(string("--osampfib"),false,
	    string("Output sampled fibres"),
	    false, no_argument),
   onewayonly(string("--onewayonly"),false,
	      string("Track in one direction only (towards the brain - only valid for surface seeds)"),
	      false, no_argument),
   opathdir(string("--opathdir"),false,
	    string("Output average local tract orientation (tangent)"),
	    false, no_argument),
   save_paths(string("--savepaths"),false,
	      string("Save path coordinates to ascii file"),
	      false, no_argument),
   locfibchoice(string("--locfibchoice"),string(""),
	      string("Local rules for fibre choice - 0=closest direction(default), 1=equal prob (f>thr), 2=equal prob with angle threshold (=40 deg), 3=sample in proportion to f"),
	      false, requires_argument),
   loccurvthresh(string("--loccurvthresh"),string(""),
	      string("Local curvature threshold"),
	      false, requires_argument),
   targetpaths(string("--otargetpaths"),false,
	      string("Output separate fdt_paths for targets (assumes --os2t is on)"),
	      false, no_argument),
   noprobinterpol(string("--noprobinterpol"),false,
	      string("Turn off probabilistic interpolation"),
	      false, no_argument),

   options("probtrackx","probtrackx2 -s <basename> -m <maskname> -x <seedfile> -o <output> --targetmasks=<textfile>\n probtrackx2 --help\n")
   {
     
    
     try {
       options.add(verbose);
       options.add(help);

       options.add(basename);
       options.add(outfile);
       options.add(logdir); 
       options.add(forcedir); 
 
       options.add(maskfile);
       options.add(seedfile); 
 
       options.add(simple);
       options.add(network);
       options.add(simpleout);
       options.add(pathdist);
       options.add(omeanpathlength);
       options.add(pathfile);
       options.add(s2tout);
       options.add(s2tastext);
       options.add(closestvertex);

       options.add(targetfile);
       options.add(waypoints);
       options.add(waycond);
       options.add(wayorder);
       options.add(onewaycondition);
       options.add(rubbishfile);
       options.add(stopfile);
       options.add(wtstopfiles);

       options.add(matrix1out);
       options.add(distthresh1);
       options.add(matrix2out);
       options.add(lrmask);
       options.add(matrix3out);
       options.add(mask3);
       options.add(lrmask3);
       options.add(distthresh3);
       options.add(matrix4out);
       options.add(mask4);
       options.add(dtimask);

       options.add(seeds_to_dti);
       options.add(dti_to_seeds);
       options.add(seedref);
       options.add(meshspace);


       options.add(nparticles);
       options.add(nsteps);
       options.add(steplength);

       options.add(distthresh);
       options.add(c_thr);
       options.add(fibthresh);
       options.add(loopcheck);
       options.add(usef);
       options.add(modeuler);

       options.add(sampvox);
       options.add(randfib);
       options.add(fibst);
       options.add(rseed);

       options.add(skipmask);
       options.add(prefdirfile);
       options.add(forcefirststep);
       options.add(osampfib);
       options.add(onewayonly);
       options.add(opathdir);
       options.add(save_paths);
       options.add(locfibchoice);
       options.add(loccurvthresh);
       options.add(targetpaths);
       options.add(noprobinterpol);

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







