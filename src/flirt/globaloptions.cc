/*  globaloptions.cc

    Mark Jenkinson, FMRIB Image Analysis Group

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

#include "globaloptions.h"

globaloptions* globaloptions::gopt = NULL;


void globaloptions::parse_command_line(int argc,char** argv,
				       const string &p_version)
{
  version = p_version;

  if(argc<2){
    print_usage(argc,argv);
    exit(1);
  }


  int n=1;
  string arg;
  char first;

  while (n<argc) {
    arg=argv[n];
    if (arg.size()<1) { n++; continue; }
    first = arg[0];
    if (first!='-') {
      cerr << "Unrecognised option " << first << endl;
      exit(-1);
    }
    
    // put options without arguments here
    if ( arg == "-help" ) {
      print_usage(argc,argv);
      exit(0);
    } else if ( arg == "-version") {
      print_version();
      exit(0);
    } else if ( arg == "-applyxfm" || arg == "-applynonisoxfm" ) {
      do_optimise = false;
      iso = false;
      nosave = false;
      n++;
      continue;
    } else if ( arg == "-i") {
      interactive = true;
      n++;
      continue;
    } else if ( arg == "-nosearch") {
      searchrx(1) = 0;
      searchrx(2) = 0;
      searchry(1) = 0;
      searchry(2) = 0;
      searchrz(1) = 0;
      searchrz(2) = 0;
      n++;
      continue;
    } else if ( arg == "-nosave") {
      nosave = true;
      n++;
      continue;
    } else if ( arg == "-noresample") {
      resample = false;
      n++;
      continue;
    } else if ( arg == "-forcescaling") {
      force_scaling = true;
      n++;
      continue;
    } else if ( arg == "-2D") {
      mode2D = true;
      n++;
      continue;
    } else if ( arg == "-noclamp") {
      clamping = false;
      n++;
      continue;
    } else if ( arg == "-noresampblur") {
      interpblur = false;
      n++;
      continue;
    } else if ( arg == "-usesqform") {
      initmatsqform = true;
      n++;
      continue;    
    } else if ( arg == "-printinit") {
      printinit = true;
      n++;
      continue;    
    } else if ( arg == "-debugsave") {
      nosave = false;
      n++;
      continue;
   } else if ( arg == "-debug" ) {
      debug = true;
      n++;
      continue;
   } else if ( arg == "-v" ) {
      verbose = 1;
      n++;
      continue;
    }

    if (n+1>=argc) 
      { 
	cerr << "Lacking argument to option " << arg << endl;
	exit(-1);
      }

    // put options with 1 argument here
    if ( (arg == "-o") || (arg == "-out") ) {
      outputfname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-ref") {
      reffname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-in") {
      inputfname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-init") {
      initmatfname = argv[n+1];
      initmatsqform = false;
      n+=2;
      continue;
    } else if ( arg == "-schedule") {
      schedulefname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-refweight") {
      refweightfname = argv[n+1];
      useweights = true;
      n+=2;
      continue;
    } else if ( arg == "-inweight") {
      testweightfname = argv[n+1];
      useweights = true;
      n+=2;
      continue;
    } else if ( arg == "-wmseg") {
      wmsegfname = argv[n+1];
      useseg = true;
      n+=2;
      continue;
    } else if ( arg == "-fieldmap") {
      fmapfname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-fieldmapmask") {
      fmapmaskfname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-wmcoords") {
      wmcoordsfname = argv[n+1];
      useseg = true;
      usecoords = true;
      n+=2;
      continue;
    } else if ( arg == "-wmnorms") {
      wmnormsfname = argv[n+1];
      useseg = true;
      usecoords = true;
      n+=2;
      continue;
    } else if ( arg == "-omat") {
      outputmatascii = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-bins") {
      no_bins = atoi(argv[n+1]);
      n+=2;
      continue;
    } else if ( arg == "-dof") {
      no_params = atoi(argv[n+1]);
      dof = atoi(argv[n+1]);
      n+=2;
      continue;
    } else if ( arg == "-sincwidth") {
      sincwidth = atoi(argv[n+1]);
      n+=2;
      continue;
    } else if ( arg == "-paddingsize") {
      paddingsize = atof(argv[n+1]);
      n+=2;
      continue;
    } else if ( arg == "-applyisoxfm" ) {
      isoscale = atof(argv[n+1]);
      do_optimise = false;
      iso = true;
      nosave = false;
      n+=2;
      continue;
    } else if ( arg == "-basescale") {
      force_basescale = true;
      basescale = atof(argv[n+1]);
      n+=2;
      continue;
    } else if ( arg == "-minsampling") {
      min_sampling = atof(argv[n+1]);
      n+=2;
      continue;
    } else if ( arg == "-setbackground") {
      backgndval = atof(argv[n+1]);
      forcebackgnd = true;
      n+=2;
      continue;
    } else if ( arg == "-coarsesearch") {
      coarsedelta = atof(argv[n+1])*M_PI/180.0;
      n+=2;
      continue;
    } else if ( arg == "-finesearch") {
      finedelta = atof(argv[n+1])*M_PI/180.0;
      n+=2;
      continue;
    } else if ( arg == "-echospacing") {
      echo_spacing = atof(argv[n+1]);
      n+=2;
      continue;
    } else if ( arg == "-pedir") {
      pe_dir = atoi(argv[n+1]);
      if (pe_dir==0) {
	cerr << "Unrecognised argument to pedir (" << argv[n+1] << ") - it should be a number between -3 and 3 (not 0)" << endl;
	exit(-1);
      }
      n+=2;
      continue;
    } else if ( arg == "-bbrtype") {
      bbr_type = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-bbrslope") {
      bbr_slope = atof(argv[n+1]);
      n+=2;
      continue;
    } else if ( arg == "-verbose") {
      verbose = atoi(argv[n+1]);
      n+=2;
      continue;
    } else if ( arg == "-datatype") {
      {
	forcedatatype = true;
	string dataarg = argv[n+1];
	if (dataarg == "double") {
	  datatype = DT_DOUBLE;
	} else if (dataarg == "float") {
	  datatype = DT_FLOAT;
	} else if (dataarg == "int") {
	  datatype = DT_SIGNED_INT;
	} else if (dataarg == "short") {
	  datatype = DT_SIGNED_SHORT;
	} else if (dataarg == "char") {
	  datatype = DT_UNSIGNED_CHAR;
	} else {
	  cerr << "Unrecognised data type: " << dataarg << endl;
	  exit(-1);
	}
      }
      n+=2;
      continue;
    } else if ( arg == "-cost") {
      {
	string costarg = argv[n+1];
	maincostfn = costfn_type(costarg);
	if (maincostfn == Unknown) {
	  cerr << "Unrecognised cost function type: " << costarg << endl;
	  exit(-1);
	}
      }
      n+=2;
      continue;
    } else if ( arg == "-searchcost") {
      {
	string costarg = argv[n+1];
	searchcostfn = costfn_type(costarg);
	if (searchcostfn == Unknown) {
	  cerr << "Unrecognised cost function type: " << costarg << endl;
	  exit(-1);
	}
      }
      n+=2;
      continue;
    } else if ( arg == "-interp") {
      {
	string interparg = argv[n+1];
	if (interparg == "trilinear") {
	  interpmethod = TriLinear;
	} else if (interparg == "nearestneighbour") {
	  interpmethod = NearestNeighbour;
	} else if (interparg == "sinc") {
	  interpmethod = NEWIMAGE::Sinc;
	} else if (interparg == "spline") {
	  interpmethod = NEWIMAGE::Spline;
	} else {
	  cerr << "Unrecognised interpolation method: " << interparg << endl;
	  exit(-1);
	}
      }
      n+=2;
      continue;
    } else if ( arg == "-sincwindow") {
      {
	string winarg = argv[n+1];
	if (winarg == "rectangular") {
	  sincwindow = Rect;
	} else if (winarg == "hanning") {
	  sincwindow = Hanning;
	} else if (winarg == "blackman") {
	  sincwindow = Blackman;
	} else {
	  cerr << "Unrecognised sinc window: " << winarg << endl;
	  exit(-1);
	}
      }
      n+=2;
      continue;
    } else if ( arg == "-anglerep" ) {
      {
	string anglearg = argv[n+1];
	if (anglearg == "quaternion") {
	  anglerep = Quaternion;
	} else if (anglearg == "euler") {
	  anglerep = Euler;
	} else {
	  cerr << "Unrecognised angle representation: " << anglearg << endl;
	  exit(-1);
	}
      }
      n+=2;
      continue;
    }

    if (n+2>=argc) 
      { 
	cerr << "Lacking argument to option " << arg << endl;
	exit(-1);
      }
    
    
    // put options with 2 arguments here
    if ( arg == "-searchrx" ) {
      searchrx(1) = Min(atof(argv[n+1]),atof(argv[n+2]))*M_PI/180.0;
      searchrx(2) = Max(atof(argv[n+1]),atof(argv[n+2]))*M_PI/180.0;
      n+=3;
      continue;
    } else if ( arg == "-searchry" ) {
      searchry(1) = Min(atof(argv[n+1]),atof(argv[n+2]))*M_PI/180.0;
      searchry(2) = Max(atof(argv[n+1]),atof(argv[n+2]))*M_PI/180.0;
      n+=3;
      continue;
    } else if ( arg == "-searchrz" ) {
      searchrz(1) = Min(atof(argv[n+1]),atof(argv[n+2]))*M_PI/180.0;
      searchrz(2) = Max(atof(argv[n+1]),atof(argv[n+2]))*M_PI/180.0;
      n+=3;
      continue;
    } else { 
      cerr << "Unrecognised option " << arg << endl;
      exit(-1);
    } 

    

  }  // while (n<argc)

  if (inputfname.size()<1) {
    cerr << "ERROR:: Input volume filename not found\n\n";
    print_usage(argc,argv);
    exit(2);
  }

  if (reffname.size()<1) {
    cerr << "ERROR:: Reference volume filename not found\n\n";
    print_usage(argc,argv);
    exit(2);
  }
}

void globaloptions::print_usage(int argc, char *argv[])
{
  print_version();
  cout << endl;
  cout << "Usage: " << argv[0] << " [options] -in <inputvol> -ref <refvol> -out <outputvol>\n"
       << "       " << argv[0] << " [options] -in <inputvol> -ref <refvol> -omat <outputmatrix>\n"
       << "       " << argv[0] << " [options] -in <inputvol> -ref <refvol> -applyxfm -init <matrix> -out <outputvol>\n\n"
       << "  Available options are:\n"
       << "        -in  <inputvol>                    (no default)\n"
       << "        -ref <refvol>                      (no default)\n"
       << "        -init <matrix-filname>             (input 4x4 affine matrix)\n"
       << "        -omat <matrix-filename>            (output in 4x4 ascii format)\n"
       << "        -out, -o <outputvol>               (default is none)\n"
       << "        -datatype {char,short,int,float,double}                    (force output data type)\n"
       << "        -cost {mutualinfo,corratio,normcorr,normmi,leastsq,labeldiff,bbr}        (default is corratio)\n"
       << "        -searchcost {mutualinfo,corratio,normcorr,normmi,leastsq,labeldiff,bbr}  (default is corratio)\n"
       << "        -usesqform                         (initialise using appropriate sform or qform)\n"
       << "        -displayinit                       (display initial matrix)\n"
       << "        -anglerep {quaternion,euler}       (default is euler)\n"
       << "        -interp {trilinear,nearestneighbour,sinc,spline}  (final interpolation: def - trilinear)\n"
       << "        -sincwidth <full-width in voxels>  (default is 7)\n"
       << "        -sincwindow {rectangular,hanning,blackman}\n"
       << "        -bins <number of histogram bins>   (default is "
                                        << no_bins << ")\n"
       << "        -dof  <number of transform dofs>   (default is "
                                        << dof << ")\n"
       << "        -noresample                        (do not change input sampling)\n"
       << "        -forcescaling                      (force rescaling even for low-res images)\n"
       << "        -minsampling <vox_dim>             (set minimum voxel dimension for sampling (in mm))\n"
    //       << "        -basescale <scale>                 (sets the scaling for the final scaling - default is 1.0)\n"
       << "        -applyxfm                          (applies transform (no optimisation) - requires -init)\n"
       << "        -applyisoxfm <scale>               (as applyxfm but forces isotropic resampling)\n"
       << "        -paddingsize <number of voxels>    (for applyxfm: interpolates outside image by size)\n"
       << "        -searchrx <min_angle> <max_angle>  (angles in degrees: default is -90 90)\n" 
       << "        -searchry <min_angle> <max_angle>  (angles in degrees: default is -90 90)\n" 
       << "        -searchrz <min_angle> <max_angle>  (angles in degrees: default is -90 90)\n" 
       << "        -nosearch                          (sets all angular search ranges to 0 0)\n" 
       << "        -coarsesearch <delta_angle>        (angle in degrees: default is 60)\n" 
       << "        -finesearch <delta_angle>          (angle in degrees: default is 18)\n" 
       << "        -schedule <schedule-file>          (replaces default schedule)\n"
       << "        -refweight <volume>                (use weights for reference volume)\n"
       << "        -inweight <volume>                 (use weights for input volume)\n"
       << "        -wmseg <volume>                    (white matter segmentation volume needed by BBR cost function)\n"
       << "        -wmcoords <text matrix>            (white matter boundary coordinates for BBR cost function)\n"
       << "        -wmnorms <text matrix>             (white matter boundary normals for BBR cost function)\n"
       << "        -fieldmap <volume>                 (fieldmap image in rads/s - must be already registered to the reference image)\n"
       << "        -fieldmapmask <volume>             (mask for fieldmap image)\n"
       << "        -pedir <index>                     (phase encode direction of EPI - 1/2/3=x/y/z & -1/-2/-3=-x/-y/-z)\n"
       << "        -echospacing <value>               (value of EPI echo spacing - units of seconds)\n"
       << "        -bbrtype <value>                   (type of bbr cost function: signed [default], global_abs, local_abs)\n"
       << "        -bbrslope <value>                  (value of bbr slope)\n"
       << "        -setbackground <value>             (use specified background value for points outside FOV)\n"
       << "        -noclamp                           (do not use intensity clamping)\n"
       << "        -noresampblur                      (do not use blurring on downsampling)\n"
       << "        -2D                                (use 2D rigid body mode - ignores dof)\n"
       << "        -verbose <num>                     (0 is least and default)\n"
       << "        -v                                 (same as -verbose 1)\n"
       << "        -i                                 (pauses at each stage: default is off)\n"
       << "        -version                           (prints version number)\n"
       << "        -help\n";
}


void globaloptions::print_version()
{
  cout << "FLIRT version " << version << endl;
}


