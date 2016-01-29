
/*  overlay.c - combine two images for colour overlay

    Stephen Smith, Christian  Beckmann and Matthew Webster, FMRIB Image Analysis Group

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

#include "libvis/miscpic.h"

using namespace NEWIMAGE;
using namespace MISCPIC;

void usage(void)
{
  printf("Usage: overlay <colour_type> <output_type> [-c] <background_image> <bg_min> <bg_max> <stat_image_1> <s1_min> <s1_max> [stat_image_2 s2min s2max] <output_image> [cbartype] [cbarfilename]\n");
  printf("colour_type: 0=solid 1=transparent colours\n");
  printf("output_type: 0=floating point (32 bit real) 1=integer (16 bit signed integer)\n");
  printf("-c : use checkerboard mask for overlay\n");
  printf("<bg_min> <bg_max> can be replaced by -a for automatic estimation of background display range or -A to use the full image range\n");
  printf("valid cbartypes colours are: ybg, valid cbartypes options are: s (stack) \n");
  exit(1);
}

//template <class T>
int fmrib_main(int argc, char* argv[], bool out_int)
{
  int colour_type, argindex=1, checker=0;
  float bgmin, bgmax, s1min, s1max, s2min, s2max;
  bool debug = false;

  volume<float> bg, s1, s2;
  string cbarfname = "";
  string cbartype = "";

  colour_type=atoi(argv[argindex++]);
  argindex++;

  if (!strcmp(argv[argindex],"-c")) {
    checker=1;
    argindex++;
  }

  if (!strcmp(argv[argindex],"-d")) {
    debug=true;
    argindex++;
  }

  read_volume(bg,string(argv[argindex++]));

  if (!strcmp(argv[argindex],"-a")) {
    bgmax = bg.percentile(0.98);
    bgmin = bg.percentile(0.02);
    argindex++;
  } else if (!strcmp(argv[argindex],"-A")) {
    bgmax = bg.max();
    bgmin = bg.min();
    argindex++;
  } else {
    bgmin=atof(argv[argindex++]);
    bgmax=atof(argv[argindex++]);
  }

  read_volume(s1,string(argv[argindex++]));
  s1min=atof(argv[argindex++]);
  s1max=atof(argv[argindex++]);
  
  if (argc-argindex-1>2){  
    read_volume(s2,string(argv[argindex++]));
    s2min=atof(argv[argindex++]);
    s2max=atof(argv[argindex++]);
  }
  else{
    s2 = s1; s2min= 0.0; s2max = 0.0;
  }

  if (argc-argindex-1==2){
    cbarfname = string(argv[argindex++]);
    cbartype = string(argv[argindex++]);
  }

  if(!argv[argindex]){
    cerr << "ERROR: Please specify an output filename " << endl << endl;
    exit(2);
  }
  else{
    //    miscpic<T> newpic;
    // volume<T> newvol;
    miscpic newpic;
    volume<float> newvol;

    newpic.overlay(newvol, bg, s1, s2, bgmin, bgmax, s1min, s1max, 
		   s2min, s2max, colour_type, checker, 
		   cbarfname, cbartype, out_int, debug);
    if(out_int)
      save_volume_dtype(newvol,string(argv[argindex]), DT_SIGNED_SHORT);
    else
      save_volume(newvol,string(argv[argindex]));

    return 0;
  } 
}


int main(int argc,char *argv[])
{

  if (argc<9)
    usage();

  int otype;
  otype = atoi(argv[2]);
  
  // if (otype == 0)
//     return call_fmrib_main(DT_FLOAT,argc,argv);
//   else if (otype == 1)
//     return call_fmrib_main(DT_SIGNED_SHORT,argc,argv);  
//   else
//     usage();

  if (otype == 0)
    return fmrib_main(argc,argv,false);
  else if (otype == 1)
    return fmrib_main(argc,argv,true);  
  else
    usage();

}

