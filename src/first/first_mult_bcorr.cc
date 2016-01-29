/*  first_mult_bcorr.cc
    Mark Jenkinson
    Copyright (C) 2009 University of Oxford  */

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

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include "math.h"

#include "utils/options.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"

using namespace std;
using namespace NEWIMAGE;
using namespace Utilities;
using namespace MISCMATHS;

string title="first_mult_bcorr (Version 1.0) University of Oxford (Mark Jenkinson)\nConverts ";
string examples="first_mult_bcorr [options] -i <T1_image> -c <4D_corrected_labels> -u <4D_uncorrected_labels> -o <output_image>";


Option<bool> verbose(string("-v,--verbose"), false, 
		     string("output F-stats to standard out"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<string> inname(string("-i,--in"), string(""),
		      string("filename of original T1 input image"),
		      true, requires_argument);
Option<string> ucorrname(string("-u,--uncorrected"), string(""),
		       string("filename of 4D image of uncorrected labels (with boundaries)"),
		       true, requires_argument);
Option<string> corrname(string("-c,--corrected"), string(""),
		       string("filename of 4D image of individually corrected labels"),
		       true, requires_argument);
Option<string> outname(string("-o,--out"), string(""),
		       string("output image name (3D label image)"),
		       true, requires_argument);

int nonoptarg;

////////////////////////////////////////////////////////////////////////////

// global variable for holding internal intensities of each structure

vector<volume<float> > internalints;

////////////////////////////////////////////////////////////////////////////


int set_up_globals(const volume<float>& t1im, const volume4D<int>& ucorrim)
{
  // store the set of internal intensities for each structure (for efficiency)
  volume<float> tmpvol;
  internalints.clear();
  for (int s=0; s<ucorrim.tsize(); s++) {
    int nvox=0;
    for (int z=ucorrim.minz(); z<=ucorrim.maxz(); z++) {
      for (int y=ucorrim.miny(); y<=ucorrim.maxy(); y++) {
	for (int x=ucorrim.minx(); x<=ucorrim.maxx(); x++) {
	  if ((ucorrim(x,y,z,s)>0) && (ucorrim(x,y,z,s)<100)) {
	    nvox++;
	  }
	}
      }
    }
    tmpvol.reinitialize(nvox,1,1);
    nvox=0;
    for (int z=ucorrim.minz(); z<=ucorrim.maxz(); z++) {
      for (int y=ucorrim.miny(); y<=ucorrim.maxy(); y++) {
	for (int x=ucorrim.minx(); x<=ucorrim.maxx(); x++) {
	  if ((ucorrim(x,y,z,s)>0) && (ucorrim(x,y,z,s)<100)) {
	    tmpvol(nvox++,0,0)=t1im(x,y,z);
	  }
	}
      }
    }
    internalints.push_back(tmpvol);
  }
  return 0;
}



int labelval(const volume<int>& lvol) 
{
  int val=200;
  for (int z=lvol.minz(); z<=lvol.maxz(); z++) {
    for (int y=lvol.miny(); y<=lvol.maxy(); y++) {
      for (int x=lvol.minx(); x<=lvol.maxx(); x++) {
	if ((lvol(x,y,z)>0) && (lvol(x,y,z)<val)) {
	  val=lvol(x,y,z);
	}
      }
    }
  }
  return val;
}

float scoreval(float fval, int ulab, int clab, int shapenum)
{
  // score based on reflected cdf around the median
  //  (i.e. closer to the median = better score)
  if (clab==0) { return 0.0; }
  int nint=internalints[shapenum].xsize(), nlower=0;
  for (int x=0; x<nint; x++) {
    if (internalints[shapenum](x,0,0)<fval) { nlower++; }
  }
  float score=0;
  if (nint>0) {
    score=((float) nlower)/((float) nint);
    if (score>0.5) { score=1.0-score; }
    if ((ulab>0) && (ulab<100)) { score += 1; }  // favour interior points
  }
  if (ulab<=0) { cerr << "Inconsistency in labelings passed to score" << endl; }
  return score;
}

int do_work() 
{
  volume<float> t1im;
  volume4D<int> ucorrim, corrim;
  read_volume(t1im,inname.value());
  read_volume4D(ucorrim,ucorrname.value());
  read_volume4D(corrim,corrname.value());
  if (!samesize(ucorrim,corrim)) {
    cerr << "ERROR: uncorrected and corrected images must be the same size" << endl;
    exit(EXIT_FAILURE);
  }
  if (!samesize(t1im,corrim[0])) {
    cerr << "ERROR: T1 image and corrected image must be the same size" << endl;
    exit(EXIT_FAILURE);
  }
  if (!samesize(t1im,ucorrim[0])) {
    cerr << "ERROR: T1 image and uncorrected image must be the same size"<<endl;
    exit(EXIT_FAILURE);
  }
  
  volume<int> outim(corrim[0]);
  outim*=0;

  int nsh=corrim.tsize(), nlab;

  set_up_globals(t1im,ucorrim);

  for (int z=t1im.minz(); z<=t1im.maxz(); z++) {
    for (int y=t1im.miny(); y<=t1im.maxy(); y++) {
      for (int x=t1im.minx(); x<=t1im.maxx(); x++) {
	// work out how many shapes label this voxel
	nlab=0;
	for (int t=0; t<nsh; t++) {
	  if (corrim(x,y,z,t)>0) {
	    nlab++;
	    outim(x,y,z)=corrim(x,y,z,t);  // useful if no conflict
	  }
	}
	if (nlab>1) {  // only do overlapping voxels (conflicts)
	  // calculate the cdf of this voxel's intensity wrt the internal
	  //  intensity dist for each shape
	  if (verbose.value()) { cout << "Found overlapping labels at voxel (" << x << "," << y << "," << z << ")" << endl; }
	  int bestlab=-1;
	  float bestscore=-1.0, score;
	  for (int lab=0; lab<nsh; lab++) {
	    score = scoreval(t1im(x,y,z),ucorrim(x,y,z,lab),corrim(x,y,z,lab),lab);
	    if (score>bestscore) { bestscore=score; bestlab=lab; }
	    //if (verbose.value()) { cout << "  score for label " << lab << " is " << score << endl; }
	  }
	  if (bestlab>=0) {
	    outim(x,y,z)=labelval(ucorrim[bestlab]);
	  } // else it failed, so do nothing
	}
      }
    }
  }
  outim.setAuxFile("MGH-Subcortical");  // set nice colormap
  save_volume(outim,outname.value());
  return 0;
}

////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{
	
  Tracer tr("main");
  OptionParser options(title, examples);
	
  try {
    // must include all wanted options here (the order determines how
    //  the help message is printed)
    options.add(inname);
    options.add(outname);
    options.add(ucorrname);
    options.add(corrname);
    options.add(verbose);
    options.add(help);
    nonoptarg = options.parse_command_line(argc, argv);
		
    // line below stops the program if the help was requested or 
    //  a compulsory option was not set
    if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
      {
	options.usage();
	exit(EXIT_FAILURE);
      }
		
    // Call the local functions
    do_work();
		
  } catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } catch(...) {
    cerr << "Aborted" << endl;
  } 
	
  return 0;// do_work(argc,argv);
}

