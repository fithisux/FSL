/*  globaloptions.h

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

#ifndef __GLOBALOPTIONS_
#define __GLOBALOPTIONS_

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "newimage/costfns.h"


namespace NEWIMAGE {
  enum anglereps { Euler, Quaternion };
  enum interps { TriLinear, NearestNeighbour, Sinc, Spline };
  enum windowtype { Rect, Hanning, Blackman };
}

using namespace NEWIMAGE;

  typedef std::vector<RowVector> MatVec;
  typedef MatVec* MatVecPtr;


class globaloptions {
 public:
  static globaloptions& get();
  ~globaloptions() { delete gopt; }
  
  string version;
  
  std::vector<MatVec> usrmat;
  MatVec searchoptmat;
  MatVec preoptsearchmat;

  string inputfname;
  string outputfname;
  string reffname;
  string outputmatascii;
  string initmatfname;
  string refweightfname;
  string testweightfname;
  string wmsegfname;
  string wmcoordsfname;
  string wmnormsfname;
  string fmapfname;
  string fmapmaskfname;
  bool initmatsqform;
  bool printinit;
  Matrix initmat;

  string schedulefname;

  Costfn *impair;
  ColumnVector refparams;
  Matrix parammask;
  int no_params;
  int dof;
  bool usrsubset;
  int searchdof;
  int no_bins;
  costfns maincostfn;
  costfns searchcostfn;
  costfns currentcostfn;
  string optimisationtype;
  anglereps anglerep;
  float isoscale;
  float min_sampling;
  float lastsampling;
  float requestedscale;
  bool force_basescale;
  float basescale;
  bool force_scaling;
  float smoothsize;
  float fuzzyfrac;
  ColumnVector tolerance;
  ColumnVector boundguess;

  ColumnVector searchrx;
  ColumnVector searchry;
  ColumnVector searchrz;
  float coarsedelta;
  float finedelta;

  short datatype;
  bool forcedatatype;
  int verbose;
  bool debug;
  bool interactive;
  bool do_optimise;
  bool nosave;
  bool iso;
  bool resample;
  bool useweights;
  bool useseg;
  bool usecoords;
  bool mode2D;
  bool clamping;
  bool forcebackgnd;
  float backgndval;
  bool interpblur;
  interps interpmethod;
  float sincwidth;
  windowtype sincwindow;
  float paddingsize;
  int pe_dir;
  float echo_spacing;
  string bbr_type;
  float bbr_slope;

  int single_param;

  void parse_command_line(int argc, char** argv, const string &);

 private:
  globaloptions();
  
  const globaloptions& operator=(globaloptions&);
  globaloptions(globaloptions&);
      
  static globaloptions* gopt;

  void print_usage(int argc, char *argv[]);  
  void print_version();
  
};



//--------------------------------------------------------------------------//


inline globaloptions& globaloptions::get(){
  if(gopt == NULL)
    gopt = new globaloptions();
  
  return *gopt;
}

inline globaloptions::globaloptions()
{
  // set up defaults

  version = "";

  searchoptmat.clear();
  preoptsearchmat.clear();
  usrmat.resize(27);
  for (unsigned int i=0; i<usrmat.size(); i++) {
    usrmat[i].clear();
  }

  outputfname = "";
  reffname = "";

  inputfname = "";
  outputmatascii = "";
  initmatfname = "";
  refweightfname = "";
  testweightfname = "";
  wmsegfname = "";
  wmcoordsfname = "";
  wmnormsfname = "";
  fmapfname = "";
  fmapmaskfname = "";
  initmat = IdentityMatrix(4);
  initmatsqform = false;
  printinit = false;

  schedulefname = "";
  
  impair = 0;
  refparams.ReSize(12);
  refparams << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 << 1.0 << 1.0 << 1.0
	    << 0.0 << 0.0 << 0.0;
  parammask=IdentityMatrix(12);
  no_params = 12;
  dof = 12;
  usrsubset = false;
  searchdof = 12;
  no_bins = 256;
  maincostfn = CorrRatio;
  searchcostfn = CorrRatio;
  currentcostfn = CorrRatio;
  optimisationtype = "brent";
  anglerep = Euler;
  isoscale = 1.0;
  min_sampling = 1.0;
  lastsampling = 8;
  requestedscale = 1.0;
  force_basescale = false;
  basescale = 1.0;
  force_scaling = false;
  smoothsize = 1.0;
  fuzzyfrac = 0.5;
  tolerance.ReSize(12);
  tolerance << 0.005 << 0.005 << 0.005 << 0.2 << 0.2 << 0.2 << 0.002 
	    << 0.002 << 0.002 << 0.001 << 0.001 << 0.001; 
  boundguess.ReSize(2);
  boundguess << 10.0 << 1.0;

  searchrx.ReSize(2);
  searchry.ReSize(2);
  searchrz.ReSize(2);
  searchrx << -M_PI/2.0 << M_PI/2.0;
  searchry << -M_PI/2.0 << M_PI/2.0;
  searchrz << -M_PI/2.0 << M_PI/2.0;
  coarsedelta = 60.0*M_PI/180.0;
  finedelta = 18.0*M_PI/180.0;

  datatype = -1;
  forcedatatype = false;
  verbose = 0;
  debug = false;
  interactive = false;
  do_optimise = true;
  nosave = true;
  iso = true;
  resample = true;
  useweights = false;
  useseg = false;
  usecoords = false;
  mode2D = false;
  clamping = true;
  forcebackgnd = false;
  backgndval = 0.0;
  interpblur = true;
  interpmethod = TriLinear;
  sincwidth = 7.0; // voxels
  sincwindow = Hanning;
  paddingsize = 0.0;
  pe_dir=0;   // 1=x, 2=y, 3=z, -1=-x, -2=-y, -3=-z, 0=none
  echo_spacing = 5e-4;  // random guess (0.5ms) - units of seconds
  bbr_type = "signed";
  bbr_slope = -0.5;

  single_param = -1;
}

#endif

