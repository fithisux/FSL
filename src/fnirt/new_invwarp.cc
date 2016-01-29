
/*  new_invwarp.cc

    Jesper Andersson, FMRIB Image Analysis Group

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

#include "utils/options.h"
#include "miscmaths/miscmaths.h"
#include "warpfns/warpfns.h"
#include "warpfns/fnirt_file_reader.h"
#include "warpfns/fnirt_file_writer.h"
#include "tetrahedron.h"
#include "dilator.h"

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

using namespace std;
using namespace Utilities;
using namespace NEWMAT;
using namespace MISCMATHS;
using namespace NEWIMAGE;

// COMMAND LINE OPTIONS

string title="new_invwarp\nCopyright(c) 2010, University of Oxford (Jesper Andersson)";
string examples="invwarp -w warpvol -o invwarpvol -r refvol";

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<bool> debug(string("--debug"), false,
		  string("turn on debugging output"),
		  false, no_argument);
Option<bool> abswarp(string("--abs"), false,
		  string("use absolute warp convention (default): x' = w(x)"),
		  false, no_argument);
Option<bool> relwarp(string("--rel"), false,
		  string("use relative warp convention: x' = x + w(x)"),
		  false, no_argument);
Option<bool> nojaccon(string("--noconstraint"), false,
		  string("do not apply the Jacobian constraint"),
		  false, no_argument);
Option<float> jmin(string("--jmin"), 0.01,
			string("minimum acceptable Jacobian value for constraint (default 0.01)"),
			false, requires_argument);
Option<float> jmax(string("--jmax"), 100.0,
			string("maximum acceptable Jacobian value for constraint (default 100.0)"),
			false, requires_argument);
Option<string> refvolname(string("-r,--ref"), string(""),
		       string("filename for new reference image, i.e., what was originally the input image (determines inverse warpvol's FOV and pixdims)"),
		       true, requires_argument);
Option<string> outname(string("-o,--out"), string(""),
		       string("filename for output (inverse warped) image"),
		       true, requires_argument);
Option<string> warpname(string("-w,--warp"), string(""),
			string("filename for warp/shiftmap transform (volume)"),
			true, requires_argument);
Option<string> extrapolation(string("-e,--extrap"), string("affine"),
			     string("Extrapolation method to use outside FOV [affine | fancy], default affine"),
			     true, requires_argument);

int new_invwarp()
{
  // Read warps
  // N.B. that latest version of FnirtFileReader
  // will separate Affine and non-linear component
  // even for files that aren't fnirt coefficient
  // files. This makes it easier to write code where
  // one truly does not need to care about where the
  // files come from,
  if (verbose.value()) cout << "Reading input files" << endl;
  AbsOrRelWarps      spec_wt=UnknownWarps;        // Specified warp convention
  if (abswarp.value()) spec_wt = AbsoluteWarps;
  else if (relwarp.value()) spec_wt = RelativeWarps;
  FnirtFileReader    fnirtfile(warpname.value(),spec_wt,verbose.set());
  boost::shared_ptr<volume4D<float> > warpvol(new volume4D<float>);
  *warpvol = fnirtfile.FieldAsNewimageVolume4D();
  NEWMAT::Matrix aff = fnirtfile.AffineMat();

  // Read reference volume and set size of invwarp
  volume4D<float>    invwarp;
  read_volume4D(invwarp,refvolname.value());
  while (invwarp.tsize()<3) { invwarp.addvolume(invwarp[0]); }
  while (invwarp.tsize()>3) { invwarp.deletevolume(invwarp.maxt()); }
  // inwarp.tsize() == 3 here
  invwarp=0.0f;
  invwarp.setDisplayMaximumMinimum(0,0);
  const volume<float>& ref=invwarp[0];   // To make subsequent code easier to understand (I hope)

  // Create mm->mm transformation matrix for
  // use within tetrahedron and for use as
  // initial guess.
  NEWMAT::Matrix B = warpvol->sampling_mat().i()*aff*ref.sampling_mat();      // For initial guess
  double b11, b12, b13, b14;
  double b21, b22, b23, b24;
  double b31, b32, b33, b34;
  b11=B(1,1); b12=B(1,2); b13=B(1,3); b14=B(1,4);
  b21=B(2,1); b22=B(2,2); b23=B(2,3); b24=B(2,4);
  b31=B(3,1); b32=B(3,2); b33=B(3,3); b34=B(3,4);

  // Rescale displacement fields mm->voxels_in_ref_space
  (*warpvol)[0] /= ref.xdim();
  (*warpvol)[1] /= ref.ydim();
  (*warpvol)[2] /= ref.zdim();

  // Do the inversion
  if (verbose.value()) {
    cout << "Inverting field" << endl;
  }
  NEWMAT::Matrix iB = B.i();
  Tetrahedron  tet(0,0,0,warpvol,iB); // The guy doing the work
  tet.SetIgnoreFOV();
  bool found_it = false;
  for (int k=0; k<invwarp.zsize(); k++) {
    if (verbose.value()) { cout << "\b\b\b\b\b\b\b\b\b\b" << "Slice: " << k+1; cout.flush(); }
    for (int j=0; j<invwarp.ysize(); j++) {
      for (int i=0; i<invwarp.xsize(); i++) {
        double ox, oy, oz;
        double ax = b11*i+b12*j+b13*k+b14;
        double ay = b21*i+b22*j+b23*k+b24;
        double az = b31*i+b32*j+b33*k+b34;
        if (!i || !found_it) {
	  tet.SetFirstPoint(int(ax),int(ay),int(az));
	}
        if (tet.FindPoint(double(i),double(j),double(k),ox,oy,oz)) {
          invwarp(i,j,k,0) = ox-ax; 
          invwarp(i,j,k,1) = oy-ay;
          invwarp(i,j,k,2) = oz-az;
          found_it = true;
	}
        else { // Indicates singularity.
          invwarp(i,j,k,0) = -999; // NaN
          invwarp(i,j,k,1) = -999; // NaN
          invwarp(i,j,k,2) = -999; // NaN
          found_it = false;
        }
      }
    }
  }
  if (verbose.value()) cout << endl;

  // Perform dilations to replace NaNs
  // with average of non-NaN neighbours.
  if (verbose.value()) cout << "Fudging values at edge of FOV" << endl;
  for (int i=0; i<3; i++) {
    Dilator dil(invwarp[i]);
    unsigned int n = dil.Dilate();
    while (n) {
      n = dil.Dilate();
    }
    invwarp[i] = dil.Get();
  }

  // Convert back to mm
  std::vector<double> vxs = fnirtfile.VoxelSize();
  invwarp[0] *= vxs[0]; invwarp[1] *= vxs[1]; invwarp[2] *= vxs[2];

  // Add affine component back in
  NEWMAT::Matrix M = ref.sampling_mat();
  double m11, m12, m13, m14;
  double m21, m22, m23, m24;
  double m31, m32, m33, m34;
  m11=M(1,1); m12=M(1,2); m13=M(1,3); m14=M(1,4);
  m21=M(2,1); m22=M(2,2); m23=M(2,3); m24=M(2,4);
  m31=M(3,1); m32=M(3,2); m33=M(3,3); m34=M(3,4);
  NEWMAT::Matrix A = aff;
  double a11, a12, a13, a14;
  double a21, a22, a23, a24;
  double a31, a32, a33, a34;
  a11=A(1,1); a12=A(1,2); a13=A(1,3); a14=A(1,4);
  a21=A(2,1); a22=A(2,2); a23=A(2,3); a24=A(2,4);
  a31=A(3,1); a32=A(3,2); a33=A(3,3); a34=A(3,4);
  for (int k=0; k<invwarp.zsize(); k++) {
    for (int j=0; j<invwarp.ysize(); j++) {
      for (int i=0; i<invwarp.xsize(); i++) {
        double mmx = m11*i+m12*j+m13*k+m14;
        double mmy = m21*i+m22*j+m23*k+m24;
        double mmz = m31*i+m32*j+m33*k+m34;
        double ax = a11*mmx+a12*mmy+a13*mmz+a14;
        double ay = a21*mmx+a22*mmy+a23*mmz+a24;
        double az = a31*mmx+a32*mmy+a33*mmz+a34;
        invwarp(i,j,k,0) += ax-mmx; 
        invwarp(i,j,k,1) += ay-mmy;
        invwarp(i,j,k,2) += az-mmz;
      }
    }
  }

  // Constrain range of Jacobians of inverse field
  if (!nojaccon.value()) {
    if (verbose.value()) cout << "Constraining range of Jacobians in inverse field" << endl;
    convertwarp_rel2abs(invwarp);
    constrain_topology(invwarp,jmin.value(),jmax.value());
    convertwarp_abs2rel(invwarp);
  }
  
  // Save inverse. I save it as a "fnirt field", which is different
  // from "old" invwarp. The advantage is that now the field is 
  // "in the system" and we can make assumptions about it.
  if (verbose.value()) cout << "Saving inverted field" << endl;
  FnirtFileWriter(outname.value(),ref,invwarp);
 
  return(EXIT_SUCCESS); 
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// The main pretty much just does some argument checking and then
// calls the new_invwarp function that does the job.
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

int main(int argc, char *argv[])
{

  Tracer tr("main");

  OptionParser options(title, examples);

  try {
    options.add(warpname);
    options.add(outname);
    options.add(refvolname);
    options.add(relwarp);
    options.add(abswarp);
    options.add(nojaccon);
    options.add(jmin);
    options.add(jmax);
    options.add(debug);
    options.add(verbose);
    options.add(help);
    
    int nparsed = options.parse_command_line(argc, argv);
    if (nparsed < argc) {
      for (; nparsed<argc; nparsed++) {
        cerr << "Unknown input: " << argv[nparsed] << endl;
      }
      exit(EXIT_FAILURE);
    }

    if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
      {
	options.usage();
	exit(EXIT_FAILURE);
      }
  } catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } 

  if (abswarp.value() && relwarp.value()) {
    cerr << "--abs and --rel flags cannot both be set" << endl;
    exit(EXIT_FAILURE);
  }

  volume4D<float>   warpvol;
  try {
    read_volume4D_hdr_only(warpvol,warpname.value());
  }
  catch(...) {
    cerr << "invwarp: Problem reading warp-file " << warpname.value() << endl;
    exit(EXIT_FAILURE);
  }
  if ((warpvol.intent_code()==FSL_CUBIC_SPLINE_COEFFICIENTS || 
       warpvol.intent_code()==FSL_QUADRATIC_SPLINE_COEFFICIENTS || 
       warpvol.intent_code()==FSL_DCT_COEFFICIENTS ||
       warpvol.intent_code()==FSL_FNIRT_DISPLACEMENT_FIELD) &&
      (abswarp.value() || relwarp.value())) {
    cout << "--abs and --rel flags ignored when reading fnirt coefficient files" << endl;
  }

  return new_invwarp();
}

