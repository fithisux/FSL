/*  fnirtfileutils.cpp

    Jesper Andersson, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

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

#ifndef EXPOSE_TREACHEROUS
#define EXPOSE_TREACHEROUS           // To allow us to use .sampling_mat()
#endif

#include <boost/shared_ptr.hpp>
#include "newmat.h"
#include "utils/options.h"
#include "miscmaths/miscmaths.h"
#include "newimage/newimageall.h"
#include "basisfield/basisfield.h"
#include "basisfield/splinefield.h"
#include "basisfield/dctfield.h"
#include "warpfns/warpfns.h"
#include "warpfns/fnirt_file_reader.h"
#include "warpfns/fnirt_file_writer.h"

using namespace boost;
using namespace Utilities;
using namespace NEWMAT;
using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace BASISFIELD;

// COMMAND LINE OPTIONS

string title="fnirtfileutils\nCopyright(c) 2012, University of Oxford (Jesper Andersson)";
string examples=string("fnirtfileutils --in=my_fnirtcoefs --ref=my_refvol --out=my_field\n") +
                string("fnirtfileutils --in=my_fnirtcoefs --ref=my_refvol --jac=my_jacobians\n") + 
                string("fnirtfileutils --in=my_fnirtcoefs --ref=my_refvol --out=my_field --withaff\n");

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<string> inname(string("-i,--in"), string(""),
		      string("\tfilename of input coefficient volume (to be converted)"),
		      true, requires_argument);
Option<string> refname(string("-r,--ref"), string(""),
		       string("filename for reference volume"),
		       false, requires_argument);
Option<string> outname(string("-o,--out"), string(""),
		       string("filename for output (field/coef) volume - uses relative warp convention"),
		       false, requires_argument);
Option<string> outformat(string("-f,--outformat"),string("field"),
			 string("Output format [field spline], default=field"),
			 false,requires_argument);

/* The DCT option hidden until bending energy has been implemented for DCT

Option<string> outformat(string("-f,--outformat"),string("field"),
			 string("Output format [field spline DCT], default=field"),
			 false,requires_argument);
*/

vector<float> warpresdefault(3,10.0);
Option<vector<float> > warpres(string("-w,--warpres"),warpresdefault,
			       string("Warp resolution (mm), only relevant when --outformat=spline"),
			       false,requires_argument);
vector<int> knotspacedefault;
Option<vector<int> > knotspace(string("-k,--knotspace"),knotspacedefault,
			         string("Knot-spacing (voxels), only relevant when --outformat=spline"),
			         false,requires_argument);
vector<int> dctorderdefault;
HiddenOption<vector<int> > dctorder(string("-d,--dctorder"),dctorderdefault,
			            string("DCT-order, only relevant when --outformat=DCT"),
			            false,requires_argument);
Option<string> jacname(string("-j,--jac"), string(""),
		       string("filename for output (jacobian map) volume"),
		       false, requires_argument);
Option<bool> withaff(string("-a,--withaff"), false,
		     string("If set, the affine transform is included in the field/jacobian"),
		     false, no_argument);


int main(int argc, char *argv[])
{
  OptionParser options(title, examples);

  //
  // Parse command line
  //
  try {
    options.add(inname);
    options.add(refname);
    options.add(outname);
    options.add(outformat);
    options.add(warpres);
    options.add(dctorder);
    options.add(knotspace);
    options.add(jacname);
    options.add(withaff);
    options.add(verbose);
    options.add(help);

    int i=options.parse_command_line(argc, argv);
    if (i < argc) {
      for (; i<argc; i++) {
        cerr << "Unknown input: " << argv[i] << endl;
      }
      exit(EXIT_FAILURE);
    }
    if (help.value()) {
      options.usage();
      exit(EXIT_SUCCESS);
    }
    else if (!options.check_compulsory_arguments(true)) {
      options.usage();
      exit(EXIT_FAILURE);
    }
  }  
  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } 
  catch(std::exception &e) {
    cerr << e.what() << endl;
    exit(EXIT_FAILURE);
  }
  //
  // Do some post-parsing not achieved by the options class
  //
  if (!outname.set() && !jacname.set()) {
    cerr << "fnirtfileutils: At least one of --out and --jac must be set" << endl;
    exit(EXIT_FAILURE);
  }
  if ((warpres.set() && dctorder.set()) ||
      (warpres.set() && knotspace.set()) ||
      (knotspace.set() && dctorder.set())) {
    cerr << "fnirtfileutils: The --warpres, --dctorder and --knotspace parameters are mutually exclusive" << endl;
    exit(EXIT_FAILURE);
  }
  if (outformat.set() && (outformat.value()!=string("field") && outformat.value()!=string("spline") && outformat.value()!=string("DCT"))) {
    cerr << "fnirtfileutils: Invalid --outformat option " << outformat.value() << endl;
    exit(EXIT_FAILURE);
  }

  //
  // Open coefficient-file and reference volume
  //
  FnirtFileReader   infile;
  volume4D<float>   coefvol;
  try {
    infile.Read(inname.value(),UnknownWarps,verbose.value());
  }
  catch(...) {
    cerr << "fnirtfileutils: Problem reading coef-file " << inname.value() << endl;
    exit(EXIT_FAILURE);
  }
  if (!infile.IsValid()) {
    cerr << "fnirtfileutils: " << inname.value() << " not a valid fnirt-coefficient file" << endl; 
    exit(EXIT_FAILURE);
  }
  volume<float>     refvol;
  if (refname.set()) { // If a reference volume was explicitly specified
    try {
      read_volume(refvol,refname.value());
    }
    catch(...) {
      cerr << "fnirtfileutils: Problem reading reference-file " << refname.value() << endl;
      exit(EXIT_FAILURE);
    }
  }
  else {
    std::vector<unsigned int> imsz = infile.FieldSize();
    std::vector<double> vxs = infile.VoxelSize();
    refvol.reinitialize(imsz[0],imsz[1],imsz[2]);
    refvol.setdims(vxs[0],vxs[1],vxs[2]);
  }

  if (outname.set()) {   // If we are to return some representation of a displacement-field
    if (outformat.value() == string("field")) { // If we are to return it in field representation
      volume4D<float> field;
      try {
        field = infile.FieldAsNewimageVolume4D(withaff.value());
      }
      catch(...) {
        cerr << "fnirtfileutils: Problem creating displacement-field " << outname.value() << endl;
        exit(EXIT_FAILURE);
      }
      try {
        FnirtFileWriter  outfile(outname.value(),refvol,field);
      }
      catch(...) {
        cerr << "fnirtfileutils: Problem writing displacement-field " << outname.value() << endl;
        exit(EXIT_FAILURE);
      }
    }
    else if (outformat.value() == string("spline")) { // If we are to return it in spline representation
      vector<unsigned int> ksp(3,0);
      if (knotspace.set()) {
        if (knotspace.value().size() == 3) for (int i=0; i<3; i++) ksp[i] = static_cast<unsigned int>(knotspace.value()[i]);
	else if (knotspace.value().size() == 1) for (int i=0; i<3; i++) ksp[i] = static_cast<unsigned int>(knotspace.value()[0]);
        else {
          cerr << "fnirtfileutils: --knotspace must be 1 or 3 elements long " << endl;
          exit(EXIT_FAILURE);
	}
      }
      else {
        vector<float> tmpwres(3,0.0);
        if (warpres.value().size() == 3) tmpwres=warpres.value();
        else if (warpres.value().size() == 1) for (int i=0; i<3; i++) tmpwres[i] = (warpres.value())[0];
        else {
          cerr << "fnirtfileutils: --warpres must be 1 or 3 elements long " << endl;
          exit(EXIT_FAILURE);
	}
        ksp[0] = static_cast<unsigned int>(floor(tmpwres[0]/refvol.xdim() + 0.5));
        ksp[1] = static_cast<unsigned int>(floor(tmpwres[1]/refvol.ydim() + 0.5));
        ksp[2] = static_cast<unsigned int>(floor(tmpwres[2]/refvol.zdim() + 0.5));
      }
      if (verbose.value()) {
        cout << "Calculating spline coefficients for ksp = " << ksp[0] << ", " << ksp[1] << ", " << ksp[2] << endl;
      }
      vector<boost::shared_ptr<basisfield> > fields(3);
      try {
        for (int i=0; i<3; i++) fields[i] = boost::shared_ptr<splinefield>(new splinefield(infile.FieldAsSplinefield(i,ksp)));
      }
      catch(...) {
        cerr << "fnirtfileutils: Problem calculating spline coefficients " << endl;
        exit(EXIT_FAILURE);
      }
      try {
        FnirtFileWriter outfile(outname.value(),fields,infile.AffineMat());
      }
      catch (...) {
        cerr << "fnirtfileutils: Problem writing spline coefficient file " << outname.value() << endl;
        exit(EXIT_FAILURE);
      }
    }
    else if (outformat.value() == string("DCT")) { // If we are to return it in DCT representation
      vector<unsigned int>  order(3);
      if (dctorder.set()) {
        if (dctorder.value().size() == 3) for (int i=0; i<3; i++) order[i] = static_cast<unsigned int>(dctorder.value()[i]);
	else if (dctorder.value().size() == 1) for (int i=0; i<3; i++) order[i] = static_cast<unsigned int>(dctorder.value()[0]);
        else {
          cerr << "fnirtfileutils: --dctorder must be 1 or 3 elements long " << endl;
          exit(EXIT_FAILURE);
	}
      }
      else {
        vector<float>  tmpwres(3);
        if (warpres.value().size() == 3) tmpwres=warpres.value();
        else if (warpres.value().size() == 1) for (int i=0; i<3; i++) tmpwres[i] = (warpres.value())[0];
        else {
          cerr << "fnirtfileutils: --warpres must be 1 or 3 elements long " << endl;
          exit(EXIT_FAILURE);
	}
        order[0] = static_cast<unsigned int>(floor(refvol.xdim()*refvol.xsize()/tmpwres[0] + 0.5));
        order[1] = static_cast<unsigned int>(floor(refvol.ydim()*refvol.ysize()/tmpwres[1] + 0.5));
        order[2] = static_cast<unsigned int>(floor(refvol.zdim()*refvol.zsize()/tmpwres[2] + 0.5));
      }
      if (verbose.value()) {
        cout << "Calculating DCT coefficients for order = " << order[0] << ", " << order[1] << ", " << order[2] << endl;
      }
      vector<boost::shared_ptr<basisfield> > fields(3);
      try {
        for (int i=0; i<3; i++) fields[i] = boost::shared_ptr<dctfield>(new dctfield(infile.FieldAsDctfield(i,order)));
      }
      catch(...) {
        cerr << "fnirtfileutils: Problem calculating DCT coefficients " << endl;
        exit(EXIT_FAILURE);
      }
      try {
        FnirtFileWriter outfile(outname.value(),fields,infile.AffineMat());
      }
      catch (...) {
        cerr << "fnirtfileutils: Problem writing DCT coefficient file " << outname.value() << endl;
        exit(EXIT_FAILURE);
      }
    }
  }  

  if (jacname.set()) { // If we are to return a Jacobian field
    if (infile.Type()!=FnirtSplineDispType && infile.Type()!=FnirtDCTDispType) {
      cerr << "fnirtfileutils: Can presently only generate Jacobian from coefficient files " << endl;
      exit(EXIT_FAILURE);
    }
    else {
      try {
        volume<float>  ojac(refvol.xsize(),refvol.ysize(),refvol.zsize());
        volume<float>  jac = infile.Jacobian(withaff.value());
        if (samesize(ojac,jac)) {
	  ojac = jac;
          ojac.copyproperties(refvol);
	}
        else {
          ojac.copyproperties(refvol);
          jac.setextrapolationmethod(extraslice);
          Matrix O2I = jac.sampling_mat().i() * ojac.sampling_mat();
          ColumnVector xo(4), xi(4);
          int zs = ojac.zsize(), ys = ojac.ysize(), xs = ojac.xsize();
          xo(4) = 1.0;
          for (int z=0; z<zs; z++) { xo(3)=double(z); for (int y=0; y<ys; y++) { xo(2)=double(y); for (int x=0; x<xs; x++) {
		xo(1)=double(x);
                xi = O2I*xo;
                ojac(x,y,z) = jac.interpolate(xi(1),xi(2),xi(3));
	  } } }
	}
        ojac.setDisplayMaximum(0.0); ojac.setDisplayMinimum(0.0);
        save_volume(ojac,jacname.value());
      }
      catch(...) {
        cerr << "fnirtfileutils: Problem creating/writing Jacobian file " << jacname.value() << endl;
        exit(EXIT_FAILURE);
      }
    }
  }
  exit(EXIT_SUCCESS);
}
