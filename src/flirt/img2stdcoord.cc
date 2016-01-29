/*  img2stdcoord.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2000 University of Oxford  */

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

#include <string>
#include <iostream>
#include <fstream>
#include <unistd.h>
#define WANT_STREAM
#define WANT_MATH

#ifndef EXPOSE_TREACHEROUS
#define EXPOSE_TREACHEROUS
#endif

#include "newmatap.h"
#include "newmatio.h"
#include "miscmaths/miscmaths.h"
#include "newimage/newimageall.h"
#include "warpfns/warpfns.h"
#include "warpfns/fnirt_file_reader.h"

#ifndef NO_NAMESPACE
 using namespace MISCMATHS;
 using namespace NEWMAT;
 using namespace NEWIMAGE;
#endif

////////////////////////////////////////////////////////////////////////////
// the real defaults are provided in the function parse_command_line

class globaloptions {
public:
  string stdfname;
  string imgfname;
  string prexfmfname;
  string coordfname;
  string warpfname;
  bool usestd;
  bool mm;
  int verbose;
public:
  globaloptions();
  ~globaloptions() {};
};

globaloptions globalopts;


globaloptions::globaloptions()
{
  // set up defaults
  stdfname = "";
  imgfname = "";
  coordfname = "";
  prexfmfname = "";
  warpfname = "";
  verbose = 0;
  usestd = false;
  mm = false;
}


////////////////////////////////////////////////////////////////////////////

// Parsing functions for command line parameters

void print_usage(int argc, char *argv[])
{
  cout << "Usage: " << argv[0] << " [options] <filename containing coordinates>\n"
       << "e.g.   " << argv[0] << " -img <invol> -std <standard image> -xfm <img2standard mat file> <coordinate file>\n"
       << "       " << argv[0] << " -img <invol> <coordinate file>\n"
       << "       " << argv[0] << " -img <invol> - \n\n"
       << "  Options are:\n"
       << "        -img <filename of input image>\n"
       << "        -std <filename of standard image>\n"
       << "        -xfm <filename of affine transform   (e.g. example_func2standard.mat)>\n"
       << "        -warp <filename of warpfield (e.g. highres2standard_warp.nii.gz)>\n"
       << "        -premat <filename of pre-warp affine transform  (e.g. example_func2highres.mat)>   (default=identity)\n"
       << "        -vox                                 (input coordinates in voxels - default)\n"
       << "        -mm                                  (input coordinates in mm)\n"
       << "        -v                                   (verbose output)\n"
       << "        -verbose                             (more verbose output)\n"
       << "        -help\n\n"
       << " Notes:\n"
       << "  (1) if '-' is used as coordinate filename then coordinates are read from standard input\n"
       << "  (2) the -img option is compulsory\n";
}


void parse_command_line(int argc, char* argv[])
{
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
      globalopts.coordfname = arg;
      n++;
      continue;
    }
    
    // put options without arguments here
    if ( arg == "-help" ) {
      print_usage(argc,argv);
      exit(0);
    } else if ( arg == "-vox" ) {
      globalopts.mm = false;
      n++;
      continue;
    } else if ( arg == "-mm" ) {
      globalopts.mm = true;
      n++;
      continue;
    } else if ( arg == "-flirt" ) {
      cerr << "WARNING::Using outdated options, please update to new usage" << endl;
      // do nothing anymore
      n++;
      continue;
    } else if ( arg == "-v" ) {
      globalopts.verbose = 1;
      n++;
      continue;
    } else if ( arg == "-verbose" ) {
      globalopts.verbose = 5;
      n++;
      continue;
    } else if ( arg == "-") {
      globalopts.coordfname = "-";
      n++;
      continue;
    }

    if (n+1>=argc) 
      { 
	cerr << "Lacking argument to option " << arg << endl;
	break; 
      }

    // put options with 1 argument here
    if ( arg == "-std") {
      globalopts.stdfname = argv[n+1];
      globalopts.usestd = true;
      n+=2;
      continue;
    } else if ( arg == "-img") {
      globalopts.imgfname = argv[n+1];
      n+=2;
      continue;
    } else if ( ( arg == "-xfm") || ( arg == "-premat") ) {
      globalopts.prexfmfname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-warp") {
      globalopts.warpfname = argv[n+1];
      n+=2;
      continue;
    } else { 
      cerr << "Unrecognised option " << arg << endl;
      exit(-1);
    } 

  }  // while (n<argc)

  if ((globalopts.imgfname.size()<1)) {
    cerr << "ERROR:: input image filename not found\n\n";
  }
  if ((globalopts.usestd) && (globalopts.stdfname.size()<1)) {
    cerr << "ERROR:: standard image filename not found\n\n";
  }
}

////////////////////////////////////////////////////////////////////////////

void print_info(const volume<float>& vol, const string& name) {
  cout << name << ":: SIZE = " << vol.xsize() << " x " << vol.ysize() 
       << " x " << vol.zsize() << endl;
  cout << name << ":: DIMS = " << vol.xdim() << " x " << vol.ydim() 
       << " x " << vol.zdim() << " mm" << endl << endl;
}  

////////////////////////////////////////////////////////////////////////////

ColumnVector NewimageCoord2NewimageCoord(const FnirtFileReader& fnirtfile, const Matrix& affmat,
					 const volume<float>& srcvol, const volume<float>& destvol,
					 const ColumnVector& srccoord)
{
  ColumnVector retvec;
  if (fnirtfile.IsValid()) {
    // in the following affmat=example_func2highres.mat, fnirtfile=highres2standard_warp.nii.gz
    static volume4D<float> fieldVolume( fnirtfile.FieldAsNewimageVolume4D(true) );
    retvec = NewimageCoord2NewimageCoord(affmat,fieldVolume,true,srcvol,destvol,srccoord);
  } else {
    retvec = NewimageCoord2NewimageCoord(affmat,srcvol,destvol,srccoord);
  }
  return retvec;
}

////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{
  parse_command_line(argc,argv);


  volume<float> imgvol, stdvol;
    // read volumes
  if (read_volume_hdr_only(imgvol,globalopts.imgfname)<0) {
    cerr << "Cannot read input image" << endl;
    return -1;
  }
  if (globalopts.usestd && (read_volume_hdr_only(stdvol,globalopts.stdfname)<0)) {
    cerr << "Cannot read standard image" << endl;
    return -1;
  }
    
  if (globalopts.verbose>3) {
    if (globalopts.usestd) { print_info(stdvol,"standard image"); }
    print_info(imgvol,"input image");
  }

  // read matrices
  Matrix affmat(4,4);
  bool use_sform=false;
  if (globalopts.prexfmfname.length()>0) {
    affmat = read_ascii_matrix(globalopts.prexfmfname);
    use_sform = false;
    if (affmat.Nrows()<4) {
      cerr << "Cannot read transform file" << endl;
      return -2;
    }
  } else {
    use_sform = true;
    affmat = IdentityMatrix(4);
  }

    
  if (globalopts.verbose>3) {
    cout << " affmat =" << endl << affmat << endl << endl;
  }

  // Read in warps from file (if specified)
  FnirtFileReader  fnirtfile;
  AbsOrRelWarps    wt = UnknownWarps;
  if (globalopts.warpfname != "") {
    try {
      fnirtfile.Read(globalopts.warpfname,wt,globalopts.verbose>3);
    }
    catch (...) {
      cerr << "An error occured while reading file: " << globalopts.warpfname << endl;
      exit(EXIT_FAILURE);
    }
  }

  /////////////// SET UP MATRICES ////////////////

  if ( (stdvol.sform_code()==NIFTI_XFORM_UNKNOWN) && 
       (stdvol.qform_code()==NIFTI_XFORM_UNKNOWN) ) { 
    cerr << "WARNING:: mm coordinate space not set in standard volume" << endl; 
  }
  if (globalopts.verbose>3) {
    cout << " stdvox2world =" << endl << stdvol.newimagevox2mm_mat() << endl << endl;
  }

  // initialise coordinate vectors

  ColumnVector imgcoord(4), stdcoord(4), oldimg(4);
  imgcoord = 0;
  stdcoord = 0;
  imgcoord(4)=1;
  stdcoord(4)=1;
  oldimg = 0;  // 4th component set to 0, so that initially oldimg -ne imgcoord


  bool use_stdin = false;
  if ( (globalopts.coordfname=="-") || (globalopts.coordfname.size()<1)) {
    use_stdin = true;
  }

  if (globalopts.verbose>0) {
    cout << "Coordinates in standard image (in mm):" << endl;
  }

  // set up coordinate reading (from file or stdin) //
  ifstream matfile(globalopts.coordfname.c_str());

  if (use_stdin) {
    if (globalopts.verbose>0) {
      cout << "Please type in input image coordinates";
      if (globalopts.mm) { 
	cout << " (in mm) :" << endl;
      } else { 
	cout << " (in voxels) :" << endl; 
      }
    } 
  } else {
    if (!matfile) { 
      cerr << "Could not open matrix file " << globalopts.coordfname << endl;
      return -1;
    }
  }
  
  // loop around reading coordinates and displaying output
  
  while ( (use_stdin && (!cin.eof())) || ((!use_stdin) && (matfile >> imgcoord(1) >> imgcoord(2) >> imgcoord(3))) ) {
    if  (use_stdin) {
      cin >> imgcoord(1) >> imgcoord(2) >> imgcoord(3);
      // this is in case the pipe continues to input a stream of zeros
      if (oldimg == imgcoord)  return 0;
      oldimg = imgcoord;
    }
     
    if (globalopts.mm) {  // in mm
      stdcoord = stdvol.newimagevox2mm_mat() * 
	NewimageCoord2NewimageCoord(fnirtfile,affmat,imgvol,stdvol,imgvol.newimagevox2mm_mat().i() * imgcoord);
    } else { // in voxels
      stdcoord = stdvol.newimagevox2mm_mat() * 
	NewimageCoord2NewimageCoord(fnirtfile,affmat,imgvol,stdvol,imgvol.niftivox2newimagevox_mat() * imgcoord); 
    }
    cout << stdcoord(1) << "  " << stdcoord(2) << "  " << stdcoord(3) << endl;
  }
  
  if (!use_stdin) { matfile.close(); }
  
  return 0;
}


