/*  makerot.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2003 University of Oxford  */

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

// Skeleton application framework for using newimage

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "utils/options.h"
#include <vector>

using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;
using namespace std;

// The two strings below specify the title and example usage that is
//  printed out as the help or usage message

string title="makerot \nCopyright(c) 2003, University of Oxford (Mark Jenkinson)";
string examples="makerot [options] --theta=angle";

// Each (global) object below specificies as option and can be accessed
//  anywhere in this file (since they are global).  The order of the
//  arguments needed is: name(s) of option, default value, help message,
//       whether it is compulsory, whether it requires arguments
// Note that they must also be included in the main() function or they
//  will not be active.

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<float> theta(string("-t,--theta"), 0.0,
		  string("angle of rotation (in degrees)"),
		  true, requires_argument);
Option<std::vector<float> > axis(string("-a,--axis"), vector<float>(3),
		  string("~<ax,ay,az>\tunnormalised axis vector (comma separated)"),
		  false, requires_argument);
Option<std::vector<float> > centre(string("-c,--centre"), vector<float>(3),
		  string("~<cx,cy,cz>\tcentre of rotation in mm (comma separated)"),
		  false, requires_argument);
Option<string> covopt(string("--cov"), string(""),
		  string("image filename used for centre of volume"),
		  false, requires_argument);
Option<string> outname(string("-o,--out"), string(""),
		  string("output filename for matrix"),
		  false, requires_argument);
int nonoptarg;

////////////////////////////////////////////////////////////////////////////

// Local functions


int do_work(int argc, char* argv[]) 
{
  ColumnVector ctr(3), ax(3);
  if (axis.unset()) {
    ax(1)=0; ax(2)=0; ax(3)=1;
  } else {
    ax(1) = axis.value()[0]; ax(2) = axis.value()[1]; ax(3) = axis.value()[2];
  }
  if (norm2(ax)==0.0) {
    cerr << "ERROR:: Axis cannot be zero" << endl;
    return 1;
  } else {
    ax = ax / sqrt(norm2(ax));
  }

  if (centre.unset()) {
    ctr(1)=0; ctr(2)=0; ctr(3)=0;
  } else {
    ctr(1) = centre.value()[0]; ctr(2) = centre.value()[1]; ctr(3) = centre.value()[2]; 
  }

  if (covopt.set()) {
    volume<float> v1;
    read_volume(v1,covopt.value());
    ctr(1) = (v1.xsize() -1.0) * v1.xdim() * 0.5;
    ctr(2) = (v1.ysize() -1.0) * v1.ydim() * 0.5;
    ctr(3) = (v1.zsize() -1.0) * v1.zdim() * 0.5;
  }

  Matrix mat(4,4);
  make_rot(ax*theta.value()*M_PI/180.0, ctr, mat);

  if (outname.set()) {
    write_ascii_matrix(mat,outname.value());
  } else {
    cout << mat << endl;
  }

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
    options.add(theta);
    options.add(axis);
    options.add(covopt);
    options.add(centre);
    options.add(outname);
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
    
  }  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } 

  // Call the local functions

  return do_work(argc,argv);
}

