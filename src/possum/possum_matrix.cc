
/*  POSSUM
    Ivana Drobnjak & Mark Jenkinson
    Copyright (C) 2005-2007 University of Oxford  */

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

// POSSUM

#include <iostream>
#include <string>
#include <fstream>
#include <unistd.h>

#ifdef USE_MPI
#include <mpi.h>
#include <unistd.h>
#endif //USE_MPI

#include "libprob.h"
#include "newmatap.h"
#include "newmatio.h"
#include "newimage/newimageall.h"
#include "possumfns.h"
#include "utils/options.h"
#include "newimage/costfns.h"
#include "miscmaths/miscmaths.h"

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

using namespace NEWIMAGE;
using namespace NEWMAT;
using namespace MISCMATHS;
using namespace Utilities;
//using namespace std;

string title="possum (Version 2.0)\nCopyright(c) 2008, University of Oxford (Ivana Drobnjak)";
string examples="possum_matrix  -p <pulse>  -m <motion file> -o <output main even matrix> [optional arguments]";

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);

Option<string> opt_motion(string("-m,--motion"), string(""),
		  string("<inputmatrix-filename> (Motion matrix [time(s) Tx(m) Ty(m) Tz(m) Rx(rad) Ry(rad) Rz(rad)]) "),
		  true, requires_argument);

Option<string> opt_pulse(string("-p,--pulse"), string(""),
		  string("<inputmatrix-basename> (Pulse sequence - all additional files .posx,.posy, etc,  expected to be in the same directory)"),
		  true, requires_argument);

Option<string> opt_mainmatrix(string("-o,--mainmatx"), string(""),
		  string("<outputmatrix-filename> (Main event matrix [t(s),rf_ang(rad),rf_freq_band(Hz),(4)=rf_cent_freq(Hz),read(1/0),Gx,Gy,Gz(T/m),Tx,Ty,Tz(m),angle_of_rot B(rad),rot_axis Bx,By,Bz(m),angle_of_rot A(rad),rot_axis Ax,Ay,Az(m)]) "),
		  true, requires_argument);

Option<bool> opt_old(string("--old"), false,
		  string("Allows for the old version of the sorter to run"),
		  false, no_argument);

Option<int> opt_seg(string("--seg"), 10000,
		  string("Seting the size of the segment of the matrix that is read in one at a time"),
		  false, requires_argument);

int nonoptarg;

/////////////////////////////////////////////////////////////////////////////////////////////////////
int compute_volume(int argc, char *argv[])
{
  ///////////////////////////////////////////////////////
  //PULSE & MOTION MATRIX SORT IN MAINMATRIX
  ///////////////////////////////////////////////////////
  cout<<"Reading the pulse sequence..."<<endl;
  RowVector pulseinfo;
  pulseinfo=read_ascii_matrix(opt_pulse.value()+".info");//[SeqType,TE,TR,TRslc,Nx,Ny,dx,dy,maxG,RiseT,BWrec, Nvol,Nslc,SlcThk,SlcDir,Gap,zstart,FlipAngle]
  cout<<"[SeqType,TE,TR,TRslc,Nx,Ny,dx,dy,maxG,RiseT,BW,Nvol,Nslc,SlcThk,SlcDir,Gap,zstart,FA]"<<endl;
  cout<<pulseinfo<<endl;
  cout<<""<<endl;
  cout<<"Reading the motion file..."<<endl;
  Matrix motion;
  int opt_test=0;
  if (verbose.value())opt_test=1;
  motion=read_ascii_matrix(opt_motion.value());
  cout<<"Motion file is "<<motion<<endl;
  if (opt_old.value()){
    Matrix pulse;
    pulse=read_binary_matrix(opt_pulse.value());
    pulse=sorter_old(pulse,motion);
    write_binary_matrix(pulse,opt_mainmatrix.value());
  }else{
    PMatrix pulse;
    read_binary_matrix(pulse,opt_pulse.value());
    cout<<"Gor here"<<endl;
    sorter(pulse,motion,opt_pulse.value(),opt_test);
    cout<<"Got here ?"<<endl;
    write_binary_matrix(pulse,opt_mainmatrix.value());
    RowVector numpoints(2);
    numpoints=0;
    numpoints(1)=pulse.Nrows();
    numpoints(2)=opt_seg.value();
    cout<<"Size of the matrix="<<numpoints(1)<<" rows"<<endl;
    cout<<"Size of the matrix segments="<<numpoints(2)<<" rows"<<endl; 
    write_ascii_matrix(numpoints,opt_mainmatrix.value()+".numpoints");
  }
  return 0;
}


int main (int argc, char *argv[])
{

  Tracer tr("main");
  OptionParser options(title, examples);

  try {
    options.add(verbose);
    options.add(help);
    options.add(opt_motion);
    options.add(opt_pulse);
    options.add(opt_old);
    options.add(opt_mainmatrix);
    options.add(opt_seg);

    nonoptarg = options.parse_command_line(argc, argv);

    // line below stops the program if there are less than 2 non-optional args
    //   or the help was requested or a compulsory option was not set
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
  compute_volume(argc, argv);
  return 0;
}
