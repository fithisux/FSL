
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

// TESTING_POSSUM

#include <iostream>
#include <string>
#include <fstream>
#include <unistd.h>

#ifdef USE_MPI
#include <mpi.h>
#endif //USE_MPI

#include "libprob.h"
#include "newmatap.h"
#include "newmatio.h"
#include "newimage/newimageall.h"
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

string title="test_possum\nCopyright(c) 2003, University of Oxford (Ivana Drobnjak)";
string examples="test_possum ";

Option<bool>   verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool>   help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);

//INPUT for the determination of testing a specific module 
Option<string> opt_mod(string("--mod"), "",
		  string("Defines which modules is to be tested: implementation, image_contrast, motion, b0, chemical_shift, rf, noise,eddy_currents, bold "),
		  true,requires_argument);

//INPUT object
Option<string> opt_object(string("-i,--inp"), string("el"),
		  string("ellipse (el) or rectangle (rec)"),
		  false, requires_argument);
Option<float>  opt_a(string("--a"),0.03,
		  string("(mm) for rec: side1 a, for el: semimajoraxis"),
		  false,requires_argument);
Option<float>  opt_b(string("--b"),0.03,
		  string("(mm) for rec: side2 b, for el: semiminoraxis"),
		  false,requires_argument);
Option<float>  opt_x0(string("--x0"),0,
		  string("(mm) center of the object"),
		  false,requires_argument);
Option<float>  opt_y0(string("--y0"),0,
		  string("(mm) center of the object"),
		  false,requires_argument);
Option<float>  opt_theta(string("--theta"),0,
		  string("(degrees) for rec: angle between side1 and x axis, for el: angle between majoraxis and x axis"),
		  false,requires_argument);

//INPUT for the kspace coordinates
Option<string> opt_kcoord(string("-k,--kcoord"), string(""),
		  string("2-raw kcoord matrix [kx;ky]"),
		  true, requires_argument);
//INPUT for the kspace coordinates
Option<string> opt_pulse(string("-p,--pulse"), string(""),
		  string("pulseseq"),
		  false, requires_argument);

//INPUT b0
Option<float>  opt_b0(string("--b0"),0,
		  string("(T) b0 perturbation = b0+b0x*x+b0y*y+b0z*z"),
		  false, requires_argument);
Option<float>  opt_b0x(string("--b0x"),0,
		  string("b0 perturbation = b0+b0x*x+b0y*y+b0z*z"),
		  false, requires_argument);
Option<float>  opt_b0y(string("--b0y"),0,
		  string("b0 perturbation = b0+b0x*x+b0y*y+b0z*z"),
		  false, requires_argument);
Option<float>  opt_b0z(string("--b0z"),0,
		  string("b0 perturbation = b0+b0x*x+b0y*y+b0z*z"),
		  false, requires_argument);

//INPUT rf
Option<float>  opt_RFrec(string("--rfr"),1,
		  string("(val  0 to 1, 1 is for perfectly homog) signal=signal*rfr"),
		  false, requires_argument);
Option<float>  opt_RFtrans(string("--rft"),1,
		  string("(val  0 to 1, 1 is for perfectly homog) flip_angle=flip_angle*rft"),
		  false, requires_argument);

//INPUT motion
Option<string> opt_motion(string("-m,--motion"), string(""),
		  string("7-col motion matrix [time(s) Tx(m) Ty(m) Tz(m) Rx(rad) Ry(rad) Rz(rad)] "),
		  true, requires_argument);

//INPUT bold
//INPUT eddys
//INPUT chemical shift


//OUTPUT signal
Option<string> opt_signal(string("-o,--out"), string(""),
		  string("2-row signal matrix, [sreal; simag]"),
		  true, requires_argument);

int nonoptarg;

///////////////////////////////////
double bessj1(const double x){
  
  double ax,z,xx,y,ans,ans1,ans2;
  
  if ((ax=fabs(x)) <8.0){
    y=x*x;
    ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1+y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
    ans2=144725228442.0+y*(2300535178.0+y*(18583304.74+y*(99447.43394+y*(376.9991397+y*1.0))));
    ans=ans1/ans2;
  } else {
    z=8.0/ax;
    y=z*z;
    xx=ax-2.356194491;
    ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4+y*(0.2457520174e-5+y*(-0.240337019e-6))));
    ans2=0.04687499995+y*(-0.2002690873e-3+y*(0.8449199096e-5+y*(-0.88228987e-6+y*0.105787412e-6)));
    ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
    if (x < 0.0) ans=-ans;
  }
  return ans;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
int implementation(){ 
  
  cout<<"Starting TESTING_implementation..."<<endl;
  //////////////////////////////////////////////////////////////////////////
  // READ IN THE OBJECT                                      
  //////////////////////////////////////////////////////////////////////////
  string obj=opt_object.value();
  cout<<obj<<endl;
  float a=opt_a.value()*1e-3;
  float b=opt_b.value()*1e-3; 
  float x0=opt_x0.value()*1e-3;
  float y0=opt_y0.value()*1e-3;
  float theta=opt_theta.value()*M_PI/180;
  ///////////////////////////////////////////////////////
  //KSPACE COORD, MOTION MATRIX
  ///////////////////////////////////////////////////////
  cout<<"Reading pulse sequence kspace ..."<<endl;
  Matrix kcoord;
  kcoord=read_binary_matrix(opt_kcoord.value());
  ///////////////////////////////////////////////////////
  //PULSE INFO
  ///////////////////////////////////////////////////////
  RowVector pulseinfo(17);
  pulseinfo=read_ascii_matrix(opt_pulse.value()+".info");
  cout<<opt_pulse.value()+"info"<<endl;
  cout<<read_ascii_matrix(opt_pulse.value()+"info")<<endl;
  //int numvol=(int) (pulseinfo(12));
  //int numslc=(int) (pulseinfo(13));
  double slcthk=pulseinfo(14);//slcthk (m)
  cout<<"slcthk  "<<slcthk<<endl;
  //double dx=pulseinfo(7);
  //double dy=pulseinfo(8);
  int resX=(int) (pulseinfo(5));
  int resY=(int) (pulseinfo(6));
  //int zstart_p=(int) (pulseinfo(17));// zstart (vox)
  //int seq=(int) (pulseinfo(1));//epi or ge
  ////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  // SIGNAL                                                             
  /////////////////////////////////////////////////////////////////////////
  int nreadp=resX*resY;
  cout<<"Number of read out points is "<<nreadp<<endl;
  Matrix signal(2,nreadp);//two rows, one real and one complex for the signal in sum
  double fix=slcthk*1e06;//a_signal is multiplied with this so that it matches difference with possum, it iwll need few more fixes though but will put them in post procesing
  signal=0;
  /////////////////////////////////////////////////////////
  //MAIN WORK
  /////////////////////////////////////////////////////////
  if (obj=="rec"){
    for (int m=1;m<=nreadp;m++){
      signal(1,m)=fix*a*b*Sinc((kcoord(1,m)*cos(theta)-kcoord(2,m)*sin(theta))*a)*Sinc((kcoord(1,m)*sin(theta)+kcoord(2,m)*cos(theta))*b)*cos(2*M_PI*(kcoord(2,m)*y0-kcoord(1,m)*x0));
      signal(2,m)=fix*a*b*Sinc((kcoord(1,m)*cos(theta)-kcoord(2,m)*sin(theta))*a)*Sinc((kcoord(1,m)*sin(theta)+kcoord(2,m)*cos(theta))*b)*sin(2*M_PI*(kcoord(2,m)*y0-kcoord(1,m)*x0));
    }

  }
  else if (obj=="el"){
    for (int m=1;m<=nreadp;m++){
      float kr=sqrt((a*kcoord(1,m)*cos(theta)+a*kcoord(2,m)*sin(theta))*(a*kcoord(1,m)*cos(theta)+a*kcoord(2,m)*sin(theta))+(b*kcoord(2,m)*cos(theta)-b*kcoord(1,m)*sin(theta))*(b*kcoord(2,m)*cos(theta)-b*kcoord(1,m)*sin(theta)));
      signal(1,m)=fix*(a*b*bessj1(2.0*M_PI*kr)/kr)*cos(2*M_PI*(-kcoord(2,m)*y0 - kcoord(1,m)*x0));
      signal(2,m)=fix*(a*b*bessj1(2.0*M_PI*kr)/kr)*sin(2*M_PI*(-kcoord(2,m)*y0 - kcoord(1,m)*x0));
      if (m==2081){
        cout.precision(20);
	cout<<"fix= "<<fix<<"a= "<<a<<"b= "<<b<<"bessj1(2.0*M_PI*kr)/kr= "<<bessj1(2.0*M_PI*kr)/kr<<endl;
        cout<<"cos(2*M_PI*(-kcoord(2,m)*y0 - kcoord(1,m)*x0))= "<<cos(2*M_PI*(-kcoord(2,m)*y0 - kcoord(1,m)*x0))<<endl;
	cout<<"sin(2*M_PI*(-kcoord(2,m)*y0 - kcoord(1,m)*x0))= "<<sin(2*M_PI*(-kcoord(2,m)*y0 - kcoord(1,m)*x0))<<endl;
	cout<<"sreal(2081)= "<<signal(1,2081)<<endl;
        cout<<"simag(2081)= "<<signal(2,2081)<<endl;
      }
    }
  }
   /////////////////////////
   //OUTPUT
   /////////////////////////
    write_binary_matrix(signal,opt_signal.value());
    cout<<"Testing_implementation finished generating the signal..."<<endl;

  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////
int compute_volume(int argc, char *argv[])
{
  
  if (opt_mod.value()=="implementation") implementation(); 
  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////
int main (int argc, char *argv[])
{

  Tracer tr("main");
  OptionParser options(title, examples);

  try {
    options.add(verbose);
    options.add(help);
    options.add(opt_object);
    options.add(opt_a);
    options.add(opt_b);
    options.add(opt_x0);
    options.add(opt_y0);
    options.add(opt_theta);
    options.add(opt_motion);
    options.add(opt_kcoord);
    options.add(opt_pulse);
    options.add(opt_b0);
    options.add(opt_b0x);
    options.add(opt_b0y);
    options.add(opt_b0z);
    options.add(opt_RFrec);
    options.add(opt_RFtrans);
    options.add(opt_mod);
    options.add(opt_signal);
   
    
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
