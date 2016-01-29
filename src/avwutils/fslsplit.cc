//     fslsplit.cc - split 4D files into 3D files for SPM
//     David Flitney and Matthew Webster, FMRIB Image Analysis Group
//     Copyright (C) 2000-2006 University of Oxford  
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

#include "newimage/newimageall.h"
#include "newimage/fmribmain.h"
#include "miscmaths/miscmaths.h"

using namespace NEWIMAGE;
using namespace MISCMATHS;

void print_usage(const string& progname) 
{
  cout << endl;
  cout << "Usage: fslsplit <input>" << endl;
  cout << "       fslsplit <input> [output_basename] [-t/x/y/z]" << endl;
  cout << "       -t : separate images in time (default behaviour)" << endl;
  cout << "       -x : separate images in the x direction"  << endl;
  cout << "       -y : separate images in the y direction"  << endl;
  cout << "       -z : separate images in the z direction" << endl;
}

template <class T>
int fmrib_main(int argc, char *argv[])
{
  volume4D<T> input_vol,output_vol;
  string input_name=string(argv[1]);
  string output_name("vol");
  int xoff(1),yoff(1),zoff(1),toff(1),nsize;
  read_volume4D(input_vol,input_name);
  if ((argc>2) && (argv[2][0]!='-')) output_name=string(argv[2]);
  if (argv[argc-1][0] == '-')  
  { 
      if (argv[argc-1][1] == 't') toff=input_vol.tsize(); 
      if (argv[argc-1][1] == 'z') zoff=input_vol.zsize();  
      if (argv[argc-1][1] == 'y') yoff=input_vol.ysize();  
      if (argv[argc-1][1] == 'x') xoff=input_vol.xsize();  
  }
  else toff=input_vol.tsize(); 
  nsize=xoff*yoff*zoff*toff;
  for(int j=0;j<nsize;j++)   
  {
    input_vol.setROIlimits(0+j*(xoff!=1),0+j*(yoff!=1),0+j*(zoff!=1),0+j*(toff!=1),input_vol.xsize()-1-(nsize-j-1)*(xoff!=1),input_vol.ysize()-1-(nsize-j-1)*(yoff!=1),input_vol.zsize()-1-(nsize-j-1)*(zoff!=1),input_vol.tsize()-1-(nsize-j-1)*(toff!=1));
  input_vol.activateROI(); 
  output_vol=input_vol.ROI();
  save_volume4D(output_vol,(output_name+num2str(j,4)));
  }
  return 0;
}


int main(int argc,char *argv[])
{

  Tracer tr("main");

  string progname=argv[0];
  if (argc <= 1 || argc >= 5) 
  { 
    print_usage(progname);
    return 1; 
  }
   
  string iname=string(argv[1]);
  return call_fmrib_main(dtype(iname),argc,argv); 
}

