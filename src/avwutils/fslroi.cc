//     fslroi.cc  extract cuboid ROI and/or timeseries from image
//     Stephen Smith, Matthew Webster and Mark Jenkinson, FMRIB Image Analysis Group
//     Copyright (C) 1999-2008 University of Oxford  
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

using namespace NEWIMAGE;

void print_usage(const string& progname) {
  cout << endl;
  cout << "Usage: fslroi <input> <output> <xmin> <xsize> <ymin> <ysize> <zmin> <zsize>" << endl;
  cout << "       fslroi <input> <output> <tmin> <tsize>\n" << endl;
  cout << "       fslroi <input> <output> <xmin> <xsize> <ymin> <ysize> <zmin> <zsize> <tmin> <tsize>" << endl;
  cout << "Note: indexing (in both time and space) starts with 0 not 1! Inputting -1 for a size will set it to the full image extent for that dimension." << endl;
}


template <class T>
int fmrib_main(int argc, char *argv[])
{
  volume4D<T> input_vol,output_vol;
  string input_name=string(argv[1]);
  string output_name=string(argv[2]);
  read_volume4D(input_vol,input_name);

  int minx = -1,miny = -1,minz = -1,mint = -1,maxx = -1,maxy = -1,maxz= -1,maxt = -1,argindex = 3;
  //Note there is currently no explicit bounds testing e.g. if maxx>inputvol.maxx() etc
  if (argc==5) //4D Timeseries only
  {
    minx=0;
    maxx=input_vol.maxx();
    miny=0;
    maxy=input_vol.maxy();
    minz=0;
    maxz=input_vol.maxz();
    mint=atoi(argv[argindex++]);
    maxt=atoi(argv[argindex++]);
    maxt+=mint-1;
  }
  else if (argc==9)  //3D Region of interest 
  {
    minx=atoi(argv[argindex++]);   //N.B. could compact some of these lines with the ones below...
    maxx=atoi(argv[argindex++]);
    miny=atoi(argv[argindex++]);
    maxy=atoi(argv[argindex++]);
    minz=atoi(argv[argindex++]);
    maxz=atoi(argv[argindex++]);
    mint=0;
    maxt=input_vol.maxt();
    maxx+=minx-1;
    maxy+=miny-1;
    maxz+=minz-1;
  }
  else if (argc==11)   //4D Timeseries and Region of interest
  {
    minx=atoi(argv[argindex++]);
    maxx=atoi(argv[argindex++]);
    miny=atoi(argv[argindex++]);
    maxy=atoi(argv[argindex++]);
    minz=atoi(argv[argindex++]);
    maxz=atoi(argv[argindex++]);
    mint=atoi(argv[argindex++]);
    maxt=atoi(argv[argindex++]); 
    maxx+=minx-1;
    maxy+=miny-1;
    maxz+=minz-1;
    maxt+=mint-1;
  }

  if ( maxx == minx-2 ) //Then the user input -1 for size
    maxx=input_vol.maxx();
  if ( maxy == miny-2 ) //Then the user input -1 for size
    maxy=input_vol.maxy();
  if ( maxz == minz-2 ) //Then the user input -1 for size
    maxz=input_vol.maxz();
  if ( maxt == mint-2 ) //Then the user input -1 for size
    maxt=input_vol.maxt();

  // sanity check
  if ((maxx<minx) || (maxy<miny) || (maxz<minz) || (maxt<mint)) { 
    imthrow("Invalid ROI dimensions",21); 
  }

  // now transform these coordinates to newimage conventions
  ColumnVector v(4);
  v << minx << miny << minz << 1.0;
  v = input_vol.niftivox2newimagevox_mat() * v;
  minx = MISCMATHS::round(v(1));
  miny = MISCMATHS::round(v(2));
  minz = MISCMATHS::round(v(3));
  v << maxx << maxy << maxz << 1.0;
  v = input_vol.niftivox2newimagevox_mat() * v;
  maxx = MISCMATHS::round(v(1));
  maxy = MISCMATHS::round(v(2));
  maxz = MISCMATHS::round(v(3));

  if (minx>maxx) { int tmpx=maxx; maxx=minx; minx=tmpx; }
  if (miny>maxy) { int tmpy=maxy; maxy=miny; miny=tmpy; }
  if (minz>maxz) { int tmpz=maxz; maxz=minz; minz=tmpz; }
  input_vol.setROIlimits(minx,miny,minz,mint,maxx,maxy,maxz,maxt);
  input_vol.activateROI();  //is this needed? //yes!
  output_vol=input_vol.ROI();
  save_volume4D(output_vol,output_name);
  return 0;
}


int main(int argc,char *argv[])
{

  Tracer tr("main");

  try {
    string progname=argv[0];
    if (argc !=5 && argc !=9 && argc !=11) 
      { 
	print_usage(progname);
	return 1; 
      }
    
    string iname=string(argv[1]);
    return call_fmrib_main(dtype(iname),argc,argv); 
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  }  catch(Exception &e) {
    exit(EXIT_FAILURE);
  } 
}


