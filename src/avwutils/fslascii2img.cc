//     fsl2ascii.cc - convert raw ASCII text to a NIFTI image
//     Mark Jenkinson, FMRIB Image Analysis Group
//     Copyright (C) 2009 University of Oxford  
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
#include "miscmaths/miscmaths.h"
#include <fstream>

using namespace NEWIMAGE;
using namespace MISCMATHS;

void print_usage(const string& progname) {
  cout << endl;
  cout << "Usage: fslascii2img <input> <xsize> <ysize> <zsize> <tsize> <xdim> <ydim> <zdim> <TR>  <output>" << endl;
  cout << "  where sizes are in voxels, dims are in mm, TR in sec " << endl;
}


int do_work(int argc, char *argv[])
{
  int xsize, ysize, zsize, tsize;
  float xdim, ydim, zdim, tr;
  xsize=atoi(argv[2]);
  ysize=atoi(argv[3]);
  zsize=atoi(argv[4]);
  tsize=atoi(argv[5]);
  xdim=atof(argv[6]);
  ydim=atof(argv[7]);
  zdim=atof(argv[8]);
  tr=atof(argv[9]);
  volume4D<float> ovol(xsize,ysize,zsize,tsize);
  ovol.setdims(xdim,ydim,zdim,tr);
  string input_name=string(argv[1]);
  string output_name=string(argv[10]);
  Matrix amat;
  amat = read_ascii_matrix(input_name);
  if (xsize*ysize*zsize*tsize != amat.Nrows() * amat.Ncols()) {
    cerr << "Sizes incompatible: " <<  xsize*ysize*zsize*tsize << " voxels vs " << amat.Nrows() * amat.Ncols() << " numbers" << endl;
    cerr << "Matrix dimensions are " << amat.Nrows() << " by " << amat.Ncols() << endl;
    exit(EXIT_FAILURE);
  }
  amat = reshape(amat.t(),tsize,xsize*ysize*zsize);
  ovol.setmatrix(amat);
  save_volume4D(ovol,output_name);
  return 0;
}


int main(int argc,char *argv[])
{

  Tracer tr("main");

  string progname=argv[0];
  if (argc != 11) 
  { 
    print_usage(progname);
    return 1; 
  }
   
  return do_work(argc,argv); 
}


