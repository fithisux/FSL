//     fslcpgeom.cc - Copy certain parts of an AVW header
//     Mark Jenkinson, Steve Smith and Matthew Webster, FMRIB Image Analysis Group
//     Copyright (C) 2001-2005 University of Oxford  
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

using namespace NEWIMAGE;

void print_usage(const string& progname) 
{
  cout << endl;
  cout << "Usage: fslcpgeom <source> <destination> [-d]" << endl;
  cout << "-d : don't copy image dimensions" << endl;
}


int avwcpgeom_main(int argc, char *argv[])
{
  FSLIO *src = NULL, *dest = NULL, *destcopy = NULL;
  short x, y, z, v, copydim=1, t, scode, qcode, dt=-1;
  float vx, vy, vz, tr;
  int filetype;
  mat44 smat, qmat;
  void *buffer = NULL;
  char desthdrname[10000];
  size_t nbytes, nsrcbytes, ndestbytes;

  if (argc>3)
    copydim=0;
    
  src = FslOpen(argv[1], "r");

  dest = FslOpen(argv[2], "r");
  destcopy = (FSLIO *)calloc(sizeof(FSLIO),1);
  FslCloneHeader(destcopy,dest);

  if ((src == NULL) || (dest == NULL)) {
    perror("Error opening files");
    return EXIT_FAILURE;
  }

  FslGetDim(src, &x, &y, &z, &v);
  nsrcbytes = x * y * z * v * (FslGetDataType(dest, &t) / 8);
  FslGetDim(dest, &x, &y, &z, &v); 
  ndestbytes = x * y * z * v * (FslGetDataType(dest, &t) / 8);
  if (nsrcbytes > ndestbytes) nbytes=nsrcbytes; else nbytes=ndestbytes;
  if( (buffer = calloc(nbytes,1)) == NULL ) {
    perror("Unable to allocate memory for copy");
    return EXIT_FAILURE;
  }
  FslReadVolumes(dest, buffer, v);


  strcpy(desthdrname,dest->niftiptr->fname);
  filetype = FslGetFileType(dest);
  FslGetDataType(dest,&dt);
  FslClose(dest);
  dest = FslXOpen(desthdrname, "wb", filetype);
  FslCloneHeader(dest,destcopy); 

  scode = FslGetStdXform(src,&smat);
  FslSetStdXform(dest,scode,smat);
  qcode = FslGetRigidXform(src,&qmat);
  FslSetRigidXform(dest,qcode,qmat);
  
  FslGetVoxDim(src, &vx, &vy, &vz, &tr);
  FslSetVoxDim(dest, vx, vy, vz, tr);
  
  if (copydim) {
    FslGetDim(src, &x, &y, &z, &v);
    FslSetDim(dest, x, y, z, v);
  }


  /* Preserve the datatype - probably unneccesary now, but left for safety */ 
  FslSetDataType(dest,dt);

 FslWriteHeader(dest);
 FslWriteVolumes(dest, buffer, v);
 FslClose(dest);

 FslClose(src);

  return 0;
}


int main(int argc,char *argv[])
{
  if (argc < 3) 
  {
    print_usage(string(argv[0]));
    return 1; 
  }
  return avwcpgeom_main(argc,argv); 
}


