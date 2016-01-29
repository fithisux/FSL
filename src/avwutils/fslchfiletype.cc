//     fslchfiletype.cc  Conversion program for different volume types
//     Matthew Webster, FMRIB Image Analysis Group
//     Copyright (C) 2008 University of Oxford  
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

int printUsage(const string& progname) 
{
  cout << "Usage: " << progname << " <filetype> <filename> [filename2]" << endl << endl;
  cout << "  Changes the file type of the image file, or copies to new file" << endl;
  cout << "  Valid values of filetype are ANALYZE, NIFTI, NIFTI_PAIR," << endl;
  cout << "                               ANALYZE_GZ, NIFTI_GZ, NIFTI_PAIR_GZ" << endl;
  return 1;
}

template <class T>
int fmrib_main(int argc, char *argv[])
{
  volume4D<T> input;
  read_volume4D(input,string(argv[0]));
  save_volume4D_dtype(input,string(argv[1]),dtype(string(argv[0])));
  return 0;
}

int main(int argc,char *argv[])
{
  if (argc < 3 || argc > 4)
    return printUsage(string(argv[0]));

  string inputFile(argv[2]);
  string outputFile(argv[2]);
  if (argc==4) outputFile=string(argv[3]);

  int outputType(-1);
  if (string(argv[1])=="ANALYZE") outputType=FSL_TYPE_ANALYZE;
  else if (string(argv[1])=="NIFTI") outputType=FSL_TYPE_NIFTI;
  else if (string(argv[1])=="NIFTI_PAIR") outputType=FSL_TYPE_NIFTI_PAIR;
  else if (string(argv[1])=="ANALYZE_GZ") outputType=FSL_TYPE_ANALYZE_GZ;
  else if (string(argv[1])=="NIFTI_GZ") outputType=FSL_TYPE_NIFTI_GZ;
  else if (string(argv[1])=="NIFTI_PAIR_GZ") outputType=FSL_TYPE_NIFTI_PAIR_GZ;
  else return printUsage(string(argv[0]));

  FslSetOverrideOutputType(outputType);

  if (dtype(inputFile)==DT_COMPLEX) 
  {
    volume4D<float> real, imaginary;
    read_complexvolume4D(real, imaginary, inputFile);
    save_complexvolume4D(real, imaginary, outputFile);
  }
  else  
  {
    argv[0]= const_cast <char*>(inputFile.c_str());
    argv[1]= const_cast <char*>(outputFile.c_str());
    call_fmrib_main(dtype(inputFile),argc,argv);
  }
  return 0;
}

