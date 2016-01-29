//     fslcorrecthd.cc - check and correct a nifti file for bad vox-offset
//     Matthew Webster, FMRIB Image Analysis Group
//     Copyright (C) 2007 University of Oxford  
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
#include "fslio/fslio.h"
#include <iostream>
using namespace NEWIMAGE;

void print_usage(const string& progname) 
{
  cout << endl;
  cout << "Usage: fslcorrecthd <input> <output>" << endl;
  cout << "       Note that fslcorrecthd only operates on uncompressed NIFTI or ANALYZE files" << endl;
}

int main(int argc,char *argv[])
{
  if (argc < 3) 
  {
    print_usage(string(argv[0]));
    return 1; 
  }
  FSLIO* fslio=NULL;
  fslio = FslOpen(FslMakeBaseName(argv[1]),"rb");
  FslClose(fslio);
  struct dsr *hdr;
  hdr = (struct dsr *)calloc(1,sizeof(struct dsr));
  FslReadRawHeader(hdr,fslio->niftiptr->fname);
  if (fslio->niftiptr->byteorder != nifti_short_order()) 
  {
    cout << "Byte swapping" << endl;
    AvwSwapHeader(hdr);
  } 
  //check nifti-libs output versus raw header info
  int offset =(int) ( fslio->niftiptr->iname_offset - hdr->dime.vox_offset );
  int minft=(int)MIN(fslio->niftiptr->iname_offset,hdr->dime.vox_offset);
  cout << "number of bytes wrong: " << offset << endl << "start at byte location: " << minft << endl;

  if (offset==0) {
    cout << "No byte correction needed, exiting." << endl;
    return 0;
  }

  if (FslIsCompressedFileType(FslGetFileType(fslio))) {
    cerr << "Error: fslcorrecthd requires uncompressed input" << endl;         
    return 1;
  }

 ifstream input_file;
 ofstream output_file;
 char *temp,*outputName,*inputName; 
 FslGetHdrImgNames(argv[2],fslio,&temp,&outputName);
 FslGetHdrImgNames(argv[1],fslio,&temp,&inputName);
 char byte[1];
 input_file.open(inputName,ios::in | ios :: binary);
 output_file.open(outputName,ofstream::out | ofstream::binary);

 for(int i=1;i<=minft;i++) //Write Header
 {
   input_file.read(byte,1);
   if (input_file.eof()) break;
   output_file.write(byte,1);
 }

 for(int i=1;i<=abs(offset) && offset>0;i++) //Pad if we have missing 4 bytes
 {
   byte[0]=0;
   output_file.write(byte,1);
 }

 for(int i=1;i<=abs(offset) && offset<0;i++) //Read past bad extensions/junk 
 {
   input_file.read(byte,1);
 }

 while(true) //Copy the data
 {
   input_file.read(byte,1);
   if (input_file.eof()) break;
   output_file.write(byte,1);
 }  
 output_file.close();
 input_file.close();

 system(("FSLOUTPUTTYPE=NIFTI; ${FSLDIR}/bin/fslmaths " + string(outputName)).c_str());  //To clean up header
 free(temp);
 free(outputName);
 free(hdr);
 return 0;
}


