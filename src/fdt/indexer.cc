/*    Copyright (C) 2012 University of Oxford  */

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


#include <iostream>
#include <fstream>
#include <iomanip>
#include "newimage/newimageall.h"
#include "utils/log.h"



using namespace std;
using namespace NEWIMAGE;
using namespace Utilities;
//using namespace NEWMAT;
//////////////////////////
/////////////////////////

void read_masks(vector<string>& masks, const string& filename){
  ifstream fs(filename.c_str());
  string tmp;
  if(fs){
    fs>>tmp;
    while(!fs.eof()){
      masks.push_back(tmp);
      fs>>tmp;
    }
  }
  else{
    cerr<<filename<<" does not exist"<<endl;
    exit(0);
  }
}


ColumnVector vox_to_vox(const ColumnVector& xyz1,const ColumnVector& dims1,const ColumnVector& dims2,const Matrix& xfm){
  ColumnVector xyz1_mm(4),xyz2_mm,xyz2(3);
  xyz1_mm<<xyz1(1)*dims1(1)<<xyz1(2)*dims1(2)<<xyz1(3)*dims1(3)<<1;
  xyz2_mm=xfm*xyz1_mm;
  xyz2_mm=xyz2_mm/xyz2_mm(4);
  xyz2<<xyz2_mm(1)/dims2(1)<<xyz2_mm(2)/dims2(2)<<xyz2_mm(3)/dims2(3);
  return xyz2;
}


int main ( int argc, char **argv ){
  if(argc<6){
    cerr<<"Usage: indexer <mask_for_going_ahead> <file_of_volume_names> <x> <y> <z>"<<endl;
    cerr<<endl;
    cerr<< "<x> <y> <z> in mni mm coords"<<endl;
    return 0;
  }
 
  vector<string> masknames;
  volume<int> inside_mask;
  read_volume(inside_mask,argv[1]);
  read_masks(masknames,argv[2]);
  vector<volume<float> > cortex_masks;
  volume<float> tmpcort;
  for( unsigned int m = 0; m < masknames.size(); m++ ){
    read_volume(tmpcort,masknames[m]);
    cortex_masks.push_back(tmpcort);
  }
  
  float x_mm=atof(argv[3]);
  float y_mm=atof(argv[4]);
  float z_mm=atof(argv[5]);
  float x_orig_mm=90;
  float y_orig_mm=124;
  float z_orig_mm=72;
  int x_roi_start_vox=45;
  int y_roi_start_vox=70;
  int z_roi_start_vox=50;
  int numsubjects=11;
  float xvox=(x_mm+x_orig_mm)/cortex_masks[0].xdim()-x_roi_start_vox;
  float yvox=(y_mm+y_orig_mm)/cortex_masks[0].ydim()-y_roi_start_vox;
  float zvox=(z_mm+z_orig_mm)/cortex_masks[0].zdim()-z_roi_start_vox;
  
  if(inside_mask((int)MISCMATHS::round(xvox),(int)MISCMATHS::round(yvox),(int)MISCMATHS::round(zvox))==0){
    cout<<"Sorry - Your input is not in our defined thalamus"<<endl;
    for(unsigned int i=0;i<masknames.size();i++){
      cout<<0<<endl;
    }
  }  
  else{
    cout<<"Input inside thalamus"<<endl;
    for(unsigned int i=0;i<masknames.size();i++){
      if(cortex_masks[i].interpolate(xvox,yvox,zvox)==0){
	cout<<std::setprecision(0)<<std::fixed;
      }
      else{
	cout<<std::setprecision(2)<<std::fixed;
      }
      cout<<cortex_masks[i].interpolate(xvox,yvox,zvox)/numsubjects<<endl;
    }
  }
  return 0;
}


 














