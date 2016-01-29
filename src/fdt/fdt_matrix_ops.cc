/*  Copyright (C) 2004 University of Oxford  */

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
#include <cmath>
#include "newimage/newimageall.h"
#include <vector>
#include <algorithm>

using namespace std;
using namespace NEWIMAGE;
using namespace NEWMAT;

string matf2coordf(string matf){
  size_t pos=matf.rfind("/");
  if(pos!=string::npos)
    matf.replace(pos,1,"/coords_for_");
  else
    matf="coords_for_"+matf;

  return matf;
}


int main ( int argc, char **argv ){
  if(argc<5){
    cout<<"usage: fdt_matrix_ops <matrix1> <seedmask1> <matrix2> <seedmask2> ... [inclusion mask] <output>"<<endl;
    cout<<"creates one big uber-matrix from lots of cute wee dinky ones"<<endl;
    cout<<"If seedmasks overlap, it just takes the values from the first of them"<<endl;
    cout<<"If you specify an incluson mask (optional), only voxels that are"<<endl;
    cout<<"inside this mask will be included in the output"<<endl;
    exit(0);
  }
  int var=argc-1;
  bool incmaskyn = (float(var)/2.0)==int(float(var)/2.0);
  int Nmats;
  if(incmaskyn) Nmats=var/2-1;
  else Nmats=(var-1)/2;
  
  vector<volume<int>* > masks;
  vector<volume<float>* > mats;
  vector<volume<int>* > coords;
  volume<int> incmask;
  volume<int> totalmask;
  for(int i=0;i<Nmats;i++){
    volume<float> *tmp= new volume<float>;
    volume<int> *tmpcoords= new volume<int>;
    volume<int> *tmpmask= new volume<int>;
    string matfile=string(argv[ 2*i + 1 ]);
    string coordfile=matf2coordf(matfile);
    read_volume(*tmp,matfile);
    read_volume(*tmpcoords,coordfile);
    read_volume(*tmpmask,argv[ 2*(i + 1) ]);
    mats.push_back(tmp);
    masks.push_back(tmpmask);
    coords.push_back(tmpcoords);
    if(i==0) totalmask=*tmpmask;
    else totalmask=totalmask+ *tmpmask;

}
  if(incmaskyn){
    read_volume(incmask,argv[var-1]);
    incmask=incmask*totalmask;
  }else{
    incmask=totalmask;
  }
  
  vector<volume<int>* >lookups;
  for(int i=0;i<Nmats;i++){
    volume<int> *lu = new volume<int>;
    *lu=*masks[i];int conrow=0;
    for(int z=0;z<(*lu).zsize();z++){
      for(int y=0;y<(*lu).ysize();y++){
	for(int x=0;x<(*lu).xsize();x++){
	  if((*masks[i])(x,y,z)>0){
	    (*lu)(x,y,z)=conrow;
	    conrow++;
	  }
	}
      }
    }
    lookups.push_back(lu);
  }
  
  int nvoxels=0;

  for(int z=0;z<incmask.zsize();z++) {
    for(int y=0;y<incmask.ysize();y++){
      for(int x=0;x<incmask.xsize();x++){
	  if(incmask(x,y,z)>0) nvoxels++;
      }
    }
  }
  
  volume<float> output(nvoxels,(*mats[0]).ysize(),1);
  volume<float> outcoords(nvoxels,3,1);
  int newrow=0;
  for(int z=0;z<incmask.zsize();z++) {
    for(int y=0;y<incmask.ysize();y++){
      for(int x=0;x<incmask.xsize();x++){
	if(incmask(x,y,z)>0){
	  bool found=false;
	  for(unsigned int i=0;i<mats.size();i++){
	    if(!found){
	      if((*masks[i])(x,y,z)>0){
		int oldrow=(*lookups[i])(x,y,z);
		for(int col=0;col<(*mats[i]).ysize();col++){
		  output(newrow,col,0)=(*mats[i])(oldrow,col,0);
		}
		//cout<<x<<" "<<y<<" "<<z<<endl;
		//cout<<newrow<<" "<<oldrow<<endl;
		//cout<<(*coords[i])(oldrow,0,0)<<" "<<(*coords[i])(oldrow,1,0)<<" "<<(*coords[i])(oldrow,2,0)<<endl;
		outcoords(newrow,0,0)=(*coords[i])(oldrow,0,0);
		outcoords(newrow,1,0)=(*coords[i])(oldrow,1,0);
		outcoords(newrow,2,0)=(*coords[i])(oldrow,2,0);
		//cout<<"yep"<<endl;
		found=true;
		newrow++;
	      }
	    }
	    
	  }
	  
	}
      }
    }
  }

  for(int i=0;i<Nmats;i++){
    delete mats[i];
    delete coords[i];
    delete masks[i];
    delete lookups[i];
}
  save_volume(output,argv[argc-1]);
  string coordout=matf2coordf(string(argv[argc-1]));
  save_volume(outcoords,coordout);
 return 0;
}
 


















