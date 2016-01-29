/*  Copyright (C) 2010 University of Oxford  */
/*  Stam Sotiropoulos  */    
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
#include <string>
#include <cmath>
#include "newimage/newimageall.h"

using namespace std;
using namespace NEWIMAGE;
using namespace NEWMAT;

//Rearranges a dataset according to the angular distance between the corresponding bvecs entry and a reference vector  
int main ( int argc, char *argv[]){
  if(argc<5 || argc>6){
    cerr<<" "<<endl;
    cerr<<"usage: rearrange data bvecs ref_dyads [brain_mask] output"<<endl<<endl;
    cerr<<"Rearranges a dataset according to the angular distance"<<endl; 
    cerr<<"between the corresponding bvecs entry and a reference vector"<<endl;
    cerr<<" "<<endl;
    exit(1);
  }
  
  ColumnVector ref_vector(3);
  Matrix bvecs;   //Read Input 
  volume4D<float> data, ref; volume<float> mask;
  
  read_volume4D(data,argv[1]);
  bvecs=read_ascii_matrix(argv[2]);
  read_volume4D(ref,argv[3]);
  if (argc==6)
    read_volume(mask,argv[4]);
  else{ 
    mask.reinitialize(data.xsize(),data.ysize(),data.zsize());
    mask=1;
  }
  
  volume4D<float> output(data.xsize(),data.ysize(),data.zsize(),data.tsize());
  copybasicproperties(data,output);output=0;

  if(bvecs.Nrows()>3) bvecs=bvecs.t();   //Make sure the bvecs entries are normalized unit vectors
  for(int i=1;i<=bvecs.Ncols();i++){
    float tmpsum=sqrt(bvecs(1,i)*bvecs(1,i)+bvecs(2,i)*bvecs(2,i)+bvecs(3,i)*bvecs(3,i));
    if(tmpsum!=0){
      bvecs(1,i)=bvecs(1,i)/tmpsum;
      bvecs(2,i)=bvecs(2,i)/tmpsum;
      bvecs(3,i)=bvecs(3,i)/tmpsum;
    }  
  }
  for(int z=data.minz();z<=data.maxz();z++){   //Rearrange data for each voxel
    for(int y=data.miny();y<=data.maxy();y++){
      for(int x=data.minx();x<=data.maxx();x++){
	if (mask(x,y,z)!=0){
	  vector<pair<float,int> > dot;  //Keep in the first entry of each pair the dot product value and the in the second the index 
	  pair<float,int> ftmp;
	  ref_vector(1)=ref(x,y,z,0); 	  ref_vector(2)=ref(x,y,z,1); 	  ref_vector(3)=ref(x,y,z,2);
	  for (int n=1; n<=bvecs.Ncols(); n++){     //Get the dot product of each bvecs entry with the reference orientation in this voxel
	    if (bvecs(1,n)==0 && bvecs(2,n)==0 && bvecs(3,n)==0)
	      ftmp.first=2.0;              //The b=0 entries are first in the sequence
	    else
	      ftmp.first=fabs(ref_vector(1)*bvecs(1,n)+ref_vector(2)*bvecs(2,n)+ref_vector(3)*bvecs(3,n));
	    ftmp.second=n-1;
	    dot.push_back(ftmp);
	    } 
	  sort(dot.begin(),dot.end());     //Sort the dot products in ascending order
	  reverse(dot.begin(),dot.end());  //Reverse the ordering, so that it is in descending order
	  for (int n=0; n<bvecs.Ncols(); n++)//Get the dot product of each bvecs entry with the reference orientation in this voxel
	    output(x,y,z,n)=data(x,y,z,dot[n].second);
	}
      }
    }
    cout<<z+1<<" slices processed"<<endl;
  }
  save_volume4D(output,argv[argc-1]);
  return 0;
}









