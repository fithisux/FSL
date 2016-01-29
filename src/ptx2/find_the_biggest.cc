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
#include "meshclass/meshclass.h"
#include "newimage/newimageall.h"
#include <vector>
#include "csv_mesh.h"
using namespace std;
using namespace NEWIMAGE;
using namespace mesh;

void biggest_from_volumes(vector<string> innames,string oname){
  vector<volume<float> > tmpvec;
  tmpvec.reserve(innames.size());
  volume<float> tmp;
  cout<<"number of inputs "<<innames.size()<<endl;
  cout<<"Indices"<<endl;
  for(unsigned int i=0;i<innames.size();i++){
    cout<<i+1<<" "<<innames[i]<<endl;
    read_volume(tmp,innames[i]);
    tmpvec.push_back(tmp);
  }
  volume<int> output(tmp.xsize(),tmp.ysize(),tmp.zsize());
  copybasicproperties(tmp,output);output=0;

  for(int z=tmp.minz();z<=tmp.maxz();z++){
    for(int y=tmp.miny();y<=tmp.maxy();y++){
      for(int x=tmp.minx();x<=tmp.maxx();x++){
	RowVector bum(innames.size());
	Matrix bugger;
	ColumnVector index;
	for(unsigned int i=0;i<innames.size();i++ ){
	    bum(i+1)=tmpvec[i](x,y,z);
	} 
	bugger=max(bum,index);
	if(bum.MaximumAbsoluteValue()!=0)
	  output(x,y,z)=(int)index.AsScalar();
	else
	  output(x,y,z)=0;
      }  
    }
  }

  output.setDisplayMaximumMinimum(innames.size(),0);
  save_volume(output,oname);

}

void biggest_from_surfaces(vector<string> innames,string oname){
  vector<CsvMesh> meshes;
  CsvMesh m;int nv=0,nt=0;

  cout<<"number of inputs "<<innames.size()<<endl;
  cout<<"Indices"<<endl;
  for(unsigned int i=0;i<innames.size();i++){
    cout<<i+1<<" "<<innames[i]<<endl;
    nt++;
    m.load(innames[i]);
    meshes.push_back(m);
    // test size
    if(i==0)nv=m.nvertices();
    else{
      if(m.nvertices()!=nv){
	cerr<<"Number of vertices incompatible amongst input"<<endl;
	exit(1);
      }
    }
  }
      
  Matrix allvals(nv,nt);
  for(unsigned int i=0;i<meshes.size();i++){
    int cnt=1;
    for(int j = 0; j< meshes[i].nvertices();j++){
      allvals(cnt,i+1)=meshes[i].get_pvalue(j);
      cnt++;	
    }
  }

  for(int i=1;i<=nv;i++){
    if(allvals.Row(i).MaximumAbsoluteValue()==0){
      m.set_pvalue(i-1,0);
      continue;
    }
    int jmax=1;float vm=0;
    for(int j=1;j<=nt;j++){
      if(j==1)vm=allvals(i,j);
      if(allvals(i,j)>=vm){vm=allvals(i,j);jmax=j;}
    }
    m.set_pvalue(i-1,jmax);
  }

  m.save(oname,meshFileType(innames[0]));		      
  
}
bool test_input(const vector<string>& filenames){
  int nsurfs=0,nvols=0;
  for(unsigned int i=0;i<filenames.size();i++){
    if(meshExists(filenames[i])){nsurfs++;}
    else if(fsl_imageexists(filenames[i])){nvols++;}
    else{
      cerr<<"Cannot open "<<filenames[i]<<" for reading"<<endl;
      exit(1);
    }
  }
  if(nvols>0 && nsurfs>0){
    cerr<<"Please use list of EITHER volumes OR surfaces"<<endl;
    exit(1);
  }
  return (nvols>0?true:false);
}
int main ( int argc, char **argv ){
  if(argc<3){
    cerr<<" "<<endl;
    cerr<<"usage: find_the_biggest <lots of volumes/surfaces> output"<<endl;
    cerr<<"  output is index in order of inputs"<<endl;
    cerr<<"  please use EITHER volumes OR surfaces"<<endl;
    cerr<<" "<<endl;
    exit(1);
  }

  vector<string> innames;
  innames.clear();
  for(int i=1;i<=argc-2;i++){
    innames.push_back(argv[i]);
  }
  bool isvols=test_input(innames);
  if(!isvols){
    biggest_from_surfaces(innames,argv[argc-1]);
  }
  else{
    biggest_from_volumes(innames,argv[argc-1]);
  }

  return(0);

}









