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
#include "newimage/newimageall.h"
#include <vector>
using namespace std;
using namespace NEWIMAGE;
int read_avg_file (vector<vector<int> >& avgs,const string fname){
  avgs.clear();
  ifstream avg_file(fname.c_str());
  string myline;
  bool nobbsize=true;

  int row = 0;
  
  if(!avg_file){return -1;}
  else{
    while(nobbsize){
      avgs.push_back(vector<int>());

      nobbsize=false;      
      getline(avg_file,myline);
      
      int pos=0;
      while(pos!=int(string::npos)) {
	pos = myline.find(",",pos);
	if(pos!=int(string::npos)){
	  myline.replace(pos,1," ");
	  pos++;
	}
      }

      istringstream mylinestr(myline.c_str());   

      while(!mylinestr.eof()){

	string startstr;
	int start;
	int length;
	mylinestr >> startstr;
	if(isNumber(startstr)){
	  nobbsize=true;
	  start = atoi(startstr.c_str());
	  mylinestr >> length;

	  //cerr<< start <<" " << length << " ";

	  for(int i=start; i<start+length; i++)
	    {
	      avgs[row].push_back(i);
	    }
	}
      }
      row++;
      cerr<<endl;
    }
  }


  row--;
  avgs.pop_back();
 
//   for(int r=0;r<row;r++)
//     {
//       cout << "size=" << avgs[r].size() << endl;
//       for(int c=0;c<avgs[r].size();c++)
// 	cout << avgs[r][c] << " ";
//       cout << endl;
//     }
 
  return 0;
}
 

int main ( int argc, char **argv ){
  if(argc<5){
    cerr<<"usage: replacevols input4D avg_file volno volno ... volno output4D"<<endl;
    cerr<<"volno starts at zero!!!"<<endl;
    exit(0);
  }

  volume4D<float> data4D;
  read_volume4D(data4D,argv[1]);
  vector<vector<int> > avgs;
  int ret;
  ret = read_avg_file(avgs, argv[2]);
  if(ret==-1){
    cerr<<"can't find "<<argv[2]<<endl;
    return -1;
  }

  vector<int> repvols;
  for(int i=3;i<=argc-2;i++){
    repvols.push_back(atoi(argv[i]));
  }


  for(unsigned int j=0;j<avgs[0].size();j++){//loop over volume numbers


    //Next loop is within volume number over averages just 
    // Working out which ones to replace and which to keep.
    
    vector<int> repthis,keepthis;
    for(unsigned int i=0;i<avgs.size();i++){ //loop over averages
      bool replaced=false;
      for(unsigned int r=0;r<repvols.size();r++){// loop over things to be replaced
	if(avgs[i][j]==repvols[r]){
	  replaced=true;
	  repthis.push_back(avgs[i][j]);
	}
      }
      if(!replaced){
	keepthis.push_back(avgs[i][j]);
      }

    }
      

    if(repthis.size()>0){
      
      cerr<<"Replacing volumes: ";
      for(unsigned int r=0;r<repthis.size();r++){  
	cerr<<repthis[r]<<" ";
      }
      cerr <<endl;
      cerr<<"with the average of volumes: ";
      for(unsigned int r=0;r<keepthis.size();r++){  
	cerr<<keepthis[r]<<" ";
      }
      cerr<<endl;
      

      if(keepthis.size()>0){
	// Next loop makes the average of all the ones we are keeping 
	volume<float> tmp;
	tmp=data4D[keepthis[0] ];
	
	for(unsigned int n=1;n<keepthis.size();n++){
	  tmp=tmp+data4D[keepthis[n] ]; 
	}
	tmp=tmp/keepthis.size(); //Average of all non-replaced ones.
	
	
	
	//Next loop replaces all the ones to be replaced with this average
	for(unsigned int n=0;n<repthis.size();n++){
	  data4D[repthis[n] ]=tmp; //replacing.
	}
	
	
      }
      else{
	cerr<<"Error: Volume number "<<j<<" has no averages to keep!!"<<endl;;
	return -1; 
      }
    }//repthis.size>0
    
    
  }// loop over volume numbers




  //   vector<volume<float> > tmpvec;
//   float test=(argc-3)/2.0f;
//   if(round(test)!=test){
//     cerr<<"Either you have different nums of volnos and vols, or you haven't specified input or output"<<endl;
//     return(0);
//   }
    
//   //  int numchanges=(argc-3)/2;
//   tmpvec.reserve((argc-3)/2);
//   vector<int> volnos;
//   volume<float> tmp;
//   tmpvec.reserve((argc-3)/2);
  
//   cout<<"number of vols to be replaced "<<(argc-3)/2<<endl;
//   for(int i=2;i<=(argc-2);i+=2){
//     volnos.push_back(atoi(argv[i]));
//     read_volume(tmp,argv[i+1]);
//     tmpvec.push_back(tmp);
//   }
  
//   for(int i=0;i<(int)volnos.size();i++){
//     int num=volnos[i];
//     data4D[num]=tmpvec[i];
//   }
    
  save_volume4D(data4D,argv[argc-1]);
}









