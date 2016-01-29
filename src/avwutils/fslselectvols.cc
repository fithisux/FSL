/*  fslselectvols.cc

    Saad Jbabdi, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford */

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

#include "miscmaths/miscmaths.h"
#include "newimage/newimageall.h"
#include "utils/options.h"

using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;


int parse_filterstring(vector<int>& comps,const string& vols){
  int ctr=0;    
  char *p;
  char t[1024];
  const char *discard = ", [];{(})abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ~!@#$%^&*_-=+|\':><./?";
  comps.clear();

  strcpy(t, vols.c_str());
  p=strtok(t,discard);
  ctr = atoi(p);
  if(ctr>=0)
    comps.push_back(ctr);
  do{
    p=strtok(NULL,discard);
    if(p){
      ctr = atoi(p);
      if(ctr>0){
	comps.push_back(ctr);
      }
    }
  }while(p);
  return 0;
}
bool file_exists (const std::string& name) {
    ifstream f(name.c_str());
    if (f.good()) {
        f.close();
        return true;
    } else {
        f.close();
        return false;
    }   
}

int get_vols(Matrix& id_vols,const string& vols){
  vector<int> comps;
  if(file_exists(vols)){
    id_vols = read_ascii_matrix(vols);
    if(id_vols.Ncols()==1){id_vols=id_vols.t();}

    comps.resize(id_vols.Ncols());
    for(int i=0;i<(int)comps.size();i++){
      comps[i]=id_vols(1,i+1);
    }
  }
  else{
    if(vols.length()>0 && parse_filterstring(comps,vols)){
      return 1;		
    }
  }

  //sort and remove duplicates 
  sort (comps.begin(), comps.end());
  vector<int>::iterator it = unique (comps.begin(), comps.end());
  comps.resize( it - comps.begin() );
    
  id_vols.ReSize(1,comps.size());
  for(int i=0;i<(int)comps.size();i++){
    id_vols(1,i+1)=comps[i];
  }
  
  
  return 0;
}



string title=string("fslselectvols\n")+
  string(" Select volumes from a 4D time series and output a subset 4D volume\n");
string examples=string("fslselectvols -i <input> -o <output> --vols=0,1,6,7\n")+
  string("fslselectvols -i <input> -o <output> --vols=volumes_list.txt");

Option<bool> help(string("-h,--help"), false,
		string("display this help text"),
		false,no_argument);
Option<string> ifilename(string("-i,--in"), string(""),
		string("\tinput file name (4D image)"),
		true, requires_argument);
Option<string> ofilename(string("-o,--out"), string(""),
		string("output file name (4D image)"),
		true, requires_argument);
Option<string> vols(string("--vols"), string(""),
		string("\tlist of volumes to extract (comma-separated list or ascii file)"),
		true, requires_argument);
Option<bool>   calc_mean(string("-m"), false,
		string("\toutput mean instead of concat"),
		false, no_argument);
Option<bool>   calc_var(string("-v"), false,
		string("\toutput variance instead of concat"),
		false, no_argument);


int main(int argc,char *argv[]){
  OptionParser options(title, examples);
  try{
    options.add(help);
    options.add(ifilename);
    options.add(ofilename);
    options.add(vols);
    options.add(calc_mean);
    options.add(calc_var);
    options.parse_command_line(argc, argv);
    if ( (help.value()) || (!options.check_compulsory_arguments(true)) ){
      options.usage();
      exit(EXIT_FAILURE);
    }
    else{
      
      volume4D<float> data;
      volume4D<float> subdata;
      read_volume4D(data,ifilename.value());

      Matrix list;
      get_vols(list,vols.value());

      bool meancalc=calc_mean.value(), varcalc=calc_var.value();

      if(!meancalc && !varcalc)
	subdata.reinitialize(data.xsize(),data.ysize(),data.zsize(),list.Ncols());
      else
	subdata.reinitialize(data.xsize(),data.ysize(),data.zsize(),1);
      subdata=0;
      copybasicproperties(data,subdata);
      
      float v,v2;
      for(int z=0;z<data.zsize();z++){
	for(int y=0;y<data.ysize();y++){
	  for(int x=0;x<data.xsize();x++){
	    v=0;v2=0;
	    for(int i=1;i<=list.Ncols();i++){	  
	      v += data[list(1,i)](x,y,z);
	      v2 += data[list(1,i)](x,y,z)*data[list(1,i)](x,y,z);
	      if(!meancalc && !varcalc)
		subdata[i-1](x,y,z)= data[list(1,i)](x,y,z);
	    }
	    if(meancalc)
	      subdata(x,y,z,0)= v/float(list.Ncols());
	    if(varcalc)
	      subdata(x,y,z,0)= v2/float(list.Ncols()) -v*v/float(list.Ncols()*list.Ncols());
	  }
	}
      }
      
      save_volume4D(subdata,ofilename.value());
      return 0;

    }
  }catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  }catch(std::exception &e) {
    cerr << e.what() << endl;
  } 
  // either mean or variance but not both
  if(calc_mean.value() && calc_var.value()){
    cerr<<" Cannot output mean *and* variance. Please choose one or the other but not both"<<endl;
    exit(1);
  }
  
  
}
