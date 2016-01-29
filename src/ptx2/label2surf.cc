/*  label2surf.cc

    Transform a set of label files into a surface

    Saad Jbabdi, FMRIB Image Analysis Group

    Copyright (C) 2010 University of Oxford */

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


#include "utils/options.h"
#include "csv_mesh.h"

using namespace Utilities;


string title="label2surf \n\t Transform a group of labels into a surface";
string examples="label2surf -s <surface> -o <outputsurface> -l <labels>";

Option<bool> verbose(string("-v,--verbose"),false,
		       string("switch on diagnostic messages"),
		       false,no_argument);
Option<bool> help(string("-h,--help"),false,
		       string("display this message"),
		       false,no_argument);
Option<string> isurf(string("-s,--surf"),"",
		     string("input surface"),
		     true,requires_argument);
Option<string> osurf(string("-o,--out"),"",
		     string("output surface"),
		    true,requires_argument);
Option<string> labels(string("-l,--labels"),"",
		      string("ascii list of label files"),
		      true,requires_argument);


void read_fnames(vector<string>& fnames,const string& filename){
  fnames.clear();
  ifstream fs(filename.c_str());
  string tmp;
  if(fs){
    fs>>tmp;
    do{
      fnames.push_back(tmp);
      fs>>tmp;
    }while(!fs.eof());
  }
  else{
    cerr<<filename<<" does not exist"<<endl;
    exit(0);
  }
}

void read_label(vector<int>& IDs,const string& labelfile){
  IDs.clear();
  ifstream fs(labelfile.c_str());
  if (!fs) { 
    cerr << "Could not open label file " << labelfile << endl;
    exit(1);
  }
  // read first line
  char str[200];
  FILE *fp;
  string firstline;
  fp = fopen(labelfile.c_str(), "r");
  fscanf(fp, "%[^\n]", str);
  firstline = str;
  string cline;
  // skip header
  cline = skip_alpha(fs);
  // read number of vertices
  string ss="";
  fs >> ss;
  float nrows = atof(ss.c_str());
  for(int r=1;r<=nrows;r++){
    for(int c=1;c<=5;c++){
      if(!fs.eof()){
	fs >> ss;
	while ( !isNumber(ss) && !fs.eof() ) {
	  fs >> ss;
	}
	if(c==1){IDs.push_back(atoi(ss.c_str()));}
      }
    }
  }
}


int main(int argc,char *argv[]){
  OptionParser options(title,examples);
 
  options.add(verbose);
  options.add(help);
  options.add(isurf);
  options.add(osurf);
  options.add(labels);

  options.parse_command_line(argc,argv);
  if ( (help.value()) || (!options.check_compulsory_arguments(true)) ){
    options.usage();
    return 1;
  }
  
  ////////
  if(verbose.value())
    cout<<"read input surface"<<endl;
  
  CsvMesh m;
  m.load(isurf.value());
  m.reset_pvalues();
  m.reset_tvalues();
  
  if(verbose.value())
    cout<<"read input labels"<<endl;
  
  vector<string> labs;
  vector<int>    IDs;
  read_fnames(labs,labels.value());
  for(unsigned int i=0;i<labs.size();i++){
    if(verbose.value())
      cout<<"   label " << i+1 << endl;
    read_label(IDs,labs[i]);
    m.set_pvalues(IDs,1);
  }

  m.save(osurf.value());

  return 0;
}
