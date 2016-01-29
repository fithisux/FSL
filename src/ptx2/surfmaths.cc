/*   surfmaths.cc 
     uniary and binary operations from fslmaths applied to surfaces

     Saad Jbabdi and Matthew Webster, FMRIB Image Analysis Group

     Copyright (C) 2012-2014 University of Oxford  */
     

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

     
#include "newmesh/newmesh.h"
#include "miscmaths/miscmaths.h"
#include "utils/fsl_isfinite.h"
#include "libprob.h"

using namespace MISCMATHS;
using namespace NEWMESH;


#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

int printUsage(const string& programName) 
{
  cout << "\nUsage: surfmaths <first_input> [operations and inputs] <output> " << endl;

  cout << "\nBasic unary operations:" << endl;
  cout << " -exp   : exponential" << endl;
  cout << " -log   : natural logarithm" << endl;
  cout << " -sin   : sine function" << endl;
  cout << " -cos   : cosine function" << endl;
  cout << " -tan   : tangent function" << endl;
  cout << " -asin  : arc sine function" << endl;
  cout << " -acos  : arc cosine function" << endl;
  cout << " -atan  : arc tangent function" << endl;
  cout << " -sqr   : square" << endl;
  cout << " -sqrt  : square root" << endl;
  cout << " -recip : reciprocal (1/current surface)" << endl;
  cout << " -abs   : absolute value" << endl;
  cout << " -bin   : use (current surface>0) to binarise" << endl;
  cout << " -nan   : replace NaNs (improper numbers) with 0" << endl;

  cout << "\nBinary operations:" << endl;
  cout << "  (some inputs can be either a surface or a number)" << endl;
  cout << " -add   : add following input to current surface" << endl;
  cout << " -sub   : subtract following input from current surface" << endl;
  cout << " -mul   : multiply current surface by following input" << endl;
  cout << " -div   : divide current surface by following input" << endl;
  cout << " -mas   : use (following surface>0) to mask current surface" << endl;
  cout << " -thr   : use following number to threshold current surface (zero anything below the number)" << endl;

  cout << "\ne.g. surfmaths inputSurface -add inputSurface2 outputSurface" << endl;
  cout << "     surfmaths inputSurface -add 2.5 outputSurface" << endl;
  cout << "     surfmaths inputSurface -add 2.5 -mul inputSurface2 outputSurface\n" << endl;

  return 1;
}

void loadSurface(const newmesh& iSurf,newmesh& tmpSurf,const string& fname){
  tmpSurf.load_gifti(fname,false);
  if(tmpSurf.npvalues()!=iSurf.npvalues()){
    cerr<<"Error: surfaces do not have the same number of vertices"<<endl;
    exit(1);
  }
}


int check_for_output_name(int i, int argc_1)
{
  if (i>argc_1) {
    cerr << "Error: no output filename specified!" << endl;
    exit(1);
  }
  return 0;
}


int inputParser(int argc, char *argv[]){
  newmesh inputSurface;
  inputSurface.load_gifti(string(argv[1]),false);
  int i=2;
  for (i = 2; i < argc-1; i++){    
    newmesh temp_surface;    
    
    /***************************************************************/
    /********************Binary Operations**************************/
    /***************************************************************/
    if (string(argv[i])=="-mas"){  
      loadSurface(inputSurface,temp_surface,string(argv[++i]));	
      for(int v=0;v<temp_surface.npvalues();v++)
	inputSurface.set_pvalue(v, temp_surface.get_pvalue(v)>0?inputSurface.get_pvalue(v):0 );
    }                                                           
    /***************************************************************/
    else if (string(argv[i])=="-add"){
      i++;
      if (isNumber(string(argv[i]))){
	for(int v=0;v<inputSurface.npvalues();v++)
	  inputSurface.set_pvalue(v, inputSurface.get_pvalue(v)+atof(argv[i]));
      }
      else{  
	loadSurface(inputSurface,temp_surface,string(argv[i]));	
	for(int v=0;v<temp_surface.npvalues();v++)
	  inputSurface.set_pvalue(v, temp_surface.get_pvalue(v)+inputSurface.get_pvalue(v) );
      }
    }
    /***************************************************************/
    else if (string(argv[i])=="-sub"){
      i++;
      if (isNumber(string(argv[i]))){
	for(int v=0;v<inputSurface.npvalues();v++)
	  inputSurface.set_pvalue(v, inputSurface.get_pvalue(v)-atof(argv[i]));
      }
      else{  
	loadSurface(inputSurface,temp_surface,string(argv[i]));	
	for(int v=0;v<temp_surface.npvalues();v++)
	  inputSurface.set_pvalue(v, temp_surface.get_pvalue(v)-inputSurface.get_pvalue(v) );
      }
    }
    /***************************************************************/
    else if (string(argv[i])=="-mul"){
      i++;
      if (isNumber(string(argv[i]))){
	for(int v=0;v<inputSurface.npvalues();v++)
	  inputSurface.set_pvalue(v, inputSurface.get_pvalue(v)*atof(argv[i]));
      }
      else{  
	loadSurface(inputSurface,temp_surface,string(argv[i]));	
	for(int v=0;v<temp_surface.npvalues();v++)
	  inputSurface.set_pvalue(v, temp_surface.get_pvalue(v)*inputSurface.get_pvalue(v) );
      }      
    }    
    /***************************************************************/
    else if (string(argv[i])=="-div"){
      i++;
      if (isNumber(string(argv[i]))){
	for(int v=0;v<inputSurface.npvalues();v++)
	  inputSurface.set_pvalue(v, inputSurface.get_pvalue(v)/atof(argv[i]));
      }
      else {  
	loadSurface(inputSurface,temp_surface,string(argv[i]));	
	for(int v=0;v<temp_surface.npvalues();v++)
	  inputSurface.set_pvalue(v, temp_surface.get_pvalue(v)/inputSurface.get_pvalue(v) );
      }   
    }    
    /***************************************************************/
    /******************** Unary Operations *************************/
    /***************************************************************/
    else if (string(argv[i])=="-thr"){
      for(int v=0;v<inputSurface.npvalues();v++){
	inputSurface.set_pvalue(v, inputSurface.get_pvalue(v)>=atof(argv[i+1])?inputSurface.get_pvalue(v):0);
      }
      i++;
    }    
    /***************************************************************/
    else if (string(argv[i])=="-sqrt"){
      for(int v=0;v<inputSurface.npvalues();v++)
	inputSurface.set_pvalue(v,inputSurface.get_pvalue(v)>=0?std::sqrt(inputSurface.get_pvalue(v)):0);
    }
    /***************************************************************/
    else if (string(argv[i])=="-pow"){
      for(int v=0;v<inputSurface.npvalues();v++)
	inputSurface.set_pvalue(v,pow(inputSurface.get_pvalue(v),atof(argv[i+1])));
      i++;
    }
    /***************************************************************/
    else if (string(argv[i])=="-sqr"){ 
    for(int v=0;v<inputSurface.npvalues();v++)
      inputSurface.set_pvalue(v,inputSurface.get_pvalue(v)*inputSurface.get_pvalue(v));
    }
    /***************************************************************/
    else if (string(argv[i])=="-recip"){
      for(int v=0;v<inputSurface.npvalues();v++)
	inputSurface.set_pvalue(v, inputSurface.get_pvalue(v)!=0?1/inputSurface.get_pvalue(v):0);
    }
    /***************************************************************/
    else if (string(argv[i])=="-exp"){
      for(int v=0;v<inputSurface.npvalues();v++)
	inputSurface.set_pvalue(v,std::exp((double)inputSurface.get_pvalue(v)));
    }
    /***************************************************************/
    else if (string(argv[i])=="-log"){
      for(int v=0;v<inputSurface.npvalues();v++)
	inputSurface.set_pvalue(v,inputSurface.get_pvalue(v)>0?std::log((double)inputSurface.get_pvalue(v)):0);
    }
    /***************************************************************/
    else if (string(argv[i])=="-cos"){
      for(int v=0;v<inputSurface.npvalues();v++)
	inputSurface.set_pvalue(v,std::cos((double)inputSurface.get_pvalue(v)));
    }
    /***************************************************************/
    else if (string(argv[i])=="-sin"){
      for(int v=0;v<inputSurface.npvalues();v++)
	inputSurface.set_pvalue(v,std::sin((double)inputSurface.get_pvalue(v)));
    }
    /***************************************************************/
    else if (string(argv[i])=="-tan"){
      for(int v=0;v<inputSurface.npvalues();v++)
	inputSurface.set_pvalue(v,std::tan((double)inputSurface.get_pvalue(v)));
    }
    /***************************************************************/
    else if (string(argv[i])=="-asin"){
      for(int v=0;v<inputSurface.npvalues();v++)
	inputSurface.set_pvalue(v,std::asin((double)inputSurface.get_pvalue(v)));
    }
    /***************************************************************/
    else if (string(argv[i])=="-acos"){
      for(int v=0;v<inputSurface.npvalues();v++)
	inputSurface.set_pvalue(v,std::acos((double)inputSurface.get_pvalue(v)));
    }
    /***************************************************************/
    else if (string(argv[i])=="-atan"){
      for(int v=0;v<inputSurface.npvalues();v++)
	inputSurface.set_pvalue(v,std::atan((double)inputSurface.get_pvalue(v)));
    }
    /***************************************************************/
    else if (string(argv[i])=="-abs"){      
      for(int v=0;v<inputSurface.npvalues();v++)
	inputSurface.set_pvalue(v,std::fabs(inputSurface.get_pvalue(v)));
    } 
    /***************************************************************/
    else if (string(argv[i])=="-bin"){
      for(int v=0;v<inputSurface.npvalues();v++)
	inputSurface.set_pvalue(v, inputSurface.get_pvalue(v)!=0?1:0);
    }
     /******************************************************/
    else if (string(argv[i])=="-nan"){
      for(int v=0;v<inputSurface.npvalues();v++)
	inputSurface.set_pvalue(v,!isfinite(inputSurface.get_pvalue(v))?0:inputSurface.get_pvalue(v));
     }
    else{
      cerr<<"unknown option "<<string(argv[i])<<endl;
      exit(1);
    }
  }
  inputSurface.save(string(argv[argc-1]));
  return 0;
}


int main(int argc,char *argv[]){
  if (argc < 2)  
    return printUsage(string(argv[0])); 
  
  if(string(argv[1]) == "-h" || string(argv[1]) == "--help") { 
    printUsage(string(argv[0])); 
    exit(0); 
  }
  return inputParser(argc,argv);
}



