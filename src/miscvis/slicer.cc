/*  slicer

    Christian Beckmann and Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 2006-2009 University of Oxford  */

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
 
#include "libvis/miscpic.h"

using namespace NEWIMAGE;
using namespace MISCPIC;

void usage(void)
{
  printf("\nUsage: slicer <input> [input2] [main options] [output options - any number of these]\n\n");

  printf("Main options: [-L] [-l <lut>] [-s <scale>] [-i <intensitymin> <intensitymax>] [-e <thr>] [-t] [-n] [-u]\n");
  printf("These must be before output options.\n");
  printf("-L       : Label slices with slice number.\n");
  printf("-l <lut> : use a different colour map from that specified in the header.\n");
  printf("-i <min> <max> : specify intensity min and max for display range.\n");
  printf("-e <thr> : use the specified threshold for edges (if >0 use this proportion of max-min, if <0, use the absolute value) \n");
  printf("-t       : produce semi-transparent (dithered) edges.\n");
  printf("-n       : use nearest-neighbour interpolation for output.\n");
  printf("-u       : do not put left-right labels in output.\n\n");
  printf("-c       : add a red dot marker to top right of image");

  printf("Output options:\n");
  printf("[-x/y/z <slice> <filename>]      : output sagittal, coronal or axial slice\n     (if <slice> >0 it is a fraction of image dimension, if <0, it is an absolute slice number)\n");
  printf("[-a <filename>]                  : output mid-sagittal, -coronal and -axial slices into one image\n");
  printf("[-A <width> <filename>]          : output _all_ axial slices into one image of _max_ width <width>\n");
  printf("[-S <sample> <width> <filename>] : as -A but only include every <sample>'th slice\n\n");

  exit(1);
}

int fmrib_main(int argc, char* argv[])
{
  //Option parsing
  string allCommands("");
  for(int token=1;token<argc;token++)
    allCommands+=string(argv[token])+" ";
  char* remainingChars = new char[allCommands.size()+1];
  strcpy(remainingChars,allCommands.c_str());

  char* tokenised;
  vector<string> tokens, miscpicOptions, nonOptionInputs;
  tokenised=strtok(remainingChars," ");
  do {
    tokens.push_back(string(tokenised));   
  } while ( (tokenised=strtok(NULL," ")) != NULL );

  for (unsigned int currentToken=0;currentToken<tokens.size();)
  {
    string currentOption;
    if ( tokens[currentToken].compare(0,1,"-")==0 ) { //Start of an option
      string optionLetter=tokens[currentToken].substr(1,1);
      unsigned int requiredArgs(0);
      switch (optionLetter.c_str()[0]) {
        case 'a' :
        case 'l' :
        case 's' :
        case 'e' :
	  requiredArgs=1;
	  break;
        case 'A' :
        case 'i' :
        case 'x' :
        case 'y' :
        case 'z' :
	  requiredArgs=2;
	  break;
        case 'S' :
	  requiredArgs=3;
	  break;
      }
      for (unsigned int requestedToken=0; requestedToken <= requiredArgs && currentToken<tokens.size(); requestedToken++) 
	currentOption+=tokens[currentToken++]+" ";
      miscpicOptions.push_back(currentOption);
    }
    else nonOptionInputs.push_back(tokens[currentToken++]);
  }
  //End of parsing

  volume<float> inputVolume, secondaryVolume(1,1,1); 
  try
  {
    read_volume(inputVolume,nonOptionInputs[0]);
  }
  catch (...)
  {
    cout << "Error in slicer input, exiting..." << endl;
    exit(1);
  }
  //Fix -x option for neurological images
  for (unsigned int option=0;option<miscpicOptions.size();option++)
    {
    if ( miscpicOptions[option].compare(0,2,"-x") == 0 ) {
      string subString(miscpicOptions[option].substr(3,string::npos));
      string filename(subString);
      filename.erase(0,filename.find(' ')+1);
      subString.erase(subString.find(' '),string::npos); //This substring contains the slice number ( in nifti format )
      if ( subString[0] == '-' ) {
	int sliceNumber=-atof(subString.c_str());
	ColumnVector v0(4);
	v0 << sliceNumber << 0 << 0 << 1.0;
	v0 = inputVolume.niftivox2newimagevox_mat() * v0;
        sliceNumber=-MISCMATHS::round(v0(1));
	miscpicOptions[option]="-x "+num2str(sliceNumber)+" "+filename;
      } else {
	double sliceNumber=atof(subString.c_str())*inputVolume.maxx();
	ColumnVector v0(4);
	v0 << sliceNumber << 0 << 0 << 1.0;
	v0 = inputVolume.niftivox2newimagevox_mat() * v0;
        sliceNumber=v0(1)/(float)inputVolume.maxx();
	sliceNumber=std::max(std::min(sliceNumber,1.0),0.0);
	miscpicOptions[option]="-x "+num2str(sliceNumber)+" "+filename;
      }
    }
  }

  if ( nonOptionInputs.size()>1 && FslFileExists(nonOptionInputs[1].c_str()) ) 
    read_volume(secondaryVolume,nonOptionInputs[1]);

  bool printDebug(false),labelSlices(false);
  allCommands="";
  for (unsigned int option=0;option<miscpicOptions.size();option++)
  {
    if ( miscpicOptions[option]=="-d " )
      printDebug=true;
    else if ( miscpicOptions[option]=="-L " )
      labelSlices=true;
    else allCommands+=miscpicOptions[option]+" ";
  }

  if ( printDebug ) {
    for (unsigned int option=0;option<miscpicOptions.size();option++)
      cerr << "Option " << option << ": " << miscpicOptions[option] << endl;
    cerr << allCommands.c_str() << endl;
  }

  delete [] remainingChars;
  miscpic newpic;
  return newpic.slicer(inputVolume, secondaryVolume, miscpicOptions, labelSlices, printDebug);
}

int main(int argc,char *argv[])
{
  if (argc<2) 
    usage();
  return fmrib_main(argc,argv); 
}
