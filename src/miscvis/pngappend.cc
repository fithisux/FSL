/*  pngappend - simple programm to append sets of png

    Christian Beckmann and Matthew Webster, FMRIB Image Analysis Group

    2007 University of Oxford  */

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

#include <string>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include "gd.h"
#include "miscmaths/miscmaths.h"

using namespace std;

void usage(void)
{
  cerr << endl << "Usage: pngappend <input 1> <+|-> [n] <input 2> [<+|-> [n] <input n>]  output>" << endl;

  cerr << endl << " + appends horizontally," << endl << " - appends vertically (i.e. works like a linebreak)" << endl;
  cerr << "[n] number ofgap pixels" <<  endl;
  cerr << "note that files with .gif extension will be input/output in GIF format" << endl; 
  exit(1);
}

void help()
{
  cerr<< endl << " pngappend  -  append PNG files horizontally and/or vertically into a new PNG (or GIF) file" << endl;
  usage();
}
    
int max (int a, int b)
{
  if (a<b) return b;
  else return a;
}

int main(int argc, char *argv[])
{
  if ((argc<2)||(strcmp(argv[1],"-help")==0)||(strcmp(argv[1],"--help")==0)) help();

  if ((argc<3)) usage();
  int argidx = 1, black, toprx = 0, topry = 0, gap = 0;
  bool appendh = true;

  gdImagePtr im_out, im_read, im_tmp;
  FILE *in;
  string filename;
  bool GIFOUT;
  // read in the first image 

  if((in = fopen(argv[argidx++], "rb"))==NULL){
    cerr << endl<< " Cannot open " << argv[argidx-1] << " for reading" << endl; 
    fclose(in);
    exit(1);
  }
  else if((im_read = gdImageCreateFromPng(in))==NULL)
  {
    rewind(in);
    if ((im_read = gdImageCreateFromGif(in))==NULL)
    { 
      cerr << endl <<argv[argidx-1] << " is not a valid png2 or gif file" << endl;
      fclose(in);
      exit(1);
    }
  };
  fclose(in);

  //copy it to im_out 

  im_out = gdImageCreateTrueColor(im_read->sx,im_read->sy); //dImageCreate rather than TrueColour for GIF - does it matter?
  gdImageCopy(im_out, im_read, 0, 0, 0, 0,im_read->sx, im_read->sy);
  gdImageDestroy(im_read);

  toprx = im_out->sx;
  topry = 0;

  while(argidx < argc-1){
    
    //copy im_out to im_tmp
    im_tmp = gdImageCreateTrueColor(im_out->sx,im_out->sy);
    gdImageCopy(im_tmp, im_out, 0, 0, 0, 0,im_out->sx, im_out->sy);
    gdImageDestroy(im_out);

    //read token
    if(strcmp(argv[argidx],"+")==0)
      appendh = true;
    else if(strcmp(argv[argidx],"-")==0)
      appendh = false;
    else{
      cerr << endl << "ERROR: use '+' or '-' to indicate horizontal or vertical concatenation" 
	   << endl;
       usage();
    }
    argidx++;

    //read in new image
    if(MISCMATHS::isNumber(argv[argidx])){
      gap = atoi(argv[argidx]);
      argidx++;
    }

    if((in = fopen(argv[argidx], "rb"))==NULL){
	cerr << endl<< " Cannot open " << argv[argidx] << " for reading" << endl; 
	fclose(in);
	exit(1);
    }
    else if((im_read = gdImageCreateFromPng(in))==NULL)
    {
      rewind(in);
      if ((im_read = gdImageCreateFromGif(in))==NULL)
      { 
        cerr << endl <<argv[argidx-1] << " is not a valid png or gif file" << endl;
        fclose(in);
        exit(1);
      }
    };
    fclose(in);

    //copy everything into img_out
    {
      int newx, newy;
  
      if(!appendh){
	toprx = 0;
	topry = im_tmp->sy + gap;
      }else
	toprx += gap;      

      newx = max(toprx + im_read->sx, im_tmp->sx);
      newy = max(topry + im_read->sy, im_tmp->sy);

      im_out = gdImageCreateTrueColor(newx,newy);
      black = gdImageColorAllocate(im_out, 0, 0, 0); 
      gdImageCopy(im_out, im_tmp, 0, 0, 0, 0,im_tmp->sx, im_tmp->sy);
      gdImageCopy(im_out, im_read, toprx, topry, 0, 0, im_read->sx, im_read->sy);       

      toprx += im_read->sx;
    }

    gdImageDestroy(im_read);
    gdImageDestroy(im_tmp);

    argidx++;
  }
	
  //output im_out
  FILE *imageout;
  filename=string(argv[argc-1]);
  GIFOUT=(filename.substr(filename.size()-4,filename.size())==string(".gif"));
  if((imageout = fopen(argv[argc-1], "wb"))==NULL){
    cerr << endl <<"ERROR: cannot open " << argv[argc-1] << "for writing" << endl;
    exit(1);
  }
  if (GIFOUT)
  {
    im_tmp=gdImageCreate(im_out->sx,im_out->sy);
    gdImageCopy(im_tmp, im_out, 0, 0, 0, 0,im_out->sx, im_out->sy);
    gdImageGif(im_tmp, imageout);
    gdImageDestroy(im_tmp);
  }
  else gdImagePng(im_out, imageout);
  gdImageDestroy(im_out);
  fclose(imageout);

  return 0;
}
