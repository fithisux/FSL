/*  libpic - collection of image display and rendering routines

    Stephen Smith, Christian Beckmann and Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 1999-2009 University of Oxford  */

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


#ifndef __MISCPIC_h
#define __MISCPIC_h

#include "newimage/newimageall.h"
#include <stdarg.h>
#include "gd.h"

namespace MISCPIC{
  
  //template <class T>
  class miscpic
    {
    public:

      //constructor
      miscpic(){
	nlut = 0;      
	compare= 0;
	writeText=false;
	LR_label_flag = true;
	markRight=false;
	trans= -10;
	edgethresh = 0.0;
	if(getenv("FSLDIR")!=0){
	  lutbase = string(getenv("FSLDIR")) + "/etc/luts/";
	}
	else{
	  lutbase = string("/");
	}
	title = string("");
	cbartype = string("");
	cbarptr = NULL;
	outim = NULL;
	picr = NULL;
	picg = NULL;
	picb = NULL;
      };

      ~miscpic(){
	if(picr!=NULL) free(picr);
	if(picg!=NULL) free(picg);
	if(picb!=NULL) free(picb);
	if(cbarptr) gdImageDestroy(cbarptr);
	if(outim) gdImageDestroy(outim);
      }

      int slicer(const NEWIMAGE::volume<float>& vol1,const NEWIMAGE::volume<float>& vol2,const char *opts, bool labelSlices=false, bool debug = false);
      int slicer(const NEWIMAGE::volume<float>& vol1,const NEWIMAGE::volume<float>& vol2,vector<string> inputOptions, bool labelSlices=false, bool debug = false);

      inline int slicer(const NEWIMAGE::volume<float>& vol1,const char *opts, bool labelSlices=false, bool debug = false)
      { NEWIMAGE::volume<float> tmp(1,1,1); 
	return this->slicer(vol1, tmp, opts, debug);}
      
      int write_png ( char *filename, int x_size, int y_size,	 
 		      unsigned char *r, unsigned char *g, unsigned char *b); 
      int write_pgm ( char *filename, int x_size, int y_size,	
		      unsigned char *i );
      int write_ppm ( char *filename, int x_size, int y_size,	
		      unsigned char *r, unsigned char *g, unsigned char *b );
      void write_pic(char *fname, int width, int height);
      
      inline void set_title(string what) { title = what;}
      inline void set_cbar(string what) { cbartype = what;}

      inline void read_lut(char *fname){
	lut = string(fname); this->read_lut();
      }

      void read_lut();
      
      inline void set_minmax(float bgmin, float bgmax, float s1min,
		  float s1max){
	this->set_minmax(bgmin, bgmax, s1min, s1max, (float)0.0, float(0.0));
      };

      void set_minmax(float bgmin, float bgmax, float s1min,
		  float s1max, float s2min, float s2max);
      
      inline int overlay(NEWIMAGE::volume<float>& newvol, NEWIMAGE::volume<float>& bg, NEWIMAGE::volume<float>& s1,
			 NEWIMAGE::volume<float>& s2, float bgmin, float bgmax, float s1min,
			 float s1max, float s2min, float s2max, int colour_type,
			 int checker, 
			 bool out_int = false, bool debug = false){
	return this->overlay(newvol, bg, s1, s2, bgmin, bgmax, s1min, 
			     s1max, s2min, s2max,
			     colour_type, checker, 
			     string(""), string(""), out_int, debug);
      }

      int overlay(NEWIMAGE::volume<float>& newvol, NEWIMAGE::volume<float>& bg, NEWIMAGE::volume<float>& s1,
		  NEWIMAGE::volume<float>& s2, float bgmin, float bgmax, float s1min,
		  float s1max, float s2min, float s2max, int colour_type,
		  int checker, 
		  string cbarfname, string cbartype, bool out_int = false, 
		  bool debug = false);

    private:

      int x_size, y_size, z_size, size, x_size_pic, y_size_pic, 
	  z_size_pic, nlut, compare, trans;

      bool debug, LR_label_flag, writeText;

      float edgethresh;

      string lut, lutbase, title, cbartype;

      gdImagePtr cbarptr, outim;

      bool markRight;

      vector<int> rlut, glut, blut; //stores the lut;
      unsigned char *picr, *picg, *picb;

      //volume<T>  inp1, inp2, imr, img, imb;	  
      NEWIMAGE::volume<float>  inp1, inp2, imr, img, imb;
      vector<float> minmax; //will store min and max values for bg and stats images
      //needed for colorbar

      void sag(float xx, int p, int width);
      void cor(float yy, int p, int width);
      void axi(float zz, int p, int width);
            
      int create_cbar(string cbartype);
      int write_cbar(string fname, string cbartype);
      int add_cbar(string cbartype);

      int add_title(int width);
      void addRlabel(int p, int width, int size_pic, int alt_size_pic, 
		     bool onleft);
      void addRlabel(unsigned char* picr, int p, int width, int size_pic, 
		     int alt_size_pic, bool onleft);
  };
}

#endif
