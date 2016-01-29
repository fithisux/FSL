/*  miscplot.cc  -- miscellaneous plotting functions

    Christian F. Beckmann and Matthew Webster, FMRIB Image Analysis Group
    
    Copyright (C) 1999-2007 University of Oxford */

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

#include <sstream>
#include <stdio.h>
#include <math.h>

#include "miscplot.h"
#include "miscmaths/miscprob.h"
#include "gdfonts.h"
#include "utils/fsl_isfinite.h"

extern "C" {
 #include "gd.h"
 #include "gdc.h"
 #include "gdchart.h"
 #include "gdcpie.h"
}

namespace MISCPLOT{
	//  MATLAB Line colors
	unsigned long sc_init[64] = {
		0x8080FF,0x007F00,0xFF0000,0x00BFBF,0xBF00BF,0xBFBF00,0x3F3F3F,
		0x0000FF,0x007F00,0xFF0000,0x00BFBF,0xBF00BF,0xBFBF00,0x3F3F3F,
		0x0000FF,0x007F00,0xFF0000,0x00BFBF,0xBF00BF,0xBFBF00,0x3F3F3F,
		0x0000FF,0x007F00,0xFF0000,0x00BFBF,0xBF00BF,0xBFBF00,0x3F3F3F,
		0x0000FF,0x007F00,0xFF0000,0x00BFBF,0xBF00BF,0xBFBF00,0x3F3F3F,
		0x0000FF,0x007F00,0xFF0000,0x00BFBF,0xBF00BF,0xBFBF00,0x3F3F3F,
		0x0000FF,0x007F00,0xFF0000,0x00BFBF,0xBF00BF,0xBFBF00,0x3F3F3F,
		0x0000FF,0x007F00,0xFF0000,0x00BFBF,0xBF00BF,0xBFBF00,0x3F3F3F,
		0x0000FF,0x007F00,0xFF0000,0x00BFBF,0xBF00BF,0xBFBF00,0x3F3F3F,0x0000FF};

	using namespace NEWMAT;
	using namespace MISCMATHS;

	string int2str(int n){
    ostringstream os;
    os.setf(ios::internal, ios::adjustfield);
    os << n;
    return os.str();
	}

string float2str(float f, int width, int prec, bool scientif){
    ostringstream os;
    int redw = int(std::abs(std::log10(std::abs(f))))+1;
    if(width>0)
      os.width(width);
    if(scientif)
      os.setf(ios::scientific);
    os.precision(redw+std::abs(prec));
    os.setf(ios::internal, ios::adjustfield);
    os << f;
    return os.str();
  }

void miscplot::add_legend(void* ptr, unsigned long cmap[], bool inside){
  
  int xsize = gdImagePtr(ptr)->sx;
  int ysize = gdImagePtr(ptr)->sy;
  int dest_x = 0;

  int linebrk = 5;
  int space = gdFontSmall->w + 2;
  int linelength = 15;

  int rlbllength = 0;
  for(int ctr= 0; ctr < int(labels.size()); ctr++)
    rlbllength = max(rlbllength,int(labels[ctr].length()));

  if(explabel.length()>0)
    ysize += space + 0*gdFontSmall->h;

  if(xlabels.size()>0)
    ysize += xlabels.size()*(linebrk + gdFontSmall->h)+linebrk;

  if(ylabels.size()>0)
    dest_x += ylabels.size()*(linebrk + gdFontSmall->h)+2*linebrk;

  if((!inside)&&(labels.size()>0))
    xsize += 2*space + linelength + rlbllength * gdFontSmall->w;
  
  xsize += dest_x;
  gdImagePtr newim;
  newim = gdImageCreate(xsize,ysize);
  gdImageCopy(newim,gdImagePtr(ptr),dest_x,0,0,0,gdImagePtr(ptr)->sx,gdImagePtr(ptr)->sy);

  int black = gdImageColorResolve(newim, 0, 0, 0);
  int ycoor, xcoor, yoffset=0;

  //make scale info
  if(explabel.length()>0){
    yoffset=space;

    xcoor=gdImagePtr(ptr)->sx - (4 + explabel.length())*gdFontSmall->w;
    ycoor=gdImagePtr(ptr)->sy+linebrk;
    char *s = (char*)string("x10").c_str();
    gdImageString(newim, gdFontSmall,xcoor,ycoor, (unsigned char*) s, black);
    xcoor+=(int)3.5*gdFontSmall->w;
    ycoor-=gdFontSmall->h/2;
    string s2 = string("-")+explabel;
    char *s3 = (char*)s2.c_str();
    gdImageString(newim, gdFontSmall,xcoor,ycoor, (unsigned char*) s3, black);
  }

  // make xlabel
  ycoor=gdImagePtr(ptr)->sy+linebrk+yoffset;
  for(int ctr=0; ctr < int(xlabels.size()); ctr++){
     xcoor = dest_x + gdImagePtr(ptr)->sx/2 - ((int)xlabels[ctr].length())/2*gdFontSmall->w;
     char *s = (char*)xlabels[ctr].c_str();
     gdImageString(newim, gdFontSmall,xcoor,ycoor, (unsigned char*) s, black);
     ycoor+=gdFontSmall->h + linebrk;
  }

  // make legend
  ycoor = 2+2*gdFontSmall->h;
  for(int ctr=0; ctr < int(labels.size()); ctr++){
    if (labels[ctr]!=""){
    	xcoor = dest_x+gdImagePtr(ptr)->sx;
    	if(inside)
      	xcoor -= 2*space + linelength + rlbllength * gdFontSmall->w;

    	int col =  gdImageColorResolve(newim, ((cmap[ctr]&PVRED)>>16),((cmap[ctr]&PVGRN)>>8), (cmap[ctr]&0x000000FF));

    	gdImageLine(newim, xcoor, ycoor + gdFontSmall->h/2, xcoor + linelength, ycoor + gdFontSmall->h/2,col);
    	gdImageLine(newim, xcoor, ycoor + gdFontSmall->h/2+1, xcoor + linelength, ycoor + gdFontSmall->h/2+1,col);
    	char *s = (char*)labels[ctr].c_str();
    	gdImageString(newim, gdFontSmall, xcoor + linelength + space, ycoor, (unsigned char*) s, black);
    	ycoor += gdFontSmall->h + linebrk;
		}
  }

  // make ylabel
  xcoor=space;
  for(int ctr=0; ctr < int(ylabels.size()); ctr++){ 
    int ycoor = 3*gdImagePtr(ptr)->sy/5 + ((int)ylabels[ctr].length())/2*gdFontSmall->w;
     char *s = (char*)ylabels[ctr].c_str();
     gdImageStringUp(newim, gdFontSmall,xcoor,ycoor, (unsigned char*) s, black);
     xcoor+=gdFontSmall->h + linebrk;
  }

  outim = newim;
}

void miscplot::timeseries(const Matrix& mat, string filename, string title, 
	  float tr, int ysize, int width, int prec, bool sci){

  int numlines=mat.Nrows();
  int numpoint=mat.Ncols();
  float* data = new float[numlines*numpoint];
 
  GDC_interpolations=TRUE;
  for(int ctr1=1;ctr1<=numlines; ctr1++)
      for(int ctr2=1;ctr2<=numpoint; ctr2++)
        if (!isfinite(data[(ctr1-1)*numpoint + ctr2-1])) 
 	   			data[(ctr1-1)*numpoint + ctr2-1] = GDC_INTERP_VALUE;
  if(ymin < ymax){
    GDC_requested_ymax = ymax;
    GDC_requested_ymin = ymin;
    //GDC_hard_size = TRUE;  //DO NOT uncomment this - it seems to causes a memory leak
    //possibly in conjunction with GDC_EXPOSE_IMAGE that will cause tsplot to stop working
    for(int ctr1=1;ctr1<=numlines; ctr1++){
      for(int ctr2=1;ctr2<=numpoint; ctr2++){ 
	    	data[(ctr1-1)*numpoint + ctr2-1] = std::min(std::max(float(mat(ctr1,ctr2)),ymin),ymax);
      }
    }
  }
	else{
    GDC_requested_ymax = minmaxscale*mat.Maximum();
    GDC_requested_ymin = (2-minmaxscale)*mat.Minimum();
    for(int ctr1=1;ctr1<=numlines; ctr1++){
      for(int ctr2=1;ctr2<=numpoint; ctr2++){
				data[(ctr1-1)*numpoint + ctr2-1] = mat(ctr1,ctr2);
      }
    }
		ymin = GDC_requested_ymin;
		ymax = GDC_requested_ymax;
  }

  for(int ctr1=1;ctr1<=numlines; ctr1++)
    for(int ctr2=1;ctr2<=numpoint; ctr2++)
      if (!isfinite(data[(ctr1-1)*numpoint + ctr2-1])) data[(ctr1-1)*numpoint + ctr2-1]=GDC_INTERP_VALUE;

  if (tr<0) tr=-tr;

  spacing =1; 
  if(numpoint>20) spacing = 5;
  if(numpoint>40) spacing = 10;
  if(numpoint>100) spacing = 20;
  if(numpoint>500) spacing = 50;
  if(numpoint>1200) spacing = 100;  
  
	int xsize = std::max(4*numpoint,750);
  if(xsize>1200) xsize=(int)floor(std::max(xsize/2.0,750.0));

  if(req_xsize>0){
    xsize=req_xsize;
    ysize=req_ysize;
  }

  string s;
  s=float2str((1)*tr,width,prec,sci);
  int maxlabellength = s.length()+1;
  s=float2str((numpoint-1) * tr,width,prec,sci);
  maxlabellength = std::max(maxlabellength,(int)s.length()+1);
 
  while ( float(xsize) * float(spacing) / float(numpoint) < 8 * maxlabellength){
    spacing = 2 * spacing;
  }
   
  char* ctl = new char[numpoint];

  if(((tr == 0.0)&&(spacing%2 == 1)) || (tr > 0.0))
    ctl[0] = TRUE;
  else
    ctl[0] = FALSE;
  for(int ctr=1;ctr<numpoint;ctr++)
    ctl[ctr]=((((ctr+1)%spacing == 0)&&(tr == 0.0))
	      ||(((ctr)%spacing == 0)&&(tr != 0.0)));

  int maxctr=0;
  for(int ctr=1;ctr<numpoint;ctr++)
    if (ctl[ctr]==true) maxctr=ctr;
  char** lbls = new char*[numpoint];
  
  lbls[0] = new char[2];
  if(tr == 0.0)
    strcpy(lbls[0],string("1").c_str());
  else
    strcpy(lbls[0],string("0").c_str());
  
  for(int ctr=1;ctr<numpoint;ctr++){
    if(ctl[ctr]){
      if(tr == 0.0)
				s = int2str(ctr+1);
      else
				s = float2str((ctr)*tr,width,prec,sci);
      if (ctr==maxctr ) //Padding for final label to avoid it being "chopped"
      { 
       int len=numpoint-maxctr;
       int len2=(s.length())-len;
       //cerr << len << " " << s << " " << s.length() << " " << len2 << endl;
        while(len2>0) { s+=" "; len2-- ; }
      }

      lbls[ctr] = new char[s.length()+1];
      strcpy(lbls[ctr],s.c_str());
    }
    else{
      lbls[ctr] = new char[1];
      strcpy(lbls[ctr],string("").c_str());
    }
  }
	
  {
    int tmp = _gdccfoo1;
    tmp++;
    long unsigned int tmp2 = _gdccfoo2;
    tmp2++;
  }
	
  GDC_xlabel_ctl = ctl;

  GDC_BGColor   = 0xFFFFFFL;                  /* backgound color (white) */
  GDC_LineColor = 0x000000L;                  /* line color      (black) */
  GDC_SetColor  = &(sc[0]);                   /* MATLAB-like colors */

  GDC_title = (char*)title.c_str();

  GDC_ticks = GDC_TICK_LABELS;
  GDC_grid  = GDC_TICK_NONE;
	if(gridswapdefault)
		GDC_grid = GDC_TICK_LABELS;
  GDC_yaxis = TRUE;
  GDC_yaxis2 = FALSE;
  GDC_xaxis = TRUE;
  GDC_xaxis_angle = 0;
  float range;
  range=abs(ymax-ymin);
  if (range<1) GDC_ylabel_density = 75;
  else GDC_ylabel_density = 70;	 

  if      (range>1e6)  GDC_ylabel_fmt = (char*)"%.1e";
  else if (range>10)   GDC_ylabel_fmt = (char*)"%.0f";
  else if (range>1)    GDC_ylabel_fmt = (char*)"%.1f";
  else if (range>0.01) GDC_ylabel_fmt = (char*)"%.3f";
  else                 GDC_ylabel_fmt = (char*)"%.5f";

	if(Ylabel_fmt > "")
		GDC_ylabel_fmt = (char*)Ylabel_fmt.c_str();

  GDC_title_size = GDC_SMALL;
  //GDC_hard_xorig = 50;
  //GDC_hard_graphwidth = xsize - 65;
  if(filename.substr(filename.size()-4,filename.size())!=string(".png"))
    filename += string(".png");
  
  FILE  *outpng1 = fopen(filename.c_str(), "wb" );
  GDC_image_type     = GDC_PNG;
  gdImagePtr outim2=NULL;

  //if (labels.size()>0||ylabels.size()>0||xlabels.size()>0)    GDC_hold_img = GDC_EXPOSE_IMAGE;
  GDC_out_graph( xsize, ysize, outpng1 , GDC_LINE, numpoint, (char**) lbls , numlines, data, NULL ); 
  fclose (outpng1 );
  if (labels.size()>0||ylabels.size()>0||xlabels.size()>0) 
  {
    outpng1=fopen(filename.c_str(), "rb" );
    outim2=gdImageCreateFromPng(outpng1);
    fclose( outpng1 );
  }
  if (labels.size()>0||ylabels.size()>0||xlabels.size()>0){
    outpng1 = fopen(filename.c_str(), "wb" );
    add_legend(outim2, sc);
    gdImagePng(outim, outpng1 );
    gdImageDestroy(outim2);
    fclose( outpng1 ); 
    if(outim) gdImageDestroy(outim);
  }
  for(int ctr=0;ctr<numpoint;ctr++){
    delete [] lbls[ctr];
  }
    
  delete [] data;
  delete [] ctl;
  delete [] lbls;
  
  }

void miscplot::add_bpdata(const NEWMAT::ColumnVector& vec){
  //add new boxplot column to the boxplot
  // cerr << bp_colctr << "   " << vec.Nrows() << " " << vec(vec.Nrows())<<endl;
  bp_median.push_back(median(vec));
  bp_iqr.push_back(iqr(vec));
  bp_q1.push_back(quantile(vec,1));
  bp_q3.push_back(quantile(vec,3));
  bp_medhi.push_back(std::min(bp_q3[bp_colctr],
			      float(bp_median[bp_colctr]+1.57*bp_iqr[bp_colctr]/std::sqrt((float)vec.Nrows()))));
  bp_medlo.push_back(std::max(bp_q1[bp_colctr],
			      float(bp_median[bp_colctr]-1.57*bp_iqr[bp_colctr]/std::sqrt((float)vec.Nrows()))));
  bp_min.push_back(vec.Minimum());
  bp_max.push_back(vec.Maximum());

  float wishi = bp_min[bp_colctr] , wislo = bp_max[bp_colctr];
  for(int i=1; i<=vec.Nrows(); i++){
    float val = float(vec(i));
    if(val > bp_q3[bp_colctr] + bp_whiskerlength * bp_iqr[bp_colctr]){
      bp_outlierindex.push_back(bp_colctr);
      bp_outliervalue.push_back(val);
    }
    else{
      wishi = std::max(wishi,val);
    }
    if(val < bp_q1[bp_colctr] - bp_whiskerlength * bp_iqr[bp_colctr]){
      bp_outlierindex.push_back(bp_colctr);
      bp_outliervalue.push_back(val);
    }
    else{
      wislo = std::min(wislo,val);
    }
  }

  if(bp_colctr<1){
    bp_minall = bp_min[bp_colctr];
    bp_maxall = bp_max[bp_colctr];
  }else{
    bp_minall = std::min(bp_minall,bp_min[bp_colctr]);
    bp_maxall = std::max(bp_maxall,bp_max[bp_colctr]);
  }

  bp_wishi.push_back(wishi);
  bp_wislo.push_back(wislo);
  bp_colctr++;  
}

void miscplot::add_bpdata(const NEWMAT::Matrix& mat){
  for (int ctr = 1;ctr <=mat.Ncols(); ctr++){
    ColumnVector vec;
    vec = mat.Column(ctr);
    add_bpdata(vec);
  }
}

void miscplot::boxplot(const Matrix& mat, string filename, string title){
  add_bpdata(mat);
  boxplot(filename,title);
}

void  miscplot::boxplot(const  ColumnVector&   vec,  string  filename,  string
title){
  add_bpdata(vec);
  boxplot(filename,title);
}

void miscplot::setscatter(Matrix &data,int width){
  deletescatter();
  GDC_scatter = new GDC_SCATTER_T[ data.Nrows() ]; 
  scat_ctr = 0;
	Matrix tmpdat = data;
	if(data.Ncols()<2){
		tmpdat = data.Column(1);
		for(int ctr=1;ctr<=tmpdat.Nrows(); ctr++)
			tmpdat(ctr,1) = ctr-1;
		tmpdat |= data.Column(1);
	}	

  for(int i=1;i<=tmpdat.Nrows();i++){
    GDC_scatter[scat_ctr].point = tmpdat(i,1);
    GDC_scatter[scat_ctr].val   = tmpdat(i,2);
    GDC_scatter[scat_ctr].width = width;
    GDC_scatter[scat_ctr].color = 0xff0000L;
    GDC_scatter[scat_ctr].ind   = GDC_SCATTER_CIRCLE;
    scat_ctr++;
  }
  GDC_num_scatter_pts = scat_ctr;
}

void miscplot::deletescatter(){
  if(scat_ctr)
  {
	  delete [] GDC_scatter;
    GDC_num_scatter_pts = scat_ctr = 0;
  }
}
 
void miscplot::GDCglobals_reset(){
		deletescatter();
		
		//reset all the globals
		GDC_BGColor   = 0xFFFFFFL;  		/* backgound color (white) */
		GDC_LineColor = GDC_DFLTCOLOR;  /* line color      (black) */
		GDC_SetColor  = &(sc[0]);				/* MATLAB-like colors 		 */
		GDC_bar_width = 75;             /* % of width              */
		GDC_grid  = GDC_TICK_LABELS;
		GDC_hard_graphwidth = 0;
		GDC_hard_size = FALSE;  
		GDC_hard_xorig = 0;
		GDC_hold_img = GDC_DESTROY_IMAGE;
		GDC_image_type = GDC_PNG;
		GDC_interpolations = FALSE;
		GDC_num_scatter_pts = 0;
		GDC_requested_ymax = GDC_NOVALUE;
		GDC_requested_ymin = GDC_NOVALUE;
		GDC_scatter         = NULL;
		GDC_num_scatter_pts = 0;
		GDC_ticks = GDC_TICK_LABELS;
		GDC_title = NULL;
		GDC_title_size = GDC_MEDBOLD;
		GDC_xaxis = TRUE;
		GDC_xaxis_angle = 90;
		GDC_xlabel_ctl = NULL;
		GDC_yaxis = TRUE;
		GDC_yaxis2 = TRUE;
		GDC_ylabel_density = 80;
		GDC_ylabel_fmt = NULL;
}

void miscplot::boxplot(string filename, string title){

		deletescatter();

  int    num_points = bp_colctr+1;

  float* data = new float[2*num_points];
  //unsigned long* extclr = new unsigned long[ num_points ];
  char** lbls = new char*[num_points];
  string s;
 
  data[0]=bp_min[0];
  data[num_points]=bp_min[0];
  for(int ctr2=1; ctr2<num_points; ctr2++){
    data[ctr2] = bp_q1[ctr2-1];
    data[ctr2+num_points] = bp_q3[ctr2-1];
  }
 
  s=string(" ");
  lbls[0]=new char[s.length()+1];
  strcpy(lbls[0],s.c_str());

  if(labels.size()>0)
    for(int ctr=1;ctr<num_points;ctr++){
      if(ctr<=int(labels.size()))
				s = labels[ctr-1];
      else
				s=string(" ");
      lbls[ctr] = new char[s.length()+1];
      strcpy(lbls[ctr],s.c_str());
    }
  else{
    for(int ctr=1;ctr<num_points;ctr++){
      s = int2str(ctr);
      lbls[ctr] = new char[s.length()+1];
      strcpy(lbls[ctr],s.c_str());
    }
  }
  labels.clear();

  // scatter points draw medians, ranges and outliers
	int scatter_points = (60 + 40 + 40 + 500 + bp_outlierindex.size())
		* bp_colctr;
	
  if(bp_notched)
    scatter_points += 120 * bp_colctr;
	
  GDC_SCATTER_T* scat = new GDC_SCATTER_T[ scatter_points ];
	int scat_ctr = 0;

  //medians
  for( int i = 1; i < num_points; i++ ){
    for ( int k = 1; k <= 60; k++ ){
  		scat[scat_ctr].point = float(i)-0.29+0.01*k;
      scat[scat_ctr].val   = bp_median[i-1];
      scat[scat_ctr].width = 1;
      scat[scat_ctr].color = 0xff0000L;
      scat[scat_ctr].ind   = GDC_SCATTER_CIRCLE;
      scat_ctr++;
    }
  }

  if(bp_notched){
    for( int i = 1; i < num_points; i++ ){
      for ( int k = 1; k <= 60; k++ ){
				scat[scat_ctr].point = float(i)-0.29+0.01*k;
				scat[scat_ctr].val   = bp_medhi[i-1];
				scat[scat_ctr].width = 1;
				scat[scat_ctr].color = 0x0000ffL;
				scat[scat_ctr].ind   = GDC_SCATTER_CIRCLE;
				scat_ctr++;
      }
    }
    
    for( int i = 1; i < num_points; i++ ){
      for ( int k = 1; k <= 60; k++ ){
				scat[scat_ctr].point = float(i)-0.29+0.01*k;
				scat[scat_ctr].val   = bp_medlo[i-1];
				scat[scat_ctr].width = 1;
				scat[scat_ctr].color = 0x0000ffL;
				scat[scat_ctr].ind   = GDC_SCATTER_CIRCLE;
				scat_ctr++;
      }
    }
  }

  //top whisker
  for( int i = 1; i < num_points; i++ ){
    for ( int k = 1; k <= 40; k++ ){
      scat[scat_ctr].point = float(i)-0.2+0.01*k;
      scat[scat_ctr].val   = bp_wishi[i-1];
      scat[scat_ctr].width = 1;
      scat[scat_ctr].color = 0xff0000L;
      scat[scat_ctr].ind   = GDC_SCATTER_CIRCLE;
      scat_ctr++;
    }
    for ( int k = 1; k <= 250; k++ ){
      scat[scat_ctr].point = float(i);
      scat[scat_ctr].val   = bp_wishi[i-1]-k*(bp_wishi[i-1]-bp_q3[i-1])/250;
      scat[scat_ctr].width = 1;
      scat[scat_ctr].color = 0xff0000L;
      scat[scat_ctr].ind   = GDC_SCATTER_CIRCLE;
      scat_ctr++;
    }
  }

  //bottom whisker
  for( int i = 1; i < num_points; i++ ){
    for ( int k = 1; k <= 40; k++ ){
      scat[scat_ctr].point = float(i)-0.2+0.01*k;
      scat[scat_ctr].val   = bp_wislo[i-1];
      scat[scat_ctr].width = 1;
      scat[scat_ctr].color = 0xff0000L;
      scat[scat_ctr].ind   = GDC_SCATTER_CIRCLE;
      scat_ctr++;
    }
    for ( int k = 1; k <= 250; k++ ){
      scat[scat_ctr].point = float(i);
      scat[scat_ctr].val   = bp_wislo[i-1]+k*(bp_q1[i-1]-bp_wislo[i-1])/250;
      scat[scat_ctr].width = 1;
      scat[scat_ctr].color = 0xff0000L;
      scat[scat_ctr].ind   = GDC_SCATTER_CIRCLE;
      scat_ctr++;
    }
  }

  //outlier
  for ( int k = 1; k <= int(bp_outlierindex.size()  ); k++ ){
      scat[scat_ctr].point = float(bp_outlierindex[k-1]+1);
      scat[scat_ctr].val   = bp_outliervalue[k-1];
      scat[scat_ctr].width = 8;
      scat[scat_ctr].color = 0xff0000L;
      scat[scat_ctr].ind   = GDC_SCATTER_CIRCLE;
      scat_ctr++;
  }

  int xsize, ysize;
  if(req_xsize>0){
    xsize=req_xsize;
  }else{
    xsize=num_points *80;
  }

  if(req_ysize>0){
    ysize=req_ysize;
  }else{
    ysize=450;
  }

  if(ymin < ymax){
    GDC_requested_ymax = ymax;
    GDC_requested_ymin = ymin;
  }else{
    float range = 0.025*(bp_maxall - bp_minall);
    GDC_requested_ymax = bp_maxall+range;
    GDC_requested_ymin = bp_minall-range;
  }

  GDC_scatter         = scat;
  GDC_num_scatter_pts = scat_ctr;
	GDC_grid = GDC_TICK_LABELS;
	if(gridswapdefault)
		GDC_grid = GDC_TICK_NONE;
		
  if(filename.substr(filename.size()-4,filename.size())!=string(".png"))
    filename += string(".png");

  FILE  *outpng1 = fopen(filename.c_str(), "wb" );
  GDC_image_type     = GDC_PNG;

  GDC_ylabel_fmt = (char*)string("%.3f").c_str();
  GDC_ylabel_density = 30;
  GDC_title = (char*)title.c_str();

  //if (ylabels.size()>0||xlabels.size()>0)
  //  GDC_hold_img = GDC_EXPOSE_IMAGE;

  GDC_BGColor   = 0xFFFFFFL;      /* backgound color (white) */
  GDC_LineColor = 0x000000L;      /* line color      (black) */
  GDC_SetColor  = &(sc[0]);    	  /* set color option */
  GDC_bar_width = 60;             /* (%)              */
  GDC_image_type = GDC_PNG;
  GDC_out_graph( xsize, ysize,    /* width, height */
		 outpng1,             /* open FILE pointer       */
		 GDC_FLOATINGBAR,     /* chart type              */
		 num_points,          /* num points per data set */
		 lbls,                /* X labels array of char* */
		 1,                   /* number of data sets     */
		 (float*)data,        /* data                    */
		 (float*)NULL );      /* right-hand-axis data    */

  fclose( outpng1 );

  gdImagePtr outim2=NULL;

 	if (labels.size()>0||ylabels.size()>0||xlabels.size()>0) 
	 {
	   outpng1=fopen(filename.c_str(), "rb" );
	   outim2=gdImageCreateFromPng(outpng1);
	   fclose( outpng1 );
	 }

  if(labels.size()>0||ylabels.size()>0||xlabels.size()>0){
    outpng1 = fopen(filename.c_str(), "wb" );
    add_legend(outim2, sc);
    gdImagePng(outim, outpng1);
 		gdImageDestroy(outim2);
    fclose( outpng1 );
    if(outim) gdImageDestroy(outim);
 	}
   

  for(int ctr=0;ctr<num_points;ctr++){
    delete [] lbls[ctr];
  }
    
  delete [] data;
  delete [] scat;
  delete [] lbls;

  GDC_scatter         = NULL;
  GDC_num_scatter_pts = 0;
}

void miscplot::histogram(const Matrix& mat, string filename, string title){

  RowVector datam = mat.Row(1);
  int numpoint=datam.Ncols();
  int i,j;
  double scale=std::max(1.0,MISCMATHS::pow((float)10.0,(double)-(std::floor(std::log10(std::min(std::abs(datam.Maximum2(i,j)),std::abs(datam.Minimum2(i,j)))/2.0)))));
  
  if (scale > 10000)
    scale=1;
  //  cerr << datam.Maximum2(i,j) << "  " << datam.Minimum2(i,j) << scale << endl;
  datam *=scale;
  if(scale>100)
    explabel=num2str(std::log10(scale));

  float tmax = datam.Maximum2(i,j);
  float tmin = datam.Minimum2(i,j); 
  float trange = tmax-tmin;

  // check to see if the data is all the same - if it is then return without plotting
  if(trange==0)
    {
      return;
    }

  int bins = (int)floor(MISCMATHS::sqrt(numpoint));
  if (histogram_bins>0)
    bins = histogram_bins;
  float intsize = trange / bins;
  int xlint = (int)ceil(trange/6);
  int binperint = (int)ceil(xlint / intsize);

  //  cerr << tmax << " " << tmin  << " " << trange  << " " << bins << " " << intsize  << " " << xlint  << " " << binperint << " :" << scale <<endl; 
  intsize = float(xlint) / float(binperint);

  float xmin = ceil(std::abs(1.02*tmin)/intsize)*intsize * sign(tmin);
  float xmax = ceil(std::abs(1.02*tmax)/intsize)*intsize * sign(tmax);

  bins = (int)((xmax-xmin)/intsize);

  Matrix bindata(1,bins);
  bindata = 0.0;

  double binsize = (xmax-xmin)/std::max(bins,1);
  for(int ctr = 1; ctr<= datam.Ncols(); ctr++){
    bindata(1,std::max(std::min(int(floor((datam(ctr) - xmin) / binsize) + 1),bins),1))++;
  }
  //  cerr << " calculated bindata " << endl;

  int factor = 1;
  numpoint = factor * bindata.Ncols();
  float* data = new float[numpoint]; 
  for(int ctr1=0;ctr1< numpoint; ctr1++){
      data[ctr1] = bindata(1,(int)floor(float(ctr1 / factor + 1)));
  }
  
  RowVector xax(numpoint);
  for(int ctr1=0;ctr1< numpoint; ctr1++){
    xax(ctr1+1) = xmin + (ctr1+ 0.5) * intsize; 
  }

  char* ctl = new char[numpoint];
  xlint = (int)ceil((xmax-xmin)/9);
  
  ctl[0]=FALSE;
  for(int ctr=1;ctr<numpoint;ctr++){
    int lblpoint = (int)(MISCMATHS::round(abs(xax(ctr))/xlint)*xlint);
    if(xax(ctr)<0) lblpoint *= -1;
     if( (xax(ctr)<lblpoint)&&(xax(ctr+1)>lblpoint) ){
       ctl[ctr]=TRUE;
     }
     else
       ctl[ctr]=FALSE;
  } 

  char** lbls = new char*[numpoint];
 
  string s;

  for(int ctr=0;ctr<numpoint;ctr++){
    if(ctl[ctr]){
      if(scale<2||scale>100)
	s = float2str((int)MISCMATHS::round(xax(ctr+1)/xlint)*xlint,3,2,false);
      else
	s = float2str((float)MISCMATHS::round(xax(ctr+1)/xlint)*xlint/scale,3,2,false);
      lbls[ctr] = new char[s.length()+1];
      strcpy(lbls[ctr],s.c_str());
    }
    else{
      lbls[ctr] = new char[1];
      strcpy(lbls[ctr],string("").c_str());
    }
  }
  
 
  GDC_xlabel_ctl = ctl;

  GDC_BGColor   = 0xFFFFFFL;                  /* backgound color (white) */
  GDC_LineColor = 0x000000L;                  /* line color      (black) */
  GDC_SetColor  = &(sc[0]);                   /* MATLAB-like colors */

  GDC_title = (char*)title.c_str();
  GDC_title_size = GDC_SMALL;

  GDC_ticks = GDC_TICK_LABELS;
  GDC_grid  = GDC_TICK_NONE;
	if(gridswapdefault)
		GDC_grid = GDC_TICK_LABELS;
 
  GDC_yaxis = ylabels.size()>0;
  GDC_yaxis2 = FALSE;
  GDC_xaxis = TRUE;
  GDC_xaxis_angle = 0;
  GDC_bar_width =  75;
  GDC_ylabel_density = 50;
  GDC_requested_ymax = 1.15*bindata.Maximum2(i,j);
  GDC_requested_ymin = 1.15*bindata.Minimum2(i,j);

  //int xsize = std::max(4*numpoint,750);
  //  if(xsize>1200) xsize=(int)floor(std::max(xsize/2.0,750.0));
  int xsize = 600;
  int ysize = 400;

  if(req_xsize>0){
    xsize=req_xsize;
    ysize=req_ysize;
  }

  if(filename.substr(filename.size()-4,filename.size())!=string(".png"))
    filename += string(".png");

  FILE  *outpng1 = fopen(filename.c_str(), "wb" );
  GDC_image_type     = GDC_PNG;

  if (labels.size()>0||ylabels.size()>0||xlabels.size()>0||explabel.length()>0)
    GDC_hold_img = GDC_EXPOSE_IMAGE;
  
  GDC_out_graph( xsize, ysize, outpng1 , GDC_BAR, numpoint, (char**) lbls , 1, data, NULL ); 
  fclose( outpng1 );

  if (labels.size()>0||ylabels.size()>0||xlabels.size()>0||explabel.length()>0){
    outpng1 = fopen(filename.c_str(), "wb" );
    add_legend(GDC_image, sc);
    
    gdImagePng(outim, outpng1); 
    fclose( outpng1 );
    GDC_destroy_image(GDC_image);
    if(outim) gdImageDestroy(outim);
  }

  for(int ctr=0;ctr<numpoint;ctr++){
    delete [] lbls[ctr];
  }
    
  delete [] data;
  delete [] lbls;
  delete [] ctl;
}

  void miscplot::gmmfit(const Matrix& mat, Matrix& mu, Matrix& sig, 
			Matrix& pi, string filename, string title, 
			bool gammamix, float meanoffset, float detailfactor){
  RowVector datam = mat.Row(1);
  int numpoint=datam.Ncols();
  int i,j;
  double scale=std::max(1.0,MISCMATHS::pow((float)10.0,(double)
		-(std::floor(std::log10(std::min(std::abs(datam.Maximum2(i,j)),
		std::abs(datam.Minimum2(i,j)))/2.0)))));

  if (scale > 10000)
    scale=1;
  //cerr << datam.Maximum2(i,j) << "  " << datam.Minimum2(i,j) << scale << endl;
  datam *=scale; mu *= scale; sig *= scale;
  if(scale>100)
    explabel=num2str(std::log10(scale));

  float tmax = datam.Maximum2(i,j);
  float tmin = datam.Minimum2(i,j); 
  float trange = tmax-tmin;
  int bins = (int)floor(MISCMATHS::sqrt(numpoint));
  if (histogram_bins>0)
    bins = histogram_bins;
  float intsize = trange / bins;
  int xlint = (int)ceil(trange/6);
  int binperint = (int)ceil(xlint / intsize);

  intsize = float(xlint) / float(binperint);

  float xmin = ceil(std::abs(1.02*tmin)/intsize)*intsize * sign(tmin);
  float xmax = ceil(std::abs(1.02*tmax)/intsize)*intsize * sign(tmax);

  bins = (int)((xmax-xmin)/intsize);

  Matrix bindata(1,bins);
  bindata = 0.0;

  double binsize = (xmax-xmin)/std::max(bins,1);
  for(int ctr = 1; ctr<= datam.Ncols(); ctr++){
    bindata(1,std::max(std::min(int(floor((datam(ctr) - xmin) / binsize) + 1),bins),1))++;
  }

  int factor = 5;

  bindata = bindata / max(double(factor*bindata.SumAbsoluteValue()),
		double(1.0));

  numpoint = factor * bindata.Ncols();
  float* histdata = new float[numpoint]; 
  for(int ctr1=0;ctr1< numpoint; ctr1++){
   	histdata[ctr1] = bindata(1,(int)floor(float(ctr1 / factor + 1)));
  }
  
  RowVector xax(numpoint);
  for(int ctr1=0;ctr1< numpoint; ctr1++){
    xax(ctr1+1) = xmin + (ctr1+ 0.5) * intsize/factor; 
  }

  char* ctl = new char[numpoint];
  xlint = (int)ceil((xmax-xmin)/(factor*2));
  
  ctl[0]=FALSE;
  for(int ctr=1;ctr<numpoint;ctr++){
    int lblpoint = (int)(MISCMATHS::round(abs(xax(ctr))/xlint)*xlint);
    if(xax(ctr)<0) 
			lblpoint *= -1;
    if( (xax(ctr)<lblpoint)&&(xax(ctr+1)>lblpoint) )	
     	ctl[ctr]=TRUE;
    else
     	ctl[ctr]=FALSE;
  } 

  char** lbls = new char*[numpoint];
 
  string s;
  for(int ctr=0;ctr<numpoint;ctr++){
    if(ctl[ctr]){
      if(scale<2||scale>100)
				s = float2str((int)MISCMATHS::round(xax(ctr+1)/xlint)*xlint,
					3,2,false);
      else
				s = float2str((float)MISCMATHS::round(xax(ctr+1)/xlint)*xlint/scale,
					3,2,false);
      lbls[ctr] = new char[s.length()+1];
      strcpy(lbls[ctr],s.c_str());
    }
    else{
      lbls[ctr] = new char[1];
      strcpy(lbls[ctr],string("").c_str());
    }
  }
  
  // Calculate lines
  int numlines=mu.Ncols()+1;
  for(int ctr1=1;ctr1< sig.Ncols(); ctr1++)
    if(sig(1,ctr1)<0.000001){
      sig(1,ctr1) = 0.000001;
      pi(1,ctr1) = 0.0;
    }

  Matrix fit;
  fit = SP(normpdf(xax,mu,sig),pi.t()*ones(1,numpoint));

  if(null_shift!=0.0)
    fit.Row(1) = pi(1,1)*gammapdf(-(xax + null_shift) + 
			meanoffset, std::abs(mu(1,1) + null_shift), sig(1,1));

  if(gammamix){
    if(pi(1,2)>0.0001)
      fit.Row(2)  = pi(1,2)*gammapdf( xax - meanoffset, 
				std::abs(mu(1,2)), sig(1,2));
    if(pi(1,3)>0.0001)
      fit.Row(3)  = pi(1,3)*gammapdf(-xax + meanoffset, 
				std::abs(mu(1,3)), sig(1,3));
  }
  fit = sum(fit,1) & fit;
  fit = fit / fit.Row(1).SumAbsoluteValue();

  float* linesdata = new float[(numlines)*numpoint];
  for(int ctr1=1;ctr1<=numlines; ctr1++)
    for(int ctr2=1;ctr2<=numpoint; ctr2++)
      linesdata[(ctr1-1)*numpoint + ctr2-1] = fit(ctr1,ctr2);	

  GDC_xlabel_ctl = ctl;

  unsigned long*   sc2 = new unsigned long[numlines];
  for(int ctr=0;ctr<numlines;ctr++) sc2[ctr] = 0xFFDD00;
  sc2[0]=0xFF0000;
  if(gammamix)
    sc2[1]=0x00aa00;

  unsigned long  sc3[64];
  for(int ctr=0;ctr<64;ctr++) sc3[ctr] = sc[ctr];
  for(int ctr=0;ctr<numlines;ctr++)
    sc3[ctr+1] = sc2[ctr];


  GDC_BGColor   = 0xFFFFFFL;                  /* backgound color (white) */
  GDC_LineColor = 0x000000L;                  /* line color      (black) */
  GDC_SetColor  = &(sc2[0]);                  /* MATLAB-like colors */

  GDC_title = (char*)title.c_str();
  GDC_title_size = GDC_SMALL;

  GDC_ticks = GDC_TICK_LABELS;
  GDC_grid  = GDC_TICK_NONE;
	if(gridswapdefault)
		GDC_grid = GDC_TICK_LABELS;
  GDC_yaxis = ylabels.size()>0;
  GDC_yaxis2 = FALSE;
  GDC_xaxis = TRUE;
  GDC_xaxis_angle = 0;
  GDC_bar_width =  75;
  GDC_ylabel_density = 50;
  GDC_requested_ymax = 1.15*bindata.Maximum2(i,j);
  GDC_requested_ymin = 1.15*bindata.Minimum2(i,j);

  int xsize = 600;
  int ysize = 400;

  if(req_xsize>0){
    xsize=req_xsize;
    ysize=req_ysize;
  }

  if((detailfactor > 0.0)&&(fit.Nrows()>2)){
    Matrix tmp(1,fit.Ncols()-2);
    tmp=sum(fit.Rows(3,fit.Nrows()),1);
    float req_max = detailfactor*tmp.Maximum2(i,j);
    if(req_max<0.0000001)
      req_max= 1.05*std::max(bindata.Maximum2(i,j),fit.Maximum2(i,j));
    for(int ctr1=1;ctr1<=numlines; ctr1++)
      for(int ctr2=1;ctr2<=numpoint; ctr2++)
				linesdata[(ctr1-1)*numpoint + ctr2-1] = std::min(req_max,
					linesdata[(ctr1-1)*numpoint + ctr2-1]);

    for(int ctr1=0;ctr1< numpoint; ctr1++)
     	histdata[ctr1] =  std::min(req_max, histdata[ctr1]);
   
    GDC_requested_ymax = req_max;
    GDC_requested_ymin = float(0.0);
  }

  if(filename.substr(filename.size()-4,filename.size())!=string(".png"))
    filename += string(".png");

  FILE  *outpng1 = fopen(filename.c_str(), "wb" );
  GDC_image_type     = GDC_PNG;

  if (labels.size()>0||ylabels.size()>0||xlabels.size()>0||explabel.length()>0)
    GDC_hold_img = GDC_EXPOSE_IMAGE;
  
  GDC_out_graph( xsize, ysize, outpng1 , GDC_COMBO_LINE_AREA, numpoint, 
		(char**) lbls , numlines, linesdata, histdata ); 
  fclose( outpng1 );

  if (labels.size()>0||ylabels.size()>0||xlabels.size()>0||
		explabel.length()>0){
    outpng1 = fopen(filename.c_str(), "wb" );
    add_legend(GDC_image, sc3, TRUE);
   
    gdImagePng(outim, outpng1); 
    fclose( outpng1 );
    GDC_destroy_image(GDC_image);
    if(outim) gdImageDestroy(outim);
    }

  for(int ctr=0;ctr<numpoint;ctr++)
    delete [] lbls[ctr];
    
  delete [] histdata;
  delete [] linesdata;
  delete [] lbls;
  delete [] ctl;
  }




void miscplot::gmmfit(const NEWMAT::Matrix& mat,const NEWMAT::ColumnVector& mu,const NEWMAT::ColumnVector& var,const NEWMAT::ColumnVector& pi,
		      string filename,string title,bool mtype,float offset, float detailfactor){

  RowVector datam = mat.Row(1);
  RowVector _mu(mu.Nrows());
  RowVector _sig(var.Nrows());
  RowVector _pi(pi.Nrows());
  

  _mu  = mu.t();
  _sig = var.t();
  _pi  = pi.t();
  

  int numpoint=datam.Ncols();
  int i,j;
  double scale=std::max(1.0,MISCMATHS::pow((float)10.0,(double)
					   -(std::floor(std::log10(std::min(std::abs(datam.Maximum2(i,j)),
									    std::abs(datam.Minimum2(i,j)))/2.0)))));


  float tmax = datam.Maximum2(i,j);
  float tmin = datam.Minimum2(i,j); 
  float trange = tmax-tmin;
  int bins = (int)floor(MISCMATHS::sqrt(numpoint));

  if (histogram_bins>0)
    bins = histogram_bins;
  float intsize = trange / bins;
  int xlint = (int)ceil(trange/6);
  int binperint = (int)ceil(xlint / intsize);

  intsize = float(xlint) / float(binperint);

  float xmin = ceil(std::abs(1.02*tmin)/intsize)*intsize * sign(tmin);
  float xmax = ceil(std::abs(1.02*tmax)/intsize)*intsize * sign(tmax);

  bins = (int)((xmax-xmin)/intsize);

  Matrix bindata(1,bins);
  bindata = 0.0;

  double binsize = (xmax-xmin)/std::max(bins,1);
  for(int ctr = 1; ctr<= datam.Ncols(); ctr++){
    bindata(1,std::max(std::min(int(floor((datam(ctr) - xmin) / binsize) + 1),bins),1))++;
  }

  int factor = 1;

  bindata = bindata / max(double(factor*bindata.SumAbsoluteValue()),double(1.0));

  numpoint = factor * bindata.Ncols();
  float* histdata = new float[numpoint]; 
  for(int ctr1=0;ctr1< numpoint; ctr1++){
   	histdata[ctr1] = bindata(1,(int)floor(float(ctr1 / factor + 1)));
  }
  
  RowVector xax(numpoint);
  for(int ctr1=0;ctr1< numpoint; ctr1++){
    xax(ctr1+1) = xmin + (ctr1+ 0.5) * intsize/factor; 
  }

  char* ctl = new char[numpoint];
  xlint = (int)ceil((xmax-xmin)/(factor*2));
  
  ctl[0]=FALSE;
  for(int ctr=1;ctr<numpoint;ctr++){
    int lblpoint = (int)(MISCMATHS::round(abs(xax(ctr))/xlint)*xlint);
    if(xax(ctr)<0) 
			lblpoint *= -1;
    if( (xax(ctr)<lblpoint)&&(xax(ctr+1)>lblpoint) )	
     	ctl[ctr]=TRUE;
    else
     	ctl[ctr]=FALSE;
  } 

  char** lbls = new char*[numpoint];
 
  string s;
  for(int ctr=0;ctr<numpoint;ctr++){
    if(ctl[ctr]){
      if(scale<2||scale>100)
	s = float2str((int)MISCMATHS::round(xax(ctr+1)/xlint)*xlint,
		      3,2,false);
      else
	s = float2str((float)MISCMATHS::round(xax(ctr+1)/xlint)*xlint/scale,
		      3,2,false);
      lbls[ctr] = new char[s.length()+1];
      strcpy(lbls[ctr],s.c_str());
    }
    else{
      lbls[ctr] = new char[1];
      strcpy(lbls[ctr],string("").c_str());
    }
  }
  
  // Calculate lines
  int numlines=mu.Nrows()+1;
  
   


  for(int i=1;i<=var.Nrows();i++){
    if(_sig(i)<0.000001){
      _sig(i) = 0.000001;
      _pi(i) = 0.0;
    }
  }

  Matrix fit;
  //OUT(_pi);
  fit = normpdf(xax,_mu,_sig);
  fit = SD(fit,MISCMATHS::repmat(MISCMATHS::sum(fit,2),1,numpoint));
  fit = SP(fit,_pi.t()*ones(1,numpoint));

  //OUT(fit);
  fit = sum(fit,1) & fit;
  //fit = fit / fit.Row(1).SumAbsoluteValue();

  float* linesdata = new float[(numlines)*numpoint];
  for(int ctr1=1;ctr1<=numlines; ctr1++)
    for(int ctr2=1;ctr2<=numpoint; ctr2++)
      linesdata[(ctr1-1)*numpoint + ctr2-1] = fit(ctr1,ctr2);	

  GDC_xlabel_ctl = ctl;

  unsigned long*   sc2 = new unsigned long[numlines];
  for(int ctr=0;ctr<numlines;ctr++) sc2[ctr] = 0xFFDD00;
  sc2[0]=0xFF0000;
  
  unsigned long  sc3[64];
  for(int ctr=0;ctr<64;ctr++) sc3[ctr] = sc[ctr];
  for(int ctr=0;ctr<numlines;ctr++)
    sc3[ctr+1] = sc2[ctr];


  GDC_BGColor   = 0xFFFFFFL;                  /* backgound color (white) */
  GDC_LineColor = 0x000000L;                  /* line color      (black) */
  GDC_SetColor  = &(sc2[0]);                  /* MATLAB-like colors */

  GDC_title = (char*)title.c_str();
  GDC_title_size = GDC_SMALL;

  GDC_ticks = GDC_TICK_LABELS;
  GDC_grid  = GDC_TICK_NONE;
	if(gridswapdefault)
		GDC_grid = GDC_TICK_LABELS;
  GDC_yaxis = ylabels.size()>0;
  GDC_yaxis2 = FALSE;
  GDC_xaxis = TRUE;
  GDC_xaxis_angle = 0;
  GDC_bar_width =  75;
  GDC_ylabel_density = 50;
  GDC_requested_ymax = 1.15*bindata.Maximum2(i,j);
  GDC_requested_ymin = 1.15*bindata.Minimum2(i,j);

  int xsize = 600;
  int ysize = 400;

  if(req_xsize>0){
    xsize=req_xsize;
    ysize=req_ysize;
  }

  

  if(filename.substr(filename.size()-4,filename.size())!=string(".png"))
    filename += string(".png");

  FILE  *outpng1 = fopen(filename.c_str(), "wb" );
  GDC_image_type     = GDC_PNG;

  if (labels.size()>0||ylabels.size()>0||xlabels.size()>0||explabel.length()>0)
    GDC_hold_img = GDC_EXPOSE_IMAGE;
  
  GDC_out_graph( xsize, ysize, outpng1 , GDC_COMBO_LINE_AREA, numpoint, 
		(char**) lbls , numlines, linesdata, histdata ); 
  fclose( outpng1 );

  if (labels.size()>0||ylabels.size()>0||xlabels.size()>0||
		explabel.length()>0){
    outpng1 = fopen(filename.c_str(), "wb" );
    add_legend(GDC_image, sc3, TRUE);
   
    gdImagePng(outim, outpng1); 
    fclose( outpng1 );
    GDC_destroy_image(GDC_image);
    if(outim) gdImageDestroy(outim);
    }

  for(int ctr=0;ctr<numpoint;ctr++)
    delete [] lbls[ctr];
    
  delete [] histdata;
  delete [] linesdata;
  delete [] lbls;
  delete [] ctl;
  }


}





















