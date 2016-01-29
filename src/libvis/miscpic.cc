/* {{{ Copyright */

/*  miscpic - collection of image display and rendering routines

    Stephen Smith, Christian Beckmann and Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 1999-2014 University of Oxford  */

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

/* }}} */

#include "miscpic.h"
#include "gdfontt.h"
#include "gdfonts.h"
#include "gdfontmb.h"
#include "gdfontl.h"
#include "gdfontg.h"

using namespace NEWIMAGE;

namespace MISCPIC{

#include <sstream>

class TextWriter
{
public:
  unsigned int x,y;
  string text;
  TextWriter(int inputX, int inputY, string inputText);
};

TextWriter::TextWriter(int inputX, int inputY, string inputText)
{
  x=inputX;
  y=inputY;
  text=inputText;
}

vector<TextWriter> textWriterVector;
/* {{{ read_lut */

//template <class T>
void miscpic::read_lut()
{
  FILE *tmpfp; 

  if ((tmpfp=fopen(lut.c_str(),"rb"))!=NULL)
    {
      char lutline[10000];
      nlut = 0;rlut.clear();glut.clear();blut.clear();
      while (fgets(lutline, 10000, tmpfp)!=NULL)
	if ( (strncmp(lutline,"<-color{",8)==0) )
	  {
	    float fr, fg, fb;
	    sscanf(lutline+8,"%f,%f,%f",&fr,&fg,&fb);
	    rlut.push_back((int)MISCMATHS::Max(0,MISCMATHS::Min(255,255*fr)));
	    glut.push_back((int)MISCMATHS::Max(0,MISCMATHS::Min(255,255*fg)));
	    blut.push_back((int)MISCMATHS::Max(0,MISCMATHS::Min(255,255*fb)));
	    /*printf("%d %d %d\n",rlut[nlut],glut[nlut],blut[nlut]);*/
	    nlut++;
	  }
      fclose(tmpfp);
    }
}

/* }}} */
/* {{{ float2str */

string float2str(float f, int prec)
  {
    ostringstream os;
    if(std::abs(f)>0.00001){
      os.precision(std::min((int)(prec+ceil(abs(log10(abs(f))))),7));
      os.setf(ios::internal, ios::adjustfield);
      os << f << '\0';
    }
    else
      os << "0.0";
    return os.str();    
  }

/* }}} */
/* {{{ set_minmax */

  //template <class T>
void miscpic::set_minmax(float bgmin, float bgmax, float s1min,
		  float s1max, float s2min, float s2max)
{
  minmax.push_back(bgmin);
  minmax.push_back(bgmax);
  minmax.push_back(s1min);
  minmax.push_back(s1max);
  minmax.push_back(s2min);
  minmax.push_back(s2max);
}

/* }}} */
/* {{{ add_title */

//template <class T>
int miscpic::add_title(int width)
{
  if(title.length()>0)
    {
      
      string tmptitle = title;
      int strlen=0, numns = 1; 

      //work out number of lines etc
      while(tmptitle.find("\n")!=string::npos && tmptitle.find("\n")<tmptitle.length()){
	strlen = max(strlen,(int)tmptitle.find("\n"));
	string tmp = tmptitle;
	tmp=tmp.erase(tmptitle.find("\n"),tmptitle.length());
	tmptitle=tmptitle.erase(0,tmptitle.find("\n")+1);
	numns++;
      }
      strlen = max(strlen,int(tmptitle.length()));

      int linebrk = 3;
      int xsize = max(width, int(strlen*gdFontSmall->w)), 
	ysize = 3*linebrk + numns*(linebrk+gdFontSmall->h);
     
      gdImagePtr tmpim;
      tmpim = gdImageCreateTrueColor(xsize, ysize);
      int fontclr   = gdImageColorResolve(tmpim,240, 240, 240);

      tmptitle = title;
      int xcoor = linebrk, ycoor = 2*linebrk; 
      while(tmptitle.find("\n")!=string::npos && tmptitle.find("\n")<tmptitle.length()){
	strlen = (int)tmptitle.find("\n");
	string tmp = tmptitle;
	tmp=tmp.erase(tmptitle.find("\n"),tmptitle.length());
	tmptitle=tmptitle.erase(0,tmptitle.find("\n")+1);
	xcoor = linebrk + tmpim->sx / 2 - (strlen * gdFontSmall->w)/2;
	char *s = (char*)tmp.c_str();
	gdImageString(tmpim, gdFontSmall, xcoor , ycoor,(unsigned char*) s, fontclr);
	ycoor += linebrk + gdFontSmall->h;
      }

      char *s = (char*)tmptitle.c_str();
      xcoor = linebrk + tmpim->sx / 2 - (tmptitle.length() * gdFontSmall->w)/2;
      gdImageString(tmpim, gdFontSmall, xcoor , ycoor,(unsigned char*) s, fontclr);

      //assemble everything into one image

      gdImagePtr all;
      xsize = tmpim->sx + outim->sx - width;
      ysize = tmpim->sy + outim->sy; 
      all = gdImageCreateTrueColor(xsize, ysize);      
      gdImageCopy(all, tmpim, 0,0,0,0, tmpim->sx, tmpim->sy);
      gdImageCopy(all, outim, 0,tmpim->sy,0,0,outim->sx, outim->sy);

      gdImageDestroy(outim);
      outim = all;
      //outim = gdImageCreateTrueColor(xsize, ysize);
      //gdImageCopy(outim, all, 0, 0 ,0,0,outim->sx, outim->sy);
      
      //gdImageDestroy(all);
      gdImageDestroy(tmpim);
    }
  return 0;
}

/* }}} */
/* {{{ colourbar stuff */

//template <class T>
int miscpic::add_cbar(string cbartype)
{
  if(outim){ 
    if(cbarptr) gdImageDestroy(cbarptr);
    
    cbarptr = gdImageCreateTrueColor(10,outim->sy);
    create_cbar(cbartype);

    gdImagePtr tmpim;
    int xsize = outim->sx + cbarptr->sx;
    int ysize = max(int(outim->sy),int(cbarptr->sy));
    tmpim = gdImageCreateTrueColor(xsize, ysize);
    gdImageCopy(tmpim, outim, 0, (ysize-outim->sy) / 2, 0, 0,outim->sx, outim->sy);
    gdImageCopy(tmpim, cbarptr, outim->sx, (ysize-cbarptr->sy)/2, 0, 0,
		cbarptr->sx, cbarptr->sy);

    gdImageDestroy(outim);
    outim = tmpim;
    
    if(cbarptr){gdImageDestroy(cbarptr); cbarptr=NULL;}

 
//     tmpim = gdImageCreateTrueColor(outim->sx,outim->sy);
//     gdImageCopy(tmpim, outim, 0, 0, 0, 0, outim->sx,outim->sy);
    
    
//     gdImageDestroy(outim);
//     int xsize = tmpim->sx + cbarptr->sx;
//     int ysize = max(int(tmpim->sy),int(cbarptr->sy));
//     outim = gdImageCreateTrueColor(xsize, ysize);
 
//     gdImageCopy(outim, tmpim, 0, (ysize-tmpim->sy) / 2, 0, 0,tmpim->sx, tmpim->sy); 
//     gdImageCopy(outim, cbarptr, tmpim->sx, (ysize-cbarptr->sy)/2, 0, 0, cbarptr->sx, cbarptr->sy);
//     if(cbarptr) gdImageDestroy(cbarptr);
//     if(tmpim) gdImageDestroy(tmpim);
  }

  return 0;
}

//template <class T>
int miscpic::write_cbar(string fname, string cbartype)
{
  
  if(cbarptr){
    gdImageDestroy(cbarptr);   
    cbarptr = NULL;
  }
  
  //output colorbar
  if((create_cbar(cbartype)==0)&&cbarptr){
    // if(debug)
    //  cerr << " write  cbar" << endl;
    FILE *pngout;
    if ((pngout = fopen(fname.c_str(), "wb"))==NULL)
      {
	cerr << "ERROR: Can't open " << fname << "for writing" << endl;
	return(1);
      }
    
    gdImagePng(cbarptr, pngout);
    fclose(pngout);
    gdImageDestroy(cbarptr); cbarptr=NULL;
    return 0;
  }
  else return -1;
}

//template <class T>
int miscpic::create_cbar(string cbartype)
{
  //if(debug)
  // cerr << "  In create_cbar" << endl;
  if(lut.find("render1t")<lut.length()) lut = lutbase + "render1.lut";
  if(lut.find("render2t")<lut.length()) lut = lutbase + "render2.lut";
  read_lut();
  if(nlut<200){
    cerr << "lut " << lut << " contains less than 200 colours" << endl;
    return -1;
  }

  if(minmax.size()<4){
    cerr << "ERROR: set minmax correctly" << endl;
    return -1;
  }

  int width = 14;
  int y_size = 100;
  gdImagePtr redyell, bluecyn, grays;
  
  redyell = gdImageCreateTrueColor(width, y_size);
  bluecyn = gdImageCreateTrueColor(width, y_size);
  grays   = gdImageCreateTrueColor(width, y_size);
  int black = gdImageColorResolve(grays, 0,0,0);
      black = gdImageColorResolve(bluecyn, 0,0,0);
      black = gdImageColorResolve(redyell, 0,0,0);

      //if(debug)
      // cerr << "  black allocated" << endl; 
  //create gray bar
  for(int y=0; y < y_size; y++){
    int col = gdImageColorResolve(grays, 
	      (int)rlut[y],(int)glut[y],(int)blut[y]);
    if(minmax[0]<minmax[1])
      for(int i = 0; i < width; i++)
	gdImageSetPixel(grays, i, y_size-(y+1), col);
    else{
      for(int i = 0; i < width; i++)
	gdImageSetPixel(grays, i, y, col);
    }
  }
  //if(debug)
  //  cerr << "  gray created" << endl;
  //create redyell bar
  for(int y=0; y < y_size; y++){
    int col = gdImageColorResolve(redyell, 
		(int)rlut[y+100],(int)glut[y+100],(int)blut[y+100]);
    if(minmax[2]<minmax[3])
      for(int i = 0; i < width; i++)
	gdImageSetPixel(redyell, i, y_size-(y+1), col);	
    else{
      for(int i = 0; i < width; i++)
	gdImageSetPixel(redyell, i, y, col);
    }
  } 
  //if(debug)
  //  cerr << "  red created" << endl;
  if((nlut>200)&&(minmax.size()>5)){
    for(int y=0; y < y_size; y++){
      int col = gdImageColorResolve(bluecyn, 
		    (int)rlut[y+200],(int)glut[y+200],(int)blut[y+200]);
      if(minmax[4]<minmax[5]){
	for(int i = 0; i < width; i++)
	  gdImageSetPixel(bluecyn, i, y_size-(y+1), col);	
      }
      else{
	for(int i = 0; i < width; i++)
	  gdImageSetPixel(bluecyn, i, y, col);
      }
    }
    //if(debug)
    //  cerr << "  blue created" << endl;
  }
 
  if(minmax[0]>minmax[1]){
    float tmp = minmax[1]; minmax[1] = minmax[0]; minmax[0]=tmp;
  }
  if(minmax[2]>minmax[3]){
    float tmp = minmax[2]; minmax[2] = minmax[3]; minmax[3]=tmp;
  }
  if(minmax[4]>minmax[5]){
    float tmp = minmax[4]; minmax[4] = minmax[5]; minmax[5]=tmp;
  }

  //if(debug){
  //  cerr << "    minmax[0] : " <<  minmax[0] << "    minmax[1] : " <<  minmax[1] << endl;
  //  cerr << "    minmax[2] : " <<  minmax[2] << "    minmax[3] : " <<  minmax[3] << endl;
  //  cerr << "    minmax[4] : " <<  minmax[4] << "    minmax[5] : " <<  minmax[5] << endl;
  // }

  //set up colorbar creation

  vector<string> lbls(minmax.size());
  vector<int> lsize(minmax.size());
  
  for(int ctr= 0; ctr< int(minmax.size()); ctr++){
    lbls[ctr]=float2str(minmax[ctr],2);
    lsize[ctr] = int((lbls[ctr].length() + 2) * gdFontSmall->w);
  }     
  //if(debug)
  //    cerr << "  labels created" << endl;
  int ignores = 0;
  int nrows=1;
  int stepy = 100;
  int topmargin = 10, leftmargin = 10;
  int xsize = 2*leftmargin, ysize = 0;
  int rctr=0, bctr=0, gctr=0;
  
  if(cbartype.length()>5){
    cerr << "ERROR: to many colorbars requested" << endl;
    return -1;
  }
  
  // parse the first letter

  switch(char(cbartype[0])){
  case 'y':  
    xsize += max(lsize[2],lsize[3]) + width;rctr++;  
    break;
  case 'b':  
    if(nlut>200){
      xsize += max(lsize[4],lsize[5]) + width;
      bctr++;
  }
  break;
  case 'g':  
    xsize += max(lsize[0],lsize[1]);gctr++;
    break; 
  default : 
    cerr << endl<< 
      "ERROR: Only use 'y','b' or 'g' as the first colorbar option"
	 << endl; return -1;
  }
  //if(debug)
  //    cerr << "  parsed first symbol" << endl;

  // parse the rest of the colorbar line

  for(int ctr = 1; ctr < int(cbartype.length()); ctr++) {
    switch(char(cbartype[ctr])) {
    case 'y':
      rctr++;
      if(char(cbartype[ctr-1])=='s')
	xsize = max(xsize,leftmargin+width + max(lsize[2],lsize[3]));
      else
	xsize += max(lsize[2],lsize[3]) + width; 
      break;
      
    case 'b':  
      if(nlut>200){
	bctr++;
	if(char(cbartype[ctr-1])=='s')
	  xsize = max(xsize,leftmargin+width + max(lsize[4],lsize[5]));
	else
	  xsize += max(lsize[4],lsize[5]) + width;
      }
      break;
    case 'g':
      gctr++;
      if(char(cbartype[ctr-1])=='s')
	xsize = max(xsize,leftmargin+ width + max(lsize[0],lsize[1]));
      else
	xsize += max(lsize[0],lsize[1]);
      break;
      
    case 's': 
      nrows++; 
      break;
    default : 
      cerr << endl<< 
	"ERROR: Only use 'y','b','g' or 's' to specify type of colorbar"
	   << endl; return -1;
    } 
  }
  //if(debug)
  //    cerr << "  completed parsing" << endl;

  if(cbarptr){
    //if(debug)
    //  cerr << "  cbarptr exists" << endl;
    ysize = cbarptr->sy;
    if((0.9*(cbarptr->sy)/nrows) > 50){
      stepy = int(0.9 * cbarptr->sy / nrows);
      topmargin = int((cbarptr->sy - nrows*stepy)/2);
    }
    else if(cbarptr->sy > 55){
      ignores = 1;
      nrows = 1;
      stepy = int(0.8 * cbarptr->sy / nrows);
      topmargin = int((cbarptr->sy - nrows*stepy)/2);
      xsize = 2*leftmargin + rctr * (max(lsize[2],lsize[3]) + width) +
	gctr * (max(lsize[0],lsize[1]) + width) + bctr * (max(lsize[4],lsize[5]) + width);
    }
    else{
      return -1;
    }
  }
  else{
    ysize = nrows * 100 + 2*topmargin;
  }
  //know now what to do
  //if(debug)
  //   cerr << "  destroy cbarptr" << endl;
  if(cbarptr) gdImageDestroy(cbarptr);
  //if(cbarptr&&debug)
  //  cerr << "  failed" << endl;
  //if(debug)
  //   cerr << "  create cbarptr" << endl;
  cbarptr = gdImageCreateTrueColor(xsize,ysize);
  //if(debug)
  //   cerr << "    allocate black" << endl;
  black = gdImageColorResolve(cbarptr, 0,0,0);

  //define the color for labels
  //int lineclr   = gdImageColorResolve(cbarptr,230,230,230);
  int fontclr   = gdImageColorResolve(cbarptr,255, 0, 0);
    
  int xcoor = leftmargin;
  int ycoor = topmargin;
  gdImagePtr src =  NULL;
  int inmin = 0;
  int inmax = 0;

  //parse colorbar option again and create plot 
  for(int ctr=0; ctr<int(cbartype.length()); ctr++){
    switch(char(cbartype[ctr])) {
    case 'y': 
      src = redyell;
      inmin = 2;
      inmax = 3;
      break;
    case 'g':
      src = grays;
      inmin = 0;
      inmax = 1;
      break;
    case 'b':
      if(nlut>200){
	src = bluecyn;
	inmin=4;
	inmax=5;
      }else src = NULL;
      break;
    case 's':
      if(ignores==0){
	src = NULL;
	xcoor = leftmargin;
	ycoor = ycoor + stepy;
      }else{
	src = NULL;
      }
      break;
    default : break;
    }

    if(src){
      if(stepy == y_size){
	gdImageCopy(cbarptr, src, xcoor,ycoor, 0, 0, width, y_size);
      }
      else{
	gdImageCopyResized(cbarptr, src, xcoor,ycoor, 0, 0, width, stepy+1, width, y_size);
      }
      //gdImageLine(cbarptr, xcoor+10, ycoor, xcoor+18, ycoor, lineclr);
      //gdImageLine(cbarptr, xcoor+10, ycoor+stepy-1, xcoor+18, 
      //		  ycoor+stepy-1, lineclr);
      {
	char *s = (char*)lbls[inmax].c_str();
	gdImageString(cbarptr, gdFontSmall, xcoor + width + 
		      gdFontSmall->w  , ycoor + 1, 
		      (unsigned char*) s, fontclr);
      }
      {
	char *s = (char*)lbls[inmin].c_str();
	gdImageString(cbarptr, gdFontSmall, xcoor + width + 
		      gdFontSmall->w  , ycoor + stepy-1 - 
		      gdFontSmall->h, (unsigned char*) s, fontclr);
      }
      xcoor += width + max(lsize[inmin],lsize[inmax]);
    }
  }

  if(redyell) gdImageDestroy(redyell);
  if(bluecyn) gdImageDestroy(bluecyn);
  if(grays)   gdImageDestroy(grays);
  
  return 0;
}

/* }}} */

/* {{{ write_pic */

//template <class T>
void miscpic::write_pic(char *fname, int width, int height)
{
  remove(fname);
  if ( (nlut>0) || (compare) )
    {
      if(strstr(fname,".png")==0)
	write_ppm(fname,width,height,picr,picg,picb);
      else
	write_png(fname,width,height,picr,picg,picb); 
    }
  else
    {
      if(strstr(fname,".png")==0)
	write_pgm(fname,width,height,picr);
      else
	write_png(fname,width,height,picr,picr,picr); 
    }
}

/* }}} */
/* {{{ write_png */

//template <class T>
int miscpic::write_png ( char *filename, int x_size, int y_size, 
			 unsigned char *r, unsigned char *g, unsigned char *b)
{  
  FILE *pngout;

  if(strstr(filename,".png")==NULL)
    strcat(filename,".png");

  outim = gdImageCreateTrueColor(x_size, y_size);

  for(int x=0; x < x_size; x++){
    for(int y=0; y < y_size; y++){
      int col = gdImageColorResolve(outim, 
				    (int) r[y*x_size+x], 
				    (int) g[y*x_size+x], 
				    (int) b[y*x_size+x]);
      gdImageSetPixel(outim, x, y, col);	
    }
  }

  for (unsigned int i=0;i<textWriterVector.size();i++)
  {
    int fontclr   = gdImageColorResolve(outim, 255, 255, 255);
    unsigned char *s = (unsigned char*)textWriterVector[i].text.c_str();
    gdImageString(outim, gdFontTiny, textWriterVector[i].x , textWriterVector[i].y, s, fontclr);
    //gdImageString(outim, gdFontSmall, textWriterVector[i].x , textWriterVector[i].y, s, fontclr);

  }
  textWriterVector.clear();

  if ( markRight ) {
    int circleWidth=x_size/20;
    gdImageFilledEllipse(outim, x_size-circleWidth, circleWidth, circleWidth, circleWidth, gdImageColorAllocate(outim, 255, 0, 0));
  }

  if(!(cbartype==string(""))) add_cbar(cbartype);
  add_title(x_size);

  if ((pngout = fopen(filename, "wb"))==NULL)
    {
      printf("Can't open %s for writing\n",filename);
      return(1);
    }

  gdImagePng(outim, pngout);
  fclose(pngout);
  gdImageDestroy(outim); outim=NULL;
  return(0);
}

/* }}} */
/* {{{ write_pgm */
 

//template <class T>
int miscpic::write_pgm ( char *filename, int x_size, int y_size,	unsigned char *i )
{
  FILE *ofp;
  int  x, y;

  if ((ofp=fopen(filename,"wb"))==NULL)
    {
      printf("Can't open %s for writing\n",filename);
      return(1);
    }

  fprintf(ofp,"P5\n");
  fprintf(ofp,"%d %d\n",x_size,y_size);
  fprintf(ofp,"255\n");

  for(y=0; y<y_size; y++)
    for(x=0; x<x_size; x++)
      fwrite(&i[y*x_size+x],1,1,ofp);

  fclose(ofp);
  return(0);
}

/* }}} */
/* {{{ write_ppm */

//template <class T>
int miscpic::write_ppm ( char *filename, int x_size, int y_size,
		unsigned char *r, unsigned char *g, unsigned char *b )
{
  FILE *ofp;
  int  x, y;

  if ((ofp=fopen(filename,"wb"))==NULL)
    {
      printf("Can't open %s for writing\n",filename);
      return(1);
    }

  fprintf(ofp,"P6\n");
  fprintf(ofp,"%d %d\n",x_size,y_size);
  fprintf(ofp,"255\n");

  for(y=0; y<y_size; y++)
    for(x=0; x<x_size; x++)
      {
	fwrite(&r[y*x_size+x],1,1,ofp);
	fwrite(&g[y*x_size+x],1,1,ofp);
	fwrite(&b[y*x_size+x],1,1,ofp);
      }

  fclose(ofp);
  return(0);
}

/* }}} */

/* {{{ add Right label */

void miscpic::addRlabel(unsigned char* picr, int p, int width, int size_pic, 
			int alt_size_pic, bool onleft) 
{
  if (!LR_label_flag) { return; }
  int xstart=0;
  if (!onleft) {
    xstart=alt_size_pic-6;
  }
  picr[p+(size_pic-1-1)*width+xstart+1]=255;
  picr[p+(size_pic-1-2)*width+xstart+1]=255;
  picr[p+(size_pic-1-3)*width+xstart+1]=255;
  picr[p+(size_pic-1-4)*width+xstart+1]=255;
  picr[p+(size_pic-1-5)*width+xstart+1]=255;
  picr[p+(size_pic-1-5)*width+xstart+2]=255;
  picr[p+(size_pic-1-5)*width+xstart+3]=255;
  picr[p+(size_pic-1-5)*width+xstart+4]=255;
  picr[p+(size_pic-1-4)*width+xstart+4]=255;
  picr[p+(size_pic-1-3)*width+xstart+4]=255;
  picr[p+(size_pic-1-3)*width+xstart+3]=255;
  picr[p+(size_pic-1-3)*width+xstart+2]=255;
  picr[p+(size_pic-1-2)*width+xstart+3]=255;
  picr[p+(size_pic-1-1)*width+xstart+4]=255;
}

void miscpic::addRlabel(int p, int width, int size_pic, int alt_size_pic, 
			bool onleft)
{
  if (!LR_label_flag) { return; }
  addRlabel(picr,p,width,size_pic,alt_size_pic,onleft);
  addRlabel(picg,p,width,size_pic,alt_size_pic,onleft);
  addRlabel(picb,p,width,size_pic,alt_size_pic,onleft);
}

/* }}} */

/* {{{ write sagittal slice */

//template <class T>
  void miscpic::sag(float xx, int p, int imageWidth)
{
  float yy, zz;
  int   y, z;
  if (xx<0)
    xx=-xx;
  else
    xx*=(float)x_size;
  xx=MISCMATHS::Min(x_size-1.0,xx);
  ostringstream tempBuffer;
  tempBuffer << (int)xx;
  //Image top left is: (p%imageWidth,p/imageWidth)
  if ( writeText ) {
    TextWriter tempTextWriter((int)(p%imageWidth),(int)(p/imageWidth),"X="+tempBuffer.str());
    textWriterVector.push_back(tempTextWriter);
  }
  for(y=0; y<y_size_pic; y++)
    for(z=0; z<z_size_pic; z++)
      {
	yy=MISCMATHS::Min(y_size-1.0,y/inp1.ydim());
	zz=MISCMATHS::Min(z_size-1.00,z/inp1.zdim());

	if (nlut==0) {
	  picr[p+(z_size_pic-1-z)*imageWidth+y] = (unsigned char)MISCMATHS::Min(255,MISCMATHS::Max(0,inp1.interpolate(xx,yy,zz)));

	  if (compare)
	    {
	      picg[p+(z_size_pic-1-z)*imageWidth+y] = picr[p+(z_size_pic-1-z)*imageWidth+y];
	      picb[p+(z_size_pic-1-z)*imageWidth+y] = picr[p+(z_size_pic-1-z)*imageWidth+y];
	      if ( (inp2((int)(xx+0.5),(int)(yy+0.5),(int)(zz+0.5))>0) &&
		   ((y+z)%2>trans) )
		{
		  picr[p+(z_size_pic-1-z)*imageWidth+y] = 255;
		  picg[p+(z_size_pic-1-z)*imageWidth+y] = picb[p+(z_size_pic-1-z)*imageWidth+y] = 0;
		}
	    }
	} else {
	  picr[p+(z_size_pic-1-z)*imageWidth+y]= (unsigned char)MISCMATHS::Min(255,MISCMATHS::Max(0,imr.interpolate(xx,yy,zz)));
	  picg[p+(z_size_pic-1-z)*imageWidth+y]= (unsigned char)MISCMATHS::Min(255,MISCMATHS::Max(0,img.interpolate(xx,yy,zz)));
	  picb[p+(z_size_pic-1-z)*imageWidth+y]= (unsigned char)MISCMATHS::Min(255,MISCMATHS::Max(0,imb.interpolate(xx,yy,zz)));
	}
      }
  // put in "R" label as appropriate
  int icode, jcode, kcode;
  get_axis_orientations(inp1,icode,jcode,kcode);
  // check if y-axis (x-axis in the figure) is left-right
  if (jcode==NIFTI_L2R) { addRlabel(p,imageWidth,z_size_pic,y_size_pic,false); }
  if (jcode==NIFTI_R2L) { addRlabel(p,imageWidth,z_size_pic,y_size_pic,true); }
}

/* }}} */
/* {{{ write coronal slice */

//template <class T>
  void miscpic::cor(float yy, int p, int imageWidth)
{
  float xx, zz;
  int   x, z;
  if (yy<0)
    yy=-yy;
  else
    yy*=(float)y_size;
  yy=MISCMATHS::Min(y_size-1.0,yy);
  ostringstream tempBuffer;
  tempBuffer << (int)yy;
  //Image top left is: (p%imageWidth,p/imageWidth)
  if ( writeText ) {
    TextWriter tempTextWriter((int)(p%imageWidth),(int)(p/imageWidth),"Y="+tempBuffer.str());
    textWriterVector.push_back(tempTextWriter);
  }
  for(x=0; x<x_size_pic; x++)
    for(z=0; z<z_size_pic; z++)
      {
	xx=MISCMATHS::Min(x_size-1.0,x/inp1.xdim());
	zz=MISCMATHS::Min(z_size-1.0,z/inp1.zdim());

	if (nlut==0) {
	  picr[p+(z_size_pic-1-z)*imageWidth+x] = (unsigned char)MISCMATHS::Min(255,MISCMATHS::Max(0,inp1.interpolate(xx,yy,zz)));
	  
	  if (compare)
	    {
	      picg[p+(z_size_pic-1-z)*imageWidth+x] = picr[p+(z_size_pic-1-z)*imageWidth+x];
	      picb[p+(z_size_pic-1-z)*imageWidth+x] = picr[p+(z_size_pic-1-z)*imageWidth+x];
	      if ( (inp2((int)(xx+0.5),(int)(yy+0.5),(int)(zz+0.5))>0) &&
		   ((x+z)%2>trans) )
		{
		  picr[p+(z_size_pic-1-z)*imageWidth+x] = 255;
		  picg[p+(z_size_pic-1-z)*imageWidth+x] = picb[p+(z_size_pic-1-z)*imageWidth+x] = 0;
		}
	    }
	} else {
	  picr[p+(z_size_pic-1-z)*imageWidth+x]= (unsigned char)MISCMATHS::Min(255,MISCMATHS::Max(0,imr.interpolate(xx,yy,zz)));
	  picg[p+(z_size_pic-1-z)*imageWidth+x]= (unsigned char)MISCMATHS::Min(255,MISCMATHS::Max(0,img.interpolate(xx,yy,zz)));
	  picb[p+(z_size_pic-1-z)*imageWidth+x]= (unsigned char)MISCMATHS::Min(255,MISCMATHS::Max(0,imb.interpolate(xx,yy,zz)));
	}
      }
  // put in "R" label as appropriate
  int icode, jcode, kcode;
  get_axis_orientations(inp1,icode,jcode,kcode);
  // check if x-axis (x-axis in the figure) is left-right
  if (icode==NIFTI_L2R) { addRlabel(p,imageWidth,z_size_pic,x_size_pic,false); }
  if (icode==NIFTI_R2L) { addRlabel(p,imageWidth,z_size_pic,x_size_pic,true); }
}

/* }}} */
/* {{{ write axial slice */

//template <class T>
  void miscpic::axi(float zCoord, int p, int imageWidth)
{
  if (zCoord<0)
    zCoord=-zCoord;
  else
    zCoord*=(float)z_size;
  zCoord=MISCMATHS::Max(0,MISCMATHS::Min(z_size-1.0,zCoord));
  ostringstream tempBuffer;
  tempBuffer << (int)zCoord;
  //Image top left is: (p%imageWidth,p/imageWidth)
  if ( writeText ) {
    TextWriter tempTextWriter((int)(p%imageWidth),(int)(p/imageWidth),"Z="+tempBuffer.str());
    textWriterVector.push_back(tempTextWriter);
  }
  for(int x=0; x<x_size_pic; x++)
    for(int y=0; y<y_size_pic; y++)
    {
      float xx(MISCMATHS::Min(x_size-1.0,x/inp1.xdim()));
      float yy(MISCMATHS::Min(y_size-1.0,y/inp1.ydim()));
      if (nlut==0) {
	picr[p+(y_size_pic-1-y)*imageWidth+x] = (unsigned char)MISCMATHS::Min(255,MISCMATHS::Max(0,inp1.interpolate(xx,yy,zCoord)));
	if (compare) {
	  picg[p+(y_size_pic-1-y)*imageWidth+x] = picr[p+(y_size_pic-1-y)*imageWidth+x];
	  picb[p+(y_size_pic-1-y)*imageWidth+x] = picr[p+(y_size_pic-1-y)*imageWidth+x];
	  if ( (inp2((int)(xx+0.5),(int)(yy+0.5),(int)(zCoord+0.5))>0) && ((x+y)%2>trans) ) {
	    picr[p+(y_size_pic-1-y)*imageWidth+x] = 255;
	    picg[p+(y_size_pic-1-y)*imageWidth+x] = picb[p+(y_size_pic-1-y)*imageWidth+x] = 0;
	  }
	}
      } else {
	picr[p+(y_size_pic-1-y)*imageWidth+x]= (unsigned char)MISCMATHS::Min(255,MISCMATHS::Max(0,imr.interpolate(xx,yy,zCoord)));
	picg[p+(y_size_pic-1-y)*imageWidth+x]= (unsigned char)MISCMATHS::Min(255,MISCMATHS::Max(0,img.interpolate(xx,yy,zCoord)));
	picb[p+(y_size_pic-1-y)*imageWidth+x]= (unsigned char)MISCMATHS::Min(255,MISCMATHS::Max(0,imb.interpolate(xx,yy,zCoord)));
      }
    }
  // put in "R" label as appropriate
  int icode, jcode, kcode;
  get_axis_orientations(inp1,icode,jcode,kcode);
  // check if x-axis (x-axis in the figure) is left-right
  if (icode==NIFTI_L2R) { addRlabel(p,imageWidth,y_size_pic,x_size_pic,false); }
  if (icode==NIFTI_R2L) { addRlabel(p,imageWidth,y_size_pic,x_size_pic,true); }
}

/* }}} */

/* {{{ slicer */
int miscpic::slicer(const volume<float>& vol1, const volume<float>& vol2,vector<string> inputOptions,  bool labelSlices, bool debug)
{
  /* {{{ vars */
int  picsetup=0, sagcorskip=0, axiskip=0, height=0, picsize=0;
/* }}} */
  /* {{{ first image stuff  */
writeText=labelSlices;
inp1=vol1;
inp1.setinterpolationmethod(trilinear);
inp1.setextrapolationmethod(zeropad);

x_size = inp1.xsize();
y_size = inp1.ysize();
z_size = inp1.zsize();

size = x_size * y_size * z_size;

if(inp1.max()==inp1.min()){
  cerr << " WARNING: input image is empty " << endl;
}

 lut = lutbase + vol1.getAuxFile() + ".lut";

float intensity_min(inp1.percentile(0.02)), intensity_max(inp1.percentile(0.98));

{
  FILE *tmpfp(NULL);
  if( strlen(vol1.getAuxFile().c_str())>0 &&
      ((tmpfp=fopen(lut.c_str(),"rb"))!=NULL) )
    {
      fclose(tmpfp);
      intensity_max=vol1.getDisplayMaximum();
      intensity_min=vol1.getDisplayMinimum(); 
    }
}

if (intensity_min==intensity_max)
{
  intensity_min=inp1.min();
  intensity_max=inp1.max();
}

float scale = MISCMATHS::Min(inp1.xdim(),MISCMATHS::Min(inp1.ydim(),inp1.zdim()));

/* }}} */
  /* {{{ read second (optional) image */

volume<float> edge(vol2.xsize(),vol2.ysize(),vol2.zsize());
float edgethreshold;

if (vol2.nvoxels() > 1 )
{
  inp2=vol2;

  int x, y, z;
  compare=1;

  if ( (x_size!=inp2.xsize()) || (y_size!=inp2.ysize()) || (z_size!=inp2.zsize()) )
    {
      cerr << "Images aren't the same size!" << endl;
      return 1;
    }

  volume<float> tmpedge(inp2.xsize(),inp2.ysize(),inp2.zsize());
  volume<char>  orient(inp2.xsize(),inp2.ysize(),inp2.zsize());
  
  inp2.setextrapolationmethod(mirror);
  tmpedge.setextrapolationmethod(zeropad);
   
  for(z=inp2.minz(); z<=inp2.maxz(); z++)
    for(y=inp2.miny(); y<=inp2.maxy(); y++)
      for(x=inp2.minx(); x<=inp2.maxx(); x++)
	{
	  float xx=inp2(x+1,y,z)-inp2(x-1,y,z),
	    yy=inp2(x,y+1,z)-inp2(x,y-1,z),
	    zz=inp2(x,y,z+1)-inp2(x,y,z-1);

	  xx*=xx; yy*=yy; zz*=zz;

	  if ( (xx>yy) && (xx>zz) )
	    orient.value(x,y,z)=0;
	  else
	    {
	      if (yy>zz)
		orient.value(x,y,z)=1;
	      else
		orient.value(x,y,z)=2;
	    }

	  tmpedge.value(x,y,z) = (float) xx+yy+zz;
	}

  edge = 0; 
  for(z=tmpedge.minz(); z<=tmpedge.maxz(); z++)
    for(y=tmpedge.miny(); y<=tmpedge.maxy(); y++)
      for(x=tmpedge.minx(); x<=tmpedge.maxx(); x++)
	{
	  if (orient.value(x,y,z)==0)
	    {
	      if ( (tmpedge.value(x,y,z)>tmpedge(x-1,y,z)) 
		   && (tmpedge.value(x,y,z)>=tmpedge(x+1,y,z)) ) 
		edge.value(x,y,z)=std::sqrt(tmpedge.value(x,y,z))/2;
	    }
	  else
	    {
	      if (orient.value(x,y,z)==1)
		{
		  if ( (tmpedge.value(x,y,z)>tmpedge(x,y-1,z)) && (tmpedge.value(x,y,z)>=tmpedge(x,y+1,z)) )
		    edge.value(x,y,z)=std::sqrt(tmpedge.value(x,y,z))/2;
		}
	      else
		{
		  if ( (tmpedge.value(x,y,z)>tmpedge(x,y,z-1)) && 
		       (tmpedge.value(x,y,z)>=tmpedge(x,y,z+1)) )
		    edge.value(x,y,z)=std::sqrt(tmpedge.value(x,y,z))/2;
		}
	    }
	}

  float edge2=edge.percentile(0.02);
  float edge98=edge.percentile(0.98);

  if (edge2==edge98)
    edgethreshold = 0.5*(edge.max()-edge.min())+edge.min();
  else
    edgethreshold = 0.5*(edge98-edge2)+edge2;
  

  for(z=inp2.minz(); z<=inp2.maxz(); z++)
    for(y=inp2.miny(); y<=inp2.maxy(); y++)
      for(x=inp2.minx(); x<=inp2.maxx(); x++)
	if(edge.value(x,y,z) > edgethreshold) inp2.value(x,y,z) = 1; 
	else inp2.value(x,y,z) = 0;
}

/* }}} */
  
 for (unsigned int current=0;current<inputOptions.size();current++ ) {
   char* input=new char[inputOptions[current].size()+1];
   strcpy(input,inputOptions[current].c_str());
   char* theopt=strtok(input," ");
   const char* discard=" ";
   if (strncmp(theopt,"-l",2)==0) { /* lut option lut */
      theopt = strtok(NULL,discard); 
      if (!theopt || strncmp(theopt,"-",1)==0 ) {
         cerr << "Error - must have a LUT name after the -l option." << endl;
	exit(1);
      }
      if(strstr(theopt,"/")==NULL){
	lut = lutbase + string(theopt)+".lut";}
      else
	lut = string(theopt);
   }
   else if (strncmp(theopt, "-s", 2)==0) { /* set scale etc */
      theopt = strtok(NULL,discard);
      if (!theopt || strncmp(theopt,"-",1)==0) {
	 cerr << "Error - must have a scaling factor after the -s option." << endl;
	 exit (1);
      }
      scale /= atof(theopt);
   } else if (strncmp(theopt, "-c", 2)==0) { //mark image
     markRight=true;
   } else if (strncmp(theopt, "-u", 2)==0) /* Do not label the image (i.e. u for unlabeled!)  */
      LR_label_flag = false;
   else if (strncmp(theopt, "-i", 2)==0) { /* set intensity range */
      theopt = strtok(NULL,discard);
      if (!theopt) {
	 cerr << "Error - must set min and max intensities after the -i option." << endl;
	 exit (1);
      }
      intensity_min=atof(theopt);
      theopt = strtok(NULL,discard);
      if (!theopt) {
	 cerr << "Error - must set min and max intensities after the -i option." << endl;
	 exit (1);
      }
      intensity_max=atof(theopt);
   }
   else if (strncmp(theopt, "-e", 2)==0) { /* {{{ user-defined edge threshold? */
      theopt = strtok(NULL,discard);
      if (!theopt) {
	 cerr << "Error - must set threshold value after the -e option." << endl;
	 exit (1);
      }
      edgethresh=atof(theopt);	
      if (edgethresh<0) 
	edgethreshold=-edgethresh;  // negative values are abs (like in -x)
      else 
	edgethreshold = edgethresh*(edge.max()-edge.min()) + edge.min();
      if (inp2.nvoxels()>0) {
	 for(int z=inp2.minz(); z<=inp2.maxz(); z++)
	    for(int y=inp2.miny(); y<=inp2.maxy(); y++)
	       for(int x=inp2.minx(); x<=inp2.maxx(); x++)
		  if ( edge.value(x,y,z) > edgethreshold ) 
		    inp2.value(x,y,z) = 1; 
		  else 
		    inp2.value(x,y,z) = 0;
      }  
   }
   else if (strncmp(theopt, "-t", 2)==0) { /* {{{ transparent edges? */
      cerr << "TRANS = 1" << endl;
      trans=1;
   }
   else if (strncmp(theopt, "-n", 2)==0) /* {{{ nearestneighbour interpolation? */
      inp1.setinterpolationmethod(nearestneighbour);
   else { 
      if ( !picsetup ) {
	 /* {{{ read lut */
	 read_lut();
	 /* {{{ setup scaling and pictures sizes */

	 inp1.setxdim(inp1.xdim()/scale);
	 inp1.setydim(inp1.ydim()/scale);
	 inp1.setzdim(inp1.zdim()/scale);
  
	 x_size_pic = (int)MISCMATHS::round(x_size*inp1.xdim());
	 y_size_pic = (int)MISCMATHS::round(y_size*inp1.ydim());
	 z_size_pic = (int)MISCMATHS::round(z_size*inp1.zdim());

	 //if(debug) printf("%f %f %f %d %d %d\n",inp1.xdim(),inp1.ydim(),inp1.zdim(),x_size_pic,y_size_pic,z_size_pic);

	 picsize = MISCMATHS::Max(x_size_pic,MISCMATHS::Max(y_size_pic,z_size_pic));
	 picsize = picsize*picsize*MISCMATHS::Max(10,2*z_size);

	 //if(debug) printf("%d %d %d %d\n",x_size_pic,y_size_pic,z_size_pic,picsize);*/

	 picr = (unsigned char*)malloc(picsize);
	 picg = (unsigned char*)malloc(picsize);
	 picb = (unsigned char*)malloc(picsize);

	 sagcorskip=0;
	 axiskip=0;
	 height=MISCMATHS::Max(y_size_pic,z_size_pic);

	 if (y_size_pic>z_size_pic)
	   sagcorskip=(y_size_pic-z_size_pic)/2;
	 else
	   axiskip=(z_size_pic-y_size_pic)/2;

	 /* {{{ rescale image intensity */

	 if(debug)
	   cerr << "Entries in LUT: " << nlut << endl;

	 if(debug)
	   cerr << " Intensity min, max : " << intensity_min << "  " << intensity_max << endl;


	 if (nlut==0) { /* otherwise (if using luts) then the display range should already have been set correctly */
	   for(int z=inp1.minz(); z<=inp1.maxz(); z++)
	     for(int y=inp1.miny(); y<=inp1.maxy(); y++)
	       for(int x=inp1.minx(); x<=inp1.maxx(); x++)
		 if(intensity_max>intensity_min)
		   inp1.value(x,y,z) =  (int)(((inp1.value(x,y,z) - intensity_min ) * 255.001 ) / (intensity_max-intensity_min));
		 else
		   inp1.value(x,y,z) = 0;
	 } else {
	   imr = inp1;
	   img = inp1;
	   imb = inp1;
	   imr.setextrapolationmethod(zeropad);
	   img.setextrapolationmethod(zeropad);
	   imb.setextrapolationmethod(zeropad);

	   for(int z=inp1.minz(); z<=inp1.maxz(); z++)
	     for(int y=inp1.miny(); y<=inp1.maxy(); y++)
	       for(int x=inp1.minx(); x<=inp1.maxx(); x++) {	       
		  int tmp; 
		  if(intensity_max>intensity_min)
		     tmp = (int)(((float)inp1.value(x,y,z) - intensity_min ) * (nlut-1) / (intensity_max-intensity_min));   
		  else
		     tmp = 0;
		  tmp = MISCMATHS::Min(nlut-1,MISCMATHS::Max(0,tmp));
		  imr.value(x,y,z) = rlut[tmp];
		  img.value(x,y,z) = glut[tmp];
		  imb.value(x,y,z) = blut[tmp];
	       }
	   if(debug){
	      save_volume(imr,"IMG_R");
	      save_volume(img,"IMG_G");
	      save_volume(imb,"IMG_B");
	   }
	 }
	 picsetup=1;
      }
      memset(picr,0,picsize);
      memset(picg,0,picsize);
      memset(picb,0,picsize);
      /* {{{ can specify a title for each plot using -T */
      if ( !strncmp(theopt, "-T", 2) ) {
	 string tmpstr = string("");
	 int otheropt;
	 do {
	   theopt = strtok(NULL,discard);
	   otheropt = 1;
	   otheropt=(strncmp(theopt, "-x", 2)&&
		    strncmp(theopt, "-y", 2)&&
		    strncmp(theopt, "-z", 2)&&
		    strncmp(theopt, "-a", 2)&&
		    strncmp(theopt, "-A", 2)&&
		    strncmp(theopt, "-S", 2));
	   if( otheropt ) 
	     tmpstr+= string(theopt) + " ";
	 } while( otheropt );
	 set_title(tmpstr);
      }
      if (!strncmp(theopt, "-x", 2)) { /* {{{ write sagittal slice */
	 theopt = strtok(NULL,discard);
	 if (theopt==NULL) {
	   cerr << "Error - must have a slice choice after the -x option." << endl;
	   exit (1);
	 }
	 float slice=atof(theopt);
	 theopt = strtok(NULL,discard);
	 if (theopt==NULL) {
	   cerr << "Error - must have a picture output filename when using the -x option." << endl;
	   exit (1);
	 }
	 sag(slice, sagcorskip*y_size_pic, y_size_pic);
	 write_pic(theopt,y_size_pic,height);
      }
      else if (!strncmp(theopt, "-y", 2)) { /* {{{ write coronal slice */
	 theopt = strtok(NULL,discard);
	 if (theopt==NULL) {
	   cerr << "Error - must have a slice choice after the -y option." << endl;
	   exit (1);
	 }
	 float slice=atof(theopt);
	 theopt = strtok(NULL,discard);
	 if (theopt==NULL) {
	   cerr << "Error - must have a picture output filename when using the -y option." << endl;
	   exit (1);
	 }
	 cor(slice, sagcorskip*x_size_pic, x_size_pic);
	 write_pic(theopt,x_size_pic,height);
      }
      else if (!strncmp(theopt, "-z", 2)) { /* {{{ write axial slice */
	 theopt = strtok(NULL,discard);
	 if (theopt==NULL) {
	   cerr << "Error - must have a slice choice after the -z option." << endl;
	   exit (1);
	 }
	 float slice=atof(theopt);

	 theopt = strtok(NULL,discard);
	 if (theopt==NULL) {
	   cerr << "Error - must have a picture output filename when using the -z option." << endl;
	   exit (1);
	 }
	 axi(slice, axiskip*x_size_pic, x_size_pic);
	 write_pic(theopt,x_size_pic,height);
      }
      else if (!strncmp(theopt, "-a", 2)) { /* {{{ write default slices */
	theopt = strtok(NULL,discard);
	if (theopt==NULL) {
	  cerr << "Error - must have a picture output filename after the -a option." << endl;
	  exit (1);
	}
	int width=y_size_pic+2*x_size_pic;
	sag(0.5,                         sagcorskip*width, width);
	cor(0.5, y_size_pic            + sagcorskip*width, width);
	axi(0.5, y_size_pic+x_size_pic + axiskip*width,    width);
	write_pic(theopt,width,height);
      }
      else if ( (!strncmp(theopt, "-A", 2)) || (!strncmp(theopt, "-S", 2)) ) { /* {{{ write all slices */
	int maxwidth, nx, ny, width, row, column, height, z=0, sample=1, slices;
	if ( !strcmp(theopt, "-S") ) {
	  theopt = strtok(NULL,discard);
	  if (theopt==NULL) {
	    cerr << "Error - must set the <sample> option choice when using the -S option." << endl;
	    exit (1);
	  }
	  sample=atoi(theopt);
	}
	theopt = strtok(NULL,discard);
	if (theopt==NULL) {
	  cerr << "Error - must set the <width> option when using the -S or -A options." << endl;
	  exit (1);
	}
	maxwidth=atoi(theopt);
	theopt = strtok(NULL,discard);
	if (theopt==NULL) {
	  cerr << "Error - must have a picture output filename after the -S or -A options." << endl;
	  exit (1);
	}

	slices=(int)(ceil(float(z_size / sample)));
	nx=MISCMATHS::Min(MISCMATHS::Max(maxwidth/x_size_pic,1),slices);
	ny=(int)ceil( ((float)slices) / nx);
	width=x_size_pic*nx;
	height=y_size_pic*ny;
	/*printf("%d %d %d %d %d\n",slices,nx,ny,width,height);*/
	for(row=0;row<ny;row++)
	  for(column=0;column<nx;column++) {
	    if (z<z_size)
	      axi(-z, row*y_size_pic*width + x_size_pic*column, width);
	    z+=sample;
	  }

	write_pic(theopt,width,height);
      }
   }
   delete [] input;
 }
 return(0);
}



int miscpic::slicer(const volume<float>& vol1, const volume<float>& vol2, const char *opts,  bool labelSlices, bool debug)
{
  /* {{{ vars */

int  picsetup=0, sagcorskip=0, axiskip=0, height=0, picsize=0;

/* }}} */
  /* {{{ first image stuff  */
writeText=labelSlices;
inp1=vol1;
inp1.setinterpolationmethod(trilinear);
inp1.setextrapolationmethod(zeropad);

x_size = inp1.xsize();
y_size = inp1.ysize();
z_size = inp1.zsize();

size = x_size * y_size * z_size;

if(inp1.max()==inp1.min()){
  cerr << " WARNING: input image is empty " << endl;
}

 lut = lutbase + vol1.getAuxFile() + ".lut";

float intensity_min(inp1.percentile(0.02)), intensity_max(inp1.percentile(0.98));

{

  FILE *tmpfp(NULL);
  if( strlen(vol1.getAuxFile().c_str())>0 &&
      ((tmpfp=fopen(lut.c_str(),"rb"))!=NULL) )
    {
      fclose(tmpfp);
      intensity_max=vol1.getDisplayMaximum();
      intensity_min=vol1.getDisplayMinimum(); 
    }
}

if (intensity_min==intensity_max)
{
  intensity_min=inp1.min();
  intensity_max=inp1.max();
}

float scale = MISCMATHS::Min(inp1.xdim(),MISCMATHS::Min(inp1.ydim(),inp1.zdim()));

/* }}} */
  /* {{{ read second (optional) image */

if (vol2.nvoxels() > 1 )
{
  inp2=vol2;

  int x, y, z;
  compare=1;

  if ( (x_size!=inp2.xsize()) || (y_size!=inp2.ysize()) || (z_size!=inp2.zsize()) )
    {
      cerr << "Images aren't the same size!" << endl;
      return 1;
    }

  volume<float> edge(inp2.xsize(),inp2.ysize(),inp2.zsize());
  volume<float> tmpedge(inp2.xsize(),inp2.ysize(),inp2.zsize());
  volume<char>  orient(inp2.xsize(),inp2.ysize(),inp2.zsize());
  
  inp2.setextrapolationmethod(mirror);
  tmpedge.setextrapolationmethod(zeropad);
   
  for(z=inp2.minz(); z<=inp2.maxz(); z++)
    for(y=inp2.miny(); y<=inp2.maxy(); y++)
      for(x=inp2.minx(); x<=inp2.maxx(); x++)
	{
	  float xx=inp2(x+1,y,z)-inp2(x-1,y,z),
	    yy=inp2(x,y+1,z)-inp2(x,y-1,z),
	    zz=inp2(x,y,z+1)-inp2(x,y,z-1);

	  xx*=xx; yy*=yy; zz*=zz;

	  if ( (xx>yy) && (xx>zz) )
	    orient.value(x,y,z)=0;
	  else
	    {
	      if (yy>zz)
		orient.value(x,y,z)=1;
	      else
		orient.value(x,y,z)=2;
	    }

	  tmpedge.value(x,y,z) = (float) xx+yy+zz;
	}

  edge = 0; 
  for(z=tmpedge.minz(); z<=tmpedge.maxz(); z++)
    for(y=tmpedge.miny(); y<=tmpedge.maxy(); y++)
      for(x=tmpedge.minx(); x<=tmpedge.maxx(); x++)
	{
	  if (orient.value(x,y,z)==0)
	    {
	      if ( (tmpedge.value(x,y,z)>tmpedge(x-1,y,z)) 
		   && (tmpedge.value(x,y,z)>=tmpedge(x+1,y,z)) ) 
		edge.value(x,y,z)=std::sqrt(tmpedge.value(x,y,z))/2;
	    }
	  else
	    {
	      if (orient.value(x,y,z)==1)
		{
		  if ( (tmpedge.value(x,y,z)>tmpedge(x,y-1,z)) && (tmpedge.value(x,y,z)>=tmpedge(x,y+1,z)) )
		    edge.value(x,y,z)=std::sqrt(tmpedge.value(x,y,z))/2;
		}
	      else
		{
		  if ( (tmpedge.value(x,y,z)>tmpedge(x,y,z-1)) && 
		       (tmpedge.value(x,y,z)>=tmpedge(x,y,z+1)) )
		    edge.value(x,y,z)=std::sqrt(tmpedge.value(x,y,z))/2;
		}
	    }
	}

  float edge2=edge.percentile(0.02);
  float edge98=edge.percentile(0.98);

  float edgethreshold;
  if (edge2==edge98)
    edgethreshold = 0.5*(edge.max()-edge.min())+edge.min();
  else
    edgethreshold = 0.5*(edge98-edge2)+edge2;

  for(z=inp2.minz(); z<=inp2.maxz(); z++)
    for(y=inp2.miny(); y<=inp2.maxy(); y++)
      for(x=inp2.minx(); x<=inp2.maxx(); x++)
	if(edge.value(x,y,z) > edgethreshold) inp2.value(x,y,z) = 1; 
	else inp2.value(x,y,z) = 0;
}

/* }}} */

  char *theopt;
  char t[10243]; /* needs to be of fixed width for strtok to work */
  //const char *discard = " ,;:'][}{";
  const char *discard = " ";

  strcpy(t, opts);
  theopt=strtok(t,discard);

  while ( theopt ) {

    if (strncmp(theopt,"-l",2)==0)
      /* {{{ lut option lut */

{
  theopt = strtok(NULL,discard); 
  if (!theopt || strncmp(theopt,"-",1)==0 )
    {
      cerr << "Error - must have a LUT name after the -l option." << endl;
      exit (1);
    }

  if(strstr(theopt,"/")==NULL){
    lut = lutbase + string(theopt)+".lut";}
  else{
    lut = string(theopt);}
}

/* }}} */
    else if (strncmp(theopt, "-s", 2)==0)
      /* {{{ set scale etc */

{
  theopt = strtok(NULL,discard);
  if (!theopt || strncmp(theopt,"-",1)==0)
    {
      cerr << "Error - must have a scaling factor after the -s option." << endl;
      exit (1);
    }
  scale /= atof(theopt);
}

/* }}} */
    else if (strncmp(theopt, "-u", 2)==0)
      /* {{{ Do not label the image (i.e. u for unlabeled!) */
{
  LR_label_flag = false;
}

/* }}} */
    else if (strncmp(theopt, "-i", 2)==0)
      /* {{{ set intensity range */

{
  theopt = strtok(NULL,discard);
  if (!theopt)
    {
      cerr << "Error - must set min and max intensities after the -i option." << endl;
      exit (1);
    }
  intensity_min=atof(theopt);

  theopt = strtok(NULL,discard);
  if (!theopt)
    {
      cerr << "Error - must set min and max intensities after the -i option." << endl;
      exit (1);
    }
  intensity_max=atof(theopt);
}

/* }}} */
    else if (strncmp(theopt, "-t", 2)==0)
      /* {{{ transparent edges? */

trans=1;

/* }}} */
    else if (strncmp(theopt, "-n", 2)==0)
      /* {{{ nearestneighbour interpolation? */

inp1.setinterpolationmethod(nearestneighbour);

/* }}} */
    else
      { 
	if ( !picsetup )
	  {
	    /* {{{ read lut */

read_lut();

/* }}} */
	    /* {{{ setup scaling and pictures sizes */

inp1.setxdim(inp1.xdim()/scale);
inp1.setydim(inp1.ydim()/scale);
inp1.setzdim(inp1.zdim()/scale);
  
x_size_pic = (int)MISCMATHS::round(x_size*inp1.xdim());
y_size_pic = (int)MISCMATHS::round(y_size*inp1.ydim());
z_size_pic = (int)MISCMATHS::round(z_size*inp1.zdim());

/*if(debug)
  printf("%f %f %f %d %d %d\n",inp1.xdim(),inp1.ydim(),inp1.zdim(),x_size_pic,y_size_pic,z_size_pic);*/

picsize = MISCMATHS::Max(x_size_pic,MISCMATHS::Max(y_size_pic,z_size_pic));
picsize = picsize*picsize*MISCMATHS::Max(10,2*z_size);

/*if(debug)
  printf("%d %d %d %d\n",x_size_pic,y_size_pic,z_size_pic,picsize);*/

picr = (unsigned char*)malloc(picsize);
picg = (unsigned char*)malloc(picsize);
picb = (unsigned char*)malloc(picsize);

sagcorskip=0;
axiskip=0;
height=MISCMATHS::Max(y_size_pic,z_size_pic);

if (y_size_pic>z_size_pic)
  sagcorskip=(y_size_pic-z_size_pic)/2;
else
  axiskip=(z_size_pic-y_size_pic)/2;

 /* }}} */
	    /* {{{ rescale image intensity */

if(debug)
  cerr << "Entries in LUT: " << nlut << endl;

if(debug)
  cerr << " Intensity min, max : " << intensity_min 
       << "  " << intensity_max << endl;


if (nlut==0) /* otherwise (if using luts) then the display range should already have been set correctly */
{ 
  for(int z=inp1.minz(); z<=inp1.maxz(); z++)
    for(int y=inp1.miny(); y<=inp1.maxy(); y++)
      for(int x=inp1.minx(); x<=inp1.maxx(); x++)
	if(intensity_max>intensity_min)
	  inp1.value(x,y,z) =  (int)(((inp1.value(x,y,z) 
				       - intensity_min ) * 255.001 ) 
				     / (intensity_max-intensity_min));
	else
	  inp1.value(x,y,z) = 0;
} else {
  imr = inp1;
  img = inp1;
  imb = inp1;
  imr.setextrapolationmethod(zeropad);
  img.setextrapolationmethod(zeropad);
  imb.setextrapolationmethod(zeropad);

  for(int z=inp1.minz(); z<=inp1.maxz(); z++)
    for(int y=inp1.miny(); y<=inp1.maxy(); y++)
      for(int x=inp1.minx(); x<=inp1.maxx(); x++)
  {
    int tmp; 
    if(intensity_max>intensity_min)
      tmp = (int)(((float)inp1.value(x,y,z) - intensity_min ) * (nlut-1) 
		  / (intensity_max-intensity_min));   
    else
      tmp = 0;

    tmp = MISCMATHS::Min(nlut-1,MISCMATHS::Max(0,tmp));
    imr.value(x,y,z) = rlut[tmp];
    img.value(x,y,z) = glut[tmp];
    imb.value(x,y,z) = blut[tmp];
  }
  if(debug){
      save_volume(imr,"IMG_R");
      save_volume(img,"IMG_G");
      save_volume(imb,"IMG_B");
  }
}

/* }}} */
	    picsetup=1;
	  }
	/* {{{ zero picture */

memset(picr,0,picsize);
memset(picg,0,picsize);
memset(picb,0,picsize);

/* }}} */
        /* {{{ can specify a title for each plot using -T */

	if(!strncmp(theopt, "-T", 2))
	  {
	    string tmpstr = string("");
	    int otheropt;
	    do{
	      theopt = strtok(NULL,discard);
	      otheropt = 1;
	      otheropt=(strncmp(theopt, "-x", 2)&&
		    strncmp(theopt, "-y", 2)&&
		    strncmp(theopt, "-z", 2)&&
		    strncmp(theopt, "-a", 2)&&
		    strncmp(theopt, "-A", 2)&&
		    strncmp(theopt, "-S", 2));
	      if(otheropt) tmpstr+= string(theopt) + " ";
	    }while(otheropt);
	    set_title(tmpstr);
	  }

/* }}} */

	if (!strncmp(theopt, "-x", 2))
	  /* {{{ write sagittal slice */

{
  theopt = strtok(NULL,discard);
  if (theopt==NULL)
    {
      cerr << "Error - must have a slice choice after the -x option." << endl;
      exit (1);
    }
  float slice=atof(theopt);

  theopt = strtok(NULL,discard);
  if (theopt==NULL)
    {
      cerr << "Error - must have a picture output filename when using the -x option." << endl;
      exit (1);
    }

  sag(slice, sagcorskip*y_size_pic, y_size_pic);
  write_pic(theopt,y_size_pic,height);
}

/* }}} */
	else if (!strncmp(theopt, "-y", 2))
	  /* {{{ write coronal slice */

{
  theopt = strtok(NULL,discard);
  if (theopt==NULL)
    {
      cerr << "Error - must have a slice choice after the -y option." << endl;
      exit (1);
    }
  float slice=atof(theopt);

  theopt = strtok(NULL,discard);
  if (theopt==NULL)
    {
      cerr << "Error - must have a picture output filename when using the -y option." << endl;
      exit (1);
    }

  cor(slice, sagcorskip*x_size_pic, x_size_pic);
  write_pic(theopt,x_size_pic,height);
}

/* }}} */
	else if (!strncmp(theopt, "-z", 2))
	  /* {{{ write axial slice */

{
  theopt = strtok(NULL,discard);
  if (theopt==NULL)
    {
      cerr << "Error - must have a slice choice after the -z option." << endl;
      exit (1);
    }
  float slice=atof(theopt);

  theopt = strtok(NULL,discard);
  if (theopt==NULL)
    {
      cerr << "Error - must have a picture output filename when using the -z option." << endl;
      exit (1);
    }

  axi(slice, axiskip*x_size_pic, x_size_pic);
  write_pic(theopt,x_size_pic,height);
}

/* }}} */
	else if (!strncmp(theopt, "-a", 2))
	  /* {{{ write default slices */

{
  theopt = strtok(NULL,discard);
  if (theopt==NULL)
    {
      cerr << "Error - must have a picture output filename after the -a option." << endl;
      exit (1);
    }

  int width=y_size_pic+2*x_size_pic;

  sag(0.5,                         sagcorskip*width, width);
  cor(0.5, y_size_pic            + sagcorskip*width, width);
  axi(0.5, y_size_pic+x_size_pic + axiskip*width,    width);

  write_pic(theopt,width,height);
}

/* }}} */
	else if ( (!strncmp(theopt, "-A", 2)) || (!strncmp(theopt, "-S", 2)) )
	  /* {{{ write all slices */

{
  int maxwidth, nx, ny, width, row, column, height, z=0, sample=1, slices;

  if (!strcmp(theopt, "-S"))
    {
      theopt = strtok(NULL,discard);
      if (theopt==NULL)
	{
	  cerr << "Error - must set the <sample> option choice when using the -S option." << endl;
	  exit (1);
	}
      sample=atoi(theopt);
    }

  theopt = strtok(NULL,discard);
  if (theopt==NULL)
    {
      cerr << "Error - must set the <width> option when using the -S or -A options." << endl;
      exit (1);
    }
  maxwidth=atoi(theopt);

  theopt = strtok(NULL,discard);
  if (theopt==NULL)
    {
      cerr << "Error - must have a picture output filename after the -S or -A options." << endl;
      exit (1);
    }

  slices=(int)(ceil(float(z_size / sample)));
  nx=MISCMATHS::Min(MISCMATHS::Max(maxwidth/x_size_pic,1),slices);
  ny=(int)ceil( ((float)slices) / nx);
  width=x_size_pic*nx;
  height=y_size_pic*ny;
  /*printf("%d %d %d %d %d\n",slices,nx,ny,width,height);*/
    
  for(row=0;row<ny;row++)
    for(column=0;column<nx;column++)
      {
	if (z<z_size)
	  axi(-z, row*y_size_pic*width + x_size_pic*column, width);
	z+=sample;
      }

  write_pic(theopt,width,height);
}

/* }}} */
      }
    theopt = strtok(NULL,discard);
}

  return(0);
}

/* }}} */
/* {{{ overlay */

//template <class T>
int miscpic::overlay(volume<float>& newvol, volume<float>& bg, 
			volume<float>& s1,volume<float>& s2, 
			float bgmin, float bgmax, float s1min,
			float s1max, float s2min, float s2max, 
			int colour_type, int checker, 
			string cbarfname, string cbartype, bool out_int, bool dbg)
{
  float latitude = 0.00001, hrange, A , B;
  int ns = 1;
  volume<float> tmpvol;
  float factor1 = 1.0, factor2 = 1.0;
  debug = false;

  set_minmax(bgmin,bgmax,s1min,s1max,s2min,s2max);

  if(debug){
    cerr << "  libvis - overlay " << endl
	 << "  bgmin/bgmax : " << bgmin << " " << bgmax<< endl
	 << "  s1min/s1max : " << s1min << " " << s1max<< endl
	 << "  s2min/s2max : " << s2min << " " << s2max<< endl
	 << "  colour_type : " << colour_type << " checker  : " << checker  << endl
	 << "  cbarfname   : " << cbarfname   << " cbartype : " << cbartype << endl
         << "  hdrinfo lut : " << newvol.getAuxFile() << endl
	 << "  hdr cal_min : " << newvol.getDisplayMinimum() << endl
	 << "  hdr cal_max : " << newvol.getDisplayMaximum() << endl
         << "  bg perc.    : " << bg.percentile(0.1) << " " << bg.percentile(0.9) << endl
	 << "  s1 perc.    : " << s1.percentile(0.1) << " " << s1.percentile(0.9) << endl;
    if(s2min != s2max)
      cerr << "  s2 perc.    : " << s2.percentile(0.1) << " " << s2.percentile(0.9) << endl;

    if(&s1==&s2)
      cerr << "  s1 and s2 point to the same volume " << endl;
    else
      cerr << "  s1 and s2 point to different volumes " << endl;
  }


  if(s1max<s1min){
    factor1 *= -1.0;
    s1min *= -1.0; s1max *= -1.0;
  }
  
  if(s2min != s2max) ns=2;

  if((ns==2)&&(s2max<s2min)){
    factor2 *= -1.0;
    s2min *= -1.0 ; s2max *= -1.0;
  }


  //  if( (dtype(newvol) == DT_FLOAT)&&(s1min != s1max) )
  if( ( !out_int )&&(s1min != s1max) )
    {
      A = (13*s1min - 7*s1max ) / 6;
      B = (20*s1min + A) / 21;
    }
  else
    {
      A = 0.0;
      B = 32000.0 / (1 + ( ns * ns ));
    }

  hrange = B - A;
  tmpvol = bg;
  tmpvol.setDisplayMaximumMinimum(A + ( hrange * ( 1 + ( ns * ns )) ) ,A);
  for(int z = bg.minz(); z<=bg.maxz();z++){
    for(int y = bg.miny(); y<=bg.maxy();y++){
      for(int x = bg.minx(); x<=bg.maxx();x++){
	
	float newbg;
	if(bgmin==bgmax)
	  newbg = float(0.0);
	else
	  newbg = (MISCMATHS::Min(MISCMATHS::Max(bg(x,y,z),bgmin),bgmax)
		       - bgmin ) / (bgmax-bgmin);

	float news1;
	if(s1max==s1min)
	  news1 = float(0.0);
	else
	  news1 = (MISCMATHS::Min(factor1*s1(x,y,z),s1max)-s1min) / (s1max-s1min); 
	int   mask1 = ((factor1*s1(x,y,z)) >= (s1min - latitude));          

	if(checker > 0)
	  mask1 *= (x+y+z) % 2;

	if (ns == 1){
	  if ( colour_type == 0 )
	    tmpvol.value(x,y,z) = A + newbg*0.95*hrange*(1-mask1) 
	      + ( news1*0.9*hrange + 1.05*hrange ) * mask1;
	  else
	    tmpvol.value(x,y,z) = A + newbg*0.95*hrange*(1-mask1) +
	      ( ((int)(news1*5*0.999)) * 0.2*hrange +
		newbg*hrange*0.5*0.2 + 1.05*hrange) * mask1;
	}
	else{
	  float news2 = (MISCMATHS::Min((factor2*s2(x,y,z)),s2max)-s2min)/(s2max-s2min); 
	  int   mask2 = ((factor2*s2(x,y,z)) >= (s2min - latitude));
	  
	  if (checker > 0) 
	    mask2 *=(x+y+z)%2;
	  
	  if ( colour_type == 0 )
	    tmpvol.value(x,y,z) = A + newbg*hrange*0.95*(1-mask1)*(1-mask2) +
	      ( news1*hrange*0.9 + hrange*1.05 ) * mask1 +
	      ( news2*hrange*0.9 + hrange*2.05 ) * mask2;
	  else
	    tmpvol.value(x,y,z) = A + newbg*hrange*0.95*(1-mask1)*(1-mask2) +
	      ( ((int)(news1*5*0.999)) * 0.2*hrange +
		newbg*hrange*0.5*0.2 + 1.05*hrange) * mask1*(1-mask2) +
	      ( ((int)(news2*5*0.999)) * 0.2*hrange +
		newbg*hrange*0.5*0.2 + 2.05*hrange) * mask2*(1-mask1) +
	      ( ((int)((news1+news2)*5*0.999)) * 0.2*hrange +
		newbg*hrange*0.5*0.4 + 3.05*hrange) * mask1*mask2;
	}  
      }
    }
  }   
    
  if (ns == 1)
    lut = string("render1");
  else
    lut = string("render2");

  if (colour_type == 1)
    lut += string("t"); 
  tmpvol.setAuxFile(lut);
  
  lut = lutbase + lut + ".lut";

  copyconvert(tmpvol,newvol);
  // if(debug)
  //  cerr << " calling write_cbar() " << endl;
  if(cbarfname.length()>1)
    write_cbar(cbarfname, cbartype);

 //  cbarptr = gdImageCreateTrueColor(200,400);
//   if(cbarfname.length()>1)
//     write_cbar(cbarfname, cbartype);
//   gdImageDestroy(cbarptr);

  return(0);
}

/* }}} */

/* {{{ template class miscpic */

//template class miscpic<char>;
//template class miscpic<short>;
//template class miscpic<int>;
//template class miscpic<unsigned char>;
//template class miscpic<unsigned short>;
//template class miscpic<unsigned int>;
//template class miscpic<float>;
//template class miscpic<double>;

/* }}} */

}


