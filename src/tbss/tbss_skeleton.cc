/*  tbss_skeleton.cc

    Stephen Smith, FMRIB Analysis Group

    Copyright (C) 2005-2008 University of Oxford  */

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

// {{{ includes and options

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "utils/options.h"

using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;

string title="tbss_skeleton (Version 1.03)\nCopyright(c) 2005-2007, University of Oxford (Stephen Smith)";
string examples="tbss_skeleton -i <inputimage> -o <skeleton>\ntbss_skeleton -i <inputimage> -p <skel_thresh> <distancemap> <search_rule_mask> <4Ddata> <projected_4Ddata> [-a <alt_4D>] [-s <alt_skeleton>]}";

Option<string> inname(string("-i,--in"), string(""),
		  string("input image"),
		  true, requires_argument);
Option<string> outname(string("-o,--out"), string(""),
		  string("output image"),
		  false, requires_argument);
Option<string> projectargs(string("-p"), "",
			   string("~<skel_thresh> <distancemap> <search_rule_mask> <4Ddata> <projected_4Ddata>"), false, requires_5_arguments);
Option<string> alt4Dname(string("-a"), string(""),
		  string("alternative 4Ddata (e.g. L1)"),
		  false, requires_argument);
Option<string> altskelname(string("-s"), string(""),
		  string("alternative skeleton"),
		  false, requires_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<bool> debugging(string("-d,--debug"), false, 
		     string("switch on debugging image outputs"), 
		     false, no_argument);
Option<string> debugging2(string("-D,--debug2"), string(""),
		  string("~<skelpoints>\tde-project <skelpoints> points on skeleton back to all_FA space"),
		  false, requires_argument);

int nonoptarg;

// }}}

#define SEARCHSIGMA 10 /* length in linear voxel dimensions */
#define MAXSEARCHLENGTH (3*SEARCHSIGMA)

int main(int argc,char *argv[])
{
  // {{{ parse options

  Tracer tr("main");
  OptionParser options(title, examples);

  try {
    options.add(inname);
    options.add(outname);
    options.add(projectargs);
    options.add(alt4Dname);
    options.add(altskelname);
    options.add(help);
    options.add(debugging);
    options.add(debugging2);

    nonoptarg = options.parse_command_line(argc, argv);

    // line below stops the program if the help was requested or 
    //  a compulsory option was not set
    if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
      {
	options.usage();
	exit(EXIT_FAILURE);
      }
    
  }  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } 

// }}}

  volume<float> im;
  read_volume(im,inname.value());

  // {{{ estimate perp from CofG and curvature, and store this in X,Y,Z

volume<short> X(im.xsize(),im.ysize(),im.zsize()),
  Y(im.xsize(),im.ysize(),im.zsize()),
  Z(im.xsize(),im.ysize(),im.zsize());

X=0; Y=0; Z=0;

for(int z=1;z<im.zsize()-1;z++) for(int y=1;y<im.ysize()-1;y++) for(int x=1;x<im.xsize()-1;x++)
  {
    float theval = im(x,y,z);

    if ( theval != 0 )
      {
	float CofGx=0, CofGy=0, CofGz=0, Sum=0, CofGl;
	int xxx=0, yyy=0, zzz=0;

	for(int zz=-1; zz<=1; zz++) for(int yy=-1; yy<=1; yy++) for(int xx=-1; xx<=1; xx++)
	  {
	    float val = im(x+xx,y+yy,z+zz);
	    Sum   += val;
	    CofGx += xx * val;  CofGy += yy * val;  CofGz += zz * val;
	  }	      
	      
	CofGx /= Sum;  CofGy /= Sum;  CofGz /= Sum;
	CofGl = sqrt(CofGx*CofGx+CofGy*CofGy+CofGz*CofGz);
	      
	if (CofGl > .1)  /* is CofG far enough away from centre voxel? */
	  {
	    // if (CofGl > 1.4) cout << CofGx << " " << CofGy << " " << CofGz << " " << Sum << endl;
	    xxx = MISCMATHS::Max( MISCMATHS::Min( MISCMATHS::round(CofGx/CofGl) , 1 ) , -1);
	    yyy = MISCMATHS::Max( MISCMATHS::Min( MISCMATHS::round(CofGy/CofGl) , 1 ) , -1);
	    zzz = MISCMATHS::Max( MISCMATHS::Min( MISCMATHS::round(CofGz/CofGl) , 1 ) , -1);
	  }
	else
	  // {{{ find direction of max curvature

	{
	  float maxcost=0, centreval=2*theval;

	  for(int zz=0; zz<=1; zz++) /* note - starts at zero as we're only searching half the voxels */
	    for(int yy=-1; yy<=1; yy++)
	      for(int xx=-1; xx<=1; xx++)
		if ( (zz==1) || (yy==1) || ((yy==0)&&(xx==1)) ) /* only search half the voxels */
		  {
		    float weighting = pow( (float)(xx*xx+yy*yy+zz*zz) , -0.7 ); /* power is arbitrary: maybe test other functions here */
		    float cost = weighting * ( centreval 
					       - (float)im(x+xx,y+yy,z+zz)
					       - (float)im(x-xx,y-yy,z-zz) );

		    if (cost>maxcost)
		      {
			maxcost=cost;
			xxx=xx;
			yyy=yy;
			zzz=zz;
		      }
		  }
	}

// }}}
							   
	  X(x,y,z)=xxx;
	  Y(x,y,z)=yyy;
	  Z(x,y,z)=zzz;
      }

  }

// {{{ save perp image

if (debugging.set())
  {
    volume4D<float>tmpim(im.xsize(),im.ysize(),im.zsize(),3);
    copybasicproperties(im,tmpim);
    tmpim=0;
    
    for(int z=0;z<im.zsize();z++) for(int y=0;y<im.ysize();y++) for(int x=0;x<im.xsize();x++)
      {
	float tmpX=X(x,y,z), tmpY=Y(x,y,z), tmpZ=Z(x,y,z),
	  tmpf=sqrt(tmpX*tmpX+tmpY*tmpY+tmpZ*tmpZ);

	if (tmpf>0)
	  {
	    tmpim(x,y,z,0)=tmpX/tmpf;
	    tmpim(x,y,z,1)=tmpY/tmpf;
	    tmpim(x,y,z,2)=tmpZ/tmpf;
	  }
      }
    
    save_volume4D(tmpim,inname.value()+"_flow");
  }

// }}}

// }}}
  // {{{ smooth X,Y,Z and store in XX,YY,ZZ

volume<short> XX(im.xsize(),im.ysize(),im.zsize()),
  YY(im.xsize(),im.ysize(),im.zsize()),
  ZZ(im.xsize(),im.ysize(),im.zsize());

XX=0; YY=0; ZZ=0;

for(int z=1;z<im.zsize()-1;z++) for(int y=1;y<im.ysize()-1;y++) for(int x=1;x<im.xsize()-1;x++)
  {
    int* localsum = new int[27];
    int localmax=0, xxx, yyy, zzz;

    for(int zz=0; zz<27; zz++) localsum[zz]=0;

    for(int zz=-1; zz<=1; zz++) for(int yy=-1; yy<=1; yy++) for(int xx=-1; xx<=1; xx++)
      {
	xxx = X(x+xx,y+yy,z+zz);
	yyy = Y(x+xx,y+yy,z+zz);
	zzz = Z(x+xx,y+yy,z+zz);
	// cout << xxx << " " << yyy << " " << zzz << " " << (1+zzz)*9+(1+yyy)*3+1+xxx << endl;
	localsum[(1+zzz)*9+(1+yyy)*3+1+xxx]++;
	localsum[(1-zzz)*9+(1-yyy)*3+1-xxx]++;
      }
    
    for(int zz=-1; zz<=1; zz++) for(int yy=-1; yy<=1; yy++) for(int xx=-1; xx<=1; xx++)
      {
	if (localsum[(1+zz)*9+(1+yy)*3+1+xx]>localmax)
	  {
	    localmax=localsum[(1+zz)*9+(1+yy)*3+1+xx];
	    XX(x,y,z)=xx;
	    YY(x,y,z)=yy;
	    ZZ(x,y,z)=zz;
	  }
      }

    delete localsum;
  }

X.destroy();
Y.destroy();
Z.destroy();

// {{{ save perp image

if (debugging.set())
  {
    volume4D<float>tmpim(im.xsize(),im.ysize(),im.zsize(),3);
    copybasicproperties(im,tmpim);
    tmpim=0;
    
    for(int z=0;z<im.zsize();z++) for(int y=0;y<im.ysize();y++) for(int x=0;x<im.xsize();x++)
      {
	float tmpX=XX(x,y,z), tmpY=YY(x,y,z), tmpZ=ZZ(x,y,z),
	  tmpf=sqrt(tmpX*tmpX+tmpY*tmpY+tmpZ*tmpZ);

	if (tmpf>0)
	  {
	    tmpim(x,y,z,0)=tmpX/tmpf;
	    tmpim(x,y,z,1)=tmpY/tmpf;
	    tmpim(x,y,z,2)=tmpZ/tmpf;
	  }
      }
    
    save_volume4D(tmpim,inname.value()+"_flowsmooth");
  }

// }}}

// }}}
  // {{{ do non-max-suppression in the direction of perp and save to file

volume<float> tmpim(im);
tmpim=0;

for(int z=1;z<im.zsize()-1;z++) for(int y=1;y<im.ysize()-1;y++) for(int x=1;x<im.xsize()-1;x++)
  {
    float theval = im(x,y,z);
    int xxx=XX(x,y,z);
    int yyy=YY(x,y,z);
    int zzz=ZZ(x,y,z);
    
    if ( ( (xxx!=0) || (yyy!=0) || (zzz!=0) ) &&
	 ( theval >= im(x+xxx,y+yyy,z+zzz) ) &&
	 ( theval >  im(x-xxx,y-yyy,z-zzz) ) &&
	 ( theval >= im(x+2*xxx,y+2*yyy,z+2*zzz) ) &&
	 ( theval >  im(x-2*xxx,y-2*yyy,z-2*zzz) ) )
      tmpim(x,y,z) = theval;
  }

if (outname.set())
  save_volume(tmpim,outname.value());

// }}}

  if (projectargs.set())
    // {{{ do search in the direction of perp

{
  float origthresh = atof(projectargs.value(0).c_str());

  volume<float> distancemap;
  read_volume(distancemap,projectargs.value(1));

  volume<int> lowercingulum;
  read_volume(lowercingulum,projectargs.value(2));

  volume4D<float> data_4d;
  read_volume4D(data_4d,projectargs.value(3));

  volume4D<float> alt_data_4d;
  if (alt4Dname.set())
    read_volume4D(alt_data_4d,alt4Dname.value());

  if (altskelname.set())
    read_volume(tmpim,altskelname.value());

  volume4D<float> data_4d_projected(data_4d);
  data_4d_projected=0;

  // {{{ debugging

  volume4D<short> tmpimFLOWx, tmpimFLOWy, tmpimFLOWz;
  if (debugging.set())
    {
      tmpimFLOWx.reinitialize(im.xsize(),im.ysize(),im.zsize(),data_4d.tsize());
      tmpimFLOWy.reinitialize(im.xsize(),im.ysize(),im.zsize(),data_4d.tsize());
      tmpimFLOWz.reinitialize(im.xsize(),im.ysize(),im.zsize(),data_4d.tsize());
      copybasicproperties(im,tmpimFLOWx);
      copybasicproperties(im,tmpimFLOWy);
      copybasicproperties(im,tmpimFLOWz);
      tmpimFLOWx=0;
      tmpimFLOWy=0;
      tmpimFLOWz=0;
    }

  volume<float> debug2in;
  volume4D<float> debug2out;
  if (debugging2.set())
    {
      read_volume(debug2in,debugging2.value());
      debug2out.reinitialize(im.xsize(),im.ysize(),im.zsize(),data_4d.tsize());
      copybasicproperties(im,debug2out);
      debug2out=0;
    }

// }}}

  for(int T=0;T<data_4d.tsize();T++) for(int z=1;z<im.zsize()-1;z++) for(int y=1;y<im.ysize()-1;y++) for(int x=1;x<im.xsize()-1;x++)
    if (tmpim(x,y,z) > origthresh)
      {
	int xxx=XX(x,y,z), yyy=YY(x,y,z), zzz=ZZ(x,y,z);
	short maxvalX=0, maxvalY=0, maxvalZ=0;
	float maxval=data_4d(x,y,z,T), maxval_weighted=maxval,
	  exponentfactor = -0.5 * (xxx*xxx+yyy*yyy+zzz*zzz) / (float)(SEARCHSIGMA*SEARCHSIGMA);
	if (alt4Dname.set()) maxval=alt_data_4d(x,y,z,T);

	if (lowercingulum(x,y,z) == 0)
	  // {{{ search perp to sheet

{
  for(int iters=0;iters<2;iters++)
    {
      float distance=0;

      for(int d=1;d<MAXSEARCHLENGTH;d++)
	{
	  int D=d;
	  if (iters==1) D=-d;
		    
	  if (distancemap(x+xxx*D,y+yyy*D,z+zzz*D)>=distance)
	    {
	      float distanceweight = exp(d * d * exponentfactor);
	      distance=distancemap(x+xxx*D,y+yyy*D,z+zzz*D);
	      if (distanceweight * data_4d(x+xxx*D,y+yyy*D,z+zzz*D,T)>maxval_weighted)
		{
		  maxval=data_4d(x+xxx*D,y+yyy*D,z+zzz*D,T);
		  maxval_weighted=maxval*distanceweight;
		  maxvalX=xxx*D;
		  maxvalY=yyy*D;
		  maxvalZ=zzz*D;
		  if (alt4Dname.set()) maxval=alt_data_4d(x+xxx*D,y+yyy*D,z+zzz*D,T);
		}
	    }
	  else
	    d=MAXSEARCHLENGTH;
	}
    }
}

// }}}
	else
	  // {{{ search all around tube

{
  for(int yyy=-MAXSEARCHLENGTH; yyy<=MAXSEARCHLENGTH; yyy++) for(int xxx=-MAXSEARCHLENGTH; xxx<=MAXSEARCHLENGTH; xxx++) 
    {
      float distanceweight = exp(-0.5 * (xxx*xxx+yyy*yyy) / (float)(SEARCHSIGMA*SEARCHSIGMA) );

      float r=sqrt((float)(xxx*xxx+yyy*yyy));

      if (r>0)
	{
	  int allok=1;
	  
	  for(float rr=1; rr<=r+0.1; rr++) /* search outwards from centre to current voxel - test that distancemap always increasing */
	    {
	      int xxx1=MISCMATHS::round(rr*xxx/r);
	      int yyy1=MISCMATHS::round(rr*yyy/r);
	      int xxx2=MISCMATHS::round((rr+1)*xxx/r);
	      int yyy2=MISCMATHS::round((rr+1)*yyy/r);
	      if ( distancemap(x+xxx1,y+yyy1,z) > distancemap(x+xxx2,y+yyy2,z) )
		allok=0;
	    }

	  if ( allok &&
	       ( distanceweight * data_4d(x+xxx,y+yyy,z,T) > maxval_weighted ) )
	    {
	      maxval=data_4d(x+xxx,y+yyy,z,T);
	      maxval_weighted=maxval*distanceweight;
	      maxvalX=xxx;
	      maxvalY=yyy;
	      maxvalZ=0;
	      if (alt4Dname.set()) maxval=alt_data_4d(x+xxx,y+yyy,z,T);
	    }
	}
    }

}

// }}}
	
	  data_4d_projected(x,y,z,T)=maxval; /* output maxsearch data */

	  // {{{ debugging search

	  if (debugging.set())
	    {
	      tmpimFLOWx(x,y,z,T)=maxvalX;
	      tmpimFLOWy(x,y,z,T)=maxvalY;
	      tmpimFLOWz(x,y,z,T)=maxvalZ;
	    }

	  if (debugging2.set())
	    if (debug2in(x,y,z)>0)
	      debug2out(x+maxvalX,y+maxvalY,z+maxvalZ,T)=debug2in(x,y,z);

// }}}
      }

  data_4d.destroy();
  alt_data_4d.destroy();
  save_volume4D(data_4d_projected,projectargs.value(4));

  // {{{ debugging search

  if (debugging.set())
    {
      save_volume4D(tmpimFLOWx,projectargs.value(4)+"_search_X");
      save_volume4D(tmpimFLOWy,projectargs.value(4)+"_search_Y");
      save_volume4D(tmpimFLOWz,projectargs.value(4)+"_search_Z");
    }

  if (debugging2.set())
    save_volume4D(debug2out,projectargs.value(4)+"_deprojected");

// }}}
}

// }}}
}

