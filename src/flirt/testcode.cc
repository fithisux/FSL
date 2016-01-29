/*  testcode.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

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
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <time.h>
#include <vector>
#include <algorithm>
#define WANT_STREAM
#define WANT_MATH

#include "newmatap.h"
#include "newmatio.h"
//#include "miscimfns.h"
#include "miscmaths/miscmaths.h"
//#include "interpolation.h"
//#include "mjimage.h"
#include "newimage/costfns.h"
#include "newimage/newimageall.h"
#include "globaloptions.h"

#ifndef NO_NAMESPACE
 using namespace MISCMATHS;
 using namespace MISCIMFNS;
 using namespace COSTFNS;
 using namespace INTERPOLATION;
 using namespace NEWMAT;
 using namespace MJIMAGE;
 using namespace GENERALIO;
#endif



int testfn_yy(int argc, char *argv[]) {
  volume gfimage;
  if (argc > 1) {
    read_volume(gfimage,argv[1]);
  } else {
    read_volume(gfimage,"gfimage.hdr");
  }
	cout << "Read File Successfully!" << endl;
	cout << "Width: " << gfimage.rows() << "\tHeight: " 
	     << gfimage.columns() << "\tDepth: " << gfimage.slices() << endl; 

	//ZEnumPixelOf<float>	pixel(&gfimage);

	//for(pixel.Reset(); pixel.FMore(); pixel.Next()) *pixel *= -0.5;

	float max=0.0, min=0.0;
	get_min_max(gfimage,min,max);
	cout << "Min and Max = " << min << " & " << max << endl;

	for (float val=0.3; val>-0.005; val-=0.0014) {
	  for (int z=0; z<gfimage.zsize(); z++) {
	    for (int y=0; y<gfimage.ysize(); y++) {
	      for (int x=0; x<gfimage.xsize(); x++) {
		gfimage(x,y,z) *= (1.0+val);
	      }
	    }
	  }
	  cerr << ".";
	}

	get_min_max(gfimage,min,max);
	cout << "Min and Max = " << min << " & " << max << endl;


	save_volume(gfimage,"result");

  return -1;
}



int main(void)
{


  //  if (testfn_yy(argc,argv)<0) exit(-1);


  volume ref, test, bindex;
  read_volume(ref,"/usr/people/flitney/src/cr/ref");
  ref.setvoxelsize(1.0,1.0,1.0);
  ref.setvoxelorigin(0,0,0);
  read_volume(test,"/usr/people/flitney/src/cr/test");
  test.setvoxelsize(1.0,1.0,1.0);
  test.setvoxelorigin(0,0,0);
  read_volume(bindex,"/usr/people/flitney/src/cr/ref");
  bindex.setvoxelsize(1.0,1.0,1.0);
  bindex.setvoxelorigin(0,0,0);

  float theta=0.0, thecr=0.0;
  Matrix xfm;

  imagepair imp(ref,test);
  imp.set_no_bins(256);
  globaloptions::get().impair = &imp;

  for(theta = 0; theta < 360.0; theta += 4.0)
    {
      Matrix tr1, tr2, rot(4,4);
      ColumnVector angl(3), trans(3);
      trans = 0;
      angl = 0;
      angl(1) = theta*M_PI/180.0;

      make_rot(angl,trans,rot);

      tr1=IdentityMatrix(4);
      tr1(1,4) = -(test.xsize())/2.0;
      tr1(2,4) = -(test.ysize())/2.0;
      tr1(3,4) = -(test.zsize())/2.0;
      tr2=IdentityMatrix(4);
      tr2(1,4) = (test.xsize())/2.0;
      tr2(2,4) = (test.ysize())/2.0;
      tr2(3,4) = (test.zsize())/2.0;
      
      xfm = tr2 * rot * tr1;

      thecr = corr_ratio(globaloptions::get().impair, xfm);
      cerr << 1.0 - thecr << endl;
    }

  return 0;
}
