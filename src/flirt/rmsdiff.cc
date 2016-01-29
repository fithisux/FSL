/*  rmsdiff.cc

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
#include <unistd.h>
#define WANT_STREAM
#define WANT_MATH

#include "newmatap.h"
#include "newmatio.h"
#include "miscmaths/miscmaths.h"
#include "newimage/newimageall.h"

#ifndef NO_NAMESPACE
 using namespace MISCMATHS;
 using namespace NEWMAT;
 using namespace NEWIMAGE;
#endif


////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{

  float rmax=80.0;

  if (argc<4) { 
    cerr << "Usage: " << argv[0] << " matrixfile1 matrixfile2 refvol [mask]" << endl; 
    cerr << "        Outputs rms deviation between matrices (in mm)" << endl;
    return -1; 
  }
  
  Matrix affmat1(4,4), affmat2(4,4);
  affmat1 = read_ascii_matrix(argv[1]);
  if (affmat1.Nrows()<4) {
    cerr << "Could not read matrix " << argv[1] << endl;
    return -2;
  }
  affmat2 = read_ascii_matrix(argv[2]);
  if (affmat2.Nrows()<4) {
    cerr << "Could not read matrix " << argv[2] << endl;
    return -2;
  }

  if (fabs(affmat1.Determinant())<0.1) {
    cerr << "WARNING:: matrix 1 has low determinant" << endl;
    cerr << affmat1 << endl;
  }

  if (fabs(affmat2.Determinant())<0.1) {
    cerr << "WARNING:: matrix 2 has low determinant" << endl;
    cerr << affmat2 << endl;
  }

  ColumnVector centre(3);
  centre = 0;

  volume<float> refvol;
  read_volume(refvol,argv[3]);

  if (argc<5) {
    // do the RMS
    // compute the centre of gravity
    centre = refvol.cog("scaledmm");
    float rms = rms_deviation(affmat1,affmat2,centre,rmax);
    cout << rms << endl;
  } else {  
    // do the extreme distance
    double maxdist=0, dist=0, sumdistsq=0;
    ColumnVector cvec(4);
    cvec=0;  cvec(4)=1;
    long int nvox=0;
    volume<float> mask;
    read_volume(mask,argv[4]);
    for (int z=mask.minz(); z<=mask.maxz(); z++) {
      for (int y=mask.miny(); y<=mask.maxy(); y++) {
	for (int x=mask.minx(); x<=mask.maxx(); x++) {
	  if (mask(x,y,z)>0.5) {
	    cvec(1)=x*refvol.xdim();  cvec(2)=y*refvol.ydim();  cvec(3)=z*refvol.zdim();
	    dist = norm2((affmat1 -affmat2)*cvec);
	    maxdist=Max(dist,maxdist);
	    sumdistsq+=dist*dist;
	    nvox++;
	  }
	}
      }
    }
    cout << maxdist << endl;
    double rms = sqrt(sumdistsq/nvox);
    cout << rms << endl;
  }

  return 0;

}





