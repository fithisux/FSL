/*  avscale.cc

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
#define WANT_STREAM
#define WANT_MATH

#include "newmatap.h"
#include "newmatio.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"

#ifndef NO_NAMESPACE
 using namespace MISCMATHS;
 using namespace NEWMAT;
 using namespace NEWIMAGE;
#endif


////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{

  try {

    if (argc<2) { 
      cerr << "Usage: " << argv[0] << " [--allparams/--inverteddies] matrixfile [non-reference-volume]\n"; 
      return -1; 
    }

    int offset=1;
    string optionarg=argv[1];
    if (optionarg[0]!='-') {  optionarg=""; offset=0; }

    Matrix affmat(4,4);
    ColumnVector params(12), cor(3);
    cor=0;  // centre of rotations

    if (argc>2 + offset) {
      string volname=argv[2 + offset];
      volume<float> testvol;
      if (volname!="") {
        //if (read_volume_hdr_only(testvol,argv[2 + offset])<0)  return -1;
	// NB: need to read all of volume to get cog() - may not need always
        if (read_volume(testvol,argv[2 + offset])<0)  return -1;
      }
      affmat = read_ascii_matrix(argv[1 + offset]);
      if (affmat.Nrows()<4)   return -2;
      cor = testvol.cog("scaled_mm");
      // the following is for when I get around to offering a cov option
      //cor(1) = (testvol.xsize() - 1.0) * testvol.xdim() / 2.0;
      //cor(2) = (testvol.ysize() - 1.0) * testvol.ydim() / 2.0;
      //cor(3) = (testvol.zsize() - 1.0) * testvol.zdim() / 2.0;
    } else {
      affmat = read_ascii_matrix(argv[1 + offset]);
      if (affmat.Nrows()<4) return -2;
    }

    //cout << affmat << endl;
    decompose_aff(params,affmat,cor,rotmat2euler);

    Matrix rotmat(4,4);
    construct_rotmat_euler(params,6,rotmat);
    
    cout << "Rotation & Translation Matrix:\n" << rotmat << endl;
    if (optionarg=="--allparams") {
      cout << "Rotation Angles (x,y,z) [rads] = " << params.SubMatrix(1,3,1,1).t() << endl;
      cout << "Translations (x,y,z) [mm] = " << params.SubMatrix(4,6,1,1).t() << endl;
    }

    cout << "Scales (x,y,z) = " << params.SubMatrix(7,9,1,1).t() << endl;
    cout << "Skews (xy,xz,yz) = " << params.SubMatrix(10,12,1,1).t() << endl;
    float avscale = (params(7) + params(8) + params(9))/3.0;
    cout << "Average scaling = " << avscale << endl << endl;

    cout << "Determinant = " << affmat.Determinant() << endl;
    cout << "Left-Right orientation: ";
    if (affmat.Determinant()>0) { cout << "preserved" << endl; }
    else { cout << "swapped" << endl; }
    cout << endl;

    Matrix swapmat;
    swapmat = IdentityMatrix(4);
    if (affmat.Determinant()<0) swapmat(1,1) = -1;
    Matrix m2 = sqrtaff(affmat * swapmat);
    Matrix m0 = m2*affmat.i();
    cout << "Forward half transform =\n" << m2 << endl;
    cout << "Backward half transform =\n" << m0 << endl;

    Matrix scale,skew;
    scale=IdentityMatrix(4);
    skew=IdentityMatrix(4);
    scale(1,1) = params(7);
    scale(2,2) = params(8);
    scale(3,3) = params(9);
    skew(1,2) = params(10);
    skew(1,3) = params(11);
    skew(2,3) = params(12);
    Matrix ans;
    ans = affmat.i() * rotmat * skew * scale;

    if (optionarg=="--inverteddies") {
      cout << endl << "Matrix check: mat^-1 * rotmat * skew * scale =\n";
      cout << ans << endl;

      // reconstitute scale and skew matrices with the inverse effect in y
      scale(1,1) = params(7);
      scale(2,2) = 1.0/params(8);
      scale(3,3) = params(9);
      skew(1,2) = -params(10);
      skew(1,3) = params(11);
      skew(2,3) = -params(12);
      Matrix newxfm;
      newxfm = rotmat * skew * scale;
      cout << "Inverted eddy matrix:" << endl << newxfm << endl;
    }
    

    return 0;
  }
  catch(Exception exc) {
    cerr << exc.what() << endl;
    throw;
  }
  catch(...) {
    cerr << "Image error" << endl;
    throw;
  } 
  return(0);
}








