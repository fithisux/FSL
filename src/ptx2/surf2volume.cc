/*  Copyright (C) 2004 University of Oxford  */

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


#include "csv.h"
#include "newimage/newimageio.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include "miscmaths/SpMat.h"
#include "meshclass/meshclass.h"

using namespace NEWIMAGE;
using namespace MISCMATHS;
using namespace mesh;


int main(int argc, char** argv){

  if(argc<2){
    cout<<"surf2vol <surf> <refvol> <outvol> <convention>"<<endl;
    exit(1);
  }

  volume<short int> refvol;
  read_volume(refvol,argv[2]);

  CSV csv(refvol);
  csv.set_convention(argv[4]);
  csv.load_rois(argv[1]);
  csv.save_surfvol(argv[3],false);
  //  csv.save_normalsAsVol(0,string(argv[3])+"_normal");


//   // Save surface normal
//   volume4D<float> odir(refvol.xsize(),
// 		       refvol.ysize(),
// 		       refvol.zsize(),
// 		       3);
//   copybasicproperties(refvol,odir);
//   odir=0;

//   ColumnVector pos(3),dir(3);
//   for(int i=0;i<csv.get_mesh(0).nvertices();i++){
//     pos=csv.get_vertex_as_vox(0,i);
//     dir=csv.get_normal_as_vox(0,i);
    
//     odir((int)round((float)pos(1)),
// 	 (int)round((float)pos(2)),
// 	 (int)round((float)pos(3)),0)=dir(1);
//     odir((int)round((float)pos(1)),
// 	 (int)round((float)pos(2)),
// 	 (int)round((float)pos(3)),1)=dir(2);
//     odir((int)round((float)pos(1)),
// 	 (int)round((float)pos(2)),
// 	 (int)round((float)pos(3)),2)=dir(3);
	 
//   }
//   odir.setDisplayMaximumMinimum(1,-1);
//   save_volume4D(odir,string(argv[3])+"_normal");

//   Matrix n(csv.get_mesh(0).nvertices(),3);
//   for(int i=0;i<csv.get_mesh(0).nvertices();i++){
//     dir=csv.get_normal(0,i);
//     n.Row(i+1)<<dir(1)<<dir(2)<<dir(3);
//   }
//   write_ascii_matrix(n,"normalsCSV.txt");

  //csv.save_roi(0,"surf.txt");

  return 0;
}
