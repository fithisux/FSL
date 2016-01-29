/*    Copyright (C) 2012 University of Oxford  */

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
#include "utils/options.h"
#include "newimage/newimageall.h"
#include "csv.h"
#include "stdlib.h"
#include "string.h"
#include "miscmaths/miscmaths.h"



using namespace Utilities;
using namespace std;
using namespace NEWIMAGE;
using namespace MISCMATHS;


string title="";
string examples="";


Option<string> surf(string("--surf"),string(""),
		    string("surface file"),
		    true,requires_argument);
Option<string> surfvolref(string("--meshref"),string(""),
			  string("surface volume ref"),
			  true,requires_argument);
Option<string> outvolref(string("--outref"),string(""),
			 string("output volume ref"),
			 true,requires_argument);
Option<string> xfm(string("--xfm"),string(""),
		   string("surf2out transform"),
		   true,requires_argument);
Option<string> meshspace(string("--meshspace"),string("freesurfer"),
			 string("meshspace (default='freesurfer')"),
			 true,requires_argument);
Option<string> out(string("--out"),string(""),
		   string("output file"),
		   true,requires_argument);
Option<float> dist(string("--step"),1,
		   string("step, >0 means towards brain (mm - default=1)"),
		   false,requires_argument);


int main(int argc,char *argv[]){

  OptionParser options(title,examples);
  
  
  options.add(surf);
  options.add(surfvolref);  
  options.add(outvolref);
  options.add(xfm);
  options.add(meshspace);
  options.add(out);
  options.add(dist);
  
  options.parse_command_line(argc,argv);
  if(!options.check_compulsory_arguments(true)){
    options.usage();
    return(1);
  }
  
  cout<<"read"<<endl;
  volume<short int> surfvol;
  read_volume(surfvol,surfvolref.value());

  CSV csv(surfvol);
  csv.set_convention(meshspace.value());
  csv.load_rois(surf.value());

  volume4D<float> dyads;
  volume<float> outvol;
  read_volume(outvol,outvolref.value());
  dyads.reinitialize(outvol.xsize(),
		     outvol.ysize(),
		     outvol.zsize(),3);
  copybasicproperties(outvol,dyads);


  Matrix surf2out = read_ascii_matrix(xfm.value());
  Matrix ixfm(3,3);
  ixfm=surf2out.SubMatrix(1,3,1,3);
  ixfm=ixfm.i();

  ColumnVector surfdim(3),outdim(3);
  surfdim << surfvol.xdim() << surfvol.ydim() << surfvol.zdim();
  outdim << dyads[0].xdim()<<dyads[0].ydim()<<dyads[0].zdim();

  float step = dist.value(); // mm

  cout<<"loop"<<endl;
  // loop over surface vertices
  // for each vertex, go along normal for a fixed distance, and sample orientations
  // save the orientation, the f1/f2, etc.
  ColumnVector x(3),xx(3),n(3),dti_vox(3),nn(3);
  int ix,iy,iz;
  int nsteps=10;
  dyads=0;
  for(int i=0;i<csv.get_mesh(0).nvertices();i++){
    x=csv.get_vertex_as_vox(0,i);
    n=csv.get_normal_as_vox(0,i);

    // rotate orientation
    n = surf2out.SubMatrix(1,3,1,3)*n;

    // flip normal in caret
    if(meshspace.value()=="caret")
      n*=-1;
    if(n.MaximumAbsoluteValue()>0)
      n/=sqrt(n.SumSquare());
    
    for(int s=0;s<nsteps;s++){
      xx<<x(1)+s/(nsteps-1)*step/surfdim(1)*n(1)
	<<x(2)+s/(nsteps-1)*step/surfdim(2)*n(2)
	<<x(3)+s/(nsteps-1)*step/surfdim(3)*n(3);
      
      dti_vox = vox_to_vox(x,surfdim,outdim,surf2out);
      ix = (int)round((float)dti_vox(1));
      iy = (int)round((float)dti_vox(2));
      iz = (int)round((float)dti_vox(3));
      
      dyads(ix,iy,iz,0) = n(1);
      dyads(ix,iy,iz,1) = n(2);
      dyads(ix,iy,iz,2) = n(3);
    }
    
  }
  cout<<"save"<<endl;

  save_volume4D(dyads,out.value());

  return 0;
}



