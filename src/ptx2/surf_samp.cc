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
#ifndef EXPOSE_TREACHEROUS
#define EXPOSE_TREACHEROUS
#endif

#include "utils/options.h"
#include "newimage/newimageall.h"
#include "csv.h"
#include "stdlib.h"
#include "string.h"
#include "miscmaths/miscmaths.h"
#include "warpfns/warpfns.h"
#include "warpfns/fnirt_file_reader.h"



using namespace Utilities;
using namespace std;
using namespace NEWIMAGE;
using namespace MISCMATHS;


string title="";
string examples="";


Option<string> datafile(string("--data"),string(""),
		   string("data"),
		   true,requires_argument);
Option<string> coordfile(string("--coordfile"),string(""),
		   string("coord file (list of vertex indices, start at 0)"),
		   true,requires_argument);
Option<string> surf(string("--surf"),string(""),
		    string("surface file"),
		    true,requires_argument);
Option<string> surfvol(string("--meshref"),string(""),
		    string("surface volume ref"),
		    true,requires_argument);
Option<string> xfm(string("--xfm"),string(""),
		   string("surf2diff transform"),
		   true,requires_argument);
Option<string> meshspace(string("--meshspace"),string("freesurfer"),
		   string("meshspace. Either 'freesurfer' [default], caret, first, vox "),
		   true,requires_argument);
Option<string> out(string("--out"),string(""),
		   string("output file"),
		   true,requires_argument);
Option<float> dist(string("--step"),1,
		   string("step (mm - default=1)"),
		   false,requires_argument);
Option<int>  dir(string("--dir"),1,
		  string("direction: 1=inwards(default), -1=outward"),
		  false,requires_argument);
Option<int> nsteps(string("--nsteps"),10,
		  string("number of steps (default=10)"),
		  false,requires_argument);


int main(int argc,char *argv[]){

  OptionParser options(title,examples);
  
  
  options.add(datafile);
  options.add(coordfile);
  options.add(surf);
  options.add(surfvol);
  options.add(xfm);
  options.add(meshspace);
  options.add(out);
  options.add(dist);
  options.add(dir);
  options.add(nsteps);

  
  options.parse_command_line(argc,argv);
  if(!options.check_compulsory_arguments(true)){
    options.usage();
    return(1);
  }
  
  cout<<"read"<<endl;
  volume<float> data;
  read_volume(data,datafile.value());
  volume<short int> refvol;
  read_volume(refvol,surfvol.value());

  CSV csv(refvol);
  csv.set_convention(meshspace.value());
  csv.load_rois(surf.value());

  Matrix list = read_ascii_matrix(coordfile.value());
  if(list.Nrows()==1){list=list.t();}

  ColumnVector surfdim(3),datadim(3);
  surfdim << refvol.xdim() << refvol.ydim() << refvol.zdim();
  datadim << data.xdim()<<data.ydim()<<data.zdim();

  bool isWarp=false;
  volume4D<float> data2surf_warp;
  Matrix          surf2data,data2surf;
  if(xfm.value() != ""){
    if(fsl_imageexists(xfm.value())){
      isWarp=true;
      FnirtFileReader ffr(xfm.value());
      data2surf_warp = ffr.FieldAsNewimageVolume4D(true);
    }
    else{
      data2surf = read_ascii_matrix(xfm.value());
      surf2data = data2surf.i();
    }
  }
  else{
    surf2data = IdentityMatrix(4);
  }


  Matrix surfdata(list.Nrows(),nsteps.value());
  surfdata=0;

  float step = dist.value(); // mm
  
  ColumnVector x(3),n(3),data_vox(3),nn(3);
  int ix,iy,iz;
  for(int ind=1;ind<=list.Nrows();ind++){
    //  for(int i=0;i<csv.get_mesh(0).nvertices();i++){
    int i=(int)list(ind,1);

    x=csv.get_vertex_as_vox(0,i);
    n=csv.get_normal_as_vox(0,i);
    if(dir.value()<0) n*=-1;

    if(n.MaximumAbsoluteValue()>0)
      n/=sqrt(n.SumSquare());

    ColumnVector xx(3);
    vector< vector<float> > vals(nsteps.value());
    for(int s=0;s<nsteps.value();s++){
      xx<<x(1)+(float)s/float(nsteps.value()-1.0)*step/surfdim(1)*n(1)
	<<x(2)+(float)s/float(nsteps.value()-1.0)*step/surfdim(2)*n(2)
	<<x(3)+(float)s/float(nsteps.value()-1.0)*step/surfdim(3)*n(3);
      if(!isWarp)
	data_vox = vox_to_vox(xx,surfdim,datadim,surf2data);
      else
	data_vox = NewimageCoord2NewimageCoord(data2surf_warp,false,refvol,data,xx);
	
      ix = (int)round((float)data_vox(1));
      iy = (int)round((float)data_vox(2));
      iz = (int)round((float)data_vox(3));
      
      
      surfdata(ind,s+1)=data(ix,iy,iz);
    }
  }
    

  cout<<"save"<<endl;

  write_ascii_matrix(surfdata,out.value());



  return 0;
}



