
/*  Copyright (C) 2007 University of Oxford  */

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

#include "ptx_meshmask.h"
#include "streamlines.h"
using namespace std;
using namespace NEWIMAGE;
using namespace TRACT;
using namespace Utilities;
using namespace PARTICLE;
using namespace mesh;



void meshmask()
{ 
  probtrackxOptions& opts =probtrackxOptions::getInstance();

  // load seed mesh
  cout<<"loading mesh files"<<endl;
  Mesh mseeds;
  mseeds.load(opts.meshfile.value());
  if(opts.meshspace.value()!="first")
    mseeds.load_fs_label(opts.seedfile.value());
  cout<<"mesh files loaded"<<endl;

  // internally create seed mask in voxel space
  volume<float> seeds;
  if(opts.seedref.value()!="")
    read_volume(seeds,opts.seedref.value());
  else
    read_volume(seeds,opts.maskfile.value());
  seeds=0;


  
  Matrix mm_to_vox(4,4);
  if(opts.meshspace.value()=="freesurfer"){
    mm_to_vox << -1 << 0 << 0 <<  (seeds.xsize())/2
	      <<  0 << 0 << -1 << (seeds.zsize())/2
	      <<  0 << 1 << 0 <<  (seeds.ysize())/2
	      <<  0 << 0 << 0 << 1;
  }
  else if(opts.meshspace.value()=="caret"){
    cerr<<"caret not supported at the moment"<<endl;
    exit(1);
    //mm_to_vox = seeds.sform_mat();
    //mm_to_vox = mm_to_vox.i();
  }
  else if(opts.meshspace.value()=="first"){
    mm_to_vox << 1.0/seeds.xdim() << 0 << 0 << 0
	      << 0 << 1.0/seeds.ydim() << 0 << 0
	      << 0 << 0 << 1.0/seeds.zdim() << 0
	      << 0 << 0 << 0 << 1;
  }
  else{
    cerr << "unkown mesh option - exit without doing anything" << endl;
    cerr << "mesh option is either 'freesurfer' or 'caret'" << endl;
    cerr << "this is because each surface vertex coordinates are relative to a specific " << endl;
    cerr << "coordinate system which needs to be transformed into voxel space" << endl;
    exit(1);
  }


  ColumnVector fs_coord_mm(4),xyz_vox,seeddim(3);
  seeddim << seeds.xdim() << seeds.ydim() << seeds.zdim();
  ColumnVector dir(3);
  int keeptotal=0;

  // first fill seed with ones
  seeds = 0;
  int numseeds=0;
  for(vector<Mpoint*>::iterator i = mseeds._points.begin();i!=mseeds._points.end();i++){
    if((*i)->get_value() > 0 || opts.meshspace.value()=="first"){
    
      fs_coord_mm<<(*i)->get_coord().X<<(*i)->get_coord().Y<<(*i)->get_coord().Z << 1.0; 
      xyz_vox = mm_to_vox*fs_coord_mm;


      
      float x=xyz_vox(1);float y=xyz_vox(2);float z=xyz_vox(3);
      Pt newPt(x,y,z);
      (*i)->_update_coord = newPt;

      seeds(int(round(x)),int(round(y)),int(round(z))) = 1;
      numseeds++;
    }
  }

  
  ////////////////////////////////
  //  Log& logger = LogSingleton::getInstance();
  Streamliner stline(seeds);
  Counter counter(seeds,stline,numseeds);
  counter.initialise();
  Seedmanager seedmanager(counter);


  for(vector<Mpoint*>::iterator i = mseeds._points.begin();i!=mseeds._points.end();i++){
    if((*i)->get_value() > 0 || opts.meshspace.value()=="first"){
    
      fs_coord_mm<<(*i)->get_coord().X<<(*i)->get_coord().Y<<(*i)->get_coord().Z << 1.0; 
      //      xyz_vox = seeds.qform_mat().i()*fs_coord_mm
            xyz_vox = mm_to_vox*fs_coord_mm;
	    //xyz_vox = fs_coord_mm;
      

      float x=xyz_vox(1);float y=xyz_vox(2);float z=xyz_vox(3);
      Pt newPt(x,y,z);
      (*i)->_update_coord = newPt;

		//      seeds(round(x),round(y),round(z)) = 1;
    
      cout <<"run"<<endl;
      dir << (*i)->local_normal().X << (*i)->local_normal().Y << (*i)->local_normal().Z;
 
      if(opts.meshspace.value()=="first")
	dir*=-1.0;

     keeptotal += seedmanager.run(x,y,z,true,-1,dir); 
	
    }
  }
  mseeds.update();
  //  mseeds.save("test.vtk",3);

  //return;
  
  //   for(int z=0;z<seeds.zsize();z++){
  //     cout <<"sl "<<z<<endl;
  //     for(int y=0;y<seeds.ysize();y++){
  //       for(int x=0;x<seeds.xsize();x++){
  // 	if(seeds(x,y,z)>0){
  // 	  cout <<"run"<<endl;
  // 	  dir << (*i)->local_normal().X << (*i)->local_normal().Y << (*i)->local_normal().Z;
  // 	  keeptotal += seedmanager.run(x,y,z,true,-1,dir); 
  // 	}
  //       }
  //     }
  //   }
  
  counter.save_total(keeptotal);  
  counter.save();
  
  cout<<"finished"<<endl;
}


