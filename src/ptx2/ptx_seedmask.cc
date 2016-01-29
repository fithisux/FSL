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

#include "streamlines.h"
#include "ptx_seedmask.h"
#include <time.h>

using namespace std;
using namespace NEWIMAGE;
using namespace TRACT;
using namespace Utilities;
using namespace PARTICLE;
using namespace mesh;



void seedmask()
{ 
  probtrackxOptions& opts =probtrackxOptions::getInstance();
  
  // we need a reference volume for CSV
  // (in case seeds are a list of surfaces)
  volume<short int> refvol;
  if(opts.seedref.value()!="")
    read_volume(refvol,opts.seedref.value());
  else
    read_volume(refvol,opts.maskfile.value());
  

  cout<<"load seeds"<<endl;
  CSV seeds(refvol);
  seeds.set_convention(opts.meshspace.value());
  seeds.load_rois(opts.seedfile.value());

  cout<<"done."<<endl;
  if(seeds.nVols()==0 && opts.seedref.value()==""){
    cerr<<"Warning: need to set a reference volume when defining a surface-based seed"<<endl;
  }

  Streamliner  stline      (seeds);
  Counter      counter     (stline);
  counter.initialise();
  Seedmanager  seedmanager (counter);

  srand(opts.rseed.value()); // need to reinitialise random seed because of GIFTI!!
  
  int keeptotal=0;

  time_t _time;
  _time=time(NULL);
  // seed from volume-like ROIs
  if(seeds.nVols()>0){
    cout << "Volume seeds" << endl;
    vector<int> triangles; //to avoid connections between vertices of the same triangle...but not used for volumes
    triangles.push_back(-1);

    for(int roi=1;roi<=seeds.nVols();roi++){
      cout<<"volume "<<roi<<endl;

      for(int z=0;z<seeds.zsize();z++){
	if(opts.verbose.value()>=1)
	  cout <<"sl "<<z<<endl;
	for(int y=0;y<seeds.ysize();y++){
	  for(int x=0;x<seeds.xsize();x++){
	    if(seeds.isInRoi(x,y,z,roi)){
	      counter.updateSeedLocation(seeds.get_volloc(roi-1,x,y,z),-1,triangles);
	      if(opts.verbose.value()>=1){
		cout <<"run"<<endl;
		cout <<x<<" "<<y<<" "<<z<<endl;
	      }
	      keeptotal += seedmanager.run((float)x,(float)y,(float)z,
					   false,-1,opts.sampvox.value());
	    }
	  }
	}
      }
      
    }

  }

  // seed from surface-like ROIs
  if(seeds.nSurfs()>0){
    cout << "Surface seeds" << endl;
    ColumnVector pos;//,dir;
    for(int i=0;i<seeds.nSurfs();i++){
      cout<<"surface "<<i<<endl;

      // inform user if whole surface is used or not
      if( seeds.nActVertices(i) != seeds.nVertices(i) ){
	cout << "  Using a subset of the vertices labelled active (i.e. non zero value)" << endl;
	cout << "   set all values to 0 or non-zero to use entire surface" << endl;
      }
      for(int p=0;p<seeds.get_mesh(i).nvertices();p++){
	// check if active point	
	if(seeds.get_mesh(i).get_pvalue(p)==0.0)
	  continue;

        //to avoid connections between vertices of the same triangle
	CsvMpoint vertex=seeds.get_mesh(i).get_point(p);
	vector<int> triangles;
	for(int t=0;t<vertex.ntriangles();t++){
	  triangles.push_back(vertex.get_trID(t));
	}
	
	counter.updateSeedLocation(seeds.get_surfloc(i,p),i,triangles);
	pos=seeds.get_vertex_as_vox(i,p);

	ColumnVector dir(3);
	dir=seeds.get_normal_as_vox(i,p);

	//if(opts.meshspace.value()=="caret")
	// dir*=-1; // normals in caret point away from the brain

	if(opts.verbose.value()>=1){
	  cout <<"run"<<endl;
	  cout <<pos(1)<<" "<<pos(2)<<" "<<pos(3)<<endl;
	}
	keeptotal += seedmanager.run(pos(1),pos(2),pos(3),
				       false,-1,opts.sampvox.value());
	

      }
    }
  }

  cout<<endl<<"time spent tracking: "<<(time(NULL)-_time)<<" seconds"<<endl<<endl;

  // save results
  cout << "save results" << endl;
  counter.save_total(keeptotal);  
  counter.save();

  cout<<"finished"<<endl;
}


