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

#include "ptx_twomasks.h"
#include "streamlines.h"

using namespace std;
using namespace NEWIMAGE;
using namespace TRACT;
using namespace Utilities;
using namespace PARTICLE;
using namespace mesh;



void twomasks()
{ 
  probtrackxOptions& opts =probtrackxOptions::getInstance();

  ////////////////////////////////
  //  Log& logger = LogSingleton::getInstance();
  volume<float> seeds,seeds2;
  read_volume(seeds,opts.seedfile.value());
  read_volume(seeds2,opts.mask2.value());

  if(opts.s2tout.value()){
    cerr << "Seed_to_target not available in multiple seed tractography" << endl;
    exit(0);
  }
  // correct for non-1 values
  seeds=NEWIMAGE::abs(seeds);
  seeds.binarise(0,seeds.max()+1,exclusive);
  seeds2=NEWIMAGE::abs(seeds2);
  seeds2.binarise(0,seeds2.max()+1,exclusive);

  int numseeds=0;
  for(int z=0;z<seeds.zsize();z++)
    for(int y=0;y<seeds.ysize();y++)
      for(int x=0;x<seeds.xsize();x++)
	if(seeds(x,y,z)!=0)numseeds++;	


  Streamliner stline(seeds);
  Counter counter(seeds,stline,numseeds);
  counter.initialise();
  Seedmanager seedmanager(counter);
  
  vector<int> keeptotal(2);

  stline.add_waymask(seeds2);
  keeptotal[0] = 0;
  cout << "tracking from first seed" << endl;
  for(int z=0;z<seeds.zsize();z++){
    for(int y=0;y<seeds.ysize();y++){
      for(int x=0;x<seeds.xsize();x++){
	if(seeds(x,y,z)>0){
	  keeptotal[0] += seedmanager.run(x,y,z,false,-1); 
	}
      }
    }
  }
  stline.pop_waymasks();
  
  stline.add_waymask(seeds);
  keeptotal[1] = 0;
  cout << "tracking from second seed" << endl;
  for(int z=0;z<seeds2.zsize();z++){
    for(int y=0;y<seeds2.ysize();y++){
      for(int x=0;x<seeds2.xsize();x++){
	if(seeds2(x,y,z)>0){	  
	  keeptotal[1] += seedmanager.run(x,y,z,false,-1); 
	  //keeptotal[1] += seedmanager.run(x,y,z,false,seeds2(x,y,z)-1); 
	}
      }
    }
  }
  

  counter.save_total(keeptotal);
  counter.save();
  cout<<"finished"<<endl;
}


