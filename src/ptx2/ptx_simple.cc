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
#include "ptx_simple.h"

using namespace std;
using namespace NEWIMAGE;
using namespace TRACT;
using namespace Utilities;
using namespace PARTICLE;
using namespace mesh;


void track(){
  probtrackxOptions& opts =probtrackxOptions::getInstance();
  
  ////////////////////////////
  Log& logger = LogSingleton::getInstance();
  if(opts.verbose.value()>1){
    logger.makeDir("particles","particle0",true,false);
  }
  
  volume<short int> seedref;
  if(opts.seedref.value()!=""){
    read_volume(seedref,opts.seedref.value());
  }
  else{
    read_volume(seedref,opts.maskfile.value());
  }
  CSV csv(seedref);//we need this for seed2target

  Matrix Seeds = read_ascii_matrix(opts.seedfile.value());
  if(Seeds.Ncols()!=3 && Seeds.Nrows()==3)
    Seeds=Seeds.t();

  Streamliner   stline(csv);
  Counter       counter(stline);
  counter.forceNumSeeds(Seeds.Nrows()); // csv useless in this case
  counter.initialise();
  Seedmanager   seedmanager(counter);

  srand(opts.rseed.value()); // need to reinitialise random seed because of GIFTI!!
  
  // convert coordinates from nifti (external) to newimage (internal)
  //   conventions - Note: for radiological files this should do nothing
  Matrix newSeeds(Seeds.Nrows(),3);
  for (int n=1; n<=Seeds.Nrows(); n++) {
    ColumnVector v(4);
    v << Seeds(n,1) << Seeds(n,2) << Seeds(n,3) << 1.0;
    v = seedref.niftivox2newimagevox_mat() * v;
    newSeeds.Row(n) << v(1) << v(2) << v(3);
  }

  time_t _time;
  _time=time(NULL);

  int keeptot=0;
  vector<int> triangles;  //to avoid connections between vertices of the same triangle...but not used for volumes
  triangles.push_back(-1);
  for(int SN=1; SN<=newSeeds.Nrows();SN++){
    counter.updateSeedLocation(SN,-1,triangles);
    float xst=newSeeds(SN,1);
    float yst=newSeeds(SN,2);
    float zst=newSeeds(SN,3);
    keeptot += seedmanager.run(xst,yst,zst,false,-1,opts.sampvox.value());
    string add="_"+num2str(Seeds(SN,1))+"_"+num2str(Seeds(SN,2))+"_"+num2str(Seeds(SN,3));
    
    if(opts.simpleout.value()||opts.omeanpathlength.value())
      counter.save_pathdist(add);

    counter.reset_prob();
  } //Close Seed number Loop
  
  cout<<endl<<"time spent tracking: "<<(time(NULL)-_time)<<" seconds"<<endl<<endl;


  cout << "save results" << endl;
  counter.save();

  cout<<"finished"<<endl;
}
