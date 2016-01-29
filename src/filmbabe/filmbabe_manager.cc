/*  filmbabe_manager.cc

    Mark Woolrich, FMRIB Image Analysis Group

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

#include "filmbabe_manager.h"
#include "utils/log.h"
#include "miscmaths/miscmaths.h"
#include "newimage/newimageall.h"
#include "utils/tracer_plus.h"

using namespace Utilities;
using namespace MISCMATHS;
using namespace NEWIMAGE;

namespace Filmbabe {

  Filmbabe_Manager::~Filmbabe_Manager() 
  { 
    delete filmbabe_vb_flobs;
  }

  void Filmbabe_Manager::setup()
  {
    Tracer_Plus trace("Filmbabe_Manager::setup");
    
    if(FilmbabeOptions::getInstance().debuglevel.value()==2)
      {
	cout << "******************************************" << endl
	     << "SETUP" << endl << "******************************************"
	     << endl;
      }

    cout << "datafile =" << opts.datafile.value() << endl;    
    read_volume4D(data, opts.datafile.value());
    
    cout << "maskfile =" << opts.maskfile.value() << endl;
    copybasicproperties(data[0],mask);
    read_volume(mask, opts.maskfile.value());

    cout << "designfile =" << opts.designfile.value() << endl;
    designmatrix = read_vest(opts.designfile.value()).t();

    cout << "flobsregressorsfile =" << opts.flobsregressorsfile.value() << endl;
    flobsregressors = read_ascii_matrix(opts.flobsregressorsfile.value());

    OUT(flobsregressors.t());
    OUT(designmatrix.Nrows());
    OUT(designmatrix.Ncols());

    localweights.reinitialize(data.xsize(),data.ysize(),data.zsize(),6);
    localweights = 0;
    int num_superthreshold = 0;
  
    for(int x = 0; x < data.xsize(); x++)
      for(int y = 0; y < data.ysize(); y++)
	for(int z = 0; z < data.zsize(); z++)
	  if(mask(x,y,z))
	    {
	      num_superthreshold++;
	      
	      int xi,yi,zi;
	      for(unsigned int i = 0; i < connected_offsets.size(); i++) 
		{
		  xi = x+connected_offsets[i].x;
		  yi = y+connected_offsets[i].y;
		  zi = z+connected_offsets[i].z;
		  
		  if(mask(xi,yi,zi))
		    {
		      localweights(x,y,z,connected_offsets[i].ind) = 1;
		    }
		}
	    }

    OUT(num_superthreshold);       
   
    filmbabe_vb_flobs = new Filmbabe_Vb_Flobs(data, mask, designmatrix, flobsregressors, localweights, connected_offsets, num_superthreshold);   
    filmbabe_vb_flobs->setup();
     
  } 
  
  void Filmbabe_Manager::run()
  {
    Tracer_Plus trace("Filmbabe_Manager::run");
    
    if(FilmbabeOptions::getInstance().debuglevel.value()==2)
      {
	cout << "******************************************" << endl
	     << "RUN" << endl << "******************************************" 
	     << endl;
      }

    filmbabe_vb_flobs->run();
   
    cout << endl << "Finished" << endl;
  }

  void Filmbabe_Manager::save()
  {
    Tracer_Plus trace("Filmbabe_Manager::save");

    filmbabe_vb_flobs->save();

    // save other stuff
    save_volume(mask, LogSingleton::getInstance().appendDir("mask"));

    save_volume4D(data, LogSingleton::getInstance().appendDir("data"));

  }

}
