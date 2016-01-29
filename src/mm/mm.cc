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
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#define WANT_STREAM
#define WANT_MATH
#include "newmatap.h"
#include "newmatio.h"
#include <string>
#include <math.h>
#include "utils/log.h"
#include "mmoptions.h"
#include "utils/tracer_plus.h"
#include "newimage/newimageall.h"
#include "mixture_model.h"

using namespace Utilities;
using namespace NEWMAT;
using namespace Mm;
using namespace NEWIMAGE;

ReturnMatrix volume2col(const volume<float>& spatial_data, const volume<int>& mask)
{
    int num_superthreshold = 0;
    for(int x = 0; x < mask.xsize(); x++)
      for(int y = 0; y < mask.ysize(); y++)
	for(int z = 0; z < mask.zsize(); z++)
	  if(mask(x,y,z))
	    {
	      num_superthreshold++;	      
	    }
    
    ColumnVector Y(num_superthreshold);
    Y = 0;
    int index=1;
    for(int z = 0; z < mask.zsize(); z++)    
      for(int y = 0; y < mask.ysize(); y++)	
	for(int x = 0; x < mask.xsize(); x++)
	  if(mask(x,y,z))
	    {
	      Y(index) = spatial_data(x,y,z);

	      index++;
	    }

    Y.Release();
    return Y;
}
int main(int argc, char *argv[])
{
  try{
    double tmp=10;
    tmp=10;    OUT(exp(tmp));
    tmp=100;   OUT(exp(tmp));
    tmp=200;   OUT(exp(tmp));
    tmp=300;   OUT(exp(tmp));
    tmp=400;   OUT(exp(tmp));
    tmp=500;   OUT(exp(tmp));
    tmp=600;   OUT(exp(tmp));
    tmp=700;   OUT(exp(tmp));

    // Setup logging:
    Log& logger = LogSingleton::getInstance();
    
    // parse command line - will output arguments to logfile
    MmOptions& opts = MmOptions::getInstance();
    opts.parse_command_line(argc, argv, logger);

    //Tracer_Plus::setinstantstackon();

    if(opts.debuglevel.value()==1)
      Tracer_Plus::setrunningstackon();

    if(opts.timingon.value())
      Tracer_Plus::settimingon();
 
    ////////////////////////////////

    volume<float> epi_example_data;
    volume<int> mask;
    volume<float> spatial_data;

    int epibt = 0;
    cout << "epibt =" << epibt << endl;

    cout << "spatialdatafile =" << opts.spatialdatafile.value() << endl;    
    read_volume(spatial_data, opts.spatialdatafile.value());

    bool overlay = MmOptions::getInstance().epiexampledatafile.value()!="";

    if(MmOptions::getInstance().epiexampledatafile.value()!="")
      {   
	cout << "epiexampledatafile =" << opts.epiexampledatafile.value() << endl;
	read_volume(epi_example_data, opts.epiexampledatafile.value());
      }
    else
      {
	epi_example_data.reinitialize(spatial_data.xsize(),spatial_data.ysize(),spatial_data.zsize());
	epi_example_data = 0;
      }
    
    cout << "maskfile =" << opts.maskfile.value() << endl;
    read_volume(mask, opts.maskfile.value());

    vector<Distribution*> dists;
    vector<volume<float> > w_means;

    ColumnVector Y = volume2col(spatial_data, mask);       
    float minmode=0.5;

    OUT(minmode);
    bool zfstatmode = opts.zfstatmode.value();
    //bool zfstatmode = false;

    if(!zfstatmode)
      {  
	// standard z SPM

  	GaussianDistribution nonactdist;
  	GaussianDistribution gaussian_actdist;      
  	GaussianDistribution gaussian_deactdist;
	
	GammaDistribution gamma_actdist;  
	FlippedGammaDistribution fgamma_deactdist;
	
	//dists.push_back(&fgamma_nonactdist);
	dists.push_back(&nonactdist);
 
	//dists.push_back(&gaussian_actdist);
 	dists.push_back(&gamma_actdist);

	//dists.push_back(&gaussian_deactdist);
  	dists.push_back(&fgamma_deactdist);

	//	ggmfit(Y.t(), dists, MmOptions::getInstance().nonspatial.value());
	ggmfit(Y.t(), dists, true);
	
	// set minmodes based on non-spatial fit:
	float nonactmean = dists[0]->getmean();
	float setminmode = Max(1.5*nonactmean+std::sqrt(dists[0]->getvar()),minmode);
	gamma_actdist.setminmode(setminmode);
	OUT(setminmode);

	setminmode = Min(nonactmean-1.5*std::sqrt(dists[0]->getvar()),-minmode);
	fgamma_deactdist.setminmode(setminmode);
	OUT(setminmode);

	Mixture_Model mm(spatial_data, mask, epi_example_data, epibt, dists, w_means, Y, MmOptions::getInstance());
	
	mm.setup();
	mm.run();
	mm.save();
	
	plot_ggm(w_means,dists,mask,Y);
	make_ggmreport(w_means,dists,mask,spatial_data,zfstatmode,overlay,epi_example_data,opts.threshold.value(),opts.nonspatial.value(), opts.updatetheta.value(),opts.spatialdatafile.value());
      }

    else
      {	
	// z SPM from an f-statistic - just two classes needed (no deactivation)		   

	GaussianDistribution gauss_nonactdist;	    
	GammaDistribution gamma_actdist;
	GaussianDistribution gauss_actdist;

	dists.push_back(&gauss_nonactdist);

	dists.push_back(&gamma_actdist);
	//dists.push_back(&gauss_actdist);

	//	ggmfit(Y.t(),dists, MmOptions::getInstance().nonspatial.value());
	ggmfit(Y.t(), dists, true);

	// set minmodes based on non-spatial fit:
	float nonactmean = dists[0]->getmean();
	gamma_actdist.setminmode(Max(1.5*nonactmean+std::sqrt(dists[0]->getvar()),minmode));	

	OUT(minmode);
	OUT(std::sqrt(dists[0]->getvar()));
	OUT(nonactmean);
	OUT(Max(nonactmean+std::sqrt(dists[0]->getvar())*1.5,minmode));

	Mixture_Model mm(spatial_data, mask, epi_example_data, epibt, dists, w_means, Y, MmOptions::getInstance());
	    
	mm.setup();	
	mm.run();
	mm.save();    
	    
	plot_ggm(w_means,dists,mask,Y);
	make_ggmreport(w_means,dists,mask,spatial_data,zfstatmode,overlay,epi_example_data,opts.threshold.value(),opts.nonspatial.value(), opts.updatetheta.value(), opts.spatialdatafile.value());
	 
      }
    
    ////////////////////////////

    if(opts.timingon.value())
      Tracer_Plus::dump_times(logger.getDir());
    
    cout << "Log directory was: " << logger.getDir() << endl;

  }
  catch(Exception& e) 
    {
      cerr << endl << e.what() << endl;
      return 1;
    }
  catch(X_OptionError& e) 
    {
      cerr << endl << e.what() << endl;
      return 1;
    }

  return 0;
}












