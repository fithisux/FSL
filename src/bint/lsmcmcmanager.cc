/*  LSMCMCManager.cc

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

#include "lsmcmcmanager.h"
#include "utils/log.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"

#include "utils/tracer_plus.h"
#include "model.h"

using namespace Utilities;
using namespace MISCMATHS;


namespace Bint {

  void LSMCMCManager::setup()
  {
    Tracer_Plus tr("LSMCMCManager::setup"); 

    ntpts = data.Nrows();
    nvoxels = data.Ncols();

    // need to get nparams by dummy call to model:
    model.setparams();
    nparams = model.getnparams();    
	    
    // initialise samples
    samples.resize(nparams);
    cout << "nparams=" << nparams << endl;
    cout << "nsamples=" << nsamples << endl;

    // initialise some more in case nvoxels is zero
    for(int p=0; p<nparams; p++) {      
      samples[p].ReSize(nsamples,nvoxels);
      samples[p] = 0;

      paramnames.push_back(model.getparam(p).getname());	  
    }       
    if(!analmargprec)
      {
	precsamples.ReSize(nsamples,nvoxels);
	precsamples = 0;
      }
  }

  void LSMCMCManager::run()
  {
    Tracer_Plus tr("LSMCMCManager::run");

    // now run each voxel
    for(int vox=1;vox<=data.Ncols();vox++)
      {
	cout << vox<< ",";
	cout.flush();

	if(debuglevel==2)
	  {	  
	    cout << endl;
	    cout << "----------------------------------" << endl;
	  }

	voxelmanager.setdata(data.Column(vox));	
	voxelmanager.setupparams(precin);	
	voxelmanager.run();	

	// Now store samples:
	for(int p=0; p<nparams; p++) {
	  if(voxelmanager.getmcmcparams()[p]->getallowtovary())
	    {
	      const vector<float>& samps = voxelmanager.getsamples(p);
	      samples[p].Column(vox)=vector2ColumnVector(samps);
	    }

	  if(!analmargprec)
	    {
	      const vector<float>& precsamps = voxelmanager.getprecsamples();	
	      precsamples.Column(vox)=vector2ColumnVector(precsamps);
	    }
	}      
      }

    cout << endl;
  }

  void LSMCMCManager::save()
  {
    Tracer_Plus tr("LSMCMCManager::save");
    
    cout << "Saving results...";

    for(int p=0; p<nparams; p++) {
      if(model.getparam(p).getallowtovary() && model.getparam(p).getsave())
	{
          volume4D<float> output(mask);
	  output.setmatrix(samples[p],mask[0]);
          save_volume4D(output,LogSingleton::getInstance().appendDir(paramnames[p]+string("_samples")));
	  samples[p].CleanUp();
	}
    }

    if(!analmargprec)
      {
        volume4D<float> output(mask);
	output.setmatrix(precsamples,mask[0]);
        save_volume4D(output,LogSingleton::getInstance().appendDir("prec_samples"));
	precsamples.CleanUp();
      }

    cout << " finished" << endl;
  }

  void McmcParameter::jump()
  {
    Tracer_Plus trace("McmcParameter::jump");

    if(debuglevel==2)
      {
	cout << param.getname() << " jump" << endl;
	OUT(param.getallowtovary());
	OUT(val);
	OUT(normrnd().AsScalar());
      }

    // store old values
    float old = val;

    // propose new values
    val += normrnd().AsScalar()*proposal_std;      

    // calculate acceptance threshold   
    float tmp = unifrnd().AsScalar();    
    float tmpold = old_energy();
    float tmpnew = new_energy();

    bool accept=false;

    if(tmpnew!=float(MAX_EN))
      {	
	accept = ((tmpold - tmpnew) > std::log(tmp));    
      }

    if(debuglevel==2)
      {
	float numer=(tmpold - tmpnew);
	float denom=std::log(tmp);
	OUT(numer);
	OUT(denom);
	OUT(tmp);
	cout << "proposal_std=" << proposal_std << endl;
	cout << "old=" << old << endl;
	cout << "val=" << val << endl;
	cout << "old_energy()=" << tmpold << endl;
	cout << "new_energy()=" << tmpnew << endl;
	cout << "accept=" << accept << endl;
      }

    if(accept)
      {	
	naccepted++;
      }
    else
      {
	// restore old values	    
	val = old;
	restore_energy();
	nrejected++;
      }

    if(jumpcount>updateproposalevery)
      {
	update_proposal_std();
	jumpcount = 0;
      }
    else
      {	
	jumpcount++;
      }
  }
  
  float LSMCMCPrecParameter::calc_extra() 
  { 
    Tracer_Plus trace("LSMCMCPrecParameter::calc_extra");

    extra_old_energy = extra_energy; 
    
    float minprec = 0;
    
    if(val <= minprec)
      {
	extra_energy = 1e16;
	impropercount++;
	
	if(impropercount==int(lsmcmc.getnsamples()/4.0))
	  cout << "too many low precisions for LSMCMCPrecParameter" << endl;
      }
    else
      {	       
	extra_energy = -(N/2.0)*std::log(val)+param.getprior().calc_energy(val);	
      }
    
    if(debuglevel==2)
      {
	cout << "extra_old_energy=" << extra_old_energy << endl;
	cout << "extra_energy=" << extra_energy << endl;
      }
      
    return extra_energy; 
  }  
  
  void LSMCMCVoxelManager::calcsumsquares()
    {
      Tracer_Plus trace("LSMCMCVoxelManager::calcsumsquares");

      sumsquares_old = sumsquares;

      ColumnVector x(nparams);
      x = 0;

      for(int p=0;p<nparams;p++)
	{
	  x(p+1) = mcmcparams[p]->value();
	}

      ColumnVector tmp = model.nonlinearfunc(x);
      sumsquares=(data-tmp).SumSquare();
    }

  void LSMCMCVoxelManager::setupparams(float prec)
  {
    Tracer_Plus trace("LSMCMCVoxelManager::setupparams");
    
    model.setparams();
    model.initialise(data);
    nparams = model.getnparams();

    // set mcmcparams
    mcmcparams.clear();
    for(int p=0; p<nparams; p++)
      {
	mcmcparams.push_back(new LSMCMCParameter(model.getparam(p),nsamples,updateproposalevery,acceptancerate,*this));
	mcmcparams[p]->setup();
      }

    sumsquares = 0;
    calcsumsquares();
    
    //    OUT(analmargprec);

    if(!analmargprec)
      {
	// setup prec
	float mean = 0.0;
	if(prec<=0)
	  {
	    updateprec = true;
	    mean = ntpts/sumsquares;
	  }
	else 
	  {
	    updateprec = false;
	    mean = prec;
	  }      

	float var = Sqr(mean)*1000000;
	float a = Sqr(mean)/var;
	float b = mean/var;
	
	precparamprior = new GammaPrior(a,b);
	precparam = new Parameter("prec", mean, mean/10.0, *precparamprior);
	precmcmcparam = new LSMCMCPrecParameter(*precparam,nsamples,updateproposalevery,acceptancerate,*this);
	precmcmcparam->setup();
      }

    likelihood = 0;

    calclikelihood();
  }

  void LSMCMCVoxelManager::setdata(const ColumnVector& pdata)
  {
    Tracer_Plus trace("LSMCMCVoxelManager::setdata");

    data = pdata;
    ntpts = data.Nrows();
  }

  void LSMCMCVoxelManager::run()
  {
    Tracer_Plus trace("LSMCMCVoxelManager::run");
    
    int samples = 0;
    int jumps = 0;
    int subsamplejumps = 0;

      ColumnVector x(nparams);
      x = 0;

      for(int p=0;p<nparams;p++)
	{
	  x(p+1) = mcmcparams[p]->value();
	}
      
      ColumnVector retstart = model.nonlinearfunc(x);

    while(true)
      {
	jumps++;
	subsamplejumps++;
	    
  	jump();

	if(subsamplejumps >= sampleevery)
	  {
	    subsamplejumps = 0;
	    
	    // sample components after burnin
	    if(jumps > burnin)
	      {	   		
		sample();
		samples++;
		
		if(samples>=nsamples)
		  break;
	      }
	  }
      }    

      x = 0;

      for(int p=0;p<nparams;p++)
	{
	  x(p+1) = mcmcparams[p]->value();
	}
      
      ColumnVector retend = model.nonlinearfunc(x);
  }

  void LSMCMCVoxelManager::sample() 
  { 
    Tracer_Plus trace("LSMCMCVoxelManager::sample");

    for(int p=0; p<nparams; p++) {
      if(mcmcparams[p]->getallowtovary())
	mcmcparams[p]->sample();
    }

    if(!analmargprec)
      precmcmcparam->sample();
  }
  
  void LSMCMCVoxelManager::jump() 
  { 
    Tracer_Plus trace("LSMCMCVoxelManager::jump");

    if(debuglevel==2)
      cout << "LSMCMCVoxelManager::jump-----------" << endl;

    for(int p=0; p<nparams; p++) {
      if(mcmcparams[p]->getallowtovary())
	mcmcparams[p]->jump();
    }    

    if(!analmargprec)
      {
	if(updateprec)
	  {
	    if(debuglevel==2)
	      cout << "prec jump" << endl;
	    precmcmcparam->jump();
	  }
      }

    if(debuglevel==2)
      cout << "-----------------------------------" << endl;


  }


}
