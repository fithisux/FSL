/*  LSLaplaceManager.cc

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

#include "lslaplacemanager.h"
#include "utils/log.h"
#include "utils/tracer_plus.h"

using namespace Utilities;
using namespace MISCMATHS;
using namespace NEWIMAGE;

// priorstd=load('stddevs');clear tmp; d=ra('data'); dataprec = 1/(1)^2; priormean=0;  for i=1:length(priorstd), priorprec = 1/(priorstd(i))^2; tmp(i) = (dataprec*sum(squeeze(d(1,1,1,:)))+priorprec*priormean)/(1000*dataprec+priorprec);end;load means;plot(priorstd,means,priorstd,tmp);legend('mcmc','theory');

namespace Bint {

  float SumSquaresEvalFunction::evaluate(const ColumnVector& x) const 
    {
      Tracer_Plus tr("SumSquaresEvalFunction::evaluate");

      ntpts = data.Nrows();

      int nparams = model.getnparams();
      float energy = 1e16;

      if(!analmargprec)
	{
	  float precision;
	  if(updateprec)
	    precision = x(nparams);
	  else
	    precision = prec;
	  	  
	  if(precision > 0)
	    {
	      energy=(data-model.nonlinearfunc(x)).SumSquare()*precision/2.0 - ntpts/2.0*std::log(precision);
	      
	      for(int p=0;p<nparams;p++)
		{      
		  energy += model.getparam(p).getprior().calc_energy(x(p+1));
		}	      
	      
	      if(debuglevel==2)
		{ 
		  OUT(ntpts);
		  OUT((data-model.nonlinearfunc(x)).SumSquare()*precision/2.0 - ntpts/2.0*std::log(precision));

		  for(int p=0;p<nparams;p++)
		    { 
		      OUT(p);
		      OUT(x(p+1));
		      OUT(model.getparam(p).getprior().calc_energy(x(p+1)));
		    }
		  
		  OUT(energy);
		  OUT(precision);
		  OUT((data-model.nonlinearfunc(x)).SumSquare());
		  OUT(x);
		}
	    }
	}
      else
	{
	  energy = ntpts/2.0*std::log((data-model.nonlinearfunc(x)).SumSquare());
	  
	  for(int p=0;p<nparams;p++)
	    {      
	      energy += model.getparam(p).getprior().calc_energy(x(p+1));
	    }
	}

      return energy;
    }

  float SumSquaresgEvalFunction::evaluate(const ColumnVector& x) const 
    {
      Tracer_Plus tr("SumSquaresgEvalFunction::evaluate");
      
      ntpts = data.Nrows();

      int nparams = model.getnparams();
      float energy = 1e16;

      if(!analmargprec)
	{
	  float precision;
	  if(updateprec)
	    precision = x(nparams);
	  else
	    precision = prec;
	  	  
	  if(precision > 0)
	    {
	      energy=(data-model.nonlinearfunc(x)).SumSquare()*precision/2.0 - ntpts/2.0*std::log(precision);
	      
	      for(int p=0;p<nparams;p++)
		{      
		  energy += model.getparam(p).getprior().calc_energy(x(p+1));
		}	      
	      
	      if(debuglevel==2)
		{ 
		  OUT(ntpts);
		  OUT((data-model.nonlinearfunc(x)).SumSquare()*precision/2.0 - ntpts/2.0*std::log(precision));
		  for(int p=0;p<nparams;p++)
		    { 
		      OUT(p);
		      OUT(x(p+1));
		      OUT(model.getparam(p).getprior().calc_energy(x(p+1)));
		    }
		  
		  OUT(energy);
		  OUT(precision);
		  OUT((data-model.nonlinearfunc(x)).SumSquare());
		  OUT(x);
		}
	    }
	}
      else
	{
	  energy = ntpts/2.0*std::log((data-model.nonlinearfunc(x)).SumSquare());
	  
	  for(int p=0;p<nparams;p++)
	    {      
	      energy += model.getparam(p).getprior().calc_energy(x(p+1));
	    }
	}

//       OUT(x.t());
//       OUT(energy);

      return energy;
    }

  ReturnMatrix SumSquaresgEvalFunction::g_evaluate(const ColumnVector& x) const 
    {
      Tracer_Plus tr("SumSquaresgEvalFunction::g_evaluate");
      
      ntpts = data.Nrows();
      int nparams = model.getnparams();
      ColumnVector ret(x.Nrows());
      ret = 0;

      if(!analmargprec)
	{
	  float precision;
	  if(updateprec)
	    precision = x(nparams);
	  else
	    precision = prec;
	  	  
	  if(precision > 0)
	    {
	      // need to write this bit to make laplace work with gradient info when analmargprec is turned off
	      
	      for(int p=0;p<nparams;p++)
		{      
		  
		}
	      if(debuglevel==2)
		{ 
		  
		}
	    }
	}
      else
	{
	  float h = (data-model.nonlinearfunc(x)).SumSquare();
	  
	  Matrix grad = model.gradient(x);
	  
	  for(int p=1;p<=nparams;p++)
	    {    

	      ret(p) = -ntpts*SP(data-model.nonlinearfunc(x),grad.Row(p).AsColumn()).Sum()/h + model.getparam(p-1).getprior().calc_gradient(x(p));

// 	      OUT(float(SP(data-model.nonlinearfunc(x),grad.Row(p).AsColumn()).Sum()));
// 	      OUT(-ntpts*SP(data-model.nonlinearfunc(x),grad.Row(p).AsColumn()).Sum()/h);
// 	      OUT(model.getparam(p-1).getprior().calc_gradient(x(p)));
	    }
	}

//       OUT(ret.t());
      ret.Release();
      return ret;
    }

  void LSLaplaceManager::setup()
  {
    Tracer_Plus tr("LSLaplaceManager::setup");
    ntpts = data.Nrows();
    nvoxels = data.Ncols();
  }

  void LSLaplaceManager::run()
  {
    Tracer_Plus tr("LSLaplaceManager::run");
//      float mint = 1;
//      float maxt = 10;

//      ColumnVector stddevs(data.nvoxels());stddevs=0;
//      ColumnVector means(data.nvoxels());means=0;

    for(int vox=1;vox<=data.Ncols();vox++)    
      {
	cout << vox<< ",";
	cout.flush();

	if(debuglevel==2)
	  {	  
	    cout << endl;
	    cout << "----------------------------------" << endl;
	  }

//  	float stddev = float((maxt-mint)*vox)/data.nvoxels()+mint;
//  	OUT(stddev);

	voxelmanager->setdata(data.Column(vox));
	voxelmanager->setupparams(precin);
	nparams = voxelmanager->getnparams();
	int nvaryingparams = voxelmanager->getnvaryingparams();
	voxelmanager->run();
//  	stddevs(vox) = stddev;

	// use first voxel to get size of results storage
	if(vox==1)
	  {  
	    covs.ReSize(nvaryingparams*nvaryingparams,nvoxels); //should this be nparams to fix error below...
	    //covs.ReSize(nparams*nparams,nvoxels);

	    covs = 0; 
	    mns.ReSize(nparams,nvoxels);
	    mns = 0;

	    if(!analmargprec)
	      {
		prec.ReSize(nvoxels);
		prec = 0;
	      }
	  }
//  	ColumnVector tmp = mns.Column(vox);
//  	OUT(tmp.Nrows());
//  	tmp = voxelmanager->getparammeans();
//  	OUT(tmp.Nrows());

	mns.Column(vox) = voxelmanager->getparammeans();
	const SymmetricMatrix& symmat = voxelmanager->getparaminvcovs();
	
	if(!analmargprec)
	  prec(vox) = voxelmanager->geterrorprecisionmean();
	
//  	OUT(mns.Column(vox));
//  	OUT(prec(vox));
  	OUT(symmat);

	ColumnVector col = reshape(symmat.i(), Sqr(symmat.Nrows()), 1).AsColumn();

	OUT(symmat.i());
	//cerr << vox << "    here5          " << col.Nrows() << " " << col.Ncols() << " " << covs.Nrows() << " " << covs.Ncols() << endl;
	covs.Column(vox) = col; //ERROR this apppears to be broken (for certain inputs ) in the volumeseries lslaplace for noamp and still needs fixing
      }

    cout << endl;
  }

  void LSLaplaceManager::save()
  {
    Tracer_Plus tr("LSLaplaceManager::save");

    volume4D<float> output(mask);
    output.setmatrix(mns,mask[0]);

    for(int p=0;p<nparams;p++)
      {
	OUT(p);
	save_volume(output[p],LogSingleton::getInstance().appendDir(voxelmanager->getparamname(p)+string("_means")));	
      }
      mns.CleanUp();

      output.setmatrix(covs,mask[0]);
      save_volume4D(output,LogSingleton::getInstance().appendDir("covs"));
      covs.CleanUp();

      if(!analmargprec)
	{
	  output.setmatrix(prec.t(),mask[0]);
	  save_volume4D(output,LogSingleton::getInstance().appendDir("prec_means"));
	  prec.CleanUp(); 
	}
  }

  void LSLaplaceVoxelManager::setdata(const ColumnVector& pdata)
  {
    Tracer_Plus trace("LSLaplaceVoxelManager::setdata");

    data = pdata;    
    ntpts = data.Nrows();
  }

  void LSLaplaceVoxelManager::setupparams(float precin)
  {
    Tracer_Plus trace("LSLaplaceVoxelManager::setupparams");

    prec = precin;
    model.setparams();
    model.initialise(data);
    nparams = model.getnparams();

    nvaryingparams = 0;
    for(int p=0;p<nparams;p++)
      {
	if(model.getparam(p).getallowtovary()) nvaryingparams++;
      }

    if(!analmargprec)
      {
	// nparams now plus one to include precision:
	parammeans.ReSize(nparams+1);	     
      }
    else
      {
	parammeans.ReSize(nparams);
      }

    parammeans = 0;     

    // initialize param values to init values
    for(int p=1; p<=nparams;p++)
      {
	parammeans(p)  = model.getparam(p-1).getinitvalue();
      }

    if(!analmargprec)
      {       
	float errorprecision = 0.0;

	if(prec<0)
	  {
	    updateprec = true;
	    ColumnVector r = data-model.nonlinearfunc(parammeans);

	    if(updateprec)
	      {
		errorprecision = ntpts/SumSquare(r);
	      }	

	    float var = Sqr(errorprecision)*1000000;
	    float a = Sqr(errorprecision)/var;
	    float b = errorprecision/var;

	    GammaPrior tmpgamprior = GammaPrior(a,b);
	    model.add_param("prec",errorprecision,errorprecision/10.0,tmpgamprior);    
	    // set precision in parammeans
	    parammeans(nparams+1) = model.getparam(nparams).getinitvalue();
	    nparams =  model.getnparams();
	  }
	else 
	  {
	    updateprec = false;
	    errorprecision = prec;
	    // remove precision in parammeans
	    parammeans = parammeans.Rows(1,nparams);
	  }
      }

  }

 void LSLaplaceVoxelManager::run()
  {
    Tracer_Plus trace("LSLaplaceVoxelManager::run");
    
    if(debuglevel==2)
      { 
	OUT(parammeans.t());
	OUT(evalfunction->evaluate(parammeans));	
      }

    ColumnVector paramstovaryflags(parammeans.Nrows());
    for(int p=0;p<nparams;p++)
      {
	paramstovaryflags(p+1) = model.getparam(p).getallowtovary();
      }
    
    evalfunction->minimize(parammeans,paramstovaryflags);

    if(debuglevel==2)
      {
	OUT(parammeans.t());
	OUT(evalfunction->evaluate(parammeans));
      }
    
    bool finished = false;

    int power=-10;
    while(!finished && power<10)
      {	
	OUT(power);
	OUT(std::pow(double(10.0),double(power)));
	paraminvcovs = hessian(parammeans, *evalfunction, std::pow(double(10.0),double(power)), 4);
	finished = true;
	
	for(int p=0;p<nparams;p++)
	  if(paramstovaryflags(p+1))
	    if(paraminvcovs(p+1,p+1) == 0)      
	      {
		OUT(p);
		finished = false;
		power++;
		break;
	      }
      }

    if(debuglevel==2)
      {
	OUT(power);
	OUT(paraminvcovs);
      }

    // prune out non varying parameters:
    SymmetricMatrix paraminvcovstmp = paraminvcovs;
    paraminvcovstmp = 0;
    int vp = 0;
    for(int p=0;p<nparams;p++)
      if(paramstovaryflags(p+1))
	{
	  vp++;
	  paraminvcovstmp(vp,vp)=paraminvcovs(p+1,p+1);
	  // for(int q=0;q<nparams;q++)
// 	    if(paramstovaryflags(q+1))
// 	      paraminvcovstmp(vp,q+1)=paraminvcovs(p+1,q+1);
	}

    paraminvcovs = paraminvcovstmp.SymSubMatrix(1,vp);
   
    if(power > 9) 
      {
	cout << "Second derivative zero in hessian calculation" << endl;
	paraminvcovs << IdentityMatrix(nparams);
	//throw Exception("Second derivative zero in hessian calculation");
      }

  }

}
