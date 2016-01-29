/*  mcmc.cc

    Mark Woolrich, Tim Behrens  FMRIB Image Analysis Group

    Copyright (C) 2002 University of Oxford  */

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

#include "mcmc.h"
#include "utils/log.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include "newimage/newimageall.h"
#include "utils/tracer_plus.h"
#include "gam.h"
#include <set>

using namespace Utilities;
using namespace MISCMATHS;
using namespace NEWIMAGE;

namespace Gs {


  void Mcmc::setup()
  {
    Tracer_Plus trace("Mcmc::setup");
        
    gamma_latest.ReSize(nevs);
    gamma_samples.ReSize(nevs,nsamples);
    beta_latest.ReSize(ngs);
    m_latest.ReSize(ntpts);
    gamma_latest = gamma_mean;
    gamma_samples = 0;
    beta_latest = SP(beta_b,beta_c);
    m_latest = 0;

    c_samples.ReSize(nsamples);

    for(int t = 1; t <= ntpts; t++)
      {
	m_latest(t) = normal.Next()*sqrt(m_var(t))+m_mean(t);
      }

  } 

  void Mcmc::jump(bool relax)
  {
    Tracer_Plus trace("Mcmc::jump");
    
    if(!switched)
      {
	m_jump(relax);
	beta_jump(relax);
	gamma_jump(relax);
      }
    else
      {
	m_jump_switched(relax);
	beta_jump_switched(relax);
	gamma_jump_switched(relax);
      }

  } 

  void Mcmc::m_jump(bool relax)
  {
    Tracer_Plus trace("Mcmc::m_jump");
    
    for(int t = 1; t <= ntpts; t++)
      {
	float sumovere = 0;
	for(int e = 1; e <= nevs; e++)
	  {
	    sumovere += design.getdm()(t,e)*gamma_latest(e)*beta_latest(design.getgroup(t));
	  }
	float preccope = 1.0/varcopedata(t);

	float invstd =  sqrt(beta_latest(design.getgroup(t))+preccope);
	float mn = (sumovere + copedata(t)*preccope)/(beta_latest(design.getgroup(t))+preccope);

	if(relax)
	  {
	    // do ordered overrelaxation
	    multiset<float> ms;
	    ms.insert(normal.Next()/invstd + mn);
	    
	    for(int k = 2; k <= 20; k++)
	      {
		ms.insert(normal.Next()/invstd + mn);	    
	      }
	    
	    int d = distance(ms.insert(m_latest(t)),ms.end())-1;
	    multiset<float>::iterator iter = ms.begin();
	    advance(iter,d);
	
	    m_latest(t) = *iter;
	  }
	else
	  {
	    m_latest(t) = normal.Next()/invstd + mn;
	  }
      }
  }
  
  void Mcmc::m_jump_switched(bool relax)
  {
    Tracer_Plus trace("Mcmc::m_jump_switched");
    
    for(int t = 1; t <= ntpts; t++)
      {
	float sumovere = 0; 
	for(int e = 1; e <= nevs; e++)
	  {
	    sumovere += design.getdm()(t,e)*gamma_latest(e)/varcopedata(t);
	  }
	float preccope = 1.0/varcopedata(t);	

	float invstd = sqrt(preccope+beta_latest(design.getgroup(t)));
	float mn = (sumovere + copedata(t)*beta_latest(design.getgroup(t)))/(preccope+beta_latest(design.getgroup(t)));

	if(relax)
	  {
	    // do ordered overrelaxation
	    multiset<float> ms;
	    ms.insert(normal.Next()/invstd + mn);
	    
	    for(int k = 2; k <= 20; k++)
	      {
		ms.insert(normal.Next()/invstd + mn);	    
	      }
	
	    int d = distance(ms.insert(m_latest(t)),ms.end())-1;
	    multiset<float>::iterator iter = ms.begin();
	    advance(iter,d);
	
	    m_latest(t) = *iter;
	  }
	else
	  {
	    m_latest(t) = normal.Next()/invstd + mn;
	  }

      }

  }

  void Mcmc::beta_jump(bool relax)
  {
    Tracer_Plus trace("Mcmc::beta_jump");

    for(int g = 1; g <= ngs; g++)
      {
	float sumovert = 0; 
	for(int t = 1; t <= ntpts; t++)
	  {
	    if(design.getgroup(t)==g)
	      {
		float sumovere = 0; 
		for(int e = 1; e <= nevs; e++)
		  {
		    sumovere += design.getdm()(t,e)*gamma_latest(e);
		  }
		sumovert += Sqr(m_latest(t)-sumovere);
	      }
	  }

	c_latest = 1e-6+0.5*sumovert;
	
	if(relax)
	  {
	    // do ordered overrelaxation
	    multiset<float> betas;
	    betas.insert(gam.rnd(1e-6+design.getntptsingroup(g)/2.0, 1e-6+0.5*sumovert));
	    
	    for(int k = 2; k <= 20; k++)
	      {
		betas.insert(gam.rnd());	    
	      }
	    
	    int d = distance(betas.insert(beta_latest(g)),betas.end())-1;
	    multiset<float>::iterator iter = betas.begin();
	    advance(iter,d);
	    
	    beta_latest(g) = *iter;
	  }
	else
	  {
	    beta_latest(g) = gam.rnd(1e-6+design.getntptsingroup(g)/2.0, 1e-6+0.5*sumovert);
	  }
      }
  }

  void Mcmc::beta_jump_switched(bool relax)
  {
    Tracer_Plus trace("Mcmc::beta_jump_switched");

    for(int g = 1; g <= ngs; g++)
      {
	float sumovert = 0; 
	for(int t = 1; t <= ntpts; t++)
	  {
	    if(design.getgroup(t)==g)
	      {
		sumovert += Sqr(copedata(t)-m_latest(t));
	      }
	  }

	c_latest = 1e-6+0.5*sumovert;

	if(relax)
	  {
	    // do ordered overrelaxation
	    multiset<float> betas;
	    betas.insert(gam.rnd(1e-6+design.getntptsingroup(g)/2.0, 1e-6+0.5*sumovert));
	
	    for(int k = 2; k <= 20; k++)
	      {
		betas.insert(gam.rnd());	    
	      }
	
	    int d = distance(betas.insert(beta_latest(g)),betas.end())-1;
	    multiset<float>::iterator iter = betas.begin();
	    advance(iter,d);

	    beta_latest(g) = *iter;

	  }
	else
	  {
	    beta_latest(g) = gam.rnd(1e-6+design.getntptsingroup(g)/2.0, 1e-6+0.5*sumovert);
	  }
      }

  }

  void Mcmc::gamma_jump(bool relax)
  {
    Tracer_Plus trace("Mcmc::gamma_jump");

    ColumnVector gammean(nevs);
    SymmetricMatrix gamcovar(nevs);

    for(int e = 1; e <= nevs; e++)
      {
	float sumovert1 = 0;
	float sumovert2 = 0;
	for(int t = 1; t <= ntpts; t++)
	  {
	    sumovert1 += design.getdm()(t,e)*beta_latest(design.getgroup(t))*m_latest(t);
	    sumovert2 += Sqr(design.getdm()(t,e))*beta_latest(design.getgroup(t));
	  }

	gammean(e) = sumovert1/(sumovert2);

	for(int e2 = e; e2 <= nevs; e2++)
	  {
	    float sumovert1 = 0;
	    for(int t = 1; t <= ntpts; t++)
		sumovert1 +=  design.getdm()(t,e)*design.getdm()(t,e2)*beta_latest(design.getgroup(t));
	      
	    gamcovar(e,e2) =  1.0/(sumovert1);
	  }
      }

    gamma_latest = mvnormrandm(gammean,gamcovar);
 
  }

  void Mcmc::gamma_jump_switched(bool relax)
  {
    Tracer_Plus trace("Mcmc::gamma_jump_switched");

    ColumnVector gammean(nevs);
    SymmetricMatrix gamcovar(nevs);

    for(int e = 1; e <= nevs; e++)
      {
	float sumovert1 = 0;
	float sumovert2 = 0;
	for(int t = 1; t <= ntpts; t++)
	  {
	    float preccope = 1.0/varcopedata(t);
	    sumovert1 += design.getdm()(t,e)*preccope*m_latest(t);
	    sumovert2 += Sqr(design.getdm()(t,e))*preccope;
	  }

	gammean(e) = sumovert1/(sumovert2);

	for(int e2 = e; e2 <= nevs; e2++)
	  {
	    float sumovert1 = 0;
	    for(int t = 1; t <= ntpts; t++)
	      {
		float preccope = 1.0/varcopedata(t);
		sumovert1 +=  design.getdm()(t,e)*design.getdm()(t,e2)*preccope;
	      }
	    gamcovar(e,e2) =  1.0/(sumovert1);
	  }
      }

    gamma_latest = mvnormrandm(gammean,gamcovar);
  }

  void Mcmc::sample(int samp)
  {
    Tracer_Plus trace("Mcmc::sample");

    for(int e = 1; e <= nevs; e++)
      {
	gamma_samples(e,samp) = gamma_latest(e);	
      }
    c_samples(samp) = c_latest;
  }

  void Mcmc::run()
  {
    Tracer_Plus trace("Mcmc::run");
    
    int samples = 1;
    int jumps = 0;
    int subsamplejumps = 0;
    int relax = false;

    while(true)
      {
	jumps++;
	subsamplejumps++;
	    
  	jump(relax);

	if(subsamplejumps >= opts.sampleevery.value())
	  {
	    if(!switched && c_latest<1e-2)
	      {
		prior_dominating = true;
		break;
	      }

	    subsamplejumps = 0;
	    
	    // sample components after burnin
	    if(jumps > opts.burnin.value())
	      {	   		
		//		if(!relax) relax = true;
		sample(samples);
		samples++;
		
		if(samples>=nsamples)
		  break;
	      }
	  }
      }    
  }

}
