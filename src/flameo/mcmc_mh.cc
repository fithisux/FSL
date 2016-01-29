/*  mcmc_mh.cc

    Mark Woolrich, Tim Behrens, FMRIB Image Analysis Group

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

#include "mcmc_mh.h"
#include "utils/log.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include "newimage/newimageall.h"
#include "utils/tracer_plus.h"
#include <set>

using namespace Utilities;
using namespace MISCMATHS;
using namespace NEWIMAGE;

namespace Gs {


  void Mcmc_Mh::setup()
  {
    Tracer_Plus trace("Mcmc_Mh::setup");
    
    beta_naccepted = 0;
    phi_naccepted = 0;
    gamma_naccepted = 0;
    beta_nrejected = 0;
    phi_nrejected = 0;
    gamma_nrejected = 0;

    gamma_latest.ReSize(nevs);
    gamma_samples.ReSize(nevs,nsamples);
    beta_latest.ReSize(ngs);
    beta_samples.ReSize(ngs,nsamples);
    phi_latest.ReSize(ntpts);
    phi_samples.ReSize(ntpts,nsamples);
    likelihood_samples.ReSize(nsamples);
    likelihood_samples = 0;

    Matrix gamC;
    reshape(gamC, gamma_S, nevs, nevs); // gamC is Covariance
    
    gamma_latest = gamma_mean;
    gamma_proposal_std = sqrt(diag(gamC))*8;
    gamma_samples = 0;

    for(int g = 1; g <= ngs; g++)
      {
//      // use when sampling from precision
// 	  beta_latest(g) = 1.0/(beta_b(g)/beta_c(g));

// 	beta_proposal_std(g) = 4.0/(sqrt(beta_b(g)/Sqr(beta_c(g))));

// 	// use when sampling from variance
	beta_latest(g) = beta_b(g)/beta_c(g);

	beta_proposal_std(g) = 4.0*(sqrt(beta_b(g)/Sqr(beta_c(g))));

	// use when sampling from log(variance)
// 	  beta_latest(g) = log((beta_b(g)/beta_c(g)));

// 	beta_proposal_std(g) = fabs(log(4.0*(sqrt(beta_b(g)/Sqr(beta_c(g))))));

      }

    for(int t = 1; t <= ntpts; t++)
      {
	phi_latest(t) = varcopedata(t);

	if(opts.dofvarcopefile.value() != string(""))
	  phi_proposal_std(t) = sqrt(2*Sqr(varcopedata(t))/dofvarcopedata(t))/2;
      }

    beta_samples = 0;
    phi_samples = 0;
 
    for(int g = 1; g <= ngs; g++)
      {
	beta_prior_energy_old(g) = beta_prior_energy(g);
      }
    
    if(uncertainty_in_varcopes)
      {
	for(int t = 1; t <= ntpts; t++)
	{
	  phi_prior_energy_old(t) = phi_prior_energy(t);
	}
      }
    sumovere = 0;

    for(int t = 1; t <= ntpts; t++)
      {
	prec_ontwo(t) = 0.5/(phi_latest(t)+beta_latest(design.getgroup(t)));
	logprec_ontwo(t) = log(2*prec_ontwo(t))/2.0;

	for(int e = 1; e <= nevs; e++)
	  {
	    sumovere(t) += design_matrix(t,e)*gamma_latest(e);
	    // 	OUT(gamma_latest(e));
	    // 	OUT(design_matrix(t,e));
	  }	
      }

    likelihood_energy_old = likelihood_energy(0,0,false);

    // if all_jump
//     likelihood_energy_old = 0;
//     for(int g = 1; g <= ngs; g++)
//       {
// 	likelihood_energy_old += beta_prior_energy(g);
//       }
    
//     likelihood_energy_old += likelihood_energy();

  }

  void Mcmc_Mh::jump()
  {
    Tracer_Plus trace("Mcmc_Mh::jump");
        
    //    all_jump();

    beta_jump();

    if(uncertainty_in_varcopes)
      phi_jump();

    gamma_jump();

    if(subsampcount>40)
      {
	  for(int g = 1; g <= ngs; g++)
	    beta_proposal_std(g) *= 0.65/((1+beta_nrejected(g))/float(1+beta_naccepted(g)+beta_nrejected(g)));
	
	if(uncertainty_in_varcopes)
	  for(int t = 1; t <= ntpts; t++)
	    phi_proposal_std(t) *= 0.65/((1+phi_nrejected(t))/float(1+phi_naccepted(t)+phi_nrejected(t)));
	
	for(int e = 1; e <= nevs; e++)
	  gamma_proposal_std(e) *= 0.65/((1+gamma_nrejected(e))/float(1+gamma_naccepted(e)+gamma_nrejected(e)));
	
	beta_naccepted = 0;
	phi_naccepted = 0;
	gamma_naccepted = 0;
	beta_nrejected = 0;
	phi_nrejected = 0;
	gamma_nrejected = 0;

	subsampcount = 0;
      }
    else
      {	
	subsampcount++;
      }

  } 

  void Mcmc_Mh::beta_jump()
  {
    Tracer_Plus trace("Mcmc_Mh::beta_jump");

//     if(sampcount>9844)
//        GsOptions::getInstance().debuglevel.set_value("2");

    for(int g = 1; g <= ngs; g++)
      {
	// store old values
	float beta_old = beta_latest(g);
	ColumnVector prec_ontwo_old = prec_ontwo;
	ColumnVector logprec_ontwo_old = logprec_ontwo;

	// propose new value	
	beta_latest(g) += normrnd().AsScalar()*beta_proposal_std(g);	  

	// use when sampling from log(variance)
	//	if(abs(beta_latest(g)) > 50) {beta_latest(g) = beta_old; beta_nrejected(g)++; return;}

	// use when sampling from variance
	if(beta_latest(g) <= 0) {beta_latest(g) = beta_old; beta_nrejected(g)++; return;}
	
	float likelihood_energy_new = likelihood_energy(0,0,true);	
	float beta_prior_energy_new = beta_prior_energy(g);

	// calculate acceptance threshold
	float tmp = unifrnd().AsScalar();
	float energy_new = likelihood_energy_new + beta_prior_energy_new;
	float energy_old = likelihood_energy_old + beta_prior_energy_old(g);

	if(GsOptions::getInstance().debuglevel.value()==2)
	  {
	    cout << "--------------" << endl;	
	    OUT(varcopedata.t());
	    OUT(copedata.t());
	    OUT(sampcount);
	    OUT(gamma_latest.t());	    
	    OUT(beta_latest(g));
	    OUT(beta_proposal_std(g));
	    OUT(beta_old);
	    OUT(beta_prior_energy_new);
	    OUT(beta_prior_energy_old(g));
	    OUT(likelihood_energy_new);
	    OUT(likelihood_energy_old);
	    OUT(energy_new);
	    OUT(energy_old);
	    OUT(tmp);
	    OUT(exp(energy_old - energy_new));
	  }

	bool accept = exp(energy_old - energy_new) > tmp;
	
	if(accept)
	  {
	    if(GsOptions::getInstance().debuglevel.value()==2)
	      {
		cout << "accepted" << endl;
	      }
	    beta_prior_energy_old(g) = beta_prior_energy_new;
	    likelihood_energy_old = likelihood_energy_new;
	    beta_naccepted(g)++;
	  }
	else
	  {
	    if(GsOptions::getInstance().debuglevel.value()==2)
	      {
		cout << "rejected" << endl;
	      }
	    // restore old values
	    beta_latest(g) = beta_old;

	    prec_ontwo = prec_ontwo_old;
	    logprec_ontwo = logprec_ontwo_old;

	    beta_nrejected(g)++;
	  }
      }
  }

  void Mcmc_Mh::phi_jump()
  {
    Tracer_Plus trace("Mcmc_Mh::phi_jump");

//     if(sampcount>9844)
//        GsOptions::getInstance().debuglevel.set_value("2");

    for(int t = 1; t <= ntpts; t++)
      {
	// store old values
	float phi_old = phi_latest(t);
	ColumnVector prec_ontwo_old = prec_ontwo;
	ColumnVector logprec_ontwo_old = logprec_ontwo;

	// propose new value	
	phi_latest(t) += normrnd().AsScalar()*phi_proposal_std(t);	  

	if(phi_latest(t) <= 0) {phi_latest(t) = phi_old; phi_nrejected(t)++; return;}
	
	float likelihood_energy_new = likelihood_energy_phichanged(t);	
	float phi_prior_energy_new = phi_prior_energy(t);
	float beta_prior_energy_new = beta_prior_energy(design.getgroup(t));

	// calculate acceptance threshold
	float tmp = unifrnd().AsScalar();
	float energy_new = likelihood_energy_new + phi_prior_energy_new + beta_prior_energy_new;
	float energy_old = likelihood_energy_old + phi_prior_energy_old(t) + beta_prior_energy_old(design.getgroup(t));

	if(GsOptions::getInstance().debuglevel.value()==2)
	  {
	    cout << "--------------" << endl;	
	    OUT(varcopedata.t());
	    OUT(copedata.t());
	    OUT(sampcount);
	    OUT(gamma_latest.t());	    
	    OUT(phi_latest(t));
	    OUT(phi_proposal_std(t));
	    OUT(phi_old);
	    OUT(phi_prior_energy_new);
	    OUT(phi_prior_energy_old(t));
	    OUT(likelihood_energy_new);
	    OUT(likelihood_energy_old);
	    OUT(energy_new);
	    OUT(energy_old);
	    OUT(tmp);
	    OUT(exp(energy_old - energy_new));
	  }

	bool accept = exp(energy_old - energy_new) > tmp;
	
	if(accept)
	  {
	    if(GsOptions::getInstance().debuglevel.value()==2)
	      {
		cout << "accepted" << endl;
	      }
	    phi_prior_energy_old(t) = phi_prior_energy_new;
	    beta_prior_energy_old(design.getgroup(t)) = beta_prior_energy_new;
	    likelihood_energy_old = likelihood_energy_new;
	    phi_naccepted(t)++;
	  }
	else
	  {
	    if(GsOptions::getInstance().debuglevel.value()==2)
	      {
		cout << "rejected" << endl;
	      }
	    // restore old values
	    phi_latest(t) = phi_old;
	    prec_ontwo = prec_ontwo_old;
	    logprec_ontwo = logprec_ontwo_old;

	    phi_nrejected(t)++;
	  }
      }
  }

//   void Mcmc_Mh::gamma_jump()
//   {
//     Tracer_Plus trace("Mcmc_Mh::gamma_jump");

//     // store old values
//     // propose new values
//     ColumnVector gamma_old(nevs);
//     for(int e = 1; e <= nevs; e++)
//       {
// 	gamma_old(e) = gamma_latest(e);       
// 	gamma_latest(e) += normal.Next()*gamma_proposal_std(e);
//       }
    
//     float likelihood_energy_new = likelihood_energy();	

//     // calculate acceptance threshold
//     float tmp = uniform.Next();
//     bool accept = exp(likelihood_energy_old - likelihood_energy_new) > tmp;
    
//     if(accept)
//       {
// 	likelihood_energy_old = likelihood_energy_new;
// 	gamma_naccepted(1)++;
//       }
//     else
//       {
// 	// restore old values
// 	for(int e = 1; e <= nevs; e++)
// 	  {
// 	    gamma_latest(e) = gamma_old(e);
// 	  }
	
// 	gamma_nrejected(1)++;
//       }
//   }

  void Mcmc_Mh::gamma_jump()
  {
    Tracer_Plus trace("Mcmc_Mh::gamma_jump");
    
    for(int e = 1; e <= nevs; e++)
      {
	// store old values
	float gamma_old = gamma_latest(e);       
	ColumnVector sumovere_old = sumovere;

	// propose new values	
	gamma_latest(e) += normrnd().AsScalar()*gamma_proposal_std(e);      

	float likelihood_energy_new = likelihood_energy(e,gamma_old,false);	

	// calculate acceptance threshold

	float tmp = unifrnd().AsScalar();

	bool accept = exp(likelihood_energy_old - likelihood_energy_new) > tmp;

	if(accept)
	  {
	    likelihood_energy_old = likelihood_energy_new;
	    gamma_naccepted(e)++;
	  }
	else
	  {
	    // restore old values	    
	    gamma_latest(e) = gamma_old;
	     sumovere = sumovere_old;
	    gamma_nrejected(e)++;
	  }
      } 
  }

//   void Mcmc_Mh::all_jump()
//   {
//     Tracer_Plus trace("Mcmc_Mh::all_jump");

//     // store old values
//     // propose new values
//     ColumnVector gamma_old(nevs);
//     for(int e = 1; e <= nevs; e++)
//       {
// 	gamma_old(e) = gamma_latest(e);       
// 	gamma_latest(e) += normrnd().AsScalar()*gamma_proposal_std(e);
//       }

//     ColumnVector beta_old(ngs);
//     float energy_new = 0.0;

//     for(int g = 1; g <= ngs; g++)
//       {
// 	beta_old(g) = beta_latest(g);
// 	beta_latest(g) += normrnd().AsScalar()*beta_proposal_std(g);	  
// 	if(beta_latest(g) <= 0) {beta_latest(g) = beta_old(g);}
// 	energy_new += beta_prior_energy(g);
//       }
    
//     energy_new += likelihood_energy();	

//     // calculate acceptance threshold
//     float tmp = unifrnd().AsScalar();
//     bool accept = exp(likelihood_energy_old - energy_new) > tmp;
    
//     if(accept)
//       {
// 	likelihood_energy_old = energy_new;
// 	gamma_naccepted(1)++;
//       }
//     else
//       {
// 	// restore old values
// 	for(int e = 1; e <= nevs; e++)
// 	  {
// 	    gamma_latest(e) = gamma_old(e);
// 	  }

// 	for(int g = 1; g <= ngs; g++)
// 	  {
// 	    beta_latest(g) = beta_old(g);
// 	  }
	
// 	gamma_nrejected(1)++;
//       }
//   }

  float Mcmc_Mh::likelihood_energy(const int echanged, const float gamma_old, const bool betachanged)
  {
    Tracer_Plus trace("Mcmc_Mh::likelihood_energy");

    float en = 0.0;

//     OUT(copedata);
//     OUT(varcopedata);
//     OUT(design_matrix);
    float gamlatest = 0;

    if(echanged>0) gamlatest= gamma_latest(echanged);

    // matlab:   n = 4;y = ones(n)*o+eye(n)*k;inv(y),det(y)    
    for(int t = 1; t <= ntpts; t++)
      {
	//float sumovere = 0.0;

// 	OUT(t);

	if(echanged>0)
	  {
	    sumovere(t) += (gamlatest-gamma_old)*design_matrix(t,echanged);
	  }

	// use when sampling from variance
	if(betachanged)
	  {
	    if(!infer_outliers)
	      {
		prec_ontwo(t) = 0.5/(phi_latest(t)+beta_latest(design.getgroup(t)));				
	      }
	    else
	      {
		float vr=Sqr(1-prob_outlier(t))*(phi_latest(t)+beta_latest(design.getgroup(t)))+Sqr(prob_outlier(t))*(phi_latest(t)+beta_latest(design.getgroup(t))+beta_outlier[design.getgroup(t)-1]);
		prec_ontwo(t) = 0.5/vr;
	      }
	    logprec_ontwo(t) = log(2*prec_ontwo(t))/2.0;
	  }

	// for inferring outliers we do not need the log(l_kf + (1-l_k)(1-f)) term as it is a constant wrt mcmc params
	en += -logprec_ontwo(t) + prec_ontwo(t)*Sqr(copedata(t) - sumovere(t));

	//float prec = 1.0/(phi_latest(t)+beta_latest(design.getgroup(t)));
	//en += -0.5*log(prec)+0.5*prec*Sqr(copedata(t) - sumovere(t));

	// use when sampling from log(variance)
// 	if(betachanged)
// 	  {
// 	    prec_ontwo(t) = 0.5/(phi_latest(t)+exp(beta_latest(design.getgroup(t))));
// 	    logprec_ontwo(t) = log(2*prec_ontwo(t))/2.0;
// 	  }
// 	en += -logprec_ontwo(t) + prec_ontwo(t)*Sqr(copedata(t) - sumovere(t));

// 	OUT(varcopedata(t));
// 	OUT(beta_latest(design.getgroup(t)));
// 	OUT(prec);
// 	OUT(design.getgroup(t));	
	  
      }

    //    OUT(energy);
    return en;
  }
  
  float Mcmc_Mh::likelihood_energy_phichanged(const int t)
  {    
    // logprec_ontwo(t) and prec_ontwo(t) will be calculated using phi_old
    float old_energy = (-logprec_ontwo(t) + prec_ontwo(t)*Sqr(copedata(t) - sumovere(t)));

    // recalculate logprec_ontwo(t) and prec_ontwo(t) using phi_latest
    prec_ontwo(t) = 0.5/(phi_latest(t)+beta_latest(design.getgroup(t)));
    logprec_ontwo(t) = log(2*prec_ontwo(t))/2.0;
    float new_energy = (-logprec_ontwo(t) + prec_ontwo(t)*Sqr(copedata(t) - sumovere(t)));

    float en = likelihood_energy_old - old_energy + new_energy;

    return en;
  }

  float Mcmc_Mh::beta_prior_energy(int g)
  {
    Tracer_Plus trace("Mcmc_Mh::beta_prior_energy");
    
    float en = 0.0;
    
    // prior is 1/beta (beta is variance)
    en = log(beta_latest(g));

    return en;
  }

  float Mcmc_Mh::phi_prior_energy(int t)
  {
    Tracer_Plus trace("Mcmc_Mh::phi_prior_energy");
    
    // p276 Lee
    float S = dofvarcopedata(t)/varcopedata(t);
    float en = -(dofvarcopedata(t)/2-1)*log(phi_latest(t)) + 0.5*S*phi_latest(t);

    return en;
  }

  void Mcmc_Mh::sample(int samp)
  {
    Tracer_Plus trace("Mcmc_Mh::sample");

    sampcount++;    

    for(int g = 1; g <= ngs; g++)
      beta_samples(g,samp) = beta_latest(g);

    for(int t = 1; t <= ntpts; t++)
      phi_samples(t,samp) = phi_latest(t);

    for(int e = 1; e <= nevs; e++)
      gamma_samples(e,samp) = gamma_latest(e);	

    likelihood_samples(samp) = likelihood_energy_old;

    //    sample_sumsquares(samp);
  }

//   void Mcmc_Mh::sample_sumsquares(int samp)
//   {
//     Tracer_Plus trace("Mcmc_Mh::sample_sumsquares");
    
//     ss_samples[0](samp) = sumsquare_residuals(design_matrix,copedata,gamma_latest);
    
//     for(int f = 1; f < design.getnumfcontrasts()+1; f++)
//       {	
// 	const Matrix& reduceddm = design.getfreduceddm(f);
	
// // 	OUT(copedata.t());
// // 	OUT(gamma_latest.t());
// // 	OUT(reduceddm);

// 		ss_samples[f](samp) = sumsquare_residuals(reduceddm,copedata,gamma_latest);
//       }
//   }
  
//   float Mcmc_Mh::sumsquare_residuals(const Matrix& pdm, const ColumnVector& pdata, const ColumnVector& ppes)
//   {
//     Tracer_Plus trace("Mcmc_Mh::sumsquare_residuals");

//     float ss = 0.0;

//     for(int t = 1; t <= ntpts; t++)
//       {
// 	float sumovere = 0.0;
	
// 	for(int e = 1; e <= nevs; e++)
// 	  {
// 	    sumovere += pdm(t,e)*ppes(e);
// 	  }
// 	float prec = 1.0/(varcopedata(t)+exp(beta_latest(design.getgroup(t))));

// 	ss += prec*Sqr(pdata(t) - sumovere);
// 	//ss += Sqr(pdata(t) - sumovere);
	
//       }

//     ss = ss/ntpts;

//     return ss;
//   }
  
//   void Mcmc_Mh::dic(float& DIC, float& pd)
//   {       
//     Tracer_Plus trace("Mcmc_Mh::dic");

    // calc Dthetabar
    // set latest params to posterior means
//     gamma_latest = mean(gamma_samples,2);
//     beta_latest = mean(beta_samples,2); 

//     float Dthetabar = 2*likelihood_energy();

//     // calc Dbar
//     float Dbar = mean(2*likelihood_samples).AsScalar();
    
//     // pd = Dbar - Dthetabar
//     pd = Dbar - Dthetabar;

//     // DIC = Dbar + pd
//     DIC = Dbar + pd;
//   }

  void Mcmc_Mh::run()
  {
    Tracer_Plus trace("Mcmc_Mh::run");
    
    int samples = 1;
    int jumps = 0;
    int subsamplejumps = 0;

    while(true)
      {
	jumps++;
	subsamplejumps++;
	    
  	jump();

	if(subsamplejumps >= opts.sampleevery.value())
	  {
	    subsamplejumps = 0;
	    
	    // sample components after burnin
	    if(jumps > opts.burnin.value())
	      {	   		
		sample(samples);
		samples++;
		
		if(samples>nsamples)
		  break;
	      }
	  }
      }    
  }

}
