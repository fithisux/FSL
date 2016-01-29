/*  fwdmodel_asl_pvc.cc - Partial Volume Correction resting state ASL model (Buxton)

    Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2009 University of Oxford  */

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

#include "fwdmodel_asl_pvc.h"

#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include "newimage/newimageall.h"
#include "miscmaths/miscprob.h"
using namespace NEWIMAGE;
#include "easylog.h"

string ASL_PVC_FwdModel::ModelVersion() const
{
  return "$Id: fwdmodel_asl_pvc.cc,v 1.6 2013/09/03 15:08:04 chappell Exp $";
}

void ASL_PVC_FwdModel::HardcodedInitialDists(MVNDist& prior, 
    MVNDist& posterior) const
{
    Tracer_Plus tr("ASL_PVC_FwdModel::HardcodedInitialDists");
    assert(prior.means.Nrows() == NumParams());

     SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;

    // Set priors
    // Tissue bolus perfusion
     if (infertiss) {
     prior.means(tiss_index()) = 0;
     precisions(tiss_index(),tiss_index()) = 1e-12;

     
     //if (!singleti) {
       // Tissue bolus transit delay
       prior.means(tiss_index()+1) = setdelt;
       precisions(tiss_index()+1,tiss_index()+1) = 10;
       // }
    
     }

    // Tissue bolus length
     if (infertau && infertiss) {
       prior.means(tau_index()) = seqtau;
       precisions(tau_index(),tau_index()) = 10;
     }

     if (infertaub)
       {
	 prior.means(taub_index()) = seqtau;
	 precisions(taub_index(),taub_index()) = 10;
       }

    // Arterial Perfusion & bolus delay

    if (inferart)
      {
	int aidx = art_index();
	prior.means(aidx) = 0;
	prior.means(aidx+1) = 0.5;
	precisions(aidx+1,aidx+1) = 10;
	precisions(aidx,aidx) = 1e-12;
      }
 
    // T1 & T1b
    if (infert1) {
      int tidx = t1_index();
      prior.means(tidx) = t1;  
      prior.means(tidx+1) = t1b; 
      precisions(tidx,tidx) = 100;
      precisions(tidx+1,tidx+1) = 100;
    }

    /* if (inferart) {
      prior.means(R_index()) = log(10);
      precisions(R_index(),R_index()) = 1;
      }*/

    if (inferwm) {
      int wmi = wm_index();
      prior.means(wmi) = 0;
      prior.means(wmi+1) = 1.2;
      precisions(wmi,wmi) = 1e-12;
      precisions(wmi+1,wmi+1) = 10;

      if (infertau) {
	prior.means(wmi+2) = seqtau;
	precisions(wmi+2,wmi+2) = 10;
      }

      if (infert1) {
	prior.means(wmi+3) = t1wm;
	precisions(wmi+3,wmi+3) = 100;
      }

      if (usepve) {
      //PV entries, the means get overwritten elsewhere if the right sort of prior is specified
      // default is to allow both (NB artifically defies sum(pve)=1)
      int pvi= pv_index();
      prior.means(pvi) = 1; //GM is first
      prior.means(pvi+1)= 1; //WM is second

      // (precisions are big as we treat PV parameters as correct
      // NB they are not accesible from the data anyway)
      //std dev of 1%
      precisions(pvi,pvi) = 1e4;
      precisions(pvi+1,pvi+1) = 1e4;
      }

    }

    /*    if (inferinveff) {
      prior.means(inveff_index()) = 0.3;
      precisions(inveff_index(),inveff_index()) = 10;
      }*/


    // Set precsions on priors
    prior.SetPrecisions(precisions);
    
    // Set initial posterior
    posterior = prior;

    // For parameters with uniformative prior chosoe more sensible inital posterior
    // Tissue perfusion
    if (infertiss) {
    posterior.means(tiss_index()) = 10;
    precisions(tiss_index(),tiss_index()) = 1;
    }
    // Arterial perfusion
    if (inferart)
      {
	posterior.means(art_index()) = 10;
	precisions(art_index(),art_index()) = 1;
      }

    if (inferwm)
      {
	posterior.means(wm_index()) = 10;
	precisions(wm_index(),wm_index()) = 1;
      }
    posterior.SetPrecisions(precisions);
    
}    
    
    

void ASL_PVC_FwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
  Tracer_Plus tr("ASL_PVC_FwdModel::Evaluate");

    // ensure that values are reasonable
    // negative check
  ColumnVector paramcpy = params;
   for (int i=1;i<=NumParams();i++) {
      if (params(i)<0) { paramcpy(i) = 0; }
      }
  
     // sensible limits on transit times
   if (infertiss) {
  if (params(tiss_index()+1)>timax-0.2) { paramcpy(tiss_index()+1) = timax-0.2; }
   }
  if (inferart) {
    if (params(art_index()+1)>timax-0.2) { paramcpy(art_index()+1) = timax-0.2; }
  }


  // parameters that are inferred - extract and give sensible names
  float ftiss;
  float delttiss;
  float tauset; //the value of tau set by the sequence (may be effectively infinite)
  float taubset;
  float fblood;
  float deltblood;
  float T_1;
  float T_1b;

  float pv_gm;
  float pv_wm;

  float fwm;
  float deltwm;
  float tauwmset;
  float T_1wm;

  //  float RR;
  //  float inveffslope;
  //float trailingperiod;

  if (infertiss) {
  ftiss=paramcpy(tiss_index());
  //if (!singleti) {
  delttiss=paramcpy(tiss_index()+1);
  //}
  //else {
    //only inferring on tissue perfusion, assume fixed value for tissue arrival time
    //delttiss = 0;
    //}
  }
  else {
    ftiss=0;
    delttiss=0;
  }

  if (infertau && infertiss) { 
    tauset=paramcpy(tau_index()); 
  }
  else { tauset = seqtau;
  }
  
  if (infertaub) {
    taubset = paramcpy(taub_index());
  }
  else
    { taubset = tauset; }

  if (inferart) {
    fblood=paramcpy(art_index());
    deltblood=paramcpy(art_index()+1);
  }
  else {
    fblood = 0;
    deltblood = 0;
  }

  if (infert1) {
    T_1 = paramcpy(t1_index());
    T_1b = paramcpy(t1_index()+1);

    //T1 cannot be zero!
    if (T_1<0.01) T_1=0.01;
    if (T_1b<0.01) T_1b=0.01;
  }
  else {
    T_1 = t1;
    T_1b = t1b;
  }

  /*if (inferart) {
    RR = exp( paramcpy(R_index()) );
    if (RR<1) RR=1;
    }*/

  if (inferwm) {
    fwm=paramcpy(wm_index());
    //fwm=20;
    deltwm=paramcpy(wm_index()+1);

    if (infertau) {
      tauwmset = paramcpy(wm_index()+2);
    }
    else tauwmset = seqtau;

    if (infert1) {
      T_1wm = paramcpy(wm_index()+3);
      if (T_1<0.01) T_1=0.01;
    }
    else T_1wm = t1wm;

    if (usepve) {
      pv_gm = paramcpy(pv_index());
      pv_wm = paramcpy(pv_index()+1);
    }
    else {
      pv_gm=1;pv_wm=1;
    }
  }
  else {
    fwm=0;
    deltwm=0;
    T_1wm=t1wm;

    pv_gm=1;
    pv_wm=1;
  }


    float lambdagm = 0.98;
    float lambdawm = 0.82;

    float T_1app = 1/( 1/T_1 + 0.01/lambdagm );
    float T_1appwm = 1/( 1/T_1wm + 0.003/lambdawm );
    float R = 1/T_1app - 1/T_1b;
    float Rwm = 1/T_1appwm - 1/T_1b;

    float tau; //bolus length as seen by kintic curve
    float taub; //bolus length of blood as seen in signal
    float tauwm;
    

    float F=0;
    float Fwm=0;
    float kctissue;
    float kcblood;
    float kcwm;


    // loop over tis
    float ti;
    result.ReSize(tis.Nrows()*repeats);

    for(int it=1; it<=tis.Nrows(); it++)
      {
	ti = tis(it) + slicedt*coord_z; //account here for an increase in the TI due to delays between slices;

	if (casl)  {
	  F = 2*ftiss;
	  Fwm = 2*fwm;
	}
	else {
	  F = 2*ftiss * exp(-ti/T_1app);
	  Fwm = 2*fwm * exp(-ti/T_1appwm);
	}
	

	  //GRASE -  deal with bolus length (see above) */

	// Deal with saturation of the bolus before the TI - defined by pretisat
	if(tauset < ti - pretisat)
	  { tau = tauset; }
	else
	  { tau = ti -  pretisat; }
	
	if(taubset < ti -  pretisat)
	  {taub = taubset; }
	else
	  {taub = ti -  pretisat; }

	if(tauwmset < ti -  pretisat)
	  {tauwm = tauwmset; }
	else
	  {tauwm = ti -  pretisat; }

	    

	// (1) tissue contribution 
	    if(ti < delttiss)
	      { kctissue = 0;}
	    else if(ti >= delttiss && ti <= (delttiss + tau))
	      {
		
		if (casl)  kctissue = F * T_1app * exp(-delttiss/T_1b) * (1 - exp(-(ti-delttiss)/T_1app));
		else	   kctissue = F/R * ( (exp(R*ti) - exp(R*delttiss)) ) ;

	      }
	    else //(ti > delttiss + tau)
	      {
		if (casl)  kctissue = F * T_1app * exp(-delttiss/T_1b) * exp(-(ti-tau-delttiss)/T_1app) * (1 - exp(-tau/T_1app));
		else       kctissue = F/R * ( (exp(R*(delttiss+tau)) - exp(R*delttiss))  );
			      }
	
	    // (2) arterial contribution
	    if(ti < deltblood)
	      { 
		//kcblood = 0;
		// use a artifical lead in period for arterial bolus to improve model fitting
		kcblood = fblood * exp(-deltblood/T_1b) * (0.98 * exp( (ti-deltblood)/0.05 ) + 0.02 * ti/deltblood );
	      }
	    else if(ti >= deltblood && ti <= (deltblood + taub))
	      { 
		if (casl)  kcblood = fblood * exp(-deltblood/T_1b);
		else       kcblood = fblood * exp(-ti/T_1b); 
	      }
	    else //(ti > deltblood + tau)
	      {
		kcblood = 0; //end of bolus
		if (casl)  kcblood = fblood * exp(-deltblood/T_1b);
		else	   kcblood = fblood * exp(-(deltblood+taub)/T_1b);
		kcblood *= (0.98 * exp( -(ti - deltblood - taub)/0.05) + 0.02 * (1-(ti - deltblood - taub)/5));
		// artifical lead out period for taub model fitting
		if (kcblood<0) kcblood=0; //negative values are possible with the lead out period equation
									   
	      }
	    // full model for arterial cpt
	    /*	    if(ti < deltblood)
	      { 
		kcblood = 0;
	      }
	    else if(ti >= deltblood && ti <= (deltblood + taub))
	      { 
		
		kcblood = fblood * exp(-ti/T_1b) * (1 - exp( -RR*(ti-deltblood) ) ); 
		
	      }
	    else //(ti > deltblood + tau)
	      {
		kcblood = 0; //end of bolus
		
									   
		}*/

	    // (3) WM contribution
	    if(ti < deltwm)
	      { kcwm = 0;}
	    else if(ti >= deltwm && ti <= (deltwm + tauwm))
	      {
		if (casl)  kcwm = Fwm * T_1appwm * exp(-deltwm/T_1b) * (1 - exp(-(ti-deltwm)/T_1appwm));
		else	   kcwm = Fwm/Rwm * ( (exp(Rwm*ti) - exp(Rwm*deltwm)) ) ;

	      }
	    else //(ti > delttiss + tau)
	      {
		if (casl)  kcwm = Fwm * T_1appwm * exp(-deltwm/T_1b) * exp(-(ti-tauwm-deltwm)/T_1appwm) * (1 - exp(-tauwm/T_1appwm));
		else       kcwm = Fwm/Rwm * ( (exp(Rwm*(deltwm+tauwm)) - exp(Rwm*deltwm))  );
		
	      }



	    if (isnan(kctissue)) { kctissue=0; LOG << "Warning NaN in tissue curve at TI:" << ti << " with f:" << ftiss << " delt:" << delttiss << " tau:" << tau << " T1:" << T_1 << " T1b:" << T_1b << endl; }
	    if (isnan(kcwm)) { kcwm=0; LOG << "Warning NaN in WM curve at TI:" << ti << " with f:" << fwm << " delt:" << deltwm << " tau:" << tauwm << " T1wm:" << T_1wm << " T1b:" << T_1b << endl; }
	    //}

	/* output */
	// loop over the repeats
	for (int rpt=1; rpt<=repeats; rpt++)
	  {
	    result( (it-1)*repeats+rpt ) = pv_gm*kctissue + kcblood + pv_wm*kcwm;
	  }

 
      }
    //cout << result.t();
    

  return;
}


ASL_PVC_FwdModel::ASL_PVC_FwdModel(ArgsType& args)
{
    string scanParams = args.ReadWithDefault("scan-params","cmdline");
    
    if (scanParams == "cmdline")
    {
      // specify command line parameters here
      repeats = convertTo<int>(args.Read("repeats")); // number of repeats in data
      t1 = convertTo<double>(args.ReadWithDefault("t1","1.3"));
      t1b = convertTo<double>(args.ReadWithDefault("t1b","1.5"));
      t1wm = convertTo<double>(args.ReadWithDefault("t1wm","1.1"));
      lambda = convertTo<double>(args.ReadWithDefault("lambda","0.9")); //NOT used - here for compatibility

      pretisat = convertTo<double>(args.ReadWithDefault("pretisat","0")); // deal with saturation of the bolus a fixed time pre TI measurement
      grase = args.ReadBool("grase"); // DEPRECEATED data has come from the GRASE-ASL sequence - therefore apply pretisat of 0.1s
      if (grase) pretisat=0.1;

      casl = args.ReadBool("casl"); //set if the data is CASL or PASL (default)
      slicedt = convertTo<double>(args.ReadWithDefault("slicedt","0.0")); // increase in TI per slice

      infertau = args.ReadBool("infertau"); // infer on bolus length?
      infert1 = args.ReadBool("infert1"); //infer on T1 values?
      inferart = args.ReadBool("inferart"); //infer on arterial compartment?
      inferwm = args.ReadBool("inferwm");
      //inferinveff = args.ReadBool("inferinveff"); //infer on a linear decrease in inversion efficiency?
      //infertrailing = args.ReadBool("infertrailing"); //infers a trailing edge bolus slope using new model
      seqtau = convertTo<double>(args.ReadWithDefault("tau","1000")); //bolus length as set by sequence (default of 1000 is effectively infinite
      setdelt = convertTo<double>(args.ReadWithDefault("bat","0.7"));

      bool ardoff = false;
      ardoff = args.ReadBool("ardoff");
      bool tauboff = false;
      tauboff = args.ReadBool("tauboff"); //forces the inference of arterial bolus off

      usepve = args.ReadBool("usepve");

      // combination options
      infertaub = false;
      if (inferart && infertau && !tauboff) infertaub = true;

      
      //special - turn off tissue cpt
      infertiss=true;
      bool tissoff = args.ReadBool("tissoff");
      if (tissoff) infertiss = false;

      
      // deal with ARD selection
      doard=false;
      tissard=false;artard=true;wmard=true; //default ARD flags
      //if (inferart==true && ardoff==false) { doard=true;}
      //if (inferwm==true && ardoff==false) {doard=true; }
      //special, individual ARD switches
      bool tissardon = args.ReadBool("tissardon");
      if (tissardon) tissard=true;
      bool artardoff = args.ReadBool("artardoff");
      if (artardoff) artard=false;
      bool wmardoff = args.ReadBool("wmardoff");
      if (wmardoff) wmard=false;

      // ** ardoff overrides all other ARD options
      if ( (tissard || artard || wmard) && !ardoff) doard = true;

      /* if (infertrailing) {
	if (!infertau) {
	  // do not permit trailing edge inference without inferring on bolus length
	  throw Invalid_option("--infertrailing has been set without setting --infertau");
	}
	else if (inferinveff)
	  //do not permit trailing edge inference and inversion efficiency inference (they are mututally exclusive)
	  throw Invalid_option("--infertrailing and --inferinveff may not both be set");
	  }*/

      // Deal with tis
      tis.ReSize(1); //will add extra values onto end as needed
      tis(1) = atof(args.Read("ti1").c_str());
      
      while (true) //get the rest of the tis
	{
	  int N = tis.Nrows()+1;
	  string tiString = args.ReadWithDefault("ti"+stringify(N), "stop!");
	  if (tiString == "stop!") break; //we have run out of tis
	 
	  // append the new ti onto the end of the list
	  ColumnVector tmp(1);
	  tmp = convertTo<double>(tiString);
	  tis &= tmp; //vertical concatenation

	}
      timax = tis.Maximum(); //dtermine the final TI
      
  // need to set the voxel coordinates to a deafult of 0 (for the times we call the model before we start handling data)
      coord_x = 0;
      coord_y = 0;
      coord_z = 0;
      
      singleti = false; //normally we do multi TI ASL
      /*if (tis.Nrows()==1) {
	//only one TI therefore only infer on CBF and ignore other inference options
	LOG << "--Single inversion time mode--" << endl;
	LOG << "Only a sinlge inversion time has been supplied," << endl;
	LOG << "Therefore only tissue perfusion will be inferred." << endl;
	LOG << "-----" << endl;
	singleti = true;
	// force other inference options to be false
	infertau = false; infert1 = false; inferart = false; //inferinveff = false;
	}*/
	
      // add information about the parameters to the log
      LOG << "Inference using development model" << endl;
      if (pretisat>0) LOG << "Saturation of" << pretisat << "s before TI has been specified" << endl;
      if (grase) LOG << "Using pre TI saturation of 0.1 for GRASE-ASL sequence" << endl;
      LOG << "    Data parameters: #repeats = " << repeats << ", t1 = " << t1 << ", t1b = " << t1b;
      LOG << ", bolus length (tau) = " << seqtau << endl ;
      if (infertau) {
	LOG << "Infering on bolus length " << endl; }
      if (doard) {
	LOG << "ARD subsystem is enabled" << endl; }
      if (infertiss) {
	LOG << "Infertting on tissue component " << endl; }
      if (doard && tissard) {
	LOG << "ARD has been set on the tissue component " << endl; }
      if (inferart) {
	LOG << "Infering on artertial compartment " << endl; }
      if (doard && artard) {
	LOG << "ARD has been set on arterial compartment " << endl; }
      if (inferwm) {
	LOG << "Inferring on white matter component" << endl; 
	if (doard && wmard) { LOG << "ARD has been set on wm component" << endl;}
      }
      if (infert1) {
	LOG << "Infering on T1 values " << endl; }
      /*if (inferinveff) {
	LOG << "Infering on Inversion Efficency slope " << endl; }
      if (infertrailing) {
      LOG << "Infering bolus trailing edge period" << endl; }*/
      LOG << "TIs: ";
      for (int i=1; i <= tis.Nrows(); i++)
	LOG << tis(i) << " ";
      LOG << endl;
	  
    }

    else
        throw invalid_argument("Only --scan-params=cmdline is accepted at the moment");    
    
 
}

void ASL_PVC_FwdModel::ModelUsage()
{ 
  cout << "\nUsage info for --model=grase:\n"
       << "Required parameters:\n"
       << "--repeats=<no. repeats in data>\n"
       << "--ti1=<first_inversion_time_in_seconds>\n"
       << "--ti2=<second_inversion_time>, etc...\n"
       << "Optional arguments:\n"
       << "--grase *DEPRECEATAED* (data collected using GRASE-ASL: same as --pretissat=0.1)"
       << "--pretisat=<presat_time> (Define that blood is saturated a specific time before TI image acquired)"
       << "--tau=<temporal_bolus_length> (default 10s if --infertau not set)\n"
       << "--t1=<T1_of_tissue> (default 1.3)\n"
       << "--t1b=<T1_of_blood> (default 1.5)\n"
       << "--infertau (to infer on bolus length)\n"
       << "--inferart (to infer on arterial compartment)\n"
       << "--infert1 (to infer on T1 values)\n"
    ;
}

void ASL_PVC_FwdModel::DumpParameters(const ColumnVector& vec,
                                    const string& indent) const
{
    
}

void ASL_PVC_FwdModel::NameParams(vector<string>& names) const
{
  names.clear();
  
  if (infertiss) {
  names.push_back("ftiss");
  //if (!singleti) 
    names.push_back("delttiss");
  }
  if (infertau && infertiss)
    {
    names.push_back("tautiss");
    }
  if (inferart) {
    names.push_back("fblood");
    names.push_back("deltblood");
  }
  if (infert1) {
    names.push_back("T_1");
    names.push_back("T_1b");
  }
  /* if (inferinveff) {
    names.push_back("Inveffslope");
  }
  if (infertrailing) {
    names.push_back("trailingperiod");
    }*/
  if (infertaub) {
    names.push_back("taublood");
  }
  /*if (inferart) {
    names.push_back("R");
    }*/

  if (inferwm) {
    names.push_back("fwm");
    names.push_back("deltwm");

    if (infertau) names.push_back("tauwm");
    if (infert1) names.push_back("T_1wm");

    if (usepve) {
      names.push_back("p_gm");
      names.push_back("p_wm");
    }
  }
}

void ASL_PVC_FwdModel::SetupARD( const MVNDist& theta, MVNDist& thetaPrior, double& Fard)
{
  Tracer_Plus tr("ASL_PVC_FwdModel::SetupARD");

  if (doard)
    {
      //sort out ARD indices
      if (tissard) ard_index.push_back(tiss_index());
      if (artard) ard_index.push_back(art_index());
      if (wmard) ard_index.push_back(wm_index());

      Fard = 0;

      int ardindex;
      for (unsigned int i=0; i<ard_index.size(); i++) {
	//iterate over all ARD parameters
	ardindex = ard_index[i];

	SymmetricMatrix PriorPrec;
	PriorPrec = thetaPrior.GetPrecisions();
	
	PriorPrec(ardindex,ardindex) = 1e-12; //set prior to be initally non-informative
	
	thetaPrior.SetPrecisions(PriorPrec);
	
	thetaPrior.means(ardindex)=0;
	
	//set the Free energy contribution from ARD term
	SymmetricMatrix PostCov = theta.GetCovariance();
	double b = 2/(theta.means(ardindex)*theta.means(ardindex) + PostCov(ardindex,ardindex));
	Fard += -1.5*(log(b) + digamma(0.5)) - 0.5 - gammaln(0.5) - 0.5*log(b); //taking c as 0.5 - which it will be!
      }
  }
  return;
}

void ASL_PVC_FwdModel::UpdateARD(
				const MVNDist& theta,
				MVNDist& thetaPrior, double& Fard) const
{
  Tracer_Plus tr("ASL_PVC_FwdModel::UpdateARD");
  
  if (doard)
    Fard=0;
    {
      int ardindex;
      for (unsigned int i=0; i<ard_index.size(); i++) {
	//iterate over all ARD parameters
	ardindex = ard_index[i];

  
      SymmetricMatrix PriorCov;
      SymmetricMatrix PostCov;
      PriorCov = thetaPrior.GetCovariance();
      PostCov = theta.GetCovariance();

      PriorCov(ardindex,ardindex) = theta.means(ardindex)*theta.means(ardindex) + PostCov(ardindex,ardindex);

      
      thetaPrior.SetCovariance(PriorCov);

      //Calculate the extra terms for the free energy
      double b = 2/(theta.means(ardindex)*theta.means(ardindex) + PostCov(ardindex,ardindex));
      Fard += -1.5*(log(b) + digamma(0.5)) - 0.5 - gammaln(0.5) - 0.5*log(b); //taking c as 0.5 - which it will be!
    }
  }

  return;

  }
