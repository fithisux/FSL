/*  fwdmodel_pcASL.cc - Implements (pseudo) continuous ASL for multi-echo time functional analysis

    Michael Chappell, IBME & FMRIB Image Analysis Group

    Copyright (C) 20011 University of Oxford  */

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

#include "fwdmodel_pcASL.h"

#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include "newimage/newimageall.h"
using namespace NEWIMAGE;
#include "easylog.h"

string pcASLFwdModel::ModelVersion() const
{
  return "$Id: fwdmodel_pcASL.cc,v 1.5 2013/01/25 10:57:32 chappell Exp $";
}

void pcASLFwdModel::HardcodedInitialDists(MVNDist& prior, 
    MVNDist& posterior) const
{
    Tracer_Plus tr("pcASLFwdModel::HardcodedInitialDists");
    assert(prior.means.Nrows() == NumParams());
    
    // Set priors
    prior.means = 0;
    SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;

    if(0) 
      {
	LOG_ERR("Hack: using 1e-10 precision instead of 1e-12");
	precisions = IdentityMatrix(NumParams()) * 1e-10;
      }

    if(0) 
      {
	LOG_ERR("Hack: using slightly informative (100x expected stdev) priors"<<endl);
	precisions(1,1) = 1e-4/350/350;
	precisions(2,2) = 1e-4/350/350;
	precisions(3,3) = 1e-4/17000/17000;
	precisions(4,4) = 1e-4/700/700;
	precisions(5,5) = 1e-4/50/50;
	precisions(6,6) = 1e-4/.05/.05;
      }
    
    // informative prior on dt
    if (stdevDt>0)
      {
        prior.means(dtIndex()) = fixedDt;
        precisions(dtIndex(),dtIndex()) = 1/(stdevDt*stdevDt);
      }
    if (stdevT1b>0)
      {
	    prior.means(T1bIndex()) = fixedT1b;
	    precisions(T1bIndex(),T1bIndex()) = 1/(stdevT1b*stdevT1b);
      }
    if (stdevInvEff>0)
      {
	    prior.means(InvEffIndex()) = fixedInvEff;
	    precisions(InvEffIndex(),InvEffIndex()) = 1/(stdevInvEff*stdevInvEff);
      }
    prior.SetPrecisions(precisions);
    
    // Set informative initial posterior
    posterior = prior;

    if (stdevDt>0)
      {
        precisions(dtIndex(),dtIndex()) = 1/(stdevDt*stdevDt);
      }
    if (stdevT1b>0)
      {
	    precisions(T1bIndex(),T1bIndex()) = 1/(stdevT1b*stdevT1b);
      }
    if (stdevInvEff>0)
      {
	    precisions(InvEffIndex(),InvEffIndex()) = 1/(stdevInvEff*stdevInvEff);
      }
    precisions = IdentityMatrix(NumParams()) * 1;
    posterior.SetPrecisions(precisions); //use quite informative precisions for the initialization of things with uniformative prior

    posterior.means(Q0index()) = 200;
//    assert(id.Q0(posterior.means) == 200);
    posterior.means(M0index()) = 1.5e4;
//    assert(id.M0(posterior.means) == 1.5e4);
    posterior.means(R0index()) = 25;
//    assert(id.R0(posterior.means) == 25);
}    

void pcASLFwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
    Tracer_Plus tr("pcASLFwdModel::Evaluate");

    double R0 = params(R0index());
    if (R0<1) R0=1; //R0 cannot be negative, or very small for the matter
    
    // Parameterization used in most recent results:
    // Absolute M and Q change (same units as M0 or Q0):
    ColumnVector StatMag = params(M0index()) - Mbasis * MnOf(params);
    ColumnVector CBF = params(Q0index()) + Qbasis * QnOf(params);
    // Fractional change in BOLD effect (at TE_2), rather than using % R2* change 
    ColumnVector R2s = -1/echoTime(2) * log( 
        Rbasis * RnOf(params) + exp(-echoTime(2)*R0));

    // The following are relative magnetizations    
    double pretag = 1; // untagged
    double T1b = (stdevT1b>0 ? params(T1bIndex()) : fixedT1b);
    double invEfficiency = (stdevInvEff>0 ? params(InvEffIndex()) : fixedInvEff);
    double dt = (stdevDt>0? params(dtIndex()) : fixedDt);
    ColumnVector bolus = 1 - (1-rho)*invEfficiency*T1b*exp(-TI/T1b)*exp(-dt/T1b)*( exp(-Tau/T1b)-1 ); // tag or control
    double posttag = 1;
    ColumnVector Sb = SP( CBF,  // SP(a,b) means a.*b 
        pretag*dt + bolus*Tau + posttag*(TI-Tau-dt) );
    ColumnVector S = StatMag + Sb;
      
  int Ntimes = R2s.Nrows();
  int Nte = echoTime.Nrows();
  if (result.Nrows() != Nte*Ntimes)
    result.ReSize(Nte*Ntimes);
    
//  result = 0.0/0.0; // pre-fill with nans to check all overwritten
  
  for (int te = 1; te <= Nte; te++)
    {
      ColumnVector nuisance = Nbasis * NnOf(te, params);
      // Will be all-zero if there are no nuisance regressors
        
      for (int i = 1; i <= Ntimes; i++)
        result( Nte*(i-1) + te ) = 
            S(i) * exp(-echoTime(te) * R2s(i)) + nuisance(i);
      // Fill order: te1 te2 te1 te2 te1 te2 te1 te2 ...
    }

  return; // answer is in the "result" vector
}

void pcASLFwdModel::ModelUsage()
{
    cout << "\nUsage info for --model=pcasl-dualecho:\n"
      << "Required options:\n"
      << "--bold-basis=<bold_design_file>\n"
      << "--cbf-basis=<cbf_design_file>\n"
      << "--statmag-basis=<statmag_design_file>\n\n"
      << "Optional options:\n"
      << "--nuisance-basis=<nuisance_regressors_design_file> (default: null)\n"
      << "--ti=<ti_in_sec>, inversion time of acquisition = PLD + Bolus duration (default: 2.0)"
      << "--tau=<tau_in_sec>, bolus duration (default: 1.0)\n"
      << "--te1=<te1_in_millisec>, "
      << "--te2=<te2_in_millisec> (default: 9.1, 30)\n"
//      << "--te3=<next_echo_time>, etc.\n"
      << "--tag-pattern=<string_of_Ts_and_Cs> (default: TC)\n"      
      << "--t1b=<T1_of_blood> (default: 1.66), --t1b-stdev=<stdev> (to add it as a parameter)\n"
      << "--dt=<bolus_arrival_time>, --dt-stdev (default: --dt=0.5 --dt-stdev=0.25)\n"
      << "--inv-eff=<inversion_efficiency>, --inv-eff-stdev=<stdev> (to add it as a parameter)\n\n"      
      ;
}

pcASLFwdModel::pcASLFwdModel(ArgsType& args)
{
    string scanParams = args.ReadWithDefault("scan-params","cmdline");
    string tagPattern;
    
    if (scanParams == "cmdline")
    {
        TI = convertTo<double>(args.ReadWithDefault("ti","2.00"));
        Tau = convertTo<double>(args.ReadWithDefault("tau","1.00"));
	    stdevT1b = convertTo<double>(args.ReadWithDefault("t1b-stdev", "0"));
        fixedT1b = convertTo<double>(args.ReadWithDefault("t1b","1.66"));
	    stdevInvEff = convertTo<double>(args.ReadWithDefault("inv-eff-stdev","0"));
        fixedInvEff = convertTo<double>(args.ReadWithDefault("inv-eff","1"));
        stdevDt = convertTo<double>(args.ReadWithDefault("dt-stdev","0.25"));
        fixedDt = convertTo<double>(args.ReadWithDefault("dt","0.5"));
        
        if (stdevInvEff < 0 || stdevDt < 0 || stdevT1b < 0)
            throw Invalid_option("standard deviations must not be negative!");
    
        tagPattern = args.ReadWithDefault("tag-pattern","TC");
        if (tagPattern.find_first_not_of("TCtc") != tagPattern.npos)
            throw Invalid_option("tagpattern string must contain only Ts and Cs!");
        
        echoTime.ReSize(2);
        echoTime(1) = convertTo<double>(args.ReadWithDefault("te1","9.1"))/1000.0;
        echoTime(2) = convertTo<double>(args.ReadWithDefault("te2","30"))/1000.0;

	while (true) 
	  {
	    int N = echoTime.Nrows()+1;
	    string teString = args.ReadWithDefault("te"+stringify(N), "stop!");
	    if (teString == "stop!") break;

            // This just isn't tested enough (at all)... remove if you dare
            throw Invalid_option(
              "Using more than two echo times is implemented but completely untested... modify the code if you really want to try it.");

	    // Append this TE to the list of TEs
	    ColumnVector tmp(1); 
	    tmp = atof(teString.c_str())/1000.0;
	    echoTime &= tmp; // vertcat

	    // Sanity checks:
	    if (echoTime(N) <= 0.001)
	      throw Invalid_option(
		"Was expecting TE > 1 ms (don't use seconds!)");
	    if (echoTime(N) > 0.500)
	      throw Invalid_option("Was expecting TE < 500 ms");
	  }
        
//        if (echoTime.Nrows() < 1)
//          throw Invalid_option("The --te1=<ms> option is mandatory for --model=quipss2");

        LOG << "    Scan parameters: --ti=" << TI << " --tau=" << Tau
            << " --t1b=" << fixedT1b << "--t1b-stdev=" << stdevT1b
            << " --inv-eff=" << fixedInvEff << " --inv-eff-stdev=" << stdevInvEff
            << " --dt=" << fixedDt << "--dt-stdev=" << stdevDt
            << " --tag-pattern=" << tagPattern
            ;
	for (int i = 1; i <= echoTime.Nrows(); i++)
	  LOG << " --te" << i << "=" << echoTime(i)*1000.0;
	LOG << endl;
    }
    // It should also be possible to parse most of this information straight out of 
    // the .HEAD file.  For example:
    // type = string-attribute
    // name = NOTE_NUMBER_001
    // count = 213
    // 'Sequence parameters for file run1brik_e01:\nspep 35.5, reps 130, nEcho 2, 
    // TR 2000000, TE 9100, TE2 30000, nIntlv 1, nCoil 8, nPix 64, FOV 240.0, 
    // nSlice 3, slThick 8.0, slGap 0.0, slDelay 54740, spdir1 0, spdir2 0~
    else
        throw Invalid_option("Only --scan-params=cmdline is accepted at the moment");    
    
    string rb = args.Read("bold-basis"); // only mandatory basis set
    // default value was: "/usr/people/woolrich/scratch/tldata/analysis_protocols/response_fromroi/cbvdesign.mat";
    
//    string qb = args.ReadWithDefault("cbf-basis", rb);  
//    string mb = args.ReadWithDefault("statmag-basis", qb);
    string qb = args.Read("cbf-basis");  
    string mb = args.Read("statmag-basis");
    string nb = args.ReadWithDefault("nuisance-basis", "null");
    
    LOG_ERR( "    Reading BOLD basis functions: " << rb << endl );
    if (rb != "null")
        Rbasis = read_vest(rb);
    else
        throw Invalid_option("Currently --bold-basis=null isn't allowed...");
        // Gotta get the data length from somewhere.  Haven't loaded the data yet.
    
    const int numTR = Rbasis.Nrows();
        
    LOG_ERR( "    Reading CBF basis functions: " << qb << endl );
    if (qb != "null")
        Qbasis = read_vest(qb);
    else
        Qbasis.ReSize(numTR, 0);
        
    LOG_ERR( "    Reading Stat. Mag. basis functions: " << mb << endl );
    if (mb != "null")
        Mbasis = read_vest(mb);
    else
        Mbasis.ReSize(numTR, 0);
        
    LOG_ERR( "    Reading Nuisance basis functions: " << nb << endl );
    if (nb != "null")
        Nbasis = read_vest(nb);
    else
        Nbasis.ReSize(numTR, 0);
    
    // Now we know how many basis functions -> define parameter vector
//    id.Define(Qbasis.Ncols(), Mbasis.Ncols(), Rbasis.Ncols(), 
//              Nbasis.Ncols(), echoTime.Nrows());
    
    // Now we can parse the TagPattern string.
        
    rho.ReSize(numTR);
    
    for (unsigned i = 1; i <= tagPattern.length(); i++)
        rho(i) = (toupper(tagPattern[i-1]) == 'T') ? -1 : 1;
        
    for (int i = tagPattern.length()+1; i<=numTR; i++)
        rho(i) = rho(i-tagPattern.length());

    LOG << "Full tag-control pattern used (" << rho.Nrows() << " TRs): ";
    for (int i = 1; i <= rho.Nrows(); i++)
	LOG << (rho(i)>0 ? "C" : "T");
    LOG << endl;
}

void pcASLFwdModel::DumpParameters(const ColumnVector& vec,
                                    const string& indent) const
{
    LOG << indent << "Baseline parameters:" << endl;
    LOG << indent << "  Q0 == " << vec(Q0index()) << " (baseline CBF)\n";
    LOG << indent << "  M0 == " << vec(M0index()) << " (baseline Stat. Mag.)\n";
    LOG << indent << "  R0 == " << vec(R0index()) << " (baseline T2*)\n";
    if (stdevT1b>0)
      LOG << indent << "  T1b == " << vec(T1bIndex()) << " (T1 of blood)\n";
    if (stdevInvEff>0)
      LOG << indent << "  inv-eff == " << vec(InvEffIndex()) << " (inversion efficiency)\n";
    if (stdevDt>0)
      LOG << indent << "  dt == " << vec(dtIndex()) << " (constant bolus arrival time)\n"; 
    LOG << indent << "Absolute change parameters (CBF, StatMag, BOLD effect):" << endl;
    LOG << indent << "  Qn == " << QnOf(vec).t();// << "]\n";
    LOG << indent << "  Mn == " << MnOf(vec).t();// << "]\n";
    LOG << indent << "  Rn == " << RnOf(vec).t();// << "]\n";
    LOG << indent << "Nuisance regressors (one line per TE):\n";
    for (int i = 1; i <= echoTime.Nrows(); i++)
        LOG << indent << "  Nn == " << NnOf(i,vec).t();
}

void pcASLFwdModel::NameParams(vector<string>& names) const
{
    names.clear();
    
    names.push_back("Q0");
    for (int i = 1; i <= Qbasis.Ncols(); i++)
      names.push_back("Q_abschg_" + stringify(i));
    
    names.push_back("M0");
    for (int i = 1; i <= Mbasis.Ncols(); i++)
      names.push_back("M_abschg_" + stringify(i));
    
    names.push_back("R0");
    for (int i = 1; i <= Rbasis.Ncols(); i++)
      names.push_back("BOLD_abschg_" + stringify(i));

    for (int i = 1; i <= Nbasis.Ncols(); i++)
      names.push_back("Nuisance_signal_" + stringify(i));

    if (stdevInvEff>0)
      names.push_back("InvEff");

    if (stdevT1b>0)
      names.push_back("T1b");

    if (stdevDt>0)
      names.push_back("dt");
  
    assert(names.size() == unsigned(NumParams())); 
}
