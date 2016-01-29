/*  fwdmodel_q2tips.cc - Implements a Q2TIPS dual-echo ASL model

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

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

#include "fwdmodel_q2tips.h"

#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include "newimage/newimageall.h"
using namespace NEWIMAGE;
#include "easylog.h"



string Q2tipsFwdModel::ModelVersion() const
{
  return "$Id: fwdmodel_q2tips.cc,v 1.2 2007/08/02 15:14:11 adriang Exp $ and "
	+ Quipss2FwdModel::ModelVersion();
}


void Q2tipsFwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
    Tracer_Plus tr("Q2tipsFwdModel::Evaluate");
    // Adapted from original_fwdmodel.m
    
    // Parameterization used in most recent results:
    // Absolute M and Q change (same units as M0 or Q0):
    ColumnVector StatMag = params(M0index()) - Mbasis * MnOf(params);
    ColumnVector CBF = params(Q0index()) + Qbasis * QnOf(params);
    // Fractional change in BOLD effect (at TE_2), rather than using % R2* change 
    ColumnVector R2s = -1/echoTime(2) * log( 
        Rbasis * RnOf(params) + exp(-echoTime(2)*params(R0index())));

    // The following are relative magnetizations    
    double pretag = 1; // untagged
    double T1b = (stdevT1b>0 ? params(T1bIndex()) : fixedT1b);
    double invEfficiency = (stdevInvEff>0 ? params(InvEffIndex()) : fixedInvEff);
    ColumnVector bolus = 1 - (1-rho)*invEfficiency*exp(-TI2/T1b); // tag or control
    //    double posttag = 1 - exp(-(TI2-TI1)/T1b); // saturated
    double dt = (stdevDt>0? params(dtIndex()) : fixedDt);
    //    ColumnVector Sb = SP( CBF,  // SP(a,b) means a.*b 
    //        pretag*dt + bolus*TI1 + posttag*(TI2-TI1-dt) );
    double posttagQ2 = ((TI2-TI1-dt) 
			+ T1b*exp(-(TI2-TI1)/T1b)
			- T1b*exp(-dt/T1b));
    ColumnVector Sb = SP( CBF,  // SP(a,b) means a.*b 
        pretag*dt + bolus*TI1 + posttagQ2 );
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
