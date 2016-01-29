/*  fwdmodel_simple.cc - Implements the simplified ASL model

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

#include "fwdmodel_simple.h"

#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include "newimage/newimageall.h"
using namespace NEWIMAGE;
#include "easylog.h"

string SimpleFwdModel::ModelVersion() const
{
  return "$Id: fwdmodel_simple.cc,v 1.20 2015/09/08 13:55:19 mwebster Exp $";
}

void SimpleFwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
//  ColumnVector Qn = id.Qn(params);
  ColumnVector tcModulated = id.Q0(params) 
	* ( 1 + basis * id.Qn(params) / 100.0);
  ColumnVector tcUnmodulated = id.U0(params) 
	* ( 1 + basis * id.Un(params) / 100.0);
  ColumnVector R2s = id.R0(params) * (1 + basis * id.Rn(params) / 100.0);
  ColumnVector S = SP(rho, tcModulated) + tcUnmodulated; // SP means .*
  
  int Ntimes = R2s.Nrows();
  if (result.Nrows() != 2*Ntimes)
    result.ReSize(2*Ntimes);
    
//  result = 0.0/0.0; // pre-fill with nans to check all overwritten
  
  for (int te = 1; te <= echoTime.Nrows(); te++)
    {
      for (int i = 1; i <= Ntimes; i++)
        result( 2*i - 2 + te ) = S(i) * exp(-echoTime(te) * R2s(i));
      // Fill order: te1 te2 te1 te2 te1 te2 te1 te2 ...
    }

/*    
  LOG << "Fwdmodel input:\n" << params;
  LOG << "tcModulated:" << endl;
  LOG << tcModulated;
  LOG << "tcUnmodulated:" << endl;
  LOG << tcUnmodulated;
  LOG << "R2s:" << endl;
  LOG << R2s;
  LOG << "S" << endl << S;
  LOG << "echoTime" << echoTime;
  LOG << "Output:\n" << result;
*/   
//    LOG << "BASIS\n" << basis.t()*basis;

  return;
}

SimpleFwdModel::SimpleFwdModel(ArgsType& args)
{
    
  string scanParams = args.ReadWithDefault("scan-params","hardcoded");
  if (scanParams == "hardcoded")
    {
      LOG << "  Loading hardcoded forward-model priors" << endl;
      echoTime.ReSize(2);
      echoTime << 9.1 << 30;
      echoTime *= 0.001;

      //      LOG << "SimpleFwdModel::SimpleFwdModel isn't implemented yet!" << endl;

      // basis.. not even going to bother trying.  Load from a file.
//      string basisFile = "/home/fs11/adriang/proj/response_fir/bolddesign.mat";
      string basisFile = "/usr/people/woolrich/scratch/tldata/analysis_protocols/response_fromroi/cbvdesign.mat";
      LOG << "    Reading basis functions from file: " << basisFile << endl;
      basis = read_vest(basisFile);
      //      LOG << "basis == \n" << basis << endl;
      // Nrows = 136, Ncols = 15
      // LOG << "BASIS: " << basis.Nrows() << basis.Ncols() << endl;

      LOG << "      Read " << basis.Ncols() << " basis functions" << endl;
      id.Define(basis.Ncols());

      rho.ReSize(68*2);
      rho(1) = -1;
      for (int i = 2; i <= 68*2; i++)
        rho(i) = rho(i-1) * -1;

      //      LOG << "echoTime:\n" << echoTime << endl;
      //      LOG << "rho:\n" << rho << endl;

    } 
  else
    throw Invalid_option("Only --scan-params=hardcoded is accepted at the moment");
}

void SimpleFwdModel::DumpParameters(const ColumnVector& vec,
                                    const string& indent) const
{
    LOG << indent << "Baseline parameters:" << endl;
    LOG << indent << "  U0 == " << id.U0(vec) << " (baseline unmodulated mag.)\n";
    LOG << indent << "  Q0 == " << id.Q0(vec) << " (baseline modulated mag.)\n";
    LOG << indent << "  R0 == " << id.R0(vec) << " (baseline T2*)\n";
    LOG << indent << "Percent change parameters:" << endl;
    LOG << indent << "  Un == " << id.Un(vec).t();// << "]\n";
    LOG << indent << "  Qn == " << id.Qn(vec).t();// << "]\n";
    LOG << indent << "  Rn == " << id.Rn(vec).t();// << "]\n";
}    

void SimpleFwdModelIdStruct::NameParams(vector<string>& names) const
{
    names.clear();
    
    for (int p = 1; p <= 3; p++)
    {
        string letter;
        switch(p) { 
            case 1: letter="Q"; break;
            case 2: letter="U"; break;
            case 3: letter="R"; break;
        } 
        names.push_back(letter + "0");
        for (int i = 1; i <= Nbasis; i++)
        {
            names.push_back(letter + "_percchg_" + stringify(i));
        }
    }
    assert(names.size() == unsigned(3+3*Nbasis)); 
}
