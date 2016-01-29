/*  fwdmodel_flobs.cc - Does FLOBS

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


#include "fwdmodel_flobs.h"

#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include "newimage/newimageall.h"
using namespace NEWIMAGE;
#include "easylog.h"

string FlobsFwdModel::ModelVersion() const
{
  return "$Id: fwdmodel_flobs.cc,v 1.6 2015/09/08 13:55:19 mwebster Exp $";
}

void FlobsFwdModel::HardcodedInitialDists(MVNDist& prior, 
    MVNDist& posterior) const
{
    Tracer_Plus tr("FlobsFwdModel::HardcodedInitialDists");
    assert(prior.means.Nrows() == NumParams());
    
    // Set priors
    prior.means = 0;
    SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;

    prior.SetPrecisions(precisions);
    
    posterior = prior;
    //    if (useSeparateScale)
    //      {
    //	// Set informative initial posterior
    //	// Shape = first basis function, magnitude 1.
    //	posterior.means(1) = 1;
    //	posterior.means(NumParams()) = 1;
    //      }
}

void FlobsFwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
  Tracer_Plus tr("FlobsFwdModel::Evaluate");

  assert(params.Nrows() == NumParams());
  
  //  if (useSeparateScale)
  //    {
  //      ColumnVector beta = params.Rows(1,basis.Ncols());
  //      //double betaBar = params(NumParams());
  //      result = basis * beta; // * betaBar; Now it's purely linear and ignore betaBar!
  //      assert(false); // TODO: update for nuisance parameters
  //    }

  if (usePolarCoordinates)
    {
      assert(basis.Ncols() == 2);
      ColumnVector beta = params.Rows(1,basis.Ncols());
      double betaBar = params(1);
      beta(1) = cos(beta(2));
      beta(2) = sin(beta(2));
      result = basis * beta * betaBar;

      // No nuisance stuff yet.
      assert(params.Nrows() == 2);
    }
  else
    { 
      ColumnVector beta = params.Rows(1,basis.Ncols());
      double betaBar = params(1);
      beta(1) = 1.0; // fixed shape=1, scale=betaBar
      result = basis * beta * betaBar;
      
      beta = params.Rows(basis.Ncols()+1,params.Nrows());
      result += nuisanceBasis * beta;
    }

  return; // answer is in the "result" vector
}

void FlobsFwdModel::ModelUsage()
{
  //  if (useSeparateScale)
    {
      cout << "Usage for --model=flobs6:\n"
	   << "  --basis=<basis_functions>\n"
	   << "  A scale parameter will automatically be added.\n"
	   << "  Always use with the --flobs-prior-adjust option!!! \n"
	   << "  (Future work: specify several scaling factors, for multi-event stimuli)\n\n";
    } 
    //  else 
    {
      cout << "Usage for --model=flobs5:\n"
	   << "  --basis=<basis_functions>\n"
	   << "  The first basis function will serve as the scaling factor (fixed shape==1)\n";
    }
    // TODO: add --nuisance= option
}

FlobsFwdModel::FlobsFwdModel(ArgsType& args, bool polar) 
//  : useSeparateScale(sepScale)
  : usePolarCoordinates(polar)
{
    string basisFile = args.Read("basis"); 
    LOG_ERR( "    Reading basis functions: " << basisFile << endl );
    basis = read_vest(basisFile);
    LOG_ERR( "    Read " << basis.Ncols() << " basis functions of length " 
	     << basis.Nrows() << endl);

    basisFile = args.ReadWithDefault("nuisance","null"); 
    if (basisFile == "null")
      {
	nuisanceBasis.ReSize(basis.Nrows(), 0);
      }
    else if (basisFile == "offset")
      {
	nuisanceBasis.ReSize(basis.Nrows(), 1);
	nuisanceBasis = 1;
      }
    else
      {
	LOG_ERR( "    Reading nuisance basis functions: " 
		 << basisFile << endl );
	nuisanceBasis = read_vest(basisFile);
	LOG_ERR( "    Read " << nuisanceBasis.Ncols() 
		 << " nuisance basis functions of length "
		 << nuisanceBasis.Nrows() << endl);
      }

    if (nuisanceBasis.Nrows() != basis.Nrows())
      throw Invalid_option("Basis length mismatch!\n");
}

void FlobsFwdModel::DumpParameters(const ColumnVector& vec,
                                    const string& indent) const
{
  //  if (useSeparateScale)
  //    LOG << indent << "Scale = " << vec(NumParams()) << ", shape:\n"
  //	<< vec.Rows(1, NumParams()-1);
  //  else
    LOG << indent << "Scale = " << vec(1) << ", shape:\n1 (fixed)\n"
	<< vec.Rows(2, NumParams());
  // TODO: should dump nuisanceBasis too
}

void FlobsFwdModel::NameParams(vector<string>& names) const
{
    names.clear();
    
    for (int i = 1; i <= basis.Ncols(); i++)
      names.push_back("basis_" + stringify(i));
    
    //    if (useSeparateScale)
    //      names.push_back("scale");

    for (int i = 1; i <= nuisanceBasis.Ncols(); i++)
      names.push_back("nuisance_"+stringify(i));
  
    assert(names.size() == (unsigned)NumParams()); 
}
