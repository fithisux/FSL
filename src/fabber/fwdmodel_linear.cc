/*  fwdmodel_linear.cc - Linear forward model and related classes

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007-2008 University of Oxford  */

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

#include "fwdmodel_linear.h"
#include <iostream>
#include "newmatio.h"
#include <stdexcept>
#include "easylog.h"
#include "newimage/newimageall.h"
using namespace NEWIMAGE;

string LinearFwdModel::ModelVersion() const
{
  return "$Id: fwdmodel_linear.cc,v 1.20 2015/09/08 13:55:19 mwebster Exp $";
}

void LinearFwdModel::ModelUsage()
{
  cout << "\nUsage info for --model=linear:\n"
       << "Required options:\n"
       << "--basis=<design_file>\n"
    ;
}

LinearFwdModel::LinearFwdModel(ArgsType& args)
{
  Tracer_Plus tr("LinearFwdModel::LinearFwdModel(args)");
  string designFile = args.Read("basis");
  LOG_ERR("    Reading design file: " << designFile << endl);
  jacobian = read_vest(designFile);

  const int Ntimes = jacobian.Nrows();
  const int Nbasis = jacobian.Ncols();

  LOG_ERR("      Loaded " << jacobian.Ncols() 
	  << " basis functions of length " << Ntimes << endl);

  centre.ReSize(Nbasis); 
  centre = 0;
  offset.ReSize(Ntimes);
  offset = 0;

  if (args.ReadBool("add-ones-regressor"))
    {
      LOG_ERR("      Plus an additional regressor of all ones\n");
      ColumnVector ones(Ntimes); ones = 1.0;
      jacobian = jacobian | ones;
      centre.ReSize(Nbasis+1);
    }
  // Warning: Nbasis is now wrong!
}

void LinearFwdModel::HardcodedInitialDists(MVNDist& prior, 
					   MVNDist& posterior) const
{
  Tracer_Plus tr("LinearFwdModel::HardcodedInitialDists");
  assert(prior.means.Nrows() == NumParams());

  prior.means = 0;
  prior.SetPrecisions(IdentityMatrix(NumParams()) * 1e-12);
  posterior = prior;

}

void LinearFwdModel::Evaluate(const ColumnVector& params,
	               		    ColumnVector& result) const
{
  result = jacobian * (params - centre) + offset;
}

void LinearizedFwdModel::ReCentre(const ColumnVector& about)
{
  Tracer_Plus tr("LinearizedFwdModel::ReCentre");
  assert(about == about); // isfinite

  // Store new centre & offset
  centre = about;
  fcn->Evaluate(centre, offset);
  if (0*offset != 0*offset) 
    {
      LOG_ERR("about:\n" << about);
      LOG_ERR("offset:\n" << offset.t());
      throw overflow_error("ReCentre: Non-finite values found in offset");
    }

  // Calculate the Jacobian numerically.  jacobian is len(y)-by-len(m)
  jacobian.ReSize(offset.Nrows(), centre.Nrows());
  // jacobian = 0.0/0.0; // fill with NaNs to check
  
  // try and get the gradient from the model first
  int gradfrommodel=false;
  gradfrommodel = fcn->Gradient(centre,jacobian);

  if (!gradfrommodel) {
  ColumnVector centre2, centre3;
  ColumnVector offset2, offset3;
  for (int i = 1; i <= centre.Nrows(); i++)
    {
      double delta = centre(i) * 1e-5;
      if (delta<0) delta = -delta;
      if (delta<1e-10) delta = 1e-10;

      // Take derivative numerically
      centre3 = centre;
      centre2 = centre;
      centre2(i) += delta;
      centre3(i) -= delta;
      fcn->Evaluate(centre2, offset2);
      fcn->Evaluate(centre3, offset3);
      jacobian.Column(i) = (offset2 - offset3) / (centre2(i) - centre3(i));

      /*
if (i==4)
{LOG << "centre2 -centre3== \n" << 1e10*(centre2-centre3) << endl;
LOG << "offset2-offset3 == \n" << offset2(33)-offset3(33) << endl;
LOG << "offset2-offset3 == \n" << float(offset2(33)-offset3(33)) << endl;
LOG << "offset2-offset3 == \n" << double(offset2(33)-offset3(33)) << endl;
LOG << "Jac 33,4 == " << jacobian(33,4) << endl;
}
//*/
    }   
  }

  if (0*jacobian != 0*jacobian) 
    {
      LOG << "jacobian:\n" << jacobian;
      LOG << "about':\n" << about.t();
      LOG << "offset':\n" << offset.t();    
      throw overflow_error("ReCentre: Non-finite values found in jacobian");
    }
}

void LinearFwdModel::DumpParameters(const ColumnVector& vec,
                                    const string& indent) const
{
    LOG << indent << "Parameters mean nothing to me!  "
                   << "I am a mere linear model." << endl;
    LOG << indent << "  Vector: " << vec.t() << endl;
}

void LinearFwdModel::NameParams(vector<string>& names) const
{
    names.clear();
 
    // Completely generic names.
    for (int i = 1; i <= NumParams(); i++)
        names.push_back("Parameter_" + stringify(i));
}

void LinearizedFwdModel::DumpParameters(const ColumnVector& vec,
                                    const string& indent) const
{
//    LOG << indent << "This is what the nonlinear model has to say:" << endl;
    fcn->DumpParameters(vec, indent);
}
