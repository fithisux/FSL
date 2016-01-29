/*  fwdmodel_custom.cc - A place for your own quick'n'dirty forward model implementations.

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2008 University of Oxford  */

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


/* fwdmodel_custom.h
 * Implementation for a custom forward model class
 * Template written by Adrian Groves, 2008
 * FMRIB Centre, University of Oxford
 *
 * Last modified: $Date: 2015/09/08 13:55:19 $ $Author: mwebster $ $Revision: 1.3 $
 */

#include "fwdmodel_custom.h"

// Constructor:
CustomFwdModel::CustomFwdModel(ArgsType& args)
{
  // Read any command-line arguments.  These should probably be stored in 

  // e.g. Mandatory option:
  // numRepeats = convertTo<int>(args.Read("reps"));
  // An error will be thrown if there's no --reps=NNN option.

  // e.g. Optional option:
  // T1 = convertTo<double>(args.ReadWithDefault("T1","1.6"))

  // e.g. a basis function:
  // string basisFile = args.Read("basis"); // error if --basis=??? option is missing
  // basis = read_vest(basisFile); // error if file is missing

  // e.g. a boolean option (present or absent);
  // wantFriesWithThat = args.ReadBool("chips");



  
  // TODO: Read any commmand-line arguments and save them in the class's member variables.





  // Before you return, all of the implementation-specific variables in fwdmodel_custom.h 
  // should have had values assigned to them.  If you didn't add any variables then you 
  // don't have to do anything here!
  return;
}

void CustomFwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
  assert(params.Nrows() == NumParams());

  // TODO:This needs to equal the number of timepoints in your data.
  // If it varies, you need to be able to calculate it from the command-line options.
  // This may mean adding a (redundant) --data-length=NNN option above.
  const int dataLength = 10; // Example value

  if (result.Nrows() != dataLength)
    result.ReSize(dataLength);

  // TODO: Your model here!

  // A very simple forward model: fitting a quadratic equation, with inputs at 1..10
  // Notice that NEWMAT vectors and matrices could from one, not zero!
  const double constantTerm = params(1);
  const double linearTerm = params(2);
  const double quadraticTerm = params(3);
  
  for (int i = 1; i <= dataLength; i++)
    {
      result(i) = constantTerm + i*linearTerm + i*i*quadraticTerm;
    }

  // Before you return, you should have assigned a predicted signal to "result".
  return;
}

int CustomFwdModel::NumParams() const
{
  // TODO: How many parameters does your model need?

  return 3; // Example function
}

void CustomFwdModel::NameParams(vector<string>& names) const
{
  // TODO: Name your model parameters here.  Should match NumParams() above!

  // Example function:
  names.push_back("Constant term");
  names.push_back("Linear term");
  names.push_back("Quadratic term");
}


string CustomFwdModel::ModelVersion() const
{ 
  return "$Id: fwdmodel_custom.cc,v 1.3 2015/09/08 13:55:19 mwebster Exp $";
}

void CustomFwdModel::HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const
{
  Tracer_Plus tr("CustomFwdModel::HardcodedInitialDists");
  // Pick a safe starting point for your model fits, if no other input is provided.
  // The default one just uses a N(0,1e12) prior, and starts with all parameters set to zeroes:

  assert(prior.means.Nrows() == NumParams());

  prior.means = 0;
  prior.SetPrecisions(IdentityMatrix(NumParams()) * 1e-12);
  posterior = prior;

  // Ask Adrian for help if you want to modify this!
  // You can override prior at the command line with the --fwd-initial-prior=VESTFILE option,
  // and you can set initialization points using the --continue-from-mvn=MVNFILE option.
}
 
void CustomFwdModel::ModelUsage()
{
  cout << "No model usage info available for --model=custom, yet." << endl;
}
