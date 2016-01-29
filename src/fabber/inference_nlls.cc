/* inference_nlls.h - Non-Linear Least Squares class declarations

   Adrian Groves Michael Chappell, FMRIB Image Analysis Group

   Copyright (C) 2007 University of Oxford */
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
#include "inference_nlls.h"

void NLLSInferenceTechnique::Setup(ArgsType& args)
{
  Tracer_Plus tr("NLLSInferenceTechnique::Setup");
  model = FwdModel::NewFromName(args.Read("model"), args);
  assert( model->NumParams() > 0 );
  LOG_ERR("    Forward Model version:\n      " 
	  << model->ModelVersion() << endl);

  //determine whether NLLS is being run in isolation or as a pre-step for VB (alters what we do if result is ill conditioned)
  vbinit = args.ReadBool("vb-init");

  // option to load a 'posterior' which will allow the setting of intial parameter estmates for NLLS
  MVNDist* loadPosterior = new MVNDist( model->NumParams() );
  MVNDist* junk = new MVNDist( model->NumParams() );
  string filePosterior = args.ReadWithDefault("fwd-inital-posterior","modeldefault");
  model->HardcodedInitialDists(*junk, *loadPosterior);
  if (filePosterior != "modeldefault") loadPosterior->Load(filePosterior);

  //assert(initialFwdPosterior == NULL);
  initialFwdPosterior = loadPosterior;
  loadPosterior = NULL;

  lm = args.ReadBool("lm"); //determine whether we use L (default) or LM converengce

}

void NLLSInferenceTechnique::DoCalculations(const DataSet& allData)
{
  Tracer_Plus tr("NLLSInferenceTechnique::DoCalculations");
  //get data for this voxel
  const Matrix& data = allData.GetVoxelData();
  const Matrix & coords = allData.GetVoxelCoords();
  unsigned int Nvoxels = data.Ncols();
  int Nsamples = data.Nrows();
  if (data.Nrows() != model->NumOutputs())
    throw Invalid_option("Data length (" 
      + stringify(data.Nrows())
      + ") does not match model's output length ("
      + stringify(model->NumOutputs())
      + ")!");

  for (unsigned int voxel = 1; voxel <= Nvoxels; voxel++)
    {
      ColumnVector y = data.Column(voxel);
      ColumnVector vcoords = coords.Column(voxel);
      // some models might want more information about the data
      model->pass_in_data( y );
      model->pass_in_coords(vcoords);
      
      LOG_ERR("  Voxel " << voxel << " of " << Nvoxels << endl);

      MVNDist fwdPosterior;
      LinearizedFwdModel linear(model);

      int Nparams = initialFwdPosterior->GetSize();
      fwdPosterior.SetSize(Nparams);
      IdentityMatrix I(Nparams);

      NLLSCF costfn(y, model);
      NonlinParam nlinpar(Nparams,NL_LM);

      if (!lm)
	{ nlinpar.SetGaussNewtonType(LM_L); }

      // set ics from 'posterior'
      ColumnVector nlinics = initialFwdPosterior->means;
      nlinpar.SetStartingEstimate(nlinics);
      nlinpar.LogPar( true);nlinpar.LogCF(true);
      
 
      try {
	 __attribute__((unused)) NonlinOut status = nonlin(nlinpar,costfn);
	 // Status is unused - unsure if nonlin has any effect so telling compiler to ignore the status variable

	/*cout << "The solution is: " << nlinpar.Par() << endl;
	cout << "and this is the process " << endl;
	for (int i=0; i<nlinpar.CFHistory().size(); i++) {
	  cout << " cf: " << (nlinpar.CFHistory())[i] <<endl;
	}
	for (int i=0; i<nlinpar.ParHistory().size(); i++) {
	  cout << (nlinpar.ParHistory())[i] << ": :";
	  }*/

	fwdPosterior.means = nlinpar.Par();

	// recenter linearized model on new parameters
	linear.ReCentre( fwdPosterior.means );
	const Matrix& J = linear.Jacobian();
	// Calculate the NLLS covariance
	/* this is inv(J'*J)*mse?*/
	double sqerr = costfn.cf( fwdPosterior.means );
	double mse = sqerr/(Nsamples - Nparams);
	
	/*	Matrix Q = J;
	UpperTriangularMatrix R;
	QRZ(Q,R);
	Matrix Rinv = R.i();
	SymmetricMatrix nllscov;
	nllscov = Rinv.t()*Rinv*mse;
	
	fwdPosterior.SetCovariance( nllscov );*/

	
	SymmetricMatrix nllsprec;
      	nllsprec << J.t()*J/mse;
	
	// look for zero diagonal elements (implies parameter is not observable) 
	//and set precision small, but non-zero - so that covariance can be calculated
	for (int i=1; i<=nllsprec.Nrows(); i++)
	  {
	    if (nllsprec(i,i) < 1e-6)
	      {
		nllsprec(i,i) = 1e-6;
	      }
	  }
	fwdPosterior.SetPrecisions( nllsprec );
	fwdPosterior.GetCovariance();

      }

catch (Exception)
	{
	  LOG_ERR("   NEWMAT Exception in this voxel:\n"
		  << Exception::what() << endl);
	  
	  //if (haltOnBadVoxel) throw;
    
	  LOG_ERR("   Estimates in this voxel may be unreliable" <<endl
		  << "(precision matrix will be set manually)" <<endl
		  << "   Going on to the next voxel" << endl);

	    // output the results where we are
	    fwdPosterior.means = nlinpar.Par();

	    // recenter linearized model on new parameters
	    linear.ReCentre( fwdPosterior.means );
	    
	    // precision matrix is probably singular so set manually
	    fwdPosterior.SetPrecisions(  I*1e-12 );
	 
	}

      resultMVNs.push_back(new MVNDist(fwdPosterior));
      assert(resultMVNs.size() == voxel);
    }
}

NLLSInferenceTechnique::~NLLSInferenceTechnique()
{

}

double NLLSCF::cf(const ColumnVector& p) const
{
  Tracer_Plus tr("NLLSCF::cf");
  ColumnVector yhat;
  model->Evaluate(p,yhat);

  double cfv = ( (y-yhat).t() * (y-yhat) ).AsScalar();

  /*double cfv = 0.0;
  for (int i=1; i<=y.Nrows(); i++) { //sum of squares cost function
    double err = y(i) - yhat(i);
    cfv += err*err;
    }*/
  return(cfv);
}

ReturnMatrix NLLSCF::grad(const ColumnVector& p) const
{
  Tracer_Plus tr("NLLSCF::grad");
  ColumnVector gradv(p.Nrows());
  gradv=0.0;

  // need to recenter the linearised model to the current parameter values
  linear.ReCentre( p );
  const Matrix& J = linear.Jacobian();
  //const ColumnVector gm = linear.Offset(); //this is g(w) i.e. model evaluated at current parameters?
  ColumnVector yhat;
  model->Evaluate(p,yhat);

  gradv = -2*J.t()*(y-yhat);

  gradv.Release();
  return(gradv);
  }

boost::shared_ptr<BFMatrix> NLLSCF::hess(const ColumnVector& p, boost::shared_ptr<BFMatrix> iptr) const
{
  Tracer_Plus tr("NLLSCF::hess");
  boost::shared_ptr<BFMatrix> hessm;

  if (iptr && iptr->Nrows()==(unsigned)p.Nrows() && iptr->Ncols()==(unsigned)p.Nrows())
    { hessm = iptr; }
  else
    {
      hessm = boost::shared_ptr<BFMatrix>(new FullBFMatrix(p.Nrows(),p.Nrows()));
    }
  
  // need to recenter the linearised model to the current parameter values
  linear.ReCentre( p );
  const Matrix& J = linear.Jacobian();
  Matrix hesstemp = 2*J.t()*J; //Make the G-N approximation to the hessian

  //(*hessm) = J.t()*J;

  for (int i=1; i<=p.Nrows(); i++) { for (int j=1; j<=p.Nrows(); j++) hessm->Set(i,j,hesstemp(i,j));}

  return(hessm);
  }
  
