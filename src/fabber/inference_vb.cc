/*  inference_vb.cc - VB inference technique class declarations

    Adrian Groves and Michael Chappell, FMRIB Image Analysis Group

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

#include "inference_vb.h"
#include "convergence.h"

#ifndef __FABBER_LIBRARYONLY
using namespace NEWIMAGE;
#endif

void VariationalBayesInferenceTechnique::Setup(ArgsType& args) 
{ 
  Tracer_Plus tr("VariationalBayesInferenceTechnique::Setup");

  // Call ancestor, which does most of the real work
  InferenceTechnique::Setup(args);

  // Load up initial prior and initial posterior
  MVNDist* loadPrior = new MVNDist( model->NumParams() );
  MVNDist* loadPosterior = new MVNDist( model->NumParams() );
  NoiseParams* loadNoisePrior = noise->NewParams();
  NoiseParams* loadNoisePosterior = noise->NewParams();
  
  string filePrior = args.ReadWithDefault("fwd-initial-prior", "modeldefault");
  string filePosterior = args.ReadWithDefault("fwd-initial-posterior", "modeldefault");
  if (filePrior == "modeldefault" || filePosterior == "modeldefault")
      model->HardcodedInitialDists(*loadPrior, *loadPosterior);
  if (filePrior != "modeldefault") loadPrior->Load(filePrior);
  if (filePosterior != "modeldefault") loadPosterior->Load(filePosterior);

  if ( (loadPosterior->GetSize() != model->NumParams()) 
    || (loadPrior->GetSize() != model->NumParams()) )
      throw Invalid_option("Size mismatch: model wants " 
      	+ stringify(model->NumParams())
	+ ", initial prior (" + filePrior + ") is " 
	+ stringify(loadPrior->GetSize()) 
	+ ", initial posterior (" + filePosterior + ") is "
        + stringify(loadPosterior->GetSize())
	+ "\n");

  filePrior = args.ReadWithDefault("noise-initial-prior", "modeldefault");
  filePosterior = args.ReadWithDefault("noise-initial-posterior", "modeldefault");
  if (filePrior == "modeldefault" || filePosterior == "modeldefault")
      noise->HardcodedInitialDists(*loadNoisePrior, *loadNoisePosterior);
  if (filePrior != "modeldefault") 
      loadNoisePrior->InputFromMVN( MVNDist(filePrior) );
  if (filePosterior != "modeldefault") 
      loadNoisePosterior->InputFromMVN( MVNDist(filePosterior) );

  // Make these distributions constant:
  assert(initialFwdPrior == NULL);
  assert(initialFwdPosterior == NULL);    
  assert(initialNoisePrior == NULL);
  assert(initialNoisePosterior == NULL);    
  initialFwdPrior = loadPrior;
  initialFwdPosterior = loadPosterior;
  initialNoisePrior = loadNoisePrior;
  initialNoisePosterior = loadNoisePosterior;
  loadPrior = loadPosterior = NULL;
  loadNoisePrior = loadNoisePosterior = NULL; // now, only accessible as consts.

  // Resume from a previous run? 
  continueFromFile = args.ReadWithDefault("continue-from-mvn", "");
  paramFilename = args.ReadWithDefault("continue-from-params",""); // optional list of parameters in MVN
  if (continueFromFile != "")
  {
    // Won't need these any more.  They don't hurt, but why leave them around?
    // They can only cause trouble (if used by mistake).
    delete initialFwdPosterior; initialFwdPosterior = NULL;
    
    if ( !args.ReadBool("continue-fwd-only") )
    {
        delete initialNoisePosterior; 
        initialNoisePosterior = NULL;
    }
  }

  // Specify prior types and sources
// Section in which spatial priors can be specified by parameter name on command line
// param_spatial_priors_bymane (PSP_byname)
// this overwrites any existing entries in the string
 vector<string> modnames; //names of model parameters
 imagepriorstr.resize(model->NumParams()); 
 model->NameParams(modnames);
 int npspnames=0;
 while (true)
   {
     npspnames++;
     string bytypeidx = args.ReadWithDefault("PSP_byname"+stringify(npspnames),"stop!");
     if (bytypeidx == "stop!") break; //no more spriors have been specified

     // deal with the specification of this sprior (by name)
     //compare name to those in list of model names
     bool found=false;
     for (int p=0; p<model->NumParams(); p++) {
       if (bytypeidx == modnames[p]) {
	 found = true;
	 char pspstype = convertTo<char>(args.Read("PSP_byname"+stringify(npspnames)+"_type"));
	 PriorsTypes[p]=pspstype;
	 PSPidx.push_back(p); //record the index at which a PSP has been defined for use in spatialvb setup
	 LOG_ERR("PSP_byname parameter " << bytypeidx << " at entry " << p << ", type: " << pspstype << endl);

	 // now read in file name for an image prior (if appropriate)
	 if (pspstype == 'I') {
	   imagepriorstr[p] = args.Read("PSP_byname"+stringify(npspnames)+"_image");
	 }
       }
     }
     if (!found) {
       throw Invalid_option("ERROR: Prior specification by name, parameter " + bytypeidx + " does not exist in the model\n");
     }
   }



  // Fix the linearization centres?
  lockedLinearFile = args.ReadWithDefault("locked-linear-from-mvn","");

  // Maximum iterations allowed:
  string maxIterations = args.ReadWithDefault("max-iterations","10");
  if (maxIterations.find_first_not_of("0123456789") != string::npos)
    throw Invalid_option("--convergence=its=?? parameter must be a positive number");    
  int its = atol(maxIterations.c_str());
  if (its<=0)
    throw Invalid_option("--convergence=its=?? paramter must be positive");
  
  // Figure out convergence-testing method:
  string convergence = args.ReadWithDefault("convergence", "maxits");
  if (convergence == "maxits")
    conv = new CountingConvergenceDetector(its);
  else if (convergence == "pointzeroone")
    conv = new FchangeConvergenceDetector(its, 0.01);
  else if (convergence == "freduce")
    conv = new FreduceConvergenceDetector(its, 0.01);
  else if (convergence == "trialmode") {
    int maxtrials = convertTo<int>(args.ReadWithDefault("max-trials","10"));
    conv = new TrialModeConvergenceDetector(its, maxtrials, 0.01);
  }
  else if (convergence == "lm")
    conv = new LMConvergenceDetector(its,0.01);
  else
    throw Invalid_option("Unrecognized convergence detector: '" 
                           + convergence + "'");

  // Figure out if F needs to be calculated every iteration
  printF = args.ReadBool("print-free-energy");
  needF = conv->UseF() || printF;

  haltOnBadVoxel = !args.ReadBool("allow-bad-voxels");
  if (haltOnBadVoxel)
    LOG << "Note: numerical errors in voxels will cause the program to halt.\n"
        << "Use --allow-bad-voxels (with caution!) to keep on calculating.\n";
  else
    LOG << "Using --allow-bad-voxels: numerical errors in a voxel will\n"
	<< "simply stop the calculation of that voxel.\n"
	<< "Check log for 'Going on to the next voxel' messages.\n"
	<< "Note that you should get very few (if any) exceptions like this;"
	<< "they are probably due to bugs or a numerically unstable model.";
  
}

void VariationalBayesInferenceTechnique::DoCalculations(const DataSet& allData) 
{
  Tracer_Plus tr("VariationalBayesInferenceTechnique::DoCalculations");
  
  // extract data (and the coords) from allData for the (first) VB run
  const Matrix& origdata = allData.GetVoxelData();
  //cerr << "Data MaxAbsValue = " << MaximumAbsoluteValue(data) << endl;
  const Matrix & coords = allData.GetVoxelCoords();
  const Matrix & suppdata = allData.GetVoxelSuppData();
  // Rows are volumes
  // Columns are (time) series
  // num Rows is size of (time) series
  // num Cols is size of volumes

  // pass in some (dummy) data/coords here just in case the model relies upon it
  // use the first voxel values as our dummies
  if (suppdata.Ncols() > 0) {
    model->pass_in_data( origdata.Column(1) , suppdata.Column(1) );
  }
  else {
    model->pass_in_data( origdata.Column(1) );
  }
  model->pass_in_coords(coords.Column(1));

       
  int Nvoxels = origdata.Ncols();
  if (origdata.Nrows() != model->NumOutputs())
    throw Invalid_option("Data length (" 
      + stringify(origdata.Nrows())
      + ") does not match model's output length ("
      + stringify(model->NumOutputs())
      + ")!");

#ifdef __FABBER_MOTION
  MCobj mcobj(allData);
#endif //__FABBER_MOTION

  Matrix data(origdata.Nrows(),Nvoxels);
  data = origdata;
  Matrix modelpred(model->NumOutputs(),Nvoxels); //use this to store the model predictions in to pass to motion correction routine

  assert(resultMVNs.empty()); // Only call DoCalculations once
  resultMVNs.resize(Nvoxels, NULL);

  assert(resultFs.empty());
  resultFs.resize(Nvoxels, 9999);  // 9999 is a garbage default value

  // If we're continuing from previous saved results, load them here:
  bool continuingFromFile = (continueFromFile != "");
  vector<MVNDist*> continueFromDists;
  if (continuingFromFile)
  {
    InitMVNFromFile(continueFromDists,continueFromFile, allData, paramFilename);
    //MVNDist::Save(continueFromDists,"temp",allData.GetMask()); // check that the MVN created is right
      //MVNDist::Load(continueFromDists, continueFromFile, allData.GetMask());
  } 

  if (lockedLinearFile != "")
    throw Invalid_option("The option --locked-linear-from-mvn doesn't work with --method=vb yet, but should be pretty easy to implement.\n");
    

  const int nFwdParams = initialFwdPrior->GetSize();
  const int nNoiseParams = initialNoisePrior->OutputAsMVN().GetSize(); 

  // sort out loading for 'I' prior
  vector<ColumnVector> ImagePrior(nFwdParams);
#ifndef __FABBER_LIBRARYONLY
  volume4D<float> imagevol;
#endif
  for (int k=1; k<=nFwdParams; k++) {
    if (PriorsTypes[k-1] == 'I') {
      LOG_ERR("Reading Image prior ("<<k<<"): " << imagepriorstr[k-1] << endl);
      if (EasyOptions::UsingMatrixIO())
	{
	  ImagePrior[k-1] = EasyOptions::InMatrix(imagepriorstr[k-1]);
	}
      else
	{
#ifdef __FABBER_LIBRARYONLY
	  throw Logic_error("Should not reach this point!");
#else
	  read_volume4D(imagevol,imagepriorstr[k-1]);
	  ImagePrior[k-1] = (imagevol.matrix(allData.GetMask())).AsColumn();
#endif //__FABBER_LIBRARYONLY
      }
  }
 }

  // main loop over motion correction iterations and VB calculations
  bool continuefromprevious = false; //indicates that we should continue from a previous run (i.e. after a motion correction step)
  for (int step = 0; step <= Nmcstep; step++) {
    if (step>0) cout << endl << "Motion correction step " << step << " of " << Nmcstep << endl;

  // loop over voxels doing VB calculations
  for (int voxel = 1; voxel <= Nvoxels; voxel++)
    {
      ColumnVector y = data.Column(voxel);
      ColumnVector vcoords = coords.Column(voxel);
      if (suppdata.Ncols() > 0) {
	ColumnVector suppy = suppdata.Column(voxel);
	model->pass_in_data( y , suppy );
      }
      else {
	model->pass_in_data( y );
      }
      model->pass_in_coords(vcoords);
      NoiseParams* noiseVox = NULL;
      
      if (continuefromprevious) {
	// noise params come from resultMVN
	noiseVox = noise->NewParams();
	noiseVox->InputFromMVN( resultMVNs.at(voxel-1)->GetSubmatrix(nFwdParams+1, nFwdParams+nNoiseParams) );
      }
      else if (initialNoisePosterior == NULL) // continuing noise params from file 
      {
        assert(continuingFromFile);
        assert(continueFromDists.at(voxel-1)->GetSize() == nFwdParams+nNoiseParams);
        noiseVox = noise->NewParams();
        noiseVox->InputFromMVN( continueFromDists.at(voxel-1)
            ->GetSubmatrix(nFwdParams+1, nFwdParams+nNoiseParams) );
      }  
      else
      {
        noiseVox = initialNoisePosterior->Clone();
	/* if (continuingFromFile)
	   assert(continueFromDists.at(voxel-1)->GetSize() == nFwdParams);*/
      }
      const NoiseParams* noiseVoxPrior = initialNoisePrior;
      NoiseParams* const noiseVoxSave = noiseVox->Clone();
      

      // give an indication of the progress through the voxels
      LOG << "  Voxel " << voxel << " of " << Nvoxels << endl;
      if (fmod(voxel,floor(Nvoxels/10))==0) {cout << ". " << flush;}

      //LOG_ERR("  Voxel " << voxel << " of " << Nvoxels << endl); 
      //  << " sumsquares = " << (y.t() * y).AsScalar() << endl;
      double F = 1234.5678;

      MVNDist fwdPrior( *initialFwdPrior );
      MVNDist fwdPosterior;
      if (continuefromprevious) {
	//use result from a previous run within fabber (presumably after motion correction)
	fwdPosterior = resultMVNs.at(voxel-1)->GetSubmatrix(1, nFwdParams);
      }
      if (continuingFromFile)
      {
	//use results from a previous run loaded from a file
        assert(initialFwdPosterior == NULL);
        fwdPosterior = continueFromDists.at(voxel-1)->GetSubmatrix(1, nFwdParams);
      }
      else
      { 
        assert(initialFwdPosterior != NULL);
        fwdPosterior = *initialFwdPosterior;
	// any voxelwise initialisation
	model->Initialise(fwdPosterior);
      }


      MVNDist fwdPosteriorSave(fwdPosterior);
      MVNDist fwdPriorSave(fwdPrior);

      
      LinearizedFwdModel linear( model );
      
      // Setup for ARD (fwdmodel will decide if there is anything to be done)
      double Fard = 0;
      model->SetupARD( fwdPosterior, fwdPrior, Fard ); // THIS USES ARD IN THE MODEL AND IS DEPRECEATED
      Fard = noise->SetupARD( model->ardindices, fwdPosterior, fwdPrior );

      // Image priors
      for (int k=1; k<=nFwdParams; k++) {
	if (PriorsTypes[k-1] == 'I') {
	  ColumnVector thisimageprior;
	  thisimageprior = ImagePrior[k-1];
	  fwdPrior.means(k) = thisimageprior(voxel);
	}
      }

      try
	{
	  linear.ReCentre( fwdPosterior.means );
	  

	  noise->Precalculate( *noiseVox, *noiseVoxPrior, y );

	  conv->Reset();

	  // START the VB updates and run through the relevant iterations (according to the convergence testing)
	  int iteration = 0; //count the iterations
	  do 
	    {
	      if ( conv-> NeedRevert() ) //revert to previous solution if the convergence detector calls for it
		{
		  *noiseVox = *noiseVoxSave;  // copy values, not pointers!
		  fwdPosterior = fwdPosteriorSave;
		  fwdPrior = fwdPriorSave; // need to revert prior too (in case ARD is in place)
		  linear.ReCentre( fwdPosterior.means );
		}
	      
	      if (needF) { 
		F = noise->CalcFreeEnergy( *noiseVox, 
					   *noiseVoxPrior, fwdPosterior, fwdPrior, linear, y);
		F = F + Fard; }
	      if (printF) 
		LOG << "      Fbefore == " << F << endl;

 
              // Save old values if called for
	      if ( conv->NeedSave() )
              {
		*noiseVoxSave = *noiseVox;  // copy values, not pointers!
                fwdPosteriorSave = fwdPosterior;
		fwdPriorSave = fwdPrior;
              }

	      // Do ARD updates (model will decide if there is anything to do here)
	      if (iteration > 0) { 
		model->UpdateARD( fwdPosterior, fwdPrior, Fard ); // THIS USES ARD IN THE MODEL AND IS DEPRECEATED
		Fard = noise->UpdateARD( model->ardindices, fwdPosterior, fwdPrior );
	      }

	      // Theta update
	      noise->UpdateTheta( *noiseVox, fwdPosterior, fwdPrior, linear, y, NULL, conv->LMalpha() );


      
	      if (needF) {
		F = noise->CalcFreeEnergy( *noiseVox, 
					   *noiseVoxPrior, fwdPosterior, fwdPrior, linear, y);
		F = F + Fard; }
	      if (printF) 
		LOG << "      Ftheta == " << F << endl;
	      
	      
	      // Alpha & Phi updates
	      noise->UpdateNoise( *noiseVox, *noiseVoxPrior, fwdPosterior, linear, y );

	      if (needF) {
		F = noise->CalcFreeEnergy( *noiseVox, 
					   *noiseVoxPrior, fwdPosterior, fwdPrior, linear, y);
		F = F + Fard; }
	      if (printF) 
	      LOG << "      Fphi == " << F << endl;

	      // Test of NoiseModel cloning:
	      // NoiseModel* tmp = noise; noise = tmp->Clone(); delete tmp;

	      // Linearization update
	      // Update the linear model before doing Free eneergy calculation (and ready for next round of theta and phi updates)
	      linear.ReCentre( fwdPosterior.means );
	      
	      
	      if (needF) {
		F = noise->CalcFreeEnergy( *noiseVox, 
					   *noiseVoxPrior, fwdPosterior, fwdPrior, linear, y);
		F = F + Fard; }
	      if (printF) 
		LOG << "      Fnoise == " << F << endl;


	      iteration++;
	    }           
	  while ( !conv->Test( F ) );
	  // END of VB updates

	  // Revert to old values at last stage if required
	  if ( conv-> NeedRevert() )
          {
	    *noiseVox = *noiseVoxSave;  // copy values, not pointers!
            fwdPosterior = fwdPosteriorSave;
	    fwdPrior = fwdPriorSave;
	    linear.ReCentre( fwdPosterior.means ); //just in case we go on to use this in motion correction
	  }
	  conv->DumpTo(LOG, "    ");
	} 
      catch (const overflow_error& e)
	{
	  LOG_ERR("    Went infinite!  Reason:" << endl
		  << "      " << e.what() << endl);
	  //todo: write garbage or best guess to memory/file
	  if (haltOnBadVoxel) throw;
	  LOG_ERR("    Going on to the next voxel." << endl);
	}
      catch (Exception)
	{
	  LOG_ERR("    NEWMAT Exception in this voxel:\n"
		  << Exception::what() << endl);
	  if (haltOnBadVoxel) throw;
	  LOG_ERR("    Going on to the next voxel." << endl);  
	}
      catch (...)
	{
	  LOG_ERR("    Other exception caught in main calculation loop!!\n");
	    //<< "    Use --halt-on-bad-voxel for more details." << endl;
	  if (haltOnBadVoxel) throw;
	  LOG_ERR("    Going on to the next voxel" << endl);
	}
      
      // now write the results to resultMVNs
      try {

	LOG << "    Final parameter estimates (" << fwdPosterior.means.Nrows() << "x" << fwdPosterior.means.Ncols() << ") are: " << fwdPosterior.means.t() << endl;
	linear.DumpParameters(fwdPosterior.means, "      ");
	
	//assert(resultMVNs.at(voxel-1) == NULL); // this is no longer a good check, since we might voerwrite previous results here
	resultMVNs.at(voxel-1) = new MVNDist(
	  fwdPosterior, noiseVox->OutputAsMVN() );
	if (needF)
	  resultFs.at(voxel-1) = F;
	modelpred.Column(voxel) = linear.Offset(); // get the model prediction which is stored within the linearized forward model

      } catch (...) {
	// Even that can fail, due to results being singular
	LOG << "    Can't give any sensible answer for this voxel; outputting zero +- identity\n";
	MVNDist* tmp = new MVNDist();
	tmp->SetSize(fwdPosterior.means.Nrows()
		    + noiseVox->OutputAsMVN().means.Nrows());
	tmp->SetCovariance(IdentityMatrix(tmp->means.Nrows()));
	resultMVNs.at(voxel-1) = tmp;

	if (needF)
	  resultFs.at(voxel-1) = F;
	modelpred.Column(voxel) = linear.Offset(); // get the model prediction which is stored within the linearized forward model
      }
      
      delete noiseVox; noiseVox = NULL;
      delete noiseVoxSave;
    } //END of voxelwise updates

  //MOTION CORRECTION
  if (step<Nmcstep) { //dont do motion correction on the last run though as that would be a waste
#ifdef __FABBER_MOTION
     mcobj.run_mc(modelpred,data);
#endif //__FABBER_MOTION
  }

  continuefromprevious = true; //we now take resultMVNs and use these as the starting point if we are to run again
  }// END of Steps that include motion correction and VB updates

    while (continueFromDists.size()>0)
    {
      delete continueFromDists.back();
      continueFromDists.pop_back();
    }
}

VariationalBayesInferenceTechnique::~VariationalBayesInferenceTechnique() 
{ 
  delete conv;
  delete initialFwdPrior;
  delete initialFwdPosterior;
}




