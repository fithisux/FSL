/*  inference_spatialvb.cc - implementation of VB with spatial priors

Adrian Groves and Matthew Webster, FMRIB Image Analysis Group

Copyright (C) 2007-2010 University of Oxford  */

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

#include "easylog.h"
using namespace Utilities;
#include "inference_spatialvb.h"
#include "convergence.h"

#define NOCACHE 1

#ifndef __FABBER_LIBRARYONLY
using namespace NEWIMAGE;
#endif

void SpatialVariationalBayes::Setup(ArgsType& args)
{
  Tracer_Plus tr("SpatialVariationalBayes::Setup");
  // Call parent to do most of the setup
  VariationalBayesInferenceTechnique::Setup(args);

  try 
    {
      spatialDims = convertTo<int>(args.ReadWithDefault("spatial-dims", "3"));
    } 
  catch (invalid_argument&)
    {
      // More meaningful error message
      throw Invalid_option("--spatial-dims= must have an integer parameter");
    }

  if (spatialDims < 0 || spatialDims > 3)
    {
      throw Invalid_option("--spatial-dims= must take 0, 1, 2, or 3");
    }
  else if (spatialDims == 1)
    {
      Warning::IssueOnce("--spatial-dims=1 is very weird... I hope you're just testing something!");
    }
  else if (spatialDims == 2)
    {
      Warning::IssueOnce("--spatial-dims=2 doesn't decompose into slices and won't help if you're using the D prior");
    }

  maxPrecisionIncreasePerIteration = convertTo<double>(args.ReadWithDefault("spatial-speed","-1"));
  assert(maxPrecisionIncreasePerIteration > 1 || maxPrecisionIncreasePerIteration == -1);



  
  distanceMeasure = args.ReadWithDefault("distance-measure","dist1");
  spatialPriorsTypes = args.ReadWithDefault("param-spatial-priors","S+");

//  if (spatialPriorsTypes == "N+")
//    {
//      Warning::IssueOnce("You should specify --param-spatial-priors... check documentation!  Default is N+.");
//    }
  
  // Some unsupported options:
  fixedDelta = 
    convertTo<double>(args.ReadWithDefault("fixed-delta","-1")); 
  fixedRho = 
    convertTo<double>(args.ReadWithDefault("fixed-rho","0"));
  updateSpatialPriorOnFirstIteration = 
    args.ReadBool("update-spatial-prior-on-first-iteration");
  newDeltaEvaluations =
    convertTo<int>(args.ReadWithDefault("new-delta-iterations","10"));
  assert(newDeltaEvaluations > 0);

  // Some depreciated options:
  useSimultaneousEvidenceOptimization = args.ReadBool("use-simultaneous-evidence-optimization");
  useFullEvidenceOptimization = useSimultaneousEvidenceOptimization || args.ReadBool("use-full-evidence-optimization");
  firstParameterForFullEO = useFullEvidenceOptimization 
    ? convertTo<int>(args.ReadWithDefault("first-parameter-for-full-eo","1")) : -999; // WARNING: May need to be set to a sensible value in other circumstances!
  useEvidenceOptimization = useFullEvidenceOptimization || args.ReadBool("use-evidence-optimization");
  useCovarianceMarginalsRatherThanPrecisions = 
    useFullEvidenceOptimization && args.ReadBool("use-covariance-marginals");
  keepInterparameterCovariances = 
    useFullEvidenceOptimization && args.ReadBool("keep-interparameter-covariances");
  alwaysInitialDeltaGuess = convertTo<double>(args.ReadWithDefault("always-initial-delta-guess", "-1"));
  assert(!(updateSpatialPriorOnFirstIteration && !useEvidenceOptimization)); // currently doesn't work, but fixable
  bruteForceDeltaSearch = args.ReadBool("brute-force-delta-search");

  // Preferred way of using these options
  if (!useFullEvidenceOptimization && 
      !args.ReadBool("no-eo") &&
      spatialPriorsTypes.find_first_of("DR") != string::npos)
    {

      useFullEvidenceOptimization = true;
      useEvidenceOptimization = true;
      // keepInterparameterCovariances = true; // hacky
      useSimultaneousEvidenceOptimization = args.ReadBool("slow-eo");
      if (!useSimultaneousEvidenceOptimization)
	Warning::IssueOnce("Defaulting to Full (non-simultaneous) Evidence Optimization");
    }

  if (spatialPriorsTypes.find("F") != string::npos) // F found
    {
      if (fixedDelta < 0)
	throw Invalid_option("If --param-spatial-priors=F, you must specify a --fixed-delta value.\n");
    }
  else
    {
      if (fixedDelta == -1)
	fixedDelta = 0.5; // Default initial value (in mm!)
    }

  //get file names for I priors
  imagepriorstr.resize(model->NumParams());
  for (int k=1; k<=model->NumParams(); k++) {
    if (spatialPriorsTypes[k-1] == 'I') {
    imagepriorstr[k-1] = args.Read("image-prior"+stringify(k));
    }
  }

// deal with the spatial prior string, expand the '+' if it has been used
const int Nparams = model->NumParams();
const string::size_type thePlus = spatialPriorsTypes.find_first_of("+",1);
if (thePlus != string::npos)
{
assert(spatialPriorsTypes.find_last_of("+") == thePlus);
string before(spatialPriorsTypes, 0, thePlus-1);
string after(spatialPriorsTypes, thePlus+1, Nparams);
char repeatme = spatialPriorsTypes[thePlus-1];
//cout << "before == " << before << "\nafter == " << after
//   << "\nrepeatme == " << repeatme << "\nthePlus == " << thePlus
//   << endl;
spatialPriorsTypes = before;
for (int k = before.length()+after.length(); k < Nparams; k++)
  {spatialPriorsTypes += repeatme;   }
spatialPriorsTypes += after;


// deal with the shifting of image prior file names from expanding the +
 int nins = Nparams-before.length()-after.length()-2; // -2 accoutns for the letter and + in the original string
  for (unsigned int k = Nparams; k> (Nparams-after.length()) ; k--)
   {
     imagepriorstr[k-1] = imagepriorstr[k-1-nins];
     if (nins>0) imagepriorstr[k-1-nins] = ""; //clear old entry
   }
 }

// Deal with priors specified on the command line using the PSP_byname syntax
// NOTE: the section the queries the command line options is now in the inferencevb steup routine
// Here we just copy across any entries from PriorsTypes to spatialPriorsTypes
 if ( !PSPidx.empty() ) {
   for (unsigned int k=0; k<PSPidx.size(); k++) {
     spatialPriorsTypes[PSPidx[k]] = PriorsTypes[PSPidx[k]];
     //cout << "Copy " << k << " for parameter " << PSPidx[k] << endl;
   }
 }

// Section in which spatial priors can be specified by parameter name on command line
// param_spatial_priors_bymane (PSP_byname)
// this overwrites any existing entries in the string
 // vector<string> modnames; //names of model parameters
 // model->NameParams(modnames);
 // int npspnames=0;
 // while (true)
 //   {
 //     npspnames++;
 //     cout << "PSP_byname"+stringify(npspnames) << endl;
 //     string bytypeidx = args.ReadWithDefault("PSP_byname"+stringify(npspnames),"stop!");
 //     LOG_ERR("PSP_byname: " << bytypeidx << endl);
 //     if (bytypeidx == "stop!") break; //no more spriors have been specified

 //     // deal with the specification of this sprior (by name)
 //     //compare name to those in list of model names
 //     bool found=false;
 //     for (int p=0; p<model->NumParams(); p++) {
 //       if (bytypeidx == modnames[p]) {
 // 	 found = true;
 // 	 char pspstype = convertTo<char>(args.Read("PSP_byname"+stringify(npspnames)+"_type"));
 // 	 spatialPriorsTypes[p]=pspstype;
 // 	 

 // 	 // now read in file name for an image prior (if appropriate)
 // 	 if (pspstype == 'I') {
 // 	   imagepriorstr[p] = args.Read("PSP_byname"+stringify(npspnames)+"_image");
 // 	 }
 //       }
 //     }
 //     if (!found) {
 //       throw Invalid_option("ERROR: Spatial prior specification by name, parameter " + bytypeidx + " does not exist in the model\n");
 //   }
 //   }


// finally check that there are the right number of spatial priors specified and write full expanded string to the log
if ((int)spatialPriorsTypes.length() != Nparams)// && !useShrinkageMethod)
{
throw Invalid_option("--param-spatial-priors=" + spatialPriorsTypes 
		   + ", but there are " + stringify(Nparams)
		   + " parameters!\n");
}
else
{
LOG_ERR("Expanded, --param-spatial-priors=" 
      << spatialPriorsTypes << endl);
}

}

void SpatialVariationalBayes::DoCalculations(const DataSet& allData)
{
Tracer_Plus tr("SpatialVariationalBayes::DoCalculations");
const Matrix& data = allData.GetVoxelData();
const Matrix & coords = allData.GetVoxelCoords();
const Matrix & suppdata = allData.GetVoxelSuppData();
const int Nvoxels = data.Ncols();
// Rows are volumes
// Columns are (time) series
// num Rows is size of (time) series
// num Cols is size of volumes       

  // pass in some (dummy) data/coords here just in case the model relies upon it
  // use the first voxel values as our dummies
  if (suppdata.Ncols() > 0) {
    model->pass_in_data( data.Column(1) , suppdata.Column(1) );
  }
  else {
    model->pass_in_data( data.Column(1) );
  }
  model->pass_in_coords(coords.Column(1));


const int Nparams = model->NumParams();

// Added to diagonal to make sure the spatial precision matrix
// doesn't become singular -- and isolated voxels behave sensibly. 
const double tiny = 0; // turns out to be no longer necessary.

// Sanity checks:

if (data.Nrows() != model->NumOutputs())
throw Invalid_option("Data length (" 
+ stringify(data.Nrows())
+ ") does not match model's output length ("
+ stringify(model->NumOutputs())
+ ")!");

assert(resultMVNs.empty()); // Only call DoCalculations once
assert(resultMVNsWithoutPrior.empty());;
assert(resultFs.empty());

// Initialization:

// Make the neighbours[] lists if required
if (spatialPriorsTypes.find_first_of("mMpPSZ") != string::npos)
{
#ifndef __FABBER_LIBRARYONLY
	if (allData.GetMask().nvoxels()>0)
  	    CalcNeighbours(allData.GetMask());
	else 
#endif // __FABBER_LIBRARYONLY
	    CalcNeighbours(allData.GetVoxelCoords());
}

// Make distance matrix if required
if (spatialPriorsTypes.find_first_of("RDF") != string::npos)
{
#ifndef __FABBER_LIBRARYONLY
	if (allData.GetMask().nvoxels()>0)
  	    covar.CalcDistances(allData.GetMask(), distanceMeasure);
	else
#endif //__FABBER_LIBRARYONLY
	    covar.CalcDistances(allData.GetVoxelCoords(), distanceMeasure); // Note: really ought to know the voxel dimensions and multiply by those, because CalcDistances expects an input in mm, not index.
}

// If we haven'd done this, then covar is invalid and it'll return a 
// zero-size matrix for GetC, etc!

// Make each voxel's distributions

vector<NoiseParams*> noiseVox; // these change
vector<NoiseParams*> noiseVoxPrior; // these may change in future
vector<MVNDist> fwdPriorVox;
vector<MVNDist> fwdPosteriorVox;
vector<LinearizedFwdModel> linearVox;

 bool alsoSaveWithoutPrior = useEvidenceOptimization; // or other reasons?
 bool alsoSaveSpatialPriors = false;
 Warning::IssueOnce("Not saving finalSpatialPriors.nii.gz -- too huge!!");

 vector<MVNDist*> fwdPosteriorWithoutPrior(Nvoxels, (MVNDist*)NULL);
  if ( alsoSaveWithoutPrior )
    for (int v = 1; v <= Nvoxels; v++)
      fwdPosteriorWithoutPrior.at(v-1) = new MVNDist;

// Locked linearizations, if requested
bool lockedLinearEnabled = (lockedLinearFile != "");
Matrix lockedLinearCentres;  // empty by default

{ Tracer_Plus tr("SpatialVariationalBayes::DoCalculations - initialization");

// If we're continuing from previous saved results, load them here:
continuingFromFile = (continueFromFile != "");

vector<MVNDist*> continueFromDists;
if (continuingFromFile)
{
  InitMVNFromFile(continueFromDists,continueFromFile, allData, paramFilename);
  //MVNDist::Load(continueFromDists, continueFromFile, allData.GetMask());
}

// Locked linearizations, if requested
if (lockedLinearEnabled)
{
LOG_ERR("Loading fixed linearization centres from the MVN '" 
      << lockedLinearFile << "'\nNOTE: This does not check if the correct number of parameters is present!\n");
vector<MVNDist*> lockedLinearDists;
#ifndef __FABBER_LIBRARYONLY
MVNDist::Load(lockedLinearDists, lockedLinearFile, allData.GetMask());
#else
throw Logic_error("lockedLinearEnabled not supported yet for fabber_library");
#endif
lockedLinearCentres.ReSize(Nparams, Nvoxels);

for (int v = 1; v <= Nvoxels; v++)
{
  lockedLinearCentres.Column(v) = 
    lockedLinearDists.at(v-1)->means.Rows(1,Nparams);
}
}

const int nFwdParams = initialFwdPrior->GetSize();
const int nNoiseParams = initialNoisePrior->OutputAsMVN().GetSize();   


noiseVox.resize(Nvoxels, NULL); // polymorphic type, so need to use pointers
noiseVoxPrior.resize(Nvoxels,NULL);
fwdPriorVox.resize(Nvoxels, *initialFwdPrior);
if (continuingFromFile)
{
  fwdPosteriorVox.resize(Nvoxels);
  for (int v = 1; v <= Nvoxels; v++)
    {
      fwdPosteriorVox[v-1] = continueFromDists.at(v-1)
	->GetSubmatrix(1, nFwdParams);
    }
 }
else
{ 
  fwdPosteriorVox.resize(Nvoxels, *initialFwdPosterior);
}

linearVox.resize(Nvoxels, LinearizedFwdModel(model) );
resultMVNs.resize(Nvoxels, NULL);

if (alsoSaveWithoutPrior)
  resultMVNsWithoutPrior.resize(Nvoxels, NULL);

resultFs.resize(Nvoxels, 9999); // 9999 is a garbage default value

for (int v = 1; v <= Nvoxels; v++)
{
linearVox[v-1].ReCentre(lockedLinearEnabled
		      ? lockedLinearCentres.Column(v)
		      : fwdPosteriorVox[v-1].means 
		      );

if (initialNoisePosterior == NULL) // continuing Noise from file
{
assert(nFwdParams + nNoiseParams == continueFromDists.at(v-1)->GetSize());
noiseVox[v-1] = noise->NewParams();
noiseVox[v-1]->InputFromMVN( continueFromDists.at(v-1)
    ->GetSubmatrix(nFwdParams+1, nFwdParams+nNoiseParams) );
}
else
{ 
noiseVox[v-1] = initialNoisePosterior->Clone();
}
noiseVoxPrior[v-1] = initialNoisePrior->Clone();
noise->Precalculate( *noiseVox[v-1], *noiseVoxPrior[v-1], data.Column(v) );
}
} // end tracer  

// Make the spatial normalization parameters
//akmean = 0*distsMaster.theta.means + 1e-8;
DiagonalMatrix akmean(Nparams); akmean = 1e-8;




DiagonalMatrix delta(Nparams); 
DiagonalMatrix rho(Nparams);
delta = fixedDelta; // Hard-coded initial value (in mm!)
rho = 0;
LOG_ERR("Using initial value for all deltas: " << delta(1) << endl);
//  delta(1) = delta(3) = .5;
//  LOG_ERR("Except delta([1 3]) (Q0,M0) = " << delta(3) << endl);
//  delta(3) = .5;
//  LOG_ERR("Except delta(3) (M0) = " << delta(3) << endl);
//delta(7) = 1e12;
//LOG_ERR("Except delta(7) (dt) = " << delta(7) << endl);
//  bool HackDelta7NextTime = true;
vector<SymmetricMatrix> Sinvs(Nparams);

SymmetricMatrix StS; // Cache for StS matrix in 'S' mode

const double globalF = 1234.5678; // no sensible updates yet
//  if (needF)
//    throw Invalid_option("Can't calculate free energy with --inference=spatialvb yet!\n");
// Allow calculation to continue even with bad voxels
// Note that this is a bad idea when spatialDims>0, because the bad voxel
// will drag its neighbours around... but can never recover!  Maybe a more
// sensible approach is to reset bad voxels to the prior on each iteration.
//  vector<bool> keepGoing;
//  keepGoing.resize(Nvoxels, true);


// sort out loading for 'I' prior
vector<ColumnVector> ImagePrior(Nparams);
#ifndef __FABBER_LIBRARYONLY
volume4D<float> imagevol;
#endif
for (int k=1; k<=Nparams; k++) {
  if (spatialPriorsTypes[k-1] == 'I') {
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


// Quick check.. which shrinkage prior to use?  For various reasons we can't 
// yet mix different shinkage priors for different parameters (and I can't
// see any good reason to implement it).
char shrinkageType = '-';
for (int k = 1; k <= Nparams; k++)
{
char type = spatialPriorsTypes[k-1];
switch (type)
{
 case 'R': case 'D': case 'N': case 'F': case 'I': case 'A':
  break;

 case 'm': case 'M': case 'p': case 'P': case 'S':
  if (type != shrinkageType && shrinkageType != '-')
    throw Invalid_option("Sorry, only one type of shrinkage prior at a time, please!\n");
  shrinkageType = type;
  break;

default:
  LOG_ERR("What the heck? spatialPriorsType[" << k << "-1] == "
	  << spatialPriorsTypes[k-1] << endl);
  assert(false);
}
}

if (StS.Nrows() == 0 && neighbours.size() > 0 && (shrinkageType == 'S' || shrinkageType == 'Z'))
  {
    Tracer_Plus tr("Creating StS matrix");
    assert((int)neighbours.size() == Nvoxels);

    const double tiny = 1e-6;
    Warning::IssueOnce("Using 'S' prior with fast-calculation method and constant diagonal weight of " + stringify(tiny));

    // NEW METHOD
    { Tracer_Plus tr("New method for generating StS matrix");
    cout << "Attempting to allocate, Nvoxels = " << Nvoxels << endl;
    StS.ReSize(Nvoxels);
    cout << "Allocated" << endl;
    StS = 0;
    for (int v = 1; v <= Nvoxels; v++)
      {
        int Nv = neighbours[v-1].size(); // Number of neighbours v has

        // Diagonal value = N + (N+tiny)^2
        StS(v,v) = Nv + (Nv+tiny)*(Nv+tiny);

	// Off-diagonal value = num 2nd-order neighbours (with duplicates) - Aij(Ni+Nj+2*tiny)
	for (vector<int>::iterator nidIt = neighbours[v-1].begin();
             nidIt != neighbours[v-1].end(); nidIt++)
          {
            if (v < *nidIt)
              StS(v,*nidIt) -= Nv + neighbours[*nidIt-1].size() + 2*tiny;
          }
	for (vector<int>::iterator nidIt = neighbours2[v-1].begin();
             nidIt != neighbours2[v-1].end(); nidIt++)
          {
            if (v < *nidIt)
              StS(v,*nidIt) += 1;
          }
      }
      cout << "Done generating StS matrix (New method)" << endl;
    } // end NEW METHOD tracer block

    /* OLD METHOD
    // Slow because connect*connect is O(N^3)	      	    
    SymmetricMatrix connect =
      IdentityMatrix(Nvoxels) * 1e-6;
    Warning::IssueOnce("Using 'S' prior with constant diagonal weight of " + stringify(connect(1,1)));
    
    for (int v = 1; v <= Nvoxels; v++)
      {
	for (vector<int>::iterator nidIt = neighbours[v-1].begin();
	     nidIt != neighbours[v-1].end(); nidIt++) 
	  {
	    int nid = *nidIt;
	    connect(v, nid) = -1;
	    connect(v,v) += 1;
	  }
      }  
    
    Tracer_Plus tr2("connect*connect");
    
    Matrix sTmp = connect*connect;
    assert(sTmp == sTmp.t());
    // Skip this, since we're just doing a validation check: StS << sTmp;

    // VALIDATE THAT OLD == NEW
    for (int i = 1; i <= Nvoxels; i++)
      for (int j = 1; j <= Nvoxels; j++)
        if (sTmp(i,j) != StS(i,j))
          cout << "Mismatch at " << i << "," << j << ": " << sTmp(i,j) << "," << StS(i,j) << endl;

    assert(sTmp == StS);

    // END OLD METHOD */
  }


conv->Reset();
bool isFirstIteration = true; // slightly different behaviour in first iteratio

//  if (!useShrinkageMethod) LOG_ERR("HACK: using --fixed-delta value on first iteration instead of automatically determining delta from priors\n");

// MAIN ITERATION LOOP
do {
Tracer_Plus tr("Main iteration loop");
conv->DumpTo(LOG);
conv->DumpTo(cout);    

// UPDATE SPATIAL SHRINKAGE PRIOR PARAMETERS

    //    if (useShrinkageMethod)
 if (shrinkageType != '-' && (!isFirstIteration || updateSpatialPriorOnFirstIteration))
      { 
	Tracer_Plus tr("SpatialVariationalBayes::DoCalculations - old spatial norm update");
	// Update spatial normalization term
       
	// Collect gk, wk, sigmak across all voxels
	DiagonalMatrix gk(Nparams); //gk = 0.0/0.0;
	for (int k = 1; k <= Nparams; k++)
	  {
	    ColumnVector wk(Nvoxels); //wk = 0.0/0.0;
	    DiagonalMatrix sigmak(Nvoxels); //sigmak=0.0/0.0;
	    for (int v = 1; v <= Nvoxels; v++)
	      {
		wk(v) = fwdPosteriorVox.at(v-1).means(k);
		sigmak(v,v) = fwdPosteriorVox.at(v-1).GetCovariance()(k,k);
	      }
	    //    wk = wk(:);
	    //    tmp1 = trace(diag(sigmak)*S'*S)
	    
	    //cout << "fwdPosteriorVox[0].GetCovariance():" << fwdPosteriorVox[0].GetCovariance();
	    

	    // Update from Penny05:
	    // 1/gk = 0.5*Trace[Sigmak*S'*S] + 0.5*wk'*S'*S*wk + 1/q1
	    // hk = N/2 + q2
	    // To do the MRF, just replace S'*S with S.
	    
	    switch (shrinkageType)
	      {
	      case 'Z': //case 'S':

		{
		  Tracer_Plus tr("S-prior akmean estimation");
		  assert(alsoSaveWithoutPrior);
		  
		  assert(StS.Nrows() == Nvoxels);

		  // Prior used by penny:
		  //		  double q1 = 10, q2 = 1;
		  // Noninformative prior:
		  double q1 = 1e12, q2 = 1e-12;

		  Warning::IssueOnce("Hyperpriors on S prior: using q1 == " + stringify(q1) + ", q2 == " + stringify(q2));

		  gk(k) = 1/( 0.5*(sigmak * StS).Trace() + (wk.t()*StS*wk).AsScalar() + 1/q1);
		  
		  akmean(k) = gk(k) * (0.5*Nvoxels + q2);
		}
		// only output is akmean(k)
		break;

	      case 'p': case 'P': case 'm': case 'M': case 'S':

		{
		  // The following calculates Tr[Sigmak*S'*S]
		  // using the fact that this == sum(diag(sigmak) .* diag(S'*S))
		  // (since sigmak is diagonal!)
		  
		  double tmp1 = 0.0;	
		  for (int v = 1; v <= Nvoxels; v++)
		    {
		      int nn = neighbours.at(v-1).size();
		      //cout << v << ": " << nn << "," << wk(v) << "," << sigmak(v,v) << endl;
		      if (shrinkageType == 'm') //useMRF)
			tmp1 += sigmak(v,v) * spatialDims*2;
		      else if (shrinkageType == 'M') //useMRF2)
			tmp1 += sigmak(v,v) * (nn+1e-8);
		      else if (shrinkageType == 'p')
			tmp1 += sigmak(v,v) * ( 4*spatialDims*spatialDims + nn );
		      else if (shrinkageType == 'S')
			tmp1 += sigmak(v,v) * ( (nn+1e-6)*(nn+1e-6) + nn );
		      else
			tmp1 += sigmak(v,v) * ( (nn+tiny)*(nn+tiny) + nn );
		    }
		  
		  //    tmp2 = wk'*S'*S*wk
		  ColumnVector Swk = tiny * wk;
		  
		  if (shrinkageType == 'S')
		    Swk = 1e-6 * wk;

		  for (int v = 1; v <= Nvoxels; v++)
		    {
		      for (vector<int>::iterator v2It = neighbours[v-1].begin();
			   v2It != neighbours.at(v-1).end(); v2It++)
			{
			  Swk(v) += wk(v) - wk(*v2It);
			}
		      //		if (useDirichletBC || useMRF) // but not useMRF2
		      if (shrinkageType == 'p' || shrinkageType == 'm')
			Swk(v) += wk(v)*(spatialDims*2 - neighbours.at(v-1).size());
		      // Do nothing for 'S'
		    }
		  double tmp2 = Swk.SumSquare(); //(Swk.t() * Swk).AsScalar();
		  
		  //	    if (useMRF || useMRF2) // overwrite this for MRF
		  if (shrinkageType == 'm' || shrinkageType == 'M')
		    tmp2 = DotProduct(Swk, wk);
		  
		  cout << "k=" << k << ", tmp1=" << tmp1 << ", tmp2=" << tmp2 << endl;
		  //cout << Swk.t();
		  
		  //  gk(k) = 1/(0.5*tmp1 + 0.5*tmp2 + 1/10)
		  gk(k,k) = 1/(0.5*tmp1 + 0.5*tmp2 + 0.1); // prior q1 == 10 (1/q1 == 0.1)
		  //  end
		  
		  akmean(k) = gk(k) * (Nvoxels*0.5 + 1.0); // prior q2 == 1.0
		}
		
		break;
		
	      default:
		assert(false);
		break;
	      }
	  }
		
	      	DiagonalMatrix akmeanMax = akmean * maxPrecisionIncreasePerIteration;
     
	for (int k = 1; k <= akmean.Nrows(); k++)
	  {
	    if (akmean(k,k) < 1e-50)
	    {
	      LOG_ERR("akmean(" << k << ") was " << akmean(k,k) << endl);
	      Warning::IssueOnce("akmean value was tiny!");
	      akmean(k,k) = 1e-50; // prevent crashes
	    }
	  }
      
	// Rate limiting section - temporary turn back on (also akmeanMax line above)
		for (int k = 1; k <= akmean.Nrows(); k++)
		  {
		    if (akmeanMax(k) < 0.5)
		      akmeanMax(k) = 0.5;  // totally the wrong scaling.. oh well
		    if (maxPrecisionIncreasePerIteration > 0 && akmean(k) > akmeanMax(k))
		      {
			LOG_ERR("Rate-limiting the increase on akmean " << k 
				<< ": was " << akmean(k));
			akmean(k) = akmeanMax(k);
			LOG_ERR(", now " << akmean(k) << endl);
		      }
		  }
	

	//  cout << "Using fixed values for akmean!!\n";
	//  akmean << 0.0011 << 9.9704e-04 << 2.1570e-07 
	//    << 5.2722 << 0.0177 << 0.5041;

	  LOG_ERR("New akmean: " << akmean.AsColumn().t());
//	for (int asdf = 1; asdf <= akmean.Nrows(); asdf++)
//	  LOG_ERR(akmean(asdf) << endl);
	
      } // tracer: spatial norm update
    



// UPDATE DELTA & RHO ESTIMATES
for (int k = 1; k <= Nparams; k++)
{
Tracer tr("SpatialVariationalBayes::DoCalculations - delta updates");

//if(k<6){LOG_ERR("Skipping parameter "<<k<<endl);continue;}
        LOG_ERR("Optimizing for parameter " << k << endl);

	char type = spatialPriorsTypes[k-1];


	// Each type should issue exactly one line to the logfile of the form
	// SpatialPrior on k type ? ?? : x y z
	// ? = single-character type
	// ?? = any other subtype info (optional, free-form but no : character)
	// x, y, z = numerical parameters (e.g. delta, rho, 0)
	switch (type)
	  {
	  case 'N': case 'I': case 'A':
	    // Nonspatial priors
	    delta(k) = 0;
	    rho(k) = 0; // no shrinkage/scaling factor either	    
	    LOG_ERR("\nSpatialPrior " << k << " type " << type << " : 0 0 0\n");
	    break;

	    //	  case 'F':
	    //	    delta(k) = fixedDelta;
	    //	    rho(k) = fixedRho;
	    //	    break;

	  case 'm': case 'M': case 'p': case 'P': case 'S':

	    assert(type == shrinkageType);
	    // fill with invalid values:
	    delta(k) = -3;
	    rho(k) = 1234.5678;
	    LOG_ERR("\nSpatialPrior " << k << " type " << type << " : " << akmean(k) << " 0 0\n");

	    break;
	    
	  default:
		throw Invalid_option(string("Invalid spatial prior type '")
			     + type + "' given to --param-spatial-priors\n");
		break;

	  case 'R': case 'D': case 'F':
	    // Reorganize data by parameter (rather than by voxel)
	    DiagonalMatrix covRatio(Nvoxels);
	    ColumnVector meanDiffRatio(Nvoxels);
	    const double priorCov = initialFwdPrior->GetCovariance()(k,k);
	    const double priorCovSqrt = sqrt(priorCov);
	    const double priorMean = initialFwdPrior->means(k);

	    
	    
	    for (int v = 1; v <= Nvoxels; v++)
	      {
		// Isolate just the dimensionless quantities we need

		// Penny:
	        covRatio(v,v) = fwdPosteriorVox.at(v-1).GetCovariance()(k,k) 
		  / priorCov;
		// Hacky:
		//LOG_ERR("WARNING: Using hacky covRatio calculation (precision rather than covariance!\n!");		
		//		covRatio(v,v) = 1 / fwdPosteriorVox.at(v-1).GetPrecisions()(k,k) / priorCov; 

		meanDiffRatio(v) = (fwdPosteriorVox.at(v-1).means(k) - priorMean) 
		  / priorCovSqrt;
	      }
	    
	    //	SymmetricMatrix covRatioSupplemented(Nvoxels);
	    // Recover the off-diagonal elements of fwdPriorVox/priorCov
	    
	    /* TURN COVRATIOSUPPLEMENTED BACK INTO COVRATIO!
	       
	    if (Sinvs.at(k-1).Nrows() != 0)
	    {
	    // Implicitly use all zeroes for first iteration
	    assert(Sinvs.at(k-1).Nrows() == Nvoxels);
	    covRatioSupplemented = Sinvs.at(k-1);
	    
	    for (int v=1; v<=Nvoxels; v++)
	    {
	    //		cout << "Should be one: " << initialFwdPrior->GetPrecisions()(k,k) * Sinvs.at(k-1)(v,v) / fwdPriorVox.at(v-1).GetPrecisions()(k,k) << endl;
	    assert(initialFwdPrior->GetPrecisions()(k,k)*Sinvs.at(k-1)(v,v) == fwdPriorVox.at(v-1).GetPrecisions()(k,k));
	    
	    covRatioSupplemented(v,v) = 0;
	    }
	    }
	    
	    covRatioSupplemented += covRatio.i();
	    covRatioSupplemented = covRatioSupplemented.i();
	    */

	    // INSTEAD DO THIS:
	    //	covRatioSupplemented = covRatio;
	    
	    if (isFirstIteration && !updateSpatialPriorOnFirstIteration)
	      {
		if (type == 'F' && bruteForceDeltaSearch)
		  LOG_ERR("Doing calc on first iteration, just because it's F and bruteForceDeltaSearch is on.  Temporary hack!\n");
		else
		  break; // skip the updates
	      }

	    DiagonalMatrix deltaMax = delta * maxPrecisionIncreasePerIteration;
	    

	    if (type == 'R')
	      {
		if (alwaysInitialDeltaGuess>0) delta(k) = alwaysInitialDeltaGuess;
		if (useEvidenceOptimization)
		  {
		    Warning::IssueAlways("Using R... mistake??");
		    delta(k) = OptimizeEvidence(fwdPosteriorWithoutPrior, k, initialFwdPrior, delta(k), true, &rho(k));
		    LOG_ERR("\nSpatialPrior " << k << " type R eo : " << delta(k) << " " << rho(k) << " 0\n");		  }
		else
		  {		    
		    Warning::IssueAlways("Using R without EO... mistake??");
		    // Spatial priors with rho & delta
		    delta(k) = OptimizeSmoothingScale( covRatio, meanDiffRatio,
						       delta(k), &rho(k), true);
		    LOG_ERR("\nSpatialPrior " << k << " type R vb : " << delta(k) << " " << rho(k) << " 0\n");

		  }
	      }
	    else if (type == 'D')
	      {
                if (alwaysInitialDeltaGuess>0) delta(k) = alwaysInitialDeltaGuess;

		// Spatial priors with only delta
		if (useEvidenceOptimization)
		  {
		    delta(k) = OptimizeEvidence(fwdPosteriorWithoutPrior, k, initialFwdPrior, delta(k));
		    LOG_ERR("\nSpatialPrior " << k << " type D eo : " << delta(k) << " 0 0\n");
		  }
		else
		  {
		    Warning::IssueAlways("Using D without EO... mistake??");
		    delta(k) = OptimizeSmoothingScale( covRatio, meanDiffRatio, 
						       delta(k), &rho(k), false );
		    LOG_ERR("\nSpatialPrior " << k << " type D vb : " << delta(k) << " 0 0\n");

		  }
		
	      }
	    else // type == 'F'
	      {	
		delta(k) = fixedDelta;
		rho(k) = fixedRho;

		// The following does nothing BUT it's neccessary to 
		// make the bruteForceDeltaEstimates work.
		double newDelta = OptimizeSmoothingScale( 
				    covRatio, meanDiffRatio,
				    delta(k), &rho(k), 
				    false,
				    false);
		assert(newDelta == fixedDelta);
		assert(rho(k) == fixedRho);
		deltaMax(k) = delta(k);
		LOG_ERR("\nSpatialPrior " << k << " type F : " << delta(k) << " " << rho(k) << " 0\n");
	      }
		
	
	    // enforce maxPrecisionIncreasePerIteration
	    if (deltaMax(k) < 0.5)
	      deltaMax(k) = 0.5;
	    if (maxPrecisionIncreasePerIteration > 0 && delta(k) > deltaMax(k))
	      {
		LOG_ERR("Rate-limiting the increase on delta " << k 
			<< ": was " << delta(k));
		delta(k) = deltaMax(k);
		LOG_ERR(", now " << delta(k) << endl);
		
		// Re-evaluate rho, for this delta
		double newDelta = OptimizeSmoothingScale( 
							 covRatio, meanDiffRatio,
							 delta(k), &rho(k), 
							 type == 'R',
							 false);
		assert(newDelta == delta(k)); // a quick check
	      }
	    // default: dealt with earlier.
	  }
    
        LOG_ERR("    delta(k) = " << delta(k) << 
		", rho(k) == " << rho(k) << endl);
      }
    
    
    // CALCULATE THE C^-1 FOR THE NEW DELTAS
    {
      Tracer_Plus tr("SpatialVariationalBayes::DoCalculations - Cinv calculations");
      // Calculate the Cinv
      for (int k = 1; k <= Nparams; k++)
        { 
	  if (delta(k) >= 0)
	    {
	      Sinvs.at(k-1) = covar.GetCinv(delta(k)) * exp(rho(k));
	      
	      if (delta(k) == 0 && alsoSaveWithoutPrior)
		{
		  assert(spatialPriorsTypes[k-1] == 'N' | spatialPriorsTypes[k-1] == 'I' | spatialPriorsTypes[k-1] == 'A');
		  Sinvs.at(k-1) = IdentityMatrix(Nvoxels); 
		}

	      assert( SP( initialFwdPrior->GetPrecisions(), IdentityMatrix(Nparams)-1 ).MaximumAbsoluteValue() == 0 );
	      Sinvs[k-1] *= initialFwdPrior->GetPrecisions()(k,k);
	    }
	  
          if (delta(k)<0 && alsoSaveWithoutPrior)
	    {
	      Tracer_Plus tr("Building spatial precision matrix for Penny prior");
	      assert(spatialPriorsTypes[k-1] == shrinkageType);
	      
	      if (shrinkageType == 'S')
		{ 
		  assert(StS.Nrows() == Nvoxels);
		  Sinvs.at(k-1) = StS * akmean(k);
		}
	      else
		{
		  assert(shrinkageType == 'p');

		  // Build up the second-order matrix directly, row-by-row
		  Matrix sTmp(Nvoxels, Nvoxels); 
		  sTmp = 0;
		  
		  for (int v = 1; v <= Nvoxels; v++)
		    {
		      
		      // self = (2*Ndim)^2 + (nn)
		      sTmp(v,v) = 4*spatialDims*spatialDims; // nn added later
		      
		      // neighbours = (2*Ndim) * -2
		      for (vector<int>::iterator nidIt = neighbours[v-1].begin();
			   nidIt != neighbours[v-1].end(); nidIt++)
			{
			  int nid = *nidIt; // neighbour ID (voxel number)
			  assert(sTmp(v,nid) == 0);
			  sTmp(v,nid) = -2 * 2 * spatialDims;
			  sTmp(v,v) += 1;
			}
		      
		      // neighbours2 = 1 (for each appearance)	    
		      for (vector<int>::iterator nidIt = neighbours2[v-1].begin();
			   nidIt != neighbours2[v-1].end(); nidIt++)
			{
			  int nid2 = *nidIt; // neighbour ID (voxel number)
			  sTmp(v,nid2) += 1; // not =1, because duplicates are ok.
			}
		    }
		  
		  // Store back into Sinvs[k-1] and apply akmean(k)
		  assert(sTmp == sTmp.t());
		  Sinvs.at(k-1) << sTmp;
		  Sinvs[k-1] *= akmean(k);
		}
	    }
	}
    }        
    

    // ITERATE OVER VOXELS
    for (int v = 1; v <= Nvoxels; v++)
      {
	// some models may want extra information about the data
	if (suppdata.Ncols() > 0) {
	  model->pass_in_data(  data.Column(v) ,  suppdata.Column(v) );
	}
	else {
	  model->pass_in_data(  data.Column(v) );
	}
	model->pass_in_coords(coords.Column(v));
	double &F = resultFs.at(v-1);  // short name

	if (!continuingFromFile) {
	  //voxelwise initialisation - only if we dont have initial values from a pre loaded MVN
	  model->Initialise(fwdPosteriorVox[v-1]);
	}

	// from simple_do_vb_ar1c_spatial.m

	// Note: this sets the priors as if all parameters were shrinkageType.
	// We overwrite the non-shrinkageType parameter priors later.

	if (shrinkageType == 'S')
	  {
	    Tracer_Plus tr("shrinkage spatial priors S");
	    Warning::IssueOnce("Using new S VB spatial thingy");

	    assert(StS.Nrows() == Nvoxels);

	    double weight = 1e-6; // weakly pulled to zero
	    ColumnVector contrib(Nparams); 
	    contrib = 0;
	    
	    for (int i = 1; i <= Nvoxels; i++)
	      {
		if (v != i)
		  {
		    weight += StS(v,i);
		    contrib += StS(v,i) * fwdPosteriorVox[i-1].means;
		  }
	      }
	    
	    DiagonalMatrix spatialPrecisions;
	    spatialPrecisions = akmean * StS(v,v);
	    
	    fwdPriorVox[v-1].SetPrecisions(spatialPrecisions);

	    fwdPriorVox[v-1].means = contrib / weight;   

	  }
	else if (shrinkageType != '-')
	  { 
	    Tracer_Plus tr("SpatialVariationalBayes::DoCalculations - shrinkage spatial priors");
	    



	    double weight8 = 0; // weighted +8
	    ColumnVector contrib8(Nparams); contrib8 = 0.0;
	    for (vector<int>::iterator nidIt = neighbours[v-1].begin();
		 nidIt != neighbours[v-1].end(); nidIt++) 
	      // iterate over neighbour ids
	      {
		int nid = *nidIt;
		const MVNDist& neighbourPost = fwdPosteriorVox[nid-1];
		contrib8 += 8 * neighbourPost.means;
		weight8 += 8;
	      }
	    
	    double weight12 = 0; // weighted -1, may be duplicated
	    ColumnVector contrib12(Nparams); contrib12 = 0.0;
	    for (vector<int>::iterator nidIt = neighbours2[v-1].begin();
		 nidIt != neighbours2[v-1].end(); nidIt++)
	      // iterate over neighbour ids
	      {
		int nid = *nidIt;
		const MVNDist& neighbourPost = fwdPosteriorVox[nid-1];
		contrib12 += -neighbourPost.means;
		weight12 += -1;
	      }
	    
	    // Set prior mean & precisions
	    
	    int nn = neighbours[v-1].size();
	    
	    //	    if (useDirichletBC)
	    if (shrinkageType == 'p')
	      {
		//	cout << nn << " -> " << 2*spatialDims << endl;
		assert(nn <= spatialDims*2);
		//nn = spatialDims*2;
		weight8 = 8*2*spatialDims;
		weight12 = -1*(4*spatialDims*spatialDims-nn);
	      }
	    
	    DiagonalMatrix spatialPrecisions;

	    if (shrinkageType == 'P')
	      spatialPrecisions = 
		akmean * ( (nn+tiny)*(nn+tiny) + nn );
	    else if (shrinkageType == 'm')
	      spatialPrecisions = 
		akmean * spatialDims*2;
	    else if (shrinkageType == 'M')
	      spatialPrecisions = 
		akmean * (nn+1e-8);
	    else if (shrinkageType == 'p')
	      spatialPrecisions = 
		akmean * (4*spatialDims*spatialDims + nn);
	    else if (shrinkageType == 'S')
	      {
		spatialPrecisions =
		  akmean * ( (nn+1e-6)*(nn+1e-6) + nn );
		Warning::IssueOnce("Using a hacked-together VB version of the 'S' prior");
	      }

	    //	    if (useDirichletBC || useMRF)
	    if (shrinkageType == 'p' || shrinkageType == 'm')
	      {
		//	LOG_ERR("Penny-style DirichletBC priors -- ignoring initialFwdPrior completely!\n");
		fwdPriorVox[v-1].SetPrecisions(spatialPrecisions);
	      }  
	    else
	      {
		fwdPriorVox[v-1].SetPrecisions(
			 initialFwdPrior->GetPrecisions() + spatialPrecisions );
	      }
	    
	    ColumnVector mTmp(Nparams);
	    
	    if (weight8 != 0)
	      mTmp = (contrib8+contrib12)/(weight8+weight12);
	    else
	      mTmp = 0;
	    
	    if (shrinkageType == 'm') //useMRF) // overwrite this for MRF
	      mTmp = contrib8 / (8*spatialDims*2); // note: Dirichlet BCs on MRF
	    if (shrinkageType == 'M') //useMRF2)
	      mTmp = contrib8 / (8*(nn+1e-8));
	    
	    // equivalent, when non-spatial priors are very weak:
	    //    fwdPriorVox[v-1].means = mTmp; 
	    
	    fwdPriorVox[v-1].means =
	      fwdPriorVox[v-1].GetCovariance() *  
	      (spatialPrecisions * mTmp 
	       + initialFwdPrior->GetPrecisions() * initialFwdPrior->means);
	    
	    //	    if (useMRF || useMRF2) // overwrite this for MRF
	    if (shrinkageType == 'm' || shrinkageType == 'M')
	      fwdPriorVox[v-1].means = 
		fwdPriorVox[v-1].GetCovariance() * spatialPrecisions * mTmp; // = mTmp;
	    


	  } 
	// else

//cout << "Sinvs[0] is:\n" << Sinvs[0] << endl; // Verified against matlab 2008-02-16

	double Fard=0;
	if (1)
	  { 
	  
	  // Use the new spatial priors
	  Tracer_Plus tr("SpatialVariationalBayes::DoCalculations - new spatial prior calculations");

	  // Marginalize out all the other voxels
	  
	  DiagonalMatrix spatialPrecisions(Nparams);
	  ColumnVector weightedMeans(Nparams);

	  ColumnVector priorMeans(Nparams);
	  priorMeans = initialFwdPrior->means; // default is to get these from intialFwdPrior
	                                       // this is overwritten for I priors 
	                                       // or ignored for spatial priors
	 
	  for (int k = 1; k <= Nparams; k++) 
	    {
	      if (spatialPriorsTypes[k-1] == shrinkageType)
		{
		  spatialPrecisions(k) = -9999;
		  weightedMeans(k) = -9999;
		  continue;
		}
	      else if (spatialPriorsTypes[k-1] == 'A')
		{
		  if (isFirstIteration)
		    {
		      spatialPrecisions(k) = initialFwdPrior->GetPrecisions()(k,k);
		      weightedMeans(k) = initialFwdPrior->means(k);
		      //Fard = 0;		
		    }		  
		  else
		    {
		      double ARDparam = 1/fwdPosteriorVox[v-1].GetPrecisions()(k,k) + 
			fwdPosteriorVox[v-1].means(k)*fwdPosteriorVox[v-1].means(k) ;
		      spatialPrecisions(k) = 1/ARDparam;
		      weightedMeans(k) = 0;
		      Fard -= 2.0*log(2.0/ARDparam);
		    }
		  continue;
		}
	      else if (spatialPriorsTypes[k-1] == 'N')
		{
		  // special case because Sinvs is 0x0, but should actually
		  // be the identity matrix.
		  spatialPrecisions(k) = initialFwdPrior->GetPrecisions()(k,k);
		  assert( SP( initialFwdPrior->GetPrecisions(), IdentityMatrix(Nparams)-1 ).MaximumAbsoluteValue() == 0);

		  // Don't worry, this is multiplied by initialFwdPrior later
		  weightedMeans(k) = 0;
		  continue;
		}
	      else if (spatialPriorsTypes[k-1] == 'I')
		{
		  // get means from image prior MVN
		  priorMeans(k) = ImagePrior[k-1](v);

		  // precisions in same way as 'N' prior (for time being!)
		  spatialPrecisions(k) = initialFwdPrior->GetPrecisions()(k,k);
		  assert( SP( initialFwdPrior->GetPrecisions(), IdentityMatrix(Nparams)-1 ).MaximumAbsoluteValue() == 0);

		  weightedMeans(k) = 0;
		  continue;
		}

	      spatialPrecisions(k) = Sinvs[k-1](v,v);
	  //	  double testWeights = 0;	      
	      weightedMeans(k) = 0;
	      for (int n = 1; n <= Nvoxels; n++)
		if (n != v)
		  {
		    weightedMeans(k) += Sinvs[k-1](n,v) * 
		      (fwdPosteriorVox[n-1].means(k) - initialFwdPrior->means(k));
		    //	      testWeights += Cinvs[k-1](n,v);
		  }
	      //	  LOG_ERR("Parameter " << k << ", testWeights == " << testWeights << ", spatialPrecisions(k) == " << spatialPrecisions(k) << ", delta(k) == " << delta(k) << ", test2 == " << test2 << endl);
	    }
	  //      LOG_ERR("--------- end of voxel " << v << endl);
	  
	  assert(initialFwdPrior->GetPrecisions().Nrows() == spatialPrecisions.Nrows());
	  // Should check that earlier!  It's possible for basis=1 and priors=2x2 to slip through.  TODO
	  
	  // Should check that initialFwdPrior was already diagonal --
	  // this will cause real problems if it isn't!!
	  // (Safe way: SP of the covariance matrices -- that'd force 
	  // diagonality while preserving individual variance.)
	 
//cout << "Spatial precisions: " << spatialPrecisions;

	  DiagonalMatrix finalPrecisions = spatialPrecisions;
	    //	    SP(initialFwdPrior->GetPrecisions(),spatialPrecisions);

//cout << "initialFwdPrior->GetPrecisions() == " << initialFwdPrior->GetPrecisions();
//cout << "initialFwdPrior->GetCovariance() == " << initialFwdPrior->GetCovariance();

	  ColumnVector finalMeans = priorMeans
	    - spatialPrecisions.i() * weightedMeans;

//cout << "Final means and precisions: " << finalMeans << finalPrecisions;

	  // Preserve the shrinkageType ones from before.
	  // They'd better be diagonal!
	  for (int k = 1; k <= Nparams; k++)
	    if (spatialPriorsTypes[k-1] == shrinkageType)
	      {
		finalPrecisions(k) = fwdPriorVox[v-1].GetPrecisions()(k,k);
		finalMeans(k) = fwdPriorVox[v-1].means(k);
	      }
	  
	  fwdPriorVox[v-1].SetPrecisions( finalPrecisions );	  
	  fwdPriorVox[v-1].means = finalMeans;
	  // Definitely a minus here.
	}
	
	
	if (needF)
	  { 
	    F = noise->CalcFreeEnergy( *noiseVox[v-1], *noiseVoxPrior[v-1], 
				       fwdPosteriorVox[v-1], fwdPriorVox[v-1],
				       linearVox[v-1], data.Column(v) );
	    F += Fard;
	  }

	if (printF) 
	  LOG << "      Fbefore == " << F << endl;
	
        // Produces heaps of output and not very useful for debugging:
	//        LOG << "Voxel " << v << " of " << Nvoxels << endl;
	
	noise->UpdateTheta( *noiseVox[v-1],  
			    fwdPosteriorVox[v-1], fwdPriorVox[v-1], 
			    linearVox[v-1], data.Column(v), 
			    fwdPosteriorWithoutPrior.at(v-1));	


	if (needF) 
	  {
	    F = noise->CalcFreeEnergy( *noiseVox[v-1], *noiseVoxPrior[v-1], 
				       fwdPosteriorVox[v-1], fwdPriorVox[v-1],
				       linearVox[v-1], data.Column(v) );
	    F += Fard;
	    // Fard does NOT change because we haven't updated fwdPriorVox yet.
	  }
	
	if (printF) 
	  LOG << "      Ftheta == " << F << endl;

	/* MOVED BELOW -- 2007-11-23
	if (!lockedLinearEnabled)
	  linearVox[v-1].ReCentre( fwdPosteriorVox[v-1].means );
	
	if (needF) 
	  F = noise->CalcFreeEnergy( *noiseVox[v-1], *noiseVoxPrior[v-1], 
				     fwdPosteriorVox[v-1], fwdPriorVox[v-1],
				     linearVox[v-1], data.Column(v) );
	if (printF) 
	  LOG << "      Flin == " << F << endl;
	*/

	
      } 
    // QUICK INTERRUPTION: Voxelwise calculations continue below.

    if (useSimultaneousEvidenceOptimization)
      {
	Tracer_Plus tr("useSimultaneousEvidenceOptimization calculations");

	Warning::IssueOnce("Using simultaneous evidence optimization");
	
	// Re-estimate fwdPriorVox for all voxels simultaneously, 
	// based on the full covariance matrix
	
	// Check it's the simple case (haven't coded up the correction 
	// factors yet)
	//	assert(initialFwdPrior->GetPrecisions() == IdentityMatrix(Nparams)); // now part of Sinvs
	//	assert(initialFwdPrior->means == -initialFwdPrior->means);  // but this still applies
	
	if (!(initialFwdPrior->means == -initialFwdPrior->means))
	  Warning::IssueAlways("Quick hack to avoid assertion with initialFwdPrior->means != 0");
	
	SymmetricMatrix SigmaInv(Nparams*Nvoxels);
	//	SymmetricMatrix Sigma(Nparams*Nvoxels);
	ColumnVector Mu(Nparams*Nvoxels);
	
	Tracer_Plus tr5("useSimultaneousEvidenceOptimization calculations -- first part");
	// These matrices consist of NxN matrices blocked together
	// so parameter k, voxel v is in row (or col): v + (k-1)*Nparams
	SymmetricMatrix Ci = -999*IdentityMatrix(Nvoxels*Nparams);
	SymmetricMatrix XXtr = -999*IdentityMatrix(Nvoxels*Nparams);
	ColumnVector XYtr(Nvoxels*Nparams); XYtr = -999;
	
	// Build Ci
	for (int k = 1; k <= Nparams; k++)
	  {
	    Ci.SymSubMatrix(Nvoxels*(k-1)+1, Nvoxels*k) = Sinvs[k-1];
	    // off-diagonal blocks are zero, by definition of the our priors
	    // (priors between parameters are independent)
	  }
	
	// Build XXtr and XYtr
	for (int v = 1; v <= Nvoxels; v++)
	  {
	    const SymmetricMatrix& tmp = fwdPosteriorWithoutPrior[v-1]->GetPrecisions();
	    ColumnVector tmp2 = tmp * (fwdPosteriorWithoutPrior[v-1]->means - initialFwdPrior->means);
	    for (int k1 = 1; k1 <= Nparams; k1++)
	      {
		XYtr(v+(k1-1)*Nvoxels) = tmp2(k1);
		for (int k2 = 1; k2 <= Nparams; k2++)
		  {
		    XXtr(v+(k1-1)*Nvoxels,v+(k2-1)*Nvoxels) = tmp(k1,k2);
		  }
	      }
	    // XXtr != 0 only where the row and col refer to the same voxel.
	  }
	
	    

	{ Tracer_Plus tr("useSimultaneousEvidenceOptimization calculations -- 1a");
	  SigmaInv = XXtr + Ci;
	}

	// OLD SLOW CODE
	//       	{ Tracer_Plus tr("useSimultaneousEvidenceOptimization calculations -- 1b");
	//	  Sigma = SigmaInv.i();
	//	}
	//	{ Tracer_Plus tr("useSimultaneousEvidenceOptimization calculations -- 1c");
	//	  Mu = Sigma * XYtr;
	//	}
	
	{ Tracer_Plus tr("useSimultaneousEvidenceOptimization calculation -- 1bc replacement");
	    Mu = SigmaInv.i() * XYtr;
	}
	



      
	for (int v = 1; v <= Nvoxels; v++)
	  {
	    Tracer_Plus tr("useSimultaneousEvidenceOptimization calculations -- second loop");
	    
	    ColumnVector muBefore = fwdPosteriorVox[v-1].means - initialFwdPrior->means;
	    
	    assert(firstParameterForFullEO == 1);

	    for (int k = 1; k <= Nparams; k++)
	      fwdPosteriorVox[v-1].means(k) = Mu(v+(k-1)*Nvoxels) + initialFwdPrior->means(k);

	    //	    if ((muBefore - fwdPosteriorVox[v-1].means).MaximumAbsoluteValue() > 1e-10)
	    //	      LOG << "mBef = " << muBefore.t() << "mAft = " << fwdPosteriorVox[v-1].means.t();

            if (useCovarianceMarginalsRatherThanPrecisions) {
		SymmetricMatrix Sigma = SigmaInv.i();

		SymmetricMatrix cov = fwdPosteriorVox[v-1].GetCovariance();
		SymmetricMatrix covOld = cov;

		Warning::IssueOnce("Full simultaneous diagonal thingy -- now in covariances!");

		for (int k1 = 1; k1 <= Nparams; k1++)
		    {
			for (int k2 = 1; k2 <= Nparams; k2++)
			    {
				cov(k1,k2) = Sigma(v+(k1-1)*Nvoxels,v+(k2-1)*Nvoxels);
			    }
		    }
		if ((cov-covOld).MaximumAbsoluteValue() > 1e-10)
		    LOG << "covBefore: " << covOld.AsColumn().t() << "covAfter: " << cov.AsColumn().t();
		fwdPosteriorVox[v-1].SetCovariance(cov);

	    } else {

	    SymmetricMatrix prec = fwdPosteriorVox[v-1].GetPrecisions();
	    
	    SymmetricMatrix precOld = prec;
	    //cout << "precBefore:\n" << prec;
	    
	    Warning::IssueOnce("Full simultaneous diagonal thingy");
	    //cout << "prec = \n" << prec << endl;
	    for (int k1 = 1; k1 <= Nparams; k1++)
	      {
		for (int k2 = 1; k2 <= Nparams; k2++)
		  {
		    prec(k1,k2) = SigmaInv(v+(k1-1)*Nvoxels,v+(k2-1)*Nvoxels);
		  }
	      }
	    
	    if ((prec-precOld).MaximumAbsoluteValue() > 1e-10)
	      LOG << "precBefore: " << precOld.AsColumn().t() << "precAfter: " << prec.AsColumn().t();
	    
	    //if ((fwdPosteriorVox[v-1].GetPrecisions() - prec).MaximumAbsoluteValue() > 1e-10)
	    //  cout << "pbef: " << fwdPosteriorVox[v-1].GetPrecisions() << "paft: " << prec;
	    
	    fwdPosteriorVox[v-1].SetPrecisions(prec);
	    }
	  }	
      }
    else if (useFullEvidenceOptimization)
      {
	Tracer_Plus tr("useFullEvidenceOptimization calculations");

//	assert(!useCovarianceMarginalsRatherThanPrecisions); 
	// Covariance marginals are broken below, and I think they're 
	// rubbish anyway


	Warning::IssueOnce("Using full evidence optimization; using " 
			   + string(useCovarianceMarginalsRatherThanPrecisions? "covariances." : "precisions."));
	
	// Re-estimate fwdPriorVox for all voxels simultaneously, 
	// based on the full covariance matrix
	
	// Check it's the simple case (haven't coded up the correction 
	// factors yet)
	//	assert(initialFwdPrior->GetPrecisions() == IdentityMatrix(Nparams)); // now part of Sinvs
	//	assert(initialFwdPrior->means == -initialFwdPrior->means);  // but this still applies
	
	vector<SymmetricMatrix> SigmaInv(Nparams);
	vector<SymmetricMatrix> Sigma(Nparams);
	vector<ColumnVector> Mu(Nparams);
	
	for (int k = 1; k <= Nparams; k++)
	  {
	    Tracer_Plus tr("useFullEvidenceOptimization calculations -- first loop");
	    const SymmetricMatrix& Ci = Sinvs[k-1];
	    SymmetricMatrix XXtr(Nvoxels);
	    ColumnVector XYtr(Nvoxels);

	    ColumnVector XXtrMuOthers(Nvoxels);
	    
	    // Initialize to junk values
	    XXtr = -999*IdentityMatrix(Nvoxels);
	    XYtr = -999;
	    
	    for (int v = 1; v <= Nvoxels; v++)
	      {
		const SymmetricMatrix& tmp = fwdPosteriorWithoutPrior[v-1]->GetPrecisions();
		XXtr(v,v) = tmp(k,k);
		
		ColumnVector tmp2 = tmp * (fwdPosteriorWithoutPrior[v-1]->means - initialFwdPrior->means); 
		XYtr(v) = tmp2(k);

		//		ColumnVector MuOthers = fwdPosteriorVox[v-1].means;
		ColumnVector MuOthers = fwdPosteriorVox[v-1].means - initialFwdPrior->means;
		MuOthers(k) = 0;
		ColumnVector tmp3 = tmp * MuOthers;
		XXtrMuOthers(v) = tmp3(k);

		Warning::IssueOnce("Corrected mistake in useFullEvidenceOptimization: initialFwdPrior->means (not k)");
		// Also notice the subtle difference above: MuOthers uses the actual posterior means, while XYtr uses
		// the priorless posterior means.
		// Also, the above XXtr do NOT include the correction for non-N(0,1) initialFwdPriors, because the 
		// Sinvs (and hence Ci etc) already include this correction.  This is different from the DerivEdDelta
		// calculation which uses Cs with 1 on the diagonal (originally to facilitate reuse across different rho
		// values, now just confusing).
	      }

	    //	    ColumnVector tmp4(Nvoxels); 
	    //	    tmp4 = initialFwdPrior->means(k);
	    //	    ColumnVector CiMu0 = Ci * tmp4;

	    { Tracer_Plus tr("useFullEvidenceOptimization calculations -- 1a");
	      SigmaInv.at(k-1) = XXtr + Ci;
	    }
	    { Tracer_Plus tr("useFullEvidenceOptimization calculations -- 1b");
	      Sigma.at(k-1) = SigmaInv[k-1].i();
	    }
	    { Tracer_Plus tr("useFullEvidenceOptimization calculations -- 1c");
	      //	      Mu.at(k-1) = Sigma[k-1] * (XYtr - XXtrMuOthers);
	      //	      Mu.at(k-1) = Sigma[k-1] * (XYtr - XXtrMuOthers + CiMu0);
	      Mu.at(k-1) = Sigma[k-1] * (XYtr - XXtrMuOthers);
	    }
	  }
	for (int v = 1; v <= Nvoxels; v++)
	  {
	    Tracer_Plus tr("useFullEvidenceOptimization calculations -- second loop");

	    ColumnVector muBefore = fwdPosteriorVox[v-1].means;

	    for (int k = firstParameterForFullEO; k <= Nparams; k++)
	      fwdPosteriorVox[v-1].means(k) = Mu[k-1](v) + initialFwdPrior->means(k);

	    //	    if ((muBefore - fwdPosteriorVox[v-1].means).MaximumAbsoluteValue() > 1e-10)
	    //	      {
	    //		LOG << "mBef = " << muBefore.t() << "mAft = " << fwdPosteriorVox[v-1].means.t();
	    //	      }

	    if (useCovarianceMarginalsRatherThanPrecisions)
	      {
		SymmetricMatrix cov = 
		  SP(fwdPosteriorVox[v-1].GetCovariance(),
		     IdentityMatrix(Nparams));
		Warning::IssueOnce("Covariance diagonal thingy");
		//cout << "cov = \n" << cov << endl;

		for (int k = firstParameterForFullEO; k <= Nparams; k++)
		  cov(k,k) = Sigma[k-1](v,v);

		fwdPosteriorVox[v-1].SetCovariance(cov);
	      }
	    else if (keepInterparameterCovariances)
	      {
		Warning::IssueOnce("Keeping inter-parameter covariances from VB!");
	      }
	    else
	      {
		SymmetricMatrix prec = 
		  SP(fwdPosteriorVox[v-1].GetPrecisions(),
		     IdentityMatrix(Nparams));

		SymmetricMatrix precOld = prec;
		//cout << "precBefore:\n" << prec;

		Warning::IssueOnce("Precision diagonal thingy");
		//cout << "prec = \n" << prec << endl;
		for (int k = firstParameterForFullEO; k <= Nparams; k++)
		  prec(k,k) = SigmaInv[k-1](v,v);

		if ((prec-precOld).MaximumAbsoluteValue() > 1e-10)
		  LOG << "precBefore: " << precOld.AsColumn().t() << "precAfter: " << prec.AsColumn().t();

		//if ((fwdPosteriorVox[v-1].GetPrecisions() - prec).MaximumAbsoluteValue() > 1e-10)
		//  cout << "pbef: " << fwdPosteriorVox[v-1].GetPrecisions() << "paft: " << prec;

		fwdPosteriorVox[v-1].SetPrecisions(prec);
	      }
	    assert(fwdPosteriorVox[v-1].GetSize() == Nparams);
	  }
      }


    // Back to your regularly-scheduled voxelwise calculations
    for (int v = 1; v <= Nvoxels; v++)
      {
	// some models may want extra information about the data
	if (suppdata.Ncols() > 0) {
	  model->pass_in_data(  data.Column(v) ,  suppdata.Column(v) );
	}
	else {
	  model->pass_in_data(  data.Column(v) );
	}
	model->pass_in_coords(coords.Column(v));

	double &F = resultFs.at(v-1);  // short name

	noise->UpdateNoise( *noiseVox[v-1], *noiseVoxPrior[v-1], 
        fwdPosteriorVox[v-1], linearVox[v-1], data.Column(v) );

	if (needF) 
	  F = noise->CalcFreeEnergy( *noiseVox[v-1], *noiseVoxPrior[v-1], 
				     fwdPosteriorVox[v-1], fwdPriorVox[v-1],
				     linearVox[v-1], data.Column(v) );
	if (printF) 
	  LOG << "      Fnoise == " << F << endl;

	//} catch (...) 
	//{keepGoing[v-1] = false; cout << "Bad Voxel! " << v << endl;} 

	//* MOVED HERE on Michael's advice -- 2007-11-23
	if (!lockedLinearEnabled)
	  linearVox[v-1].ReCentre( fwdPosteriorVox[v-1].means );
	
	if (needF) 
	  F = noise->CalcFreeEnergy( *noiseVox[v-1], *noiseVoxPrior[v-1], 
				     fwdPosteriorVox[v-1], fwdPriorVox[v-1],
				     linearVox[v-1], data.Column(v) );
	if (printF) 
	  LOG << "      Flin == " << F << endl;
	// */
	
      }


    
    
    
    /*  
    if (delta(7) == .5)
      {
	delta(7) = 1e12;
	LOG_ERR("Testing Hack: Just before delta updates, set delta(7) = " << delta(7) << endl);
      }
    */

    // Moved shrinkage updates to the beginning!!

    isFirstIteration = false;
    
    // next iteration:
  } while (!conv->Test( globalF ));
  
  // Phew!

  
  // Interesting addition: calculate "coefficient resels" from Penny et al. 2005
  for (int k = 1; k <= Nparams; k++)
    {
      ColumnVector gamma_vk(Nvoxels); // might be handy
      ColumnVector gamma_vk_eo(Nvoxels); // slightly different calculation (differs if using EO)
      gamma_vk_eo = -999;
      for (int v = 1; v <= Nvoxels; v++)
	{
	  gamma_vk(v) = 1 - fwdPosteriorVox[v-1].GetCovariance()(k,k) / fwdPriorVox[v-1].GetCovariance()(k,k);
	  if (fwdPosteriorWithoutPrior.at(v-1) != NULL)
	    {
	      gamma_vk_eo(v) = fwdPosteriorVox[v-1].GetCovariance()(k,k) 
		/ fwdPosteriorWithoutPrior[v-1]->GetCovariance()(k,k);
	    }
	}
      LOG_ERR("Coefficient resels per voxel for param " << k << ": " << gamma_vk.Sum()/Nvoxels
	      << " (vb) or " << gamma_vk_eo.Sum()/Nvoxels << " (eo)\n");
    }

  //  if (spatialPriorOutputCorrection)
  //    {
  //      Tracer_Plus tr("SpatialVariationalBayes::DoCalculations - Spatial Prior Ouput Correction");
  //
  //      // Instead of using the diagonal of the precision matrix as the prior
  //      // precision, use 1/ the diagonal of the covariance matrix!
  //      // All elements on this diagonal are exactly exp(-rho(k)).
  //      
  //      DiagonalMatrix spatialCovariance(Nparams);
  //      for (int k = 1; k <= Nparams; k++) 
  //	{
  //	  spatialCovariance(k) = exp(-rho(k));
  //	}
  //
  //      for (int v = 1; v <= Nvoxels; v++)
  //	{
  //	  assert(
  //		 fwdPosteriorVox[v-1].GetPrecisions() 
  //		 == fwdPosteriorWithoutPrior.at(v-1).GetPrecisions() 
  //	 + fwdPriorVox[v-1].GetPrecisions());
  //
  //  fwdPosteriorVox[v-1].SetPrecisions(
  //      fwdPosteriorWithoutPrior.at(v-1).GetPrecisions() 
  //      + SP(spatialCovariance.i(), initialFwdPrior->GetPrecisions()) );
  //
  //  //	  cout << "fwdPriorVox[v-1].GetPrecisions:\n" 
  //  //	       << fwdPriorVox[v-1].GetPrecisions()
  //  //	       << "spatialCovariance.i(): \n"
  //  //	       << spatialCovariance.i();
  //
  //  // Leave the mean unchanged
  //}
  //}


  for (int v = 1; v <= Nvoxels; v++)
    {
      resultMVNs[v-1] = new MVNDist(
        fwdPosteriorVox[v-1], noiseVox[v-1]->OutputAsMVN() );

      if (alsoSaveWithoutPrior)
	{
	  resultMVNsWithoutPrior.at(v-1) = new MVNDist(
	    *fwdPosteriorWithoutPrior[v-1], noiseVox[v-1]->OutputAsMVN() );
	  // Should probably save the noiseWithoutPriors, but don't need that yet (ever?)
	}
    }
 
  // resultFs are already stored as we go along.

  if (!needF)
    {
      for (int v = 1; v <= Nvoxels; v++)
	assert(resultFs.at(v-1) == 9999); 
      // check we're not throwing away anything useful

      resultFs.clear();
      // clearing resultFs here should prevent an F image from being saved.
    }

  // Save Sinvs if possible
  if (alsoSaveSpatialPriors) 
    {
#ifdef __FABBER_LIBRARYONLY
	throw Logic_error("Not implemented for fabber_library");
#else
      // Copied from MVNDist::Save.  There are enough subtle differences 
      // to justify duplicating the code here.

      Tracer_Plus tr("Saving Sinvs");
      Matrix vols;
      
      vols.ReSize(Nparams, Nvoxels*Nvoxels);
      for (int k = 1; k <= Nparams; k++)
	{
	  assert(Sinvs.at(k-1).Nrows() == Nvoxels);
	  Matrix full = Sinvs.at(k-1); // easier to visualize if in full form
	  vols.Row(k) = full.AsColumn().t();
	}
      
      volume4D<float> output(Nvoxels,Nvoxels,1,Nparams);
      output.set_intent(NIFTI_INTENT_SYMMATRIX,1,1,1);
      output.setdims(1,1,1,1);
      output.setmatrix(vols);
      // no unThresholding needed
      cout << vols.Nrows() << "," << vols.Ncols() << endl;
      
      save_volume4D(output,outputDir + "/finalSpatialPriors");
#endif
    }
      


  // Delete stuff (avoid memory leaks)
  for (int v = 1; v <= Nvoxels; v++)
    {
      delete noiseVox[v-1];
      delete noiseVoxPrior[v-1];
      delete fwdPosteriorWithoutPrior.at(v-1);
    }
}

// Binary search for data(index) == num
// Assumes data is sorted ascending!!
// Either returns an index such that data(index) == num
//   or -1 if num is not present in data. 
inline int binarySearch(const ColumnVector& data, int num)
{
  int first = 1, last = data.Nrows();

  while (first <= last)
    {
      int test = (first+last)/2;

      if (data(test) < num)
	first = test + 1;
      else if (data(test) > num)
	last = test - 1;
      else if (data(test) == num)
	return test;
      else
	assert(false); // logic error!  data wasn't sorted?
    }
  return -1;
}

bool IsCoordMatrixCorrectlyOrdered(const Matrix& voxelCoords)
{
    assert(voxelCoords.Nrows()==3);
    const int nVoxels = voxelCoords.Ncols();
    for (int v = 1; v <= nVoxels-1; v++)
    {
	ColumnVector diff = voxelCoords.Column(v+1) - voxelCoords.Column(v);
	int d = sign(diff(1)) + 10*sign(diff(2)) + 100*sign(diff(3)); // +1 = +x, +10 = +y, +100 = +z, -99 = -z+x, etc.
        if (d<=0) 
	{
	    LOG << "Found mis-ordered voxels " << v << " and " << v+1 << ": d=" << d << endl;
	    return false;
	}
    }
    return true;
}

void SpatialVariationalBayes::CalcNeighbours(const Matrix& voxelCoords)
{
    Tracer_Plus tr("SpatialVariationalBayes::CalcNeighbours from voxelCoords");
    // NOTE there's a bit of an incompatibility here: CalcDistances assumes voxelCoords are in mm, while this assumes 
    // that they're integers!

    const int nVoxels = voxelCoords.Ncols();
    if (nVoxels < 2) return;
    if (!IsCoordMatrixCorrectlyOrdered(voxelCoords))
	throw Invalid_option("Coordinate matrix must be in correct order to use adjacency-based priors.");

#ifdef __FABBER_LIBRARYONLY
    throw Logic_error("CalcNeighbours currently requires NEWIMAGE to be compiled in.");
#else
    volume<float> mask(voxelCoords.Row(1).Maximum()+1, voxelCoords.Row(2).Maximum()+1, voxelCoords.Row(3).Maximum()+1); 
    mask = 0.0; 
    for (int v=1; v<=nVoxels; v++)
    {
	float& maskVox = mask( voxelCoords(1,v), voxelCoords(2,v), voxelCoords(3,v) );
	assert(maskVox == 0);
	maskVox = v;
    }

    // Essential thing: threshold the mask!
    mask.binarise(1e-16,mask.max()+1,exclusive);

    // Okay, now we've built up a close-enough mask matrix, pass the buck to the old version:
    CalcNeighbours(mask);
#endif // __FABBER_LIBRARYONLY
}

#ifndef __FABBER_LIBRARYONLY
void SpatialVariationalBayes::CalcNeighbours(const volume<float>& mask)
{
  Tracer_Plus tr("SpatialVariationalBayes::CalcNeighbours");

  ColumnVector preThresh((int)mask.sum());
  const int nVoxels = preThresh.Nrows();

  assert(nVoxels > 0); // Probably because mask is empty, logic error

  int offset(0);
  int count(1);
  for(int z=0;z<mask.zsize();z++)
    for(int y=0;y<mask.ysize();y++)     
      for(int x=0;x<mask.xsize();x++)
      {
	if (mask(x,y,z)!=0) 
	  preThresh(count++)=offset;
	offset++;
      }

  vector<int> delta;
  if (spatialDims >= 1)
  { 
    delta.push_back(1); // next row
    delta.push_back(-1); // prev row
  }
  if (spatialDims >= 2)
  {
    delta.push_back(mask.xsize()); // next column
    delta.push_back(-mask.xsize()); // prev column
  }
  if (spatialDims >= 3)
  {
    delta.push_back(mask.xsize()*mask.ysize()); // next slice
    delta.push_back(-mask.xsize()*mask.ysize()); // prev slice
  }

  /* Tedious output and I'm pretty sure this is working correctly.  
  LOG << "Using offsets of " << delta 
      << "\nNote these have not been checked!\n"
      << "Also note that they may wrap around if both edges are in mask\n";
  */
  
  neighbours.resize(nVoxels);
  
  for(int vid = 1; vid <= nVoxels; vid++) // voxel id
    {
      int pos = int(preThresh(vid));

      for (unsigned n = 0; n < delta.size(); n++) // neighbour-search number
	{
	  
	  // is there a voxel at this neighbour position?
	  int id = binarySearch(preThresh, pos + delta[n]);
	  // indexed from 1; id == -1 if not found.
	  
	  // No such voxel: continue
	  if (id < 0)
	    continue;

	  // Check for wrap-around
	  if ( (delta.size() > n+2) &&
	       (
		( (delta[n]>0) && (pos %  delta[n+2] == 0) ) ||
		( (delta[n]<0) && (pos % -delta[n+2] == 1) ) ) 
	       )
	    {
	      /* Tedious, and I never looked at it anyway.
	      LOG_ERR("Note: avoiding wrap-around from position "
		      << pos << " to " << pos + delta[n] << endl);
	      */
	      continue;
	    }

	  // If we get this far, add it to the list
	  neighbours.at(vid-1).push_back(id);
	}
    }
  
  // Neighbours-of-neighbours, excluding self, and duplicated if there 
  // are two routes to get there (diagonally connected)
  neighbours2.resize(nVoxels);
  for(int vid = 1; vid <= nVoxels; vid++)
    {
      for (unsigned n1 = 0; n1 < neighbours.at(vid-1).size(); n1++)
	{
	  int n1id = neighbours[vid-1].at(n1);
	  int checkNofN = 0;
	  for (unsigned n2 = 0; n2 < neighbours.at(n1id-1).size(); n2++)
	    {
	      int n2id = neighbours[n1id-1].at(n2);
	      if (n2id != vid)
		neighbours2[vid-1].push_back(n2id);
	      else
		checkNofN++;
	    }
	  assert(checkNofN == 1);
	  // Each of this voxel's neighbours must have this voxel 
	  // as a neighbour.
	}
    }    

  /* Tedious, and I never looked at it anyway.
  LOG << "Neighbours are (for dims==" << spatialDims << ")\n";
  for (int v = 1; v <= nVoxels; v++)
    LOG << v << ": " << neighbours.at(v-1) 
	<< "-" << neighbours2.at(v-1) << endl;
  */
}
#endif //__FABBER_LIBRARYONLY

#if defined(__FABBER_LIBRARYONLY_TESTWITHNEWIMAGE) || !defined(__FABBER_LIBRARYONLY)
// Helper function, also used in fabber_library's test main()
void ConvertMaskToVoxelCoordinates(const volume<float>& mask, Matrix& voxelCoords)
{
    Tracer_Plus("ConvertMaskToVoxelCoordinates");

    ColumnVector preThresh((int)mask.sum()); // mask has previously been binarized to 0 or 1
    const int nVoxels = preThresh.Nrows();

    int offset(0);
    int count(1);
    for(int z=0;z<mask.zsize();z++)
      for(int y=0;y<mask.ysize();y++)     
	for(int x=0;x<mask.xsize();x++)
	{
	  if (mask(x,y,z)!=0) 
	    preThresh(count++)=offset;
	  offset++;
	}
 
    const int dims[3] = {mask.xsize(),mask.ysize(),mask.zsize()};

    assert(mask.xsize()*mask.ysize()*mask.zsize() > 0);
    assert(mask.xdim()*mask.ydim()*mask.zdim() > 0); // no zeroes!
    
    LOG_SAFE_ELSE_DISCARD( "Calculating distance matrix, using voxel dimensions: " 
	<< mask.xdim() << " by " << mask.ydim() << " by " << mask.zdim() << " mm\n" );
    
    ColumnVector positions[3];  // indices
    positions[0].ReSize(nVoxels);
    positions[1].ReSize(nVoxels);
    positions[2].ReSize(nVoxels);
    
    for (int vox = 1; vox <= nVoxels; vox++)
    {
        int pos = (int)preThresh(vox) - 1; // preThresh appears to be 1-indexed.
	assert(pos>=0);
        positions[0](vox) = pos % dims[0];
        pos = pos / dims[0];
        positions[1](vox) = pos % dims[1];
        pos = pos / dims[1];
        positions[2](vox) = pos % dims[2];
        pos = pos / dims[2];
	//cout << vox << ' ' << (int)preThresh(vox)<< ' ' << dims[0] << ' ' << dims[1] << ' ' << dims[2] << ' ' << pos << endl;
        assert(pos == 0); // fails if preThresh >= product of dims[0,1,2]
        
        //double pos = preThresh(vox);
        //positions[2](vox) = floor(pos/dims[0]/dims[1]);
        //pos -= positions[2](vox)*dims[0]*dims[1];
        //positions[1](vox) = floor(pos/dims[0]);
        //pos -= positions[1](vox)*dims[0];               
        //positions[0](vox) = pos;
        //assert(positions[2](vox) < dims[2]);
        //assert(positions[1](vox) < dims[1]);
        //assert(positions[0](vox) < dims[0]);
        assert(preThresh(vox)-1 == 
              positions[0](vox) 
            + positions[1](vox)*dims[0]
            + positions[2](vox)*dims[0]*dims[1] );

    }    
    
    // Pass on to the alternative form
    voxelCoords = ( positions[0] | positions[1] | positions[2] ).t(); // 3 x Nvox
}
#endif //__FABBER_LIBRARYONLY_TESTWITHNEWIMAGE || !__FABBER_LIBRARYONLY

#ifndef __FABBER_LIBRARYONLY
void CovarianceCache::CalcDistances(const volume<float>& mask, const string& distanceMeasure)
{
    Tracer_Plus tr("CovarianceCache::CalcDistances mask -> voxelCoords");

    Matrix voxelCoords;
    ConvertMaskToVoxelCoordinates(mask, voxelCoords);
    Matrix voxelCoordsInMm = voxelCoords;
    voxelCoordsInMm.Row(1) *= mask.xdim();
    voxelCoordsInMm.Row(2) *= mask.ydim();
    voxelCoordsInMm.Row(3) *= mask.zdim();
    CalcDistances(voxelCoords, distanceMeasure);
}
#endif //__FABBER_LIBRARYONLY

// Note: voxelCoords should really be in MM, not indices; only really matters if it's aniostropic or you're using the
// smoothness values directly.
void CovarianceCache::CalcDistances(const NEWMAT::Matrix& voxelCoords, const string& distanceMeasure)
{
    // Back to the original form
    ColumnVector positions[3];
    positions[0] = voxelCoords.Row(1).t();
    positions[1] = voxelCoords.Row(2).t();
    positions[2] = voxelCoords.Row(3).t();
    const int nVoxels = positions[0].Nrows();
    const double dimSize[3] = {1.0, 1.0, 1.0};  // dimSize is already included in voxelCoords
    
    if (nVoxels > 7500)
      LOG_SAFE_ELSE_CERR("WARNING: Over " << int(2.5*nVoxels*nVoxels*8/1e9) 
	      << " GB of memory will be used just to calculate "
	      << "the distance matrix.  Hope you're not trying to invert "
	      << "this sucker!\n" << endl);

    SymmetricMatrix relativePos[3];  // millimetres
    ColumnVector allOnes(nVoxels); allOnes = 1.0;
    
    for (int dim = 0; dim < 3; dim++)
    {
        Matrix rel = dimSize[dim] * 
            (positions[dim]*allOnes.t() - allOnes*positions[dim].t());
        assert(rel == -rel.t());
        assert(rel.Nrows() == nVoxels);
        relativePos[dim] << rel;  
            // down-convert to symmetric matrix (lower triangle so all positive??) 
    }

    distances.ReSize(nVoxels);
  
    if (distanceMeasure == "dist1") // absolute Euclidean distance
      {
	LOG_ERR("Using absolute Euclidean distance\n");
        for (int a = 1; a <= nVoxels; a++)
	  {
            for (int b = 1; b <= a; b++)
	      {
                distances(a,b) = sqrt(
                    relativePos[0](a,b)*relativePos[0](a,b) +
                    relativePos[1](a,b)*relativePos[1](a,b) +
                    relativePos[2](a,b)*relativePos[2](a,b) );
	      }
	  }
      }
    else if (distanceMeasure == "dist2") // Euclidian distance squared
      {
	LOG_ERR("Using almost-squared (^1.99) Euclidean distance\n");
        for (int a = 1; a <= nVoxels; a++)
	  {
            for (int b = 1; b <= a; b++)
	      {
                distances(a,b) = pow(
                    relativePos[0](a,b)*relativePos[0](a,b) +
                    relativePos[1](a,b)*relativePos[1](a,b) +
                    relativePos[2](a,b)*relativePos[2](a,b), 0.995 );
	      }
	  }	
      }
    else if (distanceMeasure == "mdist") // Manhattan distance (bad?)
      {
	LOG_ERR("Using Manhattan distance\n");
	LOG_ERR("WARNING: Seems to result in numerical problems down the line (not sure why)\n");
        for (int a = 1; a <= nVoxels; a++)
	  {
            for (int b = 1; b <= a; b++)
	      {
                distances(a,b) = 
		  fabs(relativePos[0](a,b)) +
		  fabs(relativePos[1](a,b)) +
                  fabs(relativePos[2](a,b));
	      }
	  }	
      }   
    else
      {
        throw Invalid_option("\nUnrecognized distance measure: " + distanceMeasure + "\n");
      }
    
}

#include "tools.h"

// /*
class DerivFdRho : public GenericFunction1D
{
public:
  virtual double Calculate(double input) const;
  DerivFdRho(const CovarianceCache& cov, const DiagonalMatrix& cr,
	     const ColumnVector& mdr, const double d)
    : covar(cov), covRatio(cr), meanDiffRatio(mdr), delta(d)
      { return; }
  
  // don't need a PickFasterGuess

private:
  const CovarianceCache& covar;
  const DiagonalMatrix& covRatio;
  const ColumnVector& meanDiffRatio;
  const double delta;
};

double DerivFdRho::Calculate(const double rho) const
{
  Tracer_Plus tr("DerivFdRho::Calculate");
  
  const int Nvoxels = covar.GetDistances().Nrows();
  const SymmetricMatrix& Cinv = covar.GetCinv(delta);  
  
  double out = 0;
  out += 0.5 * Nvoxels;
  out += -0.5 * (covRatio * exp(rho)*Cinv).Trace();
  out += -0.5 * (meanDiffRatio.t() * exp(rho)*Cinv * meanDiffRatio).AsScalar();

  /* Old version (pre Oct 10) -- actually gives almost-identical results (just the prior).

  double out = Nvoxels; // == Trace(S*Sinv)
  //  out -= exp(rho) * SP(covRatio, Cinv).Trace();
  out -= exp(rho) * (covRatioSupplemented * Cinv).Trace();
  out -= exp(rho) * (meanDiffRatio.t() * Cinv * meanDiffRatio).AsScalar();
  
  // Want a scale-free prior on exp(rho)

  //  LOG_ERR("Scale-free prior on rho... +exp(-2*rho)\n");
  out -= -1/exp(2*rho); // correct?

  */
  
  return out;
}
// */

// Evidence optimization
class DerivEdDelta : public GenericFunction1D
{
public:
  virtual double Calculate(double delta) const;
  DerivEdDelta(const CovarianceCache& c, 
	       const vector<MVNDist*>& fpwp, 
  //	       const vector<SymmetricMatrix>& Si)
	       const int kindex,
	       const MVNDist* initFwdPrior,
	       const bool r = false)
    : covar(c), fwdPosteriorWithoutPrior(fpwp),
      k(kindex), initialFwdPrior(initFwdPrior), allowRhoToVary(r)
      //Sinvs(Si) 
  { return; }

  //  virtual bool PickFasterGuess(double* guess, double lower, double upper, bool allowEndpoints = false) const
  virtual double OptimizeRho(double delta) const;

private:

  const CovarianceCache& covar; // stores C, Cinv, distance matrix?, etc.
  //  const SymmetricMatrix& XXtr; // = X*X'*precision -- replaces covRatio
  //  const ColumnVector& XYtr; // = X*Y'*precision -- replaces meanDiffRatio;
  const vector<MVNDist*>& fwdPosteriorWithoutPrior;

  // To allow multiple parameters, will need storage for k (which param to optimize) 
  // and the other prior precision matrices (may be from Penny, Nonspatial, etc.)

  const int k;
  const MVNDist* initialFwdPrior;

  //const vector<SymmetricMatrix>& Sinvs;

  bool allowRhoToVary;
};

double DerivEdDelta::OptimizeRho(double delta) const
{
  Tracer_Plus tr("DerivEdDelta::OptimizeRho");

  double rho;
  if (!allowRhoToVary)
    return 0.0;

  // This is just copy-pasted from ::Calculate.  There are more efficient 
  // ways to do this!
  const SymmetricMatrix& dist = covar.GetDistances();
  const int Nvoxels = dist.Nrows();

  assert(initialFwdPrior->GetCovariance()(k,k) == 1); // unimplemented correction factor!

  DiagonalMatrix XXtr(Nvoxels);
  ColumnVector XYtr(Nvoxels);
  { 
    Tracer_Plus tr("Populating XXtr and XYtr");
    assert(Nvoxels == (int)fwdPosteriorWithoutPrior.size());
    for (int v = 1; v <= Nvoxels; v++)
      {
	XXtr(v,v) = fwdPosteriorWithoutPrior[v-1]->GetPrecisions()(k,k);
	XYtr(v) = XXtr(v,v) * (fwdPosteriorWithoutPrior[v-1]->means(k) - initialFwdPrior->means(k));
      }
    assert(XXtr.Nrows() == Nvoxels);
    assert(XYtr.Nrows() == Nvoxels);
  }

  SymmetricMatrix Sigma;
  {
    Tracer_Plus tr("Just calculating Sigma");
    Sigma = (XXtr + covar.GetCinv(delta)).i();
  }

  const ColumnVector mu = Sigma * XYtr; 
 
  rho = -log(1.0/Nvoxels * ( (Sigma+mu*mu.t())*covar.GetCinv(delta) ).Trace() );

  LOG_ERR("rho == " << rho);

  return rho;
}

double DerivEdDelta::Calculate(double delta) const
{
  Tracer_Plus tr("DerivEdDelta::Calculate");

//  assert(delta >= 0.05); // Will be slow below this scale

  const SymmetricMatrix& dist = covar.GetDistances();
  const int Nvoxels = dist.Nrows();

  DiagonalMatrix XXtr(Nvoxels);
  ColumnVector XYtr(Nvoxels);
  
  { 
    Tracer_Plus tr("Populating XXtr and XYtr");
    assert(Nvoxels == (int)fwdPosteriorWithoutPrior.size());
    for (int v = 1; v <= Nvoxels; v++)
      {
	XXtr(v,v) = fwdPosteriorWithoutPrior[v-1]->GetPrecisions()(k,k) * initialFwdPrior->GetCovariance()(k,k);
	XYtr(v) = XXtr(v,v) * (fwdPosteriorWithoutPrior[v-1]->means(k) - initialFwdPrior->means(k))
	  * sqrt(initialFwdPrior->GetPrecisions()(k,k));
	Warning::IssueOnce("Using the new XYtr correction (*sqrt(precision))");
      }
    assert(XXtr.Nrows() == Nvoxels);
    assert(XYtr.Nrows() == Nvoxels);
  }
  
  Tracer_Plus tr2("Calculating Sigma etc...");

  double out; // Assigned by the following statement:
  const SymmetricMatrix& CiCodistCi = covar.GetCiCodistCi(delta, &out);
  SymmetricMatrix Sigma;
  {
    Tracer_Plus tr("Just calculating Sigma");
    Sigma = (XXtr + covar.GetCinv(delta)).i();
  }

  out -= (Sigma * CiCodistCi).Trace();

    // If the trace turns out to be slow, then
    // use identity: trace(a*b) == sum(sum(a.*b'))

  const ColumnVector mu = Sigma * XYtr;

  out -= (mu.t() * CiCodistCi * mu).AsScalar();
  out /= -4*delta*delta; // = -1/2 * d(1/delta)/ddelta.  Note Sahani used d(1/delta^2)/ddelta.

    if (0)
      {

// WORKS! This is a very strong prior that attracts delta towards 5.
//	LOG_ERR("Using a Ga(.0001,50000) prior on delta!\n");
//	const double c = 50000, b = .0001;
//	// DERIVATIVE OF -gammaln(c) + (c-1)*log(delta) - c*log(b) - delta/b;
//	out += (c-1)/delta - 1/b;
// Note: Matlab is gampdf(x, c, b) NOT b,c!

	// Uninformative prior should have b >>> 1.
	const double b = 1e40;
	const double c = 1/b;
	// Mean = b*c
	// Variance = b^2 * c

	// Strong prior:
	//	const double m = 3.0/4.0; // correcting for 1/4 scaling of this ASL data
	//	const double s = 1.0/4.0;
	//
	//	const double b = (s*s)/m;
	//	const double c = (m*m)/(s*s);

	Warning::IssueOnce("Using a Ga(" + stringify(b) + ", " 
			   + stringify(c) + " prior on delta!");
	// DERIVATIVE OF -gammaln(c) + (c-1)*log(delta) - c*log(b) - delta/b;
	out += (c-1)/delta - 1/b;

      }
    else
      {
	//	LOG_ERR("Not using any prior at all on delta\n");
	Warning::IssueOnce("Not using any prior at all on delta");
      }

  return out;
}

// helper function
void BlockInverse(const SymmetricMatrix& in, SymmetricMatrix& out);


// Old Free energy optimization
class DerivFdDelta : public GenericFunction1D
{
public:
    virtual double Calculate(double input) const;
    DerivFdDelta(const CovarianceCache& cov, const DiagonalMatrix& cr,
		 const ColumnVector& mdr, bool rv = true) 
      : covar(cov), covRatio(cr), meanDiffRatio(mdr), allowRhoToVary(rv)
        { return; } 

  virtual bool PickFasterGuess(double* guess, double lower, double upper, bool allowEndpoints = false) const
  { return covar.GetCachedInRange(guess, lower, upper, allowEndpoints); }
    
  virtual double OptimizeRho(double delta) const
    {
      if (!allowRhoToVary)
	return 0.0;

      // Use the following approximation, without the prior.
      // The difference due to prior is typically very small..
      // >5: identical to 6 digits
      // 1.80468 -> 1.80471
      // -7 -> -4
      // -12 -> -6
      // -14 -> -7
      // -19 -> -8
      // In my data set, these correspond to delta = 1e12 or 1e15
      // For values with a rho (not dt), typically <1000 and highest
      // observed stop value was 700,000.

      const int Nvoxels = covar.GetDistances().Nrows();
      const SymmetricMatrix& Cinv = covar.GetCinv(delta);
      //      const double tmp = SP(covRatio, Cinv).Trace() 
      const double tmp = (covRatio * Cinv).Trace() 
        + (meanDiffRatio.t() * Cinv * meanDiffRatio).AsScalar();
      // Note: tmp can be negative if there's a numerical problem.
      // this means rho2 = NaN so it'll go on to the search method,
      // which should deal with this case reasonably well...

      //      cout << Cinv.Trace()<<endl
      //	   << SP(Cinv, IdentityMatrix(Nvoxels)) << endl
	//<< covRatio 
	// << meanDiffRatio.t() 
      //	   << SP(covRatio, Cinv).Trace() << endl
      //	   << (meanDiffRatio.t() * Cinv * meanDiffRatio).AsScalar()
      //	   << "  tmp == " << tmp << endl;
      
      const double rho2 = -log(tmp/Nvoxels);
      LOG_ERR("  rho2 == " << rho2 << "\n");

      double rho;
      if (true) //(rho2 > 2.0)
	{
	  rho = rho2;
	}
      else
	{
	  // Rho is small enough that the prior might actually 
	  // make a difference -- so calculate it the more-accurate-but-slow
	  // way.

	  //	  LOG_ERR("\n--- OPTIMIZING FOR RHO ---\n");
	  DerivFdRho fcn2(covar, covRatio, meanDiffRatio, delta);
	  BisectionGuesstimator guesser;
	  rho = DescendingZeroFinder(fcn2)
	    .InitialGuess(1).TolY(0.0001).RatioTolX(1.001).Verbosity(0)
	    .SetGuesstimator(&guesser).SearchMin(-70).SearchMax(70);
	  // Observed range: -35 to +40
	  // If the inversion causes numerical problems, then the above
	  // we get tmp < 0, so rho2 = NaN, so we search, and typically end
	  // up stuck at the upper limit.  700 seems to cause problems, but it
	  // seems to recover if we only go as high as 70.

	  //	  LOG_ERR("--- OPTIMIZED, RHO == " << rho << " ---\n\n");
	  LOG_ERR(" rho == " << rho << endl);
	}

      return rho;
    }

private:
    const CovarianceCache& covar;
    const DiagonalMatrix& covRatio;
  //const SymmetricMatrix& covRatioSupplemented;
    const ColumnVector& meanDiffRatio;
    bool allowRhoToVary;
};

double DerivFdDelta::Calculate(const double delta) const
{
  Tracer_Plus tr("DerivFdDelta::Calculate");

  const double rho = OptimizeRho(delta);
  // Returns rho = 0 if !allowRhoToVary.
  // Otherwise we return the value of delta optimized over all rhos!

  // Now that we've optimized over rho, get on with the calculation:

    //    LOG_ERR("        Calculating dF/ddelta at " << delta << endl);
    
    assert(delta >= 0.05);
    //    const SymmetricMatrix& dist = covar.GetDistances();
    const SymmetricMatrix& dist = covar.GetDistances();
    const int Nvoxels = dist.Nrows();
    assert(covRatio.Nrows() == Nvoxels);
    assert(meanDiffRatio.Nrows() == Nvoxels);
    
    //SymmetricMatrix C = covar.GetC(delta);
    //const SymmetricMatrix& Ci = covar.GetCinv(delta);
    //const Matrix& CiCodist = covar.GetCiCodist(delta);
    //double out = covar.GetCiCodist(delta).Trace();

    double out; 
    const SymmetricMatrix& CiCodistCi = covar.GetCiCodistCi(delta, &out);
    // Above does: out = trace(CiCodist)

    //    cout << "The uncacheable parts... " << flush;
    //LOG_ERR("values: " << out);

    // Correct???, but slower:
    //    out -= exp(rho) * (
    //		         (covRatio.i() - diag(old Ci) + Ci).i() * CiCodistCi
    //		      ).Trace();

    // Hopefully also correct (after iterations) but faster:
// METHOD USED BEFORE 2008-03-13
    out -= exp(rho) * (
		       covRatio * CiCodistCi
		       ).Trace();

    // If the trace turns out to be slow, then
    // use identity: trace(a*b) == sum(sum(a.*b'))

    out -= exp(rho) * 
      (meanDiffRatio.t() * CiCodistCi * meanDiffRatio).AsScalar() ; // REQUIRES MEANDIFFRATIO
    out /= -4*delta*delta;
    //    cout << "done." << endl;

    // Slower (equivalent) method
    //C = exp(-0.5*dist/delta_s);
    //Ci = inv(C);
    //Codist = C.*dist;
    //out2 = -0.25/delta_s^2 * ( ...
    //    trace(Ci * Codist) ...
    //    + trace(inv(diag(precRatio)) * Ci * Codist * Ci) ...
    //    + meanDiffRatio' * Ci * Codist * Ci * meanDiffRatio );
    //[out out2];
    //assert(diff(ans)/mean(ans) < 1e-10)
    
    //    LOG_ERR("        dF/ddelta at\t" << delta << " =\t" << out << endl);

    //    LOG_ERR(" - " << SP(covRatio, CiCodistCi).Trace() 
    //    << " - " << (meanDiffRatio.t() * CiCodistCi * meanDiffRatio).AsScalar() 
    //    << " all / " <<  -4*delta*delta 
    //    << " - " << log(delta)
    //    <<" == " << out << endl);

    // Use a scale-free prior: F -= 1/delta
    //   so F' -= -1/delta^2
    
    //    out -= -1/delta/delta;
//    LOG_ERR("Using NO PRIOR on delta\n");

    if (0)
      {

// WORKS! This is a very strong prior that attracts delta towards 5.
//	LOG_ERR("Using a Ga(.0001,50000) prior on delta!\n");
//	const double c = 50000, b = .0001;
//	// DERIVATIVE OF -gammaln(c) + (c-1)*log(delta) - c*log(b) - delta/b;
//	out += (c-1)/delta - 1/b;

	Warning::IssueOnce("Using a Ga(.1,50) prior on delta!");
	const double c = 50, b = .1;
	// DERIVATIVE OF -gammaln(c) + (c-1)*log(delta) - c*log(b) - delta/b;
	out += (c-1)/delta - 1/b;

      }
    else
      {
	Warning::IssueOnce("Not using any prior at all on delta");
      }


    //    LOG_ERR("Using CORRECT scale-free prior on DerivFdDelta\n");

    return out;
}

double SpatialVariationalBayes::OptimizeEvidence(
  // const vector<MVNDist>& fwdPriorVox, // used for parameters other than k
  const vector<MVNDist*>& fwdPosteriorWithoutPrior, // used for parameter k
  //  const vector<SymmetricMatrix>& Si,
  int k, const MVNDist* initialFwdPrior, double guess, bool allowRhoToVary, double* rhoOut) const
{
  Tracer_Plus tr("SpatialVariationalBayes::OptimizeEvidence");

  assert(fwdPosteriorWithoutPrior.at(0) != NULL);
  const int Nparams = fwdPosteriorWithoutPrior[0]->GetSize();
  //const int Nvoxels = fwdPosteriorWithoutPrior.size();
  cout << Nparams << ", " << k << endl;
  assert(Nparams >= 1);
  assert(k <= Nparams);

  DerivEdDelta fcn(covar, fwdPosteriorWithoutPrior, k,  initialFwdPrior, allowRhoToVary);

  LogBisectionGuesstimator guesser;
  //LogRiddlersGuesstimator guesser;

  double hardMin = 0.05; // Much below 0.2 inversion becomes painfully slow
  double hardMax = 1e3;  // Above 1e15, exp(-0.5*1/delta) == 1 (singular)


  double delta = DescendingZeroFinder(fcn)
          .InitialGuess(guess)
    //          .InitialScale(guess/16) // In two guesses, scale = guess (anywhere down to zero).
    //          .ScaleGrowth(4)  // However, increasing from 0.05 to 1000 in one iterataion would take 10 new guesses.
          .InitialScale(guess*0.009) // So hopefully it'll only take two guesses once settled
          .ScaleGrowth(16)  // But it can escape pretty quickly!  7 guesses to get from 0.05 to 1000.
          .SearchMin(hardMin)
          .SearchMax(hardMax)
          .RatioTolX(1.01).MaxEvaluations(2+newDeltaEvaluations)
          .SetGuesstimator(&guesser);

  //Warning::IssueAlways("Increased delta search precision to 1.0001 from 1.01");
  Warning::IssueOnce("Hard limits on delta: [" + stringify(hardMin) + ", " + stringify(hardMax) + "]");

  if (rhoOut != NULL)
    *rhoOut = fcn.OptimizeRho(delta); // 0 if allowRhoToVary == false

  return delta;
}

double SpatialVariationalBayes::OptimizeSmoothingScale( 
    const DiagonalMatrix& covRatio, const ColumnVector& meanDiffRatio, 
    double guess, double* optimizedRho, bool allowRhoToVary,
    bool allowDeltaToVary) const
{
    Tracer_Plus tr("SpatialVariationalBayes::OptimizeSmoothingScale");
    
    DerivFdDelta fcn( covar, covRatio, meanDiffRatio, allowRhoToVary );
    LogBisectionGuesstimator guesser;


    if (bruteForceDeltaSearch) {
      Tracer_Plus tr("SpatialVariationalBayes::OptimizeSmoothingScale - brute force data output");
      
      LOG_ERR("BEGINNING BRUTE-FORCE DELTA SEARCH.\n");
      LOG << "PARAMETERS:\ncovRatio = ["
      	  << covRatio << endl
          << "];\nmeanDiffRatio = [" << meanDiffRatio << "];\n";

      // BELOW: changed cutoff (was 1e16 originally)
      for (double dk = 0.001; dk < 1e4; dk *= (sqrt((2))))
	{
	  cout << "dk = " << dk << endl;
	  LOG  << "BRUTEFORCE=" << dk << "\t"
	       << -0.5*covar.GetC(dk).LogDeterminant().LogValue() << "\t"
	       << -0.5*(covar.GetCinv(dk)*covRatio).Trace() << "\t"
	       << -0.5*(meanDiffRatio.t()*covar.GetCinv(dk)*meanDiffRatio).AsScalar() << endl;
	}
	LOG_ERR("END OF BRUTE-FORCE DELTA SEARCH.\n");

    }
   
    double delta;
    if (allowDeltaToVary)
      {
        delta = DescendingZeroFinder(fcn)
          .InitialGuess(guess)
          .SearchMin(0.2) // Below this, inversion becomes painfully slow
          //.SearchMin(0.01) // Below this, inversion becomes painfully slow
          .SearchMax(1e15) // Above this, exp(-0.5*1/delta) == 1 (singular)
          .RatioTolX(1.01).MaxEvaluations(2+newDeltaEvaluations)
          .SetGuesstimator(&guesser);

//LOG_ERR("HORRIBLE HACK: delta *= 0.95;\n");
//	delta *= 0.95;
      }
    else
      {
        delta = guess;
      }


    if (allowDeltaToVary && optimizedRho != NULL)
      *optimizedRho = fcn.OptimizeRho(delta);

    return delta;
}

const ReturnMatrix CovarianceCache::GetC(double delta) const
{
  Tracer_Plus tr("CovarianceCache::GetC");
  const int Nvoxels = distances.Nrows();

  if (delta == 0)
    return IdentityMatrix(Nvoxels);

  SymmetricMatrix C(Nvoxels);
  for (int a = 1; a <= Nvoxels; a++)
    for (int b = 1; b <= a; b++)
      C(a,b) = exp(-0.5*distances(a,b)/delta);

  // NOTE: when distances = squared distance, prior is equivalent to white
  // noise smoothed with a Gaussian with sigma^2 = 2*delta (haven't actually
  // double-checked this yet).
  // BEWARE: delta is measured in millimeters!! (based on NIFTI file info).

  C.Release(); return C;
}


bool CovarianceCache::GetCachedInRange(double* guess, double lower, double upper, bool allowEndpoints) const
{
  Tracer_Plus tr("CovarianceCache::GetCachedInRange");
  assert(guess != NULL);
  const double initialGuess = *guess;
  if (!(lower < initialGuess && initialGuess < upper))
    {
      cout << "Uh-oh... lower = " << lower
	   << ", initialGuess = " << initialGuess
	   << ", upper = " << upper << endl;

    }
  assert(lower < initialGuess && initialGuess < upper);
  
  Cinv_cache_type::iterator it = Cinv_cache.lower_bound(lower);
  if (it == Cinv_cache.end()) return false;
  if (it->first == lower && !allowEndpoints) it++;
  if (it == Cinv_cache.end()) return false;
  if (it->first > upper) return false;
  if (it->first == upper && !allowEndpoints) return false;

  // Success -- we have at least one fast guess!
  *guess = it->first;

  //  cout << "Found a guess! " << lower << " < " << *guess << " < " << upper << endl;

  // Can we find a better one?
  while (++it != Cinv_cache.end() && it->first <= upper)
    {
      if (it->first == upper && !allowEndpoints) break;

      //      if ( abs(it->first - initialGuess) < abs(*guess - initialGuess) )
      if (it->first < initialGuess
	  || it->first - initialGuess < initialGuess - *guess)
	*guess = it->first;

      //      cout << "Improved guess! " << lower << " < " << *guess << " < " << upper << endl;
    }

  assert(lower < *guess && *guess < upper);

  return true;
}

const SymmetricMatrix& CovarianceCache::GetCinv(double delta) const
{
  Tracer_Plus tr("CovarianceCache::GetCinv");
  if (Cinv_cache[delta].Nrows() == 0)
    {

#ifdef NOCACHE
	Warning::IssueOnce("Cache is disabled to avoid memory problems!");
	Cinv_cache.clear();
#endif


      //      cout << "[" << flush;
      //      cout << "GetCinv cache miss... " << flush;
      //      const int Nvoxels = distances.Nrows();
      //      SymmetricMatrix C(Nvoxels);
      //      for (int a = 1; a <= Nvoxels; a++)
      //	for (int b = 1; b <= a; b++)
      //	  C(a,b) = exp(-0.5*distances(a,b)/delta);
      //      Cinv_cache[delta] = C.i();
      Cinv_cache[delta] = GetC(delta).i();
      //      cout << "done." << endl;
      //      cout << "]" << flush;
    }
  else
    {
      //      cout << "GetCinv cache hit!\n";
    }

  return Cinv_cache[delta];
}

const SymmetricMatrix& CovarianceCache::GetCiCodistCi(double delta, 
		double* CiCodistTrace) const
{
  if (CiCodistCi_cache[delta].first.Nrows() == 0)
    {
#ifdef NOCACHE
	CiCodistCi_cache.clear();
#endif
      //      cout << "{" << flush;
      GetCinv(delta); // for sensible messages, make sure cache hits
      //cout << "GetCiCodistCi cache miss... " << flush;
      Matrix CiCodist = GetCinv(delta) * SP(GetC(delta), distances);
      CiCodistCi_cache[delta].second = CiCodist.Trace();
      Matrix CiCodistCi_tmp = CiCodist*GetCinv(delta);
      CiCodistCi_cache[delta].first << CiCodistCi_tmp; // Force symmetric
    
      { // check something
	double maxAbsErr = 
	  (CiCodistCi_cache[delta].first - CiCodistCi_tmp).MaximumAbsoluteValue();
	if (maxAbsErr > CiCodistCi_tmp.MaximumAbsoluteValue() * 1e-5 )
	  // If that test fails, you're probably in trouble.
	  // Reducing it to e.g. 1e-5 (to make dist2 work) 
	  //   => non-finite alpha in iteration 2
	  // Reduced it to 1e-5 to get mdist to work....
	  //   => same result.  (oops, was mdist2)
	  // Reduced it to 1e-5 to get mdist to work, again...
	  //   => 
	  {
	    LOG_ERR("In GetCiCodistCi -- matrix not symmetric!\nError = "
		    << maxAbsErr << ", maxabsvalue = " 
		    << CiCodistCi_tmp.MaximumAbsoluteValue() << endl);
	    assert(false);
	  }
      }      
      //      cout << "}" << flush;
    }

  if (CiCodistTrace != NULL) 
    (*CiCodistTrace) = CiCodistCi_cache[delta].second;
  return CiCodistCi_cache[delta].first;
}

