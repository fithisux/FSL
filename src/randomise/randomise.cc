/*  randomise.cc
    Tim Behrens & Steve Smith & Matthew Webster (FMRIB) & Tom Nichols (UMich)
    Copyright (C) 2004-2010 University of Oxford  */
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

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1
#define CLUST_CON 26

// 26 for FSL 18 for SPM
#include "miscmaths/f2z.h"
#include "newimage/newimageall.h"
#include "libprob.h"
#include "ranopts.h"

#include <algorithm>

using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;
using namespace RANDOMISE;

class VoxelwiseDesign
{
public:
  bool isSet;
  vector<Matrix> EV;
  vector<int> location;
  vector<Matrix> designAtVoxel;
  vector< vector<int> > dofAtVoxel;
  vector<Matrix> contrastAtVoxel;
  void setup(const vector<int>& voxelwise_ev_numbers,const vector<string>& voxelwise_ev_filenames,const volume<float>& mask, const int maximumLocation, const bool isVerbose);
  VoxelwiseDesign() : isSet(false), storingByVoxel(false) {}
  void storeByVoxel(size_t nVoxels) { storingByVoxel=true; EV.clear(); location.clear(); designAtVoxel.resize(nVoxels+1); dofAtVoxel.resize(nVoxels+1); contrastAtVoxel.resize(nVoxels+1);}
  Matrix adjustDesign(const Matrix& originalDesign,const int voxelNo);
private:
  bool storingByVoxel;
};

Matrix VoxelwiseDesign::adjustDesign(const Matrix& originalDesign, const int voxelNo)
{
  if ( storingByVoxel ) 
    return designAtVoxel[voxelNo];
  Matrix newDesign(originalDesign);
  for (unsigned int currentEV=0;currentEV<EV.size();currentEV++) {
    if ( location[currentEV] < 0 )
      newDesign.Columns( abs(location[currentEV]),abs(location[currentEV])+originalDesign.Nrows()-1 )=EV[currentEV].Column(voxelNo).AsDiagonal();
    else
       newDesign.Column(location[currentEV])=EV[currentEV].Column(voxelNo);
  }
  return newDesign;
}

void VoxelwiseDesign::setup(const vector<int>& EVnumbers,const vector<string>& EVfilenames,const volume<float>& mask,const int maximumLocation,const bool isVerbose)
{
  isSet=false;
  if(EVnumbers.size() != EVfilenames.size())
    throw Exception("Number of input voxelwise_ev_filenames must match number of voxelwise_ev_numbers");
  location=EVnumbers;
  EV.resize(EVfilenames.size());
  volume4D<float> input;
  for(unsigned int i=0; i<EV.size(); i++)      
  {
    if(location[i]>maximumLocation)
      throw Exception("voxelwise_ev_numbers option specifies a number greater than number of design EVs)");
    if (isVerbose) cout << "Loading voxelwise ev: " << EVfilenames.at(i) << " for EV " << location.at(i) << endl;
    read_volume4D(input,EVfilenames.at(i));
    EV.at(i)=input.matrix(mask);
  }
  isSet=true;
}

VoxelwiseDesign voxelwiseInput;

class Permuter
{ 
public:
  bool isFlipping;
  bool isRandom;
  bool isPermutingBlocks;
  int nBlocks;
  int nSubjects;
  double finalPermutation;
  vector<double>       uniquePermutations; //0 is unique for whole design, 1..nBlocks is unique per block
  vector<ColumnVector> permutedLabels;
  vector<ColumnVector> originalLabels;
  vector<ColumnVector> originalLocation;
  vector<ColumnVector> previousPermutations;
  ColumnVector truePermutation;
  ColumnVector unpermutedVector;
  Permuter();
  ~Permuter();
  void writePermutationHistory(const string& filename);
  void createPermutationScheme(const Matrix& design, ColumnVector groups, const bool oneNonZeroContrast, const long requiredPermutations, const bool detectingNullElements, const bool outputDebug, const bool permuteBLocks=false,const bool forceFlipping=false);
  void initialisePermutationBlocks(const ColumnVector& labels, const long requiredPermutations);
  ColumnVector createDesignLabels(const Matrix& design);
  void createTruePermutation(const ColumnVector& labels, ColumnVector copyOldlabels, ColumnVector& permvec);
  ColumnVector nextPermutation(const long perm);
  ColumnVector nextPermutation(const long permutationNumber, const bool printStatus, const bool isStoring);
  bool isPreviousPermutation(const ColumnVector& newPermutation);
  ColumnVector permutationVector();
  double reportRequiredPermutations(const bool printToScreen);
  ColumnVector returnPreviousTruePermutation(const long permutationNumber);
  ColumnVector returnPreviousTruePermutation(const long permutationNumber, ColumnVector& previousState);

private:
  Permuter *blockPermuter;
  double computeUniquePermutations(const ColumnVector& labels, const bool calculateFlips);
  void nextShuffle(ColumnVector& perm);
  void nextFlip(ColumnVector& mult);
};

class ParametricStatistic
{
public:
  Matrix originalStatistic,uncorrectedStatistic,maximumDistribution,sumStatMat,sumSampMat;
  bool isAveraging,storingUncorrected;
  string outputName;
  void   store(const volume<int>& clusterLabels, const ColumnVector& clusterSizes, const volume<float>& mask, const int contrastNo, const unsigned long permNo, const bool outputRaw);
  void   store(const Matrix& parametricMatrix, const unsigned long permNo,const volume<float> *mask, const bool outputRaw);
  void   setup(const int nContrasts,const unsigned long nPerms, const int nVoxels, const bool wantAverage, const bool wantUncorrected, const string outputFileName);
  void   average(const string filename, const float percentileThreshold,const volume<float>& mask);
  ParametricStatistic() { isAveraging=false; }
};

void ParametricStatistic::setup(const int nContrasts,const unsigned long nPerms, const int nVoxels, const bool wantAverage,const bool wantUncorrected=false,const string outputFileName="")
{
  outputName=outputFileName;
  isAveraging=wantAverage;
  storingUncorrected=wantUncorrected;
  maximumDistribution.ReSize(nContrasts,nPerms);
  maximumDistribution=0;
  if ( storingUncorrected ) {
    uncorrectedStatistic.ReSize(1,nVoxels);
    uncorrectedStatistic=0;
  }
  originalStatistic.ReSize(nContrasts,nVoxels);
  originalStatistic=0;
  if ( isAveraging ) {
    sumStatMat=originalStatistic;
    sumSampMat=originalStatistic;
  }
}

void ParametricStatistic::store(const volume<int>& clusterLabels, const ColumnVector& clusterSizes , const volume<float>& mask, const int contrastNo, const unsigned long permNo, const bool outputtingRaw=false)
{
  if ( clusterSizes.Nrows() > 0 ) 
    maximumDistribution(contrastNo,permNo)=clusterSizes.MaximumAbsoluteValue(); 
  if ( permNo==1 || isAveraging || outputtingRaw ) { 
    volume4D<float> parametricImage(mask.xsize(),mask.ysize(),mask.zsize(),1);
    parametricImage=0;
    for(int z=0; z<mask.zsize(); z++)
      for(int y=0; y<mask.ysize(); y++)
	for(int x=0; x<mask.xsize(); x++)
	  if( clusterLabels(x,y,z) ) 
	    parametricImage(x,y,z,0)=clusterSizes(clusterLabels(x,y,z));	  
    if (permNo==1) 
      originalStatistic.Row(contrastNo)=parametricImage.matrix(mask);
    if (isAveraging) {
      sumStatMat.Row(contrastNo)+=parametricImage.matrix(mask);
      sumSampMat.Row(contrastNo)+=SD(parametricImage.matrix(mask),parametricImage.matrix(mask));
    }
    if ( outputtingRaw )
      save_volume4D( parametricImage, outputName+"_perm"+(num2str(permNo).insert(0,"00000")).erase(0,num2str(permNo).length()) );
  }
}

void ParametricStatistic::store(const Matrix& parametricMatrix, const unsigned long permNo, const volume<float> *mask=NULL, const bool outputtingRaw=false )
{
  maximumDistribution.Column(permNo)=max(parametricMatrix.t()).t();
  if (permNo==1) 
    originalStatistic=parametricMatrix;
  if (storingUncorrected) 
    uncorrectedStatistic += gt(originalStatistic,parametricMatrix);
  if (isAveraging) {
    sumStatMat+=parametricMatrix;
    sumSampMat+=SD(parametricMatrix,parametricMatrix);
  }
  if ( outputtingRaw ) {
     volume4D<float> rawImage;
     rawImage.setmatrix( parametricMatrix, *mask );
     save_volume4D( rawImage, outputName+"_perm" + (num2str(permNo).insert(0,"00000")).erase(0,num2str(permNo).length()) );
  }
  
}

void ParametricStatistic::average(const string filename, const float percentileThreshold,const volume<float>& mask)
{
  if (isAveraging) {
    volume4D<float> temp;
    temp.setmatrix(sumStatMat,mask);
    save_volume4D(temp,filename+"sum");
    temp.setmatrix(sumSampMat,mask);
    save_volume4D(temp,filename+"samp");
    sumStatMat=SD(sumStatMat,sumSampMat);

    if (percentileThreshold>0) {
      temp.setmatrix(sumStatMat,mask);
      float min=temp.percentile(percentileThreshold,mask);
      //cerr << min << " " << percentile((Matrix)tstat_ceav.t(),percentileThreshold*100) << endl;     
      for(int i=1;i<=sumStatMat.Ncols();i++)
	if(sumStatMat(1,i)<min) sumStatMat(1,i)=min;
    }

    temp.setmatrix(sumStatMat,mask);
    save_volume4D(temp,filename+"post");
  }
}

Matrix tfce(const Matrix& tstat, const volume<float>& mask, const float delta, float height_power, float size_power, int connectivity){
  volume4D<float> spatialStatistic;
  spatialStatistic.setmatrix(tstat,mask);
  tfce(spatialStatistic[0],height_power,size_power,connectivity,0,delta);
  return(spatialStatistic.matrix(mask));
}

void clusterStatistic(ParametricStatistic& output, const Matrix& inputStatistic, const volume<float>& mask, const float threshold, const int permutationNo, const bool outputPerms)
{
ColumnVector clusterSizes;
volume4D<float> spatialStatistic;  
   spatialStatistic.setmatrix(inputStatistic,mask);
   spatialStatistic.binarise(threshold);
   volume<int> clusterLabels=connected_components(spatialStatistic[0],clusterSizes,CLUST_CON);
   output.store(clusterLabels,clusterSizes,mask,1,permutationNo,outputPerms);
}

void clusterMassStatistic(ParametricStatistic& output, const Matrix& inputStatistic, const volume<float>& mask, const float threshold, const int permutationNo, const bool outputPerms)
{
ColumnVector clusterSizes;
volume4D<float> spatialStatistic, originalSpatialStatistic;  
   spatialStatistic.setmatrix(inputStatistic,mask);
   originalSpatialStatistic=spatialStatistic;
   spatialStatistic.binarise(threshold);
   volume<int> clusterLabels=connected_components(spatialStatistic[0],clusterSizes,CLUST_CON);
   clusterSizes=0;	
   for(int z=0; z<mask.zsize(); z++)
     for(int y=0; y<mask.ysize(); y++)
       for(int x=0; x<mask.xsize(); x++)
	 if(clusterLabels(x,y,z)>0)
	   clusterSizes(clusterLabels(x,y,z))=clusterSizes(clusterLabels(x,y,z))+originalSpatialStatistic[0](x,y,z);
   output.store(clusterLabels,clusterSizes,mask,1,permutationNo,outputPerms);
}

Matrix tfceStatistic(ParametricStatistic& output, const Matrix& inputStatistic, const volume<float>& mask, float& tfceDelta, const float tfceHeight, const float tfceSize, const int tfceConnectivity, const int permutationNo, const bool isF, const int numContrasts, const vector<int>& dof, const bool outputPerms, const bool overrideDelta)
{
  Matrix tstat_ce(inputStatistic);
  if ( isF ) { 
    ColumnVector zstat, dofVector(inputStatistic.AsColumn());
    dofVector=dof[0];
    F2z::ComputeFStats( tstat_ce.AsColumn(), numContrasts, dofVector, zstat);
    tstat_ce=zstat.AsRow();
  }

  if (permutationNo==1 && !overrideDelta) {
     tfceDelta=tstat_ce.Maximum()/100.0;  // i.e. 100 subdivisions of the max input stat height
     if ( tfceDelta <= 0 )
       cout << "Warning: The unpermuted statistic image for the current image contains no positive values, and cannot be processed with TFCE. A blank output image will be created." << endl;
  }
 
  if ( tfceDelta > 0 ) 
    tstat_ce=tfce(tstat_ce,mask,tfceDelta,tfceHeight,tfceSize,tfceConnectivity);
  else 
    tstat_ce=0;

  output.store(tstat_ce, permutationNo,&mask,outputPerms);
  return (tstat_ce.Row(1));
}

void checkInput(const short st,const  Matrix& dm,const  Matrix& tc,const  Matrix& fc) {
  if (dm.Nrows()!=st) throw Exception("number of rows in design matrix doesn't match number of \"time points\" in input data!"); 
  if (tc.Ncols()!=dm.Ncols()) throw Exception("number of columns in t-contrast matrix doesn't match number of columns in design matrix!");
  if (fc.Ncols() !=0 && fc.Ncols()!=tc.Nrows()) throw Exception("number of columns in f-contrast matrix doesn't match number of rows in t-contrast matrix!");
}

volume<float> nonConstantMask(volume4D<float>& data, const bool allOnes)
{
  volume<float> nonConstantMask(data.xsize(),data.ysize(),data.zsize());
  nonConstantMask.copyproperties(data[0]);
  if ( allOnes ) {
    nonConstantMask=1;
    return nonConstantMask;
  }
  nonConstantMask=0;
  for(int z=0; z<data.zsize(); z++)
    for(int y=0; y<data.ysize(); y++)
      for(int x=0; x<data.xsize(); x++)
	for(int t=1; t<data.tsize(); t++)
	{
	  if ( data(x,y,z,t)!=data(x,y,z,0) ) {
	    nonConstantMask(x,y,z)=1;
	    break;
	  }
	}
  return nonConstantMask;
}

void Initialise(ranopts& opts, volume<float>& mask, Matrix& datam, Matrix& tc, Matrix& dm, Matrix& fc, Matrix& gp, Matrix& effectiveDesign)
{
  if (opts.tfce2D.value()) {
    opts.tfce.set_value("true");
    opts.tfce_height.set_value("2");     
    opts.tfce_size.set_value("1");     
    opts.tfce_connectivity.set_value("26");  
  }

  if ( opts.n_perm.value() < 0 ) {
    throw Exception(("Randomise requires a postive number of permutations, did you mean to type -n "+num2str(opts.n_perm.value()).erase(0,1)+"?").c_str());
  }

  if ( opts.randomSeed.set() ) srand(opts.randomSeed.value() );
  if ( opts.randomSeed.set() && opts.verbose.value() ) cout << "Seeding with " << opts.randomSeed.value() << endl;
  if (opts.verbose.value()) cout << "Loading Data: "; 
  volume4D<float> data;
  read_volume4D(data,opts.in_fileroot.value());
  
  if(opts.one_samp.value())
  {
    dm.ReSize(data.tsize(),1);
    dm=1;
    tc.ReSize(1,1);
    tc=1;
  }
  else if ( opts.dm_file.value()=="" || opts.tc_file.value()=="" ) throw Exception("Randomise requires a design matrix and contrast as input");
  if (opts.dm_file.value()!="") dm=read_vest(opts.dm_file.value());
  if (opts.tc_file.value()!="") tc=read_vest(opts.tc_file.value());
  if (opts.fc_file.value()!="") fc=read_vest(opts.fc_file.value());
  if (opts.gp_file.value()!="") gp=read_vest(opts.gp_file.value());
  else {
    gp.ReSize(dm.Nrows(),1);
    gp=1;
  }
  vector<bool> groupListed((int)gp.Maximum()+1,false);
  for ( int i=1; i<=gp.Nrows() ; i++ )
    groupListed[(int)gp(i,1)]=true;
  for ( int i=1; i<=gp.Maximum(); i++ )
    if ( !groupListed[i] )
      throw Exception(("Error: block "+num2str(i)+" must be assigned to at least one design row in the blocks file.").c_str());

  if (opts.effectiveDesignFile.value()!="") effectiveDesign=read_vest(opts.effectiveDesignFile.value());
  if ( opts.nMultiVariate.value() == 1 ) checkInput(data.tsize(),dm,tc,fc);  // should do a different check in the Multivariate case!

  if (opts.parallelData.value()) {
    int permutationsPerContrast(300); 
    if (opts.tfce.value()) permutationsPerContrast=100;
    if (opts.voxelwise_ev_numbers.set() && opts.voxelwise_ev_filenames.set()) permutationsPerContrast=100;

    cout << opts.n_perm.value() << " " << int(tc.Nrows()+fc.Nrows()) << " " << opts.out_fileroot.value() << " " << permutationsPerContrast << endl;
    exit(0);
  }

  if (opts.nMultiVariate.set()) {
    // Read in 4D data - see error message below for format details
    if ((data.ysize()!=opts.nMultiVariate.value()) || (data.zsize()!=1) || (data.tsize()!=dm.Nrows())) { 
      throw Exception("Multi-Variate input data of wrong size!\nSize must be: N x k x 1 x M\n   where N=#vertices, k=#multi-variate dims, M=#subjects");
    }
    // make data matrix of concatenated components, with all of component 1, then all of component 2, etc.
    // this way the first sx values represent a whole mesh/volume of values for a component
    // and the output stats can just be taken as the first set of values (with corresponding row indices)
    datam.ReSize(data.tsize(),data.xsize()*opts.nMultiVariate.value());
    for (int t=1; t<=data.tsize(); t++) {
      for (int n=1; n<=opts.nMultiVariate.value(); n++) {
	for (int x=1; x<=data.xsize(); x++) {
	  datam(t,x+(n-1)*data.xsize())=data(x-1,n-1,0,t-1);
	}
      }
    }
    // dummy mask (if needed) of size of output
    mask.reinitialize(data.xsize(),1,1);
    mask=1.0f;
  } else {
    if (opts.maskname.value()!="") {
      read_volume(mask,opts.maskname.value());
      if (!samesize(data[0],mask)) throw Exception("Mask dimensions do not match input data dimensions!");
    }
    else mask=nonConstantMask(data, opts.disableNonConstantMask.value() );
    if ( mask.sum() < 1 ) throw Exception("Data mask is blank.");
    datam=data.matrix(mask);
    if (opts.demean_data.value()) datam=remmean(datam);
  }

  if ( datam.Nrows() == 0 || datam.Ncols() == 0 )
    throw Exception("No data voxels present.");

  if (opts.verbose.value()) cout << endl;

  if (opts.voxelwise_ev_numbers.set() && opts.voxelwise_ev_filenames.set())
    voxelwiseInput.setup(opts.voxelwise_ev_numbers.value(),opts.voxelwise_ev_filenames.value(),mask,dm.Ncols(),opts.verbose.value());
  
  if (opts.verbose.value()) cout << "Data loaded" << endl;
}

Matrix PermutedDesign(const Matrix& originalDesign,const ColumnVector& permutation,const bool multiply){
  Matrix output=originalDesign;
  for(int row=1;row<=originalDesign.Nrows();row++)
  {
    if (multiply) output.Row(row)=originalDesign.Row(row)*permutation(row);
    else output.Row(row) << originalDesign.Row(int(permutation(row)));
  }
  return output;
}


ReturnMatrix sumSqr(const Matrix& mat) {
  Matrix res(1,mat.Ncols());
  for (int mc=1; mc<=mat.Ncols(); mc++) {
    double sumsqr(0);
    for (int mr=1; mr<=mat.Nrows(); mr++) {
      double foo(mat(mr,mc));
      sumsqr += foo*foo;
    }
    res(1,mc)=sumsqr;
  }
res.Release();
return res;
  }

ReturnMatrix sumSqr2(const Matrix& mat) {
  Matrix res(1,mat.Nrows());
  res=1.0;
  res*=SP(mat,mat);
  res.Release();
  return res;
  }


Matrix calculateTstat(const Matrix& data, const Matrix& model, const Matrix& tc, Matrix& estimate, Matrix& residuals, Matrix& sigmaSquared, const int dof)
{
  Matrix pinvModel(pinv(model)); // inverted model used several times
  estimate=pinvModel*data;
  residuals=data-model*estimate;
  estimate=tc*estimate; //estimate now is cope
  sigmaSquared=sumSqr2(residuals)/(float)dof;
  residuals=diag(tc*pinvModel*pinvModel.t()*tc.t())*sigmaSquared; //residuals now is varcope
  return(SD(estimate,sqrt(residuals)));
}

Matrix calculateFStat(const Matrix& data, const Matrix& model, const Matrix& contrast, const int dof,const int rank)
{ 
  // model is N_subject by N_ev
  // data is N_subject by N_voxels
  Matrix pinvModel(pinv(model)); // inverted model used several times
  Matrix estimate = pinvModel*data;
  Matrix residuals= data-model*estimate;
  residuals = sum(SP(residuals,residuals))/(float)dof; //residuals now hold sigmasquared
  estimate = pinv((contrast*pinvModel).t()).t()*contrast*estimate;
  estimate = sum(SP(estimate,estimate))/(float)rank;   
  return(SD(estimate,residuals));
}        

Matrix smoothTstat(const Matrix inputSigmaSquared,const volume<float>& mask,const volume<float>& smoothedMask, const float sigma_mm)
{
  volume4D<float> sigsqvol;
  sigsqvol.setmatrix(inputSigmaSquared,mask);
  sigsqvol[0]=smooth(sigsqvol[0],sigma_mm);
  sigsqvol[0]/=smoothedMask;
  Matrix newSigmaSquared=sigsqvol.matrix(mask);
  return(SD(newSigmaSquared,inputSigmaSquared));
}

void OutputStat(const ParametricStatistic input,const volume<float>& mask, const int nPerms,string statLabel,const string fileRoot,const bool outputText, const bool outputRaw=true, const bool writeCritical=true)
{ 
 volume4D<float> output(mask.xsize(),mask.ysize(),mask.zsize(),1);
 long nVoxels(input.originalStatistic.Ncols());
 Matrix currentStat(1,nVoxels);
 output.setmatrix(input.originalStatistic.Row(1),mask);
 if (outputRaw) save_volume4D(output,fileRoot+statLabel);
 RowVector distribution = input.maximumDistribution.Row(1);    
 if (outputText)
 {
   ofstream output_file((fileRoot+"_corrp"+statLabel+".txt").c_str());
   output_file << distribution.t();
   output_file.close();
 }
 SortAscending(distribution);

 int criticalLocation=(int)ceil(0.95*nPerms);
 double criticalValue=distribution(criticalLocation);
 if ( writeCritical )
   cout << "Critical Value for: " << fileRoot+"_corrp"+statLabel << " is: " << criticalValue << endl;

 currentStat=0;
 for(int i=1; i<=nVoxels; i++)
   for(int j=nPerms; j>=1; j--) //N.B. it's probably safe to start at nPerms-1, which would also help the '1's bug
     if ((float)input.originalStatistic(1,i)>(float)distribution(j))
     {
       currentStat(1,i) = float(j)/nPerms;
       j=0;
     }
 output.setmatrix(currentStat,mask);
 save_volume4D(output,fileRoot+"_corrp"+statLabel);
 if (input.storingUncorrected) {
   output.setmatrix(input.uncorrectedStatistic.Row(1)/float(nPerms),mask);
   save_volume4D(output,fileRoot+"_p"+statLabel);
 }
}           

float MVGLM_fit(const Matrix& X, const Matrix& Y, const Matrix& contrast, vector<int>& dof)
{
  // adapted by Mark Jenkinson from code in first_utils (by Brian Patenaude)
  // Y is data : N_subject x 3
  // X is design : N_subject x N_ev
  // contrast is : N_con x N_ev
  // g = Y.Ncols = 3
  // p = X.Ncols = N_ev
  // N = Y.Nrows = N_subjects

  // Calculate estimated values
  Matrix Yhat=X*(X.t()*X).i()*X.t()*Y;
  // Calculate R0 (residual) covariance matrix
  Matrix R0=Y-Yhat;
  R0=R0.t()*R0;
  // Calculate R1, the sum-square /cross square product for hypothesis test
  Matrix Yhat1= X*contrast.t()*(contrast*X.t()*X*contrast.t()).i()*contrast*X.t()*Y;
  Matrix R1=Y-Yhat1;
  // Not efficient but easy to convert to other statistics
  R1=R1.t()*R1-R0;
	
  // Calculate Pillai F
  int g=Y.Ncols();
  int p=X.Ncols();//number of dependant
  int N=Y.Nrows();//total sample size
  float F=0, df2=0,df1=0;
				
  float pillai=(R1*(R1+R0).i()).Trace();
				
  int s=1;
  if (p<(g-1)) {s=p;}
  else {s=g-1;}
  float t=(abs(p-g-1)-1)/2.0;
  float u=(N-g-p-1)/2.0;				
  F=((2*u+s+1)/(2*t+s+1))*(pillai/(s-pillai));
  df1=s*(2*t+s+1);
  df2=s*(2*u+s+1);
  if (dof.size()!=2) dof.resize(2);
  dof[0]=df1;  dof[1]=df2;
  //    cout<<"Pillai F "<<pillai<<" "<<F<<" "<<df1<<" "<<df2<<endl;
  return F;
}


Matrix calculateMultiVariateFStat(const Matrix& model, const Matrix& data, vector<int>& dof, int nMultiVariate)
{ 
  // model is N_subject by N_ev
  // data is N_subject by (N_vertex * 3)
  // dof[0] for numerator (F) and dof[1] for denominator  - but are these ever needed?
  int nvert=data.Ncols()/nMultiVariate;
  int nsubj=data.Nrows();
  int nev=model.Ncols();
  Matrix Fstat(1,nvert), datav(nsubj,3), contrast(nev,nev);
  contrast=IdentityMatrix(nev);  // is this what is needed after initial mangled of design?!?
  for (int n=1; n<=nvert; n++) {
    for (int r=1; r<=nsubj; r++) {
      datav(r,1)=data(r,n);
      datav(r,2)=data(r,n+nvert);
      datav(r,3)=data(r,n+2*nvert);
    }
    Fstat(1,n) = MVGLM_fit(model,datav,contrast,dof);
  }
  return Fstat;
}

Matrix evaluateStatistics(const Matrix& data,const Matrix& model,const Matrix& contrast, Matrix& cope, Matrix& varcope, Matrix& sigmaSquared,vector<int>& dof,const int rank,const int multiVariate, const bool doingF)
{
  if ( doingF ) {
    if ( multiVariate > 1 )
      return calculateMultiVariateFStat(model, data, dof, multiVariate);
    else
      return calculateFStat(data, model, contrast, dof[0], rank);
  } else 
    return calculateTstat(data, model, contrast, cope, varcope, sigmaSquared, dof[0]);
}

void calculatePermutationStatistics(ranopts& opts, const volume<float>& mask, Matrix& datam, const Matrix& tc, Matrix& dm,int tstatnum, vector<int>& dof, Permuter& permuter, VoxelwiseDesign& voxelwiseDesign)
{
  int nVoxels=(int)no_mask_voxels(mask);
  int rankF=MISCMATHS::rank(tc.t());
  if ( opts.isDebugging.value() ) {
    cerr << "Input Design: " << endl << dm << endl;
    cerr << "Input Contrast: " << endl << tc << endl;
    cerr << "Contrast rank: " << rankF << endl;
    cerr << "Dof: " << dof[0] << " original dof: " << ols_dof(dm) << endl;
  }
  volume4D<float> tstat4D(mask.xsize(),mask.ysize(),mask.zsize(),1);
  float tfce_delta(0), clusterThreshold(0), massThreshold(0);
  if ( opts.tfce_delta.set() )
    tfce_delta=opts.tfce_delta.value();
  if (tstatnum>=0) clusterThreshold=opts.cluster_thresh.value();
  else clusterThreshold=opts.f_thresh.value();
  if (tstatnum>=0) massThreshold=opts.clustermass_thresh.value();
  else massThreshold=opts.fmass_thresh.value();
  bool isNormalising( opts.cluster_norm.value() && tstatnum >=0 ), lowram(false);

  string statLabel;
  if (tstatnum<0) 
    statLabel="_fstat"+num2str(-tstatnum);
  else 
    statLabel="_tstat"+num2str(tstatnum);
  // prepare smoothed mask for use (as a convolution renormaliser) in variance smoothing if required
  volume<float> smoothedMask;
  if(opts.var_sm_sig.value()>0) 
    smoothedMask=smooth(binarise(mask,(float)0.5),opts.var_sm_sig.value());
  // containers for different inference distribution
  ParametricStatistic clusters, clusterMasses, clusterNormals, clusterEnhanced, clusterEnhancedNormals, voxels;
  Matrix dmperm, tstat(1,nVoxels), pe(dm.Ncols(),nVoxels), cope(1,nVoxels), varcope(1,nVoxels), sigmaSquared(1,nVoxels), previousTFCEStat;
  unsigned long nPerms=(unsigned long)permuter.reportRequiredPermutations(opts.verbose.value());

  if ( !((clusterThreshold>0) || (massThreshold>0) || opts.tfce.value() || opts.voxelwiseOutput.value()) )
  {
    cout << "Warning! No output options selected. Outputing raw tstat only" << endl;
    opts.outputRaw.set_value("true");
    nPerms=1;
  }
  // resize the containers for the relevant inference distributions
  voxels.setup(1,nPerms,nVoxels,false,opts.outputUncorr.value(),opts.out_fileroot.value()+"_vox"+statLabel);
  if ( clusterThreshold >0 )  
    clusters.setup(1,nPerms,nVoxels,isNormalising,opts.outputUncorr.value(),opts.out_fileroot.value()+"_clustere"+statLabel);
  if ( massThreshold>0 ) 
    clusterMasses.setup(1,nPerms,nVoxels,false,opts.outputUncorr.value(),opts.out_fileroot.value()+"_clusterm"+statLabel);
  if ( clusters.isAveraging ) 
    clusterNormals.setup(1,nPerms,nVoxels,false,opts.outputUncorr.value(),opts.out_fileroot.value()+"_clustern"+statLabel);
  if ( opts.tfce.value() )    
    clusterEnhanced.setup(1,nPerms,nVoxels,isNormalising,opts.outputUncorr.value(),opts.out_fileroot.value()+"_tfce"+statLabel);
  if ( clusterEnhanced.isAveraging ) { 
    clusterEnhancedNormals.setup(1,nPerms,nVoxels,false,opts.outputUncorr.value(),opts.out_fileroot.value()+"_tfcen"+statLabel);
    try { previousTFCEStat.ReSize(nPerms,clusterEnhanced.sumStatMat.Ncols());} //between 5e17 - 5e18 values for a 2gb machine
    catch (...) {cerr << "using lowram" << endl; lowram=true;}         
  }
  
  for(unsigned long perm=1; perm<=nPerms; perm++) {

    ColumnVector permvec = permuter.nextPermutation(perm,opts.verbose.value(), isNormalising || opts.outputTextPerm.value());   
    dmperm=PermutedDesign(dm,permvec,permuter.isFlipping);   

    if (voxelwiseDesign.isSet)
      for(int voxel=1;voxel<=datam.Ncols();voxel++) {
	Matrix dmtemp(voxelwiseDesign.adjustDesign(dm,voxel)), sigmaTemp(1,1), copeTemp(1,1), varcopeTemp(1,1);
	dmperm=PermutedDesign(dmtemp,permvec,permuter.isFlipping);
	tstat.Column(voxel)=evaluateStatistics(datam.Column(voxel), dmperm, voxelwiseDesign.contrastAtVoxel[voxel], copeTemp, varcopeTemp, sigmaTemp, voxelwiseDesign.dofAtVoxel[voxel], rankF, opts.nMultiVariate.value(), (tstatnum < 0) );
	sigmaSquared.Column(voxel)=sigmaTemp;
	cope.Column(voxel)=copeTemp;
	varcope.Column(voxel)=varcopeTemp;
      }
    else 
      tstat=evaluateStatistics(datam, dmperm, tc, cope, varcope, sigmaSquared, dof, rankF, opts.nMultiVariate.value(), (tstatnum < 0) );

    if ( opts.outputGlm.value() && perm == 1 && tstatnum>0 ) {
       volume4D<float> temp4D;
       Matrix estimate( pinv(dm)*datam ); //regenerate paramter estimates - make efficient later
       temp4D.setmatrix(estimate,mask);
       save_volume4D(temp4D,opts.out_fileroot.value()+"_glm_pe"+statLabel);
       temp4D.setmatrix(cope,mask);
       save_volume4D(temp4D,opts.out_fileroot.value()+"_glm_cope"+statLabel);
       temp4D.setmatrix(varcope,mask);
       save_volume4D(temp4D,opts.out_fileroot.value()+"_glm_varcope"+statLabel);
       temp4D.setmatrix(sigmaSquared,mask);
       save_volume4D(temp4D,opts.out_fileroot.value()+"_glm_sigmasqr"+statLabel);    
    }
 
    if( opts.var_sm_sig.value()>0 && tstatnum > 0 ) 
      tstat=SD(tstat,sqrt(smoothTstat(sigmaSquared,mask,smoothedMask,opts.var_sm_sig.value())));
    if ( opts.isDebugging.value() ) 
      cerr << "statistic Maximum: " << tstat.Maximum() << endl;
    voxels.store(tstat,perm,&mask,opts.output_permstat.value());
    if (opts.tfce.value())
    {
      Matrix tfceOutput=tfceStatistic(clusterEnhanced,tstat,mask,tfce_delta,opts.tfce_height.value(),opts.tfce_size.value(),opts.tfce_connectivity.value(),perm,(tstatnum<0),tc.Nrows(),dof,opts.output_permstat.value(), opts.tfce_delta.set());
      if(!lowram && clusterEnhanced.isAveraging ) previousTFCEStat.Row(perm)=tfceOutput.Row(1);
    }
    if ( clusterThreshold > 0 ) 
      clusterStatistic(clusters,tstat,mask,clusterThreshold,perm,opts.output_permstat.value());
    
    if ( massThreshold > 0 ) 
      clusterMassStatistic(clusterMasses,tstat,mask,massThreshold,perm,opts.output_permstat.value());
  }
  //End of Permutations
    
  //Rerun perms for clusternorm
  if ( isNormalising )
  { 
    volume4D<float> temp4D;
    if ( clusters.isAveraging ) {
      clusters.average(opts.out_fileroot.value()+statLabel+"_clusternorm",0,mask);
      temp4D.setmatrix(clusters.sumStatMat,mask);
    }
    if (clusterEnhanced.isAveraging)
      clusterEnhanced.average(opts.out_fileroot.value()+statLabel+"_tfcenorm",0.02,mask);

    for(unsigned long perm=1; perm<=nPerms; perm++)
    {
      if (opts.verbose.value()) cout << "Starting second-pass " << perm << endl;
      if ( clusters.isAveraging || ( clusterEnhanced.isAveraging && lowram ) ) //Regenerate stats
      { 
	ColumnVector permvec=permuter.returnPreviousTruePermutation(perm);
	dmperm=PermutedDesign(dm,permvec,permuter.isFlipping);
	tstat=calculateTstat(datam,dmperm,tc,cope,varcope,sigmaSquared,dof[0]);    
      }
      if ( clusters.isAveraging )
      { 
	ColumnVector clustersizes;
	tstat4D.setmatrix(tstat,mask);
	tstat4D.binarise(clusterThreshold);
	volume<int> clusterLabels=connected_components(tstat4D[0],clustersizes,CLUST_CON);
	ColumnVector entries,cluster(clustersizes.Nrows());
	cluster=0;
	entries=cluster;
	for(int z=0; z<mask.zsize(); z++)
	  for(int y=0; y<mask.ysize(); y++)
	    for(int x=0; x<mask.xsize(); x++)
	      if (clusterLabels(x,y,z))
	      {
		cluster(clusterLabels(x,y,z))+=temp4D(x,y,z,0);
		entries(clusterLabels(x,y,z))++;
	      }
	clustersizes=SD(clustersizes,SD(cluster,entries));
	clusterNormals.store(clusterLabels,clustersizes,mask,1,perm);
      }
      if ( clusterEnhanced.isAveraging )
      {
	if (!lowram) tstat=previousTFCEStat.Row(perm);
	else tstat=tfce(tstat,mask,tfce_delta,opts.tfce_height.value(),opts.tfce_size.value(),opts.tfce_connectivity.value());
	tstat=SD(tstat,clusterEnhanced.sumStatMat); 
	clusterEnhancedNormals.store(tstat,perm);
      }
    }
  }
   
  //OUTPUT Routines 
  tstat4D.setmatrix(voxels.originalStatistic.Row(1),mask);
  save_volume4D(tstat4D,opts.out_fileroot.value()+statLabel);
  if ( opts.voxelwiseOutput.value() ) OutputStat(voxels,mask,nPerms,statLabel,opts.out_fileroot.value()+"_vox",opts.outputTextNull.value(),false,opts.verbose.value());
  if ( clusterThreshold > 0 ) OutputStat(clusters,mask,nPerms,statLabel,opts.out_fileroot.value()+"_clustere",opts.outputTextNull.value(),opts.outputRaw.value(),opts.verbose.value());
  if ( massThreshold > 0 )    OutputStat(clusterMasses,mask,nPerms,statLabel,opts.out_fileroot.value()+"_clusterm",opts.outputTextNull.value(),opts.outputRaw.value(),opts.verbose.value());
  if ( clusters.isAveraging ) OutputStat(clusterNormals,mask,nPerms,statLabel,opts.out_fileroot.value()+"_clustern",opts.outputTextNull.value(),opts.outputRaw.value(),opts.verbose.value());  
  if ( opts.tfce.value() )    OutputStat(clusterEnhanced,mask,nPerms,statLabel,opts.out_fileroot.value()+"_tfce",opts.outputTextNull.value(),opts.outputRaw.value(),opts.verbose.value());
  if ( clusterEnhanced.isAveraging ) OutputStat(clusterEnhancedNormals,mask,nPerms,statLabel,opts.out_fileroot.value()+"_tfcen",opts.outputTextNull.value(),opts.outputRaw.value(),opts.verbose.value());
  if (opts.outputTextPerm.value()) 
    permuter.writePermutationHistory(opts.out_fileroot.value()+"_perm"+statLabel+".txt");  
}

bool convertContrast(const Matrix& inputModel,const Matrix& inputContrast,const Matrix& inputData,Matrix& outputModel,Matrix& outputContrast, Matrix& outputData, const int mode, const bool debug)
{
    DiagonalMatrix D;
    Matrix U,V,inverseContrast=pinv(inputContrast);
    Matrix W1=inputModel*inverseContrast;
    Matrix W2=inputModel-inputModel*inputContrast.t()*inverseContrast.t();
    int originalConfounds(W2.Ncols());
    Matrix W2tmp(W2.Nrows(),0);
    for(int col=1;col<=W2.Ncols();col++) {
      ColumnVector tmpCol(W2.Column(col));
      if ( !tmpCol.IsZero() )
	W2tmp |= tmpCol;
    }
    W2=W2tmp;
    bool confoundsExist( W2.Ncols() > 0 );
    if ( confoundsExist ) {
      if ( W2.Ncols() <= W2.Nrows() ) 
	SVD(W2,D,U,V);
      else 
	SVD(W2.t(),D,V,U); //Note the swap of U and V for transposed input
      float confoundMax(D.Maximum());
      while ( D(D.Ncols()) < confoundMax*1e-10 ) 
	D = D.SymSubMatrix(1,D.Ncols()-1);
      if ( debug )
	cerr << "Removed " << originalConfounds-D.Ncols() << " null confounds" << endl;
      W2=U.Columns(1,D.Ncols()); //Multiplying by D just preserved scaling, not needed as per MJ.
      SVD(W1,D,U,V);
      float interestMax(D.Maximum());
      if ( (interestMax>confoundMax) && (interestMax == (interestMax+(confoundMax*10.0)))) { //test if confound space is very small compared to interest
	confoundsExist=false;
	W2=0;
      }
    }
    if ( confoundsExist &&  ( mode == 0 || mode == 1 ) ) {
      outputData=(IdentityMatrix(W2.Nrows())-W2*W2.t())*inputData;
      if ( debug )
	cerr << "Orthogonalising data wrt to confounds." << endl;
    }
    else { 
      outputData=inputData;
      if ( debug )
	cerr << "Not orthogonalising data wrt to confounds." << endl;
    }
    outputModel=W1;                  //All modes start with "base" model and contrast of  W1 and I(r)
    outputContrast=IdentityMatrix(inputContrast.Nrows());
    if ( mode == 0 && confoundsExist ) //Kennedy  Regress Y_a on X_a
      outputModel=W1-W2*W2.t()*W1;
    if ( mode == 1 || mode == 2 || mode == 3 ) { //Regress Y_a (Freedman_Lane) or Y (No unconfounding) or Y_aFull ( ter Braak ) on X | Z 
      if ( confoundsExist ) {
	Matrix nuisanceContrast(inputContrast.Nrows(),W2.Ncols());
	nuisanceContrast=0;
	outputContrast |= nuisanceContrast;
	outputModel |= W2;
      }        
      if ( mode == 3 ) 
	outputData=(IdentityMatrix(outputModel.Nrows())-outputModel*pinv(outputModel))*inputData;
      
    }
    return(confoundsExist);
}


void analyseContrast(const Matrix& inputContrast, const Matrix& dm, const Matrix& datam, const volume<float>& mask,const Matrix& gp,const int& contrastNo,ranopts& opts, const Matrix& effectiveDesign)
{
  //-ve num for f-stat contrast
  Matrix NewModel,NewDataM, NewCon;
  VoxelwiseDesign fullVoxelwiseDesign;
  Permuter permuter;
  bool hasConfounds(false);
  if (voxelwiseInput.isSet) {
    fullVoxelwiseDesign.storeByVoxel(datam.Ncols());
    NewDataM.ReSize(datam);
    for(int voxel=1;voxel<=datam.Ncols();voxel++)
    {
      Matrix tempData;
      hasConfounds=convertContrast(voxelwiseInput.adjustDesign(dm,voxel),inputContrast,datam.Column(voxel),fullVoxelwiseDesign.designAtVoxel[voxel],fullVoxelwiseDesign.contrastAtVoxel[voxel],tempData,opts.confoundMethod.value(),opts.isDebugging.value()) || hasConfounds;
      NewDataM.Column(voxel)=tempData;
      fullVoxelwiseDesign.dofAtVoxel[voxel].push_back((int)ols_dof(fullVoxelwiseDesign.designAtVoxel[voxel])-(int)opts.demean_data.value());
    }
    fullVoxelwiseDesign.isSet=true;
    NewModel=fullVoxelwiseDesign.designAtVoxel[1]; //Arbitrary choices
    NewCon=fullVoxelwiseDesign.contrastAtVoxel[1];
  } else hasConfounds=convertContrast(dm,inputContrast,datam,NewModel,NewCon,NewDataM,opts.confoundMethod.value(),opts.isDebugging.value());
  
  if ( opts.isDebugging.value() ) {
    if ( hasConfounds ) 
      cerr << "Confounds detected." << endl;
    else 
      cerr << "No confounds detected." << endl;
  }

  bool oneRegressor( inputContrast.SumAbsoluteValue() == inputContrast.MaximumAbsoluteValue() );
  if ( opts.effectiveDesignFile.value()!="" )
    permuter.createPermutationScheme(effectiveDesign,gp.Column(1),(contrastNo>0 && oneRegressor),opts.n_perm.value(),opts.detectNullSubjects.value(),opts.isDebugging.value(),opts.permuteBlocks.value(),opts.one_samp.value()); 
  else
    permuter.createPermutationScheme(remmean(dm)*inputContrast.t(),gp.Column(1),(contrastNo>0 && oneRegressor),opts.n_perm.value(),opts.detectNullSubjects.value(),opts.isDebugging.value(),opts.permuteBlocks.value(),opts.one_samp.value());
  if( permuter.isFlipping ) cout << "One-sample design detected; sign-flipping instead of permuting." << endl;
  if( opts.verbose.value() || opts.how_many_perms.value() ) 
  {
    if( permuter.isFlipping ) cout << permuter.uniquePermutations[0] << " sign-flips required for exhaustive test";
    else cout << permuter.uniquePermutations[0] << " permutations required for exhaustive test";
    if (contrastNo>0)  cout << " of t-test " << contrastNo << endl;
    if (contrastNo==0) cout << " of all t-tests " << endl;
    if (contrastNo<0)  cout << " of f-test " << abs(contrastNo) << endl;
    if(opts.how_many_perms.value()) return;
  }  
  vector<int> dof(1,(int)ols_dof(NewModel)-(int)opts.demean_data.value()); 
  calculatePermutationStatistics(opts,mask,NewDataM,NewCon,NewModel,contrastNo,dof,permuter,fullVoxelwiseDesign); 
}


void analyseFContrast(Matrix& fc,Matrix& tc,Matrix& model,Matrix& data,volume<float>& mask,Matrix& gp,ranopts& opts, const Matrix& effectiveDesign)
{   
  int startingContrast(1);
  int finalContrast(fc.Nrows());
  if ( opts.skipTo.value() > 0 ) {
    startingContrast=opts.skipTo.value();
    finalContrast=min( fc.Nrows(), opts.skipTo.value() );
  }

    for( int fstat=startingContrast; fstat<=finalContrast ; fstat++ ) {
      Matrix fullFContrast( 0, tc.Ncols() );
      for (int tcon=1; tcon<=fc.Ncols() ; tcon++ )
	if (fc(fstat,tcon)==1) fullFContrast &= tc.Row(tcon);
      analyseContrast(fullFContrast,model,data,mask,gp,-fstat,opts,effectiveDesign);
    }
}

int main(int argc,char *argv[]) {
  Log& logger = LogSingleton::getInstance();
  ranopts& opts = ranopts::getInstance();
  opts.parse_command_line(argc,argv,logger);
  if (opts.parallelData.value()) opts.verbose.set_value("false");
  Matrix model, Tcontrasts, Fcontrasts, data, blockLabels, effectiveDesign;
  volume<float> mask;
  if ( opts.verbose.value() ) { 
    cout << "randomise options: ";
    for (int i=1;i<argc;i++) cout << argv[i] << " ";
    cout << endl;
  }
  try { 
    Initialise(opts,mask,data,Tcontrasts,model,Fcontrasts,blockLabels,effectiveDesign);
    bool needsDemean=true;
    for (int i=1;i<=model.Ncols();i++) if ( fabs( (model.Column(i)).Sum() ) > 0.0001 ) needsDemean=false;
    if (needsDemean && !opts.demean_data.value()) cerr << "Warning: All design columns have zero mean - consider using the -D option to demean your data" << endl;
    if (!needsDemean && opts.demean_data.value()) {
      cerr << "Warning: Data demeaning selected, but at least one design column has non-zero mean - therefore invoking automatic demeaning of design matrix" << endl;
      model=remmean(model);
    }
    if(opts.fc_file.value()!="") analyseFContrast(Fcontrasts,Tcontrasts,model,data,mask,blockLabels,opts,effectiveDesign); 

    int startingContrast(1);
    int finalContrast( Tcontrasts.Nrows() );
    if ( opts.skipTo.value() > 0 ) {
      startingContrast=opts.skipTo.value()-Fcontrasts.Nrows();
      finalContrast=min( Tcontrasts.Nrows(), opts.skipTo.value()-Fcontrasts.Nrows() );
    }


    for (int tstat = startingContrast; tstat <= finalContrast  && tstat > 0 && !opts.doFOnly.value(); tstat++ )  
      analyseContrast(Tcontrasts.Row(tstat),model,data,mask,blockLabels,tstat,opts,effectiveDesign); 
  }
  catch(Exception& e) 
  { 
    cerr << "ERROR: Program failed" <<  e.what() << endl << endl << "Exiting" << endl; 
    return 1;
  }
  catch(...) 
  { 
    cerr << "ERROR: Program failed, unknown exception" << endl << endl << "Exiting" << endl; 
    return 1;
  }
  if ( opts.verbose.value() ) 
    cout << "Finished, exiting." << endl;
  return 0;
}

//Permuter Class
void Permuter::createPermutationScheme(const Matrix& design, ColumnVector groups, const bool oneNonZeroContrast, const long requiredPermutations, const bool detectingNullElements, const bool outputDebug, const bool permuteBlocks, const bool forceFlipping)
{
  nBlocks=int(groups.Maximum())+1; //+1 to include the "0" block
  nSubjects=design.Nrows();
  if (detectingNullElements)
    for(int row=1;row<=nSubjects;row++)
      if (abs(design.Row(row).Sum())<1e-10 && !isFlipping) //original just checked if Sum()==0
	groups(row)=0;
  isPermutingBlocks=permuteBlocks;

  if ( isPermutingBlocks ) {
    blockPermuter=new Permuter;
      ColumnVector dummyBlocks(nBlocks-1);
      dummyBlocks=1;
      Matrix effectiveBlockDesign(0,design.Ncols()*design.Nrows()/(nBlocks-1));
    for(int group=1;group<nBlocks;group++) {
      RowVector currentRow;
      for(int row=1;row<=nSubjects;row++)
	if(groups(row)==group) currentRow |= design.Row(row);
      effectiveBlockDesign &= currentRow;
    }
    blockPermuter->createPermutationScheme(effectiveBlockDesign,dummyBlocks,false,requiredPermutations,false,outputDebug,false,forceFlipping);
  }
  
  ColumnVector labels = createDesignLabels(design|groups);
  if ( forceFlipping )
    labels=1;

  if ( isPermutingBlocks ) 
    for( int row=1;row<=nSubjects;row++ )
      labels(row)=row;

  if ( forceFlipping )
    labels=1;
  isFlipping = ( (labels.Maximum()==1) && oneNonZeroContrast ) || forceFlipping ;

  originalLocation.resize(nBlocks);
  permutedLabels.resize(nBlocks);
  originalLabels.resize(nBlocks);
  for(int group=0;group<nBlocks;group++)
  {
     int member(0);
     for(int row=1;row<=nSubjects;row++)
       if(groups(row)==group) member++;
     originalLocation[group].ReSize(member);
     permutedLabels[group].ReSize(member);
     for(int row=nSubjects;row>=1;row--) //Now work backwards to fill in the starting locations
       if(groups(row)==group) originalLocation[group](member--)=row;
  }
  initialisePermutationBlocks(labels,requiredPermutations);
  if (outputDebug) 
    cerr << "Subject | Design | group | label" << endl << ( truePermutation | design | groups | labels ) << endl;
}

ColumnVector Permuter::returnPreviousTruePermutation(const long permutationNumber)
{
  if (isFlipping) 
    return previousPermutations[permutationNumber-1];
  else {    
    ColumnVector permvec(unpermutedVector);
    for(long perm=1; perm<=permutationNumber; perm++) 
      createTruePermutation(previousPermutations[perm-1],previousPermutations[perm-1-int(perm!=1)],permvec);  
    return permvec;
  }
}

ColumnVector Permuter::returnPreviousTruePermutation(const long permutationNumber, ColumnVector& previousState)
{
  if (isFlipping) 
    return previousPermutations[permutationNumber-1];
  else {    
    createTruePermutation(previousPermutations[permutationNumber-1],previousPermutations[permutationNumber-1-int(permutationNumber!=1)],previousState);  
    return previousState;
  }
}

void Permuter::initialisePermutationBlocks(const ColumnVector& designLabels,const long requiredPermutations)
{
  truePermutation.ReSize(nSubjects);
  for(int i=1;i<=nSubjects;i++) truePermutation(i)=i;
  if (isFlipping) truePermutation=1;
  unpermutedVector=truePermutation;
  uniquePermutations.resize(nBlocks);
  uniquePermutations[0]=1;
  for(int group=0;group<nBlocks;group++)
  { 
    for(int row=1;row<=permutedLabels[group].Nrows();row++) 
      permutedLabels[group](row)=designLabels((int)originalLocation[group](row));
    if (group>0) uniquePermutations[group]=computeUniquePermutations(permutedLabels[group],isFlipping);
    uniquePermutations[0]*=uniquePermutations[group];
    originalLabels[group]=permutedLabels[group];
  }

  if ( isPermutingBlocks )
    uniquePermutations[0]=blockPermuter->uniquePermutations[0];
  isRandom=!(requiredPermutations==0 || requiredPermutations>=uniquePermutations[0]);
  if (isRandom) finalPermutation=requiredPermutations;
  else finalPermutation=uniquePermutations[0];
  previousPermutations.reserve((long)finalPermutation);
}

ColumnVector Permuter::permutationVector()
{
ColumnVector newvec(nSubjects); 
   for(int block=0;block<nBlocks;block++)
     for(int row=1;row<=permutedLabels[block].Nrows();row++) 
       newvec( (int)originalLocation[block](row) ) = permutedLabels[block](row);
   return newvec;
}

void Permuter::createTruePermutation(const ColumnVector& newLabels,ColumnVector copyOldLabels,ColumnVector& permvec)
{
  if (isFlipping) permvec=permutationVector();
  else
  {
    for(int k=1;k<=newLabels.Nrows();k++)	       
      if(newLabels(k)!=copyOldLabels(k))
        for(int l=1;l<=newLabels.Nrows();l++) 
	  if(newLabels(l)!=copyOldLabels(l) && copyOldLabels(l)==newLabels(k) )
 	   {
	     swap(permvec(l),permvec(k));
	     swap(copyOldLabels(l),copyOldLabels(k));
           }
  }
}

ColumnVector Permuter::nextPermutation(const long permutationNumber)
{
  if ( isPermutingBlocks ) {
    ColumnVector permutedBlocksFoo=blockPermuter->nextPermutation(permutationNumber,false,false);
    for(int block=1;block<nBlocks;block++) {
      if ( blockPermuter->isFlipping )
	permutedLabels[block]=originalLabels[block]*permutedBlocksFoo(block);
      else
	permutedLabels[block]=originalLabels[(int)permutedBlocksFoo(block)];
    }
    return(permutationVector());
  }

  for(int group=1;group<nBlocks;group++)
  {
      if(isFlipping) nextFlip(permutedLabels[group]);
      else nextShuffle(permutedLabels[group]);
      if (!isRandom && permutedLabels[group]!=originalLabels[group] ) //Move to next group as either "reset" has occurred or we are in random mode
	break;
  }
  return(permutationVector());
}

ColumnVector Permuter::nextPermutation(const long permutationNumber, const bool printStatus, const bool isStoring)
{
  if (permutationNumber!=1 && printStatus) cout << "Starting permutation " << permutationNumber << endl;
  else if (printStatus) cout << "Starting permutation " << permutationNumber << " (Unpermuted data)" << endl;

  ColumnVector currentLabels=permutationVector();
  ColumnVector newPermutation;
  do
  {
    if (permutationNumber!=1) newPermutation=nextPermutation(permutationNumber);
  } while(isRandom && isPreviousPermutation(newPermutation));
  if(isStoring || isRandom) previousPermutations.push_back(permutationVector());
  createTruePermutation(permutationVector(),currentLabels,truePermutation);				       
  return(truePermutation);
}

bool Permuter::isPreviousPermutation(const ColumnVector& newPermutation){
  for(int i=previousPermutations.size()-1; i>=0; i--)
    if(newPermutation==previousPermutations[i]) return true;
  return false;
  }

void Permuter::nextShuffle(ColumnVector& perm){
   vector<int> temp;
   for (int i=1;i<=perm.Nrows();i++) temp.push_back((int)perm(i));
   if (isRandom) random_shuffle(temp.begin(),temp.end());
   else next_permutation(temp.begin(),temp.end());
   for (int i=1;i<=perm.Nrows();i++) perm(i)=temp[i-1];
}

void Permuter::nextFlip(ColumnVector& flip){
       
  if (isRandom)
  {
    for(int i=1;i<=flip.Nrows();i++)
    {
      float tmp=(float)rand()/RAND_MAX;
      if(tmp > 0.5) flip(i)=1;
      else  flip(i)=-1;
    }     
  }
  else for (int n=flip.Nrows();n>0;n--)
    if(flip(n)==1)
    {
      flip(n)=-1;
      if (n<flip.Nrows()) flip.Rows(n+1,flip.Nrows())=1;
      return;
    } 
}

double Permuter::computeUniquePermutations(const ColumnVector& labels,const bool calculateFlips){
  if (calculateFlips) return std::pow(2.0,labels.Nrows());
  ColumnVector label_counts((int)labels.MaximumAbsoluteValue());
  label_counts=0;
  for(int i=1; i<=labels.Nrows(); i++) label_counts(int(labels(i)))++;
  double yo = lgam(labels.Nrows()+1);
  for(int i=1; i<=labels.MaximumAbsoluteValue(); i++)
    yo -= lgam(label_counts(i)+1);
  return std::floor(exp(yo)+0.5);
}

ColumnVector Permuter::createDesignLabels(const Matrix& design){
  ColumnVector designLabels(design.Nrows());
  vector<RowVector> knownLabels;
  for(int i=1;i<=design.Nrows();i++){
    bool wasExistingLabel(false);
    for(unsigned int l=0;l<knownLabels.size();l++){
      if(design.Row(i)==knownLabels[l]){
	designLabels(i)=l+1;
	wasExistingLabel=true;
      }
    }
    if(!wasExistingLabel){
      knownLabels.push_back(design.Row(i));
      designLabels(i)=knownLabels.size();
    }
  }
  return(designLabels);
}

double Permuter::reportRequiredPermutations(const bool printToScreen)
{
  if (printToScreen)
  {
    if(isRandom) cout<<"Doing " << finalPermutation << " random permutations"<<endl;
    else cout<<"Doing all "<< finalPermutation <<" unique permutations"<<endl;
  }
  return(finalPermutation);
}

void Permuter::writePermutationHistory(const string& fileName)
{
  ofstream output_file(fileName.c_str());
  ColumnVector state(unpermutedVector);
  for(unsigned long perm=1; perm<=finalPermutation; perm++)
    output_file << returnPreviousTruePermutation(perm,state).t(); 
  output_file.close();
}

Permuter::Permuter()
{
  isPermutingBlocks=false;
  blockPermuter=NULL;
}

Permuter::~Permuter()
{
  if ( blockPermuter != NULL )
     delete blockPermuter;
}
