//     
//     fslmaths.cc Image processing routines, some basic, some not so basic...
//     Steve Smith, David Flitney, Mark Jenkinson, Stuart Clare, Thomas Nichols and Matthew Webster, FMRIB Image Analysis Group
//     Copyright (C) 2000-2008 University of Oxford  
//     

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

//     
#define EXPOSE_TREACHEROUS
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "utils/fsl_isfinite.h"
#include "libprob.h"

using namespace MISCMATHS;
using namespace NEWIMAGE;

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

int printUsage(const string& programName) 
{
  cout << "\nUsage: fslmaths [-dt <datatype>] <first_input> [operations and inputs] <output> [-odt <datatype>]" << endl;

  cout << "\nDatatype information:" << endl;
  cout << " -dt sets the datatype used internally for calculations (default float for all except double images)" << endl;
  cout << " -odt sets the output datatype ( default is float )" << endl;
  cout << " Possible datatypes are: char short int float double input" << endl;
  cout << " \"input\" will set the datatype to that of the original image" << endl;

  cout << "\nBinary operations:" << endl;
  cout << "  (some inputs can be either an image or a number)" << endl;
  cout << " -add   : add following input to current image" << endl;
  cout << " -sub   : subtract following input from current image" << endl;
  cout << " -mul   : multiply current image by following input" << endl;
  cout << " -div   : divide current image by following input" << endl;
  cout << " -rem   : modulus remainder - divide current image by following input and take remainder" << endl;
  cout << " -mas   : use (following image>0) to mask current image" << endl;
  cout << " -thr   : use following number to threshold current image (zero anything below the number)" << endl;
  cout << " -thrp  : use following percentage (0-100) of ROBUST RANGE to threshold current image (zero anything below the number)" << endl;
  cout << " -thrP  : use following percentage (0-100) of ROBUST RANGE of non-zero voxels and threshold below" << endl;
  cout << " -uthr  : use following number to upper-threshold current image (zero anything above the number)" << endl;
  cout << " -uthrp : use following percentage (0-100) of ROBUST RANGE to upper-threshold current image (zero anything above the number)" << endl;
  cout << " -uthrP : use following percentage (0-100) of ROBUST RANGE of non-zero voxels and threshold above" << endl;
  cout << " -max   : take maximum of following input and current image" << endl;
  cout << " -min   : take minimum of following input and current image" << endl;
  cout << " -seed  : seed random number generator with following number" << endl;
  cout << " -restart : replace the current image with input for future processing operations" << endl;
  cout << " -save : save the current working image to the input filename" << endl;

  cout << "\nBasic unary operations:" << endl;
  cout << " -exp   : exponential" << endl;
  cout << " -log   : natural logarithm" << endl;
  cout << " -sin   : sine function" << endl;
  cout << " -cos   : cosine function" << endl;
  cout << " -tan   : tangent function" << endl;
  cout << " -asin  : arc sine function" << endl;
  cout << " -acos  : arc cosine function" << endl;
  cout << " -atan  : arc tangent function" << endl;
  cout << " -sqr   : square" << endl;
  cout << " -sqrt  : square root" << endl;
  cout << " -recip : reciprocal (1/current image)" << endl;
  cout << " -abs   : absolute value" << endl;
  cout << " -bin   : use (current image>0) to binarise" << endl;
  cout << " -binv  : binarise and invert (binarisation and logical inversion)" << endl;
  cout << " -fillh : fill holes in a binary mask (holes are internal - i.e. do not touch the edge of the FOV)" << endl;
  cout << " -fillh26 : fill holes using 26 connectivity" << endl;
  cout << " -index : replace each nonzero voxel with a unique (subject to wrapping) index number" << endl;
  cout << " -grid <value> <spacing> : add a 3D grid of intensity <value> with grid spacing <spacing>" << endl;
  cout << " -edge  : edge strength" << endl;
  cout << " -tfce <H> <E> <connectivity>: enhance with TFCE, e.g. -tfce 2 0.5 6 (maybe change 6 to 26 for skeletons)" << endl;
  cout << " -tfceS <H> <E> <connectivity> <X> <Y> <Z> <tfce_thresh>: show support area for voxel (X,Y,Z)" << endl;
  cout << " -nan   : replace NaNs (improper numbers) with 0" << endl;
  cout << " -nanm  : make NaN (improper number) mask with 1 for NaN voxels, 0 otherwise" << endl;
  cout << " -rand  : add uniform noise (range 0:1)" << endl;
  cout << " -randn : add Gaussian noise (mean=0 sigma=1)" << endl;
  cout << " -inm <mean> :  (-i i ip.c) intensity normalisation (per 3D volume mean)" << endl;
  cout << " -ing <mean> :  (-I i ip.c) intensity normalisation, global 4D mean)" << endl;
  cout << " -range : set the output calmin/max to full data range" << endl;

  cout << "\nMatrix operations:" << endl;
  cout << " -tensor_decomp : convert a 4D (6-timepoint )tensor image into L1,2,3,FA,MD,MO,V1,2,3 (remaining image in pipeline is FA)" << endl;

  cout << "\nKernel operations (set BEFORE filtering operation if desired):" << endl;
  cout << " -kernel 3D : 3x3x3 box centered on target voxel (set as default kernel)" << endl;
  cout << " -kernel 2D : 3x3x1 box centered on target voxel" << endl;
  cout << " -kernel box    <size>     : all voxels in a cube of width <size> mm centered on target voxel" << endl;
  cout << " -kernel boxv   <size>     : all voxels in a cube of width <size> voxels centered on target voxel, CAUTION: size should be an odd number" << endl;
  cout << " -kernel boxv3  <X> <Y> <Z>: all voxels in a cuboid of dimensions X x Y x Z centered on target voxel, CAUTION: size should be an odd number" << endl;
  cout << " -kernel gauss  <sigma>    : gaussian kernel (sigma in mm, not voxels)" << endl;
  cout << " -kernel sphere <size>     : all voxels in a sphere of radius <size> mm centered on target voxel" << endl;
  cout << " -kernel file   <filename> : use external file as kernel" << endl;

  cout << "\nSpatial Filtering operations: N.B. all options apart from -s use the default kernel or that _previously_ specified by -kernel" << endl;
  cout << " -dilM    : Mean Dilation of non-zero voxels" << endl;
  cout << " -dilD    : Modal Dilation of non-zero voxels" << endl;
  cout << " -dilF    : Maximum filtering of all voxels" << endl;
  cout << " -dilall  : Apply -dilM repeatedly until the entire FOV is covered" << endl;
  cout << " -ero     : Erode by zeroing non-zero voxels when zero voxels found in kernel" << endl;
  cout << " -eroF    : Minimum filtering of all voxels" << endl;
  cout << " -fmedian : Median Filtering " << endl;
  cout << " -fmean   : Mean filtering, kernel weighted (conventionally used with gauss kernel)" << endl;
  cout << " -fmeanu  : Mean filtering, kernel weighted, un-normalised (gives edge effects)" << endl;
  cout << " -s <sigma> : create a gauss kernel of sigma mm and perform mean filtering" << endl;
  cout << " -subsamp2  : downsamples image by a factor of 2 (keeping new voxels centred on old)" << endl;
  cout << " -subsamp2offc  : downsamples image by a factor of 2 (non-centred)" << endl;

  cout << "\nDimensionality reduction operations:" << endl;
  cout << "  (the \"T\" can be replaced by X, Y or Z to collapse across a different dimension)" << endl;
  cout << " -Tmean   : mean across time" << endl;
  cout << " -Tstd    : standard deviation across time" << endl;
  cout << " -Tmax    : max across time" << endl;
  cout << " -Tmaxn   : time index of max across time" << endl;
  cout << " -Tmin    : min across time" << endl;
  cout << " -Tmedian : median across time" << endl;
  cout << " -Tperc <percentage> : nth percentile (0-100) of FULL RANGE across time" << endl;
  cout << " -Tar1    : temporal AR(1) coefficient (use -odt float and probably demean first)" << endl;

  cout << "\nBasic statistical operations:" << endl;
  cout << " -pval    : Nonparametric uncorrected P-value, assuming timepoints are the permutations; first timepoint is actual (unpermuted) stats image" << endl;
  cout << " -pval0   : Same as -pval, but treat zeros as missing data" << endl;
  cout << " -cpval   : Same as -pval, but gives FWE corrected P-values" << endl;
  cout << " -ztop    : Convert Z-stat to (uncorrected) P" << endl;
  cout << " -ptoz    : Convert (uncorrected) P to Z" << endl;
  cout << " -rank    : Convert data to ranks (over T dim)" << endl;
  cout << " -ranknorm: Transform to Normal dist via ranks" << endl;

  cout << "\nMulti-argument operations:" << endl;
  cout << " -roi <xmin> <xsize> <ymin> <ysize> <zmin> <zsize> <tmin> <tsize> : zero outside roi (using voxel coordinates). Inputting -1 for a size will set it to the full image extent for that dimension." << endl;
  cout << " -bptf  <hp_sigma> <lp_sigma> : (-t in ip.c) Bandpass temporal filtering; nonlinear highpass and Gaussian linear lowpass (with sigmas in volumes, not seconds); set either sigma<0 to skip that filter" << endl;
  cout << " -roc <AROC-thresh> <outfile> [4Dnoiseonly] <truth> : take (normally binary) truth and test current image in ROC analysis against truth. <AROC-thresh> is usually 0.05 and is limit of Area-under-ROC measure FP axis. <outfile> is a text file of the ROC curve (triplets of values: FP TP threshold). If the truth image contains negative voxels these get excluded from all calculations. If <AROC-thresh> is positive then the [4Dnoiseonly] option needs to be set, and the FP rate is determined from this noise-only data, and is set to be the fraction of timepoints where any FP (anywhere) is seen, as found in the noise-only 4d-dataset. This is then controlling the FWE rate. If <AROC-thresh> is negative the FP rate is calculated from the zero-value parts of the <truth> image, this time averaging voxelwise FP rate over all timepoints. In both cases the TP rate is the average fraction of truth=positive voxels correctly found." << endl;

  cout << "\nCombining 4D and 3D images:" << endl;
  cout << " If you apply a Binary operation (one that takes the current image and a new image together), when one is 3D and the other is 4D," << endl;
  cout << " the 3D image is cloned temporally to match the temporal dimensions of the 4D image." << endl;

  cout << "\ne.g. fslmaths inputVolume -add inputVolume2 output_volume" << endl;
  cout << "     fslmaths inputVolume -add 2.5 output_volume" << endl;
  cout << "     fslmaths inputVolume -add 2.5 -mul inputVolume2 output_volume\n" << endl;
  cout << "     fslmaths 4D_inputVolume -Tmean -mul -1 -add 4D_inputVolume demeaned_4D_inputVolume\n" << endl;
  return 1;
}

template <class T>
void loadNewImage(volume4D<T> &oldI, volume4D<T> &newI, string filename)
{
  read_volume4D(newI, filename);
  if ( (oldI.tsize() == 1) && (oldI.tsize() < newI.tsize()) )
    {
      volume4D<T> tmpvol=oldI;
      oldI.reinitialize(oldI.xsize(),oldI.ysize(),oldI.zsize(),newI.tsize()); 
      oldI.copyproperties(tmpvol);
      for (int t=0;t<newI.tsize();t++) oldI[t]=tmpvol[0]; 
    }

  // sanity check on valid orientation info (it should be consistent)
  if ((newI.sform_code()!=NIFTI_XFORM_UNKNOWN) || (newI.qform_code()!=NIFTI_XFORM_UNKNOWN))
    {
      if ((oldI.sform_code()!=NIFTI_XFORM_UNKNOWN) || (oldI.qform_code()!=NIFTI_XFORM_UNKNOWN))
	{
	  float rms = rms_deviation(newI.newimagevox2mm_mat(),oldI.newimagevox2mm_mat());
	  if (rms>0.5) {   // arbitrary 0.5mm rms diff threshold - maybe too sensitive?
	    cerr << endl << "WARNING:: Inconsistent orientations for individual images in pipeline!" <<endl;  
	    cerr <<"          Will use voxel-based orientation which is probably incorrect - *PLEASE CHECK*!" <<endl<<endl; }
	}
    }
}

int check_for_output_name(int i, int argc_1)
{
  if (i>argc_1) {
    cerr << "Error: no output filename specified!" << endl;
    exit(EXIT_FAILURE);
  }
  return 0;
}


template <class T>
int inputParser(int argc, char *argv[], short output_dt)
{
  volume4D<T> inputVolume;
  volume<float> kernel(box_kernel(3,3,3));
  bool separable(false);
  read_volume4D(inputVolume,string(argv[1]));
  bool modifiedInput(false);
  bool setDisplayRange(false);
  float tfce_delta(0);
  float tfce_minT(0);

  int i=2;
  for (i = 2; i < argc-1; i++)  //main loop
  {    
    volume4D<T> temp_volume;
    modifiedInput=true;
    /***************************************************************/
    /******************** Dimensionality Reduction *****************/
    /***************************************************************/
    if (isupper((int)argv[i][1]) && argv[i][0] == '-')  //if first letters are -capital - dimensionality reduction...
    { 
      int xoff=1,yoff=1,zoff=1,toff=1,nsize;
      if (argv[i][1] == 'T') toff=inputVolume.tsize(); 
      if (argv[i][1] == 'Z') zoff=inputVolume.zsize();  
      if (argv[i][1] == 'Y') yoff=inputVolume.ysize();  
      if (argv[i][1] == 'X') xoff=inputVolume.xsize();  
      temp_volume=inputVolume;
      inputVolume.reinitialize(inputVolume.xsize()/xoff,inputVolume.ysize()/yoff,inputVolume.zsize()/zoff,inputVolume.tsize()/toff); 
      inputVolume.copyproperties(temp_volume);
      nsize=xoff*yoff*zoff*toff;
      volume<T> column_volume(nsize,1,1); //will be size of appropriate dimension, as only 1 arg is non-unitary
      for(int t=0;t<inputVolume.tsize();t++)           
        for(int z=0;z<inputVolume.zsize();z++)
          for(int y=0;y<inputVolume.ysize();y++)	    
	    for(int x=0;x<inputVolume.xsize();x++)
	    {
              for (int j=0;j<nsize;j++) column_volume.value(j,0,0)=temp_volume(x+j*(xoff!=1),y+j*(yoff!=1),z+j*(zoff!=1),t+j*(toff!=1));
              //This goes along the appropriate axis (non unitary offset variable) and fills a "column" volume with data
	      if (string(argv[i]+2) == "max")    inputVolume.value(x,y,z,t)=column_volume.max(); 
	      else if (string(argv[i]+2) == "min")    inputVolume.value(x,y,z,t)=column_volume.min();
              else if (string(argv[i]+2) == "mean")   inputVolume.value(x,y,z,t)=(T)column_volume.mean();
              else if (string(argv[i]+2) == "std")    inputVolume.value(x,y,z,t)=(T)column_volume.stddev();
	      else if (string(argv[i]+2) == "maxn")   inputVolume.value(x,y,z,t)=column_volume.maxcoordx();
              else if (string(argv[i]+2) == "median") inputVolume.value(x,y,z,t)=column_volume.percentile(0.5);
	      else if (string(argv[i]+2) == "perc")   inputVolume.value(x,y,z,t)=column_volume.percentile(atof(argv[i+1])/100.0);
              else if (string(argv[i]+2) == "ar1") {
                column_volume-=(T)column_volume.mean();
                double sumsq=column_volume.sumsquares();
                inputVolume(x,y,z,t)=0;
		if(sumsq!=0) for (int k=1;k<nsize;k++) inputVolume(x,y,z,t)+=(T)(column_volume(k,0,0)*column_volume(k-1,0,0)/sumsq);
	      }
	      else {
		cerr << "Error unknown Dimensionality operation: " << string(argv[i]) << endl;
		return(1);
	      }
	    }
       if (string(argv[i]+2) == "perc") i++;
    }

    /***************************************************************/
    /******************** ROC **************************************/
    /***************************************************************/
    else if (string(argv[i])=="-roc")
    {
      // {{{ variables

float aroc_thresh = atof(argv[++i]);

ofstream ofs(argv[++i]);

int border=5;

int separatenoise=1;
if (aroc_thresh<0)
  {
    separatenoise=0;
    aroc_thresh*=-1;
  }

// }}}
      // {{{ setup minval, read pure-noise and log-transform it, log-transform input

float minval=inputVolume.min();

volume4D<float> noise;
if (separatenoise)
  {
    read_volume4D(noise,string(argv[++i]));
    minval=min(noise.min(),minval);
    for(int t=0;t<inputVolume.tsize();t++)
      for(int z=0;z<inputVolume.zsize();z++)
	for(int y=0;y<inputVolume.ysize();y++)	    
	  for(int x=0;x<inputVolume.xsize();x++)
	    noise.value(x,y,z,t)=log(noise.value(x,y,z,t)-minval+1);
  }

volume4D<float> loginput(inputVolume.xsize(),inputVolume.ysize(),inputVolume.zsize(),inputVolume.tsize());
for(int t=0;t<inputVolume.tsize();t++)
  for(int z=0;z<inputVolume.zsize();z++)
    for(int y=0;y<inputVolume.ysize();y++)	    
      for(int x=0;x<inputVolume.xsize();x++)
	loginput.value(x,y,z,t)=log(inputVolume.value(x,y,z,t)-minval+1);

// }}}
      // {{{ get FP from separate noise image

{
  float aroc=0, FP=0, TP=0, TPprev=0, FPprev=0; // note that all the aroc calculations could be removed and replaced at the end with:   mean(arocT).AsScalar()/aroc_thresh
  ColumnVector TPvals(loginput.tsize()), TPprevvals(loginput.tsize()), FPuncorrected(loginput.tsize()), arocT(loginput.tsize());
  TPvals=0; TPprevvals=0; FPuncorrected=0; arocT=0;
  float maxlogval=loginput.max();
  if (separatenoise)
    maxlogval=max(maxlogval,noise.max());
  float delta=maxlogval/1000;
  if (delta==0) delta=1;
  //cout << "minval=" << minval << " maxlogval=" << maxlogval << " delta=" << delta << endl;

  // {{{ setup truth and invtruth etc, images

// read truth
volume<float> truth;
read_volume(truth,string(argv[++i]));

// mask
volume<float> maskim;
if (!separatenoise)
  {
    maskim=truth;
    maskim.binarise(-0.001);
  }

// truth
int invtruecount=0;
for(int z=0;z<loginput.zsize();z++) for(int y=0;y<loginput.ysize();y++) for(int x=0;x<loginput.xsize();x++)
  {
    if ((x<border)||(x>=loginput.xsize()-border)||(y<border)||(y>=loginput.ysize()-border)||(z<border)||(z>=loginput.zsize()-border))
      truth.value(x,y,z)=0;
    else
      invtruecount++;
  }
truth.binarise(truth.max() * 0.05);
int truecount=(int)truth.sum();

// TPim
volume4D<float> TPim;
copyconvert(loginput,TPim);
for(int t=0;t<loginput.tsize();t++)           
  TPim[t]*=truth;

// noise
if (!separatenoise)
  {
    volume<float> invtruth = ((truth*-1.0)+1)*maskim;
    for(int z=0;z<loginput.zsize();z++) for(int y=0;y<loginput.ysize();y++) for(int x=0;x<loginput.xsize();x++)
      if ((x<border)||(x>=loginput.xsize()-border)||(y<border)||(y>=loginput.ysize()-border)||(z<border)||(z>=loginput.zsize()-border))
	invtruth.value(x,y,z)=0;
    invtruecount=(int)invtruth.sum();
    copyconvert(loginput,noise);
    for(int t=0;t<loginput.tsize();t++)           
      noise[t]*=invtruth;
  }

// }}}
					  
  for (float thresh=maxlogval+delta; thresh>0 && FP<aroc_thresh; thresh-=delta)
    {
      TP=FP=0;
      float FPfwe=0;
      for(int t=0;t<loginput.tsize();t++)           
	{
	  float sigTP=0, sigFP=0;
	  for(int z=border;z<loginput.zsize()-border;z++)
	    for(int y=border;y<loginput.ysize()-border;y++)	    
	      for(int x=border;x<loginput.xsize()-border;x++)
		{
		  sigTP+=TPim.value(x,y,z,t)>=thresh;
		  sigFP+=noise.value(x,y,z,t)>=thresh;
		}
	  TPvals(t+1) = sigTP / truecount;
	  FPuncorrected(t+1) = sigFP / invtruecount;
	  FPfwe += sigFP>0;
	  if (!separatenoise)
	    ofs << FPuncorrected(t+1) << " " << TPvals(t+1) << " " << exp(thresh)-1+minval << endl;
	}
	      
      TP = mean(TPvals).AsScalar();
      if (separatenoise)
	FP = FPfwe / loginput.tsize();
      else
	FP = mean(FPuncorrected).AsScalar();

      if (FP<aroc_thresh)
	{
	  aroc += (TP+TPprev)/2 * (FP-FPprev);
	  arocT += (TPvals+TPprevvals)/2 * (FP-FPprev);
	}
      else  // case for when the latest update straddles the FP threshold
	{
	  aroc += (TPprev + (TP-TPprev)*0.5*((aroc_thresh-FPprev)/(FP-FPprev))) * (aroc_thresh-FPprev);
	  arocT += (TPprevvals + (TPvals-TPprevvals)*0.5*((aroc_thresh-FPprev)/(FP-FPprev))) * (aroc_thresh-FPprev);
	}
	      
      TPprev=TP;
      TPprevvals=TPvals;
      FPprev=FP;
      
      if (separatenoise)
	ofs << FP << " " << TP << " " << sqrt(var(FPuncorrected).AsScalar()) << " " << sqrt(var(TPvals).AsScalar()) << " " << exp(thresh)-1+minval << endl;
    }
  
  if (FP<aroc_thresh) // deal with case of when FP never reached aroc_thresh
    {
      aroc += TP*(aroc_thresh-FP);
      arocT += TPvals*(aroc_thresh-FP);
    }
  
  ofs.close();
  cout << aroc / aroc_thresh << " " << quantile(arocT,1) / aroc_thresh << " " << quantile(arocT,3) / aroc_thresh << endl;
}

// }}}
    }


    /***************************************************************/
    /********************Binary Operations**************************/
    /***************************************************************/
    else if (string(argv[i])=="-mas")
    {  
	loadNewImage(inputVolume, temp_volume, string(argv[++i]));
        temp_volume.binarise(0,temp_volume.max()+1,exclusive); // needed to binarise max() value + 1 due to
        for (int t=0;t<inputVolume.tsize();t++)                     //potential issue with exclusive binarise
          inputVolume[t]*=temp_volume[t%temp_volume.tsize()]; //this gives compatibility with 3 and 4D masks
    }                                                            //without needing to have a specific volume3D variable
    /***************************************************************/
    else if (string(argv[i])=="-add"){
      i++;
      if (isNumber(string(argv[i]))) inputVolume+=(T)atof(argv[i]); 
      else 
      {  
	loadNewImage(inputVolume, temp_volume, string(argv[i]));
        for (int t=0;t<inputVolume.tsize();t++) inputVolume[t]+=temp_volume[t%temp_volume.tsize()]; 
      }}
    /***************************************************************/
    else if (string(argv[i])=="-sub"){
      i++;
      if (isNumber(string(argv[i]))) inputVolume-=(T)atof(argv[i]);
      else 
      {  
	loadNewImage(inputVolume, temp_volume, string(argv[i]));
        for (int t=0;t<inputVolume.tsize();t++) inputVolume[t]-=temp_volume[t%temp_volume.tsize()]; 
      }}
    /***************************************************************/
    else if (string(argv[i])=="-mul"){
      i++;
      if (isNumber(string(argv[i]))) inputVolume*=(T)atof(argv[i]);
      else  
      {  
	loadNewImage(inputVolume, temp_volume, string(argv[i]));
        for (int t=0;t<inputVolume.tsize();t++) inputVolume[t]*=temp_volume[t%temp_volume.tsize()]; 
      }}
    /***************************************************************/
    else if (string(argv[i])=="-rem"){
      i++;
      if (isNumber(string(argv[i]))) 
      {
	int denom=(int)atof(argv[i]);
	for(int t=0;t<inputVolume.tsize();t++) 
          for(int z=0;z<inputVolume.zsize();z++)
	    for(int y=0;y<inputVolume.ysize();y++)	    
	      for(int x=0;x<inputVolume.xsize();x++)
		inputVolume(x,y,z,t)=(int)inputVolume(x,y,z,t)%denom;
      }
      else 
      {  
	loadNewImage(inputVolume, temp_volume, string(argv[i]));
        for(int t=0;t<inputVolume.tsize();t++) 
	  { 
            int t2=t%temp_volume.tsize();     
            for(int z=0;z<inputVolume.zsize();z++)
	      for(int y=0;y<inputVolume.ysize();y++)	    
	        for(int x=0;x<inputVolume.xsize();x++)
		  if(temp_volume(x,y,z,t2)!=0) inputVolume(x,y,z,t)=(int)inputVolume(x,y,z,t)%(int)temp_volume(x,y,z,t2); 
	  }
      }
    }
    /***************************************************************/
    else if (string(argv[i])=="-div"){
      i++;
      if (isNumber(string(argv[i]))) {if (atof(argv[i])!=0) inputVolume/=(T)atof(argv[i]);}
      else 
      {  
        loadNewImage(inputVolume, temp_volume, string(argv[i]));
        for(int t=0;t<inputVolume.tsize();t++)      
	{
          int t2=t%temp_volume.tsize();     
          for(int z=0;z<inputVolume.zsize();z++)
            for(int y=0;y<inputVolume.ysize();y++)	    
	      for(int x=0;x<inputVolume.xsize();x++)
		  if(temp_volume(x,y,z,t2)!=0) inputVolume.value(x,y,z,t) /= temp_volume.value(x,y,z,t2);
                  else inputVolume.value(x,y,z,t)=(T)0.0;
        }          
      }}
    /***************************************************************/
    else if (string(argv[i])=="-max" || string(argv[i])=="-min")
    {
      T param=0;
      bool max=false;
      bool file=false;
      if (string(argv[i])=="-max") max=true;
      i++;
      if (isNumber(string(argv[i]))) param=(T)atof(argv[i]);
      else 
      {                           
	loadNewImage(inputVolume, temp_volume, string(argv[i]));
        file=true;
      }     
      for(int t=0;t<inputVolume.tsize();t++)           
        for(int z=0;z<inputVolume.zsize();z++)
          for(int y=0;y<inputVolume.ysize();y++)	    
	    for(int x=0;x<inputVolume.xsize();x++)
	    {
              if (max && file) inputVolume.value(x,y,z,t)=MAX(inputVolume.value(x,y,z,t),temp_volume.value(x,y,z,t%temp_volume.tsize()));
              if (max && !file) inputVolume.value(x,y,z,t)=MAX(inputVolume.value(x,y,z,t),param);
              if (!max && file) inputVolume.value(x,y,z,t)=MIN(inputVolume.value(x,y,z,t),temp_volume.value(x,y,z,t%temp_volume.tsize()));
              if (!max && !file) inputVolume.value(x,y,z,t)=MIN(inputVolume.value(x,y,z,t),param);
            }
    }    
    /***************************************************************/
    else if (string(argv[i])=="-seed") {
      srand(atoi(argv[++i]));
      sdrand(rand(),rand(),rand());
    }
    /***************************************************************/
    else if (string(argv[i])=="-restart")
    {  
      read_volume4D(inputVolume, string(argv[++i]) );
    }   
    /***************************************************************/
    else if (string(argv[i])=="-save")
    {  
      save_volume4D(inputVolume, string(argv[++i]) );
    }   
    /***************************************************************/
    /******************** Unary Operations *************************/
    /***************************************************************/
    else if (string(argv[i])=="-thr") inputVolume.threshold((T)atof(argv[++i]),inputVolume.max()+1,inclusive);
    /***************************************************************/
    else if (string(argv[i])=="-range") setDisplayRange=true;
    /***************************************************************/
    else if (string(argv[i])=="-thrp") 
    {
      T lowerlimit =(T)(inputVolume.robustmin()+(atof(argv[++i])/100.0)*(inputVolume.robustmax()-inputVolume.robustmin())); 
      inputVolume.threshold(lowerlimit,inputVolume.max()+1,inclusive);
    }
    /***************************************************************/
    else if (string(argv[i])=="-thrP") 
    {
      volume4D<T> mask(inputVolume);
      mask.binarise(0,inputVolume.max()+1,exclusive);
      T lowerlimit =(T)(inputVolume.robustmin(mask)+(atof(argv[++i])/100.0)*(inputVolume.robustmax(mask)-inputVolume.robustmin(mask))); 
      inputVolume.threshold(lowerlimit,inputVolume.max()+1,inclusive);
    }
    /***************************************************************/
    else if (string(argv[i])=="-uthr") inputVolume.threshold(inputVolume.min()-1,(T)atof(argv[++i]),inclusive);
    /***************************************************************/
    else if (string(argv[i])=="-uthrp") 
    {
      T upperlimit = (T)(inputVolume.robustmin()+(atof(argv[++i])/100.0)*(inputVolume.robustmax()-inputVolume.robustmin())); 
       inputVolume.threshold(inputVolume.min()-1,upperlimit,inclusive);
    }
    /***************************************************************/
    else if (string(argv[i])=="-uthrP") 
    {
       volume4D<T> mask(inputVolume);
       mask.binarise(0,inputVolume.max()+1,exclusive);
       T upperlimit = (T)(inputVolume.robustmin(mask)+(atof(argv[++i])/100.0)*(inputVolume.robustmax(mask)-inputVolume.robustmin(mask))); 
       inputVolume.threshold(inputVolume.min()-1,upperlimit,inclusive);
    }
    /***************************************************************/
    else if (string(argv[i])=="-kernel") 
    {
       kernel.destroy();
       float xdim=inputVolume.xdim();
       float ydim=inputVolume.ydim();
       float zdim=inputVolume.zdim();
       if(string(argv[i+1])=="2D")      kernel=box_kernel(3,3,1);
       else if(string(argv[i+1])=="3D") kernel=box_kernel(3,3,3);
       else
       {
	 float size=atof(argv[i+2]);
	 separable=false;
         if(string(argv[i+1])=="box" || string(argv[i+1])=="boxv" || string(argv[i+1])=="boxv3" || string(argv[i+1])=="gauss") separable=true;       

	 if(string(argv[i+1])=="box")         kernel=box_kernel(size,xdim,ydim,zdim);
	 else if(string(argv[i+1])=="boxv")   kernel=box_kernel((int)size,(int)size,(int)size);
	 else if(string(argv[i+1])=="boxv3")  { kernel=box_kernel((int)atof(argv[i+2]),(int)atof(argv[i+3]),(int)atof(argv[i+4])); i+=2; }
         else if(string(argv[i+1])=="gauss")  kernel=gaussian_kernel3D(size,xdim,ydim,zdim);
         else if(string(argv[i+1])=="sphere") kernel=spherical_kernel(size,xdim,ydim,zdim);
	 else if(string(argv[i+1])=="file")   read_volume(kernel,string(argv[i+2]));
	 i++;
       }
       i++;
       //save_volume(kernel,"kernel");
       
    }
    /***************************************************************/
    else if (string(argv[i])=="-s") 
    {
      kernel.destroy();
      float xdim=inputVolume.xdim();
      float ydim=inputVolume.ydim();
      float zdim=inputVolume.zdim();
      kernel=gaussian_kernel3D(atof(argv[i+1]),xdim,ydim,zdim);
      separable=true;
      inputVolume=generic_convolve(inputVolume,kernel,separable,true); 
      i++;
    }
    /***************************************************************/
    else if (string(argv[i])=="-subsamp2"){
      temp_volume.clear();
      extrapolation oldex = inputVolume.getextrapolationmethod();
      inputVolume.setextrapolationmethod(extraslice);
      temp_volume = subsample_by_2(inputVolume,true);
      temp_volume.setextrapolationmethod(oldex);
      inputVolume = temp_volume;
    }
    /***************************************************************/
    else if (string(argv[i])=="-subsamp2offc"){
      temp_volume.clear();
      extrapolation oldex = inputVolume.getextrapolationmethod();
      inputVolume.setextrapolationmethod(extraslice);
      temp_volume = subsample_by_2(inputVolume,false);
      temp_volume.setextrapolationmethod(oldex);
      inputVolume = temp_volume;
    }
    /***************************************************************/
    else if (string(argv[i])=="-sqrt")
    {
      for(int t=0;t<inputVolume.tsize();t++)           
        for(int z=0;z<inputVolume.zsize();z++)
          for(int y=0;y<inputVolume.ysize();y++)	    
	    for(int x=0;x<inputVolume.xsize();x++)
	    {
              if (inputVolume.value(x,y,z,t)> 0) inputVolume.value(x,y,z,t)=(T)sqrt(inputVolume.value(x,y,z,t));
              else  inputVolume.value(x,y,z,t)= 0; 
            }
    }
    /***************************************************************/
    else if (string(argv[i])=="-pow")
    {
      for(int t=0;t<inputVolume.tsize();t++)           
        for(int z=0;z<inputVolume.zsize();z++)
          for(int y=0;y<inputVolume.ysize();y++)	    
	    for(int x=0;x<inputVolume.xsize();x++)
	      inputVolume.value(x,y,z,t)=(T)pow(inputVolume.value(x,y,z,t),atof(argv[i+1]));
      i++;
    }
    /***************************************************************/
    else if (string(argv[i])=="-invbin")
    {
      for(int t=0;t<inputVolume.tsize();t++)           
        for(int z=0;z<inputVolume.zsize();z++)
          for(int y=0;y<inputVolume.ysize();y++)	    
	    for(int x=0;x<inputVolume.xsize();x++)
	      inputVolume.value(x,y,z,t)=(T)!inputVolume.value(x,y,z,t);
    }
    /***************************************************************/
    else if (string(argv[i])=="-sqr") inputVolume*=inputVolume;
    /***************************************************************/
    else if (string(argv[i])=="-recip")
      for(int t=0;t<inputVolume.tsize();t++)      
	{
          for(int z=0;z<inputVolume.zsize();z++) for(int y=0;y<inputVolume.ysize();y++) for(int x=0;x<inputVolume.xsize();x++)
            if(inputVolume.value(x,y,z,t)!=0)
	      inputVolume.value(x,y,z,t) = (T)1.0 / inputVolume.value(x,y,z,t);
        }          
    /***************************************************************/
    else if (string(argv[i])=="-exp")
    {
      for(int t=0;t<inputVolume.tsize();t++)           
        for(int z=0;z<inputVolume.zsize();z++)
          for(int y=0;y<inputVolume.ysize();y++)	    
	    for(int x=0;x<inputVolume.xsize();x++)
               inputVolume.value(x,y,z,t)=(T)exp((double)inputVolume.value(x,y,z,t));
    }
    /***************************************************************/
    else if (string(argv[i])=="-log")
    {
      for(int t=0;t<inputVolume.tsize();t++)           
        for(int z=0;z<inputVolume.zsize();z++)
          for(int y=0;y<inputVolume.ysize();y++)	    
	    for(int x=0;x<inputVolume.xsize();x++)
              if (inputVolume.value(x,y,z,t)> 0) inputVolume.value(x,y,z,t)=(T)log((double)inputVolume.value(x,y,z,t));
    }
    /***************************************************************/
    else if (string(argv[i])=="-cos")
    {
      for(int t=0;t<inputVolume.tsize();t++)           
        for(int z=0;z<inputVolume.zsize();z++)
          for(int y=0;y<inputVolume.ysize();y++)	    
	    for(int x=0;x<inputVolume.xsize();x++)
              inputVolume.value(x,y,z,t)=(T)cos((double)inputVolume.value(x,y,z,t));
    }
    /***************************************************************/
    else if (string(argv[i])=="-sin")
    {
      for(int t=0;t<inputVolume.tsize();t++)           
        for(int z=0;z<inputVolume.zsize();z++)
          for(int y=0;y<inputVolume.ysize();y++)	    
	    for(int x=0;x<inputVolume.xsize();x++)
              inputVolume.value(x,y,z,t)=(T)sin((double)inputVolume.value(x,y,z,t));
    }
    /***************************************************************/
    else if (string(argv[i])=="-tan")
    {
      for(int t=0;t<inputVolume.tsize();t++)           
        for(int z=0;z<inputVolume.zsize();z++)
          for(int y=0;y<inputVolume.ysize();y++)	    
	    for(int x=0;x<inputVolume.xsize();x++)
              inputVolume.value(x,y,z,t)=(T)tan((double)inputVolume.value(x,y,z,t));
    }
    /***************************************************************/
    else if (string(argv[i])=="-asin")
    {
      for(int t=0;t<inputVolume.tsize();t++)           
        for(int z=0;z<inputVolume.zsize();z++)
          for(int y=0;y<inputVolume.ysize();y++)	    
	    for(int x=0;x<inputVolume.xsize();x++)
              inputVolume.value(x,y,z,t)=(T)asin((double)inputVolume.value(x,y,z,t));
    }
    /***************************************************************/
    else if (string(argv[i])=="-acos")
    {
      for(int t=0;t<inputVolume.tsize();t++)           
        for(int z=0;z<inputVolume.zsize();z++)
          for(int y=0;y<inputVolume.ysize();y++)	    
	    for(int x=0;x<inputVolume.xsize();x++)
              inputVolume.value(x,y,z,t)=(T)acos((double)inputVolume.value(x,y,z,t));
    }
    /***************************************************************/
    else if (string(argv[i])=="-atan")
    {
      for(int t=0;t<inputVolume.tsize();t++)           
        for(int z=0;z<inputVolume.zsize();z++)
          for(int y=0;y<inputVolume.ysize();y++)	    
	    for(int x=0;x<inputVolume.xsize();x++)
              inputVolume.value(x,y,z,t)=(T)atan((double)inputVolume.value(x,y,z,t));
    }
    /***************************************************************/
    /* Uncorrected nonparametric P-value, assuming t-dim is perm-dim */
    else if (string(argv[i])=="-pval" || string(argv[i])=="-pval0")
    {
      int IgnoreZeros = string(argv[i])=="-pval0";

      // Reduce to 3D from 4D
      volume4D<T> temp_volume;
      temp_volume=inputVolume;
      inputVolume.reinitialize(inputVolume.xsize(),inputVolume.ysize(),inputVolume.zsize(),1); 
      inputVolume.copyproperties(temp_volume);
      
      // Compute p-value for each voxel
      float pval; int cnt;
      for(int z=0;z<inputVolume.zsize();z++)
        for(int y=0;y<inputVolume.ysize();y++)     
          for(int x=0;x<inputVolume.xsize();x++) {
            if (IgnoreZeros && temp_volume.value(x,y,z,0)==0)
              inputVolume.value(x,y,z,0)=(T)0;
            else {
              pval=1;
              if (IgnoreZeros) {
                cnt=1;
                for(int t=1;t<temp_volume.tsize();t++)
                  if (temp_volume.value(x,y,z,t)!=0) {
                    pval+=(temp_volume.value(x,y,z,0)<=temp_volume.value(x,y,z,t));
                    cnt++;
                  }
		pval/=cnt;
              } else {
                for(int t=1;t<temp_volume.tsize();t++)
                  pval+=(temp_volume.value(x,y,z,0)<=temp_volume.value(x,y,z,t));
                pval/=temp_volume.tsize();
              }
              inputVolume.value(x,y,z,0)=(T)pval;
            }
          }
    }
    /********* Corrected Nonparametric FWE P-value, assuming t-dim is perm-dim  ***************/
    else if (string(argv[i])=="-cpval")
    {
      // Initialize max distribution to min
      volume<float> max_dist(inputVolume.tsize(),1,1);
      for(int t=0;t<max_dist.xsize();t++)
        max_dist(t,0,0)=inputVolume.min();

      // Compute max distribution (t-dim is perm-dim)
      for(int t=0;t<inputVolume.tsize();t++)
        for(int z=0;z<inputVolume.zsize();z++)
          for(int y=0;y<inputVolume.ysize();y++)           
            for(int x=0;x<inputVolume.xsize();x++)
              max_dist(t,0,0)=MAX(max_dist(t,0,0),inputVolume.value(x,y,z,t));

      // Reduce to 3D from 4D
      volume4D<T> temp_volume;
      temp_volume=inputVolume;
      inputVolume.reinitialize(inputVolume.xsize(),inputVolume.ysize(),inputVolume.zsize(),1); 
      inputVolume.copyproperties(temp_volume);

      // Compute corrected p-value for each voxel
      float cpval;
      for(int z=0;z<inputVolume.zsize();z++)
        for(int y=0;y<inputVolume.ysize();y++)     
          for(int x=0;x<inputVolume.xsize();x++) {
            cpval=0;
            for(int t=0;t<max_dist.xsize();t++)
              cpval+=(temp_volume.value(x,y,z,0)<=max_dist(t,0,0));
            cpval/=temp_volume.tsize();
            inputVolume.value(x,y,z,0)=(T)cpval;
          }
    }
    /***** Uncorrected Parametric P-value for Normal Distribution ********/
    else if (string(argv[i])=="-ztop")
    {
      for(int t=0;t<inputVolume.tsize();t++)           
        for(int z=0;z<inputVolume.zsize();z++)
          for(int y=0;y<inputVolume.ysize();y++)           
            for(int x=0;x<inputVolume.xsize();x++)
              inputVolume.value(x,y,z,t)=(T)(1-ndtr((double)inputVolume.value(x,y,z,t)));
    }
    /***** Inverse Uncorrected Parametric P-value for Normal Distribution ********/
    else if (string(argv[i])=="-ptoz")
    {
      double v;
      for(int t=0;t<inputVolume.tsize();t++)           
        for(int z=0;z<inputVolume.zsize();z++)
          for(int y=0;y<inputVolume.ysize();y++)           
            for(int x=0;x<inputVolume.xsize();x++) {
              v=inputVolume.value(x,y,z,t);
              inputVolume.value(x,y,z,t)=(T)((v<=0 || v>=1)?0:ndtri(1-v));
            }
    }
    /***** Convert to ranks ********/
    else if (string(argv[i])=="-rank")
    {
      ColumnVector val(inputVolume.tsize()),rank(inputVolume.tsize()),sortval(inputVolume.tsize());
      for(int z=0;z<inputVolume.zsize();z++)
	for(int y=0;y<inputVolume.ysize();y++)           
	  for(int x=0;x<inputVolume.xsize();x++) {
	    for(int t=0;t<inputVolume.tsize();t++)
	      val(t+1)=inputVolume.value(x,y,z,t);
	    
	    
	    /* Take sortval and 'unsort' it, finding ranks */
	    sortval=val;
	    SortAscending(sortval);
	    for(int k=1;k<=val.Nrows();k++)	       
	      rank(k)=k;
	    for(int k=1;k<=val.Nrows();k++)	       
	      if(val(k)!=sortval(k))
		for(int l=k+1;l<=val.Nrows();l++) 
		  if(sortval(l)==val(k))
		    {
		      swap(rank(l),rank(k));
		      swap(sortval(l),sortval(k));
		    }
	    for(int t=0;t<inputVolume.tsize();t++)
	      inputVolume.value(x,y,z,t)=(T)rank(t+1);
	  }

    }
    /***** Convert to normal distribution via ranks ********/
    else if (string(argv[i])=="-ranknorm")
    {
      ColumnVector val(inputVolume.tsize()),rank(inputVolume.tsize()),sortval(inputVolume.tsize());
      volume<double> valv(inputVolume.tsize(),1,1);
      for(int z=0;z<inputVolume.zsize();z++)
	for(int y=0;y<inputVolume.ysize();y++)           
	  for(int x=0;x<inputVolume.xsize();x++) {
	    for(int t=0;t<inputVolume.tsize();t++) {
	      val(t+1)=inputVolume.value(x,y,z,t); 
	      valv(t,0,0)=val(t+1);
	    }
	    
	    /* Take sortval and 'unsort' it, finding ranks */
	    sortval=val;
	    SortAscending(sortval);
	    for(int k=1;k<=val.Nrows();k++)	       
	      rank(k)=k;
	    for(int k=1;k<=val.Nrows();k++)	       
	      if(val(k)!=sortval(k))
		for(int l=k+1;l<=val.Nrows();l++) 
		  if(sortval(l)==val(k))
		    {
		      swap(rank(l),rank(k));
		      swap(sortval(l),sortval(k));
		    }
	    /* Transform to expected order statisics of a Uniform */
	    rank = (rank-0.5)/rank.Nrows();
	    /* Transform to expected order statisics of a Normal with same Mean & Sd */
	    for(int t=0;t<inputVolume.tsize();t++)
	      inputVolume.value(x,y,z,t)=(T)(valv.mean()+ndtri(rank(t+1))*valv.stddev()); 
	    
	  }

    }
    /***************************************************************/
    else if (string(argv[i])=="-abs") inputVolume=abs(inputVolume);
    /***************************************************************/
    else if (string(argv[i])=="-bin") inputVolume.binarise(0,inputVolume.max()+1,exclusive); 
    /***************************************************************/
    else if (string(argv[i])=="-binv") 
      {
	inputVolume.binarise(0,inputVolume.max()+1,exclusive);
	inputVolume*=(T)(-1);
	inputVolume+=(T)1;
	inputVolume.binarise((T)0.5,inputVolume.max()+1,exclusive);  // do this to cope with a bit of rounding error
      }
    /***************************************************************/
    else if (string(argv[i])=="-fillh") 
      {
	inputVolume.binarise(0,inputVolume.max()+1,exclusive);
	for (int t=0; t<=inputVolume.maxt(); t++) {
	  inputVolume[t] = fill_holes(inputVolume[t],6);
	}
	inputVolume.binarise(0,inputVolume.max()+1,exclusive);
      }
    /***************************************************************/
    else if (string(argv[i])=="-fillh26") 
      {
	inputVolume.binarise(0,inputVolume.max()+1,exclusive);
	for (int t=0; t<=inputVolume.maxt(); t++) {
	  inputVolume[t] = fill_holes(inputVolume[t],26);
	}
	inputVolume.binarise(0,inputVolume.max()+1,exclusive);
      }
    /***************************************************************/
    else if (string(argv[i])=="-index")
    {
      int indexval=0;
      for(int t=0;t<inputVolume.tsize();t++)           
        for(int z=0;z<inputVolume.zsize();z++)
          for(int y=0;y<inputVolume.ysize();y++)	    
	    for(int x=0;x<inputVolume.xsize();x++)
	      if (inputVolume.value(x,y,z,t)>0) 
		{
		  inputVolume.value(x,y,z,t)=(T)indexval;
		  indexval++;
		}
    }
    /***************************************************************/
    else if (string(argv[i])=="-grid")
    {
      double gridvalue = atof(argv[++i]);
      int gridspacing = atoi(argv[++i]);
      for(int t=0;t<inputVolume.tsize();t++)           
        for(int z=0;z<inputVolume.zsize();z++)
          for(int y=0;y<inputVolume.ysize();y++)	    
	    for(int x=0;x<inputVolume.xsize();x++)
	      if ( x%gridspacing==0 || y%gridspacing==0 || z%gridspacing==0 )
		inputVolume.value(x,y,z,t)=(T)gridvalue;
    }
    /*****************SPATIAL FILTERING OPTIONS*********************/
    /***********************All Dilation***************************/
    else if (string(argv[i])=="-dilall") {
        volume4D<T> mask(inputVolume);   
        mask.binarise(0,0,inclusive); // get a zero mask
	mask=(T)1-mask;  // convert to a non-zero mask
	for(int t=0;t<inputVolume.tsize();t++) dilall(inputVolume[t],mask[t]);
    }
    /***********************Mean Dilation***************************/
    else if (string(argv[i])=="-dilM")
	for(int t=0;t<inputVolume.tsize();t++) inputVolume[t]=morphfilter(inputVolume[t],kernel,"dilateM");
    /***********************Modal Dilation**************************/
    else if (string(argv[i])=="-dilD")
	for(int t=0;t<inputVolume.tsize();t++) inputVolume[t]=morphfilter(inputVolume[t],kernel,"dilateD");
    /***********************MJ Dilation**************************/
    else if (string(argv[i])=="-dilF")
	for(int t=0;t<inputVolume.tsize();t++) inputVolume[t]=morphfilter(inputVolume[t],kernel,"dilate");
    /***********************Steves Erosion**************************/
    else if (string(argv[i])=="-ero")
	for(int t=0;t<inputVolume.tsize();t++) inputVolume[t]=morphfilter(inputVolume[t],kernel,"erodeS");
    /**************************MJ Erosion**************************/
    else if (string(argv[i])=="-eroF")
	for(int t=0;t<inputVolume.tsize();t++) inputVolume[t]=morphfilter(inputVolume[t],kernel,"erode");
    /***********************Median Filtering***********************/
    else if (string(argv[i])=="-fmedian")
	for(int t=0;t<inputVolume.tsize();t++) inputVolume[t]=morphfilter(inputVolume[t],kernel,"median");
    /******************Mean Filtering*************************/
    else if (string(argv[i])=="-fmean")
       inputVolume=generic_convolve(inputVolume,kernel,separable,true);
    /******************Mean Filtering Unnormalised************/
    else if (string(argv[i])=="-fmeanu")
       inputVolume=generic_convolve(inputVolume,kernel,false,false);
    /*****************END OF FILTERING OPTIONS***************/
    else if (string(argv[i])=="-edge")
       inputVolume=edge_strengthen(inputVolume);
    else if (string(argv[i])=="-tfce_minT")
      tfce_minT = atof(argv[++i]);
    else if (string(argv[i])=="-tfce_delta")
      tfce_delta = atof(argv[++i]);
    else if (string(argv[i])=="-tfce") {
      float height_power = atof(argv[++i]);
      float size_power = atof(argv[++i]);
      int connectivity = atoi(argv[++i]);

      for(int t=0;t<inputVolume.tsize();t++) {
	try { 
	  tfce(inputVolume[t], height_power, size_power, connectivity, tfce_minT, tfce_delta);
	}
	catch(Exception& e) { 
	  cerr << "ERROR: TFCE failed, please check your file for valid sizes and voxel values." <<  e.what() << endl << endl << "Exiting" << endl; 
	  return 1;
	}
      }
    }
// }}}
    else if (string(argv[i])=="-tfceS")
      // {{{ TFCE supporting-area option

      {
	float height_power = atof(argv[++i]);
	float size_power = atof(argv[++i]);
	int connectivity = atoi(argv[++i]);
	int X = atoi(argv[++i]);
	int Y = atoi(argv[++i]);
	int Z = atoi(argv[++i]);
	float tfce_thresh = atof(argv[++i]);
	
	for(int t=0;t<inputVolume.tsize();t++)
	  tfce_support(inputVolume[t], height_power, size_power, connectivity, tfce_minT, tfce_delta, X, Y, Z, tfce_thresh);
      }

// }}}
    /******************************************************/
    else if (string(argv[i])=="-nanm")
     {   
       for(int t=0;t<inputVolume.tsize();t++)           
         for(int z=0;z<inputVolume.zsize();z++)
           for(int y=0;y<inputVolume.ysize();y++)	    
	     for(int x=0;x<inputVolume.xsize();x++)
               if ( isfinite((double)inputVolume.value(x,y,z,t))) inputVolume.value(x,y,z,t)=0;
	       else inputVolume.value(x,y,z,t)=1;
     }
     /******************************************************/
    else if (string(argv[i])=="-nan")
     {   
       for(int t=0;t<inputVolume.tsize();t++)           
         for(int z=0;z<inputVolume.zsize();z++)
           for(int y=0;y<inputVolume.ysize();y++)	    
	     for(int x=0;x<inputVolume.xsize();x++)
               if (!isfinite((double)inputVolume.value(x,y,z,t))) inputVolume.value(x,y,z,t)=0;     
     }
     /******************************************************/
    else if (string(argv[i])=="-rand")
     {   
       double rnd;
       for(int t=0;t<inputVolume.tsize();t++)           
         for(int z=0;z<inputVolume.zsize();z++)
           for(int y=0;y<inputVolume.ysize();y++)          
             for(int x=0;x<inputVolume.xsize();x++) {
               drand(&rnd);
               inputVolume.value(x,y,z,t) += (T)(rnd-1);
             }
     }
     /******************************************************/
    else if (string(argv[i])=="-randn")
     {   
       double rnd;
       for(int t=0;t<inputVolume.tsize();t++)           
         for(int z=0;z<inputVolume.zsize();z++)
           for(int y=0;y<inputVolume.ysize();y++)          
             for(int x=0;x<inputVolume.xsize();x++) {
               drand(&rnd);
               inputVolume.value(x,y,z,t) += (T)ndtri(rnd-1);
             }
     }
     /******************************************************/
    else if (string(argv[i])=="-roi")
     { 
       int x0=atoi(argv[i+1]);
       int x1=atoi(argv[i+1])+atoi(argv[i+2])-1;
       int y0=atoi(argv[i+3]);
       int y1=atoi(argv[i+3])+atoi(argv[i+4])-1;
       int z0=atoi(argv[i+5]);
       int z1=atoi(argv[i+5])+atoi(argv[i+6])-1;
       int t0=atoi(argv[i+7]);
       int t1=atoi(argv[i+7])+atoi(argv[i+8])-1;
       if ( x1 == x0-2 )  //Then the user input -1 for size
	 x1 = inputVolume.maxx();
       if ( y1 == y0-2 )  //Then the user input -1 for size
	 y1 = inputVolume.maxy();
       if ( z1 == z0-2 )  //Then the user input -1 for size
	 z1 = inputVolume.maxz();
       if ( t1 == t0-2 )  //Then the user input -1 for size
	 t1 = inputVolume.maxt();
       ColumnVector v0(4), v1(4);
       v0 << x0 << y0 << z0 << 1.0;
       v1 << x1 << y1 << z1 << 1.0;
       v0 = inputVolume.niftivox2newimagevox_mat() * v0;
       v1 = inputVolume.niftivox2newimagevox_mat() * v1;
       x0=MISCMATHS::round(v0(1));
       y0=MISCMATHS::round(v0(2));
       z0=MISCMATHS::round(v0(3));
       x1=MISCMATHS::round(v1(1));
       y1=MISCMATHS::round(v1(2));
       z1=MISCMATHS::round(v1(3));
       // swap back to restore min/max order as necessary
       int tmp;
       if (x0>x1) { tmp=x0; x0=x1;  x1=tmp; }
       if (y0>y1) { tmp=y0; y0=y1;  y1=tmp; }
       if (z0>z1) { tmp=z0; z0=z1;  z1=tmp; }
       for(int t=0;t<inputVolume.tsize();t++)           
         for(int z=0;z<inputVolume.zsize();z++)
           for(int y=0;y<inputVolume.ysize();y++)	    
	     for(int x=0;x<inputVolume.xsize();x++)
               if((x<x0) || (x>x1) || (y<y0) || (y>y1) || (z<z0) || (z>z1) || (t<t0) || (t>t1) )
                 inputVolume.value(x,y,z,t)=0;
       i+=8;
     }
     /*******************IP functions***********************/
    else if (string(argv[i])=="-inm")
     { 
       double target,tmean;
       target = atof(argv[++i]);
       volume4D<T> mask(inputVolume);   
       mask.binarise(0,mask.max()+1,exclusive); 
       for(int t=0;t<inputVolume.tsize();t++)     
       {
         tmean=target/inputVolume[t].mean(mask[t]);
         for(int z=0;z<inputVolume.zsize();z++)
           for(int y=0;y<inputVolume.ysize();y++)	    
	     for(int x=0;x<inputVolume.xsize();x++)
               inputVolume.value(x,y,z,t)=(T)(inputVolume.value(x,y,z,t)*tmean);
       }
     }
    else if (string(argv[i])=="-ing")
     { 
       double tmean,target;
       target = atof(argv[++i]);
       volume4D<T> mask(inputVolume);   
       mask.binarise(0,mask.max()+1,exclusive); 
       tmean=target/inputVolume.mean(mask);
       for(int t=0;t<inputVolume.tsize();t++)     
         for(int z=0;z<inputVolume.zsize();z++)
           for(int y=0;y<inputVolume.ysize();y++)	    
	     for(int x=0;x<inputVolume.xsize();x++)
               inputVolume.value(x,y,z,t)=(T)(inputVolume.value(x,y,z,t)*tmean);
     }
    else if (string(argv[i])=="-tensor_decomp")
     { 
       // {{{ tensor decomposition

       if ( inputVolume.tsize() != 6 )
	 {
	   cout << "Error: input to tensor option is not a tensor!" << endl;
	   exit(1);
	 }
       volume<float> dti_L1(inputVolume.xsize(),inputVolume.ysize(),inputVolume.zsize()),
	 dti_L2(inputVolume.xsize(),inputVolume.ysize(),inputVolume.zsize()),
	 dti_L3(inputVolume.xsize(),inputVolume.ysize(),inputVolume.zsize()),
	 dti_FA(inputVolume.xsize(),inputVolume.ysize(),inputVolume.zsize()),
	 dti_MD(inputVolume.xsize(),inputVolume.ysize(),inputVolume.zsize()),
	 dti_MO(inputVolume.xsize(),inputVolume.ysize(),inputVolume.zsize());
       volume4D<float> dti_V1(inputVolume.xsize(),inputVolume.ysize(),inputVolume.zsize(),3),
	 dti_V2(inputVolume.xsize(),inputVolume.ysize(),inputVolume.zsize(),3),
	 dti_V3(inputVolume.xsize(),inputVolume.ysize(),inputVolume.zsize(),3);

       dti_L1=0; dti_L2=0; dti_L3=0; dti_FA=0; dti_MD=0; dti_MO=0; dti_V1=0; dti_V2=0; dti_V3=0;
       
       copybasicproperties(inputVolume,dti_L1);
       copybasicproperties(inputVolume,dti_L2);
       copybasicproperties(inputVolume,dti_L3);
       copybasicproperties(inputVolume,dti_FA);
       copybasicproperties(inputVolume,dti_MD);
       copybasicproperties(inputVolume,dti_MO);
       copybasicproperties(inputVolume,dti_V1);
       copybasicproperties(inputVolume,dti_V2);
       copybasicproperties(inputVolume,dti_V3);

       for(int z=0;z<inputVolume.zsize();z++)
	 for(int y=0;y<inputVolume.ysize();y++)	    
	   for(int x=0;x<inputVolume.xsize();x++)
	     {
	       SymmetricMatrix dti_tensor(3);
	       dti_tensor(1,1) = (float)inputVolume.value(x,y,z,0);
	       dti_tensor(2,1) = (float)inputVolume.value(x,y,z,1);
	       dti_tensor(2,2) = (float)inputVolume.value(x,y,z,3); // note strange order because we're filling lower part of matrix
	       dti_tensor(3,1) = (float)inputVolume.value(x,y,z,2); // and the existing tensor coeffients assumed upper
	       dti_tensor(3,2) = (float)inputVolume.value(x,y,z,4);
	       dti_tensor(3,3) = (float)inputVolume.value(x,y,z,5);

	       Matrix dti_V(3,3);
	       DiagonalMatrix dti_L(3);

	       EigenValues(dti_tensor,dti_L,dti_V);

	       if ( dti_L(3) > 0 )
		 {
		   dti_L1.value(x,y,z)=dti_L(3);
		   dti_L2.value(x,y,z)=dti_L(2);
		   dti_L3.value(x,y,z)=dti_L(1);

		   float lmean=(dti_L1.value(x,y,z)+dti_L2.value(x,y,z)+dti_L3.value(x,y,z))/3;
		   float lnum=3*( (dti_L1.value(x,y,z)-lmean)*(dti_L1.value(x,y,z)-lmean) +
				  (dti_L2.value(x,y,z)-lmean)*(dti_L2.value(x,y,z)-lmean) +
				  (dti_L3.value(x,y,z)-lmean)*(dti_L3.value(x,y,z)-lmean) );
		   float lden=2*( dti_L1.value(x,y,z)*dti_L1.value(x,y,z) +
				  dti_L2.value(x,y,z)*dti_L2.value(x,y,z) +
				  dti_L3.value(x,y,z)*dti_L3.value(x,y,z) );
		   dti_FA.value(x,y,z)=sqrt(lnum/lden);

		   dti_MD.value(x,y,z)=lmean;

		   float e1=dti_L1.value(x,y,z)-lmean, e2=dti_L2.value(x,y,z)-lmean, e3=dti_L3.value(x,y,z)-lmean;
		   float n = (e1 + e2 - 2*e3)*(2*e1 - e2 - e3)*(e1 - 2*e2 + e3);
		   float d = (e1*e1 + e2*e2 + e3*e3 - e1*e2 - e2*e3 - e1*e3);
		   d = sqrt(MAX(0, d));
		   d = 2*d*d*d;
		   dti_MO.value(x,y,z) = MIN(MAX(d ? n/d : 0.0, -1),1);

		   dti_V1.value(x,y,z,0)=dti_V(1,3);
		   dti_V1.value(x,y,z,1)=dti_V(2,3);
		   dti_V1.value(x,y,z,2)=dti_V(3,3);
		   dti_V2.value(x,y,z,0)=dti_V(1,2);
		   dti_V2.value(x,y,z,1)=dti_V(2,2);
		   dti_V2.value(x,y,z,2)=dti_V(3,2);
		   dti_V3.value(x,y,z,0)=dti_V(1,1);
		   dti_V3.value(x,y,z,1)=dti_V(2,1);
		   dti_V3.value(x,y,z,2)=dti_V(3,1);
		 }
	       
	       inputVolume.value(x,y,z,0) = (T)dti_FA.value(x,y,z);
	     }
       for(int j=5;j>0;j--)
	 inputVolume.deletevolume(j);

       // if i+1>argc-1 then can save (otherwise it is a syntax error with no specific output specified)
       check_for_output_name(i+1,argc-1);
       dti_L1.setDisplayMaximumMinimum(dti_L1.max(),dti_L1.min());
       save_volume(dti_L1,string(argv[argc-1])+"_L1");
       dti_L2.setDisplayMaximumMinimum(dti_L2.max(),dti_L2.min());
       save_volume(dti_L2,string(argv[argc-1])+"_L2");
       dti_L3.setDisplayMaximumMinimum(dti_L3.max(),dti_L3.min());
       save_volume(dti_L3,string(argv[argc-1])+"_L3");
       dti_FA.setDisplayMaximumMinimum(1,0);
       save_volume(dti_FA,string(argv[argc-1])+"_FA");
       dti_MD.setDisplayMaximumMinimum(dti_MD.max(),dti_MD.min());
       save_volume(dti_MD,string(argv[argc-1])+"_MD");
       dti_MO.setDisplayMaximumMinimum(1,-1);
       save_volume(dti_MO,string(argv[argc-1])+"_MO");
       dti_V1.setDisplayMaximumMinimum(1,-1);
       save_volume4D(dti_V1,string(argv[argc-1])+"_V1");
       dti_V2.setDisplayMaximumMinimum(1,-1);
       save_volume4D(dti_V2,string(argv[argc-1])+"_V2");
       dti_V3.setDisplayMaximumMinimum(1,-1);
       save_volume4D(dti_V3,string(argv[argc-1])+"_V3");

// }}}
     }
    else if(string(argv[i])=="-bptf")
     {
       inputVolume=bandpass_temporal_filter(inputVolume,atof(argv[i+1]),atof(argv[i+2]));
       i+=2;
     }
    else { cout << "\n Error in command line: unknown option \"" << argv[i] << "\"\n" << endl; return printUsage("blah"); }
     /******************************************************/
  } 

  // if i>argc-1 then can save (otherwise it is a syntax error with no specific output specified)
  check_for_output_name(i,argc-1);

  if (dtype(inputVolume)>=DT_FLOAT && output_dt < DT_FLOAT)
  {
    for(int t=0;t<inputVolume.tsize();t++)           
        for(int z=0;z<inputVolume.zsize();z++)
          for(int y=0;y<inputVolume.ysize();y++)	    
	    for(int x=0;x<inputVolume.xsize();x++)
              inputVolume.value(x,y,z,t)=(T) MISCMATHS::round(inputVolume.value(x,y,z,t));
  }

  if (modifiedInput)
    inputVolume.setDisplayMaximumMinimum(0,0);
  if (setDisplayRange)
    inputVolume.setDisplayMaximumMinimum((float)inputVolume.max(),(float)inputVolume.min());


  save_volume4D_datatype(inputVolume,string(argv[argc-1]),output_dt);
  return 0;
}


int main(int argc,char *argv[])
{
  if (argc < 2)  
    return printUsage(string(argv[0])); 
  
  short inputType;
  if(string(argv[1]) == "-h" || string(argv[1]) == "--help") { printUsage(string(argv[0])); exit(0); }
  if(string(argv[1]) =="-datatype" || string(argv[1])== "-dt")  inputType = dtype(string(argv[3]));
  else inputType = dtype(string(argv[1]));
  short outputType = DT_FLOAT;
  if ( inputType == DT_DOUBLE )
    outputType = DT_DOUBLE; 

  if(string(argv[argc-2])=="-output_datatype" || string(argv[argc-2])== "-odt") //output datatype
  {
    if ( string(argv[argc-1]) == "input" )   outputType =  inputType;
    else if(string(argv[argc-1])=="char")   outputType =  DT_UNSIGNED_CHAR;  
    else if(string(argv[argc-1])=="short")  outputType =  DT_SIGNED_SHORT;
    else if(string(argv[argc-1])=="int")    outputType =  DT_SIGNED_INT;
    else if(string(argv[argc-1])=="float")  outputType =  DT_FLOAT;
    else if(string(argv[argc-1])=="double") outputType =  DT_DOUBLE;
    else {cout << "Error: Unknown datatype \"" << argv[argc-1] << "\" - Possible datatypes are: char short int float double input" << endl; return 1;}
    argc-=2;
  }
  if(string(argv[1])=="-datatype" || string(argv[1])== "-dt") //input datatype
  {
    if(string(argv[2])!="input") inputType=-1;
    if(string(argv[2])=="char" || inputType == DT_UNSIGNED_CHAR)      return inputParser<char>(argc-2, argv+2,outputType);
    else if(string(argv[2])=="short" || inputType == DT_SIGNED_SHORT) return inputParser<short>(argc-2, argv+2,outputType);
    else if(string(argv[2])=="int"   || inputType == DT_SIGNED_INT)   return inputParser<int>(argc-2, argv+2,outputType);
    else if(string(argv[2])=="float" || inputType == DT_FLOAT)        return inputParser<float>(argc-2, argv+2,outputType); 
    else if(string(argv[2])=="double"|| inputType == DT_DOUBLE)       return inputParser<double>(argc-2, argv+2,outputType); 
    else if (inputType==-1)
      { cout << "Error: Unknown datatype \"" << argv[2] <<  "\" - Possible datatypes are: char short int float double input" << endl; return 1;}
  }
  else if (dtype(string(argv[1]))==DT_DOUBLE) return inputParser<double>(argc,argv,outputType);
  else return inputParser<float>(argc,argv,outputType);
}



