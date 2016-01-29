/*  fslstats.cc

    Mark Jenkinson and Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 2003-2009 University of Oxford  */

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

#include "miscmaths/miscmaths.h"
#include "newimage/newimageall.h"
#include "newimage/costfns.h"
#include "utils/fsl_isfinite.h"

using namespace NEWIMAGE;

int print_usage(const string& progname) {
  cout << "Usage: fslstats [preoptions] <input> [options]" << endl << endl; 
  cout << "preoption -t will give a separate output line for each 3D volume of a 4D timeseries" << endl;
  cout << "preoption -K < indexMask > will generate seperate n submasks from indexMask, for indexvalues 1..n where n is the maximum index value in indexMask, and generate statistics for each submask" << endl;
  cout << "Note - options are applied in order, e.g. -M -l 10 -M will report the non-zero mean, apply a threshold and then report the new nonzero mean" << endl << endl;
  cout << "-l <lthresh> : set lower threshold" << endl;
  cout << "-u <uthresh> : set upper threshold" << endl;
  cout << "-r           : output <robust min intensity> <robust max intensity>" << endl;
  cout << "-R           : output <min intensity> <max intensity>" << endl;
  cout << "-e           : output mean entropy ; mean(-i*ln(i))" << endl;
  cout << "-E           : output mean entropy (of nonzero voxels)" << endl;
  cout << "-v           : output <voxels> <volume>" << endl;
  cout << "-V           : output <voxels> <volume> (for nonzero voxels)" << endl;
  cout << "-m           : output mean" << endl;
  cout << "-M           : output mean (for nonzero voxels)" << endl;
  cout << "-s           : output standard deviation" << endl;
  cout << "-S           : output standard deviation (for nonzero voxels)" << endl;
  cout << "-w           : output smallest ROI <xmin> <xsize> <ymin> <ysize> <zmin> <zsize> <tmin> <tsize> containing nonzero voxels" << endl;
  cout << "-x           : output co-ordinates of maximum voxel" << endl;
  cout << "-X           : output co-ordinates of minimum voxel" << endl;
  cout << "-c           : output centre-of-gravity (cog) in mm coordinates" << endl;
  cout << "-C           : output centre-of-gravity (cog) in voxel coordinates" << endl;
  cout << "-p <n>       : output nth percentile (n between 0 and 100)" << endl;
  cout << "-P <n>       : output nth percentile (for nonzero voxels)" << endl;
  cout << "-a           : use absolute values of all image intensities"<< endl;
  cout << "-n           : treat NaN or Inf as zero for subsequent stats" << endl;
  cout << "-k <mask>    : use the specified image (filename) for masking - overrides lower and upper thresholds" << endl;
  cout << "-d <image>   : take the difference between the base image and the image specified here" << endl;
  cout << "-h <nbins>   : output a histogram (for the thresholded/masked voxels only) with nbins" << endl; 
  cout << "-H <nbins> <min> <max>   : output a histogram (for the thresholded/masked voxels only) with nbins and histogram limits of min and max" << endl << endl;
  cout << "Note - thresholds are not inclusive ie lthresh<allowed<uthresh" << endl;
  return 1;
}

// Some specialised nonzero functions just for speedup
//  (it avoids generating masks when not absolutely necessary)

long int nonzerocount(const volume4D<float>& vol)
{
  long int totn=0;
  for (int t=vol.mint(); t<=vol.maxt(); t++) {
    for (int z=vol.minz(); z<=vol.maxz(); z++) {
      for (int y=vol.miny(); y<=vol.maxy(); y++) {
	for (int x=vol.minx(); x<=vol.maxx(); x++) {
	  if (vol(x,y,z,t)!=0.0) {
	    totn++;
	  }
	}
      }
    }
  }
  return totn;
}

double nonzeromean(const volume4D<float>& vol)
{
  double totv=0.0;
  long int totn=0;
  for (int t=vol.mint(); t<=vol.maxt(); t++) {
    for (int z=vol.minz(); z<=vol.maxz(); z++) {
      for (int y=vol.miny(); y<=vol.maxy(); y++) {
	for (int x=vol.minx(); x<=vol.maxx(); x++) {
	  if (vol(x,y,z,t)!=0.0) {
	    totv+=(double) vol(x,y,z,t);
	    totn++;
	  }
	}
      }
    }
  }
  double meanval=0.0;
  if (totn>0) {
    meanval=totv/totn;
  }
  return meanval;
}

double nonzerostddev(const volume4D<float>& vol)
{
  double totv=0.0, totvv=0.0;
  long int totn=0;
  for (int t=vol.mint(); t<=vol.maxt(); t++) {
    for (int z=vol.minz(); z<=vol.maxz(); z++) {
      for (int y=vol.miny(); y<=vol.maxy(); y++) {
	for (int x=vol.minx(); x<=vol.maxx(); x++) {
	  if (vol(x,y,z,t)!=0.0) {
	    float v=vol(x,y,z,t);
	    totvv+=(double) v*v;
	    totv+=(double) v;
	    totn++;
	  }
	}
      }
    }
  }
  double var=0.0;
  if (totn>1) {
    double meanval = totv/totn;
    var = (totvv - totn*meanval*meanval)/(totn-1);
  }
  return std::sqrt(var);
}

int generateNonZeroMask(const volume4D<float> &mask, volume4D<float> &masknz, const volume4D<float> &input)
{
  masknz.reinitialize(mask.xsize(),mask.ysize(),mask.zsize(),input.tsize());
  for (int t=input.mint(); t<=input.maxt(); t++) 
    masknz[t]=((binarise(input[t],0.0f, 0.0f)-1.0f)*-1.0f*mask[t % mask.tsize()]);
  return 0;
}

int generate_masks(volume4D<float>& mask, volume4D<float>& masknz, const volume4D<float>& input, const float& lthr, const float& uthr) 
{
  mask = binarise(input,lthr,uthr,exclusive);
  return generateNonZeroMask(mask,masknz,input);
}

int fmrib_main_float(int argc, char* argv[],const bool timeseriesMode, const string& indexMaskName) 
{
  cout.setf(ios::dec); 
  cout.setf(ios::fixed, ios::floatfield); 
  cout.setf(ios::left, ios::adjustfield); 
  cout.precision(6);  
  volume4D<float> vol, inputMaster;
  volume4D<int> indexMask;
  if ( timeseriesMode || indexMaskName != "" ) 
    read_volume4D(inputMaster,argv[1]);
  else 
    read_volume4D(vol,argv[1]);
  volume4D<float> & indexMaster = (timeseriesMode ) ? vol : inputMaster;
  if ( indexMaskName != "" ) 
    read_volume4D(indexMask,indexMaskName);
  int nTimepoints((timeseriesMode) ? inputMaster.tsize() : 1), nIndices((indexMaskName != "") ? indexMask.max() : 1);
  for ( int timepoint=0; timepoint < nTimepoints ; timepoint++ ) {
   for ( int index=1; index <= nIndices; index++ ) {
    if ( timeseriesMode )
      vol=inputMaster[timepoint];
    volume4D<float> mask, masknz;
    if ( indexMask.nvoxels() ) {
      if ( indexMask.tsize() != 1 ) {
	cerr << "Index mask must be 3D" << endl;
        return 3;
      }
      copyconvert(indexMask,mask);
      mask.binarise(index-1,index+1,exclusive);
      vol=indexMaster*mask[0];
      generateNonZeroMask(mask,masknz,vol);
    }
    float lthr(vol.min()-1);
    float uthr=(vol.max()+1);    
    int narg(2);
 
  while (narg<argc) {
    string sarg(argv[narg]);
    if (sarg=="-n") {
      for (int t=vol.mint(); t<=vol.maxt(); t++)
        for (int z=vol.minz(); z<=vol.maxz(); z++)
          for (int y=vol.miny(); y<=vol.maxy(); y++)
            for (int x=vol.minx(); x<=vol.maxx(); x++)
              if (!isfinite((double)vol(x,y,z,t)))
	        vol(x,y,z,t)=0;
    } else if (sarg=="-m") {
      if (mask.nvoxels()>0) cout <<  vol.mean(mask) << " ";
      else cout << vol.mean() << " ";
    } else if (sarg=="-M") {
      if (masknz.nvoxels()>0) cout << vol.mean(masknz) << " ";
      else {
	double nzmean=0;
	nzmean = nonzeromean(vol);
	cout << nzmean << " ";
      }
    } else if (sarg=="-X") {
      ColumnVector coord(4);
      coord(4)=1.0;
      if (mask.nvoxels()>0) {
	coord(1) = vol.mincoordx(mask);
	coord(2) = vol.mincoordy(mask);
	coord(3) = vol.mincoordz(mask);
      } else {
	coord(1) = vol.mincoordx();
	coord(2) = vol.mincoordy();
	coord(3) = vol.mincoordz();
      }
      coord = (vol[0].niftivox2newimagevox_mat()).i() * coord;
      cout << MISCMATHS::round(coord(1)) << " " << 
	MISCMATHS::round(coord(2)) << " " << MISCMATHS::round(coord(3)) << " ";
    } else if (sarg=="-x") { 
      ColumnVector coord(4);
      coord(4)=1.0;
      if (mask.nvoxels()>0) {
	coord(1) = vol.maxcoordx(mask);
	coord(2) = vol.maxcoordy(mask);
	coord(3) = vol.maxcoordz(mask);
      } else {
	coord(1) = vol.maxcoordx();
	coord(2) = vol.maxcoordy();
	coord(3) = vol.maxcoordz();
      }
      coord = (vol[0].niftivox2newimagevox_mat()).i() * coord;
      cout << MISCMATHS::round(coord(1)) << " " << 
	MISCMATHS::round(coord(2)) << " " << MISCMATHS::round(coord(3)) << " ";
    } else if (sarg=="-w") {
	if (masknz.nvoxels()<1) { //Need to generate non-zeromask 
	  generate_masks(mask,masknz,vol,lthr,uthr); 
	  vol*=mask; 
	}
      int xmin=masknz.maxx(),xmax=masknz.minx(),ymin=masknz.maxy(),ymax=masknz.miny(),zmin=masknz.maxz(),zmax=masknz.minz(),tmin=masknz.maxt(),tmax=masknz.mint();
      
      for(int t=masknz.mint();t<=masknz.maxt();t++) {
	for(int z=masknz.minz();z<=masknz.maxz();z++) {
	  for(int y=masknz.miny();y<=masknz.maxy();y++) {
	    for(int x=masknz.minx();x<=masknz.maxx();x++) {
	      if (masknz(x,y,z,t)>0.5) {
		// if (masknz(x,y,z)>0.5) {
		if (x<xmin) xmin=x;
		if (x>xmax) xmax=x;
		if (y<ymin) ymin=y;
		if (y>ymax) ymax=y;
		if (z<zmin) zmin=z;
		if (z>zmax) zmax=z;
		if (t<tmin) tmin=t;
		if (t>tmax) tmax=t;
	      }
	    }
	  }
	}
      }
      // change voxel coords from newimage to nifti convention for output
      ColumnVector v(4);
      v << xmin << ymin << zmin << 1.0;
      v = masknz.niftivox2newimagevox_mat().i() * v;
      xmin = MISCMATHS::round(v(1));
      ymin = MISCMATHS::round(v(2));
      zmin = MISCMATHS::round(v(3));
      v << xmax << ymax << zmax << 1.0;
      v = masknz.niftivox2newimagevox_mat().i() * v;
      xmax = MISCMATHS::round(v(1));
      ymax = MISCMATHS::round(v(2));
      zmax = MISCMATHS::round(v(3));
      if (xmin>xmax) { int tmp=xmax;  xmax=xmin;  xmin=tmp; }
      if (ymin>ymax) { int tmp=ymax;  ymax=ymin;  ymin=tmp; }
      if (zmin>zmax) { int tmp=zmax;  zmax=zmin;  zmin=tmp; }
      // now output nifti coords
      cout << xmin << " " << 1+xmax-xmin << " " << ymin << " " << 1+ymax-ymin << " " << zmin << " " << 1+zmax-zmin << " " << tmin << " " << 1+tmax-tmin << " ";
      } else if (sarg=="-e") {
	if (mask.nvoxels()<1) {
	  generate_masks(mask,masknz,vol,lthr,uthr); 
	  vol*=mask; 
	}
      ColumnVector hist;
      int nbins=1000;
      double entropy=0;
      hist = vol.histogram(nbins,mask);
      double ntot = hist.Sum();
      for (int j=1; j<=nbins; j++) {
	if (hist(j)>0) {
	  entropy -= (hist(j)/ntot) * log(hist(j)/ntot);	
	}
      }
      entropy /= log((double) nbins);
      cout << entropy << " ";
      } else if (sarg=="-E") { 
      ColumnVector hist;
      int nbins=1000;
      double entropy=0;
      if (mask.nvoxels()<1) {
	generate_masks(mask,masknz,vol,lthr,uthr); 
	vol*=mask; 
      }
      hist = vol.histogram(nbins,masknz);
      double ntot = hist.Sum();
      for (int j=1; j<=nbins; j++) {
	if (hist(j)>0) {
	  entropy -= (hist(j)/ntot) * log(hist(j)/ntot);	
	}
      }
      entropy /= log((double) nbins);
      cout << entropy << " ";
    } else if (sarg=="-k") {
      narg++;
      volume4D<float> mask2;
      if (narg>=argc) {
        cerr << "Must specify an argument to -k" << endl;
        exit(2);
      }
      read_volume4D(mask,argv[narg]);
      if (!samesize(mask[0],vol[0])) {
        cerr << "Mask and image must be the same size" << endl;
        exit(3);
      }
      if (timeseriesMode && mask.tsize() != 1 ) { mask2=mask[timepoint]; mask=mask2; }
      if ( mask.tsize() > vol.tsize() ) {
        cerr << "Mask and image must be the same size" << endl;
        exit(3);
      }
      if ( mask.tsize() != vol.tsize() && mask.tsize() != 1) {
        // copy the last max volume until the correct size is reached
        while (mask.tsize() < vol.tsize() ) {
          mask.addvolume(mask[mask.maxt()]);
        }
      }

      mask.binarise(0.5);
      generateNonZeroMask(mask,masknz,vol);
        if (mask.tsize()!=1) vol*=mask;
        else vol*=mask[0];
    } else if (sarg=="-d") {
      narg++;
      if (narg>=argc) {
	cerr << "Must specify an argument to -d" << endl;
	exit(2);
      }
      volume4D<float> image2,image3;
      read_volume4D(image2,argv[narg]);
      if (!samesize(image2[0],vol[0])) {
	cerr << "Image used with -d must be the same size as the base image" << endl;
	exit(3);
      }
      if (timeseriesMode) { image3=image2[timepoint]; image2=image3; }
      if ( image2.tsize() > vol.tsize() ) {
	cerr << "Image must be the same size as the base image" << endl;
	exit(3);
      }
      if ( image2.tsize() != vol.tsize() && image2.tsize() != 1) {
	// copy the last max volume until the correct size is reached
	while (image2.tsize() < vol.tsize() ) {
   	  image2.addvolume(image2[image2.maxt()]);
        }
      }
      // now substract the new image from the base image 
      vol -= image2;
    } else if (sarg=="-l") {
      narg++;
      if (narg<argc) {
        lthr = atof(argv[narg]);
      } else {
	cerr << "Must specify an argument to -l" << endl;
	exit(2);
      }
      generate_masks(mask,masknz,vol,lthr,uthr);
      vol*=mask;
    } else if (sarg=="-u") {
      narg++;
      if (narg<argc) {
        uthr = atof(argv[narg]);
      } else {
	cerr << "Must specify an argument to -u" << endl;
	exit(2);
      }
      generate_masks(mask,masknz,vol,lthr,uthr);
      vol*=mask;
    } else if (sarg=="-a") {
      vol = abs(vol);
    } else if (sarg=="-v") {
      if (mask.nvoxels()>0) {
        long int nvox = mask.sum();
        if (mask.tsize() == 1) nvox = nvox * vol.tsize();
	cout << (long int) nvox << " " 
	     << nvox * vol.xdim() * vol.ydim() * vol.zdim() << " ";
      } else {
	cout << (long int) vol.nvoxels() * vol.tsize() << " "
	     << vol.nvoxels() * vol.tsize() * vol.xdim() * vol.ydim() * vol.zdim() << " ";
      }
    } else if (sarg=="-V") {
      if (masknz.nvoxels()>0) {
	cout << (long int) masknz.sum() << " " 
	     << masknz.sum() * vol.xdim() * vol.ydim() * vol.zdim() << " ";
      } else {
	long int nzvox;
	nzvox = nonzerocount(vol);
	cout << nzvox << " " 
	     << nzvox * vol.xdim() * vol.ydim() * vol.zdim() << " ";
      }
    } else if (sarg=="-D") {
	// hidden debug option!  // now -D not -d
      cout << vol.sum() << " ";
    } else if (sarg=="-s") {
	if (mask.nvoxels()>0) cout << vol.stddev(mask) << " ";
	else cout << vol.stddev() << " ";
    } else if (sarg=="-S") {
      if (masknz.nvoxels()>0) {
	cout << vol.stddev(masknz) << " ";
      } else {
	cout << nonzerostddev(vol) << " ";
      }
    } else if (sarg=="-r") {
      if (mask.nvoxels()>0) cout << vol.robustmin(mask) << " " << vol.robustmax(mask) << " ";
        else cout << vol.robustmin() << " " << vol.robustmax() << " ";
    } else if (sarg=="-R") {
	if (mask.nvoxels()>0) cout << vol.min(mask) << " " << vol.max(mask) << " ";
	else cout << vol.min() << " " << vol.max() << " ";
    } else if (sarg=="-c") {
	ColumnVector cog(4);
	// convert from fsl mm to voxel to nifti sform coord
	cog.SubMatrix(1,3,1,1) = vol[0].cog();
	cog(4) = 1.0;
	cog = vol[0].newimagevox2mm_mat() * cog; 
	cout << cog(1) << " " << cog(2) << " " << cog(3) << " " ;
    } else if (sarg=="-C") {
    ColumnVector cog(4);
	// convert from fsl mm to fsl voxel coord to nifti voxel coord
	cog.SubMatrix(1,3,1,1) = vol[0].cog();
	cog(4) = 1.0;
	cog = (vol[0].niftivox2newimagevox_mat()).i() * cog;
	cout << cog(1) << " " << cog(2) << " " << cog(3) << " " ;
    } else if (sarg=="-p") {
      float n;
      narg++;
      if (narg<argc) {
        n = atof(argv[narg]);
      } else {
	cerr << "Must specify an argument to -p" << endl;
	exit(2);
      }
      if ( (n<0) || (n>100) ) {
    	cerr << "Percentile must be between 0 and 100" << endl;
    	exit(1);
      }
      if (mask.nvoxels()>0) cout << vol.percentile((float) n/100.0, mask) << " ";
      else cout << vol.percentile((float) n/100.0) << " ";
    } else if (sarg=="-P") { 
      float n;
      narg++;
      if (narg<argc) {
        n = atof(argv[narg]);
      } else {
	cerr << "Must specify an argument to -P" << endl;
	exit(2);
      }
      if ( (n<0) || (n>100) ) {
    	cerr << "Percentile must be between 0 and 100" << endl;
    	exit(1);
      }
      if (mask.nvoxels()<1) {
	generate_masks(mask,masknz,vol,lthr,uthr); 
	vol*=mask; 
      }
      cout << vol.percentile((float) n/100.0,masknz) << " ";
    } else if (sarg=="-h") {
      float n;
      narg++;
      if (narg<argc) {
        n = atof(argv[narg]);
      } else {
	cerr << "Must specify the number of bins" << endl;
	exit(2);
      }
      int nbins = (int) n;
      if (nbins<1) {
    	cerr << "Must specify at least 1 bin" << endl;
    	exit(1);
      }
      if (mask.nvoxels()>0) {
	cout << vol.histogram(nbins,vol.min(),vol.max(),mask) << " ";
      } else {
	cout << vol.histogram(nbins,vol.min(),vol.max()) << " ";
      }
   } else if (sarg=="-H") {
      float n;
      narg++;
      if (narg<argc) {
        n = atof(argv[narg]);
      } else {
	cerr << "Must specify the number of bins" << endl;
	exit(2);
      }
      int nbins = (int) n;
      if (nbins<1) {
    	cerr << "Must specify at least 1 bin" << endl;
    	exit(1);
      }
      float min=0;
      narg++;
      if (narg<argc) {
        min = atof(argv[narg]);
      } else {
	cerr << "Must specify the histogram minimum intensity" << endl;
	exit(2);
      }
      float max=0;
      narg++;
      if (narg<argc) {
        max = atof(argv[narg]);
      } else {
	cerr << "Must specify the histogram maximum intensity" << endl;
	exit(2);
      }
      if (mask.nvoxels()>0) {
	cout << vol.histogram(nbins,min,max,mask) << " ";
      } else {
	cout << vol.histogram(nbins,min,max) << " ";
      }
    } else {
	cerr << "Unrecognised option: " << sarg << endl;
	exit(3);
    }
  
    narg++;
  }
   }
   cout << endl;
   
  }
  return 0;
}



int main(int argc,char *argv[])
{

  Tracer tr("main");
  string progname(argv[0]);
  bool timeseriesMode(false);
  string indexMask("");
  while ( argc > 2 && ( string(argv[1])=="-t" || string(argv[1]) =="-K" ) ) {
    if ( string(argv[1])=="-t" )
      timeseriesMode=true;
    if ( string(argv[1])=="-K" ) {
      indexMask=string(argv[2]);
      argv++;
      argc--;
    }
    argv++;
    argc--;
  }
  if (argc < 3 )
    return print_usage(progname);  
  try {
    return fmrib_main_float(argc,argv,timeseriesMode,indexMask);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } catch (...) {
    // do nothing - just exit without garbage message
  }

  return -1;
  
}

