/*  robustfov.cc

    Mark Jenkinson and Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

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

// Calculates a robust FOV to help with brain extraction

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

#include <stdlib.h>
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "utils/options.h"

using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;

// The two strings below specify the title and example usage that is
//  printed out as the help or usage message

string title="robustfov \nCopyright(c) 2012, University of Oxford (Mark Jenkinson)\n\nReduce FOV of image to remove lower head and neck.";
string examples="robustfov [options] -i <image> -r <outputimage>\nrobustfov [options] -i <image>";

// Each (global) object below specificies as option and can be accessed
//  anywhere in this file (since they are global).  The order of the
//  arguments needed is: name(s) of option, default value, help message,
//       whether it is compulsory, whether it requires arguments
// Note that they must also be included in the main() function or they
//  will not be active.

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<bool> debug(string("--debug"), false,
		  string("turn on debugging output"),
		  false, no_argument);
Option<float> brainsize(string("-b"), 170.0f,
		  string("size of brain in z-dimension (default 170mm)"),
		  false, requires_argument);
Option<string> roivol(string("-r"), string(""),
		  string("ROI volume output name"),
		  false, requires_argument);
Option<string> matname(string("-m"), string(""),
		  string("matrix output name (roi to full fov)"),
		  false, requires_argument);
Option<string> inname(string("-i"), string(""),
		  string("input filename"),
		  true, requires_argument);
int nonoptarg;

////////////////////////////////////////////////////////////////////////////

// Local functions

ColumnVector calc_FOV(int q, int brainlen, int xmax, int ymax, int zmax, int dim) {
  ColumnVector FOV(6);
  // set default FOV as whole image
  FOV=0;
  FOV(2) = xmax+1;
  FOV(4) = ymax+1;
  FOV(6) = zmax+1;
  // adjust for the case where adding the brainlen would go beyond the image FOV
  int adjust=0;
  if (q-brainlen<0) adjust=brainlen-q;
  if (dim==1) { FOV(1) = Max(q-brainlen,0); FOV(2) = brainlen-adjust; }
  if (dim==2) { FOV(3) = Max(q-brainlen,0); FOV(4) = brainlen-adjust; }
  if (dim==3) { FOV(5) = Max(q-brainlen,0); FOV(6) = brainlen-adjust; }
  // for the negative dims, q=max corresponds to coord=0 (as mapping q:max->0 is actually coord:0->max)
  if (dim==-1) { FOV(1) = xmax-q; FOV(2) = brainlen-adjust; }
  if (dim==-2) { FOV(3) = ymax-q; FOV(4) = brainlen-adjust; }
  if (dim==-3) { FOV(5) = zmax-q; FOV(6) = brainlen-adjust; }
  if (debug.value()) { cout << "Internal FOV is: " << FOV.t() << endl; }
  return FOV;
}

void convert_fov_coord2nifti(int& minx, int& maxx, int& miny, int& maxy, int& minz, int& maxz, 
			     const volume<float>& input_vol)
{
  // now transform these coordinates to nifti conventions
  ColumnVector v(4);
  v << minx << miny << minz << 1.0;
  v = input_vol.niftivox2newimagevox_mat().i() * v;
  minx = MISCMATHS::round(v(1));
  miny = MISCMATHS::round(v(2));
  minz = MISCMATHS::round(v(3));
  v << maxx << maxy << maxz << 1.0;
  v = input_vol.niftivox2newimagevox_mat().i() * v;
  maxx = MISCMATHS::round(v(1));
  maxy = MISCMATHS::round(v(2));
  maxz = MISCMATHS::round(v(3));
  if (minx>maxx) { int tmp=minx; minx=maxx; maxx=tmp; }
  if (miny>maxy) { int tmp=miny; miny=maxy; maxy=tmp; }
  if (minz>maxz) { int tmp=minz; minz=maxz; maxz=tmp; }
}

template <class T>
bool is_int_image(const volume<T>& im)
{
  if ((im.max()>im.min()) && (im.max()-im.min()<1.0)) return false;  // small floating-point range
  double fracsum=0.0;
  float val=0.0, frac=0.0;
  int nvox=0;
  for (int z=im.minz(); z<=im.maxz(); z++) {
    for (int y=im.miny(); y<=im.maxy(); y++) {
      for (int x=im.minx(); x<=im.maxx(); x++) {
	val=im(x,y,z);
	if (val>1e-8) {
	  nvox++;
	  frac=fabs(val - MISCMATHS::round(val));
	  fracsum+= frac;
	}
      }
    }
  }
  if (nvox>0) fracsum /= nvox;
  return (fracsum<1e-5);   // somewhat arbitrary threshold, but non-integer images should very rarely pass this
}

////////////////////////////////////////////////////////////////////////////

// do main work

int do_work(int argc, char* argv[]) 
{
  volume<float> vol;
  read_volume(vol,inname.value());
  if (vol.min()>0) { // cope with images set in a negative range
    volume<float> mask(vol);
    mask.binarise(0,vol.max()+1,exclusive);
    vol -= vol.min();  // make everything non-negative
    vol *= mask;  // keep previous zeros as zeros
  }
  bool is_int_im = is_int_image(vol);
  int xmax, ymax, zmax;
  xmax=vol.xsize()-1;
  ymax=vol.ysize()-1;
  zmax=vol.zsize()-1;

  int dim=3;
  Matrix sqform;
  sqform=vol.newimagevox2mm_mat();  // TODO - need to deal differently when qform/sform are not set
  sqform=sqform.SubMatrix(1,3,1,3);
  ColumnVector zhat(3), zvec;
  zhat << 0.0 << 0.0 << 1.0;
  zvec = sqform.i() * zhat;
  zvec(1) *= vol.xdim();
  zvec(2) *= vol.ydim();
  zvec(3) *= vol.zdim();
  if ( (fabs(zvec(3))>fabs(zvec(2))) && (fabs(zvec(3))>fabs(zvec(1))) ) { if (zvec(3)>0) dim=3; else dim=-3; }
  if ( (fabs(zvec(2))>fabs(zvec(3))) && (fabs(zvec(2))>fabs(zvec(1))) ) { if (zvec(2)>0) dim=2; else dim=-2; }
  if ( (fabs(zvec(1))>fabs(zvec(3))) && (fabs(zvec(1))>fabs(zvec(2))) ) { if (zvec(1)>0) dim=1; else dim=-1; }
  if (dim==0) { cerr << "Cannot determine direction of S-I" << endl; return 1; }
  if (verbose.value()) { cout << "Dimension chosen for Superior is: " << dim << endl; }

  int qstart=0, qend=0, qinc=0, nslicevox=0;
  float qdim=0;
  if (dim==1) { qstart=xmax; qend=-1; qinc=-1; qdim=vol.xdim(); nslicevox=vol.ysize()*vol.zsize(); }
  if (dim==2) { qstart=ymax; qend=-1; qinc=-1; qdim=vol.ydim(); nslicevox=vol.xsize()*vol.zsize(); }
  if (dim==3) { qstart=zmax; qend=-1; qinc=-1; qdim=vol.zdim(); nslicevox=vol.xsize()*vol.ysize(); }
  if (dim==-1) { qstart=0; qend=xmax+1; qinc=1; qdim=vol.xdim(); nslicevox=vol.ysize()*vol.zsize(); }
  if (dim==-2) { qstart=0; qend=ymax+1; qinc=1; qdim=vol.ydim(); nslicevox=vol.xsize()*vol.zsize(); }
  if (dim==-3) { qstart=0; qend=zmax+1; qinc=1; qdim=vol.zdim(); nslicevox=vol.xsize()*vol.ysize(); }


  // save thresholded voxel counts for slices - done from upper end of array (as dim=3 is "natural" counting)
  int rowmax=abs(qend-qstart)+1;
  ColumnVector frac(rowmax), nzsum(rowmax), threshvals(rowmax);
  { 
    volume<float> copyv;
    for (int loopnum=1; loopnum<=2; loopnum++) {
      copyv.deactivateROI();
      copyv=vol;
      if (loopnum==1) copyv.binarise(0,copyv.max()+1,exclusive);  // set to 1 anything above 0
      copyv.activateROI();
      for (int q=qstart, row=rowmax; q!=qend; q+=qinc, row--) {
	if (abs(dim)==1) {
	  copyv.setROIlimits(q,0,0,q,vol.ysize()-1,vol.zsize()-1);
	} else if (abs(dim)==2) {
	  copyv.setROIlimits(0,q,0,vol.xsize()-1,q,vol.zsize()-1);
	} else if (abs(dim)==3) {
	  copyv.setROIlimits(0,0,q,vol.xsize()-1,vol.ysize()-1,q);
	}
	float p1, p99, thresh;
	if (loopnum==1) { nzsum(row)=copyv.sum(); }
	if (loopnum==2) {
	  if (nzsum(row)>0) {
	    p1=copyv.percentile(1-0.99*nzsum(row)/nslicevox);
	    p99=copyv.percentile(1-0.01*nzsum(row)/nslicevox);
	    thresh=0.1*(p99-p1)+p1;  // 10% of the value from 1st percentile to 99th percentile
	    copyv.binarise(thresh);
	    frac(row)=copyv.sum()/nzsum(row);   // percentage of super-threshold voxels (of non-zero voxels)
	    if ((is_int_im) && (p99-p1<5)) {  // check for very low dynamic range of int (indicates heavy quantisation)
	      frac(row)=0.85;  // use value appropriate to pure Gaussian noise
	    }
	    threshvals(row)=thresh;
	  } else {
	    frac(row)=1;
	    threshvals(row)=0;
	  }
	  if (verbose.value()) { cout << "At q="<<q<<" : frac = " << frac(row) << " : thresh = " << threshvals(row) << endl; }
	}
      }
    }
  }
  
  // now search for point where thresholded voxels dip significantly (minchange) as noisy slices
  //   contain a large number of super-threshold voxels but slices near the top of the brain contain
  //   far fewer  (if no dip is found then assume the top of the head/brain was already in the top slice)
  ColumnVector fov(6);
  float minchange=0.3;  // pick 30% : TODO - calculate this from a mm^2 area of the top of the brain
  int gap=MISCMATHS::round(3.0/qdim+0.5);  // 3mm gap  (distance between baseline slice and slice being tested)
  int minlen=gap;
  int brainlen=MISCMATHS::round(brainsize.value()/qdim);  // 150mm of brain + top scalp
  int minvox=0;    // only start counting slices with at least minvoxfraction of non-zero voxels (avoids blank areas due to gradient unwarping or similar)
  float minvoxfraction=0.02;  // 2% of non-zero voxels  (equivalent of 14% by 14% of coverage in a slice)
  if (abs(dim)==1) minvox=MISCMATHS::round(ymax*zmax*minvoxfraction);
  if (abs(dim)==2) minvox=MISCMATHS::round(xmax*zmax*minvoxfraction);
  if (abs(dim)==3) minvox=MISCMATHS::round(xmax*ymax*minvoxfraction);
  float baseval;
  int qmax=Max(qstart,qend-1), qbest=-1;
  bool foundbest=false;
  for (int q=qmax-gap; q>=brainlen*2/3; q--) {
    if (!foundbest) {
      baseval=frac(q+gap+1);
      if ( (nzsum(q+1)>minvox) && (frac(q+1)<baseval-minchange) ) { 
	//if ( (nzsum(q+gap+1)>minvox) && (nzsum(q+1)>minvox) && (frac(q+1)<baseval-minchange) ) { 
	// TODO add a test for threshold values? (avoid artefact being identified as brain in superior slices ala Marco's data)
	//   something like threshold must be > 0.5 * max threshold found   (but don't want bias field to stop proper slices being found!)
	if (verbose.value()) { cout << "Testing stability of q = " << q << " versus value " << baseval - minchange << endl; }
	// check that the frac isn't a fluke and is no worse for a continuous stretch (minlen)
	bool stable=true;
	for (int q0=q; q0>=q-minlen; q0--) {
	  if (frac(q0+1)>=baseval-minchange) stable=false;
	}
	if (stable) {
	  if (verbose.value()) { cout << "Chosen q = " << q << endl; }
	  qbest=q;
	  foundbest=true; // effectively break from loop (but retain this info)
	}
      }
    }
  }
  if (!foundbest) qbest=qmax;  // if nothing found, assume that there is brain in the top slice
  if (verbose.value()) { cout << "Found best change at: " << qbest << endl; }
  fov=calc_FOV(qbest,brainlen,xmax,ymax,zmax,dim);

  // generate a FLIRT matrix to allow user to go backwards and forwards between FOV and full image
  // NB: do this *before* converting to nifti coords, as FLIRT uses newimage coords
  Matrix aff(4,4);
  aff=IdentityMatrix(4);
  aff(1,4)=fov(1)*vol.xdim();
  aff(2,4)=fov(3)*vol.ydim();
  aff(3,4)=fov(5)*vol.zdim();
  if (matname.set()) { write_ascii_matrix(aff,matname.value()); }

  // generate ROI volume output
  if (roivol.set()) {
    vol.activateROI();
    vol.setROIlimits(MISCMATHS::round(fov(1)),MISCMATHS::round(fov(3)),MISCMATHS::round(fov(5)),
		     MISCMATHS::round(fov(1)+fov(2))-1,MISCMATHS::round(fov(3)+fov(4))-1,MISCMATHS::round(fov(5)+fov(6))-1);
    volume<float> roiv;
    roiv = vol.ROI();
    save_volume(roiv,roivol.value());
  }

  // convert to nifti conventions
  int minx, miny, minz, maxx, maxy, maxz;
  minx=MISCMATHS::round(fov(1));
  miny=MISCMATHS::round(fov(3));
  minz=MISCMATHS::round(fov(5));
  maxx=MISCMATHS::round(fov(2))+minx-1;
  maxy=MISCMATHS::round(fov(4))+miny-1;
  maxz=MISCMATHS::round(fov(6))+minz-1;
  convert_fov_coord2nifti(minx,maxx,miny,maxy,minz,maxz,vol);
  fov(1)=minx; fov(3)=miny; fov(5)=minz;
  fov(2)=maxx-minx+1;
  fov(4)=maxy-miny+1;
  fov(6)=maxz-minz+1;
  cout << "Final FOV is: " << endl << fov.t() << endl;
  
  return 0;
}

////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{

  Tracer tr("main");
  OptionParser options(title, examples);

  try {
    // must include all wanted options here (the order determines how
    //  the help message is printed)
    options.add(inname);
    options.add(brainsize);
    options.add(matname);
    options.add(roivol);
    options.add(debug);
    options.add(verbose);
    options.add(help);
    
    nonoptarg = options.parse_command_line(argc, argv);

    // line below stops the program if the help was requested or 
    //  a compulsory option was not set
    if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
      {
	options.usage();
	exit(EXIT_FAILURE);
      }
    
  }  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } 

  // Call the local functions

  return do_work(argc,argv);
}

