/*  pnm_evs.cc  

    Mark Jenkinson, FMRIB Image Analysis Group

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

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "utils/options.h"

using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;

// The two strings below specify the title and example usage that is
//  printed out as the help or usage message

string title="pnm_evs (Version 2.0)\nCopyright(c) 2011, University of Oxford (Mark Jenkinson)";
string examples="pnm_evs [options] --tr=3.0 -i -o -r -c ... TBD ...";

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
Option<int> cardorder(string("--oc"), 2,
		  string("order of basic cardiac regressors (number of Fourier pairs) - default=2"),
		  false, requires_argument);
Option<int> resporder(string("--or"), 1,
		  string("order of basic respiratory regressors (number of Fourier pairs) - default=1"),
		  false, requires_argument);
Option<int> multc(string("--multc"), 0,
		  string("order of multiplicative cardiac terms (also need to set --multr) - default=0"),
		  false, requires_argument);
Option<int> multr(string("--multr"), 0,
		  string("order of multiplicative respiratory terms (also need to set --multc) - default=0"),
		  false, requires_argument);
Option<float> tr(string("--tr"), 0.0,
		  string("TR of fMRI volumes (in seconds)"),
		  true, requires_argument);
Option<float> heartratesmooth(string("--heartratesmooth"), 15.0,
		  string("Optional smoothing of heart rate regressor (in seconds - e.g. 10)"),
		  false, requires_argument);
Option<float> rvtsmooth(string("--rvtsmooth"), 0.0,
		  string("Optional smoothing of RVT regressor (in seconds - default 0)"),
		  false, requires_argument);
Option<string> slicedir(string("--slicedir"), string("z"),
		  string("specify slice direction (x/y/z) - default is z"),
		  false, requires_argument);
Option<string> sliceordering(string("--sliceorder"), string("up"),
		  string("specify slice ordering (up/down/interleaved_up/interleaved_down)"),
		  false, requires_argument);
Option<string> slicetiming(string("--slicetiming"), string(""),
		  string("specify slice timing via an external file"),
		  false, requires_argument);
Option<string> csfmask(string("--csfmask"), string(""),
		  string("filename of csf mask image (and generate csf regressor)"),
		  false, requires_argument);
Option<string> cardname(string("-c,--cardiac"), string(""),
		  string("input filename for cardiac values (1 or 2 columns: time [phase])"),
		  false, requires_argument);
Option<string> respname(string("-r,--respiratory"), string(""),
		  string("input filename for respiratory phase values (1 or 2 columns: time [phase])"),
		  false, requires_argument);
Option<string> rvt(string("--rvt"), string(""),
		   string("input filename of RVT data (2 columns: time value)"),
		   false, requires_argument);
Option<string> heartrate(string("--heartrate"), string(""),
			 string("input filename for heartrate data (2 columns: time value)"),
			 false, requires_argument);
Option<string> inname(string("-i,--in"), string(""),
		  string("input image filename (of 4D functional/EPI data)"),
		  true, requires_argument);
Option<string> outname(string("-o,--out"), string(""),
		  string("output filename (for confound/EV matrix)"),
		  true, requires_argument);
int nonoptarg;

// GLOBAL VARIABLES (!)

Matrix slicetimingvals;

////////////////////////////////////////////////////////////////////////////

// Local functions

int sanity_check()
{
  // sanity checking
  if (cardorder.value()<0) { 
    cerr << "Invalid order for cardiac (" << cardorder.value() << ") - must be non-negative" << endl; 
    return 1;
  }
  if ((cardorder.value()>0) && cardname.unset()) {
    cerr << "Must specify cardiac phase in order to generate cardiac regressors" << endl; 
    return 1;
  }
  if (heartrate.set() && cardname.unset()) { 
    cerr << "Must specify cardiac phase in order to generate heartrate regressor" << endl; 
    return 1;
  }
  if (resporder.value()<0) { 
    cerr << "Invalid order for respiratory (" << resporder.value() << ") - must be non-negative" << endl;    return 1;
  }
  if ((resporder.value()>0) && respname.unset()) {
    cerr << "Must specify respiratory phase in order to generate respiratory regressors" << endl; 
    return 1;
  }
  if (tr.value()<=0.0) {
    cerr << "Must specify a positive TR value (specified value = " << tr.value() << ")" << endl;
    return 1;
  }
  return 0;
}

int numslices(int64_t sx, int64_t sy, int64_t sz) {
  int nslices=0;
  if (slicedir.value()=="x") nslices=sx;
  if (slicedir.value()=="y") nslices=sy;
  if (slicedir.value()=="z") nslices=sz;
  if (nslices<=0) { 
    cerr << "Cannot determine number of slices: slicedir = " << slicedir.value() << endl; 
    exit(EXIT_FAILURE); 
  }
  return nslices;
}


int moving_average(ColumnVector& vec, int ns) {
  // use the efficient circular storage (one in, one out) method
  ColumnVector tmp(ns);
  tmp=0.0;
  double movsum=0.0;
  int m=1;  // index for circular tmp buffer
  int n2=(ns-1)/2, den=0;
  for (int n=1; n<=vec.Nrows()+n2; n++) {
    if (n>ns) { movsum-=tmp(m); den--; }
    if (n<=vec.Nrows()) {
      tmp(m)=vec(n);
      movsum+=tmp(m);  den++;
    }
    if ((n-n2)>=1) vec(n-n2)=movsum/den;
    m++;  if (m>ns) m=1;
  }
  return 0;
}

ColumnVector interp1(const ColumnVector& x, const ColumnVector& y, const ColumnVector& xi) 
// Look-up function for data table defined by x, y
// Returns the values yi at xi using linear interpolation
// Assumes that x is sorted in ascending order
// Also assumes that xi is sorted in ascending order
{
  ColumnVector ans;
  ans=xi*0.0;
  int ind=2;
  float xmax, xmin;
  xmax = x.Maximum();
  xmin = x.Minimum();
  for (int idx=1; idx<=xi.Nrows(); idx++) {
    if(xi(idx) >= xmax) {
      ind=x.Nrows();
    } else if(xi(idx) <= xmin) {
      ind=2;
    } else {
      while(xi(idx) >= x(ind)) { ind++; }
    }
    float xa = x(ind-1), xb = x(ind), ya = y(ind-1), yb = y(ind);
    ans(idx) = ya + (xi(idx) - xa)/(xb - xa) * (yb - ya);
  }
  return ans;
}

ColumnVector bandpass_temporal_filter(const ColumnVector& invals, double lp_sigma) {
  ColumnVector outvals(invals.Nrows());
  volume4D<float> tmpv(1,1,1,invals.Nrows());
  for (int z=0; z<invals.Nrows(); z++) { tmpv(0,0,0,z)=invals(z+1); }
  tmpv = bandpass_temporal_filter(tmpv,0.0,lp_sigma);
  for (int z=0; z<invals.Nrows(); z++) { outvals(z+1)=tmpv(0,0,0,z); }
  return outvals;
}

int bandpass_temporal_filter(volume<float>& invals, int x, int y, double lp_sigma) {
  volume4D<float> tmpv(1,1,1,invals.zsize());
  for (int z=0; z<invals.zsize(); z++) { tmpv(0,0,0,z)=invals(x,y,z); }
  tmpv = bandpass_temporal_filter(tmpv,0.0,lp_sigma);
  for (int z=0; z<invals.zsize(); z++) { invals(x,y,z)=tmpv(0,0,0,z); }
  return 0;
}


//////////////////////////////////////////////////////////////////////////
// SLICE ORDERING INFORMATION
//////////////////////////////////////////////////////////////////////////
// Slice ordering information: ascending => acquisition starts with slices closest to the
// feet, and moves upwards towards the head. descending => topmost (closest to the head) slices
// acquired first, then slices further towards the feet. interleaved => odd slices before even
// ones, with the most inferior (bottom most) odd numbered slice acquired first, i.e. like
// an ascending acquisition. NOTE!!!! a) obviously this applies primarily to axial or "transversal"
// acquisitions, b) that *irrespective* of the slice acquisition ordering, the slices are numbered
// and stored on the scanner with the bottom-most slice as number 1 and the topmost slice as number N.
//
//
// SIEMENS WEIRDNESS!!!! If there are an odd number of slices and ordering
// is interleaved then scanner will acquire ODD then EVEN. HOWEVER, if there
// are and even number of slices it acquires EVEN then ODD. Go figure!

// Generate timing for slice #slicenum  (where slicenum starts at 1, not 0)
ColumnVector slice_timing(int slicenum, int64_t sx, int64_t sy, int64_t sz, int64_t st)
{
  ColumnVector stimes(st);
  if (slicetiming.unset()) {
    int nslices=numslices(sx,sy,sz);
    
    float firsttime=0.0, halftime=0.0;
    if (sliceordering.value()=="up") {
      firsttime = (slicenum-1)*tr.value()/nslices;
    } else if (sliceordering.value()=="down") {
      firsttime = (nslices - slicenum)*tr.value()/nslices;
    } else if (sliceordering.value()=="interleaved_up") {
      if ((nslices % 2) == 1) { // odd number of slices - odd first, then even
	if ((slicenum % 2) == 1) { firsttime = ((slicenum-1)/2)*tr.value()/nslices; }
	else { halftime = ((nslices-1)/2)*tr.value()/nslices + tr.value()/nslices;
	  firsttime = ((slicenum-2)/2)*tr.value()/nslices + halftime;}
      } else { // even number of slices - even first, then odd
	if ((slicenum % 2) == 0) { firsttime = ((slicenum-2)/2)*tr.value()/nslices; }
	else { halftime = ((nslices-2)/2)*tr.value()/nslices + tr.value()/nslices;
	  firsttime = ((slicenum-1)/2)*tr.value()/nslices + halftime; }
      }
    } else if (sliceordering.value()=="interleaved_down") {
      int snum = nslices - slicenum + 1;  // convert to pseudo-ascending number
      if ((nslices % 2) == 1) { // odd number of slices - odd first, then even
	if ((slicenum % 2) == 1) { firsttime = ((snum-1)/2)*tr.value()/nslices; }
	else { halftime = ((nslices-1)/2)*tr.value()/nslices + tr.value()/nslices;
	  firsttime = ((snum-2)/2)*tr.value()/nslices + halftime;}
      } else { // even number of slices - even first, then odd
	if ((slicenum % 2) == 0) { firsttime = ((snum-1)/2)*tr.value()/nslices; }
	else { halftime = ((nslices-2)/2)*tr.value()/nslices + tr.value()/nslices;
	  firsttime = ((snum-2)/2)*tr.value()/nslices + halftime; }
      }
    } else {
	cerr << "Unrecognised slice ordering: " << sliceordering.value() << endl;
	exit(EXIT_FAILURE);
      } 
    for (int n=1; n<=stimes.Nrows(); n++) {
      stimes(n) = firsttime + (n-1)*tr.value();
    }
  } else {
    // user provided timings
    ColumnVector usertimes = slicetimingvals.Column(Min(slicenum,slicetimingvals.Ncols()));
    if (usertimes.Ncols()==stimes.Nrows()) { stimes=usertimes; }
    else if (usertimes.Ncols()==1) {
      for (int n=1; n<=stimes.Nrows(); n++) {
	stimes(n) = usertimes(1) + (n-1)*tr.value();
      }
    }

  }
  return stimes;
}


volume<float> calc_confounds(const Matrix& cardph, const Matrix& respph, 
			     int64_t sx, int64_t sy, int64_t sz, int64_t st)
{
  if (verbose.value()) { cout << "Beginning to calculate confounds" << endl; }
  int nreg = cardorder.value()*2 + resporder.value()*2 + multc.value()*multr.value()*4;
  if (rvt.set()) nreg++;
  if (heartrate.set()) nreg++;
  if (csfmask.set()) nreg++;

  if (verbose.value()) { cout << "Number of confounds = " << nreg << endl; }

  int nslices = numslices(sx,sy,sz);
  if (verbose.value()) { cout << "Number of slices = " << nslices << endl; }
  volume<float> confoundvol(nreg,nslices,st);

  volume4D<float> vol;
  volume<float> csfm, varim, totalmask;
  if (csfmask.set()) {
    read_volume4D(vol,inname.value());
    read_volume(csfm,csfmask.value());
    varim = variancevol(vol);
    varim *= csfm;
    totalmask = csfm*0.0f;
  }

  // following few lines are just to setup tmpv (to be more efficient)
  ColumnVector tslice;
  tslice = slice_timing(1,sx,sy,sz,st);
  int nt=tslice.Nrows();

  for (int n=1; n<=nslices; n++) {
    if (verbose.value()) { cout << "Loop number " << n << endl; }
    // get slice timings
    tslice = slice_timing(n,sx,sy,sz,st);
    if (debug.value() && (n==1)) { write_ascii_matrix(tslice,"tslice"); }
    
    // resample phase values at the appropriate slice timings
    if (verbose.value()) { cout << "Resampling phase" << endl; }
    ColumnVector cph_slice(nt), rph_slice(nt);
    cph_slice = 0.0f;
    rph_slice = 0.0f;
    if (cardorder.value()>0) { cph_slice = interp1(cardph.Column(1),cardph.Column(2),tslice); }
    if (resporder.value()>0) { rph_slice = interp1(respph.Column(1),respph.Column(2),tslice); }
    if (debug.value() && (n==1)) { write_ascii_matrix(cph_slice,"cph_slice"); }
    if (debug.value() && (n==1)) { write_ascii_matrix(rph_slice,"rph_slice"); }
    
    // generate the required regressors
    int col=1;

    // cardiac-only terms
    if (verbose.value()) { cout << "Cardiac terms" << endl; }
    for (int m=1; m<=cardorder.value(); m++) {
      for (int nn=1; nn<=nt; nn++) {
	confoundvol(col-1,n-1,nn-1) = sin(cph_slice(nn)*m);
	confoundvol(col,n-1,nn-1) = cos(cph_slice(nn)*m);
      }
      col+=2;
    }

    if (verbose.value()) { cout << "Respiratory terms" << endl; }
    // respiratory-only terms
    for (int m=1; m<=resporder.value(); m++) {
      for (int nn=1; nn<=nt; nn++) {
	confoundvol(col-1,n-1,nn-1) = sin(rph_slice(nn)*m);
	confoundvol(col,n-1,nn-1) = cos(rph_slice(nn)*m);
      }
      col+=2;
    }

    if (verbose.value()) { cout << "Amplitude modulated terms" << endl; }
    // amplitude modulated terms
    for (int mc=1; mc<=multc.value(); mc++) {
      for (int mr=1; mr<=multr.value(); mr++) {
	for (int nn=1; nn<=nt; nn++) {
	  confoundvol(col-1,n-1,nn-1) = cos(cph_slice(nn)*mc + rph_slice(nn)*mr);
	  confoundvol(col,n-1,nn-1) = sin(cph_slice(nn)*mc + rph_slice(nn)*mr);
	  confoundvol(col+1,n-1,nn-1) = cos(cph_slice(nn)*mc - rph_slice(nn)*mr);
	  confoundvol(col+2,n-1,nn-1) = sin(cph_slice(nn)*mc - rph_slice(nn)*mr);
	}
	col+=4;
      }
    }

    if (heartrate.set()) {  // heartrate
      if (verbose.value()) { cout << "HR term" << endl; }
      Matrix hrmat;
      hrmat = read_ascii_matrix(heartrate.value());
      if (heartratesmooth.value()>0.0) {  // Smooth with moving average over higher sampling (~10Hz)
	int nsec = MISCMATHS::round((tslice(2)-tslice(1))/0.1);   // divide TR into N sections (close to 0.1s)
	float tsamp=(tslice(2)-tslice(1))/nsec;
	ColumnVector tvals(nt*nsec);
	for (int nn=1; nn<=nt*nsec; nn++) {
	  tvals(nn) = interp1(hrmat.Column(1),hrmat.Column(2),tslice(1)+(nn-1)*tsamp); 
	}
	// smooth
	int nsmooth  = MISCMATHS::round(heartratesmooth.value()/tsamp);
	moving_average(tvals,nsmooth);
	// resample
	for (int nn=1; nn<=nt; nn++) {
	  confoundvol(col-1,n-1,nn-1) = tvals(1+(nn-1)*nsec); 
	}
      } else {
	for (int nn=1; nn<=nt; nn++) {
	  confoundvol(col-1,n-1,nn-1) = interp1(hrmat.Column(1),hrmat.Column(2),tslice(nn)); 
	}
      }
      col++;
    }

    if (rvt.set()) {  // rvt
      if (verbose.value()) { cout << "RVT term" << endl; }    // MJ CHECK - SEEMS TO BE TOO LITTLE CHANGE WRT SLICE
      Matrix rvtmat;
      rvtmat = read_ascii_matrix(rvt.value());
      ColumnVector rvt_slicevals(nslices);
      if (rvtsmooth.value()>0.0) {  // Smooth with moving average over higher sampling (~10Hz)
	int nsec = MISCMATHS::round((tslice(2)-tslice(1))/0.1);   // divide TR into N sections (close to 0.1s)
	float tsamp=(tslice(2)-tslice(1))/nsec;
	ColumnVector tvals(nt*nsec);
	for (int nn=1; nn<=nt*nsec; nn++) {
	  tvals(nn) = interp1(rvtmat.Column(1),rvtmat.Column(2),tslice(1)+(nn-1)*tsamp); 
	}
	// smooth
	int nsmooth  = MISCMATHS::round(rvtsmooth.value()/tsamp);
	moving_average(tvals,nsmooth);
	// resample
	for (int nn=1; nn<=nt; nn++) {
	  confoundvol(col-1,n-1,nn-1) = tvals(1+(nn-1)*nsec); 
	}
      } else {
	for (int nn=1; nn<=nt; nn++) {
	  confoundvol(col-1,n-1,nn-1) = interp1(rvtmat.Column(1),rvtmat.Column(2),tslice(nn)); 
	}
      }
      col++;
    }

    if (csfmask.set()) {  // csf
      if (verbose.value()) { cout << "CSF term" << endl; }
      // take top 10% of variance voxels and generate mean timeseries from this
      // Define a single slice mask
      volume<float> slicemask(varim);
      slicemask=0.0f;
      for (int x=slicemask.minx(); x<=slicemask.maxx(); x++) {
	for (int y=slicemask.miny(); y<=slicemask.maxy(); y++) {
	  slicemask(x,y,n-1)=1.0f;
	}
      }
      slicemask *= csfm;
      float varthresh = varim.percentile(0.9f,slicemask);
      slicemask *= varim;  // put the variance values into this volume
      slicemask.binarise(varthresh);
      totalmask += slicemask;
      for (int mm=1; mm<=nt; mm++) {
	confoundvol(col-1,n-1,mm-1) = mean(vol[mm-1],slicemask);
      }
      col++;
    }

  }

  if (csfmask.set()) {
    save_volume(totalmask,outname.value()+"_csfmask");
  }

  return confoundvol;
}


int do_work(int argc, char* argv[]) 
{
  if (sanity_check()!=0) { exit(EXIT_FAILURE); }

  if (verbose.value()) { cout << "Calculating slice timings" << endl; }
  if (slicetiming.set()) {
    slicetimingvals = read_ascii_matrix(slicetiming.value());
  }

  // read in cardiac phase values
  Matrix cardph;
  if (cardname.set()) {
    if (verbose.value()) { cout << "Reading in cardiac phase" << endl; }
    Matrix cardpht;
    cardpht = read_ascii_matrix(cardname.value());
    if (cardpht.Ncols()==1) {
      cardph.ReSize(cardpht.Nrows(),2);
      for (int n=1; n<=cardpht.Nrows(); n++) { cardph(n,1)=cardpht(n,1); cardph(n,2)=(n-1)*2.0*M_PI; }
    } else {
      cardph=cardpht;
    }
  }

  // read in respiratory phase values
  Matrix respph;
  if (respname.set()) {
    if (verbose.value()) { cout << "Reading in respiratory phase" << endl; }
    respph = read_ascii_matrix(respname.value());
    // decide whether input is just timing of zero phase points, or full histogram normalised phases
    if (respph.Ncols()==1) {
      Matrix resppht;
      resppht=respph;
      respph.ReSize(resppht.Nrows(),2);
      for (int n=1; n<=resppht.Nrows(); n++) { respph(n,1)=resppht(n,1); respph(n,2)=(n-1)*2.0*M_PI; }
    }
  }

  if (verbose.value()) { cout << "Reading in functional image" << endl; }
  int64_t vsx,vsy,vsz,vst,vs5;
  read_volume_size(inname.value(),vsx,vsy,vsz,vst,vs5);

  int nslices, nvols;
  nvols = vst;
  nslices = numslices(vsx,vsy,vsz);
  
  // confound storage
  volume<float> confoundvol;
  if (verbose.value()) { cout << "Calculating confounds" << endl; }
  confoundvol = calc_confounds(cardph,respph,vsx,vsy,vsz,vst);
  int numevs = confoundvol.xsize();

  // save output
  if (verbose.value()) { cout << "Calculating and saving EVs" << endl; }
  double tmean=0;
  string outbase=fslbasename(outname.value());
  int nnx=1, nny=1, nnz=1;
  if (slicedir.value()=="x") nnx=nslices;
  if (slicedir.value()=="y") nny=nslices;
  if (slicedir.value()=="z") nnz=nslices;
  volume4D<float> ev_vol(nnx,nny,nnz,nvols);
  for (int ev=1; ev<=numevs; ev++) {
    for (int ns=1; ns<=nslices; ns++) {
      tmean=0;
      for (int nv=1; nv<=nvols; nv++) {
	ev_vol(Min(ns,nnx)-1,Min(ns,nny)-1,Min(ns,nnz)-1,nv-1)=confoundvol(ev-1,ns-1,nv-1);
	tmean+=ev_vol(Min(ns,nnx)-1,Min(ns,nny)-1,Min(ns,nnz)-1,nv-1);
      }
      // demean EV
      tmean/=nvols;
      for (int nv=1; nv<=nvols; nv++) {
	ev_vol(Min(ns,nnx)-1,Min(ns,nny)-1,Min(ns,nnz)-1,nv-1) -= tmean;
      }      
    }
    
    save_volume4D(ev_vol,outbase+"ev"+num2str(ev,3));
  }

  return 0;

}

////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{


  Tracer tr_main("main");
  OptionParser options(title, examples);

  try {
    // must include all wanted options here (the order determines how
    //  the help message is printed)
    options.add(inname);
    options.add(outname);
    options.add(tr);
    options.add(cardname);
    options.add(respname);
    options.add(cardorder);
    options.add(resporder);
    options.add(multc);
    options.add(multr);
    options.add(csfmask);
    options.add(rvt);
    options.add(heartrate);
    options.add(rvtsmooth);
    options.add(heartratesmooth);
    options.add(slicedir);
    options.add(sliceordering);
    options.add(slicetiming);
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
    
    // Call the local functions
    
    return do_work(argc,argv);
    
  } catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } catch(Exception &e) {
    cerr << e.what() << endl;
  } catch(...) {
    cerr << "Fatal Error" << endl;
  }
  return 0;
}

