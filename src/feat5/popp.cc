/*  popp.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2011 University of Oxford  */

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

// Pulse Oximeter Peak Picker  (also includes REPP - Respiratory Peak Picker)
//   adapted from matlab code originally written by MJ and modified by Jon Brooks and Yazhuo Kong

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

#include <vector>
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "utils/options.h"

using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;

// The two strings below specify the title and example usage that is
//  printed out as the help or usage message

string title="popp (Version 1.0)\nCopyright(c) 2011, University of Oxford (Mark Jenkinson)";
string examples="popp [options] -i <input data file> -o <output data file>";

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
		  string("output debugging information"),
		  false, no_argument);
Option<bool> noclean1(string("--noclean1"), false,
		  string("turn off cleanup phase 1"),
		  false, no_argument);
Option<bool> noclean2(string("--noclean2"), false,
		  string("turn off cleanup phase 2"),
		  false, no_argument);
Option<bool> invertcardiac(string("--invertcardiac"), false,
		  string("invert cardiac trace"),
		  false, no_argument);
Option<bool> invertresp(string("--invertresp"), false,
		  string("invert respiratory trace"),
		  false, no_argument);
Option<bool> rvt(string("--rvt"), false,
		  string("generate RVT data"),
		  false, no_argument);
Option<bool> heartrate(string("--heartrate"), false,
		  string("generate heartrate data"),
		  false, no_argument);
Option<bool> pulseox_trigger(string("--pulseox_trigger"), false,
		  string("specify that cardiac data is a trigger"),
		  false, no_argument);
Option<int> cardiac(string("--cardiac"), 0,
		    string("specify column number of cardiac input"),
		    false, requires_argument);
Option<int> resp(string("--resp"), 0,
		 string("specify column number of respiratory input"),
		 false, requires_argument);
Option<int> trigger(string("--trigger"), 0,
		 string("specify column number of trigger input"),
		 false, requires_argument);
Option<int> startingsample(string("--startingsample"), 0,
		 string("set sample number of the starting time (t=0)"),
		 false, requires_argument);
Option<float> initthreshc(string("--initthreshc"), 90,
		  string("initial threshold percentile for cardiac (default 90)"),
		  false, requires_argument);  
Option<float> neighbourthreshc(string("--nthreshc"), 60,
		  string("neighbourhood exclusion threshold percentile for cardiac (default 60)"),
		  false, requires_argument);
Option<float> initthreshr(string("--initthreshr"), 80,
		  string("initial threshold percentile for respiratory (default 80)"),
		  false, requires_argument);  
Option<float> neighbourthreshr(string("--nthreshr"), 80,
		  string("neighbourhood exclusion threshold percentile for respiratory (default 80)"),
		  false, requires_argument);
Option<float> smoothcard(string("--smoothcard"), 0.1, 
		  string("specify smoothing amount for cardiac (in seconds)"),
		  false, requires_argument);
Option<float> smoothresp(string("--smoothresp"), 0.1, 
		  string("specify smoothing amount for respiratory (in seconds)"),
		  false, requires_argument);
Option<float> highfreqcutoff(string("--highfreqcutoff"), 5.0, 
		  string("high frequency cut off for respiratory filter in Hz (default is 5Hz)"),
		  false, requires_argument);
Option<float> lowfreqcutoff(string("--lowfreqcutoff"), 0.1, 
		  string("low frequency cut off for respiratory filter in Hz (default is 0.1Hz)"),
		  false, requires_argument);
Option<float> samplingrate(string("-s,--samplingrate"), 100.0, 
		  string("sampling rate in Hz (default is 100Hz)"),
		  false, requires_argument);
Option<float> plotsampling(string("--plotsampling"), 40.0, 
		  string("Sampling rate in Hz (webpage plot)"),
		  false, requires_argument);
Option<float> tr(string("--tr"), 1.0, 
		  string("TR value in seconds"),
		  false, requires_argument);
Option<std::vector<float> > respadd(string("--respadd"), vector<float>(),
		  string("comma separated list of times (in sec) to add to respiratory peak list (uses nearest local max)"),
		  false, requires_argument);
Option<std::vector<float> > respdel(string("--respdel"), vector<float>(),
		  string("comma separated list of times (in sec) to delete from respiratory peak list (uses nearest existing peak)"),
		  false, requires_argument);
Option<std::vector<float> > cardadd(string("--cardadd"), vector<float>(),
		  string("comma separated list of times (in sec) to add to cardiac peak list (uses nearest local max)"),
		  false, requires_argument);
Option<std::vector<float> > carddel(string("--carddel"), vector<float>(),
		  string("comma separated list of times (in sec) to delete from cardiac peak list (uses nearest existing peak)"),
		  false, requires_argument);
Option<string> loadcardphase(string("--loadcardphase"), string(""),
		  string("input cardiac phase for reprocessing (text format)"),
		  false, requires_argument);
Option<string> loadrespphase(string("--loadrespphase"), string(""),
		  string("input respiratory phase for reprocessing (text format)"),
		  false, requires_argument);
Option<string> inname(string("-i,--in"), string(""),
		  string("input physiological data filename (text format)"),
		  true, requires_argument);
Option<string> volname(string("--vol"), string(""),
		  string("input volumetric image (EPI) filename"),
		  false, requires_argument);
Option<string> outname(string("-o,--out"), string(""),
		  string("output basename for physiological data and timing/triggers (no extensions)"),
		  true, requires_argument);
int nonoptarg;

// GLOBAL variables

ColumnVector times;
Matrix dataplot;

////////////////////////////////////////////////////////////////////////////

// Local functions (some of which should go into libraries!)

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

int mode(const vector<int>& vv, int binWidth)
{
  // TODO - Introduce a min/max limitation ?  or maybe a bin size?
  // mode function - for discrete data with not much spread
  //     (e.g. interval times in a long physio trace)
  // Uses an iterative approach of merging cells until max bin contains >10% of samples
  int maxvv=*max_element(vv.begin(),vv.end());
  vector<double> hist(maxvv+1,0); //The zero bin is not used but present to make indexing easier 
  for (int n=0; n<(int) vv.size(); n++) 
    hist[vv[n]]++;
  float totsamp=vv.size();
  int modev=0, maxh=0, niter=1;
  while ((binWidth<maxvv/2) && (maxh<0.1*totsamp)) {
    if (verbose.value()) { cout << "Mode loop: binsize = " << binWidth << " maxh/totsamp = " << maxh/totsamp << endl; }
    modev=0; maxh=0;
    for (int n=1; n<=maxvv; n+=binWidth) {
      float val(0);
      for (int m=n; m<MISCMATHS::Min(n+binWidth,maxvv); m++) { val+=hist[m]; }
      if (val>maxh) {
	maxh=val;
	modev=n + binWidth/2;
      }
    }
    // for next time around while loop (only needed if <10% samples in max bin)
    (niter++)==1 ? binWidth=2 : binWidth++;
  }
  return modev;
}


ColumnVector diff(const ColumnVector& data) {
  ColumnVector dv(data.Nrows()-1);
  for (int n=1; n<data.Nrows(); n++) {
    dv(n)=data(n+1)-data(n);
  }
  return dv;
}


vector<int> calc_diffs(const ColumnVector& data, float thresh) {
  vector<int> diffs;
  diffs.clear();
  int idx0=0;
  for (int n=1; n<=data.Nrows(); n++) {
    if (data(n)>thresh) {
      if ((idx0>0) && ((n-idx0)>1)) {
	diffs.push_back(n-idx0);
      }
      idx0=n;
    }
  }
  return diffs;
}

ColumnVector stl2newmat(const vector<int>& svec) {
  ColumnVector nvec(svec.size());
  for (int n=1; n<=(int) svec.size(); n++) { nvec(n) = svec[n-1]; }
  return nvec;
}

ColumnVector find(const ColumnVector& vec, float thresh=0.0f) {
  ColumnVector fvec(vec);
  int m=1;
  for (int n=1; n<=vec.Nrows(); n++) {
    if (vec(n)>thresh) { fvec(m++) = (float) n; }
  }
  fvec=fvec.Rows(1,m-1);
  return fvec;
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

int sort2(ColumnVector& vec, ColumnVector& idx) {
  // does a standard sort on the first argument and also returns the index of the sort
  idx=vec*0.0;  // initial size only
  // set up an stl sorted container - multimap (only using first element 'key' to sort)
  multimap<float, int> mm;
  int n=1;
  for (n=1; n<=vec.Nrows(); n++) {
    mm.insert(pair<float,int>(vec(n),n));
  }
  n=1;
  for (multimap<float, int>::iterator it = mm.begin(); it != mm.end() ; ++it) {
    vec(n)=(*it).first;
    idx(n)=(float) (*it).second;
    n++;
  }
  return 0;
}


// More specific support functions



ColumnVector get_times(const ColumnVector& tp) {
  ColumnVector timepts(tp.Nrows());
  for (int n=1; n<=tp.Nrows(); n++) { timepts(n) = times(MISCMATHS::round(tp(n))); }
  return timepts;
}

int nearest_time_index(float tval)
{
  // find closest index to requested time
  int bestidx=0;
  float tdiff=-1;
  for (int n=1; n<=times.Nrows(); n++) {
    if ((tdiff<0) || (fabs(times(n)-tval)<tdiff)) {
      tdiff=fabs(times(n)-tval);
      bestidx=n;
    }
  }
  return bestidx;
}

ColumnVector delete_times(const ColumnVector& peakind, const std::vector<float>& del_times)
{
  // peakind is an ordered list of indices for the peaks
  // for each delete time choose the closest in times(peakind(...)) and deletes it (can be multiple)
  if (verbose.value()) { cout << "Deleting points" << endl; }
  ColumnVector flags;
  flags=peakind;
  flags=1.0;
  float mindiff0=times(MISCMATHS::round(peakind.Maximum()));
  for (unsigned int n=1; n<=del_times.size(); n++) {
    int minloc=0;
    float mindiff=mindiff0;
    for (int m=1; m<=peakind.Nrows(); m++) {
      if (fabs(del_times[n-1]-times(MISCMATHS::round(peakind(m))))<mindiff) {
	mindiff=fabs(del_times[n-1]-times(MISCMATHS::round(peakind(m))));
	minloc=m;
      }
    }
    if (minloc>0) { flags(minloc)=0; }
    if (verbose.value()) { cout << "Deleted point " << times(MISCMATHS::round(peakind(minloc))) << " as closest to " << del_times[n-1] << endl; }
  }
  ColumnVector new_ind(MISCMATHS::round(flags.Sum()));
  int mm=1;
  for (int n=1; n<=peakind.Nrows(); n++) {
    if (flags(n)>0.5) {
      new_ind(mm++) = peakind(n);
    }
  }
  return new_ind;
}

int add1ind(ColumnVector& allind, double nind) {
  // only add indices (at end) that are not already present and different
  ColumnVector tt(1);
  if (nind>0.0) {
    tt=nind;
    if (allind.Nrows()==0) { allind = tt; }
    else if (fabs(allind(allind.Nrows())-nind)>0.5) { allind &= tt; }
  }
  return 0;
}

int localmaxidx(int startidx, const ColumnVector& vals, int dir=1) 
{
  int lmax;
  for (lmax=startidx; ((lmax>1) && (lmax<vals.Nrows())); lmax+=dir) {
    if ((vals(lmax-1)<=vals(lmax)) && (vals(lmax+1)<=vals(lmax))) { return lmax; }
  }
  return -1;
}

ColumnVector add_times(const ColumnVector& peakind, const std::vector<float>& add_times, 
		       const ColumnVector& vals)
{
  // peakind is an ordered list of indices for the peaks
  // choose nearest local max to add (irrespective of whether it is already an indicated point or not)
  if (verbose.value()) { cout << "Adding points" << endl; }
  ColumnVector new_ind(add_times.size());
  for (int n=1; n<=new_ind.Nrows(); n++) {
    int idx=nearest_time_index(add_times[n-1]);
    if (idx>0) {
      // find closest local max (forward and back) within vals
      int lmaxplus=localmaxidx(idx,vals,+1);
      int lmaxminus=localmaxidx(idx,vals,-1);
      int maxidx=lmaxplus;
      if ((abs(lmaxplus-idx)>abs(lmaxminus-idx)) || (lmaxplus<0)) maxidx=lmaxminus;
      if (maxidx>0) {
	new_ind(n)=maxidx;
      } else {
	new_ind(n)=-1.0;
      }
    }
  }
  SortAscending(new_ind);
  // insert values back into peakind (add at end and then sort)
  // peakind must be sorted upon entry
  ColumnVector allind;
  int n=1;
  for (int m=1; m<=peakind.Nrows(); m++) {
    while ((n<=new_ind.Nrows()) && (new_ind(n)<peakind(m))) { add1ind(allind,new_ind(n)); n++; }
    add1ind(allind,peakind(m));
  }
  for (;n<=new_ind.Nrows(); n++) { add1ind(allind,new_ind(n)); }  // do any stragglers
  return allind;
}


ColumnVector detect_triggers(const ColumnVector& datavals, float trigLevel) {
  ColumnVector trigs;
  trigs = datavals*0.0f;
  for (int ind=1; ind<=datavals.Nrows(); ind++) {
    if (datavals(ind)>trigLevel) {
      if ((ind==1) || (datavals(ind-1)<trigLevel)) {
	trigs(ind)=1.0;
      }
    }
  }
  return trigs;
}


// MAIN PEAK FINDING ALGORITHM (returns vector with non-peaks = 0 and peaks = 1)
//   minsep specifies a fixed minimum separation, meaninterval is an estimate only,
//   nthresh (between 0 and 100) is the percentile a peak must exceed its neighbours (in interval)
ColumnVector calc_peaks(const ColumnVector& datavals, int minsep, int meaninterval, float nthresh) {
  ColumnVector pp, datavalscopy;
  pp = datavals*0.0f;
  datavalscopy=datavals;
  int nn = datavals.Nrows();

  float maxv;
  int idx;
  while (datavalscopy.Maximum()>0) {
    maxv=datavalscopy.Maximum1(idx);
    pp(idx)=datavals(idx);
    for (int n=Max(1,idx-minsep); n<=Min(nn,idx+minsep); n++) {
      datavalscopy(n)=0;
    }
  }

  int stdivl=MISCMATHS::round(0.9*meaninterval);  // 90 percent of meaninterval (HBI)
 
  // Cleanup peaks
  ColumnVector data;
  for (int idx=1; idx<=pp.Nrows(); idx++) {
    if (pp(idx)>0) {
      // check that it is a local maxima
      if (datavals(Max(1,idx-1))>datavals(idx)) { pp(idx)=0; }
      if (datavals(Min(nn,idx+1))>datavals(idx)) { pp(idx)=0; }
      if (pp(idx)>0) {
        data=datavals.Rows(Max(1,idx-stdivl),Min(nn,idx+stdivl));
        SortAscending(data);
        // check that it is greater than the Nth percentile of the neighbourhood values
        if (datavals(idx)<data(MISCMATHS::round(data.Nrows()*nthresh/100))) { pp(idx)=0; }
      }
    }
  }

  for (int idx=1; idx<=pp.Nrows(); idx++) { if (pp(idx)>0) pp(idx)=1; }

  return pp;
}


int repp_cleanup(ColumnVector& pp, const ColumnVector& datavals , float thresh) {

  ColumnVector tp, dt, short_dt;

  // TODO - split this into cleanup1 and cleanup2 functions and allow cardiac to be processed too

  // PART 1:
  if (noclean1.value()) {
    // Delete peaks with very short interval <2s
    float short_int = 2;  // in seconds
    // MJ - TODO - it would be good to replace this hard threshold with one based on the median interval over a window (say 60 seconds) about the point of interest, and check for ones closer to half the this median interval
    int nshort=1;
    while (nshort>=1) {
      tp=find(pp,0.5);
      dt=diff(tp); 
      short_dt=find(dt*-1.0f,-1*short_int*samplingrate.value());  // checks dt<n  (i.e. -dt>-n)
      nshort=short_dt.Nrows();
      if (verbose.value()) { cout << "Found " << nshort << " short periods" << endl; }
      for (int ind=1; ind<=nshort; ind++) {
	// ORIGINALLY: "if pp(tp(short_dt)) < pp(tp(short_dt+1))"  // MJ: which would never work! 
	if (datavals(tp(short_dt(ind))) < datavals(tp(short_dt(ind)+1))) { 
	  pp(MISCMATHS::round(tp(short_dt(ind))))=0;
	} else {
	  pp(MISCMATHS::round(tp(short_dt(ind)+1)))=0;
	}
      }
    }
  }

  // PART 2:
  if (noclean2.value()) {
    // Delete artifacts (huge variance between peaks) and fake peaks
    // use dt_mean as interval and thresh1 as magnitude
    
    // MJ - TODO - again this would be better to do using a local threshold based on a median interval within a window (say 60 seconds)

    float dt_mean=mean(dt).AsScalar();
    ColumnVector dt_std, data;
    for (int m=1; m<tp.Nrows(); m++) {
      data=datavals.Rows(MISCMATHS::round(tp(m)),MISCMATHS::round(tp(m+1)));
      dt_std = dt_std & stdev(data);  // MJ: this is the only output of this loop: stddev in each interval
    }
  
    // Large variance greater than the 90th percentile of all peak to peak variance
    float thresh2=percentile(dt_std,99); 
    ColumnVector big_std=find(dt_std,thresh2);  // guaranteed 1% found here  (returns interval numbers)
    ColumnVector bigstd_diff = diff(big_std);  // how close they are in intervals
  
    int ind = 1; 
    ColumnVector bigstd_tag1; // big variance starts
    ColumnVector bigstd_tag2; // big variance ends
    while (ind <= bigstd_diff.Nrows()) {
      if (bigstd_diff(ind) > 3) {  // Only do <= 3 intervals separation (large stddev close together!)
	// do nothing
	ind++;
      } else {
	ColumnVector tmp(1);
	tmp = (float) ind;
	bigstd_tag1 = bigstd_tag1 & tmp;  // accumulate interval nums into bigstd_diff
	ind++;
	if (ind <= bigstd_diff.Nrows()) {
	  while ((ind<=bigstd_diff.Nrows()) && (bigstd_diff(ind) <= 3)) {
	    ind++;
	  }
	}
	tmp = (float) (ind-1);
	bigstd_tag2 &= tmp;
      }
    }
  
    // The result of the above is to get the start and end interval numbers of sections of the data that represent short patches of high variance (highest 1% of variance in intervals and intervals not separated by more than 2 intervals from each other, although these can be strung together)

    for (ind=1; ind<=bigstd_tag1.Nrows(); ind++) {
      if (bigstd_tag1(ind) != bigstd_tag2(ind)) {
	int tp1=MISCMATHS::round(big_std(bigstd_tag1(ind)));  // interval number of start of sequence
	int tp2=MISCMATHS::round(big_std(bigstd_tag2(ind)+1)+1); // interval number of end of sequence (plus 1)
	for (int n=MISCMATHS::round(tp(tp1)+1); n<=MISCMATHS::round(tp(tp2)-1); n++) {
	  // tp gives timepoint for start of each interval, which means that this 
	  // deletes all peaks within this noisy region (but not the end ones)
	  if (verbose.value()) { cout << "Deleting peak at timepoint " << n << endl; }
	  pp(n)=0; 
	}
	// put in fake peaks at an average separation instead!
	int peak_number=MISCMATHS::round((tp(tp2)-tp(tp1))/dt_mean);
	float dt_new=MISCMATHS::round((tp(tp2)-tp(tp1))/peak_number);
	for (int n=1; n<peak_number; n++) {
	  if (verbose.value()) { cout << "Adding peak at timepoint " << MISCMATHS::round(tp(tp1)+dt_new*n) << endl; }
	  pp(MISCMATHS::round(tp(tp1)+dt_new*n))=thresh; 
	}
      }
    }
  }
  return 0;
}

bool check_valid_column(int cn, const Matrix& rpm) {
  if ((cn<1) || (cn>rpm.Ncols())) {
    cerr << "Cannot access column number " << cn << " in matrix " << inname.value() << endl;
    exit(1);
  }
  return true;
}

Matrix calc_rvt(const ColumnVector& tp, const ColumnVector& datavals) {
  // input tp is a list of the index values (not in seconds)
  Matrix rvt(tp.Nrows()-1,2);
  rvt = 0.0f;
  if (verbose.value()) { cout << "Pre rvt loop" << endl; }
  for (int n=1; n<=tp.Nrows()-1; n++) {
    float tav=(times(MISCMATHS::round(tp(n)))+times(MISCMATHS::round(tp(n+1))))/2.0;
    rvt(n,1)=tav;
    float tdiff=tp(n+1)-tp(n);
    float rmax=MISCMATHS::Max(datavals(MISCMATHS::round(tp(n))),datavals(MISCMATHS::round(tp(n+1))));
    ColumnVector subv;
    float rmin;
    if (tdiff>0.0) { 
      subv=datavals.Rows(MISCMATHS::round(tp(n))+1,MISCMATHS::round(tp(n+1))-1); 
      rmin=subv.Minimum(); 
    } else { 
      rmin=0.0; 
    }
    rvt(n,2)=(rmax-rmin)/tdiff;
  }
  if (verbose.value()) { cout << "Post rvt loop" << endl; }
  return rvt;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
int set_dataplot(const ColumnVector& datavals, const ColumnVector& tp, int colno)
{
      // write out text file for javascript/browser editor
      // calculate timepoints for 40Hz sampling
      int npts=dataplot.Nrows();
      float datavalsmax, datavalsmin;
      datavalsmax=datavals.Maximum();
      datavalsmin=datavals.Minimum();
      for (int nn=1; nn<=npts; nn++) {
	float tval=(nn-1)/plotsampling.value();
	// interpolate respiratory trace (datavals) at these times
	int datavalsidx=MISCMATHS::round(tval*samplingrate.value());
	if (datavalsidx<=0) datavalsidx=1;
	if (datavalsidx>datavals.Nrows()) datavalsidx=datavals.Nrows();
	dataplot(nn,colno)=(datavals(datavalsidx)-datavalsmin)/(datavalsmax-datavalsmin);
      }
      // create a peak trace and fill the appropriate timepoints with appropriate values
      for (int nn=1; nn<=tp.Nrows(); nn++) {
	int tidx=MISCMATHS::round(tp(nn)/samplingrate.value()*plotsampling.value());
	if (tidx>0) dataplot(tidx,colno+1)=1;
      }
      return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////

int process_scanner_triggers(const ColumnVector& datavals) {
  // TRIGGERS
  // detect slice/volume triggers and mark initial detection in Etrig

  // loop is to progressively lower the threshold in case the initial one is too high
  ColumnVector Etrig(1);
  Etrig=0.0;
  float frac=1.0, trigLevel;
  while (Etrig.Sum()<=0.0) {
    frac-=0.1;
    if (frac<=0.0) { imthrow("Cannot find any valid triggers",60);  }
    // threshold starts at 90% of the 6th largest value (assume up to 5 outliers!)  
    trigLevel = frac*percentile(datavals,100.0*(1.0-5.0/datavals.Nrows())) + (1-frac)*percentile(datavals,5.0);
    if (debug.value()) { cout << "trigLevel is " << trigLevel << endl; }
    Etrig = detect_triggers(datavals,trigLevel);
  }

  ColumnVector trigpts;
  trigpts=find(Etrig,0.5);
  int ntrigs=trigpts.Nrows();
  float sec_per_trig=0.0;
  //   if ntrigs = number of time points then sec_per_trig=TR
  //   else if ntrigs = number of time points * number of slices, sec_per_trig=TR/nslices
  //   estimate how many slices per trigger based on TR / average time between trigger points
  float nslice_est = tr.value() / ((trigpts(ntrigs)-trigpts(1))/(samplingrate.value()*(ntrigs-1)));
  if ((nslice_est>0.9) && (nslice_est<1.1)) { sec_per_trig = tr.value(); }
  else if ((nslice_est>1.0) && (fabs(nslice_est-(MISCMATHS::round(nslice_est)))<0.2)) { sec_per_trig = tr.value()/MISCMATHS::round(nslice_est); }
  else { cerr << "ERROR:: Problem with trigger timing" << endl << "Estimated time per trigger does not match with integer number of slices" << endl << "Time per trigger / TR = " << 1/nslice_est << endl; }

  ColumnVector timepts(ntrigs);
  for (int n=1; n<=ntrigs; n++) { timepts(n) = (n-1)*sec_per_trig; }
  // set the GLOBAL variable "times" (to be the acquisition time of each timept)
  times.ReSize(datavals.Nrows());
  dataplot.ReSize(datavals.Nrows(),5);
  for (int n=1; n<=datavals.Nrows(); n++) { 
    if ( (n>=trigpts(1)) && (n<=trigpts(ntrigs)) ) {
      times(n) = interp1(trigpts,timepts,n);  
    } else {
      // need to extrapolate at either end - use sampling rate for this
      if (n<trigpts(1)) {
	times(n) = (n-trigpts(1))*(1.0/samplingrate.value()) + timepts(1);
      } else {
	times(n) = (n-trigpts(ntrigs))*(1.0/samplingrate.value()) + timepts(ntrigs);
      }
    }
    // also save into dataplot Matrix for webpage report
    dataplot(n,1)=times(n);
  }
  
  write_ascii_matrix(times,outname.value()+"_time.txt");
  return 0;
}

int set_default_sample_times(int ntimepts) {
  // WARNING MESSAGE ?
  times.ReSize(ntimepts);
  for (int n=1; n<=ntimepts; n++) {
    times(n) = (n-startingsample.value())*(1.0/samplingrate.value());
  }
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////

int process_cardiac(const ColumnVector& datavals_orig) {

  ColumnVector peakind, peakmask, timepts, datavals;
  datavals=datavals_orig;
  if (loadcardphase.unset()) {
    // PULSE-OX TRIGGER FILE 
    if (pulseox_trigger.value()) {
      // calculate the pulse ox triggers here (leave out the top 10, as potential outliers)
      float PtrigLevel = 0.8*percentile(datavals,100.0*(1.0-10.0/datavals.Nrows())) + 0.2*percentile(datavals,50);
      peakmask = detect_triggers(datavals,PtrigLevel);   // binary vector same length as datavals
      peakind = find(peakmask,0.5);
    } else {
      // Optional Smoothing - (originally done in popp)
      if (smoothcard.value()>0) {
	// moving average filter
	moving_average(datavals,MISCMATHS::round(smoothcard.value() * samplingrate.value()/2)*2 + 1);
      }
      if (debug.value()) { write_ascii_matrix(datavals,"debug_smooth_card"); }
      
      // initialise sensible HBI (heart beat interval)
      // loop over progressively decreasing percentage thresholds, in case of clipping near top
      vector<int> diffs;
      int hbi, minsep;
      float thresh, pthresh=initthreshc.value();
      while (diffs.size()<=0) {
	// find percentile for threshold
	thresh = percentile(datavals,pthresh); // empirical val (not too sensitive to this)
	if (verbose.value()) { cout << "Cardiac threshold = " << thresh << endl; }
	// Do the peak calculation
	diffs = calc_diffs(datavals,thresh);
	if (pthresh<10) { imthrow("Cannot find any valid cardiac peaks",60);  }
	pthresh-=10;
      }


      if (debug.value()) { ColumnVector tmp(diffs.size()); for (int n=1; n<=tmp.Nrows(); n++) {tmp(n)=diffs[n-1];} write_ascii_matrix(tmp,"diffs4mode"); }
      hbi=mode(diffs,MISCMATHS::round(0.1*samplingrate.value()));   // 0.1 seconds as a rough bin size guess
      if (verbose.value()) { cout << "Mean Heart Beat Interval = " << hbi/samplingrate.value() << endl; }
      minsep = MISCMATHS::round(hbi/2);
      peakmask = calc_peaks(datavals,minsep,hbi,neighbourthreshc.value());   // binary vector same length as datavals
      
      if (debug.value()) { 
	ColumnVector diffvec; 
	diffvec = diff(find(peakmask,0.5)); 
	write_ascii_matrix(diffvec,"card_peak_diffs"); 
      }
    }
    
    peakind=find(peakmask,0.5);

    // delete and add requested points
    if (carddel.set()) { peakind=delete_times(peakind,carddel.value()); }
    if (cardadd.set()) { peakind=add_times(peakind,cardadd.value(),datavals); }

    timepts = get_times(peakind);
    set_dataplot(datavals,peakind,2);
    write_ascii_matrix(timepts,outname.value()+"_card.txt");
  } else {
    // Load external (e.g. manually fixed) data
    timepts = read_ascii_matrix(loadcardphase.value());
  }

  // HEARTRATE
  if (heartrate.value()) {
    Matrix hr(timepts.Nrows()-1,2);
    for (int n=1; n<timepts.Nrows(); n++) { 
      float tav=(timepts(n)+timepts(n+1))/2.0;
      float bpm=60.0/(timepts(n+1)-timepts(n));
      hr(n,1) = tav;
      hr(n,2) = bpm;
    }
    write_ascii_matrix(hr,outname.value()+"_hr.txt");
  }

  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////

int process_respiratory(const ColumnVector &datavals_orig) {
  
  ColumnVector datavals, peakmask, peakind, timepts;
  datavals=datavals_orig;

  if (loadrespphase.unset()) {
    
  // band pass filter (0.1<f<2Hz) the respiratory data to account for slow
  // drift and noise (air leaks?) from the respiratory bellows
  
  // Use 5Hz rather than 2Hz for lowpass filter - more conservative (moving_av vs butter)
  // Note:  round(blah/2)*2 + 1  guarantees the closest (ceil) odd number
  // 0.1Hz to 5Hz is the default

    int n_lp_filt = MISCMATHS::round(1/highfreqcutoff.value() * samplingrate.value()/2)*2 + 1;
    moving_average(datavals,n_lp_filt);
    int n_hp_filt = MISCMATHS::round(1/lowfreqcutoff.value() * samplingrate.value()/2)*2 + 1;
    ColumnVector datavalslow;
    datavalslow = datavals;
    moving_average(datavalslow,n_hp_filt);
    datavals -= datavalslow;
    if (debug.value()) { write_ascii_matrix(datavals,"datavals_pre_norm"); }

    // Optional Smoothing - (originally done in popp)
    if (smoothresp.value()>0) {
      if (verbose.value()) { cout << "Smoothing with " << MISCMATHS::round(smoothresp.value() * samplingrate.value()/2)*2 + 1 << endl; }
      // moving average filter
      moving_average(datavals,MISCMATHS::round(smoothresp.value() * samplingrate.value()/2)*2 + 1);
    }
    if (debug.value()) { write_ascii_matrix(datavals,"debug_smooth_resp"); }
    
    ColumnVector datavalsnorm;
    // other output is histogram normalised respiratory data (uses sort)
    {
      ColumnVector idx, datavalscopy;
      datavalscopy = datavals;
      sort2(datavals,idx);
      datavals = datavalscopy;
      sort2(idx,datavalsnorm);
      datavalsnorm /= datavals.Nrows();
    }
    
    if (debug.value()) { write_ascii_matrix(datavalsnorm,"datavals_post_norm"); }
    // Convert to phase (use gradient sign to distinguish inhale from exhale)
    ColumnVector grad(datavalsnorm.Nrows());
    for (int n=2; n<grad.Nrows(); n++) { grad(n)=datavalsnorm(n+1)-datavalsnorm(n-1); }
    grad(1)=grad(2);
    grad(grad.Nrows())=grad(grad.Nrows()-1);
    moving_average(grad,20);
    
    Matrix respmat(datavals.Nrows(),2);
    for (int n=1; n<=datavals.Nrows(); n++) {
      respmat(n,1) = times(n);
      int factor=0;  // if grad(n)==0 
      if (grad(n)>0) factor=1;
      if (grad(n)<0) factor=-1;
      respmat(n,2) = M_PI * factor * datavalsnorm(n);
    }
    write_ascii_matrix(respmat,outname.value()+"_resp.txt");

    if (rvt.set()) { 
      // used to pass non-normalised data to repp(), but I think it makes little
      //  difference, so might as well pass the normalised one
            
      // initialise sensible BBI (breath-to-breath interval)
      // loop over progressively decreasing percentage thresholds, in case of clipping near top
      vector<int> diffs;
      int bbi, minsep;
      float thresh, pthresh=initthreshr.value();
      while (diffs.size()<=0) {
	// find percentile for threshold
	float thresh = percentile(datavals,pthresh); // empirical val (not too sensitive to this)
	if (verbose.value()) { cout << "Threshold = " << thresh << endl; }
	if (verbose.value()) { cout << "Min / Max = " << datavals.Minimum() << " & " << datavals.Maximum() << endl; }
	// Do the peak calculation
	if (verbose.value()) { cout << "Calculating differences" << endl; }
	if (debug.value()) { write_ascii_matrix(datavals,"datavals_pre_diffs"); }
	diffs = calc_diffs(datavals,thresh);
	if (pthresh<10) { imthrow("Cannot find any valid cardiac peaks",60);  }
	pthresh-=10;
      }

      if (verbose.value()) { cout << "Calculating mode: number of diffs = " << diffs.size() << endl; }
      bbi=mode(diffs,MISCMATHS::round(0.2*samplingrate.value()));  // 0.2 seconds as a rough bin size
      if (verbose.value()) { cout << "Mean Breath to Breath Interval = " << bbi/samplingrate.value() << endl; }
      minsep = MISCMATHS::round(bbi/2);
      peakmask = calc_peaks(datavals,minsep,bbi,neighbourthreshr.value());
      
      if (debug.value()) { 
	ColumnVector diffvec; 
	diffvec = diff(find(peakmask,0.5)); 
	write_ascii_matrix(diffvec,"resp_peak_diffs"); 
      }
      
      repp_cleanup(peakmask,datavals,thresh);

      peakind=find(peakmask,0.5);
      
      // delete and add requested points
      if (respdel.set()) { peakind=delete_times(peakind,respdel.value()); }
      if (respadd.set()) { peakind=add_times(peakind,respadd.value(),datavals); }

      timepts=get_times(peakind);
      set_dataplot(datavals,peakind,4);
    } 
  } else {
    // Load external (e.g. manually fixed) data
    timepts = read_ascii_matrix(loadrespphase.value());
  }
  
  if (rvt.set()) {
    // feed calc_rvt the peak times and raw respiratory trace
    Matrix rvtmat;
    if (verbose.value()) { cout << "Pre RVT calculation" << endl; }
    rvtmat = calc_rvt(peakind,datavals);
    write_ascii_matrix(rvtmat,outname.value()+"_rvt.txt");
  }

  return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// Main function that does the work
int do_work(int argc, char* argv[]) 
{
  Matrix rpm;
  rpm = read_ascii_matrix(inname.value());
  if (verbose.value()) { cout << "Input matrix is " << rpm.Nrows() << " by " << rpm.Ncols() << endl; }
  if (verbose.value() && respdel.set()) { cout << "Respiratory delete times are: " << respdel.value()[0] << " to " << respdel.value()[respdel.value().size()-1] << endl; }
  int nn=rpm.Nrows();
  // make sure minimum in each column is set to 0 (and hence all values are positive)
  for (int cn=1; cn<=rpm.Ncols(); cn++) {
    rpm.Column(cn)=rpm.Column(cn)-rpm.Column(cn).Minimum();
  }
  
  if (verbose.value()) { cout << "Checking options" << endl; }
  // Sanity checking and pre-processing (inverting)
  if (trigger.set()) { check_valid_column(trigger.value(),rpm); }
  if (resp.set()) {
    check_valid_column(resp.value(),rpm);
    if (invertresp.value()) { rpm.Column(resp.value()) = -rpm.Column(resp.value()); } 
  }
  if (cardiac.set()) {
    check_valid_column(cardiac.value(),rpm);
    if (invertcardiac.value()) { rpm.Column(cardiac.value()) = -rpm.Column(cardiac.value()); } 
  }

  if (verbose.value()) { cout << "Main Processing" << endl; }
  // Do the main processing
  ColumnVector datavals;
  if (trigger.set()) {
    if (verbose.value()) { cout << "Processing triggers in column " << trigger.value() << endl; }
    datavals = rpm.Column(trigger.value());
    process_scanner_triggers(datavals);
  } else {
    set_default_sample_times(nn);
  }

  // set webpage plotting (but in a global, so modified in other places too)
  dataplot.ReSize(MISCMATHS::round(times.Nrows()*plotsampling.value()/samplingrate.value()),5);
  dataplot=0.0;
  for (int nn=1; nn<=dataplot.Nrows(); nn++) { 
    int tidx=nn*samplingrate.value()/plotsampling.value();
    dataplot(nn,1)=times(Max(Min(tidx,times.Nrows()),1));
  }

  if (resp.set()) {
    if (verbose.value()) { cout << "Processing resp in column " << resp.value() << endl; }
    datavals = rpm.Column(resp.value());
    process_respiratory(datavals);
 }
  if (cardiac.set()) {
    if (verbose.value()) { cout << "Processing cardiac in column " << cardiac.value() << endl; }
    datavals = rpm.Column(cardiac.value());
    process_cardiac(datavals);
  }

  // save webpage plotting data
  //write_ascii_matrix(dataplot,outname.value()+"_web.txt");
  ofstream output_file((outname.value()+"_pnm.js").c_str());
  output_file.setf(ios::floatfield);  // use fixed or scientific notation as appropriate
  output_file << "function pnm_data() {" << endl;
  output_file << "return \"\" +" << endl;
  for (int nn=1; nn<=dataplot.Nrows(); nn++) {
    if (nn>1) { output_file << " + " << endl; }
    output_file << "\"" << dataplot(nn,1) << "," << dataplot(nn,2) << "," << dataplot(nn,3) 
		<< "," << dataplot(nn,4) << "," << dataplot(nn,5) << "\\n\"";
  }
  output_file << ";" << endl << "}" << endl;
  output_file.close();
  return 0;
}

////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{

  Tracer trcr("main");
  OptionParser options(title, examples);

  try {
    // must include all wanted options here (the order determines how
    //  the help message is printed)
    options.add(inname);
    options.add(outname);
    options.add(samplingrate);
    options.add(tr);
    options.add(resp);
    options.add(cardiac);
    options.add(trigger);
    options.add(rvt);
    options.add(heartrate);
    options.add(pulseox_trigger);
    options.add(smoothcard);
    options.add(smoothresp);
    options.add(highfreqcutoff);
    options.add(lowfreqcutoff);
    options.add(initthreshc);
    options.add(neighbourthreshc);
    options.add(initthreshr);
    options.add(neighbourthreshr);
    options.add(invertresp);
    options.add(invertcardiac);
    options.add(noclean1);
    options.add(noclean2);
    options.add(loadcardphase);
    options.add(loadrespphase);
    options.add(volname);
    options.add(startingsample);
    options.add(respadd);
    options.add(respdel);
    options.add(cardadd);
    options.add(carddel);
    options.add(verbose);
    options.add(debug);
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

  }  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } catch(Exception &e) {
    cerr << e.what() << endl;
  } catch (...) {
    cerr << "Fatal error" << endl;
  }
}

