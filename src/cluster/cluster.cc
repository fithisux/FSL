/*  cluster.cc

    Mark Jenkinson & Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 2000-2008 University of Oxford  */

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

#ifndef EXPOSE_TREACHEROUS
#define EXPOSE_TREACHEROUS
#endif    

#include <vector>
#include <algorithm>
#include <iomanip>
#include "newimage/fmribmain.h"
#include "newimage/newimageall.h"
#include "utils/options.h"
#include "infer.h"
#include "warpfns/warpfns.h"
#include "warpfns/fnirt_file_reader.h"

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1        

using namespace NEWIMAGE;
using std::vector;
using namespace Utilities;
using namespace NEWMAT;
using namespace MISCMATHS;

string title="cluster \nCopyright(c) 2000-2013, University of Oxford (Mark Jenkinson, Matthew Webster)";
string examples="cluster --in=<filename> --thresh=<value> [options]";

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<bool> mm(string("--mm"), false,
		  string("use mm, not voxel, coordinates"),
		  false, no_argument);
Option<bool> minv(string("--min"), false,
		  string("find minima instead of maxima"),
		  false, no_argument);
Option<bool> fractional(string("--fractional"), false,
			string("interprets the threshold as a fraction of the robust range"),
			false, no_argument);
Option<bool> no_table(string("--no_table"), false,
		      string("suppresses printing of the table info"),
		      false, no_argument);
Option<bool> minclustersize(string("--minclustersize"), false,
		      string("prints out minimum significant cluster size"),
		      false, no_argument);
Option<int> sizethreshold(string("--minextent"), 1,
		      string("do not report clusters with less than extent voxels"),
		      false, requires_argument);
Option<int> voxvol(string("--volume"), 0,
		   string("number of voxels in the mask"),
		   false, requires_argument);
Option<int> numconnected(string("--connectivity"), 26,
		   string("the connectivity of voxels (default 26)"),
		   false, requires_argument);
Option<int> mx_cnt(string("-n,--num"), 6,
		   string("no of local maxima to report"),
		   false, requires_argument);
Option<float> dLh(string("-d,--dlh"), 4.0,
		  string("smoothness estimate = sqrt(det(Lambda))"),
		  false, requires_argument);
Option<float> thresh(string("-t,--thresh,--zthresh"), 2.3,
		     string("threshold for input volume"),
		     true, requires_argument);
Option<float> pthresh(string("-p,--pthresh"), 0.01,
		      string("p-threshold for clusters"),
		      false, requires_argument);
Option<float> peakdist(string("--peakdist"), 0,
		      string("minimum distance between local maxima/minima, in mm (default 0)"),
		      false, requires_argument);
Option<string> inputname(string("-i,--in,-z,--zstat"), string(""),
			 string("filename of input volume"),
			 true, requires_argument);
Option<string> copename(string("-c,--cope"), string(""),
			string("filename of input cope volume"),
			false, requires_argument);
Option<string> outpvals(string("--opvals"), string(""),
			string("filename for image output of log pvals"),
			false, requires_argument);
Option<string> outindex(string("-o,--oindex"), string(""),
			string("filename for output of cluster index (in size order)"),
			false, requires_argument);
Option<string> outlmax(string("--olmax"), string(""),
			string("filename for output of local maxima text file"),
			false, requires_argument);
Option<string> outlmaxim(string("--olmaxim"), string(""),
			string("filename for output of local maxima volume"),
			false, requires_argument);
Option<string> outthresh(string("--othresh"), string(""),
			 string("filename for output of thresholded image"),
			 false, requires_argument);
Option<string> outsize(string("--osize"), string(""),
		       string("filename for output of size image"),
		       false, requires_argument);
Option<string> outmax(string("--omax"), string(""),
		      string("filename for output of max image"),
		      false, requires_argument);
Option<string> outmean(string("--omean"), string(""),
		       string("filename for output of mean image"),
		       false, requires_argument);
Option<string> transformname(string("-x,--xfm"), string(""),
		       string("filename for Linear: input->standard-space transform. Non-linear: input->highres transform"),
		       false, requires_argument);
Option<string> stdvolname(string("--stdvol"), string(""),
		       string("filename for standard-space volume"),
		       false, requires_argument);
Option<string> warpname(string("--warpvol"), string(""),
		       string("filename for warpfield"),
		       false, requires_argument);
Option<string> scalarname(string("--scalarname"), string(""),
		       string("give name of scalars (e.g. Z) - to be used in printing output tables"),
		       false, requires_argument);


int num(const char x) { return (int) x; }
short int num(const short int x) { return x; }
int num(const int x) { return x; }
float num(const float x) { return x; }
double num(const double x) { return x; }


template <class T>
struct triple {
T x,y,z;
triple() {}
triple(const T ix,const T iy,const T iz) : x(ix),y(iy),z(iz) {}
};

template <class T>
struct cluster {
cluster();
int originalLabel;
int size;
T maxval;
float meanval;
triple<float> maxpos; //float as may be mm or vox
triple<float> cog;
float pval;
float logpval;
};

template <class T>
cluster<T>::cluster() : originalLabel(0), size(0), maxval(0), meanval(0), pval(1),logpval(0) {
  maxpos.x=maxpos.y=maxpos.z=cog.x=cog.y=cog.z=0;
}

template <class T>
bool operator< (const cluster<T> &c1, const cluster<T> &c2)
{
    return c1.size < c2.size;
}

template <class T>
bool operator< (const pair<T, triple<float> > &p1, const pair<T, triple<float> > &p2)
{
    return p1.first < p2.first;
}

template <class T, class S>
void copyconvert(const vector<triple<T> >& oldcoords, 
		 vector<triple<S> >& newcoords)
{
  newcoords.erase(newcoords.begin(),newcoords.end());
  newcoords.resize(oldcoords.size());
  for (unsigned int n=0; n<oldcoords.size(); n++) {
    newcoords[n].x = (S) oldcoords[n].x;
    newcoords[n].y = (S) oldcoords[n].y;
    newcoords[n].z = (S) oldcoords[n].z;
  }
}

template <class T>
void MultiplyCoordinateVector(triple<T> & coords, const Matrix& mat)
{
  ColumnVector vec(4);
  vec << coords.x << coords.y << coords.z << 1.0;
  vec = mat * vec;     // apply voxel xfm
  coords.x = vec(1);
  coords.y = vec(2);
  coords.z = vec(3);
}

template <class T, class S>
void TransformToReference(triple<T> & coordlist, const Matrix& affine, 
			  const volume<S>& source, const volume<S>& dest, const volume4D<float>& warp,bool doAffineTransform, bool doWarpfieldTransform)
{
  ColumnVector coord(4);
  coord << coordlist.x << coordlist.y << coordlist.z << 1.0;
  if ( doAffineTransform && doWarpfieldTransform ) coord = NewimageCoord2NewimageCoord(affine,warp,true,source,dest,coord);
  if ( doAffineTransform && !doWarpfieldTransform) coord = NewimageCoord2NewimageCoord(affine,source,dest,coord);
  if ( !doAffineTransform && doWarpfieldTransform) coord = NewimageCoord2NewimageCoord(warp,true,source,dest,coord);
  coordlist.x = coord(1);
  coordlist.y = coord(2);
  coordlist.z = coord(3);
}

template <class T>
bool checkIfLocalMaxima(const int& index, const volume<int>& labelim, const volume<T>& zvol, const int& connectivity, const int& x, const int& y, const int& z )
{	       
  if (connectivity==6)
    return ( index==labelim(x,y,z) &&
	     zvol(x,y,z)>zvol(x,  y,  z-1) &&
	     zvol(x,y,z)>zvol(x,  y-1,z) &&
	     zvol(x,y,z)>zvol(x-1,y,  z) &&
	     zvol(x,y,z)>=zvol(x+1,y,  z) &&
	     zvol(x,y,z)>=zvol(x,  y+1,z) &&
	     zvol(x,y,z)>=zvol(x,  y,  z+1) );

  else 
    return ( index==labelim(x,y,z) &&
	     zvol(x,y,z)>zvol(x-1,y-1,z-1) &&
	     zvol(x,y,z)>zvol(x,  y-1,z-1) &&
	     zvol(x,y,z)>zvol(x+1,y-1,z-1) &&
	     zvol(x,y,z)>zvol(x-1,y,  z-1) &&
	     zvol(x,y,z)>zvol(x,  y,  z-1) &&
	     zvol(x,y,z)>zvol(x+1,y,  z-1) &&
	     zvol(x,y,z)>zvol(x-1,y+1,z-1) &&
	     zvol(x,y,z)>zvol(x,  y+1,z-1) &&
	     zvol(x,y,z)>zvol(x+1,y+1,z-1) &&
	     zvol(x,y,z)>zvol(x-1,y-1,z) &&
	     zvol(x,y,z)>zvol(x,  y-1,z) &&
	     zvol(x,y,z)>zvol(x+1,y-1,z) &&
	     zvol(x,y,z)>zvol(x-1,y,  z) &&
	     zvol(x,y,z)>=zvol(x+1,y,  z) &&
	     zvol(x,y,z)>=zvol(x-1,y+1,z) &&
	     zvol(x,y,z)>=zvol(x,  y+1,z) &&
	     zvol(x,y,z)>=zvol(x+1,y+1,z) &&
	     zvol(x,y,z)>=zvol(x-1,y-1,z+1) &&
	     zvol(x,y,z)>=zvol(x,  y-1,z+1) &&
	     zvol(x,y,z)>=zvol(x+1,y-1,z+1) &&
	     zvol(x,y,z)>=zvol(x-1,y,  z+1) &&
	     zvol(x,y,z)>=zvol(x,  y,  z+1) &&
	     zvol(x,y,z)>=zvol(x+1,y,  z+1) &&
	     zvol(x,y,z)>=zvol(x-1,y+1,z+1) &&
	     zvol(x,y,z)>=zvol(x,  y+1,z+1) &&
	     zvol(x,y,z)>=zvol(x+1,y+1,z+1) );

}

template <class T>
void get_stats(const volume<int>& labelim, const volume<T>& origim,vector<cluster<T> >& clusters,bool minv) 
{
  int labelnum = labelim.max();
  clusters.resize(labelnum);
  for (int z=labelim.minz(); z<=labelim.maxz(); z++) {
    for (int y=labelim.miny(); y<=labelim.maxy(); y++) {
      for (int x=labelim.minx(); x<=labelim.maxx(); x++) {
	int idx = labelim(x,y,z);
	if ( idx-- ) {
	  clusters[idx].originalLabel=idx+1; //slightly wasteful, but doesn't really matter
	  T oxyz = origim(x,y,z);
	  clusters[idx].size++;
	  clusters[idx].cog.x+=((float) oxyz)*x;
	  clusters[idx].cog.y+=((float) oxyz)*y;
	  clusters[idx].cog.z+=((float) oxyz)*z;
	  clusters[idx].meanval+=(float) oxyz;
	  if ((clusters[idx].size==1) || ((oxyz>clusters[idx].maxval) && !minv ) || ((oxyz<clusters[idx].maxval) && minv )) {
	    clusters[idx].maxval = oxyz;
	    clusters[idx].maxpos.x = x;
	    clusters[idx].maxpos.y = y;
	    clusters[idx].maxpos.z = z;
	  }
	}
      }
    }
  }
  for (unsigned int n=0; n<clusters.size(); n++) {
    if (clusters[n].size) {
      clusters[n].cog.x /= clusters[n].meanval; //meanval is currently just sum
      clusters[n].cog.y /= clusters[n].meanval;
      clusters[n].cog.z /= clusters[n].meanval;
      clusters[n].meanval /= clusters[n].size;
    }
  }
}

template <class T, class S>
void relabel_image(const volume<int>& labelim, volume<T>& relabelim,
		   const vector<S>& newlabels)
{
  copyconvert(labelim,relabelim);
  for (int z=relabelim.minz(); z<=relabelim.maxz(); z++) 
    for (int y=relabelim.miny(); y<=relabelim.maxy(); y++) 
      for (int x=relabelim.minx(); x<=relabelim.maxx(); x++) 
	relabelim(x,y,z) = (T) newlabels[labelim(x,y,z)];
}

template <class T>
void print_results(vector<cluster<T> >& clusters, 
		   vector<cluster<T> >& clustersCope,
		   const volume<T>& zvol, const volume<T>& cope, 
		   const volume<int> &labelim)
{
  bool doAffineTransform=false;
  bool doWarpfieldTransform=false;
  volume4D<float>  full_field;
  volume<T> stdvol;
  Matrix trans;
  const volume<T> *refvol = &zvol;
  if ( transformname.set() && stdvolname.set() ) {
    read_volume(stdvol,stdvolname.value());
    trans = read_ascii_matrix(transformname.value());
    if (verbose.value()) { 
      cout << "Transformation Matrix filename = "<<transformname.value()<<endl;
      cout << trans.Nrows() << " " << trans.Ncols() << endl; 
      cout << "Transformation Matrix = " << endl;
      cout << trans << endl;
    }
    doAffineTransform=true;
  }

  FnirtFileReader   reader;
  if (warpname.value().size()) 
  {
    reader.Read(warpname.value());
    full_field = reader.FieldAsNewimageVolume4D(true);
    doWarpfieldTransform=true;
  }

  if ( doAffineTransform || doWarpfieldTransform ) {
    for (unsigned int n=0; n<clusters.size(); n++) {
      TransformToReference(clusters[n].maxpos,trans,zvol,stdvol,full_field,doAffineTransform,doWarpfieldTransform);
      TransformToReference(clusters[n].cog,trans,zvol,stdvol,full_field,doAffineTransform,doWarpfieldTransform);
      if (copename.set()) TransformToReference(clustersCope[n].maxpos,trans,zvol,stdvol,full_field,doAffineTransform,doWarpfieldTransform);
    }
  }

  if ( doAffineTransform ) refvol = &stdvol;
  Matrix toDisplayCoord(refvol->niftivox2newimagevox_mat().i());
  if (mm.value()) {
    if (verbose.value()) 
      cout << "Using matrix : " << endl << refvol->newimagevox2mm_mat() << endl;
    toDisplayCoord=refvol->newimagevox2mm_mat();
  }
  for (unsigned int n=0; n<clusters.size(); n++) {
    MultiplyCoordinateVector(clusters[n].maxpos,toDisplayCoord);
    MultiplyCoordinateVector(clusters[n].cog,toDisplayCoord);
    if (copename.set()) MultiplyCoordinateVector(clustersCope[n].maxpos,toDisplayCoord);   // used cope before
  }
  

  string units(mm.value() ? " (mm)" : " (vox)");
  string tablehead;
  tablehead = "Cluster Index\tVoxels";
  if (pthresh.set()) tablehead += "\tP\t-log10(P)";
  string z=scalarname.value()+"-";
  if (z=="-") { z=""; }
  tablehead += "\t"+z+"MAX\t"+z+"MAX X" + units + "\t"+z+"MAX Y" + units + "\t"+z+"MAX Z" + units
    + "\t"+z+"COG X" + units + "\t"+z+"COG Y" + units + "\t"+z+"COG Z" + units;
  if (copename.set()) {
    tablehead+= "\tCOPE-MAX\tCOPE-MAX X" + units + "\tCOPE-MAX Y" + units + "\tCOPE-MAX Z" + units 
                 + "\tCOPE-MEAN";
  }

  if (!no_table.value()) cout << tablehead << endl;
  for (int n=clusters.size()-1; n>=0 && !no_table.value(); n--) {
      cout << setprecision(3) << num(n+1) << "\t" << clusters[n].size << "\t"; 
      if (!pthresh.unset()) { cout << num(clusters[n].pval) << "\t" << num(-clusters[n].logpval) << "\t"; }
        cout << num(clusters[n].maxval) << "\t" 
	   << num(clusters[n].maxpos.x) << "\t" << num(clusters[n].maxpos.y) << "\t" 
	   << num(clusters[n].maxpos.z) << "\t"
	   << num(clusters[n].cog.x) << "\t" << num(clusters[n].cog.y) << "\t" 
	   << num(clusters[n].cog.z);
      if (!copename.unset()) { 
	  cout   << "\t" << num(clustersCope[n].maxval) << "\t"
	       << num(clustersCope[n].maxpos.x) << "\t" << num(clustersCope[n].maxpos.y) << "\t" 
	       << num(clustersCope[n].maxpos.z) << "\t" << num(clustersCope[n].meanval);
	}
        cout << endl;
    
  }
  // output local maxima
  if (outlmax.set() || outlmaxim.set()) {
    string outlmaxfile="/dev/null";
    if (outlmax.set()) { outlmaxfile=outlmax.value(); }
    ofstream lmaxfile(outlmaxfile.c_str());
    if (!lmaxfile)
      cerr << "Could not open file " << outlmax.value() << " for writing" << endl;
    string scalarnm=scalarname.value();
    if (scalarnm=="") { scalarnm="Value"; }
    lmaxfile << "Cluster Index\t"+scalarnm+"\tx\ty\tz\t" << endl;
    volume<char> lmaxvol;
    copyconvert(zvol,lmaxvol);
    lmaxvol=0;
    zvol.setextrapolationmethod(zeropad);
    for (int n=clusters.size()-1; n>=0; n--) {
      vector<pair<T, triple<float> > > maxima;
      for (int z=labelim.minz(); z<=labelim.maxz(); z++)
	for (int y=labelim.miny(); y<=labelim.maxy(); y++)
	  for (int x=labelim.minx(); x<=labelim.maxx(); x++)
	    if ( checkIfLocalMaxima(clusters[n].originalLabel,labelim,zvol,numconnected.value(),x,y,z)) 
	      maxima.push_back(make_pair(zvol(x,y,z),triple<float>(x,y,z)));
      sort(maxima.rbegin(),maxima.rend());
      if (peakdist.value()>0) {
	for(unsigned int source=0;source<maxima.size();source++)
	  {
	    triple<float> sourcecoord(maxima[source].second);
	    MultiplyCoordinateVector(sourcecoord,refvol->newimagevox2mm_mat());
	    for(typename vector<pair<T, triple<float> > >::iterator clust=maxima.begin()+source+1; clust !=maxima.end();)
	    {
	      triple<float> clustcoord((*clust).second);
	      MultiplyCoordinateVector(clustcoord,refvol->newimagevox2mm_mat());
	      float dist(sqrt((sourcecoord.x-clustcoord.x)*(sourcecoord.x-clustcoord.x) + (sourcecoord.y-clustcoord.y )*(sourcecoord.y-clustcoord.y) + (sourcecoord.z-clustcoord.z)*(sourcecoord.z-clustcoord.z)));
	      clust = dist<peakdist.value() ? maxima.erase(clust) : clust+1;
	    }
	  }
      }
      maxima.resize(std::min(maxima.size(),(size_t)mx_cnt.value()));
      for(typename vector<pair<T, triple<float> > >::iterator point=maxima.begin(); point !=maxima.end(); ++point) { //output results
	lmaxvol(MISCMATHS::round((*point).second.x),
		MISCMATHS::round((*point).second.y),
		MISCMATHS::round((*point).second.z))=1;
	if ( doAffineTransform || doWarpfieldTransform ) TransformToReference((*point).second,trans,zvol,stdvol,full_field,doAffineTransform,doWarpfieldTransform);
	MultiplyCoordinateVector((*point).second, toDisplayCoord);
	lmaxfile << setprecision(3) << n+1 << "\t" << (*point).first << "\t" << 
	  (*point).second.x << "\t" << (*point).second.y << "\t" << (*point).second.z << endl;
      }
    }
    lmaxfile.close();
    if (outlmaxim.set()) {
      lmaxvol.setDisplayMaximumMinimum(0.0f,0.0f);
      save_volume(lmaxvol,outlmaxim.value());
    }
  }
}



template <class T>
int fmrib_main(int argc, char *argv[])
{
  volume<int> labelim;
  float th = thresh.value();

  // read in the volume
  volume<T> zvol, mask, cope;
  read_volume(zvol,inputname.value());
  if (verbose.value())  print_volume_info(zvol,"Zvol");
  
  if ( fractional.value() ) {
    float frac = th;
    th = frac*(zvol.robustmax() - zvol.robustmin()) + zvol.robustmin();
  }
  mask = zvol;
  mask.binarise((T) th);
  if (minv.value()) { mask = ((T) 1) - mask; }
  if (verbose.value())  print_volume_info(mask,"Mask");
  
  // Find the connected components
  labelim = connected_components(mask,numconnected.value());
  if (verbose.value())  print_volume_info(labelim,"Labelim");
  
  // process according to the output format
  vector<cluster<T> > clusters;
  get_stats(labelim,zvol,clusters,minv.value());
  int nOriginalLabels(clusters.size()+1); //0 is also a label of sorts
  if (verbose.value()) cout<<"Number of labels = "<<clusters.size()<<endl;

  // process the cope image if entered
  vector<cluster<T> > clustersCope;
  if (!copename.unset()) {
    read_volume(cope,copename.value());
    get_stats(labelim,cope,clustersCope,minv.value());
  }

 sort(clusters.rbegin(),clusters.rend());        //Sort descending for threshold purposes
 sort(clustersCope.rbegin(),clustersCope.rend());
 

  // re-threshold for p
  int nozeroclust=0;
  if (pthresh.set()) {
    if (verbose.value()) 
      cout<<"Re-thresholding with p-value"<<endl;
    Infer infer(dLh.value(), th, voxvol.value());
    if (labelim.zsize()<=1) 
      infer.setD(2); // the 2D option
    if (minclustersize.value()) {
      float pmin=1.0;
      int nmin=0;
      while (pmin>=pthresh.value()) pmin=exp(infer(++nmin)); 
      cout << "Minimum cluster size under p-threshold = " << nmin << endl;
    }
    for (unsigned int n=0; n<clusters.size(); n++) {
      clusters[n].logpval = infer(clusters[n].size)/log(10);
      clusters[n].pval = exp(clusters[n].logpval*log(10));
      if (clusters[n].pval>pthresh.value()) 
	nozeroclust++;
    }
  }
  if (verbose.value()) cout<<"Number of sub-p clusters = "<<nozeroclust<<endl;


  unsigned int n=0;
  for (; n<clusters.size(); n++) { //find first sub-p/size index
    if ( pthresh.set() && clusters[n].pval>pthresh.value() ) break;
    if ( sizethreshold.set() && clusters[n].size < sizethreshold.value() ) break;
  }
  clusters.resize(n); //remove all clusters which are either sub p or sub sizethresh
  clustersCope.resize(n);
  reverse(clusters.begin(),clusters.end());        //Ascending for output
  reverse(clustersCope.begin(),clustersCope.end());
  if (verbose.value()) {cout<<clusters.size()<<" labels in sortedidx"<<endl;}

  // print table
  print_results(clusters, clustersCope, zvol, cope, labelim);
  
  labelim.setDisplayMaximumMinimum(0,0);
  // save relevant volumes
  if ( outindex.set() ) {
    volume<int> relabeledim;
    vector<int> indexMap(nOriginalLabels,0);
    for (unsigned int n=0; n<clusters.size(); n++)
      indexMap[clusters[n].originalLabel]=n+1;
    relabel_image(labelim,relabeledim,indexMap);
    save_volume(relabeledim,outindex.value());  
  }
  if (outsize.set()) {
    volume<int> relabeledim;
    vector<int> sizeMap(nOriginalLabels,0);
    for (unsigned int n=0; n<clusters.size(); n++)
      sizeMap[clusters[n].originalLabel]=clusters[n].size;
    relabel_image(labelim,relabeledim,sizeMap);
    save_volume(relabeledim,outsize.value());
  }
  if (outmax.set()) {
    volume<T> relabeledim;
    vector<T> maxMap(nOriginalLabels,0);
    for (unsigned int n=0; n<clusters.size(); n++)
      maxMap[clusters[n].originalLabel]=clusters[n].maxval;
    relabel_image(labelim,relabeledim,maxMap);
    save_volume(relabeledim,outmax.value());
  }
  if (outmean.set()) {
    volume<float> relabeledim;
    vector<float> meanMap(nOriginalLabels,0);
    for (unsigned int n=0; n<clusters.size(); n++)
      meanMap[clusters[n].originalLabel]=clusters[n].meanval;
    relabel_image(labelim,relabeledim,meanMap);
    save_volume(relabeledim,outmean.value());
  }
  if (outpvals.set()) {
    volume<float> relabeledim;
    vector<float> pMap(nOriginalLabels,0);
    for (unsigned int n=0; n<clusters.size(); n++)
      pMap[clusters[n].originalLabel]=clusters[n].logpval;
    relabel_image(labelim,relabeledim,pMap);
    save_volume(relabeledim,outpvals.value());
    }
  if (!outthresh.unset()) {
    // Threshold the input volume st it is 0 for all non-clusters
    //   and maintains the same values otherwise
    volume<T> lcopy;
    vector<int> indexMap(nOriginalLabels,0);
    for (unsigned int n=0; n<clusters.size(); n++)
      indexMap[clusters[n].originalLabel]=n+1;
    relabel_image(labelim,lcopy,indexMap);
    lcopy.binarise(1);
    save_volume(lcopy*zvol,outthresh.value());
  }

  return 0;
}



int main(int argc,char *argv[])
{

  Tracer tr("main");

  OptionParser options(title, examples);

  try {
    options.add(inputname);
    options.add(thresh);
    options.add(outindex);
    options.add(outthresh);
    options.add(outlmax);
    options.add(outlmaxim);
    options.add(outsize);
    options.add(outmax);
    options.add(outmean);
    options.add(outpvals);
    options.add(pthresh);
    options.add(peakdist);
    options.add(copename);
    options.add(voxvol);
    options.add(dLh);
    options.add(fractional);
    options.add(numconnected);
    options.add(mm);
    options.add(minv);
    options.add(no_table);
    options.add(minclustersize);
    options.add(sizethreshold);
    options.add(transformname);
    options.add(stdvolname);
    options.add(scalarname);
    options.add(mx_cnt);
    options.add(verbose);
    options.add(help);
    options.add(warpname);
    
    options.parse_command_line(argc, argv);

    if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
      {
	options.usage();
	exit(EXIT_FAILURE);
      }
    
    if ( (!pthresh.unset()) && (dLh.unset() || voxvol.unset()) ) 
      {
	options.usage();
	cerr << endl 
	     << "Both --dlh and --volume MUST be set if --pthresh is used." 
	     << endl;
	exit(EXIT_FAILURE);
      }
    
    if ( ( !transformname.unset() && stdvolname.unset() ) ||
	 ( transformname.unset() && (!stdvolname.unset()) ) ) 
      {
	options.usage();
	cerr << endl 
	     << "Both --xfm and --stdvol MUST be set together." 
	     << endl;
	exit(EXIT_FAILURE);
      }
  }  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } 

  return call_fmrib_main(dtype(inputname.value()),argc,argv);

}

