/*  General IO functions (images and transformation files)

    Mark Jenkinson and Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 1999-2008 University of Oxford  */

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

#include "newimageio.h"

using namespace MISCMATHS;

namespace NEWIMAGE {



////////////////////////////////////////////////////////////////////////////


int  handle_read_error(int errorflag, const string& filename) {
  if ((errorflag & 1) == 1) { imthrow("ERROR:: Could not open file " + filename,22); }
  if ((errorflag & 2) == 2) { imthrow("ERROR:: Illegal NIfTI file! Inconsistent sform and qform information set in " + filename,40); }
  if ((errorflag & 4) == 4) { imthrow("ERROR:: Illegal NIfTI file! Zero determinant for sform and/or qform set in  " + filename,41); }
  return errorflag;
}



// VOLUME I/O
template <class T>
void set_volume_properties(FSLIO* IP1, volume<T>& target)
{
  float x,y,z,tr;
  FslGetVoxDim(IP1,&x,&y,&z,&tr);
  target.setdims(x,y,z);

  int sform_code, qform_code;
  mat44 smat, qmat;
  sform_code = FslGetStdXform(IP1,&smat);
  qform_code = FslGetRigidXform(IP1,&qmat);
  Matrix snewmat(4,4), qnewmat(4,4);
  for (int i=1; i<=4; i++) {
    for (int j=1; j<=4; j++) {
      snewmat(i,j) = smat.m[i-1][j-1];
      qnewmat(i,j) = qmat.m[i-1][j-1];
    }
  }
  target.set_sform(sform_code,snewmat);
  target.set_qform(qform_code,qnewmat);
  target.RadiologicalFile = (FslGetLeftRightOrder(IP1)==FSL_RADIOLOGICAL);

  short intent_code;
  float p1, p2, p3;
  FslGetIntent(IP1, &intent_code, &p1, &p2, &p3);
  target.set_intent(intent_code,p1,p2,p3);
  FslGetCalMinMax(IP1,&p1,&p2);
  target.setDisplayMinimum(p1);
  target.setDisplayMaximum(p2);
  char fileName[24];
  FslGetAuxFile(IP1,fileName);
  target.setAuxFile(string(fileName));
}

template void set_volume_properties(FSLIO* IP1, volume<char>& target);
template void set_volume_properties(FSLIO* IP1, volume<short>& target);
template void set_volume_properties(FSLIO* IP1, volume<int>& target);
template void set_volume_properties(FSLIO* IP1, volume<float>& target);
template void set_volume_properties(FSLIO* IP1, volume<double>& target);

template <class T>
int read_volumeROI(volume<T>& target, const string& filename, 
		   short& dtype, bool read_img_data,
		   int x0, int y0, int z0, int x1, int y1, int z1,
		   bool swap2radiological)
{
  // to get the whole volume use x0=y0=z0=0 and x1=y1=z1=-1
  // NB: coordinates are in "radiological" convention when swapping (i.e.
  ///    *not* the same as nifti/fslview), or untouched otherwise
  Tracer trcr("read_volumeROI");

  FSLIO *IP1;
  IP1 = NewFslOpen(filename.c_str(), "r");
  int errorflag=FslGetErrorFlag(IP1);
  if (errorflag==1) { imthrow("Failed to read volume "+filename,22); }
  short sx,sy,sz,st;
  FslGetDim(IP1,&sx,&sy,&sz,&st);
  size_t volsize=sx*sy*sz;

  T* tbuffer;
  if (read_img_data) {
    tbuffer = new T[volsize];
    if (tbuffer==0) { imthrow("Out of memory",99); }
    FslReadBuffer(IP1,tbuffer);
  } else {
    tbuffer  = new T[volsize];  // a hack to stop reinitialize from allocating memory //originally 1
  }
  target.reinitialize(sx,sy,sz,tbuffer,true);

  FslGetDataType(IP1,&dtype);
  set_volume_properties(IP1,target);

  FslClose(IP1);

  // swap to radiological if necessary
  if (swap2radiological && !target.RadiologicalFile) {
    target.makeradiological();
  }

  // now get the ROI (if necessary)
  // this is a hack until the disk reading functions integrated here
  // use -1 to signify end point
  if (x1<0) { x1=sx-1; }
  if (y1<0) { y1=sy-1; }
  if (z1<0) { z1=sz-1; }
  // truncate to known limits
  if (x0<0) { x0=0; }
  if (y0<0) { y0=0; }
  if (z0<0) { z0=0; }
  if (x1>sx-1) { x1=sx-1; }
  if (y1>sy-1) { y1=sy-1; }
  if (z1>sz-1) { z1=sz-1; }
  if (x0>x1) { x0=x1; }
  if (y0>y1) { y0=y1; }
  if (z0>z1) { z0=z1; }
  if ((x0!=0) || (y0!=0) || (z0!=0) || (x1!=sx-1) || (y1!=sy-1) || (z1!=sz-1))
    {
      target.setROIlimits(x0,y0,z0,x1,y1,z1);
      target.activateROI();
      target = target.ROI();
    }

  return errorflag;
}

template int read_volumeROI(volume<char>& target, const string& filename, 
		   short& dtype, bool read_img_data,
		   int x0, int y0, int z0, int x1, int y1, int z1,
		   bool swap2radiological);
template int read_volumeROI(volume<short>& target, const string& filename, 
		   short& dtype, bool read_img_data,
		   int x0, int y0, int z0, int x1, int y1, int z1,
		   bool swap2radiological);
template int read_volumeROI(volume<int>& target, const string& filename, 
		   short& dtype, bool read_img_data,
		   int x0, int y0, int z0, int x1, int y1, int z1,
		   bool swap2radiological);
template int read_volumeROI(volume<float>& target, const string& filename, 
		   short& dtype, bool read_img_data,
		   int x0, int y0, int z0, int x1, int y1, int z1,
		   bool swap2radiological);
template int read_volumeROI(volume<double>& target, const string& filename, 
		   short& dtype, bool read_img_data,
		   int x0, int y0, int z0, int x1, int y1, int z1,
		   bool swap2radiological);


int read_volume_size(const string& filename, 
		     int64_t& sx, int64_t& sy, int64_t& sz, int64_t& st, int64_t& s5)
{
  // read in sizes only
  Tracer trcr("read_volume_size");

  FSLIO *IP1;
  IP1 = NewFslOpen(filename.c_str(), "r");
  int errorflag=FslGetErrorFlag(IP1);
  if (errorflag==1) { imthrow("Failed to read volume "+filename,22); }

  short ssx,ssy,ssz,sst,ss5;
  FslGetDim5(IP1,&ssx,&ssy,&ssz,&sst,&ss5);
  if (sst<1) sst=1;  // make it robust to dim4=0
  sst*=ss5;  // in newimage the time dimension is used to store both dim4 and dim5 (so these are not raw)
  sx=ssx;
  sy=ssy;
  sz=ssz;
  st=sst;
  s5=ss5;
  return errorflag;
}

template <class T>
int read_volume4DROI(volume4D<T>& target, const string& filename, 
		     short& dtype, bool read_img_data,
		     int x0, int y0, int z0, int t0, 
		     int x1, int y1, int z1, int t1,
		     bool swap2radiological)
{
  // to get the whole volume use x0=y0=z0=t0=0 and x1=y1=z1=t1=-1
  // NB: coordinates are in "radiological" convention when swapping (i.e.
  ///    *not* the same as nifti/fslview), or untouched otherwise

  Tracer trcr("read_volume4DROI");

  target.destroy();

  FSLIO *IP1;
  IP1 = NewFslOpen(filename.c_str(), "r");
  int errorflag=FslGetErrorFlag(IP1);
  if (errorflag==1) { imthrow("Failed to read volume "+filename,22); }

  short sx,sy,sz,st,s5;
  FslGetDim5(IP1,&sx,&sy,&sz,&st,&s5);
  if (st<1) st=1;  // make it robust to dim4=0
  if (s5<1) s5=1;  // make it robust to dim5=0
  st*=s5;  // in newimage the time dimension is used to store both dim4 and dim5
  size_t volsize=sx*sy*sz;
  
  // use -1 to signify end point
  if (t1<0) { t1=st-1; }
  // truncate to known limits
  if (t0<0) { t0=0; }
  if (t1>st-1) { t1=st-1; }
  if (t0>t1) { t0=t1; }
  // use -1 to signify end point
  if (x1<0) { x1=sx-1; }
  if (y1<0) { y1=sy-1; }
  if (z1<0) { z1=sz-1; }
  // truncate to known limits
  if (x0<0) { x0=0; }
  if (y0<0) { y0=0; }
  if (z0<0) { z0=0; }
  if (x1>sx-1) { x1=sx-1; }
  if (y1>sy-1) { y1=sy-1; }
  if (z1>sz-1) { z1=sz-1; }
  if (x0>x1) { x0=x1; }
  if (y0>y1) { y0=y1; }
  if (z0>z1) { z0=z1; }

  volume<T> dummyvol(sx,sy,sz), tmpvol;
  // now take ROI (if necessary)
  if ((x0!=0) || (y0!=0) || (z0!=0) || (x1!=sx-1) || (y1!=sy-1) || (z1!=sz-1))
    {
      tmpvol = dummyvol;
      dummyvol.setROIlimits(x0,y0,z0,x1,y1,z1);
      dummyvol.activateROI();
      dummyvol = dummyvol.ROI();
    }
  if (t0>0) {
    if (t0>st-1) t0=st-1;
    FslSeekVolume(IP1,t0);
  }
  for (int t=t0; t<=t1; t++) {
    target.addvolume(dummyvol);
    T* tbuffer;
    if (read_img_data) {
      tbuffer = new T[volsize];
      if (tbuffer==0) { imthrow("Out of memory",99); }
      FslReadBuffer(IP1,tbuffer);
    } else {
      tbuffer = new T[volsize];  // set 1 as a bad hack to stop reinitialize from allocating memory // 
    }
    // Note that the d_owner flag = true in the following so that the
    //  control for delete is passed to the volume class
    // now take ROI (if necessary)
    if ((x0!=0) || (y0!=0) || (z0!=0) || (x1!=sx-1) || (y1!=sy-1) || (z1!=sz-1))
      {
	tmpvol.reinitialize(sx,sy,sz,tbuffer,true);
	tmpvol.setROIlimits(x0,y0,z0,x1,y1,z1);
	tmpvol.activateROI();
	target[t-t0] = tmpvol.ROI();
      } else {
	target[t-t0].reinitialize(sx,sy,sz,tbuffer,true);
      }
    set_volume_properties(IP1,target[t-t0]);
  }

  target.setROIlimits(target.limits());
  float x,y,z,tr;
  FslGetVoxDim(IP1,&x,&y,&z,&tr);
  target.setdims(x,y,z,tr);
  target.setsize5(s5);

  FslGetDataType(IP1,&dtype);

  float maximum,minimum;
  FslGetCalMinMax(IP1,&minimum,&maximum);
  target.setDisplayMinimum(minimum);
  target.setDisplayMaximum(maximum);
  char fileName[24];
  FslGetAuxFile(IP1,fileName);
  target.setAuxFile(string(fileName));
  FslClose(IP1);

  // swap to radiological if necessary
  if (swap2radiological && !target[0].RadiologicalFile) target.makeradiological();

  return errorflag;
}

template int read_volume4DROI(volume4D<char>& target, const string& filename, 
		     short& dtype, bool read_img_data,
		     int x0, int y0, int z0, int t0, 
		     int x1, int y1, int z1, int t1,
			      bool swap2radiological);
template int read_volume4DROI(volume4D<short>& target, const string& filename, 
		     short& dtype, bool read_img_data,
		     int x0, int y0, int z0, int t0, 
		     int x1, int y1, int z1, int t1,
			       bool swap2radiological);
template int read_volume4DROI(volume4D<int>& target, const string& filename, 
		     short& dtype, bool read_img_data,
		     int x0, int y0, int z0, int t0, 
		     int x1, int y1, int z1, int t1,
			      bool swap2radiological);
template int read_volume4DROI(volume4D<float>& target, const string& filename, 
		     short& dtype, bool read_img_data,
		     int x0, int y0, int z0, int t0, 
		     int x1, int y1, int z1, int t1,
			      bool swap2radiological);
template int read_volume4DROI(volume4D<double>& target, const string& filename, 
		     short& dtype, bool read_img_data,
		     int x0, int y0, int z0, int t0, 
		     int x1, int y1, int z1, int t1,
			      bool swap2radiological);

template <class T>
int save_basic_volume(const volume<T>& source, const string& filename, 
		      int filetype, bool save_orig)
{
  // if filetype < 0 then it is ignored, otherwise it overrides everything
  Tracer tr("save_basic_volume");

  bool currently_rad = source.left_right_order()==FSL_RADIOLOGICAL;
  if (!save_orig && !source.RadiologicalFile && currently_rad) const_cast< volume <T>& > (source).makeneurological();
  FSLIO *OP = NewFslOpen(filename.c_str(),"wb",filetype);
  if (OP==0) { imthrow("Failed to open volume "+filename+" for writing",23); }
  set_fsl_hdr(source,OP,1,1);
  FslWriteAllVolumes(OP,&(source(0,0,0)));
  FslClose(OP);
  if (!save_orig && !source.RadiologicalFile && currently_rad) const_cast< volume <T>& > (source).makeradiological();
  return 0;
}

template int save_basic_volume(const volume<char>& source, const string& filename, 
			       int filetype, bool save_orig);
template int save_basic_volume(const volume<short>& source, const string& filename, 
			       int filetype, bool save_orig);
template int save_basic_volume(const volume<int>& source, const string& filename, 
			       int filetype, bool save_orig);
template int save_basic_volume(const volume<float>& source, const string& filename, 
			       int filetype, bool save_orig);
template int save_basic_volume(const volume<double>& source, const string& filename, 
			       int filetype, bool save_orig);

template <class T>
int save_basic_volume4D(const volume4D<T>& source, const string& filename,
			int filetype, bool save_orig)
{
  Tracer tr("save_basic_volume4D");
  if (source.tsize()<1) return -1;
  bool currently_rad = source.left_right_order()==FSL_RADIOLOGICAL;
  if (!save_orig && !source[0].RadiologicalFile && currently_rad)  const_cast< volume4D <T>& > (source).makeneurological();
  // if filetype < 0 then it is ignored, otherwise it overrides everything
  FSLIO *OP = NewFslOpen(filename.c_str(),"wb",filetype);
  if (OP==0) { imthrow("Failed to open volume "+filename+" for writing",23); }
  set_fsl_hdr(source[0],OP,source.tsize(),source.tdim(),source.size5());
  if (filetype>=0) FslSetFileType(OP,filetype);
  FslWriteHeader(OP);
  if (source.nvoxels()>0) {
    for (int t=0; t<source.tsize(); t++) {
      FslWriteVolumes(OP,&(source[t](0,0,0)),1);
    }
  }
  FslClose(OP); 
  if (!save_orig && !source[0].RadiologicalFile && currently_rad)  const_cast< volume4D <T>& > (source).makeradiological();
  return 0;
}

template int save_basic_volume4D(const volume4D<char>& source, const string& filename,
				 int filetype, bool save_orig);
template int save_basic_volume4D(const volume4D<short>& source, const string& filename,
				 int filetype, bool save_orig);
template int save_basic_volume4D(const volume4D<int>& source, const string& filename,
				 int filetype, bool save_orig);
template int save_basic_volume4D(const volume4D<float>& source, const string& filename,
				 int filetype, bool save_orig);
template int save_basic_volume4D(const volume4D<double>& source, const string& filename,
				 int filetype, bool save_orig);

void WriteClonedHeader(FSLIO *dest, const FSLIO *src)
{
  FslCloneHeader(dest,src);
  FslSetIntensityScaling(dest,1.0,0.0); //set to (dest,0.0,0.0) for binary comparison with older versions
}

mat44 newmat2mat44(const Matrix& nmat)
{
  mat44 ret;
  for (int i=1; i<=4; i++) {
    for (int j=1; j<=4; j++) {
      ret.m[i-1][j-1] = nmat(i,j);
    }
  }
  return ret;
}


string fslbasename(const string& filename)
{
  // does string() copy safely and dispose of temporary c-string storage?
  return string(FslMakeBaseName(filename.c_str()));
}


int make_basename(string& filename)
{
  char *tmpname;
  tmpname = FslMakeBaseName(filename.c_str());
  if (tmpname==NULL) return -1;
  // this is now the basename
  filename = string(tmpname);
  // free(tmpname);  // is this safe to do?
  return 0;
}


int find_pathname(string& filename)
{
  Tracer tr("find_pathname");
  if (filename.size() < 1) return -1;
  string pathname = filename;
  int fsize = pathname.length(), indx;

  // working backwards, find '/' and remove everything after it

  indx = fsize-1;
  while ((pathname[indx] != '/') && (indx != 0))
    indx--;
  
  if (indx<fsize-1)
    pathname.erase(indx+1);
  
  filename = pathname;
  return 0;
}


bool fsl_imageexists(const string& filename) {
  return (FslFileExists(filename.c_str())!=0);
}


void check_filename(const string& basename)
{
  FSLIO* OP=FslOpen(basename.c_str(),"r");
  if (OP==NULL) {
    cerr << "ERROR: Cannot open volume " << basename << " for reading!\n";
    exit(1);
  }
}

FSLIO* NewFslOpen(const string& filename, const string& permissions, 
		  int filetype)
{
  string basename = filename;
  make_basename(basename);
  if ( basename.size()<1 ) {
    return 0;
  }

  bool writemode=false;
  if ( (permissions.find('w')!=string::npos) || 
       (permissions.find('+')!=string::npos) )  { writemode=true; }

  FSLIO* OP=FslXOpen(basename.c_str(),permissions.c_str(),filetype);
  int errorflag=FslGetErrorFlag(OP);
  if (errorflag==1) {
    imthrow("ERROR: Could not open image "+basename,22);
  }

  return OP;
}

FSLIO* NewFslOpen(const string& filename, const string& permissions)
{
  return NewFslOpen(filename,permissions,-1);
}

short closestTemplatedType(const short inputType)
{
  switch (inputType) {
  case DT_UNSIGNED_CHAR:
  case DT_INT8:
    return DT_UNSIGNED_CHAR;
  case DT_SIGNED_SHORT:
    return DT_SIGNED_SHORT;
  case DT_SIGNED_INT:
  case DT_UINT16:
    return DT_SIGNED_INT;
  case DT_FLOAT:
  case DT_UINT32:
  case DT_INT64:
  case DT_UINT64:
    return DT_FLOAT;
  case DT_DOUBLE:
  case DT_FLOAT128:
    return DT_DOUBLE;
  case DT_COMPLEX:
    cerr << "COMPLEX not supported as an independent type" << endl;
    return -1;
  default:
    cerr << "Datatype " << inputType << " is NOT supported - please check your image" << endl;
    return -1;
  }
}

short dtype(const char* T)   { return DT_UNSIGNED_CHAR; }
short dtype(const short* T)  { return DT_SIGNED_SHORT; }
short dtype(const int* T)    { return DT_SIGNED_INT; }
short dtype(const float* T)  { return DT_FLOAT; }
short dtype(const double* T) { return DT_DOUBLE; }

short dtype(const volume<char>& vol)   { return DT_UNSIGNED_CHAR; }
short dtype(const volume<short>& vol)  { return DT_SIGNED_SHORT; }
short dtype(const volume<int>& vol)    { return DT_SIGNED_INT; }
short dtype(const volume<float>& vol)  { return DT_FLOAT; }
short dtype(const volume<double>& vol) { return DT_DOUBLE; }

short dtype(const volume4D<char>& vol)   { return DT_UNSIGNED_CHAR; }
short dtype(const volume4D<short>& vol)  { return DT_SIGNED_SHORT; }
short dtype(const volume4D<int>& vol)    { return DT_SIGNED_INT; }
short dtype(const volume4D<float>& vol)  { return DT_FLOAT; }
short dtype(const volume4D<double>& vol) { return DT_DOUBLE; }

short dtype(const string& filename) 
{
  Tracer trcr("dtype");
  if ( filename.size()<1 ) return -1;
  string basename = fslbasename(filename);

  FSLIO* IP1;
  IP1 = FslOpen(basename.c_str(),"rb");
  if (IP1==NULL) {
    cerr << "Cannot open volume " << basename << " for reading!\n";
    exit(1);
  }

  short dtype;
  FslGetDataType(IP1,&dtype);
  float slope, intercept;
  int doscaling;
  doscaling = FslGetIntensityScaling(IP1,&slope,&intercept);
  if (doscaling==1) { dtype = DT_FLOAT; }
  FslClose(IP1);
  free(IP1);

  return dtype;
}

short fslFileType(const string& filename) 
{
  Tracer trcr("fslFileType");
  if ( filename.size()<1 ) return -1;
  string basename = fslbasename(filename);

  FSLIO* IP1;
  IP1 = FslOpen(basename.c_str(),"rb");
  if (IP1==NULL) {
    cerr << "Cannot open volume " << basename << " for reading!\n";
    exit(1);
  }

  short fileType(-1);
  fileType=FslGetFileType(IP1);
  FslClose(IP1);
  free(IP1);

  return fileType;
}


//////////////////////////////////////////////////////////////////////////

// COMPLEX IMAGE I/O

void FslReadComplexBuffer(FSLIO* IP, float* realbuffer, float* imagbuffer) 
{
  short sx,sy,sz,st;
  FslGetDim(IP,&sx,&sy,&sz,&st);
  size_t imagesize=sx*sy*sz;
  short type;
  FslGetDataType(IP,&type);
  switch(type)
    {
    case DT_COMPLEX:
      {
	float* sbuffer=new float[2*imagesize];
	if (sbuffer==0) { imthrow("Out of memory",99); }
	FslReadVolumes(IP,sbuffer,1);
	float *sptr=sbuffer, *rptr=realbuffer, *iptr=imagbuffer;
	for (size_t poff=0; poff<imagesize; poff++) {
	  *rptr++ = *sptr++;
	  *iptr++ = *sptr++;
	}
	delete[] sbuffer;
      }
      break;
    default:
      {  // read in the real part instead
	FslReadBuffer<float>(IP,realbuffer);
	float *iptr=imagbuffer;
	for (size_t poff=0; poff<imagesize; poff++) {
	  *iptr++ = 0.0;
	}
      }
    }
}
  
void FslWriteComplexVolume(FSLIO* OP, const float* realbuffer,
			    const float* imagbuffer) 
{
  short sx,sy,sz,st;
  FslGetDim(OP,&sx,&sy,&sz,&st);
  size_t imagesize=sx*sy*sz*1;
  float* sbuffer=new float[2*imagesize];
  if (sbuffer==0) { imthrow("Out of memory",99); }
  float *sptr=sbuffer;
  const float *rptr=realbuffer, *iptr=imagbuffer;
  for (size_t poff=0; poff<imagesize; poff++) {
    *sptr++ = *rptr++;
    *sptr++ = *iptr++;
  }
  FslWriteVolumes(OP,sbuffer,1);
  delete[] sbuffer;
}


int read_complexvolume(volume<float>& realvol, volume<float>& imagvol,
		       const string& filename,
		       bool read_img_data)
{
  Tracer trcr("read_complexvolume");
  if ( filename.size()<1 ) return -1;
  string basename = filename;
  make_basename(basename);

  FSLIO* IP1 = FslOpen(basename.c_str(), "r");
  int errorflag=FslGetErrorFlag(IP1);
  if (errorflag==1) { imthrow("Failed to read volume "+basename,22); }
  short sx,sy,sz,st;
  FslGetDim(IP1,&sx,&sy,&sz,&st);
  size_t volsize=sx*sy*sz;

  float* realbuffer=new float[volsize];
  if (realbuffer==0) { imthrow("Out of memory",99); }
  float* imagbuffer=new float[volsize];
  if (imagbuffer==0) { imthrow("Out of memory",99); }
  if (read_img_data)  FslReadComplexBuffer(IP1,realbuffer,imagbuffer);
  realvol.reinitialize(sx,sy,sz,realbuffer,true);
  imagvol.reinitialize(sx,sy,sz,imagbuffer,true);

  float x,y,z,tr;
  FslGetVoxDim(IP1,&x,&y,&z,&tr);
  realvol.setdims(x,y,z);
  imagvol.setdims(x,y,z);

  // swap to Radiological when necessary
  if (FslGetLeftRightOrder(IP1)!=FSL_RADIOLOGICAL) {
    realvol.RadiologicalFile = false;
    realvol.makeradiological();
    imagvol.RadiologicalFile = false;
    imagvol.makeradiological();
  } else {
    realvol.RadiologicalFile = true;
    imagvol.RadiologicalFile = true;
  }
  FslClose(IP1);
  return errorflag;
}


int read_complexvolume(volume<float>& realvol, volume<float>& imagvol, const string& filename)
{
  int retval = read_complexvolume(realvol,imagvol,filename,true);
  return retval;
}


int read_complexvolume(complexvolume& vol, const string& filename)
{
  return read_complexvolume(vol.re(),vol.im(),filename,true);
}


int read_complexvolume4D(volume4D<float>& realvols, volume4D<float>& imagvols,
			 const string& filename, bool read_img_data)
{
  Tracer trcr("read_complexvolume4D");
  if ( filename.size()<1 ) return -1;
  string basename = filename;
  make_basename(basename);

  FSLIO* IP1 = FslOpen(basename.c_str(), "r");
  int errorflag=FslGetErrorFlag(IP1);
  if (errorflag==1) { imthrow("Failed to read volume "+basename,22); }
  short sx,sy,sz,st;
  FslGetDim(IP1,&sx,&sy,&sz,&st);
  size_t volsize=sx*sy*sz;
  if (st<1) st=1;   //make it robust to dim4<1

  volume<float> dummyvol(sx,sy,sz);
  for (int t=0; t<st; t++) {
    realvols.addvolume(dummyvol);
    imagvols.addvolume(dummyvol);
    float* rbuffer=new float[volsize];
    if (rbuffer==0) { imthrow("Out of memory",99); }
    float* ibuffer=new float[volsize];
    if (ibuffer==0) { imthrow("Out of memory",99); }
    if (read_img_data)  FslReadComplexBuffer(IP1,rbuffer,ibuffer);
    // Note that the d_owner flag = true in the following so that the
    //  control for delete is passed to the volume class
    realvols[t].reinitialize(sx,sy,sz,rbuffer,true);
    imagvols[t].reinitialize(sx,sy,sz,ibuffer,true);
  }

  float x,y,z,tr;
  FslGetVoxDim(IP1,&x,&y,&z,&tr);
  realvols.setdims(x,y,z,tr);
  imagvols.setdims(x,y,z,tr);
  // swap to Radiological when necessary
  if (FslGetLeftRightOrder(IP1)!=FSL_RADIOLOGICAL) {
    realvols[0].RadiologicalFile = false;
    realvols.makeradiological();
    imagvols[0].RadiologicalFile = false;
    imagvols.makeradiological();
  } else {
    realvols[0].RadiologicalFile = true;
    imagvols[0].RadiologicalFile = true;
  }
  FslClose(IP1);
  return errorflag;
}


int read_complexvolume4D(volume4D<float>& realvol, volume4D<float>& imagvol, const string& filename)
{
  int retval = read_complexvolume4D(realvol,imagvol,filename,true);
  return retval;
}


int save_complexvolume(const volume<float>& realvol, const volume<float>& imagvol, const string& filename)
{
  Tracer tr("save_complexvolume");
  string basename = filename;
  make_basename(basename);
  if ( basename.size()<1 ) return -1;

  // convert back to Neurological if necessary
  if (!realvol.RadiologicalFile) { const_cast< volume <float>& > (realvol).makeneurological(); }
  if (!imagvol.RadiologicalFile) { const_cast< volume <float>& > (imagvol).makeneurological(); }

  FSLIO* OP=FslOpen(basename.c_str(),"w");
  if (OP==0) return -1;

  set_fsl_hdr(realvol,OP,1,1);

  FslSetDataType(OP, DT_COMPLEX);
  
  FslWriteHeader(OP);
  FslWriteComplexVolume(OP,&(realvol(0,0,0)),&(imagvol(0,0,0)));

  FslClose(OP); 

  // restore to original ?
  if (!realvol.RadiologicalFile) { const_cast< volume <float>& > (realvol).makeradiological(); }
  if (!imagvol.RadiologicalFile) { const_cast< volume <float>& > (imagvol).makeradiological(); }
  return 0;
}


int save_complexvolume(const complexvolume& vol, const string& filename)
{
  return save_complexvolume(vol.re(),vol.im(),filename);
}


int save_complexvolume4D(const volume4D<float>& realvols, 
			 const volume4D<float>& imagvols, 
			 const string& filename)
{
  Tracer tr("save_complexvolume4D");

  if (realvols.tsize()<=0) return -1;

  string basename = filename;
  make_basename(basename);
  if ( basename.size()<1 ) return -1;

  // convert back to Neurological if necessary
  if (!realvols[0].RadiologicalFile) { const_cast< volume4D <float>& > (realvols).makeneurological(); }
  if (!imagvols[0].RadiologicalFile) { const_cast< volume4D <float>& > (imagvols).makeneurological(); }

  FSLIO* OP=FslOpen(basename.c_str(),"w");
  if (OP==0) return -1;
  set_fsl_hdr(realvols[0],OP,realvols.tsize(),realvols.tdim(),realvols.size5());
  FslSetDataType(OP, DT_COMPLEX);

  FslWriteHeader(OP);
  
  for (int t=0; t<realvols.tsize(); t++) {
    FslWriteComplexVolume(OP,&(realvols[t](0,0,0)),&(imagvols[t](0,0,0)));
  }

  FslClose(OP); 

  // restore to original ?
  if (!realvols[0].RadiologicalFile) { const_cast< volume4D <float>& > (realvols).makeradiological(); }
  if (!imagvols[0].RadiologicalFile) { const_cast< volume4D <float>& > (imagvols).makeradiological(); }
  return 0;
}

//////////////////////////////////////////////////////////////////////////

int load_complexvolume(volume<float>& realvol, volume<float>& imagvol,
		       const string& filename)
  { return read_complexvolume(realvol,imagvol,filename); }
int load_complexvolume(complexvolume& vol, const string& filename)
  { return read_complexvolume(vol,filename); }
int load_complexvolume4D(volume4D<float>& realvol, volume4D<float>& imagvol,const string& filename)
  { return read_complexvolume4D(realvol,imagvol,filename); }

int write_complexvolume(const volume<float>& realvol, 
		       const volume<float>& imagvol, const string& filename)
  { return save_complexvolume(realvol,imagvol,filename); }
int write_complexvolume(const complexvolume& vol, const string& filename)
  { return save_complexvolume(vol,filename); }
int write_complexvolume4D(const volume4D<float>& realvol, 
			 const volume4D<float>& imagvol, 
			 const string& filename)
  { return save_complexvolume4D(realvol,imagvol,filename); }

}








