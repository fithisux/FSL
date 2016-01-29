/*  imfft.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2001 University of Oxford  */

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


// FFT routines for 3D images


#include "imfft.h"

using namespace MISCMATHS;
using namespace NEWMAT;


namespace NEWIMAGE {

//////////////////////////////////////////////////////////////////////////////

template <class T>
int ifft(volume<T>& revol, volume<T>& imvol, bool transformz=true)
{
  ColumnVector fvr, fvi, vecr, veci;
  // do the transform in x
  int xoff = revol.minx()-1;
  vecr.ReSize(revol.maxx()-xoff);
  veci.ReSize(revol.maxx()-xoff);
  for (int z=revol.minz(); z<=revol.maxz(); z++) {
    for (int y=revol.miny(); y<=revol.maxy(); y++) {
      for (int x=revol.minx(); x<=revol.maxx(); x++) {
	vecr(x-xoff) = revol(x,y,z);
	veci(x-xoff) = imvol(x,y,z);
      }
      FFTI(vecr,veci,fvr,fvi);
      for (int x=revol.minx(); x<=revol.maxx(); x++) {
	revol(x,y,z) = fvr(x-xoff);
	imvol(x,y,z) = fvi(x-xoff);
      }
    }
  }
  // do the transform in y
  int yoff = revol.miny()-1;
  vecr.ReSize(revol.maxy()-yoff);
  veci.ReSize(revol.maxy()-yoff);
  for (int z=revol.minz(); z<=revol.maxz(); z++) {
    for (int x=revol.minx(); x<=revol.maxx(); x++) {
      for (int y=revol.miny(); y<=revol.maxy(); y++) {
	vecr(y-yoff) = revol(x,y,z);
	veci(y-yoff) = imvol(x,y,z);
      }
      FFTI(vecr,veci,fvr,fvi);
      for (int y=revol.miny(); y<=revol.maxy(); y++) {
	revol(x,y,z) = fvr(y-yoff);
	imvol(x,y,z) = fvi(y-yoff);
      }
    }
  }

  if (transformz && ((revol.maxz()-revol.minz())>0)) {
    // do the transform in z
    int zoff = revol.minz()-1;
    vecr.ReSize(revol.maxz()-zoff);
    veci.ReSize(revol.maxz()-zoff);
    for (int x=revol.minx(); x<=revol.maxx(); x++) {
      for (int y=revol.miny(); y<=revol.maxy(); y++) {
	for (int z=revol.minz(); z<=revol.maxz(); z++) {
	  vecr(z-zoff) = revol(x,y,z);
	  veci(z-zoff) = imvol(x,y,z);
	}
	FFTI(vecr,veci,fvr,fvi);
	for (int z=revol.minz(); z<=revol.maxz(); z++) {
	  revol(x,y,z) = fvr(z-zoff);
	  imvol(x,y,z) = fvi(z-zoff);
	}
      }
    }
  }
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int ifft3(complexvolume& vol)
{
  return ifft(vol.re(),vol.im(),true);
}

int ifft2(complexvolume& vol)
{
  return ifft(vol.re(),vol.im(),false);
}

int ifft3(const complexvolume& vin, complexvolume& vout)
{
  // set up dimensions
  vout = vin;
  return ifft(vout.re(),vout.im(),true);
}

int ifft2(const complexvolume& vin, complexvolume& vout)
{
  // set up dimensions
  vout = vin;
  return ifft(vout.re(),vout.im(),false);
}

int ifft2(volume<float>& realvol, volume<float>& imagvol)
{  return ifft(realvol,imagvol,false); }

int ifft2(const volume<float>& realvin, const volume<float>& imagvin,
	  volume<float>& realvout, volume<float>& imagvout)
{
  // set up dimensions
  realvout = realvin;
  imagvout = imagvin;
  return ifft(realvout,imagvout,false);
}

int ifft3(volume<float>& realvol, volume<float>& imagvol)
{  return ifft(realvol,imagvol,true); }

int ifft3(const volume<float>& realvin, const volume<float>& imagvin,
	  volume<float>& realvout, volume<float>& imagvout)
{
  // set up dimensions
  realvout = realvin;
  imagvout = imagvin;
  return ifft(realvout,imagvout,true);
}

int ifft2(volume<double>& realvol, volume<double>& imagvol)
{  return ifft(realvol,imagvol,false); }

int ifft2(const volume<double>& realvin, const volume<double>& imagvin,
	  volume<double>& realvout, volume<double>& imagvout)
{
  // set up dimensions
  realvout = realvin;
  imagvout = imagvin;
  return ifft(realvout,imagvout,false);
}

int ifft3(volume<double>& realvol, volume<double>& imagvol)
{  return ifft(realvol,imagvol,true); }

int ifft3(const volume<double>& realvin, const volume<double>& imagvin,
	  volume<double>& realvout, volume<double>& imagvout)
{
  // set up dimensions
  realvout = realvin;
  imagvout = imagvin;
  return ifft(realvout,imagvout,true);
}


///////////////////////////////////////////////////////////////////////////////

template <class T>
int fft(volume<T>& revol, volume<T>& imvol, bool transformz=true)
{
  ColumnVector fvr, fvi, vecr, veci;
  // do the transform in x
  int xoff = revol.minx()-1;
  vecr.ReSize(revol.maxx()-xoff);
  veci.ReSize(revol.maxx()-xoff);
  for (int z=revol.minz(); z<=revol.maxz(); z++) {
    for (int y=revol.miny(); y<=revol.maxy(); y++) {
      for (int x=revol.minx(); x<=revol.maxx(); x++) {
	vecr(x-xoff) = revol(x,y,z);
	veci(x-xoff) = imvol(x,y,z);
      }
      FFT(vecr,veci,fvr,fvi);
      for (int x=revol.minx(); x<=revol.maxx(); x++) {
	revol(x,y,z) = fvr(x-xoff);
	imvol(x,y,z) = fvi(x-xoff);
      }
    }
  }
  // do the transform in y
  int yoff = revol.miny()-1;
  vecr.ReSize(revol.maxy()-yoff);
  veci.ReSize(revol.maxy()-yoff);
  for (int z=revol.minz(); z<=revol.maxz(); z++) {
    for (int x=revol.minx(); x<=revol.maxx(); x++) {
      for (int y=revol.miny(); y<=revol.maxy(); y++) {
	vecr(y-yoff) = revol(x,y,z);
	veci(y-yoff) = imvol(x,y,z);
      }
      FFT(vecr,veci,fvr,fvi);
      for (int y=revol.miny(); y<=revol.maxy(); y++) {
	revol(x,y,z) = fvr(y-yoff);
	imvol(x,y,z) = fvi(y-yoff);
      }
    }
  }

  if (transformz && ((revol.maxz()-revol.minz())>0)) {
    // do the transform in z
    int zoff = revol.minz()-1;
    vecr.ReSize(revol.maxz()-zoff);
    veci.ReSize(revol.maxz()-zoff);
    for (int x=revol.minx(); x<=revol.maxx(); x++) {
      for (int y=revol.miny(); y<=revol.maxy(); y++) {
	for (int z=revol.minz(); z<=revol.maxz(); z++) {
	  vecr(z-zoff) = revol(x,y,z);
	  veci(z-zoff) = imvol(x,y,z);
	}
	FFT(vecr,veci,fvr,fvi);
	for (int z=revol.minz(); z<=revol.maxz(); z++) {
	  revol(x,y,z) = fvr(z-zoff);
	  imvol(x,y,z) = fvi(z-zoff);
	}
      }
    }
  }
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

int fft3(complexvolume& vol)
{
  return fft(vol.re(),vol.im(),true);
}

int fft2(complexvolume& vol)
{
  return fft(vol.re(),vol.im(),false);
}

int fft3(const complexvolume& vin, complexvolume& vout)
{
  // set up dimensions
  vout = vin;
  return fft(vout.re(),vout.im(),true);
}

int fft2(const complexvolume& vin, complexvolume& vout)
{
  // set up dimensions
  vout = vin;
  return fft(vout.re(),vout.im(),false);
}

int fft2(volume<float>& realvol, volume<float>& imagvol)
{  return fft(realvol,imagvol,false); }

int fft2(const volume<float>& realvin, const volume<float>& imagvin,
	 volume<float>& realvout, volume<float>& imagvout)
{
  // set up dimensions
  realvout = realvin;
  imagvout = imagvin;
  return fft(realvout,imagvout,false);
}

int fft3(volume<float>& realvol, volume<float>& imagvol)
{  return fft(realvol,imagvol,true); }

int fft3(const volume<float>& realvin, const volume<float>& imagvin,
	 volume<float>& realvout, volume<float>& imagvout)
{
  // set up dimensions
  realvout = realvin;
  imagvout = imagvin;
  return fft(realvout,imagvout,true);
}

int fft2(volume<double>& realvol, volume<double>& imagvol)
{  return fft(realvol,imagvol,false); }

int fft2(const volume<double>& realvin, const volume<double>& imagvin,
	 volume<double>& realvout, volume<double>& imagvout)
{
  // set up dimensions
  realvout = realvin;
  imagvout = imagvin;
  return fft(realvout,imagvout,false);
}

int fft3(volume<double>& realvol, volume<double>& imagvol)
{  return fft(realvol,imagvol,true); }

int fft3(const volume<double>& realvin, const volume<double>& imagvin,
	 volume<double>& realvout, volume<double>& imagvout)
{
  // set up dimensions
  realvout = realvin;
  imagvout = imagvin;
  return fft(realvout,imagvout,true);
}

//////////////////////////////////////////////////////////////////////////////

template <class T>
void fftshift(volume<T>& vol, bool transformz) {
  if (transformz) {
    cerr << "WARNING:: fftshift not implemented in 3D - doing 2D instead"<<endl;
  }
  // does the fftshift for each 2D (z) plane separately
  volume<T> volb;
  volb = vol;
  int Na, Nb, mida, midb;
  Na = vol.xsize();
  Nb = vol.ysize();
  mida = (Na+1)/2;  // effectively a ceil()
  midb = (Nb+1)/2;  // effectively a ceil()

  for (int z=vol.minz(); z<=vol.maxz(); z++) { 

    for (int a=0; a<Na; a++) {
      for (int b=midb; b<=Nb-1; b++) {
	vol(a,b-midb,z) = volb(a,b,z);
      }
      for (int b=0; b<=midb-1; b++) {
	vol(a,b+Nb-midb,z) = volb(a,b,z);
      }
    }

    volb = vol;

    for (int b=0; b<Nb; b++) {
      for (int a=mida; a<=Na-1; a++) {
	vol(a-mida,b,z) = volb(a,b,z);
      }
      for (int a=0; a<=mida-1; a++) {
	vol(a+Na-mida,b,z) = volb(a,b,z);
      }
    }

  }
}


void fftshift(complexvolume& vol) {
  fftshift(vol.re(),false);
  fftshift(vol.im(),false);
}

void fftshift(volume<double>& vol) {
  return fftshift(vol,false);
}

void fftshift(volume<float>& vol) {
  return fftshift(vol,false);
}

}
