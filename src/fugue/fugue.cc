/*  fugue.cc

    Mark Jenkinson and Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 2000-2012 University of Oxford  */

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

#include "unwarpfns.h"
#include "utils/options.h"

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

using namespace Utilities;

string title="fugue \nFMRIB's Utility for Geometric Unwarping of EPIs\nCopyright(c) 2000-2008, University of Oxford (Mark Jenkinson)";
string examples="fugue -i <epi> -p <unwrapped phase map> -d <dwell-to-asym-ratio> -u <result> [options]\nfugue  -i <unwarped-image> -p <unwrapped phase map> -d <dwell-to-asym-ratio> -w <epi-like-result> [options]\nfugue -p <unwrapped phase map> -d <dwell-to-asym-ratio> --saveshift=<shiftmap> [options]";

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<bool> nocheck(string("--nocheck"), false,
		  string("turn off all checking"),
		  false, no_argument);
Option<bool> medianfilter(string("-m,--median"), false,
		  string("apply 2D median filtering"),
		  false, no_argument);
Option<bool> despike(string("--despike"), false,
		  string("apply a 2D de-spiking filter"),
		  false, no_argument);
Option<bool> pavafilter(string("--pava"), false,
		  string("apply monotonic enforcement via PAVA"),
		  false, no_argument);
Option<bool> nofill(string("--nofill"), false,
		  string("do not apply gap-filling measure to the fieldmap"),
		  false, no_argument);
Option<bool> norigidextend(string("--noextend"), false,
		  string("do not apply rigid-body extrapolation to the fieldmap"),
		  false, no_argument);
Option<bool> unmaskedfmap(string("--unmaskfmap"), false,
			  string("saves the unmasked fieldmap when using --savefmap"),
			  false, no_argument);
Option<bool> unmaskedshift(string("--unmaskshift"), false,
			  string("saves the unmasked shiftmap when using --saveshift"),
			  false, no_argument);
Option<bool> icorr(string("--icorr"), false,
		  string("apply intensity correction to unwarping (pixel shift method only)"),
		  false, no_argument);
Option<bool> icorronly(string("--icorronly"), false,
		  string("apply intensity correction only (must specify output with -u,--unwarp)"),
		  false, no_argument);
Option<bool> phaseconj(string("--phaseconj"), false,
		  string("apply phase conjugate method of unwarping"),
		  false, no_argument);
Option<bool> matrixinverse(string("--matrixinverse"), false,
			   string("apply matrix inverse method of unwarping"),
			   false, no_argument);
Option<bool> nokspace(string("--nokspace"), false,
			   string("do not use k-space forward warping"),
			   false, no_argument);
Option<string> unwarpdir(string("--unwarpdir"), "y",
			 string("specifies direction of warping [x/y/z/x-/y-/z-] (default y)"),
			 false, requires_argument);
Option<int> polynomial(string("--poly"), 0,
		   string("apply polynomial fitting of order N"),
		   false, requires_argument);
Option<int> fourier(string("--fourier"), 0,
		   string("apply Fourier (sinusoidal) fitting of order N"),
		   false, requires_argument);
Option<float> smooth2(string("--smooth2"), 0.0,
		      string("apply 2D Gaussian smoothing of sigma N (in mm)"),
		      false, requires_argument);
Option<float> smooth3(string("-s,--smooth3"), 0.0,
		      string("apply 3D Gaussian smoothing of sigma N (in mm)"),
		      false, requires_argument);
Option<string> inputname(string("-i,--in"), string(""),
			 string("filename of input volume"),
			 false, requires_argument);
Option<string> loadshift(string("--loadshift"), string(""),
			string("filename for reading pixel shift volume"),
			false, requires_argument);
Option<string> warpfname(string("-w,--warp"), string(""),
			 string("apply forward warp and save as filename"),
			 false, requires_argument);
Option<string> unwarpfname(string("-u,--unwarp"), string(""),
			   string("apply unwarping and save as filename"),
			   false, requires_argument);
Option<string> inphase(string("-p,--phasemap"), string(""),
		       string("filename for input phase image"),
		       false, requires_argument);
Option<float> dwelltoasym(string("-d,--dwelltoasym"), 0.295,
			  string("set the dwell to asym time ratio"),
			  false, requires_argument);
Option<float> dwell(string("--dwell"), 0.0,
		    string("set the EPI dwell time per phase-encode line - same as echo spacing - (sec)"),
		    false, requires_argument);
Option<float> asym(string("--asym"), 1.0,
		   string("set the fieldmap asymmetric spin echo time (sec)"),
		   false, requires_argument);
Option<float> despikethreshold(string("--despikethreshold"), 3.0,
		      string("specify the threshold for de-spiking (default=3.0)"),
		      false, requires_argument);
Option<string> saveshift(string("--saveshift"), string(""),
			 string("filename for saving pixel shift volume"),
			 false, requires_argument);
Option<string> maskfname(string("--mask"), string(""),
			 string("filename for loading valid mask"),
			 false, requires_argument);
Option<string> savefmap(string("--savefmap"), string(""),
			 string("filename for saving fieldmap (rad/s)"),
			 false, requires_argument);
Option<string> loadfmap(string("--loadfmap"), string(""),
			 string("filename for loading fieldmap (rad/s)"),
			 false, requires_argument);

string unwarpdir_value;

////////////////////////////////////////////////////////////////////////////

void pava_regularise(volume<float>& pixshift, const volume<float>& mask) 
{
  volume<float> regshift = pixshift;
  int ysize = pixshift.ysize();
  ColumnVector weights(ysize), data(ysize), pavadata(ysize);
  for (int z=pixshift.minz(); z<=pixshift.maxz(); z++) {
    for (int x=pixshift.minx(); x<=pixshift.maxx(); x++) {
      for (int y=0; y<ysize; y++) {
	// NB: data = -(shift+y) as want to enforce monotonic INCREASE
	//     but PAVA enforces DECREASE
	data(y+1) = -pixshift(x,y,z) - y;
	weights(y+1) = mask(x,y,z);
      }
      pavadata = pava(data,weights);
//        pavadata = pava(data);
      for (int y=0; y<ysize; y++) {
	pixshift(x,y,z) = -pavadata(y+1) - y;
      }
    }
  }
}


void fftshift(volume<float>& vol) {
  // does the fftshift for each 2D (z) plane separately
  volume<float> volb;
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


volume<float> warplikeepi(const volume<float>& absfmap, 
			  const volume<float>& fmap, float dwelltime)
{
  volume<float> kre, kim, epir, epii;
  kre = absfmap*0;
  kim = absfmap*0;
  epir = absfmap*0;
  epii = absfmap*0;
  float tab, Nxsize = absfmap.xsize(), Nysize=absfmap.ysize();
  // convert to k-space
  cout << "Z-limits are "<<absfmap.minz()<<" and "<<absfmap.maxz()<<endl;
  for (int z=absfmap.minz(); z<=absfmap.maxz(); z++) { 
    for (int b=0; b<Nysize; b++) {
      for (int a=0; a<Nxsize; a++) {
	tab = b + Nysize/2;  
	if (tab>Nysize)  tab-=Nysize;
	for (int n=0; n<Nysize; n++) {
	  for (int m=0; m<Nxsize; m++) {
	    float angle = tab*dwelltime*fmap(m,n,z) 
	      + 2.0*M_PI*(a*m/Nxsize + b*n/Nysize);
	    kre(a,b,z) += absfmap(m,n,z)*cos(angle);
	    kim(a,b,z) += absfmap(m,n,z)*sin(angle);
	  }
	}
      }
    }
    if (verbose.value()) { cout << "."; }
  }
  if (verbose.value()) { cout << endl; }

//    fftshift(kre);
//    fftshift(kim);

  // convert back to image space 
  for (int z=absfmap.minz(); z<=absfmap.maxz(); z++) { 
    for (int q=0; q<Nysize; q++) {
      for (int p=0; p<Nxsize; p++) {
	for (int b=0; b<Nysize; b++) {
	  for (int a=0; a<Nxsize; a++) {
	    float angle = -2.0*M_PI*(a*p/Nxsize + b*q/Nysize);
	    float cosa = cos(angle), sina = sin(angle);
	    epir(p,q,z) += kre(a,b,z)*cosa - kim(a,b,z)*sina;
	    epii(p,q,z) += kre(a,b,z)*sina + kim(a,b,z)*cosa;
	  }
	}
      }
    }
    if (verbose.value()) { cout << "."; }
  }
  if (verbose.value()) { cout << endl; }
  
  float normalisefac = 1.0/((float) Nxsize * Nysize);
  return abs(epir,epii) * normalisefac;
}



float partial_intensity(float ydash, float y0, float y1, float I0) 
{
  // calculate the part of the intensity that maps into distorted 
  //  voxel at coordinate ydash, given the half voxel end-point coords 
  //  of y0 and y1
  // That is, undistorted voxel centre maps to y0, edge maps to y1 in
  //  distorted space and undistorted intensity, I0, is mapped across
  //  distorted space with Idash in this particular voxel (centred at ydash)
  // NB: everything done in voxel coordinates
  float yL=0, yR=0, yA, yB;
  yA = Min(y0,y1);
  yB = Max(y0,y1);
  // determine yL and yR - the left and right boundaries of y0-y1 interval
  //  within the voxel with coords: ydash +/- 0.5
  if ( (yA > ydash + 0.5) || (yB < ydash - 0.5) ) { yL=0; yR =0; }
  else {
    if (yB <= ydash + 0.5) { yR = yB; } else { yR = ydash + 0.5; }
    if (yA >= ydash - 0.5) { yL = yA; } else { yL = ydash - 0.5; }
  }
  float L = yR - yL;
  // distribute intensity according to L / (yB-yA) => intensity preservation
  //  also note that 0.5 is used since there are two halves to the original
  //  undistorted voxel, each with 0.5 * I0 as intensity (additive in output)
  float Idash=0;
  if (yB - yA>0) Idash = 0.5 * I0 * L / (yB - yA);
  return Idash;
}


volume<float> warplikeepi1D_imspace(const volume<float>& undistvol, 
				   const volume<float>& fmap, float dwelltime)
{
  // MJ NOTE: ASSUMING FMAP AND UNDISTVOL ARE SAME SIZE (FOR NOW!) - TEMP!!!
  // Use the same dimensions for the output as for the undistvol
  // And only warp in the Y direction (do a swapdim before passing in to 
  //  this function for other directions)
  volume<float> distvol(undistvol);
  distvol = 0.0;
  float f2vox = fmap2pixshift_factor(undistvol,dwelltime,"y");
  float ys1, ys2, ys0, ys15, ys05, I0;

  // x,y,z in undistorted space : x,ydash,z in distorted space
  for (int z=undistvol.minz(); z<=undistvol.maxz(); z++) {
    for (int x=undistvol.minx(); x<=undistvol.maxx(); x++) {
      for (int y=undistvol.miny(); y<=undistvol.maxy(); y++) {
	// coordinates of undistorted voxel centre + left and right edges
	ys1 = y + fmap(x,y,z) * f2vox;
	if (y<undistvol.maxy()) { ys2 = y + 1 + fmap(x,y+1,z) * f2vox; }
	else { ys2 = ys1 + 1; }  
	if (y>0) { ys0 = y - 1 + fmap(x,y-1,z) * f2vox; }
	else { ys0 = ys1 - 1; }  
	ys15 = 0.5 * ( ys1 + ys2 );
	ys05 = 0.5 * ( ys1 + ys0 );
	// undistorted voxel intensity (to be distributed in output image)
	I0 = undistvol(x,y,z);
	for (int ydash=undistvol.miny(); ydash<=undistvol.maxy(); ydash++) {
	  // get partial intensities from left and right undistorted voxel parts
	  float Idash0 = partial_intensity(ydash,ys05,ys1,I0);
	  float Idash1 = partial_intensity(ydash,ys1,ys15,I0);
	  distvol(x,ydash,z) += Idash0 + Idash1;
	}
      }
    }
    if (verbose.value()) { cout << "."; }
  }
  if (verbose.value()) { cout << endl; }
  return distvol;
}


volume<float> warplikeepi1D_kspace(const volume<float>& absfmap, 
				   const volume<float>& fmap, float dwelltime)
{
  // 1D form of warp (where tab MUST ONLY DEPEND ON b)
  volume<float> kre, kim, epir, epii, fre, fim;
  kre = absfmap*0.0f;
  kim = absfmap*0.0f;
  epir = absfmap*0.0f;
  epii = absfmap*0.0f;
  fre = kre;
  fim = kim;
  float tb, Nxsize = absfmap.xsize(), Nysize=absfmap.ysize();
  // convert to k-space
  for (int z=absfmap.minz(); z<=absfmap.maxz(); z++) { 

    for (int b=0; b<Nysize; b++) {
      // assume dwell-time is 1 (factored into fmap already)
      tb = b + Nysize/2;  
      if (tb>Nysize)  tb-=Nysize;
      for (int n=0; n<Nysize; n++) {
	for (int m=0; m<Nxsize; m++) {
	  float angle = tb*dwelltime*fmap(m,n,z) + 2.0*M_PI*(b*n)/Nysize;
	  fre(m,b,z) += absfmap(m,n,z)*cos(angle);
	  fim(m,b,z) += absfmap(m,n,z)*sin(angle);
	}
      }
    }
    
    for (int b=0; b<Nysize; b++) {
      for (int a=0; a<Nxsize; a++) {
	for (int m=0; m<Nxsize; m++) {
	  float angle = 2.0*M_PI*(a*m)/Nxsize;
	  float cosa = cos(angle), sina = sin(angle);
	  kre(a,b,z) += fre(m,b,z)*cosa - fim(m,b,z)*sina;
	  kim(a,b,z) += fre(m,b,z)*sina + fim(m,b,z)*cosa;
	}
      }
    }
    
    if (verbose.value()) { cout << "."; }
  }
  if (verbose.value()) { cout << endl; }

  fre *= 0.0f;
  fim *= 0.0f;
  // convert back to image space 
  for (int z=absfmap.minz(); z<=absfmap.maxz(); z++) { 

    for (int q=0; q<Nysize; q++) {
      for (int b=0; b<Nysize; b++) {
	for (int a=0; a<Nxsize; a++) {
	  float angle = -2.0*M_PI*(b*q)/Nysize;
	  float cosa = cos(angle), sina = sin(angle);
	  fre(a,q,z) += kre(a,b,z)*cosa - kim(a,b,z)*sina;
	  fim(a,q,z) += kre(a,b,z)*sina + kim(a,b,z)*cosa;
	}
      }
    }

    for (int q=0; q<Nysize; q++) {
      for (int p=0; p<Nxsize; p++) {
	for (int a=0; a<Nxsize; a++) {
	  float angle = -2.0*M_PI*(a*p)/Nxsize;
	  float cosa = cos(angle), sina = sin(angle);
	  epir(p,q,z) += fre(a,q,z)*cosa - fim(a,q,z)*sina;
	  epii(p,q,z) += fre(a,q,z)*sina + fim(a,q,z)*cosa;
	}
      }
    }

    if (verbose.value()) { cout << "."; }
  }
  if (verbose.value()) { cout << endl; }

  float normalisefac = 1.0/((float) Nxsize * Nysize);
  return abs(epir,epii) * normalisefac;
}



volume<float> warplikeepi1D(const volume<float>& absfmap, 
			    const volume<float>& fmap, float dwelltime,
			    bool usekspace=true)
{
  if (usekspace) {
    return warplikeepi1D_kspace(absfmap,fmap,dwelltime);
  }
  // else
  return warplikeepi1D_imspace(absfmap,fmap,dwelltime);
}




volume<float> unwarpPhaseConj1D(const volume<float>& absfmap, 
			 const volume<float>& fmap, float dwelltime)
{
  // 1D form of warp (where tab MUST ONLY DEPEND ON b)
  volume<float> kre, kim, uepir, uepii, fre, fim;
  kre = absfmap*0.0f;
  kim = absfmap*0.0f;
  uepir = absfmap*0.0f;
  uepii = absfmap*0.0f;
  fre = kre;
  fim = kim;
  float tb, Nxsize=absfmap.xsize(), Nysize=absfmap.ysize();
  // convert to k-space
  for (int z=absfmap.minz(); z<=absfmap.maxz(); z++) { 

    for (int b=0; b<Nysize; b++) {
      for (int n=0; n<Nysize; n++) {
	for (int m=0; m<Nxsize; m++) {
	  float angle = 2.0*M_PI*(b*n)/Nysize;
	  fre(m,b,z) += absfmap(m,n,z)*cos(angle);
	  fim(m,b,z) += absfmap(m,n,z)*sin(angle);
	}
      }
    }
    
    for (int b=0; b<Nysize; b++) {
      for (int a=0; a<Nxsize; a++) {
	for (int m=0; m<Nxsize; m++) {
	  float angle = 2.0*M_PI*(a*m)/Nxsize;
	  float cosa = cos(angle), sina = sin(angle);
	  kre(a,b,z) += fre(m,b,z)*cosa - fim(m,b,z)*sina;
	  kim(a,b,z) += fre(m,b,z)*sina + fim(m,b,z)*cosa;
	}
      }
    }
    
    if (verbose.value()) { cout << "."; }
  }
  if (verbose.value()) { cout << endl; }

  fre *= 0.0f;
  fim *= 0.0f;
  // convert back to image space 
  for (int z=absfmap.minz(); z<=absfmap.maxz(); z++) { 

    for (int p=0; p<Nxsize; p++) {
      for (int b=0; b<Nysize; b++) {
	for (int a=0; a<Nxsize; a++) {
	  float angle = -2.0*M_PI*(a*p)/Nxsize;
	  float cosa = cos(angle), sina = sin(angle);
	  fre(p,b,z) += kre(a,b,z)*cosa - kim(a,b,z)*sina;
	  fim(p,b,z) += kre(a,b,z)*sina + kim(a,b,z)*cosa;
	}
      }
    }

    for (int q=0; q<Nysize; q++) {
      for (int p=0; p<Nxsize; p++) {
	for (int b=0; b<Nysize; b++) {
	  tb = b + Nysize/2;  
	  if (tb>Nysize)  tb-=Nysize;
	  float angle = -tb*dwelltime*fmap(p,q,z) -2.0*M_PI*(b*q)/Nysize;
	  float cosa = cos(angle), sina = sin(angle);
	  uepir(p,q,z) += fre(p,b,z)*cosa - fim(p,b,z)*sina;
	  uepii(p,q,z) += fre(p,b,z)*sina + fim(p,b,z)*cosa;
	}
      }
    }

    if (verbose.value()) { cout << "."; }
  }
  if (verbose.value()) { cout << endl; }

  float normalisefac = 1.0/((float) Nxsize * Nysize);
  return abs(uepir,uepii) * normalisefac;
}


void ComplexInvertMatrix(const Matrix& Ar, const Matrix& Ai,
			 Matrix& Br, Matrix& Bi)
{
  // simplistic version assumes that Ar is invertible
  Br = (Ar + Ai*Ar.i()*Ai).i();
  Bi = -Br*Ai*Ar.i();
}


volume<float> unwarpMatrixInverse1D(const volume<float>& warpedvol, 
				    const volume<float>& fmap, float dwelltime)
{
  float tb;
  int Nxsize=warpedvol.xsize(), Nysize=warpedvol.ysize();
  Matrix Fwarpr(Nysize,Nysize), Fwarpi(Nysize,Nysize);
  Matrix Bwarpr(Nysize,Nysize), Bwarpi(Nysize,Nysize);
  ColumnVector Warped(Nysize), Unwarpedr(Nysize), Unwarpedi(Nysize);

  volume<float> unwarpedvol;
  unwarpedvol = warpedvol;

  for (int z=0; z<warpedvol.zsize(); z++) {
    for (int x=0; x<Nxsize; x++) {
      for (int q=0; q<Nysize; q++) {
	Warped(q+1) = warpedvol(x,q,z);
	for (int n=0; n<Nysize; n++) {
	  Fwarpr(q+1,n+1)=0.0; Fwarpi(q+1,n+1)=0.0;
	  for (int b=0; b<Nysize; b++) {
	    tb = b + Nysize/2;  
	    if (tb>Nysize)  tb-=Nysize;
	    float angle = tb*dwelltime*fmap(x,n,z) + 2.0*M_PI*b*(n-q)/Nysize;
	    Fwarpr(q+1,n+1) += cos(angle);
	    Fwarpi(q+1,n+1) += sin(angle);
	  }
	}
      }
      ComplexInvertMatrix(Fwarpr,Fwarpi,Bwarpr,Bwarpi);
      Unwarpedr = Bwarpr * Warped;
      Unwarpedi = Bwarpi * Warped;
      for (int n=0; n<Nysize; n++) {
	// take the absolute value of the complex result (?)
	unwarpedvol(x,n,z) = sqrt(Sqr(Unwarpedr(n+1))+Sqr(Unwarpedi(n+1)));
      }
    }
    if (verbose.value()) { cout << "."; }
  }
  if (verbose.value()) { cout << endl; }
  return unwarpedvol;
}


///////////////////////////////////////////////////////////////////////////////

template <class T> 
void swapdirections(T& invol, const string& dir) 
{
  // This is used for swapping the required dimension to the +y axis
  // and for swapping it back
  // Each case must preserve LR order and when repeated give the identity 
      if (dir=="x")  { invol.swapdimensions(2,1,-3); }
      if (dir=="y")  { /* do nothing */ }
      if (dir=="z")  { invol.swapdimensions(-1,3,2); }
      if (dir=="x-")  { invol.swapdimensions(-2,-1,-3); }
      if (dir=="y-")  { invol.swapdimensions(-1,-2,3); }
      if (dir=="z-")  { invol.swapdimensions(-1,-3,-2); }
}


void rigidextend_y(volume<float>& pixshift, const volume<float>& mask)
{
  // Perform rigid extrapolation in the y direction
  float defval=0;

  for (int z=pixshift.minz(); z<=pixshift.maxz(); z++) {
    for (int x=pixshift.minx(); x<=pixshift.maxx(); x++) {

      int y=0, y0=0, y1=0, p;
      float i0=defval, i1=defval;
      while (y<pixshift.ysize()) {
	if (mask(x,y,z)<0.5) {   // found a region to fill/extend
	  // find end of unmasked segment (p)
	  for (p=y+1; (p<pixshift.ysize()) && (mask(x,p,z)<0.5); p++) { }
	  // set initial y coord and intensity value
	  y0 = y;
	  if (y0>0) {
	    i0 = pixshift(x,y0-1,z);
	  }
	  // set final y coord and intensity value
	  if (p<pixshift.ysize()) {
	    i1 = pixshift(x,p,z);
	    y1 = p;
	  } else {
	    i1 = i0;
	    y1 = pixshift.ysize()-1;
	  }
	  // fix initial intensity value if starting at y0==0
	  if (y0==0) {
	    i0 = i1;
	  }
	  // go through this area and set the values in the unmasked region
	  for (p=y0; p<=y1; p++) {
	    pixshift(x,p,z) = (i1-i0)*((float) (p-y0))/((float) (y1-y0)) + i0;
	  }
	  y=y1;
	} 
	y++;
      }
    }
  }
}



void regularise_pixshift(volume<float>& pixshift, 
			 volume<float>& validmask,
			 volume<float>& filledmask) 
{
  if (unwarpdir.set()) swapdirections(pixshift,unwarpdir_value);
  if (unwarpdir.set()) swapdirections(validmask,unwarpdir_value);
  if (unwarpdir.set()) swapdirections(filledmask,unwarpdir_value);

  if (!nofill.value()) {
    pixshift = extrapolate_volume(pixshift,validmask,filledmask);
  }

  if (despike.value()) {
    pixshift = masked_despike_filter2D(pixshift,filledmask,
				       despikethreshold.value());
  }

  if (medianfilter.value()) {
    pixshift = masked_median_filter2D(pixshift,filledmask);
    // pixshift = median_filter2D(pixshift);
  }

  if (polynomial.set()) {
    pixshift = polynomial_extrapolate(pixshift,validmask,polynomial.value());
    pixshift *= filledmask;
  }

  if (fourier.set()) {
    pixshift = fourier_extrapolate(pixshift,validmask,fourier.value());
    pixshift *= filledmask;
  }

  if (smooth3.set()) {
    pixshift = smooth(pixshift,smooth3.value());
  }

  if (smooth2.set()) {
    pixshift = smooth2D(pixshift,smooth2.value(),3);
  }

  if (pavafilter.value()) {
    //      mask2*=mask1;
    //      pixshift = limit_pixshift(pixshift,mask2,-0.9);
    // USE PAVA
    //  Let mask be either 0.1 or 1.0 and use it as the weighting
    volume<float> pavamask;
    pavamask = validmask * 0.9;
    pavamask += 0.1;
    pava_regularise(pixshift,pavamask);
  }

  if (!norigidextend.value()) {
    rigidextend_y(pixshift,filledmask);
  }

  if (unwarpdir.set()) swapdirections(pixshift,unwarpdir_value);       
  if (unwarpdir.set()) swapdirections(validmask,unwarpdir_value);
  if (unwarpdir.set()) swapdirections(filledmask,unwarpdir_value);

  return;
}

int get_dwell_and_asym(float& dwellval, float& asymval) {
  asymval = asym.value();
  dwellval = dwell.value();
  bool needasym = inphase.set() && savefmap.set();
  bool needdwell = loadfmap.set() && ( unwarpfname.set() || warpfname.set() 
				       || saveshift.set() );
  needdwell = needdwell || ( loadshift.set() && savefmap.set() );
  bool needratio = inphase.set() && ( unwarpfname.set() || warpfname.set() 
				       || saveshift.set() );

  if (needasym && asym.unset() && (dwelltoasym.unset() || dwell.unset())) {
    cerr << "Must set a value for asym time (or both dwell time and ratio)" 
	 << endl;
    return -1;
  }
  if (needdwell && dwell.unset() && (dwelltoasym.unset() || asym.unset())) {
    cerr << "Must set a value for dwell time (or both asym time and ratio)" 
	 << endl;
    return -1;
  }
  if (needratio && dwelltoasym.unset() && (dwell.unset() || asym.unset())) {
    cerr << "Must set the dwell to asym ratio (or both dwell and asym times)" 
	 << endl;
    return -1;
  }

  if (dwelltoasym.set()) {
    if ( asym.unset() && dwell.unset() ) {
      asymval = 1.0;
      dwellval = dwelltoasym.value();
    } else if ( asym.unset() ) {
      dwellval = dwell.value();
      asymval = dwellval / dwelltoasym.value();
    } else if ( dwell.unset() ) {
      asymval = asym.value();
      dwellval = dwelltoasym.value() * asymval;
    } else {
      cerr << "Warning: ignoring dwelltoasym value as both dwell "
	   << "and asym are separately set." << endl;
    }
  }

  if (dwellval==0.0) {
    if (needdwell) {
      cerr << "Cannot have zero dwell value!" << endl;
      return -1;
    } else {
      dwellval = 1.0;  // arbitrary value to allow division by fmap2pixshift
    }
  }

  if (needdwell) {
    if (dwellval>0.2) {
      if (nocheck.unset()) {
	cerr << "ERROR:: dwell time should be in seconds but the value of " << dwellval << "is unusually large and maybe incorrectly specified in units of milliseconds." << endl;
	cerr << "Try running with the value " << dwellval/1000.0 << endl;
	cerr << "Alternatively, to force the code to use the exact value as specified, re-run with the --nocheck option on" << endl;
	exit(EXIT_FAILURE);
      }
    }
  }

  // otherwise all is OK
  return 0;
}



int do_unwarping()
{
  volume4D<float> invol, resvol;
  volume<float> mask1, mask2;
  volume<float> pixshift, fmap;

  if (inputname.set()) {
    if (verbose.value()) { cout << "Reading input volume" << endl; }
    read_volume4D(invol,inputname.value());
  } else if (unwarpfname.set() || warpfname.set()) {
    cerr << "Must specify an input volume (-i or --in) to use (un)warping." 
	 << endl;
    return -1;
  }

  unwarpdir_value = unwarpdir.value();

  if (maskfname.set()) {
    if (verbose.value()) { cout << "Reading mask volume" << endl; }
    read_volume(mask1,maskfname.value());
    mask2 = mask1;
  }

  float dwellval, asymval;
  int retval = get_dwell_and_asym(dwellval,asymval);
  if (retval!=0)  return retval;

  if (inphase.set()) { // PROCESS UNWRAPPED PHASE MAP
    volume4D<float> uphase;
    if (verbose.value()) { cout << "Reading unwrapped phasemap" << endl; }
    read_volume4D(uphase,inphase.value());
    if (uphase.tsize()<2) {
      cerr << "Unwrappedphasemaps must contain at least two volumes!" << endl;
      return -1;
    }

    if (!maskfname.set()) {
      if (verbose.value()) { cout << "Calculating masks from phase" << endl; }
      // if no combined mask is given, must estimate from the phase images
      mask1 = 1.0f - binarise(uphase[0],-1e-6f,1e-6f);
      mask2 = 1.0f - binarise(uphase[1],-1e-6f,1e-6f);
    }
    
    if (verbose.value()) { cout << "Calculating pixel-shift map" << endl; }
    fmap = calc_fmap(uphase[0],uphase[1],mask1,mask2,asymval);
    if (inputname.unset()) { invol.addvolume(fmap); }  // needed for voxdims
    pixshift = fmap * fmap2pixshift_factor(invol[0],dwellval,unwarpdir_value);
  } 
  else { // LOAD SHIFTMAP or FIELDMAP
    if (loadshift.set()) {
      if (verbose.value()) { cout << "Reading pixel-shift map" << endl; }
      read_volume(pixshift,loadshift.value());
      if (inputname.unset()) { invol.addvolume(pixshift); }// needed for voxdims
      fmap = pixshift/fmap2pixshift_factor(invol[0],dwellval,unwarpdir_value);
    } else if (loadfmap.set()) {
      if (verbose.value()) { cout << "Reading fieldmap" << endl; }
      read_volume(fmap,loadfmap.value());
      if (inputname.unset()) { invol.addvolume(fmap); }  // needed for voxdims
      pixshift = fmap * fmap2pixshift_factor(invol[0],dwellval,unwarpdir_value);
    } else {
      cerr << "Must be able to get pixshift map from file or phasemap" << endl;
      cerr << "Re-run with either --loadshift or --phasemap set" << endl;
      return -1;
    }

    if (!maskfname.set()) {
      if (verbose.value()) { cout << "Calculating masks" << endl; }
      if (loadfmap.set()) {
	mask1 = 1.0f - binarise(fmap,-1e-6f,1e-6f);
      } else {
	mask1 = 1.0f - binarise(pixshift,-1e-6f,1e-6f);
      }
      mask2 = mask1;
    }
  }

  // Get valid and filled masks
  if (verbose.value()) { cout << "Calculating valid and filled masks" << endl; }
  volume<float> validmask, filledmask;
  validmask = mask1*mask2;
  
  filledmask = mask1;
  filledmask += mask2;
  filledmask.binarise(0.5);
  filledmask = fill_head_mask(filledmask);
  
//    if (verbose.value()) {
//      save_volume(validmask,inputname.value()+"_valid");
//      save_volume(filledmask,inputname.value()+"_mask");
//    }
  
  
  // regularise the pixshift map if required
  if (verbose.value()) { cout << "Regularising the fieldmap" << endl; }
  regularise_pixshift(pixshift,validmask,filledmask);
  fmap = pixshift / fmap2pixshift_factor(invol[0],dwellval,unwarpdir.value());

  // pixshift is now fully filtered
  if (saveshift.set()) {
    if (unmaskedshift.value()) {
      save_volume(pixshift,saveshift.value());
    } else {
      save_volume(pixshift * filledmask,saveshift.value());
    }
  }
  if (savefmap.set()) {
    if (asym.unset() && ( dwell.unset() || dwelltoasym.unset() )
	&& loadfmap.unset()) {
      if ( !(dwell.set() && asym.unset() && dwelltoasym.unset() 
	  && loadshift.set()) ) {
	cerr << "Warning in Save Fieldmap: no asym value was set, so fieldmap "
	     << "has undetermined scale factor." << endl;
      }
    }
    if (!unmaskedfmap.value()) {
      save_volume(fmap * filledmask,savefmap.value());
    } else {
      save_volume(fmap,savefmap.value());
    }
  }
  

  // THE UNWARPING PART
  if (unwarpfname.set()) {
    if (verbose.value()) { cout << "Applying unwarping" << endl; }
    resvol = invol;
    for (int t0=invol.mint(); t0<=invol.maxt(); t0++) {
      if (unwarpdir.set()) { 
	swapdirections(invol[t0],unwarpdir_value); 
	swapdirections(pixshift,unwarpdir_value);       
      }      
      if (phaseconj.value()) {
	// Phase Conjugate method
	resvol[t0] = unwarpPhaseConj1D(invol[t0],fmap,dwellval);
      } else if (matrixinverse.value()) {
	// Matrix Inverse method
	resvol[t0] = unwarpMatrixInverse1D(invol[t0],fmap,dwellval);
      } else {
	// Pixel Shifting method (default)
	if (!icorronly.value()) {
	  resvol[t0] = apply_pixshift(invol[t0],pixshift);
	} else {
	  resvol[t0] = invol[t0];
	}
	if (icorr.value() || icorronly.value()) {
	  volume<float> derivshift;
	  derivshift = yderiv(pixshift);
	  derivshift += 1.0f;
	  resvol[t0] *= derivshift;
	}
      }
      if (unwarpdir.set()) { 
	swapdirections(resvol[t0],unwarpdir_value); 
	swapdirections(invol[t0],unwarpdir_value); 
	swapdirections(pixshift,unwarpdir_value);       
      }      
    }
    // save_volume(derivshift,unwarpfname.value()+"_deriv");
    save_volume4D(resvol,unwarpfname.value());
  }


  // THE WARPING PART
  if (warpfname.set()) {
    if (verbose.value()) { cout << "Applying forward warping" << endl; }
    resvol = invol;
    for (int t0=invol.mint(); t0<=invol.maxt(); t0++) {
      // warp invol to look like an EPI acquired volume    
      if (unwarpdir.set()) { 
	swapdirections(invol[t0],unwarpdir_value); 
	swapdirections(fmap,unwarpdir_value);       
      }      
      resvol[t0] = warplikeepi1D(invol[t0],fmap,dwellval,!nokspace.value());
      if (unwarpdir.set()) { 
	swapdirections(resvol[t0],unwarpdir_value); 
	swapdirections(invol[t0],unwarpdir_value); 
	swapdirections(fmap,unwarpdir_value);       
      }      
    }    
    save_volume4D(resvol,warpfname.value());
  }
    
  return 0;
}


int main(int argc,char *argv[])
{

  Tracer tr("main");

  OptionParser options(title, examples);

  try {
    options.add(inputname);
    options.add(unwarpfname);
    options.add(warpfname);
    options.add(inphase);
    options.add(dwelltoasym);
    options.add(dwell);
    options.add(asym);
    options.add(savefmap);
    options.add(loadfmap);
    options.add(saveshift);
    options.add(loadshift);
    options.add(medianfilter);
    options.add(despike);
    options.add(nofill);
    options.add(norigidextend);
    options.add(smooth2);
    options.add(smooth3);
    options.add(polynomial);
    options.add(fourier);
    options.add(pavafilter);
    options.add(despikethreshold);
    options.add(unwarpdir);
    options.add(phaseconj);
    //options.add(matrixinverse);
    options.add(icorr);
    options.add(icorronly);
    options.add(maskfname);
    options.add(unmaskedfmap);
    options.add(unmaskedshift);
    options.add(nokspace);
    options.add(nocheck);
    options.add(verbose);
    options.add(help);
    
    options.parse_command_line(argc, argv);

    if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
      {
	options.usage();
	exit(EXIT_FAILURE);
      }
        
    if ( inphase.unset() && loadshift.unset() && loadfmap.unset() ) 
      {
	options.usage();
	cerr << endl 
	     << "Either --phasemap , --loadshift or --loadfmap MUST be used." 
	     << endl;
	exit(EXIT_FAILURE);
      }

    if (unwarpdir.value()!="x" && unwarpdir.value()!="y" && 
	unwarpdir.value()!="z" && unwarpdir.value()!="x-" && 
	unwarpdir.value()!="y-" && unwarpdir.value()!="z-") 
     {
	cerr << "Illegal value for unwarpdir!" << endl
   	     << "Use x, y, z, x-, y- or z- only." << endl;
  	exit(EXIT_FAILURE);
     }
    
  }  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } 

  int retval = do_unwarping();
  if (retval!=0) {
    cerr << endl << endl << "Error detected: try -h for help" << endl;
  }
  return retval;
}

