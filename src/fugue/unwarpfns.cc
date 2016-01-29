/*  unwarpfns.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2000 University of Oxford  */

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

#include "unwarpfns.h"

int maximum_position(const ColumnVector& vec)
{
  if (vec.Nrows()<=0) return -1;
  int idx=1;
  float max=vec(1);
  for (int n=2; n<=vec.Nrows(); n++) {
    if (vec(n)>max) { max = vec(n); idx = n; }
  }
  return idx;
}

int maximum_position(const Matrix& vec)
{
  if (vec.Nrows()<=0) return -1;
  int idx=1;
  float max=vec(1,1);
  for (int n=2; n<=vec.Nrows(); n++) {
    if (vec(n,1)>max) { max = vec(n,1); idx = n; }
  }
  return idx;
}


ColumnVector pava(const ColumnVector& data, const ColumnVector& weight0) 
{
  Tracer tr("pava");
  
  int length = (int) data.Nrows();
  ColumnVector values = data;
  ColumnVector weights = weight0;
  ColumnVector index(length);
  for(int j=1; j<=length; j++)  index(j) = j;
  
  bool anyviolators = true;
  // the PAVA
  while(anyviolators) {
    anyviolators = false;
    
    for(int k=2; k<=values.Nrows(); k++) {
      if (values(k) > values(k-1)) {
	anyviolators = true;
	values(k-1) = (values(k-1)*weights(k-1) + values(k)*weights(k))
	  / (weights(k-1) + weights(k));
	values = values.Rows(1,k-1) & values.Rows(k+1,values.Nrows());
	weights(k-1) = weights(k) + weights(k-1);
	weights = weights.Rows(1,k-1) & weights.Rows(k+1,weights.Nrows());
	
	for(int j=1; j<=length; j++) {
	  if (index(j)>=k)  index(j) = index(j)-1;
	}
	break;
      }
    }
  }

  // put the results back
  ColumnVector result = data;
  for(int j=1; j<=length; j++) {
    result(j) = values(MISCMATHS::round(double(index(j))));
  }
  return result;
}


ColumnVector pava(const ColumnVector& data) 
{
  ColumnVector weights=data;
  weights=1.0;
  return pava(data,weights);
}


/////////////////////////////////////////////////////////////////////////

volume<float> despike_filter2D(const volume<float>& vol, float threshold)
{
  volume<float> outvol(vol);
  for (int z=outvol.minz(); z<=outvol.maxz(); z++) {
    for (int y=outvol.miny(); y<=outvol.maxy(); y++) {
      for (int x=outvol.minx(); x<=outvol.maxx(); x++) {
	ColumnVector vals(8);
	int count=1;
	for (int y0=y-1; y0<=y+1; y0++) {
	  for (int x0=x-1; x0<=x+1; x0++) {
	    if ( (y0!=y) || (x0!=x) ) vals(count++) = vol(x0,y0,z);
	  }
	}
	SortAscending(vals);
	float range = vals(8) - vals(1);
	float median = vals(5);
	float cval = vol(x,y,z);
	if ( (fabs(range)<1e-9) || 
	     (fabs((cval - median)/range) > threshold) ) {
	  outvol(x,y,z) = median;
	}
      }
    }
  }
  return outvol;
}


volume<float> median_filter2D(const volume<float>& vol)
{
  volume<float> outvol(vol);
  for (int z=outvol.minz(); z<=outvol.maxz(); z++) {
    for (int y=outvol.miny(); y<=outvol.maxy(); y++) {
      for (int x=outvol.minx(); x<=outvol.maxx(); x++) {
	ColumnVector vals(9);
	int count=1;
	for (int y0=y-1; y0<=y+1; y0++) {
	  for (int x0=x-1; x0<=x+1; x0++) {
	    vals(count++) = vol(x0,y0,z);
	  }
	}
	SortAscending(vals);
	outvol(x,y,z) = vals(5);
      }
    }
  }
  return outvol;
}

volume<float> masked_despike_filter2D(const volume<float>& vol, 
				      const volume<float>& mask,
				      float threshold)
{
  volume<float> outvol(vol);
  for (int z=outvol.minz(); z<=outvol.maxz(); z++) {
    for (int y=outvol.miny(); y<=outvol.maxy(); y++) {
      for (int x=outvol.minx(); x<=outvol.maxx(); x++) {
	ColumnVector vals(9);
	int count=1;
	for (int y0=y-1; y0<=y+1; y0++) {
	  for (int x0=x-1; x0<=x+1; x0++) {
	    if (mask(x0,y0,z)>0.5) {
	      if ( (y0!=y) || (x0!=x) ) vals(count++) = vol(x0,y0,z);
	    }
	  }
	}
	if (count>1) {
	  ColumnVector littlevals;
	  littlevals = vals.SubMatrix(1,count-1,1,1);
	  SortAscending(littlevals);
	  float median = littlevals(Max(1,(count-1)/2));
	  float range = littlevals(count-1) - littlevals(1);
	  float cval = vol(x,y,z);
	  if ( (fabs(range)<1e-9) || 
	       (fabs((cval - median)/range) > threshold) ) {
	    outvol(x,y,z) = median;
	  }
	}
      }
    }
  }
  return outvol;
}


volume<float> masked_median_filter2D(const volume<float>& vol, 
				     const volume<float>& mask)
{
  volume<float> outvol(vol);
  for (int z=outvol.minz(); z<=outvol.maxz(); z++) {
    for (int y=outvol.miny(); y<=outvol.maxy(); y++) {
      for (int x=outvol.minx(); x<=outvol.maxx(); x++) {
	ColumnVector vals(9);
	int count=1;
	for (int y0=y-1; y0<=y+1; y0++) {
	  for (int x0=x-1; x0<=x+1; x0++) {
	    if (mask(x0,y0,z)>0.5) {
	      vals(count++) = vol(x0,y0,z);
	    }
	  }
	}
	if (count>1) {
	  ColumnVector littlevals;
	  littlevals = vals.SubMatrix(1,count-1,1,1);
	  SortAscending(littlevals);
	  outvol(x,y,z) = littlevals(Max(1,(count-1)/2));
	} else {
	  outvol(x,y,z) = 0.0;  // THE DEFAULT CASE (what is most sensible?)
	}
      }
    }
  }
  return outvol;
}


/////////////////////////////////////////////////////////////////////////

float basic_mask_threshold(const volume<float>& absmap)
{
  float thresh, t2, t98;
  t2 = absmap.percentile(0.02);
  t98 = absmap.percentile(0.98);
  thresh = 0.1*(t98 - t2) + t2;
  return thresh;
}


volume<float> make_basic_head_mask(const volume<float>& absmap, float thresh)
{
  volume<float> mask;
  mask = binarise(absmap,thresh);
  return mask;
}


volume<float> largest_connected_component(const volume<float>& mask)
{ 
  mask.setextrapolationmethod(extraslice);
  volume<int> label;
  ColumnVector hist;
  label = connected_components(mask,6);
  // check to see if there are ANY valid voxels
  if ( (label.max()<=0) && (mask(0,0,0)==0.0) )  return mask*0.0f;
  // only keep component with most voxels
  hist = label.histogram(label.max()+1,0,label.max()+1);
  hist(1) = 0;  // ignore background
  int brainidx = maximum_position(hist) - 1;
  label.binarise(brainidx,brainidx);
  volume<float> singlecomp;
  copyconvert(label,singlecomp);
  return singlecomp;
}


volume<float> make_head_mask(const volume<float>& absmap, float thresh) 
{ 
  volume<float> mask;
  mask = make_basic_head_mask(absmap,thresh);
  mask = largest_connected_component(mask);
  return mask;
}


volume<float> make_head_mask2D(const volume<float>& absmap, float thresh) 
{ 
  volume<float> mask(absmap);
  mask = 0.0;
  // ONLY USES 2D CONNECTIVITY

  int zmin = absmap.minz(), zmax = absmap.maxz();
  absmap.activateROI();
  mask.activateROI();

  for (int z=zmin; z<=zmax; z++) {
    volume<float> aslice, mslice;
    // get slices from absmap
    absmap.setROIlimits(absmap.minx(),absmap.miny(),z,
			absmap.maxx(),absmap.maxy(),z);
    aslice = absmap.ROI();
    // find the mask
    mslice = make_head_mask(aslice,thresh);
    // set the corresponding slice of mask to mslice
    mask.setROIlimits(mask.minx(),mask.miny(),z,
		      mask.maxx(),mask.maxy(),z);
    mask.copyROIonly(mslice);
  }
  
  // restore previous ROI status...
  absmap.deactivateROI();
  
  return mask;
}


volume<float> fill_head_mask(const volume<float>& origmask)
{
  volume<int> label, label2;
  ColumnVector hist;
  volume<float> mask = origmask;
  simple_dilate(mask);
  simple_erode(mask);
  mask.setextrapolationmethod(extraslice);
  // invert mask so that is keeps all the bits NOT connected to the background
  mask = ((float) 1.0) - mask;
  // some of the bits are false positives from the initial thresholding
  //   so get rid of them by only retaining the largest connected component
  label = connected_components(mask,6);
  // background = component with most voxels (ignoring brain)
  hist = label.histogram(label.max()+1,0,label.max()+1);
  hist(1) = 0;  // ignore brain
  int background = maximum_position(hist) - 1;
  label.binarise(background,background);
  // invert the mask so background is zero
  copyconvert(1 - label,mask);
  // now eliminate all but the single biggest component (the brain)
  label = connected_components(mask,6);
  hist = label.histogram(label.max()+1,0,label.max()+1);
  hist(1) = 0;  // ignore background
  int brainidx = maximum_position(hist) - 1;
  label.binarise(brainidx,brainidx);
  copyconvert(label,mask);
  return mask;
}


volume<float> make_filled_head_mask(const volume<float>& absmap) 
{ 
  volume<float> mask;
  float thresh = basic_mask_threshold(absmap);
  mask = make_basic_head_mask(absmap, thresh);
  mask = fill_head_mask(mask);
  return mask;
}

void fill_holes(volume<float>& vol, const volume<float>& holemask, 
		const volume<float>& mask) 
{
  // initially just find average of nearest valid entries and fill
  //  invalid voxels with these
  mask.setextrapolationmethod(zeropad);
  holemask.setextrapolationmethod(zeropad);
  int maxradius = Min(Min(vol.maxx()-vol.minx(),vol.maxy()-vol.miny()),
		      vol.maxz()-vol.minz());
  float defaultval = vol.percentile(0.5);
  for (int z=vol.minz(); z<=vol.maxz(); z++) {
    for (int y=vol.miny(); y<=vol.maxy(); y++) {
      for (int x=vol.minx(); x<=vol.maxx(); x++) {
	if (holemask(x,y,z)>0.5) {
	  int radius=1, count=0;
	  float sum=0.0;
	  // keep searching in greater and greater boxes for valid neighbours
	  while ((radius<maxradius) && (count==0) ) {
	    for (int z0=z-radius; z0<=z+radius; z0++) {
	      for (int y0=y-radius; y0<=y+radius; y0++) {
		for (int x0=x-radius; x0<=x+radius; x0++) {
		  if ( (holemask(x0,y0,z0)<0.5) && (mask(x0,y0,z0)>0.5) ) {
		    sum+=vol(x0,y0,z0);
		    count++;
		  }
		}
	      }
	    }
	    radius++;
	  } // end while
	  if (count==0) {
	    vol(x,y,z) = defaultval;
	  } else {
	    vol(x,y,z) = sum/((float) count);
	  }
	}
      }
    }
  }
  // now smooth the whole image, but only include the smoothed result in
  //  the holes
  ColumnVector blurmask(3);
  blurmask << 0.2 << 0.6 << 0.2;
  volume<float> vol2, volnorm;
  vol2 = convolve_separable(vol,blurmask,blurmask,blurmask);
  volnorm = convolve_separable(mask + holemask,blurmask,blurmask,blurmask);
  volnorm += 0.0001;  // avoid division by zeros
  vol2 /= volnorm;
  // now calculate : mask*vol2 + (1-mask)*vol
  vol2 -= vol;
  vol2 *= holemask;
  vol += vol2;
}


volume<float> extrapolate_volume(const volume<float>& datavol,
				 const volume<float>& origmask,
				 const volume<float>& enlargedmask)
{
  // create a mask that is one for all the within-head holes
  volume<float> holemask;
  holemask = enlargedmask*(1.0f - origmask);
  // now fill in the holes in uphase
  volume<float> extrapvol = datavol;
  fill_holes(extrapvol,holemask,origmask);
  return extrapvol;
}


ColumnVector extrapolate_coefficients(const Matrix& BB,
				      const ColumnVector& IB)
{
  float thresh = 1e-12;
  DiagonalMatrix D, Dinv;
  Matrix U, V;
  SVD(BB, D, U, V);
  Dinv = D;
  thresh *= D.MaximumAbsoluteValue();
  for (int n=1; n<=D.Nrows(); n++) {
    if (fabs(D(n))>thresh) {
      Dinv(n) = 1.0/D(n);
    } else{
      Dinv(n) = 0.0;
    }
  }
  ColumnVector C;
  C = U.t() * IB;
  C = Dinv * C;
  C = V * C;
  return C;
}


volume<float> polynomial_extrapolate(const volume<float>& datavol,
				     const volume<float>& validmask,
				     int n, bool verbose)
{
  int numb = (n+1)*(n+1)*(n+1);
  vector<volume<float> > basisx(n+1), basisy(n+1), basisz(n+1);
  vector<int> ix(numb), iy(numb), iz(numb);
  Matrix BB(numb,numb);
  ColumnVector IB(numb);

  // generate basis images
  if (verbose) cout << "Generating basis images" << endl;

  int idx = 0;
  for (int mz=0; mz<=n; mz++) {
    for (int my=0; my<=n; my++) {
      for (int mx=0; mx<=n; mx++) {
	ix[idx]=mx; iy[idx]=my; iz[idx]=mz;
	idx++;
      }
    }
  }

  if (verbose) cout << "Generating basis set" << endl;
  volume<float> temp(datavol);
  float midx = (temp.maxx() + temp.minx()) * 0.5;
  float midy = (temp.maxy() + temp.miny()) * 0.5;
  float midz = (temp.maxz() + temp.minz()) * 0.5;
  for (int dir=1; dir<=3; dir++) {
    for (int m=0; m<=n; m++) {
      temp = datavol;
      for (int z=temp.minz(); z<=temp.maxz(); z++) {
	for (int y=temp.miny(); y<=temp.maxy(); y++) {
	  for (int x=temp.minx(); x<=temp.maxx(); x++) {
	    if (dir==1)  temp(x,y,z) = MISCMATHS::pow(x-midx,int(m))/
			   MISCMATHS::pow(16.0,int(m));  
	    if (dir==2)  temp(x,y,z) = MISCMATHS::pow(y-midy,int(m))/
			   MISCMATHS::pow(16.0,int(m));
	    if (dir==3)  temp(x,y,z) = MISCMATHS::pow(z-midz,int(m))/
			   MISCMATHS::pow(16.0,int(m));
	  }
	}
      }
      if (dir==1) basisx[m] = temp;
      if (dir==2) basisy[m] = temp;
      if (dir==3) basisz[m] = temp;
      if (verbose) cout << ".";
    }
  }
  if (verbose) cout << endl << "Forming dot products" << endl;

  // form basis dot products
  for (int m=0; m<numb; m++) {
    temp = datavol * basisx[ix[m]];
    temp *= basisy[iy[m]];
    temp *= basisz[iz[m]];
    temp *= validmask;
    IB(m+1) = temp.sum();
    for (int m2=m; m2<numb; m2++) {
      temp = basisx[ix[m]] * basisx[ix[m2]];
      temp *= basisy[iy[m]];
      temp *= basisy[iy[m2]];
      temp *= basisz[iz[m]];
      temp *= basisz[iz[m2]];
      temp *= validmask;
      BB(m+1,m2+1) = temp.sum();
      BB(m2+1,m+1) = temp.sum();
    }
    if (verbose) cout << ".";
  }
  if (verbose) cout << endl;
  
  // find least-squares solution
  if (verbose) cout << "Finding least-squares solution" << endl;
  ColumnVector C;
  C = extrapolate_coefficients(BB,IB);
  
  // generate the resultant extrapolation
  if (verbose) cout << "Regenerating result" << endl;
  volume<float> resvol;
  resvol = datavol * 0.0f;
  for (int m=0; m<numb; m++) {
    resvol += basisx[ix[m]] * basisy[iy[m]] * basisz[iz[m]] * C(m+1);
  }
//    resvol -= resvol * validmask;
//    resvol += datavol * validmask;

  return resvol;
}



volume<float> fourier_extrapolate(const volume<float>& datavol,
				  const volume<float>& validmask,
				  int n, bool verbose)
{
  int numb = (2*n+1)*(2*n+1)*(2*n+1);
  vector<volume<float> >  basisx(2*n+1), basisy(2*n+1), basisz(2*n+1);
  vector<int> ix(numb), iy(numb), iz(numb);
  Matrix BB(numb,numb);
  ColumnVector IB(numb);

  // generate basis images
  if (verbose) cout << "Generating basis images" << endl;
  int idx = 0;
  for (int mz=0; mz<=2*n; mz++) {
    for (int my=0; my<=2*n; my++) {
      for (int mx=0; mx<=2*n; mx++) {
	ix[idx]=mx; iy[idx]=my; iz[idx]=mz;
	idx++;
      }
    }
  }


  volume<float> temp(datavol);
  float onlenx = 1.0/(temp.maxx() - temp.minx() + 1.0);
  float onleny = 1.0/(temp.maxy() - temp.miny() + 1.0);
  float onlenz = 1.0/(temp.maxz() - temp.minz() + 1.0);

  for (int dir=1; dir<=3; dir++) {
    for (int m=0; m<=2*n; m++) {
      temp = datavol;
      for (int z=temp.minz(); z<=temp.maxz(); z++) {
	for (int y=temp.miny(); y<=temp.maxy(); y++) {
	  for (int x=temp.minx(); x<=temp.maxx(); x++) {
	    if (dir==1) {
	      if (m<=n) { temp(x,y,z) = cos(2.0*M_PI*x*m*onlenx); }
	      else { temp(x,y,z) = sin(2.0*M_PI*x*(m-n)*onlenx); }
	    }
	    if (dir==2) {
	      if (m<=n) { temp(x,y,z) = cos(2.0*M_PI*y*m*onleny); }
	      else { temp(x,y,z) = sin(2.0*M_PI*y*(m-n)*onleny); }
	    }
	    if (dir==3) {
	      if (m<=n) { temp(x,y,z) = cos(2.0*M_PI*z*m*onlenz); }
	      else { temp(x,y,z) = sin(2.0*M_PI*z*(m-n)*onlenz); }
	    }
	  }
	}
      }
      if (dir==1) basisx[m] = temp;
      if (dir==2) basisy[m] = temp;
      if (dir==3) basisz[m] = temp;
      if (verbose) cout << ".";
    }
  }
  if (verbose) cout << endl << "Forming dot products" << endl;

  // form basis dot products
  for (int m=0; m<numb; m++) {
    temp = datavol * basisx[ix[m]];
    temp *= basisy[iy[m]];
    temp *= basisz[iz[m]];
    temp *= validmask;
    IB(m+1) = temp.sum();
    for (int m2=m; m2<numb; m2++) {
      temp = basisx[ix[m]] * basisx[ix[m2]];
      temp *= basisy[iy[m]];
      temp *= basisy[iy[m2]];
      temp *= basisz[iz[m]];
      temp *= basisz[iz[m2]];
      temp *= validmask;
      BB(m+1,m2+1) = temp.sum();
      BB(m2+1,m+1) = temp.sum();
    }
    if (verbose) cout << ".";
  }
  if (verbose) cout << endl;
  
  // find least-squares solution
  if (verbose) cout << "Finding least-squares solution" << endl;
  ColumnVector C;
  C = extrapolate_coefficients(BB,IB);
  
  // generate the resultant extrapolation
  if (verbose) cout << "Regenerating result" << endl;
  volume<float> resvol;
  resvol = datavol * 0.0f;
  for (int m=0; m<numb; m++) {
    resvol += basisx[ix[m]] * basisy[iy[m]] * basisz[iz[m]] * C(m+1);
  }
//    resvol -= resvol * validmask;
//    resvol += datavol * validmask;

  return resvol;
}

///////////////////////////////////////////////////////////////////////////////

float wrap(float theta)
{  // returns wrapped phase between -pi and +pi
  const float oneontwopi = 1.0/(2.0*M_PI);
  const float twopi = 2.0*M_PI;
  return (theta - MISCMATHS::round(theta*oneontwopi)*twopi);
}


void restore_linear_ramps(const ColumnVector& ramp,
			  volume<float>& uph, const volume<float>& mask)
{
  float xav, yav, zav;
  xav = (uph.maxx() + uph.minx())/2.0;
  yav = (uph.maxy() + uph.miny())/2.0;
  zav = (uph.maxz() + uph.minz())/2.0;
  for (int z=uph.minz(); z<uph.maxz(); z++) {
    for (int y=uph.miny(); y<uph.maxy(); y++) {
      for (int x=uph.minx(); x<uph.maxx(); x++) {
	if (mask(x,y,z)>0.5) {
	  uph(x,y,z) += ramp(1)*(x-xav) + ramp(2)*(y-yav) + ramp(3)*(z-zav);
	}	
      }
    }
  }
}

void remove_linear_ramps(const ColumnVector& ramp,
			 volume<float>& ph, const volume<float>& mask)
{
  float xav, yav, zav;
  xav = (ph.maxx() + ph.minx())/2.0;
  yav = (ph.maxy() + ph.miny())/2.0;
  zav = (ph.maxz() + ph.minz())/2.0;
  for (int z=ph.minz(); z<ph.maxz(); z++) {
    for (int y=ph.miny(); y<ph.maxy(); y++) {
      for (int x=ph.minx(); x<ph.maxx(); x++) {
	if (mask(x,y,z)>0.5) {
	  ph(x,y,z) = wrap(ph(x,y,z) - ramp(1)*(x-xav) - ramp(2)*(y-yav) - ramp(3)*(z-zav));
	}	
      }
    }
  }
}

ColumnVector estimate_linear_ramps(volume<float>& ph, const volume<float>& mask)
{
  ColumnVector ramp(3);
  ramp = 0.0;
  long int nx=0, ny=0, nz=0;
  for (int z=ph.minz(); z<ph.maxz(); z++) {
    for (int y=ph.miny(); y<ph.maxy(); y++) {
      for (int x=ph.minx(); x<ph.maxx(); x++) {
	if (mask(x,y,z)>0.5) {
	  if (ph.in_bounds(x+1,y,z) && (mask(x+1,y,z)>0.5)) {
	    ramp(1) += wrap(ph(x+1,y,z)-ph(x,y,z));
	    nx++;
	  }
	  if (ph.in_bounds(x,y+1,z) && (mask(x,y+1,z)>0.5)) {
	    ramp(2) += wrap(ph(x,y+1,z)-ph(x,y,z));
	    ny++;
	  }
	  if (ph.in_bounds(x,y,z+1) && (mask(x,y,z+1)>0.5)) {
	    ramp(3) += wrap(ph(x,y,z+1)-ph(x,y,z));
	    nz++;
	  }
	}
      }
    }
  }
  if (nx>0) ramp(1) /= (float) nx;
  if (ny>0) ramp(2) /= (float) ny;
  if (nz>0) ramp(3) /= (float) nz;
  return ramp;
}


///////////////////////////////////////////////////////////////////////////////

struct sconstraint { float C; float N; float P; float K; };
typedef sconstraint constraint;
typedef map<pair<int, int> , constraint > constraintmap;

inline float getP(constraintmap::iterator& it, int r, int s)
 { if (r<s) return (*it).second.P; else return -(*it).second.P; }

inline void setP(constraintmap::iterator& it, int r, int s, float Pval)
 { if (r<s) (*it).second.P = Pval; else (*it).second.P = -Pval; }

constraintmap::iterator find(constraintmap& conmap, int r, int s)
{
  return conmap.find(pair<int,int>(Min(r,s),Max(r,s)));
}

constraintmap::iterator find_maxC(constraintmap& conmap) 
{
  if (conmap.empty()) return conmap.begin();
  constraintmap::iterator itm = conmap.begin(), it;
  float MaxC = (*itm).second.C;
  for (it = conmap.begin(); it!=conmap.end(); ++it) {
    if ((*it).second.C > MaxC) {
      MaxC = (*it).second.C;
      itm = it;
    }
  }
  return itm;
}

void calc_constraint(constraint& con)
{
  float N, P, K;
  N = (float) con.N;
  P = con.P;
  if (N>0) {
    K = -P/(2.0*M_PI*N);
  } else {
    K = 0.0;
  }
  con.K = K;
  con.C = N*(0.5 - fabs(MISCMATHS::round(K) - K));
}

void make_constraints(const volume<float>& phasemap, const volume<int>& label, 
		      constraintmap& conmap)
{
  label.setextrapolationmethod(zeropad);
  int l1=0,l2=0;
  float p1=0.0, p2=0.0;
  constraintmap::iterator it;
  for (int z=label.minz(); z<=label.maxz(); z++) {
    for (int y=label.miny(); y<=label.maxy(); y++) {
      for (int x=label.minx(); x<=label.maxx(); x++) {
	l1 = label(x,y,z);
	p1 = phasemap(x,y,z);
	if (l1!=0) {
	  for (int dir=1; dir<=3; dir++) {
	    if (dir==1) { l2=label(x-1,y,z); p2=phasemap(x-1,y,z); }
	    if (dir==2) { l2=label(x,y-1,z); p2=phasemap(x,y-1,z); }
	    if (dir==3) { l2=label(x,y,z-1); p2=phasemap(x,y,z-1); }
	    if ( (l2!=0) && (l1!=l2) ) {
	      int v1=l1, v2=l2;
	      float pp1=p1, pp2=p2;
	      if (v2<v1) { swap(v1,v2); swap(pp1,pp2); }
	      it = conmap.find(pair<int,int>(v1,v2));
	      if (it != conmap.end()) { 
		(*it).second.N++;
		(*it).second.P += pp1-pp2;
	      } else { 
		// insert this new value
		constraint newcon;
	        newcon.C=0;
		newcon.K=0;
		newcon.N = 1;
		newcon.P = pp1-pp2;
		conmap[pair<int,int>(v1,v2)] = newcon;
	      }
	    }
	  }
	}
      }
    }
  }

  // form K which is the matrix of local offset estimates
  // also form C which are the delta cost estimates
  for (it = conmap.begin(); it!=conmap.end(); ++it) {
    calc_constraint((*it).second);
  }  
  return;
}


////////////////////////////////////////////////////////////////////////////

volume<float> apply_unwrapping(const volume<float>& phasemap, 
			       const volume<int>& label, 
			       const Matrix& label_offsets)
{
  volume<float> uphase(phasemap);
  if (label_offsets.Nrows() < label.max()) {
    return uphase;
  }
  int lab;
  for (int z=label.minz(); z<=label.maxz(); z++) {
    for (int y=label.miny(); y<=label.maxy(); y++) {
      for (int x=label.minx(); x<=label.maxx(); x++) {
	lab = label(x,y,z);
	if (lab>0) {
	  uphase(x,y,z) += 2.0*M_PI*label_offsets(lab,1);
	} else {
	  uphase(x,y,z) = 0.0;
	}
      }
    }
  }
  return uphase;
}


volume<int> find_phase_labels(const volume<float>& phasemap, 
			      const volume<float>& mask,
			      int n) 
{
  // divides the phasemap into n regions and labels the connected
  //  components in each
  // Assumes that the max/min of phase is +/- PI

  volume<int> label;
  float binsize = 2.0*M_PI/((float) n);

  volume<float> phase = phasemap;
  for (int z=phase.minz(); z<=phase.maxz(); z++) {
    for (int y=phase.miny(); y<=phase.maxy(); y++) {
      for (int x=phase.minx(); x<=phase.maxx(); x++) {
	if (mask(x,y,z)>0.5) {
	  int binno = (int) ((phase(x,y,z) + M_PI)/binsize) + 1;
	  if (binno<1) binno=1;
	  if (binno>n) binno=n;
	  phase(x,y,z) = (float) binno;
	} else {
	  phase(x,y,z) = 0.0;
	}
      }
    }
  }

  label = connected_components(phase,6);
  return label;
}


volume<int> find_phase_labels2D(const volume<float>& phasemap, 
				const volume<float>& mask,
				int n, bool unique3Dlabels) 
{
  // divides the phasemap into n regions and labels the connected
  //  components in each
  // ONLY USES 2D CONNECTIVITY
  // Assumes that the max/min of phase is +/- PI

  phasemap.activateROI();
  mask.activateROI();

  volume<int> label;
  copyconvert(phasemap,label);
  label = 0;
  label.activateROI();
  int max_label = 0;

  int zmin = phasemap.minz(), zmax = phasemap.maxz();
  for (int z=zmin; z<=zmax; z++) {
    volume<float> pslice, mslice;
    volume<int> lslice, lsmask;
    // get slices from phasemap and mask
    phasemap.setROIlimits(phasemap.minx(),phasemap.miny(),z,
			  phasemap.maxx(),phasemap.maxy(),z);
    pslice = phasemap.ROI();
    mask.setROIlimits(mask.minx(),mask.miny(),z,
		      mask.maxx(),mask.maxy(),z);
    mslice = mask.ROI();
    // get phase labels for this slice
    lslice = find_phase_labels(pslice,mslice,n);
    if (unique3Dlabels) {
      // offset labels so they start at one above the previous slice
      lsmask = binarise(lslice,1);
      lslice += max_label;
      lslice *= lsmask;  // reset zero to zero
      max_label = lslice.max();
    }
    // set the corresponding slice of label to lslice
    label.setROIlimits(label.minx(),label.miny(),z,
		       label.maxx(),label.maxy(),z);
    label.copyROIonly(lslice);
  }

  // restore previous ROI status...
  phasemap.deactivateROI();
  mask.deactivateROI();
  label.deactivateROI();

  return label;
}



volume<float> unwrap(const volume<float>& phasemap, 
		     const volume<int>& label, bool verbose)
{

  //////////////////////////////////////
  // New consistent constraint method //
  //////////////////////////////////////

  if (verbose) cout << "Calculating starting matrices ("<<label.max()<<" by " 
       <<label.max()<<")"<<endl;

  constraintmap conmap;
  make_constraints(phasemap, label, conmap);
  if (verbose) cout << "Finished connection_matrices" << endl;

  int nclass = label.max();

  // find the area with the biggest interface 
  vector<float> Q(nclass+1,0.0);
  for (constraintmap::iterator it=conmap.begin(); it!=conmap.end(); ++it)
    {
      Q[(*it).first.first] += (*it).second.N;
      Q[(*it).first.second] += (*it).second.N;
    }
  int bigclassnum = 1;
  for (int n=1; n<=nclass; n++) {
    if (Q[bigclassnum] < Q[n])  bigclassnum = n;
  }

  vector<int> Mvec(nclass+1), Dvec(nclass+1,0);
  for (int k=0; k<=nclass; k++)  Mvec[k] = k;

  int iterations=0;
  while (!conmap.empty()) {  // should execute nclass times
    iterations++;
    if ((iterations % 1000) == 1)
      if (verbose) cout << conmap.size() << " constraints left" << endl;
    constraintmap::iterator it;
    it = find_maxC(conmap);
    int r,s;
    r = (*it).first.first;
    s = (*it).first.second;
    int Krs = MISCMATHS::round((*it).second.K);
    conmap.erase(it);
    // merge s into r (note that r<s due to insertion into constraintmap)
    for (int p=1; p<=nclass; p++) {
      if (Mvec[p]==s) {
	Mvec[p] = r;
	Dvec[p] -= Krs;
      }
    }
    for (int t=1; t<=nclass; t++) {
      constraintmap::iterator itst = find(conmap,s,t);
      if (itst != conmap.end()) {
	constraintmap::iterator itrt = find(conmap,r,t);
	float Pst = getP(itst,s,t);
	float Nst = (*itst).second.N;
	float Prt=0.0, Nrt=0.0;
	if (itrt != conmap.end()) {
	  if (itrt == itst) { cerr << "ERROR!  itrt==itst" << endl; }
	  Prt = getP(itrt,r,t);
	  Nrt = (*itrt).second.N;
	  conmap.erase(itrt);
	}
	conmap.erase(itst);
	constraint newcon;
	newcon.N = Nrt + Nst;
	float Prt2 = Prt + Pst - 2*M_PI*Nst*Krs;
	if (r<t) newcon.P = Prt2;
	else     newcon.P = -Prt2;
	calc_constraint(newcon);
	conmap[pair<int,int>(Min(r,t),Max(r,t))] = newcon;
      }
    }
  }
  if (verbose) cout << "Did while loop " << iterations << " times" << endl;
  
  Matrix Mbest(nclass,1);
  for (int n=1; n<=nclass; n++)  
      Mbest(n,1) = (float) Dvec[n] - Dvec[bigclassnum];

  // now unwrap the phase
  volume<float> uphase;
  uphase = apply_unwrapping(phasemap,label,Mbest);

  return uphase;
}



volume<float> unwrap2D(const volume<float>& phasemap, 
		       const volume<int>& label, bool verbose)
{
  // ONLY USES 2D CONNECTIVITY

  int zmin = phasemap.minz(), zmax = phasemap.maxz();
  phasemap.activateROI();
  label.activateROI();

  volume<float> uphase;
  uphase = phasemap;
  uphase = 0;
  uphase.activateROI();

  volume<float> pslice, uslice, udiff, mslice;
  volume<int> lslice;

  for (int z=zmin; z<=zmax; z++) {
    if (verbose) cout << "SLICE NUMBER " << z << endl;

    // get slices from phasemap and label
    phasemap.setROIlimits(phasemap.minx(),phasemap.miny(),z,
			  phasemap.maxx(),phasemap.maxy(),z);
    pslice = phasemap.ROI();
    label.setROIlimits(label.minx(),label.miny(),z,
		       label.maxx(),label.maxy(),z);
    lslice = label.ROI();
    copyconvert(lslice,mslice);
    mslice.binarise(0.5);
    // unwrap 
    uslice = unwrap(pslice,lslice,verbose);
    // set the corresponding slice of uphase to uslice
    uphase.setROIlimits(uphase.minx(),uphase.miny(),z,
			uphase.maxx(),uphase.maxy(),z);
    uphase.copyROIonly(uslice);
  }
  
  // restore previous ROI status...
  phasemap.deactivateROI();
  label.deactivateROI();
  uphase.deactivateROI();
  
// now fix slice-offsets using median difference

// //  NB: used to try the same unwrap() call with a single label per slice but it did not work well!
//   volume<int> slicelabels(label);
//   int labelnum=0;
//   for (int z=zmin; z<=zmax; z++) {
//     bool empty=true;
//     for (int y=slicelabels.miny(); y<=slicelabels.maxy(); y++) {
//       for (int x=slicelabels.minx(); x<=slicelabels.maxx(); x++) {
// 	if (slicelabels(x,y,z)>0) {
// 	  if (empty) { empty=false; labelnum++; }
// 	  slicelabels(x,y,z)=labelnum;
// 	}
//       }
//     }
//   }
//   save_volume(slicelabels,"slicelabels");
//   unwrap(uphase,slicelabels,true);

  // instead, use a median difference (best when combined with --removeramps)
  vector<float> diffvec;
  for (int z=zmin; z<zmax; z++) {
    float phasediff=0.0;
    diffvec.clear();
    for (int y=label.miny(); y<=label.maxy(); y++) {
      for (int x=label.minx(); x<=label.maxx(); x++) {
	if ((label(x,y,z)>0) && (label(x,y,z+1)>0))  {
	  phasediff = uphase(x,y,z+1)-uphase(x,y,z);
	  diffvec.push_back(phasediff);
	}
      }
    }
    if (diffvec.size()>0) {
      sort(diffvec.begin(),diffvec.end());
      float median_diff=diffvec[MISCMATHS::round((int) diffvec.size()/2)];
      int m=MISCMATHS::round(median_diff/(2.0*M_PI));
      for (int y=label.miny(); y<=label.maxy(); y++) {
	for (int x=label.minx(); x<=label.maxx(); x++) {
	  if (label(x,y,z+1)>0) { uphase(x,y,z+1) -= m*2.0*M_PI; }
	}
      }
    }
  }

  return uphase;
}


////////////////////////////////////////////////////////////////////////////

volume<float> calc_fmap(const volume<float>& phase1, 
			const volume<float>& phase2,
			const volume<float>& mask1,
			const volume<float>& mask2,
			float asym_time)
{
  volume<float> fmap, fmask;
  fmap = phase2 - phase1;
  fmask = mask1 * mask2;
  fmap *= fmask;
  // check to see if there is an overall phase offset, due to the
  //  arbitrary phase chosen in the phase maps
  float sum, n;
  sum = fmap.sum();
  n = fmask.sum();
  int offset = MISCMATHS::round(sum/(n*2.0*M_PI));
  if (offset!=0) {
    fmap -= (float) (2.0*M_PI*offset) * fmask;
  }
  float phase_to_rads_factor = 1.0/asym_time;
  fmap *= phase_to_rads_factor;
  return fmap;
}


float fmap2pixshift_factor(const volume<float>& invol, float pe_dwell_time,
			   const string& dir)
{
  int n_pe=64;
  if ((dir=="x") || (dir=="x-")) { n_pe = invol.xsize(); }
  if ((dir=="y") || (dir=="y-")) { n_pe = invol.ysize(); }
  if ((dir=="z") || (dir=="z-")) { n_pe = invol.zsize(); }
  return (n_pe*pe_dwell_time)/(2.0*M_PI);
}


volume<float> yderiv(const volume<float>& vol) 
{
  vol.setextrapolationmethod(zeropad);
  ColumnVector derivmask(3), unitmask(1);
  unitmask = 1.0;
  derivmask << -0.5 << 0.0 << 0.5;
  return convolve_separable(vol,unitmask,derivmask,unitmask);
}


////////////////////////////////////////////////////////////////////////////

volume<float> apply_pixshift(const volume<float>& epivol,
			     const volume<float>& pixshiftmap)
{
  volume<float> outvol(epivol);
  for (int z=0; z<outvol.zsize(); z++) {
    for (int y=0; y<outvol.ysize(); y++) {
      for (int x=0; x<outvol.xsize(); x++) {
	outvol(x,y,z) = epivol.interpolate((float) x, 
					   ((float) y)+pixshiftmap(x,y,z), 
					   (float) z);
      }
    }
  }
  return outvol;
}


volume<float> limit_pixshift(const volume<float>& pixshift, 
			     const volume<float>& mask, float mindiff)
{
  volume<float> outvol(pixshift);
  for (int z=0; z<pixshift.zsize(); z++) {
    for (int x=0; x<pixshift.xsize(); x++) {

      // find COM for this line
      int ycom = pixshift.ysize()/2, n=0;
      float sum=0.0;
      for (int y0=0; y0<pixshift.ysize(); y0++) {
	if (mask(x,y0,z)>0.5) {
	  n++;  sum+=y0;
	}
      }
      if (n>2) {
	ycom = MISCMATHS::round(sum/((float) n));
      }
	
      for (int y=ycom; y<pixshift.ysize()-1; y++) {
	outvol(x,y+1,z) = Max(pixshift(x,y+1,z),outvol(x,y,z) + mindiff);
      }
      for (int y=ycom; y>=1; y--) {
	outvol(x,y-1,z) = Min(pixshift(x,y-1,z),outvol(x,y,z) - mindiff);
      }
    }
  }
  return outvol;
} 

////////////////////////////////////////////////////////////////////////////

