/*  costfns.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

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

// Interpolation functions
//  Written by Mark Jenkinson  18/2/99

#include <iostream>
#include <cassert>
#include "costfns.h"
#include "newmatio.h"
#include "miscmaths/miscmaths.h"
#include "newimage.h"
#include "newimagefns.h"

using namespace MISCMATHS;

namespace NEWIMAGE {

// make SAFE_FLIRT the default now (bounds checking on all interpolations)
#define SAFE_FLIRT 1

  costfns costfn_type(const string& cname) 
  {
    costfns cfn=Unknown;
    if (cname == "mutualinfo") {
      cfn = MutualInfo;
    } else if (cname == "corratio") {
      cfn = CorrRatio;
    } else if (cname == "woods") {
      cfn = Woods;
    } else if (cname == "normcorr") {
      cfn = NormCorr;
    } else if (cname == "normmi") {
      cfn = NormMI;
    } else if (cname == "leastsq") {
      cfn = LeastSq;
    } else if (cname == "labeldiff") {
      cfn = LabelDiff;
    } else if (cname == "bbr") {
      cfn = BBR;
    }
    return cfn;
  }

   /////////////////////////////////////////////////////////////////////////

   // quick, efficient interpolation calls

  inline bool in_interp_bounds(const volume<float>& v,
			       const float x, const float y, const float z)
  {
      int ix, iy, iz;
      ix=(int) x; iy=(int) y; iz=(int) z;
      if ( (!v.in_bounds(ix,iy,iz)) || (!v.in_bounds(ix+1,iy+1,iz+1)) )
	 return false;
      else
	return true;
  }


  inline float q_tri_interpolation(const volume<float>& v, 
				   const float x, const float y, const float z)
    {
      int ix, iy, iz;
      ix=(int) x; iy=(int) y; iz=(int) z;
      float dx=x-ix, dy=y-iy, dz=z-iz;
      float v000, v001, v010, v011, v100, v101, v110, v111;
      float retval=0;
#ifdef SAFE_FLIRT
      if ( (ix>=0) && (iy>=0) && (iz>=0) 
	   && (ix<v.maxx()) && (iy<v.maxy()) && (iz<v.maxz()) ) {
#endif
	v.getneighbours(ix,iy,iz, v000,v001,v010,v011,v100,v101,v110,v111);
	{
	  float temp1, temp2, temp3, temp4, temp5, temp6;
	  temp1 = (v100 - v000)*dx + v000;
	  temp2 = (v101 - v001)*dx + v001;
	  temp3 = (v110 - v010)*dx + v010;
	  temp4 = (v111 - v011)*dx + v011;
	  // second order terms
	  temp5 = (temp3 - temp1)*dy + temp1;
	  temp6 = (temp4 - temp2)*dy + temp2;
	  // final third order term
	  retval = (temp6 - temp5)*dz + temp5;
	}
#ifdef SAFE_FLIRT
      } else {
	retval=v.getpadvalue();
      }
#endif
      return retval;
    }


  int q_get_neighbours(const volume<float>& v, 
		       const float x, const float y, const float z,
		       float& v000, float& v001, float& v010, float& v011,
		       float& v100, float& v101, float& v110, float& v111,
		       float& dx, float& dy, float& dz)
  {
    int ix, iy, iz;
    ix=(int) x; iy=(int) y; iz=(int) z;
    dx=x-ix; dy=y-iy; dz=z-iz;
#ifdef SAFE_FLIRT
    if ( (ix>=0) && (iy>=0) && (iz>=0) 
	 && (ix<v.maxx()) && (iy<v.maxy()) && (iz<v.maxz()) ) {
#endif
      v.getneighbours(ix,iy,iz, v000,v001,v010,v011,v100,v101,v110,v111);
#ifdef SAFE_FLIRT
    } else {
      v000=v001=v010=v011=v100=v101=v110=v111=v.getpadvalue();
    }
#endif
    return 0;
  }


   /////////////////////////////////////////////////////////////////////////

   float q_sinc(float x)
   {
     if (fabs(x)<1e-7) { return 1.0-fabs(x); }
     float y=M_PI*x;
     return sin(y)/y;
   }
   
   float q_hanning(float x, int w)
   {
     if (fabs(x)>w) 
       return 0.0;
     else
       return (0.5 + 0.5 *cos(M_PI*x/w));
   }

   // somewhat cheating global statics
   static int Globalkernelwidth;
   static float Globalsinckernel[201];
   
   void q_setupkernel()
   {
     Globalkernelwidth=3;
     // set x between +/- kernelwidth
     for (int n=0; n<=200; n++) {
       float x=(n-100)/100.0*Globalkernelwidth;
       Globalsinckernel[n] = q_sinc(x)*q_hanning(x,Globalkernelwidth);
     }
   }
   
   float q_kernelval(float x, int w)
   {
     // effectively returns  sinc(x)*hanning(x,w);
     if (fabs(x)>w) return 0.0;
     float dn = x/w*100.0 + 100;
     int n = (int) floor(dn);
     dn -= n;
     if (n>=200) return 0.0;
     if (n<0) return 0.0;
     
     return Globalsinckernel[n]*(1.0-dn) + Globalsinckernel[n+1]*dn;
   }
   
   
   float q_sinc_interpolation(const volume<float>& v, 
			      const float x, const float y, const float z)
   {
     // kernel half-width  (i.e. range is +/- w)
     int w=Globalkernelwidth;  
     if (w<1) { 
       q_setupkernel(); 
       w=Globalkernelwidth;
     }
     
     int ix0, iy0, iz0;
     ix0 = (int) floor(x);
     iy0 = (int) floor(y);
     iz0 = (int) floor(z);
     
     float convsum=0.0, interpval=0.0, kersum=0.0;
     static float sincz[201], sincy[201], sincx[201];  // limits width to 200
     
     for (int d=-w; d<=w; d++) {
       sincz[d+w] = q_kernelval((z-iz0+d),w);
       sincy[d+w] = q_kernelval((y-iy0+d),w);
       sincx[d+w] = q_kernelval((x-ix0+d),w);
     }
     
     int xj, yj, zj;
     int x1a, x1b, y1a, y1b, z1a, z1b;
     x1a = Max(0,ix0-w);  x1b=Min(v.xsize()-1,ix0+w);
     y1a = Max(0,iy0-w);  y1b=Min(v.ysize()-1,iy0+w);
     z1a = Max(0,iz0-w);  z1b=Min(v.zsize()-1,iz0+w);
     for (int z1=z1a; z1<=z1b; z1++) {
       zj=iz0-z1+w;
       for (int y1=y1a; y1<=y1b; y1++) {
	 yj=iy0-y1+w;
	 for (int x1=x1a; x1<=x1b; x1++) {
	     xj=ix0-x1+w;
	     float sincfac = sincx[xj] * sincy[yj] * sincz[zj];
	     convsum += v.value(x1,y1,z1) * sincfac;
	     kersum += sincfac;
	 }
       }
     }
     
     if (fabs(kersum)>1e-9) {
       interpval = convsum / kersum;
     } else {
       return v.backgroundval();
     }
     return interpval;
   }
   
   /////////////////////////////////////////////////////////////////////////


   Costfn::Costfn(const volume<float>& refv, const volume<float>& inputv) :
     refvol(refv), testvol(inputv), rweight(refv), tweight(inputv), 
     wmseg(refv), fmap(), fmap_mask(), debugvol(),
     bindex(0), no_bins(0),jointhist(0), marghist1(0), marghist2(0), 
     fjointhist(0), fmarghist1(0), fmarghist2(0), p_count(0), 
     p_costtype(CorrRatio), validweights(false), bin_a0(0), bin_a1(1),
     bbr_dist(2.0), bbr_offset(0.0), bbr_slope(-0.5), gm_coord_x(0), 
     gm_coord_y(0), gm_coord_z(0), wm_coord_x(0), wm_coord_y(0), wm_coord_z(0),
     no_coords(0), vertex_step(1), pe_dir(0), bbr_type("signed"), debug_mode(false), 
     smoothsize(1.0), fuzzyfrac(0.5)
   { 
   }

   Costfn::Costfn(const volume<float>& refv, const volume<float>& inputv,
		  const volume<float>& refweight, 
		  const volume<float>& inweight) :
     refvol(refv), testvol(inputv), rweight(refweight), tweight(inweight), 
     wmseg(refv), fmap(), fmap_mask(), debugvol(),
     bindex(0), no_bins(0),jointhist(0), marghist1(0), marghist2(0), 
     fjointhist(0), fmarghist1(0), fmarghist2(0), p_count(0), 
     p_costtype(CorrRatio), validweights(true), bin_a0(0), bin_a1(1),
     bbr_dist(2.0), bbr_offset(0.0), bbr_slope(-0.5), gm_coord_x(0), 
     gm_coord_y(0), gm_coord_z(0), wm_coord_x(0), wm_coord_y(0), wm_coord_z(0),
     no_coords(0), vertex_step(1), pe_dir(0), bbr_type("signed"), debug_mode(false), 
     smoothsize(1.0), fuzzyfrac(0.5)
   { 
   }

    Costfn::~Costfn() { 
     if (jointhist)  delete [] jointhist; 
     if (marghist1)  delete [] marghist1; 
     if (marghist2)  delete [] marghist2; 
     if (fjointhist)  delete [] fjointhist; 
     if (fmarghist1)  delete [] fmarghist1; 
     if (fmarghist2)  delete [] fmarghist2; 
     if (bindex)     delete [] bindex;
     if (gm_coord_x) delete [] gm_coord_x;
     if (gm_coord_y) delete [] gm_coord_y;
     if (gm_coord_z) delete [] gm_coord_z;
     if (wm_coord_x) delete [] wm_coord_x;
     if (wm_coord_y) delete [] wm_coord_y;
     if (wm_coord_z) delete [] wm_coord_z;
   }

   // General cost function calls


//   // affmat is voxel to voxel and non-linear parameters are arbitrary
//   float Costfn::cost(const Matrix& affmat, const ColumnVector& nonlin_params) const
//   {
//   }

   float Costfn::cost(const Matrix& affmat) const // affmat is voxel to voxel
    {
      if (validweights) {
	return this->cost(affmat,rweight,tweight);
      }

      float retval = 0.0;
      switch (p_costtype)
	{
	case NormCorr:  // MAXimise corr
	  if (smoothsize > 0.0) { 
	    retval = 1.0 - 
	      fabs(this->normcorr_smoothed(affmat));
	  } else {
	    retval = 1.0 - fabs(this->normcorr(affmat));
	  }
	  break;
	case NormCorrSinc:  // MAXimise corr
	  retval = 1.0 - fabs(this->normcorr_smoothed_sinc(affmat));
	  break;
	case LeastSq:  // Minimise square
	  if (smoothsize > 0.0) {
	    retval = this->leastsquares_smoothed(affmat);
	  } else {
	    retval = this->leastsquares(affmat);
	  }
	  break;
	case LabelDiff:  // Minimise label difference (any mismatch is equal)
	  if (smoothsize > 0.0) {
	    retval = this->labeldiff_smoothed(affmat);
	  } else {
	    retval = this->labeldiff(affmat);
	  }
	  break;
	case BBR:  // Minimise appropriate BBR metric tanh(intensity difference)
	  retval = this->bbr(affmat);
	  break;
	case CorrRatio:  // MAXimise corr
	  if (smoothsize > 0.0) {
	    retval = 1.0 - this->corr_ratio_smoothed(affmat);
	  } else {
	    retval = 1.0 - this->corr_ratio(affmat); 
	  }
	  break;
	case Woods:  // minimise variance/mean
	  retval = this->woods_fn(affmat); 
	  break;
	case MutualInfo:  // MAXimise info
	  if ((smoothsize > 0.0) || (fuzzyfrac > 0.0)) {
	    retval = -this->mutual_info_smoothed(affmat); 
	  } else {
	    retval = -this->mutual_info(affmat); 
	  }
	  break;
	case NormMI:  // MAXimise
	  if ((smoothsize > 0.0) || (fuzzyfrac > 0.0)) {
	    retval = -this->normalised_mutual_info_smoothed(affmat); 
	  } else {
	    retval = -this->normalised_mutual_info(affmat); 
	  }
	  break;
	default:
	  cerr << "Invalid cost function type" << endl;
	  break;
	}
      return retval;
    }


  float Costfn::cost(const Matrix& affmat, const ColumnVector& nonlin_params) const // affmat is voxel to voxel
    {
//       if (validweights) {
// 	return this->cost(affmat,nonlin_params,rweight,tweight);
//       }

      float retval = 0.0;
      switch (p_costtype)
	{
	case BBR:  // Minimise appropriate BBR metric tanh(intensity difference)
	  retval = this->bbr(affmat,nonlin_params);
	  break;
	default:
	  cerr << "Invalid cost function type" << endl;
	  break;
	}
      return retval;
    }


//    float Costfn::cost(const Matrix& affmat,   // affmat is voxel to voxel
// 		      const ColumnVector& nonlin_params,  // arbitrary units
// 		      const volume<float>& refweight, 
// 		      const volume<float>& testweight) const
//    {
//    }

   float Costfn::cost(const Matrix& affmat,   // affmat is voxel to voxel
		      const volume<float>& refweight, 
		      const volume<float>& testweight) const
   {
     float retval = 0.0;
     switch (p_costtype) 
       {
       case NormCorr:  // MAXimise corr
	 retval = 1.0-this->normcorr_fully_weighted(affmat,refweight,
						    testweight);
	 break;
	case NormCorrSinc:  // MAXimise corr
	 cerr<<"WARNING: NormCorrSinc is not implemented with cost function weighting"
	     << endl;
	  retval = 1.0 - fabs(this->normcorr_smoothed_sinc(affmat));
	  break;
       case LeastSq:  // Minimise square
	 retval = this->leastsquares_fully_weighted(affmat,refweight,
							testweight);
	 break;
       case LabelDiff:  // Minimise label difference (any mismatch is equal)
	 retval = this->labeldiff_fully_weighted(affmat,refweight,
							testweight);
	 break;
       case BBR:  // Minimise appropriate BBR metric tanh(intensity difference)
	 retval = this->bbr(affmat);
	 break;
       case CorrRatio:  // MAXimise corr
	 retval = 1.0-this->corr_ratio_fully_weighted(affmat,refweight,
						      testweight);
	 break;
       case Woods:  // minimise variance/mean
	 cerr<<"WARNING: Woods is not implemented with cost function weighting"
	     << endl;
	 retval = this->woods_fn(affmat); 
	 break;
       case MutualInfo:  // MAXimise info
	retval = -this->mutual_info_fully_weighted(affmat,refweight,
						      testweight);
	break;
       case NormMI:  // MAXimise
	retval = -this->normalised_mutual_info_fully_weighted(affmat,
								 refweight,
								 testweight);
	break;
       default:
	 cerr << "Invalid cost function type" << endl;
	 break;
       }
     return retval;
   }



  float Costfn::cost(const volume4D<float>& warp,  // warp is mm to mm
		     const volume<float>& refweight, 
		     const volume<float>& testweight) const

    {
      float retval = 0.0;
      switch (p_costtype)
	{
	case CorrRatio:  // MAXimise corr
	  retval = 1.0 - this->corr_ratio_fully_weighted(warp, refweight, 
							 testweight);
	    break;
	default:
	  cerr << "Invalid cost function type" << endl;
	  break;
	}
      return retval;
    }




   float Costfn::cost(const volume4D<float>& warp) const // warp is mm to mm
    {
      if (validweights) {
	return this->cost(warp,rweight,tweight);
      }

      float retval = 0.0;
      switch (p_costtype)
	{
	case CorrRatio:  // MAXimise corr
	  cerr << "Non-weighted Correlation Ratio not yet available" << endl;
	  break;
	default:
	  cerr << "Invalid cost function type" << endl;
	  break;
	}
      return retval;
    }


 
  float Costfn::cost_gradient(volume4D<float>& gradvec,
			      const volume4D<float>& warp,  // warp is mm to mm
			      const volume<float>& refweight, 
			      const volume<float>& testweight, 
			      bool nullbc) const

  {
      float retval = 0.0;
      switch (p_costtype)
	{
	case CorrRatio:  // MAXimise corr
	  retval = 1.0 - this->corr_ratio_gradient_fully_weighted(gradvec,
								  warp, 
								  refweight, 
								  testweight,
								  nullbc);
	  gradvec *= -1.0f;
	  break;
	default:
	  cerr << "Invalid cost function type" << endl;
	  break;
	}
      return retval;
    }

  
  float Costfn::cost_gradient(volume4D<float>& gradvec, // warp is mm to mm
			      const volume4D<float>& warp, bool nullbc) const
   
    {
      if (validweights) {
	return this->cost_gradient(gradvec,warp,rweight,tweight,nullbc);
      }
      float retval = 0.0;
      switch (p_costtype)
	{
	case CorrRatio:  // MAXimise corr
	  cerr << "Non-weighted Correlation Ratio not yet available" << endl;
	  break;
	default:
	  cerr << "Invalid cost function type" << endl;
	  break;
	}
      return retval;
    }

  

   // Helper functions
   inline int *get_bindexptr(unsigned int x, unsigned int y, unsigned int z,
			     const volume<float>& refvol, int *bindex) 
   { return bindex + ( z * refvol.ysize() + y ) * refvol.xsize() + x; }
   

   void findrangex(unsigned int &xmin1 , unsigned int &xmax1,
		   float o1, float o2, float o3,
		   float a11, float a21, float a31,
		   unsigned int xb1, unsigned int yb1, unsigned int zb1,
		   float xb2, float yb2, float zb2) {
     
     float x1, x2, xmin, xmax, xmin0, xmax0;
     
     xmin0 = 0;
     xmax0 = xb1;
      
     if (fabs(a11)<1.0e-8) {
       if ((0.0<=o1) && (o1<=xb2)) {
	 x1 = -1.0e8; x2 = 1.0e8;
       } else {
	 x1 = -1.0e8; x2 = -1.0e8;
       }
     } else {
       x1 = -o1/a11;
       x2 = (xb2-o1)/a11;
     }
     xmin = Min(x1,x2);
     xmax = Max(x1,x2);
     // intersect ranges
     xmin0 = Max(xmin0,xmin);
     xmax0 = Min(xmax0,xmax);
	  
     if (fabs(a21)<1.0e-8) {
       if ((0.0<=o2) && (o2<=yb2)) {
	 x1 = -1.0e8; x2 = 1.0e8;
       } else {
	 x1 = -1.0e8; x2 = -1.0e8;
       }
     } else {
       x1 = -o2/a21;
       x2 = (yb2-o2)/a21;
     }
     xmin = Min(x1,x2);
     xmax = Max(x1,x2);
     // intersect ranges
     xmin0 = Max(xmin0,xmin);
     xmax0 = Min(xmax0,xmax);

     if (fabs(a31)<1.0e-8) {
       if ((0.0<=o3) && (o3<=zb2)) {
	 x1 = -1.0e8; x2 = 1.0e8;
       } else {
	 x1 = -1.0e8; x2 = -1.0e8;
       }
     } else {
       x1 = -o3/a31;
       x2 = (zb2-o3)/a31;
     }
     xmin = Min(x1,x2);
     xmax = Max(x1,x2);
     // intersect ranges
     xmin0 = Max(xmin0,xmin);
     xmax0 = Min(xmax0,xmax);
    
     //assert(xmin0>=0.0);
     //assert(xmax0<=xb1);

     if (xmax0<xmin0) {
       xmax1=0;
       xmin1=1;
     } else {
       xmin1 = (unsigned int) ceil(xmin0);
       xmax1 = (unsigned int) floor(xmax0);
     }

     // brute force check to see that the floating point accumulation will
     //  not go wrong later
     unsigned int xminT = xmin1, xmaxT = xmax1;
     float x=o1, y=o2, z=o3;
     x += xminT * a11;
     y += xminT * a21;
     z += xminT * a31;
     for (unsigned int xc=xminT; xc<=xmaxT; xc++) {
       bool inbounds=( (x<=xb2) && (x>=0.0) && (y<=yb2) && (y>=0.0) && 
		       (z<=zb2) && (z>=0.0) );
       if ( (xc==xmin1) && (!inbounds) ) {
	 xmin1++;
       } else {
	 if (!inbounds) {
	   xmax1=xc-1;  return;
	 }
       }
       x += a11;
       y += a21;
       z += a31;
     }
     // if it gets here then there was no problem so just return xmin1,xmax1
   }

   //--------------------------------------------------------------------//

   float p_corr_ratio_fully_weighted(const volume<float>& vref, 
				   const volume<float>& vtest,
				   const volume<float>& refweight, 
				   const volume<float>& testweight,
				   int *bindex, const Matrix& aff,
				   const int no_bins, const float smoothsize)
    {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1

      Matrix iaffbig = vtest.sampling_mat().i() * aff.i() *
	                     vref.sampling_mat();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.xsize()-1, yb1=vref.ysize()-1, zb1=vref.zsize()-1;
      float  xb2 = ((float) vtest.xsize())-1.0001,
	yb2=((float) vtest.ysize())-1.0001, zb2=((float) vtest.zsize())-1.0001;

      float *sumy, *sumy2;
      sumy = new float[no_bins+1];
      sumy2 = new float[no_bins+1];
      float *numy;
      numy = new float[no_bins+1];
      int b=0;
 
      for (int i=0; i<=no_bins; i++) {
	numy[i]=0.0; sumy[i]=0.0;  sumy2[i]=0.0;
      }

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4);
      float wval=0,val=0,o1,o2,o3;

      float smoothx, smoothy, smoothz, weight;
      smoothx = smoothsize / vtest.xdim();
      smoothy = smoothsize / vtest.ydim();
      smoothz = smoothsize / vtest.zdim();

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      unsigned int xmin, xmax;
      int *bptr;

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  // assume that this is always OK
	  bptr = get_bindexptr(xmin,y,z,vref,bindex);

	  for (unsigned int x=xmin; x<=xmax; x++) {

	    if ( !((x==xmin) || (x==xmax)) 
		 || in_interp_bounds(vtest,o1,o2,o3) )
	      {
		val = q_tri_interpolation(vtest,o1,o2,o3);
		wval = q_tri_interpolation(testweight,o1,o2,o3);
		
		// do the cost function record keeping...
		b=*bptr;
		weight=wval*refweight(x,y,z);
		if (o1<smoothx)  weight*=o1/smoothx;
		else if ((xb2-o1)<smoothx) weight*=(xb2-o1)/smoothx;
		if (o2<smoothy)  weight*=o2/smoothy;
		else if ((yb2-o2)<smoothy) weight*=(yb2-o2)/smoothy;
		if (o3<smoothz)  weight*=o3/smoothz;
		else if ((zb2-o3)<smoothz) weight*=(zb2-o3)/smoothz;
		if (weight<0.0)  weight=0.0;
		numy[b]+=weight;
		sumy[b]+=weight*val;
		sumy2[b]+=weight*val*val;
	      }

	    bptr++;
	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	}
      }


      float corr_ratio=0.0, var=0.0, totsumy=0.0, totsumy2=0.0;
      float numtoty=0.0;

      // correct for occasion lapses into the last bin
      numy[no_bins-1] += numy[no_bins];
      sumy[no_bins-1] += sumy[no_bins];
      sumy2[no_bins-1] += sumy2[no_bins];
      numy[no_bins]=0.0;
      sumy[no_bins]=0.0;
      sumy2[no_bins]=0.0;

      // now calculate the individual variances for each iso-set
      //  weighting them by the number of pixels from Image x that contribute
      for (b=0; b<no_bins; b++) {
	if (numy[b]>2.0) {
	  numtoty += numy[b];
	  totsumy += sumy[b];
	  totsumy2 += sumy2[b];
	  // the following should be the variance of the bth iso-subset
	  var = (sumy2[b] - sumy[b]*sumy[b]/numy[b] ) / ( numy[b]-1);
	  // cerr << "Set #" << b << " has " << numy[b] << " elements and " 
	  //   << var << " variance" << endl;
	  corr_ratio += var * ((float) numy[b]);
	}
      }
      delete [] numy; delete [] sumy; delete [] sumy2;

      // normalise the weighting of numy[]
      if (numtoty>0)  corr_ratio/=((float) numtoty);
      // calculate the total variance of Image y and then normalise by this
      if (numtoty>1)
	var = ( totsumy2 - totsumy*totsumy/numtoty ) / (numtoty - 1);
      //cerr << "TOTALS are:" << endl 
      //   << " numerator variance is : " << corr_ratio << endl
      //   << " and denominator variance is: " << var << " from " << numtoty 
      //   << " valid elements" << endl;
      if (var>0.0)  corr_ratio/=var;
      // the above is actually 1 - correlation ratio, so correct this now
      if ( (numtoty<=1) || (var<=0.0) )
	return 0.0;   // the totally uncorrelated condition
      else
	return (1.0 - corr_ratio);

      // an alternative is to return 1.0/corr_ratio (=1/(1-correlation ratio))
      //  which may be better at rewarding gains near the best solution

      return 0;

    }

  ///////////////////////////////////////////////////////////////////////


   float p_corr_ratio_fully_weighted(const volume<float>& vref, 
				   const volume<float>& vtest,
				   const volume<float>& refweight, 
				   const volume<float>& testweight,
				   int *bindex, const volume4D<float>& warpvol,
				   const int no_bins, const float smoothsize)
   {

      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1
      // This assumes that warpvol is in the same space (and same sampling)
      //  as vref

      float  xb2 = ((float) vtest.xsize())-1.0001,
	yb2=((float) vtest.ysize())-1.0001, zb2=((float) vtest.zsize())-1.0001;

      float *sumy, *sumy2;
      sumy = new float[no_bins+1];
      sumy2 = new float[no_bins+1];
      float *numy;
      numy = new float[no_bins+1];
      int b=0;
 
      for (int i=0; i<=no_bins; i++) {
	numy[i]=0.0; sumy[i]=0.0;  sumy2[i]=0.0;
      }

      float wval=0,val=0,o1,o2,o3;

      float smoothx, smoothy, smoothz, weight;
      smoothx = smoothsize / vtest.xdim();
      smoothy = smoothsize / vtest.ydim();
      smoothz = smoothsize / vtest.zdim();

      int *bptr;

      for (int z=vref.minz(); z<=vref.maxz(); z++) {
	for (int y=vref.miny(); y<=vref.maxy(); y++) {
	  for (int x=vref.minx(); x<=vref.maxx(); x++) {
	    //   assume vref and warpvol have same voxel size
	    //   look up the warp dest coordinate (in mm)
	    if (warpvol[0].in_bounds(x,y,z)) {
	      o1 = warpvol[0](x,y,z) / vtest.xdim();
	      o2 = warpvol[1](x,y,z) / vtest.ydim();
	      o3 = warpvol[2](x,y,z) / vtest.zdim();

	      if (in_interp_bounds(vtest,o1,o2,o3)) {
		val = q_tri_interpolation(vtest,o1,o2,o3);
		wval = q_tri_interpolation(testweight,o1,o2,o3);
		
		// do the cost function record keeping...
		bptr = get_bindexptr(x,y,z,vref,bindex);
		b=*bptr;
		weight=wval*refweight(x,y,z);
		if (o1<smoothx)  weight*=o1/smoothx;
		else if ((xb2-o1)<smoothx) weight*=(xb2-o1)/smoothx;
		if (o2<smoothy)  weight*=o2/smoothy;
		else if ((yb2-o2)<smoothy) weight*=(yb2-o2)/smoothy;
		if (o3<smoothz)  weight*=o3/smoothz;
		else if ((zb2-o3)<smoothz) weight*=(zb2-o3)/smoothz;
		if (weight<0.0)  weight=0.0;
		numy[b]+=weight;
		sumy[b]+=weight*val;
		sumy2[b]+=weight*val*val;
	      }
	    }
	  }
	}
      }
      
      
      float corr_ratio=0.0, var=0.0, totsumy=0.0, totsumy2=0.0;
      float numtoty=0.0;

      // correct for occasion lapses into the last bin
      numy[no_bins-1] += numy[no_bins];
      sumy[no_bins-1] += sumy[no_bins];
      sumy2[no_bins-1] += sumy2[no_bins];
      numy[no_bins]=0.0;
      sumy[no_bins]=0.0;
      sumy2[no_bins]=0.0;

      // now calculate the individual variances for each iso-set
      //  weighting them by the number of pixels from Image x that contribute
      for (b=0; b<no_bins; b++) {
	if (numy[b]>2.0) {
	  numtoty += numy[b];
	  totsumy += sumy[b];
	  totsumy2 += sumy2[b];
	  // the following should be the variance of the bth iso-subset
	  var = (sumy2[b] - sumy[b]*sumy[b]/numy[b] ) / ( numy[b]-1);
	  // cerr << "Set #" << b << " has " << numy[b] << " elements and " 
	  //   << var << " variance" << endl;
	  corr_ratio += var * ((float) numy[b]);
	}
      }
      delete [] numy; delete [] sumy; delete [] sumy2;

      // normalise the weighting of numy[]
      if (numtoty>0)  corr_ratio/=((float) numtoty);
      // calculate the total variance of Image y and then normalise by this
      if (numtoty>1)
	var = ( totsumy2 - totsumy*totsumy/numtoty ) / (numtoty - 1);
      //cerr << "TOTALS are:" << endl 
      //   << " numerator variance is : " << corr_ratio << endl
      //   << " and denominator variance is: " << var << " from " << numtoty 
      //   << " valid elements" << endl;
      if (var>0.0)  corr_ratio/=var;
      // the above is actually 1 - correlation ratio, so correct this now
      if ( (numtoty<=1) || (var<=0.0) )
	return 0.0;   // the totally uncorrelated condition
      else
	return (1.0 - corr_ratio);

      // an alternative is to return 1.0/corr_ratio (=1/(1-correlation ratio))
      //  which may be better at rewarding gains near the best solution

      return 0;
    }

   // slow but general warp function... (doesn't need to be the same size as vref)
   /*
	xout(1) = x;  xout(2) = y;  xout(3) = z;
	xout = outvol.sampling_mat() * xout;
	//   apply inverse of postmat
	xout = postmat.i() * xout;
	//   convert xout to warpvol voxel coords
	xout = warpvol[0].sampling_mat().i() * xout;
	//   look up the warp dest coordinate (in mm)
	if (warpvol[0].in_bounds(MISCMATHS::round(xout(1)),
				 MISCMATHS::round(xout(2)),
				 MISCMATHS::round(xout(3)))) {
	  xin(1) = warpvol[0].interpolate(xout(1),xout(2),xout(3));
	  xin(2) = warpvol[1].interpolate(xout(1),xout(2),xout(3));
	  xin(3) = warpvol[2].interpolate(xout(1),xout(2),xout(3));
	  //   apply inverse of premat
	  xin = premat.i() * xin;
	  //   convert xin from mm to voxel coords
	  xin = invol.sampling_mat().i() * xin;
	  I_in = invol.interpolate(xin(1),xin(2),xin(3));
   */


  ///////////////////////////////////////////////////////////////////////

  float p_corr_ratio_gradient_fully_weighted(volume4D<float>& gradvec,
					     const volume<float>& vref, 
					     const volume<float>& vtest,
					     const volume<float>& refweight, 
					     const volume<float>& testweight,
					     int *bindex, 
					     const volume4D<float>& warpvol,
					     const int no_bins, 
					     const float smoothsize, 
					     bool nullbc)
    {

      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1
      // This assumes that warpvol is in the same space (and same sampling)
      //  as vref

      float  xb2 = ((float) vtest.xsize())-1.0001,
	yb2=((float) vtest.ysize())-1.0001, zb2=((float) vtest.zsize())-1.0001;

      float *sumy, *sumy2;
      sumy = new float[no_bins+1];
      sumy2 = new float[no_bins+1];
      float *numy;
      numy = new float[no_bins+1];
      int b=0;
 
      for (int i=0; i<=no_bins; i++) {
	numy[i]=0.0; sumy[i]=0.0;  sumy2[i]=0.0;
      }

      float wval=0,val=0,o1=0,o2=0,o3=0;

      float smoothx, smoothy, smoothz, weight;
      smoothx = smoothsize / vtest.xdim();
      smoothy = smoothsize / vtest.ydim();
      smoothz = smoothsize / vtest.zdim();

      int *bptr;

      for (int z=vref.minz(); z<=vref.maxz(); z++) {
	for (int y=vref.miny(); y<=vref.maxy(); y++) {
	  for (int x=vref.minx(); x<=vref.maxx(); x++) {
	    //   assume vref and warpvol have same voxel size
	    //   look up the warp dest coordinate (in mm)
	    if (warpvol[0].in_bounds(x,y,z)) {
	      o1 = warpvol[0](x,y,z) / vtest.xdim();
	      o2 = warpvol[1](x,y,z) / vtest.ydim();
	      o3 = warpvol[2](x,y,z) / vtest.zdim();

	      if (in_interp_bounds(vtest,o1,o2,o3)) {
		val = q_tri_interpolation(vtest,o1,o2,o3);
		wval = q_tri_interpolation(testweight,o1,o2,o3);
		
		// do the cost function record keeping...
		bptr = get_bindexptr(x,y,z,vref,bindex);
		b=*bptr;
		weight=wval*refweight(x,y,z);
		if (o1<smoothx)  weight*=o1/smoothx;
		else if ((xb2-o1)<smoothx) weight*=(xb2-o1)/smoothx;
		if (o2<smoothy)  weight*=o2/smoothy;
		else if ((yb2-o2)<smoothy) weight*=(yb2-o2)/smoothy;
		if (o3<smoothz)  weight*=o3/smoothz;
		else if ((zb2-o3)<smoothz) weight*=(zb2-o3)/smoothz;
		if (weight<0.0)  weight=0.0;
		numy[b]+=weight;
		sumy[b]+=weight*val;
		sumy2[b]+=weight*val*val;
	      }
	    }
	  }
	}
      }
      
      
      float corr_ratio=0.0, var=0.0, totsumy=0.0, totsumy2=0.0;
      float numtoty=0.0;

      // correct for occasion lapses into the last bin
      numy[no_bins-1] += numy[no_bins];
      sumy[no_bins-1] += sumy[no_bins];
      sumy2[no_bins-1] += sumy2[no_bins];
      numy[no_bins]=0.0;
      sumy[no_bins]=0.0;
      sumy2[no_bins]=0.0;

      // now calculate the individual variances for each iso-set
      //  weighting them by the number of pixels from Image x that contribute
      for (b=0; b<no_bins; b++) {
	if (numy[b]>2.0) {
	  numtoty += numy[b];
	  totsumy += sumy[b];
	  totsumy2 += sumy2[b];
	  // the following should be the variance of the bth iso-subset
	  var = (sumy2[b] - sumy[b]*sumy[b]/numy[b] ) / ( numy[b]-1);
	  corr_ratio += var * ((float) numy[b]);
	}
      }

      // normalise the weighting of numy[]
      if (numtoty>0)  corr_ratio/=((float) numtoty);
      // calculate the total variance of Image y and then normalise by this
      if (numtoty>1)
	var = ( totsumy2 - totsumy*totsumy/numtoty ) / (numtoty - 1);
      if (var>0.0)  corr_ratio/=var;
      // the above is actually 1 - correlation ratio, so correct this now
      if ( (numtoty<=1) || (var<=0.0) )
	corr_ratio = 0.0;   // the totally uncorrelated condition
      else
	corr_ratio = 1.0 - corr_ratio;
      
      ////////////////////////////////
      // now calculate the gradient //
      ////////////////////////////////
      gradvec = warpvol * 0.0f;   
      float delta_cr[2], delta_i, i0=0, i1=0;
      float do1=0, do2=0, do3=0;
      float nn, voxdim=1;
      bool prein=false, postin=false;

      for (int z=vref.minz(); z<=vref.maxz(); z++) {
	for (int y=vref.miny(); y<=vref.maxy(); y++) {
	  for (int x=vref.minx(); x<=vref.maxx(); x++) {
	    
	    bptr = get_bindexptr(x,y,z,vref,bindex);
	    b=*bptr;
	    nn = numy[b];
	    if (nn<1) nn=1;

	    if (warpvol[0].in_bounds(x,y,z)) {
	      o1 = warpvol[0](x,y,z) / vtest.xdim();
	      o2 = warpvol[1](x,y,z) / vtest.ydim();
	      o3 = warpvol[2](x,y,z) / vtest.zdim();
	    }
	    prein = in_interp_bounds(vtest,o1,o2,o3);
	    if (prein) { i0 = q_tri_interpolation(vtest,o1,o2,o3); }
	    
	    for (int dirn=0; dirn<=2; dirn++) {
	      // set the perturbed direction (in voxel units)
	      if (dirn==0) { do1=0.5; do2=0; do3=0; voxdim=vtest.xdim(); }
	      if (dirn==1) { do1=0; do2=0.5; do3=0; voxdim=vtest.ydim(); }
	      if (dirn==2) { do1=0; do2=0; do3=0.5; voxdim=vtest.zdim(); }
	      for (int diridx=0; diridx<=1; diridx++) {
		delta_cr[diridx] = 0;  // default
		// set the perturbed direction
		if (diridx==1) { do1=-do1; do2=-do2; do3=-do3; }
		postin = in_interp_bounds(vtest,o1+do1,o2+do2,o3+do3);
		if (postin) { 
		  i1 = q_tri_interpolation(vtest,o1+do1,o2+do2,o3+do3); 
		}
		if (prein && postin) {
		  // for both pre- and post-perturbed point inside the FOV
		  delta_i = i1 - i0;
		  delta_cr[diridx] = -delta_i * ( 2 * (i0 - sumy[b]/nn) + 
						  (1 - 1.0/nn)*delta_i ) / nn;
		}
		if (prein && !postin && !nullbc) {
		  // if the post-perturbed point is outside the FOV, but pre- is inside
		  delta_cr[diridx] = (sumy2[b] - sumy[b]*sumy[b]/nn) / (nn-1) 
		    - (sumy2[b] - i0*i0 - Sqr(sumy[b] - i0)/(nn-1)) / (nn-2);
		}
		if (postin && !prein && !nullbc) {
		  // if the pre-perturbed point is outside the FOV, but post- is inside
		  delta_cr[diridx] =  (sumy2[b] - sumy[b]*sumy[b]/nn) / (nn-1) 
		    - ( sumy2[b] + i1*i1 - Sqr(sumy[b] + i1) / (nn+1) ) / nn;
		}
	      }
	      
	      // check if the behaviour about this point is linear or quadratic
	      if ( (delta_cr[0]<0) && (delta_cr[1]<0) ) {
		// currently at a maximum, so pick the direction of biggest decrease
		if (delta_cr[0]<delta_cr[1]) {
		  gradvec[dirn](x,y,z) = delta_cr[0] / (0.5 * voxdim);
		} else {
		  // as diridx 1 is -ve mvmt
		  gradvec[dirn](x,y,z) = -delta_cr[1] / (0.5 * voxdim);  
		}
	      } else if ( (delta_cr[0]>0) && (delta_cr[1]>0) ) {
		// currently at a minimum, so zero gradient (don't want to move)
		gradvec[dirn](x,y,z) = 0;
	      } else {
		// average two values (but note that 0 = +ve mvmt, 1 = -ve mvmt)
		gradvec[dirn](x,y,z) = (delta_cr[0]-delta_cr[1]) / voxdim;
	      }
	    }
	  }
	}
      }
      
      return corr_ratio;
    }

  ///////////////////////////////////////////////////////////////////////

   float p_corr_ratio_smoothed(const volume<float>& vref, 
			     const volume<float>& vtest,
			     int *bindex, const Matrix& aff,
			     const int no_bins, const float smoothsize)
    {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1

      Matrix iaffbig = vtest.sampling_mat().i() * aff.i() *
	                     vref.sampling_mat();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.xsize()-1, yb1=vref.ysize()-1, zb1=vref.zsize()-1;
      float  xb2 = ((float) vtest.xsize())-1.0001,
	yb2=((float) vtest.ysize())-1.0001, zb2=((float) vtest.zsize())-1.0001;

      float *sumy, *sumy2;
      sumy = new float[no_bins+1];
      sumy2 = new float[no_bins+1];
      float *numy;
      numy = new float[no_bins+1];
      int b=0;

      for (int i=0; i<=no_bins; i++) {
	numy[i]=0.0; sumy[i]=0.0;  sumy2[i]=0.0;
      }

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4);
      float val=0,o1,o2,o3;

      float smoothx, smoothy, smoothz, weight;
      smoothx = smoothsize / vtest.xdim();
      smoothy = smoothsize / vtest.ydim();
      smoothz = smoothsize / vtest.zdim();

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      unsigned int xmin, xmax;
      int *bptr;

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  bptr = get_bindexptr(xmin,y,z,vref,bindex);

	  for (unsigned int x=xmin; x<=xmax; x++) {

	    if ( !((x==xmin) || (x==xmax)) 
		 || in_interp_bounds(vtest,o1,o2,o3) )
	      {
		val = q_tri_interpolation(vtest,o1,o2,o3);
	    
	        // do the cost function record keeping...
		b=*bptr;
		weight=1.0;
		if (o1<smoothx)  weight*=o1/smoothx;
		else if ((xb2-o1)<smoothx) weight*=(xb2-o1)/smoothx;
		if (o2<smoothy)  weight*=o2/smoothy;
		else if ((yb2-o2)<smoothy) weight*=(yb2-o2)/smoothy;
		if (o3<smoothz)  weight*=o3/smoothz;
		else if ((zb2-o3)<smoothz) weight*=(zb2-o3)/smoothz;
		if (weight<0.0)  weight=0.0;
		numy[b]+=weight;
		sumy[b]+=weight*val;
		sumy2[b]+=weight*val*val;
	      }

	    bptr++;
	    o1 += a11;
	    o2 += a21;
	    o3 += a31;

	  }
	}
      }

      float corr_ratio=0.0, var=0.0, totsumy=0.0, totsumy2=0.0;
      float numtoty=0.0;

      // correct for occasion lapses into the last bin
      numy[no_bins-1] += numy[no_bins];
      sumy[no_bins-1] += sumy[no_bins];
      sumy2[no_bins-1] += sumy2[no_bins];
      numy[no_bins]=0.0;
      sumy[no_bins]=0.0;
      sumy2[no_bins]=0.0;

      // now calculate the individual variances for each iso-set
      //  weighting them by the number of pixels from Image x that contribute
      for (b=0; b<no_bins; b++) {
	if (numy[b]>2.0) {
	  numtoty += numy[b];
	  totsumy += sumy[b];
	  totsumy2 += sumy2[b];
	  // the following should be the variance of the bth iso-subset
	  var = (sumy2[b] - sumy[b]*sumy[b]/numy[b] ) / ( numy[b]-1);
	  // cerr << "Set #" << b << " has " << numy[b] << " elements and " 
	  //   << var << " variance" << endl;
	  corr_ratio += var * ((float) numy[b]);
	}
      }
      delete [] numy; delete [] sumy; delete [] sumy2;

      // normalise the weighting of numy[]
      if (numtoty>0)  corr_ratio/=((float) numtoty);
      // calculate the total variance of Image y and then normalise by this
      if (numtoty>1)
	var = ( totsumy2 - totsumy*totsumy/numtoty ) / (numtoty - 1);
      //cerr << "TOTALS are:" << endl 
      //   << " numerator variance is : " << corr_ratio << endl
      //   << " and denominator variance is: " << var << " from " << numtoty 
      //   << " valid elements" << endl;
      if (var>0.0)  corr_ratio/=var;
      // the above is actually 1 - correlation ratio, so correct this now

      if ( (numtoty<=1) || (var<=0.0) )
	return 0.0;   // the totally uncorrelated condition
      else
	return (1.0 - corr_ratio);

      // an alternative is to return 1.0/corr_ratio (=1/(1-correlation ratio))
      //  which may be better at rewarding gains near the best solution

      return 0;

    }

  ///////////////////////////////////////////////////////////////////////

   float p_corr_ratio(const volume<float>& vref, const volume<float>& vtest,
		    int *bindex, const Matrix& aff,
		    const int no_bins)
    {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1

      Matrix iaffbig = vtest.sampling_mat().i() * aff.i() *
	                     vref.sampling_mat();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.xsize()-1, yb1=vref.ysize()-1, zb1=vref.zsize()-1;
      float  xb2 = ((float) vtest.xsize())-1.0001,
	yb2=((float) vtest.ysize())-1.0001, zb2=((float) vtest.zsize())-1.0001;

      float *sumy, *sumy2;
      sumy = new float[no_bins+1];
      sumy2 = new float[no_bins+1];
      int *numy;
      numy = new int[no_bins+1];
      int b=0;
 
      for (int i=0; i<=no_bins; i++) {
	numy[i]=0; sumy[i]=0.0;  sumy2[i]=0.0;
      }

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4);
      float val=0,o1,o2,o3;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      unsigned int xmin, xmax;
      int *bptr;

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  bptr = get_bindexptr(xmin,y,z,vref,bindex);

	  for (unsigned int x=xmin; x<=xmax; x++) {

	    if ( !((x==xmin) || (x==xmax)) 
		 || in_interp_bounds(vtest,o1,o2,o3) )
	      {
		val = q_tri_interpolation(vtest,o1,o2,o3);
		
		// do the cost function record keeping...
		b=*bptr;
		numy[b]++;
		sumy[b]+=val;
		sumy2[b]+=val*val;
	      }

	    bptr++;
	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	}
      }


      float corr_ratio=0.0, var=0.0, totsumy=0.0, totsumy2=0.0;
      int numtoty=0;

      // correct for occasion lapses into the last bin
      numy[no_bins-1] += numy[no_bins];
      sumy[no_bins-1] += sumy[no_bins];
      sumy2[no_bins-1] += sumy2[no_bins];
      numy[no_bins]=0;
      sumy[no_bins]=0.0;
      sumy2[no_bins]=0.0;

      // now calculate the individual variances for each iso-set
      //  weighting them by the number of pixels from Image x that contribute
      for (b=0; b<no_bins; b++) {
	if (numy[b]>2) {
	  numtoty += numy[b];
	  totsumy += sumy[b];
	  totsumy2 += sumy2[b];
	  // the following should be the variance of the bth iso-subset
	  var = (sumy2[b] - sumy[b]*sumy[b]/((float) numy[b]) ) /
	    ((float) (numy[b]-1));
	  // cerr << "Set #" << b << " has " << numy[b] << " elements and " 
	  //   << var << " variance" << endl;
	  corr_ratio += var * ((float) numy[b]);
	}
      }
      delete [] numy; delete [] sumy; delete [] sumy2;

      // normalise the weighting of numy[]
      if (numtoty>0)  corr_ratio/=((float) numtoty);
      // calculate the total variance of Image y and then normalise by this
      if (numtoty>1)
	var = ( totsumy2 - totsumy*totsumy/((float) numtoty) ) /
	  ((float) (numtoty - 1));
      //cerr << "TOTALS are:" << endl 
      //   << " numerator variance is : " << corr_ratio << endl
      //   << " and denominator variance is: " << var << " from " << numtoty 
      //   << " valid elements" << endl;
      if (var>0.0)  corr_ratio/=var;
      // the above is actually 1 - correlation ratio, so correct this now
      if ( (numtoty<=1) || (var<=0.0) )
	return 0.0;   // the totally uncorrelated condition
      else
	return (1.0 - corr_ratio);

      // an alternative is to return 1.0/corr_ratio (=1/(1-correlation ratio))
      //  which may be better at rewarding gains near the best solution

      return 0;

    }

  ///////////////////////////////////////////////////////////////////////


  float p_woods_fn(const volume<float>& vref, 
		 const volume<float>& vtest, int *bindex, 
		 const Matrix& aff, const int no_bins)
    {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1
      Matrix iaffbig = vtest.sampling_mat().i() * aff.i() *
	                     vref.sampling_mat();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.xsize()-1, yb1=vref.ysize()-1, zb1=vref.zsize()-1;
      float  xb2 = ((float) vtest.xsize())-1.0001,
	yb2=((float) vtest.ysize())-1.0001, zb2=((float) vtest.zsize())-1.0001;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      float *sum, *sum2;
      sum = new float[no_bins+1];
      sum2 = new float[no_bins+1];
      int *num;
      num = new int[no_bins+1];
      int b=0;

      for (int i=0; i<=no_bins; i++) {
	num[i]=0; sum[i]=0.0;  sum2[i]=0.0;
      }
  
      float val=0.0;
      unsigned int xmin, xmax;
      int *bptr;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  bptr = get_bindexptr(xmin,y,z,vref,bindex);

	  for (unsigned int x=xmin; x<=xmax; x++) {

	    if ( !((x==xmin) || (x==xmax)) 
		 || in_interp_bounds(vtest,o1,o2,o3) )
	      {
		val = q_tri_interpolation(vtest,o1,o2,o3);
		
		// do the cost function record keeping...
		b=*bptr;
		num[b]++;
		sum[b]+=val;
		sum2[b]+=val*val;
	      }
	    
	    bptr++;
	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	}
      }

      // now calculate  W = sum_j (n_j/N)*(sigma_j / mu_j)
      //  where n_j = num[j], N = sum_j n_j, mu_j = sum[j]/num[j]
      //        sigma_j^2 = sum_j 1/(n_j - 1) * (sum2[j] - mu_j^2 * n_j)
      float woods=0.0, stdev=0.0, var=0.0;
      int numtot=0;
      for (b=0; b<=no_bins; b++) {
	if (num[b]>2) {
	  numtot += num[b];
	  // the following should be the variance of the bth subset
	  var = (sum2[b] - sum[b]*sum[b]/((float) num[b]) ) /
	    ((float) (num[b]-1));
	  if (var>0.0)
	    stdev = sqrt(var);
	  else
	    stdev = 0.0;
	  if (sum[b]>0)
	    woods += Sqr((float) num[b])*stdev/sum[b];
	  else
	    woods += Sqr((float) num[b])*stdev;
	}
      }
      delete [] num; delete [] sum; delete [] sum2;
      if (numtot>0) {
	woods/=((float) numtot);
	return woods;
      } else {
	return 1e+10;
      }
    }


  ///////////////////////////////////////////////////////////////////////


  float p_woods_fn_smoothed(const volume<float>& vref, 
			  const volume<float>& vtest, int *bindex, 
			  const Matrix& aff, const int no_bins, 
			  const float smoothsize)
    {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1
      Matrix iaffbig = vtest.sampling_mat().i() * aff.i() *
	                     vref.sampling_mat();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.xsize()-1, yb1=vref.ysize()-1, zb1=vref.zsize()-1;
      float  xb2 = ((float) vtest.xsize())-1.0001,
	yb2=((float) vtest.ysize())-1.0001, zb2=((float) vtest.zsize())-1.0001;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      float *sum, *sum2;
      sum = new float[no_bins+1];
      sum2 = new float[no_bins+1];
      float *num;
      num = new float[no_bins+1];
      int b=0;

      for (int i=0; i<=no_bins; i++) {
	num[i]=0.0; sum[i]=0.0;  sum2[i]=0.0;
      }
  
      float smoothx, smoothy, smoothz, weight;
      smoothx = smoothsize / vtest.xdim();
      smoothy = smoothsize / vtest.ydim();
      smoothz = smoothsize / vtest.zdim();

      float val=0.0;
      unsigned int xmin, xmax;
      int *bptr;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  bptr = get_bindexptr(xmin,y,z,vref,bindex);

	  for (unsigned int x=xmin; x<=xmax; x++) {

	    if ( !((x==xmin) || (x==xmax)) 
		 || in_interp_bounds(vtest,o1,o2,o3) )
	      {
		val = q_tri_interpolation(vtest,o1,o2,o3);
		
		// do the cost function record keeping...
		b=*bptr;
		weight=1.0;
		if (o1<smoothx)  weight*=o1/smoothx;
		else if ((xb2-o1)<smoothx) weight*=(xb2-o1)/smoothx;
		if (o2<smoothy)  weight*=o2/smoothy;
		else if ((yb2-o2)<smoothy) weight*=(yb2-o2)/smoothy;
		if (o3<smoothz)  weight*=o3/smoothz;
		else if ((zb2-o3)<smoothz) weight*=(zb2-o3)/smoothz;
		if (weight<0.0)  weight=0.0;
		num[b]+=weight;
		sum[b]+=weight*val;
		sum2[b]+=weight*val*val;
	      }

	    bptr++;
	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	}
      }

      // now calculate  W = sum_j (n_j/N)*(sigma_j / mu_j)
      //  where n_j = num[j], N = sum_j n_j, mu_j = sum[j]/num[j]
      //        sigma_j^2 = sum_j 1/(n_j - 1) * (sum2[j] - mu_j^2 * n_j)
      float woods=0.0, stdev=0.0, var=0.0;
      float numtot=0.0;
      for (b=0; b<=no_bins; b++) {
	if (num[b]>2.0) {
	  numtot += num[b];
	  // the following should be the variance of the bth subset
	  var = (sum2[b] - sum[b]*sum[b]/(num[b])) / (num[b]-1.0);
	  if (var>0.0)
	    stdev = sqrt(var);
	  else
	    stdev = 0.0;
	  if (sum[b]>0)
	    woods += Sqr(num[b])*stdev/sum[b];
	  else
	    woods += Sqr(num[b])*stdev;
	}
      }
      delete [] num; delete [] sum; delete [] sum2;
      if (numtot>0) {
	woods/=numtot;
	return woods;
      } else {
	return 1e10;
      }
    }


  ///////////////////////////////////////////////////////////////////////


  float p_normcorr(const volume<float>& vref, const volume<float>& vtest,
		 const Matrix& aff)
    {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1
      Matrix iaffbig = vtest.sampling_mat().i() * aff.i() *
	                     vref.sampling_mat();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.xsize()-1, yb1=vref.ysize()-1, zb1=vref.zsize()-1;
      float  xb2 = ((float) vtest.xsize())-1.0001,
	yb2=((float) vtest.ysize())-1.0001, zb2=((float) vtest.zsize())-1.0001;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      float corr=0.0;
      float sumx=0.0, sumy=0.0, sumx2=0.0, sumy2=0.0, sumxy=0.0;
      float sumxA=0.0, sumyA=0.0, sumx2A=0.0, sumy2A=0.0, sumxyA=0.0;
      float sumxB=0.0, sumyB=0.0, sumx2B=0.0, sumy2B=0.0, sumxyB=0.0;
      float varx=0.0, vary=0.0, varxy=0.0;
      float valx=0.0, valy=0.0, val=0.0;
      long int num=0;

      unsigned int xmin, xmax;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  for (unsigned int x=xmin; x<=xmax; x++) {

	    if ( !((x==xmin) || (x==xmax)) 
		 || in_interp_bounds(vtest,o1,o2,o3) )
	      {
		val = q_tri_interpolation(vtest,o1,o2,o3);
		
		// do the cost function record keeping...
		num++;
		valx = vref(x,y,z);
		valy = val;
		sumx += valx;
		sumx2 += valx*valx;
		sumy += valy;
		sumy2 += valy*valy;
		sumxy += valx*valy;
	      }

	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	  sumxA+=sumx; sumyA+=sumy; 
	  sumx2A+=sumx2; sumy2A+=sumy2; sumxyA+=sumxy; 
	  sumx=0.0; sumy=0.0; sumx2=0.0; sumy2=0.0; sumxy=0.0;
	}
	sumxB+=sumxA; sumyB+=sumyA; 
	sumx2B+=sumx2A; sumy2B+=sumy2A; sumxyB+=sumxyA; 
	sumxA=0.0; sumyA=0.0; sumx2A=0.0; sumy2A=0.0; sumxyA=0.0;
      }
      assert(fabs(sumxA+sumx)<1e-9);
      sumx=sumxB; sumy=sumyB; sumx2=sumx2B; sumy2=sumy2B; sumxy=sumxyB;

      corr = 0.0;  // uncorrelated (worst) case
      if (num>2) {
	float numsq = ((float) num)*((float) num);
	varxy = sumxy/((float) num-1) - (sumx*sumy)/numsq;
	varx = sumx2/((float) num-1) - (sumx*sumx)/numsq;
	vary = sumy2/((float) num-1) - (sumy*sumy)/numsq;
	if ((varx>0.0) && (vary>0.0)) {
	  corr = varxy/sqrt(varx)/sqrt(vary);
	} 
      }
      return corr;
    }

  
  ///////////////////////////////////////////////////////////////////////


  float p_normcorr_smoothed(const volume<float>& vref, 
			    const volume<float>& vtest,
			    const Matrix& aff, const float smoothsize)
    {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1
      Matrix iaffbig = vtest.sampling_mat().i() * aff.i() *
	                     vref.sampling_mat();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.xsize()-1, yb1=vref.ysize()-1, zb1=vref.zsize()-1;
      float  xb2 = ((float) vtest.xsize())-1.0001,
	yb2=((float) vtest.ysize())-1.0001, zb2=((float) vtest.zsize())-1.0001;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      float corr=0.0;
      float sumx=0.0, sumy=0.0, sumx2=0.0, sumy2=0.0, sumxy=0.0;
      float sumxA=0.0, sumyA=0.0, sumx2A=0.0, sumy2A=0.0, sumxyA=0.0;
      float sumxB=0.0, sumyB=0.0, sumx2B=0.0, sumy2B=0.0, sumxyB=0.0;
      float varx=0.0, vary=0.0, varxy=0.0;
      float valx=0.0, valy=0.0, val=0.0;
      float num=0.0, numA=0.0, numB=0.0;

      float smoothx, smoothy, smoothz, weight;
      smoothx = smoothsize / vtest.xdim();
      smoothy = smoothsize / vtest.ydim();
      smoothz = smoothsize / vtest.zdim();

      unsigned int xmin, xmax;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      //cerr << "debug - normcorr_smoothed (1)" << endl;
    
      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 
	  
	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  //cerr << "xmin: " << xmin << endl;
	  //cerr << "xmax: " << xmax << endl;

	  for (unsigned int x=xmin; x<=xmax; x++) {

	    if ( !((x==xmin) || (x==xmax)) 
		 || in_interp_bounds(vtest,o1,o2,o3) )
	      {
		val = q_tri_interpolation(vtest,o1,o2,o3);
		
		// do the cost function record keeping...
		weight=1.0;
		if (o1<smoothx)  weight*=o1/smoothx;
		else if ((xb2-o1)<smoothx) weight*=(xb2-o1)/smoothx;
		if (o2<smoothy)  weight*=o2/smoothy;
		else if ((yb2-o2)<smoothy) weight*=(yb2-o2)/smoothy;
		if (o3<smoothz)  weight*=o3/smoothz;
		else if ((zb2-o3)<smoothz) weight*=(zb2-o3)/smoothz;
		if (weight<0.0)  weight=0.0;
		
		valx = vref(x,y,z);
		valy = val;
		num += weight;
		sumx += weight*valx;
		sumx2 += weight*valx*valx;
		sumy += weight*valy;
		sumy2 += weight*valy*valy;
		sumxy += weight*valx*valy;
	      }

	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }

	  numA+=num; sumxA+=sumx; sumyA+=sumy; 
	  sumx2A+=sumx2; sumy2A+=sumy2; sumxyA+=sumxy; 
	  sumx=0.0; sumy=0.0; sumx2=0.0; sumy2=0.0; sumxy=0.0;
	}

	numB+=numA; sumxB+=sumxA; sumyB+=sumyA; 
	sumx2B+=sumx2A; sumy2B+=sumy2A; sumxyB+=sumxyA; 
	sumxA=0.0; sumyA=0.0; sumx2A=0.0; sumy2A=0.0; sumxyA=0.0;
      }

      assert(fabs(sumxA+sumx)<1e-9);
      num = numB;
      sumx=sumxB; sumy=sumyB; sumx2=sumx2B; sumy2=sumy2B; sumxy=sumxyB;
  
      corr = 0.0;  // uncorrelated (worst) case
      if (num>2.0) {
	varxy = sumxy/(num-1.0) - (sumx*sumy)/(num*num);
	varx = sumx2/(num-1.0) - (sumx*sumx)/(num*num);
	vary = sumy2/(num-1.0) - (sumy*sumy)/(num*num);
	if ((varx>0.0) && (vary>0.0)) {
	  corr = varxy/sqrt(varx)/sqrt(vary);
	} 
      }
      return corr;
    }

  
  ///////////////////////////////////////////////////////////////////////

 
  float p_normcorr_smoothed_sinc(const volume<float>& vref, 
				 const volume<float>& vtest,
				 const Matrix& aff, const float smoothsize)
  {
    // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1

      ///vtest.setinterpolationmethod(sinc);

      Matrix iaffbig = vtest.sampling_mat().i() * aff.i() *
	                     vref.sampling_mat();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.xsize()-1, yb1=vref.ysize()-1, zb1=vref.zsize()-1;
      float  xb2 = ((float) vtest.xsize())-1.0001,
	yb2=((float) vtest.ysize())-1.0001, zb2=((float) vtest.zsize())-1.0001;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      float corr=0.0;
      float sumx=0.0, sumy=0.0, sumx2=0.0, sumy2=0.0, sumxy=0.0;
      float sumxA=0.0, sumyA=0.0, sumx2A=0.0, sumy2A=0.0, sumxyA=0.0;
      float sumxB=0.0, sumyB=0.0, sumx2B=0.0, sumy2B=0.0, sumxyB=0.0;
      float varx=0.0, vary=0.0, varxy=0.0;
      float valx=0.0, valy=0.0, val=0.0;
      float num=0.0, numA=0.0, numB=0.0;

      float smoothx, smoothy, smoothz, weight;
      smoothx = smoothsize / vtest.xdim();
      smoothy = smoothsize / vtest.ydim();
      smoothz = smoothsize / vtest.zdim();

      unsigned int xmin, xmax;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  for (unsigned int x=xmin; x<=xmax; x++) {

	    if ( !((x==xmin) || (x==xmax)) 
		 || in_interp_bounds(vtest,o1,o2,o3) )
	      {
		val = q_sinc_interpolation(vtest,o1,o2,o3);
		///val = vtest.interpolate(o1,o2,o3);
		
		// do the cost function record keeping...
		weight=1.0;
		if (o1<smoothx)  weight*=o1/smoothx;
		else if ((xb2-o1)<smoothx) weight*=(xb2-o1)/smoothx;
		if (o2<smoothy)  weight*=o2/smoothy;
		else if ((yb2-o2)<smoothy) weight*=(yb2-o2)/smoothy;
		if (o3<smoothz)  weight*=o3/smoothz;
		else if ((zb2-o3)<smoothz) weight*=(zb2-o3)/smoothz;
		if (weight<0.0)  weight=0.0;
		
		valx = vref(x,y,z);
		valy = val;
		num += weight;
		sumx += weight*valx;
		sumx2 += weight*valx*valx;
		sumy += weight*valy;
		sumy2 += weight*valy*valy;
		sumxy += weight*valx*valy;
	      }

	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	  numA+=num; sumxA+=sumx; sumyA+=sumy; 
	  sumx2A+=sumx2; sumy2A+=sumy2; sumxyA+=sumxy; 
	  sumx=0.0; sumy=0.0; sumx2=0.0; sumy2=0.0; sumxy=0.0;
	}
	numB+=numA; sumxB+=sumxA; sumyB+=sumyA; 
	sumx2B+=sumx2A; sumy2B+=sumy2A; sumxyB+=sumxyA; 
	sumxA=0.0; sumyA=0.0; sumx2A=0.0; sumy2A=0.0; sumxyA=0.0;
      }
      assert(fabs(sumxA+sumx)<1e-9);
      num = numB;
      sumx=sumxB; sumy=sumyB; sumx2=sumx2B; sumy2=sumy2B; sumxy=sumxyB;
  
      corr = 0.0;  // uncorrelated (worst) case
      if (num>2.0) {
	varxy = sumxy/(num-1.0) - (sumx*sumy)/(num*num);
	varx = sumx2/(num-1.0) - (sumx*sumx)/(num*num);
	vary = sumy2/(num-1.0) - (sumy*sumy)/(num*num);
	if ((varx>0.0) && (vary>0.0)) {
	  corr = varxy/sqrt(varx)/sqrt(vary);
	} 
      }
      return corr;
    }

  


  ///////////////////////////////////////////////////////////////////////


  float p_normcorr_fully_weighted(const volume<float>& vref, 
				const volume<float>& vtest,
				const volume<float>& refweight, 
				const volume<float>& testweight,
				const Matrix& aff, const float smoothsize)
    {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1
      Matrix iaffbig = vtest.sampling_mat().i() * aff.i() *
	                     vref.sampling_mat();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.xsize()-1, yb1=vref.ysize()-1, zb1=vref.zsize()-1;
      float  xb2 = ((float) vtest.xsize())-1.0001,
	yb2=((float) vtest.ysize())-1.0001, zb2=((float) vtest.zsize())-1.0001;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      float corr=0.0;
      float sumx=0.0, sumy=0.0, sumx2=0.0, sumy2=0.0, sumxy=0.0;
      float sumxA=0.0, sumyA=0.0, sumx2A=0.0, sumy2A=0.0, sumxyA=0.0;
      float sumxB=0.0, sumyB=0.0, sumx2B=0.0, sumy2B=0.0, sumxyB=0.0;
      float varx=0.0, vary=0.0, varxy=0.0;
      float valx=0.0, valy=0.0, val=0.0;
      float num=0.0, numA=0.0, numB=0.0;

      float wval=0, smoothx, smoothy, smoothz, weight;
      smoothx = smoothsize / vtest.xdim();
      smoothy = smoothsize / vtest.ydim();
      smoothz = smoothsize / vtest.zdim();

      unsigned int xmin, xmax;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  for (unsigned int x=xmin; x<=xmax; x++) {

	    if ( !((x==xmin) || (x==xmax)) 
		 || in_interp_bounds(vtest,o1,o2,o3) )
	      {
		val = q_tri_interpolation(vtest,o1,o2,o3);
		wval = q_tri_interpolation(testweight,o1,o2,o3);
		
		// do the cost function record keeping...
		weight=wval*refweight(x,y,z);
		if (o1<smoothx)  weight*=o1/smoothx;
		else if ((xb2-o1)<smoothx) weight*=(xb2-o1)/smoothx;
		if (o2<smoothy)  weight*=o2/smoothy;
		else if ((yb2-o2)<smoothy) weight*=(yb2-o2)/smoothy;
		if (o3<smoothz)  weight*=o3/smoothz;
		else if ((zb2-o3)<smoothz) weight*=(zb2-o3)/smoothz;
		if (weight<0.0)  weight=0.0;
		
		valx = vref(x,y,z);
		valy = val;
		num += weight;
		sumx += weight*valx;
		sumx2 += weight*valx*valx;
		sumy += weight*valy;
		sumy2 += weight*valy*valy;
		sumxy += weight*valx*valy;
	      }
	    
	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	  numA+=num; sumxA+=sumx; sumyA+=sumy; 
	  sumx2A+=sumx2; sumy2A+=sumy2; sumxyA+=sumxy; 
	  sumx=0.0; sumy=0.0; sumx2=0.0; sumy2=0.0; sumxy=0.0;
	}
	numB+=numA; sumxB+=sumxA; sumyB+=sumyA; 
	sumx2B+=sumx2A; sumy2B+=sumy2A; sumxyB+=sumxyA; 
	sumxA=0.0; sumyA=0.0; sumx2A=0.0; sumy2A=0.0; sumxyA=0.0;
      }
      assert(fabs(sumxA+sumx)<1e-9);
      num = numB;
      sumx=sumxB; sumy=sumyB; sumx2=sumx2B; sumy2=sumy2B; sumxy=sumxyB;
  
      corr = 0.0;  // uncorrelated (worst) case
      if (num>2.0) {
	varxy = sumxy/(num-1.0) - (sumx*sumy)/(num*num);
	varx = sumx2/(num-1.0) - (sumx*sumx)/(num*num);
	vary = sumy2/(num-1.0) - (sumy*sumy)/(num*num);
	if ((varx>0.0) && (vary>0.0)) {
	  corr = varxy/sqrt(varx)/sqrt(vary);
	} 
      }
      return corr;
    }

  
  ///////////////////////////////////////////////////////////////////////


  float p_leastsquares(const volume<float>& vref, const volume<float>& vtest,
		     const Matrix& aff)
    {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1
      Matrix iaffbig = vtest.sampling_mat().i() * aff.i() *
	                     vref.sampling_mat();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.xsize()-1, yb1=vref.ysize()-1, zb1=vref.zsize()-1;
      float  xb2 = ((float) vtest.xsize())-1.0001,
	yb2=((float) vtest.ysize())-1.0001, zb2=((float) vtest.zsize())-1.0001;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      float lsq=0.0;
      float sum=0.0, sumA=0.0, sumB=0.0;
      float valx=0.0, valy=0.0, val=0.0;
      long int num=0;

      unsigned int xmin, xmax;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  for (unsigned int x=xmin; x<=xmax; x++) {

	    if ( !((x==xmin) || (x==xmax)) 
		 || in_interp_bounds(vtest,o1,o2,o3) )
	      {
		val = q_tri_interpolation(vtest,o1,o2,o3);
		
		// do the cost function record keeping...
		num++;
		valx = vref(x,y,z);
		valy = val;
		sum += (valx-valy)*(valx-valy);
	      }

	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	  sumA+=sum; sum=0.0;
	}
	sumB+=sumA; sumA=0.0;
      }
      assert(fabs(sumA+sum)<1e-9);
      sum = sumB;

      if (num>1) {
	lsq = sum/((float) num);
      } else {
	  // return the worst cost = (max-min)^2
	lsq = (Max(vref.max(),vtest.max())-Min(vref.min(),vtest.min()));
	lsq = lsq*lsq;
      }
      
      return lsq;
    }


  ///////////////////////////////////////////////////////////////////////


  float p_leastsquares_smoothed(const volume<float>& vref, 
			      const volume<float>& vtest,
			      const Matrix& aff, const float smoothsize)
    {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1
      Matrix iaffbig = vtest.sampling_mat().i() * aff.i() *
	                     vref.sampling_mat();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.xsize()-1, yb1=vref.ysize()-1, zb1=vref.zsize()-1;
      float  xb2 = ((float) vtest.xsize())-1.0001,
	yb2=((float) vtest.ysize())-1.0001, zb2=((float) vtest.zsize())-1.0001;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      float smoothx, smoothy, smoothz, weight;
      smoothx = smoothsize / vtest.xdim();
      smoothy = smoothsize / vtest.ydim();
      smoothz = smoothsize / vtest.zdim();

      float lsq=0.0;
      float sum=0.0, sumA=0.0, sumB=0.0;
      float valx=0.0, valy=0.0, val=0.0;
      float num=0.0, numA=0.0, numB=0.0;

      unsigned int xmin, xmax;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  for (unsigned int x=xmin; x<=xmax; x++) {

	    if ( !((x==xmin) || (x==xmax)) 
		 || in_interp_bounds(vtest,o1,o2,o3) )
	      {
		val = q_tri_interpolation(vtest,o1,o2,o3);
		
		// do the cost function record keeping...
		weight=1.0;
		if (o1<smoothx)  weight*=o1/smoothx;
		else if ((xb2-o1)<smoothx) weight*=(xb2-o1)/smoothx;
		if (o2<smoothy)  weight*=o2/smoothy;
		else if ((yb2-o2)<smoothy) weight*=(yb2-o2)/smoothy;
		if (o3<smoothz)  weight*=o3/smoothz;
		else if ((zb2-o3)<smoothz) weight*=(zb2-o3)/smoothz;
		if (weight<0.0)  weight=0.0;
		
		valx = vref(x,y,z);
		valy = val;
		num+=weight;
		sum += weight*(valx-valy)*(valx-valy);
	      }

	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	  sumA+=sum; sum=0.0;
	  numA+=num; num=0.0;
	}
	sumB+=sumA; sumA=0.0;
	numB+=numA; numA=0.0;
      }
      assert(fabs(sumA+sum)<1e-9);
      sum = sumB;  num = numB;
  
      
      if (num>1.0) {
	lsq = sum/num;
      } else {
	  // return the worst cost = (max-min)^2
	lsq = (Max(vref.max(),vtest.max())-Min(vref.min(),vtest.min()));
	lsq = lsq*lsq;
      }
      
      return lsq;
    }


  ///////////////////////////////////////////////////////////////////////


  float p_leastsquares_fully_weighted(const volume<float>& vref, 
				    const volume<float>& vtest,
				    const volume<float>& refweight, 
				    const volume<float>& testweight,
				    const Matrix& aff, const float smoothsize)
  {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1
      Matrix iaffbig = vtest.sampling_mat().i() * aff.i() *
	                     vref.sampling_mat();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.xsize()-1, yb1=vref.ysize()-1, zb1=vref.zsize()-1;
      float  xb2 = ((float) vtest.xsize())-1.0001,
	yb2=((float) vtest.ysize())-1.0001, zb2=((float) vtest.zsize())-1.0001;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      float wval=0, smoothx, smoothy, smoothz, weight;
      smoothx = smoothsize / vtest.xdim();
      smoothy = smoothsize / vtest.ydim();
      smoothz = smoothsize / vtest.zdim();

      float lsq=0.0;
      float sum=0.0, sumA=0.0, sumB=0.0;
      float valx=0.0, valy=0.0, val=0.0;
      float num=0.0, numA=0.0, numB=0.0;

      unsigned int xmin, xmax;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  for (unsigned int x=xmin; x<=xmax; x++) {

	    if ( !((x==xmin) || (x==xmax)) 
		 || in_interp_bounds(vtest,o1,o2,o3) )
	      {
		val = q_tri_interpolation(vtest,o1,o2,o3);
		wval = q_tri_interpolation(testweight,o1,o2,o3);
		
		// do the cost function record keeping...
		weight=wval*refweight(x,y,z);
		if (o1<smoothx)  weight*=o1/smoothx;
		else if ((xb2-o1)<smoothx) weight*=(xb2-o1)/smoothx;
		if (o2<smoothy)  weight*=o2/smoothy;
		else if ((yb2-o2)<smoothy) weight*=(yb2-o2)/smoothy;
		if (o3<smoothz)  weight*=o3/smoothz;
		else if ((zb2-o3)<smoothz) weight*=(zb2-o3)/smoothz;
		if (weight<0.0)  weight=0.0;
		
		valx = vref(x,y,z);
		valy = val;
		num+=weight;
		sum += weight*(valx-valy)*(valx-valy);
	      }

	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	  sumA+=sum; sum=0.0;
	  numA+=num; num=0.0;
	}
	sumB+=sumA; sumA=0.0;
	numB+=numA; numA=0.0;
      }
      assert(fabs(sumA+sum)<1e-9);
      sum = sumB;  num = numB;
  
      
      if (num>1.0) {
	lsq = sum/num;
      } else {
	  // return the worst cost = (max-min)^2
	lsq = (Max(vref.max(),vtest.max())-Min(vref.min(),vtest.min()));
	lsq = lsq*lsq;
      }
      
      return lsq;
    }


  ///////////////////////////////////////////////////////////////////////

  float p_labeldiff(const volume<float>& vref, const volume<float>& vtest,
		    const Matrix& aff)
    {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1
      Matrix iaffbig = vtest.sampling_mat().i() * aff.i() *
	                     vref.sampling_mat();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.xsize()-1, yb1=vref.ysize()-1, zb1=vref.zsize()-1;
      float  xb2 = ((float) vtest.xsize())-1.0001,
	yb2=((float) vtest.ysize())-1.0001, zb2=((float) vtest.zsize())-1.0001;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      float v000, v001, v010, v011, v100, v101, v110, v111;

      float lsq=0.0;
      float sum=0.0, sumA=0.0, sumB=0.0, lsum=0.0;
      float valx=0.0, dx, dy, dz;
      long int num=0;

      unsigned int xmin, xmax;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  for (unsigned int x=xmin; x<=xmax; x++) {

	    if ( !((x==xmin) || (x==xmax)) 
		 || in_interp_bounds(vtest,o1,o2,o3) )
	      {
		// do the cost function record keeping...
		num++;
		valx = vref(x,y,z);
		q_get_neighbours(vtest,o1,o2,o3,v000,v001,v010,v011,v100,v101,v110,v111,
				 dx, dy, dz);
		lsum=0.0;
		if (fabs(v000-valx)>0.5) lsum += (1.0-dx)*(1.0-dy)*(1.0-dz);
		if (fabs(v001-valx)>0.5) lsum += (1.0-dx)*(1.0-dy)*dz;
		if (fabs(v011-valx)>0.5) lsum += (1.0-dx)*dy*dz;
		if (fabs(v010-valx)>0.5) lsum += (1.0-dx)*dy*(1.0-dz);
		if (fabs(v110-valx)>0.5) lsum += dx*dy*(1.0-dz);
		if (fabs(v100-valx)>0.5) lsum += dx*(1.0-dy)*(1.0-dz);
		if (fabs(v101-valx)>0.5) lsum += dx*(1.0-dy)*dz;
		if (fabs(v111-valx)>0.5) lsum += dx*dy*dz;
		sum += lsum;
	      }

	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	  sumA+=sum; sum=0.0;
	}
	sumB+=sumA; sumA=0.0;
      }
      assert(fabs(sumA+sum)<1e-9);
      sum = sumB;

      if (num>1) {
	lsq = sum/((float) num);
      } else {
	  // return the worst cost = (max-min)^2
	lsq = (Max(vref.max(),vtest.max())-Min(vref.min(),vtest.min()));
	lsq = lsq*lsq;
      }
      
      return lsq;
    }


  ///////////////////////////////////////////////////////////////////////


  float p_labeldiff_smoothed(const volume<float>& vref, 
			     const volume<float>& vtest,
			     const Matrix& aff, const float smoothsize)
    {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1
      Matrix iaffbig = vtest.sampling_mat().i() * aff.i() *
	                     vref.sampling_mat();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.xsize()-1, yb1=vref.ysize()-1, zb1=vref.zsize()-1;
      float  xb2 = ((float) vtest.xsize())-1.0001,
	yb2=((float) vtest.ysize())-1.0001, zb2=((float) vtest.zsize())-1.0001;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      float smoothx, smoothy, smoothz, weight;
      smoothx = smoothsize / vtest.xdim();
      smoothy = smoothsize / vtest.ydim();
      smoothz = smoothsize / vtest.zdim();

      float v000, v001, v010, v011, v100, v101, v110, v111;
      float lsq=0.0;
      float sum=0.0, sumA=0.0, sumB=0.0, lsum=0.0;
      float valx=0.0, dx, dy, dz;
      float num=0.0, numA=0.0, numB=0.0;

      unsigned int xmin, xmax;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  for (unsigned int x=xmin; x<=xmax; x++) {

	    if ( !((x==xmin) || (x==xmax)) 
		 || in_interp_bounds(vtest,o1,o2,o3) )
	      {

		valx = vref(x,y,z);
	
		// do the cost function record keeping...
		weight=1.0;
		if (o1<smoothx)  weight*=o1/smoothx;
		else if ((xb2-o1)<smoothx) weight*=(xb2-o1)/smoothx;
		if (o2<smoothy)  weight*=o2/smoothy;
		else if ((yb2-o2)<smoothy) weight*=(yb2-o2)/smoothy;
		if (o3<smoothz)  weight*=o3/smoothz;
		else if ((zb2-o3)<smoothz) weight*=(zb2-o3)/smoothz;
		if (weight<0.0)  weight=0.0;
		
		num+=weight;
		q_get_neighbours(vtest,o1,o2,o3,v000,v001,v010,v011,v100,v101,v110,v111,
				 dx, dy,dz);
		lsum=0.0;
		if (fabs(v000-valx)>0.5) lsum += (1.0-dx)*(1.0-dy)*(1.0-dz);
		if (fabs(v001-valx)>0.5) lsum += (1.0-dx)*(1.0-dy)*dz;
		if (fabs(v011-valx)>0.5) lsum += (1.0-dx)*dy*dz;
		if (fabs(v010-valx)>0.5) lsum += (1.0-dx)*dy*(1.0-dz);
		if (fabs(v110-valx)>0.5) lsum += dx*dy*(1.0-dz);
		if (fabs(v100-valx)>0.5) lsum += dx*(1.0-dy)*(1.0-dz);
		if (fabs(v101-valx)>0.5) lsum += dx*(1.0-dy)*dz;
		if (fabs(v111-valx)>0.5) lsum += dx*dy*dz;

		sum += weight*lsum;
	      }

	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	  sumA+=sum; sum=0.0;
	  numA+=num; num=0.0;
	}
	sumB+=sumA; sumA=0.0;
	numB+=numA; numA=0.0;
      }
      assert(fabs(sumA+sum)<1e-9);
      sum = sumB;  num = numB;
  
      
      if (num>1.0) {
	lsq = sum/num;
      } else {
	  // return the worst cost = (max-min)^2
	lsq = (Max(vref.max(),vtest.max())-Min(vref.min(),vtest.min()));
	lsq = lsq*lsq;
      }
      
      return lsq;
    }


  ///////////////////////////////////////////////////////////////////////


  float p_labeldiff_fully_weighted(const volume<float>& vref, 
				   const volume<float>& vtest,
				   const volume<float>& refweight, 
				   const volume<float>& testweight,
				   const Matrix& aff, const float smoothsize)
  {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1
      Matrix iaffbig = vtest.sampling_mat().i() * aff.i() *
	                     vref.sampling_mat();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.xsize()-1, yb1=vref.ysize()-1, zb1=vref.zsize()-1;
      float  xb2 = ((float) vtest.xsize())-1.0001,
	yb2=((float) vtest.ysize())-1.0001, zb2=((float) vtest.zsize())-1.0001;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      float wval=0, smoothx, smoothy, smoothz, weight;
      smoothx = smoothsize / vtest.xdim();
      smoothy = smoothsize / vtest.ydim();
      smoothz = smoothsize / vtest.zdim();

      float v000, v001, v010, v011, v100, v101, v110, v111;
      float lsq=0.0;
      float sum=0.0, sumA=0.0, sumB=0.0, lsum=0.0;
      float valx=0.0, dx, dy, dz;
      float num=0.0, numA=0.0, numB=0.0;

      unsigned int xmin, xmax;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  for (unsigned int x=xmin; x<=xmax; x++) {

	    if ( !((x==xmin) || (x==xmax)) 
		 || in_interp_bounds(vtest,o1,o2,o3) )
	      {
		wval = q_tri_interpolation(testweight,o1,o2,o3);
		
		// do the cost function record keeping...
		weight=wval*refweight(x,y,z);
		if (o1<smoothx)  weight*=o1/smoothx;
		else if ((xb2-o1)<smoothx) weight*=(xb2-o1)/smoothx;
		if (o2<smoothy)  weight*=o2/smoothy;
		else if ((yb2-o2)<smoothy) weight*=(yb2-o2)/smoothy;
		if (o3<smoothz)  weight*=o3/smoothz;
		else if ((zb2-o3)<smoothz) weight*=(zb2-o3)/smoothz;
		if (weight<0.0)  weight=0.0;
		
		valx = vref(x,y,z);

		num+=weight;
		q_get_neighbours(vtest,o1,o2,o3,v000,v001,v010,v011,v100,v101,v110,v111,
				 dx,dy,dz);
		lsum=0.0;
		if (fabs(v000-valx)>0.5) lsum += (1.0-dx)*(1.0-dy)*(1.0-dz);
		if (fabs(v001-valx)>0.5) lsum += (1.0-dx)*(1.0-dy)*dz;
		if (fabs(v011-valx)>0.5) lsum += (1.0-dx)*dy*dz;
		if (fabs(v010-valx)>0.5) lsum += (1.0-dx)*dy*(1.0-dz);
		if (fabs(v110-valx)>0.5) lsum += dx*dy*(1.0-dz);
		if (fabs(v100-valx)>0.5) lsum += dx*(1.0-dy)*(1.0-dz);
		if (fabs(v101-valx)>0.5) lsum += dx*(1.0-dy)*dz;
		if (fabs(v111-valx)>0.5) lsum += dx*dy*dz;

		sum += weight*lsum;
	      }

	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	  sumA+=sum; sum=0.0;
	  numA+=num; num=0.0;
	}
	sumB+=sumA; sumA=0.0;
	numB+=numA; numA=0.0;
      }
      assert(fabs(sumA+sum)<1e-9);
      sum = sumB;  num = numB;
  
      
      if (num>1.0) {
	lsq = sum/num;
      } else {
	  // return the worst cost = (max-min)^2
	lsq = (Max(vref.max(),vtest.max())-Min(vref.min(),vtest.min()));
	lsq = lsq*lsq;
      }
      
      return lsq;
    }

  ///////////////////////////////////////////////////////////////////////

  float approx1tanh(const float x)
  {
    float val=0.0;
    if (x<-4.0) { val=0.0; } 
    else if (x<-2.0) { val=0.4 + 0.1*x; }
    else if (x<2.0) { val=1.0 + 0.4*x; }
    else if (x<4.0) { val=1.6 + 0.1*x; }
    else { val=2.0; }
    return val;
  }


  float Costfn::fmap_extrap(const double& x_vox, const double& y_vox, const double& z_vox, const ColumnVector& v_pe) const
  {
    double fmap_max_size=Max(Max(fmap.xsize()*fmap.xdim(),fmap.ysize()*fmap.ydim()),fmap.zsize()*fmap.zdim());
    for (double dist=0.0; dist<=fmap_max_size; dist+=1.0) {
      for (int dir=-1; dir<=1; dir+=2) {
	float x_ref=x_vox+dir*dist*v_pe(1);
	float y_ref=y_vox+dir*dist*v_pe(2);
	float z_ref=z_vox+dir*dist*v_pe(3);
	if (fmap_mask.in_bounds(x_ref,y_ref,z_ref) && fmap_mask.interpolate(x_ref,y_ref,z_ref)>0.95) {
	  return fmap.interpolate(x_ref,y_ref,z_ref);
	}
      }
    }
    return 0.0;
  }



  int Costfn::vox_coord_calc(ColumnVector& tvc, ColumnVector& rvc, const Matrix& aff, const ColumnVector& nonlin_params, 
			     const Matrix& iaffbig, const Matrix& mm2vox, const ColumnVector& pe_dir_vec) const
  {
    // NOTE! input coordinates rvc are in mm **but get converted to vox** during this call  (tvc output as vox)
    tvc = iaffbig * rvc;  // vox
    rvc = mm2vox * rvc;  // vox
    // tvc is currently the undistorted coords in the EPI (test image)
    if (pe_dir!=0) {
      // fmap gives the distance required to shift from undistorted to distorted coords
      //   and fmap will be in the reference frame
      // add in the fmap transformations (locked along the EPI PE direction)
      // multiply fmap by a free parameter - nonlin_params(1)
      // convert rvc to voxel coords
      if (fmap_mask.interpolate(rvc(1),rvc(2),rvc(3)) < 0.95) {
	tvc(abs(pe_dir)) += nonlin_params(1)*fmap_extrap(rvc(1),rvc(2),rvc(3),pe_dir_vec);
      } else {
	tvc(abs(pe_dir)) += nonlin_params(1)*fmap.interpolate(rvc(1),rvc(2),rvc(3));
      }
    }
    return 0;
  }


  float Costfn::bbr(const Matrix& aff, const ColumnVector& nonlin_params) const
  {
    // use the wmseg image as a dummy image (it will not be changed, as resample_required is set to false)
    volume<float> dummy;
    return bbr(aff,nonlin_params,dummy,false);
  }

  float Costfn::bbr_resamp(const Matrix& aff, const ColumnVector& nonlin_params, volume<float>& resampvol) const
  {
    // use the wmseg image as a dummy image (it will not be changed, as resample_required is set to false)
    return bbr(aff,nonlin_params,resampvol,true);
  }

  float Costfn::bbr(const Matrix& aff, const ColumnVector& nonlin_params, 
		    volume<float>& resampvol, bool resampling_required) const
    {
      p_count++;
      // Input GM and WM coordinates are in (flirt) mm coords in the ref space
      Matrix iaffbig = testvol.sampling_mat().i() * aff.i();  // ref mm to test vox
      Matrix mm2vox = refvol.sampling_mat().i();  // ref mm to ref vox
      ColumnVector pe_dir_vec(4);
      pe_dir_vec(1)=0.0;  pe_dir_vec(2)=1.0;  pe_dir_vec(3)=0.0;  pe_dir_vec(4)=1.0;
      pe_dir_vec = refvol.sampling_mat().i() * aff * pe_dir_vec;
      pe_dir_vec /= norm2(pe_dir_vec);   // "unit" vector, but in voxel coords

      double bbrval=0.0, gv=0.0, wv=0.0, qv=0.0, bbrupdate=0.0, bbrwsum=0.0, weight1=1.0, weight2=1.0, w12=1.0;
      //float echo_spacing = 5e-4;   // in seconds - a reasonable guess for modern scanners
      //fmapscale = echo_spacing * tsize(pe_dir) / (2.0*M_PI);  // converts rad/s to shift in voxels
      
      // rvc = ref voxel coord, tvc = test voxel coord
      ColumnVector rvc1(4), rvc2(4), tvc1(4), tvc2(4);
      
      if (debug_mode) { debugvol = refvol; debugvol.addvolume(refvol); debugvol=0.0f; }
	  
      if (!resampling_required) {
	for (int row=0; row<no_coords; row+=vertex_step) {
	  // convert to voxel coords (rvc=ref voxel coord, tvc=test vox coord)
	  rvc1(1)=gm_coord_x[row]; rvc1(2)=gm_coord_y[row]; rvc1(3)=gm_coord_z[row]; rvc1(4)=1;  // mm
	  rvc2(1)=wm_coord_x[row]; rvc2(2)=wm_coord_y[row]; rvc2(3)=wm_coord_z[row]; rvc2(4)=1;
	  
	  vox_coord_calc(tvc1, rvc1, aff, nonlin_params, iaffbig, mm2vox, pe_dir_vec);
	  vox_coord_calc(tvc2, rvc2, aff, nonlin_params, iaffbig, mm2vox, pe_dir_vec);
	  
	  gv = testvol.interpolate(tvc1(1),tvc1(2),tvc1(3));
	  wv = testvol.interpolate(tvc2(1),tvc2(2),tvc2(3));
	  
	  if (validweights) {
	    weight1 = tweight.interpolate(tvc1(1),tvc1(2),tvc1(3)) * rweight.interpolate(rvc1(1),rvc1(2),rvc1(3));
	    weight2 = tweight.interpolate(tvc2(1),tvc2(2),tvc2(3)) * rweight.interpolate(rvc2(1),rvc2(2),rvc2(3));
	    w12 = weight1 * weight2;
	  } else {
	    w12=1.0;  // need this if w12 can be reset after this (but not when validweights is on)
	  }
	  
	  
	  if (fabs(gv+wv)>1e-6) {
	    qv = 200.0 * (gv-wv) / (gv+wv);
	  } else {
	    qv=0.0;
	  }
	  
	  bbrupdate = 1 + tanh(bbr_slope*(qv-bbr_offset));
	  if (bbr_type=="local_abs") {
	    bbrval += w12*Min(bbrupdate,2.0-bbrupdate);
	  } else {
	    bbrval += w12*bbrupdate;
	  }
	  bbrwsum += w12;
	  if (debug_mode) { 
	    debugvol(MISCMATHS::round((rvc1(1)+rvc2(1))/2),MISCMATHS::round((rvc1(2)+rvc2(2))/2),MISCMATHS::round((rvc1(3)+rvc2(3))/2),0)=w12*bbrupdate+(1.0-w12); 
	    debugvol(MISCMATHS::round((rvc1(1)+rvc2(1))/2),MISCMATHS::round((rvc1(2)+rvc2(2))/2),MISCMATHS::round((rvc1(3)+rvc2(3))/2),1)=w12; 
	  }
	  
	  if (std::isnan(bbrval)) { cerr << "WARNING:: Found NaN in BBR cost" << endl; }
	}
	
	bbrval /= bbrwsum;
	if (bbr_type=="global_abs") {
	  bbrval = Min(bbrval,2.0-bbrval);
	}
	
	if (std::isnan(bbrval)) { cerr << "WARNING:: Found NaN in BBR cost (2)" << endl; }
	
	if (debug_mode) { save_volume4D(debugvol,"debug"); }
      }

      // optional resampling
      if (resampling_required) {
	resampvol=refvol;
	ColumnVector rvc(4), tvc(4);
	rvc(4)=1; tvc(4)=1;
	for (int refx=refvol.minx(); refx<=refvol.maxx(); refx++) {
	  for (int refy=refvol.miny(); refy<=refvol.maxy(); refy++) {
	    for (int refz=refvol.minz(); refz<=refvol.maxz(); refz++) {
	      // convert coords to mm (for the vox_coord_calc function)
	      rvc(1)=refx*refvol.xdim();
	      rvc(2)=refy*refvol.ydim();
	      rvc(3)=refz*refvol.zdim();
	      vox_coord_calc(tvc, rvc, aff, nonlin_params, iaffbig, mm2vox, pe_dir_vec);
	      resampvol(refx,refy,refz)=testvol.interpolate(tvc(1),tvc(2),tvc(3));
	    }
	  }
	}
      }

      return bbrval;
    }

  ///////////////////////////////////////////////////////////////////////


  void calc_entropy(const volume<float>& vref, const volume<float>& vtest,
		    int *bindex,  const Matrix& aff,
		    const float mintest, const float maxtest,
		    const int no_bins, const ColumnVector& plnp, 
		    int *jointhist, int *marghist1, int *marghist2,
		    float& jointentropy, float& margentropy1,
		    float& margentropy2)
    {
      // the joint and marginal entropies between the two images are
      //  calculated here and returned
      // the last parameter, plnp, is a vector containing values of -p*log(p)
      //  which makes the calculation significantly more efficient

      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1

      Matrix iaffbig = vtest.sampling_mat().i() * aff.i() *
	                     vref.sampling_mat();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.xsize()-1, yb1=vref.ysize()-1, zb1=vref.zsize()-1;
      float  xb2 = ((float) vtest.xsize())-1.0001,
	yb2=((float) vtest.ysize())-1.0001, zb2=((float) vtest.zsize())-1.0001;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      for (long int i=0; i<((no_bins+1)*(no_bins+1)); i++) {
	  jointhist[i]=0;
      }
      for (int i=0; i<=no_bins; i++) {
	  marghist1[i]=0;
	  marghist2[i]=0;
      }

      long int a,b;
      float b1=no_bins/(maxtest-mintest), b0=-mintest*no_bins/(maxtest-mintest);
      float val=0.0;
 
      unsigned int xmin, xmax;
      int *bptr;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  bptr = get_bindexptr(xmin,y,z,vref,bindex);

	  for (unsigned int x=xmin; x<=xmax; x++) {

	    if ( !((x==xmin) || (x==xmax)) 
		 || in_interp_bounds(vtest,o1,o2,o3) )
	      {
		val = q_tri_interpolation(vtest,o1,o2,o3);
		
		// do the cost function record keeping...
		a=*bptr;
		b=(long int) (val*b1 + b0);
		if (b>=no_bins) b=no_bins-1;
		if (b<0) b=0;
		(jointhist[a*(no_bins+1) + b])++;
		(marghist1[a])++;
		(marghist2[b])++;
	      }

	    bptr++;
	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	}
      }

      // note that the plnp values indexed by integers such that:
      //   plnp(n) = n/N * log(n/N)
      float p=0.0;
      int n=0, psize=plnp.Nrows();
      long int nvoxels = (long int) (vref.ysize() * vref.xsize() * 
				     vref.zsize());
      for (long int i=0; i<((no_bins+1)*(no_bins+1)); i++) {
	n = jointhist[i];
	if (n>0) {
	  if (n<=psize)
	    jointentropy+=plnp(n);
	  else {
	    p = ((float) n) / ((float) nvoxels);
	    jointentropy+= - p*log(p);
	  }
	}
      }
      for (int i=0; i<=no_bins; i++) {
	n = marghist1[i];
	if (n>0) {
	  if (n<=psize)
	    margentropy1+=plnp(n);
	  else {
	    p = ((float) n) / ((float) nvoxels);
	    margentropy1+= - p*log(p);
	  }
	}
      }
      long int noverlap=0;
      for (int i=0; i<=no_bins; i++) {
	n = marghist2[i];
	if (n>0) {
	  noverlap += n;
	  if (n<=psize)
	    margentropy2+=plnp(n);
	  else {
	    //cerr << ":";
	    p = ((float) n) / ((float) nvoxels);
	    margentropy2+= - p*log(p);
	  }
	}
      }

      // correct for difference in total histogram size
      //  that is: noverlap vs nvoxels
      // H_1 = N_0/N_1 * H_0 + log(N_1/N_0)
      //     = N_0/N_1 * H_0 - log(N_0/N_1)
      if (noverlap > 0) {
	float nratio = ((float) nvoxels) / ((float) noverlap);
	jointentropy = nratio * jointentropy - log(nratio);
	margentropy1 = nratio * margentropy1 - log(nratio);
	margentropy2 = nratio * margentropy2 - log(nratio);
      } else {
	// Put in maximum entropy values as base cases = BAD registration
	jointentropy = 2.0*log(no_bins);
	margentropy1 = log(no_bins);
	margentropy2 = log(no_bins);
      }
      return;
    }



  float p_test_entropy(const volume<float>& vref, const volume<float>& vtest,
		    int *bindex, const Matrix& aff,
		    const float mintest, const float maxtest,
		    const int no_bins, const ColumnVector& plnp, 
		    int *jointhist, int *marghist1, int *marghist2)
    {
      float jointentropy=0.0, margentropy1=0.0, margentropy2=0.0;
      calc_entropy(vref,vtest,bindex,aff,mintest,maxtest,no_bins,
		   plnp,jointhist,marghist1,marghist2,
		   jointentropy,margentropy1,margentropy2);
      return margentropy1;
    }

  float p_ref_entropy(const volume<float>& vref, const volume<float>& vtest,
		    int *bindex, const Matrix& aff,
		    const float mintest, const float maxtest,
		    const int no_bins, const ColumnVector& plnp, 
		    int *jointhist, int *marghist1, int *marghist2)
    {
      float jointentropy=0.0, margentropy1=0.0, margentropy2=0.0;
      calc_entropy(vref,vtest,bindex,aff,mintest,maxtest,no_bins,
		   plnp,jointhist,marghist1,marghist2,
		   jointentropy,margentropy1,margentropy2);
      return margentropy2;
    }

  float p_joint_entropy(const volume<float>& vref, const volume<float>& vtest,
		    int *bindex, const Matrix& aff,
		    const float mintest, const float maxtest,
		    const int no_bins, const ColumnVector& plnp, 
		    int *jointhist, int *marghist1, int *marghist2)
    {
      float jointentropy=0.0, margentropy1=0.0, margentropy2=0.0;
      calc_entropy(vref,vtest,bindex,aff,mintest,maxtest,no_bins,
		   plnp,jointhist,marghist1,marghist2,
		   jointentropy,margentropy1,margentropy2);
      return jointentropy;
    }

  float p_mutual_info(const volume<float>& vref, const volume<float>& vtest,
		    int *bindex, const Matrix& aff,
		    const float mintest, const float maxtest,
		    const int no_bins, const ColumnVector& plnp, 
		    int *jointhist, int *marghist1, int *marghist2)
    {
      float jointentropy=0.0, margentropy1=0.0, margentropy2=0.0;
      calc_entropy(vref,vtest,bindex,aff,mintest,maxtest,no_bins,
		   plnp,jointhist,marghist1,marghist2,
		   jointentropy,margentropy1,margentropy2);
      float mutualinformation = margentropy1 + margentropy2 - jointentropy;
      return mutualinformation;
    }



  float p_normalised_mutual_info(const volume<float>& vref, 
			       const volume<float>& vtest,
			       int *bindex, const Matrix& aff,
			       const float mintest, const float maxtest,
			       const int no_bins, const ColumnVector& plnp, 
			       int *jointhist, int *marghist1, int *marghist2)
    {
      float jointentropy=0.0, margentropy1=0.0, margentropy2=0.0;
      calc_entropy(vref,vtest,bindex,aff,mintest,maxtest,no_bins,
		   plnp,jointhist,marghist1,marghist2,
		   jointentropy,margentropy1,margentropy2);
      float normmi;
      if (fabs(jointentropy)<1e-9) {
	normmi = 0.0;  // BAD registration result
      } else {
	normmi = (margentropy1 + margentropy2)/jointentropy;
      }
      return normmi;
    }


  ///////////////////////////////////////////////////////////////////////

   void calc_smoothed_entropy(const volume<float>& vref, 
			      const volume<float>& vtest,
			      int *bindex,  const Matrix& aff,
			      const float mintest, const float maxtest,
			      const int no_bins,
			      float *jointhist, float *marghist1, 
			      float *marghist2,
			      float& jointentropy, float& margentropy1,
			      float& margentropy2,
			      const float smoothsize, const float fuzzyfrac)
   {
      // the joint and marginal entropies between the two images are
      //  calculated here and returned

      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1

      Matrix iaffbig = vtest.sampling_mat().i() * aff.i() *
	                     vref.sampling_mat();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.xsize()-1, yb1=vref.ysize()-1, zb1=vref.zsize()-1;
      float  xb2 = ((float) vtest.xsize())-1.0001,
	yb2=((float) vtest.ysize())-1.0001, zb2=((float) vtest.zsize())-1.0001;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      for (long int i=0; i<((no_bins+1)*(no_bins+1)); i++) {
	  jointhist[i]=0;
      }
      for (int i=0; i<=no_bins; i++) {
	  marghist1[i]=0;
	  marghist2[i]=0;
      }

      long int a;
      float b1=no_bins/(maxtest-mintest), b0=-mintest*no_bins/(maxtest-mintest);
      float val=0.0;
 
      float smoothx, smoothy, smoothz, geomweight, wcentre, wplus, wminus, bidx;
      long int bcentre, bplus, bminus;
      smoothx = smoothsize / vtest.xdim();
      smoothy = smoothsize / vtest.ydim();
      smoothz = smoothsize / vtest.zdim();

      unsigned int xmin, xmax;
      int *bptr;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  bptr = get_bindexptr(xmin,y,z,vref,bindex);

	  for (unsigned int x=xmin; x<=xmax; x++) {

	    if ( !((x==xmin) || (x==xmax)) 
		 || in_interp_bounds(vtest,o1,o2,o3) )
	      {
		val = q_tri_interpolation(vtest,o1,o2,o3);
		
		// do the cost function record keeping...
		geomweight=1.0;
		if (o1<smoothx)  geomweight*=o1/smoothx;
		else if ((xb2-o1)<smoothx) geomweight*=(xb2-o1)/smoothx;
		if (o2<smoothy)  geomweight*=o2/smoothy;
		else if ((yb2-o2)<smoothy) geomweight*=(yb2-o2)/smoothy;
		if (o3<smoothz)  geomweight*=o3/smoothz;
		else if ((zb2-o3)<smoothz) geomweight*=(zb2-o3)/smoothz;
		if (geomweight<0.0)  geomweight=0.0;
		
		// do the cost function record keeping...
		a=*bptr;
		bidx=val*b1 + b0;
		bcentre=(long int) (bidx);
		bplus = bcentre + 1;
		bminus = bcentre - 1;
		if (bcentre>=no_bins) {
		  bcentre=no_bins-1;
		  bplus = bcentre;
		}
		if (bcentre<0) {
		  bcentre=0;
		  bminus = bcentre;
		}
		if (bplus>=no_bins) { bplus = no_bins-1; }
		if (bminus<0) { bminus = 0; }
		// Fuzzy binning weights
		bidx = fabs(bidx - (int) bidx);  // get fractional component : [0,1]
		if (bidx<fuzzyfrac) {
		  wcentre = 0.5 + 0.5*(bidx/fuzzyfrac);
		  wminus = 1 - wcentre;
		  wplus = 0;
		} else if (bidx>(1.0-fuzzyfrac)) {
		  wcentre = 0.5 + 0.5*((1.0-bidx)/fuzzyfrac);
		  wplus = 1 - wcentre;
		  wminus=0;
		} else {
		  wcentre = 1;
		  wplus =0;
		  wminus=0;
		}
		(jointhist[a*(no_bins+1) + bcentre])+=geomweight*wcentre;
		(marghist2[bcentre])+=geomweight*wcentre;
		(jointhist[a*(no_bins+1) + bplus])+=geomweight*wplus;
		(marghist2[bplus])+=geomweight*wplus;
		(jointhist[a*(no_bins+1) + bminus])+=geomweight*wminus;
		(marghist2[bminus])+=geomweight*wminus;
		(marghist1[a])+=geomweight;
	      }

	    bptr++;
	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	}
      }

      float p=0.0, n=0.0;
      long int nvoxels = (long int) (vref.ysize() * vref.xsize() * 
				     vref.zsize());
      for (long int i=0; i<((no_bins+1)*(no_bins+1)); i++) {
	n = jointhist[i];
	if (n>0) {
	  p = n / ((float) nvoxels);
	  jointentropy+= - p*log(p);
	}
      }
      for (int i=0; i<=no_bins; i++) {
	n = marghist1[i];
	if (n>0) {
	  p = n / ((float) nvoxels);
	  margentropy1+= - p*log(p);
	}
      }
      float noverlap=0;
      for (int i=0; i<=no_bins; i++) {
	n = marghist2[i];
	if (n>0) {
	  noverlap += n;
	  p = n / ((float) nvoxels);
	  margentropy2+= - p*log(p);
	}
      }

      // correct for difference in total histogram size
      //  that is: noverlap vs nvoxels
      // H_1 = N_0/N_1 * H_0 + log(N_1/N_0)
      //     = N_0/N_1 * H_0 - log(N_0/N_1)
      if (noverlap > 0) {
	float nratio = ((float) nvoxels) / ((float) noverlap);
	jointentropy = nratio * jointentropy - log(nratio);
	margentropy1 = nratio * margentropy1 - log(nratio);
	margentropy2 = nratio * margentropy2 - log(nratio);
      } else {
	// Put in maximum entropy values as base cases = BAD registration
	jointentropy = 2.0*log(no_bins);
	margentropy1 = log(no_bins);
	margentropy2 = log(no_bins);
      }
      return;
    }



  float p_mutual_info_smoothed(const volume<float>& vref, 
			     const volume<float>& vtest,
			     int *bindex, const Matrix& aff,
			     const float mintest, const float maxtest,
			     const int no_bins, 
			     float *jointhist, float *marghist1, 
			     float *marghist2,
			     const float smoothsize, const float fuzzyfrac)
    {
      float jointentropy=0.0, margentropy1=0.0, margentropy2=0.0;
      calc_smoothed_entropy(vref,vtest,bindex,aff,mintest,maxtest,no_bins,
			    jointhist,marghist1,marghist2,
			    jointentropy,margentropy1,margentropy2,
			    smoothsize,fuzzyfrac);
      float mutualinformation = margentropy1 + margentropy2 - jointentropy;
      return mutualinformation;
    }



  float p_normalised_mutual_info_smoothed(const volume<float>& vref, const volume<float>& vtest,
			       int *bindex, const Matrix& aff,
			       const float mintest, const float maxtest,
			       const int no_bins, float *jointhist, 
			       float *marghist1, float *marghist2,
			       const float smoothsize, const float fuzzyfrac)
    {
      float jointentropy=0.0, margentropy1=0.0, margentropy2=0.0;
      calc_smoothed_entropy(vref,vtest,bindex,aff,mintest,maxtest,no_bins,
			    jointhist,marghist1,marghist2,
			    jointentropy,margentropy1,margentropy2,
			    smoothsize,fuzzyfrac);
      float normmi;
      if (fabs(jointentropy)<1e-9) {
	normmi = 0.0;  // BAD registration result
      } else {
	normmi = (margentropy1 + margentropy2)/jointentropy;
      }
      return normmi;
    }

  ///////////////////////////////////////////////////////////////////////

   void calc_fully_weighted_entropy(const volume<float>& vref, 
				    const volume<float>& vtest,
				    const volume<float>& refweight, 
				    const volume<float>& testweight,
				    int *bindex,  const Matrix& aff,
				    const float mintest, const float maxtest,
				    const int no_bins,
				    float *jointhist, float *marghist1, 
				    float *marghist2,
				    float& jointentropy, float& margentropy1,
				    float& margentropy2,
				    const float smoothsize, 
				    const float fuzzyfrac)
   {
      // the joint and marginal entropies between the two images are
      //  calculated here and returned

      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1

      Matrix iaffbig = vtest.sampling_mat().i() * aff.i() *
	                     vref.sampling_mat();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.xsize()-1, yb1=vref.ysize()-1, zb1=vref.zsize()-1;
      float  xb2 = ((float) vtest.xsize())-1.0001,
	yb2=((float) vtest.ysize())-1.0001, zb2=((float) vtest.zsize())-1.0001;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      for (long int i=0; i<((no_bins+1)*(no_bins+1)); i++) {
	  jointhist[i]=0;
      }
      for (int i=0; i<=no_bins; i++) {
	  marghist1[i]=0;
	  marghist2[i]=0;
      }

      long int a;
      float b1=no_bins/(maxtest-mintest), b0=-mintest*no_bins/(maxtest-mintest);
      float val=0.0;
 
      float wval=0, smoothx, smoothy, smoothz;
      float geomweight, wcentre, wplus, wminus, bidx;
      long int bcentre, bplus, bminus;
      smoothx = smoothsize / vtest.xdim();
      smoothy = smoothsize / vtest.ydim();
      smoothz = smoothsize / vtest.zdim();

      unsigned int xmin, xmax;
      int *bptr;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  bptr = get_bindexptr(xmin,y,z,vref,bindex);

	  for (unsigned int x=xmin; x<=xmax; x++) {

	    if ( !((x==xmin) || (x==xmax)) 
		 || in_interp_bounds(vtest,o1,o2,o3) )
	      {
		val = q_tri_interpolation(vtest,o1,o2,o3);
		wval = q_tri_interpolation(testweight,o1,o2,o3);
		
		// do the cost function record keeping...
		geomweight=wval*refweight(x,y,z);
		if (o1<smoothx)  geomweight*=o1/smoothx;
		else if ((xb2-o1)<smoothx) geomweight*=(xb2-o1)/smoothx;
		if (o2<smoothy)  geomweight*=o2/smoothy;
		else if ((yb2-o2)<smoothy) geomweight*=(yb2-o2)/smoothy;
		if (o3<smoothz)  geomweight*=o3/smoothz;
		else if ((zb2-o3)<smoothz) geomweight*=(zb2-o3)/smoothz;
		if (geomweight<0.0)  geomweight=0.0;
		
		// do the cost function record keeping...
		a=*bptr;
		bidx=val*b1 + b0;
		bcentre=(long int) (bidx);
		bplus = bcentre + 1;
		bminus = bcentre - 1;
		if (bcentre>=no_bins) {
		  bcentre=no_bins-1;
		  bplus = bcentre;
		}
		if (bcentre<0) {
		  bcentre=0;
		  bminus = bcentre;
		}
		if (bplus>=no_bins) { bplus = no_bins-1; }
		if (bminus<0) { bminus = 0; }
		// Fuzzy binning weights
		bidx = fabs(bidx - (int) bidx);  // get fractional component : [0,1]
		if (bidx<fuzzyfrac) {
		  wcentre = 0.5 + 0.5*(bidx/fuzzyfrac);
		  wminus = 1 - wcentre;
		  wplus = 0;
		} else if (bidx>(1.0-fuzzyfrac)) {
		  wcentre = 0.5 + 0.5*((1.0-bidx)/fuzzyfrac);
		  wplus = 1 - wcentre;
		  wminus=0;
		} else {
		  wcentre = 1;
		  wplus =0;
		  wminus=0;
		}
		(jointhist[a*(no_bins+1) + bcentre])+=geomweight*wcentre;
		(marghist2[bcentre])+=geomweight*wcentre;
		(jointhist[a*(no_bins+1) + bplus])+=geomweight*wplus;
		(marghist2[bplus])+=geomweight*wplus;
		(jointhist[a*(no_bins+1) + bminus])+=geomweight*wminus;
		(marghist2[bminus])+=geomweight*wminus;
		(marghist1[a])+=geomweight;
	      }

	    bptr++;
	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	}
      }

      float p=0.0, n=0.0;
      long int nvoxels = (long int) (vref.ysize() * vref.xsize() * 
				     vref.zsize());
      for (long int i=0; i<((no_bins+1)*(no_bins+1)); i++) {
	n = jointhist[i];
	if (n>0) {
	  p = n / ((float) nvoxels);
	  jointentropy+= - p*log(p);
	}
      }
      for (int i=0; i<=no_bins; i++) {
	n = marghist1[i];
	if (n>0) {
	  p = n / ((float) nvoxels);
	  margentropy1+= - p*log(p);
	}
      }
      float noverlap=0;
      for (int i=0; i<=no_bins; i++) {
	n = marghist2[i];
	if (n>0) {
	  noverlap += n;
	  p = n / ((float) nvoxels);
	  margentropy2+= - p*log(p);
	}
      }

      // correct for difference in total histogram size
      //  that is: noverlap vs nvoxels
      // H_1 = N_0/N_1 * H_0 + log(N_1/N_0)
      //     = N_0/N_1 * H_0 - log(N_0/N_1)
      if (noverlap > 0) {
	float nratio = ((float) nvoxels) / ((float) noverlap);
	jointentropy = nratio * jointentropy - log(nratio);
	margentropy1 = nratio * margentropy1 - log(nratio);
	margentropy2 = nratio * margentropy2 - log(nratio);
      } else {
	// Put in maximum entropy values as base cases = BAD registration
	jointentropy = 2.0*log(no_bins);
	margentropy1 = log(no_bins);
	margentropy2 = log(no_bins);
      }
      return;
    }



  float p_mutual_info_fully_weighted(const volume<float>& vref, 
				   const volume<float>& vtest,
				    const volume<float>& refweight, 
				    const volume<float>& testweight,
			     int *bindex, const Matrix& aff,
			     const float mintest, const float maxtest,
			     const int no_bins, 
			     float *jointhist, float *marghist1, 
			     float *marghist2,
			     const float smoothsize, const float fuzzyfrac)
    {
      float jointentropy=0.0, margentropy1=0.0, margentropy2=0.0;
      calc_fully_weighted_entropy(vref,vtest,refweight,testweight,
				  bindex,aff,mintest,maxtest,no_bins,
			    jointhist,marghist1,marghist2,
			    jointentropy,margentropy1,margentropy2,
			    smoothsize,fuzzyfrac);
      float mutualinformation = margentropy1 + margentropy2 - jointentropy;
      return mutualinformation;
    }



  float p_normalised_mutual_info_fully_weighted(const volume<float>& vref, 
					      const volume<float>& vtest,
					      const volume<float>& refweight, 
					      const volume<float>& testweight,
			       int *bindex, const Matrix& aff,
			       const float mintest, const float maxtest,
			       const int no_bins, float *jointhist, 
			       float *marghist1, float *marghist2,
			       const float smoothsize, const float fuzzyfrac)
    {
      float jointentropy=0.0, margentropy1=0.0, margentropy2=0.0;
      calc_fully_weighted_entropy(vref,vtest,refweight,testweight,
				  bindex,aff,mintest,maxtest,no_bins,
			    jointhist,marghist1,marghist2,
			    jointentropy,margentropy1,margentropy2,
			    smoothsize,fuzzyfrac);
      float normmi;
      if (fabs(jointentropy)<1e-9) {
	normmi = 0.0;  // BAD registration result
      } else {
	normmi = (margentropy1 + margentropy2)/jointentropy;
      }
      return normmi;
    }


  ///////////////////////////////////////////////////////////////////////

  void Costfn::set_debug_mode(bool debug_flag) 
  {
    debug_mode = debug_flag;
  }


  void Costfn::set_no_bins(int n_bins) 
    {
      no_bins=n_bins;
      //hist.ReSize(no_bins+1,no_bins+1);
      jointhist = new int[(no_bins+1)*(no_bins+1)];
      marghist1 = new int[no_bins+1];
      marghist2 = new int[no_bins+1];
      fjointhist = new float[(no_bins+1)*(no_bins+1)];
      fmarghist1 = new float[no_bins+1];
      fmarghist2 = new float[no_bins+1];
      unsigned long int N = this->refvol.nvoxels();
      float p=0.0;
      try {
        plnp.ReSize(Min((unsigned long int) 10000, (unsigned long int) (10*N/(no_bins+1))));
      } catch(...) { 
	cerr<<"ERROR: failed on plnp.Resize("
	    <<Min((unsigned long int) 10000, (unsigned long int) (10*N/(no_bins+1)))
	    <<")"<<endl;
	cerr <<"Probably out of memory" << endl;
      }
      for (int num=1; num<=plnp.Nrows(); num++) {
	p = ((float) num) / ((float) N);
	plnp(num) = -p*log(p);
      }

      if (bindex) delete [] bindex;
      bindex = new int[refvol.nvoxels()];
      float refmin, refmax;
      refmin = this->refvol.min();
      refmax = this->refvol.max();
      if (refmax-refmin==0.0) refmax+=1;
      bin_a1=no_bins/(refmax-refmin);
      bin_a0=-refmin*no_bins/(refmax-refmin);
      int *bptr = bindex;
      for (int z=0; z<refvol.zsize(); z++) { 
	for (int y=0; y<refvol.ysize(); y++) { 
	  for (int x=0; x<refvol.xsize(); x++) {
	    *bptr = (int) get_bin_number(refvol(x,y,z));
	    if (*bptr >= no_bins)  *bptr = no_bins-1;
	    if (*bptr < 0.0)  *bptr = 0;
	    bptr++;
	  }
	}
      }
    }

  float Costfn::get_bin_intensity(int bin_number) const
  {
    return (bin_number +0.5 - bin_a0)/bin_a1;
  }

  float Costfn::get_bin_number(float intensity) const
  {
    return (intensity*bin_a1 + bin_a0);
  }

  bool Costfn::is_bbr_set(void) const
  {
    if (no_coords>0) { return true; }
    return false;
  }

  int Costfn::set_bbr_seg(const volume<float>& bbr_seg) 
    {
      this->wmseg = bbr_seg;
      // pre-calculate points on the interface and normals to these
      int npts=0;
      volume<float> wmsum, wmsmooth;
      volume4D<float> wmgrad;
      ColumnVector coord(4), box3(3);
      box3=1.0;
      wmsum=convolve_separable(wmseg,box3,box3,box3);
      wmsmooth=smooth(wmseg,2.0);  // 2.0mm smoothing sigma (is this good?) - just used for getting gradients/normals
      gradient(wmsmooth,wmgrad);
      //save_volume(wmsmooth,"wmsmooth");
      //save_volume(wmsum,"wmsum");
      //save_volume4D(wmgrad,"wmgrad");
      for (int z=wmseg.minz(); z<=wmseg.maxz(); z++) {
	for (int y=wmseg.miny(); y<=wmseg.maxy(); y++) {
	  for (int x=wmseg.minx(); x<=wmseg.maxx(); x++) {
	    // if this point is WM but some neighbours are not
	    if (wmseg(x,y,z)>0.5) {
	      if (wmsum(x,y,z)<26.5) { npts++; }
	    }
	  }
	}
      }
      //bbr_pts.ReSize(npts,3);
      //bbr_norms.ReSize(npts,3);
      no_coords = npts;
      // cerr << "Number of points on surface for BBR is " << no_coords << endl; // DEBUG DEBUG DEBUG
      gm_coord_x = new float[npts];
      gm_coord_y = new float[npts];
      gm_coord_z = new float[npts];
      wm_coord_x = new float[npts];
      wm_coord_y = new float[npts];
      wm_coord_z = new float[npts];
      npts=0;
      for (int z=wmseg.minz(); z<=wmseg.maxz(); z++) {
	for (int y=wmseg.miny(); y<=wmseg.maxy(); y++) {
	  for (int x=wmseg.minx(); x<=wmseg.maxx(); x++) {
	    if (wmseg(x,y,z)>0.5) {
	      if (wmsum(x,y,z)<26.5) { 
		coord(1)=x; coord(2)=y; coord(3)=z; coord(4)=1;  // vox coord
		coord = wmseg.sampling_mat() * coord;  // mm coord
		//bbr_pts(npts,1)=coord(1);
		//bbr_pts(npts,2)=coord(2);
		//bbr_pts(npts,3)=coord(3);
		// turn gradient into mm gradients
		float normx=-wmgrad(x,y,z,0)/refvol.xdim();
		float normy=-wmgrad(x,y,z,1)/refvol.ydim();
		float normz=-wmgrad(x,y,z,2)/refvol.zdim();
		float normnorm=std::sqrt(norm2sq(normx,normy,normz));
		normx /= normnorm;
		normy /= normnorm;
		normz /= normnorm;
		// stored coordinates are in mm
		gm_coord_x[npts] = coord(1) + normx * bbr_dist; 
		gm_coord_y[npts] = coord(2) + normy * bbr_dist; 
		gm_coord_z[npts] = coord(3) + normz * bbr_dist; 
		wm_coord_x[npts] = coord(1) - normx * bbr_dist;
		wm_coord_y[npts] = coord(2) - normy * bbr_dist; 
		wm_coord_z[npts] = coord(3) - normz * bbr_dist; 
		npts++; 
	      }
	    }
	  }
	}
      }
      if (no_coords==0) {
	cerr << "ERROR::set_bbr_seg: could not find any boundary points!" << endl;
	return 1;
      }
      return 0;
    }

  int Costfn::set_bbr_coords(const Matrix& coords, const Matrix& norms) 
    {
      //bbr_pts = coords;
      //bbr_norms = norms;
      if ( (coords.Nrows()==0) || (norms.Nrows()==0) ||
	   (coords.Nrows() != norms.Nrows()) ) {
	cerr << "ERROR::set_bbr_coords: coords and norms are different sizes or zero size" << endl;
	return 1;
      }
      no_coords = coords.Nrows();
      gm_coord_x = new float[no_coords];
      gm_coord_y = new float[no_coords];
      gm_coord_z = new float[no_coords];
      wm_coord_x = new float[no_coords];
      wm_coord_y = new float[no_coords];
      wm_coord_z = new float[no_coords];
      for (int npts=0; npts<no_coords; npts++) {
	gm_coord_x[npts] = coords(npts+1,1) + norms(npts+1,1) * bbr_dist; 
	gm_coord_y[npts] = coords(npts+1,2) + norms(npts+1,2) * bbr_dist; 
	gm_coord_z[npts] = coords(npts+1,3) + norms(npts+1,3) * bbr_dist; 
	wm_coord_x[npts] = coords(npts+1,1) - norms(npts+1,1) * bbr_dist;
	wm_coord_y[npts] = coords(npts+1,2) - norms(npts+1,2) * bbr_dist; 
	wm_coord_z[npts] = coords(npts+1,3) - norms(npts+1,3) * bbr_dist; 
      }
      return 0;
    }

  int Costfn::set_bbr_step(int step) 
  {
    vertex_step = step;
    return 0;
  }

  int Costfn::set_bbr_slope(float slope) 
  {
    bbr_slope = slope;
    return 0;
  }

  int Costfn::set_bbr_type(const string& typenm) 
  {
    if ((typenm=="signed") || (typenm=="local_abs") || (typenm=="global_abs")) {
      bbr_type=typenm;
    } else {
      imthrow("Unrecognised BBR type: " + typenm + "\nValid types are: signed, global_abs, local_abs",30);
    }
    return 0;
  }


  int Costfn::set_bbr_fmap(const volume<float>& fieldmap, int phase_encode_direction) 
  {
    fmap=fieldmap;
    fmap_mask=fmap*0.0f+1.0f;
    pe_dir=phase_encode_direction;  // TODO: some sanity checking on this first
    return 0;
  }

  int Costfn::set_bbr_fmap(const volume<float>& fieldmap, const volume<float>& fieldmap_mask, int phase_encode_direction)
  {
    fmap=fieldmap;
    fmap_mask=fieldmap_mask;   // TODO: check they are the same dimensions
    pe_dir=phase_encode_direction;  // TODO: some sanity checking on this first
    return 0;

  }


  // Member function interfaces


  float Costfn::normcorr(const Matrix& aff) const
    {
      p_count++;
      return p_normcorr(this->refvol,this->testvol,aff);
    }

  float Costfn::normcorr_smoothed(const Matrix& aff) const
    {
      p_count++;
      return p_normcorr_smoothed(this->refvol,this->testvol,aff, 
				 this->smoothsize);
    }


  float Costfn::normcorr_smoothed_sinc(const Matrix& aff) const
    {
      p_count++;
      return p_normcorr_smoothed_sinc(this->refvol,this->testvol,aff, 
				      this->smoothsize);
    }

  float Costfn::normcorr_fully_weighted(const Matrix& aff,
				const volume<float>& refweight, 
				const volume<float>& testweight) const
   {
     p_count++;
     return p_normcorr_fully_weighted(this->refvol,this->testvol,
				    refweight,testweight,
				    aff, this->smoothsize);
   }


  float Costfn::leastsquares(const Matrix& aff) const
    {
      p_count++;
     return p_leastsquares(this->refvol,this->testvol,aff);
    }

  float Costfn::leastsquares_smoothed(const Matrix& aff) const
    {
      p_count++;
      return p_leastsquares_smoothed(this->refvol,this->testvol,aff, 
				   this->smoothsize);
    }

   float Costfn::leastsquares_fully_weighted(const Matrix& aff, 
				     const volume<float>& refweight, 
				     const volume<float>& testweight) const
   {
     p_count++;
     return p_leastsquares_fully_weighted(this->refvol,this->testvol,
					refweight,testweight,
					aff,this->smoothsize);
   }


  float Costfn::labeldiff(const Matrix& aff) const
    {
      p_count++;
     return p_labeldiff(this->refvol,this->testvol,aff);
    }

  float Costfn::labeldiff_smoothed(const Matrix& aff) const
    {
      p_count++;
      return p_labeldiff_smoothed(this->refvol,this->testvol,aff, 
				   this->smoothsize);
    }

   float Costfn::labeldiff_fully_weighted(const Matrix& aff, 
				     const volume<float>& refweight, 
				     const volume<float>& testweight) const
   {
     p_count++;
     return p_labeldiff_fully_weighted(this->refvol,this->testvol,
					refweight,testweight,
					aff,this->smoothsize);
   }


  float Costfn::bbr(const Matrix& aff) const
    {
      ColumnVector unit_value(1);
      unit_value=1.0;
      return bbr(aff, unit_value);
    }


  float Costfn::woods_fn(const Matrix& aff) const
    {
      p_count++;
      return p_woods_fn(this->refvol,this->testvol,this->bindex,aff,
		      this->no_bins);
    }

  float Costfn::woods_fn_smoothed(const Matrix& aff) const
    {
      p_count++;
      return p_woods_fn_smoothed(this->refvol,this->testvol,this->bindex,aff,
		      this->no_bins, this->smoothsize);
    }


  float Costfn::corr_ratio(const Matrix& aff) const
    {
      p_count++;
      return p_corr_ratio(this->refvol,this->testvol,this->bindex,aff,
			this->no_bins);
    }
  
  float Costfn::corr_ratio_smoothed(const Matrix& aff) const
    {
      p_count++;
      return p_corr_ratio_smoothed(this->refvol,this->testvol,this->bindex,aff,
			this->no_bins, this->smoothsize);
    }
  
  float Costfn::corr_ratio_fully_weighted(const Matrix& aff,
				  const volume<float>& refweight, 
				  const volume<float>& testweight) const
    {
      p_count++;
      return p_corr_ratio_fully_weighted(this->refvol,this->testvol,
				       refweight,testweight,
				       this->bindex,aff,
				       this->no_bins, this->smoothsize);
    }
  
  float Costfn::corr_ratio_fully_weighted(const volume4D<float>& warpvol,
				  const volume<float>& refweight, 
				  const volume<float>& testweight) const
    {
      p_count++;
      return p_corr_ratio_fully_weighted(this->refvol,this->testvol,
				       refweight,testweight,
				       this->bindex,warpvol,
				       this->no_bins, this->smoothsize);
    }
  

  float Costfn::corr_ratio_gradient_fully_weighted(volume4D<float>& gradvec,
					  const volume4D<float>& warpvol,
					  const volume<float>& refweight, 
					  const volume<float>& testweight, 
						   bool nullbc) const
  {
      p_count++;
      return p_corr_ratio_gradient_fully_weighted(gradvec,
						  this->refvol,this->testvol,
						  refweight,testweight,
						  this->bindex,warpvol,
						  this->no_bins, 
						  this->smoothsize, nullbc);
    }
  

  float Costfn::ref_entropy(const Matrix& aff) const
    {
      p_count++;
      return p_ref_entropy(this->refvol,this->testvol,this->bindex,aff,
			 this->testvol.min(),this->testvol.max(),
			 this->no_bins,this->plnp,this->jointhist,
			 this->marghist1,this->marghist2);
    }

  float Costfn::test_entropy(const Matrix& aff) const
    {
      p_count++;
      return p_test_entropy(this->refvol,this->testvol,this->bindex,aff,
			 this->testvol.min(),this->testvol.max(),
			 this->no_bins,this->plnp,this->jointhist,
			 this->marghist1,this->marghist2);
    }

  float Costfn::joint_entropy(const Matrix& aff) const
    {
      p_count++;
      return p_joint_entropy(this->refvol,this->testvol,this->bindex,aff,
			 this->testvol.min(),this->testvol.max(),
			 this->no_bins,this->plnp,this->jointhist,
			 this->marghist1,this->marghist2);
    }


  float Costfn::mutual_info(const Matrix& aff) const
    {
      p_count++;
      return p_mutual_info(this->refvol,this->testvol,this->bindex,aff,
			 this->testvol.min(),this->testvol.max(),
			 this->no_bins,this->plnp,this->jointhist,
			 this->marghist1,this->marghist2);
    }


  float Costfn::mutual_info_smoothed(const Matrix& aff) const
    {
      p_count++;
      return p_mutual_info_smoothed(this->refvol,this->testvol,
				  this->bindex,aff,
				  this->testvol.min(),this->testvol.max(),
				  this->no_bins,this->fjointhist,
				  this->fmarghist1,this->fmarghist2,
				  this->smoothsize, this->fuzzyfrac);
    }


  float Costfn::mutual_info_fully_weighted(const Matrix& aff,
				   const volume<float>& refweight, 
				   const volume<float>& testweight) const
    {
      p_count++;
      return p_mutual_info_fully_weighted(this->refvol,this->testvol,
				  refweight,testweight,
				  this->bindex,aff,
				  this->testvol.min(),this->testvol.max(),
				  this->no_bins,this->fjointhist,
				  this->fmarghist1,this->fmarghist2,
				  this->smoothsize, this->fuzzyfrac);
    }


  float Costfn::normalised_mutual_info(const Matrix& aff) const
    {
      p_count++;
      return p_normalised_mutual_info(this->refvol,this->testvol,this->bindex,aff,
			 this->testvol.min(),this->testvol.max(),
			 this->no_bins,this->plnp,this->jointhist,
			 this->marghist1,this->marghist2);
    }

  float Costfn::normalised_mutual_info_smoothed(const Matrix& aff) const
    {
      p_count++;
      return p_normalised_mutual_info_smoothed(this->refvol,this->testvol,
					     this->bindex,aff,
					     this->testvol.min(),
					     this->testvol.max(),
					     this->no_bins,this->fjointhist,
					     this->fmarghist1,this->fmarghist2,
					     this->smoothsize, this->fuzzyfrac);
    }

  float Costfn::normalised_mutual_info_fully_weighted(const Matrix& aff,
					      const volume<float>& refweight, 
					      const volume<float>& testweight) const
    {
      p_count++;
      return p_normalised_mutual_info_fully_weighted(this->refvol,this->testvol,
						   refweight,testweight,
						   this->bindex,aff,
						   this->testvol.min(),
						   this->testvol.max(),
						   this->no_bins,
						   this->fjointhist,
						   this->fmarghist1,
						   this->fmarghist2,
						   this->smoothsize, 
						   this->fuzzyfrac);
    }


  ///////////////////////////////////////////////////////////////////////////


  int Costfn::p_corr_ratio_image_mapper(volume<float>& vout,
				Matrix& mappingfn,
				const volume<float>& vref, 
				const volume<float>& vtest,
				const volume<float>& refweight, 
				const volume<float>& testweight,
				int *bindex, const Matrix& aff,
				const int no_bins, const float smoothsize) const
    {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1

      // NB: only maps the image if it has non-zero size on entry

      bool mapimage=true;
      if (vout.nvoxels()==0) {
	mapimage=false;
      }

      Matrix iaffbig = vtest.sampling_mat().i() * aff.i() *
	                     vref.sampling_mat();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.xsize()-1, yb1=vref.ysize()-1, zb1=vref.zsize()-1;
      float  xb2 = ((float) vtest.xsize())-1.0001,
	yb2=((float) vtest.ysize())-1.0001, zb2=((float) vtest.zsize())-1.0001;

      float *sumy, *sumy2;
      sumy = new float[no_bins+1];
      sumy2 = new float[no_bins+1];
      float *numy;
      numy = new float[no_bins+1];
      int b=0;
 
      for (int i=0; i<=no_bins; i++) {
	numy[i]=0.0; sumy[i]=0.0;  sumy2[i]=0.0;
      }

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4);
      float wval=0,val=0,o1,o2,o3;

      float smoothx, smoothy, smoothz, weight;
      smoothx = smoothsize / vtest.xdim();
      smoothy = smoothsize / vtest.ydim();
      smoothz = smoothsize / vtest.zdim();

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      unsigned int xmin, xmax;
      int *bptr;

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  // assume that this is always OK
	  bptr = get_bindexptr(xmin,y,z,vref,bindex);

	  for (unsigned int x=xmin; x<=xmax; x++) {

	    if ( !((x==xmin) || (x==xmax)) 
		 || in_interp_bounds(vtest,o1,o2,o3) )
	      {
		val = q_tri_interpolation(vtest,o1,o2,o3);
		wval = q_tri_interpolation(testweight,o1,o2,o3);
		
		// do the cost function record keeping...
		b=*bptr;
		weight=wval*refweight(x,y,z);
		if (o1<smoothx)  weight*=o1/smoothx;
		else if ((xb2-o1)<smoothx) weight*=(xb2-o1)/smoothx;
		if (o2<smoothy)  weight*=o2/smoothy;
		else if ((yb2-o2)<smoothy) weight*=(yb2-o2)/smoothy;
		if (o3<smoothz)  weight*=o3/smoothz;
		else if ((zb2-o3)<smoothz) weight*=(zb2-o3)/smoothz;
		if (weight<0.0)  weight=0.0;
		numy[b]+=weight;
		sumy[b]+=weight*val;
		sumy2[b]+=weight*val*val;
	      }

	    bptr++;
	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	}
      }

      // correct for occasion lapses into the last bin
      numy[no_bins-1] += numy[no_bins];
      sumy[no_bins-1] += sumy[no_bins];
      sumy2[no_bins-1] += sumy2[no_bins];
      numy[no_bins]=0.0;
      sumy[no_bins]=0.0;
      sumy2[no_bins]=0.0;

      if (mapimage) {
	for (unsigned int z=0; z<=zb1; z++) { 
	  for (unsigned int y=0; y<=yb1; y++) { 
	    xmin=0; xmax=xb1;
	    bptr = get_bindexptr(xmin,y,z,vref,bindex);
	    for (unsigned int x=xmin; x<=xmax; x++) {
	      b=*bptr;
	      if (numy[b]>0.01) {
		vout(x,y,z) = sumy[b]/numy[b];
	      } else {
		vout(x,y,z) = 0.0;
	      }
	      bptr++;
	    }
	  }
	}
      }

     // Generate mapping function
     mappingfn.ReSize(no_bins,2);
     for (int b=0; b<=no_bins-1; b++) {
       mappingfn(b+1,1) = get_bin_intensity(b);
   	if (numy[b]>0) {
	  mappingfn(b+1,2) = sumy[b]/numy[b];
	} else { 
	  mappingfn(b+1,2) = 0;
	}
     }

     delete [] numy; delete [] sumy; delete [] sumy2;
     
     return 0;
     
    }


  volume<float> Costfn::image_mapper(const Matrix& affmat) const // affmat is voxel to voxel
  {
    volume<float> vnew(this->refvol);
    Matrix mappingfn;
    p_corr_ratio_image_mapper(vnew,mappingfn,
			      this->refvol,this->testvol,
			      this->rweight,
			      this->tweight,
			      this->bindex,affmat,
			      this->no_bins,
			      this->smoothsize);
    return vnew;
  }

  Matrix Costfn::mappingfn(const Matrix& affmat) const // affmat is voxel to voxel
  {
    volume<float> vnew;
    Matrix mappingmat;
    p_corr_ratio_image_mapper(vnew,mappingmat,
			      this->refvol,this->testvol,
			      this->rweight,
			      this->tweight,
			      this->bindex,affmat,
			      this->no_bins,
			      this->smoothsize);
    return mappingmat;
  }

}



