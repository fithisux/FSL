/*  b0calc.cc

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


// B0 Calculation program - intended for use with the Virtual Scanner


#include <iostream>
#include <string>
#include "newimage/newimageall.h"
#include "utils/options.h"

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

using namespace Utilities;  
using namespace NEWIMAGE;

// Global constants

string title="b0calc\nB0 field calculation program\nCopyright(c) 2001, University of Oxford (Mark Jenkinson)";
string examples="b0calc -i <input> -o <output> [options]";

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<bool> calcxyz(string("--xyz"), false,
		  string("calculate and save all 3 field components (i.e. x,y,z)"),
		  false, no_argument);
Option<bool> direct_conv(string("--directconv"), false,
		  string("use direct (image space) convolution, not FFT"),
		  false, no_argument);
Option<float> delta(string("-d"), -9.45e-6,   // chi_tissue - chi_air
		  string("Delta value (chi_tissue - chi_air): default=-9.45e-6"),
		  false, requires_argument);
Option<float> chi0(string("--chi0"), +4e-7,   // chi of air
		  string("Value for susceptibility of air: default=+4e-7"),
		  false, requires_argument);
Option<float> b0z(string("--b0"), 1.0,
		  string("Value for zeroth-order b0 field (z-component): default=1"),
		  false, requires_argument);
Option<float> b0x(string("--b0x"), 0.0,
		  string("Value for zeroth-order b0 field (x-component): default=0"),
		  false, requires_argument);
Option<float> b0y(string("--b0y"), 0.0,
		  string("Value for zeroth-order b0 field (y-component): default=0"),
		  false, requires_argument);
Option<float> gx(string("--gx"), 0.0,
		  string("Value for zeroth-order x-gradient field (per mm): default=0"),
		  false, requires_argument);
Option<float> gy(string("--gy"), 0.0,
		  string("Value for zeroth-order y-gradient field (per mm): default=0"),
		  false, requires_argument);
Option<float> gz(string("--gz"), 0.0,
		  string("Value for zeroth-order z-gradient field (per mm): default=0"),
		  false, requires_argument);
Option<float> extendboundary(string("--extendboundary"), 1.0,
		  string("Relative proportion to extend voxels at boundary: default=1"),
		  false, requires_argument);
Option<string> inname(string("-i,--in"), string(""),
			 string("filename of input image (usually a tissue/air segmentation)"),
			 true, requires_argument);
Option<string> outname(string("-o,--out"), string(""),
			 string("filename of B0 output volume"),
			 true, requires_argument);

/////////////////////////////////////////////////////////////////////////////

// % Sphere simulation (including Lorentz Correction)
// % cx=100; cy=100; cz=100;
// % a = 60;
// % [xs,ys,zs]=ndgrid((-(cx-1)):(cx-1),(-(cy-1)):(cy-1),(-(cz-1)):(cz-1));
// % rs=sqrt(xs.^2+ys.^2+zs.^2);
// % mu=(rs<=a);
// % mus=(round(rs)==a);      % surface voxels
// % musdotn=mus.*zs./max(rs,1e-8);         % surface normals (in z)
// % bztheory = (1-mu)/3.*(a^3).*(3*zs.^2 - rs.^2)./max((rs.^5) + mu,1e-8);

/////////////////////////////////////////////////////////////////////////////

float arctan(float num, float den)
{
  if (fabs(den)>1e-8) { 
    return atan(num/den);
  } else {
    return atan2(num,den);
  }
}


float arcsinh(float num, float den)
{
  if (fabs(den)>1e-8) { 
    return asinh(num/den);
  } else {
    return log(num + sqrt(num*num + den*den)) - log(den);
    // uses asinh(x) = log(x + sqrt(x*x+1))
  }
}

volume<float> partbzkernel(const volume<float>& obj, 
			   float xoff, float yoff, float zoff,
			   const string& b0type)
{
  // This function returns F(x';x) as described in the technical report
  //  where xoff = x' - x , yoff = y' - y , etc
  int nx, ny, nz;
  nx = obj.xsize();  ny = obj.ysize();  nz = obj.zsize();
  float dx, dy, dz;
  dx = obj.xdim();  dy = obj.ydim();  dz = obj.zdim();
  volume<float> kernel(2*nx-1,2*ny-1,2*nz-1);
  for (int z1=-(nz-1); z1<=(nz-1); z1++) {
    for (int y1=-(ny-1); y1<=(ny-1); y1++) {
      for (int x1=-(nx-1); x1<=(nx-1); x1++) {
	// the following are x' in F(x';x)  (as always use F(x+offset;x))
	float x = x1*dx+xoff, y = y1*dy+yoff, z = z1*dz+zoff;
	float r = sqrt(x*x + y*y + z*z);
	if (b0type=="0,0,1") {
	  kernel(x1+nx-1,y1+ny-1,z1+nz-1) = arctan(x*y,z*r);
	} else if (b0type=="1,0,0") {
	  kernel(x1+nx-1,y1+ny-1,z1+nz-1) = -arcsinh(y,sqrt(x*x+z*z));
	} else if (b0type=="0,1,0") {
	  kernel(x1+nx-1,y1+ny-1,z1+nz-1) = -arcsinh(x,sqrt(y*y+z*z));
	} else if (b0type=="z,0,x") {
	  kernel(x1+nx-1,y1+ny-1,z1+nz-1) = (x-xoff)*arctan(x*y,z*r)
	    + x*arctan(z*y,x*r) + zoff*arcsinh(y,sqrt(x*x+z*z)) 
	    - y*arcsinh(z,sqrt(y*y+x*x));
	} else if (b0type=="0,z,y") {
	  kernel(x1+nx-1,y1+ny-1,z1+nz-1) = (y-yoff)*arctan(x*y,z*r)
	    + y*arctan(z*x,y*r) + zoff*arcsinh(x,sqrt(y*y+z*z)) 
	    - x*arcsinh(z,sqrt(y*y+x*x));
	} else if (b0type=="-x,0,z") {
	  kernel(x1+nx-1,y1+ny-1,z1+nz-1) = -xoff*arcsinh(y,sqrt(x*x+z*z))
	    -zoff*arctan(x*y,z*r);
	} else {
	  cerr << "Unknown b0 type requested:: " << b0type << endl;
	  exit(-1);
	}
      }
    }
  }
  cerr << "." << endl;
  // apply 1/(4*pi) factor - as it is more efficient here
  kernel *= 1.0/(4.0*M_PI);
  return kernel;
}


volume<float> calculate_kernel(const volume<float>& obj, 
			       const string& b0type="0,0,1")
{
  // Combines the partial voxel Bz solutions to get a full solution
  // Answer is for position x,y,z (relative to voxel centre) 
  //   given a voxel of size (dx,dy,dz)
  cout << "Calculating kernel" << endl;
  float dx2, dy2, dz2;
  dx2 = obj.xdim()/2.0;  dy2 = obj.ydim()/2.0;  dz2 = obj.zdim()/2.0;
  volume<float> kernel;
  kernel = partbzkernel(obj,dx2,dy2,dz2,b0type);
  kernel -= partbzkernel(obj,-dx2,dy2,dz2,b0type);
  kernel -= partbzkernel(obj,dx2,-dy2,dz2,b0type);
  kernel -= partbzkernel(obj,dx2,dy2,-dz2,b0type);
  kernel += partbzkernel(obj,dx2,-dy2,-dz2,b0type);
  kernel += partbzkernel(obj,-dx2,dy2,-dz2,b0type);
  kernel += partbzkernel(obj,-dx2,-dy2,dz2,b0type);
  kernel -= partbzkernel(obj,-dx2,-dy2,-dz2,b0type);
  return kernel;
}
  
/////////////////////////////////////////////////////////////////////////////

volume<float> gradientfield(const volume<float>& chi, const string& dir)
{
  int nx, ny, nz;
  nx = chi.xsize()/2;  ny = chi.ysize()/2;  nz = chi.zsize()/2;
  float dx, dy, dz;
  dx = chi.xdim();  dy = chi.ydim();  dz = chi.zdim();
  volume<float> gfield;
  gfield = chi * 0.0f;
  for (int z1=0; z1<chi.zsize(); z1++) {
    for (int y1=0; y1<chi.ysize(); y1++) {
      for (int x1=0; x1<chi.xsize(); x1++) {
	if (dir=="x") {
	  gfield(x1,y1,z1)=(x1-nx)*dx;
	} else if (dir=="y") {
	  gfield(x1,y1,z1)=(y1-ny)*dy;
	} else if (dir=="z") {
	  gfield(x1,y1,z1)=(z1-nz)*dz;
	} else {
	  cerr << "Unknown direction '"<<dir<<"' specified in gradientfield" 
	       << endl;
	  return gfield;
	}
      }
    }
  }
  return gfield;
}

/////////////////////////////////////////////////////////////////////////////


volume<float> convolve3_complex(const volume<float>& chi, 
				const volume<float>& kernel,
				complexvolume& kernelc)
{
  // assumes that kernel has size >= (2*nx-1,2*ny-1,2*nz-1)
  // also uses kernelc unless it has the wrong size
  int nx,ny,nz;
  nx = chi.xsize();  ny=chi.ysize();  nz=chi.zsize();
  complexvolume chi1c;
  cout << "zero pad" << endl;
  chi1c.re() = kernel*0.0f;
  chi1c.im() = chi1c.re();
  chi1c.re().setROIlimits(0,0,0,nx-1,ny-1,nz-1);
  chi1c.re().activateROI();
  chi1c.re().copyROIonly(chi);
  chi1c.re().deactivateROI();
  cout << "Forward FFT" << endl;
  fft3(chi1c);
  {
    if (!samesize(kernelc.re(),kernel)) {
      cout << "Forward FFT" << endl;
      kernelc.re() = kernel;
      kernelc.im() = kernel*0.0f;
      fft3(kernelc);
    }
    cout << "Kernel multiplication" << endl;
    {
      volume<float> tmp;
      tmp = chi1c.re() * kernelc.re() - chi1c.im() * kernelc.im();
      chi1c.im() = chi1c.re() * kernelc.im() + chi1c.im() * kernelc.re();
      chi1c.re() = tmp;
    } // destroy tmp
  } // destroy kernelc
  cout << "Inverse FFT" << endl;
  ifft3(chi1c);
  cout << "select ROI" << endl;
  chi1c.re().setROIlimits(nx-1,ny-1,nz-1,2*nx-2,2*ny-2,2*nz-2);
  chi1c.re().activateROI();
  volume<float> bz;
  bz = chi*0.0f;
  bz.copyROIonly(chi1c.re());
  return bz;
}


volume<float> direct_convolve(const volume<float>& source, 
			      const volume<float>& kernel)
{ 
  extrapolation oldex = source.getextrapolationmethod();
  if ((oldex==boundsassert) || (oldex==boundsexception)) 
    { source.setextrapolationmethod(constpad); }
  if (    (( (kernel.maxz() - kernel.minz()) % 2)==1) || 
	  (( (kernel.maxy() - kernel.miny()) % 2)==1) ||
	  (( (kernel.maxx() - kernel.minx()) % 2)==1) ) 
    {
      cerr << "WARNING:: Off-centre convolution being performed as kernel"
	   << " has even dimensions" << endl;
    }
  volume<float> result(source);
  result = 0.0f;
  int midx, midy, midz;
  midz=(kernel.maxz() - kernel.minz())/2;
  midy=(kernel.maxy() - kernel.miny())/2;
  midx=(kernel.maxx() - kernel.minx())/2;
  for (int z=source.minz(); z<=source.maxz(); z++) {
    for (int y=source.miny(); y<=source.maxy(); y++) {
      for (int x=source.minx(); x<=source.maxx(); x++) {
	if (source(x,y,z)>0.0) {
	  for (int mz=kernel.minz(); mz<=kernel.maxz(); mz++) {
	    for (int my=kernel.miny(); my<=kernel.maxy(); my++) {
	      for (int mx=kernel.minx(); mx<=kernel.maxx(); mx++) {
		result(x+mx-midx,y+my-midy,z+mz-midz) += kernel(mx,my,mz);
	      }
	    }
	  }
	}
      }
    }
  }
  source.setextrapolationmethod(oldex);
  return result;
}


volume<float> convolve3(const volume<float>& chi, 
			const volume<float>& kernel,
			complexvolume& kernelc)
{
  if (direct_conv.value()) {
    return direct_convolve(chi,kernel);
  }
  return convolve3_complex(chi,kernel,kernelc);
}


/////////////////////////////////////////////////////////////////////////////

volume4D<float> BzField(const volume4D<float>& chi1, 
			float gx, float gy, float gz,
			float b0x, float b0y, float b0z)
{
  // Implements the first order perturbation solution which is:
  //    Bz(1) = delta/(1+chi0) . ( (1+chi0)/(3+chi0) . chi1.Bz(0) - { ...
  //                        (d2G/dxdz)*(chi1.Bx(0)) + ...
  //                        (d2G/dydz)*(chi1.By(0)) + (d2G/dz2)*(chi1.Bz(0)) } )
  // NOTE: 1/3 in first term comes from using the Lorentz Correction
  //       which is -2/3 . (chi / (1 + chi)) . Bz
  cout << "Calculating Bz field" << endl;
  volume4D<float> bz;
  volume<float> kernel, chi1b0;
  complexvolume fftkernel, dummy;
  string kernelstr="";
  float tol = 0.001*fabs(delta.value())*fabs(b0z);
  for (int t=chi1.mint(); t<=chi1.maxt(); t++) {
    // Start with the (1/(3+chi0) * chi1 * B_0(z)) term (non-convolved)
    bz.addvolume(b0z * chi1[t] / (3.0 + chi0.value()) * (1 + chi0.value()));
    if (fabs(gx)>tol) bz[t] += gx * gradientfield(chi1[t],"x")* chi1[t];
    if (fabs(gy)>tol) bz[t] += gy * gradientfield(chi1[t],"y")* chi1[t];
    if (fabs(gz)>tol) bz[t] += gz * gradientfield(chi1[t],"z")* chi1[t];
    // convolve appropriate kernels and B_0*chi1 terms
    //  x-gradient term
    if (fabs(gx)>tol) {
      if (kernelstr != "z,0,x") {
	kernelstr = "z,0,x";
	kernel = calculate_kernel(chi1[t],kernelstr);
	fftkernel = dummy;
      }
      bz[t] -=  gx * convolve3(chi1[t],kernel,fftkernel);
    }
    //  y-gradient term
    if (fabs(gy)>tol) {
      if (kernelstr != "0,z,y") {
	kernelstr = "0,z,y";
	kernel = calculate_kernel(chi1[t],kernelstr);
	fftkernel = dummy;
      }
      bz[t] -=  gy * convolve3(chi1[t],kernel,fftkernel);
    }
    //  z-gradient term
    if (fabs(gz)>tol) {
      if (kernelstr != "-x,0,z") {
	kernelstr = "-x,0,z";
	kernel = calculate_kernel(chi1[t],kernelstr);
	fftkernel = dummy;
      }
      bz[t] -=  gz * convolve3(chi1[t],kernel,fftkernel);
    }
    //  constant B_x term
    if (Max(fabs(b0x),Max(fabs(gx),fabs(gz)))>tol) {
      chi1b0 = b0x * chi1[t];
      if (fabs(gx)>tol) chi1b0 += gx * gradientfield(chi1[t],"z")* chi1[t];
      if (fabs(gz)>tol) chi1b0 -= gz * gradientfield(chi1[t],"x")* chi1[t];
      if (kernelstr != "1,0,0") {
	kernelstr = "1,0,0";
	kernel = calculate_kernel(chi1[t],kernelstr);
	fftkernel = dummy;
      }
      bz[t] -= convolve3(chi1b0,kernel,fftkernel);
    }
    //  constant B_y term
    if (Max(fabs(b0y),fabs(gy))>tol) {
      chi1b0 = b0y * chi1[t];
      if (fabs(gy)>tol) chi1b0 += gy * gradientfield(chi1[t],"z")* chi1[t];
      if (kernelstr != "0,1,0") {
	kernelstr = "0,1,0";
	kernel = calculate_kernel(chi1[t],kernelstr);
	fftkernel = dummy;
      }
      bz[t] -= convolve3(chi1b0,kernel,fftkernel);
    }
    //  constant B_z term
    if (Max(fabs(b0z),Max(fabs(gx),Max(fabs(gy),fabs(gz))))>tol) {
      chi1b0 = b0z * chi1[t];
      if (fabs(gx)>tol) chi1b0 += gx * gradientfield(chi1[t],"x")* chi1[t];
      if (fabs(gy)>tol) chi1b0 += gy * gradientfield(chi1[t],"y")* chi1[t];
      if (fabs(gz)>tol) chi1b0 += gz * gradientfield(chi1[t],"z")* chi1[t];
      if (kernelstr != "0,0,1") {
	kernelstr = "0,0,1";
	kernel = calculate_kernel(chi1[t],kernelstr);
	fftkernel = dummy;
      }
      bz[t] -= convolve3(chi1b0,kernel,fftkernel);
    }
  }
  // multiply everything by for delta/(1+chi0) (first order general coefficient)
  bz *= (delta.value()/(1.0 + chi0.value()));
  return bz;
}


volume4D<float> BzField(const volume4D<float>& chi1)
{
  return BzField(chi1,gx.value(),gy.value(),gz.value(),
		 b0x.value(),b0y.value(),b0z.value());
}

volume4D<float> ByField(volume4D<float>& chi1)  // chi1 is const (but needs swapping)
{
  volume4D<float> by;
  // flip object axes
  chi1.swapdimensions(3,1,2);
  by = BzField(chi1,gz.value(),gx.value(),gy.value(),
		 b0z.value(),b0x.value(),b0y.value());
  // restore original axes
  chi1.swapdimensions(2,3,1);
  by.swapdimensions(2,3,1);
  return by;
}

volume4D<float> BxField(volume4D<float>& chi1)  // chi1 is const (but needs swapping)
{
  volume4D<float> bx;
  // flip object axes
  chi1.swapdimensions(2,3,1);  
  bx = BzField(chi1,gy.value(),gz.value(),gx.value(),
		 b0y.value(),b0z.value(),b0x.value());
  // restore original axes
  chi1.swapdimensions(3,1,2);
  bx.swapdimensions(3,1,2);
  return bx;
}


/////////////////////////////////////////////////////////////////////////////

int do_calculation()
{
  volume4D<float> bz, chi;
  read_volume4D(chi,inname.value());
  if (calcxyz.value()) { bz.addvolume(BxField(chi)); }
  if (calcxyz.value()) { bz.addvolume(ByField(chi)); }
  bz.addvolume(BzField(chi));
  save_volume4D(bz,outname.value());
  return 0;
}

/////////////////////////////////////////////////////////////////////////////


int main(int argc, char* argv[])
{
  OptionParser options(title, examples);

  try {
    options.add(inname);
    options.add(outname);
    options.add(gx);
    options.add(gy);
    options.add(gz);
    options.add(b0x);
    options.add(b0y);
    options.add(b0z);
    options.add(delta);
    options.add(chi0);
    options.add(calcxyz);
    options.add(extendboundary);
    options.add(direct_conv);
    options.add(verbose);
    options.add(help);
    
    options.parse_command_line(argc, argv);

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
  return do_calculation();
}


