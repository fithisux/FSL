/*  smoothest.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2000-2004 University of Oxford  */

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

#include <cmath>
#include <iostream>
#include <string>
#include <map>

#include "smoothest.h"

#include "utils/options.h"

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

using namespace Utilities;
using namespace NEWIMAGE;

Option<bool> verbose(string("-V,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<float> dof(string("-d,--dof"), 100.0,
		  string("number of degrees of freedom"),
		  true, requires_argument);
Option<string> maskname(string("-m,--mask"), "mask",
			string("brain mask volume"),
			true, requires_argument);
Option<string> zstatname(string("-z,--zstat"), "zstat",
			 string("filename of zstat image (not with -d)"),
			 true, requires_argument);
Option<string> residname(string("-r,--res"), "res4d",
			 string("filename of `residual-fit' image (use -d)"),
			 true, requires_argument);

namespace SMOOTHEST {

class Interpolate {
public:
  Interpolate() {
    lut[5]   = 1.5423138; lut[6]   = 1.3757105; lut[7]   = 1.2842680;
    lut[8]   = 1.2272151; lut[9]   = 1.1885232; lut[10]  = 1.1606988;
    lut[11]  = 1.1398000; lut[12]  = 1.1235677; lut[13]  = 1.1106196;
    lut[14]  = 1.1000651; lut[15]  = 1.0913060; lut[16]  = 1.0839261;
    lut[17]  = 1.0776276; lut[18]  = 1.0721920; lut[19]  = 1.0674553;
    lut[20]  = 1.0632924; lut[25]  = 1.0483053; lut[30]  = 1.0390117;
    lut[40]  = 1.0281339; lut[50]  = 1.0219834; lut[60]  = 1.0180339;
    lut[70]  = 1.0152850; lut[80]  = 1.0132621; lut[90]  = 1.0117115;
    lut[100] = 1.0104851; lut[150] = 1.0068808; lut[200] = 1.0051200;
    lut[300] = 1.0033865; lut[500] = 1.0020191;
  }

  inline float operator()(float v) {

    float retval = 0;

    if (v<6) return 1.1; // ?? no idea - steve ??

    map<int, float>::iterator i = lut.lower_bound(int(v));

    if(i != lut.end()) {
      if(i != lut.begin()) {
	map<int, float>::iterator j = i--;

	retval = (j->second - i->second)/(j->first - i->first)*(v - i->first)
	  + j->second;
      } 
    } else {
      retval = 1.0321/v + 1;
    }
    
    retval = pow(retval, 0.5);

    return retval;
  }

private:
  map<int, float> lut;
};

Interpolate interpolate;

}

//////////////////////////////////////////////////////////////////////////////
// Standardise the residual field (assuming gaussianity)
unsigned long standardise(volume<float>& mask, 
			  volume4D<float>& R)
{
  unsigned long count = 0;
  int M=R.tsize();

  for (int z=mask.minz(); z<=mask.maxz(); z++) {
    for (int y=mask.miny(); y<=mask.maxy(); y++) {
      for (int x=mask.minx(); x<=mask.maxx(); x++) {

	if( mask(x,y,z) > 0.5) {
	  
	  count ++;
	  
	  if( M > 2 ) {
	    
	    // For each voxel 
	    //    calculate mean and standard deviation...
	    double Sx = 0.0, SSx = 0.0;
	    
	    for ( int t = 0; t < M; t++ ) {
	      double R_it = R(x,y,z,t);
	      
	      Sx += R_it;
	      SSx += Sqr(R_it);
	    }
	    
	    double mean = Sx / M;
	    double sdsq = (SSx - (Sqr(Sx) / M)) / (M - 1) ;
	    
	    if (sdsq<=0) {
	      // trap for differences between mask and invalid data
	      mask(x,y,z)=0;
	      count--;
	    } else {
	      //    ... and use them to standardise to N(0, 1).
	      for ( unsigned short t = 0; t < M; t++ ) {
		R(x,y,z,t) = (R(x,y,z,t) - mean) / sqrt(sdsq);
	      }
	    } 
	  }
	}
      }
    }
  }  
  return count;
}


string title = "\
smoothest \nCopyright(c) 2000-2002, University of Oxford (Dave Flitney and Mark Jenkinson)";

string examples = "\
\tsmoothest -d <number> -r <filename> -m <filename>\n\
\tsmoothest -z <filename> -m <filename>";

int main(int argc, char **argv) {

  OptionParser options(title, examples);

  options.add(verbose);
  options.add(help);
  options.add(dof);
  options.add(maskname);
  options.add(residname);
  options.add(zstatname);

  options.parse_command_line(argc, argv);

//  if(verbose.value()) options.check_compulsory_arguments(true);

  if(help.value()) {
    options.usage();
    exit(EXIT_SUCCESS);
  }

  if( !((zstatname.set() && residname.unset()) || (residname.set() && zstatname.unset())) ||
      (zstatname.set() && (residname.set() || dof.set())) || 
      (residname.set() && dof.unset()) ||
      (maskname.unset()) ) {
    options.usage();
    cerr << endl;
    cerr << "***************************************************************************" << endl;
    cerr << "You must specify either a zstat image OR a 4d residual image." << endl;
    cerr << "If processing a zstat image then you should not set degrees of freedom." << endl; 
    cerr << "You must specify a mask volume image filename" << endl; 
    cerr << "You must specify the degrees of freedom for processing a 4d residual image." << endl; 
    cerr << "***************************************************************************" << endl;
    cerr << endl;
    exit(EXIT_FAILURE);    
  }

  if(verbose.value()) {
    cout << "verbose = " << verbose.value() << endl;
    cout << "help = " << help.value() << endl;
    cout << "dof = " << dof.value() << endl;
    cout << "maskname = " << maskname.value() << endl;
    cout << "residname = " << residname.value() << endl;
    cout << "zstatname = " << zstatname.value() << endl;
  }

  // Read the AVW mask image (single volume)
  if(verbose.value()) cerr << "Reading mask....";
  volume<float> mask;
  read_volume(mask,maskname.value());
  if(verbose.value()) cerr << "done" << endl;
  if (verbose.value()) print_volume_info(mask,"mask");

  if(verbose.value()) cerr << "Reading datafile....";
  string datafilename;
  if(residname.set()) {
    // Read the AVW residual images (array of volumes)
    datafilename = residname.value();
  } else {
    // Read the AVW zstat image (array of one volume)
    datafilename = zstatname.value();
  }
  if(verbose.value()) cerr << "done" << endl;
  

  volume4D<float> R;
  read_volume4D(R,datafilename);
  if (verbose.value()) print_volume_info(R,"Data (residuals/zstat)");

  if (!samesize(R[0],mask)) {
    cerr << "Mask and Data (residuals/zstat) volumes MUST be the same size!"
	 << endl;
    exit(EXIT_FAILURE);
  }

  if(verbose.value()) cerr << "Standardising....";
  unsigned long mask_volume = standardise(mask, R);
  if(verbose.value()) cerr << "done" << endl;
  
  if(verbose.value()) cerr << "Masked-in voxels = " << mask_volume << endl;

  unsigned long N = 0;

  // MJ additions to make it cope with 2D images
  bool usez = true;
  if (R.zsize() <= 1) { usez = false; }
  if ((!usez) && verbose.value()) {
    cout << "Using 2D image mode." << endl;
  }

  // Estimate the smoothness of the normalised residual field
  // see TR00DF1 for mathematical description of the algorithm.
  enum {X = 0, Y, Z};
  double SSminus[3] = {0, 0, 0}, S2[3] = {0, 0, 0};

  int zstart=1;
  if (!usez) zstart=0;
  for ( unsigned short z = zstart; z < R.zsize() ; z++ )
    for ( unsigned short y = 1; y < R.ysize() ; y++ )
      for ( unsigned short x = 1; x < R.xsize() ; x++ )
	// Sum over N
	if( (mask(x, y, z)>0.5) &&
	    (mask(x-1, y, z)>0.5) && 
	    (mask(x, y-1, z)>0.5) && 
	    ( (!usez) || (mask(x, y, z-1)>0.5) ) ) {
	  
	  N++;
	  
	  for ( unsigned short t = 0; t < R.tsize(); t++ ) {
	    // Sum over M
	    SSminus[X] += R(x, y, z, t) * R(x-1, y, z, t);
	    SSminus[Y] += R(x, y, z, t) * R(x, y-1, z, t);
	    if (usez) SSminus[Z] += R(x, y, z, t) * R(x, y, z-1, t);

	    S2[X] += 0.5 * (Sqr(R(x, y, z, t)) + Sqr(R(x-1, y, z, t)));
	    S2[Y] += 0.5 * (Sqr(R(x, y, z, t)) + Sqr(R(x, y-1, z, t)));
	    if (usez) S2[Z] += 0.5 * (Sqr(R(x, y, z, t)) + Sqr(R(x, y, z-1, t)));
	  }
	}

  double norm = 1.0/(double) N;
  double v = dof.value();	// v - degrees of freedom (nu)  
  if(R.tsize() > 1) {
    if(verbose.value()) {
      cerr << "Non-edge voxels = " << N << endl;
      cerr << "(v - 2)/(v - 1) = " << (v - 2)/(v - 1) << endl;
    }
    norm = (v - 2) / ((v - 1) * N * R.tsize());
  }
 
//    SSminus[X] *= norm;
//    SSminus[Y] *= norm;
//    SSminus[Z] *= norm;
  
//    S2[X] *= norm;
//    S2[Y] *= norm;
//    S2[Z] *= norm;

  if(verbose.value()) {
    cout << "SSminus[X] = " << SSminus[X] << ", SSminus[Y] = " << SSminus[Y] << ", SSminus[Z] = " << SSminus[Z] 
	 << ", S2[X] = " << S2[X] << ", S2[Y] = " << S2[Y] << ", S2[Z] = " << S2[Z]
	 << endl;
  }

  // for extreme smoothness 
  if (SSminus[X]>=0.99999999*S2[X]) {
    SSminus[X]=0.99999*S2[X];  
    cerr << "WARNING: Extreme smoothness detected in X - possibly biased"
	 << " global estimate." << endl; } 
  if (SSminus[Y]>=0.99999999*S2[Y]) {
    SSminus[Y]=0.99999*S2[Y];
    cerr << "WARNING: Extreme smoothness detected in Y - possibly biased"
	 << " global estimate." << endl; } 
  if (usez) {
    if (SSminus[Z]>=0.99999999*S2[Z]) {
      SSminus[Z]=0.99999*S2[Z];
      cerr << "WARNING: Extreme smoothness detected in Z - possibly biased"
	   << " global estimate." << endl; } 
  }

  // Convert to sigma squared
  double sigmasq[3];

  sigmasq[X] = -1.0 / (4 * log(fabs(SSminus[X]/S2[X])));
  sigmasq[Y] = -1.0 / (4 * log(fabs(SSminus[Y]/S2[Y])));
  if (usez) { sigmasq[Z] = -1.0 / (4 * log(fabs(SSminus[Z]/S2[Z]))); }
  else { sigmasq[Z]=0; }
  // the following is determininant of Lambda to the half 
  //   i.e. dLh = | Lambda |^(1/2)
  // Furthermore, W_i = 1/(2.lambda_i) = sigma_i^2 => 
  //   det(Lambda) = det( lambda_i ) = det ( (2 W_i)^-1 ) = (2^D det(W))^-1
  //   where D = number of dimensions (2 or 3)
  double dLh;
  if (usez) { dLh=pow(sigmasq[X]*sigmasq[Y]*sigmasq[Z], -0.5)*pow(8, -0.5); }
  else { dLh = pow(sigmasq[X]*sigmasq[Y], -0.5)*pow(4, -0.5); }

  if(verbose.value()) {
    cout << "DLH " << dLh << " voxels^-3 before correcting for temporal DOF" << endl;
  }
  if(R.tsize() > 1) dLh *= SMOOTHEST::interpolate(v);

  // Convert to full width half maximum
  double FWHM[3];
  FWHM[X] = sqrt(8 * log(2) * sigmasq[X]);
  FWHM[Y] = sqrt(8 * log(2) * sigmasq[Y]);
  if (usez) { FWHM[Z] = sqrt(8 * log(2) * sigmasq[Z]); }
  else { FWHM[Z]=0; }
  double resels = FWHM[X] * FWHM[Y];
  if (usez) resels *= FWHM[Z];

  if(verbose.value()) {
    cout << "FWHMx = " << FWHM[X] << " voxels, " 
	 << "FWHMy = " << FWHM[Y] << " voxels"; 
    if (usez) cout << ", FWHMz = " << FWHM[Z] << " voxels";
    cout << endl;
  }
  
  FWHM[X] *= R.xdim(); FWHM[Y] *= R.ydim(); if (usez) FWHM[Z] *= R.zdim(); 

  if(verbose.value()) {
    cout << "FWHMx = " << FWHM[X] << " mm, "
	 << "FWHMy = " << FWHM[Y] << " mm";
    if (usez) cout << ", FWHMz = " << FWHM[Z] << " mm";
    cout << endl;
    cout << "DLH " << dLh << " voxels^-3" << endl;
    cout << "VOLUME " << mask_volume << " voxels" << endl;
    cout << "RESELS " << resels << " voxels per resel" << endl;
  }
  
  cout << "DLH " << dLh << endl;
  cout << "VOLUME " << mask_volume << endl;
  cout << "RESELS " << resels << endl;

  return EXIT_SUCCESS;
}


