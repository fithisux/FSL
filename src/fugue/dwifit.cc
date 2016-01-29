/*  dwifit.cc

    Mark Jenkinson

    Copyright (C) 2002 FSL  */

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

#include "newimage/newimageall.h"
#include "utils/options.h"

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

using namespace NEWIMAGE;
using namespace Utilities;

string title="dwifit\nCopyright(c) 2002, FSL";
string examples="dwifit -i <dwivolumes> -m <brainmask> -o <diffusionvol> -b <bmatrix> [options]";

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<bool> noweight(string("-n,--noweight"), false,
		  string("turn off SNR weighting"),
		  false, no_argument);
Option<string> involname(string("-i,--in"), string(""),
			   string("filename for input 4D volume"),
			   true, requires_argument);
Option<string> outvolname(string("-o,--out"), string(""),
			   string("filename for output parameter fits"),
			   true, requires_argument);
Option<string> maskvolname(string("-m,--mask"), string(""),
		      string("filename of mask input volume"),
		      true, requires_argument);
Option<string> bmatrix(string("-b,--matrix"), string(""),
		      string("filename of (ordered) values: b gx gy gz (as text)"),
		      false, requires_argument);
Option<string> designmatrix(string("--designmatrix"), string(""),
		      string("read in the design matrix directly - ignoring b-values"),
		      false, requires_argument);

////////////////////////////////////////////////////////////////////////////


int do_dwifit()
{
  volume4D<float> invol;
  volume<float> mask;
  read_rad_volume4D(invol,involname.value());
  read_rad_volume(mask,maskvolname.value());

  // make numbers better suited for calculations (only affects S0 term)
  float scalefactor = invol[0].mean();
  invol /= scalefactor;

  int ntime, nparams;
  ntime=invol.tsize();

  ColumnVector invnoise2est(ntime);
  {
    volume<float> invmask;
    invmask = 1.0f - mask;  // this could be improved - want only air!
    float noise;
    for (int n=1; n<=ntime; n++) {
      // estimate noise from oppositely masked images
      volume<float> tmpvol = invmask * invol[n-1];
      noise = tmpvol.sumsquares();
      noise /= invmask.sum();
      if (verbose.value()) {cout << "noise =  " << noise << endl;}
      invnoise2est(n) = 1.0/Sqr(noise);
    }
  }
  if ((invnoise2est.MinimumAbsoluteValue())<1e-8) {
    cerr << "WARNING:: near zero noise estimate - using 1 instead" << endl;
    invnoise2est = 1.0;
  }
  if (verbose.value()) { 
    cout << "Inverse noise is " << invnoise2est.t() << endl;
  }


  Matrix outmat;
  { // scope for (most) matrices

    // convert to simple data matrix form (ntime x nvox)
    Matrix Y, logY;
    Y = invol.matrix(mask);
    int nvox= Y.Ncols();
    if (ntime != Y.Nrows()) {
      cerr << "Error in converting image series to matrix form!" << endl;
      exit(-1);
    }
    
    // convert all values to log() and demean
    if (verbose.value()) { cout << "Taking logs of data" << endl; }
    logY=log(Y);
    // estimate noise for each input image
    
    
    // read in b value matrix and convert to design matrix, X
    Matrix X;
    if (!designmatrix.set()) {
      // read in simple set of values [ b gx gy gz ] and construct the
      //   design matrix from them
      Matrix bvals;
      if (!bmatrix.set()) {
	cerr << "Must specify either b-value matrix or design matrix!" << endl;
	exit(-1);
      }
      bvals = read_ascii_matrix(bmatrix.value());
      if (bvals.Nrows() != ntime) {
	cerr << "Number of rows in b-value matrix must be the same as number of images in DWI series" << endl;
	exit(-1);
      }
      X.ReSize(bvals.Nrows(),7);
      // set up rows as: B_11 B_22 B_33 2*B_12 2*B_13 2*B_23 1 
      //    where B_ij = b.g_i.g_j     st. log(S) = log(S0) - sum_ij B_ij D_ij
      // NB: last column models the constant log(S0) component
      for (int bv=1; bv<=bvals.Nrows(); bv++) {
	float bb, gx, gy, gz;
	bb = -bvals(bv,1);
	gx=bvals(bv,2);  gy=bvals(bv,3);  gz=bvals(bv,4);
	float gnorm = sqrt(gx*gx + gy*gy + gz*gz);
	gx /= gnorm;  gy /= gnorm;  gz /= gnorm;
	X(bv,1) = bb * gx * gx;
	X(bv,2) = bb * gy * gy;
	X(bv,3) = bb * gz * gz;
	X(bv,4) = 2.0 * bb * gx * gy;
	X(bv,5) = 2.0 * bb * gx * gz;
	X(bv,6) = 2.0 * bb * gy * gz;
	X(bv,7) = 1.0;
      }
    } else {
      // read in design matrix directly (good for simpler fits)
      X = read_ascii_matrix(designmatrix.value());
      if (X.Nrows() != ntime) {
	cerr << "Number of rows in design matrix must be the same as number of images in DWI series" << endl;
	exit(-1);
      }
    }
    
    if (verbose.value()) { cout << "Design matrix = " << endl << X << endl; }
    nparams=X.Ncols();
    
    // set up the matrices
    Matrix W(ntime,ntime), XtW(ntime,ntime), beta(nparams,1);
    outmat.ReSize(nparams,nvox);

    // fit the model separately for each voxel (as the weighting changes...)
    for (int n=1; n<=nvox; n++) {
      // set up weighting matrix
      W = 0.0;
      for (int m=1; m<=ntime; m++) { 
	if (noweight.value()) {
	  W(m,m) = 1.0;
	} else {
	  W(m,m) = Sqr(Y(m,n)) * invnoise2est(m); 
	}
      }
      // calculate fit (pseudo-inverse method)
      XtW = X.t() * W;
      Matrix tmp = XtW*X;
      beta = pinv(XtW * X) * XtW * logY.SubMatrix(1,ntime,n,n);
      if (!designmatrix.set()) {
	beta(7,1) = exp(beta(7,1))*scalefactor;  // reconstitute S0
      }
      outmat.SubMatrix(1,nparams,n,n)=beta;
      if (verbose.value()) { cout << "."; }
    }
    
  } // end scope for matrices
  
  if (verbose.value()) { cout << endl << "Finished Calculations" << endl; }

  volume4D<float> outvol(invol.xsize(),invol.ysize(),invol.zsize(),nparams);
  outvol.setdims(invol.xdim(),invol.ydim(),invol.zdim(),invol.tdim());
  outvol.setmatrix(outmat,mask,0.0f);
  save_volume4D(outvol,outvolname.value());

  return 0;
}


int main(int argc,char *argv[])
{

  OptionParser options(title, examples);

  try {
    options.add(involname);
    options.add(outvolname);
    options.add(maskvolname);
    options.add(bmatrix);
    options.add(designmatrix);
    options.add(noweight);
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

  return do_dwifit();
}






