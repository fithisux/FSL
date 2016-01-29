
/*  POSSUM
    Ivana Drobnjak & Mark Jenkinson
    Copyright (C) 2005-2007 University of Oxford  */

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

#include <iostream>
#include <string>
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "utils/options.h"
  
using namespace NEWIMAGE;
using namespace MISCMATHS;
using namespace Utilities;

// The two strings below specify the title and example usage that is
//  printed out as the help or usage message

string title="spharm_rm\nCopyright(c) 2006, University of Oxford (Mark Jenkinson)";
string examples="spharm_rm [options] -i <input_image> -o <output_image>";

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
Option<bool> unmasked(string("--unmasked"), false,
		  string("do not mask final removal"),
		  false, no_argument);
Option<int> numterms(string("-n"), 9,
		  string("number of terms to remove (order is 1,x,y,z,z^2+(x^2+y^2)/2,zx,zy,xy,x^2-y^2)"),
		  false, requires_argument);
Option<string> inname(string("-i,--in"), string(""),
		  string("input filename"),
		  true, requires_argument);
Option<string> outname(string("-o,--out"), string(""),
		  string("output filename"),
		  true, requires_argument);
Option<string> maskname(string("-m"), string(""),
		  string("mask filename"),
		  false, requires_argument);
int nonoptarg;


/////////////////////////////////////////////////////////////////////////////

volume4D<float> generate_spherical_harmonics(const volume<float>& ref, int n=9)
{
  if (n<=0) {
    cerr << "WARNING::generate_spherical_harmonics:: n must be > 0" << endl;
    volume4D<float> rubbish;
    return rubbish;
  }
  volume4D<float> confounds(ref.xsize(),ref.ysize(),ref.zsize(),n);
  float xdim = ref.xdim(), ydim = ref.ydim(),zdim = ref.zdim();
  confounds.setdims(xdim,ydim,zdim,1.0);
  int xc = ref.xsize()/2;
  int yc = ref.ysize()/2;
  int zc = ref.zsize()/2;
  for (int z=0; z<ref.zsize(); z++) {
    for (int y=0; y<ref.ysize(); y++) {
      for (int x=0; x<ref.xsize(); x++) {
	// Set the spherical harmonics: x,y,z,z^2+(x^2+y^2)/2,zx,zy,xy,x^2-y^2
	float x0=(x-xc)*xdim, y0=(y-yc)*ydim, z0=(z-zc)*zdim;
	if (n>0)  confounds(x,y,z,0) = 1.0;
	if (n>1)  confounds(x,y,z,1) = x0;
	if (n>2)  confounds(x,y,z,2) = y0;
	if (n>3)  confounds(x,y,z,3) = z0;
	if (n>4)  confounds(x,y,z,4) = z0*z0 + 0.5*(x0*x0+y0*y0);
	if (n>5)  confounds(x,y,z,5) = z0*x0;
	if (n>6)  confounds(x,y,z,6) = z0*y0;
	if (n>7)  confounds(x,y,z,7) = x0*y0;
	if (n>8)  confounds(x,y,z,8) = x0*x0-y0*y0;
      }
    }
  }
  return confounds;
}


ColumnVector fit_functions(const volume<float>& invol, 
			   const volume4D<float>& confounds)
{			   
  int num = confounds.tsize();
  ColumnVector params(num), xty(num);
  Matrix xtx(num,num);
  volume<float> dotprod;
  // Form the quantities (X' * X) and (X' * Y)
  for (int n=0; n<confounds.tsize(); n++) {
    for (int m=n; m<confounds.tsize(); m++) {
      dotprod = confounds[n] * confounds[m];
      xtx(m+1,n+1) = dotprod.sum();
      xtx(n+1,m+1) = xtx(m+1,n+1);
    }
    dotprod = confounds[n] * invol;
    xty(n+1) = dotprod.sum();
  }
  // Can now get the parameters by solving the GLM
  //   That is:   B = (X' * X)^{-1} * (X' * Y)
  //   Use the pseudo-inverse just in case there is linear dependency
  //    as it should still give a sensible (though not unique) answer
  params = pinv(xtx) * xty;
  return params;
}


volume<float> remove_fitted_functions(const volume<float>& invol,
				      const volume4D<float>& confounds,
				      const ColumnVector& params)
{
  if (confounds.tsize() != params.Nrows() ) {
    cerr << "WARNING:: Do not have the correct number of parameters required" 
	 << endl;
    cerr << "  Removing only a subset of functions" << endl;
  }
  volume<float> outvol = invol;
  int num = Min(confounds.tsize(),params.Nrows());
  for (int n=0; n<num; n++) { 
    outvol -= ((float) params(n+1)) * confounds[n];
  }
  return outvol;
}


/////////////////////////////////////////////////////////////////////////////


int do_work(int argc, char* argv[])
{
  volume<float> invol, outvol, mask;
  volume4D<float> confounds;

  read_volume(invol,inname.value());
  if (maskname.set()) {
    read_volume(mask,maskname.value());
  } else {
    mask = invol*0.0f + 1.0f;
  }
  outvol = invol;


  if (verbose.value()) {cout << "Generating Spherical Harmonic Terms" << endl;}
  confounds = generate_spherical_harmonics(invol,numterms.value());
  for (int n=0; n<confounds.tsize(); n++) { confounds[n] *= mask; }
  if (verbose.value()) {cout << "Fitting Spherical Harmonics" << endl;}
  ColumnVector fits;
  fits = fit_functions(invol*mask,confounds);
  if (verbose.value()) {
    cout << "Spherical Harmonic Amplitudes: ";
    for (int n=1; n<=fits.Nrows(); n++) { cout << fits(n) << "  "; }  
    cout << endl;
  }
  if (unmasked.value()) {
    if (verbose.value()) {cout<<"Regenerating Spherical Harmonic Terms"<<endl;}
    confounds = generate_spherical_harmonics(invol,numterms.value());
  }
  if (verbose.value()) {cout << "Removing Spherical Harmonics" << endl;}
  outvol = remove_fitted_functions(invol,confounds,fits);

  if (!unmasked.value()) { outvol *= mask; }
  save_volume(outvol,outname.value());
  return 0;
}


/////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{

  Tracer tr("main");
  OptionParser options(title, examples);

  try {
    // must include all wanted options here (the order determines how
    //  the help message is printed)
    options.add(inname);
    options.add(outname);
    options.add(maskname);
    options.add(numterms);
    options.add(verbose);
    options.add(help);
    
    nonoptarg = options.parse_command_line(argc, argv);

    // line below stops the program if the help was requested or 
    //  a compulsory option was not set
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

  // Call the local functions

  return do_work(argc,argv);
}

