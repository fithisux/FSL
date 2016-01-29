/*  retroicor_evs.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2008 University of Oxford  */

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

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "utils/options.h"

using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;

// The two strings below specify the title and example usage that is
//  printed out as the help or usage message

string title="retroicor_evs (Version 1.0)\nCopyright(c) 2008, University of Oxford (Mark Jenkinson)";
string examples="retroicor_evs [options] -c <cardiac phase file> -r <respiratory phase file> -o <confound matrix> --oc=3 --or=2 --oac=2 --oar=1 --tcsamp=0.1 --trsamp=0.5 --nt=100 --tr=3.0";

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
Option<int> ntimepoints(string("--nt"), 0,
		  string("number of timepoints for EVs (equal to number of volumes in fMRI)"),
		  true, requires_argument);
Option<int> cardorder(string("--oc"), 2,
		  string("order of cardiac regressors (number of Fourier pairs) - default=2"),
		  true, requires_argument);
Option<int> resporder(string("--or"), 1,
		  string("order of respiratory regressors (number of Fourier pairs) - default=1"),
		  true, requires_argument);
Option<int> amodcorder(string("--oac"), 1,
		  string("order of amplitude modulated cardiac regressors - default=1"),
		  true, requires_argument);
Option<int> amodrorder(string("--oar"), 1,
		  string("order of amplitude modulated respiratory regressors - default=1"),
		  true, requires_argument);
Option<float> tcsamp(string("--tcsamp"), 0.0,
		  string("sampling time for cardiac phase input (in seconds)"),
		  true, requires_argument);
Option<float> trsamp(string("--trsamp"), 0.0,
		  string("sampling time for respiratory phase input (in seconds)"),
		  true, requires_argument);
Option<float> tr(string("--tr"), 0.0,
		  string("TR of fMRI volumes (in seconds)"),
		  true, requires_argument);
Option<float> toffset(string("--toffset"), 0.0,
		  string("offset time for EVs; time of first acquisition relative to phase inputs -  (in seconds)"),
		  false, requires_argument);
Option<string> cardname(string("-c,--cardiac"), string(""),
		  string("input filename for cardiac phase values"),
		  true, requires_argument);
Option<string> respname(string("-r,--respiratory"), string(""),
		  string("input filename for respiratory phase values"),
		  true, requires_argument);
Option<string> outname(string("-o,--out"), string(""),
		  string("output filename (for confound/EV matrix)"),
		  true, requires_argument);
int nonoptarg;

////////////////////////////////////////////////////////////////////////////

// Local functions

int sanity_check()
{
  // sanity checking
  if (cardorder.value()<0) { 
    cerr << "Invalid order for cardiac (" << cardorder.value() << ") - must be non-negative" << endl; 
    return 1;
  }
  if (resporder.value()<0) { 
    cerr << "Invalid order for respiratory (" << resporder.value() << ") - must be non-negative" << endl; 
    return 2;
  }
  if (amodcorder.value()<0) { 
    cerr <<"Invalid order for amplitude modulated cardiac terms (" << amodcorder.value() << ") - must be non-negative"<<endl; 
    return 3;
  }
  if (amodrorder.value()<0) { 
    cerr <<"Invalid order for amplitude modulated respiratory terms (" << amodrorder.value() << ") - must be non-negative"<<endl; 
    return 4;
  }
  int nt=ntimepoints.value();
  if (nt<2) { 
    cerr << "Invalid number of timepoints (" << nt << ")" << endl;
    return 5;
  }
  return 0;
}


Matrix calc_confoundmat(int cardorder, int resporder, int amodcorder, int amodrorder, 
			int nt, float tcsmap, float trsamp, float tr, float toffset)
{
  Matrix confoundmat(nt,cardorder*2 + resporder*2 + amodcorder*amodrorder*4);

  // create timing for cardiac phase, respiratory phase and slice timing of EPI
  ColumnVector tcard, tresp, tslice(nt);
  tcard = cardph;
  for (int n=1; n<=cardph.Nrows(); n++) { tcard(n) = (n-1)*tcsamp; }
  for (int n=1; n<=respph.Nrows(); n++) { tresp(n) = (n-1)*trsamp; }
  for (int n=1; n<=nt; n++) { tslice(n) = tr*n + toffset; }
  // resample phase values at the appropriate slice timings
  ColumnVector cph_slice(nt), rph_slice(nt);
  for (int n=1; n<=nt; n++) {
    cph_slice(n) = interp1(tcard,cardph,tslice(n));
    rph_slice(n) = interp1(tresp,respph,tslice(n));
  }

  // generate the required regressors
  int col=1;
  // cardiac-only terms
  for (int m=1; m<=cardorder; m++) {
    for (int n=1; n<=nt; n++) {
      confoundmat(col,n) = cos(cph_slice(n)*m);
      confoundmat(col+1,n) = sin(cph_slice(n)*m);
    }
    col+=2;
  }
  // respiratory-only terms
  for (int m=1; m<=resporder; m++) {
    for (int n=1; n<=nt; n++) {
      confoundmat(col,n) = cos(rph_slice(n)*m);
      confoundmat(col+1,n) = sin(rph_slice(n)*m);
    }
    col+=2;
  }
  // amplitude modulated terms
  for (int mc=1; mc<=amodcorder; mc++) {
    for (int mr=1; mr<=amodrorder; mr++) {
      for (int n=1; n<=nt; n++) {
	confoundmat(col,n)   = cos(cph_slice(n)*mc + rph_slice(n)*mr);
	confoundmat(col+1,n) = sin(cph_slice(n)*mc + rph_slice(n)*mr);
	confoundmat(col+2,n) = cos(cph_slice(n)*mc - rph_slice(n)*mr);
	confoundmat(col+3,n) = sin(cph_slice(n)*mc - rph_slice(n)*mr);
      }
      col+=4;
    }
  }
  return confoundmat;
}


int do_work(int argc, char* argv[]) 
{
  if (sanity_check()!=0) { exit(EXIT_FAILURE); }

  // read in phase values
  ColumnVector cardph, respph;
  cardph = read_ascii_matrix(cardname.value());
  respph = read_ascii_matrix(respname.value());

  // setup matrix
  Matrix confoundmat;

  confoundmat = calc_confoundmat(cardorder.value(), resporder.value(), 
				 amodcorder.value(), amodrorder.value(), 
				 ntimepoints.value(), tcsmap.value(), trsamp.value(), 
				 tr.value(), toffset.value());

  // save output
  write_ascii_matrix(confoundmat,outname.value());
  return 0;

}

////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{

  Tracer tr_main("main");
  OptionParser options(title, examples);

  try {
    // must include all wanted options here (the order determines how
    //  the help message is printed)
    options.add(cardname);
    options.add(respname);
    options.add(outname);
    options.add(cardorder);
    options.add(resporder);
    options.add(amodrorder);
    options.add(tcsamp);
    options.add(trsamp);
    options.add(ntimepoints);
    options.add(tr);
    options.add(toffset);
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

