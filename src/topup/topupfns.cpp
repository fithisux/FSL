// Definitions of global functions used by fnirt
//
// fnirtfns.cpp
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2012 University of Oxford  */
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

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/shared_ptr.hpp>
#include "newmat.h"
#include "newmatio.h"
#include "utils/options.h"
#include "miscmaths/miscmaths.h"
#include "newimage/newimageall.h"
// #include "warpfns/fnirt_file_reader.h"
#include "topupfns.h"

using namespace std;
using namespace BASISFIELD;

namespace TOPUP {

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Creates a topup_clp object after having done consistency checking on the arguments
//
///////////////////////////////////////////////////////////////////////////////////////////////

topup_clp::topup_clp(const Utilities::Option<string>&            imain,
                     const Utilities::Option<string>&            datain,
                     const Utilities::Option<string>&            out,
                     const Utilities::Option<string>&            fout,
                     const Utilities::Option<string>&            iout,
                     const Utilities::Option<string>&            dfout,
                     const Utilities::Option<string>&            rbmout,
                     const Utilities::Option<string>&            jacout,
                     const Utilities::Option<string>&            lout,
                     const Utilities::Option<vector<float> >&    warpres,
                     const Utilities::Option<vector<int> >&      subsamp,
                     const Utilities::Option<vector<float> >&    fwhm,
                     const Utilities::Option<vector<int> >&      miter,
                     const Utilities::Option<vector<float> >&    lambda,
                     const Utilities::Option<int>&               ssqlambda,
                     const Utilities::Option<vector<int> >&      estmov,
                     const Utilities::Option<vector<int> >&      minmet,
                     const Utilities::Option<string>&            regmod,
                     const Utilities::Option<int>&               sporder,
                     const Utilities::Option<string>&            numprec,
                     const Utilities::Option<string>&            interp,
                     const Utilities::Option<int>&               indscale,
                     const Utilities::Option<int>&               regrid,
                     const Utilities::Option<bool>&              verbose,
                     const Utilities::Option<int>&               debug,
                     const Utilities::Option<bool>&              trace)
  : _imain(imain.value()), _out(out.value()), _fout(fout.value()), _iout(iout.value()), _dfout(dfout.value()), 
    _rbmout(rbmout.value()), _jacout(jacout.value()), _lout(lout.value()), _datafile(datain.value()), 
    _warpres(warpres.value()), _subsamp(subsamp.value()), _fwhm(fwhm.value()), _miter(miter.value()), 
    _lambda(lambda.value()), _ssqlambda((ssqlambda.value()==0) ? false : true), _indscale((indscale.value()==0) ? false : true), 
    _regrid((regrid.value()==0) ? false : true), _estmov(estmov.value()), _sporder(static_cast<unsigned int>(sporder.value())), 
    _verbose(verbose.value()), _trace(trace.value()), _debug(debug.value())
{
  // Parse and assert input
  // Check input volume
  NEWIMAGE::volume4D<float> invol;
  NEWIMAGE::read_volume4D_hdr_only(invol,_imain);   // Will throw if there is a problem
  int nvols = invol.tsize();

  // Check file with Phase-encode vectors and read-out times
  if (int(_datafile.N()) != nvols) throw TopupException(string("topup_clp::topup_clp: Mismatch between ")+imain.value()+string(" and ")+datain.value());
  
  // Go through all parameters that are defined
  // Once per level and make sure they agree

  _nlev = std::max(_subsamp.size(),std::max(_warpres.size(),std::max(_fwhm.size(),std::max(_miter.size(),std::max(_lambda.size(),std::max(_estmov.size(),minmet.value().size()))))));
  if (_nlev > 1) {
    if (_subsamp.size() == 1) _subsamp.resize(_nlev,_subsamp[0]);
    else if (_subsamp.size() != _nlev) throw TopupException("topup_clp::topup_clp: Mismatch between --subsamp and other parameters ");
    if (_warpres.size() == 1) _warpres.resize(_nlev,_warpres[0]);
    else if (_warpres.size() != _nlev) throw TopupException("topup_clp::topup_clp: Mismatch between --warpres and other parameters ");
    if (_fwhm.size() == 1) _fwhm.resize(_nlev,_fwhm[0]);
    else if (_fwhm.size() != _nlev) throw TopupException("topup_clp::topup_clp: Mismatch between --fwhm and other parameters ");
    if (_miter.size() == 1) _miter.resize(_nlev,_miter[0]);
    else if (_miter.size() != _nlev) throw TopupException("topup_clp::topup_clp: Mismatch between --miter and other parameters ");
    if (_lambda.size() == 1) _lambda.resize(_nlev,_lambda[0]);
    else if (_lambda.size() != _nlev) throw TopupException("topup_clp::topup_clp: Mismatch between --lambda and other parameters ");
    if (_estmov.size() == 1) _estmov.resize(_nlev,_estmov[0]);
    else if (_estmov.size() != _nlev) throw TopupException("topup_clp::topup_clp: Mismatch between --estmov and other parameters ");

    if (minmet.value().size() != _nlev && minmet.value().size() != 1) throw TopupException("topup_clp::topup_clp: Mismatch between --minmet and other parameters ");
  }

  // And validate all parameters individually
  // --subsamp  
  for (unsigned int i=0; i<_nlev; i++) {
    if (_subsamp[i] < 1 || invol.xsize()%_subsamp[i] || invol.ysize()%_subsamp[i] || invol.zsize()%_subsamp[i]) {
      throw TopupException("topup_clp::topup_clp: Subsampling levels incompatible with image data.");
    }
  }
  // --minmet, should be either 0 (Levenberg-Marquardt) or 1 (Scaled conjugate gradient)
  _optmet.resize(_nlev,MISCMATHS::NL_LM);
  for (unsigned int i=0; i<_nlev; i++) {
    int tmp_met = (minmet.value().size() == 1) ? minmet.value()[0] : minmet.value()[i];
    if (tmp_met == 0) _optmet[i] = MISCMATHS::NL_LM;
    else if (tmp_met == 1) _optmet[i] = MISCMATHS::NL_SCG;
    else throw TopupException("topup_clp::topup_clp: Unrecognized minimization method code.");
  }
  // Warn if they try to estimate movements with scaled conjugate gradient
  for (unsigned int i=0; i<_nlev; i++) {
    if (_estmov[i] && _optmet[i] == MISCMATHS::NL_SCG) {
      cout << "topup_clp::topup_clp: WARNING, it is generally NOT a good idea to include movements when using the Scaled Conjugate Gradient method." << endl;
    }
  }
  // --numprec
  if (numprec.value() == "float") _numprec = MISCMATHS::BFMatrixFloatPrecision;
  else if (numprec.value() == "double") _numprec = MISCMATHS::BFMatrixDoublePrecision;
  else throw TopupException(string("topup_clp::topup_clp: Unrecognized value ")+numprec.value()+string(" as argument to --numprec."));
    
  // --regmod
  if (regmod.value() == "membrane_energy") _regtype = MembraneEnergy;
  else if (regmod.value() == "bending_energy") _regtype = BendingEnergy;
  else throw TopupException(string("topup_clp::topup_clp: Unrecognized value ")+regmod.value()+string(" as argument to --regmod."));
    
  // --sporder
  if (_sporder < 2 || _sporder > 3) throw TopupException("topup_clp::topup_clp: Only quadratic or cubic splines allowed.");

  // --interp
  if (interp.value() == "linear") _interp = TOPUP::LinearInterp;
  else if (interp.value() == "spline") _interp = TOPUP::SplineInterp;
  else throw TopupException(string("topup_clp::topup_clp: Unrecognized value ")+interp.value()+string(" as argument to --interp."));
}

// Returns suitable regridded matrix size given original matrix size and levels of subsampling
std::vector<unsigned int> topup_clp::Regridding(const NEWIMAGE::volume4D<float>& invols) const
{
  std::vector<unsigned int> ims(3,0);

  if (!_regrid) {
    ims[0] = invols.xsize(); ims[1] = invols.ysize(); ims[2] = invols.zsize();
  }
  else {
    unsigned int maxss = 1;
    for (unsigned int i=0; i<NoOfLevels(); i++) maxss = (SubSampling(i+1)>maxss) ? SubSampling(i+1) : maxss;
    ims[0] = invols.xsize() + maxss; ims[1] = invols.ysize() + maxss; ims[2] = invols.zsize() + maxss;
  }
  return(ims);
}

boost::shared_ptr<topup_clp> parse_topup_command_line(unsigned int   narg,
                                                      char           *args[])
{
  // Define the names and types of all input arguments

  Utilities::Option<string> imaname(string("--imain"), string(""),
      string("\tname of 4D file with images"), true, Utilities::requires_argument);

  Utilities::Option<string> dataname(string("--datain"), string(""),
      string("name of text file with PE directions/times"), true, Utilities::requires_argument);
  
  Utilities::Option<string> outname(string("--out"), string(""),
      string("\tbase-name of output files (spline coefficients (Hz) and movement parameters)"), false, Utilities::requires_argument);
  
  Utilities::Option<string> foutname(string("--fout"), string(""),
      string("\tname of image file with field (Hz)"), false, Utilities::requires_argument);
 
  Utilities::Option<string> ioutname(string("--iout"), string(""),
      string("\tname of 4D image file with unwarped images"), false, Utilities::requires_argument);
 
  Utilities::HiddenOption<string> dfoutname(string("--dfout"), string(""),
      string("\tbasename of 4D image files with displacement fields"), false, Utilities::requires_argument);
 
  Utilities::HiddenOption<string> rbmoutname(string("--rbmout"), string(""),
      string("\tbasename of ascii files with RB matrices"), false, Utilities::requires_argument);
 
  Utilities::HiddenOption<string> jacoutname(string("--jacout"), string(""),
      string("\tbasename of image files with jacobian determinants"), false, Utilities::requires_argument);
 
  Utilities::Option<string> logoutname(string("--logout"),string(""),
			    string("Name of log-file"),false,Utilities::requires_argument);
 
  vector<float> warpresdefault(1,10.0);
  Utilities::Option<vector<float> > warpres(string("--warpres"), warpresdefault,
      string("(approximate) resolution (in mm) of warp basis for the different sub-sampling levels, default 10"), false, Utilities::requires_argument);

  vector<int> subsamplingdefault(1,1);
  Utilities::Option<vector<int> > subsampling(string("--subsamp"), subsamplingdefault,
      string("sub-sampling scheme, default 1"), false, Utilities::requires_argument);

  vector<float> smoothdefault(1,8.0);
  Utilities::Option<vector<float> > smoothing(string("--fwhm"),smoothdefault,
      string("\tFWHM (in mm) of gaussian smoothing kernel, default 8"), false, Utilities::requires_argument);
  
  Utilities::Option<string> configname(string("--config"),string(""),
      string("Name of config file specifying command line arguments"),false,Utilities::requires_argument);

  vector<int> maxiterdefault(1,5);
  Utilities::Option<vector<int> > maxiter(string("--miter"), maxiterdefault,
      string("\tMax # of non-linear iterations, default 5"), false, Utilities::requires_argument);

  vector<float> lambdadefault(1,0.0);
  Utilities::Option<vector<float> > lambda(string("--lambda"),lambdadefault,
      string("Weight of regularisation, default depending on --ssqlambda and --regmod switches. See user documetation."), 
      false, Utilities::requires_argument);

  Utilities::Option<int> ssqlambda(string("--ssqlambda"),1,string("If set (=1), lambda is weighted by current ssq, default 1"),
      false,Utilities::requires_argument);

  vector<int> estmovdefault(1,100);
  Utilities::Option<vector<int> > estimatemovements(string("--estmov"),estmovdefault,
      string("Estimate movements if set, default 1 (true)"),false,Utilities::requires_argument);

  vector<int> minmetdefault(1,0); 
  Utilities::Option<vector<int> > minimisationmethod(string("--minmet"),minmetdefault,
      string("Minimisation method 0=Levenberg-Marquardt, 1=Scaled Conjugate Gradient, default 0 (LM)"),false,Utilities::requires_argument);

  Utilities::Option<string> regularisationmodel(string("--regmod"),string("bending_energy"),
      string("Model for regularisation of warp-field [membrane_energy bending_energy], default bending_energy"),
      false,Utilities::requires_argument); 

  Utilities::Option<int> splineorder(string("--splineorder"),3,
      string("Order of spline, 2->Qadratic spline, 3->Cubic spline. Default=3"),false, Utilities::requires_argument);

  Utilities::Option<string> numprec(string("--numprec"),string("double"),
      string("Precision for representing Hessian, double or float. Default double"),false,Utilities::requires_argument);

  Utilities::Option<string> interpolation(string("--interp"),string("spline"),
      string("Image interpolation model, linear or spline. Default spline"),false,Utilities::requires_argument);

  Utilities::Option<int> individualscaling(string("--scale"),0,
      string("\tIf set (=1), the images are individually scaled to a common mean, default 0 (false)"),false,Utilities::requires_argument);

  Utilities::Option<int> regridding(string("--regrid"),1,
      string("\tIf set (=1), the calculations are done in a different grid, default 1 (true)"),false,Utilities::requires_argument);

  Utilities::Option<bool> help(string("-h,--help"), false,
      string("display help info"), false, Utilities::no_argument);

  Utilities::Option<bool> verbose(string("-v,--verbose"), false,
      string("Print diagonostic information while running"), false, Utilities::no_argument);

  Utilities::HiddenOption<int> debug(string("--debug"),0,
      string("Save debug information while running, levels 0 (no info), 1 (some info), 2 (little more info) or 3 (LOTS of info)"), false, Utilities::requires_argument);
      
  Utilities::HiddenOption<bool> traceprint(string("--trace"),0,
      string(""), false, Utilities::no_argument);
      
  // Some explanatory text

  string title = "topup";
  string examples = string("topup --imain=<some 4D image> --datain=<text file> --config=<text file with parameters> --out=my_topup_results\n");

  // Create and load options-parser

  Utilities::OptionParser options(title, examples);

  try {
    // Load parser

    options.add(imaname);
    options.add(dataname);
    options.add(outname);
    options.add(foutname);
    options.add(ioutname);
    options.add(dfoutname);
    options.add(rbmoutname);
    options.add(jacoutname);
    options.add(logoutname);
    options.add(warpres);
    options.add(subsampling);
    options.add(smoothing);
    options.add(configname);
    options.add(maxiter);
    options.add(lambda);
    options.add(ssqlambda);
    options.add(regularisationmodel);
    options.add(estimatemovements);
    options.add(minimisationmethod);      
    options.add(splineorder);
    options.add(numprec);
    options.add(interpolation);
    options.add(individualscaling);
    options.add(regridding);
    options.add(help);
    options.add(verbose);
    options.add(debug);
    options.add(traceprint);
    options.add(help);

    // Parse command line

    unsigned int i = options.parse_command_line(narg,args);
    if (i < narg) {
      for (; i<narg; i++) {
        cerr << "Unknown input: " << args[i] << endl;
      }
      exit(EXIT_FAILURE);
    }

    // Then parse config file if one present

    string   final_configname;
    if (configname.set()) {
      final_configname = TOPUP::existing_conf_file(configname.value());
      if (!final_configname.length()) {
        cerr << "Cannot find config-file " << configname.value() << endl;
        exit(EXIT_FAILURE);
      }
      options.parse_config_file(final_configname);
    }

    if (help.value() || !options.check_compulsory_arguments()) {
      options.usage();
      exit(EXIT_FAILURE);
    }
  }
  catch(Utilities::X_OptionError& e) {
    options.usage();
    cerr << "\n" << e.what() << "\n";
    exit(EXIT_FAILURE);
  }
  catch(std::exception &e) {
    cerr << "\n" << e.what() << "\n";
    exit(EXIT_FAILURE);
  }

  // Do the higher-level parsing, i.e. make sure this particular combination
  // of options make sense. If it does, create a fnirt_clp object that can
  // be interrogated about the different options.

  boost::shared_ptr<topup_clp>   clp;
  try {
    clp = boost::shared_ptr<topup_clp>(new topup_clp(imaname,dataname,outname,foutname,ioutname,dfoutname,rbmoutname,
						     jacoutname,logoutname,warpres,subsampling,smoothing,maxiter,lambda,ssqlambda,
						     estimatemovements,minimisationmethod,regularisationmodel,splineorder,numprec,
						     interpolation,individualscaling,regridding,verbose,debug,traceprint));
  }
  catch(std::exception& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  }

  // Write out a log-file containing all the parsed/expanded parameters
  try {
    std::ofstream logfs(clp->LogFname().c_str());
    logfs << imaname << endl;
    logfs << dataname << endl;
    logfs << outname << endl;
    if (foutname.set()) logfs << foutname << endl;
    if (ioutname.set()) logfs << ioutname << endl;
    logfs << warpres << endl;
    logfs << subsampling << endl;
    logfs << smoothing << endl;
    logfs << maxiter << endl;
    logfs << lambda << endl;
    logfs << ssqlambda << endl;
    logfs << estimatemovements << endl;
    logfs << minimisationmethod << endl;
    logfs << regularisationmodel << endl;
    logfs << splineorder << endl;
    logfs << numprec << endl;
    logfs << interpolation << endl;
    logfs << individualscaling << endl;
    logfs << regridding << endl;

    logfs.close();
  }
  catch (...) {
    cerr << "Error when writing log-file: " << clp->LogFname() << endl;
    throw;
  }

  return(clp);  
}

////////////////////////////////////////////////////////////////////////////
//
// Try to find existing file matching the name cfname. It will look for
// (in turn) i) cfname , ii) cfname + ".cnf" , 
// iii) $FSLDIR + "/etc/flirtsch/" + cfname and
// iv) $FSLDIR + "/etc/flirtsch/" + cfname + ".cnf"
// and returns the first one that exists.
//
////////////////////////////////////////////////////////////////////////////

string existing_conf_file(const string& cfname)
{
  string    ecfname;
  ifstream  ins;

  ecfname = cfname;
  if (TOPUP::check_exist(ecfname)) return(ecfname);
  if (!TOPUP::extension(ecfname).length()) {        // If no extension explicitly given
    ecfname += string(".cnf");
    if (TOPUP::check_exist(ecfname)) return(ecfname);
  }
  if (!TOPUP::path(cfname).length()) {              // If no path explicitly given
    const char *fsldir_ptr = getenv("FSLDIR");
    ecfname = string(fsldir_ptr) + string("/etc/flirtsch/") + cfname;
    if (TOPUP::check_exist(ecfname)) return(ecfname);
    else if (!TOPUP::extension(ecfname).length()) { // If no path _and_ no extension given
      ecfname += string(".cnf");
      if (TOPUP::check_exist(ecfname)) return(ecfname);
    }
  }
  return(string(""));
}

// New version that _might_ help a little

bool check_exist(const string& fname)
{
  // cout << "Attempting to open file named " << fname << endl;

  std::ifstream  ins(fname.c_str(),std::ios::in);
  return((ins) ? true : false);
}

string path(const string& fullname)
{
  string             pathv;
  string::size_type  idx = fullname.find_last_of("/");
  if (idx == string::npos) pathv = "";
  else pathv = fullname.substr(0,idx+1);

  return(pathv);
}

string filename(const string& fullname)
{
  string             fnamev;

  string::size_type  idx = fullname.find_last_of("/");
  if (idx == string::npos) fnamev = fullname;
  else fnamev = fullname.substr(idx+1);

  idx = fnamev.find_first_of(".");
  if (idx != string::npos) fnamev = fnamev.substr(0,idx);

  return(fnamev);
}

string extension(const string& fullname)
{
  string   extv;

  string::size_type  dotidx = fullname.find_last_of(".");
  if (dotidx != string::npos) {    // If there is a dot
    string::size_type  eopidx = fullname.find_last_of("/");
    if (eopidx != string::npos) {  // If there is an explicit path
      if (dotidx > eopidx) extv = fullname.substr(dotidx);
    }
    else extv = fullname.substr(dotidx);
  }
  else extv = "";

  return(extv);
}


} // End namespace TOPUP
