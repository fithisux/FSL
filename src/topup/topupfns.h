// Declarations of utility functions/classes used by topup
//
// topupfns.h
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2012 University of Oxford 

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

#ifndef topupfns_h
#define topupfns_h

#include <string>
#include <vector>
#include "newmat.h"
#include "newimage/newimage.h"
#include "basisfield/basisfield.h"
#include "basisfield/splinefield.h"
#include "basisfield/dctfield.h"
#include "utils/options.h"
#include "miscmaths/nonlin.h"
#include "topup_costfunctions.h"
#include "topup_file_io.h"

namespace TOPUP {

std::string path(const std::string& fullname);
std::string filename(const std::string& fullname);
std::string extension(const std::string& fullname);

enum PrecisionType         {FloatPrecision, DoublePrecision};

///////////////////////////////////////////////////////////////////////////////////////////////
//
// topup_clp is a glorified struct that holds the command line parameters of topup
//
///////////////////////////////////////////////////////////////////////////////////////////////

class topup_clp
{
public:
  topup_clp(const Utilities::Option<string>&            imain,
	    const Utilities::Option<string>&            datain,
	    const Utilities::Option<string>&            cout,
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
	    const Utilities::Option<bool>&              trace);
  const std::string& ImaFname()  const { return(_imain); }
  std::string CoefFname() const {
    if (!_out.length()) _out = TOPUP::path(_imain) + TOPUP::filename(_imain);
    return(_out + string("_fieldcoef"));
  }
  std::string MovParFname() const {
    if (!_out.length()) _out = TOPUP::path(_imain) + TOPUP::filename(_imain);
    return(_out + string("_movpar.txt"));
  }
  const std::string& LogFname() const {
    if (!_lout.length()) _lout = TOPUP::path(_imain) + TOPUP::filename(_imain) + string(".topup_log");
    return(_lout);
  }
  const std::string& FieldFname() const { return(_fout); }
  const std::string& ImaOutFname() const { return(_iout); }
  const std::string& DisplacementFieldBaseFname() const { return(_dfout); }
  const std::string& RigidBodyBaseFname() const { return(_rbmout); }
  const std::string& JacobianBaseFname() const { return(_jacout); }
  const NEWMAT::Matrix PhaseEncodeVectors() const { return(_datafile.PhaseEncodeVectors()); }
  const NEWMAT::ColumnVector ReadoutTimes() const { return(_datafile.ReadOutTimes()); }
  unsigned int NoOfLevels() const { return(_nlev); }
  TopupInterpolationType InterpolationModel() const { return(_interp); }
  MISCMATHS::BFMatrixPrecisionType HessianPrecision() const { return(_numprec); }
  RegularisationType RegularisationModel() const { return(_regtype); }
  unsigned int SplineOrder() const { return(_sporder); }
  unsigned int DebugLevel() const { return(_debug); }
  bool Verbose() const { return(_verbose); }
  bool Trace() const { return(_trace); }
  bool SSQLambda() const { return(_ssqlambda); }
  bool IndividualScaling() const { return(_indscale); }
  std::vector<unsigned int> Regridding(const NEWIMAGE::volume4D<float>& invols) const; 
  unsigned int SubSampling(unsigned int level) const {
    if (level < 1 || level > _nlev) throw TopupException("topup_clp::SubSampling: Out-of-range value of level");
    return(static_cast<unsigned int>(_subsamp[level-1]));
  }
  double WarpRes(unsigned int level) const {
    if (level < 1 || level > _nlev) throw TopupException("topup_clp::WarpRes: Out-of-range value of level");
    return(_warpres[level-1]);
  }
  double FWHM(unsigned int level) const {
    if (level < 1 || level > _nlev) throw TopupException("topup_clp::FWHM: Out-of-range value of level");
    return(_fwhm[level-1]);
  }
  double Lambda(unsigned int level) const {
    if (level < 1 || level > _nlev) throw TopupException("topup_clp::Lambda: Out-of-range value of level");
    return(_lambda[level-1]);
  }
  unsigned int MaxIter(unsigned int level) const {
    if (level < 1 || level > _nlev) throw TopupException("topup_clp::MaxIter: Out-of-range value of level");
    return(static_cast<unsigned int>(_miter[level-1]));
  }
  bool EstimateMovements(unsigned int level) const {
    if (level < 1 || level > _nlev) throw TopupException("topup_clp::EstimateMovements: Out-of-range value of level");
    return(static_cast<bool>(_estmov[level-1] != 0));
  }
  MISCMATHS::NLMethod OptimisationMethod(unsigned int level) const { 
    if (level < 1 || level > _nlev) throw TopupException("topup_clp::OptimisationMethod: Out-of-range value of level");
    return(_optmet[level-1]); 
  }

private:
  unsigned int                           _nlev;
  std::string                            _imain;
  mutable std::string                    _out;
  std::string                            _fout;
  std::string                            _iout;
  std::string                            _dfout;
  std::string                            _rbmout;
  std::string                            _jacout;
  mutable std::string                    _lout;
  TopupDatafileReader                    _datafile;
  std::vector<float>                     _warpres;
  std::vector<int>                       _subsamp;
  std::vector<float>                     _fwhm;
  std::vector<int>                       _miter;
  std::vector<float>                     _lambda;
  bool                                   _ssqlambda;
  bool                                   _indscale;
  bool                                   _regrid;
  RegularisationType                     _regtype;
  TopupInterpolationType                 _interp;
  MISCMATHS::BFMatrixPrecisionType       _numprec;
  std::vector<MISCMATHS::NLMethod>       _optmet;
  std::vector<int>                       _estmov;
  unsigned int                           _sporder;
  bool                                   _verbose;
  bool                                   _trace;
  unsigned int                           _debug;
};

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Declarations of global functions used by topup
//
///////////////////////////////////////////////////////////////////////////////////////////////

boost::shared_ptr<topup_clp> parse_topup_command_line(unsigned int   narg,
                                                      char           *args[]);
bool check_exist(const std::string& fname);
std::string existing_conf_file(const std::string& cfname);

} // End namespace TOPUP

#endif // End #ifndef topupfns_h
