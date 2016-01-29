// Declarations of global utility functions used by fnirt
//
// fnirtfns.h
//
// Jesper Andersson, FMRIB Image Analysis Group
//
//
/*    Copyright (C) 2012 University of Oxford  */

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


#ifndef fnirtfns_h
#define fnirtfns_h

#include <string>
#include <vector>
#include "newmat.h"
#include "newimage/newimage.h"
#include "basisfield/basisfield.h"
#include "basisfield/splinefield.h"
#include "basisfield/dctfield.h"
#include "utils/options.h"
#include "miscmaths/nonlin.h"
#include "fnirt_costfunctions.h"
#include "intensity_mappers.h"

namespace FNIRT {

std::string path(const std::string& fullname);
std::string filename(const std::string& fullname);
std::string extension(const std::string& fullname);

class fnirt_error: public std::exception
{
public:
  fnirt_error(const string& pmsg) throw() : msg(pmsg) {}
  const char *what() const throw() {return(string("fnirt::" + msg).c_str());}
  ~fnirt_error() throw() {}
private:
  string msg;
};

///////////////////////////////////////////////////////////////////////////////////////////////
//
// fnirt_clp is a glorified struct that holds the command line parameters of fnirt
//
///////////////////////////////////////////////////////////////////////////////////////////////

enum MaskType {InclusiveMask, ExclusiveMask, IgnoreMask}; 
enum CostFunctionType {SSD};
enum BasisFieldType   {Spline, DCT};
enum IntensityMappingType {NONE, GLOBAL_LINEAR, GLOBAL_NON_LINEAR, LOCAL_LINEAR, 
			   LOCAL_BIAS_WITH_GLOBAL_NON_LINEAR, LOCAL_NON_LINEAR};

class fnirt_clp
{
private:
  unsigned int                                 nlev;
  std::string                                  ref;
  std::string                                  obj;
  std::string                                  inwarp;
  std::string                                  in_int;
  mutable std::string                          coef;
  std::string                                  objo;
  std::string                                  fieldo;
  std::string                                  jaco;
  std::string                                  refo;
  std::string                                  into;
  std::string                                  logo;
  std::string                                  refm;
  std::string                                  objm;
  std::string                                  ref_pl;
  std::string                                  obj_pl;
  bool                                         rimf;
  double                                       rimv;
  bool                                         oimf;
  double                                       oimv;
  CostFunctionType                             cf;
  BasisFieldType                               bf;
  NEWMAT::Matrix                               aff;
  std::vector<unsigned int>                    ss;
  std::vector<unsigned int>                    mi;
  unsigned int                                 spordr;
  std::vector<float>                           ofwhm;
  std::vector<float>                           rfwhm;
  std::vector<std::vector<unsigned int> >      ksp;
  std::vector<std::vector<unsigned int> >      dco;
  RegularisationType                           regmod;
  std::vector<float>                           lambda;
  bool                                         ssqlambda;
  std::vector<float>                           jacrange;
  bool                                         userefderiv;
  IntensityMappingType                         imt;
  std::vector<bool>                            estint;
  std::vector<bool>                            userefmask;
  std::vector<bool>                            useobjmask;
  unsigned int                                 intord;
  RegularisationType                           biasregmod;
  std::vector<float>                           biaslambda;
  std::vector<float>                           mpl_lambda;
  std::vector<std::vector<unsigned int> >      bias_ksp;
  std::vector<std::vector<unsigned int> >      bias_dco;
  MISCMATHS::NLMethod                          nlm;
  bool                                         verbose;
  unsigned int                                 debug;
  BFMatrixPrecisionType                        hess_prec;
  FnirtInterpolationType                       interp_type;

public:
  fnirt_clp(const Utilities::Option<std::string>&                     pref,
            const Utilities::Option<std::string>&                     pobj,
            const Utilities::Option<std::string>&                     paff,
            const Utilities::Option<std::string>&                     pinwarp,
            const Utilities::Option<std::string>&                     pin_int,
            const Utilities::Option<std::string>&                     pcoef,
            const Utilities::Option<std::string>&                     pobjo,
            const Utilities::Option<std::string>&                     pfieldo,
            const Utilities::Option<std::string>&                     pjaco,
            const Utilities::Option<std::string>&                     prefo,
            const Utilities::Option<std::string>&                     pinto,
            const Utilities::Option<std::string>&                     plogo,
            const Utilities::Option<std::string>&                     prefm,
            const Utilities::Option<std::string>&                     pobjm,
            const Utilities::Option<std::string>&                     pref_pl,
            const Utilities::Option<std::string>&                     pobj_pl,
            const Utilities::Option<std::vector<int> >&               puserefm,
            const Utilities::Option<std::vector<int> >&               puseobjm,
	    const Utilities::Option<int>&                             primf,
	    const Utilities::Option<int>&                             poimf,
	    const Utilities::Option<float>&                           primv,
	    const Utilities::Option<float>&                           poimv,
            const Utilities::Option<std::string>&                     pcf,
            const Utilities::Option<std::string>&                     pbf,
            const Utilities::Option<std::string>&                     pnlm,
            const Utilities::Option<std::vector<int> >&               pmi,
            const Utilities::Option<std::vector<int> >&               pss,
            const Utilities::Option<std::vector<float> >&             pwres,
            const Utilities::Option<int>&                             pspordr,
	    const Utilities::Option<std::vector<float> >&             pofwhm,
            Utilities::Option<std::vector<float> >&                   prfwhm,
	    const Utilities::Option<std::string>&                     pregmod,
            Utilities::Option<std::vector<float> >&                   plambda,
            const Utilities::Option<int>&                             pssqlambda,
            const Utilities::Option<float>&                           pmpl_lambda,
	    const Utilities::Option<std::vector<float> >&             pjacrange,
            const Utilities::Option<int>&                             puserefderiv,
            const Utilities::Option<std::string>&                     pintmod,
	    const Utilities::Option<std::vector<int> >&               pestint,
            const Utilities::Option<int>&                             pintord,
            const Utilities::Option<std::vector<float> >&             pbiasres,
            const Utilities::Option<std::string>&                     pbiasregmod,
            const Utilities::Option<float>&                           pbiaslambda,
            const Utilities::Option<bool>&                            pverbose,
            const Utilities::Option<int>&                             pdebug,
            const Utilities::Option<std::string>&                     p_hess_prec,
            const Utilities::Option<std::string>&                     p_interp_type);
  ~fnirt_clp() {}
  const std::string& Obj() const {return(obj);}
  const std::string& Ref() const {return(ref);}
  const std::string& InWarp() const {return(inwarp);}
  const std::string& InInt() const {return(in_int);}
  const std::string& CoefFname() const {
    // If no name has been supplied we will derive it from the obj image name.
    if (!coef.length()) coef = FNIRT::path(obj) + FNIRT::filename(obj) + string("_warpcoef");
    return(coef);
  }
  const std::string& LogFname() {
    // If not specified, the name of the log-file will be derived from the --ref and --in names
    if (!logo.length()) logo = path(obj) + filename(obj) + string("_to_") + filename(ref) + string(".log");
    return(logo);
  }
  const std::string& ObjOutFname() const {return(objo);}
  const std::string& RefOutFname() const {return(refo);}
  const string& FieldFname() const {return(fieldo);}
  const string& JacFname() const {return(jaco);}
  const string& IntensityMappingFname() const {return(into);}
  const string& ObjMask() const {return(objm);}
  const string& RefMask() const {return(refm);}
  const string& RefPointList() const {return(ref_pl);}
  const string& ObjPointList() const {return(obj_pl);}
  bool UseImplicitRefMask() const {return(rimf);}
  double ImplicitRefValue() const {return(rimv);}
  bool UseImplicitObjMask() const {return(oimf);}
  double ImplicitObjValue() const {return(oimv);}
  int NoLevels() const {return(nlev);}
  CostFunctionType CostFunction() const {return(cf);}
  BasisFieldType Basis() const {return(bf);}
  BFMatrixPrecisionType HessianPrecision() const {return(hess_prec);}
  FnirtInterpolationType InterpolationModel() const {return(interp_type);}
  unsigned int SplineOrder() const {return(spordr);}
  MISCMATHS::NLMethod MinimisationMethod() const {return(nlm);}
  const NEWMAT::Matrix& Affine() const {return(aff);}
  RegularisationType RegularisationModel() const {return(regmod);}
  double JacLowerBound() const {return(jacrange[0]);}
  double JacUpperBound() const {return(jacrange[1]);}
  bool WeightLambdaBySSD() const {return(ssqlambda);}
  bool UseRefDeriv() const {return(userefderiv);}
  std::vector<unsigned int> RefSize() const {
    std::vector<unsigned int>   rsz(3,0);
    NEWIMAGE::volume<float>     vref;
    NEWIMAGE::read_volume_hdr_only(vref,ref);
    rsz[0] = vref.xsize(); rsz[1] = vref.ysize(); rsz[2] = vref.zsize();
    return(rsz);
  } 
  std::vector<double> RefVxs() const {
    std::vector<double>         vxs(3,0);
    NEWIMAGE::volume<float>     vref;
    NEWIMAGE::read_volume_hdr_only(vref,ref);
    vxs[0] = vref.xdim(); vxs[1] = vref.ydim(); vxs[2] = vref.zdim();
    return(vxs);
  }
  std::vector<unsigned int> FullResKsp() const {
    if (bf != Spline) throw fnirt_error("fnirt_clp::FullResKsp: Knot-spacing not defined for DCT basis");
    std::vector<unsigned int>  full_ksp = KnotSpacing(NoLevels());
    for (unsigned int i=0; i<full_ksp.size(); i++) full_ksp[i] /= SubSampling(NoLevels());
    return(full_ksp);
  }
  std::vector<unsigned int> FullResIntKsp() const {
    if (bf != Spline) throw fnirt_error("fnirt_clp::FullResIntKsp: Knot-spacing not defined for DCT basis");
    std::vector<unsigned int>  full_ksp = IntKnotSpacing(NoLevels());
    // for (unsigned int i=0; i<full_ksp.size(); i++) full_ksp[i] /= SubSampling(NoLevels());
    for (unsigned int i=0; i<full_ksp.size(); i++) full_ksp[i] *= SubSampling(NoLevels());
    return(full_ksp);
  }
  unsigned int SubSampling(unsigned int lev) const {
    if (lev < 1 || lev > nlev) throw fnirt_error("fnirt_clp::Subsampling: Out-of-range value of lev");
    return(ss[lev-1]);
  }
  double RefFWHM(unsigned int lev) const {
    if (lev < 1 || lev > nlev) throw fnirt_error("fnirt_clp::RefFWHM: Out-of-range value of lev");
    return(static_cast<double>(rfwhm[lev-1]));
  }
  double ObjFWHM(unsigned int lev) const {
    if (lev < 1 || lev > nlev) throw fnirt_error("fnirt_clp::ObjFWHM: Out-of-range value of lev");
    return(static_cast<double>(ofwhm[lev-1]));
  }
  double Lambda(unsigned int lev) const {
    if (lev < 1 || lev > nlev) throw fnirt_error("fnirt_clp::Lambda: Out-of-range value of lev");
    return(static_cast<double>(lambda[lev-1]));
  }
  unsigned int MaxIter(unsigned int lev) const {
    if (lev < 1 || lev > nlev) throw fnirt_error("fnirt_clp::MaxIter: Out-of-range value of lev");
    return(mi[lev-1]);
  }
  const vector<unsigned int>& KnotSpacing(unsigned int lev) const {
    if (bf != Spline) throw fnirt_error("fnirt_clp::KnotSpacing: Knot spacing not defined for DCT basis");
    if (lev < 1 || lev > nlev) throw fnirt_error("fnirt_clp::KnotSpacing: Out-of-range value of lev");
    return(ksp[lev-1]);
  }
  const vector<unsigned int>& DCTOrder(unsigned int lev) const {
    if (bf != DCT) throw fnirt_error("fnirt_clp::DCTOrder: Order not defined for spline basis");
    if (lev < 1 || lev > nlev) throw fnirt_error("fnirt_clp::DCTOrder: Out-of-range value of lev");
    return(dco[lev-1]);
  }
  IntensityMappingType IntMapType() const {return(imt);}
  bool EstimateIntensity(unsigned int lev) const {
    if (lev < 1 || lev > nlev) throw fnirt_error("fnirt_clp::EstimateIntensity: Out-of-range value of lev");
    return(estint[lev-1]);
  }
  bool UseRefMask(unsigned int lev) const {
    if (lev < 1 || lev > nlev) throw fnirt_error("fnirt_clp::UseRefMask: Out-of-range value of lev");
    return(userefmask[lev-1]);
  }
  bool UseObjMask(unsigned int lev) const {
    if (lev < 1 || lev > nlev) throw fnirt_error("fnirt_clp::UseObjMask: Out-of-range value of lev");
    return(useobjmask[lev-1]);
  }
  unsigned int IntMapOrder() const {return(intord);}
  RegularisationType IntRegType() const {return(biasregmod);}
  double IntLambda(unsigned int lev=1) const {
    if (lev < 1 || lev > nlev) throw fnirt_error("fnirt_clp::IntLambda: Out-of-range value of lev");
    return(static_cast<double>(biaslambda[lev-1]));
  }
  double MatchingPointsLambda(unsigned int lev=1) const {
    if (lev < 1 || lev > nlev) throw fnirt_error("fnirt_clp::MatchingPointsLambda: Out-of-range value of lev");
    return(static_cast<double>(mpl_lambda[lev-1]));
  }
  const vector<unsigned int> IntKnotSpacing(unsigned int lev) const {
    if (bf != Spline) throw fnirt_error("fnirt_clp::IntKnotSpacing: Knot spacing not defined for DCT basis");    
    if (lev < 1 || lev > nlev) throw fnirt_error("fnirt_clp::IntKnotSpacing: Out-of-range value of lev");
    return(bias_ksp[lev-1]);
  }
  const vector<unsigned int> IntDCTOrder(unsigned int lev) const {
    if (bf != DCT) throw fnirt_error("fnirt_clp::IntDCTOrder: Order not defined for Spline basis");    
    if (lev < 1 || lev > nlev) throw fnirt_error("fnirt_clp::IntDCTOrder: Out-of-range value of lev");
    return(bias_dco[lev-1]);
  }
  bool Verbose() const {return(verbose);}
  unsigned int Debug() const {return(debug);}
};

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Declarations of templated global functions
//
///////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////
//
// Declarations of global functions used by fnirt
//
///////////////////////////////////////////////////////////////////////////////////////////////

boost::shared_ptr<fnirt_clp> parse_fnirt_command_line(unsigned int   narg,
                                                      char           *args[]);

bool constrain_warpfield(const SSD_fnirt_CF&   cf,
                         const fnirt_clp&      clp,
                         unsigned int          max_try);

std::vector<boost::shared_ptr<BASISFIELD::basisfield> > init_warpfield(const fnirt_clp&  clp);

boost::shared_ptr<IntensityMapper> init_intensity_mapper(const fnirt_clp&  clp);

void make_field_with_same_fov(const BASISFIELD::splinefield&                   in,
                              boost::shared_ptr<BASISFIELD::basisfield>        out);

boost::shared_ptr<NEWIMAGE::volume<char> > make_mask(const string&                   fname, 
                                                     MaskType                        mt, 
                                                     const NEWIMAGE::volume<float>&  ima, 
                                                     bool                            impf, 
                                                     double                          impv);

bool trying_to_register_to_self(const string&                    ref_fname,
				const NEWIMAGE::volume<float>&   ref,
				const string&                    obj_fname,
				const NEWIMAGE::volume<float>&   obj,
                                const NEWMAT::Matrix&            aff);

void write_self_results(const fnirt_clp&                clp,
                        const NEWIMAGE::volume<float>&  ref);

void set_nlpars(NonlinParam&  nlp);
double spmlike_mean(NEWIMAGE::volume<float>&  ima);
bool is_identity(const NEWMAT::Matrix& A, double prec=1e-8);
bool check_exist(const std::string& fname);
std::string existing_conf_file(const std::string& cfname);
std::string existing_ref_fname(const std::string& ref_fname);
  
} // End namespace FNIRT

#endif // end #ifndef fnirtfns_h
