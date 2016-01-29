// Definitions of global functions used by fnirt
//
// fnirtfns.cpp
//
// Jesper Andersson and Matthew Webster, FMRIB Image Analysis Group
//
//
/*    Copyright (C) 2012-2015 University of Oxford  */

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
#include "warpfns/fnirt_file_reader.h"
#include "warpfns/point_list.h"
#include "intensity_mappers.h"
#include "fnirtfns.h"

using namespace std;
using namespace BASISFIELD;

namespace FNIRT {

fnirt_clp::fnirt_clp(const Utilities::Option<string>&                     pref,
                     const Utilities::Option<string>&                     pobj,
                     const Utilities::Option<string>&                     paff,
                     const Utilities::Option<string>&                     pinwarp,
                     const Utilities::Option<string>&                     pin_int,
                     const Utilities::Option<string>&                     pcoef,
                     const Utilities::Option<string>&                     pobjo,
                     const Utilities::Option<string>&                     pfieldo,
                     const Utilities::Option<string>&                     pjaco,
                     const Utilities::Option<string>&                     prefo,
                     const Utilities::Option<string>&                     pinto,
                     const Utilities::Option<string>&                     plogo,
                     const Utilities::Option<string>&                     prefm,
                     const Utilities::Option<string>&                     pobjm,
                     const Utilities::Option<string>&                     pref_pl,
                     const Utilities::Option<string>&                     pobj_pl,
                     const Utilities::Option<vector<int> >&               puserefm,
                     const Utilities::Option<vector<int> >&               puseobjm,
		     const Utilities::Option<int>&                        primf,
		     const Utilities::Option<int>&                        poimf,
		     const Utilities::Option<float>&                      primv,
		     const Utilities::Option<float>&                      poimv,
                     const Utilities::Option<string>&                     pcf,
                     const Utilities::Option<string>&                     pbf,
                     const Utilities::Option<string>&                     pnlm,
                     const Utilities::Option<vector<int> >&               pmi,
                     const Utilities::Option<vector<int> >&               pss,
                     const Utilities::Option<vector<float> >&             pwres,
                     const Utilities::Option<int>&                        pspordr,
	             const Utilities::Option<vector<float> >&             pofwhm,
                     Utilities::Option<vector<float> >&                   prfwhm,
		     const Utilities::Option<string>&                     pregmod,
                     Utilities::Option<vector<float> >&                   plambda,
                     const Utilities::Option<int>&                        pssqlambda,
                     const Utilities::Option<float>&                      pmpl_lambda,
		     const Utilities::Option<vector<float> >&             pjacrange,
                     const Utilities::Option<int>&                        puserefderiv,
                     const Utilities::Option<string>&                     pintmod,
		     const Utilities::Option<vector<int> >&               pestint,
                     const Utilities::Option<int>&                        pintord,
                     const Utilities::Option<vector<float> >&             pbiasres,
                     const Utilities::Option<string>&                     pbiasregmod,
                     const Utilities::Option<float>&                      pbiaslambda,
                     const Utilities::Option<bool>&                       pverbose,
                     const Utilities::Option<int>&                        pdebug,
                     const Utilities::Option<string>&                     p_hess_prec,
                     const Utilities::Option<string>&                     p_interp_type)
  : ref(pref.value()), obj(pobj.value()), inwarp(pinwarp.value()), in_int(pin_int.value()), coef(pcoef.value()), objo(pobjo.value()), 
    fieldo(pfieldo.value()), jaco(pjaco.value()), refo(prefo.value()), into(pinto.value()), logo(plogo.value()), 
    refm(prefm.value()), objm(pobjm.value()), ref_pl(pref_pl.value()), obj_pl(pobj_pl.value()), rimf((primf.value()==0) ? false : true), 
    rimv(primv.value()), oimf((poimf.value()==0) ? false : true), oimv(poimv.value()), spordr(static_cast<unsigned int>(pspordr.value())), 
    ssqlambda((pssqlambda.value()==0) ? false : true), jacrange(pjacrange.value()), userefderiv((puserefderiv.value()==0) ? false : true), 
    verbose(pverbose.value()), debug(static_cast<unsigned int>(pdebug.value()))
{
  // Parse and assert input

  // Check that image files exist and are valid
  NEWIMAGE::volume<float> vref, vobj, vany;
  if (!(ref = existing_ref_fname(ref)).size()) throw fnirt_error(string("fnirt_clp: No file matching --ref=")+ref+string(" found"));
  read_volume_hdr_only(vref,ref);
  if (refm.length()) {
    if (!(refm = existing_ref_fname(refm)).size()) throw fnirt_error(string("fnirt_clp: No file matching --refmask=")+refm+string(" found"));
    read_volume_hdr_only(vany,refm);
    if (vref.xsize() != vany.xsize() || vref.ysize() != vany.ysize() || vref.zsize() != vany.zsize() ||
        vref.xdim() != vany.xdim() || vref.ydim() != vany.ydim() || vref.zdim() != vany.zdim()) {
      throw fnirt_error(string("fnirt_clp: --ref=")+ref+string(" and --refmask=")+refm+string(" have different dimensions."));
    }
  }
  NEWIMAGE::read_volume_hdr_only(vobj,obj);     // Throws if there is a problem
  if (objm.length()) {
    NEWIMAGE::read_volume_hdr_only(vany,objm);  // Throws if there is a problem
    if (vobj.xsize() != vany.xsize() || vobj.ysize() != vany.ysize() || vobj.zsize() != vany.zsize() ||
        vobj.xdim() != vany.xdim() || vobj.ydim() != vany.ydim() || vobj.zdim() != vany.zdim()) {
      throw fnirt_error(string("fnirt_clp: --in=")+obj+string(" and --inmask: ")+objm+string(" have different dimensions."));
    }
  }

  // Check that either two or zero point-lists have been defined.
  if ((ref_pl.length() && !obj_pl.length()) || (!ref_pl.length() && obj_pl.length())) {
    throw fnirt_error("fnirt_clp: Must define zero or two point lists");
  }
  if (ref_pl.length()) {  // This means there are two point lists
    try {
      NEWIMAGE::PointList rpl(ref_pl,ref);
      NEWIMAGE::PointList opl(obj_pl,obj);
      if (rpl.NPoints() != opl.NPoints()) {
        throw fnirt_error(string("fnirt_clp: Mismatch in number of points in lists ")+ref_pl+string(" and ")+obj_pl);
      }
    }
    catch (std::exception& error) {
      throw fnirt_error(string("fnirt_clp: Error when reding point-list. Message was: ")+string(error.what()));
    }
  }

  // Read and assert affine matrix
  if (paff.value().empty() && pinwarp.value().empty()) aff = NEWMAT::IdentityMatrix(4);
  else if (!paff.value().empty()) {
    if (!pinwarp.value().empty()) { // If we have both affine and non-linear starting guess
      NEWIMAGE::FnirtFileReader    reader;
      try { reader.Read(pinwarp.value()); }
      catch (...) { throw fnirt_error(string("fnirt_clp: Problems reading initial warp file ")+pinwarp.value()); }
      NEWMAT::Matrix tmp_aff = reader.AffineMat();
      if ((tmp_aff-NEWMAT::IdentityMatrix(4)).MaximumAbsoluteValue() > 1e-6) { // This means we have a potential conflict
        throw fnirt_error("fnirt_clp: Dual affine transform (in --aff and --inwarp)");
      }
    }
    ifstream  ifs(paff.value().c_str());
    if (!ifs.fail()) {
      aff = MISCMATHS::read_ascii_matrix(ifs);
    }
    if (ifs.fail() || aff.Nrows() != 4 || aff.Ncols() != 4) {
      throw fnirt_error(string("fnirt_clp: Invalid affine matrix ")+paff.value()+string(" passed as argument to --aff"));
    }
  }
  else { // This means we (only) have a non-linear initial guess
    NEWIMAGE::FnirtFileReader    reader;
    try { reader.Read(pinwarp.value()); }
    catch (...) { throw fnirt_error(string("fnirt_clp: Problems reading initial warp file ")+pinwarp.value()); }
    aff = reader.AffineMat();
  }


  // Assert the categorical parameters
  if (p_hess_prec.value() == "float") hess_prec = BFMatrixFloatPrecision;
  else if (p_hess_prec.value() == "double") hess_prec = BFMatrixDoublePrecision;
  else throw fnirt_error("fnirt_clp: --numprec takes values float or double");
  if (p_interp_type.value() == "linear") interp_type = LinearInterp;
  else if (p_interp_type.value() == "spline") interp_type = SplineInterp;
  else throw fnirt_error("fnirt_clp: --interp takes values linear or spline");
  if (debug > 3) throw fnirt_error("fnirt_clp: --debug takes values 0, 1, 2 or 3");
  if (pcf.value() == "ssd") cf = SSD;
  else throw fnirt_error("fnirt_clp: Invalid cost-function option");
  if (pbf.value() == "spline") bf = Spline;
  else if (pbf.value() == "dct") bf = DCT;
  else throw fnirt_error("fnirt_clp: Invalid basis option");
  if (pnlm.value() == "lm") nlm = NL_LM;
  else if (pnlm.value() == "scg") nlm = NL_SCG;
  else throw fnirt_error("fnirt_clp: Invalid minimisation option");
  if (pregmod.value() == "membrane_energy") regmod = MembraneEnergy;
  else if (pregmod.value() == "bending_energy") regmod = BendingEnergy;
  else throw fnirt_error("fnirt_clp: Invalid regularisation type "+pregmod.value()+string(" passed as argument to --regmod"));
  if (spordr < 2 || spordr > 3) throw fnirt_error("fnirt_clp: Spline order (--splineorder) must be 2 or 3");

  // Make sure jacobian range isn't silly
  if (jacrange.size() == 1) {
    if (jacrange[0] == -1) { // -1 used to indicate no Jacobian constraints
      jacrange[0] = -1e6;
      jacrange.push_back(1e6);
    }
    else throw fnirt_error("fnirt_clp: --jacrange takes 2 arguments, or a single -1");
  }
  else if (jacrange.size() == 2 && (jacrange[0] >= jacrange[1])) throw fnirt_error("fnirt_clp: Invalid jacobian range");
  else if (jacrange.size() != 2) throw fnirt_error("fnirt_clp: --jacrange takes 2 arguments, or a single -1");

  // Make sure subsampling is by a power of 2
  for (int i=0; i<int(pss.value().size()); i++) {
    bool valid = false;
    for (int j=1; j<64; j*=2) {
      if ((pss.value())[i] == j) valid = true;
    }
    if (!valid) throw fnirt_error("fnirt_clp: Subsampling must be by a power of 2 at all levels");
  }
  // Make sure we have correspondence between subsampling
  // and the various parameters that should be defined for
  // each level of subsampling
  ss.resize(nlev = pss.value().size());
  for (unsigned int i=0; i<pss.value().size(); i++) ss[i] = static_cast<unsigned int>(pss.value()[i]);
  if (pmi.value().size() != nlev && pmi.value().size() != 1) {
    throw fnirt_error("fnirt_clp: Mismatch between --subsamp and --miter options");
  }
  if (pofwhm.value().size() != nlev && pofwhm.value().size() != 1) {
    throw fnirt_error("fnirt_clp: Mismatch between --subsamp and --infwhm options");
  }
  if (prfwhm.value().size() != nlev && prfwhm.value().size() != 1 && prfwhm.value().size() != 0) {
    throw fnirt_error("fnirt_clp: Mismatch between --subsamp and --reffwhm options");
  }
  if (plambda.value().size() != nlev && plambda.value().size() != 1 && plambda.value().size() != 0) {
    throw fnirt_error("fnirt_clp: Mismatch between --subsamp and --lambda options");
  }
  if (pestint.value().size() != nlev && pestint.value().size() != 1) {
    throw fnirt_error("fnirt_clp: Mismatch between --subsamp and --estint options");
  }
  if (puserefm.value().size() != nlev && puserefm.value().size() != 1) {
    throw fnirt_error("fnirt_clp: Mismatch between --subsamp and --applyrefmask options");
  }
  if (puseobjm.value().size() != nlev && puseobjm.value().size() != 1) {
    throw fnirt_error("fnirt_clp: Mismatch between --subsamp and --applyinmask options");
  }

  // Assign
  if (pmi.value().size() == nlev) {
    mi.resize(pmi.value().size()); 
    for (unsigned int i=0; i<pmi.value().size(); i++) mi[i] = static_cast<unsigned int>(pmi.value()[i]);
  }
  else {mi = vector<unsigned int>(nlev); for (unsigned int i=0; i<nlev; i++) mi[i] = (pmi.value())[0];}

  if (pofwhm.value().size() == nlev) ofwhm = pofwhm.value();
  else {ofwhm = vector<float>(nlev); for (unsigned int i=0; i<nlev; i++) ofwhm[i] = (pofwhm.value())[0];}

  if (prfwhm.value().size() == nlev) rfwhm = prfwhm.value();
  else if (prfwhm.value().size() == 1) {rfwhm = vector<float>(nlev); for (unsigned int i=0; i<nlev; i++) rfwhm[i] = (prfwhm.value())[0];}
  else {
    rfwhm = ofwhm;
    prfwhm.set_T(pofwhm.value());
  }

  if (plambda.value().size() == nlev) lambda = plambda.value();
  else {
    lambda = vector<float>(nlev);
    if (plambda.value().size() == 1) {for (unsigned int i=0; i<nlev; i++) lambda[i] = (plambda.value())[0];}
    else {
      // Here we have to assign values for lambda
      // based on what regularisation model has been
      // specified and wether to weight by ssq
      if (ssqlambda) {
        if (regmod == BendingEnergy) {for (unsigned int i=0; i<nlev; i++) lambda[i] = ss[i] * 30.0;}
	else if (regmod == MembraneEnergy) {for (unsigned int i=0; i<nlev; i++) lambda[i] = ss[i] * 2.5;}
      }
      else {
	if (regmod == BendingEnergy) {for (unsigned int i=0; i<nlev; i++) lambda[i] = ss[i] * 15000.0;} 
	else if (regmod == MembraneEnergy) {for (unsigned int i=0; i<nlev; i++) lambda[i] = ss[i] * 1500;}
      }
      plambda.set_T(lambda);
    }
  }

  estint.resize(nlev);
  if (pestint.value().size() == nlev) for (unsigned int i=0; i<nlev; i++) estint[i] = (pestint.value()[i] ? true : false);
  else if (pestint.value().size() == 1) for (unsigned int i=0; i<nlev; i++) estint[i] = (pestint.value()[0] ? true : false);

  userefmask.resize(nlev);
  if (puserefm.value().size() == nlev) for (unsigned int i=0; i<nlev; i++) userefmask[i] = (puserefm.value()[i] ? true : false);
  else if (puserefm.value().size() == 1) for (unsigned int i=0; i<nlev; i++) userefmask[i] = (puserefm.value()[0] ? true : false);

  useobjmask.resize(nlev);
  if (puseobjm.value().size() == nlev) for (unsigned int i=0; i<nlev; i++) useobjmask[i] = (puseobjm.value()[i] ? true : false);
  else if (puseobjm.value().size() == 1) for (unsigned int i=0; i<nlev; i++) useobjmask[i] = (puseobjm.value()[0] ? true : false);

  //
  // Make sure we are not trying to estimate intensities with Scaled Conjugate-gradient.
  //
  if (nlm == NL_SCG) for (unsigned int i=0; i<nlev; i++) if (estint[i]) throw fnirt_error("fnirt_clp: Cannot estimate intensity mapping with Scaled Conjugate Gradient method");

  //
  // Assert and parse parameters pertaining to intensity mapping.
  //
  if (!pintmod.value().compare(string("none"))) imt = NONE;
  else if (!pintmod.value().compare(string("global_linear"))) imt = GLOBAL_LINEAR;
  else if (!pintmod.value().compare(string("global_non_linear"))) imt = GLOBAL_NON_LINEAR;
  else if (!pintmod.value().compare(string("local_linear"))) imt = LOCAL_LINEAR;
  else if (!pintmod.value().compare(string("global_non_linear_with_bias"))) imt = LOCAL_BIAS_WITH_GLOBAL_NON_LINEAR;
  else if (!pintmod.value().compare(string("local_non_linear"))) imt = LOCAL_NON_LINEAR;
  else throw fnirt_error(string("fnirt_clp: Argument to --intmod is unknown scaling option: ")+pintmod.value());

  if (pbiasregmod.value() == string("bending_energy")) biasregmod = BendingEnergy;
  else throw fnirt_error(string("fnirt_clp: Argument to --biasregmod is unknown option: ")+pbiasregmod.value());

  biaslambda = vector<float>(nlev,pbiaslambda.value());
  mpl_lambda = vector<float>(nlev,pmpl_lambda.value());

  if (imt == GLOBAL_NON_LINEAR && (intord = static_cast<unsigned int>(pintord.value())) > 10) {
    throw fnirt_error("fnirt_clp: Argument to --intorder, cannot use polynomial of order > 10 for intensity mapping");
  }
  else if (imt == LOCAL_BIAS_WITH_GLOBAL_NON_LINEAR && (intord = static_cast<unsigned int>(pintord.value())) > 5) {
    throw fnirt_error("fnirt_clp: Argument to --intorder, cannot use polynomial of order > 5 for intensity mapping when also modelling intensity bias");
  }
  else if (imt == LOCAL_NON_LINEAR && (intord = static_cast<unsigned int>(pintord.value())) > 5) {
    throw fnirt_error("fnirt_clp: Argument to --intorder, cannot use polynomial of order > 5 for intensity mapping when also modelling local non-linear mappings");
  }
  
  double pxs[] = {vref.xdim(), vref.ydim(), vref.zdim()};
  int    dim[] = {vref.xsize(), vref.ysize(), vref.zsize()};

  //
  // Convert mm to voxels/order for splines/DCT for displacement fields
  //
  vector<float>  tmpwres(3,0);
  if (pwres.value().size() == 3) tmpwres = pwres.value();
  else if (pwres.value().size() == 1) {for (int i=0; i<3; i++) tmpwres[i] = (pwres.value())[0];}
  else throw fnirt_error("fnirt_clp: --warpres option must be 3x1 vector or scalar");
  if (bf == Spline) {
    vector<unsigned int>     tmpksp(3,0);
    for (int i=0; i<3; i++) {
      tmpksp[i] = max(static_cast<unsigned int>(1),static_cast<unsigned int>(floor(tmpwres[i]/pxs[i] + 0.5)));
      // cout << "tmpwres[" << i << "]=" << tmpwres[i] << ", pxs[" << i << "]=" << pxs[i] << ", tmpksp[" << i << "]=" << tmpksp[i] << endl;
    }
    ksp = vector<vector<unsigned int> >(nlev);
    for (unsigned int i=0; i<nlev; i++) ksp[i] = tmpksp;
    if (!(tmpksp[0]/ss[nlev-1])) throw fnirt_error("fnirt_clp: Incompatible combination of --subsamp, --warpres and voxel size in x-direction");
    if (!(tmpksp[1]/ss[nlev-1])) throw fnirt_error("fnirt_clp: Incompatible combination of --subsamp, --warpres and voxel size in y-direction");
    if (!(tmpksp[2]/ss[nlev-1])) throw fnirt_error("fnirt_clp: Incompatible combination of --subsamp, --warpres and voxel size in z-direction");
  }
  else if (bf == DCT) {
    vector<int> tmpdco(3,0);
    for (int i=0; i<3; i++) tmpdco[i] = static_cast<int>(floor(pxs[i]*dim[i]/tmpwres[i] + 0.5));
    dco = vector<vector<unsigned int> >(nlev);
    for (unsigned int i=0; i<nlev; i++) {
      for (int j=0; j<3; j++) {
        dco[i][j] = static_cast<unsigned int>(floor(float(tmpdco[j])/float(ss[i]/ss[nlev-1]) + 0.5));
      }
    }
  }

  //
  // Convert mm to voxels/order for splines/DCT for bias-fields,
  // taking into account that we want to keep the resolution (in mm)
  // constant across the subsamplings.
  //
  vector<float>  tmpbres(3,0);
  if (pbiasres.value().size() == 3) tmpbres = pbiasres.value();
  else if (pbiasres.value().size() == 1) {for (int i=0; i<3; i++) tmpbres[i] = (pbiasres.value())[0];}
  else throw fnirt_error("fnirt_clp: --biasres option must be 3x1 vector or scalar");
  if (bf == Spline) {
    vector<unsigned int>  tmpksp(3,0);       // Knot-spacing in voxels for lowest resolution (most sub-sampling)
    for (int i=0; i<3; i++) tmpksp[i] = static_cast<unsigned int>(floor(double(tmpbres[i] / double(pxs[i] * ss[0])) + 0.5));
    bias_ksp.resize(nlev);
    for (unsigned int l=0; l<nlev; l++) {
      bias_ksp[l].resize(3);
      for (int i=0; i<3; i++) {
        bias_ksp[l][i] = (ss[0]/ss[l]) * tmpksp[i];
      }
      // cout << "level = " << l << ", bias_ksp[0] = " << bias_ksp[l][0] << ", bias_ksp[1] = " << bias_ksp[l][1] << ", bias_ksp[2] = " << bias_ksp[l][2] << endl;
    }
  }
  else if (bf == DCT) {
    throw fnirt_error("Sorry, not yet implemented");
  }  
}

boost::shared_ptr<fnirt_clp> parse_fnirt_command_line(unsigned int   narg,
                                                      char           *args[])
{
  // These are the specifications of all the command-line options

  Utilities::Option<string> refname(string("--ref"), string(""),
      string("\tname of reference image"), true, Utilities::requires_argument);

  Utilities::Option<string> objname(string("--in"), string(""),
      string("\tname of input image"), true, Utilities::requires_argument);

  Utilities::Option<string> affname(string("--aff"), string(""),
      string("\tname of file containing affine transform"), false, Utilities::requires_argument);

  Utilities::Option<string> inwarpname(string("--inwarp"), string(""),
      string("name of file containing initial non-linear warps"), false, Utilities::requires_argument);

  Utilities::Option<string> in_intname(string("--intin"), string(""),
      string("\tname of file/files containing initial intensity mapping"), false, Utilities::requires_argument);

  Utilities::Option<string> coefname(string("--cout"), string(""),
      string("\tname of output file with field coefficients"), false, Utilities::requires_argument);

  Utilities::Option<string> imoutname(string("--iout"), string(""),
      string("\tname of output image"), false, Utilities::requires_argument);

  Utilities::Option<string> fieldoutname(string("--fout"), string(""),
      string("\tname of output file with field"), false, Utilities::requires_argument);

  Utilities::Option<string> jacname(string("--jout"), string(""),
      string("\tname of file for writing out the Jacobian of the field (for diagnostic or VBM purposes)"),
      false,Utilities::requires_argument);

  Utilities::Option<string> refoutname(string("--refout"), string(""),
      string("name of file for writing out intensity modulated --ref (for diagnostic purposes)"),
      false,Utilities::requires_argument);

  Utilities::Option<string> intoutname(string("--intout"), string(""),
      string("name of files for writing information pertaining to intensity mapping"),false,Utilities::requires_argument);

  Utilities::Option<string> logoutname(string("--logout"),string(""),
			    string("Name of log-file"),false,Utilities::requires_argument);

  Utilities::Option<string> configname(string("--config"),string(""),
			    string("Name of config file specifying command line arguments"),false,Utilities::requires_argument);

  Utilities::Option<string> costfunction(string("--cf"), string("ssd"),
      string("cost-function to minimise"), false, Utilities::requires_argument);

  Utilities::Option<string> basis(string("--basis"), string("spline"),
      string("\tbasis used to model field [spline | dct]"), false, Utilities::requires_argument);

  Utilities::Option<string> refmaskname(string("--refmask"), string(""),
      string("name of file with mask in reference space"), false, Utilities::requires_argument);

  Utilities::Option<string> objmaskname(string("--inmask"), string(""),
      string("name of file with mask in input image space"),false, Utilities::requires_argument);

  // The point-lists are hidden for the time being
  Utilities::HiddenOption<string> refpointlistname(string("--refpointlist"), string(""),
      string("name of file with points/coordinates in reference space"), false, Utilities::requires_argument);

  Utilities::HiddenOption<string> objpointlistname(string("--inpointlist"), string(""),
      string("name of file with points/coordinates in input image space"),false, Utilities::requires_argument);

  Utilities::HiddenOption<float> mpl_lambda(string("--pointlistlambda"),1.0,
      string("Weight of landmark distances relative to intensities, default 1"), 
      false, Utilities::requires_argument);

  vector<int> applyrefmaskdefault(4,1);
  Utilities::Option<vector<int> > applyrefmask(string("--applyrefmask"),applyrefmaskdefault,
      string("Use specified refmask if set, default 1 (true)"),false,Utilities::requires_argument);

  vector<int> applyobjmaskdefault(4,1);
  Utilities::Option<vector<int> > applyobjmask(string("--applyinmask"),applyobjmaskdefault,
      string("Use specified inmask if set, default 1 (true)"),false,Utilities::requires_argument);

  Utilities::Option<int> imprefflag(string("--imprefm"), 1,
      string("If =1, use implicit masking based on value in --ref image. Default =1"),false, Utilities::requires_argument);

  Utilities::Option<int> impobjflag(string("--impinm"), 1,
      string("If =1, use implicit masking based on value in --in image, Default =1"),false, Utilities::requires_argument);

  Utilities::Option<float> imprefval(string("--imprefval"), 0.0,
      string("Value to mask out in --ref image. Default =0.0"),false, Utilities::requires_argument);

  Utilities::Option<float> impobjval(string("--impinval"), 0.0,
      string("Value to mask out in --in image. Default =0.0"),false, Utilities::requires_argument);

  Utilities::Option<string> minimisationmethod(string("--minmet"), string("lm"),
      string("non-linear minimisation method [lm | scg] (Levenberg-Marquardt or Scaled Conjugate Gradient)"), 
      false, Utilities::requires_argument);

  vector<int> maxiterdefault(4,0);
  maxiterdefault[0] = 5; maxiterdefault[1] = 5; maxiterdefault[2] = 5; maxiterdefault[3] = 5; 
  Utilities::Option<vector<int> > maxiter(string("--miter"), maxiterdefault,
      string("\tMax # of non-linear iterations, default 5,5,5,5"), false, Utilities::requires_argument);

  vector<int> subsamplingdefault(4,0);
  subsamplingdefault[0] = 4; subsamplingdefault[1] = 2; subsamplingdefault[2] = 1; subsamplingdefault[3] = 1; 
  Utilities::Option<vector<int> > subsampling(string("--subsamp"), subsamplingdefault,
      string("sub-sampling scheme, default 4,2,1,1"), false, Utilities::requires_argument);

  vector<float> warpresdefault(3,10.0);
  Utilities::Option<vector<float> > warpres(string("--warpres"), warpresdefault,
      string("(approximate) resolution (in mm) of warp basis in x-, y- and z-direction, default 10,10,10"), false, Utilities::requires_argument);

  Utilities::Option<int> splineorder(string("--splineorder"), 3,
      string("Order of spline, 2->Quadratic spline, 3->Cubic spline. Default=3"),false, Utilities::requires_argument);

  vector<float> objsmoothdefault(4,0);
  objsmoothdefault[0] = 6.0; objsmoothdefault[1] = 4.0; objsmoothdefault[2] = 2.0; objsmoothdefault[3] = 2.0;
  Utilities::Option<vector<float> > objsmoothing(string("--infwhm"),objsmoothdefault,
      string("FWHM (in mm) of gaussian smoothing kernel for input volume, default 6,4,2,2"), false, Utilities::requires_argument);

  vector<float> refsmoothdefault(4,0);
  refsmoothdefault[0] = 4.0; refsmoothdefault[1] = 2.0; refsmoothdefault[2] = 0.0; refsmoothdefault[3] = 0.0;
  Utilities::Option<vector<float> > refsmoothing(string("--reffwhm"),refsmoothdefault,
      string("FWHM (in mm) of gaussian smoothing kernel for ref volume, default 4,2,0,0"), false, Utilities::requires_argument);

  Utilities::Option<string> regularisationmodel(string("--regmod"),string("bending_energy"),
				                string("Model for regularisation of warp-field [membrane_energy bending_energy], default bending_energy"),
				                false,Utilities::requires_argument); 

  vector<float> lambdadefault(0);
  Utilities::Option<vector<float> > lambda(string("--lambda"),lambdadefault,
      string("Weight of regularisation, default depending on --ssqlambda and --regmod switches. See user documentation."), 
      false, Utilities::requires_argument);

  Utilities::Option<int> ssqlambda(string("--ssqlambda"),true,string("If set (=1), lambda is weighted by current ssq, default 1"),
      false,Utilities::requires_argument);

  vector<float> jacrangedefault(2,0);
  jacrangedefault[0] = 0.01; jacrangedefault[1] = 100.0;   // Basically non-negativity
  Utilities::Option<vector<float> > jacrange(string("--jacrange"),jacrangedefault,
      string("Allowed range of Jacobian determinants, default 0.01,100.0"), 
      false, Utilities::requires_argument);

  Utilities::Option<int> userefderiv(string("--refderiv"),0,string("If =1, ref image is used to calculate derivatives. Default =0"),
      false,Utilities::requires_argument);

  Utilities::Option<string> intensitymodel(string("--intmod"),string("global_non_linear_with_bias"),
      string("Model for intensity-mapping [none global_linear global_non_linear local_linear global_non_linear_with_bias local_non_linear]"),
      false,Utilities::requires_argument); 

  int intensityorderdefault = 5;
  Utilities::Option<int> intensityorder(string("--intorder"),intensityorderdefault,
      string("Order of polynomial for mapping intensities, default 5"),false,Utilities::requires_argument);
   
  vector<float> biasresdefault(3,50.0);
  Utilities::Option<vector<float> > biasfieldres(string("--biasres"),biasresdefault,
      string("Resolution (in mm) of bias-field modelling local intensities, default 50,50,50"),false,Utilities::requires_argument);

  Utilities::Option<float> biasfieldlambda(string("--biaslambda"),1e4,
      string("Weight of regularisation for bias-field, default 10000"),false,Utilities::requires_argument);

  Utilities::Option<string> biasfieldregmod(string("--biasregmod"),string("bending_energy"),
      string("Model for regularisation of biasfield, default bending_energy"),false,Utilities::requires_argument);

  vector<int> estintdefault(4,1);
  estintdefault[3] = 0;
  Utilities::Option<vector<int> > estimateintensity(string("--estint"),estintdefault,
      string("Estimate intensity-mapping if set, default 1 (true)"),false,Utilities::requires_argument);

  Utilities::Option<string> numprec(string("--numprec"),string("double"),
      string("Precision for representing Hessian, double or float. Default double"),false,Utilities::requires_argument);

  Utilities::Option<string> interpolation(string("--interp"),string("linear"),
      string("Image interpolation model, linear or spline. Default linear"),false,Utilities::requires_argument);

  Utilities::Option<string> configfile(string("--config"),string(""),
      string("Name of configuration field with settings for some/all fnirt parameters"),false,Utilities::requires_argument);

  Utilities::Option<bool> help(string("-h,--help"), false,
      string("display help info"), false, Utilities::no_argument);

  Utilities::Option<bool> verbose(string("-v,--verbose"), false,
      string("Print diagnostic information while running"), false, Utilities::no_argument);

  Utilities::HiddenOption<int> debug(string("--debug"), 0,
      string("Save debug information while running, levels 0 (no info), 1 (some info), 2 (little more info) or 3 (LOTS of info)"), false, Utilities::requires_argument);

  // Some explanatory text

  string title = "fnirt";
  string examples = string("fnirt --ref=<some template> --in=<some image>\n") +
                    string("fnirt --ref=<some template> --in=<some image> --infwhm=8,4,2 --subsamp=4,2,1 --warpres=8,8,8");

  // Create and load options-parser

  Utilities::OptionParser options(title, examples);

  try {
    // Load parser

    options.add(refname);
    options.add(objname);
    options.add(affname);
    options.add(inwarpname);
    options.add(in_intname);
    options.add(coefname);
    options.add(imoutname);
    options.add(fieldoutname);
    options.add(jacname);
    options.add(refoutname);
    options.add(intoutname);
    options.add(logoutname);
    options.add(configname);
    options.add(refmaskname);
    options.add(objmaskname);
    options.add(refpointlistname);
    options.add(objpointlistname);
    options.add(applyrefmask);
    options.add(applyobjmask);
    options.add(imprefflag);
    options.add(impobjflag);
    options.add(imprefval);
    options.add(impobjval);
    // options.add(costfunction);         Cost-function option hidden for now
    // options.add(basis);                Basis-set hidden for now
    options.add(minimisationmethod);      
    options.add(maxiter);
    options.add(subsampling);
    options.add(warpres);
    options.add(splineorder);
    options.add(objsmoothing);
    options.add(refsmoothing);
    options.add(regularisationmodel);
    options.add(lambda);
    options.add(ssqlambda);
    options.add(mpl_lambda);
    options.add(jacrange);
    options.add(userefderiv);
    options.add(intensitymodel);
    options.add(intensityorder);
    options.add(biasfieldres);
    // options.add(biasfieldregmod);      Regularisation model for bias field hidden for now
    options.add(biasfieldlambda);
    options.add(estimateintensity);
    options.add(numprec);
    options.add(interpolation);
    options.add(verbose);
    options.add(debug);
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
      final_configname = existing_conf_file(configname.value());
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

  boost::shared_ptr<fnirt_clp>   clp;
  try {
    clp = boost::shared_ptr<fnirt_clp>(new fnirt_clp(refname,objname,affname,inwarpname,in_intname,coefname,imoutname,fieldoutname,
                                                     jacname,refoutname,intoutname,logoutname,refmaskname,objmaskname,refpointlistname,
                                                     objpointlistname,applyrefmask,applyobjmask,imprefflag,impobjflag,imprefval,impobjval,costfunction,
                                                     basis,minimisationmethod,maxiter,subsampling,warpres,splineorder,objsmoothing,
                                                     refsmoothing,regularisationmodel,lambda,ssqlambda,mpl_lambda,jacrange,userefderiv,intensitymodel,
                                                     estimateintensity,intensityorder,biasfieldres,biasfieldregmod,
                                                     biasfieldlambda,verbose,debug,numprec,interpolation));
  }
  catch(fnirt_error& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  }

  // Write out a log-file containing all the parsed/expanded parameters
  try {
    std::ofstream logfs(clp->LogFname().c_str());
    logfs << refname << endl;
    logfs << objname << endl;
    if (affname.set()) logfs << affname << endl;
    if (inwarpname.set()) logfs << inwarpname << endl;
    if (in_intname.set()) logfs << in_intname << endl;
    if (coefname.set()) logfs << coefname << endl;
    if (imoutname.set()) logfs << imoutname << endl;
    if (fieldoutname.set()) logfs << fieldoutname << endl;
    if (jacname.set()) logfs << jacname << endl;
    if (refoutname.set()) logfs << refoutname << endl;
    if (intoutname.set()) logfs << intoutname << endl;
    if (logoutname.set()) logfs << logoutname << endl;
    if (refmaskname.set()) logfs << refmaskname << endl;
    if (objmaskname.set()) logfs << objmaskname << endl;
    // Uncomment these when pointlist functionality is released
    // if (refpointlistname.set()) logfs << refpointlistname << endl;
    // if (objpointlistname.set()) logfs << objpointlistname << endl;
    // logfs << mpl_lambda << endl;    
    logfs << imprefflag << endl;
    logfs << impobjflag << endl;
    logfs << imprefval << endl;
    logfs << impobjval << endl;
    logfs << subsampling << endl;
    logfs << maxiter << endl;
    logfs << objsmoothing << endl;
    logfs << refsmoothing << endl;
    logfs << lambda << endl;
    logfs << estimateintensity << endl;
    logfs << warpres << endl;
    logfs << splineorder << endl;
    logfs << ssqlambda << endl;
    logfs << jacrange << endl;
    logfs << regularisationmodel << endl;
    logfs << intensitymodel << endl;
    logfs << intensityorder << endl;
    logfs << biasfieldres << endl;
    logfs << biasfieldlambda << endl;
    logfs << numprec << endl;
    logfs << interpolation << endl;
    logfs << userefderiv << endl;    
    logfs.close();
  }
  catch (...) {
    cerr << "Error when writing log-file: " << clp->LogFname() << endl;
    throw;
  }

  return(clp);  
}

/////////////////////////////////////////////////////////////////////
//
// This routines will use the constrain_topology function from
// warpfns to make sure the Jacobians of the field is within
// the prescribed range. The break statement will be executed if
// subsequent calls to ForceJacobianRange produces identical
// results. This indicates a disagreement between the continous
// Jacobian calculated within fnirt and the discrete representation
// used by constrain_topology. Effectively, constrain_topology
// considers the field to already be within the Jacobian bounds
// and will make no further changes. This will/may happen for very
// high resolution fields with very small regions (single voxels)
// outside the permitted range.
//
/////////////////////////////////////////////////////////////////////

bool constrain_warpfield(const SSD_fnirt_CF&   cf,
                         const fnirt_clp&      clp,
                         unsigned int          max_try)
{
  std::pair<double,double> range = cf.JacobianRange();
  std::pair<double,double> last_range = range;
  if (clp.Verbose()) cout << "Jacobian range is " << range.first << " -- " << range.second << endl;

  unsigned int n_try=0;
  while ((range.first < clp.JacLowerBound() || range.second > clp.JacUpperBound()) && n_try < max_try) {
    if (clp.Verbose()) cout << "Forcing Jacobian range to " << clp.JacLowerBound() << " -- " << clp.JacUpperBound() << endl;
    cf.ForceJacobianRange(clp.JacLowerBound(),clp.JacUpperBound());
    range = cf.JacobianRange();
    if (clp.Verbose()) cout << "Jacobian range is " << range.first << " -- " << range.second << endl;
    if (std::fabs(last_range.first-range.first) < 1e-6 || std::fabs(last_range.second-range.second) < 1e-6) break;
    last_range = range;
    n_try++;
  } 

  if (range.first < clp.JacLowerBound() || range.second > clp.JacUpperBound()) return(false);

  return(true);  
}

vector<boost::shared_ptr<basisfield> > init_warpfield(const fnirt_clp&  clp)
{
  vector<boost::shared_ptr<basisfield> >    field(3);

  if (clp.InWarp().size()) {  // If there is a starting guess for the field
    NEWIMAGE::FnirtFileReader    reader;
    try {
      reader.Read(clp.InWarp());
    }
    catch (...) {
      throw fnirt_error(string("Problems reading initial warp file ")+clp.InWarp());
    }
    if (clp.RefSize() != reader.FieldSize() || clp.RefVxs() != reader.VoxelSize()) {
      // Different size fields. Check if it maybe is up/down-sampling.
      if ((std::fabs(reader.FieldSize()[0]*reader.VoxelSize()[0] - clp.RefSize()[0]*clp.RefVxs()[0]) > 1e-3) || 
          (std::fabs(reader.FieldSize()[1]*reader.VoxelSize()[1] - clp.RefSize()[1]*clp.RefVxs()[1]) > 1e-3) || 
          (std::fabs(reader.FieldSize()[2]*reader.VoxelSize()[2] - clp.RefSize()[2]*clp.RefVxs()[2]) > 1e-3)) {
        // This means FOV is different, i.e. not just up/down-sampling
	if (clp.RefSize() != reader.FieldSize()) throw fnirt_error(string("Field in file ")+clp.InWarp()+string(" not of same size as image in ")+clp.Ref());
        if (clp.RefVxs() != reader.VoxelSize()) throw fnirt_error(string("Field in file ")+clp.InWarp()+string(" has different voxel size from image in ")+clp.Ref());
      } 
      else { // We will assume this means that it is just a matter of up/down-sampling
        if (clp.Basis() == FNIRT::Spline) {
          if (clp.Verbose()) cout << "Resampling --inwarp field" << endl;
          for (unsigned int fi=0; fi<3; fi++) {
            field[fi] = boost::shared_ptr<splinefield>(new splinefield(clp.RefSize(),clp.RefVxs(),clp.FullResKsp(),clp.SplineOrder()));
            BASISFIELD::splinefield original = reader.FieldAsSplinefield(fi);
            make_field_with_same_fov(original,field[fi]);
	  }
	}
	else throw fnirt_error(string("FNIRT cannot zoom DCT --inwarp field ")+clp.InWarp());
      }
    }
    else if (clp.Basis() == FNIRT::Spline) {
      for (int i=0; i<3; i++) field[i] = boost::shared_ptr<splinefield>(new splinefield(reader.FieldAsSplinefield(i,clp.FullResKsp())));
    }
    else if (clp.Basis() == FNIRT::DCT) {  
      for (int i=0; i<3; i++) field[i] = boost::shared_ptr<dctfield>(new dctfield(reader.FieldAsDctfield(i,clp.DCTOrder(clp.NoLevels()))));
    }      
  }
  else { // If we are starting from scratch
    if (clp.Basis() == FNIRT::Spline) {
      for (int i=0; i<3; i++) field[i] = boost::shared_ptr<splinefield>(new splinefield(clp.RefSize(),clp.RefVxs(),clp.FullResKsp(),clp.SplineOrder()));
    }
    else if (clp.Basis() == FNIRT::DCT) {  
      for (int i=0; i<3; i++) field[i] = boost::shared_ptr<dctfield>(new dctfield(clp.RefSize(),clp.RefVxs(),clp.DCTOrder(clp.NoLevels())));
    }
  }
  return(field);
}

boost::shared_ptr<IntensityMapper> init_intensity_mapper(const fnirt_clp&  clp)
{
  boost::shared_ptr<IntensityMapper>  mymapper;  // Null pointer
  IntensityMapperReader               mapreader;

  if (clp.InInt().size()) mapreader.Read(clp.InInt());    // If initialisation file name has been specified

  if (clp.IntMapType() == NONE) {
    mymapper = boost::shared_ptr<IntensityMapper>(new SSDIntensityMapper);
  }
  else if (clp.IntMapType() == GLOBAL_LINEAR) {
    double se = 1.0;
    if (mapreader.HasGlobal()) {
      vector<double>  old_global = mapreader.GetGlobal();
      if (old_global.size()==1) se = old_global[0];
      else if (old_global.size()>1) se = old_global[1];
    }
    mymapper = boost::shared_ptr<IntensityMapper>(new SSDIntensityMapper(se));
  }
  else if (clp.IntMapType() == GLOBAL_NON_LINEAR) {
    vector<double>  se(clp.IntMapOrder(),0.0);       // Set all components to zero
    se[1] = 1.0;                                     // Set linear component to 1
    if (mapreader.HasGlobal()) {
      vector<double> old_global = mapreader.GetGlobal();
      if (old_global.size()==1) se[1] = old_global[0];
      else if (old_global.size()>1) for (unsigned int i=0; i<Min(se.size(),old_global.size()); i++) se[i] = old_global[i];
    }
    mymapper = boost::shared_ptr<IntensityMapper>(new SSDIntensityMapper(se));
  }
  else if (clp.IntMapType() == LOCAL_LINEAR) {
    boost::shared_ptr<basisfield>  fld;  // Null-pointer
    if (clp.Basis() == Spline) {
      if (mapreader.HasLocal()) {
        if (mapreader.LocalFieldSize()!=clp.RefSize() || mapreader.LocalFieldVoxelSize()!=clp.RefVxs()) { // This means dimensions are not identical
          if ((std::fabs(mapreader.LocalFieldSize()[0]*mapreader.LocalFieldVoxelSize()[0] - clp.RefSize()[0]*clp.RefVxs()[0]) > 1e-3) || 
              (std::fabs(mapreader.LocalFieldSize()[1]*mapreader.LocalFieldVoxelSize()[1] - clp.RefSize()[1]*clp.RefVxs()[1]) > 1e-3) || 
              (std::fabs(mapreader.LocalFieldSize()[2]*mapreader.LocalFieldVoxelSize()[2] - clp.RefSize()[2]*clp.RefVxs()[2]) > 1e-3)) { // This means FOV is different
            if (mapreader.LocalFieldSize()!=clp.RefSize()) throw fnirt_error("init_intensity_mapper: Intensity Mapper file has different matrix-size from --ref");
            if (mapreader.LocalFieldVoxelSize()!=clp.RefVxs()) throw fnirt_error("init_intensity_mapper: Intensity Mapper file has different voxel-size from --ref");
	  }
	  else { // This means FOV is the same, so we'll assume up/down-sampling
            fld = boost::shared_ptr<BASISFIELD::basisfield>(new BASISFIELD::splinefield(clp.RefSize(),clp.RefVxs(),clp.FullResIntKsp()));
	    std::vector<unsigned int>  ksp(3); // To avoid possible loss of resolution
            ksp[0] = static_cast<unsigned int>((clp.RefVxs()[0]/mapreader.LocalFieldVoxelSize()[0])*clp.FullResIntKsp()[0]);
            ksp[1] = static_cast<unsigned int>((clp.RefVxs()[1]/mapreader.LocalFieldVoxelSize()[1])*clp.FullResIntKsp()[1]);
            ksp[2] = static_cast<unsigned int>((clp.RefVxs()[2]/mapreader.LocalFieldVoxelSize()[2])*clp.FullResIntKsp()[2]);
            if (mapreader.NLocalFields()==1) {
	      if (clp.Verbose()) cout << "Resampling --intin field" << endl;
	      BASISFIELD::splinefield original = mapreader.GetLocalAsSingleSplinefield(ksp);
	      make_field_with_same_fov(original,fld);
	    }
	    else if (mapreader.NLocalFields()>1) {
	      if (clp.Verbose()) cout << "Resampling --intin field" << endl;
              BASISFIELD::splinefield original = mapreader.GetLocalAsSingleSplinefield(ksp,1);
	      make_field_with_same_fov(original,fld);
	    }
	  }
	}
	else { // Dimensions are identical
	  if (mapreader.NLocalFields()==1) {
            fld = boost::shared_ptr<basisfield>(new splinefield(mapreader.GetLocalAsSingleSplinefield(clp.FullResIntKsp())));
          }
          else if (mapreader.NLocalFields()>1) {
            fld = boost::shared_ptr<basisfield>(new splinefield(mapreader.GetLocalAsSingleSplinefield(clp.FullResIntKsp(),1)));
	  }
	}
      }
      else {
        fld = boost::shared_ptr<basisfield>(new splinefield(clp.RefSize(),clp.RefVxs(),clp.FullResIntKsp()));
        fld->SetToConstant(1.0);  // Make it fields of ones to start out with
      }
    }
    else if (clp.Basis() == DCT) {
      cout << "NYI" << endl;
    }
    mymapper = boost::shared_ptr<IntensityMapper>(new SSDIntensityMapper(fld,clp.IntLambda()));
  }
  else if (clp.IntMapType() == LOCAL_BIAS_WITH_GLOBAL_NON_LINEAR) {
    boost::shared_ptr<BASISFIELD::basisfield>  fld;  // Null-pointer
    if (clp.Basis() == Spline) {
      if (mapreader.HasLocal()) {
        if (mapreader.LocalFieldSize()!=clp.RefSize() || mapreader.LocalFieldVoxelSize()!=clp.RefVxs()) { // This means dimensions are not identical
          if ((std::fabs(mapreader.LocalFieldSize()[0]*mapreader.LocalFieldVoxelSize()[0] - clp.RefSize()[0]*clp.RefVxs()[0]) > 1e-3) || 
              (std::fabs(mapreader.LocalFieldSize()[1]*mapreader.LocalFieldVoxelSize()[1] - clp.RefSize()[1]*clp.RefVxs()[1]) > 1e-3) || 
              (std::fabs(mapreader.LocalFieldSize()[2]*mapreader.LocalFieldVoxelSize()[2] - clp.RefSize()[2]*clp.RefVxs()[2]) > 1e-3)) { // This means FOV is different
            if (mapreader.LocalFieldSize()!=clp.RefSize()) throw fnirt_error("init_intensity_mapper: Intensity Mapper file has different matrix-size from --ref");
            if (mapreader.LocalFieldVoxelSize()!=clp.RefVxs()) throw fnirt_error("init_intensity_mapper: Intensity Mapper file has different voxel-size from --ref");
	  }
	  else { // This means FOV is the same, so we'll assume up/down-sampling
            fld = boost::shared_ptr<BASISFIELD::basisfield>(new BASISFIELD::splinefield(clp.RefSize(),clp.RefVxs(),clp.FullResIntKsp()));
	    std::vector<unsigned int>  ksp(3); // To avoid possible loss of resolution
            ksp[0] = static_cast<unsigned int>((clp.RefVxs()[0]/mapreader.LocalFieldVoxelSize()[0])*clp.FullResIntKsp()[0]);
            ksp[1] = static_cast<unsigned int>((clp.RefVxs()[1]/mapreader.LocalFieldVoxelSize()[1])*clp.FullResIntKsp()[1]);
            ksp[2] = static_cast<unsigned int>((clp.RefVxs()[2]/mapreader.LocalFieldVoxelSize()[2])*clp.FullResIntKsp()[2]);
            if (mapreader.NLocalFields()==1) {
	      if (clp.Verbose()) cout << "Resampling --intin field" << endl;
	      BASISFIELD::splinefield original = mapreader.GetLocalAsSingleSplinefield(ksp);
	      make_field_with_same_fov(original,fld);
	    }
	    else if (mapreader.NLocalFields()>1) {
	      if (clp.Verbose()) cout << "Resampling --intin field" << endl;
              BASISFIELD::splinefield original = mapreader.GetLocalAsSingleSplinefield(ksp,1);
	      make_field_with_same_fov(original,fld);
	    }
	  }
	}
	else { // Dimensions are identical
	  if (mapreader.NLocalFields()==1) {
            fld = boost::shared_ptr<basisfield>(new splinefield(mapreader.GetLocalAsSingleSplinefield(clp.FullResIntKsp())));
          }
          else if (mapreader.NLocalFields()>1) {
            fld = boost::shared_ptr<basisfield>(new splinefield(mapreader.GetLocalAsSingleSplinefield(clp.FullResIntKsp(),1)));
	  }
	}
      }
      else {
        fld = boost::shared_ptr<BASISFIELD::basisfield>(new BASISFIELD::splinefield(clp.RefSize(),clp.RefVxs(),clp.FullResIntKsp()));
        fld->SetToConstant(1.0);  // Make it fields of ones to start out with
      }
    }
    else if (clp.Basis() == DCT) {
      cout << "NYI" << endl;
    }
    vector<double> se(clp.IntMapOrder(),0.0);
    se[1] = 1.0;                // Set linear component to 1
    if (mapreader.HasGlobal()) {
      vector<double> old_global = mapreader.GetGlobal();
      if (old_global.size()==1) se[1] = old_global[0];
      else if (old_global.size()>1) for (unsigned int i=0; i<Min(se.size(),old_global.size()); i++) se[i] = old_global[i];
    }
    mymapper = boost::shared_ptr<IntensityMapper>(new SSDIntensityMapper(se,fld,clp.IntLambda()));
  }
  else if (clp.IntMapType() == LOCAL_NON_LINEAR) {
    std::vector<boost::shared_ptr<basisfield> >  fld(clp.IntMapOrder());
    if (clp.Basis() == Spline) {
      if (mapreader.HasLocal()) {
      }
      else {  // Brand new fields
        for (unsigned int i=0; i<clp.IntMapOrder(); i++) {
          fld[i] = boost::shared_ptr<BASISFIELD::basisfield>(new BASISFIELD::splinefield(clp.RefSize(),clp.RefVxs(),clp.FullResIntKsp()));
          if (i==1) fld[i]->SetToConstant(1.0);  // Make it linear to start out with
        }
      }
    }
    else if (clp.Basis() == DCT) {
      cout << "NYI" << endl;
    }
    mymapper = boost::shared_ptr<IntensityMapper>(new SSDIntensityMapper(fld,clp.IntLambda()));
  }
  return(mymapper);
}

void make_field_with_same_fov(const BASISFIELD::splinefield&                   in,
                              boost::shared_ptr<BASISFIELD::basisfield>        out)
{
  NEWIMAGE::volume<float>  volout(out->FieldSz_x(),out->FieldSz_y(),out->FieldSz_z());
  volout.setdims(out->Vxs_x(),out->Vxs_y(),out->Vxs_z());

  double zf = out->Vxs_z() / in.Vxs_z();
  double yf = out->Vxs_y() / in.Vxs_y();
  double xf = out->Vxs_x() / in.Vxs_x();
  for (unsigned int k=0; k<out->FieldSz_z(); k++) {
    for (unsigned int j=0; j<out->FieldSz_y(); j++) {
      for (unsigned int i=0; i<out->FieldSz_x(); i++) {
        volout(i,j,k) = in.Peek(static_cast<double>(xf*i),static_cast<double>(yf*j),static_cast<double>(zf*k));
      }
    }
  }
  out->Set(volout);

  return;
}

void set_nlpars(NonlinParam&  nlp)
{
  if (nlp.Method() == NL_LM) {
    nlp.SetEquationSolverMaxIter(500);
    nlp.SetEquationSolverTol(1.0e-3);
  }
  else if (nlp.Method() == NL_CG) {
  }
  else if (nlp.Method() == NL_SCG) {
  }
}

// Combine an inclusive/exclusive explicit mask with an implicit mask

boost::shared_ptr<NEWIMAGE::volume<char> > make_mask(const string&                  mfname, 
                                                     MaskType                       mt, 
                                                     const NEWIMAGE::volume<float>& ima, 
                                                     bool                           impf, 
                                                     double                         impv)
{
  boost::shared_ptr<NEWIMAGE::volume<char> >  maskp;
  if (!((mfname.length() && mt!=IgnoreMask) || impf)) return(maskp); // Return null-pointer
  
  if (mfname.length() && mt!=IgnoreMask) { // If there is an explicit mask
    maskp = boost::shared_ptr<NEWIMAGE::volume<char> >(new NEWIMAGE::volume<char>);
    NEWIMAGE::volume<char>&  mask = *maskp;
    read_volume(mask,mfname);
    for (int k=0; k<mask.zsize(); k++) {
      for (int j=0; j<mask.ysize(); j++) {
        for (int i=0; i<mask.xsize(); i++) {
	  if (mask(i,j,k)!=1 && mask(i,j,k)!=0) throw fnirt_error("make_mask: Mask must be binary");
	  if (mt == ExclusiveMask) mask(i,j,k) = (mask(i,j,k)) ? 0 : 1;  // Turn exclusive to inclusive
	}
      }
    }
  }
  else { // No explicit mask then, but we know impf is true
    maskp = boost::shared_ptr<NEWIMAGE::volume<char> >(new NEWIMAGE::volume<char>(ima.xsize(),ima.ysize(),ima.zsize()));
    maskp->setdims(ima.xdim(),ima.ydim(),ima.zdim());
    *maskp = 1;
  }

  if (impf) { // Add implicit mask 
    NEWIMAGE::volume<char>&  mask = *maskp;
    for (int k=0; k<mask.zsize(); k++) {
      for (int j=0; j<mask.ysize(); j++) {
        for (int i=0; i<mask.xsize(); i++) {
          mask(i,j,k) = (fabs(double(ima(i,j,k))-impv) < 1e-16) ? 0 : mask(i,j,k);
	}
      }
    }
  }
  
  return(maskp);  
}

double spmlike_mean(NEWIMAGE::volume<float>&  ima)
{
  double mean = 0.0;
  // First pass to get mean of all voxels
  for (int k=0; k<ima.zsize(); k++) {
    for (int j=0; j<ima.ysize(); j++) {
      for (int i=0; i<ima.xsize(); i++) mean += ima(i,j,k);
    }
  }
  mean /= (ima.xsize()*ima.ysize()*ima.zsize());
  double thres = 0.125*mean;  // 1/8
  mean = 0.0;
  // Now take the mean of all voxels > 1/8 * mean
  int n = 0;
  for (int k=0; k<ima.zsize(); k++) {
    for (int j=0; j<ima.ysize(); j++) {
      for (int i=0; i<ima.xsize(); i++) {
        if (ima(i,j,k) > thres) {mean += ima(i,j,k); n++;}
      }
    }
  }
  mean /= n;
  
  return(mean);
}


////////////////////////////////////////////////////////////////////////////
//
// Check for attempt to register to self
//
////////////////////////////////////////////////////////////////////////////

bool trying_to_register_to_self(const string&                    ref_fname,
				const NEWIMAGE::volume<float>&   ref,
				const string&                    obj_fname,
				const NEWIMAGE::volume<float>&   obj,
                                const NEWMAT::Matrix&            aff)
{
  if (is_identity(aff)) {
    if (ref_fname == obj_fname) return(true);
    if (ref == obj) return(true);
  }
  
  return(false);
}

bool is_identity(const NEWMAT::Matrix&   A,
                 double                  prec)
{
  if (A.Nrows() != A.Ncols()) return(false);
  if ((A - IdentityMatrix(A.Nrows())).MaximumAbsoluteValue() < prec) return(true);
  return(false);
}

////////////////////////////////////////////////////////////////////////////
//
// Write results relevant to self-registration
//
////////////////////////////////////////////////////////////////////////////

void write_self_results(const fnirt_clp&                clp,
                        const NEWIMAGE::volume<float>&  ref)
{
  // Create a cf-object to do the job for us
  std::vector<boost::shared_ptr<basisfield> >   field = init_warpfield(clp);
  boost::shared_ptr<IntensityMapper>            intmap = init_intensity_mapper(clp);
  boost::shared_ptr<SSD_fnirt_CF>   cf = boost::shared_ptr<SSD_fnirt_CF>(new SSD_fnirt_CF(ref,ref,IdentityMatrix(4),field,intmap));

  cf->SaveDefCoefs(clp.CoefFname());                                      // Coefficients
  if (clp.FieldFname().length()) cf->SaveDefFields(clp.FieldFname());     // Field
  if (clp.JacFname().length()) cf->SaveJacobian(clp.JacFname());          // Jacobian
  if (clp.RefOutFname().length()) cf->SaveScaledRef(clp.RefOutFname());   // Intensity modulated ref scan
  if (clp.ObjOutFname().length()) cf->SaveScaledRef(clp.ObjOutFname());   // Warped object image
  // Intensity-mapping
  if (clp.IntensityMappingFname().length()) cf->SaveIntensityMapping(clp.IntensityMappingFname());  
}

////////////////////////////////////////////////////////////////////////////
//
// Try to find existing file matching the name ref_fname
// It is needed because when we read the filenames from 
// the configuration file the shell will not be able to
// do variable substitutions for us.
//
// It will
// 1. Check to see if there is an absolute path,
//    if so it will check for existence of that.
// 2. If no path explicitly given it will look 
//    in the current directory.
// 3. If no path given and not in the current 
//    directory it will look in ${FSLDIR}/data/standard
//
////////////////////////////////////////////////////////////////////////////

string existing_ref_fname(const string& ref_fname)
{
  if (FNIRT::path(ref_fname).length()) {      // If there is an explicit path
    if (NEWIMAGE::fsl_imageexists(ref_fname)) return(string(ref_fname));
    else return(string(""));
  }
  else {        // If no explicit path
    if (NEWIMAGE::fsl_imageexists(ref_fname)) { // Try to open it in current directory
      return(string(ref_fname));
    }
    else {
      const char *fsldir_ptr = getenv("FSLDIR");
      string eref_fname = string(fsldir_ptr) + string("/data/standard/") + ref_fname;
      if (NEWIMAGE::fsl_imageexists(eref_fname)) return(eref_fname);
      else return(string(""));
    }
  }
}
/*
string existing_ref_fname(const string& ref_fname)
{
  string   eref_fname;
  if (FNIRT::path(ref_fname).length()) {      // If there is an explicit path
    try {
      NEWIMAGE::volume<float> vref;
      NEWIMAGE::read_volume_hdr_only(vref,ref_fname); // Throws if file dont exist
      eref_fname = ref_fname;
    }
    catch(...) {
      return(string(""));     // Return empty string
    }
  }
  else {        // If no explicit path
    NEWIMAGE::volume<float> vref;
    try {       // Try to open it in current directory
      NEWIMAGE::read_volume_hdr_only(vref,ref_fname); // Throws if file dont exist
      eref_fname = ref_fname;
    }
    catch(...) { // Didn't exist in current directory, try in ${FSLDIR}/data/standard
      const char *fsldir_ptr = getenv("FSLDIR");
      eref_fname = string(fsldir_ptr) + string("/data/standard/") + ref_fname;
      try {
        cout << "Could not find " << ref_fname << ", now checking " << eref_fname << endl;
        NEWIMAGE::read_volume_hdr_only(vref,eref_fname); // Throws if file dont exist
      }
      catch(...) {
        return(string(""));
      }
    }
  }
  return(eref_fname);
}
*/
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
  if (check_exist(ecfname)) return(ecfname);
  if (!FNIRT::extension(ecfname).length()) {        // If no extension explicitly given
    ecfname += string(".cnf");
    if (check_exist(ecfname)) return(ecfname);
  }
  if (!FNIRT::path(cfname).length()) {              // If no path explicitly given
    const char *fsldir_ptr = getenv("FSLDIR");
    ecfname = string(fsldir_ptr) + string("/etc/flirtsch/") + cfname;
    if (check_exist(ecfname)) return(ecfname);
    else if (!FNIRT::extension(ecfname).length()) { // If no path _and_ no extension given
      ecfname += string(".cnf");
      if (check_exist(ecfname)) return(ecfname);
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

/* Old version (used for first release that had a problem 
bool check_exist(const string& fname)
{
  std::ifstream  ins;

  cout << "Attempting to open file named " << fname << endl;
  ins.open(fname.c_str(),std::ios::in);
  ins.close();
  return((ins.fail()) ? false : true);
}
*/

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


} // End namespace FNIRT
