// Declarations of cost-function classes used by FNIRT
//
// fnirt_costfunctions.h
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2007 University of Oxford 
//
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


#ifndef fnirt_costfunctions_h
#define fnirt_costfunctions_h

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
#include "newmat.h"
#include "miscmaths/bfmatrix.h"
#include "miscmaths/nonlin.h"
#include "basisfield/basisfield.h"
#include "basisfield/splinefield.h"
#include "basisfield/dctfield.h"
#include "intensity_mappers.h"
#include "matching_points.h"

namespace FNIRT {

enum FnirtInterpolationType {LinearInterp, SplineInterp, UnknownInterp};
enum RegularisationType {MembraneEnergy, BendingEnergy, LinearEnergy};

class FnirtException: public std::exception
{
private:
  std::string m_msg;
public:
  FnirtException(const std::string& msg) throw(): m_msg(msg) {}

  virtual const char * what() const throw() {
    return string("Fnirt: msg=" + m_msg).c_str();
  }

  ~FnirtException() throw() {}
};


///////////////////////////////////////////////////////////////////////////////////////////////
//
// fnirt_CF is the base-class for cost functions used by fnirt. It implements
// general functions like initialisation, masking etc, but is still not a 
// proper cost-function and is still pure virtual. It should be sub-classed
// to implement an actual cost-function.
//
// Most objects in the class are represented as smart pointers, which should
// allow both safe and memory efficient use. It allows the user to request and
// get handles to the internal objects without risk of memory leak or of 
// objects suddenly passing out of scope.
//
// Note that there are two constructors. One which takes references to the input
// image volumes and copies these, making it safe albeit a little wasteful. The
// other takes smart pointers and just copies these. That means no copying of
// images take place, but also that there is another set of pointers out there
// that might potentially change our image from under our feet.
//
///////////////////////////////////////////////////////////////////////////////////////////////

class fnirt_CF: public MISCMATHS::NonlinCF
{
public:
  // Constructors and destructors

  // Safe (and slightly wasteful) constructor
  fnirt_CF(const NEWIMAGE::volume<float>&                                 pvref, 
           const NEWIMAGE::volume<float>&                                 pvobj,
           const NEWMAT::Matrix&                                          aff,
           std::vector<boost::shared_ptr<BASISFIELD::basisfield> >&       pdf_field,
           boost::shared_ptr<IntensityMapper>                             im = boost::shared_ptr<IntensityMapper>());

  // A little unsafe, but faster
  fnirt_CF(const boost::shared_ptr<NEWIMAGE::volume<float> >              pvref, 
           const boost::shared_ptr<NEWIMAGE::volume<float> >              pvobj,
           const NEWMAT::Matrix&                                          paff,
           std::vector<boost::shared_ptr<BASISFIELD::basisfield> >&       pdf_field,
           boost::shared_ptr<IntensityMapper>                             im = boost::shared_ptr<IntensityMapper>());

  virtual ~fnirt_CF() {}

  // Find out # of parameters
  virtual int NPar() const = 0;
  // Get current set of parameters
  virtual NEWMAT::ReturnMatrix Par() const = 0;
  // Set model for interpolation
  virtual void SetInterpolationModel(FnirtInterpolationType it) {interp = it; vobj->setinterpolationmethod(translate_interp_type(it)); svobj->setinterpolationmethod(translate_interp_type(it)); }
  // Set model to use for regularisation of field
  virtual void SetRegularisationModel(RegularisationType  pregmod) {regmod = pregmod;} 
  // Set relative weight of (membrane energy) regularisation
  virtual void SetLambda(double plambda) {lambda = plambda; if (Verbose()) cout << "New Lambda: " << lambda << endl;}
  // Specify that lambda should be weighted by latest ssd
  virtual void WeightLambdaBySSD(bool  wl=true) {ssd_lambda=wl;}
  // Specify that derivatives of reference image should be used (fast approximation)
  virtual void UseRefDerivs() {use_ref_derivs = true;}
  // Specify that derivatives of object image should be used ("exact", but slower)
  virtual void UseObjDerivs() {use_ref_derivs = false;}
  // Specify how chatty we should be while running
  virtual void SetVerbose(bool status=true) {verbose=status;}
  // Save debug information or not
  virtual void SetDebug(unsigned int level=1) {debug=level;}
  // Set what level we are at (used for debug info)
  virtual void SetLevel(unsigned int plevel=0) {level=plevel; iter=0; attempt=0;}
  // Set precision for representation of Hessian
  virtual void SetHessianPrecision(MISCMATHS::BFMatrixPrecisionType prec) {hess_prec=prec;}
  // Set list of matching points/landmarks to be included in matching.
  virtual void SetMatchingPoints(const MatchingPoints& pmpl) {mpl = boost::shared_ptr<MatchingPoints>(new MatchingPoints(pmpl));}
  virtual void SetMatchingPointsLambda(double pl) {mpl_lambda = pl;}

  // Routines for re-scaling the intensity of ref or obj images
  virtual void IntensityScaleObj(double sfac);  

  // Routines for applying Gaussian smoothing to ref or obj images
  virtual void SmoothObj(double  fwhm);
  virtual void SmoothRef(double  fwhm);
  
  // Subsampling reference volume
  virtual void SubsampleRef(unsigned int ss) {std::vector<unsigned int> ssv(3,ss); SubsampleRef(ssv);}
  virtual void SubsampleRef(const std::vector<unsigned int>& ss);

  // Setting optional masks in the space of the reference and/or the object image
  // Safe versions
  virtual void SetRefMask(const NEWIMAGE::volume<char>& prefmask);
  virtual void SetObjMask(const NEWIMAGE::volume<char>& pobjmask);

  // Potentially dodgy, but faster versions
  virtual void SetRefMask(const boost::shared_ptr<NEWIMAGE::volume<char> > prefmask);
  virtual void SetObjMask(const boost::shared_ptr<NEWIMAGE::volume<char> > pobjmask);

  // Routines to save fields as "fields" or as coefficient volumes
  virtual void SaveDefFields(std::string fname) const;
  virtual void SaveDefCoefs(const std::string &fname) const;

  // Save Jacobian of displacement field. For diagnostic purposes
  virtual void SaveJacobian(const std::string &fname) const;

  // Save a copy of the resampled object image.
  virtual void SaveRobj(const std::string &fname) const;

  // Save a copy of "working copy" of reference volume. For debugging mainly
  virtual void SaveRef(const std::string &fname) const;

  // Save a copy of "working copy" of reference volume after intensity scaling. For debugging mainly
  virtual void SaveScaledRef(const std::string &fname) const;

  // Save a copy of "working copy" of mask volume. For debugging mainly
  virtual void SaveMask(const std::string &fname) const;

  // Save a copy of "working copy" of ScaledRef-Robj masked by mask
  virtual void SaveMaskedDiff(const std::string &fname) const;

  // Save a copy of "working copy" of various components of the mask volume. For debugging mainly
  virtual void SaveDataMask(const std::string &fname) const;
  virtual void SaveRobjMask(const std::string &fname) const;
  virtual void SaveObjMask(const std::string &fname) const;

  // Set intensity mapping to be fixed, i.e. not part of the 
  // parameters we are trying to estimate.
  virtual void SetIntensityMappingFixed(bool status=true) {imap->SetFixed(status);}

  // Save local and global parameters/fields for intensity mapping
  virtual void SaveLocalIntensityMapping(std::string fname = string("fnirt_local_intensity_fields")) const;
  virtual void SaveGlobalIntensityMapping(std::string fname = string("fnirt_global_intensity_parameters.txt")) const;
  virtual void SaveIntensityMapping(std::string fname = string("fnirt_intensity_parameters")) const
  {
    SaveGlobalIntensityMapping(fname+string(".txt"));
    SaveLocalIntensityMapping(fname);
  }

  // Access displacement field as volume
  virtual const NEWIMAGE::volume<float>& DefFieldAsVol(int index) const
  {
    if (index<0 || index >2) {throw FnirtException("fnirt_CF::DefFieldAsVol: index must be 0, 1 or 2");}
    return((*df)[index]);
  }

  // Access displacement field as field class
  virtual const BASISFIELD::basisfield& DefField(int index) const
  {
    if (index<0 || index >2) {throw FnirtException("fnirt_CF::DefField: index must be 0, 1 or 2");}
    return(*(df_field[index]));
  }

  // Get range of Jacobian values for non-linear part
  std::pair<double,double> JacobianRange() const;

  // Set range of Jacobians to within specified bounds
  void ForceJacobianRange(double minj, double maxj) const;

  // Get some general info
  virtual unsigned int OrigRefSz_x() const {return(static_cast<unsigned int>(vref->xsize()));}
  virtual unsigned int OrigRefSz_y() const {return(static_cast<unsigned int>(vref->ysize()));}
  virtual unsigned int OrigRefSz_z() const {return(static_cast<unsigned int>(vref->zsize()));}
  virtual unsigned int OrigRefSz() const {return(OrigRefSz_x()*OrigRefSz_y()*OrigRefSz_z());}
  virtual double OrigRefVxs_x() const {return(vref->xdim());}
  virtual double OrigRefVxs_y() const {return(vref->ydim());}
  virtual double OrigRefVxs_z() const {return(vref->zdim());}
  virtual unsigned int RefSz_x() const {return(static_cast<unsigned int>(ssvref->xsize()));}
  virtual unsigned int RefSz_y() const {return(static_cast<unsigned int>(ssvref->ysize()));}
  virtual unsigned int RefSz_z() const {return(static_cast<unsigned int>(ssvref->zsize()));}
  virtual unsigned int RefSz() const {return(RefSz_x()*RefSz_y()*RefSz_z());}
  virtual double RefVxs_x() const {return(ssvref->xdim());}
  virtual double RefVxs_y() const {return(ssvref->ydim());}
  virtual double RefVxs_z() const {return(ssvref->zdim());}
  virtual unsigned int ObjSz_x() const {return(static_cast<unsigned int>(vobj->xsize()));}
  virtual unsigned int ObjSz_y() const {return(static_cast<unsigned int>(vobj->ysize()));}
  virtual unsigned int ObjSz_z() const {return(static_cast<unsigned int>(vobj->zsize()));}
  virtual unsigned int ObjSz() const {return(ObjSz_x()*ObjSz_y()*ObjSz_z());}
  virtual double ObjVxs_x() const {return(vobj->xdim());}
  virtual double ObjVxs_y() const {return(vobj->ydim());}
  virtual double ObjVxs_z() const {return(vobj->zdim());}
  virtual unsigned int DefCoefSz_x() const {return(df_field[0]->CoefSz_x());}
  virtual unsigned int DefCoefSz_y() const {return(df_field[0]->CoefSz_y());}
  virtual unsigned int DefCoefSz_z() const {return(df_field[0]->CoefSz_z());}
  virtual unsigned int DefCoefSz() const {return(df_field[0]->CoefSz());}
  virtual unsigned int ScaleCoefSz() const {return(imap->NPar());}

protected:
  // Utility functions for internal use in this and derived classes.

  // Get read access to initial affine guess.
  virtual const NEWMAT::Matrix& AffineMat() const {return(aff); }

  // Find out what interpolation model to use
  virtual FnirtInterpolationType InterpolationType() const {return(interp);}

  // Find out what regularisation model to use
  virtual RegularisationType RegularisationModel() const {return(regmod);}
  
  // Get weight of regularisation
  virtual double Lambda() const {if (ssd_lambda) return(latest_ssd*lambda); else return(lambda);}

  // Get relative weight of matching points/landmarks
  virtual double MatchingPointsLambda() const {return(mpl_lambda);}

  // See if we are using derivatives from reference image
  virtual bool UsingRefDeriv() const {return(use_ref_derivs);}

  // Set/Get parameters describing the field
  virtual void SetDefFieldParams(const NEWMAT::ColumnVector& p) const;
  virtual NEWMAT::ColumnVector GetDefFieldParams(int  findx = -1) const; 
 
  // Set deformation field directly (i.e. not via params)
  virtual void SetDefField(const NEWIMAGE::volume4D<float>&   vfield) const;

  // Set/Get parameters describing intensity scaling
  virtual void SetScaleParams(const NEWMAT::ColumnVector& p) const;
  virtual NEWMAT::ColumnVector GetScaleParams() const; 

  // Calculate resulting mask
  virtual const NEWIMAGE::volume<char>& Mask() const;

  // Get resampled object image.
  virtual const NEWIMAGE::volume<float>& Robj() const;

  // Get resampled derivatives of object image
  virtual const NEWIMAGE::volume4D<float>& RobjDeriv() const;

  // Get reference image
  virtual const NEWIMAGE::volume<float>& Ref() const;

  // Get intensity-scaled reference image
  virtual void ScaledRef(NEWIMAGE::volume<float>&  sref) const;

  // Get derivatives of reference image in original space
  virtual const NEWIMAGE::volume4D<float>& RefDeriv() const;

  // Get derivatives according to setting of use_ref_derivs
  virtual const NEWIMAGE::volume4D<float>& Deriv() const;

  // Set last observed value for ssd
  virtual void SetLatestSSD(double ssd) const {latest_ssd = ssd;}

  // Allow read-access to intensity mapper
  virtual const boost::shared_ptr<IntensityMapper>& IMap() const {return(imap);}

  // Allow read-access to point-list mapper
  virtual const boost::shared_ptr<MatchingPoints>& MPL() const {return(mpl);}

  // Set new value for intensity mapper
  virtual void SetIMap(const boost::shared_ptr<IntensityMapper>&  new_imap) {imap = new_imap;}

  // Find out how chatty to be
  virtual bool Verbose() const {return(verbose);}

  // "Local" change of verbose flag
  virtual void LocalSetVerbose(bool status=true) const {verbose=status;}
  
  // Find out what precision to use for Hessian
  virtual BFMatrixPrecisionType HessianPrecision() const {return(hess_prec);}

  // Find out if debug info should be saved
  virtual unsigned int Debug() const {return(debug);}

  // Some routines used for debug info
  virtual void SetIter(unsigned int piter=0) const {iter=piter; attempt=0;}
  virtual void SetAttempt(unsigned int pattempt) const {attempt=pattempt;}
  virtual unsigned int Level() const {return(level);}
  virtual unsigned int Iter() const {return(iter);}
  virtual unsigned int Attempt() const {return(attempt);}
  virtual std::string DebugString() const {char c_str[256]; sprintf(c_str,"level%02d_iter%02d_attempt%02d",Level(),Iter(),Attempt()); return(std::string(c_str));}

private:
  // Hide auto generated functions
  fnirt_CF();
  fnirt_CF(const fnirt_CF& inf);
  virtual fnirt_CF& operator=(const fnirt_CF& inf) {throw FnirtException("fnirt_CF:: Assignment explicitly disallowed"); return(*this);}

  const boost::shared_ptr<NEWIMAGE::volume<float> >        vref;               // Reference volume
  boost::shared_ptr<NEWIMAGE::volume<float> >              svref;              // Smoothed reference volume
  double                                                   ref_fwhm;           // FWHM of smoothed reference volume
  boost::shared_ptr<NEWIMAGE::volume<float> >              ssvref;             // Subsampled smoothed reference volume
  std::vector<unsigned int>                                subsamp;            // Subsampling of subsampled reference volume

  const boost::shared_ptr<NEWIMAGE::volume<float> >        vobj;       // Object volume, i.e. volume to be warped to vref
  boost::shared_ptr<NEWIMAGE::volume<float> >              svobj;      // Smoothed object volume
  double                                                   obj_fwhm;   // FWHM of smoothed object volume

  const NEWMAT::Matrix                                     aff;        // Affine transformation matrix

  // The field is represented as a 3-element vector of fields, one per displacement direction  

  mutable std::vector<boost::shared_ptr<BASISFIELD::basisfield> >  df_field;   // Displacement field implementing the warp

  // The IntensityMapper will perform any calculations associated with mapping intensitites of ref to obj

  mutable boost::shared_ptr<IntensityMapper>               imap;

  // The matching-points list will calculate contribution from (optional) anatomical landmark information

  boost::shared_ptr<MatchingPoints>                        mpl;

  mutable double                                           latest_ssd;         // For use when weighting lambda by ssd
  RegularisationType                                       regmod;             // Model for regularisation
  double                                                   lambda;             // Relative weight of regularisation
  double                                                   mpl_lambda;         // Relative weight of landmark-list
  bool                                                     ssd_lambda;         // Set if lambda is to be multiplied by latest ssd
  bool                                                     use_ref_derivs;     // Use derivatives of reference image
  BFMatrixPrecisionType                                    hess_prec;          // Can be float or double
  FnirtInterpolationType                                   interp;             // Interpolation model trilinear/spline
  mutable bool                                             verbose;            // Print diagnostic information
  unsigned int                                             debug;              // Level of debug info to save

  // State variables used when outputting debug information

  mutable unsigned int                                     level;              // Level (of subsampling, fwhm etc)
  mutable unsigned int                                     iter;               // Iteration, decided by # of times we calculate hessian
  mutable unsigned int                                     attempt;            // Attempt (within an iteration)

  // Optional masks in ref and/or object space
  
  boost::shared_ptr<NEWIMAGE::volume<char> >    refmask;                       // Mask in space of reference volume
  boost::shared_ptr<NEWIMAGE::volume<char> >    ssrefmask;                     // Mask in sub-sampled reference volume
  boost::shared_ptr<NEWIMAGE::volume<char> >    objmask;                       // Mask in space of object volume
  
  // Scratch pads

  mutable boost::shared_ptr<NEWIMAGE::volume<float> >   robj;                // Resampled version of object volume
  mutable boost::shared_ptr<NEWIMAGE::volume<char> >    datamask;            // One where robj reflects "true" interpolated intensity
  mutable bool                                          robj_updated;        // True if robj accurately reflects df
  mutable boost::shared_ptr<NEWIMAGE::volume4D<float> > df;                  // Volume representation of displacement field
  mutable boost::shared_ptr<NEWIMAGE::volume<char> >    robjmask;            // Resampled version of objmask
  mutable bool                                          robjmask_updated;    // True if robjmask accurately reflects df
  mutable boost::shared_ptr<NEWIMAGE::volume<char> >    totmask;             // Total mask, only used if refmask and/or objmask is set.
  mutable bool                                          totmask_updated;   
  mutable boost::shared_ptr<NEWIMAGE::volume4D<float> > robj_deriv;          // Derivative of object image at resampled points.
  mutable bool                                          robj_deriv_updated;  // True if derivatives accurately reflects current df
  mutable boost::shared_ptr<NEWIMAGE::volume4D<float> > ref_deriv;           // Derivative of reference image in native space

  // Utility functions for use only in this class

  NEWIMAGE::interpolation translate_interp_type(FnirtInterpolationType it) const { if (it==LinearInterp) return(NEWIMAGE::trilinear); else return(NEWIMAGE::spline); }
  FnirtInterpolationType translate_interp_type(NEWIMAGE::interpolation it) const { if (it==NEWIMAGE::trilinear) return(LinearInterp); else if (it==NEWIMAGE::spline) return(SplineInterp); else return(UnknownInterp); }

  // Convert subsampling-factor to new matrix size
  std::vector<unsigned int> ss2ms(const std::vector<unsigned int>& ss) const
  {
    std::vector<unsigned int> ms(3);
    ms[0]=OrigRefSz_x(); ms[1]=OrigRefSz_y(); ms[2]=OrigRefSz_z();
    return(df_field[0]->SubsampledMatrixSize(ss,ms));
  }
  // Convert subsampling-factor to new voxel size
  std::vector<double> ss2vxs(const std::vector<unsigned int> ss) const
  {
    std::vector<unsigned int> ms(3);
    ms[0]=OrigRefSz_x(); ms[1]=OrigRefSz_y(); ms[2]=OrigRefSz_z();
    std::vector<double> vxs(3);
    vxs[0]=OrigRefVxs_x(); vxs[1]=OrigRefVxs_y(); vxs[2]=OrigRefVxs_z();
    return(df_field[0]->SubsampledVoxelSize(ss,vxs,ms));
  }
	   
  void subsample_ref(const std::vector<unsigned int>  nms,
                     const std::vector<double>        nvxs);

  void subsample_refmask(const std::vector<unsigned int>  nms,
                         const std::vector<double>        nvxs);

  boost::shared_ptr<NEWIMAGE::volume<float> > masked_smoothing(const NEWIMAGE::volume<float>&              vol,
                                                               double                                      fwhm,
                                                               boost::shared_ptr<NEWIMAGE::volume<char> >  mask);

  unsigned int good_fft_size(unsigned int isz) const;

  // Tasks common to both constructors
  void common_construction(std::vector<boost::shared_ptr<BASISFIELD::basisfield> >& pdf_field);

  // Update robj using current field
  void update_robj() const;

  // Update robjmask using current field
  void update_robjmask() const;

};


///////////////////////////////////////////////////////////////////////////////////////////////
//
// Here comes the first proper cost function. It inherits all the functionality from
// fnirt_CF regarding maskin, initialisation etc. In addition it implements the sum-of-squared
// differences cost-function, its gradient and its hessian. The parameters are the parameters
// for the x-, y- and z-displacement fields, in that order, followed by intensity scaling
// parameters that are also estimated.
//
///////////////////////////////////////////////////////////////////////////////////////////////

class SSD_fnirt_CF: public fnirt_CF
{
public:
  // Constructors and destructors

  // Safe and wasteful
  SSD_fnirt_CF(const NEWIMAGE::volume<float>&                            pvref, 
               const NEWIMAGE::volume<float>&                            pvobj,
               const NEWMAT::Matrix                                      aff,
               std::vector<boost::shared_ptr<BASISFIELD::basisfield> >&  pdf_field,
               boost::shared_ptr<IntensityMapper>                        im = boost::shared_ptr<IntensityMapper>())
    : fnirt_CF(pvref,pvobj,aff,pdf_field,im)
  {
    last_hess = boost::shared_ptr<BFMatrix>();   // NULL
    if (!im) {                                   // Linear scaling default for ssd
      boost::shared_ptr<IntensityMapper>  def_imap(new SSDIntensityMapper(1.0));
      SetIMap(def_imap);
    }
  }

  // Unsafe (if you misbehave), but faster.
  SSD_fnirt_CF(const boost::shared_ptr<NEWIMAGE::volume<float> >         pvref, 
               const boost::shared_ptr<NEWIMAGE::volume<float> >         pvobj,
               const NEWMAT::Matrix                                      aff,
               std::vector<boost::shared_ptr<BASISFIELD::basisfield> >&  pdf_field,
               boost::shared_ptr<IntensityMapper>                        im = boost::shared_ptr<IntensityMapper>())
    : fnirt_CF(pvref,pvobj,aff,pdf_field,im)
  {
    last_hess = boost::shared_ptr<BFMatrix>();   // NULL
    if (!im) {                                   // Linear scaling default for ssd
      boost::shared_ptr<IntensityMapper>  def_imap(new SSDIntensityMapper(1.0));
      SetIMap(def_imap);
    }
  }
  
  virtual ~SSD_fnirt_CF() {}

  // Find out # of parameters
  virtual int NPar() const {return(3*DefCoefSz()+ScaleCoefSz());}
  // Read access to current set of parameters
  virtual NEWMAT::ReturnMatrix Par() const; 

  // Cost function, gradient and hessian
  // For SSD_fnirt_CF the first parameter is a scale-factor between the reference
  // and the object and the remaining parameters describe the x-, y- and z-displacement
  // fields in that order.
  virtual double cf(const NEWMAT::ColumnVector& p) const;
  virtual NEWMAT::ReturnMatrix grad(const NEWMAT::ColumnVector& p) const;
  virtual boost::shared_ptr<MISCMATHS::BFMatrix> hess(const NEWMAT::ColumnVector& p,
                                                      boost::shared_ptr<MISCMATHS::BFMatrix>  iptr=boost::shared_ptr<MISCMATHS::BFMatrix>()) const;
  
private:
  // I can't see a need for default constructor, assignment or copy construction, so I'll hide them.
  SSD_fnirt_CF();
  SSD_fnirt_CF(const SSD_fnirt_CF& inf);
  virtual SSD_fnirt_CF& operator=(const SSD_fnirt_CF& inf) {throw FnirtException("SSD_fnirt_CF:: Assignment explicitly disallowed"); return(*this);}
  
  // Internal copy of ptr-> last calculated Hessian. Useful when basing hessian on reference image derivatives
  mutable boost::shared_ptr<MISCMATHS::BFMatrix>                   last_hess;
  // Internal copy of ptr-> regularisation part of Hessian. Saves us a little time
  mutable boost::shared_ptr<MISCMATHS::BFMatrix>                   reg;
};

  
} // End namespace FNIRT

#endif // end #ifndef fnirt_costfunctions_h
