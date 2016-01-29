// Declarations of intensity-mapper classes used by FNIRT
//
// intensity_mappers.h
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


#ifndef intensity_mappers_h
#define intensity_mappers_h

#include <string>
#include <vector>
#include "newmat.h"
#include "miscmaths/bfmatrix.h"
#include "basisfield/basisfield.h"
#include "basisfield/splinefield.h"
#include "basisfield/dctfield.h"

namespace FNIRT {

class IntensityMapperException: public std::exception
{
private:
  std::string m_msg;
public:
  IntensityMapperException(const std::string& msg) throw(): m_msg(msg) {}

  virtual const char * what() const throw() {
    return string("IntensityMapper:: msg=" + m_msg).c_str();
  }

  ~IntensityMapperException() throw() {}
};

///////////////////////////////////////////////////////////////////////////////////////////////
//
// The IntensityMapper is a class that takes care of the intensity-mapping between a
// reference-volume and an object-colume within a firt_CF class. The IntensityMapper
// may map the intensity in the following ways depending on the respective sizes of
// _sfac and _sfld.
//
// 1. _sfac.size() == 1 && _sfld.size() == 0 :
// In this case there is no intensity scaling and it is assumed that
// obj(x,y,z) = ref(x,y,z) + err 
// Probably never useful.
//
// 2. _sfac.size() == 1 && _sfld.size() == 0 :
// This is the "normal" case where we have a single scale-factor between the ref and
// obj images such that 
// obj(x,y,z) = _sfac[0] * ref(x,y,z) + err
//
// 3. _sfac.size() > 1 && _sfld.size() == 0 :
// This is used to accomplish a "global" non-linear mapping between the intensities
// in the ref and obj such that 
// obj(x,y,z) = _sfac[0] + _sfac[1]*ref(x,y,z) + _sfac[2]*ref(x,y,z).^2 + ... + err
// This gives a "correlation-ratio like" cost-function
//
// 4. _sfac.size() == 0 && _sfld.size() == 1 :
// This is used to accomplish a "local" intensity modulation to account for e.g.
// different bias-fields between the ref and end the obj volume such that
// obj(x,y,z) = _sfld[0].Peek(x,y,z) * ref(x,y,z) + err
//
// 5. _sfac.size() > 1 && _sfld.size() == 1 :
// This combines a global non-linear intensity mapping between ref and obj with a
// a local bias-field modulation such that
// obj(x,y,z) = _sfld[0].Peek(x,y,z) * (_sfac[0] + _sfac[1]*ref(x,y,z) + ...) + err
//
// 6. _sfac.size() == 0 && _sfld.size() > 1
// This is the most general case and and allows for a local non-linear intensity
// mapping between ref and obj such that
// obj(x,y,z) = _sfld[0].Peek(x,y,z) + _sfld[1].Peek(x,y,z)*ref(x,y,z) + ... + err
//
///////////////////////////////////////////////////////////////////////////////////////////////

class IntensityMapper
{
public:
  // Constructors and destructors

  // Constructor for no scaling
  IntensityMapper() : _sfac(), _presc(), _sfld(), _lambda(0.0), _mt(NO_SCALING), _fixed(false) {}

  // Constructor for global linear mapping
  IntensityMapper(double sfac) : _sfac(1,sfac), _presc(1,1.0), _sfld(), _lambda(0.0), _mt(GLOBAL_LINEAR), _fixed(false) {}

  // Constructor for global non-linear mapping  
  IntensityMapper(const std::vector<double>&  sfac) : _sfac(sfac), _presc(_sfac.size(),1.0), _sfld(), _lambda(0.0), _mt((_sfac.size() > 1) ? GLOBAL_NON_LINEAR : GLOBAL_LINEAR), _fixed(false)
  {
    if (!_sfac.size()) throw IntensityMapperException("IntensityMapper::IntensityMapper: Cannot pass zero length vector");
  }

  // Constructor for local linear mapping
  IntensityMapper(const boost::shared_ptr<BASISFIELD::basisfield>&  sfld,
                  double                                            lambda=10000.0) : _sfac(), _presc(), _sfld(1,sfld), _lambda(lambda), _mt(LOCAL_LINEAR), _fixed(false) {}

  // Constructor for global non-linear mapping with local bias-field
  IntensityMapper(const std::vector<double>&                        sfac,
                  const boost::shared_ptr<BASISFIELD::basisfield>&  sfld,
                  double                                            lambda=10000.0) : _sfac(sfac), _presc(_sfac.size(),1.0), _sfld(1,sfld), _lambda(lambda), _mt(LOCAL_BIAS_WITH_GLOBAL_NON_LINEAR), _fixed(false)
  {
    if (_sfac.size() < 2) throw IntensityMapperException("IntensityMapper::IntensityMapper: Global polynomial must have order > 1 when combined with bias field");
  }

  // Constructor for local non-linear mapping
  IntensityMapper(const std::vector<boost::shared_ptr<BASISFIELD::basisfield> >&   sfld,
                  double                                                           lambda=10000.0) : _sfac(), _presc(), _sfld(sfld), _lambda(lambda), _mt(LOCAL_NON_LINEAR), _fixed(false) {}

  virtual ~IntensityMapper() {} // Thanks to boost

  // General utility functions

  // Number of "components" the the mapping consists of 
  unsigned int NComp() const {return(_sfac.size() + _sfld.size());}

  // Number of parameters needed to completely determine the mapping
  unsigned int NPar() const {
    if (_fixed) return(0);   // No parameters, just constants
    else return(_sfac.size() + ((_sfld.size()) ? _sfld.size() * _sfld[0]->CoefSz() : 0));
  }

  // Returns the parameters of the mapping
  NEWMAT::ReturnMatrix GetPar() const;

  // Returns status of fixed-flag
  bool Fixed() const {return(_fixed);}

  // Sets the status of the fixed-flag
  void SetFixed(bool status=true) {_fixed = status;}

  // Returns the bending-energy weighting
  double Lambda() const {return(_lambda);}

  // Sets the bending-energy weighting
  void SetLambda(double lambda) {_lambda = lambda;}

  // Displays debug info
  void DebugPrint() const;

  // Sets new parameters for the mapping
  void SetPar(const NEWMAT::ColumnVector&  par);

  // SCales the input volume
  void ScaleMe(const NEWIMAGE::volume<float>&  ivol,
               NEWIMAGE::volume<float>&        ovol) const;

  // Returns the scaled membrane energy of any local mapping
  double CFContribution() const;

  // Saves the global intensity mapping vector
  // in NEWMAT ColumnVector style.
  void SaveGlobalMapping(const std::string& fname) const;

  // Saves the local intensity mapping in
  // image format.
  void SaveLocalMapping(const std::string& fname) const;

  // Updates the fields (if any) in accordance with
  // the new parameters.
  void NewSubSampling(const std::vector<unsigned int>&    ms,       // New matrix size
                      const std::vector<double>&          vxs,      // New voxel-size
                      const std::vector<unsigned int>&    new_ss,   // New level of sub-sampling
                      const std::vector<unsigned int>&    old_ss);  // Old level of sub-sampling

  // Calculates the gradient of the cost-function with 
  // respect to the parameters determining the scaling.
  virtual NEWMAT::ReturnMatrix Gradient(const NEWIMAGE::volume<float>&   ref,
                                        const NEWIMAGE::volume<float>&   diff,
                                        const NEWIMAGE::volume<char>     *mask);

  // Calculates the Hessian of the cost-function with 
  // respect to the parameters determining the scaling.
  virtual boost::shared_ptr<MISCMATHS::BFMatrix> Hessian(const NEWIMAGE::volume<float>&    ref,
                                                         const NEWIMAGE::volume<char>      *mask,
                                                         MISCMATHS::BFMatrixPrecisionType  prec=BFMatrixDoublePrecision) const;

  // Calculates the cross-terms (off-diagonal) of the Hessian
  // of the cost-function. I.e. the mixed partial derivatives
  // where the differentiation is w.r.t. one displacement parameter
  // and one scaling parameter.
  virtual boost::shared_ptr<MISCMATHS::BFMatrix> CrossHessian(const BASISFIELD::basisfield&        dfield,
                                                              const NEWIMAGE::volume<float>&       dima,
                                                              const NEWIMAGE::volume<float>&       ref,
                                                              const NEWIMAGE::volume<char>         *mask,
                                                              MISCMATHS::BFMatrixPrecisionType     prec) const;
  
protected:
  enum MappingType {NO_SCALING, GLOBAL_LINEAR, GLOBAL_NON_LINEAR,
                    LOCAL_LINEAR, LOCAL_BIAS_WITH_GLOBAL_NON_LINEAR,
                    LOCAL_NON_LINEAR};
  std::vector<double>                                       _sfac;   // Vector of coefficients for global mapping
  std::vector<double>                                       _presc;  // "pre-scaling" to improve condition # of Hessian
                                                                     // _presc is _only_ used by global_non_linear model
  std::vector<boost::shared_ptr<BASISFIELD::basisfield> >   _sfld;   // Vector of fields for local mapping
  double                                                    _lambda; // Weighting of membrane-energy of field/fields
  const MappingType                                         _mt;     // Flag to determine the current mapping type
  bool                                                      _fixed;  // True if mapping should be kept fixed
};

///////////////////////////////////////////////////////////////////////////////////////////////
//
// The SSDIntensityMapper is a sub-class of IntensityMapper. It provides instances of the
// gradient and hessian methods that are relevent/suitable for sum-of-squared differences
// type cost-functions. 
//
///////////////////////////////////////////////////////////////////////////////////////////////

class SSDIntensityMapper : public IntensityMapper
{
public:
  // Constructors and destructor

  // Constructor for no scaling
  SSDIntensityMapper() : IntensityMapper() {}
  // Constructor for global linear mapping
  SSDIntensityMapper(double sfac) : IntensityMapper(sfac) {}
  // Constructor for global non-linear mapping  
  SSDIntensityMapper(const std::vector<double>&  sfac) : IntensityMapper(sfac) {} 
  // Constructor for local linear mapping
  SSDIntensityMapper(boost::shared_ptr<BASISFIELD::basisfield>&  sfld, 
                     double                                      lambda=1000.0) : IntensityMapper(sfld,lambda) {}
  // Constructor for global non-linear mapping with local bias-field
  SSDIntensityMapper(const std::vector<double>&                        sfac,
                     const boost::shared_ptr<BASISFIELD::basisfield>&  sfld,
                     double                                            lambda=1000.0) : IntensityMapper(sfac,sfld,lambda) {}
  // Constructor for local non-linear mapping
  SSDIntensityMapper(const std::vector<boost::shared_ptr<BASISFIELD::basisfield> >&   sfld,
                     double                                                           lambda=1000.0) : IntensityMapper(sfld,lambda) {}

  ~SSDIntensityMapper() {} // Thanks to boost

  // Calculates the gradient of the cost-function with 
  // respect to the parameters determining the scaling.
  virtual NEWMAT::ReturnMatrix Gradient(const NEWIMAGE::volume<float>&   ref,
                                        const NEWIMAGE::volume<float>&   diff,
                                        const NEWIMAGE::volume<char>     *mask);

  // Calculates the Hessian of the cost-function with 
  // respect to the parameters determining the scaling.
  virtual boost::shared_ptr<MISCMATHS::BFMatrix> Hessian(const NEWIMAGE::volume<float>&    ref,
                                                         const NEWIMAGE::volume<char>      *mask,
                                                         MISCMATHS::BFMatrixPrecisionType  prec) const;

  // Calculates the cross-terms (off-diagonal) of the Hessian
  // of the cost-function. I.e. the mixed partial derivatives
  // where the differentiation is w.r.t. one displacement parameter
  // and one scaling parameter.
  virtual boost::shared_ptr<MISCMATHS::BFMatrix> CrossHessian(const BASISFIELD::basisfield&        dfield,
                                                              const NEWIMAGE::volume<float>&       dima,
                                                              const NEWIMAGE::volume<float>&       ref,
                                                              const NEWIMAGE::volume<char>         *mask,
                                                              MISCMATHS::BFMatrixPrecisionType     prec) const;

private:
  // Utility function that returns mean of a vector, possibly ignoring zeros.
  double vector_mean(const NEWMAT::ColumnVector&  vec, bool exclude_zeros=true) const;

  // Utility function that calculates some sums of elementwise products
  // of images. If you consider vol1, vol2 and vol3 below as column vectors
  // then it returns ((v1.^n).*v2)'*v3
  double funny_dotproduct(const NEWIMAGE::volume<float>&   vol1,
                          unsigned int                     n,
                          const NEWIMAGE::volume<float>&   vol2,
                          const NEWIMAGE::volume<float>&   vol3,
                          const NEWIMAGE::volume<char>&    mask) const;
  double funny_dotproduct(const NEWIMAGE::volume<float>&   vol1,
                          unsigned int                     n,
                          const NEWIMAGE::volume<float>&   vol2,
                          const NEWIMAGE::volume<float>&   vol3) const;
  double funny_dotproduct(const NEWIMAGE::volume<float>&   vol1,
                          unsigned int                     n,
                          const NEWIMAGE::volume<float>&   vol2,
                          const NEWIMAGE::volume<float>&   vol3,
                          const NEWIMAGE::volume<char>     *mask) const;


};

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Class used for reading intensity mapping info from file.
//
///////////////////////////////////////////////////////////////////////////////////////////////

class IntensityMapperReaderException: public std::exception
{
private:
  std::string m_msg;
public:
  IntensityMapperReaderException(const std::string& msg) throw(): m_msg(msg) {}

  virtual const char * what() const throw() {
    return string("IntensityMapperReader:: msg=" + m_msg).c_str();
  }

  ~IntensityMapperReaderException() throw() {}
};

class IntensityMapperReader
{
public:
  IntensityMapperReader() : _has_global(false), _has_local(false) {}
  IntensityMapperReader(const std::string& fname) : _has_global(false), _has_local(false) {common_read(fname);}
  void Read(const std::string& fname) {common_read(fname);}
  bool HasGlobal() const {return(_has_global);}
  bool HasLocal() const {return(_has_local);}
  std::vector<unsigned int> LocalFieldSize() const;
  std::vector<double> LocalFieldVoxelSize() const;
  unsigned int NLocalFields() const {return(_local.tsize());}
  unsigned int GlobalN() const {return(_global.size());}
  std::vector<double> GetGlobal() const {return(_global);}
  std::vector<boost::shared_ptr<BASISFIELD::basisfield> > GetLocalAsSplinefieldVector(const std::vector<unsigned int>&  ksp) const;
  BASISFIELD::splinefield GetLocalAsSingleSplinefield(const std::vector<unsigned int>&  ksp,
                                                      unsigned int                      indx=0) const;  
protected:
  void common_read(const std::string& fname);
private:
  std::string                                                _fname;
  std::vector<double>                                        _global;
  bool                                                       _has_global;
  NEWIMAGE::volume4D<float>                                  _local;
  bool                                                       _has_local;
};

} // End namespace FNIRT

#endif // end #ifndef intensity_mappers_h
