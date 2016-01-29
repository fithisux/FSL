// Declarations of classes and functions for the
// gps utility. The gps utility will determine a 
// set of isotropic directions on the sphere.
//
// gps.h
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2011 University of Oxford 
//

#ifndef gps_h
#define gps_h

#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include "newmat.h"
#include "utils/options.h"
#include "miscmaths/miscmaths.h"

namespace GPS {

class GpsException: public std::exception
{
private:
  std::string m_msg;
public:
  GpsException(const std::string& msg) throw(): m_msg(msg) { cout << what() << endl; }

  virtual const char * what() const throw() {
    return string("gps: msg=" + m_msg).c_str();
  }

  ~GpsException() throw() {}
};

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class Bvec
//
// This class manages a single diffusion direction. It utilies a dual
// representation in cartesian and spherical coordinates and these
// are kept in accordance at all times. It also supplies the partial
// derivatives of the cartesion represenation w.r.t. to the angles 
// of the spherical representation.
// scan. 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class Bvec
{
public:
  Bvec(); // Initialises to a random vector on the unity sphere
  Bvec(const NEWMAT::ColumnVector& vec);
  void SetAngles(const NEWMAT::ColumnVector& angl);
  void SetVector(const NEWMAT::ColumnVector& vec);
  void SetRandomVector();
  void SignSwap();
  const NEWMAT::ColumnVector& Vector() const { return(_vr); }
  const NEWMAT::ColumnVector& Angles() const { return(_ar); }
  const NEWMAT::ColumnVector& dVec_dtheta() const { if (!_dutd) update_derivatives(); return(_dvdt); }
  const NEWMAT::ColumnVector& dVec_dphi() const { if (!_dutd) update_derivatives(); return(_dvdp); }
  const NEWMAT::ColumnVector& d2Vec_dtheta2() const { if (!_dutd) update_derivatives(); return(_d2vdt2); }
  const NEWMAT::ColumnVector& d2Vec_dtheta_dphi() const { if (!_dutd) update_derivatives(); return(_d2vdtdp); }
  const NEWMAT::ColumnVector& d2Vec_dphi2() const { if (!_dutd) update_derivatives(); return(_d2vdp2); }

private:
  NEWMAT::ColumnVector _vr;                // Vector representation
  NEWMAT::ColumnVector _ar;                // Spherical angles representation
  mutable NEWMAT::ColumnVector _dvdt;      // Partial derivative of vector representation w.r.t. theta
  mutable NEWMAT::ColumnVector _dvdp;      // Partial derivative of vector representation w.r.t. phi
  mutable NEWMAT::ColumnVector _d2vdt2;    // 2nd derivative, see above
  mutable NEWMAT::ColumnVector _d2vdp2;    // 2nd derivative, see above
  mutable NEWMAT::ColumnVector _d2vdtdp;   // 2nd derivative, see above
  mutable bool                 _dutd;      // Indicates if Derivatives are Up To Date

  static double PI;

  void set_ar_from_vr();
  void set_vr_from_ar();
  void update_derivatives() const; 
};

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class BvecCollection
//
// This class manages a set of Bvecs. The first of these should 
// always be [1 0 0] and the remaining can take any (unity norm)
// values. Its methods supplies the total Coulomb forces (cf) of a 
// system of positive charges at both ends of all directions. It 
// also  supplies the gradient and the Hessian of the cf w.r.t. 
// to all "free" angles, noting that the first direction is fixed to 
// [1 0 0]. 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class BvecCollection
{
public:
  BvecCollection(unsigned int n);
  BvecCollection(NEWMAT::Matrix mat);
  unsigned int NBvec() const { return(_bv.size()); }
  unsigned int NPar() const { return(2*(_bv.size()-1)); }
  NEWMAT::ColumnVector GetPar() const;
  void SetPar(const NEWMAT::ColumnVector& p);
  void ShakeIt(unsigned int ns=50);
  void OptimiseOnWholeSphere(unsigned int niter=50, bool verbose=false);
  const Bvec& GetBvec(unsigned int i) const { if (i>=NBvec()) GpsException("BvecCollection::GetBvec: index out of range"); return(_bv[i]); }
  void SetBvec(unsigned int i, const NEWMAT::ColumnVector& vec) {
    if (!i || i>=NBvec()) GpsException("BvecCollection::GetBvec: index out of range");
    _bv[i].SetVector(vec);
  }
  void SetBvecToRandom(unsigned int i) {
    if (!i || i>=NBvec()) GpsException("BvecCollection::GetBvec: index out of range");
    _bv[i].SetRandomVector();
  }
  NEWMAT::ReturnMatrix GetAllBvecs() const;
  double CoulombForces() const;
  NEWMAT::ReturnMatrix CoulombForcesGradient() const;
  NEWMAT::ReturnMatrix CoulombForcesNumericalGradient();
  NEWMAT::ReturnMatrix CoulombForcesHessian() const;
  NEWMAT::ReturnMatrix CoulombForcesNumericalHessian();
  double SingleChargeCoulombForces() const;
private:
  std::vector<Bvec> _bv;

  double sqr(double a) const { return(a*a); }
};

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class GpsCF
//
// This class is effectively a wrapper class for BvecCollection that
// is a sub-class of NonlinCF which means that it can be used with
// the non-linear optimisation methods in MISCMATHS.
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class GpsCF: public MISCMATHS::NonlinCF
{
public:
  GpsCF(const BvecCollection& bc, bool verbose) : _bc(bc), _verbose(verbose) {}
  unsigned int NPar() const { return(_bc.NPar()); }
  BvecCollection GetBvecs() const { return(_bc); }

  // cost-function (Coloumb forces), gradient and Hessian as specified in NonlinCF class

  double cf(const NEWMAT::ColumnVector& p) const {
    if (p.Nrows() != int(NPar())) throw GpsException("GpsCF::cf: Mismatch between p and number of parameters");
    _bc.SetPar(p); 
    double cfv = _bc.CoulombForces();
    if (_verbose) cout << "cf = " << cfv << endl;
    return(cfv);
  }
  NEWMAT::ReturnMatrix grad(const NEWMAT::ColumnVector& p) const {
    if (p.Nrows() != int(NPar())) throw GpsException("GpsCF::grad: Mismatch between p and number of parameters");
    _bc.SetPar(p); return(_bc.CoulombForcesGradient());
  }
  boost::shared_ptr<MISCMATHS::BFMatrix> hess(const NEWMAT::ColumnVector& p,
					      boost::shared_ptr<MISCMATHS::BFMatrix>  iptr=boost::shared_ptr<MISCMATHS::BFMatrix>()) const
  {
    if (p.Nrows() != int(NPar())) throw GpsException("GpsCF::hess: Mismatch between p and number of parameters");
    _bc.SetPar(p); return(boost::shared_ptr<MISCMATHS::BFMatrix>(new MISCMATHS::FullBFMatrix(_bc.CoulombForcesHessian())));
  }
private:
  mutable BvecCollection _bc;
  bool                   _verbose;  
};

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class GpsCommandLineOptions
//
// This class parses the input to the gps utility.
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class GpsCommandLineOptions {
public:
  GpsCommandLineOptions(int argc, char *argv[]);
  std::string OutFname() const { return(_out.value()); }
  unsigned int NDir() const { return(static_cast<unsigned int>(_ndir.value())); }
  bool OptimiseOnWholeSphere() const { return(_optws.value()); }
  bool Verbose() const { return(_verbose.value()); }
  bool Debug() const { return(_debug.value()); }
  bool HasInitMatrix() const { return((_init.value().size()) ? true : false); }
  NEWMAT::Matrix InitMatrix() const { return(_init_mat); }
  bool Report() const { return(_report.value()); }
private:
  std::string                   _title;
  std::string                   _examples;
  Utilities::Option<bool>       _verbose;
  Utilities::Option<bool>       _help;
  Utilities::Option<string>     _out;
  Utilities::Option<int>        _ndir;
  Utilities::Option<bool>       _optws;
  Utilities::Option<int>        _ranseed;
  Utilities::HiddenOption<bool> _debug;
  Utilities::Option<string>     _init;
  Utilities::Option<bool>       _report;
  NEWMAT::Matrix                _init_mat;
};

} // End namespace GPS

#endif // End #ifndef gps_h 
