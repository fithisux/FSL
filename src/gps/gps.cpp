// gps - FMRIB's Tool for getting directions
//
// gps.cpp
//
// Jesper Andersson, FMRIB Image Analysis Group
//
//

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/shared_ptr.hpp>
#include "newmat.h"
#include "newmatio.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/nonlin.h"
#include "gps.h"

using namespace GPS;

double Bvec::PI = 3.141592653589793;

Bvec::Bvec() : _vr(3), _ar(2), _dvdt(3), _dvdp(3), _d2vdt2(3), _d2vdp2(3), _d2vdtdp(3), _dutd(false)
{
  for (int i=1; i<=3; i++) _vr(i) = -1.0 + 2.0*(double(rand()) / double(RAND_MAX));
  _vr /= _vr.NormFrobenius();
  set_ar_from_vr();
}

Bvec::Bvec(const NEWMAT::ColumnVector& vec) : _vr(3), _ar(2), _dvdt(3), _dvdp(3), _d2vdt2(3), _d2vdp2(3), _d2vdtdp(3), _dutd(false)
{
  if (vec.Nrows() != 3) throw GpsException("Bvec::Bvec: vec must have 3 elements");
  _vr = vec/vec.NormFrobenius(); set_ar_from_vr();
}

void Bvec::SetAngles(const NEWMAT::ColumnVector& angl) 
{
  if (angl.Nrows() != 2) throw GpsException("Bvec::SetAngles: angl must have 2 elements");
  if (_ar != angl) { _ar = angl; set_vr_from_ar(); _dutd = false; }
}

void Bvec::SetVector(const NEWMAT::ColumnVector& vec) 
{
  if (vec.Nrows() != 3) throw GpsException("Bvec::SetVector: vec must have 3 elements");
  if (_vr != vec) { _vr = vec/vec.NormFrobenius(); set_ar_from_vr(); _dutd = false; }
}

void Bvec::SetRandomVector() 
{
  for (int i=1; i<=3; i++) _vr(i) = -1.0 + 2.0*(double(rand()) / double(RAND_MAX));
  _vr /= _vr.NormFrobenius();
  set_ar_from_vr();
  _dutd = false;
}

void Bvec::SignSwap() 
{
  for (int i=1; i<=3; i++) _vr(i) = -_vr(i);
  set_ar_from_vr();
  _dutd = false;
}

void Bvec::set_ar_from_vr() {
  double theta = asin(_vr(3));
  double cost = cos(theta);
  std::vector<double> phi1(2);
  phi1[0] = acos(_vr(1)/cost); phi1[1] = 2*PI-phi1[0];
  std::vector<double> phi2(2);
  double tmp = asin(_vr(2)/cost);
  phi2[1] = PI-tmp;
  if (tmp >= 0.0) phi2[0] = tmp;
  else phi2[0] = 2.0*PI+tmp;
  double phi = 0.0;
  if (fabs(phi1[0]-phi2[0])<1e-6 || fabs(phi1[0]-phi2[1])<1e-6) phi = phi1[0];
  else phi = phi1[1];
  _ar(1) = theta;
  _ar(2) = phi;
}

void Bvec::set_vr_from_ar() {
  double cost = cos(_ar(1));
  _vr(1) = cost*cos(_ar(2)); _vr(2) = cost*sin(_ar(2)); _vr(3) = sin(_ar(1)); 
}

void Bvec::update_derivatives() const {
  double sint = sin(_ar(1)); double cost = cos(_ar(1));
  double sinp = sin(_ar(2)); double cosp = cos(_ar(2));
  _dvdt(1) = -sint*cosp; _dvdt(2) = -sint*sinp; _dvdt(3) = cost;
  _dvdp(1) = -cost*sinp; _dvdp(2) = cost*cosp; _dvdp(3) = 0.0;
  _d2vdt2(1) = -cost*cosp; _d2vdt2(2) = -cost*sinp; _d2vdt2(3) = -sint;
  _d2vdp2(1) = -cost*cosp; _d2vdp2(2) = -cost*sinp; _d2vdp2(3) = 0.0;
  _d2vdtdp(1) = sint*sinp; _d2vdtdp(2) = -sint*cosp; _d2vdtdp(3) = 0.0;
  _dutd = true;
}

BvecCollection::BvecCollection(unsigned int n) : _bv(n)
{
  NEWMAT::ColumnVector first(3);
  first << 1 << 0 << 0;
  _bv[0].SetVector(first);
  for (unsigned int i=1; i<n; i++) _bv[i].SetRandomVector(); 
}

BvecCollection::BvecCollection(NEWMAT::Matrix mat)
{
  if (mat.Nrows() != 3 && mat.Ncols() != 3) throw GpsException("BvecCollection::BvecCollection: mat must be Nx3 or 3xN");
  if (mat.Nrows() != 3) mat = mat.t();
  if (mat(1,1) != 1 || mat(2,1) != 0 || mat(3,1) != 0) throw GpsException("BvecCollection::BvecCollection: the first direction must be [1 0 0]");
  _bv.resize(mat.Ncols());
  for (unsigned int i=0; i<_bv.size(); i++) {
    _bv[i].SetVector(mat.Column(i+1));
  }
}

NEWMAT::ColumnVector BvecCollection::GetPar() const
{
  NEWMAT::ColumnVector par(NPar());
  for (unsigned int i=1; i<NBvec(); i++) {
    par.Rows(2*i-1,2*i) = _bv[i].Angles();
  }
  return(par);
}

void BvecCollection::SetPar(const NEWMAT::ColumnVector& p)
{
  if (p.Nrows() != int(NPar())) throw GpsException("BvecCollection::Set: p of wrong size");
  for (unsigned int i=1; i<NBvec(); i++) {
    _bv[i].SetAngles(p.Rows(2*i-1,2*i));
  }
}

void BvecCollection::ShakeIt(unsigned int n) 
{
  for (unsigned int i=0; i<n; i++) { 
    NEWMAT::ColumnVector grad = CoulombForcesGradient();
    double mgrad = grad.SumAbsoluteValue()/double(grad.Nrows());
    for (unsigned int j=1; j<NBvec(); j++) {
      if (grad(2*j-1) > 5*mgrad || grad(2*j) > 5*mgrad) {
        double oldcf = CoulombForces();
	NEWMAT::ColumnVector oldvec = _bv[j].Vector();
        SetBvecToRandom(j);
        if (CoulombForces() > oldcf) SetBvec(j,oldvec);
      }
    }
  }
}

void BvecCollection::OptimiseOnWholeSphere(unsigned int niter, bool verbose) 
{  
  double cf = SingleChargeCoulombForces();
  for (unsigned int i=0; i<niter; i++) {
    double itercf = cf;
    if (verbose) cout << "iter = " << i << ", cf = " << cf << endl;
    for (unsigned int j=1; j<NBvec(); j++) {
      _bv[j].SignSwap();
      double ncf = SingleChargeCoulombForces();
      if (ncf > cf) _bv[j].SignSwap(); // Swap back if it was no good
      else cf = ncf;
    }
    if (cf == itercf) break; // If it got no better
  }
}

NEWMAT::ReturnMatrix BvecCollection::GetAllBvecs() const 
{
  NEWMAT::Matrix mat(NBvec(),3);
  for (unsigned int i=0; i<NBvec(); i++) {
    mat.Row(i+1) = _bv[i].Vector().t();
  }
  mat.Release();
  return(mat);
}

double BvecCollection::CoulombForces() const
{
  double cf = 0.0;
  for (unsigned int i=0; i<NBvec()-1; i++) {
    for (unsigned int j=i+1; j<NBvec(); j++) {
      cf += 1.0 / (_bv[i].Vector()-_bv[j].Vector()).SumSquare();
      cf += 1.0 / (_bv[i].Vector()+_bv[j].Vector()).SumSquare();
    }
  }
  return(cf);  
}

NEWMAT::ReturnMatrix BvecCollection::CoulombForcesGradient() const
{
  NEWMAT::ColumnVector grad(2*(NBvec()-1));
  grad = 0.0;  
  for (unsigned int i=1; i<NBvec(); i++) {
    for (unsigned int j=0; j<NBvec(); j++) {
      if (i != j) {
	NEWMAT::ColumnVector tmp = _bv[i].Vector() - _bv[j].Vector();
        grad(2*i-1) -= 2.0*NEWMAT::DotProduct(_bv[i].dVec_dtheta(),tmp)/sqr(tmp.SumSquare()); // df_dtheta
        grad(2*i) -= 2.0*NEWMAT::DotProduct(_bv[i].dVec_dphi(),tmp)/sqr(tmp.SumSquare());     // df_dphi
        tmp = _bv[i].Vector() + _bv[j].Vector();
        grad(2*i-1) -= 2.0*NEWMAT::DotProduct(_bv[i].dVec_dtheta(),tmp)/sqr(tmp.SumSquare()); // df_dtheta
        grad(2*i) -= 2.0*NEWMAT::DotProduct(_bv[i].dVec_dphi(),tmp)/sqr(tmp.SumSquare());     // df_dphi
      }
    }
  }
  grad.Release();
  return(grad);  
}

NEWMAT::ReturnMatrix BvecCollection::CoulombForcesNumericalGradient()
{
  NEWMAT::ColumnVector grad(2*(NBvec()-1));
  grad = 0.0;
  double cf = CoulombForces();
  for (unsigned int i=1; i<NBvec(); i++) {
    NEWMAT::ColumnVector angles = _bv[i].Angles();
    angles(1) += 1e-6;
    _bv[i].SetAngles(angles);
    grad(2*i-1) = (CoulombForces()-cf)/1e-6;
    angles(1) -= 1e-6;
    angles(2) += 1e-6;
    _bv[i].SetAngles(angles);
    grad(2*i) = (CoulombForces()-cf)/1e-6;
    angles(2) -= 1e-6;
    _bv[i].SetAngles(angles);    
  }
  grad.Release();
  return(grad);
}

NEWMAT::ReturnMatrix BvecCollection::CoulombForcesHessian() const
{
  NEWMAT::Matrix hess(2*(NBvec()-1),2*(NBvec()-1));
  hess = 0.0;
  for (unsigned int i=1; i<NBvec(); i++) {
    for (unsigned int j=i; j<NBvec(); j++) {
      if (i == j) {
        for (unsigned int k=0; k<NBvec(); k++) {
          if (i != k) {
	    NEWMAT::ColumnVector tmp1 = _bv[i].Vector() - _bv[k].Vector();
            double tmp2 = tmp1.SumSquare();
            double tmp3 = tmp2*tmp2*tmp2;
            hess(2*i-1,2*j-1) += (8*sqr(NEWMAT::DotProduct(_bv[i].dVec_dtheta(),tmp1)) 
				  - 2*(NEWMAT::DotProduct(_bv[i].d2Vec_dtheta2(),tmp1) + NEWMAT::DotProduct(_bv[i].dVec_dtheta(),_bv[i].dVec_dtheta()))*tmp2) / tmp3;
            hess(2*i,2*j) += (8*sqr(NEWMAT::DotProduct(_bv[i].dVec_dphi(),tmp1)) 
			      - 2*(NEWMAT::DotProduct(_bv[i].d2Vec_dphi2(),tmp1) + NEWMAT::DotProduct(_bv[i].dVec_dphi(),_bv[i].dVec_dphi()))*tmp2) / tmp3;
            hess(2*i-1,2*j) += (8*NEWMAT::DotProduct(_bv[i].dVec_dtheta(),tmp1)*NEWMAT::DotProduct(_bv[i].dVec_dphi(),tmp1) 
				- 2*(NEWMAT::DotProduct(_bv[i].d2Vec_dtheta_dphi(),tmp1) + NEWMAT::DotProduct(_bv[i].dVec_dtheta(),_bv[i].dVec_dphi()))*tmp2) / tmp3;
	    tmp1 = _bv[i].Vector() + _bv[k].Vector();
            tmp2 = tmp1.SumSquare();
            tmp3 = tmp2*tmp2*tmp2;
            hess(2*i-1,2*j-1) += (8*sqr(NEWMAT::DotProduct(_bv[i].dVec_dtheta(),tmp1)) 
				  - 2*(NEWMAT::DotProduct(_bv[i].d2Vec_dtheta2(),tmp1) + NEWMAT::DotProduct(_bv[i].dVec_dtheta(),_bv[i].dVec_dtheta()))*tmp2) / tmp3;
            hess(2*i,2*j) += (8*sqr(NEWMAT::DotProduct(_bv[i].dVec_dphi(),tmp1)) 
			      - 2*(NEWMAT::DotProduct(_bv[i].d2Vec_dphi2(),tmp1) + NEWMAT::DotProduct(_bv[i].dVec_dphi(),_bv[i].dVec_dphi()))*tmp2) / tmp3;
            hess(2*i-1,2*j) += (8*NEWMAT::DotProduct(_bv[i].dVec_dtheta(),tmp1)*NEWMAT::DotProduct(_bv[i].dVec_dphi(),tmp1) 
				- 2*(NEWMAT::DotProduct(_bv[i].d2Vec_dtheta_dphi(),tmp1) + NEWMAT::DotProduct(_bv[i].dVec_dtheta(),_bv[i].dVec_dphi()))*tmp2) / tmp3;
	    hess(2*i,2*j-1) = hess(2*i-1,2*j);
	  }
	} 
      }
      else {
	NEWMAT::ColumnVector tmp1 = _bv[i].Vector() - _bv[j].Vector();
	double tmp2 = tmp1.SumSquare();
	double tmp3 = tmp2*tmp2*tmp2;
        hess(2*i-1,2*j-1) = (2*NEWMAT::DotProduct(_bv[i].dVec_dtheta(),_bv[j].dVec_dtheta())*tmp2 - 
			     8*NEWMAT::DotProduct(_bv[i].dVec_dtheta(),tmp1)*NEWMAT::DotProduct(_bv[j].dVec_dtheta(),tmp1)) / tmp3;
        hess(2*i,2*j) = (2*NEWMAT::DotProduct(_bv[i].dVec_dphi(),_bv[j].dVec_dphi())*tmp2 - 
			 8*NEWMAT::DotProduct(_bv[i].dVec_dphi(),tmp1)*NEWMAT::DotProduct(_bv[j].dVec_dphi(),tmp1)) / tmp3;
        hess(2*i-1,2*j) = (2*NEWMAT::DotProduct(_bv[i].dVec_dtheta(),_bv[j].dVec_dphi())*tmp2 - 
			   8*NEWMAT::DotProduct(_bv[i].dVec_dtheta(),tmp1)*NEWMAT::DotProduct(_bv[j].dVec_dphi(),tmp1)) / tmp3;
        hess(2*i,2*j-1) = (2*NEWMAT::DotProduct(_bv[i].dVec_dphi(),_bv[j].dVec_dtheta())*tmp2 - 
			   8*NEWMAT::DotProduct(_bv[i].dVec_dphi(),tmp1)*NEWMAT::DotProduct(_bv[j].dVec_dtheta(),tmp1)) / tmp3;
	tmp1 = _bv[i].Vector() + _bv[j].Vector();
	tmp2 = tmp1.SumSquare();
	tmp3 = tmp2*tmp2*tmp2;
        hess(2*i-1,2*j-1) -= (2*NEWMAT::DotProduct(_bv[i].dVec_dtheta(),_bv[j].dVec_dtheta())*tmp2 - 
			     8*NEWMAT::DotProduct(_bv[i].dVec_dtheta(),tmp1)*NEWMAT::DotProduct(_bv[j].dVec_dtheta(),tmp1)) / tmp3;
        hess(2*i,2*j) -= (2*NEWMAT::DotProduct(_bv[i].dVec_dphi(),_bv[j].dVec_dphi())*tmp2 - 
			 8*NEWMAT::DotProduct(_bv[i].dVec_dphi(),tmp1)*NEWMAT::DotProduct(_bv[j].dVec_dphi(),tmp1)) / tmp3;
        hess(2*i-1,2*j) -= (2*NEWMAT::DotProduct(_bv[i].dVec_dtheta(),_bv[j].dVec_dphi())*tmp2 - 
			   8*NEWMAT::DotProduct(_bv[i].dVec_dtheta(),tmp1)*NEWMAT::DotProduct(_bv[j].dVec_dphi(),tmp1)) / tmp3;
        hess(2*i,2*j-1) -= (2*NEWMAT::DotProduct(_bv[i].dVec_dphi(),_bv[j].dVec_dtheta())*tmp2 - 
			   8*NEWMAT::DotProduct(_bv[i].dVec_dphi(),tmp1)*NEWMAT::DotProduct(_bv[j].dVec_dtheta(),tmp1)) / tmp3;
	hess(2*j-1,2*i-1) = hess(2*i-1,2*j-1);
	hess(2*j,2*i-1) = hess(2*i-1,2*j);
	hess(2*j-1,2*i) = hess(2*i,2*j-1);
	hess(2*j,2*i) = hess(2*i,2*j);
      }
    }
  }
  hess.Release();
  return(hess);
}

NEWMAT::ReturnMatrix BvecCollection::CoulombForcesNumericalHessian()
{
  NEWMAT::Matrix hess(2*(NBvec()-1),2*(NBvec()-1));
  hess = 0.0;
  NEWMAT::ColumnVector grad = CoulombForcesNumericalGradient();
  for (unsigned int i=1; i<NBvec(); i++) {
    NEWMAT::ColumnVector angles = _bv[i].Angles();
    angles(1) += 1e-6;
    _bv[i].SetAngles(angles);
    hess.Column(2*i-1) = (CoulombForcesNumericalGradient()-grad)/1e-6;
    angles(1) -= 1e-6;
    angles(2) += 1e-6;
    _bv[i].SetAngles(angles);
    hess.Column(2*i) = (CoulombForcesNumericalGradient()-grad)/1e-6;
    angles(2) -= 1e-6;
    _bv[i].SetAngles(angles);    
  }
  hess.Release();
  return(hess);
}

double BvecCollection::SingleChargeCoulombForces() const
{
  double cf = 0.0;
  for (unsigned int i=0; i<NBvec()-1; i++) {
    for (unsigned int j=i+1; j<NBvec(); j++) {
      cf += 1.0 / (_bv[i].Vector()-_bv[j].Vector()).SumSquare();
    }
  }
  return(cf);  
}

GpsCommandLineOptions::GpsCommandLineOptions(int argc, char *argv[]) :
  _title("gps (Version 1.0)\nCopyright(c) 2011, University of Oxford (Jesper Andersson)"),
  _examples("gps --ndir=128 --optws --out=my_128_directions.txt"),
  _verbose(string("-v,--verbose"),false,string("switch on diagnostic messages"),false, Utilities::no_argument), 
  _help(string("-h,--help"),false,string("display this message"),false, Utilities::no_argument),
  _out(string("--out"),string(""),string("Name of output file (default: bvecs#.txt)"),false,Utilities::requires_argument), 
  _ndir(string("--ndir"),0,string("Number of directions"),true,Utilities::requires_argument),
  _optws(string("--optws"),false,string("Perform additional optimisation on the whole sphere (needed for eddy)"),false,Utilities::no_argument),
  _ranseed(string("--ranseed"),0,string("Seed random generator with supplied number"),false,Utilities::requires_argument),
  _debug(string("--debug"),false,string("Writes out some debug info and terminates"),false,Utilities::no_argument),
  _init(string("--init"),string(""),string("File with bvecs to use as initialisation"),false,Utilities::requires_argument),
  _report(string("--report"),false,string("Report coulomb forces for inital configuration"),false,Utilities::no_argument)
{
  // Parse arguments
  Utilities::OptionParser options(_title,_examples);
  try {
    options.add(_out);
    options.add(_ndir);
    options.add(_optws);
    options.add(_ranseed);
    options.add(_init);
    options.add(_report);
    options.add(_verbose);
    options.add(_help);
    options.add(_debug);

    int i=options.parse_command_line(argc, argv);
    if (i < argc) {
      for (; i<argc; i++) {
        cerr << "Unknown input: " << argv[i] << endl;
      }
      exit(EXIT_FAILURE);
    }

    if (_help.value() || !options.check_compulsory_arguments(true)) {
      options.usage();
      exit(EXIT_FAILURE);
    }
  }  
  catch(Utilities::X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } 
  catch(std::exception &e) {
    cerr << e.what() << endl;
    exit(EXIT_FAILURE);
  } 
  // Do a bit of sanity checking
  if (_ndir.value() > 2048) throw GpsException("");
  if (!_out.set()) {
    char ofname[128];
    sprintf(ofname,"bvecs%d.txt",_ndir.value());
    _out.set_value(std::string(ofname));
  }
  if (_init.set()) {
    try { _init_mat = MISCMATHS::read_ascii_matrix(_init.value()); }
    catch(...) { throw GpsException("GpsCommandLineOptions::GpsCommandLineOptions: Failed to read --init file"); }
    if (_init_mat.Nrows() != 3 && _init_mat.Ncols() != 3) throw GpsException("GpsCommandLineOptions::GpsCommandLineOptions: --init matrix must be 3xN or Nx3");
    if (_init_mat.Nrows() != 3) _init_mat = _init_mat.t();
  }
    
  // Perform some global tasks
  srand(_ranseed.value());
}

int main(int argc, char *argv[])
{
  BvecCollection *bvecp;
  // Parse command line
  GpsCommandLineOptions clo(argc,argv);
  

  if (clo.HasInitMatrix()) {
    NEWMAT::Matrix mat_bvecs = clo.InitMatrix();
    bvecp = new BvecCollection(mat_bvecs);
  }
  else {
    // Get a random collection of b-vectors
    bvecp = new BvecCollection(clo.NDir());
    // Do a little initial shaking
    bvecp->ShakeIt(50);
  }

  // If debug, write some debug info and exit
  if (clo.Debug()) {
    cout << "Coulomb Forces = " << bvecp->CoulombForces() << endl;
    cout << "Single Charge Coulomb Forces = " << bvecp->SingleChargeCoulombForces() << endl;
    NEWMAT::Matrix mat_bvecs = bvecp->GetAllBvecs();
    MISCMATHS::write_ascii_matrix(mat_bvecs,clo.OutFname()+string("_debug_vectors"));
    NEWMAT::ColumnVector grad = bvecp->CoulombForcesGradient();
    MISCMATHS::write_ascii_matrix(grad,clo.OutFname()+string("_debug_analytic_gradient"));
    grad = bvecp->CoulombForcesNumericalGradient();
    MISCMATHS::write_ascii_matrix(grad,clo.OutFname()+string("_debug_numerical_gradient"));
    NEWMAT::Matrix hess = bvecp->CoulombForcesHessian();
    MISCMATHS::write_ascii_matrix(hess,clo.OutFname()+string("_debug_analytic_hessian"));
    hess = bvecp->CoulombForcesNumericalHessian();
    MISCMATHS::write_ascii_matrix(hess,clo.OutFname()+string("_debug_numerical_hessian"));
    exit(EXIT_SUCCESS);    
  }

  // If report, report forces and exit
  if (clo.Report()) {
    cout << "Coulomb Forces = " << bvecp->CoulombForces() << endl;
    cout << "Single Charge Coulomb Forces = " << bvecp->SingleChargeCoulombForces() << endl;
    exit(EXIT_SUCCESS);
  }

  // Set up a cost-function object
  GpsCF cf(*bvecp,clo.Verbose());

  // A nonlinear object to guide the optimisation
  MISCMATHS::NonlinParam nlpar(bvecp->NPar(),MISCMATHS::NL_LM,bvecp->GetPar());
  nlpar.SetGaussNewtonType(MISCMATHS::LM_L);
  nlpar.SetMaxIter(2000*clo.NDir());

  // Do the minimisation
  try {
    nonlin(nlpar,cf);
  }
  catch (const std::exception& error) {
    cerr << "Error occurred when attempting to optimise directions" << endl;
    cerr << "Exception thrown with message: " << error.what() << endl; 
    exit(EXIT_FAILURE);
  }
  BvecCollection opt_bvecs = cf.GetBvecs();

  // At this stage they should be optimised on a half sphere, and
  // distributed "randomly" on the whole sphere. If the --optws
  // flag is set we will sign swap directions and attempt do find
  // a more even configuration on the whole sphere. To do so we
  // use the coulomb forces resulting from putting charges only at
  // the positive ends of the gradients.
  if (clo.OptimiseOnWholeSphere()) opt_bvecs.OptimiseOnWholeSphere(50,clo.Verbose());

  // Write the results
  NEWMAT::Matrix mat_bvecs = opt_bvecs.GetAllBvecs();
  MISCMATHS::write_ascii_matrix(mat_bvecs,clo.OutFname());

  exit(EXIT_SUCCESS);
}
