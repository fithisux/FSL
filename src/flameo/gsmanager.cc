/*  gsmanager.cc

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

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

#include "gsmanager.h"
#include "utils/log.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include "newimage/newimageall.h"
#include "utils/tracer_plus.h"
//#include "mcmc.h"
#include "mcmc_mh.h"
#include "miscmaths/t2z.h"
#include "miscmaths/f2z.h"
#include <algorithm>
#include <limits.h>

using namespace Utilities;
using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace std;

namespace Gs {
  
  bool compare(const pair<float,int> &r1,const pair<float,int> &r2){
    return (r1.first<r2.first);
  }

  void Gsmanager::randomise(vector< pair<float,int> >& r){
    for(unsigned int i=0;i<r.size();i++){
    pair<float,int> p(rand()/float(RAND_MAX),i);
    r[i]=p;
    }
    sort(r.begin(),r.end(),compare);
  }

  vector< pair<float,int> > Gsmanager::randomise(const int n){
    vector< pair<float,int> > v(n);
    randomise(v);
    return v;
  }
  
  void Gsmanager::do_kmeans(const Matrix& data,vector<int>& z,const int k,Matrix& means){
    Tracer tr("Gsmanager::do_kmeans");  

    // note that first class is restricted to having a mean of zero

    int numiter=100;
    if(data.Nrows() != (int)z.size())
      z.resize(data.Nrows());
    
    // k is number of classes
    int n = data.Nrows(); // number of observations
    int d = data.Ncols(); // dimension
    
    //    Matrix means(d,k),newmeans(d,k);
    Matrix newmeans(d,k);
    ColumnVector nmeans(k);
       
    nmeans=0;
	
 //    //    cout<<"inside kmeans"<<endl;
//     // initialise random
//     vector< pair<float,int> > rindex(n);
//     randomise(rindex);
//     vector<pair<float,int> >::iterator riter;
//     int nn=0,cl=1,nperclass=(int)(float(n)/float(k));
//     for(riter=rindex.begin();riter!=rindex.end();++riter){
//       means.Column(cl) += data.Row((*riter).second+1).t();
//       nmeans(cl) += 1;
//       z[(*riter).second]=cl;
//       nn++;
//       if(nn>=nperclass && cl<k){
// 	nn=0;
// 	cl++;
//       }
//     }
//     for(int m=1;m<=k;m++)
//       means.Column(m) /= nmeans(m);


    // iterate
    for(int iter=0;iter<numiter;iter++){
      // loop over datapoints and attribute z for closest mean
      newmeans=0;
    nmeans=0;
    for(int i=1;i<=n;i++){
      float mindist=1E20,dist=0;
      int mm=1;
      for(int m=1;m<=k;m++){
	dist = (means.Column(m)-data.Row(i).t()).SumSquare();
	if( dist<mindist){
	  mindist=dist;
	  mm = m;
	}
      }
      z[i-1] = mm;
      newmeans.Column(mm) += data.Row(i).t();
      nmeans(mm) += 1;
    }
    
    // compute means
    for(int m=1;m<=k;m++){
      // if(nmeans(m)==0){
      // 	do_kmeans(data,z,k);
      // 	return;
      //       }
      newmeans.Column(m) /= nmeans(m);
    }
    means = newmeans;
  }
  }


  float Gsmanager::log_likelihood(float beta, const ColumnVector& gam, const ColumnVector& Y, const Matrix& z, const ColumnVector& S)
  {
    Tracer tr("Gsmanager::log_likelihood");  

    int nts=Y.Nrows();

    DiagonalMatrix U(nts);

    float ret=0;

    for(int t=1;t<=nts;t++)
      { 
	if(beta+S(t) == 0) return 1e32;
	U(t) = S(t)+beta;
	float vr = S(t)+beta;
	ret+=lognormpdf(Y(t),(z.Row(t)*gam).AsScalar(),vr);
      }

//     OUT(beta);
//     OUT(Y.t());
//     OUT(gam);
//     OUT(ret);
    
    //float ret=std::log(mvnpdf(Y.t(),(z*gam).t(),U));

    return ret;

  }
 
  float Gsmanager::marg_posterior_energy(float x, const ColumnVector& y, const Matrix& z, const ColumnVector& S)
  {
    //Tracer tr("Gsmanager::marg_posterior_energy");  

    // x is log(variance)

    float ex=std::exp(x);

    // ex is variance

    if(ex <= 0)
      {
	float ret =1e32;
	return ret;
      }

    int nts=y.Nrows();

    DiagonalMatrix iU(nts);

    for(int t=1;t<=nts;t++)
      { 
	if(ex+S(t) == 0) return 1e32;
	iU(t) = 1.0/(S(t)+ex);
      }


    Matrix tmp = z.t()*iU;
    SymmetricMatrix ziUz;
    ziUz << tmp*z;
    ColumnVector gam=ziUz.i()*tmp*y;  

    //   OUT(iU);
    //   OUT(log(priorsum));

//     float logprior = -x;
//     float logdfxbydx = x;

    float ret = -(0.5*iU.LogDeterminant().LogValue()-0.5*ziUz.LogDeterminant().LogValue()-0.5*(y.t()*iU*y-gam.t()*ziUz*gam).AsScalar());//+logprior+logdfxbydx);

    //   OUT(ret);
    //   exit(1);
    return ret;

  }
  
  float Gsmanager::log_likelihood_outlier(float beta, float beta_outlier, const ColumnVector& gam, const ColumnVector& Y, const Matrix& z, const ColumnVector& S, float global_prob_outlier, const ColumnVector& prob_outlier)
  {
    //Tracer tr("Gsmanager::log_likelihood_outlier");  


    int nts=Y.Nrows();

    DiagonalMatrix U1(nts);
    DiagonalMatrix U2(nts);
    
    float ret=0;
    for(int t=1;t<=nts;t++)
      { 
	float vr1=S(t)+beta;
	if(vr1 <= 0) return 1e32;
	U1(t) = vr1;
	//	float vr2=S(t)+beta_outlier;
	float vr2=beta_outlier;
	if(vr2 <= 0) return 1e32;

	U2(t) = vr2;
	if(prob_outlier(t)>0.999)
	  U2(t)=0;

	ret+=std::log((1-global_prob_outlier)*normpdf(Y(t),(z.Row(t)*gam).AsScalar(),vr1)+global_prob_outlier*normpdf(Y(t),(z.Row(t)*gam).AsScalar(),vr2));
      }

//     OUT(global_prob_outlier);
//     OUT(beta);
//     OUT(beta_outlier);
//     OUT(Y.t());
//     OUT(gam);
//     OUT(ret);

    //float ret=std::log((1-global_prob_outlier)*mvnpdf(Y.t(),(z*gam).t(),U1)+global_prob_outlier*mvnpdf(Y.t(),(z*gam).t(),U2));

    return ret;

  }

  float Gsmanager::marg_posterior_energy_outlier(float log_beta, float log_beta_outlier, const ColumnVector& y, const Matrix& z, const ColumnVector& S, const ColumnVector& prob_outlier)
  {
    //   Tracer tr("Gsmanager::marg_posterior_energy_outlier");  

    float beta=std::exp(log_beta);    // beta is a variance
    float beta_outlier=std::exp(log_beta_outlier);

    if(beta <= 0 || beta_outlier <=0)
      {
	return 1e32;
      }

    DiagonalMatrix iU(y.Nrows());
    for(int t=1;t<=y.Nrows();t++)
      { 
	//double vr=Sqr(1-prob_outlier(t))*(S(t)+beta)+Sqr(prob_outlier(t))*(S(t)+beta_outlier);
	double vr=Sqr(1-prob_outlier(t))*(S(t)+beta)+Sqr(prob_outlier(t))*(beta_outlier);
// 	double vr2=(S(t)+beta);

// 	if(vr!=vr2)
// 	  {
// 	    OUT(t);
// 	    OUT(vr);
// 	    OUT(vr2);
// 	    OUT(vr-vr2);
// 	  }

	if(vr <= 0) return 1e32;
	iU(t) = 1.0/vr;

	if(prob_outlier(t)>0.999)
	  iU(t)=0;

      }

    Matrix tmp = z.t()*iU;
    SymmetricMatrix ziUz;
    ziUz << tmp*z;
    ColumnVector gam=pinv(ziUz)*tmp*y;  

//     float logprior = -log_beta-log_beta_outlier;
//     float logdfxbydx = log_beta+log_beta_outlier;

//     float ret = -(0.5*iU.LogDeterminant().LogValue()-0.5*ziUz.LogDeterminant().LogValue()-0.5*(y.t()*iU*y-gam.t()*ziUz*gam).AsScalar()+logprior+logdfxbydx);
 
    float ret = -(0.5*iU.LogDeterminant().LogValue()-0.5*ziUz.LogDeterminant().LogValue()-0.5*(y.t()*iU*y-gam.t()*ziUz*gam).AsScalar());

//    if(log_beta>2.5 && log_beta_outlier>6.06)
//      {    
//     OUT(beta);
//     OUT(beta_outlier);
//     OUT(prob_outlier);
//     OUT(gam);
//     OUT(iU);
//     OUT(ret);
//     exit(0);
//      }
    return ret;

  }

  float Gsmanager::solveforbeta(const ColumnVector& y, const Matrix& z, const ColumnVector& S)
  { 
    //    Tracer tr("Gsmanager::solveforbeta");

    // Value of golden section (1 + sqrt(5))/2.0
    double phi = 1.6180339887499;
    double cphi = 1.0 - 1.0/phi;
    double TOL = 1.0e-6;	// Maximal fractional precision
    double TINY = 1.0e-10;         // Can't use fractional precision when minimum is at 0
  
    // Bracket the minimum
    double br_min=log(1e-10);
    double br_mid=0;
    double br_max=log(1e10);

    int dir=1;
    double pt=1e-10;

    // Use Brent's algorithm to find minimum
    // Initialise the points and function values
    double w = br_mid;   	// Where second from minimum is
    double v = br_mid;   	// Previous value of w
    double x = v;   	// Where current minimum is
    double e = 0.0; 	// Distance moved on step before last

    // x is log(variance)
    double fx = marg_posterior_energy(pt+x*dir, y, z, S);

    double fv = fx;
    double fw = fx;
    double d = 0.0;
    double tol1 = TOL*abs(x) + TINY;
   
    float prec = 0.00001;
    int niters=100;
    double u = 0.0;
    for(int n = 1;n<=niters;n++)
      {
	double xm = 0.5*(br_min+br_max);
	// Make sure that tolerance is big enough
	tol1 = TOL*abs(x) + TINY;

	// Decide termination on absolute precision required by options(2)
	if (abs(x - xm) <= prec && br_max-br_min < 4*prec)
	  break;
      
	// Check if step before last was big enough to try a parabolic step.
	// Note that this will fail on first iteration, which must be a golden
	// section step.
	if (abs(e) > tol1)
	  {
	    // Construct a trial parabolic fit through x, v and w
	    double r = (fx - fv) * (x - w);
	    double q = (fx - fw) * (x - v);

	    double p = (x - v)*q - (x - w)*r;

	    q = 2.0 * (q - r);

	    if (q > 0.0 && p != 0)
	      p = -p;
	    q = abs(q);
	    // Test if the parabolic fit is OK
	  
	    if (abs(p) >= abs(0.5*q*e) || p <= q*(br_min-x) || p >= q*(br_max-x))
	      {

		// No it isn't, so take a golden section step
		if (x >= xm)
		  e = br_min-x;
		else
		  e = br_max-x;

		d = cphi*e;
	      }
	    else
	      {
	      
		// Yes it is, so take the parabolic step
		e = d;
		d = p/q;
		u = x+d;
		if (u-br_min < 2*tol1 || br_max-u < 2*tol1)
		  d = sign(xm-x)*tol1;
	      
	      }
	  }
	else
	  {
	  
	    // Step before last not big enough, so take a golden section step
	    if (x >= xm)
	      {
		e = br_min - x;
	      }
	    else
	      {
		e = br_max - x;
	      }

	    d = cphi*e;
	  }
      
	// Make sure that step is big enough
	if (abs(d) >= tol1)
	  {
	    u = x+d;
	  }
	else
	  {
	    u = x + sign(d)*tol1;
	  }


	// Evaluate function at u
	double fu = marg_posterior_energy(pt+u*dir, y, z, S);

	// Reorganise bracket
	if (fu <= fx)
	  {
	    if (u >= x)
	      br_min = x;
	    else
	      br_max = x;
	  
	    v = w; w = x; x = u;
	    fv = fw; fw = fv; fx = fu;
	  }
	else
	  {
	    if (u < x)
	      br_min = u;   
	    else
	      br_max = u;
	  
	    if (fu <= fw || w == x)
	      {
		v = w; w = u;
		fv = fw; fw = fu;
	      }
	    else if (fu <= fv || v == x || v == w)
	      {
		v = u;
		fv = fu;
	      }
	  }      
      }

    //    OUT(exp(x));

    return Max(1e-10,exp(x));
  }

  void Gsmanager::multitfit(const Matrix& x, ColumnVector& m, SymmetricMatrix& covar, float& v, bool fixmean/*=false*/) const
  {
    Tracer_Plus trace("Gsmanager::multitfit");

    int n = x.Ncols();
    int P = x.Nrows();

    v = 10;

    Matrix mx;
    if(!fixmean)
      {	
	Matrix mn;
	remmean(x,mx,mn,2);
	m = mn.t().AsColumn();
      }
    else
      {
	mx=x;
	for(int ctr = 1; ctr <= x.Nrows(); ctr++) {
	  for (int ctr2 =1; ctr2 <= x.Ncols(); ctr2++) {
	    mx(ctr,ctr2)-=m(ctr);
	  }
	}
      }

    covar = cov(mx.t());
    float tmp = pow(covar.Determinant(),(1.0/P));    

    covar=covar/tmp;

    // xsq(i) = x(i)'*inv(c)*x(i)
    ColumnVector xsq(n);
    SymmetricMatrix invcovar = covar.i();

    float cosi = 0.0;    
    for(int i = 1; i <= n; i++)
      {
	xsq(i) = (mx.Column(i).t()*invcovar*mx.Column(i)).AsScalar();
	cosi += xsq(i);
      }
    cosi /= n-1;

    float phi = tmp;

    for(int i = 1; i <= 50; i++)
      {
	float newphi = 0.0;
	for(int j = 1; j <= n; j++)
	  {
	    float tau = phi*(v+P)/(v*phi+xsq(j));
	    newphi += tau*xsq(j);
	  }
	phi = newphi/float(n*P);

	// 	write_ascii_matrix(xsq,'xsq');
	// 	OUT(P);
	// 	OUT(phi);
	//	v = mdofls(xsq,phi,P);
	v =  2.0/(1.0-phi/cosi);
	//	OUT(v);
      } 

    covar = covar*phi;
  }

  void Gsmanager::setup()
  {
    Tracer_Plus trace("Gsmanager::setup");
    
    if(opts.debuglevel.value()==2)
      {
	cout << "******************************************" << endl
	     << "SETUP" << endl << "******************************************"
	     << endl;
      }

    dofpassedin = opts.dofvarcopefile.value() != string("");

    // setup design
    design.setup();    
        
    mcmc_mask=design.getmask();

    volume<float> tmp_mask=mcmc_mask;
    tmp_mask.binarise(1e-32);
    nmaskvoxels = int(tmp_mask.sum());

    ngs = design.getngs();
    nevs = design.getnevs();
    ntpts = design.getntpts();
    xsize = design.getmask().xsize();
    ysize = design.getmask().ysize();
    zsize = design.getmask().zsize();
    
    cout << "nevs=" << nevs << endl;
    cout << "ntpts=" << ntpts << endl;
    cout << "ngs=" << ngs << endl;
    cout << "nvoxels=" << nmaskvoxels << endl;

    initialise();

  }

  void Gsmanager::initialise()
  {
    Tracer_Plus trace("Gsmanager::initialise");

    pes.resize(design.getnevs());      

    ts.resize(design.getnumtcontrasts());
    tdofs.resize(design.getnumtcontrasts());
    zts.resize(design.getnumtcontrasts());
    tcopes.resize(design.getnumtcontrasts());
    tvarcopes.resize(design.getnumtcontrasts());
    fs.resize(design.getnumfcontrasts());
    fdof1s.resize(design.getnumfcontrasts());
    fdof2s.resize(design.getnumfcontrasts());
    zfs.resize(design.getnumfcontrasts());

    zflame1lowerts.resize(design.getnumtcontrasts());
    zflame1upperts.resize(design.getnumtcontrasts());
    zflame1lowerfs.resize(design.getnumfcontrasts());
    zflame1upperfs.resize(design.getnumfcontrasts());

    beta_b.resize(design.getngs());
    beta_c.resize(design.getngs());
    beta_mean.resize(design.getngs());
    
    weights.resize(design.getngs());

    if(opts.infer_outliers.value())
      {
	beta_outlier_mean.resize(design.getngs());
	global_prob_outlier_mean.resize(design.getngs());
	prob_outlier_mean.resize(design.getngs());

      }

    if(nevs>100)
      {
	if(!(opts.runmode.value()==string("fe") || (opts.runmode.value()==string("ols") && !opts.infer_outliers.value())))
	  {
	    cout << "WARNING: Do anything other than straight ols or fixed effects with a large number of EVs will be VERY slow in the current implementation." << endl;
	  }
	else if(!opts.no_pe_output.value())
	  {
	    cout << "WARNING: Number of EVs in the design matrix is large. Consider turning the no PE output option on (--npo)." << endl;
	  }
      }

    if(!(opts.runmode.value()==string("fe") || (opts.runmode.value()==string("ols") && !opts.infer_outliers.value())))
      {
	cov_pes.reinitialize(xsize,ysize,zsize,nevs*nevs);   
	cov_pes = 0;
      }

    for(int g = 0; g < ngs; g++)
      {
	beta_mean[g].reinitialize(xsize,ysize,zsize);		
	beta_mean[g].copyproperties(design.getmask());
	beta_mean[g] = 0;
	beta_b[g].reinitialize(xsize,ysize,zsize);		
	beta_b[g] = 0;
	beta_c[g].reinitialize(xsize,ysize,zsize);		
	beta_c[g] = 0;
	weights[g].reinitialize(xsize,ysize,zsize,design.getntptsingroup(g+1));
	weights[g].copyproperties(design.getmask());
	weights[g] = 0;

	if(opts.infer_outliers.value())
	  {
	    beta_outlier_mean[g].reinitialize(xsize,ysize,zsize);
	    beta_outlier_mean[g].copyproperties(design.getmask());
	    beta_outlier_mean[g] = 0;
	    global_prob_outlier_mean[g].reinitialize(xsize,ysize,zsize);
	    global_prob_outlier_mean[g].copyproperties(design.getmask());
	    global_prob_outlier_mean[g] = 0;
	    prob_outlier_mean[g].reinitialize(xsize,ysize,zsize,design.getntptsingroup(g+1));
	    prob_outlier_mean[g].copyproperties(design.getmask());
	    prob_outlier_mean[g] = 0;

	  }
    }    

    for(int e = 0; e < nevs; e++)
      {
	pes[e].reinitialize(xsize,ysize,zsize);
	pes[e].copyproperties(design.getmask());
	pes[e] = 0;
      }

    for(int t = 0; t < design.getnumtcontrasts(); t++)
      {
	ts[t].reinitialize(xsize,ysize,zsize);
	ts[t].copyproperties(design.getmask());
	ts[t] = 0;

	tdofs[t].reinitialize(xsize,ysize,zsize);
	tdofs[t].copyproperties(design.getmask());
	tdofs[t] = 0;

	zts[t].reinitialize(xsize,ysize,zsize);
	zts[t].copyproperties(design.getmask());
	zts[t] = 0;

// 	OUT(zts[t].xdim());
// 	OUT(zts[t].ydim());
// 	OUT(zts[t].zdim());

	zflame1upperts[t].reinitialize(xsize,ysize,zsize);
	zflame1upperts[t].copyproperties(design.getmask());
	zflame1upperts[t] = 0;

	zflame1lowerts[t].reinitialize(xsize,ysize,zsize);
	zflame1lowerts[t].copyproperties(design.getmask());
	zflame1lowerts[t] = 0;

	tcopes[t].reinitialize(xsize,ysize,zsize);
	tcopes[t].copyproperties(design.getmask());
	tcopes[t] = 0;

	tvarcopes[t].reinitialize(xsize,ysize,zsize);
	tvarcopes[t].copyproperties(design.getmask());
	tvarcopes[t] = 0;
      }

    for(int f = 0; f < design.getnumfcontrasts(); f++)
      {
	fs[f].reinitialize(xsize,ysize,zsize);
	fs[f].copyproperties(design.getmask());
	fs[f] = 0;
	fdof1s[f].reinitialize(xsize,ysize,zsize);
	fdof1s[f].copyproperties(design.getmask());
	fdof1s[f] = 0;
	fdof2s[f].reinitialize(xsize,ysize,zsize);
	fdof2s[f].copyproperties(design.getmask());
	fdof2s[f] = 0;
	zfs[f].reinitialize(xsize,ysize,zsize);
	zfs[f].copyproperties(design.getmask());
	zfs[f] = 0;
	zflame1upperfs[f].reinitialize(xsize,ysize,zsize);
	zflame1upperfs[f].copyproperties(design.getmask());
	zflame1upperfs[f] = 0;	
	zflame1lowerfs[f].reinitialize(xsize,ysize,zsize);
	zflame1lowerfs[f].copyproperties(design.getmask());
	zflame1lowerfs[f] = 0;
      }
    
  }
  
  void Gsmanager::run()
  {
    Tracer_Plus trace("Gsmanager::run");

    if(opts.runmode.value()==string("fe"))
      {
	// check that varcope data is available
	if(GsOptions::getInstance().varcopefile.value() == string(""))	
	  {
	    cout << "Fixed effects requires lower level variance, this must be passed in when doing fixed effects using the --varcope option." << endl;
	    throw Exception("Fixed effects requires lower level variance, this must be passed in when doing fixed effects using the --varcope option.");
	  }
	else
	  {
	    fixed_effects();
	  }
      }    
    else if(opts.runmode.value()==string("ols"))
      {
	if(opts.infer_outliers.value())
	  {
	    // set varcopedata (var_filtered_func_data) to zero and run flame stage 1
	    design.setvarcopedata(0);
	    flame_stage1();	 
	  }
	else
	  {
	    ols();
	  }
	
      }      
    else if(opts.runmode.value()==string("flame1"))
      {	  	
	flame_stage1();	    
      }
    else if(opts.runmode.value()==string("flame12"))	
      {
	flame_stage1();	
	flame_stage2();	    
      }      
    else
      {	
	throw Exception("Invalid run mode");
      }
  }

  void Gsmanager::save()
  {
    Tracer_Plus trace("Gsmanager::save");

    // need to save residuals:
    volume4D<float> res = design.getcopedata();
    res = 0;

    for(int x = 0; x < xsize; x++)
      for(int y = 0; y < ysize; y++)
	for(int z = 0; z < zsize; z++)
	  {
	    if(design.getmask()(x,y,z))
	      {
		Matrix dm = design.getdm(x,y,z);
		
		ColumnVector petmp(nevs);
		for(int e = 0; e < nevs; e++)
		  {
		    petmp(e+1) = pes[e](x,y,z);
		  }	      
		
		// remove any zero EV PEs
		petmp = design.remove_zeroev_pes(x,y,z,petmp);
		
		res.setvoxelts(design.getcopedata().voxelts(x,y,z)-dm*petmp,x,y,z);
		  
	      }
	  }
    

    if ( opts.outputDof.value() ) {
      ColumnVector dofVec(1);
      dofVec = design.dof();
      write_ascii_matrix( LogSingleton::getInstance().appendDir("dof"), dofVec);
    }

    res.set_intent(NIFTI_INTENT_ESTIMATE, 0, 0, 0);
    res.setDisplayMaximumMinimum(res.max(),res.min());
    save_volume4D(res, LogSingleton::getInstance().appendDir("res4d"));      

    design.getmask().set_intent(NIFTI_INTENT_NONE, 0, 0, 0);
    const_cast< volume <float>& >(design.getmask()).setDisplayMaximumMinimum(design.getmask().max(),design.getmask().min());
    save_volume(design.getmask(), LogSingleton::getInstance().appendDir("mask"));

    for(int t = 0; t < design.getnumtcontrasts(); t++)
      {	
	ts[t].set_intent(NIFTI_INTENT_TTEST, 0, 0, 0);
	ts[t].setDisplayMaximumMinimum(ts[t].max(),ts[t].min());
	save_volume(ts[t],LogSingleton::getInstance().appendDir("tstat"+num2str(t+1)));

	tdofs[t].set_intent(NIFTI_INTENT_NONE, 0, 0, 0);
	tdofs[t].setDisplayMaximumMinimum(tdofs[t].max(),tdofs[t].min());
	save_volume(tdofs[t],LogSingleton::getInstance().appendDir("tdof_t"+num2str(t+1)));

	zts[t].set_intent(NIFTI_INTENT_ZSCORE, 0, 0, 0);
	zts[t].setDisplayMaximumMinimum(zts[t].max(),zts[t].min());
// 	float min; float max;
// 	OUT(zts[t].min());
// 	OUT(zts[t].max());
// 	OUT(min); OUT(max);
	save_volume(zts[t],LogSingleton::getInstance().appendDir("zstat"+num2str(t+1)));

	zflame1upperts[t].set_intent(NIFTI_INTENT_ZSCORE, 0, 0, 0);
	zflame1upperts[t].setDisplayMaximumMinimum(zflame1upperts[t].max(),zflame1upperts[t].min());
	save_volume(zflame1upperts[t],LogSingleton::getInstance().appendDir("zflame1uppertstat"+num2str(t+1)));

	zflame1lowerts[t].set_intent(NIFTI_INTENT_ZSCORE, 0, 0, 0);
	zflame1lowerts[t].setDisplayMaximumMinimum(zflame1lowerts[t].max(),zflame1lowerts[t].min());
	save_volume(zflame1lowerts[t],LogSingleton::getInstance().appendDir("zflame1lowertstat"+num2str(t+1)));

	tcopes[t].set_intent(NIFTI_INTENT_ESTIMATE, 0, 0, 0);
	tcopes[t].setDisplayMaximumMinimum(tcopes[t].max(),tcopes[t].min());
	save_volume(tcopes[t],LogSingleton::getInstance().appendDir("cope"+num2str(t+1)));

	tvarcopes[t].set_intent(NIFTI_INTENT_ESTIMATE, 0, 0, 0);
	tvarcopes[t].setDisplayMaximumMinimum(tvarcopes[t].max(),tvarcopes[t].min());
	save_volume(tvarcopes[t],LogSingleton::getInstance().appendDir("varcope"+num2str(t+1)));

      }

    for(int f = 0; f < design.getnumfcontrasts(); f++)
      {
	fs[f].set_intent(NIFTI_INTENT_FTEST, 0, 0, 0);
	fs[f].setDisplayMaximumMinimum(fs[f].max(),fs[f].min());
	save_volume(fs[f],LogSingleton::getInstance().appendDir("fstat"+num2str(f+1)));

	zflame1upperfs[f].set_intent(NIFTI_INTENT_ZSCORE, 0, 0, 0);
	zflame1upperfs[f].setDisplayMaximumMinimum(zflame1upperfs[f].max(),zflame1upperfs[f].min());
	save_volume(zflame1upperfs[f],LogSingleton::getInstance().appendDir("zflame1upperfstat"+num2str(f+1)));

	zflame1lowerfs[f].set_intent(NIFTI_INTENT_ZSCORE, 0, 0, 0);
	zflame1lowerfs[f].setDisplayMaximumMinimum(zflame1lowerfs[f].max(),zflame1lowerfs[f].min());
	save_volume(zflame1lowerfs[f],LogSingleton::getInstance().appendDir("zflame1lowerfstat"+num2str(f+1)));

	fdof1s[f].set_intent(NIFTI_INTENT_NONE, 0, 0, 0);
	fdof1s[f].setDisplayMaximumMinimum(fdof1s[f].max(),fdof1s[f].min());
	save_volume(fdof1s[f],LogSingleton::getInstance().appendDir("fdof1_f"+num2str(f+1)));

	fdof2s[f].set_intent(NIFTI_INTENT_NONE, 0, 0, 0);
	fdof2s[f].setDisplayMaximumMinimum(fdof2s[f].max(),fdof2s[f].min());
	save_volume(fdof2s[f],LogSingleton::getInstance().appendDir("fdof2_f"+num2str(f+1)));

	zfs[f].set_intent(NIFTI_INTENT_ZSCORE, 0, 0, 0);
	zfs[f].setDisplayMaximumMinimum(zfs[f].max(),zfs[f].min());
	save_volume(zfs[f],LogSingleton::getInstance().appendDir("zfstat"+num2str(f+1)));

      }

    if(!opts.no_pe_output.value())
      for(int e = 0; e < nevs; e++)
	{
	  pes[e].set_intent(NIFTI_INTENT_ESTIMATE, 0, 0, 0);
	  pes[e].setDisplayMaximumMinimum(pes[e].max(),pes[e].min());
	  save_volume(pes[e],LogSingleton::getInstance().appendDir("pe"+num2str(e+1)));
	}

    for(int g = 0; g < ngs; g++)
      {
	beta_mean[g].set_intent(NIFTI_INTENT_ESTIMATE, 0, 0, 0);
	beta_mean[g].setDisplayMaximumMinimum(beta_mean[g].max(),beta_mean[g].min());
	save_volume(beta_mean[g],LogSingleton::getInstance().appendDir("mean_random_effects_var"+num2str(g+1)));

	weights[g].set_intent(NIFTI_INTENT_ESTIMATE, 0, 0, 0);
	weights[g].setDisplayMaximumMinimum(weights[g].max(),weights[g].min());
	save_volume4D(weights[g],LogSingleton::getInstance().appendDir("weights"+num2str(g+1)));
	      
	if(opts.infer_outliers.value())
	  {
	    beta_outlier_mean[g].set_intent(NIFTI_INTENT_ESTIMATE, 0, 0, 0);
	    beta_outlier_mean[g].setDisplayMaximumMinimum(beta_outlier_mean[g].max(),beta_outlier_mean[g].min());
	    save_volume(beta_outlier_mean[g],LogSingleton::getInstance().appendDir("mean_outlier_random_effects_var"+num2str(g+1)));

	    global_prob_outlier_mean[g].set_intent(NIFTI_INTENT_ESTIMATE, 0, 0, 0);
	    global_prob_outlier_mean[g].setDisplayMaximumMinimum(global_prob_outlier_mean[g].max(),global_prob_outlier_mean[g].min());
	    save_volume(global_prob_outlier_mean[g],LogSingleton::getInstance().appendDir("global_prob_outlier"+num2str(g+1)));

	    prob_outlier_mean[g].set_intent(NIFTI_INTENT_ESTIMATE, 0, 0, 0);
	    prob_outlier_mean[g].setDisplayMaximumMinimum(prob_outlier_mean[g].max(),prob_outlier_mean[g].min());
	    save_volume4D(prob_outlier_mean[g],LogSingleton::getInstance().appendDir("prob_outlier"+num2str(g+1)));
	   
	  }
      }

  }

  void Gsmanager::ols()
  {
    Tracer_Plus trace("Gsmanager::ols");

    Matrix dm = design.getdm();
    Matrix pinvdm = (dm.t()*dm).i();
    Matrix pinvdmdm = pinvdm*dm.t();

    if(nevs >= ntpts)
      {
	throw Exception("Singular design. Number of EVs > number of time points. ");
      }
 
    int vox2=0;
    int voxout=0;
    for(int x = 0; x < xsize; x++)
      for(int y = 0; y < ysize; y++)
	for(int z = 0; z < zsize; z++)
	  {
	    if(design.getmask()(x,y,z))
	      {

		vox2++;
		if(vox2 > voxout*nmaskvoxels/100.0)
		  {
		    //cout<<(voxout+1)<<'%' <<'\r';		
		    cout << " " << (voxout+1);
		    cout.flush();
		    voxout++;
		  }

		if(design.is_voxelwise_dm())
		  {
		    dm = design.getdm(x,y,z);
		    pinvdm = (dm.t()*dm).i();
		    pinvdmdm = pinvdm*dm.t();
// 		    OUT(x);OUT(y);OUT(z);
// 		    OUT(dm);
		  }

		ColumnVector Y = design.getcopedata().voxelts(x,y,z);

		ColumnVector petmp = pinvdmdm*Y;

		ColumnVector r = Y-dm*petmp;
		float r2 = (r.t()*r).AsScalar();
		SymmetricMatrix covariance;
		covariance << (pinvdm*r2/(ntpts-nevs));

		// insert any zero EV PEs back in
		petmp = design.insert_zeroev_pes(x,y,z,petmp);

		if(!opts.no_pe_output.value())
		  {
		    // store pe
		    for(int e = 0; e < nevs; e++)
		      {
			pes[e](x,y,z) = petmp(e+1);
		      }
		  }

		// insert any zero EV PEs back in
		covariance = design.insert_zeroev_covpes(x,y,z,covariance);
			
// 		if(!opts.no_pe_output.value())
// 		  {
// 		    // store cov
// 		    ColumnVector covariance2;
// 		    reshape(covariance2, covariance, nevs*nevs, 1);
// 		    cov_pes.setvoxelts(covariance2,x,y,z);
// 		  }

		ols_contrasts(petmp,covariance,x,y,z);
	      }
	  }

    cout << endl;

  }

  void Gsmanager::fixed_effects_onvoxel(const ColumnVector& Y, const Matrix& z, const ColumnVector& S, ColumnVector& gam, SymmetricMatrix& gamcovariance)
  {
 
    Tracer_Plus trace("Gsmanager::fixed_effects_onvoxel");

    // calc gam
    DiagonalMatrix iU(ntpts);
    iU = 0;
   
    for(int t=1;t<=ntpts;t++)
      {		    
	iU(t) = 1.0/S(t);
      }
   
    // calc gam covariance
    gamcovariance << (z.t()*iU*z).i();

    gam=gamcovariance*z.t()*iU*Y;          	

  }

  void Gsmanager::fixed_effects()
  {
    Tracer_Plus trace("Gsmanager::fixed_effects");
  
    Matrix dm = design.getdm();    

    // loop through voxels calling fixed effects inference on each
    OUT(nmaskvoxels);
    
    int vox2=0;
    int voxout=0;
    for(int x = 0; x < xsize; x++)
      for(int y = 0; y < ysize; y++)
	for(int z = 0; z < zsize; z++)
	  {
	    if(design.getmask()(x,y,z))
	      {
		vox2++;
		if(vox2 > voxout*nmaskvoxels/100.0)
		  {
		    //cout<<(voxout+1)<<'%' <<'\r';		
		    cout << " " << (voxout+1);
		    cout.flush();
		    voxout++;
		  }

		if(design.is_voxelwise_dm())
		  {
		    dm = design.getdm(x,y,z);
		  }	       

		//	cout << x << "," << y << "," << z << endl;

		// setup data
		ColumnVector Y = design.getcopedata().voxelts(x,y,z);
		ColumnVector S = design.getvarcopedata().voxelts(x,y,z);

		ColumnVector gam;
		SymmetricMatrix gamcovariance;

		fixed_effects_onvoxel(Y, dm, S, gam, gamcovariance);
  
		// insert any zero EV PEs back in
		gam = design.insert_zeroev_pes(x,y,z,gam);		    

		if(!opts.no_pe_output.value())
		  {
		    // store results for gam:
		    for(int e = 0; e < nevs; e++)
		      {
			pes[e](x,y,z) = gam(e+1);
		      }
		  }
		   
		// insert any zero EV PEs back in
		gamcovariance = design.insert_zeroev_covpes(x,y,z,gamcovariance);
		    		
// 		if(!opts.no_pe_output.value())
// 		  {
// 		    // store results for gam covariance
// 		    ColumnVector gamcovariance2;
// 		    reshape(gamcovariance2, gamcovariance, nevs*nevs, 1);		
// 		    cov_pes.setvoxelts(gamcovariance2,x,y,z);
// 		  }
		
		fe_contrasts(gam,gamcovariance,x,y,z);
	      }
	  }
    cout << endl;
  }	

  void Gsmanager::flame_stage1_onvoxel(const vector<ColumnVector>& Yg, const ColumnVector& Y, const vector<Matrix>& zg, const Matrix& z, const vector<ColumnVector>& Sg, const ColumnVector& S, ColumnVector& beta, ColumnVector& gam, SymmetricMatrix& gamcovariance, vector<float>& marg, vector<ColumnVector>& weights_g, int& nparams, int px, int py, int pz)
  {
 
    Tracer_Plus trace("Gsmanager::flame_stage1_onvoxel");

    ////////////////////////////////////////////////////////////
    // runs flame stage 1 on a voxel without outlier inference

    // calc betas and weights
    beta.ReSize(ngs);
    beta=0;
 
   for(int g = 1; g <= ngs; g++)
      {
	beta(g) = solveforbeta(Yg[g-1],zg[g-1],Sg[g-1]);

	weights_g[g-1].ReSize(design.getntptsingroup(g));
	
	for(int t=1;t<=design.getntptsingroup(g);t++)
	  {		 		
	    weights_g[g-1](t)=1.0/(Sg[g-1](t)+beta(g));	
	  }
      }
      
    // calc gam
    DiagonalMatrix iU(ntpts);
    iU = 0;
   
    for(int t=1;t<=ntpts;t++)
      {		    
	//	iU(t) = 1.0/(S(t)+beta(design.getgroup(t)));
// 	OUT(t);
// 	OUT(design.getgroup(t));
// 	OUT(design.getglobalindex(design.getgroup(t),t));

	iU(t)=weights_g[design.getgroup(t)-1](design.getindexingroup(design.getgroup(t),t));
      }

    // calc gam covariance
    gamcovariance << (z.t()*iU*z).i();

//     OUT(S);
//     OUT(beta);
//     OUT(gamcovariance);
//     OUT(Y);
//     OUT(z);
//     OUT(iU);
    gam=gamcovariance*z.t()*iU*Y;    	

    // calc log marg posterior
//     loglik=0;
//     float loglik2=0;
    
    ColumnVector gamtmp=gam;
    gamtmp = design.insert_zeroev_pes(px,py,pz,gamtmp);
    for(int g = 1; g <= ngs; g++)
      {	
	// need to extract gamg from gam
	ColumnVector gamg(zg[g-1].Ncols());
	int ind=1;

	for(int e = 1; e <= nevs; e++)
	  {

	    if(design.get_evs_group(e,px,py,pz)==g)
	      gamg(ind++)=gamtmp(e);
	  }

// 	loglik+=log_likelihood(beta(g), gamg, Yg[g-1], zg[g-1], Sg[g-    

	marg[g-1]=marg_posterior_energy(log(beta(g)),Yg[g-1],zg[g-1],Sg[g-1]);       
   
// 	loglik2+=log_likelihood_outlier(beta(g), 0, gam, Yg[g-1], zg[g-1], Sg[g-1],0.3);
      }

    nparams=ngs+gam.Nrows();
  }

  void Gsmanager::init_flame_stage1_inferoutliers()
  {
    Tracer_Plus trace("Gsmanager::init_flame_stage1_inferoutliers");

    ////////////////////////////////////////////////////////////
    // init flame stage 1 for outlier inference 

    for(int x = 0; x < xsize; x++)
      for(int y = 0; y < ysize; y++)
	for(int z = 0; z < zsize; z++)
	  {
	    if(design.getmask()(x,y,z))
	      {
		// extract gam
		ColumnVector gam(nevs);
		for(int e = 0; e < nevs; e++)
		  {
		    gam(e+1)=pes[e](x,y,z);	
		  }
		// remove any zero EV PEs
		gam = design.remove_zeroev_pes(x,y,z,gam);

		//		cout << x << "," << y << "," << z << endl;
		// setup data
		ColumnVector Y = design.getcopedata().voxelts(x,y,z);
	        float norm=sqrt(var(Y).AsScalar());
		OUT(norm);
		Y=Y/norm;
		ColumnVector S = design.getvarcopedata().voxelts(x,y,z);
		S=S/Sqr(norm);
		// setup data for groups		
		vector<Matrix> zg(ngs);
		vector<ColumnVector> Yg(ngs);
		vector<ColumnVector> Sg(ngs);
		for(int g = 1; g <= ngs; g++)
		  {
		    zg[g-1]=design.getdm(x,y,z,g);
		    Yg[g-1]=design.getcopedata(x,y,z,g);
		    Yg[g-1]=Yg[g-1]/norm;
		    Sg[g-1]=design.getvarcopedata(x,y,z,g);
		    Sg[g-1]=Sg[g-1]/Sqr(norm);
		  }
		Matrix dm=design.getdm(x,y,z);

		ColumnVector gamtmp=gam;
		gamtmp = design.insert_zeroev_pes(x,y,z,gamtmp);
		for(int g = 1; g <= ngs; g++)
		  {
		    // need to extact gamg from gam
		    ColumnVector gamg(zg[g-1].Ncols());
		    int ind=1;
		    for(int e = 1; e <= nevs; e++)
		      {
			if(design.get_evs_group(e,x,y,z)==g)
			  gamg(ind++)=gamtmp(e);
		      }
		    
		    // calc residuals
		    ColumnVector resids=Yg[g-1]-zg[g-1]*gamg;
		    // normalise
		    resids=abs(resids-mean(resids).AsScalar());

		    // call k-means on resids
		    vector<int> classif;
		    Matrix data(resids.Nrows(),1);
		    data.Column(1)=resids;
		    Matrix means(1,2);
		    means(1,1)=mean(data).AsScalar();
		    means(1,2)=4*means(1,1);
		    
		    do_kmeans(data,classif,2,means);
		    
		    int count_outliers=0;
		    ColumnVector zout(classif.size());
		    for(unsigned int cock=0; cock<classif.size(); cock++)
		      {
			if(classif[cock]==2) count_outliers++;
			zout(cock+1)=classif[cock];
		      }
		    // 	OUT(zout.t());
		    // 	OUT(means);

		    //	global_prob_outlier(g)=(1.0/design.getntptsingroup(g))/2.0;
		    float gpo=count_outliers/float(classif.size());
		    
		    // 	OUT(gpo);
		    
		    if(gpo>0.25) gpo=0.25;
		    
		    // set prob(outlier) parameters with values found from k-means
		    ColumnVector tmp(design.getntptsingroup(g));
		    tmp=gpo;
		    global_prob_outlier_mean[g-1](x,y,z)=gpo;		   

		    prob_outlier_mean[g-1].setvoxelts(tmp,x,y,z);

		  }
	      }
	  }
  }

//   void Gsmanager::flame_stage1_inferoutliers()
//   {
//     Tracer_Plus trace("Gsmanager::flame_stage1_inferoutliers");
    
//     ////////////////////////////////////////////////////////////
//     // runs flame stage 1  with outlier inference 
//     init_flame_stage1_inferoutliers();

//     int niters=opts.io_niters.value();
//     ColumnVector loglik_io(niters);
//     loglik_io=0;

//     // store no outlier results in case we want to revert to them later		    
//     vector<volume<float> > pes_nooutliers(nevs);
//     vector<volume<float> > beta_mean_nooutliers(ngs);
//     volume4D<float> cov_pes_nooutliers=cov_pes;    
//     for(int e = 0; e < nevs; e++)
//       {
// 	pes_nooutliers[e] = pes[e];
//       }
//     for(int g = 0; g < ngs; g++)
//       {
// 	beta_mean_nooutliers[g]=beta_mean[g];
//       }	  
    
//     ColumnVector num_bins_log_beta_col(niters);
//     ColumnVector num_bins_log_beta_outlier_col(niters);
//     for(int it = 1; it <= niters; it++)
//       {
// 	num_bins_log_beta_col(it)=12+2*(it-1);
// 	num_bins_log_beta_outlier_col(it)=floor(5+1*(it-1));
//       }

// //     num_bins_log_beta_col(niters)=100;
// //     num_bins_log_beta_outlier_col(niters)=30;

//     // iterate
//     for(int it = 1; it <= niters; it++)
//       {

// 	// smooth global prob outlier
// 	if(opts.sigma_smooth_globalproboutlier.value() > 0)
// 	  {
// 	    // smooth global prob(outlier)
// 	    for(int g = 0; g < ngs; g++)
// 	      {
// 		float sig=opts.sigma_smooth_globalproboutlier.value();
// 		volume<float> kern=gaussian_kernel3D(sig, int(4*sig/design.getmask().xdim()));
// 		global_prob_outlier_mean[g]=convolve(global_prob_outlier_mean[g], kern, design.getmask());
// 		//	    global_prob_outlier_mean[g]=smooth2D(global_prob_outlier_mean[g],4);
// 	      }
// 	  }

// 	for(int x = 0; x < xsize; x++)
// 	  for(int y = 0; y < ysize; y++)
// 	    for(int z = 0; z < zsize; z++)
// 	      {
// 		if(design.getmask()(x,y,z))
// 		  {
// 		    /////////// extract params for this voxel
// 		    // extract gam
// 		    ColumnVector gam(nevs);
// 		    for(int e = 0; e < nevs; e++)
// 		      {
// 			gam(e+1)=pes[e](x,y,z);	
// 		      }
// 		    // remove any zero EV PEs
// 		    gam = design.remove_zeroev_pes(x,y,z,gam);

// 		    ColumnVector global_prob_outlier(ngs); // one for each group
// 		    vector<ColumnVector> prob_outlier_g(ngs); // one for each "subject" within each group
// 		    ColumnVector prob_outlier(ntpts); // one for each "subject" (just a restructuring of prob_outlier_g)
// 		    ColumnVector beta_outlier(ngs); // one for each group
// 		    ColumnVector beta(ngs); // one for each group

// 		    for(int g = 0; g < ngs; g++)
// 		      {
// 			beta(g+1)=beta_mean[g](x,y,z);
// 			beta_outlier(g+1)=beta_outlier_mean[g](x,y,z);
// 			global_prob_outlier(g+1)=global_prob_outlier_mean[g](x,y,z);
// 			prob_outlier_g[g]=prob_outlier_mean[g].voxelts(x,y,z);
// 		      }
// 		    ///////////

// 		    /////////// 
// 		    // setup data
// 		    ColumnVector Y = design.getcopedata().voxelts(x,y,z);
// 		    ColumnVector S = design.getvarcopedata().voxelts(x,y,z);
// 		    Matrix dm=design.getdm(x,y,z);

// 		    // setup data for groups		
// 		    vector<Matrix> zg(ngs);
// 		    vector<ColumnVector> Yg(ngs);
// 		    vector<ColumnVector> Sg(ngs);
// 		    for(int g = 1; g <= ngs; g++)
// 		      {
// 			zg[g-1]=design.getdm(x,y,z,g);
// 			Yg[g-1]=design.getcopedata(x,y,z,g);
// 			Sg[g-1]=design.getvarcopedata(x,y,z,g);
// 		      }
		    
// 		    /////////////
// 		    // setup grid for marg
// 		    ColumnVector min_log_beta(ngs);
// 		    min_log_beta=log(1e-9);
// 		    ColumnVector max_log_beta(ngs);
// 		    max_log_beta=log(beta*1e2);
		    
// 		    ColumnVector min_log_beta_outlier(ngs);
// 		    min_log_beta_outlier=log((mean(S)+Minimum(beta))/100.0).AsScalar();
// 		    ColumnVector max_log_beta_outlier(ngs);
// 		    max_log_beta_outlier=log((mean(S)+Minimum(beta))*1e9).AsScalar();
		    
// 		    vector<ColumnVector> log_beta(ngs);
// 		    vector<ColumnVector> log_beta_outlier(ngs);
		    
// 		    // need to extract gamg from gam
// 		    ColumnVector gamtmp=gam;
// 		    gamtmp = design.insert_zeroev_pes(x,y,z,gamtmp);
//  		    vector<ColumnVector> gamg(ngs);
// 		    for(int g = 1; g <= ngs; g++)
// 		      {
// 			gamg[g-1].ReSize(zg[g-1].Ncols());
// 			int ind=1;
// 			for(int e = 1; e <= nevs; e++)
// 			  {
// 			    if(design.get_evs_group(e,x,y,z)==g)
// 			      gamg[g-1](ind++)=gamtmp(e);
// 			  }
// 		      }
		    
// 		    // calc betas for each group
// 		    for(int g = 1; g <= ngs; g++)
// 		      {
			
// 			// setup values for doing 2D grid integration later
// 			int num_bins_log_beta=int(num_bins_log_beta_col(it));

// 			float res_log_beta=(max_log_beta(g)-min_log_beta(g))/(num_bins_log_beta-1);
			
// 			log_beta[g-1].ReSize(num_bins_log_beta);
// 			for(int i=1;i<=log_beta[g-1].Nrows();i++){
// 			  log_beta[g-1](i)=min_log_beta(g)+(i-1)*res_log_beta;
// 			} 
			
// 			// setup values for doing 2D grid integration later		
// 			int num_bins_log_beta_outlier=int(num_bins_log_beta_outlier_col(it));

// 			float res_log_beta_outlier=(max_log_beta_outlier(g)-min_log_beta_outlier(g))/(num_bins_log_beta_outlier-1);
// 			log_beta_outlier[g-1].ReSize(num_bins_log_beta_outlier);
// 			for(int i=1;i<=log_beta_outlier[g-1].Nrows();i++){
// 			  log_beta_outlier[g-1](i)=min_log_beta_outlier(g)+(i)*res_log_beta_outlier;
// 			} 
// 			log_beta_outlier[g-1](log_beta_outlier[g-1].Nrows())=std::log(1e9);

// 			// compute membership given gam, beta, beta_outlier and global_prob_outlier
// 			ColumnVector Uin=Sg[g-1]+beta(g); //variances
// 			ColumnVector Uout=Sg[g-1]+beta(g)+beta_outlier(g); //variances

// 			ColumnVector mu=zg[g-1]*gamg[g-1];

// 			for(int i = 1; i <= mu.Nrows(); i++)
// 			  {
// 			    double pin=(1-global_prob_outlier(g))*normpdf(Yg[g-1](i),mu(i),Uin(i));
// 			    double pout=global_prob_outlier(g)*normpdf(Yg[g-1](i),mu(i),Uout(i));

// 			    if(pin+pout==0)
// 			      {
// 				pin=0.5;
// 				pout=0.5;
// 			      }
			    
// 			    prob_outlier_g[g-1](i)=pout/(pin+pout);

// 			    prob_outlier(design.getglobalindex(g,i))=prob_outlier_g[g-1](i);
// 			  }
			
// 			float tmp=mean(prob_outlier_g[g-1]).AsScalar();

// 			if(isnan(tmp))
// 			  {
// 			    OUT(x); OUT(y); OUT(z);
// 			    OUT(g);
// 			    OUT(gam);
// 			    OUT(gamg[g-1]);
// 			    OUT(zg[g-1]);
// 			    OUT(global_prob_outlier(g));
// 			    OUT(prob_outlier_g[g-1].t());
// 			    OUT(Uin.t());
// 			    OUT(Uout.t());
// 			    OUT(Sg[g-1].t());
// 			    OUT(beta(g));
// 			    OUT(beta_outlier(g));
// 			    OUT(gamg[g-1]);
// 			    OUT(mu.t());
// 			    OUT(Yg[g-1].t());
// 			    for(int i = 1; i <= mu.Nrows(); i++)
// 			      {
// 				OUT(normpdf(Yg[g-1](i),mu(i),Uin(i)));
// 				OUT(normpdf(Yg[g-1](i),mu(i),Uout(i)));
// 			      }
// 			    OUT("Error in Gsmanager::flame_stage1_onvoxel");
// 			    exit(0);
// 			  }
 
// 			// update g with membership
// 			global_prob_outlier(g)=tmp;
			
// 			if(global_prob_outlier(g)>0.25) global_prob_outlier(g)=0.25; 
			
// 			// marg out gam and estimate beta and beta_outlier from numerical integration on a 2D grid
// 			Matrix margpost(log_beta[g-1].Nrows(),log_beta_outlier[g-1].Nrows());
// 			margpost=0;
			
// 			ColumnVector marg_beta(log_beta[g-1].Nrows());
// 			marg_beta=0;
// 			ColumnVector marg_beta_outlier(log_beta_outlier[g-1].Nrows());
// 			marg_beta_outlier=0;
// 			for(int j=1;j<=log_beta_outlier[g-1].Nrows();j++){
// 			  for(int i=1;i<=log_beta[g-1].Nrows();i++){
			    
// 			    margpost(i,j) = marg_posterior_energy_outlier(log_beta[g-1](i),log_beta_outlier[g-1](j),Yg[g-1],zg[g-1],Sg[g-1],prob_outlier_g[g-1]);
// 			    marg_beta_outlier(j)+=margpost(i,j);
// 			    marg_beta(i)+=margpost(i,j);
// 			  }	
// 			}     
			
// 			int ind;
			
// 			marg_beta.Minimum1(ind);
// 			beta(g)=std::exp(log_beta[g-1](ind));
			
// 			marg_beta_outlier.Minimum1(ind);
// 			OUT(ind);
// 			beta_outlier(g)=std::exp(log_beta_outlier[g-1](ind));
			
// 		      }  
		    
// 		    // calc gam
// 		    DiagonalMatrix iU(ntpts);
// 		    iU = 0;
		    
// 		    for(int t=1;t<=ntpts;t++)
// 		      {		   
// 			float vr=Sqr(1-prob_outlier(t))*(S(t)+beta(design.getgroup(t)))+Sqr(prob_outlier(t))*(S(t)+beta(design.getgroup(t))+beta_outlier(design.getgroup(t)));
// 			iU(t) = 1.0/vr;
// 		      }
		    
// 		    // calc gam covariance
// 		    SymmetricMatrix gamcovariance;
// 		    gamcovariance << (dm.t()*iU*dm).i();		    
// 		    gam=gamcovariance*dm.t()*iU*Y;		   
		    
// 		    // calc likelihood
// 		    // first need to extract gamg from gam
// 		    gamtmp=gam;
// 		    gamtmp = design.insert_zeroev_pes(x,y,z,gamtmp);
//  		    for(int g = 1; g <= ngs; g++)
// 		      {
// 			gamg[g-1].ReSize(zg[g-1].Ncols());
// 			int ind=1;
// 			for(int e = 1; e <= nevs; e++)
// 			  {
// 			    if(design.get_evs_group(e,x,y,z)==g)
// 			      gamg[g-1](ind++)=gamtmp(e);
// 			  }
// 		      }
// 		    for(int g = 1; g <= ngs; g++)
// 		      {
// 			loglik_io(it)+=log_likelihood_outlier(beta(g),beta_outlier(g),gamg[g-1],Yg[g-1],zg[g-1],Sg[g-1],global_prob_outlier(g));	
// 		      }
		    		    
// 		    // insert any zero EV PEs back in
// 		    gam = design.insert_zeroev_pes(x,y,z,gam);	      
		
// 		    // store latest results for this voxel
// 		    for(int e = 0; e < nevs; e++)
// 		      {
// 			pes[e](x,y,z)=gam(e+1);	
// 		      }
// 		    for(int g = 0; g < ngs; g++)
// 		      {
// 			beta_mean[g](x,y,z) = beta(g+1);
// 			beta_outlier_mean[g](x,y,z) = beta_outlier(g+1);
// 			global_prob_outlier_mean[g](x,y,z) = global_prob_outlier(g+1);
// 			prob_outlier_mean[g].setvoxelts(prob_outlier_g[g],x,y,z);
// 		      }

// 		    // insert any zero EV PEs back in
// 		    gamcovariance = design.insert_zeroev_covpes(x,y,z,gamcovariance);
	
// 		    // store results for gam covariance
// 		    ColumnVector gamcovariance2;
// 		    reshape(gamcovariance2, gamcovariance, nevs*nevs, 1);		
// 		    cov_pes.setvoxelts(gamcovariance2,x,y,z);		      

// 		  }
// 	      }
//       } // end iterations

//     // write out likelihood history
//     write_ascii_matrix(LogSingleton::getInstance().appendDir("loglik_io"),loglik_io);

//     // Loop through seeing if there were any outliers and calc contrasts if there was
//     for(int x = 0; x < xsize; x++)
//       for(int y = 0; y < ysize; y++)
// 	for(int z = 0; z < zsize; z++)
// 	  {
// 	    if(design.getmask()(x,y,z))
// 	      {	
// 		bool no_outliers=true;
// 		for(int g = 0; g < ngs; g++)
// 		  {		   	
// 		    if(global_prob_outlier_mean[g](x,y,z)>(1.0/ntpts)/10.0)
// 		      no_outliers=false;
// 		  }

// 		if(!no_outliers)
// 		  {
// 		    /////////// extract params for this voxel
// 		    // extract gam
// 		    ColumnVector gam(nevs);
// 		    for(int e = 0; e < nevs; e++)
// 		      {
// 			gam(e+1)=pes[e](x,y,z);	
// 		      }

// 		    // extract gamcovariance
// 		    SymmetricMatrix gamcovariance;
// 		    Matrix tmp_mat;
// 		    ColumnVector tmp=cov_pes.voxelts(x,y,z);
// 		    reshape(tmp_mat, tmp, nevs, nevs);		
// 		    gamcovariance << tmp_mat;

// 		    //////////
// 		    // calc contrasts
// 		    flame1_contrasts_with_outliers(gam,gamcovariance,x,y,z);

// 		  }
// 		else
// 		  {
// 		    // restore no outlier results
// 		    for(int e = 0; e < nevs; e++)
// 		      {
// 			pes[e](x,y,z) = pes_nooutliers[e](x,y,z);
// 		      }
// 		    cov_pes.setvoxelts(cov_pes_nooutliers.voxelts(x,y,z),x,y,z);
		    
// 		    for(int g = 0; g < ngs; g++)
// 		      {
// 			beta_mean[g](x,y,z) = beta_mean_nooutliers[g](x,y,z);
// 		      }	      
// 		  }

// 		for(int g = 1; g <= ngs; g++)
// 		  for(int i = 1; i <= prob_outlier_mean[g-1].tsize(); i++)
// 		    {
// 		      if(prob_outlier_mean[g-1](x,y,z,i-1) < 1e-4) prob_outlier_mean[g-1](x,y,z,i-1)=0;		    
// 		    }
// 	      }
// 	  }			    
//   }

  void Gsmanager::flame_stage1_onvoxel_inferoutliers(const vector<ColumnVector>& Yg, const ColumnVector& Y, const vector<Matrix>& zg, const Matrix& z, const vector<ColumnVector>& Sg, const ColumnVector& S, ColumnVector& beta, ColumnVector& gam, SymmetricMatrix& gamcovariance, ColumnVector& global_prob_outlier, vector<ColumnVector>& prob_outlier_g,  ColumnVector& prob_outlier, ColumnVector& beta_outlier, vector<float>& marg, vector<ColumnVector>& weights_g, int& nparams, vector<bool>& no_outliers, int px, int py, int pz)
  {
    Tracer_Plus trace("Gsmanager::flame_stage1_onvoxel_inferoutliers");

    ////////////////////////////////////////////////////////////
    // runs flame stage 1 on a voxel with outlier inference 

    global_prob_outlier.ReSize(ngs);
    prob_outlier_g.resize(ngs);

    prob_outlier.ReSize(ntpts);
    prob_outlier=0.05;

    beta_outlier.ReSize(ngs);
    beta_outlier=Sqr(100);       

//     OUT("===================");

 //    OUT(beta);
//     OUT(gam);
//     OUT(prob_outlier);

    // setup grid for marg
    ColumnVector min_log_beta(ngs);
    min_log_beta=log(beta*1e-6);
    ColumnVector max_log_beta(ngs);
    max_log_beta=log(beta*1e2);

    ColumnVector min_log_rho_outlier(ngs); // beta_outlier=rho*beta
    min_log_rho_outlier=log(100);
    ColumnVector max_log_rho_outlier(ngs);
    max_log_rho_outlier=log(1e9);

    vector<ColumnVector> log_beta(ngs);
    vector<ColumnVector> log_rho_outlier(ngs);
    
    ColumnVector gamtmp=gam;
    ColumnVector maxsg(ngs);

    gamtmp = design.insert_zeroev_pes(px,py,pz,gamtmp);

 //    OUT(gamtmp);

    for(int g = 1; g <= ngs; g++)
      {
	// calc max gs:
	maxsg(g)=100*Maximum(Sg[g-1]);

	// 	OUT(maxsg);

	// need to extact gamg from gam
	ColumnVector gamg(zg[g-1].Ncols());
	int ind=1;
	for(int e = 1; e <= nevs; e++)
	  {
	    if(design.get_evs_group(e,px,py,pz)==g)
	      gamg(ind++)=gamtmp(e);
	  }

	// calc residuals
	ColumnVector resids=Yg[g-1]-zg[g-1]*gamg;
	// normalise
	resids=abs(resids-mean(resids).AsScalar());

	//////////////////
	// call k-means on resids
	vector<int> zc; // classification labels
	Matrix data(resids.Nrows(),1);
	data.Column(1)=resids;
	Matrix means(1,2);
	means(1,1)=mean(data).AsScalar();
	means(1,2)=4*means(1,1);

// 	OUT(data.t());

	do_kmeans(data,zc,2,means);

	int count_outliers=0;
	ColumnVector zout(zc.size());
	prob_outlier_g[g-1].ReSize(design.getntptsingroup(g));
	prob_outlier_g[g-1]=0.2;
	for(unsigned int cock=0; cock<zc.size(); cock++)
	  {
	    if(zc[cock]==2) count_outliers++;
	    zout(cock+1)=zc[cock];

	    if(zout(cock+1)==2)
	      {
		prob_outlier_g[g-1](cock+1)=0.8;
	      }
	  }

//  	OUT(zout.t());
// 	OUT(means);

	//	global_prob_outlier(g)=(1.0/design.getntptsingroup(g))/2.0;
	float gpo=count_outliers/float(zc.size());

	// 	OUT(gpo);
	
	if(gpo>0.25) 
	  {
	    gpo=0.25;

	    global_prob_outlier(g)=gpo;
	    prob_outlier_g[g-1]=global_prob_outlier(g);
	  }
	else
	  {
	    global_prob_outlier(g)=gpo;
	  }

// 	OUT(prob_outlier_g[g-1]);    
      }    

    int niters=opts.io_niters.value();
    vector<bool> converged(ngs,false);
    vector<float> marg_io_old(ngs,0);

    float tol=1e-3;

    ColumnVector num_bins_log_beta_col(niters);
    ColumnVector num_bins_log_rho_outlier_col(niters);
    for(int it = 1; it <= niters; it++)
      {
//    	num_bins_log_beta_col(it)=10+2*(it-1);
//    	num_bins_log_rho_outlier_col(it)=5+1*(it-1);
   	num_bins_log_beta_col(it)=20+2*(it-1);
   	num_bins_log_rho_outlier_col(it)=20+1*(it-1);
      }

    num_bins_log_beta_col(niters)=100;
    num_bins_log_rho_outlier_col(niters)=30;

    for(int it = 1; it <= niters; it++)
      {
	// need to extact gamg from gam	
	vector<ColumnVector> gamg(ngs);
	ColumnVector gamtmp=gam;
	gamtmp = design.insert_zeroev_pes(px,py,pz,gamtmp);
 	for(int g = 1; g <= ngs; g++)
	  {
	    gamg[g-1].ReSize(zg[g-1].Ncols());
	    int ind=1;
	    for(int e = 1; e <= nevs; e++)
	      {
		if(design.get_evs_group(e,px,py,pz)==g)
		  gamg[g-1](ind++)=gamtmp(e);
	      }
	  }

	// calc betas for each group
	for(int g = 1; g <= ngs; g++)
	  {
	     if(!converged[g-1])
	      {
	    // setup values for doing 2D grid integration later
	    //	    int num_bins_log_beta=30;
	    int num_bins_log_beta=int(num_bins_log_beta_col(it));

	    float res_log_beta=(max_log_beta(g)-min_log_beta(g))/float(num_bins_log_beta-1);

	    log_beta[g-1].ReSize(num_bins_log_beta);
	    for(int i=1;i<=log_beta[g-1].Nrows();i++){
	      log_beta[g-1](i)=min_log_beta(g)+(i-1)*res_log_beta;
	    } 

	    // setup values for doing 2D grid integration later		
	    int num_bins_log_rho_outlier=int(num_bins_log_rho_outlier_col(it));

	    float res_log_rho_outlier=(max_log_rho_outlier(g)-min_log_rho_outlier(g))/float(num_bins_log_rho_outlier-1);
	    log_rho_outlier[g-1].ReSize(num_bins_log_rho_outlier);
	    for(int i=1;i<=log_rho_outlier[g-1].Nrows();i++){
	      log_rho_outlier[g-1](i)=min_log_rho_outlier(g)+(i-1)*res_log_rho_outlier;
	    } 

// 	    OUT(res_log_rho_outlier);
// 	    OUT(exp(min_log_rho_outlier(g)));
// 	    OUT(exp(max_log_rho_outlier(g)));
// 	    OUT(exp(log_rho_outlier[g-1]));

	    // compute membership given gam, beta, beta_outlier and global_prob_outlier
 	    ColumnVector Uin=Sg[g-1]+beta(g); //variances
	    // 	    ColumnVector Uout=Sg[g-1]+beta_outlier(g); //variances
 	    ColumnVector Uout(Uin.Nrows()); Uout=beta_outlier(g); //variances

// 	    OUT(g)
//  	    OUT(gam);
//  	    OUT(gamg[g-1]);
//  	    OUT(zg[g-1]);

	    ColumnVector mu=zg[g-1]*gamg[g-1];

	    for(int i = 1; i <= mu.Nrows(); i++)
	      {
		double pin=(1-global_prob_outlier(g))*normpdf(Yg[g-1](i),mu(i),Uin(i));
		double pout=global_prob_outlier(g)*normpdf(Yg[g-1](i),mu(i),Uout(i));

	// 	if(it==niters){
// 		  OUT(Yg[g-1](i));
// 		  OUT(mu(i));
// 		  OUT(Yg[g-1](i)-mu(i));	
// 		  OUT(Sg[g-1](i));
// 		  OUT(beta(g));		  
// 		  OUT(beta_outlier(g));
  
// 		  OUT(pin);
// 		  OUT(pout);
// 		  OUT(Uout(i));
// 		  OUT(Uin(i));
// 		  OUT(normpdf(Yg[g-1](i),mu(i),Uout(i)));
// 		  OUT(normpdf(Yg[g-1](i),mu(i),Uin(i)));
// // 		}

		if(pin+pout==0)
		  {
		    pin=0.5;
		    pout=0.5;
		  }

 		prob_outlier_g[g-1](i)=pout/(pin+pout);	

		// this will encourage a sparse solution
		if(prob_outlier_g[g-1](i)<0.1)
		  prob_outlier_g[g-1](i)=0;

		prob_outlier(design.getglobalindex(g,i))=prob_outlier_g[g-1](i);
 
	      }
	
// 	    prob_outlier_g[g-1]=0;
// 	    prob_outlier_g[g-1](1)=1;
// 	    prob_outlier=0;
// 	    prob_outlier(1)=1;
// 	    global_prob_outlier(g)=0.05;

	    float tmp=mean(prob_outlier_g[g-1]).AsScalar();
   	    if(isnan(tmp))
	      {
		OUT(px); OUT(py); OUT(pz);
		OUT(g);
		OUT(gam);
		OUT(gamg[g-1]);
		OUT(zg[g-1]);
		OUT(global_prob_outlier(g));
		OUT(prob_outlier_g[g-1].t());
		OUT(Uin.t());
		OUT(Uout.t());
		OUT(Sg[g-1].t());
		OUT(beta(g));
		OUT(beta_outlier(g));
		OUT(gamg[g-1]);
		OUT(mu.t());
		OUT(Yg[g-1].t());
		for(int i = 1; i <= mu.Nrows(); i++)
		  {
		    OUT(normpdf(Yg[g-1](i),mu(i),Uin(i)));
		    OUT(normpdf(Yg[g-1](i),mu(i),Uout(i)));
		  }
		OUT("Error in Gsmanager::flame_stage1_onvoxel_inferoutliers");
		exit(0);
	      }
 
	    // update g with membership
	    global_prob_outlier(g)=tmp;

	    if(global_prob_outlier(g)>0.25) global_prob_outlier(g)=0.25; 
	    
	    // marg out gam and estimate beta and beta_outlier from numerical integration on a 2D grid
	    Matrix margpost(log_beta[g-1].Nrows(),log_rho_outlier[g-1].Nrows());
	    margpost=0;

	    ColumnVector marg_beta(log_beta[g-1].Nrows());
	    marg_beta=0;
	    ColumnVector marg_beta_outlier(log_rho_outlier[g-1].Nrows());
	    marg_beta_outlier=0;
	    for(int j=1;j<=log_rho_outlier[g-1].Nrows();j++){
	      for(int i=1;i<=log_beta[g-1].Nrows();i++){
		
 		float log_beta_outlier_tmp=std::log(std::exp(log_rho_outlier[g-1](j))*(std::exp(log_beta[g-1](i))+maxsg(g)));
		//float log_beta_outlier_tmp=std::log(std::exp(log_rho_outlier[g-1](j))*(std::exp(log_beta[g-1](i))));
		margpost(i,j) = marg_posterior_energy_outlier(log_beta[g-1](i),log_beta_outlier_tmp,Yg[g-1],zg[g-1],Sg[g-1],prob_outlier_g[g-1]);
		marg_beta_outlier(j)+=margpost(i,j);
		marg_beta(i)+=margpost(i,j);
	      }	
	    }     

// 	    if(it==1)
// 	      {
// 		OUT(prob_outlier_g[g-1]);
// 	      }

// 	    if(it==niters)
// 	      {
// 		OUT(it);
// 		write_ascii_matrix(log_rho_outlier[0],LogSingleton::getInstance().appendDir("log_rho_outlier"));
// 		write_ascii_matrix(log_beta[0],LogSingleton::getInstance().appendDir("log_beta"));
// 		write_ascii_matrix(Sg[0],LogSingleton::getInstance().appendDir("Sg"));
// 		write_ascii_matrix(Yg[0],LogSingleton::getInstance().appendDir("Yg"));
// 		write_ascii_matrix(zg[0],LogSingleton::getInstance().appendDir("zg"));
// 		write_ascii_matrix(prob_outlier_g[0],LogSingleton::getInstance().appendDir("prob_outlier_g"));
// 		write_ascii_matrix(margpost,LogSingleton::getInstance().appendDir("margpost"));
// 		write_ascii_matrix(marg_beta,LogSingleton::getInstance().appendDir("marg_beta"));
// 		write_ascii_matrix(marg_beta_outlier,LogSingleton::getInstance().appendDir("marg_beta_outlier"));
		
	
// 		//X=load('zg'); y=load('Yg'); s=load('Sg'); x2=X(2:end); y2=y(2:end); s2=s(2:end); res=ols(y2,x2); [gam, beta, covgam] = flame1(y2,x2,s2); gam, covgam, beta, gam/sqrt(covgam)
// 		// p=load('prob_outlier_g');
// 		//mp=load('margpost');mb=load('marg_beta');mbo=load('marg_beta_outlier');lb=load('log_beta');[i,j]=min(mb);exp(lb(j))
// 	      }
	    
	    int ind;	    
	    marg_beta.Minimum1(ind);
	    beta(g)=std::exp(log_beta[g-1](ind));

	    int ind2;
	    marg_beta_outlier.Minimum1(ind2);
	    beta_outlier(g)=std::exp(log_rho_outlier[g-1](ind2))*(std::exp(log_beta[g-1](ind))+maxsg(g));
	    //beta_outlier(g)=std::exp(log_rho_outlier[g-1](ind2))*(std::exp(log_beta[g-1](ind)));	  

	    // 	    if(it==niters)
	    // 	      {
	    // 		OUT(ind);
	    // 		OUT(ind2);
	    // 		OUT(log_rho_outlier[0].Nrows());
	    // 		OUT(std::exp(log_rho_outlier[0](ind2)));
	    // 		OUT(beta_outlier(g));
	    // 		OUT(std::exp(log_beta[0](ind)));
	    // 		OUT(log_beta[0](ind));
	    // 		OUT(marg_posterior_energy_outlier(log(beta(g)),log(beta_outlier(g)),Yg[0],zg[0],Sg[0],prob_outlier_g[0]));
	    // 		OUT(marg_posterior_energy(log(beta(g)),Yg[0],zg[0],Sg[0]));

	    // 	      }

	    weights_g[g-1].ReSize(design.getntptsingroup(g));

	    for(int t=1;t<=design.getntptsingroup(g);t++)
	      {		 		
		weights_g[g-1](t)=1.0/(Sqr(1-prob_outlier_g[g-1](t))*(Sg[g-1](t)+beta(g))+Sqr(prob_outlier_g[g-1](t))*(beta_outlier(g)));		
	      
		if(prob_outlier_g[g-1](t)>0.999)
		  weights_g[g-1](t)=0;	

	      }

	      }
	  }
   	
	// calc gam
	DiagonalMatrix iU(ntpts);
	// 	DiagonalMatrix iU2(ntpts);
	
 	for(int t=1;t<=ntpts;t++)
 	  {		 

	    // 	    OUT(design.getgroup(t));
	    // 	    OUT(design.getglobalindex(design.getgroup(t),t));

	    iU(t)=weights_g[design.getgroup(t)-1](design.getindexingroup(design.getgroup(t),t));

	    // 	    // 	    float vr=Sqr(1-prob_outlier(t))*(S(t)+beta(design.getgroup(t)))+Sqr(prob_outlier(t))*(S(t)+beta_outlier(design.getgroup(t)));
	    // 	    float vr=Sqr(1-prob_outlier(t))*(S(t)+beta(design.getgroup(t)))+Sqr(prob_outlier(t))*(beta_outlier(design.getgroup(t)));
	    
	    // 	    iU(t) = 1.0/vr;
	      
	    // 	      if(prob_outlier(t)>0.999)
	    // 		iU(t)=0;

 	  }

	// calc gam covariance
	gamcovariance << pinv(z.t()*iU*z);

	gam=gamcovariance*z.t()*iU*Y; 
	       
	  // 	if(it==niters)	      	
	  // 	  {
	  //  	    OUT(Y.t());
	  //  	    OUT(z.t());
	  //  	    OUT(S.t());
	  // // 	    OUT(iU);
	  //  	    OUT(gam);
	  // 	    OUT(prob_outlier_g[0].t());
	  // 	    OUT(prob_outlier.t());
	  // 	    OUT(global_prob_outlier(1));
	  // 	    OUT(beta);
	  // 	    OUT(beta_outlier);
	  // 	    OUT(max_log_beta);
	  // 	    OUT(min_log_beta);
	  // // 	    OUT(gamcovariance);
	  // 	    OUT(gam(1)/sqrt(gamcovariance(1,1)));	  	
	  // 	  }

	  // 	gamtmp=gam;
	  // 	gamtmp = design.insert_zeroev_pes(px,py,pz,gamtmp);

	  for(int g = 1; g <= ngs; g++)
 	  {

	    if(it==niters && mean(prob_outlier_g[g-1]).AsScalar() > 0.5)	      	
	      {
		OUT(px); OUT(py); OUT(pz);
		OUT(g);
		OUT(Y.t());
		OUT(z.t());
		OUT(S.t());
		// 	    OUT(iU);
		OUT(gam);
		OUT(prob_outlier_g[g-1].t());
		OUT(prob_outlier.t());
		OUT(global_prob_outlier(g));
		OUT(beta);
		OUT(beta_outlier);
		OUT(max_log_beta);
		OUT(min_log_beta);
		// 	    OUT(gamcovariance);
		OUT(gam(1)/sqrt(gamcovariance(1,1)));	  
		throw Exception("Excessive number of outliers detected. Something has gone wrong.");	
	      }

// 	    // need to extract gamg from gam
// 	    ColumnVector gamg(zg[g-1].Ncols());
// 	    int ind=1;
// 	    for(int e = 1; e <= nevs; e++)
// 	      {
// 		if(design.get_evs_group(e,px,py,pz)==g)
// 		  gamg(ind++)=gamtmp(e);
// 	      }
	    
// 	    loglik_io+=log_likelihood_outlier(beta(g),beta_outlier(g),gamg,Yg[g-1],zg[g-1],Sg[g-1],global_prob_outlier(g),prob_outlier_g[g-1]);	
	
	    // calc log marg posterior to monitor convergence
	   
 	float marg_io=marg_posterior_energy_outlier(log(beta(g)),log(beta_outlier(g)),Yg[g-1],zg[g-1],Sg[g-1],prob_outlier_g[g-1]);	


	if(abs((marg_io - marg_io_old[g-1])/marg_io_old[g-1]) < tol && it>6)
	      {
		converged[g-1]=true;
	       
	      }
	    
	    marg_io_old[g-1]=marg_io;

	  }
    

      } // end iterations

    // check to see if there are any outliers
    no_outliers.resize(ngs);
    for(int g = 0; g < ngs; g++)
      {		   	
	no_outliers[g]=true;
	if(global_prob_outlier(g+1)>(1.0/ntpts)/10.0)
	  no_outliers[g]=false;
      }
 
//    // due to finite variance for outlier class, set prob_outlier to zero if less than 1e-4:
//     for(int g = 1; g <= ngs; g++)
//       for(int i = 1; i <= prob_outlier_g[g-1].Nrows(); i++)
// 	{
// 	  if(prob_outlier_g[g-1](i) < 1e-4) prob_outlier_g[g-1](i)=0;
	  
// 	  prob_outlier(design.getglobalindex(g,i))=prob_outlier_g[g-1](i);
// 	}
    

//     float loglik_io=0;
//     gamtmp=gam;
//     gamtmp = design.insert_zeroev_pes(px,py,pz,gamtmp);
//     for(int g = 1; g <= ngs; g++)
//       {
// 	// need to extact gamg from gam
// 	ColumnVector gamg(zg[g-1].Ncols());
// 	int ind=1;
// 	for(int e = 1; e <= nevs; e++)
// 	  {
// 	    if(design.get_evs_group(e,px,py,pz)==g)
// 	      gamg(ind++)=gamtmp(e);
// 	  }

// 	loglik_io+=log_likelihood_outlier(beta(g),beta_outlier(g),gamg,Yg[g-1],zg[g-1],Sg[g-1],global_prob_outlier(g),prob_outlier_g[g-1]);
//       }
//     int nparams_io=ngs+gam.Nrows();

//     nparams=gam.Nrows();

//  //    OUT(loglik_io);
// //     OUT(loglik);  
// //     OUT(nparams_io);
// //     OUT(nparams);

//     no_outliers=true;
//     for(int g = 0; g < ngs; g++)
//       {		   	
// 	if(global_prob_outlier(g+1)>(1.0/ntpts)/10.0)
// 	  no_outliers=false;
//       }
    
//     float complexity_term=0.0;
//     float complexity_term_io=0.0;
    
//     if(opts.model_select_mode.value()==string("aic"))
//       {
// 	complexity_term=nparams*2;
// 	complexity_term_io=nparams_io*2;	
//       } 
//     else if(opts.model_select_mode.value()==string("loglik"))
//       {
// 	complexity_term=0;
// 	complexity_term_io=0;	
//       }
//     else
//       {
// 	cout<< "Invalid model select mode. Using bic."<< endl;
// 	complexity_term=nparams*std::log(ntpts);
// 	complexity_term_io=nparams_io*std::log(ntpts);
//       }

//     float bic=-2*loglik+complexity_term;
//     float bic_io=-2*loglik_io+complexity_term_io;

//     // check best BIC between with and without outliers
//     if(!no_outliers)
//        {
// 	 //	OUT("Oh no");
// 	 //if(loglik>=loglik_io)
// 	 if(bic<=bic_io)
// 	   {
// 	     // 	     OUT(loglik);
// 	     // 	     OUT(loglik_io);
// 	     OUT("Oh yes");
// 	     no_outliers=true;
// 	   }
//        }
    
    ///////////////////////////

    // double check best energy between with and without outliers
    // calc log marg posterior
    
    for(int g = 1; g <= ngs; g++)
      {	
	// calc log marg posterior 
	float marg_io=marg_posterior_energy_outlier(log(beta(g)),log(beta_outlier(g)),Yg[g-1],zg[g-1],Sg[g-1],prob_outlier_g[g-1]); 
	
	if(!no_outliers[g-1])
	  {	
	    if(marg[g-1]<marg_io)
	      {
		// 	    OUT(marg_io);
		// 	    OUT(marg);

		//		OUT("No outliers");
		// 		OUT(beta(g));
		// 		OUT(beta_outlier(g));
		// 		OUT(g);
		
		no_outliers[g-1]=true;
	      }
	  }
      }
  }
  
  void Gsmanager::flame_stage1()
  {
    Tracer_Plus trace("Gsmanager::flame_stage1");
  
    if(nevs >= ntpts)
      {
	throw Exception("nevs >= ntpts, Singular matrix.");
      }       

    // loop through voxels calling flame stage 1 on each
    OUT(nmaskvoxels);
    int vox2=0;
    int voxout=0;
    for(int x = 0; x < xsize; x++)
      for(int y = 0; y < ysize; y++)
	for(int z = 0; z < zsize; z++)
	  {
	    if(design.getmask()(x,y,z))
	      //           if(x==20 && y==36 && z==0)		      
	      {
		vox2++;
		if(vox2 > voxout*nmaskvoxels/100.0)
		  {
		    //cout<<(voxout+1)<<'%' <<'\r';		
		    cout << " " << (voxout+1);
		    cout.flush();
		    voxout++;
		  }
		//		cout << x << "," << y << "," << z << endl;       
		// setup data
	// 	ColumnVector Y = design.getcopedata().voxelts(x,y,z);
// 		ColumnVector S = design.getvarcopedata().voxelts(x,y,z);
// 		Matrix dm = design.getdm(x,y,z);    
   
// 		// setup data for groups		
// 		vector<Matrix> zg(ngs);
// 		vector<ColumnVector> Yg(ngs);
// 		vector<ColumnVector> Sg(ngs);
// 		for(int g = 1; g <= ngs; g++)
// 		  {
// 		    zg[g-1]=design.getdm(x,y,z,g);
// 		    Yg[g-1]=design.getcopedata(x,y,z,g);
// 		    Sg[g-1]=design.getvarcopedata(x,y,z,g);
// 		  }

 		Matrix dm = design.getdm(x,y,z);    
		ColumnVector Y = design.getcopedata().voxelts(x,y,z);
	        float norm=sqrt(var(Y).AsScalar());
		//		OUT(norm);
		Y=Y/norm;
		ColumnVector S = design.getvarcopedata().voxelts(x,y,z);
		S=S/Sqr(norm);
		// setup data for groups		
		vector<Matrix> zg(ngs);
		vector<ColumnVector> Yg(ngs);
		vector<ColumnVector> Sg(ngs);
		for(int g = 1; g <= ngs; g++)
		  {
		    zg[g-1]=design.getdm(x,y,z,g);
		    Yg[g-1]=design.getcopedata(x,y,z,g);
		    Yg[g-1]=Yg[g-1]/norm;
		    Sg[g-1]=design.getvarcopedata(x,y,z,g);
		    Sg[g-1]=Sg[g-1]/Sqr(norm);
		  }

		// containers for output
		ColumnVector beta;
		ColumnVector gam;
		SymmetricMatrix gamcovariance;

		vector<float> marg(ngs);
		int nparams;
		vector<ColumnVector> weights_g(ngs);

		flame_stage1_onvoxel(Yg, Y, zg, dm, Sg, S, beta, gam, gamcovariance,marg,weights_g,nparams,x,y,z);

		if(opts.infer_outliers.value() && opts.sigma_smooth_globalproboutlier.value()<=0)
		  {
		    vector<bool> no_outliers(ngs);

		    ColumnVector global_prob_outlier; // one for each group
		    vector<ColumnVector> prob_outlier_g; // one for each "subject" within each group
		    ColumnVector prob_outlier; // one for each "subject" (just a restructuring of prob_outlier_g)
		    ColumnVector beta_outlier; // one for each group

		    ColumnVector beta_io=beta;
		    ColumnVector gam_io=gam;
		    SymmetricMatrix gamcovariance_io=gamcovariance;
		    vector<ColumnVector> weights_g_io(ngs);

		    flame_stage1_onvoxel_inferoutliers(Yg, Y, zg, dm, Sg, S, beta_io, gam_io, gamcovariance_io, global_prob_outlier, prob_outlier_g, prob_outlier, beta_outlier,marg, weights_g_io,nparams,no_outliers,x,y,z);
		    
		    // if outliers then replace flame1 results with outlier results		    					

		    // OUT(no_outliers[0]);

		    for(int g = 0; g < ngs; g++)
		      {
			if(!no_outliers[g])
			  {
			    beta(g+1)=beta_io(g+1);
			    weights_g[g] = weights_g_io[g];
			  }
		      }

		    // recalc gam and gamcovariance with blend of outlier and no outlier betas
		    DiagonalMatrix iU(ntpts);		    
		    for(int t=1;t<=ntpts;t++)
		      {		 
			iU(t)=weights_g[design.getgroup(t)-1](design.getindexingroup(design.getgroup(t),t));
		      }
		    gamcovariance << pinv(dm.t()*iU*dm);		    
		    gam=gamcovariance*dm.t()*iU*Y;     

		    // store outlier results		    
		    for(int g = 0; g < ngs; g++)
		      {
			if(!no_outliers[g])
			  {
			    beta_outlier_mean[g](x,y,z) = beta_outlier(g+1)*Sqr(norm);
			    global_prob_outlier_mean[g](x,y,z) = global_prob_outlier(g+1);
			    prob_outlier_mean[g].setvoxelts(prob_outlier_g[g],x,y,z);
			  }
		      }		    		 
		
		  }

		for(int g = 0; g < ngs; g++)
		  {
		    float betamean = beta(g+1);
		    float betavar = 1;
		    beta_b[g](x,y,z) = Sqr(betamean)/betavar;
		    beta_c[g](x,y,z) = betamean/betavar;
		    weights[g].setvoxelts(weights_g[g],x,y,z);
		  }

		// insert any zero EV PEs back in
		gam = design.insert_zeroev_pes(x,y,z,gam)*norm;

		// store results for gam:
		for(int e = 0; e < nevs; e++)
		  {
		    pes[e](x,y,z) = gam(e+1);			
		  }

		// store results for beta:
		for(int g = 0; g < ngs; g++)
		  {
		    beta_mean[g](x,y,z) = beta(g+1)*Sqr(norm);
		  }	      
	
		// insert any zero EV PEs back in
		gamcovariance = design.insert_zeroev_covpes(x,y,z,gamcovariance)*Sqr(norm);

		// store results for gam covariance
		ColumnVector gamcovariance2;
		reshape(gamcovariance2, gamcovariance, nevs*nevs, 1);
		cov_pes.setvoxelts(gamcovariance2,x,y,z);

		if(opts.infer_outliers.value() && opts.sigma_smooth_globalproboutlier.value()<=0)
		  flame1_contrasts_with_outliers(gam,gamcovariance,x,y,z);
		else
		  flame1_contrasts(gam,gamcovariance,x,y,z);	  
	      }
	  }
    cout << endl;

    if(opts.infer_outliers.value() && opts.sigma_smooth_globalproboutlier.value()>0)
      {
	throw Exception("flame_stage1_inferoutliers() unsupported");
	// flame_stage1_inferoutliers();	              
      }

    ////////////////////////////////////////////////////////////////
    // if doing flame stage 2 then find which voxels need processing
    if(opts.zlowerthreshold.value() > 0 && opts.runmode.value()==string("flame12"))
      {
	for(int x = 0; x < xsize; x++)
	  for(int y = 0; y < ysize; y++)
	    for(int z = 0; z < zsize; z++)
	      {
		if(design.getmask()(x,y,z))
		  {		    
		    mcmc_mask(x,y,z) = pass_through_to_mcmc(opts.zlowerthreshold.value(),opts.zupperthreshold.value(),x,y,z);		      
		  }

	      }

	  }

    volume<float> tmp_mask=mcmc_mask;
    tmp_mask.binarise(1e-32);
    nmaskvoxels = int(tmp_mask.sum());
    OUT(nmaskvoxels);
  }

  void Gsmanager::flame_stage2()
  {
    // flame stage 2
    // MCMC
    vector<volume4D<float> > gamma_samples(nevs);
    vector<volume<float> > gamma_naccepted(nevs);
    vector<volume<float> > gamma_nrejected(nevs);
    vector<volume4D<float> > beta_samples(ngs);
    vector<volume<float> > beta_naccepted(ngs);
    vector<volume<float> > beta_nrejected(ngs);
    vector<volume4D<float> > phi_samples(ntpts);
    vector<volume<float> > phi_naccepted(ntpts);
    vector<volume<float> > phi_nrejected(ntpts);
    vector<volume4D<float> > ss_samples(design.getnumfcontrasts()+1);
    
    int nsamples = (opts.njumps.value()-opts.burnin.value())/opts.sampleevery.value();
    cout << "njumps = " << opts.njumps.value() << endl;
    cout << "burnin = " << opts.burnin.value() << endl;
    cout << "sampleevery = " << opts.sampleevery.value() << endl;
    cout << "nsamples = " << nsamples << endl << endl;        

    if(opts.verbose.value())
      {
	for(int e = 0; e < nevs; e++)
	  {	
	    gamma_samples[e].reinitialize(xsize,ysize,zsize,nsamples);
	    gamma_samples[e] = 0;
	    gamma_naccepted[e].reinitialize(xsize,ysize,zsize);
	    gamma_naccepted[e] = 0;
	    gamma_nrejected[e].reinitialize(xsize,ysize,zsize);
	    gamma_nrejected[e] = 0;
	  }

	if(dofpassedin)
	  for(int t = 0; t < ntpts; t++)
	    {
	      phi_samples[t].reinitialize(xsize,ysize,zsize,nsamples);
	      phi_samples[t] = 0;
	      phi_naccepted[t].reinitialize(xsize,ysize,zsize);
	      phi_naccepted[t] = 0;
	      phi_nrejected[t].reinitialize(xsize,ysize,zsize);
	      phi_nrejected[t] = 0;
	    }
	
	for(int g = 0; g < ngs; g++)
	  {
	    beta_samples[g].reinitialize(xsize,ysize,zsize,nsamples);
	    beta_samples[g] = 0;
	    beta_naccepted[g].reinitialize(xsize,ysize,zsize);
	    beta_naccepted[g] = 0;
	    beta_nrejected[g].reinitialize(xsize,ysize,zsize);
	    beta_nrejected[g] = 0;
	  }

	for(int f = 0; f < design.getnumfcontrasts()+1; f++)
	  {
	    ss_samples[f].reinitialize(xsize,ysize,zsize,nsamples);
	    ss_samples[f] = 0;
	  }

      }

    Matrix gamsamples(nevs, nsamples);
    Matrix betasamples(ngs, nsamples);
    Matrix phisamples(ntpts, nsamples);
    ColumnVector likelihood_samples(nsamples);
    vector<ColumnVector> sssamples(design.getnumfcontrasts()+1);

    for(int f = 0; f < design.getnumfcontrasts()+1; f++)
      {
	sssamples[f].ReSize(nsamples);
	sssamples[f] = 0;
      }


    cout << "Metropolis Hasting Sampling" << endl;
    cout << "Number of voxels=" << nmaskvoxels << endl;
    cout << "Percentage done:" << endl;
    int vox2=0;
    int voxout=0;

    for(int x = 0; x < xsize; x++)
      for(int y = 0; y < ysize; y++)
	for(int z = 0; z < zsize; z++)
	  {
	    if(design.getmask()(x,y,z))
	      {		
		if(mcmc_mask(x,y,z))
		  {
		    vox2++;
		    gamsamples = 0;
		    betasamples = 0;
		    phisamples = 0;

		    if(vox2 > voxout*nmaskvoxels/100.0)
		      {
			//cout<<(voxout+1)<<'%' <<'\r';		
			cout << " " << (voxout+1);
			cout.flush();
			voxout++;
		      }

		    if(opts.debuglevel.value()==2)
		      {
			cout << "--------------" << endl;
			cout << "x=" << x << "y=" << y << "z=" << z << endl;
		      }


		    ColumnVector gammamean(nevs);
		    for(int e = 0; e < nevs; e++)
		      {
			gammamean(e+1) = pes[e](x,y,z);
		      }

		    // remove any zero EV PEs
		    gammamean = design.remove_zeroev_pes(x,y,z,gammamean);

		    // extract gamcovariance
		    SymmetricMatrix gamcovariance;
		    Matrix tmp_mat;
		    ColumnVector tmp=cov_pes.voxelts(x,y,z);
		    reshape(tmp_mat, tmp, nevs, nevs);		
		    gamcovariance << tmp_mat;
		    gamcovariance=design.remove_zeroev_covpes(x,y,z,gamcovariance);

		    ColumnVector betab(ngs);
		    ColumnVector betac(ngs);
		    for(int g = 0; g < ngs; g++)
		      {
			betab(g+1) = beta_b[g](x,y,z);
			betac(g+1) = beta_c[g](x,y,z);
		      }

		    srand(opts.seed.value());

		    // get ready the outlier results		    
		    ColumnVector prob_outlier_mcmc(ntpts);
		    vector<float> global_prob_outlier_mcmc(ngs);
		    vector<float> beta_outlier_mcmc(ngs);		    		    

		    if(opts.infer_outliers.value())
		      {
			for(int g = 0; g < ngs; g++)
			  {
			    beta_outlier_mcmc[g]=beta_outlier_mean[g](x,y,z);
			    global_prob_outlier_mcmc[g]=global_prob_outlier_mean[g](x,y,z);
			    for(int t = 1; t <= design.getntptsingroup(g+1); t++)
			      {
				prob_outlier_mcmc(design.getglobalindex(g+1,t))=prob_outlier_mean[g](x,y,z,t-1);
			      }
			  }
		      }

// 		    write_ascii_matrix(LogSingleton::getInstance().appendDir("copedata"),design.getcopedata().voxelts(x,y,z));
// 			    write_ascii_matrix(LogSingleton::getInstance().appendDir("varcopedata"),design.getvarcopedata().voxelts(x,y,z));
// 			    OUT(gammamean.t());
// 			    OUT(gamcovariance);
// 			    OUT(betab);
// 			    OUT(betac);
		    Mcmc_Mh mcmc_mh(design.getcopedata().voxelts(x,y,z), design.getvarcopedata().voxelts(x,y,z),design.getdofvarcopedata().voxelts(x,y,z), design, gammamean, gamcovariance, betab, betac, gamsamples, betasamples, phisamples, likelihood_samples, sssamples, nsamples,x,y,z,prob_outlier_mcmc,global_prob_outlier_mcmc,beta_outlier_mcmc,opts.infer_outliers.value());
		    mcmc_mh.setup();
		    mcmc_mh.run();

		    if(!opts.fixmeanfortfit.value())
		      {
			// get mean gamma from samples
			ColumnVector tmp_gam(gammamean.Nrows());
			for(int e = 0; e < gammamean.Nrows(); e++)
			  {
			    tmp_gam(e+1)=mean(gamsamples.Row(e+1).t()).AsScalar();			    
			  }

			// insert any zero EV PEs
			tmp_gam = design.insert_zeroev_pes(x,y,z,tmp_gam);

			for(int e = 0; e < nevs; e++)
			  {
			    pes[e](x,y,z) = tmp_gam(e+1);
			  }

		      }

		    for(int g = 0; g < ngs; g++)
		      {
			beta_mean[g](x,y,z) = mean(betasamples.Row(g+1).t()).AsScalar();
		      }

		    // insert any zero EV PE samples
		    gamsamples=design.insert_zeroev_pemcmcsamples(x,y,z,gamsamples);
		    
		    //		    write_ascii_matrix(LogSingleton::getInstance().appendDir("gamsamples"),gamsamples);
		    flame2_contrasts(gamsamples,x,y,z);

		    if((std::abs(zts[0](x,y,z)-zflame1lowerts[0](x,y,z))>3))
		      {
			cout << endl << "WARNING: FLAME stage 2 has given an abnormally large difference to stage 1" << endl;
			cout << "x=" << x << ",y=" << y << ",z=" << z << endl;
			OUT(zts[0](x,y,z));
			OUT(zflame1lowerts[0](x,y,z));
			OUT(beta_mean[0](x,y,z));
			OUT(beta_b[0](x,y,z)/beta_c[0](x,y,z));

// 			for(int e = 0; e < nevs; e++)
// 			  {
// 			    write_ascii_matrix(LogSingleton::getInstance().appendDir("gamma_"+num2str(e+1)),gamsamples.Row(e+1).t());
// 			  }
// 			for(int g = 0; g < ngs; g++)
// 			  {
// 			    write_ascii_matrix(LogSingleton::getInstance().appendDir("beta_"+num2str(g+1)),betasamples.Row(g+1).t());
// 			  }
		      }

		    if(opts.verbose.value())
		      {
			for(int e = 0; e < nevs; e++)
			  {
			    gamma_samples[e].setvoxelts(gamsamples.Row(e+1).t(),x,y,z);
			    gamma_naccepted[e](x,y,z) = mcmc_mh.getgamma_naccepted()(e+1);
			    gamma_nrejected[e](x,y,z) = mcmc_mh.getgamma_nrejected()(e+1);
			  }	

			for(int f = 0; f < design.getnumfcontrasts()+1; f++)
			  {
			    ss_samples[f].setvoxelts(sssamples[f],x,y,z);
			  }
		    
			for(int g = 0; g < ngs; g++)
			  {
			    beta_samples[g].setvoxelts(betasamples.Row(g+1).t(),x,y,z);
			    beta_naccepted[g](x,y,z) = mcmc_mh.getbeta_naccepted()(g+1);
			    beta_nrejected[g](x,y,z) = mcmc_mh.getbeta_nrejected()(g+1);
			  }

			if(dofpassedin)
			  for(int t = 0; t < ntpts; t++)
			    {
			      phi_samples[t].setvoxelts(phisamples.Row(t+1).t(),x,y,z);
			      phi_naccepted[t](x,y,z) = mcmc_mh.getphi_naccepted()(t+1);
			      phi_nrejected[t](x,y,z) = mcmc_mh.getphi_nrejected()(t+1);
			    }

		      }
		  }	      	    
	      }
	  }
	  
    if(opts.verbose.value())
      {
	for(int g = 0; g < ngs; g++)
	  {
	    save_volume4D(beta_samples[g], LogSingleton::getInstance().appendDir("beta"+num2str(g+1)+"_samples"));	
	    save_volume(beta_naccepted[g], LogSingleton::getInstance().appendDir("beta"+num2str(g+1)+"_naccepted"));
	    save_volume(beta_nrejected[g], LogSingleton::getInstance().appendDir("beta"+num2str(g+1)+"_nrejected"));	
	  }

	if(dofpassedin)
	  for(int t = 0; t < ntpts; t++)
	    {
	      save_volume4D(phi_samples[t], LogSingleton::getInstance().appendDir("phi"+num2str(t+1)+"_samples"));	
	      save_volume(phi_naccepted[t], LogSingleton::getInstance().appendDir("phi"+num2str(t+1)+"_naccepted"));
	      save_volume(phi_nrejected[t], LogSingleton::getInstance().appendDir("phi"+num2str(t+1)+"_nrejected"));	
	    }
	
	for(int e = 0; e < nevs; e++)
	  {
	    save_volume4D(gamma_samples[e], LogSingleton::getInstance().appendDir("gamma"+num2str(e+1)+"_samples"));
	    save_volume(gamma_naccepted[e], LogSingleton::getInstance().appendDir("gamma"+num2str(e+1)+"_naccepted"));
	    save_volume(gamma_nrejected[e], LogSingleton::getInstance().appendDir("gamma"+num2str(e+1)+"_nrejected"));
	    
	  }
	
	// 	for(int f = 0; f < design.getnumfcontrasts()+1; f++)
	// 	  {
	// 	    save_volume4D(ss_samples[f], LogSingleton::getInstance().appendDir("ss"+num2str(f+1)+"_samples"));
	// 	  }

      }

    regularise_flame2_contrasts();
  
    cout << endl;
  }

  void Gsmanager::regularise_flame2_contrasts()
  {
    Tracer_Plus trace("Gsmanager::regularise_flame2_contrasts");

    float sig=opts.sigma_smooth_flame2_dofs.value();

    if(sig != -1)
      {
	volume<float> kern=gaussian_kernel3D(sig, int(4*sig/design.getmask().xdim()));
	for(int t = 0; t < design.getnumtcontrasts(); t++)
	  {	   	
	    tdofs[t]=convolve(tdofs[t], kern, design.getmask());
	  } 
	
	for(int f = 0; f < design.getnumfcontrasts(); f++)
	  {
	    fdof2s[f]=convolve(fdof2s[f], kern, design.getmask());		    
	  }
      }

    for(int x = 0; x < xsize; x++)
      for(int y = 0; y < ysize; y++)
	for(int z = 0; z < zsize; z++)
	  {
	    if(design.getmask()(x,y,z))
	      {		
		if(mcmc_mask(x,y,z))
		  {
		    for(int t = 0; t < design.getnumtcontrasts(); t++)
		      {
			if(int(tdofs[t](x,y,z))<=0)
			  {
			    cout << "x=" << x << ",y=" << y << ",z=" << z << endl;	
			    OUT(zflame1lowerts[0](x,y,z));
			    OUT(tdofs[t](x,y,z));
			    OUT(ts[t](x,y,z));
			    throw Exception("Error. DOF<=0 in Gsmanager::regularise_flame2_contrasts. ");
			  }
			zts[t](x,y,z) = T2z::getInstance().convert(ts[t](x,y,z),int(tdofs[t](x,y,z)));
		      }

		    for(int f = 0; f < design.getnumfcontrasts(); f++)
		      {
			zfs[f](x,y,z) = F2z::getInstance().convert(fs[f](x,y,z),int(fdof1s[f](x,y,z)),int(fdof2s[f](x,y,z)));
		      }

		  }
	      }
	  }
  }

  bool Gsmanager::pass_through_to_mcmc(float zlowerthresh, float zupperthresh, int px, int py, int pz)
  {
    Tracer_Plus trace("Gsmanager::pass_through_to_mcmc");   
    
    bool ret = false;
    
    // both emupper and emlower need to be jointly above
    // or below the threshold region to avoid the need 
    // for MCMC
    for(int t = 0; !ret && t < design.getnumtcontrasts(); t++)
      {		
	ret = !(((zflame1upperts[t](px,py,pz)) > (zupperthresh) &&
		 (zflame1lowerts[t](px,py,pz)) > (zupperthresh)) ||
		((zflame1upperts[t](px,py,pz)) < (zlowerthresh) &&
		 (zflame1lowerts[t](px,py,pz)) < (zlowerthresh)));
      }

    for(int f = 0; !ret && f < design.getnumfcontrasts(); f++)
      {
	ret = !(((zflame1upperfs[f](px,py,pz)) > (zupperthresh) &&
		 (zflame1lowerfs[f](px,py,pz)) > (zupperthresh)) ||
		((zflame1upperfs[f](px,py,pz)) < (zlowerthresh) &&
		 (zflame1lowerfs[f](px,py,pz)) < (zlowerthresh)));
      }

//     OUT(zupperthresh);
//     OUT(zflame1upperts[0](px,py,pz))
//     OUT(zlowerthresh);
//     OUT(zflame1lowerts[0](px,py,pz))
//     OUT(ret);

    return ret;
  }

  void Gsmanager::ols_contrasts(const ColumnVector& mn, const SymmetricMatrix& covariance, int px, int py, int pz)
  {
    Tracer_Plus trace("Gsmanager::ols_contrasts");    
    
    for(int t = 0; t < design.getnumtcontrasts(); t++)
      {
	tdofs[t](px,py,pz) = ntpts - nevs;

	t_ols_contrast(mn,covariance,design.gettcontrast(t+1),tcopes[t](px,py,pz), tvarcopes[t](px,py,pz), ts[t](px,py,pz), tdofs[t](px,py,pz), zts[t](px,py,pz), px,py,pz);
      }

    for(int f = 0; f < design.getnumfcontrasts(); f++)
      {
	fdof1s[f](px,py,pz) = float(design.getfcontrast(f+1).Nrows());
	fdof2s[f](px,py,pz) = float(ntpts - nevs);

	f_ols_contrast(mn,covariance,design.getfcontrast(f+1),fs[f](px,py,pz), fdof1s[f](px,py,pz), fdof2s[f](px,py,pz), zfs[f](px,py,pz), px,py,pz);
      }    
  }

  void Gsmanager::fe_contrasts(const ColumnVector& mn, const SymmetricMatrix& covariance, int px, int py, int pz)
  {
    Tracer_Plus trace("Gsmanager::fe_contrasts");    

    // contrasts for fixed effects

    float sumdof=1000;
    sumdof=design.getsumdofvarcopedata()(px,py,pz);    
    if(sumdof>1000) sumdof=1000; 

    for(int t = 0; t < design.getnumtcontrasts(); t++)
      {
	if(dofpassedin)
	  tdofs[t](px,py,pz) = sumdof - nevs;
	else
	  tdofs[t](px,py,pz) = 1000;

	t_ols_contrast(mn,covariance,design.gettcontrast(t+1),tcopes[t](px,py,pz), tvarcopes[t](px,py,pz), ts[t](px,py,pz), tdofs[t](px,py,pz), zts[t](px,py,pz), px,py,pz);		
      }


    for(int f = 0; f < design.getnumfcontrasts(); f++)
      {
	fdof1s[f](px,py,pz) = float(design.getfcontrast(f+1).Nrows());
	if(dofpassedin)
	  fdof2s[f](px,py,pz) = sumdof - nevs;
	else
	  fdof2s[f](px,py,pz) = 1000;

	f_ols_contrast(mn,covariance,design.getfcontrast(f+1),fs[f](px,py,pz), fdof1s[f](px,py,pz), fdof2s[f](px,py,pz), zfs[f](px,py,pz), px,py,pz);
      }

  }

  void Gsmanager::flame1_contrasts(const ColumnVector& mn, const SymmetricMatrix& covariance, int px, int py, int pz)
  {
    Tracer_Plus trace("Gsmanager::flame1_contrasts");    
    
    for(int t = 0; t < design.getnumtcontrasts(); t++)
      {
	float tdofupper = 1000;

	// double tdoflower = float(ntpts - nevs);
	// eventually work out effective DOF based on which variance groups interact with the contrast
	// for now just take dof from variance group with fewest members (to be conservative)
	float tdoflower = INT_MAX/10.0;

	for(int g=0; g<ngs; g++)
	  {
	    float tmpdof = design.getntptsingroup(g+1) - design.getnevsingroup(g+1);
	    if(design.is_group_in_tcontrast(g+1, design.gettcontrast(t+1)) && tmpdof>0 && tmpdof<tdoflower)
	      tdoflower=tmpdof;
	  }       

	// call first with highest possible DOF
	t_ols_contrast(mn,covariance, design.gettcontrast(t+1), tcopes[t](px,py,pz), tvarcopes[t](px,py,pz), ts[t](px,py,pz), tdofupper, zflame1upperts[t](px,py,pz), px,py,pz);
	
	// call with OLS DOF
	t_ols_contrast(mn,covariance, design.gettcontrast(t+1), tcopes[t](px,py,pz), tvarcopes[t](px,py,pz), ts[t](px,py,pz), tdoflower, zflame1lowerts[t](px,py,pz), px,py,pz);

	// set z-score to lower z bound	
	tdofs[t](px,py,pz) = tdoflower;
	zts[t](px,py,pz) = zflame1lowerts[t](px,py,pz);
      }
	      
    for(int f = 0; f < design.getnumfcontrasts(); f++)
      {
	fdof1s[f](px,py,pz) = float(design.getfcontrast(f+1).Nrows());
	float fdof2upper = 1000;       

	// double fdof2lower = float(ntpts - nevs);
	// eventually work out effective DOF based on which variance groups interact with the contrast
	// for now just take dof from variance group with fewest members (to be conservative)
	float fdof2lower = INT_MAX/10.0;
	for(int g=0; g<ngs; g++)
	  {
	    float tmpdof = design.getntptsingroup(g+1) - design.getnevsingroup(g+1);
	    if(design.is_group_in_fcontrast(g+1, design.getfcontrast(f+1)) && tmpdof>0 && tmpdof<fdof2lower)
	      fdof2lower=tmpdof;
	  }

	// call first with highest possible DOF
	f_ols_contrast(mn,covariance,design.getfcontrast(f+1),fs[f](px,py,pz), fdof1s[f](px,py,pz), fdof2upper, zflame1upperfs[f](px,py,pz), px,py,pz);

	// call with OLS DOF
	f_ols_contrast(mn,covariance,design.getfcontrast(f+1),fs[f](px,py,pz), fdof1s[f](px,py,pz), fdof2lower, zflame1lowerfs[f](px,py,pz), px,py,pz);

	// set z-score to lower z bound
	fdof2s[f](px,py,pz) = fdof2lower;
	zfs[f](px,py,pz) = zflame1lowerfs[f](px,py,pz);	
      }
  }

  void Gsmanager::flame1_contrasts_with_outliers(const ColumnVector& mn, const SymmetricMatrix& covariance, int px, int py, int pz)
  {
    Tracer_Plus trace("Gsmanager::flame1_contrasts_with_outliers");    
    
    for(int t = 0; t < design.getnumtcontrasts(); t++)
      {
	float tdofupper = 1000;

	// eventually work out effective DOF based on which variance groups interact with the contrast
	// for now just take dof from variance group with fewest members (to be conservative)
	float tdoflower = 1000;
	for(int g=0; g<ngs; g++)
	  {
	    float tmpdof = (1-global_prob_outlier_mean[g](px,py,pz))*design.getntptsingroup(g+1) - design.getnevsingroup(g+1);
	    if(design.is_group_in_tcontrast(g+1, design.gettcontrast(t+1)) && tmpdof>0 && tmpdof<tdoflower)
	      tdoflower=tmpdof;
	  }

//  	OUT(global_prob_outlier_mean[0](px,py,pz));
//  	OUT(design.getntptsingroup(1));
//  	OUT(design.getnevsingroup(1));
// 	OUT(nevs);
// 	OUT(ngs);
//  	OUT(tdoflower);
// 	OUT(beta_mean[0](px,py,pz));
// 	OUT(beta_outlier_mean[0](px,py,pz));


	// call first with highest possible DOF
	t_ols_contrast(mn,covariance, design.gettcontrast(t+1), tcopes[t](px,py,pz), tvarcopes[t](px,py,pz), ts[t](px,py,pz), tdofupper, zflame1upperts[t](px,py,pz), px,py,pz);
       
	// call with OLS DOF
	t_ols_contrast(mn,covariance, design.gettcontrast(t+1), tcopes[t](px,py,pz), tvarcopes[t](px,py,pz), ts[t](px,py,pz), tdoflower, zflame1lowerts[t](px,py,pz), px,py,pz);

	// set z-score to lower z bound	
	tdofs[t](px,py,pz) = tdoflower;
	zts[t](px,py,pz) = zflame1lowerts[t](px,py,pz);

// 	OUT(zts[t](px,py,pz));
// 	OUT(tdoflower);
// 	exit(1);
      }	        

    for(int f = 0; f < design.getnumfcontrasts(); f++)
      {
	fdof1s[f](px,py,pz) = float(design.getfcontrast(f+1).Nrows());
	float fdof2upper = 1000;

	// eventually work out effective DOF based on which variance groups interact with the contrast
	// for now just take dof from variance group with fewest members (to be conservative)
	float fdof2lower = INT_MAX/10.0;

	for(int g=0; g<ngs; g++)
	  {
	    float tmpdof = (1-global_prob_outlier_mean[g](px,py,pz))*design.getntptsingroup(g+1) - design.getnevsingroup(g+1);
	    if(design.is_group_in_fcontrast(g+1, design.getfcontrast(f+1)) && tmpdof>0 && tmpdof<fdof2lower)
	      fdof2lower=tmpdof;

// 	    OUT(nevs);
// 	    OUT(ngs);
// 	    OUT(global_prob_outlier_mean[g](px,py,pz)*design.getntptsingroup(g+1));

// 	    OUT(tmpdof);
// 	    OUT(design.is_group_in_fcontrast(g+1, design.getfcontrast(f+1)));

// 	    OUT(g+1);
// 	    OUT(design.getfcontrast(f+1));
	      
	  }

	// call first with highest possible DOF
	f_ols_contrast(mn,covariance,design.getfcontrast(f+1),fs[f](px,py,pz), fdof1s[f](px,py,pz), fdof2upper, zflame1upperfs[f](px,py,pz), px,py,pz);

	// call with OLS DOF
	f_ols_contrast(mn,covariance,design.getfcontrast(f+1),fs[f](px,py,pz), fdof1s[f](px,py,pz), fdof2lower, zflame1lowerfs[f](px,py,pz), px,py,pz);

	// set z-score to lower z bound
	fdof2s[f](px,py,pz) = fdof2lower;
	zfs[f](px,py,pz) = zflame1lowerfs[f](px,py,pz);	

      }

  }

  void Gsmanager::flame2_contrasts(const Matrix& gamsamples, int px, int py, int pz)
  {
    Tracer_Plus trace("Gsmanager::flame2_contrasts");   
    
    for(int t = 0; t < design.getnumtcontrasts(); t++)
      {
	RowVector tcon = design.gettcontrast(t+1);
	t_mcmc_contrast(gamsamples, tcon, tcopes[t](px,py,pz), tvarcopes[t](px,py,pz), ts[t](px,py,pz), tdofs[t](px,py,pz), zts[t](px,py,pz), px,py,pz);		    
      }

    for(int f = 0; f < design.getnumfcontrasts(); f++)
      {
	f_mcmc_contrast(gamsamples, design.getfcontrast(f+1), fs[f](px,py,pz), fdof1s[f](px,py,pz), fdof2s[f](px,py,pz), zfs[f](px,py,pz), px,py,pz);		    
      }
       
  }

  inline void Gsmanager::t_ols_contrast(const ColumnVector& mn, const SymmetricMatrix& covariance, const RowVector& tcontrast, float& cope, float& varcope, float& t, float& dof, float& z, int px, int py, int pz)
  {
    //    Tracer_Plus trace("Gsmanager::t_ols_contrast");

    if(design.tcontrast_has_zeroevs(px,py,pz,tcontrast))
      {
// 	cout << endl << "WARNING: FLAME has detected a contrast that includes a voxelwise regressor which is all zeros at voxel:" << endl;
// 	cout << "x=" << px << ",y=" << py << ",z=" << pz << endl;	
// 	OUT(tcontrast);
	// set values to something sensible and continue
	t=0;
	z=0;
	cope=0;
	varcope=1e32;
	dof=ntpts;
      }

    varcope = (tcontrast*covariance*tcontrast.t()).AsScalar();

    cope = (tcontrast*mn).AsScalar();

    t = cope/sqrt(varcope);
    z = T2z::getInstance().convert(t,int(floor(dof+0.5)));

  }

  inline void Gsmanager::f_ols_contrast(const ColumnVector& mn, const SymmetricMatrix& covariance, const Matrix& fcontrast, float& f, float& dof1, float& dof2, float& z, int px, int py, int pz)
  {
    //    Tracer_Plus trace("Gsmanager::f_ols_contrast");
    if(design.fcontrast_has_zeroevs(px,py,pz,fcontrast))
      {
// 	cout << endl << "WARNING: FLAME has detected a contrast that includes a voxelwise regressor which is all zeros at voxel:" << endl;
// 	cout << "x=" << px << ",y=" << py << ",z=" << pz << endl;	
// 	OUT(fcontrast);
	// set values to something sensible and continue
	f=0;
	z=0;
	dof2=ntpts;
      }

    f = (mn.t()*fcontrast.t()*(fcontrast*covariance*fcontrast.t()).i()*fcontrast*mn/dof1).AsScalar();
    z = F2z::getInstance().convert(f,int(floor(dof1+0.5)),int(floor(dof2+0.5)));

//     OUT(dof1);
//     OUT(dof2);
//     OUT(int(floor(dof1+0.5)));
//     OUT(int(floor(dof2+0.5)));
//     OUT(f);
//     OUT(z);
  }
 
  void Gsmanager::t_mcmc_contrast(const Matrix& gamsamples, const RowVector& tcontrast, float& cope, float& varcope, float& t, float& dof, float& z, int px, int py, int pz)
  {
    Tracer_Plus trace("Gsmanager::t_mcmc_contrast");
 
    //gamsamples(nevs, nsamples);
  
    ColumnVector m;
    SymmetricMatrix covar;

    Matrix tcsamples = tcontrast*gamsamples;
   
    float tmpdof;

    ColumnVector gammean(nevs);gammean=0;
    for(int e = 0; e < nevs; e++)
      gammean(e+1) = pes[e](px,py,pz);

    m = tcontrast*gammean;
    multitfit(tcsamples, m, covar, tmpdof, true); 
       
    dof = float(tmpdof);
 
    varcope = covar(1,1);
    cope = m(1);

    t = cope/sqrt(varcope);
   
    if(design.tcontrast_has_zeroevs(px,py,pz,tcontrast))
      {
// 	cout << endl << "WARNING: FLAME stage 2 has detected a contrast that includes a voxelwise regressor which is all zeros at voxel:" << endl;
// 	cout << "x=" << px << ",y=" << py << ",z=" << pz << endl;	
// 	OUT(tcontrast);
	// set values to something sensible and continue
	t=0;
	z=0;
	cope=0;
	varcope=1e32;
	dof=ntpts;
      }
    else if(int(floor(dof+0.5))<=0)
      {
	cout << "x=" << px << ",y=" << py << ",z=" << pz << endl;	
	OUT(cope);
	OUT(varcope);
	OUT(zflame1lowerts[0](px,py,pz));
	OUT(dof);
	OUT(gammean);
	OUT(m);
	OUT(tcontrast);
	OUT(design.getdm());
	OUT(design.getdm(px,py,pz));   

	throw Exception("Error. DOF<=0 in Gsmanager::t_mcmc_contrast. ");
      }
    else
      {
	z = T2z::getInstance().convert(t,int(floor(dof+0.5)));
      }
  }

  void Gsmanager::f_mcmc_contrast(const Matrix& gamsamples, const Matrix& fcontrast, float& f, float& dof1, float& dof2, float& z, int px, int py, int pz)
  {
    Tracer_Plus trace("Gsmanager::f_mcmc_contrast");
    
    //gamsamples(nevs, nsamples);
    
    ColumnVector m;
    SymmetricMatrix covar;
    
    float tmpdof2;
    ColumnVector gammean(nevs);gammean=0;
    for(int e = 0; e < nevs; e++)
      gammean(e+1) = pes[e](px,py,pz);

    m = gammean;
    multitfit(gamsamples, m, covar, tmpdof2, true);

    dof2 = float(tmpdof2);
    dof1 = float(fcontrast.Nrows());

    f = (m.t()*fcontrast.t()*(fcontrast*covar*fcontrast.t()).i()*fcontrast*m/dof1).AsScalar();

    if(design.fcontrast_has_zeroevs(px,py,pz,fcontrast))
      {
// 	cout << endl << "WARNING: FLAME stage 2 has detected a contrast that includes a voxelwise regressor which is all zeros at voxel:" << endl;
// 	cout << "x=" << px << ",y=" << py << ",z=" << pz << endl;	
// 	OUT(fcontrast);
	// set values to something sensible and continue
	f=0;
	z=0;
	dof2=ntpts;
      }
    else if(int(floor(dof2+0.5))<=0)
      {
	cout << "x=" << px << ",y=" << py << ",z=" << pz << endl;	
	OUT(dof2);
	OUT(gammean);
	OUT(m);
	OUT(fcontrast);
	OUT(design.getdm());
	OUT(design.getdm(px,py,pz));   

	throw Exception("Error. DOF<=0 in Gsmanager::f_mcmc_contrast. ");
      }
    else
      {	       
	z = F2z::getInstance().convert(f,int(floor(dof1+0.5)),int(floor(dof2+0.5)));
      }
  }

}
