/*  mixture_model.cc

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

#include "mixture_model.h"
#include "utils/log.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include "libvis/miscplot.h"
#include "libvis/miscpic.h"
#include "newmat.h"
#include "utils/tracer_plus.h"
#include "mmoptions.h"
#include "newimage/newimagefns.h"
#include "miscmaths/sparsefn.h"
#include <utility>
#include <iomanip>

using namespace NEWMAT;
using namespace MISCPLOT;
using namespace MISCPIC;

void matout(const Matrix& mat, const string& name)
{
  cout << name << "=[";
  cout.setf(ios::scientific);
  cout.width(10);
  cout.precision(10);
  for(int z = 1; z <= mat.Nrows(); z++)
    {
      for(int y = 1; y <= mat.Ncols(); y++)
	{      
	  cout<< mat(z,y);
	  if(y<mat.Ncols())
	    cout << " ";
	}
      if(z<mat.Nrows())
	cout << ";" << endl;
    }
  cout << "]" << endl;
  cout.setf(ios::fixed);
}

// VOXEL ANISOTROPY

namespace Mm {

    string float2str(float f,int width, int prec, bool scientif)
  {
    ostringstream os;
    int redw = int(std::abs(std::log10(std::abs(f))))+1;
    if(width>0)
      os.width(width);
    if(scientif)
      os.setf(ios::scientific);
    os.precision(redw+std::abs(prec));
    os.setf(ios::internal, ios::adjustfield);
    os << f;
    return os.str();
  } 

  bool GammaDistribution::validate() {
    
    if(var<=0)
      return false;

    // want mn to be positive
      if(mn <= minmode)
	{
	  return false;
	}

      // want mode to be positive
      float mode = mn - var/mn;

      if(mode <= minmode)
	{
	  return false;
	}

      return true;
    }

  bool FlippedGammaDistribution::validate() {
        
    if(var<=0)
      return false;

    // want mn to be negative
    if(mn >= -std::abs(minmode))
	  {
	    return false;
	  }
	
	// want mode to be !positive!
	float mode = abs(mn) - var/abs(mn);
	
	if(mode <= std::abs(minmode))
	  {
	    return false;	    
	  }

	return true;
    }

  float SmmVoxelFunction::evaluate(const ColumnVector& wtilde) const
  {
    Tracer_Plus trace("SmmVoxelFunction::evaluate");

    float h = 0;

//      OUT(nclasses);
//     OUT(data);
//      OUT(log_bound);
//     OUT(wtilde.t());

    RowVector w = logistic_transform(wtilde.t(),lambda,log_bound);
//     OUT(w);

    for(int c=1; c <= nclasses; c++)
      {     
	h += w(c)*dists[c-1]->pdf(data);
      }
   
    if(h<=0) {
//       OUT(w);
//       OUT(h);
//      for(int c=1; c <= nclasses; c++)OUT(dists[c-1]->pdf(data));
      return 1e32;
    }
    return -std::log(h);

  }

  SmmFunction::SmmFunction(const ColumnVector& pdata, vector<Distribution*>& pdists, const float& pmrf_precision, const volume<int>& pmask, const vector<Connected_Offset>& pconnected_offsets, const volume<int>& pindices, const SparseMatrix& pD, float plambda, float plog_bound)
      : gEvalFunction(),
	data(pdata),
	dists(pdists),
	mrf_precision(pmrf_precision),
	mask(pmask),
	connected_offsets(pconnected_offsets),
	indices(pindices),
	D(pD),
	num_superthreshold(pdata.Nrows()),
	nclasses(pdists.size()),
	lambda(plambda),
	log_bound(plog_bound)
    {}

  float SmmFunction::evaluate(const ColumnVector& m_tildew) const
  {
    Tracer_Plus trace("SmmFunction::evaluate");

    float ret = 0.0;

    // stuff from MRF term:
    ret = mrf_precision/2.0*quadratic(m_tildew,D);

//     OUT(ret);
//     ret = 0;
//     for(int z = 0; z < mask.zsize(); z++)
//       for(int y = 0; y < mask.ysize(); y++)
// 	for(int x = 0; x < mask.xsize(); x++)
// 	  if(mask(x,y,z))
// 	    {
// 	      int xi=0,yi=0,zi=0;
// 	      for(unsigned int i = 0; i < connected_offsets.size(); i++) 
// 		{
// 		  xi = x+connected_offsets[i].x;
// 		  yi = y+connected_offsets[i].y;
// 		  zi = z+connected_offsets[i].z;
		  		  
// 		  if(mask(xi,yi,zi))
// 		    {
// 		      for(int c = 0; c < nclasses; c++)
// 			{
// 			  ret +=  Sqr(m_tildew(c*num_superthreshold+indices(x,y,z)) - m_tildew(c*num_superthreshold+indices(xi,yi,zi)));
// 			}
// 		    }
// 		}
// 	    }

//     ret *= mrf_precision/4.0;
//     OUT(ret);

    // likelihood terms stuff:
    //    float A = 1.0/std::sqrt(2*M_PI);
    
    for(int r = 1; r<=num_superthreshold; r++)
      {       
	float h = 0;
       
	RowVector wtilde(nclasses);
	wtilde = 0;

	for(int c=0; c < nclasses; c++)
	  {
	    wtilde(c+1) = m_tildew(c*num_superthreshold+r);
	  }

	RowVector w = logistic_transform(wtilde,lambda,log_bound);

	for(int c=1; c <= nclasses; c++)
	  {
	    h += w(c)*dists[c-1]->pdf(data(r));
// 	    OUT(c);
// 	    OUT(dists[c-1]->pdf(data(r)));
	  }
       
	if(h<=0) {
// 	  OUT(w);
// 	  OUT(h);
//	  for(int c=1; c <= nclasses; c++)OUT(dists[c-1]->pdf(data(r)));
	  return 1e32;
	}
	ret += -std::log(h);
      }

    return ret;
  }

  ReturnMatrix SmmFunction::g_evaluate(const ColumnVector& m_tildew) const
  {
    Tracer_Plus trace("SmmFunction::g_evaluate");
        
    //    float A = 1.0/std::sqrt(2*M_PI);        
    ColumnVector derivative_tilde(num_superthreshold*nclasses);
    derivative_tilde = 0;

    multiply(D,m_tildew,derivative_tilde);
    derivative_tilde *= mrf_precision;

//     OUT(derivative_tilde.Rows(1,10).t());
//     derivative_tilde = 0;
//     // first derivative stuff from MRF term:
//     for(int z = 0; z < mask.zsize(); z++)
//       for(int y = 0; y < mask.ysize(); y++)
// 	for(int x = 0; x < mask.xsize(); x++)
// 	  if(mask(x,y,z))
// 	    {
// 	      ColumnVector derivsum(nclasses);
// 	      derivsum = 0;

// 	      int xi=0,yi=0,zi=0;
// 	      for(unsigned int i = 0; i < connected_offsets.size(); i++) 
// 		{
// 		  xi = x+connected_offsets[i].x;
// 		  yi = y+connected_offsets[i].y;
// 		  zi = z+connected_offsets[i].z;
		  		  
// 		  if(mask(xi,yi,zi))
// 		    {
// 		      for(int c = 0; c < nclasses; c++)
// 			{  
// 			  derivsum(c+1) +=  m_tildew(c*num_superthreshold+indices(x,y,z)) - m_tildew(c*num_superthreshold+indices(xi,yi,zi)); 
// 			}
// 		    }
// 		}

// 	      for(int c = 0; c < nclasses; c++)
// 		{		  
// 		  derivative_tilde(c*num_superthreshold+indices(x,y,z)) += derivsum(c+1)*mrf_precision;
// 		}
// 	    }
// OUT(derivative_tilde.Rows(1,10).t());

    // likelihood terms stuff:
    for(int r = 1; r<=num_superthreshold; r++)
      {
	RowVector wtildetmp(nclasses);
	wtildetmp = 0;
	for(int c=1; c <= nclasses; c++)
	  {	
	    // calculate terms for use in calculating hessian and deriv of likelihood:
	    wtildetmp(c) = m_tildew((c-1)*num_superthreshold+r);
	  }

	// LT of y = LT of demean(y)
	    
	const RowVector wtilde = wtildetmp - mean(wtildetmp,2).AsScalar();
		
	RowVector w = logistic_transform(wtilde,lambda,log_bound);
	
	vector<double> R(nclasses,0);
	
	double P = 0;
	float h = 0;
	//float Q = 0;

	for(int c=1; c <= nclasses; c++)
	  {	
	    // calculate terms for use in deriv of likelihood:
	    h += w(c)*dists[c-1]->pdf(data(r));
	    R[c-1] = boundexp(wtilde(c)/(lambda*log_bound));
	    P += R[c-1];
	  }

//  	OUT(r);
//   	OUT(wtilde);
//   	OUT(w);
//  	OUT(P);

	// calculate dw/dy and dy/dx
	vector<ColumnVector> dwdy(nclasses);
	for(int k=1; k <= nclasses; k++)
	  {	
	    dwdy[k-1].ReSize(nclasses);
	    dwdy[k-1] = 0;
//  	    dydx[k-1].ReSize(nclasses);
//  	    dydx[k-1] = 0;
	   
	    for(int c2=1; c2 <= nclasses; c2++)
	      {		
		if(c2==k)
		  {
		    dwdy[k-1](k) = R[k-1]*(1-R[k-1]/P)/(lambda*log_bound*P);
		  }
		else
		  {
		    dwdy[k-1](c2) = -R[k-1]*R[c2-1]/(lambda*log_bound*Sqr(P));
		  }	       
	      }
	  }

	// calculate df/dw
	ColumnVector dfdw(nclasses);
	dfdw = 0;       

	for(int c=1; c <= nclasses; c++)
	  {
	    dfdw(c) = -dists[c-1]->pdf(data(r))/h;
	  }
           
    	// Now fill up derivative vector for tildew (aka x) for this voxel
	ColumnVector derivanal(nclasses);
	derivanal = 0;
	
	for(int k=1; k <= nclasses; k++)
	  {
	    float sum_l = 0;
	    for(int l=1; l <= nclasses; l++)
	      {
		sum_l += dfdw(l)*dwdy[l-1](k);
	      }
	    
	    derivanal(k) = sum_l;
	  }

	// now fill up return vector
	for(int k=1; k <= nclasses; k++)
	  {
	    derivative_tilde((k-1)*num_superthreshold+r) += derivanal(k);
	  }

      }

//     OUT(derivative_tilde.Rows(1,10).t());
    derivative_tilde.Release();
    
    return derivative_tilde;
  }


  SmmFunctionDists::SmmFunctionDists(const ColumnVector& pdata, vector<Distribution*>& pdists, const float& pmrf_precision, const volume<int>& pmask, const vector<Connected_Offset>& pconnected_offsets, const volume<int>& pindices, float plambda, float plog_bound, const ColumnVector& pm_tildew)
    : gEvalFunction(),
      //EvalFunction(),
      data(pdata),
      dists(pdists),
      mrf_precision(pmrf_precision),
      mask(pmask),
      connected_offsets(pconnected_offsets),
      indices(pindices),
      w(pdata.Nrows()),
      num_superthreshold(pdata.Nrows()),
      nclasses(pdists.size()),
      lambda(plambda),
      log_bound(plog_bound),
      m_tildew(pm_tildew)
  {

    for(int r = 1; r<=num_superthreshold; r++)
      {
	RowVector wtilde(nclasses);
	wtilde = 0;
	
	for(int c=0; c < nclasses; c++)
	  {
	    wtilde(c+1) = m_tildew(c*num_superthreshold+r);
	  }
	
	//	OUT(wtilde);
	w[r-1] = logistic_transform(wtilde,lambda,log_bound);
	//OUT(w);
      }

  }

  float SmmFunctionDists::evaluate(const ColumnVector& x) const
  {
    Tracer_Plus trace("SmmFunctionDists::evaluate");

    float ret = 0.0;

    // get dists params from passed in vector
    for(int c = 0; c < nclasses; c++)
      {
	if(!dists[c]->setparams(x(c*2+1),x(c*2+2),1))
	  {
	    return 1e32;
	  }
      }

    // likelihood terms stuff:
    //    float A = 1.0/std::sqrt(2*M_PI);
    
    for(int r = 1; r<=num_superthreshold; r++)
      {       
	float h = 0;
       
	for(int c=1; c <= nclasses; c++)
	  {
	    h += w[r-1](c)*dists[c-1]->pdf(data(r));
	  }


	ret += -std::log(h);
      }

 
    return ret;
  }

  ReturnMatrix SmmFunctionDists::g_evaluate(const ColumnVector& x) const
  {
    Tracer_Plus trace("SmmFunctionDists::g_evaluate");

    // get dists params from passed in vector
    for(int c = 0; c < nclasses; c++)
      {
	dists[c]->setparams(x(c*2+1),x(c*2+2),1);	
      }      

    // likelihood terms stuff:
    ColumnVector meansum(nclasses); 
    meansum = 0.0;
    ColumnVector varsum(nclasses); 
    varsum = 0.0;
    for(int r = 1; r<=num_superthreshold; r++)
      {		
	float h = 0;

	//float Q = 0;

	for(int c=1; c <= nclasses; c++)
	  {	
	    // calculate terms for use in deriv of likelihood:
	    float mult=1;
	    //	    if(dists[c-1]->flipped) mult=-1;
	    h += w[r-1](c)*mult*dists[c-1]->pdf(data(r));	    
	  }	

	for(int c=1; c <= nclasses; c++)
	  {
	    meansum(c) += -w[r-1](c)*dists[c-1]->dpdfdmn(data(r))/h;
	    varsum(c) += -w[r-1](c)*dists[c-1]->dpdfdvar(data(r))/h;
	  }
      }
	
    ColumnVector derivative_tilde(nclasses*2);
    derivative_tilde = 0;
    
    for(int c = 0; c < nclasses; c++)
      {	
	derivative_tilde(c*2+1) = meansum(c+1);
	derivative_tilde(c*2+2) = varsum(c+1);
      }

    
    derivative_tilde.Release();
    
    return derivative_tilde;
  }

  Mixture_Model::Mixture_Model(const volume<float>& pspatial_data, const volume<int>& pmask, const volume<float>& pepi_example_data, float pepibt, vector<Distribution*>& pdists, vector<volume<float> >& pw_means, ColumnVector& pY, MmOptions& popts) :
    xsize(pmask.xsize()),
    ysize(pmask.ysize()),
    zsize(pmask.zsize()),
    num_superthreshold(0),
    nclasses(pdists.size()),
    spatial_data(pspatial_data),
    mask(pmask),
    epi_example_data(pepi_example_data),
    epibt(pepibt),
    localweights(),
    connected_offsets(),
    indices(pmask.xsize(),pmask.ysize(),pmask.zsize()),    
    Y(pY),
    D(),
    m_tildew(),
    mrf_precision(popts.mrfprecstart.value()),
    nonspatial(popts.nonspatial.value()),
    niters(popts.niters.value()),  
    stopearly(),
    updatetheta(popts.updatetheta.value()),
    debuglevel(popts.debuglevel.value()),
    lambda(popts.phi.value()),
    log_bound(10),
    dists(pdists),
    w_means(pw_means),
    mrfprecmultiplier(popts.mrfprecmultiplier.value()),
    initmultiplier(popts.initmultiplier.value()),
    fixmrfprec(popts.fixmrfprec.value())
  {
    if(nonspatial)
      {
	mrf_precision=1e-10;
	fixmrfprec=true;
	updatetheta=false;
	initmultiplier=0.6;
	niters=1;
      }

    connected_offsets.push_back(Connected_Offset(-1,0,0,0,1));
    connected_offsets.push_back(Connected_Offset(1,0,0,1,0));
    connected_offsets.push_back(Connected_Offset(0,-1,0,2,3));
    connected_offsets.push_back(Connected_Offset(0,1,0,3,2));
    connected_offsets.push_back(Connected_Offset(0,0,-1,4,5));
    connected_offsets.push_back(Connected_Offset(0,0,1,5,4));
  }

  Mixture_Model::Mixture_Model(const volume<float>& pspatial_data, const volume<int>& pmask, const volume<float>& pepi_example_data, float pepibt, vector<Distribution*>& pdists, vector<volume<float> >& pw_means, ColumnVector& pY, bool pnonspatial, int pniters, bool pupdatetheta, int pdebuglevel, float pphi, float pmrfprecstart, int pntracesamps, float pmrfprecmultiplier, float pinitmultiplier, bool pfixmrfprec) :
    xsize(pmask.xsize()),
    ysize(pmask.ysize()),
    zsize(pmask.zsize()),
    num_superthreshold(0),
    nclasses(pdists.size()),
    spatial_data(pspatial_data),
    mask(pmask),
    epi_example_data(pepi_example_data),
    epibt(pepibt),
    localweights(),
    connected_offsets(),
    indices(pmask.xsize(),pmask.ysize(),pmask.zsize()),
    Y(pY),
    D(),
    m_tildew(),
    mrf_precision(pmrfprecstart),
    nonspatial(pnonspatial),
    niters(pniters),
    stopearly(),
    updatetheta(pupdatetheta),
    debuglevel(pdebuglevel),
    lambda(pphi),
    log_bound(10.0),
    dists(pdists),
    w_means(pw_means),
    mrfprecmultiplier(pmrfprecmultiplier),
    initmultiplier(pinitmultiplier),
    fixmrfprec(pfixmrfprec)
  {
    connected_offsets.push_back(Connected_Offset(-1,0,0,0,1));
    connected_offsets.push_back(Connected_Offset(1,0,0,1,0));
    connected_offsets.push_back(Connected_Offset(0,-1,0,2,3));
    connected_offsets.push_back(Connected_Offset(0,1,0,3,2));
    connected_offsets.push_back(Connected_Offset(0,0,-1,4,5));
    connected_offsets.push_back(Connected_Offset(0,0,1,5,4));
  }

  void Mixture_Model::run()
  {
    Tracer_Plus trace("Mixture_Model::run");    
    
//  	float mrf_precision_saved = mrf_precision;
	//	mrf_precision = 10;

    save_weights(m_tildew,"_init",false);

	for(it=1; it<=niters; it++)
	  {
	    OUT(it);
	    
	    calculate_taylor_lik();
	    update_voxel_tildew_vb();

// 	    if(it==20) mrf_precision = mrf_precision_saved;

	    if(!fixmrfprec)// && it>20 )
	      {		
		OUT("Calculating trace");
		calculate_trace_tildew_D();

		OUT("Updating MRF precision");
		update_mrf_precision();

		OUT(mrf_precision);
	      }

	    if(updatetheta)	    
	      {
		OUT("Updating distribution params");
		update_theta();
	      }	    	    	  

 //  	    save_weights(volinfo,m_tildew,num2str(it).c_str(),true);

	    cout << "Iterations=" << it << endl;
	  }        
      
  }


  void Mixture_Model::setup()
  {
    Tracer_Plus trace("Mixture_Model::setup");
    
    trace_tol = 0.0001;
    scg_tol = 0.001;  

    if(niters<0)
      {
	stopearly=true;
	niters=120;
      }
    else
      stopearly=false;

    // setup lattice:
    localweights.reinitialize(mask.xsize(),mask.ysize(),mask.zsize(),6);
    localweights = 0;
  
    for(int x = 0; x < mask.xsize(); x++)
      for(int y = 0; y < mask.ysize(); y++)
	for(int z = 0; z < mask.zsize(); z++)
	  if(mask(x,y,z))
	    {
	      num_superthreshold++;
	      
	      int xi,yi,zi;
	      for(unsigned int i = 0; i < connected_offsets.size(); i++) 
		{
		  xi = x+connected_offsets[i].x;
		  yi = y+connected_offsets[i].y;
		  zi = z+connected_offsets[i].z;
		  
		  if(mask(xi,yi,zi))
		    {
		      if(epibt > 0)
			{
			  localweights(x,y,z,connected_offsets[i].ind) = exp(-Sqr(epi_example_data(xi,yi,zi)-epi_example_data(x,y,z))/(2*epibt*epibt));			  
			}
		      else
			{
			  localweights(x,y,z,connected_offsets[i].ind) = 1;
			}
		    }
		}
	    }   
      
    OUT(nclasses);
    OUT(num_superthreshold);
    OUT(niters);
    OUT(mrf_precision);
    OUT(fixmrfprec);
    OUT(mrfprecmultiplier);
    OUT(initmultiplier);
    OUT(updatetheta);
    OUT(nonspatial);
    OUT(lambda);
    OUT(log_bound);

    w_means.resize(nclasses);    
    indices = 0;
    int index=1;
    Y.ReSize(num_superthreshold);
    Y = 0;
    D.ReSize(num_superthreshold*nclasses,num_superthreshold*nclasses);

    m_tildew.ReSize(num_superthreshold*nclasses);
    m_tildew = 0;

    prec_tildew.resize(num_superthreshold);
    cov_tildew.resize(num_superthreshold);
    for(int c = 0; c < nclasses; c++)
      {
	prec_tildew[c].ReSize(nclasses);
	prec_tildew[c] = 0;
	cov_tildew[c].ReSize(nclasses);
	cov_tildew[c] = 0;

      }
    
    precision_lik.ReSize(num_superthreshold*nclasses,num_superthreshold*nclasses);
 
//     derivative_tildew.ReSize(num_superthreshold*nclasses);
//     derivative_tildew = 0;
    derivative_lik.ReSize(num_superthreshold*nclasses);
    derivative_lik = 0;

    ColumnVector num_neigbours(num_superthreshold);
    num_neigbours = 0;

    volume<int> maxcsmap(mask.xsize(),mask.ysize(),mask.zsize());
    maxcsmap=0;

    for(int c = 1; c <= nclasses; c++)
      dists[c-1]->setuseprop(true);
    
    for(int z = 0; z < mask.zsize(); z++)    
      for(int y = 0; y < mask.ysize(); y++)	
	for(int x = 0; x < mask.xsize(); x++)
	  if(mask(x,y,z))
	    {
	      int xi=0,yi=0,zi=0;
	      for(unsigned int i = 0; i < connected_offsets.size(); i++) 
		{
		  xi = x+connected_offsets[i].x;
		  yi = y+connected_offsets[i].y;
		  zi = z+connected_offsets[i].z;
		  
		  if(mask(xi,yi,zi))
		    {
		      num_neigbours(index) += localweights(x,y,z,connected_offsets[i].ind);
		    }
		}

	      indices(x,y,z) = index;
	      Y(index) = spatial_data(x,y,z);

	      ///////////
	      // setup classifications
	      ColumnVector probs(nclasses);
	      probs = 0;

	      //float sum = 0;
	      float tmp = Y(index);

	      int maxc = 1;
	      float maxprob = 0;

	      int minc = 1;
	      float minprob = 1e32;

	      float sumprob = 0;
	
	      for(int c = 1; c <= nclasses; c++)
		{	    
		  probs(c) = dists[c-1]->pdf(tmp);
		  sumprob += probs(c);

		  if(probs(c) > maxprob) {
		    maxc = c;
		    maxprob = probs(c);
		  }	    	    

		  if(probs(c) < minprob) {
		    minc = c;
		    minprob = probs(c);
		  }	    	    
		}

	      maxcsmap(x,y,z)=maxc;
	      ////////////////////

	      index++;
	    }

    save_volume(maxcsmap, LogSingleton::getInstance().appendDir("maxcsmap"));
	
    volume<int> maxcsmapnew;
    maxcsmapnew=maxcsmap;

    for(int z = 0; z < mask.zsize(); z++)    
      for(int y = 0; y < mask.ysize(); y++)	
	for(int x = 0; x < mask.xsize(); x++)
	  if(mask(x,y,z))
	    {
	      // count up classes in neighbourhood
	      RowVector count(nclasses);
	      count=0;
	      int xi=0,yi=0,zi=0;
	      for(unsigned int i = 0; i < connected_offsets.size(); i++) 
		{
		  xi = x+connected_offsets[i].x;
		  yi = y+connected_offsets[i].y;
		  zi = z+connected_offsets[i].z;
		  
		  if(mask(xi,yi,zi))
		    {
		      count(int(maxcsmap(xi,yi,zi)))++;
		    }
		}

	      maxcsmapnew(x,y,z)=maxcsmap(x,y,z);

	      // assign voxel to its own class if number of voxels in 
	      // neighbourhood is gt 1
	      // otherwise assign to class with greatest no of voxels
	      //  if(count(int(maxcsmap(x,y,z)))<2 && !nonspatial)		 
// 	      if(count(int(maxcsmap(x,y,z)))<2 && !nonspatial)	  
// 		{
// 		  int tmpmax=1;
// 		  int tmpmaxindex=1;
// 		  for(int c = 1; c <= nclasses; c++)
// 		    {
// 		      if(count(c)>tmpmax) 
// 			{
// 			  tmpmaxindex=c;
// 			  tmpmax=int(count(c));
// 			}
// 		    }

// 		  maxcsmapnew(x,y,z)=tmpmaxindex;
// 		}

	    }
 
    save_volume(maxcsmapnew, LogSingleton::getInstance().appendDir("maxcsmapnew"));

    // second erosion
    ColumnVector maxcs(num_superthreshold);
    maxcs = 0;
 
    for(int z = 0; z < mask.zsize(); z++)    
      for(int y = 0; y < mask.ysize(); y++)	
	for(int x = 0; x < mask.xsize(); x++)
	  if(mask(x,y,z))
	    {
	      // count up classes in neighbourhood
	      RowVector count(nclasses);
	      count=0;
	      int xi=0,yi=0,zi=0;
	      for(unsigned int i = 0; i < connected_offsets.size(); i++) 
		{
		  xi = x+connected_offsets[i].x;
		  yi = y+connected_offsets[i].y;
		  zi = z+connected_offsets[i].z;
		  
		  if(mask(xi,yi,zi))
		    {
		      count(int(maxcsmapnew(xi,yi,zi)))++;
		    }
		}

	      maxcs(int(indices(x,y,z)))=maxcsmapnew(x,y,z);

	      // assign voxel to its own class if number of voxels in 
	      // neighbourhood is gt 1
	      // otherwise assign to class with greatest no of voxels
	      //  if(count(int(maxcsmap(x,y,z)))<2 && !nonspatial)		 
// 	      if(count(int(maxcsmapnew(x,y,z)))<2 && !nonspatial)	  
// 		{
// 		  int tmpmax=1;
// 		  int tmpmaxindex=1;
// 		  for(int c = 1; c <= nclasses; c++)
// 		    {
// 		      if(count(c)>tmpmax) 
// 			{
// 			  tmpmaxindex=c;
// 			  tmpmax=int(count(c));
// 			}
// 		    }


// 		  // final erosion result
// 		  maxcs(int(indices(x,y,z)))=tmpmaxindex;
		  
// 		}
	    }
 
    for(int z = 0; z < mask.zsize(); z++)
      for(int y = 0; y < mask.ysize(); y++)
	for(int x = 0; x < mask.xsize(); x++)
	  if(mask(x,y,z))
	    {
	      int xi=0,yi=0,zi=0;
	      for(unsigned int i = 0; i < connected_offsets.size(); i++) 
		{
		  xi = x+connected_offsets[i].x;
		  yi = y+connected_offsets[i].y;
		  zi = z+connected_offsets[i].z;
		  
		  if(mask(xi,yi,zi))
		    {
		      for(int c = 0; c < nclasses; c++)
			{
			  // this gives same as sumpairs:
//  			  D.insert(c*num_superthreshold+indices(x,y,z),c*num_superthreshold+indices(xi,yi,zi), -1.0);
//  			  D.insert(c*num_superthreshold+indices(xi,yi,zi),c*num_superthreshold+indices(x,y,z), -1.0);	 

			  float val = std::sqrt(num_neigbours(indices(x,y,z))*num_neigbours(indices(xi,yi,zi)));
			  D.insert(c*num_superthreshold+indices(x,y,z),c*num_superthreshold+indices(xi,yi,zi), -1.0/val);
 			  D.insert(c*num_superthreshold+indices(xi,yi,zi),c*num_superthreshold+indices(x,y,z), -1.0/val);	 
			  
			}
		    }
		}

	      for(int c = 0; c < nclasses; c++)
		//this gives same as sumpairs:
// 		D.insert(c*num_superthreshold+indices(x,y,z),c*num_superthreshold+indices(x,y,z),num_neigbours(indices(x,y,z)));
		
		D.insert(c*num_superthreshold+indices(x,y,z),c*num_superthreshold+indices(x,y,z),1);
	      
	    }
    
    // initialise tildew
    for(int r = 1; r<=num_superthreshold; r++)
      {

	//   	OUT(w);
	// 	for(int c = 1; c <= nclasses; c++)
	// 	  {
	// 	    w(c) = probs(c)/sumprob;
	// 	  }
	//   	OUT(w);
	RowVector w(nclasses);
	w = 0;
	w(int(maxcs(r))) = 1;

	RowVector wtilde = inv_transform(w,lambda,log_bound,initmultiplier);

	if(r<4 || (it>5 && r==-1))
	  {
	    OUT(r);
	    OUT(Y(r));
	    matout(w,"w");
	    matout(wtilde,"wtilde");
	    RowVector w2 = logistic_transform(wtilde,lambda,log_bound);
	    matout(w2,"w");
	  }

	for(int c=0; c < nclasses; c++)
	  {
	    m_tildew(c*num_superthreshold+r) = wtilde(c+1);
	    //m_tildew(c*num_superthreshold+r) = normrnd().AsScalar()*0.1;
	  }

      }

    for(int c = 1; c <= nclasses; c++)
      dists[c-1]->setuseprop(false);
    
  }

  void Mixture_Model::update_theta()
  {
    Tracer_Plus trace("Mixture_Model::update_theta");
   
    // update means and vars using ML estimates, given weights are known
    
    if(true)
      {
	SmmFunctionDists smmfunc(Y,dists,mrf_precision,mask,connected_offsets,indices,lambda,log_bound,m_tildew);
	// set dists params in vector
	ColumnVector x(nclasses*2);
	x=0;    	
	for(int c = 0; c < nclasses; c++)
	  {	    
	    x(c*2+1) = dists[c]->getmean();
	    x(c*2+2) = dists[c]->getvar();
	  }
	
// 	ColumnVector tmp2(x.Nrows());tmp2=1;	
// 	scg(x, smmfunc, tmp2, scg_tol);
// 	//smmfunc.minimize(x,tmp2);

      	float tmp = smmfunc.evaluate(x);
      	OUT(tmp);
	ColumnVector tmp2(x.Nrows());tmp2=1;
	scg(x,smmfunc, tmp2, 0.01);
      	tmp = smmfunc.evaluate(x);
      	OUT(tmp);

	// get dists params from vector
	for(int c = 0; c < nclasses; c++)
	  {
	    dists[c]->setparams(x(c*2+1),x(c*2+2),1);
	  }
      }
    else
      {    
	vector<ColumnVector> weights;
	vector<vector<vector<float> > > weights_samps;
	vector<vector<vector<float> > > tildew_samps;
	//get_weights2(weights, weights_samps, tildew_samps, 100, m_tildew);
	get_weights(weights, m_tildew);
  	
	for(int c = 0; c < nclasses; c++)
	  {	
	    float sumweights = 0;
	    float sumy = 0;

	    for(int r = 1; r<=num_superthreshold; r++)
	      {
		float wtmp = weights[c](r);
		sumweights += wtmp;
		sumy += wtmp*Y(r);	    
	      }
	
	    float mu = sumy/sumweights;
	
	    float ssy = 0;
	    for(int r = 1; r<=num_superthreshold; r++)
	      {
		float wtmp = weights[c](r);	    
		ssy += wtmp*Sqr(Y(r)-mu);
	      }

	    //  	OUT(c);
	    //  	OUT(sumweights);
	    //  	OUT(sumy/sumweights);
	    //  	OUT(ssy/sumweights);

	    if(updatetheta)	
	      dists[c]->setparams(mu,ssy/sumweights,sumweights/num_superthreshold);
	    else
	      dists[c]->setparams(dists[c]->getmean(),dists[c]->getvar(),sumweights/num_superthreshold);
      
	  }    
      }

    meanhist.push_back(dists[0]->getmean());

    OUT(dists[0]->getmean());
  }

  void Mixture_Model::update_tildew_scg()
  {
    Tracer_Plus trace("Mixture_Model::update_tildew_scg");

    OUT("Doing tildew SCG");
    SmmFunction smmfunc(Y,dists,mrf_precision,mask,connected_offsets,indices,D,lambda,log_bound);
    
    float tmp = smmfunc.evaluate(m_tildew);
    OUT(tmp);
    ColumnVector tmp2(m_tildew.Nrows());tmp2=1;    
    scg(m_tildew,smmfunc, tmp2, 0.01);
    tmp = smmfunc.evaluate(m_tildew);
    OUT(tmp); 
  }

  void Mixture_Model::update_voxel_tildew_vb()
  {
    Tracer_Plus trace("Mixture_Model::update_voxel_tildew_vb");

    cout << "Doing voxel-wise tildew VB" << endl;

    ColumnVector m_tildew_new = m_tildew;

    SparseMatrix Lambda;
    Lambda = precision_lik;
    symmetric_addto(Lambda,D,mrf_precision);

    ColumnVector beta;
    multiply(precision_lik,m_tildew,beta);
    beta -= derivative_lik;

    float count = 0;
    for(int z = 0; z < mask.zsize(); z++)
      for(int y = 0; y < mask.ysize(); y++)
	for(int x = 0; x < mask.xsize(); x++)
	  if(mask(x,y,z))
	    {
	      ColumnVector sumneighs(nclasses);
	      sumneighs = 0;

	      int r = indices(x,y,z);   
	     
// 	      if(it>0 && r==17)
// 	      {
// 		OUT(r);
// 		RowVector cv(nclasses);
// 		cv = 0;		
// 		for(int c = 0; c < nclasses; c++)
// 		  {
// 		    cv(c+1) = m_tildew(c*num_superthreshold+r);
// 		  }
		
// 		matout(cv,"cv");
// 		RowVector w = logistic_transform(cv,lambda,log_bound);
// 		OUT(w);
// 	      }
	   
	      int xi=0,yi=0,zi=0;
	      for(unsigned int i = 0; i < connected_offsets.size(); i++) 
		{
		  xi = x+connected_offsets[i].x;
		  yi = y+connected_offsets[i].y;
		  zi = z+connected_offsets[i].z;
		  
		  if(mask(xi,yi,zi))
		    {
		      int r2 = indices(xi,yi,zi);

		      ColumnVector ck(nclasses);
		      ck = 0;
		     
		      DiagonalMatrix lam_vk(nclasses);
		      lam_vk = 0;

		      for(int c = 0; c < nclasses; c++)
			{
			  ck(c+1) = m_tildew(c*num_superthreshold+r2);
			  lam_vk(c+1,c+1) =  Lambda(c*num_superthreshold+r2,c*num_superthreshold+r);
			  
			}

		      sumneighs += lam_vk*ck;
// 		      if(it>0 && r==17)
// 			{
// 			  OUT(lam_vk);
// 			  OUT(ck);
// 			  OUT(lam_vk*ck);
// 			}
		    }
		}	    

	      ColumnVector betav(nclasses);
	      for(int c=0; c < nclasses; c++)
		{
		  betav(c+1) = beta(c*num_superthreshold+r);
		}

	      // get wtildecov
	      SymmetricMatrix wtildeprec(nclasses);
	      wtildeprec = 0;
	      SymmetricMatrix wtildecov(nclasses);
	      wtildecov = 0;

	      for(int c=0; c < nclasses; c++)
		{
		  wtildeprec(c+1,c+1) = Lambda(c*num_superthreshold+r,c*num_superthreshold+r);
		  for(int k=c+1; k < nclasses; k++)
		    {
		      wtildeprec(c+1,k+1) = Lambda(c*num_superthreshold+r,k*num_superthreshold+r); 
		    }
		}

	      try{
		wtildecov = wtildeprec.i();
	      }
	      catch(Exception& exp)
		{
		  OUT(exp.what());
		  matout(wtildeprec,"wtildeprec");
		  matout(betav.t(),"betav");
		  matout(sumneighs.t(),"sumneighs");
		  RowVector cv(nclasses);
		  cv = 0;		
		  for(int c = 0; c < nclasses; c++)
		    {
		      cv(c+1) = m_tildew(c*num_superthreshold+r);
		    }		
		  matout(cv,"cv");
		  RowVector w = logistic_transform(cv,lambda,log_bound);
		  matout(w,"w");
		  OUT(r);
		  exit(0);
		}
	    

	      ColumnVector wtilde = wtildecov*(betav - sumneighs);

	      wtilde = wtilde - mean(wtilde).AsScalar();
	
	      bool valid = true;
	      for(int c=0; c < nclasses; c++)
		{
		  if(std::fabs(wtilde(c+1)) > 10)
		    {
		      // 		      OUT(r);
		      valid = false;
		    }

		  // 		  if(wtilde(c+1) > 10)
		  // 		    {
		  // 		      wtilde(c+1) = 10;
		  // 		      count++;
		  
		  // 		    }
		  // 		  if(wtilde(c+1) < -10)
		  // 		    {
		  // 		      wtilde(c+1) = -10;
		  // 		      count++;
		  // 		    }
		}
	  
	      //if(it>0 && r==11335)
	      // 	      if(!valid)
// 		{	
// 		  OUT(r);		
// 		  RowVector w = logistic_transform(wtilde.t(),lambda,log_bound);
// 		  matout(wtildeprec,"wtildeprec");
// 		  matout(wtildecov,"wtildecov");
// 		  matout(betav.t(),"betav");
// 		  matout(sumneighs.t(),"sumneighs");
// 		  matout((betav-sumneighs).t(),"betav-sumneighs");
// 		  matout(wtilde.t(),"wtilde");
// 		  matout(w,"w");

// 		  RowVector cv(nclasses);
// 		  cv = 0;		
// 		  for(int c = 0; c < nclasses; c++)
// 		    {
// 		      cv(c+1) = m_tildew(c*num_superthreshold+r);
// 		    }		
// 		  matout(cv,"cv");
// 		  w = logistic_transform(cv,lambda,log_bound);
// 		  matout(w,"w");

// 		}	
    
	      if(valid || it<2)
		{
		  count++;
		  prec_tildew[r-1] = wtildeprec;
		  cov_tildew[r-1] = wtildecov;
		  
		  for(int c=0; c < nclasses; c++)
		    {
		      //m_tildew_new(c*num_superthreshold+r) = wtilde(c+1);
		      m_tildew(c*num_superthreshold+r) = wtilde(c+1);
		    }
 		}
	    }

    OUT(num_superthreshold - count);

    //    m_tildew = m_tildew_new;
  }

  void Mixture_Model::update_mrf_precision()
  {
    Tracer_Plus trace("Mixture_Model::update_mrf_precision"); 
    
    mrf_precision_hist.push_back(mrf_precision);
    
    float errorprecision = 1;
    float var = Sqr(errorprecision)*10;
    float c_0 = Sqr(errorprecision)/var;
    float b_0 = errorprecision/var;
    float c_g = nclasses*num_superthreshold/2.0 + c_0;
    float b_g = 1.0/(0.5*(quadratic(m_tildew,D) + trace_covariance_tildew_D)+1.0/b_0);
//     OUT(1/b_0);
//     OUT(0.5*(quadratic(m_tildew,D)));
//     OUT(0.5*trace_covariance_tildew_D);
//     OUT(1/b_g);
//     OUT(c_0);
//     OUT(c_g);
    float mrf_precision_new = std::exp(std::log(b_g)+MISCMATHS::lgam(c_g+1)-MISCMATHS::lgam(c_g));

    if(mrfprecmultiplier>0 && it>2)
      {
	float mrf_precision_old = mrf_precision_hist[mrf_precision_hist.size()-1];
	float mrf_precision_oldold = mrf_precision_hist[mrf_precision_hist.size()-2];     
	
	// 	OUT(mrfprecmultiplier);
	// 	OUT(mrf_precision_oldold);
	// 	OUT(mrf_precision_old);
	// 	OUT(mrf_precision_new);
	// 	OUT(mrf_precision);

	if(sign(mrf_precision_oldold - mrf_precision_old) != sign(mrf_precision_old-mrf_precision_new))
	  {
	    mrfprecmultiplier *= 0.5;
	  }

	mrf_precision = mrfprecmultiplier*(mrf_precision_new-mrf_precision_old)+mrf_precision_old;

	if(mrf_precision<=0)
	  {
	    mrf_precision = 1;
	    mrfprecmultiplier *= 0.5;
	  } 

	if(mrfprecmultiplier<1)mrfprecmultiplier=1;
	OUT(mrfprecmultiplier);
	
      }
    else
      {
	mrf_precision = mrf_precision_new;

	// check for convergence:
	if(it>10 && stopearly)
	  {
	    float mrf_precision_old = mrf_precision_hist[mrf_precision_hist.size()-1];
	    float mrf_precision_oldold = mrf_precision_hist[mrf_precision_hist.size()-2];
	    
// 	    OUT(mrf_precision);
// 	    OUT(mrf_precision_old);
// 	    OUT(mrf_precision_oldold);
//  	    OUT(std::fabs((mrf_precision-mrf_precision_old)/mrf_precision_old));
//  	    OUT(std::fabs((mrf_precision-mrf_precision_oldold)/mrf_precision_oldold));

	    float precision = 0.005;
	    
	    if(std::fabs((mrf_precision-mrf_precision_old)/mrf_precision_old) < precision && std::fabs((mrf_precision-mrf_precision_oldold)/mrf_precision_oldold) < precision)
	      {
		it=niters;
	      }
	  }
      }
    
  }

  void Mixture_Model::calculate_taylor_lik()
  {
    Tracer_Plus trace("Mixture_Model::calculate_taylor_lik");

    cout << "Doing 2nd Order Taylor Expansion of Likelihood" << endl; 

    // build up precision/hessian matrix for tildew (aka x)
    //    float A = 1.0/std::sqrt(2*M_PI);
    float lamsqr = Sqr(lambda*log_bound);

    derivative_lik.ReSize(num_superthreshold*nclasses);
    derivative_lik = 0;
    precision_lik.ReSize(num_superthreshold*nclasses,num_superthreshold*nclasses);

    for(int r = 1; r<=num_superthreshold; r++)
      {
	RowVector wtildetmp(nclasses);
	wtildetmp = 0;
	for(int c=1; c <= nclasses; c++)
	  {	
	    // calculate terms for use in calculating hessian and deriv of likelihood:
	    wtildetmp(c) = m_tildew((c-1)*num_superthreshold+r);
	  }

	// LT of y = LT of demean(y)
	    
	const RowVector wtilde = wtildetmp - mean(wtildetmp,2).AsScalar();

	RowVector w = logistic_transform(wtilde,lambda,log_bound);

	if(it>0 && r==-1)
	  {
	    OUT(wtildetmp);
	    OUT(mean(wtildetmp,2).AsScalar());
	    OUT(wtilde);
	    matout(w,"w");
	  }

	RowVector gamma(nclasses);
	gamma = 0;
	vector<double> R(nclasses,0);
	    
	double P = 0;
	double h = 0;

	for(int c=1; c <= nclasses; c++)
	  {	
	    // calculate terms for use in calculating hessian and deriv of likelihood:	    
	    h += w(c)*dists[c-1]->pdf(Y(r));
	    R[c-1] = boundexp(wtilde(c)/(lambda*log_bound));
	    P += R[c-1];
	    if(it>0 && r==-1)
	      {
		OUT(dists[c-1]->pdf(Y(r)));
	      }
	  }

	if(it>0 && r==-1)
	  {
	    OUT(r);
	    OUT(Y(r));
	    matout(w,"w");
	    matout(wtilde,"wtilde");
	    for(int c=1; c <= nclasses; c++)
	      OUT(R[c-1]);
	    OUT(P);
	    OUT(h);
	    OUT(Sqr(h));
	  }

	// calculate dw_k/dx and d^2w_k/dx^2 
	vector<SymmetricMatrix> dwdydy(nclasses);
	vector<ColumnVector> dwdy(nclasses);
	for(int k=1; k <= nclasses; k++)
	  {	
	    dwdydy[k-1].ReSize(nclasses);
	    dwdydy[k-1] = 0;
	    dwdy[k-1].ReSize(nclasses);
	    dwdy[k-1] = 0;
	    
	    for(int c2=1; c2 <= nclasses; c2++)
	      {		
		
		if(c2==k)
		  {
		    dwdydy[k-1](k,k) = R[k-1]/(lamsqr*P)*(1-3*R[k-1]/P+2*Sqr(R[k-1]/P));
		    dwdy[k-1](k) = R[k-1]*(1-R[k-1]/P)/(lambda*log_bound*P);
		  }
		else
		  {
		    dwdydy[k-1](c2,c2) = R[k-1]*R[c2-1]/(lamsqr*Sqr(P))*(2*R[c2-1]/P-1);	
		    dwdy[k-1](c2) = -R[k-1]*R[c2-1]/(lambda*log_bound*Sqr(P));
		  }

		for(int c3=c2+1; c3 <= nclasses; c3++)
		  {
		    if(c2==k) // we know that c3 != c2
		      {
			dwdydy[k-1](k,c3) = R[k-1]*R[c3-1]/(lamsqr*Sqr(P))*(2*R[k-1]/P-1);
		      }
		    else if(c3==k) // we know that c3 != c2
		      {
			dwdydy[k-1](k,c2) = R[k-1]*R[c2-1]/(lamsqr*Sqr(P))*(2*R[k-1]/P-1);
		      }
		    else
		      {
			dwdydy[k-1](c2,c3) = 2*R[k-1]*R[c2-1]*R[c3-1]/(lamsqr*Sqr(P)*P);
		      }
		  }
	      }
	  }

	// calculate d^2f/dw^2 and df/dw
	SymmetricMatrix dfdwdw(nclasses);
	dfdwdw = 0;
	ColumnVector dfdw(nclasses);
	dfdw = 0;

	//	float premult = Sqr(A/h);

	for(int c=1; c <= nclasses; c++)
	  {
	    dfdwdw(c,c) =  Sqr(dists[c-1]->pdf(Y(r))/h);
	    dfdw(c) = -dists[c-1]->pdf(Y(r))/h;

	    for(int c2=c+1; c2 <= nclasses; c2++)
	      {
		dfdwdw(c,c2) = dists[c-1]->pdf(Y(r))*dists[c2-1]->pdf(Y(r))/Sqr(h);
		//dfdwdw(c2,c) = dists[c-1]->pdf(Y(r))*dists[c2-1]->pdf(Y(r))/Sqr(h);
	      }
	  }

	if(it>0 && r==-1)
	  {
	    matout(dfdwdw,"dfdwdw");
	    matout(dfdw,"dfdw");
	  }

	// Now fill up precision/Hessian for tildew (aka x) for this voxel
	SymmetricMatrix hessanal(nclasses);
	hessanal = 0;
	ColumnVector derivanal(nclasses);
	derivanal = 0;

	for(int k=1; k <= nclasses; k++)
	  {
	    // diagonal terms k=j
	    float sum_l = 0;
	    float sum_l2 = 0;
	    for(int l=1; l <= nclasses; l++)
	      {
		float sum_m = 0;
		for(int m=1; m <= nclasses; m++)
		  {
		    sum_m += dfdwdw(m,l)*dwdy[m-1](k);
		  }
		if(it==-1 && r==40 && k==1)
		  {
		    matout(dfdwdw,"dfdwdw");
		    matout(dwdy[l-1],"dwdy");
		    OUT(l);
		    OUT(sum_m);
		    OUT(dwdydy[l-1](k,k));
		    OUT(dwdy[l-1](k));
		    OUT(dfdw(l));
		  }

		sum_l += sum_m*dwdy[l-1](k) + dfdw(l)*dwdydy[l-1](k,k);

		sum_l2 += dfdw(l)*dwdy[l-1](k);
	      }

	    derivanal(k) = sum_l2;

	    if(it==-1 && r==40 && k==1)
	      {
		OUT(sum_l);
	      }
	    hessanal(k,k) = sum_l;

	    // off-diagonal terms jk (j is called n here)
	    for(int n=k+1; n <= nclasses; n++)
	      {
		float sum_l = 0;
		for(int l=1; l <= nclasses; l++)
		  {
		    float sum_m = 0;
		    for(int m=1; m <= nclasses; m++)
		      {
			sum_m += dfdwdw(m,l)*dwdy[m-1](k);
		      }
		    sum_l += sum_m*dwdy[l-1](n) + dfdw(l)*dwdydy[l-1](n,k);		 
		  }

		hessanal(n,k) = sum_l;
	      }
	  }

	for(int k=1; k <= nclasses; k++)
	  {
	    derivative_lik((k-1)*num_superthreshold+r) += derivanal(k);

	    precision_lik.addto((k-1)*num_superthreshold+r,(k-1)*num_superthreshold+r,hessanal(k,k));	    
	    for(int l=k+1; l <= nclasses; l++)
	      {
		precision_lik.addto((k-1)*num_superthreshold+r,(l-1)*num_superthreshold+r,hessanal(k,l));
		precision_lik.addto((l-1)*num_superthreshold+r,(k-1)*num_superthreshold+r,hessanal(k,l));
	      }	
	  }      

	////////////////////
	if(it>0 && r==-1)  
	  {
	    matout(derivanal,"derivanal");
	    matout(hessanal,"hessanal");
	  }

// 	if(it>5 && r==-1)  
// 	  {
// 	    SmmVoxelFunction smmfunc(Y(r),dists,lambda,log_bound);
// 	    ColumnVector deriv1 = gradient(wtilde.t(),smmfunc,1e-2,1);
// 	    ColumnVector deriv2 = gradient(wtilde.t(),smmfunc,1e-2,2);
// 	    ColumnVector deriv4 = gradient(wtilde.t(),smmfunc,1e-2,4);
// 	    OUT(r);
// 	    matout(derivanal,"derivanal");
// 	    matout(deriv2, "deriv2");
// 	    matout(deriv4, "deriv4");
// 	    SymmetricMatrix hess1 = hessian(wtilde.t(),smmfunc,1e-2,1);
// 	    SymmetricMatrix hess2 = hessian(wtilde.t(),smmfunc,1e-2,2);
// 	    SymmetricMatrix hess4 = hessian(wtilde.t(),smmfunc,1e-2,4);
// 	    matout(hessanal,"hessanal");
// 	    matout(hess2,"hess2");
// 	    matout(hess4,"hess4");
// 	  }
	//////////////
      }
  }

  void Mixture_Model::calculate_trace_tildew_D()
  {
    Tracer_Plus trace("Mixture_Model::calculate_trace_tildew_D");

    // Now calculate trace    
    float trace_new = 0; 
    DiagonalMatrix tmp1(num_superthreshold*nclasses);
    tmp1 = 0;

    for(int r=1; r <= num_superthreshold; r++)
      {	
// 	OUT(r);
// 	matout(cov_tildew[r-1],"cov_tildew[r-1]");
	for(int c=1; c <= nclasses; c++)
	  {
	    tmp1((c-1)*num_superthreshold+r,(c-1)*num_superthreshold+r) = cov_tildew[r-1](c,c);
	  }
      }

    SparseMatrix tmpres(num_superthreshold*nclasses,num_superthreshold*nclasses);
    multiply(tmp1,D,tmpres);
    trace_new = tmpres.trace();

    OUT(trace_new);

    trace_covariance_tildew_D = trace_new;	

    OUT(trace_covariance_tildew_D);  
  }

  ReturnMatrix sum_transform(const RowVector& wtilde, float log_bound)
  {
    //Tracer_Plus trace("sum_transform");

    // converts from wtilde (aka x) to y
    RowVector ret_y  = log_bound*wtilde/std::sqrt(wtilde.SumSquare());    
    
    ret_y.Release();
    return ret_y;    
  }

  ReturnMatrix logistic_transform(const RowVector& py,float lambda,float log_bound)
  {
    //    Tracer_Plus trace("logistic_transform");

    // LT of y = LT of demean(y)
    const RowVector y = py - mean(py,2).AsScalar();

    // converts from y to w
    
    int nclasses = y.Ncols();
    double sum = 0.0;
    RowVector ret_weights(nclasses);
    ret_weights = 0;

    double phi = lambda*log_bound;

    for(int c=0; c < nclasses; c++)
      {
	sum += boundexp(y(c+1)/phi);
      }
    
    for(int c=0; c < nclasses; c++)
      {
	ret_weights(c+1) = boundexp(y(c+1)/phi)/sum;	
      }

    if(ret_weights(2)>1)
      {
	OUT(phi);
	OUT(y);
	OUT(sum);
	OUT(ret_weights);
	OUT(boundexp(y(2)/phi));

      }
    
    ret_weights.Release();
    return ret_weights;
    
  }

  ReturnMatrix inv_transform(const RowVector& w,float lambda,float log_bound,float initmultiplier)
  {
    Tracer_Plus trace("inv_transform");

    // converts from w to wtilde (aka x)

    int nclasses = w.Ncols();
    
    RowVector ret_wtilde(nclasses);
    ret_wtilde = 0;


    ///////////////////
//     double phi = lambda*log_bound;

//     float noise = 0;
//     int max = 1;
//     float maxret = phi*initmultiplier+normrnd().AsScalar()*phi*noise;

//     for(int c=0; c < nclasses; c++)
//       {
// 	if(w(c+1)>w(max)) max = c+1;
// 	ret_wtilde(c+1) = -phi*initmultiplier+normrnd().AsScalar()*phi*noise;

// 	while(ret_wtilde(c+1) >= maxret)
// 	{	  
// 	  ret_wtilde(c+1) = -phi*initmultiplier+normrnd().AsScalar()*phi*noise;
// 	}
//       }
//     ret_wtilde(max) = maxret;
///////////////////

    for(int c=0; c < nclasses; c++)
      {
	if(w(c+1)==1)
	  ret_wtilde(c+1) = log_bound*initmultiplier;
// 	else if(w(c+1)==0.5)
// 	  ret_wtilde(c+1) = 0;
	else
	  ret_wtilde(c+1) = -log_bound*initmultiplier;
      }

    ret_wtilde.Release();
    return ret_wtilde; 
  }

  void Mixture_Model::save_weights(const ColumnVector& pmtildew, const char* affix, bool usesamples) 
  {
    Tracer_Plus trace("Mixture_Model::save_weights");

    //    vector<volume<float> > w_means(nclasses);   
    vector<volume<float> > logistic_w_means(nclasses);   
    vector<ColumnVector> weights;
    vector<vector<vector<float> > > weights_samps;
    vector<vector<vector<float> > > tildew_samps;

    int nsamps = 50;

    OUT("Calculating weights");
    if(nonspatial || !usesamples)
      {
	get_weights(weights, pmtildew);
      }
    else
      {	
	get_weights2(weights,weights_samps,tildew_samps,nsamps,pmtildew);
      }

    vector<volume4D<float> > w_samples(nclasses);
    vector<volume4D<float> > tildew_samples(nclasses);

    for(int c=0; c < nclasses; c++)
      {
	logistic_w_means[c].reinitialize(xsize,ysize,zsize);
	logistic_w_means[c] = 0.0;
	w_means[c].reinitialize(xsize,ysize,zsize);
	w_means[c] = 0.0;
	w_samples[c].reinitialize(xsize,ysize,zsize,nsamps);
	w_samples[c] = 0.0;
	tildew_samples[c].reinitialize(xsize,ysize,zsize,nsamps);
	tildew_samples[c] = 0.0;

	for(int z = 0; z < mask.zsize(); z++)    
	  for(int y = 0; y < mask.ysize(); y++)	
	    for(int x = 0; x < mask.xsize(); x++)
	      if(mask(x,y,z))
		{
		  w_means[c](x,y,z) = weights[c](indices(x,y,z));
		  logistic_w_means[c](x,y,z) = pmtildew(c*num_superthreshold+indices(x,y,z));
		  //		  if(w_means[c](x,y,z)>0.5) {OUT(indices(x,y,z)); OUT(w_means[c](x,y,z));}
		  
		  if(!nonspatial && usesamples)
		    for(int s = 0; s < nsamps; s++) 
		      {
			w_samples[c](x,y,z,s) = weights_samps[indices(x,y,z)-1][s][c];
			tildew_samples[c](x,y,z,s) = tildew_samps[indices(x,y,z)-1][s][c];
		      }
		}
	
	copybasicproperties(spatial_data,logistic_w_means[c]);		
	save_volume(logistic_w_means[c], LogSingleton::getInstance().appendDir("logistic_w"+num2str(c+1)+"_mean"+affix));	

	copybasicproperties(spatial_data,w_means[c]);		
	save_volume(w_means[c], LogSingleton::getInstance().appendDir("w"+num2str(c+1)+"_mean"+affix));	
	
	copybasicproperties(spatial_data,w_samples[c]);
	save_volume4D(w_samples[c], LogSingleton::getInstance().appendDir("w"+num2str(c+1)+"_samples"+affix));

	copybasicproperties(spatial_data,tildew_samples[c]);
	save_volume4D(tildew_samples[c], LogSingleton::getInstance().appendDir("logistic_w"+num2str(c+1)+"_samples"+affix));
      }

  }

  void Mixture_Model::save() 
  {
    Tracer_Plus trace("Mixture_Model::save");

    save_volume(spatial_data, LogSingleton::getInstance().appendDir("spatial_data"));
    save_volume(mask, LogSingleton::getInstance().appendDir("mask"));

    // save weights
    save_weights(m_tildew, "", true);

//     update_tildew_scg();
//     save_weights(volinfo,m_tildew,"_scg",false);

    calculate_props(w_means,dists,mask);

    // save distribution params
    ColumnVector means(nclasses);
    ColumnVector vars(nclasses);
    ColumnVector props(nclasses);
    means = 0;
    vars = 0;
    props = 0;

    for(int c = 0; c < nclasses; c++)
      {
	means(c+1) = dists[c]->getmean();
	vars(c+1) = dists[c]->getvar();
	props(c+1) = dists[c]->getprop();
      }

    for(int c=0; c < nclasses; c++)
      {
	write_ascii_matrix(means,LogSingleton::getInstance().appendDir("mu_mean"));
	write_ascii_matrix(vars,LogSingleton::getInstance().appendDir("var_mean"));
	write_ascii_matrix(props,LogSingleton::getInstance().appendDir("prop_mean"));
      }

    if(!nonspatial && !fixmrfprec)
      {
	miscplot newplot;
	newplot.add_xlabel("Iterations");    
	newplot.set_xysize(610,300);
	newplot.timeseries(vector2ColumnVector(mrf_precision_hist).t(), LogSingleton::getInstance().appendDir("mrfprechist"), "MRF Precision", 0,400,3,0,false);
      }

    if(updatetheta)
      {
	miscplot newplot;
	newplot.add_xlabel("Iterations");    
	newplot.set_xysize(610,300);
	newplot.timeseries(vector2ColumnVector(meanhist).t(), LogSingleton::getInstance().appendDir("meanhist"), "class 1 mean", 0,400,3,0,false);
      }


    write_vector(mrf_precision_hist, LogSingleton::getInstance().appendDir("mrf_precision_hist"));

  }

  void Mixture_Model::get_weights(vector<ColumnVector>& weights, const ColumnVector& pmtildew)
  {
    weights.resize(nclasses);
    for(int c = 0; c < nclasses; c++)
      {
	weights[c].ReSize(num_superthreshold);
	weights[c] = 0;
      } 

    for(int r = 1; r<=num_superthreshold; r++)
      {       
	RowVector wtilde(nclasses);
	wtilde = 0;
	
	for(int c=0; c < nclasses; c++)
	  wtilde(c+1) = pmtildew(c*num_superthreshold+r);

	RowVector w = logistic_transform(wtilde,lambda,log_bound);

	//  	OUT(r);
	//  	OUT(wtilde);
	//  	OUT(w);

	for(int c=0; c < nclasses; c++)
	  weights[c](r) = w(c+1);
      }
  }


  void Mixture_Model::get_weights2(vector<ColumnVector>& weights, vector<vector<vector<float> > >& weights_samps, vector<vector<vector<float> > >& tildew_samps, int nsamps, const ColumnVector& pmtildew)
  {

    Tracer_Plus trace("Mixture_Model::get_weights2");

    weights.resize(nclasses);
    for(int c = 0; c < nclasses; c++)
      {
	weights[c].ReSize(num_superthreshold);
	weights[c] = 0;
      } 
    weights_samps.reserve(num_superthreshold);
    tildew_samps.reserve(num_superthreshold);

    for(int r = 1; r<=num_superthreshold; r++)
      {	
	//	OUT(r);
	RowVector wtilde(nclasses);	
	wtilde = 0;	
	for(int c=0; c < nclasses; c++)
	  {
	    wtilde(c+1) = pmtildew(c*num_superthreshold+r);	    
	  }
	SymmetricMatrix wtildecov = cov_tildew[r-1];

// 	if(it>0 && r==-1)
// 	  {
// 	    OUT(r);
// 	    matout(wtilde,"wtilde");
// 	    matout(wtildecov,"wtildecov");
// 	  }

	// now sample
	Matrix wtildesamps_mat = mvnrnd(wtilde,wtildecov,nsamps);

	vector<vector<float> > wsamps;
	wsamps.reserve(nsamps);

	vector<vector<float> > wtildesamps;
	wtildesamps.reserve(nsamps);

	RowVector w(nclasses);
	w = 0;
	for(int s=1; s <= nsamps; s++)
	  {
	    RowVector wrow = logistic_transform(wtildesamps_mat.Row(s),lambda,log_bound);
	    w += wrow;

// 	    if(r<2 && s<2)
// 	      {
// 		OUT(s);
// 		OUT(wrow);
// 		OUT(wtildesamps_mat.Row(s));
// 	      }

	    vector<float> wvec(nclasses);
	    for(int c=0; c < nclasses; c++)
	      wvec[c] = wrow(c+1);
	    wsamps.push_back(wvec);

	    for(int c=0; c < nclasses; c++)
	      wvec[c] = wtildesamps_mat(s,c+1);

	    wtildesamps.push_back(wvec);
	  }

	weights_samps.push_back(wsamps);
	tildew_samps.push_back(wtildesamps);

	for(int c=0; c < nclasses; c++)
	  weights[c](r) = w(c+1)/float(nsamps);
      }
  }

  void ggmfit(const RowVector& dat, vector<Distribution*>& pdists, bool useprops)
  {// fit a mixture of a Gaussian and multiple Gamma functions to the histogram
  
    //normalise data 
    float datamean = mean(dat,2).AsScalar();
    float datastdev= stdev(dat,2).AsScalar();
    RowVector data=(dat - datamean)/datastdev;

    int numdata = dat.Ncols();
    int nummix = 3;
    RowVector means(nummix);
    RowVector vars(nummix);
    RowVector props(nummix);
    means = 0;
    vars = 0;
    props = 0;
//      Params=zeros(1,nummix);
    //    float logprobY = 1.0;

    props = std::pow(float(nummix),float(-1.0));

    Matrix tmp1;
    tmp1 = data * data.t() / numdata;
    vars = tmp1.AsScalar();

    float Dmin1, Dmax1, IntSize;
    Dmin1 =  min(data).AsScalar(); Dmax1 = max(data).AsScalar();
    IntSize = Dmax1 / nummix;

    means(1)=mean(data,2).AsScalar(); 
    for (int ctr=2; ctr <= means.Ncols(); ctr++){
      means(ctr) =  Dmax1 - (ctr - 1.5) * IntSize; 
    }
    means(2)=means(1)+2*sqrt(vars(1));
    //means(2)=means(1)+ 0.6*(Dmax-means(1));
    if(nummix>2)
      //means(3)=means(1)-0.6*(means(1)-Dmin);
            means(3)=means(1)-2*sqrt(vars(1));

//      epsilon = eps;
//      if(epsilon >=0 && epsilon < 0.0000001)
//        {epsilon = std::log(float(dat.Ncols()))/1000 ;}


    // ggmfit
    float log_p_y_theta = 1.0;
    float old_ll = 2.0;
    float g_eps = 0.000001;
    int it_ctr = 0;
    double Dmax, Dmin;
   
    Dmax = 2 * data.Maximum();
    Dmin = -2 * data.Minimum();

    //resize means, vars and props
    if(nummix > 3)
      nummix = 3;
    means = means.Columns(1,nummix);
    vars  = vars.Columns(1,nummix);
    props = props.Columns(1,nummix);

    means(1) = -2*mean(data,2).AsScalar();

    Matrix prob_K__y_theta;
    Matrix kdata;
    RowVector prob_Y__theta;RowVector Nbar;
    RowVector mubar;RowVector sigmabar;RowVector pibar;
    Matrix p_ygx(nummix,numdata);
    //    float offset = 0.0;
    float const2;
    Matrix negdata(data);
    negdata = -data;

    while((it_ctr<10) ||
	  ((std::abs(old_ll - log_p_y_theta)>g_eps) && (it_ctr<100)))
      { // fit GGM
	
 	it_ctr++;
	//offset = (std::min(0.2,1-props(1)))*std::sqrt(vars(1));

// 	//make sure all variances are acceptable
 	for(int ctr1=1; ctr1 <=nummix; ctr1++)
 	  if(vars(ctr1)<0.0001){
 	    vars(ctr1) = 0.0001;
 	  }

 	p_ygx = 0.0;
 	p_ygx.Row(1) << normpdf(data,means(1),vars(1));
       
 	const2 = (2.6-props(1))*sqrt(vars(1))+means(1); //min. nht level
 
	means(2) = (std::max(means(2), std::max(0.001,
	   0.5 * ( const2 + std::sqrt( const2*const2 + 4*vars(2))))));
	vars(2)  = std::max(std::min(vars(2), 0.5*std::pow(means(2),2)),0.0001);
	p_ygx.Row(2) << gammapdf(data,means(2),vars(2));
   
 	if(nummix>2){
	  const2 = (2.6-props(1))*sqrt(vars(1))-means(1);
	
	  means(3) = -(std::max(-means(3), std::max(0.001,
	      0.5 * ( const2 + std::sqrt( const2*const2 + 4*vars(3))))));
	  vars(3)  = std::max(std::min(vars(3), 0.5*std::pow(means(3),2)),0.0001);
 	  p_ygx.Row(3) << gammapdf(negdata,-means(3),vars(3));
	}

 	tmp1 = SP(props.t()*ones(1,numdata),p_ygx);
 	prob_Y__theta << sum(tmp1,1);
	
	//deal with non-modelled voxels
	for(int ctr=1; ctr<=tmp1.Ncols(); ctr++)
	  if(prob_Y__theta(ctr) < 0.0001)
	    prob_Y__theta(ctr) = 0.0001;

 	old_ll = log_p_y_theta;
 	log_p_y_theta = log(prob_Y__theta).Sum();
	// 	cerr << "calculated log_prob_Y__theta" <<endl;
	// 	cerr << old_ll << "   " << log_p_y_theta << "   " 
	//	cerr << float(std::abs(old_ll - log_p_y_theta)) << endl;
 	if((it_ctr<10) ||
	   ((std::abs(old_ll - log_p_y_theta)>g_eps) && (it_ctr<100))){//update
	  
 	  prob_K__y_theta = SP(tmp1,pow(ones(nummix,1)*prob_Y__theta,-1));
 	  Nbar << sum(prob_K__y_theta,2).t();
	  for(int ctr=1; ctr<=Nbar.Ncols(); ctr++)
	    if(Nbar(ctr) < 0.0001 * numdata)
	      Nbar = Nbar + 0.0001;
 	  pibar= Nbar / numdata;
	  // 	  cerr << "pibar :" << pibar << endl;
 	  kdata = ones(nummix,1)*data;
 	  mubar <<SP(sum(SP(kdata,prob_K__y_theta),2).t(),pow(Nbar,-1)); 
	  // 	  cerr << "mubar :" << mubar << endl;

 	  kdata -= mubar.t()*ones(1,numdata);
 	  kdata = pow(kdata,2);
 	  sigmabar << SP(sum(SP(kdata,prob_K__y_theta),2).t(),pow(Nbar,-1));
      
 	  means = mubar;
 	  vars  = sigmabar;

	  if(useprops)
	    props = pibar;
 	}//update
      } //while loop


    data = data*datastdev+datamean;
    means = means*datastdev+datamean;

    for(int c=0; c < nummix; c++)
      vars(c+1) = Sqr(std::sqrt(vars(c+1))*datastdev);

    if(((it_ctr%1)==0)||(it_ctr==1)){
      string pad = "-";
      miscplot newplot;

      newplot.ggmfit(data,means,vars,props,
		     LogSingleton::getInstance().appendDir("init_mmfit.png"),string("Initial Fit"),0.0);
    }    

    for(unsigned int c=0; c < pdists.size(); c++)
      {	
	pdists[c]->setparams(means(c+1),vars(c+1),props(c+1));
	  
 	OUT(c);
 	OUT(pdists[c]->getmean());
 	OUT(pdists[c]->getvar());
      }
  }

  void calculate_props(const vector<volume<float> >& w_means, vector<Distribution*>& dists, const volume<int>& mask)
  {
    int nclasses = w_means.size();
    for(int c = 0; c < nclasses; c++)
      {	
	OUT(c);
	float sumweights = 0;
	int num_superthreshold = 0;

	for(int z = 0; z < mask.zsize(); z++)    
	  for(int y = 0; y < mask.ysize(); y++)	
	    for(int x = 0; x < mask.xsize(); x++)
	      if(mask(x,y,z))
		{
		  float wtmp = w_means[c](x,y,z);
		  sumweights += wtmp;
		  num_superthreshold++;
		}

	OUT(num_superthreshold);

	dists[c]->setparams(dists[c]->getmean(),dists[c]->getvar(),sumweights/num_superthreshold);
      }
}

  void plot_ggm(const vector<volume<float> >& w_means, const vector<Distribution*>& dists, const volume<int>& mask, const ColumnVector& Y)
  {
    OUT("plot_ggm");
    // update means and vars using ML estimates, given weights are known    
    int nclasses = w_means.size();

    OUT(nclasses);

    int nummix = 3;
    RowVector means(nummix);
    RowVector vars(nummix);
    RowVector props(nummix);
    means = 0;
    vars = 0;
    props = 0;

    for(int c = 0; c < nclasses; c++)
      {
	means(c+1) = dists[c]->getmean();
	vars(c+1) = dists[c]->getvar();
	props(c+1) = dists[c]->getprop();
      }

    if(nclasses==2)
      {
	means(3) = 0;
	vars(3) = 0.1;
	props(3) = 0;
      }

    OUT(means);
    OUT(vars);
    OUT(props);

    miscplot newplot;
 
    newplot.ggmfit(Y.t(),means,vars,props,LogSingleton::getInstance().appendDir("final_mmfit.png"),string("Final Fit"),0.0);
    
  }

  void make_ggmreport(const vector<volume<float> >& w_means, const vector<Distribution*>& dists, const volume<int>& mask, const volume<float>& spatial_data, bool zfstatmode, bool overlay, const volume<float>& epivol, float thresh, bool nonspatial, bool updatetheta, const string& data_name)
  {
    int nclasses = dists.size();

    RowVector means(nclasses);
    RowVector vars(nclasses);
    RowVector props(nclasses);
    means = 0;
    vars = 0;
    props = 0;

    for(int c = 0; c < nclasses; c++)
      {
	means(c+1) = dists[c]->getmean();
	vars(c+1) = dists[c]->getvar();
	props(c+1) = dists[c]->getprop();
      }

    LogSingleton::getInstance().setLogFile(string("MM.html"));

    Log& htmllog = LogSingleton::getInstance();
      {	
//  	OUT(htmllog.getDir());
//  	OUT(htmllog.getLogFileName());

	htmllog << "<HTML> " << endl
		<< "<TITLE>Mixture Model fit for" << data_name << "</TITLE>" << endl
		<< "<BODY BACKGROUND=\"file:" << getenv("FSLDIR") 
		<< "/doc/images/fsl-bg.jpg\">" << endl 
		<< "<hr><CENTER><H1>Mixture Model fit for<br>" << data_name << " </H1>"<< endl;
     	
	htmllog << "<hr><p>" << endl;
      }

      {
//        volume<float> map1;
//        volume<float> map2;

//        map1 = threshold(spatial_data,float(0.0), 
//  		       spatial_data.max());
//        map2 = threshold(spatial_data,spatial_data.min(), 
//  		       float(-0.0));
      
        volume<float> newvol; 
        miscpic newpic;
      
	newvol = spatial_data;

 	char instr[10000];
	
 	sprintf(instr," ");
 	strcat(instr,"-s 2 ");	
	strcat(instr,(string("-i ")+num2str(spatial_data.min())+" "+num2str(spatial_data.max())+" ").c_str());
 	strcat(instr,"-A 950 ");
 	strcat(instr,LogSingleton::getInstance().appendDir("spatial_data.png").c_str());

 	newpic.set_title(string("Raw spatial map"));
	
 	newpic.slicer(newvol, instr);
      

	htmllog << "<img BORDER=0 SRC=\"spatial_data.png\"><p>" << endl;
      }

      if(overlay)
	{	
	  volume<float> map;	  

	  map = w_means[1];
	
 	  if(thresh>0)
 	    map.threshold(0.5);	 
	  else
	    thresh=0;

	  volume<float> newvol; 
	  miscpic newpic;

	  volume<float> epivoltmp = epivol;

	  save_volume(map, LogSingleton::getInstance().appendDir("map"));
	  save_volume(epivoltmp, LogSingleton::getInstance().appendDir("epivol"));
	 
	  newpic.overlay(newvol, epivoltmp, map, map, 
			 epivol.percentile(0.01),
			 epivol.percentile(0.99),
			 float(thresh), float(1.0), float(0.0), float(0.0),
			 0, 0);
		
	  char instr[10000];
	
	  sprintf(instr," ");
	  strcat(instr,"-l render1 -s 2 ");
	  strcat(instr,"-A 950 ");
	  strcat(instr,LogSingleton::getInstance().appendDir("actprobmap.png").c_str());     

	  string tit = "Activation prob map";
	  if(thresh>0)
	    tit += " thresholded at p>" + num2str(thresh);
	  newpic.set_title(tit);	    
	
	  newpic.set_cbar(string("y"));
	  newpic.slicer(newvol, instr); 
	
	  htmllog << "<img BORDER=0 SRC=\"actprobmap.png\">" << endl;
	  htmllog << "<p>" << endl;

	}
      else
	{			
	  volume<float> newvol; 
	  miscpic newpic;

	  newvol = w_means[1];

	  if(thresh>0)
	    newvol.threshold(0.5);

	  char instr[10000];
	
	  	      //volume<float> tmp = newvol/newvol.max();
	      volume<float> tmp = spatial_data;
// 	      newpic.overlay(newvol, newvol, tmp, tmp,
// 			     float(0),float(1),
// 			     float(2), float(3), float(0.0), float(0.0),
// 			     0, 0, &volinfo);
	      newpic.overlay(newvol, tmp,newvol, newvol,
			     spatial_data.percentile(0.01),
			     spatial_data.percentile(0.99),
			     float(thresh), float(1.0), float(0.0), float(0.0),
			     0, 0);

	  sprintf(instr," ");
	  strcat(instr,"-s 2 ");
	  //	strcat(instr,"-i 0 1 ");
	  strcat(instr,"-A 950 ");
	  strcat(instr,LogSingleton::getInstance().appendDir("actprobmap.png").c_str());      

	  string tit = "Activation prob map";
	  if(thresh>0)
	    tit += " thresholded at p>" + num2str(thresh);
	  newpic.set_title(tit);	    
	  
	  newpic.set_cbar(string("y"));
	  newpic.slicer(newvol, instr); 
		
	  htmllog << "<img BORDER=0 SRC=\"actprobmap.png\">" << endl;
	  htmllog << "<p>" << endl;
	}
      
      if(!zfstatmode)
	{	
	  // Output deactivation
	  if(overlay)
	    {	
	      volume<float> map;
	      map = w_means[2];
	
	      if(thresh>0)
		map.threshold(0.5);
	      else
		thresh=0;

	      volume<float> newvol; 
	      
	      miscpic newpic;

	      volume<float> epivoltmp = epivol;

	      newpic.overlay(newvol, epivoltmp, map, map, 
			     epivol.percentile(0.01),
			     epivol.percentile(0.99),
			     float(thresh), float(1.0), float(0.0), float(0.0),
			     0, 0);
		
	      char instr[10000];
	
	      sprintf(instr," ");
	      strcat(instr,"-s 2 ");
	      strcat(instr,"-A 950 ");
	      strcat(instr,LogSingleton::getInstance().appendDir("deactprobmap.png").c_str());

	      string tit = "Deactivation prob map";
	      if(thresh>0)
		tit += " thresholded at p>" + num2str(thresh);
	      newpic.set_title(tit);
	
	      newpic.set_cbar(string("y"));
	      newpic.slicer(newvol, instr); 
	
	      htmllog << "<img BORDER=0 SRC=\"deactprobmap.png\">" << endl;
	      htmllog << "<p>" << endl;

	    }
	  else
	    {	      
	    	    
	      volume<float> newvol; 
	      miscpic newpic;
	
	      newvol = w_means[2];
	      if(thresh>0)
		newvol.threshold(0.5);

	      char instr[10000];
	
	      //volume<float> tmp = newvol/newvol.max();
	      volume<float> tmp = spatial_data;
// 	      newpic.overlay(newvol, newvol, tmp, tmp,
// 			     float(0),float(1),
// 			     float(2), float(3), float(0.0), float(0.0),
// 			     0, 0, &volinfo);
	      newpic.overlay(newvol, tmp,newvol, newvol,
			     spatial_data.percentile(0.01),
			     spatial_data.percentile(0.99),
			     float(thresh), float(1.0), float(0.0), float(0.0),
			     0, 0);

	      sprintf(instr," ");
	      strcat(instr,"-s 2 ");
	      //	      strcat(instr,"-i 0 1 ");
	      strcat(instr,"-A 950 ");
	      strcat(instr,LogSingleton::getInstance().appendDir("deactprobmap.png").c_str());      
	      string tit = "Deactivation prob map";
	      if(thresh>0)
		tit += " thresholded at p>" + num2str(thresh);
	      newpic.set_title(tit);

	      newpic.set_cbar(string("y"));
	      newpic.slicer(newvol, instr); 
	      htmllog << "<img BORDER=0 SRC=\"deactprobmap.png\">" << endl;
	      htmllog << "<p>" << endl;
	    }
	}
      
      {//Output GGM fit
	miscplot newplot;

	htmllog << "<A><img BORDER=0 SRC=\"final_mmfit.png\"></A><p>" << endl;
      }
      
      {
	htmllog << "<br> Mixture Model fit <br>" << endl
		<< "<br> &nbsp; Means :  " << means << endl
		<< "<br> &nbsp;  Vars  :  " << vars  << endl
		<< "<br> &nbsp;  Prop. :  " << props    << endl;	
      }

//      if(updatetheta)
// 	{
// 	  htmllog << "<p><A><img BORDER=0 SRC=\"meanhist.png\"></A><p>" << endl;
// 	}

      if(!nonspatial)
	{
	  // output MRF precision convergence plot
	  htmllog << "<p><A><img BORDER=0 SRC=\"mrfprechist.png\"></A><p>" << endl;
	  htmllog << "</CENTER> Spatial mixture modelling based on: <br> Mixture Models with Adaptive Spatial Regularisation for Segmentation with an Application to FMRI Data; Woolrich, M., Behrens, T., Beckmann, C., and Smith, S.; IEEE TMI In Press 2004.";

	}
      else
	{
	  htmllog << "</CENTER> Non-spatial mixture modelling based on: <br> Mixture Models with Adaptive Spatial Regularisation for Segmentation with an Application to FMRI Data; Woolrich, M., Behrens, T., Beckmann, C., and Smith, S.; IEEE TMI In Press 2004.";
	
	}
 
      { 
	htmllog<< "<HR><FONT SIZE=1>This page produced automatically by "
	       << "mm </A>" 
	       << " - a part of <A HREF=\"http://www.fmrib.ox.ac.uk/fsl\">FSL - "
	       << "FMRIB Software Library</A>.</FONT>" << endl
	       << "</BODY></HTML>" << endl;
      } 
  }

}
