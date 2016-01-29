/*  mixture_model.h

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 2004 University of Oxford  */

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

#if !defined(mixture_model_h)
#define mixture_model_h

#include <iostream>
#include <fstream>
#include "newmatap.h"
#include "newmatio.h"
#include "newimage/newimageall.h"
#include "utils/tracer_plus.h"
#include "mmoptions.h"
#include "connected_offset.h"
#include "miscmaths/sparse_matrix.h"
#include "miscmaths/minimize.h"
#include "libprob.h"

using namespace NEWMAT;
using namespace NEWIMAGE;
using namespace MISCMATHS;

namespace Mm{

  inline double boundexp(double x)
  {
    //    OUT(x);

    double bound=700;
    if(x>bound) x=bound;
    else if(x<-bound) x=-bound;
    double ret = std::exp(x);

    //    OUT(ret);
    return ret;
  }

  class Distribution
  {
  public:
    Distribution() : useprop(false){}
    virtual float pdf(float val) const = 0;
    virtual float dpdfdmn(float val) const = 0;
    virtual float dpdfdvar(float val) const  = 0;

    virtual ~Distribution(){}

    float getmean() const {return mn;}
    float getvar() const {return var;}
    float getprop() const {return prop;}

    virtual bool setparams(float pmn, float pvar, float pprop){mn=pmn;var=pvar;prop=pprop;return true;}    
    void setuseprop(bool puseprop) {useprop = puseprop;}

  protected:
    float mn;
    float var;
    float prop;
    float useprop;   
  };

  class GaussianDistribution : public Distribution
  {
  public:

    GaussianDistribution() : Distribution() {
    }

    virtual float pdf(float val) const {
      float ret = premult*std::exp(-0.5/var*Sqr(val-mn));
      if(useprop) ret *= prop;
      return ret;
    }

    virtual float dpdfdmn(float val) const {      
      float ret = premult*(val-mn)/var*std::exp(-0.5/var*Sqr(val-mn));
      return ret;
    }

    virtual float dpdfdvar(float val) const {
      float ret = premult*0.5*(Sqr(val-mn)-var)/std::pow(var,2)*std::exp(-0.5/var*Sqr(val-mn));
      
      return ret;
    }

    virtual ~GaussianDistribution(){}

    virtual bool setparams(float pmean, float pvar, float pprop) 
    {
      Distribution::setparams(pmean, pvar, pprop);

      if(pvar<=0)
	return false;

      useprop ? premult=pprop/std::sqrt(2*M_PI*pvar) : premult=1.0/std::sqrt(2*M_PI*pvar);

      return true;
    }

  private:
    float premult;

  };

  class GammaDistribution : public Distribution
  {
  public:
    GammaDistribution(float pminmode = 0) : Distribution(), minmode(pminmode) {}

    virtual float pdf(float val) const 
    {
      float ret = 1e-32;
      if(val > 0){
	ret = boundexp(preadd + (a-1) * std::log(val) - b*val);
	if(useprop) ret *= prop;
      }

      return ret;
    }

    virtual float dpdfdmn(float val) const {

      float ret = 0;

      if(val > 0)
	ret = dpdfda(val)*(2*mn/var)+dpdfdb(val)*(1/var);

      return ret;
    }

    virtual float dpdfdvar(float val) const {
      
      float ret = 0;

      if(val > 0)
	ret = dpdfda(val)*(-Sqr(mn)/Sqr(var))+dpdfdb(val)*(-mn/Sqr(var));

      return ret;
    }

    virtual ~GammaDistribution(){}

    virtual bool setparams(float pmn, float pvar, float pprop) {
      Distribution::setparams(pmn, pvar, pprop);
      bool ret = validate();

      a = std::pow(mn,2)/var;
      b = mn/var;

      useprop ? preadd= log(prop)+a*std::log(b)-lgam(a) : preadd=a*std::log(b)-lgam(a);
      useprop ? premult=pprop/std::sqrt(2*M_PI*pvar) : premult=1.0/std::sqrt(2*M_PI*pvar);
      digama = digamma(a);

      if(!ret) OUT("invalid gamma");
      return ret;
    }

    void setminmode(float pminmode) {minmode=pminmode; validate();}

  private:

    float dpdfdb(float val) const { 
      return std::pow(b,a-1)*std::pow(val,a-1)/std::exp(lgam(a))*std::exp(-b*val)*(a-b*val);
    }

    float dpdfda(float val) const { 
      return std::pow(b,a)*std::pow(val,a-1)/std::exp(lgam(a))*std::exp(-b*val)*(std::log(b)+std::log(val)-digama);
    }

    bool validate();

    float digama;
    float preadd;
    float a;
    float b;
    float minmode;    
    float premult;

  };

  class FlippedGammaDistribution : public Distribution
  {
  public:
    FlippedGammaDistribution(float pminmode = 0) : Distribution(), minmode(pminmode) {}

    virtual float pdf(float val) const
    {
      float ret = 1e-32;
      val = -val;

      if(val > 0) {
	ret = boundexp(preadd + (a-1) * std::log(val) - b*val);
	//ret = boundexp(preadd + (a-1) * std::log(std::abs(val)) - b*std::abs(val));

	if(useprop) ret *= prop;
      }
      return ret;
    }
    
    virtual float dpdfdmn(float val) const {

      // flip val
      val = -val;
      float pmn = -mn;
 
      float ret = 0;
      if(val > 0)
	ret = dpdfda(val)*(2*pmn/var)+dpdfdb(val)*(1/var);
      
      return -ret;
    }

    virtual float dpdfdvar(float val) const {

      // flip val
      val = -val;
      float pmn = -mn;

      float ret = 0;

      if(val > 0)
	ret = dpdfda(val)*(-Sqr(pmn)/Sqr(var))+dpdfdb(val)*(-pmn/Sqr(var));

      return ret;
    }

    virtual ~FlippedGammaDistribution(){}

    virtual bool setparams(float pmn, float pvar, float pprop)
    {
      Distribution::setparams(pmn, pvar, pprop);   
      bool ret = validate();

      a = std::pow(mn,2)/var;      
      b = -(mn)/var;
//        OUT(a);
//        OUT(b);
//        OUT(prop);
      useprop ? preadd= log(prop)+a*std::log(b)-lgam(a) : preadd=a*std::log(b)-lgam(a);
      useprop ? premult=pprop/std::sqrt(2*M_PI*pvar) : premult=1.0/std::sqrt(2*M_PI*pvar);
      digama = digamma(a);

      if(!ret) OUT("invalid gamma");
      return ret;
    }

    void setminmode(float pminmode) {minmode=pminmode; validate();}

  private:
    bool validate();

    float dpdfdb(float val) const { 
      return std::pow(b,a-1)*std::pow(val,a-1)/std::exp(lgam(a))*std::exp(-b*val)*(a-b*val);
    }

    float dpdfda(float val) const { 
      return std::pow(b,a)*std::pow(val,a-1)/std::exp(lgam(a))*std::exp(-b*val)*(std::log(b)+std::log(val)-digama);
    }

    float digama;
    float preadd;
    float a;
    float b;
    float minmode;
    float premult;
  };
 
  class SmmVoxelFunction : public EvalFunction
  {
  public:
    SmmVoxelFunction(float pdata, vector<Distribution*>& pdists, float plambda, float plog_bound)
      : EvalFunction(),
	data(pdata),
	dists(pdists),
	nclasses(pdists.size()),
	lambda(plambda),
	log_bound(plog_bound)
    {}

    float evaluate(const ColumnVector& x) const; //evaluate the function
    
    virtual ~SmmVoxelFunction(){};

  private:
    SmmVoxelFunction();
    const SmmVoxelFunction& operator=(SmmVoxelFunction& par);
    SmmVoxelFunction(const SmmVoxelFunction&);

    float data;
    vector<Distribution*>& dists;
    int nclasses;
    float lambda;
    float log_bound;
  };

  class SmmFunction : public gEvalFunction
  {
  public:
    SmmFunction(const ColumnVector& pdata, vector<Distribution*>& pdists, const float& pmrf_precision, const volume<int>& pmask, const vector<Connected_Offset>& pconnected_offsets, const volume<int>& pindices, const SparseMatrix& pD, float plambda, float plog_bound);

    float evaluate(const ColumnVector& x) const; //evaluate the function
    
    ReturnMatrix g_evaluate(const ColumnVector& x) const; //evaluate the gradient function
    
    virtual ~SmmFunction(){};

  private:
    SmmFunction();
    const SmmFunction& operator=(SmmFunction& par);
    SmmFunction(const SmmFunction&);

    const ColumnVector& data; 
    vector<Distribution*>& dists;
    const float& mrf_precision;
    const volume<int>& mask;
    const vector<Connected_Offset>& connected_offsets;
    const volume<int>& indices;
    const SparseMatrix& D;

    int num_superthreshold;
    int nclasses;
    float lambda;
    float log_bound;
    
  }; 

  class SmmFunctionDists : public gEvalFunction
  //class SmmFunctionDists : public EvalFunction
  {
  public:
    SmmFunctionDists(const ColumnVector& pdata, vector<Distribution*>& pdists, const float& pmrf_precision, const volume<int>& pmask, const vector<Connected_Offset>& pconnected_offsets, const volume<int>& pindices, float plambda, float plog_bound, const ColumnVector& m_tildew);

    float evaluate(const ColumnVector& x) const; //evaluate the function
    
    ReturnMatrix g_evaluate(const ColumnVector& x) const; //evaluate the gradient function
    
    virtual ~SmmFunctionDists(){};

  private:
    SmmFunctionDists();
    const SmmFunctionDists& operator=(SmmFunctionDists& par);
    SmmFunctionDists(const SmmFunctionDists&);

    const ColumnVector& data; 
    vector<Distribution*>& dists;
    const float& mrf_precision;
    const volume<int>& mask;
    const vector<Connected_Offset>& connected_offsets;
    const volume<int>& indices;
    vector<RowVector> w;

    int num_superthreshold;
    int nclasses;
    float lambda;
    float log_bound;
    
    const ColumnVector& m_tildew;
  };

  class Mixture_Model
    {
    public:

      // Constructor
      Mixture_Model(const volume<float>& pspatial_data, const volume<int>& pmask, const volume<float>& pepi_example_data, float pepibt, vector<Distribution*>& pdists, vector<volume<float> >& pw_means, ColumnVector& pY, MmOptions& popts);

      Mixture_Model(const volume<float>& pspatial_data, const volume<int>& pmask, const volume<float>& pepi_example_data, float pepibt, vector<Distribution*>& pdists, vector<volume<float> >& pw_means, ColumnVector& pY, bool pnonspatial=false, int pniters=10, bool pupdatetheta=true, int pdebuglevel=0, float pphi=0.015, float pmrfprecstart=10.0, int pntracesamps=10, float pmrfprecmultiplier=10.0, float pinitmultiplier=6.0, bool pfixmrfprec=false);

      // setup
      void setup();

      // run
      void run();

      // save data to logger dir
      void save() ;
       
      // Destructor
      virtual ~Mixture_Model(){}

    private:
    
      Mixture_Model();
      const Mixture_Model& operator=(Mixture_Model&);     
      Mixture_Model(Mixture_Model&);
 
      void update_theta();
      void update_mrf_precision();
      void update_tildew_scg();
      void update_voxel_tildew_vb();
      void calculate_taylor_lik();     
      void calculate_trace_tildew_D();

      void get_weights(vector<ColumnVector>& weights, const ColumnVector& pmtildew);

      void get_weights2(vector<ColumnVector>& weights, vector<vector<vector<float> > >& weights_samps, vector<vector<vector<float> > >& tildew_samps, int nsamps, const ColumnVector& pmtildew);

      void save_weights(const ColumnVector& pmtildew, const char* affix, bool usesamples = true);

      int xsize; 
      int ysize;
      int zsize;
      int num_superthreshold;
      int nclasses;

      const volume<float>& spatial_data;
      const volume<int>& mask;
      const volume<float>& epi_example_data;
      float epibt;

      volume4D<float> localweights;      
      vector<Connected_Offset> connected_offsets;

      volume<int> indices;

      ColumnVector& Y;
      SparseMatrix D;

      ColumnVector m_tildew;
      vector<SymmetricMatrix> prec_tildew;
      vector<SymmetricMatrix> cov_tildew;

      SparseMatrix precision_lik;
      ColumnVector derivative_lik;

      float mrf_precision;
      //      float mrf_precision_old;

      bool nonspatial;
      int niters;
      bool stopearly;

      bool updatetheta;
      int debuglevel;

      // logistic transform params:
      float lambda; 
      float log_bound;

      float trace_covariance_tildew_D;

      int it;

      vector<Distribution*>& dists;
      vector<volume<float> >& w_means;

      int ntracesamps;
      float mrfprecmultiplier;
      float initmultiplier;
      bool fixmrfprec;

      float trace_tol;
      float scg_tol;

      vector<float> meanhist;
      vector<float> mrf_precision_hist;
    }; 
  
  ReturnMatrix sum_transform(const RowVector& wtilde, float log_bound);
  ReturnMatrix logistic_transform(const RowVector& wtilde,float lambda,float log_bound);
  ReturnMatrix inv_transform(const RowVector& w,float lambda,float log_bound,float initmultiplier);
  
  void ggmfit(const RowVector& data, vector<Distribution*>& pdists, bool useprops);
  void plot_ggm(const vector<volume<float> >& w_means, const vector<Distribution*>& dists, const volume<int>& mask, const ColumnVector& Y);

  void make_ggmreport(const vector<volume<float> >& w_means, const vector<Distribution*>& dists, const volume<int>& mask, const volume<float>& spatial_data, bool zfstatmode, bool overlay, const volume<float>& epivol, float thresh, bool nonspatial, bool updatetheta, const string& data_name);
  void calculate_props(const vector<volume<float> >& w_means, vector<Distribution*>& dists, const volume<int>& mask);
}
#endif
