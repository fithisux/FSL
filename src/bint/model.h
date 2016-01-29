/*  model.h

    Mark Woolrich FMRIB Image Analysis Group

    Copyright (C) 2002 University of Oxford  */

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

#if !defined(model_h)
#define model_h

#include <string>
#include <vector>
#include <math.h>

#include "utils/tracer_plus.h"
#include "miscmaths/miscmaths.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef MAX_EN
#define MAX_EN 1e16
#endif

using namespace Utilities;
using namespace MISCMATHS;

namespace Bint {
    
  class Prior
  {
  public:
    
    Prior(){}
    
    virtual ~Prior(){}
    virtual float calc_energy(float value) const  = 0;
    virtual float calc_gradient(float value) const  = 0;

  private:
    
  };
  
  class GaussPrior : public Prior
  {
  public:
    GaussPrior(float pmn,float pstd, float pmin = -1e16, float pmax = 1e16) : 
      Prior(),
      mn(pmn),std(pstd),min(pmin),max(pmax) {}

    const GaussPrior& operator=(const GaussPrior& pin)
    {
      mn = pin.mn;
      std = pin.std;
      min = pin.min;
      max = pin.max;

      return *this;
    }

    GaussPrior(const GaussPrior& pin)
    {
      operator=(pin);
    }
    
    virtual ~GaussPrior(){}

    virtual float calc_energy(float value) const { 
      float energy = MAX_EN;
      if(value > min && value < max)
	energy = pow(double((value-mn)/std),2.0)/2.0;
      
      return energy; }

    virtual float calc_gradient(float value) const {
      float grad = 0;
      
      if(value > min && value < max)
	grad = ((mn-value)/(std*std));
      
      return grad;
    }

  private:

    GaussPrior();    
    float mn;
    float std;
    float min;
    float max;
  };
  
  
  class UnifPrior : public Prior
  {
  public:
    UnifPrior(float pmin = -1e16, float pmax = 1e16) : 
      Prior(),
      min(pmin),max(pmax) {
      
      en=-std::log( 1/(max-min) );
    }
    
    const UnifPrior& operator=(const UnifPrior& pin)
    {
      min = pin.min;
      max = pin.max;
      en=-std::log( 1/(max-min) );
      return *this;
    }

    UnifPrior(const UnifPrior& pin)
    {
      operator=(pin);
    }

    virtual ~UnifPrior(){}
    

    virtual float calc_energy(float value) const { 
      float energy=MAX_EN;
      if(value >= min && value <= max )
	energy=en;
      return energy; }

    virtual float calc_gradient(float value) const {
      float grad = 0;
     
      return grad;
    }

  private:

    UnifPrior();    
    float min;
    float max;
    float en;
  };


  class GammaPrior : public Prior
  {
  public:
    // m=a/b, v=a/b^2, a=m^2/v, b=m/v
    GammaPrior(float pa,float pb, float pmin = 0, float pmax = 1e16) : 
      Prior(),
      a(pa),b(pb),min(pmin),max(pmax) {}

    virtual ~GammaPrior(){}

    const GammaPrior& operator=(const GammaPrior& pin)
    {
      a = pin.a;
      b = pin.b;
      min = pin.min;
      max = pin.max;

      return *this;
    }

    GammaPrior(const GammaPrior& pin)
    {
      operator=(pin);
    }

    virtual float calc_energy(float value) const { 
     
      float energy = MAX_EN;
      if(value > min && value < max)
	{
	  energy = -(a-1)*std::log(value)+b*value;
	}

      return energy; 
    }

    virtual float calc_gradient(float value) const {
      float grad = 0;
      if(value > min && value < max)
	{
	  grad = -(a-1)/value+b;
	}

      return grad;
    }

  private:

    GammaPrior();

    float a;
    float b;
    float min;
    float max;
  };



  class SinPrior : public Prior
  {
  public:
    SinPrior(float pscale=1,float pmin=-1e16,float pmax=1e16) : 
      Prior(),
      scale(pscale),min(pmin),max(pmax){}
    
    const SinPrior& operator=(const SinPrior& pin)
    {
      scale = pin.scale;
      min=pin.min;
      max=pin.max;
      return *this;
    }

    SinPrior(const SinPrior& pin)
    {
      operator=(pin);
    }

    virtual ~SinPrior(){}
    

    virtual float calc_energy(float value) const { 
      
      float energy= MAX_EN;
      
      if(value > min && value < max){
	if(!value==0){
	  energy=-std::log(fabs(std::sin(value/scale)/2));
	}
      }
      
      return energy; 
    }
    
    virtual float calc_gradient(float value) const {
      float grad = 0;
      if(value > min && value < max)
	{
	  grad = sign(std::sin(value/scale))*(std::cos(value/scale)/scale)/(std::sin(value/scale));

	}

      return grad;
    }

  private:

    SinPrior();    
    float scale;
    float min;
    float max;
  };

 class GaussARDPrior : public Prior
  {
  public:
    GaussARDPrior() : 
      Prior(){}
    
    const GaussARDPrior& operator=(const GaussARDPrior& pin)
    {
      return *this;
    }

    GaussARDPrior(const GaussARDPrior& pin)
    {
      operator=(pin);
    }

    virtual ~GaussARDPrior(){}
    

    virtual float calc_energy(float value) const { 
      
      float energy= MAX_EN;
      if(value!=0){
	energy=std::log(fabs(value));
      }
      
      return energy; 
    }
    
    virtual float calc_gradient(float value) const {
      float grad = 1; //Wooly thinks this is better than 0
      if(value!=0){ grad =sign(value)/fabs(value);}
      return grad;
    }

  private:
 
 };



  class Parameter
  {
  public:
        
    Parameter(const string& pname, float pinitvalue, float pinitstd, Prior& pprior, bool pallowtovary = true, bool psave = true) :
      name(pname),
      init_value(pinitvalue),
      init_std(pinitstd),
      priorobj(pprior),
      allowtovary(pallowtovary),
      save(psave)
    {      
    }

    virtual ~Parameter(){}
    const string& getname() const {return name;}
    const Prior& getprior() const {return priorobj;}
    float getinitvalue() const {return init_value;}
    float getinitstd() const {return init_std;}
    bool getallowtovary() const {return allowtovary;}
    bool getsave() const { return save; }
    
    void setinitvalue(float pinitvalue) {init_value = pinitvalue;}

  protected:

    string name;
    float init_value;    
    float init_std;    
    Prior& priorobj;
    bool allowtovary;
    bool save;
  
  private:

    Parameter();
    const Parameter& operator=(Parameter& par);
    Parameter(const Parameter&);
  };
  
  class ForwardModel
  {
  public:
    ForwardModel(int pdebuglevel) : 
      debuglevel(pdebuglevel),
      paramcount(0)
    {}

    virtual ~ForwardModel(){params.clear();priors.clear();}

    virtual  ReturnMatrix nonlinearfunc(const ColumnVector& paramvalues) const = 0;

    virtual void setparams() = 0;
    virtual void initialise(const ColumnVector& data) = 0;
    // indexes from zero:
    Parameter& getparam(int p) {return *params[p];}

    int getnparams() const {return paramcount;}
  
    void clear_params() {params.clear();paramcount = 0;}

    void add_param(const string& pname, float pinit_value, float pinit_std, Prior* tmp, bool pallowtovary, bool psave)
    {
      paramcount++;

      priors.push_back(tmp);
      params.push_back(new Parameter(pname,pinit_value,pinit_std,*tmp,pallowtovary,psave));
    }   

    void add_param(const string& pname, float pinit_value, float pinit_std, GaussPrior& pprior, bool pallowtovary = true, bool psave = true)
    {   
      add_param(pname,pinit_value,pinit_std,new GaussPrior(pprior),pallowtovary,psave);
    }  
    
    void add_param(const string& pname, float pinit_value, float pinit_std, UnifPrior& pprior, bool pallowtovary = true, bool psave = true)
    {
      add_param(pname,pinit_value,pinit_std,new UnifPrior(pprior),pallowtovary,psave);
    }  

    void add_param(const string& pname, float pinit_value, float pinit_std, GammaPrior& pprior, bool pallowtovary = true, bool psave = true)
    {
      add_param(pname,pinit_value,pinit_std,new GammaPrior(pprior),pallowtovary,psave);
    }  

    void add_param(const string& pname, float pinit_value, float pinit_std, SinPrior& pprior, bool pallowtovary = true, bool psave = true)
    {
      add_param(pname,pinit_value,pinit_std,new SinPrior(pprior),pallowtovary,psave);
    }  
  
    void add_param(const string& pname, float pinit_value, float pinit_std, GaussARDPrior& pprior, bool pallowtovary = true, bool psave = true)
    {
      add_param(pname,pinit_value,pinit_std,new GaussARDPrior(pprior),pallowtovary,psave);
    } 
  protected:
    
    int debuglevel;   
    
    vector<Parameter*> params;

    int paramcount;

    vector<Prior*> priors;

  private:

    ForwardModel();
    const ForwardModel& operator=(ForwardModel& par);
    ForwardModel(const ForwardModel&);
  };

  class gForwardModel : public ForwardModel
  {
  public:
    gForwardModel(int pdebuglevel) : 
      ForwardModel(pdebuglevel)
    {}

    virtual ~gForwardModel(){}

    // gradient
    virtual  ReturnMatrix gradient(const ColumnVector& paramvalues) const = 0;

  private:

    gForwardModel();
    const gForwardModel& operator=(gForwardModel& par);
    gForwardModel(const gForwardModel&);
  };
}
   
#endif

