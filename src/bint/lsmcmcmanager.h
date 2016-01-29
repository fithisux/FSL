/*  lsmcmcmanager.h

    Mark Woolrich - FMRIB Image Analysis Group

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

#if !defined(lsmcmcmanager_h)
#define lsmcmcmanager_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "newimage/newimageall.h"
#include "model.h"
#include "bintoptions.h"

using namespace NEWMAT;
using namespace MISCMATHS;
using namespace NEWIMAGE;

namespace Bint {

  class LSMCMCParameter;
  class LSMCMCPrecParameter;

  class LSMCMCVoxelManager
  {
  public:

    LSMCMCVoxelManager(ForwardModel& pmodel,int pburnin, int pnjumps, int psampleevery, int pupdateproposalevery, float pacceptancerate, int pdebuglevel, bool panalmargprec, int pnsamples) : model(pmodel), burnin(pburnin), njumps(pnjumps), sampleevery(psampleevery), updateproposalevery(pupdateproposalevery), acceptancerate(pacceptancerate), nsamples(pnsamples), nparams(0), sumsquares(0), likelihood(0), debuglevel(pdebuglevel), analmargprec(panalmargprec) {}
  
    virtual ~LSMCMCVoxelManager();

    void run();
    void jump();
    void sample();

    void setupparams(float prec);
    void setdata(const ColumnVector& pdata);

    void calcsumsquares();
    float calclikelihood();

    void restoresumsquares(){ sumsquares = sumsquares_old; }

    void restorelikelihood(){ likelihood = likelihood_old; }

    int getnsamples() const {return nsamples;}
    const vector<float>& getsamples(int p);
    const vector<float>& getprecsamples();
    const string& getparamname(int p);
    int getnparams() const {return nparams;}

    int getdebuglevel() const {return debuglevel;}
    float getlikelihood() const {return likelihood;}
    int getntpts() const {return ntpts;}
    const vector<LSMCMCParameter*>& getmcmcparams() const {return mcmcparams;}

  protected:

    ForwardModel& model;

    vector<LSMCMCParameter*> mcmcparams;

    Parameter* precparam;
    LSMCMCPrecParameter* precmcmcparam;
    GammaPrior* precparamprior;

    int burnin;
    int njumps;
    int sampleevery;
    int updateproposalevery; 
    float acceptancerate;

    int nsamples;
    int ntpts;
    int nparams;
    float sumsquares;
    float likelihood;

    int debuglevel;
    bool analmargprec;

    float sumsquares_old;
    float likelihood_old;
    bool updateprec;
    ColumnVector data;    

  private:

    LSMCMCVoxelManager();
    const LSMCMCVoxelManager& operator=(LSMCMCVoxelManager& par);     
    LSMCMCVoxelManager(LSMCMCVoxelManager& des);
  };

class McmcParameter
  {
  public:
    
    McmcParameter(Parameter& pparam, int pnsamples, int pupdateproposalevery, float pacceptancerate, int pdebuglevel) :       
      param(pparam),
      val(pparam.getinitvalue()),
      naccepted(0),
      nrejected(0),
      proposal_std(pparam.getinitstd()),
      jumpcount(0),
      debuglevel(pdebuglevel),
      updateproposalevery(pupdateproposalevery),
      acceptancerate(pacceptancerate)
    {
      samples.reserve(pnsamples);
    }
    
    // get new energy taking into account that this parameter's value has changed
    virtual float new_energy() = 0;

    // get energy assuming current value has not changed
    virtual float old_energy() = 0;

    // jump has been rejected: restore energies
    virtual void restore_energy() = 0;

    void update_proposal_std() {
      proposal_std *= acceptancerate/((1+nrejected)/float(1+naccepted+nrejected));
      //cout <<param.getname()<<endl;
      //cout<<"Nacc "<<naccepted<<endl;
      //cout<<"Nrej "<<nrejected<<endl;
      naccepted = 0;
      nrejected = 0;
    }

    virtual void setup() = 0;

//      virtual void reset(float pvalue, float pproposal_std)
//      {
//        naccepted = 0;
//        nrejected = 0;
//        proposal_std = pproposal_std;
//        val = pvalue;
//        jumpcount = 0;
//      }

    void jump();
    void sample(){samples.push_back(val);}    

    virtual ~McmcParameter(){}

    const float value() const {return val;}

    float& value() {return val;}

    const vector<float>& getsamples() const {return samples;}

    bool getallowtovary() const {return param.getallowtovary();}
    bool getsave() const {return param.getsave();}
  protected:

    Parameter& param;
    float val;
    int naccepted;
    int nrejected;
    float proposal_std;
    int jumpcount;
    vector<float> samples;
    int debuglevel;
    int updateproposalevery;
    float acceptancerate;

  private:

    McmcParameter();
    const McmcParameter& operator=(McmcParameter& par);     
    McmcParameter(McmcParameter& des);
  };


  class LSMCMCParameter : public McmcParameter
  {
  public:
    LSMCMCParameter(Parameter& pparam,int pnsamples, int pupdateproposalevery, float pacceptancerate,LSMCMCVoxelManager& plsmcmc) : 
      McmcParameter(pparam,pnsamples,pupdateproposalevery,pacceptancerate,plsmcmc.getdebuglevel()), lsmcmc(plsmcmc), prior_energy(0.0), prior_old_energy(0.0)
    {}

    ~LSMCMCParameter(){}
    
    void setup()
    {
      calc_prior();
    }

    float new_energy() 
    {      
      float energy = calc_prior();
      if(energy != float(MAX_EN)) 
	{
	  lsmcmc.calcsumsquares();
	  energy += lsmcmc.calclikelihood();
	}
	
      return energy;
    }

    float old_energy() 
    {
      float energy = prior_energy;
      
      if(energy != float(MAX_EN)) 
	{
	  energy += lsmcmc.getlikelihood();
	}

      return energy;
    }

    void restore_energy()
    {
      restoreprior();
      lsmcmc.restoresumsquares();
      lsmcmc.restorelikelihood();
    }

    float calc_prior() { 

      prior_old_energy = prior_energy; 
      prior_energy = param.getprior().calc_energy(val); 

      if(debuglevel==2)
	{
	  cout << "prior_old_energy=" << prior_old_energy << endl;
	  cout << "prior_energy=" << prior_energy << endl;
	}

      return prior_energy; 
    }  

    void restoreprior() { prior_energy = prior_old_energy; }

  protected:
    LSMCMCVoxelManager& lsmcmc;
    float prior_energy;
    float prior_old_energy;

  private:

    LSMCMCParameter();
    const LSMCMCParameter& operator=(LSMCMCParameter& par);     
    LSMCMCParameter(LSMCMCParameter& des);

  };

  class LSMCMCPrecParameter : public McmcParameter
  {
  public:
    LSMCMCPrecParameter(Parameter& pparam,int pnsamples, int pupdateproposalevery, float pacceptancerate,LSMCMCVoxelManager& plsmcmc) : 
      McmcParameter(pparam,pnsamples,pupdateproposalevery,pacceptancerate,plsmcmc.getdebuglevel()), lsmcmc(plsmcmc), extra_energy(0.0), extra_old_energy(0.0), N(plsmcmc.getntpts()), priormean(pparam.getinitvalue()), impropercount(0)
    {}

    ~LSMCMCPrecParameter(){}
    
    void setup()
    {      
      calc_extra();
    }

//     void reset(float pvalue, float pproposal_std)
//     { 
//       McmcParameter::reset(pvalue,pproposal_std);
//       extra_energy = 0.0;
//       impropercount = 0;
//       calc_extra();
//     }

    float new_energy() 
    {      
      return calc_extra() + lsmcmc.calclikelihood();
    }

    float old_energy() 
    {
      return extra_energy + lsmcmc.getlikelihood();
    }

    void restore_energy()
    {
      restoreextra();      
      lsmcmc.restorelikelihood();
    }

    float calc_extra();

    void restoreextra() { extra_energy = extra_old_energy; }

  protected:

    LSMCMCVoxelManager& lsmcmc;

    float extra_energy;
    float extra_old_energy;

    int N;

    float priormean;

    int impropercount;

 private:

    LSMCMCPrecParameter();
    const LSMCMCPrecParameter& operator=(LSMCMCPrecParameter& par);     
    LSMCMCPrecParameter(LSMCMCPrecParameter& des);

  };

  inline LSMCMCVoxelManager::~LSMCMCVoxelManager() { 
    mcmcparams.clear();
    if(!analmargprec) {
      delete precparam; delete precmcmcparam; delete precparamprior;
    }
  }
  
  inline const vector<float>& LSMCMCVoxelManager::getsamples(int p) {return mcmcparams[p]->getsamples();}
  inline const vector<float>& LSMCMCVoxelManager::getprecsamples() {return precmcmcparam->getsamples();}
  inline const string& LSMCMCVoxelManager::getparamname(int p) {return model.getparam(p).getname();}
  
  inline float LSMCMCVoxelManager::calclikelihood() 
  { 
    // calculates -log(likelihood):

    likelihood_old = likelihood;

    if(!analmargprec)
      {               
	likelihood = precmcmcparam->value()*sumsquares/2.0;
      }
    else
      {
	likelihood = ntpts/2.0*std::log(sumsquares);
      }

    if(debuglevel==2)
      {
	cout << "likelihood_old="<< likelihood_old << endl; 
	cout << "likelihood="<< likelihood << endl; 
      }
    
    return likelihood;   
  }

  class LSMCMCManager
  {
  public:
    // constructor
    LSMCMCManager(int pnjumps, int pnburnin, int psampleevery, int pupdateproposalevery, int pacceptancerate, int pdebuglevel, float pprecin, bool panalmargprec, ForwardModel& pmodel, const Matrix& pdata, const volume4D<float>& pmask) : 
      data(pdata),
      mask(pmask),      
      debuglevel(pdebuglevel),
      precin(pprecin),
      analmargprec(panalmargprec),
      nsamples((pnjumps-pnburnin)/psampleevery),
      voxelmanager(pmodel,pnburnin,pnjumps,psampleevery,pupdateproposalevery,pacceptancerate,pdebuglevel,analmargprec,nsamples),
      model(pmodel)
    {
    }

    LSMCMCManager(BintOptions& opts, ForwardModel& pmodel,const Matrix& pdata, const volume4D<float>& pmask) : 
      data(pdata),
      mask(pmask),      
      debuglevel(opts.debuglevel.value()),
      precin(opts.prec.value()),
      analmargprec(opts.analmargprec.value()),
      nsamples((opts.njumps.value()-opts.burnin.value())/opts.sampleevery.value()),
      voxelmanager(pmodel,opts.burnin.value(),opts.njumps.value(),opts.sampleevery.value(),opts.updateproposalevery.value(),opts.acceptancerate.value(),opts.debuglevel.value(),opts.analmargprec.value(),nsamples),
      model(pmodel)      
    {       
    }

    // load data from file in from file and set up starting values
    void setup();

    void run();

    // saves results in logging directory
    void save();     

    int getntpts() const {return ntpts;}
    int getnvoxels() const {return nvoxels;}
    
    const Matrix& getsamples(int paramnum) const {return samples[paramnum];}
    Matrix& getsamples(int paramnum) {return samples[paramnum];}
   
    // Destructor
    virtual ~LSMCMCManager() {}
 
  protected:     

    Matrix data;
    volume4D<float> mask;

    int ntpts;    
    int nvoxels;
    int nparams;

    vector<Matrix> samples;

    Matrix precsamples;

    vector<string> paramnames;

    int debuglevel;
    float precin;
    bool analmargprec;

    int nsamples;

    LSMCMCVoxelManager voxelmanager;
    
    ForwardModel& model;

  private:

    LSMCMCManager();
    const LSMCMCManager& operator=(LSMCMCManager& par);     
    LSMCMCManager(LSMCMCManager& des);

  };
}

#endif







