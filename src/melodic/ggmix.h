/*  MELODIC - Multivariate exploratory linear optimized decomposition into 
              independent components
    
    ggmix.h - class for Gaussian/Gamma Mixture Model

    Christian F. Beckmann, FMRIB Analysis Group
    
    Copyright (C) 1999-2013 University of Oxford */

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

#ifndef __GGMIX_h
#define __GGMIX_h

#include "newimage/newimageall.h"
/*#include "utils/log.h"
#include "melodic.h"
#include "utils/options.h"
#include "meloptions.h"
*/

//using namespace Utilities;
using namespace NEWIMAGE;

namespace GGMIX{
  
  class ggmix
    {
    public:
    
      ggmix(){}
      //MelodicOptions &popts, Log &plogger):
      //	opts(popts),
      //	logger(plogger)
      //	{
      //	}

      ~ggmix() { 
	  }

      void save();

      void setup(const RowVector& dat, const string dirname,
		 int here, volume<float> themask, 
		 volume<float> themean, int num_mix = 3, 
		 float eps = 0.0, bool fixdim = false);
      
      void gmmfit();
      void ggmfit();

      inline void fit(string mtype = string("GGM"))
	{
	  mmtype = mtype;
	  if(mmtype==string("GGM")) 
	    this->ggmfit(); 
	  else 
	    this->gmmfit();

	  //re-insert mean and stdev

	  data = data*datastdev + datamean;
	  //threshmaps = threshmaps*datastdev + datamean;
	  means = means*datastdev + datamean;
	  vars = vars*datastdev*datastdev;
	}

      inline Matrix threshold(string levels)
	{return this->threshold(data, levels);}
      inline Matrix threshold(RowVector& levels) 
	{return this->threshold(data, levels);}
      Matrix threshold(const RowVector& dat, Matrix& levels);
      Matrix threshold(const RowVector& dat, string levels);

      void status(const string &txt);
 
      inline RowVector& get_means() {return means;}
      inline void set_means(RowVector& Arg) {means = Arg;}
    
      inline RowVector& get_vars() {return vars;}
      inline void set_vars(RowVector& Arg) {vars = Arg;}
      
      inline RowVector& get_pi() {return props;}
      inline void set_pi(RowVector& Arg) {props = Arg;}
      
      inline RowVector& get_data() {return data;}
      inline void set_data(RowVector& Arg) {data = Arg;}

      inline RowVector& get_prob() {return probmap;}

      inline float get_eps() {return epsilon;}
      inline void set_eps(float Arg) {epsilon = Arg;}

      inline Matrix& get_threshmaps() {return threshmaps;}
      inline void set_threshmaps(Matrix& Arg) {threshmaps = Arg;}

      inline bool isfitted(){return fitted;}

      inline int mixtures(){return nummix;}
     
      inline string get_type() { return mmtype;}
      inline void set_type(string Arg) { mmtype = Arg;}
      
      inline string get_prefix() { return prefix;}
      inline void  set_prefix(string Arg) { prefix = Arg;}

      inline RowVector get_probmap() {return probmap;}

      inline float get_offset() {return offset;}
      inline void set_offset(float Arg) {offset = Arg;}

      inline void flipres(int num){
	means = -means;
	data = -data;
	threshmaps = -threshmaps;
	if(mmtype=="GGM"){
	  float tmp;
	  tmp= means(2);means(2)=means(3);means(3)=tmp;
	  tmp=vars(2);vars(2)=vars(3);vars(3)=tmp;
	  tmp=props(2);props(2)=props(3);props(3)=tmp;
	}
      }      

      inline void add_infstr(string what){
	threshinfo.push_back(what);
      }

      inline string get_infstr(int num){
	if((threshinfo.size()<(unsigned int)(num-1))||(num<1))
	  return string("");
	else
	  return threshinfo[num-1];
      }

      inline int size_infstr(){
      return threshinfo.size();
      }

      inline void clear_infstr(){
	threshinfo.clear();
      }

      inline void smooth_probs(float howmuch){
	volume4D<float> tempVol;
	tempVol.setmatrix(probmap,Mask);
        tempVol[0]= smooth(tempVol[0],howmuch);
        probmap = tempVol.matrix(Mask);
      }



      double datamean;
      double datastdev;

    private:
      //    MelodicOptions &opts;     
      // Log &logger; //global log file

      //Log mainhtml;

      void gmmupdate();
      float gmmevidence();
      void gmmreducemm();
      void add_params(Matrix& mu, Matrix& sig, Matrix& pi, 
		      float logLH, float MDL, float Evi, bool advance = false);
      void get_params(int index, Matrix& mu, Matrix& sig, Matrix& pi, 
		      float logLH, float MDL, float Evi);

      Matrix Params;
      Matrix threshmaps;

      RowVector means;
      RowVector vars;
      RowVector props;
      RowVector data;
      RowVector probmap;

      volume<float> Mean;
      volume<float> Mask;

      float epsilon;
      float logprobY;
      float MDL;
      float Evi;
      float offset;

      int nummix;
      int numdata;
      int cnumber;

      bool fitted;
      bool fixdim;

      string prefix;
      string mmtype;
      string dirname;

      vector<string> threshinfo;

    };
}

#endif
