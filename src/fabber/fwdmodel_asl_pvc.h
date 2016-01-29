/*  fwdmodel_asl_pvc.h - Partial Volume Correction resting state ASL model (Buxton)

    Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2009 University of Oxford  */

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

#include "fwdmodel.h"
#include "inference.h"
#include <string>
using namespace std;

class ASL_PVC_FwdModel : public FwdModel {
public: 
  // Virtual function overrides
  virtual void Evaluate(const ColumnVector& params, 
			      ColumnVector& result) const;
  static void ModelUsage();
  virtual string ModelVersion() const;
                  
  virtual void DumpParameters(const ColumnVector& vec,
                                const string& indents = "") const;
                                
  virtual void NameParams(vector<string>& names) const;     
  virtual int NumParams() const 
  { return (infertiss?2:0) - (singleti?1:0) + (infertiss?(infertau?1:0):0) + (inferart?2:0) + (infert1?2:0) + (infertaub?1:0)  + (inferwm?(2+(infertau?1:0)+(infert1?1:0)+(usepve?2:0)):0);  

    //return 2 - (singleti?1:0) + (infertau?1:0) + (inferart?2:0) + (infert1?2:0) + (inferinveff?1:0) + (infertrailing?1:0) + (infertaub?1:0); 
  } 

  virtual ~ASL_PVC_FwdModel() { return; }

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

  using FwdModel::SetupARD;
  virtual void SetupARD(const MVNDist& posterior, MVNDist& prior, double& Fard);
  virtual void UpdateARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const;

  // Constructor
  ASL_PVC_FwdModel(ArgsType& args);


protected: // Constants

  // Lookup the starting indices of the parameters
  int tiss_index() const {return (infertiss?1:0);} //main tissue parameters: ftiss and delttiss alway come first

  int tau_index() const {  return (infertiss?2:0) + (infertiss?(infertau?1:0):0);  }

  int art_index() const {  return (infertiss?2:0) + (infertiss?(infertau?1:0):0) + (inferart?1:0); }

  int t1_index() const { return (infertiss?2:0) + (infertiss?(infertau?1:0):0) + (inferart?2:0) + (infert1?1:0); }
  
  //int inveff_index() const { return 2 + (infertau?1:0) + (inferart?2:0) + (infert1?2:0) +(inferinveff?1:0); }

  //int trailing_index() const { return 2 + (infertau?1:0) + (inferart?2:0) + (infert1?2:0) + (infertrailing?1:0); }

  //int taub_index() const { return 2 + (infertau?1:0) + (inferart?2:0) + (infert1?2:0) + (inferinveff?1:0) + (infertrailing?1:0) + (infertaub?1:0);}

  int taub_index() const { return (infertiss?2:0) + (infertiss?(infertau?1:0):0) + (inferart?2:0) + (infert1?2:0) + (infertaub?1:0);}

  //int R_index() const { return 2 + (infertau?1:0) + (inferart?2:0) + (infert1?2:0) + (infertaub?1:0) + (inferart?1:0);}
  int wm_index() const { return (infertiss?2:0) + (infertiss?(infertau?1:0):0) + (inferart?2:0) + (infert1?2:0) + (infertaub?1:0)  + (inferwm?1:0); }

  int pv_index() const { return (infertiss?2:0) + (infertiss?(infertau?1:0):0) + (inferart?2:0) + (infert1?2:0) + (infertaub?1:0)  + (inferwm?(2 + (infertau?1:0) + (infert1?1:0) ):0) + (usepve?1:0); }

  // vector indices for the parameters to expereicne ARD
  vector<int> ard_index;


  // scan parameters
  double seqtau; //bolus length as set by the sequence
  double setdelt;

  int repeats;
  double t1;
  double t1b;
  double t1wm;
  double lambda;
  double pretisat;
  bool grase; //to indicate data was collected with GRASE-ASL
 double slicedt;
 bool casl;

  bool infertiss;
  bool singleti; //specifies that only tissue perfusion should be inferred
  bool infertau;
  bool infertaub;
  bool inferart;
  bool infert1;
  bool inferwm;
  bool usepve;
  //bool inferinveff;
  //bool infertrailing;

  // ard flags
  bool doard;
  bool tissard;
  bool artard;
  bool wmard;

  ColumnVector tis;
  Real timax;


};
