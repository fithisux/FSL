/*  FAST4 - FMRIB's Automated Segmentation Tool v4

    John Vickers, Mark Jenkinson, Steve Smith and Matthew Webster
    FMRIB Image Analysis Group

    Copyright (C) 2005-2009 University of Oxford  */

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

#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "utils/options.h"
#include <vector>

class ZMRISegmentation
{
 public:
   ZMRISegmentation();
  ~ZMRISegmentation()
    {
    }
  void TanakaCreate(const NEWIMAGE::volume<float>& image, float fbeta, int nclasses, float nblowpass, bool bbias, int biterationspve, float mixeltypeMRF, int nbiter, int initinitfixed, int winitfixed, int bapused, float hyp, bool verb,string mansegfle,int typeoffile);
  int TanakaMain(NEWIMAGE::volume<float>& pcsf, NEWIMAGE::volume<float>& pgm, NEWIMAGE::volume<float>& pwm);

  NEWIMAGE::volume4D<float> members;
  NEWIMAGE::volume4D<float> m_post;
  NEWIMAGE::volume4D<float> m_pve;
  NEWIMAGE::volume<float> m_BiasField;
  NEWIMAGE::volume<int> m_Segment;
  NEWIMAGE::volume<int> m_pveSegment;
  NEWIMAGE::volume<int>hardPV;


private:
  void BiasRemoval();
  void Classification();
  void Classification(int x, int y, int z);
  NEWIMAGE::volume<float> Convolve(NEWIMAGE::volume<float>& resfieldimage);
  void Dimensions();
  void Initclass();
  float MRFWeightsTotal();
  double MRFWeightsInner(const int x, const int y, const int z,const int c);
  double MRFWeightsAM(const int l, const int m, const int n);
  float logGaussian(float y_i, float mu, float sigma);
  float PVEnergy(int x, int y, int z, float mu, float sigma);
  void qsort();
  void TanakaIterations();
  void TanakaHyper();
  void TanakaPriorHyper();
  void Initialise();
  void InitSimple(const NEWIMAGE::volume<float>& pcsf, const NEWIMAGE::volume<float>& pgm, const NEWIMAGE::volume<float>& pwm);
  void UpdateMembers(NEWIMAGE::volume4D<float>& probability);
  void WeightedKMeans();
  void pveIterations();
  void pveInitialize();
  void takeexpo();
  void InitKernel(float kernalsize);
  void pveClassification(int x, int y, int z);
  void pveClassification();
  void Probabilities(int index);
  void calculateVolumes(const NEWIMAGE::volume4D<float>& probs);
  void printVolumeTotals();
  void PVClassificationStep();
  void ICMPV();
  void PVEnergyInit();
  void PVestimation();
  void MeansVariances(int numberofclasses, NEWIMAGE::volume4D<float>& probability);

  NEWIMAGE::volume4D<float>PVprob;
  NEWIMAGE::volume4D<float> pvprobsinit;
  NEWIMAGE::volume4D<float> talpriors;
  NEWIMAGE:: volume4D<float> m_prob;
  NEWIMAGE::volume <float>  m_Mricopy;
  NEWIMAGE::volume<float> m_Mri;
  NEWIMAGE::volume<float> m_mask;
  ColumnVector kernelx, kernely, kernelz;
  vector<double> volumequant;
  vector<float> m_mean;
  vector<float> m_variance;
  vector<float> weight;
  vector<float> rhs;
  double m_nxdim, m_nydim, m_nzdim;
  float m_nbLowpass;
  float beta;
  float pveBmixeltype;
  float nvoxel;
  float Hyper;
  float amx, amy, amz, amxy, amzy, amzx;
  int imagetype,m_nWidth, m_nHeight, m_nDepth;
  int neighbour[27];
  int m_nbIter;
  int bapusedflag;
  int nClasses;
  int inititerations;
  int usingmeanfield;
  int iterationspve;
  int initfixed;
  bool biasfieldremoval;
  bool verboseusage;
  string mansegfile;

  class kmeansexception: public exception
  {
  public:
    virtual const char* what() const throw ()
    {
      return "Exception: Not enough classes detected to init KMeans";
    }
  } kmeansexc;

};
