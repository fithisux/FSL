/*  FAST4 - FMRIB's Automated Segmentation Tool v4

    John Vickers, Mark Jenkinson and Steve Smith
    FMRIB Image Analysis Group

    Copyright (C) 2005-2007 University of Oxford  */

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


using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;
using namespace std;

class ZMRIMULTISegmentation
{
 public:

  ZMRIMULTISegmentation();
  ~ZMRIMULTISegmentation()
    {
    }
  void TanakaCreate(const NEWIMAGE::volume<float>* images, int nclasses, bool b2d, int nbiter,float nblowpass, float fbeta, int bapused, bool pveboolean, int nchan, bool bbias, int initfixity, bool verb,int biterationspve, int winit, float mixelpB, float Hyp,string mansegfle,int finalT2file);
  int TanakaMain(NEWIMAGE::volume<float>& pcsf, NEWIMAGE::volume<float>& pgm, NEWIMAGE::volume<float>& pwm);

 
  volume4D<float> m_post;
  volume4D<float> members;
  volume4D<float> m_pve;
  volume<float>* m_Finalbias;
  volume<int> m_Segment;
  volume<int> m_pveSegment;
  volume<int> hardPV;
  Matrix coord_minus_mean;
  float maximum, minimum;


private:
  NEWIMAGE::volume<float> Convolve(NEWIMAGE::volume<float>& resfieldimage);
  NEWIMAGE::volume4D<float> InitclassAlt(int);
  NEWIMAGE::volume4D<float>  Initclass(int noclasses);
  Matrix covariancematrix(int classid, volume4D<float> probability);
  float M_2PI(int numberofchan);
  float logpveGaussian(int x, int y, int z, Matrix mu, Matrix sig, float detsig);
  float PVEnergy(int x, int y, int z, Matrix mu, Matrix sigma, float sigdet);
  float logGaussian(int classnumber, int x, int y, int z);
  float pvmeans(int clas);
  float pvvar(int clas);
  void Volumesquant(const NEWIMAGE::volume4D<float>& probs);
  void Initialise();
  void takeexpo();
  float  MRFWeightsTotal();
  double MRFWeightsInner(const int x, const int y, const int z,const int c);
  double MRFWeightsAM(const int l, const int m, const int n);
  void InitWeights();
  void UpdateWeights();
  void InitSimple(const NEWIMAGE::volume<float>& pcsf, const NEWIMAGE::volume<float>& pgm, const NEWIMAGE::volume<float>& pwm);
  void InitAprioriKMeans();
  void BiasRemoval();
  void MeansVariances(int numberofclasses);
  void Dimensions();
  void WeightedKMeans();
  void mresmatrix();
  void UpdateMembers(NEWIMAGE::volume4D<float>& probability);
  void PVClassificationStep();
  void ICMPV();
  void PVMoreSophisticated();
  void PVestimation();
  void PVMeansVar(const NEWIMAGE::volume4D<float>& m_posttemp);
  void Classification(int x, int y, int z);
  void Classification();
  void InitKernel();
  void pveClassification(int x, int y, int z);
  void pveClassification();
  int qsort();
  void TanakaHyper();
  void TanakaPriorHyper();
  void TanakaIterations();
  Matrix*  m_inv_co_variance;
  Matrix* m_co_variance;
  Matrix m_mean;
  volume4D<float> m_prob;
  volume4D<float>PVprob;
  volume4D<float> talpriors;
  volume<float>* m_Mri;
  volume<float>* m_Mricopy;
  volume<float>* p_bias;
  volume<float>* m_meaninvcov;
  volume<float>* m_resmean;
  volume<float> pve_eng;
  volume<float> m_maskc;
  volume<int> m_mask;
  ColumnVector kernelx, kernely, kernelz;
  float amx, amy, amz, amxy, amzx, amzy;
  double* volumequant;
  float* rhs;
  float* weight;
  int imagetype;
  float m_nbLowpass;
  float m_nxdim, m_nydim, m_nzdim;
  float beta, Hyper;
  float pveBmixeltype;
  int noclasses;
  int numberofchannels;
  int m_nSlicesize, m_nWidth, m_nHeight, m_nDepth, m_nbIter, initfixed, inititerations, iterationspve;
  int bapusedflag;
  bool verboseusage;
  bool biasfieldremoval;
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

