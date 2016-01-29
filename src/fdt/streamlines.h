/*    Copyright (C) 2012 University of Oxford  */

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
#include <fstream>
#ifndef EXPOSE_TREACHEROUS
#define EXPOSE_TREACHEROUS
#endif
#include "newimage/newimageall.h"
#include "utils/log.h"
#include "meshclass/meshclass.h"
#include "probtrackxOptions.h"
#include "particle.h"
#include "tractvolsx.h"
#include "miscmaths/SpMat.h" 



using namespace std;
using namespace NEWIMAGE;
using namespace Utilities;
using namespace TRACTVOLSX;
using namespace mesh;
using namespace PARTICLE;

namespace TRACT{
  void read_masks(vector<string>& masks, const string& filename);


  class Streamliner{
    //Everything in DTI space is done INSIDE this class and lower level classes (particle and tractvolsx)
    //This class communicates with higher level classes in Seed voxels.
    //
    probtrackxOptions& opts;
    Log& logger;
    Particle m_part;
    vector<ColumnVector> m_path;
    volume<int> m_mask;


    volume<int> m_skipmask;
    volume<int> m_rubbish;
    volume<int> m_stop;
    volume4D<float> m_prefdir;
    volume4D<float> m_loopcheck;
    vector<volume<float>* > m_waymasks;
    vector<bool> m_passed_flags;
    vector<bool> m_own_waymasks;

    Matrix m_Seeds_to_DTI;
    Matrix m_DTI_to_Seeds;
    volume4D<float> m_Seeds_to_DTI_warp;
    volume4D<float> m_DTI_to_Seeds_warp;
    volume4D<float> m_jacx;
    volume4D<float> m_jacy;
    volume4D<float> m_jacz;
    bool m_IsNonlinXfm;
    Matrix m_rotdir;

    Tractvolsx vols;
    float m_lcrat;
    float m_x_s_init;
    float m_y_s_init;
    float m_z_s_init;

    // Streamliner needs to know about matrix3 
    volume<int>  m_mask3;
    volume<int>  m_beenhere3;
    vector<ColumnVector> m_inmask3;

    // we need this class to know about seed space
    const volume<float>& m_seeds;

  public:
    //Constructors
    Streamliner(const volume<float>&);
    ~Streamliner(){
      for(unsigned int i=0;i<m_waymasks.size();i++)
	if(m_own_waymasks[i]) delete m_waymasks[i];
    }
    void add_waymask(volume<float>& myway,const bool& ownership=false){
      //make sure that the waymask to add will not be deleted before
      //this streamliner goes out of scope!!
      m_waymasks.push_back(&myway);
      m_own_waymasks.push_back(ownership);
      m_passed_flags.push_back(false);
    }
    void pop_waymasks(){
      volume<float>* tmpptr=m_waymasks[m_waymasks.size()-1];
      m_waymasks.pop_back();
      m_passed_flags.pop_back();
      if(m_own_waymasks[m_own_waymasks.size()-1]){
	delete tmpptr;
      }
      m_own_waymasks.pop_back();
      
    }
    void clear_waymasks(){
      // clear all waymasks
      for(unsigned int i=0;i<m_waymasks.size();i++)
	pop_waymasks();
    }

    const Tractvolsx& get_vols() const {return vols;}
    inline int nfibres() const {return vols.nfibres();}
    inline const float get_x_seed() const {return m_x_s_init;}
    inline const float get_y_seed() const {return m_y_s_init;}
    inline const float get_z_seed() const {return m_z_s_init;}
    inline const vector<ColumnVector>& get_path_ref() const{return m_path;}
    inline vector<ColumnVector> get_path() const{return m_path;}
    inline void reset(){
      m_part.reset();
      vols.reset(opts.fibst.value());
      for(unsigned int i=0;i<m_passed_flags.size();i++)
	m_passed_flags[i]=false;
    }
    inline void reverse(){
      m_part.restart_reverse();
    }
    int streamline(const float& x_init,const float& y_init, const float& z_init,const ColumnVector& dim_seeds,const int& fibst,const ColumnVector& dir);

    void rotdir(const ColumnVector& dir,ColumnVector& rotdir,const float& x,const float& y,const float& z);

    const volume<int>& get_stop()const{return m_stop;}


    // matrix3 methods
    void clear_beenhere3(){
      for(unsigned int i=0;i<m_inmask3.size();i++){
	m_beenhere3((int)MISCMATHS::round(float(m_inmask3[i](1))),
		    (int)MISCMATHS::round(float(m_inmask3[i](2))),
		    (int)MISCMATHS::round(float(m_inmask3[i](3))))=0;
      }
    }
    void                  clear_inmask3(){m_inmask3.clear();}
    vector<ColumnVector>& get_inmask3(){return m_inmask3;}
    volume<int>&          get_mask3(){return m_mask3;}
    volume<int>&          get_beenhere3(){return m_beenhere3;}
  };


  class Counter{
    probtrackxOptions& opts;
    Log& logger;
    volume<int> m_prob; 
    volume<int> m_beenhere;
    Matrix m_I;
    vector<ColumnVector> m_path;
    
    vector<ColumnVector> m_seedcounts;
    Matrix m_SeedCountMat;
    int    m_SeedRow;
    int    m_numseeds;

    Matrix m_targetmasks;
    vector<string> m_targetmasknames;
    vector<int> m_targflags;


    volume<int> m_seeds_vol2mat;
    Matrix      m_seeds_mat2vol;
    volume<int> m_targets_vol2mat;
    Matrix      m_targets_mat2vol;

    
    volume<int> m_ConMat;
    volume<int> m_CoordMat;
    int m_Conrow;

    volume<int> m_ConMat2;
    volume<int> m_CoordMat2;
    volume<int> m_CoordMat_tract2;
    volume<int> m_lrmask;
    volume4D<int> m_lookup2;
    volume<int> m_beenhere2;
    int m_Conrow2;
    ColumnVector m_lrdim;

    //volume<int>  m_ConMat3;
    SpMat<int>  *m_ConMat3;    // Use sparse representation to allow for big NxN matrices 
    volume<int>  m_Lookup3;
    Matrix       m_CoordMat3;
    
    const volume<float>& m_seeds;
    ColumnVector m_seedsdim;
    Streamliner& m_stline;
    Streamliner& m_nonconst_stline;
    
    Matrix m_pathlengths;
    //int m_nsamp;

  public:
    Counter(const volume<float>& seeds,Streamliner& stline,int numseeds):opts(probtrackxOptions::getInstance()),
									 logger(LogSingleton::getInstance()),
									 m_numseeds(numseeds),
									 m_seeds(seeds),m_stline(stline),
									 m_nonconst_stline(stline){
      //are they initialised to zero?
      m_beenhere.reinitialize(m_seeds.xsize(),m_seeds.ysize(),m_seeds.zsize());
      m_seedsdim.ReSize(3);
      m_seedsdim << m_seeds.xdim() <<m_seeds.ydim() <<m_seeds.zdim();
      m_I=IdentityMatrix(4);
      
    }
    
    void initialise();
    
    void initialise_path_dist(){
      m_prob.reinitialize(m_seeds.xsize(),m_seeds.ysize(),m_seeds.zsize());
      copybasicproperties(m_seeds,m_prob);
      m_prob=0;
    }
    void initialise_seedcounts();
    
    void initialise_matrix1(); //Need to make sure that initialise_path_dist is run first
    void initialise_matrix2();
    void initialise_matrix3();
    
    void initialise_maskmatrix(){} //not written yet
    
    void store_path(){ m_path=m_stline.get_path();}
    void append_path(){
      for(unsigned int i=0;i<m_stline.get_path_ref().size();i++)
	m_path.push_back(m_stline.get_path_ref()[i]);
    }
    void clear_path(){ m_path.clear(); };

    void count_streamline();
    void count_seed();
    void clear_streamline();
    
    
    void update_pathdist();
    void reset_beenhere();
    
    void reset_prob(){m_prob=0;}
    void update_seedcounts();
    void reset_targetflags(){
      for(unsigned int i=0;i<m_targflags.size();i++) m_targflags[i]=0;
    }
    
    
    void update_matrix1(); //update path_dist after each streamline, only run this after each voxel!!
    
    void update_matrix2_row(); //but run this one every streamline as with the others
    void next_matrix2_row(){m_Conrow2++;}//and then run this after each voxel..
    void reset_beenhere2();

    void update_matrix3();
    void reset_beenhere3();

    void update_maskmatrix(){} //not written yet
    
    void save_total(const int& keeptotal);
    void save_total(const vector<int>& keeptotal);
    void save();
    void save_pathdist();
    void save_pathdist(string add);
    void save_seedcounts();
    void save_matrix1();
    void save_matrix2();
    void save_matrix3();
    void save_maskmatrix(){}//not written yet
    

    inline const Streamliner& get_streamline() const {return m_stline;}
    inline Streamliner& get_nonconst_streamline() const {return m_nonconst_stline;}
    inline const volume<float>& get_seeds() const {return m_seeds;}

    
  };
  
  class Seedmanager{
    probtrackxOptions& opts;
    Log& logger;
    Counter& m_counter;    
    Streamliner& m_stline;
    const volume<float>& m_seeds;
    ColumnVector m_seeddims;
  public:
    Seedmanager(Counter& counter):opts(probtrackxOptions::getInstance()),
				  logger(LogSingleton::getInstance()),
				  m_counter(counter),
				  m_stline(m_counter.get_nonconst_streamline()),
				  m_seeds(m_counter.get_seeds()){
      m_seeddims.ReSize(3);
      m_seeddims<<m_seeds.xdim()<<m_seeds.ydim()<<m_seeds.zdim();
    }
    int run(const float& x,const float& y,const float& z,bool onewayonly=false, int fibst=-1);
    int run(const float& x,const float& y,const float& z,bool onewayonly, int fibst,const ColumnVector& dir);
  };

}
