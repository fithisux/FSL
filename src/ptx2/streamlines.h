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
#include "csv.h"
#include "utils/tracer_plus.h"

using namespace std;
using namespace NEWIMAGE;
using namespace Utilities;
using namespace TRACTVOLSX;
using namespace mesh;
using namespace PARTICLE;

struct infoVertex{ // to avoid connections between vertices of the same triangle
  float value;
  int mesh;
  int triangle;
};
struct infoVertexSeed{ // to avoid connections between vertices of the same triangle ... for seeds
  int loc;
  int mesh;
  vector<int> triangles;
};

  class MatCell_cmpr{
    // This class contains compressed information on entries for matrix4 format
  public:
    MatCell_cmpr():code2(0) {}
    MatCell_cmpr(double val):code2(0) { }
    MatCell_cmpr(const MatCell_cmpr& rhs){ *this=rhs; }
    void add_one(float dist,int fib); 
    void add_n(float dist,vector<float> props,int n);
    void add_n(int64_t newcode2);
    int64_t getcode2() const { return code2; }
    MatCell_cmpr& operator=(const MatCell_cmpr& rhs){
      code2=rhs.code2;
      return *this;
    }
    MatCell_cmpr& operator+=(const MatCell_cmpr& rhs){      
      add_n(rhs.code2);
      return *this;
    }
    void Print(){
      int nsamples, fibcnt1, fibcnt2; float length_tot;
      decode(code2,nsamples, fibcnt1, fibcnt2, length_tot);
      cout<<"code2 = "<<code2<<endl;
      cout<<"nsamples = "<<nsamples<<endl;
      cout<<"fibcnt1  = "<<fibcnt1<<endl;
      cout<<"fibcnt2  = "<<fibcnt2<<endl;
      cout<<"fibcnt3  = "<<nsamples-fibcnt1-fibcnt2<<endl;
      cout<<"length_tot = "<<length_tot<<endl;
      cout<<"avg_length = "<<(nsamples!=0?length_tot/(float)nsamples:0.0)<<endl;
      cout<<"----------------------"<<endl;
    }
 
  private:
    int64_t code2; //Compressed info on fibre_count, fibre1_prop, fibre1_prop, avg_length: code2 = two32*fibre_count + mult*mult*fibre_prop1 + mult*fibre_prop2 + length_val
    void decode(int64_t incode, int& nsamples, int& fibcnt1, int& fibcnt2, float& length_tot) const;
    void encode(const int nsamples, const int fibcnt1, const int fibcnt2, const float length_tot);
  };
  inline MatCell_cmpr operator*(const double& x, const MatCell_cmpr& rhs){ return rhs;}
  


  class SpMatHCPException: public std::exception{
  private:
    std::string m_msg;
  public:
    SpMatHCPException(const std::string& msg) throw(): m_msg(msg){ }
    virtual const char * what() const throw() { return string("SpMat_HCP::" + m_msg).c_str(); }
    ~SpMatHCPException() throw() {}
  };



  class SpMat_HCP : public SpMat<MatCell_cmpr>{
    public:
    SpMat_HCP():SpMat<MatCell_cmpr>::SpMat(){}
    SpMat_HCP(unsigned int m, unsigned int n):SpMat<MatCell_cmpr>::SpMat(m,n){}
    SpMat_HCP(unsigned int m, unsigned int n,const string& basename);
    ~SpMat_HCP(){}
    int SaveTrajFile(const string& basename)const;
    void AddToTraj(unsigned int r,unsigned int c, float dist,int fib){ here(r,c).add_one(dist,fib); } 
    void AddToTraj(unsigned int r,unsigned int c, float dist, vector<float> props, int n){ here(r,c).add_n(dist,props,n); }
    void AddToTraj(unsigned int r,unsigned int c, int64_t Newcode2){ here(r,c).add_n(Newcode2); }
    
    void Print(){
      for(unsigned int c=0;c<Ncols();c++){
	const std::vector<unsigned int>&    ri = get_ri(c);
	for (unsigned int r=0; r<ri.size(); r++) { 	  
	  cout<<"Element " << ri[r]+1 <<","<<c+1<<endl;
	  Peek(ri[r]+1,c+1).Print();
	}
      }
    }
    void Print(int r,int c){
      Peek(r,c).Print();
    }
  };


/*
  //Old SpMat_HCP with MattCell (less rounding errors, but takes ~5 times as much memory and ~40% more execution time

  class MatCell{
    // This class contains information on entries for matrix4 format
  public:
    MatCell():nsamples(0),length_tot(0.0){fibcnt.clear();fibcnt.resize(3,0);}
    MatCell(double val):nsamples(0),length_tot(0.0){fibcnt.clear();fibcnt.resize(3,0);}
    float   get_avg_length()const{
      return (nsamples!=0?length_tot/(float)nsamples:0.0);
    }
    float   get_length_tot()const{return length_tot;}
    int     get_nsamples()const{return nsamples;}
    float   get_fibprop(const int& f)const{return float(fibcnt[f-1])/float(nsamples);}
    float   get_fibcnt(const int& f)const{return fibcnt[f-1];}
    void    add_one(float dist,int fib){
      fibcnt[fib-1]+=1;
      length_tot+=dist;
      nsamples+=1;      
    }
    void    add_n(float dist,vector<float> props,int n){
      if(fibcnt.size()!=3){
	cerr<<"MatCell::add_n:Only valid with 3 fibres"<<endl;
	exit(1);
      }
      int n0=(int)round(props[0]*n);
      int n1=(int)round(props[1]*n);
      fibcnt[0]+=n0;
      fibcnt[1]+=n1;
      fibcnt[2]+=(n-n0-n1);
      length_tot+=(dist*n);
      nsamples+=n;      
    }
    void Print(){
      cout<<"nsamples   = "<<nsamples<<endl;
      cout<<"fibcnt[0]  = "<<fibcnt[0]<<endl;
      cout<<"fibcnt[1]  = "<<fibcnt[1]<<endl;
      cout<<"fibcnt[2]  = "<<fibcnt[2]<<endl;
      cout<<"length_tot = "<<length_tot<<endl;
      cout<<"avg_length = "<<length_tot/float(nsamples)<<endl;
      cout<<"----------------------"<<endl;
    }
    MatCell(const MatCell& rhs){
      *this=rhs;
    }
    MatCell& operator=(const MatCell& rhs){
      fibcnt=rhs.fibcnt;
      nsamples=rhs.nsamples;
      length_tot=rhs.length_tot;
      return *this;
    }
    MatCell& operator+=(const MatCell& rhs){      
      length_tot += rhs.get_length_tot();
      nsamples   += rhs.get_nsamples();
      for(unsigned int i=0;i<fibcnt.size();i++)
	fibcnt[i] += rhs.get_fibcnt(i+1);     
      return *this;
    }
  private:
    vector<int> fibcnt;
    int         nsamples;
    float       length_tot;
  };
  inline MatCell operator*(const double& x, const MatCell& rhs){return rhs;}

  class SpMat_HCP : public SpMat<MatCell>{
    public:
    SpMat_HCP():SpMat<MatCell>::SpMat(){}
    SpMat_HCP(unsigned int m, unsigned int n):SpMat<MatCell>::SpMat(m,n){}
    SpMat_HCP(unsigned int m, unsigned int n,const string& basename);
    ~SpMat_HCP(){}
    // HCP Trajectory-file writer (MJ+SJ)
    int SaveTrajFile(const string& basename)const;
    void AddToTraj(unsigned int r,unsigned int c, float dist,int fib){ here(r,c).add_one(dist,fib); } 
    void AddToTraj(unsigned int r,unsigned int c, float dist, vector<float> props, int n){ here(r,c).add_n(dist,props,n); }
    
    void Print(){
      for(unsigned int c=0;c<Ncols();c++){
	const std::vector<unsigned int>&    ri = get_ri(c);
	for (unsigned int r=0; r<ri.size(); r++) { 	  
	  cout<<"Element " << ri[r]+1 <<","<<c+1<<endl;
	  here(ri[r]+1,c+1).Print();
	}
      }
    }
    void Print(int r,int c){
      here(r,c).Print();
    }
  };

*/


namespace TRACT{

  void read_ascii_files(const string& filename,vector<string>& content);  

  //  the following is a helper function for save_matrix*
  //  to convert between nifti coords (external) and newimage coord (internal)
  void applycoordchange(volume<int>& coordvol, const Matrix& old2new_mat);
  void applycoordchange(Matrix& coordvol, const Matrix& old2new_mat);

  class Streamliner{
    //Everything in DTI space is done INSIDE this class and lower level classes (particle and tractvolsx)
    //This class communicates with higher level classes in Seed voxels.
    //
    probtrackxOptions&            opts;
    Log&                          logger;

    Particle                      m_part;
    vector<ColumnVector>          m_path;
    vector<ColumnVector>          m_diff_path;
    int                           m_tracksign;

    volume<float>                 m_mask;


    // prior masks
    volume<int>                   m_skipmask;
    CSV                           m_rubbish;
    CSV                           m_stop; 
    CSV                           m_waymasks;
    CSV                           m_netmasks;
    CSV                           m_wtstopmasks;

    vector<int>                   m_way_passed_flags;
    
    // for network mode
    int                           m_seed_id;
    Matrix                        m_network_mat;
    ColumnVector                  m_net_passed_flags;
    ColumnVector                  m_net_passed;

    vector< vector<ColumnVector> > m_crossedvox;
    bool                           m_surfexists;

    string                        m_waycond;
    
    volume4D<float>               m_prefdir;
    volume4D<float>               m_loopcheck;
    volume<float>                 m_loccurvthresh;

    // transform seed<->diff space
    Matrix                        m_Seeds_to_DTI;
    Matrix                        m_DTI_to_Seeds;
    volume4D<float>               m_Seeds_to_DTI_warp;
    volume4D<float>               m_DTI_to_Seeds_warp;
    bool                          m_IsNonlinXfm;
    // rotdir stuff
    Matrix                        m_rotdir;
    volume4D<float>               m_jacx,m_jacy,m_jacz;

    Tractvolsx                    vols;
    float                         m_lcrat;
    float                         m_x_s_init;
    float                         m_y_s_init;
    float                         m_z_s_init;

    // Streamliner needs to know about matrix3 
    CSV                           m_mask3;
    CSV                           m_lrmask3;
    vector< pair<int,infoVertex> > m_inmask3; // knows which node in mask3 and how far from seed (signed distance)
    vector< pair<int,infoVertex> > m_inlrmask3;
    vector< pair<int,infoVertex> > m_inmask3_aux; // write here and update m_inmask3 only if a part is accepted
    vector< pair<int,infoVertex> > m_inlrmask3_aux;	

    // we need this class to know about seed space
    const CSV&                    m_seeds;

  public:
    //Constructors
    Streamliner  (const CSV& seeds);
    ~Streamliner (){}
    const CSV&         get_seeds()  const {return m_seeds;}
    const Tractvolsx&  get_vols()   const {return vols;}
    inline int         nfibres()    const {return vols.get_nfibres();}
    inline const float get_x_seed() const {return m_x_s_init;}
    inline const float get_y_seed() const {return m_y_s_init;}
    inline const float get_z_seed() const {return m_z_s_init;}
    const vector<ColumnVector>& get_path_ref() const{return m_path;}
    vector<ColumnVector>        get_path()     const{return m_path;}
    const vector<ColumnVector>& get_diff_path_ref() const{return m_diff_path;}
    vector<ColumnVector>        get_diff_path()     const{return m_diff_path;}

    vector< vector<ColumnVector> >        get_crossedvox()const{return m_crossedvox;}
    const vector< vector<ColumnVector> >& get_crossedvox_ref()const{return m_crossedvox;}

    void surfexists(){m_surfexists=true;}

    inline void reset(){
      m_part.reset();
      vols.reset(opts.fibst.value());
      for(unsigned int i=0;i<m_way_passed_flags.size();i++)
	m_way_passed_flags[i]=0;
      if(opts.network.value()){
	m_net_passed_flags=0;
	m_net_passed=0;
      }
      m_tracksign=1;
    }
    inline void reverse(){
      m_part.restart_reverse();
      m_tracksign=-1;
    }
    void rotdir(const ColumnVector& dir,ColumnVector& rotdir,
		const float& x,const float& y,const float& z);

    int streamline(const float& x_init,const float& y_init, const float& z_init,
		   const ColumnVector& dim_seeds,const int& fibst);

    
    // separate masks loading from class constructor
    void load_netmasks(const string& filename,const int& excl){
      vector<int> vexcl;vexcl.push_back(excl);
      load_netmasks(filename,vexcl);
    }
    void load_netmasks(const string& filename,const vector<int>& exclude){

      string tmpfilename=logger.appendDir("tmpnetmaskfile");
      ofstream of(tmpfilename.c_str());
      vector<string> filelist;

      read_ascii_files(filename,filelist);

      int nfiles=0;
      for(unsigned int i=0;i<filelist.size();i++){
	bool iselmnt=false;
	for(unsigned int j=0;j<exclude.size();j++){
	  if(i==(unsigned int)exclude[j]){iselmnt=true;break;}
	}
	if(!iselmnt){
	  nfiles++;
	  of<<filelist[i]<<endl;
	}
      }
      
      m_netmasks.reinitialize(m_seeds.get_refvol());
      m_netmasks.set_convention(opts.meshspace.value());
      m_netmasks.load_rois(tmpfilename);    
      if(m_netmasks.nSurfs()>0){surfexists();}
      m_net_passed_flags.ReSize(m_netmasks.nRois());
      m_net_passed_flags=0;
      m_net_passed.ReSize(m_netmasks.nRois());
      m_net_passed=0;

    }
    void set_seed_id(const int i){m_seed_id=i;}
    void init_network_mat(const int n){
      m_network_mat.ReSize(n,n);
      m_network_mat=0;
    }
    void update_mat(){
      for(int i=1;i<=m_net_passed.Nrows();i++){
	if(m_net_passed(i)==0){continue;}
	if(m_seed_id+1>i){m_network_mat(m_seed_id+1,i)++;}
	else{m_network_mat(m_seed_id+1,i+1)++;}    
      }   
    }
    void save_network_mat(){
      write_ascii_matrix(m_network_mat,logger.appendDir("fdt_network_matrix"));
    }

    void load_waymasks(const string& filename){
      m_waymasks.reinitialize(m_seeds.get_refvol());
      m_waymasks.set_convention(opts.meshspace.value());
      m_waymasks.load_rois(filename);    
      if(m_waymasks.nSurfs()>0){surfexists();}
      m_way_passed_flags.clear(); 
      for(int i=0;i<m_waymasks.nRois();i++)
	m_way_passed_flags.push_back(0);
    }

    void load_wtstopmasks(const string& filename){
      m_wtstopmasks.reinitialize(m_seeds.get_refvol());
      m_wtstopmasks.set_convention(opts.meshspace.value());
      m_wtstopmasks.load_rois(filename);    
      if(m_wtstopmasks.nSurfs()>0){surfexists();} //cerr<<"Surface has been provided as a walk-through stopping mask! Currently unsupported!"; exit(1); }
    }

    void    set_waycond(const string& cond){m_waycond=cond;}
    string  get_wacond()const{return m_waycond;}

    void load_stop(const string& filename){    
      m_stop.reinitialize(m_seeds.get_refvol());
      m_stop.set_convention(opts.meshspace.value());
      m_stop.load_rois(filename);     
      if(m_stop.nSurfs()>0){surfexists();}
    }
    void load_rubbish(const string& filename){      
      m_rubbish.reinitialize(m_seeds.get_refvol());
      m_rubbish.set_convention(opts.meshspace.value());
      m_rubbish.load_rois(filename);     
      if(m_rubbish.nSurfs()>0){surfexists();}
    }

    // //////    matrix3 methods
    void init_mask3(){
      m_mask3.reinitialize(m_seeds.get_refvol());
      m_mask3.set_convention(opts.meshspace.value());
      m_mask3.load_rois(opts.mask3.value());
      if(m_mask3.nSurfs()>0){surfexists();}

      if(opts.lrmask3.value()!=""){
	m_lrmask3.reinitialize(m_seeds.get_refvol());
	m_lrmask3.set_convention(opts.meshspace.value());      
	m_lrmask3.load_rois(opts.lrmask3.value());
	if(m_lrmask3.nSurfs()>0){surfexists();}
      }
    }
    void                       clear_inmask3()   {m_inmask3.clear();m_inmask3_aux.clear();}
    void                       clear_inlrmask3() {m_inlrmask3.clear();m_inlrmask3_aux.clear();}
    void		       reset_m_inmask3_aux() {m_inmask3_aux.clear();m_inlrmask3_aux.clear();}
    vector< pair<int,infoVertex> >& get_inmask3()     {return m_inmask3;}
    vector< pair<int,infoVertex> >& get_inlrmask3()   {return m_inlrmask3;}
    CSV                        get_mask3()       {return m_mask3;}
    CSV                        get_lrmask3()     {return m_lrmask3;}
    void fill_inmask3(const vector<int>& crossedlocs3,vector< pair<int,int> >& surf_Triangle,const float& pathlength){
      vector< pair<int,infoVertex> > inmask3;
      for(unsigned int iter=0;iter<crossedlocs3.size();iter++){
	pair<int,infoVertex> mypair;
	mypair.first=crossedlocs3[iter];
	mypair.second.value=m_tracksign*pathlength;
	mypair.second.mesh=surf_Triangle[iter].first;
	mypair.second.triangle=surf_Triangle[iter].second;
	inmask3.push_back(mypair);
      }
      m_inmask3_aux.insert(m_inmask3_aux.end(),inmask3.begin(),inmask3.end());
    }
    void fill_inlrmask3(const vector<int>& crossedlocs3,vector< pair<int,int> >& surf_Triangle,const float& pathlength){
      vector< pair<int,infoVertex> > inmask3;
      for(unsigned int iter=0;iter<crossedlocs3.size();iter++){
	pair<int,infoVertex> mypair;
	mypair.first=crossedlocs3[iter];
	mypair.second.value=m_tracksign*pathlength;
	mypair.second.mesh=surf_Triangle[iter].first;
	mypair.second.triangle=surf_Triangle[iter].second;
	inmask3.push_back(mypair);
      }
      m_inlrmask3_aux.insert(m_inlrmask3_aux.end(),inmask3.begin(),inmask3.end());
    }
    void copy_inmask3(){
	m_inmask3.insert(m_inmask3.end(),m_inmask3_aux.begin(),m_inmask3_aux.end());
	m_inlrmask3.insert(m_inlrmask3.end(),m_inlrmask3_aux.begin(),m_inlrmask3_aux.end());
    }	
    // /////////////////////////////////////////////////////////////////
  };


  class Counter{
    probtrackxOptions&           opts;
    Log&                         logger;

    volume<float>                m_prob;      // spatial histogram of tract location within brain mask (in seed space)
    volume<float>                m_prob2;     // for mean path length	
    volume4D<float>              m_localdir;
    volume<int>                  m_beenhere;
    Matrix                       m_I;
    vector<ColumnVector>         m_path;
    vector<ColumnVector>         m_diff_path;
    vector< vector<ColumnVector> > m_crossedvox;
    CSV                          m_prob_alt;  // spatial histogram of tracts with alternative user-defined mask
    CSV                          m_prob_alt2; // for mean path length
    CSV                          m_beenhere_alt;

    // same as m_prob and m_localdir but split into the different
    // target masks if the option opts.targetpaths is ON
    vector< volume<float> >      m_prob_multi;
    vector< volume<float> >      m_prob_multi2; // for mean path length
    vector< volume4D<float> >    m_localdir_multi;
    

    // temp 
    volume<float>                m_lastpoint; // store last point in trajectory

    vector< vector<ColumnVector> > m_save_paths;

    // do we still need these?
    vector<ColumnVector>         m_seedcounts;
    Matrix                       m_SeedCountMat;
    int                          m_numseeds;

    // know where we are in seed space/counts (because seeds are now CSV)
    string                       m_curtype;
    infoVertexSeed	         m_curloc;

    // classification targets
    CSV                          m_targetmasks;
    vector<bool>                 m_targflags;
    CSV                          m_s2t_count;
    CSV                          m_s2t_count2;  // for mean path length
    Matrix                       m_s2tastext;
    Matrix                       m_s2tastext2;  // for mean path length
    int                          m_s2trow;
    volume4D<float>              m_targetpaths;

    // MATRIX 1
    SpMat<float>                *m_ConMat1; // using sparse
    SpMat<float>                *m_ConMat1b; // for mean path length
    int                          m_Conrow1;

    // MATRIX 2
    SpMat<float>                 *m_ConMat2; // using sparse
    SpMat<float>                 *m_ConMat2b; // for mean path length
    volume<int>                  m_lrmask;
    volume<int>                  m_lookup2;
    volume<int>                  m_beenhere2;
    ColumnVector                 m_lrdim;

    // MATRIX 3
    SpMat<float>                 *m_ConMat3; // using sparse
    SpMat<float>                 *m_ConMat3b; // for mean path length

    // MATRIX 4 - columns are seed space, rows are diffusion space
    SpMat_HCP                   *m_ConMat4;     
    volume<int>                  m_dtimask;
    CSV                          m_mask4;
    volume<int>                  m_lookup4;
    volume<int>                  m_beenhere4;
    ColumnVector                 m_dtidim;

    // misc    
    ColumnVector                 m_seedsdim;
    Streamliner&                 m_stline;
    
  public:
    Counter(Streamliner& stline):opts(probtrackxOptions::getInstance()),
				 logger(LogSingleton::getInstance()),						  
				 m_stline(stline)
    {
      m_numseeds = m_stline.get_seeds().nLocs();

      m_beenhere.reinitialize(m_stline.get_seeds().xsize(),
			      m_stline.get_seeds().ysize(),
			      m_stline.get_seeds().zsize());
      m_beenhere=0;
      m_seedsdim.ReSize(3);
      m_seedsdim << m_stline.get_seeds().xdim() 
		 << m_stline.get_seeds().ydim() 
		 << m_stline.get_seeds().zdim();
      m_I=IdentityMatrix(4);      

    }
    ~Counter(){}

    Streamliner& get_stline(){return m_stline;}

    void initialise();
    
    void initialise_path_dist(){
      if(opts.verbose.value()>0)
	cout<<"Initialise pathdist"<<endl;
      m_prob.reinitialize(m_stline.get_seeds().xsize(),
			  m_stline.get_seeds().ysize(),
			  m_stline.get_seeds().zsize());
      copybasicproperties(m_stline.get_seeds().get_refvol(),m_prob);
      m_prob=0;   
      if(opts.omeanpathlength.value()){    
        m_prob2.reinitialize(m_stline.get_seeds().xsize(),
			     m_stline.get_seeds().ysize(),
			     m_stline.get_seeds().zsize());
        copybasicproperties(m_stline.get_seeds().get_refvol(),m_prob2);
        m_prob2=0;
      }      
      if(opts.opathdir.value()){
	m_localdir.reinitialize(m_stline.get_seeds().xsize(),
				m_stline.get_seeds().ysize(),
				m_stline.get_seeds().zsize(),6);
	copybasicproperties(m_stline.get_seeds().get_refvol(),m_localdir);
	m_localdir=0;
      }
      

      if(opts.pathfile.set()){
	m_prob_alt.reinitialize(m_stline.get_seeds().get_refvol());
	m_prob_alt.set_convention(opts.meshspace.value());
	m_prob_alt.load_rois(opts.pathfile.value());
	m_prob_alt.reset_values();
        if(opts.omeanpathlength.value()){
          m_prob_alt2.reinitialize(m_stline.get_seeds().get_refvol());
	  m_prob_alt2.set_convention(opts.meshspace.value());
	  m_prob_alt2.load_rois(opts.pathfile.value());
	  m_prob_alt2.reset_values();
	}
	m_beenhere_alt.reinitialize(m_stline.get_seeds().get_refvol());
	m_beenhere_alt.set_convention(opts.meshspace.value());
	m_beenhere_alt.load_rois(opts.pathfile.value());
	m_beenhere_alt.set_vol_values(1);
	if(m_prob_alt.nSurfs()>0){m_stline.surfexists();}
      }
      if(opts.verbose.value()>0)
	cout<<"....done"<<endl;
    }
    void initialise_seedcounts();
    
    void initialise_matrix1(); 
    void initialise_matrix2();
    void initialise_matrix3();
    void initialise_matrix4();
    
    void forceNumSeeds(const int& n) {m_numseeds=n;}
    void updateSeedLocation(int loc, int roi, vector<int>& triangles) {
      m_curloc.triangles.clear();
      m_curloc.loc=loc;
      m_curloc.mesh=roi;
      for(unsigned int i=0;i<triangles.size();i++){
        m_curloc.triangles.push_back(triangles[i]);
      }
    }

    void store_path(){ 
      m_path=m_stline.get_path();
      if(opts.matrix4out.value())
	m_diff_path=m_stline.get_diff_path();
      m_crossedvox=m_stline.get_crossedvox();
    }
    void append_path(){      
      for(unsigned int i=0;i<m_stline.get_path_ref().size();i++){
	m_path.push_back(m_stline.get_path_ref()[i]);
	if(opts.matrix4out.value())
	  m_diff_path.push_back(m_stline.get_diff_path_ref()[i]);
      }
      for(unsigned int i=0;i<m_stline.get_crossedvox_ref().size();i++){
	m_crossedvox.push_back(m_stline.get_crossedvox_ref()[i]);
      }
    }
    float calc_pathlength(const int& redund=0){
      return( float(m_path.size()-redund)*opts.steplength.value() );
    }

    void clear_path(){ 
      m_path.clear(); 
      if(opts.matrix4out.value()){
	m_diff_path.clear(); 
      }      
      m_crossedvox.clear();
      
    };

    void count_streamline();
    void count_seed();
    void clear_streamline();
    
    
    void update_pathdist();
    void reset_beenhere();
    void update_pathdist_multi();
    
    void reset_prob(){m_prob=0;m_prob2=0;}
    void update_seedcounts();
    void reset_targetflags(){
      for(unsigned int i=0;i<m_targflags.size();i++) m_targflags[i]=false;
    }
    
    
    void update_matrix1(); //update path_dist after each streamline, only run this after each voxel!!
    
    void update_matrix2_row(); //but run this one every streamline as with the others
    void reset_beenhere2();

    void update_matrix3();
    void reset_beenhere3();

    void update_matrix4_col(); 
    void reset_beenhere4();
    
    void save_total(const int& keeptotal);
    void save_total(const vector<int>& keeptotal);
    void save();
    void save_pathdist();
    void save_pathdist(string add);
    void save_seedcounts();
    void save_matrix1();
    void save_matrix2();
    void save_matrix3();
    void save_matrix4();

    void add_path();
    void save_paths();
    
  };
  
  class Seedmanager{
    probtrackxOptions&        opts;
    Log&                      logger;

    Counter&                  m_counter;    
    ColumnVector              m_seeddims;

  public:
    Seedmanager(Counter& counter):opts(probtrackxOptions::getInstance()),
				  logger(LogSingleton::getInstance()),
				  m_counter(counter){
      m_seeddims.ReSize(3);
      m_seeddims<<m_counter.get_stline().get_seeds().xdim()
		<<m_counter.get_stline().get_seeds().ydim()
		<<m_counter.get_stline().get_seeds().zdim();
    }
    ~Seedmanager(){}
    Streamliner& get_stline(){return m_counter.get_stline();}

    int run(const float& x,const float& y,const float& z,
	    bool onewayonly, int fibst,float sampvox);



  };

}
