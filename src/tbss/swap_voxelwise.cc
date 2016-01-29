/*  tbss_reorder_vectors.cc

    Saad Jbabdi, FMRIB Image Analysis Group

    Copyright (C) 2009 University of Oxford */

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

#include "utils/options.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscprob.h"
#include "stdlib.h"

#define AP 2
#define NB 1
#define FA 0


using namespace Utilities;
using namespace std;
using namespace NEWIMAGE;
using namespace MISCMATHS;


string title="swap_voxelwise \nReordering of vectors with direction preservation";
string examples="swap_voxelwise -v <vectorsFileList> [-s <scalarsFileList> -b <outputBaseName> -m <mask>]";

Option<bool> help(string("-h,--help"),false,
		       string("display this message"),
		       false,no_argument);
Option<bool> verbose(string("-V,--verbose"),false,
		       string("switch on diagnostic messages"),
		       false,no_argument);
Option<string> vectors(string("-v,--vectors"),string(""),
		       string("vectors list (textfile)"),
		       true,requires_argument);
Option<string> scalars(string("-s,--scalars"),string(""),
		       string("scalars list (textfile)"),
		       true,requires_argument);
Option<string> mode(string("--mode"),string("voxels"),
		       string("reordering mode - choose between 'voxels' (default) or 'volumes' (in which the volumes must be aligned)"),
		       false,requires_argument);
Option<string> obasename(string("-b,--obasename"),string("reordered"),
		       string("output obasename - default='reordered'"),
		       false,requires_argument);
Option<string> bmask(string("-m,--mask"),string(""),
		       string("filename of brain mask or skeleton"),
		       true,requires_argument);
Option<string> initmask(string("--initmask"),string(""),
		       string("filename of initialisation mask"),
		       false,requires_argument);
Option<float> xthresh(string("--xthresh"),0.1,
		       string("threshold for considering a crosing fibre region - default=0.1"),
		       false,requires_argument);

////////////////////////////////////////////////////////

void read_masks(vector<string>& masks, const string& filename){
  ifstream fs(filename.c_str());
  string tmp;
  if(fs){
    fs>>tmp;
    do{
      masks.push_back(tmp);
      fs>>tmp;
    }while(!fs.eof());
  }
  else{
    cerr<<filename<<" does not exist"<<endl;
    exit(1);
  }
} 

// HEAP STRUCTURE FOR FAST MARCHING
/*---------------------*/
/* HEAP DATA STRUCTURE */
/*---------------------*/
class Heap{
 private:
  double *m_val;
  volume<int> m_bpr;
  int *m_x;
  int *m_y;
  int *m_z;
  int m_N;
  int m_size;

 public:

  Heap(const int& nx,const int& ny,const int& nz){

  int siz=MAX(MAX(nx,ny),MAX(ny,nz));
  siz*=siz*10;
  //  OUT(siz);

  m_size=siz;

  m_val=(double *)malloc(siz*sizeof(double));
  m_x=(int *)malloc(siz*sizeof(int));
  m_y=(int *)malloc(siz*sizeof(int));
  m_z=(int *)malloc(siz*sizeof(int));
  
  m_bpr.reinitialize(nx,ny,nz) ;
  
  m_N=0;

  m_val[0]=-1;

  //  cout<<"end of Heap constructor"<<endl;
  }

  inline int get_N()const{return m_N;}
  inline int get_bpr(int i,int j,int k)const{return m_bpr(i,j,k);}
  inline float get_val(int pos)const{return m_val[pos];}
  void set_val(int pos,float v){
    m_val[pos]=v;
  }

  void heapFree(){
    free(m_val);
    free(m_x);
    free(m_y);
    free(m_z); 
  }

void heapUp(int k)
{
  double v;
  int x,y,z;
  v = m_val[k];
  x = m_x[k];
  y = m_y[k];
  z = m_z[k];
  /*h->bpr[z][y][x]=k;*/
  while (m_val[k/2] > v & (k/2)>=1){
    m_val[k]=m_val[k/2];
    m_x[k]=m_x[k/2];
    m_y[k]=m_y[k/2];
    m_z[k]=m_z[k/2];
    m_bpr(m_x[k],m_y[k],m_z[k])=k;
    
    m_val[k/2]=v;
    m_x[k/2]=x;
    m_y[k/2]=y;
    m_z[k/2]=z;
    m_bpr(x,y,z)=k/2;
    
    k = k/2;
  }
  m_val[k]=v;
  m_x[k]=x;
  m_y[k]=y;
  m_z[k]=z;
  m_bpr(x,y,z)=k;
  
}

void heapDown(int k)
{
  int j;
  double v;
  int x,y,z;
  v = m_val[k];
  x = m_x[k];
  y = m_y[k];
  z = m_z[k];
  m_bpr(x,y,z)=k;
  while (k <= m_N/2)
    {
      j=k+k;
      if ((j<m_N) && (m_val[j]>m_val[j+1])) j++;
      if (v <= m_val[j]) break;
      m_val[k]=m_val[j];
      m_x[k]=m_x[j];
      m_y[k]=m_y[j];
      m_z[k]=m_z[j];
      m_bpr(m_x[k],m_y[k],m_z[k])=k;
      k = j;
    }
  m_val[k] = v;
  m_x[k] = x;
  m_y[k] = y;
  m_z[k] = z;
  m_bpr(x,y,z)=k;
  
}

void heapInsert(float v,int x,int y,int z)
{  
  m_val[++(m_N)] = v;
  m_x[(m_N)] = x;
  m_y[(m_N)] = y;
  m_z[(m_N)] = z;
  m_bpr(x,y,z)=m_N;
 
  heapUp(m_N);
}

void heapRemove(float& v, int& x, int& y, int& z)
{
  v = m_val[1];
  x = m_x[1];
  y = m_y[1];
  z = m_z[1];
  
  m_x[1] = m_x[(m_N)];
  m_y[1] = m_y[(m_N)];
  m_z[1] = m_z[(m_N)];
  m_val[1] = m_val[(m_N)--];
  m_bpr(m_x[1],m_y[1],m_z[1])=1;
  
  heapDown(1);
  
}

void heapRemovePoint(int x, int y, int z)
{
  int i;

  i=m_bpr(x,y,z);

  m_x[i]=m_x[(m_N)];
  m_y[i]=m_y[(m_N)];
  m_z[i]=m_z[(m_N)];
  m_val[i]=m_val[(m_N)--];
  m_bpr(m_z[i],m_y[i],m_x[i])=i;
  
  heapDown(i);
  
}

void print()
{
  
  cout<<"heap structure size = "<<m_N<<endl<<endl;
  cout<<"values : "<<endl<<endl;
  
  for(int i=1;i<=m_N;i++){
    cout<<i<<"\t"<<m_val[i]<<"\t i="<<m_x[i]<<",j="<<m_y[i]<<",k="<<m_z[i]<<endl;
    cout<<"bpr of "<<i<<" is "<<m_bpr(m_x[i],m_y[i],m_z[i])<<endl;
  }
}
};

void make_lut(const volume<float>& mask,volume<int>& vol2mat,Matrix& mat2vol){
  vol2mat.reinitialize(mask.xsize(),mask.ysize(),mask.zsize());
  vol2mat=0;
  int nnz=0;
  for(int z=0;z<mask.zsize();z++){
    for(int y=0;y<mask.ysize();y++){
      for(int x=0;x<mask.xsize();x++){
	if(mask(x,y,z)==0)continue;
	nnz++;
	vol2mat(x,y,z)=nnz;
      }
    }
  }
  mat2vol.ReSize(nnz,3);
  nnz=0;
  for(int z=0;z<mask.zsize();z++){
    for(int y=0;y<mask.ysize();y++){
      for(int x=0;x<mask.xsize();x++){
	if(mask(x,y,z)==0)continue;
	nnz++;
	mat2vol.Row(nnz) << x << y << z;
      }
    }
  }
}

void fill_matrix(Matrix& mat,const volume4D<float>& vol,const Matrix& mat2vol){
  mat.ReSize(mat2vol.Nrows(),vol.tsize());
  for(int i=1;i<=mat.Nrows();i++){
    for(int j=1;j<=mat.Ncols();j++){
      mat(i,j) = vol((int)mat2vol(i,1),(int)mat2vol(i,2),(int)mat2vol(i,3),j-1);
    }
  }
}

void fill_vector(ColumnVector& mat,const volume<float>& vol,const Matrix& mat2vol){
  mat.ReSize(mat2vol.Nrows());
  for(int i=1;i<=mat.Nrows();i++){
    mat(i) = vol((int)mat2vol(i,1),(int)mat2vol(i,2),(int)mat2vol(i,3));
  }
}

void fill_volume4D(volume4D<float>& vol,const Matrix& mat,const Matrix& mat2vol){
  for(int i=1;i<=mat.Nrows();i++){
    for(int j=1;j<=mat.Ncols();j++){
      vol((int)mat2vol(i,1),(int)mat2vol(i,2),(int)mat2vol(i,3),j-1)=mat(i,j);
    }
  }
}
void fill_volume(volume<float>& vol,const ColumnVector& mat,const Matrix& mat2vol){
  for(int i=1;i<=mat.Nrows();i++){
    vol((int)mat2vol(i,1),(int)mat2vol(i,2),(int)mat2vol(i,3))=mat(i);
  }
}



class FM{
 private:
  volume<float> mask,tmap,imask;
  volume<int> state;
  vector< Matrix > dyads;
  vector< ColumnVector >  meanf;
  int nx,ny,nz;
  Matrix P;
  Matrix neighbours;
  
  volume<int> vol2mat;
  Matrix mat2vol;

 public:
  FM(){
    // read the mask
    read_volume(mask,bmask.value());
    make_lut(mask,vol2mat,mat2vol);


    if(verbose.value())
      cout<<"reading vectors"<<endl;

    
    vector<string> vectornames;
    vector<string> scalarnames;

    {
      read_masks(vectornames,vectors.value());
      volume4D<float> tmpvol;
      Matrix tmpmat;
      for(unsigned int i=0;i<vectornames.size();i++){
	read_volume4D(tmpvol,vectornames[i]);
	fill_matrix(tmpmat,tmpvol,mat2vol);
	dyads.push_back(tmpmat);
      }
      
    }
    {
      read_masks(scalarnames,scalars.value());
      if(scalarnames.size() != vectornames.size()){
	cerr << "vector and scalar inputs do not have the same number of items" << endl;
	exit(1);
      }
      ColumnVector tmpmat;
      volume<float> tmpvol;
      for(unsigned int i=0;i<scalarnames.size();i++){
	read_volume(tmpvol,scalarnames[i]);
	fill_vector(tmpmat,tmpvol,mat2vol);
	meanf.push_back(tmpmat);
      }
    }
  
  
    /////////////////////////////////

    if(verbose.value())
      cout<<"done reading"<<endl;

    // read input volumes
    nx=mask.xsize();
    ny=mask.ysize();
    nz=mask.zsize();

    state.reinitialize(nx,ny,nz);
    tmap.reinitialize(nx,ny,nz);

    state=FA;
    tmap=0;

    P=MISCMATHS::perms(dyads.size());
    //    OUT(P);

    neighbours.ReSize(3,26);
    neighbours << 1 << 0 << 0 << -1 << 0 << 0 << 1 << 1 <<-1 <<-1 << 1 <<-1 << 1 <<-1 << 0 << 0 << 0 << 0 << 1 <<-1 << 1 << 1 <<-1 <<-1 << 1 <<-1 
	       << 0 << 1 << 0 << 0  <<-1 << 0 << 1 <<-1 << 1 <<-1 << 0 << 0 << 0 << 0 << 1 <<-1 << 1 <<-1 << 1 << 1 <<-1 << 1 <<-1 << 1 <<-1 <<-1
	       << 0 << 0 << 1 << 0  << 0 <<-1 << 0 << 0 << 0 << 0 << 1 << 1 <<-1 <<-1 << 1 << 1 <<-1 <<-1 << 1 << 1 << 1 <<-1 << 1 <<-1 <<-1 <<-1;
  }
  
  ////////////////////// functions
  void do_fast_marching(){
    int i0=0,j0=0,k0=0;           /* current point coords */
    float tval; 

     /* create heap sort structure */
    Heap h(nx,ny,nz);

    
    /********************************/
    /*** performing fmt algorithm ***/
    /********************************/
    
    /*** initialization ***/
    if(verbose.value())
      cout<<"initialise"<<endl;

    if(initmask.value() == ""){
      /* look for bigger f1+f2 point as a seed */
      float maxf=0,curf;
      //OUT(nx);OUT(ny);OUT(nz);
      for(int z=0;z<nz;z++)
	for(int y=0;y<ny;y++)
	  for(int x=0;x<nx;x++){
	    if(mask(x,y,z)==0)
	      continue;
	    curf=0;
	    for(unsigned int f=0;f<meanf.size();f++)
	      curf += meanf[f](vol2mat(x,y,z));
	    if(curf>maxf){
	      i0=x;j0=y;k0=z;
	      maxf=curf;
	    }
	    state(x,y,z)=FA;
	    tmap(x,y,z)=1;
	  }
      state(i0,j0,k0)=AP;
      tmap(i0,j0,k0)=0;
      updateNBvalue(i0,j0,k0,h);
    }
    else{
      /* make all voxels in init mask as AP */
      read_volume(imask,initmask.value());
      for(int z=0;z<nz;z++)
	for(int y=0;y<ny;y++)
	  for(int x=0;x<nx;x++){
	    if(mask(x,y,z)==0)
	      continue;
	    
	    if(imask(x,y,z)==0){
	      state(x,y,z)=FA;
	      tmap(x,y,z)=1;
	    }
	    else{
	      //	      OUT(z);
	      state(x,y,z)=AP;
	      tmap(x,y,z)=0;
	      updateNBvalue(x,y,z,h);
	      i0=x;j0=y;k0=z;
	    }
	  }
    }


    //cout<<"nbvalue1"<<endl;
    /* and all points of the ROIs as Alive Points */
    
    //h.print();

    //return;
    /*--------------------------------------------------------*/
    /*** big loop ***/
    if(verbose.value())
      cout<<"start FM"<<endl;
    int STOP = 0;
    int counter = 0;
    while(STOP==0){
      /*break;*/
      counter++;
      //OUT(counter);
      
      /*** uses the heap sort structure to find the NB point with the smallest T-value ***/
      h.heapRemove(tval,i0,j0,k0);
      //cout << i0 << " " << j0 << " " << k0 << endl;
      
      /*** add this point to the set of alive points ***/
      state(i0,j0,k0)=AP;
      
      if(h.get_N()==0)
	break;
      
      /*** update narrow band's T-value ***/
      updateNBvalue(i0,j0,k0,h);
      
    }
    
  }
  bool isInside(int i,int j,int k){
    //cout<<"je suis dans isInside"<<endl;
    //OUT(i);OUT(j);OUT(k);
    //OUT(mask(i,j,k));
    return((i>=0) && (i<nx)
	   && (j>=0) && (j<ny) 
	   && (k>=0) && (k<nz) 
	   && (mask(i,j,k)!=0));
  }
  void computeT(int x,int y,int z){
    Matrix V(P.Ncols(),3);
    Matrix nV(P.Ncols(),3);
    ColumnVector mF(P.Ncols());

    V=get_vector(x,y,z);
    mF=get_f(x,y,z);

    //cout << "NB voxel: ";
    //cout << x << " " << y << " " << z << endl;

    int nbx,nby,nbz;
    int opt_perm=1;
    float opt_sperm=0;
    for(int per=1;per<=P.Nrows();per++){
      //OUT(per);
      float opt_nb=0;int nnb=0;
      for(int nb=1;nb<=neighbours.Ncols();nb++){
	//OUT(nb);
	nbx=x+(int)neighbours(1,nb);nby=y+(int)neighbours(2,nb);nbz=z+(int)neighbours(3,nb);
	if(!isInside(nbx,nby,nbz))
	  continue;
	if(state(nbx,nby,nbz)==AP){
	  //cout<<"this neighbour is an AP"<<endl;
	  nV=get_vector(nbx,nby,nbz);
	  
	  float opt_s=0,s;
	  for(int f=1;f<=V.Nrows();f++){
	    s=abs(dot(nV.Row(f).t(),V.Row((int)P(per,f)).t()));
	    if(s>opt_s){
	      opt_s=s;
	    }
	  }//f
	  opt_nb+=opt_s;
	  nnb++;
	}//endif
	
      }//nb

      opt_nb/=nnb;
      //OUT(opt_nb);
	
      if(opt_nb>opt_sperm){
	opt_sperm=opt_nb;
	opt_perm=per;
      }

    }//perm
    tmap(x,y,z)=1-opt_sperm; // store optimal mean scalar product

    int nfibres = 0;
    for(int i=1;i<=V.Nrows();i++){
      if(mF(i) >= xthresh.value()){
	nfibres++;
      }
    }
    if(nfibres<2) opt_perm=2;


    // reorder dyads
    if(initmask.value()!="")
      if(imask(x,y,z)!=0){
	OUT(state(x,y,z));
	cerr<<"WARNING!!!!! Changing Init Mask!!!!"<<endl;
      }

    for(unsigned int f=0;f<dyads.size();f++){
      dyads[f](vol2mat(x,y,z),1)=V((int)P(opt_perm,f+1),1);
      dyads[f](vol2mat(x,y,z),2)=V((int)P(opt_perm,f+1),2);
      dyads[f](vol2mat(x,y,z),3)=V((int)P(opt_perm,f+1),3);

      meanf[f](vol2mat(x,y,z))=mF((int)P(opt_perm,f+1));
    }
  }

  ReturnMatrix get_vector(int x,int y,int z){
    Matrix V(P.Ncols(),3);
    for(unsigned int f=0;f<dyads.size();f++)
      for(int i=1;i<=3;i++)
	V(f+1,i)=dyads[f](vol2mat(x,y,z),i);

    V.Release();
    return V;
  }

  ReturnMatrix get_f(int x,int y,int z){
    ColumnVector V(P.Ncols());
    for(unsigned int f=0;f<dyads.size();f++)
      V(f+1)=meanf[f](vol2mat(x,y,z));

    V.Release();
    return V;
  }

  void updateNBvalue(int i,int j,int k,Heap& h){    
    int ni,nj,nk;
    int pos;
    double val;
   
    for(int nb=1;nb<=neighbours.Ncols();nb++){
      ni=i+(int)neighbours(1,nb);nj=j+(int)neighbours(2,nb);nk=k+(int)neighbours(3,nb);
      if(isInside(ni,nj,nk)){
	if(state(ni,nj,nk)!=AP && imask(ni,nj,nk)==0){
	  computeT(ni,nj,nk);
	  val=tmap(ni,nj,nk);
	  /* update value if in NB */
	  if(state(ni,nj,nk)==NB){
	    pos=h.get_bpr(ni,nj,nk);
	    h.set_val(pos,MIN(h.get_val(pos),val));
	    h.heapUp(pos);
	  }
	  /* insert value if in FA */
	  if(state(ni,nj,nk)==FA){
	    state(ni,nj,nk)=NB;
	    h.heapInsert(val,ni,nj,nk);
	  }
	}
      }
      
    }
  }
  
  void save_results(){
    if(verbose.value())
      cout<<"saving results"<<endl;

    volume4D<float> tmpd;
    volume<float> tmpf;
    tmpd.reinitialize(mask.xsize(),mask.ysize(),mask.zsize(),3);
    tmpf.reinitialize(mask.xsize(),mask.ysize(),mask.zsize());

    copybasicproperties(mask,tmpd[0]);
    copybasicproperties(mask,tmpf);

    int fib=1;
    for(unsigned int f=0;f<dyads.size();f++){
      tmpd=0;tmpf=0;
      fill_volume4D(tmpd,dyads[f],mat2vol);
      fill_volume(tmpf,meanf[f],mat2vol);

      save_volume4D(tmpd,obasename.value()+"_vectors"+num2str(fib));
      save_volume(tmpf,obasename.value()+"_scalars"+num2str(fib));
      fib++;
    }
    //cout << "saved safely" << endl;
  }

};

void do_crossvolume_reorder(){
  
  cerr << "This mode is not implemented!" << endl;   

}



int main(int argc,char *argv[]){

  Tracer tr("main");
  OptionParser options(title,examples);

  try{
    options.add(help);
    options.add(verbose);
    options.add(vectors);
    options.add(scalars);
    options.add(mode);
    options.add(obasename);
    options.add(bmask);
    options.add(initmask);
    options.add(xthresh);

    options.parse_command_line(argc,argv);

    
    if ( (help.value()) || (!options.check_compulsory_arguments(true)) ){
      options.usage();
      exit(EXIT_FAILURE);
    }

    if( mode.value()=="voxels" ) {
      if(verbose.value())
	cout<<"call for fast marching"<<endl;
      FM fm;
      if(verbose.value())
	cout<<"perform fast marching"<<endl;
      fm.do_fast_marching();
      if(verbose.value())
	cout<<"save results"<<endl;
      fm.save_results();
    }
    else if( mode.value()=="volumes"){
      cout << "volumes mode not implemented yet!" << endl;
      exit(0);

      if(verbose.value())
	cout<<"reorder across volumes"<<endl;

      do_crossvolume_reorder();


    }
    else{
      cerr << "unkown mode: " << mode.value() << endl;
      exit(1);
    }
    



  }
  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } 
  catch(std::exception &e) {
    cerr << e.what() << endl;
  } 

  return 0;
  
  
}
