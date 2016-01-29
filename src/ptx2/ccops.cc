/*  ccops.cc

    Tim Behrens, Saad Jbabdi, FMRIB Image Analysis Group

    Copyright (C) 1999-2010 University of Oxford  */

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

#include <iostream>
#include <fstream>
#include <cmath>
#include "newimage/newimageall.h"
#include "ccopsOptions.h"
#include <vector>
#include <algorithm>
#include "csv.h"
#include "miscmaths/SpMat.h"


using namespace std;
using namespace NEWIMAGE;
using namespace NEWMAT;
using namespace MISCMATHS;
using namespace CCOPS;

///////  RE-ORDERING //////////////////////

void spect_reord(SymmetricMatrix& A,ColumnVector& r,ColumnVector& y){
  SymmetricMatrix Q=-A;
  DiagonalMatrix t(Q.Nrows());
  t=0;
  DiagonalMatrix D;
  Matrix V;
  for(int i=1;i<=Q.Nrows();i++){
    float rowsum1=0, rowsum2=0;
    for(int j=1;j<=Q.Ncols();j++){
      if(i!=j) rowsum1+=Q(i,j);
      rowsum2+=A(i,j);
    }
    Q(i,i)=-rowsum1;
    t(i)=1/sqrt(rowsum2);
  }
  Q << t*Q*t;
  EigenValues(Q,D,V);
  vector<pair<float,int> > myvec;
  vector<pair<float,int> > myvec2;
  
  for(int i=1;i<=D.Nrows();i++){
    pair<float,int> mypair;
    mypair.first=D(i);
    mypair.second=i;
    myvec.push_back(mypair);
  }
  sort(myvec.begin(),myvec.end());
  int ind=myvec[1].second; // index for second eigenval
  
  ColumnVector v2scale(V.Nrows());
  for(int i=1;i<=V.Nrows();i++){
    v2scale(i)=V(i,ind); //second eigvec
  }
  v2scale=t*v2scale; //scale it
  
  
  for(int i=1;i<=D.Nrows();i++){
    pair<float,int> mypair;
    mypair.first=v2scale(i);
    mypair.second=i;
    myvec2.push_back(mypair);
  }
  //myvec2 contains scaled second eigenvector and index for sorting.
  
  sort(myvec2.begin(),myvec2.end());
  r.ReSize(D.Nrows());
  y.ReSize(D.Nrows());
  
  for(int i=1;i<=D.Nrows();i++){
    y(i)=myvec2[i-1].first;
    r(i)=myvec2[i-1].second;
  }
  
} 

bool compare(const pair<float,int> &r1,const pair<float,int> &r2){
  return (r1.first<r2.first);
}
void randomise(vector< pair<float,int> >& r){
  for(unsigned int i=1;i<=r.size();i++){
    pair<float,int> p(rand()/float(RAND_MAX),i);
    r[i-1]=p;
  }
  sort(r.begin(),r.end(),compare);
  
}
void do_kmeans(const Matrix& data,ColumnVector& y,const int k){
  int numiter=50; // hard-coded number of iterations in kmeans
  if(data.Nrows() != (int)y.Nrows()){
    y.ReSize(data.Nrows());
  }
  int n = data.Nrows();
  int d = data.Ncols();
  Matrix means(d,k),newmeans(d,k);
  ColumnVector nmeans(k);
  means=0;
  nmeans=0;
  

// initialise with far-away trick
// start with a random class centre. then each new class centre is
  // as far as possible from the cog of the previous classes
  means.Column(1) = data.Row(round(rand()/float(RAND_MAX)*float(n-1))+1).t();
  ColumnVector cog(d);
  for(int cl=2;cl<=k;cl++){
    cog = sum(means.SubMatrix(1,d,1,cl-1),2);
    
    int maxi=1;float dist=0,maxdist=0;
    for(int i=1;i<=n;i++){
      float cdist=0,mindist=-1;int minc=1;
      for(int prevcl=cl-1;prevcl>=1;prevcl--){
	cdist = (means.Column(prevcl)-data.Row(i).t()).SumSquare();
	if(mindist==-1 || cdist<mindist){mindist=cdist;minc=prevcl;}
      }
      dist = mindist;
      if(dist>=maxdist){maxdist=dist;maxi=i;}
    }
    means.Column(cl)=data.Row(maxi).t();
  }


  
  // iterate
  for(int iter=0;iter<numiter;iter++){
    // loop over datapoints and attribute z for closest mean
    newmeans=0;
    nmeans=0;
    for(int i=1;i<=n;i++){
      float mindist=1E20,dist=0;
      int mm=1;
      for(int m=1;m<=k;m++){
	dist = (means.Column(m)-data.Row(i).t()).SumSquare();
	if( dist<mindist){
	  mindist=dist;
	  mm = m;
	}
      }
      y(i) = mm;
      newmeans.Column(mm) += data.Row(i).t();
      nmeans(mm) += 1;
    }
    
    // compute means
    for(int m=1;m<=k;m++){
      if(nmeans(m)==0){
	cout << "Only found " << k-1 << " clusters!!!" << endl;
	do_kmeans(data,y,k-1);
	return;
      }
      newmeans.Column(m) /= nmeans(m);
    }
    means = newmeans;
  }
  
}


void kmeans_reord(const Matrix& A,ColumnVector& r,ColumnVector& y,const int k){
  do_kmeans(A,y,k);
 
  vector< pair<float,int> > myvec2;
  for(int i=1;i<=A.Nrows();i++){
    pair<int,int> mypair;
    mypair.first=(int)y(i);
    mypair.second=i;
    myvec2.push_back(mypair);
  }  
  
  sort(myvec2.begin(),myvec2.end());
  r.ReSize(A.Nrows());
  y.ReSize(A.Nrows());
  
  for(int i=1;i<=A.Nrows();i++){
    y(i)=myvec2[i-1].first;
    r(i)=myvec2[i-1].second;
  } 
}

void do_fuzzy(const Matrix& data,Matrix& u,const int k){
  int numiter     = 50; // hard-coded #iterations
  float fuzziness = 2;  // hard-coded fuzziness factor

  int n = data.Nrows();
  int d = data.Ncols();
  Matrix means(d,k),newmeans(d,k);
  ColumnVector nmeans(k);
  means=0;
  nmeans=0;
  
  // initialise with far-away trick
  // start with a random class centre. then each new class centre is
  // as far as possible from the cog of the previous classes
  means.Column(1) = data.Row(round(rand()/float(RAND_MAX)*float(n-1))+1).t();
  ColumnVector cog(d);
  for(int cl=2;cl<=k;cl++){
    cog = sum(means.SubMatrix(1,d,1,cl-1),2);    
    int maxi=1;float dist=0,maxdist=0;
    for(int i=1;i<=n;i++){
      float cdist=0,mindist=-1;int minc=1;
      for(int prevcl=cl-1;prevcl>=1;prevcl--){
	cdist = (means.Column(prevcl)-data.Row(i).t()).SumSquare();
	if(mindist==-1 || cdist<mindist){mindist=cdist;minc=prevcl;}
      }
      dist = mindist;
      if(dist>=maxdist){maxdist=dist;maxi=i;}
    }
    means.Column(cl)=data.Row(maxi).t();
  }

  // iterate
  for(int iter=0;iter<numiter;iter++){
    // loop over datapoints and attribute z for closest mean
    newmeans=0.0;
    nmeans=0;
    for(int i=1;i<=n;i++){
      float mindist=-1,dist=0;
      int mm=1;
      for(int m=1;m<=k;m++){
	dist = (means.Column(m)-data.Row(i).t()).SumSquare();
	if( mindist==-1 || dist<mindist){
	  mindist=dist;
	  mm = m;
	}
      }
      //y(i) = mm;
      newmeans.Column(mm) += data.Row(i).t();
      nmeans(mm) += 1;
    }
    
    // compute means
    for(int m=1;m<=k;m++){
      if(nmeans(m)!=0)
	newmeans.Column(m) /= nmeans(m);
      means.Column(m) = newmeans.Column(m);
    }
  }

  // now use this to calculate u
  u.ReSize(n,k);
  u=0.0;
  for(int i=1;i<=n;i++){
    for(int j=1;j<=k;j++){
      float xi_cj = (data.Row(i) - means.Column(j).t()).SumSquare();
      if(xi_cj==0){u.Row(i)=0.01/float(k-1);u(i,j)=0.99;break;}
      float xi_cl;
      for(int l=1;l<=k;l++){
	xi_cl = (data.Row(i) - means.Column(l).t()).SumSquare();
	u(i,j) += std::exp(std::log(xi_cj/xi_cl)/(fuzziness-1));
      }
      u(i,j) = 1/u(i,j);
    }
  }
  

}


void fuzzy_reord(const Matrix& A,Matrix& u,ColumnVector& r,ColumnVector& y,const int k){
  do_fuzzy(A,u,k);

  float junk;
  int index;
  y.ReSize(A.Nrows());
  r.ReSize(A.Nrows());
  for(int i=1;i<=A.Nrows();i++){
    junk = u.Row(i).Maximum1(index);
    y(i) = index;
  }
 
  vector< pair<float,int> > myvec2;
  for(int i=1;i<=A.Nrows();i++){
    pair<int,int> mypair;
    mypair.first=(int)y(i);
    mypair.second=i;
    myvec2.push_back(mypair);
  }  
  
  sort(myvec2.begin(),myvec2.end());
  r.ReSize(A.Nrows());
  y.ReSize(A.Nrows());
  
  
  for(int i=1;i<=A.Nrows();i++){
    y(i)=myvec2[i-1].first;
    r(i)=myvec2[i-1].second;
  } 
  
}

///////////// PRE-PROCESS MATRIX /////////////

void rem_zrowcol(const Matrix& myOM3Col,vector<int>& excl_cols,
		 const Matrix& coordmat,const Matrix& tractcoordmat,
		 const bool coordbool,const bool tractcoordbool,
		 Matrix& newOM3Col,Matrix& newcoordmat, Matrix& newtractcoordmat)
{
  // first pass to determine zero rows/columns to keep
  int nrows = (int)myOM3Col(myOM3Col.Nrows(),1);
  int ncols = (int)myOM3Col(myOM3Col.Nrows(),2);

  vector<int> keep_cols(ncols,0),keep_rows(nrows,0);
  for(int i=1;i<myOM3Col.Nrows();i++){
    if(myOM3Col(i,3)!=0){
      keep_cols[ (int)myOM3Col(i,2)-1 ] =1;
      keep_rows[ (int)myOM3Col(i,1)-1 ] =1;
    }
  }
  for(unsigned int i=0;i<excl_cols.size();i++)
    keep_cols[i]=0;

  // new indices for these rows/columns
  int newncols=0,newnrows=0;
  vector<int> lu_r(nrows,0);
  for(unsigned int i=0;i<keep_rows.size();i++){
    if(keep_rows[i]==0)continue;
    newnrows++;
    lu_r[i]=newnrows;
  }
  vector<int> lu_c(ncols,0);
  for(unsigned int i=0;i<keep_cols.size();i++){
    if(keep_cols[i]==0)continue;
    newncols++;
    lu_c[i]=newncols;
  }
  ////////////////////////////
  vector<int> _r,_c,_old_r,_old_c;
  vector<double> _v;

  int r,c;
  double v;
  for(int i=1;i<myOM3Col.Nrows();i++){
    r=(int)myOM3Col(i,1);
    c=(int)myOM3Col(i,2);
    v=(double)myOM3Col(i,3);
    
    //column is in excl_cols?
    if(keep_cols[c-1]==0){continue;}
    if(keep_rows[r-1]==0){continue;}
    
    _r.push_back(lu_r[r-1]);
    _c.push_back(lu_c[c-1]);
    _v.push_back(v);

  }

  // now make new matrix with remaining values
  newOM3Col.ReSize(_r.size()+1,3);
  for(unsigned int i=0;i<_r.size();i++){
    newOM3Col.Row(i+1) << _r[i] << _c[i] << _v[i];
  }
  newOM3Col.Row(_r.size()+1) << newnrows << newncols << 0;

  // newcoords
  if(coordbool){
    newcoordmat.ReSize(newnrows,coordmat.Ncols());
    for(int i=1;i<=nrows;i++){
      if(keep_rows[i-1]==0)continue;
      newcoordmat.Row(lu_r[i-1]) << coordmat.Row(i);
    }
  }
  if(tractcoordbool){
    newtractcoordmat.ReSize(newncols,tractcoordmat.Ncols());
    for(int i=1;i<=ncols;i++){
      if(keep_cols[i-1]==0)continue;
      newtractcoordmat.Row(lu_c[i-1]) << tractcoordmat.Row(i);
    }
  }
  
}




void add_connexity(SymmetricMatrix& CtCt,const Matrix& coord,const float p=.5){
  // compute CtCt range
  float r=CtCt.Minimum();
  float R=CtCt.Maximum();

  // compute distance matrix
  SymmetricMatrix D(coord.Nrows());
  for(int i=1;i<=coord.Nrows();i++)
    for(int j=1;j<=i;j++){
      

      D(i,j) = std::sqrt(
			 (coord(i,1)-coord(j,1))*(coord(i,1)-coord(j,1))+
			 (coord(i,2)-coord(j,2))*(coord(i,2)-coord(j,2))+
			 (coord(i,3)-coord(j,3))*(coord(i,3)-coord(j,3))
			 );
    }
  D=D.MaximumAbsoluteValue()-D;

  // change distance range
  float m=D.Minimum(),M=D.Maximum();
  D=(D-m)*(R-r)/(M-m)+r;

  // add distance to CtCt matrix
  for(int i=1;i<=coord.Nrows();i++)
    for(int j=1;j<=i;j++){
      CtCt(i,j)=std::sqrt((1-p)*CtCt(i,j)*CtCt(i,j)+p*D(i,j)*D(i,j));
    }
 
}



// corrcoef for a Nx3-represented sparse matrix
//  (like the matlab)
ReturnMatrix mycorrcoef(Matrix& M){
  int nrows = (int)M(M.Nrows(),1);
  int ncols = (int)M(M.Nrows(),2);

  // first pass to calculate means
  ColumnVector m(nrows),m2(nrows),s(nrows);
  m=0;m2=0;
  for(int i=1;i<M.Nrows();i++){
    m( (int)M(i,1) )  += M(i,3);
    m2( (int)M(i,1) ) += M(i,3)*M(i,3);
  }
  m /= float(ncols);
  s = MISCMATHS::sqrt(m2 - float(ncols)*NEWMAT::SP(m,m));
  
  // second pass to calculate corrcoef(M')
  SymmetricMatrix C(nrows);
  C=0;

  vector<int> is,rows;
  for(int i=1;i<M.Nrows();i++){
    if(i>=2){
      if(M(i,2)>M(i-1,2)){//we are still in the same column
	for(unsigned int ii=0;ii<rows.size();ii++){
	  for(unsigned int jj=0;jj<ii;jj++){
	    C(rows[jj],rows[ii]) += M(is[ii],3)*M(is[jj],3);
	  }
	}
	is.clear();
	rows.clear();
      }
    }
    is.push_back(i);
    rows.push_back((int)M(i,1));
  }
  for(int i=1;i<=nrows;i++){
    C(i,i)=1;
    for(int j=1;j<i;j++){
      C(i,j) -= (float)ncols * m(i)*m(j);
      if(s(i)!=0)
	C(i,j) /= s(i);
      if(s(j)!=0)
	C(i,j) /= s(j);
    }
  }

  C.Release();
  return C;
}


int main ( int argc, char **argv ){
  ccopsOptions& opts = ccopsOptions::getInstance();
  int success=opts.parse_command_line(argc,argv);
  if(!success) return -1;

  string ip=opts.inmatrix2.value();
  make_basename(ip);

 
  ColumnVector y1,r1,y2,r2;
  Matrix U;

  Matrix myOMmat = read_ascii_matrix(opts.ptxdir.value()+"/"+opts.inmatrix2.value()+".dot");
  // binarise
  if(opts.bin.value()>=0){
    for(int i=1;i<myOMmat.Nrows();i++){
      if(myOMmat(i,3)<(int)opts.bin.value()) myOMmat(i,3)=0;
      else myOMmat(i,3)=1;
    }
  }
  

  //Checking for and loading up Seed Coordinates
  string coordname=opts.ptxdir.value()+"/coords_for_"+ip;
  Matrix coordmat;
  bool   coordbool=false;
  ifstream fs1(coordname.c_str());
  if(fs1){
    coordmat  = read_ascii_matrix(coordname);
    coordbool = true;
  }
  else{
    cout<<"Seed Space Coordinate File Not present - Ignoring"<<endl;
  }

  //Checking For and Loading Up Tract coordinates
  string tractcoordname=opts.ptxdir.value()+"/tract_space_coords_for_"+ip;
  Matrix tractcoordmat;
  bool   tractcoordbool=false;
  ifstream fs2(coordname.c_str());
  if(fs2){
    tractcoordmat  = read_ascii_matrix(tractcoordname);
    tractcoordbool = true;
  }
  else{
    cout<<"Tract Space Coordinate File Not present - Ignoring"<<endl;
  }
 

  // If user specifies an exclusion mask in tract space. 
  // work out which columns in the matrix to remove 
  // This only works if there is a lookup matrix available
  vector<int> excl_cols;
  if(opts.excl_mask.value()!=""){
    volume<int>   lookup_tract;
    volume<float> excl;
    string        exname=opts.excl_mask.value();

    read_volume(lookup_tract,opts.ptxdir.value()+"/lookup_tractspace_"+ip);
    make_basename(exname);
    read_volume(excl,exname);
    if(!samesize(excl,lookup_tract)){
      cerr<<"Whoops - your exlusion mask does not appear to be "
	  <<"in the same space as your original low resolution mask - sorry"<<endl;
      return(-1);
    }

    for(int z=0;z<=excl.zsize();z++){
      for(int y=0;y<=excl.ysize();y++){
	for(int x=0;x<=excl.xsize();x++){
	  if(excl(x,y,z)==0)continue;
	  if(lookup_tract(x,y,z)==0)continue;

	  if(lookup_tract(x,y,z)<=tractcoordmat.Nrows()){
	    excl_cols.push_back(lookup_tract(x,y,z));
	  }
	  else{
	    cerr<<"Something a bit dodgy has happened here"<<endl;
	    cerr<<"Have you already run a reord_OM on this matrix"<<endl;
	    cerr<<"If so you can't use an exclusion mask as the"<<endl;
	    cerr<<"tractspace_lookup volume is not valid for this matrix"<<endl;
	    return(-1);
	  }

	}
      }
    }

  }

  if(opts.verbose.value())
    cout<<"remove zero rows and columns"<<endl;
  Matrix newcoordmat,newtractcoordmat;
  Matrix newOMmat;
  rem_zrowcol(myOMmat,excl_cols,
  	      coordmat,tractcoordmat,
	      coordbool,tractcoordbool,
	      newOMmat,newcoordmat,newtractcoordmat);
  
  int newnrows = (int)newOMmat( newOMmat.Nrows(),1 );

  string base=opts.obasename.value();
  make_basename(base);
  volume<float> outCCvol  (newnrows,newnrows,1);
  Matrix        outcoords (newcoordmat.Nrows(),coordmat.Ncols());

  if(opts.verbose.value())
    cout<<"Computing correlation"<<endl;
  SymmetricMatrix CtCt;
  CtCt = mycorrcoef(newOMmat);
  CtCt << CtCt + 1;


  // adding connexity constraint
  if(opts.verbose.value())
    cout<<"Adding Connexity constraint"<<endl;
  if(!coordbool){
    cerr<<"Requested connexitiy constraint, but couldn't locate coordinate file - Ignoring"<<endl;
  }
  else{
    add_connexity(CtCt,newcoordmat,opts.connexity.value());
  }

  if(opts.power.value()!=1){
    CtCt << pow(CtCt,opts.power.value());
  }
 
  if(!opts.reord1.value()){
     
    for(int j=0;j<outCCvol.ysize();j++){
      for(int i=0;i<outCCvol.xsize();i++){
	outCCvol(i,j,0)=CtCt(i+1,j+1);
      }
    }
      
    save_volume(outCCvol,opts.ptxdir.value()+"/CC_"+base);
    if(coordbool)
      write_ascii_matrix(newcoordmat,opts.ptxdir.value()+"/coords_for_"+base);

  }
  else{
    if(opts.verbose.value())
      cout<<"Starting First Reordering"<<endl;

    if(opts.scheme.value()=="spectral")
      spect_reord(CtCt,r1,y1);
    else if(opts.scheme.value()=="kmeans")
      kmeans_reord(CtCt,r1,y1,opts.nclusters.value());
    else if(opts.scheme.value()=="fuzzy")
      fuzzy_reord(CtCt,U,r1,y1,opts.nclusters.value());
    else{
      cerr << "unkown reordering scheme" << endl;
      return(-1);
    }
   
    if(opts.verbose.value())
      cout<<"Permuting seed CC matrix"<<endl;
    for(int j=0;j<outCCvol.ysize();j++){
      for(int i=0;i<outCCvol.xsize();i++){
	outCCvol(i,j,0)=CtCt((int)r1(i+1),(int)r1(j+1));
      }
    }
   
   
    if(coordbool){
      if(opts.verbose.value())
	cout<<"Permuting Seed Coordinates"<<endl;

      for(int i=0;i<outcoords.xsize();i++){
	outcoords.Row(i) << newcoordmat.Row(int(r1(i+1)));
      } 
    }

    if(opts.verbose.value())
      cout<<"Saving results"<<endl;
    write_ascii_matrix(r1,opts.ptxdir.value()+"/"+base+"r1");
    write_ascii_matrix(y1,opts.ptxdir.value()+"/"+base+"y1");
    save_volume(outCCvol,opts.ptxdir.value()+"/reord_CC_"+base);
    if(coordbool)
      write_ascii_matrix(outcoords,opts.ptxdir.value()+"/coords_for_reord_"+base);

    // propagate seed clustering onto tract space if requested
    if(opts.reord3.value()){
      if(opts.scheme.value() == "spectral"){
	cerr << "Warning: cannot propagate to tract space under this scheme." << endl;
	cerr << "will carry on ignoring this option"<<endl;
      }
      else{
	if(opts.verbose.value())
	  cout << "Propagating seed clustering into tract space" << endl;
	volume<float> opaths;
	read_volume(opaths,opts.ptxdir.value()+"/lookup_tractspace_"+ip);
	opaths = 0.0;

	// firstly, determine sum of matrix2 over different classes
	Matrix sumOverSeeds((int)y1.Maximum(),newtractcoordmat.Nrows());
	sumOverSeeds = 0.0;
	for(int i=1;i<newOMmat.Nrows();i++){
	  int r=(int)newOMmat(i,1),c=(int)newOMmat(i,2);
	  float v=newOMmat(i,3);
	  sumOverSeeds((int)y1(r),c) += v;
	}
	float minval = sumOverSeeds.Minimum();
	for(int j=1;j<=newtractcoordmat.Nrows();j++){
	  float maxval;int maxind;
	  maxval = sumOverSeeds.Column(j).Maximum1(maxind);
	  if(minval != maxval)
	    opaths((int)newtractcoordmat(j,1),
		   (int)newtractcoordmat(j,2),
		   (int)newtractcoordmat(j,3)) = maxind;
	}
	
	opaths.setDisplayMaximumMinimum(opaths.max(),0);
	save_volume(opaths,opts.ptxdir.value()+"/tract_space_propagated_"+base);
      }
    }

    // save clustering if kmeans used
    if(opts.scheme.value() == "kmeans" || opts.scheme.value()=="fuzzy"){
      if(opts.mask.value()!=""){
	volume<unsigned int> refvol;
	read_volume(refvol,opts.ptxdir.value()+"/fdt_paths");
	CSV mask(refvol);
	mask.load_rois(opts.mask.value());
	mask.reset_values();
	
	string type;int roiind,roiloc;
	for(int i=1;i<=outcoords.Nrows();i++){
	  roiind=outcoords(i,4);
	  roiloc=outcoords(i,5);
	  if(mask.isVol(roiind))
	    type="volume";
	  else
	    type="surface";
	  
	  mask.set_value(type,roiind,roiloc,(int)y1(i+1));
	}

	if(mask.Nrois()>1){
	  for(int i=0;i<mask.Nrois();i++)
	    mask.save_roi(i,opts.ptxdir.value()+"/reord_mask_"+num2str(i)+"_"+base);
	}
	else
	  mask.save_roi(i,opts.ptxdir.value()+"/reord_mask_"+base);



	// save memberships if fuzzy clustering used
	if(opts.scheme.value() == "fuzzy"){
	  CSV umask(refvol);
	  umask.load_rois(opts.mask.value());
	  umask.reset_values();

	  volume<float> umask;
	  umask.reinitialize(mask.xsize(),mask.ysize(),mask.zsize());
	  //OUT(U);
	  for(int cl=1;cl<=opts.nclusters.value();cl++){
	    umask=0;
	    for(int i=0;i<outcoords.xsize();i++){
	      //if((int)y1(i+1) == cl){
	      umask(outcoords(i,0,0),
		    outcoords(i,1,0),
		    outcoords(i,2,0)) = U((int)r1(i+1),cl);
	      //}
	    }
	    umask.setDisplayMaximumMinimum(1,0);
	    save_volume(umask,opts.ptxdir.value()+"/reord_membership_class"+num2str(cl)+"_"+base);
	  }
	}
      }
    }
    
  }

  if(opts.reord2.value()){
    if(opts.verbose.value())
      cout<<"Starting Second Reordering"<<endl;
    // need to change all this....
    SymmetricMatrix CC;
    CC << corrcoef(newOMmat);
    CC<<CC+1;

    if(opts.scheme.value()=="spectral")
      spect_reord(CC,r2,y2);
    else if(opts.scheme.value()=="kmeans")
      kmeans_reord(CC,r2,y2,opts.nclusters.value());
    else{
      cerr << "unkown reordering scheme" << endl;
      return(-1);
    }
 
    write_ascii_matrix(r2,opts.ptxdir.value()+"/"+base+"r2");
    write_ascii_matrix(y2,opts.ptxdir.value()+"/"+base+"y2");

    volume<int> outvol(newOMmat.Nrows(),newOMmat.Ncols(),1);
    volume<int> outtractcoords(newtractcoordmat.Nrows(),3,1);


    if(opts.verbose.value())
      cout<<"Permuting Matrix"<<endl;
    for(int j=0;j<outvol.ysize();j++){
      for(int i=0;i<outvol.xsize();i++){
	outvol(i,j,0)=(int)newOMmat((int)r1(i+1),(int)r2(j+1));
      }
    }

    if(tractcoordbool){
      if(opts.verbose.value())
	cout<<"Permuting Tract Coordinates"<<endl;
      for(int i=0;i<outtractcoords.xsize();i++){
	outtractcoords(i,0,0)=(int)newtractcoordmat(int(r2(i+1)),1);
	outtractcoords(i,1,0)=(int)newtractcoordmat(int(r2(i+1)),2);
	outtractcoords(i,2,0)=(int)newtractcoordmat(int(r2(i+1)),3);
      } 
    }
    save_volume(outvol,opts.ptxdir.value()+"/reord_"+base);
    save_volume(outtractcoords,opts.ptxdir.value()+"/tract_space_coords_for_reord_"+base);
  }

  
  if(opts.verbose.value())
    cout << "Done." << endl;
  return 0;
}
 


















