/*  Copyright (C) 2004 University of Oxford  */

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

using namespace std;
using namespace NEWIMAGE;
using namespace NEWMAT;
using namespace CCOPS;
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
  int numiter=200;
  if(data.Nrows() != (int)y.Nrows()){
    y.ReSize(data.Nrows());
  }
  int n = data.Nrows();
  int d = data.Ncols();
  Matrix means(d,k),newmeans(d,k);
  ColumnVector nmeans(k);
  means=0;
  nmeans=0;
  
  // // initialise random
//   vector< pair<float,int> > rindex(n);
//   randomise(rindex);
//   vector<pair<float,int> >::iterator riter;
//   int nn=0,cl=1,nperclass=(int)(float(n)/float(k));
//   for(riter=rindex.begin();riter!=rindex.end();++riter){
//     means.Column(cl) += data.Row((*riter).second).t();
//     nmeans(cl) += 1;
//     y((*riter).second)=cl;
//     nn++;
//     if(nn>=nperclass && cl<k){
//       nn=0;
//       cl++;
//     }
//   }
//   for(int m=1;m<=k;m++)
//     means.Column(m) /= nmeans(m);


// initialise with far-away trick
// start with a random class centre. then each new class centre is
  // as far as possible from the cog of the previous classes
  means.Column(1) = data.Row(MISCMATHS::round(rand()/float(RAND_MAX)*float(n-1))+1).t();
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
    mypair.first=y(i);
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
  int numiter     = 200;
  float fuzziness = 2;

  int n = data.Nrows();
  int d = data.Ncols();
  Matrix means(d,k),newmeans(d,k);
  ColumnVector nmeans(k);
  means=0;
  nmeans=0;
  
  // initialise with far-away trick
  // start with a random class centre. then each new class centre is
  // as far as possible from the cog of the previous classes
  means.Column(1) = data.Row(MISCMATHS::round(rand()/float(RAND_MAX)*float(n-1))+1).t();
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
  //OUT(u);

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
    mypair.first=y(i);
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



void rem_zrowcol(const Matrix& myOMmat,const Matrix& coordmat,const Matrix& tractcoordmat,const bool coordbool,const bool tractcoordbool,Matrix& newOMmat,Matrix& newcoordmat, Matrix& newtractcoordmat)
{
 
 vector<int> zerorows;
 vector<int> zerocols;
 int dimsum =0;



 cout<< "Checking for all zero rows"<<endl;
 for(int i=1;i<=myOMmat.Nrows();i++){
   dimsum=0;
   for(int j=1;j<=myOMmat.Ncols();j++){
     dimsum+=(int)myOMmat(i,j);
   }
   if(dimsum==0){zerorows.push_back(i);}
 }



 cout<< "Checking for all zero cols"<<endl;
 for(int j=1;j<=myOMmat.Ncols();j++){
   dimsum=0;
   for(int i=1;i<=myOMmat.Nrows();i++){
     dimsum+=(int)myOMmat(i,j);
   }
   if(dimsum==0){zerocols.push_back(j);}
 }


 newOMmat.ReSize(myOMmat.Nrows()-zerorows.size(),myOMmat.Ncols()-zerocols.size());
 if(coordbool){
   newcoordmat.ReSize(coordmat.Nrows()-zerorows.size(),3);
 }
 if(tractcoordbool){
   newtractcoordmat.ReSize(tractcoordmat.Nrows()-zerocols.size(),3);
 }

 int zrowcounter=0,zcolcounter=0,nzrowcounter=1,nzcolcounter=1;
 cout<<"Forming New Matrix"<<endl;
 for(int j=1;j<=myOMmat.Ncols();j++){
   zrowcounter=0;
   nzrowcounter=1;

   if(zerocols.size()>0){ //Are there any Zero Columns
     if(zerocols[zcolcounter]!=j){  // Only add a col if it's not the next zcol
       for(int i=1;i<=myOMmat.Nrows();i++){
	 if(zerorows.size()>0){ //Are there any zero rows?
	   if(zerorows[zrowcounter]!=i){ //Only add a row if it's not the next zrow
	     newOMmat(nzrowcounter,nzcolcounter)=myOMmat(i,j);
	     nzrowcounter++;
	   }
	   else{zrowcounter++;}
	 }
	 else{newOMmat(nzrowcounter,nzcolcounter)=myOMmat(i,j);nzrowcounter++;}//No Zero Rows
       }
       nzcolcounter++;
     }
     else{zcolcounter++;} //move onto next z col
   }
   
   else{  //No zero Columns
     for(int i=1;i<=myOMmat.Nrows();i++){
       if(zerorows.size()>0){
	 if(zerorows[zrowcounter]!=i){ //Only add a row if it's not the next zrow
	   newOMmat(nzrowcounter,nzcolcounter)=myOMmat(i,j);
	   nzrowcounter++;
	 }
	 else{zrowcounter++;}
       }
       else{newOMmat(nzrowcounter,nzcolcounter)=myOMmat(i,j);nzrowcounter++;}
     }
     nzcolcounter++;
   }

}


 if(coordbool){
   cout<<"Updating Seed Coordinates"<<endl;
   zrowcounter=0;nzrowcounter=1;
   if(zerorows.size()>0){//Are there any zero rows?
     for(int i=1;i<=coordmat.Nrows();i++){
       if(zerorows[zrowcounter]!=i){
	 newcoordmat(nzrowcounter,1)=coordmat(i,1);
	 newcoordmat(nzrowcounter,2)=coordmat(i,2);
	 newcoordmat(nzrowcounter,3)=coordmat(i,3);
	 nzrowcounter++;
       }
       else{zrowcounter++;}
     }
   } 
   else{//No zero Rows
     newcoordmat=coordmat;
   }
}

 if(tractcoordbool){
   cout<<"Updating Tract Coordinates"<<endl;
   zcolcounter=0;nzcolcounter=1;
   if(zerocols.size()>0){//Are there any zero cols?
     for(int i=1;i<=tractcoordmat.Nrows();i++){
       if(zerocols[zcolcounter]!=i){
	 newtractcoordmat(nzcolcounter,1)=tractcoordmat(i,1);
	 newtractcoordmat(nzcolcounter,2)=tractcoordmat(i,2);
	 newtractcoordmat(nzcolcounter,3)=tractcoordmat(i,3);
	 nzcolcounter++;
       }
       else{zcolcounter++;}
     }
   }
   else{//No zero Cols 
     newtractcoordmat=tractcoordmat;
   }
 }

}









void rem_cols(Matrix& myOMmat,Matrix& tractcoordmat,const bool tractcoordbool,const vector<int>& excl_cols)
{
 

 Matrix newOMmat,newtractcoordmat;
 newOMmat.ReSize(myOMmat.Nrows(),myOMmat.Ncols()-excl_cols.size());

 if(tractcoordbool){
   newtractcoordmat.ReSize(tractcoordmat.Nrows()-excl_cols.size(),3);
 }

 int zrowcounter=0,zcolcounter=0,nzcolcounter=1,nzrowcounter=1;
 vector<int> zerorows;

 for(int j=1;j<=myOMmat.Ncols();j++){
   zrowcounter=0;
   nzrowcounter=1;   

   if(excl_cols.size()>0){ //Are there any excl Columns
     if(excl_cols[zcolcounter]!=j){  // Only add a col if it's not the next zcol
       for(int i=1;i<=myOMmat.Nrows();i++){
	 if(zerorows.size()>0){ //Are there any excl rows?
	   if(zerorows[zrowcounter]!=i){ //Only add a row if it's not the next zrow
	     newOMmat(nzrowcounter,nzcolcounter)=myOMmat(i,j);
	     nzrowcounter++;
	   }
	   else{zrowcounter++;}
	 }
	 else{newOMmat(nzrowcounter,nzcolcounter)=myOMmat(i,j);nzrowcounter++;}//No Zero Rows
       }
       nzcolcounter++;
     }
     else{zcolcounter++;} //move onto next z col
   }
   else{  //No zero Columns
     for(int i=1;i<=myOMmat.Nrows();i++){
       if(zerorows.size()>0){
	 if(zerorows[zrowcounter]!=i){ //Only add a row if it's not the next zrow
	   newOMmat(nzrowcounter,nzcolcounter)=myOMmat(i,j);
	   nzrowcounter++;
	 }
	 else{zrowcounter++;}
       }
       else{newOMmat(nzrowcounter,nzcolcounter)=myOMmat(i,j);nzrowcounter++;}
     }
     nzcolcounter++;
   }
}

 
 if(tractcoordbool){
   zcolcounter=0;nzcolcounter=1;
   if(excl_cols.size()>0){//Are there any zero cols?
     for(int i=1;i<=tractcoordmat.Nrows();i++){
       if(excl_cols[zcolcounter]!=i){
	 newtractcoordmat(nzcolcounter,1)=tractcoordmat(i,1);
	 newtractcoordmat(nzcolcounter,2)=tractcoordmat(i,2);
	 newtractcoordmat(nzcolcounter,3)=tractcoordmat(i,3);
	 nzcolcounter++;
       }
       else{zcolcounter++;}
     }
   }
   else{//No zero Cols 
     newtractcoordmat=tractcoordmat;
   }
 }

 myOMmat = newOMmat;
 tractcoordmat = newtractcoordmat;
  
}

void add_connexity(SymmetricMatrix& CtCt,const Matrix& coord,const float p=.5){
  // compute CtCt range
  float r=CtCt.Minimum();
  float R=CtCt.Maximum();

  // compute distance matrix
  SymmetricMatrix D(coord.Nrows());
  for(int i=1;i<=coord.Nrows();i++)
    for(int j=1;j<=i;j++){
      D(i,j)=std::sqrt((coord.Row(i)-coord.Row(j)).SumSquare());
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



int main ( int argc, char **argv ){
  ccopsOptions& opts = ccopsOptions::getInstance();
  int success=opts.parse_command_line(argc,argv);
  if(!success) return 1;
  if(opts.inmatrix.value()=="" && opts.directory.value()==""){
    cerr << "Specify either input matrix or tractography directory" << endl;
    return 1;
  }
  string ip=opts.inmatrix.value();
  make_basename(ip);
  
  ColumnVector y1,r1,y2,r2;
  Matrix U;
  volume<int> myOM;
  volume<int> coordvol;
  volume<int> tractcoordvol;
  bool coordbool=false,tractcoordbool=false;
  //  read_volume(myOM,ip);
  read_volume(myOM,opts.directory.value()+"/"+opts.inmatrix.value());

  Matrix myOMmat(myOM.xsize(),myOM.ysize());
  Matrix mycoordmat,mytractcoordmat;
  Matrix newOMmat,newcoordmat,newtractcoordmat;
  
  for(int j=0;j<myOM.ysize();j++){
    for(int i=0;i<myOM.xsize();i++){
      if(opts.bin.value()==0)
	myOMmat(i+1,j+1)=float(myOM(i,j,0)); 
      else{
	if(myOM(i,j,0)>opts.bin.value()){
	  myOMmat(i+1,j+1)=1.0f;
	}
	else{
	  myOMmat(i+1,j+1)=0.0f;
	}
      }
    }
  }
  // write_ascii_matrix(newOMmat.t(),"preprecock");
  //Checking for and loading up Seed Coordinates

  string coordname=opts.directory.value()+"/coords_for_"+ip;
  if(fsl_imageexists(coordname)){
    read_volume(coordvol,coordname);
    coordbool=true;
    mycoordmat.ReSize(coordvol.xsize(),coordvol.ysize());
    for(int j=0;j<coordvol.ysize();j++){
      for(int i=0;i<coordvol.xsize();i++){
	mycoordmat(i+1,j+1)=coordvol(i,j,0);
      }
    }
  }
  else{
    cout<<"Seed Space Coordinate File Not present - Ignoring"<<endl;
  }

  //Checking For and Loading Up Tract coordinates
  string trcoordname=opts.directory.value()+"/tract_space_coords_for_"+ip;
  if(fsl_imageexists(trcoordname)){
    read_volume(tractcoordvol,trcoordname);
    tractcoordbool=true;
    mytractcoordmat.ReSize(tractcoordvol.xsize(),tractcoordvol.ysize());
    for(int j=0;j<tractcoordvol.ysize();j++){
      for(int i=0;i<tractcoordvol.xsize();i++){
	mytractcoordmat(i+1,j+1)=tractcoordvol(i,j,0);
      }
    }
  }
  else{
    cout<<"Tract Space Coordinate File Not present - Ignoring"<<endl;
  }
 


  // If user specifies an exclusion mask in tract space. 
  // work out which columns in the matrix to remove 
  // This only works if there is a lookup matrix available


  if(opts.excl_mask.value()!=""){
    volume<int> lookup_tract;
    volume<int> excl;
    read_volume(lookup_tract,opts.directory.value()+"/lookup_tractspace_"+ip);
    string exname=opts.excl_mask.value();
    make_basename(exname);
    read_volume(excl,exname);
    if(!samesize(excl,lookup_tract)){
      cerr<<"Whoops - your exlusion mask does not appear to be in the same space as your original low resolution mask - sorry"<<endl;
      return(-1);
    }
    vector<int> excl_cols;
    for(int k=0;k<=excl.zsize();k++){
      for(int j=0;j<=excl.ysize();j++){
	for(int i=0;i<=excl.xsize();i++){
	  if(excl(i,j,k)==1){
	    if(lookup_tract(i,j,k)!=0){
	     
	      if(lookup_tract(i,j,k)<=mytractcoordmat.Nrows()){
		excl_cols.push_back(lookup_tract(i,j,k)+1);
	      }
	      else{
		cerr<<"Something a bit dodgy has happened here"<<endl;
		cerr<<"Have you already run a reord_OM on this matrix"<<endl;
		cerr<<"If so you can't use an exclusion mask as the"<<endl;
		cerr<<"tractspace_lookup volume is not valid for this matrix"<<endl;
	      }
	     
	     
	    }
	  }
	}
      }
    }

    rem_cols(myOMmat,mytractcoordmat,tractcoordbool,excl_cols);
  }

 
  rem_zrowcol(myOMmat,mycoordmat,mytractcoordmat,coordbool,tractcoordbool,newOMmat,newcoordmat,newtractcoordmat);
  //   cerr<<"NOW"<<endl;
  //   cerr<<myOMmat.MaximumAbsoluteValue()<<endl;
  //   cerr<<newOMmat.MaximumAbsoluteValue()<<endl;
 
  //write_ascii_matrix("ncm",newcoordmat);
  // write_ascii_matrix("nctm",newtractcoordmat);



  string base=opts.basename.value();
  make_basename(base);
  volume<float> outCCvol(newOMmat.Nrows(),newOMmat.Nrows(),1);
  volume<int> outcoords(newcoordmat.Nrows(),3,1);



  cout<<"Computing correlation"<<endl;
  SymmetricMatrix CtCt;
  CtCt << corrcoef(newOMmat.t());
  CtCt << CtCt+1;

  // adding connexity constraint
  if(!coordbool){
    cout<<"WARNING !! No coordinates provided. I cannot apply any connexity constraint."<<endl;
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
   
   
    if(coordbool){
      for(int i=0;i<outcoords.xsize();i++){
	outcoords(i,0,0)=(int)newcoordmat(i+1,1);
	outcoords(i,1,0)=(int)newcoordmat(i+1,2);
	outcoords(i,2,0)=(int)newcoordmat(i+1,3);
      } 
    }
   
    save_volume(outCCvol,opts.directory.value()+"/CC_"+base);
    save_volume(outcoords,opts.directory.value()+"/coords_for_"+base);

  }
  else{
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
   
    cout<<"Permuting seed CC matrix"<<endl;
    for(int j=0;j<outCCvol.ysize();j++){
      for(int i=0;i<outCCvol.xsize();i++){
	outCCvol(i,j,0)=CtCt((int)r1(i+1),(int)r1(j+1));
      }
    }
   
   
    if(coordbool){
      cout<<"Permuting Seed Coordinates"<<endl;
//       OUT(r1.Maximum());
//       OUT(r1.Minimum());
//       OUT(newcoordmat.Ncols());
//       OUT(newcoordmat.Nrows());

      for(int i=0;i<outcoords.xsize();i++){
	outcoords(i,0,0)=(int)newcoordmat(int(r1(i+1)),1);
	outcoords(i,1,0)=(int)newcoordmat(int(r1(i+1)),2);
	outcoords(i,2,0)=(int)newcoordmat(int(r1(i+1)),3);
      } 
    }
    cout<<"Saving results"<<endl;
    write_ascii_matrix(r1,opts.directory.value()+"/"+base+"r1");
    write_ascii_matrix(y1,opts.directory.value()+"/"+base+"y1");
    save_volume(outCCvol,opts.directory.value()+"/reord_CC_"+base);
    save_volume(outcoords,opts.directory.value()+"/coords_for_reord_"+base);

    // propagate seed clustering onto tract space if requested
    if(opts.reord3.value()){
      if(opts.scheme.value() == "spectral"){
	cerr << "Warning: cannot propagate to tract space under this scheme." << endl;
	cerr << "will carry on ignoring this option"<<endl;
      }
      else{
	cout << "Propagating seed clustering into tract space" << endl;
	volume<float> opaths;
	read_volume(opaths,opts.directory.value()+"/lookup_tractspace_"+ip);
	opaths = 0.0;

	// firstly, determine sum of matrix2 over different classes
	Matrix sumOverSeeds(y1.Maximum(),newtractcoordmat.Nrows());
	sumOverSeeds = 0.0;
	for(int j=1;j<=newtractcoordmat.Nrows();j++){
	  for(int i=1;i<=y1.Nrows();i++){
	    sumOverSeeds(y1(i),j) += newOMmat(r1(i),j);
	  }

	  float maxval,minval;
	  int maxind;
	  minval = sumOverSeeds.Minimum();
	  maxval = sumOverSeeds.Column(j).Maximum1(maxind);
	  if(minval != maxval)
	    opaths(newtractcoordmat(j,1),
		   newtractcoordmat(j,2),
		   newtractcoordmat(j,3)) = maxind;
	}

	
	opaths.setDisplayMaximumMinimum(opaths.max(),0);
	save_volume(opaths,opts.directory.value()+"/tract_space_propagated_"+base);
      }
    }

    // save clustering if kmeans used
    if(opts.scheme.value() == "kmeans" || opts.scheme.value()=="fuzzy"){
      volume<int> mask;
      if(opts.mask.value() == "")
	read_volume(mask,opts.directory.value()+"/fdt_paths");
      else
	read_volume(mask,opts.mask.value());
      mask = 0;
      for(int i=0;i<outcoords.xsize();i++){
	mask(outcoords(i,0,0),
	     outcoords(i,1,0),
	     outcoords(i,2,0)) = (int)y1(i+1);
      }
      mask.setDisplayMaximumMinimum(y1.Maximum(),0);
      save_volume(mask,opts.directory.value()+"/reord_mask_"+base);
      
      // save memberships if fuzzy clustering used
      if(opts.scheme.value() == "fuzzy"){
	volume<float> umask;
	umask.reinitialize(mask.xsize(),mask.ysize(),mask.zsize());
	copybasicproperties(mask,umask);
	//OUT(U);
	for(int cl=1;cl<=opts.nclusters.value();cl++){
	  umask=0;
	  for(int i=0;i<outcoords.xsize();i++){
	    //if((int)y1(i+1) == cl){
	      umask(outcoords(i,0,0),
		    outcoords(i,1,0),
		    outcoords(i,2,0)) = U(r1(i+1),cl);
	      //}
	  }
	  umask.setDisplayMaximumMinimum(1,0);
	  save_volume(umask,opts.directory.value()+"/reord_membership_class"+num2str(cl)+"_"+base);
	}
	
      }
	
    }
    
  }

  if(opts.reord2.value()){
    cout<<"Starting Second Reordering"<<endl;
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
 
    write_ascii_matrix(r2,opts.directory.value()+"/"+base+"r2");
    write_ascii_matrix(y2,opts.directory.value()+"/"+base+"y2");

    volume<int> outvol(newOMmat.Nrows(),newOMmat.Ncols(),1);
    volume<int> outtractcoords(newtractcoordmat.Nrows(),3,1);


    cout<<"Permuting Matrix"<<endl;
    for(int j=0;j<outvol.ysize();j++){
      for(int i=0;i<outvol.xsize();i++){
	outvol(i,j,0)=(int)newOMmat((int)r1(i+1),(int)r2(j+1));
      }
    }

    if(tractcoordbool){
      cout<<"Permuting Tract Coordinates"<<endl;
      for(int i=0;i<outtractcoords.xsize();i++){
	outtractcoords(i,0,0)=(int)newtractcoordmat(int(r2(i+1)),1);
	outtractcoords(i,1,0)=(int)newtractcoordmat(int(r2(i+1)),2);
	outtractcoords(i,2,0)=(int)newtractcoordmat(int(r2(i+1)),3);
      } 
    }
    save_volume(outvol,opts.directory.value()+"/reord_"+base);
    save_volume(outtractcoords,opts.directory.value()+"/tract_space_coords_for_reord_"+base);
  }

  cout << "Done." << endl;
  return 0;
}
 


















