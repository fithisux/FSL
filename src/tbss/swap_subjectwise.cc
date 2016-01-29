/*  swap_dyadic_vectors.cc

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
#include "string.h"
#include "miscmaths/miscmaths.h"



using namespace Utilities;
using namespace std;
using namespace NEWIMAGE;
using namespace MISCMATHS;


string title="swap_subjectwise\n Reordering of the dyadic vectors and fsamples according to average inter-subject modal orientations";
string examples="swap_subjectwise -r <ListOfListOfDyads> -f <ListOfListOfFsamples>[-m mask -b <obase>]\n";


Option<bool> help(string("-h,--help"),false,
		       string("display this message"),
		       false,no_argument);
Option<bool> verbose(string("-v,--verbose"),false,
		       string("switch on diagnostic messages"),
		       false,no_argument);
Option<string> obasename(string("-b,--obasename"),string("swapped"),
		       string("output obasename [default=swapped]"),
		       false,requires_argument);
Option<string> maskfile(string("-m,--mask"),string(""),
		       string("filename of brain mask"),
		       true,requires_argument);
Option<string> dlist(string("-r,--dyads"),string(""),
		       string("list of list of dyads"),
		       true,requires_argument);
Option<string> flist(string("-f,--fmean"),string(""),
		       string("list of list of meanfsamples"),
		       true,requires_argument);
Option<float> xthresh(string("--xthresh"),0.1,
		      string("A.R.D. threshold - default=0.1"),
		      false,requires_argument);
Option<bool> avg(string("--averageonly"),false,
		      string("average only?"),
		      false,no_argument);

////////////////////////////////////////////////////////

void get_filenames(vector<string>& masks, const string& filename){
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



int main(int argc,char *argv[]){

  Tracer tr("main");
  OptionParser options(title,examples);

  try{
    options.add(help);
    options.add(verbose);
    options.add(obasename);
    options.add(maskfile);
    options.add(dlist);
    options.add(flist);
    options.add(xthresh);
    options.add(avg);

    options.parse_command_line(argc,argv);

    
    if ( (help.value()) || (!options.check_compulsory_arguments(true)) ){
      options.usage();
      exit(EXIT_FAILURE);
    }

    if(verbose.value())
      cout << "reading file names" << endl;

    vector<string> dyadsList,fmeanList;
  
    get_filenames(dyadsList,dlist.value());
    get_filenames(fmeanList,flist.value());

    if(dyadsList.size() != fmeanList.size()){
      cerr << "number of dyads and fmeans not equal" << endl;
      exit(1);
    }
    unsigned int ndyads = dyadsList.size();
    unsigned int nsubj = 0;

    vector< vector<string> > dyadspersubject(ndyads),fmeanpersubject(ndyads);
    for(unsigned int i=0;i<dyadsList.size();i++){
      get_filenames(dyadspersubject[i],dyadsList[i]);
      get_filenames(fmeanpersubject[i],fmeanList[i]);

      if(dyadspersubject[i].size() != fmeanpersubject[i].size()){
	cerr << "number of subjects not equal for " << i+1 << "-th entry" << endl;
	exit(1);
      }
      if(i==0)nsubj=dyadspersubject[i].size();
      else{
	if(dyadspersubject[i].size() != nsubj){
	  cerr << "number of subjects different for " << i+1 << "-th entry (at least)" << endl;
	  exit(1);
	}
      }
    }
    
    if(verbose.value())
      cout << "reading data" << endl;

    volume<float> mask;
    volume<int> vol2mat;
    Matrix mat2vol;
    read_volume(mask,maskfile.value());
    make_lut(mask,vol2mat,mat2vol);

    vector< vector< ColumnVector > >   fmeans(ndyads);
    vector< vector< Matrix > > dyads(ndyads);
    
    {
      volume<float>   tmp3D;ColumnVector tmpvec;
      volume4D<float> tmp4D;Matrix tmpmat;
      for(unsigned int i=0;i<ndyads;i++){
	if(verbose.value())
	  cout << "Vector " << i+1 << endl;
	for(unsigned int j=0;j<nsubj;j++){
	  if(verbose.value())
	    cout << "subject " << j+1 << endl;
	  read_volume(tmp3D,fmeanpersubject[i][j]);
	  read_volume4D(tmp4D,dyadspersubject[i][j]);
	  fill_vector(tmpvec,tmp3D,mat2vol);
	  fill_matrix(tmpmat,tmp4D,mat2vol);
	  
	  dyads[i].push_back(tmpmat);
	  fmeans[i].push_back(tmpvec);
	}
      }
    }

    vector< Matrix > mean_dyads(ndyads);


    if(verbose.value())
      cout << "FIRST PASS" << endl;
    
    // FIRST DETERMINE MEAN ORIENTATIONS ACROSS SUBJECTS
    ColumnVector V(3);
    Matrix U;
    SymmetricMatrix T(3);
    DiagonalMatrix D;
    for(unsigned int i=0;i<ndyads;i++){
      mean_dyads[i] = dyads[i][0];
      mean_dyads[i] = 0;


      for(int z=0;z<mask.zsize();z++)
	for(int y=0;y<mask.ysize();y++)
	  for(int x=0;x<mask.xsize();x++){
	    if(mask(x,y,z)==0)continue;

	    // read dyadic vectors from all subjects
	    T=0;
	    for(unsigned int j=0;j<nsubj;j++){
	      V << dyads[i][j](vol2mat(x,y,z),1) 
		<< dyads[i][j](vol2mat(x,y,z),2) 
		<< dyads[i][j](vol2mat(x,y,z),3);

	      T << T + V*V.t();
	    }

	    // extract principal component
	  
	    EigenValues(T,D,U);

	    mean_dyads[i](vol2mat(x,y,z),1) = U(1,3);
	    mean_dyads[i](vol2mat(x,y,z),2) = U(2,3);
	    mean_dyads[i](vol2mat(x,y,z),3) = U(3,3);
	  }
    }

    if(verbose.value())
      cout<<"saving modes"<<endl;
    {
      volume4D<float> tmp4D;
      for(unsigned int i=0;i<ndyads;i++){
	tmp4D.reinitialize(mask.xsize(),mask.ysize(),mask.zsize(),3);
	tmp4D=0;
	fill_volume4D(tmp4D,mean_dyads[i],mat2vol);
      save_volume4D(tmp4D,obasename.value()+"_meandyads"+num2str(i+1));
      }
    }

    if(avg.value()){
      return(0);
    }

    if(verbose.value())
      cout << "SECOND PASS" << endl;
    // NOW RE-ATTRIBUTE F-VALUES TO APPROPRIATE ORIENTATIONS
    Matrix Perm;
    Perm = MISCMATHS::perms(ndyads);


    for(int z=0;z<mask.zsize();z++)
      for(int y=0;y<mask.ysize();y++)
	for(int x=0;x<mask.xsize();x++){
	  if(mask(x,y,z)==0)continue;

	  // for each subject, calculate the score of each permutation and keep the strongest one
	  ColumnVector vsubj(3),vgrp(3);
	  for(unsigned int j=0;j<nsubj;j++){

	    // only do something when there are crossing fibres!
	    int nfibres = 0;
	    for(unsigned int i=0;i<ndyads;i++){
	      if(fmeans[i][j](vol2mat(x,y,z)) >= xthresh.value()){
	    nfibres++;
	    }
	    }
	    if(nfibres<2) continue;

	    float dotprod=0,tmpdot;//,fsubj;
	    int optperm=1;
	    for(int p=1;p<=Perm.Nrows();p++){
	      tmpdot=0;
	      for(unsigned int i=0;i<ndyads;i++){
		//fsubj = fmeans[i][j](x,y,z);
		vsubj << dyads[i][j](vol2mat(x,y,z),1) 
		      << dyads[i][j](vol2mat(x,y,z),2) 
		      << dyads[i][j](vol2mat(x,y,z),3);
		int ii = (int)Perm(p,i+1)-1;
		vgrp << mean_dyads[ii](vol2mat(x,y,z),1) 
		     << mean_dyads[ii](vol2mat(x,y,z),2) 
		     << mean_dyads[ii](vol2mat(x,y,z),3);
		tmpdot = max(tmpdot,abs(dot(vsubj,vgrp)));
	      }
	      //	      tmpdot /= float(ndyads);
	      if(tmpdot > dotprod){
		dotprod = tmpdot;
		optperm = p;
	      }

	    }

	    // Use those permutations!
	    Matrix newV(ndyads,3);
	    ColumnVector newF(ndyads);
	    for(unsigned int i=0;i<ndyads;i++){
	      newV(i+1,1) = dyads[(int)Perm(optperm,i+1)-1][j](vol2mat(x,y,z),1);
	      newV(i+1,2) = dyads[(int)Perm(optperm,i+1)-1][j](vol2mat(x,y,z),2);
	      newV(i+1,3) = dyads[(int)Perm(optperm,i+1)-1][j](vol2mat(x,y,z),3);
	      newF(i+1) = fmeans[(int)Perm(optperm,i+1)-1][j](vol2mat(x,y,z));
	      
	    }
	    for(unsigned int i=0;i<ndyads;i++){
	      dyads[i][j](vol2mat(x,y,z),1) = newV(i+1,1);
	      dyads[i][j](vol2mat(x,y,z),2) = newV(i+1,2);
	      dyads[i][j](vol2mat(x,y,z),3) = newV(i+1,3);
	      fmeans[i][j](vol2mat(x,y,z)) = newF(i+1);
	    }


	  }//end of subject loop


	}//end of volume loop


    if(verbose.value())
      cout<<"saving results"<<endl;
    volume<float> tmp3D;
    tmp3D.reinitialize(mask.xsize(),mask.ysize(),mask.zsize());
    copybasicproperties(mask,tmp3D);
    volume4D<float> tmp4D;
    tmp4D.reinitialize(mask.xsize(),mask.ysize(),mask.zsize(),3);
    copybasicproperties(mask,tmp4D[0]);
    for(unsigned int i=0;i<ndyads;i++){
      for(unsigned int j=0;j<nsubj;j++){
	tmp3D=0;
	fill_volume(tmp3D,fmeans[i][j],mat2vol);
	save_volume(tmp3D,obasename.value()+"_"+fmeanpersubject[i][j]);
	tmp4D=0;
	fill_volume4D(tmp4D,dyads[i][j],mat2vol);
	save_volume4D(tmp4D,obasename.value()+"_"+dyadspersubject[i][j]);
      }
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
