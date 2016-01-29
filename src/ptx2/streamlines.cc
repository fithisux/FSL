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

#include "streamlines.h"
#include "warpfns/fnirt_file_reader.h"
#include "warpfns/warpfns.h"


  //***************************************** 
  //************* MatCell_cmpr ************** 
  void MatCell_cmpr::decode(int64_t incode, int& nsamples, int& fibcnt1, int& fibcnt2, float& length_tot) const{
    //undo the coding for incode
    //code = two32*fibre_count + mult*mult*fibre_prop1 + mult*fibre_prop2 + length_val;
    int64_t fibre_count, fibre_prop1, fibre_prop2, length_val, two32=(1LL<<32), mult=1001, multmult=mult*mult;    
    
    fibre_count = incode / two32;
    incode = incode % two32;
    fibre_prop1 = incode / multmult;
    incode = incode % multmult;
    fibre_prop2 = incode / mult;
    incode = incode % mult;
    length_val = incode;

    fibcnt1=(int)MISCMATHS::round((float)(fibre_prop1*0.001*fibre_count));
    fibcnt2=(int)MISCMATHS::round((float)(fibre_prop2*0.001*fibre_count));
    length_tot=float(length_val*fibre_count);
    nsamples=(int)fibre_count;
  } 

  void MatCell_cmpr::encode(const int nsamples, const int fibcnt1, const int fibcnt2, const float length_tot){
    //store new coding in code2
    int64_t  fibre_count, fibre_prop1, fibre_prop2, length_val, two32=(1LL<<32), mult=1001;    
    // fibre_prop1, fibre_prop2 and length_val ***MUST*** BE WITHIN 0 and 1000 INCLUSIVE

    fibre_count = (int64_t)MIN((int64_t)nsamples,two32-1);
    //Notice that every time we decode and encode there are rounding errors because of the following functions.
    //We considered using the absolute values instead of proportions or average distances, but then the dynamic range of the values we can code is reduced. 
    fibre_prop1 = (int64_t)MIN((int64_t)(MISCMATHS::round(float(fibcnt1)/float(nsamples)*1000)),1000); 
    fibre_prop2 = (int64_t)MIN((int64_t)(MISCMATHS::round(float(fibcnt2)/float(nsamples)*1000)),1000);
    fibre_prop2 = (int64_t)MIN(fibre_prop2, (int64_t)(1000-fibre_prop1)); //Avoid sum over 1000 due to rounding errors
    length_val  = (int64_t)MIN((int64_t)(MISCMATHS::round((nsamples!=0?length_tot/float(nsamples):0.0))),1000); 
    code2 = two32*fibre_count + mult*mult*fibre_prop1 + mult*fibre_prop2 + length_val;
  }

    
  void MatCell_cmpr::add_one(float dist,int fib){
    int nsamples, fibcnt1, fibcnt2; float length_tot;
    
    //Undo the coding  for current value
    decode(code2,nsamples, fibcnt1, fibcnt2, length_tot);
    
    //Update Values
    if (fib==1)
      fibcnt1+=1;
    if (fib==2)
      fibcnt2+=1;
    length_tot+=dist;
    nsamples+=1;   

    //Code again
    encode(nsamples, fibcnt1, fibcnt2, length_tot);
  }


  void MatCell_cmpr::add_n(float dist,vector<float> props,int n){
    int nsamples, fibcnt1, fibcnt2; float length_tot;

    if(props.size()!=2){
      cerr<<"MatCell_cmpr::add_n:Only valid with 3 fibres for now"<<endl;
      exit(1);
    }
        
    //Undo the coding  for current value
    decode(code2, nsamples, fibcnt1, fibcnt2, length_tot);
    
    //Update Values
    fibcnt1+=(int)MISCMATHS::round(props[0]*n);
    fibcnt2+=(int)MISCMATHS::round(props[1]*n);
    length_tot+=(dist*n);
    nsamples+=n;   

    //Code again
    encode(nsamples, fibcnt1, fibcnt2, length_tot);
  }
  

  void MatCell_cmpr::add_n(int64_t newcode2){
    int nsamples, fibcnt1, fibcnt2; float length_tot;
    int Newnsamples, Newfibcnt1, Newfibcnt2; float Newlength_tot;
    
    //Undo the coding  for current value
    decode(code2, nsamples, fibcnt1, fibcnt2, length_tot);      
    
    //Undo the coding  for new value
    decode(newcode2, Newnsamples, Newfibcnt1, Newfibcnt2, Newlength_tot);      
    
    //Update the values by adding the New ones
    fibcnt1+=Newfibcnt1;
    fibcnt2+=Newfibcnt2;
    length_tot+=Newlength_tot;
    nsamples+=Newnsamples;

    //Perform coding again 
    encode(nsamples, fibcnt1, fibcnt2, length_tot);
  }



  //************************************** 
  //************* SpMat_HCP ************** 
  //Constructor from a file
  SpMat_HCP::SpMat_HCP(unsigned int m, unsigned int n,const string& basename):SpMat<MatCell_cmpr>::SpMat(m,n){
      string extension="mtx";
      string file1=basename+"1."+extension;
      string file2=basename+"2."+extension;
      string file3=basename+"3."+extension;
      ifstream fs1,fs2,fs3;
      // FIRST FILE (text file of matrix dimensions)
      try{ fs1.open(file1.c_str()); }
      catch(...) {
         throw SpMatHCPException("SpMat_HCP::SpMat_HCP(unsigned int m, unsigned int n,const string& basename): cannot open file1 for reading");
      }
      int nrows, ncols;
      fs1 >> nrows;
      fs1 >> ncols;
      fs1.close();

      if(nrows!=(int)Nrows() || ncols!=(int)Ncols()){
      	cerr<<"SpMat_HCP::SpMat_HCP(unsigned int m, unsigned int n,const string& basename): Incompatible matrix dimensions"<<endl;
	exit(1);
      }

      // SECOND FILE (binary file of lengths)
      try{ fs2.open(file2.c_str(), ios::in | ios::binary); }
      catch(...) {
         throw SpMatHCPException("SpMat_HCP::SpMat_HCP(unsigned int m, unsigned int n,const string& basename): cannot open file2 for reading");
      }
      ColumnVector colsize(ncols);
      for(int c=1; c <= ncols; c++) {
	int64_t sz;
	fs2.read((char*)&sz,sizeof(sz));
	colsize(c)=sz;    
      }
      fs2.close();
      
      // THIRD FILE (binary file of contents, integer coded)
      try{ fs3.open(file3.c_str(), ios::in | ios::binary); }
      catch(...) {
         throw SpMatHCPException("SpMat_HCP::SpMat_HCP(unsigned int m, unsigned int n,const string& basename): cannot open file3 for reading");
      }
      int64_t code1, code2;
      
      for(int c=0; c < ncols; c++) {
	for (unsigned int s=1; s <= colsize(c+1); s++) {
	  fs3.read((char*)&code1,sizeof(code1));
	  fs3.read((char*)&code2,sizeof(code2));
	  unsigned int r = code1; 
	  AddToTraj(r+1,c+1,code2);
	}
      }
      fs3.close();
    }
 

  // Trajectory file writer
  int SpMat_HCP::SaveTrajFile(const string& basename)const
  {
    if ( (basename.size()<1) ) return -1;
    string extension="mtx";
    string file1=basename+"1."+extension;
    string file2=basename+"2."+extension;
    string file3=basename+"3."+extension;
    ofstream fs1, fs2, fs3;
    // FIRST FILE (text file of matrix dimensions)
    try{ fs1.open(file1.c_str(), ios::out | ios::binary); }
    catch(...) {
      throw SpMatHCPException("SpMat_HCP::SaveTrajFile(const string& basename): cannot open file1 for writing");
    } 
    fs1 << Nrows() << endl;
    fs1 << Ncols() << endl;
    fs1.close();

    // SECOND FILE (binary file of lengths)
    try{ fs2.open(file2.c_str(), ios::out | ios::binary); }
    catch(...) {
      throw SpMatHCPException("SpMat_HCP::SaveTrajFile(const string& basename): cannot open file2 for writing");
    }
    for(unsigned int c=0; c < Ncols(); c++) {
      int64_t sz = get_ri(c).size();  
      fs2.write((char*)&sz,sizeof(sz));
    }
    fs2.close();

    // THIRD FILE (binary file of contents, integer coded)
    try{ fs3.open(file3.c_str(), ios::out | ios::binary); }
    catch(...) {
      throw SpMatHCPException("SpMat_HCP::SaveTrajFile(const string& basename): cannot open file3 for writing");
    }
    int64_t code1, code2;
    
    for(unsigned int c=0; c < Ncols(); c++) {
      if(get_ri(c).size()){
	const std::vector<unsigned int>&    ri = get_ri(c);
	const std::vector<MatCell_cmpr>&   val = get_val(c);
	for (unsigned int r=0; r<ri.size(); r++) { 	
	  code1 = ri[r];
	  code2 = val[r].getcode2();
	  fs3.write((char*)&code1,sizeof(code1));
	  fs3.write((char*)&code2,sizeof(code2));
	}
      }
    }  
    fs3.close();
    return 0;
  }

/*

  //Old SpMat_HCP with MattCell (less rounding errors, but takes ~5 times as much memory and ~40% more execution time
  //Constructor from a file
  SpMat_HCP::SpMat_HCP(unsigned int m, unsigned int n,const string& basename):SpMat<MatCell>::SpMat(m,n){
      string extension="mtx";
      string file1=basename+"1."+extension;
      string file2=basename+"2."+extension;
      string file3=basename+"3."+extension;
      ifstream fs1,fs2,fs3;
      // FIRST FILE (text file of matrix dimensions)
      try{ fs1.open(file1.c_str()); }
      catch(...) {
         throw SpMatHCPException("SpMat_HCP::SpMat_HCP(unsigned int m, unsigned int n,const string& basename): cannot open file1 for reading");
      }
      int nrows, ncols;
      fs1 >> nrows;
      fs1 >> ncols;
      fs1.close();

      if(nrows!=(int)Nrows() || ncols!=(int)Ncols()){
      	cerr<<"SpMat_HCP::SpMat_HCP(unsigned int m, unsigned int n,const string& basename): Incompatible matrix dimensions"<<endl;
	exit(1);
      }

      // SECOND FILE (binary file of lengths)
      try{ fs2.open(file2.c_str(), ios::in | ios::binary); }
      catch(...) {
         throw SpMatHCPException("SpMat_HCP::SpMat_HCP(unsigned int m, unsigned int n,const string& basename): cannot open file2 for reading");
      }
      ColumnVector colsize(ncols);
      for(int c=1; c <= ncols; c++) {
	int64_t sz;
	fs2.read((char*)&sz,sizeof(sz));
	colsize(c)=sz;    
      }
      fs2.close();
      
      // THIRD FILE (binary file of contents, integer coded)
      try{ fs3.open(file3.c_str(), ios::in | ios::binary); }
      catch(...) {
         throw SpMatHCPException("SpMat_HCP::SpMat_HCP(unsigned int m, unsigned int n,const string& basename): cannot open file3 for reading");
      }
      int64_t code1, code2;
      int64_t two32=(1LL<<32), mult=1001;    
      
      for(int c=0; c < ncols; c++) {
	for (unsigned int s=1; s <= colsize(c+1); s++) {
	  // fibre_prop1, fibre_prop2 and length_val ***MUST*** BE WITHIN 0 and 1000 INCLUSIVE	  
	  int64_t fibre_count, fibre_prop1, fibre_prop2, length_val;
	  fs3.read((char*)&code1,sizeof(code1));
	  fs3.read((char*)&code2,sizeof(code2));
	  //Undo this coding: code2 = two32*fibre_count + mult*mult*fibre_prop1 + mult*fibre_prop2 + length_val;
	  fibre_count = code2 / two32;
	  code2 = code2 % two32;
	  fibre_prop1 = code2 / (mult*mult);
	  code2 = code2 % (mult*mult);
	  fibre_prop2 = code2 / mult;
	  code2 = code2 % mult;
	  length_val = code2;
	  // Now undo the integer coding of fibre proportions
	  float fprop1=fibre_prop1 / 1000.0;
	  float fprop2=fibre_prop2 / 1000.0;

	  // Uncode the row 
	  unsigned int r = code1;
	  // update matrix
	  vector<float> props(2);
	  props[0]=fprop1;props[1]=fprop2;
	  AddToTraj(r+1,c+1,(float)length_val,props,fibre_count);
	}
      }
      fs3.close();
    }
 

  // Trajectory file writer
  int SpMat_HCP::SaveTrajFile(const string& basename)const
  {
    if ( (basename.size()<1) ) return -1;
    string extension="mtx";
    string file1=basename+"1."+extension;
    string file2=basename+"2."+extension;
    string file3=basename+"3."+extension;
    ofstream fs1, fs2, fs3;
    // FIRST FILE (text file of matrix dimensions)
    try{ fs1.open(file1.c_str(), ios::out | ios::binary); }
    catch(...) {
      throw SpMatHCPException("SpMat_HCP::SaveTrajFile(const string& basename): cannot open file1 for writing");
    } 
    fs1 << Nrows() << endl;
    fs1 << Ncols() << endl;
    fs1.close();

    // SECOND FILE (binary file of lengths)
    try{ fs2.open(file2.c_str(), ios::out | ios::binary); }
    catch(...) {
      throw SpMatHCPException("SpMat_HCP::SaveTrajFile(const string& basename): cannot open file2 for writing");
    }
    for(unsigned int c=0; c < Ncols(); c++) {
      int64_t sz = get_ri(c).size();  
      fs2.write((char*)&sz,sizeof(sz));
    }
    fs2.close();

    // THIRD FILE (binary file of contents, integer coded)
    try{ fs3.open(file3.c_str(), ios::out | ios::binary); }
    catch(...) {
      throw SpMatHCPException("SpMat_HCP::SaveTrajFile(const string& basename): cannot open file3 for writing");
    }
    int64_t code1, code2;
    int64_t two32=(1LL<<32), mult=1001;    
    
    for(unsigned int c=0; c < Ncols(); c++) {
      if(get_ri(c).size()){
	const std::vector<unsigned int>&    ri = get_ri(c);
	const std::vector<MatCell>&         val = get_val(c);
	for (unsigned int r=0; r<ri.size(); r++) { 	
	  code1 = ri[r];
	  int64_t fibre_count, fibre_prop1, fibre_prop2, length_val;
	  // fibre_prop1, fibre_prop2 and length_val ***MUST*** BE WITHIN 0 and 1000 INCLUSIVE	  
	  fibre_count = (int64_t)MIN((int64_t)val[r].get_nsamples(),two32-1);
	  fibre_prop1 = (int64_t)MIN((int64_t)(MISCMATHS::round(val[r].get_fibprop(1)*1000)),1000);
	  fibre_prop2 = (int64_t)MIN((int64_t)(MISCMATHS::round(val[r].get_fibprop(2)*1000)),1000);
	  length_val  = (int64_t)MIN((int64_t)(MISCMATHS::round(val[r].get_avg_length())),1000); 

	  code2 = two32*fibre_count + mult*mult*fibre_prop1 + mult*fibre_prop2 + length_val;
	  fs3.write((char*)&code1,sizeof(code1));
	  fs3.write((char*)&code2,sizeof(code2));
	}
      }
    }  
    fs3.close();
     
    return 0;
    
  }
*/

namespace TRACT{

  void read_ascii_files(const string& filename,vector<string>& content){
    ifstream fs(filename.c_str());
    string tmp;
    if(fs){
      fs>>tmp;
      do{
	content.push_back(tmp);
	fs>>tmp;
      }while(!fs.eof());
    }
    else{
      cerr<<"TRACT::read_ascii_files: "<<filename<<" does not exist"<<endl;
      exit(0);
    }
  }
  void write_matrix_as_volume(const Matrix& M,const string& filename){
    volume<float> tmp(M.Nrows(),M.Ncols(),1);
    for(int i=1;i<=M.Nrows();i++)
      for(int j=1;j<=M.Ncols();j++)
	tmp(i-1,j-1,0)=M(i,j);
    save_volume(tmp,filename);    
  }

  ColumnVector mean_sph_pol(ColumnVector& A, ColumnVector& B){
    // A and B contain th, ph f. 
    float th,ph;
    ColumnVector rA(3), rB(3);
    
    rA << (sin(A(1))*cos(A(2))) << (sin(A(1))*sin(A(2))) << (cos(A(1)));
    rB << (sin(B(1))*cos(B(2))) << (sin(B(1))*sin(B(2))) << (cos(B(1)));
    
    if(sum(SP(rA,rB)).AsScalar()>0)
      cart2sph((rA+rB)/2,th,ph);
    else
      cart2sph((rA-rB)/2,th,ph);
    
    A(1)=th; A(2)=ph;
    return A;
  }
  
  void imgradient(const volume<float>& im,volume4D<float>& grad){
    
    grad.reinitialize(im.xsize(),im.ysize(),im.zsize(),3);
    copybasicproperties(im,grad[0]);
    
    int fx,fy,fz,bx,by,bz;
    float dx,dy,dz; 
    for (int z=0; z<grad.zsize(); z++){
      fz = z ==(grad.zsize()-1) ? 0 :  1;
      bz = z == 0              ? 0 : -1;
      dz = (fz==0 || bz==0)    ? 1.0 :  2.0;
      for (int y=0; y<grad.ysize(); y++){
	fy = y ==(grad.ysize()-1) ? 0 :  1;
	by = y == 0              ? 0 : -1;
	dy = (fy==0 || by==0)    ? 1.0 :  2.0;
	for (int x=0; x<grad.xsize(); x++){
	  fx = x ==(grad.xsize()-1) ? 0 :  1;
	  bx = x == 0              ? 0 : -1;
	  dx = (fx==0 || bx==0)    ? 1.0 :  2.0;
	  grad[0](x,y,z) = (im(x+fx,y,z) - im(x+bx,y,z))/dx;
	  grad[1](x,y,z) = (im(x,y+fy,z) - im(x,y+by,z))/dy;
	  grad[2](x,y,z) = (im(x,y,z+fz) - im(x,y,z+bz))/dz;
	}
      }
    }
    
  }
  bool compare_Vertices(const pair<int,infoVertex> &x, const pair<int,infoVertex> &y){
    if(x.first<y.first) return true;
    if(x.first>y.first) return false;
    if(x.second.mesh<y.second.mesh) return true;
    if(x.second.mesh>y.second.mesh) return false;
    if(x.second.triangle<y.second.triangle) return true;
    if(x.second.triangle>y.second.triangle) return false;
    return (x.second.value<y.second.value);
  }
  void make_unique(vector<int>& x,vector<float>& y){
    // sort
    vector< pair<int,int> > p;
    for(unsigned int i=0;i<x.size();i++){
      pair<int,int> pp;
      pp.first=x[i];
      pp.second=i;
      p.push_back(pp);
    }
    sort(p.begin(),p.end());
    vector<float> ycp(y.begin(),y.end());
    for(unsigned int i=0;i<x.size();i++){
      x[i]=p[i].first;
      y[i]=ycp[ p[i].second ];
    }
    // unique
    vector<int> xx;vector<float> yy;
    for(unsigned int i=0;i<x.size();i++){
      if(i>0){
	if( x[i]==x[i-1] )continue;
      }
      xx.push_back(x[i]);
      yy.push_back(y[i]);
    }
    x=xx;
    y=yy;
  }

  void make_unique(vector<int>& x){
    sort(x.begin(),x.end());
    vector<int> xx;
    for(unsigned int i=0;i<x.size();i++){
      if(i>0){
	if( x[i]==x[i-1] )continue;
      }
      xx.push_back(x[i]);      
    }
    x=xx;
  }

  void make_unique(vector< pair<int,float> >&x){
    sort(x.begin(),x.end());
    vector< pair<int,float> > xx;
    for(unsigned int i=0;i<x.size();i++){
      if(i>0){
	if( x[i].first==x[i-1].first )continue;
      }
      xx.push_back(x[i]);
    }
    x=xx;
  }
 
  void make_unique(vector< pair<int,infoVertex> >&x){
    sort(x.begin(),x.end(),compare_Vertices);
    vector< pair<int,infoVertex> > xx;
    for(unsigned int i=0;i<x.size();i++){
      if(i>0){
	if( x[i].first==x[i-1].first && x[i].second.mesh==x[i-1].second.mesh && x[i].second.triangle==x[i-1].second.triangle )continue;
      }
      xx.push_back(x[i]);
    }
    x=xx;
  }

  Streamliner::Streamliner(const CSV& seeds):opts(probtrackxOptions::getInstance()),
					     logger(LogSingleton::getInstance()),					     					     
					     vols(opts.usef.value()),
					     m_seeds(seeds)
  {    

    // the tracking mask - no tracking outside of this
    read_volume(m_mask,opts.maskfile.value());
    
    // particle
    m_part.initialise(0,0,0,0,0,0,opts.steplength.value(),
		      m_mask.xdim(),m_mask.ydim(),m_mask.zdim(),
		      false);
    
    // no-integrity
    if(opts.skipmask.value()!="") 
      read_volume(m_skipmask,opts.skipmask.value());
    
    // loop checking box
    m_lcrat=5; // hard-coded! box is five times smaller than brain mask
    if(opts.loopcheck.value()){
      m_loopcheck.reinitialize(int(ceil(m_mask.xsize()/m_lcrat)+1),
			       int(ceil(m_mask.ysize()/m_lcrat)+1),
			       int(ceil(m_mask.zsize()/m_lcrat)+1),3);
      m_loopcheck=0;
    }

    // by default, assume there is no surface (speed)
    m_surfexists=false;
    if(m_seeds.nSurfs()>0){surfexists();}

    // exclusion mask
    // now in CSV format
    if(opts.rubbishfile.value()!=""){
      load_rubbish(opts.rubbishfile.value());   
    }
    
    // stopping mask
    // now in CSV format
    if(opts.stopfile.value()!=""){
      load_stop(opts.stopfile.value());
    }

    // Walk-through stopping mask
    if(opts.wtstopfiles.value()!=""){
      load_wtstopmasks(opts.wtstopfiles.value());
    }

     
    // waymasks in CSV format
    if(opts.waypoints.set()){
      load_waymasks(opts.waypoints.value());      
      m_waycond=opts.waycond.value();
    }
    
    // choose preferred fibre within this mask (prior knowledge)
    if(opts.prefdirfile.value()!=""){
      read_volume4D(m_prefdir,opts.prefdirfile.value());
    }    
    
    // local curvature threshold
    if(opts.loccurvthresh.value()!=""){
      read_volume(m_loccurvthresh,opts.loccurvthresh.value());
    }

    // Allow for either matrix transform (12dof affine) or nonlinear (warpfield)
    m_Seeds_to_DTI = IdentityMatrix(4);
    m_DTI_to_Seeds = IdentityMatrix(4);
    m_rotdir       = IdentityMatrix(3);
    
    m_IsNonlinXfm = false;   
    if(opts.seeds_to_dti.value()!=""){
      if(!fsl_imageexists(opts.seeds_to_dti.value())){//presumably ascii file provided
	m_Seeds_to_DTI = read_ascii_matrix(opts.seeds_to_dti.value());
	m_DTI_to_Seeds = m_Seeds_to_DTI.i();
	m_rotdir       = m_Seeds_to_DTI.SubMatrix(1,3,1,3);	
      }
      else{
	m_IsNonlinXfm = true;
	FnirtFileReader ffr(opts.seeds_to_dti.value());
	m_Seeds_to_DTI_warp = ffr.FieldAsNewimageVolume4D(true);
	if(opts.dti_to_seeds.value()==""){
	  cerr << "TRACT::Streamliner:: DTI -> Seeds transform needed" << endl;
	  exit(1);
	}
	FnirtFileReader iffr(opts.dti_to_seeds.value());
	m_DTI_to_Seeds_warp = iffr.FieldAsNewimageVolume4D(true);

	// now calculate the jacobian of this transformation (useful for rotating vectors)
	imgradient(m_Seeds_to_DTI_warp[0],m_jacx);
	imgradient(m_Seeds_to_DTI_warp[1],m_jacy);
	imgradient(m_Seeds_to_DTI_warp[2],m_jacz);
      }
    }
    
    vols.initialise(opts.basename.value(),m_mask);
    m_path.reserve(opts.nsteps.value());
    m_diff_path.reserve(opts.nsteps.value());
    if(m_surfexists)
      m_crossedvox.reserve(opts.nsteps.value());

    m_x_s_init=0;
    m_y_s_init=0;
    m_z_s_init=0;

    m_inmask3.reserve(opts.nsteps.value());
    m_inlrmask3.reserve(opts.nsteps.value());
    m_inmask3_aux.reserve(opts.nsteps.value());
    m_inlrmask3_aux.reserve(opts.nsteps.value()); 	
  }
  
  void Streamliner::rotdir(const ColumnVector& dir,ColumnVector& rotdir,
			   const float& x,const float& y,const float& z){
    if(!m_IsNonlinXfm){
      rotdir = m_rotdir*dir;
      if(rotdir.MaximumAbsoluteValue()>0)
	rotdir = rotdir/std::sqrt(rotdir.SumSquare());
    }
    else{
      ColumnVector xyz_dti(3),xyz_seeds(3);
      xyz_seeds << x << y << z;
      xyz_dti = NewimageCoord2NewimageCoord(m_DTI_to_Seeds_warp,
					    false,m_seeds.get_refvol(),m_mask,xyz_seeds);

      Matrix F(3,3),Jw(3,3);
      int x=(int)MISCMATHS::round((float)xyz_dti(1));
      int y=(int)MISCMATHS::round((float)xyz_dti(2));
      int z=(int)MISCMATHS::round((float)xyz_dti(3));
      Jw << m_jacx(x,y,z,0) << m_jacx(x,y,z,1) << m_jacx(x,y,z,2)
	 << m_jacy(x,y,z,0) << m_jacy(x,y,z,1) << m_jacy(x,y,z,2)
	 << m_jacz(x,y,z,0) << m_jacz(x,y,z,1) << m_jacz(x,y,z,2);
      F = (IdentityMatrix(3) + Jw).i();
      rotdir = F*dir;
      if(rotdir.MaximumAbsoluteValue()>0)
	rotdir = rotdir/std::sqrt(rotdir.SumSquare());

    }
  }
  
  int Streamliner::streamline(const float& x_init,const float& y_init,const float& z_init, 
			      const ColumnVector& dim_seeds,const int& fibst){ 
    //Tracer_Plus tr("Streamliner::streamline");
    //fibst tells tractvolsx which fibre to start with if there are more than one..
    //x_init etc. are in seed space...
    vols.reset(fibst);
    m_x_s_init=x_init; // seed x position in voxels
    m_y_s_init=y_init; // and y
    m_z_s_init=z_init; // and z
    ColumnVector xyz_seeds(3);
    xyz_seeds<<x_init<<y_init<<z_init;
    ColumnVector xyz_dti;
    ColumnVector th_ph_f;
    float xst,yst,zst,x,y,z,tmp2,cthr=opts.c_thr.value();
    float pref_x=0,pref_y=0,pref_z=0;
    int x_s,y_s,z_s;
    int x_p,y_p,z_p;
    int sampled_fib=fibst;

    // find xyz in dti space
    if(!m_IsNonlinXfm)
      xyz_dti = vox_to_vox(xyz_seeds,dim_seeds,vols.dimensions(),m_Seeds_to_DTI);    
    else{
      xyz_dti = NewimageCoord2NewimageCoord(m_DTI_to_Seeds_warp,false,m_seeds.get_refvol(),m_mask,xyz_seeds);
    }


    xst=xyz_dti(1);yst=xyz_dti(2);zst=xyz_dti(3);
    m_path.clear();
    m_diff_path.clear();
    x=xst;y=yst;z=zst;
    m_part.change_xyz(x,y,z);
    x_p=(int)MISCMATHS::round(m_part.x());
    y_p=(int)MISCMATHS::round(m_part.y());
    z_p=(int)MISCMATHS::round(m_part.z());
    //    x_p=(int)MISCMATHS::round(x);y_p=(int)MISCMATHS::round(y);z_p=(int)MISCMATHS::round(z);

    float          pathlength=0;
    bool           rubbish_passed=false;
    bool           wayorder=true;
    vector<int>    waycrossed;
    
    bool            wtstop=false;
    //Counter that determines propagation using wtstop masks. 
    //In case of volume wt_stop masks, this means that a streamline can enter the wt volume, propagate within it, but it will stop upon exiting it.
    //In case of surface wt_stop masks, a streamline can cross the surface, but propagation will stop when it crosses the surface again.
    //It treats each volume and each surface independently
    // wtstop_flags: 
    // 0-> seed is inside the roi, let it go out
    // 1-> not inside roi yet
    // 2-> inside the roi
    // 3-> going outside the roi ... STOP

    vector<int>    wtstop_flags;  //temporary structures for wtstop propagation
    vector<int>    wtstopcrossed;

    vector<int>    crossedlocs3;
    vector< pair<int,int> > surf_Triangle3;
    vector<ColumnVector> crossedvox;

    int            cnt=-1;
    // loopcheck stuff
    float oldrx,oldry,oldrz;
    int   lcx,lcy,lcz;
    bool forcedir=false;
    //NB - this only goes in one direction!!

    vector<int> m_way_passed_flags_updated;	// it keeps a record to know what positions of m_way_passed_flags have been updated
    for(unsigned int i=0; i<m_way_passed_flags.size();i++) // to undo them in case rubbish_passed=1
	m_way_passed_flags_updated.push_back(0);
    if(opts.onewaycondition.value()){
      for(unsigned int i=0; i<m_way_passed_flags.size();i++)
	m_way_passed_flags[i]=0;
      if(opts.network.value()){
	m_net_passed_flags=0;
      }
    }
    
    if(m_surfexists){m_crossedvox.clear();}

    if (opts.wtstopfiles.value()!=""){ //use wtstop masks
      for (int i=0; i<m_wtstopmasks.nRois(); i++)
	wtstop_flags.push_back(1);	// not inside roi yet

      wtstopcrossed.clear();
      // if seed is inside roi, then set to 0 (only for volumes, not for surfaces)
      m_wtstopmasks.has_crossed_roi_vols(xyz_seeds,wtstopcrossed);
      for(unsigned int wm=0;wm<wtstopcrossed.size();wm++) 
          wtstop_flags[wtstopcrossed[wm]]=0;  // seed is inside the roi, let it go out                              
    }

    for(int it=1;it<=opts.nsteps.value()/2;it++){
      
      if((m_mask(x_p,y_p,z_p)!=0)){
	///////////////////////////////////
	//loopchecking
	///////////////////////////////////
	if(opts.loopcheck.value()){
	  lcx=(int)MISCMATHS::round(m_part.x()/m_lcrat);
	  lcy=(int)MISCMATHS::round(m_part.y()/m_lcrat);
	  lcz=(int)MISCMATHS::round(m_part.z()/m_lcrat);
	  oldrx=m_loopcheck(lcx,lcy,lcz,0);
	  oldry=m_loopcheck(lcx,lcy,lcz,1);
	  oldrz=m_loopcheck(lcx,lcy,lcz,2);
	  if(m_part.rx()*oldrx+m_part.ry()*oldry+m_part.rz()*oldrz<0){break;}
	  m_loopcheck(lcx,lcy,lcz,0)=m_part.rx();
	  m_loopcheck(lcx,lcy,lcz,1)=m_part.ry();
	  m_loopcheck(lcx,lcy,lcz,2)=m_part.rz();
	}
	
	
	x=m_part.x();y=m_part.y();z=m_part.z();
	xyz_dti <<x<<y<<z;
	x_p=(int)MISCMATHS::round(x);y_p=(int)MISCMATHS::round(y);z_p=(int)MISCMATHS::round(z);

	// now find xyz in seeds space
	if(cnt>=0){
	  if(!m_IsNonlinXfm)
	    xyz_seeds = vox_to_vox(xyz_dti,vols.dimensions(),dim_seeds,m_DTI_to_Seeds);    
	  else{
	    xyz_seeds = NewimageCoord2NewimageCoord(m_Seeds_to_DTI_warp,false,m_mask,m_seeds.get_refvol(),xyz_dti);
	    }	  
	}
	
	x_s =(int)MISCMATHS::round((float)xyz_seeds(1));
	y_s =(int)MISCMATHS::round((float)xyz_seeds(2));
	z_s =(int)MISCMATHS::round((float)xyz_seeds(3));
	

	// how prefdir works:
	// if a vector field is given (tested using 4th dim==3) then set prefdir to that orientation
	// eventually will be flipped to match the current direction of travel
	// if dim4<3, then :
	//    - use dim1 scalar as fibre sampling method (i.e. randfib 1 2 3 )
	//    - use dim2 scalar (when given) as which fibre to follow blindly, i.e. without flipping
	pref_x=0;pref_y=0;pref_z=0;
	int sample_fib = 0;
	if(opts.prefdirfile.value()!=""){
	  if(m_prefdir.tsize()==3){
	    pref_x = m_prefdir(x_p,y_p,z_p,0);
	    pref_y = m_prefdir(x_p,y_p,z_p,1);
	    pref_z = m_prefdir(x_p,y_p,z_p,2);
	  }
	  else{
	    if(m_prefdir(x_s,y_s,z_s,0)!=0)
	      sample_fib = (int)m_prefdir(x_p,y_p,z_p,0);
	    if(sample_fib>3)sample_fib=3;	    
	  }
	}

	// // attractor mask 
	// if(opts.pullmask.value()!=""){
	//   ColumnVector pulldir(3);float mindist;
	//   if(m_pullmask.is_near_surface(xyz_seeds,opts.pulldist.value(),pulldir,mindist)){
	//     m_prefx = pulldir(1);
	//     m_prefy = pulldir(2);
	//     m_prefz = pulldir(3);
	//   }	  
	// }


	// augment the path	
	m_path.push_back(xyz_seeds);cnt++;
	if(cnt>0)
	  pathlength += opts.steplength.value();

	// Do this only once as it is the same for all surfaces
	if(cnt>0 && m_surfexists){
	  float line[2][3]={{m_path[cnt-1](1),m_path[cnt-1](2),m_path[cnt-1](3)},
			    {m_path[cnt](1),m_path[cnt](2),m_path[cnt](3)}};
	  m_seeds.line_crossed_voxels(line,crossedvox); //crossed in voxels
	  m_crossedvox.push_back(crossedvox);	      
	}
		
	// only test exclusion after at least one step
	if(opts.rubbishfile.value()!="" && cnt>0){
	  if(m_rubbish.has_crossed(m_path[cnt-1],m_path[cnt],crossedvox)){
	    rubbish_passed=1;
	    for(unsigned int i=0; i<m_way_passed_flags_updated.size();i++){ // undo updated in this part
		if(m_way_passed_flags_updated[i])
			m_way_passed_flags[i]=0;
	    }
	    break;
	  }
	}

	// update every passed_flag
	if(m_way_passed_flags.size()>0){
	  if(cnt>0){
	    waycrossed.clear();
	    m_waymasks.has_crossed_roi(m_path[cnt-1],m_path[cnt],crossedvox,waycrossed);

	    for(unsigned int wm=0;wm<waycrossed.size();wm++){
              if(m_way_passed_flags[waycrossed[wm]]==0){ 	// not crossed yet
		 m_way_passed_flags_updated[waycrossed[wm]]=1;
		 m_way_passed_flags[waycrossed[wm]]=1;
	      }
	    }
	    // check if order is respected
	    if(opts.wayorder.value() && wayorder){
	      for(unsigned int wm=0;wm<waycrossed.size();wm++){
		for(int wm2=0;wm2<waycrossed[wm];wm2++){
		  if(m_way_passed_flags[wm2]==0){ wayorder=false; break;}
		}
		if(!wayorder)break;
	      }
	    }	    
	  }
	}
	
	if(opts.network.value()){
	  if(cnt>0){
	    waycrossed.clear();
	    m_netmasks.has_crossed_roi(m_path[cnt-1],m_path[cnt],crossedvox,waycrossed);

	    for(unsigned int wm=0;wm<waycrossed.size();wm++){
	      m_net_passed_flags(waycrossed[wm]+1)=1;
	      m_net_passed(waycrossed[wm]+1)=1;
	    }
	  }	  
	}

	// Test WT-stopping after at least one step
	if(opts.wtstopfiles.value()!="" && cnt>0){
	  wtstopcrossed.clear();
	  m_wtstopmasks.has_crossed_roi(m_path[cnt-1],m_path[cnt],crossedvox,wtstopcrossed);

	  for (unsigned int r=0; r<wtstop_flags.size(); r++){
 	    bool crossed=false;
	    for(unsigned int wm=0;wm<wtstopcrossed.size();wm++){
              if(r==(unsigned int)wtstopcrossed[wm]) crossed=true;
	    }

	    if ((wtstop_flags[r]==0) && (!crossed)) wtstop_flags[r]=1;		// going outside the roi, but we started inside
	    else if ((wtstop_flags[r]==1) && (crossed)) wtstop_flags[r]=2; 	// going into the roi
	    else if ((wtstop_flags[r]==2) && (!crossed) & (m_wtstopmasks.get_roitype(r)==VOLUME)){ 	// going outside the roi, stop, works for volumes
		wtstop=true;	
	    }
	    else if ((wtstop_flags[r]==2) && (crossed) & (m_wtstopmasks.get_roitype(r)==SURFACE)){ 	// second time that a surface is crossed, stop
		wtstop=true;	
	    }
          }
	  if (wtstop){ //stop
	    // remove last path point? Yes, so that counters don't count 2nd crossings	
	    m_path.pop_back();cnt--;
	    if(cnt>0) pathlength -= opts.steplength.value();
	    break;
	  }
	}

	// //////////////////////////////
	// update locations for matrix3
	if(opts.matrix3out.value() && cnt>0){
	  waycrossed.clear();crossedlocs3.clear();surf_Triangle3.clear();
	  if(m_mask3.has_crossed_roi(m_path[cnt-1],m_path[cnt],crossedvox,waycrossed,crossedlocs3,surf_Triangle3,opts.closestvertex.value())){	    
	    fill_inmask3(crossedlocs3,surf_Triangle3,pathlength);
	  }
	  if(opts.lrmask3.value()!=""){
	    waycrossed.clear();crossedlocs3.clear();surf_Triangle3.clear();
	    if(m_lrmask3.has_crossed_roi(m_path[cnt-1],m_path[cnt],crossedvox,waycrossed,crossedlocs3,surf_Triangle3,opts.closestvertex.value())){	    
	      fill_inlrmask3(crossedlocs3,surf_Triangle3,pathlength);
	    }	 
	  }
	}


	// //////////////////////////////

	// sample a new fibre orientation
	int newx,newy,newz;	
	if(opts.skipmask.value() == ""){	  
	  th_ph_f = vols.sample(m_part.x(),m_part.y(),m_part.z(),    // sample at this location
				m_part.rx(),m_part.ry(),m_part.rz(), // choose closest sample to this
				pref_x,pref_y,pref_z,                // unless we have this prefered direction 
				sample_fib,                          // controls probabilities of sampling a given fibre
				sampled_fib,                         // which fibre has been sampled
				newx,newy,newz);                     // which voxel has been considered
	}
	else{
	  if(m_skipmask(x_s,y_s,z_s)==0)//only sample outside skipmask
	    th_ph_f=vols.sample(m_part.x(),m_part.y(),m_part.z(),
				m_part.rx(),m_part.ry(),m_part.rz(),
				pref_x,pref_y,pref_z,
				sample_fib,
				sampled_fib,
				newx,newy,newz);
	}
	
	ColumnVector voxfib(4);
	voxfib<<newx<<newy<<newz<<sampled_fib;
	m_diff_path.push_back(voxfib);
	

	// only test stopping after at least one step
	if(opts.stopfile.value()!="" && m_path.size()>1){
	  if(m_path.size()==2 && opts.forcefirststep.value()){
	    // do nothing
	  }
	  else if(m_stop.has_crossed(m_path[cnt-1],m_path[cnt],crossedvox)){
	    break;	    
	  }	  
	}

	// jump
	tmp2=0;
	if(opts.usef.value()){tmp2=(float)rand()/(float)RAND_MAX;}

	if(th_ph_f(3)>tmp2){ //volume fraction criterion
	  if(opts.loccurvthresh.value()!=""){
	    cthr=m_loccurvthresh(x_p,y_p,z_p);	    
	  }
	  if(!m_part.check_dir(th_ph_f(1),th_ph_f(2),cthr,forcedir)){ // curvature threshold
	    break;
	  }

	  if((th_ph_f(1)!=0 && th_ph_f(2)!=0)){
	    if( (m_mask(x_p,y_p,z_p) != 0) ){

	      if(!opts.modeuler.value())
		m_part.jump(th_ph_f(1),th_ph_f(2),forcedir);

	      else{
		  ColumnVector test_th_ph_f;		  
		  m_part.testjump(th_ph_f(1),th_ph_f(2),forcedir);

		  test_th_ph_f=vols.sample(m_part.testx(),m_part.testy(),m_part.testz(),
					   m_part.rx(),m_part.ry(),m_part.rz(),
					   pref_x,pref_y,pref_z,sample_fib,sampled_fib,newx,newy,newz);
		  test_th_ph_f=mean_sph_pol(th_ph_f,test_th_ph_f);
		  m_part.jump(test_th_ph_f(1),test_th_ph_f(2),forcedir);
		}
	    }
	    else{
	      break; // outside mask	    
	    }
	  }
	  else{
	    break; // direction=zero
	  }
	}
	else{
	  break; // volume fraction threshold
	}
	
      }
      else{
	break; // outside mask
      }
      
    } // Close Step Number Loop (done tracking sample)


    // reset loopcheck box
    if(opts.loopcheck.value()){
      m_loopcheck=0;
    }
    
    // rejflag = 0 (accept), 1 (reject) or 2 (wait for second direction)
    int rejflag=0;

    if(rubbish_passed) return(1);
    if(pathlength<opts.distthresh.value()) return(1);

    if(opts.network.value()){
      unsigned int numpassed=0;
      for(int i=1; i<=m_net_passed_flags.Nrows();i++){
	if(m_net_passed_flags(i))numpassed++;
      }
      if(numpassed==0)rejflag=1;
      else if((int)numpassed<m_net_passed_flags.Nrows()){
	if(m_waycond=="AND"){
	  rejflag=opts.onewaycondition.value()?1:2;
	}
	else{
	  rejflag=0;
	}	
      }
      else rejflag=0;
    }
    if(rejflag==1) return(rejflag);
    
    // Waypoints
    if(m_way_passed_flags.size()!=0){
      unsigned int numpassed=0;
      for(unsigned int i=0; i<m_way_passed_flags.size();i++){
	if(m_way_passed_flags[i]>0)numpassed++;
      }
      
      if(numpassed==0)rejflag=1;
      else if(numpassed<m_way_passed_flags.size()){
	if(m_waycond=="AND"){
	  if(opts.onewaycondition.value()){rejflag=1;}
	  else{
	    if(opts.wayorder.value()){rejflag=wayorder?2:1;}
	    else{rejflag=2;}
	  }
	}
	else{
	  rejflag=0;
	}
      }
      else{
	if(m_waycond=="OR" || !opts.wayorder.value())
	  rejflag=0;
	else
	  rejflag=wayorder?0:1;
      }
    }

    return rejflag;    
  }
  
  void Counter::initialise(){
   
    if(opts.simpleout.value()||opts.omeanpathlength.value()){
      initialise_path_dist();
    }
    if(opts.s2tout.value()){
      initialise_seedcounts();
    }
    if(opts.matrix1out.value()){
      initialise_matrix1();
    }
    if(opts.matrix2out.value()){
      initialise_matrix2();
    }
    if(opts.matrix3out.value()){
      initialise_matrix3();
    }
    if(opts.matrix4out.value()){
      initialise_matrix4();
    }
  }

  
  void Counter::initialise_seedcounts(){
    // now the CSV class does all the work
    m_targetmasks.reinitialize(m_stline.get_seeds().get_refvol());
    m_targetmasks.set_convention(opts.meshspace.value());
    m_targetmasks.load_rois(opts.targetfile.value());
    m_targetmasks.reset_values();

    if(m_targetmasks.nSurfs()>0){m_stline.surfexists();}

    // seeds are CSV-format
    if(!opts.simple.value()){
      m_s2t_count=m_stline.get_seeds();
      m_s2t_count.reset_values();
      for(int m=0;m<m_targetmasks.nRois();m++){
	m_s2t_count.add_map();
      }
      if(opts.omeanpathlength.value()){
	m_s2t_count2=m_stline.get_seeds();
        m_s2t_count2.reset_values();
        for(int m=0;m<m_targetmasks.nRois();m++){
	  m_s2t_count2.add_map();
        }
      }
    }
    // seeds are text
    if(opts.simple.value() || opts.s2tastext.value() || opts.omeanpathlength.value()){
      m_s2tastext.ReSize(m_numseeds,m_targetmasks.nRois());
      m_s2tastext=0;
      if(opts.omeanpathlength.value()){
	m_s2tastext2.ReSize(m_numseeds,m_targetmasks.nRois());
        m_s2tastext2=0;
      }
    }

    m_s2trow=1;
    m_targflags.resize(m_targetmasks.nRois());

    if(opts.targetpaths.value()){
      m_targetpaths.reinitialize(m_prob.xsize(),m_prob.ysize(),m_prob.zsize(),m_targetmasks.nRois());
      m_targetpaths=0;
    }

    // init separate tract counters
    if(opts.targetpaths.value()){
      m_prob_multi.resize(m_targetmasks.nRois());
      for(int t=0;t<m_targetmasks.nRois();t++){
	m_prob_multi[t].reinitialize(m_stline.get_seeds().xsize(),
				     m_stline.get_seeds().ysize(),
				     m_stline.get_seeds().zsize());
	copybasicproperties(m_stline.get_seeds().get_refvol(),m_prob_multi[t]);
	m_prob_multi[t]=0;  
      }
      if(opts.omeanpathlength.value()){
        m_prob_multi2.resize(m_targetmasks.nRois());
        for(int t=0;t<m_targetmasks.nRois();t++){
	  m_prob_multi2[t].reinitialize(m_stline.get_seeds().xsize(),
				       m_stline.get_seeds().ysize(),
				       m_stline.get_seeds().zsize());
	  copybasicproperties(m_stline.get_seeds().get_refvol(),m_prob_multi2[t]);
	  m_prob_multi2[t]=0;  
        }
      }
      if(opts.opathdir.set()){
	m_localdir_multi.resize(m_targetmasks.nRois());
	for(int t=0;t<m_targetmasks.nRois();t++){
	  m_localdir_multi[t].reinitialize(m_stline.get_seeds().xsize(),
					   m_stline.get_seeds().ysize(),
					   m_stline.get_seeds().zsize(),6);
	  copybasicproperties(m_stline.get_seeds().get_refvol(),m_localdir_multi[t]);
	  m_localdir_multi[t]=0;	
	}
      }
    }
  }
  

  // matrix1 is nseeds X nseeds
  void Counter::initialise_matrix1(){
    m_ConMat1 = new SpMat<float> (m_numseeds,m_numseeds);
    if(opts.omeanpathlength.value()) m_ConMat1b = new SpMat<float> (m_numseeds,m_numseeds);
    m_Conrow1 = 1;

    vector<ColumnVector> coords = m_stline.get_seeds().get_locs_coords();
    vector<int> roicind         = m_stline.get_seeds().get_locs_coord_index();
    vector<int> roiind          = m_stline.get_seeds().get_locs_roi_index();
    
    Matrix CoordMat1(m_numseeds,5); 
    for (unsigned int i=0;i<coords.size();i++)
      CoordMat1.Row(i+1) << (float)coords[i](1) 
			 << (float)coords[i](2)
			 << (float)coords[i](3)
			 << roiind[i]
			 << roicind[i];
    
    applycoordchange(CoordMat1, m_stline.get_seeds().get_refvol().niftivox2newimagevox_mat().i());
    write_ascii_matrix(CoordMat1,logger.appendDir("coords_for_fdt_matrix1"));
  }
  
  // matrix2 is nseeds X nlrmask
  void Counter::initialise_matrix2(){
    // init lrmask-related 
    read_volume(m_lrmask,opts.lrmask.value());
    m_beenhere2.reinitialize(m_lrmask.xsize(),m_lrmask.ysize(),m_lrmask.zsize());
    m_beenhere2=0;
    m_lookup2.reinitialize(m_lrmask.xsize(),m_lrmask.ysize(),m_lrmask.zsize());
    copybasicproperties(m_lrmask,m_lookup2);
    m_lookup2=0;
    m_lrdim.ReSize(3);
    m_lrdim<<m_lrmask.xdim()<<m_lrmask.ydim()<<m_lrmask.zdim();

    int numnz=0;    
    for(int Wz=m_lrmask.minz();Wz<=m_lrmask.maxz();Wz++)
      for(int Wy=m_lrmask.miny();Wy<=m_lrmask.maxy();Wy++)
	for(int Wx=m_lrmask.minx();Wx<=m_lrmask.maxx();Wx++)
	  if(m_lrmask.value(Wx,Wy,Wz)!=0){
	    numnz++;
	    m_lookup2(Wx,Wy,Wz)=numnz;
	  }
    
    Matrix CoordMat_tract2(numnz,3);
    int mytrow=1;
    for(int Wz=m_lrmask.minz();Wz<=m_lrmask.maxz();Wz++)
      for(int Wy=m_lrmask.miny();Wy<=m_lrmask.maxy();Wy++)
	for(int Wx=m_lrmask.minx();Wx<=m_lrmask.maxx();Wx++)
	  if(m_lrmask(Wx,Wy,Wz)!=0){
	    CoordMat_tract2(mytrow,1)=Wx;
	    CoordMat_tract2(mytrow,2)=Wy;
	    CoordMat_tract2(mytrow,3)=Wz;
	    mytrow++;
	  }

    applycoordchange(CoordMat_tract2, m_lrmask.niftivox2newimagevox_mat().i());
    write_ascii_matrix(CoordMat_tract2,logger.appendDir("tract_space_coords_for_fdt_matrix2"));

    save_volume(m_lookup2,logger.appendDir("lookup_tractspace_fdt_matrix2"));

    
    // init matrix2-related
    m_ConMat2 = new SpMat<float> (m_numseeds,numnz);
    if(opts.omeanpathlength.value()) m_ConMat2b = new SpMat<float> (m_numseeds,numnz);

    if( !opts.simple.value()){
      vector<ColumnVector> coords = m_stline.get_seeds().get_locs_coords();
      vector<int> roicind         = m_stline.get_seeds().get_locs_coord_index();
      vector<int> roiind          = m_stline.get_seeds().get_locs_roi_index();

      Matrix CoordMat2(m_numseeds,5);
      for (unsigned int i=0;i<coords.size();i++)
	CoordMat2.Row(i+1) << (float)coords[i](1) 
			   << (float)coords[i](2)
			   << (float)coords[i](3)
			   << roiind[i]
			   << roicind[i];
      
      applycoordchange(CoordMat2, m_stline.get_seeds().get_refvol().niftivox2newimagevox_mat().i());
      write_ascii_matrix(CoordMat2,logger.appendDir("coords_for_fdt_matrix2"));            

    }


  }
  
  // the following will use the voxels of a given mask to initialise 
  // an NxN matrix. This matrix will store the number of samples from 
  // each seed voxel that have made it to each pair of voxels in the mask
  // now mask3 is in CSV format
  void Counter::initialise_matrix3(){
    m_stline.init_mask3();

    int nmask3 = m_stline.get_mask3().nLocs();
    int nlrmask3;
    if(opts.lrmask3.value()!="")
      nlrmask3 = m_stline.get_lrmask3().nLocs();
    else
      nlrmask3 = nmask3;

    // recalculate nmask3 if lowres surface provided
    m_ConMat3  = new SpMat<float> (nmask3,nlrmask3); 
    if(opts.omeanpathlength.value()) m_ConMat3b = new SpMat<float> (nmask3,nlrmask3); 
    //OUT(m_ConMat3->Nrows());
    //OUT(m_ConMat3->Ncols());
    //exit(1);

    // save lookup tables...
    
    CSV  mask3(m_stline.get_mask3());
    vector<ColumnVector> coords = mask3.get_locs_coords();
    vector<int> roicind = mask3.get_locs_coord_index();
    vector<int> roiind = mask3.get_locs_roi_index();
    
    Matrix mat(mask3.nLocs(),5);
    for (unsigned int i=0;i<coords.size();i++)
      mat.Row(i+1) << coords[i](1) 
		   << coords[i](2)
		   << coords[i](3)
		   << roiind[i]
		   << roicind[i];
    
    applycoordchange(mat, m_stline.get_seeds().get_refvol().niftivox2newimagevox_mat().i());

    write_ascii_matrix(mat,logger.appendDir("coords_for_fdt_matrix3"));

    if(opts.lrmask3.value()!=""){
      CSV lrmask3(m_stline.get_lrmask3());
      vector<ColumnVector> lrcoords = lrmask3.get_locs_coords();
      vector<int> lrroicind = lrmask3.get_locs_coord_index();
      vector<int> lrroiind = lrmask3.get_locs_roi_index();
    
      Matrix lrmat(lrmask3.nLocs(),5);
      for (unsigned int i=0;i<lrcoords.size();i++)
	lrmat.Row(i+1) << lrcoords[i](1) 
		     << lrcoords[i](2)
		     << lrcoords[i](3)
		     << lrroiind[i]
		     << lrroicind[i];
    
      applycoordchange(lrmat, m_stline.get_seeds().get_refvol().niftivox2newimagevox_mat().i());
      MISCMATHS::write_ascii_matrix(lrmat,logger.appendDir("tract_space_coords_for_fdt_matrix3")); 
    }

  }


 
  // matrix4 is ndtimask X nseeds
  // unless mask4 is set
  void Counter::initialise_matrix4(){
    if(opts.simple.value()){
      cerr<<"Matrix4 output not compatible with --simple mode"<<endl;
      exit(1);
    }

    // columns are brain mask in diffusion space
    read_volume(m_dtimask,opts.dtimask.value());
    m_beenhere4.reinitialize(m_dtimask.xsize(),m_dtimask.ysize(),m_dtimask.zsize());
    m_beenhere4=0;
    m_lookup4.reinitialize(m_dtimask.xsize(),m_dtimask.ysize(),m_dtimask.zsize());
    copybasicproperties(m_dtimask,m_lookup4);
    m_lookup4=0;
    m_dtidim.ReSize(3);
    m_dtidim<<m_dtimask.xdim()<<m_dtimask.ydim()<<m_dtimask.zdim();
    int numnz=0;    
    for(int Wz=m_dtimask.minz();Wz<=m_dtimask.maxz();Wz++)
      for(int Wy=m_dtimask.miny();Wy<=m_dtimask.maxy();Wy++)
	for(int Wx=m_dtimask.minx();Wx<=m_dtimask.maxx();Wx++)
	  if(m_dtimask.value(Wx,Wy,Wz)!=0){
	    numnz++;
	    m_lookup4(Wx,Wy,Wz)=numnz;
	  }
    Matrix CoordMat_tract4(numnz,3);
    int mytrow=1;
    for(int Wz=m_dtimask.minz();Wz<=m_dtimask.maxz();Wz++)
      for(int Wy=m_dtimask.miny();Wy<=m_dtimask.maxy();Wy++)
	for(int Wx=m_dtimask.minx();Wx<=m_dtimask.maxx();Wx++)
	  if(m_dtimask(Wx,Wy,Wz)!=0){
	    CoordMat_tract4(mytrow,1)=Wx;
	    CoordMat_tract4(mytrow,2)=Wy;
	    CoordMat_tract4(mytrow,3)=Wz;
	    mytrow++;
	  }

    applycoordchange(CoordMat_tract4, m_dtimask.niftivox2newimagevox_mat().i());
    write_ascii_matrix(CoordMat_tract4,logger.appendDir("tract_space_coords_for_fdt_matrix4"));
    save_volume(m_lookup4,logger.appendDir("lookup_tractspace_fdt_matrix4"));

    
    // init matrix4-related
    if(opts.mask4.value()==""){
      m_ConMat4 = new SpMat_HCP(numnz,m_numseeds);
    }
    else{
      m_mask4.reinitialize(m_stline.get_seeds().get_refvol());
      m_mask4.set_convention(opts.meshspace.value());
      m_mask4.load_rois(opts.mask4.value());
      if(m_mask4.nSurfs()>0){m_stline.surfexists();}

      m_ConMat4 = new SpMat_HCP(numnz,m_mask4.nLocs());
    }

    // COMMENTED OUT SOME OUTPUTS [NOT BEING USED AT THIS STAGE]
    // vector<ColumnVector> coords = m_stline.get_seeds().get_locs_coords();
    // vector<int> roicind         = m_stline.get_seeds().get_locs_coord_index();
    // vector<int> roiind          = m_stline.get_seeds().get_locs_roi_index();
    // Matrix CoordMat4(m_numseeds,5);
    // for (unsigned int i=0;i<coords.size();i++)
    //   CoordMat4.Row(i+1) << (float)coords[i](1) 
    // 			 << (float)coords[i](2)
    // 			 << (float)coords[i](3)
    // 			 << roiind[i]
    // 			 << roicind[i];
    // applycoordchange(CoordMat4, m_stline.get_seeds().get_refvol().niftivox2newimagevox_mat().i());
    // write_ascii_matrix(CoordMat4,logger.appendDir("coords_for_fdt_matrix4"));                  
    
  }
  

  void Counter::count_streamline(){
    if(opts.save_paths.value()){
      add_path();
    }
    if(opts.simpleout.value()||opts.matrix1out.value()||opts.omeanpathlength.value()){
      update_pathdist();
    }
    if(opts.s2tout.value()){
      update_seedcounts();
    }
    if(opts.matrix2out.value()){
      update_matrix2_row();
    }
    if(opts.matrix4out.value()){
      update_matrix4_col();
    }
  }

  void Counter::count_seed(){
    if(opts.s2tout.value()){
      m_s2trow++;
    }
    if(opts.matrix1out.value()){
      m_Conrow1++;;
    }
  }

  void Counter::clear_streamline(){
    if(opts.simpleout.value()||opts.matrix1out.value()||opts.omeanpathlength.value()){
      reset_beenhere();
    }
    if(opts.s2tout.value()){
      reset_targetflags();
    }
    if(opts.matrix2out.value()){
      reset_beenhere2();
    }
    if(opts.matrix3out.value()){
      reset_beenhere3();
    }
    if(opts.matrix4out.value()){
      reset_beenhere4();
    }
  }

  void Counter::update_pathdist(){
    //Tracer_Plus tr("Counter::update_pathdist");
    if(m_path.size()<2){return;}
    int x_s,y_s,z_s;
    float pathlength=0;
    vector<int> crossedrois,crossedlocs;
    vector< pair<int,int> > surf_Triangle;
    vector<ColumnVector> crossedvox;
    int nlocs=0;
    int offset=-1;
    bool restarted=false;
    for(unsigned int i=0;i<m_path.size();i++){
      x_s=(int)MISCMATHS::round((float)m_path[i](1));
      y_s=(int)MISCMATHS::round((float)m_path[i](2));
      z_s=(int)MISCMATHS::round((float)m_path[i](3));
      // check here if back to seed
      if(i>0 && (m_path[i]-m_path[0]).MaximumAbsoluteValue()==0){
	//m_lastpoint(x_s,y_s,z_s)+=1;  
	pathlength=0;
	offset-=1;
	restarted=true;
      }
      if(m_beenhere(x_s,y_s,z_s)==0){
	if(!opts.pathdist.value())
	  m_prob(x_s,y_s,z_s)+=1; 
	else
	  m_prob(x_s,y_s,z_s)+=pathlength;
        if(opts.omeanpathlength.value()){
	  if(!opts.pathdist.value())
	    m_prob2(x_s,y_s,z_s)+=pathlength; 
	  else
	    m_prob2(x_s,y_s,z_s)+=1;
        }
	m_beenhere(x_s,y_s,z_s)=1;
	
	if(opts.opathdir.value() && !restarted){
	  ColumnVector v(3);
	  if(i==0)
	    v=m_path[i+1]-m_path[i];	  
	  else
	    v=m_path[i]-m_path[i-1];	  
	  float ss=v.SumSquare();
	  if(ss>0){
	    v/=std::sqrt(ss);
	  }

	  // Add direction (needs to account for the current direction and flip if necessary)
	  
	  // lower diagonal rows (because of the way SymmetricMatrix works)
	  m_localdir(x_s,y_s,z_s,0)+=v(1)*v(1);
	  m_localdir(x_s,y_s,z_s,1)+=v(1)*v(2);
	  m_localdir(x_s,y_s,z_s,2)+=v(2)*v(2);
	  m_localdir(x_s,y_s,z_s,3)+=v(1)*v(3);
	  m_localdir(x_s,y_s,z_s,4)+=v(2)*v(3);
	  m_localdir(x_s,y_s,z_s,5)+=v(3)*v(3);
	  
	}
      }
      
      // Fill alternative mask
      // This mask's values are:
      //  0: location not to be considered
      //  1: location has not been visited 
      //  2: location has been visited

      if(opts.pathfile.set()){
	if(pathlength>0 && !restarted){
	  if(m_beenhere_alt.nSurfs()>0){
	    crossedvox=m_crossedvox[i+offset];
	  }

	  crossedrois.clear();crossedlocs.clear();
	  if(m_beenhere_alt.has_crossed_roi(m_path[i-1],m_path[i],crossedvox,crossedrois,crossedlocs,surf_Triangle,opts.closestvertex.value())){
	    nlocs+=crossedlocs.size();
	    for(unsigned int i=0;i<crossedlocs.size();i++){
	      if(m_beenhere_alt.get_value(crossedlocs[i])==1){	      
		if(!opts.pathdist.value())
		  m_prob_alt.add_value(crossedlocs[i],1); 
		else
		  m_prob_alt.add_value(crossedlocs[i],pathlength);
                if(opts.omeanpathlength.value()){
		  if(!opts.pathdist.value())
		    m_prob_alt2.add_value(crossedlocs[i],pathlength); 
		  else
	            m_prob_alt2.add_value(crossedlocs[i],1);
                }
		m_beenhere_alt.set_value(crossedlocs[i],2);
	      }
	    }
	  }
	}
      }
      pathlength+=opts.steplength.value();
      restarted=false;
      
    }

    // // Fill last point
    // int i=m_path.size()-1;
    // x_s=(int)MISCMATHS::round((float)m_path[i](1));
    // y_s=(int)MISCMATHS::round((float)m_path[i](2));
    // z_s=(int)MISCMATHS::round((float)m_path[i](3));
    // m_lastpoint(x_s,y_s,z_s)+=1;  

    // In network mode, update network matrix
    if(opts.network.value()){
      m_stline.update_mat();      
    }

  }

  void Counter::reset_beenhere(){
    int x_s,y_s,z_s,offset=-1;
    vector<int> crossedlocs,crossedrois;
    vector< pair<int,int> > surf_Triangle;
    vector<ColumnVector> crossedvox;
    bool hascrossed=false;
    for(unsigned int i=0;i<m_path.size();i++){
      x_s=(int)MISCMATHS::round((float)m_path[i](1));
      y_s=(int)MISCMATHS::round((float)m_path[i](2));
      z_s=(int)MISCMATHS::round((float)m_path[i](3));
      m_beenhere(x_s,y_s,z_s)=0;

      // back to first point? keep going
      if(i>0 && (m_path[i]-m_path[0]).MaximumAbsoluteValue()==0){
	offset-=1;
	continue;
      }

      // alternative roi
      if(opts.pathfile.set()){
	if(i>0){
	  if(m_beenhere_alt.nSurfs()>0){
	    crossedvox=m_crossedvox[i+offset];
	  }
	  if(m_beenhere_alt.has_crossed_roi(m_path[i-1],m_path[i],crossedvox,crossedrois,crossedlocs,surf_Triangle,opts.closestvertex.value()))
	    hascrossed=true;
	}
      }
    }
    if(opts.pathfile.set() && hascrossed){
      //cout<<"reset"<<endl;
      //OUT(crossedlocs.size());
      for(unsigned int i=0;i<crossedlocs.size();i++){
	if(m_beenhere_alt.get_value(crossedlocs[i])>1)
	  m_beenhere_alt.set_value(crossedlocs[i],1);
      }      
    }

  }

  void Counter::add_path(){
    m_save_paths.push_back(m_path);
  }
  void Counter::save_paths(){
    string filename=logger.appendDir("saved_paths.txt");
    ofstream of(filename.c_str());
    if (of.is_open()){
      for(unsigned int i=0;i<m_save_paths.size();i++){
	stringstream flot;
	flot << "# " << m_save_paths[i].size()<<endl;
	for(unsigned int j=0;j<m_save_paths[i].size();j++){
	  flot << m_save_paths[i][j](1) << " " << m_save_paths[i][j](2) << " " << m_save_paths[i][j](3) << endl;
	}
	of<<flot.str();
      }
      of.close();
    }
    else{
      cerr<<"Counter::save_paths:error opening file for writing: "<<filename<<endl;
    }
  }
  

  void Counter::update_seedcounts(){
    //Tracer_Plus tr("Counter::update_seedcounts");
    vector<int> crossed;
    vector<ColumnVector> crossedvox;
    float       pathlength=0;
    int         cnt=0,offset=-1;

    for(unsigned int i=1;i<m_path.size();i++){
      pathlength+=opts.steplength.value();
      // check here if back to seed
      if((m_path[i]-m_path[0]).MaximumAbsoluteValue()==0){
	offset-=1;
	pathlength=0;
	continue;
      }

      // this would be more efficient if we only checked 
      // masks that haven't been crossed yet...
      crossed.clear();
      if(m_targetmasks.nSurfs()>0){
	// below gives a warning because addition of an unsigned int and a signed int is always an unsigned quantity
	if( (( i+offset ) < 0 ) || ((i+offset)>=m_crossedvox.size())){
	  cout<<"-----------------------"<<endl;
	  OUT(m_path.size());
	  OUT(m_crossedvox.size());
	  OUT(i);
	  OUT(offset);
	  for(unsigned int ii=0;ii<m_path.size();ii++)
	    cout<<m_path[ii](1)<<" "<<m_path[ii](2)<<" "<<m_path[ii](3)<<endl;
	  cout<<"-----------------------"<<endl;
	  exit(1);
	}
	crossedvox=m_crossedvox[i+offset];//stepcnt++;
      }

      m_targetmasks.has_crossed_roi(m_path[i-1],m_path[i],crossedvox,crossed);
      for(unsigned int t=0;t<crossed.size();t++){
	if(m_targflags[crossed[t]])continue;

	if(!opts.simple.value()){
	  if(!opts.pathdist.value())
	    m_s2t_count.add_map_value(m_curloc.loc,1,crossed[t]);
	  else
	    m_s2t_count.add_map_value(m_curloc.loc,pathlength,crossed[t]);
	  if(opts.omeanpathlength.value()){
	    if(!opts.pathdist.value())
	      m_s2t_count2.add_map_value(m_curloc.loc,pathlength,crossed[t]);
	    else
	      m_s2t_count2.add_map_value(m_curloc.loc,1,crossed[t]);
	  }
	}
	if(opts.simple.value() || opts.s2tastext.value()){
	  if(!opts.pathdist.value())
	    m_s2tastext(m_s2trow,crossed[t]+1)+=1;
	  else
	    m_s2tastext(m_s2trow,crossed[t]+1)+=pathlength;
          if(opts.omeanpathlength.value()){
	    if(!opts.pathdist.value())
	      m_s2tastext2(m_s2trow,crossed[t]+1)+=pathlength;
	    else
	      m_s2tastext2(m_s2trow,crossed[t]+1)+=1;
          }
	}
	m_targflags[crossed[t]]=true;cnt++;

      }
      // stop if they all have been crossed
      if(cnt==m_targetmasks.nRois())break;
      
    }
    // save separate paths
      if(opts.targetpaths.value()){
	// WHAT TO DO HERE???
	// 1st Loop = streamlines (like in m_prob update)
	// beenhere should be used
	// when updating, loop over crossed targets
	reset_beenhere();
	update_pathdist_multi();
      }
  }

  void Counter::update_pathdist_multi(){
    if(m_path.size()<1){return;}
    int x_s,y_s,z_s;
    float pathlength=0;
    for(unsigned int i=0;i<m_path.size();i++){
      pathlength+=opts.steplength.value();
      x_s=(int)MISCMATHS::round((float)m_path[i](1));
      y_s=(int)MISCMATHS::round((float)m_path[i](2));
      z_s=(int)MISCMATHS::round((float)m_path[i](3));
      // check here if back to seed
      if(i>0 && (m_path[i]-m_path[0]).MaximumAbsoluteValue()==0){
	pathlength=0;     
      }
      if(m_beenhere(x_s,y_s,z_s)==0){
	for(unsigned int t=0;t<m_targflags.size();t++){
	  if(!m_targflags[t]){continue;}
	  if(!opts.pathdist.value())
	    m_prob_multi[t](x_s,y_s,z_s)+=1; 
	  else
	    m_prob_multi[t](x_s,y_s,z_s)+=pathlength;
	  if(opts.omeanpathlength.value()){
	    if(!opts.pathdist.value())
	      m_prob_multi2[t](x_s,y_s,z_s)+=pathlength; 
	    else
	      m_prob_multi2[t](x_s,y_s,z_s)+=1;
 	  }
	  if(opts.opathdir.value() && i>0){
	    ColumnVector v(3);
	    v=m_path[i]-m_path[i-1];
	    v/=std::sqrt(v.SumSquare());
	    m_localdir_multi[t](x_s,y_s,z_s,0)+=v(1)*v(1);
	    m_localdir_multi[t](x_s,y_s,z_s,1)+=v(1)*v(2);
	    m_localdir_multi[t](x_s,y_s,z_s,2)+=v(2)*v(2);
	    m_localdir_multi[t](x_s,y_s,z_s,3)+=v(1)*v(3);
	    m_localdir_multi[t](x_s,y_s,z_s,4)+=v(2)*v(3);
	    m_localdir_multi[t](x_s,y_s,z_s,5)+=v(3)*v(3);
	  }	  	
	}
	m_beenhere(x_s,y_s,z_s)=1;
	
	
      }
    }
  }

  void Counter::update_matrix1(){
    //Tracer_Plus tr("Counter::update_matrix1");
    // use path and has_crossed
    float pathlength=opts.steplength.value(),val=1;
    vector<int> crossedseeds,crossedlocs;
    vector< pair<int,int> > surf_Triangle;
    vector< pair<int,infoVertex> > locs; //locs,values,rois and triangles
    vector<ColumnVector> crossedvox;
    int offset=-1;
    for(unsigned int i=1;i<m_path.size();i++){
      // check here if back to seed (not a robust way to do it...)
      if((m_path[i]-m_path[0]).MaximumAbsoluteValue()==0){
	pathlength=opts.steplength.value();
	val = opts.pathdist.value()?pathlength:1;
	if(opts.omeanpathlength.value()) val = pathlength;
	offset-=1;
	continue;
      }
      crossedseeds.clear();crossedlocs.clear();surf_Triangle.clear();
      if(m_stline.get_seeds().nSurfs()>0){
	crossedvox=m_crossedvox[i+offset];       
      }

      if(m_stline.get_seeds().has_crossed_roi(m_path[i-1],m_path[i],crossedvox,crossedseeds,crossedlocs,surf_Triangle,opts.closestvertex.value())){
	for(unsigned int j=0;j<crossedlocs.size();j++){
	  pair<int,infoVertex> mypair;
	  mypair.first=crossedlocs[j];
	  mypair.second.value=val;
	  mypair.second.mesh=surf_Triangle[j].first;
	  mypair.second.triangle=surf_Triangle[j].second;
	  locs.push_back(mypair);
	}
      }
      
      val = opts.pathdist.value()?pathlength:1;      
      if(opts.omeanpathlength.value()) val = pathlength;  
      pathlength+=opts.steplength.value();
    }
    // fill matrix1
    {
      //Tracer_Plus tr("make_unique");
      make_unique(locs);
    }

    vector< pair<int,int> > mytrianglesj; //list with the roi and triangles of an individual vertex
    for(unsigned int j=0;j<locs.size();j++){
      mytrianglesj.clear();
      pair<int,int> mypair;
      int index=locs[j].first;
      mypair.first=locs[j].second.mesh;
      mypair.second=locs[j].second.triangle;
      mytrianglesj.push_back(mypair);
      for(;(j+1)<locs.size() && locs[j+1].first==index;j++){  //same vertix - different roi-triangle
        pair<int,int> mypair2;
	mypair2.first=locs[j+1].second.mesh;
	mypair2.second=locs[j+1].second.triangle;
	mytrianglesj.push_back(mypair2); 
      }
      bool connect=false;
      if(m_curloc.loc!=locs[j].first){
        for(unsigned int ii=0;ii<m_curloc.triangles.size()&&!connect;ii++){
          for(unsigned int jj=0;jj<mytrianglesj.size()&&!connect;jj++){
            if(m_curloc.mesh!=mytrianglesj[jj].first || m_curloc.triangles[ii]!=mytrianglesj[jj].second ||
              m_curloc.triangles[ii]==-1 || mytrianglesj[jj].first==-1){
	      connect=true;
	    }
	  }
        }
      }
      if(connect){  
        m_ConMat1->AddTo(m_curloc.loc+1,locs[j].first+1,locs[j].second.value);
        if(opts.omeanpathlength.value()) m_ConMat1b->AddTo(m_curloc.loc+1,locs[j].first+1,1);
      }
    }
  }
  
  void Counter::update_matrix2_row(){
        //Tracer_Plus tr("Counter::update_matrix2");
    //run this one every streamline - not every voxel..
    float d=opts.steplength.value();
    int x_lr,y_lr,z_lr,Concol2;
    ColumnVector xyz_lr(3);
    for(unsigned int i=0;i<m_path.size();i++){
      // check here if back to seed
      if(i>0 && (m_path[i]-m_path[0]).MaximumAbsoluteValue()==0)
	d=opts.steplength.value();

      xyz_lr=vox_to_vox(m_path[i],m_seedsdim,m_lrdim,m_I);	
      x_lr=(int)MISCMATHS::round((float)xyz_lr(1));
      y_lr=(int)MISCMATHS::round((float)xyz_lr(2));
      z_lr=(int)MISCMATHS::round((float)xyz_lr(3));
      Concol2=m_lookup2(x_lr,y_lr,z_lr);

      if(Concol2>0){
	if(m_beenhere2(x_lr,y_lr,z_lr)==0){
	  if(opts.omeanpathlength.value()){
            m_ConMat2->AddTo(m_curloc.loc+1,Concol2,d);
	    m_ConMat2b->AddTo(m_curloc.loc+1,Concol2,1);
          }else{  
	    if(!opts.pathdist.value())
	      m_ConMat2->AddTo(m_curloc.loc+1,Concol2,1);
	    else
	      m_ConMat2->AddTo(m_curloc.loc+1,Concol2,d);
          }
	  m_beenhere2(x_lr,y_lr,z_lr)=1;
	  d+=opts.steplength.value();
	}
      }
    } 
  }

  void Counter::update_matrix3(){
    //Tracer_Plus tr("Counter::update_matrix3");
    vector< pair<int,infoVertex> >& inmask3 = m_stline.get_inmask3();
    vector< pair<int,infoVertex> >& inlrmask3 = m_stline.get_inlrmask3();
    bool uselr = (opts.lrmask3.value()!="");
    if(!uselr){
      if(inmask3.size()<2)return;
    }
    else{
      if(inmask3.size()<1 || inlrmask3.size()<1)return;
    }

    // remove duplicates
    make_unique(inmask3);

    vector< pair<int,int> > mytrianglesi; //list with the roi and triangles of an individual vertex
    vector< pair<int,int> > mytrianglesj; //list with the roi and triangles of an individual vertex
    if(!uselr){
      // where we update NxN matrix
      for(unsigned int i=0;i<inmask3.size();i++){
	mytrianglesi.clear();
	pair<int,int> mypair;
	int index=inmask3[i].first;
	mypair.first=inmask3[i].second.mesh;
	mypair.second=inmask3[i].second.triangle;
	mytrianglesi.push_back(mypair);
	for(;(i+1)<inmask3.size() && inmask3[i+1].first==index;i++){  //same vertix - different roi-triangle
	    pair<int,int> mypair;
	    mypair.first=inmask3[i+1].second.mesh;
	    mypair.second=inmask3[i+1].second.triangle;
	    mytrianglesi.push_back(mypair);
	}
	unsigned int j=i+1;
	for(;j<inmask3.size();j++){
	  mytrianglesj.clear();
	  pair<int,int> mypair;
          index=inmask3[j].first;
	  mypair.first=inmask3[j].second.mesh;
	  mypair.second=inmask3[j].second.triangle;
	  mytrianglesj.push_back(mypair);
	  for(;(j+1)<inmask3.size() && inmask3[j+1].first==index;j++){  //same vertix - different roi-triangle
	    pair<int,int> mypair;
	    mypair.first=inmask3[j+1].second.mesh;
	    mypair.second=inmask3[j+1].second.triangle;
	    mytrianglesj.push_back(mypair); 
	  }

	  bool connect=false;
	  for(unsigned int ii=0;ii<mytrianglesi.size()&&!connect;ii++){
	    for(unsigned int jj=0;jj<mytrianglesj.size()&&!connect;jj++){
              if(mytrianglesi[ii].first!=mytrianglesj[jj].first || mytrianglesi[ii].second!=mytrianglesj[jj].second ||
		mytrianglesi[ii].first==-1 || mytrianglesj[jj].first==-1){
		//if first is -1 is because is not a vertex, it is a voxel 
	        connect=true;
	      }
	    }
	  }
	  if(connect){
	    if(opts.omeanpathlength.value()){
	      float val = fabs(inmask3[i].second.value-inmask3[j].second.value);
	      m_ConMat3->AddTo(inmask3[i].first+1,inmask3[j].first+1,val);
              m_ConMat3b->AddTo(inmask3[i].first+1,inmask3[j].first+1,1); 
	    }else{
	      if(!opts.pathdist.value()){
	        m_ConMat3->AddTo(inmask3[i].first+1,inmask3[j].first+1,1);
	      }else{
	        float val = fabs(inmask3[i].second.value-inmask3[j].second.value);
	        m_ConMat3->AddTo(inmask3[i].first+1,inmask3[j].first+1,val);
	      }
	    }
          }
	}
      }
    }
    else{ // where we update Nxn matrix
      make_unique(inlrmask3);
      for(unsigned int i=0;i<inmask3.size();i++){
        mytrianglesi.clear();
	pair<int,int> mypair;
	int index=inmask3[i].first;
	mypair.first=inmask3[i].second.mesh;
	mypair.second=inmask3[i].second.triangle;
	mytrianglesi.push_back(mypair);
	for(;(i+1)<inmask3.size() && inmask3[i+1].first==index;i++){  //same vertix - different roi-triangle
	    pair<int,int> mypair;
	    mypair.first=inmask3[i+1].second.mesh;
	    mypair.second=inmask3[i+1].second.triangle;
	    mytrianglesi.push_back(mypair);
	}
	for(unsigned int j=0;j<inlrmask3.size();j++){
	  mytrianglesj.clear();
	  pair<int,int> mypair;
          index=inlrmask3[j].first;
	  mypair.first=inlrmask3[j].second.mesh;
	  mypair.second=inlrmask3[j].second.triangle;
	  mytrianglesj.push_back(mypair);
 	  for(;(j+1)<inlrmask3.size() && inlrmask3[j+1].first==index;j++){  //same vertix - different roi-triangle
	    pair<int,int> mypair;
	    mypair.first=inlrmask3[j+1].second.mesh;
	    mypair.second=inlrmask3[j+1].second.triangle;
	    mytrianglesj.push_back(mypair); 
	  }
	  bool connect=false;
	  for(unsigned int ii=0;ii<mytrianglesi.size()&&!connect;ii++){
	    for(unsigned int jj=0;jj<mytrianglesj.size()&&!connect;jj++){
              if(mytrianglesi[ii].first!=mytrianglesj[jj].first || mytrianglesi[ii].second!=mytrianglesj[jj].second ||
		mytrianglesi[ii].first==-1 || mytrianglesj[jj].first==-1){
	        connect=true;
	      }
	    }
	  }
	  if(connect){
	    if(opts.omeanpathlength.value()){
	      float val = fabs(inmask3[i].second.value-inlrmask3[j].second.value);
	      m_ConMat3->AddTo(inmask3[i].first+1,inlrmask3[j].first+1,val);
	      m_ConMat3b->AddTo(inmask3[i].first+1,inlrmask3[j].first+1,1);
	    }else{
	      if(!opts.pathdist.value())
	        m_ConMat3->AddTo(inmask3[i].first+1,inlrmask3[j].first+1,1);
	      else{
	        float val = fabs(inmask3[i].second.value-inlrmask3[j].second.value);
	        m_ConMat3->AddTo(inmask3[i].first+1,inlrmask3[j].first+1,val);
	      }
	    }
	  }
	}
      }
    }
  }  

  void Counter::reset_beenhere3(){
    m_stline.clear_inmask3();
    m_stline.clear_inlrmask3();
  }  
  
  void Counter::reset_beenhere2(){
    ColumnVector xyz_lr(3);
    for(unsigned int i=0;i<m_path.size();i++){
      xyz_lr=vox_to_vox(m_path[i],m_seedsdim,m_lrdim,m_I);
      m_beenhere2((int)MISCMATHS::round((float)xyz_lr(1)),
		  (int)MISCMATHS::round((float)xyz_lr(2)),
		  (int)MISCMATHS::round((float)xyz_lr(3)))=0;
    }    
  }

  void Counter::reset_beenhere4(){
    ColumnVector xyz(3);
    for(unsigned int i=0;i<m_diff_path.size();i++){
      xyz<<m_diff_path[i](1)<<m_diff_path[i](2)<<m_diff_path[i](3);
      m_beenhere4((int)MISCMATHS::round((float)xyz(1)),
		  (int)MISCMATHS::round((float)xyz(2)),
		  (int)MISCMATHS::round((float)xyz(3)))=0;
    }    
  }


  void Counter::update_matrix4_col(){
    //Tracer_Plus tr("Counter::update_matrix4");
    float d=opts.steplength.value();
    int x,y,z,Conrow4;
    ColumnVector xyz(3);
    if(opts.mask4.value()==""){
      for(unsigned int i=0;i<m_diff_path.size();i++){
	// check here if back to seed
	if(i>0 && (m_path[i]-m_path[0]).MaximumAbsoluteValue()==0)
	  d=opts.steplength.value();
	
	xyz<<m_diff_path[i](1)<<m_diff_path[i](2)<<m_diff_path[i](3);
	x=(int)MISCMATHS::round((float)xyz(1));
	y=(int)MISCMATHS::round((float)xyz(2));
	z=(int)MISCMATHS::round((float)xyz(3));
	Conrow4=m_lookup4(x,y,z);
	
	if(Conrow4>0){
	  if(m_beenhere4(x,y,z)==0){
	    m_ConMat4->AddToTraj(Conrow4,m_curloc.loc+1,d,(int)m_diff_path[i](4));  
	    m_beenhere4(x,y,z)=1;
	    d+=opts.steplength.value();
	  }
	}
      } 
    }
    else{ //case where mask4 is set
      vector<int> rows,cols,vals,crossedrois,crossedlocs;
      vector<float> ds;
      bool restarted=false;int offset=-1;
      vector<ColumnVector> crossedvox;
      vector< pair<int,int> > surf_Triangle;
      for(unsigned int i=0;i<m_diff_path.size();i++){
	// check here if back to seed
	if(i>0 && (m_path[i]-m_path[0]).MaximumAbsoluteValue()==0){
	  d=opts.steplength.value();restarted=true;
	  offset-=1;
	}	
	x=(int)MISCMATHS::round(float(m_diff_path[i](1)));
	y=(int)MISCMATHS::round(float(m_diff_path[i](2)));
	z=(int)MISCMATHS::round(float(m_diff_path[i](3)));

	if(m_beenhere4(x,y,z)==0){
	  rows.push_back(m_lookup4(x,y,z));
	  vals.push_back((int)m_diff_path[i](4));
	  ds.push_back(d);
	  m_beenhere4(x,y,z)=1;
	  d+=opts.steplength.value();
	}

	if(i>0 && !restarted){
	  crossedlocs.clear();crossedrois.clear();
	  if(m_mask4.nSurfs()>0){
	    crossedvox=m_crossedvox[i+offset];
	  }
	  if(m_mask4.has_crossed_roi(m_path[i-1],m_path[i],crossedvox,crossedrois,crossedlocs,surf_Triangle,opts.closestvertex.value())){
	    for(unsigned int j=0;j<crossedlocs.size();j++)
	      cols.push_back(crossedlocs[j]+1);
	  }
	}
	restarted=false;
      }
      // now fill in the matrix
      make_unique(cols);      
      for(unsigned int i=0;i<rows.size();i++){
	if(rows[i]>0){
	  for(unsigned int j=0;j<cols.size();j++){	 
	    m_ConMat4->AddToTraj(rows[i],cols[j],ds[i],vals[i]);
	  }
	}
      }
      
    }

  }


  void Counter::save_total(const int& keeptotal){
    // save total number of particles that made it through the streamlining
    ColumnVector keeptotvec(1);
    keeptotvec(1)=keeptotal;
    write_ascii_matrix(keeptotvec,logger.appendDir("waytotal"));

  }
  void Counter::save_total(const vector<int>& keeptotal){
    // save total number of particles that made it through the streamlining
    ColumnVector keeptotvec(keeptotal.size());
    for (int i=1;i<=(int)keeptotal.size();i++)
      keeptotvec(i)=keeptotal[i-1];
    write_ascii_matrix(keeptotvec,logger.appendDir("waytotal"));
  }

  void Counter::save(){
    if((opts.simpleout.value()||opts.omeanpathlength.value()) && !opts.simple.value()){
      save_pathdist();
    }
    if(opts.network.value()){
      m_stline.save_network_mat();
    }
    if(opts.save_paths.value())
      save_paths();
    if(opts.s2tout.value()){
      save_seedcounts();
    }
    if(opts.matrix1out.value()){
      save_matrix1();
    }
    if(opts.matrix2out.value()){
      save_matrix2();
    }
    if(opts.matrix3out.value()){
      save_matrix3();
    }
    if(opts.matrix4out.value()){
      save_matrix4();
    }
  }
  
  void Counter::save_pathdist(){  
    m_prob.setDisplayMaximumMinimum(m_prob.max(),m_prob.min());
    save_volume(m_prob,logger.appendDir(opts.outfile.value()));

    //m_lastpoint.setDisplayMaximumMinimum(m_lastpoint.max(),m_lastpoint.min());
    //save_volume(m_lastpoint,logger.appendDir("lastpoint"));

    if(opts.pathfile.set()){
      m_prob_alt.save_rois(logger.appendDir(opts.outfile.value())+"_alt");
      //m_prob_alt.save_as_volume(logger.appendDir(opts.outfile.value())+"_alt_vol");
      //m_beenhere_alt.save_rois(logger.appendDir(opts.outfile.value())+"_beenhere");
    }
    if(opts.opathdir.value()){
      volume4D<float> tmplocdir(m_prob.xsize(),m_prob.ysize(),m_prob.zsize(),3);
      copybasicproperties(m_prob,tmplocdir);
      tmplocdir=0;
      SymmetricMatrix Tens(3);
      DiagonalMatrix D;Matrix V;
      for(int z=0;z<m_prob.zsize();z++){
	for(int y=0;y<m_prob.ysize();y++){
	  for(int x=0;x<m_prob.xsize();x++){
	    if(m_prob(x,y,z)==0)continue;
	    Tens<<m_localdir(x,y,z,0)
		<<m_localdir(x,y,z,1)
		<<m_localdir(x,y,z,2)
		<<m_localdir(x,y,z,3)
		<<m_localdir(x,y,z,4)
		<<m_localdir(x,y,z,5);
	    if(m_prob(x,y,z)!=0)Tens=Tens/m_prob(x,y,z);
	    EigenValues(Tens,D,V);
	    for(int t=0;t<3;t++)
	      tmplocdir(x,y,z,t)=V(t+1,3);
	  }
	}
      }
      tmplocdir.setDisplayMaximumMinimum(1,-1);
      save_volume4D(tmplocdir,logger.appendDir(opts.outfile.value()+"_localdir"));
    }
    if(opts.omeanpathlength.value()){
      if(!opts.pathdist.value()){
	for (int z=0; z<m_prob.zsize(); z++) {
          for (int y=0; y<m_prob.ysize(); y++) {
	    for (int x=0; x<m_prob.xsize(); x++) {
	      if(m_prob(x,y,z)){
	        m_prob(x,y,z)=m_prob2(x,y,z)/m_prob(x,y,z);
              }else{
	        m_prob(x,y,z)=0;
	      }
	    }
	  }
        }
      }else{
	for (int z=0; z<m_prob.zsize(); z++) {
          for (int y=0; y<m_prob.ysize(); y++) {
	    for (int x=0; x<m_prob.xsize(); x++) {
	      if(m_prob2(x,y,z)){
	        m_prob(x,y,z)=m_prob(x,y,z)/m_prob2(x,y,z);
              }else{
	        m_prob(x,y,z)=0;
	      }
	    }
	  }
        }
      }
      m_prob.setDisplayMaximumMinimum(m_prob.max(),m_prob.min());
      save_volume(m_prob,logger.appendDir(opts.outfile.value())+"_lengths");
      if(opts.pathfile.set()){
        if(!opts.pathdist.value()){
	  m_prob_alt2.divide_rois(m_prob_alt);
	  m_prob_alt2.save_rois(logger.appendDir(opts.outfile.value())+"_alt_lengths");
	}else{
	  m_prob_alt.divide_rois(m_prob_alt2);
          m_prob_alt.save_rois(logger.appendDir(opts.outfile.value())+"_alt_lengths");
	}
      }
    }
  }
  
  void Counter::save_pathdist(string add){  //for simple mode
    string thisout=opts.outfile.value();
    make_basename(thisout);
    thisout+=add;
    m_prob.setDisplayMaximumMinimum(m_prob.max(),m_prob.min());
    save_volume(m_prob,logger.appendDir(thisout));
    if(opts.pathfile.set()){
      m_prob_alt.save_rois(logger.appendDir(thisout)+"_alt");
    }
    if(opts.omeanpathlength.value()){
      if(!opts.pathdist.value()){
        for (int z=0; z<m_prob.zsize(); z++) {
          for (int y=0; y<m_prob.ysize(); y++) {
	    for (int x=0; x<m_prob.xsize(); x++) {
	      if(m_prob(x,y,z)){
	        m_prob(x,y,z)=m_prob2(x,y,z)/m_prob(x,y,z);
              }else{
	        m_prob(x,y,z)=0;
	      }
	    }
	  }
        }
      }else{
    	for (int z=0; z<m_prob.zsize(); z++) {
          for (int y=0; y<m_prob.ysize(); y++) {
	    for (int x=0; x<m_prob.xsize(); x++) {
	      if(m_prob2(x,y,z)){
	        m_prob(x,y,z)=m_prob(x,y,z)/m_prob2(x,y,z);
              }else{
	        m_prob(x,y,z)=0;
	      }
	    }
	  }
        }
      } 
      m_prob.setDisplayMaximumMinimum(m_prob.max(),m_prob.min());
      save_volume(m_prob,logger.appendDir(thisout)+"_lengths");
      if(opts.pathfile.set()){
	if(!opts.pathdist.value()){
	  m_prob_alt2.divide_rois(m_prob_alt);
	  m_prob_alt2.save_rois(logger.appendDir(thisout)+"_alt_lengths");
	}else{
	  m_prob_alt.divide_rois(m_prob_alt2);
          m_prob_alt.save_rois(logger.appendDir(thisout)+"_alt_lengths");
	}
      }
    }
  }

  void Counter::save_seedcounts(){
    // get target names 
    vector<string> targetnames;
    for(int m=0;m<m_targetmasks.nRois();m++){
      string tmpname=m_targetmasks.get_name(m);
      int pos=tmpname.find("/",0);
      int lastpos=pos;
      while(pos>=0){
	lastpos=pos;
	pos=tmpname.find("/",pos);
	// replace / with _
	tmpname[pos]='_';
      }      
      //only take things after the last pos
      tmpname=tmpname.substr(lastpos+1,tmpname.length()-lastpos-1);
      targetnames.push_back(tmpname);
    }

    for(int m=0;m<m_targetmasks.nRois();m++){
      if(m_s2t_count.nRois()>1){
	for(int i=0;i<m_s2t_count.nRois();i++)
	  m_s2t_count.save_map(i,m,logger.appendDir("seeds_"+num2str(i)+"_to_"+targetnames[m]));
      }
      else{// keep this nomenclature for backward compatibility
	m_s2t_count.save_map(0,m,logger.appendDir("seeds_to_"+targetnames[m]));
      }	
    }

    if(opts.s2tastext.value()){
      write_ascii_matrix(m_s2tastext,logger.appendDir("matrix_seeds_to_all_targets"));
      //write_matrix_as_volume(m_s2tastext,logger.appendDir("matrix_seeds_to_all_targets"));
    }

    if(opts.targetpaths.value()){
      for(unsigned int t=0;t<m_prob_multi.size();t++){
	m_prob_multi[t].setDisplayMaximumMinimum(m_prob_multi[t].max(),m_prob_multi[t].min());
	save_volume(m_prob_multi[t],logger.appendDir("target_paths_"+targetnames[t]));
      }
      if(opts.opathdir.value()){
	volume4D<float> tmplocdir(m_prob.xsize(),m_prob.ysize(),m_prob.zsize(),3);
	copybasicproperties(m_prob,tmplocdir);
	for(unsigned int t=0;t<m_prob_multi.size();t++){
	  tmplocdir=0;
	  SymmetricMatrix Tens(3);
	  DiagonalMatrix D;Matrix V;
	  for(int z=0;z<m_prob.zsize();z++){
	    for(int y=0;y<m_prob.ysize();y++){
	      for(int x=0;x<m_prob.xsize();x++){
		if(m_prob_multi[t](x,y,z)==0)continue;
		Tens<<m_localdir(x,y,z,0)
		    <<m_localdir(x,y,z,1)
		    <<m_localdir(x,y,z,2)
		    <<m_localdir(x,y,z,3)
		    <<m_localdir(x,y,z,4)
		    <<m_localdir(x,y,z,5);
		if(m_prob_multi[t](x,y,z)!=0)Tens=Tens/m_prob(x,y,z);
		EigenValues(Tens,D,V);
		for(int tt=0;tt<3;tt++)
		  tmplocdir(x,y,z,tt)=V(tt+1,3);
	      }
	    }
	  }
	  tmplocdir.setDisplayMaximumMinimum(1,-1);
	  save_volume4D(tmplocdir,logger.appendDir("target_localdir_"+targetnames[t]));
	}
      }
    }
    if(opts.omeanpathlength.value()){
      if(!opts.pathdist.value()){
        m_s2t_count2.divide_maps(m_s2t_count);
        for(int m=0;m<m_targetmasks.nRois();m++){
          if(m_s2t_count2.nRois()>1){
	    for(int i=0;i<m_s2t_count2.nRois();i++)
	      m_s2t_count2.save_map(i,m,logger.appendDir("seeds_"+num2str(i)+"_to_"+targetnames[m]+"_lengths"));
          }else{// keep this nomenclature for backward compatibility
	    m_s2t_count2.save_map(0,m,logger.appendDir("seeds_to_"+targetnames[m]+"_lengths"));
          }	
        }
      }else{
	m_s2t_count.divide_maps(m_s2t_count2);
        for(int m=0;m<m_targetmasks.nRois();m++){
          if(m_s2t_count.nRois()>1){
	    for(int i=0;i<m_s2t_count.nRois();i++)
	      m_s2t_count.save_map(i,m,logger.appendDir("seeds_"+num2str(i)+"_to_"+targetnames[m]+"_lengths"));
          }else{// keep this nomenclature for backward compatibility
	    m_s2t_count.save_map(0,m,logger.appendDir("seeds_to_"+targetnames[m]+"_lengths"));
          }	
        }
      }
      if(opts.s2tastext.value()){
	if(!opts.pathdist.value()){
          for(int i=1;i<=m_s2tastext.Nrows();i++){
	    for(int j=1;j<=m_s2tastext.Ncols();j++){
	      if(m_s2tastext(i,j))
              	m_s2tastext(i,j)=m_s2tastext2(i,j)/m_s2tastext(i,j);
	      else
		m_s2tastext(i,j)=0;
            }
          }
	}else{
	  for(int i=1;i<=m_s2tastext.Nrows();i++){
	    for(int j=1;j<=m_s2tastext.Ncols();j++){
	      if(m_s2tastext2(i,j))
                m_s2tastext(i,j)=m_s2tastext(i,j)/m_s2tastext2(i,j);
	      else
		m_s2tastext(i,j)=0;
            }
          }
	}
	write_ascii_matrix(m_s2tastext,logger.appendDir("matrix_seeds_to_all_targets_lengths"));
      }
      if(opts.targetpaths.value()){
	if(!opts.pathdist.value()){
          for(unsigned int t=0;t<m_prob_multi.size();t++){
            for (int z=0; z<m_prob_multi[t].zsize(); z++) {
              for (int y=0; y<m_prob_multi[t].ysize(); y++) {
	        for (int x=0; x<m_prob_multi[t].xsize(); x++) {
	          if(m_prob_multi[t](x,y,z)){
	            m_prob_multi[t](x,y,z)=m_prob_multi2[t](x,y,z)/m_prob_multi[t](x,y,z);
                  }else{
	            m_prob_multi[t](x,y,z)=0;
	          }
	        }
	      }
            }
	    m_prob_multi[t].setDisplayMaximumMinimum(m_prob_multi[t].max(),m_prob_multi[t].min());
	    save_volume(m_prob_multi[t],logger.appendDir("target_paths_"+targetnames[t]+"_lengths"));
          }
	}else{
	  for(unsigned int t=0;t<m_prob_multi.size();t++){
            for (int z=0; z<m_prob_multi[t].zsize(); z++) {
              for (int y=0; y<m_prob_multi[t].ysize(); y++) {
	        for (int x=0; x<m_prob_multi[t].xsize(); x++) {
	          if(m_prob_multi2[t](x,y,z)){
	            m_prob_multi[t](x,y,z)=m_prob_multi[t](x,y,z)/m_prob_multi2[t](x,y,z);
                  }else{
	            m_prob_multi[t](x,y,z)=0;
	          }
	        }
	      }
            }
	    m_prob_multi[t].setDisplayMaximumMinimum(m_prob_multi[t].max(),m_prob_multi[t].min());
	    save_volume(m_prob_multi[t],logger.appendDir("target_paths_"+targetnames[t]+"_lengths"));
          }
	}
      }
    }
  }
    
  // the following is a helper function for save_matrix*
  //  to convert between nifti coords (external) and newimage coord (internal)
  void applycoordchange(volume<int>& coordvol, const Matrix& old2new_mat)
  {
    for (int n=0; n<=coordvol.maxx(); n++) {
      ColumnVector v(4);
      v << coordvol(n,0,0) << coordvol(n,1,0) << coordvol(n,2,0) << 1.0;
      v = old2new_mat * v;
      coordvol(n,0,0) = MISCMATHS::round(v(1));
      coordvol(n,1,0) = MISCMATHS::round(v(2));
      coordvol(n,2,0) = MISCMATHS::round(v(3));
    }
  }
  void applycoordchange(Matrix& coordvol, const Matrix& old2new_mat)
  {
    for (int n=1; n<=coordvol.Nrows(); n++) {
      ColumnVector v(4);
      v << coordvol(n,1) << coordvol(n,2) << coordvol(n,3) << 1.0;
      v = old2new_mat * v;
      coordvol(n,1) = MISCMATHS::round(v(1));
      coordvol(n,2) = MISCMATHS::round(v(2));
      coordvol(n,3) = MISCMATHS::round(v(3));
    }
  }

  void Counter::save_matrix1(){
    if(!opts.omeanpathlength.value()){
      m_ConMat1->Print(logger.appendDir("fdt_matrix1.dot"));
    }else{
      if(!opts.pathdist.value()){
	m_ConMat1b->Print(logger.appendDir("fdt_matrix1.dot"));
      }else{
	m_ConMat1->Print(logger.appendDir("fdt_matrix1.dot"));
      }
      for(unsigned int i=1;i<=m_ConMat1->Nrows();i++){
	for(unsigned int j=1;j<=m_ConMat1->Ncols();j++){
	  if(m_ConMat1b->Peek(i,j)){
	    m_ConMat1->Set(i,j,m_ConMat1->Peek(i,j)/m_ConMat1b->Peek(i,j));
          }
        }
      }
      m_ConMat1->Print(logger.appendDir("fdt_matrix1_lengths.dot"));
    }
  }

  void Counter::save_matrix2(){
    if(!opts.omeanpathlength.value()){
      m_ConMat2->Print(logger.appendDir("fdt_matrix2.dot"));
    }else{
      if(!opts.pathdist.value()){
	m_ConMat2b->Print(logger.appendDir("fdt_matrix2.dot")); 
      }else{
	m_ConMat2->Print(logger.appendDir("fdt_matrix2.dot")); 
      }
      for(unsigned int i=1;i<=m_ConMat2->Nrows();i++){
	for(unsigned int j=1;j<=m_ConMat2->Ncols();j++){
	  if(m_ConMat2b->Peek(i,j)){
	    m_ConMat2->Set(i,j,m_ConMat2->Peek(i,j)/m_ConMat2b->Peek(i,j));
          }
        }
      }
      m_ConMat2->Print(logger.appendDir("fdt_matrix2_lengths.dot"));
    }
  }

  void Counter::save_matrix3(){
    if(!opts.omeanpathlength.value()){
      m_ConMat3->Print(logger.appendDir("fdt_matrix3.dot"));
    }else{
      if(!opts.pathdist.value()){
	m_ConMat3b->Print(logger.appendDir("fdt_matrix3.dot"));
      }else{
	m_ConMat3->Print(logger.appendDir("fdt_matrix3.dot"));
      }
      if(opts.omeanpathlength.value()){
        for(unsigned int i=1;i<=m_ConMat3->Nrows();i++){
	  for(unsigned int j=1;j<=m_ConMat3->Ncols();j++){
	    if(m_ConMat3b->Peek(i,j)){
	      m_ConMat3->Set(i,j,m_ConMat3->Peek(i,j)/m_ConMat3b->Peek(i,j));
            }
          }
        }
      }
      m_ConMat3->Print(logger.appendDir("fdt_matrix3_lengths.dot"));
    }
  }
  void Counter::save_matrix4(){
    m_ConMat4->SaveTrajFile(logger.appendDir("fdt_matrix4_"));    
  }

  // this function returns the total number of pathways that survived streamlining 
  int Seedmanager::run(const float& x,const float& y,const float& z,
		       bool onewayonly, int fibst,float sampvox){
    //Tracer_Plus tr("Seedmanager::run");

    //onewayonly for mesh things..
    if(opts.fibst.set()){
      fibst=opts.fibst.value()-1;      
    }
    else{
      if(fibst == -1){
	fibst=0;
      }
      //TB moved randfib option inside tractvols.h 28/10/2009
      // This means that we have access to fsamples when figuring out fibst
      // so we can choose to seed in proportion to f in that voxel. 
    }
  
  
    int nlines=0;
    for(int p=0;p<opts.nparticles.value();p++){
      if(opts.randfib.value()>2){ 
	//This bit of code just does random sampling from all fibre populations - even those whose f value is less than f-thresh. 
	//3 other possibilities - randfib==0 -> use fibst (default first fibre but can be set)
	// randfib==1 - random sampling of fibres bigger than fthresh
	// randfib==2 random sampling of fibres bigger than fthresh in proporthion to their f-values. 
	double tmp=(rand()/(double(RAND_MAX)+1)) * m_counter.get_stline().nfibres();
	fibst = (int)floor(tmp);
      }
    
      // random jitter of seed point inside a sphere
      float newx=x,newy=y,newz=z;    
      if(sampvox>0){
	bool rej=true;float dx,dy,dz;float r2=sampvox*sampvox;
	while(rej){
	  dx=2.0*sampvox*((float)rand()/float(RAND_MAX)-.5);
	  dy=2.0*sampvox*((float)rand()/float(RAND_MAX)-.5);
	  dz=2.0*sampvox*((float)rand()/float(RAND_MAX)-.5);
	  if( dx*dx+dy*dy+dz*dz <= r2 )
	    rej=false;
	}
	newx+=dx/m_seeddims(1);
	newy+=dy/m_seeddims(2);
	newz+=dz/m_seeddims(3);
      }
    
      if(opts.verbose.value()>1)
	logger.setLogFile("particle"+num2str(p));
    
      m_counter.get_stline().reset(); //This now includes a vols.reset() in order to get fibst right. 

      bool forwardflag=false,backwardflag=false;    
      int  rejflag1=1,rejflag2=1; // 0:accept, 1:reject, 2:wait
      m_counter.clear_path();
    
      // track in one direction
      if(!onewayonly || opts.matrix3out.value()){//always go both ways in matrix3 mode
	rejflag1 = m_counter.get_stline().streamline(newx,newy,newz,m_seeddims,fibst);

	if(rejflag1==0 || rejflag1==2){ 
	  forwardflag=true;
	  m_counter.append_path();
	  if(opts.save_paths.value())
	    m_counter.add_path();
	}
        if(rejflag1==0 && opts.matrix3out.value()){ 
	  m_counter.get_stline().copy_inmask3();
	  m_counter.get_stline().reset_m_inmask3_aux();
	}else if(rejflag1==1 && opts.matrix3out.value()){
	  m_counter.get_stline().reset_m_inmask3_aux();
	}
	m_counter.get_stline().reverse();
      }
      

      // track in the other direction
      rejflag2=m_counter.get_stline().streamline(newx,newy,newz,m_seeddims,fibst);

      if(rejflag2==0){	
	backwardflag=true;
        if(opts.matrix3out.value())
	  m_counter.get_stline().copy_inmask3();
      }
      if(rejflag2>0){
	backwardflag=false;
	if(rejflag1>0)
	  forwardflag=false;
      }
    
      if(!forwardflag)
	m_counter.clear_path();
      if(backwardflag){
	m_counter.append_path();   	
	if(opts.save_paths.value())
	  m_counter.add_path();
      }
      if(forwardflag || backwardflag){
	nlines++; 
	m_counter.count_streamline();
      
	if(opts.matrix3out.value()||opts.matrix1out.value()){
	  float pathlength;
	  if( !forwardflag || !backwardflag)
	    pathlength = m_counter.calc_pathlength(0);
	  else // if the paths are appended, one step has zero length
	    pathlength = m_counter.calc_pathlength(1);
	
	  if(opts.matrix3out.value() && pathlength>=opts.distthresh3.value())
	    m_counter.update_matrix3();
	  
	  if(opts.matrix1out.value() && pathlength>=opts.distthresh1.value())
	    m_counter.update_matrix1();
	}
      }

      m_counter.clear_streamline(); 
    }
  
    m_counter.count_seed();
  
  
    return nlines;
    
  }

}

  
  

