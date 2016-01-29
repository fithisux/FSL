/*  fdt_matrix_split.cc

    Saad Jbabdi, FMRIB Image Analysis Group

    Copyright (C) 2010 University of Oxford */

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
#include "newimage/newimageio.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/SpMat.h"

using namespace Utilities;
using namespace NEWIMAGE;
using namespace MISCMATHS;

string title="Splits matrix3 output into seeds_to_target files";
string examples="fdt_matrix_split [options]";

Option<bool> verbose(string("-v,--verbose"),false,
		       string("switch on diagnostic messages"),
		       false,no_argument);
Option<bool> help(string("-h,--help"),false,
		       string("display this message"),
		       false,no_argument);
Option<string> ptxdir(string("-d,--ptxdir"),string(""),
		      string("probtrackx output directory"),
		      true,requires_argument);
Option<string> oBase(string("-o,--output"),string(""),
		       string("output basename (results will be <outbase>_<seed>_to_<target>)"),
		       true,requires_argument);
Option<string> seedsFileName(string("-s,--seeds"),string(""),
		       string("filename for ascii list of seeds"),
		       true,requires_argument);
Option<string> targetsFileName(string("-t,--targets"),string(""),
		       string("filename for ascii list of targets"),
		       true,requires_argument);
Option<string> refFileName(string("--refvol"),string(""),
			   string("reference volume (required if only using surfaces)"),
			   false,requires_argument);
Option<string> mask3FileName(string("-m,--mask3"),string(""),
		       string("filename for mask3 used in probtrackx"),
		       true,requires_argument);
Option<string> meshspace(string("--meshspace"),string("freesurfer"),
			  string("mesh space for surfaces. Either 'freesurfer' (default) or 'first' or 'caret'"),
			  false,requires_argument);


////////////////////////////////////////////////////////
void readMasks(vector<string>& masks, const string& filename){
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

void vol2vol(const volume<float>& mask,volume<float>& lu){
  lu.reinitialize(mask.xsize(),mask.ysize(),mask.zsize());
  lu=-1;
  int counter=0;
  for(int z=mask.minz();z<=mask.maxz();z++)
    for(int y=mask.miny();y<=mask.maxy();y++)
      for(int x=mask.minx();x<=mask.maxx();x++){
	if(mask(x,y,z)==0)continue;
	lu(x,y,z)=counter;
	counter++;
      }
}
void vol2mat(const volume<float>& mask,Matrix& coords,const volume<float>& lu,const volume<float>& mask3,const Matrix& xfm){
  int counter=0;
  ColumnVector seed_xyz(3),mask3_xyz(3),seed_dims(3),mask3_dims(3);
  seed_dims << mask.xdim() << mask.ydim() << mask.zdim();
  mask3_dims << mask3.xdim() << mask3.ydim() << mask3.zdim();
  for(int z=mask.minz();z<=mask.maxz();z++)
    for(int y=mask.miny();y<=mask.maxy();y++)
      for(int x=mask.minx();x<=mask.maxx();x++){
	if(mask(x,y,z)==0)continue;
	counter++;
      }
  coords.ReSize(counter,4);
  counter=0;
  int xx,yy,zz;
  for(int z=mask.minz();z<=mask.maxz();z++)
    for(int y=mask.miny();y<=mask.maxy();y++)
      for(int x=mask.minx();x<=mask.maxx();x++){
	if(mask(x,y,z)==0)continue;
	counter++;
	seed_xyz << x << y << z;
	mask3_xyz = vox_to_vox(seed_xyz,seed_dims,mask3_dims,xfm);	
	xx=(int)round(float(mask3_xyz(1)));
	yy=(int)round(float(mask3_xyz(2)));
	zz=(int)round(float(mask3_xyz(3)));
	if(mask3(xx,yy,zz)!=0)
	  coords(counter,1) = lu(xx,yy,zz);
	else
	  coords(counter,1) = -1; // seed voxel not found in mask3
	coords(counter,2) = x;
	coords(counter,3) = y;
	coords(counter,4) = z;
      }  
}

void getName(const string& istring,string& ostring){
  ostring=istring;
  int pos=istring.find("/",0);
  int lastpos=pos;
  while(pos>=0){
    lastpos=pos;
    pos=ostring.find("/",pos);
    ostring[pos]='_';
  }
  ostring=ostring.substr(lastpos+1,ostring.length()-lastpos-1);
}

int do_split_matrix(){
  if(verbose.value())
    cout<<"read all inputs"<<endl;
  
  SpMat<int> matrix3(ptxdir.value()+"/fdt_matrix3.dot");
  Matrix     coords3 = read_ascii_matrix(ptxdir.value()+"/coords_for_fdt_matrix3");

  Matrix seed2mask3 = IdentityMatrix(4);
  if(xfmFileName.value()!="")
    seed2mask3 = read_ascii_matrix(xfmFileName.value());
  

  CSV *mask3;
  {
    if(refFileName.value()!=0){
      volume<int> tmprefvol;
      read_volume(tmprefvol,refFileName.value());
      mask3 = new CSV(tmprefvol);
    }
    else{
      mask3 = new CSV();
    }
    mask3->load_rois(mask3FileName.value());
    mask3->set_convention(meshspace.value());
    if(mask3->nVols()==0 && refFileName.value()==""){
      cerr<<"Must provide a refvol for mask3 if only using surfaces"<<endl;
      exit(1);
    }
  }

  CSV *seeds,*targets;
  {
    if(refFileName.value()!=0){
      volume<int> tmprefvol;
      read_volume(tmprefvol,refFileName.value());
      seeds   = new CSV(tmprefvol);
      targets = new CSV(tmprefvol);
    }
    else{
      seeds   = new CSV();
      targets = new CSV();
    }
    seeds->load_rois(seedsFileName.value());
    targets->load_rois(seedsFileName.value());

    if(refFileName.value()==""){
      if(seeds->nVols()!=0){
	if(targets->nVols()==0)
	  targets->change_refvol(seeds->get_refvol());
      }
      else{
	if(targets->nVols()!=0)
	  seeds->change_refvol(targets->get_refvol());
	else{
	  cerr<<"Must provide a refvol for seeds/targets if only using surfaces"<<endl;
	  exit(1);
	}
      }
    }
    seeds->set_convention(meshspace.value());
    targets->set_convention(meshspace.value());
  }


  if(verbose.value())
    cout<<"split the matrix"<<endl;


  // iterate over seeds
  // iterate over targets
  // find coords for seed and target
  // find corresponding coords in matrix3? Argh, they must be in the same space! that's fine
  // and seed and target need to have correspondants in mask3 (when they are surfaces)
  // split the matrix
  
  vector<CSV*> out(targets->nRois());
  for(int i=0;i<targets->nRois();i++){
    out[i] = new CSV(seeds->get_refvol());
    out[i]->reset_values();
  }


  ColumnVector coords(3);vector<int> tlocs;int x,y,z;
  for(int si=0;si<seeds->nRois();si++){
    for(int st=0;st<targets->nRois();st++){
      string sname = seeds->get_name(si);
      string tname = targets->get_name(st);

      // argh! i need to loop through locations within that volume/surface, not all rois
      for(int sloc=1;sloc<=seeds->nLocs();sloc++){
	int isloc=intersectRoi(mask3,seeds,sloc);
	if(isloc<0)continue;
	tlocs.clear();
	for(int tloc=1;tloc<=targets->nLocs();tloc++){
	  int itloc = intersectRoi(mask3,targets,tloc);
	  if(iloc<0)continue;
	  tlocs.push_back(iloc);
	}
	// fill with matrix values! (remember the matrix is symmetric but only half filled)
	for(unsigned int i=0;i<tlocs.size();i++){
	  out->add_value(type,si,sloc,matrix3->Peek(isloc,tlocs[i]));
	  out->add_value(type,si,sloc,matrix3->Peek(tlocs[i],isloc));
	}
      }
    }

  }


  for(int i=0;i<targets->nRois();i++){
    tlocs.clear();
    for(int loc=1;loc<=targets->nLocs();loc++){
      if(targets->isVol(i)){
	coords = targets->get_loc_coords(loc);
	x=(int)round((float)coords(1));
	y=(int)round((float)coords(2));
	z=(int)round((float)coords(3));

	if(!mask3->isInRoi(x,y,z))continue;
	for(int j=0;j<mask3->nVols();j++)
	  tlocs.push_back(mask3->get_loc(j,x,y,z));

      }
      else{
      }
    }

  }
  
  copybasicproperties(mask3,tmpvol);
  string sname,tname;
  for(unsigned int i=0;i<seedCoords.size();i++){
    if(verbose.value())
      cout<<"Seed: "<<seedFiles[i]<<endl;
    getName(seedFiles[i],sname);
    for(unsigned int j=0;j<targetCoords.size();j++){
      if(verbose.value())
	cout<<"Target: "<<targetFiles[j]<<endl;
      getName(targetFiles[j],tname);

      tmpvol=0;
      for(int sr=1;sr<=seedCoords[i].Nrows();sr++){
	if(seedCoords[i](sr,1)==-1)continue;
	for(int tr=1;tr<=targetCoords[j].Nrows();tr++){
	  if(targetCoords[j](tr,1)==-1)continue;
	  int sc=(int)seedCoords  [i](sr,1)+1;
	  int tc=(int)targetCoords[j](tr,1)+1;
	  if(!isVol){
	    tmpvol((int)seedCoords[i](sr,2),
		   (int)seedCoords[i](sr,3),
		   (int)seedCoords[i](sr,4))+=matrix3->Peek(sc,tc);
	    tmpvol((int)seedCoords[i](sr,2),
		   (int)seedCoords[i](sr,3),
		   (int)seedCoords[i](sr,4))+=matrix3->Peek(tc,sc);
	  }
	  else{
	    tmpvol((int)seedCoords[i](sr,2),
		   (int)seedCoords[i](sr,3),
		   (int)seedCoords[i](sr,4))+=matrix3vol(tc-1,sc-1,0);	  
	  }
	}
      }
      // save
      tmpvol.setDisplayMaximumMinimum(tmpvol.max(),0);
      save_volume(tmpvol,oBase.value()+"_"+sname+"_to_"+tname);
    }
  }
  if(verbose.value())
    cout<<"Done!"<<endl;

  return 0;
}

int main(int argc,char *argv[]){

  Tracer tr("main");
  OptionParser options(title,examples);

  try{
    options.add(verbose);
    options.add(help);
    options.add(matFileName);
    options.add(oBase);
    options.add(seedsFileName);
    options.add(targetsFileName);
    options.add(refFileName);
    options.add(mask3FileName);
    options.add(meshspace);

    options.parse_command_line(argc,argv);

    
    if ( (help.value()) || (!options.check_compulsory_arguments(true)) ){
      options.usage();
      exit(EXIT_FAILURE);
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
  
  return do_split_matrix();
  
  
}
