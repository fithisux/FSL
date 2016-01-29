/*  fsl_glm - 

    Christian F. Beckmann, FMRIB Analysis Group

    Copyright (C) 2008-2013 University of Oxford */

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

#include "libvis/miscplot.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include "utils/options.h"
#include <vector>
#include "newimage/newimageall.h"
#include "melhlprfns.h"

using namespace MISCPLOT;
using namespace MISCMATHS;
using namespace Utilities;
using namespace std;

// The two strings below specify the title and example usage that is
// printed out as the help or usage message

  string title=string("fsl_sbca (Version 1.0)")+
		string("\nAuthor: Christian F. Beckmann \nCopyright(C) 2008-2013 University of Oxford \n")+
		string(" \n Performs seed-based correlation analysis on FMRI data\n")+
		string(" using either a single seed coordinate or a seed mask \n")+
		string(" ");
  string examples="fsl_sbca -i <input> -s <seed> -t <target> -o <basename> [options]";

//Command line Options {
  Option<string> fnin(string("-i,--in"), string(""),
		string("        input file name (4D image file)"),
		true, requires_argument);
  Option<string> fnout(string("-o,--out"), string(""),
		string("output file base name"),
		true, requires_argument);
  Option<string> fnseed(string("-s,--seed"), string(""),
		string("seed voxel coordinate or file name of seed mask (3D/4D file)"),
		true, requires_argument);
  Option<string> fntarget(string("-t,--target"), string(""),
		string("file name of target mask(s) (3D or 4D file)"),
		true, requires_argument);
  Option<bool> regress_only(string("-r,--reg"), false,
		string("perform time series regression rather than classification to targets"),
		false, no_argument);
  Option<string> fnconf(string("--conf"), string(""),
		string("        file name (or comma-separated list of file name) for confound ascii txt files"),
		false, requires_argument);
  Option<string> fnseeddata(string("--seeddata"), string(""),
		string("file name of 4D data file for the seed"),
		false, requires_argument);		
  Option<bool> map_bin(string("--bin"), false,
		string("        binarise spatial maps prior to calculation of time courses"),
		false, no_argument);
  Option<bool> verbose(string("-v,--verbose"), false,
		string("switch on diagnostic messages"),
		false, no_argument);
  Option<bool> tc_mean(string("--mean"), false,
		string("        use mean instead of Eigenvariates for calculation of time courses"),
		false, no_argument);
  Option<int> tc_order(string("--order"), 1,
		string("        number of Eigenvariates (default 1)"),
		false, requires_argument);	
  Option<bool> abscc(string("--abscc"), false,
		string("        use maximum absolute value instead of of maximum value of the cross-correlations"),
		false, no_argument);			
  Option<bool> out_seeds(string("--out_seeds"), false,
		string("output seed mask image as <basename>_seeds"),
		false, no_argument);
  Option<bool> out_seedmask(string("--out_seedmask"), false,
		string("output seed mask image as <basename>_seedmask"),
		false, no_argument);
  Option<bool> out_ttcs(string("--out_ttcs"), false,
		string("output target time courses as <basename>_ttc<X>.txt"),
		false, no_argument);
  Option<bool> out_conf(string("--out_conf"), false,
		string("output confound time courses as <basename>_confounds.txt"),
		false, no_argument);
  Option<bool> out_tcorr(string("--out_tcorr"), false,
		string("output target correlations as <basename>_tcorr.txt"),
		false, no_argument, false);
  Option<int> help(string("-h,--help"), 0,
		string("display this help text"),
		false,no_argument);
  Option<bool> debug(string("--debug"), false,
		string("        switch on debug messages"),
		false, no_argument, false);
		/*
}
*/
//Globals {
	Matrix data, confounds;
	volume4D<float> orig_data;
	volume<float> maskS, maskT; 
	int voxels = 0;
	Matrix seeds, coords;
	vector<Matrix> ttcs;
	
	Matrix out1, out2;
	 /*
}
*/
////////////////////////////////////////////////////////////////////////////

// Local functions
void save4D(Matrix what, volume<float>& msk, string fname){
  if(debug.value()) 
	cerr << "DBG: in save4D" << endl;
  volume4D<float> tempVol;
  tempVol.setmatrix(what,msk);
  save_volume4D(tempVol,fname);
}

void save4D(volume<float>& in, string fname){
  if(debug.value()) 
	cerr << "DBG: in save4D" << endl;
  volume4D<float> tempVol;
  tempVol.addvolume(in);
  save_volume4D(tempVol,fname);
}

ReturnMatrix create_coords(string what){
  if(debug.value()) 
	cerr << "DBG: in create_coords" << endl;
  Matrix res;
  ifstream fs(what.c_str());
  if (!fs) { 
	Matrix tmp(1,3);
	char *p;
	char t[1024];
	const char *discard = ", [];{(})abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ~!@#$%^&*_-=+|\':><./?";	  
	strcpy(t, what.c_str());
	p=strtok(t,discard);
	tmp(1,1) = atoi(p);
	p=strtok(NULL,discard);
	tmp(1,2) = atoi(p);
	p=strtok(NULL,discard);
	tmp(1,3) = atoi(p);
	res = tmp;
	
	do{
	  p=strtok(NULL,discard);
	  if(p){
		tmp(1,1) = atoi(p);
		p=strtok(NULL,discard);
		tmp(1,2) = atoi(p);
		p=strtok(NULL,discard);
		tmp(1,3) = atoi(p);
		res &= tmp;   	
      }
    }while(p);			
  }else{
    res = read_ascii_matrix(fs);
    fs.close();
  }

  if(res.Ncols()!=3){
	cerr << "ERROR: incorrect format " << what << endl;
	exit(1);
  }

 // if(verbose.value())
//	cout << " Created seed coordinates (size: " << res.Nrows() << " x " << res.Ncols() << ")" << endl;
  res.Release();
	
  return res;	
}

void create_mask(string what){
  if(debug.value()) 
	cerr << "DBG: in create_mask" << endl;
  
  coords = create_coords(what);
  maskS = orig_data[0] * 0.0;
  for(int ctr = 1; ctr <= coords.Nrows(); ctr++)
	maskS(coords(ctr,1),coords(ctr,2),coords(ctr,3)) = 1.0;
  maskS.binarise(1e-8);
}	

void create_seeds(string what){
  if(debug.value()) 
	cerr << "DBG: in create_seeds" << endl;


  volume4D<float> tmp_vol;
  if(fsl_imageexists(what)){
	read_volume4D(tmp_vol,what);
	maskS = tmp_vol[0];
    if(!samesize(orig_data[0],maskS)){
	  cerr << "ERROR: Seed mask image does not match input image" << endl;
      exit(1);
    }
  }  
  else
    create_mask(what);	

  if(tmp_vol.tsize() > 1 && tmp_vol.tsize() == orig_data.tsize()){
	maskS *= tmp_vol[0] / tmp_vol.tsize();
	for(int ctr=1; ctr < tmp_vol.tsize(); ctr++)
	  maskS += tmp_vol[ctr] * tmp_vol[ctr] / tmp_vol.tsize(); 
    maskS.binarise(1e-8);
	seeds = remmean(tmp_vol.matrix(maskS),1);
  }
  else{    
    volume4D<float> tmp_mask;
    tmp_mask.addvolume(maskS);
    maskS.binarise(1e-8);

    if(fnseeddata.value()>"" && fsl_imageexists(fnseeddata.value())){
      volume4D<float> seed_data;
	  if(verbose.value())
		cout << " Reading input data for seeds " << fnseeddata.value() << endl;
      read_volume4D(seed_data,fnseeddata.value());
      seeds = remmean(seed_data.matrix(maskS),1);
    }else{
	  seeds = remmean(orig_data.matrix(maskS),1);	
      if(!map_bin.value()){
    	Matrix scales = tmp_mask.matrix(maskS);	
        seeds = SP(seeds, ones(seeds.Nrows(),1) * scales);
      }
    }
  }
	
  voxels = seeds.Ncols();

  if(debug.value()){
	cerr << "DBG: " << voxels << " voxels" << endl;
	cerr << "DBG: seeds matrix is " << seeds.Nrows() << " x " << seeds.Ncols() << endl;
  }

  if(verbose.value())
	cout << " Created seed time courses "  << endl;

}

ReturnMatrix create_confs(string what){
  if(debug.value()) 
	cerr << "DBG: in create_confs" << endl;
  
  Matrix res, tmp;    
  char *p;
  char t[1024];
  const char *discard = ",";

  strcpy(t, what.c_str());
  p=strtok(t,discard);
  res = remmean(read_ascii_matrix(string(p)),1);
  do{
    p=strtok(NULL,discard);
    if(p){
		tmp = read_ascii_matrix(string(p));
		if(tmp.Nrows()!=res.Nrows()){
			cerr << "ERROR: confound matrix" << string(p) << " is of wrong size "<< endl;
			exit(1);
		}
    	res |= remmean(tmp,1);	
    }
  }while(p);

  if(verbose.value())
	cout << " Created confound matrix (size: " << res.Nrows() << " x " << res.Ncols() << ")" << endl;
  res.Release();
  return res;
}

ReturnMatrix calc_ttc(volume<float>& in){
	if(debug.value()) 
	cerr << "DBG: in calc_ttc" << endl;
  
	Matrix res, tmp, scales;
	
	volume<float> tmp1;
	volume4D<float> tmp2;

	tmp1 = in;
	tmp1.binarise(1e-8);
	maskT += tmp1;
	
	tmp2.addvolume(in);	
	scales = tmp2.matrix(tmp1);
	tmp = remmean(orig_data.matrix(tmp1),1);
	
	if(!map_bin.value())
		tmp = SP(tmp, ones(tmp.Nrows(),1) * scales);
	
	if(tc_mean.value())
		res = mean(tmp,2);
	else{
	   	SymmetricMatrix Corr;
		Corr << tmp * tmp.t() / tmp.Ncols();
		DiagonalMatrix tmpD;
	    EigenValues(Corr,tmpD,res);	
		res = fliplr(res.Columns(res.Ncols()-tc_order.value()+1 , res.Ncols())) * std::sqrt(tmp.Nrows());	
		
		Matrix res2 = mean(tmp,2);

		if(debug.value())
			cerr << "DBG: mean size is " << res2.Nrows() << " x " << res2.Ncols() << endl;
		res2 = res2.Column(1).t() * res.Column(1);
		
		if((float)res2.AsScalar() < 0){
			res = -1.0 * res;
			if(debug.value())
				cerr << "DBG: flipping first eigenvariates" << endl;
		}
	}
	
	if(debug.value())
		cerr << "DBG: size is " << res.Nrows() << " x " << res.Ncols() << endl;
	res.Release();
	return res;
}

void create_target_tcs(){
	if(debug.value()) 
	cerr << "DBG: in create_target_tcs" << endl;
  
	volume4D<float> tmptarg;
	read_volume4D(tmptarg,fntarget.value());	
    maskT = orig_data[0] * 0.0;

	for(int ctr=0; ctr < tmptarg.tsize(); ctr++){
	   ttcs.push_back(calc_ttc(tmptarg[ctr]));
	}	
	
	if(debug.value()) {
		cerr << "DBG: " << ttcs.size() << " target matrices created " << endl;
	}
	if(verbose.value())
	  cout << " Created target mask time courses "  << endl;
	
}

int setup(){
  if(debug.value()) 
	cerr << "DBG: in setup" << endl;
  if(fsl_imageexists(fnin.value())){ //read data
	if(verbose.value())
      cout << " Reading input file " << fnin.value() << endl;
    read_volume4D(orig_data,fnin.value());
  }
  else{
	cerr << "ERROR: Invalid input file " << fnin.value() << endl;       
	exit(1);
    }

  create_seeds(fnseed.value());

  if(!regress_only.value())
    create_target_tcs();
  else{
   	volume4D<float> tmptarg;
	read_volume4D(tmptarg,fntarget.value());	
	maskT = tmptarg[0];
	maskT.binarise(1e-8);
	data = orig_data.matrix(maskT);
	data = remmean(data,1);
  }

  if(fnconf.value()>"")
    confounds = create_confs(fnconf.value());

  return 0;	
}

ReturnMatrix calc_tcorr(int in){
  if(debug.value()) 
	cerr << "DBG: in calc_tcorr" << endl;

	Matrix res = zeros(1,seeds.Ncols()), partial_conf, targetcol;
	
	for(int ctr = 0; ctr < (int)ttcs.size(); ctr++)
	  if(ctr != in){
	  	if(partial_conf.Storage() == 0)
	    	partial_conf = ttcs.at(ctr);
	  	else	
    		partial_conf |= ttcs.at(ctr);
      }

	if(ttcs.at(in).Ncols()>1) {
	  if(partial_conf.Storage()>0)
	    partial_conf = ttcs.at(in).Columns(2,ttcs.at(in).Ncols()) | partial_conf;
	  else
	    partial_conf = ttcs.at(in).Columns(2,ttcs.at(in).Ncols());
	}
 	if(confounds.Storage() > 0) {
	  if(partial_conf.Storage()>0)
	    partial_conf |= confounds;	
	  else
	    partial_conf = confounds;	
	}
    if(debug.value() && partial_conf.Storage()>0) 
      cerr << "DBG: partial_conf " << partial_conf.Nrows() << " x " << partial_conf.Ncols() << endl;

	targetcol = ttcs.at(in).Column(1);
    if(debug.value()) 
      cerr << "DBG: targetcol " << targetcol.Nrows() << " x " << targetcol.Ncols() << endl;

	for(int ctr = 1; ctr <= seeds.Ncols(); ctr++)
      res(1,ctr) = Melodic::corrcoef(targetcol, seeds.Column(ctr), partial_conf).AsScalar();
    
	res.Release();
	return res;	
}

void calc_res(){
  if(debug.value()) 
	cerr << "DBG: in calc_res" << endl;

  out2 = zeros(1,seeds.Ncols());	

  if(!regress_only.value()){
	//Target TCs exist
	if(verbose.value())
	  cout << " Calculating partial correlation scores between seeds and targets "  << endl;

	Matrix tmp;
	int tmp2;
	out1=zeros(ttcs.size(),seeds.Ncols());
	for(int ctr = 0 ;ctr < (int)ttcs.size(); ctr++)
	  out1.Row(ctr+1) = calc_tcorr(ctr);					

    for(int ctr = 1 ;ctr <= out1.Ncols(); ctr++){
	  if(!abscc.value()){
	  	out1.Column(ctr).Maximum1(tmp2);
      	out2(1,ctr) = tmp2;
	  }else
	  {
	  	out1.Column(ctr).MaximumAbsoluteValue1(tmp2);
      	out2(1,ctr) = tmp2;		
	  }
    }
      
	if(debug.value()){ 
      cerr << "DBG: out1 " << out1.Nrows() << " x " << out1.Ncols() << endl;
      cerr << "DBG: out2 " << out2.Nrows() << " x " << out2.Ncols() << endl;
	}
  }
  else{
	//no Target TCs
	if(verbose.value())
	  cout << " Calculating partial correlation maps "  << endl;

	out1 = zeros(seeds.Ncols(), data.Ncols());
	if(confounds.Storage()>0){
	  data = data - confounds * pinv(confounds) * data;
      seeds = seeds - confounds * pinv(confounds) * seeds;
    }

	if(debug.value()){ 
      cerr << "DBG: seeds " << seeds.Nrows() << " x " << seeds.Ncols() << endl;
      cerr << "DBG: data " << data.Nrows() << " x " << data.Ncols() << endl;
    }

	for(int ctr = 1 ;ctr <= seeds.Ncols(); ctr++){
	  Matrix tmp;
	  if(coords.Storage()>0){
	    tmp = orig_data.voxelts(coords(ctr,1), coords(ctr,2), coords(ctr,3));
		volume4D<float> tmpVol;
		tmpVol.setmatrix(out2,maskS);
		tmpVol( coords(ctr,1), coords(ctr,2), coords(ctr,3), 0) = ctr; 
		out2 = tmpVol.matrix(maskS); 
     	if(confounds.Storage()>0)
		  tmp = tmp - confounds * pinv(confounds) * tmp;
	  }
	  else{
		tmp = seeds.Column(ctr);
     	out2(1,ctr) = ctr;
	  }
	  for(int ctr2 =1; ctr2 <= data.Ncols(); ctr2++)
        out1(ctr,ctr2) = Melodic::corrcoef(tmp,data.Column(ctr2)).AsScalar();      
    }

	if(debug.value()){ 
      cerr << "DBG: out1 " << out1.Nrows() << " x " << out1.Ncols() << endl;
      cerr << "DBG: out2 " << out2.Nrows() << " x " << out2.Ncols() << endl;
    }
  }
}

void write_res(){
  if(verbose.value())
    cout << " Saving results " << endl;
		
  if(debug.value()) 
	cerr << "DBG: in write_res" << endl;
    
  if(regress_only.value()){
	save4D(out2,maskS, fnout.value()+"_index");
	save4D(out1,maskT, fnout.value()+"_corr");
  }
  else{
	save4D(out1,maskS, fnout.value()+"_corr");
	save4D(out2,maskS, fnout.value()+"_index");
  }

  if(out_ttcs.value() && ttcs.size()>0)
    for(int ctr = 0 ;ctr < (int)ttcs.size(); ctr++)
      write_ascii_matrix(ttcs.at(ctr),fnout.value()+"_ttc"+num2str(ctr+1)+".txt");	
  
  if(out_conf.value() && confounds.Storage()>0)
	write_ascii_matrix(confounds, fnout.value()+"_confounds.tx");
	
  if(out_seeds.value())   	
    save4D(seeds, maskS, fnout.value()+"_seeds");
  if(out_seedmask.value())   
    save4D(maskS,fnout.value()+"_seedmask");
  	
}

int do_work(int argc, char* argv[]) {
  if(setup())
	exit(1);
  calc_res();
  write_res();
  return 0;
}

////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[]){
	  Tracer tr("main");
	  OptionParser options(title, examples);
	  try{
	    // must include all wanted options here (the order determines how
	    //  the help message is printed)
		options.add(fnin);
		options.add(fnseed);
		options.add(fnout);
		options.add(fntarget);
		options.add(regress_only);
		options.add(fnconf);
		options.add(fnseeddata);
		options.add(map_bin);
		options.add(tc_mean);
		options.add(abscc);
		options.add(tc_order);
		options.add(out_seeds);
		options.add(out_seedmask);
		options.add(out_ttcs);
		options.add(out_conf);
		options.add(out_tcorr);
		options.add(verbose);
		options.add(help);
		options.add(debug);
	    options.parse_command_line(argc, argv);

	    // line below stops the program if the help was requested or 
	    //  a compulsory option was not set
	    if ( (help.value()) || (!options.check_compulsory_arguments(true)) ){
			options.usage();
			exit(EXIT_FAILURE);
	    }else{
	  		// Call the local functions
	  		return do_work(argc,argv);
		}
	  }catch(X_OptionError& e) {
		options.usage();
	  	cerr << endl << e.what() << endl;
	    exit(EXIT_FAILURE);
	  }catch(std::exception &e) {
	    cerr << e.what() << endl;
	  } 
}

