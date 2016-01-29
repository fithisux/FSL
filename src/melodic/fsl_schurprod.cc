/*  fsl_glm - 

    Christian F. Beckmann, FMRIB Analysis Group

    Copyright (C) 2008-2013 University of Oxford  */

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

  string title=string("fsl_schurprod (Version 1.0)")+
		string("\nAuthor: Christian F. Beckmann \nCopyright(C) University of Oxford\n")+
		string(" \n Generates element-wise matrix products or product of matrices against vectors from 4D data\n")+
		string(" ");
  string examples="fsl_schurprod -i <input> -d <time series> -o <basename> [options]";

//Command line Options {
  Option<string> fnin(string("-i,--in"), string(""),
		string("        input file name (4D image file)"),
		true, requires_argument);
  Option<string> fnout(string("-o,--out"), string(""),
		string("output file base name"),
		true, requires_argument);
  Option<string> fndes(string("-d,--design"), string(""),
		string("ASCII text matrix of time series to be correlated"),
		true, requires_argument);
  Option<int> indx(string("-i,--index"), 0,
		string("index of column in the design to be used for matrix product calculation"),
		false, requires_argument);
  Option<string> fnmask(string("-m,--mask"), string(""),
		string("mask image file name"),
		false, requires_argument);
  Option<bool> regress_only(string("-r,--regression"), TRUE,
		string("use regression rather than correlation"),
		false, no_argument);
  Option<bool> verbose(string("-v,--verbose"), false,
		string("switch on diagnostic messages"),
		false, no_argument);
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
	int voxels = 0;
	Matrix data;
	Matrix design;
	Matrix meanR, meanC;
	Matrix newData;
	volume<float> mask;
	volume<float> Mean;
	
	 /*
}
*/
////////////////////////////////////////////////////////////////////////////

// Local functions
void save4D(Matrix what, string fname){
		if(what.Ncols()==data.Ncols()||what.Nrows()==data.Nrows()){
			volume4D<float> tempVol;
			if(what.Nrows()>what.Ncols())
				tempVol.setmatrix(what.t(),mask);
			else
				tempVol.setmatrix(what,mask);
			save_volume4D(tempVol,fname);
		}
}

bool isimage(Matrix what){
	if((voxels > 0)&&(what.Ncols()==voxels || what.Nrows()==voxels))
		return TRUE;
	else
		return FALSE;
}

void saveit(Matrix what, string fname){
	if(isimage(what))
		save4D(what,fname);
	else
		write_ascii_matrix(what,fname);
}

int setup(){
  	if(fsl_imageexists(fnin.value())){//read data
		//input is 3D/4D vol
		volume4D<float> tmpdata;
		read_volume4D(tmpdata,fnin.value());
		
		// create mask
		if(fnmask.value()>""){
			read_volume(mask,fnmask.value());
			if(!samesize(tmpdata[0],mask)){
				cerr << "ERROR: Mask image does not match input image" << endl;
				return 1;
			};
		}else{
			if(verbose.value())
				cout << "  Creating mask image "  << endl;
			Mean = meanvol(tmpdata);
			float Mmin, Mmax;
			Mmin = Mean.min(); Mmax = Mean.max();
			mask = binarise(Mean,float(Mmin + 0.01* (Mmax-Mmin)),Mmax);
		}
		
		data = tmpdata.matrix(mask);
		voxels = data.Ncols();
		if(verbose.value())
			cout << " Data matrix size : " << data.Nrows() << " x " << voxels << endl;
		
	}else{
		cerr << "ERROR: cannot read input image " << fnin.value()<<endl;
		return 1;
	}
	
	design = read_ascii_matrix(fndes.value());

	if(!isimage(data)){
		cerr << "ERROR: need to specify 4D input file" << endl;
		return 1;
	}
	
	if(indx.value()>design.Ncols()){
		cerr << "ERROR: index value too high" << endl;
		return 1;
	}
	
	meanR=mean(data,1);
	data = remmean(data,1);
	meanC=mean(design,1);
	design = remmean(design,1);
	if(regress_only.value()){
		Matrix tmpscales;
		tmpscales = stdev(data);
		data=SP(data,ones(data.Nrows(),1)*pow(tmpscales,-1));
		tmpscales = stdev(design);
		design=SP(data,ones(design.Nrows(),1)*pow(tmpscales,-1));
	}

    if(indx.value()>0)
		design = design.Columns(indx.value(),indx.value());


	if(debug.value()) cerr << " data: " << data.Nrows() << " x " << data.Ncols() << endl;	
    if(debug.value()) cerr << " design: " << design.Nrows() << " x " << design.Ncols() << endl;	

  return 0;	
}


void calc_res(int id=0){
  if(debug.value()) 
	cerr << "DBG: in calc_res" << endl;

	newData = SP(data,design.Column(id)*ones(1,data.Ncols()));
}

void write_res(int id=0){
  if(verbose.value())
    cout << " Saving results " << endl;
		
  if(debug.value()) 
	cerr << "DBG: in write_res" << endl;  	

  if(indx.value()>0)
  	saveit(newData,fnout.value()+num2str(id));
  else
	saveit(newData,fnout.value());
}

int do_work(int argc, char* argv[]) {
  if(setup())
	exit(1);
  for (int ctr=1; ctr<= design.Ncols(); ctr++){
  	calc_res(ctr);
  	write_res(ctr);
  }	
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
		options.add(fndes);
		options.add(fnout);
		options.add(regress_only);
		options.add(indx);
		options.add(fnmask);
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

