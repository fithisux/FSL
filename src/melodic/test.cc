/*  test.cc
  
    Christian F. Beckmann, FMRIB Analysis Group
  
    Copyright (C) 1999-20013 University of Oxford */

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
#include <ctime>
#include "newimage/newimageall.h"
#include "melhlprfns.h"
#include <iostream>

#ifdef __APPLE__
#include <mach/mach.h>
#define memmsg(msg) { \
  struct task_basic_info t_info; \
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT; \
  if (KERN_SUCCESS == task_info(mach_task_self(), TASK_BASIC_INFO, (task_info_t) &t_info, &t_info_count)) \
	{ \
		cout << msg << " res: " << t_info.resident_size/1000000 << " virt: " << t_info.virtual_size/1000000 << "\n"; \
		cout.flush(); \
	} \
}
#else
#define memmsg(msg) { \
   cout << msg; \
}
#endif

// a simple message macro that takes care of cout and log
#define message(msg) { \
  cout << msg; \
  cout.flush(); \
}

#define outMsize(msg,Mat) { \
    cerr << "     " << msg << "  " <<Mat.Nrows() << " x " << Mat.Ncols() << endl;	\
}


using namespace MISCPLOT;
using namespace MISCMATHS;
using namespace Utilities;
using namespace std;
using namespace Melodic;

// GLOBALS

clock_t tictime;


// The two strings below specify the title and example usage that is
// printed out as the help or usage message

  string title=string("fsl_BLAH")+
		string("\nAuthor: Christian F. Beckmann \nCopyright(c) 2008-2013 University of Oxford\n")+
		string(" \n \n")+
		string(" \n");
  string examples="fsl_BLAH [options]";

//Command line Options {
    Option<string> fnin(string("-i,--in"), string(""),
		string("input file name (matrix 3D or 4D image)"),
		false, requires_argument);
    Option<string> fnmask(string("-m"), string(""),
			string("mask file name "),
			false, requires_argument);
	Option<int> help(string("-h,--help"), 0,
		string("display this help text"),
		false,no_argument);
	Option<int> xdim(string("-x,--xdim"), 0,
		string("xdim"),
		false,requires_argument);
	Option<int> ydim(string("-y,--ydim"), 0,
		string("ydim"),
		false,requires_argument);
	Option<int> econ(string("-e,--econ"), 0,
		string("econ: how to liump stuff"),
		false,requires_argument);
		/*
}
*/
////////////////////////////////////////////////////////////////////////////

// Local functions

void tic(){
	tictime = clock();
}

void toc(){
	cerr << endl << "TOC: " << float(clock()-tictime)/CLOCKS_PER_SEC << " seconds" << endl<<endl;
}

Matrix calccorr(const Matrix& in, int econ)
  { 
    Matrix Res;
	int nrows=in.Nrows();
	int ncols=in.Ncols();    
    Res = zeros(nrows,nrows);

    if(econ>0){
      RowVector colmeans(ncols);
	  for (int n=1; n<=ncols; n++) {
        colmeans(n)=0;
        for (int m=1; m<=nrows; m++) {
          colmeans(n)+=in(m,n);
        }
        colmeans(n)/=nrows;
      }
      int dcol = econ;
	  Matrix suba; 

      for(int ctr=1; ctr <= in.Ncols(); ctr+=dcol){
	    suba=in.SubMatrix(1,nrows,ctr,Min(ctr+dcol-1,ncols));
		int scolmax = suba.Ncols();

		for (int n=1; n<=scolmax; n++) {
	        double cmean=colmeans(ctr + n - 1);
	        for (int m=1; m<=nrows; m++) {
	          suba(m,n)-=cmean;
	        }
	    }
		
	    Res += suba*suba.t() / ncols;
      }
    }
    else
      Res = cov(in.t());
    return Res;
  }  //Matrix calccorr

int do_work(int argc, char* argv[]) {

	tic();
	Matrix MatrixData;
	volume<float> Mean;
  	
	if(xdim.value()==0 && ydim.value()==0)
    {
		volume4D<float> RawData;
		volume<float> theMask;
		toc();
    	//read data
    	message("Reading data file " << (string)fnin.value() << "  ... ");
    	read_volume4D(RawData,fnin.value());
    	message(" done" << endl);
 
		Mean = meanvol(RawData);
		toc();
    	message("Reading mask file " << (string)fnmask.value() << "  ... ");
      	read_volume(theMask,fnmask.value());

		memmsg(" Before reshape ");
    	MatrixData = RawData.matrix(theMask);
    }
	else{
		Matrix data = unifrnd(xdim.value(),ydim.value());
	    outMsize("data", data);
		tic(); calccorr(data,econ.value()); toc();
		
		data = unifrnd(10,100);
		outMsize("data", data);
		tic(); calccorr(data,econ.value()); toc();
		
		data = unifrnd(100,1000);
		outMsize("data", data);
		tic(); calccorr(data,econ.value()); toc();
		
		data = unifrnd(100,10000);
		outMsize("data", data);
		tic(); calccorr(data,econ.value()); toc();

		data = unifrnd(300,200000);
		outMsize("data", data);
		tic(); calccorr(data,econ.value()); toc();		
		
		data = unifrnd(500,20000);
		outMsize("data", data);
		tic(); calccorr(data,econ.value()); toc();
		
		data = unifrnd(500,200000);
		outMsize("data", data);
		tic(); calccorr(data,econ.value()); toc();
		
		
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
			options.add(fnmask);		
			options.add(help);
			options.add(xdim);
			options.add(ydim);
			options.add(econ);
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

