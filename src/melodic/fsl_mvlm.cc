/*  fsl_glm - 

    Christian F. Beckmann, FMRIB Analysis Group

    Copyright (C) 2006-2013 University of Oxford  */

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


//Header & includes
#include "libvis/miscplot.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include "utils/options.h"
#include <vector>
#include "newimage/newimageall.h"
#include "melhlprfns.h"

#define message(msg) { \
  if(verbose.value()) \
    { \
      cout << msg; \
    } \
}

#define dbgmsg(msg) { \
  if(debug.value()) {\
	cerr << msg << endl; } \
}

#define outMsize(msg,Mat) { \
  if(debug.value())						\
    cerr << "     " << msg << "  " <<Mat.Nrows() << " x " << Mat.Ncols() << endl;	\
}

#define outM(msg,Mat) { \
  if(verbose.value())						\
    cout << "     " << msg << "  " <<Mat.Nrows() << " x " << Mat.Ncols() << endl;	\
}


using namespace MISCPLOT;
using namespace MISCMATHS;
using namespace Utilities;
using namespace std;


//Command-line output
// The two strings below specify the title and example usage that is
// printed out as the help or usage message

  string title=string("fsl_mvlm (Version 1.0)")+
		string("\nAuthor: Christian F. Beckmann \nCopyright(C) 2006-2013 University of Oxford\n")+
		string(" \n Multivariate Linear Model regression on\n")+
		string(" time courses and/or 3D/4D imges using  SVD (PCA), PLS, normalised PLS, \n")+
		string(" CVA, SVD-CVA or  MLM\n\n");
  string examples="fsl_mvlm -i <input> -o <output> [options]";


//Command line Options 
    Option<string> fnin(string("-i,--in"), string(""),
		string("        input file name (text matrix or 3D/4D image file)"),
		true, requires_argument);
    Option<string> fnout(string("-o,--out"), string(""),
		string("basename for output files "),
		true, requires_argument);
	Option<string> approach(string("-a,--alg"), string("PCA"),
		string("algorithm for decomposition: PCA (or SVD; default), PLS, orthoPLS, CVA, SVD-CVA, MLM, NMF"),
		false, requires_argument);	
    Option<string> fndesign(string("-d,--design"), string(""),
		string("file name of the GLM design matrix (time courses or spatial maps)"),
		false, requires_argument);
    Option<string> fnmask(string("-m,--mask"), string(""),
		string("mask image file name if input is image"),
		false, requires_argument);
	Option<bool> normdes(string("--des_norm"),FALSE,
		string("switch on normalisation of the design matrix columns to unit std. deviation"),
		false, no_argument);
	Option<bool> perfvn(string("--vn"),FALSE,
		string("        perform MELODIC variance-normalisation on data"),
		false, no_argument);
	Option<bool> perf_demean(string("--demean"),FALSE,
		string("switch on de-meaning of design and data"),
		false, no_argument);
	Option<int> nmfdim(string("--nmf_dim"), 0,
		string(" Number of underlying factors for NMF"),
		false,requires_argument);
	Option<int> nmfitt(string("--nmfitt"), 100,
		string("number of NMF itterations (default 100)"),
		false,requires_argument);
	Option<int> help(string("-h,--help"), 0,
		string("display this help text"),
		false,no_argument);
	Option<bool> verbose(string("-v,--verbose"),FALSE,
		string("switch on verbose output"),
		false, no_argument);
	Option<bool> debug(string("--debug"),FALSE,
		string("switch on debug output"),
		false, no_argument, false);
	// Output options	
	Option<string> outres(string("--out_res"),string(""),
		string("output file name for residuals"),
		false, requires_argument, false);
	Option<string> outdata(string("--out_data"),string(""),
		string("output file name for pre-processed data"),
		false, requires_argument);
	Option<string> outvnscales(string("--out_vnscales"),string(""),
		string("output file name for scaling factors for variance normalisation"),
		false, requires_argument);

//Globals 
	Melodic::basicGLM glm;
	int voxels = 0;
	Matrix data, tmpdata;
	Matrix design;
	Matrix meanR;
	Matrix svd_X_U, svd_X_V, svd_Y_U, svd_Y_V;
	DiagonalMatrix svd_X_D, svd_Y_D;	
	RowVector vnscales;
	volume<float> mask;  

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
	else if(fsl_imageexists(fndesign.value()))
		write_ascii_matrix(what.t(),fname);
	else
		write_ascii_matrix(what,fname);
}

int setup(){
	
	dbgmsg(" In <setup>");
	message(" Reading data " << fnin.value() << " ... ");
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
			mask = tmpdata[0]*0.0+1.0;	
		}
		
		data = tmpdata.matrix(mask);
		voxels = data.Ncols();
		if(perf_demean.value())
			data = remmean(data,1);
		if(perfvn.value())
			vnscales = Melodic::varnorm(data);
	}
	else
		data = read_ascii_matrix(fnin.value());	

	message("done" << endl;);
	if(fndesign.value().length()>0){
		message(" Reading design " << fndesign.value() << " ... ");
		if(fsl_imageexists(fndesign.value())){//read design
			volume4D<float> tmpdata;
			read_volume4D(tmpdata,fndesign.value());
			if(!samesize(tmpdata[0],mask)){
				cerr << "ERROR: GLM design does not match input image in size" << endl;
				return 1;
			}
			design = tmpdata.matrix(mask).t();
			data = data.t();
		}else{
			design = read_ascii_matrix(fndesign.value());
		}
		message("done" << endl;);
	}else
		design = ones(data.Nrows(),1);

	meanR=mean(data,1);
	if(perf_demean.value()){
		data = remmean(data,1);
		design = remmean(design,1);
	}
	
	if(normdes.value())
		design =  SP(design,ones(design.Nrows(),1)*pow(stdev(design,1),-1));
	
	SVD( design, svd_X_D, svd_X_U, svd_X_V );
	if(approach.value()!=string("NMF")){
		if(data.Nrows()>=data.Ncols())
			SVD ( data, svd_Y_D, svd_Y_U, svd_Y_V );
		else{
			SVD ( data.t(), svd_Y_D, svd_Y_V, svd_Y_U );
		}
	}		

	if(fnout.value().length() == 0){
		string basename = fnin.value();
		basename = make_basename(basename);
		fnout.set_T(basename+string("_mvlm_"));	
	}
	
	outM("Data matrix : ", data);
	outM("Design matrix : ", design);

	dbgmsg(" initial SVD : ");
	outMsize("svd_Y_U",svd_Y_U);
	outMsize("svd_Y_V",svd_Y_V);
	outMsize("svd_Y_D",svd_Y_D);

	outMsize("svd_X_U",svd_X_U);
	outMsize("svd_X_V",svd_X_V);
	outMsize("svd_X_D",svd_X_D);

	dbgmsg(" Leaving <setup>");
	return 0;	
}

void write_res(){
	dbgmsg(" In <write_res>");
	
	message(" Writing results ... ")
	if(isimage(svd_Y_V)){
		saveit(svd_Y_V,fnout.value()+string("maps"));
		saveit(svd_Y_U,fnout.value()+string("tcs"));
	}
	else{
		saveit(svd_Y_V.t(),fnout.value()+string("tcs"));
		saveit(svd_Y_U,fnout.value()+string("maps"));
	}
	saveit(svd_Y_D.AsColumn(),fnout.value()+string("scales"));
	if(outres.value()>"")
		
	if(outdata.value()>"")
		saveit(data,outdata.value());
	if(outvnscales.value()>"")
		saveit(vnscales,outvnscales.value());
	
	message("done" << endl;);
	dbgmsg(" Leaving <write_res>");
		
}

int do_work(int argc, char* argv[]) {
	dbgmsg(" In <do_work>");
	
    if(setup())
      exit(1);

	//modify data
	
	//X = svd_X_U * svd_X_D * svd_X_V.t();
	//Y = svd_Y_U * svd_Y_D * svd_Y_V.t();
	//X'X = svd_X_V *pow(svd_X_D,2) * svd_X_V.t();
	//(X'X)^(-1) = svd_X_V *pow(svd_X_D,-2) * svd_X_V.t()
 	//(X'X)^(-1/2) = svd_X_V *pow(svd_X_D,-1) * svd_X_V.t()
			
	if(approach.value()==string("PLS")) {
		message(" Using method : " << approach.value() << endl;);
		data = design.t() * data;
	} 
	if(approach.value()==string("orthoPLS")) {
		message(" Using method : " << approach.value() << endl;);
		data = (svd_X_V * svd_X_D.i() * svd_X_V.t()) * design.t() * data;
	} 
	if(approach.value()==string("CVA")) {
		message(" Using method : " << approach.value() << endl;);
		data = design.t() * svd_Y_U * svd_Y_V.t();
		data = (svd_X_V * svd_X_D.i() * svd_X_V.t() ) * data;
	}
	if(approach.value()==string("SVD-CVA")) {
		message(" Using method : " << approach.value() << endl;);
		tmpdata = data;
		data = design.t() * svd_Y_U;
		data = (svd_X_V * svd_X_D.i() * svd_X_V.t() ) * data;
	}
	if(approach.value()==string("MLM")) {
		message(" Using method : " << approach.value() << endl;);

		Matrix RE;
	    DiagonalMatrix RD;
	    SymmetricMatrix tmp;
	    tmp << cov(data.t());
	    EigenValues(tmp,RD,RE);
//		S = RE * RD * RE.t()
		
		tmp << sqrtm(svd_X_V * svd_X_D * svd_X_U.t() * RE * RD * RE.t() *svd_X_U * svd_X_D * svd_X_V.t());
		data =  tmp.i()*design.t() * data;
	}

    if( approach.value()!=string("MLM") && approach.value()!=string("CVA") && approach.value()!=string("PLS") &&
		approach.value()!=string("SVD-CVA") && approach.value()!=string("orthoPLS") && approach.value()!=string("NMF"))
		message(" Using method : PCA" << endl;);
	
	//perform an SVD on data	
	outMsize(" New Data ", data);
	
	if(approach.value()!=string("NMF")){	
		if(data.Nrows()>=data.Ncols())
			SVD ( data, svd_Y_D, svd_Y_U, svd_Y_V );
		else{
			SVD ( data.t(), svd_Y_D, svd_Y_V, svd_Y_U );
			svd_Y_U = svd_Y_U.t();
			svd_Y_V = svd_Y_V.t();		
		}	

		dbgmsg(" Finished SVD : ");
		outMsize("svd_Y_U",svd_Y_U);
		outMsize("svd_Y_V",svd_Y_V);
		outMsize("svd_Y_D",svd_Y_D);
	
		svd_Y_V = sqrtm(svd_Y_D) * svd_Y_V;	
		svd_Y_U = svd_Y_U * sqrtm(svd_Y_D);	

		if(approach.value()==string("SVD-CVA"))
			svd_Y_V = svd_Y_V *tmpdata;
	}
	else{ //NMF
		float err, err_old;
		Matrix Ratio, Diff;

		if(nmfdim.value()==0)
			nmfdim.set_T(data.Nrows());
		
		message("Using "<< nmfdim.value() << " dimensions" << endl;); 	
		svd_Y_U =  unifrnd(data.Nrows(), nmfdim.value());
		svd_Y_V =  unifrnd(nmfdim.value(), data.Ncols());
		// re-scale columns of svd_Y_U to unit amplitude
		Ratio = pow(stdev(svd_Y_U),-1);
		svd_Y_U = SP(svd_Y_U,ones(svd_Y_U.Nrows(),1)*Ratio);
				
		Diff = data - svd_Y_U * svd_Y_V;
		err = Diff.SumAbsoluteValue()/(data.Ncols()*data.Nrows());
			
		for(int k=1; k< nmfitt.value(); k++)
		{
		//	Ratio = SP(data,pow(svd_Y_U * svd_Y_V,-1));
		//	svd_Y_U  =  SP(svd_Y_U, Ratio * svd_Y_V.t());
		//	svd_Y_U  =  SP(svd_Y_U, pow( ones(svd_Y_U.Nrows(),1) * sum(svd_Y_U,1),-1));
		//	svd_Y_V  =  SP(svd_Y_V, svd_Y_U.t()* Ratio);
		//	
		// Lee & Seung multiplicatice updates
		
			Ratio = SP(svd_Y_U.t() * data, pow(svd_Y_U.t() * svd_Y_U * svd_Y_V ,-1));
			svd_Y_V  = SP(svd_Y_V,Ratio);

			Ratio = SP(data * svd_Y_V.t(),pow(svd_Y_U * (svd_Y_V * svd_Y_V.t()),-1));
			svd_Y_U  = SP(svd_Y_U,Ratio);

			// re-scale columns of svd_Y_U to unit amplitude
			Ratio = pow(stdev(svd_Y_U),-1);
			svd_Y_U = SP(svd_Y_U,ones(svd_Y_U.Nrows(),1)*Ratio);
			
			Diff = data - svd_Y_U * svd_Y_V;
			err_old = err;
			err = Diff.SumSquare()/(data.Ncols()*data.Nrows());			
			message(" Error " << err << endl;);
		}	
	}
	
	write_res();
	
	dbgmsg(" Leaving <do_work>");
	
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
			options.add(fnout);
			options.add(approach);
			options.add(fndesign);
			options.add(fnmask);
			options.add(normdes);
			options.add(perfvn);
			options.add(perf_demean);
			options.add(nmfdim);
			options.add(nmfitt);
			options.add(help);
			options.add(verbose);
			options.add(debug);
			options.add(outres);
			options.add(outdata);
			options.add(outvnscales);

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
      }
      catch(X_OptionError& e){
		  options.usage();
	      cerr << endl << e.what() << endl;
	      exit(EXIT_FAILURE);
	  }
	  catch(std::exception &e){
	    cerr << e.what() << endl;
	  } 
}

