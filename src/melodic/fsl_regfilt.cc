/*  fsl_regfilt - 

    Christian F. Beckmann, FMRIB Analysis Group

    Copyright (C) 2006-2013 University of Oxford */
 
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

  string title=string("fsl_regfilt")+
		string("\nAuthor: Christian F. Beckmann \n Copyright(C) 2016-2013 University of Oxford\n")+
		string(" Data de-noising by regressing out part of a design matrix\n")+
		string(" using simple OLS regression on 4D images");
  string examples="fsl_regfilt -i <input> -d <design> -f <component numbers or filter threshold> -o <out> [options]";

//Command line Options {
  Option<string> fnin(string("-i,--in"), string(""),
		string("        input file name (4D image)"),
		true, requires_argument);
  Option<string> fnout(string("-o,--out"), string(""),
		string("output file name for the filtered data"),
		true, requires_argument);
  Option<string> fndesign(string("-d,--design"), string(""),
		string("file name of the matrix with time courses (e.g. GLM design or MELODIC mixing matrix)"),
		true, requires_argument);
  Option<string> fnmask(string("-m,--mask"), string(""),
		string("mask image file name"),
		false, requires_argument);		
   Option<string> filter(string("-f,--filter"),string(""),
		string("filter out part of the regression model, e.g. -f \"1,2,3\" "),
		false, requires_argument);
	Option<bool> freqfilt(string("-F,--freqfilt"),false,
		string("filter out components based on high vs. low frequency content "),
		false, no_argument);
	Option<bool> freq_ic(string("--freq_ic"),true,
		string("switch off IC Z-stats filtering as part of frequency filtering"),
		false, no_argument);		
	Option<float> freq_ic_smooth(string("--freq_ic_smooth"),5.0,
		string("smoothing width for IC Z-stats filtering as part of frequency filtering"),
		false, no_argument);		
	Option<float> freqthresh(string("--fthresh"),0.15,
		string("frequency threshold ratio - default: 0.15"),
		false,requires_argument);		
	Option<float> freqthresh2(string("--fthresh2"),0.02,
		string("frequency filter score threshold - default: 0.02"),
		false,requires_argument);		
	Option<bool> verbose(string("-v"),FALSE,
		string("        switch on diagnostic messages"),
		false, no_argument);
	Option<bool> aggressive(string("-a"),FALSE,
		string("        switch on aggressive filtering (full instead of partial regression)"),
		false, no_argument);
	Option<bool> perfvn(string("--vn"),FALSE,
		string("        perform variance-normalisation on data"),
		false, no_argument);
	Option<int> help(string("-h,--help"), 0,
		string("display this help text"),
		false,no_argument);
	Option<bool> debug(string("--debug"), false,
		string("switch on debug messages"),
		false,no_argument,false);
	// Output options	
	Option<string> outdata(string("--out_data"),string(""),
		string("output file name for pre-processed data (prior to denoising)"),
		false, requires_argument);
	Option<string> outmix(string("--out_mix"),string(""),
		string("output file name for new mixing matrix"),
		false, requires_argument);
	Option<string> outvnscales(string("--out_vnscales"),string(""),
		string("output file name for scaling factors from variance normalisation"),
		false, requires_argument);
		/*
}
*/
//Globals {
	int voxels = 0;
	float TR;
	Matrix data;
	Matrix design;
	Matrix fdesign;
	Matrix meanR, meanC;
	Matrix newData, newMix;
	RowVector vnscales;
	volume<float> mask;
	volume<float> Mean;
	vector<int> comps, ind;
	vector<int>::iterator it;
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
			tempVol.setTR(TR);
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

Matrix smooth_map(Matrix what, float howmuch){
	volume4D<float> tempVol;
	tempVol.setmatrix(what,mask);
    tempVol= smooth(tempVol,howmuch);
	Matrix out;
	out = tempVol.matrix(mask);
	return out;
}

int parse_filterstring(){
	int ctr=0;    
	char *p;
	char t[1024];
	const char *discard = ", [];{(})abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ~!@#$%^&*_-=+|\':><./?";

	strcpy(t, filter.value().c_str());
	p=strtok(t,discard);
	ctr = atoi(p);
	if(ctr>0 && ctr<=design.Ncols())
		comps.push_back(ctr);
	do{
	  	p=strtok(NULL,discard);
	  	if(p){
			ctr = atoi(p);
	    	if(ctr>0 && ctr<=design.Ncols())
				comps.push_back(ctr);
	    }
	}while(p);
	return 0;
}

int calc_freqindex(){
	if(debug.value()) cerr << " In calc_freqindex " << endl;
	
	fdesign = Melodic::calc_FFT(design);
	if(debug.value()) cerr << " fdesign: " << fdesign.Nrows() << " x " << fdesign.Ncols() << endl;	
	
	int Nps = fdesign.Nrows();
	float MAXf = 1/(2*TR);
	float Nthresh = ceil(Nps * freqthresh.value()/MAXf);	
	if(debug.value()) cerr << " Nps: " << Nps << " MAXf: " << MAXf << " Nthresh: " << Nthresh << endl;	
	
	Matrix sum_ratio;	
	sum_ratio = SP(sum(fdesign.Rows(1,Nthresh),1),pow(sum(sum(fdesign.Rows(Nthresh+1,Nps))),-1));
	sum_ratio /= (float)sum_ratio.MaximumAbsoluteValue();
	if(debug.value()) cerr << " sum_ratio: " << sum_ratio << endl;		

	if(freq_ic.value()){		
		Matrix scores = zeros(1,design.Ncols());
		{
			Matrix ICs, noisestddev, stdNoisei,unmixMatrix;
	
			unmixMatrix = pinv(design);
			ICs = unmixMatrix * data;
			noisestddev = stdev(data-design*ICs);	
			stdNoisei = pow(noisestddev*
				std::sqrt((float)(data.Nrows()-1))/
				std::sqrt((float)(data.Nrows()-ICs.Nrows())),-1);
			
			ColumnVector diagvals;
			diagvals = pow(diag( unmixMatrix*unmixMatrix.t()),-0.5);
			
			ICs=smooth_map(SP(ICs,diagvals*stdNoisei),freq_ic_smooth.value());
			ICs= SP(ICs,ones(ICs.Nrows(),1)*meanR);
			volume4D<float> tempVol;
			tempVol.setmatrix(ICs,mask);
			tempVol.threshold(0.0);
			for(int ctr = 0; ctr < design.Ncols(); ctr++ )
				scores(1,ctr+1) = tempVol[ctr].percentile(0.99,mask);
			scores/=scores.MaximumAbsoluteValue();
			scores-=scores.MinimumAbsoluteValue();
			if(debug.value()) cerr << " initial scores: " << scores << endl;
		}
		scores = SP(scores,sum_ratio);
		scores /= scores.Maximum();
		if(debug.value()) cerr << " scores: " << scores << endl;
		for(int ctr = 1; ctr <= design.Ncols(); ctr++ )
			if(scores(1,ctr) < freqthresh2.value())	
			comps.push_back(ctr);
	}	
	return 0;
}

int get_comp(){
	if(filter.value().length()>0 && parse_filterstring())
		return 1;

	if(freqfilt.value() && calc_freqindex())
		return 1;
		
	//sort and remove duplicates 
	sort (comps.begin(), comps.end());
	it = unique (comps.begin(), comps.end());
	comps.resize( it - comps.begin() );
	
	if(debug.value()){
		for (it=comps.begin(); it!=comps.end(); ++it)
	    	cout << " " << *it;
	  	cout << endl; 
	}
	
	return 0;
}

int dofilter(){

	if(verbose.value())
		cout << "  Calculating maps " << endl;  

  	Matrix unmixMatrix = pinv(design);
  	Matrix maps = unmixMatrix * data;

	Matrix noisedes;
 	Matrix noisemaps;

   	noisedes = design.Column(comps.at(0));
    noisemaps  = maps.Row(comps.at(0)).t();    
 	
  	for(int ctr = 1; ctr < (int)comps.size();++ctr){
    	noisedes  |= design.Column(comps.at(ctr));
  		noisemaps |= maps.Row(comps.at(ctr)).t();
	}
	if(debug.value()) cerr << " noisedes " << noisedes.Nrows() << " x " << noisedes.Ncols() << endl;	

	if(verbose.value())
		cout << "  Calculating filtered data " << endl;

  	if(aggressive.value())
		newData = data - noisedes * (pinv(noisedes)*data);
	else
		newData = data - noisedes * noisemaps.t();
		
	if(perfvn.value())
		newData = SP(newData,ones(newData.Nrows(),1)*vnscales);	
  	newData = newData + ones(newData.Nrows(),1)*meanR;

	for(int ctr = 1; ctr <= design.Ncols();++ctr)
		ind.push_back(ctr);
	for(int ctr = 0; ctr < (int)comps.size();++ctr)
		it=remove(ind.begin(),ind.end(),comps.at(ctr));
	ind.resize(design.Ncols()-comps.size());
	if(debug.value()){
		for (it=ind.begin(); it!=ind.end(); ++it)
	    	cout << " " << *it;
	  	cout << endl; 
	}	
	
	if(ind.size()>0){
		newMix=design.Column(ind.at(0));
		for(int ctr = 1; ctr < (int)ind.size();++ctr)
			newMix |= design.Column(ind.at(ctr));
		newMix = newMix - noisedes * (pinv(noisedes)*newMix);
		
		if(debug.value())
			cerr << " newMix " << newMix.Nrows() << " x " << newMix.Ncols() << endl;
	}
		
  	return 0;	
}

int setup(){
	if(fsl_imageexists(fnin.value())){//read data
		//input is 3D/4D vol
		volume4D<float> tmpdata;
		read_volume4D(tmpdata,fnin.value());
		TR=tmpdata.TR();
		
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
	
	design = read_ascii_matrix(fndesign.value());

	if(!isimage(data)){
		cerr << "ERROR: need to specify 4D input to use filtering" << endl;
		return 1;
	}
	
	meanR=mean(data,1);
	data = remmean(data,1);
	meanC=mean(design,1);
	design = remmean(design,1);
	if(perfvn.value())
		vnscales = Melodic::varnorm(data);
		
	if(debug.value()) cerr << " data: " << data.Nrows() << " x " << data.Ncols() << endl;	
    if(debug.value()) cerr << " design: " << design.Nrows() << " x " << design.Ncols() << endl;	
				
	return 0;	
}

void write_res(){	
	saveit(newData,fnout.value());
	if(outdata.value()>"")
		saveit(data,outdata.value());
	if(outvnscales.value()>"")
		saveit(vnscales,outvnscales.value());
	if(outmix.value()>"" && newMix.Storage()>0)
		saveit(newMix,outmix.value());
}

int do_work(int argc, char* argv[]) {
  	if(setup())
		exit(1);
	if(get_comp())
		exit(1);
	if(dofilter())
		exit(1);	
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
			options.add(fnout);
			options.add(fndesign);
			options.add(fnmask);
			options.add(filter);
			options.add(freqfilt);
			options.add(freq_ic);
			options.add(freq_ic_smooth);
			options.add(freqthresh);
			options.add(freqthresh2);
			options.add(perfvn);
			options.add(verbose);
			options.add(aggressive);
			options.add(help);
			options.add(debug);
			options.add(outdata);
			options.add(outmix);
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
		}catch(X_OptionError& e) {
			options.usage();
	  	cerr << endl << e.what() << endl;
	    exit(EXIT_FAILURE);
	  }catch(std::exception &e) {
	    cerr << e.what() << endl;
	  } 
}

