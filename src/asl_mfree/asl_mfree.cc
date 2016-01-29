/*   asl_mfree.cc 'Model-free' deconvolution for multi-TI ASL data

      Michael Chappell and Matthew Webster - IBME & FMIRB Image Analysis Group

      Copyright (C) 2011 University of Oxford */

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
#include <math.h>
#include <string>
#include "newmatap.h"
#include "newmatio.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "utils/tracer_plus.h"
#include "stdlib.h"

#include "readoptions.h"
#include "asl_mfree_functions.h"

using namespace Utilities;
using namespace NEWMAT;
using namespace NEWIMAGE;
using namespace MISCMATHS;

using namespace OXASL;

int main(int argc, char *argv[])
{
  try {

    cout << "ASL_MFREE (1.0)" << endl;

    //parse command line (puts these into the log file)
    ReadOptions& opts = ReadOptions::getInstance();
    opts.parse_command_line(argc,argv);

    cout << "Loading data" << endl;

    // load data    
    volume4D<float> data;
    read_volume4D(data,opts.datafile.value());

    // load mask
    volume<float> mask(data.xsize(),data.ysize(),data.zsize());
    mask.setdims(data.xdim(),data.ydim(),data.zdim());
    read_volume(mask,opts.maskfile.value());

    Matrix asldata;
    asldata = data.matrix(mask);
    //data.setmatrix(asldata,mask);
    int nvox=asldata.Ncols();
    
    //load aif
    volume4D<float> aif;
    read_volume4D(aif,opts.aif.value());

    
    // select AIF based on metric
    // load metric image if it exists
    volume<float> metric;
    if (opts.metric.set()) {
      cout << "Preparing AIFs" << endl;
      read_volume(metric,opts.metric.value());

      Prepare_AIF(aif, metric, mask, opts.mthresh.value());
      
      //volume4D<float> aifout;
      //aifout.setmatrix(aif,mask);
      save_volume4D(aif,opts.outname.value()+"_aifs");
    }

    Matrix aifmtx;
    aifmtx = aif.matrix(mask);

    // do deconvolution
    cout << "Performing deconvolution" << endl;
    ColumnVector mag;
    Matrix resid;
    Deconv(asldata,aifmtx,opts.dt.value(),mag,resid);

    // estimate BAT (of tissue)
    ColumnVector batt;
    if (opts.batout.set() | (opts.tcorrect.set() & !opts.batt.set())) {
      cout << "Estimating BAT" << endl;
      Estimate_onset(asldata,batt,opts.dt.value());

      if (opts.batout.set()) {
	//output the BAT image (from the tissue)
	volume4D<float> batout;
	batout.setmatrix(batt.AsMatrix(1,nvox),mask);
	save_volume4D(batout,opts.outname.value()+"_bat");
      }
    }

    // correct aif magntiude for differences in arrival time between aif and tissue 
    ColumnVector batd;
    if(opts.tcorrect.set()) {
      cout << "Performing timing correction" << endl;
      

      if (opts.bata.set()) {
	//calculate difference between tissue and AIF curves using suppled BAT images
	volume4D<float> bat_art;
	read_volume4D(bat_art,opts.bata.value());
	
	if (opts.batt.set()) {
	  //load supplied tissue BAT
	volume4D<float> bat_tiss;
	read_volume4D(bat_tiss,opts.batt.value());
	
	batt = (bat_tiss.matrix(mask)).AsColumn();
	}
	
	if (opts.metric.set()) {
	  // correct the AIF bat image to match the AIFs where a metric image has been supplied
	  Prepare_AIF(bat_art,metric,mask,opts.mthresh.value());
	}
	
	ColumnVector bata;
	bata = (bat_art.matrix(mask)).AsColumn();
	
	batd = batt-bata;
      }
      else {
	//otherwise estimate BAT difference using the peak in the residue function
	//Estimate_BAT_difference(resid,batd,opts.dt.value());

	// Estiamte BAT difference using edge detection
	ColumnVector bata;
	Estimate_onset(aifmtx,bata,opts.dt.value());

	batd = batt-bata;
      }

      for (int i=1; i<=batd.Nrows(); i++) { if (batd(i)<0.0) batd(i)=0.0; }
      Correct_magnitude(mag,batd,opts.T1.value(),opts.dt.value(),opts.fa.value());
    }

    if(opts.std.set()) {
      // do wild boostrapping std dev estimation for cbf
      cout << "Performing wild bootstrapping for precision estimation" << endl;
      ColumnVector magstd;
      BootStrap(aifmtx, asldata, opts.dt.value(), mag, resid, opts.nwb.value(), magstd);

      if (opts.tcorrect.set()) {
	// if needed we should correct the std dev for timing discrpancies between aif and ctc
	Correct_magnitude(magstd,batd,opts.T1.value(),opts.dt.value(),opts.fa.value());
      }

      // save it
      volume4D<float> stdoutVol; //stdout is a reserved name - can cause weird compile errors
      stdoutVol.setmatrix(magstd.AsMatrix(1,nvox),mask);
      save_volume4D(stdoutVol,opts.outname.value()+"_stddev");
    }
    
    cout << "Saving results" << endl;
    //output 
    volume4D<float> residout;
    residout.setmatrix(resid,mask);
    save_volume4D(residout,opts.outname.value()+"_residuals");

    volume4D<float> magout;
    magout.setmatrix(mag.AsMatrix(1,nvox),mask);
    save_volume4D(magout,opts.outname.value()+"_magntiude");

    cout << "ASL_MFREE - Done!" << endl;

  }

catch(Exception& e)
    {
      cerr << endl << e.what() << endl;
      return 1;
    }
  catch(X_OptionError& e)
    {
      cerr << endl << e.what() << endl;
      return 1;
    }

  cout << "Done." << endl;

  return 0;


}
