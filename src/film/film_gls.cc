/*  film_gls.cc

    Mark Woolrich and Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 1999-2008 University of Oxford  */

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
#define WANT_STREAM
#define WANT_MATH

#include "newimage/newimageall.h"
#include "fslsurface/fslsurface.h"
#include "utils/log.h"
#include "AutoCorrEstimator.h"
#include "paradigm.h"
#include "FilmGlsOptions.h"
#include "glimGls.h"
#include <vector>
#include <string>

using namespace NEWMAT;
using namespace FILM;
using namespace Utilities;
using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace fslsurface_name;


int main(int argc, char *argv[])
{
  try{
    
    Log& logger = LogSingleton::getInstance();
    FilmGlsOptions& globalopts = FilmGlsOptions::getInstance();
    globalopts.parse_command_line(argc, argv, logger);
    Matrix tContrasts,fContrasts;
    if (globalopts.contrastFile.value() != "") {
      tContrasts=read_vest(globalopts.contrastFile.value());
    }
    if (globalopts.fContrastFile.value() != "") {
      fContrasts=read_vest(globalopts.fContrastFile.value());
    }
    // load data
    string epifname("epivolume");
    volume4D<float> input_data, reference;
    volume<float> mask, variance;
    fslSurface<float, unsigned int> surfaceData;
    Matrix datam;
    vector<long> labels;


    if ( globalopts.analysisMode.value()=="surface" ) {
      read_surface(surfaceData,globalopts.inputDataName.value());
      datam.ReSize(surfaceData.getNumberOfScalarData(),surfaceData.getNumberOfVertices());
      for(unsigned int vertex=0; vertex < surfaceData.getNumberOfVertices(); vertex++) {
	for(unsigned int timepoint=0; timepoint < surfaceData.getNumberOfScalarData(); timepoint++) 
	   datam(timepoint+1,vertex+1)=surfaceData.getScalar(timepoint,vertex);
	datam.Column(vertex+1)-=(datam.Column(vertex+1).Sum()/datam.Nrows());
      }
    } else {
      read_volume4D(input_data,globalopts.inputDataName.value());
      reference=input_data[int(input_data.tsize()/2)-1];
      copybasicproperties(input_data,reference);
      mask=meanvol(input_data);
      variance=variancevol(input_data);
      input_data-=mask;
      mask.binarise(globalopts.thresh.value(),mask.max()+1,exclusive);
      variance.binarise(1e-10,variance.max()+1,exclusive); //variance mask needed if thresh is -ve to remove background voxels (0 variance)
      mask*=variance; //convolved mask ensures that only super-threshold non-background voxels pass
      datam=input_data.matrix(mask,labels);
    }
    int sizeTS(datam.Nrows()), numTS(datam.Ncols());

    // Load paradigm:
    Paradigm paradigm;
    if(!globalopts.ac_only.value())
      paradigm.load(globalopts.paradigmfname.value(), "", "", false, sizeTS);
    else
      paradigm.setDesignMatrix(sizeTS); // set design matrix to be one ev with all ones:

    if(globalopts.verbose.value())
      write_vest(logger.appendDir("Gc"), paradigm.getDesignMatrix());

    if (globalopts.voxelwise_ev_numbers.value().size()>0 && globalopts.voxelwiseEvFilenames.value().size()>0)
      paradigm.loadVoxelwise(globalopts.voxelwise_ev_numbers.value(),globalopts.voxelwiseEvFilenames.value(),mask);
    

    OUT(paradigm.getDesignMatrix().Nrows());
    OUT(paradigm.getDesignMatrix().Ncols());
    OUT(sizeTS);
    OUT(numTS);

    // Setup GLM:
    int numParams = paradigm.getDesignMatrix().Ncols();
    GlimGls glimGls(numTS, sizeTS, numParams, tContrasts.Nrows(),fContrasts.Nrows());

    // Residuals container:
    Matrix residuals(sizeTS, numTS);

    // Setup autocorrelation estimator:
    AutoCorrEstimator acEst(residuals);

    acEst.mask=mask;

    if(!globalopts.noest.value()) {
	cout << "Calculating residuals..." << endl; 
	for(int i = 1; i <= numTS; i++)
	  {						    
            glimGls.setData(datam.Column(i), paradigm.getDesignMatrix(i,mask,labels), i, tContrasts, fContrasts);
	    residuals.Column(i)=glimGls.getResiduals();
	  }
	cout << "Completed" << endl; 
	
	cout << "Estimating residual autocorrelation..." << endl; 
		
	if( globalopts.fitAutoRegressiveModel.value() ) {
          glimGls.saveData(logger.getDir() + "/betas",acEst.fitAutoRegressiveModel(),input_data,mask,true,true,reference.tdim(),false,-1,surfaceData,globalopts.analysisMode.value());
	}
	else if( globalopts.multitapersize.set() ) {
	  acEst.calcRaw();
	  acEst.multitaper(int(globalopts.multitapersize.value()));
	}
	else if(globalopts.pava.value()) {
	  acEst.calcRaw();
	  if(globalopts.smoothACEst.value()) {
	    if ( globalopts.analysisMode.value()=="surface" ) {
	      fslSurface<float, unsigned int> topologyData;
	      read_surface(topologyData,globalopts.inputDataName2.value());
	      acEst.spatiallySmooth(topologyData,globalopts.epith.value(),globalopts.ms.value());
	    } else
	      acEst.spatiallySmooth(reference.matrix(mask).t(), globalopts.ms.value(), globalopts.epith.value(), reference[0]);
	  }
	    acEst.pava();
	}
	else {    
	  if(globalopts.tukeysize.value() == 0)
	     globalopts.tukeysize.set_value(num2str((int)(2*sqrt(sizeTS))/2));
	   acEst.calcRaw();
	   if(globalopts.smoothACEst.value()) {
	     if ( globalopts.analysisMode.value()=="surface" ) {
	        fslSurface<float, unsigned int> topologyData;
	        read_surface(topologyData,globalopts.inputDataName2.value());
	        acEst.spatiallySmooth(topologyData,globalopts.epith.value(),globalopts.ms.value(),globalopts.tukeysize.value());
	     } else
	        acEst.spatiallySmooth(reference.matrix(mask).t(), globalopts.ms.value(), globalopts.epith.value(), reference[0], globalopts.tukeysize.value());		
	   }
	   acEst.tukey(globalopts.tukeysize.value());
	}
    }
    cout << "Completed" << endl; 

    cout << "Prewhitening and Computing PEs..." << endl;
    cout << "Percentage done:" << endl;
    int co = 1;

    Matrix mean_prewhitened_dm(paradigm.getDesignMatrix().Nrows(),paradigm.getDesignMatrix().Ncols());
    mean_prewhitened_dm=0;
    for(int i = 1; i <= numTS; i++)
    {	
      Matrix effectiveDesign(paradigm.getDesignMatrix(i,mask,labels));
      if ( (100.0*i)/numTS > co )
        cout << co++ << "," << flush;	   
      if(!globalopts.noest.value()) {
	acEst.setDesignMatrix(effectiveDesign);
	// Use autocorr estimate to prewhiten data and design:
	ColumnVector xw;
	acEst.preWhiten(datam.Column(i), xw, i, effectiveDesign);   
        datam.Column(i)=xw;
      }
      glimGls.setData(datam.Column(i), effectiveDesign, i, tContrasts, fContrasts);
      residuals.Column(i)=glimGls.getResiduals();
      if(globalopts.output_pwdata.value() || globalopts.verbose.value())
        mean_prewhitened_dm+=effectiveDesign;	
    }
     
    cout << "Completed" << endl << "Saving results... " << endl;

    if (globalopts.meanInputFile.value()=="" || globalopts.minimumTimepointFile.value()=="") {
      glimGls.saveData(logger.getDir() + "/res4d",residuals,input_data,mask,true,true,reference.tdim(),false,-1,surfaceData,globalopts.analysisMode.value());
    } else {
      int minimumTimepoint(0);
      ifstream inputTextFile(globalopts.minimumTimepointFile.value().c_str());
      if(inputTextFile.is_open()) {
	inputTextFile >> minimumTimepoint;
	inputTextFile.close();
      }
      cout << "Calculating new mean functional image using timepoint " << minimumTimepoint << endl;
      volume4D<float> residualsImage;
      volume<float> meanInput;
      read_volume4D(input_data,globalopts.inputDataName.value());
      read_volume(meanInput,globalopts.meanInputFile.value());
      residualsImage.setmatrix(residuals,mask);
      residualsImage.setDisplayMaximumMinimum(residualsImage.max(),residualsImage.min());
      save_volume4D(residualsImage,logger.getDir() + "/res4d");
      input_data-=residualsImage;
      input_data-=meanInput;
      save_volume(input_data[minimumTimepoint],"mean_func2"); 
    }

    if(globalopts.output_pwdata.value() || globalopts.verbose.value())
      {
	// Write out whitened data
        input_data.setmatrix(datam,mask);
	input_data.setDisplayMaximumMinimum(input_data.max(),input_data.min());
        save_volume4D(input_data,logger.getDir() + "/prewhitened_data");
	// Write out whitened design matrix
	write_vest(logger.appendDir("mean_prewhitened_dm.mat"), mean_prewhitened_dm/numTS);
		
      }

    // Write out threshac:
    Matrix& threshacm = acEst.getEstimates();
    int cutoff = sizeTS/2;
    if(globalopts.tukeysize.value()>0) cutoff = globalopts.tukeysize.value();
    if( globalopts.noest.value() ) cutoff=1;	
    threshacm = threshacm.Rows(1,MISCMATHS::Max(1,cutoff)); 


    glimGls.saveData(logger.getDir() + "/threshac1",threshacm,input_data,mask,true,true,reference.tdim(),true,NIFTI_INTENT_ESTIMATE,surfaceData,globalopts.analysisMode.value());
    threshacm.CleanUp();

    // save gls results:
    glimGls.Save(mask,input_data,surfaceData,globalopts.analysisMode.value(),reference.tdim());
    glimGls.CleanUp();

    cout << "Completed" << endl;
  }  
  catch(Exception p_excp) 
  {
      cerr << p_excp.what() << endl;
      return 1;
  }
  catch(...)
  {
      cerr << "Uncaught exception!" << endl;
      return 1;
  }
  return 0;
}

