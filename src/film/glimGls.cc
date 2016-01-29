/*  glimGls.cc

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

#include <sstream>
#include "glimGls.h"
#include "miscmaths/miscmaths.h"
#include "utils/log.h"
#include "fslsurface/fslsurface.h"
#include "miscmaths/t2z.h"
#include "miscmaths/f2z.h"

using namespace MISCMATHS;
using namespace Utilities;
using namespace NEWIMAGE;
using namespace fslsurface_name;


namespace FILM {

  GlimGls::GlimGls(const int pnumTS, const int psizeTS, const int pnumParams, const int pnContrasts, const int pnFContrasts) : 
    numTS(pnumTS),
    sizeTS(psizeTS),
    nContrasts(pnContrasts),
    nFContrasts(pnFContrasts),
    numParams(pnumParams),
    b(numParams, numTS),
    copes(nContrasts, numTS),
    varcopes(nContrasts, numTS),  
    fstats(nFContrasts,numTS),
    sigmaSquareds(numTS),
    dof(sizeTS - numParams)
    {
      I=IdentityMatrix(sizeTS);
      fDof.resize(nFContrasts,0);
    }
  
  void GlimGls::CleanUp()
  {
    copes.CleanUp();
    varcopes.CleanUp();
    sigmaSquareds.CleanUp();
    b.CleanUp();
    r.CleanUp();
  }

 void GlimGls::setData(const ColumnVector& y, const Matrix& x, const int ind, const Matrix& tContrasts, const Matrix& fContrasts)
 {
      // compute b
      Matrix inv_xx = pinv(x.t()*x);
      b.Column(ind) = inv_xx*x.t()*y;
      // compute r
      Matrix R = I-x*inv_xx*x.t();
      r = R*y;

      sigmaSquareds(ind) = (r.t()*r).AsScalar()/R.Trace();
      for ( int tContrast=1;tContrast <= tContrasts.Nrows();tContrast++) {
	Matrix con=tContrasts.Row(tContrast);
	copes(tContrast,ind)=(con*b.Column(ind)).AsScalar();
	double scale((con*inv_xx*con.t()).AsScalar());
	if ( scale <= 0 && con.MaximumAbsoluteValue() > 0 )
	  cerr << "Neff Error" << endl;
	varcopes(tContrast,ind)=sigmaSquareds(ind)*scale;
      }
      for ( int fContrast=1;fContrast <= fContrasts.Nrows();fContrast++) {
	Matrix fullFContrast( 0, tContrasts.Ncols() );
	for (int tContrast=1; tContrast<=fContrasts.Ncols() ; tContrast++ )
	  if (fContrasts(fContrast,tContrast)==1) fullFContrast &= tContrasts.Row(tContrast);
	fDof[fContrast-1]=fullFContrast.Nrows();
	fstats(fContrast,ind) = (b.Column(ind).t()*fullFContrast.t()*(fullFContrast*inv_xx*fullFContrast.t()*sigmaSquareds(ind)).i()*fullFContrast*b.Column(ind)).AsScalar()/(float)fullFContrast.Nrows();
      }
      // set corrections  SetCorrection(inv_xx, ind);
 }


    void GlimGls::Save(const volume<float>& mask, volume4D<float>& saveVolume, fslSurface<float, unsigned int>& saveSurface, const string& saveMode, const float reftdim)
    {
      // Need to save b, sigmaSquareds, corrections and dof 
      Log& logger = LogSingleton::getInstance();
      for(int i = 1; i <= numParams; i++) //Beta
	  saveData(logger.getDir() + "/pe" + num2str(i),b.Row(i),saveVolume,mask,true,false,-1,true,NIFTI_INTENT_ESTIMATE,saveSurface,saveMode);



      ColumnVector zstat(numTS);
      for ( int tContrast=1;tContrast <= nContrasts;tContrast++) {
	  saveData(logger.getDir() + "/cope" + num2str(tContrast),copes.Row(tContrast),saveVolume,mask,true,false,-1,true,NIFTI_INTENT_ESTIMATE,saveSurface,saveMode);
	  saveData(logger.getDir() + "/varcope" + num2str(tContrast),varcopes.Row(tContrast),saveVolume,mask,true,false,-1,true,NIFTI_INTENT_ESTIMATE,saveSurface,saveMode);
	  RowVector tstat(SD(copes.Row(tContrast),sqrt(varcopes.Row(tContrast))));
	  T2z::ComputeZStats(varcopes.Row(tContrast).AsColumn(), copes.Row(tContrast).AsColumn(), dof, zstat);
	  saveData(logger.getDir() + "/tstat" + num2str(tContrast),tstat,saveVolume,mask,true,false,-1,true,NIFTI_INTENT_TTEST,saveSurface,saveMode);
	  saveData(logger.getDir() + "/zstat" + num2str(tContrast),zstat.AsRow(),saveVolume,mask,true,false,-1,true,NIFTI_INTENT_ZSCORE,saveSurface,saveMode);
      }
      for ( int fContrast=1;fContrast <= nFContrasts;fContrast++) {
	  saveData(logger.getDir() + "/fstat" + num2str(fContrast),fstats.Row(fContrast),saveVolume,mask,true,false,-1,true,NIFTI_INTENT_FTEST,saveSurface,saveMode);
	  F2z::ComputeFStats(fstats.Row(fContrast).AsColumn(), fDof[fContrast-1], dof, zstat);
	  saveData(logger.getDir() + "/zfstat" + num2str(fContrast),zstat.AsRow(),saveVolume,mask,true,false,-1,true,NIFTI_INTENT_ZSCORE,saveSurface,saveMode);
      }
      saveData(logger.getDir() + "/sigmasquareds",sigmaSquareds,saveVolume,mask,true,false,-1,false,-1,saveSurface,saveMode);
      ColumnVector dofVec(1);
      dofVec = dof;
      write_ascii_matrix(logger.appendDir("dof"), dofVec);    
    }

  void GlimGls::saveData(const string& outputName, const Matrix& data, volume4D<float>& saveVolume, const volume<float>& volumeMask, const  bool setVolumeRange, const bool setVolumeTdim, const int outputTdim, const bool setIntent, const int intentCode, fslSurface<float, unsigned int>& saveSurface, const string& saveMode)
  {
    if ( saveMode=="surface" ) {
      saveSurface.reinitialiseScalars(data.Nrows());
      for(unsigned int vertex=0; vertex < saveSurface.getNumberOfVertices(); vertex++)
	for(unsigned int timepoint=0; timepoint < saveSurface.getNumberOfScalarData(); timepoint++)
           saveSurface.setScalar(timepoint,vertex,data(timepoint+1,vertex+1));
      writeGIFTI(saveSurface,outputName+".func.gii",GIFTI_ENCODING_B64GZ);
    } else { 
      saveVolume.setmatrix(data,volumeMask);
      if ( setVolumeTdim ) 
	saveVolume.settdim(outputTdim);
      if ( setIntent )
	saveVolume.set_intent(intentCode,0,0,0);
      if ( setVolumeRange )
	saveVolume.setDisplayMaximumMinimum(saveVolume.max(),saveVolume.min());
      save_volume4D(saveVolume,outputName);
    }
  }


}







