/*  ContrastMgr.cc

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
#include <fstream>
#include "ContrastMgr.h"
#include "ContrastMgrOptions.h"
#include "newimage/newimageall.h"
#include "utils/log.h"
#include "gaussComparer.h"
#include "miscmaths/t2z.h"
#include "miscmaths/f2z.h"
#include "paradigm.h"
#include "utils/tracer_plus.h"

using namespace NEWIMAGE;
using namespace Utilities;
using namespace MISCMATHS;
using namespace std;

namespace FILM {

  ContrastMgr::ContrastMgr() :    
    tc(),
    fc(),
    c_counter(0),
    numParams(0),
    num_Ccontrasts_in_Fcontrast(0),
    contrast_valid(false),
    contrast_num(0),
    parad(),
    corrections(),
    b(),
    dof(),
    sigmaSquareds(),
    varcb(),
    cb(),
    neff(),
    numTS(0)
  {}

  void ContrastMgr::SetFContrast(const int p_num, const int p_c_counter)
  {
    fc.ReSize(parad.getTContrasts().Nrows(),numParams);
    int count = 0;

    for(int c = 1; c <= parad.getTContrasts().Nrows(); c++)
      {	
	if(parad.getFContrasts()(p_num,c) == 1)
	  {	    
	    count++;
	    fc.Row(count) = parad.getTContrasts().Row(c);	
	  }
      }

    fc = fc.Rows(1,count);
    contrast_num = p_num;
    contrast_valid = true;

    num_Ccontrasts_in_Fcontrast = fc.Nrows();

    c_counter = p_c_counter;
  }

  void ContrastMgr::run()
  {
    Tracer_Plus ts("ContrastMgr::run");
    Load();

    // Loop through tcontrasts:
    for(int c = 1; c <= parad.getTContrasts().Nrows(); c++)
      {
	if(ContrastMgrOptions::getInstance().verbose)
	  {
	    cerr << "T contrast no. " << c << endl;
	    cerr << parad.getTContrasts().Row(c) << endl;
	  }

	SetTContrast(c, c+ContrastMgrOptions::getInstance().copenumber-1);
	ComputeNeff();    
	ComputeCope();
	ComputeVarCope();
	ComputeZStat();
	
	SaveTContrast(ContrastMgrOptions::getInstance().suffix);
      }

    // Loop through fcontrasts:
    for(int c = 1; c <= parad.getFContrasts().Nrows(); c++)
      {
	SetFContrast(c, c+ContrastMgrOptions::getInstance().copenumber-1);

	if(ContrastMgrOptions::getInstance().verbose)
	  {
	    cerr << "F contrast no." << c << endl;
	    cerr << parad.getFContrasts().Row(c) << endl;
	    cerr << fc << endl;
	  }

	ComputeFStat();
	
	if(contrast_valid)
	  SaveFContrast(ContrastMgrOptions::getInstance().suffix);

      }
  }

  void ContrastMgr::Load()
  {
    Tracer_Plus ts("ContrastMgr::Load");
    // Need to read in b, sigmaSquareds, corrections and dof 
    Log& logger = LogSingleton::getInstance();
    // Load contrasts:     
    parad.load("", ContrastMgrOptions::getInstance().contrastfname, ContrastMgrOptions::getInstance().fcontrastfname, false, 0);
       
    numParams = parad.getTContrasts().Ncols(); 
    if(ContrastMgrOptions::getInstance().verbose)
      {
	cerr << "T Contrasts:" << endl << parad.getTContrasts();
	cerr << "F Contrasts:" << endl << parad.getFContrasts();
      }
    volume4D<float> input_data;
    read_volume4D(input_data,logger.getDir() + "/sigmasquareds");

    mask=input_data[0];
    mask.binarise(0.0,mask.max()+1,exclusive);
    sigmaSquareds=input_data.matrix(mask).AsColumn();
    numTS = sigmaSquareds.Nrows();    
    // b:
    ColumnVector peVol;
    b.ReSize(numTS, numParams);
    for(int i = 1; i <= numParams; i++)
      { 
	// Add param number to "pe" to create filename:
	ostringstream osc;
	osc << i;
	read_volume4D(input_data,logger.getDir() + "/pe" + osc.str().c_str());
	peVol=input_data.matrix(mask).AsColumn();
	b.Column(i) = peVol; 
      }

    // dof: - maybe single value ASCII or avw file:
    ifstream in;
    in.open(string(logger.getDir() + "/dof").c_str(), ios::in);    
    if(!in) //avw format
      {
	read_volume4D(input_data,logger.getDir() + "/dof");
	dof=input_data.matrix(mask).AsColumn();
      }
    else
      {
	// single value ascii format
	in.close();
	ColumnVector dofVec = MISCMATHS::read_ascii_matrix(logger.getDir() + "/dof");
	dof = sigmaSquareds;
	dof = dofVec(1);
      }      

    // corrections - maybe ASCII (old version) or avw file:    
    // corrections are the correlation matrix of the pes.
    in.open(string(logger.getDir() + "/corrections").c_str(), ios::in);    
    if(!in)
      {
	// avw format
	is_avw_corrections = true;
	read_volume4D(input_data,logger.getDir() + "/corrections");
	corrections=input_data.matrix(mask);
      }
    else
      {
	// old ascii format
	in.close();      
	is_avw_corrections = false;
	corrections.ReSize(numTS,numParams*numParams);
	corrections=read_vest(logger.getDir() + "/corrections");
	corrections = corrections.t();
      }
  }

  void ContrastMgr::SaveFContrast(const string& suffix)
  {
    Tracer_Plus ts("ContrastMgr::SaveFContrast");
    Log& logger = LogSingleton::getInstance();

    // prepare contrast number:
    ostringstream osc;
    osc << suffix << c_counter;
    volume4D<float> output;
    output.setmatrix(fstat.AsRow(),mask);
    output.setDisplayMaximumMinimum(output.max(),output.min());
    output.set_intent(NIFTI_INTENT_FTEST,0,0,0);
    save_volume4D(output,logger.getDir() + "/fstat" + osc.str().c_str());
    output.setmatrix(zstat.AsRow(),mask);
    output.setDisplayMaximumMinimum(output.max(),output.min());
    output.set_intent(NIFTI_INTENT_ZSCORE,0,0,0);
    save_volume4D(output,logger.getDir() + "/zfstat" + osc.str().c_str());
  }

  void ContrastMgr::SaveTContrast(const string& suffix)
  {
    Tracer_Plus ts("ContrastMgr::SaveTContrast");
    Log& logger = LogSingleton::getInstance();
    // prepare contrast number:
    ostringstream osc;
    osc << suffix << c_counter;
    volume4D<float> output;
    output.setmatrix(neff.AsRow(),mask);
    output.setDisplayMaximumMinimum(output.max(),output.min());
    output.set_intent(NIFTI_INTENT_NONE,0,0,0);
    save_volume4D(output,logger.getDir() + "/neff" + osc.str().c_str());
    output.setmatrix(cb.AsRow(),mask);
    output.setDisplayMaximumMinimum(output.max(),output.min());
    output.set_intent(NIFTI_INTENT_ESTIMATE,0,0,0);
    save_volume4D(output,logger.getDir() + "/cope" + osc.str().c_str());
    output.setmatrix(varcb.AsRow(),mask);
    output.setDisplayMaximumMinimum(output.max(),output.min());
    output.set_intent(NIFTI_INTENT_ESTIMATE,0,0,0);
    save_volume4D(output,logger.getDir() + "/varcope" + osc.str().c_str());
    output.setmatrix(tstat.AsRow(),mask);
    output.setDisplayMaximumMinimum(output.max(),output.min());
    output.set_intent(NIFTI_INTENT_TTEST,0,0,0);
    save_volume4D(output,logger.getDir() + "/tstat" + osc.str().c_str());
    output.setmatrix(zstat.AsRow(),mask);
    output.setDisplayMaximumMinimum(output.max(),output.min());
    output.set_intent(NIFTI_INTENT_ZSCORE,0,0,0);
    save_volume4D(output,logger.getDir() + "/zstat" + osc.str().c_str());
  }


  void ContrastMgr::GetCorrection(Matrix& corr, const int ind)
  {
    Tracer_Plus ts("ContrastMgr::GetCorrection");

    // puts ColumnVector of length p*p from correction
    // into Matrix corr which is p*p:
    corr.ReSize(numParams, numParams);
    for (int i = 1; i <= numParams; i++)
      {
	for (int j = 1; j <= numParams; j++)
	  {
	    corr(i,j) = corrections((i-1)*numParams + j, ind); 
	  }
      }
  }

  void ContrastMgr::ComputeZStat()
  {
    Tracer_Plus ts("ContrastMgr::ComputeZStat");

    Log& logger = LogSingleton::getInstance();

    // calulate Zstat:
    tstat.ReSize(numTS);
    for(int i = 1; i <= numTS; i++)
      {
	if(varcb(i) > 0 && neff(i) > 0)
	  {
	    tstat(i) = cb(i)/sqrt(varcb(i));
	  }
	else
	  tstat(i) = 0.0;
      }
      
    // Calculate tstat:
    zstat.ReSize(numTS);
    T2z::ComputeZStats(varcb, cb, dof, zstat);

    // Compare with theory:
    GaussComparer gaussComp(zstat);
    gaussComp.setup();

    ColumnVector ratios(5);
    ColumnVector probs(5);
    int co = 1;
    for(float p = 0.05; p >= 0.0005; p=p/sqrt(10))
      {
	float temp = gaussComp.computeRatio(p, logger.str());
	logger.str() << "p = " << p << ": ratio = " << temp << endl;
	ratios(co) = temp;
	probs(co) = p;
	co++;
      }
      
    write_ascii_matrix(logger.appendDir("ratios"), ratios);
    write_ascii_matrix(logger.appendDir("probs"), probs);
      
  }

  void ContrastMgr::ComputeCope()
  { 
    Tracer_Plus ts("ContrastMgr::ComputeCope");
    cb.ReSize(numTS);

    for(int i = 1; i <= numTS; i++)
      {	
	cb(i) = (tc.t()*b.Row(i).t()).AsScalar();
      }
  }

  void ContrastMgr::ComputeNeff()
  {
    Tracer_Plus ts("ContrastMgr::ComputeNeff");
   
    //Log& logger = LogSingleton::getInstance();
    Matrix corr;
      
    neff.ReSize(numTS);
    neff = 0;

    int numNegs = 0;
    int maxNeff = 0;

    for(int i = 1; i <= numTS; i++)
      {	
	GetCorrection(corr, i);

	neff(i) = 1/(tc.t()*corr*tc).AsScalar();

	if(maxNeff < neff(i))
	  maxNeff = (int)neff(i);
	  
	if(neff(i) < 0.0 && i > 1)
	  {
	    neff(i) = neff(i-1);
	    numNegs++;
	  }  
      }
  }

  void ContrastMgr::ComputeFStat()
  {
    Tracer_Plus ts("ContrastMgr::ComputeFStat");
 
    //Log& logger = LogSingleton::getInstance();
    Matrix corr;

    fstat.ReSize(numTS);
    fstat = 1;

    for(int i = 1; i <= numTS; i++)
      {
	GetCorrection(corr, i);

	try
	  {	    	    
	    fstat(i) = (b.Row(i)*fc.t()*(fc*corr*fc.t()*sigmaSquareds(i)).i()*fc*b.Row(i).t()).AsScalar()/num_Ccontrasts_in_Fcontrast;
	  }
	catch(SingularException& ex)
	  {	    
	    cerr << ex.what() << endl;
	    cerr << "F contrast no. " << contrast_num << " produces singular variance matrix." << endl; 
	    cerr << "No results will be produced for this contrast"  << endl;
	    contrast_valid = false;
	    break;	    
	  }
      }
   
    // Calculate zstat:
    zstat.ReSize(numTS);
    F2z::ComputeFStats(fstat, num_Ccontrasts_in_Fcontrast, dof, zstat);
  }

  void ContrastMgr::ComputeVarCope()
  { 
    Tracer_Plus ts("ContrastMgr::ComputeVarCope");
    varcb.ReSize(numTS);

    for(int i = 1; i <= numTS; i++)
      {
	if(neff(i) > 0)
	  varcb(i) = sigmaSquareds(i)/neff(i);
	else
	  varcb(i) = 0;
      }
  }

}










