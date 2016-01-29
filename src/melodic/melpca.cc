/*  MELODIC - Multivariate exploratory linear optimized decomposition into 
              independent components
    
    melpca.cc - PCA and whitening 

    Christian F. Beckmann, FMRIB Analysis Group
    
    Copyright (C) 1999-2013 University of Oxford */

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

#include "utils/log.h"
#include "meloptions.h"
#include "meldata.h"
#include "melodic.h"
#include "newmatap.h"
#include "newmatio.h"
#include "melpca.h"
#include "melhlprfns.h"
#include "libprob.h"

using namespace Utilities;

namespace Melodic{

  void MelodicPCA::perf_pca(Matrix& in, Matrix& weights){    
  	message("Starting PCA  ... ");

    SymmetricMatrix Corr;
    Matrix tmpE;
    RowVector tmpD, AdjEV, PercEV;
   
	if(opts.paradigmfname.value().length()>0)
	{
		basicGLM tmpglm;
		tmpglm.olsfit(in,melodat.get_param(),IdentityMatrix(melodat.get_param().Ncols()));
		std_pca(tmpglm.get_residu(),weights,Corr,tmpE,tmpD);
	}
	else{
 		std_pca(in,weights,Corr,tmpE,tmpD); 
	}	
	
    if(opts.tsmooth.value()){
      message(endl << "  temporal smoothing of Eigenvectors " << endl);
      tmpE=smoothColumns(tmpE);
    }
   
    adj_eigspec(tmpD,AdjEV,PercEV);
    melodat.set_pcaE(tmpE);
    melodat.set_pcaD(tmpD);
    melodat.set_EVP(PercEV); 
    melodat.set_EV((AdjEV));
    write_ascii_matrix(logger.appendDir("eigenvalues_percent"),PercEV);
   
    message("done" << endl);
  }
  
  void MelodicPCA::perf_white(){
    int N = melodat.get_pcaE().Ncols();    
    if(opts.pca_dim.value() > N){
      message("dimensionality too large - using -dim " << N << 
              " instead " << endl);
      opts.pca_dim.set_T(N);
    }
    if(opts.pca_dim.value() < 0){
      if(opts.remove_meanvol.value()){
				opts.pca_dim.set_T(N-2);
      }else{
   			opts.pca_dim.set_T(N-1);
      }
    }
    if(opts.pca_dim.value() ==0){
      opts.pca_dim.set_T(pcadim());
      if(melodat.get_Data().Nrows()<20){
				opts.pca_dim.set_T(N-2);
				message("too few data points for full estimation, using "
					<< opts.pca_dim.value() << " instead"<< endl);
      }
    }
    if(opts.approach.value()==string("jade") && opts.pca_dim.value() > 30){
      message("dimensionality too large for jade estimation - using --dim 30 instead" << endl);
      opts.pca_dim.set_T(30);
    }

    message("Start whitening using  "<< opts.pca_dim.value()<<" dimensions ... " << endl);
    Matrix tmpWhite;
    Matrix tmpDeWhite;

    float varp = 1.0;

	if (opts.paradigmfname.value().length()>0)
	    varp = calc_white(melodat.get_pcaE(),melodat.get_pcaD(), 
    		melodat.get_EVP(),opts.pca_dim.value(),melodat.get_param(),melodat.get_paramS(),tmpWhite,tmpDeWhite);
	else
		varp = calc_white(melodat.get_pcaE(),melodat.get_pcaD(), 
    		melodat.get_EVP(),opts.pca_dim.value(),tmpWhite,tmpDeWhite);

    melodat.set_white(tmpWhite);
    melodat.set_dewhite(tmpDeWhite);
    message("  retaining "<< varp*100 <<" percent of the variability " << endl);
    message(" ... done"<< endl << endl);
  }

  int MelodicPCA::pcadim()
  { 
    message("Estimating the number of sources from the data (PPCA) ..." << endl);

    ColumnVector PPCA; RowVector AdjEV, PercEV;   
    int res = ppca_dim(melodat.get_Data(),melodat.get_RXweight(), PPCA, 
			AdjEV, PercEV, melodat.get_resels(), opts.pca_est.value());
     
    write_ascii_matrix(logger.appendDir("PPCA"),PPCA);			      
    write_ascii_matrix(logger.appendDir("eigenvalues_adjusted"),AdjEV.t());
    write_ascii_matrix(logger.appendDir("eigenvalues_percent"),PercEV.t());

    melodat.set_EVP(PercEV);
    melodat.set_EV(AdjEV);
    melodat.set_PPCA(PPCA);
    
    return res;
  }

}


