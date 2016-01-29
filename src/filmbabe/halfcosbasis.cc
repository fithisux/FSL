/*  halfcosbasis.cc

    Mark Woolrich - FMRIB Image Analysis Group

    Copyright (C) 2004 University of Oxford  */

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
#include <fstream>
#include <iomanip>
#include <sstream>
#define WANT_STREAM
#define WANT_MATH
#include "newmatap.h"
#include "newmatio.h"
#include <string>
#include <math.h>
#include "utils/log.h"
#include "halfcosbasisoptions.h"
#include "utils/tracer_plus.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include "libvis/miscplot.h"

using namespace Utilities;
using namespace NEWMAT;
using namespace MISCMATHS;
using namespace MISCPLOT;
using namespace Halfcosbasis;

inline float halfcos(float t,float ybot,float ytop,float xleft,float xright,int flipud)
{
  float y;

  if(xright>xleft && t>xleft && t<=xright)
    y=((ytop-ybot)/2.0*flipud*std::cos(2*M_PI*(t-xleft)/((xright-xleft)*2))+(ytop-ybot)/2+ybot);
  else
    y=0;

  return y;
}

float halfcos_hrf(float t,  const ColumnVector& halfcosparams)
  //float m1,float m2,float m3,float m4,float m5,float c1,float c2,float c3)
{
  float m1 = halfcosparams(1);
  float m2 = halfcosparams(2);
  float m3 = halfcosparams(3);
  float m4 = halfcosparams(4);
  float m5 = halfcosparams(5);
  float c1 = halfcosparams(6);
  float c2 = halfcosparams(7);
  float c3 = halfcosparams(8);

  float y=halfcos(t,-1*c1,0,0,m1,1);
  y+=halfcos(t,-1*c1,1,m1,m1+m2,-1);
  y+=halfcos(t,-c2,1,m1+m2,m1+m2+m3,1);
  y+=halfcos(t,-c2,c3,m1+m2+m3,m1+m2+m3+m4,-1);
  y+=halfcos(t,0,c3,m1+m2+m3+m4,m1+m2+m3+m4+m5,1);

  return y;
}


ReturnMatrix halfcos_hrf(const ColumnVector& t, const ColumnVector& halfcosparams)
{
  Tracer_Plus trace("halfcos_hrf");

  int T = t.Nrows();
  ColumnVector y(T);
  y = 0;
  
  for(int i = 1; i <= T; i++)
    {
      //	OUT(i);
      y(i) = halfcos_hrf(t(i),halfcosparams);
    }
  
  y.Release();
  return y;
}

ReturnMatrix samplehalfcosparams(const Matrix& halfcosparamsranges)
{
  Tracer_Plus trace("samplehalfcosparams");

  int nparams = halfcosparamsranges.Nrows();
  //OUT(nparams);

  ColumnVector halfcosparams(nparams);
  halfcosparams = 0;
  
  for(int i = 1; i <= nparams; i++)
    {
      //	OUT(i);
      halfcosparams(i) = unifrnd(1,1,halfcosparamsranges(i,1),halfcosparamsranges(i,2)).AsScalar();
    }
  
  halfcosparams.Release();
  return halfcosparams;
}

int main(int argc, char *argv[])
{
  try{
  
    // Setup logging:
    Log& logger = LogSingleton::getInstance();
    
    // parse command line - will output arguments to logfile
    HalfcosbasisOptions& opts = HalfcosbasisOptions::getInstance();
    opts.parse_command_line(argc, argv, logger);

    //Tracer_Plus::setinstantstackon();

    Matrix halfcosparamranges = read_ascii_matrix(opts.halfcosparamrangesfile.value());
    
    if(opts.verbose.value())
      OUT(halfcosparamranges);

    if(halfcosparamranges.Nrows() != 8 || halfcosparamranges.Ncols() != 2)
      {
	throw Exception("Invalid halfcosparamranges file");
      }

    int T = int(opts.nsecs.value()/opts.res.value());

    if(opts.verbose.value()) 
      {
	OUT(T);
	OUT(opts.nsecs.value());
	OUT(opts.res.value());
	OUT(opts.nhrfsamps.value());
      }

    // setup time vector
    ColumnVector t(T);
    t = 0;
    for(int i = 1; i <= T; i++)
      {
	//	OUT(i);
	t(i) = i*opts.res.value();
      }

    if(opts.verbose.value())
      OUT("Generating HRF samples");
    Matrix hrfsamps(T,opts.nhrfsamps.value());
    hrfsamps = 0;

    for(int i = 1; i <= opts.nhrfsamps.value(); i++)
      {
	//OUT(i);
	ColumnVector halfcosparams = samplehalfcosparams(halfcosparamranges);
	hrfsamps.Column(i) = halfcos_hrf(t,halfcosparams);
      }

    if(opts.verbose.value())
	OUT("Finished generating HRF samples");

    if(opts.verbose.value())
	OUT("Performing SVD");

    DiagonalMatrix eigenvals;
    Matrix eigenvecs;    

    try
      {
        //SymmetricMatrix corr;
        //corr << hrfsamps*hrfsamps.t()/opts.nhrfsamps.value();
        //EigenValues(corr, eigenvals, eigenvecs);

        // use svd instead of eigenvalues - uses much less memory
        if (hrfsamps.Ncols()>hrfsamps.Nrows())  // Ncols has to be <= Nrows for SVD()
          {
            Matrix grot;
            SVD(hrfsamps.t(), eigenvals, grot, eigenvecs);
          }
        else
          SVD(hrfsamps, eigenvals, eigenvecs);
        eigenvals = (eigenvals * eigenvals).Reverse(); // square and reverse order
        for(int i = 1; i <= eigenvecs.Nrows(); i++)
          eigenvecs.Row(i)=eigenvecs.Row(i).Reverse(); // reverse order of columns
      }
    catch(Exception& e) 
      {
	if(opts.verbose.value())
	  cerr << endl << e.what() << endl;
	
	throw e;
      }
	
    if(opts.verbose.value())
      {
	OUT("Finished SVD");
	write_ascii_matrix(eigenvecs, LogSingleton::getInstance().appendDir("hrfeigenvecs.txt"));
	write_ascii_matrix(diag(eigenvals), LogSingleton::getInstance().appendDir("hrfeigenvals.txt"));
      }

    // need to flip cos of the way EigenValues does its output
    eigenvecs = eigenvecs.Reverse();
    eigenvals = eigenvals.Reverse();

    Matrix basisfns = eigenvecs.Columns(1,opts.nbfs.value());

    // need to flip (in time) basis function time series cos of the way EigenValues does its output
    for(int i=1; i<=opts.nbfs.value(); i++)
      {
	basisfns.Column(i) = basisfns.Column(i).Reverse();
      }

//     OUT(basisfns.Column(1).Minimum());
//     OUT(basisfns.Column(1).Maximum());

    // find max and min of biggest component and flip all (in amp) if abs(min)>abs(max)
    if(abs(basisfns.Column(1).Minimum())>abs(basisfns.Column(1).Maximum()))
      {
	basisfns = -basisfns;
      }

    // Output hrf samples
    write_ascii_matrix(hrfsamps, LogSingleton::getInstance().appendDir("hrfsamps.txt"));
    miscplot newplot;
    newplot.add_xlabel("time (secs)");    
    newplot.set_xysize(610,300);
    newplot.timeseries(hrfsamps.t(), LogSingleton::getInstance().appendDir("hrfsamps"), "HRF Samples", opts.res.value(),400,3,0,false);

    // output basis functions
    write_ascii_matrix(basisfns, LogSingleton::getInstance().appendDir("hrfbasisfns.txt"));
    miscplot newplot2;
    for(int i = 1; i <= opts.nbfs.value(); i++)
      {
	newplot2.add_label(string("Basis fn ")+num2str(i));
      }
    newplot2.add_xlabel("time (secs)");    
    newplot2.set_xysize(600,300);
    newplot2.timeseries(basisfns.t(), LogSingleton::getInstance().appendDir("hrfbasisfns"), "HRF Basis Functions", opts.res.value(),400,3,0,false);    

    // find number of eigenvalues to display:    
    float sumeigs = eigenvals.Sum();    
    float runsum = 0.0;
    int numeigs = 1;    
    
    if(opts.verbose.value())
      OUT(sumeigs);

    for(; numeigs < 100; numeigs++)
      {
	runsum += eigenvals(numeigs,numeigs);
	
	if(runsum>0.995*sumeigs) break;
      }
    
    if(opts.verbose.value())
      {
	OUT(runsum);
	OUT(numeigs);
      }

    miscplot newplot3;
    newplot3.add_xlabel("basis function number");    
    newplot3.set_xysize(300,300);
    newplot3.timeseries((diag(eigenvals).Rows(1,numeigs)/eigenvals(1,1)).t(), LogSingleton::getInstance().appendDir("eigenvalues"), "Normalised Eigenvalues (99.5% of variance)", 0,400,3,0,false);    

    if(opts.verbose.value())
      OUT("Normalising HRF samples and basis set");

//     // remove mean and variance from basis functions and hrf samples
//     for(int i = 1; i <= opts.nbfs.value(); i++)
//       {
// 	ColumnVector tmpcol = basisfns.Column(i);
// 	basisfns.Column(i) = (tmpcol - mean(tmpcol).AsScalar())/stdev(tmpcol).AsScalar();
//       }

//     for(int i = 1; i <= opts.nhrfsamps.value(); i++)
//       {
// 	ColumnVector tmpcol = hrfsamps.Column(i);
// 	hrfsamps.Column(i) = (tmpcol - mean(tmpcol).AsScalar())/stdev(tmpcol).AsScalar();
//       }

    // remove mean and variance from basis functions and hrf samples
    for(int i = 1; i <= opts.nbfs.value(); i++)
      {
	ColumnVector tmpcol = basisfns.Column(i);
	basisfns.Column(i) = (tmpcol - mean(tmpcol).AsScalar());
      }

    for(int i = 1; i <= opts.nhrfsamps.value(); i++)
      {
	ColumnVector tmpcol = hrfsamps.Column(i);
	hrfsamps.Column(i) = (tmpcol - mean(tmpcol).AsScalar());
      }

    if(opts.verbose.value())
      OUT("Regressing HRF samples onto basis set");

    Matrix paramsamps = pinv(basisfns)*hrfsamps;
    //write_ascii_matrix(paramsamps, LogSingleton::getInstance().appendDir("paramsamps"));
    
    if(opts.verbose.value())
      OUT("Fitting MVN parameter constraints");

    ColumnVector priormeans = mean(paramsamps.t()).t();
    Matrix priorcovars = cov(paramsamps.t());    

    // output constraints
    write_vest(priormeans, LogSingleton::getInstance().appendDir("priormeans.mat"));
    write_vest(priorcovars, LogSingleton::getInstance().appendDir("priorcovars.mat"));

    if(opts.verbose.value())
      {
	OUT(priormeans);
	OUT(priorcovars);
      }

    if(opts.debuglevel.value()==1)
      Tracer_Plus::setrunningstackon();

    if(opts.timingon.value())
      Tracer_Plus::settimingon();
    
    if(opts.timingon.value())
      Tracer_Plus::dump_times(logger.getDir());
    
    cout << "Log directory was: " << logger.getDir() << endl;

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

  return 0;
}












