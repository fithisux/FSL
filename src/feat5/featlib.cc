/* {{{ Copyright etc. */

/*  featlib - FMRI time series and model plotting

    Stephen Smith and Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2002 University of Oxford  */

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

/* }}} */
/* {{{ defines, includes and typedefs */
#include "featlib.h"

/* }}} */
/* {{{ find_key */

char *find_key(FILE *fd, char *control_line, const char *key)
{
  while (fgets(control_line, 1000, fd)!=NULL)
    if ( (strncmp(control_line,key,strlen(key))==0) )
      return control_line+strlen(key);
	  
  printf("Error: key \"%s\" not found.\n",key);
  exit(1);
}

/* }}} */
/* {{{ establish_pwfilter */

void establish_pwfilter(const ColumnVector& ac, ColumnVector& pwfilter, int zeropad, int npts)
{
  // FFT auto corr estimate
  ColumnVector dummy(zeropad);
  dummy = 0;
  ColumnVector vrow(zeropad);
  vrow = 0;

  // ac maybe cutoff to be smaller than npts
  int sizeTS = ac.Nrows();
  if (sizeTS > npts/2)
    sizeTS = npts/2;
  
  vrow.Rows(1,sizeTS) = ac.Rows(1,sizeTS);
  vrow.Rows(zeropad - sizeTS + 2, zeropad) = ac.Rows(2, sizeTS).Reverse();

  ColumnVector ac_fft_imag;
    
  FFT(vrow, dummy, pwfilter, ac_fft_imag);      

  // inverse auto corr to give prewhitening filter
  // no DC component so set first value to 0
  pwfilter(1) = 0.0;
  
  for(int j = 2; j <= zeropad; j++)
    {
      if (pwfilter(j)<0)
	{
	  cout << "Warning: possible high autocorrelation in time series" << endl;
	  pwfilter(j)=0;
	}
      else
	pwfilter(j) = 1.0/sqrt(pwfilter(j));	      
    }  

  // normalise pwfilter such that sum(j)((pwfilter)^2/zeropad)) = 1
  pwfilter /= sqrt(pwfilter.SumSquare()/zeropad);
}

/* }}} */
/* {{{ prewhiten_data */

void prewhiten_data(const ColumnVector& data, ColumnVector& pwdata, ColumnVector& pwfilter, int zeropad, int npts)
{ 
  ColumnVector data_fft_real, data_fft_imag, realifft, dummy;
  dummy.ReSize(zeropad);
  dummy = 0;

  // Remove and store mean
  float mn = MISCMATHS::mean(data).AsScalar();
  pwdata.ReSize(zeropad);
  pwdata = 0;
  pwdata.Rows(1,npts) = data - mn;

  // FFT data
  FFT(pwdata, dummy, data_fft_real, data_fft_imag);
  FFTI(SP(pwfilter, data_fft_real), SP(pwfilter, data_fft_imag), realifft, dummy); 
  // take first npts and restore mean
  pwdata = realifft.Rows(1,npts) + mn;  
}

/* }}} */
/* {{{ prewhiten_model */

void prewhiten_model(const ColumnVector& ac,vector<double>& model,vector<double>& pwmodel, int nevs, int npts)
{
  Matrix dm(npts, nevs), pwdm(npts, nevs);

  // Get values out of model
  for(int ev=0; ev<nevs; ev++)
    for(int t=0; t<npts; t++)
      dm(t+1,ev+1) = model[t*nevs+ev];

  pwdm = dm;

  int zeropad = (int)pow(2,ceil(log(npts)/log(2)));
  
  // get prewhitening filter from ac
  ColumnVector pwfilter;
  establish_pwfilter(ac, pwfilter, zeropad, npts);

  // prewhiten each of the evs
  for(int ev = 1; ev <= nevs; ev++)
    {
      ColumnVector pwdata;
      prewhiten_data(dm.Column(ev), pwdata, pwfilter, zeropad, npts);
      pwdm.Column(ev) = pwdata;
    }

  // setup and return result
  for(int ev=0; ev<nevs; ev++)
    for(int t=0; t<npts; t++)
      pwmodel[t*nevs+ev] = pwdm(t+1,ev+1);
}

/* }}} */
/* {{{ prewhiten_timeseries */

void prewhiten_timeseries(const ColumnVector& ac, const ColumnVector& ts, ColumnVector& pwts, int npts)
{
  int zeropad = (int)pow(2,ceil(log(npts)/log(2)));
  
  // get prewhitening filter from ac
  ColumnVector pwfilter;
  establish_pwfilter(ac, pwfilter, zeropad, npts);

  // prewhiten
  prewhiten_data(ts, pwts, pwfilter, zeropad, npts);
}

/* }}} */
/* {{{ read_model */

vector<double> read_model(const string& filename, int *nevs, int *npts)
{
  FILE *designFile;
  char control_line[1010];
  vector<double> model;
  int i;

  if((designFile=fopen(filename.c_str(),"r"))==NULL)
    {
      *npts=0;
      return model;
    }

  *nevs=(int)atof(find_key(designFile,control_line,"/NumWaves"));
  *npts=(int)atof(find_key(designFile,control_line,"/NumPoints"));

  atof(find_key(designFile,control_line,"/Matrix"));

  model.resize( *nevs * *npts );

  for(i=0; i<*nevs * *npts; i++)
    fscanf(designFile,"%lf",&model[i]);

  return model;
}

/* }}} */
/* {{{ read_ftests */

void read_ftests(const string& filename, int *nftests)
{
  FILE *designfp;
  char control_line[1010];

  *nftests=0;

  if((designfp=fopen(filename.c_str(),"r"))!=NULL)
    *nftests=(int)atof(find_key(designfp,control_line,"/NumContrasts"));
}

/* }}} */
/* {{{ read_triggers */

bool read_triggers(const string& filename,vector<double>& triggers, int nevs, int npts)
{
  ifstream triggerFile(filename.c_str());
  if ( triggerFile.is_open() )
  {
    triggers.resize(nevs*2*npts);
    for(int ev=0;ev<nevs;ev++)
    {
      int i=1;
      string input;
      getline(triggerFile,input);
      stringstream inputStream(input);
      while (inputStream >> input)
	triggers[(i++)*nevs+ev]=atof(input.c_str());
      triggers[ev]=i-2;
    }
    triggerFile.close();
    return true;
  }
  return false;
}

/* }}} */
