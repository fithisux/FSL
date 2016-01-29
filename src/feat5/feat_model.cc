// {{{ Copyright etc. 

/*  feat_model - create FEAT design matrix, contrasts etc.

    Stephen Smith and Matthew Webster, FMRIB Analysis Group

    Copyright (C) 1999-2010 University of Oxford  */
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

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

#include "featlib.h"

#include "miscmaths/miscprob.h"
#include "miscmaths/t2z.h"
#include "libprob.h"

#include <cstring>
#include <sstream>

using namespace MISCMATHS;
using namespace NEWIMAGE;

#define DT                0.05 /* temporal sampling for initial model, in seconds */
#define NEGSECS           30.0 /* amount of negative modelling allowed, for custom 3 etc, in seconds */

class Contrasts {
public:
  Matrix C;
  vector<string> name;
  ColumnVector realHeights;
  Contrasts( const int nRealEVs, const int nContrasts );
  int nC() const { return C.Ncols(); };
  int nEV() const { return C.Nrows(); };
  ColumnVector  RE;
};

Contrasts::Contrasts( const int nRealEVs, const int nContrasts )
{
  C.ReSize(nRealEVs,nContrasts); //this is the transpose of FEAT's design.con format
  realHeights.ReSize(nContrasts);
  C=0;
  name.resize(nContrasts);
  RE.ReSize(nContrasts);
}

class Regressors {
public:
  int nReal;
  int nRealNoMotion;
  int nOriginal;
  int nOriginalNoMotion;
  int nMotion;
  int nTimepoints;
  vector<int> nRealPerOrig;
  vector<bool> usesDeriv;
  vector<bool> usesFilter;
  Matrix orthogonals;
  Matrix copyOfRealDesign;
  ColumnVector realHeights;
  DiagonalMatrix eigenvals;
  vector<int> basisorth; // whether or not to orth basis functions wrt each other
  int equivalentRealLocation(const int orginalLocation);
};

int Regressors::equivalentRealLocation(const int originalLocation)
{
  int currentRealEV(0);
  for (int ev=1; ev < originalLocation; ev++)
    currentRealEV+=nRealPerOrig[ev-1];
  return(++currentRealEV);
}

Regressors peakAndFilter(const bool estimateRealHeights, Regressors& evs, Contrasts& contrasts, const double hp_sigma, const double lp_sigma, const int nTimepoints, const Matrix& preWhitening, const float critical_t, const float noise, const string filename);
void writeSubmodel( const string filename, const Regressors& evs);
void writeDesignMatrix(const string filenameRoot, const ColumnVector& realDesignHeights, const Matrix& realDesign);
void writeTriggers(const string filename, const vector<float>& triggers, const int* shape, const Regressors& evs, const double factor);
void writeContrastMatrix(const string filename, const Contrasts& contrasts, const ColumnVector& RE );
void writeFMatrix(const string filename, const Contrasts& contrasts, const Matrix& F, const int nFTests, const Matrix& realDesign );
void writeGroupFile(const string filename, const vector<int> G, const int nTimepoints, const Matrix& realDesign);

ColumnVector filterRegressor(const ColumnVector& inputRegressor, const double hp_sigma, const double lp_sigma)
{
  ColumnVector outputRegressor(inputRegressor);
  volume4D <float> im(1,1,1,inputRegressor.Nrows());
  for(int t=0;t<inputRegressor.Nrows();t++)
    im.value(0,0,0,t) = inputRegressor(t+1);
  im=bandpass_temporal_filter(im, hp_sigma, lp_sigma);
  for(int t=0;t<inputRegressor.Nrows();t++)
    outputRegressor(t+1) = im.value(0,0,0,t);
  return( outputRegressor );
}


float z2t(float z, int dof)
{
  float tmin=0, tmax=1e10, absz=fabs(z), old_diff=1e11;
  
  while (old_diff-(tmax-tmin)>1e-5) {
      old_diff=tmax-tmin;
      
      float t=(tmax+tmin)/2;
      float z=T2z::getInstance().convert(t,dof);
      if (z>absz)
	tmax=t;
      else
	tmin=t;
    }

  if (z<0)
    tmin=-tmin;

  return tmin;
}

// }}}
// {{{ mygammapdf 

ReturnMatrix mygammapdf(const int npts, const float mult, const float delay, const float sigma)
{
  ColumnVector grot(npts);
  
  for(int i=0; i<npts; i++)
    grot(i+1)=i/mult;

  grot = gammapdf(grot.t(),delay,sigma*sigma).t();

  grot.Release();
  return grot;
}

// }}}
// {{{ estimate_X_heights 

ReturnMatrix estimate_X_heights(const Matrix real_X)
{
  ColumnVector real_X_heights(real_X.Ncols());
  for(int real_ev=0; real_ev<real_X.Ncols(); real_ev++) // make sure that 0 is included in the X2 min:max range this is so that (eg) all-1s EVs have height 1 (etc)
    real_X_heights(real_ev+1) = Max( real_X.Column(real_ev+1).Maximum() , 0.0f ) - Min( real_X.Column(real_ev+1).Minimum() , 0.0f );
  real_X_heights.Release();
  return real_X_heights;
}

// }}}
// {{{ feat_svd 
ReturnMatrix feat_svd(const Matrix real_X)
{
  DiagonalMatrix eigenvals(1);
      
  if (real_X.Ncols()>1)
    {
      if ( real_X.Ncols() > real_X.Nrows() )
	SVD(real_X.t(),eigenvals);
      else 
	SVD(real_X,eigenvals);	

      SortDescending(eigenvals);
      float inv_condition = eigenvals.Minimum() / eigenvals.Maximum();

      if ( real_X.Ncols() > real_X.Nrows() ) {
	DiagonalMatrix eigenvalues( real_X.Ncols() );
	eigenvalues=0;
	for ( int i=1; i<= eigenvals.Nrows() ; i++ )
	  eigenvalues(i)=eigenvals(i);
	eigenvals=eigenvalues;

      }
      //cout << "eigenvals = " << eigenvals << endl;
      //cout << "inv_condition = " << inv_condition << endl;

      if (inv_condition<1e-7)
	{
	  cout << "Warning: at least one EV is (close to) a linear combination of the others. You probably need to alter your design.\n(Design matrix is rank deficient - ratio of min:max eigenvalues in SVD of matrix is " << inv_condition << ")\n Contrasts involving these combinations will be set to zero.\n" << endl;

	}

    }

  eigenvals.Release();
  return eigenvals;
}
// }}}
// {{{ renorm kernel 

// void renorm_kernel(ColumnVector &X)
// {
//   double kernelsum=0;

//   for(int i=0; i<X.Nrows(); i++)
//     kernelsum += X(i+1); /* maybe this should be fabs(cX[i]), but looking at double-gamma output, I've not done this */

//   X /= kernelsum;
// }

// }}}
// {{{ orth_i_wrt_j 

// B -= B.A * A / A.A

void orth_i_wrt_j(Matrix &X, int i, int j)
{
  //cout << "orthogonalising real EV " << i << " wrt " << j << endl;

  float thedenom = ( X.Column(j).t() * X.Column(j) ).AsScalar();

  if (thedenom != 0)     // do nothing if Column(j) is constant
    X.Column(i) -= ( X.Column(i).t() * X.Column(j) ).AsScalar() * X.Column(j) / thedenom;
}

// }}}
// {{{ carry out the convolution 

ReturnMatrix do_convolve(const ColumnVector input, const ColumnVector kernel, const int phase, const int renorm)
{
  ColumnVector output(input);
  
  output=0;

  for(int t=0; t<input.Nrows(); t++)
    {
      float kernel_norm=0;
      for(int i=MAX(0,1+t+phase-input.Nrows()); i<MIN(kernel.Nrows(),t+phase+1); i++)
	{
	  output(t+1) += input(t+phase-i+1) * kernel(i+1);
	  kernel_norm += kernel(i+1);
	}
      if (renorm)
	output(t+1) /= kernel_norm;
    }

  output.Release();
  return output;
}

// }}}
// {{{ resample down in time 

/* sample in the MIDDLE of the upsampled period */

void do_resample(const ColumnVector input, Matrix &output, const int real_ev, float trmult, int negpts)
{
  for(int t=0;t<output.Nrows();t++)
    output(t+1,real_ev+1) = input(((int)((t+0.5)*trmult))+negpts+1);
}

// }}}
// {{{ find_line 

/* finds LAST matching entry in setup file */

char *find_line(const string filename, const string key, char *fl)
{
  FILE *fd;
  char *return_ptr=NULL, tmp_fl[10000];
  
  fd=fopen( filename.c_str() , "rb" );
  while ( fgets(tmp_fl, 1000, fd) != NULL ) 
    if (strncmp(tmp_fl,"set ",4)==0)
      for(int j=4; tmp_fl[j]!=0; j++)
	if (strncmp(key.c_str(),tmp_fl+j,key.size())==0)
	  {
	    strcpy(fl,tmp_fl+j);
	    return_ptr = fl+1+key.size();
	    if (fl[strlen(fl)-1]==10) fl[strlen(fl)-1]=0;
	    if (fl[strlen(fl)-1]==13) fl[strlen(fl)-1]=0;
	    if (fl[strlen(fl)-1]=='"') fl[strlen(fl)-1]=0;
	    if (fl[1+key.size()]=='"') return_ptr++;
	  }

  fclose(fd);

  if (return_ptr==NULL)
    {
      cerr << "Can't find key " << key << endl;
      exit(1);
    }
  else
    return return_ptr;
}

string tclKey(const string keyName, const int number1, const int number2, const bool debug=false)
{
  stringstream stream1,stream2;
  stream1 << number1;
  stream2 << number2;
  string str1(stream1.str()), str2(stream2.str());
  if ( debug )
    cerr << "fmri("+keyName+str1+"."+str2+")" << endl;
  return "fmri("+keyName+str1+"."+str2+")";
}

string tclKey(const string keyName, const int number, const bool debug=false)
{
  stringstream stream1;
  stream1 << number;
  string str1(stream1.str());
  if ( debug )
    cerr << "fmri("+keyName+str1+")" << endl;
  return "fmri("+keyName+str1+")";  
}

bool testForZeroContrasts( const Contrasts& contrasts )
{
  for(int con=1; con<=contrasts.nC(); con++) {
    if ( contrasts.C.Column(con).Maximum()==0 && contrasts.C.Column(con).Minimum()==0 ) {
      cerr << "Contrast " << con << " is empty!" << endl;
      return true;
    }
  }
    return false;
}
// }}}
// {{{ setup_font 

/* taken (I think) from pbmtext etc by Jef Poskanzer and George Phillips */

#define DEFAULTFONT_ROWS 155
#define DEFAULTFONT_COLS 112

#define FONT_WIDTH       7
#define FONT_HEIGHT      13

#define FONT_Y_SIZE     30

typedef struct { unsigned char font[DEFAULTFONT_ROWS][DEFAULTFONT_COLS];
                 int           char_row0[95], char_col0[95],
                               char_width, char_height; } FONT_DATA;

void writeCovarianceImage(const string filename, const Contrasts& contrasts, const Matrix& F, const int nFTests, const Matrix& realDesign, const int level, const DiagonalMatrix eigenvals, FONT_DATA *font_data, const ColumnVector& RE );
void writeImagePreview( const string filename, const Contrasts& contrasts, const Matrix& F, const int nFTests, Matrix realDesign, const int level, const Regressors& evs, FONT_DATA *font_data, const vector<string>& titles, const float tr, const float nltffwhm, const int nTimepoints, const vector<int>& G );

void error_exit(const char *outkey )
{
  printf("%s",outkey);
  exit(1);
}

void setup_font(FONT_DATA *font_data)
{
    // {{{ Default Font 

/* The default font, packed in hex so this source file doesn't get huge.
   You can replace this with your own font using pbm_dumpfont().
*/

static unsigned long defaultfont_bits[DEFAULTFONT_ROWS][(DEFAULTFONT_COLS+31)/32] = {
    {0x00000000,0x20000c00,0x10000000,0x00000000},
    {0xc600a000,0x42000810,0x00000002,0x00000063},
    {0x6c00a000,0x45000810,0x00000002,0x00000036},
    {0x6c00a000,0x88800808,0xf2e1dee2,0x00000036},
    {0x54000000,0x80000800,0x11122442,0x0000002a},
    {0x54000001,0x00000800,0x11122442,0x0000002a},
    {0x54000001,0x00000800,0x11122282,0x0000002a},
    {0x44000102,0x00000800,0x11122382,0x00000022},
    {0xee000102,0x00000800,0x11e1e102,0x00000077},
    {0x00000204,0x00000800,0x11002102,0x00000000},
    {0x00000000,0x00000c00,0x11002102,0x00000000},
    {0x00000000,0x003f8000,0xe3807600,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x02000080,0x00040000,0x00120000,0x00000001},
    {0x04000082,0x828e1838,0x20210100,0x00000002},
    {0x04000082,0x82912448,0x20210100,0x00000002},
    {0x08000082,0x8fd01940,0x404087c2,0x00000004},
    {0x08000080,0x050c0622,0x00408102,0x00000004},
    {0x10000080,0x05061874,0x0040828f,0x00008008},
    {0x10000080,0x1f912688,0x00408002,0x00000008},
    {0x20000000,0x0a11098c,0x00408002,0x00000010},
    {0x20000080,0x0a0e0672,0x00210000,0x00000010},
    {0x40000000,0x00040000,0x00210000,0x00000020},
    {0x00000000,0x00000000,0x00120000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x004e0838,0x7023e1cf,0x00008000},
    {0x00000000,0x00913844,0x88620208,0x00008000},
    {0x08000000,0x00910844,0x08a20401,0x00000004},
    {0x10000000,0x01110844,0x08a20401,0x00000008},
    {0x20000000,0x01110808,0x3123c781,0x00000010},
    {0x400003e0,0x02110810,0x0a202441,0x00000020},
    {0x20000000,0x02110820,0x0bf02442,0x00000010},
    {0x10008000,0x04110844,0x88242442,0x00000008},
    {0x08008002,0x040e3e7c,0x7073c382,0x00000004},
    {0x00010000,0x08000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x0000e1c0,0x00000000,0x00000000,0x00000000},
    {0x00011220,0x00000000,0x70e38f87,0x00000000},
    {0x20011220,0x00020020,0x89108448,0x00008010},
    {0x10011220,0x00040010,0x09314448,0x00008008},
    {0x0800e221,0x02083e08,0x11514788,0x00000004},
    {0x040111e0,0x00100004,0x2153e448,0x00000002},
    {0x08011020,0x00083e08,0x213a2448,0x00008004},
    {0x10011040,0x02040010,0x01022448,0x00008008},
    {0x2000e381,0x02020020,0x20e77f87,0x00000010},
    {0x00000000,0x04000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x3803e7ef,0xc73bbe3d,0xdb863ce7,0x0000001c},
    {0x44011224,0x48910808,0x91036648,0x00008022},
    {0x4c011285,0x48910808,0xa1036648,0x00008026},
    {0x54011387,0x081f0808,0xc102a548,0x0000802a},
    {0x54011285,0x09910808,0xe102a548,0x0000802a},
    {0x4e011204,0x08910848,0x9112a4c8,0x00008027},
    {0x40011224,0x08910848,0x891224c8,0x00008020},
    {0x3803e7ef,0x073bbe31,0xcff77e47,0x0000001c},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000003,0x00000000},
    {0x0003e1cf,0x87bff7ef,0xdfbf77c2,0x00000000},
    {0x00013224,0x48a4a244,0x89122442,0x00000000},
    {0x00011224,0x4824a244,0xa8a14482,0x00000000},
    {0x00013227,0x8e04226c,0xa8414102,0x00000000},
    {0x0001e224,0x83842228,0xa8a08102,0x00000000},
    {0x00010224,0x40842228,0xd8a08242,0x00000000},
    {0x00010224,0x48843638,0x51108442,0x00000000},
    {0x0003c1ce,0x6f1f1c10,0x53b9c7c2,0x00000000},
    {0x00000060,0x00000000,0x00000002,0x00000000},
    {0x00000000,0x00000000,0x00000003,0x00000000},
    {0xfe000000,0x00000000,0x00000000,0x0000007f},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00010180,0x000000c0,0x003001c0,0x00000000},
    {0x08008081,0x00040040,0x00100200,0x00000004},
    {0x10008082,0x80040040,0x00100200,0x00000008},
    {0x10004084,0x40023c78,0x70f1c7c7,0x00004008},
    {0x10004080,0x00000244,0x89122208,0x00008008},
    {0x20002080,0x00001e44,0x8113e208,0x00008010},
    {0x10002080,0x00002244,0x81120208,0x00008008},
    {0x10001080,0x00002244,0x89122208,0x00008008},
    {0x10001080,0x00001db8,0x70e9c787,0x00008008},
    {0x10000880,0x00000000,0x00000000,0x00008008},
    {0x08000180,0x00000000,0x00000000,0x00008004},
    {0x00000000,0x1fc00000,0x00000007,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00030080,0x981c0000,0x00000000,0x00000000},
    {0x20010000,0x08040000,0x00000000,0x00000010},
    {0x10010000,0x08040000,0x00000000,0x00000008},
    {0x10016387,0x898474b8,0x72e1d5c7,0x00000008},
    {0x10019080,0x8a042a64,0x89122208,0x00008008},
    {0x08011080,0x8c042a44,0x89122207,0x00000004},
    {0x10011080,0x8a042a44,0x89122200,0x00008008},
    {0x10011080,0x89042a44,0x89122208,0x00008008},
    {0x1003bbe0,0x98dfebe6,0x71e1e787,0x00000008},
    {0x10000000,0x80000000,0x01002000,0x00000008},
    {0x20000000,0x80000000,0x01002000,0x00000010},
    {0x00000007,0x00000000,0x03807000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00008000,0x00000000,0x10410000,0x00000000},
    {0x00008000,0x00000000,0x20408000,0x00000000},
    {0x0001f66e,0xfdfbf77c,0x20408000,0x00000000},
    {0x24008224,0x488a2248,0x20408240,0x00000012},
    {0x54008224,0x4a842210,0x40404540,0x0000002a},
    {0x48008222,0x8a8a1420,0x20408480,0x00000024},
    {0x00008a23,0x85111c44,0x20408000,0x00000000},
    {0x000071d1,0x0531887c,0x20408000,0x00000000},
    {0x00000000,0x00000800,0x20408000,0x00000000},
    {0x00000000,0x00000800,0x10410000,0x00000000},
    {0x00000000,0x00003000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x20000c00,0x10000000,0x00000000},
    {0xc600a000,0x42000810,0x00000002,0x00000063},
    {0x6c00a000,0x45000810,0x00000002,0x00000036},
    {0x6c00a000,0x88800808,0xf2e1dee2,0x00000036},
    {0x54000000,0x80000800,0x11122442,0x0000002a},
    {0x54000001,0x00000800,0x11122442,0x0000002a},
    {0x54000001,0x00000800,0x11122282,0x0000002a},
    {0x44000102,0x00000800,0x11122382,0x00000022},
    {0xee000102,0x00000800,0x11e1e102,0x00000077},
    {0x00000204,0x00000800,0x11002102,0x00000000},
    {0x00000000,0x00000c00,0x11002102,0x00000000},
    {0x00000000,0x003f8000,0xe3807600,0x00000000}
    };

// }}}
    unsigned char b;
    int rows, cols,
        scol, brow, d,
        bcol, ch;
    unsigned long l;

    // {{{ text explanation 

    /*
    ** This routine expects a font bitmap representing the following text:
    **
    ** (0,0)
    **    M ",/^_[`jpqy| M
    **
    **    /  !"#$%&'()*+ /
    **    < ,-./01234567 <
    **    > 89:;<=>?@ABC >
    **    @ DEFGHIJKLMNO @
    **    _ PQRSTUVWXYZ[ _
    **    { \]^_`abcdefg {
    **    } hijklmnopqrs }
    **    ~ tuvwxyz{|}~  ~
    **
    **    M ",/^_[`jpqy| M
    **
    ** The bitmap must be cropped exactly to the edges.
    **
    ** The dissection works by finding the first blank row and column; that
    ** gives the height and width of the maximum-sized character, which is
    ** not too useful.  But the distance from there to the opposite side is
    ** an integral multiple of the cell size, and that's what we need.  Then
    ** it's just a matter of filling in all the coordinates.
    */

// }}}
    // {{{ create font bits 

  for ( rows = 0; rows < DEFAULTFONT_ROWS; ++rows )
  {
    for ( cols = 0; cols < DEFAULTFONT_COLS; cols += 32 )
    {
      l = defaultfont_bits[rows][cols / 32];
      if (cols + 32 < DEFAULTFONT_COLS)
        ch = cols + 32;
      else
        ch = DEFAULTFONT_COLS;
      for ( scol = ch - 1; scol >= cols; --scol )
      {
        if ( l & 1 )
          font_data -> font[rows][scol] = 1;
        else
          font_data -> font[rows][scol] = 0;
        l >>= 1;
      }
    }
  }

// }}}
    // {{{ Find first blank row 

    for ( brow = 0; brow < DEFAULTFONT_ROWS / 6; ++brow )
	{
	b = font_data -> font[brow][0];
	for ( cols = 1; cols < DEFAULTFONT_COLS; ++ cols )
	    if ( font_data -> font[brow][cols] != b )
		goto nextrow;
	goto gotblankrow;
    nextrow: ;
	}

    error_exit("Couldn't find blank row in font\n");

// }}}
    // {{{ gotblankrow: 

gotblankrow:
    /* Find first blank col. */
    for ( bcol = 0; bcol < DEFAULTFONT_COLS / 8; ++bcol )
	{
	b = font_data -> font[0][bcol];
	for ( rows = 1; rows < DEFAULTFONT_ROWS; ++ rows )
	    if ( font_data -> font[rows][bcol] != b )
		goto nextcol;
	goto gotblankcol;
    nextcol: ;
	}
    error_exit("Couldn't find blank col in font.\n");

// }}}
    // {{{ gotblankcol: 

gotblankcol:
    /* Now compute character cell size. */
    d = DEFAULTFONT_ROWS - brow;
    font_data -> char_height = d / 11;
    if ( font_data -> char_height * 11 != d )
	error_exit("problem computing character cell height");
    d = DEFAULTFONT_COLS - bcol;
    font_data -> char_width = d / 15;
    if ( font_data -> char_width * 15 != d )
	error_exit("problem computing character cell width");

    /*printf("height=%d width=%d\n",font_data -> char_height,
           font_data -> char_width);*/

    /* Now fill in the 0,0 coords. */
    rows = font_data -> char_height * 2;
    cols = font_data -> char_width * 2;
    for ( ch = 0; ch < 95; ++ch )
	{
	font_data -> char_row0[ch] = rows;
	font_data -> char_col0[ch] = cols;
	cols += font_data -> char_width;
	if ( cols >= font_data -> char_width * 14 )
	    {
	    cols = font_data -> char_width * 2;
	    rows += font_data -> char_height;
	    }
	}

// }}}
}

// }}}
// {{{ write_string 

void write_string(unsigned char *in, int x, int y, const char *the_string,
		  FONT_DATA *font_data, int colour, int x_size, int y_size)
{
  char string_bit;
  int  i, r, c, X, Y;

  for (i=0; the_string[i] != '\0'; i++)
  {
    string_bit = the_string[i] - ' ';
    for ( r = 0; r < font_data->char_height; ++r )
      for ( c = 0; c < font_data->char_width; ++c )
      {
        X = (i*font_data->char_width)+c+x;
        Y = r+y;
        if ( (X>-1) && (X<x_size) && (Y>-1) && (Y<y_size) &&
             (font_data->font[font_data->char_row0[(int)string_bit] + r][font_data->char_col0[(int)string_bit] + c]==1) )
          in[Y*x_size+X]=colour;
      }
  }
}

void write_string_rgb(unsigned char *r, unsigned char *g, unsigned char *b, 
		      int x, int y, const char *the_string,
		      FONT_DATA *font_data, int cr, int cg, int cb,
		      int x_size, int y_size)
{
  write_string(r, x, y, the_string, font_data, cr, x_size, y_size);
  write_string(g, x, y, the_string, font_data, cg, x_size, y_size);
  write_string(b, x, y, the_string, font_data, cb, x_size, y_size);
}

// }}}
double RequiredEffect(const Matrix& realDesign, const Matrix& Q, const ColumnVector& singleContrast, const Matrix& preWhitening, const float critical_t, double& height, const float noise, const bool usePreWhitening=true)
{
  Matrix X2( realDesign * Q * singleContrast * pinv(singleContrast.t() * Q * singleContrast) );
  if ( usePreWhitening ) // whiten X2
    X2 = preWhitening*X2;
  if ( height == -1 ) // make sure that 0 is included in the X2 min:max range this is so that (eg) all-1s EVs have height 1 (etc)
    height= Max( X2.Maximum() , 0.0f ) - Min( X2.Minimum() , 0.0f );
  float D = height / sqrt( (X2.t() * X2).AsScalar() );
  return ( critical_t * D * noise );
}

int main(int argc, char **argv)
{
  // {{{ variables 
  int    i, t, mnpts, level,
    orig_ev, real_ev,
    shape[10000], nftests, f;
  vector<int> G;
  vector<string> titles;
  float  tr, mult, trmult, nltffwhm=0, maxconvwin=0;
  char   fl[10000], *FSLDIR;
  string fn, filename;
  FONT_DATA *font_data = new FONT_DATA[1];


  Regressors evs;

  // }}}
  // {{{ read arguments and prepare variables 

  if (argc<2) {
    cout << "Usage: feat_model <design_name_root> [confound matrix text file]" << endl;
    exit(1);
  }

  Matrix motionparams(0,0);
  if (argc==3)
    motionparams=remmean(read_ascii_matrix(argv[2]));

  FSLDIR=getenv("FSLDIR");

  fn = string(argv[1])+".fsf";

  level      = atoi(find_line(fn, "fmri(level)", fl));
  tr         = atof(find_line(fn, "fmri(tr)", fl));
  int nTimepoints  = atoi(find_line(fn, "fmri(npts)", fl));
  evs.nTimepoints=nTimepoints;
  int ndelete      = atoi(find_line(fn, "fmri(ndelete)", fl));
  evs.nOriginalNoMotion = atoi(find_line(fn, "fmri(evs_orig)", fl));
  evs.nOriginal=evs.nOriginalNoMotion + (motionparams.Ncols()>0);
  evs.nRealNoMotion = atoi(find_line(fn, "fmri(evs_real)", fl));
  evs.nReal= evs.nRealNoMotion + motionparams.Ncols();
  int nVoxelwiseEVs= atoi(find_line(fn, "fmri(evs_vox)", fl));
  int nContrasts = atoi(find_line(fn, "fmri(ncon_real)", fl));
  nftests    = atoi(find_line(fn, "fmri(nftests_real)", fl));
  titles.resize(evs.nOriginal);
  evs.usesFilter.resize(evs.nOriginal,true);
  evs.usesDeriv.resize(evs.nOriginalNoMotion,false);
  evs.orthogonals.ReSize(evs.nOriginalNoMotion,evs.nOriginalNoMotion+1); //the '0' column (1st) is whether any ortho's are present
  evs.orthogonals=0;

  for(int orig_ev=0; orig_ev<evs.nOriginalNoMotion; orig_ev++) {
      titles[orig_ev]=string(find_line(fn, tclKey("evtitle",orig_ev+1), fl));
      if ( atoi ( find_line( fn, tclKey("deriv_yn",orig_ev+1), fl) ) ) 
	evs.usesDeriv[orig_ev]=true;
      if ( atoi( find_line( fn, tclKey("tempfilt_yn",orig_ev+1), fl ) ) == 0)
	evs.usesFilter[orig_ev]=false;
      evs.orthogonals(orig_ev+1,1) = atoi( find_line( fn, tclKey("ortho",orig_ev+1,0), fl ) ); 
      for(int tmpEV=0; tmpEV<evs.nOriginalNoMotion; tmpEV++)
	if (tmpEV!=orig_ev) 
	  evs.orthogonals(orig_ev+1,tmpEV+2) = atoi( find_line( fn, tclKey("ortho",orig_ev+1,tmpEV+1), fl) );
      
  }

  nTimepoints    -= ndelete;     /* subtract off images which have already been deleted */
  mult     = 1/DT;        /* points in high-res model per second */
  trmult   = tr*mult;     /* points in high-res model per TR */
  mnpts    = (int)(((float)(nTimepoints+10))*trmult+NEGSECS*mult); /* number of points in high-res model - add a few on for safety */

  Matrix originalHighResDesign(mnpts,evs.nOriginal);
  originalHighResDesign=0;

  ColumnVector orig_level(evs.nOriginal); // records the difference between min and first point in an originalHighResDesign column

  Matrix realDesign(nTimepoints,evs.nReal);
  realDesign=0;

  Contrasts contrasts(evs.nReal,nContrasts);
  Matrix F(nftests,contrasts.nC());

  evs.nRealPerOrig.resize(evs.nOriginal,1);

  setup_font(font_data);

  vector<float> triggers(evs.nOriginal*2*10*nTimepoints,-1e7);/* factor of 10 for safety */
  vector<int> negpts(evs.nOriginal,0);
  vector<int> convolve_interaction(evs.nOriginal,0);
  evs.basisorth.resize(evs.nOriginal,0);

  float critical_t = z2t( atof(find_line(fn, "fmri(critical_z)", fl)) , MAX(nTimepoints-evs.nReal,1) );
  float noise      = atof(find_line(fn, "fmri(noise)", fl));
  float noisear    = atof(find_line(fn, "fmri(noisear)", fl));

  /* setup prewhitening stuff */
  SymmetricMatrix pwV(nTimepoints);
  for(int i=1; i<=nTimepoints; i++)
    for(int j=1; j<=i; j++)
      pwV(i,j)=pow(noisear,(double)abs(i-j));
  Matrix pwA = (Cholesky(pwV)).i(); // note newmat definition opposite of matlab
  //cout << pwA * pwV * pwA.t() << endl; // should be identity

  // }}}
  // {{{ read contrasts 
  

  // modified from find_line
  FILE *fd;
  char tmp_fl[10000];
  fd=fopen(fn.c_str(),"rb");
  while ( fgets(tmp_fl, 1000, fd) != NULL ) {
    if (strncmp(tmp_fl,"set fmri(conname_real.",22)==0)
      {
	//tokenize tmp_vals...
	int con=atoi(strtok(tmp_fl+22,") "));
	if ( con <= nContrasts ) {
	  contrasts.name[con-1]=string(strtok(NULL,"\n\r\f"));
	  size_t found( contrasts.name[con-1].find_first_of('\"') );
	  if ( found != string::npos ) //These lines are to match the output of the previous feat_model ( spacing the png and removing quotes etc )
	    contrasts.name[con-1].erase(contrasts.name[con-1].begin(),contrasts.name[con-1].begin()+found+1);
	  found = contrasts.name[con-1].find_last_of('\"');
	  if ( found != string::npos )
	    contrasts.name[con-1].replace(found,1,1,' ');
	}
      }
    if (strncmp(tmp_fl,"set fmri(con_real",17)==0)
      {
	//tokenize tmp_vals...
	int con=atoi(strtok(tmp_fl+17,"."));// before period is contrastnum
	if ( con <= nContrasts ) {
	  real_ev=atoi(strtok(NULL,") "));//between period and ") " is real_ev num
	  contrasts.C(real_ev,con)=atof(strtok(NULL," \0\n\t")); // rest is value
	}
      }
    if (nftests>0)
      if (strncmp(tmp_fl,"set fmri(ftest_real",19)==0)
	{
	  //tokenize tmp_vals...
	  f=atoi(strtok(tmp_fl+19,"."));// before period is f-testnum
	  int con=atoi(strtok(NULL,") "));//between period and ") " is contrast num
	  if ( con <= nContrasts ) 
	    F(f,con)=atof(strtok(NULL," \0\n\t")); // rest is value
	}
  }
  fclose(fd);

  if ( testForZeroContrasts(contrasts) )
    exit(1);

  // }}}
  // {{{ read and create design matrix 

  if (level==1) {

      bool temphp_yn( atoi(find_line(fn, "fmri(temphp_yn)", fl)));
      bool templp_yn( atoi(find_line(fn, "fmri(templp_yn)", fl)));
      double hp_sigma(-1), lp_sigma(-1);
      nltffwhm  = atof(find_line(fn, "fmri(paradigm_hp)", fl)); 
      if ( templp_yn ) 
	lp_sigma=2.8/tr;
      if ( temphp_yn ) 
	hp_sigma=0.5*nltffwhm/tr;
      // {{{ create basic shape 

      for(orig_ev=0; orig_ev<evs.nOriginalNoMotion; orig_ev++)
	{
	  evs.nRealPerOrig[orig_ev]=1;

	  shape[orig_ev]=atoi(find_line(fn, tclKey("shape",orig_ev+1), fl));

	  switch (shape[orig_ev])
	    {
	    case 0:
	      // {{{ square wave 

	      {
		float skip=atof(find_line(fn, tclKey("skip",orig_ev+1), fl))*mult;
		float off=atof(find_line(fn, tclKey("off",orig_ev+1), fl))*mult;
		float on=atof(find_line(fn, tclKey("on",orig_ev+1), fl))*mult;
		float phase=atof(find_line(fn, tclKey("phase",orig_ev+1), fl))*mult;
		float stop=atof(find_line(fn, tclKey("stop",orig_ev+1), fl))*mult;

		if ( (stop<0) || (stop+skip>mnpts) ) stop=mnpts-skip;
  
		for(t=(int)skip;t<(int)(skip+stop);t++)
		  { // do modulo maths in float not int - necessary for very short TR and block length
		    if ( t+phase-skip - ((int)((t+phase-skip)/(off+on)))*(off+on) >= off )
		      originalHighResDesign(t+1,orig_ev+1)=1.0;
		    else
		      originalHighResDesign(t+1,orig_ev+1)=0.0;
		  }
	      }
	      break;

	    case 1:
	      // {{{ sinusoid 

	      { 
		int skip=(int)(atof(find_line(fn, tclKey("skip",orig_ev+1), fl))*mult); 
		float period=atof(find_line(fn, tclKey("period",orig_ev+1), fl)) * mult * 0.5 / M_PI;
		float phase=atof(find_line(fn, tclKey("phase",orig_ev+1), fl)) * mult;
		int stop=(int)(atof(find_line(fn, tclKey("stop",orig_ev+1), fl))*mult);

		if ( (stop<0) || (stop+skip>mnpts) ) stop=mnpts-skip;

		for(t=skip;t<skip+stop;t++)
		  originalHighResDesign(t+1,orig_ev+1) = 0.5 * ( 1.0 + sin( M_PI + (t+phase-skip) / period ) );
	      }
	      break;

	      // }}}
	    case 2:
	      // {{{ custom single column 

	      {
		FILE *ifp2;

		filename=string(find_line(fn, tclKey("custom",orig_ev+1) , fl));

		if ((ifp2=fopen(filename.c_str(),"rb"))==NULL)
		  {
		    cout << "Can't open " << filename << " for reading" << endl;
		    exit(1);
		  }

		for(i=0,t=0; t<nTimepoints; t++)
		  {
		    float tmpf;
	  
		    if (fscanf(ifp2,"%f",&tmpf)!=1)
		      {
			cout << "Not enough data in " << filename << endl;
			exit(1);
		      }

		    for(; i<trmult*(t+1); i++)
		      originalHighResDesign(i+1,orig_ev+1)=tmpf;
		  }

		fclose(ifp2);
	      }
	      break;

	      // }}}
	    case 3:
	      // {{{ custom 3 columns 

	      {
		FILE *ifp2;
		float start, stop, value;
		int   success=0;

		filename=string(find_line(fn, tclKey("custom",orig_ev+1) , fl));

		if ((ifp2=fopen(filename.c_str(),"rb"))==NULL)
		  {
		    cout << "Can't open " << filename << " for reading" << endl;
		    exit(1);
		  }

		negpts[orig_ev] = (int)(NEGSECS*mult);

		while( fscanf(ifp2,"%f %f %f",&start,&stop,&value) == 3 )
		  {
		    start=start*mult+negpts[orig_ev]; stop=start+stop*mult;

		    if ( (stop>0) && (stop-start<1) ) stop=start+1.1; /* if very short stim time make sure it happens */

		    if (start<0)     start=0;    /* don't enter stuff before t=-negpts */
		    if (start>mnpts) stop=0;     /* don't do anything */
		    if (stop>mnpts)  stop=mnpts; /* don't overrun */

		    for(t=(int)start; t<(int)stop; t++)
		      {
			originalHighResDesign(t+1,orig_ev+1)=value;
			success=1;
		      }
		  }

		fclose(ifp2);

		if ( success==0 )
		  {
		    cout << "No valid\n[onset duration strength]\ntriplets found in " << filename << endl;
		    exit(1);
		  }

	      }
	      break;

	      // }}}
	    case 4:
	      // {{{ interactions

	      // complication; if the EVs used to feed into the interaction EV don't
	      // have the same convolution settings, then we cannot do this
	      // interaction before convolving the interaction - we must do the convolutions first....

	      {
		originalHighResDesign.Column(orig_ev+1)=1;
		int tmp, convEV=-1, other_convolves=-1;
		for(tmp=0; tmp<orig_ev; tmp++)
		  {
		    if ( atoi(find_line(fn, tclKey("interactions",orig_ev+1,tmp+1), fl)) )
		      {
			float subtraction = originalHighResDesign(1,tmp+1) - orig_level(tmp+1);

                        int convolve=atoi(find_line(fn, tclKey("convolve",tmp+1), fl)); convEV=tmp;

			if (other_convolves==-1)
			  other_convolves=convolve;
			else
			  if ( other_convolves != convolve )
			    other_convolves=0;

			ColumnVector tmpv(mnpts);
			tmpv=0;
			tmpv.Rows(1, mnpts-negpts[tmp]) = originalHighResDesign.SubMatrix(1+negpts[tmp], mnpts, tmp+1, tmp+1);

			int zero_how = atoi(find_line(fn, tclKey("interactionsd",orig_ev+1,tmp+1), fl));
			if ( zero_how == 1 )
			  subtraction = ( originalHighResDesign.Column(tmp+1).Maximum() + originalHighResDesign.Column(tmp+1).Minimum() ) / 2.0;
			else if ( zero_how == 2 )
			  subtraction = originalHighResDesign.Column(tmp+1).Sum()/originalHighResDesign.Nrows();

			// cout << "subtraction=" << subtraction <<endl;
			originalHighResDesign.Column(orig_ev+1) = SP( originalHighResDesign.Column(orig_ev+1) , tmpv - subtraction );
		      }
		  }

		if (other_convolves>0)
		  convolve_interaction[orig_ev]=convEV;
		//cout << orig_ev << " " << convolve_interaction[orig_ev] << endl;
	      }
	      break;
	    case 9:
	      // {{{ voxelwise input

	      {
		volume4D<float> ev_image;
		if ( read_volume4D(ev_image,find_line(fn, tclKey("evs_vox_",orig_ev+1), fl)) )
		  cout << "Warning: voxelwise EV " << orig_ev+1 << " isn't readable" << endl;
		for(i=0,t=0; t<nTimepoints; t++)
		  {
		    float tmpf(ev_image[t].mean());
		    for(; i<trmult*(t+1); i++)
		      originalHighResDesign(i+1,orig_ev+1)=tmpf;
		  }
		save_volume4D(ev_image,string(argv[1])+"InputVoxelwiseEV"+num2str(orig_ev+1));
		ev_image=bandpass_temporal_filter(ev_image, hp_sigma, lp_sigma);
		save_volume4D(ev_image,string(argv[1])+"VoxelwiseEV"+num2str(orig_ev+1));
	      }
	      // }}}
	    }

	  //cout << originalHighResDesign << endl;

	  orig_level(orig_ev+1) = originalHighResDesign(1,orig_ev+1) - originalHighResDesign.Column(orig_ev+1).Minimum();
	  //cout << "EV" << orig_ev+1 << " min=" << originalHighResDesign.Column(orig_ev+1).Minimum() << " t[0]=" << originalHighResDesign(1,orig_ev+1) << " level=" << orig_level(orig_ev+1) << endl;

	  if (shape[orig_ev]>1) // demean this orig EV
	    originalHighResDesign.Column(orig_ev+1)=remmean(originalHighResDesign.Column(orig_ev+1));

	  // {{{ triggers 

#define TRIGGER_THRESH 0.001

	  {
	    float previous=0;
	    int trigger_count=0;
  
	    for(t=0;t<nTimepoints*trmult;t++) {
		if ( (originalHighResDesign(t+negpts[orig_ev]+1,orig_ev+1)>TRIGGER_THRESH) && (previous<=TRIGGER_THRESH) && (trigger_count<nTimepoints*2) )
		  triggers[trigger_count++*evs.nOriginal+orig_ev]=t/trmult - 0.5; /* the 0.5 because sampling is halfway through TR */
		if ( (originalHighResDesign(t+negpts[orig_ev]+1,orig_ev+1)<=TRIGGER_THRESH) && (previous>TRIGGER_THRESH) && (trigger_count<nTimepoints*2) )
		  triggers[trigger_count++*evs.nOriginal+orig_ev]=t/trmult - 0.5; /* the 0.5 because sampling is halfway through TR */
		previous=originalHighResDesign(t+negpts[orig_ev]+1,orig_ev+1);
	      }

	  }

	  // }}}
	}
      //Find first minimum timepoint for perfusion subtraction
      Matrix downsampledOriginalModel(nTimepoints,evs.nOriginal);
      for(orig_ev=0; orig_ev<evs.nOriginalNoMotion; orig_ev++)
	do_resample(originalHighResDesign.Column(orig_ev+1),downsampledOriginalModel,orig_ev, trmult, negpts[orig_ev]);
      int minimumTimepoint=1;
      for (int row=1;row<=downsampledOriginalModel.Nrows();row++)
	{
	  if ( downsampledOriginalModel.Row(row).Sum() < downsampledOriginalModel.Row(minimumTimepoint).Sum() )
	    minimumTimepoint=row;
	}  
      ofstream outputFile((string(argv[1])+".min").c_str());
      if(outputFile.is_open()) {
	outputFile << minimumTimepoint;
	outputFile.close();
      }
      // }}}
      // {{{ convolve and resample down in time 

      for(orig_ev=real_ev=0; orig_ev<evs.nOriginalNoMotion; orig_ev++, real_ev++)
	if (shape[orig_ev]<10) {
	    int convolve=atoi(find_line(fn, tclKey("convolve",orig_ev+1), fl));
	    int c_orig_ev=orig_ev; // normally the same as orig_ev except for interactions
	    if (shape[orig_ev] == 4) {// i.e. interactions
	      if (convolve_interaction[orig_ev]==0)
		convolve=0;
	      else {
		c_orig_ev=convolve_interaction[orig_ev];
		convolve=atoi(find_line(fn, tclKey("convolve",c_orig_ev+1), fl));
	      }
	    }
	    if (convolve>0) {
	      int convolve_phase=(int)(atof(find_line(fn, tclKey("convolve_phase",c_orig_ev+1), fl))*mult);

	      if ( (convolve>3) && (convolve<8) ) {
		evs.basisorth[orig_ev]=atoi(find_line(fn, tclKey("basisorth",c_orig_ev+1), fl));
	      }

	      switch (convolve) {
	         case 1: // {{{ Gaussian 
		 {
		   float sigma=atof(find_line(fn, tclKey("gausssigma",c_orig_ev+1), fl))*mult;
		   float delay=atof(find_line(fn, tclKey("gaussdelay",c_orig_ev+1), fl))*mult;
		   
		   int fw = (int)(delay + sigma*5);
		   maxconvwin=MAX(fw,maxconvwin);
		   ColumnVector cX(fw);
		   
		   for(i=0; i<fw; i++) {
		     float tmpf=(((float)i)-delay)/sigma;
		     tmpf=exp(-0.5*tmpf*tmpf);
		     if (i<1/DT) tmpf *= exp( -pow( 1-0.5*(((float)i)*DT) , 6.0 ) ); /* between 0 and 1s, window start of kernel */
		     cX(i+1)=tmpf;
		   }

		   ColumnVector oX = do_convolve(originalHighResDesign.Column(orig_ev+1),cX,convolve_phase,1);
		   do_resample(oX, realDesign, real_ev, trmult, negpts[orig_ev]);
		 }
		 break;
	         case 2: // {{{ Gamma 
		 {
		   float sigma=atof(find_line(fn, tclKey("gammasigma",c_orig_ev+1), fl));
		   float delay=atof(find_line(fn, tclKey("gammadelay",c_orig_ev+1), fl));

		   int fw = (int)((delay + sigma*5)*mult);
		   maxconvwin=MAX(fw,maxconvwin);

		   ColumnVector cX = mygammapdf(fw,mult,delay,sigma);

		   ColumnVector oX=do_convolve(originalHighResDesign.Column(orig_ev+1),cX,convolve_phase,1);
		   do_resample(oX, realDesign, real_ev, trmult, negpts[orig_ev]);
		 }
		 break;
		 case 3: // {{{ double-gamma HRF 
		 {
		   float sigma1=2.449, delay1=6, // first gamma
		     sigma2=4, delay2=16,    // second gamma
		     ratio=6;                // hrf = gammapdf1 - gammapdf2/ratio;
		   
		   int fw = (int)((delay2 + sigma2*5)*mult);
		   maxconvwin=MAX(fw,maxconvwin);
		   
		   ColumnVector cX = mygammapdf(fw,mult,delay1,sigma1) - mygammapdf(fw,mult,delay2,sigma2)/ratio;
		   
		   ColumnVector oX=do_convolve(originalHighResDesign.Column(orig_ev+1),cX,convolve_phase,1);
		   do_resample(oX, realDesign, real_ev, trmult, negpts[orig_ev]);
		 }
		 break;
		 case 4: // {{{ Gamma basis functions 
                 #define WINDOW_FRAC 0.25
		 {
		   int fnumber=atoi(find_line(fn, tclKey("basisfnum",c_orig_ev+1), fl));	
		   int window=(int)(atof(find_line(fn, tclKey("basisfwidth",c_orig_ev+1), fl))*mult);
		   evs.nRealPerOrig[orig_ev]+=fnumber-1;
		   
		   int fw = window;
		   maxconvwin=MAX(fw,maxconvwin);
		   
		   float delay = fw/(2*mult);

		   for(int fnum=1; fnum<=fnumber; fnum++) {
		     float sigma = delay/2;
		     ColumnVector cX = mygammapdf(fw,mult,delay,sigma);
		     for(i=(int)((1-WINDOW_FRAC)*fw); i<fw; i++)
		       cX(i+1) *= 1 - (i-(1-WINDOW_FRAC)*fw)/(WINDOW_FRAC*fw);

		     ColumnVector oX=do_convolve(originalHighResDesign.Column(orig_ev+1),cX,convolve_phase,1);
		     do_resample(oX, realDesign, real_ev, trmult, negpts[orig_ev]);

		     delay *= 0.5;
		     real_ev++;
		   }

		   real_ev--;
		 }
		 break;
		 case 5: // {{{ Sine basis functions 
		 {
		   int fnumber=atoi(find_line(fn, tclKey("basisfnum",c_orig_ev+1), fl));
		   int window=(int)(atof(find_line(fn, tclKey("basisfwidth",c_orig_ev+1), fl))*mult);
		   evs.nRealPerOrig[orig_ev]+=fnumber-1;
		   int fw = window;
		   maxconvwin=MAX(fw,maxconvwin);
		   ColumnVector cX(fw);

		   for(int fnum=1; fnum<=fnumber; fnum++) {
		     for(i=0; i<fw; i++) {
		       if (i<0) cout << "gcc is broken" << endl; // bizzarely, without this pointless line, the loop crashes....
		       cX(i+1)= sin(M_PI*i*fnum/fw);
		     }

		     ColumnVector oX=do_convolve(originalHighResDesign.Column(orig_ev+1),cX,convolve_phase,0);
		     do_resample(oX, realDesign, real_ev, trmult, negpts[orig_ev]);
		     
		     real_ev++;
		   }

		   real_ev--;
		 }
		 break;
		 case 6: // {{{ FIR basis functions 
		 {
		   int fnumber=atoi(find_line(fn, tclKey("basisfnum",c_orig_ev+1), fl));
		   int window=(int)(atof(find_line(fn, tclKey("basisfwidth",c_orig_ev+1), fl))*mult);
		   evs.nRealPerOrig[orig_ev]+=fnumber-1;
		   int fw = window;
		   maxconvwin=MAX(fw,maxconvwin);
		   ColumnVector cX(fw);

		   for(int fnum=0; fnum<fnumber; fnum++) {
		     for(i=0; i<fw; i++)
		       if ( (fw*fnum/fnumber<=i) && (i<fw*(fnum+1)/fnumber) )
			 cX(i+1)=1;
		       else
			 cX(i+1)=0;

		     ColumnVector oX=do_convolve(originalHighResDesign.Column(orig_ev+1),cX,convolve_phase,0);
		     do_resample(oX, realDesign, real_ev, trmult, negpts[orig_ev]);
		     real_ev++;
		   }

		   real_ev--;
		 }
		 break;
		 case 7: // {{{ Custom basis functions 
		 {
		   int fnumber=atoi(find_line(fn, tclKey("basisfnum",c_orig_ev+1), fl));
		   string bfcustomname(string(find_line(fn, tclKey("bfcustom",c_orig_ev+1), fl)));
		   Matrix icX(mnpts,fnumber);
		   evs.nRealPerOrig[orig_ev]+=fnumber-1;
 
		   FILE *ifp2;
		   if ((ifp2=fopen(bfcustomname.c_str(),"rb"))==NULL) {
		     cout << "Can't open " << bfcustomname << " for reading" << endl;
		     exit(1);
		   }

		   int fw, carryon=1;
		   for(fw=0; carryon; fw++)
		     for(int fnum=0; fnum<fnumber; fnum++) {
		       float tmpval;
		       if( fscanf(ifp2,"%f",&tmpval) != 1 )
			 carryon=0;
		       icX(fw+1,fnum+1)=tmpval;
		     }

		   fw--;
		   fclose(ifp2);
		   maxconvwin=MAX(fw,maxconvwin);
  
		   for(int fnum=0; fnum<fnumber; fnum++) {
		     ColumnVector cX=icX.SubMatrix(1,fw,fnum+1,fnum+1);
		     ColumnVector oX=do_convolve(originalHighResDesign.Column(orig_ev+1),cX,convolve_phase,0);
		     do_resample(oX, realDesign, real_ev, trmult, negpts[orig_ev]);

		     real_ev++;
		   }

		   real_ev--;
		 }
		 break;
	      }
	    } else {
	      do_resample(originalHighResDesign.Column(orig_ev+1), realDesign, real_ev, trmult, negpts[orig_ev]);
	      if ( shape[orig_ev]==1 )
		  // {{{ create sinusoid harmonics 

		  /* this is treated like a convolution to keep the structure simple -
		     thus the file key value searching is just a duplicate of that from the
		     principal sinusoid EV already created */

		  {
		    int nharmonics;

		    nharmonics=atoi(find_line(fn, tclKey("nharmonics",c_orig_ev+1), fl));

		    if (nharmonics>0)
		      {
			int skip=(int)(atof(find_line(fn, tclKey("skip",c_orig_ev+1), fl))*mult);
			float period=atof(find_line(fn, tclKey("period",c_orig_ev+1), fl)) * mult * 0.5 / M_PI;
			float phase=atof(find_line(fn, tclKey("phase",c_orig_ev+1), fl)) * mult;
			int stop=(int)(atof(find_line(fn, tclKey("stop",c_orig_ev+1), fl))*mult);

			if ( (stop<0) || (stop+skip>mnpts) ) stop=mnpts-skip;

			evs.nRealPerOrig[orig_ev]+=nharmonics;

			for(int harm=1; harm<=nharmonics; harm++)
			  {
			    real_ev++;

			    for(t=skip;t<skip+stop;t++)
			      originalHighResDesign(t+1,orig_ev+1) = 0.5 * ( 1.0 + sin( M_PI + (harm+1)*(t+phase-skip) / period ) );
	
			    do_resample(originalHighResDesign.Column(orig_ev+1), realDesign, real_ev, trmult, negpts[orig_ev]);
			  }
		      }
		  }

		// }}}
	      }
	  }

      // }}}
      // {{{ recomputed interactions

      for(orig_ev=real_ev=0; orig_ev<evs.nOriginalNoMotion; real_ev+=evs.nRealPerOrig[orig_ev], orig_ev++)
	if ( (shape[orig_ev]==4) && (convolve_interaction[orig_ev]==0) ) {
	  realDesign.Column(real_ev+1)=1;
	  for(int tmp(0),tmp_real(0); tmp<orig_ev; tmp_real+=evs.nRealPerOrig[tmp], tmp++) {
	    if ( atoi(find_line(fn,  tclKey("interactions",orig_ev+1,tmp+1), fl)) ) {
		float subtraction = originalHighResDesign(1,tmp+1) - orig_level(tmp+1);
		int zero_how = atoi(find_line(fn, tclKey("interactionsd",orig_ev+1,tmp+1), fl));
		if ( zero_how == 1 )
		  subtraction = ( realDesign.Column(tmp_real+1).Maximum() + realDesign.Column(tmp_real+1).Minimum() ) / 2.0;
		else if ( zero_how == 2 )
		  subtraction = realDesign.Column(tmp_real+1).Sum()/realDesign.Nrows();
		realDesign.Column(real_ev+1) = SP( realDesign.Column(real_ev+1) , realDesign.Column(tmp_real+1) - subtraction );
	    }
	  }  
	}
      // }}}
      // {{{ add motion params to model if required 

      if ( motionparams.Ncols() > 0 ) {
	for(orig_ev=real_ev=0; orig_ev<evs.nOriginalNoMotion; real_ev+=evs.nRealPerOrig[orig_ev], orig_ev++)
	  ; //Fast forward through normal EVs
	  //cout << "inserting motion parameters starting at real EV " << real_ev+1 << endl;
	  realDesign.Columns(real_ev+1,real_ev+motionparams.Ncols()) = motionparams;
	  evs.nRealPerOrig[evs.nOriginal-1] = motionparams.Ncols();
	  shape[evs.nOriginal-1]=2;
      }

      // two passes through the rest of the design setup:
      // first pass runs without temporal filtering, in order to get accurate peak-peak EV height estimation
      // second pass creates the final design, including temporal filtering

 
      vector<int> tempNReal= evs.nRealPerOrig;

      
      evs.copyOfRealDesign=realDesign;
      evs.nRealPerOrig=tempNReal;
      peakAndFilter(true, evs, contrasts, -1, -1, nTimepoints, pwA,  critical_t, noise, string(argv[1])+".frf");
      /*for(double mod=90.0/(2.0*tr);mod<=3*nTimepoints;mod+=1) {
	evs.copyOfRealDesign=realDesign;
	evs.nRealPerOrig=tempNReal;
	peakAndFilter(false, evs, contrasts, mod, -1, nTimepoints, pwA,  critical_t, noise, string(argv[1])+".frf");
	for(int con=1; con<=contrasts.nC(); con++) 
	  cerr << contrasts.RE(con) << endl;
	cerr << "****************" << endl;
	}*/
      evs.copyOfRealDesign=realDesign;
      evs.nRealPerOrig=tempNReal;
      peakAndFilter(false, evs, contrasts, hp_sigma, lp_sigma, nTimepoints, pwA,  critical_t, noise, string(argv[1])+".frf");
      realDesign=evs.copyOfRealDesign;

      //NEW
  vector<string> fileList;
  ifstream inputFile("vef.dat");
  if ( inputFile.is_open() ) {
    string input;
    while (getline(inputFile,input,',')) {
      if ( input.size() != 0 )
	fileList.push_back(input);
    }
    inputFile.close();
  }

  for( unsigned int currentRegressor=0; currentRegressor < fileList.size(); currentRegressor++ )
  {
    volume4D<float> inputImage;
    read_volume4D(inputImage,fileList[currentRegressor]);
    save_volume4D(inputImage,"Input"+fileList[currentRegressor]);
    inputImage=bandpass_temporal_filter(inputImage, hp_sigma, lp_sigma);
    save_volume4D(inputImage,fileList[currentRegressor]);
  }
  //END NEW
    } else { // group level design 

    for(real_ev=0; real_ev<evs.nReal; real_ev++) {
      if (real_ev<evs.nReal-nVoxelwiseEVs)
	for(t=0;t<nTimepoints;t++) {
	  realDesign(t+1,real_ev+1)=atof(find_line(fn, tclKey("evg",t+1,real_ev+1), fl));
	} else {
        volume4D<float> ev_image;
	if ( read_volume4D(ev_image,find_line(fn, tclKey("evs_vox_",1+real_ev-(evs.nReal-nVoxelwiseEVs)), fl)) )
	  cout << "Warning: voxelwise EV " << 1+real_ev-nVoxelwiseEVs << " isn't readable" << endl;
	for(t=0;t<nTimepoints;t++)
	  realDesign(t+1,real_ev+1)=ev_image[t].mean();
      }
    }
    // {{{ orthogonalisation 

    for(real_ev=0; real_ev<evs.nReal-nVoxelwiseEVs; real_ev++)
      for(int tmp_real_ev=0; tmp_real_ev<evs.nReal-nVoxelwiseEVs; tmp_real_ev++)
	if (tmp_real_ev!=real_ev)
	  {
	    if (atoi(find_line(fn, tclKey("ortho",real_ev+1,tmp_real_ev+1), fl)))
	      orth_i_wrt_j(realDesign,real_ev+1,tmp_real_ev+1);
	  }

    // }}}
    evs.realHeights = estimate_X_heights(realDesign);
    // {{{ check rank of DM and do efficiency test 

    // first do real rank deficiency test
    evs.eigenvals = feat_svd(realDesign);
    //eigenvals = feat_svd(realDesign.Columns(1,evs.nReal-nVoxelwiseEVs));

    // now do "meaningful" rank deficiency test
    Matrix Q = pinv(realDesign.t() * realDesign);
    for(int con=1; con<=contrasts.nC(); con++) {
      double height(-1);
      contrasts.RE(con) = RequiredEffect(realDesign, Q, contrasts.C.Column(con), IdentityMatrix(1), critical_t, height, noise, false);
      contrasts.realHeights(con)=height;
     }

 
  }

  writeDesignMatrix(string(argv[1])+".mat", evs.realHeights, realDesign);

  if (level==1)
    writeTriggers(string(argv[1])+".trg",triggers,shape,evs,maxconvwin/trmult);
   
  writeContrastMatrix(string(argv[1])+".con", contrasts, contrasts.RE );
  writeFMatrix(string(argv[1])+".fts", contrasts, F, nftests , realDesign);

  if (level==2) {
    G.resize(nTimepoints);
      for(t=0;t<nTimepoints;t++)
	  G[t]=atoi(find_line(fn, tclKey("groupmem.",t+1), fl));
    writeGroupFile(string(argv[1])+".grp",G,nTimepoints,realDesign);
  }

  writeCovarianceImage(string(argv[1])+"_cov.ppm", contrasts, F, nftests, realDesign, level, evs.eigenvals, font_data, contrasts.RE);
  writeImagePreview(string(argv[1])+".ppm", contrasts, F, nftests, realDesign, level, evs, font_data, titles, tr, nltffwhm, nTimepoints, G);

  filename=string(getenv("FSLDIR"))+"/bin/wpng -q -overwrite  "+string(argv[1])+".ppm ";
  system(filename.c_str());

  return(0);
}

Regressors peakAndFilter(const bool estimateRealHeights,Regressors& evs, Contrasts& contrasts, const double hp_sigma, const double lp_sigma, const int nTimepoints, const Matrix& preWhitening, const float critical_t, const float noise, const string filename)
{
  if ( !estimateRealHeights ) {
    //apply temporal filtering
    for(int orig_ev(0), real_ev(0); orig_ev< evs.nOriginal && ( ( hp_sigma != -1 ) || ( lp_sigma != -1 ) ); orig_ev++) 
      for(int i=0; i<evs.nRealPerOrig[orig_ev]; i++, real_ev++) //Increment real EV to start of new orig ev
	if ( evs.usesFilter[orig_ev] ) evs.copyOfRealDesign.Column(real_ev+1)=filterRegressor(evs.copyOfRealDesign.Column(real_ev+1),hp_sigma,lp_sigma);   
  } //End of second pass temporal filtering
  evs.copyOfRealDesign=remmean(evs.copyOfRealDesign);
  // basis function orthogonalisation (within orig_ev)
  for(int orig_ev(0), real_ev(0); orig_ev<evs.nOriginalNoMotion; real_ev+=evs.nRealPerOrig[orig_ev], orig_ev++)
    if ( evs.basisorth[orig_ev] ) {
      for(int tmp=1; tmp<evs.nRealPerOrig[orig_ev]; tmp++) {
	Matrix tmp_orth = evs.copyOfRealDesign.Columns(real_ev+1,real_ev+tmp);
	evs.copyOfRealDesign.Column(real_ev+tmp+1) = evs.copyOfRealDesign.Column(real_ev+tmp+1) - tmp_orth*(pinv(tmp_orth)*evs.copyOfRealDesign.Column(real_ev+tmp+1));
      }
    }
  // main orthogonalisation
  for(int orig_ev(0), real_ev(0); orig_ev<evs.nOriginalNoMotion; real_ev+=evs.nRealPerOrig[orig_ev], orig_ev++) {
    if ( evs.orthogonals(orig_ev+1,1) ) {
      Matrix tmp_orth(nTimepoints,0);
      for(int tmp_orig_ev(0), tmp_real_ev(0); tmp_orig_ev<evs.nOriginalNoMotion; tmp_real_ev+=evs.nRealPerOrig[tmp_orig_ev], tmp_orig_ev++)
	if (tmp_orig_ev!=orig_ev) {
	  if ( evs.orthogonals(orig_ev+1,tmp_orig_ev+2) )  // i.e. should we orth orig_EV(orig_ev) wrt orig_EV(tmp_orig_ev)
	    for(int tmp=0; tmp<evs.nRealPerOrig[tmp_orig_ev]; tmp++)
	      tmp_orth = tmp_orth | evs.copyOfRealDesign.Column(tmp_real_ev+tmp+1);
	}
      if (tmp_orth.Ncols()>0)
	for(int tmp=0; tmp<evs.nRealPerOrig[orig_ev]; tmp++)  // ORTH evs.copyOfRealDesign.Column(real_ev+tmp+1) WRT tmp_orth		
	  evs.copyOfRealDesign.Column(real_ev+tmp+1) = evs.copyOfRealDesign.Column(real_ev+tmp+1) - tmp_orth*(pinv(tmp_orth)*evs.copyOfRealDesign.Column(real_ev+tmp+1));
    }
  }
  for(int orig_ev(0),real_ev(0); orig_ev<evs.nOriginalNoMotion; orig_ev++) {
    if ( evs.usesDeriv[orig_ev] ) {
      evs.nRealPerOrig[orig_ev]++;
      // shift columns to the right to make space for tempderiv (if this isn't the right-most EV)
      if (orig_ev<evs.nOriginal-1) {
	//cout << "shifting real EVs " << real_ev+2 << ":" << evs.nReal-1 << " to " << real_ev+3 << ":" << evs.nReal << endl;
	evs.copyOfRealDesign.Columns(real_ev+3,evs.nReal) = evs.copyOfRealDesign.Columns(real_ev+2,evs.nReal-1);
      }
      /* do end points explicitly */
      evs.copyOfRealDesign(1,real_ev+2) = evs.copyOfRealDesign(2,real_ev+1) - evs.copyOfRealDesign(1,real_ev+1);
      evs.copyOfRealDesign(nTimepoints,real_ev+2) = evs.copyOfRealDesign(nTimepoints,real_ev+1) - evs.copyOfRealDesign(nTimepoints-1,real_ev+1);
      for(int t=1; t<nTimepoints-1 ; t++)
	evs.copyOfRealDesign(t+1,real_ev+2) = (evs.copyOfRealDesign(t+2,real_ev+1) - evs.copyOfRealDesign(t,real_ev+1))*0.5;
      orth_i_wrt_j(evs.copyOfRealDesign,real_ev+2,real_ev+1);	// well, Timmy wanted it but doesn't seem to have much effect.....
    }
    real_ev+=evs.nRealPerOrig[orig_ev];
  }

  evs.copyOfRealDesign=remmean(evs.copyOfRealDesign);

  if ( estimateRealHeights ) {
    evs.realHeights = estimate_X_heights(evs.copyOfRealDesign);
    // {{{ check rank of DM and get "real" contrast h2 heights (ie before HP filtering) 
    // first do real rank deficiency test
    // actually, no - we don't need this now we've switched to using pinv() below
    //eigenvals = feat_svd(evs.copyOfRealDesign);
    // now do "meaningful" rank deficiency test
    Matrix Q = pinv(evs.copyOfRealDesign.t() * evs.copyOfRealDesign);
    for(int con=1; con<=contrasts.nC(); con++) {
      ColumnVector contrast=contrasts.C.Column(con);
      Matrix X2 = evs.copyOfRealDesign * Q * contrast * pinv(contrast.t() * Q * contrast);
      contrasts.realHeights(con) = X2.Maximum() - X2.Minimum();
      if( fabs(contrasts.realHeights(con) ) > 1e6 ) 
	contrasts.realHeights(con) = 0; // try to catch dodgy heights, e.g. from empty EVs
    }
  } else { 
    writeSubmodel( filename, evs);
    // check rank of DM and do final contrast estimability test 
    // first do real rank deficiency test
    evs.eigenvals = feat_svd(evs.copyOfRealDesign);
    // now do "meaningful" rank deficiency test
    for(int con=1; con<=contrasts.nC(); con++) {
      ColumnVector dumbregressor( evs.copyOfRealDesign * contrasts.C.Column(con) ); /* just in order to test if it's empty */
      if ( dumbregressor.Maximum() - dumbregressor.Minimum() > 1e-10) { /* i.e. not an 'empty' contrast, e.g. from empty EVs */
	contrasts.RE(con) = RequiredEffect(evs.copyOfRealDesign, pinv(evs.copyOfRealDesign.t() * evs.copyOfRealDesign), contrasts.C.Column(con), preWhitening, critical_t, contrasts.realHeights(con), noise);
      } else {
	contrasts.realHeights(con)=0;
	contrasts.C.Column(con)=0;
	contrasts.RE(con) = 0;
      }
    }
  }
  return evs;
}

void writeDesignMatrix(const string filename, const ColumnVector& realDesignHeights, const Matrix& realDesign)
{
  FILE *outputFile;
  if ((outputFile=fopen(filename.c_str(),"wb"))==NULL) {
    cout << "Can't open " << filename << " for writing" << endl;
    exit(1);
  }
  
  fprintf(outputFile,"/NumWaves	%d\n", realDesign.Ncols() );
  fprintf(outputFile,"/NumPoints	%d\n", realDesign.Nrows() );

  fprintf(outputFile,"/PPheights	");
  for(int currentEV=1; currentEV <= realDesign.Ncols(); currentEV++)
    fprintf(outputFile,"	%e",realDesignHeights(currentEV));
  fprintf(outputFile,"\n");

  fprintf(outputFile,"\n/Matrix\n");

  for(int t=1; t <= realDesign.Nrows(); t++)
    {
      for(int currentEV=1; currentEV <= realDesign.Ncols(); currentEV++)
	fprintf(outputFile,"%e	",realDesign(t,currentEV));
      fprintf(outputFile,"\n");
    }

  fclose(outputFile);
}

void writeTriggers(const string filename, const vector<float>& triggers, const int* shape, const Regressors& evs, const double factor)
{
  FILE* triggerFile;
  if ((triggerFile=fopen(filename.c_str(),"wb"))==NULL) {
    cout << "Can't open " << filename << " for writing" << endl;
    exit(1);
  }

  for( int orig_ev=0; orig_ev<evs.nOriginal; orig_ev++ )
    for( int i=0; i<evs.nRealPerOrig[orig_ev]; i++ ) {
      if ( (shape[orig_ev]!=4) && (triggers[evs.nOriginal+orig_ev]>-1e6) ) { /* check for at least one trigger pair */
		
	ColumnVector diffs_list(2*evs.nTimepoints);
	int j;
	
	for(j=0;triggers[j*evs.nOriginal+orig_ev]>-1e6;j+=2) {
	  fprintf(triggerFile,"%e ",triggers[j*evs.nOriginal+orig_ev]);
	  if (triggers[(j+1)*evs.nOriginal+orig_ev]>-1e6)
	    diffs_list(j/2+1)=triggers[(j+1)*evs.nOriginal+orig_ev]-triggers[j*evs.nOriginal+orig_ev];
	}

	diffs_list=diffs_list.Rows(1,j/2);
	fprintf(triggerFile,"%f\n",median(diffs_list) + factor);
      }
      else
	fprintf(triggerFile,"0\n");
    }
  fclose(triggerFile);
}

void writeSubmodel( const string filename, const Regressors& evs)
{
  FILE* outputFile;
  if ( (outputFile=fopen(filename.c_str(),"wb") )==NULL) {
    cout << "Can't open " << filename << " for writing" << endl;
    exit(1);
  }
  for(int orig_ev=0; orig_ev<evs.nOriginal; orig_ev++)
    for(int tmp=0; tmp<evs.nRealPerOrig[orig_ev]; tmp++)
      fprintf(outputFile,"%d\n",orig_ev+1);
  fclose(outputFile);
}

void writeContrastMatrix(const string filename, const Contrasts& inputContrasts, const ColumnVector& RE )
{
  FILE * contrastFile;
  if ( (contrastFile=fopen(filename.c_str(),"wb")) ==NULL) {
    cout << "Can't open " << filename << " for writing" << endl;
    exit(1);
  }
  Matrix contrasts(inputContrasts.C.t());
  for(int con=1; con<=contrasts.Nrows(); con++)
    fprintf(contrastFile,"/ContrastName%d	%s\n",con,inputContrasts.name[con-1].c_str());
  fprintf(contrastFile,"/NumWaves	%d\n", contrasts.Ncols() );
  fprintf(contrastFile,"/NumContrasts	%d\n", contrasts.Nrows() );
  /* for contrasts, estimate "regressor height" as the mean of the heights of the relevant EVs */
  fprintf(contrastFile,"/PPheights	");
  for(int con=1; con<=contrasts.Nrows(); con++)
    fprintf(contrastFile,"	%e",inputContrasts.realHeights(con));
    
  fprintf(contrastFile,"\n");

  fprintf(contrastFile,"/RequiredEffect	");
  for(int con=1; con<=contrasts.Nrows(); con++)
    fprintf(contrastFile,"	%.3f",RE(con));
  fprintf(contrastFile,"\n");


  fprintf(contrastFile,"\n/Matrix\n");
  for(int con=1; con<=contrasts.Nrows(); con++) {
    for(int real_ev=1; real_ev<=contrasts.Ncols(); real_ev++) 
      fprintf(contrastFile,"%e ", contrasts(con,real_ev));
    fprintf(contrastFile,"\n");
  }

  fclose(contrastFile);
}

void writeFMatrix(const string filename, const Contrasts& contrasts, const Matrix& F, const int nFTests, const Matrix& realDesign )
{
  if (nFTests>0) {
    for(int f=1; f<=nFTests; f++) {
      Matrix Fmat(contrasts.nC(),contrasts.nEV());
      int Fmat_rows(0);

      for(int con=1; con<=contrasts.nC(); con++)
	if (F(f,con)) {
	  Fmat_rows++;
	  Fmat.Row(Fmat_rows)=( contrasts.C.Column(con)).t();
	}
    

      Fmat=Fmat.Rows(1,Fmat_rows);

      // test that F(X'X)^-1F' is invertible, i.e. of full rank
      if ( (Fmat.Nrows()==0) || ( MISCMATHS::rank(Fmat*pinv(realDesign.t()*realDesign)*Fmat.t()) < Fmat.Nrows() ) ) {
	cout << "F-test " << f << " isn't valid - each included contrast cannot be a linear combination of the others." << endl;
	exit(1);
      }
    }

    FILE *outputFile;
    if ((outputFile=fopen(filename.c_str(),"wb"))==NULL) {
      cout << "Can't open " << filename << " for writing" << endl;
      exit(1);
    }

    fprintf(outputFile,"/NumWaves	%d\n",contrasts.nC());
    fprintf(outputFile,"/NumContrasts	%d\n",nFTests);
    fprintf(outputFile,"\n/Matrix\n");

    for(int f=1; f<=nFTests; f++) {
      for(int con=1; con<=contrasts.nC(); con++)
	fprintf(outputFile,"%d ",(int)F(f,con));
      fprintf(outputFile,"\n");
    }
    fclose(outputFile);
  }	    
}

void writeGroupFile(const string filename, const vector<int> G, const int nTimepoints, const Matrix& realDesign)
{
int maxG(0), isok, isnotzero;
  for(int t=0;t<nTimepoints;t++) 
    maxG=MAX(G[t],maxG);

  FILE *outputFile;
  if ((outputFile=fopen(filename.c_str(),"wb"))==NULL) {
    cout << "Can't open " << filename << " for writing" << endl;
    exit(1);
  }
  /* if different group memberships check orthogonality of submatrices */
  if (maxG>1) {
    vector<float> sub_X(nTimepoints*maxG,0);
    vector<int> n_sub_X(maxG,0);
    for( int real_ev=1; real_ev<=realDesign.Ncols(); real_ev++) {
      for(int i=0; i < maxG; i++) 
	n_sub_X[i]=0;
      for(int t=0; t < nTimepoints; t++)
	sub_X[(G[t]-1)*nTimepoints+n_sub_X[G[t]-1]++] = realDesign(t+1,real_ev);

      isok=2;
      for(int i=0; i < maxG; i++) {
	isnotzero=0;
	for(int t=0; t < n_sub_X[i]; t++)
	  if ( sub_X[i*nTimepoints+t] != 0 )
	    isnotzero=1;
	isok-=isnotzero;
      }
      if ( isok < 1 )
	printf("Warning - design matrix uses different groups (for different variances), but these do not contain \"separable\" EVs for the different groups (it is necessary that, for each EV, only one of the groups has non-zero values). This message can be ignored if you are intending to use the groups file to define exchangeability blocks for randomise.\n");
    }
  }
  
  fprintf(outputFile,"/NumWaves	1\n");
  fprintf(outputFile,"/NumPoints	%d\n",nTimepoints);
  fprintf(outputFile,"\n/Matrix\n");

  for(int t=0; t<nTimepoints; t++)
    fprintf(outputFile,"%d\n",G[t]);

  fclose(outputFile);
}

void writeImagePreview( const string filename, const Contrasts& contrasts, const Matrix& F, const int nFTests, Matrix realDesign, const int level, const Regressors& evs, FONT_DATA *font_data, const vector<string>& titles, const float tr, const float nltffwhm, const int nTimepoints, const vector<int>& G )
{
  FILE *outputFile;
  unsigned char *r,*g,*b;
  int border=5, temp=15, fsize=11, nameLength=temp;
  char the_string[10000];

  /* how much space for contrast names? */
  for(int con =0; con <contrasts.nC(); con++) {
    int tmp =( strlen(contrasts.name[con].c_str() ) +3 ) * FONT_WIDTH;
    if ( tmp > nameLength ) nameLength=tmp;
  }

  int xmag = MIN( MAX( 1 , 600/contrasts.nEV()  ) , 50);
  int ymag = MIN( MAX( 1 , 400/nTimepoints ) , 20);

  int xsize = contrasts.nEV()*xmag + border*(contrasts.nEV()+3+nFTests+(nFTests>0)) + nameLength + nFTests*fsize;
  int ysize = nTimepoints*ymag + (contrasts.nC()+1)*FONT_HEIGHT + border*(4+contrasts.nC());

  // }}}
  // {{{ reset X[] range (but don't change offset) 

  for(int real_ev=1; real_ev <= contrasts.nEV(); real_ev++)
    if ( realDesign.Column(real_ev).MaximumAbsoluteValue() > 0 )
      realDesign.Column(real_ev) /= realDesign.Column(real_ev).MaximumAbsoluteValue();

  // }}}
  // {{{ malloc images and fill in background 

  r=(unsigned char *)malloc(xsize*ysize);
  g=(unsigned char *)malloc(xsize*ysize);
  b=(unsigned char *)malloc(xsize*ysize);

  memset((void *)r,(unsigned char)180,xsize*ysize);
  memset((void *)g,(unsigned char)215,xsize*ysize);
  memset((void *)b,(unsigned char)255,xsize*ysize);

  // }}}
  // {{{ time 

  if (level==1)
    for(int t=0; t<nTimepoints; t++)
      for(int y=0; y<ymag; y++)
	for(int x=0; x<temp; x++) {
	  int intensity(64), rintensity;
	  if ( t%10 == 0 ) 
	    intensity=255;
	  rintensity=intensity;
	  if ( abs(t-nTimepoints/2) <= (0.5*nltffwhm/tr) ) {
	    if (x>1 && x<temp-2)
	      rintensity=intensity=0;
	    if (x>3 && x<temp-4) {
	      rintensity=255;
	      intensity=0;
	    }
	  }
	  r[(border+t*ymag+y)*xsize+border+x] = rintensity;
	  g[(border+t*ymag+y)*xsize+border+x] = intensity;
	  b[(border+t*ymag+y)*xsize+border+x] = intensity;
	}

  // }}}
  // {{{ DM: grey 

  for(int t=0; t<nTimepoints; t++)
    for(int xx=0; xx<contrasts.nEV(); xx++)
      for(int y=0; y<ymag; y++)
	for(int x=0; x<xmag; x++)
	  {
	    float XX = realDesign(t+1,xx+1);

	    r[(border+t*ymag+y)*xsize+border*(xx+3)+xx*xmag+x+nameLength]=(int)((XX * 100.0)+128);
	    g[(border+t*ymag+y)*xsize+border*(xx+3)+xx*xmag+x+nameLength]=(int)((XX * 100.0)+128);
	    b[(border+t*ymag+y)*xsize+border*(xx+3)+xx*xmag+x+nameLength]=(int)((XX * 100.0)+128);
	  }

  // }}}
  // {{{ DM: red 

  for(int t=0; t<nTimepoints; t++)
    for(int xx=0; xx<contrasts.nEV(); xx++)
      for(int y=0; y<ymag; y++)
	{
	  float XX = realDesign(t+1,xx+1);
	  
	  int x = (int)((XX*0.8+1.0)*xmag/2);

	  r[(border+t*ymag+y)*xsize+border*(xx+3)+xx*xmag+x+nameLength]=255;
	  g[(border+t*ymag+y)*xsize+border*(xx+3)+xx*xmag+x+nameLength]=0;
	  b[(border+t*ymag+y)*xsize+border*(xx+3)+xx*xmag+x+nameLength]=0;
	  r[(border+t*ymag+y)*xsize+border*(xx+3)+xx*xmag+x+1+nameLength]=255;
	  g[(border+t*ymag+y)*xsize+border*(xx+3)+xx*xmag+x+1+nameLength]=0;
	  b[(border+t*ymag+y)*xsize+border*(xx+3)+xx*xmag+x+1+nameLength]=0;
	  r[(border+t*ymag+y)*xsize+border*(xx+3)+xx*xmag+x-1+nameLength]=255;
	  g[(border+t*ymag+y)*xsize+border*(xx+3)+xx*xmag+x-1+nameLength]=0;
	  b[(border+t*ymag+y)*xsize+border*(xx+3)+xx*xmag+x-1+nameLength]=0;
	  r[(border+t*ymag+y)*xsize+border*(xx+3)+xx*xmag+x+2+nameLength]=0;
	  g[(border+t*ymag+y)*xsize+border*(xx+3)+xx*xmag+x+2+nameLength]=0;
	  b[(border+t*ymag+y)*xsize+border*(xx+3)+xx*xmag+x+2+nameLength]=0;
	  r[(border+t*ymag+y)*xsize+border*(xx+3)+xx*xmag+x-2+nameLength]=0;
	  g[(border+t*ymag+y)*xsize+border*(xx+3)+xx*xmag+x-2+nameLength]=0;
	  b[(border+t*ymag+y)*xsize+border*(xx+3)+xx*xmag+x-2+nameLength]=0;
	}

  // }}}
  // {{{ EV names 
  int real_ev;
  for(int orig_ev=real_ev=0; orig_ev<evs.nOriginal; orig_ev++)
    {
      if ( orig_ev < evs.nOriginalNoMotion )
	 strcpy( the_string , titles[orig_ev].c_str() );
      else
	sprintf(the_string,"%s","conf");

      for(int i=0; i<evs.nRealPerOrig[orig_ev]; i++)
	{
	  write_string_rgb(r, g, b, 
			   border*(real_ev+3)+real_ev*xmag+xmag/2+nameLength-FONT_WIDTH*strlen(the_string)/2,
			   3*border+nTimepoints*ymag,
			   the_string, font_data, 150, 50, 50, xsize, ysize);
	  real_ev++;
	}

    }

  // }}}
  // {{{ contrasts 

  for(int yy=0; yy<contrasts.nC(); yy++)
    {
      sprintf(the_string,"C%d",yy+1);
      write_string_rgb(r, g, b, 
		       border,
		       (4+yy)*border+(yy+1)*FONT_HEIGHT+nTimepoints*ymag,
		       the_string, font_data, 0, 0, 0, xsize, ysize);

      strcpy(the_string,contrasts.name[yy].c_str());
      write_string_rgb(r, g, b, 
		       border+4*FONT_WIDTH,
		       (4+yy)*border+(yy+1)*FONT_HEIGHT+nTimepoints*ymag,
		       the_string, font_data, 150, 50, 50, xsize, ysize);

      for(int xx=0; xx<contrasts.nEV(); xx++)
	{
	  if ( contrasts.C(xx+1,yy+1) == (int)contrasts.C(xx+1,yy+1))
	    sprintf(the_string,"%d",(int)contrasts.C(xx+1,yy+1));
	  else
	    sprintf(the_string,"%.3f",contrasts.C(xx+1,yy+1));
      
	  write_string_rgb(r, g, b, 
			   border*(xx+3)+xx*xmag+xmag/2+nameLength-FONT_WIDTH*strlen(the_string)/2,
			   (4+yy)*border+(yy+1)*FONT_HEIGHT+nTimepoints*ymag,
			   the_string, font_data, 0, 0, 0, xsize, ysize);
	}
    }

  // }}}
  // {{{ ftests 

  for(int f=0; f<nFTests; f++)
    {
      write_string_rgb(r, g, b,
		       border*(contrasts.nEV()+4)+contrasts.nEV()*xmag+nameLength+f*(border+fsize),
		       3*border-FONT_HEIGHT+nTimepoints*ymag,
		       "F", font_data, 0, 0, 0, xsize, ysize);

      sprintf(the_string,"%d",f+1);
      write_string_rgb(r, g, b,
		       border*(contrasts.nEV()+4)+contrasts.nEV()*xmag+nameLength+f*(border+fsize),
		       3*border+nTimepoints*ymag,
		       the_string, font_data, 0, 0, 0, xsize, ysize);

      for(int con=0; con<contrasts.nC(); con++) {
	for(int y=0; y<fsize; y++)
	  for(int x=0; x<fsize; x++) {
	    long offset( ((4+con)*border+(con+1)*FONT_HEIGHT+nTimepoints*ymag+y)*xsize + border*(contrasts.nEV()+4)+contrasts.nEV()*xmag+nameLength+f*(border+fsize)+x );
	    r[offset] = g[offset] = b[offset] = 0; //initialise to black
	      if ( x > 0 && y > 0 && x < fsize-1 && y < fsize-1 ) {
		if (F(f+1,con+1)==1)
		  r[offset] = 175;
		else {
		  r[ offset ]=180;
		  g[ offset ]=215;
		  b[ offset ]=255;
		}
	      }
	  }
      }
    }
  // }}}
  // {{{ second-level group memberships 

  if (level==2)
    for(int t=0;t<nTimepoints;t++)
      {
	sprintf(the_string,"%d",G[t]);
	write_string_rgb(r, g, b,
			 border,
			 border+t*ymag,
			 the_string, font_data, 0, 0, 0, xsize, ysize);
      }

  // }}}
  // {{{ output image 


  if ((outputFile=fopen(filename.c_str(),"wb"))==NULL)
    {
      cout << "Can't open " << filename << " for writing" << endl;
      exit(1);
    }

  fprintf(outputFile,"P6\n");
  fprintf(outputFile,"%d %d\n",xsize,ysize);
  fprintf(outputFile,"255\n");

  for(int y=0; y < ysize; y++)
    for(int x=0; x < xsize; x++)
      {
	fwrite(&r[y*xsize+x],1,1,outputFile);
	fwrite(&g[y*xsize+x],1,1,outputFile);
	fwrite(&b[y*xsize+x],1,1,outputFile);
      }

  fclose(outputFile);
}

void writeCovarianceImage(string filename, const Contrasts& contrasts, const Matrix& F, const int nFTests, const Matrix& realDesign, const int level, const DiagonalMatrix eigenvals, FONT_DATA *font_data, const ColumnVector& RE)
{
FILE *outputfile;
unsigned char *r,*g,*b;
int border=5, mag=0, size=0, xsize=200, ysize;
float *cov=NULL, *covnorm=NULL;
char the_string[10000];

   if ( contrasts.nEV() > 1 ) {
	cov = (float *)malloc(contrasts.nEV()*contrasts.nEV()*sizeof(float));
	covnorm = (float *)malloc(contrasts.nEV()*sizeof(float));
	mag = MIN( MAX( 1 , 300/contrasts.nEV()  ) , 50);
	size = contrasts.nEV()*mag;
	xsize = size*2 + border*3;
    }

    ysize = size   + border*2 + (level==1)*(border*( contrasts.nC()+2 ) + FONT_HEIGHT*( contrasts.nC()+1 ));

    // }}}

    // {{{ malloc images and fill in background 

    r=(unsigned char *)malloc(xsize*ysize);
    g=(unsigned char *)malloc(xsize*ysize);
    b=(unsigned char *)malloc(xsize*ysize);

    memset((void *)r,(unsigned char)180,xsize*ysize);
    memset((void *)g,(unsigned char)215,xsize*ysize);
    memset((void *)b,(unsigned char)255,xsize*ysize);

    // }}}

    if (contrasts.nEV()>1) {						  
	// {{{ cov matrix orig 
	for(int evx=0; evx<contrasts.nEV(); evx++)
	  for(int evy=0; evy<contrasts.nEV(); evy++) {
	    cov[evy*contrasts.nEV()+evx]= (realDesign.Column(evx+1).t() * realDesign.Column(evy+1)).AsScalar();
	    if (evx==evy)
	      covnorm[evx]=cov[evx*contrasts.nEV()+evx];
	  }


	for(int evx=0; evx<contrasts.nEV(); evx++)
	  for(int evy=0; evy<contrasts.nEV(); evy++) {
	    float tmpf( covnorm[evy]*covnorm[evx] );
	    if (tmpf<=0)
	      cov[evy*contrasts.nEV()+evx]=0;
	    else
	      cov[evy*contrasts.nEV()+evx]=abs(cov[evy*contrasts.nEV()+evx])/sqrt(tmpf);
	  }

	for(int evx=0; evx<contrasts.nEV(); evx++)
	  for(int evy=0; evy<contrasts.nEV(); evy++)
	    for(int y=0; y<mag; y++)
	      for(int x=0; x<mag; x++)
		r[(border+evy*mag+y)*xsize+border+evx*mag+x] =
		  g[(border+evy*mag+y)*xsize+border+evx*mag+x] =
		  b[(border+evy*mag+y)*xsize+border+evx*mag+x] = (int)(255*cov[evy*contrasts.nEV()+evx]);

	// }}}
	// {{{ cov matrix SVD 

	for(int evx=0; evx<contrasts.nEV(); evx++) 
	  for(int evy=0; evy<contrasts.nEV(); evy++) {
	    int tmp=0;

	    if (evx==evy)
	      tmp=(int)(255*(eigenvals(evx+1)/eigenvals.Maximum()));

	    for(int y=0; y<mag; y++)
	      for(int x=0; x<mag; x++) {
		long offset = (border+evy*mag+y)*xsize+2*border+size+evx*mag+x;
		r[offset] = g[offset] = b[offset] = tmp;
	      }
	  }

	for(int evx=0; evx<contrasts.nEV(); evx++) 
	  for(int evy=0; evy<contrasts.nEV(); evy++)
	    if (evx==evy) {
	      sprintf(the_string,"%.3f",eigenvals(evx+1)/eigenvals.Maximum());
	      write_string_rgb(r, g, b, 
			       3*border+size+evx*mag + evx*(mag-2*border-5*FONT_WIDTH)/(contrasts.nEV()-1),
			       border+evy*mag+mag/2-FONT_HEIGHT/2,
			       the_string, font_data, 255, 0, 0, xsize, ysize);
	    }
    }

    if (level==1) { // {{{ contrast efficiencies 

	sprintf(the_string,"   Effect required (%%)");
	write_string_rgb(r, g, b, 
			 border,
			 size+3*border,
			 the_string, font_data, 0, 0, 0, xsize, ysize);

	for(int yy=0; yy<contrasts.nC(); yy++) {
	    sprintf(the_string,"C%d %.3f",yy+1,RE(yy+1));
      
	    int red_colour=0;
	    /*if (RE(yy+1)>5)
	      red_colour=255;*/

	    write_string_rgb(r, g, b, 
			     border,
			     size+(4+yy)*border+(yy+1)*FONT_HEIGHT,
			     the_string, font_data, red_colour, 0, 0, xsize, ysize);
	}
    }

    // }}}

    // {{{ output image 

    if ((outputfile=fopen(filename.c_str(),"wb"))==NULL) {
      cout << "Can't open " << filename << " for writing" << endl;
      exit(1);
    }

    fprintf(outputfile,"P6\n");
    fprintf(outputfile,"%d %d\n",xsize,ysize);
    fprintf(outputfile,"255\n");

    for(int y=0; y<ysize; y++)
      for(int x=0; x<xsize; x++) {
	fwrite(&r[y*xsize+x],1,1,outputfile);
	fwrite(&g[y*xsize+x],1,1,outputfile);
	fwrite(&b[y*xsize+x],1,1,outputfile);
      }

    fclose(outputfile);

    filename=string(getenv("FSLDIR")) + "/bin/wpng -q -overwrite  " + filename;
    system(filename.c_str());
}
