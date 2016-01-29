/* {{{ Copyright etc. */

/*  prewhiten - apply prehitening

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
/* {{{ background theory */

/*

consider the GLM:

Y = Xb + E

consider a single voxel, in which case Y (data) and X (model) and E
(residuals) are column vectors (ie time series)

we generate a prewhitening matrix S which makes the residuals white:

SY = SXb + SE

(SPM also uses this form, but in the case of SPM, S is not a whitening
matrix but a colouring matrix, so that instead of the resulting
residuals being white, they contain an amount of autcorrelation which
is approximately known and can be corrected for in the inference.)

S is a Toeplitz matrix

S =   / r0  r1  r2  r3  r4  ....... \
     |  r-1 r0  r1  r2  r3  ....... |
     |  r-2 r-1 r0  r1  r2  ....... |
     |  r-3 r-2 r-1 r0  r1  ....... |
        ...........................
        ...........................
        ...........................
     |  ....... r-2 r-1 r0  r1  r2  |
     |  ....... r-3 r-2 r-1 r0  r1  |
     \  ....... r-4 r-3 r-2 r-1 r0  /

this is in theory equivalent to convolving the time series in a column
vector (Y or X or E) with 1D convolution kernel r[i];

SY[i] = SUM(r[j-i] * Y[j])

the difference in practice is that it is easier to sort out
end-effects (padding of the data to allow simple processing with a
constant width convolution kernel) if this is seen as a convolution
rather than matrix multiplication. In fact, the convolution is
normally carried out in fourier space (ie a fourier space
multiplication rather than real space convolution) for computational
efficiency.


***********************************************************

generating prewhitening matrix S:

E = N(0,sigma^2 V)

where V is the autocorrelation matrix

use Cholesky to get K where V=KK'

so if we set S = K^-1

there are more details in Mark's report which we've just put on the web:
http://www.fmrib.ox.ac.uk/analysis/techrep/tr01mw1/tr01mw1/node3.html

in detail:

1) fit model to data and remove that part of signal which correlates.
2) take resulting residuals and estimate autocorrelation parameters
3) regularise these by applying a Tukey taper and the smoothing spatially with the estimates at neighbouring voxels
4) create S=K^-1 by inverting the autocorrelation estimates in the spectral domain
5) Use S to prewhiten data and model, and refit.

again - there is more detail in the paper (and the source code ;-)


*/

/* }}} */
/* {{{ defines, includes and typedefs */

#include "featlib.h"
using namespace NEWMAT;
using namespace NEWIMAGE;

/* }}} */
/* {{{ usage */

void usage(void)
{
  printf("Usage: prewhiten <feat_directory.feat> [options]\n");
  printf("[-o <output_directory>] change output directory from default of input feat directory\n");
  exit(1);
}

/* }}} */

int main(int argc, char **argv)
{
  /* {{{ variables */

int          argi=1, x, y, z;
char         fmridata[10000], featdir[10000], outputdir[10000], thestring[10000];
ColumnVector pwts;

/* }}} */

  /* {{{ process arguments */

if (argc<2) usage();

strcpy(featdir,argv[argi++]);
sprintf(fmridata,"%s/filtered_func_data",featdir);

strcpy(outputdir,featdir);

for (;argi<argc;argi++)
{
  if (!strcmp(argv[argi], "-o"))
    /* {{{ output dir */

{
  argi++;
  if (argc<argi+1)
    {
      printf("Error: no value given following -o\n");
      usage();
    }
  strcpy(outputdir,argv[argi]);
}

/* }}} */
}

/* }}} */
  /* {{{ read filtered_func_data */

volume4D<int> im;
read_volume4D(im, fmridata);

/* }}} */
  /* {{{ read auto correlation estimates for prewhitening */

volume4D<float> acs;
sprintf(thestring,"%s/stats/threshac1",featdir);

if(fsl_imageexists(thestring)) 
  read_volume4D(acs, thestring);

/* }}} */

  for(z=0;z<im.zsize();z++) for(y=0;y<im.ysize();y++) for(x=0;x<im.xsize();x++)
    {
      prewhiten_timeseries(acs.voxelts(x,y,z), im.voxelts(x,y,z), pwts, im.tsize());
      im.voxelts(x,y,z)=pwts;
    }

  sprintf(thestring,"%s/pw_res",outputdir);
  save_volume4D(im,thestring);

  return 0;
}

