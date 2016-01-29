/* {{{ Copyright etc. */

/*  hist2prob - convert histogram into probability values

    Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

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
/* {{{ includes and defines */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

/* }}} */
/* {{{ main() */

int main(argc,argv)
  int argc;
  char *argv[];
{
/* {{{ variables */

FILE   *ifp;
int    i, image_size, count, n_points=1000, total;
float  *hist_data,low_thresh,high_thresh,incr,tmpf;

/* }}} */

  /* {{{ setup files, etc */

  /*printf("hist2prob by Steve Smith\n");*/

  if (argc<5)
  {
    printf("Usage: hist2prob <image> <size> <low_thresh> <high_thresh>\n");
    exit(1);
  }

  if ((ifp=fopen(argv[1],"rb")) == NULL)
  {
    printf("File %s isn't readable\n",argv[1]);
    exit(-1);
  }
  image_size=atoi(argv[2]);
  low_thresh=atof(argv[3]);
  high_thresh=atof(argv[4]);
  incr=(high_thresh-low_thresh)/(float)n_points;
  hist_data=(float*)malloc(image_size*sizeof(float));
  fread(hist_data,sizeof(float),image_size,ifp);

  count=0;
  for(i=0; i<image_size; i++)
    if ( (hist_data[i]>low_thresh) && (hist_data[i]<high_thresh) )
      count++;

  /*printf("count=%d\n",count);*/

  for(tmpf=high_thresh;tmpf>low_thresh;tmpf-=incr)
  {
    total=0;
    for(i=0; i<image_size; i++)
      if ( (hist_data[i]>tmpf) && (hist_data[i]<high_thresh) )
        total++;
    total /= 2; /* this is to make it a one-tailed test? */
    printf("%f %f\n",tmpf,log10(((float)total)/((float)count)));
  }

/* }}} */
  return (0);
}

/* }}} */
