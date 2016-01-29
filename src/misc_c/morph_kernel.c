/* {{{ Copyright etc. */

/*  morph_kernel - create 3D 0/1 ASCII morphology kernels in a very proper way

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
#include <math.h>

/* }}} */

int main(argc,argv)
  int argc;
  char *argv[];
{
  /* {{{ variables */

float a,b,c,d,r,delta=0.05;
int   w,x,y,z;

/* }}} */

  /* {{{ usage and read arg */

if (argc<3)
{
  printf("Usage: morph_kernel <cube_side_length> <sphere_radius>\n");
  printf("e.g.: morph_kernel 11 5.5 > sphere111111.ker\n");
  exit(1);
}

w=atoi(argv[1]);
r=atof(argv[2]);

d=(w-1.0)/2;
r+=0.1; /* this so that x.5 values do actually include x.5 voxels */
r=r*r;   /* this to save having to do a sqrt */

/* }}} */
  /* {{{ print header */

printf("%%!XC-Volume\n\n/PixelEncoding  /Unsigned       def\n/Width  %d       def\n/Height %d       def\n/Slices %d       def\n/Depth  8       def\n\nBeginVolumeSamples\n\n",w,w,w);

/* }}} */
  /* {{{ output 0s and 1s */

for (z=0; z<w; z++)
{
  for (y=0; y<w; y++)
    {
      for (x=0; x<w; x++)
	{
	  int count=0, sum=0;

	  for (c=-d-0.5; c<-d+0.5+delta/2; c+=delta)
	    for (b=-d-0.5; b<-d+0.5+delta/2; b+=delta)
	      for (a=-d-0.5; a<-d+0.5+delta/2; a+=delta)
		{
		  count++;
		  if ( (z+c)*(z+c) + (y+b)*(y+b) + (x+a)*(x+a) <= r )
		    sum++;
		}

	  printf("%d ",(int)(sum>(int)(count/2)));
	}
      printf("\n");
    }

  printf("\n");
}

/* }}} */
  /* {{{ print footer */

  printf("EndVolumeSamples\n");

/* }}} */

  return 0;
}
