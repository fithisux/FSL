/* {{{ Copyright etc. */

/*  create_lut - create MEDx LUTs and material maps for colour rendering

    Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 1999-2001 University of Oxford  */

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

#define NLUTS 3

/* note that for number 2 (the third map) only the render3.lut makes sense */

int main(argc,argv)
  int argc;
  char *argv[];
{
  /* {{{ variables */

FILE   *lut[NLUTS], *map[NLUTS], *lutt[NLUTS], *mapt[NLUTS];
int    i1=0, i2=0, j, k, l;
double a, intensity=0, red, green, blue;
char   filename[1000];

/* }}} */

  /* {{{ usage */

  if (argc<2)
  {
    printf("Usage: create_lut <output_file_root>\n");
    exit(1);
  }

/* }}} */
  /* {{{ prepare output files and print starting texts */

for(l=0;l<NLUTS;l++)
{
  sprintf(filename,"%s%d.lut",argv[1],l+1);  lut[l]=fopen(filename,"wb");
  sprintf(filename,"%s%d.map",argv[1],l+1);  map[l]=fopen(filename,"wb");

  sprintf(filename,"%s%dt.lut",argv[1],l+1); lutt[l]=fopen(filename,"wb");
  sprintf(filename,"%s%dt.map",argv[1],l+1); mapt[l]=fopen(filename,"wb");

  fprintf(lut[l],"%%!VEST-LUT\n%%%%BeginInstance\n<<\n/SavedInstanceClassName /ClassLUT \n/PseudoColorMinimum 0.00 \n/PseudoColorMaximum 1.00 \n/PseudoColorMinControl /Low \n/PseudoColorMaxControl /High \n/PseudoColormap [\n");
  fprintf(map[l],"%%!VEST-MaterialMap\n%%%%BeginInstance\n<<\n/SavedInstanceClassName /ClassMaterialMap \n/Name (Render Map) \n/Colors [\n");
  fprintf(lutt[l],"%%!VEST-LUT\n%%%%BeginInstance\n<<\n/SavedInstanceClassName /ClassLUT \n/PseudoColorMinimum 0.00 \n/PseudoColorMaximum 1.00 \n/PseudoColorMinControl /Low \n/PseudoColorMaxControl /High \n/PseudoColormap [\n");
  fprintf(mapt[l],"%%!VEST-MaterialMap\n%%%%BeginInstance\n<<\n/SavedInstanceClassName /ClassMaterialMap \n/Name (Render Map) \n/Colors [\n");
}

/* }}} */
  /* {{{ colour maps for luts and mmaps */

  /* {{{ type 0 */

for (a=0; a<1; a+=0.01)
{
  fprintf(lut[0],"<-color{%f,%f,%f}->\n",a,a,a);
  fprintf(map[0],"<-color{%f,%f,%f}->\n",a,a,a);
  i1++;
}
for (a=0; a<1; a+=0.01)
{
  fprintf(lut[0],"<-color{1.0,%f,0.0}->\n",a);
  fprintf(map[0],"<-color{1.0,%f,0.0}->\n",a);
  i1++;
}

/* }}} */
  /* {{{ type 1 */

for (a=0; a<1; a+=0.01)
{
  fprintf(lut[1],"<-color{%f,%f,%f}->\n",a,a,a);
  fprintf(map[1],"<-color{%f,%f,%f}->\n",a,a,a);
  i2++;
}
for (a=0; a<1; a+=0.01)
{
  fprintf(lut[1],"<-color{1.0,%f,0.0}->\n",a);
  fprintf(map[1],"<-color{1.0,%f,0.0}->\n",a);
  i2++;
}
for (a=0; a<1; a+=0.01)
{
  fprintf(lut[1],"<-color{0.0,%f,1.0}->\n",a);
  fprintf(map[1],"<-color{0.0,%f,1.0}->\n",a);
  i2++;
}
for (a=0; a<1; a+=0.005)
{
  fprintf(lut[1],"<-color{0.0,%f,0.0}->\n",a/2+0.5);
  fprintf(map[1],"<-color{0.0,%f,0.0}->\n",a/2+0.5);
  i2++;
}

/* }}} */
  /* {{{ type 2 */

for (a=0; a<1; a+=0.01) /* LB-B */
{
  fprintf(lut[2],"<-color{0.0,%f,1.0}->\n",1-a);
  fprintf(map[2],"<-color{0.0,%f,1.0}->\n",1-a);
  i2++;
}
for (a=0; a<1; a+=0.02) /* B-G */
{
  fprintf(lut[2],"<-color{%f,%f,%f}->\n",a/2,a/2,1-a/2);
  fprintf(map[2],"<-color{%f,%f,%f}->\n",a/2,a/2,1-a/2);
  i2++;
}
for (a=0; a<1; a+=0.02) /* G-R */
{
  fprintf(lut[2],"<-color{%f,%f,%f}->\n",(1+a)/2,(1-a)/2,(1-a)/2);
  fprintf(map[2],"<-color{%f,%f,%f}->\n",(1+a)/2,(1-a)/2,(1-a)/2);
  i2++;
}
for (a=0; a<1; a+=0.01) /* R-Y */
{
  fprintf(lut[2],"<-color{1.0,%f,0.0}->\n",a);
  fprintf(map[2],"<-color{1.0,%f,0.0}->\n",a);
  i2++;
}

/* }}} */

#define STARTINT 0.25
#define K1 5
#define K2 10
#define K3 20

  /* {{{ type t0 */

for (a=0; a<1; a+=0.01)
{
  fprintf(lutt[0],"<-color{%f,%f,%f}->\n",a,a,a);
  fprintf(mapt[0],"<-color{%f,%f,%f}->\n",a,a,a);
}

for(j=0; j<5; j++)
  for(k=0; k<K3; k++)
  {
    if (k<K1)              intensity = STARTINT;
    if ((k>=K1) && (k<K2)) intensity = ((double)k-K1)/((K2-K1-1)/(1-STARTINT)) + STARTINT;
    if (k>=K2)             intensity = 1;

    red=1.0 * intensity;
    green=((double)j)/4.0 * intensity;

    fprintf(lutt[0],"<-color{%f,%f,0.0}->\n",red,green);
    fprintf(mapt[0],"<-color{%f,%f,0.0}->\n",red,green);
  }

/* }}} */
  /* {{{ type t1 */

for (a=0; a<1; a+=0.01)
{
  fprintf(lutt[1],"<-color{%f,%f,%f}->\n",a,a,a);
  fprintf(mapt[1],"<-color{%f,%f,%f}->\n",a,a,a);
}

for(j=0; j<5; j++)
  for(k=0; k<K3; k++)
  {
    if (k<K1)              intensity = STARTINT;
    if ((k>=K1) && (k<K2)) intensity = ((double)k-K1)/((K2-K1-1)/(1-STARTINT)) + STARTINT;
    if (k>=K2)             intensity = 1;

    red=1.0 * intensity;
    green=((double)j)/4.0 * intensity;

    fprintf(lutt[1],"<-color{%f,%f,0.0}->\n",red,green);
    fprintf(mapt[1],"<-color{%f,%f,0.0}->\n",red,green);
  }

for(j=0; j<5; j++)
  for(k=0; k<K3; k++)
  {
    if (k<K1)              intensity = STARTINT;
    if ((k>=K1) && (k<K2)) intensity = ((double)k-K1)/((K2-K1-1)/(1-STARTINT)) + STARTINT;
    if (k>=K2)             intensity = 1;

    blue=1.0 * intensity;
    green=((double)j)/4.0 * intensity;

    fprintf(lutt[1],"<-color{0.0,%f,%f}->\n",green,blue);
    fprintf(mapt[1],"<-color{0.0,%f,%f}->\n",green,blue);
  }

for(j=0; j<5; j++)
  for(k=0; k<K3*2; k++)
  {
    if (k<K1*2)                intensity = STARTINT;
    if ((k>=K1*2) && (k<K2*2)) intensity = ((double)k-K1*2)/((K2*2-K1*2-1)/(1-STARTINT)) + STARTINT;
    if (k>=K2*2)               intensity = 1;

    green=(((double)j)/8.0+0.5) * intensity;

    fprintf(lutt[1],"<-color{0.0,%f,0.0}->\n",green);
    fprintf(mapt[1],"<-color{0.0,%f,0.0}->\n",green);
  }

/* }}} */

/* }}} */
  /* {{{ rest of material map */

for(l=0;l<NLUTS;l++)
{
  fprintf(map[l],"]\n/Names [\n");
  fprintf(mapt[l],"]\n/Names [\n");
}

for (j=0; j<i1; j++)
{
  fprintf(map[0],"(Matrl%d)\n",j);
  fprintf(mapt[0],"(Matrl%d)\n",j);
}
for (j=0; j<i2; j++)
{
  fprintf(map[1],"(Matrl%d)\n",j);
  fprintf(mapt[1],"(Matrl%d)\n",j);
} 
 
for(l=0;l<NLUTS;l++)
{
  fprintf(map[l],"]\n/MinValues [\n");
  fprintf(mapt[l],"]\n/MinValues [\n");
}

for (j=0; j<i1; j++)
{
  fprintf(map[0],"%f\n",((double)j)*32000.0/((double)i1));
  fprintf(mapt[0],"%f\n",((double)j)*32000.0/((double)i1));
}
for (j=0; j<i2; j++)
{
  fprintf(map[1],"%f\n",((double)j)*32000.0/((double)i2));
  fprintf(mapt[1],"%f\n",((double)j)*32000.0/((double)i2));
} 
 
for(l=0;l<NLUTS;l++)
{
  fprintf(map[l],"]\n/MaxValues [\n");
  fprintf(mapt[l],"]\n/MaxValues [\n");
}

for (j=1; j<=i1; j++)
{
  fprintf(map[0],"%f\n",((double)j)*32000.0/((double)i1));
  fprintf(mapt[0],"%f\n",((double)j)*32000.0/((double)i1));
}
for (j=1; j<=i2; j++)
{
  fprintf(map[1],"%f\n",((double)j)*32000.0/((double)i2));
  fprintf(mapt[1],"%f\n",((double)j)*32000.0/((double)i2));
}

for(l=0;l<NLUTS;l++)
{
  fprintf(map[l],"]\n/Opacities [\n");
  fprintf(mapt[l],"]\n/Opacities [\n");
}

for (j=0; j<i1/2; j++)
{
  fprintf(map[0],"%f\n",((double)j)*0.05/(0.5*(double)i1)+0.05);
  fprintf(mapt[0],"%f\n",((double)j)*0.05/(0.5*(double)i1)+0.05);
}
for (j=i1/2; j<i1; j++)
{
  fprintf(map[0],"0.1\n");
  fprintf(mapt[0],"0.1\n");
}

for (j=0; j<i2/5; j++)
{
  fprintf(map[1],"%f\n",((double)j)*0.05/(0.2*(double)i2)+0.05);
  fprintf(mapt[1],"%f\n",((double)j)*0.05/(0.2*(double)i2)+0.05);
}
for (j=i2/5; j<i2; j++)
{
  fprintf(map[1],"0.1\n");
  fprintf(mapt[1],"0.1\n");
}

for(l=0;l<NLUTS;l++)
{
  fprintf(map[l],"]\n/RenderMaterial? [\n");
  fprintf(mapt[l],"]\n/RenderMaterial? [\n");
}

for (j=1; j<=i1; j++)
{
  fprintf(map[0],"true\n");
  fprintf(mapt[0],"true\n");
}
for (j=1; j<=i2; j++)
{
  fprintf(map[1],"true\n");
  fprintf(mapt[1],"true\n");
}

/* }}} */
  /* {{{ end texts for luts and mmaps */

for(l=0;l<NLUTS;l++)
{
  fprintf(lut[l],"]\n>>\n\n%%%%EndInstance\n%%%%EOF\n");
  fprintf(map[l],"]\n>>\n\n%%%%EndInstance\n%%%%EOF\n");
  fprintf(lutt[l],"]\n>>\n\n%%%%EndInstance\n%%%%EOF\n");
  fprintf(mapt[l],"]\n>>\n\n%%%%EndInstance\n%%%%EOF\n");
}

/* }}} */

  return(0);
}

/* }}} */
