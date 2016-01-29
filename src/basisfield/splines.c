/*     CCOPYRIGHT     */
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
#ifdef MEX
#include "mex.h"
#endif

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include "splines.h"
 
/*
This is used to ensure that memory allocated C-style is
also freed C-style. This is intended to use when this
module is called from C++.
*/
void please_free(void *ptr)
{
  if (ptr) {my_free(ptr);}
}
 
/* 
Calculates a 2- or 3D cubic B-spline given an array
of 1D splines (obtained by get_1D_spline).
*/
int spline_kron(/* Input */
                int           ndim,     /* Dimensionality of spline. */
                int           dim[3],   /* Size of spline in the different directions. */
                double        *sp1d[3], /* Set of 1D splines. */
                /* Output */
                double        *spline)  /* nD (n>0 & n<4) spline. */
{
   int     i=0, j=0, k=0;
   int     ldim[3];
   int     n=1;
   double  tmp = 1.0;
   double  *lsp1d[3];

   for (i=0; i<ndim; i++) {ldim[i] = dim[i]; n *= dim[i]; lsp1d[i] = sp1d[i];}
   for (i=ndim; i<3; i++) {ldim[i] = 1; lsp1d[i] = &tmp;}
   
   for (k=0; k<ldim[2]; k++)
   {
      for (j=0; j<ldim[1]; j++)
      {
         for (i=0; i<ldim[0]; i++)
	 {
	    spline[index(i,j,k,ldim)] = lsp1d[0][i]*lsp1d[1][j]*lsp1d[2][k];
         }
      }
   }

   return(n);
}

/* 
Calculates all non-zero values of a 1D cubic B-spline
given a specific knot spacing.
*/
int get_1D_spline(/* Input */
                  int     knsp,      /* Knot-spacing (in # of voxels) */
                  /* Output */
                  double  **spline)  /* 1D spline function */
{
   unsigned int    i=0;
   double          tmp=0.0;

   *spline = (double *) my_calloc(ksp2sz(knsp),sizeof(double));
   for (i=0; i<(knsp-1); i++)
   {
      tmp = ((double) i+1) / ((double) knsp);
      (*spline)[i] = (1.0/6.0) * tmp*tmp*tmp;
   }
   for (i=(knsp-1); i<(2*knsp); i++)
   {
      tmp = ((double) i+1) / ((double) knsp);
      (*spline)[i] = (2.0/3.0) - 0.5*(2.0-tmp)*(2.0-tmp)*tmp;
   }
   for (i=(2*knsp); i<(4*knsp-1); i++)
   {
     (*spline)[i] = (*spline)[(4*knsp-2)-i];
   }

   return(ksp2sz(knsp));
}

/* 
Calculates all non-zero values of the first spatial derivative of a
1D cubic B-spline given a specific knot spacing.
*/
int get_1D_spline_d(/* Input */
                    int     knsp,      /* Knot-spacing (in # of voxels) */
                    /* Output */
                    double  **spline)  /* 1D spline function */
{
   unsigned int    i=0;
   double          tmp=0.0;

   *spline = (double *) my_calloc(ksp2sz(knsp),sizeof(double));
   for (i=0; i<(knsp-1); i++)
   {
      tmp = ((double) i+1) / ((double) knsp);
      (*spline)[i] = (1.0 / ((double) knsp)) * 0.5 * tmp*tmp;
   }
   for (i=(knsp-1); i<(2*knsp); i++)
   {
      tmp = ((double) i+1) / ((double) knsp);
      (*spline)[i] = (1.0 / ((double) knsp)) *(-2.0*(tmp-2.0) - 1.5*(tmp-2.0)*(tmp-2.0));
   }
   for (i=(2*knsp); i<(4*knsp-1); i++)
   {
     (*spline)[i] = -(*spline)[(4*knsp-2)-i];
   }

   return(ksp2sz(knsp));
}

/* 
Calculates all non-zero values of the second spatial derivative of a
1D cubic B-spline given a specific knot spacing.
*/
int get_1D_spline_dd(/* Input */
                     int     knsp,      /* Knot-spacing (in # of voxels) */
                     /* Output */
                     double  **spline)  /* 1D spline function */
{
   unsigned int    i=0;
   double          tmp=0.0;

   *spline = (double *) my_calloc(ksp2sz(knsp),sizeof(double));
   for (i=0; i<(knsp-1); i++)
   {
      tmp = ((double) i+1) / ((double) knsp);
      (*spline)[i] = pow(1.0/((double) knsp),2.0) * tmp;
   }
   for (i=(knsp-1); i<(2*knsp); i++)
   {
      tmp = ((double) i+1) / ((double) knsp);
      (*spline)[i] = pow(1.0/((double) knsp),2.0) * (4.0 - 3.0*tmp);
   }
   for (i=(2*knsp); i<(3*knsp-1); i++)
   {
      tmp = ((double) i+1) / ((double) knsp);
     (*spline)[i] = pow(1.0/((double) knsp),2.0) * (3.0*tmp - 8.0);
   }
   for (i=(3*knsp-1); i<(4*knsp-1); i++)
   {
      tmp = ((double) i+1) / ((double) knsp);
     (*spline)[i] = pow(1.0/((double) knsp),2.0) * (4.0 - tmp);
   }

   return(ksp2sz(knsp));
}

/*
Calculates the new spline-coefficients when a field is zoomed
by a factor that is a power of 2 in one or more directions.
*/
int zoom_field(/* input */
               int      ndim,
               int      oksp[3],
               int      nksp[3],
               int      idim[3],
               double   *oc,
               /* Output */
               double   *nc)
{
   int      i=0;
   int      dim=0;
   int      sf=0;
   int      sz=0;
   int      tksp[3];
   double   *tmpc1=NULL;
   double   *tmpc2=NULL;

   for (i=ndim; i<3; i++) {oksp[i]=1; nksp[i]=1; idim[i]=1;}
   for (i=0, sz=1; i<3; i++) {tksp[i] = oksp[i]; sz *= no_of_knots(oksp[i],idim[i]);}

   tmpc1 = oc;
   
   for (dim=0; dim<3; dim++)
   {
      sf = oksp[dim]/nksp[dim];
      while (sf > 1)
      {
	 zoom_field_by2(3,tksp,idim,dim,tmpc1,&tmpc2);
         tksp[dim] /= 2;
         sf /= 2;
         if (tmpc1 != oc)
	 {
            my_free(tmpc1);
         }
         tmpc1 = tmpc2;
      }
      if (sf != 1)
      {
	 printf("\nzoom_field: zooming must be by power of 2.");
	 return(-1);
      }
   }

   for (i=0, sz=1; i<3; i++) {sz *= no_of_knots(nksp[i],idim[i]);}

   memcpy(nc,tmpc2,sz*sizeof(double));
   my_free(tmpc2);
 
   return(1);
}
/*
Calculates the new set of spline-coefficients resulting when the 
knot-spacing of an existing field is cut by half in one direction.
*/
int zoom_field_by2(/* Input */
                   int     ndim,
                   int     ksp[3],
                   int     idim[3],
                   int     zdim,
                   double  *oc,
                   /* Output */
                   double  **nc)
{
   int      i=0, j=0, k=0;
   int      nsz=0;
   int      ocdim[3], ncdim[3];
   double   a0=1.0/8.0, a4=a0;
   double   a1=1.0/2.0, a3=a1;
   double   a2=3.0/4.0;

   for (i=ndim; i<3; i++) {ocdim[i] = 1; idim[i] = 1; ksp[i] = 1;}
   for (i=0, nsz=1; i<3; i++) 
   {
      ocdim[i] = no_of_knots(ksp[i],idim[i]);
      ncdim[i] = (i==zdim) ? no_of_knots(ksp[i]/2,idim[i]) : ocdim[i]; 
      nsz *= ncdim[i];
   }

   *nc = (double *) my_calloc(nsz,sizeof(double));

   if (zdim==0)
   {
      for (k=0; k<ncdim[2]; k++)
      {
	 for (j=0; j<ncdim[1]; j++)
	 {
	    for (i=0; i<ncdim[0]; i++) 
	    {
	       if (i % 2)
	       {
		  (*nc)[index(i,j,k,ncdim)] = a4*oc[index((i+1)/2-1,j,k,ocdim)] + 
		                              a2*oc[index((i+1)/2,j,k,ocdim)];
                  if (((i+1)/2+1) < ocdim[0]) {
                    (*nc)[index(i,j,k,ncdim)] += a0*oc[index((i+1)/2+1,j,k,ocdim)];
		  }
               }
               else
	       {
		  (*nc)[index(i,j,k,ncdim)] = a3*oc[index(i/2,j,k,ocdim)];
                  if ((i/2+1) < ocdim[0]) { 
		    (*nc)[index(i,j,k,ncdim)] += a1*oc[index(i/2+1,j,k,ocdim)];
		  } 
               }
	    }
	 }
      }
   }
   else if (zdim==1)
   {
      for (k=0; k<ncdim[2]; k++)
      {
	 for (j=0; j<ncdim[1]; j++)
	 {
	    for (i=0; i<ncdim[0]; i++) 
	    {
	       if (j % 2)
	       {
		  (*nc)[index(i,j,k,ncdim)] = a4*oc[index(i,(j+1)/2-1,k,ocdim)] + 
		                              a2*oc[index(i,(j+1)/2,k,ocdim)]; 
                  if (((j+1)/2+1) < ocdim[1]) {
                    (*nc)[index(i,j,k,ncdim)] += a0*oc[index(i,(j+1)/2+1,k,ocdim)];
		  }
               }
               else
	       {
		  (*nc)[index(i,j,k,ncdim)] = a3*oc[index(i,j/2,k,ocdim)];
                  if ((j/2+1) < ocdim[1]) { 
		    (*nc)[index(i,j,k,ncdim)] += a1*oc[index(i,j/2+1,k,ocdim)]; 
		  } 
               }
	    }
	 }
      }
   }
   else if (zdim==2)
   {
      for (k=0; k<ncdim[2]; k++)
      {
	 for (j=0; j<ncdim[1]; j++)
	 {
	    for (i=0; i<ncdim[0]; i++) 
	    {
	       if (k % 2)
	       {
		  (*nc)[index(i,j,k,ncdim)] = a4*oc[index(i,j,(k+1)/2-1,ocdim)] + 
		                              a2*oc[index(i,j,(k+1)/2,ocdim)];
                  if (((k+1)/2+1) < ocdim[2]) {
                    (*nc)[index(i,j,k,ncdim)] += a0*oc[index(i,j,(k+1)/2+1,ocdim)];
		  }
               }
               else
	       {
		  (*nc)[index(i,j,k,ncdim)] = a3*oc[index(i,j,k/2,ocdim)];
                  if ((k/2+1) < ocdim[2]) { 
		    (*nc)[index(i,j,k,ncdim)] += a1*oc[index(i,j,k/2+1,ocdim)]; 
		  } 
               }
	    }
	 }
      }
   }
   else
   {
     return(-1);
   }

   return(1);
}

int get_field(/* Input */
              int     ndim,
              int     cdim[3],
              double  *c,
              int     sdim[3],
              double  *sp,
              int     fdim[3],
              /* Output */
              double  *f)
{
   int     i=0;
   int     ck=0, cj=0, ci=0;
   int     fk=0, fj=0, fi=0;
   int     sk=0, sj=0, si=0;
   int     fks=0, fjs=0, fis=0;
   int     fke=0, fje=0, fie=0;
   int     sks=0, sjs=0, sis=0;
   double  cv = 0.0;

   for (i=ndim; i<3; i++) {cdim[i]=1; sdim[i]=1; fdim[i]=1;}
   memset(f,0,fdim[0]*fdim[1]*fdim[2]*sizeof(double));

   for (ck=0; ck<cdim[2]; ck++)
   {
      get_range(ck,sdim[2],fdim[2],&fks,&fke,&sks);
      for (cj=0; cj<cdim[1]; cj++)
      {
	 get_range(cj,sdim[1],fdim[1],&fjs,&fje,&sjs);
	 for (ci=0; ci<cdim[0]; ci++)
	 {
            get_range(ci,sdim[0],fdim[0],&fis,&fie,&sis);
	    if ((cv = c[index(ci,cj,ck,cdim)]))
	    {
               for (fk=fks,sk=sks; fk<fke; fk++,sk++)
	       {
	          for (fj=fjs, sj=sjs; fj<fje; fj++, sj++)
	          {
		     for (fi=fis, si=sis; fi<fie; fi++, si++)
		     {
		        f[index(fi,fj,fk,fdim)] += cv*sp[index(si,sj,sk,sdim)];
                     }
                  }
               }
            }
         }
      }
   }
   return(1);
}


int get_range(/* Input */
              int      k,
              int      ssz,
              int      fsz,
              /* Output */
              int      *fs,
              int      *fe,
              int      *ks)
{
   if (ssz == 1) {(*fs)=0; (*fe)=fsz; (*ks)=0; return(1);}
   (*fs) = 1 + (k-3)*sz2ksp(ssz);
   (*fe) = (k+1)*sz2ksp(ssz) - 1;
   (*ks) = 0;
   if ((*fs) < 0)
   {
      (*ks) -= (*fs);
      (*fs) = 0;
   }
   if ((*fe) > fsz) {(*fe) = fsz;}

   return(1);
} 
  
int no_of_knots(int    ksp,
                int    msz)
{
   if (msz == 1) return(1); /* Collapsed dimension. */ 
   else return(((int) ceil((((double) msz) + 1.0) / ((double) ksp))) + 2);
}

int make_A(/* Input. */
           int      ndim,
           int      kdim[3],
           int      sdim[3],
           double   *spl,
           int      idim[3],
           double   *ima,
           /* Output. */
           int      *irp,
           int      *jcp,
           double   *sp)
{
   int     cntr = 0;
   int     i=0;
   int     ci=0;
   int     kk=0, kj=0, ki=0;
   int     ik=0, ij=0, ii=0;
   int     sk=0, sj=0, si=0;
   int     iks=0, ijs=0, iis=0;
   int     ike=0, ije=0, iie=0;
   int     sks=0, sjs=0, sis=0;

   for (i=ndim; i<3; i++) {kdim[i]=1; sdim[i]=1; idim[i]=1;}

   for (kk=0,ci=0; kk<kdim[2]; kk++)
   {
      get_range(kk,sdim[2],idim[2],&iks,&ike,&sks);
      for (kj=0; kj<kdim[1]; kj++)
      {
	 get_range(kj,sdim[1],idim[1],&ijs,&ije,&sjs);
	 for (ki=0; ki<kdim[0]; ki++)
	 {
            get_range(ki,sdim[0],idim[0],&iis,&iie,&sis);
            jcp[ci] = cntr; 
            for (ik=iks,sk=sks; ik<ike; ik++,sk++)
	    {
	       for (ij=ijs, sj=sjs; ij<ije; ij++,sj++)
	       {
		  for (ii=iis, si=sis; ii<iie; ii++,si++)
		  {
                     *sp++ = ima[index(ii,ij,ik,idim)] * spl[index(si,sj,sk,sdim)];
                     cntr++;
		     *irp++ = index(ii,ij,ik,idim);
                  }
               }
            }
            ci++;
         }
      }
   }
   jcp[ci] = cntr;
   
   return(cntr);
}

int make_Aty(/* Input. */
             int      ndim,      /* Actual dimensionality of problem (1, 2 or 3). */
             int      cdim[3],   /* # of knots in the three dimensions. */
             int      sdim[3],   /* Size of spline kernel in the three dimensions. */
             double   *spl,      /* Spline kernel. */
             int      idim[3],   /* Size of image matrix. */
             double   *ima,      /* Image. */
             double   *y,        /* y-vector. */
             /* Output. */
             double   *aty)      /* Resulting vector */
{
   int     i=0, j=0, k=0;
   int     ri = 0;              /* Row index for Aty. */
   double  *s_by_i = NULL;      /* Spline multiplied by appurtenant image intensities. */
   double  *y_by_i = NULL;      /* Elementwise product of image and y. */

   for (i=ndim; i<3; i++) {cdim[i]=1; sdim[i]=1; idim[i]=1;}

   s_by_i = (double *) my_calloc(sdim[0]*sdim[1]*sdim[2],sizeof(double));
   y_by_i = (double *) my_calloc(idim[0]*idim[1]*idim[2],sizeof(double));

   /* First multiply image and y-vector. */
   for (i=0; i<idim[0]*idim[1]*idim[2]; i++) {y_by_i[i] = ima[i]*y[i];}
   
   /* First three loops over splines determine the row of Aty. */
   ri = 0;
   for (k=0; k<cdim[2]; k++)
   {
      for (j=0; j<cdim[1]; j++) 
      {
	 for (i=0; i<cdim[0]; i++)
	 {
            /* Multiply this spline with product of image and y. */
	    aty[ri++] = get_s_by_i(i,j,k,sdim,idim,spl,y_by_i,s_by_i);
	 }
      }
   }

   my_free(s_by_i);
   my_free(y_by_i);

   return(ri);      
}

/*
Calculates the gradient of the membrane energy for the field given by the
spline-coefficients beta, without taking the detour via the H-matrix.
This way the user won't need to store the (potentially) quite big matrix H. 
The membrane energy can be expressed as b'*A_x'*A_x*b + b'*A_y'*A_y*b + b'*A_z'*A_z*b 
where A_i is a matrix with fsize rows and nbas columns, and where the jth column 
contains an unravelled version of a spline-kernel once differentiated in the ith 
direction.
If we denote the matrix (A_x'*A_x + A_y'*A_y + A_z'*A_z) by H, then the gradient
of the membrane energy is given by 2*H*b.
The helper function memen_AtAb simply calculates, depending on the what-flag, either A*b
or A'*A*b. 
*/
int get_memen_grad(/* Input. */
                   int           ndim,    /* # of dimensions (1,2 or 3) */
                   const int     cdim[3], /* Size of coefficient array. */
                   const int     ksp[3],  /* Knot-spacings. */
                   const double  *beta,   /* Coefficients. */
                   /* Output.*/
                   double        *grad)
{
   int     i=0, j=0;
   int     sdim[3] = {0,0,0}; /* Size of spline kernel. */
   int     ssz = 1;           /* Total size of spline kernel. */
   int     csz = 1;           /* Total # of coefficients. */
   int     AtAb_sz = 0;       /* # of rows in A'*A*b. */
   double  *spline = NULL;    /* Spline kernel. */
   double  *AtAb = NULL;   
   double  *sp_1D[3] = {NULL, NULL, NULL};
   double  *sp_1D_d[3] = {NULL, NULL, NULL};
   double  *tmp[3] = {NULL, NULL, NULL};


   /* First get 1D spline kernels and their derivatives. */
   for (i=0; i<ndim; i++) {
      csz *= cdim[i];
      ssz *= (sdim[i] = get_1D_spline(ksp[i],&(sp_1D[i])));
      get_1D_spline_d(ksp[i],&(sp_1D_d[i]));
   }

   memset(grad,0,csz*sizeof(double));
   spline = (double *) my_calloc(ssz,sizeof(double));

   /* Then do one dimension at a time. */
   for (i=0; i<ndim; i++) {
      memset(spline,0,ssz*sizeof(double));
      for (j=0; j<ndim; j++) {
	 if (i==j) {tmp[j] = sp_1D_d[j];}
         else {tmp[j] = sp_1D[j];}
      }
      spline_kron(ndim,sdim,tmp,spline);
      memen_AtAb(ndim,cdim,ksp,sdim,spline,beta,1,&AtAb_sz,&AtAb);
      for (j=0; j<csz; j++) {grad[j] += 2.0 * AtAb[j];}
      my_free(AtAb);
   }

   my_free(spline);
   for (i=0; i<ndim; i++) {my_free(sp_1D[i]); my_free(sp_1D_d[i]);}
                  
   return(csz);
}

/*
Calculates the membrane energy for the field given by the spline coefficients b, without
taking the detour via the H-matrix (see make_memen_H()). This way the user won't need
to store the (potentially) quite big matrix H. The membrane energy can be expressed as
b'*A_x'*A_x*b + b'*A_y'*A_y*b + b'*A_z'*A_z*b where A_i is a matrix with fsize rows and
nbas columns, and where the jth column contains an unravelled version of a spline-kernel
once differentiated in the ith direction.
Note that fsize is larger than nvox. This is because we base the membrane energy on the
full extent of the field, not just the part that fits into the image FOV.
The helper function memen_AtAb simply calculates, depending on the what-flag, either A*b
or A'*A*b. 
*/
double get_memen(/* Input. */
                 int           ndim,    /* # of dimensions (1,2 or 3) */
                 const int     cdim[3], /* Size of coefficient array. */
                 const int     ksp[3],  /* Knot-spacings. */
                 const double  *beta)   /* Coefficients. */
{
   int     i=0, j=0;
   int     sdim[] = {1,1,1}; /* Size of spline kernel. */
   int     ssz = 1;          /* Total size of spline kernel. */
   int     Ab_sz = 0;        /* # of rows in A*b. */
   double  memen = 0.0;      /* Membrane energy. */
   double  *spline = NULL;   /* Spline kernel. */
   double  *Ab = NULL;   
   double  *sp_1D[3] = {NULL, NULL, NULL};
   double  *sp_1D_d[3] = {NULL, NULL, NULL};
   double  *tmp[3] = {NULL, NULL, NULL};

   /* First get 1D spline kernels and their derivatives. */
   for (i=0; i<ndim; i++) {
      ssz *= (sdim[i] = get_1D_spline(ksp[i],&(sp_1D[i])));
      get_1D_spline_d(ksp[i],&(sp_1D_d[i]));
   }

   spline = (double *) my_calloc(ssz,sizeof(double));
   
   /* Then do one dimension at a time. */
   memen = 0.0;
   for (i=0; i<ndim; i++) {
      memset(spline,0,ssz*sizeof(double));
      for (j=0; j<ndim; j++) {
	 if (i==j) {tmp[j] = sp_1D_d[j];}
         else {tmp[j] = sp_1D[j];}
      }
      spline_kron(ndim,sdim,tmp,spline);
      memen_AtAb(ndim,cdim,ksp,sdim,spline,beta,0,&Ab_sz,&Ab);
      for (j=0; j<Ab_sz; j++) {memen += Ab[j]*Ab[j];}
      my_free(Ab);
   }

   my_free(spline);
   for (i=0; i<ndim; i++) {my_free(sp_1D[i]); my_free(sp_1D_d[i]);}
                  
   return(memen);
}
double *memen_AtAb(/* Input. */
                   int           ndim,     /* # of dimensions (1,2 or 3). */
                   const int     cdim[3],  /* Size of coeffient array. */
                   const int     ksp[3],   /* Knot-spacing. */
                   const int     sdim[3],  /* Kernel dimensions. */
                   const double  *spl,     /* Spline kernel. */
                   const double  *beta,    /* Spline coefficients. */
                   int           what,     /* 0->A*b, 1->A'*A*b */
                   /* Output. */
                   int           *sz,      /* Length of output vector */
                   double        **ovec)   /* Output vector */
{
   int    x=0, y=0, z=0;
   int    xs=0, ys=0, zs=0;
   int    xoff=0, yoff=0, zoff=0;
   int    c=0, i=0;
   int    offs=0;
   int    ksz = 1;       /* Total size of spline kernel. */
   int    eidim[3];      /* Size of full field. */
   int    eisz = 1;      /* Total size of field. */
   int    csz = 1;       /* Total size of coefficient matrix. */
   int    *indx = NULL;  /* Indicies into field for first spline. */
   double *Ab = NULL;
   double *AtAb = NULL;

   for (i=0, eisz=1, ksz=1; i<ndim; i++) 
   {
      eidim[i] = ksp[i]*(cdim[i]+3) - 1;
      eisz *= eidim[i];
      ksz *= sdim[i];
      csz *= cdim[i];
   }
   for (i=ndim; i<3; i++) {eidim[i] = 1;}

   Ab = (double *) my_calloc(eisz,sizeof(double));
   indx = (int *) my_calloc(ksz,sizeof(int));

   /* 
   First get a range of indicies into Abeta
   for the first ((0,0,0) in coef-space)
   spline kernel.
   */
   for (z=0, i=0; z<(ndim>2 ? 4*ksp[2]-1 : 1); z++) {
      for (y=0; y<(ndim>1 ? 4*ksp[1]-1 : 1); y++) {
	 for (x=0; x<4*ksp[0]-1; x++, i++) {
	    indx[i] = index(x,y,z,eidim);
         }
      }
   }

   /*
   Then calculate the index offset as we go
   one step in the x-, y- and z-direction
   respectively on the coefficient space.
   */
   xs = ksp[0];
   ys = ksp[1]*eidim[0];
   zs = ksp[2]*eidim[0]*eidim[1];

   /*
   Then build A*b as weighted sum
   of the columns of A.
   */
   for (z=0, zoff=0, c=0; z<cdim[2]; z++, zoff+=zs) {
      for (y=0, yoff=0; y<cdim[1]; y++, yoff+=ys) {
	 for (x=0, xoff=0; x<cdim[0]; x++, xoff+=xs, c++) {
	    offs = xoff+yoff+zoff;
            for (i=0; i<ksz; i++) {
	       Ab[indx[i]+offs] += beta[c] * spl[i];
            }
         }
      }
   }

   /* 
   Go on to build A'*A'b if that is what
   we are asked to do.
   */
   if (what == 1)
   {
      *sz = csz;
      AtAb = (double *) my_calloc(csz,sizeof(double));
      for (z=0, zoff=0, c=0; z<cdim[2]; z++, zoff+=zs) {
         for (y=0, yoff=0; y<cdim[1]; y++, yoff+=ys) {
	    for (x=0, xoff=0; x<cdim[0]; x++, xoff+=xs, c++) {
	       offs = xoff+yoff+zoff;
               for (i=0; i<ksz; i++) {
	          AtAb[c] += spl[i] * Ab[indx[i]+offs];
               }
            }
         }
      }
      *ovec = AtAb;
      my_free(Ab);
   }
   else {*sz = eisz; *ovec = Ab;}
               
   my_free(indx);

   return(*ovec);
}
/* 
Calculates the H-matrix such that b'*H*b gives the membrane energy for the field given
by the spline coefficients b. The reason we are not using make_AtA for this is that that
would give the energy for a field extending exactly over the image FOV, which would give
very low weight to splines at edges and corners. What we really want is the energy of the
field that extends all the way out to the edges of the outermost splines.

The routine as it stands is ludicrously inefficient. This is an effect of lazily reusing
code from make_AtA. It would be possible to speed it up a lot.
*/
int make_memen_H(/* Input. */
                 const int         ndim,      /* Actual dimensionality of problem (1, 2 or 3). */
                 int               cdim[3],   /* # of knots in the three dimensions. */
                 const int         ksp[3],    /* Knot spacing in the three dimensions. */
                 /* Output. */
                 int               *irp,            /* Array of row-indicies. */
                 int               *jcp,            /* Array of pointers into column-starts in irp. */
                 double            *sp)             /* Array of non-zero values in sparse matric. */
{
   int     i = 0, j = 0;
   int     ci = 0;              /* Column index for H. */
   int     ri = 0;              /* Row index for H. */
   int     cntr = 0;
   int     kns=0, jns=0, ins=0; /* Start of neighbours in k, j and i directions. */
   int     kne=0, jne=0, ine=0; /* End of neighbours in k, j and i directions. */
   int     s1k=0, s1j=0, s1i=0; /* index for "first" spline */
   int     s2k=0, s2j=0, s2i=0; /* index for "second" spline */
   int     sdim[3]={0,0,0};     /* Size of spline kernel in the three dimensions. */
   int     memsize = 0;
   double  *sp_1D[3] = {NULL, NULL, NULL};
   double  *splines[3] = {NULL, NULL, NULL};

   /* Create the spline kernels that we'll need. */
   for (i=0; i<ndim; i++)
   {
      /* Make as many 1D splines as we need. */
      memsize = 1;
      for (j=0; j<ndim; j++)
      {
	 if (i==j) {sdim[j] = get_1D_spline_d(ksp[j],&(sp_1D[j])); memsize *= sdim[j];} /* Get first derivative. */
	 else {sdim[j] = get_1D_spline(ksp[j],&(sp_1D[j])); memsize *= sdim[j];}        /* Get straight spline. */
      }
      splines[i] = (double *) my_calloc(memsize,sizeof(double));
      spline_kron(ndim,sdim,sp_1D,splines[i]);
      for (j=0; j<ndim; j++) 
      {
         my_free(sp_1D[j]);
      }
   }

   for (i=ndim; i<3; i++) {cdim[i]=1; sdim[i]=1;}

   
   /* First three loops over splines determine the column of H. */
   ci = 0;
   cntr = 0;
   for (s1k=0; s1k<cdim[2]; s1k++)
   {
      get_nabos(s1k,cdim[2],sdim[2],&kns,&kne);
      for (s1j=0; s1j<cdim[1]; s1j++) 
      {
	 get_nabos(s1j,cdim[1],sdim[1],&jns,&jne);
	 for (s1i=0; s1i<cdim[0]; s1i++)
	 {
  	    get_nabos(s1i,cdim[0],sdim[0],&ins,&ine);
            jcp[ci] = cntr;
            /* Loop over all neighbours. */
            for (s2k=kns; s2k<kne; s2k++)
	    {
	       for (s2j=jns; s2j<jne; s2j++)
	       {
		  for (s2i=ins; s2i<ine; s2i++)
		  {
		     ri = index(s2i,s2j,s2k,cdim);
                     irp[cntr] = ri;
                     if (ri < ci) /* If above diagonal, utilise symmetry. */
		     {
		        sp[cntr++] = find_val(&(irp[jcp[ri]]),jcp[ri+1]-jcp[ri],ci,&(sp[jcp[ri]]));
                     }
                     else
		     {
		        sp[cntr] = 0.0;
                        for (i=0; i<ndim; i++)
			{
                           sp[cntr] += dot_prod_H(s1i,s1j,s1k,s2i,s2j,s2k,sdim,splines[i],splines[i]);
                        }
                        cntr++;
                     }
                  }
	       }
	    }
            ci++;
         }
      }
   }
   jcp[ci] = cntr;   
   
   for (i=0; i<ndim; i++) {my_free(splines[i]);}

   return(cntr);
}

int make_AtA(/* Input. */
             int               ndim,      /* Actual dimensionality of problem (1, 2 or 3). */
             int               cdim[3],   /* # of knots in the three dimensions. */
             int               sdim[3],   /* Size of spline kernel in the three dimensions. */
             double            *spl,      /* Spline kernel. */
             int               idim[3],   /* Size of image matrix. */
             double            *ima,      /* Image. */
             /* Output. */
             int               *irp,      /* Array of row-indicies. */
             int               *jcp,      /* Array of pointers into column-starts in irp. */
             double            *sp)       /* Array of non-zero values in sparse matric. */
{
   int     i = 0;
   int     ci = 0;              /* Column index for AtA. */
   int     ri = 0;              /* Row index for AtA. */
   int     cntr = 0;
   int     kns=0, jns=0, ins=0; /* Start of neighbours in k, j and i directions. */
   int     kne=0, jne=0, ine=0; /* End of neighbours in k, j and i directions. */
   int     s1k=0, s1j=0, s1i=0; /* index for "first" spline */
   int     s2k=0, s2j=0, s2i=0; /* index for "second" spline */
   int     is[3]={0,0,0};       /* Start indicies in image for first spline. */
   int     ie[3]={0,0,0};       /* End indicies in image for first spline. */
   int     ss[3]={0,0,0};       /* Offset for first spline. */
   double  *s_by_i = NULL;      /* Spline multiplied by appurtenant image intensities. */
   double  tmp = 0.0;

   for (i=ndim; i<3; i++) {cdim[i]=1; sdim[i]=1; idim[i]=1;}

   s_by_i = (double *) my_calloc(sdim[0]*sdim[1]*sdim[2],sizeof(double));
   
   /* First three loops over splines determine the column of AtA. */
   ci = 0;
   cntr = 0;
   for (s1k=0; s1k<cdim[2]; s1k++)
   {
      get_nabos(s1k,cdim[2],sdim[2],&kns,&kne);
      for (s1j=0; s1j<cdim[1]; s1j++) 
      {
	 get_nabos(s1j,cdim[1],sdim[1],&jns,&jne);
	 for (s1i=0; s1i<cdim[0]; s1i++)
	 {
  	    get_nabos(s1i,cdim[0],sdim[0],&ins,&ine);
            jcp[ci] = cntr;
            /* Pre-multiply this spline with image intensities. */
	    get_s_by_i(s1i,s1j,s1k,sdim,idim,spl,ima,s_by_i);
            /* And get indicies into image array. */
            get_range(s1i,sdim[0],idim[0],&(is[0]),&(ie[0]),&(ss[0]));
            get_range(s1j,sdim[1],idim[1],&(is[1]),&(ie[1]),&(ss[1]));
            get_range(s1k,sdim[2],idim[2],&(is[2]),&(ie[2]),&(ss[2]));
            /* Loop over all neighbours. */
            for (s2k=kns; s2k<kne; s2k++)
	    {
	       for (s2j=jns; s2j<jne; s2j++)
	       {
		  for (s2i=ins; s2i<ine; s2i++)
		  {
		     ri = index(s2i,s2j,s2k,cdim);
                     irp[cntr] = ri;
                     if (ri < ci) /* If above diagonal, utilise symmetry. */
		     {
		        sp[cntr++] = find_val(&(irp[jcp[ri]]),jcp[ri+1]-jcp[ri],ci,&(sp[jcp[ri]]));
                     }
                     else
		     {
                        tmp = dot_prod(s1i,s1j,s1k,s2i,s2j,s2k,sdim,
			               idim,spl,ima,s_by_i,is,ie,ss);
                        if (tmp) {sp[cntr++] = tmp;}
                        else {sp[cntr++] = DBL_MIN;}  /* Really really silly. */
                     }
                  }
	       }
	    }
            ci++;
         }
      }
   }
   jcp[ci] = cntr;   

   my_free(s_by_i);

   return(cntr);
}

double find_val(int     *a,
                int     n,
                int     key,
                double  *val)
{
   int   j = 0;
   int   jlo = -1;     /* Just had to! */
   int   jup = n;

   if (key < a[0] || key > a[n-1]) {return(0.0);}

   while (jup-jlo > 1)
   {
      j = (jlo+jup) >> 1;
      if (key >= a[j]) {jlo = j;}
      else {jup = j;}
   }

   if (a[jlo] == key) {return(val[jlo]);}
   else return(0.0);
}

int make_AtB(/* Input. */
             int               ndim,      /* Actual dimensionality of problem (1, 2 or 3). */
             int               cdim[3],   /* # of knots in the three dimensions. */
             int               sdim[3],   /* Size of spline kernel in the three dimensions. */
             double            *splB,     /* Spline kernel for B. */
             double            *splA,     /* Spline kernel for A. */
             int               idim[3],   /* Size of image matrix. */
             double            *imaB,     /* Image for B. */
             double            *imaA,     /* Image for A. */
             /* Output. */
             int               *irp,      /* Array of row-indicies. */
             int               *jcp,      /* Array of pointers into column-starts in irp. */
             double            *sp)       /* Array of non-zero values in sparse matric. */
{
   int     i = 0;
   int     ci = 0;              /* Column index for AtB. */
   int     ri = 0;              /* Row index for AtB. */
   int     cntr = 0;
   int     kns=0, jns=0, ins=0; /* Start of neighbours in k, j and i directions. */
   int     kne=0, jne=0, ine=0; /* End of neighbours in k, j and i directions. */
   int     s1k=0, s1j=0, s1i=0; /* index for "first" spline */
   int     s2k=0, s2j=0, s2i=0; /* index for "second" spline */
   int     is[3]={0,0,0};       /* Start indicies in image for first spline. */
   int     ie[3]={0,0,0};       /* End indicies in image for first spline. */
   int     ss[3]={0,0,0};       /* Offset for first spline. */
   double  *s_by_i = NULL;      /* Spline multiplied by appurtenant image intensities. */
   double  tmp = 0.0;

   for (i=ndim; i<3; i++) {cdim[i]=1; sdim[i]=1; idim[i]=1;}

   s_by_i = (double *) my_calloc(sdim[0]*sdim[1]*sdim[2],sizeof(double));
   
   /* First three loops over splines determine the column of AtB. */
   ci = 0;
   cntr = 0;
   for (s1k=0; s1k<cdim[2]; s1k++)
   {
      get_nabos(s1k,cdim[2],sdim[2],&kns,&kne);
      for (s1j=0; s1j<cdim[1]; s1j++) 
      {
	 get_nabos(s1j,cdim[1],sdim[1],&jns,&jne);
	 for (s1i=0; s1i<cdim[0]; s1i++)
	 {
  	    get_nabos(s1i,cdim[0],sdim[0],&ins,&ine);
            jcp[ci] = cntr;
            /* Pre-multiply spline B with image B intensities. */
	    get_s_by_i(s1i,s1j,s1k,sdim,idim,splB,imaB,s_by_i);
            /* And get indicies into image array. */
            get_range(s1i,sdim[0],idim[0],&(is[0]),&(ie[0]),&(ss[0]));
            get_range(s1j,sdim[1],idim[1],&(is[1]),&(ie[1]),&(ss[1]));
            get_range(s1k,sdim[2],idim[2],&(is[2]),&(ie[2]),&(ss[2]));
            /* Loop over all neighbours. */
            for (s2k=kns; s2k<kne; s2k++)
	    {
	       for (s2j=jns; s2j<jne; s2j++)
	       {
		  for (s2i=ins; s2i<ine; s2i++)
		  {
		     ri = index(s2i,s2j,s2k,cdim);
                     irp[cntr] = ri;
                     tmp = dot_prod(s1i,s1j,s1k,s2i,s2j,s2k,sdim,idim,splA,imaA,s_by_i,is,ie,ss);
                     if (tmp) {sp[cntr++] = tmp;}
                     else {sp[cntr++] = DBL_MIN;}
                  }
	       }
	    }
            ci++;
         }
      }
   }
   jcp[ci] = cntr;   

   my_free(s_by_i);

   return(cntr);
}

double dot_prod(int       i1,         /* [i1 j1 k1] is index of first spline kernel/coef. */ 
                int       j1,
                int       k1,
                int       i2,         /* [i2 j2 k2] is index of second spline kernel/coef. */ 
                int       j2,
                int       k2,
                int       sdim[3],    /* Size of spline kernel (common to s1 and s2). */
                int       idim[3],    /* Size of image. */
                double    *s2,        /* Second spline kernel. */
                double    *ima,       /* Image. */
                double    *s1,        /* First spline kernel (pre-multiplied with image. */
                int       is1[3],     /* Start indices of first spline kernel in image. */
                int       ie1[3],     /* End indices of first spline kernel in image. */
                int       ss1[3])     /* Offset into first spline kernel (to handle edges). */
{
   int   is2i=0, is2j=0, is2k=0;    /* Start indicies into ima. */
   int   ie2i=0, ie2j=0, ie2k=0;    /* End indicies into ima. */
   int   ss2i=0, ss2j=0, ss2k=0;    /* Start indicies into s2 (second spline kernel). */
   int   ss1i=0, ss1j=0, ss1k=0;    /* Start indicies into s1 (first spline kernel pre-multiplied with ima). */
   int   i2i=0, i2j=0, i2k=0;       /* Indicies into ima. */
   int   s2i=0, s2j=0, s2k=0;       /* Indicies into s2. */
   int   s1i=0, s1j=0, s1k=0;       /* Indicies into s1. */
   double   dp = 0.0;               /* Dot product. */  

   get_range(i2,sdim[0],idim[0],&is2i,&ie2i,&ss2i);
   get_range(j2,sdim[1],idim[1],&is2j,&ie2j,&ss2j);
   get_range(k2,sdim[2],idim[2],&is2k,&ie2k,&ss2k);

   /* 
      Get indicies for overlap 
   */

   if (is2i>is1[0]) {ss1i = ss1[0] + (is2i-is1[0]);}
   else if (is1[0]>is2i) {ss1i = ss1[0]; ss2i += (is1[0]-is2i); is2i = is1[0];}
   else {ss1i = ss1[0];} 
   ie2i = MIN(ie2i,ie1[0]); 
   
   if (is2j>is1[1]) {ss1j = ss1[1] + (is2j-is1[1]);}
   else if (is1[1]>is2j) {ss1j = ss1[1]; ss2j += (is1[1]-is2j); is2j = is1[1];}
   else {ss1j = ss1[1];} 
   ie2j = MIN(ie2j,ie1[1]); 
   
   if (is2k>is1[2]) {ss1k = ss1[2] + (is2k-is1[2]);}
   else if (is1[2]>is2k) {ss1k = ss1[2]; ss2k += (is1[2]-is2k); is2k = is1[2];}
   else {ss1k = ss1[2];} 
   ie2k = MIN(ie2k,ie1[2]); 

   dp = 0.0;
   for (i2k=is2k,s2k=ss2k,s1k=ss1k; i2k<ie2k; i2k++,s2k++,s1k++)
   {
      for (i2j=is2j,s2j=ss2j,s1j=ss1j; i2j<ie2j; i2j++,s2j++,s1j++)
      {
         for (i2i=is2i,s2i=ss2i,s1i=ss1i; i2i<ie2i; i2i++,s2i++,s1i++)
         {
	    dp += s2[index(s2i,s2j,s2k,sdim)] * ima[index(i2i,i2j,i2k,idim)] * s1[index(s1i,s1j,s1k,sdim)];
         }
      }
   }

   return(dp);
}
                
double dot_prod_H(int       i1,         /* [i1 j1 k1] is index of first spline kernel/coef. */ 
                  int       j1,
                  int       k1,
                  int       i2,         /* [i2 j2 k2] is index of second spline kernel/coef. */ 
                  int       j2,
                  int       k2,
                  int       sdim[3],    /* Size of spline kernel (common to s1 and s2). */
                  double    *s1,        /* First spline kernel. */
                  double    *s2)        /* Second spline kernel. */
{
   int     ss1i=0, ss1j=0, ss1k=0;    /* Start indicies into s1 (first spline kernel). */
   int     se1i=0, se1j=0, se1k=0;    /* End indicies into s1. */
   int     ss2i=0, ss2j=0, ss2k=0;    /* Start indicies into s2 (second spline kernel). */
   int     s2i=0, s2j=0, s2k=0;       /* Indicies into s2. */
   int     s1i=0, s1j=0, s1k=0;       /* Indicies into s1. */
   double  dp = 0.0;                  /* Dot product. */  

   if (i1<i2) {se1i=sdim[0]; ss2i=0; ss1i=(i2-i1)*sz2ksp(sdim[0]); }
   else if (i2<i1) {ss1i=0; se1i=sdim[0]-(i1-i2)*sz2ksp(sdim[0]); ss2i=(i1-i2)*sz2ksp(sdim[0]);}
   else {ss1i=ss2i=0; se1i=sdim[0];}
   
   if (j1<j2) {se1j=sdim[1]; ss2j=0; ss1j=(j2-j1)*sz2ksp(sdim[1]); }
   else if (j2<j1) {ss1j=0; se1j=sdim[1]-(j1-j2)*sz2ksp(sdim[1]); ss2j=(j1-j2)*sz2ksp(sdim[1]);}
   else {ss1j=ss2j=0; se1j=sdim[1];}
   
   if (k1<k2) {se1k=sdim[2]; ss2k=0; ss1k=(k2-k1)*sz2ksp(sdim[2]); }
   else if (k2<k1) {ss1k=0; se1k=sdim[2]-(k1-k2)*sz2ksp(sdim[2]); ss2k=(k1-k2)*sz2ksp(sdim[2]);}
   else {ss1k=ss2k=0; se1k=sdim[2];}
   

   dp = 0.0;
   for (s2k=ss2k,s1k=ss1k; s1k<se1k; s2k++,s1k++)
   {
      for (s2j=ss2j,s1j=ss1j; s1j<se1j; s2j++,s1j++)
      {
         for (s2i=ss2i,s1i=ss1i; s1i<se1i; s2i++,s1i++)
         {
	    dp += s2[index(s2i,s2j,s2k,sdim)] * s1[index(s1i,s1j,s1k,sdim)];
         }
      }
   }

   return(dp);
}
                
double get_s_by_i(/* Input */
                  int      i,         /* [i j k] index of spline kernel/coef. */
                  int      j,
                  int      k,
                  int      sdim[3],   /* Size of spline kernel. */
                  int      idim[3],   /* Size of image. */
                  double   *spl,      /* Spline. */
                  double   *ima,      /* Image. */
                  /* Output */
                  double   *sbyi)     /* Spline multiplied with appuretenant values in image. */
{    
   int     ii=0, ij=0, ik=0;     /* Indicies into image. */     
   int     si=0, sj=0, sk=0;     /* Indicies into spline. */
   int     sindex=0;             /* linear index into spline. */
   int     iks=0, ijs=0, iis=0;  /* Start indicies into image. */
   int     ike=0, ije=0, iie=0;  /* End indicies into image. */
   int     sks=0, sjs=0, sis=0;  /* Start indicies into spline. */
   double  sum=0.0;              /* Sum over valid range of spline. */

   get_range(i,sdim[0],idim[0],&iis,&iie,&sis);
   get_range(j,sdim[1],idim[1],&ijs,&ije,&sjs);
   get_range(k,sdim[2],idim[2],&iks,&ike,&sks);
   memset(sbyi,0,sdim[0]*sdim[1]*sdim[2]*sizeof(double));

   sum = 0.0;
   for (ik=iks,sk=sks; ik<ike; ik++,sk++)
   {
      for (ij=ijs,sj=sjs; ij<ije; ij++,sj++)
      {
         for (ii=iis,si=sis; ii<iie; ii++,si++)
	 {
	    sindex = index(si,sj,sk,sdim);
	    sbyi[sindex] = ima[index(ii,ij,ik,idim)] * spl[sindex];
            sum += sbyi[sindex];
         }
      }
   }
			      
   return(sum);
}

int get_nabos(/* Input. */
              int   i,       /* Index (in one dimension) of spline-coef. */
              int   csz,     /* Total # of spline-coef. */
              int   ssz,     /* Size of spline kernel. */
              /* Output. */
              int   *ns,     /* Lowest index of overlapping neighbour. */
              int   *ne)     /* Highest index of overlapping neighbour. */
{
   int   ol=0;

   ol = MIN(sz2ksp(ssz),3);
   if ((*ns=i-ol) < 0) {*ns = 0;}
   if ((*ne=i+ol+1) > csz) {*ne = csz;}

   return(*ns); /* As an extra courtesy. */  
}

int get_A_nzmax(/* Input. */
                int      ndim,
                int      kdim[3],
                int      sdim[3],
                int      idim[3])
{
   int     nzmax = 0;
   int     i=0;
   int     kk=0, kj=0, ki=0;
   int     iks=0, ijs=0, iis=0;
   int     ike=0, ije=0, iie=0;
   int     sks=0, sjs=0, sis=0;

   for (i=ndim; i<3; i++) {kdim[i]=1; sdim[i]=1; idim[i]=1;}

   for (kk=0; kk<kdim[2]; kk++)
   {
      get_range(kk,sdim[2],idim[2],&iks,&ike,&sks);
      for (kj=0; kj<kdim[1]; kj++)
      {
	 get_range(kj,sdim[1],idim[1],&ijs,&ije,&sjs);
	 for (ki=0; ki<kdim[0]; ki++)
	 {
            get_range(ki,sdim[0],idim[0],&iis,&iie,&sis);
            nzmax += (iie-iis)*(ije-ijs)*(ike-iks);
         }
      }
   }
   
   return(nzmax);
}

int get_AtA_nzmax(int          ndim,
                  const int    nknot[3],
                  const int    ksp[3])
{
   int    lnknot[3], lksp[3];
   int    nzmax=0;
   int    n_knabo=0, n_jnabo=0, n_inabo=0;
   int    i=0, j=0, k=0;

   for (i=0; i<ndim; i++) {lnknot[i]=nknot[i]; lksp[i]=ksp[i];}
   for (i=ndim; i<3; i++) {lnknot[i]=1; lksp[i]=1;}

   nzmax = 0;
   for (k=0; k<lnknot[2]; k++)
   {
      n_knabo = n_nabo(k,lnknot[2],lksp[2]);
      for (j=0; j<lnknot[1]; j++)
      {
	 n_jnabo = n_nabo(j,lnknot[1],lksp[1]);
         for (i=0; i<lnknot[0]; i++)
	 {
	    n_inabo = n_nabo(i,lnknot[0],lksp[0]);
            nzmax += n_knabo*n_jnabo*n_inabo;
         }
      }
   }   
   return(nzmax);
}
 
int n_nabo(int   i,
           int   n,
           int   ksp)
{
   int   nn = 0;

   ksp = MIN(ksp,3);
   nn = 2*ksp + 1;

   if (i < ksp) nn -= (ksp - i);
   if ((n-1-i) < ksp) nn -= (ksp - (n-1-i));

   return(nn);     
}


int AtimesB(/* Input. */
            int      *ir_A,
            int      *jc_A,
            double   *s_A,
            int      mA,
            int      *ir_B,
            int      *jc_B,
            double   *s_B,
            int      nB,
            int      nzmax,
            /* Output. */
            int      **ir_out_orig, /* These have to be pointers */
            int      *jc_out,       /* to pointers to allow for  */
            double   **s_out_orig)  /* realloc.                  */
{
   int      i=0, j=0;
   int      si=0, ei=0;
   int      si2=0, ei2=0;
   int      ndx = 0;         /* Scratch. */
   int      ci = 0;          /* Column index for A*B */
   int      cnt = 0;         /* Total # of non-zero elements. */ 
   int      nc = 0;          /* # of non-zero elements in present column. */
   double   bval = 0.0;      /* Scratch. */
   int      *ir_out = NULL; 
   double   *s_out = NULL; 
   int      *ir_tmp = NULL;  /* ir_tmp, full _ir and full_s plays the role of SPA   */
   int      *full_ir = NULL; /* as described in Gilbert, Moler & Schreiber 1992.    */
   double   *full_s = NULL;

   ir_out = *ir_out_orig;
   s_out = *s_out_orig;
   full_s = (double *) my_calloc(mA,sizeof(double));
   full_ir = (int *) my_calloc(mA,sizeof(int));
   ir_tmp = (int *) my_calloc(mA,sizeof(int));

   cnt = 0;
   for (ci=0; ci<nB; ci++)  /* Loop over all columns of A*B (and hence B). */
   {
      jc_out[ci] = cnt;
      nc = 0; 
      si = jc_B[ci]; ei = jc_B[ci+1];
      for (i=si; i<ei; i++) /* Loop over all non-zero values in B column. */
      {
	 bval = s_B[i];
	 si2 = jc_A[ir_B[i]]; ei2 = jc_A[ir_B[i]+1];
         for (j=si2; j<ei2; j++) /* And multiply with all non-zero values in A column. */
	 {
	    ndx = ir_A[j];
	    full_s[ndx] += bval * s_A[j];
            if (!full_ir[ndx])
	    {
	       full_ir[ndx] = 1;
               ir_tmp[nc++] = ndx;
            }
         }
      }
      /* Sort ir_tmp prior to copying to ir_out. */
      qsort(ir_tmp,nc,sizeof(int),cmpf); 

      /* Check we have room for new column. */
      if ((cnt+nc) > nzmax)
      {
	 nzmax = ((int) 1.1 * nzmax);  /* Increase memory by 10% */
         ir_out = (int *) my_realloc(ir_out,nzmax*sizeof(int));
         s_out = (double *) my_realloc(s_out,nzmax*sizeof(double));
         printf("\nWarning, non-optimal nzmax passed to AtranspA.");
      }
 
      /* Put values in ir_out and s_out. */
      memcpy(&(ir_out[cnt]),ir_tmp,nc*sizeof(int)); 
      for (i=cnt; i<(cnt+nc); i++) 
      {
	 ndx = ir_out[i];
         s_out[i] = full_s[ndx];
         full_s[ndx] = 0.0;
         full_ir[ndx] = 0;
      }
      cnt += nc;
   }
   jc_out[ci] = cnt;

   *ir_out_orig = ir_out;
   *s_out_orig = s_out;
 
   my_free(full_s);
   my_free(full_ir);
   my_free(ir_tmp);

   return(nzmax);             
}


int cmpf(const void    *el1,
         const void    *el2)
{
   if (*((int *)el1) < *((int *)el2)) return(-1);
   else if (*((int *)el1) > *((int *)el2)) return(1);
   else return(0); 
}


int AtranspA(/* Input. */
             int      *ir_in,
             int      *jc_in,
             double   *s_in,
             int      m,
             int      n,
             int      nzmax,
             /* Output. */
             int      **ir_out_orig,/* These have to be pointers */
             int      *jc_out,      /* to pointers to allow for  */
             double   **s_out_orig) /* realloc.                  */
{
   int     i=0, j=0;
   int     ri=0;          /* Row index of AtA. */
   int     ci=0;          /* Column index of AtA. */
   int     cnt=0;         /* Count of non-zero elements in AtA. */
   int     si=0, ei=0;
   int     si2=0, ei2=0;
   double  *s_tmp = NULL; /* Temprary non-sparse column of A. */
   int     *flags = NULL; /* Temprary non-sparse column of A. */
   int     *ir_out = NULL;
   double  *s_out = NULL;

   /* 
      We are working with copies of the input pointers
      to avoid double dereferencing at each access.
   */

   ir_out = *ir_out_orig;
   s_out = *s_out_orig;

   s_tmp = s_in - 1;
   flags = (int *) my_calloc(m,sizeof(int));

   for (ci=0,cnt=0; ci<n; ci++) /* Loop over all columns of AtA. */
   {
      jc_out[ci] = cnt;
      si = jc_in[ci];
      ei = jc_in[ci+1];
      if (ei > si) /* If anything non-zero at all in this column. */
      {
	 /* Fill in dense copy of column of A. */
	 for (i=si; i<ei; i++) {flags[ir_in[i]] = i+1;}

         for (ri=0; ri<n; ri++) /* Loop over all rows of AtA for this column. */
         {
	    if (cnt >= nzmax) /* Oh oh, didn't really want that. */
	    {
	       nzmax = ((int) 1.1 * nzmax);  /* Increase memory by 10% */
               ir_out = (int *) my_realloc(ir_out,nzmax*sizeof(int));
               s_out = (double *) my_realloc(s_out,nzmax*sizeof(double));
               printf("\nWarning, non-optimal nzmax passed to AtranspA.");
            }

	    /* Divide into diagonal, sub-diagonal and super-diagonal cases. */
	    if (ri==ci) /* If diagonal. */
	    {

	       for (i=si; i<ei; i++) /* Dont need to use dense copy for this case. */
	       {
		  s_out[cnt] += s_in[i]*s_in[i];
	       }
               ir_out[cnt++] = ri;
            }
            else if (ri > ci) /* If subdiagonal, calculate it. */
	    {
	       si2 = jc_in[ri];
               ei2 = jc_in[ri+1];
               for (i=si2; i<ei2; i++) /* Use dense column to calculate dot-product. */
	       {
		  if ((j = flags[ir_in[i]]))
		  {
		     s_out[cnt] += s_in[i] * s_tmp[j];
                  }
               }
               if (s_out[cnt]) {ir_out[cnt++] = ri;}
            }
            else /* Super-diagonal, reflect it. */
	    {
	       si2 = jc_out[ri];
	       ei2 = jc_out[ri+1];
               if ((s_out[cnt] = find_val(&(ir_out[si2]),ei2-si2,ci,&(s_out[si2]))))
	       {
		  ir_out[cnt++] = ri;
               }
            }
         } 
	 /* Reset dense copy of column of A. */
	 for (i=si; i<ei; i++) {flags[ir_in[i]] = 0;}
      }
   }

   jc_out[ci] = cnt;

   *ir_out_orig = ir_out;
   *s_out_orig = s_out;

   my_free(flags);

   return(nzmax);
}

int AtranspB(/* Input. */
             int      *ir_A,
             int      *jc_A,
             double   *s_A,
             int      mA,
             int      nA,
             int      *ir_B,
             int      *jc_B,
             double   *s_B,
             int      mB,
             int      nB,
             int      nzmax,
             /* Output. */
             int      **ir_out_orig, /* These have to be pointers */
             int      *jc_out,       /* to pointers to allow for  */
             double   **s_out_orig)  /* realloc.                  */
{
   int     i=0, j=0;
   int     ri=0;          /* Row index of AtB. */
   int     ci=0;          /* Column index of AtB. */
   int     cnt=0;         /* Count of non-zero elements in AtB. */
   int     si=0, ei=0;
   int     si2=0, ei2=0;
   double  *s_tmp = NULL; /* Temporary non-sparse column of A. */
   int     *flags = NULL; /* Temporary non-sparse column of A. */
   int     *ir_out = NULL;
   double  *s_out = NULL;

   if (mA != mB) {return(-1);} /* mA is n of A' */

   /* 
      We are working with copies of the input pointers
      to avoid double dereferencing at each access.
   */

   ir_out = *ir_out_orig;
   s_out = *s_out_orig;

   s_tmp = s_B - 1;
   flags = (int *) my_calloc(mB,sizeof(int));

   for (ci=0,cnt=0; ci<nB; ci++) /* Loop over all columns of AtB. */
   {
      jc_out[ci] = cnt;
      si = jc_B[ci];
      ei = jc_B[ci+1];
      if (ei > si) /* If anything non-zero at all in this column. */
      {
	 /* Fill in dense copy of column of B. */
	 for (i=si; i<ei; i++) {flags[ir_B[i]]=i+1;}

         for (ri=0; ri<nA; ri++) /* Loop over all rows of AtB for this column. */
         {
	    if (cnt >= nzmax) /* Oh oh, didn't really want that. */
	    {
	       nzmax = ((int) 1.1 * nzmax);  /* Increase memory by 10% */
               ir_out = (int *) my_realloc(ir_out,nzmax*sizeof(int));
               s_out = (double *) my_realloc(s_out,nzmax*sizeof(double));
               printf("\nWarning, non-optimal nzmax passed to AtranspA.");
            }

	    si2 = jc_A[ri];
            ei2 = jc_A[ri+1];
	    if ((ir_B[si] <= ir_A[si2] && ir_B[ei-1] >= ir_A[si2]) ||
		(ir_A[si2] <= ir_B[si] && ir_A[ei2-1] >= ir_B[si]))
	    {
               for (i=si2; i<ei2; i++) /* Calculate dot-product. */
	       {
		  if ((j = flags[ir_A[i]]))
		  {
                     s_out[cnt] += s_A[i] * s_tmp[j];
                  }
	       }
               if (s_out[cnt]) {ir_out[cnt++] = ri;}
	    } 
         } 
         for (i=si; i<ei; i++) {flags[ir_B[i]] = 0;}
      }
   }

   jc_out[ci] = cnt;

   *ir_out_orig = ir_out;
   *s_out_orig = s_out;

   my_free(flags);

   return(nzmax);
}

/*

These routines "shadow" zoom_field and zoom_field_by2,
but take the "true" old coefficient matrix size as input
parameter. 

*/

/*
Calculates the new spline-coefficients when a field is zoomed
by a factor that is a power of 2 in one or more directions.
*/
int fnirt_zoom_field(/* input */
                     int      ndim,
                     int      oksp[3],
                     int      ocdim[3],
                     int      nksp[3],
                     int      idim[3],
                     double   *oc,
                     /* Output */
                     double   *nc)
{
   int      i=0;
   int      dim=0;
   int      sf=0;
   int      sz=0;
   int      tksp[3];
   int      tcdim[3];
   double   *tmpc1=NULL;
   double   *tmpc2=NULL;

   for (i=ndim; i<3; i++) {oksp[i]=1; ocdim[i] = 1; nksp[i]=1; idim[i]=1;}
   for (i=0; i<3; i++) {tksp[i] = oksp[i]; tcdim[i] = ocdim[i];}

   tmpc1 = oc;
   
   for (dim=0; dim<3; dim++)
   {
      sf = oksp[dim]/nksp[dim];
      while (sf > 1)
      {
	 fnirt_zoom_field_by2(3,tksp,tcdim,idim,dim,tmpc1,&tmpc2);
         tksp[dim] /= 2;
         tcdim[dim] = no_of_knots(tksp[dim],idim[dim]);
         sf /= 2;
         if (tmpc1 != oc)
	 {
            my_free(tmpc1);
         }
         tmpc1 = tmpc2;
      }
      if (sf != 1)
      {
	 printf("\nzoom_field: zooming must be by power of 2.");
	 return(-1);
      }
   }

   for (i=0, sz=1; i<3; i++) {sz *= no_of_knots(nksp[i],idim[i]);}

   memcpy(nc,tmpc2,sz*sizeof(double));
   my_free(tmpc2);
 
   return(1);
}
/*
Calculates the new set of spline-coefficients resulting when the 
knot-spacing of an existing field is cut by half in one direction.
*/
int fnirt_zoom_field_by2(/* Input */
                         int     ndim,
                         int     ksp[3],
                         int     ocdim[3],
                         int     idim[3],
                         int     zdim,
                         double  *oc,
                         /* Output */
                         double  **nc)
{
   int      i=0, j=0, k=0;
   int      nsz=0;
   int      ncdim[3];
   double   a0=1.0/8.0, a4=a0;
   double   a1=1.0/2.0, a3=a1;
   double   a2=3.0/4.0;

   for (i=ndim; i<3; i++) {ocdim[i] = 1; idim[i] = 1; ksp[i] = 1;}
   for (i=0, nsz=1; i<3; i++) 
   {
      ncdim[i] = (i==zdim) ? no_of_knots(ksp[i]/2,idim[i]) : ocdim[i]; 
      nsz *= ncdim[i];
   }

   *nc = (double *) my_calloc(nsz,sizeof(double));

   if (zdim==0)
   {
      for (k=0; k<ncdim[2]; k++)
      {
	 for (j=0; j<ncdim[1]; j++)
	 {
	    for (i=0; i<ncdim[0]; i++) 
	    {
	       if (i % 2)
	       {
		  (*nc)[index(i,j,k,ncdim)] = a4*oc[index((i+1)/2-1,j,k,ocdim)] + 
		                              a2*oc[index((i+1)/2,j,k,ocdim)];
                  if (((i+1)/2+1) < ocdim[0]) {
                    (*nc)[index(i,j,k,ncdim)] += a0*oc[index((i+1)/2+1,j,k,ocdim)];
		  }
               }
               else
	       {
		  (*nc)[index(i,j,k,ncdim)] = a3*oc[index(i/2,j,k,ocdim)];
                  if ((i/2+1) < ocdim[0]) { 
		    (*nc)[index(i,j,k,ncdim)] += a1*oc[index(i/2+1,j,k,ocdim)];
		  } 
               }
	    }
	 }
      }
   }
   else if (zdim==1)
   {
      for (k=0; k<ncdim[2]; k++)
      {
	 for (j=0; j<ncdim[1]; j++)
	 {
	    for (i=0; i<ncdim[0]; i++) 
	    {
	       if (j % 2)
	       {
		  (*nc)[index(i,j,k,ncdim)] = a4*oc[index(i,(j+1)/2-1,k,ocdim)] + 
		                              a2*oc[index(i,(j+1)/2,k,ocdim)]; 
                  if (((j+1)/2+1) < ocdim[1]) {
                    (*nc)[index(i,j,k,ncdim)] += a0*oc[index(i,(j+1)/2+1,k,ocdim)];
		  }
               }
               else
	       {
		  (*nc)[index(i,j,k,ncdim)] = a3*oc[index(i,j/2,k,ocdim)];
                  if ((j/2+1) < ocdim[1]) { 
		    (*nc)[index(i,j,k,ncdim)] += a1*oc[index(i,j/2+1,k,ocdim)]; 
		  } 
               }
	    }
	 }
      }
   }
   else if (zdim==2)
   {
      for (k=0; k<ncdim[2]; k++)
      {
	 for (j=0; j<ncdim[1]; j++)
	 {
	    for (i=0; i<ncdim[0]; i++) 
	    {
	       if (k % 2)
	       {
		  (*nc)[index(i,j,k,ncdim)] = a4*oc[index(i,j,(k+1)/2-1,ocdim)] + 
		                              a2*oc[index(i,j,(k+1)/2,ocdim)];
                  if (((k+1)/2+1) < ocdim[2]) {
                    (*nc)[index(i,j,k,ncdim)] += a0*oc[index(i,j,(k+1)/2+1,ocdim)];
		  }
               }
               else
	       {
		  (*nc)[index(i,j,k,ncdim)] = a3*oc[index(i,j,k/2,ocdim)];
                  if ((k/2+1) < ocdim[2]) { 
		    (*nc)[index(i,j,k,ncdim)] += a1*oc[index(i,j,k/2+1,ocdim)]; 
		  } 
               }
	    }
	 }
      }
   }
   else
   {
     return(-1);
   }

   return(1);
}
