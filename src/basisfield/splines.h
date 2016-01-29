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
#ifndef __splines_h
#define __splines_h

#ifdef __cplusplus
   extern "C" {
#endif

/* Silly little macros. */

/*
Get index into C-type array given 3D 
subscript and array of given dimensions.
*/
#ifndef index
#define index(A,B,C,DIM) ((C)*DIM[0]*DIM[1] + (B)*DIM[0] + (A))
#endif

#ifndef MIN
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#endif

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

/* 
The following 2 macros are used to obtain start and
end indicies into (C-type 0 offset) image-matrix 
given knot spacing and spline-coefficient index.
*/

#ifndef k2strt
#define k2strt(K,KSP) (1+(KSP)*((K)-4)))
#endif

#ifndef k2end
#define k2end(K,KSP) ((K)*(KSP)-1)
#endif

/* 
The following 2 macros are used to go between 
knot-spacing and spline kernel size.
*/

#ifndef sz2ksp
#define sz2ksp(SZ) (((SZ)+1)/4)
#endif

#ifndef ksp2sz
#define ksp2sz(KSP) (4*(KSP)-1)
#endif

/* 
These are to ensure that right memory allocation
routines are used depending on wether the routines
are linked to mex-files or to C++.
*/

#ifdef MEX
#define my_calloc  mxCalloc
#define my_realloc mxRealloc
#define my_free    mxFree
#else
#define my_calloc  calloc
#define my_realloc realloc
#define my_free    free
#endif

void please_free(void *ptr);

int spline_kron(/* Input */
                int           ndim,     /* Dimensionality of spline. */
                int           dim[3],   /* Size of spline in the different directions. */
                double        *sp1d[3], /* Set of 1D splines. */
                /* Output */
                double        *spline); /* nD (n>0 & n<4) spline. */

int get_1D_spline(/* Input */
                  int     knsp,      /* Knot-spacing (in # of voxels) */
                  /* Output */
                  double  **spline); /* 1D spline function */

int get_1D_spline_d(/* Input */
                    int     knsp,       /* Knot-spacing (in # of voxels) */
                    /* Output */
                    double  **spline);  /* 1D spline function */

int get_1D_spline_dd(/* Input */
                     int     knsp,       /* Knot-spacing (in # of voxels) */
                     /* Output */
                     double  **spline);  /* 1D spline function */

int zoom_field(/* input */
               int      ndim,
               int      oksp[3],
               int      nksp[3],
               int      idim[3],
               double   *oc,
               /* Output */
	       double   *nc);

int zoom_field_by2(/* Input */
                   int     ndim,
                   int     ksp[3],
                   int     idim[3],
                   int     zdim,
                   double  *oc,
                   /* Output */
		   double  **nc);

int get_field(/* Input */
              int     ndim,
              int     cdim[3],
              double  *c,
              int     sdim[3],
              double  *sp,
              int     fdim[3],
              /* Output */
	      double  *f);

int get_range(/* Input */
              int      k,
              int      ssz,
              int      fsz,
              /* Output */
              int      *fs,
              int      *fe,
	      int      *ks);

int no_of_knots(int    ksp,
                int    msz);

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
	   double   *sp);

int make_Aty(/* Input. */
             int      ndim,      /* Actual dimensionality of problem (1, 2 or 3). */
             int      cdim[3],   /* # of knots in the three dimensions. */
             int      sdim[3],   /* Size of spline kernel in the three dimensions. */
             double   *spl,      /* Spline kernel. */
             int      idim[3],   /* Size of image matrix. */
             double   *ima,      /* Image. */
             double   *y,        /* y-vector. */
             /* Output. */
             double   *aty);      /* Resulting vector */

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
             double            *sp);      /* Array of non-zero values in sparse matric. */

int make_AtB(/* Input. */
             int      ndim,      /* Actual dimensionality of problem (1, 2 or 3). */
             int      cdim[3],   /* # of knots in the three dimensions. */
             int      sdim[3],   /* Size of spline kernel in the three dimensions. */
             double   *splB,     /* Spline kernel for B. */
             double   *splA,     /* Spline kernel for A. */
             int      idim[3],   /* Size of image matrix. */
             double   *imaB,     /* Image for B. */
             double   *imaA,     /* Image for A. */
             /* Output. */
             int      *irp,      /* Array of row-indicies. */
             int      *jcp,      /* Array of pointers into column-starts in irp. */
             double   *sp);      /* Array of non-zero values in sparse matric. */

int get_memen_grad(/* Input. */
                   int           ndim,    /* # of dimensions (1,2 or 3) */
                   const int     cdim[3], /* Size of coefficient array. */
                   const int     ksp[3],  /* Knot-spacings. */
                   const double  *beta,   /* Coefficients. */
                   /* Output.*/
                   double        *grad);

double get_memen(/* Input. */
                 int           ndim,    /* # of dimensions (1,2 or 3) */
                 const int     cdim[3], /* Size of coefficient array. */
                 const int     ksp[3],  /* Knot-spacings. */
                 const double  *beta);  /* Coefficients. */

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
                   double        **ovec);  /* Output vector */

int make_memen_H(/* Input. */
                 const int      ndim,      /* Actual dimensionality of problem (1, 2 or 3). */
                 int            cdim[3],   /* # of knots in the three dimensions. */
                 const int      ksp[3],    /* Knot spacing in the three dimensions. */
                 /* Output. */
                 int            *irp,      /* Array of row-indicies. */
                 int            *jcp,      /* Array of pointers into column-starts in irp. */
                 double         *sp);      /* Array of non-zero values in sparse matric. */

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
                int       ss1[3]);    /* Offset into first spline kernel (to handle edges). */

double dot_prod_H(int       i1,         /* [i1 j1 k1] is index of first spline kernel/coef. */ 
                  int       j1,
                  int       k1,
                  int       i2,         /* [i2 j2 k2] is index of second spline kernel/coef. */ 
                  int       j2,
                  int       k2,
                  int       sdim[3],    /* Size of spline kernel (common to s1 and s2). */
                  double    *s1,        /* First spline kernel. */
                  double    *s2);       /* Second spline kernel. */

double get_s_by_i(/* Input */
                  int      i,         /* [i j k] index of spline kernel/coef. */
                  int      j,
                  int      k,
                  int      sdim[3],   /* Size of spline kernel. */
                  int      idim[3],   /* Size of image. */
                  double   *spl,      /* Spline. */
                  double   *ima,      /* Image. */
                  /* Output */
                  double   *sbyi);    /* Spline multiplied with appuretenant values in image. */

int get_nabos(/* Input. */
              int   i,       /* Index (in one dimension) of spline-coef. */
              int   csz,     /* Total # of spline-coef. */
              int   ssz,     /* Size of spline kernel. */
              /* Output. */
              int   *ns,     /* Lowest index of overlapping neighbour. */
              int   *ne);    /* Highest index of overlapping neighbour. */

int get_A_nzmax(/* Input. */
                int      ndim,
                int      kdim[3],
                int      sdim[3],
	        int      idim[3]);

int get_AtA_nzmax(int          ndim,
                  const int    nknot[3],
                  const int    ksp[3]);

int n_nabo(int   i,
           int   n,
           int   ksp);

int AtranspA(/* Input. */
             int      *ir_in,
             int      *jc_in,
             double   *s_in,
             int      m,
             int      n,
             int      nzmax,
             /* Output. */
             int      **ir_out_orig, /* These have to be pointers */
             int      *jc_out, /* to pointers to allow for  */
             double   **s_out_orig); /* realloc.                  */

int AtranspB(/* Input. */
             int      *ir_inA,
             int      *jc_inA,
             double   *s_inA,
             int      mA,
             int      nA,
             int      *ir_inB,
             int      *jc_inB,
             double   *s_inB,
             int      mB,
             int      nB,
             int      nzmax,
             /* Output. */
             int      **ir_out_orig, /* These have to be pointers */
             int      *jc_out,       /* to pointers to allow for  */
             double   **s_out_orig); /* realloc.                  */

double find_val(int     *a,
                int     n,
                int     key,
                double  *val);

int cmpf(const void    *el1,
         const void    *el2);


int fnirt_zoom_field(/* input */
                     int      ndim,
                     int      oksp[3],
                     int      ocdim[3],
                     int      nksp[3],
                     int      idim[3],
                     double   *oc,
                     /* Output */
                     double   *nc);

int fnirt_zoom_field_by2(/* Input */
                         int     ndim,
                         int     ksp[3],
                         int     ocdim[3],
                         int     idim[3],
                         int     zdim,
                         double  *oc,
                         /* Output */
                         double  **nc);

#ifdef __cplusplus
   }
#endif

#endif
