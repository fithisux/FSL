/* GDCHART 0.11.0dev  GDCHART.C  28 Feb 2002 */
/* Copyright Bruce Verderaime 1998-2002 */

/* vi:set tabstop=4 */

#define GDC_INCL
#define GDC_LIB
#include "gdc.h"	/* gdc.h before system includes to pick up features */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
/*#include <string.h>*/
#include <stdarg.h>
#include <errno.h>

#include "gdchart.h"

#define HIGHSET		0
#define LOWSET		1
#define CLOSESET	2

/* scaled translation onto graph */
#define PX( x )		(int)( xorig + (setno*xdepth_3D) + (x)*xscl )
#define PY( y )		(int)( yorig - (setno*ydepth_3D) + (y)*yscl )
/*#define PV( y )		(int)( vyorig - (setno*ydepth_3D) + (y)*vyscl ) CFB */
#define PV( y )		(int)( yorig - (setno*ydepth_3D) + (y)*yscl )

#define SET_RECT( gdp, x1, x2, y1, y2 )	gdp[0].x = gdp[3].x = x1,	\
										gdp[0].y = gdp[1].y = y1,	\
										gdp[1].x = gdp[2].x = x2,	\
										gdp[2].y = gdp[3].y = y2


#define SET_3D_POLY( gdp, x1, x2, y1, y2, xoff, yoff )						\
								gdp[0].x  = x1,        gdp[0].y = y1,		\
								gdp[1].x  = x1+(xoff), gdp[1].y = y1-yoff,	\
								gdp[2].x  = x2+(xoff), gdp[2].y = y2-yoff,	\
								gdp[3].x  = x2,        gdp[3].y = y2
/* ------------------------------------------------------------------------- */
/* vals in pixels */
/* ref is front plane */
/* allows for intersecting 3D lines      */
/*  (also used for single 3D lines >:-Q  */
struct YS { int y1; int y2; float slope; int lnclr; int shclr; };
static int qcmpr( const void *a, const void *b )
{ if( ((struct YS*)a)->y2 < ((struct YS*)b)->y2 ) return 1;
 if( ((struct YS*)a)->y2 > ((struct YS*)b)->y2 ) return -1;
 return 0; }
void
draw_3d_line( gdImagePtr	im,
	      int			y0,
	      int			x1,
	      int			x2,
	      int			y1[],
	      int			y2[],
	      int			xdepth,
	      int			ydepth,
	      int			num_sets,
	      int			clr[],
	      int			clrshd[] )
{
#define F(x,i)	(int)( (float)((x)-x1)*slope[i]+(float)y1[i] )
  float		depth_slope  = xdepth==0? FLT_MAX: (float)ydepth/(float)xdepth;
  CREATE_ARRAY1( slope, float,     num_sets );		/* float		slope[num_sets] */
  CREATE_ARRAY1( lnclr, int,       num_sets );		/* int			slope[num_sets] */
  CREATE_ARRAY1( shclr, int,       num_sets );		/* int			slope[num_sets] */
  CREATE_ARRAY1( ypts,  struct YS, num_sets );		/* struct YS	slope[num_sets] */
  int			i;
  int			x;
  gdPoint		poly[4];

  for( i=0; i<num_sets; ++i )
    {
      /* lnclr[i] = clr[i]; */
      /* shclr[i] = clrshd[i]; */
      slope[i] = x2==x1? FLT_MAX: (float)(y2[i]-y1[i])/(float)(x2-x1);
    }

  for( x=x1+1; x<=x2; ++x )
    {
      for( i=0; i<num_sets; ++i )						/* load set of points */
	{
	  ypts[i].y1    = F(x-1,i);
	  ypts[i].y2    = F(x,i);
	  ypts[i].lnclr = clr[i];
	  ypts[i].shclr = clrshd[i];
	  ypts[i].slope = slope[i];
	}											/* sorted "lowest" first */
      qsort( ypts, num_sets, sizeof(struct YS), qcmpr );
      /* put out in that order */
      for( i=0; i<num_sets; ++i )
	{											/* top */
	  SET_3D_POLY( poly, x-1, x, ypts[i].y1, ypts[i].y2, xdepth, ydepth );
	  gdImageFilledPolygon( im, poly, 4,			/* depth_slope ever < 0 ? */
				-ypts[i].slope>depth_slope? ypts[i].shclr: ypts[i].lnclr );
	  if( x == x1+1 )								/* edging */
	    gdImageLine( im,
			 x-1, ypts[i].y2,
			 x-1+xdepth, ypts[i].y2-ydepth,
			 -ypts[i].slope<=depth_slope? ypts[i].shclr: ypts[i].lnclr );
	}
    }
  FREE_ARRAY1( slope );
  FREE_ARRAY1( lnclr );
  FREE_ARRAY1( shclr );
  FREE_ARRAY1( ypts );
}

/* ------------------------------------------------------------------------- */
/* vals in pixels */
/* ref is front plane */
void
draw_3d_area( gdImagePtr		im,
	      int				x1,
	      int				x2,
	      int				y0,			/* drawn from 0 */
	      int				y1,
	      int				y2,
	      int				xdepth,
	      int				ydepth,
	      int				clr,
	      int				clrshd )
{

  gdPoint     poly[4];
  int			y_intercept = 0;									/* if xdepth || ydepth */

  if( xdepth || ydepth )
    {
      float		line_slope   = x2==x1?    FLT_MAX: (float)-(y2-y1) / (float)(x2-x1);
      float		depth_slope  = xdepth==0? FLT_MAX: (float)ydepth/(float)xdepth;

      y_intercept = (y1 > y0 && y2 < y0) ||						/* line crosses y0 */
	(y1 < y0 && y2 > y0)?
	(int)((1.0/ABS(line_slope))*(float)(ABS(y1-y0)))+x1:
	0;										/* never */

      /* edging along y0 depth */
      gdImageLine( im, x1+xdepth, y0-ydepth, x2+xdepth, y0-ydepth, clrshd );

      SET_3D_POLY( poly, x1, x2, y1, y2, xdepth, ydepth );		/* top */
      gdImageFilledPolygon( im, poly, 4, line_slope>depth_slope? clrshd: clr );

      SET_3D_POLY( poly, x1, x2, y0, y0, xdepth, ydepth+1 );	/* along y axis */
      gdImageFilledPolygon( im, poly, 4, clr );

      SET_3D_POLY( poly, x2, x2, y0, y2, xdepth, ydepth );		/* side */
      gdImageFilledPolygon( im, poly, 4, clrshd );

      if( y_intercept )
	gdImageLine( im, y_intercept,        y0,
		     y_intercept+xdepth, y0-ydepth, clrshd );	/* edging */
      gdImageLine( im, x1, y0, x1+xdepth, y0-ydepth, clrshd );	/* edging */
      gdImageLine( im, x2, y0, x2+xdepth, y0-ydepth, clrshd );	/* edging */

      /* SET_3D_POLY( poly, x2, x2, y0, y2, xdepth, ydepth );	// side */
      /* gdImageFilledPolygon( im, poly, 4, clrshd ); */

      gdImageLine( im, x1, y1, x1+xdepth, y1-ydepth, clrshd );	/* edging */
      gdImageLine( im, x2, y2, x2+xdepth, y2-ydepth, clrshd );	/* edging */
    }

  if( y1 == y2 )												/* bar rect */
    SET_RECT( poly, x1, x2, y0, y1 );							/* front */
  else
    {
      poly[0].x = x1;	poly[0].y = y0;
      poly[1].x = x2;	poly[1].y = y0;
      poly[2].x = x2;	poly[2].y = y2;
      poly[3].x = x1;	poly[3].y = y1;
    }
  gdImageFilledPolygon( im, poly, 4, clr );

  gdImageLine( im, x1, y0, x2, y0, clrshd );			/* edging along y0 */

  if( (xdepth || ydepth) &&								/* front edging only on 3D */
      (y1<y0 || y2<y0) )									/* and only above y0 */
    {
      if( y1 > y0 && y2 < y0 )							/* line crosses from below y0 */
	gdImageLine( im, y_intercept, y0, x2, y2, clrshd );
      else
	if( y1 < y0 && y2 > y0 )							/* line crosses from above y0 */
	  gdImageLine( im, x1, y1, y_intercept, y0, clrshd );
	else												/* completely above */
	  gdImageLine( im, x1, y1, x2, y2, clrshd );
    }
}

/* ------------------------------------------------------------------------- */
/* vals in pixels */
/* ref is front plane */
void
draw_3d_bar( gdImagePtr			im,
	     int				x1,
	     int				x2,
	     int				y0,
	     int				yhigh,
	     int				xdepth,
	     int				ydepth,
	     int				clr,
	     int				clrshd )
{
#define SET_3D_BAR( gdp, x1, x2, y1, y2, xoff, yoff )						\
								gdp[0].x  = x1,        gdp[0].y = y1,		\
								gdp[1].x  = x1+(xoff), gdp[1].y = y1-yoff,	\
								gdp[2].x  = x2+(xoff), gdp[2].y = y2-yoff,	\
								gdp[3].x  = x2,        gdp[3].y = y2

  gdPoint     poly[4];
  int			usd = MIN( y0, yhigh );								/* up-side-down bars */


  if( xdepth || ydepth )
    {
      if( y0 != yhigh )											/* 0 height? */
	{
	  SET_3D_BAR( poly, x2, x2, y0, yhigh, xdepth, ydepth );	/* side */
	  gdImageFilledPolygon( im, poly, 4, clrshd );
	}

      SET_3D_BAR( poly, x1, x2, usd, usd, xdepth, ydepth );		/* top */
      gdImageFilledPolygon( im, poly, 4, clr );
    }

  SET_RECT( poly, x1, x2, y0, yhigh );							/* front */
  gdImageFilledPolygon( im, poly, 4, clr );

  if( xdepth || ydepth )
    gdImageLine( im, x1, usd, x2, usd, clrshd );
}

/* ------------------------------------------------------------------------- */
struct BS { float y1; float y2; int clr; int shclr; };
static int barcmpr( const void *a, const void *b )
{ if( ((struct BS*)a)->y2 < ((struct BS*)b)->y2 ) return -1;
 if( ((struct BS*)a)->y2 > ((struct BS*)b)->y2 ) return 1;
 return 0; }

/* ------------------------------------------------------------------------- */
/* simple two-point linear interpolation */
/* attempts between first, then nearest */
void
do_interpolations( int		num_points,
		   int		interp_point,
		   float	vals[] )
{
  int		i, j;
  float	v1 = GDC_NOVALUE,
    v2 = GDC_NOVALUE;
  int		p1 = -1,
    p2 = -1;

  /* find backwards */
  for( i=interp_point-1; i>=0 && p1==-1; --i )
    if( vals[i] != GDC_NOVALUE && vals[i] != GDC_INTERP_VALUE )
      {
	v1 = vals[i];
	p1 = i;
      }
  /* find forwards */
  for( j=interp_point+1; j<num_points && p2==-1; ++j )
    if( vals[j] != GDC_NOVALUE && vals[j] != GDC_INTERP_VALUE )
      {
	v2 = vals[j];
	p2 = j;
      }
  /* no forward sample, find backwards */
  for( ; i>=0 && p2==-1; --i )
    if( vals[i] != GDC_NOVALUE && vals[i] != GDC_INTERP_VALUE )
      {
	v2 = vals[i];
	p2 = i;
      }
  /* no backwards sample, find forwards */
  for( ; j<num_points && p1==-1; ++j )
    if( vals[j] != GDC_NOVALUE && vals[j] != GDC_INTERP_VALUE )
      {
	v1 = vals[j];
	p1 = j;
      }
  if( p1==-1 || p2==-1 ||							/* need both */
      p1 == p2 )									/* idiot */
    {
      vals[interp_point] = GDC_NOVALUE;
      return;
    }

  /* Point-slope formula */
  vals[interp_point] = ((v2-v1)/(float)(p2-p1)) * (float)(interp_point-p1) + v1;
  return;
}

/* ========================================================================= */
/* little error checking  0: ok,                      */
/*                     -ret: error no graph output    */
/*                      ret: error graph out          */
/* watch out for # params and array sizes==num_points */
/* ------------------------------------------------------------------------- */
/* original var arg interface */
int
out_graph( short		IMGWIDTH,		/* no check for a image that's too small to fit */
	   short		IMGHEIGHT,		/* needed info (labels, etc), could core dump */
	   FILE			*img_fptr,		/* open file pointer (img out) */
	   GDC_CHART_T	type,
	   int			num_points,     /* points along x axis (even iterval) */
	   /*	all arrays dependant on this */
	   char			*xlbl[],		/* array of xlabels */
	   int			num_sets,
	   ... )
{
  char	do_hlc = ( type == GDC_HILOCLOSE        ||
		   type == GDC_3DHILOCLOSE      ||
		   type == GDC_3DCOMBO_HLC_BAR  ||
		   type == GDC_3DCOMBO_HLC_AREA ||
		   type == GDC_COMBO_HLC_BAR    ||
		   type == GDC_COMBO_HLC_AREA );

  char	do_fb  = ( type == GDC_FLOATINGBAR ||
		   type == GDC_3DFLOATINGBAR );

  char	do_vol = ( type == GDC_COMBO_HLC_BAR   ||
		   type == GDC_COMBO_HLC_AREA  ||
		   type == GDC_COMBO_LINE_BAR  ||
		   type == GDC_COMBO_LINE_AREA ||
		   type == GDC_COMBO_LINE_LINE ||
		   type == GDC_3DCOMBO_HLC_BAR ||
		   type == GDC_3DCOMBO_HLC_AREA||
		   type == GDC_3DCOMBO_LINE_BAR||
		   type == GDC_3DCOMBO_LINE_AREA ||
		   type == GDC_3DCOMBO_LINE_LINE );

  int		num_arrays = num_sets * (do_hlc? 3:
					 do_fb?  2: 1);

  CREATE_ARRAY1( data, float, num_arrays*num_points );	/* float data[num_arrays*num_points]  */
  float	*combo_data = (float*)NULL;

  va_list	ap;
  int		i,
    rtn;

  va_start( ap, num_sets );
  for( i=0; i<num_arrays; ++i )
    memcpy( data+i*num_points, va_arg(ap, float*), num_points*sizeof(float) );
  if( do_vol )
    combo_data = va_arg(ap, float*);
  va_end(ap);
  rtn =  GDC_out_graph( IMGWIDTH,
			IMGHEIGHT,
			img_fptr,
			type,
			num_points,
			xlbl,
			num_sets,
			data,
			combo_data );
  FREE_ARRAY1( data );

  return rtn;
}

/* ------------------------------------------------------------------------- */
/* multi array interface */
int
GDC_out_graph( short		IMGWIDTH,		/* no check for a img that's too small to fit */
	       short		IMGHEIGHT,		/* needed info (labels, etc), could core dump */
	       FILE			*img_fptr,		/* open file pointer (img out) */
	       GDC_CHART_T	type,
	       int			num_points,     /* points along x axis (even iterval) */
	       /*	all arrays dependant on this */
	       char			*xlbl[],		/* array of xlabels */
	       int			num_sets,
	       float		*data,			/* (float*) cast on multi-dim array (num_sets > 1) */
	       float		*combo_data )	/* only used on COMBO chart types */
{
  int			i, j;
  /*int k; */

  int			graphwidth;
  int			grapheight;
  gdImagePtr	im;
  gdImagePtr	bg_img = NULL;

  float		xorig, yorig, vyorig = 0.0;
  float		yscl     = 0.0;
  float		vyscl    = 0.0;
  float		xscl     = 0.0;
  float		vhighest = -FLT_MAX;
  float		vlowest  = FLT_MAX;
  float		highest  = -FLT_MAX;
  float		lowest   = FLT_MAX;
  /*	gdPoint     volpoly[4];*/

  char		do_vol = ( type == GDC_COMBO_HLC_BAR   ||		/* aka: combo */
			   type == GDC_COMBO_HLC_AREA  ||
			   type == GDC_COMBO_LINE_BAR  ||
			   type == GDC_COMBO_LINE_AREA ||
			   type == GDC_COMBO_LINE_LINE ||
			   type == GDC_3DCOMBO_HLC_BAR ||
			   type == GDC_3DCOMBO_HLC_AREA||
			   type == GDC_3DCOMBO_LINE_BAR||
			   type == GDC_3DCOMBO_LINE_AREA ||
			   type == GDC_3DCOMBO_LINE_LINE );
  char		threeD = ( type == GDC_3DAREA          ||
			   type == GDC_3DLINE          ||
			   type == GDC_3DBAR           ||
			   type == GDC_3DFLOATINGBAR   ||
			   type == GDC_3DHILOCLOSE     ||
			   type == GDC_3DCOMBO_HLC_BAR ||
			   type == GDC_3DCOMBO_HLC_AREA||
			   type == GDC_3DCOMBO_LINE_BAR||
			   type == GDC_3DCOMBO_LINE_AREA ||
			   type == GDC_3DCOMBO_LINE_LINE );
  char		num_groups = num_sets;							/* set before num_sets gets adjusted */
  char		set_depth = ( GDC_stack_type == GDC_STACK_DEPTH )? num_groups:
    1;
  char		do_bar = ( type == GDC_3DBAR           ||		/* offset X objects to leave */
			   type == GDC_BAR             ||		/*  room at X(0) and X(n) */
			   type == GDC_3DFLOATINGBAR   ||		/*  i.e., not up against Y axes */
			   type == GDC_FLOATINGBAR);
  char		do_ylbl_fractions = 					      	/* %f format not given, or */
    ( !GDC_ylabel_fmt ||					/*  format doesn't have a %,g,e,E,f or F */
      strlen(GDC_ylabel_fmt) == strcspn(GDC_ylabel_fmt,"%geEfF") );
  float		ylbl_interval  = 0.0;
  int			xlbl_hgt       = 0;
  int			xdepth_3Dtotal = 0;
  int			ydepth_3Dtotal = 0;
  int			xdepth_3D      = 0;		/* affects PX() */
  int			ydepth_3D      = 0;		/* affects PY() and PV() */
  int			hlf_barwdth	   = 0;		/* half bar widths */
  int			hlf_hlccapwdth = 0;		/* half cap widths for HLC_I_CAP and DIAMOND */
  int			annote_len     = 0,
    annote_hgt     = 0;

  /* args */
  int			setno = 0;				/* affects PX() and PY() */
  CREATE_ARRAY1( uvals, float *, type == GDC_HILOCLOSE        ||
		 type == GDC_3DHILOCLOSE      ||
		 type == GDC_3DCOMBO_HLC_BAR  ||
		 type == GDC_3DCOMBO_HLC_AREA ||
		 type == GDC_COMBO_HLC_BAR    ||
		 type == GDC_COMBO_HLC_AREA?  num_sets *= 3:	/* 1 more last set is vol */
		 type == GDC_FLOATINGBAR      ||
		 type == GDC_3DFLOATINGBAR?   num_sets *= 2:
		 type == GDC_COMBO_LINE_BAR   ||
		 type == GDC_3DCOMBO_LINE_BAR ||
		 type == GDC_3DCOMBO_LINE_AREA||
		 type == GDC_3DCOMBO_LINE_LINE||
		 type == GDC_COMBO_LINE_AREA  ||
		 type == GDC_COMBO_LINE_LINE? num_sets:		/* 1 more last set is vol */
		 num_sets );
  CREATE_ARRAY1( ExtVolColor, int, num_points );				/* int	ExtVolColor[num_points],           */
  CREATE_ARRAY2( ExtColor,    int, num_sets, num_points );	/*		ExtColor[num_sets][num_points],    */
  CREATE_ARRAY2( ExtColorShd, int, threeD?num_sets:1,			/*		ExtColorShd[num_sets][num_points]; */
		 threeD?num_points:1 );		/* shade colors only with 3D */
  float		*uvol = 0;

  int			BGColor,
    LineColor,
    PlotColor,
    GridColor,
    VolColor = 0,
    /*ThumbDColor,*/
    /*    ThumbLblColor,*/
    /*ThumbUColor,*/
    /*				ArrowDColor, */
    /*				ArrowUColor, */
    AnnoteColor = 0;
#ifdef HAVE_LIBFREETYPE
  char		*gdc_title_font	      = GDC_title_font;			/* for convienience  */
  char		*gdc_ytitle_font      = GDC_ytitle_font;		/* in func calls */
  char		*gdc_xtitle_font      = GDC_xtitle_font;
  /*	char		*gdc_yaxis_font	      = GDC_yaxis_font; */
  char		*gdc_xaxis_font	      = GDC_xaxis_font;
  double		gdc_title_ptsize      = GDC_title_ptsize;
  double		gdc_ytitle_ptsize     = GDC_ytitle_ptsize;
  double		gdc_xtitle_ptsize     = GDC_xtitle_ptsize;
  /*	double		gdc_yaxis_ptsize      = GDC_yaxis_ptsize; */
  double		gdc_xaxis_ptsize      = GDC_xaxis_ptsize;
  double		gdc_xaxis_rad         = TO_RAD( GDC_xaxis_angle );
  char		*gdc_annotation_font  = GDC_annotation_font;
  double		gdc_annotation_ptsize = GDC_annotation_ptsize;

#else
  char		*gdc_title_font		 = NULL;
  char		*gdc_ytitle_font	 = NULL;
  char		*gdc_xtitle_font	 = NULL;
  /*	char		*gdc_yaxis_font		 = NULL; */
  char		*gdc_xaxis_font		 = NULL;
  double		gdc_title_ptsize	 = 0.0;
  double		gdc_ytitle_ptsize	 = 0.0;
  double		gdc_xtitle_ptsize	 = 0.0;
  /*	double		gdc_yaxis_ptsize	 = 0.0; */
  double		gdc_xaxis_ptsize	 = 0.0;
  double		gdc_xaxis_rad        = GDC_xaxis_angle==90.0? M_PI/2.0: 0.0;
  char		*gdc_annotation_font = NULL;
  double		gdc_annotation_ptsize=0.0;
#endif
  double		sin_xangle = 1.0,								/* calc only when&if needed */
    cos_xangle = 0.0;
  /* idiot checks */
  if( IMGWIDTH<=0 || IMGHEIGHT<=0 || (!img_fptr && GDC_generate_img) )
    {
      FREE_ARRAY1( uvals );
      FREE_ARRAY1( ExtVolColor );
      FREE_ARRAY2( ExtColor );
      FREE_ARRAY2( ExtColorShd );
      return -1;
    }
  if( num_points <= 0 )
    {
      out_err( IMGWIDTH, IMGHEIGHT, img_fptr, GDC_BGColor, GDC_LineColor, "No Data Available" );
      FREE_ARRAY1( uvals );
      FREE_ARRAY1( ExtVolColor );
      FREE_ARRAY2( ExtColor );
      FREE_ARRAY2( ExtColorShd );
      return 1;
    }

  load_font_conversions();
  if( GDC_thumbnail )
    {
      GDC_grid = FALSE;
      GDC_xaxis = FALSE;
      GDC_yaxis = FALSE;
    }

  /* ----- get args ----- */
  for( i=0; i<num_sets; ++i )
    uvals[i] = data+i*num_points;
  if( do_vol ){
    if( !combo_data )
      {
	out_err( IMGWIDTH, IMGHEIGHT, img_fptr, GDC_BGColor, GDC_LineColor, "No Combo Data Available" );
	FREE_ARRAY1( uvals );
	FREE_ARRAY1( ExtVolColor );
	FREE_ARRAY2( ExtColor );
	FREE_ARRAY2( ExtColorShd );
	return -2;
      }
    else
      uvol = combo_data;
  }
  /* ----- calculate interpretations first ----- */
  if( GDC_interpolations )
    {
      for( i=0; i<num_sets; ++i )
	for( j=0; j<num_points; ++j )
	  if( uvals[i][j] == GDC_INTERP_VALUE )
	    {
	      do_interpolations( num_points, j, uvals[i] );
	    }
      if( do_vol )
	for( j=0; j<num_points; ++j )
	  if( uvol[j] == GDC_INTERP_VALUE )
	    {
	      do_interpolations( num_points, j, uvol );
	    }
    }

  /* ----- highest & lowest values ----- */
  if( GDC_stack_type == GDC_STACK_SUM ) 		/* need to walk sideways */
    for( j=0; j<num_points; ++j )
      {
	float	set_sum = 0.0;
	for( i=0; i<num_sets; ++i )
	  if( uvals[i][j] != GDC_NOVALUE )
	    {
	      set_sum += uvals[i][j];
	      highest = MAX( highest, set_sum );
	      lowest  = MIN( lowest,  set_sum );
	    }
      }
  else
    if( GDC_stack_type == GDC_STACK_LAYER )		/* need to walk sideways */
      for( j=0; j<num_points; ++j )
	{
	  float	neg_set_sum = 0.0,
	    pos_set_sum = 0.0;
	  for( i=0; i<num_sets; ++i )
	    if( uvals[i][j] != GDC_NOVALUE ){
	      if( uvals[i][j] < 0.0 )
		neg_set_sum += uvals[i][j];
	      else
		pos_set_sum += uvals[i][j];
	    }
	  lowest  = MIN( lowest,  MIN(neg_set_sum,pos_set_sum) );
	  highest = MAX( highest, MAX(neg_set_sum,pos_set_sum) );
	}
    else
      for( i=0; i<num_sets; ++i )
	for( j=0; j<num_points; ++j )
	  if( uvals[i][j] != GDC_NOVALUE )
	    {
	      highest = MAX( uvals[i][j], highest );
	      lowest  = MIN( uvals[i][j], lowest );
	    }
  if( GDC_scatter )
    for( i=0; i<GDC_num_scatter_pts; ++i )
      {
	highest = MAX( (GDC_scatter+i)->val, highest );
	lowest  = MIN( (GDC_scatter+i)->val, lowest  );
      }
  if( do_vol )						  /* for now only one combo set allowed */
    {
      /* vhighest = 1.0; */
      /* vlowest  = 0.0; */
      for( j=0; j<num_points; ++j )
	if( uvol[j] != GDC_NOVALUE )
	  {
	    vhighest = MAX( uvol[j], vhighest );
	    vlowest  = MIN( uvol[j], vlowest );
	  }
      if( vhighest == -FLT_MAX )				/* no values */
	vhighest = 1.0;						/* for scaling, need a range */
      if( vlowest == FLT_MAX )
	vlowest = 0.0;
      if( type == GDC_COMBO_LINE_BAR    ||
	  type == GDC_COMBO_HLC_BAR     ||
	  type == GDC_COMBO_LINE_AREA   ||
	  type == GDC_COMBO_HLC_AREA    ||
	  type == GDC_3DCOMBO_LINE_BAR  ||
	  type == GDC_3DCOMBO_LINE_AREA ||
	  type == GDC_3DCOMBO_HLC_BAR   ||
	  type == GDC_3DCOMBO_HLC_AREA ){
	if( vhighest < 0.0 )
	  vhighest = 0.0;
	else
	  if( vlowest > 0.0 )
	    vlowest = 0.0;						/* bar, area should always start at 0 */
      }
      /*vhighest = highest;  CFB */
      /*vlowest = lowest;    CFB */
    }

  if( lowest == FLT_MAX )
    lowest = 0.0;
  if( highest == -FLT_MAX )
    highest = 1.0;							/* need a range */
  if( type == GDC_AREA  ||					/* bars and area should always start at 0 */
      type == GDC_BAR   ||
      type == GDC_3DBAR ||
      type == GDC_3DAREA ){
    if( highest < 0.0 )
      highest = 0.0;
    else
      if( lowest > 0.0 )						/* negs should be drawn from 0 */
	lowest = 0.0;
  }
  if( GDC_requested_ymin != GDC_NOVALUE && GDC_requested_ymin < lowest ){
    lowest = GDC_requested_ymin;
    vlowest = GDC_requested_ymin;
  }
  if( GDC_requested_ymax != GDC_NOVALUE && GDC_requested_ymax > highest ){
    highest = GDC_requested_ymax;
    vhighest = GDC_requested_ymax;
  }

  /* ----- graph height and width within the img height width ----- */
  /* grapheight/height is the actual size of the scalable graph */
  {
    int	title_hgt  = GDC_title? 2				/* title? horizontal text line(s) */
      + GDCfnt_sz(GDC_title,GDC_title_size,gdc_title_font,gdc_title_ptsize,0.0,NULL).h
      + 2:
      2;
    int	xlabel_hgt = 0;
    int	xtitle_hgt = GDC_xtitle? 1+GDCfnt_sz(GDC_xtitle,GDC_xtitle_size,gdc_xtitle_font,gdc_xtitle_ptsize,0.0,NULL).h+1: 0;
    int	ytitle_hgt = GDC_ytitle? 1+GDCfnt_sz(GDC_ytitle,GDC_ytitle_size,gdc_ytitle_font,gdc_ytitle_ptsize,M_PI/2.0,NULL).h+1: 0;
    int	vtitle_hgt = do_vol&&GDC_ytitle2? 1+GDCfnt_sz(GDC_ytitle2,GDC_ytitle_size,gdc_ytitle_font,gdc_ytitle_ptsize,M_PI/2.0,NULL).h+1: 0;
    int	ylabel_wth = 0;
    int	vlabel_wth = 0;

    int	xtics      = GDC_ticks && (GDC_grid||GDC_xaxis)? 1+2: 0;
    int	ytics      = GDC_ticks && (GDC_grid||GDC_yaxis)? 1+3: 0;
    int	vtics      = GDC_ticks && (GDC_yaxis&&do_vol)? 3+1: 0;


#define	HYP_DEPTH	( (double)((IMGWIDTH+IMGHEIGHT)/2) * ((double)GDC_3d_depth)/100.0 )
#define RAD_DEPTH	( (double)GDC_3d_angle*2*M_PI/360 )
    xdepth_3D      = threeD? (int)( cos(RAD_DEPTH) * HYP_DEPTH ): 0;
    ydepth_3D      = threeD? (int)( sin(RAD_DEPTH) * HYP_DEPTH ): 0;
    xdepth_3Dtotal = xdepth_3D*set_depth;
    ydepth_3Dtotal = ydepth_3D*set_depth;
    annote_hgt     = GDC_annotation && *(GDC_annotation->note)?
      1 +											/* space to note */
      (1+GDCfnt_sz( GDC_annotation->note,GDC_annotation_font_size,
		    gdc_annotation_font,gdc_annotation_ptsize,0.0,NULL ).h) +
      1 +											/* space under note */
      2: 0;										/* space to chart */
    annote_len     = GDC_annotation && *(GDC_annotation->note)?
      GDCfnt_sz( GDC_annotation->note,GDC_annotation_font_size,
		 gdc_annotation_font,gdc_annotation_ptsize,0.0,NULL ).w:
      0;

    /* find length of "longest" (Y) xaxis label */
    /* find the average "width" (X) xaxis label */
    /*	avg method fails when 2 or 3 very wide are consecutive, with the rest being thin */
    /*	this is most evident with horizontal (0deg) xlabels */
    /*	assume in this case they are quite uniform, e.g., dates */
    /* find affects on graphwidth/xorig of wildly overhanging angled labels */
    if( GDC_xaxis && xlbl )
      {
	int		biggest = -INT_MAX,
	  widest  = -INT_MAX;
#ifdef HAVE_LIBFREETYPE
	if( gdc_xaxis_rad!=M_PI/2.0 && gdc_xaxis_font && gdc_xaxis_ptsize )
	  {
	    sin_xangle = sin( gdc_xaxis_rad ),
	      cos_xangle = cos( gdc_xaxis_rad );
	  }
#endif

	for( i=0; i<num_points; ++i )
	  {
	    int	len = 0,
	      wdth = 0;
	    if( !GDC_xlabel_ctl ||
		(GDC_xlabel_ctl && GDC_xlabel_ctl[i]) )
	      {
		char	*sts;
		struct fnt_sz_t	lftsz = GDCfnt_sz( xlbl[i], GDC_xaxisfont_size,
						   gdc_xaxis_font, gdc_xaxis_ptsize, gdc_xaxis_rad,
						   &sts );

		if( gdc_xaxis_rad == M_PI/2.0 ||			/* no need to do the floating point math */
		    (sts && *sts) )							/* FT fail status, used default gdfont */
		  {
#ifdef DEBUG
		    fprintf( stderr, "TTF/FT failure: %s\n", sts );
#endif
		    len  = lftsz.w;
		    wdth = lftsz.h;
		  }
		else
		  if( gdc_xaxis_rad == 0.0 )			/* protect /0 */
		    {								/* reverse when horiz. */
		      len  = lftsz.h;
		      wdth = lftsz.w;
		    }
		  else									/* length & width due to the angle */
		    {
		      len  = (int)( (double)lftsz.w * sin_xangle + (double)lftsz.h * cos_xangle );
		      wdth = (int)( (double)lftsz.h / sin_xangle );
		    }
	      }
	    biggest = MAX( len, biggest );				/* last seg */
	    widest  = MAX( wdth, widest );				/* last seg */
	  }
	xlbl_hgt = 1+ widest +1;
	xlabel_hgt   = 1+ biggest +1;
      }

    grapheight = IMGHEIGHT - ( xtics          +
			       xtitle_hgt     +
			       xlabel_hgt     +
			       title_hgt      +
			       annote_hgt     +
			       ydepth_3Dtotal +
			       2 );
    if( GDC_hard_size && GDC_hard_grapheight )				/* user wants to use his */
      grapheight = GDC_hard_grapheight;
    GDC_hard_grapheight = grapheight;
    /* before width can be known... */
    /* ----- y labels intervals ----- */
    {
      float	tmp_highest;
      /* possible y gridline points */
#define	NUM_YPOINTS	(sizeof(ypoints_2f) / sizeof(float))
      float	ypoints_2f[] = { 1.0/1024,1.0/512,1.0/256.0,1.0/128.0,1.0/64.0, 1.0/32.0, 1.0/16.0, 1.0/8.0, 1.0/4.0, 1.0/3.0,1.0/2.0,
				 1.0,      2.0,      3.0, 4.0,     5.0, 6.0,    10.0,    25.0,
				 50.0,     100.0,    250.0,    500.0,   1000.0,  2500,    5000.0,
				 10000.0,  25000.0,  50000.0,  100000.0,500000.0,1000000, 5000000,
				 10000000, 50000000, 100000000 , 5e8,1e9,5e9 },
	ypoints_dec[NUM_YPOINTS] = 						/* "pretty" points for dec (non-fraction) this is used in normal (non-fraction) y labels*/
	  { 0.00001,0.0001,0.00025,0.0005, 0.001,0.0025,0.005,    0.01,     0.025,    0.05,     0.1,     0.2,     0.25,    0.5,
	  1.0,      2.0,      3.0,      5.0,     10.0,    25.0,
	  50.0,     100.0,    250.0,    500.0,   1000.0,  2500,    5000.0,
	    10000.0,  25000.0,  50000.0,  100000.0,500000.0,1000000, 5000000 , 1e7,2.5e7, 5e7, 1e8, 1e9 },
	  *ypoints = do_ylbl_fractions? ypoints_2f: ypoints_dec; 

	  int		max_num_ylbls;
	  int		longest_ylblen = 0;
	  /* maximum y lables that'll fit... */
	  max_num_ylbls = grapheight / (3+GDC_fontc[GDC_yaxisfont_size==GDC_TINY? GDC_yaxisfont_size+1:
						   GDC_yaxisfont_size].h);
	  if( max_num_ylbls < 3 )
	    {
	      /* gdImageDestroy(im);		haven't yet created it */
	      out_err( IMGWIDTH, IMGHEIGHT,
		       img_fptr,
		       GDC_BGColor, GDC_LineColor,
		       "Insificient Height" );
	      FREE_ARRAY1( uvals );
	      FREE_ARRAY1( ExtVolColor );
	      FREE_ARRAY2( ExtColor );
	      FREE_ARRAY2( ExtColorShd );
	      return 2;
	    }

	  {													/* one "space" interval above + below */
	    float	ylbl_density_space_intvl = ((float)max_num_ylbls-(1.0+1.0)) * (float)GDC_ylabel_density/100.0;
	    for( i=1; i<NUM_YPOINTS; ++i )
	      /* if( ypoints[i] > ylbl_interval ) */
	      /*	break; */
	      if( (highest-lowest)/ypoints[i] < ylbl_density_space_intvl )
		break;
	    /* gotta go through the above loop to catch the 'tweeners :-| */
	  }

	  ylbl_interval = GDC_requested_yinterval != GDC_NOVALUE &&
	    GDC_requested_yinterval > ypoints[i-1]?	  GDC_requested_yinterval:
	ypoints[i-1];

	  /* perform floating point remainders */
	  /* gonculate largest interval-point < lowest */
	  if( lowest != 0.0 &&
	      lowest != GDC_requested_ymin )
	    {
	      if( lowest < 0.0 )
		lowest -= ylbl_interval;
	      /* lowest = (lowest-ypoints[0]) - */
	      /* 			( ( ((lowest-ypoints[0])/ylbl_interval)*ylbl_interval ) - */
	      /* 			   ( (float)((int)((lowest-ypoints[0])/ylbl_interval))*ylbl_interval ) ); */
	      lowest = ylbl_interval * (float)(int)((lowest-ypoints[0])/ylbl_interval);
	    }
	  /* find smallest interval-point > highest */
	  tmp_highest = lowest;
	  do	/* while( (tmp_highest += ylbl_interval) <= highest ) */
	    {
	      int		nmrtr, dmntr, whole;
	      char	*price_to_str( float, int*, int*, int*, char* );
	      int		lbl_len;
	      char	foo[32];

	      if( GDC_yaxis )
		{											/* XPG2 compatibility */
		  sprintf( foo, do_ylbl_fractions? "%.0f": GDC_ylabel_fmt, tmp_highest );
		  lbl_len = ylbl_interval<1.0? strlen( price_to_str(tmp_highest,
								    &nmrtr,
								    &dmntr,
								    &whole,
								    do_ylbl_fractions? NULL: GDC_ylabel_fmt) ):
		    strlen( foo );
		  longest_ylblen = MAX( longest_ylblen, lbl_len );
		}
	    } while( (tmp_highest += ylbl_interval) <= highest );
	  ylabel_wth = longest_ylblen * GDC_fontc[GDC_yaxisfont_size].w;
	  highest = GDC_requested_ymax==GDC_NOVALUE? tmp_highest:
	MAX( GDC_requested_ymax, highest );

	  if( do_vol )
	    {
	      float	num_yintrvls = (highest-lowest) / ylbl_interval;
	      /*no skyscrapers  CFB*/
	      if( vhighest != 0.0 ) 
		vhighest += (vhighest-vlowest) / (num_yintrvls*2.0);
	      if( vlowest != 0.0 ) 
		vlowest -= (vhighest-vlowest) / (num_yintrvls*2.0); 
	      
	      if( GDC_yaxis2 )
		{
		  char	svlongest[32];
		  int		lbl_len_low  = sprintf( svlongest, GDC_ylabel2_fmt? GDC_ylabel2_fmt: "%.0f", vlowest );
		  int		lbl_len_high = sprintf( svlongest, GDC_ylabel2_fmt? GDC_ylabel2_fmt: "%.0f", vhighest );
		  vlabel_wth = 1
		    + MAX( lbl_len_low,lbl_len_high ) * GDC_fontc[GDC_yaxisfont_size].w;
		}
	    }
    }

    graphwidth = IMGWIDTH - ( ( (GDC_hard_size && GDC_hard_xorig)? GDC_hard_xorig:
				( ytitle_hgt +
				  ylabel_wth +
				  ytics ) )
			      + vtics
			      + vtitle_hgt
			      + vlabel_wth
			      + xdepth_3Dtotal );
    if( GDC_hard_size && GDC_hard_graphwidth )				/* user wants to use his */
      graphwidth = GDC_hard_graphwidth;
    GDC_hard_graphwidth = graphwidth;

    /* ----- scale to img size ----- */
    /* offset to 0 at lower left (where it should be) */
    xscl = (float)(graphwidth-xdepth_3Dtotal) / (float)(num_points + (do_bar?2:0));
    yscl = -((float)grapheight) / (float)(highest-lowest);
    if( do_vol )
      {
	float	hilow_diff = vhighest-vlowest==0.0? 1.0: vhighest-vlowest;

	vyscl = -((float)grapheight) / hilow_diff;
	vyorig = (float)grapheight
	  + ABS(vyscl) * MIN(vlowest,vhighest)
	  + ydepth_3Dtotal
	  + title_hgt
	  + annote_hgt;
      }
    xorig = (float)( IMGWIDTH - ( graphwidth +
				  vtitle_hgt +
				  vtics      +
				  vlabel_wth ) );
    if( GDC_hard_size && GDC_hard_xorig )
      xorig = GDC_hard_xorig;
    GDC_hard_xorig = xorig;
    /*	yorig = (float)grapheight + ABS(yscl * lowest) + ydepth_3Dtotal + title_hgt; */
    yorig = (float)grapheight
      + ABS(yscl) * MIN(lowest,highest)
      + ydepth_3Dtotal
      + title_hgt
      + annote_hgt;
    /*????	if( GDC_hard_size && GDC_hard_yorig )					vyorig too? */
    /*????		yorig = GDC_hard_yorig;	FRED - check email */
    GDC_hard_yorig = yorig;

    hlf_barwdth     = (int)( (float)(PX(2)-PX(1)) * (((float)GDC_bar_width/100.0)/2.0) );	/* used only for bars */
    hlf_hlccapwdth  = (int)( (float)(PX(2)-PX(1)) * (((float)GDC_HLC_cap_width/100.0)/2.0) );
  }
  /* scaled, sized, ready */


  /* ----- OK start the graphic ----- */
  if( (GDC_hold_img & GDC_REUSE_IMAGE) &&
      GDC_image != (void*)NULL )
    im = GDC_image;
  else
    im = gdImageCreate( IMGWIDTH, IMGHEIGHT );


  BGColor        = gdImageColorAllocate( im, l2gdcal(GDC_BGColor) );
  LineColor      = clrallocate( im, GDC_LineColor );
  PlotColor      = clrallocate( im, GDC_PlotColor );
  GridColor      = clrallocate( im, GDC_GridColor );
  if( do_vol )
    {
      VolColor     = clrallocate( im, GDC_VolColor );
      for( i=0; i<num_points; ++i )
	if( GDC_ExtVolColor )
	  ExtVolColor[i] = clrallocate( im, GDC_ExtVolColor[i] );
	else
	  ExtVolColor[i] = VolColor;
    }
  /*	ArrowDColor    = gdImageColorAllocate( im, 0xFF,    0, 0 ); */
  /*	ArrowUColor    = gdImageColorAllocate( im,    0, 0xFF, 0 ); */
  if( GDC_annotation )
    AnnoteColor = clrallocate( im, GDC_annotation->color );

  /* attempt to import optional background image */
  if( GDC_BGImage )
    {
      FILE	*in = fopen(GDC_BGImage, "rb");
      if( !in )
	{
	  ; /* Cant load background image, drop it */
	}
      else
	{
	  /* assume GIF */
	  /* should determine type by file extension, option, ... */
	  bg_img = gdImageCreateFromGif(in);
	  if( bg_img )					/* = */
	    {
	      int	bgxpos = gdImageSX(bg_img)<IMGWIDTH?  IMGWIDTH/2 - gdImageSX(bg_img)/2:  0,
		bgypos = gdImageSY(bg_img)<IMGHEIGHT? IMGHEIGHT/2 - gdImageSY(bg_img)/2: 0;


	      if( gdImageSX(bg_img) > IMGWIDTH ||				/* resize only if too big */
		  gdImageSY(bg_img) > IMGHEIGHT )				/*  [and center] */
		{
		  gdImageCopyResized( im, bg_img,				/* dst, src */
				      bgxpos, bgypos,			/* dstX, dstY */
				      0, 0,					/* srcX, srcY */
				      IMGWIDTH, IMGHEIGHT,	/* dstW, dstH */
				      IMGWIDTH, IMGHEIGHT );	/* srcW, srcH */
		}
	      else											/* just center */
		gdImageCopy( im, bg_img,					/* dst, src */
			     bgxpos, bgypos,				/* dstX, dstY */
			     0, 0,							/* srcX, srcY */
			     IMGWIDTH, IMGHEIGHT );			/* W, H */
	    }
	  fclose(in);
	}
    }

  for( j=0; j<num_sets; ++j )
    for( i=0; i<num_points; ++i )
      if( GDC_ExtColor )
	{
	  unsigned long	ext_clr = *(GDC_ExtColor+num_points*j+i);

	  ExtColor[j][i]            = clrallocate( im, ext_clr );
	  if( threeD )
	    ExtColorShd[j][i]     = clrshdallocate( im, ext_clr );
	}
      else if( GDC_SetColor )
	{
	  int	set_clr = GDC_SetColor[j];
	  ExtColor[j][i]     = clrallocate( im, set_clr );
	  if( threeD )
	    ExtColorShd[j][i] = clrshdallocate( im, set_clr );
	}
      else
	{
	  ExtColor[j][i]     = PlotColor;
	  if( threeD )
	    ExtColorShd[j][i] = clrshdallocate( im, GDC_PlotColor );
	}
			

  if( GDC_transparent_bg )
    gdImageColorTransparent( im, BGColor );

  if( GDC_title )
    {
      struct fnt_sz_t	tftsz      = GDCfnt_sz( GDC_title, GDC_title_size, gdc_title_font, gdc_title_ptsize, 0.0, NULL );
      int				titlecolor = clrallocate( im, GDC_TitleColor );

      GDCImageStringNL( im,
			&GDC_fontc[GDC_title_size],
			gdc_title_font, gdc_title_ptsize,
			0.0,
			IMGWIDTH/2 - tftsz.w/2,
			1,
			GDC_title,
			titlecolor,
			GDC_JUSTIFY_CENTER,
			NULL );
    }
  if( GDC_xtitle )
    {
      struct fnt_sz_t	xtftsz     = GDCfnt_sz( GDC_xtitle, GDC_xtitle_size, gdc_xtitle_font, gdc_xtitle_ptsize, 0.0, NULL );
      int				titlecolor = GDC_XTitleColor==GDC_DFLTCOLOR? PlotColor:
	clrallocate( im, GDC_XTitleColor );
      GDCImageStringNL( im,
			&GDC_fontc[GDC_xtitle_size],
			gdc_xtitle_font, gdc_xtitle_ptsize,
			0.0,
			IMGWIDTH/2 - xtftsz.w/2,
			IMGHEIGHT-1-xtftsz.h-1,
			GDC_xtitle,
			titlecolor,
			GDC_JUSTIFY_CENTER,
			NULL );
    }


  /* ----- start drawing ----- */
  /* ----- backmost first - border, grid & labels ----- */
  /* if no grid, on 3D, border needs to handle it */
  if( !GDC_grid && threeD &&
      ((GDC_border == GDC_BORDER_ALL) || (GDC_border & GDC_BORDER_X) || (GDC_border & GDC_BORDER_Y)) )
    {
      int	x1, x2,
	y1, y2;

      x1 = PX(0);
      y1 = PY(lowest);

      setno = set_depth;
      x2 = PX(0);
      y2 = PY(lowest);

      gdImageLine( im, x1, y1, x2, y2, LineColor );			/* depth at origin */
      if( (GDC_border == GDC_BORDER_ALL) || (GDC_border & GDC_BORDER_X) )
	gdImageLine( im, x2, y2, PX(num_points-1+(do_bar?2:0)), y2, LineColor );
      if( (GDC_border == GDC_BORDER_ALL) || (GDC_border & GDC_BORDER_Y) )
	gdImageLine( im, x2, PY(highest), x2, y2, LineColor );
      setno = 0;
    }
  if( GDC_grid || GDC_ticks || GDC_yaxis )
    {	/* grid lines & y label(s) */
      float	tmp_y = lowest;
      int		labelcolor = GDC_YLabelColor==GDC_DFLTCOLOR? 
	LineColor: clrallocate( im, GDC_YLabelColor );
      int		label2color = GDC_YLabel2Color==GDC_DFLTCOLOR? 
	VolColor: clrallocate( im, GDC_YLabel2Color );

      /* step from lowest to highest puting in labels and grid at interval points */
      /* since now "odd" intervals may be requested, try to step starting at 0,   */
      /* if lowest < 0 < highest                                                  */
      for( i=-1; i<=1; i+=2 )		        /* -1, 1 */
	{
	  if( i == -1 ){	
	    if( lowest >= 0.0 )	/*	all pos plotting */
	      continue;
	    else
	      tmp_y = MIN( 0, highest );        /*	step down to lowest */
	  }
	  if( i == 1 ){	
	    if( highest <= 0.0 )    /*	all neg plotting */
	      continue;
	    else
	      tmp_y = MAX( 0, lowest );		/*	step up to highest */
	  }
	  
	  do	/* while( (tmp_y (+-)= ylbl_interval) < [highest,lowest] ) */
	    {
	      int		n, d, w;
	      char	*price_to_str( float, int*, int*, int*, char* );
	      char	nmrtr[3+1], dmntr[3+1], whole[8];
	      char	all_whole = ylbl_interval<1.0? FALSE: TRUE;
	      
	      char	*ylbl_str = price_to_str( tmp_y,&n,&d,&w,
						  do_ylbl_fractions? NULL: GDC_ylabel_fmt );
	      if( do_ylbl_fractions )
		{
		  sprintf( nmrtr, "%d", n );
		  sprintf( dmntr, "%d", d );
		  sprintf( whole, "%d", w );
		}
	      
	      if( GDC_grid || GDC_ticks )
		{
		  int	x1, x2, y1, y2;
		  /* int	gridline_clr = tmp_y == 0.0? LineColor: GridColor; */
		  /* tics */
		  x1 = PX(0);		y1 = PY(tmp_y);
		  if( GDC_ticks )
		    gdImageLine( im, x1-2, y1, x1, y1, GridColor );
		  if( GDC_grid )
		    {
		      setno = set_depth;
		      x2 = PX(0);		y2 = PY(tmp_y);						/* w/ new setno */
		      gdImageLine( im, x1, y1, x2, y2, GridColor );		/* depth for 3Ds */
		      gdImageLine( im, x2, y2, PX(num_points-1+(do_bar?2:0)), y2, GridColor );
		      setno = 0;											/* set back to foremost */
		    }
		}
	      if( GDC_yaxis ){
		if( do_ylbl_fractions ){
		  if( w || (!w && !n && !d) )
		    {
		      gdImageString( im,
				     GDC_fontc[GDC_yaxisfont_size].f,
				     PX(0)-2-strlen(whole)*GDC_fontc[GDC_yaxisfont_size].w
				     - ( (!all_whole)?
					 (strlen(nmrtr)*GDC_fontc[GDC_yaxisfont_size-1].w +
					  GDC_fontc[GDC_yaxisfont_size].w                 +
					  strlen(nmrtr)*GDC_fontc[GDC_yaxisfont_size-1].w) :
					 1 ),
				     PY(tmp_y)-GDC_fontc[GDC_yaxisfont_size].h/2,
				     (unsigned char*)whole,labelcolor );
		    }
		  if( n )
		    {
		      gdImageString( im,
				     GDC_fontc[GDC_yaxisfont_size-1].f,
				       PX(0)-2-strlen(nmrtr)*GDC_fontc[GDC_yaxisfont_size-1].w
				     -GDC_fontc[GDC_yaxisfont_size].w
				     -strlen(nmrtr)*GDC_fontc[GDC_yaxisfont_size-1].w + 1,
				     PY(tmp_y)-GDC_fontc[GDC_yaxisfont_size].h/2 + 1,
				     (unsigned char*)nmrtr,
				     labelcolor );
		      gdImageString( im,
				     GDC_fontc[GDC_yaxisfont_size].f,
				     PX(0)-2-GDC_fontc[GDC_yaxisfont_size].w
				     -strlen(nmrtr)*GDC_fontc[GDC_yaxisfont_size-1].w,
				     PY(tmp_y)-GDC_fontc[GDC_yaxisfont_size].h/2,
				     (unsigned char*)"/",
				     labelcolor );
		      gdImageString( im,
				     GDC_fontc[GDC_yaxisfont_size-1].f,
				       PX(0)-2-strlen(nmrtr)*GDC_fontc[GDC_yaxisfont_size-1].w - 2,
				     PY(tmp_y)-GDC_fontc[GDC_yaxisfont_size].h/2 + 3,
				       (unsigned char*)dmntr,
				     labelcolor );
		    }
		}
		else
		  gdImageString( im,
				 GDC_fontc[GDC_yaxisfont_size].f,
				 PX(0)-2-strlen(ylbl_str)*GDC_fontc[GDC_yaxisfont_size].w,
				 PY(tmp_y)-GDC_fontc[GDC_yaxisfont_size].h/2,
				 (unsigned char*)ylbl_str,
				 labelcolor );
		
	      }
	      if( do_vol && GDC_yaxis2 )
		{
		  char	vylbl[16];
		  /* opposite of PV(y) */
		  sprintf( vylbl,
			   GDC_ylabel2_fmt? GDC_ylabel2_fmt: "%.0f",
			   ((float)(PY(tmp_y)+(setno*ydepth_3D)-vyorig))/vyscl );
		  
		  setno = set_depth;
		  if( GDC_ticks )
		    gdImageLine( im, PX(num_points-1+(do_bar?2:0)), PY(tmp_y),
				 PX(num_points-1+(do_bar?2:0))+3, PY(tmp_y), GridColor );
		  if( atof(vylbl) == 0.0 )									/* rounding can cause -0 */
		    strcpy( vylbl, "0" );
		  gdImageString( im,
				 GDC_fontc[GDC_yaxisfont_size].f,
				 PX(num_points-1+(do_bar?2:0))+6,
				 PY(tmp_y)-GDC_fontc[GDC_yaxisfont_size].h/2,
				 (unsigned char*)vylbl,
				 label2color );
		  setno = 0;
		}
	    }
	  while( ((i>0) && ((tmp_y += ylbl_interval) < highest)) ||
		 ((i<0) && ((tmp_y -= ylbl_interval) > lowest)) );
	}
      
      /* catch last (bottom) grid line - specific to an "off" requested interval */
      if( GDC_grid && threeD )
	{
	  setno = set_depth;
	  gdImageLine( im, PX(0), PY(lowest), PX(num_points-1+(do_bar?2:0)), PY(lowest), GridColor );
	  setno = 0;											/* set back to foremost */
	}
      
      /* vy axis title */
      if( do_vol && GDC_ytitle2 )
	{
	  struct fnt_sz_t	ytftsz     = GDCfnt_sz( GDC_ytitle2, GDC_ytitle_size, gdc_ytitle_font, gdc_ytitle_ptsize, 0.0, NULL );
	  int				titlecolor = GDC_YTitle2Color==GDC_DFLTCOLOR? VolColor:
	    clrallocate( im, GDC_YTitle2Color );
	  GDCImageStringNL( im,
			    &GDC_fontc[GDC_ytitle_size],
			    gdc_ytitle_font, gdc_ytitle_ptsize,
			    M_PI/2.0,
			    IMGWIDTH-(1+ytftsz.h),
			    yorig/2+ytftsz.w/2,
			    GDC_ytitle2,
			    titlecolor,
			    GDC_JUSTIFY_CENTER,
			    NULL );
	}

      /* y axis title */
      if( GDC_yaxis && GDC_ytitle )
	{
	  struct fnt_sz_t	ytftsz     = GDCfnt_sz( GDC_ytitle, GDC_ytitle_size, gdc_ytitle_font, gdc_ytitle_ptsize, 0.0, NULL );
	  int				titlecolor = GDC_YTitleColor==GDC_DFLTCOLOR? PlotColor:
	    clrallocate( im, GDC_YTitleColor );
	  GDCImageStringNL( im,
			    &GDC_fontc[GDC_ytitle_size],
			    gdc_ytitle_font, gdc_ytitle_ptsize,
			    M_PI/2.0,
			    1,
			    yorig/2+ytftsz.w/2,
			    GDC_ytitle,
			    titlecolor,
			    GDC_JUSTIFY_CENTER,
			    NULL );
	}
    }

  /* interviening set grids */
  /*  0 < setno < num_sets   non-inclusive, they've already been covered */
  if( GDC_grid && threeD )
    {
      for( setno=set_depth - 1;
	   setno > 0;
	   --setno )
	{
	  gdImageLine( im, PX(0), PY(lowest), PX(0), PY(highest), GridColor );
	  gdImageLine( im, PX(0), PY(lowest), PX(num_points-1+(do_bar?2:0)), PY(lowest), GridColor );
	}
      setno = 0;
    }

  if( ( GDC_grid || GDC_0Shelf ) &&							/* line color grid at 0 */
      ( (lowest < 0.0 && highest > 0.0) ||
	( (lowest == 0.0 || highest == 0.0) && !(GDC_border&GDC_BORDER_X) ) ) )
    {
      int	x1, x2, y1, y2;
      /* tics */
      x1 = PX(0);		y1 = PY(0);
      if( GDC_ticks )
	gdImageLine( im, x1-2, y1, x1, y1, LineColor );
      setno = set_depth;
      x2 = PX(0);		y2 = PY(0);								/* w/ new setno */
      gdImageLine( im, x1, y1, x2, y2, LineColor );			/* depth for 3Ds */
      gdImageLine( im, x2, y2, PX(num_points-1+(do_bar?2:0)), y2, LineColor );
      setno = 0;												/* set back to foremost */
    }

  /* x ticks and xlables */
  if( GDC_grid || GDC_xaxis )
    {
      int		num_xlbls  = graphwidth / 						/* maximum x lables that'll fit */
	( (GDC_xlabel_spacing==SHRT_MAX?0:GDC_xlabel_spacing) + xlbl_hgt );
      int		labelcolor = GDC_XLabelColor==GDC_DFLTCOLOR? LineColor:
	clrallocate( im, GDC_XLabelColor );
      for( i=0; i<num_points+(do_bar?2:0); ++i )
	{
	  int	xi = do_bar? i-1: i;
	  int	x1, x2, y1, y2, yh;									/* ticks & grids */
#define DO_TICK(x,y)			if( GDC_ticks )													\
												gdImageLine( im, x, y, x,  y+2, GridColor );				\
											else
#define DO_GRID(x1,y1,x2,y2)	if( GDC_grid )													\
												{															\
												gdImageLine( im, x1, y1, x2,  y2, GridColor ); /* depth */	\
												gdImageLine( im, x2, y2, x2,  yh, GridColor );				\
												}															\
											else

	  x1 = PX(i);		y1 = PY(lowest);
	  setno = set_depth;
	  x2 = PX(i);		y2 = PY(lowest);	yh = PY(highest);
	  setno = 0;												/* reset to foremost */

	  if( i == 0 ){					/* catch 3D Y back corner */
	    DO_GRID(x1,y1,x2,y2);
	  }
	  /* labeled points */
	  if( (!GDC_xlabel_ctl && ( (i%(1+num_points/num_xlbls) == 0) ||	/* # x labels are regulated */
				    num_xlbls >= num_points       ||
				    GDC_xlabel_spacing == SHRT_MAX ))
	      ||
	      (GDC_xlabel_ctl && xi>=0 && *(GDC_xlabel_ctl+xi)) )
	    {
	      DO_TICK(x1,y1);										/* labeled points tick & grid */
	      DO_GRID(x1,y1,x2,y2);

	      if( !do_bar || (i>0 && xi<num_points) )
		if( GDC_xaxis && xlbl && xlbl[xi] && *(xlbl[xi]) )
		  {
		    char			*sts;
		    struct fnt_sz_t	lftsz = GDCfnt_sz( xlbl[xi], GDC_xaxisfont_size,
							   gdc_xaxis_font, gdc_xaxis_ptsize,
							   gdc_xaxis_rad,
							   &sts );
		    if( gdc_xaxis_rad == M_PI/2.0 ||
			(sts && *sts) )							/* FT fail status, used default gdfont */
		      {
#ifdef DEBUG
			fprintf( stderr, "TTF/FT failure: %s\n", sts );
#endif
			GDCImageStringNL( im,
					  &GDC_fontc[GDC_xaxisfont_size],
					  gdc_xaxis_font, gdc_xaxis_ptsize,
					  M_PI/2.0,
					  PX(i)-1 - (lftsz.h/2),
					  PY(lowest) + 2 + 1 + lftsz.w,
					  xlbl[xi],
					  labelcolor,
					  GDC_JUSTIFY_RIGHT,
					  NULL );
		      }
		    else
		      if( gdc_xaxis_rad == 0.0 )
			GDCImageStringNL( im,
					  &GDC_fontc[GDC_xaxisfont_size],
					  gdc_xaxis_font, gdc_xaxis_ptsize,
					  0.0,
					  PX(i)-1 - (lftsz.w/2),
					  PY(lowest) + 2 + 1,
					  xlbl[xi],
					  labelcolor,
					  GDC_JUSTIFY_CENTER,
					  NULL );
		      else
			GDCImageStringNL( im,
					  &GDC_fontc[GDC_xaxisfont_size],
					  gdc_xaxis_font, gdc_xaxis_ptsize,
					  gdc_xaxis_rad,
					  PX(i)-1 - (int)((double)lftsz.w*cos_xangle
							  + (double)lftsz.h*gdc_xaxis_rad/(M_PI/2.0)/2.0),
					  PY(lowest) + 2 + 1 + (int)((double)lftsz.w*sin_xangle),
					  xlbl[xi],
					  labelcolor,
					  GDC_JUSTIFY_RIGHT,
					  NULL );
		  }
	    }
	  /* every point, on-point */
	  if( i>0 )
	    {
	      if( GDC_grid == GDC_TICK_POINTS )				/* --- GRID --- */
		DO_GRID( x1, y1, x2, y2 );
	      else if( GDC_grid > GDC_TICK_NONE )
		{
		  int k;
		  int	xt;
		  int	prevx      = PX(i-1);
		  int	intrv_dist = (x1-prevx)/(GDC_grid+1);
		  DO_GRID( x1, y1, x2, y2 );
		  for( k=0,          xt=prevx + intrv_dist;
		       k<GDC_grid && xt<x1;
		       ++k,          xt += intrv_dist )
		    DO_GRID( xt, y1, xt+xdepth_3Dtotal, y2 );
		}

	      if( GDC_ticks == GDC_TICK_POINTS )				/* --- TICKS --- */
		DO_TICK(x1,y1);
	      else if( GDC_ticks > GDC_TICK_NONE )
		{
		  int k;
		  int	xt;
		  int	prevx=PX(i-1);
		  int	intrv_dist = (x1-prevx)/(GDC_ticks+1);
		  DO_TICK( x1, y1 );
		  for( k=0,           xt=prevx + intrv_dist;
		       k<GDC_ticks && xt<x1;
		       ++k,           xt += intrv_dist )
		    DO_TICK( xt, y1 );
		}
	    }
	}
    }

  /* ----- secondard data plotting (volume) ----- */
  /*  so that grid lines appear under vol */
  if( do_vol )
    {
      setno = set_depth;
      if( type == GDC_COMBO_HLC_BAR    ||
	  type == GDC_COMBO_LINE_BAR   ||
	  type == GDC_3DCOMBO_LINE_BAR ||
	  type == GDC_3DCOMBO_HLC_BAR )
	{
	  if( uvol[0] != GDC_NOVALUE )
	    draw_3d_bar( im, PX(0), PX(0)+hlf_barwdth,
			 PV(0), PV(uvol[0]),
			 0, 0,
			 ExtVolColor[0],
			 ExtVolColor[0] );
	  for( i=1; i<num_points-1; ++i )
	    if( uvol[i] != GDC_NOVALUE )
	      draw_3d_bar( im, PX(i)-hlf_barwdth, PX(i)+hlf_barwdth,
			   PV(0), PV(uvol[i]),
			   0, 0,
			   ExtVolColor[i],
			   ExtVolColor[i] );
	  if( uvol[i] != GDC_NOVALUE )
	    draw_3d_bar( im, PX(i)-hlf_barwdth, PX(i),
			 PV(0), PV(uvol[i]),
			 0, 0,
			 ExtVolColor[i],
			 ExtVolColor[i] );
	}
      else
	if( type == GDC_COMBO_HLC_AREA   ||
	    type == GDC_COMBO_LINE_AREA  ||
	    type == GDC_3DCOMBO_LINE_AREA||
	    type == GDC_3DCOMBO_HLC_AREA )
	  {
	    for( i=1; i<num_points; ++i )
	      if( uvol[i-1] != GDC_NOVALUE && uvol[i] != GDC_NOVALUE )
		draw_3d_area( im, PX(i-1), PX(i),
			      PV(0), PV(uvol[i-1]), PV(uvol[i]),
			      0, 0,
			      ExtVolColor[i],
			      ExtVolColor[i] );
	  }
	else
	  if( type == GDC_COMBO_LINE_LINE ||
	      type == GDC_3DCOMBO_LINE_LINE )
	    {
	      for( i=1; i<num_points; ++i )
		if( uvol[i-1] != GDC_NOVALUE && uvol[i] != GDC_NOVALUE )
		  gdImageLine( im, PX(i-1), PV(uvol[i-1]),
			       PX(i),   PV(uvol[i]),
			       ExtVolColor[i] );
	    }
      setno = 0;
    }		/* volume polys done */

  if( GDC_annotation && threeD )		/* back half of annotation line */
    {
      int	x1 = PX(GDC_annotation->point+(do_bar?1:0)),
	y1 = PY(lowest);
      setno = set_depth;
      gdImageLine( im, x1, y1, PX(GDC_annotation->point+(do_bar?1:0)), PY(lowest), AnnoteColor );
      gdImageLine( im, PX(GDC_annotation->point+(do_bar?1:0)), PY(lowest),
		   PX(GDC_annotation->point+(do_bar?1:0)), PY(highest)-2, AnnoteColor );
      setno = 0;
    }

  /* ---------- start plotting the data ---------- */
  switch( type )
    {
    case GDC_3DBAR:					/* depth, width, y interval need to allow for whitespace between bars */
    case GDC_BAR:
      /* --------- */
      switch( GDC_stack_type )
	{
	case GDC_STACK_DEPTH:
	  for( setno=num_sets-1; setno>=0; --setno )		/* back sets first   PX, PY depth */
	    for( i=0; i<num_points; ++i )
	      if( uvals[setno][i] != GDC_NOVALUE )
		draw_3d_bar( im, PX(i+(do_bar?1:0))-hlf_barwdth, PX(i+(do_bar?1:0))+hlf_barwdth,
			     PY(0), PY(uvals[setno][i]),
			     xdepth_3D, ydepth_3D,
			     ExtColor[setno][i],
			     threeD? ExtColorShd[setno][i]: ExtColor[setno][i] );
	  setno = 0;
	  break;

	case GDC_STACK_LAYER:
	  {
	    int			j = 0;
	    CREATE_ARRAY1( barset, struct BS, num_sets );

	    /*				float	lasty[ num_points ]; */
	    /*				for( i=0; i<num_points; ++i ) */
	    /*					if( uvals[j][i] != GDC_NOVALUE ) */
	    /*						{ */
	    /*						lasty[i] = uvals[j][i]; */
	    /*						draw_3d_bar( im, PX(i+(do_bar?1:0))-hlf_barwdth, PX(i+(do_bar?1:0))+hlf_barwdth, */
	    /*										 PY(0), PY(uvals[j][i]), */
	    /*										 xdepth_3D, ydepth_3D, */
	    /*										 ExtColor[j][i], */
	    /*										 threeD? ExtColorShd[j][i]: ExtColor[j][i] ); */
	    /*						} */
	    for( i=0; i<num_points; ++i )
	      {
		float		lasty_pos = 0.0;
		float		lasty_neg = 0.0;
		int			k;

		for( j=0, k=0; j<num_sets; ++j )
		  {
		    if( uvals[j][i] != GDC_NOVALUE )
		      {
			if( uvals[j][i] < 0.0 )
			  {
			    barset[k].y1 = lasty_neg;
			    barset[k].y2 = uvals[j][i] + lasty_neg;
			    lasty_neg    = barset[k].y2;
			  }
			else
			  {
			    barset[k].y1 = lasty_pos;
			    barset[k].y2 = uvals[j][i] + lasty_pos;
			    lasty_pos    = barset[k].y2;
			  }
			barset[k].clr   = ExtColor[j][i];
			barset[k].shclr = threeD? ExtColorShd[j][i]: ExtColor[j][i];
			++k;
		      }
		  }
		qsort( barset, k, sizeof(struct BS), barcmpr );

		for( j=0; j<k; ++j )
		  {
		    draw_3d_bar( im,
				 PX(i+(do_bar?1:0))-hlf_barwdth, PX(i+(do_bar?1:0))+hlf_barwdth,
				 PY(barset[j].y1), PY(barset[j].y2),
				 xdepth_3D, ydepth_3D,
				 barset[j].clr,
				 barset[j].shclr );
		  }
	      }
	    FREE_ARRAY1( barset );
	  }
	  break;

	case GDC_STACK_SUM : break;
	case GDC_STACK_BESIDE:
	  {												/* h/.5, h/1, h/1.5, h/2, ... */
	    int	new_barwdth = (int)( (float)hlf_barwdth / ((float)num_sets/2.0) );
	    for( i=0; i<num_points; ++i )
	      for( j=0; j<num_sets; ++j )
		if( uvals[j][i] != GDC_NOVALUE )
		  draw_3d_bar( im, PX(i+(do_bar?1:0))-hlf_barwdth+new_barwdth*j+1,
			       PX(i+(do_bar?1:0))-hlf_barwdth+new_barwdth*(j+1),
			       PY(0), PY(uvals[j][i]),
			       xdepth_3D, ydepth_3D,
			       ExtColor[j][i],
			       threeD? ExtColorShd[j][i]: ExtColor[j][i] );
	  }
	  break;
	}
      break;

    case GDC_3DFLOATINGBAR:
    case GDC_FLOATINGBAR:
      /* --------- */
      switch( GDC_stack_type )
	{
	case GDC_STACK_DEPTH:
	  for( setno=num_groups-1; setno>=0; --setno )	/* back sets first   PX, PY depth */
	    for( i=0; i<num_points; ++i )
	      if( uvals[0+setno*2][i] != GDC_NOVALUE &&
		  uvals[1+setno*2][i] != GDC_NOVALUE &&
		  uvals[1+setno*2][i] > uvals[0+setno*2][i] )
		draw_3d_bar( im, PX(i+(do_bar?1:0))-hlf_barwdth, PX(i+(do_bar?1:0))+hlf_barwdth,
			     PY(uvals[0+setno*2][i]), PY(uvals[1+setno*2][i]),
			     xdepth_3D, ydepth_3D,
			     ExtColor[setno][i],
			     threeD? ExtColorShd[setno][i]: ExtColor[setno][i] );
	  setno = 0;
	  break;
	case GDC_STACK_SUM: break;
	case GDC_STACK_LAYER: break;

	case GDC_STACK_BESIDE:
	  {												/* h/.5, h/1, h/1.5, h/2, ... */
	    int	new_barwdth = (int)( (float)hlf_barwdth / ((float)num_groups/2.0) );
	    for( i=0; i<num_points; ++i )
	      for( j=0; j<num_groups; ++j )
		if( uvals[0+j*2][i] != GDC_NOVALUE &&
		    uvals[1+j*2][i] != GDC_NOVALUE &&
		    uvals[1+j*2][i] > uvals[0+j*2][i] )
		  draw_3d_bar( im, PX(i+(do_bar?1:0))-hlf_barwdth+new_barwdth*j+1,
			       PX(i+(do_bar?1:0))-hlf_barwdth+new_barwdth*(j+1),
			       PY(uvals[0+j*2][i]), PY(uvals[1+j*2][i]),
			       xdepth_3D, ydepth_3D,
			       ExtColor[j][i],
			       threeD? ExtColorShd[j][i]: ExtColor[j][i] );
	  }
	  break;
	}
      break;

    case GDC_LINE:
    case GDC_COMBO_LINE_BAR:
    case GDC_COMBO_LINE_AREA:
    case GDC_COMBO_LINE_LINE:
      for( j=num_sets-1; j>=0; --j )
	for( i=1; i<num_points; ++i )
	  if( uvals[j][i-1] != GDC_NOVALUE && uvals[j][i] != GDC_NOVALUE )
	    {
	      gdImageLine( im, PX(i-1), PY(uvals[j][i-1]), PX(i), PY(uvals[j][i]), ExtColor[j][i] );
	      gdImageLine( im, PX(i-1), PY(uvals[j][i-1])+1, PX(i), PY(uvals[j][i])+1, ExtColor[j][i] );
	    }
	  else
	    {
	      if( uvals[j][i-1] != GDC_NOVALUE )
		gdImageSetPixel( im, PX(i-1), PY(uvals[j][i-1]), ExtColor[j][i] );
	      if( uvals[j][i] != GDC_NOVALUE )
		gdImageSetPixel( im, PX(i), PY(uvals[j][i]), ExtColor[j][i] );
	    }
      break;

    case GDC_3DLINE:
    case GDC_3DCOMBO_LINE_BAR:
    case GDC_3DCOMBO_LINE_AREA:
    case GDC_3DCOMBO_LINE_LINE:
      {
	CREATE_ARRAY1( y1, int, num_sets );
	CREATE_ARRAY1( y2, int, num_sets );
	CREATE_ARRAY1( clr,    int, num_sets );
	CREATE_ARRAY1( clrshd, int, num_sets );

	for( i=1; i<num_points; ++i )
	  {
	    if( GDC_stack_type == GDC_STACK_DEPTH )
	      {
		for( j=num_sets-1; j>=0; --j )
		  if( uvals[j][i-1] != GDC_NOVALUE &&
		      uvals[j][i]   != GDC_NOVALUE )
		    {
		      setno = j;
		      y1[j] = PY(uvals[j][i-1]);
		      y2[j] = PY(uvals[j][i]);

		      draw_3d_line( im,
				    PY(0),
				    PX(i-1), PX(i), 
				    &(y1[j]), &(y2[j]),
				    xdepth_3D, ydepth_3D,
				    1,
				    &(ExtColor[j][i]),
				    &(ExtColorShd[j][i]) );
		      setno = 0;
		    }
	      }
	    else
	      if( GDC_stack_type == GDC_STACK_BESIDE ||
		  GDC_stack_type == GDC_STACK_SUM )			/* all same plane */
		{
		  int		set;
		  float	usey1 = 0.0,
		    usey2 = 0.0;
		  for( j=0,set=0; j<num_sets; ++j )
		    if( uvals[j][i-1] != GDC_NOVALUE &&
			uvals[j][i]   != GDC_NOVALUE )
		      {
			if( GDC_stack_type == GDC_STACK_SUM )
			  {
			    usey1 += uvals[j][i-1];
			    usey2 += uvals[j][i];
			  }
			else
			  {
			    usey1 = uvals[j][i-1];
			    usey2 = uvals[j][i];
			  }
			y1[set]     = PY(usey1);
			y2[set]     = PY(usey2);
			clr[set]    = ExtColor[j][i];
			clrshd[set] = ExtColorShd[j][i];	/* fred */
			++set;
		      }
		  draw_3d_line( im,
				PY(0),
				PX(i-1), PX(i), 
				y1, y2,
				xdepth_3D, ydepth_3D,
				set,
				clr,
				clrshd );
		}
	  }
	FREE_ARRAY1( clr );
	FREE_ARRAY1( clrshd );
	FREE_ARRAY1( y1 );
	FREE_ARRAY1( y2 );
      }
      break;

    case GDC_AREA:
    case GDC_3DAREA:
      switch( GDC_stack_type )
	{
	case GDC_STACK_SUM:
	  {
	    CREATE_ARRAY1( lasty, float, num_points );
	    int		j = 0;
	    for( i=1; i<num_points; ++i )
	      if( uvals[j][i] != GDC_NOVALUE )
		{
		  lasty[i] = uvals[j][i];
		  if( uvals[j][i-1] != GDC_NOVALUE )
		    draw_3d_area( im, PX(i-1), PX(i),
				  PY(0), PY(uvals[j][i-1]), PY(uvals[j][i]),
				  xdepth_3D, ydepth_3D,
				  ExtColor[j][i],
				  threeD? ExtColorShd[j][i]: ExtColor[j][i] );
		}
	    for( j=1; j<num_sets; ++j )
	      for( i=1; i<num_points; ++i )
		if( uvals[j][i] != GDC_NOVALUE && uvals[j][i-1] != GDC_NOVALUE )
		  {
		    draw_3d_area( im, PX(i-1), PX(i),
				  PY(lasty[i]), PY(lasty[i-1]+uvals[j][i-1]), PY(lasty[i]+uvals[j][i]),
				  xdepth_3D, ydepth_3D,
				  ExtColor[j][i],
				  threeD? ExtColorShd[j][i]: ExtColor[j][i] );
		    lasty[i] += uvals[j][i];
		  }
	    FREE_ARRAY1( lasty );
	  }
	  break;

	case GDC_STACK_BESIDE:								/* behind w/o depth */
	  for( j=num_sets-1; j>=0; --j )					/* back sets 1st  (setno = 0) */
	    for( i=1; i<num_points; ++i )
	      if( uvals[j][i-1] != GDC_NOVALUE && uvals[j][i] != GDC_NOVALUE )
		draw_3d_area( im, PX(i-1), PX(i),
			      PY(0), PY(uvals[j][i-1]), PY(uvals[j][i]),
			      xdepth_3D, ydepth_3D,
			      ExtColor[j][i],
			      threeD? ExtColorShd[j][i]: ExtColor[j][i] );
	  break;

	case GDC_STACK_DEPTH:
	default:
	  for( setno=num_sets-1; setno>=0; --setno )		/* back sets first   PX, PY depth */
	    for( i=1; i<num_points; ++i )
	      if( uvals[setno][i-1] != GDC_NOVALUE && uvals[setno][i] != GDC_NOVALUE )
		draw_3d_area( im, PX(i-1), PX(i),
			      PY(0), PY(uvals[setno][i-1]), PY(uvals[setno][i]),
			      xdepth_3D, ydepth_3D,
			      ExtColor[setno][i],
			      threeD? ExtColorShd[setno][i]: ExtColor[setno][i] );
	  setno = 0;
	}
      break;

    case GDC_3DHILOCLOSE:
    case GDC_3DCOMBO_HLC_BAR:
    case GDC_3DCOMBO_HLC_AREA:
      {
	gdPoint     poly[4];
	for( j=num_groups-1; j>=0; --j )
	  {
	    for( i=1; i<num_points+1; ++i )
	      if( uvals[CLOSESET+j*3][i-1] != GDC_NOVALUE )
		{
		  if( (GDC_HLC_style & GDC_HLC_I_CAP) &&			/* bottom half of 'I' */
		      uvals[LOWSET+j*3][i-1] != GDC_NOVALUE )
		    {
		      SET_3D_POLY( poly, PX(i-1)-hlf_hlccapwdth, PX(i-1)+hlf_hlccapwdth,
				   PY(uvals[LOWSET+j*3][i-1]), PY(uvals[LOWSET+j*3][i-1]),
				   xdepth_3D, ydepth_3D );
		      gdImageFilledPolygon( im, poly, 4, ExtColor[LOWSET+j*3][i-1] );
		      gdImagePolygon( im, poly, 4, ExtColorShd[LOWSET+j*3][i-1] );
		    }
		  /* all HLC have vert line */
		  if( uvals[LOWSET+j*3][i-1] != GDC_NOVALUE )
		    {											/* bottom 'half' */
		      SET_3D_POLY( poly, PX(i-1), PX(i-1),
				   PY(uvals[LOWSET+j*3][i-1]), PY(uvals[CLOSESET+j*3][i-1]),
				   xdepth_3D, ydepth_3D );
		      gdImageFilledPolygon( im, poly, 4, ExtColor[LOWSET+j*3][i-1] );
		      gdImagePolygon( im, poly, 4, ExtColorShd[LOWSET+j*3][i-1] );
		    }
		  if( uvals[HIGHSET+j*3][i-1] != GDC_NOVALUE )
		    {											/* top 'half' */
		      SET_3D_POLY( poly, PX(i-1), PX(i-1),
				   PY(uvals[CLOSESET+j*3][i-1]), PY(uvals[HIGHSET+j*3][i-1]),
				   xdepth_3D, ydepth_3D );
		      gdImageFilledPolygon( im, poly, 4, ExtColor[HIGHSET+j*3][i-1] );
		      gdImagePolygon( im, poly, 4, ExtColorShd[HIGHSET+j*3][i-1] );
		    }
		  /* line at close */
		  gdImageLine( im, PX(i-1),           PY(uvals[CLOSESET+j*3][i-1]),
			       PX(i-1)+xdepth_3D, PY(uvals[CLOSESET+j*3][i-1])-ydepth_3D,
			       ExtColorShd[CLOSESET+j*3][i-1] );
		  /* top half 'I' */
		  if( !( (GDC_HLC_style & GDC_HLC_DIAMOND) &&
			 (PY(uvals[HIGHSET+j*3][i-1]) > PY(uvals[CLOSESET+j*3][i-1])-hlf_hlccapwdth) ) &&
		      uvals[HIGHSET+j*3][i-1] != GDC_NOVALUE )
		    if( GDC_HLC_style & GDC_HLC_I_CAP )
		      {
			SET_3D_POLY( poly, PX(i-1)-hlf_hlccapwdth, PX(i-1)+hlf_hlccapwdth,
				     PY(uvals[HIGHSET+j*3][i-1]), PY(uvals[HIGHSET+j*3][i-1]),
				     xdepth_3D, ydepth_3D );
			gdImageFilledPolygon( im, poly, 4, ExtColor[HIGHSET+j*3][i-1] );
			gdImagePolygon( im, poly, 4, ExtColorShd[HIGHSET+j*3][i-1] );
		      }

		  if( i < num_points &&
		      uvals[CLOSESET+j*3][i] != GDC_NOVALUE )
		    {
		      if( GDC_HLC_style & GDC_HLC_CLOSE_CONNECTED )	/* line from prev close */
			{
			  SET_3D_POLY( poly, PX(i-1), PX(i),
				       PY(uvals[CLOSESET+j*3][i-1]), PY(uvals[CLOSESET+j*3][i-1]),
				       xdepth_3D, ydepth_3D );
			  gdImageFilledPolygon( im, poly, 4, ExtColor[CLOSESET+j*3][i] );
			  gdImagePolygon( im, poly, 4, ExtColorShd[CLOSESET+j*3][i] );
			}
		      else	/* CLOSE_CONNECTED and CONNECTING are mutually exclusive */
			if( GDC_HLC_style & GDC_HLC_CONNECTING )	/* thin connecting line */
			  {
			    int	y1 = PY(uvals[CLOSESET+j*3][i-1]),
			      y2 = PY(uvals[CLOSESET+j*3][i]);
			    draw_3d_line( im,
					  PY(0),
					  PX(i-1), PX(i),
					  &y1, &y2,					/* rem only 1 set */
					  xdepth_3D, ydepth_3D,
					  1,
					  &(ExtColor[CLOSESET+j*3][i]),
					  &(ExtColorShd[CLOSESET+j*3][i]) );
			    /* edge font of it */
			    gdImageLine( im, PX(i-1), PY(uvals[CLOSESET+j*3][i-1]),
					 PX(i), PY(uvals[CLOSESET+j*3][i]),
					 ExtColorShd[CLOSESET+j*3][i] );
			  }
		      /* top half 'I' again */
		      if( PY(uvals[CLOSESET+j*3][i-1]) <= PY(uvals[CLOSESET+j*3][i]) &&
			  uvals[HIGHSET+j*3][i-1] != GDC_NOVALUE  )
			if( GDC_HLC_style & GDC_HLC_I_CAP )
			  {
			    SET_3D_POLY( poly, PX(i-1)-hlf_hlccapwdth, PX(i-1)+hlf_hlccapwdth,
					 PY(uvals[HIGHSET+j*3][i-1]), PY(uvals[HIGHSET+j*3][i-1]),
					 xdepth_3D, ydepth_3D );
			    gdImageFilledPolygon( im, poly, 4, ExtColor[HIGHSET+j*3][i-1] );
			    gdImagePolygon( im, poly, 4, ExtColorShd[HIGHSET+j*3][i-1] );
			  }
		    }
		  if( GDC_HLC_style & GDC_HLC_DIAMOND )
		    {									/* front */
		      poly[0].x = PX(i-1)-hlf_hlccapwdth;
		      poly[0].y = PY(uvals[CLOSESET+j*3][i-1]);
		      poly[1].x = PX(i-1);
		      poly[1].y = PY(uvals[CLOSESET+j*3][i-1])+hlf_hlccapwdth;
		      poly[2].x = PX(i-1)+hlf_hlccapwdth;
		      poly[2].y = PY(uvals[CLOSESET+j*3][i-1]);
		      poly[3].x = PX(i-1);
		      poly[3].y = PY(uvals[CLOSESET+j*3][i-1])-hlf_hlccapwdth;
		      gdImageFilledPolygon( im, poly, 4, ExtColor[CLOSESET+j*3][i-1] );
		      gdImagePolygon( im, poly, 4, ExtColorShd[CLOSESET+j*3][i-1] );
		      /* bottom side */
		      SET_3D_POLY( poly, PX(i-1), PX(i-1)+hlf_hlccapwdth,
				   PY(uvals[CLOSESET+j*3][i-1])+hlf_hlccapwdth,
				   PY(uvals[CLOSESET+j*3][i-1]),
				   xdepth_3D, ydepth_3D );
		      gdImageFilledPolygon( im, poly, 4, ExtColorShd[CLOSESET+j*3][i-1] );
		      /* gdImagePolygon( im, poly, 4, ExtColor[CLOSESET+j*3][i-1] ); */
		      /* top side */
		      SET_3D_POLY( poly, PX(i-1), PX(i-1)+hlf_hlccapwdth,
				   PY(uvals[CLOSESET+j*3][i-1])-hlf_hlccapwdth,
				   PY(uvals[CLOSESET+j*3][i-1]),
				   xdepth_3D, ydepth_3D );
		      gdImageFilledPolygon( im, poly, 4, ExtColor[CLOSESET+j*3][i-1] );
		      gdImagePolygon( im, poly, 4, ExtColorShd[CLOSESET+j*3][i-1] );
		    }
		}
	  }
      }
      break;

    case GDC_HILOCLOSE:
    case GDC_COMBO_HLC_BAR:
    case GDC_COMBO_HLC_AREA:
      for( j=num_groups-1; j>=0; --j )
	{
	  for( i=0; i<num_points; ++i )
	    if( uvals[CLOSESET+j*3][i] != GDC_NOVALUE )
	      {											/* all HLC have vert line */
		if( uvals[LOWSET+j*3][i] != GDC_NOVALUE )
		  gdImageLine( im, PX(i), PY(uvals[CLOSESET+j*3][i]),
			       PX(i), PY(uvals[LOWSET+j*3][i]),
			       ExtColor[LOWSET+(j*3)][i] );
		if( uvals[HIGHSET+j*3][i] != GDC_NOVALUE )
		  gdImageLine( im, PX(i), PY(uvals[HIGHSET+j*3][i]),
			       PX(i), PY(uvals[CLOSESET+j*3][i]),
			       ExtColor[HIGHSET+j*3][i] );

		if( GDC_HLC_style & GDC_HLC_I_CAP )
		  {
		    if( uvals[LOWSET+j*3][i] != GDC_NOVALUE )
		      gdImageLine( im, PX(i)-hlf_hlccapwdth, PY(uvals[LOWSET+j*3][i]),
				   PX(i)+hlf_hlccapwdth, PY(uvals[LOWSET+j*3][i]),
				   ExtColor[LOWSET+j*3][i] );
		    if( uvals[HIGHSET+j*3][i] != GDC_NOVALUE )
		      gdImageLine( im, PX(i)-hlf_hlccapwdth, PY(uvals[HIGHSET+j*3][i]),
				   PX(i)+hlf_hlccapwdth, PY(uvals[HIGHSET+j*3][i]),
				   ExtColor[HIGHSET+j*3][i] );
		  }
		if( GDC_HLC_style & GDC_HLC_DIAMOND )
		  {
		    gdPoint         cd[4];

		    cd[0].x = PX(i)-hlf_hlccapwdth;	cd[0].y = PY(uvals[CLOSESET+j*3][i]);
		    cd[1].x = PX(i);	cd[1].y = PY(uvals[CLOSESET+j*3][i])+hlf_hlccapwdth;
		    cd[2].x = PX(i)+hlf_hlccapwdth;	cd[2].y = PY(uvals[CLOSESET+j*3][i]);
		    cd[3].x = PX(i);	cd[3].y = PY(uvals[CLOSESET+j*3][i])-hlf_hlccapwdth;
		    gdImageFilledPolygon( im, cd, 4, ExtColor[CLOSESET+j*3][i] );
		  }
	      }
	  for( i=1; i<num_points; ++i )
	    if( uvals[CLOSESET+j*3][i-1] != GDC_NOVALUE && uvals[CLOSESET+j*3][i] != GDC_NOVALUE )
	      {
		if( GDC_HLC_style & GDC_HLC_CLOSE_CONNECTED )	/* line from prev close */
		  gdImageLine( im, PX(i-1), PY(uvals[CLOSESET+j*3][i-1]),
			       PX(i), PY(uvals[CLOSESET+j*3][i-1]),
			       ExtColor[CLOSESET+j*3][i] );
		else	/* CLOSE_CONNECTED and CONNECTING are mutually exclusive */
		  if( GDC_HLC_style & GDC_HLC_CONNECTING )		/* thin connecting line */
		    gdImageLine( im, PX(i-1), PY(uvals[CLOSESET+j*3][i-1]),
				 PX(i), PY(uvals[CLOSESET+j*3][i]),
				 ExtColor[CLOSESET+j*3][i] );
	      }
	}
      break;
    }
  setno = 0;

  /* ---------- scatter points  over all other plots ---------- */
  /* scatters, by their very nature, don't lend themselves to standard array of points */
  /* also, this affords the opportunity to include scatter points onto any type of chart */
  /* drawing of the scatter point should be an exposed function, so the user can */
  /*  use it to draw a legend, and/or add their own */
  if( GDC_scatter )
    {
      CREATE_ARRAY1( scatter_clr, int, GDC_num_scatter_pts );
      gdPoint	ct[3];

      for( i=0; i<GDC_num_scatter_pts; ++i )
	{
	  int		hlf_scatterwdth = (int)( (float)(PX(2)-PX(1))
						 * (((float)((GDC_scatter+i)->width)/100.0)/2.0) );
	  int	scat_x = PX( (GDC_scatter+i)->point + (do_bar?1:0) ),
	    scat_y = PY( (GDC_scatter+i)->val );

	  if( (GDC_scatter+i)->point >= num_points ||				/* invalid point */
	      (GDC_scatter+i)->point <  0 )
	    continue;
	  scatter_clr[i] = clrallocate( im, (GDC_scatter+i)->color );

	  switch( (GDC_scatter+i)->ind )
	    {
	    case GDC_SCATTER_CIRCLE:
	      {
		long	uniq_clr = get_uniq_color( im );
		int		s        = 0,
		  e        = 360,
		  fo       = 0;

		if( !do_bar ){
		  if( (GDC_scatter+i)->point == 0 )
		    { s = 270; e = 270+180; fo = 1; }
		  else
		    if( (GDC_scatter+i)->point == num_points-1 )
		      { s = 90; e = 90+180; fo = -1; }
		}
		if( uniq_clr != -1L )							/* the safe way */
		  {
		    int	uc = gdImageColorAllocate( im, l2gdcal(uniq_clr) );
		    gdImageArc( im, scat_x, scat_y,
				hlf_scatterwdth*2, hlf_scatterwdth*2,
				s, e,
				uc );
		    if( fo )									/* close off  semi-circle case */
		      gdImageLine( im, scat_x, scat_y+hlf_scatterwdth,
				   scat_x, scat_y-hlf_scatterwdth,
				   uc );
		    gdImageFillToBorder( im, scat_x+fo, scat_y, uc, scatter_clr[i] );
		    gdImageArc( im, scat_x, scat_y,
				hlf_scatterwdth*2, hlf_scatterwdth*2,
				s, e,
				scatter_clr[i] );
		    if( fo )
		      gdImageLine( im, scat_x, scat_y+hlf_scatterwdth,
				   scat_x, scat_y-hlf_scatterwdth,
				   scatter_clr[i] );
		    gdImageColorDeallocate( im, uc );
		  }
		else											/* chance it */
		  {
		    gdImageArc( im, scat_x, scat_y,
				hlf_scatterwdth*2, hlf_scatterwdth*2,
				s, e,
				scatter_clr[i] );
		    if( fo )
		      gdImageLine( im, scat_x, scat_y+hlf_scatterwdth,
				   scat_x, scat_y-hlf_scatterwdth,
				   scatter_clr[i] );
		    gdImageFillToBorder( im, scat_x+fo, scat_y,
					 scatter_clr[i], scatter_clr[i] );
		  }
	      }
	      break;
	    case GDC_SCATTER_TRIANGLE_UP:
	      ct[0].x = scat_x;
	      ct[0].y = scat_y;
	      ct[1].x = scat_x - hlf_scatterwdth;
	      ct[1].y = scat_y + hlf_scatterwdth;;
	      ct[2].x = scat_x + hlf_scatterwdth;
	      ct[2].y = scat_y + hlf_scatterwdth;
	      if( !do_bar ){
		if( (GDC_scatter+i)->point == 0 )
		  ct[1].x = scat_x;
		else
		  if( (GDC_scatter+i)->point == num_points-1 )
		    ct[2].x = scat_x;
	      }
	      gdImageFilledPolygon( im, ct, 3, scatter_clr[i] );
	      break;
	    case GDC_SCATTER_TRIANGLE_DOWN:
	      ct[0].x = scat_x;
	      ct[0].y = scat_y;
	      ct[1].x = scat_x - hlf_scatterwdth;
	      ct[1].y = scat_y - hlf_scatterwdth;;
	      ct[2].x = scat_x + hlf_scatterwdth;
	      ct[2].y = scat_y - hlf_scatterwdth;
	      if( !do_bar ){
		if( (GDC_scatter+i)->point == 0 )
		  ct[1].x = scat_x;
		else
		  if( (GDC_scatter+i)->point == num_points-1 )
		    ct[2].x = scat_x;
	      }
	      gdImageFilledPolygon( im, ct, 3, scatter_clr[i] );
	      break;
	    }
	}
      FREE_ARRAY1( scatter_clr );
    }

  /* box it off */
  /*  after plotting so the outline covers any plot lines */
  {
    if( GDC_border == GDC_BORDER_ALL || (GDC_border & GDC_BORDER_X) )
      gdImageLine( im,          PX(0),   PY(lowest), PX(num_points-1+(do_bar?2:0)),  PY(lowest), LineColor );

    if( GDC_border == GDC_BORDER_ALL || (GDC_border & GDC_BORDER_TOP) )
      {
	setno = set_depth;
	gdImageLine( im,          PX(0),   PY(highest), PX(num_points-1+(do_bar?2:0)),  PY(highest), LineColor );
	setno = 0;
      }
  }
  if( GDC_border )
    {
      int	x1, y1, x2, y2;

      x1 = PX(0);
      y1 = PY(highest);
      x2 = PX(num_points-1+(do_bar?2:0));
      y2 = PY(lowest);
      if( GDC_border == GDC_BORDER_ALL || (GDC_border & GDC_BORDER_Y) )
	gdImageLine( im, x1, PY(lowest), x1, y1, LineColor );

      setno = set_depth;
      if( GDC_border == GDC_BORDER_ALL || (GDC_border & GDC_BORDER_Y) || (GDC_border & GDC_BORDER_TOP) )
	gdImageLine( im, x1, y1, PX(0), PY(highest), LineColor );
      /* if( !GDC_grid || do_vol || GDC_thumbnail )					// grid leaves right side Y open */
	{
	  if( GDC_border == GDC_BORDER_ALL || (GDC_border & GDC_BORDER_X) || (GDC_border & GDC_BORDER_Y2) )
	    gdImageLine( im, x2, y2, PX(num_points-1+(do_bar?2:0)), PY(lowest), LineColor );
	  if( GDC_border == GDC_BORDER_ALL || (GDC_border & GDC_BORDER_Y2) )
	    gdImageLine( im, PX(num_points-1+(do_bar?2:0)), PY(lowest),
			 PX(num_points-1+(do_bar?2:0)), PY(highest), LineColor );
	}
      setno = 0;
    }

  if( GDC_0Shelf && threeD &&								/* front of 0 shelf */
      ( (lowest < 0.0 && highest > 0.0) ||
	( (lowest == 0.0 || highest == 0.0) && !(GDC_border&GDC_BORDER_X) ) ) )
    {
      int	x2 = PX( num_points-1+(do_bar?2:0) ),
	y2 = PY( 0 );

      gdImageLine( im, PX(0), PY(0), x2, y2, LineColor );		/* front line */
      setno = set_depth;
      /* depth for 3Ds */
      gdImageLine( im, x2, y2, PX(num_points-1+(do_bar?2:0)), PY(0), LineColor );
      setno = 0;												/* set back to foremost */
    }

  if( GDC_annotation )			/* front half of annotation line */
    {
      int		x1 = PX(GDC_annotation->point+(do_bar?1:0)),
	y1 = PY(highest);
      int		x2;
      /* front line */
      gdImageLine( im, x1, PY(lowest)+1, x1, y1, AnnoteColor );
      if( threeD )
	{												/* on back plane */
	  setno = set_depth;
	  x2 = PX(GDC_annotation->point+(do_bar?1:0));
	  /* prspective line */
	  gdImageLine( im, x1, y1, x2, PY(highest), AnnoteColor );
	}
      else												/* for 3D done with back line */
	{
	  x2 = PX(GDC_annotation->point+(do_bar?1:0));
	  gdImageLine( im, x1, y1, x1, y1-2, AnnoteColor );
	}
      /* line-to and note */
      if( *(GDC_annotation->note) )						/* any note? */
	{
	  if( GDC_annotation->point >= (num_points/2) )		/* note to the left */
	    {
	      gdImageLine( im, x2,              PY(highest)-2,
			   x2-annote_hgt/2, PY(highest)-2-annote_hgt/2,
			   AnnoteColor );
	      GDCImageStringNL( im,
				&GDC_fontc[GDC_annotation_font_size],
				gdc_annotation_font, gdc_annotation_ptsize,
				0.0,
				x2-annote_hgt/2-1-annote_len - 1,
				PY(highest)-annote_hgt+1,
				GDC_annotation->note,
				AnnoteColor,
				GDC_JUSTIFY_RIGHT,
				NULL );
	    }
	  else												/* note to right */
	    {
	      gdImageLine( im, x2,              PY(highest)-2,
			   x2+annote_hgt/2, PY(highest)-2-annote_hgt/2,
			   AnnoteColor );
	      GDCImageStringNL( im,
				&GDC_fontc[GDC_annotation_font_size],
				gdc_annotation_font, gdc_annotation_ptsize,
				0.0,
				x2+annote_hgt/2+1 + 1,
				PY(highest)-annote_hgt+1,
				GDC_annotation->note,
				AnnoteColor,
				GDC_JUSTIFY_LEFT,
				NULL );
	    }
	}
      setno = 0;
    }


  /* usually GDC_generate_img is used in conjunction with hard or hold options */
  if( GDC_generate_img )
    {
      fflush(img_fptr);			/* clear anything buffered  */
      switch( GDC_image_type )
	{
	case GDC_PNG:	gdImagePng( im, img_fptr );						break;
#ifdef HAVE_JPEG
	case GDC_JPEG:	gdImageJpeg( im, img_fptr, GDC_jpeg_quality );	break;
#endif
	case GDC_WBMP:	gdImageWBMP( im, PlotColor, img_fptr );			break;
	case GDC_GIF:
	default:		gdImageGif( im, img_fptr);
	}
    }

  if( bg_img )
    gdImageDestroy(bg_img);
  if( GDC_hold_img & GDC_EXPOSE_IMAGE )
    GDC_image = (void*)im;
  else
    gdImageDestroy(im);
  FREE_ARRAY1( uvals );
  FREE_ARRAY1( ExtVolColor );
  FREE_ARRAY2( ExtColor );
  FREE_ARRAY2( ExtColorShd );
  return 0;
}
