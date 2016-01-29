/* GDCHART 0.10.0dev  GDC.C  2 Nov 2000 */
/* Copyright Bruce Verderaime 1998, 1999, 2000 */

#define GDC_INCL
#define GDC_LIB
#include <math.h>
#include "gdc.h"

struct	GDC_FONT_T	GDC_fontc[GDC_numfonts] = { { (gdFontPtr)NULL, 8,  5 } ,
						    { (gdFontPtr)NULL, 8,  5 } ,
						    { (gdFontPtr)NULL, 12, 6 } ,
						    { (gdFontPtr)NULL, 13, 7 } ,
						    { (gdFontPtr)NULL, 16, 8 } ,
						    { (gdFontPtr)NULL, 15, 9 } };

/* ------------------------------------------------------------------- *\ 
 * convert from enum GDC_font_size to gd fonts
 * for now load them all
 *	#defines and #ifdefs might enable loading only needed fonts
 *	gd2.0 is to be built as a shared obj.
\* ------------------------------------------------------------------- */
void
load_font_conversions()
{
	GDC_fontc[GDC_pad].f     = gdFontTiny;
	GDC_fontc[GDC_TINY].f    = gdFontTiny;
	GDC_fontc[GDC_SMALL].f   = gdFontSmall;
	GDC_fontc[GDC_MEDBOLD].f = gdFontMediumBold;
	GDC_fontc[GDC_LARGE].f   = gdFontLarge;
	GDC_fontc[GDC_GIANT].f   = gdFontGiant;
}

/* ------------------------------------------------------------------ *\ 
 * count (natural) substrings (new line sep)
\* ------------------------------------------------------------------ */
short
cnt_nl( char	*nstr,
		int		*len )			/* strlen - max seg */
{
	short	c           = 1;
	short	max_seg_len = 0;
	short	tmplen      = 0;

	if( !nstr )
		{
		if( len )
			*len = 0;
		return 0;
		}
	while( *nstr )
		{
		if( *nstr == '\n' )
			{
			++c;
			max_seg_len = MAX( tmplen, max_seg_len );
			tmplen = 0;
			}
		else
			++tmplen;
		++nstr;
		}

	if( len )
		*len = MAX( tmplen, max_seg_len );		/* don't forget last seg */
	return c;
}

/* ------------------------------------------------------------------ *\ 
 * gd out a string with '\n's
 * handle FTs (TTFs) and gd fonts
 * gdImageString() draws from the upper left;
 * gdImageStringFT() draws from lower left (one font height, even with '\n's)! >:-|
\* ------------------------------------------------------------------ */
int
GDCImageStringNL( gdImagePtr		im,
				  struct GDC_FONT_T	*f,
				  char				*ftfont,
				  double			ftptsz,
				  double			rad,
				  int				x,
				  int				y,
				  char				*str,
				  int				clr,
				  GDC_justify_t		justify,
				  char				**sts )
{
	int		retval = 0;
	char	*err   = NULL;

#ifdef HAVE_LIBFREETYPE
	/* TODO: honor justifies */
	if( ftfont && ftptsz ) 
		{
		/* need one line height */
		/* remember last one (will likely be the same) */
		/*	is this needed? */
		/*	gdImageStringFT() utilizes some caching */
		/*	saves a couple floating point trig calls */
		static int		f1hgt = 0;
		static double	xs,
						ys;
		static double	lftptsz = 0.0;
		static char		*lftfont = (char*)-1;

		if( !f1hgt ||
			( lftfont != ftfont || lftptsz != ftptsz ) )
			{
			f1hgt = GDCfnt_sz( "Aj",
							   0,
							   ftfont,
							   ftptsz,
							   rad,
							   NULL ).h;
			xs = (double)f1hgt * sin(rad);
			ys = (double)(f1hgt-1) * cos(rad);
			}
		x += (int)xs;
		y += (int)ys;
		if( (err = gdImageStringFT( im,
									(int*)NULL,
									clr,
									ftfont,
									ftptsz,
									rad,
									x,
									y,
									str)) == NULL )
			{
			if( sts )	*sts = err;
			return 0;
			}
		else
			{
			/* TTF failed */
			retval = 1;
			/* fall through - default to gdFonts */
			/* reinstate upper left reference */
			x -= (int)xs;
			y -= (int)ys;
			}
		}
#endif

	{
	int		i;
	int		len;
	int     max_len;
	short   strs_num = cnt_nl( str, &max_len );
	CREATE_ARRAY1( sub_str, unsigned char, max_len+1 );	/* char sub_str[max_len+1]; */

	len      = -1;
	strs_num = -1;
	i = -1;
	do
		{
		++i;
		++len;
		sub_str[len] = *(str+i);
		if( *(str+i) == '\n' ||
			*(str+i) == '\0' )
			{
			int	xpos;

			sub_str[len] = '\0';
			++strs_num;
			switch( justify )
			  {
			  case GDC_JUSTIFY_LEFT:	xpos = 0;					break;
			  case GDC_JUSTIFY_RIGHT:	xpos = f->w*(max_len-len);	break;
			  case GDC_JUSTIFY_CENTER:
			  default:					xpos = f->w*(max_len-len)/2;
			  }
			if( rad == 0.0 )
				gdImageString( im,
							   f->f,
							   x + xpos,
							   y + (f->h-1)*strs_num,
							   sub_str,
							   clr );
			else /* if( rad == M_PI/2.0 ) */
				gdImageStringUp( im,
								 f->f,
								 x + (f->h-1)*strs_num,
								 y - xpos,
								 sub_str,
								 clr );
			len = -1;
			}
		}
	while( *(str+i) );
	}

	if( sts )	*sts = err;
	return retval;
}

/* ------------------------------------------------------------------------ *\ 
 * TODO:                                                                    *
 * really get a unique color from the color map                             *
\* ------------------------------------------------------------------------ */
long
get_uniq_color( gdImagePtr	im )
{
	return 0x123454;
}

/* ------------------------------------------------------------------------ */
struct fnt_sz_t
GDCfnt_sz( char					*s,
		   enum GDC_font_size	gdfontsz,
		   char					*ftfont,
		   double				ftfptsz,
		   double				rad,		/* w,h still relative to horiz. */
		   char					**sts )
{
	struct fnt_sz_t	rtnval;
	int				len;
	char			*err = NULL;

#ifdef HAVE_LIBFREETYPE
	if( ftfont && ftfptsz )
		{
		int		brect[8];

		/* obtain brect so that we can size the image */
		if( (err = gdImageStringFT( (gdImagePtr)NULL,
									&brect[0],
									0,
									ftfont,
									ftfptsz,
									0.0,	/* rad, */	/* always match angled size??? */
									0,
									0,
									s)) == NULL )
			{
			rtnval.h = brect[1] - brect[7];
			rtnval.w = brect[2] - brect[0];
			if( sts )	*sts = err;
			return rtnval;
			}
		}
#endif

	rtnval.h = cnt_nl(s,&len) * GDC_fontc[gdfontsz].h;
	rtnval.w = len * GDC_fontc[gdfontsz].w;
	if( sts )	*sts = err;
	return rtnval;
}

/* ------------------------------------------------------------------------ */
void
GDC_destroy_image(void *im)
{
	if( im )
		gdImageDestroy( (gdImagePtr)im );
}

/* ------------------------------------------------------------------------ */
void
out_err( int			IMGWIDTH,
		 int			IMGHEIGHT,
		 FILE			*fptr,
		 unsigned long	BGColor,
		 unsigned long	LineColor,
		 char			*err_str )
{

	gdImagePtr	im;
	int			lineclr;
	int			bgclr;


	if( (GDC_hold_img & GDC_REUSE_IMAGE) &&
		GDC_image != (void*)NULL )
		im = GDC_image;
	else
		im = gdImageCreate( IMGWIDTH, IMGHEIGHT );

	bgclr    = gdImageColorAllocate( im, l2gdcal(BGColor) );
	lineclr  = gdImageColorAllocate( im, l2gdcal(LineColor) );

	gdImageString( im,
				   gdFontMediumBold,
				   IMGWIDTH/2 - GDC_fontc[GDC_MEDBOLD].w*strlen(err_str)/2,
				   IMGHEIGHT/3,
				   (unsigned char*)err_str,
				   lineclr );

	/* usually GDC_generate_img is used in conjunction with hard or hold options */
	if( GDC_generate_img )
		{
		fflush(fptr);			/* clear anything buffered  */
		switch( GDC_image_type )
			{
			case GDC_PNG:	gdImagePng( im, fptr );						break;
#ifdef HAVE_JPEG
			case GDC_JPEG:	gdImageJpeg( im, fptr, GDC_jpeg_quality );	break;
#endif
			case GDC_WBMP:	gdImageWBMP( im, lineclr, fptr );			break;
			case GDC_GIF:
			default:		gdImageGif( im, fptr);
			}
		}

	if( GDC_hold_img & GDC_EXPOSE_IMAGE )
		GDC_image = (void*)im;
	else
		gdImageDestroy(im);
	return;
}
