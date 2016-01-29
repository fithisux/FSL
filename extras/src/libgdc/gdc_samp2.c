/* GDCHART 0.10.0dev  2nd CHART SAMPLE  2 Nov 2000 */
/* Copyright Bruce Verderaime 1998, 1999, 2000 */

/*
** vi note  :set tabstop=4 **

 a more complicated example
 High Low Close Combo (Volume)  with annotation

 produces a file: g2.gif

 Until a README is ready, see gdchart.h for options
	All options are defaulted, no need to set any
*/

#include <stdio.h>

#include "gdc.h"
#include "gdchart.h"


main()
{
	/* set some sample data points */
	float	h[12]  = {	17.8,  17.1,  17.3,  GDC_NOVALUE,  17.2,  17.1,
						17.3,  17.3,  17.3,  17.1,         17.5,  17.4 };

	float	c[12]  =  { 17.0,  16.8,  16.9,  GDC_NOVALUE,  16.9,  16.8,
					    17.2,  16.8,  17.0,  16.9,         16.4,  16.1 };

	float	l[12]  = {  16.8,  16.8,  16.7,  GDC_NOVALUE,  16.5,  16.0,
						16.1,  16.8,  16.5,  16.9,         16.2,  16.0 };

	float	v[12]  = {  150.0, 100.0, 340.0,  GDC_NOVALUE, 999.0, 390.0,
						420.0, 150.0, 100.0,  340.0,       1590.0, 700.0 };
	char	*t[12] = {	"May", "Jun", "Jul",  "Aug",       "Sep",  "Oct",
						"Nov", "Dec", "Jan", "Feb",        "Mar",  "Apr" };

	/* set color RGB as ulong array */
	unsigned long	setcolor[3]    = { 0xC0C0FF, 0xFF4040, 0xFFFFFF }; 

	GDC_ANNOTATION_T	anno;

	/* need an open FILE pointer  - can be stdout */
	FILE				*outpng1 = fopen( "g2.png", "wb" );	/* rem: test open() fail */

	anno.color = 0x00FF00;
	strncpy( anno.note, "Did Not\nTrade", MAX_NOTE_LEN );	/* don't exceed MAX_NOTE_LEN */
	anno.point = 3;											/* first is 0 */
	GDC_annotation_font_size = GDC_TINY;
	GDC_annotation = &anno;									/* set annote option */

	GDC_HLC_style = GDC_HLC_I_CAP | GDC_HLC_CLOSE_CONNECTED;
	GDC_HLC_cap_width = 45;

	GDC_bar_width     = 15;									/* % */

	GDC_thumbnail = TRUE;
	/*	GDC_BGImage = "W.gif"; */

	GDC_title = "Widget Corp.";
	GDC_ytitle = "Price ($)";
	GDC_ytitle2 = "Volume (K)";
	GDC_ytitle_size = GDC_SMALL;
	GDC_VolColor = 0xAAAAFFL;								/* aka combo */
	GDC_3d_depth  = 4.0;									/* % entire gif */

	/*	GDC_SetColor  = setcolor;						*/		/* see README */ GDC_image_type     = GDC_PNG;
	GDC_PlotColor = 0xFFFFFFL;
	GDC_grid = FALSE;

/*	GDC_xtitle="fy.1998"; */

/*	fprintf( stdout, "Content-Type: image/gif\n\n" );*/		/* rem: for web use */
															/* finally: make the call */
	out_graph( 200, 275,									/* overall width, height */
			   outpng1,										/* open FILE pointer */
			   GDC_COMBO_LINE_BAR,							/* chart type */
			   12,											/* number of points */
		   t,											/* X axis label array */
			   2,											/* number of sets (see README) */
			   											/* set 1 (high) */
			   l,											/* low */
			   c,											/* close */
			   v );											/* combo/volume */

	fclose( outpng1 );
	exit(0);
}
