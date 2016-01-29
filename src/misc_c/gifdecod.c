/*
 * gifdecode.c
 *
 * Copyright (c) 1997,1998,1999 by Hans Dinsen-Hansen
 * Partly inspired by Helene Schulerud's inversrl_imageb.c of 13.2.1991.
 * Copyright (c) 1995,1996 by Kevin Kadow
 * Copyright (c) 1990,1991,1992,1993 by Mark Podlipec
 * All rights reserved.
 *
 * This software may be freely copied, modified and redistributed
 * without fee provided that this copyright notice is preserved
 * intact on all copies and modified copies.
 *
 * There is no warranty or other guarantee of fitness of this software.
 * It is provided solely "as is". The author(s) disclaim(s) all
 * responsibility and liability with respect to this software's usage
 * or its effect upon hardware or computer systems.
 *
 * The Graphics Interchange format (c) is the Copyright property of
 * Compuserve Incorporated.  Gif(sm) is a Service Mark property of
 * Compuserve Incorporated.
 *
 */

#include <stdio.h>
#include "whirlgif.h"

extern unsigned int debugFlag, verbose;
extern int count;

ULONG codeSize, expected, imgsize, mask, old, picI, rootCodeSize, first[MAXVAL];
UBYTE *topGifBuff, *picture, gifBuff[MAXVAL], last[MAXVAL];
int imgheight, imgwidth, interlaced, pass,
	step[5]={7,7,3,1,0}, start[5]= {0,4,2,1,0};

void GifDecode(FILE *fp, UBYTE *pix, GifImageHdr gifimage)
{ UBYTE *chPos, firstCodeOut = 0, charBuff[256];
  ULONG CLEAR, EOI, bits = 0, code = 0, codeFeched, buffCount = 0;
  int need = 0;

  interlaced = gifimage.i;
  imgwidth = gifimage.width;
  imgheight = gifimage.height;
  imgsize = imgwidth * imgheight;
  picture = pix;
  pass = picI = 0;
  if ( debugFlag > 1 )
     fprintf(stderr, "(interlaced,imgwidth,imgheight,imgsize)=(%d,%d,%d,%d)\n",
		interlaced, imgwidth, imgheight, imgsize);

  rootCodeSize = Xgetc(fp);
  CLEAR = 1 << rootCodeSize;
  EOI = CLEAR + 1;

  GifClearTable();

  if ( (buffCount = Xgetc(fp)) == 0 ) {
    sprintf( (char *)charBuff, "End of image # %d before it began!\n", count);
    TheEnd1( (char *)charBuff);
  }

  while(buffCount > 0) {
    if ( fread(charBuff, 1, buffCount, fp) != buffCount ) { 
      sprintf( (char *)charBuff, "Premature end of file; Image # %d\n", count);
      TheEnd1( (char *)charBuff);
    }
    for(chPos = charBuff; buffCount-- > 0; chPos++) {
      need += (int) *chPos << bits;
      bits += 8;
      while (bits >= codeSize) {
	code = need & mask;
	need >>= codeSize;
	bits -= codeSize;
	if(code > expected)
	   TheEnd1("Neither LZW nor RunLength (new code != expected)\n");
	if (code == EOI) {
	  if (debugFlag > 1) fprintf(stderr, "Image # %d ends; buffCount=%d\n",
			count, buffCount);
	  goto skipRest;
	}

	if(code == CLEAR) {
	GifClearTable();
	continue;
	}
	if(old == MAXVAL) {  /* i.e. first code after clear table */
	  pix = GifSendData(pix, 1, &last[code]);
	  firstCodeOut = last[code];
	  old = code;
	  continue;
	}
	codeFeched = code;
	if(code == expected) {
	  *topGifBuff++ = firstCodeOut;
	  code = old;
	}
	while(code > CLEAR) {
	  *topGifBuff++ = last[code];
	  code = first[code];
	}

	*topGifBuff++ = firstCodeOut = last[code];
	first[expected] = old;
	last[expected] = firstCodeOut;
	if(expected < MAXVAL) expected++;
	if(((expected & mask) == 0) && (expected < MAXVAL)) {
	  codeSize++;
	  mask += expected;
	}
	old = codeFeched;
	pix = GifSendData(pix, topGifBuff - gifBuff, gifBuff);
	topGifBuff = gifBuff;

      }   /* end of extracting codes */
    }   /* end of reading a block of data */
    if ( (buffCount = Xgetc(fp)) == 0 ) {
      sprintf( (char *)charBuff, "End of image # %d without EOI\n", count);
      TheEnd1( (char *)charBuff);
    }
  }

skipRest: 
  if (debugFlag) fprintf(stderr, "Ending GifDecode, written: %ld=%d\n",
	  interlaced && (pix-picture == 0) ? imgsize : pix - picture, imgsize);
  return ;
}

UBYTE *GifSendData(UBYTE *pix, int bytes, UBYTE source[])

{ int j=0;
  for(j = bytes - 1; j >= 0; j--) {
    picI++;
    *pix = source[j]; pix++;
    if ( interlaced && (picI % imgwidth == 0) ) {
      picI += ( imgwidth * step[pass]);
      if (picI >= imgsize) {
	picI = start[++pass] * imgwidth;
	if ( debugFlag > 1 )
	fprintf(stderr, "De-interlacing: (picI,pass,start[pass])=(%d,%d,%d)\n",
				picI, pass, start[pass]);
      }
      pix = &picture[picI];
    }
  }
  return(pix);
}

void GifClearTable()
{ int i, maxi;
  maxi = 1 << rootCodeSize;
  expected = maxi + 2;
  if (debugFlag > 1 ) fprintf(stderr, "Initing Table...");
  old = MAXVAL;
  codeSize = rootCodeSize + 1;
  mask = (1<<codeSize) - 1;

  for(i = 0; i < maxi; i++) {
    first[i] = MAXVAL;
    last[i] = i & 0xff;
  }
  topGifBuff = gifBuff;
}
