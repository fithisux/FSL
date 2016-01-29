/*
 * whirlgif.h
 *
 * Copyright (c) 1997,1998,1999 by Hans Dinsen-Hansen
 * Copyright (c) 1995,1996 by Kevin Kadow
 * Copyright (c) 1990,1991,1992 by Mark Podlipec.
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

#define DA_REV  3.04

/* common includes */
#include <stdio.h>
#include <stdlib.h>

#ifdef _USE_STRINGS_H
#include <strings.h>
#else
#include <string.h>
#endif

#ifdef _FOPEN_TXT_OR_BIN
#define WRIBIN	"wb"
#define REATXT	"rt"
#define REABIN	"rb"
#else
/* Usually there is no need to distinguish between binary and txt */
#define WRIBIN	"w"
#define REATXT	"r"
#define REABIN	"r"
#endif

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/* define constants and defaults */
    /* Default amount of inter-frame time */
#define DEFAULT_TIME 10
    /* If set to 1, Netscape 'loop' code will be added by default */
#define DEFAULT_LOOP 0
    /* If set to 1, use the colormaps from all images, not just the first */
#define DEFAULT_USE_COLORMAP 0

    /* Used in calculating the transparent color */
#define TRANS_NONE 1
#define TRANS_RGB 2
#define TRANS_MAP 3

#define DISP_NONE 0
#define DISP_NOT  1
#define DISP_BACK 2
#define DISP_PREV 3
#define DEFAULT_DISPOSAL DISP_NONE
    /* set default disposal method here to any of the DISP_XXXX values */

#define BIGSTRING 256
#define MAXVAL  4096        /* maxval of lzw coding size */
#define MAXVALP 4096
#define TERMIN 'T'
#define LOOKUP 'L'
#define SEARCH 'S'
#define noOfArrays 20
/* defines the amount of memory set aside in the encoding for the
 * LOOKUP type nodes; for a 256 color GIF, the number of LOOKUP
 * nodes will be <= noOfArrays, for a 128 color GIF the number of
 * LOOKUP nodes will be <= 2 * noOfArrays, etc.  */

/* define shorthand for various types */
#define LONG int
#define ULONG unsigned int
#define BYTE char
#define UBYTE unsigned char
#define SHORT short
#define USHORT unsigned short
#define WORD short int
#define UWORD unsigned short int


/* definition of various structures */
typedef struct Transparency {
  int type;
  UBYTE valid;
  UBYTE map;
  UBYTE red;
  UBYTE green;
  UBYTE blue;
  } Transparency;

typedef struct Global {
  Transparency trans;
  int left;
  int top;
  unsigned int time;
  unsigned short disposal;
  } Global;

typedef struct GifScreenHdr {
  int width;
  int height;
  UBYTE m;
  UBYTE cres;
  UBYTE pixbits;
  UBYTE bc;
  UBYTE aspect;
 } GifScreenHdr;

typedef union GifColor {
  struct cmap {
    UBYTE red;
    UBYTE green;
    UBYTE blue;
    UBYTE pad;
   } cmap;
  ULONG pixel;
 } GifColor;

typedef struct GifImageHdr {
  int left;
  int top;
  int width;
  int height;
  UBYTE m;
  UBYTE i;
  UBYTE pixbits;
  UBYTE reserved;
 } GifImageHdr;

typedef struct GifTree {
  char typ;             /* terminating, lookup, or search */
  int code;             /* the code to be output */
  UBYTE ix;             /* the color map index */
  struct GifTree **node, *nxt, *alt;
} GifTree;

/* define inline functions */
#define GifPutShort(i, fout)    {fputc(i&0xff, fout); fputc(i>>8, fout);}
#define GifGetShort(fin)        (Xgetc(fin) | Xgetc(fin)<<8)

/* forward declaration of the functions  */
char *AddCodeToBuffer(int, short, char *);
void CalcTrans(char *);
void ClearTree(int, GifTree *);
void GifClearTable();
void GifComment(FILE *, char *);
void GifDecode(FILE *, UBYTE *, GifImageHdr);
void GifEncode(FILE *, UBYTE *, int, int);
void GifLoop(FILE *, unsigned int);
void GifReadFile(FILE *, char *, int);
void GifScreenHeader(FILE *, FILE *, int);
UBYTE *GifSendData(UBYTE *, int, UBYTE *);
void ReadImageHeader(FILE *);
void SetOffset(char *);
long sq(UBYTE, UBYTE);
void TheEnd();
void TheEnd1(char *);
void Usage();
void WriteImageHeader(FILE *);
UBYTE Xgetc(FILE *);
