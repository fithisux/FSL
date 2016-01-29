/*
 * gifencode.c
 *
 * Copyright (c) 1997,1998 by Hans Dinsen-Hansen
 * The algorithms are inspired by those of gifcode.c
 * Copyright (c) 1995,1996 Michael A. Mayer
 * All rights reserved.
 *
 * This software may be freely copied, modified and redistributed
 * without fee provided that above copyright notices are preserved
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
 *
 *           Implements GIF encoding by means of a tree search.
 *           --------------------------------------------------
 *
 *  - The string table may be thought of being stored in a "b-tree of
 * steroids," or more specifically, a {256,128,...,4}-tree, depending on
 * the size of the color map.
 *  - Each (non-NULL) node contains the string table index (or code) and
 * {256,128,...,4} pointers to other nodes.
 *  - For example, the index associated with the string 0-3-173-25 would be
 * stored in:
 *       first->node[0]->node[3]->node[173]->node[25]->code
 *
 *  - Speed and effectivity considerations, however, have made this
 * implementation somewhat obscure, because it is costly to initialize
 * a node-array where most elements will never be used.
 *  - Initially, a new node will be marked as terminating, TERMIN.
 * If this node is used at a later stage, its mark will be changed.
 *  - Only nodes with several used nodes will be associated with a
 * node-array.  Such nodes are marked LOOKUP.
 *  - The remaining nodes are marked SEARCH.  They are linked together
 * in a search-list, where a field, NODE->alt, points at an alternative
 * following color.
 *  - It is hardly feasible exactly to predict which nodes will have most
 * used node pointers.  The theory here is that the very first node as
 * well as the first couple of nodes which need at least one alternative
 * color, will be among the ones with many nodes ("... whatever that
 * means", as my tutor in Num. Analysis and programming used to say).
 *  - The number of possible LOOKUP nodes depends on the size of the color 
 * map.  Large color maps will have many SEARCH nodes; small color maps
 * will probably have many LOOKUP nodes.
*/

#include "whirlgif.h"

#define BLOKLEN 255
#define BUFLEN 1000


int chainlen = 0, maxchainlen = 0, nodecount = 0, lookuptypes = 0, nbits;
short need = 8;
GifTree *empty[256], GifRoot = {LOOKUP, 0, 0, empty, NULL, NULL},
	*topNode, *baseNode, **nodeArray, **lastArray;

extern unsigned int debugFlag, verbose;
extern int count;

void GifEncode(FILE *fout, UBYTE *pixels, int depth, int siz)
{
  GifTree *first = &GifRoot, *newNode, *curNode;
  UBYTE   *end;
  int     cc, eoi, next, tel=0;
  short   cLength;

  char    *pos, *buffer;

  empty[0] = NULL;
  need = 8;

  nodeArray = empty;
  memmove(++nodeArray, empty, 255*sizeof(GifTree **));
  if (( buffer = (char *)malloc((BUFLEN+1)*sizeof(char))) == NULL )
	 TheEnd1("No memory for writing");
  buffer++;


  pos = buffer;
  buffer[0] = 0x0;

  cc = (depth == 1) ? 0x4 : 1<<depth;
  fputc((depth == 1) ? 2 : depth, fout);
  eoi = cc+1;
  next = cc+2;

  cLength = (depth == 1) ? 3 : depth+1;

  if (( topNode = baseNode = (GifTree *)malloc(sizeof(GifTree)*4094)) == NULL )
	 TheEnd1("No memory for GIF-code tree");
  if (( nodeArray = first->node = (GifTree **)malloc(256*sizeof(GifTree *)*noOfArrays)) == NULL )
	 TheEnd1("No memory for search nodes");
  lastArray = nodeArray + ( 256*noOfArrays - cc);
  ClearTree(cc, first);

  pos = AddCodeToBuffer(cc, cLength, pos);

  end = pixels+siz;
  curNode = first;
  while(pixels < end) {

    if ( curNode->node[*pixels] != NULL ) {
      curNode = curNode->node[*pixels];
      tel++;
      pixels++;
      chainlen++;
      continue;
    } else if ( curNode->typ == SEARCH ) {
      newNode = curNode->nxt;
      while ( newNode->alt != NULL ) {
	if ( newNode->ix == *pixels ) break;
	newNode = newNode->alt;
      }
      if (newNode->ix == *pixels ) {
	tel++;
	pixels++;
	chainlen++;
	curNode = newNode;
	continue;
      }
    }

/* ******************************************************
 * If there is no more thread to follow, we create a new node.  If the
 * current node is terminating, it will become a SEARCH node.  If it is
 * a SEARCH node, and if we still have room, it will be converted to a
 * LOOKUP node.
*/
  newNode = ++topNode;
  switch (curNode->typ ) {
   case LOOKUP:
     newNode->nxt = NULL;
     newNode->alt = NULL,
     curNode->node[*pixels] = newNode;
   break;
   case SEARCH:
     if ( nodeArray != lastArray ) {
       nodeArray += cc;
       curNode->node = nodeArray;
       curNode->typ = LOOKUP;
       curNode->node[*pixels] = newNode;
       curNode->node[(curNode->nxt)->ix] = curNode->nxt;
       lookuptypes++;
       newNode->nxt = NULL;
       newNode->alt = NULL,
       curNode->nxt = NULL;
       break;
     }
/*   otherwise do as we do with a TERMIN node  */
   case TERMIN:
     newNode->alt = curNode->nxt;
     newNode->nxt = NULL,
     curNode->nxt = newNode;
     curNode->typ = SEARCH;
     break;
   default:
     fprintf(stderr, "Silly node type: %d\n", curNode->typ);
  }
  newNode->code = next;
  newNode->ix = *pixels;
  newNode->typ = TERMIN;
  newNode->node = empty;
  nodecount++;
/*
* End of node creation
* ******************************************************
*/
  if (debugFlag) {
    if (curNode == newNode) fprintf(stderr, "Wrong choice of node\n");
    if ( curNode->typ == LOOKUP && curNode->node[*pixels] != newNode ) fprintf(stderr, "Wrong pixel coding\n");
    if ( curNode->typ == TERMIN ) fprintf(stderr, "Wrong Type coding; frame no = %d; pixel# = %d; nodecount = %d\n", count, tel, nodecount);
  }
    pos = AddCodeToBuffer(curNode->code, cLength, pos);
    if ( chainlen > maxchainlen ) maxchainlen = chainlen;
    chainlen = 0;
    if(pos-buffer>BLOKLEN) {
      buffer[-1] = (char)BLOKLEN;
      fwrite(buffer-1, 1, BLOKLEN+1, fout);
      buffer[0] = buffer[BLOKLEN];
      buffer[1] = buffer[BLOKLEN+1];
      buffer[2] = buffer[BLOKLEN+2];
      buffer[3] = buffer[BLOKLEN+3];
      pos -= BLOKLEN;
    }
    curNode = first;

    if(next == (1<<cLength)) cLength++;
    next++;

    if(next == 0xfff) {
      ClearTree(cc,first);
      pos = AddCodeToBuffer(cc, cLength, pos);
      if(pos-buffer>BLOKLEN) {
	buffer[-1] = (char)BLOKLEN;
	fwrite(buffer-1, 1, BLOKLEN+1, fout);
	buffer[0] = buffer[BLOKLEN];
	buffer[1] = buffer[BLOKLEN+1];
	buffer[2] = buffer[BLOKLEN+2];
	buffer[3] = buffer[BLOKLEN+3];
	pos -= BLOKLEN;
      }
      next = cc+2;
      cLength = (depth == 1)?3:depth+1;
    }
  }

  pos = AddCodeToBuffer(curNode->code, cLength, pos);
  if(pos-buffer>BLOKLEN-3) {
    buffer[-1] = (char)(BLOKLEN-3);
    fwrite(buffer-1, 1, BLOKLEN-2, fout);
    buffer[0] = buffer[BLOKLEN-3];
    buffer[1] = buffer[BLOKLEN-2];
    buffer[2] = buffer[BLOKLEN-1];
    buffer[3] = buffer[BLOKLEN];
    buffer[4] = buffer[BLOKLEN+1];
    pos -= BLOKLEN-3;
  }
  pos = AddCodeToBuffer(eoi, cLength, pos);
  pos = AddCodeToBuffer(0x0, -1, pos);
  buffer[-1] = pos-buffer;

  fwrite(buffer-1, pos-buffer+1, 1, fout);
  free(buffer-1); free(first->node); free(baseNode);
  if (debugFlag) fprintf(stderr, "pixel count = %d; nodeCount = %d lookup nodes = %d\n", tel, nodecount, lookuptypes);
  return;

}

void ClearTree(int cc, GifTree *root)
{
  int i;
  GifTree *newNode, **xx;

  if (debugFlag>1) fprintf(stderr, "Clear Tree  cc= %d\n", cc);
  if (debugFlag>1) fprintf(stderr, "nodeCount = %d lookup nodes = %d\n", nodecount, lookuptypes);
  maxchainlen=0; lookuptypes = 1;
  nodecount = 0;
  nodeArray = root->node;
  xx= nodeArray;
  for (i = 0; i < noOfArrays; i++ ) {
    memmove (xx, empty, 256*sizeof(GifTree **));
    xx += 256;
  }
  topNode = baseNode;
  for(i=0; i<cc; i++) {
    root->node[i] = newNode = ++topNode;
    newNode->nxt = NULL;
    newNode->alt = NULL;
    newNode->code = i;
    newNode->ix = i;
    newNode->typ = TERMIN;
    newNode->node = empty;
    nodecount++;
  }
}

char *AddCodeToBuffer(int code, short n, char *buf)
{
  int    mask;

  if(n<0) {
    if(need<8) {
      buf++;
      *buf = 0x0;
    }
    need = 8;
    return buf;
  }

  while(n>=need) {
    mask = (1<<need)-1;
    *buf += (mask&code)<<(8-need);
    buf++;
    *buf = 0x0;
    code = code>>need;
    n -= need;
    need = 8;
  }
  if(n) {
    mask = (1<<n)-1;
    *buf += (mask&code)<<(8-need);
    need -= n;
  }
  return buf;
}
