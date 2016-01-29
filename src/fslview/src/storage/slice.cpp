/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

// slice.cpp: implementation of the Slice class.
//
//////////////////////////////////////////////////////////////////////

#include "slice.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

#include "slice.hpp"

template SliceStore<unsigned char>;
template SliceStore<unsigned short>;
template SliceStore<float>;
