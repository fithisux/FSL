/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined (SLICEVISITOR_H)
#define SLICEVISITOR_H

#include "slice.h"

class SliceVisitor
{
public:
	SliceVisitor() {} 
	virtual ~SliceVisitor() {}

	virtual void visit(SliceB::WeakHandle target) = 0;
	virtual void visit(SliceUS::WeakHandle target) = 0;
	virtual void visit(SliceF::WeakHandle target) = 0;
};

#endif
