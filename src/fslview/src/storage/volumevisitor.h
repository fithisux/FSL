/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined (VOLUMEVISITOR_H)
#define VOLUMEVISITOR_H

#include "volume.h"

class VolumeVisitor
{
public:
  VolumeVisitor() {} 
  virtual ~VolumeVisitor() {}
  
  virtual void visit(VolumeB::WeakHandle target) = 0;
  virtual void visit(VolumeS::WeakHandle target) = 0;
  virtual void visit(VolumeF::WeakHandle target) = 0;
};

#endif
