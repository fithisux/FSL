/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Rama Aravind Vorray
		James Saunders
		David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(VTKTOOLBAR_H)
#define VTKTOOLBAR_H

#include "vtktoolbarbase.h"

class VTKProperties;

class VTKToolbar: public VTKToolbarBase
{
public:
  VTKToolbar(QWidget *, VTKProperties&);

private slots:
  void thresholdValueChanged(int);
  void clippingStateChanged(int);

private:
  VTKProperties &m_props;
};

#endif
