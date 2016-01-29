/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Rama Aravind Vorray
		James Saunders
		David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "vtktoolbar.h"
#include "vtkwidget.h"

#include <qspinbox.h>

VTKToolbar::VTKToolbar(QWidget *parent, VTKProperties& p): 
  VTKToolbarBase(parent), m_props(p) 
{
  m_threshold->setValue(m_props.inqLowerThreshold());
}

void VTKToolbar::thresholdValueChanged(int v)
{
  m_props.setLowerThreshold(v);
}  

void VTKToolbar::clippingStateChanged(int state)
{
  m_props.setClipping(state);
}
