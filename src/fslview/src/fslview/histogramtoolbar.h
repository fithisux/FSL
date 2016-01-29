/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Rama Aravind Vorray
		James Saunders
		David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(HISTOGRAMTOOLBAR_H)
#define HISTOGRAMTOOLBAR_H

#include "histogramtoolbarbase.h"

//! @brief Extends HistogramToolbarBase to provide HistogramToolbar
//! behaviour.
//!
//! @author Dave Flitney
class HistogramToolbar: public HistogramToolbarBase
{
public:
  HistogramToolbar(QWidget *parent): HistogramToolbarBase(parent) {}

private:
};

#endif
