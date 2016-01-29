/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "viewwidget.h"

#include <qgroupbox.h>
#include <qlayout.h>

#include "tracker.h"

ViewWidget::ViewWidget(QWidget *parent)
  : QMainWindow(parent, 0, WDestructiveClose)
{
  TRACKER("ViewWidget::ViewWidget");
}

ViewWidget::~ViewWidget()
{
}

void ViewWidget::setImageCursor(int x, int y, int z)
{
  emit imageCursorChanged(x, y, z);
}

void ViewWidget::closeEvent(QCloseEvent* e)
{
  emit windowClose(e);
}

