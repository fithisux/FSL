/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include <qcheckbox.h>
#include <qspinbox.h>
#include "maintoolbar.h"

MainToolBarWidget::MainToolBarWidget(QWidget *parent, int min, int max): 
  MainToolBarWidgetBase(parent)
{
  m_jumpToMaxButton->hide();
  connect(m_zoomSpinBox,       SIGNAL(valueChanged(int)), SIGNAL(zoomValueChanged(int)));
  connect(m_crossHairsButton,  SIGNAL(stateChanged(int)), SIGNAL(crossHairStateChanged(int)));

//   connect(m_cursorModeButton,  SIGNAL(stateChanged(int)), SLOT(setCursorState(int)));
//   connect(m_panModeButton,     SIGNAL(stateChanged(int)), SLOT(setPanState(int)));
//   connect(m_maskModeButton,    SIGNAL(stateChanged(int)), SLOT(setMaskState(int)));
//   connect(m_zoomModeButton,    SIGNAL(stateChanged(int)), SLOT(setZoomState(int)));
  connect(m_viewResetButton,   SIGNAL(clicked()),         SIGNAL(resetZoomClicked()));
}

MainToolBarWidget::~MainToolBarWidget()
{
}

void MainToolBarWidget::setCrossHairsMode(bool checked)
{ 
  m_crossHairsButton->blockSignals(true); 
  m_crossHairsButton->setOn(checked); 
  m_crossHairsButton->blockSignals(false); 
}

void MainToolBarWidget::setZoomValue(int f)
{
  m_zoomSpinBox->blockSignals(true);
  m_zoomSpinBox->setValue(f);
  m_zoomSpinBox->blockSignals(false);
}

void MainToolBarWidget::setPanState(int state) 
{ 
  if(state == QButton::On)
    emit modeChanged(SliceWidget::Pan);
}

void MainToolBarWidget::setCursorState(int state) 
{ 
  if(state == QButton::On)
    emit modeChanged(SliceWidget::Cursing); 
}

void MainToolBarWidget::setZoomState(int state) 
{ 
  if(state == QButton::On)
    emit modeChanged(SliceWidget::Zoom);
}

void MainToolBarWidget::setMaskState(int state) 
{
  if(state == QButton::On)
    emit modeChanged(SliceWidget::Masking); 
}

