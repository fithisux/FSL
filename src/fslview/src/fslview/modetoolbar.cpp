/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "modetoolbar.h"

#include <qtooltip.h>

using namespace std;

ModeToolBarWidget::ModeToolBarWidget(QWidget *parent): 
  ModeToolbarBase(parent)
{
  connect(m_movieModeButton,   SIGNAL(stateChanged(int)), SIGNAL(movieStateChanged(int)));

  connect(m_sliceRollButton,   SIGNAL(stateChanged(int)), SIGNAL(sliceRollStateChanged(int)));
  connect(m_switchViewsButton, SIGNAL(clicked()),         SIGNAL(switchViewsClicked()));
}

ModeToolBarWidget::~ModeToolBarWidget()
{
}

void
ModeToolBarWidget::setSwitchHelpText(const string& s)
{
  QToolTip::add(m_switchViewsButton, s);
}

