/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(MODETOOLBAR_H)
#define MODETOOLBAR_H

#include <qtoolbutton.h>
#include <string>

#include "slicewidget.h"
#include "modetoolbarbase.h"

class ModeToolBarWidget : public ModeToolbarBase
{
  Q_OBJECT

public:
  ModeToolBarWidget(QWidget *parent);
  virtual ~ModeToolBarWidget();

  void enableMovieMode(bool on)     { m_movieModeButton->setEnabled(on);   }
  void enableSliceRollMode(bool on) { m_sliceRollButton->setEnabled(on);   }
  void enableSwitchViews(bool on)   { m_switchViewsButton->setEnabled(on); }
  bool inqMovieMode() { return m_movieModeButton->isOn(); }
  void setSwitchHelpText(const std::string&);

signals:
  void sliceRollStateChanged(int);
  void movieStateChanged(int);
  void switchViewsClicked();
  void optionsClicked();
  void printClicked();
};

#endif
