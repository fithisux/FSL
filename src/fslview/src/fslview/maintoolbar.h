/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(MAINTOOLBAR_H)
#define MAINTOOLBAR_H

#include <qtoolbutton.h>
#include <qcheckbox.h>
#include <qspinbox.h>

#include "slicewidget.h"
#include "maintoolbarbase.h"

class MainToolBarWidget : public MainToolBarWidgetBase
{
  Q_OBJECT

public:
  MainToolBarWidget(QWidget *parent, int, int);
  virtual ~MainToolBarWidget();

  void enableMaskMode(bool on)      { m_maskModeButton->setEnabled(on);    }
 
  bool inqMaskMode()  { return m_maskModeButton->isOn();  }
  bool inqCrossHairState() { return m_crossHairsButton->isOn(); }

  void setCursorMode() { m_cursorModeButton->setOn(true); }
  void setPanMode()    { m_panModeButton->setOn(true);    }
  void setMaskMode()   { m_maskModeButton->setOn(true);   }
  void setZoomMode()   { m_zoomModeButton->setOn(true);   }
  void setCrossHairsMode(bool on); // { m_crossHairsButton->setOn(true); }

public slots: 
  void setZoomValue(int);
  
private slots:
  void setCursorState(int);
  void setPanState(int);
  void setZoomState(int);
  void setMaskState(int);
  
signals:
  void modeChanged(SliceWidget::Mode);
  void zoomValueChanged(int);
  void crossHairStateChanged(int);
  void resetZoomClicked();
};

#endif
