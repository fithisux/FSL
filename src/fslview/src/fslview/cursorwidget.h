
/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(CURSORWIDGET_H)
#define CURSORWIDGET_H

#include "cursor.h"
#include "overlaylist.h"
#include "cursorwidgetbase.h"

#include <qspinbox.h>

class QLineEdit;
class QDoubleValidator;

//! @brief User interface behaviour for displaying and controlling a Cursor object.
class CursorWidget : public CursorWidgetBase, CursorObserver, OverlayListObserver
{
  Q_OBJECT

public:
  CursorWidget(QWidget *parent,const Cursor::Handle& c,OverlayList::Handle ol);
  virtual ~CursorWidget();

  virtual void update(const Cursor::Handle& c);
  virtual void update(const OverlayList* ol, OverlayListMsg msg);
  
  inline void enableVolumeSpinBox(bool on) { m_volumeBox->setEnabled(on);     }

  void setValBoxState(bool);

signals:
  void volumeValueChanged(int);

public slots:
  void setVolumeValue(int);

private:
  Cursor::Handle m_cursor;
  OverlayList::Handle   m_overlayList;

  QDoubleValidator* m_xBoxValidator;
  QDoubleValidator* m_yBoxValidator;
  QDoubleValidator* m_zBoxValidator;
  void blockBoxSignals(bool);
  void updateValBox();
  void setInputValidators();
  
  QString fixMmBoxVal(QDoubleValidator*,QString &);
  int m_xVox,m_yVox,m_zVox;
  int m_xMm, m_yMm, m_zMm;
  bool m_valBoxState;

private slots:
  void internalVolumeValueChanged(int);
  void voxBoxChanged(int);
  void mmBoxChanged();
};

#endif
