/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(SINGLEWIDGET_H)
#define SINGLEWIDGET_H

#include "imagewidget.h"
#include "slicewidget.h"

class QTimer;

class SingleWidget : public ImageWidget  
{
  Q_OBJECT
public:
  SingleWidget(QWidget *parent, ImageGroup::Handle i, 
              OverlayList::Handle ol,
              Cursor::Handle& c);
  virtual ~SingleWidget();
  //  virtual void update(const Cursor::Handle& c);

signals:
 
  void  volChanged(int);

private slots:
  void nextSlice();
  void toggleSliceRoll(int);
  void changeView();
  void setMovieFrameRate(int);
  void print();

private:
  void newSlice(int orient, int mode);
  SliceWidget::Handle  m_slice;
  ImageGroup::Handle   m_image;
  QTimer              *m_sliceRollTimer;  
  QToolButton         *m_cursorModeButton;
  int                  m_viewNumber;
};

#endif
