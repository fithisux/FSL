/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(LIGHTBOXWIDGET_H)
#define LIGHTBOXWIDGET_H

#include "imagewidget.h"
#include "slicewidget.h"

class OverlayWidget;
class QScrollView;

class LightboxWidget : public ImageWidget  
{
  Q_OBJECT
public:
  LightboxWidget(QWidget *parent,ImageGroup::Handle i,
                 OverlayList::Handle ol, Cursor::Handle& c);
  virtual ~LightboxWidget();

  //  virtual void update(const Cursor::Handle& c);

  virtual void resizeEvent(QResizeEvent*);

signals:

  void  volChanged(int); 

public slots:
  void scrolled(int);
  void repaintSlices();
  void setZoom(int);
  void print();

private:
  void layoutSlices() const;

  SliceListHandle    m_slices;
  ImageGroup::Handle m_image;
  QScrollView       *m_sv;
  QToolButton       *m_cursorModeButton;
  float              m_zoom;
};

#endif
