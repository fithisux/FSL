/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(VIEWWIDGET_H)
#define VIEWWIDGET_H

#include <boost/shared_ptr.hpp>
#include <qmainwindow.h>

#include "imagegroup.h"
#include "cursor.h"
#include "overlaylist.h"

class QGroupBox;
class OverlayList;

/** 
 * An abstract base class for all widgets which can display
 * image views such as ortho projections or lightboxes.
 */
class ViewWidget : public QMainWindow, public CursorObserver
{
  Q_OBJECT
		
public:
  typedef boost::shared_ptr< ViewWidget > Handle;

  ViewWidget(QWidget *parent);
  virtual ~ViewWidget();
  virtual OverlayList::Handle getOverlayList(){OverlayList::Handle h;return h;}
  virtual void update(const Cursor::Handle& c){}

signals:
  void message(const QString&, int);
  void imageCursorChanged(int, int, int);
  void volumeChanged(int);
  void modeChanged(int);
  void addLookUpTable();
  void windowClose(QCloseEvent*);
  void overlayEvent();

public slots:
  virtual void setImageCursor(int, int, int);

protected:
  virtual void closeEvent(QCloseEvent*);
};

#endif 
