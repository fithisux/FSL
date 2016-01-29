/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(IMAGEWIDGET_H)
#define IMAGEWIDGET_H

#include <qtimer.h>
#include "drawwidget.h"
#include "viewwidget.h"
#include "slicewidget.h"
#include "overlaylist.h"
#include "shape.h"
#include "cursor.h"
#include "viewoptions.h"

class QToolBar;
class QToolButton;
class QSpinBox;
class QCheckBox;
class OverlayWidget;
class BriConWidget;
class CursorWidget;
class TalairachWidget;
class ImageGroup;
class OverlayInfoDialog;
class QTimer;
class MainToolBarWidget;
class ModeToolBarWidget;

class ImageWidget : public ViewWidget,
                    public OverlayListObserver
{
  Q_OBJECT

public:
  ImageWidget(QWidget* parent,ImageGroup::Handle i,
              OverlayList::Handle,Cursor::Handle c);
  ~ImageWidget();  
  OverlayList::Handle getOverlayList(); 
  void update(const DrawWidget* w);
  virtual void update(const OverlayList* i, OverlayListMsg msg);
  virtual void update(const Cursor::Handle& c);
  virtual void setLabels(const OverlayList*) {}

  inline int inqX();
  inline int inqY();
  inline int inqZ();

  QSize sizeHint() const;

public slots:
  void setZoomValue(int);
  void setVolumeValue(int);
  void openOverlayDialog();
  void windowActivated(QWidget *);

  virtual void print() = 0;

private slots:
  void toggleMovie(int);
  void nextFrame();
  void crossHairModeChanged(int);
  void changeMode(SliceWidget::Mode);
  void undoGraphics();
  void redoGraphics();

protected slots:
  virtual void options();

signals:
  virtual void volumeValueChanged(int n);
  virtual void zoomValueChanged(int n);
  virtual void resetZoom();
  virtual void crossHairModeChanged(bool);
  virtual void modeChanged(SliceWidget::Mode);

protected:
  virtual void setMovieFrameRate(int);

private:
  void connectControls();
  void constructToolBar(); 
  void loadOverlaysList();
  void clearUndoList();
  void dtiDisplayMode(int);

  QSpinBox            *m_speedSpinBox;
  QTimer              *m_movieTimer;
  QToolButton         *m_movieButton;
  QToolButton         *m_maskModeButton;
  QSpinBox            *m_volSpinBox;

  ImageGroup::Handle   m_imageGroup; 

  int m_xDim;
  int m_yDim; 
  int m_zDim;

  int m_movieVols;

protected:
  QWidget             *m_centralWidget;
  QToolButton         *m_noneModeButton;
  Cursor::Handle       m_globalCursor;
  Cursor::Handle       m_cursor;
  OverlayList::Handle  m_overlayList;  
  DrawSettings::Handle m_drawSettings;
  MainToolBarWidget   *m_mainToolbarWidget;
  ModeToolBarWidget   *m_modeWidget;
  OverlayWidget       *m_overlayWidget;
  CursorWidget        *m_cursorWidget;
  TalairachWidget     *m_talairachWidget;
  BriConWidget        *m_briconWidget;
  DrawWidget          *m_drawWidget;
  OverlayInfoDialog   *m_overlayDialog;

  QToolBar            *m_toolbar; 
  QToolBar            *m_modebar; 
  QToolBar            *m_briconToolbar; 
  QToolBar            *m_cursorToolbar; 
  QToolBar            *m_drawToolbar; 
  QDockWindow         *m_overlayDock; 

  std::list<Shape::Handle>        m_undoList;
  std::list<Shape::Handle>        m_redoList;
  QCheckBox           *m_crossHairCheckBox;

  ViewOptions         m_opts;
};

int ImageWidget::inqX(){return m_xDim;}
int ImageWidget::inqY(){return m_yDim;}
int ImageWidget::inqZ(){return m_zDim;}
#endif
