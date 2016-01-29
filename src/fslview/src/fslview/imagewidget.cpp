/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "imagewidget.h"
#include "slicewidget.h"
#include <qwidget.h>
#include <qcheckbox.h>
#include <qtoolbar.h>
#include <qdockwindow.h>
#include <qspinbox.h>
#include <qslider.h>
#include <qtooltip.h>
#include <qlayout.h>
#include <qpixmap.h>
#include <qtoolbutton.h>
#include <qbuttongroup.h>
#include <qmessagebox.h>
#include <qstatusbar.h>

#include "maintoolbar.h"
#include "modetoolbar.h"
#include "overlaylist.h"
#include "overlaywidget.h"
#include "overlayinfodialog.h"
#include "briconwidget.h"
#include "drawwidget.h"
#include "cursorwidget.h"
#include "talairachwidget.h"
#include "viewoptionsdialog.h"

//#define DEBUGGING
#include "tracker.h"

#include <iostream>
using namespace std;

ImageWidget::ImageWidget(QWidget* parent,        ImageGroup::Handle i,
                         OverlayList::Handle ol, Cursor::Handle c): 
  ViewWidget(parent), m_imageGroup(i), m_cursor(c->clone()), m_globalCursor(c), 
  m_overlayList(ol)
{
  setFocusPolicy(ClickFocus);

  font().setPointSize(8);

  if(!m_overlayList)
      m_overlayList = OverlayList::create(m_imageGroup);
  m_drawSettings = DrawSettings::create();
  
  m_overlayDialog = NULL;
 
  m_cursor->attach(this);
  m_globalCursor->attach(this);

  constructToolBar();

  m_overlayList->attach(this);

  connect(parent, SIGNAL(windowActivated(QWidget*)), this, SLOT(windowActivated(QWidget*)));

  m_mainToolbarWidget->setCursorMode();
}

ImageWidget::~ImageWidget()
{
  TRACKER("ImageWidget::~ImageWidget()");
  m_overlayWidget->close();
  m_overlayList->detach(this); 
  m_cursor->detach(this);
  m_globalCursor->detach(this);
  m_movieTimer->stop();
}

void ImageWidget::options()
{
  TRACKER("ImageWidget::options()");
  ViewOptionsDialog optionsDialog(this, m_opts);

  optionsDialog.m_showLabels->setEnabled(true);

  if(optionsDialog.exec() == QDialog::Accepted)
    {
       m_opts = optionsDialog.getOptions();
       setLabels(m_overlayList.get());
       setMovieFrameRate(m_opts.inqMovieFrameRate());
       short volume( m_overlayList->getActiveMetaImage()->getDs()->inqCurrentVolume() );
       MESSAGE(QString("volume = %1").arg(volume));
       setVolumeValue(volume);
       m_cursor->setVolume(volume);
       m_cursor->repaint();
    }
}

void ImageWidget::windowActivated(QWidget *w) 
{
  TRACKER("ImageWidget::windowActivated(QWidget *w)");
  if(w == this) {
    setFocusPolicy(QWidget::StrongFocus);
    MESSAGE("setFocusPolicy(QWidget::StrongFocus)");
  } else {
    setFocusPolicy(QWidget::NoFocus);
    MESSAGE("setFocusPolicy(QWidget::NoFocus)");
  }
  emit overlayEvent();
}

QSize ImageWidget::sizeHint() const
{
  return QSize(600, 600);
}

void ImageWidget::constructToolBar()
{
  m_toolbar = new QToolBar(this, "Display mode tools");
  m_modebar = new QToolBar(this, "Advanced mode tools");
  QDockWindow *cursordock = new QDockWindow(QDockWindow::InDock, this);
  cursordock->setCloseMode(QDockWindow::Always);
  cursordock->setResizeEnabled(true);
  cursordock->setCaption( tr("Cursor tools") );
  setAppropriate(cursordock, true);

  QDockWindow *talairachdock = new QDockWindow(QDockWindow::InDock, this);
  talairachdock->setCloseMode(QDockWindow::Always);
  talairachdock->setResizeEnabled(true);
  talairachdock->setCaption( tr("Atlas tools") );
  setAppropriate(talairachdock, true);

  m_briconToolbar = new QToolBar(this, "Brightness/contrast controls");

  m_overlayDock = new QDockWindow(QDockWindow::InDock, this);
  m_overlayDock->setCloseMode(QDockWindow::Always);
  m_overlayDock->setResizeEnabled(true);
  m_overlayDock->setCaption( tr("Overlay settings") );
  setAppropriate(m_overlayDock, true);

  m_drawToolbar = new QToolBar(this, "Draw toolbar");

  ImageInfo::Handle info(m_imageGroup->getMainImage()->getInfo());

  m_modeWidget     = new ModeToolBarWidget(m_modebar);
  m_briconWidget   = new BriConWidget(m_briconToolbar,m_overlayList);
  m_overlayWidget  = new OverlayWidget(m_overlayDock, m_overlayList); 
  m_drawWidget     = new DrawWidget(m_drawToolbar, m_overlayList, m_drawSettings); 
  m_drawWidget->setEnabled(false);
  m_cursorWidget   = new CursorWidget(cursordock, m_cursor, m_overlayList);
  m_talairachWidget= new TalairachWidget(talairachdock, m_cursor, m_overlayList);

  m_overlayDock->setWidget(m_overlayWidget);
  cursordock->setWidget(m_cursorWidget);
  talairachdock->setWidget(m_talairachWidget);

  connect(m_overlayWidget,SIGNAL(infoButtonAction()),
          this,           SLOT(openOverlayDialog()));

  m_mainToolbarWidget = new MainToolBarWidget(m_toolbar, 0, info->inqNumVolumes() - 1);
  m_mainToolbarWidget->enableMaskMode(false);

  connect(m_mainToolbarWidget,  SIGNAL(zoomValueChanged(int)),      SIGNAL(zoomValueChanged(int)));
  connect(m_mainToolbarWidget,  SIGNAL(crossHairStateChanged(int)), SLOT(crossHairModeChanged(int)));
  connect(m_mainToolbarWidget,  SIGNAL(modeChanged(SliceWidget::Mode)), 
	  SLOT(changeMode(SliceWidget::Mode)));
  connect(m_mainToolbarWidget, SIGNAL(resetZoomClicked()), SIGNAL(resetZoom()));
  connect(m_modeWidget, SIGNAL(movieStateChanged(int)), SLOT(toggleMovie(int)));
  connect(m_modeWidget, SIGNAL(optionsClicked()), SLOT(options()));
  connect(m_modeWidget, SIGNAL(printClicked()), SLOT(print()));
  connect(m_drawWidget, SIGNAL(undoButtonClicked()), SLOT(undoGraphics()));
  connect(m_drawWidget, SIGNAL(redoButtonClicked()), SLOT(redoGraphics()));

  connect(m_cursorWidget,  SIGNAL(volumeValueChanged(int)),    SLOT(setVolumeValue(int)));

  addDockWindow(m_toolbar,       tr("Main mode tools"), Top, FALSE);
  addDockWindow(m_modebar,       tr("View toolbar"), Top, FALSE);
  addDockWindow(m_briconToolbar, tr("Brightness/contrast tools"), Top, FALSE);  
  addDockWindow(m_drawToolbar,   tr("Pen/drawing palette"), Top, FALSE);
  addDockWindow(cursordock,      tr("Cursor tools"), Bottom, FALSE);
  addDockWindow(talairachdock,   tr("Talairach tools"), Bottom, FALSE);
  addDockWindow(m_overlayDock,   tr("Overlay settings"), Bottom, FALSE);

  m_drawToolbar->hide();
  talairachdock->hide();

  m_movieTimer = new QTimer(this);
  connect( m_movieTimer, SIGNAL(timeout()), SLOT(nextFrame()) );
   
  m_modeWidget->enableMovieMode(info->inqNumVolumes() > 1);
}

void ImageWidget::changeMode(SliceWidget::Mode m)
{
  switch(m) {
  case SliceWidget::Masking:
    m_drawToolbar->show();
    m_drawWidget->setEnabled(true);
    m_drawWidget->updateControls();
    break;
  default:
    m_drawToolbar->hide();
    m_drawWidget->setEnabled(false);
    break;
  }

  emit modeChanged(m);
}

void ImageWidget::crossHairModeChanged(int state) 
{ 
  emit crossHairModeChanged(state == QButton::On); 
}

struct SetVolume
{
  SetVolume(int v): m_vol(v) {}

  void operator()(const MetaImage::Handle& mi)
  {
    mi->getDs()->setCurrentVolume(m_vol);
  }

  int m_vol;
};

void ImageWidget::setVolumeValue(int n)
{
  TRACKER("ImageWidget::setVolumeValue()");

  if( m_opts.inqVolumeIndexingWithinView() )
    for_each( m_overlayList->begin(), m_overlayList->end(), SetVolume(n) );

  m_cursor->setCursor(m_cursor->inqX(),m_cursor->inqY(),m_cursor->inqZ(), n);
}

void ImageWidget::setZoomValue(int f)
{
   TRACKER("ImageWidget::setZoomValue");
   m_mainToolbarWidget->setZoomValue(f);
}

void ImageWidget::openOverlayDialog()
{
  TRACKER("ImageWidget::openOverlayDialog");
  if(m_overlayDialog)delete m_overlayDialog;
  m_overlayDialog = new OverlayInfoDialog(this,m_overlayList,m_imageGroup);
  m_overlayDialog->show();
  connect(m_overlayDialog, SIGNAL(message(const QString&, int)), SIGNAL(message(const QString&, int)));
}

void ImageWidget::update(const Cursor::Handle& c)
{
  TRACKER("ImageWidget::update(const Cursor::Handle& c)");
  MESSAGE( QString("c = %1, %2, %3, %4").arg(c->inqX()).arg(c->inqY()).arg(c->inqZ()).arg(c->inqV()) );

  // Important - cross talk between cursor events should be co-ordinated here.
  //             If you can't explain behaviour by looking at this method
  //             then you probably are the proud owner of a nasty bug.
  MetaImage::Handle mi = m_overlayList->getActiveMetaImage();

  bool local(c == m_cursor);

  // 0 0 0 0
  if( local && m_opts.inqUseSharedVolume() && m_opts.inqUseSharedLocation() && m_opts.inqVolumeIndexingWithinView() ) {
    for_each( m_overlayList->begin(), m_overlayList->end(), SetVolume(c->inqV()) );
    m_globalCursor->detach(this);
    m_globalCursor->setCursor(c);
    m_globalCursor->attach(this);    
  }
  // 0 0 0 1
  if( local && m_opts.inqUseSharedVolume() && m_opts.inqUseSharedLocation() && !m_opts.inqVolumeIndexingWithinView() ) {
    mi->getDs()->setCurrentVolume(c->inqV());
    m_globalCursor->detach(this);
    m_globalCursor->setCursor(c);
    m_globalCursor->attach(this);   
  }
  // 0 0 1 0
  if( local && m_opts.inqUseSharedVolume() && !m_opts.inqUseSharedLocation() && m_opts.inqVolumeIndexingWithinView() ) {
    for_each( m_overlayList->begin(), m_overlayList->end(), SetVolume(c->inqV()) );
    m_globalCursor->detach(this);
    m_globalCursor->setVolume(c->inqV());
    m_globalCursor->attach(this);   
  }
  // 0 0 1 1
  if( local && m_opts.inqUseSharedVolume() && !m_opts.inqUseSharedLocation() && !m_opts.inqVolumeIndexingWithinView() ) {
    mi->getDs()->setCurrentVolume(c->inqV());
    m_globalCursor->detach(this);
    m_globalCursor->setVolume(c->inqV());
    m_globalCursor->attach(this);   
  }
  // 0 1 0 0
  if( local && !m_opts.inqUseSharedVolume() && m_opts.inqUseSharedLocation() && m_opts.inqVolumeIndexingWithinView() ) {
    for_each( m_overlayList->begin(), m_overlayList->end(), SetVolume(mi->getDs()->inqCurrentVolume()) );
    m_globalCursor->detach(this);
    m_globalCursor->setCursor(c->inqX(), c->inqY(), c->inqZ());
    m_globalCursor->attach(this);    
  }
  // 0 1 0 1
  if( local && !m_opts.inqUseSharedVolume() && m_opts.inqUseSharedLocation() && !m_opts.inqVolumeIndexingWithinView() ) {
    m_globalCursor->detach(this);
    m_globalCursor->setCursor(c->inqX(), c->inqY(), c->inqZ());
    m_globalCursor->attach(this);    
  }
  // 0 1 1 0
  if( local && !m_opts.inqUseSharedVolume() && !m_opts.inqUseSharedLocation() && m_opts.inqVolumeIndexingWithinView() ) {
    for_each( m_overlayList->begin(), m_overlayList->end(), SetVolume(mi->getDs()->inqCurrentVolume()) );
  }
  // 0 1 1 1
  if( local && !m_opts.inqUseSharedVolume() && !m_opts.inqUseSharedLocation() && m_opts.inqVolumeIndexingWithinView() ) {
    mi->getDs()->setCurrentVolume(c->inqV());
  }
  // 1 0 0 0
  if( !local && m_opts.inqUseSharedVolume() && m_opts.inqUseSharedLocation() && m_opts.inqVolumeIndexingWithinView() ) {
    for_each( m_overlayList->begin(), m_overlayList->end(), SetVolume(c->inqV()) );
    m_cursor->detach(this);
    m_cursor->setCursor(c);
    m_cursor->attach(this);    
  }
  // 1 0 0 1
  if( !local && m_opts.inqUseSharedVolume() && m_opts.inqUseSharedLocation() && !m_opts.inqVolumeIndexingWithinView() ) {
    mi->getDs()->setCurrentVolume(c->inqV());
    m_cursor->detach(this);
    m_cursor->setCursor(c);
    m_cursor->attach(this);    
  }
  // 1 0 1 0
  if( !local && m_opts.inqUseSharedVolume() && !m_opts.inqUseSharedLocation() && m_opts.inqVolumeIndexingWithinView() ) {
    for_each( m_overlayList->begin(), m_overlayList->end(), SetVolume(c->inqV()) );
    m_cursor->detach(this);
    m_cursor->setVolume(c->inqV());
    m_cursor->attach(this);
  }
  // 1 0 1 1
  if( !local && m_opts.inqUseSharedVolume() && !m_opts.inqUseSharedLocation() && !m_opts.inqVolumeIndexingWithinView() ) {
    mi->getDs()->setCurrentVolume(c->inqV());
    m_cursor->detach(this);
    m_cursor->setVolume(c->inqV());
    m_cursor->attach(this);
  }
  // 1 1 0 0
  if( !local && !m_opts.inqUseSharedVolume() && m_opts.inqUseSharedLocation() && m_opts.inqVolumeIndexingWithinView() ) {
    for_each( m_overlayList->begin(), m_overlayList->end(), SetVolume(mi->getDs()->inqCurrentVolume()) );
    m_cursor->detach(this);
    m_cursor->setCursor(c->inqX(), c->inqY(), c->inqZ());
    m_cursor->attach(this);    
  }
  // 1 1 0 1
  if( !local && !m_opts.inqUseSharedVolume() && m_opts.inqUseSharedLocation() && !m_opts.inqVolumeIndexingWithinView() ) {
    m_cursor->detach(this);
    m_cursor->setCursor(c->inqX(), c->inqY(), c->inqZ());
    m_cursor->attach(this);    
  }
  // 1 1 1 0
  if( !local && !m_opts.inqUseSharedVolume() && !m_opts.inqUseSharedLocation() && m_opts.inqVolumeIndexingWithinView() ) {
    for_each( m_overlayList->begin(), m_overlayList->end(), SetVolume(mi->getDs()->inqCurrentVolume()) );
  }
  // 1 1 1 1
  // Do nothing
}

void ImageWidget::update(const OverlayList* i, OverlayListMsg msg)
{
  TRACKER("ImageWidget::update(const OverlayList* i, OverlayListMsg msg)");

  if(OverlayListMsg(Select) == msg)
    {
      clearUndoList();  
    }

  if(OverlayListMsg(DtiMode) == msg)
    {
      MetaImage::Handle m = i->getActiveMetaImage();
      if(m)
      {
        dtiDisplayMode(m->getDs()->inqDtiDisplay());
      }
    }
  if(OverlayListMsg(Select) == msg || OverlayListMsg(DtiMode) == msg || 
     OverlayListMsg(Visibility) == msg || OverlayListMsg(Security) == msg)
    {
   
    bool state(i->getActiveMetaImage()->inqVisibility() && !i->getActiveMetaImage()->inqReadOnly());
    m_mainToolbarWidget->enableMaskMode(state);
    m_modeWidget->enableMovieMode(i->getActiveMetaImage()->getInfo()->inqNumVolumes() > 1);
    if(!state)
      {
	m_drawWidget->setEnabled(state);
      }
    else if(m_mainToolbarWidget->inqMaskMode())
      {
        m_drawWidget->setEnabled(true);
	m_drawSettings->setMode(DrawSettings::FreeHand);
      }
    else 
      {
	m_drawWidget->setEnabled(false);
	m_mainToolbarWidget->setCursorMode();
      }
    }

  emit overlayEvent();
}  


void ImageWidget::clearUndoList()
{
  TRACKER("ImageWidget::clearUndoList()");
  m_undoList.clear();
}


void ImageWidget::undoGraphics()
{
  TRACKER("ImageWidget::undoGraphics()");
  if(!m_undoList.empty())
    {
      Shape::Handle s = m_undoList.back();
      m_undoList.pop_back();
      m_redoList.push_back(s->getBuffer());
      s->commit();
    } 
  m_cursor->repaint();
}

void ImageWidget::redoGraphics()
{
  TRACKER("ImageWidget::redoGraphics()");
  if(!m_redoList.empty())
    {
      Shape::Handle s = m_redoList.back();
      m_redoList.pop_back();
      m_undoList.push_back(s->getBuffer());
      s->commit();
    } 
  m_cursor->repaint();
}

OverlayList::Handle ImageWidget::getOverlayList()
{
 return  m_overlayList;
}

void ImageWidget::nextFrame()
{
  MetaImage::Handle mi(m_overlayList->getActiveMetaImage());
  int n(mi->getDs()->inqCurrentVolume());

  if(++n >= m_movieVols) n = 0;

  setVolumeValue(n);  
}

void ImageWidget::toggleMovie(int state)
{
  MetaImage::Handle mi(m_overlayList->getActiveMetaImage());

  m_movieVols = mi->getInfo()->inqNumVolumes();

  if(!m_movieTimer->isActive())
    m_movieTimer->start(m_opts.inqMovieFrameRate(), false);
  else
    m_movieTimer->stop();
}  

void ImageWidget::setMovieFrameRate(int ms)
{
  if(m_movieTimer->isActive())
    m_movieTimer->changeInterval(ms);
}


void ImageWidget::dtiDisplayMode(int dtiMode)
{
  if(DtiDisplay(None) == dtiMode)
    {  
      ImageInfo::Handle info(m_imageGroup->getMainImage()->getInfo());
      m_modeWidget->enableMovieMode(info->inqNumVolumes() > 1);
      m_cursorWidget->enableVolumeSpinBox(true); 
    }
  else
    { 
      m_cursorWidget->setVolumeValue(0);
      m_cursorWidget->enableVolumeSpinBox(false); 
      m_modeWidget->enableMovieMode(false); 
      m_movieTimer->stop();
    }
}

