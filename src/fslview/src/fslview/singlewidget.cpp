/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if defined(WIN32)
#pragma warning(disable:4786)
#endif

#include "singlewidget.h"
#include "maintoolbar.h"
#include "modetoolbar.h"

#include <qcheckbox.h>
#include <qtimer.h>
#include <qtoolbar.h>
#include <qtoolbutton.h>
#include <qlayout.h>
#include "tracker.h"

// #include "sliceroll.xpm"
// #include "view.xpm"

// static const char * sliceRollText = "Slice roll mode.<br><hr>"
// "Click to automatically move through the slices.";
// static const char * viewText = "View Button.<br><hr>"
// "Press to move through axial, coronal and sagittal views.";

SingleWidget::SingleWidget(QWidget *parent, ImageGroup::Handle i,OverlayList::Handle ol, Cursor::Handle& c) :  
  ImageWidget(parent,i,ol,c), m_image(i),m_viewNumber(SliceWidget::Axial)
{
  TRACKER("SingleWidget::SingleWidget");

  newSlice(SliceWidget::Axial,SliceWidget::Cursing);

  ImageInfo::Handle info(m_image->getMainImage()->getInfo());

  m_modeWidget->enableSliceRollMode(true);
  m_sliceRollTimer = new QTimer( this );
  connect(m_sliceRollTimer, SIGNAL(timeout()), this, SLOT(nextSlice()));
  connect(m_modeWidget, SIGNAL(sliceRollStateChanged(int)), 
	  SLOT(toggleSliceRoll(int)));

  m_modeWidget->enableSwitchViews(true);
  connect(m_modeWidget, SIGNAL(switchViewsClicked()), SLOT(changeView()));
}

SingleWidget::~SingleWidget()
{
  TRACKER("SingleWidget::~SingleWidget");
  m_sliceRollTimer->stop();

}

void SingleWidget::toggleSliceRoll(int state)
{
  if(state)
    m_sliceRollTimer->start(m_opts.inqMovieFrameRate(), false);
  else
    m_sliceRollTimer->stop();
}

void SingleWidget::setMovieFrameRate(int ms)
{
  ImageWidget::setMovieFrameRate(ms);
  if(m_sliceRollTimer->isActive())
    m_sliceRollTimer->changeInterval(ms);
}

// void SingleWidget::update(const Cursor::Handle& c)
// {
//   TRACKER("SingleWidget::update(Cursor)");
  
//   m_slice->setImageCursor(c->inqX(), c->inqY(), c->inqZ(), c->inqV());
// }

void SingleWidget::nextSlice()
{
  int n;

  switch(m_slice->inqOrient())
    {
    case SliceWidget::Axial:
      n = (m_cursor->inqZ() + 1) % m_image->inqZ(); 
      m_cursor->setCursorRepaint(m_cursor->inqX(),m_cursor->inqY(),n);
      break;
    
    case SliceWidget::Coronal:
      n = (m_cursor->inqY() + 1) % m_image->inqY(); 
      m_cursor->setCursorRepaint(m_cursor->inqX(), n ,m_cursor->inqZ());
      break;
    
    case SliceWidget::Sagittal:
      n = (m_cursor->inqX() + 1) % m_image->inqX(); 
      m_cursor->setCursorRepaint(n, m_cursor->inqY(),m_cursor->inqZ());
      break;
   }
}


void SingleWidget::changeView()
{
  if(++m_viewNumber > 2) m_viewNumber = 0;
  newSlice(m_viewNumber,m_slice->inqMode());
}


void SingleWidget::newSlice(int orient, int mode)
{

 if(orient == SliceWidget::Coronal)
   m_slice = SliceWidget::Handle(new CoronalWidget(this, "coronal", m_cursor, 
						   m_overlayList, m_drawSettings, m_undoList,
						   m_opts)); 

  if(orient == SliceWidget::Sagittal)
   m_slice = SliceWidget::Handle(new SagittalWidget(this, "sagittal", m_cursor, 
						    m_overlayList, m_drawSettings, m_undoList,
						    m_opts)); 

  if(orient == SliceWidget::Axial)
   m_slice = SliceWidget::Handle(new AxialWidget(this, "axial", m_cursor, 
                                                 m_overlayList, m_drawSettings, m_undoList,
						 m_opts)); 

  setCentralWidget(m_slice.get());
  
  connect(m_mainToolbarWidget, SIGNAL(modeChanged(SliceWidget::Mode)),
	  m_slice.get(), SLOT(setMode(SliceWidget::Mode)));
  connect(this, SIGNAL(zoomValueChanged(int)),          m_slice.get(), SLOT(setZoom(int)));
  connect(this, SIGNAL(resetZoom()),                    m_slice.get(), SLOT(resetZoom()));
  connect(this, SIGNAL(crossHairModeChanged(bool)),     m_slice.get(), SLOT(crossHairMode(bool)));  
  connect(m_slice.get(),SIGNAL(message(const QString&, int )), 
          this ,SIGNAL(message(const QString&, int )));

  m_slice->crossHairMode(m_mainToolbarWidget->inqCrossHairState());
  
  m_slice->setImageCursor(m_cursor->inqX(), m_cursor->inqY(), m_cursor->inqZ(), m_cursor->inqV());  

  //  emit modeChanged(mode);

  m_slice->show(); 
}

#include <qfiledialog.h>
#include <qpixmap.h>

void SingleWidget::print()
{
  QString fn = QFileDialog::getSaveFileName("screenshot.png", 
					    "PNG files (*.png)", this,
					    "Screenshot dialog",
					    "Select a filename for saving");
  if(!fn.isNull()) 
    {
      QPixmap pm(centralWidget()->size());
      bitBlt(&pm, 0, 0, centralWidget());

//       QImage im = pm.convertToImage();
//       int dpm( (72.0 / 2.54) * 100.0 );
//       im.setDotsPerMeterX(dpm);
//       im.setDotsPerMeterY(dpm);
      pm.save(fn, "PNG", 100);
    }

}
