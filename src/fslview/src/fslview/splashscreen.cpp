/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if defined(WIN32)
#pragma warning(disable:4786)
#endif

#include <qtimer.h>
#include <qpixmap.h>
#include <qapplication.h>

#include "application.h"
#include "splashscreen.h"
#include "version.h"

#if !defined(WIN32)
#include "fslstart.xpm"
#endif

/**
 * Display a splash screen for a few seconds and then invoke the app.
 */
SplashScreen::SplashScreen(QWidget *parent, ApplicationOptions& opts, const char *name) 
  : QFrame(parent, name, QWidget::WStyle_NoBorder | QWidget::WStyle_Customize), m_options(opts)
{

  #if !defined(WIN32)

  QPixmap pm( fslstart_xpm );
  setBackgroundPixmap( pm );

  int  w = pm.width()/2;
  int  h = pm.height()/2;
  int dw = QApplication::desktop()->width()/2;
  int dh = QApplication::desktop()->height()/2;

  setGeometry( (dw) - (w),
      	       (dh) - (h),
	       pm.width(), pm.height() );
  setFrameStyle( QFrame::Box | QFrame::Raised );


  show();
  #endif
}

SplashScreen::~SplashScreen()
{
}

void SplashScreen::showEvent(QShowEvent *e)
{

    QTimer *timer = new QTimer( this );
    connect( timer, SIGNAL(timeout()), this, SLOT(runApplication()) );
#if defined(WIN32)    
	timer->start( 1, TRUE );
#else
	timer->start( 1000, TRUE );
#endif
}

void SplashScreen::runApplication()
{
  this->hide();

  ApplicationWindow *w = new ApplicationWindow(m_options);
  w->setCaption(QString("FSLView (%1.%2)").arg(Version).arg(Release));
  w->show();
}
