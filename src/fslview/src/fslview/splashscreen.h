/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(SPLASHSCREEN_H)
#define SPLASHSCREEN_H

#include <qframe.h>

#include "options.h"

class SplashScreen : public QFrame
{
  Q_OBJECT
public:
  SplashScreen(QWidget *parent, ApplicationOptions& opts, const char *name=0);
  virtual ~SplashScreen();
  
  void showEvent(QShowEvent *);

public slots:
  void runApplication();

private:
  ApplicationOptions& m_options;
};

#endif
