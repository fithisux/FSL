/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(TIMESERIESWIDGET_H)
#define TIMESERIESWIDGET_H

#include <qtoolbar.h>
#include <qstringlist.h>
#include <qcheckbox.h>
#include "viewwidget.h"
#include "singleserieswidget.h"
#include "cursor.h"
#include "qcombobox.h"
#include "qstring.h"

#include "timeserieswindowbase.h"

class QSpinBox;
class QToolButton;
class TimeSeriesPlot;
class TimeSeriesToolbar;

class TimeSeriesWidget : public TimeSeriesWindowBase
{
  Q_OBJECT
public:
  TimeSeriesWidget(QWidget *parent,
                   Image::Handle& image,
                   Cursor::Handle& cursor);  

  TimeSeriesWidget(QWidget *parent,
                   Image::Handle& image,

                   Cursor::Handle& cursor,
                   ModelFit::Handle& modelFit);

  virtual ~TimeSeriesWidget();
  
  void addFeatComboBox(QToolBar *);
  
private:

  void constructor();

  Image::Handle    m_image;
  Cursor::Handle   m_cursor;
  PlotOptions::Handle m_options;
  //int              m_viewNumber;
  int              m_contrListIndex;
  
  TimeSeriesDisplay::Handle  m_displayWidget;

public  slots:
  void closeEvent(QCloseEvent*);  

  void addPressed();
  void removePressed();
  void demeanToggled(bool);
  void percentToggled(bool);
  void modelComboActivated(int);
  void featModeToggled(bool);
  void showAxesToggled(bool);
  void printPressed();
  void intensityChanged(float,float);

signals:
  void windowClose(QCloseEvent*);
};

#endif
