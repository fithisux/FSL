/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(GRIDSERIESWIDGET_H)
#define GRIDSERIESWIDGET_H

#include "singleserieswidget.h"
#include "graphmanager.h"
#include "cursor.h"
#include "storage/image.h"
#include <qwidget.h>
#include <qlayout.h>
#include <list>

class GridSeriesWidget : public  TimeSeriesDisplay
{
  Q_OBJECT
public:
  GridSeriesWidget(QWidget *parent,
                   Image::Handle& i,
                   Cursor::Handle& c,
                   PlotOptions::Handle&,
                   short zOffset);
  virtual ~GridSeriesWidget();

private:
  std::list<TimeSeriesDisplay::Handle> m_plots;
  QGridLayout*        m_grid;
  

public  slots:
   void addTimeSeries();  
   void remTimeSeries();
   void demeanButtonToggle(bool);
   void setEnabled(bool);
   void axisDisplay();
signals:
   void addTimeSeriesSignal();
   void remTimeSeriesSignal();
   void demeanButtonToggleSignal(bool);
   void setEnabledSignal(bool);
   void axisDisplaySignal();
};



#endif
