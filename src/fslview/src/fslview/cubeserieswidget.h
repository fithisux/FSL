/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(CUBESERIESWIDGET_H)
#define CUBESERIESWIDGET_H

#include "gridserieswidget.h"
#include "cursor.h"
#include "storage/image.h"
#include <qwidget.h>
#include <qlayout.h>
#include <qtabwidget.h>
#include <list>

class CubeSeriesWidget : public  TimeSeriesDisplay
{
  Q_OBJECT
public:
  CubeSeriesWidget(QWidget *parent,Image::Handle& i,Cursor::Handle& c,PlotOptions::Handle& options);
  virtual ~CubeSeriesWidget();

private:
  QTabWidget*         m_tabWidget;
  

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
protected:
  virtual void resizeEvent( QResizeEvent* );  
 

};



#endif
