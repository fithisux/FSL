/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(SINGLESERIESWIDGET_H)
#define SINGLESERIESWIDGET_H

#include <boost/shared_ptr.hpp>
#include "qwt_plot.h"
#include "curvedatalist.h"
#include "modelfit.h"
#include "storage/timeseries.h"
#include "storage/image.h"
#include <vector>
#include "cursor.h"
#include <qwidget.h>
#include "graphmanager.h"
//#include "plotoptions.h"

class PlotOptions
{
public:  
  typedef boost::shared_ptr< PlotOptions > Handle;

  static Handle create();

  virtual ~PlotOptions();

  bool inqTitle();
  bool inqXGrid();
  bool inqYGrid();
  bool inqXNums();
  bool inqYNums();
  bool inqXLabel();
  bool inqYLabel();
  int  inqXOffset();
  int  inqYOffset();
  int  inqZOffset();
  bool inqFeedback();
  bool inqFeatMode();
  bool inqAddRemEnabled();
  ModelFit::Handle & getModelFit();

  void setTitle(bool);
  void setGrids(bool x, bool y);
  void setNums(bool x, bool y);
  void setLabels(bool x, bool y);
  void setOffsets(int x,int y,int z);
  void setFeedBack(bool);
  void setFeatMode(bool);
  void setModelFit(ModelFit::Handle &);

  void showFull(bool);
  bool showFull(void) const;

  void showPartial(bool);
  bool showPartial(void) const;

private:
  PlotOptions();
  struct Implementation;
  const std::auto_ptr<Implementation> m_impl;
};

class TimeSeriesDisplay: public QwtPlot
{
  Q_OBJECT
public:  
  typedef boost::shared_ptr< TimeSeriesDisplay > Handle;
  TimeSeriesDisplay(QWidget* parent):QwtPlot(parent){};
  virtual ~TimeSeriesDisplay(){};

  virtual void redraw() = 0;

  virtual void setDemean(bool) = 0;
  virtual void setPercent(bool) = 0;
  virtual void addTimeSeries() = 0;  
  virtual void remTimeSeries() = 0;
  virtual void remAllTimeSeries() = 0;
  virtual void axisDisplay(bool) = 0;

// public  slots:
//    virtual void setEnabled(bool) = 0;
//    virtual void axisDisplay() = 0;
//    virtual void demeanButtonToggle(bool) = 0;

};


class SingleSeriesWidget: public  TimeSeriesDisplay, CursorObserver, GraphManagerObserver, ModelFitObserver
{	
  Q_OBJECT
public:

  SingleSeriesWidget(QWidget *parent, Image::Handle ,
                     Cursor::Handle,PlotOptions::Handle);

  SingleSeriesWidget(QWidget *parent, Image::Handle ,
                     Cursor::Handle,GraphManager::Handle,
                     PlotOptions::Handle);

  QSizePolicy sizePolicy();

  virtual ~SingleSeriesWidget();
  bool addTimeSeries(const TimeSeries::Handle &timeSeries, bool browse);
  void remTimeSeries(bool browse);
  bool addFeatSeries(const TimeSeries::Handle &, int index);
  

  void setLastCurveActive(bool setCursor);
  void setAllInActive();

  void setEnabled(bool);
  void axisDisplay(bool);
  void addTimeSeries();  
  void remTimeSeries();
  void remAllTimeSeries();
  void setDemean(bool state) {m_demean = state;}
  void setPercent(bool state) {m_percent = state;}
  void redraw();

  bool inqCausedCursorUpdate(){return m_causedCursorUpdate;} 
  virtual void update(const Cursor::Handle& c);
  virtual void update(const GraphManager::Handle& c);
  virtual void update(ModelFit *m);

private:
  CurveDataList::Handle m_curveDataList;  
  Image::Handle         m_image;
  Cursor::Handle        m_cursor; 
  GraphManager::Handle  m_graphManager;
  PlotOptions::Handle   m_options;
  //ModelFit              *m_modelFit;

  double m_range;
  bool m_enabled;
  bool m_percent;
  bool m_axisDisplay;
  
  bool m_demean;
  bool m_causedCursorUpdate;
  void startPlotProcess();
  void plotAllTimeSeries();
  void plotTimeSeries(CurveData::Handle);  
  void setActiveCurve(CurveData::Handle,bool setCursor);
  void constructor();
  void setGraphOptions();
  void drawMarker(int mouseXPos,int mouseYPos);
  void setCursorVolume(short vol);

private slots:
    virtual void mousePressed(const QMouseEvent &);
    virtual void mouseMoved(const QMouseEvent &);
    virtual void mouseReleased(const QMouseEvent  &);

signals:
  void intensityChanged(float, float);
};



#endif
