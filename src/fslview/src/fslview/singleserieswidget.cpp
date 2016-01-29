/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

//#define DEBUGGING
#include "tracker.h"

#include <boost/shared_ptr.hpp>
#include "singleserieswidget.h"
#include "timeserieswidget.h"

SingleSeriesWidget::SingleSeriesWidget(QWidget *parent,
                                       Image::Handle i,
                                       Cursor::Handle c,
                                       PlotOptions::Handle p): 
  TimeSeriesDisplay(parent),m_image(i),m_cursor(c),
  m_options(p),m_enabled(true),
  m_percent(false),m_axisDisplay(true),m_demean(false),
  m_causedCursorUpdate(false)
{
  constructor();
}

SingleSeriesWidget::SingleSeriesWidget(QWidget *parent,
                                       Image::Handle i,
                                       Cursor::Handle c,
                                       GraphManager::Handle g,
                                       PlotOptions::Handle p): 
  TimeSeriesDisplay(parent),m_image(i),m_cursor(c),m_graphManager(g),
  m_options(p),m_enabled(true),
  m_percent(false),m_axisDisplay(true),m_demean(false),
  m_causedCursorUpdate(false)
{
  constructor();
}

QSizePolicy SingleSeriesWidget::sizePolicy()
{
  return QSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);
}

void SingleSeriesWidget::constructor()
{ 
  TRACKER("SingleSeriesWidget::constructor");

  if(m_graphManager){m_graphManager->attach(this);}
  
  m_curveDataList = CurveDataList::create();
  
  setMinimumSize(50, 50);
  setMargin(10);

  m_cursor->attach(this);
  if(m_options->getModelFit())
    m_options->getModelFit()->attach(this);
  
  setGraphOptions();

  connect(this, SIGNAL(plotMousePressed(const QMouseEvent &)),
          SLOT(mousePressed( const QMouseEvent&)));
  connect(this, SIGNAL(plotMouseMoved(const QMouseEvent &)),
          SLOT(mouseMoved( const QMouseEvent&)));
  connect(this, SIGNAL(plotMouseReleased(const QMouseEvent &)),
          SLOT(mouseReleased( const QMouseEvent&)));
}

SingleSeriesWidget::~SingleSeriesWidget()
{  
  TRACKER("SingleSeriesWidget::~SingleSeriesWidget");
  m_cursor->detach(this);
  if(m_options->getModelFit()) m_options->getModelFit()->detach(this);
  if(m_graphManager) m_graphManager->detach(this);
}


bool SingleSeriesWidget::addTimeSeries(const TimeSeries::Handle &timeSeries,
                                       bool browse)
{  
  TRACKER("SingleSeriesWidget::addTimeSeries");
  CurveData::Handle curveData = CurveData::create(timeSeries,browse);
  return m_curveDataList->push_back(curveData);
}

void SingleSeriesWidget::setLastCurveActive(bool setCursor)
{   
  TRACKER("SingleSeriesWidget::setLastCurveActive"); 
  m_curveDataList->setAllInActive();
  setActiveCurve(m_curveDataList->back(),setCursor);
}

void SingleSeriesWidget::setAllInActive()
{  
  TRACKER("SingleSeriesWidget::setAllInActive");
  m_curveDataList->setAllInActive();
}

void SingleSeriesWidget::remTimeSeries(bool browse)
{ 
  TRACKER("SingleSeriesWidget::remTimeSeries");
  if(browse){ m_curveDataList->removeBrowse();}
  else      { m_curveDataList->removeActive();}
}

void SingleSeriesWidget::remAllTimeSeries()
{
  m_curveDataList->removeAll();
}

void SingleSeriesWidget::setEnabled(bool state)
{  
  TRACKER("SingleSeriesWidget::setEnabled");
  m_enabled = state;
}

void SingleSeriesWidget::axisDisplay(bool y)
{ 
  TRACKER("SingleSeriesWidget::axisDisplay");
  if(m_enabled)
    {
      m_axisDisplay = y;
      m_options->setGrids(m_axisDisplay, m_axisDisplay);
      m_options->setNums(m_axisDisplay, m_axisDisplay);
      setGraphOptions();
      redraw();
    }
  
}

void SingleSeriesWidget::startPlotProcess()
{  
  TRACKER("SingleSeriesWidget::startPlotProcess");
  if(m_graphManager)
    {
      m_graphManager->submitRange(m_curveDataList->inqMinCurveValue(),
                                  m_curveDataList->inqMaxCurveValue());
    }
  else
    {
      setAxisOptions(QwtPlot::yLeft, QwtAutoScale::None);
      plotAllTimeSeries();
    }
}

void SingleSeriesWidget::update(const GraphManager::Handle& gm)
{  
  TRACKER("SingleSeriesWidget::update GraphManager");
  setAxisScale(QwtPlot::yLeft,gm->inqMin(),gm->inqMax()); 
  plotAllTimeSeries();
}

void SingleSeriesWidget::plotAllTimeSeries()
{   
  TRACKER("SingleSeriesWidget::plotAllTimeSeries");
    clear();
    
    for(CurveDataList::It i = m_curveDataList->begin();
        i != m_curveDataList->end();
        i++)
      {
        plotTimeSeries(*i);
      }   
    replot();
}

void SingleSeriesWidget::plotTimeSeries(CurveData::Handle cd)
{  
  TRACKER("SingleSeriesWidget::plotTimeSeries");
    TimeSeries::Handle ts = cd->inqTimeSeries();  

    int volCount = ts->inqVolCount();
    long curve = insertCurve("Timeseries");
    cd->setCurve(curve);
    if(cd->inqIsActive()){    setCurvePen(curve, QPen(green));}
    else                 {    setCurvePen(curve, QPen(blue));}
    
    switch(cd->inqIndex())
      {
      case CurveData::Null:                                            break;
      case CurveData::FiltFunc:  setCurvePen(curve,QPen(red));         break;      
      case CurveData::Full:      setCurvePen(curve,QPen(blue));        break;
      case CurveData::Cope1:     setCurvePen(curve,QPen(green));       break;
      case CurveData::Cope2:     setCurvePen(curve,QPen(cyan));        break;      
      case CurveData::Cope3:     setCurvePen(curve,QPen(darkMagenta)); break;      
      case CurveData::Cope4:     setCurvePen(curve,QPen(darkGreen));   break;
      case CurveData::PE:        setCurvePen(curve,QPen(green));       break;
      }

    double *x = new double[volCount];
    double *y = new double[volCount];
    double mean = ts->mean();
    
    if(m_graphManager.get())
      {
        setAxisScale(QwtPlot::yLeft,
                     mean - (m_graphManager->inqRange()/2.0),
                     mean + (m_graphManager->inqRange()/2.0)); 
      }
    else
      {
        setAxisOptions(QwtPlot::yLeft, QwtAutoScale::None);
      }

    if(!m_demean)
      {
        for(int n = 0;n < volCount;++n)
          {
            *(x + n) = n;
            *(y + n) = ts->value(n);
          }
      }
    else
      {
	for(int n = 0;n < volCount;++n)
	  {
	    *(x + n) = n;
	    if (!m_percent)
	      *(y + n) = ts->value(n) - mean;
	    else
	      *(y + n) = ((ts->value(n) - mean) / mean) * 100.0;
	  } 
      }

    setCurveData(curve, x, y, volCount);

    delete [] x;
    delete [] y;
}

void SingleSeriesWidget::drawMarker(int mouseXPos, int mouseYPos)
{ 
  TRACKER("SingleSeriesWidget::drawMarker");
  CurveData::Handle activeCurve = m_curveDataList->getActiveData();
  if(isValidCurveData(activeCurve))
    {
      removeMarkers();
      long m = insertMarker();

      float x = invTransform(QwtPlot::xBottom, mouseXPos);
      float y = activeCurve->inqYValue(short(x));
      if(m_demean) {
	float mean(activeCurve->inqTimeSeries()->mean());
	if(!m_percent)
	  y = y - mean;
	else
	  y = ((y - mean) / mean) * 100.0;
      }

      emit intensityChanged(x, y);

      QwtSymbol s;
      s.setStyle(QwtSymbol::Cross);
      s.setSize(30);
      
      setMarkerSymbol(m, s);
      setMarkerXPos(m, short(x));
      setMarkerYPos(m, y);

      replot();
    }
}

void SingleSeriesWidget::mouseMoved(const QMouseEvent & e)
{   
  TRACKER("SingleSeriesWidget::mouseMoved");  
  if(m_options->inqFeedback())
    {
      drawMarker(e.x(),e.y());
      setCursorVolume(short(invTransform(QwtPlot::xBottom,e.x())));
    }
}

void SingleSeriesWidget::mouseReleased(const QMouseEvent & e)
{   
  TRACKER("SingleSeriesWidget::mouseReleased");
  if(m_options->inqFeedback())
    {
      drawMarker(e.x(),e.y());
      setCursorVolume(short(invTransform(QwtPlot::xBottom,e.x())));
    }
}

void SingleSeriesWidget::mousePressed(const QMouseEvent & e)
{  
  TRACKER("SingleSeriesWidget::mousePressed");
  if(m_options->inqFeedback())
    {
      int dist;
      long closeCurve = closestCurve(e.x(),e.y(),dist);
      if(dist < 10)
        {
          m_curveDataList->setAllInActive();
          CurveData::Handle foundCurve = m_curveDataList->getCurveData(closeCurve);

          setActiveCurve(foundCurve,true);
          startPlotProcess();
        }
      else
        {
          drawMarker(e.x(),e.y());
          setCursorVolume(short(invTransform(QwtPlot::xBottom,e.x())));
        }
    }
}

void SingleSeriesWidget::setActiveCurve(CurveData::Handle curve,bool setCursor)
{  
  TRACKER("SingleSeriesWidget::setActiveCurve");
 if(isValidCurveData(curve))
    {
      curve->setIsActive(true);
      TimeSeries::Handle ts = curve->inqTimeSeries();
      if(setCursor)
        {
          m_causedCursorUpdate = true;
          m_cursor->setCursor(ts->inqX()-m_options->inqXOffset(),
                              ts->inqY()-m_options->inqYOffset(),
                              ts->inqZ()-m_options->inqZOffset());

                  
          m_causedCursorUpdate = false;
        }
    }
}

void SingleSeriesWidget::setCursorVolume(short vol)
{   
  TRACKER("SingleSeriesWidget::setCursorVolume"); 
  CurveData::Handle active = m_curveDataList->getActiveData();      
  
  if(isValidCurveData(active))
    {
      TimeSeries::Handle ts = active->inqTimeSeries();
      m_causedCursorUpdate = true;
      m_cursor->setCursor(ts->inqX()-m_options->inqXOffset(),
                          ts->inqY()-m_options->inqYOffset(),
                          ts->inqZ()-m_options->inqZOffset(),
                          vol);
      m_causedCursorUpdate = false;
    }
}

void SingleSeriesWidget::redraw()
{  
  TRACKER("SingleSeriesWidget::redraw");
  startPlotProcess();
}

void SingleSeriesWidget::addTimeSeries()
{ 
  TRACKER("SingleSeriesWidget::addTimeSeries");
 
  if(m_enabled)
    {
      if (m_image->getInfo()->isValidCoordinate(
                  m_cursor->inqX() + m_options->inqXOffset(),
                  m_cursor->inqY() + m_options->inqYOffset(),
                  m_cursor->inqZ() + m_options->inqZOffset()))
        {

          if(!m_options->inqFeatMode())
           {
              TimeSeries::Handle timeSeries =  
                m_image->getTimeSeries(m_cursor->inqX() + m_options->inqXOffset(),
                                       m_cursor->inqY() + m_options->inqYOffset(),
                                       m_cursor->inqZ() + m_options->inqZOffset());
              addTimeSeries(timeSeries,false);
            }

          setAllInActive();
          redraw();
        }
    }
}

void SingleSeriesWidget::remTimeSeries()
{  
  TRACKER("SingleSeriesWidget::remTimeSeries");  
  if(m_enabled)
    {
      remTimeSeries(false);
      setLastCurveActive(true);
      redraw();
    }
}

void SingleSeriesWidget::setGraphOptions()
{  
  TRACKER("SingleSeriesWidget::setGraphOptions");

  std::string title =  m_image->getInfo()->inqImageName();
  if(m_options->inqTitle())
    {
      QString preTitle("Timeseries - "); 
      setTitle(preTitle + title.c_str());
    }
  else
    {
      setTitle("");
    }

  if(m_options->inqXLabel()){setAxisTitle(xBottom, "Time");}
  else                      {setAxisTitle(xBottom,"");}
  if(m_options->inqYLabel()){setAxisTitle(yLeft,   "Value");}
  else                      {setAxisTitle(yLeft,"");}
    
  setCanvasBackground(QColor(white));

  setAxisMargins(QwtPlot::yLeft, 0, 0);  

  enableAxis(xBottom,m_options->inqXNums());
  enableAxis(yLeft,m_options->inqYNums());
  enableGridX(m_options->inqXGrid());
  enableGridY(m_options->inqYGrid());

  if(!m_options->inqXNums() && !m_options->inqYNums() &&
     !m_options->inqXGrid() && !m_options->inqYGrid())
    {
      m_axisDisplay = false;
    }
  else
    {
      m_axisDisplay = true;
    }
}

void SingleSeriesWidget::update(ModelFit *m)
{
  update(m_cursor);
}

void SingleSeriesWidget::update(const Cursor::Handle& c)
{  
  TRACKER("SingleSeriesWidget::update Cursor");   
  if(m_enabled)
    {
      unsigned short x(m_cursor->inqX() + m_options->inqXOffset());
      unsigned short y(m_cursor->inqY() + m_options->inqYOffset());
      unsigned short z(m_cursor->inqZ() + m_options->inqZOffset());

      MESSAGE(QString("Location x(%1) y(%2) z(%3)").arg(x).arg(y).arg(z));

      bool validTimeSeries = m_image->getInfo()->isValidCoordinate(x, y, z);
      
      if(!inqCausedCursorUpdate())
        {
          if(validTimeSeries)
            { 
              if(m_options->inqFeatMode())
                {          
                  remTimeSeries(true);
                  remTimeSeries(true);
		  remTimeSeries(true);
		  //remTimeSeries(true);
		  //remTimeSeries(true);
 
		  ModelFit::Handle model(m_options->getModelFit());

		  TimeSeries::Handle timeSeries(m_image->getTimeSeries(x, y, z));
                  addFeatSeries(timeSeries, CurveData::FiltFunc);
		  
		  if(m_options->showFull())
		    addFeatSeries(model->fullModel(x, y, z, timeSeries->mean()), CurveData::Full);
		  
		  if(m_options->showPartial())
		    if(model->copePe())
		      addFeatSeries(model->CopeCurve(x, y, z, timeSeries->mean()), CurveData::Cope1);
		    else
		      addFeatSeries(model->peCurve(x, y, z, timeSeries->mean()), CurveData::PE);

// 		  addFeatSeries(model->perCentChange(x, y, z, timeSeries->mean()), CurveData::Cope3);
                }
              else
                {          
                  remTimeSeries(true);
                  if(addTimeSeries(m_image->getTimeSeries(x, y, z), true))
                    setLastCurveActive(false);
                }     
            }
          redraw();
        }   
    }
}

bool SingleSeriesWidget::addFeatSeries(const TimeSeries::Handle & ts, 
                                       int index)
{  
  TRACKER("SingleSeriesWidget::addFeatSeries");  
  bool browse(true);
  CurveData::Handle curveData = CurveData::create(ts,browse,index);
  return m_curveDataList->push_back(curveData);
}

struct PlotOptions::Implementation
{
  Implementation()
  {
    m_title = true;
    m_xGrid = true;
    m_yGrid = true;
    m_xNums = true;
    m_yNums = true;
    m_xLabel = true;
    m_yLabel = true;
    m_xOffset = 0;
    m_yOffset = 0;
    m_zOffset = 0;
    m_feedback = true;
    m_addRemEnabled = true;
    m_modelEnabled = true;
  };  

  bool m_title;
  bool m_xGrid;
  bool m_yGrid;
  bool m_xNums;
  bool m_yNums;
  bool m_xLabel;
  bool m_yLabel;
  int  m_xOffset;
  int  m_yOffset;
  int  m_zOffset;
  bool m_feedback;
  bool m_addRemEnabled;
  ModelFit::Handle m_modelFit;
  bool m_modelEnabled;
  bool m_showFull;
  bool m_showPartial;
};
  
bool PlotOptions::inqTitle() {return m_impl->m_title;}
bool PlotOptions::inqXGrid() {return m_impl->m_xGrid;}
bool PlotOptions::inqYGrid() {return m_impl->m_yGrid;}
bool PlotOptions::inqXNums() {return m_impl->m_xNums;}
bool PlotOptions::inqYNums() {return m_impl->m_yNums;}
bool PlotOptions::inqXLabel()  {return m_impl->m_xLabel;}
bool PlotOptions::inqYLabel()  {return m_impl->m_yLabel;}
int  PlotOptions::inqXOffset() {return m_impl->m_xOffset;}
int  PlotOptions::inqYOffset() {return m_impl->m_yOffset;}
int  PlotOptions::inqZOffset() {return m_impl->m_zOffset;}
bool PlotOptions::inqFeedback() {return m_impl->m_feedback;}
bool PlotOptions::inqFeatMode() {return m_impl->m_modelFit && m_impl->m_modelEnabled;}
bool PlotOptions::inqAddRemEnabled(){return m_impl->m_addRemEnabled;}
  
ModelFit::Handle& PlotOptions::getModelFit(){return m_impl->m_modelFit;}

void PlotOptions::showFull(bool y) {  m_impl->m_showFull=y; /*notify();*/ }
bool PlotOptions::showFull(void) const { return m_impl->m_showFull; }

void PlotOptions::showPartial(bool y) {  m_impl->m_showPartial=y; /*notify();*/ }
bool PlotOptions::showPartial(void) const { return m_impl->m_showPartial; }

void PlotOptions::setTitle(bool state)    
{
  m_impl->m_title = state;
}

void PlotOptions::setGrids(bool x, bool y)
{
  m_impl->m_xGrid = x;
  m_impl->m_yGrid = y;
}
void PlotOptions::setNums(bool x, bool y) 
{
  m_impl->m_xNums = x;
  m_impl->m_yNums = y;
}
void PlotOptions::setLabels(bool x, bool y)
{
  m_impl->m_xLabel = x;
  m_impl->m_yLabel = y;
}
void PlotOptions::setOffsets(int x,int y,int z)
{
  m_impl->m_xOffset = x;
  m_impl->m_yOffset = y;
  m_impl->m_zOffset = z;
}
void PlotOptions::setFeedBack(bool state)
{
  m_impl->m_feedback = state;
}

void PlotOptions::setFeatMode(bool y)
{
  m_impl->m_modelEnabled = y;
}


void PlotOptions::setModelFit(ModelFit::Handle& model)
{
  m_impl->m_modelFit = model;
  m_impl->m_addRemEnabled = false;
}

PlotOptions::PlotOptions():
  m_impl(new Implementation)
{ 
 
}

PlotOptions::~PlotOptions(){}

PlotOptions::Handle PlotOptions::create()
{
  Handle dst(new PlotOptions());
  return dst;
}

