/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "cubeserieswidget.h"
#include "singleserieswidget.h"
#include "cursor.h"
#include "storage/timeseries.h"
#include <qlayout.h>
#include <qpushbutton.h>
#include <qtoolbar.h>
#include <qspinbox.h>
#include <qtoolbutton.h>
#include <qtooltip.h>
#include <list>


CubeSeriesWidget::CubeSeriesWidget(QWidget *parent, Image::Handle& i, Cursor::Handle& c,PlotOptions::Handle& options):
  TimeSeriesDisplay(parent)
{   
  m_tabWidget = new QTabWidget(this);

  GridSeriesWidget* gridFront  = new GridSeriesWidget(this,i,c,options,1);
  GridSeriesWidget* gridMiddle = new GridSeriesWidget(this,i,c,options,0);
  GridSeriesWidget* gridBack   = new GridSeriesWidget(this,i,c,options,-1);

  m_tabWidget->addTab(gridFront,  "z+1");  
  m_tabWidget->addTab(gridMiddle, "z+0");  
  m_tabWidget->addTab(gridBack,   "z-1");  
 
 
  m_tabWidget->showPage(gridMiddle);  


  connect(this,SIGNAL(addTimeSeriesSignal()),
          gridFront,SLOT(addTimeSeries()));
  connect(this,SIGNAL(remTimeSeriesSignal()),
          gridFront,SLOT(remTimeSeries()));      
  connect(this,SIGNAL(demeanButtonToggleSignal(bool)),
          gridFront,SLOT(demeanButtonToggle(bool)));  
  connect(this,SIGNAL(setEnabledSignal(bool)),
          gridFront,SLOT(setEnabled(bool)));  
  connect(this,SIGNAL(axisDisplaySignal()),
          gridFront,SLOT(axisDisplay())); 
 

  connect(this,SIGNAL(demeanButtonToggleSignal(bool)),
          gridMiddle,SLOT(demeanButtonToggle(bool)));  
  connect(this,SIGNAL(addTimeSeriesSignal()),
          gridMiddle,SLOT(addTimeSeries()));
  connect(this,SIGNAL(remTimeSeriesSignal()),
          gridMiddle,SLOT(remTimeSeries()));
  connect(this,SIGNAL(setEnabledSignal(bool)),
          gridMiddle,SLOT(setEnabled(bool)));  
  connect(this,SIGNAL(axisDisplaySignal()),
          gridMiddle,SLOT(axisDisplay())); 


  connect(this,SIGNAL(addTimeSeriesSignal()),
          gridBack,SLOT(addTimeSeries()));
  connect(this,SIGNAL(remTimeSeriesSignal()),
          gridBack,SLOT(remTimeSeries()));      
  connect(this,SIGNAL(demeanButtonToggleSignal(bool)),
          gridBack,SLOT(demeanButtonToggle(bool)));
  connect(this,SIGNAL(setEnabledSignal(bool)),
          gridBack,SLOT(setEnabled(bool)));  
  connect(this,SIGNAL(axisDisplaySignal()),
          gridBack,SLOT(axisDisplay()));  

}

CubeSeriesWidget::~CubeSeriesWidget()
{
};



void CubeSeriesWidget::addTimeSeries()
{
  emit addTimeSeriesSignal();
}

void CubeSeriesWidget::remTimeSeries()
{
  emit remTimeSeriesSignal();
}

void CubeSeriesWidget::setEnabled(bool state)
{
  emit setEnabledSignal(state);
}

void CubeSeriesWidget::demeanButtonToggle(bool state)
{
  emit demeanButtonToggleSignal(state);
}

void CubeSeriesWidget::axisDisplay()
{
  emit axisDisplaySignal();
}

void CubeSeriesWidget::resizeEvent( QResizeEvent* )
{
    m_tabWidget->resize(this->size());
}
 
