/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "gridserieswidget.h"
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


GridSeriesWidget::GridSeriesWidget(QWidget *parent, Image::Handle& image, Cursor::Handle& cursor, 
				   PlotOptions::Handle& opts, short zOffset):
  TimeSeriesDisplay(parent)
{  
  int numRows(3),numCols(3);

  m_grid = new QGridLayout(this,numRows,numCols);
  m_grid->setAutoAdd(false);
  GraphManager::Handle graphManager = GraphManager::create();

  for (int row = 0; row < numRows; ++row){
     for (int col = 0; col < numCols; ++col){

       PlotOptions::Handle options = PlotOptions::create();

       options->setNums(false,false);
       options->setTitle(false);
       options->setLabels(false,false);
       options->setOffsets(col - 1, -(row - 1),zOffset);
       options->setGrids(false,false);
       options->setFeedBack(true);

       TimeSeriesDisplay::Handle graph = 
         TimeSeriesDisplay::Handle(new SingleSeriesWidget(this,
                                                          image,
                                                          cursor,
                                                          graphManager,
                                                          options));
      m_plots.push_back(graph);
      m_grid->addWidget(graph.get(),row,col);

      connect(this,SIGNAL(addTimeSeriesSignal()),
              graph.get(),SLOT(addTimeSeries()));
      connect(this,SIGNAL(remTimeSeriesSignal()),
              graph.get(),SLOT(remTimeSeries()));      
      connect(this,SIGNAL(demeanButtonToggleSignal(bool)),
              graph.get(),SLOT(demeanButtonToggle(bool)));
      connect(this,SIGNAL(setEnabledSignal(bool)),
              graph.get(),SLOT(setEnabled(bool)));     
      connect(this,SIGNAL(axisDisplaySignal()),
              graph.get(),SLOT(axisDisplay()));
     }
  }
  
  m_grid->activate();
}

GridSeriesWidget::~GridSeriesWidget()
{
};

void GridSeriesWidget::addTimeSeries()
{
  emit addTimeSeriesSignal();
}

void GridSeriesWidget::remTimeSeries()
{
  emit remTimeSeriesSignal();
}

void GridSeriesWidget::demeanButtonToggle(bool state)
{
  emit demeanButtonToggleSignal(state);
}

void GridSeriesWidget::setEnabled(bool state)
{
  emit setEnabledSignal(state);
}

void GridSeriesWidget::axisDisplay()
{
  emit axisDisplaySignal();
}


