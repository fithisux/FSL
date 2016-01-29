/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "timeserieswidget.h"
#include "cursor.h"
#include "tracker.h"
#include "modelfit.h"
#include "storage/timeseries.h"
#include <qlayout.h>
#include <qpushbutton.h>
#include <qspinbox.h>
#include <qbuttongroup.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qtooltip.h>
#include <list>
#include <qlabel.h>
#include <qcheckbox.h>
#include <qprinter.h>
#include <qregexp.h>

#include "singleserieswidget.h"

#include <qwt_plot.h>

class Singleserieswidget;


TimeSeriesWidget::TimeSeriesWidget(QWidget *parent, 
                                   Image::Handle& image, 
                                   Cursor::Handle& cursor):
  TimeSeriesWindowBase(parent, 0, WDestructiveClose), m_image(image),
  m_cursor(cursor)
{  
  m_options = PlotOptions::create();

  constructor();
}

TimeSeriesWidget::TimeSeriesWidget(QWidget *parent, 
                                   Image::Handle& image, 
                                   Cursor::Handle& cursor,
                                   ModelFit::Handle& modelFit):
  TimeSeriesWindowBase(parent, 0, WDestructiveClose), m_image(image),
  m_cursor(cursor)
{  
  m_options = PlotOptions::create();
  m_options->setModelFit(modelFit); 
  m_options->setFeedBack(true);
  constructor();
}

void TimeSeriesWidget::closeEvent(QCloseEvent* e)
{
  emit windowClose(e);
}

void TimeSeriesWidget::constructor()
{  
  TRACKER("TimeSeriesWidget::constructor");  

  m_contrListIndex=0;

  layout()->remove(m_pixmap);	 // This section replaces the m_pixmap placeholder with a TSD
  delete m_pixmap;
  m_displayWidget = TimeSeriesDisplay::Handle(new SingleSeriesWidget(this,m_image,m_cursor,m_options));
  layout()->add(m_displayWidget.get());
  layout()->remove(m_intensity); // Fudge to get the positioning correct
  layout()->add(m_intensity);

  if(m_options->inqFeatMode()) {
      m_add->setEnabled(false);
      m_remove->setEnabled(false);
      m_modelCombo->setEnabled(true);
      m_modelCombo->clear();
      m_modelCombo->insertItem("No model");
      m_modelCombo->insertItem("Full model only");
      for(unsigned int i=0; i<m_options->getModelFit()->numFits(); i++)
	m_modelCombo->insertItem(m_options->getModelFit()->getConName(i));
      m_featMode->setChecked(true);
  } else {
      m_featMode->setChecked(false);
      m_featMode->setEnabled(false);
      m_modelCombo->setEnabled(false);
  }

  m_showAxes->setOn(true);

  connect(m_displayWidget.get(), SIGNAL(intensityChanged(float, float)), 
	  SLOT(intensityChanged(float, float)));
 
  m_displayWidget->show();
  m_cursor->setCursor(m_cursor->inqX(),m_cursor->inqY(),m_cursor->inqZ());
}

TimeSeriesWidget::~TimeSeriesWidget()
{
}

void TimeSeriesWidget::intensityChanged(float x, float y)
{
  m_intensity->setText(QString("Time: %1, Intensity: %2").arg(short(x)).arg(y));
}

class PrintFilter: public QwtPlotPrintFilter
{
public:
  PrintFilter() {};

  virtual QFont font(const QFont &f, Item, int) const
  {
    QFont f2 = f;
    f2.setPointSizeFloat(f.pointSize() * 0.75);
    return f2;
  }
};

void TimeSeriesWidget::printPressed()
{
  QPrinter printer;

  QString docName = m_displayWidget->title();
  if ( docName.isEmpty() )
    {
      docName.replace (QRegExp (QString::fromLatin1 ("\n")), tr (" -- "));
      printer.setDocName (docName);
    }

  printer.setCreator("fslview");
  printer.setOrientation(QPrinter::Landscape);

  if (printer.setup())
    m_displayWidget->print(printer, PrintFilter());
}

void TimeSeriesWidget::featModeToggled(bool y)
{
  m_options->setFeatMode(y);
  if(m_options->inqFeatMode())
    {
      m_add->setEnabled(false);
      m_remove->setEnabled(false);
      m_modelCombo->setEnabled(true);
      m_modelCombo->clear();
      m_modelCombo->insertItem("No model");
      m_modelCombo->insertItem("Full model only");
      for(unsigned int i=0; i<m_options->getModelFit()->numFits(); i++)
	m_modelCombo->insertItem(m_options->getModelFit()->getConName(i));
      m_modelCombo->setCurrentItem(1); // Full model only
      modelComboActivated(1);
    }
  else
    {
      m_add->setEnabled(true);
      m_remove->setEnabled(true);
    }
  m_displayWidget->remAllTimeSeries();
  m_cursor->setCursor(m_cursor->inqX(),m_cursor->inqY(),m_cursor->inqZ());
}

void TimeSeriesWidget::modelComboActivated(int item)
{
  TRACKER("TimeSeriesWidget::modelComboActivated(int curItem)");
  
  MESSAGE(QString("curItem = %1").arg(item));

  ModelFit::Handle m(m_options->getModelFit());
  m_contrListIndex=item-2;
  switch(item) {
  case 0:
    m_options->showFull(false);
    m_options->showPartial(false);
    break;
  case 1:
    m_options->showFull(true);
    m_options->showPartial(false);
    break;
  default:
    m_options->showFull(true);
    m_options->showPartial(true);
    m->curFit(m_contrListIndex);
    break;
  }
  m_cursor->setCursor(m_cursor->inqX(),m_cursor->inqY(),m_cursor->inqZ());
}

void TimeSeriesWidget::removePressed()
{
  m_displayWidget->remTimeSeries();
}

void TimeSeriesWidget::addPressed()
{
  m_displayWidget->addTimeSeries();
}

void TimeSeriesWidget::demeanToggled(bool y)
{
  m_displayWidget->setDemean(y);
  if( !y ) {
    m_displayWidget->setPercent(y);
    m_percent->setOn(false);
  }
  m_displayWidget->redraw();
}

void TimeSeriesWidget::percentToggled(bool y)
{
  m_displayWidget->setPercent(y); 
  m_displayWidget->redraw();
}

void TimeSeriesWidget::showAxesToggled(bool y)
{
  m_displayWidget->axisDisplay(y);
}
