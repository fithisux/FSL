/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */


#include "briconwidget.h"

#include <qlineedit.h>
#include <qvalidator.h>
#include <qtoolbutton.h>
#include <qlayout.h>
#include <qpixmap.h>
#include <qlabel.h>
#include <qtooltip.h>

#include <qwt_wheel.h>

#include <algorithm>
#include <stdio.h>

//#define DEBUGGING
#include "tracker.h"

// class BriConWheel: public QwtWheel
// {
// public:
//   BriConWheel(QWidget *parent, const char *name = 0): QwtWheel(parent, name)
//   {
//     setOrientation(Qt::Horizontal);
//     setTickCnt(20);
//     setViewAngle(170);
//     setTotalAngle(6 * 360);
//   }
// };

// class LimitBox: public QLineEdit
// {
//  public:
//   LimitBox(QWidget *parent): QLineEdit(parent)
//   {
//     //    setFont(QFont("Helvetica", 10));
//     setMinimumWidth( QFontMetrics(font()).width(QString("-0.000000e-00")) );
//     setMaximumWidth( QFontMetrics(font()).width(QString("-0.000000e-00")) );
//     setValidator(new QDoubleValidator(this));
//     setAlignment(AlignLeft);
//   }
// };

BriConWidget::BriConWidget(QWidget* w, OverlayList::Handle list): 
  BriConWidgetBase(w), m_list(list)//, m_blockEvents(false)
{
  MetaImage::Handle mi = m_list->getActiveMetaImage();
  if(!mi.get())
    mi = m_list->getMainMetaImage();
  m_bricon = mi->getDs()->inqBriCon();

  m_bricon->attach(this);

  m_list->attach(this);

  m_minBox->setText( QString("%1").arg(m_bricon->inqMin(),7,'g',5) );
  m_minBox->setValidator(new QDoubleValidator(m_minBox));

  m_maxBox->setText( QString("%1").arg(m_bricon->inqMax(),7,'g',5) );
  m_maxBox->setValidator(new QDoubleValidator(m_minBox));

  m_briSlider->setRange(-49, 49, 2, 1);
  m_conSlider->setRange(-49, 49, 2, 1);

  connect(m_resetBriCon, SIGNAL(pressed()),this, SLOT(reset())); 

  connect(m_minBox, SIGNAL(lostFocus()), SLOT(minChanged()));
  connect(m_maxBox, SIGNAL(lostFocus()), SLOT(maxChanged()));
  connect(m_minBox, SIGNAL(returnPressed()), SLOT(minChanged()));
  connect(m_maxBox, SIGNAL(returnPressed()), SLOT(maxChanged()));

  connect(m_briSlider, SIGNAL(valueChanged(double)),this, SLOT(briSliderChanged(double))); 
  connect(m_conSlider, SIGNAL(valueChanged(double)),this, SLOT(conSliderChanged(double)));
  connect(m_briSlider, SIGNAL(sliderReleased()), this, SLOT(sliderReleased()));
  connect(m_conSlider, SIGNAL(sliderReleased()), this, SLOT(sliderReleased()));
  
}

BriConWidget::~BriConWidget()
{
  m_list->detach(this);
  m_bricon->detach(this);
}

void BriConWidget::setMinMaxBoxesState(bool state)
{
  m_minBox->setEnabled(state);
  m_maxBox->setEnabled(state);
}

void BriConWidget::setBriSliderState(bool state)
{
  m_briSlider->blockSignals(true);
  m_briSlider->setValid(state);
  m_briSlider->blockSignals(false);
}

void BriConWidget::setConSliderState(bool state)
{
  m_conSlider->blockSignals(true);
  m_conSlider->setValid(state);
  m_conSlider->blockSignals(false);
}

void BriConWidget::reset()
{
  m_bricon->detach(this);
  m_bricon->reset();
  m_bricon->attach(this);
  sliderReleased();
}

void BriConWidget::update(const BriCon* bricon)
{
  TRACKER("BriConWidget::update(const BriCon* bricon)");
  updateMinMaxBoxes();
}

void BriConWidget::update(const OverlayList* list, OverlayListMsg msg)
{
  TRACKER("BriConWidget::update(const OverlayList* list, OverlayListMsg msg)");

  if(OverlayListMsg(Select) == msg) 
    {
      MESSAGE("Select");
      MetaImage::Handle mi = list->getActiveMetaImage();
      if(!mi.get())
        mi = list->getMainMetaImage();
      m_bricon->detach(this);
      m_bricon = mi->getDs()->inqBriCon();
      m_bricon->attach(this);
      m_minBox->setText( QString("%1").arg(m_bricon->inqMin(),7,'g',5) );
      m_maxBox->setText( QString("%1").arg(m_bricon->inqMax(),7,'g',5) );
    }
  if((OverlayListMsg(Select) == msg) || (OverlayListMsg(Visibility) == msg))
    {
      MESSAGE("Select || Visibility");
      bool state(list->getActiveMetaImage()->inqVisibility());
      setMinMaxBoxesState(state);
      setBriSliderState(state);
      setConSliderState(state);
    }
}

void BriConWidget::minChanged()
{
  TRACKER("BriConWidget::minChanged");

  float value = m_minBox->text().toFloat();
  m_bricon->setMin(value);
}

void BriConWidget::maxChanged()
{
  TRACKER("BriConWidget::maxChanged");

  float value = m_maxBox->text().toFloat();
  m_bricon->setMax(value);
}

void BriConWidget::briSliderChanged(double value)
{
  TRACKER("BriConWidget::briSliderChanged");

  m_bricon->modifyRange(value / 100.0, 0.0);
  updateMinMaxBoxes();
}

void BriConWidget::conSliderChanged(double value)
{
  TRACKER("BriConWidget::conSliderChanged");

  m_bricon->modifyRange(0.0, value / 100.0);
  updateMinMaxBoxes();
}

void BriConWidget::sliderReleased()
{
  TRACKER("BriConWidget::sliderReleased");

  m_bricon->updateRange();

  m_minBox->setText( QString("%1").arg(m_bricon->inqMin(),7,'g',5) );
  m_maxBox->setText( QString("%1").arg(m_bricon->inqMax(),7,'g',5) );
  
  m_briSlider->setValue(0);
  m_conSlider->setValue(0);
}

void BriConWidget::updateMinMaxBoxes()
{
  TRACKER("BriConWidget::updateMinMaxBoxes()");
  m_minBox->blockSignals(true);
  m_minBox->setText( QString("%1").arg(m_bricon->inqAdjustedMin(),7,'g',5) );      
  m_minBox->blockSignals(false);

  m_maxBox->blockSignals(true);
  m_maxBox->setText( QString("%1").arg(m_bricon->inqAdjustedMax(),7,'g',5) );
  m_maxBox->blockSignals(false);
}
