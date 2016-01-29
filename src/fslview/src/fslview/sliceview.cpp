/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Rama Aravind Vorray
		James Saunders
		David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2002-2005 University of Oxford  */

/*  CCOPYRIGHT */

#include "qlabel.h"
#include "qlayout.h"

#include "sliceview.h"
#include "slicewidget.h"

SliceView::SliceView(QWidget* parent, const char *name):
  SliceViewBase(parent, name)
{
  //  setBackgroundColor(QColor(0, 0, 0));
  delete m_pixmapLabel;
}

void SliceView::setLabelsState(LabelState s)
{
  switch(s)
    {
    case Enabled:
      m_northLabel->setEnabled(true);
      m_southLabel->setEnabled(true);
      m_eastLabel->setEnabled(true);
      m_westLabel->setEnabled(true);
      break;
    case Greyed:
      m_northLabel->setEnabled(false);
      m_southLabel->setEnabled(false);
      m_eastLabel->setEnabled(false);
      m_westLabel->setEnabled(false);
      break;
    case Disabled:
      m_northLabel->setText("");
      m_southLabel->setText("");
      m_eastLabel->setText("");
      m_westLabel->setText("");
   default:
      break;
    }      
}

void SliceView::setSliceWidget(SliceWidget* slice)
{
  m_slice = slice;
  m_gridLayout->addWidget(slice, 1, 1);
  slice->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
}

QPixmap SliceView::getPixmap() const
{
  QPixmap pm(size());

  bitBlt(&pm,0,0,this);

  return pm;
}

void SliceView::setNorthText(const std::string& s)
{
  m_northLabel->setText(s);
}

void SliceView::setSouthText(const std::string& s)
{
  m_southLabel->setText(s);
}

void SliceView::setEastText(const std::string& s)
{
  m_eastLabel->setText(s);
}

void SliceView::setWestText(const std::string& s)
{
  m_westLabel->setText(s);
}
