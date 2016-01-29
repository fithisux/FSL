/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

//#define DEBUGGING

#include <qspinbox.h>
#include <qtoolbutton.h>
#include <qpixmap.h>
#include <qtooltip.h>
#include <qlayout.h>
#include <qlabel.h>
#include <qvbox.h>
#include <qcolor.h>
#include <qcstring.h>

#include <algorithm>
#include <stdio.h>
#include <qcanvas.h>
#include <qtimer.h>

//#include "fill.xpm"
//#include "undo.xpm"
//#include "pen.xpm"
//#include "eraser.xpm"
//#include "lutcolor.xpm"

#include "drawwidget.h"

//#define DEBUGGING
#include "tracker.h"

DrawWidget::DrawWidget(QWidget* w, OverlayList::Handle ol, DrawSettings::Handle ds):
  DrawToolbarBase(w),
  m_overlayList(ol),
  m_drawSettings(ds)
{
  connect(m_valBox,      SIGNAL(valueChanged(int)), SLOT(valBoxChanged(int)));
  connect(m_sizeBox,     SIGNAL(valueChanged(int)), SLOT(sizeBoxChanged(int)));
 
  connect(m_linkCursorButton,   
                         SIGNAL(toggled(bool)),     SLOT(linkButtonToggled(bool)));
  connect(m_penButton,   SIGNAL(toggled(bool)),     SLOT(penButtonToggled(bool)));
  connect(m_eraseButton, SIGNAL(toggled(bool)),     SLOT(eraseButtonToggled(bool))); 
  connect(m_fillButton,  SIGNAL(toggled(bool)),     SLOT(fillButtonToggled(bool)));

  connect(m_undoButton,  SIGNAL(clicked()),         SIGNAL(undoButtonClicked()));
  connect(m_redoButton,  SIGNAL(clicked()),         SIGNAL(redoButtonClicked()));

  m_undoButton->setAccel(CTRL+Key_Z);

  m_overlayList->attach(this);
  MetaImage::Handle mi = ol->getActiveMetaImage();
  if(!mi.get())
    mi = ol->getMainMetaImage();
  m_bricon = mi->getDs()->inqBriCon();
  m_bricon->attach(this);
  m_drawSettings->attach(this);

  updateControls();  
}

DrawWidget::~DrawWidget()
{
  m_bricon->detach(this);
  m_overlayList->detach(this);
  m_drawSettings->detach(this);
}

void DrawWidget::update(const DrawSettings* d)
{
  updateControls();
}

void DrawWidget::update(const OverlayList* ol, OverlayListMsg msg)
{
  TRACKER("DrawWidget::update(const OverlayList* ol, OverlayListMsg msg)");

  MESSAGE("updateControls");
  updateControls();
}

void DrawWidget::sizeBoxChanged(int value)
{
  m_drawSettings->setPenSize(value);
}

void DrawWidget::valBoxChanged(int value)
{
  m_drawSettings->setPenValue(value);
}

void DrawWidget::linkButtonToggled(bool state)
{
  m_drawSettings->setLinkCursor(state);
}

void DrawWidget::fillButtonToggled(bool state)
{
  if(state){m_drawSettings->setMode(DrawSettings::Fill);updateControls();}
}

void DrawWidget::penButtonToggled(bool state)
{
  if(state){m_drawSettings->setMode(DrawSettings::FreeHand);updateControls();}
}

void DrawWidget::eraseButtonToggled(bool state)
{
  if(state){m_drawSettings->setMode(DrawSettings::Erase);updateControls();}
}

void DrawWidget::updateControls()
{  
  TRACKER("DrawWidget::updateControls");

  blockControlSignals(true);

  m_penButton->setOn(m_drawSettings->inqMode() == DrawSettings::FreeHand);
  m_eraseButton->setOn(m_drawSettings->inqMode() == DrawSettings::Erase);
  m_fillButton->setOn(m_drawSettings->inqMode() == DrawSettings::Fill);
  m_linkCursorButton->setOn(m_drawSettings->linkCursorOn());
  m_valBox->setValue(m_drawSettings->inqPenValue());
  m_sizeBox->setValue(m_drawSettings->inqPenSize()); 
  LookUpTable::Handle lut = m_overlayList->getActiveMetaImage()->getDs()->inqLookUpTable();
  if(lut) {
    MESSAGE(QString("setting color %1").arg(m_drawSettings->inqPenValue()));
    BriCon::Handle briCon = m_overlayList->getActiveMetaImage()->getDs()->inqBriCon();

    QPixmap colorPixmap = QPixmap( *(m_colorPixmapLabel->pixmap()) );
    if(lut->isIndexLut()) {
      const LutElement& penColor = lut->inqValueIndex(m_drawSettings->inqPenValue());
      MESSAGE(QString("...to: %1 %2 %3").arg(penColor.red()).arg(penColor.green()).arg(penColor.blue()) );
      colorPixmap.fill( QColor(penColor.red(), penColor.green(), penColor.blue()) );
    } else {
      const LutElement& penColor = lut->inqValue(briCon->adjust(m_drawSettings->inqPenValue()));
      colorPixmap.fill( QColor(penColor.red(), penColor.green(), penColor.blue()) );
    }
    m_colorPixmapLabel->setPixmap(colorPixmap);
  }

  blockControlSignals(false);
}

void DrawWidget::blockControlSignals(bool state)
{
  m_valBox->blockSignals(state); 
  m_sizeBox->blockSignals(state);
  m_undoButton->blockSignals(state);
  m_penButton->blockSignals(state);
  m_eraseButton->blockSignals(state);
  m_fillButton->blockSignals(state);
}

void DrawWidget::update(const BriCon* b)
{
  if(b->inqMin()  < b->inqMax())
  {  
    if(m_valBox->value() < b->inqMin()
    || m_valBox->value() > b->inqMax())
      {
	m_valBox->setValue(int(b->inqMin()));
      }
  }
  else return; 
  
  updateControls();
}
