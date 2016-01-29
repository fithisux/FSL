/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "orthowidget.h"
#include <qlayout.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <qwidget.h>
#include <qtimer.h>
#include <qpixmap.h>
#include <qimage.h>
#include <qspinbox.h>
#include <qstatusbar.h>
#include <qpushbutton.h>
#include <qtoolbar.h>
#include <qfiledialog.h>
#include <qtoolbutton.h>
#include <qlayout.h>
#include "tracker.h"
#include "ortho.xpm"
#include "overlaywidget.h"
#include "maintoolbar.h"
#include "modetoolbar.h"

#include <iostream>
using namespace std;

OrthoWidget::OrthoWidget(QWidget *parent, ImageGroup::Handle& i,
                         OverlayList::Handle ol, Cursor::Handle& c): 
  ImageWidget(parent, i, ol, c), m_image(i), m_layout(Traditional)
{
  TRACKER("OrthoWidget::OrthoWidget(QWidget *, Cursor*, ImageGroup::Handle&)");

  m_centralWidget = new QWidget(this); 

  m_grid = new QGridLayout(m_centralWidget, 2, 2, 6, -1);
  m_centralWidget->setBackgroundColor(QColor(128, 128, 128));

  setIcon( QPixmap(ortho_xpm) );

  //m_overlayList->attach(this);

  //m_slices = SliceListHandle(new SliceList);
 

//   SliceWidget::Handle m_coronal = 
//     SliceWidget::Handle(new CoronalWidget(m_centralWidget, "coronal", m_cursor,
// 					  m_overlayList, m_drawSettings, m_undoList));

//   SliceWidget::Handle m_axial = 
//     SliceWidget::Handle(new AxialWidget(m_centralWidget, "axial", m_cursor, 
// 					m_overlayList, m_drawSettings, m_undoList));

//   SliceWidget::Handle m_sagittal = 
//     SliceWidget::Handle(new SagittalWidget(m_centralWidget, "sagittal", m_cursor, 
// 					   m_overlayList, m_drawSettings, m_undoList));

  m_coronal  = new SliceView(m_centralWidget, "corornal");
  m_axial    = new SliceView(m_centralWidget, "axial");
  m_sagittal = new SliceView(m_centralWidget, "sagittal");

  SliceWidget* coronal =  new CoronalWidget(m_coronal, "coronal", m_cursor,
					    m_overlayList, m_drawSettings, m_undoList, m_opts);
  
  SliceWidget* axial =    new AxialWidget(m_axial, "axial", m_cursor, 
					  m_overlayList, m_drawSettings, m_undoList, m_opts);
  
  SliceWidget* sagittal = new SagittalWidget(m_sagittal, "sagittal", m_cursor, 
					     m_overlayList, m_drawSettings, m_undoList, m_opts);
  m_coronal->setSliceWidget(coronal);
  m_axial->setSliceWidget(axial);
  m_sagittal->setSliceWidget(sagittal);

  ImageInfo::Handle info(m_image->getMainImage()->getInfo());

  setLayout();

  connect(this,           SIGNAL(modeChanged(SliceWidget::Mode)), 
	  coronal,  SLOT(setMode(SliceWidget::Mode)));
  connect(this,           SIGNAL(modeChanged(SliceWidget::Mode)), 
	  sagittal, SLOT(setMode(SliceWidget::Mode)));
  connect(this,           SIGNAL(modeChanged(SliceWidget::Mode)),
	  axial,    SLOT(setMode(SliceWidget::Mode)));

  connect(this, SIGNAL(zoomValueChanged(int)), coronal, SLOT(setZoom(int)));
  connect(this, SIGNAL(zoomValueChanged(int)), sagittal,SLOT(setZoom(int)));
  connect(this, SIGNAL(zoomValueChanged(int)), axial,   SLOT(setZoom(int)));
 
  connect(this, SIGNAL(resetZoom()), coronal, SLOT(resetZoom()));
  connect(this, SIGNAL(resetZoom()), sagittal,SLOT(resetZoom()));
  connect(this, SIGNAL(resetZoom()), axial,   SLOT(resetZoom()));
  
  connect(coronal,  SIGNAL(message(const QString&, int )), SIGNAL(message(const QString&, int )));
  connect(sagittal, SIGNAL(message(const QString&, int )), SIGNAL(message(const QString&, int )));
  connect(axial,    SIGNAL(message(const QString&, int )), SIGNAL(message(const QString&, int )));

  connect(this, SIGNAL(crossHairModeChanged(bool)), coronal,  SLOT(crossHairMode(bool)));
  connect(this, SIGNAL(crossHairModeChanged(bool)), sagittal, SLOT(crossHairMode(bool)));
  connect(this, SIGNAL(crossHairModeChanged(bool)), axial,    SLOT(crossHairMode(bool)));

  connect(coronal, SIGNAL(emitZoomFactor(int)), SLOT(setZoomValue(int)));
  connect(sagittal,SIGNAL(emitZoomFactor(int)), SLOT(setZoomValue(int)));
  connect(axial,   SIGNAL(emitZoomFactor(int)), SLOT(setZoomValue(int)));

  m_modeWidget->enableSwitchViews(true);
  m_modeWidget->setSwitchHelpText("Switch View<hr>Toggles between different layouts.");
  connect(m_modeWidget, SIGNAL(switchViewsClicked()), SLOT(changeView()));

  //m_cursor->attach(this);

  //std::for_each(m_slices->begin(), m_slices->end(), SetImageCursor(m_cursor));

  setLabels(m_overlayList.get());

  setCentralWidget(m_centralWidget);
}

OrthoWidget::~OrthoWidget()
{
  TRACKER("OrthoWidget::~OrthoWidget");  
  m_cursor->detach(this);
}

std::string axisCodeToString(int code, bool lower)
{
  std::string s;
  switch(code)
    {
    case NIFTI_L2R: lower ? s = "L" : s = "R"; break;
    case NIFTI_R2L: lower ? s = "R" : s = "L"; break;
    case NIFTI_P2A: lower ? s = "P" : s = "A"; break;
    case NIFTI_A2P: lower ? s = "A" : s = "P"; break;
    case NIFTI_I2S: lower ? s = "I" : s = "S"; break;
    case NIFTI_S2I: lower ? s = "S" : s = "I"; break;
    default: s = ""; break;
    }

  return s;
}

void OrthoWidget::print()
{
  QString fn = QFileDialog::getSaveFileName("screenshot.png", 
					    "PNG files (*.png)", this,
					    "Screenshot dialog",
					    "Select a filename for saving");
  if(!fn.isNull()) 
    {
//       QPixmap axial(m_axial->getPixmap());
//       QPixmap coronal(m_coronal->getPixmap());
//       QPixmap sagittal(m_sagittal->getPixmap());
  
//       int width  = sagittal.width() + coronal.width();
//       int height = axial.height() + coronal.height();
      
//       QPixmap composite(width, height );
//       composite.fill(QColor(128,128,128));

//       int ax(0), ay(coronal.height()), cx(0), cy(0), sx(axial.width()), sy(0);

//       copyBlt(&composite, ax, ay, &axial, 0, 0);
//       copyBlt(&composite, cx, cy, &coronal, 0, 0);
//       copyBlt(&composite, sx, sy, &sagittal, 0, 0);

      QPixmap pm(m_centralWidget->size());
      bitBlt(&pm, 0, 0, m_centralWidget);

//       QImage im = pm.convertToImage();
//       int dpm( (72.0 / 2.54) * 100.0 );
//       im.setDotsPerMeterX(dpm);
//       im.setDotsPerMeterY(dpm);
      pm.save(fn, "PNG", 100);
    }
}

void OrthoWidget::update(const OverlayList* o, OverlayListMsg msg)
{
  ImageWidget::update(o, msg);

  setLabels(o);
}

void OrthoWidget::setLayout()
{
  if(m_grid)
    delete m_grid;
  switch(m_layout)
    {
    case Traditional:
      m_grid = new QGridLayout(m_centralWidget, 2, 2, 6, -1);
      m_grid->addWidget(m_coronal,  0, 0);
      m_grid->addWidget(m_sagittal, 0, 1);
      m_grid->addWidget(m_axial,    1, 0);
      break;
    case InRow:
      m_grid = new QGridLayout(m_centralWidget, 1, 3, 6, -1);
      m_grid->addWidget(m_sagittal,  0, 0);
      m_grid->addWidget(m_coronal, 0, 1);
      m_grid->addWidget(m_axial,    0, 2);
      break;
    case InColumn:
      m_grid = new QGridLayout(m_centralWidget, 3, 1, 6, -1);
      m_grid->addWidget(m_sagittal,  0, 0);
      m_grid->addWidget(m_coronal, 1, 0);
      m_grid->addWidget(m_axial,    2, 0);
      break;
    default:
      break;
    }
  m_grid->activate();
}

void OrthoWidget::changeView()
{
  m_layout = Layout( (m_layout+1) % (InColumn+1) );
  setLayout();
}

void OrthoWidget::setLabels(const OverlayList* o)
{
  // coronal=ik
  // axial=ij
  // sagittal=jk

  int icode(0), jcode(0), kcode(0);

  ImageInfo::Handle i(o->getActiveMetaImage()->getImage()->getInfo());
  
  i->inqAxisOrientations(icode, jcode, kcode);

  m_sagittal->setWestText(axisCodeToString(jcode, true));
  m_sagittal->setEastText(axisCodeToString(jcode, false));
  m_sagittal->setNorthText(axisCodeToString(kcode, false));
  m_sagittal->setSouthText(axisCodeToString(kcode, true));

  m_axial->setWestText(axisCodeToString(icode, i->isStoredRadiological()));
  m_axial->setEastText(axisCodeToString(icode, !i->isStoredRadiological()));
  m_axial->setNorthText(axisCodeToString(jcode, false));
  m_axial->setSouthText(axisCodeToString(jcode, true));

  m_coronal->setWestText(axisCodeToString(icode, i->isStoredRadiological()));
  m_coronal->setEastText(axisCodeToString(icode, !i->isStoredRadiological()));
  m_coronal->setNorthText(axisCodeToString(kcode, false));
  m_coronal->setSouthText(axisCodeToString(kcode, true));

  if(m_opts.inqShowLabels()) {
    if(i->hasValidXfms() ) {
      m_sagittal->setLabelsState(SliceView::Enabled);
      m_coronal->setLabelsState(SliceView::Enabled);
      m_axial->setLabelsState(SliceView::Enabled);
    } else {
      m_sagittal->setLabelsState(SliceView::Disabled);
      m_coronal->setLabelsState(SliceView::Disabled);
      m_axial->setLabelsState(SliceView::Disabled);
    }
  } else {
    m_sagittal->setLabelsState(SliceView::Disabled);
    m_coronal->setLabelsState(SliceView::Disabled);
    m_axial->setLabelsState(SliceView::Disabled);
  }
}
