/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "lightboxwidget.h"
#include "maintoolbar.h"

#include <qscrollview.h>
#include <qpixmap.h>
#include <qcheckbox.h>
#include <qtoolbutton.h>
#include <qlayout.h>

#include "tracker.h"

#include "lightbox.xpm"

LightboxWidget::LightboxWidget(QWidget *parent,ImageGroup::Handle i,OverlayList::Handle ol, Cursor::Handle& c ) :  
  ImageWidget(parent,i,ol,c), m_image(i), m_zoom(1.0)
{
  TRACKER("LightboxWidget::LightboxWidget");
  m_sv = new QScrollView(this);
  
  setIcon( QPixmap(lightbox_xpm) );

  setCentralWidget(m_sv);
  
  m_sv->viewport()->setBackgroundColor(QColor(128, 128, 128));

  connect(m_sv->verticalScrollBar(), SIGNAL(valueChanged(int)), this, SLOT(scrolled(int)));
  connect(m_sv->verticalScrollBar(), SIGNAL(sliderReleased()),  this, SLOT(repaintSlices()));

  connect(m_mainToolbarWidget, SIGNAL(zoomValueChanged(int)), this, SLOT(setZoom(int)));

  m_slices = SliceListHandle(new SliceList);

  for(int n = 0; n < m_image->inqZ(); ++n)
  {
    SliceWidget::Handle axial =
      SliceWidget::Handle(new AxialWidget(m_sv->viewport(), "axial", 
					  m_cursor, m_overlayList, m_drawSettings, m_undoList,
					  m_opts));
    
    m_sv->addChild(axial.get());
    m_slices->push_back(axial);

    connect(this, SIGNAL(crossHairModeChanged(bool)),     axial.get(), SLOT(crossHairMode(bool)));
    connect(this, SIGNAL(modeChanged(SliceWidget::Mode)), axial.get(), SLOT(setMode(SliceWidget::Mode)));
    connect(this, SIGNAL(crossHairModeChanged(bool)),     axial.get(), SLOT(crossHairMode(bool)));
    connect(this, SIGNAL(resetZoom()),                    axial.get(), SLOT(resetZoom()));
    connect(axial.get(), SIGNAL(message(const QString&, int )), SIGNAL(message(const QString&, int )));

//    m_mainToolbarWidget->setCursorMode();

    axial->setSlice(n,m_cursor->inqV());
    axial->setSliceIsFixed(true);
  }
  
  m_mainToolbarWidget->setCrossHairsMode(true);
}

LightboxWidget::~LightboxWidget()
{
  TRACKER("LightboxWidget::~LightboxWidget");
}

struct SliceHider
{
  SliceHider(const QRect &b): m_boundingBox(b) {}

  void operator() (SliceWidget::Handle s)
  {
    QRect sRect = s->geometry().normalize();

    if(m_boundingBox.intersects(sRect))
      s->enableUpdates(true);
    else
      s->enableUpdates(false);

	s->show();//Added to Ensure repaint at correct place
  }

  const QRect &m_boundingBox;
};

struct SliceRepainter
{
  SliceRepainter(){}

  void operator() (SliceWidget::Handle s)
  {
	s->repaint();
  }
};

void LightboxWidget::scrolled(int v)
{
  TRACKER("LightboxWidget::scrolled");
  
  std::for_each(m_slices->begin(), m_slices->end(), SliceHider(m_sv->viewport()->geometry()));
  if(!m_sv->verticalScrollBar()->draggingSlider())repaintSlices();
}

struct SlicePlacer
{
  SlicePlacer(QScrollView* sv, unsigned int border, float zoom):
    m_sv(sv), m_borderPixels(border), m_zoom(zoom),
    m_x(border), m_y(border), 
    m_availableWidth(sv->visibleWidth()), m_heightIncrement(0),
    m_width(0) {}

  void operator()(SliceWidget::Handle s) 
  {
    unsigned int height = unsigned(s->inqHeight() * m_zoom);
    unsigned int width  = unsigned(s->inqWidth() * m_zoom);

    if((m_x + width) > m_availableWidth) {
	  m_width = std::max(m_width, m_x);
      m_x = m_borderPixels;
      m_y += m_heightIncrement + m_borderPixels;
      m_heightIncrement = 0;
    }

    int vx, vy;
    m_sv->contentsToViewport(m_x, m_y, vx, vy);
    s->hide(); //Added to Ensure repaint at correct place
    s->setGeometry(vx, vy, int(m_zoom * s->inqWidth()), int(m_zoom * s->inqHeight()));

    m_x += width + m_borderPixels;
    m_heightIncrement = std::max(m_heightIncrement, height);
  }

  unsigned int width() const { return m_width; }
  unsigned int height() const { return m_y + m_heightIncrement; }

  QScrollView* m_sv;
  unsigned int m_borderPixels;
  float        m_zoom;

  unsigned int m_x, m_y;
  unsigned int m_availableWidth;
  unsigned int m_heightIncrement;
  unsigned int m_width;
};

void LightboxWidget::layoutSlices() const
{
  TRACKER("LightboxWidget::layoutSlices");
  SlicePlacer sp = std::for_each(m_slices->begin(), m_slices->end(), SlicePlacer(m_sv, 6, m_zoom));
  m_sv->resizeContents(sp.width(), sp.height());
  std::for_each(m_slices->begin(), m_slices->end(), SliceHider(m_sv->viewport()->geometry()));
}

void LightboxWidget::setZoom(int factor)
{
  m_zoom = factor / 100.0;

  layoutSlices();
}

void LightboxWidget::resizeEvent(QResizeEvent *e)
{
  TRACKER("LightboxWidget::resizeEvent");
  layoutSlices();
}

// void LightboxWidget::update(const Cursor::Handle& c)
// {
//   TRACKER("LightboxWidget::update");
//   std::for_each(m_slices->begin(), m_slices->end(), SetImageCursor(c));  
// }

void LightboxWidget::repaintSlices()
{
#ifdef WIN32
	std::for_each(m_slices->begin(), m_slices->end(), SliceRepainter());
#endif
}

#include <qfiledialog.h>
#include <qpixmap.h>
 
void LightboxWidget::print()
{
  QString fn = QFileDialog::getSaveFileName("screenshot.png", 
					    "PNG files (*.png)", this,
					    "Screenshot dialog",
					    "Select a filename for saving");
  if(!fn.isNull()) 
    {
      QPixmap pm(centralWidget()->size());
      bitBlt(&pm, 0, 0, centralWidget());

//       QImage im = pm.convertToImage();
//       int dpm( (72.0 / 2.54) * 100.0 );
//       im.setDotsPerMeterX(dpm);
//       im.setDotsPerMeterY(dpm);
      pm.save(fn, "PNG", 100);
    }

}
