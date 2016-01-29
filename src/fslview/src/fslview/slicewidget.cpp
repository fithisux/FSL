
/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */
/*  CCOPYRIGHT */


/****************************************************************************
** $Id: slicewidget.cpp,v 1.204.2.2 2010/04/19 15:46:30 flitney Exp $
**
** Copyright (C) 2002 University of Oxford.  All rights reserved.
**
** FSLView
**
*****************************************************************************/

#if defined(WIN32)
#pragma warning (disable:4786)
#endif

#include <qbitmap.h>
#include <qlabel.h>
#include <qimage.h>
#include <algorithm>
#include <list>
#include "slicewidget.h"
#include "imagewidget.h"
#include "overlaylist.h"
#include "rect.h"
#include <assert.h>
#include <math.h>
#include "shape.h"
#include "qtimer.h"
#include "imagebuffer.h"
#include "drawsettings.h"

//#define DEBUGGING
#include "tracker.h"

#include "icons/crosscursor.xpm"
#include "icons/crossmask.xpm"
#include "icons/zoomcursor.xpm"
#include "icons/zoommask.xpm"
#include "icons/pancursor.xpm"
#include "icons/pencursor.xpm"
#include "icons/penmask.xpm"
#include "icons/panmask.xpm"
#include "icons/erasercursor.xpm"
#include "icons/erasermask.xpm"
#include "icons/fillcursor.xpm"
#include "icons/fillmask.xpm"

using namespace std;

QCursor* SliceWidget::m_crossCursor = 0;
QCursor* SliceWidget::m_zoomCursor = 0;
QCursor* SliceWidget::m_panCursor = 0;
QCursor* SliceWidget::m_penCursor = 0;
QCursor* SliceWidget::m_eraserCursor = 0;
QCursor* SliceWidget::m_fillCursor = 0;

struct Detach
{
  Detach(SliceWidget* s): m_sliceWidget(s) {}
  void operator()(const BriCon::Handle& bh) 
  {
    bh->detach(m_sliceWidget);
  }

  SliceWidget* m_sliceWidget;
};

struct Attach
{
  Attach(SliceWidget* s, std::list<BriCon::Handle>& bl): m_sliceWidget(s), m_briconList(bl) {}
  void operator()(MetaImage::Handle& mi)
  {
    BriCon::Handle bh(mi->getDs()->inqBriCon());
    bh->attach(m_sliceWidget);
    m_briconList.push_back(bh);
  }

  SliceWidget* m_sliceWidget;
  std::list<BriCon::Handle>& m_briconList;
};

//! @brief Class constructor
//!
//! @param parent The parent widget to be passed to QWidget constructor
//! @param name   The widget name to be passed to QWidget constructor
//! @param orient This slices @ref Orientation
//! @param c A Cursor which this SliceWidget will observe and track
//! @param l An OverlayList containing the images to be rendered
//! @param d The DrawSettings object of the parent object. Used to indicate
//!          mask drawing options.
//! @param u List of shapes drawn in mask mode. Used for the undo functionality
//!
SliceWidget::SliceWidget(QWidget* parent, const char *name, Orientation orient, 
			 Cursor::Handle c, OverlayList::Handle l, 
			 DrawSettings::Handle d, 
			 std::list<Shape::Handle>& u,
			 const ViewOptions& vo): 
  QWidget(parent, name),
  m_sliceIsFixed(false), m_updatesEnabled(true), m_crossHairsOn(true),
  m_imagesEnabled(true), 
  m_slice(0), m_noSliceSet(true),m_volume(0), m_cursor(c),
  m_zoom(1.0),m_origX(0),m_origY(0),m_scaleX(1.0),m_scaleY(1.0),
  m_overlayList(l), m_drawSettings(d),m_undoList(u),
  m_orient(orient),m_mode(None),m_zooming(false),m_trueScale(true),
  m_forceRender(false), 
  m_imageBuffersValid(false),m_displayPixmapValid(false),
  m_displayBufferValid(false),
  m_opts(vo)
{
  TRACKER("SliceWidget::SliceWidget");  

  m_cursor->attach(this);
  m_overlayList->attach(this);
  std::for_each(m_overlayList->begin(), m_overlayList->end(), Attach(this, m_briconList));
  m_drawSettings->attach(this);
  
  setFocusPolicy(StrongFocus);
  loadStore();

  m_prevCursor.reset();

  //  setMinimumSize(200,200);

  if(!m_crossCursor)
    {
      QPixmap pixmap = QPixmap(crosscursor_xpm);
      QBitmap mask;
      mask = QPixmap(crossmask_xpm);
      pixmap.setMask(mask);

      m_crossCursor = new QCursor(pixmap, 9, 9);
    }

  if(!m_zoomCursor)
    {
      QPixmap zoomPixmap = QPixmap(zoomcursor_xpm);
      QBitmap zoomMask;
      zoomMask = QPixmap(zoommask_xpm);
      zoomPixmap.setMask(zoomMask);

      m_zoomCursor = new QCursor(zoomPixmap, 9, 9);
    }

  if(!m_panCursor)
    {
      QPixmap panPixmap = QPixmap(pancursor_xpm);
      QBitmap panMask;
      panMask = QPixmap(panmask_xpm);
      panPixmap.setMask(panMask);

      m_panCursor = new QCursor(panPixmap, 9, 9);
    }
 
  if(!m_penCursor)
    {
      QPixmap penPixmap = QPixmap(pencursor_xpm);
      QBitmap penMask;
      penMask = QPixmap(penmask_xpm);
      penPixmap.setMask(penMask);

      m_penCursor = new QCursor(penPixmap, 1, 1);
    }
  
  if(!m_fillCursor)
    {
      QPixmap fillPixmap = QPixmap(fillcursor_xpm);
      QBitmap fillMask;
      fillMask = QPixmap(fillmask_xpm);
      fillPixmap.setMask(fillMask);

      m_fillCursor = new QCursor(fillPixmap, 2, 10);
    }
  
  if(!m_eraserCursor)
    {
      QPixmap eraserPixmap = QPixmap(erasercursor_xpm);
      QBitmap eraserMask;
      eraserMask = QPixmap(erasermask_xpm);
      eraserPixmap.setMask(eraserMask);

      m_eraserCursor = new QCursor(eraserPixmap, 5, 2);
    }
  
//   setCursor(*m_crossCursor);
  setMode(m_mode);

  setBackgroundMode(NoBackground);
   
  //Need this or lightbox will crash
  m_viewRect = Rect::createRect(0,0,64,64);
  m_zoomRect = Rect::createRect(0,0,0,0);
  m_dataRect = Rect::createRect(0,0,64,64);
   
  m_displayPixmap = new QPixmap();
  
  m_timer = new QTimer(this);
  connect( m_timer, SIGNAL(timeout()), this, SLOT(showSlice()));
  m_timer->start(50,false);
  
  emitZoomFactor(100);

  m_dtiColors.push_back(QColor(255, 0,   0));
  m_dtiColors.push_back(QColor(  0, 0, 255));
  m_dtiColors.push_back(QColor(255, 0, 255));
}

SliceWidget::~SliceWidget()
{
  TRACKER("SliceWidget::~SliceWidget");
  m_cursor->detach(this);
  m_overlayList->detach(this);
  std::for_each(m_briconList.begin(), m_briconList.end(), Detach(this));
  m_drawSettings->detach(this);
  delete m_displayPixmap;
}

void SliceWidget::polish()
{  
  TRACKER("SliceWidget::polish");
  m_viewRect = Rect::createRect(0,0,inqWidth(),inqHeight());
  m_zoomRect = Rect::createRect(0,0,0,0);
  m_dataRect = Rect::createRect(0,0,inqWidth(),inqHeight());           
  m_displayPixels = ColorRGBAHandle( new ColorRGBA[inqWidth() * inqHeight()] );
}

//! @brief 
//!
//!
void SliceWidget::resizeEvent(QResizeEvent* e)
{  
  TRACKER("SliceWidget::resizeEvent");

  initZoom();
}


void SliceWidget::initZoom()
{
  int l  = std::min(width(), height());
  m_zoom = std::min( (l/(float)inqWidth()), 
		     std::min( (l/((float)inqHeight()*(float)inqRatio())), 
			       (l/((float)inqDepth()*(float)depthRatio()))) );
}

//! @brief Turn on/off response to the update method.
//!
//! @param enabled Use the value false to disable updates.
//!
void SliceWidget::enableUpdates(bool enabled)
{
  m_updatesEnabled = enabled;
}

QPixmap SliceWidget::getPixmap() const
{
  QPixmap pm(size());

  bitBlt(&pm,0,0,this);

  return pm;
}

void SliceWidget::paintEvent(QPaintEvent*)
{
  if(!m_updatesEnabled) { m_imagesEnabled = true; return; }
  TRACKER("SliceWidget::paintEvent");

  if(!imageBuffersValid()) 
    renderBuffer();  

  CHECKPOINT();
  QPixmap pm(size());

  if(!m_imagesEnabled){bitBlt(&pm,0,0,this);} // Why!?
  else{pm.fill( black );}

  m_paint.begin(&pm,this);

  float fit;

  if(m_trueScale)
    {
      fit = 1.0;
    }
  else
    {
      fit = std::min(width()/(m_viewRect->width()*m_zoom), height()/(m_viewRect->height()*m_zoom*inqRatio()));
    }

  float zoom(m_zoom * fit);

  m_origX  = int(width() -  (m_viewRect->width()            *zoom))/2;
  m_origY  = int(height() - (m_viewRect->height()*inqRatio()*zoom))/2; 
  m_scaleX = (float)(m_viewRect->width()*zoom)/(float)m_viewRect->width();
  m_scaleY = (float)(m_viewRect->height()*inqRatio()*zoom)/(float)(m_viewRect->height());
  

  m_paint.setViewport(m_origX, m_origY,
		      (int)(m_viewRect->width()*zoom),
		      (int)(m_viewRect->height()*inqRatio()*zoom));  
  
  m_paint.setWindow(m_viewRect->left(),
                    m_viewRect->top(),
                    m_viewRect->width(),
                    -m_viewRect->height());
  
  m_paint.setClipRect(m_paint.viewport().left(),
                       m_paint.viewport().top(),
                       m_paint.viewport().width() + 1,
                       m_paint.viewport().height()+ 1);
 
  m_paint.setClipping(true);

  if(m_imagesEnabled)
    {
      m_prevCursor.reset();
      paintImages();
      
      m_paint.setWindow(m_viewRect->left()   *256,
                        m_viewRect->top()    *256,
                        m_viewRect->width()  *256,
                        -m_viewRect->height()*256);
      drawDtiLines();
    } 
 
  m_paint.setWindow(m_viewRect->left(),
                    m_viewRect->top(),
                    m_viewRect->width(),
                    -m_viewRect->height());
  
  if(m_zooming) drawZoomRectangle();   
  if(m_mode == Masking) paintGraphics();
  m_paint.setWindow(m_viewRect->left()   *2,
                    m_viewRect->top()    *2,
                    m_viewRect->width()  *2,
                    -m_viewRect->height()*2);
  
//    if(hasMouseTracking())
//    m_paint.setBrush(QColor(255, 255, 0));
//    else m_paint.setBrush(QColor(255, 0, 0));
//    m_paint.drawRoundRect(m_viewRect->left(),m_viewRect->top(),10, 10);
  
  m_paint.setPen(QColor(128, 128, 128));

  if(m_crossHairsOn)
    {
      if(m_mode == Masking)
        {
          if (!m_shape) { drawCrossHairs(); }
          else if (m_shape->empty()) { drawCrossHairs(); }
        }
      else
        {
          drawCrossHairs();
        }
    }

  m_forceRender = false;
  m_imagesEnabled = true; 

  if(m_opts.inqShowSliceLabels())
    {
      QFont font = m_paint.font();
      font.setStyleHint(QFont::SansSerif, 
			QFont::NoAntialias);
      font.setPointSize(12);
      m_paint.setPen(QColor(255, 255, 255));
      //m_paint.setFont(font);
      m_paint.setViewport(geometry());
      m_paint.setClipping(false);
      m_paint.setWindow(geometry());
      m_paint.drawText(1, 10, inqLocationText());
    }

  m_paint.end();
  
  bitBlt(this, 0, 0, &pm);
}

QString AxialWidget::inqLocationText() const
{
  QString str;
  if(m_opts.inqUnitsAreVoxels())
    str = QString("z=%1").arg(m_slice);
  else {
    float x(0), y(0), z(0);
    ImageInfo::Handle im_info(m_overlayList->getActiveMetaImage()->getInfo());

    im_info->voxToMMCoord(0, 0, m_slice, x, y, z);
    str = QString("z=%1 mm").arg(z);
  }
  return str;
}

QString SagittalWidget::inqLocationText() const
{
  QString str;
  if(m_opts.inqUnitsAreVoxels())
    str = QString("x=%1").arg(m_slice);
  else {
    float x(0), y(0), z(0);
    int radiogX(m_slice);
    ImageInfo::Handle im_info(m_overlayList->getActiveMetaImage()->getInfo());
    if(!im_info->isStoredRadiological())
      radiogX = im_info->inqX()-1-radiogX;
    im_info->voxToMMCoord(radiogX, 0, 0, x, y, z);
    str = QString("x=%1 mm").arg(x);
  }
  return str;
}

QString CoronalWidget::inqLocationText() const
{
  QString str;
  if(m_opts.inqUnitsAreVoxels())
    str = QString("y=%1").arg(m_slice);
  else {
    float x(0), y(0), z(0);
    ImageInfo::Handle im_info(m_overlayList->getActiveMetaImage()->getInfo());
    im_info->voxToMMCoord(0, m_slice, 0, x, y, z);
    str = QString("y=%1 mm").arg(y);
  }
  return str;
}

void SliceWidget::paintGraphics()
{
  TRACKER("SliceWidget::paintGraphics()");
  if(m_shape)
    m_shape->draw();
}


void SliceWidget::paintImages()
{ 
  TRACKER("SliceWidget::paintImages()");

  bool isBottomImage(true);

  if(!displayBufferValid())
    {
      m_store->resetPos();
      
      while(!m_store->currentEmpty())
	{
	  ImageData::Handle i = m_store->current();
	  if(i->inqVisibility() && (i->inqDtiDisplay() != DtiDisplay(Lines)) && 
	     (i->inqDtiDisplay() != DtiDisplay(LinesRGB)) )
            {
	      ImageBuffer::blendBuffers(m_displayPixels, i->getBuffer(), i->inqTransparency(), 
					isBottomImage, inqWidth()*inqHeight());    
              isBottomImage = false;
            }

	  m_store->next();
	}
    }
  if(isBottomImage)setToZero(m_displayPixels);
    
  if(!displayPixmapValid())
    {  
      if(QImage::BigEndian == QImage::systemByteOrder())reorderBytes(m_displayPixels);
      
      QImage image((unsigned char*)m_displayPixels.get(),                   
		   inqWidth(),
		   inqHeight(),
		   32,NULL,0,QImage::IgnoreEndian);         

      m_displayPixmap->convertFromImage(image);
    }
 
  m_paint.drawPixmap(m_dataRect->left(), m_dataRect->bottom(),
		     *m_displayPixmap,
                     m_dataRect->left(), m_dataRect->bottom(),
                     m_dataRect->width(),m_dataRect->height());
   
  validateDisplayPixmap();
  validateDisplayBuffer();
}

void SliceWidget::crossHairMode(bool mode)
{
  m_crossHairsOn = mode;

  QWidget::repaint();
}

/** 
 * A fixed slice will not change with the cursor. Use this if
 * you want to force the slice to stay the same regardless of the
 * cursor "depth" value.
 * 
 * @param fixed true if the slice is fixed.
 */
void SliceWidget::setSliceIsFixed(bool fixed)
{
  m_sliceIsFixed = fixed;
}

void SliceWidget::setSlice(int s, int v)
{  
  TRACKER("SliceWidget::setSlice");
  CHECKPOINT();

  ImageDisplaySetting::Handle ds = m_overlayList->getMainMetaImage()->getDs();

  m_imagesEnabled = false;
  if(( m_slice  != s) ||
     ((m_volume != v) && (ds->inqDtiDisplay() == DtiDisplay(None))) || 
     ( m_noSliceSet) || 
     ( m_forceRender))
  {
      if(!m_sliceIsFixed || (m_volume != v)|| (m_forceRender)) 
      {
      if(!m_sliceIsFixed)m_slice = s;
      m_volume = v;
      m_noSliceSet = false;
      invalidateImageBuffers();
      m_imagesEnabled = true;
      }
  }

  QWidget::repaint();
}
    



void CoronalWidget::setImageCursor(int x, int y, int z, int v)
{  
  TRACKER("CoronalWidget::setImageCursor");
  CHECKPOINT();
  setSlice(y,v);
}

void AxialWidget::setImageCursor(int x, int y, int z, int v)

{  
  TRACKER("AxialWidget::setImageCursor");
  CHECKPOINT();
  setSlice(z,v);
}

//! @brief Cursor objects update the SliceWidget via this method
void SliceWidget::update(const Cursor::Handle& c)
{
  TRACKER("SliceWidget::update(const Cursor* c)");
  CHECKPOINT();

  if(c->inqRepaint()){m_forceRender = true;invalidateImageBuffers();}
 
  setImageCursor(c->inqX(), c->inqY(), c->inqZ(), c->inqV());
}

//! @brief BriCon objects update the SliceWidget via this method
void SliceWidget::update(const BriCon* b)
{
  TRACKER("SliceWidget::update(const BriCon* b)");

  invalidateImageBuffers();
  QWidget::repaint();
}

//! @brief DrawSettings objects update the SliceWidget via this method
void SliceWidget::update(const DrawSettings* d)
{
  TRACKER("SliceWidget::update(const DrawSettings* d)");
  if(m_mode == Masking) 
    switch(d->inqMode())
      {
      case DrawSettings::FreeHand:
	setCursor(*m_penCursor);
	break;
      case DrawSettings::Erase:
	setCursor(*m_eraserCursor);
	break;    
      case DrawSettings::Fill:
	setCursor(*m_fillCursor);
	break;
      }
}

//! @brief OverlayList objects update the SliceWidget via this method
void SliceWidget::update(const OverlayList* i, OverlayListMsg msg)
{  
  TRACKER("SliceWidget::update(const OverlayList* i, OverlayListMsg msg)");

  if( msg != Select )
    {
      invalidateDisplayBuffer();
      if( (msg == Add) || (msg == Rem) ) {
	std::for_each(m_briconList.begin(), m_briconList.end(), 
		      Detach(this));
	std::for_each(m_overlayList->begin(), m_overlayList->end(), 
		      Attach(this, m_briconList));
      }
      loadStore();
      invalidateImageBuffers();
      QWidget::repaint();
    }

//   switch(msg)
//   {
//    case OverlayListMsg(Select):      MESSAGE("Select");break;     
//    case OverlayListMsg(Visibility):  MESSAGE("Visibility");
//                                      invalidateDisplayBuffer();
//                                      QWidget::repaint();break;
//    case OverlayListMsg(Transparency):MESSAGE("Transparency");
//                                      invalidateDisplayBuffer();
//                                      QWidget::repaint(); break;
//    case OverlayListMsg(Order):       MESSAGE("Order");
//                                      loadStore();invalidateImageBuffers(); 
//                                      QWidget::repaint(); break;
//    case OverlayListMsg(Add):         MESSAGE("Add");
//    case OverlayListMsg(Rem):         MESSAGE("Rem");
//      std::for_each(m_briconList.begin(), m_briconList.end(), Detach(this));
//      std::for_each(m_overlayList->begin(), m_overlayList->end(), Attach(this, m_briconList));
//      loadStore();
//    case OverlayListMsg(LookUpTable): MESSAGE("LookUpTable");invalidateImageBuffers(); 
//                                      QWidget::repaint();
//                                      break;   
//    case OverlayListMsg(ModImage):    MESSAGE("ModImage");invalidateImageBuffers(); 
//                                      QWidget::repaint();
//                                      break;
//    case OverlayListMsg(DtiMode):     MESSAGE("DtiMode");
//                                      invalidateImageBuffers(); 
//                                      QWidget::repaint();
//                                      break;
//   }
}

void SliceWidget::mouseMoveEvent(QMouseEvent *e)
{ 
  QPoint w = convMouseToWorld(QPoint(e->x(),e->y()));

  bool lButton(e->state() & LeftButton);
  bool mButton(e->state() & MidButton);
  bool rButton(e->state() & RightButton);

  int dx = m_startX - w.x();
  int dy = m_startY - w.y();
  
  if(lButton) {
    switch(m_mode)
      {
      case None:
      case Cursing: 
        cursorEvent(w.x(), w.y()); 
  
        break;
      case Pan:
        {
          transEvent(dx, dy);
          QPoint v = convMouseToWorld(QPoint(e->x(),e->y()));
          setStartMove(v.x(),v.y()); 
        }
        break;
      case Zoom:
        zoomEvent(w.x(),w.y());
        break;
      case Masking:
        if(m_drawSettings->inqMode() != DrawSettings::Fill)
	  drawEvent(w.x(),w.y());
        break;
      default:
        break;
      }
  } else if(mButton) {
    transEvent(dx, dy);
    QPoint v = convMouseToWorld(QPoint(e->x(),e->y()));
    setStartMove(v.x(),v.y()); 
  } else if(rButton) {
    zoomEvent(w.x(),w.y());
  }
}

void SliceWidget::mousePressEvent(QMouseEvent *e)
{
  QPoint w = convMouseToWorld(QPoint(e->x(),e->y()));

  bool lButton(e->button() == LeftButton);
  bool mButton(e->button() == MidButton);
  bool rButton(e->button() == RightButton);

  setStartMove(w.x(),w.y());

  if(lButton) {
    switch(m_mode)
      {
      case None:
      case Cursing: 
        cursorEvent(w.x(), w.y()); break;
      case Pan: 	  
        setCursor(*m_panCursor);
        break;
      case Zoom: 
        if(e->state() & ControlButton) 
          zoomOut(w.x(), w.y());
        else {
          setCursor(*m_zoomCursor);
          m_zooming = true;
        }
        break;

      case Masking:       

        if(layerValidForDrawing())
          {  
            MetaImage::Handle mi = m_overlayList->getActiveMetaImage();
            if(mi)
              {
                m_shape = Shape::create(&m_paint,
                                        mi->getImage()->getVolume(m_volume),
                                        m_orient,
                                        m_slice);
                
		switch(m_drawSettings->inqMode())
		  {
		  case DrawSettings::FreeHand:
		  case DrawSettings::Erase:
                    drawEvent(w.x(),w.y());
		    break;
		  case DrawSettings::Fill:
                    floodEvent(w.x(),w.y());
		    break;
		  }
                mi->getInfo()->setTarnished(true);
              }
          }
        
        break;

      default:
        break;
      }
  } else if(mButton) {
  } else if(rButton) {
    if(e->state() & ControlButton) 
      zoomOut(w.x(), w.y());
    else
      m_zooming = true;
  }
}

void SliceWidget::mouseReleaseEvent(QMouseEvent *e)
{

  bool lButton(e->button() == LeftButton);
  bool mButton(e->button() == MidButton);
  bool rButton(e->button() == RightButton);

  if(lButton) {
    switch(m_mode)
      {
      case None:
      case Cursing: 
        break;

      case Pan:  
        setCursor(*m_panCursor);
        break;

      case Zoom:    
        if(!(e->state() & ControlButton))
          setViewRect(m_zoomRect->left(),m_zoomRect->bottom(),
                      m_zoomRect->right(),m_zoomRect->top());
        
        m_zoomRect->setRect(0,0,0,0);
        setCursor(*m_zoomCursor);
        m_zooming = false;
        break;

      case Masking:
        if(m_drawSettings->inqMode() != DrawSettings::Fill)
	  commitGraphics();
        break;

      default:
        break;
      }
  } else if(mButton) {
  } else if(rButton && !(e->state() & ControlButton)) {
    setViewRect(m_zoomRect->left(),m_zoomRect->bottom(),
		m_zoomRect->right(),m_zoomRect->top());
    m_zoomRect->setRect(0,0,0,0);
    m_zooming = false;    
  }
}

void SliceWidget::commitGraphics()
{
  TRACKER("SliceWidget::commitGraphics()");

  MetaImage::Handle mi = m_overlayList->getActiveMetaImage();
  if(mi.get() && (m_shape.get() != NULL))
    {
      Shape::Handle undoBuffer = m_shape->getBuffer();
      m_undoList.push_back(undoBuffer);
      if(m_undoList.size() > 20) m_undoList.pop_front();
       
      m_shape->commit();

      mi->getInfo()->setTarnished(true);
    }
 
     m_cursor->repaint();
}


void SliceWidget::enterEvent( QEvent *e )
{
  TRACKER("SliceWidget::enterEvent(QEvent)");
//   topLevelWidget()->setCaption("Here I am!");
//   if(topLevelWidget()->isFocusEnabled()) {
//     MESSAGE("topLevelWidget()->isFocusEnabled()");
//     setFocus();
//   }
//  grabKeyboard();
}

void SliceWidget::leaveEvent( QEvent *e )
{
  TRACKER("SliceWidget::leaveEvent(QEvent)");
//   clearFocus();
//   QWidget::repaint();
//  releaseKeyboard();
}

const QPoint  SliceWidget::convMouseToWorld(const QPoint & p) const
{   
  QPoint world;

  world.setX ((int)((float)(p.x() - m_origX)/(float)m_scaleX) + m_viewRect->left());
  world.setY ((int)((float)((height() - p.y()) - m_origY)/(float)m_scaleY) +m_viewRect->bottom());
  return world;
}

void SliceWidget::transEvent(int dx, int dy)
{
  TRACKER("SliceWidget::transEvent");
  m_viewRect->translate(dx,dy);
  setDataRect();
  QWidget::repaint();
}

void SliceWidget::zoomEvent(int x, int y)
{
  m_zoomRect->setRect(m_startX,m_startY,x,y);
  m_zoomRect->setUnion(m_viewRect); 

  QWidget::repaint();
}

void SliceWidget::briconEvent(int dx, int dy)
{
}

void AxialWidget::setZoom(int factor)
{  
  TRACKER("AxialWidget::setZoom(int)");

  float f = (factor/100.0);
  int dx = (int)(inqWidth()/f/2);
  int dy = (int)(inqHeight()/f/2);

  setViewRect(m_cursor->inqX()-dx,m_cursor->inqY()-dy,
	      m_cursor->inqX()+dx,m_cursor->inqY()+dy);
  QWidget::repaint();
}

void CoronalWidget::setZoom(int factor)
{  
  TRACKER("CoronalWidget::setZoom(int)");  

  float f = (factor/100.0);
  int dx = (int)(inqWidth()/f/2);
  int dy = (int)(inqHeight()/f/2);

  setViewRect(m_cursor->inqX()-dx,m_cursor->inqZ()-dy,
              m_cursor->inqX()+dx,m_cursor->inqZ()+dy);
  QWidget::repaint();
}

void SagittalWidget::setZoom(int factor)
{  
  TRACKER("SagittalWidget::setZoom(int)");  

  float f = (factor/100.0);
  int dx = (int)(inqWidth()/f/2);
  int dy = (int)(inqHeight()/f/2);

  setViewRect(m_cursor->inqY()-dx,m_cursor->inqZ()-dy,
              m_cursor->inqY()+dx,m_cursor->inqZ()+dy);
  QWidget::repaint();
}

void SliceWidget::zoomOut(int x,int y)
{
  if(!m_zoomHistory.empty())
    {
      m_viewRect = m_zoomHistory.top().first;
      m_trueScale = m_zoomHistory.top().second;
      m_zoomHistory.pop();
    }

  setDataRect();
  QWidget::repaint(); 

}

void SliceWidget::resetZoom()
{
  TRACKER("SliceWidget::resetZoom()");

  while(!m_zoomHistory.empty()){m_zoomHistory.pop();}

  m_viewRect = Rect::createRect(0,0,inqWidth(),inqHeight());
  m_zoomRect = Rect::createRect(0,0,inqWidth(),inqHeight());
  m_dataRect = Rect::createRect(0,0,inqWidth(),inqHeight());
  
  m_trueScale = true;
  
  emitZoomFactor(100);

  initZoom();

  QWidget::repaint();
}

void SliceWidget::setViewRect(int startX,int startY,int curX,int curY)
{ 
  TRACKER("SliceWidget::setViewRect(int startX,int startY,int curX,int curY)");


  if(!(startX == curX && startY == curY))
    { 
      m_zoomHistory.push(std::make_pair(m_viewRect->clone(),m_trueScale));

      float winRatio = height()/(float)width();

      m_viewRect->setRect(startX,startY,curX,curY);
      m_trueScale = false;
      
      if( m_viewRect->height() >= m_viewRect->width()*winRatio)
        m_viewRect->setWidth((int)((m_viewRect->height()*inqRatio())/winRatio));
      else if( m_viewRect->height() < m_viewRect->width()*winRatio)
        m_viewRect->setHeight((int)((m_viewRect->width()/inqRatio())*winRatio));  

      //fix zero sizes to avoid divide by zeros later on
      if(m_viewRect->width() == 0){m_viewRect->setWidth(m_viewRect->width() +1);}
      if(m_viewRect->height() == 0){m_viewRect->setHeight(m_viewRect->height() +1);}

      setDataRect();  
    }
  else if((startX == curX) && (startY == curY))
    {  
      //if right button just pressed not moved
      if(m_trueScale)
        {
          m_zoomHistory.push(std::make_pair(m_viewRect->clone(),m_trueScale));  
          m_trueScale = false;
        }  
    }  

  QWidget::repaint(); 
}

void SliceWidget::setDataRect()
{
  TRACKER("SliceWidget::setDataRect");
 
  m_dataRect->setRect(0,0,inqWidth(),inqHeight());
  m_dataRect->setUnion(m_viewRect);
}

void SliceWidget::loadStore()
{  
  TRACKER("SliceWidget::loadStore");
  m_store = ImageDataStore::create(m_overlayList);
}

void SliceWidget::renderBuffer()
{   
  TRACKER("SliceWidget::renderBuffer");
  CHECKPOINT();

  m_store->resetPos();
      
  bool isBottomImage(true);

  while(!m_store->currentEmpty())
    {
      ImageData::Handle  i(m_store->current());
      MetaImage::Handle mi(i->getMetaImage());

      if(i->inqDtiDisplay() == DtiDisplay(None))
	{
	  ColorRGBAHandle buffer(bufferVolume(mi));

	  LookUpTable::Handle slut(mi->getDs()->inqSecondaryLookUpTable());
	  if(mi->getDs()->inqUseSecondaryLookUpTable()) {
	    MetaImage::Handle si(mi->clone());
	    BriCon::Handle    sb(si->getDs()->inqBriCon());
	    sb->setRange(-(sb->inqMin()), -(sb->inqMax()));
	    si->getDs()->setLookUpTable(slut);
	    ColorRGBAHandle sbuffer(bufferVolume(si));
 	    ImageBuffer::blendBuffers(buffer, sbuffer, 1, false, inqWidth() * inqHeight());
	  }
	  i->setBuffer(buffer);
	}
      else
        i->setBuffer(dtiVolume(i->getMetaImage()));

      isBottomImage = false;
      m_store->next();
    }

  validateImageBuffers();
}

void SliceWidget::setMode(SliceWidget::Mode m)
{  
  TRACKER("SliceWidget::setMode");
  m_mode = Mode(m);
  switch(m_mode)
    {
    case None:
    case Cursing:
      setCursor(*m_crossCursor);
      break;
    case Pan:
      setCursor(*m_panCursor);
      break;
    case Zoom:
      setCursor(*m_zoomCursor);
      break;
    case Masking:
      switch (m_drawSettings->inqMode())
        {
        case DrawSettings::FreeHand:
          setCursor(*m_penCursor);
          break;
        case DrawSettings::Erase:
          setCursor(*m_eraserCursor);
          break;
        case DrawSettings::Fill:
          setCursor(*m_fillCursor);
          break;
        }
      break;
    default:
      break;
    }
  QWidget::repaint();
}

void SliceWidget::drawZoomRectangle()
{  
  TRACKER("SliceWidget::drawZoomRectangle");
  m_paint.setPen(QColor(0,255,0));
  m_paint.drawRect(m_zoomRect->left(),
                   m_zoomRect->bottom(),
                   m_zoomRect->width()+1,
                   m_zoomRect->height()+1);
}

void SliceWidget::drawCrossHairs()
{
  if(m_prevCursor.get() != NULL)
      drawCrossHairLines(m_prevCursor,m_prevSlice);
 
  drawCrossHairLines(m_cursor,m_slice);
 
  m_prevCursor = m_cursor->clone();
  m_prevSlice = m_slice;
}

void SliceWidget::reorderBytes(ColorRGBAHandle buffer)
{
  ImageBuffer::reorderBytes(buffer,inqWidth() * inqHeight());
}

void SliceWidget::setToZero(ColorRGBAHandle buffer)
{
  ImageBuffer::setToZero(buffer,inqWidth() * inqHeight());
}

void SliceWidget::showSlice()
{
  m_timer->stop();
  resize(width(),height()+1);
}

void SliceWidget::keyPressEvent(QKeyEvent* e)
{
  if(e->key() == Key_PageUp)  pageUpPressed();
  else if(e->key() == Key_PageDown)pageDownPressed();
  else if(e->key() == Key_Up)      moveCursor( 0, 1);
  else if(e->key() == Key_Down)    moveCursor( 0,-1);
  else if(e->key() == Key_Left)    moveCursor(-1, 0);
  else if(e->key() == Key_Right)   moveCursor( 1, 0);
  else if(e->key() == Key_Control) setEraseMode(true);
  else if(e->key() == Key_Shift)   setFillMode(true);
}

void SliceWidget::keyReleaseEvent(QKeyEvent* e)
{
  if(e->key() == Key_Control) setEraseMode(false);
  if(e->key() == Key_Shift)   setFillMode(false);
}

void SliceWidget::setEraseMode(bool state)
{
  if(m_mode == Masking)
    {
      if(state){m_drawSettings->setMode(DrawSettings::Erase);}
      else{     m_drawSettings->setPrevMode();}
    }
}

void SliceWidget::setFillMode(bool state)
{
  if(m_mode == Masking)
    {
      if(state){m_drawSettings->setMode(DrawSettings::Fill);}
      else{     m_drawSettings->setPrevMode();}
    }
}

void SliceWidget::pageUpPressed()
{
  int s = m_slice + 1;
  if(s >= inqDepth()) s = 0;
  setCursorSlice(s);
  m_cursor->repaint();
}

void SliceWidget::pageDownPressed()
{
  int s = m_slice - 1;
  if(s < 0) s = (inqDepth() -1);
  setCursorSlice(s);
  m_cursor->repaint();
}

void SliceWidget::drawEvent(int x, int y)
{  
  int   size  = m_drawSettings->inqPenSize();
  int   value = m_drawSettings->inqPenValue();

  if(m_drawSettings->linkCursorOn())
    cursorEvent(x, y);

  if (m_shape)
    {
      if(layerValidForDrawing()){m_shape->addVertex(x, y, size,(float)value);}
      if(m_shape->size() > 1)m_imagesEnabled = false;
    }
  QWidget::repaint();
  m_imagesEnabled = true;
}

void SliceWidget::floodEvent(int x, int y)
{
  int   value = m_drawSettings->inqPenValue();

  if (m_shape)
    {
      if(layerValidForDrawing()){m_shape->floodFill(x, y,(float)value);}

      Shape::Handle undoBuffer = m_shape->getFloodBuffer();
      m_undoList.push_back(undoBuffer);
      if(m_undoList.size() > 20)m_undoList.pop_front();
    }
  m_cursor->repaint();
}

bool SliceWidget::layerValidForDrawing()
{  
  bool result(false);

  MetaImage::Handle mi = m_overlayList->getActiveMetaImage();
  
  if(mi)
    { 
      if(!mi->inqReadOnly())
        {result = true;}
      else
        {emit message("Warning: Drawing disabled. The currently selected layer is locked.", 4000);}
    }
  else
    {
      emit message("Warning: No valid overlay selected",2000);
    }

  return result;
}

/***************************************
*Sagittal Widget
*
****************************************/
void SagittalWidget::setImageCursor(int x, int y, int z, int v)
{  
  TRACKER("SagittalWidget::setImageCursor");
  CHECKPOINT();
  setSlice(x,v);
}

//! @brief Renders the cross hairs.
//! @param c the Cursor::Handle for the location at which the cross hairs will
//!        be drawn
//! @param slice the highlighted slice (rendered brighter than the 
//! crosshairs of the other slices)
void SagittalWidget::drawCrossHairLines(const Cursor::Handle c,int slice)
{ 
  if(c->inqX() == slice) {
    m_paint.setRasterOp(Qt::XorROP);
    m_paint.setPen( QColor(0,255,0));
    if (m_opts.inqShowCursorGap())
      drawBrokenCrossHairLines( (c->inqY()*2)+1, (c->inqZ()*2)+1 );
    else
      drawSimpleCrossHairLines( (c->inqY()*2)+1, (c->inqZ()*2)+1 );
    m_paint.setRasterOp(Qt::CopyROP);
  }
}

void SagittalWidget::cursorEvent(int x, int y)
{
  m_cursor->setCursor(m_slice, x, y);
}

int SagittalWidget::inqWidth() const 
{
  return m_overlayList->inqY(); 
}

int SagittalWidget::inqHeight() const 
{ 
  return m_overlayList->inqZ(); 
}

int SagittalWidget::inqDepth() const
{
  return m_overlayList->inqX();
}

float SagittalWidget::depthRatio() const
{
  ImageInfo::Handle info(m_overlayList->getMainImage()->getInfo());  
  return info->inqNoDimensions()? 1.0 : fabs(info->inqXDim() / info->inqYDim());
}


float SagittalWidget::inqRatio() const
{
  ImageInfo::Handle info(m_overlayList->getMainImage()->getInfo());  
  return info->inqNoDimensions()? 1.0 : fabs(info->inqZDim() / info->inqYDim());
}

ColorRGBAHandle SagittalWidget::dtiVolume(MetaImage::Handle mi)
{
  return ImageBuffer::sagittalDtiBuffer(mi,m_slice);
}

ColorRGBAHandle SagittalWidget::bufferVolume(MetaImage::Handle mi)
{
  TRACKER("SagittalWidget::bufferVolume");
  CHECKPOINT();
  return ImageBuffer::sagittalBuffer(mi,m_slice,mi->getDs()->inqCurrentVolume());
}

void SagittalWidget::setCursorSlice(short s)
{
  m_cursor->setCursor(s, m_cursor->inqY(), m_cursor->inqZ());
}  

// QColor SliceWidget::getDTIVectorColor(const Image::Handle& i,
// 				      int x, int y, int z)
// {
//   Volume::Handle vR(i->getImage()->getVolume(0));
//   Volume::Handle vG(i->getImage()->getVolume(1));
//   Volume::Handle vB(i->getImage()->getVolume(2));

//   int r = fabs(vR->value(m_slice, y, z)) * 255;
//   int g = fabs(vG->value(m_slice, y, z)) * 255;
//   int b = fabs(vB->value(m_slice, y, z)) * 255;

//   return QColor(r, g, b);  
// }

void SagittalWidget::drawDtiLines()
{
  m_store->resetPos();

  unsigned int c = 0;

  while(!m_store->currentEmpty())
    {
      ImageData::Handle i(m_store->current());
      ImageInfo::Handle info(i->getImage()->getInfo());

      if( i->inqVisibility() && ( (i->inqDtiDisplay() == DtiDisplay(Lines)) ||
				  (i->inqDtiDisplay() == DtiDisplay(LinesRGB)) ) )
	{
	  int yVec,zVec;

	  Volume::Handle vG(i->getImage()->getVolume(1));
	  Volume::Handle vB(i->getImage()->getVolume(2));

	  unsigned int width  = vG->inqY();
	  unsigned int height = vG->inqZ();
	  float minDimension = std::min(info->inqXDim(), std::min(info->inqYDim(), info->inqZDim()));

	  m_paint.setPen(m_dtiColors[c]);

	  for(unsigned int z = 0; z < height; ++z)
	    {
	      for( unsigned int y = 0; y < width; ++y)
		{
		  yVec = int(vG->value(m_slice, y, z)*( 255.0/info->inqYDim() )*minDimension);
		  zVec = int(vB->value(m_slice, y, z)*( 255.0/info->inqZDim() )*minDimension);
		  if(yVec != 0 || zVec != 0)
		    {
		      if( i->inqDtiDisplay() == DtiDisplay(LinesRGB) )
			m_paint.setPen(getDTIVectorColor(i->getImage(), m_slice, y, z));
		      m_paint.drawLine((y*256) + 128 - int( 0.5 * yVec),(z*256) + 128 - int( 0.5 * zVec),
				       (y*256) + 128 + int( 0.5 * yVec),(z*256) + 128 + int( 0.5 * zVec));
		    }
		}
	    }
	}
      if(info->isDtiCompatible())
	c = (c + 1) % m_dtiColors.size();
      m_store->next();
    }
}

void SagittalWidget::moveCursor(short dx, short dy)
{
  m_cursor->setCursor(m_slice,
                      m_cursor->inqY() + dx,
                      m_cursor->inqZ() +dy);
}



/***********************************
*Axial Widget
*
************************************/
void AxialWidget::cursorEvent(int x, int y)
{  
  m_cursor->setCursor(x, y, m_slice);
}

int AxialWidget::inqWidth() const 
{ 
  return m_overlayList->inqX(); 
}

int AxialWidget::inqHeight() const 
{
  return m_overlayList->inqY(); 
}

int AxialWidget::inqDepth() const
{
  return m_overlayList->inqZ();
}

float AxialWidget::depthRatio() const
{
  ImageInfo::Handle info(m_overlayList->getMainImage()->getInfo());
  return info->inqNoDimensions()? 1.0 : fabs(info->inqZDim() / info->inqXDim());
}

float AxialWidget::inqRatio() const
{
  ImageInfo::Handle info(m_overlayList->getMainImage()->getInfo());  
  return info->inqNoDimensions()? 1.0 : fabs(info->inqYDim() / info->inqXDim());
}

ColorRGBAHandle AxialWidget::dtiVolume(MetaImage::Handle mi)
{
  return ImageBuffer::axialDtiBuffer(mi,m_slice);
}

ColorRGBAHandle AxialWidget::bufferVolume(MetaImage::Handle mi)
{
  TRACKER("AxialWidget::bufferVolume");
  CHECKPOINT();
  return ImageBuffer::axialBuffer(mi,m_slice,mi->getDs()->inqCurrentVolume());
}

void SliceWidget::drawSimpleCrossHairLines(int x, int y)
{
  int h = inqHeight() * 2;
  int w = inqWidth() * 2;

  m_paint.drawLine(x, 0, x, h);
  m_paint.drawLine(0, y, w, y);
}

void SliceWidget::drawBrokenCrossHairLines(int x, int y)
{
  int gapsz = m_opts.inqCursorGapSize();
  int h = inqHeight() * 2;
  int w = inqWidth() * 2;

  // vertical bar
  m_paint.drawLine( x, 0,         x, y - gapsz );
  m_paint.drawLine( x, y + gapsz, x, h         );

  // horizontal bar
  m_paint.drawLine( 0,         y, x - gapsz, y );
  m_paint.drawLine( x + gapsz, y, w,         y );
}

//! @brief Renders the cross hairs.
//! @param c the Cursor::Handle for the location at which the cross hairs will
//!        be drawn
//! @param slice the highlighted slice (rendered brighter than the 
//! crosshairs of the other slices)
void AxialWidget::drawCrossHairLines(const Cursor::Handle c,int slice)
{  
  if(c->inqZ() == slice) {
    m_paint.setRasterOp(Qt::XorROP);
    m_paint.setPen( QColor(0,255,0));
    if (m_opts.inqShowCursorGap())
      drawBrokenCrossHairLines( (c->inqX()*2)+1, (c->inqY()*2)+1 );
    else
      drawSimpleCrossHairLines( (c->inqX()*2)+1, (c->inqY()*2)+1 );
    m_paint.setRasterOp(Qt::CopyROP);
  }
}

void AxialWidget::moveCursor(short dx, short dy)
{
  m_cursor->setCursor(m_cursor->inqX() + dx, 
                      m_cursor->inqY() + dy, 
                      m_slice);
}

void AxialWidget::setCursorSlice(short s)   
{      
  m_cursor->setCursor(m_cursor->inqX(), m_cursor->inqY(),s); 
}

void AxialWidget::drawDtiLines()
{
  m_store->resetPos();
  
  unsigned int c = 0;

  while(!m_store->currentEmpty())
    {
      ImageData::Handle i(m_store->current());
      ImageInfo::Handle info(i->getImage()->getInfo());

      if( i->inqVisibility() && ( (i->inqDtiDisplay() == DtiDisplay(Lines)) || 
				  (i->inqDtiDisplay() == DtiDisplay(LinesRGB)) ) )
	{
	  int xVec,yVec;
	  
	  Volume::Handle vR(i->getImage()->getVolume(0));
	  Volume::Handle vG(i->getImage()->getVolume(1));
	  
	  unsigned int width  = vR->inqX();
	  unsigned int height = vR->inqY();
	  float minDimension = std::min(info->inqXDim(), std::min(info->inqYDim(), info->inqZDim()));
	  
	  //qDebug("XDim = %f, YDim = %f, ZDim = %f, minDimension = %f", 
	  //      info->inqXDim(), info->inqYDim(), info->inqZDim(), minDimension);

	  m_paint.setPen(m_dtiColors[c]);
	  
	  for(unsigned int y = 0; y < height; ++y)
	    { 
	      for( unsigned int x = 0; x < width; ++x) 
		{          
		  xVec = int(vR->value(x, y, m_slice)*( 255.0/info->inqXDim() )*minDimension);
		  yVec = int(vG->value(x, y, m_slice)*( 255.0/info->inqYDim() )*minDimension);
		  //qDebug("xVec = %d, yVec = %d", xVec, yVec);
		  if(xVec != 0 || yVec != 0)
		    {
		      if( i->inqDtiDisplay() == DtiDisplay(LinesRGB) )
			m_paint.setPen(getDTIVectorColor(i->getImage(), x, y, m_slice));
		      m_paint.drawLine((x*256) + 128 - int( 0.5 * xVec),(y*256) + 128 - int( 0.5 * yVec),
				       (x*256) + 128 + int( 0.5 * xVec),(y*256) + 128 + int( 0.5 * yVec));
		    }
		}
	    }
	}
      if(info->isDtiCompatible())
	c = (c + 1) % m_dtiColors.size();
      m_store->next();
    }
}




/***************************************
*Coronal Widget
*
****************************************/

//! @brief Renders the cross hairs.
//! @param c the Cursor::Handle for the location at which the cross hairs will
//!        be drawn
//! @param slice the highlighted slice (rendered brighter than the 
//! crosshairs of the other slices)
void CoronalWidget::drawCrossHairLines(const Cursor::Handle c,int slice)
{ 
  if(c->inqY() == slice) {
    m_paint.setRasterOp(Qt::XorROP);
    m_paint.setPen( QColor(0,255,0));
    if (m_opts.inqShowCursorGap())
      drawBrokenCrossHairLines( (c->inqX()*2)+1, (c->inqZ()*2)+1 );
    else
      drawSimpleCrossHairLines( (c->inqX()*2)+1, (c->inqZ()*2)+1 );
    m_paint.setRasterOp(Qt::CopyROP);
  }
}

void CoronalWidget::cursorEvent(int x, int y)
{
  m_cursor->setCursor(x, m_slice, y); 
}

int CoronalWidget::inqWidth() const 
{ 
  return m_overlayList->inqX();
}

int CoronalWidget::inqHeight() const 
{
  return m_overlayList->inqZ();
}

int CoronalWidget::inqDepth() const
{
  return m_overlayList->inqY();
}

float CoronalWidget::depthRatio() const
{
  ImageInfo::Handle info(m_overlayList->getMainImage()->getInfo());  
  return info->inqNoDimensions()? 1.0 : fabs(info->inqYDim() / info->inqXDim());
}

float CoronalWidget::inqRatio() const
{
  ImageInfo::Handle info(m_overlayList->getMainImage()->getInfo());  
  return info->inqNoDimensions()? 1.0 : fabs(info->inqZDim() / info->inqXDim());
}

ColorRGBAHandle CoronalWidget::dtiVolume(MetaImage::Handle mi)
{
  return ImageBuffer::coronalDtiBuffer(mi,m_slice);
}

ColorRGBAHandle CoronalWidget::bufferVolume(MetaImage::Handle mi)
{
  TRACKER("CoronalWidget::bufferVolume");
  CHECKPOINT();
  return ImageBuffer::coronalBuffer(mi,m_slice,mi->getDs()->inqCurrentVolume());
}

void CoronalWidget::moveCursor(short dx, short dy)
{
  m_cursor->setCursor(m_cursor->inqX() + dx ,
                      m_slice,
                      m_cursor->inqZ() + dy);
}

void CoronalWidget::setCursorSlice(short s) 
{
  m_cursor->setCursor(m_cursor->inqX(), s, m_cursor->inqZ());  
}

void CoronalWidget::drawDtiLines()
{
  m_store->resetPos();

  unsigned int c = 0;
  while(!m_store->currentEmpty())
    {
    ImageData::Handle i(m_store->current());
    ImageInfo::Handle info(i->getImage()->getInfo());

    if( i->inqVisibility() && ( (i->inqDtiDisplay() == DtiDisplay(Lines)) || 
				(i->inqDtiDisplay() == DtiDisplay(LinesRGB)) ) )
      {
        int xVec,zVec;

        Volume::Handle vR(i->getImage()->getVolume(0));
        Volume::Handle vB(i->getImage()->getVolume(2));
 	float minDimension = std::min(info->inqXDim(), std::min(info->inqYDim(), info->inqZDim()));

        unsigned int width  = vR->inqX();
        unsigned int height = vR->inqZ();
       
	m_paint.setPen(m_dtiColors[c]);
        
        for(unsigned int z = 0; z < height; ++z)
        { 
          for( unsigned int x = 0; x < width; ++x) 
            {          
              xVec = int(vR->value(x, m_slice, z)*( 255.0/info->inqXDim() )*minDimension);
              zVec = int(vB->value(x, m_slice, z)*( 255.0/info->inqZDim() )*minDimension);
	      //qDebug("xVec = %d, zVec = %d", xVec, zVec);
              if(xVec != 0 || zVec != 0)
                {
		  if( i->inqDtiDisplay() == DtiDisplay(LinesRGB) )
		    m_paint.setPen(getDTIVectorColor(i->getImage(), x, m_slice, z));
                  m_paint.drawLine((x*256) + 128 - int( 0.5 * xVec),(z*256) + 128 - int( 0.5 * zVec),
                                   (x*256) + 128 + int( 0.5 * xVec),(z*256) + 128 + int( 0.5 * zVec));
                }
            }
        }
      }
    if(info->isDtiCompatible())
      c = (c + 1) % m_dtiColors.size();

    m_store->next();
    }
}
