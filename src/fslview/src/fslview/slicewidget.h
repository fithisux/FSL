/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */
/*  CCOPYRIGHT */

#if !defined (SLICEWIDGET_H)
#define SLICEWIDGET_H

#include <qcursor.h>
#include <qcolor.h>

#include <boost/shared_ptr.hpp>
#include <vector>
#include <list>
#include <stack>

#include "storage/image.h"
#include "cursor.h"
#include "overlaylist.h"
#include "briconwidget.h"
#include "drawsettings.h"
#include "lookuptable.h"
#include "metaimage.h"
#include "imagedata.h"
#include "imagedatastore.h"
#include "viewoptions.h"
#include "rect.h"
#include "shape.h"
#include <qpainter.h>

class QTimer;
class QImage;

//! 
//! @author Dave Flitney <flitney@fmrib.ox.ac.uk>
//! @author James Saunders
//!
//! @brief  Manage the display of slices through a set of volumes from a given OverlayList
//!
//! A SliceWidget observes given Cursor, BriCon and DrawSettings objects to display
//! the group of Images referred to in its OverlayList. The Images are combined and
//! displayed according to the current settings of the observed objects.
//!
class SliceWidget : public QWidget,
		    public OverlayListObserver,
		    public CursorObserver,
		    public DrawSettingsObserver,
		    public BriConObserver
{
  Q_OBJECT
public:
  typedef boost::shared_ptr<SliceWidget> Handle;
  typedef std::list<ImageData::Handle> ImageDataList;

  typedef std::pair<Rect::Handle, bool> ZoomParams;

  typedef enum {Sagittal, Axial, Coronal} Orientation;
  typedef enum {None, Cursing, Pan, Zoom, TimeSeries, Masking} Mode;

  SliceWidget(QWidget* parent, const char* name, 
              Orientation orient, Cursor::Handle c, OverlayList::Handle l,
	      DrawSettings::Handle d, std::list<Shape::Handle>& u,
	      const ViewOptions& vo);
  virtual ~SliceWidget();
  virtual void polish();

  //! @brief Returns the current display mode
  //! @return The current SliceWidget::Mode 
  int inqMode()const {return m_mode;}
  

  //! @brief Returns the current display orientation
  //! @return the current SliceWidget::Orientation 
  int inqOrient() const {return m_orient;}

  void renderBuffer();
  void setToZero(ColorRGBAHandle buffer);
  void reorderBytes(ColorRGBAHandle buffer);
  virtual void update(const Cursor::Handle& c);  
  virtual void update(const BriCon* b);
  virtual void update(const OverlayList* l, OverlayListMsg msg);
  virtual void update(const DrawSettings* d);
  virtual ColorRGBAHandle bufferVolume(MetaImage::Handle) = 0;
  virtual ColorRGBAHandle dtiVolume(MetaImage::Handle) = 0;

  QPixmap getPixmap() const;

  QColor getDTIVectorColor(const Image::Handle& i,
			   int x, int y, int z)
  {
    Volume::Handle vR(i->getVolume(0));
    Volume::Handle vG(i->getVolume(1));
    Volume::Handle vB(i->getVolume(2));
    
    int r = int(fabs(vR->value(x, y, z)) * 255);
    int g = int(fabs(vG->value(x, y, z)) * 255);
    int b = int(fabs(vB->value(x, y, z)) * 255);

    return QColor(r, g, b);  
  }

  void setNorthText(const std::string& s);
  void setSouthText(const std::string& s);
  void setWestText(const std::string& s);
  void setEastText(const std::string& s);
signals:

  void imageCursorChanged(int, int, int);
  void volumeChanged(int);
  void emitZoomFactor(int);
  void message(const QString&, int );


public slots:
  virtual void setImageCursor(int, int, int, int) = 0;

  virtual void setZoom(int) = 0;
  void zoomOut(int,int);  
  void setViewRect(int,int,int,int);
  void resetZoom();
  void setMode(SliceWidget::Mode);
  void setSlice(int,int);
  void setSliceIsFixed(bool);
  virtual int inqWidth() const = 0;
  virtual int inqHeight() const = 0;
  void enableUpdates(bool);
  void crossHairMode(bool);

protected:
  virtual void resizeEvent( QResizeEvent* );  
  virtual void paintEvent( QPaintEvent* );
  virtual void mousePressEvent(QMouseEvent*);
  virtual void mouseMoveEvent(QMouseEvent*);
  virtual void mouseReleaseEvent(QMouseEvent*);
  virtual void keyPressEvent(QKeyEvent*);
  virtual void keyReleaseEvent(QKeyEvent*);

  virtual void enterEvent(QEvent*);
  virtual void leaveEvent(QEvent*);  
  inline unsigned char briConAdjust(float) const;
 

  bool               m_sliceIsFixed;
  bool               m_updatesEnabled;
  bool               m_crossHairsOn;
  bool               m_imagesEnabled;
  int                m_slice;  
  bool               m_noSliceSet;
  int                m_volume;  
  int                m_startX, m_startY;
  Cursor::Handle     m_cursor;
  Cursor::Handle     m_prevCursor;
  int                m_prevSlice; 
  float              m_zoom; 
  int                m_origX,m_origY;
  float              m_scaleX,m_scaleY;
  Rect::Handle       m_viewRect, m_dataRect, m_zoomRect, m_winRect;

  ColorRGBAHandle        m_displayPixels;
  OverlayList::Handle    m_overlayList;
  DrawSettings::Handle   m_drawSettings;  
  QPainter               m_paint;
  ImageDataStore::Handle m_store;

  std::list<BriCon::Handle> m_briconList;
  std::stack<ZoomParams>    m_zoomHistory;
  std::vector<QColor>       m_dtiColors; 

  const ViewOptions& m_opts;

protected:
  void drawSimpleCrossHairLines(int, int);
  void drawBrokenCrossHairLines(int, int);
  
private slots:

  void showSlice();
  void setEraseMode(bool state);  
  void setFillMode(bool state);
  void pageUpPressed();
  void pageDownPressed(); 

private:
  virtual int inqDepth() const = 0;
  virtual float depthRatio() const = 0;
  virtual float inqRatio() const = 0;
  void drawZoomRectangle();
  void paintImages();
  void paintGraphics();  
  void drawCrossHairs(); 
  virtual void moveCursor(short dx, short dy) = 0;
  virtual void setCursorSlice(short) = 0;
  virtual void drawDtiLines() = 0;
  virtual QString inqLocationText() const = 0;

  //! @brief Renders the cross hairs.
  //! @param c the Cursor::Handle for the location at which the cross hairs will
  //!        be drawn
  //! @param slice the highlighted slice (rendered brighter than the 
  //! crosshairs of the other slices)
  virtual void drawCrossHairLines(const Cursor::Handle c,int slice) = 0;

  virtual void cursorEvent(int x, int y) = 0; 
  virtual void drawEvent(int x, int y);
  virtual void floodEvent(int x, int y);
  virtual void commitGraphics();
  virtual void loadStore();

  void undoGraphics();
  void briconEvent(int x, int y);
  void zoomEvent(int x, int y);
  void transEvent(int dx,int dy);
  void setStartMove(int x, int y) { m_startX = x; m_startY = y; }
  void setDataRect();
  void setWinRect();
  void initZoom();
  bool layerValidForDrawing();

  const QPoint  convMouseToWorld(const QPoint &) const;
  const QPoint  convWorldToViewport(const QPoint &) const;
 
  inline void invalidateImageBuffers()  const { m_imageBuffersValid  = false;invalidateDisplayBuffer();}
  inline bool imageBuffersValid()       const { return m_imageBuffersValid; }
  inline void validateImageBuffers()    const { m_imageBuffersValid = true; }
  
  inline void invalidateDisplayPixmap() const { m_displayPixmapValid = false;}
  inline void validateDisplayPixmap()   const { m_displayPixmapValid = true;}
  inline bool displayPixmapValid()      const { return m_displayPixmapValid;}

  inline void invalidateDisplayBuffer() const { m_displayBufferValid = false;invalidateDisplayPixmap();}
  inline void validateDisplayBuffer()   const { m_displayBufferValid = true;}
  inline bool displayBufferValid()      const { return m_displayBufferValid;}

  static QCursor    *m_crossCursor;
  static QCursor    *m_zoomCursor;
  static QCursor    *m_panCursor;
  static QCursor    *m_penCursor; 
  static QCursor    *m_eraserCursor;
  static QCursor    *m_fillCursor;
  std::list<Shape::Handle>     &m_undoList;
  Orientation        m_orient;  
  SliceWidget::Mode  m_mode;
  //  int                m_maskMode;
  bool               m_zooming; 
  bool               m_trueScale;
  bool               m_forceRender;
  mutable bool       m_imageBuffersValid; 
  mutable bool       m_displayPixmapValid;
  mutable bool       m_displayBufferValid;
  Shape::Handle      m_shape;
  QTimer            *m_timer;
  QPixmap           *m_displayPixmap;
};

/**
 * @brief Specialisation of SliceWidget for displaying sagittally sliced orientations
 */
class SagittalWidget : public SliceWidget
{
public:
  SagittalWidget(QWidget* parent, const char* name, 
		 Cursor::Handle c, OverlayList::Handle l,
		 DrawSettings::Handle d,
		 std::list<Shape::Handle>& u,
		 const ViewOptions& vo):
    SliceWidget(parent, name, SliceWidget::Sagittal, c, l, d, u, vo) 
    {};

  virtual void setImageCursor(int, int, int, int);
  virtual void setZoom(int);
private:  
  virtual int inqWidth() const;
  virtual int inqHeight() const;
  virtual int inqDepth() const;
  virtual float depthRatio() const;
  virtual float inqRatio() const;  
  virtual void setCursorSlice(short);  
  virtual void moveCursor(short dx, short dy);
  virtual void drawCrossHairLines(const Cursor::Handle c, int slice);
  virtual ColorRGBAHandle bufferVolume(MetaImage::Handle);  
  virtual ColorRGBAHandle dtiVolume(MetaImage::Handle);  
  virtual void drawDtiLines();
  virtual void cursorEvent(int x, int y);
  virtual QString inqLocationText() const;
};

/**
 * @brief Specialisation of SliceWidget for displaying coronally sliced orientations
 */
class CoronalWidget : public SliceWidget
{
public:
  CoronalWidget(QWidget* parent, const char* name, 
		Cursor::Handle c, OverlayList::Handle l,
		DrawSettings::Handle d,
		std::list<Shape::Handle>& u,
		const ViewOptions& vo):
    SliceWidget(parent, name, SliceWidget::Coronal, c, l, d, u, vo) 
    {};
  virtual void setImageCursor(int, int, int, int);  
  virtual void setZoom(int);
private:  
  virtual int inqWidth() const;
  virtual int inqHeight() const;
  virtual int inqDepth() const;
  virtual float depthRatio() const;
  virtual float inqRatio() const;  
  virtual void setCursorSlice(short);  
  virtual void moveCursor(short dx, short dy);
  virtual void drawCrossHairLines(const Cursor::Handle c,int slice);
  virtual ColorRGBAHandle bufferVolume(MetaImage::Handle);  
  virtual ColorRGBAHandle dtiVolume(MetaImage::Handle);  
  virtual void drawDtiLines();
  virtual void cursorEvent(int x, int y);
  virtual QString inqLocationText() const;
};

/**
 * @brief Specialisation of SliceWidget for displaying axialally sliced orientations
 */
class AxialWidget : public SliceWidget
{
public:
  AxialWidget(QWidget* parent, const char* name, 
	      Cursor::Handle c, OverlayList::Handle l, 
	      DrawSettings::Handle d,std::list<Shape::Handle>& u,
	      const ViewOptions& vo):
    SliceWidget(parent, name, SliceWidget::Axial, c, l, d, u, vo) 
    {};
  virtual void setImageCursor(int, int, int, int);  
  virtual void setZoom(int);
private:
    
  virtual int inqWidth() const;
  virtual int inqHeight() const;
  virtual int inqDepth() const;
  virtual float depthRatio() const;
  virtual float inqRatio() const;
  virtual void setCursorSlice(short);  
  virtual void moveCursor(short dx, short dy);
  virtual void drawCrossHairLines(const Cursor::Handle c, int slice);
  virtual ColorRGBAHandle bufferVolume(MetaImage::Handle); 
  virtual ColorRGBAHandle dtiVolume(MetaImage::Handle);  
  virtual void drawDtiLines();
  virtual void cursorEvent(int x, int y);
  virtual QString inqLocationText() const;
};

typedef std::vector<SliceWidget::Handle> SliceList;
typedef boost::shared_ptr<SliceList> SliceListHandle;

/* SetImageCursor
 *
 * @author David Flitney <flitney@fmrib.ox.ac.uk>
 *
 * @brief Provides access to SliceWidget::setImageCursor for use with std::for_each etc
 */
class SetImageCursor {
public:
  SetImageCursor(const Cursor::Handle c) : m_c(c) {}
  void operator()(SliceList::value_type& e) const
  { 
    e->setImageCursor(m_c->inqX(), m_c->inqY(), m_c->inqZ(), m_c->inqV()); 
  }
private:
  const Cursor::Handle m_c;
};

#endif
