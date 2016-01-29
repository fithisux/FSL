/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include <qstatusbar.h>
#include <qtoolbar.h>
#include <qtoolbutton.h>
#include <qcheckbox.h>
#include <qprinter.h>
#include <qregexp.h>

#include "histogramwidget.h"
#include "histogramoptionsdialogimpl.h"
#include "qwt_plot.h"

#include <math.h>

#include <algorithm>

//#define HAVE_QWTSTDVECTORDATA

#if !defined(HAVE_QWTSTDVECTORDATA)

class QWT_EXPORT QwtStdVectorData: public QwtData
{
public:
    /*!
      Constructor
      
      \sa QwtCurve::setData and QwtPlot::setCurveData.
     */
    QwtStdVectorData(const std::vector<double> &x,
                     const std::vector<double> &y);
    QwtStdVectorData &operator=(const QwtStdVectorData &);
    virtual QwtData *copy() const;

    virtual size_t size() const;
    virtual double x(size_t i) const;
    virtual double y(size_t i) const;

    virtual QwtDoubleRect boundingRect() const;

private:
    void initCache();
    std::vector<double> d_x;
    std::vector<double> d_y;
    QwtDoubleRect d_cache;
};

QwtStdVectorData::QwtStdVectorData(
    const std::vector<double> &x,
    const std::vector<double> &y)
    : d_x(x), d_y(y)
{
    initCache();
}

QwtStdVectorData& QwtStdVectorData::operator=(const QwtStdVectorData &data)
{
    if (this != &data)
    {
        d_x = data.d_x;
        d_y = data.d_y;
        initCache();
    }
    return *this;
}

size_t QwtStdVectorData::size() const 
{ 
    return QMIN(d_x.size(), d_y.size()); 
}

double QwtStdVectorData::x(size_t i) const 
{ 
    return d_x[int(i)]; 
}

double QwtStdVectorData::y(size_t i) const 
{ 
    return d_y[int(i)]; 
}

QwtData *QwtStdVectorData::copy() const 
{ 
    return new QwtStdVectorData(d_x, d_y); 
}

/*!
  Returns the bounding rectangle of the data. If there is
  no bounding rect, like for empty data the rectangle is invalid:
  QwtDoubleRect::isValid() == FALSE
*/
QwtDoubleRect QwtStdVectorData::boundingRect() const
{
    return QwtDoubleRect(
        d_cache.x1(), d_cache.x2(), d_cache.y1(), d_cache.y2());
}

void QwtStdVectorData::initCache()
{
    const size_t sz = size();

    if ( sz <= 0 )
    {
        d_cache = QwtDoubleRect(1.0, -1.0, 1.0, -1.0); // invalid
        return;
    }

    double minX, maxX, minY, maxY;
    std::vector<double>::const_iterator xIt = d_x.begin();
    std::vector<double>::const_iterator yIt = d_y.begin();
    std::vector<double>::const_iterator end = d_x.begin() + sz;
    minX = maxX = *xIt++;
    minY = maxY = *yIt++;

    while ( xIt < end )
    {
        const double xv = *xIt++;
        if ( xv < minX )
            minX = xv;
        if ( xv > maxX )
            maxX = xv;

        const double yv = *yIt++;
        if ( yv < minY )
            minY = yv;
        if ( yv > maxY )
            maxY = yv;
    }
    
    d_cache.setRect(minX, maxX, minY, maxY);
}
#endif

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

/** Automatically calculate, min, max and bin size for histogram display
 *
 * @return The number of bins that will be need
 *
 * @param min         The minimum value in your data set
 * @param max         The maximum value in your data set
 * @param adjustedMin A new minimum which will include all your data
 * @param adjustedMax A new maximum which will include all your data
 * @param binSize     The calculated size for each bin
 * @param isInteger   Set to true if you want only integer binning
 */
unsigned int AutoBin(double min, 
		     double max, 
		     double &adjustedMin, 
		     double &adjustedMax, 
		     double &binSize,
		     bool   isInteger)
{
  double range(max - min);
  
  //
  // Find a natural bin size by rounding log10 of the range.
  // 
  binSize = pow(10, int(ceil(log10(range) - 1) - 1));

  //
  // Half the bin size if there's going to be less than 100 bins.
  // 
  float nbins(range / binSize);
  while( nbins < 100.0 )  { binSize /= 2.0; nbins = range / binSize; }

  // 
  // We don't want to use fractional binSizes unless absolutely
  // necessary so clamp the size for indicated integer images.
  // 
  if( isInteger ) binSize = std::max(1.0, ceil(binSize));

  adjustedMin = floor(min / binSize) * binSize;
  adjustedMax =  ceil(max / binSize) * binSize;

  return unsigned((adjustedMax - adjustedMin) / binSize) + 1;
}

HistogramOptions::HistogramOptions():
  m_bins(100),
  m_min(0), m_max(255),
  m_intensityRange(false),
  m_logScale(false),
  m_ignoreZeros(true),
  m_specifyBins(false)
{
}

HistogramOptions::HistogramOptions(const HistogramOptions& options):
  m_bins(options.m_bins),
  m_min(options.m_min), m_max(options.m_max),
  m_intensityRange(options.m_intensityRange),
  m_logScale(options.m_logScale),
  m_ignoreZeros(options.m_ignoreZeros),
  m_specifyBins(options.m_specifyBins)
{
}

void HistogramOptions::Swap(HistogramOptions& other)
{
  std::swap(m_bins, other.m_bins);
  std::swap(m_min, other.m_min);
  std::swap(m_max, other.m_max);
  std::swap(m_intensityRange, other.m_intensityRange);
  std::swap(m_logScale, other.m_logScale);
  std::swap(m_ignoreZeros, other.m_ignoreZeros);
  std::swap(m_specifyBins, other.m_specifyBins);
}

HistogramOptions& HistogramOptions::operator=(const HistogramOptions& rhs)
{
  HistogramOptions temp(rhs);
  Swap(temp);
  return *this;
}

class Histogram
{
public:
  Histogram(Volume::Handle v, bool isIntegerData):
    m_volume(v), m_isIntegerData(isIntegerData)
  {
    m_volume->calculateMinMax();

    m_options.setMin(m_volume->inqMin());
    m_options.setMax(m_volume->inqMax());

    calculate();
  }

  void calculate(void)
  {
    unsigned int bins(0U);

    if(!m_options.inqIntensityRange()) {
      bins = AutoBin(m_volume->inqMin(), m_volume->inqMax(), 
		     m_adjustedMin, m_adjustedMax, m_delta, m_isIntegerData);
    } else {
      bins = AutoBin(m_options.inqMin(), m_options.inqMax(),
		     m_adjustedMin, m_adjustedMax, m_delta, m_isIntegerData);
    }

    if(m_options.inqSpecifyBins())
      {
	bins = m_options.inqBins();
	m_delta = (m_adjustedMax - m_adjustedMin) / bins;
      }

    m_x.resize(bins);
    m_y.resize(bins);

    for(unsigned int n = 0; n < bins; n++)
      {
        m_x[n] = m_adjustedMin + (n * m_delta);
        m_y[n] = 0.1;
      }

    unsigned int nVoxels = m_volume->inqX() * m_volume->inqY() * m_volume->inqZ();

    for(unsigned int voxel = 0; voxel < nVoxels; voxel++)
      {
        unsigned int binNumber = (int)floor((m_volume->value(voxel) - m_adjustedMin) / m_delta);

	if((binNumber >= 0) && (binNumber < bins)) {
	  if(! (m_options.inqIgnoreZeros() && (fabs(m_volume->value(voxel)) < 0.0001)) )
	    m_y[binNumber]++;
	}
      }

    m_options.setBins(bins);
  }

  double inqYValue(double x)
  {
    x = std::max(m_adjustedMin, std::min(x, m_adjustedMax - m_delta));
    unsigned int binNumber = (int)floor((x - m_adjustedMin) / m_delta);
    return m_y[binNumber] - 0.1;
  }

  double inqXValue(double x)
  {
    x = std::max(m_adjustedMin, std::min(x, m_adjustedMax - m_delta));
    unsigned int binNumber = (int)floor((x - m_adjustedMin) / m_delta);
    return m_x[binNumber];
  }

  QwtStdVectorData inqData()
  {
    return QwtStdVectorData(m_x, m_y);
  }

  void options(QWidget* parent)
  {
    HistogramOptionsDialogImpl optionsDialog(parent, m_options);

    if(optionsDialog.exec() == QDialog::Accepted)
      {
	m_options = optionsDialog.getOptions();
	calculate();
      }
  }

  bool inqLogScale()
  {
    return m_options.inqLogScale();
  }

private:
  HistogramOptions m_options;
  Volume::Handle m_volume;

  bool m_isIntegerData;
  std::vector<double> m_x;
  std::vector<double> m_y;

  double m_adjustedMin, m_adjustedMax, m_delta;

};

class Plot: public QwtPlot
{
public:
  Plot(QWidget *parent, const QString& title): 
    QwtPlot(parent), 
    m_curve(0)
  {
    setTitle(title);
    
    setAxisTitle(xBottom, "Intensity");
    setAxisTitle(yLeft,   "#voxels");

    m_curve = insertCurve("Intensity histogram");
    setCanvasBackground(QColor(white));
    setCurvePen(m_curve, QPen(blue));

    setTitleFont(parent->font());
    setAxisTitleFont(QwtPlot::yLeft,   parent->font());
    setAxisTitleFont(QwtPlot::xBottom, parent->font());
//     enableOutline(TRUE);
//     setOutlineStyle(Qwt::VLine);
//     setOutlinePen(QPen(green));

    replot();
  }

  void redrawHistogram(const QwtStdVectorData& data, bool logScale)
  {
    setCurveData(m_curve, data);

    if(logScale)
      setAxisOptions(QwtPlot::yLeft, QwtAutoScale::Logarithmic);
    else
      setAxisOptions(QwtPlot::yLeft, QwtAutoScale::Floating);
    setAxisMargins(QwtPlot::yLeft, 0, 0);

    replot();
  }

  void showMarker(double x, double y)
  {
    removeMarkers();
    long m = insertMarker();
    QwtSymbol s;
    s.setStyle(QwtSymbol::Cross);
    s.setSize(30);
    setMarkerSymbol(m, s);
    setMarkerXPos(m, x); setMarkerYPos(m, y + 0.1);
    replot();
 }

private:
  long m_curve;
};
 
#include "histogramtoolbar.h"

/** 
 * Constructor
 * 
 * @param parent Parent widget.
 * @param vol The volume to be analysed.
 * @param isInteger If set the binning will use integer bins sizes.
 */
HistogramWidget::HistogramWidget(QWidget *parent, Volume::Handle vol, const std::string& name, unsigned int n, bool isInteger): 
  ViewWidget(parent), m_zoom(false)
{
  QString title = QString("Histogram of %1 volume %2").arg(name.c_str()).arg(n);
  m_graphWidget = new Plot(this, title);
  m_histogram = new Histogram(vol, isInteger);

  m_graphWidget->redrawHistogram(m_histogram->inqData(), m_histogram->inqLogScale());

  setCaption(title);
  setMinimumSize(400,400);
  m_graphWidget->setMargin(10);

  setCentralWidget(m_graphWidget);

  QToolBar *t = new QToolBar(this);
  HistogramToolbar *ht = new HistogramToolbar(t);

  addToolBar(t, Top, FALSE);  

  statusBar()->addWidget(new QLabel(statusBar()), 1, FALSE);

  connect(ht->m_zoomButton, SIGNAL(toggled(bool)), SLOT(toggleZoom(bool)));
  connect(ht->m_printButton, SIGNAL(clicked()), SLOT(print()));
  connect(ht->m_optionsButton, SIGNAL(clicked()), SLOT(options()));
  connect(m_graphWidget, SIGNAL(plotMousePressed(const QMouseEvent&)),
	  SLOT(plotMousePressed(const QMouseEvent&)));
  connect(m_graphWidget, SIGNAL(plotMouseMoved(const QMouseEvent&)),
	  SLOT(plotMouseMoved(const QMouseEvent&)));
  connect(m_graphWidget, SIGNAL(plotMouseReleased(const QMouseEvent&)),
	  SLOT(plotMouseReleased(const QMouseEvent&)));
}

HistogramWidget::~HistogramWidget() { delete m_histogram; } 

void HistogramWidget::showInfo(QString text)
{
  statusBar()->message(text);
}

void HistogramWidget::options()
{
  m_histogram->options(this);
  m_graphWidget->redrawHistogram(m_histogram->inqData(), m_histogram->inqLogScale());
}

void HistogramWidget::print()
{
  QPrinter printer;

  QString docName = m_graphWidget->title();
  if ( docName.isEmpty() )
    {
      docName.replace (QRegExp (QString::fromLatin1 ("\n")), tr (" -- "));
      printer.setDocName (docName);
    }

  printer.setCreator("fslview");
  printer.setOrientation(QPrinter::Landscape);

  if (printer.setup())
    m_graphWidget->print(printer, PrintFilter());
}

void HistogramWidget::toggleZoom(bool on)
{
  if (on)
    {
      m_zoom = true;
    }
  else
    {
      // Disable Zooming.
      m_graphWidget->setAxisAutoScale(QwtPlot::yLeft);
      m_graphWidget->setAxisAutoScale(QwtPlot::yRight);
      m_graphWidget->setAxisAutoScale(QwtPlot::xBottom);
      m_graphWidget->replot();
      m_zoom = false;
    }
  
//   if (m_zoom)
//     showInfo(zoomInfo);
//   else
//     showInfo(cursorInfo);
}
  
void HistogramWidget::plotMousePressed(const QMouseEvent &e)
{
  m_p1 = e.pos();

  plotMouseMoved(e);

  if (m_zoom)
    {
      m_graphWidget->enableOutline(true);
      m_graphWidget->setOutlineStyle(Qwt::Rect);
    } 
  else
    m_graphWidget->enableOutline(false);
}

void HistogramWidget::plotMouseMoved(const QMouseEvent &e)
{
  QString info;
  float x = m_graphWidget->invTransform(QwtPlot::xBottom, e.pos().x());
  
  if(!m_zoom)
    {  
      info.sprintf("Intensity %g #voxels %g",  m_histogram->inqXValue(x), m_histogram->inqYValue(x));
      m_graphWidget->showMarker(m_histogram->inqXValue(x), m_histogram->inqYValue(x));
      showInfo(info);
    }
  
}

void HistogramWidget::plotMouseReleased(const QMouseEvent &e)
{
  // some shortcuts
  int axl= QwtPlot::yLeft, axb= QwtPlot::xBottom;

  if (m_zoom)
    {
      int x1 = std::min(m_p1.x(), e.pos().x());
      int x2 = std::max(m_p1.x(), e.pos().x());
      int y1 = std::min(m_p1.y(), e.pos().y());
      int y2 = std::max(m_p1.y(), e.pos().y());
        
      // limit selected area to a minimum of 11x11 points
      int lim = 5 - (y2 - y1) / 2;
      if (lim > 0)
        {
	  y1 -= lim;
	  y2 += lim;
        }
      lim = 5 - (x2 - x1 + 1) / 2;
      if (lim > 0)
        {
	  x1 -= lim;
	  x2 += lim;
        }

      // Set fixed scales
      m_graphWidget->setAxisScale(axl, m_graphWidget->invTransform(axl,y1), 
				  m_graphWidget->invTransform(axl,y2));
      m_graphWidget->setAxisScale(axb, m_graphWidget->invTransform(axb,x1), 
				  m_graphWidget->invTransform(axb,x2));
      m_graphWidget->replot();
        
      m_graphWidget->setOutlineStyle(Qwt::Triangle);

      m_zoom = false;
    }
}


