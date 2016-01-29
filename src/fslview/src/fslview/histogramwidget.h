/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(HISTOGRAMWIDGET_H)
#define HISTOGRAMWIDGET_H

#include "viewwidget.h"

class QwtPlot;
class Plot;
class QString;
class QPoint;

class Histogram;

class HistogramOptions
{
public:
  HistogramOptions();
  HistogramOptions(const HistogramOptions& options);
  HistogramOptions& operator=(const HistogramOptions& rhs);

  void Swap(HistogramOptions& other);

  unsigned int inqBins() const   { return m_bins; } 
  bool inqIntensityRange() const { return m_intensityRange; } 
  bool inqLogScale() const       { return m_logScale; } 
  bool inqIgnoreZeros() const    { return m_ignoreZeros; } 
  bool inqSpecifyBins() const    { return m_specifyBins; } 
  double inqMin() const          { return m_min; }
  double inqMax() const          { return m_max; }

  void setBins(unsigned int i)    { m_bins = i; } 
  void setIntensityRange(bool on) { m_intensityRange = on; } 
  void setLogScale(bool on)       { m_logScale = on; } 
  void setIgnoreZeros(bool on)    { m_ignoreZeros = on; } 
  void setSpecifyBins(bool on)    { m_specifyBins = on; } 
  void setMin(double v)           { m_min = v; }
  void setMax(double v)           { m_max = v; }

private:
  unsigned int m_bins;
  double       m_min;
  double       m_max;
  bool         m_intensityRange;
  bool         m_logScale;
  bool         m_ignoreZeros;
  bool         m_specifyBins;
};

/**
 * @author David Flitney <flitney@fmrib.ox.ac.uk>
 * @date   Thu Jan  2 14:36:21 2003
 * 
 * @brief  Customises @ref ViewWidget for displaying a @ref Histogram graph.
 * 
 */
class HistogramWidget : public ViewWidget
{
  Q_OBJECT
public:
  HistogramWidget(QWidget *parent, Volume::Handle v, const std::string& name, unsigned int n, bool isInteger);
  virtual ~HistogramWidget();

private slots:
  void print();
  void options();

  void toggleZoom(bool);
  
  void plotMousePressed(const QMouseEvent &e);
  void plotMouseReleased(const QMouseEvent &e);
  void plotMouseMoved(const QMouseEvent &e);

private:
  void showInfo(QString text);

  Histogram *m_histogram;

  Plot*  m_graphWidget;
  bool   m_zoom;
  QPoint m_p1;
};

#endif
