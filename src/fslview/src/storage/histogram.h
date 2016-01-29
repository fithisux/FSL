/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(HISTOGRAM_H)
#define HISTOGRAM_H

#include <boost/shared_ptr.hpp>
#include <vector>

#include "volume.h"

class Histogram
{
public: 
  typedef boost::shared_ptr< Histogram > Handle;

  static Handle getHistogram(Volume::Handle v);

  unsigned int inqMin() const;
  unsigned int inqMax() const;

  float inqFrequency(unsigned int n) const;
  unsigned int inqBins() const;
  float inqDelta() const { return m_delta; }

  ~Histogram();

private:
  Histogram(Volume::Handle v, unsigned int n);

  Volume::Handle m_volume;
  unsigned int m_nBins;
  float m_delta;
  unsigned int m_max;
  unsigned int m_min;
  std::vector<unsigned int> m_bins;
};

#endif
