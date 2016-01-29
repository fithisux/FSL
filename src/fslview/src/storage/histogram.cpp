/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "histogram.h"

#include <algorithm>
#include <cmath>

Histogram::Histogram(Volume::Handle v, unsigned int n):
  m_volume(v), m_nBins(n), m_delta(0), m_max(0), m_min(0)
{
  m_bins.reserve(m_nBins);
  m_bins.resize(m_nBins, 0);

  float offset = m_volume->inqMin();
  float range = m_volume->inqMax() - m_volume->inqMin();
  m_delta = range / m_nBins;

  unsigned int nVoxels = m_volume->inqX() * m_volume->inqY() * m_volume->inqZ();

  for(unsigned int voxel = 0; voxel < nVoxels; voxel++)
    {
      unsigned int binNumber = (int)floor((m_volume->value(voxel) - offset) / m_delta);
      if(binNumber < 0) binNumber = 0;
      if(binNumber >= m_nBins) binNumber = m_nBins - 1;
      m_bins[binNumber]++;
    }

  m_min = *std::min_element(m_bins.begin(), m_bins.end());
  m_max = *std::max_element(m_bins.begin(), m_bins.end());
}

Histogram::~Histogram()
{
}

Histogram::Handle Histogram::getHistogram(Volume::Handle v)
{
  return Histogram::Handle(new Histogram(v, 200));
}

float Histogram::inqFrequency(unsigned int n) const { return m_bins[n]; }

unsigned int Histogram::inqBins() const { return m_nBins; }
unsigned int Histogram::inqMin() const { return m_min; }
unsigned int Histogram::inqMax() const { return m_max; }
