#include "timeseries.h"
#include "fslio/fslio.h"

template <class VoxelType>
TimeSeriesStore<VoxelType>::TimeSeriesStore(VoxelType* buf)
  :  m_buffer(buf), m_doscaling(false), m_slope(1.0), m_intercept(0.0)
{
  m_count = 0;
}

template <class VoxelType>
TimeSeriesStore<VoxelType>::~TimeSeriesStore()
{
  delete [] m_buffer;
}

template <class VoxelType>
typename TimeSeriesStore<VoxelType>::Handle TimeSeriesStore<VoxelType>::create(short x, short y, short z, int n)
{
  VoxelType* buf = new VoxelType[n];

  Handle dst(new TimeSeriesStore(buf));
  dst->setMin(0);
  dst->setMax(0);  
  dst->setVolCount(n);
  dst->calculateMinVal();
  dst->calculateMaxVal();
  dst->setCoordinates(x,y,z);

  return dst;
}

template <class VoxelType>
typename TimeSeriesStore<VoxelType>::Handle TimeSeriesStore<VoxelType>::getTimeSeries(FSLIO* avw, short x, short y, short z)
{
  short xDim, yDim, zDim, vDim;
  FslGetDim(avw, &xDim, &yDim, &zDim, &vDim);
  VoxelType* buf = new VoxelType[vDim];
  FslSeekVolume(avw, 0);

  short storedX(x);
  if(FslGetLeftRightOrder(avw) != FSL_RADIOLOGICAL)
    storedX = xDim-1-storedX;
  FslReadTimeSeries(avw, buf, storedX, y, z, vDim);


  Handle dst(new TimeSeriesStore(buf));

  dst->m_doscaling = FslGetIntensityScaling(avw,&(dst->m_slope),&(dst->m_intercept));

  dst->setMin(VoxelType(0));
  dst->setMax(VoxelType(0));  
  dst->setVolCount(vDim);
  dst->calculateMinVal();
  dst->calculateMaxVal();
  dst->setCoordinates(x,y,z);

  return dst;
}

template <class VoxelType>
void TimeSeriesStore<VoxelType>::calculateMaxVal()
{
  m_maxVal = *std::max_element(m_buffer,m_buffer + m_count);
  if (m_doscaling) {
    m_maxVal = m_slope * m_maxVal + m_intercept;
  }
}

template <class VoxelType>
void TimeSeriesStore<VoxelType>::calculateMinVal()
{
  m_minVal = *std::min_element(m_buffer,m_buffer + m_count);
  if (m_doscaling) {
    m_minVal = m_slope * m_minVal + m_intercept;
  }
}

