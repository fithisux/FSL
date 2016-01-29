
#include "volume.h"
#include "histogramfns.hpp"
#include "fslio/fslio.h"
#include <algorithm>
#include <stdio.h>

template <class VoxelType>
VolumeStore<VoxelType>::VolumeStore(short x, short y, short z, VoxelType* buf)
  : m_x(x), m_y(y), m_z(z), m_doscaling(false), m_slope(1.0), m_intercept(0.0), m_buffer(buf)
{
}

template <class VoxelType>
VolumeStore<VoxelType>::~VolumeStore()
{
  delete [] m_buffer;
}

template <class VoxelType>
typename VolumeStore<VoxelType>::Handle VolumeStore<VoxelType>::getVolume(FSLIO* avw, float min, float max, size_t n)
{
  short xDim, yDim, zDim, vDim;
  FslGetDim(avw, &xDim, &yDim, &zDim, &vDim);

  VoxelType* buf = new VoxelType[xDim * yDim * zDim];

  FslSeekVolume(avw, n);
  FslReadVolumes(avw, buf, 1);

  if(FslGetLeftRightOrder(avw) != FSL_RADIOLOGICAL) 
    {
      for(short z = 0; z < zDim; z++)
	for(short y = 0; y < yDim; y++)
	  for(short x = 0; x < xDim/2; x++)
	    std::swap(buf[((z * yDim) + y) * xDim + x], 
		      buf[((z * yDim) + y) * xDim + xDim - 1 - x]);

    }

  Handle dst(new VolumeStore(xDim, yDim, zDim, buf));

  dst->m_doscaling = FslGetIntensityScaling(avw,&(dst->m_slope),&(dst->m_intercept));
  
  if((max - min) != 0.0) {
    dst->setMin(min);
    dst->setMax(max);
  } else
    dst->calculateRobustMinMax(false);

  return dst;
}

template <class VoxelType>
bool VolumeStore<VoxelType>::saveVolume(FSLIO* avw)
{
  if(FslGetLeftRightOrder(avw) != FSL_RADIOLOGICAL) 
  {
    for(short z = 0; z < inqZ(); z++)
      for(short y = 0; y < inqY(); y++)
	for(short x = 0; x < inqX()/2; x++)
	  std::swap(m_buffer[((z * inqY()) + y) * inqX() + x], 
		    m_buffer[((z * inqY()) + y) * inqX() + inqX() - 1 - x]);

  }
  return FslWriteVolumes(avw, m_buffer, 1);
  if(FslGetLeftRightOrder(avw) != FSL_RADIOLOGICAL) 
  {
    for(short z = 0; z < inqZ(); z++)
      for(short y = 0; y < inqY(); y++)
	for(short x = 0; x < inqX()/2; x++)
	  std::swap(m_buffer[((z * inqY()) + y) * inqX() + x], 
		    m_buffer[((z * inqY()) + y) * inqX() + inqX() - 1 - x]);

  }
}

template <class VoxelType>
bool VolumeStore<VoxelType>::inRange(short x, short y, short z) const
{
  return ( (x >= 0) && (x < inqX()) &&
	   (y >= 0) && (y < inqY()) &&
	   (z >= 0) && (z < inqZ()) );
}

//! @brief Clone a volume
template <class VoxelType>
typename VolumeStore<VoxelType>::Handle VolumeStore<VoxelType>::clone(const Handle& src)
{
  VoxelType* buf = new VoxelType[src->inqX() * src->inqY() * src->inqZ() * sizeof(VoxelType)];
	
  // Copy src's buffer here
  Handle dst(new VolumeStore(src->inqX(), src->inqY(), src->inqZ(), buf));

  return dst;
}

template <class VoxelType>
typename VolumeStore<VoxelType>::Handle VolumeStore<VoxelType>::create(short x, short y, short z, VoxelType* buf)
{
  Handle dst(new VolumeStore(x, y, z, buf));

  return dst;
}

template <class VoxelType>
float VolumeStore<VoxelType>::minValue() const
{
  VoxelType min = *std::min_element(&m_buffer[0], &m_buffer[m_x * m_y * m_z]);
   // scale these values for now (want to redo this as robust in the future)
  if (m_doscaling) 
    min = VoxelType(m_slope * min + m_intercept);
  return float(min);
}

template <class VoxelType>
float VolumeStore<VoxelType>::maxValue() const
{
  VoxelType max = *std::max_element(&m_buffer[0], &m_buffer[m_x * m_y * m_z]);
  // scale these values for now (want to redo this as robust in the future)
  if (m_doscaling)
    max = VoxelType(m_slope * max + m_intercept);
  return float(max);
}

template <class VoxelType>
void VolumeStore<VoxelType>::calculateMinMax()
{
  m_min = *std::min_element(&m_buffer[0], &m_buffer[m_x * m_y * m_z]);
  m_max = *std::max_element(&m_buffer[0], &m_buffer[m_x * m_y * m_z]);
  // scale these values for now (want to redo this as robust in the future)
  if (m_doscaling) { 
    m_min = m_slope * m_min + m_intercept;
    m_max = m_slope * m_max + m_intercept;
  }
}

//! @brief Calculate robust min and max of the volume
template <class VoxelType>
void VolumeStore<VoxelType>::calculateRobustMinMax(const bool ignoreZeros)
{
  unsigned int n(m_x * m_y * m_z);

  std::valarray<float> data(n);

  for(unsigned int i = 0; i < n; ++i)
    data[i] = value(i);
 
  std::valarray<int> histogram;
  
  //  find_thresholds(std::valarray<float>(data[data != float(0)]), histogram, 1000, m_min, m_max); 
  m_min=data.min(); m_max=data.max();

  if(ignoreZeros)
    find_thresholds(std::valarray<float>(data[data != float(0)]), histogram, 1000, m_min, m_max); 
  else
    find_thresholds(data, histogram, 1000, m_min, m_max); 
   
  if (m_doscaling) { 
    m_min = m_slope * m_min + m_intercept;
    m_max = m_slope * m_max + m_intercept;
  }
  
  if(minValue() == 0)
    m_min = 0;
  if(maxValue() == 0)
    m_max = 0;
}

//! @brief Create a blank image
//! @param x X dimension
//! @param y Y dimension
//! @param z Z dimension
template <class VoxelType>
typename VolumeStore<VoxelType>::Handle VolumeStore<VoxelType>::createBlank(short x,short y,short z)
{
  VoxelType* buf = new VoxelType[x * y * z];


  Handle dst(new VolumeStore(x, y, z, buf));

  unsigned int width  = dst->inqX();
  unsigned int height = dst->inqY();
  unsigned int depth  = dst->inqZ();

  for(unsigned int zl = 0; zl< depth; ++zl)  
    {
      for(unsigned int yl = 0; yl < height; ++yl)
        {            
        for(unsigned int xl = 0; xl < width; ++xl) 
          {
            (*dst)(xl,yl,zl) = 0;
          }
        }
    }

  int min, max;
  min = 0;
  max = 4;
  dst->setMin(min);
  dst->setMax(max);

  return dst;
}



