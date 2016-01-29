

template <class VoxelType>
inline float TimeSeriesStore<VoxelType>::normalized(short nVol)
{
  return value(nVol) / (float)(inqMax() - inqMin());
}

template <class VoxelType>
inline float TimeSeriesStore<VoxelType>::value(short nVol)
{
  float tmp;
  if (m_doscaling) { tmp=m_slope * this->operator()(nVol) + m_intercept; }
  else { tmp=this->operator()(nVol); }
  return tmp;
}

template <class VoxelType>
inline void TimeSeriesStore<VoxelType>::setValue(int index, float val)
{
  VoxelType tmp;
  if(m_doscaling) 
    tmp=VoxelType((val-m_intercept) / m_slope);
  else 
    tmp=VoxelType(val);
  m_buffer[index] = tmp;
}

template <class VoxelType>
inline float TimeSeriesStore<VoxelType>::mean() const
{
  float mean;
  double accm(0.0);
  for (int n=0;n <m_count;++n)
      accm = accm + m_buffer[n];
  
  mean = accm/m_count;
  return mean;
}

template <class VoxelType>
inline int TimeSeriesStore<VoxelType>::inqVolCount()
{
  return m_count;
}


template <class VoxelType>
inline VoxelType& TimeSeriesStore<VoxelType>::operator()(short nVol)
{
	return m_buffer[nVol];
}

template <class VoxelType>
inline float TimeSeriesStore<VoxelType>::inqMaxVal()
{
  return	m_maxVal;
}

template <class VoxelType>
inline float TimeSeriesStore<VoxelType>::inqMinVal()
{
  return	m_minVal;
}
