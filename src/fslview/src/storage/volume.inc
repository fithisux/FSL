
template <class VoxelType>
inline short VolumeStore<VoxelType>::inqX() const { return m_x; }

template <class VoxelType>
inline short VolumeStore<VoxelType>::inqY() const { return m_y; }

template <class VoxelType>
inline short VolumeStore<VoxelType>::inqZ() const { return m_z; }

//! @brief Normalized value of volume at given location
//!
//! @param x X voxel coordinate
//! @param y Y voxel coordinate
//! @param z Z voxel coordinate
//! @return v(x, y, z) / (max - min)
template <class VoxelType>
inline float VolumeStore<VoxelType>::normalized(short x, short y, short z) const
{
	return value(x,y,z) / (float)(inqMax() - inqMin());
}

//! @brief Value of volume at given location
template <class VoxelType>
inline float VolumeStore<VoxelType>::value(unsigned int offset) const
{
  float tmp;
  if (m_doscaling) { tmp=m_slope *  m_buffer[offset] + m_intercept; }
  else { tmp= m_buffer[offset]; }
  return tmp;
}

//! @brief Value of volume at given location
//!
//! @param x X voxel coordinate
//! @param y Y voxel coordinate
//! @param z Z voxel coordinate
//! @return v(x, y, z)
template <class VoxelType>
inline float VolumeStore<VoxelType>::value(short x, short y, short z) const
{
  float tmp;
  if (m_doscaling) { tmp=m_slope * this->operator()(x,y,z) + m_intercept; }
  else { tmp=this->operator()(x,y,z); }
  return tmp;
}

//! @param x X voxel coordinate
//! @param y Y voxel coordinate
//! @param z Z voxel coordinate
//! @return v(x, y, z)
template <class VoxelType>
inline VoxelType& VolumeStore<VoxelType>::operator()(short x, short y, short z)
{
	return m_buffer[((z * m_y) + y) * m_x + x];
}

//! @param offset
//! @return v[offset]
template <class VoxelType>
inline VoxelType& VolumeStore<VoxelType>::operator()(unsigned int offset)
{
	return m_buffer[offset];
}

//! @param x X voxel coordinate
//! @param y Y voxel coordinate
//! @param z Z voxel coordinate
//! @return v(x, y, z)
template <class VoxelType>
inline const VoxelType& VolumeStore<VoxelType>::operator()(short x, short y, short z) const
{
	return m_buffer[((z * m_y) + y) * m_x + x];
}

//! @param offset
//! @return v[offset]
template <class VoxelType>
inline const VoxelType& VolumeStore<VoxelType>::operator()(unsigned int offset) const
{
	return m_buffer[offset];
}

//! @param x X voxel coordinate
//! @param y Y voxel coordinate
//! @param z Z voxel coordinate
//! @param value the value to set v(x, y, z) to be
template <class VoxelType>
inline void VolumeStore<VoxelType>::setValue(short x, short y, short z, float value)
{
  if(x >= 0 && x < m_x &&
     y >= 0 && y < m_y &&
     z >= 0 && z < m_z)
    {
      m_buffer[((z * m_y) + y) * m_x + x] = (VoxelType)value;  
    }

}
