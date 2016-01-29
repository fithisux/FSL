
template <class VoxelType>
inline short SliceStore<VoxelType>::inqX() const { return m_x; }

template <class VoxelType>
inline short SliceStore<VoxelType>::inqY() const { return m_y; }

template <class VoxelType>
inline VoxelType& SliceStore<VoxelType>::operator()(short x, short y)
{
	return m_buffer[(y * m_x) + x];
}

template <class VoxelType>
inline VoxelType& SliceStore<VoxelType>::operator()(unsigned int offset)
{
	return m_buffer[offset];
}

