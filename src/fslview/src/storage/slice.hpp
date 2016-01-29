#include "slice.h"
#include "slicevisitor.h"

template <class VoxelType>
SliceStore<VoxelType>::SliceStore(short x, short y, VoxelType* buf)
: m_x(x), m_y(y), m_buffer(buf)
{
}

template <class VoxelType>
SliceStore<VoxelType>::Handle SliceStore<VoxelType>::create(short x, short y, VoxelType* buf)
{
	Handle dst(new SliceStore(x, y, buf));

	dst->m_handle = dst;

	return dst;
}

template <class VoxelType>
SliceStore<VoxelType>::~SliceStore()
{
	delete [] m_buffer;
}

template <class VoxelType> 
void SliceStore<VoxelType>::accept(SliceVisitor& v)
{
	v.visit(m_handle); 
}


