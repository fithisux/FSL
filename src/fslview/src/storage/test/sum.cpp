// sum.cpp: implementation of the Sum class.
//
//////////////////////////////////////////////////////////////////////

#include "sum.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Sum::Sum() : total_(0), count_(0)
{
}

template <typename VoxelType>
void Sum::calculate(VolumeStore<VoxelType>::WeakHandle target)
{
	total_ = 0;
	count_ = 0;
	
	for(unsigned short z = 0; z < target->inqZ(); ++z) 
	{
		unsigned int offset = z * target->inqY();

		for(unsigned short y = 0; y < target->inqY(); ++y)
		{
			unsigned int xy = (offset + y) * target->inqX();

			for(unsigned short x = 0; x < target->inqX(); x++)
			{
				unsigned int xyz = xy + x;

				total_ += (*target)(xyz);
				++count_;
			}
		}
	}
}

