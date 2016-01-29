/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

// slice.h: interface for the Slice class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SLICE_H)
#define SLICE_H

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

class SliceVisitor;

class Slice  
{
public:
	typedef boost::shared_ptr< Slice > Handle;

	Slice();
	virtual ~Slice();

	virtual short inqX() const = 0;
	virtual short inqY() const = 0;

	virtual void accept(SliceVisitor& v) = 0;
};

template <class VoxelType>
class SliceStore: public Slice
{
public:
	typedef boost::shared_ptr< SliceStore<VoxelType> > Handle;
	typedef boost::weak_ptr< SliceStore<VoxelType> > WeakHandle;

	virtual ~SliceStore();
	
	VoxelType& operator()(short x, short y);
	VoxelType& operator()(unsigned int offset);

	virtual short inqX() const;
	virtual short inqY() const;

	static Handle create(short x, short y, VoxelType *buffer);

	virtual void accept(SliceVisitor &v);

private:
	SliceStore(short x, short y, VoxelType *buffer);

	short m_x, m_y;
	VoxelType *m_buffer;
	WeakHandle m_handle;
};

typedef SliceStore<unsigned char> SliceB;
typedef SliceStore<unsigned short> SliceUS;
typedef SliceStore<float> SliceF;

#include "slice.inc"

#endif
