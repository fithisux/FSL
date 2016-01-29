// testvolume.cpp: implementation of the testVolume class.
//
//////////////////////////////////////////////////////////////////////

#include <assert.h>

#include "testvolume.h"
#include "sum.h"

#include "../volume.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

class TestVolumeStoreImpl
{
public:
	TestVolumeStoreImpl() {}
	~TestVolumeStoreImpl() {}

	void testReferenceCounting();
	void testValidRead();
	void testAccept();
};

void TestVolumeStoreImpl::testReferenceCounting()
{
	Volume::Handle v1 = VolumeUS::load("structural");

	Volume::Handle v2 = v1;
	{
		Volume::Handle v3 = v1;
		Volume::Handle v4 = v2;

		unsigned int x = v1.use_count();

		assert(v1.use_count() == 4);
	}

	assert(v1.use_count() == 2);
}

void TestVolumeStoreImpl::testValidRead()
{
	Volume::Handle v = VolumeUS::load("structural");

	assert(v->inqX() == 256);
	assert(v->inqY() == 256);
	assert(v->inqZ() == 128);

//	assert((*v)(10, 20, 15) == 497);
}


void TestVolumeStoreImpl::testAccept()
{
	Volume::Handle v = VolumeUS::load("structural");
	Volume::Handle c = v;
	Sum sum;
	v->accept(sum);
	
	float total = sum.total();
	int count = sum.count();

	assert(count == v->inqX() * v->inqY() * v->inqZ());

	Volume::Handle d = c;
}

TestVolumeStore::TestVolumeStore()
{
	pImpl_ = new TestVolumeStoreImpl;
}

TestVolumeStore::~TestVolumeStore()
{
	delete pImpl_;
}

void TestVolumeStore::run()
{
//	pImpl_->testReferenceCounting();
//	pImpl_->testValidRead();
	pImpl_->testAccept();
}
