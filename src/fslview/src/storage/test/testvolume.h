#if !defined (TESTVOLUME_H)
#define TESTVOLUME_H

class TestVolumeStoreImpl;

class TestVolumeStore  
{
public:
	TestVolumeStore();
	~TestVolumeStore();

	void run();

private:
	TestVolumeStoreImpl *pImpl_;
};

#endif