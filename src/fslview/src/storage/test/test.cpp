#include "testvolume.h"
#include "testimage.h"
#include "../volume.h"

int main()
{
	TestVolumeStore volumeTest;
	volumeTest.run();

//  The following line should fail to compile. You aren't allowed to have raw VolumeStores!
//	VolumeStore<float> notAllowed(2, 2, 2, new float[8]);  // Error: must use a Handle!

	TestImage imageTest;
	imageTest.run();

	return 1;
}

