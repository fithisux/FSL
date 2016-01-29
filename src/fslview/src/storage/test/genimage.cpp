#include "../volume.h"
#include "../image.h"
#include "../error.h"

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include <iostream>

class testImageInfo: public CppUnit::TestCase
{
  CPPUNIT_TEST_SUITE( testImageInfo );
  CPPUNIT_TEST( testReadOnly );
  CPPUNIT_TEST( testDti );
  CPPUNIT_TEST( testDimensions );
  CPPUNIT_TEST( testValidCoordinate );
  CPPUNIT_TEST( testTarnished );
  CPPUNIT_TEST( testImageName );
  CPPUNIT_TEST( testStdMat );
  CPPUNIT_TEST_SUITE_END();

public:

  void setUp()
  {
    m_i = ImageInfo::Handle(new ImageInfo(5,5,5,5,DT_SIGNED_SHORT,
					  -2.0,2.0,2.0,1.0, 10,
					  "CleanTestImage", false));
  }

  void testReadOnly()
  {
    CPPUNIT_ASSERT(m_i->inqReadOnly() == false);
    m_i->setReadOnly(true);
    CPPUNIT_ASSERT(m_i->inqReadOnly() == true);
  }

  void testDti()
  {
    ImageInfo::Handle a = ImageInfo::create(5,5,5,2,DT_SIGNED_SHORT,
					    -2.0,2.0,2.0,1.0, 10,
					    "a", false);
    ImageInfo::Handle b = ImageInfo::create(5,5,5,3,DT_SIGNED_SHORT,
					    -2.0,2.0,2.0,1.0, 10,
					    "b", true);
    ImageInfo::Handle c = ImageInfo::create(5,5,5,3,DT_SIGNED_SHORT,
					    -2.0,2.0,2.0,1.0, 10,
					    "c", false);
    CPPUNIT_ASSERT(a->isDtiCompatible() == false);
    CPPUNIT_ASSERT(a->isDtiImage()      == false);
    CPPUNIT_ASSERT(b->isDtiCompatible() == true);
    CPPUNIT_ASSERT(b->isDtiImage()      == true);
    CPPUNIT_ASSERT(c->isDtiCompatible() == true);
    CPPUNIT_ASSERT(c->isDtiImage()      == false);
  }

  void testDimensions()
  {
    CPPUNIT_ASSERT((m_i->m_x == 5) && (m_i->m_y == 5) && (m_i->m_z == 5) &&
		   (m_i->m_v == 5) && (m_i->m_dt == DT_SIGNED_SHORT));
    CPPUNIT_ASSERT((m_i->m_xDim == -2.0) && 
		   (m_i->m_yDim == 2.0) &&
		   (m_i->m_zDim == 2.0));		   
  }

  void testValidCoordinate()
  {
    CPPUNIT_ASSERT(m_i->isValidCoordinate(2, 2, 2) == true);
    CPPUNIT_ASSERT(m_i->isValidCoordinate(4, 4, 4) == true);
    CPPUNIT_ASSERT(m_i->isValidCoordinate(5, 5, 5) == false);
    CPPUNIT_ASSERT(m_i->isValidCoordinate(-2, -2, -2) == false);
  }

  void testTarnished()
  {
    CPPUNIT_ASSERT(m_i->inqTarnished() == false);
    m_i->setTarnished(true);
    CPPUNIT_ASSERT(m_i->inqTarnished() == true);
    m_i->setTarnished(false);
    CPPUNIT_ASSERT(m_i->inqTarnished() == false);
  }
    
  void testImageName()
  {
    CPPUNIT_ASSERT(m_i->inqImageName() == "CleanTestImage");
  }

  void testStdMat()
  {
    mat44 stdmat = m_i->inqStdMat();
    
    for (int i=0; i<4; i++)
      for (int j=0; j<4; j++) {
	if ((i == 0) && (j == 0))
	  CPPUNIT_ASSERT(stdmat.m[i][j]==fabs(m_i->m_xDim));
	else if ((i == 1) && (j == 1))
	  CPPUNIT_ASSERT(stdmat.m[i][j]==fabs(m_i->m_yDim));
	else if ((i == 2) && (j == 2))
	  CPPUNIT_ASSERT(stdmat.m[i][j]==fabs(m_i->m_zDim));
	else if ((i == 3) && (j == 3))
	  CPPUNIT_ASSERT(stdmat.m[i][j]==1);
	else
	  CPPUNIT_ASSERT(stdmat.m[i][j]==0);
      }
  }

private:
  ImageInfo::Handle m_test;
  ImageInfo::Handle m_i;
};

class testImage: public CppUnit::TestCase
{
  CPPUNIT_TEST_SUITE( testImage );
  CPPUNIT_TEST( testSave );
  CPPUNIT_TEST( testData );
  CPPUNIT_TEST( testMask );  
  CPPUNIT_TEST( testStat );
  CPPUNIT_TEST( testSXForm );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp()
  {
    ImageInfo::Handle info = ImageInfo::create(5,5,5,5,DT_SIGNED_SHORT,
					       -2.0,2.0,2.0,1.0, 10,
					       "StorageTestImage", false);

    m_i = Image::Handle(new Image(info));

    for(int v = 0; v < info->inqNumVolumes(); v++) 
      {

	Volume::Handle vol = m_i->getVolume(v);

	for(int z = 0; z < info->inqZ(); z++)
	  for(int y = 0; y < info->inqY(); y++)
	    for(int x = 0; x < info->inqX(); x++)
	      if(x < info->inqX()/2)
		vol->setValue(x, y, z, z+y+x+v);
	      else if(x > info->inqX()/2)
		vol->setValue(x, y, z, x+y+z-v);
	      else
		vol->setValue(x, y, z, x+y+z);
      }
  }

  void testData()
  {
    ImageInfo::Handle info(m_i->getInfo());

    for(int v = 0; v < info->inqNumVolumes(); v++) 
      {

	Volume::Handle vol = m_i->getVolume(v);

	for(short z = 0; z < info->inqZ(); z++)
	  for(short y = 0; y < info->inqY(); y++)
	    for(short x = 0; x < info->inqX(); x++)
	      if(x < info->inqX()/2)
		CPPUNIT_ASSERT(vol->value(x, y, z) == z+y+x+v);
	      else if(x > info->inqX()/2)
		CPPUNIT_ASSERT(vol->value(x, y, z) == x+y+z-v);
	      else
		CPPUNIT_ASSERT(vol->value(x, y, z) == x+y+z);
      }
  }

  void testSave()
  { 
    bool success = true;

    try {
      m_i->save("data/testimage");
    } catch (FileError &e) {
      std::cout << "FileError::" << e.inqMessage() << std::endl;
      success = false;
    }

    CPPUNIT_ASSERT( success == true );
  }

  void testMask()
  {
    m_i->save("data/mask");
    m_i->save("data/stat_mask");
    Image::Handle i = Image::load("data/mask");
    Image::Handle j = Image::load("data/stat_mask");
    
    CPPUNIT_ASSERT(!i->getInfo()->isStatImage());
    CPPUNIT_ASSERT( i->getInfo()->isMaskImage());
    CPPUNIT_ASSERT(!j->getInfo()->isStatImage());
    CPPUNIT_ASSERT( j->getInfo()->isMaskImage());
  }

  void testStat()
  {
    m_i->save("data/stat");
    Image::Handle i = Image::load("data/stat");
    
    CPPUNIT_ASSERT( i->getInfo()->isStatImage());
    CPPUNIT_ASSERT(!i->getInfo()->isMaskImage());
  }

  void testSXForm()
  {
    Image::Handle a(new Image(ImageInfo::create(5,5,5,2,DT_SIGNED_SHORT,
						-2.0,2.0,2.0,1.0, 10,
						"a", false)));
    FslSetOverrideOutputType(FSL_TYPE_ANALYZE_GZ);
    a->save("data/analyze_no_orig");
    FslSetOverrideOutputType(FSL_TYPE_NIFTI_GZ);
    a->save("data/nifti_no_orig");

    ImageInfo::Handle info(a->getInfo());

    mat44 orig = info->inqStdMat();

    orig.m[0][3] = -((5 - 1) * -2.0);
    orig.m[1][3] = -((5 - 1) * 2.0);
    orig.m[2][3] = -((5 - 1) * 2.0);

    info->setStdMat(orig);
    info->setIntent(NIFTI_XFORM_ALIGNED_ANAT);

    FslSetOverrideOutputType(FSL_TYPE_ANALYZE_GZ);
    a->save("data/analyze_aligned");
    FslSetOverrideOutputType(FSL_TYPE_NIFTI_GZ);
    a->save("data/nifti_aligned");
  }

private:
  Image::Handle m_i;
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(testImageInfo, "testImageInfo");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(testImage, "testImage");

CppUnit::Test *suite(const std::string& name)
{
  CppUnit::TestFactoryRegistry &reg = 
    CppUnit::TestFactoryRegistry::getRegistry();

  reg.registerFactory(&CppUnit::TestFactoryRegistry::getRegistry(name));

  return reg.makeTest();
}

int main(int argc, char **argv)
{
  CppUnit::TextUi::TestRunner runner;

  runner.addTest( suite("testImageInfo") );
  runner.addTest( suite("testImage") );
  
  bool wasSuccessful = runner.run( "testImageInfo", false, true );
  wasSuccessful = wasSuccessful && runner.run( "testImage", false, true );

  return wasSuccessful ? 0 : 1;
}
