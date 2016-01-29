#if !defined (TESTIMAGE_H)
#define TESTIMAGE_H

class TestImageImpl;

class TestImage 
{
public:
	TestImage();
	virtual ~TestImage();

	void run();

private:
	TestImageImpl *pImpl_;
};

#endif
