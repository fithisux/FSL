// testimage.cpp: implementation of the TestImage class.
//
//////////////////////////////////////////////////////////////////////

#include <assert.h>
#include <math.h>

#include "testimage.h"

#include "../image.h"
#include "../volume.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

#include "utils/options.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

class TestImageImpl
{
public:
	TestImageImpl() {}
	~TestImageImpl() {}

	void test();
};

TestImage::TestImage()
{
	pImpl_ = new TestImageImpl;
}

TestImage::~TestImage()
{
	delete pImpl_;
}

void TestImage::run()
{
	pImpl_->test();
}

Utilities::Option<int> volume(string("-v,--volume"), 1,
			      "Volume number",
			      false, Utilities::requires_argument);
Utilities::Option<bool> assertion(string("-a"), false,
			      "Assert expected result - use during automated tests",
				  false);
Utilities::Option<string> imagename(string("-i,--image"), "testimage",
				    "Image name",
				    false, Utilities::requires_argument);
Utilities::Option<string> cursorfile(string("-c,--cursor"), "testimage.txt",
				       "List of cursor locations to look up",
				       false, Utilities::requires_argument);

ostream& operator<<(ostream &os, vector<int>& iv);

istream& operator>>(istream &is, vector<int>& iv)
{
  int v, x, y, z;

  is >> v; is >> x; is >> y; is >> z;

  iv.clear();
  iv.push_back(v);
  iv.push_back(x);
  iv.push_back(y);
  iv.push_back(z);


  return is;
}

ostream& operator<<(ostream &os, vector<int>& iv)
{
  os << "Volume " << setw(3) << iv.at(0) << 
    " at (" << setw(3) << iv.at(1) << "," << setw(3) << iv.at(2) << "," << setw(3) << iv.at(3) << ")";

  return os;
}

void TestImageImpl::test()
{
  Image::Handle im = Image::load(imagename.value());

  int vol = volume.value();

  ifstream f;
  f.open(cursorfile.value().c_str());

  vector<int> loc;

  f >> loc;

  do  {
    float result;
    f >> result;

    float v = im->getVolume(loc.at(0))->value(loc.at(1), loc.at(2), loc.at(3));

    cout << "should be " << result <<  " got " << v << endl;
    if(assertion.value())
      assert( (fabs(v - result) < 1e-6) );

    f >> loc;
  } while( !f.eof() );
}

string title = 
"testimage (Version 1.0)\n\n\
Copyright(c) 2000, University of Oxford\n\
Dave Flitney";

string examples = "testimage --cursor=<number>,<number>,<number> --volume=<number> --image=<filename>";

int main(int argc, char **argv)
{
  Utilities::OptionParser options(title, examples);

  options.add(volume);
  options.add(assertion);
  options.add(imagename);
  options.add(cursorfile);

  for(unsigned int a = options.parse_command_line(argc, argv); 
      a < argc; ) ;

  TestImage t;

  t.run();
}
