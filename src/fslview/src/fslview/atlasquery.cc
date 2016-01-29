
/*  FSLView - 2D/3D Interactive Image Viewer

    David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2006 University of Oxford  */

/*  CCOPYRIGHT */

#include "preferences.h"
#include "filemanager.h"
#include "storage/image.h"
#include <iostream>
#include <string>

#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;

namespace Utilities {
  extern bool string_to_T(Image::Handle& im, const string& s);
}
#include "options.h"
namespace Utilities {

  bool string_to_T(Image::Handle& im, const string& s) {
    try {
      im = Image::load(s);
    } catch (Image::Exception& e) {
      throw Utilities::X_OptionError("Loading image", e.what());
    }
    return true;
  }
}

AtlasGroup::Handle atlases;

int main(int argc, char **argv)
{
  using namespace Utilities;
  
  try {
    atlases = AtlasGroup::create();

    Option<bool> verbose(string("-V,--verbose"), false, 
			 string("switch on diagnostic messages"), 
			 false, no_argument);
    Option<bool> help(string("-h,--help"), false,
		      string("display this message"),
		      false, no_argument);
    Option<bool> dumpatlases(string("--dumpatlases"), false,
			   string("Dump a list of the available atlases"),
			   false, no_argument);
    Option<Image::Handle> mask(string("-m,--mask"), Image::Handle(),
			string("a mask image to use during structural lookups"),
			false, requires_argument);
    Option< vector<float> > coord(string("-c"), vector<float>(),
				string(""),
				false, requires_argument);
    Option<string> atlasname(string("-a,--atlas"), "",
			 string("name of atlas to use"),
			 true, requires_argument);
    string title("atlasquery (version 1.0)\n\nCopyright(c) 2005, University of Oxford\nDave Flitney");
    string usage("atlasquery [-a \"<atlasname>\"] [-m <maskimage>] [-c <X>,<Y>,<Z>]");
    
    OptionParser options(title, usage);
    options.add(verbose);
    options.add(help);
    options.add(mask);
    options.add(coord);
    options.add(atlasname);
    options.add(dumpatlases);
      
    for(int a = options.parse_command_line(argc, argv); a < argc; ++a) {
    }
    
    if(dumpatlases.value()) {
      for(AtlasGroup::ConstIterator it = atlases->begin(); it != atlases->end(); ++it)
	cout << it->second->inqName() << endl;
      return 0;
    }

    if(help.value() || !options.check_compulsory_arguments(true)) {
      options.usage();
      return 1;
    }

    if(verbose.value())
      cout << "Using atlas: " << atlasname.value() << endl;

    Atlas::Handle atlas = atlases->getAtlasByName(atlasname.value());
    if(!atlas) {
      cout << "Invalid atlas name. Try one of:" << endl;
      dumpatlases.set_T(true);
    }
    if(dumpatlases.value()) {
      for(AtlasGroup::ConstIterator it = atlases->begin(); it != atlases->end(); ++it)
	cout << it->second->inqName() << endl;
      return 0;
    }

    if(mask.set())
      {
	ImageInfo::Handle mi(mask.value()->getInfo());
	if(verbose.value())
	  cout << "Working from mask : " << mi->inqFileName() << endl;

	for(unsigned int index=0; index < atlas->inqNumLabels(); ++index) {
	  float prob(atlas->getAverageProbability(mask.value(), index));
	  string structure(atlas->inqStructureNameByIndex(index));
	  if(verbose.value())
	    cout << index << endl;
	  if(prob > 0)
	    cout << structure << ":" << prob << endl;
	}
      }
    else
      {
	if(verbose.value())
	  cout << "Working from coord: "
	       << coord.value()[0] << "," 
	       << coord.value()[1] << "," 
	       << coord.value()[2] << endl;

	cout << atlas->getDescription(coord.value()[0], coord.value()[1], coord.value()[2]) << endl;
      }

  } catch(X_OptionError& e) {
    cout << "X_OptionError: " << e.what() << endl;
    return 1;
  } catch(bad_lexical_cast& e) {
    cout << "Problem converting coordinate: " << e.what() << endl;
    return 1;
  } catch(AtlasGroup::Exception& e) {
    cout << "Opps! AtlasGroup::Exception: " << e.what() << endl;
    return 1;
  } catch(...) {
    cout << "Opps! Unknown exception! " << endl;
    return 1;
  }    
  return 0;
}
