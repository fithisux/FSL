/*  Copyright (C) 1999-2004 University of Oxford  */

/*  CCOPYRIGHT */

#include "options.h"

namespace Utilities {

  using namespace std;

  bool OptionParser::check_compulsory_arguments(bool verbose)
  {
    bool okay = true;

    for(Options::iterator option = options_.begin();
	option != options_.end();
	option++) {
    
      if((*option)->compulsory() && (*option)->unset()) {
	if(okay) {
	  if(verbose) {
	    cerr << "***************************************************" << endl;
	    cerr << "The following COMPULSORY options have not been set:" << endl;
	  }
	  okay = false;
	}
	if(verbose)
	  (*option)->usage(cerr); cerr << endl;
      }
    }
    if(!okay && verbose)
      cerr << "***************************************************" << endl; 

    return okay;
  }

}
