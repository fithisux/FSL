/*  Copyright (C) 1999-2004 University of Oxford  */

/*  CCOPYRIGHT */

#include "options.h"
#include "buildno.h"

namespace Utilities {

  using namespace std;

  void OptionParser::describe_options()
  {
    for(Options::iterator option = options_.begin(); option != options_.end(); 
	option++)
      {
	if((*option)->compulsory() && (*option)->visible()) {
	  static bool banner = true;
	  if(banner) {
	    cerr << endl << "Compulsory arguments (You MUST set one or more of):" << endl;
	    banner = false;
	  }
	  (*option)->usage(cerr); cerr << endl;
	}
      }


    for(Options::iterator optionx = options_.begin(); optionx != options_.end(); 
	optionx++)
      {
	if(!(*optionx)->compulsory() && (*optionx)->visible()) {
	  static bool banner = true;
	  if(banner) {
	    cerr << endl << "Optional arguments (You may optionally specify one or more of):" << endl;
	    banner = false;
	  }
	  (*optionx)->usage(cerr); cerr << endl;
	}
      }

    cerr << endl;
    cerr << endl;
  }

  void OptionParser::brief_usage()
  {
    cerr << progname_ << endl << endl;
    cerr << "Usage: " << endl << example_ << endl;

    describe_options();
  }

  void OptionParser::usage()
  {
    cerr << endl << "Part of FSL (build " << build << ")"<< endl;
    cerr << progname_ << endl << endl;
    cerr << "Usage: " << endl << example_ << endl;

    describe_options();
  }
}
