/*  Copyright (C) 1999-2004 University of Oxford  */

/*  CCOPYRIGHT */

#include "options.h"

namespace Utilities {

  using namespace std;

  bool BaseOption::matches(const string& arg)
  {
    string::size_type pos = 0, np;
    while((np = key_.find(",", pos)) != string::npos) {
      if(arg == key_.substr(pos, np - pos))
	return true;
      pos = np + 1;
    }
    if(arg == key_.substr(pos, string::npos))
      return true;
    return false;
  }

}
