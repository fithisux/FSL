/*  Copyright (C) 1999-2004 University of Oxford  */

/*  CCOPYRIGHT */

#include "options.h"

namespace Utilities {

  OptionParser *OptionParser::instance_ = 0;

  OptionParser *OptionParser::Instance()
  {
    if(instance_ == 0)
      instance_ = new OptionParser;
    return instance_;
  }

}
