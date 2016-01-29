/*  Time_Tracer.cc

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

#include "time_tracer.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <set>

namespace Utilities {

  bool Time_Tracer::instantstack = false;
  bool Time_Tracer::runningstack = false;
  bool Time_Tracer::timingon = false;
  unsigned int Time_Tracer::pad = 0;
  set<TimingFunction*, TimingFunction::comparer_name> Time_Tracer::timingFunctions;
  stack<string> Time_Tracer::stk;
}






