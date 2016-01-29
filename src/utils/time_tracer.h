/*  Time_Tracer.h

    Mark Woolrich and Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 1999-2010 University of Oxford  */

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk
    
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
    
    LICENCE
    
    FMRIB Software Library, Release 5.0 (c) 2012, The University of
    Oxford (the "Software")
    
    The Software remains the property of the University of Oxford ("the
    University").
    
    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.
    
    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.
    
    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.
    
    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. If you are
    interested in using the Software commercially, please contact Isis
    Innovation Limited ("Isis"), the technology transfer company of the
    University, to negotiate a licence. Contact details are:
    innovation@isis.ox.ac.uk quoting reference DE/9564. */

#if !defined(Time_Tracer_h)
#define Time_Tracer_h

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <time.h>
#include <set>
#include <stack>
#include <iterator>

using namespace std;

namespace Utilities{

  class TimingFunction
    {
    public:
      TimingFunction(const char * pstr):
	str(pstr),
	time_taken(0),
	times_called(0)
	{}

      class comparer_name
	{
	public:
	  bool operator()(const TimingFunction* t1, const TimingFunction* t2) const
	    {
	      return strcmp(t1->str, t2->str) < 0;
	    }
	};

      class comparer_time_taken
	{
	public:
	  bool operator()(const TimingFunction* t1, const TimingFunction* t2) const
	    {
	      return t1->time_taken > t2->time_taken;
	    }
	};

      void start() {start_time = clock();}
      void end() {time_taken += clock()-start_time; times_called++;}

      friend class comparer_name;
      friend class comparer_time_taken;
      friend std::ostream& operator<<(std::ostream& ostr, const TimingFunction* t);

    protected:
      const char* str;
      clock_t time_taken;
      int times_called;
      clock_t start_time;

    private:
      TimingFunction();
      const TimingFunction& operator=(TimingFunction&);
      TimingFunction(TimingFunction&);
    };

  inline std::ostream& operator<<(std::ostream& ostr, const TimingFunction* t)
    {
      ostr << "<tr><td>" << t->str;
      ostr.setf(std::ios::fmtflags(0),ios::floatfield);
      ostr << "<td align=center>" << float(t->time_taken)/CLOCKS_PER_SEC;
      ostr.setf(ios::scientific, ios::floatfield);
      ostr <<  "<td align=center>" << t->times_called <<  "<td align=center>" << (t->time_taken/float(t->times_called))/CLOCKS_PER_SEC;
      ostr << "</tr>";
      return ostr;
    }

  // Non Newmat Tracer:
  class Time_Tracer
    {
    public:
      Time_Tracer(const char* str)
	{
	  construct(str);
	}

      Time_Tracer(char* str)
	{		  	  
	  construct(str);
	}

      void construct(const char* str)
	{
	  if(instantstack || runningstack)
	    {
	      stk.push(string(str));

	      if(runningstack)
		{
		  tmp = "";
		  pad++;
		  for(unsigned int i = 0; i < pad; i++)
		    tmp = tmp + "  ";
		  
		  std::cout << tmp << str << std::endl;
		}
	    }
	  if(timingon)
	    {
	      // see if already in list:
	      timingFunction = new TimingFunction(str);
	      set<TimingFunction*, TimingFunction::comparer_name>::iterator it = timingFunctions.find(timingFunction);
	      if(it== timingFunctions.end())
		{		  
		  timingFunctions.insert(timingFunction);
		}
	      else
		{
		  delete timingFunction;
		  timingFunction = *it;
		}
		
	      timingFunction->start();
	    }
	}

      virtual ~Time_Tracer() 
	{ 
	  if(instantstack)
	    {
	      stk.pop();
	    }

	  if(runningstack && pad > 0) 
	    {
		  std::cout << tmp << "finished" << std::endl;
	      pad--;
	    }
	  if(timingon)
	    {
	      timingFunction->end();
	    }
	  
	}

      static void dump_times(const string& dir)
	{
	  multiset<TimingFunction*, TimingFunction::comparer_time_taken> timingFunctionsByTimeTaken(timingFunctions.begin(), timingFunctions.end());
	  //copy(timingFunctions.begin(), timingFunctions.end(), timingFunctionsByTimeTaken.begin());

	  ofstream out;
	  out.open((dir + "/timings.html").c_str(), ios::out);	  
	  out << "<HTML><TITLE>Tracer Timings</TITLE><BODY><table border=3 cellspacing=5>" << endl;
	  out << "<tr><td>Function<td align=center>Total Time(secs)<td align=center>Num of calls<td align=center>Time per call(secs)</tr>" << endl;	  
	  copy(timingFunctionsByTimeTaken.begin(), timingFunctionsByTimeTaken.end(), ostream_iterator<TimingFunction*>(out, "\n"));	
	  out << "</table></BODY></HTML>" << endl;
	  out.close();
	}

      static void dump_instant_stack()
	{
	  // tmp stack to put values into for restoring stack after outputting
	  stack<string> tmpstk;

	  while(!stk.empty())
	    {
	  
		  std::cout << stk.top() << std::endl;
	      tmpstk.push(stk.top());
	      stk.pop();
	    }

	  while(!tmpstk.empty())
	    {
	      stk.push(tmpstk.top());
	      tmpstk.pop();
	    }
	}

      static void setinstantstackon() {instantstack = true;}
      static void setrunningstackon() {runningstack = true;}
      static void settimingon() {timingon = true;}

    protected:
      static bool instantstack;
      static bool runningstack;
      static bool timingon;
      static unsigned int pad;
      static set<TimingFunction*, TimingFunction::comparer_name> timingFunctions;
      static stack<string> stk;

      string tmp;
      TimingFunction* timingFunction;

    private:
      Time_Tracer();
      const Time_Tracer& operator=(Time_Tracer&);
      Time_Tracer(Time_Tracer&);
    };

}
#endif

