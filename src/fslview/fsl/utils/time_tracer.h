/*  Time_Tracer.h

    Mark Woolrich and Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 1999-2010 University of Oxford  */

/*  CCOPYRIGHT  */

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

