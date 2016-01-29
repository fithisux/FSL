/*  log.cc

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

#include "log.h"

namespace Utilities {

  Log* LogSingleton::logger = NULL;
  int LogSingleton::count = 0;

  void Log::makeDir(const string& pdirname, const string& plogfilename, bool pstream_to_logfile, bool pstream_to_cout) 
    {
      if(logEstablished)
	{
	  logfileout.close();
	}

      dir = pdirname;
      logfilename = plogfilename;
      stream_to_logfile = pstream_to_logfile;
      stream_to_cout = pstream_to_cout;

      // make directory to place results into:
      // keep adding "+" until directory is made:      
      int count = 0;
      while(true)
	{
	  if(count >= 20)
	    {
	      string s("Cannot create directory " + dir);
	      throw Exception(s.c_str());
	    }

	  int ret = system(("mkdir "+ dir + " 2>/dev/null").c_str());
	  if(ret == 0)
	    {
	      break;
	    }
	  dir = dir + "+";
	  count++;
	}
      
      // setup logfile
      if(stream_to_logfile)
	{
	  logfileout.open((dir + "/" + logfilename).c_str(), ios::app);
	  if(logfileout.bad())
	    {
	      throw Exception(string(string("Unable to setup logfile ")+logfilename+string(" in directory ")+dir).c_str());
	    }
	}

      logEstablished = true;
    }

  void Log::setDir(const string& pdirname, const string& plogfilename, bool pstream_to_logfile, bool pstream_to_cout, ios_base::openmode mode) 
    {
      if(logEstablished)
	{
	  logfileout.close();
	}

      dir = pdirname;
      logfilename = plogfilename;
      stream_to_logfile = pstream_to_logfile;
      stream_to_cout = pstream_to_cout;

      // setup logfile
      if(stream_to_logfile)
	{
	  logfileout.open((dir + "/" + logfilename).c_str(), mode);
	   if(logfileout.bad())
	     {
	      throw Exception(string(string("Unable to setup logfile ")+logfilename+string(" in directory ")+dir).c_str());
	    }
	}
      
      logEstablished = true;
    }

  void Log::setthenmakeDir(const string& pdirname, const string& plogfilename, bool pstream_to_logfile, bool pstream_to_cout) 
    {
      if(logEstablished)
	{
	  logfileout.close();
	}

      dir = pdirname;
      logfilename = plogfilename;
      stream_to_logfile = pstream_to_logfile;
      stream_to_cout = pstream_to_cout;

      // make directory
      int ret = system(("mkdir -p "+ dir + " 2>/dev/null").c_str());
      if(ret == -1)
	{
	  throw Exception(string(string("Unable to make directory ")+dir).c_str());
	}

      // setup logfile
      if(stream_to_logfile)
	{
	  logfileout.open((dir + "/" + logfilename).c_str(), ios::app);
	   if(logfileout.bad())
	     {
	      throw Exception(string(string("Unable to setup logfile ")+logfilename+string(" in directory ")+dir).c_str());
	    }
	}
      
      logEstablished = true;
    }

}






