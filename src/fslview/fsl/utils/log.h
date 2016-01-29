/*  log.h

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/* The Log class allows for instantiation of more than one Log
either sharing directories or not. However, Logs can not share log files.
Or you can work with the LogSIngleton class.

A Log can open new logfiles in the same log directory or start on an
entirely new directory. You can stream directly to a Log with flags   
determining streaming to the Logfile and/or cout. */

/*  CCOPYRIGHT  */

#if !defined(log_h)
#define log_h

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>
#include "newmatap.h"
#include "newmatio.h"

using namespace std;

// for the Exception class:
using namespace NEWMAT;

namespace Utilities{

  template<class t> string tostring(const t obj)
  {
    ostringstream str;
    str << obj;
    return str.str();
  }

  class Log
    {
    public:
      Log():logEstablished(false) {}

      Log(const string& pdirname, const string& plogfilename = "logfile", bool pstream_to_logfile = true, bool pstream_to_cout = false, bool makedir = true):logEstablished(false) 
	{
	  if(makedir)
	    makeDir(pdirname, plogfilename, pstream_to_logfile, pstream_to_cout);
	  else
	    setDir(pdirname, plogfilename, pstream_to_logfile, pstream_to_cout);
	}

      ~Log() { logfileout.close(); }

      /** Need to call makeDir or setDir before Log can be used */

      /** Makes a directory to place results into:
	  keeps adding "+" to pdirname until unique directory is made. */
      /** The stream_to* variables define the streaming behaviour */
      void makeDir(const string& pdirname, const string& plogfilename = "logfile", bool pstream_to_logfile = true, bool pstream_to_cout = false);

      /** Sets an existing directory to place results into. */
      /** The stream_to* variables define the streaming behaviour */
      void setDir(const string& pdirname, const string& plogfilename = "logfile", bool pstream_to_logfile = true, bool pstream_to_cout = false, ios_base::openmode mode=ios::app); 

      /** Sets an existing directory to place results into. */
      /** If does not exist then makes it. */
      /** The stream_to* variables define the streaming behaviour */
      void setthenmakeDir(const string& pdirname, const string& plogfilename = "logfile", bool pstream_to_logfile = true, bool pstream_to_cout = false);

      /** Closes old logfile buffer and attempts to open new one with name specified and sets streaming to logfile on */
      void setLogFile(const string& plogfilename, ios_base::openmode mode=ios::app);

      const string& getDir() const { if(!logEstablished)throw Exception("Log not setup");return dir; }

      const string& getLogFileName() const { if(!logEstablished)throw Exception("Log not setup");return logfilename; }

      /** returns passed in filename appended onto the end of the dir name */
      const string appendDir(const string& filename) const;

      ofstream& get_logfile_ofstream() { return  logfileout;}     

      inline void flush() { 
	if(stream_to_logfile)
	  logfileout.flush();
	
	if(stream_to_cout)
	  cout.flush();
      } 

      /** allows streaming into cout and/or logfile depending upon the */
      /** stream_to_cout and stream_to_logfile respectively */
      /** use like a normal ostream, e.g. log.str() << "hello" << endl */
      /** NOTE: can simply stream straight to Log instead, e.g. log << "hello" << endl */
      Log& str();
            
      /** sets whether or not you stream to cout */
      void set_stream_to_cout(bool in = true) { stream_to_cout = in; }

      /** sets whether or not you stream to logfile */
      void set_stream_to_logfile(bool in = true) { 
	if(!stream_to_logfile && in)
	  {
	    if(logfileout.bad())
	      {
		cerr << "Warning: Unable to stream to logfile " << logfilename << ". Need to have called log.setLogFile. Therefore, no streaming to logfile will be performed" << endl;
	      }
	  }
	else stream_to_logfile = in; 
      }

    private:      
      
      const Log& operator=(Log&);
      Log(Log&);
      
      string dir;
      ofstream logfileout;
      string logfilename;

      bool logEstablished;

      bool stream_to_logfile;
      bool stream_to_cout;

      friend Log& operator<<(Log& log, ostream& (*obj) (ostream &));

      template<class t> 
	friend Log& operator<<(Log& log, const t& obj); 

      template<class t> 
	friend Log& operator<<(Log& log, t& obj); 
    };

  template<class t> Log& operator<<(Log& log, const t& obj)
    {
      if(log.stream_to_logfile)
	log.logfileout << obj;

      if(log.stream_to_cout)
	cout << obj;

      return log;
    }

  template<class t> Log& operator<<(Log& log, t& obj)
    {
      if(log.stream_to_logfile)
	log.logfileout << obj;

      if(log.stream_to_cout)
	cout << obj;

      return log;
    }
   
  class LogSingleton
    {
    public:

      static Log& getInstance();

      ~LogSingleton() { delete logger; }      

      /** hacked in utility provides a global counter for general use: */
      static int counter() { return count++; }

    private:
      LogSingleton() {}
      
      const LogSingleton& operator=(LogSingleton&);
      LogSingleton(LogSingleton&);
      
      static Log* logger;
      
      static int count;
      
    };

  inline Log& LogSingleton::getInstance(){
    if(logger == NULL)
      logger = new Log();
  
    return *logger;
  }
  
  inline void Log::setLogFile(const string& plogfilename, ios_base::openmode mode) {

    if(!logEstablished)
      {
	throw Exception("Log not setup");
      }

    logfileout.close();

    logfilename = plogfilename;
    
    // setup logfile
    logfileout.open((dir + "/" + logfilename).c_str(), mode);
    if(logfileout.bad())
      {
	throw Exception(string(string("Unable to setup logfile ")+logfilename+string(" in directory ")+dir).c_str());
      }
    
    stream_to_logfile = true;

    logEstablished = true;
  }

  inline Log& Log::str() { 
    
    if(!logEstablished)
      {
	throw Exception("Log not setup");
      }

    return *this; 
  }
 
  inline const string Log::appendDir(const string& filename) const { 

    if(!logEstablished)
      {
	throw Exception("Log not setup");
      }
  
    return dir + "/" + filename;
  }


  inline Log& operator<<(Log& log, ostream& (*obj)(ostream &))
    {
      if(log.stream_to_logfile)
 	log.logfileout << obj;
      
      if(log.stream_to_cout)
  	cout << obj;
      
      return log;
    }
  
}

#endif

