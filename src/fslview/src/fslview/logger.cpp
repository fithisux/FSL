/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Rama Aravind Vorray
		James Saunders
		David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2002-2005 University of Oxford  */

/*  CCOPYRIGHT */

#include "config.h"
#if defined(HAS_SYSLOG)
# include <syslog.h>
#endif

#include <fstream>
#include "logger.h"

using namespace std;

struct Logger::Implementation
{
  Implementation(const string& filename):
    m_logfile(filename.c_str())
  {
  }
 
  void message(const string& facility, const string& message, const Logger::Level level)
  {
#if defined(HAS_SYSLOG)
    syslog(LOG_PID | level, "%s", message.c_str());
#else
    m_logfile << message << endl;
#endif
  }

  ofstream m_logfile;
};

Logger::Logger(const string& filename): m_impl(new Implementation(filename))
{
}

Logger::Handle Logger::create(const string& filename)
{
  return Logger::Handle(new Logger(filename));
}

void Logger::emergency(const string& facility, const string& message)
{
  m_impl->message(facility, message, Emergency);
}

void Logger::alert(const string& facility, const string& message)
{
  m_impl->message(facility, message, Alert);
}

void Logger::critical(const string& facility, const string& message)
{
  m_impl->message(facility, message, Critical);
}

void Logger::error(const string& facility, const string& message)
{
  m_impl->message(facility, message, Error);
}

void Logger::warning(const string& facility, const string& message)
{
  m_impl->message(facility, message, Warning);
}

void Logger::notice(const string& facility, const string& message)
{
  m_impl->message(facility, message, Notice);
}

void Logger::info(const string& facility, const string& message)
{
  m_impl->message(facility, message, Info);
}

void Logger::debug(const string& facility, const string& message)
{
  m_impl->message(facility, message, Debug);
}


