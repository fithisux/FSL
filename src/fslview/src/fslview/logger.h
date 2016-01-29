/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:	David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2002-2005 University of Oxford  */

/*  CCOPYRIGHT */

#include <string>
#include <boost/shared_ptr.hpp>

//////////////////////////////////////////////////
//! @author David Flitney <flitney@fmrib.ox.ac.uk>
//!
//! @date   Mon Dec 23 17:25:52 2002
//!
//! @brief Logger provides facilities to abstract the syslogd(8)
//! service found on most unix flavours.
//!
//! On systems without a syslog facility Logger will output to "facility".log
//! in the user's home directory. Otherwise it will send output via the syslog(3)
//! protocols.
class Logger
{
public:
  typedef enum {Emergency, Alert, Critical, Error, Warning, Notice, Info, Debug} Level;

  typedef boost::shared_ptr<Logger> Handle;

  static Handle create(const std::string& filename = std::string("fslview.log"));

  void emergency(const std::string& facility, const std::string& message);
  void alert(const std::string& facility, const std::string& message);
  void critical(const std::string& facility, const std::string& message);
  void error(const std::string& facility, const std::string& message);
  void warning(const std::string& facility, const std::string& message);
  void notice(const std::string& facility, const std::string& message);
  void info(const std::string& facility, const std::string& message);
  void debug(const std::string& facility, const std::string& message);
  
private:
  struct Implementation;
  const std::auto_ptr<Implementation> m_impl;

  Logger(const std::string& filename);
};
