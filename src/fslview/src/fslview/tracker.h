/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(TRACKER_H)
#define TRACKER_H

#include <qobject.h>

#include <boost/shared_ptr.hpp>
#include <string>
#include <memory>

//! @brief Disgnostic instrumentation class
//!
//! @author Dave Flitney
//!
//! Use TRACKER, MESSAGE and CHECKPOINT methods liberally throughout
//! your code for diagnostic output. Define the DEBUGGING symbol to
//! enable inclusion of the MACROS a-la assert().
class Tracker
{
public:
  typedef boost::shared_ptr< Tracker > Handle;

  static Handle create(const void *, const std::string &);

  void trace();
  void checkpoint();
  void message(const std::string &) const;

  virtual ~Tracker();

private:
  Tracker(const void *, const std::string &);
  
  const std::string message() const;
  unsigned int count() const;

  struct Implementation; 
  const std::auto_ptr<Implementation> m_impl;
};

#if defined(DEBUGGING)

//! @brief Use in static member functions
#define STATIC_TRACKER(s) Tracker::Handle t = Tracker::create(0x0, s)
//! @brief Add to methods to enable Tracker functionality
#define TRACKER(s) Tracker::Handle t = Tracker::create(this, s)
#define TRACE() t->trace()
//! @brief Output automatic checkpoint messages
#define CHECKPOINT() t->checkpoint()
//! @brief Output diagnostic message
#define MESSAGE(s) t->message(s)


#else

//! @brief Use in static member functions
#define STATIC_TRACKER(s)
//! @brief Add to methods to enable Tracker functionality
#define TRACKER(s)
#define TRACE()
//! @brief Output automatic checkpoint messages
#define CHECKPOINT()
//! @brief Output diagnostic message
#define MESSAGE(s)

#endif

#endif
