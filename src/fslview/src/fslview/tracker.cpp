/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "tracker.h"

#include <stdio.h>

struct Tracker::Implementation {

  Implementation(const void *o, const std::string &msg): 
    m_object(o), m_msg(msg), m_count(0) { }

  ~Implementation() {}

  const void *m_object;
  static std::string m_indentation;
  const std::string m_msg;
  unsigned int m_count;
};

std::string Tracker::Implementation::m_indentation = ""; 

Tracker::Tracker(const void *o, const std::string &msg):m_impl(new Implementation(o, msg))
{
  m_impl->m_indentation = m_impl->m_indentation + " ";
  message("Begin");
}

Tracker::Handle Tracker::create(const void *object, const std::string &msg)
{ 
  return Tracker::Handle(new Tracker(object, msg)); 
}

Tracker::~Tracker()
{
  message("End");
  m_impl->m_indentation = m_impl->m_indentation.substr(0, m_impl->m_indentation.length() - 1);
}

// #include <execinfo.h>

void Tracker::trace()
{
//   void *array[10];
//   size_t size;
//   char **strings;
//   size_t i;
     
//   size = backtrace (array, 10);
//   strings = backtrace_symbols (array, size);
     
//   cout << "Obtained " << size << "stack frames." << endl;
  
//   for (i = 0; i < size; i++)
//     cout << strings[i] << endl;
     
//   free (strings);
}

void Tracker::checkpoint()
{
  m_impl->m_count++;
  warning("%s(%p): %d", message().c_str(), m_impl->m_object, count());
}

void Tracker::message(const std::string &msg) const
{
  warning("%s(%p): %s", message().c_str(), m_impl->m_object, msg.c_str());  
}

const std::string Tracker::message() const { return m_impl->m_indentation + m_impl->m_msg; }

unsigned int Tracker::count() const { return m_impl->m_count; }

//
// Some test code
//
#if defined(TESTING)
void anotherFunc()
{
  TRACKER("anotherFunc()");

  CHECKPOINT();
  MESSAGE("D'oh!");
}

void testTracker()
{
  TRACKER("testTracker()");

  CHECKPOINT();
  anotherFunc();
  CHECKPOINT();
}

int main()
{
  TRACKER("main");

  t->checkpoint();
  t->checkpoint();
  t->checkpoint();

  testTracker();
}
#endif
