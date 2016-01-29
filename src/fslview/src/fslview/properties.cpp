#include "properties.h"

struct Properties::Implementation
{
  Implementation(): m_create4dMask(false), m_askCreate4dMask(false) {}

  bool m_create4dMask;
  bool m_askCreate4dMask;
};

Properties::Handle Properties::create()
{
  return Properties::Handle(new Properties);
}

Properties::Properties(): m_impl(new Properties::Implementation())
{
}

Properties::~Properties()
{
}

bool Properties::inqAskCreate4dMask() const      { return m_impl->m_askCreate4dMask; }
bool Properties::inqCreate4dMask() const         { return m_impl->m_create4dMask; }
/*
 * param ask true if the system should ask the user what he/she wants to do
 */
void Properties::setAskCreate4dMask(bool ask)    { m_impl->m_askCreate4dMask = ask; }
/*
 * param create4d true if we should be allowing 4d masks
 */
void Properties::setCreate4dMask(bool create4d)  { m_impl->m_create4dMask = create4d; }
