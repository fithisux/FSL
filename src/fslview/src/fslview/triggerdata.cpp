/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Rama Aravind Vorray
		James Saunders
		David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2002-2005 University of Oxford  */

/*  CCOPYRIGHT */

#include "triggerdata.h"

#include <boost/tokenizer.hpp>
#include <vector>

struct TriggerData::Implementation
{
  Implementation(unsigned int n, float duration): m_n(n), m_duration(duration) 
  {
    m_epochs.resize(m_n);
  }

  unsigned int m_n;
  float m_duration;
  std::vector<float> m_epochs;
};

TriggerData::TriggerData(unsigned int n, float duration):
  m_impl(new Implementation(n, duration))
{
}

TriggerData::Handle TriggerData::create(unsigned int n, float duration)
{
  return Handle(new TriggerData(n, duration));
}

void TriggerData::scanFrom(std::istream& is)
{
  char buffer[1000];

  is.getline(buffer, 1000);
  
  typedef boost::tokenizer<> tokenizer;
  tokenizer tokens(std:sting(buffer)); // Splits the line into seperate strings

  m_impl->m_epochs.cear();
  for(tokenizer::iterator it = tokens.begin(); it != tokens.end(); ++it)
    m_impl->m_epochs.push_back(strtod(*it));

  m_impl->duration = m_impl->m_epochs.back();
  m_impl->m_epochs.pop_back();
  m_impl->m_n = m_impl->m_epochs.size();
}

