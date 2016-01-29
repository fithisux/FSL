/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Rama Aravind Vorray
		James Saunders
		David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "drawsettings.h"

DrawSettings::DrawSettings():
  m_linkCursor(false), m_mode(FreeHand), m_prevMode(FreeHand),
  m_pen(1, 1), m_eraser(1, 0), m_filler(1, 1), m_currentPen(&m_pen),
  m_colourIndex(1)
{
}

DrawSettings::~DrawSettings()
{
}

DrawSettings::Handle DrawSettings::create()
{
  return DrawSettings::Handle(new DrawSettings);
}

//! @brief Attach a viewer to this cursor
//! @param o handle of a viewer which requires notification of any 
//!          changes to this cursor
void DrawSettings::attach(DrawSettingsObserver *o)
{
  m_observers.push_back(o);
}

//! @brief Detach a viewer from this cursor
//! @param o handle of the viewer to be removed from the notification list
void DrawSettings::detach(DrawSettingsObserver *o)
{
  m_observers.remove(o);
}

struct Update 
{
  Update(const DrawSettings* s): m_settings(s) {}

  void operator()(DrawSettingsObserver *v)
  {
    v->update(m_settings);
  }

  const DrawSettings *m_settings;
};

void DrawSettings::notify() const
{
  std::for_each(m_observers.begin(), m_observers.end(), Update(this));
}

void DrawSettings::switchPen()
{
  switch(m_mode)
    {
    case DrawSettings::FreeHand: m_currentPen = &m_pen;    break;
    case DrawSettings::Erase:    m_currentPen = &m_eraser; break;
    case DrawSettings::Fill:     m_currentPen = &m_filler; break;
    }
}

void DrawSettings::setMode(Mode m)
{
  m_prevMode = m_mode; 
  m_mode = m; 
  switchPen();
  notify(); 
}

void DrawSettings::setPrevMode()
{ 
  m_mode = m_prevMode; 
  switchPen();
  notify(); 
}

void DrawSettings::setPenSize(int s)     { m_currentPen->setSize(s);  notify(); }
void DrawSettings::setPenValue(int s)    { m_currentPen->setValue(s); notify(); }
int DrawSettings::inqPenSize() const     { return m_currentPen->inqSize();  }
int DrawSettings::inqPenValue() const    { return m_currentPen->inqValue(); }
