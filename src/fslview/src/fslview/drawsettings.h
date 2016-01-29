/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(DRAWSETTINGS_H)
#define DRAWSETTINGS_H

#if defined(WIN32)
#pragma warning (disable:4786)
#endif

#include <boost/shared_ptr.hpp>
#include <list>

class DrawSettingsObserver;

class Pen
{
public:
  Pen(): m_size(1), m_value(1) {}
  Pen(int s, int v): m_size(s), m_value(v) {}

  virtual ~Pen() {}

  void setSize(int s)   { m_size  = s; }
  void setValue(int v)  { m_value = v; }
  int  inqSize() const  { return m_size; }
  int  inqValue() const { return m_value; }

private:
  int m_size;
  int m_value;
};

//! @brief Stores current pen mode, value and size
//!
//! Used by mask drawing code to store the current pen settings. Stored along
//! with the Image::Handle in an OverlayList object one can track per overlay
//! mask drawing settings.
class DrawSettings
{
public:
  typedef boost::shared_ptr< DrawSettings > Handle;

  typedef enum {FreeHand, Erase, Fill} Mode;

  static Handle create();
  virtual ~DrawSettings();

  void setMode(Mode);
  void setPrevMode();

  void setPenSize(int);
  void setPenValue(int);
  int inqPenSize() const;
  int inqPenValue() const;

  void setLinkCursor(bool on) { m_linkCursor = on; notify(); }
  bool linkCursorOn() const   { return m_linkCursor; }

  void setColourIndex(int c) { m_colourIndex = c;notify(); }
  int inqColourIndex() const { return m_colourIndex;       }

  Mode inqMode() const       { return m_mode; }

  void attach(DrawSettingsObserver*);
  void detach(DrawSettingsObserver*);

private:

  DrawSettings();

  void switchPen();
  void notify() const;

  bool m_linkCursor;

  Mode m_mode;
  Mode m_prevMode;

  Pen m_pen;
  Pen m_eraser;
  Pen m_filler;
  Pen *m_currentPen;

  int m_colourIndex;

  std::list< DrawSettingsObserver* > m_observers;
};

//! @brief interface for any class wishing to observe DrawSettings objects
//!
//! A class which wants to implement DrawSettingsObserver should
//! subclass itself from DrawSettingsObserver and implement the
//! DrawSettingsObserver::update method.
class DrawSettingsObserver
{
public:
  virtual ~DrawSettingsObserver() {}

  virtual void update(const DrawSettings*) = 0;

  DrawSettingsObserver() {}
};

#endif
