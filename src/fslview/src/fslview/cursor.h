/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(CURSOR_H)
#define CURSOR_H

#include <list>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

#include <iostream>

class CursorObserver;

//! @brief An observable cursor class @see Design Patterns p293 Observer
//! @author Dave Flitney
class Cursor
{
public:
  typedef boost::shared_ptr< Cursor > Handle;
  typedef boost::weak_ptr< Cursor > WeakHandle;
  
  Handle clone();

  void setCursor(const Cursor::Handle c);
  void setCursor(short x, short y, short z, short v);
  void setCursor(short x, short y, short z);
  void setVolume(short v);
  void setCursorRepaint(short x, short y, short z);
  void repaint();
  bool  inqRepaint() const;
  short inqX() const;
  short inqY() const;
  short inqZ() const;
  short inqV() const;
  void setVMax(short);
  short inqVMax() const;

  static Handle create(CursorObserver *o, 
                       short xMax, short yMax, short zMax,short vMax);
  static Handle create(short xMax, short yMax, short zMax,short vMax);

  void attach(CursorObserver *o);
  void detach(CursorObserver *o);

  virtual ~Cursor() {}

protected:
  Handle countedThis() const { return Handle(m_countedThis); }

private:
  Cursor(CursorObserver *o,short xMax, short yMax, short zMax, short vMax);
  Cursor(short x,    short y,    short z,    short v,
         short xMax, short yMax, short zMax, short vMax);

  void setCountedThis(const Handle c) { m_countedThis = WeakHandle(c); }

  void notify() const;
  inline void setXYZ(short,short,short);
  inline void setV(short);

  short m_x;
  short m_y;
  short m_z;
  short m_v;
  short m_xMax;
  short m_yMax;
  short m_zMax;
  short m_vMax;

  bool  m_repaint;

  WeakHandle m_countedThis;

  std::list< CursorObserver * > m_observers;
};

inline std::ostream& operator << (std::ostream& os, const Cursor& c)
{
  os << c.inqX() << " " << c.inqY() << " " << c.inqZ() << " " << c.inqV();
  return os;
}

inline short Cursor::inqX() const { return m_x; }
inline short Cursor::inqY() const { return m_y; }
inline short Cursor::inqZ() const { return m_z; }
inline short Cursor::inqV() const { return m_v; }
inline short Cursor::inqVMax() const { return m_vMax; }
inline bool  Cursor::inqRepaint() const {return m_repaint;}

//! @brief interface for any class wishing to observe Cursor objects
//!
//! A class which wants to implement CursorObserver should subclass itself from
//! CursorObserver and implement the CursorObserver::update method.
class CursorObserver
{
public:
  virtual ~CursorObserver() {}

  //! @brief Essential API method for CursorObserver objects
  //! Re-implement in derived classes to recieve notifications
  //! of Cursor updates.
  virtual void update(const Cursor::Handle& c) = 0;

  CursorObserver() {}
};


#endif
