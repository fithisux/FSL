/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

//#define DEBUGGING

#include <algorithm>
#include "cursor.h"

#include "tracker.h"

Cursor::Cursor(CursorObserver *o,
               short xMax, short yMax, short zMax, short vMax): 
               m_xMax(xMax), m_yMax(yMax), m_zMax(zMax),m_vMax(vMax),
               m_repaint(false)
{
  TRACKER("Cursor::Cursor(...) first c'tor");

  attach(o);
  setXYZ(0,0,0);
  setV(0);
  notify();
}

Cursor::Cursor(short x, short y, short z, short v,
               short xMax, short yMax, short zMax, short vMax):
  m_xMax(xMax),m_yMax(yMax),m_zMax(zMax),m_vMax(vMax),
  m_repaint(false)
{
  TRACKER("Cursor::Cursor(...) second c'tor");

  setXYZ(x,y,z);
  setV(v);
}

//! @brief Create a new cursor
//! @param o the initial observer for this cursor
//! @param maxX the extent of the cursor [0,maxX] in x
//! @param maxY the extent of the cursor [0,maxY] in y
//! @param maxZ the extent of the cursor [0,maxZ] in z
//! @param maxV the extent of the cursor [0,maxV] in volumes
Cursor::Handle Cursor::create(CursorObserver *o, 
                              short maxX, short maxY, short maxZ, short maxV) 
{
  STATIC_TRACKER("Cursor::create(...) first form");

  Handle c = Handle(new Cursor(o,maxX,maxY,maxZ,maxV));
  c->setCountedThis(c);
      
  return c;
}

//! @brief Create a new cursor
//! @param maxX the extent of the cursor [0,maxX] in x
//! @param maxY the extent of the cursor [0,maxY] in y
//! @param maxZ the extent of the cursor [0,maxZ] in z
//! @param maxV the extent of the cursor [0,maxV] in volumes
Cursor::Handle Cursor::create(short maxX, short maxY, short maxZ, short maxV) 
{
  STATIC_TRACKER("Cursor::create(...) second form");

  Handle c = Handle(new Cursor(0, 0, 0, 0, maxX,maxY,maxZ,maxV));
  c->setCountedThis(c);
      
  return c;
}

//! @brief Set the cursor position
//! @param x the x coordinate
//! @param y the y coordinate
//! @param z the z coordinate
void Cursor::setCursor(short x, short y, short z)
{
  TRACKER("Cursor::setCursor(short x, short y, short z)");
  MESSAGE(QString("x = %1, y = %2, z = %3").arg(x).arg(y).arg(z));

  setXYZ(x,y,z);
  notify();
}

//! @brief Set the cursor position
//! @param x the x coordinate
//! @param y the y coordinate
//! @param z the z coordinate
void Cursor::setCursorRepaint(short x, short y, short z)
{
  TRACKER("Cursor::setCursorRepaint(short x, short y, short z)");
  MESSAGE(QString("x = %1, y = %2, z = %3").arg(x).arg(y).arg(z));

  setXYZ(x,y,z);
  m_repaint = true;
  notify();
  m_repaint = false;
}

void Cursor::repaint()
{
  TRACKER("Cursor::repaint()");  

  m_repaint = true;
  notify();
  m_repaint = false;
}

Cursor::Handle Cursor::clone()
{
  TRACKER("Cursor::clone()");

  Handle c = Handle(new Cursor(m_x, m_y, m_z, m_v,
                           m_xMax, m_yMax, m_zMax, m_vMax));
  c->setCountedThis(c);

  return c;
}

//! @brief Set the cursor position
//! @param c a handle to an existing cursor whose locayion is to be copied
void Cursor::setCursor(const Cursor::Handle c)
{
  setCursor(c->inqX(), c->inqY(), c->inqZ(), c->inqV());
}

//! @brief Set the cursor position
//! @param x the x coordinate
//! @param y the y coordinate
//! @param z the z coordinate
//! @param v the volume number
void Cursor::setCursor(short x, short y, short z, short v)
{
  TRACKER("Cursor::setCursor(short x, short y, short z, short v)");
  MESSAGE(QString("x = %1, y = %2, z = %3, v = %4").arg(x).arg(y).arg(z).arg(v));

  setXYZ(x, y, z);
  setV(v);
  notify();
}

//! @brief Set the cursor position
//! @param v the volume number
void Cursor::setVolume(short v)
{
  TRACKER("Cursor::setVolume(short v)");
  MESSAGE(QString("v = %1").arg(v));

  setV(v);
  notify();
}

//! @brief Attach a view to this cursor
//! @param o handle of a view which requires notification of any 
//!          changes to this cursor
void Cursor::attach(CursorObserver *o)
{
  TRACKER("Cursor::attach(CursorObserver *o)");

  m_observers.push_back(o);
}

//! @brief Detach a view from this cursor
//! @param o handle of the view to be removed from the notification list
void Cursor::detach(CursorObserver *o)
{
  TRACKER("Cursor::dettach(CursorObserver *o)");

  m_observers.remove(o);
}

struct Update
{
  Update(Cursor::Handle c): m_cursor(c) {}

  void operator()(CursorObserver *v)
  {
    if(m_cursor) v->update(m_cursor);
  }

  const Cursor::Handle m_cursor;
};

void Cursor::notify() const
{
  TRACKER("Cursor::notify()");
  MESSAGE(QString("Notifying %1 observers").arg(m_observers.size()));

  std::for_each(m_observers.begin(), m_observers.end(), 
		Update(countedThis()));
}

void Cursor::setXYZ(short x,short y,short z)
{
  TRACKER("Cursor::setXYZ()");
  MESSAGE(QString("x = %1, y = %2, z = %3").arg(x).arg(y).arg(z));

  m_x = std::max(0, (int)x); m_x = std::min((int)m_xMax - 1, (int)m_x);
  m_y = std::max(0, (int)y); m_y = std::min((int)m_yMax - 1, (int)m_y);
  m_z = std::max(0, (int)z); m_z = std::min((int)m_zMax - 1, (int)m_z);
}

void Cursor::setV(short v)
{
  TRACKER("Cursor::setV(short v)");
  MESSAGE(QString("v = %1").arg(v));

  //  m_v = std::max(0, (int)v); m_v = std::min((int)m_vMax - 1, (int)m_v);
  m_v = v;
}

void Cursor::setVMax(short max)
{
  TRACKER("Cursor::setVMax(short max)");
  MESSAGE(QString("max -= %1").arg(max));

  if(max > m_vMax)
    m_vMax = max;

}

