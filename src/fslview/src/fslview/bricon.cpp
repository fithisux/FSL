/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "bricon.h"
#include <qobject.h>

#include "tracker.h"

struct BriCon::Implementation 
{
  Implementation(float min, float max):
    m_originalMin(min), m_originalMax(max),
    m_deltaBrightness(0),
    m_deltaContrast(0),
    m_minVal(min), m_maxVal(max) {}
    
  Implementation(float origMin, float origMax,
                 float min,     float max):
    m_originalMin(origMin), m_originalMax(origMax),
    m_deltaBrightness(0),
    m_deltaContrast(0),
    m_minVal(min), m_maxVal(max) {}
  
  void reset() { m_minVal = m_originalMin; m_maxVal = m_originalMax; }
  
  std::list<BriConObserver*> m_observers;
  
  float      m_originalMin;
  float      m_originalMax;
  float      m_deltaBrightness;
  float      m_deltaContrast;
  float      m_minVal;
  float      m_maxVal;
};

BriCon::BriCon(float min, float max):m_impl(new Implementation(min,max))
{
}

BriCon::BriCon(float origMin, float origMax, 
               float min,     float max):
               m_impl(new Implementation(origMin,origMax,min,max))
{
}

BriCon::~BriCon()
{
}

/** 
 * Attach an observer to this BriCon
 * 
 * @param o new BriConObserver which wants to recieve notifications
 */
void BriCon::attach(BriConObserver* o)
{
  m_impl->m_observers.remove(o);
  m_impl->m_observers.push_back(o);
}

/** 
 * Dettach an observer from this BriCon
 * 
 * @param o old BriConObserver which no longer wants to recieve notifications
 */
void BriCon::detach(BriConObserver* o)
{
  m_impl->m_observers.remove(o);
}

struct Update
{
  Update(BriCon* l): m_ol(l) {}

  void operator()(BriConObserver* v)
  {
    v->update(m_ol);
  }

  BriCon* m_ol;
};

/** 
 * Notify all current observers of a change to this BriCon
 */
void BriCon::notify()
{
  TRACKER("BriCon::notify()");
  MESSAGE(QString("Notifying %1 observers").arg(m_impl->m_observers.size()));

  std::for_each(m_impl->m_observers.begin(), m_impl->m_observers.end(), Update(this));
}

/** 
 * Set the min and max displayable intensities.
 * 
 * @param min new minimum displayable intensity.
 * @param max new maximum displayable intensity.
 */
void BriCon::setRange(float min, float max) { m_impl->m_minVal = min; m_impl->m_maxVal = max; notify(); }

/** 
 * Set the min displayable intensity.
 * 
 * @param min new minimum intensity.
 */
void BriCon::setMin(float min)              { m_impl->m_minVal = min; notify(); }

/** 
 * Set the max displayable intensity.
 * 
 * @param max new maximum intensity.
 */
void BriCon::setMax(float max)              { m_impl->m_maxVal = max; notify(); }

/** 
 * Get the current min intensity.
 * 
 * @return minimum displayable intensity.
 */
float BriCon::inqMin() const                { return m_impl->m_minVal; }

/** 
 * Get the current max intensity.
 * 
 * @return maximum displayable intensity.
 */
float BriCon::inqMax() const                { return m_impl->m_maxVal; }

/** 
 * Reset the range to original values
 * 
 */
void BriCon::reset() { m_impl->reset(); notify(); }

/** 
 * Used to perform temporary updates to the current range, e.g. while
 * a bricon control is being manipulated.
 * 
 * @param deltaBri fractional modifier for the brightness [-1,1]
 * @param deltaCon fractional modifier for the contrast [-1,1]
 */
void BriCon::modifyRange(float deltaBri, float deltaCon)
{
  m_impl->m_deltaBrightness = deltaBri;
  m_impl->m_deltaContrast = deltaCon;
  notify(); 
}

float BriCon::inqAdjustedMin() const
{
  float range = m_impl->m_maxVal - m_impl->m_minVal;
  
  return m_impl->m_minVal - (m_impl->m_deltaBrightness * range) + (m_impl->m_deltaContrast * range);
}

float BriCon::inqAdjustedMax() const
{
  float range = m_impl->m_maxVal - m_impl->m_minVal;
    return m_impl->m_maxVal - (m_impl->m_deltaBrightness * range) - (m_impl->m_deltaContrast * range);

}

/** 
 * Finalise any changes stored via modifyRange.
 * 
 */
void BriCon::updateRange()
{
  float range = m_impl->m_maxVal - m_impl->m_minVal;
  
  m_impl->m_minVal -= (m_impl->m_deltaBrightness * range) - (m_impl->m_deltaContrast * range);
  m_impl->m_maxVal -= (m_impl->m_deltaBrightness * range) + (m_impl->m_deltaContrast * range);
}

/** 
 * Calculate the adjusted intensity value. Do this before looking up your screen
 * value from say a look up table.
 * 
 * @param x the voxels intensity value
 * 
 * @return the adjusted value
 */
float BriCon::adjust(float x) const
{
  float min = m_impl->m_minVal;
  float max = m_impl->m_maxVal;
  float range = max - min;

  min -= (m_impl->m_deltaBrightness * range) - (m_impl->m_deltaContrast * range);
  max -= (m_impl->m_deltaBrightness * range) + (m_impl->m_deltaContrast * range);
  range = max-min;

  float m = ( range != 0.0) ? 1.0/range : 0.0; 
  float c = -min * m;
  float y = m * x + c;
  return y;
}

BriCon::Handle BriCon::clone()
{
  return  Handle(new BriCon(m_impl->m_originalMin, m_impl->m_originalMax,
                            m_impl->m_minVal,      m_impl->m_maxVal));
}
