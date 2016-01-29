/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(BRICON_H)
#define BRICON_H

#include <memory>
#include <boost/shared_ptr.hpp>

class BriConObserver;

/**
 * @author David Flitney <flitney@fmrib.ox.ac.uk>
 * @date   Thu Jan  2 15:18:22 2003
 * 
 * @brief  Class for managing brightness and contrast changes
 * 
 * A BriCon object tracks modifications to the brightness and contrast settings
 * by updating the range of voxel intensities that are displayable. The adjust
 * method can be called for a given image intensity to find the bricon adjusted
 * intensity.
 */
class BriCon {
public:
  typedef boost::shared_ptr< BriCon > Handle;

  static Handle create(float min, float max) { return Handle(new BriCon(min, max)); }

  float adjust(float v) const;

  void setRange(float min, float max);
  void modifyRange(float deltaBri, float deltaCon);
  void updateRange();
  void reset();
  void setMin(float min);
  void setMax(float max);
  float inqMin() const;
  float inqMax() const;
  float inqAdjustedMin() const;
  float inqAdjustedMax() const;
  Handle clone();

  void attach(BriConObserver* o);
  void detach(BriConObserver* o);
  void notify();

  virtual ~BriCon();

private:
  BriCon(float min, float max);
  BriCon(float origMin, float origMax, float min, float max);

  struct Implementation;  
  const std::auto_ptr< Implementation > m_impl;  
};

//! @brief interface for any class wishing to observe BriCon objects
//!
//! A class which wants to implement BriConObserver should subclass itself from
//! CursorObserver and implement the BriConObserver::update method.
class BriConObserver {
public:
  BriConObserver() {}

  virtual void update(const BriCon *) = 0;

  virtual ~BriConObserver() {}
};

#endif
