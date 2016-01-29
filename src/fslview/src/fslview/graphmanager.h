/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(GRAPHMANAGER_H)
#define GRAPHMANAGER_H

#include <list>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

class GraphManagerObserver;


class GraphManager
{
public:
  typedef boost::shared_ptr< GraphManager > Handle;
  typedef boost::weak_ptr< GraphManager > WeakHandle;

  static Handle create();
  void attach(GraphManagerObserver *o);
  void detach(GraphManagerObserver *o);
  virtual ~GraphManager() {}
  void notify() const;
  void submitRange(double min, double max);
  double inqMin(){return m_min;}
  double inqMax(){return m_max;}
  double inqRange(){return m_range;}

protected:
  Handle countedThis() const { return Handle(m_countedThis); }

private:
  GraphManager();
  void setCountedThis(const Handle c) { m_countedThis = WeakHandle(c); }
  WeakHandle m_countedThis;
  std::list< GraphManagerObserver * > m_observers;
  double m_min;
  double m_max;
  double m_range;
  unsigned int m_submittedCount;
};



class GraphManagerObserver
{
public:
  virtual ~GraphManagerObserver() {}

  virtual void update(const GraphManager::Handle& c) = 0;

  GraphManagerObserver() {}
};


#endif
