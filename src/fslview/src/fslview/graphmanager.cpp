/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include <algorithm>
#include "graphmanager.h"
#include <qobject.h>

GraphManager::GraphManager()
{
  m_submittedCount = 0;
}

GraphManager::Handle GraphManager::create() 
{
   Handle g = Handle(new GraphManager());  
   g->setCountedThis(g);
   return g;
}

void GraphManager::attach(GraphManagerObserver *o)
{
  m_observers.push_back(o);
}

void GraphManager::detach(GraphManagerObserver *o)
{
  m_observers.remove(o);
}

struct Update
{
  Update(GraphManager::Handle g): m_gm(g) {}

  void operator()(GraphManagerObserver *v)
  {
    v->update(m_gm);
  }

  const GraphManager::Handle m_gm;
};

void GraphManager::notify() const
{
  std::for_each(m_observers.begin(), m_observers.end(), 
                Update(countedThis()));
}

void GraphManager::submitRange(double min, double max)
{
  double range = max - min;

  if (m_submittedCount == 0)
    {
      m_min = min;m_max = max; m_range = range;
    }
  else
    {
      if(min < m_min)m_min = min;
      if(max > m_max)m_max = max;
      if(range > m_range)m_range = range;
    }

  ++m_submittedCount;

  if(m_submittedCount >= m_observers.size())
    {
      m_submittedCount = 0;
      notify();
    }
}

