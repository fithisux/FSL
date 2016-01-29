/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "curvedatalist.h"
#include <qobject.h>

class SetAllInActive
{
public:
  SetAllInActive(){}
  void operator()(CurveData::Handle cd) 
  {
    cd->setIsActive(false);
  }
};

class ActiveSearch 
{
public:
  ActiveSearch(){}
  void operator()(CurveData::Handle cd) 
  {
    if(cd->inqIsActive())m_curve = cd;
  }
  CurveData::Handle m_curve;
};

class BrowseSearch 
{
public:
  BrowseSearch(){}
  void operator()(CurveData::Handle cd) 
  {
    if(cd->inqBrowse())m_curve = cd;
  }
  CurveData::Handle m_curve;
};

class CurveDataSearch 
{
public:
  CurveDataSearch(long curve) : m_found(false),
                                m_curve(curve) {}
  void operator()(CurveData::Handle cd) 
  {
    if(cd->inqCurve() ==  m_curve){m_found = true;m_curveData = cd;}
  }
  bool m_found;
  CurveData::Handle m_curveData;
private:
  long m_curve;
};



class CoordinatesSearch 
{
public:
  CoordinatesSearch(CurveData::Handle newCd) : 
    m_found(false),m_newCd(newCd){}
  void operator()(CurveData::Handle cd) 
  {
    if ((cd->inqX() == m_newCd->inqX()) &&
        (cd->inqY() == m_newCd->inqY()) &&
        (cd->inqZ() == m_newCd->inqZ()) &&
        (cd->inqIndex() == m_newCd->inqIndex()))
      {
        m_found = true; m_curveData = cd;
      }
  }
  bool m_found;
  CurveData::Handle m_curveData;
private:
  CurveData::Handle m_newCd;
};

class MaxMinSearch 
{
public:
  MaxMinSearch() : 
    m_max(0),m_min(0),m_firstItem(true){}
  void operator()(CurveData::Handle cd) 
  {
    double min =  cd->inqTimeSeries()->inqMinVal();
    double max =  cd->inqTimeSeries()->inqMaxVal();

    if(m_firstItem)
      {
        m_min = min;
        m_max = max;
        m_firstItem = false;
      }
    else
      {
        if(min < m_min)m_min = min;
        if(max > m_max)m_max = max;
      }
  }

  double m_max;
  double m_min;
private:
  bool m_firstItem;
};


CurveData::Handle CurveData::create(const TimeSeries::Handle &ts,bool browse)
  {
    Handle dst(new CurveData(ts,browse,0));
    return dst;
  }
  
CurveData::Handle CurveData::create(const TimeSeries::Handle &ts,
                                    bool browse,int index)
{
    Handle dst(new CurveData(ts,browse,index));
    return dst;
}


CurveData::CurveData(const TimeSeries::Handle &timeSeries,bool browse,
                     int index):
  m_timeSeries(timeSeries),m_isActive(false),m_isBrowseCurve(browse),
  m_index(index)
  {
    
    
  }

CurveDataList::CurveDataList(){}

CurveDataList::Handle CurveDataList::create()
{
    Handle dst(new CurveDataList());
    return dst;
}

bool CurveDataList::push_back(CurveData::Handle cd)
{  
  CoordinatesSearch search = std::for_each(m_list.begin(),
                                           m_list.end(),
                                           CoordinatesSearch(cd));
  if(!search.m_found)
    {m_list.push_back(cd);}
  else
    {search.m_curveData->setBrowse(false);}

  return !search.m_found;
}

CurveDataList::It CurveDataList::begin()
{
  return m_list.begin();
}

CurveDataList::It CurveDataList::end()
{
  return m_list.end();
}



void CurveDataList::setAllInActive()
{
  std::for_each(m_list.begin(),
                m_list.end(),
                SetAllInActive());
}
  
CurveData::Handle CurveDataList::getCurveData(long curve)
{
  CurveDataSearch search = std::for_each(m_list.begin(),
                                         m_list.end(),
                                         CurveDataSearch(curve));
  return search.m_curveData;
}

CurveData::Handle CurveDataList::getActiveData()
{    
  ActiveSearch search = std::for_each(m_list.begin(),
                                      m_list.end(),
                                      ActiveSearch());

  return search.m_curve;
}


void CurveDataList::removeActive()
{
    ActiveSearch search = std::for_each(m_list.begin(),
                                         m_list.end(),
                                         ActiveSearch());

    if(isValidCurveData(search.m_curve))
      {
        m_list.erase(std::remove(m_list.begin(),m_list.end(),search.m_curve),
                   m_list.end());
      }
}
void CurveDataList::removeBrowse()
{
    BrowseSearch search = std::for_each(m_list.begin(),
                                         m_list.end(),
                                         BrowseSearch());

    if(isValidCurveData(search.m_curve))
      {
        m_list.erase(std::remove(m_list.begin(),m_list.end(),search.m_curve),
                   m_list.end());
      }
}
void CurveDataList::removeAll()
{
  m_list.clear();
}

CurveData::Handle CurveDataList::back()
{
  CurveData::Handle curve;

  if(!m_list.empty()){curve = m_list.back();}

  return curve;
}

double CurveDataList::inqMaxCurveValue() const
{
  MaxMinSearch search = std::for_each(m_list.begin(),
                                      m_list.end(),
                                      MaxMinSearch());
  
  return search.m_max;
}

double CurveDataList::inqMinCurveValue() const
{
  MaxMinSearch search = std::for_each(m_list.begin(),
                                      m_list.end(),
                                      MaxMinSearch());
  
  return search.m_min;
}
