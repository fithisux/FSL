/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(CURVEDATALIST_H)
#define CURVEDATALIST_H

#include "storage/timeseries.h"
#include <boost/shared_ptr.hpp>
#include <vector>

//! @brief Record details of curve data.
//! @author James Saunders
//!
//! Stores references to timeseries location, plot id and associated
//! properties.
class CurveData
{
public:  
  typedef boost::shared_ptr< CurveData > Handle;
  typedef enum {Null,FiltFunc,Full,Cope1,Cope2,Cope3,Cope4, PE} Feat;

  static Handle create(const TimeSeries::Handle &ts,bool browse);  
  static Handle create(const TimeSeries::Handle &ts,bool browse,
                       int index);
  TimeSeries::Handle inqTimeSeries(){return m_timeSeries;}
  long               inqCurve(){return m_curve;}
  void               setCurve(const long curve){m_curve = curve;}
  void               setBrowse(bool state){m_isBrowseCurve = state;}
  bool               inqBrowse(){return m_isBrowseCurve;}
  bool               inqIsActive(){return m_isActive;}
  void               setIsActive(bool state){m_isActive = state;}
  int                inqX(){return m_timeSeries->inqX();}
  int                inqY(){return m_timeSeries->inqY();}
  int                inqZ(){return m_timeSeries->inqZ();}
  float              inqYValue(short x){return m_timeSeries->value(x);}
  int                inqIndex(){return m_index;}

private: 
  
  CurveData(const TimeSeries::Handle &timeSeries,
            bool browse, int index);
  TimeSeries::Handle m_timeSeries;
  long m_curve;
  bool m_isActive;
  bool m_isBrowseCurve;
  int  m_index;
};
  
//! @brief Manage a list of CurveData objects
class CurveDataList
{
public:  
  typedef boost::shared_ptr< CurveDataList > Handle;
  typedef std::vector<CurveData::Handle>::iterator It;
  static Handle create();
  bool  push_back(CurveData::Handle);
  void setAllInActive();
  void removeActive();
  void removeBrowse();
  void removeAll();
  double inqMaxCurveValue()const;
  double inqMinCurveValue()const;

  CurveData::Handle back();
  CurveData::Handle getCurveData(long curve);
  CurveData::Handle getActiveData();
  It begin();
  It end();

private: 
  
  CurveDataList();
  std::vector<CurveData::Handle> m_list;
};

inline bool isValidCurveData(const CurveData::Handle cd)
{
  if(!cd.get()){return false;}
  else         {return true;}
}

#endif
