/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */


// timeseries.h: interface for the TimeSeries class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(TIMESERIES_H)
#define TIMESERIES_H

#include <fslio/fslio.h>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <string>

//! @brief Templated time series storage class
//!
//! @author James Saunders
//! @author Dave Flitney
class TimeSeries
{
public:
  typedef boost::shared_ptr< TimeSeries > Handle;

  virtual float normalized (short nVol) = 0;
  virtual float value      (short nVol) = 0;
  virtual void  setValue   (int index, float val) = 0;
  virtual float inqMaxVal  () = 0;
  virtual float inqMinVal  () = 0;
  virtual int   inqVolCount() = 0;
  virtual short inqX() = 0;
  virtual short inqY() = 0;
  virtual short inqZ() = 0;
  virtual float mean()const = 0;

  TimeSeries() {}
  virtual ~TimeSeries() {}
};

template <class VoxelType>
class TimeSeriesStore: public TimeSeries
{
  // friend class TestVolumeStore;

public:
  typedef boost::shared_ptr< TimeSeriesStore<VoxelType> > Handle;

  virtual ~TimeSeriesStore();

  VoxelType& operator()    (short nVol);
  virtual float normalized (short nVol);
  virtual float value      (short nVol);
  virtual void  setValue   (int index, float val);
  virtual int   inqVolCount();
  virtual float inqMaxVal  ();
  virtual float inqMinVal  ();
  
  void setMin(VoxelType min) { m_min = min; }
  void setMax(VoxelType max) { m_max = max; }
  void setCoordinates(short x, short y, short z){m_x = x; m_y = y; m_z = z;}
  VoxelType inqMin() const   { return  m_min; }
  VoxelType inqMax() const   { return  m_max; }
  void setVolCount(int n)    { m_count = n; }
  virtual short inqX(){return m_x;}
  virtual short inqY(){return m_y;}
  virtual short inqZ(){return m_z;}
  virtual float mean() const;
  
  static Handle create(short x, short y, short z, int n);
  static Handle getTimeSeries(FSLIO* avw, short x, short y, short z);

private:
  TimeSeriesStore(VoxelType* buffer);
  void calculateMinVal();
  void calculateMaxVal();

  int m_count;
  VoxelType  m_min, m_max;
  float  m_minVal, m_maxVal;
  VoxelType* m_buffer;
  short m_x,m_y,m_z;
  bool m_doscaling;
  float m_slope, m_intercept;
};

typedef TimeSeriesStore<char> 			TimeSeriesB;
typedef TimeSeriesStore<unsigned char> 	TimeSeriesUB;
typedef TimeSeriesStore<short>			TimeSeriesS;
typedef TimeSeriesStore<unsigned short>	TimeSeriesUS;
typedef TimeSeriesStore<int>			TimeSeriesI;
typedef TimeSeriesStore<unsigned int>	TimeSeriesUI;
typedef TimeSeriesStore<long long>		TimeSeriesI64;
typedef TimeSeriesStore<unsigned long long>	TimeSeriesUI64;
typedef TimeSeriesStore<float>			TimeSeriesF;
typedef TimeSeriesStore<double>			TimeSeriesD;
typedef TimeSeriesStore<long double>	TimeSeriesF128;

#include "timeseries.inc"

#endif
