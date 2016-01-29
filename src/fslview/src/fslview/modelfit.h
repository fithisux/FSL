/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2005 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(MODELFIT_H)
#define MODELFIT_H
#include "newmatap.h"
#include <qstringlist.h>
#include <boost/shared_ptr.hpp>
#include "miscmaths/miscmaths.h"
#include "storage/image.h"
#include <qstring.h>
#include <vector>
//#include "stdio.h"
//#include <stdlib.h>

#include <stdexcept>

using namespace NEWMAT;
using namespace MISCMATHS;

class ModelFitObserver;

//! @brief Evaluate PE and COPE images to display fitted model data.
//!
//! @author James Saunders
//! @author V Rama Aravind
//! @author Dave Flitney
//!
//! Read: design.con; design.fsf; parameter estimates and contrasts to
//! calculate model fits for a given voxel timeseries.
class ModelFit
{
public:
  typedef boost::shared_ptr< ModelFit > Handle;
  static Handle create(const QString & featDir);
  typedef class Exception;

  virtual ~ModelFit(); 

  QString& featDir() const;
  static bool isFeatDir(const QString& featDir);

  //  TimeSeries::Handle getModelTimeSeries(short x, short y, short z, float offset);
  TimeSeries::Handle getDataTimeSeries(short x, short y, short z);
  TimeSeries::Handle fullModel(short x, short y, short z, float mean); 
  TimeSeries::Handle CopeCurve(short x, short y, short z, float mean);
  TimeSeries::Handle peCurve(short x, short y, short z, float mean); 
  TimeSeries::Handle perCentChange(short x, short y, short z, float mean); 
  
  Image::Handle getFilteredFuncImage();

  virtual void attach(ModelFitObserver *);
  virtual void detach(ModelFitObserver *);
  virtual void notify();
  
  QString getConName(unsigned int i) const;
  
  unsigned int curFit(void) const;
  void curFit(unsigned int i);
  
  unsigned int numEVs(void) const;
  unsigned int numFits(void) const;

  unsigned int numContrasts(void) const;
  unsigned int numFtests(void) const;

  //  void copePe(bool);
  bool copePe(void) const;

private:
  ModelFit(const QString & featDir); 

  typedef struct Implementation;
  const std::auto_ptr<Implementation> m_impl;
};

class ModelFit::Exception: public std::runtime_error
{
public:
  Exception(const std::string& s): std::runtime_error(s) {}
};

class ModelFitObserver
{
public:
  virtual ~ModelFitObserver(){}
  virtual void update(ModelFit *m)=0;
  
protected:
  ModelFitObserver(){}
};

#endif
