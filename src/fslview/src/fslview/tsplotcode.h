/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(TSPLOTCODE_H)
#define TSPLOTCODE_H
#include "newmatap.h"
#include <boost/shared_ptr.hpp>
#include "miscmaths/miscmaths.h"
#include "storage/timeseries.h"
#include <qstring.h>
#include <vector>
#include "stdio.h"

using namespace NEWMAT;
using namespace MISCMATHS;

class TsPlotCode
{
public:
  static void     preWhitenModel(const ColumnVector& ac,
                                 const Matrix& designMatrix,
                                 Matrix& preWhitenedMatrix);
  static void     establishPwFilter(const ColumnVector& ac, 
                                    ColumnVector& pwfilter, 
                                    int zeropad, int npts);
  static void     preWhitenData(const ColumnVector& data, 
                                ColumnVector& pwdata, 
                                ColumnVector& pwfilter, 
                                int zeropad, int npts);

  static TimeSeries::Handle     preWhitenTimeseries(const ColumnVector& ac, 
                                             TimeSeries::Handle& ts);  
  static ColumnVector convertTimeSeries(const TimeSeries::Handle&);
  static TimeSeries::Handle convertMatrix(const Matrix&, int col, 
                                   short x, short y, short z); 
  static TimeSeries::Handle convertColumnVector(const ColumnVector&, 
                                   short x, short y, short z);
};

#endif
