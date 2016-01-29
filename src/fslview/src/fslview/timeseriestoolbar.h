#if !defined(_timeseriestoolbar_h)
#define _timeseriestoolbar_h

#include "modelfit.h"
#include "timeseriestoolbarbase.h"

class QWidget;

class TimeSeriesToolbar: public TimeSeriesToolbarBase
{
public:
  TimeSeriesToolbar(QWidget *parent = 0);

  void populateFeatComboBox(const ModelFit::Handle&);
};

#endif
