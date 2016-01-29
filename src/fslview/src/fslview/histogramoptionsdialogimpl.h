
#ifndef HISTOGRAMOPTIONSDIALOGIMPL_H
#define HISTOGRAMOPTIONSDIALOGIMPL_H

#include "histogramoptionsdialog.h"

class HistogramOptions;

class HistogramOptionsDialogImpl : public HistogramOptionsDialog
{
public:
  HistogramOptionsDialogImpl( QWidget* parent, HistogramOptions& options );
  ~HistogramOptionsDialogImpl();

  HistogramOptions& getOptions(void) const;

  void setMin(unsigned int);
  void setMax(unsigned int);
  unsigned int inqMin(void) const;
  unsigned int inqMax(void) const;

private:
  HistogramOptions& m_options;
};

#endif // HISTOGRAMOPTIONSDIALOG_H
