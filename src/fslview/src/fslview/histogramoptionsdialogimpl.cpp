#include "histogramoptionsdialogimpl.h"
#include "histogramwidget.h"

#include <qvalidator.h>
#include <qlineedit.h>
#include <qcheckbox.h>
#include <qspinbox.h>

HistogramOptionsDialogImpl::HistogramOptionsDialogImpl(QWidget *parent, HistogramOptions &options):
  HistogramOptionsDialog(parent), m_options(options)
{
  minIntensity->setValidator(new QDoubleValidator(this));
  intensityRange->setChecked(m_options.inqIntensityRange());
  minIntensity->setText(tr("%1").arg(options.inqMin()));
  maxIntensity->setText(tr("%1").arg(options.inqMax()));

  logScale->setChecked(m_options.inqLogScale());
  
  ignoreZeros->setChecked(m_options.inqIgnoreZeros());

  specifyBins->setChecked(m_options.inqSpecifyBins());
  numberOfBins->setValue(m_options.inqBins());
}

HistogramOptionsDialogImpl::~HistogramOptionsDialogImpl(void) {}

HistogramOptions& HistogramOptionsDialogImpl::getOptions(void) const
{
  m_options.setIntensityRange(intensityRange->isChecked());
  m_options.setMin(minIntensity->text().toDouble());
  m_options.setMax(maxIntensity->text().toDouble());

  m_options.setLogScale(logScale->isChecked());

  m_options.setIgnoreZeros(ignoreZeros->isChecked());

  m_options.setSpecifyBins(specifyBins->isChecked());
  m_options.setBins(numberOfBins->value());

  return m_options;
}
