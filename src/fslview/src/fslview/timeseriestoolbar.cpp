#include "timeseriestoolbar.h"

#include <qcombobox.h>

TimeSeriesToolbar::TimeSeriesToolbar(QWidget *parent):
    TimeSeriesToolbarBase(parent)
{
}

void TimeSeriesToolbar::populateFeatComboBox(const ModelFit::Handle &m)
{
  m_contrastComboBox->insertItem("No model");
  m_contrastComboBox->insertItem("Full model only");
  for(unsigned int i=0; i<m->numFits(); i++)
    m_contrastComboBox->insertItem(m->getConName(i));
}
