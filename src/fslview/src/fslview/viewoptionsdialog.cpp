#include "viewoptionsdialog.h"
#include "viewoptions.h"

#include <qspinbox.h>
#include <qcheckbox.h>
#include <qradiobutton.h>

ViewOptionsDialog::ViewOptionsDialog(QWidget *p, ViewOptions& v):
  ViewOptionsDialogBase(p)
{
  m_linkLocalVolume->setChecked(v.inqVolumeIndexingWithinView());

  m_linkGlobalVolume->setChecked(v.inqUseSharedVolume());
  m_linkGlobalPosition->setChecked(v.inqUseSharedLocation());

  m_showLabels->setChecked(v.inqShowLabels());

  m_showSliceLabels->setChecked(v.inqShowSliceLabels());
  m_voxButton->setChecked(v.inqUnitsAreVoxels());
  m_mmButton->setChecked(!v.inqUnitsAreVoxels());

  m_showCursorGap->setChecked(v.inqShowCursorGap());
  m_gapSize->setValue(v.inqCursorGapSize());

  m_movieFrameRate->setValue(v.inqMovieFrameRate());
}

ViewOptionsDialog::~ViewOptionsDialog()
{
}

ViewOptions ViewOptionsDialog::getOptions() const
{
  ViewOptions v;

  v.setVolumeIndexingWithinView(m_linkLocalVolume->isChecked());

  v.setUseSharedVolume(m_linkGlobalVolume->isChecked());
  v.setUseSharedLocation(m_linkGlobalPosition->isChecked());

  v.setShowLabels(m_showLabels->isChecked());

  v.setShowSliceLabels(m_showSliceLabels->isChecked());
  v.setUnitsAreVoxels(m_voxButton->isChecked());

  v.setShowCursorGap(m_showCursorGap->isChecked());
  v.setCursorGapSize(m_gapSize->value());

  v.setMovieFrameRate(m_movieFrameRate->value());

  return v;
}
