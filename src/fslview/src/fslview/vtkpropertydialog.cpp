/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Rama Aravind Vorray
		James Saunders
		David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "vtkpropertydialog.h"
#include "vtkwidget.h"

#include <qspinbox.h>
#include <qlineedit.h>
#include <qvalidator.h>
#include <qcolordialog.h>
#include <qpushbutton.h>

VTKPropertyDialog::VTKPropertyDialog(QWidget* parent, VTKProperties& props): 
  VTKPropertyDialogBase(parent), m_props(props)
{
//   m_lowerThreshold->setValue(m_props.inqLowerThreshold());
//   m_upperThreshold->setValue(m_props.inqUpperThreshold());
  m_mcThreshold->setValidator(new QDoubleValidator(this));
  m_mcThreshold->setText(tr("%1").arg(m_props.inqMcThreshold()));

  m_sd->setValidator(new QDoubleValidator(this));
  m_sd->setText( tr("%1").arg(m_props.inqStdDev()) );
  m_radius->setValidator(new QDoubleValidator(this));
  m_radius->setText( tr("%1").arg(m_props.inqRadius()) );

  m_iterations->setValue(m_props.inqIterations());
  m_relaxFactor->setValidator(new QDoubleValidator(this));
  m_relaxFactor->setText( tr("%1").arg(m_props.inqRelaxationFactor()) );

  m_interpMode->setCurrentItem(m_props.inqInterpMode());
  m_ambient->setValidator(new QDoubleValidator(this));
  m_diffuse->setValidator(new QDoubleValidator(this));
  m_opacity->setValidator(new QDoubleValidator(this));
  m_specular->setValidator(new QDoubleValidator(this));
  m_specularPower->setValidator(new QDoubleValidator(this));
  m_ambient->setText( tr("%1").arg(m_props.inqAmbient()) );
  m_diffuse->setText( tr("%1").arg(m_props.inqDiffuse()) );
  m_opacity->setText( tr("%1").arg(m_props.inqOpacity()) );
  m_specular->setText( tr("%1").arg(m_props.inqSpecular()) );
  m_specularPower->setText( tr("%1").arg(m_props.inqSpecularPower()) );

  m_featureAngle->setValidator(new QDoubleValidator(this));
  m_featureAngle->setText( tr("%1").arg(m_props.inqFeatureAngle()) );

  float cr, cg, cb;
  m_props.inqColor(cr, cg, cb);
  QColor color(int(cr * 255), int(cg * 255), int(cb * 255));
  m_colorSwatch->setPalette(QPalette(color));
  m_colorR->setValidator(new QDoubleValidator(this));
  m_colorR->setText( tr("%1").arg(cr) );
  m_colorG->setValidator(new QDoubleValidator(this));
  m_colorG->setText( tr("%1").arg(cg) );
  m_colorB->setValidator(new QDoubleValidator(this));
  m_colorB->setText( tr("%1").arg(cb) );
}

VTKProperties& VTKPropertyDialog::getProperties()
{
//   m_props.setLowerThreshold(m_lowerThreshold->value());
//   m_props.setUpperThreshold(m_upperThreshold->value());
  m_props.setMcThreshold(m_mcThreshold->text().toDouble());
  m_props.setInterpMode(m_interpMode->currentItem());
  m_props.setIterations(m_iterations->value());
  m_props.setAmbient(m_ambient->text().toDouble());
  m_props.setDiffuse(m_diffuse->text().toDouble());
  m_props.setOpacity(m_opacity->text().toDouble());
  m_props.setSpecular(m_specular->text().toDouble());
  m_props.setSpecularPower(m_specularPower->text().toDouble());
  m_props.setFeatureAngle(m_featureAngle->text().toDouble());
  m_props.setStdDev(m_sd->text().toDouble());
  m_props.setRadii(m_radius->text().toDouble());
  m_props.setColor(m_colorR->text().toDouble(), 
		   m_colorG->text().toDouble(), 
		   m_colorB->text().toDouble());
  return m_props;
}

void VTKPropertyDialog::selectColor()
{
  float cr, cg, cb;
  m_props.inqColor(cr, cg, cb);
  QColor orig = QColor(int(cr * 255), int(cg * 255), int(cb * 255));
  QColor color = QColorDialog::getColor(orig, this);
  if (color.isValid()) {
    m_colorSwatch->setPalette(QPalette(color));
    m_colorR->setText( tr("%1").arg(color.red()/255.0) );
    m_colorG->setText( tr("%1").arg(color.green()/255.0) );
    m_colorB->setText( tr("%1").arg(color.blue()/255.0) );
  }
}

#include "assistantclient.h"

void VTKPropertyDialog::help()
{
  AssistantClient::getInstance()->showPage("3D-dialog.html");
}
