#include "propertiesdialogimpl.h"
#include "preferences.h"

#include "qlineedit.h"
#include "qtextedit.h"

PropertiesDialogImpl::PropertiesDialogImpl(QWidget *parent)
{
  Preferences prefs;
  
  m_fslDir->setText(prefs.inqFSLDir());
  m_mniImage->setText(prefs.inqMni152());
  m_atlasPath->setText(prefs.inqAtlasPath());
  m_assistantPath->setText(prefs.inqAssistantPath());
}

void PropertiesDialogImpl::commit()
{
  Preferences prefs;
  
  prefs.setFSLDir(m_fslDir->text());
  prefs.setMni152(m_mniImage->text());
  prefs.setAssistantPath(m_assistantPath->text());
  prefs.setAtlasPath(m_atlasPath->text());
}

PropertiesDialogImpl::~PropertiesDialogImpl()
{
}

void PropertiesDialogImpl::getProperties(QWidget *parent)
{
  PropertiesDialogImpl propertiesDialog(parent);

  if(propertiesDialog.exec() == QDialog::Accepted)
    {
      propertiesDialog.commit();
    }
}
