/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    David Flitney 

    FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include <iostream>

using namespace std;

#include <qstring.h>
#include <qmessagebox.h>
#include <qapplication.h>

#include "assistantclient.h"
#include "preferences.h"

AssistantClient* AssistantClient::m_instance=0;

AssistantClient* AssistantClient::getInstance()
{
  if(!m_instance) {
    Preferences p;
    QString path( p.inqAssistantPath() );
    m_instance = new AssistantClient(path);
  }
  return m_instance;
}

AssistantClient::AssistantClient(const QString& path):
  QAssistantClient( path, qApp->mainWidget())
{
  Preferences p;
  m_docPath = p.inqFSLDir() + QString("/doc/fslview");
  QStringList args;

  args << "-profile" << QString("%1/fslview.adp").arg(m_docPath);
  setArguments(args);
  connect( this, SIGNAL(error(const QString&)), SLOT(showError(const QString&)) );
}

void AssistantClient::showPage(const QString& page)
{
  QString fullPath(m_docPath + "/" + page);
  QAssistantClient::showPage(fullPath);
}

void AssistantClient::showError(const QString& msg)
{
  QMessageBox::warning(qApp->mainWidget(), "AssistantClient", msg);
}
