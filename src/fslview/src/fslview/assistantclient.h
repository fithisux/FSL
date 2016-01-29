/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    David Flitney 

    FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(ASSISTANTCLIENT_H)
#define ASSISTANTCLIENT_H

#include <qassistantclient.h>

class QString;

class AssistantClient: public QAssistantClient
{
Q_OBJECT
public:
  static AssistantClient* getInstance();

protected:
  AssistantClient(const QString&);

private slots:
  void showError(const QString&);
public slots:
  void showPage(const QString&);

private:
  QString m_docPath;
  static AssistantClient* m_instance;
};

#endif
