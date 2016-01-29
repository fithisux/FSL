/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Rama Aravind Vorray
		James Saunders
		David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2002-2005 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(CLUSTERBROWSER_H)
#define CLUSTERBROWSER_H

#include "clusterbrowserbase.h"

#include "cursor.h"
#include "modelfit.h"
#include "clusterdata.h"

#include "storage/image.h"

#include <vector>

class ClusterBrowser: public ClusterBrowserBase
{
  Q_OBJECT

public:
  ClusterBrowser(QWidget*, Image::Handle, 
		 Cursor::Handle, ModelFit::Handle);

  class Exception: public std::runtime_error
  {
  public:
    Exception(const std::string& s): std::runtime_error(s) {}
  };

signals:
  void windowClose(QCloseEvent*);

private slots:
  void selectStatistic(int);
  void clusterListSelectionChanged(QListViewItem *);
  void showTalairach(bool);
protected:
  virtual void closeEvent(QCloseEvent*);

private:

  bool m_showTalairach;
  unsigned short m_currentSelection;
  std::vector<ClusterTable> m_clusterTables;
  ImageInfo::Handle m_imageInfo;
  Cursor::Handle m_cursor;
  ModelFit::Handle m_model;
};

#endif

