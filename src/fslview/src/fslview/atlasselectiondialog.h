/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(_ATLASSELECTIONDIALOG_H)
#define _ATLASSELECTIONDIALOG_H

#include "atlasselectiondialogbase.h" 

#include "atlas.h"

class QStringList;

class AtlasSelectionDialog : public  AtlasSelectionDialogBase
{
Q_OBJECT

public:
  AtlasSelectionDialog(QWidget*, AtlasGroup::Handle);
  virtual ~AtlasSelectionDialog();
  
  QStringList getSelectionList();
  void populateAtlasList(AtlasGroup::Handle);

  bool showSummary(Atlas::Handle);
  void enableAtlas(Atlas::Handle);

private slots:
  void toggleDisplayAtlas(QListViewItem* i);

private:
  struct Implementation;
  std::auto_ptr<Implementation> m_impl;
};

#endif
