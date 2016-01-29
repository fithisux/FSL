/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Rama Aravind Vorray
		James Saunders
		David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */



#ifndef ATLASOPTIONSDIALOG_H
#define ATLASOPTIONSDIALOG_H

#include "atlasoptionsdialogbase.h"
#include "atlas.h"

#include "overlaylist.h"

class AtlasOptions;

class AtlasOptionsDialog : public AtlasOptionsDialogBase
{
Q_OBJECT

public:
  AtlasOptionsDialog(QWidget*, AtlasOptions&, OverlayList::Handle, Cursor::Handle);

  void setAtlas(Atlas::Handle);

  void show(AtlasGroup::Handle);

private slots:
  void selectAtlas(const QString&);
  void structureSelected(int);
  void superimpose(bool);
  void locate(bool);
  void addStructure();
  void accept();

private:
  Image::Handle getStructureImage(int);
  
  void show() {}

  AtlasOptions& m_options;
  AtlasGroup::Handle m_atlases;
  OverlayList::Handle m_overlayList;
  Cursor::Handle m_cursor;

  Atlas::Handle m_atlas;
  Image::Handle m_probImage;
};

#endif // ATLASOPTIONSDIALOG_H
