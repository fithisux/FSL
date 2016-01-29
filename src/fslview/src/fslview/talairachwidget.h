/*  FSLView - 2D/3D Interactive Image Viewer

    David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2006 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(TALAIRACHWIDGET_H)
#define TALAIRACHWIDGET_H

#include "cursor.h"
#include "overlaylist.h"
#include "talairachwidgetbase.h"
#include "atlas.h"
#include "atlasoptionsdialog.h"
#include "atlasselectiondialog.h"

#include <vector>

class QStringList;

class TalairachWidget : public TalairachWidgetBase, CursorObserver
{
  Q_OBJECT

public:
  TalairachWidget(QWidget *parent, const Cursor::Handle& c, const OverlayList::Handle &ol);
  virtual ~TalairachWidget();

  virtual void update(const Cursor::Handle& c);

private slots:
  void showSettingsDialog();
  void showInspector();
  void help();

private:
  void readAtlas(const std::string&, const std::string&);
  Image::Handle getStructureImage(int);

  Cursor::Handle m_cursor;
  OverlayList::Handle m_overlayList;
  std::vector<MetaImage::Handle> m_imageVector;
  
  QStringList m_selectedAtlases;
  AtlasOptions m_options;
  AtlasOptionsDialog m_optionsDialog;
  AtlasSelectionDialog m_selectionDialog;

  static AtlasGroup::Handle m_atlasGroup;
};

#endif
