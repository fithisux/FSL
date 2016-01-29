/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "atlasselectiondialog.h"

#include <qlistview.h>
#include <qvaluelist.h>
#include <qstringlist.h>
#include <qheader.h>

#include "eye.xpm"

//#define DEBUGGING
#include "tracker.h"

struct SelectionListItem: public QCheckListItem
{
  SelectionListItem(QListView* parent, const QString& text): 
    QCheckListItem(parent, text, CheckBox), m_preview(false) { refresh(); }

  void refresh();

  bool m_preview;
};

void SelectionListItem::refresh()
{
  if(m_preview)
    setPixmap(1, QPixmap(eye));
  else
    setPixmap(1, NULL);
}

struct AtlasSelectionDialog::Implementation
{
  Implementation(AtlasGroup::Handle ag):
    atlases(ag) { TRACKER("AtlasSelectionDialog::Implementation"); TRACE(); }
  ~Implementation() { TRACKER("AtlasSelectionDialog::~Implementation"); TRACE(); }

  void populateAtlasList(QListView* l)
  {
    for(AtlasGroup::ConstIterator it = atlases->begin(); it != atlases->end(); ++it)
      parentList.append( new SelectionListItem(l, it->second->inqName()) );
  }

  Image::Handle probImage;
  AtlasGroup::Handle atlases;
  QValueList<QListViewItem *> parentList;
};

AtlasSelectionDialog::~AtlasSelectionDialog() {}

AtlasSelectionDialog::AtlasSelectionDialog(QWidget* p, const AtlasGroup::Handle ag):
  AtlasSelectionDialogBase(p),
//   AtlasSelectionDialogBase(p, "AtlasSelectionDialog", true, 
// 			   WStyle_Customize|WStyle_DialogBorder), 
  m_impl(new Implementation(ag))
{
  m_atlasList->clear();
}

void AtlasSelectionDialog::toggleDisplayAtlas(QListViewItem* i)
{
  SelectionListItem *item = dynamic_cast<SelectionListItem*>(i);
  if(item) {
    item->m_preview ? item->m_preview = false : item->m_preview = true;
    item->refresh();
  }
}

bool AtlasSelectionDialog::showSummary(Atlas::Handle ah)
{
  SelectionListItem *item = 
    dynamic_cast<SelectionListItem*>( m_atlasList->findItem(ah->inqName(), 0) );
  if(item)
    return item->m_preview;
  else
    return false;
}

void AtlasSelectionDialog::enableAtlas(Atlas::Handle ah)
{
  if(ah)
    {
      SelectionListItem *i =
	dynamic_cast<SelectionListItem*>( m_atlasList->findItem(ah->inqName(), 0) );
      if(i) i->setOn(true);
    }
}

void AtlasSelectionDialog::populateAtlasList(AtlasGroup::Handle ag)
{
  m_impl->atlases=ag;
  m_impl->populateAtlasList(m_atlasList);
}

QStringList AtlasSelectionDialog::getSelectionList()
{
  QStringList atlasNames;
  QListViewItemIterator it(m_atlasList, QListViewItemIterator::Checked);
  while( it.current() ) {
    atlasNames.append(it.current()->text(0));
    ++it;
  }

  return atlasNames;
}

