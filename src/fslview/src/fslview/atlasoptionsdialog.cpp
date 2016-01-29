/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Rama Aravind Vorray
		James Saunders
		David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "atlasoptionsdialog.h"

#include <qlistbox.h>
#include <qcheckbox.h>
#include <qcombobox.h>

#include "atlas.h"

#include <iostream>
using namespace std;

AtlasOptionsDialog::AtlasOptionsDialog(QWidget *p, AtlasOptions &o,
				       OverlayList::Handle oh, Cursor::Handle ch):
  AtlasOptionsDialogBase(p),
//   AtlasOptionsDialogBase(p, "AtlasOptionsDialog", false, WStyle_Customize | WStyle_NormalBorder),
  m_options(o), m_overlayList(oh), m_cursor(ch)
{
}

Image::Handle AtlasOptionsDialog::getStructureImage(int index)
{
  Atlas::Handle atlas( m_atlases->getAtlasByName(m_atlasSelection->currentText()) );
  Image::Handle aimage( atlas->inqCurrentImage() );
  Image::Handle image;
  if(index >= 0) {
    image = aimage->clone3dStructure();
    *image->getVolume(0) = *aimage->getVolume(index);
    image->getInfo()->setReadOnly(true);
    image->getInfo()->setImageName(m_options.structureName());
  }
  return image;
}

void AtlasOptionsDialog::show(AtlasGroup::Handle atlases)
{
  m_atlases = atlases;
  m_atlasSelection->clear();
  for(AtlasGroup::ConstIterator it = atlases->begin(); it != atlases->end(); ++it)
    m_atlasSelection->insertItem(it->second->inqName());

  AtlasOptionsDialogBase::show();
  m_atlasSelection->setCurrentItem(1);
  selectAtlas(m_atlasSelection->currentText());
}

void AtlasOptionsDialog::selectAtlas(const QString& name)
{
  if( name != "None" )
    setAtlas( m_atlases->getAtlasByName(name) );
}

void AtlasOptionsDialog::setAtlas(Atlas::Handle atlas)
{
  m_atlas = atlas;
  m_structureList->clear();
  //  m_structureList->insertItem("None");
  for(Atlas::ConstLabelIterator it = m_atlas->begin(); it != m_atlas->end(); ++it)
    m_structureList->insertItem(it->second);
  if( atlas->inqType() != Atlas::Probabilistic ) {
    m_superimpose->setChecked(false);
    m_superimpose->setEnabled(false);   
  } else {
    m_superimpose->setEnabled(true);   
  }
}

void AtlasOptionsDialog::accept()
{
  m_superimpose->setChecked(false);
  AtlasOptionsDialogBase::accept();
}

void AtlasOptionsDialog::addStructure()
{
  if( Image::Handle simage = getStructureImage(m_options.structure()) ) {
    m_overlayList->getImageGroup()->addOverlay(simage);
    m_overlayList->getActiveMetaImage()->getDs()->inqBriCon()->setRange(5, 100);
    m_overlayList->setLookUpTable(LookUpTable::redYellow());
  }
}

void AtlasOptionsDialog::structureSelected(int i)
{
  m_options.structure(i);
  m_options.structureName(m_structureList->currentText());

  if(m_options.locate()) {
    Atlas::Handle atlas( m_atlases->getAtlasByName(m_atlasSelection->currentText()) );
    Image::Handle refimage( m_overlayList->getActiveMetaImage()->getImage() );
    
    m_cursor->setCursor( atlas->getCursor(refimage, m_options.structure()) );
  }

  if(m_options.superimpose()) {
    if(m_probImage)
      m_overlayList->getImageGroup()->remOverlay(m_probImage);
    if( (m_probImage = getStructureImage(m_options.structure())) &&
	m_overlayList->getMainImage()->getInfo()->isCompatible(m_probImage->getInfo()) ) {
      m_overlayList->getImageGroup()->addOverlay(m_probImage);
      m_overlayList->getActiveMetaImage()->getDs()->inqBriCon()->setRange(5, 100);
      m_overlayList->setLookUpTable(LookUpTable::redYellow());
    }
  } else {
    if(m_probImage) {
      m_overlayList->getImageGroup()->remOverlay(m_probImage);
      m_probImage.reset();
    }      
  }
}

void AtlasOptionsDialog::superimpose(bool y)
{
  m_options.superimpose(y);
  structureSelected(m_options.structure());
}

void AtlasOptionsDialog::locate(bool y)
{
  m_options.locate(y);
  structureSelected(m_options.structure());
}
