
/*  FSLView - 2D/3D Interactive Image Viewer

    David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2006 University of Oxford  */

/*  CCOPYRIGHT */

#include "preferences.h"
#include "filemanager.h"
#include "talairachwidget.h"
#include "assistantclient.h"

#include "application.h"

#include "fslio/fslio.h"

#include <qapplication.h>
#include <qcursor.h>
#include <qtextbrowser.h>
#include <qstringlist.h>
#include <qcstring.h>
#include <qmessagebox.h>
#include <qcheckbox.h>

#include <iostream>
#include <string>

using namespace std;

AtlasGroup::Handle TalairachWidget::m_atlasGroup;

TalairachWidget::TalairachWidget(QWidget *parent, const Cursor::Handle& c, 
				 const OverlayList::Handle &ol):
  TalairachWidgetBase(parent), m_cursor(c), m_overlayList(ol),
  m_optionsDialog(this, m_options, ol, c),
  m_selectionDialog(this, m_atlasGroup)
{
  try {
    if( !m_atlasGroup )
      m_atlasGroup = AtlasGroup::create();

   if(m_atlasGroup) {
      m_selectionDialog.populateAtlasList(m_atlasGroup);
      m_selectionDialog.enableAtlas( m_atlasGroup->getAtlasByName("Harvard-Oxford Subcortical Structural Atlas") );
      m_selectionDialog.enableAtlas( m_atlasGroup->getAtlasByName("Harvard-Oxford Cortical Structural Atlas") );
      m_selectedAtlases = m_selectionDialog.getSelectionList();
      m_atlasGroup->selectCompatibleAtlases(ol->getMainImage());
    }
    update(m_cursor);

    m_cursor->attach(this);
  } catch (...) {
    QMessageBox::warning(this, "AtlasWidget", "Failed to initialise AtlasWidget");
  }
}

void TalairachWidget::help()
{
  AssistantClient::getInstance()->showPage( QString("./atlas.html") );
}

void TalairachWidget::readAtlas(const string& dirname, const string& fname)
{
  try {

    m_atlasGroup->readAtlas(dirname, fname);

  } catch (ios::failure &e) {
    QMessageBox::warning(this, "AtlasWidget", 
			 string("XML error while parsing atlas: ") + fname + "<br><br>" +
			 e.what());    
  } catch (Image::Exception &e) {
    QMessageBox::warning(this, "AtlasWidget", 
			 string("Exception while parsing atlas: ") + fname + "<br><br>" +
			 e.what());    
  }    
}

TalairachWidget::~TalairachWidget()
{
  m_cursor->detach(this);
}

void TalairachWidget::showSettingsDialog()
{
  if(m_selectionDialog.exec() == QDialog::Accepted) {
    m_selectedAtlases = m_selectionDialog.getSelectionList();
    for(AtlasGroup::ConstIterator it = m_atlasGroup->begin(); 
	it != m_atlasGroup->end(); ++it) {
      it->second->inqCurrentImage()->clearCache();
      Image::Handle im(it->second->inqCurrentSummaryImage());
      if( m_selectionDialog.showSummary(it->second) && 
	  im->getInfo()->isCompatible(m_overlayList->getMainImage()->getInfo()) ) {
	im->getInfo()->setPurpose(ImageIntent::Label);
	if(m_overlayList->getImageGroup()->addUniqueOverlay(im))
	  m_overlayList->setTransparency(0.5);
      }
      else
	m_overlayList->getImageGroup()->remOverlay(im);
      ApplicationWindow *w = dynamic_cast<ApplicationWindow*>(qApp->mainWidget());
      if(w)
	w->setFileMenuItemsState();
    }
  }

  QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );
  qApp->processEvents();
  update(m_cursor);
  QApplication::restoreOverrideCursor();
}

void TalairachWidget::showInspector()
{
  m_optionsDialog.show(m_atlasGroup);
}

void TalairachWidget::update(const Cursor::Handle& c)
{
  QStringList labels;
  try {
    for(QStringList::Iterator it = m_selectedAtlases.begin(); it != m_selectedAtlases.end(); ++it) {
      Atlas::Handle atlas( m_atlasGroup->getAtlasByName(*it) );
      Image::Handle image( m_overlayList->getActiveMetaImage()->getImage() );

      short radiogX(c->inqX());
      if(!image->getInfo()->isStoredRadiological()) radiogX = image->getInfo()->inqX()-1-radiogX;
      float x(0), y(0), z(0);
      image->getInfo()->voxToMMCoord(radiogX, c->inqY(), c->inqZ(),
				     x, y, z);

      labels.append( atlas->getDescription(x, y, z) );
    }

    m_text->setText( labels.join("<br>") );
  } catch (...) {
    QMessageBox::warning( this, "FSLView",
			  "AtlasWidget: An unexpected exception has occured!" );
  }
}

