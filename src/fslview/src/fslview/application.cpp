/*  FSLView - 2D/3D Interactive Image Viewer
 
    V Rama Aravind, James Saunders, David Flitney, Mark Jenkinson,

    Christian Beckmann and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

/****************************************************************************
** $Id: application.cpp,v 1.146.2.3 2010/04/19 15:46:29 flitney Exp $
**
** Copyright (C) 2002 University of Oxford.  All rights reserved.
**
** FSLView
**
*****************************************************************************/

#if defined(WIN32)
#include <strstream>
//using namespace std;
#endif

#include "viewwidget.h"
#include "orthowidget.h"
#include "lightboxwidget.h"
#include "singlewidget.h"
#include "vtkwidget.h"
#include "timeserieswidget.h"
#include "clusterbrowser.h"
#include "histogramwidget.h"
#include "propertiesdialogimpl.h"
#include "createmaskdialog.h"
#include "application.h"
#include "assistantclient.h"

#include <qaction.h>
#include <qworkspace.h>
#include <qpixmap.h>
#include <qlabel.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qkeycode.h>
#include <qfiledialog.h>
#include <qstatusbar.h>
#include <qmessagebox.h>
#include <qapplication.h>
#include <qaccel.h>
#include <qwhatsthis.h>
#include <qcheckbox.h>
#include <qtoolbutton.h>
#include <qtooltip.h>

#include "version.h"
#include "tracker.h"
#include "preferences.h"
#include "storage/error.h"
#include "storage/image.h"
#include "fslio/fslio.h"

#include "modelfit.h"

#include <algorithm>
#include <functional>
#include <boost/bind.hpp>


bool ComparePaths(const Image::Handle im1, const Image::Handle im2)
{
  return (im1->getInfo()->inqFileName() == im2->getInfo()->inqFileName());
}

// // Concrete command class for FileOpen menu item
// class FileOpen : public Command
// {
// public:
//   FileOpen(ApplicationWindow *);
//   virtual void execute(void);
// private:
//   ApplicationWindow *m_applicationWindow;
// };

// // Concrete FileOpen command subroutines
// FileOpen::FileOpen(ApplicationWindow *app)
// {
// m_applicationWindow = app;
// }

// void FileOpen::execute(void)
// {
//   QString fn = QFileDialog::getOpenFileName( QDir::currentDirPath(), "Image files (*.hdr *.hdr.gz *.nii *.nii.gz)", 
// 					     m_applicationWindow );

//   if(ModelFit::isFeatDir(fn))
//   {
//     m_applicationWindow->loadFeat(fn);
//   }
//   else
//   {
//     if(!fn.isEmpty())
//       m_applicationWindow->loadFile(fn);
//   }
//   if(!fn.isNull())m_applicationWindow->setCurrentDir(fn);
// }

// const char * fileOpenText = "Click this button to open a <em>folder</em>.<br><hr>"
// "You can also select the <b>Open command</b> from the File menu.";

ApplicationWindow::ApplicationWindow(ApplicationOptions& options):
  ApplicationWindowBase( 0, "FslView", WDestructiveClose ), 
  m_properties(Properties::create()), 
  m_toolbarMenuId(0),
  m_options(options)
{  
  TRACKER("ApplicationWindow::ApplicationWindow");
  qApp->setMainWidget(this);

#if !defined(Q_OS_MACX)
# include "icons/fsllogo.xpm"
  setIcon( QPixmap(fsllogo_xpm) );
#endif

  int dw(QApplication::desktop()->width());
  int dh(QApplication::desktop()->height());

  setMinimumSize(200,200);
  Preferences p;
  setGeometry(p.inqGeometry(dw, dh));

  setupStatusBar(); 

  m_ws = new QWorkspace(this);
  m_ws->setBackgroundMode( PaletteMidlight );

  m_ws->setScrollBarsEnabled(true);

  setCentralWidget(m_ws);
  buildMenus(); 
  connectControls(); 
  show();
  
  emit message( "Ready", 2000 ); 
  
  // do not want fslview to crash if multiple files found!
  FslSetIgnoreMFQ(1);
  if (getenv("FSLOUTPUTTYPE")==NULL) {
    FslSetOverrideOutputType(FSL_TYPE_ANALYZE);
  }
  
  if(!m_options.empty())
    {
      OverlayOptionList::const_iterator it = m_options.begin();

      QString fn(it->fileInfo().absFilePath());
      if( loadFile(fn) ) {

	OverlayList::Handle ol(activeOverlayList());
	ImageDisplaySetting::Handle ds(ol->getActiveMetaImage()->getDs());

	if( ModelFit::isFeatDir(fn) )
	  loadFeat(fn);

	if(it->lutSpecified())
	  ol->setLookUpTable(m_imageGroup->getLut(it->lutname()));
	if(it->transparencySpecified())
	  ol->setTransparency(it->transparency());
	if(it->briconSpecified())
	  ds->inqBriCon()->setRange(it->min(), it->max());

	while(++it != m_options.end())
	  {
	    loadOverlay(it->fileInfo().absFilePath());

	    ImageDisplaySetting::Handle ds = ol->getActiveMetaImage()->getDs();

	    if(it->lutSpecified())
	      ol->setLookUpTable(m_imageGroup->getLut(it->lutname()));
	    if(it->transparencySpecified())
	      ol->setTransparency(it->transparency());
	    if(it->briconSpecified())
	      ds->inqBriCon()->setRange(it->min(), it->max());
	  }

	if(m_modelFit) {
	  Image::Handle image = m_modelFit->getFilteredFuncImage();

	  // Need to add this image iff not already in image group
	  ImageGroup::ImageList::iterator it = 
	    std::find_if( m_imageGroup->begin(), m_imageGroup->end(),
			  boost::bind(ComparePaths, _1, image) );
      
	  if( it ==  m_imageGroup->end() ) {
	    m_imageGroup->addOverlay(image);
	    activeOverlayList()->setVisibility(false);
	  }
	}

	std::list<ApplicationOptions::Mode> modes = m_options.inqModes();
	if(modes.size()) {
	  while(!modes.empty()) {
	    switch(modes.front()) {
	    case ApplicationOptions::Ortho :
	      viewOrthographic();
	      break;
	    case ApplicationOptions::Lightbox :
	      viewLightbox();
	      break;
	    case ApplicationOptions::Single :
	      viewSingle();
	      break;
	    case ApplicationOptions::ThreeD :
	      view3d();
	      break;
	    default:
	      break;
	    }
	    modes.pop_front();
	  }
	} else
	  viewOrthographic();
	
	if(m_modelFit) {
	  Image::Handle image = m_modelFit->getFilteredFuncImage();
	  // And display it as a timeseries with model viewing capabilities
	  TimeSeriesWidget* timeseries = new TimeSeriesWidget(m_ws, image, m_cursor, m_modelFit);
	  connect( timeseries, SIGNAL(windowClose(QCloseEvent*)), 
		   this, SLOT(childWindowClose(QCloseEvent*)));
	  timeseries->setCaption("Feat data");
	  timeseries->resize(250,250);
	  viewShow(timeseries);
	}
    
	/*
	 * set menu items to appropriate states after image has been 
	 * loaded, as they were set to some
	 * states, initially, when no images were loaded
	 */
	setFileMenuItemsState();
	setViewMenuItemsState();

	setCurrentDir(m_options.begin()->fileInfo().absFilePath());      
      } else
	QMessageBox::warning( this, "FSLView", "Failed to load base image. Command line processing aborted." );
    }
}


void ApplicationWindow::assistantError(const QString& msg)
{
  QMessageBox::warning(this, "FSLView - error while invoking help client", msg);
}


ApplicationWindow::~ApplicationWindow()
{
  TRACKER("ApplicationWindow::~ApplicationWindow");
}

//! @brief load an image file
//!
//! This method is used to load a given file into the application. 
//! It can be called from event handlers as well as directly from client code.
//! 
//! @param absFilePath the path to the file to be loaded
//!
//! @return true if operation succeeded
bool ApplicationWindow::loadFile(const QString & absFilePath)
{ 
  bool status(true);

  try {
    QString fn(absFilePath);
    QFileInfo fi(absFilePath);

    removeExtensions(fn);

    if ( !fn.isEmpty() && checkFilesExist(fn,false))
    {
      m_masterOverlayList.reset();
      m_modelFit.reset();
      m_imageGroup.reset();

      QApplication::setOverrideCursor(Qt::waitCursor);
      emit message(QString("Loading.... %1").arg(fn), 5000);  
      qApp->processEvents();
 
      m_imageGroup = ImageGroup::create( Image::load((const char *)fn) );
      ImageInfo::Handle info(m_imageGroup->getMainImage()->getInfo());
      m_cursor = Cursor::create(info->inqX(),info->inqY(),info->inqZ(), 
                                info->inqNumVolumes());

      OverlayList::Handle ol(activeOverlayList());
      Image::Handle im(m_imageGroup->getMainImage());

      LookUpTable::Handle lut;
      // Choose lut in image file
      if(info->inqLutName() != "")
	lut = m_imageGroup->getLut(info->inqLutName());
      // look for "stat/mask" in the filename and choose colour lut if true
      if ( info->isStatImage() || info->isMaskImage() )
	{
	  if(!lut) { lut = m_imageGroup->getNextLut(); }
	}
      if(!lut) lut = LookUpTable::greyScale();

      ol->getActiveMetaImage()->getDs()->setLookUpTable(lut);

      setCaption( fn );

      m_cursor->setCursor( info->inqX()/2, info->inqY()/2, info->inqZ()/2 ); 

      QApplication::restoreOverrideCursor();      
 
      if(checkForDuplicates(fn))
        {
          QMessageBox::warning(this,"FSLView",
          "Warning: Multiple versions of the image files exist!");
        }
      }  
    else 
      {
        if(!checkFilesExist(fn,false))
          {
	    QMessageBox::warning( this, "FSLView",
                                  "Missing header/image file" );
	    emit message( QString("Could not open %1").arg(fn), 2000 );
            if(!checkFilesExist(fn,true))
              {
                emit message( QString("Missing image file %1").arg(fn), 2000 );
              }
	    status = false;
          }
        else
          {
            emit message( "Loading aborted", 2000 );
          }
      }
  }
  catch(FileError& f)
    {
      catchFileError(f);
    }

  initFileMenuItems(status);
  initViewMenuItems(status);

  return status;
}

void ApplicationWindow::fileCloseWindow()
{
  QWidget* w = m_ws->activeWindow();
  if(w){w->close(true);}
}

/* Return whether the workspace (or main window) is empty or not. i.e., are any images opened
 */
bool ApplicationWindow::windowListEmpty()
{
  return m_ws->windowList().isEmpty();
}

int ApplicationWindow::windowListCount()
{
  return m_ws->windowList().count();
}

void ApplicationWindow::helpAbout()
{
  QMessageBox::about( this, "FSLView",
		      QString("<h2>FSLView</h2><hr><b>Version %1.%2</b><br>"
			      "Written by:"
			      "<ul>"
			      "<li>Dave Flitney</li>"
			      "<li>James Saunders</li>"
			      "<li>Mark Jenkinson</li>"
			      "<li>Steve Smith</li>"
			      "<li>V Rama Aravind</li>"
			      "</ul><br>"
			      "Copyright(c) 2004-2009 University of Oxford<br><hr>"
			      "<b>Help pages: http://www.fmrib.ox.ac.uk/fsl/fslview</b><hr>"
			      "Please report  bugs to:<br><b>fslview-bugs@fmrib.ox.ac.uk<b>")
		      .arg(Version).arg(Release));
}

void ApplicationWindow::helpAboutQt()
{
  QMessageBox::aboutQt( this, "fslview" );
}

void ApplicationWindow::windowMenuAboutToShow()
{ 
  TRACKER("ApplicationWindow::windowMenuAboutToShow()");

  // Remove any sub-window items
  while(Window->count() > 3)
    Window->removeItemAt(3);

  // Re-create the sub-window menu and enable the Cascade and Tile items
  // if appropriate.
#if (QT_VERSION < 0x030300)
  QWidgetList wl = m_ws->windowList();
#else
  QWidgetList wl = m_ws->windowList(QWorkspace::CreationOrder);
#endif
  windowCascadeAction->setEnabled(!wl.isEmpty());
  windowTileAction->setEnabled(!wl.isEmpty());
  if (!wl.isEmpty()) {
    for (unsigned int i = 0; i < wl.count(); ++i ) {
      int id = Window->insertItem(wl.at(i)->caption(),
				  this, SLOT( windowMenuActivated( int ) ) );
      Window->setItemParameter( id, i );
      Window->setItemChecked( id, m_ws->activeWindow() == wl.at(i) );
    }
  }
}

void ApplicationWindow::windowMenuActivated( int id )
{
  QWidget* w = m_ws->windowList().at( id );
  if ( w ) {
    w->showNormal();
    w->setFocus();
  }
}

void ApplicationWindow::helpOnlineHelp()
{
   AssistantClient::getInstance()->showPage( QString("index.html") );
}

void ApplicationWindow::help3DRendering()
{
   AssistantClient::getInstance()->showPage( QString("3D.html") );
}

void ApplicationWindow::menusUpdate()
{
  setFileMenuItemsState();
  setViewMenuItemsState();
}

void ApplicationWindow::fileMenuAboutToShow()
{
  setFileMenuItemsState();
}

void ApplicationWindow::viewMenuAboutToShow()
{  
  TRACKER("ApplicationWindow::viewMenuAboutToShow");
  setViewMenuItemsState();

  if(m_toolbarMenuId)
    Tools->removeItem(m_toolbarMenuId);

  if(QMainWindow *view = dynamic_cast<QMainWindow*>(m_ws->activeWindow()))
    {
      m_toolbarMenuId = Tools->insertItem("Toolbars", view->createDockWindowMenu());
    }
}

void ApplicationWindow::setupStatusBar()
{
  m_statusBar = statusBar();
  m_statusBar->addWidget(new QLabel(m_statusBar), 1, FALSE);
}

void ApplicationWindow::setFileMenuItemsState(void)
{
  bool remEnabled(false);
  bool state(windowListEmpty());  
  ViewWidget *view = dynamic_cast<ViewWidget*>(m_ws->activeWindow());

  if (!state && view)
  {
    OverlayList::Handle ol = view->getOverlayList();
    if(ol)
    {
       if(ol->getActiveMetaImage())
         remEnabled = !ol->getActiveMetaImage()->getInfo()->isMainImage();
    }
  }

  initFileMenuItems(state);
  fileRemoveAction->setEnabled(remEnabled);
}

void ApplicationWindow::setViewMenuItemsState(void)
{
  TRACKER("ApplicationWindow::setViewMenuItemsState(void)");
  bool state(!windowListEmpty());      
  bool multiVolume(false); 
  bool validImage(false);

  if(ViewWidget *view = dynamic_cast<ViewWidget*>(m_ws->activeWindow()))
    {
      OverlayList::Handle ol = view->getOverlayList();
      if(ol)
	{
	  Image::Handle image = ol->inqActiveImage();
	  if(isValidImage(image))
	    {
	      MESSAGE("Is valid image");
	      validImage= true;
	      multiVolume = (image->getInfo()->inqNumVolumes()>1);
	    }
	}
    }
  else
    state = false;

  /**
   * Initialize the Menuitems with requitred state, after user clicked on the menu button.
   *  original code above commented; Rama 3/11/04
   */
  initViewMenuItems(state);
  viewImageHistogramAction->setEnabled(state && validImage);  
  viewTimeseriesAction->setEnabled(state && multiVolume && validImage);
  viewClusterBrowserAction->setEnabled(state && m_modelFit);
}
  
void ApplicationWindow::initFileMenuItems(bool state)
{
  fileCreateMaskAction->setEnabled(!state);
  fileSaveAsAction->setEnabled(!state );
  fileOpen152Action->setEnabled(state);
  fileAdd152Action->setEnabled(!state);
  fileOpenAction->setEnabled(state);
  fileAddAction->setEnabled(!state);
  fileCloseAction->setEnabled(!state);
}

void ApplicationWindow::initViewMenuItems(bool state)
{
  viewOrthographicAction->setEnabled(state);
  viewLightboxAction->setEnabled(state);
  viewSingleAction->setEnabled(state);
  view3DViewerAction->setEnabled(state);
}

void ApplicationWindow::setMenuItems_NoImages(void)
{
  setFileMenuItemsState();
  setViewMenuItemsState();
}

void ApplicationWindow::fileOpen(void)
{
  QString fn = QFileDialog::getOpenFileName(QDir::currentDirPath(),
					    "Image files (*.hdr *.hdr.gz *.nii *.nii.gz)", 
					    this );
  if(!fn.isEmpty()) {
    loadFile(fn);
    if(ModelFit::isFeatDir(fn))
      loadFeat(fn);
    setCurrentDir(fn);

    viewOrthographic();

    if(m_modelFit) {
      Image::Handle image = m_modelFit->getFilteredFuncImage();
  
      // Need to add this image iff not already in image group
      ImageGroup::ImageList::iterator it = 
	std::find_if( m_imageGroup->begin(), m_imageGroup->end(),
		      boost::bind(ComparePaths, _1, image) );
      
      if( it ==  m_imageGroup->end() ) {
	m_imageGroup->addOverlay(image);
	activeOverlayList()->setVisibility(false);
      }
      // And display it as a timeseries with model viewing capabilities      
      TimeSeriesWidget* timeseries = new TimeSeriesWidget(m_ws, image, m_cursor, m_modelFit);
      connect( timeseries, SIGNAL(windowClose(QCloseEvent*)), 
	       this, SLOT(childWindowClose(QCloseEvent*)));
      timeseries->setCaption("Feat data");
      timeseries->resize(250,250);
      viewShow(timeseries);
    }

    /*
     * set menu items to appropriate states after image has been 
     * loaded, as they were set to some
     * states, initially, when no images were loaded
     */
    setFileMenuItemsState();
    setViewMenuItemsState();
  }
}

void ApplicationWindow::fileAdd()
{
  QString fn = QFileDialog::getOpenFileName( QDir::currentDirPath(), "Image files (*.hdr *.hdr.gz *.nii *.nii.gz)", this );

  if(!fn.isNull())setCurrentDir(fn);

  if(!fn.isEmpty()) {
    loadOverlay(fn);

    /*
     * set menu items to appropriate states after image has been 
     * loaded, as they were set to some
     * states, initially, when no images were loaded
     */
    setFileMenuItemsState();
    setViewMenuItemsState();
  }
}

void ApplicationWindow::fileOpen152(void)
{
  Preferences p;
  QString fn(QFileDialog::getOpenFileName(p.inqMni152().c_str(),
					  "Image files (*.hdr *.hdr.gz *.nii *.nii.gz)", 
					  this ));

  if(!fn.isEmpty()) {
    if( loadFile(fn) )
      viewOrthographic();

    /*
     * set menu items to appropriate states after image has been 
     * loaded, as they were set to some
     * states, initially, when no images were loaded
     */
    setFileMenuItemsState();
    setViewMenuItemsState();
  }
}

void ApplicationWindow::fileAdd152(void)
{
  Preferences p;
  QString fn(QFileDialog::getOpenFileName(p.inqMni152().c_str(),
					  "Image files (*.hdr *.hdr.gz *.nii *.nii.gz)", 
					  this ));

  if(!fn.isEmpty()) {
    loadOverlay(fn);

    /*
     * set menu items to appropriate states after image has been 
     * loaded, as they were set to some
     * states, initially, when no images were loaded
     */
    setFileMenuItemsState();
    setViewMenuItemsState();
  }
}

void ApplicationWindow::buildMenus()
{
  bool state(windowListEmpty());

  //
  // Initially disable the menu item-RemOverlay
  //
  //fileMenu->setItemEnabled(fileRemoveAction,false);

  connect( fileMenu, SIGNAL( aboutToShow() ),
           this, SLOT( fileMenuAboutToShow() ) );
  connect( Window, SIGNAL( aboutToShow() ),
	   this, SLOT( windowMenuAboutToShow() ) );

  connect( windowCascadeAction, SIGNAL( activated() ),
	   m_ws, SLOT( cascade() ));
  connect( windowTileAction, SIGNAL( activated() ),
	   m_ws, SLOT( tile() ));

  Window->setCheckable(true);
  
  connect( Tools, SIGNAL( aboutToShow() ),
	   this, SLOT( viewMenuAboutToShow() ) );

  state=!windowListEmpty();

  Tools->setCheckable(true);
  //
  // Initially disable TimeSeries and Histogram options
  //
  viewImageHistogramAction->setEnabled(false);
  viewTimeseriesAction->setEnabled(false);

  initViewMenuItems(state);
}

void ApplicationWindow::filePreferences()
{
  PropertiesDialogImpl::getProperties(this);
}

void ApplicationWindow::connectControls()
{  
  connect( this, SIGNAL(      message(const QString&, int)),
           this, SLOT( displayMessage(const QString&, int)));
  connect( this, SIGNAL(workSpaceEmpty(void)),
  	   this, SLOT(setMenuItems_NoImages(void)));
}

void ApplicationWindow::update(const Cursor::Handle& c)
{
  TRACKER("ApplicationWindow::update(const Cursor::Handle& c)");
  MESSAGE("Updating");
}

void ApplicationWindow::fileSaveAs()
{
  OverlayList::Handle ol = activeOverlayList();
  QString fn;
	
  if(!ol)
    {
      QMessageBox::warning(this, "FSLView",
			   "Active window must be an OrthoView, Lightbox or Single View");
      return;
    }

  Image::Handle image = ol->inqActiveImage();
    
  if(!isValidImage(image))
    {
      QMessageBox::warning(this,"FSLView",
			   "No image selected in Active window");
      return;
    }
   
  // Okay can go ahead and try to save it then
  QString initFileName;
  // use basename only
  initFileName = QString(image->getInfo()->inqImageName().c_str());

  fn = QFileDialog::getSaveFileName(initFileName, 
                                    "Image files (*.hdr *.hdr.gz *.nii *.nii.gz)", this,
                                    "save file dialog",
                                    "Select a filename for saving");
  	

  if(!fn.isNull())setCurrentDir(fn);    

  if(checkSpecificFilesExist(fn))   
    {
      if(QMessageBox::warning( this,
                               "FSLView",
                               "File already exists. Do you want to overwrite it?",
                               "Cancel","OK","",1,0) == 0) 
	{      
	  return;
	}
    }
  
  if ( !fn.isEmpty() ) 
    {
      if(image->save((const char *)fn))
	{
	  image->getInfo()->setTarnished(false);  
	  if(checkForDuplicates(fn))
	    {
	      QMessageBox::warning(this,"FSLView",
				   "Warning: Multiple versions of the image files exist!");
	    }
	}
      else
	{
	  QMessageBox::warning(this, "FSLView",
			       QString("<h2>Save failed!</h2>") +
			       "<p>Perhaps you're out of disk space or " +
			       "you don't have permission to write to " +
			       "this directory." +
			       "<p>Please see console for more details." +
			       "<p>Click Ok and try again.",
			       QMessageBox::Ok | QMessageBox::Default | QMessageBox::Escape,
			       QMessageBox::NoButton);
	}
    }
}

void ApplicationWindow::fileRemove()
{  
  
 OverlayList::Handle ol = activeOverlayList();
 if(!ol)
   {
     QMessageBox::warning(this, "FSLView",
                          "Active window must be either OrthoView, Lightbox or Single View");
   }
 else
  {
    Image::Handle mainImg = ol->getMainImage();
    Image::Handle image = ol->inqActiveImage();
    
    if(!isValidImage(image))
      {
        QMessageBox::warning(this,
                             "FSLView",
                             "No image selected in Active window");
      }
    else if(mainImg == image)
      {
        QMessageBox::warning(this,"FSLView","Cannot remove main image");
      }
    else
      {
        if(tarnishCheck(image)){ m_imageGroup->remOverlay(image);}
      }
  }
}

bool ApplicationWindow::loadOverlay(const QString & absFilePath)
{
  bool status(true);

  try{
    QString fn(absFilePath);
    QFileInfo fi(absFilePath);
    
    removeExtensions(fn);

    if ( !fn.isEmpty() && checkFilesExist(fn,false) ) 
      {
        QApplication::setOverrideCursor(Qt::waitCursor);
        emit message( QString("Loading.... %1").arg(fn), 5000);
	qApp->processEvents();

        Image::Handle overlay = Image::load((const char *)fn);

        if(m_imageGroup->getMainImage()->getInfo()->isCompatible(overlay->getInfo()))
          {
            m_imageGroup->addOverlay( overlay );
	    m_cursor->setVMax(overlay->getInfo()->inqNumVolumes());
          }
        else
          {
            QMessageBox::warning( this, "FSLView",
                                  QString("Unable to load incompatible overlay!<br><br>") +
				  "All overlays <b>must</b> have same dimensions as the base image!");
            emit message( "Loading aborted", 2000 );
	    status = false;
          }
        
	QApplication::restoreOverrideCursor();

        if(checkForDuplicates(fn))
        {
          QMessageBox::warning(this,"FSLView",
          "Warning: Multiple versions of the image files exist!");
        }
      } 
    else 
      {        
        if(!checkFilesExist(fn,false))
          {
            QMessageBox::warning( this, "FSLView",
                                  "Missing header/image file" );

	    emit message( QString("Could not open %1").arg(fn), 2000 );            
	    if(!checkFilesExist(fn,true))
	      {
		emit message( QString("Missing image file %1").arg(fn), 2000 );
	      }
	    status = false;
          }
        else
          {
            emit message( "Loading aborted", 2000 );
          }
	status = false;
      }
  }
  catch(FileError& f)
    {
      catchFileError(f);
    }

  setFileMenuItemsState();

  return status;
}

void ApplicationWindow::addLookUpTable()
{
  QString fn = QFileDialog::getOpenFileName( QDir::currentDirPath(), "LUTs (*.lut *.rgb)", this );
  if ( !fn.isEmpty() ) {
    QApplication::setOverrideCursor(Qt::waitCursor);
    emit message( QString("Loading lookup table.... %1").arg(fn), 2000 );
    LookUpTable::Handle lookUpTable = LookUpTable::load((const char *)fn);
    m_imageGroup->addLookUpTable( lookUpTable );

    QApplication::restoreOverrideCursor();

  }  else {
    emit message( "Loading aborted", 2000 );
  }
}


void ApplicationWindow::fileCreateMask()
{
 TRACKER("ApplicationWindow::createMask");      
 OverlayList::Handle ol = activeOverlayList();

 if(!ol)
   {
     QMessageBox::warning(this,
                          "FSLView",
                          "Active window must be either OrthoView, Lightbox or Single View");
   }
 else
  {
    Image::Handle mainImg = ol->getMainImage();
    Image::Handle image = ol->inqActiveImage(), ci;
    
    if(!isValidImage(image))
      {
        QMessageBox::warning(this,
                             "FSLView",
                             "No image selected in Active window");
      }
    else
    {
      /*********************************
      Check if the main image is a 4D or 3D image. If the image is a 4D image,
      i.e., if (info->inqNumVolumes()>1) is TRUE, then ask for creating a 3D
      mask or 4D mask. If the main image is 3D, then dont ask for creating
      4D mask (point less), i.e., ELSE block is executed
      *********************************/
      ImageInfo::Handle info(m_imageGroup->getMainImage()->getInfo());
      if(info->inqNumVolumes()>1)
      {
        std::auto_ptr<CreateMaskDialog> cd(new CreateMaskDialog(this));
        
        cd->setFixedSize(300, 100);

        cd->m_create4dMask->setChecked(m_properties->inqCreate4dMask());
        cd->m_dontAsk->setChecked(m_properties->inqAskCreate4dMask());

        if(!m_properties->inqAskCreate4dMask())
          cd->exec();
        if(cd->m_create4dMask->isChecked())
          ci = image->cloneStructure();
        else
          ci = image->clone3dStructure();
        m_properties->setCreate4dMask(cd->m_create4dMask->isChecked());
        m_properties->setAskCreate4dMask(cd->m_dontAsk->isChecked());
      }
      else
      {
        ci = image->clone3dStructure();
      }

      m_imageGroup->addOverlay(ci);

      MetaImage::Handle cmi = ol->getMetaImage(ci);
      cmi->getDs()->setLookUpTable(LookUpTable::redYellow());

      emit message( QString("Created new mask: %1")
                    .arg(ci->getInfo()->inqImageName().c_str()), 2000 );
      /**
      needed for the first time execution. Even though redundant after wards, no side affect.
      **/
      fileRemoveAction->setEnabled(true);
    }
  }

}

void ApplicationWindow::viewOrthographic()
{  
  ViewWidget* view = new OrthoWidget(m_ws, m_imageGroup,
                                     copyActiveOverlayList(), m_cursor);

  connect( view, SIGNAL(message(const QString&, int)), 
           this, SIGNAL(message(const QString&, int)) );
  connect( view, SIGNAL(addLookUpTable()),
           this, SLOT(addLookUpTable()));
  connect( view, SIGNAL(windowClose(QCloseEvent*)),
           this, SLOT(childWindowClose(QCloseEvent*)));

  connect( view, SIGNAL(overlayEvent()),
           this, SLOT(menusUpdate()) );

  view->setCaption("Ortho view");

  view->resize(580,500);

  viewShow(view);
  m_cursor->repaint();
}

void ApplicationWindow::viewLightbox()
{
  ViewWidget* view = new LightboxWidget(m_ws, m_imageGroup,
                                        copyActiveOverlayList(), m_cursor); 

  connect( view, SIGNAL(message(const QString&, int)), 
           this, SIGNAL(message(const QString&, int)));   
  connect( view, SIGNAL(addLookUpTable()),
           this, SLOT(addLookUpTable()));  
  connect( view, SIGNAL(windowClose(QCloseEvent*)),
           this, SLOT(childWindowClose(QCloseEvent*)));

  connect( view, SIGNAL(overlayEvent()),
           this, SLOT(menusUpdate()) );

  view->setCaption("Lightbox view");  

  view->resize(580,500);

  viewShow(view);
  m_cursor->repaint();
}

void ApplicationWindow::viewSingle()
{
  ViewWidget* view = new SingleWidget(m_ws, m_imageGroup,
                                     copyActiveOverlayList(), m_cursor); 

  connect( view, SIGNAL(message(const QString&, int)), 
           this, SIGNAL(message(const QString&, int)) );   
  connect( view, SIGNAL(addLookUpTable()),
           this, SLOT(addLookUpTable()));  
  connect( view, SIGNAL(windowClose(QCloseEvent*)),
           this, SLOT(childWindowClose(QCloseEvent*)));

  connect( view, SIGNAL(overlayEvent()),
           this, SLOT(menusUpdate()) );

  view->setCaption("Single view");  

  view->resize(210,345);

  viewShow(view);  
}

void ApplicationWindow::view3d()
{
  ViewWidget* view = new VTKWidget(m_ws, m_imageGroup,
				   copyActiveOverlayList(),
				   m_cursor); 

  connect( view, SIGNAL(message(const QString&, int)), 
           this, SIGNAL(message(const QString&, int)) );   
  connect( view, SIGNAL(addLookUpTable()),
           this, SLOT(addLookUpTable()));  
  connect( view, SIGNAL(windowClose(QCloseEvent*)),
           this, SLOT(childWindowClose(QCloseEvent*)));

  connect( view, SIGNAL(overlayEvent()),
           this, SLOT(menusUpdate()) );

  view->setCaption("3D view");  

  view->resize(210,345);

  viewShow(view);  
  m_cursor->repaint();
}

void ApplicationWindow::viewClusterBrowser()
{
  try {
    if(m_modelFit) {
      ClusterBrowser *cb = 
	new ClusterBrowser(m_ws, m_imageGroup->getMainImage(), 
			   m_cursor, m_modelFit);
      connect( cb,   SIGNAL(windowClose(QCloseEvent*)),
	       this, SLOT(childWindowClose(QCloseEvent*)));
      cb->setCaption("Cluster Browser");
      cb->showNormal();
    } else {
      QMessageBox::warning( this, "Cluster Browser", "Unable to open browser: no valid model!"); 
    }
  } catch(const ClusterBrowser::Exception& e) {
    QMessageBox::warning( this, "Cluster Browser", e.what());
  }
}

void ApplicationWindow::viewTimeseries()
{
  if(m_imageGroup.use_count())
  {  
    OverlayList::Handle ol = activeOverlayList();
    
    if(ol)
    {
      Image::Handle image = ol->inqActiveImage();
    
      if(isValidImage(image))
        {
          if(image->getAvw() == NULL)
            {
              QMessageBox::warning( this,
                                    "FSLView",
                                    "New masks must be edited, saved and reloaded before viewing their timeseries. This is a known error.");
            }
          else
            {
	      TimeSeriesWidget* timeseries = new TimeSeriesWidget(m_ws, image, m_cursor, m_modelFit);
              connect( timeseries, SIGNAL(windowClose(QCloseEvent*)),
                       this, SLOT(childWindowClose(QCloseEvent*)));
              timeseries->setCaption("Timeseries");    
              timeseries->resize(250,250);
              viewShow(timeseries);
            }
        }
    }
  }
}


void ApplicationWindow::viewImageHistogram()
{
  if(m_imageGroup.use_count())
  {  
    if(ViewWidget *view = dynamic_cast<ViewWidget*>(m_ws->activeWindow())) {
      if(OverlayList::Handle ol = view->getOverlayList()) {
	if(Image::Handle image = ol->inqActiveImage()) {
	  unsigned int v = m_cursor->inqV();
	  HistogramWidget* histogram = 
	    new HistogramWidget(m_ws, image->getVolume(v),
				image->getInfo()->inqImageName(),
				v, image->getInfo()->isInteger());   
          connect( histogram, SIGNAL(windowClose(QCloseEvent*)),
                   this, SLOT(childWindowClose(QCloseEvent*))); 
          viewShow(histogram);
        }
      }
    }
  }  
}

void ApplicationWindow::catchFileError(FileError f)
{
  emit message( QString("Error loading %1, %2")
                        .arg(f.inqFileName().c_str())
                        .arg(f.inqMessage().c_str()),
                         3000 ); 
}

bool ApplicationWindow::checkAbsFilePath(const QString & absFilePath,
                                         QString ext)
{
  return QFile::exists(absFilePath + "." + ext);
}

OverlayList::Handle ApplicationWindow::activeOverlayList()
{
  OverlayList::Handle ol;

  ol = m_masterOverlayList;

  if (!windowListEmpty()) {
    ViewWidget* view = dynamic_cast<ViewWidget*>(m_ws->activeWindow());
    if(view)
      ol = view->getOverlayList();
  }

  if (!ol && m_imageGroup ) {
    ol = OverlayList::create(m_imageGroup);
    m_masterOverlayList = ol;
  }

  return ol;
}

OverlayList::Handle ApplicationWindow::copyActiveOverlayList()
{
  OverlayList::Handle null;
  OverlayList::Handle ol(activeOverlayList());
  
  if(ol) {return ol->clone();}
  else   {return null;}
}
 
void ApplicationWindow::displayMessage(const QString & msg, int time)
{
  if(time == -1)  statusBar()->message(msg);
  else            statusBar()->message( msg, time );
}

void ApplicationWindow::viewShow(QWidget* v)
{
  if(windowListEmpty())
    { 
      v->showMaximized();
    }
  else if(m_ws->windowList().count() == 1) 
    {
      m_ws->windowList().at(0)->showNormal();
      v->showNormal();
      m_ws->tile();
    } 
  else
    {
      v->show();
    }
}

bool ApplicationWindow::tarnishCheck(Image::Handle& image)
{
  bool result(true);

  if(image->getInfo()->inqTarnished())
    {
        QString imageMessage = QString("Continuing this action will lose changes to %1?")
          .arg(image->getInfo()->inqImageName().c_str());
        
        switch(QMessageBox::warning( this,
                               "FSLView",
                               imageMessage,"Cancel","OK","",0,0))
          {
          case 0:result = false;break;
          case 1:result = true ;break;
          }
    }
  
  return result;
}

bool ApplicationWindow::tarnishCheck()
{

 bool result(true);

  if(m_imageGroup->inqTarnished())
    {
         switch(QMessageBox::warning( this,
				      "FSLView", 
				      QString("Continuing this action will lose unsaved data."),
				      "Cancel","OK","",0,0))
          {
          case 0:result = false;break;
          case 1:result = true;break;
          }
    }
  
  return result;
}

void ApplicationWindow::childWindowClose(QCloseEvent* e)
{
  if(windowListCount() > 1)
    {
      e->accept();
    }
  else
    {
      if(tarnishCheck())e->accept();
    }
  // schedule a menu item check till the child window is actually closed by the application.
  QTimer::singleShot(1, this, SLOT(setMenuItems_NoImages()));
}

void ApplicationWindow::closeEvent(QCloseEvent* e)
{
  Preferences p;
  p.setGeometry(geometry());

  if(windowListEmpty())
    {
      e->accept();
    }
  else
    {
      if(tarnishCheck())
	e->accept();
    }
}

void ApplicationWindow::setCurrentDir(QString path)
{
  QFileInfo fi(path);
  if(!QDir::setCurrent(fi.dirPath(true)))
    warning("Unable to set Current Directory path");
}

void ApplicationWindow::removeExtensions(QString & fileName)
{
  //Strips off the gz and then the hdr and img or nii extensions
  // MJ NOTE: not using FSLIO code as I don't know to integrate it with Q* calls

  QFileInfo fiA(fileName);
  if(fiA.extension(false) == "gz")  
    fileName.remove( fileName.findRev("."), 3 );

  QFileInfo fiB(fileName);
  if((fiB.extension(false) == "img") || (fiB.extension(false) == "hdr") || 
     (fiB.extension(false) == "nii"))
    fileName.remove( fileName.findRev("."), 4 );
}
  
bool ApplicationWindow::checkFilesExist(const QString & fn, bool justImg)// why is bool necessary here?
{
  //Checks that files exist
  return (FslFileExists(fn.ascii())>0);
}
  

bool ApplicationWindow::checkSpecificFilesExist(const QString & fn)
{
  // MJ NOTE: why is this different from checkFilesExist ?
  return (FslFileExists(fn.ascii())>0);
}

bool ApplicationWindow::checkForDuplicates(const QString & fn)
{
  // for now just disable this and see if things can work
  return (FslCheckForMultipleFileNames(fn.ascii())>0);
}

void ApplicationWindow::loadFeat(const QString &fn)
{
  QFileInfo  fi(fn);
  
  try {
    m_modelFit = ModelFit::create(fi.dirPath(true));
  
    emit message("*File Opened is a FEAT analysis Directory*", 2000);
  } catch (std::ios::failure& e) {
    QMessageBox::warning(this, "Attempting to load FEAT model", e.what());
  } catch (Image::Exception& e) {
    QMessageBox::warning(this, "Attempting to load FEAT model", e.what());
  } catch (...) {
    QMessageBox::warning(this, "Unhandled exception!", "The program may become unstable.");
    //    throw;
  }
}
