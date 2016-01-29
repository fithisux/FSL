/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "overlaywidget.h"

#include <qapplication.h>
#include <qlayout.h>
#include <qcheckbox.h>
#include <qlistview.h>
#include <qslider.h>
#include <qtoolbutton.h>
#include <qobject.h>
#include <qheader.h>
#include <qlabel.h>
#include <qtooltip.h>

#include "tracker.h"

#include "mainimage.xpm"
//#include "info.xpm"
#include "eye.xpm"
#include "padlock.xpm"
//#include "uparrow.xpm"
//#include "downarrow.xpm"

const char * toggleVisibiltyText = "Visibility on/off.<br><hr>Hides/shows the currently highlighted overlay.";

const char * toggleLockText = "Lock on/off.<br><hr>Images must be unlocked before they are edited.";

const char * transparencySliderText = "Transparency slider.<br><hr>Controls the transparency of the currently selected overlay.";

const char * modTransparencySliderText = "Modulation transparency slider.<br><hr>Controls the effects of image modulation on transparency.";

const char * infoButtonText = "Open information dialog.<br><hr>Shows details of the currently selected image and allows the look up table to be changed.";

const char * upButtonText = "Move overlay up.<br><hr>Moves the currently selected overlay so that it appears on top of other overlays.";

const char * downButtonText = "Move overlay down.<br><hr>Moves the currently selected overlay so that it appears below other overlays.";

class LayerListItem : public QListViewItem
{
public:
  LayerListItem(QListView* parent, const  MetaImage::Handle mi):
    QListViewItem(parent),
    m_mi(mi)
  {
    refresh();
  }
  
  void refresh()
  {
    if(inqVisibility())
      setPixmap(0,QPixmap(eye));
    else 
      setPixmap(0,NULL);
    if(inqReadOnly())
      setPixmap(1,QPixmap(padlock));
    else 
      setPixmap(1,NULL);
    setText(2,inqImageName().c_str());
    setText(3,QString::number(inqTransparency()));
    if(inqIsMainImage())
      setPixmap(4,QPixmap(mainimage));
    else 
      setPixmap(4,NULL);
  }
  
  const MetaImage::Handle getMetaImage() { return m_mi;}
  bool inqIsMainImage()          { return m_mi->getInfo()->isMainImage();}
  std::string inqImageName()     { return m_mi->inqImageName();}
  float inqTransparency()        { return m_mi->getDs()->inqTransparency();}  
  bool  inqVisibility()          { return m_mi->inqVisibility();}
  bool  inqReadOnly()            { return m_mi->inqReadOnly();}
  bool  inqTransMod()            { return m_mi->getDs()->inqTransMod();}
  float inqModTransparency()     { return m_mi->getDs()->inqModTransparency();}

private:
  const MetaImage::Handle m_mi;

};

class InsertLayerItem {
public:
  InsertLayerItem(QListView* lv) : m_lv(lv) {}
  void operator()(MetaImage::Handle mi);

private:
  QListView* m_lv;
};

void InsertLayerItem::operator()(MetaImage::Handle mi) 
{
   new LayerListItem(m_lv, mi);
}

OverlayWidget::OverlayWidget(QWidget* w, OverlayList::Handle l):
  OverlayWidgetBase(w, 0), m_overlayList(l), m_blockEvents(false), m_blockSliderUpdate(false)
{
  TRACKER("OverlayWidget::OverlayWidget(QWidget* w, OverlayList::Handle l)");
  m_overlayList->attach(this);

  m_transSlider->setMinValue(0);
  m_transSlider->setMaxValue(10);
  m_transSlider->setPageStep(1); 

  m_modTransSlider->setMinValue(0);
  m_modTransSlider->setMaxValue(10);
  m_modTransSlider->setPageStep(1);
  m_modTransSlider->hide();
  
  m_addButton->hide();
  m_removeButton->hide();

  m_overlayListView->header()->hide();
  m_overlayListView->setSorting(-1);
  m_overlayListView->setAllColumnsShowFocus(true);

  connect(m_upButton, SIGNAL(pressed()), this, SLOT(upButtonPressed()));
  connect(m_downButton, SIGNAL(pressed()), this, SLOT(downButtonPressed()));
  connect(m_detailsButton, SIGNAL(pressed()), this, SLOT(detailsButtonPressed()));

  connect(m_overlayListView,SIGNAL(selectionChanged()),          
          this ,SLOT(listSelectChanged()));
  connect(m_overlayListView,SIGNAL(doubleClicked(QListViewItem*)),
        this,SLOT(listDoubleClicked(QListViewItem*)));

  connect(m_visibleButton,SIGNAL(toggled(bool)),
        this,SLOT(visibleButtonChanged(bool)));  
  connect(m_lockedButton,SIGNAL(toggled(bool)),
        this,SLOT(lockedButtonChanged(bool)));

  connect(m_transSlider, SIGNAL(valueChanged(int)),
        this,SLOT(transSliderChanged(int)));  
  connect(m_modTransSlider, SIGNAL(valueChanged(int)),
        this,SLOT(modTransSliderChanged(int)));

//   connect(m_addButton, SIGNAL(pressed()), qApp, SLOT(addOverlay()));
//   connect(m_removeButton, SIGNAL(pressed()), qApp, SLOT(remOverlay()));

  QToolTip::add(m_visibleButton,tr(toggleVisibiltyText));  
  QToolTip::add(m_lockedButton, tr(toggleLockText));
  QToolTip::add(m_detailsButton,  tr(infoButtonText));
  QToolTip::add(m_upButton,     tr(upButtonText));
  QToolTip::add(m_downButton,   tr(downButtonText));
  QToolTip::add(m_transSlider,  tr(transparencySliderText));
  QToolTip::add(m_modTransSlider,  tr(modTransparencySliderText));

  updateListView();
  updateUpDownButtons();
  updateControls();
}

OverlayWidget::~OverlayWidget()
{
  TRACKER("OverlayWidget::~OverlayWidget()");
  m_overlayList->detach(this);
}

QSize OverlayWidget::sizeHint() const
{
  TRACKER("OverlayWidget::sizeHint()");
  return QSize(225, 65);
}

void OverlayWidget::addButtonPressed()
{
  TRACKER("OverlayWidget::addButtonPressed()");
  //  qApp->addOverlay();
}

void OverlayWidget::removeButtonPressed()
{
  TRACKER("OverlayWidget::removeButtonPressed()");
  //  qApp->remOverlay();
}

void OverlayWidget::updateListView()
{
  TRACKER("OverlayWidget::updateListView()");

  m_overlayListView->clear();

  std::for_each(m_overlayList->begin(),
                m_overlayList->end(),
                InsertLayerItem(m_overlayListView));
  
  QListViewItem* curItem = getLayerItem(m_overlayList->getActiveMetaImage());

  if(curItem)
    m_overlayListView->setSelected(curItem, true);
}

void OverlayWidget::update(const OverlayList* l, OverlayListMsg msg)
{
  TRACKER("OverlayWidget::update(const OverlayList* l, OverlayListMsg msg)");

  m_blockEvents = true;

  switch(msg)
  {
   case OverlayListMsg(Select):       updateListItem();updateControls();break;
   case OverlayListMsg(Visibility):   
   case OverlayListMsg(Transparency): updateListItem();updateControls();break;
   case OverlayListMsg(Security):
   case OverlayListMsg(Order):
   case OverlayListMsg(Add):
   case OverlayListMsg(Rem):          updateListView();updateControls();
                                      updateUpDownButtons();break;
   case OverlayListMsg(LookUpTable):  break;   
   case OverlayListMsg(ModImage):     updateListView();updateControls();break;
   case OverlayListMsg(ImageName):    updateListView();break;
  default: break;
  }
  m_blockEvents = false;
}

void OverlayWidget::listSelectChanged()
{
  TRACKER("OverlayWidget::listSelectChanged");
  if(m_blockEvents) return;

  LayerListItem* i = (LayerListItem*)m_overlayListView->selectedItem();
  MetaImage::Handle mi;

  if(i){

    mi = i->getMetaImage();

    updateUpDownButtons();

    m_overlayList->setActiveMetaImage(mi);
    
  }
  else 
  {
    m_overlayList->setActiveMetaImage(m_overlayList->getMainMetaImage());
  }
}

QListViewItem* OverlayWidget::getLayerItem(const MetaImage::Handle mi) const
{
  TRACKER("OverlayWidget::getLayerItem(const MetaImage::Handle mi)");
  bool success(false);
  QListViewItem* itemFound = NULL;
  QListViewItemIterator it( m_overlayListView );
    
  while(it.current() && !success)   
  {
    if(((LayerListItem*)it.current())->getMetaImage() == mi)
      {
      itemFound = it.current(); 
      success = true;
      }

    ++it;
  }

  return itemFound;
}

void OverlayWidget::transSliderChanged(int value)
{
  TRACKER("OverlayWidget::transSliderChanged(int value)");
  if(m_blockEvents) return;
  m_blockSliderUpdate = true;  
  m_overlayList->setTransparency((float)value/10.0f);
  m_blockSliderUpdate = false;
}

void OverlayWidget::modTransSliderChanged(int value)
{  
  TRACKER("OverlayWidget::modTransSliderChanged(int value)");
  if(m_blockEvents) return;
  m_blockSliderUpdate = true;
  m_overlayList->setModTransparency((float)value/10.0f);
  m_blockSliderUpdate = false;
}

void OverlayWidget::visibleButtonChanged(bool state)
{  
  TRACKER("OverlayWidget::visibleButtonChanged(bool state)");
  if(m_blockEvents) return;
  m_overlayList->setVisibility(state);
  //if(!state)
  //m_overlayListView->clearSelection();
}

void OverlayWidget::lockedButtonChanged(bool state)
{  
  TRACKER("OverlayWidget::lockedButtonChanged(bool state)");
  if(m_blockEvents) return;
  m_overlayList->setReadOnly(state);
}

void OverlayWidget::listDoubleClicked(QListViewItem* item)
{
  TRACKER("OverlayWidget::listDoubleClicked(QListViewItem* item)");
  if(m_blockEvents) return; 
 
  bool state;
   
  state = m_overlayList->getActiveMetaImage()->inqVisibility();
  m_overlayList->setVisibility(!state);

  // clears the list view item (overlays) selection thereby forcing user to select an overlay
  if(state)
    m_overlayListView->clearSelection();
}

void OverlayWidget::upButtonPressed()
{
  TRACKER("OverlayWidget::upButtonPressed()");
   m_overlayList->moveOverlayUp();
}

void OverlayWidget::downButtonPressed()
{
  TRACKER("OverlayWidget::downButtonPressed()");
   m_overlayList->moveOverlayDown();
}

void OverlayWidget::updateUpDownButtons()
{   
  TRACKER("OverlayWidget::updateUpDownButtons()");
  LayerListItem* i = (LayerListItem*)m_overlayListView->selectedItem();
 
  if(i)
    {
      if(i->itemAbove()) 
        m_upButton->setEnabled(true);
      else
        m_upButton->setEnabled(false);

      if(i->itemBelow()) 
        m_downButton->setEnabled(true);
      else
        m_downButton->setEnabled(false);
    }
  else
    {
      m_upButton->setEnabled(false);
      m_downButton->setEnabled(false);
    }

}

void OverlayWidget::updateListItem()
{
  TRACKER("OverlayWidget::updateListItem()");
  LayerListItem* curOverlay = (LayerListItem*)getLayerItem(m_overlayList->getActiveMetaImage()); 
  if(curOverlay) 
    curOverlay->refresh();
}

void OverlayWidget::updateControls()
{
  TRACKER("OverlayWidget::updateControls()");    
  
  LayerListItem* curOverlay = (LayerListItem*)getLayerItem(m_overlayList->getActiveMetaImage()); 

       if(curOverlay)
         {
           if(!m_blockSliderUpdate)
             m_transSlider->setValue((int)(curOverlay->inqTransparency() * 10));
           
           m_visibleButton->setChecked(curOverlay->inqVisibility());
           m_lockedButton->setChecked(curOverlay->inqReadOnly());
          
           if(!m_blockSliderUpdate) 
             {
               m_modTransSlider->setValue((int)(curOverlay->inqModTransparency() * 10));

               if(curOverlay->inqTransMod()){m_modTransSlider->show();}
               else                         {m_modTransSlider->hide();}
             }
           
         }
}

void OverlayWidget::detailsButtonPressed()
{
  TRACKER("OverlayWidget::detailsButtonPressed()");
  emit infoButtonAction();
}
