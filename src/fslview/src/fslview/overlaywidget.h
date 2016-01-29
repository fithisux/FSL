/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined (OVERLAYWIDGET_H)
#define OVERLAYWIDGET_H

#include <algorithm>
#include <qlistview.h>
#include "overlaywidgetbase.h"
#include "overlaylist.h"
#include "metaimage.h"

class LayerListItem;
class ImageWidget;

class OverlayWidget : public OverlayWidgetBase, OverlayListObserver
{
 Q_OBJECT

 public:
  //depending upon the selected overlaylist item state, set states of certain toolbars
  OverlayWidget(QWidget* w, OverlayList::Handle l);
  virtual ~OverlayWidget();
  void update(const OverlayList* l, OverlayListMsg msg);

  QSize sizeHint() const;

 private:  
  void updateListView();
  void updateListItem();
  void setAllInActive(); 
  void updateUpDownButtons();
  void setVisibility(bool state);
  void updateControls();
  QListViewItem* getLayerItem(const MetaImage::Handle mi)const;

  OverlayList::Handle   m_overlayList; 

  bool           m_blockEvents;
  bool           m_blockSliderUpdate;
  QListViewItem *m_prevItem;
  
 signals:
  virtual void infoButtonAction();

 private slots:

  void listSelectChanged();
  void visibleButtonChanged(bool state);
  void lockedButtonChanged(bool state);
 
  void listDoubleClicked(QListViewItem* item);

  void transSliderChanged(int value);     
  void modTransSliderChanged(int value);
  void upButtonPressed();
  void downButtonPressed();
  void addButtonPressed();
  void removeButtonPressed();
  void detailsButtonPressed();
};


#endif
