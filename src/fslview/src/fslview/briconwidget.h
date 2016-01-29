/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(BRICONWIDGET_H)
#define BRICONWIDGET_H

#if defined(WIN32)
#pragma warning (disable:4786)
#endif

#include <qwidget.h>
#include <list>

#include "bricon.h"
#include "overlaylist.h"

#include "briconwidgetbase.h"

//! @brief Implentation of BriCon toolbar behaviour
//!
//! @author Dave Flitney
class BriConWidget : public BriConWidgetBase, BriConObserver, OverlayListObserver
{
  Q_OBJECT
public:
  typedef boost::shared_ptr< BriConWidget > Handle;
  
  BriConWidget(QWidget *parent, OverlayList::Handle list);
  virtual ~BriConWidget();

  virtual void update(const BriCon *);
  virtual void update(const OverlayList* list, OverlayListMsg message);
  
  void setMinMaxBoxesState(bool state);
  void setBriSliderState(bool state);
  void setConSliderState(bool state);

private:
  
  BriCon::Handle        m_bricon;
  OverlayList::Handle   m_list;
  //  bool         m_blockEvents;
  float        m_originalMin;
  float        m_originalMax;
  void         updateMinMaxBoxes();

private slots:

  void reset();
  void minChanged();
  void maxChanged();
  void briSliderChanged(double value);
  void conSliderChanged(double value);
  void sliderReleased();
};

#endif
