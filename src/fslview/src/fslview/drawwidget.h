/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(DRAWWIDGET_H)
#define DRAWWIDGET_H

#if defined(WIN32)
#pragma warning (disable:4786)
#endif

#include <qwidget.h>
#include <list>
#include <qcombobox.h>
#include <qpopupmenu.h>
#include "overlaylist.h"
#include "imagedata.h"
#include "imagedatastore.h"
#include "metaimage.h"
#include "storage/image.h"
#include "lookuptable.h"
#include "bricon.h"
#include "drawtoolbarbase.h"
#include "drawsettings.h"

class QSpinBox;
class DrawWidget;
class QToolButton;

//! @brief Behavioural  implementation of drawing palette toolbar
//! @author Dave Flitney
class DrawWidget : public DrawToolbarBase, 
		   public DrawSettingsObserver,
		   public BriConObserver,
		   public OverlayListObserver
{
  Q_OBJECT
public:

  DrawWidget(QWidget *parent, OverlayList::Handle ol, DrawSettings::Handle ds);
  virtual ~DrawWidget();
 
  void update(const OverlayList* ol, OverlayListMsg msg);
  void update(const BriCon* b);
  void update(const DrawSettings*);

  void updateControls();
  void setCurComboBoxColor(void);
  
signals:
  void undoButtonClicked();
  void redoButtonClicked();

private:
  
  QToolButton *m_lutColInd;
  
  QPopupMenu *popUpMenu;
    
  OverlayList::Handle  m_overlayList;
  BriCon::Handle       m_bricon;
  DrawSettings::Handle m_drawSettings;

  void initialiseLutComboBox();
  void blockControlSignals(bool);
//   void syncLutComboFromValBox();
//   void syncValBoxFromLutCombo();

private slots:
  void sizeBoxChanged(int);
  void valBoxChanged(int);
  void linkButtonToggled(bool);
  void fillButtonToggled(bool);
  void penButtonToggled(bool);
  void eraseButtonToggled(bool);
//   void lutComboActivated(int);
};


#endif
