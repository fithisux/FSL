/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined (OVERLAYINFODIALOG_H)
#define OVERLAYINFODIALOG_H

#include "overlaylist.h"
#include "overlayinfodialogbase.h"

class OverlayInfoDialog : public OverlayInfoDialogBase, OverlayListObserver
{
Q_OBJECT

public:
  OverlayInfoDialog(QWidget* w,OverlayList::Handle l,ImageGroup::Handle i);
  virtual ~OverlayInfoDialog();
  void update(const OverlayList* l, OverlayListMsg msg);

private:
  void clearDialog();
  void synchronizeDialog();
  OverlayList::Handle m_overlayList; 
  ImageGroup::Handle  m_imageGroup;

  bool m_blockOverlayListUpdate;
  

protected slots:

   void lutComboChanged(int n);
   void slutComboChanged(int n);
   void slutBoxChecked(bool);
   void dtiComboChanged(int n);
   void modComboChanged(int n);
   void intentChanged(int n);
   void lutButtonPressed(); 
   void overlayTextChanged( const QString & );
   void help();

signals: 

   virtual void openLookUpTable();
   virtual void message(const QString&, int);
};


#endif
