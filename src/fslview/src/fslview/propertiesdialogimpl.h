#if !defined(PROPERTIESDIALOGIMPL_H)
#define PROPERTIESDIALOGIMPL_H

/* CCOPYRIGHT */

#include "propertiesdialog.h"

class PropertiesDialogImpl : public PropertiesDialog
{
public:
  static void getProperties(QWidget *);

private:
  void commit();

  PropertiesDialogImpl(QWidget *);
  ~PropertiesDialogImpl();
};
  
#endif
