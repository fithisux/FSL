#ifndef VIEWOPTIONSDIALOG_H
#define VIEWOPTIONSDIALOG_H

#include "viewoptionsdialogbase.h"

class ViewOptions;

class ViewOptionsDialog: public ViewOptionsDialogBase
{
Q_OBJECT

public:
  ViewOptionsDialog(QWidget *, ViewOptions&);

  ViewOptions getOptions(void) const;

  ~ViewOptionsDialog();
};

#endif //VIEWOPTIONSDIALOG_H
