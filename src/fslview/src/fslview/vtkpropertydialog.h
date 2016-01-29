/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Rama Aravind Vorray
		James Saunders
		David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(VTKPROPERTYDIALOG_H)
#define VTKPROPERTYDIALOG_H

#include "vtkpropertydialogbase.h"

class VTKProperties;

class VTKPropertyDialog: public VTKPropertyDialogBase
{
public:
  VTKPropertyDialog(QWidget* parent, VTKProperties& props);
  VTKPropertyDialog(const VTKPropertyDialog& options);
  VTKProperties& operator=(const VTKProperties& rhs);
  virtual ~VTKPropertyDialog() {}

  VTKProperties& getProperties();

private:
  VTKProperties& m_props;

private slots:
  void selectColor();
  void help();
};

#endif
