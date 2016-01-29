/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    David Flitney 

    FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#ifndef MESHOPTIONSDIALOG_H
#define MESHOPTIONSDIALOG_H

#include "meshoptionsdialogbase.h"
#include "vtkmeshsurface.h"

#include <memory>
#include <vector>

class MeshOptionsDialog : public MeshOptionsDialogBase
{
Q_OBJECT

public:
  MeshOptionsDialog(QWidget*);
  virtual ~MeshOptionsDialog();

  void populateMeshList(std::vector<VTKMeshSurface::Handle>& meshes);

private slots:
  void showGlyphs(bool);
  void showCellData(bool);
  void setMesh(int);
  void setWarpFactor(int);
  void setOpacity(int);
  void setBounds();

private:
  struct Implementation;
  std::auto_ptr<Implementation> m_impl;
};

#endif // MESHOPTIONSDIALOG_H

