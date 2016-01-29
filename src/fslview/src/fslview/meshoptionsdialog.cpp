/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    David Flitney 

    FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "meshoptionsdialog.h"

#include <qlineedit.h>
#include <qcheckbox.h>
#include <qvalidator.h>

//#define DEBUGGING
#include "tracker.h"

using namespace std;

struct MeshOptionsDialog::Implementation
{
  Implementation(): m_selectedMesh(0) { TRACKER("MeshOptionsDialog::Implementation"); }
  ~Implementation() { TRACKER("MeshOptionsDialog::~Implementation"); CHECKPOINT(); }
  
  vector<VTKMeshSurface::Handle> m_meshes;
  unsigned int m_selectedMesh;
};

MeshOptionsDialog::MeshOptionsDialog(QWidget *p):
  MeshOptionsDialogBase(p), m_impl(new Implementation)
{
  TRACKER("MeshOptionsDialog::MeshOptionsDialog");
  m_lower->setValidator(new QDoubleValidator(m_lower));
  m_upper->setValidator(new QDoubleValidator(m_upper));
  m_scaleGlyphFactor->setValidator(new QDoubleValidator(m_upper));
}

MeshOptionsDialog::~MeshOptionsDialog()
{
  TRACKER("MeshOptionsDialog::~MeshOptionsDialog");
}

void MeshOptionsDialog::populateMeshList(vector<VTKMeshSurface::Handle>& m)
{
  m_impl->m_meshes.clear();
  m_impl->m_meshes.reserve(m.size());
  copy(m.begin(), m.end(), back_inserter(m_impl->m_meshes));
  if(m_impl->m_meshes.size())
    {
      VTKMeshSurface::Handle mh(m_impl->m_meshes.at(m_impl->m_selectedMesh));
      if(mh) {
	m_lower->setText( QString("%1").arg(mh->inqLower()) );
	m_upper->setText( QString("%1").arg(mh->inqUpper()) );
	m_showGlyphs->setChecked( mh->inqShowGlyphs() );
	m_showCellData->setChecked( mh->inqShowCellData() );
      }
    }
}

void MeshOptionsDialog::showGlyphs(bool y)
{
  m_impl->m_meshes.at(m_impl->m_selectedMesh)->showGlyphs(y);
}

void MeshOptionsDialog::showCellData(bool y)
{
  m_impl->m_meshes.at(m_impl->m_selectedMesh)->showCellData(y);
}

void MeshOptionsDialog::setMesh(int meshId)
{
  m_impl->m_selectedMesh = meshId;
  VTKMeshSurface::Handle mh(m_impl->m_meshes.at(m_impl->m_selectedMesh));
  if(mh) {
    m_lower->setText( QString("%1").arg(mh->inqLower()) );
    m_upper->setText( QString("%1").arg(mh->inqUpper()) );
    m_showGlyphs->setChecked( mh->inqShowGlyphs() );
    m_showCellData->setChecked( mh->inqShowCellData() );
  }
}

void MeshOptionsDialog::setOpacity(int f)
{
  m_impl->m_meshes.at(m_impl->m_selectedMesh)->setOpacity(f / 100.0);
}

void MeshOptionsDialog::setBounds()
{
  float lower(m_lower->text().toFloat());
  float upper(m_upper->text().toFloat());
  float scale(m_scaleGlyphFactor->text().toFloat());
  m_impl->m_meshes.at(m_impl->m_selectedMesh)->setScalarRange(lower, upper);
  m_impl->m_meshes.at(m_impl->m_selectedMesh)->setGlyphScaleFactor(scale);
}

void MeshOptionsDialog::setWarpFactor(int f)
{
  m_impl->m_meshes.at(m_impl->m_selectedMesh)->setWarpFactor(f / 100.0);
}
