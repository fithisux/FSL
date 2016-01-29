/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Brian Patenaude
                David Flitney 

    FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "vtkmeshsurface.h"

//
// VTKMeshSurface
//
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkPointData.h>
#include <vtkVectorNorm.h>
#include <vtkWarpVector.h>
#include <vtkGlyph3D.h>
#include <vtkScalarBarActor.h>
#include <vtkScalarBarWidget.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataMapper.h>
#include <vtkArrowSource.h>
#include <vtkThresholdPoints.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

#include <string>
#include <iostream>
using namespace std;

VTKMeshSurface::VTKMeshSurface(vtkRenderer *ren, const string& filename, float xmm):
  m_warpFactor(0), m_glyphScaleFactor(1),
  m_opacity(1.0),
  m_lower(2), m_upper(5),
  m_showGlyphs(true),
  m_showCellData(false),
  m_renderer(ren)
{
  m_model = vtkPolyDataReader::New();
  m_model->SetFileName(filename.c_str());

  vtkTransform* xform = vtkTransform::New();
  xform->Identity();
  xform->Translate(xmm, 0, 0);
  xform->Scale(-1, 1, 1);

  vtkTransformPolyDataFilter* trans = vtkTransformPolyDataFilter::New();
  trans->SetInput(m_model->GetOutput());
  trans->SetTransform(xform);

  m_warp = vtkWarpVector::New();
  m_warp->SetInput(trans->GetOutput());
  m_warp->SetScaleFactor(m_warpFactor);
  
  m_modelNormals = vtkPolyDataNormals::New();
  m_modelNormals->SetInput(m_warp->GetPolyDataOutput());
  m_modelNormals->SetFeatureAngle(60);

  m_modelMap=vtkPolyDataMapper::New();
  m_modelMap->SetInput(m_modelNormals->GetOutput());
  m_modelMap->ScalarVisibilityOn();
  m_modelMap->SetScalarRange(m_lower, m_upper);
  //  m_modelMap->SetScalarModeToUseCellData();
  m_actor = vtkActor::New();
  m_actor->SetMapper(m_modelMap);
  m_actor->GetProperty()->SetInterpolationToGouraud();

  ///************This deals with creation of vector arrows******************//

  m_thresh = vtkThresholdPoints::New();
  m_thresh->SetInput(trans->GetOutput());
  m_thresh->ThresholdByUpper(m_lower);

  //Create GLyphs
  m_arrow = vtkArrowSource::New();
  //Something to turn on off and scale .....
  m_glyph = vtkGlyph3D::New();
  m_glyph->SetInput(m_thresh->GetOutput());
  m_glyph->SetSource(m_arrow->GetOutput());
  m_glyph->SetScaleModeToScaleByVector();
  m_glyph->SetColorModeToColorByScalar();
  m_glyph->SetVectorModeToUseVector();
  //scales overall magintude
  m_glyph->SetScaleFactor(m_glyphScaleFactor);
  m_glyph->ScalingOn();
  m_glyph->SetRange(0,1);
 
  m_glyphMap = vtkPolyDataMapper::New();
  m_glyphMap->SetInput(m_glyph->GetOutput());
  m_glyphMap->ScalarVisibilityOn();	
  ///This should me the same as for the surface
  m_glyphMap->SetScalarRange(m_lower, m_upper);

  m_glyphActor = vtkActor::New();
  m_glyphActor->SetMapper(m_glyphMap);

  //************Creates a colour bar******************//
  m_scalarWidget = vtkScalarBarWidget::New();
  m_scalarWidget->SetInteractor(ren->GetRenderWindow()->GetInteractor());
  m_scalarWidget->GetScalarBarActor()->SetLookupTable(m_modelMap->GetLookupTable());
  m_scalarWidget->GetScalarBarActor()->SetWidth(0.05);
  m_scalarWidget->EnabledOn();
  m_scalarWidget->SetCurrentRenderer(ren);
}

VTKMeshSurface::Handle VTKMeshSurface::create(vtkRenderer *ren, const string& filename, float xmm)
{
  return Handle(new VTKMeshSurface(ren, filename, xmm));
}

VTKMeshSurface::~VTKMeshSurface()
{
  m_actor->Delete();
  m_model->Delete();
  m_modelNormals->Delete();
  m_modelMap->Delete();
  m_glyph->Delete();
  m_glyphMap->Delete();
  m_glyphActor->Delete();
  m_scalarWidget->Delete();
}

void VTKMeshSurface::forceRedraw()
{
  m_renderer->GetRenderWindow()->Render();
}

void VTKMeshSurface::showCellData(bool y)
{
  m_showCellData = y;
  if(m_showCellData)
    m_modelMap->SetScalarModeToUseCellData();
  else
    m_modelMap->SetScalarModeToUsePointData();
  forceRedraw();
}

void VTKMeshSurface::showGlyphs(bool y)
{
  m_showGlyphs = y;
  if(m_showGlyphs)
    m_glyphActor->VisibilityOn();
  else
    m_glyphActor->VisibilityOff();
  forceRedraw();
}

void VTKMeshSurface::setOpacity(float f)
{
  m_actor->GetProperty()->SetOpacity(f);
  m_opacity = f;
  // Add sorting if f < 1.0
  forceRedraw();
}

void VTKMeshSurface::setWarpFactor(float f)
{
  m_warpFactor = f;
  m_warp->SetScaleFactor(m_warpFactor);
  forceRedraw();
}

void VTKMeshSurface::setGlyphScaleFactor(float f)
{
  m_glyphScaleFactor = f;
  m_glyph->SetScaleFactor(m_glyphScaleFactor);
  forceRedraw();
}

void VTKMeshSurface::setScalarRange(float l, float u)
{
  m_upper=u; m_lower=l;

  m_modelMap->SetScalarRange(m_lower, m_upper);
  m_glyphMap->SetScalarRange(m_lower, m_upper);

  m_thresh->ThresholdByUpper(m_lower);

  forceRedraw();
}

