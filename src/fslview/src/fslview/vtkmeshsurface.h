#if !defined(_VTKMESHSURFACE_H)
#define _VTKMESHSURFACE_H

/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Brian Patenaude
                David Flitney 

    FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

class vtkPolyDataReader;
class vtkPolyDataNormals;
class vtkVectorNorm;
class vtkWarpVector;
class vtkPolyDataMapper;
class vtkActor;
class vtkArrowSource;
class vtkGlyph3D;
class vtkScalarBarWidget;
class vtkThresholdPoints;
class vtkRenderer;
class vtkActor;

#include <boost/shared_ptr.hpp>

#include <string>

class VTKMeshSurface
{
public:
  typedef boost::shared_ptr<VTKMeshSurface> Handle;

  static Handle create(vtkRenderer *, const std::string&, float);

  void  showCellData(bool);
  bool  inqShowCellData() const {return m_showCellData; }
  void  showGlyphs(bool);
  bool  inqShowGlyphs() const { return m_showGlyphs; }
  void  setScalarRange(float, float);
  float inqLower() const { return m_lower; }
  float inqUpper() const { return m_upper; }
  void  setWarpFactor(float f);
  float inqWarpFactor() const { return m_warpFactor; }
  void  setOpacity(float f);
  float inqOpacity() const { return m_lower; }
  void  setGlyphScaleFactor(float f);
  float inqGlyphScaleFactor() const { return m_glyphScaleFactor; }
  
  virtual ~VTKMeshSurface();

  vtkActor* getActor() { return m_actor; }
  vtkActor* getGlyphActor() { return m_glyphActor; }

private:
  VTKMeshSurface(vtkRenderer *, const std::string&, float);
  
  void forceRedraw();

  float m_warpFactor, m_glyphScaleFactor;
  float m_opacity;
  float m_lower, m_upper;
  bool m_showGlyphs;
  bool m_showCellData;
  vtkRenderer *m_renderer;

  vtkPolyDataReader *m_model;
  vtkPolyDataNormals *m_modelNormals;
  vtkWarpVector *m_warp;
  vtkPolyDataMapper *m_modelMap, *m_glyphMap;
  vtkActor *m_actor, *m_glyphActor;
  vtkArrowSource *m_arrow;
  vtkGlyph3D *m_glyph;
  vtkScalarBarWidget *m_scalarWidget;
  vtkThresholdPoints *m_thresh;
};

#endif
