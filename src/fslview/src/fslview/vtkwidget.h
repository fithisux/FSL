/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Rama Aravind Vorray
		James Saunders
		David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(VTKWIDGET_H)
#define VTKWIDGET_H

#include "imagewidget.h"

#include <boost/shared_ptr.hpp>
#include <vtkProperty.h>

#include <list>

class QVTKWidget;
class VTKPropertiesObserver;

class VTKProperties
{
public:
  typedef boost::shared_ptr< VTKProperties > Handle;

  VTKProperties();
  VTKProperties(const VTKProperties& rhs);

  VTKProperties& operator=(const VTKProperties& rhs);

  void Swap(VTKProperties& other);

  int inqLowerThreshold() const { return m_lowerThreshold; }
  void setLowerThreshold(int t) { m_lowerThreshold = t; notify(); }

  int inqUpperThreshold() const { return m_upperThreshold; }
  void setUpperThreshold(int t) { m_upperThreshold = t; notify(); }

  float inqMcThreshold() const { return m_mcThreshold; }
  void setMcThreshold(float f) { m_mcThreshold = f; notify(); }

  int inqInterpMode() const { return m_interpMode; }
  void setInterpMode(int t) { m_interpMode = t; notify(); }

  float inqStdDev() const 
  { return m_stdDev; }

  void setStdDev(float sd) 
  { m_stdDev = sd; notify(); }

  float inqRadius() const
  { return m_radius; }

  void setRadii(float r)
  { m_radius = r; notify(); }

  void inqColor(float& r, float& g, float& b)
  { r = m_colorR; g = m_colorG; b = m_colorB; }

  void setColor(float r, float g, float b)
  { m_colorR = r; m_colorG = g; m_colorB = b; notify(); }

  float inqRelaxationFactor() const { return m_relaxationFactor; }
  void setRelaxationFactor(float f) { m_relaxationFactor = f; notify(); }

  int inqIterations() const { return m_iterations; }
  void setIterations(int i) { m_iterations = i; notify(); }

  float inqAmbient() const { return m_ambient; }
  void setAmbient(float f) { m_ambient = f; notify(); }

  float inqDiffuse() const { return m_diffuse; }
  void setDiffuse(float f) { m_diffuse = f; notify(); }

  float inqOpacity() const { return m_opacity; }
  void setOpacity(float f) { m_opacity = f; notify(); }

  float inqSpecular() const { return m_specular; }
  void setSpecular(float f) { m_specular = f; notify(); }

  float inqSpecularPower() const { return m_specularPower; }
  void setSpecularPower(float f) { m_specularPower = f; notify(); }

  float inqFeatureAngle() const { return m_featureAngle; }
  void setFeatureAngle(float fa) { m_featureAngle = fa; notify(); }

  void setClipping(bool y) { m_clipping = y; notify(); }
  bool inqClipping() const { return m_clipping; }

  void attach(VTKPropertiesObserver* o);
  void detach(VTKPropertiesObserver* o);
  void notify();

private:
  int m_upperThreshold, m_lowerThreshold;
  float m_mcThreshold;
  int m_interpMode;
  int m_iterations;
  float m_relaxationFactor;
  float m_ambient, m_diffuse, m_opacity, m_specular, m_specularPower;
  float m_featureAngle;
  float m_stdDev;
  float m_radius;
  float m_colorR, m_colorG, m_colorB;
  bool m_clipping;

  std::list<VTKPropertiesObserver*> m_observers;
};

//! @brief interface for any class wishing to observe VTKProperties objects
//!
//! A class which wants to implement VTKPropertiesObserver should subclass itself from
//! VTKPropertiesObserver and implement the VTKPropertiesObserver::update method.
class VTKPropertiesObserver
{
public:
  VTKPropertiesObserver() {}

  virtual void update(const VTKProperties*) = 0;

  virtual ~VTKPropertiesObserver() {}
};

class VTKWidget : public ImageWidget, public VTKPropertiesObserver, public BriConObserver
{
  Q_OBJECT
public:
  VTKWidget(QWidget *parent, 
	    ImageGroup::Handle i, 
	    OverlayList::Handle ol,
	    Cursor::Handle c);
  virtual ~VTKWidget();
  virtual void update(const VTKProperties*);
  virtual void update(const BriCon*);
  virtual void update(const Cursor::Handle);
  virtual void update(const OverlayList*, OverlayListMsg);

private slots:
  void print();
  void options();
  void addMesh();
  void meshOptions();
  

private:
  struct Implementation;  
  const std::auto_ptr<Implementation> m_impl;

  ImageGroup::Handle m_image;
  QVTKWidget        *m_vtkwidget;
};

#endif
