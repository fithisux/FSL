/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Brian Patenaude
		David Flitney 
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#if defined(WIN32)
#pragma warning(disable:4786)
#endif

#include <vector>
#include <map>

#include <qapplication.h>
#include <qfiledialog.h>

#include "vtkwidget.h"
#include "maintoolbar.h"
#include "modetoolbar.h"
#include "metaimage.h"

//#define DEBUGGING
#include "tracker.h"

#include <vtkImageActor.h>
#include <vtkVolume.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkLookupTable.h>
#include <vtkImageBlend.h>
#include <vtkImagePlaneWidget.h>
#include <vtkImageMapToColors.h>
#include <vtkImageOpenClose3D.h>
#include <vtkVolumeRayCastCompositeFunction.h>
#include <vtkVolumeRayCastMapper.h>
#include <vtkColorTransferFunction.h>
#include <vtkMarchingCubes.h>
#include <vtkDiscreteMarchingCubes.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkClipPolyData.h>
#include <vtkClipVolume.h>
#include <vtkCutter.h>
#include <vtkDecimatePro.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataWriter.h>
#include <vtkDataSetWriter.h>
#include <vtkBoxWidget.h>
#include <vtkPlaneSource.h>
#include <vtkPlanes.h>
#include <vtkPlane.h>
#include <vtkImageThreshold.h>
#include <vtkCellDataToPointData.h>

//#include <vtkImageReslice.h>
#include <vtkProperty.h>
#include <vtkVolumeProperty.h>
#include <vtkPiecewiseFunction.h>
#include <vtkScalarsToColors.h>
#include <vtkTextureMapToPlane.h>
#include <vtkTransform.h>
#include <vtkDataSetMapper.h>
#include <vtkImageMapper.h>
#include <vtkProbeFilter.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkImageStencil.h>
#include <vtkImplicitDataSet.h>
#include <vtkDataSetWriter.h>
#include <QVTKWidget.h>

#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkDataSetAttributes.h>
#include <vtkPolyDataNormals.h>
#include <vtkTensorGlyph.h>
#include <vtkSphereSource.h>

#include <vtkImageGaussianSmooth.h>

#include <vtkSmoothPolyDataFilter.h>
#include <vtkImageDilateErode3D.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkDepthSortPolyData.h>
#include <vtkPolyDataSource.h>
#include <vtkWindowToImageFilter.h>
#include <vtkTIFFWriter.h>

#include <vtkPolyDataConnectivityFilter.h>
#include <vtkAppendPolyData.h>

#include "vtktoolbar.h"
#include "vtkpropertydialog.h"
#include "meshoptionsdialog.h"
using namespace std;

class LookUpTableFactory
{
public:
  static vtkLookupTable* convert(LookUpTable::Handle lh)
  {
    vtkLookupTable* ctfun = vtkLookupTable::New();

    ctfun->SetNumberOfTableValues(lh->size());
    LookUpTable::SizeType i(0);
    for(LookUpTable::ConstIterator it=lh->begin(); it!=lh->end(); ++it, ++i) {
      float r = it->red() / 255.0;
      float g = it->green() / 255.0;
      float b = it->blue() / 255.0;
      float a = (i == 0) ? 0 : it->alpha() / 255.0;
      ctfun->SetTableValue(int(i), r, g, b, a);
    }

    ctfun->SetTableRange(0, lh->size());
    ctfun->Build();

    return ctfun;
  }

  static vtkLookupTable* GetLookupTableByName(const std::string& name)
  {
    LookUpTable::Handle lh;

    if(name == "Greyscale")      lh = LookUpTable::greyScale();
    if(name == "Red-Yellow")     lh = LookUpTable::redYellow();
    if(name == "Blue-Lightblue") lh = LookUpTable::blueLightblue();
    if(name == "Red")            lh = LookUpTable::red();
    if(name == "Green")          lh = LookUpTable::green();
    if(name == "Blue")           lh = LookUpTable::blue();
    if(name == "Pink")           lh = LookUpTable::pink();
    if(name == "Hot")            lh = LookUpTable::hot();
    if(name == "Copper")         lh = LookUpTable::copper();
    if(name == "Cool")           lh = LookUpTable::cool();
    if(name == "MGH-Cortical")   lh = LookUpTable::cortical();
    if(name == "MGH-Subcortical")lh = LookUpTable::subcortical();
    if(name == "Random-Rainbow") lh = LookUpTable::rainbow();

    if(!lh)
      lh = LookUpTable::greyScale();

    return convert(lh);
  }

private:
  mutable unsigned int m_index;
  
  std::vector<vtkLookupTable*> m_luts;
};

#include "vtkmeshsurface.h"

class InteractorCallback: public vtkCommand
{
public:
  static InteractorCallback* New()
  { return new InteractorCallback; }
  virtual void Execute(vtkObject* caller, unsigned long eventid, void*)
  {
    //cout << "An event!" << endl;
    switch(eventid) {
    case vtkCommand::StartInteractionEvent:
      //cout << "Start event" << endl;
      if(m_props->inqOpacity() != 1.0)
	{
	  cout << "Disconnect sorter" << endl;
	  m_mapper->SetInput(m_clipper->GetOutput());
	}
      break;
    case vtkCommand::EndInteractionEvent:
      //cout << "End event" << endl;
      if(m_props->inqOpacity() != 1.0)
	{
	  cout << "Connect sorter" << endl;
	  m_mapper->SetInput(m_sorter->GetOutput());
	}
      break;
    default:
      break;
    }
  }
  void SetProperties(VTKProperties* p) { m_props = p; }
  void SetClipper(vtkClipPolyData* c) { m_clipper = c; }
  void SetSorter(vtkDepthSortPolyData* s) { m_sorter = s; }
  void SetMapper(vtkPolyDataMapper* m) { m_mapper = m; }

private:
  
  vtkDepthSortPolyData* m_sorter;
  vtkClipPolyData* m_clipper;
  vtkPolyDataMapper* m_mapper;
  VTKProperties* m_props;
};

class ClippingBoxCallback: public vtkCommand
{
public:
  static ClippingBoxCallback* New()
  { return new ClippingBoxCallback; }
  virtual void Execute(vtkObject* caller, unsigned long eventid, void*)
  {
    vtkBoxWidget* box = reinterpret_cast<vtkBoxWidget*>(caller);
    box->GetPlanes(m_planes);
    box->GetTransform(m_xform);
  }
  void SetTransform(vtkTransform* t) { m_xform = t; }
  void SetPlanes(vtkPlanes* p) { m_planes = p; }

private:
  
  vtkPlanes* m_planes;
  vtkTransform* m_xform;
};

VTKProperties::VTKProperties(): 
  m_upperThreshold(32000), 
  m_lowerThreshold(80), 
  m_mcThreshold(0.999),
  m_interpMode(VTK_GOURAUD),
  m_iterations(5),
  m_relaxationFactor(0.2), 
  m_ambient(0.1), 
  m_diffuse(0.9), 
  m_opacity(1.0),
  m_specular(0.1),
  m_specularPower(20),
  m_featureAngle(179.0),
  m_stdDev(0.1), m_radius(1.0),
  m_colorR(0.80), m_colorG(0.80), m_colorB(0.80),
  m_clipping(false)
{ 
} 

VTKProperties::VTKProperties(const VTKProperties& rhs): 
  m_upperThreshold(rhs.m_upperThreshold),
  m_lowerThreshold(rhs.m_lowerThreshold),
  m_mcThreshold(rhs.m_mcThreshold),
  m_interpMode(rhs.m_interpMode),
  m_iterations(rhs.m_iterations),
  m_relaxationFactor(rhs.m_relaxationFactor),
  m_ambient(rhs.m_ambient),
  m_diffuse(rhs.m_diffuse),
  m_opacity(rhs.m_opacity),
  m_specular(rhs.m_specular),
  m_specularPower(rhs.m_specularPower),
  m_featureAngle(rhs.m_featureAngle),
  m_stdDev(rhs.m_stdDev), m_radius(rhs.m_radius),
  m_colorR(rhs.m_colorR), m_colorG(rhs.m_colorG), m_colorB(rhs.m_colorB),
  m_clipping(rhs.m_clipping)
{
}

VTKProperties& VTKProperties::operator=(const VTKProperties& rhs)
{
  VTKProperties temp(rhs);
  Swap(temp);
  return *this;
}

void VTKProperties::Swap(VTKProperties& other)
{
  std::swap(m_lowerThreshold, other.m_lowerThreshold);
  std::swap(m_upperThreshold, other.m_upperThreshold);
  std::swap(m_mcThreshold, other.m_mcThreshold);
  std::swap(m_interpMode, other.m_interpMode);
  std::swap(m_relaxationFactor, other.m_relaxationFactor);
  std::swap(m_iterations, other.m_iterations);
  std::swap(m_ambient, other.m_ambient);
  std::swap(m_diffuse, other.m_diffuse);
  std::swap(m_opacity, other.m_opacity);
  std::swap(m_specular, other.m_specular);
  std::swap(m_specularPower, other.m_specularPower);
  std::swap(m_featureAngle, other.m_featureAngle);
  std::swap(m_colorR, other.m_colorR);
  std::swap(m_colorG, other.m_colorG);
  std::swap(m_colorB, other.m_colorB);
  std::swap(m_stdDev, other.m_stdDev);
  std::swap(m_radius, other.m_radius);
  std::swap(m_clipping, other.m_clipping);
}

void VTKProperties::attach(VTKPropertiesObserver* o)
{
  m_observers.remove(o);
  m_observers.push_back(o);
}

void VTKProperties::detach(VTKPropertiesObserver* o)
{
  m_observers.remove(o);
}

struct Update
{
  Update(VTKProperties* p): m_p(p) {}

  void operator()(VTKPropertiesObserver* v)
  {
    v->update(m_p);
  }

  VTKProperties* m_p;
};

void VTKProperties::notify()
{
  TRACKER("VTKProperties::notify()");
  MESSAGE(QString("Notifying %1 observers").arg(m_observers.size()));

  std::for_each(m_observers.begin(), m_observers.end(), Update(this));
}

template <typename T>
struct DeleteVTKObject:
  public std::unary_function<T*, void> {

  void operator()(T* ptr) const
  {
    ptr->Delete();
  }
};

class ImagePipelineObject: public BriConObserver
{
public:
  typedef boost::shared_ptr<ImagePipelineObject> Handle;

  ImagePipelineObject(MetaImage::Handle mi, vtkRenderer *ren): 
    m_smooth(vtkSmoothPolyDataFilter::New()),
//     m_sorter(vtkDepthSortPolyData::New()),
    m_normals(vtkPolyDataNormals::New()),
    m_bricon(mi->getDs()->inqBriCon()), 
    m_metaimage(mi), m_layerMapper(0),
    m_ren(ren)
  {

    Image::Handle im(mi->getImage());
    ImageInfo::Handle ii(mi->getInfo());

    m_bricon->attach(this);

    m_imageData = vtkImageData::New();
    m_imageData->SetDimensions(ii->inqX(), ii->inqY(), ii->inqZ());
    m_imageData->SetSpacing(ii->inqXDim(), ii->inqYDim(), ii->inqZDim());
    m_imageData->SetScalarTypeToFloat();
    m_imageData->SetNumberOfScalarComponents(1);
    m_imageData->AllocateScalars();

    float *ptr = (float *)m_imageData->GetScalarPointer();
    for(int z = 0; z < ii->inqZ(); ++z)
      for(int y = 0; y < ii->inqY(); ++y)
	for(int x = 0; x < ii->inqX(); ++x)
	  // 
	  // 	  NB. X inversion to view 3D in neurological convention
	  // 	  since data is in radiological order
	  // 
	  *ptr++ = im->getVolume(0)->value(ii->inqX() - 1 - x, y, z);

    m_surfLut = LookUpTableFactory::convert(LookUpTable::greyScale());
    m_imageLut = LookUpTableFactory::convert(LookUpTable::greyScale());

    m_layerRGBA = vtkImageMapToColors::New();
    m_layerRGBA->SetOutputFormatToRGBA();
    m_layerRGBA->SetInput(m_imageData);
    m_layerRGBA->SetLookupTable(m_imageLut);
  }
  
  virtual void setSurfaceLut(vtkLookupTable* lut) { m_surfLut = lut; }
  virtual vtkLookupTable* getSurfaceLut() const { return m_surfLut; }
  virtual void setImageLut(vtkLookupTable* lut) 
  { 
    m_imageLut = lut;
    m_layerRGBA->SetLookupTable(m_imageLut);
    update(m_bricon.get());
  }
  virtual vtkLookupTable* getImageLut() const { return m_imageLut; }

  virtual ~ImagePipelineObject()
  {
    m_bricon->detach(this);
    m_smooth->Delete();
    m_normals->Delete();
//     m_sorter->Delete();
  }

  virtual void update(const BriCon* bc) {}

  virtual vtkMapper *getLayerMapper()
  { 
    if(!m_layerMapper)
      m_layerMapper = vtkPolyDataMapper::New();
    m_layerMapper->SetInput(GetSurface());
    m_layerMapper->SetColorModeToMapScalars();
    m_layerMapper->UseLookupTableScalarRangeOn();
    m_layerMapper->SetLookupTable(getSurfaceLut());
    m_layerMapper->ScalarVisibilityOn();
    return m_layerMapper;
  }

//   virtual void setSorting(bool sorting)
//   {
//     if(sorting)
//       m_normals->SetInput(m_smooth->GetOutput());
//     else
//       m_normals->SetInput(m_sorter->GetOutput());
//   }

  void Render() { m_ren->GetRenderWindow()->Render(); }

  //virtual vtkImageData *GetThreshOutput() = 0;
  virtual vtkPolyData* GetSurface() { return m_normals->GetOutput(); }
  //virtual vtkPolyData* GetSortedSurface() = 0;

  virtual void setProperties(const VTKProperties& props) {}

  virtual vtkImageData *GetOutputRGBA() { return m_layerRGBA->GetOutput(); }
  virtual vtkImageData *GetOutput() { return m_imageData; }

protected:
  vtkSmoothPolyDataFilter *m_smooth;
//   vtkDepthSortPolyData *m_sorter;
  vtkPolyDataNormals *m_normals;
  BriCon::Handle m_bricon;
  MetaImage::Handle m_metaimage;
  vtkLookupTable *m_surfLut, *m_imageLut;
  vtkPolyDataMapper *m_layerMapper;
  vtkImageData *m_imageData;
  vtkImageMapToColors *m_layerRGBA;
  vtkRenderer *m_ren;
};

class StatsImage: public ImagePipelineObject
{
public:
  StatsImage(MetaImage::Handle mi, vtkRenderer* ren):
    ImagePipelineObject(mi, ren),
    m_th(vtkImageThreshold::New())
  {
    LookUpTable::Handle lh(mi->getDs()->inqLookUpTable());
    
    if(lh) {
      MESSAGE(std::string("Searching for LUT:") + lh->inqLutName());
      setImageLut(LookUpTableFactory::convert(lh));
      setSurfaceLut(LookUpTableFactory::convert(lh));
    } else {
      MESSAGE(std::string("Using LUT: GreyScale"));
      setImageLut(LookUpTableFactory::convert(LookUpTable::greyScale()));
      setSurfaceLut(LookUpTableFactory::convert(LookUpTable::greyScale()));
    }
    BriCon::Handle bc(mi->getDs()->inqBriCon());
    m_imageLut->SetRange(bc->inqMin(), bc->inqMax());

    m_layerRGBA->SetLookupTable(m_imageLut);
 
    // And IsoSurface
    m_th->SetInput(m_imageData);
    m_th->ThresholdByUpper(bc->inqMin());
    m_th->SetOutValue(0.0);
    m_th->SetInValue(1.0);
    m_th->ReplaceOutOn();
    m_th->ReplaceInOn();
    m_th->SetOutputScalarTypeToFloat();
    m_th->UpdateInformation();

//     vtkImageDilateErode3D *mcdil = vtkImageDilateErode3D::New();
//     mcdil->SetInput(m_th->GetOutput());
//     mcdil->SetKernelSize(2, 2, 2);
//     mcdil->SetDilateValue(1);
//     mcdil->SetErodeValue(0);
	  
    vtkMarchingCubes *mclayer = vtkMarchingCubes::New();
    mclayer->SetInput(m_th->GetOutput());
    mclayer->SetValue(0,0.5);
    mclayer->ComputeNormalsOn();

    vtkDecimatePro *mcdecim = vtkDecimatePro::New();
    mcdecim->SetInput(mclayer->GetOutput());
    mcdecim->PreserveTopologyOn();
    mcdecim->SplittingOn();
    mcdecim->BoundaryVertexDeletionOn();

    m_smooth->SetInput(mcdecim->GetOutput());
    m_smooth->SetNumberOfIterations(5);
    m_smooth->SetFeatureAngle(150);
    m_smooth->SetFeatureEdgeSmoothing(true);
    m_smooth->SetBoundarySmoothing(true);
    m_smooth->SetRelaxationFactor(0.5);

    m_normals->SetInput(m_smooth->GetOutput());
    m_normals->SplittingOn();
    m_normals->ConsistencyOn();
    m_normals->NonManifoldTraversalOff();
    m_normals->SetFeatureAngle(150);

//     m_sorter->SetInput(m_smooth->GetOutput());
//     m_sorter->SetCamera(ren->GetActiveCamera());
//     m_sorter->SortScalarsOn();
//     m_sorter->SetDirectionToBackToFront();
  }

  virtual vtkMapper *getLayerMapper()
  { 
    if(!m_layerMapper)
      m_layerMapper = vtkPolyDataMapper::New();
    m_layerMapper->SetInput(GetSurface());
    m_layerMapper->SetColorModeToMapScalars();
    m_layerMapper->UseLookupTableScalarRangeOn();
    getSurfaceLut()->SetRange(0, 1);
    m_layerMapper->SetLookupTable(getSurfaceLut());
    m_layerMapper->ScalarVisibilityOn();
    
    return m_layerMapper;
  }

  virtual void update(const BriCon* bc)
  {
    m_th->ThresholdByUpper(bc->inqMin());    
    m_imageLut->SetRange(bc->inqMin(), bc->inqMax());
    Render();
  }

private:
  vtkImageThreshold *m_th;
};

class MaskImage: public ImagePipelineObject
{
public:
  MaskImage(MetaImage::Handle mi, vtkRenderer* ren):
    ImagePipelineObject(mi, ren) 
  {
    LookUpTable::Handle lh(mi->getDs()->inqLookUpTable());
    
    if(lh) {
      MESSAGE(std::string("Searching for LUT:") + lh->inqLutName());
      setImageLut(LookUpTableFactory::convert(lh));
      setSurfaceLut(LookUpTableFactory::convert(lh));
    } else {
      MESSAGE(std::string("Using LUT: GreyScale"));
      setImageLut(LookUpTableFactory::convert(LookUpTable::greyScale()));
      setSurfaceLut(LookUpTableFactory::convert(LookUpTable::greyScale()));
    }
    BriCon::Handle bc(mi->getDs()->inqBriCon());
    m_imageLut->SetRange(bc->inqMin(), bc->inqMax());

    m_layerRGBA->SetLookupTable(m_imageLut);
 
    // And IsoSurface
    m_th = vtkImageThreshold::New();
    m_th->SetInput(m_imageData);
    m_th->ThresholdByUpper(0.5);
    m_th->ReplaceOutOn();
    m_th->SetOutputScalarTypeToFloat();
    m_th->SetOutValue(0.0);
    m_th->UpdateInformation();

    vtkDiscreteMarchingCubes *mclayer = vtkDiscreteMarchingCubes::New();
    mclayer->SetInput(m_th->GetOutput());
    unsigned int nlayers(int(bc->inqMax()) - int(bc->inqMin()));
    mclayer->GenerateValues(nlayers, 1, nlayers);

    vtkCellDataToPointData *c2p = vtkCellDataToPointData::New();
    c2p->SetInput(mclayer->GetOutput());
    c2p->PassCellDataOn();
    c2p->UpdateInformation();

    vtkDecimatePro *mcdecim = vtkDecimatePro::New();
    mcdecim->SetInput(c2p->GetPolyDataOutput());
    mcdecim->PreserveTopologyOn();
    mcdecim->SplittingOn();
    mcdecim->BoundaryVertexDeletionOn();

    m_smooth->SetInput(mcdecim->GetOutput());
    m_smooth->SetNumberOfIterations(2);
//     m_smooth->SetFeatureAngle(150);
//     m_smooth->SetFeatureEdgeSmoothing(true);
//     m_smooth->SetBoundarySmoothing(false);
//     m_smooth->SetRelaxationFactor(0.1);

    m_normals->SetInput(m_smooth->GetOutput());
    m_normals->SplittingOn();
    m_normals->ConsistencyOn();
    m_normals->NonManifoldTraversalOff();
    m_normals->SetFeatureAngle(150);
    
//     m_sorter->SetInput(m_smooth->GetOutput());
//     m_sorter->SetCamera(ren->GetActiveCamera());
//     m_sorter->SortScalarsOn();
//     m_sorter->SetDirectionToBackToFront();
  }

  vtkImageData* GetThreshOutput() { return m_th->GetOutput(); }

private:
  vtkImageThreshold *m_th;
};

class MainImage: public ImagePipelineObject
{
public:
  typedef boost::shared_ptr<MainImage> Handle;

  MainImage(MetaImage::Handle mi, const VTKProperties& props, vtkRenderer* ren):
    ImagePipelineObject(mi, ren)
  {
    m_thresh = vtkImageThreshold::New();
    m_thresh->SetInput(m_imageData);
    m_thresh->ThresholdByUpper(props.inqLowerThreshold());
    m_thresh->SetOutputScalarTypeToFloat();
    m_thresh->ReplaceOutOn();
    m_thresh->SetOutValue(0);
    m_thresh->ReplaceInOff();
    m_thresh->UpdateInformation();

    float sd(props.inqStdDev());
    float r(props.inqRadius());

    m_gaussian = vtkImageGaussianSmooth::New();
    m_gaussian->SetInput(m_thresh->GetOutput());
    m_gaussian->SetDimensionality(3);
    m_gaussian->SetStandardDeviations(sd, sd, sd);
    m_gaussian->SetRadiusFactors(r, r, r);

    m_mc = vtkMarchingCubes::New();
    m_mc->SetInput(m_gaussian->GetOutput());
    m_mc->SetValue(0, props.inqMcThreshold());
 
    vtkPolyDataConnectivityFilter *connect = vtkPolyDataConnectivityFilter::New();
    connect->SetInput(m_mc->GetOutput());
    connect->SetExtractionModeToLargestRegion();

    vtkDecimatePro *decim = vtkDecimatePro::New();
    decim->SetInput(connect->GetOutput());
    decim->PreserveTopologyOn();
    decim->SplittingOn();

    m_smooth->SetInput(decim->GetOutput());
    m_smooth->SetNumberOfIterations(props.inqIterations());
    m_smooth->SetRelaxationFactor(props.inqRelaxationFactor());

    m_normals->SetInput(m_smooth->GetOutput());
    m_normals->SplittingOn();
    m_normals->ConsistencyOn();
    m_normals->NonManifoldTraversalOff();
    m_normals->SetFeatureAngle(props.inqFeatureAngle());

//     m_sorter->SetInput(m_smooth->GetOutput());
//     m_sorter->SortScalarsOn();
//     m_sorter->SetCamera(ren->GetActiveCamera());
//     m_sorter->SetDirectionToBackToFront();

    connect->Delete();
    decim->Delete();
  }

  vtkImageData *GetThreshOutput() { return m_thresh->GetOutput(); }
  
  virtual void update(const BriCon* bc) 
  { 
    m_imageLut->SetRange(bc->inqMin(), bc->inqMax());
    Render();
  }

  void setProperties(const VTKProperties& props)
  {
    float sd(props.inqStdDev());
    float r(props.inqRadius());

    m_gaussian->SetStandardDeviations(sd, sd, sd);
    m_gaussian->SetRadiusFactors(r, r, r);
    m_mc->SetValue(0, props.inqMcThreshold());
    m_normals->SetFeatureAngle(props.inqFeatureAngle());
    m_smooth->SetNumberOfIterations(props.inqIterations());
    m_smooth->SetRelaxationFactor(props.inqRelaxationFactor());
    m_thresh->ThresholdByUpper(props.inqLowerThreshold());
  }

private:
  vtkImageThreshold *m_thresh;
  vtkImageGaussianSmooth *m_gaussian;
  vtkMarchingCubes *m_mc;
};

struct VTKWidget::Implementation
{
  Implementation(QWidget* parent, OverlayList::Handle ol) :   
    m_brainActor(vtkActor::New()),
    m_sorter(vtkDepthSortPolyData::New()),
    m_surfMapper(vtkPolyDataMapper::New()),
    m_clipper(vtkClipPolyData::New()),
    m_erode(vtkImageDilateErode3D::New()),
    m_meshOptions(parent),
    m_ol(ol)
  {
    for(int p = 0; p < 6; ++p)
      m_cutActor[p] = vtkActor::New();
  }

  ~Implementation()
  {
    TRACKER("VTKWidget::Implementation::~Implementation");

    MESSAGE("destroying implementation elements");

    std::for_each(m_actors.begin(), m_actors.end(), 
		  DeleteVTKObject<vtkActor>());
    std::for_each(m_gaussians.begin(), m_gaussians.end(), 
		  DeleteVTKObject<vtkImageGaussianSmooth>());

    m_sorter->Delete();
    m_surfMapper->Delete();
    m_clipper->Delete();
    m_erode->Delete();
    for(int p = 0; p < 6; ++p)
      m_cutActor[p]->Delete();
  }

  vtkActor *m_brainActor;
  vtkActor *m_cutActor[6];
  vtkDepthSortPolyData *m_sorter;
  vtkPolyDataMapper *m_surfMapper;
  vtkClipPolyData *m_clipper;
  vtkImageDilateErode3D *m_erode;
  vtkRenderer *m_renderer;

  std::vector<vtkActor *> m_actors;
  std::vector<ImagePipelineObject::Handle> m_pipelineObjects;
  std::vector<VTKMeshSurface::Handle> m_meshes;
  std::vector<MetaImage::Handle> m_metaImages;
  std::vector<vtkImageGaussianSmooth *> m_gaussians;

  MainImage::Handle m_mainImage;
  VTKProperties m_props;
  MeshOptionsDialog m_meshOptions;
  OverlayList::Handle m_ol;
};

VTKWidget::VTKWidget(QWidget *parent, 
		     ImageGroup::Handle i,
		     OverlayList::Handle ol, 
		     Cursor::Handle c) :  
  ImageWidget(parent,i,ol,c), m_impl(new Implementation(this, ol)), m_image(i)
{
  TRACKER("VTKWidget::VTKWidget");
  
  QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );

  m_vtkwidget = new QVTKWidget( this, "vtkWidget" );
  setCentralWidget(m_vtkwidget);

  QToolBar *tb = new QToolBar(this);
  VTKToolbar *vt = new VTKToolbar(tb, m_impl->m_props);

  addDockWindow(tb, tr("VTK Rendering Tools"), Top, FALSE);  
  m_toolbar->hide();
  m_modebar->hide();

  connect(vt->m_printButton,   SIGNAL(clicked()), SLOT(print()));
  connect(vt->m_optionsButton, SIGNAL(clicked()), SLOT(options()));
  connect(vt->m_addMeshButton, SIGNAL(clicked()), SLOT(addMesh()));
  connect(vt->m_meshOptionsButton, SIGNAL(clicked()), SLOT(meshOptions()));

  vtkRenderer *ren = vtkRenderer::New();
  m_impl->m_renderer=ren;
  ren->SetBackground(0.4, 0.4, 0.4);
  m_vtkwidget->GetRenderWindow()->AddRenderer(ren);

  //  LookUpTableFactory *lut = new LookUpTableFactory();
    
  m_impl->m_mainImage = MainImage::Handle(new MainImage(ol->getMainMetaImage(), m_impl->m_props, ren));
  m_impl->m_mainImage->update(ol->getMainMetaImage()->getDs()->inqBriCon().get());

  vtkImageBlend *blend = vtkImageBlend::New();
  blend->SetBlendModeToNormal();
  blend->AddInput(m_impl->m_mainImage->GetOutputRGBA());

  unsigned int count(1);
  MetaImageListIt it = ol->begin();
  ++it;
  for(; it != ol->end(); ++it)
    {
      // Process the overlays for possible rendering/blending along
      // with the main image surface.
      MetaImage::Handle mi = (*it);

      if(mi->inqVisibility()) {

	Image::Handle im(mi->getImage());
	ImageInfo::Handle info(mi->getInfo());
	BriCon::Handle bc(mi->getDs()->inqBriCon());

	vtkImageData *layer = vtkImageData::New();
	layer->SetDimensions(info->inqX(), info->inqY(), info->inqZ());
	layer->SetSpacing(info->inqXDim(), info->inqYDim(), info->inqZDim());
	layer->SetScalarTypeToFloat();
	layer->SetNumberOfScalarComponents(1);
	layer->AllocateScalars();
	
// 	vtkPoints *points = vtkPoints::New();
// 	vtkFloatArray *tensors = vtkFloatArray::New();
// 	tensors->SetNumberOfComponents(9);

// 	unsigned int offset(0);

	if(info->isMaskImage()) {
	  
	  MaskImage::Handle maskim(new MaskImage(mi, ren));

	  vtkActor *layerActor = vtkActor::New();
	  layerActor->SetMapper(maskim->getLayerMapper());
	  layerActor->GetProperty()->SetOpacity(0.4);

	  ren->AddViewProp(layerActor);
	  
	  m_impl->m_actors.push_back(layerActor);
	  m_impl->m_pipelineObjects.push_back(maskim);
	  m_impl->m_metaImages.push_back(mi);
	}
	if(info->isStatImage()) {
	  // Accumulate this layer into composite "blend"
	  // image for rendering onto cut surfaces
	  
	  StatsImage::Handle si(new StatsImage(mi, ren));

	  MESSAGE("Blending in stat image");
	  blend->AddInput(si->GetOutputRGBA());
	  blend->SetOpacity(count, 1.0); count++;
	  
	  vtkActor *layerActor = vtkActor::New();
	  layerActor->SetMapper(si->getLayerMapper());
	  layerActor->GetProperty()->SetOpacity(0.4);

	  ren->AddViewProp(layerActor);
	  
	  m_impl->m_actors.push_back(layerActor);
	  m_impl->m_pipelineObjects.push_back(si);
	  m_impl->m_metaImages.push_back(mi);
	  
	} else if(info->isDtiImage()) {
// 	  // Create glyphs for each tensor
// 	  // in the data set.
// 	  MESSAGE("Showing tensor image");

// 	  //	  float *ptr = (float *)layer->GetScalarPointer();
// 	  for(int z = 0; z < info->inqZ(); ++z)
// 	    for(int y = 0; y < info->inqY(); ++y)
// 	      for(int x = 0; x < info->inqX(); ++x)
// 		{
// 		  float vx, vy, vz, mmx, mmy, mmz;

// 		  vx = im->getVolume(0)->value(x, y, z);
// 		  vy = im->getVolume(1)->value(x, y, z);
// 		  vz = im->getVolume(2)->value(x, y, z);

// 		  if((vx != 0) && (vy != 0) && (vz != 0)) {
// 		    FslGetMMCoord(info->inqStdMat(), 
// 				  x, y, z, &mmx, &mmy, &mmz);
		  
// 		    float tensor[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
// 		    for(unsigned short j = 0; j < 9; ++j)
// 		      tensor[j] = im->getVolume(j)->value(x, y, z);

// 		    points->InsertPoint(offset, mmx, mmy, mmz);
// 		    tensors->InsertTuple(offset, tensor);
// 		    ++offset;
// 		  }
// 		}

// 	  vtkSphereSource *sphere = vtkSphereSource::New();
// 	  sphere->SetCenter(0.0, 0.0, 0.0);
// 	  sphere->SetRadius(1.0);
// 	  sphere->SetThetaResolution(20);
// 	  sphere->SetPhiResolution(20);

// 	  vtkPolyData *tensorData = vtkPolyData::New();
// 	  tensorData->SetPoints(points); points->Delete();
// 	  tensorData->GetPointData()->SetTensors(tensors); tensors->Delete();

// 	  vtkTensorGlyph *glyph = vtkTensorGlyph::New();
// 	  glyph->ExtractEigenvaluesOff();
// // 	  glyph->ThreeGlyphsOn();
// 	  glyph->SetInput(tensorData);
// 	  glyph->SetSource(sphere->GetOutput());

// 	  vtkPolyDataNormals *normals = vtkPolyDataNormals::New();
// 	  normals->SetInput(glyph->GetOutput());
// 	  vtkPolyDataMapper *ellipseMapper = vtkPolyDataMapper::New();
// 	  ellipseMapper->SetInput(normals->GetOutput());
	  
// 	  LookUpTable::Handle lh(mi->getDs()->inqLookUpTable());
// 	  vtkLookupTable *thisLut;
// 	  if(lh)
// 	    thisLut = LookUpTableFactory::convert(lh);
// 	  else
// 	    thisLut = LookUpTableFactory::convert(LookUpTable::greyScale());
	       

// 	  vtkActor *ellipseActor = vtkActor::New();
// 	  ellipseActor->SetMapper(ellipseMapper);
// 	  ellipseActor->GetProperty()->SetOpacity(0.8);
// 	  ellipseActor->GetProperty()->SetColor(thisLut->GetTableValue(16000));

// 	  ren->AddViewProp(ellipseActor);
// 	  m_impl->m_actors.push_back(ellipseActor);
// 	  m_impl->m_metaImages.push_back(mi);
	}
      }
    }

  vtkBoxWidget *box = vtkBoxWidget::New();
  box->SetKeyPressActivationValue('b');
  box->SetInteractor(ren->GetRenderWindow()->GetInteractor());
  box->PlaceWidget(10,200,10,200,10,200);
  vtkPlanes *planes = vtkPlanes::New();
  box->GetPlanes(planes);
  vtkPolyData *pd = vtkPolyData::New();
  box->GetPolyData(pd);

  m_impl->m_clipper->SetInput(m_impl->m_mainImage->GetSurface());
  m_impl->m_clipper->SetClipFunction(planes);

  m_impl->m_sorter->SetInput(m_impl->m_clipper->GetOutput());
  m_impl->m_sorter->SetCamera(ren->GetActiveCamera());
  m_impl->m_sorter->SortScalarsOn();
  m_impl->m_sorter->SetDirectionToBackToFront();

  vtkPolyData *ds = m_impl->m_mainImage->GetSurface();
  if(m_impl->m_props.inqClipping())
    ds = m_impl->m_clipper->GetOutput();
  m_impl->m_sorter->SetInput(ds);
  if(m_impl->m_props.inqOpacity() != 1.0)
    m_impl->m_surfMapper->SetInput(m_impl->m_sorter->GetOutput());
  else
    m_impl->m_surfMapper->SetInput(ds);
  m_impl->m_surfMapper->ScalarVisibilityOff();

//   vtkTransform *boxXForm = vtkTransform::New();
//   box->GetTransform(boxXForm);
//   xform->SetTransform(boxXForm);

  vtkTransform *boxXForm = vtkTransform::New();
  box->GetTransform(boxXForm);

  ClippingBoxCallback *callback = ClippingBoxCallback::New();
  callback->SetPlanes(planes);
  callback->SetTransform(boxXForm);
  box->AddObserver(vtkCommand::EndInteractionEvent, callback);

  InteractorCallback *icallback = InteractorCallback::New();
  icallback->SetClipper(m_impl->m_clipper);
  icallback->SetMapper(m_impl->m_surfMapper);
  icallback->SetSorter(m_impl->m_sorter);
  icallback->SetProperties(&(m_impl->m_props));
  m_vtkwidget->GetInteractor()->AddObserver(vtkCommand::StartInteractionEvent, icallback);
  m_vtkwidget->GetInteractor()->AddObserver(vtkCommand::EndInteractionEvent, icallback);

  m_impl->m_brainActor->SetMapper(m_impl->m_surfMapper);
  m_impl->m_brainActor->GetProperty()->SetInterpolation(m_impl->m_props.inqInterpMode());
  m_impl->m_brainActor->GetProperty()->SetAmbient(m_impl->m_props.inqAmbient());
  m_impl->m_brainActor->GetProperty()->SetDiffuse(m_impl->m_props.inqDiffuse());
  m_impl->m_brainActor->GetProperty()->SetOpacity(m_impl->m_props.inqOpacity());
  m_impl->m_brainActor->GetProperty()->SetSpecular(m_impl->m_props.inqSpecular());
  m_impl->m_brainActor->GetProperty()->SetSpecularPower(m_impl->m_props.inqSpecularPower());
  float cr, cg, cb;
  m_impl->m_props.inqColor(cr, cg, cb);
  m_impl->m_brainActor->GetProperty()->SetColor(cr, cg, cb);
  m_impl->m_actors.push_back(m_impl->m_brainActor);
  m_impl->m_pipelineObjects.push_back(m_impl->m_mainImage);
  m_impl->m_metaImages.push_back(ol->getMainMetaImage());

  // Add Actor to renderer
  ren->AddViewProp(m_impl->m_brainActor);

  unsigned short verts[][3] = { {1, 5, 0}, 
				{0, 4, 3}, {1, 2, 0}, {2, 6, 3},
				{5, 6, 4}, {1, 2, 5}
  };

  MESSAGE("Adding clipped planes");
  for(unsigned short p = 0; p < 6; ++p)
    {
      vtkPlaneSource *ps = vtkPlaneSource::New();

      ps->SetOrigin(pd->GetPoints()->GetPoint(verts[p][0]));
      ps->SetPoint1(pd->GetPoints()->GetPoint(verts[p][1]));
      ps->SetPoint2(pd->GetPoints()->GetPoint(verts[p][2]));

      ps->SetXResolution(100);
      ps->SetYResolution(100);
  
      vtkTransformPolyDataFilter *xform = vtkTransformPolyDataFilter::New();
      xform->SetInput(ps->GetOutput());
      xform->SetTransform(boxXForm);

      vtkImplicitDataSet *clipFun = vtkImplicitDataSet::New();
      clipFun->SetDataSet(m_impl->m_mainImage->GetThreshOutput());
      clipFun->SetOutValue(0.0);

      vtkClipPolyData *clipPoly = vtkClipPolyData::New();
      clipPoly->SetClipFunction(clipFun);
      clipPoly->SetInput(xform->GetOutput());

      vtkProbeFilter *probe = vtkProbeFilter::New();
      probe->SetInput(clipPoly->GetOutput());
      probe->SetSource(blend->GetOutput());

//       vtkPolyDataWriter *pdw = vtkPolyDataWriter::New();
//       pdw->SetInput(probe->GetPolyDataOutput());
//       pdw->SetFileName("/tmp/Goop");
//       pdw->Write();

      vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
      mapper->SetInput(probe->GetPolyDataOutput());

      m_impl->m_cutActor[p]->SetMapper(mapper);
      m_impl->m_cutActor[p]->GetProperty()->SetAmbient(1.0);
      m_impl->m_cutActor[p]->GetProperty()->SetDiffuse(0.2);
      m_impl->m_cutActor[p]->GetProperty()->SetInterpolationToGouraud();
      m_impl->m_cutActor[p]->GetProperty()->SetOpacity(1);

      ren->AddViewProp(m_impl->m_cutActor[p]);
      m_impl->m_cutActor[p]->VisibilityOff();
   }

  vtkImagePlaneWidget *planex = vtkImagePlaneWidget::New();
  planex->SetInput(m_impl->m_mainImage->GetThreshOutput());
  planex->SetPlaneOrientationToXAxes();
  planex->GetColorMap()->GetLookupTable()->SetAlpha(0.6);
  planex->SetKeyPressActivationValue('x');
  planex->GetTexturePlaneProperty()->SetOpacity(1);
  planex->DisplayTextOn();
  vtkImagePlaneWidget *planey = vtkImagePlaneWidget::New();
  planey->SetInput(m_impl->m_mainImage->GetThreshOutput());
  planey->SetPlaneOrientationToYAxes();
  planey->GetColorMap()->GetLookupTable()->SetAlpha(0.6);
  planey->SetKeyPressActivationValue('y');
  planey->GetTexturePlaneProperty()->SetOpacity(1);
  planey->DisplayTextOn();
  vtkImagePlaneWidget *planez = vtkImagePlaneWidget::New();
  planez->SetInput(m_impl->m_mainImage->GetThreshOutput());
  planez->SetPlaneOrientationToZAxes();
  planez->GetColorMap()->GetLookupTable()->SetAlpha(0.6);
  planez->SetKeyPressActivationValue('z');
  planez->GetTexturePlaneProperty()->SetOpacity(1);
  planez->DisplayTextOn();

  planey->SetLookupTable(planex->GetLookupTable());
  planez->SetLookupTable(planex->GetLookupTable());

  planex->SetInteractor(ren->GetRenderWindow()->GetInteractor());
  planey->SetInteractor(ren->GetRenderWindow()->GetInteractor());
  planez->SetInteractor(ren->GetRenderWindow()->GetInteractor());

  MESSAGE("Rendering");
  
  // Reset camera
  ren->ResetCamera();
  m_vtkwidget->GetRenderWindow()->Render();
  QApplication::restoreOverrideCursor();

  m_impl->m_ol->attach(this);
  m_impl->m_props.attach(this);
}

VTKWidget::~VTKWidget()
{
  TRACKER("VTKWidget::~VTKWidget");
  m_impl->m_ol->detach(this);
  m_impl->m_props.detach(this);
}

void VTKWidget::update(const BriCon* bc)
{
  TRACKER("VTKWidget::update(const BriCon*)");

  m_impl->m_surfMapper->SetScalarRange(bc->inqMin(), bc->inqMax());  
}

void VTKWidget::addMesh()
{
  //FileDialog to find mesh file

  QString fn = QFileDialog::getOpenFileName(QDir::currentDirPath(),
					    "Mesh files (*.vtk)", 
					    this );
  if(!fn.isEmpty()) 
    {
      QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );
      ImageInfo::Handle info(m_impl->m_ol->getMainMetaImage()->getInfo());
      float xoff(info->inqXDim() * (info->inqX() - 1));
      VTKMeshSurface::Handle mesh = VTKMeshSurface::create(m_impl->m_renderer, fn, xoff);
      m_impl->m_meshes.push_back(mesh);
      m_impl->m_renderer->AddViewProp(mesh->getActor());
      m_impl->m_renderer->AddViewProp(mesh->getGlyphActor());
      m_impl->m_meshOptions.populateMeshList(m_impl->m_meshes);
      QApplication::restoreOverrideCursor();
    }
}

void VTKWidget::meshOptions()
{
  if(m_impl->m_meshes.size())
    {
      m_impl->m_meshOptions.populateMeshList(m_impl->m_meshes);
      m_impl->m_meshOptions.show();
    }
}

void VTKWidget::update(const VTKProperties* p)
{
  m_impl->m_mainImage->setProperties(*p);

  vtkPolyData *ds = m_impl->m_mainImage->GetSurface();

  if(m_impl->m_props.inqClipping()) {
    ds = m_impl->m_clipper->GetOutput();
    for(int p = 0; p < 6; ++p)
      m_impl->m_cutActor[p]->VisibilityOn();
  } else {
    for(int p = 0; p < 6; ++p)
      m_impl->m_cutActor[p]->VisibilityOff();
  }
  m_impl->m_sorter->SetInput(ds);

  if(m_impl->m_ol->getMainMetaImage()->inqTransparency() != 1.0)
    m_impl->m_surfMapper->SetInput(m_impl->m_sorter->GetOutput());
  else
    m_impl->m_surfMapper->SetInput(ds);

  QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );
  m_vtkwidget->GetRenderWindow()->Render();
  QApplication::restoreOverrideCursor();
}

void VTKWidget::update(const Cursor::Handle c)
{
  TRACKER("VTKWidget::update(const Cursor::Handle)");
}

void VTKWidget::update(const OverlayList *ol, OverlayListMsg msg)
{
  TRACKER("VTKWidget::update(const OverlayList *ol, OverlayListMsg msg)");
  
  vtkPolyData *ds = m_impl->m_mainImage->GetSurface();
  if(m_impl->m_props.inqClipping())
    ds = m_impl->m_clipper->GetOutput();
  m_impl->m_sorter->SetInput(ds); // AddInputConnection()
  if(m_impl->m_ol->getMainMetaImage()->inqTransparency() != 1.0)
    m_impl->m_surfMapper->SetInput(m_impl->m_sorter->GetOutput());
  else
    m_impl->m_surfMapper->SetInput(ds);

  unsigned int count(0);

  for(std::vector<vtkActor *>::iterator it = m_impl->m_actors.begin(); 
      it != m_impl->m_actors.end(); ++it, ++count)
    {
      vtkActor *thisActor = m_impl->m_actors.at(count);
      ImagePipelineObject::Handle po = m_impl->m_pipelineObjects.at(count);
      MetaImage::Handle mi = m_impl->m_metaImages.at(count);

      thisActor->SetVisibility(mi->inqVisibility());
      thisActor->GetProperty()->SetOpacity(mi->inqTransparency());
      //thisActor->SetMapper(po->getLayerMapper());
      //po->setSorting(mi->inqTransparency() != 1.0);

      LookUpTable::Handle lh(mi->getDs()->inqLookUpTable());
      vtkLookupTable *lut;
      if(lh)
	lut = LookUpTableFactory::convert(lh);
      else
	lut = LookUpTableFactory::convert(LookUpTable::greyScale());
      po->setImageLut(lut);
      po->setSurfaceLut(lut);
    }
  
  QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );
  m_vtkwidget->GetRenderWindow()->Render();
  QApplication::restoreOverrideCursor();
}

void VTKWidget::print() 
{
  QString fn = QFileDialog::getSaveFileName("screenshot.tiff", 
					    "TIFF files (*.tif, *.tiff)", this,
					    "Screenshot dialog",
					    "Select a filename for saving");

  if(!fn.isNull()) {
    vtkWindowToImageFilter *w2i = vtkWindowToImageFilter::New();
    vtkTIFFWriter *writer = vtkTIFFWriter::New();
    
    w2i->SetInput(m_vtkwidget->GetRenderWindow());
    w2i->Update();
    writer->SetInput(w2i->GetOutput());
    writer->SetFileName((const char *)fn);
    QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );
    m_vtkwidget->GetRenderWindow()->Render();
    QApplication::restoreOverrideCursor();
    writer->Write();

    writer->Delete();
    w2i->Delete();
  }
}

struct SetGaussianParams
{
  SetGaussianParams(float sd, float  r):
    m_sd(sd), m_r(r)
  {}

  void operator()(vtkImageGaussianSmooth *g) const
  {
    g->SetStandardDeviations(m_sd, m_sd, m_sd);
    g->SetRadiusFactors(m_r, m_r, m_r);
  }

  float m_sd, m_r;
};

void VTKWidget::options()
{
  VTKPropertyDialog optionsDialog(this, m_impl->m_props);
  
  if(optionsDialog.exec() == QDialog::Accepted)
    {
      m_impl->m_props = optionsDialog.getProperties();
      float r(m_impl->m_props.inqRadius());
      float sd(m_impl->m_props.inqStdDev());
      m_impl->m_mainImage->setProperties(m_impl->m_props);

      std::for_each(m_impl->m_gaussians.begin(), m_impl->m_gaussians.end(),
		    SetGaussianParams(sd, r));

      m_impl->m_brainActor->GetProperty()->SetInterpolation(m_impl->m_props.inqInterpMode());
      m_impl->m_brainActor->GetProperty()->SetAmbient(m_impl->m_props.inqAmbient());
      m_impl->m_brainActor->GetProperty()->SetDiffuse(m_impl->m_props.inqDiffuse());
      m_impl->m_brainActor->GetProperty()->SetOpacity(m_impl->m_props.inqOpacity());
      m_impl->m_brainActor->GetProperty()->SetSpecular(m_impl->m_props.inqSpecular());
      m_impl->m_brainActor->GetProperty()->SetSpecularPower(m_impl->m_props.inqSpecularPower());
      float cr, cg, cb;
      m_impl->m_props.inqColor(cr, cg, cb);
      m_impl->m_brainActor->GetProperty()->SetColor(cr, cg, cb);
      if(m_impl->m_props.inqOpacity() != 1.0)
	m_impl->m_surfMapper->SetInput(m_impl->m_sorter->GetOutput());
      else if(m_impl->m_props.inqClipping())
	m_impl->m_surfMapper->SetInput(m_impl->m_clipper->GetOutput());
      else
	m_impl->m_surfMapper->SetInput(m_impl->m_mainImage->GetSurface());
	
      QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );
      m_vtkwidget->GetRenderWindow()->Render();
      QApplication::restoreOverrideCursor();
    }
}

