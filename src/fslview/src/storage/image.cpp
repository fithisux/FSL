/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include <algorithm>
#include "image.h"
#include <stdio.h>
#include "error.h"
#include <assert.h>

using namespace std;

Image::Image(const std::string& filename):
  m_avw(NULL)
{
  if(!(m_avw = FslOpen(filename.c_str(), "r")))
    throw Exception("Failed to open file " + filename);
  m_imageInfo = ImageInfo::init(m_avw,filename); 
  m_cachedVolumes.clear();
  m_cachedVolumes.resize(m_imageInfo->inqNumVolumes());
}

Image::Image(ImageInfo::Handle info):
  m_avw(NULL), m_imageInfo(info)
{
  m_cachedVolumes.resize(info->inqNumVolumes());
 
  int x = m_imageInfo->inqX();
  int y = m_imageInfo->inqY();
  int z = m_imageInfo->inqZ();
  int v = m_imageInfo->inqNumVolumes();

  m_imageInfo->setMin(0);
  m_imageInfo->setMax(4);

  for(int n = 0; n < v; ++n)
    {
      switch(m_imageInfo->inqDt()) 
      {
       case DT_UNSIGNED_CHAR: m_cachedVolumes[n] = VolumeUB::createBlank(x,y,z); break;
       case DT_SIGNED_SHORT:  m_cachedVolumes[n] = VolumeS::createBlank(x,y,z); break;
       case DT_SIGNED_INT:    m_cachedVolumes[n] = VolumeI::createBlank(x,y,z); break;
       case DT_FLOAT:         m_cachedVolumes[n] = VolumeF::createBlank(x,y,z); break;
       case DT_DOUBLE:        m_cachedVolumes[n] = VolumeD::createBlank(x,y,z); break;
       
       // Support for "New" NIFTI data types
       
       case DT_INT8:		  m_cachedVolumes[n] = VolumeB::createBlank(x,y,z); break;
       case DT_UINT16:		  m_cachedVolumes[n] = VolumeUS::createBlank(x,y,z); break;
       case DT_UINT32:		  m_cachedVolumes[n] = VolumeUI::createBlank(x,y,z); break;
       case DT_INT64:		  m_cachedVolumes[n] = VolumeI64::createBlank(x,y,z); break;
       case DT_UINT64:		  m_cachedVolumes[n] = VolumeUI64::createBlank(x,y,z); break;
       case DT_FLOAT128:	  m_cachedVolumes[n] = VolumeF128::createBlank(x,y,z); break;
      }
    }
}

Image::~Image()
{  
  if(m_avw != NULL){FslClose(m_avw);}
}

Image::Handle Image::load(const std::string& filename, bool calc)
{
	Image::Handle im(new Image(filename));
	//   if(im->getInfo()->inqLutName() == "")
	//     im->getVolume(0)->calculateMinMax();
	//   else if(im->getInfo()->inqMax() <= im->getInfo()->inqMin())
	//     im->getVolume(0)->calculateMinMax();

	float min(im->getInfo()->inqMin());
	float max(im->getInfo()->inqMax());

	if((min == max) && calc) {
		switch(im->getInfo()->inqDt()) 
		{
		case DT_INT8:
		case DT_UNSIGNED_CHAR:
		case DT_SIGNED_SHORT:
		case DT_UINT16:
		case DT_SIGNED_INT:
		case DT_UINT32:
		case DT_INT64:
		case DT_UINT64:
			if((im->getVolume(0)->minValue() >= -100) &&
					(im->getVolume(0)->maxValue() <=  100))
			{
				// This call causes volume to cache min/max values
				im->getVolume(0)->calculateMinMax();
				break;
			}
		case DT_FLOAT:
		case DT_DOUBLE:
		case DT_FLOAT128:
			if((min == 0) || (max == 0))
				im->getVolume(0)->calculateRobustMinMax();
			else
				im->getVolume(0)->calculateRobustMinMax(true);
			break;
		default:
			throw FileError(filename,"Data type not supported.");
		}
	}

	return im;
}

bool Image::save(const std::string& filename)
{
  FSLIO *avw = FslOpen(filename.c_str(), "w");

  if( !avw ) return false;

  m_imageInfo->saveAvwHeader(avw);

  int v = m_imageInfo->inqNumVolumes();
  for(int n = 0; n < v; ++n)
    {
      Volume::Handle v = getVolume(n);

      v->saveVolume(avw);
    }

  FslClose(avw);

  return true;
}

Image::Handle Image::cloneStructure()
{
  Image::Handle h = Image::Handle(new Image(m_imageInfo->clone()));
  h->getInfo()->setLutName("Red-Yellow");

  return h;
}

Image::Handle Image::clone3dStructure()
{
  ImageInfo::Handle i(m_imageInfo->clone());

  i->setNumVolumes(1);

  Image::Handle h = Image::Handle(new Image(i));
  h->getInfo()->setLutName("Red-Yellow");

  return h;
}

void Image::clearCache()
{
  vector<Volume::Handle>::size_type sz(m_cachedVolumes.size());
  m_cachedVolumes.clear();
  m_cachedVolumes.resize(sz);
}

Volume::Handle Image::getVolume(short n, bool cache) const
{
  n = std::max(0, std::min(int(n), m_imageInfo->inqNumVolumes() - 1));

  float min = m_imageInfo->inqMin();
  float max = m_imageInfo->inqMax();

  Volume::Handle v(m_cachedVolumes[n]);

  if(!v)
    {
      std::string name;
      
      switch(m_imageInfo->inqDt()) 
      {
      case DT_UNSIGNED_CHAR: 
    	  name="unsigned char"; 
    	  v = VolumeUB::getVolume(m_avw, min, max, n);
    	  break;
      case DT_SIGNED_SHORT:
    	  name="signed short"; 
    	  v = VolumeS::getVolume(m_avw, min, max, n);
    	  break;
      case DT_SIGNED_INT: 
    	  name="signed int"; 
    	  v = VolumeI::getVolume(m_avw, min, max, n);
    	  break;
      case DT_FLOAT:
    	  name="float";
    	  v = VolumeF::getVolume(m_avw, min, max, n);
    	  break;
      case DT_DOUBLE:
    	  name="double";
    	  v = VolumeD::getVolume(m_avw, min, max, n);
    	  break;
      case DT_INT8: 
    	  name="char"; 
    	  v = VolumeB::getVolume(m_avw, min, max, n);
    	  break;
      case DT_UINT16:
    	  name="unsigned short"; 
    	  v = VolumeUS::getVolume(m_avw, min, max, n);
    	  break;
      case DT_UINT32: 
    	  name="unsigned int"; 
    	  v = VolumeUI::getVolume(m_avw, min, max, n);
    	  break;
      case DT_INT64: 
    	  name="long long"; 
    	  v = VolumeI64::getVolume(m_avw, min, max, n);
    	  break;
      case DT_UINT64: 
    	  name="unsigned long long"; 
    	  v = VolumeUI64::getVolume(m_avw, min, max, n);
    	  break;
      case DT_FLOAT128:
    	  name="float";
    	  v = VolumeF128::getVolume(m_avw, min, max, n);
    	  break;
      default:
    	  throw FileError(m_imageInfo->inqFileName(),"Data type not supported.");
      }
    }
  if(cache) 
    m_cachedVolumes[n] = v;
  return v;
}

TimeSeries::Handle Image::getTimeSeries(short x, short y, short z) const
{
  short xDim,yDim,zDim;
  long  voxel;
  TimeSeriesMap::iterator i;
  
  xDim = m_imageInfo->inqX();
  yDim = m_imageInfo->inqY();
  zDim = m_imageInfo->inqZ();

  voxel = xDim*(y + z * yDim) + x;
  i =  m_cachedTimeSeries.find(voxel);

  if(m_cachedTimeSeries.end() == i)
  {     
	  switch(m_imageInfo->inqDt()) 
	  {
	  case DT_INT8:			 m_cachedTimeSeries.insert(TimeSeriesMap::value_type(voxel, TimeSeriesB::getTimeSeries(m_avw,x,y,z))); break;
	  case DT_UNSIGNED_CHAR: m_cachedTimeSeries.insert(TimeSeriesMap::value_type(voxel, TimeSeriesUB::getTimeSeries(m_avw,x,y,z))); break;
	  case DT_SIGNED_SHORT:  m_cachedTimeSeries.insert(TimeSeriesMap::value_type(voxel, TimeSeriesS::getTimeSeries(m_avw,x,y,z))); break;
	  case DT_UINT16:		 m_cachedTimeSeries.insert(TimeSeriesMap::value_type(voxel, TimeSeriesUS::getTimeSeries(m_avw,x,y,z))); break;
	  case DT_SIGNED_INT:    m_cachedTimeSeries.insert(TimeSeriesMap::value_type(voxel, TimeSeriesI::getTimeSeries(m_avw,x,y,z))); break;
	  case DT_UINT32:		 m_cachedTimeSeries.insert(TimeSeriesMap::value_type(voxel, TimeSeriesUI::getTimeSeries(m_avw,x,y,z))); break;
	  case DT_INT64:		 m_cachedTimeSeries.insert(TimeSeriesMap::value_type(voxel, TimeSeriesI64::getTimeSeries(m_avw,x,y,z))); break;
	  case DT_UINT64:		 m_cachedTimeSeries.insert(TimeSeriesMap::value_type(voxel, TimeSeriesUI64::getTimeSeries(m_avw,x,y,z))); break;
	  case DT_FLOAT:         m_cachedTimeSeries.insert(TimeSeriesMap::value_type(voxel, TimeSeriesF::getTimeSeries(m_avw,x,y,z))); break;
	  case DT_DOUBLE:        m_cachedTimeSeries.insert(TimeSeriesMap::value_type(voxel, TimeSeriesD::getTimeSeries(m_avw,x,y,z))); break;
	  case DT_FLOAT128:		 m_cachedTimeSeries.insert(TimeSeriesMap::value_type(voxel, TimeSeriesF128::getTimeSeries(m_avw,x,y,z))); break;
	  }
  }
 return m_cachedTimeSeries[voxel];
}

Volume::Handle Image::blankDraw()
{
  VolumeS::Handle v = VolumeS::create(m_imageInfo->inqX(),
                                      m_imageInfo->inqY(),
                                      m_imageInfo->inqZ(),
                                      new short[m_imageInfo->inqX()*
                                                m_imageInfo->inqY()*
                                                m_imageInfo->inqZ()]);
  v->setMin(0);
  v->setMax(255);
  

  unsigned int width  = v->inqX();
  unsigned int height = v->inqY();
  unsigned int depth  = v->inqZ();

  for(unsigned int z = 0; z< depth; ++z)    
  {
    for(unsigned int y = 0; y < height; ++y)
    {      
      for(unsigned int x = 0; x < width; ++x)
        {
          (*v)(x,y,z) = 0;
        }
    }
  }
  
  return v;
}
