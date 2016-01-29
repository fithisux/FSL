/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */


// Image.h: interface for the Image class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(IMAGE_H)
#define IMAGE_H

#include <boost/shared_ptr.hpp>
#include <string>
#include <vector>
#include <list>
#include <map>
#include "fslio/fslio.h"
#include "volume.h"
#include "timeseries.h"
#include "imageinfo.h"

#include <stdexcept>

//! @brief Container class for image volumes.
//!
//! @author Dave Flitney
//!
//! Provides interface for getting access to image data either by Volume or
//! TimeSeries.
class Image
{
  friend class testImage;

public:
  typedef boost::shared_ptr< Image > Handle;
  typedef std::map<long,TimeSeries::Handle> TimeSeriesMap;
  class Exception;

  virtual ~Image();

  static Image::Handle load(const std::string& filename, bool calc=true);
  bool valid() const;

  bool save(const std::string& filename);

  Handle cloneStructure();
  Handle clone3dStructure();
  Volume::Handle getVolume(short n, bool cache=true) const;
  TimeSeries::Handle getTimeSeries(short x,short y,short z) const;
  void clearCache();
  const ImageInfo::Handle getInfo() const;
  void setAvw(FSLIO* avw){m_avw = avw;}
  const FSLIO* getAvw()const {return m_avw;}
  
  //  FslNiftiExtension * getExtension(int n) { return FslGetExtension(m_avw, n); }

private:
  Image(ImageInfo::Handle);
  Image(const std::string& filename);

  FSLIO* m_avw;
  ImageInfo::Handle m_imageInfo;
  
  mutable std::vector<Volume::Handle> m_cachedVolumes;
  mutable TimeSeriesMap               m_cachedTimeSeries;
  Volume::Handle blankDraw();
};

class Image::Exception: public std::runtime_error
{
public:
  Exception(const std::string& s): std::runtime_error(s) {}
};

inline const ImageInfo::Handle Image::getInfo()     const { return m_imageInfo; }
inline bool isValidImage(const Image::Handle img)
{
  if(!img.get()){return false;}
  else           {return true;}
}
 
#endif
