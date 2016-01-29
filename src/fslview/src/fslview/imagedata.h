/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(IMAGEDATA_H)
#define IMAGEDATA_H

#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>
#include "metaimage.h"
#include "storage/image.h"
#include "imagedisplaysetting.h"

typedef boost::shared_array< ColorRGBA > ColorRGBAHandle;

class ImageData
{
public:
  typedef boost::shared_ptr< ImageData > Handle;
  static Handle create(MetaImage::Handle,ColorRGBAHandle);
  MetaImage::Handle getMetaImage();
  ColorRGBAHandle   getBuffer();
  void              setBuffer(ColorRGBAHandle);
  Image::Handle     getImage();   
  ImageInfo::Handle getInfo();
  ImageDisplaySetting::Handle getDs();
  bool inqVisibility() const;
  int  inqDtiDisplay() const;
  float inqTransparency() const;
  virtual ~ImageData();
 
private:
  ImageData(MetaImage::Handle, ColorRGBAHandle);

  struct Implementation;  
  const std::auto_ptr<Implementation> m_impl;
};


#endif
