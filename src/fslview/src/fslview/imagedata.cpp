/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "imagedata.h"

struct ImageData::Implementation
{
  Implementation(MetaImage::Handle mi, ColorRGBAHandle rgba):
    m_metaImage(mi),m_rgba(rgba){};
  MetaImage::Handle m_metaImage;
  ColorRGBAHandle   m_rgba;
};

ImageData::ImageData(MetaImage::Handle mi, ColorRGBAHandle rgba):
  m_impl(new Implementation(mi,rgba)){}

ImageData::~ImageData(){}

ImageData::Handle ImageData::create(MetaImage::Handle mi, 
                                    ColorRGBAHandle rgba)
{
  Handle dst(new ImageData(mi,rgba));
  return dst;
}  

MetaImage::Handle ImageData::getMetaImage()
{
  return m_impl->m_metaImage;
}

ColorRGBAHandle ImageData::getBuffer()
{
  return m_impl->m_rgba;
}

bool ImageData::inqVisibility()const
{
  return m_impl->m_metaImage->inqVisibility();
}

int   ImageData::inqDtiDisplay()const
{
  return m_impl->m_metaImage->getDs()->inqDtiDisplay();
}

float ImageData::inqTransparency() const
{
  return m_impl->m_metaImage->inqTransparency();
}

ImageInfo::Handle ImageData::getInfo()
{
  return m_impl->m_metaImage->getInfo();
}

Image::Handle ImageData::getImage()
{
  return m_impl->m_metaImage->getImage();
}

ImageDisplaySetting::Handle ImageData::getDs()
{
  return m_impl->m_metaImage->getDs();
}

void ImageData::setBuffer(ColorRGBAHandle rgba)
{
  m_impl->m_rgba = rgba; 
}
