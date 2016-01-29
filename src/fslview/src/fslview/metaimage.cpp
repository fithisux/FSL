/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "metaimage.h"

struct MetaImage::Implementation
{
  Implementation(Image::Handle i,ImageDisplaySetting::Handle d):
    m_image(i),m_ds(d){};
  Image::Handle m_image;
  ImageDisplaySetting::Handle m_ds;
};

MetaImage::MetaImage(Image::Handle i, ImageDisplaySetting::Handle d):
  m_impl(new Implementation(i,d)){}

MetaImage::~MetaImage(){}

MetaImage::Handle MetaImage::create(Image::Handle i, 
                                    ImageDisplaySetting::Handle d)
{
  Handle dst(new MetaImage(i,d));
  return dst;
}  

Image::Handle MetaImage::getImage()const
{
  return m_impl->m_image;
}

ImageDisplaySetting::Handle MetaImage::getDs()
{
  return m_impl->m_ds;
}

ImageInfo::Handle MetaImage::getInfo()
{
  return m_impl->m_image->getInfo();
}

MetaImage::Handle MetaImage::clone()
{
  return Handle(new MetaImage(getImage(),getDs()->clone()));
}

short MetaImage::inqX()const{return m_impl->m_image->getVolume(0)->inqX();}
short MetaImage::inqY()const{return m_impl->m_image->getVolume(0)->inqY();}
short MetaImage::inqZ()const{return m_impl->m_image->getVolume(0)->inqZ();}
  
bool  MetaImage::inqVisibility()const  { return m_impl->m_ds->inqVisibility();}
bool  MetaImage::inqReadOnly()  const  { return m_impl->m_image->getInfo()->inqReadOnly();}
float MetaImage::inqTransparency()   const { return m_impl->m_ds->inqTransparency();}
std::string MetaImage::inqImageName()const { return m_impl->m_image->getInfo()->inqImageName();}

void MetaImage::setVisibility(bool state){m_impl->m_ds->setVisibility(state);}
void MetaImage::setReadOnly(bool state)  {m_impl->m_image->getInfo()->setReadOnly(state);}
void MetaImage::setTransparency(float trans){m_impl->m_ds->setTransparency(trans);}

void MetaImage::setImageName(std::string imageName)
{
  m_impl->m_image->getInfo()->setImageName(imageName);
}
