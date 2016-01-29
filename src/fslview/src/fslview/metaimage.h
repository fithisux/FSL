/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(METAIMAGE_H)
#define METAIMAGE_H

#include <boost/shared_ptr.hpp>
#include "storage/image.h"
#include "imagedisplaysetting.h"

/*
 * @author Dave Flitney
 *
 * @date   Dec 2002
 *
 * @brief  MetaImage objects associate an Image with an ImageDisplaySetting.
 * 
 * Use MetaImage to record the display options for each image.
 */ 
class MetaImage
{
public:
  typedef boost::shared_ptr< MetaImage > Handle;
  static Handle create(Image::Handle,ImageDisplaySetting::Handle);
  Image::Handle getImage()const;
  ImageDisplaySetting::Handle getDs();
  MetaImage::Handle clone();
  virtual ~MetaImage();
  short inqX()const;
  short inqY()const;
  short inqZ()const;
  bool  inqVisibility()const; 
  bool  inqReadOnly()const;
  float inqTransparency()const;
  std::string inqImageName()const;  
  void  setVisibility(bool); 
  void  setReadOnly(bool);
  void  setTransparency(float);
  void setImageName(std::string);

  ImageInfo::Handle getInfo();

private:
  MetaImage(Image::Handle i, ImageDisplaySetting::Handle d);

  struct Implementation;  
  const std::auto_ptr<Implementation> m_impl;
};

//typedef std::pair<ImageDisplaySetting::Handle,Image::Handle> MetaImage;
typedef std::list<MetaImage::Handle> MetaImageList;
typedef MetaImageList::iterator MetaImageListIt;














#endif
