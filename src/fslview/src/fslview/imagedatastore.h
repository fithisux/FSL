/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(IMAGEDATASTORE_H)
#define IMAGEDATASTORE_H

#include <boost/shared_ptr.hpp>
#include "metaimage.h"
#include "imagedata.h"
#include "overlaylist.h"

class ImageDataStore
{
public:
  typedef boost::shared_ptr< ImageDataStore > Handle;
  static Handle create(OverlayList::Handle);
  bool isDtiLineOverlay();
  ImageData::Handle getDtiLineOverlay();
  void resetPos();        
  bool currentEmpty();
  void next();
  ImageData::Handle current();
  virtual ~ImageDataStore();
 
private:
  ImageDataStore(OverlayList::Handle);

  struct Implementation;  
  const std::auto_ptr<Implementation> m_impl;
};


#endif
