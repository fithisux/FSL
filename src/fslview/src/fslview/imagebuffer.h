/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(IMAGEBUFFER_H)
#define IMAGEBUFFER_H

#include "metaimage.h"
#include "imagedata.h"


//! @brief Functions intended for rendering slices into ColorRGBAHandle
//! objects
//!
//! @author James Saunders
class ImageBuffer
{
public:

  static ColorRGBAHandle axialBuffer(MetaImage::Handle, int slice, int vol);
  static ColorRGBAHandle coronalBuffer(MetaImage::Handle, int slice, int vol);
  static ColorRGBAHandle sagittalBuffer(MetaImage::Handle, int slice, int vol);

  static ColorRGBAHandle axialDtiBuffer(MetaImage::Handle, int slice);
  static ColorRGBAHandle coronalDtiBuffer(MetaImage::Handle, int slice);
  static ColorRGBAHandle sagittalDtiBuffer(MetaImage::Handle, int slice);

  static void blendBuffers(ColorRGBAHandle dest, ColorRGBAHandle source,
                           float trans, bool bottomLayer, 
                           unsigned int length);

  static void setToZero(ColorRGBAHandle,unsigned int length);
  static void reorderBytes(ColorRGBAHandle, unsigned int length);

  static inline unsigned char clamp(float c)
  {return (unsigned char)( std::max(0.0, std::min(255.0, c * 255.0)) );}

};


#endif
