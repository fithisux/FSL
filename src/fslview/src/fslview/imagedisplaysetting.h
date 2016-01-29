/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(IMAGEDISPLAYSETTING_H)
#define IMAGEDISPLAYSETTING_H

#if defined(WIN32)
#pragma warning(disable:4786)
#endif

#include <boost/shared_ptr.hpp>
#include <memory>
#include "storage/image.h"
#include "lookuptable.h"
#include "bricon.h"
  
typedef enum {None, Lines, RGB, LinesRGB} DtiDisplay;

class ImageDisplaySetting  
{
public:
  typedef boost::shared_ptr< ImageDisplaySetting > Handle;
  static Handle create(Image::Handle image, LookUpTable::Handle lut,float trans = 0.5, bool visible = true);
  
  void  setTransparency(float trans);
  float inqTransparency() const;

  void setVisibility(bool visible); 
  bool inqVisibility() const;

  void                setLookUpTable(LookUpTable::Handle);
  LookUpTable::Handle inqLookUpTable() const;
  void                setSecondaryLookUpTable(LookUpTable::Handle);
  LookUpTable::Handle inqSecondaryLookUpTable() const;
  void setUseSecondaryLookUpTable(bool);
  bool inqUseSecondaryLookUpTable() const;

  void setDtiDisplay(DtiDisplay);
  int  inqDtiDisplay() const;

  void setTransMod(bool);
  bool inqTransMod() const;

  void  setModTransparency(float);
  float inqModTransparency() const;

  void          setCurrentVolume(unsigned int);
  unsigned int  inqCurrentVolume() const;

  void          setModImage(Image::Handle);
  Image::Handle inqModImage() const;

  BriCon::Handle inqBriCon();

  Handle clone();

  virtual ~ImageDisplaySetting();

private:
  ImageDisplaySetting(Image::Handle image, LookUpTable::Handle lut,
                      float trans, bool visible);
  ImageDisplaySetting(BriCon::Handle bricon, 
		      LookUpTable::Handle lut, LookUpTable::Handle slut,
                      float trans, bool visible,
                      int dti, Image::Handle mod,
                      bool transMod, float modTransVal, int vol=0);

  struct Implementation;  
  const std::auto_ptr<Implementation> m_impl;
};



#endif
