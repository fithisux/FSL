/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if defined(WIN32)
#pragma warning(disable:4786)
#endif

//#define DEBUGGING

#include "lookuptable.h"
#include "imagedisplaysetting.h"

#include "tracker.h"

#include <iostream>
using namespace std;

struct ImageDisplaySetting::Implementation
{
  Implementation(LookUpTable::Handle lut,float trans, bool visible, Image::Handle image): 
    m_transparency(trans), 
    m_isVisible(visible), 
    m_lookUpTable(lut),
    m_dtiDisplay(None),
    m_transMod(false),
    m_useSlut(false),
    m_modTransparency(0.0f),
    m_currentVolume(0)
  {
    Volume::Handle v = image->getVolume(0);

    m_bricon = BriCon::create(v->inqMin(), v->inqMax());
  }

  Implementation(LookUpTable::Handle lut,
		 LookUpTable::Handle slut,
                 float trans, 
                 bool visible, 
                 BriCon::Handle bricon,
                 int dti,
                 Image::Handle mod,
                 bool transMod,
                 float modTransVal,
		 int vol=0): 
    m_transparency(trans), 
    m_isVisible(visible), 
    m_lookUpTable(lut),
    m_secondaryLookUpTable(slut),
    m_bricon(bricon),
    m_dtiDisplay(dti),
    m_modImage(mod),
    m_transMod(transMod),
    m_useSlut(false),
    m_modTransparency(modTransVal),
    m_currentVolume(vol)
  { }

  float               m_transparency;
  bool                m_isVisible;
  LookUpTable::Handle m_lookUpTable;
  LookUpTable::Handle m_secondaryLookUpTable;
  BriCon::Handle      m_bricon;  
  int                 m_dtiDisplay;
  Image::Handle       m_modImage;
  bool                m_transMod;
  bool                m_useSlut;
  float               m_modTransparency;
  unsigned int        m_currentVolume;
};

ImageDisplaySetting::Handle 
ImageDisplaySetting::create(Image::Handle image,
			    LookUpTable::Handle lut,
			    float trans, 
			    bool visible)
{
  return Handle(new ImageDisplaySetting(image, lut, trans, visible));
}

ImageDisplaySetting::ImageDisplaySetting(Image::Handle image, 
                                         LookUpTable::Handle lut, 
                                         float trans, 
                                         bool visible) :
  m_impl(new Implementation(lut, trans, visible, image))
{
}

ImageDisplaySetting::ImageDisplaySetting(BriCon::Handle bricon, 
                                         LookUpTable::Handle lut, 
                                         LookUpTable::Handle slut, 
                                         float trans, 
                                         bool visible,
                                         int dti,
                                         Image::Handle mod,
                                         bool transMod,
                                         float modTransVal,
					 int vol) :
  m_impl(new Implementation(lut, slut, trans, visible, 
                            bricon,dti,mod,
                            transMod,modTransVal, vol))
{
}

ImageDisplaySetting::~ImageDisplaySetting()
{
}

void 
ImageDisplaySetting::setTransparency(float trans)
{
  m_impl->m_transparency = trans;
}

void 
ImageDisplaySetting::setModTransparency(float trans)
{
  m_impl->m_modTransparency = trans;
}

void 
ImageDisplaySetting::setVisibility(bool visible)
{
  m_impl->m_isVisible = visible;
}

void 
ImageDisplaySetting::setTransMod(bool state)
{
  m_impl->m_transMod = state;
}

void 
ImageDisplaySetting::setLookUpTable(LookUpTable::Handle lut)
{
  TRACKER("ImageDisplaySetting::setLookUpTable(LookUpTable::Handle lut)");

  m_impl->m_lookUpTable = lut;
}

void 
ImageDisplaySetting::setSecondaryLookUpTable(LookUpTable::Handle lut)
{
  m_impl->m_secondaryLookUpTable = lut;
}

void 
ImageDisplaySetting::setModImage(Image::Handle img)
{
  m_impl->m_modImage = img;
}

bool
ImageDisplaySetting::inqTransMod() const
{
  bool result(false);
  if (m_impl->m_modImage.get())result = true;
  return result;
}

Image::Handle 
ImageDisplaySetting::inqModImage() const
{
  return m_impl->m_modImage;
}

unsigned int 
ImageDisplaySetting::inqCurrentVolume() const
{
  return m_impl->m_currentVolume;
}

void 
ImageDisplaySetting::setCurrentVolume(unsigned int v)
{
  m_impl->m_currentVolume = v;
}

float 
ImageDisplaySetting::inqTransparency() const
{
  return m_impl->m_transparency;
}

float 
ImageDisplaySetting::inqModTransparency() const
{
  return m_impl->m_modTransparency;
}

bool 
ImageDisplaySetting::inqVisibility() const
{
  return m_impl->m_isVisible;
}

LookUpTable::Handle  
ImageDisplaySetting::inqLookUpTable() const
{
  return  m_impl->m_lookUpTable;
}

void 
ImageDisplaySetting::setUseSecondaryLookUpTable(bool yn)
{
  m_impl->m_useSlut = yn;
}

bool 
ImageDisplaySetting::inqUseSecondaryLookUpTable() const
{
  return m_impl->m_useSlut;
}

LookUpTable::Handle  
ImageDisplaySetting::inqSecondaryLookUpTable() const
{
  return m_impl->m_secondaryLookUpTable;
}

BriCon::Handle  
ImageDisplaySetting::inqBriCon()
{
  return  m_impl->m_bricon;
}

int  
ImageDisplaySetting::inqDtiDisplay() const
{  
  return  m_impl->m_dtiDisplay;
}

void  
ImageDisplaySetting::setDtiDisplay(DtiDisplay mode)
{  
  m_impl->m_dtiDisplay = mode;
}

ImageDisplaySetting::Handle  
ImageDisplaySetting::clone()
{
  Handle clone(new ImageDisplaySetting(m_impl->m_bricon->clone(),
				       m_impl->m_lookUpTable,
				       m_impl->m_secondaryLookUpTable,
				       m_impl->m_transparency,
				       m_impl->m_isVisible,
				       m_impl->m_dtiDisplay,
				       m_impl->m_modImage,
				       m_impl->m_transMod,
				       m_impl->m_modTransparency,
				       m_impl->m_currentVolume));

  return clone;
}
