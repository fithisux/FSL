/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

//#define DEBUGGING

#include <math.h>
#include <iostream>

#include "imagebuffer.h"
#include "tracker.h"

#include <qobject.h>

struct ColorPicker
{
  ColorPicker(const MetaImage::Handle mi)
  {
    m_bricon = mi->getDs()->inqBriCon();  
    m_lut = mi->getDs()->inqLookUpTable();
  }

  inline void byIndex(float v)
  {
    const LutElement& col = m_lut->inqValueIndex(v);
    m_rd = col.red(); m_gr = col.green(); m_bl = col.blue();
    m_al = (v > 0.0) ? 255 : 0;
  }

  inline void byLookUp(float v)
  {
    float adjustedValue = m_bricon->adjust(v);
    const LutElement& col = m_lut->inqValue(adjustedValue);
    m_rd = col.red(); m_gr = col.green(); m_bl = col.blue();
    m_al = (adjustedValue > 1e-15) ? 255 : 0;     
  }

  inline void byIntensity(float v)
  {
    float adjustedValue = m_bricon->adjust(v);
    m_rd = m_bl = m_gr = ImageBuffer::clamp(adjustedValue);
    m_al = (m_rd > 0) ? 255 : 0;    
  }

  inline int red()   const { return m_rd; }
  inline int green() const { return m_gr; }
  inline int blue()  const { return m_bl; }
  inline int alpha() const { return m_al; }

  unsigned short m_rd, m_gr, m_bl, m_al; 
  BriCon::Handle m_bricon;
  LookUpTable::Handle m_lut;
};

ColorRGBAHandle ImageBuffer::axialBuffer(MetaImage::Handle mi, 
                                         int slice, int vol)
{
  STATIC_TRACKER("ImageBuffer::axialBuffer");

  Volume::Handle v = mi->getImage()->getVolume(vol);
  unsigned int width  = v->inqX();
  unsigned int height = v->inqY();

  ColorRGBAHandle buffer = 
    ColorRGBAHandle( new ColorRGBA[width*height] );

  //ColorPicker color(mi);
  LookUpTable::Handle lut = mi->getDs()->inqLookUpTable();
  BriCon::Handle   bricon = mi->getDs()->inqBriCon();  
  unsigned short rd, gr, bl, al;

  if(lut)
    if(lut->isIndexLut()) {
      MESSAGE("Index LUT");
      for(unsigned int y = 0; y < height; ++y)
	{ 
	  unsigned int offset = width * y;
	  
	  for( unsigned int x = 0; x < width; ++x) 
	    {              
	      unsigned int xloc = offset + x;
	      
	      float origValue =  v->value(x, y, slice);              
	      
	      //color.byIndex(origValue);
	      const LutElement& col = lut->inqValueIndex(origValue);
	      rd = col.red(); gr = col.green(); bl = col.blue();

	      al = (origValue > 0.0) ? 255 : 0;
	      
	      buffer[xloc][2] = rd;
	      buffer[xloc][1] = gr;
	      buffer[xloc][0] = bl;
	      buffer[xloc][3] = al;
	    }
	}
    } else {
      MESSAGE("RGB LUT");
      for(unsigned int y = 0; y < height; ++y)
	{ 
	  unsigned int offset = width * y;
	  
	  for( unsigned int x = 0; x < width; ++x) 
	    {              
	      unsigned int xloc = offset + x;
	      
	      float origValue =  v->value(x, y, slice);              
	      
	      //color.byLookUp(origValue);
	      float adjustedValue = bricon->adjust(origValue);
	      const LutElement& col = lut->inqValue(adjustedValue);
	      rd = col.red(); gr = col.green(); bl = col.blue();

	      al = (adjustedValue > 1e-15) ? 255 : 0;     
	      
	      buffer[xloc][2] = rd;
	      buffer[xloc][1] = gr;
	      buffer[xloc][0] = bl;
	      buffer[xloc][3] = al;
	    }
	}
    }
  else
    {
      MESSAGE("Greyscale");
      for(unsigned int y = 0; y < height; ++y)
	{ 
	  unsigned int offset = width * y;
	
	  for( unsigned int x = 0; x < width; ++x) 
	    {              
	      unsigned int xloc = offset + x;
	    
	      float origValue =  v->value(x, y, slice);              
	      
	      //color.byIntensity(origValue);
	      unsigned short i = ImageBuffer::clamp(origValue);
	      
	      buffer[xloc][2] = i;
	      buffer[xloc][1] = i;
	      buffer[xloc][0] = i;
	      buffer[xloc][3] = (i > 0) ? 255 : 0;
	    }
	}
    }
  
  return buffer;
}

ColorRGBAHandle ImageBuffer::coronalBuffer(MetaImage::Handle mi,
                                           int slice, int vol)
{
  STATIC_TRACKER("ImageBuffer::coronalBuffer");

  Volume::Handle v = mi->getImage()->getVolume(vol); 
  unsigned int width  = v->inqX();
  unsigned int height = v->inqZ();
    
  ColorRGBAHandle buffer = 
    ColorRGBAHandle( new ColorRGBA[width*height] );

  //ColorPicker color(mi);
  LookUpTable::Handle lut = mi->getDs()->inqLookUpTable();
  BriCon::Handle   bricon = mi->getDs()->inqBriCon();  
  unsigned short rd, gr, bl, al;

  if(lut)
    if(lut->isIndexLut()) {
      for(unsigned int z = 0; z < height; ++z)
	{
	  unsigned int offset = width * z;
          
	  for(unsigned int x = 0; x < width; ++x)
	    {
	      unsigned int xloc = offset + x;
	      
	      float origValue =  v->value(x, slice, z);
	      
	      //color.byIndex(origValue);
	      const LutElement& col = lut->inqValueIndex(origValue);
	      rd = col.red(); gr = col.green(); bl = col.blue();

	      al = (origValue > 0.0) ? 255 : 0;
	      
	      buffer[xloc][2] = rd;
	      buffer[xloc][1] = gr;
	      buffer[xloc][0] = bl;
	      buffer[xloc][3] = al;
	    }
	}
    } else {
      for(unsigned int z = 0; z < height; ++z)
	{
	  unsigned int offset = width * z;
          
	  for(unsigned int x = 0; x < width; ++x)
	    {
	      unsigned int xloc = offset + x;
	      
	      float origValue =  v->value(x, slice, z);
	      
	      //color.byLookUp(origValue);
	      float adjustedValue = bricon->adjust(origValue);
	      const LutElement& col = lut->inqValue(adjustedValue);
	      rd = col.red(); gr = col.green(); bl = col.blue();

	      al = (adjustedValue > 1e-15) ? 255 : 0;     
	      
	      buffer[xloc][2] = rd;
	      buffer[xloc][1] = gr;
	      buffer[xloc][0] = bl;
	      buffer[xloc][3] = al;
	    }
	}
    }
  else
    for(unsigned int z = 0; z < height; ++z)
      {
	unsigned int offset = width * z;
        
	for(unsigned int x = 0; x < width; ++x)
	  {
	    unsigned int xloc = offset + x;
	    
	    float origValue =  v->value(x, slice, z);
	      
	    //color.byIntensity(origValue);
	    unsigned char i = ImageBuffer::clamp(origValue);

	    buffer[xloc][2] = i;
	    buffer[xloc][1] = i;
	    buffer[xloc][0] = i;
	    buffer[xloc][3] = (i > 0) ? 255 : 0;
	  }
      }

  return buffer;
}

ColorRGBAHandle ImageBuffer::sagittalBuffer(MetaImage::Handle mi,
                                            int slice, int vol)
{
  STATIC_TRACKER("ImageBuffer::sagittalBuffer");

  Volume::Handle v = mi->getImage()->getVolume(vol);
  unsigned int width  = v->inqY();
  unsigned int height = v->inqZ();
    
  ColorRGBAHandle buffer = 
    ColorRGBAHandle( new ColorRGBA[width*height] );

  //ColorPicker color(mi);
  LookUpTable::Handle lut = mi->getDs()->inqLookUpTable();
  BriCon::Handle   bricon = mi->getDs()->inqBriCon();  
  unsigned short rd, gr, bl, al;

  if(lut)
    if(lut->isIndexLut()) {
      for(unsigned int z = 0; z < height; ++z)
	{
	  unsigned int offset = width * z;
      
	  for(unsigned int y = 0; y < width; ++y)
	    {
	      unsigned int xloc = offset + y;
	      
	      float origValue =  v->value(slice, y, z);
	      
	      //color.byIndex(origValue);
	      const LutElement& col = lut->inqValueIndex(origValue);
	      rd = col.red(); gr = col.green(); bl = col.blue();
	      al = (origValue > 0.0) ? 255 : 0;
	      
	      buffer[xloc][2] = rd;
	      buffer[xloc][1] = gr;
	      buffer[xloc][0] = bl;
	      buffer[xloc][3] = al;
	    }
	}
    } else {
      for(unsigned int z = 0; z < height; ++z)
	{
	  unsigned int offset = width * z;
      
	  for(unsigned int y = 0; y < width; ++y)
	    {
	      unsigned int xloc = offset + y;

	      float origValue =  v->value(slice, y, z);
	  
	      //color.byLookUp(origValue);
	      float adjustedValue = bricon->adjust(origValue);
	      const LutElement& col = lut->inqValue(adjustedValue);
	      rd = col.red(); gr = col.green(); bl = col.blue();
	      al = (adjustedValue > 1e-15) ? 255 : 0;     
	      
	      buffer[xloc][2] = rd;
	      buffer[xloc][1] = gr;
	      buffer[xloc][0] = bl;
	      buffer[xloc][3] = al;
	    }
	}
    }
  else
    for(unsigned int z = 0; z < height; ++z)
      {
	unsigned int offset = width * z;
      
	for(unsigned int y = 0; y < width; ++y)
	  {
	    unsigned int xloc = offset + y;

	    float origValue =  v->value(slice, y, z);
	  
	    //color.byIntensity(origValue);
	    unsigned char i = ImageBuffer::clamp(origValue);

	    buffer[xloc][2] = i;
	    buffer[xloc][1] = i;
	    buffer[xloc][0] = i;
	    buffer[xloc][3] = (i > 0) ? 255 : 0;
	  }
      }
  
  return buffer;
}

ColorRGBAHandle ImageBuffer::axialDtiBuffer(MetaImage::Handle mi,
                                            int slice)
{
  Volume::Handle vR = mi->getImage()->getVolume(0);
  Volume::Handle vG = mi->getImage()->getVolume(1);
  Volume::Handle vB = mi->getImage()->getVolume(2);
  
  unsigned int width  = vR->inqX();
  unsigned int height = vR->inqY();
  BriCon::Handle bricon   = mi->getDs()->inqBriCon();

  ColorRGBAHandle buffer = 
    ColorRGBAHandle( new ColorRGBA[width*height] );
  if(!mi->getDs()->inqModImage().get())
    {
      for(unsigned int y = 0; y < height; ++y)
        { 
          unsigned int offset = width * y;

          for( unsigned int x = 0; x < width; ++x) 
            {
              unsigned int xloc = offset + x;
              buffer[xloc][2] = clamp(bricon->adjust(fabs(vR->value(x, y, slice))));
              buffer[xloc][1] = clamp(bricon->adjust(fabs(vG->value(x, y, slice))));
              buffer[xloc][0] = clamp(bricon->adjust(fabs(vB->value(x, y, slice))));
              buffer[xloc][3] = 255;
            }
        }
    }
  else
    {  
      float modTransparency = mi->getDs()->inqModTransparency();

      Volume::Handle vMod = mi->getDs()->inqModImage()->getVolume(0);

      for(unsigned int y = 0; y < height; ++y)
        { 
          unsigned int offset = width * y;

          for( unsigned int x = 0; x < width; ++x) 
            {
              unsigned int xloc = offset + x;
              float modVal = vMod->normalized(x, y, slice);
              
              buffer[xloc][2] = clamp(bricon->adjust(fabs(vR->value(x, y, slice))));
              buffer[xloc][1] = clamp(bricon->adjust(fabs(vG->value(x, y, slice))));
              buffer[xloc][0] = clamp(bricon->adjust(fabs(vB->value(x, y, slice))));
              buffer[xloc][3] = clamp(modVal + modTransparency);
            }
        } 

    }
  return buffer;


}
ColorRGBAHandle ImageBuffer::coronalDtiBuffer(MetaImage::Handle mi,
                                              int slice)
{
 Volume::Handle vR = mi->getImage()->getVolume(0);
  Volume::Handle vG = mi->getImage()->getVolume(1);
  Volume::Handle vB = mi->getImage()->getVolume(2);
  
  unsigned int width  = vR->inqX();
  unsigned int height = vR->inqZ();
  BriCon::Handle bricon   = mi->getDs()->inqBriCon();

  ColorRGBAHandle buffer = 
    ColorRGBAHandle( new ColorRGBA[width*height] );
  if(!mi->getDs()->inqModImage().get())
  {       
    for(unsigned int z = 0; z < height; ++z)
        { 
          unsigned int offset = width * z;

          for( unsigned int x = 0; x < width; ++x) 
            {
              unsigned int xloc = offset + x;

              buffer[xloc][2] = clamp(bricon->adjust(fabs(vR->value(x, slice, z))));
              buffer[xloc][1] = clamp(bricon->adjust(fabs(vG->value(x, slice, z))));
              buffer[xloc][0] = clamp(bricon->adjust(fabs(vB->value(x, slice, z))));
              buffer[xloc][3] = 255;
            }
        } 
  }
  else
  {      
    float modTransparency = mi->getDs()->inqModTransparency();

    Volume::Handle vMod = mi->getDs()->inqModImage()->getVolume(0);
 
    for(unsigned int z = 0; z < height; ++z)
      { 
        unsigned int offset = width * z;

        for( unsigned int x = 0; x < width; ++x) 
          {
            unsigned int xloc = offset + x;
            float modVal = vMod->normalized(x, slice,z);
        
            buffer[xloc][2] = clamp(bricon->adjust(fabs(vR->value(x, slice, z))));
            buffer[xloc][1] = clamp(bricon->adjust(fabs(vG->value(x, slice, z))));
            buffer[xloc][0] = clamp(bricon->adjust(fabs(vB->value(x, slice, z))));
            buffer[xloc][3] = clamp(modVal + modTransparency);
          }
      }
  }

  return buffer;

}

ColorRGBAHandle ImageBuffer::sagittalDtiBuffer(MetaImage::Handle mi,
                                               int slice)
{
 Volume::Handle vR = mi->getImage()->getVolume(0);
  Volume::Handle vG = mi->getImage()->getVolume(1);
  Volume::Handle vB = mi->getImage()->getVolume(2);
  
  unsigned int width  = vR->inqY();
  unsigned int height = vR->inqZ();
  BriCon::Handle bricon   = mi->getDs()->inqBriCon();

  ColorRGBAHandle buffer = 
    ColorRGBAHandle( new ColorRGBA[width*height] );
  
  if(!mi->getDs()->inqModImage().get())
  {         
  for(unsigned int z = 0; z < height; ++z)
    { 
      unsigned int offset = width * z;

      for( unsigned int y = 0; y < width; ++y) 
        {
          unsigned int xloc = offset + y;
          
          buffer[xloc][2] = clamp(bricon->adjust(fabs(vR->value(slice, y, z))));
          buffer[xloc][1] = clamp(bricon->adjust(fabs(vG->value(slice, y, z))));
          buffer[xloc][0] = clamp(bricon->adjust(fabs(vB->value(slice, y, z))));
          buffer[xloc][3] = 255;
        }
    }
  }
  else
    {      
      float modTransparency = mi->getDs()->inqModTransparency();
     
      Volume::Handle vMod = mi->getDs()->inqModImage()->getVolume(0);

      for(unsigned int z = 0; z < height; ++z)
        { 
          unsigned int offset = width * z;

          for( unsigned int y = 0; y < width; ++y) 
            {
              unsigned int xloc = offset + y;
              
              float modVal = vMod->normalized(slice,y,z);

              buffer[xloc][2] = clamp(bricon->adjust(fabs(vR->value(slice, y, z))));
              buffer[xloc][1] = clamp(bricon->adjust(fabs(vG->value(slice, y, z))));
              buffer[xloc][0] = clamp(bricon->adjust(fabs(vB->value(slice, y, z))));
              buffer[xloc][3] = clamp(modVal + modTransparency);
            }
        }

    }
  
  return buffer;

}

void ImageBuffer::blendBuffers(ColorRGBAHandle dest, 
                               ColorRGBAHandle source,
                               float trans, bool bottomLayer, 
                               unsigned int length)
{
  float alpha(trans);
  float oneMinusAlpha(1 - alpha);
  
       if(bottomLayer)
        {
          for(unsigned int x = 0; x < length; ++x)
            {
              if(source[x][3] != 0)
                {  
                  alpha = trans * (source[x][3]/255.0);
                  
                  dest[x][0] = int(source[x][0] * alpha);
                  dest[x][1] = int(source[x][1] * alpha);          
                  dest[x][2] = int(source[x][2] * alpha);
                }
              else
                {
                  dest[x][0] = 0;
                  dest[x][1] = 0;
                  dest[x][2] = 0;
                }
            }
          
        }
      else
        {
          for(unsigned int x = 0; x < length; ++x)
            {
              if(source[x][3] != 0)
                {
                  alpha = trans * source[x][3]/255.0;
                  oneMinusAlpha = 1 - alpha;
                  dest[x][0] = int((source[x][0] * alpha) + (dest[x][0] * oneMinusAlpha));
                  dest[x][1] = int((source[x][1] * alpha) + (dest[x][1] * oneMinusAlpha));          
                  dest[x][2] = int((source[x][2] * alpha) + (dest[x][2] * oneMinusAlpha));
		  dest[x][3] = source[x][3];
                }
            }
        }     

}
  
void ImageBuffer::setToZero(ColorRGBAHandle buffer,unsigned int length)
{
  for(unsigned int x = 0; x < length; ++x)
  {
    buffer[x][0] = 0;         
    buffer[x][1] = 0;         
    buffer[x][2] = 0;         
    buffer[x][3] = 0; 
  }
	 
}

void ImageBuffer::reorderBytes(ColorRGBAHandle buffer, unsigned int length)
{
  unsigned int A,B,C,D;
  for(unsigned int x = 0; x < length; ++x)
  {
    A = buffer[x][0];
    B = buffer[x][1];
    C = buffer[x][2];
    D = buffer[x][3];
    
    buffer[x][0] = D;
    buffer[x][1] = C;         
    buffer[x][2] = B;         
    buffer[x][3] = A; 
  }
}

