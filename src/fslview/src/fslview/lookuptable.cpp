/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "lookuptable.h"
#include "filemanager.h"
#include <qfileinfo.h>
#include <qcolor.h>

#include <assert.h>
#include <math.h>

#include "luts/render1.ldt"
#include "luts/render1t.ldt"
#include "luts/render2.ldt"
#include "luts/render2t.ldt"
#include "luts/render3.ldt"
#include "luts/cortical.ldt"
#include "luts/subcortical.ldt"

#include <stdexcept>

using namespace std;

LookUpTable::LookUpTable(const std::string& filename):
  m_fileName(filename), m_isIndexLut(false), m_isAutoSelectable(false), m_isVisible(true)
{
  m_black.setColours(0, 0, 0);
  allocateMemory(255);
  for(unsigned char i = 0; i < 255; ++i)
    {
      int value = int((i / (255.0 - 1.0)) * 255);
      pushValue(value, value, value, i);
    }
  m_lutName = std::string("Greyscale");
  m_isAutoSelectable = false;

  QFileInfo fi(QString::fromLatin1(filename.c_str()));
  if(!fi.exists())
    throw std::runtime_error("LookUpTable: no such file: " + filename);
  
  std::string ext( fi.extension().ascii() );

  if(ext == "lml") {
    FileManager::readLMLFile(filename, this);
  } else {
    if((ext != "rgb") && (ext != "lut"))
      throw std::runtime_error("LookUpTable: lut file type unknown extension: " + fi.extension());

    if(fi.extension() == "rgb")
      m_isIndexLut = true;
 
    if(m_isIndexLut)
      FileManager::readRgbFile(m_fileName,this);
    else
      FileManager::readLutFile(m_fileName,this);

    m_lutName = extractName(m_fileName);
  }

}

LookUpTable::LookUpTable():
  m_fileName(""),
  m_isIndexLut(false), m_isAutoSelectable(true), m_isVisible(true)
{
  m_black.setColours(0, 0, 0);
}

/** 
 * Read a lut from a file.
 * 
 * @param filename The file containing the MEDx style lut definition.
 * 
 * @return Handle to the new lut.
 */
LookUpTable::Handle LookUpTable::load(const std::string& filename)
{
  return LookUpTable::Handle(new LookUpTable(filename)); 
} 

/** 
 * Generate a grey-scale look up table.
 * 
 * 
 * @return A lut populated from black to white.
 */
LookUpTable::Handle LookUpTable::greyScale()
{
  LookUpTable::Handle lut = LookUpTable::Handle(new LookUpTable());
  
  lut->allocateMemory(255);

  for(unsigned char i = 0; i < 255; ++i)
    {
      int value = int((i / (255.0 - 1.0)) * 255);
      lut->pushValue(value, value, value, i);
    }
  
  lut->m_lutName = std::string("Greyscale");
  lut->m_isAutoSelectable = false;
  return lut;
}

/** 
 * Generate a red look up table.
 * 
 * 
 * @return A lut populated from black to red.
 */
LookUpTable::Handle LookUpTable::red()
{
  LookUpTable::Handle lut = LookUpTable::Handle(new LookUpTable());
  
  lut->allocateMemory(255);

  float offset = 100.0;
  float step = (255.0 - offset) / 255.0;

  for(unsigned char i = 0; i < 255; ++i)
    {
      int red = int(((i + 1) * step) + offset);
      lut->pushValue(red, 0, 0, i);
    }
  
  lut->m_lutName = std::string("Red");

  return lut;
}

/** 
 * Generate a green look up table.
 * 
 * 
 * @return a lut populated from black to green.
 */
LookUpTable::Handle LookUpTable::green()
{
  LookUpTable::Handle lut = LookUpTable::Handle(new LookUpTable());
  
  lut->allocateMemory(255);

  float offset = 100.0;
  float step = (255.0 - offset) / 255.0;

  for(unsigned char i = 0; i < 255; ++i)
    {
      int green = int(((i + 1) * step) + offset);
      lut->pushValue(0, green, 0, i);
    }

  lut->m_lutName = std::string("Green");

  return lut;
}

/** 
 * Generate a blue look up table.
 * 
 * 
 * @return a lut populated from black to blue.
 */
LookUpTable::Handle LookUpTable::blue()
{
  LookUpTable::Handle lut = LookUpTable::Handle(new LookUpTable());
  
  float offset = 100.0;
  float step = (255.0 - offset) / 255.0;

  lut->allocateMemory(255);

  for(unsigned char i = 0; i < 255; ++i)
    {
      int blue = int(((i + 1) * step) + offset);
      lut->pushValue(0, 0, blue, i);
    }

  lut->m_lutName = std::string("Blue");

  return lut;
}

/** 
 * Generate a yellow look up table.
 * 
 * 
 * @return a lut populated from black to yellow.
 */
LookUpTable::Handle LookUpTable::yellow()
{
  LookUpTable::Handle lut = Handle(new LookUpTable());
  
  lut->allocateMemory(255);

  float offset = 100.0;
  float step = (255.0 - offset) / 255.0;

  for(unsigned char i = 0; i < 255; ++i)
    {
      int red = int((i * step) + offset);
      int green = int((i * step) + offset);
      lut->pushValue(red, green, 0, i);
    }

  lut->m_lutName = std::string("Yellow");

  return lut;
}

/** 
 * Generate a red-yellow look up table.
 * 
 * 
 * @return a lut populated from red to yellow.
 */
LookUpTable::Handle LookUpTable::redYellow()
{
  LookUpTable::Handle lut = Handle(new LookUpTable());
  
  lut->allocateMemory(255);

  for(unsigned char i = 0; i < 255; ++i)
    {
      int green = int((i / (255.0 - 1.0)) * 255);
      lut->pushValue(255, green, 0, i);
    }

  lut->m_lutName = std::string("Red-Yellow");

  return lut;
}

/** 
 * Generate a blue-lightblue look up table.
 * 
 * 
 * @return a lut populated from blue to lightblue
 */
LookUpTable::Handle LookUpTable::blueLightblue()
{
  LookUpTable::Handle lut = Handle(new LookUpTable());
  
  lut->allocateMemory(255);

  for(unsigned char i = 0; i < 255; ++i)
    {
      int green = int((i / (255.0 - 1.0)) * 255);
      lut->pushValue(0, green, 255, i);
    }

  lut->m_lutName = std::string("Blue-Lightblue");

  return lut;
}

/** 
 * Generate a pink look up table.
 * 
 * 
 * @return a lut populated with shades of pink
 */
LookUpTable::Handle LookUpTable::pink()
{
  LookUpTable::Handle lut = LookUpTable::Handle(new LookUpTable());
  
  lut->allocateMemory(255);

  float offset = 100.0;
  float step = (255.0 - offset) / 255.0;

  for(unsigned char i = 0; i < 255; ++i)
    {
      int v = int((i * step) + offset);
      lut->pushValue(255, v, v, i);
    }

  lut->m_lutName = std::string("Pink");

  return lut;
}

/** 
 * Generate a "hot" look up table.
 * 
 * 
 * @return a lut
 */
LookUpTable::Handle LookUpTable::hot()
{
  LookUpTable::Handle lut = LookUpTable::Handle(new LookUpTable());
  
  lut->allocateMemory(255);

  float offset = 100.0;
  float step = (255.0 - offset) / 85.0;

  for(unsigned char i = 0; i < 85; ++i)
    {
      int v = int((i * step) + offset);
      lut->pushValue(  v,   0,   0,     i);
    }
  for(unsigned char i = 0; i < 85; ++i)
    {
      int v = int((i * step) + offset);
      lut->pushValue(255,   v,   0,  i+85);
    }
  for(unsigned char i = 0; i < 85; ++i)
    {
      int v = int((i * step) + offset);
      lut->pushValue(255, 255,   v, i+170);
    }

  lut->m_lutName = std::string("Hot");

  return lut;
}

/** 
 * Generate a "copper" look up table.
 * 
 * 
 * @return a lut
 */
LookUpTable::Handle LookUpTable::copper()
{
  LookUpTable::Handle lut = LookUpTable::Handle(new LookUpTable());
  
  lut->allocateMemory(255);

  float rstep = 1.2;
  float gstep = 0.8;
  float bstep = 0.5;

  for(unsigned char i = 0; i < 255; ++i)
    {
      int r = std::min(255, int(i * rstep));
      int g = int(i * gstep);
      int b = int(i * bstep);
      lut->pushValue(r, g, b, i);
    }

  lut->m_lutName = std::string("Copper");

  return lut;
}

/** 
 * Generate a "cool" look up table.
 * 
 * 
 * @return a lut
 */
LookUpTable::Handle LookUpTable::cool()
{
  LookUpTable::Handle lut = LookUpTable::Handle(new LookUpTable());
  
  lut->allocateMemory(255);

  for(unsigned char i = 0; i < 255; ++i)
    {
      lut->pushValue(i, 255-i, 255, i);
    }

  lut->m_lutName = std::string("Cool");

  return lut;
}

class Random {
public:
  Random() { srand(1); }

  ptrdiff_t operator() (ptrdiff_t max) {
    double tmp;
    tmp = static_cast<double>(rand())
      / static_cast<double>(RAND_MAX);
    return static_cast<ptrdiff_t>(tmp * max);
  }
};

LookUpTable::Handle LookUpTable::rainbow()
{
  LookUpTable::Handle lut = Handle(new LookUpTable());

  lut->allocateMemory(1200);

  float step = 360.0/1200.0;
  for(unsigned int i = 0; i < 1200; ++i) 
    {
      int h(i * step), r, g, b;
      QColor c(h, 255, 255, QColor::Hsv);
      c.getRgb(&r, &g, &b);
      lut->pushValue(r, g, b, i);
    }
  Random rd;
  random_shuffle(lut->m_lookUpData.begin(), lut->m_lookUpData.end(), rd);

  lut->m_lutName = string("Random-Rainbow");
  lut->m_isVisible = true;
  lut->m_isIndexLut = true;
  return lut;
}

LookUpTable::Handle LookUpTable::cortical()
{
  LookUpTable::Handle lut = Handle(new LookUpTable());
  lut->LoadStdLut(cortical_lut, lut, 58, "MGH-Cortical");
  lut->m_isVisible = true;
  lut->m_isIndexLut = true;
  return lut;
}

LookUpTable::Handle LookUpTable::subcortical()
{
  LookUpTable::Handle lut = Handle(new LookUpTable());
  lut->LoadStdLut(subcortical_lut, lut, 89, "MGH-Subcortical");
  lut->m_isVisible = true;
  lut->m_isIndexLut = true;
  return lut;
}

LookUpTable::Handle LookUpTable::render1()
{
  LookUpTable::Handle lut = Handle(new LookUpTable());
  lut->LoadStdLut(render1_lut, lut, 200, "render1");
  return lut;
}

LookUpTable::Handle LookUpTable::render1t()
{
  LookUpTable::Handle lut = Handle(new LookUpTable());
  lut->LoadStdLut(render1t_lut, lut, 200, "render1t");
  return lut;
}

LookUpTable::Handle LookUpTable::render2()
{
  LookUpTable::Handle lut = Handle(new LookUpTable());
  lut->LoadStdLut(render2_lut, lut, 500, "render2");
  return lut;
}

LookUpTable::Handle LookUpTable::render2t()
{
  LookUpTable::Handle lut = Handle(new LookUpTable());
  lut->LoadStdLut(render2t_lut, lut, 500, "render2t");
  return lut;
}

LookUpTable::Handle LookUpTable::render3()
{
  LookUpTable::Handle lut = Handle(new LookUpTable());
  lut->LoadStdLut(render3_lut, lut, 300, "render3");
  return lut;
}


LookUpTable::Handle LookUpTable::LoadStdLut(const char** lutData,
                                            LookUpTable::Handle lut,
                                            int elementCount,
                                            const std::string name)
{
  float red, green, blue;

  for(int i = 0; i < elementCount; ++i)
    {
      sscanf(lutData[i],"%f,%f,%f",&red,&green,&blue);
      lut->pushValue((unsigned char)(red*255.0),
                     (unsigned char)(green*255.0),
                     (unsigned char)(blue*255.0),i);
    }

  lut->m_lutName = name;
  //  lut->m_isVisible = false;
  lut->m_isAutoSelectable = false;
  return lut;
}




/** 
 * Get a value from the lut.
 * 
 * @param f The value to be converted [0..1]
 * 
 * @return A color-triplet value from the look up table.
 */
const LutElement& LookUpTable::inqValue(float f)
{
  int size(m_lookUpData.size());
  int index = std::min(size - 1, std::max(0, int (ceil(f * (size-1)))));

  return m_lookUpData.at(index);
}

const LutElement& LookUpTable::inqValueIndex(float f)
{
  SizeType index((int)f);
  if (index >= m_lookUpData.size())
      return m_black;

  return m_lookUpData.at(index);
}

void LookUpTable::allocateMemory(int size)
{
  m_lookUpData.clear();
}

void LookUpTable::pushValue(const LutElement& e)
{
  m_lookUpData.push_back(e);
}

void LookUpTable::pushValue(unsigned char red, unsigned char green, unsigned char blue,
                            const std::string& label, int index)
{
  LutElement e(index, label);
  e.setColours(red, green, blue);

  m_lookUpData.push_back(e);
}

void LookUpTable::pushValue(unsigned char red, unsigned char green, unsigned char blue,
                            int index)
{
  LutElement e(index);
  e.setColours(red, green, blue);

  m_lookUpData.push_back(e);
}

std::string LookUpTable::getLabelByIndex(int n) const
{
  if(isIndexLut() && (n >= 0) && (n < int(m_lookUpData.size())))
    return m_lookUpData.at(n).label();
  else
    return std::string("No labels available!");
}

void LookUpTable::setVolumeName(int n, const std::string& l)
{
//   m_volLabels[n] = l;
  m_volLabels.insert(std::make_pair(n, l));
}

std::string LookUpTable::inqVolumeName(int n) const
{
  std::map<int, std::string>::const_iterator pos;

  pos = m_volLabels.find(n);
  if(pos != m_volLabels.end())
    return pos->second;
  else
    return std::string("No label found!");
}

std::string LookUpTable::extractName(std::string filename)
{
  std::string result;
  std::string ext;

  result.erase();
  result = filename.substr(filename.rfind('/')+ 1,filename.length()- 1);
  result = result.substr(0,result.rfind('.'));

  return result;
}

void LookUpTable::setLutName(const std::string& name)
{
  m_lutName = name;
}

