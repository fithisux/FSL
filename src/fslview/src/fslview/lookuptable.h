
/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(LOOKUPTABLE_H)
#define LOOKUPTABLE_H

#include <boost/shared_ptr.hpp>
#include <string>

#include <vector>
#include <map>
#include "stdio.h"

typedef unsigned char ColorTriplet[3];
typedef unsigned char ColorRGBA[4];

class LutElement
{
public:
  LutElement(unsigned int i = -1, const std::string& s = "Unknown") : m_index(i), m_label(s),
								      m_r(0), m_g(0), m_b(0), m_a(255)
  {
  }

  void setColours(unsigned char r, unsigned char g, unsigned char b) {
    m_r = r; m_g = g; m_b = b; m_a = 255;
  }
  void setColours(unsigned char r, unsigned char g, unsigned char b, 
		  unsigned char a) {
    m_r = r; m_g = g; m_b = b; m_a = a;
  }

  const unsigned int  index() const { return m_index; }
  const std::string&  label() const { return m_label; }
  const unsigned char red()   const { return m_r; }
  const unsigned char green() const { return m_g; }
  const unsigned char blue()  const { return m_b; }
  const unsigned char alpha() const { return m_a; }

private:
  unsigned int  m_index;
  std::string   m_label;
  unsigned char m_r, m_g, m_b;
  unsigned char m_a;
};

/**
 * @author James Saunders <jim@fmrib.ox.ac.uk>
 * @date   Mon Dec 23 17:59:36 2002
 * 
 * @brief  Implementation of a color look up table (lut).
 * 
 * Provides a color lut suitable for GL rendering along with methods
 * to create standard luts and read custom ones from file.
 */
class LookUpTable
{
public:
  typedef boost::shared_ptr< LookUpTable > Handle;
  typedef std::vector<LutElement>::const_iterator ConstIterator;
  typedef std::vector<LutElement>::size_type SizeType;

  static LookUpTable::Handle load(const std::string& filename); 

  static LookUpTable::Handle greyScale();
  static LookUpTable::Handle red();
  static LookUpTable::Handle blue();
  static LookUpTable::Handle green();
  static LookUpTable::Handle yellow();
  static LookUpTable::Handle redYellow();
  static LookUpTable::Handle blueLightblue();
  static LookUpTable::Handle pink();
  static LookUpTable::Handle hot();
  static LookUpTable::Handle cool();
  static LookUpTable::Handle copper();
  static LookUpTable::Handle spectrum();
  static LookUpTable::Handle render1();
  static LookUpTable::Handle render1t();
  static LookUpTable::Handle render2();
  static LookUpTable::Handle render2t();
  static LookUpTable::Handle render3();
  static LookUpTable::Handle cortical();
  static LookUpTable::Handle subcortical();
  static LookUpTable::Handle rainbow();

  void pushValue(const LutElement&);
  void pushValue(unsigned char red, unsigned char green, unsigned char blue, int index);
  void pushValue(unsigned char red, unsigned char green, unsigned char blue, 
		 const std::string& label, int index);
  const LutElement& inqValue(float f);  
  const LutElement& inqValueIndex(float f);

  void allocateMemory(int size);

  void setLutName(const std::string&);
  std::string inqLutName() const;
  std::string getLabelByIndex(int n) const;

  void setVolumeName(int, const std::string&);
  std::string inqVolumeName(int) const;

  bool isIndexLut() const;
  void isIndexLut(bool);

  bool isVisible() const;
  bool isAutoSelectable() const;
  
  ConstIterator begin() const { return m_lookUpData.begin(); }
  ConstIterator end() const { return m_lookUpData.end(); }
  SizeType size() const { return m_lookUpData.size(); }
  virtual ~LookUpTable() {}

private:   
  LookUpTable(const std::string& filename);
  LookUpTable();
  LookUpTable::Handle LoadStdLut(const char ** lutData,
                                 LookUpTable::Handle lut,
                                 int elementCount,
                                 const std::string  name);
  std::string extractName(std::string filename);

  std::vector<LutElement>  m_lookUpData;
  std::map<int, std::string> m_volLabels;
  LutElement               m_black;
  std::string              m_fileName;
  std::string              m_lutName;
  bool                     m_isIndexLut;
  bool                     m_isAutoSelectable;
  bool                     m_isVisible;
};


inline std::string LookUpTable::inqLutName() const { return m_lutName; }
inline bool LookUpTable::isIndexLut() const { return m_isIndexLut; }
inline void LookUpTable::isIndexLut(bool y) { m_isIndexLut = y; }
inline bool LookUpTable::isVisible() const { return m_isVisible; }
inline bool LookUpTable::isAutoSelectable() const { return m_isAutoSelectable; }

#endif
