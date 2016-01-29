/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Rama Aravind Vorray
		James Saunders
		David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2002-2005 University of Oxford  */

/*  CCOPYRIGHT */

#include "preferences.h"
#include <stdlib.h>
#include <iostream>

#include <qvariant.h>

using namespace std;

Preferences::Preferences()
{
  setPath("fmrib.ox.ac.uk", "fslview", User);
}

string Preferences::inqAtlasPath() const
{
  bool ok;
  if(m_atlaspath == "") {
    QStringList l = readListEntry("/fsl/atlaspath", &ok);
    if(ok)
      m_atlaspath = l.join(":").ascii();
  }

  if(m_atlaspath == "")
    m_atlaspath = string(getenv("FSLATLASPATH") ? 
			 getenv("FSLATLASPATH") :
			 inqFSLDir() + "/data/atlases");

  return m_atlaspath;
}

//! @brief Prefered value for window geometry
//!
//! Returns prefered values from Qt prefs or heuristically determines
//! the default placement from passed in desktop width and height measures
//!
//! @param dw width of desktop
//! @param dh height of desktop
//!
//! @return QRect The require geometry
QRect Preferences::inqGeometry(int dw, int dh) const
{
  QRect result;
    
  int x( readNumEntry("/fslview/geometry/x", -1) );
  int y( readNumEntry("/fslview/geometry/y", -1) );
  int w( readNumEntry("/fslview/geometry/width", -1) );
  int h( readNumEntry("/fslview/geometry/height", -1) );

  if( (x == -1) || (y == -1) || (w == -1) || (h == -1) ) {
    if(dh <= 800) {
      h = (dh * 85) / 100;
      w = (dw * 85) / 100;
      x = (dw - w)/2;
      y = (dh - h)/2;
    } else {
      h = (dh * 90) / 100;
      w = (h  * 80) / 100;
      x = y = (dh * 5) / 100;
    }
  }

  return QRect(x, y, w, h);
}

void Preferences::setGeometry(const QRect& r)
{
  writeEntry( "/fslview/geometry/x", r.x() );
  writeEntry( "/fslview/geometry/y", r.y() );
  writeEntry( "/fslview/geometry/width", r.width() );
  writeEntry( "/fslview/geometry/height", r.height() );
}

//! @brief Prefered value of FSLATLASPATH
//!
//! Returns vector of the prefered values of FSLATLASPATH or FSLDIR/lib/atlases
//!
//! @return The locations to look for atlas data sets
vector<string> Preferences::inqAtlasPathElements() const
{
  vector<string> result;
  string delimiters(":");

  string str(inqAtlasPath());
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  string::size_type pos = str.find_first_of(delimiters, lastPos);

  while (string::npos != pos || string::npos != lastPos)
    {
      result.push_back( str.substr(lastPos, pos - lastPos) );
      lastPos = str.find_first_not_of(delimiters, pos);
      pos = str.find_first_of(delimiters, lastPos);
    }
  
  return result;
}
   
//! @brief Prefered value of FSLDIR
//!
//! Returns the prefered value of FSLDIR
//!
//! @return The users prefered value of FSLDIR
string Preferences::inqFSLDir() const 
{
  if(m_fsldir == "")
    m_fsldir = readEntry("/fsl/fsldir", "").ascii();
  if(m_fsldir == "")
    m_fsldir = string(getenv("FSLDIR") ? getenv("FSLDIR") : "/usr/local/fsl");

  return m_fsldir; 
}

//! @brief Prefered location of MNI152 T1 brain
//!
//! Returns the prefered location where we can find MNI152 T1 brain image
//!
//! @return The users prefered location for the MNI152 T1 brain image
string Preferences::inqMni152() const
{
  if(m_mni == "")
    m_mni = readEntry("/fsl/mni","").ascii();
  if(m_mni == "")
    m_mni = inqFSLDir() + "/data/standard/MNI152_T1_2mm_brain.nii.gz";

  return m_mni;
}

string Preferences::inqAssistantPath() const
{
  if(m_assistantpath == "")
    m_assistantpath = readEntry("/qt/assistantpath","").ascii();
  if(m_assistantpath == "")
    m_assistantpath = string(getenv("FSLQTASSISTANTPATH") ? 
			     getenv("FSLQTASSISTANTPATH") : "");
  if(m_assistantpath == "")
    m_assistantpath = string(getenv("QTDIR") ? 
			     string(getenv("QTDIR")) + "/bin" : 
			     inqFSLDir() + "/bin");
  
  return m_assistantpath;
}

void Preferences::setFSLDir(const std::string& dir)
{
  m_fsldir = dir;
  writeEntry("/fsl/fsldir", dir);
}

void Preferences::setMni152(const std::string& filename)
{
  m_mni = filename;
  writeEntry("/fsl/mni", filename);
}

void Preferences::setAssistantPath(const std::string& path)
{
  m_assistantpath = path;
  writeEntry("/qt/assistantpath", path);
}

void Preferences::setAtlasPath(const std::string& path)
{
  m_atlaspath = path;
  writeEntry("/fsl/atlaspath", QStringList::split(":", path));
}
