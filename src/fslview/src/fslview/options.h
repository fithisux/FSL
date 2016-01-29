/*  FSLView - 2D/3D Interactive Image Viewer

    David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(FSLVIEWOPTIONS_H)
#define FSLVIEWOPTIONS_H

#include <exception>
#include <string>

namespace Utilities {

  extern bool string_to_T(std::pair<float,float>&, const std::string&);

}

#include "utils/options.h"
#include "fslio/fslio.h"

#include <qfileinfo.h>

class OverlayOption
{
public:
  OverlayOption(const string& filename, 
		const Utilities::Option<string>& lutname, 
		const Utilities::Option<float>& trans, 
		const Utilities::Option< std::pair<float,float> >& bricon):
  m_filename(filename), m_lutname(lutname), m_trans(trans), m_bricon(bricon)
  {
  }
  
  QFileInfo fileInfo() const
  {
    return QFileInfo(FslMakeBaseName(m_filename.c_str()));
  }

  bool lutSpecified() const { return m_lutname.set(); }
  const string& lutname() const { return m_lutname.value(); }

  bool briconSpecified() const { return m_bricon.set(); }
  float min() const { return m_bricon.value().first; }
  float max() const { return m_bricon.value().second; }

  bool transparencySpecified() const { return m_trans.set(); }
  float transparency() const { return m_trans.value(); }
  
  virtual ~OverlayOption() {}

private:
  //  OverlayOption() {}

  string m_filename;
  Utilities::Option<string> m_lutname;
  Utilities::Option<float> m_trans;
  Utilities::Option< std::pair<float,float> > m_bricon;
};

typedef std::list<OverlayOption> OverlayOptionList;

class ApplicationOptions
{
public:
  typedef enum {Ortho = 0, Single, Lightbox, ThreeD} Mode;

  ApplicationOptions() {}

  OverlayOptionList::const_iterator begin() const { return m_overlays.begin(); }
  OverlayOptionList::const_iterator end() const { return m_overlays.end(); }
  bool empty() const { return m_overlays.empty(); }
  void push_back(const OverlayOption& o) { m_overlays.push_back(o); }

  void setModes(const std::vector<string> modeStrings)
  {
    for(std::vector<string>::const_iterator it = modeStrings.begin();
	it != modeStrings.end(); ++it)
      {
	if(!it->compare("ortho")) { 
	  m_modes.push_back(Ortho); 
	} else if(!it->compare("lightbox")) { 
	  m_modes.push_back(Lightbox); 
	} else if(!it->compare("single")) { 
	  m_modes.push_back(Single); 
	} else if(!it->compare("3d")) { 
	  m_modes.push_back(ThreeD); 
	} else {
	  std::string msg(*it + ": bad mode!");
	  throw std::runtime_error(msg.c_str());
	}
      }
  }

  std::list<Mode>& inqModes()
  {
    return m_modes;
  }

private:
  std::list<Mode>   m_modes;
  OverlayOptionList m_overlays;
};
#endif
