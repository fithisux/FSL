/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Rama Aravind Vorray
		James Saunders
		David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2002-2005 University of Oxford  */

/*  CCOPYRIGHT */

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>

#include <qsettings.h>
#include <qrect.h>

class Preferences: public QSettings
{
public:
  std::string inqFSLDir() const;
  std::string inqMni152() const;
  std::string inqAssistantPath() const;
  std::string inqAtlasPath() const;
  QRect inqGeometry(int, int) const;

  std::vector<std::string> inqAtlasPathElements() const;

  void setFSLDir(const std::string&);
  void setMni152(const std::string&);
  void setAssistantPath(const std::string&);
  void setAtlasPath(const std::string&);
  void setGeometry(const QRect&);

  Preferences();
  virtual ~Preferences() {}
private:

  mutable std::string m_assistantpath;
  mutable std::string m_atlaspath;
  mutable std::string m_fsldir;
  mutable std::string m_mni;
  int m_w, m_h, m_x, m_y;
};

