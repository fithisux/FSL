/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Rama Aravind Vorray
		James Saunders
		David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2002-2005 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(_SLICEVIEW_H)
#define _SLICEVIEW_H

#include "sliceviewbase.h"

#include <string>

class OverlayList;
class SliceWidget;

class SliceView: public SliceViewBase
{
  Q_OBJECT
public:
  typedef enum { Enabled, Greyed, Disabled } LabelState;

  SliceView(QWidget *, const char *);
  
  void setSliceWidget(SliceWidget *);

  void setNorthText(const std::string& s);
  void setSouthText(const std::string& s);
  void setEastText(const std::string& s);
  void setWestText(const std::string& s);

  void setLabelsState(LabelState s);

  QPixmap getPixmap() const;

private:
  SliceWidget *m_slice;
};

#endif

