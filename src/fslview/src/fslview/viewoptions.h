#ifndef VIEWOPTIONS_H
#define VIEWOPTIONS_H

class ViewOptions {

public:
  ViewOptions(): 
    m_localVolIdx(true), m_globalVolIdx(true), m_globalLocIdx(true),
    m_labels(true), m_sliceLabels(false), m_voxels(true),
    m_showGap(true), m_gapSize(1), m_movieFrameRate(50)
  {
  }
  
  ViewOptions& operator=(const ViewOptions& rhs)
  {
    if (this != &rhs) {
      m_localVolIdx = rhs.m_localVolIdx;
      m_globalVolIdx = rhs.m_globalVolIdx;
      m_globalLocIdx =  rhs.m_globalLocIdx;
      m_labels = rhs.m_labels;
      m_sliceLabels = rhs.m_sliceLabels;
      m_voxels = rhs.m_voxels;
      m_showGap = rhs.m_showGap;
      m_gapSize = rhs.m_gapSize;
      m_movieFrameRate = rhs.m_movieFrameRate;
    }
    return *this;
  }

  int  inqMovieFrameRate() const { return m_movieFrameRate; }
  void setMovieFrameRate(int fr) { m_movieFrameRate = fr;   }

  int  inqCursorGapSize() const { return m_gapSize; }
  void setCursorGapSize(int sz) { m_gapSize = sz;   }
  bool inqShowCursorGap() const { return m_showGap; }
  void setShowCursorGap(bool y) { m_showGap = y;    }

  bool inqVolumeIndexingWithinView() const
  { 
    return m_localVolIdx; 
  }
  void setVolumeIndexingWithinView(bool y)
  {
    m_localVolIdx = y;
  }

  bool inqUseSharedVolume() const
  { 
    return m_globalVolIdx; 
  }
  void setUseSharedVolume(bool y)
  {
    m_globalVolIdx = y;
  }

  bool inqUseSharedLocation() const
  {
    return m_globalLocIdx;
  }
  void setUseSharedLocation(bool y)
  {
    m_globalLocIdx = y;
  }

  bool inqShowLabels() const
  {
    return m_labels;
  }
  void setShowLabels(bool y)
  {
    m_labels = y;
  }

  bool inqShowSliceLabels() const
  {
    return m_sliceLabels;
  }
  void setShowSliceLabels(bool y)
  {
    m_sliceLabels = y;
  }

  bool inqUnitsAreVoxels() const
  {
    return m_voxels;
  }
  void setUnitsAreVoxels(bool y)
  {
    m_voxels = y;
  }

private:
  bool m_localVolIdx;
  bool m_globalVolIdx;
  bool m_globalLocIdx;
  bool m_labels;
  bool m_sliceLabels;
  bool m_voxels;
  bool m_showGap;
  int  m_gapSize;
  int  m_movieFrameRate;
};

#endif
