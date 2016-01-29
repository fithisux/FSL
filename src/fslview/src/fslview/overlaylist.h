/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(OVERLAYLIST_H)
#define OVERLAYLIST_H

#include "metaimage.h"
#include "imagegroup.h"
#include <boost/shared_ptr.hpp>

enum OverlayListMsg {Select, Visibility, Transparency, Order,  
              Add,    Rem ,       LookUpTable,  Security, 
              ImageName,  DtiMode, ModImage};

class OverlayListObserver;

/**
 * @author David Flitney <flitney@fmrib.ox.ac.uk>
 *
 * @date   Dec 2002
 * 
 * @brief OverlayList groups images according to their display order
 * and provides convenience functions for accessing display properties.
 *
 * Links to a main image and its overlays with associated look up
 * tables and display settings can be stored in an OverlayList. It has responsibility
 * for maintaining the display order and settings. Display widgets, deciding what to
 * display and how, should iterate over an OverlayList objects contents and
 * behave accordingly.
 */
class OverlayList :  public ImageGroupObserver
{
public:
  typedef boost::shared_ptr< OverlayList > Handle;
  typedef boost::weak_ptr< OverlayList > WeakHandle;
  
  Handle clone();
  static Handle create(ImageGroup::Handle i);

  ~OverlayList();
  const MetaImage::Handle getMainMetaImage() const;
  const MetaImage::Handle getMetaImage(Image::Handle i) const;
  const Image::Handle     getMainImage() const; 
  const ImageGroup::Handle getImageGroup() const { return m_imageGroup; }

  const MetaImage::Handle getActiveMetaImage() const;
  void setActiveMetaImage(const MetaImage::Handle mi);
  const MetaImage::Handle getLabelMetaImage() const;
  void setLabelMetaImage(const MetaImage::Handle mi);

  void setTransparency(float trans);
  void setLookUpTable(LookUpTable::Handle lut);
  void setSecondaryLookUpTable(LookUpTable::Handle lut);
  void setModTransparency(float);
  void setModImage(Image::Handle img);
  void setVisibility(bool state);
  void setReadOnly(bool state);
  void moveOverlayUp();
  void moveOverlayDown();

  const MetaImageListIt begin();
  const MetaImageListIt end();

  void update(const ImageGroup* i);
  
  void attach(OverlayListObserver* o);
  void detach(OverlayListObserver* o);
  void notify(OverlayListMsg message);
  
  Image::Handle inqActiveImage();
  inline int inqX();
  inline int inqY();
  inline int inqZ();

private:  
  
  OverlayList(ImageGroup::Handle i);
  void loadOverlaysList();
  void assignLatestLUT();
  void setOverlays(std::list<MetaImage::Handle>);
  
  std::list<MetaImage::Handle> m_overlays;
  ImageGroup::Handle   m_imageGroup;    
  int m_xDim;
  int m_yDim; 
  int m_zDim;
  int m_currentLut;
  MetaImage::Handle    m_activeMetaImage;
  MetaImage::Handle    m_labelMetaImage;
  std::list<OverlayListObserver*> m_observers;

};

int OverlayList::inqX(){return m_xDim;}
int OverlayList::inqY(){return m_yDim;}
int OverlayList::inqZ(){return m_zDim;}

class OverlayListObserver
{
 
public:

  virtual ~OverlayListObserver() {}
  virtual void update(const OverlayList* list, OverlayListMsg message) = 0;

  OverlayListObserver() {}
};

inline bool isValidOverlayList(const OverlayList::Handle ol)
{
  if(!ol.get()){return false;}
  else         {return true;}
}

#endif
