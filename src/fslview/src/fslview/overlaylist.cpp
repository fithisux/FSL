/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */


#include "overlaylist.h"
#include <assert.h>
#include <algorithm>
#include <qobject.h>

//#define DEBUGGING
#include "tracker.h"

class CopyOverlay {
public:
  CopyOverlay(){}
  void operator()(MetaImage::Handle mi) 
  {
    MetaImage::Handle newMi = mi->clone();
    m_newList.push_back(newMi);
  }
  std::list<MetaImage::Handle> getOverlays(){return m_newList;}

private:
  std::list<MetaImage::Handle> m_newList;
};

class ImageSearch {
public:
  ImageSearch(Image::Handle i) : m_found(false),m_image(i) {}
  void operator()(MetaImage::Handle mi) 
  {
    if(mi->getImage() ==  m_image){m_found = true;m_metaImage = mi;}
  }
  bool m_found;
  MetaImage::Handle m_metaImage;
private:
  Image::Handle m_image;
};

class AddOverlay {
public:
  AddOverlay(std::list<MetaImage::Handle> ol,
             ImageGroup::Handle ig,int &n) 
    : m_oldlist(ol),
      m_imgGrp(ig),m_curLut(n) {}
  void operator()(Image::Handle i) 
  {
    ImageSearch search = std::for_each(m_oldlist.begin(),m_oldlist.end(),ImageSearch(i));

    if(search.m_found == true)
      {
        m_newlist.push_back(search.m_metaImage);
      }
    else
      {
        LookUpTable::Handle lut; 

        if(!m_newlist.empty())
          {
            lut = m_imgGrp->getLut(i->getInfo()->inqLutName());
            if(lut) lut = m_imgGrp->getLut(m_curLut++);        
            if (m_curLut >= m_imgGrp->getInitialLutCount())m_curLut = 0;
          }

        m_newlist.push_back(MetaImage::create(i,ImageDisplaySetting::create(i, lut, 1.0f, true)));  
      }
  }
  std::list<MetaImage::Handle> m_newlist;
private:
  std::list<MetaImage::Handle> m_oldlist;  
  ImageGroup::Handle m_imgGrp;  
  int &m_curLut;

};

class MainSearch 
{
public:
  MainSearch(){}
  void operator()(MetaImage::Handle mi) 
  {
    if(mi->getInfo()->isMainImage())
      m_mainMetaImage = mi;
  }
  MetaImage::Handle m_mainMetaImage;
};

OverlayList::OverlayList(ImageGroup::Handle i):m_imageGroup(i),m_currentLut(0)
{
  m_imageGroup->attach(this);
  
  loadOverlaysList();
}

OverlayList::~OverlayList(){m_imageGroup->detach(this);}

const MetaImage::Handle OverlayList::getMainMetaImage() const
{
  MainSearch search = std::for_each(m_overlays.begin(),
                                    m_overlays.end(),
                                    MainSearch());
  return search.m_mainMetaImage;
}

const MetaImage::Handle OverlayList::getMetaImage(Image::Handle i) const
{
  ImageSearch search = std::for_each(m_overlays.begin(),
				     m_overlays.end(),
				     ImageSearch(i));
  return search.m_metaImage;
}

OverlayList::Handle OverlayList::create(ImageGroup::Handle i)
{
  return OverlayList::Handle(new OverlayList(i));
}


const Image::Handle OverlayList::getMainImage() const 
{  
  MainSearch search = std::for_each(m_overlays.begin(),
                                    m_overlays.end(),
                                    MainSearch());
  return search.m_mainMetaImage->getImage();
}
  
void OverlayList::loadOverlaysList()
{ 
  AddOverlay add = std::for_each(m_imageGroup->begin(),m_imageGroup->end(),
				 AddOverlay(m_overlays,m_imageGroup,m_currentLut));
  m_overlays.clear();
  m_overlays = add.m_newlist;

  setActiveMetaImage(m_overlays.back());
  
  m_xDim = m_overlays.front()->inqX();
  m_yDim = m_overlays.front()->inqY();
  m_zDim = m_overlays.front()->inqZ();
}

struct ImageNotIn {
  ImageNotIn(ImageGroup::Handle ig): m_imageGroup(ig) {}
  
  bool operator()(const MetaImage::Handle mi) const
  {
	Image::Handle i = mi->getImage();
    ImageGroup::ImageList::iterator it = std::find(m_imageGroup->begin(), m_imageGroup->end(), i);

    return (it == m_imageGroup->end());
  }

  ImageGroup::Handle m_imageGroup;
};

const MetaImageListIt OverlayList::begin()
{
  return m_overlays.begin();
}

const MetaImageListIt OverlayList::end()
{
  return m_overlays.end();
}

void OverlayList::attach(OverlayListObserver* o)
{
  m_observers.remove(o);
  m_observers.push_back(o);
}

void OverlayList::detach(OverlayListObserver* o)
{
  m_observers.remove(o);
}

struct Update
{
  Update(OverlayList* l, OverlayListMsg m): m_ol(l),m_msg(m) {}

  void operator()(OverlayListObserver* v)
  {
    TRACKER("Update::operator()(OverlayListObserver* v)");
    v->update(m_ol,m_msg);
  }

  OverlayList*   m_ol;
  OverlayListMsg m_msg;
};

void OverlayList::notify(OverlayListMsg message)
{
  TRACKER("OverlayList::notify"); 
  MESSAGE(QString("notifying %1 observers").arg(m_observers.size()));

  std::for_each(m_observers.begin(), m_observers.end(), 
               Update(this,message));
}

//
//! Returns the currently selected @ref MetaImage for this overlay list.
//!
//! Use to find out which overlay is currently selected in the overlay
//! list displayed in @ref ImageWindow instances.
//!
//! @return Handle to the MetaImage
//
const MetaImage::Handle OverlayList::getActiveMetaImage() const
{
  return m_activeMetaImage;
}

//
//! Returns the currently disignated label @ref MetaImage for this overlay list.
//!
//! @return Handle to the label MetaImage
//
const MetaImage::Handle OverlayList::getLabelMetaImage() const
{
  return m_labelMetaImage;
}

void OverlayList::assignLatestLUT()
{
  TRACKER("OverlayList::assignLatestLUT()");
  MetaImage::Handle mi = getActiveMetaImage();
  
  if(m_activeMetaImage)
    {
      m_activeMetaImage->getDs()->setLookUpTable(m_imageGroup->getLatestLUT());
      notify(OverlayListMsg(LookUpTable));  
    }
}


void OverlayList::update(const ImageGroup* i)
{
  TRACKER("OverlayList::update(ImageGroup* i, ImageGroupMsg message");

  switch(i->inqMessage())
  {
  case ImageGroup::NewLookUpTable:  assignLatestLUT(); notify(OverlayListMsg(LookUpTable)); break;
  case ImageGroup::NewOverlay:    
    {
      Image::Handle newImage = m_imageGroup->getLatestImage();

      LookUpTable::Handle lut;

      if(!newImage->getInfo()->inqLutName().empty())
	lut = m_imageGroup->getLut(newImage->getInfo()->inqLutName());
      if(!lut) {
	if (newImage->getInfo()->inqImageName().find("stat") != std::string::npos)
	  lut = m_imageGroup->getNextLut();
	else if (newImage->getInfo()->inqImageName().find("mask") != std::string::npos)
	  lut = m_imageGroup->getNextLut();
	else
	  lut = LookUpTable::greyScale();
      }
      
      if (m_currentLut >= i->getInitialLutCount()) m_currentLut = 0;
      m_overlays.push_back(MetaImage::create(newImage,ImageDisplaySetting::create(newImage, lut, 1.0f, true)));
      setActiveMetaImage(m_overlays.back());
      notify(OverlayListMsg(Add));
    }
    break;
  case ImageGroup::RemOverlay:
    m_overlays.erase(std::remove_if(m_overlays.begin(),m_overlays.end(),ImageNotIn(m_imageGroup)),
					 m_overlays.end());
    setActiveMetaImage(m_overlays.back());
    notify(OverlayListMsg(Rem));
    break;
  case ImageGroup::Lock:          notify(OverlayListMsg(Security));  break;
  case ImageGroup::NameChange:    notify(OverlayListMsg(ImageName)); break;
  case ImageGroup::None: break;
  }
}

void OverlayList::setActiveMetaImage(const MetaImage::Handle mi)
{    
  TRACKER("OverlayList::setActiveMetaImage");
  m_activeMetaImage = mi;
  notify(OverlayListMsg(Select));
}

void OverlayList::setLabelMetaImage(const MetaImage::Handle mi)
{    
  TRACKER("OverlayList::setLabelMetaImage");
  m_labelMetaImage = mi;
  notify(OverlayListMsg(Select));
}

void OverlayList::setTransparency(float trans)
{
  TRACKER("OverlayList::setTransparency");
  if(m_activeMetaImage)
    {
      m_activeMetaImage->getDs()->setTransparency(trans);
      notify(OverlayListMsg(Transparency));
    }
}

void OverlayList::setModTransparency(float trans)
{
  if(m_activeMetaImage)
    {
      m_activeMetaImage->getDs()->setModTransparency(trans);
      notify(OverlayListMsg(ModImage));
    }
}

void OverlayList::setVisibility(bool state)
{
  if(m_activeMetaImage)
    {    
      m_activeMetaImage->setVisibility(state);
      notify(OverlayListMsg(Visibility));
    }
}  

void OverlayList::setReadOnly(bool state)
{
  TRACKER("OverlayList::setReadOnly(bool state)");
  if(m_activeMetaImage)
    {    
      m_activeMetaImage->getImage()->getInfo()->setReadOnly(state);
      notify(OverlayListMsg(Security));
      m_imageGroup->notify(ImageGroup::Lock);
    }
}

void OverlayList::moveOverlayUp()
{
  TRACKER("OverlayList::moveOverlayUp");

  if(m_activeMetaImage)   
  {
     MetaImageListIt cur = std::find(m_overlays.begin(),
                                     m_overlays.end(),
                                     m_activeMetaImage);
     assert(cur != m_overlays.end());
  
     MetaImageListIt next(cur);
          
     if(++next != m_overlays.end())
       {
         std::swap(*cur,*next);              
         notify(OverlayListMsg(Order));
       }
    }
}

void OverlayList::moveOverlayDown()
{
  TRACKER("OverlayList::moveOverlayDown");

  if(m_activeMetaImage)
  {
     MetaImageListIt cur = std::find(m_overlays.begin(),
                                     m_overlays.end(),
                                     m_activeMetaImage);
     assert(cur != m_overlays.end());
        
     MetaImageListIt prev(cur);
          
     if(cur != m_overlays.begin())
     {       
       --prev;
       std::swap(*cur,*prev);              
       notify(OverlayListMsg(Order));
     }
   }
}

void OverlayList::setLookUpTable(LookUpTable::Handle lut)
{
  TRACKER("OverlayList::setLookUpTable(LookUpTable::Handle lut)");
  if(m_activeMetaImage)
  {
    m_activeMetaImage->getDs()->setLookUpTable(lut);
    notify(OverlayListMsg(LookUpTable));
  }
}

void OverlayList::setSecondaryLookUpTable(LookUpTable::Handle lut)
{
  TRACKER("OverlayList::setSecondaryLookUpTable(LookUpTable::Handle lut)");
  if(m_activeMetaImage)
  {
    m_activeMetaImage->getDs()->setSecondaryLookUpTable(lut);
    BriCon::Handle bc(m_activeMetaImage->getDs()->inqBriCon());
    bc->setRange(bc->inqMax()/10, bc->inqMax());
    notify(OverlayListMsg(LookUpTable));
  }
}

void OverlayList::setModImage(Image::Handle img)
{
  if(m_activeMetaImage)
  {
     ImageSearch search = std::for_each(m_overlays.begin(),
                                        m_overlays.end(),
                                        ImageSearch(img));

     if(search.m_found == true)
      {
        search.m_metaImage->setVisibility(false);      
      }

    m_activeMetaImage->getDs()->setModImage(img);
    notify(OverlayListMsg(ModImage));
  }
}

Image::Handle OverlayList::inqActiveImage()
{      
  Image::Handle img;

  if(m_activeMetaImage)
      img = m_activeMetaImage->getImage();

  return img;
}

OverlayList::Handle OverlayList::clone()
{
  OverlayList::Handle hnd =  Handle(new OverlayList(m_imageGroup));
  hnd->setOverlays(m_overlays);
  return hnd;
}

void OverlayList::setOverlays(std::list<MetaImage::Handle> list)
{    
  CopyOverlay copyOverlay = std::for_each(list.begin(),list.end(),CopyOverlay());
  m_overlays = copyOverlay.getOverlays(); 
 
 if(!m_overlays.empty()){setActiveMetaImage(m_overlays.back());}
}

