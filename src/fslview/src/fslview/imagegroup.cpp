/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */


#include "imagegroup.h"
#include "imagedisplaysetting.h"
#include "tracker.h"
#include <qobject.h>
#include <qmessagebox.h>
#include <qapplication.h>

#include <list>

using namespace std;

class CheckImageTarnish{
public:
  CheckImageTarnish():m_tarnished(false){}
  void operator()(Image::Handle image) 
  {
    m_tarnished = image->getInfo()->inqTarnished();
  }
  bool result(){return m_tarnished;}
private:
  bool m_tarnished;
};

class LutSearch {
public:
  LutSearch(std::string const & name) : m_name(name) {}
  void operator()(LookUpTable::Handle lut)
  {
    if(lut->inqLutName() ==  m_name){m_lut = lut;}
  }

  LookUpTable::Handle result(){return m_lut;}
private:
  std::string m_name;  
  LookUpTable::Handle m_lut;
};

struct ImageGroup::Implementation
{
  Implementation(Image::Handle image)
  {
    image->getInfo()->setAsMainImage();
    m_imageList.push_back(image);
    m_messages.push_back(None);

    m_xDim = m_imageList.front()->getVolume(0)->inqX();
    m_yDim = m_imageList.front()->getVolume(0)->inqY();
    m_zDim = m_imageList.front()->getVolume(0)->inqZ();
    m_min  = m_imageList.front()->getInfo()->inqMin();
    m_max  = m_imageList.front()->getInfo()->inqMax();
  }

  void addLookUpTable(LookUpTable::Handle l)
  {
    m_lookUpTables.push_back(l);
  }

   ~Implementation(){TRACKER("ImageGroup::~Implementation()");}

  std::list<ImageGroupObserver* >  m_observers;
  ImageList m_imageList;
  LutList m_lookUpTables;
  int  m_currentLut;
  int  m_initLutSize;
  int  m_xDim;
  int  m_yDim;
  int  m_zDim;
  float  m_min;
  float  m_max;
  
  std::list<ImageGroup::Msg> m_messages;
};

//! @brief Class constructor
//!
//! @param image an initial image
//!
//! Creates a new ImageGroup consisting of image and registers the basic look up tables.
ImageGroup::ImageGroup(Image::Handle image):
  m_impl(new Implementation(image))
{
  TRACKER("ImageGroup::ImageGroup(Image::Handle& image)");

  addLookUpTable(LookUpTable::greyScale());
  addLookUpTable(LookUpTable::redYellow());
  addLookUpTable(LookUpTable::blueLightblue());
  addLookUpTable(LookUpTable::red());
  addLookUpTable(LookUpTable::blue());
  addLookUpTable(LookUpTable::green());
  addLookUpTable(LookUpTable::yellow());
  addLookUpTable(LookUpTable::pink());
  addLookUpTable(LookUpTable::hot());
  addLookUpTable(LookUpTable::cool());
  addLookUpTable(LookUpTable::copper());
  addLookUpTable(LookUpTable::render1());
  addLookUpTable(LookUpTable::render1t());
  addLookUpTable(LookUpTable::render2());
  addLookUpTable(LookUpTable::render2t());
  addLookUpTable(LookUpTable::render3());
  addLookUpTable(LookUpTable::cortical());
  addLookUpTable(LookUpTable::subcortical());
  addLookUpTable(LookUpTable::rainbow());

  m_impl->m_initLutSize = m_impl->m_lookUpTables.size();
  m_impl->m_currentLut = 0;
}

ImageGroup::~ImageGroup()
{
  TRACKER("ImageGroup::~ImageGroup()");
}

/** 
 * @brief Generate a new image group with a given image. 
 * 
 * @param image Main image for new image group
 * 
 * @return Reference counted handle to a new ImageGroup.
 *
 * The create methods are in-lieu of conventional constructors. See
 * @ref ImageGroup::ImageGroup for documentation of the constructor
 * behaviour
 */
ImageGroup::Handle ImageGroup::create(Image::Handle image)
{
  return ImageGroup::Handle(new ImageGroup(image));
}

/** 
 * @brief Add an overlay to this image group
 * 
 * @param image New image to be added to the image groups overlay list.
 * 
 * @return true if operation succeeds.
 */
bool ImageGroup::addOverlay(Image::Handle image)
{
  m_impl->m_imageList.push_back(image);

  ImageInfo::Handle info(image->getInfo());

  if(info->inqLutName().empty()) {
    if(info->isStatImage() || info->isMaskImage()) {
      m_impl->m_currentLut = (m_impl->m_currentLut + 1) % m_impl->m_lookUpTables.size();

      while(!m_impl->m_lookUpTables[m_impl->m_currentLut]->isAutoSelectable())
	m_impl->m_currentLut = (m_impl->m_currentLut + 1) % m_impl->m_lookUpTables.size();
    }
  }

  if(info->inqPurpose() == ImageIntent::Unknown) {
    if(info->isStatImage()) info->setPurpose(ImageIntent::Statistic);
    else if(info->isMaskImage()) info->setPurpose(ImageIntent::Label);
  }

  notify(NewOverlay);
  return true;
}

/** 
 * @brief Add an overlay to this image group iff it isn't already there.
 * 
 * @param image New image to be added to the image groups overlay list.
 * 
 * @return true if operation succeeds.
 */
bool ImageGroup::addUniqueOverlay(Image::Handle image)
{
  bool result(false);

  ImageList::iterator cur = std::find(m_impl->m_imageList.begin(),
				      m_impl->m_imageList.end(),
				      image);
  if(cur == m_impl->m_imageList.end())
    result = addOverlay(image);

  return result;
}

/** 
 * @brief Add a lut to the image groups look up table list.
 * 
 * @param lut New lut to be added to the lut list.
 */
void ImageGroup::addLookUpTable(LookUpTable::Handle lut)
{
  m_impl->addLookUpTable(lut);
  notify(NewLookUpTable);
}

//! @brief Get a handle to the latest Image added to this ImageGroup
//!
//! @return A handle referencing the latest(last) Image in this ImageGroup
Image::Handle ImageGroup::getLatestImage() const
{
  return m_impl->m_imageList.back();
}

//! @brief Get a handle to the latest LookUpTable added to this ImageGroup
//!
//! @return A handle referencing the latest(last) LUT in this ImageGroup
LookUpTable::Handle ImageGroup::getLatestLUT()
{
  LookUpTable::Handle latest;
  
  if(!m_impl->m_lookUpTables.empty())
  {
    latest = m_impl->m_lookUpTables.back();
  }

  return latest;
}

Image::Handle ImageGroup::getImage(int n) const
{
  Image::Handle img;
  int size = m_impl->m_imageList.size();

  if((size > 0) && (n < size))
    img = m_impl->m_imageList[n];

  return img;
}

//! @brief Get the next auto-selectable LUT
//! @return A handle referencing the next auto-selectable LookUpTable
//!
//! Use this method to select a LUT deemed suitable for stats/mask images. 
LookUpTable::Handle ImageGroup::getNextLut() const
{
  return m_impl->m_lookUpTables[m_impl->m_currentLut];
}

//! @brief Get the nth LUT
//! @return A handle referencing the nth LUT stored in this ImageGroup 
LookUpTable::Handle ImageGroup::getLut(int n) const
{
  LookUpTable::Handle lut;
  int size = m_impl->m_lookUpTables.size();

  if((size > 0) && (n < size))
    lut = m_impl->m_lookUpTables[n];

  return lut;
}

LookUpTable::Handle ImageGroup::getLut(std::string const &name) const
{
  LookUpTable::Handle r = LookUpTable::greyScale();

  try {
    if((0 == name.compare(0, 4, "none")) || (name.empty())) {
      r = LookUpTable::greyScale();
    } else {
      LutSearch search = std::for_each(m_impl->m_lookUpTables.begin(),
				     m_impl->m_lookUpTables.end(),
				     LutSearch(name));

      r = search.result();
    }

    if(!r)
      r = LookUpTable::load(name);
    
  } catch (std::exception& e) {
//     QMessageBox::warning(NULL, "ImageGroup::getLut: Exception while trying to open LUT!", 
// 			 QString("Error message was: \"%1\"! \nReverting to Greyscale colormap.").arg(e.what()));
    r = LookUpTable::greyScale();
  }
  return r;
}

int ImageGroup::getInitialLutCount() const
{
  return m_impl->m_initLutSize;
}

/** 
 * @brief Remove a given overlay from the image group.
 * 
 * @param image The image handle for the overlay image to be removed.
 * 
 * @return true if operation succeeded.
 */
bool ImageGroup::remOverlay(Image::Handle image)
{
  bool result(false);

  if (!m_impl->m_imageList.empty())
    {     

      ImageList::iterator cur = std::find(m_impl->m_imageList.begin(),
                                          m_impl->m_imageList.end(),
                                          image);
      if(cur != m_impl->m_imageList.end())
        {
          m_impl->m_imageList.erase(cur);
          notify(RemOverlay);
          result = true;
        }
    }

  return result;
}

/** 
 * @brief Access the groups main image.
 * 
 * @return Image::Handle for the groups main image.
 */
Image::Handle ImageGroup::getMainImage()
{
  return m_impl->m_imageList.front();
}

ImageGroup::ImageList::iterator  ImageGroup::begin()
{
  return m_impl->m_imageList.begin();
}

ImageGroup::ImageList::iterator  ImageGroup::end()
{
  return m_impl->m_imageList.end();
}

ImageGroup::ImageList::size_type ImageGroup::size()
{
  return m_impl->m_imageList.size();
}

ImageGroup::LutList::iterator ImageGroup::beginLutList()
{
  return m_impl->m_lookUpTables.begin();
}

ImageGroup::LutList::iterator ImageGroup::endLutList()
{
  return m_impl->m_lookUpTables.end();
}

int ImageGroup::inqX()
{
  return m_impl->m_xDim;
}

int ImageGroup::inqY()
{ 
   return m_impl->m_yDim; 
}

int ImageGroup::inqZ()
{   
  return m_impl->m_zDim;   
}

float ImageGroup::inqMin() { return m_impl->m_min; }

float ImageGroup::inqMax() { return m_impl->m_max; }

void ImageGroup::attach(ImageGroupObserver* o)
{
   m_impl->m_observers.remove(o);
   m_impl->m_observers.push_back(o);
}

void ImageGroup::detach(ImageGroupObserver* o)
{
   m_impl->m_observers.remove(o);
}

ImageGroup::Msg ImageGroup::inqMessage() const { return m_impl->m_messages.back(); }

struct Update
{
  Update(ImageGroup* i): m_imageGroup(i) {}

  void operator()(ImageGroupObserver* v)
  {
    v->update(m_imageGroup);
  }

  ImageGroup*   m_imageGroup;
};

void ImageGroup::notify(ImageGroup::Msg message)
{
  m_impl->m_messages.push_back(message);

  std::for_each(m_impl->m_observers.begin(), m_impl->m_observers.end(), 
               Update(this));

  m_impl->m_messages.pop_back();
}

bool ImageGroup::inqTarnished()
{
  CheckImageTarnish check = std::for_each(m_impl->m_imageList.begin(),
                                          m_impl->m_imageList.end(),
                                          CheckImageTarnish());
  return check.result();
}

