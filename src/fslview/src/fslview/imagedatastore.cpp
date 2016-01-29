/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "imagedatastore.h"
#include <vector>
#include <algorithm>

struct ImageDataStore::Implementation
{
  Implementation():m_listPos(0){};
  std::vector<ImageData::Handle> m_imgDataList;
  unsigned int m_listPos;
};

ImageData::Handle   createImageData(MetaImage::Handle mi)
{
  ColorRGBAHandle rgba;
  ImageData::Handle id = ImageData::create(mi,rgba);
  return id;
}

ImageDataStore::ImageDataStore(OverlayList::Handle ol):
  m_impl(new Implementation())
{

  m_impl->m_imgDataList.clear();

  std::transform(ol->begin(),
                 ol->end(),
                 std::back_inserter(m_impl->m_imgDataList),
                 createImageData);  
}

ImageDataStore::~ImageDataStore(){}

ImageDataStore::Handle ImageDataStore::create(OverlayList::Handle ol)
{
  Handle dst(new ImageDataStore(ol));
  return dst;
}  


void ImageDataStore::resetPos()
{
  m_impl->m_listPos = 0;
}

bool ImageDataStore::currentEmpty()
{
  bool result(false);

  if(m_impl->m_listPos >= m_impl->m_imgDataList.size())
    result = true;
  else if (m_impl->m_imgDataList[m_impl->m_listPos]->getBuffer())
    result = false;
  
  return result;
}

ImageData::Handle ImageDataStore::current()
{
  return m_impl->m_imgDataList.at(m_impl->m_listPos);
}

void ImageDataStore::next()
{
  m_impl->m_listPos++;
}

class ImageDataSearch {
public:
  ImageDataSearch() : m_dtiLinesFound(false) {}
  void operator()(ImageData::Handle i) 
  {
    if(i->inqDtiDisplay()== DtiDisplay(Lines))
      {
        m_dtiLinesFound = true;
        m_linesImageData = i;
      }
  }
  bool m_dtiLinesFound;
  ImageData::Handle m_linesImageData;
};

bool ImageDataStore::isDtiLineOverlay()
{
    ImageDataSearch search = std::for_each(m_impl->m_imgDataList.begin(),
                                           m_impl->m_imgDataList.end(),
                                           ImageDataSearch());


    return search.m_dtiLinesFound;
}

ImageData::Handle ImageDataStore::getDtiLineOverlay()
{
    ImageDataSearch search = std::for_each(m_impl->m_imgDataList.begin(),
                                           m_impl->m_imgDataList.end(),
                                           ImageDataSearch());

    return search.m_linesImageData;
}
