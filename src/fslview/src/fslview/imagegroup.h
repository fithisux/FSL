/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(IMAGEGROUP_H)
#define IMAGEGROUP_H

#include <boost/shared_ptr.hpp>
#include <memory>
#include "storage/image.h"
#include "imagedisplaysetting.h"
#include "lookuptable.h"

class ImageGroup;

class ImageGroupObserver
{

public:
  virtual ~ImageGroupObserver() {}

  virtual void update(const ImageGroup* i) = 0;

  ImageGroupObserver() {}
};

/**
 * @author David Flitney <flitney@fmrib.ox.ac.uk>
 * @author James Saunders <jim@fmrib.ox.ac.uk>
 *
 * @date   Mon Dec 23 17:25:52 2002
 * 
 * @brief ImageGroup groups a collection of images in a common
 * object. 
 *
 * A main image and several overlays with associated look up
 * tables can be stored in an ImageGroup object. It has responsibility
 * for maintaining the links between compatible and associated images.
 */
class ImageGroup
{
public:
  typedef std::vector<Image::Handle> ImageList;
  typedef std::vector<LookUpTable::Handle> LutList;
  typedef boost::shared_ptr< ImageGroup > Handle;
  typedef boost::weak_ptr< ImageGroup > WeakHandle;
  typedef enum {None, NewLookUpTable, NewOverlay, RemOverlay, Lock, NameChange} Msg;
 
  static Handle create(Image::Handle image);
  static Handle createNullImage();

  bool addOverlay(Image::Handle image);
  bool addUniqueOverlay(Image::Handle image);
  void addLookUpTable(LookUpTable::Handle lut);
  bool remOverlay(Image::Handle image);
  ImageList::iterator begin();
  ImageList::iterator end();
  ImageList::size_type size(); 

  Image::Handle getImage(int n) const;
  LutList::iterator beginLutList();
  LutList::iterator endLutList();
  LookUpTable::Handle getLut(int n) const;
  LookUpTable::Handle getNextLut() const;
  LookUpTable::Handle getLut(std::string const &name) const;

  int  getInitialLutCount() const;
  void attach(ImageGroupObserver* o);
  void detach(ImageGroupObserver* o);
  void notify(ImageGroup::Msg message);
  ImageGroup::Msg inqMessage() const;

  Image::Handle getMainImage();
  Image::Handle getLatestImage() const;
  LookUpTable::Handle getLatestLUT();
  bool inqTarnished();

  int inqX();
  int inqY();
  int inqZ();
  float inqMin();
  float inqMax();
  virtual ~ImageGroup();

private:
  ImageGroup(Image::Handle image);

  struct Implementation;
  const std::auto_ptr<Implementation> m_impl;
};


#endif
