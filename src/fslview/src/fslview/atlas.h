/*  FSLView - 2D/3D Interactive Image Viewer

    David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(ATLAS_H)
#define ATLAS_H

#include <boost/shared_ptr.hpp>
#include "storage/image.h"

#include "cursor.h"

class AtlasOptions
{
public:
  AtlasOptions(int structure=0, bool superimpose=false, bool locate=false);

  int structure() { return m_structure; }
  void structure(int i) { m_structure = i; }
  
  std::string structureName() { return m_structureName; }
  void structureName(std::string name) { m_structureName = name; }

  bool locate() { return m_locate; }
  void locate(bool y) { m_locate = y; }

  bool superimpose() { return m_superimpose; }
  void superimpose(bool y) { m_superimpose = y; }

private:
  int m_structure;
  std::string m_structureName;
  bool m_locate;
  bool m_superimpose;
};

class Atlas
{
public:
  typedef boost::shared_ptr< Atlas > Handle;
  typedef std::map<int, std::string>::const_iterator ConstLabelIterator;
  typedef std::vector<Image::Handle> ImageStore;

  enum Type { Unknown, Label, Probabilistic };

  void addLabel(int, const std::string&);
  void addCentre(int, short, short, short, short);

  std::string inqName() const;
  ImageStore& inqImages() const;
  ImageStore& inqSummaryImages() const;
  Image::Handle inqCurrentImage() const;
  Image::Handle inqCurrentSummaryImage() const;
  virtual Type inqType() const { return Atlas::Unknown; }
  std::string inqStructureNameByIndex(unsigned int) const;

  virtual unsigned int inqNumLabels() const = 0;

  void selectCompatibleImages(const Image::Handle&);

  virtual std::string getDescription(float, float, float) const = 0;
  virtual unsigned int getProbability(unsigned int, float, float, float) const = 0;
  virtual float getAverageProbability(Image::Handle, unsigned int) const = 0;
  
  Cursor::Handle getCursor(const Image::Handle&, int) const;

  ConstLabelIterator begin();
  ConstLabelIterator end();

  static Handle create(const ImageStore&, const ImageStore&, const std::string&);
  virtual ~Atlas();
 
protected:
  Atlas(const ImageStore&, const ImageStore&, const std::string&);

  struct Implementation;  
  const std::auto_ptr<Implementation> m_impl;
};

class ProbabilisticAtlas: public Atlas
{
public:
  virtual unsigned int inqNumLabels() const;

  virtual std::string getDescription(float, float, float) const;
  virtual unsigned int getProbability(unsigned int, float, float, float) const;
  virtual float getAverageProbability(Image::Handle, unsigned int) const;

  virtual Type inqType() const { return Atlas::Probabilistic; }

  static Handle create(const ImageStore&, const ImageStore&, const std::string&);
private:
  ProbabilisticAtlas(const ImageStore&, const ImageStore&, const std::string&);
};

class LabelAtlas: public Atlas
{
public:
  virtual unsigned int inqNumLabels() const;

  virtual std::string getDescription(float, float, float) const;
  virtual unsigned int getProbability(unsigned int, float, float, float) const;
  virtual float getAverageProbability(Image::Handle, unsigned int) const;

  virtual Type inqType() const { return Atlas::Label; }

  static Handle create(const ImageStore&, const ImageStore&, const std::string&);
private:
  LabelAtlas(const ImageStore&, const ImageStore&, const std::string&);
};

class AtlasGroup
{
public:
  typedef boost::shared_ptr<AtlasGroup> Handle;
  typedef std::exception Exception;

  typedef std::pair<std::string, Atlas::Handle> AtlasItem;
  typedef std::map<std::string, Atlas::Handle> AtlasContainer;
  typedef AtlasContainer::const_iterator ConstIterator;

  void readAtlas(const std::string&, const std::string&);

  Atlas::Handle getAtlasByName(const std::string&);
  void selectCompatibleAtlases(const Image::Handle&);

  static AtlasGroup::Handle create();

  ConstIterator begin() const { return m_atlases.begin(); }
  ConstIterator end() const { return m_atlases.end(); }

  virtual ~AtlasGroup() {}
private:
  AtlasGroup();

  static AtlasContainer m_atlases;
};

#endif
