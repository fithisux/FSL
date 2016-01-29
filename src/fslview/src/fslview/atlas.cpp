/*  FSLView - 2D/3D Interactive Image Viewer

    David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "atlas.h"
#include "filemanager.h"
#include "preferences.h"

#include <sstream>
#include <map>

using namespace std;

AtlasOptions::AtlasOptions(int structure, bool superimpose, bool locate):
  m_structure(structure), m_locate(locate), m_superimpose(superimpose)
{
}

struct Atlas::Implementation
{
  Implementation(const Atlas::ImageStore& i, const Atlas::ImageStore& s, const string& name):
    m_images(i), m_summaries(s), m_name(name) {}

  typedef map<int, string> Labels;
  typedef Labels::const_iterator LabelsConstIterator;
  typedef map<int, Cursor::Handle> Cursors;
  typedef Cursors::const_iterator CursorsConstIterator;
  typedef map<int, string> References;
  typedef References::const_iterator ReferencesConstIterator;

  Volume::Handle getVolume(short n) { return m_images.at(0)->getVolume(n); }
  LabelsConstIterator findLabel(int n) { return m_labels.find(n); }
  CursorsConstIterator findCursor(int i) { return m_cogs.find(i); }

  Atlas::ImageStore m_images;
  Atlas::ImageStore m_summaries;
  Image::Handle m_image;
  Image::Handle m_summary;
  Labels m_labels;
  Cursors m_cogs;
  References m_references;
  string m_name;
  Atlas::Type m_type;
};

Atlas::Handle ProbabilisticAtlas::create(const ImageStore& i, const ImageStore& s, const string& n)
{
  Handle atlas(new ProbabilisticAtlas(i, s, n));
  return atlas;
}

Atlas::Handle LabelAtlas::create(const ImageStore& i, const ImageStore& s, const string& n)
{
  Handle atlas(new LabelAtlas(i, s, n));
  return atlas;
}

Atlas::Atlas(const ImageStore& i, const ImageStore& s, const string& n):
  m_impl(new Implementation(i, s, n))
{
  m_impl->m_image =   i.at(0);  
  m_impl->m_summary = s.at(0);
}

ProbabilisticAtlas::ProbabilisticAtlas(const ImageStore& i, const ImageStore& s, const string& n):
  Atlas(i, s, n)
{
  m_impl->m_type=Probabilistic;
}

LabelAtlas::LabelAtlas(const ImageStore& i, const ImageStore& s, const string& n):
  Atlas(i, s, n)
{
  m_impl->m_type=Label;
}

Atlas::~Atlas() {}

void Atlas::addLabel(int n, const std::string& l)
{
  m_impl->m_labels.insert(make_pair(n, l));
}

void Atlas::addCentre(int n, short x, short y, short z, short v)
{
  Cursor::Handle c = Cursor::create(x+1, y+1, z+1, v+1);
  c->setCursor(x, y, z, v);
  m_impl->m_cogs.insert( make_pair(n, c) );
}

string Atlas::inqName() const
{
  return m_impl->m_name;
}

Image::Handle Atlas::inqCurrentImage() const
{
  return m_impl->m_image;
}

Image::Handle Atlas::inqCurrentSummaryImage() const
{
  return m_impl->m_summary;
}

Atlas::ImageStore& Atlas::inqImages() const
{
  return m_impl->m_images;
}

Atlas::ImageStore& Atlas::inqSummaryImages() const
{
  return m_impl->m_summaries;
}

Atlas::ConstLabelIterator Atlas::begin()
{
  return m_impl->m_labels.begin();
}

Atlas::ConstLabelIterator Atlas::end()
{
  return m_impl->m_labels.end();
}

Cursor::Handle Atlas::getCursor(const Image::Handle& image, int structureIndex) const
{
  ImageInfo::Handle atlasInfo(m_impl->m_images.at(0)->getInfo());
  ImageInfo::Handle imageInfo(image->getInfo());
  short tx(10), ty(10), tz(10), v(1);

  Implementation::CursorsConstIterator cit = m_impl->findCursor(structureIndex);
  if( cit != m_impl->m_cogs.end() ) {
    Cursor::Handle h = cit->second;
    float x, y, z;
    atlasInfo->voxToMMCoord(h->inqX(), h->inqY(), h->inqZ(), x, y, z);
    imageInfo->mmToVoxCoord(x, y, z, tx, ty, tz);
    v = h->inqV();
  }

  Cursor::Handle c = Cursor::create(tx+1, ty+1, tz+1, v+1);
  c->setCursor(tx, ty, tz, v);
  return c;
}

string LabelAtlas::getDescription(float x, float y, float z) const
{
  ostringstream text;

  ImageInfo::Handle info(m_impl->m_image->getInfo());
  short tx(0), ty(0), tz(0);

  info->mmToVoxCoord(x, y, z, tx, ty, tz);

  Volume::Handle vol(m_impl->m_image->getVolume(0));
  int index( info->isValidCoordinate(tx, ty, tz) ? 
	     int(vol->value(tx, ty, tz)) : 0 );

  text << "<b>" << m_impl->m_name << "</b><br>";

  map<int, string>::const_iterator pos = m_impl->m_labels.find(index);
  if(pos != m_impl->m_labels.end())
    text << pos->second;

  return text.str();
}

unsigned int LabelAtlas::getProbability(unsigned int structure, float x, float y, float z) const
{
  ImageInfo::Handle info(m_impl->m_image->getInfo());
  short tx(0), ty(0), tz(0);
  info->mmToVoxCoord(x, y, z, tx, ty, tz);
  Volume::Handle vol(m_impl->m_image->getVolume(0));
  return info->isValidCoordinate(tx, ty, tz) ? 
    ( (vol->value(tx, ty, tz) == structure) ? 100 : 0 ) : 0;
}

string ProbabilisticAtlas::getDescription(float x, float y, float z) const
{
  ImageInfo::Handle info(m_impl->m_image->getInfo());
  short tx(0), ty(0), tz(0);
  int nvols(info->inqNumVolumes());

  info->mmToVoxCoord(x, y, z, tx, ty, tz);

  multimap<int,string> labels;

  for(int v = 0; v < nvols; ++v) {
    Volume::Handle vol(m_impl->m_image->getVolume(v));
    int prob( info->isValidCoordinate(tx, ty, tz) ? 
	      int(vol->value(tx, ty, tz)) : 0 );
 
    if( prob > 0 ) {
      Implementation::LabelsConstIterator pos = m_impl->findLabel(v);
      if(pos != m_impl->m_labels.end())
	labels.insert(make_pair(prob, pos->second));
    }
  }
  
  unsigned int count(0);

  ostringstream text;
  
  text << "<b>" << m_impl->m_name << "</b><br>";

  for(map<int,string>::reverse_iterator it = labels.rbegin();
      it != labels.rend(); ++it) {
    if(count++)
      text << ", ";
    text << it->first << "% " << it->second;
  }

  if(!count)
    // May want to display nearest object here
    text << "No label found!";

  return text.str();
}

void Atlas::selectCompatibleImages(const Image::Handle& refim)
{
  for(ImageStore::const_iterator it = m_impl->m_images.begin();
      it != m_impl->m_images.end(); ++it) {
    if((*it)->getInfo()->isCompatible(refim->getInfo()))
      m_impl->m_image = *it;
  }
  for(ImageStore::const_iterator it = m_impl->m_summaries.begin();
      it != m_impl->m_summaries.end(); ++it) {
    if((*it)->getInfo()->isCompatible(refim->getInfo()))
      m_impl->m_summary = *it;
  }
}

string Atlas::inqStructureNameByIndex(unsigned int index) const
{
  string name("Unknown structure");
  map<int, string>::const_iterator pos = m_impl->m_labels.find(index);
  if(pos != m_impl->m_labels.end())
    name = pos->second;
  return name;
}

unsigned int ProbabilisticAtlas::getProbability(unsigned int index, float x, float y, float z) const
{
  ImageInfo::Handle info(m_impl->m_image->getInfo());
  short tx(0), ty(0), tz(0);
  info->mmToVoxCoord(x, y, z, tx, ty, tz);
  Volume::Handle vol(m_impl->m_image->getVolume(index));
  return info->isValidCoordinate(tx, ty, tz) ? 
    int(vol->value(tx, ty, tz)) : 0;
}

float ProbabilisticAtlas::getAverageProbability(Image::Handle mask, unsigned int index) const
{
  Volume::Handle m(mask->getVolume(0));
  Volume::Handle p(m_impl->m_image->getVolume(index));

  ImageInfo::Handle maskinfo(mask->getInfo());
  ImageInfo::Handle probinfo(m_impl->m_image->getInfo());

  float total(0); float sum(0);

  for(short z = 0; z < m->inqZ(); ++z)
    for(short y = 0; y < m->inqY(); ++y)
      for(short x = 0; x < m->inqX(); ++x) 
	{
	  float tx(0), ty(0), tz(0);
	  short i, j, k;

	  maskinfo->voxToMMCoord(x, y, z, tx, ty, tz);
	  probinfo->mmToVoxCoord(tx, ty, tz, i, j, k);
	  
	  float weight(m->value(x, y, z));
	  
	  total += p->value(i, j, k) * weight;
	  sum += weight;
	}

  return total / sum;
}

unsigned int ProbabilisticAtlas::inqNumLabels() const
{
  return m_impl->m_labels.size();
}

float LabelAtlas::getAverageProbability(Image::Handle mask, unsigned int index) const
{
  Volume::Handle l(m_impl->m_image->getVolume(0));
  Volume::Handle m(mask->getVolume(0));

  ImageInfo::Handle maskinfo(mask->getInfo());
  ImageInfo::Handle labelinfo(m_impl->m_image->getInfo());

  float total(0); float sum(0);

  for(short z = 0; z < m->inqZ(); ++z)
    for(short y = 0; y < m->inqY(); ++y)
      for(short x = 0; x < m->inqX(); ++x) 
	{
	  float tx(0), ty(0), tz(0);
	  short i, j, k;

	  maskinfo->voxToMMCoord(x, y, z, tx, ty, tz);
	  labelinfo->mmToVoxCoord(tx, ty, tz, i, j, k);
	  
	  float weight(m->value(x, y, z));

	  total += (l->value(i, j, k) == index ? 100 : 0) * weight;
	  sum += weight;
	}

  return total / sum;
}

unsigned int LabelAtlas::inqNumLabels() const
{
  return m_impl->m_labels.size();
}

AtlasGroup::AtlasContainer AtlasGroup::m_atlases;

AtlasGroup::AtlasGroup()
{
  if( m_atlases.size() == 0 ) {

    Preferences p;
    vector<string> atlasdirs(p.inqAtlasPathElements());
 
    for(vector<string>::iterator it = atlasdirs.begin(); 
	it != atlasdirs.end(); ++it){
      if( FileManager::checkFileExists(*it) ) {
	vector<string> fv = FileManager::getFilenames(*it, "*.xml");
	for(vector<string>::size_type i=0; i < fv.size(); i++)
	  readAtlas(*it, fv[i]);
      }
    }
  }
}

AtlasGroup::Handle AtlasGroup::create()
{
  return AtlasGroup::Handle(new AtlasGroup());
}

void AtlasGroup::selectCompatibleAtlases(const Image::Handle& refim)
{
  for(ConstIterator it = m_atlases.begin(); it != m_atlases.end(); ++it)
    it->second->selectCompatibleImages(refim);
}

Atlas::Handle AtlasGroup::getAtlasByName(const std::string& name)
{
  Atlas::Handle atlas;

  std::map<std::string, Atlas::Handle>::const_iterator pos = 
    m_atlases.find(name);

  if(pos != m_atlases.end())
    atlas = pos->second;

  return atlas;
}

void AtlasGroup::readAtlas(const std::string& path, const std::string& filename)
{
  Atlas::Handle atlas(FileManager::readXMLAtlas(path, filename));
  
  m_atlases.insert(make_pair(atlas->inqName(), atlas));
}
