/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Rama Aravind Vorray
		James Saunders
		David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2002-2005 University of Oxford  */

/*  CCOPYRIGHT */

#include "clusterdata.h"
#include <sstream>

//#define DEBUGGING
#include "tracker.h"

struct BaseCluster::Implementation
{
  Implementation(): 
    m_index(0), m_voxels(0),
    m_p(0.0), m_minusLog10P(0.0),
    m_maxZ(0.0), m_maxCOPE(0.0), m_meanCOPE(0.0),
    m_initialised(false)
  {}

  unsigned int m_index, m_voxels;
  float m_p, m_minusLog10P;
  float m_maxZ;
  float m_maxCOPE;
  float m_meanCOPE;
  bool m_initialised;
};

BaseCluster::BaseCluster(const ColumnList& cl):
  m_impl(new Implementation), m_columns(cl)
{
  TRACKER("BaseCluster::BaseCluster(const ColumnList& cl)");
  m_headings.insert(std::make_pair("Cluster Index", Index));
  m_headings.insert(std::make_pair("Voxels", Voxels));
  m_headings.insert(std::make_pair("P", P));
  m_headings.insert(std::make_pair("-log10(P)", MinusLog10P));
  m_headings.insert(std::make_pair("Z-MAX", ZMax));
  m_headings.insert(std::make_pair("COPE-MAX", COPEMax));
  m_headings.insert(std::make_pair("COPE-MEAN", COPEMean));
}

BaseCluster::~BaseCluster() {}

bool BaseCluster::initialised() const { return m_impl->m_initialised; }
void BaseCluster::initialised(bool y) { m_impl->m_initialised = y; }

bool BaseCluster::readColumn(std::istream& is, const std::string& col)
{
  bool flag(true);

  HeadingMap::iterator it = m_headings.find(col);
  if(it != m_headings.end()) {
    switch(it->second) {
    case Index:       is >> m_impl->m_index;       break;
    case Voxels:      is >> m_impl->m_voxels;      break;
    case P:           is >> m_impl->m_p;           break;
    case MinusLog10P: is >> m_impl->m_minusLog10P; break;
    case ZMax:        is >> m_impl->m_maxZ;        break;
    case COPEMax:     is >> m_impl->m_maxCOPE;     break;
    case COPEMean:    is >> m_impl->m_meanCOPE;    break;
    default: flag = false; break;
    }
  } else
    flag = false;

  return flag;
}

std::string BaseCluster::inqIndex() const
{ std::ostringstream s; s << m_impl->m_index; return s.str(); }
std::string BaseCluster::inqSize() const
{ std::ostringstream s; s << m_impl->m_voxels; return s.str(); }
std::string BaseCluster::inqP() const
{ std::ostringstream s; s << m_impl->m_p; return s.str(); }
std::string BaseCluster::inqMinusLog10P() const
{ std::ostringstream s; s << m_impl->m_minusLog10P; return s.str(); }
std::string BaseCluster::inqMaxZ() const
{ std::ostringstream s; s << m_impl->m_maxZ; return s.str(); }
std::string BaseCluster::inqMaxCOPE() const
{ std::ostringstream s; s << m_impl->m_maxCOPE; return s.str(); }
std::string BaseCluster::inqMeanCOPE() const
{ std::ostringstream s; s << m_impl->m_meanCOPE; return s.str(); }

struct Cluster::Implementation
{
  Implementation():
    m_maxZx(0), m_maxZy(0), m_maxZz(0),
    m_COGx(0.0), m_COGy(0.0), m_COGz(0.0),
    m_maxCOPEx(0), m_maxCOPEy(0), m_maxCOPEz(0)
  {}

  unsigned short m_maxZx, m_maxZy, m_maxZz;
  float m_COGx, m_COGy, m_COGz;
  unsigned short m_maxCOPEx, m_maxCOPEy, m_maxCOPEz;
};

Cluster::Cluster(const ColumnList& cl):
  BaseCluster(cl),
  m_impl(new Implementation)
{
  TRACKER("Cluster::Cluster(const ColumnList& cl)");

  m_headings.insert(std::make_pair("Z-MAX X (vox)",    ZMaxX));
  m_headings.insert(std::make_pair("Z-MAX Y (vox)",    ZMaxY));
  m_headings.insert(std::make_pair("Z-MAX Z (vox)",    ZMaxZ));
  m_headings.insert(std::make_pair("Z-COG X (vox)",    ZCOGX));
  m_headings.insert(std::make_pair("Z-COG Y (vox)",    ZCOGY));
  m_headings.insert(std::make_pair("Z-COG Z (vox)",    ZCOGZ));
  m_headings.insert(std::make_pair("COPE-MAX X (vox)", COPEMaxX));
  m_headings.insert(std::make_pair("COPE-MAX Y (vox)", COPEMaxY));
  m_headings.insert(std::make_pair("COPE-MAX Z (vox)", COPEMaxZ));
}

Cluster::Handle Cluster::create()
{
  STATIC_TRACKER("Cluster::create()");

  ColumnList cl;
  return Cluster::Handle(new Cluster(cl));
}

Cluster::Handle Cluster::create(const ColumnList& cl)
{
  STATIC_TRACKER("Cluster::create(const ColumnList& cl)");

  return Cluster::Handle(new Cluster(cl));
}

void Cluster::setCursorToMaxZ(Cursor::Handle& c) const
{
  return c->setCursor(m_impl->m_maxZx, m_impl->m_maxZy, m_impl->m_maxZz); 
}

void Cluster::setCursorToCOG(Cursor::Handle& c) const
{
  return c->setCursor(m_impl->m_COGx, m_impl->m_COGy, m_impl->m_COGz); 
}

void Cluster::setCursorToMaxCOPE(Cursor::Handle& c) const
{
  return c->setCursor(m_impl->m_maxCOPEx, m_impl->m_maxCOPEy, m_impl->m_maxCOPEz); 
}

bool Cluster::readColumn(std::istream& is, const std::string& col)
{
  bool flag(true);
  
  if(!is.eof())
    if(!BaseCluster::readColumn(is, col)) {

      HeadingMap::iterator it = m_headings.find(col);
      if(it != m_headings.end()) {
	switch(it->second) {
	case ZMaxX:    is >> m_impl->m_maxZx;    break;
	case ZMaxY:    is >> m_impl->m_maxZy;    break;
	case ZMaxZ:    is >> m_impl->m_maxZz;    break;
	case ZCOGX:    is >> m_impl->m_COGx;     break;
	case ZCOGY:    is >> m_impl->m_COGy;     break;
	case ZCOGZ:    is >> m_impl->m_COGz;     break;
	case COPEMaxX: is >> m_impl->m_maxCOPEx; break;
	case COPEMaxY: is >> m_impl->m_maxCOPEy; break;
	case COPEMaxZ: is >> m_impl->m_maxCOPEz; break;
	default: break;
	}
	flag = true;
      } else
	flag = false;
    }

  return flag;
}

void Cluster::scanFrom(std::istream& is)
{
  STATIC_TRACKER("Cluster::scanFrom(std::istream& is)");
  bool unknownFields(false);

  for(ColumnList::iterator it = m_columns.begin();
      it != m_columns.end(); ++it) {
    
    if(!readColumn(is, *it)) {
      std::string dummy;
      is >> dummy;
      unknownFields = true;
    } else
      initialised(true);
  }

  if(unknownFields)
    throw std::ios::failure("File contained unknown fields! Try re-running cluster.");

//   if(is.fail())
//     throw std::ios::failure("Cluster unexpected input!");
}

void Cluster::outputTo(std::ostream& os) const
{
  os << BaseCluster::m_impl->m_index << " ";
  os << BaseCluster::m_impl->m_voxels << " ";
  os << BaseCluster::m_impl->m_p << " ";
  os << BaseCluster::m_impl->m_minusLog10P << " ";
  os << BaseCluster::m_impl->m_maxZ << " ";
  os << m_impl->m_maxZx << " ";
  os << m_impl->m_maxZy << " ";
  os << m_impl->m_maxZz << " ";
  os << m_impl->m_COGx << " ";
  os << m_impl->m_COGy << " ";
  os << m_impl->m_COGz << " ";
  os << BaseCluster::m_impl->m_maxCOPE << " ";
  os << m_impl->m_maxCOPEx << " ";
  os << m_impl->m_maxCOPEy << " ";
  os << m_impl->m_maxCOPEz << " ";
  os << BaseCluster::m_impl->m_meanCOPE;
}

std::string Cluster::inqMaxZx() const
{ std::ostringstream s; s << m_impl->m_maxZx; return s.str(); }
std::string Cluster::inqMaxZy() const
{ std::ostringstream s; s << m_impl->m_maxZy; return s.str(); }
std::string Cluster::inqMaxZz() const
{ std::ostringstream s; s << m_impl->m_maxZz; return s.str(); }

std::string Cluster::inqMaxCOGx() const
{ std::ostringstream s; s << m_impl->m_COGx; return s.str(); }
std::string Cluster::inqMaxCOGy() const
{ std::ostringstream s; s << m_impl->m_COGy; return s.str(); }
std::string Cluster::inqMaxCOGz() const
{ std::ostringstream s; s << m_impl->m_COGz; return s.str(); }

std::string Cluster::inqMaxCOPEx() const
{ std::ostringstream s; s << m_impl->m_maxCOPEx; return s.str(); }
std::string Cluster::inqMaxCOPEy() const
{ std::ostringstream s; s << m_impl->m_maxCOPEy; return s.str(); }
std::string Cluster::inqMaxCOPEz() const
{ std::ostringstream s; s << m_impl->m_maxCOPEz; return s.str(); }

struct TalairachCluster::Implementation
{
  Implementation():
    m_maxZx(0.0), m_maxZy(0.0), m_maxZz(0.0),
    m_COGx(0.0), m_COGy(0.0), m_COGz(0.0),
    m_maxCOPEx(0.0), m_maxCOPEy(0.0), m_maxCOPEz(0.0)
  {}

  float m_maxZx, m_maxZy, m_maxZz;
  float m_COGx, m_COGy, m_COGz;
  float m_maxCOPEx, m_maxCOPEy, m_maxCOPEz;
};

TalairachCluster::TalairachCluster(const ColumnList& cl):
  BaseCluster(cl),
  m_impl(new Implementation)
{
  m_headings.insert(std::make_pair("Z-MAX X (mm)", ZMaxX));
  m_headings.insert(std::make_pair("Z-MAX Y (mm)", ZMaxY));
  m_headings.insert(std::make_pair("Z-MAX Z (mm)", ZMaxZ));
  m_headings.insert(std::make_pair("Z-COG X (mm)", ZCOGX));
  m_headings.insert(std::make_pair("Z-COG Y (mm)", ZCOGY));
  m_headings.insert(std::make_pair("Z-COG Z (mm)", ZCOGZ));
  m_headings.insert(std::make_pair("COPE-MAX X (mm)", COPEMaxX));
  m_headings.insert(std::make_pair("COPE-MAX Y (mm)", COPEMaxY));
  m_headings.insert(std::make_pair("COPE-MAX Z (mm)", COPEMaxZ));
}

TalairachCluster::Handle TalairachCluster::create()
{
  ColumnList cl;
  return TalairachCluster::Handle(new TalairachCluster(cl));
}

TalairachCluster::Handle TalairachCluster::create(const ColumnList& cl)
{
  return TalairachCluster::Handle(new TalairachCluster(cl));
}

void TalairachCluster::setCursorToMaxZ(ImageInfo::Handle& i, 
				       Cursor::Handle& c) const
{
  short x, y, z;
  i->mmToVoxCoord(m_impl->m_maxZx, m_impl->m_maxZy, m_impl->m_maxZz, 
		  x, y, z);
  return c->setCursor(x, y, z);
}

bool TalairachCluster::readColumn(std::istream& is, const std::string& col)
{
  bool flag(true);

  if(!BaseCluster::readColumn(is, col)) {

    HeadingMap::iterator it = m_headings.find(col);
    if(it != m_headings.end()) {
      switch(it->second) {
      case ZMaxX:is >> m_impl->m_maxZx; break;
      case ZMaxY:is >> m_impl->m_maxZy; break;
      case ZMaxZ:is >> m_impl->m_maxZz; break;
      case ZCOGX:is >> m_impl->m_COGx; break;
      case ZCOGY:is >> m_impl->m_COGy; break;
      case ZCOGZ:is >> m_impl->m_COGz; break;
      case COPEMaxX:is >> m_impl->m_maxCOPEx; break;
      case COPEMaxY:is >> m_impl->m_maxCOPEy; break;
      case COPEMaxZ:is >> m_impl->m_maxCOPEz; break;
      default: break;
      }
      flag = true;
    } else
      flag = false;
  }
  return flag;
}

void TalairachCluster::scanFrom(std::istream& is)
{
  STATIC_TRACKER("TalairachCluster::scanFrom(std::istream& is)");
  bool unknownFields(false);

  for(ColumnList::iterator it = m_columns.begin();
      it != m_columns.end(); ++it) {
    
    if(!readColumn(is, *it)) {
      std::string dummy;
      is >> dummy;
      unknownFields = true;
    } else
      initialised(true);
  }

  if(unknownFields)
    throw std::ios::failure("File contained unknown fields! Try re-running cluster.");

//   if(is.fail())
//     throw std::ios::failure("TalairachCluster unexpected input!");
}

void TalairachCluster::outputTo(std::ostream& os) const
{
  os << BaseCluster::m_impl->m_index << " ";
  os << BaseCluster::m_impl->m_voxels << " ";
  os << BaseCluster::m_impl->m_p << " ";
  os << BaseCluster::m_impl->m_minusLog10P << " ";
  os << BaseCluster::m_impl->m_maxZ << " ";
  os << m_impl->m_maxZx << " ";
  os << m_impl->m_maxZy << " ";
  os << m_impl->m_maxZz << " ";
  os << m_impl->m_COGx << " ";
  os << m_impl->m_COGy << " ";
  os << m_impl->m_COGz << " ";
  os << BaseCluster::m_impl->m_maxCOPE << " ";
  os << m_impl->m_maxCOPEx << " ";
  os << m_impl->m_maxCOPEy << " ";
  os << m_impl->m_maxCOPEz << " ";
  os << BaseCluster::m_impl->m_meanCOPE;
}

std::istream& operator>>(std::istream& is, BaseCluster::Handle& c)
{
  STATIC_TRACKER("std::istream& operator>>(std::istream& is, BaseCluster::Handle& c)");
  c->scanFrom(is);
  return is;
}

std::ostream& operator<<(std::ostream& os, const BaseCluster::Handle& c)
{
  c->outputTo(os);
  return os;
}

std::string TalairachCluster::inqMaxZx() const
{ std::ostringstream s; s << m_impl->m_maxZx; return s.str(); }
std::string TalairachCluster::inqMaxZy() const
{ std::ostringstream s; s << m_impl->m_maxZy; return s.str(); }
std::string TalairachCluster::inqMaxZz() const
{ std::ostringstream s; s << m_impl->m_maxZz; return s.str(); }

std::string TalairachCluster::inqMaxCOGx() const
{ std::ostringstream s; s << m_impl->m_COGx; return s.str(); }
std::string TalairachCluster::inqMaxCOGy() const
{ std::ostringstream s; s << m_impl->m_COGy; return s.str(); }
std::string TalairachCluster::inqMaxCOGz() const
{ std::ostringstream s; s << m_impl->m_COGz; return s.str(); }

std::string TalairachCluster::inqMaxCOPEx() const
{ std::ostringstream s; s << m_impl->m_maxCOPEx; return s.str(); }
std::string TalairachCluster::inqMaxCOPEy() const
{ std::ostringstream s; s << m_impl->m_maxCOPEy; return s.str(); }
std::string TalairachCluster::inqMaxCOPEz() const
{ std::ostringstream s; s << m_impl->m_maxCOPEz; return s.str(); }
