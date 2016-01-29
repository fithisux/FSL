
/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Rama Aravind Vorray
		James Saunders
		David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2002-2005 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(CLUSTERDATA_H)
#define CLUSTERDATA_H

#include <boost/shared_ptr.hpp>
#include <iostream>
#include <list>
#include <map>

#include "storage/imageinfo.h"
#include "cursor.h"

class BaseCluster
{
public:
  typedef boost::shared_ptr<BaseCluster> Handle;

  typedef enum {Index, Voxels, P, MinusLog10P, 
		ZMax,
		ZMaxX, ZMaxY, ZMaxZ, 
		ZCOGX, ZCOGY, ZCOGZ, 
		COPEMax,
		COPEMaxX, COPEMaxY, COPEMaxZ, 
		COPEMean} ColumnValue;

  typedef std::map<std::string, ColumnValue> HeadingMap;
  typedef std::list<std::string> ColumnList;
  
  virtual void scanFrom(std::istream&) = 0;
  virtual void outputTo(std::ostream&) const = 0;

  std::string inqIndex() const;
  std::string inqSize() const;
  std::string inqP() const;
  std::string inqMinusLog10P() const;

  std::string inqMaxZ() const;
  std::string inqMaxCOPE() const;
  std::string inqMeanCOPE() const;

  bool initialised() const;

  virtual std::string inqMaxZx() const = 0;
  virtual std::string inqMaxZy() const = 0;
  virtual std::string inqMaxZz() const = 0;

  virtual std::string inqMaxCOGx() const = 0;
  virtual std::string inqMaxCOGy() const = 0;
  virtual std::string inqMaxCOGz() const = 0;

  virtual std::string inqMaxCOPEx() const = 0;
  virtual std::string inqMaxCOPEy() const = 0;
  virtual std::string inqMaxCOPEz() const = 0;

  friend std::istream& operator>>(std::istream&, Handle&);
  friend std::ostream& operator<<(std::ostream&, const Handle&);

  virtual ~BaseCluster();
protected:
  BaseCluster(const ColumnList&);

  void initialised(bool);
  bool readColumn(std::istream& is, const std::string& col);

  struct Implementation;
  std::auto_ptr<Implementation> m_impl;

  HeadingMap m_headings;
  ColumnList m_columns;
};

class Cluster: public BaseCluster
{
public:
  typedef boost::shared_ptr<Cluster> Handle;

  static Handle create();
  static Handle create(const ColumnList& cl);

  virtual void scanFrom(std::istream&);
  virtual void outputTo(std::ostream&) const;

  virtual std::string inqMaxZx() const;
  virtual std::string inqMaxZy() const;
  virtual std::string inqMaxZz() const;

  virtual std::string inqMaxCOGx() const;
  virtual std::string inqMaxCOGy() const;
  virtual std::string inqMaxCOGz() const;
  
  virtual std::string inqMaxCOPEx() const;
  virtual std::string inqMaxCOPEy() const;
  virtual std::string inqMaxCOPEz() const;
  
  void setCursorToMaxZ(Cursor::Handle&) const;
  void setCursorToCOG(Cursor::Handle&) const;
  void setCursorToMaxCOPE(Cursor::Handle&) const;

private:
  Cluster(const ColumnList& cl);

  bool readColumn(std::istream& is, const std::string& col);

  struct Implementation;
  std::auto_ptr<Implementation> m_impl;
};

class TalairachCluster: public BaseCluster
{
public:
  typedef boost::shared_ptr<TalairachCluster> Handle;

  static Handle create();
  static Handle create(const ColumnList& cl);

  virtual void scanFrom(std::istream&);
  virtual void outputTo(std::ostream&) const;

  virtual std::string inqMaxZx() const;
  virtual std::string inqMaxZy() const;
  virtual std::string inqMaxZz() const;

  virtual std::string inqMaxCOGx() const;
  virtual std::string inqMaxCOGy() const;
  virtual std::string inqMaxCOGz() const;

  virtual std::string inqMaxCOPEx() const;
  virtual std::string inqMaxCOPEy() const;
  virtual std::string inqMaxCOPEz() const;

  void setCursorToMaxZ(ImageInfo::Handle&, Cursor::Handle&) const;
private:
  TalairachCluster(const ColumnList& cl);

  bool readColumn(std::istream& is, const std::string& col);

  struct Implementation;
  std::auto_ptr<Implementation> m_impl;
};

typedef std::list<BaseCluster::Handle> ClusterList;
typedef std::pair<ClusterList, ClusterList> ClusterListPair;
typedef std::pair<std::string, ClusterListPair> ClusterTable;

#endif
