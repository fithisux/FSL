/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(SHAPE_H)
#define SHAPE_H


#include <boost/shared_ptr.hpp>
#include <list>
#include <vector>

#include <stack>
#include "storage/volume.h"

class QPainter;
class PixelRead;


struct Location
{
  int column;
  int row;
};


class Voxel
{

public:   
  typedef boost::shared_ptr< Voxel > Handle;  

  static Handle create(int, int, int, float);   
  
  virtual ~Voxel(){};
  int inqX() const {return m_x;}
  int inqY() const {return m_y;}
  int inqZ() const {return m_z;}
  float inqVal() const {return m_val;}
  void setValue(float v) { m_val = v; }
  void setDrawn(bool state) {m_drawn = state;}
  bool inqDrawn() const {return m_drawn;}

private: 
  Voxel(int, int, int, float);
  int m_x;
  int m_y;
  int m_z;
  float m_val;
  bool m_drawn;
};

class Shape
{
public:

  typedef boost::shared_ptr< Shape > Handle;
  static Handle create(QPainter * p, Volume::Handle, int orient, int slice);
  virtual  ~Shape();
  void draw();   
  void commit();
  void list();
  bool empty(){return m_commitVoxels.empty();}
  int  size() {return m_commitVoxels.size();}
  Shape::Handle getBuffer();  
  Shape::Handle getFloodBuffer();
  virtual void addVertex(int, int, int, float);
  virtual void floodFill(int, int, float newVal);

private:
  Shape() {}
  Shape(QPainter* p, Volume::Handle vol, int orient, int slice);  
  virtual void addSurroundingVoxels(Voxel::Handle &,int size, float val);
  void push_check(Voxel::Handle & pix, int size);
  virtual float readVoxel(int, int, int);
  virtual float readPixel(int, int);
  virtual void writeVoxel(int, int, int, float);
  virtual void writePixel(int, int, float);
  virtual void pushFloodUndoVoxel(int, int, int, float oldVal);
  virtual void pushFloodUndoPixel(int, int, float oldVal);
  virtual bool inRange(int, int);
  Voxel::Handle pixelToVoxel(int, int, float);

  void drawVertex(const Voxel::Handle&);
  void commitVertex(const Voxel::Handle&);

  std::vector<Voxel::Handle> m_voxels;
  std::vector<Voxel::Handle> m_commitVoxels;
  std::vector<Voxel::Handle> m_floodUndoVoxels;
  QPainter* m_paint;
  Volume::Handle m_volume;
  int m_orient;
  int m_slice;
  int m_counter;
  std::stack<Location> m_seedStack;
};

class Sphere: public Shape
{
private:
  virtual void addSurroundingVoxels(Voxel::Handle &,int size, float val);
};

class Cube: public Shape
{
private:
  virtual void addSurroundingVoxels(Voxel::Handle &,int size, float val);
};
#endif 
