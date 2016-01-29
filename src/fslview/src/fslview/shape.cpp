/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "shape.h"
#include "slicewidget.h"
#include <algorithm>
#include <math.h>
#include <qobject.h>
#include <qpainter.h>

//#define DEBUGGING
#include "tracker.h"

#include<boost/functional.hpp>

using namespace std;

Voxel::Handle Voxel::create(int x, int y, int z, float val)
{
  Handle dst(new Voxel(x, y, z, val));
  return dst;
}

Voxel::Voxel(int x, int y, int z, float val):
  m_x(x), m_y(y), m_z(z), m_val(val), m_drawn(false)
{
}

struct VoxelRead {
  VoxelRead(Volume::Handle vol, std::vector<Voxel::Handle>& v):
    m_volume(vol), m_voxels(v) {}

  void operator()(Voxel::Handle vo)
  {

    if( m_volume->inRange(vo->inqX(), vo->inqY(), vo->inqZ()) )
      {
	Voxel::Handle nv = 
	  Voxel::create( vo->inqX(), vo->inqY(), vo->inqZ(),
			 m_volume->value(vo->inqX(), vo->inqY(), vo->inqZ()) );

	m_voxels.push_back(nv);
      }
  }

  Volume::Handle m_volume;
  std::vector<Voxel::Handle>& m_voxels;
};


struct VoxelSearch 
{
  VoxelSearch(Voxel::Handle &v):
    m_found(false), m_v(v) {}

  void operator()(Voxel::Handle v) 
  {
    if ( (v->inqX() == m_v->inqX()) &&
	 (v->inqY() == m_v->inqY()) &&
	 (v->inqZ() == m_v->inqZ()) )
      {
        m_found = true;
      }
  }

  bool m_found;
  Voxel::Handle m_v;
};

Shape::Shape(QPainter* p, Volume::Handle vol, int orient, int slice):
             m_paint(p), m_volume(vol),
             m_orient(orient), m_slice(slice)
{
  m_voxels.clear();
  m_commitVoxels.clear();
  m_floodUndoVoxels.clear();
}

Shape::~Shape(){}

Shape::Handle Shape::create(QPainter* painter,
			    Volume::Handle vol, int orient, int slice)
{
  Handle dst(new Shape(painter, vol,orient,slice));
  return dst;
}

void Shape::drawVertex(const Voxel::Handle& vo) 
{ 
  m_paint->setPen(QColor(255,255,255) );

  int x(0), y(0);
  if(m_orient == SliceWidget::Axial)         { x = vo->inqX(); y = vo->inqY(); }
  else if(m_orient == SliceWidget::Sagittal) { x = vo->inqY(); y = vo->inqZ(); }
  else if(m_orient == SliceWidget::Coronal)  { x = vo->inqX(); y = vo->inqZ(); }
  
  m_paint->drawPoint(x, y);
  m_paint->drawRect(x, y, 2, 2);
}

void Shape::draw()
{
  TRACKER("Shape::draw()");

  std::copy(m_voxels.begin(), m_voxels.end(), back_inserter(m_commitVoxels));
  std::for_each(m_voxels.begin(), m_voxels.end(),
		boost::bind1st(boost::mem_fun(&Shape::drawVertex), this));

  m_voxels.clear();
}

void Shape::commitVertex(const Voxel::Handle& vo)
{
  m_volume->setValue(vo->inqX(),vo->inqY(),vo->inqZ(),vo->inqVal());
}

void Shape::commit()
{
  TRACKER("Shape::commit()");

  std::for_each(m_commitVoxels.begin(), m_commitVoxels.end(), 
		boost::bind1st(boost::mem_fun(&Shape::commitVertex), this));

  m_commitVoxels.clear();
}

void Shape::floodFill(int x, int y, float newVal)
{
  /*
    This flood algorithm has been taken from "CVu The Journal of the ACCU" 
    (August 2003 Volume 15 No 4) See www.accu.org for more details.
    The article was by James Holland.  The algorithm is very similar to the
    spans fill algorithm mentioned in "Foley and van Dam"(Computer Graphics
    Section 19.5.2 The Basic Filling Algorithms.
  */
  
  if(inRange(x,y))
    {
      float oldVal = readPixel(x,y);
      if(oldVal != newVal)
	{
	  Location seed_location = {x,y};
	  m_seedStack.push(seed_location);
    
	  while(!m_seedStack.empty())
	    {
	      Location location = m_seedStack.top();
	      Location locationOrig = m_seedStack.top();
	      //Push pixel so that fill can be undone
	      pushFloodUndoPixel(location.column, location.row, oldVal);
	      writePixel(location.column, location.row, newVal);
	      --location.column;
	      while(location.column >= 0 && 
		    readPixel(location.column,location.row) == oldVal)
		{            
		  pushFloodUndoPixel(location.column, location.row, oldVal);
		  writePixel(location.column, location.row, newVal);
		  --location.column;
		}
	      int extreme_left = location.column + 1;

	      location.column = m_seedStack.top().column + 1;
	      m_seedStack.pop();

	      while(inRange(location.column,location.row) &&
		    readPixel(location.column, location.row) == oldVal)
		{
		  pushFloodUndoPixel(location.column, location.row, oldVal);
		  writePixel(location.column, location.row, newVal);
		  ++location.column;
		}
	      int extreme_right = location.column - 1;
        
	      //Scan above the seed row

	      if(inRange(locationOrig.column, locationOrig.row + 1))
		{        
		  location.row = locationOrig.row + 1;

		  bool previous_pixel_is_border = true;
		  for (location.column = extreme_right;
		       location.column>= extreme_left;
		       --location.column)
		    {
		      if(previous_pixel_is_border &&
			 readPixel(location.column, location.row) == oldVal)
			{
			  m_seedStack.push(location);
			  previous_pixel_is_border = false;
			}
		      else
			if(readPixel(location.column, location.row) != oldVal)
			  {
			    previous_pixel_is_border = true;
			  }
		    }
		}
	      //Scan below the seed row
	      location.row = locationOrig.row - 1;
        
	      if(location.row >= 0)
		{
		  bool previous_pixel_is_border = true;
		  for (location.column = extreme_right;
		       location.column>= extreme_left;
		       --location.column)
		    {
		      if(previous_pixel_is_border &&
			 readPixel(location.column, location.row) == oldVal)
			{
			  m_seedStack.push(location);
			  previous_pixel_is_border = false;
			}
		      else
			if(readPixel(location.column, location.row) != oldVal)
			  {
			    previous_pixel_is_border = true;
			  }
		    }
		}
	    }
	}
    } 
}

Voxel::Handle Shape::pixelToVoxel(int x, int y, float val)
{
  TRACKER("Shape::pixelToVoxel(int x, int y, float val)");

  MESSAGE(QString("val = %1").arg(val));

  Voxel::Handle v;

  if(m_orient == SliceWidget::Axial)
    v = Voxel::create(x, y, m_slice, val);
  else if(m_orient == SliceWidget::Sagittal)
    v = Voxel::create(m_slice, x, y, val);
  else if(m_orient == SliceWidget::Coronal)
    v = Voxel::create(x, m_slice, y, val);

  return v;
}

void Shape::addVertex(int x, int y, int size, float val)
{  
  /*
    When a vertex is added a line of pixels is automatically drawn
    between the last vertex and the current vertex.  The basic concept
    of how this is acheived is roughly described in Foley and van Dam 
    section 3.2.1 "The Basic Incremental Algoritm"

    Different processes occur depending on wether the gradient of the line
    is less than or bigger than 1. This test is used to decide wether to
    increment x and calculate y or increment y and calculate x.
  */
  TRACKER("Shape::addVertex(int x, int y, int size, float val)");
  
  Voxel::Handle cur = pixelToVoxel(x, y, val);

  if(m_commitVoxels.empty()) {
    addSurroundingVoxels(cur, size, val);
  } else {
    Voxel::Handle prev = m_commitVoxels.back();
    
    int diffX(0), diffY(0), prevX(0), prevY(0);

    if(m_orient == SliceWidget::Axial) {
      diffX = cur->inqX() - prev->inqX(); prevX = prev->inqX();
      diffY = cur->inqY() - prev->inqY(); prevY = prev->inqY();
    } else if(m_orient == SliceWidget::Sagittal) {
      diffX = cur->inqY() - prev->inqY(); prevX = prev->inqY();
      diffY = cur->inqZ() - prev->inqZ(); prevY = prev->inqZ();
    } else if(m_orient == SliceWidget::Coronal) {
      diffX = cur->inqX() - prev->inqX(); prevX = prev->inqX();
      diffY = cur->inqZ() - prev->inqZ(); prevY = prev->inqZ();
    }

    if((diffX != 0)||(diffY != 0))
      {
        float grad,gradrecip;
        if(diffY == 0)
          {
            grad = 0;gradrecip = 999;
          }
        else if(diffX == 0)
          {
            grad = 999;gradrecip = 0;
          }
        else
          {
            grad = float(diffY)/float(diffX);
            gradrecip = 1.0/grad;
          }
        //increment x if grad is less than 1.0
        if(fabs(grad) <= 1.0)
          {
            float yValue;
            int xStep(0),yValueInt,yValueSign;
            while(abs(xStep)<= abs(diffX))
              {
                yValue     = grad * xStep + prevY;
                if(yValue < 0){yValueSign = -1;}else{yValueSign = 1;}
                yValueInt  = int(floor(fabs(yValue + 0.5))) * yValueSign;
                Voxel::Handle mid = pixelToVoxel(xStep + prevX,
						 yValueInt,
						 val);
		MESSAGE(QString("<=1 - val = %1").arg(mid->inqVal()));
                addSurroundingVoxels(mid,size,val);
                if(diffX < 0){--xStep;}
                else         {++xStep;}
              }

          }
        else
          {
            //increment y if grad is bigger than 1.0
            float xValue;
            int yStep(0),xValueInt,xValueSign;
            while(abs(yStep) <= abs(diffY))
              {                
                xValue     = gradrecip * yStep + prevX;  
                if(xValue < 0){xValueSign = -1;}else{xValueSign = 1;}
                xValueInt  = int(floor(fabs(xValue + 0.5))) * xValueSign;                        
                Voxel::Handle mid = pixelToVoxel(xValueInt,
						 yStep + prevY,
						 val); 
		MESSAGE(QString(">1 - val = %1").arg(mid->inqVal()));
                addSurroundingVoxels(mid,size,val);
                if(diffY < 0){--yStep;}
                else         {++yStep;}
              }

          }
      }
  }
}

Shape::Handle Shape::getBuffer()
{
  TRACKER("Shape::getBuffer()");

  Shape::Handle hnd; 
  hnd = Handle(new Shape(m_paint,m_volume,m_orient,m_slice));
  std::for_each(m_commitVoxels.begin(),
                m_commitVoxels.end(),
                VoxelRead(m_volume, hnd->m_commitVoxels));
  
  return hnd; 
}

Shape::Handle Shape::getFloodBuffer()
{
  Shape::Handle hnd; 
  hnd = Handle(new Shape(m_paint,m_volume,m_orient,m_slice));
  hnd->m_commitVoxels = m_floodUndoVoxels;

  return hnd; 
}

float Shape::readVoxel(int x, int y, int z)
{
  float value( m_volume->value(x,y,z) );

  return value;
}

float Shape::readPixel(int x, int y)
{
  float v(0.0);

  if(m_orient == SliceWidget::Axial)
    v = readVoxel(x, y, m_slice);
  else if(m_orient == SliceWidget::Sagittal)
    v = readVoxel(m_slice, x, y);
  else if(m_orient == SliceWidget::Coronal)
    v = readVoxel(x, m_slice, y);

  return v;
}

void Shape::writeVoxel(int x, int y, int z, float newVal)
{ 
  m_volume->setValue(x,y,z,newVal);
}

void Shape::writePixel(int x, int y, float v)
{
  if(m_orient == SliceWidget::Axial)
    writeVoxel(x, y, m_slice, v);
  else if(m_orient == SliceWidget::Sagittal)
    writeVoxel(m_slice, x, y, v);
  else if(m_orient == SliceWidget::Coronal)
    writeVoxel(x, m_slice, y, v);
}

void Shape::pushFloodUndoVoxel(int x, int y, int z, float oldVal)
{
  Voxel::Handle cur = Voxel::create(x,y,z,oldVal);
  m_floodUndoVoxels.push_back(cur);
}

void Shape::pushFloodUndoPixel(int x, int y, float v)
{
  if(m_orient == SliceWidget::Axial)
    pushFloodUndoVoxel(x, y, m_slice, v);
  else if(m_orient == SliceWidget::Sagittal)
    pushFloodUndoVoxel(m_slice, x, y, v);
  else if(m_orient == SliceWidget::Coronal)
    pushFloodUndoVoxel(x, m_slice, y, v);
}

bool Shape::inRange(int x, int y)
{
  bool  result(false);

  if(m_orient == SliceWidget::Axial)
    result = m_volume->inRange(x, y, m_slice);
  else if(m_orient == SliceWidget::Sagittal)
    result = m_volume->inRange(m_slice, x, y);
  else if(m_orient == SliceWidget::Coronal)
    result = m_volume->inRange(x, m_slice, y);

  return result;
}

void Shape::push_check(Voxel::Handle& vox, int size)
{
  std::vector<Voxel::Handle>::iterator start;
  int listsize = m_voxels.size();
  int checkAmount = size * size;
  if(listsize < checkAmount)
    start = m_voxels.begin();
  else
    start = m_voxels.end()-checkAmount;

  VoxelSearch search = std::for_each(start,
                                     m_voxels.end(),
                                     VoxelSearch(vox));
  if(!search.m_found) 
    m_voxels.push_back(vox);
}

void Shape::addSurroundingVoxels(Voxel::Handle & v, int size, float val)
{
  TRACKER("Shape::addSurroundingVoxels(Voxel::Handle & v, int size, float val)");

  MESSAGE(QString("val = %1").arg(val));

  int x(0), y(0);
  if(m_orient == SliceWidget::Axial)         { x = v->inqX(); y = v->inqY(); }
  else if(m_orient == SliceWidget::Sagittal) { x = v->inqY(); y = v->inqZ(); }
  else if(m_orient == SliceWidget::Coronal)  { x = v->inqX(); y = v->inqZ(); }

  for(int nx = x; nx < x + size;nx++)
    for(int ny = y; ny < y + size;ny++)
      {
	Voxel::Handle loc;

	if(m_orient == SliceWidget::Axial)
	  loc = Voxel::create(nx-size/2, ny-size/2, m_slice, val);
	else if(m_orient == SliceWidget::Sagittal)
	  loc = Voxel::create(m_slice, nx-size/2, ny-size/2, val);	
	else if(m_orient == SliceWidget::Coronal)
	  loc = Voxel::create(nx-size/2, m_slice, ny-size/2, val);	
	
	push_check(loc, size);
      }

  Voxel::Handle mid(Voxel::create(v->inqX(),v->inqY(),v->inqZ(),val));
  m_voxels.push_back(mid);
}

void Sphere::addSurroundingVoxels(Voxel::Handle & pix, int size, float val)
{
}

void Cube::addSurroundingVoxels(Voxel::Handle & pix, int size, float val)
{
}
