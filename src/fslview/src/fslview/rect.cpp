/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "rect.h"
#include <assert.h>
struct Rect::Implementation
{
  Implementation(){};
  int m_top;
  int m_bottom;
  int m_left;
  int m_right;
};



Rect::Rect(int blX, int blY, int trX ,int trY):
  m_impl(new Implementation)
{ 
  setRect(blX,blY,trX,trY);
}

Rect::~Rect(){}

Rect::Handle Rect::createRect(int blX, int blY, int trX ,int trY)
{
  Handle dst(new Rect(blX,blY,trX,trY));
  return dst;
}

Rect::Handle Rect::clone()
{
  return Handle(new Rect(left(), bottom(), right(), top()));
}

void Rect::setRect(int botLeftX,int botLeftY,int topRightX,int topRightY)
{
  m_impl->m_bottom = std::min(topRightY,botLeftY);
  m_impl->m_top    = std::max(topRightY,botLeftY);
  m_impl->m_left   = std::min(topRightX,botLeftX);
  m_impl->m_right  = std::max(topRightX,botLeftX);  
}

int Rect::top(){return m_impl->m_top;}
int Rect::bottom(){return m_impl->m_bottom;}
int Rect::left(){return m_impl->m_left;}
int Rect::right(){return m_impl->m_right;}

int Rect::width()
{
  return m_impl->m_right - m_impl->m_left;
}

int Rect::height()
{
  return m_impl->m_top - m_impl->m_bottom;
}

void Rect::setWidth(int w)
{
  assert(w >= 0);
  float middle = (m_impl->m_left + m_impl->m_right)/2.0;
  m_impl->m_left  = (int)(middle - w/2.0);
  m_impl->m_right = (int)(middle + w/2.0);
}

void Rect::setHeight(int h)
{
  assert(h >=0);
  float middle = (m_impl->m_bottom + m_impl->m_top)/2.0;
  m_impl->m_bottom = (int)(middle - h/2.0);
  m_impl->m_top    = (int)(middle + h/2.0);
}

void Rect::translate(int x,int y)
{
  m_impl->m_left +=  x;
  m_impl->m_right += x;
  m_impl->m_bottom += y;
  m_impl->m_top += y;
}


void Rect::setUnion(Rect::Handle r)
{
  int left,bottom,right,top;

  int viewLeft   = r->left();
  int viewRight  = r->right();
  int viewTop    = r->top();
  int viewBottom = r->bottom();

  left   = std::max(viewLeft,m_impl->m_left);
  left   = std::min(left,viewRight);
  right  = std::max(viewLeft,m_impl->m_right);
  right  = std::min(right,viewRight);
  bottom = std::max(viewBottom,m_impl->m_bottom);
  bottom = std::min(bottom,viewTop);
  top    = std::max(viewBottom,m_impl->m_top);
  top    = std::min(top,viewTop);
  
  setRect(left, bottom, right, top); 
}
