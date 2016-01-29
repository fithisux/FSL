/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(RECT_H)
#define RECT_H

#include <boost/shared_ptr.hpp>
#include <memory>
class Rect 
{
public:
  typedef boost::shared_ptr< Rect > Handle;
  static Handle createRect(int,int,int,int);

  Handle clone();
  void setRect(int,int,int,int);
  void setUnion(Rect::Handle);
  void setHeight(int);
  void setWidth(int);
  void translate(int,int);
  int top();
  int bottom();
  int left();
  int right();
  int height();
  int width();

  virtual ~Rect();

private:

  Rect(int,int,int,int); 
  struct Implementation;  
  const std::auto_ptr<Implementation> m_impl;
};















#endif
