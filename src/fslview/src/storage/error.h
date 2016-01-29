/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

/****************************************************************************
** $Id: error.h,v 1.2 2003/07/10 14:37:46 jim Exp $
**
** Copyright (C) 2002 University of Oxford.  All rights reserved.
**
*****************************************************************************/

#if !defined (ERROR_H)
#define ERROR_H

#include <string>

class FileError
{
public:
  FileError(const std::string& filename, const std::string& message);
  ~FileError();
  const std::string& inqFileName(){return m_filename;}
  const std::string& inqMessage() {return m_message;}
private:
  std::string m_filename;
  std::string m_message;
};

#endif
