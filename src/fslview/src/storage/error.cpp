/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "error.h"

FileError::FileError(const std::string& filename, const std::string& message):m_filename(filename),m_message(message)
{

}

FileError::~FileError()
{

}
