/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(VERSION_H)
#define VERSION_H

//
// $Id: version.h,v 1.4 2003/12/01 15:49:27 jim Exp $
// $Log: version.h,v $
// Revision 1.4  2003/12/01 15:49:27  jim
// Merges from fslview-2_1_1 put in main trunk here
//
// Revision 1.3.8.1  2003/11/13 16:48:11  jim
// Tweeking changes
// VS: ----------------------------------------------------------------------
//
// Revision 1.3  2003/07/10 14:39:35  jim
// Added CopyRight Notices
//
// Revision 1.2  2002/12/15 18:31:04  flitney
// no message
//

extern const char *Version;
extern const char *Release;

class Version
{
public:   
  Version();
  ~Version();
};
#endif
