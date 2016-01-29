
/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "version.h"

const char *Version = "3";
const char *Release = "2.0";

Version::Version(){}
Version::~Version(){}

//! \mainpage FSLView
//! 
//! \section Installation
//!
//! \subsection Sources
//! The sources can be obtained from: www.fmrib.ox.ac.uk/fsldownloads

//
// Revision 1.13  2003/09/16 08:47:47  jim
// Updated version.cpp file for beta 2 release
//
// Revision 1.12  2003/09/03 14:40:43  jim
// Changed version.cpp and added version comments.
//
// Revision 1.11  2003/07/30 12:34:52  jim
// Image names can now be edited.
//
// Revision 1.10  2003/07/28 15:24:05  jim
// No changes except version.cpp
//
// Revision 1.9  2003/07/10 14:39:35  jim
// Added CopyRight Notices
//
// Revision 1.8  2003/07/07 13:22:50  jim
// *** empty log message ***
//
// Revision 1.7  2003/05/30 08:41:23  jim
// Checking in after merge with stable branch.
//
// Revision 1.6.2.2  2003/05/25 17:16:45  flitney
// Fixes problems with aux_file not working and OverlayList updates failing
// when locking/unlocking for mask editting.
//
// Revision 1.6.2.1  2003/05/23 17:49:25  flitney
// *** empty log message ***
//
//
// Revision 1.5  2003/04/03 12:38:34  jim
// Added graphics card diagnostics and disabling of lightbox view
//
// Revision 1.4  2003/03/13 14:17:01  jim
// last check in before build.
//

// 1)Command line feature added
// 2)Remove overlay now removes highlighted overlay
// 3)Graphs now have autoranging and respond to variations in fonts
// 4)Automatically assigns known LUTS
// 5)Can work with -ve scales of images

// New for version 1.0 release 4

// 1) Graphics driver diagnostics, "About Graphics" dialog
// 2) Disables lightbox view if not supported

// New for version 1.0 beta 8

// 1) Cursor key control of cursor on imagewidgets.
// 2) Page Up/Down keys control slice depth.
// 3) Warns user that unsaved data may be lost.
// 4) Pen , Erase and Undo but
//	3.  	file overrides above if LUT set in AuxFile field
//
//Improved font selection

//Version 2.3.5
//Fixes 
// Override fatal error if FSLOUTPUTTYPE not present.
// Poor design of OverlayWidget causing bad display under MacOSX Aqua.
// RT1327: masks in fslview.
// RT1304: fslview and dtiimages
// RT1128: Timeseries cross
//tons added.

// New for version 2.0 beta 1

// 1) Non-openGL implementation.
// 2) Lightbox view now fully functional.
// 3) Single widget can now display sag/cor/axi views.
// 4) Movie through time button and Slice Roll buttons introduced.
// 5) Ctrl key turns pen mode to erase mode.
// 6) Ctrl-Z causes undo.
// 7) mm to Voxel coordinate transform corrected.
// 8) Rounding error prominant in diffusion images fixed.

//Version 2.0 beta 2

// 1) Fixed a number of floating point issues on alpha.
// 2) Resolved display issues on Mac.

//Version 2.0 beta 3

// 1) Implemented DTI functionality: RGB, Lines and modulation.
// 2) Changed blendImage to copy with transparency modulation.
// 3) Got rid of references to mainImage. All images now in lists.
// 4) Made metaimage , imagedata and imagestore classes to 
//    cope with collections of data.
// 5) Moved all pixel manipulation code into imagebuffer.cpp
// 6) Second slider now appears on overlaywidget for DTI mod images.

//Version 2.1.0

//Latest bugfixes:
//Ticket #700 Floating point exception when loading images.
//Ticket #600 x mm box behaves incorrectly with certain images.
//Ticket #710 Crash with timeseries of new mask
//Ticket #654 Disable Mod combo when in Dti Lines mode
//Ticket #581 renaming images, cursor always moves to the end of line.

//Latest Featuers:
//Ticket #490 Masking Fill In Closed Loop shape

//Version 2.1.1

//Latest bugfixes:
//738 Filler tool can't be used as an eraser
//771 Pen size tool has been removed

//Comments:
//Worked to reintroduce pen size control. Improved mechanism so that
//only pixels that haven't been drawn before are drawn in mask mode.


//Version 2.2.0

//Features
//656 TimeSeries Demeaning added
//Removed timeseries browse button, feature now constantly enabled.
//406 Help Information
//653 Turn off modulation layer
//727 Timeseries graph title
//777 Disable timeseries if only one volume
//Redesigned timeseries widgets. Now single, grid and cube exists.
//New Histogram features introduced including print and zoom.
//New bricon controls introduced.

//Version 2.2.1

//Is able to read and write zip files.

//Version 2.2.2

//Is a version with .dsw file altered so that it will build on windows.  
//This version is the first one to use the licenced Qt3 for windows!!!

//Version 2.3.0

//409 Orthagonal View Sizes
//Hdr and img files no longer have to be both zipped or unzipped.
//Warnings given when zipped and unzipped versions of a file found in dir.
//Clicking on timeseries plot moves image to relevant volume.

//Version 2.3.1

//Floating point errors on alpha due to new zooming fixed.
//Zoomhistory code reworked.

//Version 2.3.3

//Niftiio support
//Talaraich coordinate fix
//Rotating colourmaps
//Multiple command line args

//Version 2.3.4

//New features 
//ImageLeftRightEncoding preserved on screen (radiological view)
// Better support for voxel-mm translations
// Improved LookUpTable selection on startup:
//	1.  	stats/mask images load RedYellow
//	2.  	stats add as next LUT in list
//	3.  	file overrides above if LUT set in AuxFile field
//
//Improved font selection

//Version 2.3.5
//Fixes 
// Override fatal error if FSLOUTPUTTYPE not present.
// Poor design of OverlayWidget causing bad display under MacOSX Aqua.
// RT1327: masks in fslview.
// RT1304: fslview and dtiimages
// RT1128: Timeseries cross

//Version 2.3.6
//Fixes
//  RT1367: DTI - RGB DISPLAY RANGE.
//  RT1407: Histogram binning.
//  RT1450: luts and aux_file etc.
//  RT1474: Shortcut options cause crash.
//  RT1504: Open/Add doesn't update working directory
//New Features
//  Option to create a 3D/4D overlay mask

//Version 3.0
//New Features
//  3D VTK-based viewer

//Version 3.0.1
//Fixes
//  Mistakes in 3D viewer - thresholding of stats images

//Version 3.0.3
//Fixes
//  3D viewer - more thresholding errors

//Version 3.0.4/5
//Feature change
//  Mesh view - facets coloured from Cell Data, not interpolated Point Data.

//Version 3.1.0
//Feature change
//  OrthoView now has Traditional, Row and Column views
//  Allows you to turn off negative lut
//  Volume selection now a per image display setting

//Version 3.1.6
//Feature change
//  Added DTI Lines(RGB) mode; lines colored to indictate orientation

//Version 3.1.7
//Fixes
//  VTK mode crashing on startup - back out some changes to vtkwidget.cpp from 3.1.6 

