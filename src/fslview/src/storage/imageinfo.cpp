/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "imageinfo.h"
#include "miscmaths/miscmaths.h"

#include <iostream>
#include <map>

using namespace MISCMATHS;

ImageInfo::Handle ImageInfo::clone()
{
  ImageInfo::Handle clone = ImageInfo::Handle(new ImageInfo(m_x,m_y,m_z,m_v,DT_SIGNED_SHORT,
							    m_xDim,m_yDim,m_zDim,
							    m_tr,m_bitsPerVoxel,m_imageName + "-mask",m_dtiImage));

  clone->setStdMat(inqStdMat());
  clone->m_sformcode = m_sformcode;
  clone->setRigidMat(inqRigidMat());
  clone->m_qformcode = m_qformcode;
  
//   if(m_sformcode != NIFTI_XFORM_UNKNOWN)
//     // use sform
//   else
//     // use qform

  clone->m_isStoredRadiological = m_isStoredRadiological;
  clone->m_intent = NIFTI_INTENT_LABEL;

  clone->setReadOnly(false);
  clone->setTarnished(false);

  return ImageInfo::Handle(clone);
}

//! @brief Set the default sform matrix
//!
//! @param xDim size of a voxel in X direction
//! @param yDim size of a voxel in Y direction
//! @param zDim size of a voxel in Z direction
//!
//! @returns The resultant 4x4 transformation matrix
//!
//! Sets up a default sform matrix with the scale factors necessary to achieve
//! voxel to mm coordinate conversion.
mat44 defaultStdMat(float xDim, float yDim, float zDim)
{
  mat44 stdmat;
  int i,j;
  for (i=0; i<4; i++) {
    for (j=0; j<4; j++) {
      stdmat.m[i][j]=0;
    }
  }
  stdmat.m[0][0]=fabs(xDim);
  stdmat.m[1][1]=fabs(yDim);
  stdmat.m[2][2]=fabs(zDim);
  stdmat.m[3][3]=1;
  return stdmat;
}

ImageInfo::ImageInfo(short x,short y,short z,short v,short dt,
                     float xDim, float yDim, float zDim,float tr,
                     size_t bitsPerVoxel,std::string name,
                     bool dtiImage):
  m_x(x),m_y(y),m_z(z),m_v(v),m_dt(dt),
  m_xDim(xDim),m_yDim(yDim),m_zDim(zDim),
  m_tr(tr),m_imageName(name),m_sformcode(NIFTI_XFORM_UNKNOWN),m_qformcode(NIFTI_XFORM_UNKNOWN),
  m_intent(NIFTI_INTENT_NONE),m_intentDescriptor("unset"),m_bitsPerVoxel(bitsPerVoxel),
  m_dtiImage(dtiImage), m_mainImage(false), m_isStoredRadiological(true),
  m_purpose(ImageIntent::Unknown)
{
  m_stdmat   = defaultStdMat(xDim,yDim,zDim);
  m_rigidmat = defaultStdMat(xDim,yDim,zDim);
}


ImageInfo::ImageInfo(FSLIO *avw,std::string filename):m_fileName(filename)
{
  float p1, p2, p3;
  size_t dim;
  m_bitsPerVoxel = FslGetDataType(avw, &m_dt);
  FslGetDim(avw, &m_x, &m_y, &m_z, &m_v);
  FslGetDimensionality(avw, &dim);
  if (dim == 3)
	  m_v = 1;
  FslGetVoxDim(avw, &m_xDim, &m_yDim, &m_zDim, &m_tr);
  FslGetAuxFile(avw,&m_auxFile[0]);
//   if(strlen(m_auxFile) != 0) {
//     FslGetCalMinMax(avw, &m_min, &m_max);
//   } else {
//     m_min = 0.0;
//     m_max = 0.0;
//   }
  FslGetCalMinMax(avw, &m_min, &m_max);

  m_sformcode = FslGetStdXform(avw,&m_stdmat);
  m_qformcode = FslGetRigidXform(avw,&m_rigidmat);
  if (m_sformcode==NIFTI_XFORM_UNKNOWN)  m_stdmat   = defaultStdMat(m_xDim,m_yDim,m_zDim);
  if (m_qformcode==NIFTI_XFORM_UNKNOWN)  m_rigidmat = defaultStdMat(m_xDim,m_yDim,m_zDim);

  FslGetIntent(avw,&m_intent,&p1,&p2,&p3);
  if (m_intent==NIFTI_INTENT_NONE) {
				// Try to guess intent from the file name
    string::size_type pos = m_fileName.find_last_of('/');
    string basename(m_fileName);
    if(pos != string::npos)
      basename = m_fileName.substr(pos);
    if (strstr(basename.c_str(),"stat")!=NULL) { m_intent=NIFTI_INTENT_ZSCORE; }
    if (strstr(basename.c_str(),"mask")!=NULL) { m_intent=NIFTI_INTENT_LABEL; }
    if (strstr(basename.c_str(),"tensor")!=NULL) 
      { m_intent=NIFTI_INTENT_GENMATRIX; m_intentDescriptor="DTI"; setDtiImage(true); }
  }
  m_purpose = ImageIntent::Unknown;

  m_isStoredRadiological = FslGetLeftRightOrder(avw) == FSL_RADIOLOGICAL;

  m_imageName = extractName(m_fileName);
  m_readOnly  = true;
  m_tarnished = false;
  m_dtiImage = false;
  m_mainImage = false;
}

ImageInfo::~ImageInfo()
{
}

void ImageInfo::setPurpose(ImageIntent::Code c)
{
  m_purpose = c;
  m_intent = NIFTI_INTENT_NONE;

  if(c == ImageIntent::Label)
    m_intent = NIFTI_INTENT_LABEL;
  if(c == ImageIntent::Statistic)
    m_intent = NIFTI_INTENT_ZSCORE;
  if(c == ImageIntent::Diffusion)
    m_intent=NIFTI_INTENT_GENMATRIX;
}

void ImageInfo::inqAxisOrientations(int& icode, int& jcode, int& kcode) const
{
  Matrix stand(Mat44ToNewmat(m_stdmat)), rigid(Mat44ToNewmat(m_rigidmat));

  get_axis_orientations(stand, m_sformcode, rigid, m_qformcode, icode, jcode, kcode);
}

std::string ImageInfo::inqDtAsString() const
{
  std::string result("unset");

  switch( inqDt() )
    {
    case 0:   result = std::string("None");          break;
    case 1:   result = std::string("Binary");        break;
    case 2:   result = std::string("Unsigned char"); break;
    case 4:   result = std::string("Signed short");  break;
    case 8:   result = std::string("Signed int");    break;
    case 16:  result = std::string("Float");         break;
    case 32:  result = std::string("Complex");       break;
    case 64:  result = std::string("Double");        break;
    case 128: result = std::string("RGB");           break;
    case 255: result = std::string("All");           break;
    case 256: result = std::string("Char");          break;
    case 512: result = std::string("Unsigned short");break;
    case 768: result = std::string("Unsigned int");  break;
    case 1024:result = std::string("Long long");     break;
    case 1280:result = std::string("Unsigned long"); break;
    case 1536:result = std::string("Long double");   break;
    
    default: result = std::string("Unknown DT: " + inqDt()); break;
    }

  return result;
}

/*
 * Commit the Header information and write it back to disc
 */
void ImageInfo::saveAvwHeader(FSLIO *avw)
{
  FslSetDataType(avw, m_dt);
  FslSetDim(avw, m_x, m_y, m_z, m_v);
  FslSetVoxDim(avw, m_xDim, m_yDim, m_zDim, m_tr);
  FslSetAuxFile(avw, &m_auxFile[0]);
  FslSetStdXform(avw,m_sformcode,inqStdMat());
  FslSetRigidXform(avw,m_qformcode,inqRigidMat());
  FslSetCalMinMax(avw, m_min, m_max);
  FslSetIntent(avw,m_intent,0.0,0.0,0.0);
  FslWriteHeader(avw);
}

void ImageInfo::setMin(float min) { m_min = min; }
void ImageInfo::setMax(float max) { m_max = max; }
float ImageInfo::inqMin() const { return m_min; }
float ImageInfo::inqMax() const { return m_max; }

std::string ImageInfo::inqLutName() const 
{ 
  return std::string(m_auxFile);
}

void ImageInfo::setLutName(std::string lutName)
{
  lutName.copy(m_auxFile,24);
}

void ImageInfo::setImageName(std::string imageName)
{
  m_imageName = imageName;
}

bool ImageInfo::isCompatible(ImageInfo::Handle other) const
{
  return ((other->inqX() == m_x) &&
          (other->inqY() == m_y) &&
          (other->inqZ() == m_z));
}

bool ImageInfo::isValidCoordinate(short x, short y, short z) const
{
  return(((x >= 0)&&(x < m_x)) &&
         ((y >= 0)&&(y < m_y)) &&
         ((z >= 0)&&(z < m_z)));
}

void ImageInfo::voxToMMCoord(short x, short y, short z, float& xmm, float& ymm, float& zmm) const
{
  if(m_sformcode != NIFTI_XFORM_UNKNOWN)
    FslGetMMCoord(inqStdMat(), x, y, z, &xmm, &ymm, &zmm);
  else
    FslGetMMCoord(inqRigidMat(), x, y, z, &xmm, &ymm, &zmm);
}

void ImageInfo::mmToVoxCoord(float xmm, float ymm, float zmm, short& x, short& y, short& z) const
{
  float xf, yf, zf;
  if(m_sformcode != NIFTI_XFORM_UNKNOWN)
    FslGetVoxCoord(inqStdMat(), xmm, ymm, zmm, &xf, &yf, &zf);
  else
    FslGetVoxCoord(inqRigidMat(), xmm, ymm, zmm, &xf, &yf, &zf);
  x = short(xf);
  y = short(yf);
  z = short(zf);
}

std::string ImageInfo::extractName(std::string filename)
{
  std::string result;
  std::string ext;

  result.erase();
  result = filename.substr(filename.rfind('/')+ 1,filename.length()- 1);

  return result;
} 

bool ImageInfo::isDtiCompatible() const
{
  bool retval=false;
  if ( (m_intent==NIFTI_INTENT_VECTOR) || (m_v > 2) )
    retval = true;
  if( m_purpose == ImageIntent::Diffusion ) retval = true;

  return retval;
}

bool ImageInfo::isInteger() const
{
  return (m_dt < 16);
}

bool  ImageInfo::inqNoDimensions() const
{
  return (m_xDim == 0.0 || m_yDim == 0.0 || m_zDim == 0.0);
}

bool ImageInfo::isStatImage() const
{
  bool retval=false;
  if ( (m_intent>=NIFTI_FIRST_STATCODE) && (m_intent<=NIFTI_LAST_STATCODE) )
    retval = true;
  if( m_purpose == ImageIntent::Statistic ) retval = true;

  return retval;
}

bool ImageInfo::isMaskImage() const
{
  bool retval=false;
  if (m_intent == NIFTI_INTENT_LABEL) retval = true;
  if (m_purpose == ImageIntent::Label) retval = true;
  return retval;
}

bool ImageInfo::hasValidXfms() const
{ 
  return ! ( (m_sformcode == NIFTI_XFORM_UNKNOWN) && 
	     (m_qformcode == NIFTI_XFORM_UNKNOWN) );
}

ImageCoordSystem::Code ImageInfo::inqCoordSystem() const
{
  return (m_sformcode != NIFTI_XFORM_UNKNOWN) ?
    ImageCoordSystem::Code(m_sformcode) :
    ImageCoordSystem::Code(m_qformcode);
}
