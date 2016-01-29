/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(IMAGEINFO_H)
#define IMAGEINFO_H

#include <boost/shared_ptr.hpp>
#include <string>

#include "fslio/fslio.h"

struct ImageIntent
{
  typedef enum { Unknown   = 0,
		 Label     = 1,
		 Statistic = 2,
		 Diffusion = 3 } Code;
};

struct ImageCoordSystem
{
  typedef enum { Unknown = 0,
		 ScannerAnatomical = 1,
		 AlignedAnatomical = 2,
		 Talairach = 3,
		 MNI_152 = 4 } Code;
};

class ImageInfo
{
  friend class testImageInfo;

public:
  typedef boost::shared_ptr< ImageInfo > Handle;

  ~ImageInfo();

  short inqX() const;
  short inqY() const;
  short inqZ() const;
  short inqNumVolumes() const;
  void setNumVolumes(short n);
  bool  inqNoDimensions() const;
  float inqXDim() const;
  float inqYDim() const;
  float inqZDim() const;

  // NEWIMAGE compatability interface
  inline int xsize() const { return inqX(); }
  inline int ysize() const { return inqY(); }
  inline int zsize() const { return inqZ(); }
  inline float xdim() const { return inqXDim(); }
  inline float ydim() const { return inqYDim(); }
  inline float zdim() const { return inqZDim(); }
 
  float inqTr() const;
  short inqDt() const;
  std::string inqDtAsString() const;
  size_t inqBytesPerVoxel() const;
  size_t inqBitsPerVoxel() const;
  std::string inqImageName() const;
  std::string inqFileName() const;
  std::string inqLutName() const;
  bool inqReadOnly() const;
  bool isStoredRadiological() const;
  bool hasValidXfms() const;
  bool isCompatible(ImageInfo::Handle) const;
  bool isValidCoordinate(short x, short y, short z) const;
  bool isMainImage() const;
  bool isInteger() const;
  bool isDtiCompatible() const;
  bool isDtiImage() const;
  bool isStatImage() const;
  bool isMaskImage() const;
  void setPurpose(ImageIntent::Code);
  ImageIntent::Code inqPurpose() const { return m_purpose; }
  ImageCoordSystem::Code inqCoordSystem() const;
  void setDtiImage(bool);
  void setImageName(std::string);
  void setReadOnly(bool);
  void setAsMainImage();
  void setLutName(std::string);
  void setTarnished(bool state);
  bool inqTarnished() const;     
  void inqAxisOrientations(int& icode, int& jcode, int& kcode) const;

  /** Image min instensity
   *  @return The stored minimum value or the minimum from the header if not yet set.
   */
  float inqMin() const;
  float inqMax() const;
  void setMin(float);
  void setMax(float);

  mat44 inqStdMat() const;
  void setStdMat(mat44 stdmat);

  mat44 inqRigidMat() const;
  void setRigidMat(mat44 rigidmat);

  void voxToMMCoord(short x, short y, short z, float& xmm, float& ymm, float& zmm) const;
  void mmToVoxCoord(float xmm, float ymm, float zmm, short& x, short& y, short& z) const;

  static ImageInfo::Handle create(short x,short y,short z,short v,short dt,
				  float xDim, float yDim, float zDim,float tr,
				  size_t bitsPerVoxel,std::string name, bool dtiImage) {
    return ImageInfo::Handle( new ImageInfo(x, y, z, v, dt,
					    xDim, yDim, zDim, tr,
					    bitsPerVoxel, name,dtiImage) );
  }
  static ImageInfo::Handle init(FSLIO *avw, std::string filename) {
    return ImageInfo::Handle( new ImageInfo(avw,filename) );
  }

  void saveAvwHeader(FSLIO *avw);

  Handle clone();

private:
  ImageInfo(FSLIO *avw,std::string filename);
  ImageInfo(short x,short y,short z,short v,short dt,
            float xDim, float yDim, float zDim,float tr,
            size_t bitsPerVoxel,std::string name, bool dtiImage);
  std::string extractName(std::string filename);


  bool m_minSet;
  bool m_maxSet;
  float m_min;
  float m_max;

  short m_x;
  short m_y;
  short m_z;
  short m_v;
  short m_dt;
  float m_xDim;
  float m_yDim;
  float m_zDim;
  float m_tr;
  char  m_auxFile[24];
  std::string m_fileName;
  std::string m_imageName;
 
  mat44 m_stdmat, m_rigidmat;
  short m_sformcode, m_qformcode;
  short m_intent;
  std::string m_intentDescriptor;

  size_t m_bitsPerVoxel;
  
  short m_type;  
  bool  m_readOnly;
  bool  m_tarnished;
  bool  m_dtiImage;
  bool  m_mainImage;
  bool  m_isStoredRadiological;
  ImageIntent::Code m_purpose;
};

inline void  ImageInfo::setTarnished(bool state) {m_tarnished = state;}
inline bool  ImageInfo::inqTarnished() const     {return m_tarnished;}

inline void  ImageInfo::setReadOnly(bool state) {m_readOnly = state;}
inline bool  ImageInfo::inqReadOnly() const {return m_readOnly;}

inline void  ImageInfo::setDtiImage(bool state) {m_dtiImage = state;}
inline bool  ImageInfo::isDtiImage() const {return m_dtiImage ; } //|| (m_purpose == ImageIntent::Diffusion) ;}

inline short ImageInfo::inqX() const { return m_x; }
inline short ImageInfo::inqY() const { return m_y; }
inline short ImageInfo::inqZ() const { return m_z; }
inline short ImageInfo::inqNumVolumes() const { return m_v; }
inline  void ImageInfo::setNumVolumes(short n) { m_v = n; m_tr = 0; m_dtiImage = false; }

inline float ImageInfo::inqXDim() const { return m_xDim; }
inline float ImageInfo::inqYDim() const { return m_yDim; }
inline float ImageInfo::inqZDim() const { return m_zDim; }
inline float ImageInfo::inqTr() const { return m_tr; }
inline short ImageInfo::inqDt() const { return m_dt; }

inline void ImageInfo::setAsMainImage(){m_mainImage = true;}
inline bool ImageInfo::isMainImage()const {return m_mainImage;}

inline std::string ImageInfo::inqImageName() const { return m_imageName; }
inline std::string ImageInfo::inqFileName() const { return m_fileName; }
inline size_t ImageInfo::inqBytesPerVoxel() const { return m_bitsPerVoxel / 8; }
inline size_t ImageInfo::inqBitsPerVoxel() const { return m_bitsPerVoxel; }

inline mat44 ImageInfo::inqStdMat() const { return m_stdmat; }
inline mat44 ImageInfo::inqRigidMat() const { return m_rigidmat; }
inline void ImageInfo::setStdMat(mat44 stdmat) { m_stdmat = stdmat; }
inline void ImageInfo::setRigidMat(mat44 rigidmat) { m_rigidmat = rigidmat; }

inline bool ImageInfo::isStoredRadiological() const { return m_isStoredRadiological; }

#endif
