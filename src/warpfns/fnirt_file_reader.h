// Declarations of class used to decode and
// read files written by fnirt, and potentially
// by other pieces of software as long as they
// are valid displacement-field files.
//
// fnirt_file_reader.h
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2007-2012 University of Oxford 
/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk
    
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
    
    LICENCE
    
    FMRIB Software Library, Release 5.0 (c) 2012, The University of
    Oxford (the "Software")
    
    The Software remains the property of the University of Oxford ("the
    University").
    
    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.
    
    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.
    
    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.
    
    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. If you are
    interested in using the Software commercially, please contact Isis
    Innovation Limited ("Isis"), the technology transfer company of the
    University, to negotiate a licence. Contact details are:
    innovation@isis.ox.ac.uk quoting reference DE/9564. */

#ifndef fnirt_file_reader_h
#define fnirt_file_reader_h

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
#include "newmat.h"
#include "basisfield/basisfield.h"
#include "basisfield/splinefield.h"
#include "basisfield/dctfield.h"

namespace NEWIMAGE {

class FnirtFileReaderException: public std::exception
{
private:
  std::string m_msg;
public:
  FnirtFileReaderException(const std::string& msg) throw(): m_msg(msg) {}

  virtual const char * what() const throw() {
    return string("FnirtFileReader:: msg=" + m_msg).c_str();
  }

  ~FnirtFileReaderException() throw() {}
};

//
// Non-standard nifti intent-codes used for FSL-derived
// displacement/coefficient files.
//
#define FSL_FNIRT_DISPLACEMENT_FIELD       2006
#define FSL_CUBIC_SPLINE_COEFFICIENTS      2007
#define FSL_DCT_COEFFICIENTS               2008
#define FSL_QUADRATIC_SPLINE_COEFFICIENTS  2009

enum DisplacementFileType {FnirtSplineDispType, FnirtDCTDispType, FnirtFieldDispType, UnknownDispType, InvalidDispType};
enum AbsOrRelWarps {RelativeWarps, AbsoluteWarps, UnknownWarps};

// These are global routines that I have chosen to leave outside FnirtFileReader to
// allow them to be called independently. They probably don't really belong here, but
// I can't really come up with a better place for them.

void deffield2jacobian(const BASISFIELD::basisfield&   dx,
                       const BASISFIELD::basisfield&   dy,
                       const BASISFIELD::basisfield&   dz,
                       NEWIMAGE::volume<float>&        jac);

void deffield2jacobian(const BASISFIELD::basisfield&   dx,
                       const BASISFIELD::basisfield&   dy,
                       const BASISFIELD::basisfield&   dz,
                       const NEWMAT::Matrix&           aff,
                       volume<float>&                  jac);

void add_affine_part(const NEWMAT::Matrix&     aff,
                     unsigned int              indx,
                     NEWIMAGE::volume<float>&  warps);

void remove_affine_part(const NEWMAT::Matrix&     aff,
                        unsigned int              indx,
                        NEWIMAGE::volume<float>&  warps);

NEWMAT::Matrix estimate_affine_part(NEWIMAGE::volume4D<float>&  warps,
                                    unsigned int                every=1);


///////////////////////////////////////////////////////////////////////////////////////////////
//
// The FnirtFileReader is a class that reads, checks and decodes a displacement-field file.
// This can be a file created by fnirt, but it could also have been created for example
// from a fieldmap. For files created by other applications than fnirt it only reads 
// "proper" displacement fields. For fnirt-files it also reads coefficient-files.
//
// The ultimate purpose of FnirtFileReader and FnirtFileWriter is to provide an 
// interface that makes it easy and transparent not only to read and write displacement
// files, buat also to convert between the different representations. For example we
// want it to be easy to take a file with a direct representation of the field 
// (no basis-set) and write that out as a spline-coefficient file. Therefore we
// have internal space for all representations. When a file is read in it is stored
// in its "native" representation. When/if asked for another representation that
// representation is calculated and returned. 
//
///////////////////////////////////////////////////////////////////////////////////////////////

class FnirtFileReader
{
public:
  // Default constructor
  FnirtFileReader() : _type(InvalidDispType), _aor(UnknownWarps), _aff(IdentityMatrix(4)), _coef_rep(3) {}
  // Copy constructor
  FnirtFileReader(const FnirtFileReader& src);
  // Construct from file
  FnirtFileReader(const std::string& fname, AbsOrRelWarps wt=UnknownWarps, bool verbose=false)
  : _fname(fname)
  {
    common_read(fname,wt,verbose);
  }
  // Destructor. No clean up thanks to boost.
  ~FnirtFileReader() {}

  void Read(const std::string& fname, AbsOrRelWarps wt=UnknownWarps, bool verbose=false) {common_read(fname,wt,verbose);}
  bool IsValid() const {if (_type==InvalidDispType) return(false); else return(true);}
  std::string FileName() const {return(_fname);}
  DisplacementFileType Type() const {return(_type);}
  AbsOrRelWarps AbsOrRel() const {return(_aor);}
  std::vector<unsigned int> FieldSize() const;
  std::vector<double> VoxelSize() const;
  std::vector<unsigned int> KnotSpacing() const;
  unsigned int SplineOrder() const;
  std::vector<unsigned int> DCTOrder() const;
  NEWMAT::Matrix AffineMat() const {return(_aff);}
  NEWMAT::ReturnMatrix FieldAsNewmatMatrix(int indx=-1, bool inc_aff=false) const;
  NEWIMAGE::volume4D<float> FieldAsNewimageVolume4D(bool inc_aff=false) const;
  NEWIMAGE::volume<float> FieldAsNewimageVolume(unsigned int indx, bool inc_aff=false) const;
  NEWIMAGE::volume<float> Jacobian(bool inc_aff=false) const;
  BASISFIELD::splinefield FieldAsSplinefield(unsigned int indx, std::vector<unsigned int> ksp=std::vector<unsigned int>(0), unsigned int order=0) const;
  BASISFIELD::dctfield FieldAsDctfield(unsigned int indx, std::vector<unsigned int> order=std::vector<unsigned int>(0)) const;

protected:
  void common_read(const std::string& fname, AbsOrRelWarps wt, bool verbose);
  std::vector<boost::shared_ptr<BASISFIELD::basisfield> > read_coef_file(const NEWIMAGE::volume4D<float>&   vcoef,
                                                                         bool                               verbose) const;
  // void add_affine_part(NEWMAT::Matrix aff, unsigned int indx, NEWIMAGE::volume<float>& warps) const;

private:
  std::string                                                         _fname;
  DisplacementFileType                                                _type;
  AbsOrRelWarps                                                       _aor;
  NEWMAT::Matrix                                                      _aff;
  mutable std::vector<boost::shared_ptr<BASISFIELD::basisfield> >     _coef_rep;
  mutable boost::shared_ptr<NEWIMAGE::volume4D<float> >               _vol_rep;  
};

} // End namespace NEWIMAGE

#endif // end #ifndef fnirt_file_reader_h
