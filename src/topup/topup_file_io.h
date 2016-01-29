// Declarations of classes used to write and
// read files written by topup, and potentially
// by other pieces of software as long as they
// are valid magnetic-field files.
//
// topup_file_io.h
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2009 University of Oxford 
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

#ifndef topup_file_io_h
#define topup_file_io_h

#include <string>
#include <vector>
#include "newmat.h"
#include "newimage/newimage.h"
#include "basisfield/basisfield.h"
#include "basisfield/splinefield.h"
#include "basisfield/dctfield.h"

namespace TOPUP {

class TopupFileIOException: public std::exception
{
private:
  std::string m_msg;
public:
  TopupFileIOException(const std::string& msg) throw(): m_msg(msg) {}

  virtual const char * what() const throw() {
    return string("TopupFileIO:: msg=" + m_msg).c_str();
  }

  ~TopupFileIOException() throw() {}
};

//
// Non-standard nifti intent-codes used for FSL-derived
// displacement/coefficient files.
//
#define FSL_TOPUP_CUBIC_SPLINE_COEFFICIENTS      2016
#define FSL_TOPUP_QUADRATIC_SPLINE_COEFFICIENTS  2017
#define FSL_TOPUP_FIELD                          2018

enum FieldFileType {SplineField, FieldField, UnknownField, InvalidField};
enum TopupAbsOrRelWarps {RelativeWarps, AbsoluteWarps, UnknownWarps};

// Global function

NEWMAT::Matrix MovePar2Matrix(const NEWMAT::ColumnVector&     mp,
                              const NEWIMAGE::volume<float>&  vol);

NEWMAT::ColumnVector Matrix2MovePar(const NEWMAT::Matrix&           M,
				    const NEWIMAGE::volume<float>&  vol);

///////////////////////////////////////////////////////////////////////////////////////////////
//
// The TopupFileWrite is a class that writes off-resonance magnetic-field coefficient-files
// and movement parameter files produced by topup.
//
///////////////////////////////////////////////////////////////////////////////////////////////

class TopupFileWriter
{
public:
  // Constructor for coefficient file
  TopupFileWriter(const std::string&                 fname, 
                  const BASISFIELD::splinefield&     field);
  // Constructor for field file
  TopupFileWriter(const std::string&                   fname,
                  const NEWIMAGE::volume<float>&       ref,
                  const BASISFIELD::splinefield&       field);
  // Constructor for movement parameter file
  TopupFileWriter(const std::string&       fname,
                  const NEWMAT::Matrix&    mp);
  // Constructor for rigid body matrix files
  TopupFileWriter(const std::string&                    fname,
                  const std::vector<NEWMAT::Matrix>&    mp);
};                         

///////////////////////////////////////////////////////////////////////////////////////////////
//
// The TopupFileReader is a class that reads, checks and decodes an off-resonance 
// magnetic-field coefficient file.
// This can be a file created by topup, but it could also have been created for example
// from a fieldmap. For files created by other applications than topup it only reads 
// "proper" displacement fields. For topup-files it also reads coefficient-files.
//
// The ultimate purpose of TopupFileReader and TopupFileWriter is to provide an 
// interface that makes it easy and transparent to read and write off-resonance 
// field files. 
//
///////////////////////////////////////////////////////////////////////////////////////////////

class TopupFileReader
{
public:
  // Default constructor
  TopupFileReader() : _mp_valid(false) {}
  // Constructor for coeffcient-field and movement parameters
  TopupFileReader(const std::string& fname);
  // Destructor
  ~TopupFileReader() {}

  void ReadField(const std::string& fname) { common_read(fname); }
  void ReadMovements(const std::string& fname);
  NEWIMAGE::volume<float> FieldAsVolume() const { ensure_volume(); return(*_vol_rep); }
  boost::shared_ptr<NEWIMAGE::volume<float> > FieldAsVolumePtr() const { ensure_volume(); return(_vol_rep); }
  BASISFIELD::splinefield Field() const { ensure_field(); return(*_coef_rep); }
  boost::shared_ptr<BASISFIELD::splinefield> FieldPtr() const { ensure_field(); return(_coef_rep); }
  NEWMAT::Matrix MoveMatrix(unsigned int i) const { if (i>=_move.size()) TopupFileIOException("TopupFileReader::MoveMatrix: Index i out of bound"); return(_move[i-1]); }
  NEWMAT::ColumnVector MovePar(unsigned int i) const { if (i>=_mp.size()) TopupFileIOException("TopupFileReader::MovePar: Index i out of bound"); return(_mp[i-1]); }

private:
  bool                                                 _mp_valid;
  std::vector<NEWMAT::ColumnVector>                    _mp;
  std::vector<NEWMAT::Matrix>                          _move;
  mutable boost::shared_ptr<BASISFIELD::splinefield>   _coef_rep;
  mutable boost::shared_ptr<NEWIMAGE::volume<float> >  _vol_rep;

  void common_read(const std::string& fname);
  boost::shared_ptr<BASISFIELD::splinefield>  read_coef_file(const NEWIMAGE::volume<float>& vcoef);
  void read_movement(const std::string& fname, const NEWIMAGE::volume<float>& vol);
  void ensure_volume() const;
  void ensure_field() const;
  
};
 
///////////////////////////////////////////////////////////////////////////////////////////////
//
// The TopupDatafileReader is a class that reads text files containing phase-encode 
// vectors and readout times.
//
///////////////////////////////////////////////////////////////////////////////////////////////
 
class TopupDatafileReader
{
public:
  // Only one constructor
  TopupDatafileReader(const std::string& fname)
  {
    NEWMAT::Matrix tmp = read_ascii_matrix(fname);
    if (tmp.Ncols() != 4) {
      if (tmp.Nrows() == 4) tmp = tmp.t();
      else throw TopupFileIOException(string("TopupDatafileReader:: error reading file ")+fname);
    }
    _pv.ReSize(3,tmp.Nrows());
    _rt.ReSize(tmp.Nrows());
    for (int i=0; i<tmp.Nrows(); i++) {
      _pv.Column(i+1) = (tmp.SubMatrix(i+1,i+1,1,3)).t();
      if (fabs(1.0 - (_pv.Column(i+1)).NormFrobenius()) > 0.01) throw TopupFileIOException(string("TopupDatafileReader:: phase-encode vectors must be unity length")+fname);
      _rt(i+1) = tmp(i+1,4);
    }
  }

  // Destructor
  ~TopupDatafileReader() {}

  unsigned int          N() const { return(_pv.Ncols()); }
  NEWMAT::Matrix        PhaseEncodeVectors() const { return(_pv); }
  NEWMAT::ColumnVector  PhaseEncodeVector(unsigned int i) const { if (i<1 || i>N()) throw TopupFileIOException("TopupDatafileReader::PhaseEncodeVector: Invalid index"); else return(_pv.Column(i)); }
  NEWMAT::ColumnVector  ReadOutTimes() const { return(_rt); }
  double                ReadOutTime(unsigned int i) const { if (i<1 || i>N()) throw TopupFileIOException("TopupDatafileReader::ReadOutTime: Invalid index"); else return(static_cast<double>(_rt(i))); }
private:
  NEWMAT::Matrix        _pv;
  NEWMAT::ColumnVector  _rt;
};

} // End namespace TOPUP

#endif // end #ifndef topup_file_io_h
