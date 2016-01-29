// Utility for resampling of DTI data using fields/data
// estimated by topup.
//
// applytopup.h
//
// Jesper Andersson, FMRIB Image Analysis Group
//
//Copyright (C) 2012 University of Oxford  

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

#ifndef applytopup_h
#define applytopup_h

#include <cstring>
#include "newmat.h"
#include "miscmaths/miscmaths.h"
#include "newimage/newimageall.h"
#include "warpfns/fnirt_file_reader.h"
#include "topup_file_io.h"

namespace TOPUP {

class ApplyTopupException: public std::exception
{
private:
  std::string m_msg;
public:
  ApplyTopupException(const std::string& msg) throw(): m_msg(msg) { cout << what() << endl; }

  virtual const char * what() const throw() {
    return string("ApplyTopup:: msg=" + m_msg).c_str();
  }

  ~ApplyTopupException() throw() {}
};

// Helper-class for divding a set of scans into collections that can be used for least squares restoration

class TopupCollections
{
public:
  TopupCollections(const std::vector<unsigned int>& indx,
                   const TopupDatafileReader&       dfile);
  ~TopupCollections() {}
  unsigned int NCollections() { return(indxs.size()); }
  unsigned int NIndices(unsigned int c) { if (c>=NCollections()) throw ApplyTopupException("TopupCollections::NIndicies: Index out of range"); else return(indxs[c].size()); }
  unsigned int NScans(unsigned int c) { return(NIndices(c)); }
  unsigned int IndexAt(unsigned int c, unsigned int i);
  unsigned int ScanAt(unsigned int c, unsigned int i);
private:
  std::vector<std::vector<unsigned int> >  indxs;
  std::vector<std::vector<unsigned int> >  scans;
};

// Global functions used by applytopup

// Does the job
int applytopup();
// Helper functions
NEWIMAGE::volume4D<float> vb_resample_4D(const TopupDatafileReader&                       datafile, 
                                         const TopupFileReader&                           topupfile, 
                                         const std::vector<unsigned int>&                 inindices, 
                                         NEWIMAGE::volume<char>&                          mask,
                                         std::vector<NEWIMAGE::volume4D<float> >&         scans);
NEWIMAGE::volume4D<float> vb_resample_3D(const TopupDatafileReader&                       datafile, 
                                         const TopupFileReader&                           topupfile, 
                                         const std::vector<unsigned int>&                 inindices, 
                                         NEWIMAGE::volume<char>&                          mask,
                                         std::vector<NEWIMAGE::volume4D<float> >&         scans);
double tau_update_contrib(const std::vector<std::vector<NEWMAT::ColumnVector> >&        mu,
                          const std::vector<std::vector<MISCMATHS::SpMat<double> > >&   Lambda,
		          const std::vector<std::vector<MISCMATHS::SpMat<double> > >&   Lijk);
double tau_update(const std::vector<std::vector<std::vector<NEWMAT::ColumnVector> > >&  mu,
                  const std::vector<std::vector<MISCMATHS::SpMat<double> > >&           Lambda,
		  const std::vector<std::vector<MISCMATHS::SpMat<double> > >&           Lijk,
                  double                                                                tau_0);
double tau_update(const std::vector<std::vector<NEWMAT::ColumnVector> >&        mu,
                  const std::vector<std::vector<MISCMATHS::SpMat<double> > >&   Lambda,
		  const std::vector<std::vector<MISCMATHS::SpMat<double> > >&   Lijk,
                  double                                                        tau_0);
NEWMAT::ColumnVector mu_update(const std::vector<std::vector<NEWMAT::ColumnVector> >&        mu,
                               const MISCMATHS::SpMat<double>&                               K,
		               const NEWMAT::ColumnVector&                                   y,
                               const NEWMAT::CroutMatrix&                                    Lambda,
                               const std::vector<std::vector<MISCMATHS::SpMat<double> > >&   Lijk,
                               double                                                        phi,
                               double                                                        lambda,
                               int                                                           sl,
                               int                                                           row);
NEWIMAGE::volume4D<float> vb_resample_2D(const TopupDatafileReader&                       datafile, 
                                         const TopupFileReader&                           topupfile, 
                                         const std::vector<unsigned int>&                 inindices, 
                                         NEWIMAGE::volume<char>&                          mask,
                                         std::vector<NEWIMAGE::volume4D<float> >&         scans);
std::vector<std::vector<MISCMATHS::SpMat<double> > > GetLijk(int                        sz, 
                                                             const std::vector<double>  sf);
std::vector<MISCMATHS::SpMat<double> > GetLij(int                                       sz, 
                                   double                                               sf);
MISCMATHS::SpMat<double> GetLii(int                                                     sz);
NEWIMAGE::volume4D<float> jac_resample(const TopupDatafileReader&                       datafile, 
                                       const TopupFileReader&                           topupfile, 
                                       const std::vector<unsigned int>&                 inindices, 
                                       NEWIMAGE::volume<char>&                          mask,
                                       std::vector<NEWIMAGE::volume4D<float> >&         scans,
                                       NEWIMAGE::interpolation                          interp);
NEWIMAGE::volume4D<float> lsr_resample(const TopupDatafileReader&                       datafile, 
                                       const TopupFileReader&                           topupfile, 
                                       const std::vector<unsigned int>&                 inindices, 
                                       NEWIMAGE::volume<char>&                          mask,
                                       std::vector<NEWIMAGE::volume4D<float> >&         scans);
void resample_using_movement_parameters(const TopupDatafileReader&                      datafile,
                                        const TopupFileReader&                          topupfile,
                                        const std::vector<unsigned int>&                inindices,
                                        NEWIMAGE::volume<char>&                         mask,
                                        std::vector<NEWIMAGE::volume4D<float> >&        scans,
                                        NEWIMAGE::interpolation                         interp);
std::vector<unsigned int> parse_commaseparated_numbers(const std::string& list);
std::vector<std::string> parse_commaseparated_list(const std::string&  list);
bool row_col_is_alright(const NEWIMAGE::volume<char>&   mask,
                        unsigned int                    k,
                        unsigned int                    ij,
                        bool                            row);
NEWMAT::ReturnMatrix extract_row_col(const NEWIMAGE::volume<float>&  vol,
                                     unsigned int                    k,
                                     unsigned int                    ij,
                                     bool                            row);
void add_to_slice(NEWIMAGE::volume<float>&                  vol,
                  const std::vector<NEWMAT::ColumnVector>&  mu,
                  unsigned int                              sl,
                  bool                                      row);
void add_to_rows_cols(NEWIMAGE::volume4D<float>&  vol,
                      const NEWMAT::Matrix&       B,
                      unsigned int                k,
                      unsigned int                ij,
                      bool                        row);
void zero_out_rows_cols(NEWIMAGE::volume4D<float>&   vol,
                        const NEWMAT::Matrix&        map,
                        bool                         row);


}      // End namespace TOPUP

#endif // End ifndef applytopup_h
