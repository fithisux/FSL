/*  Sparse_Matrix.h

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

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

#if !defined(Sparse_Matrix_h)
#define Sparse_Matrix_h

#include <iostream>
#include "newmat.h"
#include <math.h>
#include <map>
#include <vector>
#include "newmatio.h"

using namespace NEWMAT;
using namespace std;

namespace MISCMATHS {  
  
  class SparseMatrix
    {
    public:

      typedef map<int,double> Row;
 
      SparseMatrix() : nrows(0), ncols(0) {}

      SparseMatrix(int pnrows, int pncols);

      SparseMatrix(const SparseMatrix& psm) 
	{
	  operator=(psm);
	}

      const SparseMatrix& operator=(const SparseMatrix& psm) 
	{
	  nrows = psm.nrows;
	  ncols = psm.ncols;
	  data = psm.data;

	  return *this;
	}

      SparseMatrix(const Matrix& pmatin) 
	{
	  operator=(pmatin);
	}

      const SparseMatrix& operator=(const Matrix& pmatin);

      //      void ReSize(int pnrows, int pncols)
      void ReSize(int pnrows, int pncols);

      void clear()
	{
	  ReSize(0,0);
	}
      
      void transpose(SparseMatrix& ret);
      
      ReturnMatrix RowAsColumn(int r) const;      

      int maxnonzerosinrow() const;

      void permute(const ColumnVector& p, SparseMatrix& pA);

      const double operator()(int x, int y) const
	{
	  double ret = 0.0;
	  map<int,double>::const_iterator it=data[x-1].find(y-1);
	  if(it != data[x-1].end())
	    ret = (*it).second;

	  return ret;
	}
      
      void set(int x, int y, double val) 
	{
	  data[x-1][y-1] = val;
	}

      void update(int x, int y, double val) 
	{
	  data[x-1][y-1] = val;
	}

      void insert(int x, int y, double val) 
	{
	  data[x-1].insert(Row::value_type(y-1,val));
	}

      void addto(int x, int y, double val) 
	{
	  if(val!=0)
	    data[x-1][y-1] += val;
	}

      void multiplyby(int x, int y, double val) 
	{
	  if((*this)(x,y)!=0)
	    data[x-1][y-1] *= val;
	}
            
      float trace() const;      

      Row& row(int r) { return data[r-1]; }

      const Row& row(int r) const { return data[r-1]; }

      ReturnMatrix AsMatrix() const;

      int Nrows() const { return nrows; }
      int Ncols() const { return ncols; }

      void multiplyby(double S);

      void vertconcatbelowme(const SparseMatrix& B); // me -> [me; B]
      void vertconcataboveme(const SparseMatrix& A); // me -> [A; me]
      void horconcat2myright(const SparseMatrix& B); // me -> [me B]
      void horconcat2myleft(const SparseMatrix& A);  // me -> [A me]

    private:
      
      int nrows;
      int ncols;

      vector<map<int,double> > data;

    };   

  void multiply(const SparseMatrix& lm, const SparseMatrix& rm, SparseMatrix& ret);
  void multiply(const DiagonalMatrix& lm, const SparseMatrix& rm, SparseMatrix& ret);

  void multiply(const SparseMatrix& lm, const ColumnVector& rm, ColumnVector& ret);

  void multiply(const SparseMatrix& lm, const SparseMatrix::Row& rm, ColumnVector& ret);

  void add(const SparseMatrix& lm, const SparseMatrix& rm, SparseMatrix& ret);

  void colvectosparserow(const ColumnVector& col, SparseMatrix::Row& row);

  void vertconcat(const SparseMatrix& A, const SparseMatrix& B, SparseMatrix& ret);

  void horconcat(const SparseMatrix& A, const SparseMatrix& B, SparseMatrix& ret);
}

#endif
