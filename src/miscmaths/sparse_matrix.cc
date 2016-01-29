/*  sparse_matrix.h

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

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#define WANT_STREAM
#define WANT_MATH

#include "sparse_matrix.h"
#include "newmatio.h"
#include "newmat.h"
#include "miscmaths.h"
#include "utils/tracer_plus.h"

using namespace std;
using namespace Utilities;
using namespace NEWMAT;
using namespace MISCMATHS;

namespace MISCMATHS {

  SparseMatrix::SparseMatrix(int pnrows, int pncols) : 
	nrows(pnrows),
	ncols(pncols),
	data(pnrows) 
  {
  }

  void SparseMatrix::ReSize(int pnrows, int pncols)
  {
    nrows = pnrows;
    ncols = pncols;
    
    data.clear();
    data.resize(nrows);
  }

  const SparseMatrix& SparseMatrix::operator=(const Matrix& pmatin) 
  {
    data.clear();
    data.resize(pmatin.Nrows());
    nrows = pmatin.Nrows();
    ncols = pmatin.Ncols();

    for(int r=1; r <= pmatin.Nrows(); r++)
      {	
	for(int c=1; c <= pmatin.Ncols(); c++)
	  {
	    if(pmatin(r,c)!=0)
	      insert(r,c,pmatin(r,c));
	  }
      }
    return *this;
  }     
 
  void SparseMatrix::transpose(SparseMatrix& ret)
  {
    Tracer_Plus tr("SparseMatrix::transpose");

    ret.ReSize(ncols,nrows);
    
    for(int r=1; r <= nrows; r++)
      for(map<int,double>::const_iterator it=data[r-1].begin(); it!=data[r-1].end(); it++)
	ret.insert((*it).first+1, r, (*it).second);
    
  }
  
  int SparseMatrix::maxnonzerosinrow() const
  {
    int mx = 0;
    for(int r=1; r <= nrows; r++)
      {
	int si = data[r-1].size();
	if(si > mx)
	  mx = si;
      }
    return mx;
  }
  
  void SparseMatrix::permute(const ColumnVector& p, SparseMatrix& pA)
  {
    Tracer_Plus tr("SparseMatrix::permute");

    pA.ReSize(nrows,ncols);
    
    ColumnVector ip(p.Nrows());
    for(int r=1; r <= nrows; r++)
      ip(int(p(r))) = r;
    
    for(int r=1; r <= nrows; r++)
      for(map<int,double>::const_iterator it=data[r-1].begin(); it!=data[r-1].end(); it++)
	{
	  pA.insert(int(ip(r)), int(ip((*it).first+1)), (*it).second); 
	}
  }
  
  ReturnMatrix SparseMatrix::AsMatrix() const
  {
    Matrix ret(nrows,ncols);
    ret = 0;	  
    
    for(int r=1; r <= nrows; r++)
      for(map<int,double>::const_iterator it=data[r-1].begin(); it!=data[r-1].end(); it++)
	ret(r,(*it).first+1) = (*it).second;
    
    ret.Release();
    return ret;
  }

  float SparseMatrix::trace() const
  {   
    float tr = 0.0;
    for(int k = 1; k<=Nrows(); k++)
      {
	tr += (*this)(k,k);
      }

    return tr;
  }

 void SparseMatrix::vertconcatbelowme(const SparseMatrix& B)
 {
   Tracer_Plus tr("SparseMatrix::vertconcatbelowme");

   if (Ncols() != B.Ncols()) {throw Exception("Cols don't match in SparseMatrix::vertconcatbelowme");}

   data.resize(Nrows()+B.Nrows());
   for (int i=1; i<=B.Nrows(); i++) {
     this->row(Nrows()+i) = B.row(i);
   }

   nrows += B.Nrows();
 }

 void SparseMatrix::vertconcataboveme(const SparseMatrix& A)
 {
   Tracer_Plus tr("SparseMatrix::vertconcataboveme");

   if (Ncols() != A.Ncols()) {throw Exception("Cols don't match in SparseMatrix::vertconcataboveme");}

   data.resize(Nrows()+A.Nrows());
   for (int i=Nrows(); i>=1; i--) {
     this->row(i+A.Nrows()) = this->row(i);
   }
   for (int i=1; i<=A.Nrows(); i++) {
     this->row(i) = A.row(i);
   }

   nrows += A.Nrows();
 }
  
 void SparseMatrix::horconcat2myright(const SparseMatrix& B)
 {
   Tracer_Plus tr("SparseMatrix::horconcat2myright");

   if (Nrows() != B.Nrows()) {throw Exception("Rows don't match in SparseMatrix::vertconcat2myright");}

   for (int i=1; i<=Nrows(); i++) {
     const SparseMatrix::Row& tmpRow = B.row(i);
     for (SparseMatrix::Row::const_iterator it=tmpRow.begin(); it!=tmpRow.end(); it++) {
       this->insert(i,Ncols()+int(it->first)+1,double(it->second));
     }
   }
   ncols += B.Ncols();
 }

 void SparseMatrix::horconcat2myleft(const SparseMatrix& A)
 {
   Tracer_Plus tr("SparseMatrix::horconcat2myright");

   if (Nrows() != A.Nrows()) {throw Exception("Rows don't match in SparseMatrix::vertconcat2myleft");}

   for (int i=1; i<=Nrows(); i++) {
     SparseMatrix::Row oldRow = this->row(i);
     this->row(i) = SparseMatrix::Row();  // Empty row.
     const SparseMatrix::Row& tmpRow = A.row(i);
     for (SparseMatrix::Row::const_iterator it=tmpRow.begin(); it!=tmpRow.end(); it++) {
       this->insert(i,int(it->first)+1,double(it->second));
     }
     for (SparseMatrix::Row::const_iterator it=oldRow.begin(); it!=oldRow.end(); it++) {
       this->insert(i,A.Ncols()+int(it->first)+1,double(it->second));
     }
   }     
   ncols += A.Ncols();
 }

 ReturnMatrix SparseMatrix::RowAsColumn(int r) const
  {
    Tracer_Plus tr("SparseMatrix::RowAsColumn");

    ColumnVector ret;
    ret.ReSize(ncols);
    ret = 0;
    
    const SparseMatrix::Row& rowtmp = row(r);
    for(SparseMatrix::Row::const_iterator it=rowtmp.begin();it!=rowtmp.end();it++)
      {
	int c = (*it).first+1;	     	      
	double val = (*it).second;
	ret(c) = val;
      }
    
    ret.Release();
    return ret;
  }
  
  void colvectosparserow(const ColumnVector& col, SparseMatrix::Row& row)
  {
    Tracer_Plus tr("SparseMatrix::colvectosparserow");
    for(int j = 1; j<=col.Nrows(); j++)
      {
	if(std::abs(col(j))>1e-4)
	  row[j-1] = col(j);
      }
  }

  void SparseMatrix::multiplyby(double S)
  {
    Tracer_Plus tr("SparseMatrix::multiplyby");

    for(int j = 1; j<=Nrows(); j++)
      {
	SparseMatrix::Row& row = (*this).row(j);
	for(SparseMatrix::Row::iterator it=row.begin();it!=row.end();it++)
	  {
	    (*it).second *= S;
	  }
      }
  }
  
  void multiply(const SparseMatrix& lm, const SparseMatrix& rm, SparseMatrix& ret)
  {
    Tracer_Plus tr("SparseMatrix::multiply");

    int nrows = lm.Nrows();
    int ncols = rm.Ncols();

    if(lm.Ncols() != rm.Nrows()) throw Exception("Rows and cols don't match in SparseMatrix::multiply");

    ret.ReSize(nrows,ncols);

    for(int j = 1; j<=nrows; j++)
      {
	const SparseMatrix::Row& row = lm.row(j);	
	for(SparseMatrix::Row::const_iterator it=row.begin();it!=row.end();it++)
	  {
	    int c = (*it).first+1;
	    double val = (*it).second;
	    for(int k = 1; k<=ncols; k++)
	      {
		ret.addto(j,k,val*rm(c,k));
	      }
	  }
      }

  }
  
  void multiply(const SparseMatrix& lm, const ColumnVector& rm, ColumnVector& ret)
  {
    Tracer_Plus tr("SparseMatrix::multiply2");

    int nrows = lm.Nrows();   
    
    if(lm.Ncols() != rm.Nrows()) throw Exception("Rows and cols don't match in SparseMatrix::multiply");
    
    ret.ReSize(nrows);
    
    for(int j = 1; j<=nrows; j++)
      {
	float sum = 0.0;
	const SparseMatrix::Row& row = lm.row(j);	
	for(SparseMatrix::Row::const_iterator it=row.begin();it!=row.end();it++)
	  {
	    int c = (*it).first+1;
	    double val = (*it).second;
	    sum += val*rm(c);
	  }

	ret(j) = sum;
      }
  }

  void multiply(const DiagonalMatrix& lm, const SparseMatrix& rm, SparseMatrix& ret)
  {
    Tracer_Plus tr("SparseMatrix::multiply");

    int nrows = lm.Nrows();
    int ncols = rm.Ncols();

    if(lm.Ncols() != rm.Nrows()) throw Exception("Rows and cols don't match in SparseMatrix::multiply");

    ret.ReSize(nrows,ncols);

    for(int j = 1; j<=nrows; j++)
      {
	const SparseMatrix::Row& row = rm.row(j);	
	for(SparseMatrix::Row::const_iterator it=row.begin();it!=row.end();it++)
	  {
	    int c = (*it).first+1;
	    double val = (*it).second;
	    ret.insert(j,c,val*lm(j,j));
	  }
      }

  }

  void multiply(const SparseMatrix& lm, const SparseMatrix::Row& rm, ColumnVector& ret)
  {
    Tracer_Plus tr("SparseMatrix::multiply3");

    int nrows = lm.Nrows();   
       
    ret.ReSize(nrows);

    for(int j = 1; j<=nrows; j++)
      {
	float sum = 0.0;
	const SparseMatrix::Row& row = lm.row(j);

	SparseMatrix::Row::const_iterator it=row.begin();
	SparseMatrix::Row::const_iterator itrm=rm.begin();

	while(it!=row.end() && itrm!=rm.end())
	  {
	    int crm = (*itrm).first;
	    int c = (*it).first;
	    if(c==crm)
	      {
		sum += ((*itrm).second)*((*it).second);
		it++;
		itrm++;
	      }
	    else if(c < crm)
	      {
		it++;
	      }
	    else
	      {
		itrm++;
	      }
	  }

	ret(j) = sum;
      }
  }

  void add(const SparseMatrix& lm, const SparseMatrix& rm, SparseMatrix& ret)
  {
    Tracer_Plus tr("SparseMatrix::add");

    int nrows = lm.Nrows();
    int ncols = lm.Ncols();

    if(lm.Ncols() != rm.Ncols() || lm.Nrows() != rm.Nrows()) throw Exception("Rows and cols don't match in SparseMatrix::add"); 

    ret.ReSize(nrows,ncols);

    for(int j = 1; j<=nrows; j++)
      {
	const SparseMatrix::Row& lmrow = lm.row(j);
	const SparseMatrix::Row& rmrow = rm.row(j);

	SparseMatrix::Row::const_iterator lmit = lmrow.begin();
	SparseMatrix::Row::const_iterator rmit = rmrow.begin();
	int lmc = (*lmit).first+1;
	int rmc = (*rmit).first+1;

	while(lmit!=lmrow.end() || rmit!=rmrow.end())
	  {
	    if((lmc<rmc && lmit!=lmrow.end()) || rmit==rmrow.end())
	      {		
		ret.insert(j,lmc,(*lmit).second+rm(j,lmc));
		lmit++;
		lmc = (*lmit).first+1;
	      }
	    else if((rmc<lmc && rmit!=rmrow.end()) || lmit==lmrow.end())
	      {
		ret.insert(j,rmc,(*rmit).second+lm(j,rmc));
		rmit++;
		rmc = (*rmit).first+1;
	      }
	    else
	      {
		//lmc==rmc
		ret.insert(j,rmc,(*lmit).second+(*rmit).second);
		lmit++;
		lmc = (*lmit).first+1;
		rmit++;
		rmc = (*rmit).first+1;
	      }
	  }
      }
  }

  // Concatenation. Note that these work also for concatenating sparse and full (newmat) 
  // matrices by virtue of the Matrix->SparseMatrix constructor/converter.

  // ret = [A; B]; % Matlab lingo
  void vertconcat(const SparseMatrix& A, const SparseMatrix& B, SparseMatrix& ret)
  {
    if (A.Ncols() != B.Ncols()) {throw Exception("Cols don't match in SparseMatrix::vertconcat");}     

    ret.ReSize(A.Nrows()+B.Nrows(),A.Ncols());

    for (int i=1; i<=A.Nrows(); i++) {ret.row(i) = A.row(i);}
    for (int i=1; i<=B.Nrows(); i++) {ret.row(i+A.Nrows()) = B.row(i);}
  }

  // ret = [A B]; % Matlab lingo
  void horconcat(const SparseMatrix& A, const SparseMatrix& B, SparseMatrix& ret)
  {
    if (A.Nrows() != B.Nrows()) {throw Exception("Rows don't match in SparseMatrix::horconcat");}     

    ret.ReSize(A.Nrows(),A.Ncols()+B.Ncols());

    for (int i=1; i<=A.Nrows(); i++) {
      ret.row(i) = A.row(i);
      const SparseMatrix::Row& tmpRow = B.row(i);
      for (SparseMatrix::Row::const_iterator it=tmpRow.begin(); it!=tmpRow.end(); it++) {
        ret.insert(i,A.Ncols()+int(it->first)+1,double(it->second));
      }
    }         
  }

}










