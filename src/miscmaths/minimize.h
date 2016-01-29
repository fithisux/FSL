/*  minimize
 
    Tim Behrens, FMRIB Image Analysis Group

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

#if !defined(minimize_h)
#define minimize_h

#include <string>
#include <iostream>
#include <fstream>
//#include <unistd.h>
#include <vector>
#include <algorithm>
#include "newmatap.h"
#include "newmatio.h"
#include "miscmaths.h"


#define WANT_STREAM
#define WANT_MATH

using namespace MISCMATHS;
using namespace NEWMAT;
using namespace std;
///////////////////////////////////////////////////////

//fminsearch.m

namespace MISCMATHS {

class pair_comparer
{
public:
  bool operator()(const pair<float,ColumnVector>& p1,const pair<float,ColumnVector>& p2) const
  {
    return p1.first < p2.first;
  }    
};

  class EvalFunction;
  class gEvalFunction;

float diff1(const ColumnVector& x, const EvalFunction& func, int i,float h,int errorord=4);// finite diff derivative

float diff2(const ColumnVector& x, const EvalFunction& func, int i,float h,int errorord=4);// finite diff 2nd derivative

float diff2(const ColumnVector& x, const EvalFunction& func, int i,int j,float h,int errorord=4);// finite diff cross derivative

ReturnMatrix gradient(const ColumnVector& x, const EvalFunction& func,float h,int errorord=4);// finite diff derivative vector 

ReturnMatrix hessian(const ColumnVector& x, const EvalFunction& func,float h,int errorord=4);// finite diff hessian

void minsearch(ColumnVector& x, const EvalFunction& func, ColumnVector& paramstovary);

void scg(ColumnVector& x, const gEvalFunction& func, ColumnVector& paramstovary, float tol = 0.0000001, float eps=1e-16, int niters=500);

class EvalFunction
{//Function where gradient is not analytic (or you are too lazy to work it out) (required for fminsearch)
public:
  EvalFunction(){}
  virtual float evaluate(const ColumnVector& x) const = 0; //evaluate the function
  virtual ~EvalFunction(){};

  virtual void minimize(ColumnVector& x)
  {
    ColumnVector paramstovary(x.Nrows());
    paramstovary = 1;
    minsearch(x,*this,paramstovary);
  }

  virtual void minimize(ColumnVector& x, ColumnVector& paramstovary)
  {
    minsearch(x,*this,paramstovary);
  }

private:
  const EvalFunction& operator=(EvalFunction& par);
  EvalFunction(const EvalFunction&);
};

class gEvalFunction : public EvalFunction
{//Function where gradient is analytic (required for scg)
public:
  gEvalFunction() : EvalFunction(){}
  // evaluate is inherited from EvalFunction
  
  virtual ReturnMatrix g_evaluate(const ColumnVector& x) const = 0; //evaluate the gradient
  virtual ~gEvalFunction(){};

  virtual void minimize(ColumnVector& x)
  {
    ColumnVector paramstovary(x.Nrows());
    paramstovary = 1;
    scg(x,*this,paramstovary);
  }

  virtual void minimize(ColumnVector& x, ColumnVector& paramstovary)
  {
    scg(x,*this,paramstovary);
  }

private:
    
  const gEvalFunction& operator=(gEvalFunction& par);
  gEvalFunction(const gEvalFunction&);
};

}
   

#endif







