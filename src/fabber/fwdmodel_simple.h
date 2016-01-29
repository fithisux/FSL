/*  fwdmodel_simple.h - Implements the simplified ASL model

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

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

#include "fwdmodel.h"
#include "inference.h"
#include <string>
using namespace std;

class SimpleFwdModelIdStruct {
public:
  SimpleFwdModelIdStruct() 
        { Nbasis = 0; }
  void Define( int NumBasisFcns ) 
        { assert(Nbasis==0); assert(NumBasisFcns>0); Nbasis = NumBasisFcns; }

//  virtual void DumpVector(const ColumnVector& vec, const string& indent = "");
//  virtual ~SimpleFwdModelIdStruct() { return; }
  
// N=3: 1=Q0 2=Q1 3=Q2 4=Q3 5=U0 6=U1 7=U2 8=U3 9=R0 10=R1 11=R2 12=R3
  double Q0(const ColumnVector& p) const
        { return p(1); }
  ReturnMatrix Qn(const ColumnVector& p) const
        { AssertValid(p); ColumnVector tmp = p.Rows(2,Nbasis+1); return tmp; }
  double U0(const ColumnVector& p) const
	{ return p(2+Nbasis); }
  ReturnMatrix Un(const ColumnVector& p) const
	{ AssertValid(p); ColumnVector tmp = p.Rows(3+Nbasis,2+2*Nbasis); return tmp; }
  double R0(const ColumnVector& p) const
	{ AssertValid(p); return p(3+2*Nbasis); }
  ReturnMatrix Rn(const ColumnVector& p) const
	{ AssertValid(p); ColumnVector tmp = p.Rows(4+2*Nbasis,3+3*Nbasis); return tmp; }
  void NameParams(vector<string>& names) const;

private:
  int Nbasis;
  void AssertValid(const ColumnVector& p) const
    { assert(Nbasis>0); assert(p.Nrows() == 3 + 3*Nbasis); }
};


class SimpleFwdModel : public FwdModel {
public: 
  // Virtual function overrides
  virtual void Evaluate(const ColumnVector& params, 
			      ColumnVector& result) const;
                  
  virtual void DumpParameters(const ColumnVector& vec,
                                const string& indents = "") const;
  virtual void NameParams(vector<string>& names) const
    { id.NameParams(names); }                              
  virtual int NumParams() const { return 3 + 3*basis.Ncols(); }
  string ModelVersion() const;

  virtual ~SimpleFwdModel() { return; }

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const
    {
      prior.means = 0;
      prior.SetPrecisions(IdentityMatrix(NumParams()) * 1e-12);
      posterior = prior;
      // Informative starting points
      // I'm not very happy with this whole setup... the physical arrangment
      // of elements in the vector should only appear once.
      posterior.means(1) = 50; 
      posterior.means(2+basis.Ncols()) = 1.5e4; 
      posterior.means(3+2*basis.Ncols()) = 25;
      assert(id.Q0(posterior.means) == 50);
      assert(id.U0(posterior.means) == 1.5e4);
      assert(id.R0(posterior.means) == 25);
    }



  // Constructor
  SimpleFwdModel(ArgsType& args);

  //  SimpleFwdModel(const SimpleFwdModel& from); // copy constructor - default ok?

protected:
  // Constants
  ColumnVector echoTime;
  Matrix basis; // basis in columns?
  ColumnVector rho;
  SimpleFwdModelIdStruct id;
};

