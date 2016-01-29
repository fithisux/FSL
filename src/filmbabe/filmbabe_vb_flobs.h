/*  filmbabe_vb_flobs.h

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

#if !defined(filmbabe_vb_flobs_h)
#define filmbabe_vb_flobs_h

#include <iostream>
#include <fstream>
#include "newmatap.h"
#include "newmatio.h"
#include "newimage/newimageall.h"
#include "utils/tracer_plus.h"
#include "filmbabeoptions.h"
#include "connected_offset.h"
#include "miscmaths/sparse_matrix.h"

using namespace NEWMAT;
using namespace NEWIMAGE;

namespace Filmbabe {

  class Voxel
    {
    public:
      int x;
      int y;
      int z;

      Voxel(int px, int py, int pz) : x(px),y(py),z(pz){}
	
      bool operator==(const Voxel& r) const 
      {return (x==r.x && y==r.y && z==r.z);}
	
      bool is_neighbour(const Voxel& r) const 
	{return ((x==r.x-1 && y==r.y && z==r.z)||
		 (x==r.x+1 && y==r.y && z==r.z)||
		 (x==r.x && y==r.y-1 && z==r.z)||
		 (x==r.x && y==r.y+1 && z==r.z)||
		 (x==r.x && y==r.y && z==r.z-1)||
		 (x==r.x && y==r.y && z==r.z+1)) ;}
    };

  class Filmbabe_Vb_Flobs
    {
    public:
	
      // Constructor
      Filmbabe_Vb_Flobs(const volume4D<float>& pdata, const volume<int>& pmask, const Matrix& pdesignmatrix, const ColumnVector& pflobsregressors, const volume4D<float>& plocalweights, const vector<Connected_Offset>& pconnected_offsets, int pnum_superthreshold);
	
      // setup
      void setup();
	
      // run
      void run();
	
      // save data to logger dir
      void save() ;
	
      // Destructor
      virtual ~Filmbabe_Vb_Flobs(){}
	
    private:
	
      Filmbabe_Vb_Flobs();
      const Filmbabe_Vb_Flobs& operator=(Filmbabe_Vb_Flobs&);     
      Filmbabe_Vb_Flobs(Filmbabe_Vb_Flobs&);

      void process_flobsregressors();
      void update_Beta();
      void update_A();
      void update_phiA();
      void update_phie();
      void update_phiBeta();

      int xsize; 
      int ysize;
      int zsize;
      int ntpts;

      // num flobs Original evs
      int nflobsevs;

      // num basis functions for each flobs Original evs (is the same for all Original evs)
      int nbfs;

      // num of non flobs evs
      int nnonflobsevs;

      // total number of actual real evs in design matrix
      int nrealevs;

      const volume4D<float>& data;
      const volume<int>& mask;
      Matrix designmatrix;
      const ColumnVector& flobsregressors;

      Matrix QtXXt;

      const volume4D<float>& localweights;
      const vector<Connected_Offset>& connected_offsets;

      volume<int> indices;
      vector<Voxel> voxels;

      Matrix Y;

      // Prior distribution params
      // m_Beta_0(e)
      ColumnVector m_Beta_0_global;
      SymmetricMatrix lambda_Beta_0_global;

      vector<ColumnVector> m_Beta_0;
      vector<SymmetricMatrix> lambda_Beta_0;
      
      SparseMatrix D;

      ColumnVector trace_ilambdaDA;
      vector<ColumnVector> trace_ilambdaBeta;
      ColumnVector trace_ilambdaXXt;
      vector<ColumnVector> diag_ilambdaA;

      // Posterior distribution params

      // q(Beta)~MVN(m_Beta,lambda_Beta);
      // m_Beta[i]((bf-1)*nevs+e)
      vector<ColumnVector> m_Beta;
      vector<SymmetricMatrix> ilambda_Beta;

      // q(A_p)~MVN(m_A_p,lambda_A_p);
      // m_A[p](i)
      vector<ColumnVector> m_A;
      
      // q(phi_Beta)~Ga(b_Beta,c_Beta);
      // gam_Beta = gam(c_Beta+1)/gam(c_Beta)
      // gam_Beta[i](e)
      vector<ColumnVector> gam_Beta;

      // q(phi_A_p)~Ga(b_A_p,c_A_p);
      // gam_A_p = gam(c_A_p+1)/gam(c_A_p)
      ColumnVector gam_A;

      // q(phi_e)~Ga(b_e,c_e);
      // gam_e = gam(c_e+1)/gam(c_e)
      ColumnVector gam_e;

      int num_superthreshold;

      int maxnumneighs;
      int ntar;

      int count;

      vector<SparseMatrix> ilambdaDA;
      vector<SparseMatrix> ilambdaA;

      int niters;

      Matrix realevcoord;

      vector<Matrix> QemReCquad;
      vector<RowVector> Re;
      vector<Matrix> Qe;
      Matrix Q;
      Matrix designmatrixQ;

      vector<float> gamAhist;
    }; 
  
}
#endif
