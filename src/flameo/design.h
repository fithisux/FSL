/*  design.h

    Mark Woolrich, Tim Behrens and Matthew Webster - FMRIB Image Analysis Group

    Copyright (C) 2009 University of Oxford  */

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

#if !defined(design_h)
#define design_h
  
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "newmatap.h"
#include "newmatio.h"
#include "newimage/newimageall.h"

using namespace NEWIMAGE;
using namespace NEWMAT;
using namespace MISCMATHS;

namespace Gs{

  class Design
    {
    public:

      // constructor
      Design() :	
	nevs(0),
	ntpts(0),
	ngs(0),
	voxelwise_dm(false)
	{ 
	}

      // load design matrix in from file and set it up
      void setup(bool loadcontrasts = true);

      // getters
      const int getnevs() const { return nevs; }
      const int getntpts() const { return ntpts; }
      const int getngs() const { return ngs; }
      const int dof() const { return int(ols_dof(dm)); }
      void setvarcopedata(float val){varcopedata=val;}

     // returns full, global design matrix
      const Matrix& getdm() const { return dm; }
 
      // returns voxelwise design matrix with any voxelwise zero evs removed
      ReturnMatrix getdm(int x, int y, int z) const;

      // returns voxelwise design matrix for group g with any voxelwise zero evs removed
      ReturnMatrix getdm(int x, int y, int z, int g) const;

      const Matrix& getcovsplit() const { return covsplit; }
      
      ReturnMatrix gettcontrast(const int p_num) const { RowVector ret = tcontrasts.Row(p_num); ret.Release(); return ret; }

      const Matrix& getfcontrast(const int p_num) const { return fc[p_num-1]; }

      int getgroup(int t) const { return int(group_index(t)); } 
      int getglobalindex(int g, int within_t) const { return global_index[g-1][within_t-1]; } 
      int getindexingroup(int g, int global_t) const { return index_in_group[g-1][global_t-1]; } 

      int getntptsingroup(int g) const { return int(ntptsing(g)); }
      int getnevsingroup(int g) const { return int(nevsing(g)); }
 
      bool tcontrast_has_zeroevs(int x, int y, int z,const RowVector& tcontrast)
      {
// 	cout<< "evs_group" << evs_group<< endl;
// 	cout<< "tcontrast" << tcontrast<< endl;
	bool in=false;
	for(int e=1; e<=tcontrast.Ncols() && !in; e++)
	  {	    
	    in = (tcontrast(e)!=0 && zero_evs(x,y,z,e-1));
	  }
	return in;
      }

      bool is_group_in_tcontrast(int g, const RowVector& tcontrast)
      {
// 	cout<< "evs_group" << evs_group<< endl;
// 	cout<< "tcontrast" << tcontrast<< endl;
	bool in=false;
	for(int e=1; e<=tcontrast.Ncols() && !in; e++)
	  {	    
	    in = (tcontrast(e)!=0 && evs_group(e)==g);
	  }
	return in;
      }

     bool fcontrast_has_zeroevs(int x, int y, int z,const Matrix& fcontrast)
      {
	bool in=false;
	for(int f=1; f<=fcontrast.Nrows() && !in; f++)
	  {
	    for(int e=1; e<=fcontrast.Ncols() && !in; e++)
	      {	    
		in = (fcontrast(f,e)!=0 && zero_evs(x,y,z,e-1));
	      }	    	    
	  }
	return in;
      }

      bool is_group_in_fcontrast(int g, const Matrix& fcontrast)
      {
	bool in=false;
	for(int f=1; f<=fcontrast.Nrows() && !in; f++)
	  {
	    for(int e=1; e<=fcontrast.Ncols() && !in; e++)
	      {	    
		in = (fcontrast(f,e)!=0 && evs_group(e)==g);
	      }	    	    
	  }
	return in;
      }
      bool is_voxelwise_dm() const {return voxelwise_dm;}

      // indexed by ev number. Returns group that that ev belongs to
      int get_evs_group(int e, int x, int y, int z) {
	int ret = int(evs_group(e));
	if(zero_evs(x,y,z,e-1))
	  ret=-1;
	return ret;
      }

      int getnumfcontrasts() const { return numFcontrasts; }
      int getnumtcontrasts() const { return numTcontrasts; }
      //const Matrix& getfreduceddm(const int p_num) const { return reduceddms[p_num-1]; }
      
      const volume4D<float>& getcopedata() const {return copedata;}
      const volume4D<float>& getvarcopedata() const {return varcopedata;}
      const volume4D<float>& getdofvarcopedata() const {return dofvarcopedata;}
      const volume<float>& getsumdofvarcopedata() const {return sum_dofvarcopedata;}
      const volume<float>& getmask() const {return mask;}

      ReturnMatrix remove_zeroev_pes(int x, int y, int z, const ColumnVector& petmp);
      ReturnMatrix remove_zeroev_covpes(int x, int y, int z, const SymmetricMatrix& covpetmp);
      ReturnMatrix insert_zeroev_pemcmcsamples(int x, int y, int z, const Matrix& mcmcin);
      ReturnMatrix insert_zeroev_pes(int x, int y, int z, const ColumnVector& petmp);
      ReturnMatrix insert_zeroev_covpes(int x, int y, int z, const SymmetricMatrix& covpetmp);

      ReturnMatrix getcopedata(int x, int y, int z, int g);
      ReturnMatrix getvarcopedata(int x, int y, int z, int g);

    private:

      const Design& operator=(Design& par);     
      Design(Design& des) { operator=(des); }
      
      void setupfcontrasts();

      int nevs;
      int ntpts;
      int ngs;

      Matrix dm;
      Matrix covsplit;

      vector<volume4D<float> > voxelwise_evs;
      vector<int> voxelwise_ev_numbers;

      Matrix tcontrasts;
      Matrix fcontrasts;
      vector<Matrix> fc;

      //      vector<Matrix> reduceddms;

      ColumnVector group_index; // stores group for given global index
      vector<vector<int> > global_index; // stores global index for given group and within-group index

      ColumnVector ntptsing;
      ColumnVector nevsing;
      vector<vector<int> > index_in_group; // stores within-group index for given group and global index

      int numFcontrasts;
      int numTcontrasts;

      // indexed by ev number. Returns group that ev belongs to
      ColumnVector evs_group;

      bool voxelwise_dm;

      // inputs
      volume4D<float> copedata;
      volume4D<float> varcopedata;
      volume4D<float> dofvarcopedata;
      volume<float> sum_dofvarcopedata;
      volume<float> mask;

      volume4D<int> zero_evs;

    };
} 
#endif


