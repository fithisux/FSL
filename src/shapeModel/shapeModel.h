/*
 *  shapeModel.h
 *  
 *
 *  Created by Brian Patenaude on 23/06/2008.
 *  Copyright 2008 University of Oxford All rights reserved.
 *
 */
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
 #ifndef SHAPEMODEL_H
 #define SHAPEMODEL_H
 #include <vector>
#include <newimage/newimageall.h>

namespace SHAPE_MODEL_NAME{


class shapeModel
{

public:
	
	std::vector< std::vector<unsigned int> > cells;
	std::vector< std::vector<unsigned int> > localTri;
	int NumberOfSubjects;
	shapeModel();
	shapeModel( const std::vector<float> & mshape, const std::vector< std::vector<float> > & modesshape, const std::vector<float> & se, \
				const std::vector<float> & ishape, const std::vector< std::vector<float> > & modesint, const std::vector<float> & ie, \
				const int & M, const std::vector<float> & errs);
				
	shapeModel( const std::vector<float> & mshape, const std::vector< std::vector<float> > & modesshape, const std::vector<float> & se, \
				const std::vector<float> & ishape, const std::vector< std::vector<float> > & modesint, const std::vector< std::vector<float> > & Iprec, const std::vector<float> & ie, \
				const int & M, const std::vector<float> & errs, const std::vector< std::vector<unsigned int> > & cellsin, const std::vector<int> & vlabels);
	
	shapeModel( const std::vector<float> &  mshape, const std::vector< std::vector<float> > & modesshape, const std::vector<float> & se, \
				const std::vector<float> & ishape, const std::vector< std::vector<float> > & modesint, const std::vector<float> & ie, \
				const int & M, const std::vector<float> & errs, const std::vector<short> & vmaskin);
	
	std::vector<float> getDeformedGrid( const std::vector<float>  & vars ) const ;
	
	std::vector<float> getDeformedIGrid( const std::vector<float> & vars) const ;
	 void printLabel(const unsigned int & i) const ;
	int getLabel( const unsigned int & i) const { return labels.at(i); }
	
	void registerModel(const std::vector< std::vector<float> > & flirtmat);
	std::vector<float> getOrigSpaceBvars(const std::vector<float> & bvars ) const;
	
	std::vector< std::vector<float> > registerModeVectors( const std::vector< std::vector<float> >& vmodes, const std::vector< std::vector<float> >& flirtmat);

	
	void setFoundMode(bool found) const { MODE_FOUND=found; }
	bool getFoundMode() const { return MODE_FOUND; } 
	void setMode(float m) const { mode=m; }
	float getMode() const { return mode; }
	bool getCondSet() const { return USE_COND; }
	void setCondSet( const bool & b ) const { USE_COND=b; }
	void setCondMats( const std::vector< std::vector<float> > & mat1, const std::vector< std::vector<float> > & mat2 , const unsigned int & k ){ kpred = k; condMat1 =mat1 ; condMat2=mat2; }
	void setWarpField( const NEWIMAGE::volume4D<float> & w ) { warpField=w; }

	
 std::vector< std::vector<float> > getCondMat1() const { return condMat1; }	
		std::vector< std::vector<float> > getCondMat2() const { return condMat2; }	
unsigned int getKPred() const { return kpred; }

	std::vector<float> smean;
	std::vector< std::vector<float> > smodes;
	std::vector< std::vector<float> > u_xfm;
	std::vector<float> d_xfm;
	
	std::vector< float > seigs;
	std::vector< float > sqrtseigs;

	std::vector< float > sqrtseigsi;

	std::vector< int > labels;

	std::vector<float> imean;
	std::vector< std::vector<float> > imodes;
	std::vector< float > ieigs;
	std::vector< std::vector<float> > i_precision;
	std::vector< float > Errs;
	std::vector<short> stmask;
	
	
	NEWIMAGE::volume4D<float> warpField;
	static shapeModel* loadAndCreateShapeModel( const string & modelname, const bool & verbose);

	
	private:
		
		mutable std::vector< std::vector<float> > condMat1;
		mutable std::vector< std::vector<float> >  condMat2;
		mutable unsigned int kpred;
		mutable bool USE_COND;
	mutable bool STORE_REG_XFM;

		mutable bool MODE_FOUND;
		mutable float mode;
	
};

}

#endif

