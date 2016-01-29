/*    Copyright (C) 2012 University of Oxford  */

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

#ifndef _MVdisc
#define _MVdisc

#include <iostream>

#include <string>
#include <fstream>
#include <stdio.h>
#include <algorithm>

#include "newimage/newimageall.h"
#include "meshclass/meshclass.h"

//using namespace std;
//using namespace NEWIMAGE;
//using namespace mesh;


namespace mvdisc{
		class mvdiscException : public std::exception{
		
public:
		const char* errmesg;
		mvdiscException(const char* msg)
		{
			errmesg=msg;
		}
		
private:
			virtual const char* what() const throw()
		{
				return errmesg;
		}
	};

	
	
	class MVdisc {
		
public:
		MVdisc();
		~MVdisc();
		
		vector<unsigned int> applyLDA(const vector<ColumnVector> & Data, const float & eigThresh) const ;
		short applyLDA(ColumnVector & Data, const float & eigThresh) const; 
		ColumnVector run_LOO_LDA(const NEWMAT::Matrix & Data, const ColumnVector & target);

		void estimateLDAParams(const NEWMAT::Matrix & Data, const ColumnVector & target) ;
				void estimateAndAppendLDAParams(const NEWMAT::Matrix & Data, const ColumnVector & target) ;

		ReturnMatrix getGroupMeans(const Matrix & Data, const ColumnVector & target){  return calculateClassMeans(Data, target, LDAnsubs); }	
						
		//---------------------------------I/O FUNCTION-----------------------------------------//

		void saveLDAParams(const string & outname ) const;
				void saveLDAParams(const string & outname, const Matrix & polygons ) const;

		//---------------------------------SETTING PARAMETERS--------------------------------------//
		void set_LDA_Params(const NEWMAT::Matrix & means, const NEWMAT::Matrix & cov_vecs,  const vector<float> & cov_eigs,const vector<unsigned int> n) 
				{ LDAmu=means ; LDAcov_Vecs=cov_vecs; LDAcov_Eigs=cov_eigs; LDAnsubs=n; };
		
		const NEWMAT::Matrix* getLDAcov_Vecs_ptr(){ return &LDAcov_Vecs; }
		vector<float>::const_iterator getLDAcov_Eigs_iter(){ return LDAcov_Eigs.begin(); }

		
		
private:
		
		//------------------------------LDA/QDA parameters----------------------------//
		Matrix LDAmu; //A column per group
		Matrix LDAcov_Vecs;//can be shared by with QDA
		vector<float> LDAcov_Eigs;//can be shared by with QDA
		vector<unsigned int> LDAnsubs;
	

		void quickSVD(const Matrix & data,  DiagonalMatrix &D,Matrix &U, Matrix &V ) const ;
		void quickSVD(const Matrix & data,  DiagonalMatrix &D,Matrix &U ) const ;
		
		ReturnMatrix calculateClassMeans(const Matrix & Data, const ColumnVector & target, vector<unsigned int> & nk) const ;
		ReturnMatrix sortAndDemeanByClass(const Matrix & Data, const ColumnVector & target, const Matrix & muM, const vector<unsigned int> & nk,ColumnVector & targSorted) const;
		void estimateCommonCov(const Matrix & DeMean, const vector<unsigned int> & nK, Matrix & U, vector<float> & D) const;
		void estimateClassCovs(const Matrix & DeMean, const ColumnVector & target, const vector< int > & nK, vector< Matrix > & vU, vector< DiagonalMatrix > & vD) const;
		
		
		template<class T>
		vector<T> threshInvert(const vector<T> & D, const T & p) const ;	
		//discrimiant function
			
	};
}
#endif
