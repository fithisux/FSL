/*
 *  first_newmat_vec.h
 *  
 *
 *  Created by Brian Patenaude on 12/08/2008.
 *  Copyright 2008 University of Oxford. All rights reserved.
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

#ifndef NEWMAT_VEC_H
#define NEWMAT_VEC_H

#include <vector>
#include <newmat.h>

namespace FIRST_LIB{
	
	class first_newmat_vector{
		public:
			
	first_newmat_vector();
	~first_newmat_vector();
	
	
template<class T>
static std::vector<T> vectorToVector( const NEWMAT::Matrix & sm, const int & MaxModes);

template<class T>
static std::vector<T> vectorToVector( const NEWMAT::Matrix & sm);

template<class T>
static NEWMAT::ReturnMatrix vectorOfVectorsToMatrix( const std::vector< std::vector<T> > & vec);

template<class T>
static NEWMAT::ReturnMatrix vectorToDiagonalMatrix( const std::vector<T> & vec);

 template<class T>
   static std::vector<T> DiagonalMatrixToVector( const  NEWMAT::DiagonalMatrix & vec);


template<class T>
static std::vector< std::vector<T> > matrixToVector( const NEWMAT::Matrix & sm, const int & MaxModes);

template<class T>
static std::vector< std::vector<T> > matrixToVector( const NEWMAT::Matrix & sm);

template<class T>
static std::vector< std::vector<T> > matrixToVectorOfVectors(const NEWMAT::Matrix & m);

inline
static NEWMAT::ReturnMatrix unwrapMatrix(const NEWMAT::Matrix & m) 
{
	NEWMAT::ColumnVector munwrap(m.Nrows()*m.Ncols());
	unsigned int count=0;
	for (int i =0; i<m.Nrows() ; i++)
		for (int j =0; j<m.Ncols() ; j++,count++)
			munwrap.element(count)=m.element(i,j);
	return munwrap;
}

inline
static NEWMAT::ReturnMatrix wrapMatrix(const NEWMAT::ColumnVector & m) 
{
	
	NEWMAT::Matrix wrapped(m.Nrows()/3,3);
	unsigned int count=0;
	for (int i =0; i<m.Nrows() ; i+=3,count++)
	{
		wrapped.element(count,0)=m.element(i);
		wrapped.element(count,1)=m.element(i+1);
		wrapped.element(count,2)=m.element(i+2);
	}
	return wrapped;
}


	
	};
	
}

#endif
