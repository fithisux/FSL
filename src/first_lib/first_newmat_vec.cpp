/*
 *  first_newmat_vec.cpp
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

#include "first_newmat_vec.h"

using namespace std;
using namespace NEWMAT;

namespace FIRST_LIB{

first_newmat_vector::first_newmat_vector(){

}

first_newmat_vector::~first_newmat_vector(){

}

template<class T>
std::vector<T> first_newmat_vector::vectorToVector( const Matrix & sm, const int & MaxModes)
{
	std::vector<T> vecM;	
	if (sm.Nrows()==1){
		for (int i=0;i<  MaxModes ; i++){
			vecM.push_back(static_cast<T>(sm.element(0,i)));
		}
	}else{
		for (int i=0;i<  MaxModes ; i++){
			vecM.push_back(static_cast<T>(sm.element(i,0)));
		}
		
		
	}
	
	return vecM;
}

template std::vector<unsigned int> first_newmat_vector::vectorToVector<unsigned int>( const Matrix & sm, const int & MaxModes);
template std::vector<int> first_newmat_vector::vectorToVector<int>( const Matrix & sm, const int & MaxModes);
template std::vector<float> first_newmat_vector::vectorToVector<float>( const Matrix & sm, const int & MaxModes);
template std::vector<double> first_newmat_vector::vectorToVector<double>( const Matrix & sm, const int & MaxModes);

template<class T>
std::vector<T> first_newmat_vector::vectorToVector( const Matrix & sm)
{
	vector<T> vecM;	
	if (sm.Nrows()==1){
		for (int i=0;i< sm.Ncols() ; i++){
			vecM.push_back(static_cast<T>(sm.element(0,i)));
		}
	}else{
		for (int i=0;i< sm.Nrows() ; i++){
			vecM.push_back(static_cast<T>(sm.element(i,0)));
		}
		
		
	}
	
	return vecM;
}

template std::vector<unsigned int> first_newmat_vector::vectorToVector<unsigned int>( const Matrix & sm);
template std::vector<int> first_newmat_vector::vectorToVector<int>( const Matrix & sm);
template std::vector<float> first_newmat_vector::vectorToVector<float>( const Matrix & sm);
template std::vector<double> first_newmat_vector::vectorToVector<double>( const Matrix & sm);


template<class T>
ReturnMatrix first_newmat_vector::vectorOfVectorsToMatrix( const vector< vector<T> > & vec)
{
	Matrix out(vec.size(), vec.at(0).size());
	unsigned int row=0;

	for (typename vector< vector<T> >::const_iterator i=vec.begin() ; i!=vec.end();i++, row++)
	{
		unsigned int col=0;
		for (typename vector<T>::const_iterator j=i->begin() ; j!=i->end();j++,col++)
			out.element(row,col)=*j;
	}		
	out.Release();
	return out;
}

template ReturnMatrix first_newmat_vector::vectorOfVectorsToMatrix<unsigned int>( const vector< vector<unsigned int> > & vec);
template ReturnMatrix first_newmat_vector::vectorOfVectorsToMatrix<int>( const vector< vector<int> > & vec);
template ReturnMatrix first_newmat_vector::vectorOfVectorsToMatrix<float>( const vector< vector<float> > & vec);
template ReturnMatrix first_newmat_vector::vectorOfVectorsToMatrix<double>( const vector< vector<double> > & vec);

template<class T>
ReturnMatrix first_newmat_vector::vectorToDiagonalMatrix( const vector<T> & vec)
{
	DiagonalMatrix out(vec.size());
	unsigned int row=0;

	for (typename vector<T> ::const_iterator i=vec.begin() ; i!=vec.end();i++, row++)
		out.element(row)=*i;
		
	out.Release();
	return out;
}

template ReturnMatrix first_newmat_vector::vectorToDiagonalMatrix<unsigned int>( const vector<unsigned int> & vec);
template ReturnMatrix first_newmat_vector::vectorToDiagonalMatrix<int>( const vector<int> & vec);
template ReturnMatrix first_newmat_vector::vectorToDiagonalMatrix<float>( const vector<float> & vec);
template ReturnMatrix first_newmat_vector::vectorToDiagonalMatrix<double>( const vector<double> & vec);

	
	
	template<class T>
	vector<T> first_newmat_vector::DiagonalMatrixToVector(  const  NEWMAT::DiagonalMatrix & M )
	{
		vector<T> vec;
		for ( int i=0; i<M.Nrows();i++)
		{
			vec.push_back(static_cast<T>(M.element(i)));
		}
		
		return vec;
	}
	
	template vector<unsigned int> first_newmat_vector::DiagonalMatrixToVector<unsigned int>(  const  NEWMAT::DiagonalMatrix & M );
	template vector<int> first_newmat_vector::DiagonalMatrixToVector<int>(  const  NEWMAT::DiagonalMatrix & M );
	template vector<float> first_newmat_vector::DiagonalMatrixToVector<float>(  const  NEWMAT::DiagonalMatrix & M );
	template vector<double> first_newmat_vector::DiagonalMatrixToVector<double>(  const  NEWMAT::DiagonalMatrix & M );
	

template<class T>
vector< vector<T> > first_newmat_vector::matrixToVector( const Matrix & sm, const int & MaxModes)
{
	vector< vector<T> > vecM;	
    for (int j=0;j< MaxModes ; j++){
		vector<T> mode;
	    for (int i=0;i< sm.Nrows() ; i++){
			mode.push_back(sm.element(i,j));
		}
		vecM.push_back(mode);
	}   
	
	return vecM;
}

template vector< vector<unsigned int> > first_newmat_vector::matrixToVector<unsigned int>( const Matrix & sm, const int & MaxModes);
template vector< vector<int> > first_newmat_vector::matrixToVector<int>( const Matrix & sm, const int & MaxModes);
template vector< vector<float> > first_newmat_vector::matrixToVector<float>( const Matrix & sm, const int & MaxModes);
template vector< vector<double> > first_newmat_vector::matrixToVector<double>( const Matrix & sm, const int & MaxModes);


template<class T>
vector< vector<T> > first_newmat_vector::matrixToVector( const Matrix & sm)
{
	vector< vector<T> > vecM;	
    for (int j=0;j< sm.Ncols() ; j++)
	{
		vector<T> mode;
	    for (int i=0;i< sm.Nrows() ; i++)
			mode.push_back(static_cast<T>(sm.element(i,j)));
		vecM.push_back(mode);
	}   

	return vecM;
}

template vector< vector<unsigned int> > first_newmat_vector::matrixToVector<unsigned int>( const Matrix & sm);
template vector< vector<int> > first_newmat_vector::matrixToVector<int>( const Matrix & sm);
template vector< vector<float> > first_newmat_vector::matrixToVector<float>( const Matrix & sm);
template vector< vector<double> > first_newmat_vector::matrixToVector<double>( const Matrix & sm);

template<class T>
std::vector< std::vector<T> > first_newmat_vector::matrixToVectorOfVectors(const NEWMAT::Matrix & m){
	vector< vector<T> > all;
	for (int i=0;i<m.Nrows();i++)
	{
		vector<T> row;
		for (int j=0; j<m.Ncols();j++)
			row.push_back(static_cast<T>(m.element(i,j)));
		all.push_back(row);
	}
	return all;
}
template std::vector< std::vector<short> > first_newmat_vector::matrixToVectorOfVectors<short>(const NEWMAT::Matrix & m);
template std::vector< std::vector<unsigned int> > first_newmat_vector::matrixToVectorOfVectors<unsigned int>(const NEWMAT::Matrix & m);
template std::vector< std::vector<int> > first_newmat_vector::matrixToVectorOfVectors<int>(const NEWMAT::Matrix & m);
template std::vector< std::vector<float> > first_newmat_vector::matrixToVectorOfVectors<float>(const NEWMAT::Matrix & m);
template std::vector< std::vector<double> > first_newmat_vector::matrixToVectorOfVectors<double>(const NEWMAT::Matrix & m);




}
