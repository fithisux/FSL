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


#ifndef _fslvtkIO
#define _fslvtkIO

#include <iostream>

#include <string>
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include "meshclass/meshclass.h"
#include "first_lib/first_newmat_vec.h"

//using namespace std;
//using namespace mesh;


namespace fslvtkio{
		//using namespace mesh;
	
	class fslvtkIOException : public std::exception{
		
public:
		const char* errmesg;
		fslvtkIOException(const char* msg)
		{
			errmesg=msg;
		}
		
private:
			virtual const char* what() const throw()
		{
				cout<<errmesg<<endl;
				return errmesg;
		}
	};
	
	class fslvtkIO {
public:
		enum DataType{ POLYDATA, UNSTRUCTURED_GRID };
		fslvtkIO();
		fslvtkIO(const string & filename,const DataType i);
		
		~fslvtkIO();

		
		//----------------------SET COMMANDS-------------------------//
		void setPoints(const Matrix & pts);
		void setPoints(const vector<float> & pts);
		void setCells(const vector< vector<unsigned int> > & c){ Cells=c; }
		void setCells(const vector< vector<unsigned int> > & c, const vector<short> & c_t ){ Cells=c; Cell_Types=c_t; }
		
		void appendPointsAndPolygons(const Matrix & pts, const Matrix & Polygons);


		void setPolygons(const Matrix& m){ Polygons=m; }
		void setPolygons(const vector< vector<unsigned int> >& vm);
		void setMesh(const mesh::Mesh & m);
		//void setScalars(const Matrix & sc,const string & name);
		void setScalars(const Matrix& m){ Scalars=m; }
		
		void setVectors(const Matrix & vecs,const string &name);
		void setVectors(const Matrix& m){ Vectors=m; }

		template<class T>
		void setScalars(const vector<T> & sc);
		
			
		void addPointFieldData(const Matrix& m,const string & name, const string & type, const string & vtkAttType);			
		void addCellFieldData(const Matrix& m,const string & name, const string & type, const string & vtkAttType);
		
                void addFieldData(const Matrix& m,const string & name, const string & type);
                void addFieldData(const ReturnMatrix& m,const string & name, const string & type);

		template< class T > 
		void addFieldData(const vector<T> & m,const string & name, const string & type);
		
		template< class T >
		void addFieldData(const T & m,const string & name, const string & type);

		void replaceFieldData(const Matrix& m,const string & name);

		void addFieldData(vector< string >,string name);
		
		
		void printFieldDataNames() { for (vector<string>::iterator i=fieldDataNumName.begin(); i!=fieldDataNumName.end();i++) cout<<"field "<<*i<<endl; }

		
		//----------------------GET COMMANDS----------------------------//
		Matrix getPointsAsMatrix() const { return Points; }
		ColumnVector getPointsAsColumnVector() const ;



		template<class T>
		std::vector<T> getPointsAsVector();
		
		template<class T>
		std::vector< std::vector<T> > getPointsAsVectorOfVectors(){ return FIRST_LIB::first_newmat_vector::matrixToVectorOfVectors<T>(Points); }
			
		template<class T>
		std::vector< std::vector<T> > getPolygonsAsVectorOfVectors(){ return FIRST_LIB::first_newmat_vector::matrixToVectorOfVectors<T>(Polygons); }

		ReturnMatrix getPolygons() const { return Polygons; };
		vector< vector<unsigned int> > getCells() const { return Cells; }
		vector<short> getCellTypes() const { return Cell_Types; }
								
		string getFieldName(const int & ind) const { return fieldDataNumName.at(ind); };
		unsigned int getNumberOfFields() const { return fieldDataNumName.size(); };
		Matrix getField(const string & name); 
		Matrix getField(const string & name, unsigned int & ind); 
		void setField(const string & name, const Matrix & data);

		ReturnMatrix getScalars(){ return Scalars; }

		template<class T>
		  vector<T> getScalars();


		ReturnMatrix getVectors(){ return Vectors; }
		
		
		//----------------------I/O COMMANDS----------------------------//
		template<class T>
			void writePoints(ofstream & fshape, const string & str_typename);
		void writePolygons(ofstream & fshape);
		void writeCells(ofstream & fshape);
		//void writeUnstructuredGridCells(ofstream & fshape);
		void writeUnstructuredGridCellTypes(ofstream & fshape);
		void save(string s);
		
		void readUnstructuredGrid(string fname);
		void readPolyData(string fname);
		template<class T>
			void writePointData(ofstream & fshape, const string & str_typename );
		
		template <class T>
			void writeNumericField(ofstream & fvtk, const string & name, const string & type, const Matrix & Data);
		void writeStringField(ofstream & fvtk, const string & name, const vector<string> & v_string);
		
		void readFieldData(ifstream & fvtk);
		void readPointData(ifstream & fvtk, string & nextData);
		bool  readPoints(ifstream & fvtk);
		bool  readPolygons(ifstream & fvtk);
		
		template <class T>
			ReturnMatrix readField(ifstream & fvtk, const int & nrows,const int & mcols);
		
		void displayNumericFieldDataNames(); 				
				void displayNumericField(const string & name);

		//----------------------SETTING/ACCESS OF STATE VARIABLES----------------------------//
		void setSwitchRowsCols(bool n) { SWITCH_ROWS_COLS=n; }

		void setMAX(bool b) { MAX_SET=b; }
		void setMAX_Val( unsigned int n ) { MAX=n; }
		void setDataType(const DataType & dtype ){ dt=dtype; }
		void setBinaryWrite(bool state) { BINARY=state; }
		bool getBinaryWrite() { return BINARY; }
		static bool myExceptions(int e);
		
protected:
		Matrix Scalars;
		Matrix Vectors;
		Matrix Points;
		//cell data
		Matrix Polygons;
		
		
		
private:
#define BINFLAG 42
			//----------------------STATE VARIABLES/CONSTANTS----------------------------//
				
#ifdef PPC64
	int m_n;
#endif
			bool BINARY;//state variable for read write
			bool SWAP_BYTES;//only used for binary read
				bool MAX_SET; //is a max number of columns imposed
				bool SWITCH_ROWS_COLS;//SWITCHES THE COLUMN/ROWS , USED FOR FIXING OLD FILES

				unsigned int ST_COUNT;
	
					unsigned int MAX;
					DataType dt;//i.e. POLYDATA
						
						
						
						
						//----------------------READ COMMANDS----------------------------//
							
						//----------------------DATA----------------------------//
						
						//point data
						string scalarsName, vectorsName;
						
						vector< vector<unsigned int> > Cells;
						vector<short> Cell_Types;
						
						string cd_scalarsName, cd_vectorsName;
						Matrix cd_Scalars;
						Matrix cd_Vectors;
						
						
						//field data
						vector< Matrix > fieldDataNum;
						vector< string > fieldDataNumName;
						vector< string > fieldDataNumType;
						
						vector< vector<string> > fieldDataStr;
						vector< string > fieldDataStrName;
						
						//defaults FIELDData
						vector< string > pd_list;
						vector< string > pd_type;
						
						vector< string > cd_list;
						vector< string > cd_type;
						
	};
}
#endif
