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


#include "fslvtkio.h"
#include "meshclass/meshclass.h"
#include "newimage/newimageall.h"
#include <sstream>


using namespace std;
using namespace NEWIMAGE;
using namespace mesh;


namespace fslvtkio {
	
	fslvtkIO::fslvtkIO(){
		dt=POLYDATA;
		scalarsName="Scalars";
		vectorsName="Vectors";
		SWITCH_ROWS_COLS=false;
		SWAP_BYTES=false;
		BINARY=false;
		MAX_SET=false;
		MAX=0;
				ST_COUNT=0;		
#ifdef PPC64
	m_n=0;
#endif
	}
	
	fslvtkIO::fslvtkIO(const string & filename,const fslvtkIO::DataType i)
{
		scalarsName="Scalars";
		vectorsName="Vectors";
		SWITCH_ROWS_COLS=false;
		SWAP_BYTES=false;
		BINARY=false;
		MAX_SET=false;
		MAX=0;
		#ifdef PPC64
	m_n=0;
#endif
		switch (i)
		{
			case 0:
				setDataType(i);
				readPolyData(filename);
				ST_COUNT=1;
				break;
			case 1:
				setDataType(i);
				readUnstructuredGrid(filename);
				break;
			default:
				throw fslvtkIOException("Invalid data type. Cannot create object.");
		}
}

	fslvtkIO::~fslvtkIO(){
	}
	
	
	void fslvtkIO::setMesh(const  Mesh &  m){
		//set points
		//assume 3 dimensions
		ST_COUNT=1;
		Points.ReSize(m._points.size(),3);
		int  count=0;
		for (vector<Mpoint*>::iterator i =const_cast<Mesh &>(m)._points.begin(); i!=const_cast<Mesh &>(m)._points.end(); i++ )
		{ 
			Points.element(count,0)=(*i)->get_coord().X;	
			Points.element(count,1)=(*i)->get_coord().Y;	
			Points.element(count,2)=(*i)->get_coord().Z;	
			
			count++;
		}
		
		//set polygons, assumes eacxh vertex has 3 connections
		Polygons.ReSize(m._triangles.size(),3);
		
		count=0;
		for ( list<Triangle*>::const_iterator i=m._triangles.begin(); i!=m._triangles.end(); i++) {
			
			Polygons.element(count,0)=(*i)->get_vertice(0)->get_no();
			Polygons.element(count,1)=(*i)->get_vertice(1)->get_no();
			Polygons.element(count,2)=(*i)->get_vertice(2)->get_no();	  
			count++;
		}
		
		
	}
	
	void fslvtkIO::setPoints(const Matrix& m)
	{
		if (m.Ncols()==3)
			Points=m;
		else if ( (m.Ncols()==1) && ( (m.Nrows()%3) == 0) )
		{
			Points.ReSize(m.Nrows()/3,3);
		
			unsigned int count=0;
			for (int i =0 ; i < m.Nrows() ;i++,count++)
			{
				Points.element(count,0)=m.element(i,0);
				i++;
				Points.element(count,1)=m.element(i,0);
				i++;
				Points.element(count,2)=m.element(i,0);
			}
		}else
			throw fslvtkIOException("incompatible dimensions when setting points");
	}
	
	void fslvtkIO::setPoints(const vector<float> & m)
	{
		Points.ReSize(m.size()/3,3);
		
		unsigned int count=0;
		for (vector<float>::const_iterator i=m.begin();i!=m.end();i++,count++)
		{
			Points.element(count,0)=*i;
			i++;
			Points.element(count,1)=*i;
			i++;
			Points.element(count,2)=*i;
		}
	}
	
		void fslvtkIO::appendPointsAndPolygons(const Matrix & pts, const Matrix & polys)
		{
			cout<<"begin append"<<endl;

			//if  ( (pts.Ncols()==1) && ( (pts.Nrows()%3) == 0) && (Points.Ncols()==3) )
			//{
			//}else 
			if (pts.Ncols() != Points.Ncols())
					throw fslvtkIOException("incompatible dimensions when appending points");

			unsigned int Nprev=Points.Nrows();
		//	unsigned int N=pts.Nrows();
			Matrix Sc(pts.Nrows(),1);
			Sc=ST_COUNT;
			ST_COUNT++;
			cout<<"append points "<<endl;
			Points=Points & pts;
			cout<<"append polys "<<polys.Nrows()<<" "<<polys.Ncols()<<endl;

			Polygons=Polygons & (polys+Nprev);
			if (ST_COUNT==1)	
				Scalars=Sc;
			else
				Scalars=Scalars & Sc;
			//shift polygons by N
		//	for (unsigned int i=Nprev; i < (Nprev+N) ;i++)
		//		for (unsigned int j=0; j < static_cast<unsigned int>(Polygons.Ncols()) ; j++)
		//			Polygons.element(i,j)+=Nprev;
			cout<<"end append"<<endl;
		}

	
	template<class T>
	vector<T> fslvtkIO::getPointsAsVector()
{
		vector<T> allpoints;
		for (int i=0;i<Points.Nrows();i++)
			for (int j=0;j<Points.Ncols();j++)
				allpoints.push_back(static_cast<T>(Points.element(i,j)));
		
		return allpoints;
}

template vector<float> fslvtkIO::getPointsAsVector<float>();
template vector<double> fslvtkIO::getPointsAsVector<double>();



        ColumnVector fslvtkIO::getPointsAsColumnVector() const 
{

  ColumnVector pts(Points.Nrows() * Points.Ncols() );
	
  for (int i=0;i<Points.Nrows();i++)
    for (int j=0;j<Points.Ncols();j++)
	{
		pts.element( i*Points.Ncols() + j ) = Points.element(i,j);
	}
  return pts;
}







Matrix fslvtkIO::getField(const string & name){
	//search for field index
	
	int ind=-1;
	for (unsigned int i=0;i<fieldDataNumName.size();i++)
		if (!strcmp(fieldDataNumName.at(i).c_str(),name.c_str()))
			ind=i;
	if (ind==-1)
		throw fslvtkIOException("No field data of that name."); 
	
	return fieldDataNum.at(ind);
}
	

	
	
	void fslvtkIO::setField(const string & name, const Matrix & data){
		//search for field index
		
		int ind=-1;
		for (unsigned int i=0;i<fieldDataNumName.size();i++)
			if (!strcmp(fieldDataNumName.at(i).c_str(),name.c_str()))
				ind=i;
		if (ind==-1)
			throw fslvtkIOException("No field data of that name."); 
		
		fieldDataNum.at(ind)=data;
	}
	
	Matrix fslvtkIO::getField(const string & name, unsigned int & indout ){
		//search for field index
		
		int ind=-1;
		for (unsigned int i=0;i<fieldDataNumName.size();i++)
			if (!strcmp(fieldDataNumName.at(i).c_str(),name.c_str()))
				ind=i;
		if (ind==-1)
			throw fslvtkIOException("No field data of that name."); 
		
		indout=ind;
		
		return fieldDataNum.at(ind);
	}
	



template<class T>
vector<T> fslvtkIO::getScalars()
{
  vector<T> v_sc;
  for (int i = 0 ; i < Scalars.Nrows(); ++i)
    v_sc.push_back(Scalars.element(i,0));
  return v_sc;
}
template vector<int> fslvtkIO::getScalars();
template vector<unsigned int> fslvtkIO::getScalars();
template vector<float> fslvtkIO::getScalars();
template vector<double> fslvtkIO::getScalars();

template<class T>
void fslvtkIO::setScalars(const vector<T> & sc){
	Scalars.ReSize(sc.size(),1);
	for (unsigned int i =0; i<sc.size();i++)
		Scalars.element(i,0)=sc.at(i);
}
template void fslvtkIO::setScalars<unsigned int>(const vector<unsigned int> & sc);
template void fslvtkIO::setScalars<int>(const vector<int> & sc);
template void fslvtkIO::setScalars<float>(const vector<float> & sc);


void fslvtkIO::addPointFieldData(const Matrix & M, const string & name, const string & type, const string & vtkAttType){
	
	addFieldData(M, name, type);
	pd_list.push_back(name);
	pd_type.push_back(vtkAttType);
}


void fslvtkIO::addCellFieldData(const Matrix & M, const string & name, const string & type, const string & vtkAttType){
	addFieldData(M, name, type);
	cd_list.push_back(name);
	cd_type.push_back(vtkAttType);
}
			
void fslvtkIO::replaceFieldData(const Matrix& M,const string & name)
{
  unsigned int ind;
  getField(name, ind);
  fieldDataNum.at(ind)=M;
}
	
void fslvtkIO::addFieldData(const ReturnMatrix & M, const string & name, const string & type){
  fieldDataNum.push_back(M);
  fieldDataNumName.push_back(name);
  fieldDataNumType.push_back(type);
}

void fslvtkIO::addFieldData(const Matrix & M, const string & name, const string & type){
  fieldDataNum.push_back(M);
  fieldDataNumName.push_back(name);
  fieldDataNumType.push_back(type);
}
	
template< class T >
void fslvtkIO::addFieldData(const vector<T> & vM, const string & name, const string & type){
	ColumnVector M(vM.size());
	for (unsigned int i=0; i<vM.size();i++)
		M.element(i)=vM.at(i);
	
	fieldDataNum.push_back(M);
	fieldDataNumName.push_back(name);
	fieldDataNumType.push_back(type);
}

template void fslvtkIO::addFieldData<float>(const vector<float> & vM, const string & name, const string & type);
template void fslvtkIO::addFieldData<double>(const vector<double> & vM, const string & name, const string & type);
template void fslvtkIO::addFieldData<unsigned int>(const vector<unsigned int> & vM, const string & name, const string & type);
template void fslvtkIO::addFieldData<short>(const vector<short> & vM, const string & name, const string & type);
template void fslvtkIO::addFieldData<int>(const vector<int> & vM, const string & name, const string & type);

template< class T >
void fslvtkIO::addFieldData(const T & m,const string & name, const string & type){
  ColumnVector M(1);
  M.element(0)=m;
  fieldDataNum.push_back(M);
  fieldDataNumName.push_back(name);
  fieldDataNumType.push_back(type);
}

template void fslvtkIO::addFieldData<int>(const int & M, const string & name, const string & type);
template void fslvtkIO::addFieldData<float>(const float & M, const string & name, const string & type);

void fslvtkIO::addFieldData(vector< string > str, string name){
	fieldDataStr.push_back(str);
	fieldDataStrName.push_back(name);
	
}
	

template<class T>
void  fslvtkIO::writePointData(ofstream & fshape, const string & str_typename )
{
	//handle point data 
	if ( (Scalars.Nrows()>0) || (Vectors.Nrows()>0)  ) 
	{
		fshape<<"POINT_DATA "<<Points.Nrows()<<endl;
		if (Scalars.Nrows()>0)
		{
			fshape<<"SCALARS "<<scalarsName<<" "<<str_typename<<endl;
			fshape<<"LOOKUP_TABLE default"<<endl;
			for (int i=0;i<Scalars.Nrows();i++){
				for (int j=0;j<Scalars.Ncols();j++){
					if (BINARY)
					{
						T val=Scalars.element(i,j);
						Swap_Nbytes(1,sizeof(val),&val);	
						fshape.write(reinterpret_cast<char *>(&val),sizeof(val));
					}else
					{
						if (j==(Scalars.Ncols()-1)){
							fshape<<Scalars.element(i,j)<<endl;
						}else{
							fshape<<Scalars.element(i,j)<<" ";
						}
					}
				}
#ifdef PPC64
				if ((m_n++ % 20) == 0) fshape.flush();
#endif
			}
			
		}
		if (Vectors.Nrows()>0)
		{
			fshape<<"VECTORS "<<vectorsName<<" "<<str_typename<<endl;
			for (int i=0;i<Vectors.Nrows();i++){
				for (int j=0;j<Vectors.Ncols();j++){
					if (BINARY)
					{
						T val=Vectors.element(i,j);
						
						Swap_Nbytes(1,sizeof(val),&val);	
						fshape.write(reinterpret_cast<char *>(&val),sizeof(val));
					}else
					{
						if (j==(Vectors.Ncols()-1)){
							fshape<<Vectors.element(i,j)<<endl;
						}else{
							fshape<<Vectors.element(i,j)<<" ";
						}
					}
				}
#ifdef PPC64
				if ((m_n++ % 20) == 0) fshape.flush();
#endif
			}
			
		}
		
		
	}
	
	
}


template<class T>
void  fslvtkIO::writePoints(ofstream & fshape, const string & str_typename ){
	//only handle polydata currently
	//writing points to vtk file
	if (Points.Nrows()>0)
	{
		fshape<<"POINTS "<<Points.Nrows()<<" "<<str_typename<<endl;
		if (Points.Ncols()!=3)
			throw fslvtkIOException("Points does not have 3 columns");
		
		for (int i=0;i<Points.Nrows();i++)
		{
			//cols should always be three
			if (BINARY)
			{
				T tX=Points.element(i,0);
				T tY=Points.element(i,1);
				T tZ=Points.element(i,2);

				Swap_Nbytes(1,sizeof(tX),&tX);	
				Swap_Nbytes(1,sizeof(tY),&tY);	
				Swap_Nbytes(1,sizeof(tZ),&tZ);	
				
				fshape.write(reinterpret_cast<char *>(&tX),sizeof(tX));
				fshape.write(reinterpret_cast<char *>(&tY),sizeof(tY));
				fshape.write(reinterpret_cast<char *>(&tZ),sizeof(tZ));
			}else
				fshape<<Points.element(i,0)<<" "<<Points.element(i,1)<<" "<<Points.element(i,2)<<endl;
		}
#ifdef PPC64
		if ((m_n++ % 20) == 0) fshape.flush();
#endif
	}
	
	
}

template void fslvtkIO::writePoints<float>(ofstream & fshape, const string & str_typename );
template void fslvtkIO::writePoints<double>(ofstream & fshape, const string & str_typename );


bool  fslvtkIO::readPoints(ifstream & fvtk)
{
	int Npts;
	string stemp;
	fvtk>>stemp>>Npts;

	if ((strcmp(stemp.c_str(),"POINTS")) || (Npts<=0)) throw fslvtkIOException("POINTS not found");
	
	fvtk>>stemp;
	
	Points.ReSize(Npts,3);
	
	if (BINARY) 
		getline(fvtk,stemp); //gets rid of newline		
	
	for (int i=0 ; i < Npts ; i++)
	{
		float x,y,z;
		if (!BINARY)
			fvtk>>x>>y>>z;
		else
		{
			fvtk.read(reinterpret_cast<char*>(&x),sizeof(float));
			fvtk.read(reinterpret_cast<char*>(&y),sizeof(float));
			fvtk.read(reinterpret_cast<char*>(&z),sizeof(float));
			
			if (SWAP_BYTES)
			{
				Swap_Nbytes(1,sizeof(x),&x);
				Swap_Nbytes(1,sizeof(y),&y);	
				Swap_Nbytes(1,sizeof(z),&z);
			}
			
			
		}
		Points.element(i,0)=x;
		Points.element(i,1)=y;
		Points.element(i,2)=z;
				
	}
	
	return true;
}

void fslvtkIO::setPolygons(const vector< vector<unsigned int> > & vm)
{ 
	Matrix m(vm.size(),vm.at(0).size());
	for (unsigned int i = 0; i<vm.size();i++)
		for (unsigned int j=0;j<vm.at(0).size();j++)
			m.element(i,j)=vm.at(i).at(j);
	Polygons=m; 
}


bool  fslvtkIO::readPolygons(ifstream & fvtk)
{
	int NPolys;
	string stemp;
	fvtk>>stemp>>NPolys;
	if (strcmp(stemp.c_str(),"POLYGONS")) throw fslvtkIOException("POLYGONS not found");
	fvtk>>stemp;
	
	Polygons.ReSize(NPolys,3);
	
	if (BINARY) getline(fvtk,stemp); 
	for (int i=0 ; i < NPolys ; i++)
	{
		unsigned int x,y,z;
		
		if (!BINARY)
			fvtk>>x>>x>>y>>z;//just ignore number of connections assumed to be 3 for polydata
		else
		{
			fvtk.read(reinterpret_cast<char*>(&x),sizeof(unsigned int)); //just ignore number of connections assumed to be 3 for polydata
			fvtk.read(reinterpret_cast<char*>(&x),sizeof(unsigned int));
			fvtk.read(reinterpret_cast<char*>(&y),sizeof(unsigned int));
			fvtk.read(reinterpret_cast<char*>(&z),sizeof(unsigned int));
			if (SWAP_BYTES)
			{
				Swap_Nbytes(1,sizeof(x),&x);
				Swap_Nbytes(1,sizeof(y),&y);	
				Swap_Nbytes(1,sizeof(z),&z);
			}
		}
		Polygons.element(i,0)=x;
		Polygons.element(i,1)=y;
		Polygons.element(i,2)=z;
	}
	return true;
}


void  fslvtkIO::writePolygons(ofstream & fshape)
{
	
	if (Polygons.Nrows()>0){
		fshape<<"POLYGONS "<<Polygons.Nrows()<<"  "<<Polygons.Nrows()*(Polygons.Ncols()+1)<<endl;
		for (int i=0;i<Polygons.Nrows();i++)
		{
			for (int j=0;j<Polygons.Ncols();j++)
			{
				if (BINARY)
				{
					if (j==0)
					{
						unsigned int val=static_cast<unsigned int>(Polygons.Ncols());
						Swap_Nbytes(1,sizeof(val),&val);	
						fshape.write(reinterpret_cast<char *>(&val),sizeof(val));
					}
					unsigned int val=static_cast<unsigned int>(Polygons.element(i,j));
					Swap_Nbytes(1,sizeof(val),&val);	
					fshape.write(reinterpret_cast<char *>(&val),sizeof(val));
				}else
				{
					if (j==0)
						fshape<<Polygons.Ncols()<<" ";
					
					if (j==(Polygons.Ncols()-1))
						fshape<<Polygons.element(i,j)<<endl;
					else
						fshape<<Polygons.element(i,j)<<" ";
				}
			}
#ifdef PPC64
			if ((m_n++ % 20) == 0) fshape.flush();
#endif
			
		}
	}
	
}
void  fslvtkIO::writeCells(ofstream & fshape)
{
	//cout<<"still need to implement the cell writer";
	//count connection
	int N=0;
	for (unsigned int i=0;i<Cells.size();i++)
		N+=Cells.at(i).size();
	
	fshape<<"Cells "<<Cells.size()<<" "<<N<<endl;
	for (unsigned int i=0;i<Cells.size();i++)
	{
		for (unsigned int j=0; j<Cells.at(i).size();j++)
			fshape<<Cells.at(i).at(j)<<" ";
		fshape<<endl;
	}
}

void  fslvtkIO::writeUnstructuredGridCellTypes(ofstream & fshape)
{
	fshape<<"CELL_TYPES "<<Cell_Types.size()<<endl;
	for (unsigned int i=0;i<Cell_Types.size();i++)
		fshape<<Cell_Types.at(i)<<endl;
}
void fslvtkIO::save(string s)
{
	//add required field data 
	if (pd_list.size()>0){
		addFieldData(pd_list, "PointFieldNames");
		addFieldData(pd_type, "PointFieldAttTypes");
	}
	if (cd_list.size()>0){
		addFieldData(cd_list, "CellFieldNames");
		addFieldData(cd_type, "CellFieldAttTypes");
	}
	
	cout<<"open file "<<s<<" to save."<<endl;
	ofstream fshape;  
	fshape.open(s.c_str());
	cout<<"succesfully opened file "<<s<<" to save."<<endl;

	//calculate total number of points
	
	
	fshape<<"# vtk DataFile Version 3.0"<<endl;
	if (BINARY)
	{
		int test=42;
		Swap_Nbytes(1,sizeof(test),&test);
		fshape.write(reinterpret_cast<char *>(&test),sizeof(test));
		fshape<<"this file was written using fslvtkio"<<endl<<"BINARY"<<endl<<"DATASET ";
	}else
		fshape<<"this file was written using fslvtkio"<<endl<<"ASCII"<<endl<<"DATASET ";
	
	switch (dt) {
		case POLYDATA:
			fshape<<"POLYDATA"<<endl;
			writePoints<float>(fshape,"float");
			writePolygons(fshape);
			break;
		case UNSTRUCTURED_GRID:
			fshape<<"UNSTRUCTURED_GRID"<<endl;
			writePoints<float>(fshape,"float");
			writeCells(fshape);
			writeUnstructuredGridCellTypes(fshape);
			
			break;
		default:
			cerr<<"Invalid Data Type"<<endl;
	}
	
	writePointData<float>(fshape,"float");
	
	//Ignored cell data for now 
	
	//write field data
	if ( ( fieldDataStr.size()>0) || (fieldDataNum.size()>0) )
	{
		fshape<<"FIELD FieldData"<<" "<<fieldDataNum.size()+fieldDataStr.size()<<endl;
		//write out numerixc field data
		if (fieldDataNum.size()>0) {
			for (unsigned int i=0; i<fieldDataNum.size();i++)
				if ((MAX_SET) && (static_cast<unsigned int>(fieldDataNum.at(i).Ncols())>MAX) ) //Maximum limits the number of allowable columns in the field data (used to limit saved modes of variation)
					writeNumericField<float>(fshape, fieldDataNumName.at(i), "float", fieldDataNum.at(i).SubMatrix(1,fieldDataNum.at(i).Nrows(),1,MAX));
				else
					writeNumericField<float>(fshape, fieldDataNumName.at(i), "float", fieldDataNum.at(i));
		}
	}
	
				vector<string>::iterator name_i=fieldDataStrName.begin();
				for (vector< vector<string> >::iterator i=fieldDataStr.begin(); i!=fieldDataStr.end();i++, name_i++)
					writeStringField(fshape, *name_i, *i);
				
				fshape.close();
				
}

void fslvtkIO::writeStringField(ofstream & fvtk, const string & name, const vector<string> & v_string)
{	
	fvtk<<name<<" "<<1<<" "<<v_string.size()<<" string"<<endl;
	for (vector<string>::const_iterator i=v_string.begin(); i!=v_string.end();i++)
	{
		if (i==v_string.begin())
			fvtk<<i->c_str();
		else
			fvtk<<" "<<i->c_str();
#ifdef PPC64
		if ((m_n++ % 50) == 0) fvtk.flush();
#endif
	}
}




template <class T>
void fslvtkIO::writeNumericField(ofstream & fvtk, const string & name, const string & type, const Matrix & Data)
{
	unsigned int nrows=Data.Nrows();
	unsigned int ncols=Data.Ncols();
	fvtk<<name<<" "<<nrows<<" "<<ncols<<" "<<type<<endl;
	
	for (unsigned int i=0; i<nrows ;i++)
		for (unsigned int j=0;j<ncols;j++)
		{
			if (!BINARY)
			{
				if (j==(ncols-1))
					fvtk<<Data.element(i,j)<<endl;
				else
					fvtk<<Data.element(i,j)<<" ";
			}else
			{
				T val=static_cast<T>(Data.element(i,j));
				Swap_Nbytes(1,sizeof(val),&val);	
				fvtk.write(reinterpret_cast<char*>(&val),sizeof(val));
			}
			
#ifdef PPC64
			if ((m_n++ % 20) == 0) fvtk.flush();
#endif
		}
}


template <class T>
ReturnMatrix fslvtkIO::readField(ifstream & fvtk, const int & nrows,const int & mcols)
{
	Matrix fieldM(nrows,mcols);
	T val;
	for (int i=0; i<nrows ;i++){
		for (int j=0;j<mcols;j++){
			
			if (!BINARY)
				fvtk>>val;
			else
			{
				fvtk.read(reinterpret_cast<char*>(&val),sizeof(T));
				if (SWAP_BYTES)
					Swap_Nbytes(1,sizeof(val),&val);	
			}
			fieldM.element(i,j)=val;
			//	cout<<"val "<<val<<endl;
		}
	} 
	
	fieldM.Release();
	return fieldM;
}

void fslvtkIO::readFieldData(ifstream & fvtk){
	//blanks out previous field data
	fieldDataNumName.clear();
	fieldDataNum.clear();
	fieldDataNumType.clear();
	
	string stemp;
	int N;//number of fields
		fvtk>>stemp>>N;
		
		if (SWITCH_ROWS_COLS){ N-=1; }//only used to fix old model!!!!!!!!
		
		//read in field
		for ( int fieldnum=0; fieldnum< N ; fieldnum++){
			string fieldname;
			fvtk>>fieldname;
			long unsigned int nrows,mcols;
			
			if (SWITCH_ROWS_COLS)
			{//only used to fix old model!!!!!!!!
				fvtk>>mcols>>nrows;
				if (mcols==1){ mcols=nrows; nrows=1; } 
			}else
				fvtk>>nrows>>mcols;
			fvtk>>stemp;
			
			if (BINARY) { string stemp2 ; getline(fvtk,stemp2); };
			if (! ((strcmp(stemp.c_str(),"float")) && (strcmp(stemp.c_str(),"unsigned int")) && (strcmp(stemp.c_str(),"double")) && (strcmp(stemp.c_str(),"int"))) )  {
				//deal with number field data
				fieldDataNumType.push_back(stemp);
				fieldDataNumName.push_back(fieldname);
				Matrix fieldM;
				if (!(strcmp(stemp.c_str(),"float")))  
					fieldM=readField<float>(fvtk,nrows,mcols); 
				else if (!(strcmp(stemp.c_str(),"double")))  
					fieldM=readField<double>(fvtk,nrows,mcols); 
				else if (!(strcmp(stemp.c_str(),"unsigned int")))  
					fieldM=readField<unsigned int>(fvtk,nrows,mcols); 
				else if (!(strcmp(stemp.c_str(),"int"))) 
					fieldM=readField<int>(fvtk,nrows,mcols); 
				fieldDataNum.push_back(fieldM);
				
			}else if (!(strcmp(stemp.c_str(),"string")))
				fieldDataStrName.push_back(fieldname);	//deal with strings separately 
			else if(fvtk.eof())
				return;
			else
				throw fslvtkIOException(("Data type for field data not supported..."+stemp).c_str()) ;
		}
		
}
void fslvtkIO::readPointData( ifstream & fvtk, string & nextData){
	//only handles SCALARS and VECTORS currently
	//assumes POINT_DATA is found
        // WARNING!!  NEED TO DEAL WITH nextData BEING SET, AS IT INDICATES THAT A STRING HAS BEEN 
        //   REMOVED FROM THE STREAM BUT NOT YET PROCESSED!!
	string stemp,stype;
	int N;//number of Points
		fvtk>>N;
		
		if (N<=0)
			throw fslvtkIOException("no points in structure") ;
		
		//this will not read past any unsupported data
		bool still_PD=true;
		while (still_PD){
			fvtk>>stemp;		
			
			if (!strcmp(stemp.c_str(),"SCALARS")){
				fvtk>>stemp>>stype;
				string lut;
				fvtk>>lut>>lut;
				
				//there is now an options fourth number of component argument not implemented here
				int mcols=1;
				
				if (BINARY) 
				{ 
					string stemp2 ; 
					getline(fvtk,stemp2); 
				}
				
				if (! ((strcmp(stype.c_str(),"float") ) && (strcmp(stype.c_str(),"unsigned int")) && \
					   (strcmp(stype.c_str(),"double")) && (strcmp(stype.c_str(),"int")) ))  
				{
					//deal with number field data
					Matrix fieldM;
					if (!(strcmp(stype.c_str(),"float")))				fieldM=readField<float>(fvtk,N,mcols); 
					else if (!(strcmp(stype.c_str(),"double")))			fieldM=readField<double>(fvtk,N,mcols); 
					else if (!(strcmp(stype.c_str(),"unsigned int")))	fieldM=readField<unsigned int>(fvtk,N,mcols); 
					else if (!(strcmp(stype.c_str(),"int")))			fieldM=readField<int>(fvtk,N,mcols); 
					Scalars=fieldM;
				}else
					throw fslvtkIOException("Data type for points not supported.");
				
			}else if (!strcmp(stemp.c_str(),"VECTORS")){
				
				fvtk>>stemp>>stype;
				
				if (! ((strcmp(stype.c_str(),"float") ) && (strcmp(stype.c_str(),"unsigned int")) && \
					   (strcmp(stype.c_str(),"double")) && (strcmp(stype.c_str(),"int")) ))  
				{
					//deal with number field data
					Matrix fieldM;
					if (!(strcmp(stype.c_str(),"float")))				fieldM=readField<float>(fvtk,N,3); 
					else if (!(strcmp(stype.c_str(),"double")))			fieldM=readField<double>(fvtk,N,3); 
					else if (!(strcmp(stype.c_str(),"unsigned int")))	fieldM=readField<unsigned int>(fvtk,N,3); 
					else if (!(strcmp(stype.c_str(),"int")))			fieldM=readField<int>(fvtk,N,3); 
					Vectors=fieldM;
					
				}else
					throw fslvtkIOException("Data type for vectors not supported.");
				
			}else
			{
				nextData=stemp;
				still_PD=false;
			}
			
			
		}	
		
}

void fslvtkIO::readUnstructuredGrid(string fname)
{
	Cells.clear();
	Cell_Types.clear();
	
	ifstream fvtk;
	fvtk.open(fname.c_str());
	
	string stemp;
	getline(fvtk,stemp);
	getline(fvtk,stemp);//used to get rid of second ascii line to avoid "modes"
		
		fvtk>>stemp;
		fvtk>>stemp>>stemp;
		
		readPoints(fvtk);
		
		fvtk>>stemp;
		
		int Ncells, Csize;
		fvtk>>Ncells>>Csize;
		
		for ( int i =0 ; i<Ncells ; i++)
		{
			unsigned int temp1, temp2;
			fvtk>>temp1;
			vector<unsigned int> vcell;
			vcell.push_back(temp1);
			for (unsigned int j=0; j<temp1 ;j++)
			{
				fvtk>>temp2;
				vcell.push_back(temp2);
			}
			Cells.push_back(vcell);
		}
		fvtk>>stemp>>stemp;
		for ( int i =0 ; i<Ncells ; i++)
		{
			unsigned int temp1;
			fvtk>>temp1;
			Cell_Types.push_back(temp1);
		}
		while (fvtk>>stemp)
		{
			if (!strcmp(stemp.c_str(),"POINT_DATA"))
			{
				readPointData(fvtk,stemp);
			}else if(!strcmp(stemp.c_str(),"FIELD"))
			{	
				readFieldData(fvtk);
			}			
		}
}

void fslvtkIO::readPolyData(string name)
{
	ifstream fvtk;
	fvtk.open(name.c_str());
	
	if (fvtk.is_open())
	{		
		string stemp;
		getline(fvtk,stemp); //read first line 
		
		if (strcmp(stemp.substr(0,14).c_str(),"# vtk DataFile")) 
			throw fslvtkIOException("Not a vtk file (error in line 1).");
		
		getline(fvtk,stemp); 
		getline(fvtk,stemp);
		
		if (((strcmp(stemp.c_str(),"ASCII")) && (strcmp(stemp.c_str(),"BINARY")))) 
			throw fslvtkIOException("ASCII or Binary not specified (line 3)"); 
		else if (!strcmp(stemp.c_str(),"BINARY"))//if binary file test endianess, test number embedded on second line
		{
			BINARY=true;
			//reopen file and read binary flag from first line
			unsigned int testval;
			{//this bit just gets the test val
				ifstream* fvtktemp = new ifstream;
				fvtktemp->open(name.c_str());
				getline(*fvtktemp,stemp);//throw away first line
					fvtktemp->read(reinterpret_cast<char *>(&testval),sizeof(testval));
					fvtktemp->close();
					delete fvtktemp;
			}
			if (testval!=BINFLAG) 
			{
				SWAP_BYTES = true;
				Swap_Nbytes(1,sizeof(testval),&testval);
			}
			if (testval!=BINFLAG)
				throw fslvtkIOException("Unrecognised binary matrix file format");
		}
		//once this point is reached the endianess if okay
		//SWAP_BYTE and BINARY state variables have been set such that readPoints, Polygons, etc... will work properly
		
		getline(fvtk,stemp);
		if (strcmp(stemp.c_str(),"DATASET POLYDATA")) 
			throw fslvtkIOException("Is not specified as Polydata (line 4");
		
		readPoints(fvtk);
		readPolygons(fvtk);
		
		bool already_read_stemp=false;
		while (already_read_stemp || (fvtk>>stemp))
		{
		  already_read_stemp=false;
		  if (!strcmp(stemp.c_str(),"POINT_DATA")) {
		    readPointData(fvtk,stemp);
		    if (stemp.length()>0) {already_read_stemp=true;}
		  } else { 
		    if(!strcmp(stemp.c_str(),"FIELD"))
		      readFieldData(fvtk);
		  }
		}
	}else{
		throw fslvtkIOException("Cannot open file.");
	}
}

	void fslvtkIO::displayNumericFieldDataNames()
	{
		for (vector<string>::iterator i=fieldDataNumName.begin();i!=fieldDataNumName.end();i++)
			cout<<*i<<endl;
	}
	
	void fslvtkIO::displayNumericField(const string & name)
	{
		cout<<getField(name)<<endl;
		
	}
	
}


