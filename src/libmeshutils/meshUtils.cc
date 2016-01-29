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
//#include <math.h>
//#define USE_VTK
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <algorithm>
#include "math.h"
#include "meshUtils.h"
#include "fslvtkio/fslvtkio.h"

#include "utils/options.h"
#include "newimage/newimageall.h"
#include "meshclass/meshclass.h"
#include "first_lib/first_mesh.h"
#include "first_lib/first_newmat_vec.h"

//------------------------VTK STUFF------------------------//
#ifdef USE_VTK

#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"

#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"

#include "vtkPoints.h"
#include "vtkProperty.h"
#include "vtkImagePlaneWidget.h"
#include "vtkImageData.h"

#include <fslvtkconv/fslvtkconv.h>

using namespace fslvtkconvname;

#endif

using namespace std;
using namespace NEWIMAGE;

using namespace mesh;
using namespace fslvtkio;
using namespace FIRST_LIB;
// Local functions

struct float3
{
	float x,y,z;
};

struct float2
{
	float x,y;
};

//static float2 make_float2(float x, float y)
//{
//	float2 t; t.x = x; t.y = y; return t;
//
//}

static float3 make_float3(float x, float y, float z)
{
	float3 t; t.x = x; t.y = y; t.z = z; return t;

}
inline float dot_prod(const float3 & v1 , const float3 & v2)
{

	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

inline float3 sub_vec(const float3 & v1 , const float3 & v2)
{
	
	return make_float3(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z); 
}

inline float3 add_vecs(const float3 & v1 , const float3 & v2, const float3 & v3)
{
	
	return make_float3(v1.x+v2.x+v3.x, v1.y+v2.y+v3.y, v1.z+v2.z+v3.z); 
}
inline float3 mul_vec(const float & sc , const float3 & v1)
{
	
	return make_float3(sc*v1.x, sc*v1.y, sc*v1.z); 
}


namespace meshutils {
	
	meshUtils::meshUtils(){
		
	}
	
	meshUtils::meshUtils(const string & filename,const fslvtkIO::DataType i)
{
		
		try 
	{
		switch (i)
		{
			case 0:
				setDataType(i);
				readPolyData(filename);
				break;
			case 1:
				setDataType(i);
				cout<<"read unstructured grid meshutils"<<endl;
				readUnstructuredGrid(filename);
								cout<<"end read unstructured grid"<<endl;

				break;
			default:
				throw -8;
		}
	}catch (int e)
	{
		cout<<"invalid data type"<<endl;
		return;
	}
		
}
	void meshUtils::loadMesh(const string & meshname){
		//	mesh1=new fskvtkIO;
			setDataType(static_cast<fslvtkIO::DataType>(0));
			readPolyData(meshname);
	}
	
	void meshUtils::getBounds(int *bounds,const float & xdim, const float & ydim,const float & zdim)
	{
		float xmin=1000,xmax=-1000,ymin=1000,ymax=-1000,zmin=1000,zmax=-1000;
	
		for (int i = 0; i<Points.Nrows(); i++ ){
			float tempx=Points.element(i,0);
			float tempy=Points.element(i,1);
			float tempz=Points.element(i,2);
			if (tempx<xmin){
				xmin=tempx;
			}
			if (tempx>xmax){
				xmax=tempx;
			}
			if (tempy<ymin){
				ymin=tempy;
			}
			if (tempy>ymax){
				ymax=tempy;
			}
			if (tempz<zmin){
				zmin=tempz;
			}
			if (tempz>zmax){
				zmax=tempz;
			}
		}
		*bounds=static_cast<int>(floor(xmin/xdim)-1);
		*(bounds+1)=static_cast<int>(ceil(xmax/xdim)+1);
		*(bounds+2)=static_cast<int>(floor(ymin/ydim)-1);
		*(bounds+3)=static_cast<int>(ceil(ymax/ydim)+1);
		*(bounds+4)=static_cast<int>(floor(zmin/zdim)-1);
		*(bounds+5)=static_cast<int>(ceil(zmax/zdim)+1);
		
	}
	
	void meshUtils::generateRandomMeshUsingScalar(const Mesh & min, const string & outname, const vector<bool> & scal, const int & N)
{

	srand ( time(NULL) );
	for (int i=0; i<N;i++)
	{
	
	stringstream sout;
		sout<<i;
		string sind;
		sout>>sind;
		
		Mesh m=min;
		
		int count=0;
	    for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ , count++ )
		{
			cout<<count<<" "<<scal.at(count)<<endl;
			if (scal.at(count)){
				//	cout<<static_cast<float>(rand()-RAND_MAX/2.0)<<" "<<RAND_MAX<<" "<<static_cast<float>(rand()-RAND_MAX/2.0)/RAND_MAX<<endl;
				float rn=static_cast<float>(rand()-RAND_MAX/2.0)/RAND_MAX*1.0+0.5;
				//	cout<<"rn "<<rn<<endl;
				
				(*i)->_update_coord = (*i)->get_coord() + rn * ((*i)->local_normal());
			}else
			{	
				
				//						cout<<static_cast<float>(rand()-RAND_MAX/2.0)<<" "<<RAND_MAX<<" "<<static_cast<float>(rand()-RAND_MAX/2.0)/RAND_MAX<<endl;
				
				float rn2=static_cast<float>(rand()-RAND_MAX/2.0)/RAND_MAX*1.0+0;
				//	cout<<"rn2 "<<rn2<<endl;
				
				(*i)->_update_coord = (*i)->get_coord() + rn2 * ((*i)->local_normal());
			}
		}
		m.update();
		m.save(outname+sind+".vtk",3);
	}
}
	float meshUtils::drawTriangleScalars(volume<float>& image, volume<int>& count, const unsigned int & tri_index)
	{
		
		float3 p0,p1,p2;
		float sc0,sc1,sc2;

		for (int pivot=0;pivot<3;pivot++)
		{
		switch (pivot)
			{
				case 0:
		p0 =  make_float3( Points.element(Polygons.element(tri_index,0),0), Points.element(Polygons.element(tri_index,0),1), Points.element(Polygons.element(tri_index,0),2) );
		 p1 =  make_float3( Points.element(Polygons.element(tri_index,1),0), Points.element(Polygons.element(tri_index,1),1), Points.element(Polygons.element(tri_index,1),2) );
		 p2 =  make_float3( Points.element(Polygons.element(tri_index,2),0), Points.element(Polygons.element(tri_index,2),1), Points.element(Polygons.element(tri_index,2),2) );
		
					sc0=Scalars.element(Polygons.element(tri_index,0),0);	
					sc1=Scalars.element(Polygons.element(tri_index,1),0);	
					sc2=Scalars.element(Polygons.element(tri_index,2),0);	

					break;
				case 1:
					p1 =  make_float3( Points.element(Polygons.element(tri_index,0),0), Points.element(Polygons.element(tri_index,0),1), Points.element(Polygons.element(tri_index,0),2) );
					p0 =  make_float3( Points.element(Polygons.element(tri_index,1),0), Points.element(Polygons.element(tri_index,1),1), Points.element(Polygons.element(tri_index,1),2) );
					p2 =  make_float3( Points.element(Polygons.element(tri_index,2),0), Points.element(Polygons.element(tri_index,2),1), Points.element(Polygons.element(tri_index,2),2) );
					
					sc1=Scalars.element(Polygons.element(tri_index,0),0);	
					sc0=Scalars.element(Polygons.element(tri_index,1),0);	
					sc2=Scalars.element(Polygons.element(tri_index,2),0);	
					
					
					break;
				case 2:
					p2 =  make_float3( Points.element(Polygons.element(tri_index,0),0), Points.element(Polygons.element(tri_index,0),1), Points.element(Polygons.element(tri_index,0),2) );
					p1 =  make_float3( Points.element(Polygons.element(tri_index,1),0), Points.element(Polygons.element(tri_index,1),1), Points.element(Polygons.element(tri_index,1),2) );
					p0 =  make_float3( Points.element(Polygons.element(tri_index,2),0), Points.element(Polygons.element(tri_index,2),1), Points.element(Polygons.element(tri_index,2),2) );
					
					sc2=Scalars.element(Polygons.element(tri_index,0),0);	
					sc1=Scalars.element(Polygons.element(tri_index,1),0);	
					sc0=Scalars.element(Polygons.element(tri_index,2),0);	
					
					break;
			}
		
		
		float3 vu=sub_vec(p1,p0);
		float3 vv=sub_vec(p2,p0);
			float sc_u=sc1-sc0;
			float sc_v=sc2-sc0;
			
			

		float xdim=image.xdim();
		float ydim=image.ydim();
		float zdim=image.zdim();
		
		for (float u=0;u<=1;u+=0.01)
			for (float v=0;v<=1;v+=0.01)
			{
				if ( (u+v)<=1 )
				{
					float3 P=add_vecs(p0, mul_vec(u,vu), mul_vec(v,vv));
					//if ((sc0 + u*sc_u +v*sc_v)<0)
					//	cout<<"scalar "<<sc0<<" "<<sc1<<" "<<sc2<<" "<<sc0 + u*sc_u +v*sc_v<<endl;
					image.value(static_cast<int>(P.x/xdim+0.5), static_cast<int>(P.y/ydim+0.5), static_cast<int>(P.z/zdim+0.5) )+=  sc0 + u*sc_u +v*sc_v;
					count.value(static_cast<int>(P.x/xdim+0.5), static_cast<int>(P.y/ydim+0.5), static_cast<int>(P.z/zdim+0.5) )++;

				
				}	
			}
		
		
		}
		// Brian's dodgy code - return a number
		return -1;
	  }

	
	
	void meshUtils::generateRandom6DOFMatrices( const string & outname, const int & N)
{
	
	const float PI=3.14159;
	//allow to go plus minus pi/4 in each direction
	
	srand ( time(NULL) );
	for (int i=0; i<N;i++)
	{
	
	stringstream sout;
		sout<<i;
		string sind;
		sout>>sind;
	
		
		float alpha=static_cast<float>(rand()-RAND_MAX/2.0)*2*PI/(4.0*RAND_MAX) ;
		float beta=static_cast<float>(rand()-RAND_MAX/2.0)*2*PI/(4.0*RAND_MAX) ;
		float gam=static_cast<float>(rand()-RAND_MAX/2.0)*2*PI/(4.0*RAND_MAX) ;

		Matrix Rx(4,4);
		Matrix Ry(4,4);
		Matrix Rz(4,4);
		
		Rx=0;
		Rx.element(0,0)=1;
		Rx.element(1,1)=cos(alpha);
		Rx.element(1,2)=sin(alpha);
		Rx.element(2,1)=-Rx.element(1,2);
		Rx.element(2,2)=Rx.element(1,1);
		Rx.element(3,3)=1;

		Ry=0;
		Ry.element(1,1)=1;
		Ry.element(0,0)=cos(beta);
		Ry.element(0,2)=-sin(beta);
		Ry.element(2,0)=-Ry.element(0,2);
		Ry.element(2,2)=Ry.element(1,1);
		Ry.element(3,3)=1;

		Rz=0;
		Rz.element(2,2)=1;
		Rz.element(0,0)=cos(gam);
		Rz.element(0,1)=sin(gam);
		Rz.element(1,0)=-Rz.element(0,1);
		Rz.element(1,1)=Rz.element(0,0);
		Rz.element(3,3)=1;


		Matrix final=Rx*Ry*Rz;
		
		float tx=static_cast<float>(rand()-RAND_MAX/2.0)*2*50/(RAND_MAX) ;
		float ty=static_cast<float>(rand()-RAND_MAX/2.0)*2*50/(RAND_MAX) ;
		float tz=static_cast<float>(rand()-RAND_MAX/2.0)*2*50/(RAND_MAX) ;

		final.element(0,3)=tx;
		final.element(1,3)=ty;
		final.element(2,3)=tz;
		ofstream fout;
		fout.open((outname+"_"+sind+".mat").c_str());
		fout<<final;
		fout.close();
	}
}

/*
		void meshUtils::getConditionalMeanAndVariance(shapemodel::shapeModel * min, volume4D<float> & iCondMean, volume4D<float> & iCondVar , const volume<float> & im, const int & mode, const float & bmin, const float & bmax, const float & res, const float & mean_offset)
{

		vector<float> vars;
		for (int i=0;i<=mode;i++)
			vars.push_back(0);
		
		volume4D<float> meanI4D;
		volume4D<float> varI4D;
		
		
		cout<<"calculate  base covariance"<<endl;
		vector< vector<float> > prec=min->getShape(0)->getICondPrec();
		vector<float> eigs=min->getShape(0)->getICondEigs();
		unsigned int N=min->getDeformedIprof(vars,mode,vars.size()).size();
			cout<<"convert vector to MAtrix"<<endl;
			Matrix Mprec=vectorOfVectorsToMatrix(prec);
cout<<"converted vector to MAtrix "<<Mprec.Nrows()<<" "<<Mprec.Ncols()<<endl;

time_t start,end;
		cout<<"time begin "<<time(&start)<<endl;
		Matrix Mcov=Mprec*Mprec.t();
		cout<<"time end"<<time(&end)<<endl;
		cout<<"time diff "<<difftime(end,start)<<endl;
		cout<<"done ultiply"<<endl;

				cout<<"time begin "<<time(&start)<<endl;
		
		for ( float b=bmin; b<bmax;b++)
			  {
				  cout<<"set mode "<<mode<<" to "<<b<<" "<<vars.size()<<endl; 
					vars.at(mode)=b;
				  cout<<"vars adjusted "<<endl;
				  vector<float> iprof=min->getDeformedIprof(vars,mode,vars.size());
				 				  
				  vector< vector<float> > vcov;
				  vector<float> var;
				
				  if ( !(eigs.size()==prec.size()) )
					  throw fslvtkIOException("number of eigenvalues does not match number of number modes");
				  
				
				volume<float> imean=im;
				volume<float> ivar=im;
				imean=0;
				ivar=0;
				float xdim=im.xdim(),ydim=im.ydim(),zdim=im.zdim() ;
				
				Mesh m= min->getDeformedMeshTrans(vars,mode,vars.size());
		
		vector<float>::iterator iprof_iter=iprof.begin();
		int ipp=min->getIPP(0);
				cout<<"ipp "<<ipp<<endl;

		cout<<"iprof "<<iprof.at(0)<<" "<<iprof.size()<<endl;
		int count=0;
		for (vector<Mpoint*>::iterator k = m._points.begin(); k!=m._points.end(); k++ )
		{
			Pt vertNew;
			Pt vert = (*k)->get_coord();
			Vec normal;
			normal = (*k)->local_normal();
			for (int j=0;j<ipp;j++,iprof_iter++,count++){
				int sc=j-(ipp-1)/2;
				vertNew=vert + normal*0.5*sc;
				//cout<<vertNew.X/xdim<<" "<<vertNew.Y/ydim<<" "<<vertNew.Z/zdim<<" "<<(*iprof_iter)+mean_offset<<endl;
				imean.value(vertNew.X/xdim,vertNew.Y/ydim,vertNew.Z/zdim)=(*iprof_iter)+mean_offset;
				ivar.value(vertNew.X/xdim,vertNew.Y/ydim,vertNew.Z/zdim)=Mprec.element(count,count);

			}
			
		}
		meanI4D.addvolume(imean);
		varI4D.addvolume(ivar);
		}

		iCondMean=meanI4D;
		iCondVar=varI4D;
}

*/
/*	void meshUtils::addModesToModelUsingMask(shapeModel * min, const vector<bool> & scal)
{

//	float mag=0;
//	for (vector<float>::const_iterator i=scal.begin();i!=scal.end();i++)
//	{
//		if (*i) mag++;
//	}
//	mag=sqrt(mag);
		
		Mesh m =min->getShapeMesh(0);
		unsigned int count=0;
	    for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ , count++ )
		{
			vector< Vec > shapeMode;
			
			if (scal.at(count)){
				vector<float> eigNew;
				eigNew.push_back(1);
				for (int j=0;j<min->getNumberOfModes();j++)
				{
					eigNew.push_back(min->getEigenValue(j));
				}	
				min->setEigenValues(eigNew);
				
				for (unsigned int j=0;j<scal.size();j++)
				{
					if (j==count) 
						shapeMode.push_back((*i)->local_normal());
					else
					{
						Vec temp(0,0,0);
						shapeMode.push_back(temp);
					}
				}
				min->getShape(0)->insertModeVector(shapeMode,0);
			}
		}
}

*/
//static function that takes in mesh class mesh
	void meshUtils::getBounds(const Mesh & m, int *bounds, const float & xdim, const float & ydim,const float & zdim){
		
		float xmin=1000,xmax=-1000,ymin=1000,ymax=-1000,zmin=1000,zmax=-1000;
		for (vector<Mpoint*>::const_iterator i = m._points.begin(); i!=m._points.end(); i++ ){
			float tempx=(*i)->get_coord().X;
			float tempy=(*i)->get_coord().Y;
			float tempz=(*i)->get_coord().Z;
			if (tempx<xmin){
				xmin=tempx;
			}
			if (tempx>xmax){
				xmax=tempx;
			}
			if (tempy<ymin){
				ymin=tempy;
			}
			if (tempy>ymax){
				ymax=tempy;
			}
			if (tempz<zmin){
				zmin=tempz;
			}
			if (tempz>zmax){
				zmax=tempz;
			}
		}
		*bounds=static_cast<int>(floor(xmin/xdim)-1);
		*(bounds+1)=static_cast<int>(ceil(xmax/xdim)+1);
		*(bounds+2)=static_cast<int>(floor(ymin/ydim)-1);
		*(bounds+3)=static_cast<int>(ceil(ymax/ydim)+1);
		*(bounds+4)=static_cast<int>(floor(zmin/zdim)-1);
		*(bounds+5)=static_cast<int>(ceil(zmax/zdim)+1);
		
	}
	
	//static
void meshUtils::fileToVector(const string & fname, vector<string> meshlist)
{
		meshlist.clear();
		ifstream fin;
		fin.open(fname.c_str());
		string meshname;
		while (fin>>meshname)
		{
			meshlist.push_back(meshname);
		}
}
vector<string> meshUtils::fileToVector(const string & fname)
{
	vector<string> meshlist;
	ifstream fin;
	fin.open(fname.c_str());
	string meshname;
	while (fin>>meshname)
	{
		meshlist.push_back(meshname);
	}
	return meshlist;
}

ReturnMatrix meshUtils::readFlirtMat(const string & fname)
{
	ifstream fin;
	fin.open(fname.c_str());
	float ftemp;
	Matrix fmat(4,4);
	for (int i =0 ; i<4;i++)
		for (int j =0 ; j<4;j++)
		{
			fin>>ftemp;
			fmat.element(i,j)=ftemp;
		}
	fmat.Release();
	return fmat;
	
}
void meshUtils::writeFlirtMatrix(const Matrix & fmat, const string & fname)
{
	ofstream fout;
	fout.open(fname.c_str());
	fout.flags( ios::fixed );
	for (int i =0 ; i<4;i++)
	{
		for (int j =0 ; j<4;j++)
			fout<<fmat.element(i,j)<<" ";
		fout<<endl;
	}
	fout.close();
}
//static for convenience
void meshUtils::meshReg(Mesh & m, const Matrix & fmat)
{
	//refsize is actually target image
	for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ )
	{
		float x=(*i)->get_coord().X;
		float y=(*i)->get_coord().Y;
		float z=(*i)->get_coord().Z;

		Pt newPt(fmat.element(0,0)*x+fmat.element(0,1)*y+fmat.element(0,2)*z+fmat.element(0,3),\
			fmat.element(1,0)*x+fmat.element(1,1)*y+fmat.element(1,2)*z+fmat.element(1,3), \
			fmat.element(2,0)*x+fmat.element(2,1)*y+fmat.element(2,2)*z+fmat.element(2,3));
		(*i)->_update_coord = newPt;
	}
	m.update();
}
void meshUtils::meshReg(const Matrix & fmat){
	ColumnVector ones(Points.Nrows());
	ones=1;
	Points=Points | ones;
	Points=(fmat*(Points.t())).t();
	Points=Points.SubMatrix(1,Points.Nrows(),1,3); 
}


	void meshUtils::shift3DVertexMatrix(Matrix & mat, const  float & tx, const float & ty, const float & tz ){
		for (int i=0; i<mat.Nrows();i++){
			mat.element(i,0)+=tx;
			mat.element(i,1)+=ty;
			mat.element(i,2)+=tz;
		}
	}
	void meshUtils::shift3DVertexColumnVector(ColumnVector & mat, const  float & tx, const float & ty, const float & tz ){
		for (int i=0; i<mat.Nrows();i+=3){
			mat.element(i)+=tx;
			mat.element(i+1)+=ty;
			mat.element(i+2)+=tz;
		}
	}
void meshUtils::shift3DMesh(Mesh & m, const float & tx, const float & ty, const float & tz )
{
	for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ )
	{
		(*i)->_update_coord.X=(*i)->get_coord().X+tx;
		(*i)->_update_coord.Y=(*i)->get_coord().Y+ty;
		(*i)->_update_coord.Z=(*i)->get_coord().Z+tz;
	}
	m.update();
	
}

	void meshUtils::shiftPoints( const  float & tx, const float & ty, const float & tz ){
		for (int i=0; i<Points.Nrows();i++){
			Points.element(i,0)+=tx;
			Points.element(i,1)+=ty;
			Points.element(i,2)+=tz;
			
		}
	}
	void meshUtils::scalePoints( const  float & sx, const float & sy, const float & sz ){
		for (int i=0; i<Points.Nrows();i++)
		{
			Points.element(i,0)*=sx;
			Points.element(i,1)*=sy;
			Points.element(i,2)*=sz;
		}
	}
	
	void meshUtils::shiftPolygonMatrix(const int & shift ){

		for (int i=0; i<Polygons.Nrows();i++)
			for (int j=0; j<Polygons.Ncols();j++)
				Polygons.element(i,j)+=shift;
		
	}
	
	ReturnMatrix meshUtils::shiftPolygonMatrix(const Matrix & mat, const int & shift ){

		Matrix polyShift(mat.Nrows(),mat.Ncols());
		for (int i=0; i<mat.Nrows();i++)
			for (int j=0; j<mat.Ncols();j++)
				polyShift.element(i,j)=mat.element(i,j)+shift;
		
		polyShift.Release();
		return polyShift;
	}
	ReturnMatrix meshUtils::meshPointsToMatrix(const Mesh & m1)
{
//this mimics the other except the reference should already be loaded 
	Matrix pointsRef(m1._points.size(),3);
	int count=0;
	for (vector<Mpoint*>::const_iterator i=m1._points.begin(); i!=m1._points.end();i++,count++)
	{
		pointsRef.element(count,0)=(*i)->get_coord().X;		
		pointsRef.element(count,1)=(*i)->get_coord().Y;		
		pointsRef.element(count,2)=(*i)->get_coord().Z;		
	}
	pointsRef.Release();
	return pointsRef;
}
	
bool meshUtils::checkLine(const float & p1, const float & p2, const float & test){
	return (((p1>test)&&(p2<test)) ||  ((p1<test)&&(p2>test)));
}


bool meshUtils::checkTriangleNeighbour(const short & tri0, const short & tri1, const short & tri2 , const short & ind0, const short & ind1, short & ind0new , short & ind1new){
	
	if (((tri0==ind0) || (tri1==ind0) || (tri2==ind0)) && ((tri0==ind1) || (tri1==ind1) || (tri2==ind1))){
		if (tri0==ind0) { ind0new=0 ; }
		else if (tri1==ind0) { ind0new=1 ;  }
		else if (tri2==ind0) { ind0new=2 ;  }
		
		
		if (tri0==ind1) { ind1new=0 ; }
		else if (tri1==ind1) { ind1new=1 ; }
		else if (tri2==ind1) { ind1new=2 ; }
		
		return true;
	}
	return false;
	
}

void meshUtils::intersectionPoint(const float & ycut, const float & px0, const float & py0, const float & pz0, const  float & dx, const float & dy, const float & dz, vector<float> & px, vector<float> & py, vector<float> & pz){
	float t=(ycut-py0)/dy;
	px.push_back(px0+dx*t);
	py.push_back(ycut);
	pz.push_back(pz0+dz*t);
}


float meshUtils::myatan2(const float & y, const float & x)
{
	const float PI=3.14159265358979323846;
	float f=atan2(y,x);
	/*float f=atan(y/x);
	if ((y<0)&&(x<0)) 
		f+=PI;
	else if ((y>0)&&(x<0))
		f+=PI;
	else if ((y<0)&&(x>0))
		f+=2*PI;
	*/	
	if (f<0)
		f+=2*PI;
	
	return f;
}

//assume unwrapped coordinates
ReturnMatrix  meshUtils::addSphericalCorrdinates( const Matrix & m1, const Matrix & m2 )
{
	//solves problem with avreage around 2pi/0 boundary
	//assume theta ranges from 0 to 2pi and is the second coordinate
	const float TPI=2*3.14159265358979323846;
		const float PI=3.14159265358979323846;

	Matrix Sum(m1.Nrows(),1);
	
	//we assume that problem are only with quadrant 1 and 4 operations
	for (int i=0;i<m1.Nrows();i+=3)
	{
		Sum.element(i,0)=m1.element(i,0)+m2.element(i,0);
		
		int NcircM1= static_cast<int>(floor(m1.element(i+1,0) / (TPI) ));
		if (NcircM1<0) NcircM1++;
	
		int NcircM2= static_cast<int>(floor(m2.element(i+1,0) / (TPI) ));
		if (NcircM2<0) NcircM2++;



		float remM1=m1.element(i+1,0) -NcircM1*TPI; 
		float remM2=m2.element(i+1,0) - NcircM2*TPI; 

		if (( (remM1>=0)&&(remM1<=PI/2) ) && ( (remM2>=3*PI/2)&&(remM2<=2*PI) ) )
		{
			Sum.element(i+1,0)=m1.element(i+1,0)+(m2.element(i+1,0) - (abs(NcircM2)+1)*2*PI);
			cout<<"case 1 "<<NcircM1<<" "<<NcircM2<<" "<<m1.element(i+1,0)<<" "<<m2.element(i+1,0)<<" "<<m1.element(i+1,0)+m2.element(i+1,0)<<" "<<m1.element(i+1,0)+(m2.element(i+1,0) - 2*PI)<<endl;
//		}else if ( ( (remM2>=0)&&(remM2<=PI/2) ) && ( (remM1>=3*PI/2)&&(remM1<=2*PI) ) ){
		//	Sum.element(i+1,0)=m1.element(i+1,0)-((abs(NcircM1)+1)*2*PI)+m2.element(i+1,0);
	//		cout<<"case 2 "<<NcircM1<<" "<<NcircM2<<" "<<m1.element(i+1,0)<<" "<<m2.element(i+1,0)<<" "<<m1.element(i+1,0)+m2.element(i+1,0)<<" "<<m1.element(i+1,0)+(m2.element(i+1,0) - 2*PI)<<endl;
			
		}else{
			cout<<"case 3 "<<m1.element(i+1,0)<<" "<<m2.element(i+1,0)<<" "<<m1.element(i+1,0)+m2.element(i+1,0)<<endl;
			
			Sum.element(i+1,0)=m1.element(i+1,0)+m2.element(i+1,0);//+m2.element(i+1,0);
		}
	if (Sum.element(i+1,0)<0) 
			Sum.element(i+1,0)+=2*PI;
		
		Sum.element(i+2,0)=m1.element(i+2,0)+m2.element(i+2,0);
	}
	
	Sum.Release();
	return Sum;
	
}

ReturnMatrix meshUtils::subSampleMatrix(const Matrix & m, const vector<bool> & vmask)
{
	unsigned int N=0;
	for (vector<bool>::const_iterator i=vmask.begin();i!=vmask.end();i++)
		if (*i) N++;
		
		
	Matrix subM(N,m.Ncols());
	unsigned int count=1,count2=1;
	for (vector<bool>::const_iterator i=vmask.begin();i!=vmask.end();i++,count++)
	{
		if (*i)
		{
			subM.Row(count2)=m.Row(count);
			count2++;
		}
	}
	subM.Release();
	return subM;
}


ReturnMatrix meshUtils::subSample_Nby1_3D_Matrix(const Matrix & m, const vector<bool> & vmask )
{
	unsigned int N=0;
	for (vector<bool>::const_iterator i=vmask.begin();i!=vmask.end();i++)
		if (*i) N++;
		
	Matrix subM(3*N,m.Ncols());
	unsigned int count=1,count2=1;
	for (vector<bool>::const_iterator i=vmask.begin();i!=vmask.end();i++,count+=3)
		if (*i)
		{		
		cout<<"N "<<count2<<" "<<count<<" "<<m.Nrows()<<endl;
			subM.Row(count2)=m.Row(count);
			subM.Row(count2+1)=m.Row(count+1);
			subM.Row(count2+2)=m.Row(count+2);
			count2+=3;
		}
	subM.Release();
	return subM;
}


/*
ReturnMatrix  meshUtils::averageSphericalCorrdinates( const Matrix & m1, const Matrix & m2 , int & N1, const int & N2)
{
	//theta should NEVER be outside [0,2PI[
	//solves problem with avreage around 2pi/0 boundary
	//assume theta ranges from 0 to 2pi and is the second coordinate
	const float TPI=2*3.14159265358979323846;
		const float PI=3.14159265358979323846;

	Matrix Sum(m1.Nrows(),1);
	
	//we assume that problem are only with quadrant 1 and 4 operations
	int Nnew=N1+N2;
	for (int i=0;i<m1.Nrows();i+=3)
	{
		Sum.element(i,0)=(m1.element(i,0)*N1+m2.element(i,0)*N2)/Nnew;
		Sum.element(i+2,0)=(m1.element(i+2,0)*N1+m2.element(i+2,0)*N2)/Nnew;
		cout<<"PHI "<<m1.element(i+2,0)<<" "<<m2.element(i+2,0)<<endl;


//		int NcircM1= floor(m1.element(i+1,0) / (TPI) );
//		if (NcircM1<0) NcircM1++;
	
//		int NcircM2= floor(m2.element(i+1,0) / (TPI) );
//		if (NcircM2<0) NcircM2++;
		//assumign between 0 and 2PI
		float remM1=m1.element(i+1,0);// -NcircM1*TPI; 
		float remM2=m2.element(i+1,0);// - NcircM2*TPI; 

		if (( (remM1>=0)&&(remM1<=PI/2) ) && ( (remM2>=3*PI/2.0)&&(remM2<2*PI) ) )
		{
			cout<<"HMMMMM "<<3*PI/2.0<<" "<<3*PI/2<<endl;
						Sum.element(i+1,0)=m1.element(i+1,0)*N1+(m2.element(i+1,0) -2*PI)*N2;
			cout<<"case 1 "<<NcircM1<<" "<<NcircM2<<" "<<m1.element(i+1,0)<<" "<<m2.element(i+1,0)<<" "<<m1.element(i+1,0)+m2.element(i+1,0)<<" "<<m1.element(i+1,0)+(m2.element(i+1,0) - 2*PI)<<endl;
		}else if ( ( (remM2>=0)&&(remM2<=PI/2) ) && ( (remM1>=3*PI/2.0)&&(remM1<2*PI) ) ){
		Sum.element(i+1,0)=(m1.element(i+1,0)-2*PI)*N1+m2.element(i+1,0)*N2;
			cout<<"case 2 "<<NcircM1<<" "<<NcircM2<<" "<<m1.element(i+1,0)<<" "<<m2.element(i+1,0)<<" "<<m1.element(i+1,0)+m2.element(i+1,0)<<" "<<m1.element(i+1,0)+(m2.element(i+1,0) - 2*PI)<<endl;
			
		}else{
			cout<<"case 3 "<<m1.element(i+1,0)<<" "<<m2.element(i+1,0)<<" "<<m1.element(i+1,0)+m2.element(i+1,0)<<endl;
			
			Sum.element(i+1,0)=m1.element(i+1,0)*N1+m2.element(i+1,0)*N2;//+m2.element(i+1,0);
		}
		
		Sum.element(i+1,0)/=Nnew;
		if (Sum.element(i+1,0)<0) 
			Sum.element(i+1,0)+=2*PI;
		
	}
	N1=Nnew;

	Sum.Release();
	return Sum;
	
}
*/

 void meshUtils::SVDarcSpherical( Matrix & m1, DiagonalMatrix & D, Matrix & U, Matrix & V)
{


	//takes deviattion from mean as input, arc lengths need be directional
	//converts to arc lengths
	 for (int j=0; j<m1.Ncols();j++){
		 for (int i=0; i<m1.Nrows();i+=3){
			//sign for thetha contain within theta
			m1.element(i+1,j)=m1.element(i,j)*sin(abs(m1.element(i+2,j)))*m1.element(i+1,j);
			m1.element(i+2,j)=m1.element(i,j)*m1.element(i+2,j);
		}
	}
	SVD(m1,D,U,V);
	 
	 for (int j=0; j<U.Ncols();j++){
		 for (int i=0; i<U.Nrows();i+=3){
			 //sign for thetha contain within theta
			 U.element(i+1,j)=U.element(i+1,j)/asin(abs(U.element(i+2,j)))/U.element(i,j);
			 U.element(i+2,j)/=U.element(i,j);//U.element(i+2,j);
		 }
	 }
	 
}
ReturnMatrix  meshUtils::averageSphericalCorrdinates( const Matrix & m1, const Matrix & m2 , int & N1, const int & N2)
{

	//implement the +/- PI deviation, there always a continuous average
	//theta should NEVER be outside [0,2PI[
	//solves problem with avreage around 2pi/0 boundary
	//assume theta ranges from 0 to 2pi and is the second coordinate
//	const float TPI=2*3.14159265358979323846;
		const float PI=3.14159265358979323846;

	Matrix Sum(m1.Nrows(),1);
	
	//we assume that problem are only with quadrant 1 and 4 operations
	int Nnew=N1+N2;
	for (int i=0;i<m1.Nrows();i+=3)
	{
		Sum.element(i,0)=(m1.element(i,0)*N1+m2.element(i,0)*N2)/Nnew;
		Sum.element(i+2,0)=(m1.element(i+2,0)*N1+m2.element(i+2,0)*N2)/Nnew;
		cout<<"PHI "<<m1.element(i+2,0)<<" "<<m2.element(i+2,0)<<endl;


//		int NcircM1= floor(m1.element(i+1,0) / (TPI) );
//		if (NcircM1<0) NcircM1++;
	
//		int NcircM2= floor(m2.element(i+1,0) / (TPI) );
//		if (NcircM2<0) NcircM2++;
		//assumign between 0 and 2PI
		float remM1=m1.element(i+1,0);// -NcircM1*TPI; 
		float remM2=m2.element(i+1,0);// - NcircM2*TPI; 



		if ( ((remM1>=0)&&(remM1<PI)) && (remM2>(remM1+PI)) )
		{
			Sum.element(i+1,0)=m1.element(i+1,0)*N1+(m2.element(i+1,0) -2*PI)*N2;
			cout<<"case 1 "<<" "<<m1.element(i+1,0)<<" "<<m2.element(i+1,0)<<" "<<m1.element(i+1,0)+m2.element(i+1,0)<<" "<<m1.element(i+1,0)+(m2.element(i+1,0) - 2*PI)<<endl;
		}else if ( ((remM2>=0)&&(remM2<PI)) && (remM1>(remM2+PI)) ){
		Sum.element(i+1,0)=(m1.element(i+1,0)-2*PI)*N1+m2.element(i+1,0)*N2;
			cout<<"case 2 "<<" "<<m1.element(i+1,0)<<" "<<m2.element(i+1,0)<<" "<<m1.element(i+1,0)+m2.element(i+1,0)<<" "<<m1.element(i+1,0)+(m2.element(i+1,0) - 2*PI)<<endl;
			
		}else{
			cout<<"case 3 "<<m1.element(i+1,0)<<" "<<m2.element(i+1,0)<<" "<<m1.element(i+1,0)+m2.element(i+1,0)<<endl;
			
			Sum.element(i+1,0)=m1.element(i+1,0)*N1+m2.element(i+1,0)*N2;//+m2.element(i+1,0);
		}
		
		Sum.element(i+1,0)/=Nnew;
		if (Sum.element(i+1,0)<0) 
			Sum.element(i+1,0)+=2*PI;
		
	}
	N1=Nnew;

	Sum.Release();
	return Sum;
	
}

//assume unwrapped coordinates
ReturnMatrix meshUtils::subtractSphericalCoordinates( const Matrix & m1, const Matrix & m2 )
{
	//solves problem with avreage around 2pi/0 boundary
	//solves problem with avreage around 2pi/0 boundary
	//assume theta ranges from 0 to 2pi and is the second coordinate
	const float PI=3.14159265358979323846;
	Matrix Sum(m1.Nrows(),1);
	
	//we assume that problem are only with quadrant 1 and 4 operations
	for (int i=0;i<m1.Nrows();i+=3)
	{
		Sum.element(i,0)=m1.element(i,0)+m2.element(i,0);
		if ( ( (m1.element(i+1,0)>=0)&&(m1.element(i+1,0)<=PI/2) ) && ( (m2.element(i+1,0)>=3*PI/2)&&(m2.element(i+1,0)<=2*PI) ) )
			Sum.element(i+1,0)=m1.element(i+1,0)-(m2.element(i+1,0) - 2*PI);
		else if ( ( (m2.element(i+1,0)>=0)&&(m2.element(i+1,0)<=PI/2) ) && ( (m1.element(i+1,0)>=3*PI/2)&&(m1.element(i+1,0)<=2*PI) ) )
			Sum.element(i+1,0)=(m1.element(i+1,0)-2*PI)-m2.element(i+1,0);
		else
			Sum.element(i+1,0)=m1.element(i+1,0)-m2.element(i+1,0);
		
		if (Sum.element(i+1,0)<0) 
			Sum.element(i+1,0)+=2*PI;
		
		Sum.element(i+2,0)=m1.element(i+2,0)+m2.element(i+2,0);
	}
	
	Sum.Release();
	return Sum;
}



template<class T>
void meshUtils::cartesianToSphericalCoord(vector<T> & verts)
{
	for (unsigned int i=0;i<verts.size();i+=3)
	{
		float r=sqrt(verts.at(i)*verts.at(i)+verts.at(i+1)*verts.at(i+1)+verts.at(i+2)*verts.at(i+2));
		float theta=myatan2(verts.at(i+1),verts.at(i));
		float phi=acos(verts.at(i+2)/r);
		
		verts.at(i)=r;
		verts.at(i+1)=theta;
		verts.at(i+2)=phi;
	}
}
template void meshUtils::cartesianToSphericalCoord<float>(vector<float> & verts);



template<class T>
void meshUtils::sphericalToCartesianCoord(vector<T> & verts)
{
	for (unsigned int i=0;i<verts.size();i+=3)
	{
		verts.at(i)=verts.at(i)*cos(verts.at(i+1))*sin(verts.at(i+2));
		verts.at(i+1)=verts.at(i)*sin(verts.at(i+1))*sin(verts.at(i+2));
		verts.at(i+2)=verts.at(i)*cos(verts.at(i+2));
	}
}
template void meshUtils::sphericalToCartesianCoord<float>(vector<float> & verts);


//these assumer unwrapped coordinates
void meshUtils::getSphericalCoordFromCart(NEWMAT::Matrix & Mr, NEWMAT::Matrix & Mtheta, NEWMAT::Matrix & Mphi)
{
	Mr.ReSize(Points.Nrows(),1);
	Mtheta.ReSize(Points.Nrows(),1);
	Mphi.ReSize(Points.Nrows(),1);

		for (int i=0;i<Points.Nrows();i++)
		{
			float r=sqrt(Points.element(i,0)*Points.element(i,0)+Points.element(i,1)*Points.element(i,1)+Points.element(i,2)*Points.element(i,2));
		//put theta in x/z plane
			float theta=myatan2(Points.element(i,2),Points.element(i,2));
			float phi=acos(Points.element(i,1)/r);
	//		float theta=myatan2(Points.element(i+1,j),Points.element(i,j));
	//		float phi=acos(Points.element(i+2,j)/r);
	
			
			Mr.element(i,0)=r;
			Mtheta.element(i,0)=theta;
			Mphi.element(i,0)=phi;
		}

}


void meshUtils::cartesianToSphericalCoord(Matrix & verts)
{
	for (int j=0;j<verts.Ncols();j++)
		for (int i=0;i<verts.Nrows();i+=3)
		{
			float r=sqrt(verts.element(i,j)*verts.element(i,j)+verts.element(i+1,j)*verts.element(i+1,j)+verts.element(i+2,j)*verts.element(i+2,j));
		//put theta in x/z plane
	//		float theta=myatan2(verts.element(i+2,j),verts.element(i,j));
	//		float phi=acos(verts.element(i+1,j)/r);
			float theta=myatan2(verts.element(i+1,j),verts.element(i,j));
			float phi=acos(verts.element(i+2,j)/r);
			//			cout<<"cartTOspehre "<<verts.element(i,j)<<" "<<verts.element(i+1,j)<<" "<<verts.element(i+2,j)<<" "<<theta<<" "<<theta*180/3.14<<endl;

			verts.element(i,j)=r;
			verts.element(i+1,j)=theta;
			verts.element(i+2,j)=phi;
		}
}
void meshUtils::sphericalToCartesianCoord(Matrix & verts)
{
	for (int j=0;j<verts.Ncols();j++)
		for (int i=0;i<verts.Nrows();i+=3)
		{
			verts.element(i,j)=verts.element(i,j)*cos(verts.element(i+1,j))*sin(verts.element(i+2,j));
			verts.element(i+1,j)=verts.element(i,j)*sin(verts.element(i+1,j))*sin(verts.element(i+2,j));
			verts.element(i+2,j)=verts.element(i,j)*cos(verts.element(i+2,j));
		}
}

template<class T>
void meshUtils::cartesianToSphericalCoord(Mesh & m)
{
	for (vector<Mpoint*>::iterator i=m._points.begin();i!=m._points.end();i++)
	{
		
		
		T x=(*i)->get_coord().X;
		T y=(*i)->get_coord().Y;
		T z=(*i)->get_coord().Z;
		
		T r=sqrt(x*x+y*y+z*z);
		T theta=myatan2(y,x);
		T phi=cos(z/r);
		
		(*i)->_update_coord.X=r;
		(*i)->_update_coord.Y=theta;
		(*i)->_update_coord.Z=phi;
		
		
	}
				m.update();
	
}
template void meshUtils::cartesianToSphericalCoord<float>(Mesh & m);
template void meshUtils::cartesianToSphericalCoord<double>(Mesh & m);



template<class T>
void meshUtils::sphericalToCartesianCoord(Mesh & m)
{
	for (vector<Mpoint*>::iterator i=m._points.begin();i!=m._points.end();i++)
	{
		T r=(*i)->get_coord().X;
		T theta=(*i)->get_coord().Y;
		T phi=(*i)->get_coord().Z;
		//put theta in x/z plane

//		(*i)->_update_coord.X=r*cos(theta)*sin(phi);
//		(*i)->_update_coord.Z=r*sin(theta)*sin(phi);
//		(*i)->_update_coord.Y=r*cos(phi);
		
		(*i)->_update_coord.X=r*cos(theta)*sin(phi);
		(*i)->_update_coord.Y=r*sin(theta)*sin(phi);
		(*i)->_update_coord.Z=r*cos(phi);
	}
	m.update();
}
template void meshUtils::sphericalToCartesianCoord<float>(Mesh & m);
template void meshUtils::sphericalToCartesianCoord<double>(Mesh & m);

	//-----------------------------Transformation Matrix Utilities Begin---------------------------//
	void meshUtils::preMultiplyGlobalRotation(Matrix & fmat, const Matrix & R)
	{
		//assume a 3by 3 rotation matrix
		//SubMatrix command has indices starting from 1
		fmat.SubMatrix(1,3,1,3)=R*fmat.SubMatrix(1,3,1,3);
		fmat.SubMatrix(1,3,4,4)=R*fmat.SubMatrix(1,3,4,4);
	}
	void meshUtils::preMultiplyGlobalScale(Matrix & fmat, const float & s)
	{
		for (int col=0;col<4;col++)
			for (int row=0;row<3;row++)
				fmat.element(row,col)*=s;
	}
	void meshUtils::preMultiplyGlobalScale(Matrix & fmat, const float & sx,const float & sy, const float & sz)
	{
		for (int col=0;col<4;col++)
				fmat.element(0,col)*=sx;
		for (int col=0;col<4;col++)
				fmat.element(1,col)*=sy;
		for (int col=0;col<4;col++)
				fmat.element(2,col)*=sz;
	}
	void meshUtils::preMultiplyTranslation(Matrix & fmat, const  float & tx, const float & ty, const float & tz )
	{
		//must be a 4 by 4 matrix 
		fmat.element(0,3)+=tx;
		fmat.element(1,3)+=ty;
		fmat.element(2,3)+=tz;
	}
	ReturnMatrix meshUtils::getIdentityMatrix(const short N)
	{
		Matrix fmat(N,N);
		fmat=0;
		for (int i=0;i<N;i++)
			fmat.element(i,i)=1;
		fmat.Release();
		return fmat;
	}

		//-----------------------------Tranformation Matrix Utilities End---------------------------//

void meshUtils::combineMeshesWithVectorsAndScalars(const vector<string> & meshlist)
{
	Matrix AllPoints,AllScalars,AllVectors, AllPolygons;
	for (unsigned int i=0;i<meshlist.size();i++)
	{
		loadMesh(meshlist.at(i));
		if (i==0)
		{
			AllPoints=Points;
			AllScalars=Scalars;
			AllVectors=Vectors;
			AllPolygons=Polygons;
		}else
		{
			//shift polygons call must precede pts call
			shiftPolygonMatrix(AllPoints.Nrows());
			AllPolygons=AllPolygons & Polygons;
			AllPoints=AllPoints & Points;
			AllScalars=AllScalars & Scalars;
			AllVectors=AllVectors & Vectors;
		}
	}
	Points=AllPoints;
	Scalars=AllScalars;
	Vectors=AllVectors;
	Polygons=AllPolygons;
	
	AllPoints.Release();
	AllScalars.Release();
	AllVectors.Release();
	AllPolygons.Release();
	
}






	
	void meshUtils::removeRow(Matrix & mat, int ind ){
		//operate with ind from 0 to end-1
		ind++;
		mat=mat.SubMatrix(1,ind-1,1, mat.Ncols()) & mat.SubMatrix(ind+1,mat.Nrows(),1,mat.Ncols());
		
	}


	
	
	ColumnVector meshUtils::sample_image_at_vertices(string meshname, string imagename) { 
		//setup objects
		
		
		Mesh m;
		volume<float> im;
		
		//load data
		m.load(meshname);
		read_volume(im,imagename);
		float dx=im.xdim(), dy=im.ydim(), dz=im.zdim();
		
		//variables
		int nverts=m._points.size();
		ColumnVector Samples(nverts);
		
		int count=0;
		for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ )
		{
			Samples.element(count)=im.interpolate(((*i)->get_coord().X)/dx,((*i)->get_coord().Y)/dy, ((*i)->get_coord().Z)/dz);//im.interpolate((*i)->get_coord().X, (*i)->get_coord().Y, (*i)->get_coord().Z);
			count++;
		}
		
		return Samples;
	}
	
	void meshUtils::sample_image_at_verticesNN(const string & meshname,const string & imagename,const string & outname) 
{ 
	//setup objects
	
	
	Mesh m;
	volume<float> im;
	
	//load data
	m.load(meshname);
	read_volume(im,imagename);
	float dx=im.xdim(), dy=im.ydim(), dz=im.zdim();
	
	//variables
	int nverts=m._points.size();
	ColumnVector Samples(nverts);
	
	int count=0;
	for (vector<Mpoint*>::iterator iv = m._points.begin(); iv!=m._points.end(); iv++ )
	{
		
		float val=im.interpolate(((*iv)->get_coord().X)/dx,((*iv)->get_coord().Y)/dy, ((*iv)->get_coord().Z)/dz);
		if (val!=0){
			Samples.element(count)=val;
		}else{
			float dist=1000;
			for (int i=-3;i<3;i++){
				for (int j=-3;j<3;j++){
					for (int k=-3;k<3;k++){
						float val2=im.interpolate( ((*iv)->get_coord().X+i)/dx,((*iv)->get_coord().Y+j)/dy, ((*iv)->get_coord().Z+k)/dz);
						if ((val2!=0)&&(dist<(i*i+j*j+k*k))){
							Samples.element(count)=val2;
						}
					}
				}
			}	
			
			
		}//im.interpolate((*i)->get_coord().X, (*i)->get_coord().Y, (*i)->get_coord().Z);
		count++;
	}
	fslvtkIO* output = new fslvtkIO;
	output->setMesh(m);
	output->setScalars(Samples);
	output->save(outname);
	
}



void meshUtils::sampleSumAndWrite(const string & inname, const string & inmeshname, const string & outname){
	fslvtkIO* output = new fslvtkIO;
	
	Mesh m;
	m.load(inmeshname);
	output->setMesh(m);
	
	
	
	ifstream fin;
	fin.open(inname.c_str());
	int nSub;
	fin>>nSub;
	ColumnVector meanCV;
	
	
	for (int i=0;i<nSub;i++){
		string mname;
		string iname;
		fin>>mname>>iname;
		cout<<"mname "<<mname<<" "<<iname<<endl;
		ColumnVector scalarsCV;
		scalarsCV=sample_image_at_vertices(mname, iname);
		if (i==0){
			meanCV=scalarsCV;
		}else{
			meanCV+=scalarsCV;
		}
	}
	
	//meanCV/=nSub;
	
	
	output->setScalars(meanCV);
	output->save(outname);
	
	//	ofstream scout;
	//	scout.open(outname.value().c_str());
	//	scout<<"SCALARS"<<" "<<meanCV.Nrows()<<endl;
	//	for (int i=0;i<meanCV.Nrows();i++){
	//		
	//		scout<<meanCV.element(i)<<endl;
	
	//	}
	
	//	scout.close();
	
	
}
template<class T>
void meshUtils::sampleImageAtPoints(const volume<T> & image, vector<T> & vsamples)
{
	vector<short> vmask;
	float dx=image.xdim(), dy=image.ydim(), dz=image.zdim();
//	cout<<"sample image at points "<<Points.Nrows()<<" "<<Points.Ncols()<<endl;
	for (unsigned int i=0; i<static_cast<unsigned int>(Points.Nrows());i++)
	{
			//	vsamples.push_back(static_cast<T>(image.value(Points.element(i,0)/dx, Points.element(i,1)/dy, Points.element(i,2)/dz )));
//cout<<"i "<<i<<endl;
		vsamples.push_back(static_cast<T>(image.interpolate(Points.element(i,0)/dx, Points.element(i,1)/dy, Points.element(i,2)/dz )));
	}
//	cout<<"done sampe image"<<endl;

}


template void meshUtils::sampleImageAtPoints<short>(const volume<short> & image, vector<short> & vsamples);
template void meshUtils::sampleImageAtPoints<float>(const volume<float> & image, vector<float> & vsamples);
template void meshUtils::sampleImageAtPoints<double>(const volume<double> & image, vector<double> & vsamples);


void meshUtils::LQSurfaceReg(const Matrix & refPoints, Matrix & fmat, const int & dof)
{
	//demenad vertices
//	cout<<"POINTS "<<endl;
//	cout<<Points<<endl;
	Matrix PointsDM(Points.Nrows(),Points.Ncols());
	Matrix refPointsDM(refPoints.Nrows(),refPoints.Ncols());

	//calculate centroid to determine translation 
	float mean1x=0,mean1y=0,mean1z=0;
	for (int i =0; i<Points.Nrows();i++){
	  mean1x+=Points.element(i,0);
	  mean1y+=Points.element(i,1);
	  mean1z+=Points.element(i,2);
	}
	mean1x/=Points.Nrows();
	mean1y/=Points.Nrows();
	mean1z/=Points.Nrows();


	//calculate centroid to determine translation 
	float meanRefx=0,meanRefy=0,meanRefz=0;
	for (int i =0; i<refPoints.Nrows();i++){
	  meanRefx+=refPoints.element(i,0);
	  meanRefy+=refPoints.element(i,1);
	  meanRefz+=refPoints.element(i,2);
	}

	meanRefx/=refPoints.Nrows();
	meanRefy/=refPoints.Nrows();
	meanRefz/=refPoints.Nrows();
	//we now have enough to calculate translation

	//demean data 
	//both meshes must have equal number odf vertices
	for (int i =0; i<Points.Nrows();i++){
	  PointsDM.element(i,0)=Points.element(i,0)-mean1x;
	  PointsDM.element(i,1)=Points.element(i,1)-mean1y;
	  PointsDM.element(i,2)=Points.element(i,2)-mean1z;

	  refPointsDM.element(i,0)= refPoints.element(i,0) -meanRefx;
	  refPointsDM.element(i,1)= refPoints.element(i,1) -meanRefy;
	  refPointsDM.element(i,2)= refPoints.element(i,2) -meanRefz;

	}

	//calculate translation component
	float tx=meanRefx - mean1x ;
	float ty=meanRefy - mean1y ;
	float tz=meanRefz - mean1z ;
	
	//calculate scale
	float ssq1[4]={0,0,0,0};
	float ssqRef[4]={0,0,0,0};
	float scale[4]={1,1,1,1};
	
	if (dof==7)
	{
		for (int i =0; i<Points.Nrows();i++)
			ssq1[0]+=PointsDM.element(i,0)*PointsDM.element(i,0) + PointsDM.element(i,1)*PointsDM.element(i,1) + PointsDM.element(i,2)*PointsDM.element(i,2);
		for (int i =0; i<refPoints.Nrows();i++)
			ssqRef[0]+=refPointsDM.element(i,0)*refPointsDM.element(i,0) + refPointsDM.element(i,1)*refPointsDM.element(i,1) + refPointsDM.element(i,2)*refPointsDM.element(i,2) ;
		scale[0]=sqrt(ssqRef[0]/ssq1[0]);
	}else if (dof==9)
	{
		for (int i =0; i<Points.Nrows();i++)
		{
			ssq1[1]+=PointsDM.element(i,0)*PointsDM.element(i,0);
			ssq1[2]+=PointsDM.element(i,1)*PointsDM.element(i,1);
			ssq1[3]+=PointsDM.element(i,2)*PointsDM.element(i,2);
		}
		for (int i =0; i<refPoints.Nrows();i++)
		{
			ssqRef[1]+=refPointsDM.element(i,0)*refPointsDM.element(i,0);
			ssqRef[2]+=refPointsDM.element(i,1)*refPointsDM.element(i,1);
			ssqRef[3]+=refPointsDM.element(i,2)*refPointsDM.element(i,2);
		}
		scale[1]=sqrt(ssqRef[1]/ssq1[1]);
		scale[2]=sqrt(ssqRef[2]/ssq1[2]);
		scale[3]=sqrt(ssqRef[3]/ssq1[3]);

	}
	cout<<"scale calculated "<<scale[0]<<" "<<scale[1]<<" "<<scale[2]<<" "<<scale[3]<<endl;
	//calculate rotation
	Matrix M=refPointsDM.t()*PointsDM;
	cout<<"translation "<<mean1x<<" "<<mean1y<<" "<<mean1z<<" "<<tx<<" "<<ty<<" "<<tz<<endl;
		cout<<"translation "<<meanRefx<<" "<<meanRefy<<" "<<meanRefz<<" "<<tx<<" "<<ty<<" "<<tz<<endl;

	Matrix U;
	DiagonalMatrix D;
	SVD(M.t()*M,D,U);
	//M should always be a 3x3 matrix
	for (int i=0;i<D.Nrows();i++){
	  D.element(i)=1/sqrt(D.element(i));
	}
	cout<<"Rotation Matrix ..."<<endl;
	Matrix R(3,3);
	cout<<D.Nrows()<<" "<<U.Nrows()<<" "<<U.Ncols()<<" "<<M.Nrows()<<" "<<refPointsDM.Nrows()<<" " <<refPointsDM.Ncols()<<endl; 
	R=M*(U*D*U.t());
	//R=R.t();
	cout<<"Rotation Matrix "<<R<<endl;


	fmat=getIdentityMatrix(4);
		cout<<"fmat1 "<<endl;
	cout<<fmat<<endl;
	preMultiplyTranslation(fmat,-mean1x,-mean1y,-mean1z);
	cout<<"fmat2"<<endl;
	cout<<fmat<<endl;
	if (dof==7)
		preMultiplyGlobalScale(fmat,scale[0]);
	else if (dof==9)
		preMultiplyGlobalScale(fmat,scale[1],scale[2],scale[3]);

	cout<<"fmat3"<<endl;
	cout<<fmat<<endl;
	if (dof>3)
		preMultiplyGlobalRotation(fmat,R);	
	cout<<"fmat4 "<<mean1x<<" "<<tx<<endl;
	cout<<fmat<<endl;
	preMultiplyTranslation(fmat,mean1x+tx,mean1y+ty,mean1z+tz);
		cout<<"fmat4"<<endl;
		cout<<fmat<<endl;

}

void meshUtils::findMidPointOfMidSlice(const volume<char> & im, const Matrix & fmat, float & cx, float & cy, float & cz)
{
	int bounds[6]={0,0,0,0,0,0};
	getBounds(bounds,im.xdim(),im.ydim(),im.zdim());

	cy=(bounds[2]+bounds[3])/2 * im.ydim();
	vector<float> verts=sliceMesh(cy); 
	
	cx=0,cz=0;
	
	float xmin=1000000,xmax=-10000, zmin=100000, zmax=-100000;
	vector<float>::iterator i=verts.begin();
	while (i!=verts.end())
	{
		if ( (*i)>xmax ) xmax=(*i);
		if ( (*i)<xmin ) xmin=(*i);
		i++;
		i++;//skip y
	
		if ( (*i)>zmax ) zmax=(*i);
		if ( (*i)<zmin ) zmin=(*i);
		i++;

	}
	cout<<xmin<<" "<<xmax<<" "<<zmin<<" "<<zmax<<endl;
	cx=(xmax+xmin)/2;
	cz=(zmax+zmin)/2;
	cout<<cx<<" "<<cy<<endl;
//this would be centroid		
//		cx+=*i; 
//		i++;
//		i++;//skip y coordinate
//		cz+=*i;
//		i++;
//	}
//	cx/=verts.size()/3;
// cz/=verts.size()/3;

}

vector<float> meshUtils::sliceMesh(const float & ycut)
{
	
	short ind0=-1, ind1=-1, indP=-1;
	vector<short> mask;
	for (int i=0; i<Polygons.Nrows();i++) //for each polygons, check if it intersect cut plane (in z)
	{ //check each line
		if ( checkLine(Points.element(static_cast<int>(Polygons.element(i,0)),1), Points.element(static_cast<int>(Polygons.element(i,1)),1),ycut) ) 
			mask.push_back(i);
		else if ( checkLine(Points.element(static_cast<int>(Polygons.element(i,0)),1), Points.element(static_cast<int>(Polygons.element(i,2)),1),ycut) )
			mask.push_back(i);
		else if ( checkLine(Points.element(static_cast<int>(Polygons.element(i,1)),1), Points.element(static_cast<int>(Polygons.element(i,2)),1),ycut) )
			mask.push_back(i);
	}
	
	

	//this deterimes which line of first polygon is cut
	if (  checkLine(Points.element(static_cast<int>(Polygons.element(mask.at(0),0)),1), Points.element(static_cast<int>(Polygons.element(mask.at(0),1)),1),ycut) ) { 
		indP=mask.at(0); ind0=0; ind1=1;
	}
	else if (  checkLine(Points.element(static_cast<int>(Polygons.element(mask.at(0),0)),1), Points.element(static_cast<int>(Polygons.element(mask.at(0),2)),1),ycut) ) { 
		indP=mask.at(0); ind0=0; ind1=2; 
	}
	else if (  checkLine(Points.element(static_cast<int>(Polygons.element(mask.at(0),1)),1), Points.element(static_cast<int>(Polygons.element(mask.at(0),2)),1),ycut) ) {
		indP=mask.at(0); ind0=1; ind1=2;
	}
	
	//these describe the line
	float px0=Points.element(static_cast<int>(Polygons.element(indP,ind0)),0);
	float dx=Points.element(static_cast<int>(Polygons.element(indP,ind1)),0)-px0;
	float py0=Points.element(static_cast<int>(Polygons.element(indP,ind0)),1);
	float dy=Points.element(static_cast<int>(Polygons.element(indP,ind1)),1)-py0;
	float pz0=Points.element(static_cast<int>(Polygons.element(indP,ind0)),2);
	float dz=Points.element(static_cast<int>(Polygons.element(indP,ind1)),2)-pz0;
	
	//keep track of polygon id
	vector<short> vP;
	vP.push_back(indP);
	
	vector<float> vpx,vpy,vpz; //keep track of intersection points
	intersectionPoint(ycut, px0,py0,pz0,dx,dy,dz,vpx,vpy,vpz);
	
	//found start point now find second in triangle this will define the direction of rotation
	if ((ind0==0)&&(ind1==1)) {
		if ( checkLine(Points.element(static_cast<int>(Polygons.element(indP,0)),1), Points.element(static_cast<int>(Polygons.element(indP,2)),1),ycut)) { 
			ind0=0; ind1=2;
		}
		else if ( checkLine(Points.element(static_cast<int>(Polygons.element(indP,1)),1), Points.element(static_cast<int>(Polygons.element(indP,2)),1),ycut)) { 
			ind0=1; ind1=2;
		}
	}
	else if ((ind0==0)&&(ind1==2)) {
		if (  checkLine(Points.element(static_cast<int>(Polygons.element(indP,0)),1), Points.element(static_cast<int>(Polygons.element(indP,1)),1),ycut)) {
			ind0=0; ind1=1;  
		}
		else if (  checkLine(Points.element(static_cast<int>(Polygons.element(indP,1)),1), Points.element(static_cast<int>(Polygons.element(indP,2)),1),ycut)) {
			ind0=1; ind1=2;
		}
	}
	else if ((ind0==1)&&(ind1==2)) {
		if (  checkLine(Points.element(static_cast<int>(Polygons.element(indP,0)),1), Points.element(static_cast<int>(Polygons.element(indP,1)),1),ycut) ) {  
			ind0=0; ind1=1; 
		}
		else if (  checkLine(Points.element(static_cast<int>(Polygons.element(indP,0)),1), Points.element(static_cast<int>(Polygons.element(indP,2)),1),ycut) ) {  
			ind0=0; ind1=2;  
		}  
	}
	
	px0=Points.element(static_cast<int>(Polygons.element(indP,ind0)),0);
	py0=Points.element(static_cast<int>(Polygons.element(indP,ind0)),1);
	pz0=Points.element(static_cast<int>(Polygons.element(indP,ind0)),2);

	intersectionPoint(ycut, px0,py0,pz0,Points.element(static_cast<int>(Polygons.element(indP,ind1)),0)-px0,Points.element(static_cast<int>(Polygons.element(indP,ind1)),1)-py0,Points.element(static_cast<int>(Polygons.element(indP,ind1)),2)-pz0,vpx,vpy,vpz);
	
	//we now have the point associate with the first triangle now find negihbouring triangle
	int n=0;
	do{
		
		for (vector<short>::iterator i=(mask.begin()+1); i!=mask.end(); i++){
			short ind0new=-1, ind1new=-1;

			if  ( (checkTriangleNeighbour(static_cast<int>(Polygons.element(*i,0)),static_cast<int>(Polygons.element(*i,1)),static_cast<int>(Polygons.element(*i,2)), static_cast<int>(Polygons.element(indP,ind0)),static_cast<int>(Polygons.element(indP,ind1)) , ind0new, ind1new ) ) ){//&& ( notused)){

				indP=*i;
				vP.push_back(indP);				
				mask.erase(i);

				if  ( ((ind0new==0)&&(ind1new==1)) || ((ind0new==1)&&(ind1new==0)) )
				{
					if (  checkLine(Points.element(static_cast<int>(Polygons.element(indP,0)),1), Points.element(static_cast<int>(Polygons.element(indP,2)),1),ycut) )
					{ 
						ind0=0; ind1=2; 
					}
					else if (  checkLine(Points.element(static_cast<int>(Polygons.element(indP,1)),1), Points.element(static_cast<int>(Polygons.element(indP,2)),1),ycut) )
					{ 
						ind0=1; ind1=2; 
					}
				}else if ( ((ind0new==0)&&(ind1new==2)) || ((ind0new==2)&&(ind1new==0)) ) 
				{
					if (  checkLine(Points.element(static_cast<int>(Polygons.element(indP,0)),1), Points.element(static_cast<int>(Polygons.element(indP,1)),1),ycut) )
					{
						ind0=0; ind1=1;  
					}
					else if (  checkLine(Points.element(static_cast<int>(Polygons.element(indP,1)),1), Points.element(static_cast<int>(Polygons.element(indP,2)),1),ycut) )
					{
						ind0=1; ind1=2;  
					}
				}else if ( ((ind0new==1)&&(ind1new==2)) ||  ((ind0new==2)&&(ind1new==1))  )
				{
					if (  checkLine(Points.element(static_cast<int>(Polygons.element(indP,0)),1), Points.element(static_cast<int>(Polygons.element(indP,1)),1),ycut) )
					{  
						ind0=0; ind1=1;  
					}
					else if (  checkLine(Points.element(static_cast<int>(Polygons.element(indP,0)),1), Points.element(static_cast<int>(Polygons.element(indP,2)),1),ycut) )
					{  
						ind0=0; ind1=2;  
					}
				}			

				px0=Points.element(static_cast<int>(Polygons.element(indP,ind0)),0);
				py0=Points.element(static_cast<int>(Polygons.element(indP,ind0)),1);
				pz0=Points.element(static_cast<int>(Polygons.element(indP,ind0)),2);
				
				intersectionPoint(ycut, px0,py0,pz0,Points.element(static_cast<int>(Polygons.element(indP,ind1)),0)-px0,Points.element(static_cast<int>(Polygons.element(indP,ind1)),1)-py0,Points.element(static_cast<int>(Polygons.element(indP,ind1)),2)-pz0,vpx,vpy,vpz);
				break;
			}
		}

		n++;


	}while (mask.size()>1);
	
	vector<float> vVerts;
	vector<float>::iterator j=vpy.begin(); 
	vector<float>::iterator k=vpz.begin();

	for (vector<float>::iterator i=vpx.begin(); i!=vpx.end();i++,j++,k++)
	{
		vVerts.push_back(*i);
		vVerts.push_back(*j);
		vVerts.push_back(*k);
		
	}
	return vVerts;
}

	void meshUtils::sampleMeshProfilesFromImage(const volume<float> & image, const float & sample_interval, const unsigned int & ipp)
	{
		
		const float xdim=image.xdim();
		const float ydim=image.ydim();
		const float zdim=image.zdim();
		
		//read mesh	
		//calculate normals for vertices
		vector<float> shape=getPointsAsVector<float>();
		vector<float> nx,ny,nz, dif;
		first_mesh::normal<float, unsigned int>(shape,first_mesh::findNeighbourTriangles<unsigned int>(getPolygonsAsVectorOfVectors<unsigned int>(),getPointsAsMatrix().Nrows()),getPolygonsAsVectorOfVectors<unsigned int>(),nx,ny,nz);//could possibly optimize further
		
		float inc_x1= sample_interval/xdim;
		float inc_y1= sample_interval/ydim;
		float inc_z1= sample_interval/zdim;
		
		vector<float>::iterator nx_i=nx.begin();
		vector<float>::iterator ny_i=ny.begin();
		vector<float>::iterator nz_i=nz.begin();
		
		for (vector<float>::iterator k = shape.begin(); k!= shape.end(); k+=3,nx_i++,ny_i++, nz_i++)
		{
			float inc_x = (*nx_i) * inc_x1;
			float inc_y = (*ny_i) * inc_y1;
			float inc_z = (*nz_i) * inc_z1;
			
			
			(*k)=(*k)/xdim - (ipp-1)*0.5 * inc_x ;
			(*(k+1))=*(k+1)/ydim - (ipp-1)*0.5 * inc_y;
			(*(k+2))= *(k+2)/zdim - (ipp-1)*0.5 * inc_z;
			
			for (int j=0;j<ipp;j++,(*k)+=inc_x, (*(k+1))+=inc_y, (*(k+2))+=inc_z)
			{
				cout<<"sample "<<*nx_i<<" "<<*ny_i<<" "<<*nz_i<<" "<<*k<<" "<<*(k+1)<<" "<<*(k+2)<<" "<<image.interpolate(*k,*(k+1),*(k+2))<<endl;
				dif.push_back(image.interpolate(*k,*(k+1),*(k+2)));
			}
		}				
		
		cout<<"adding intensityies "<<dif.size()<<endl;
		addFieldData(dif,"intensity","float");
		
		
	}

/*
void meshUtils::meshToContours(const string & imname,const string & meshname,const string & flirtmatname,const string & outputname){
	volume<char> im;
	read_volume(im,imname);
	float xdim=im.xdim(),ydim=im.ydim(),zdim=im.zdim();
	
	fslvtkIO* min = new fslvtkIO;
	min->setDataType(static_cast<fslvtkIO::DataType>(0));
	min->readPolyData(meshname);
	
	
	//register mesh 
	//create contour in MNI space
	Matrix fmat=readFlirtMat(flirtmatname);
	meshReg(min,fmat);
	
	Matrix Points=min->getPointsAsMatrix();
	Matrix Polys=min->getPolygons();
	
	float minX=Points.Column(1).Minimum();
	float minY=Points.Column(2).Minimum();
	float minZ=Points.Column(3).Minimum();
	
	float maxX=Points.Column(1).Maximum();
	float maxY=Points.Column(2).Maximum();
	float maxZ=Points.Column(3).Maximum();
	
	cout<<"points bounds "<<minX<<" "<<minY<<" "<<minZ<<" "<<maxX<<" "<<maxY<<" "<<maxZ<<" "<<endl;
	//will create a contour for each image slices.
	//for (short i=floor(minY/ydim+0.5); i<floor(maxY/ydim+0.5);i++){
	
	//}
	
float testplane=(minY+maxY)/2;
cout<<"test plane "<<testplane<<" "<<Polys.Nrows()<<" "<<Polys.Ncols()<<endl;
//find first polygon to intersect
float ycut=testplane;
short ind0=-1, ind1=-1, indP=-1;
vector<short> mask;
for (int i=0; i<Polys.Nrows();i++){
	//check each line
	//	cout<<"i "<<i<<endl;
	
	if (  checkLine(Points.element(Polys.element(i,0),1), Points.element(Polys.element(i,1),1),ycut) ){ mask.push_back(i) ; }
	else if (  checkLine(Points.element(Polys.element(i,0),1), Points.element(Polys.element(i,2),1),ycut) ){  mask.push_back(i) ; }
	else if (  checkLine(Points.element(Polys.element(i,1),1), Points.element(Polys.element(i,2),1),ycut) ){ mask.push_back(i)  ;  }
	//	cout<<"mask "<<mask.size()<<endl;
	//	if (mask.size()>0) cout<<mask.at(mask.size()-1)<<endl;
}
for (vector<short>::iterator i=mask.begin(); i!=mask.end(); i++){
				cout<<"M "<<Polys.element(*i,0)<<" "<<Polys.element(*i,1)<<" "<<Polys.element(*i,2)<<endl;
	
}

int count=0;
//	for (vector<short>::iterator i=mask.begin(); i!=mask.end(); i++){
//check each line
cout<<"coutn "<<count<<endl;	

if (  checkLine(Points.element(Polys.element(mask.at(0),0),1), Points.element(Polys.element(mask.at(0),1),1),ycut) ){ indP=mask.at(0); ind0=0; ind1=1;  }
else if (  checkLine(Points.element(Polys.element(mask.at(0),0),1), Points.element(Polys.element(mask.at(0),2),1),ycut) ){ indP=mask.at(0); ind0=0; ind1=2; }
else if (  checkLine(Points.element(Polys.element(mask.at(0),1),1), Points.element(Polys.element(mask.at(0),2),1),ycut) ){indP=mask.at(0); ind0=1; ind1=2;  }
//	}

float px0=Points.element(Polys.element(indP,ind0),0);
float dx=Points.element(Polys.element(indP,ind1),0)-px0;
float py0=Points.element(Polys.element(indP,ind0),1);
float dy=Points.element(Polys.element(indP,ind1),1)-py0;
float pz0=Points.element(Polys.element(indP,ind0),2);
float dz=Points.element(Polys.element(indP,ind1),2)-pz0;
cout<<"point "<<px0<<" "<<py0<<" "<<pz0<<" "<<dx<<" "<<dy<<" "<<dz<<" "<<ind0<<" "<<ind1<<" "<<indP<<endl;

vector<short> vP;
vP.push_back(indP);

vector<float> vpx,vpy,vpz;

intersectionPoint(ycut, px0,py0,pz0,dx,dy,dz,vpx,vpy,vpz);


//found start point no find second in triangle this will define the direction of rotation

if ((ind0==0)&&(ind1==1)){
	if (  checkLine(Points.element(Polys.element(indP,0),1), Points.element(Polys.element(indP,2),1),ycut) ){ ind0=0; ind1=2; }
	else if (  checkLine(Points.element(Polys.element(indP,1),1), Points.element(Polys.element(indP,2),1),ycut) ){ ind0=1; ind1=2; }
}else if ((ind0==0)&&(ind1==2)){
	if (  checkLine(Points.element(Polys.element(indP,0),1), Points.element(Polys.element(indP,1),1),ycut) ){ind0=0; ind1=1;  }
	else if (  checkLine(Points.element(Polys.element(indP,1),1), Points.element(Polys.element(indP,2),1),ycut) ){ind0=1; ind1=2;  }
}else if ((ind0==1)&&(ind1==2)){
	if (  checkLine(Points.element(Polys.element(indP,0),1), Points.element(Polys.element(indP,1),1),ycut) ){  ind0=0; ind1=1;  }
	else if (  checkLine(Points.element(Polys.element(indP,0),1), Points.element(Polys.element(indP,2),1),ycut) ){  ind0=0; ind1=2;  }
}	


px0=Points.element(Polys.element(indP,ind0),0);
//	 dx=Points.element(Polys.element(indP,ind1),0)-px0;
py0=Points.element(Polys.element(indP,ind0),1);
//	 dy=Points.element(Polys.element(indP,ind1),1)-py0;
pz0=Points.element(Polys.element(indP,ind0),2);
//	 dz=Points.element(Polys.element(indP,ind1),2)-pz0;
cout<<"point "<<px0<<" "<<py0<<" "<<pz0<<" "<<ind0<<" "<<ind1<<endl;
intersectionPoint(ycut, px0,py0,pz0,Points.element(Polys.element(indP,ind1),0)-px0,Points.element(Polys.element(indP,ind1),1)-py0,Points.element(Polys.element(indP,ind1),2)-pz0,vpx,vpy,vpz);

//intersectionPoint(ycut, px0,py0,pz0,dx,dy,dz,vpx,vpy,vpz);
cout<<"maks: "<<mask.size()<<endl;

//	mask.erase(mask.begin());
cout<<"maks: "<<mask.size()<<endl;





//we now have the point associate with the first triangle now find negihbouring triangle
int n=0;
do{
	cout<<"n "<<indP<<" "<<vP.front()<<" "<<mask.size()<<endl;
	cout<<ycut<<" "<<Points.element(Polys.element(indP,ind1),1)<<" "<<Points.element(Polys.element(indP,ind0),1)<<" "<<Polys.element(indP,ind1)<<" "<<Polys.element(indP,ind0)<<" "<<endl;
	
	for (vector<short>::iterator i=mask.begin(); i!=mask.end(); i++){
		//cout<<"HMM "<<ind0<<" "<<ind1<<endl;
		//		cout<<ycut<<" "<<*i<<" "<<Points.element(Polys.element(indP,ind1),1)<<" "<<Points.element(Polys.element(indP,ind0),1)<<" "<<Polys.element(indP,ind1)<<" "<<Polys.element(indP,ind0)<<" "<<Polys.element(*i,0)<<" "<<Polys.element(*i,1)<<" "<<Polys.element(*i,2)<<endl;
		cout<<ycut<<" "<<*i<<" "<<Polys.element(indP,ind1)<<" "<<Polys.element(indP,ind0)<<" "<<Polys.element(*i,0)<<" "<<Polys.element(*i,1)<<" "<<Polys.element(*i,2)<<endl;
		
		short ind0new=-1, ind1new=-1;
		bool notused=true;
		for (vector<short>::iterator pit=vP.begin(); pit!=vP.end(); pit++){
			if ((*pit)==*i) notused=false;
		}
								cout<<"INDP "<<indP<<" "<<*i<<" "<<vP.back()<<" "<<notused<<endl;
		
		if  ( (checkTriangleNeighbour(Polys.element(*i,0),Polys.element(*i,1),Polys.element(*i,2), Polys.element(indP,ind0),Polys.element(indP,ind1) , ind0new, ind1new ) )
			  && ( notused)){
			cout<<"ENTERED "<<endl;
			cout<<"mask "<<*i<<" "<<mask.size()<<" "<<Polys.element(indP,ind0)<<" "<<Polys.element(indP,ind1)<<endl;
			indP=*i;
			vP.push_back(indP);
			
			//mask.erase(i);
			
			//							cout<<"mask "<<mask.size()<<endl;
			
			for (vector<short>::iterator i=mask.begin(); i!=mask.end(); i++){
				cout<<"M "<<Polys.element(*i,0)<<" "<<Polys.element(*i,1)<<" "<<Polys.element(*i,2)<<endl;
				
			}
			
			
			//need to get right edge of triangle
			//	ind0=ind0new;
			//	ind1=ind1new;
			if  ( ((ind0new==0)&&(ind1new==1)) || ((ind0new==1)&&(ind1new==0)) ){
				if (  checkLine(Points.element(Polys.element(indP,0),1), Points.element(Polys.element(indP,2),1),ycut) ){ ind0=0; ind1=2; }
				else if (  checkLine(Points.element(Polys.element(indP,1),1), Points.element(Polys.element(indP,2),1),ycut) ){ ind0=1; ind1=2; }
			}else if ( ((ind0new==0)&&(ind1new==2)) || ((ind0new==2)&&(ind1new==0)) ) {
				if (  checkLine(Points.element(Polys.element(indP,0),1), Points.element(Polys.element(indP,1),1),ycut) ){ind0=0; ind1=1;  }
				else if (  checkLine(Points.element(Polys.element(indP,1),1), Points.element(Polys.element(indP,2),1),ycut) ){ind0=1; ind1=2;  }
			}else if ( ((ind0new==1)&&(ind1new==2)) ||  ((ind0new==2)&&(ind1new==1))  ){
				if (  checkLine(Points.element(Polys.element(indP,0),1), Points.element(Polys.element(indP,1),1),ycut) ){  ind0=0; ind1=1;  }
				else if (  checkLine(Points.element(Polys.element(indP,0),1), Points.element(Polys.element(indP,2),1),ycut) ){  ind0=0; ind1=2;  }
			}				
			cout<<"ind01 "<<ind0<<" "<<ind1<<endl;
			px0=Points.element(Polys.element(indP,ind0),0);
			py0=Points.element(Polys.element(indP,ind0),1);
			pz0=Points.element(Polys.element(indP,ind0),2);
			
			intersectionPoint(ycut, px0,py0,pz0,Points.element(Polys.element(indP,ind1),0)-px0,Points.element(Polys.element(indP,ind1),1)-py0,Points.element(Polys.element(indP,ind1),2)-pz0,vpx,vpy,vpz);
			break;
		}
	}
	//	cout<<"Indp "<<indP<<endl;
				n++;
	
	cout<<"VP "<<vP.size()<<endl;
	
}while (n<70);//(indP!=vP.front());


//	px.push_back();





}
*/



void meshUtils::draw_segment(volume<short>& image, const Pt& p1, const Pt& p2, int label)
{
	double xdim = (double) image.xdim();
	double ydim = (double) image.ydim();
	double zdim = (double) image.zdim();
	
	//in new version of bet2
	double mininc = min(xdim,min(ydim,zdim)) * .5;
	
	
	
	Vec n = (p1 - p2);
	double d = n.norm();
	n.normalize();
	//	double l = d*4;
	for (double i=0; i<=d; i+=mininc)
    {
		Pt p = p2 + i* n;
		image((int) floor((p.X)/xdim +.5),(int) floor((p.Y)/ydim +.5),(int) floor((p.Z)/zdim +.5)) = label;
    }
}


volume<short> meshUtils::draw_mesh(const volume<short>& image, const Mesh &m, int label)
{
	
	double xdim = (double) image.xdim();
	double ydim = (double) image.ydim();
	double zdim = (double) image.zdim();
	
	//in new version of bet2
	double mininc = min(xdim,min(ydim,zdim)) * .5;
	
	
	volume<short> res = image;
	for (list<Triangle*>::const_iterator i = m._triangles.begin(); i!=m._triangles.end(); i++)
    {
		Vec n = (*(*i)->get_vertice(0) - *(*i)->get_vertice(1));
		double d = n.norm();
		n.normalize();
	       
		
		for (double j=0; j<=d ;  j+=mininc)
		{
			Pt p = (*i)->get_vertice(1)->get_coord()  + (double)j* n;
			draw_segment(res, p, (*i)->get_vertice(2)->get_coord(),label);
		} 
    }
	return res;
}

volume<short> meshUtils::make_mask_from_mesh(const volume<float> & image, const Mesh& m, int label, int* bounds, const bool & sep_boundary)
{
	
	float xdim = (float) image.xdim();
	float ydim = (float) image.ydim();
	float zdim = (float) image.zdim();
	
	volume<short> mask;
	copyconvert(image,mask);
	
	
	mask = 0;
	if (sep_boundary)
		mask = draw_mesh(mask, m,label+100);
	else
		mask = draw_mesh(mask, m,label);

	
	// THIS EXCLUDEDS THE ACTUAL MESH
	volume<short> otl=mask;
	getBounds(m,bounds,xdim,ydim,zdim);
	vector<Pt> current;
	current.clear();
	Pt c(bounds[0]-2, bounds[2]-2, bounds[4]-2);
	
	mask.value(static_cast<int>(c.X),static_cast<int>(c.Y),static_cast<int>(c.Z)) = label;
	current.push_back(c);
	int fillCount=0;
	while (!current.empty())
	{
		Pt pc = current.back();
		int x, y, z;
		x=(int) pc.X; y=(int) pc.Y; z=(int) pc.Z;
		
		current.pop_back();
		fillCount++;
		
		
		if (bounds[0]<=x-1 && mask.value(x-1, y, z)==0) {
			mask.value(x-1, y, z) = label;
			current.push_back(Pt(x-1, y, z));
		}
		if (bounds[2]<=y-1 && mask.value(x, y-1, z)==0) {
			mask.value(x, y-1, z) = label;
			current.push_back(Pt(x, y-1, z));
		}
		if (bounds[4]<=z-1 && mask.value(x, y, z-1)==0) {
			mask.value(x, y, z-1) = label;
			current.push_back(Pt(x, y, z-1));
		}
		if (bounds[1]>=x+1 && mask.value(x+1, y, z)==0){
			mask.value(x+1, y, z) = label;
			current.push_back(Pt(x+1, y, z));
		}
		if (bounds[3]>=y+1 && mask.value(x, y+1, z)==0){
			mask.value(x, y+1, z) = label;
			current.push_back(Pt(x, y+1, z));
		}
		if (bounds[5]>=z+1 && mask.value(x, y, z+1)==0){
			mask.value(x, y, z+1) = label;
			current.push_back(Pt(x, y, z+1)); 
		}
		
	}
	for (int i=bounds[0];i<bounds[1];i++){
		for (int j=bounds[2];j<bounds[3];j++){
			for (int k=bounds[4];k<bounds[5];k++){
				if (mask.value(i,j,k)==0){
					otl.value(i,j,k)=label;
				}
			}
		}
	}
	
	return otl;
	//return mask;
	
}



void meshUtils::fillMesh(const string & imname, const string & meshname, const string & outname, const int & label, const bool & sep_bound ){
	volume<float> im;
	read_volume(im,imname);
	
	Mesh m;
	m.load(meshname);
	
	int bounds[6];
	
	volume<short> imout;
	imout=make_mask_from_mesh(im, m, label, bounds,sep_bound);
	
	save_volume(imout,outname);
	
}


double meshUtils::maxScalar()
{
	double max=-10e-6;
	for (unsigned int i=0; i<static_cast<unsigned int>(Scalars.Nrows());i++)
	{
	//  cout<<"i "<<i<<" "<<Scalars.Nrows()<<" "<<Scalars.Ncols()<<endl;
	  if (Scalars.element(i,0)>max)
			max=Scalars.element(i,0);
	//		cout<<"max "<<max<<endl;

	}	
	return max;		
}


double meshUtils::meanScalar()
{
	double mean=0;
	for (unsigned int i=0; i<(unsigned)Scalars.Nrows();i++)
	{
		mean+=Scalars.element(i,0);
	
	//	cout<<"mean "<<mean<<endl;
	}
	mean/=Scalars.Nrows();
			//cout<<mean<<endl;
	return mean;
}



bool meshUtils::findVertex(const Matrix & vert1,const Matrix & vert2, int ind1 ){
	for (int i=0;i<vert2.Nrows();i++){
		if ((vert1.element(ind1,0)==vert2.element(i,0)) && \
			(vert1.element(ind1,1)==vert2.element(i,1)) && \
			(vert1.element(ind1,2)==vert2.element(i,2)) ) {
			cout<<"found common point "<<" "<<endl;
			return true;
		}
	}
	
	return false;
}
void meshUtils::combinedSharedBoundaries(const string & mname1, const string & mname2 ){
	//create mesh reader
	fslvtkIO* mesh1 = new fslvtkIO;
	mesh1->setDataType(static_cast<fslvtkIO::DataType>(0));
	mesh1->readPolyData(mname1);
	
	fslvtkIO* mesh2 = new fslvtkIO;
	mesh2->setDataType(static_cast<fslvtkIO::DataType>(0));
	mesh2->readPolyData(mname2);
	
	Matrix vert1=mesh1->getPointsAsMatrix();
	Matrix vert2=mesh2->getPointsAsMatrix();
	
	Matrix Poly1=mesh1->getPolygons();
	Matrix Poly2=mesh2->getPolygons();
	vector<short> vind1,vind2;
	bool foundM12=false;
	for (int i=0; i<vert1.Nrows();i++){
		//finds the common points but need to keep common points that are connect to 2 separate vertices
		for (int j=0; j<vert2.Nrows();j++){
			if ((vert1.element(i,0)==vert2.element(j,0)) && \
				(vert1.element(i,1)==vert2.element(j,1)) && \
				(vert1.element(i,2)==vert2.element(j,2))) {
				cout<<"found common point "<<" "<<endl;
				
				for (int k=0;k<Poly1.Nrows();k++){
					if (Poly1.element(k,0)==i){
						if ( (findVertex(vert1,vert2,static_cast<int>(Poly1.element(k,1)))) || (findVertex(vert1,vert2,static_cast<int>(Poly1.element(k,2))))){
							foundM12=true;
						}
					}else if(Poly1.element(k,1)==i){
						if ( (findVertex(vert1,vert2,static_cast<int>(Poly1.element(k,0)))) || (findVertex(vert1,vert2,static_cast<int>(Poly1.element(k,2))))){
							foundM12=true;
						}
					}else if(Poly1.element(k,2)==i){
						if ( (findVertex(vert1,vert2,static_cast<int>(Poly1.element(k,0)))) || (findVertex(vert1,vert2,static_cast<int>(Poly1.element(k,1))))){
							foundM12=true;
						}
					}
				}
				if (!foundM12){
					vind1.push_back(i);
					vind2.push_back(j);
				}
			}
			
		}
	}
	
	vector<short> vindP1;
	vector<short> vindP2;
	int vindc=0;
	for ( int i=0;i<Poly1.Nrows();i++){
		for (int j=0;j<3;j++){
			if (vind1.at(vindc)==Poly1.element(i,j)){ vindP1.push_back(i); break; }
		}
		vindc++;
	}
	
	vindc=0;
	for ( int i=0;i<Poly2.Nrows();i++){
		for (int j=0;j<3;j++){
			if (vind2.at(vindc)==Poly2.element(i,j)){ vindP2.push_back(i); break; }
		}
		vindc++;
	}
	
	Matrix vert1new(vert1.Nrows()-vind1.size(),3);
	Matrix vert2new(vert2.Nrows()-vind2.size(),3);
	Matrix Poly1new(Poly1.Nrows()-vindP1.size(),3);
	Matrix Poly2new(Poly2.Nrows()-vindP2.size(),3);
	
	
	int count=0;
	vindc=0;
	for (int i=0; i<vert1.Nrows();i++){
		if (vind1.at(vindc)!=i){
			vert1new.element(count,0)=vert1.element(i,0);
			vert1new.element(count,1)=vert1.element(i,1);
			vert1new.element(count,2)=vert1.element(i,2);
		}else{
			vindc++;
		}
	}
	count=0;
	vindc=0;
	for (int i=0; i<Poly1new.Nrows();i++){
		if (vindP1.at(vindc)!=i){
			Poly1new.element(count,0)=vert1.element(i,0);
			Poly1new.element(count,1)=vert1.element(i,1);
			Poly1new.element(count,2)=vert1.element(i,2);
		}else{
			vindc++;
		}
	}
	
	count=0;
	vindc=0;
	for (int i=0; i<vert2.Nrows();i++){
		if (vind2.at(vindc)!=i){
			vert2new.element(count,0)=vert2.element(i,0);
			vert2new.element(count,1)=vert2.element(i,1);
			vert2new.element(count,2)=vert2.element(i,2);
		}else{
			vindc++;
		}
	}
	
	fslvtkIO* mout1 = new fslvtkIO;
	mout1->setDataType(static_cast<fslvtkIO::DataType>(0));
	mout1->setPoints(vert1new);
	mout1->setPolygons(vert2new);
	
	
}

void meshUtils::labelAndCombineSharedBoundaries(const string & mname1, const string & mname2, const string & mout1name ){
	//create mesh reader
	fslvtkIO* mesh1 = new fslvtkIO;
	mesh1->setDataType(static_cast<fslvtkIO::DataType>(0));
	mesh1->readPolyData(mname1);
	
	fslvtkIO* mesh2 = new fslvtkIO;
	mesh2->setDataType(static_cast<fslvtkIO::DataType>(0));
	mesh2->readPolyData(mname2);
	
	Matrix vert1=mesh1->getPointsAsMatrix();
	Matrix vert2=mesh2->getPointsAsMatrix();
	
	Matrix Poly1=mesh1->getPolygons();
	Matrix Poly2=mesh2->getPolygons();
	vector<short> vind1,vind2;
	//bool foundM12=false;
	Matrix Sc1(vert1.Nrows(),1);
	int count=0;
	for (int i=0; i<vert1.Nrows();i++){
		//finds the common points but need to keep common points that are connect to 2 separate vertices
		bool found=false;
		for (int j=0; j<vert2.Nrows();j++){
			if ((vert1.element(i,0)==vert2.element(j,0)) && \
				(vert1.element(i,1)==vert2.element(j,1)) && \
				(vert1.element(i,2)==vert2.element(j,2))) {
				cout<<"found common point "<<count<<endl;
				found=true;
				vind1.push_back(i);
				vind2.push_back(j);
				count++;
				break;
			}
		}
		if (found) { Sc1.element(i,0)=2; } else { Sc1.element(i,0)=0; }
	}
	Matrix Sc2(vert2.Nrows()-vind2.size(),1);
	Sc2=1;
	cout<<"sixes "<<vind1.size()<<" "<<vind2.size()<<endl;
	//remove vertices
	//change reference in polygonsand shift polygons matuirx
	int shift=vert1.Nrows();
	shiftPolygonMatrix(Poly2, shift);
	for (unsigned int i=0;i<vind2.size();i++){
		for ( int k=0;k<Poly2.Nrows();k++){
			for ( int l=0;l<Poly2.Ncols();l++){
				if ( Poly2.element(k,l) == vind2.at(i)+shift) {Poly2.element(k,l)=vind1.at(i); }
			}
		}
	}
	
	
	
	for (unsigned int i=0;i<vind2.size();i++){
		removeRow(vert2,vind2.at(i));
		for (unsigned int j=0; j<vind2.size();j++)
			if (vind2.at(j) > vind2.at(i)) { vind2.at(j)--; }//vind2 will store the indices of the altered vertices
				
				for ( int k=0;k<Poly2.Nrows();k++){
					for ( int l=0;l<Poly2.Ncols();l++){
						if ( Poly2.element(k,l) > (vind2.at(i)+shift)) { Poly2.element(k,l)--; }
					}
				}
				
	}
		
		/*	for (unsigned int i=0;i<vind2.size();i++){
			for ( int j=0;j<Poly2.Nrows();j++){
				for ( int k=0;k<Poly2.Ncols();k++){
					if ( Poly2.element(j,k) == vind2.at(i)+shift) { Poly2.element(j,k)=vind1.at(i); }
				}
			}
		}
*/	
		//	mesh1->setPoints(vert2);
		//	mesh1->setPolygons(Poly2);
		
		
		
		mesh1->setPoints(vert1 & vert2);
		mesh1->setPolygons(Poly1 & Poly2);
		mesh1->setScalars(Sc1 & Sc2);
		mesh1->save(mout1name);
		
}


 ReturnMatrix meshUtils::appendSharedBoundaryMask(const Matrix & Points2)
 {
	//create mesh reader
	//m1 is the template first styrcture
	//m2 is the boundary fro which we are searching for shared boundaries
	//mbase is of the same topology as m1, on which we are placing the scalars
	
	Matrix Sc1(Points.Nrows(),1);
	for (int i=0; i<Points.Nrows();i++){
		//finds the common points but need to keep common points that are connect to 2 separate vertices
		bool found=false;
		short index=-1;
		for (int j=0; j<Points2.Nrows();j++)
		{
			if ((Points.element(i,0)==Points2.element(j,0)) && (Points.element(i,1)==Points2.element(j,1)) && (Points.element(i,2)==Points2.element(j,2))) 
			{
				found=true;
				index=j;
				break;
			}
		}
		Sc1.element(i,0)=index;
	}	
	Sc1.Release();
	return Sc1;
}

void meshUtils::sampleSharedBoundaryByMask(const Matrix & Points2)
{
	for (int i=0;i<Scalars.Nrows();i++)
	{
		if (Scalars.element(i,0)>=0)
			Points.Row(i+1)=Points2.Row(static_cast<int>(Scalars.element(i,0))+1);
	}
}

/*
void meshUtils::appendSharedBoundaryMask(const string & mname1, const string & mname2,const string & mbase, const string & mout1name, const bool & indexed, const bool & useSc2 ){
	//create mesh reader
	//m1 is the template first styrcture
	//m2 is the boundary fro which we are searching for shared boundaries
	//mbase is of the same topology as m1, on which we are placing the scalars
	
	fslvtkIO* mesh1 = new fslvtkIO;
	mesh1->setDataType(static_cast<fslvtkIO::DataType>(0));
	mesh1->readPolyData(mname1);
	
	fslvtkIO* mesh2 = new fslvtkIO;
	mesh2->setDataType(static_cast<fslvtkIO::DataType>(0));
	mesh2->readPolyData(mname2);
	
	fslvtkIO* meshB = new fslvtkIO;
	meshB->setDataType(static_cast<fslvtkIO::DataType>(0));
	meshB->readPolyData(mbase);
	
	
	
	Matrix vert1=mesh1->getPointsAsMatrix();
	Matrix vert2=mesh2->getPointsAsMatrix();
	
	Matrix Poly1=mesh1->getPolygons();
	Matrix Poly2=mesh2->getPolygons();
	vector<short> vind1,vind2;
	//bool foundM12=false;
	Matrix Sc1(vert1.Nrows(),1);
	Matrix Sc2(vert2.Nrows(),1);
	Sc2=0;
	int count=0;
	for (int i=0; i<vert1.Nrows();i++){
		//finds the common points but need to keep common points that are connect to 2 separate vertices
		bool found=false;
		for (int j=0; j<vert2.Nrows();j++){
			if ((vert1.element(i,0)==vert2.element(j,0)) && \
				(vert1.element(i,1)==vert2.element(j,1)) && \
				(vert1.element(i,2)==vert2.element(j,2))) {
				cout<<"found common point "<<count<<endl;
				found=true;
				vind1.push_back(i);
				vind2.push_back(j);
				if (indexed){ count++; }
				break;
			}
		}
		if (found) { Sc1.element(i,0)=count; Sc2.element(vind2.back(),0)=count; } else { Sc1.element(i,0)=0; }
	}
	
	if (useSc2){
		meshB->setScalars(Sc2);
	}else{
		meshB->setScalars(Sc1);
	}
	
	meshB->save(mout1name);
	
	mesh2->setScalars(Sc2);
	mesh2->save("st2_"+mout1name);
	
	
}
*/
/*
void meshUtils::applyFlirtThenSBmask(const string & mname1, const string & mname2,const string & mflirtname, const string & mout1name){
	//create mesh reader
	//m1 is the template first styrcture
	//m2 is the boundary fro which we are searching for shared boundaries
	//mbase is of the same topology as m1, on which we are placing the scalars
	
	fslvtkIO* mesh1 = new fslvtkIO;
	mesh1->setDataType(static_cast<fslvtkIO::DataType>(0));
	mesh1->readPolyData(mname1);
	
	Matrix fmat=readFlirtMat(mflirtname);
	fmat=fmat.i();
	meshReg(mesh1,fmat);
	
	
	fslvtkIO* mesh2 = new fslvtkIO;
	mesh2->setDataType(static_cast<fslvtkIO::DataType>(0));
	mesh2->readPolyData(mname2);
	
	
	
	Matrix vert1=mesh1->getPointsAsMatrix();
	Matrix vert2=mesh2->getPointsAsMatrix();
	
	//	Matrix Poly1=mesh1->getPolygons();
	//	Matrix Poly2=mesh2->getPolygons();
	//	vector<short> vind1,vind2;
	//	bool foundM12=false;
	Matrix Sc1 = mesh1->getScalars();
	Matrix Sc2 = mesh2->getScalars();
	
	//	Matrix Sc2(vert2.Nrows(),1);
	//	Sc2=0;
	//int count=0;
	for (int i=0; i<vert1.Nrows();i++){
		//finds the common points but need to keep common points that are connect to 2 separate vertices
		if (Sc1.element(i,0)>0){//search for coresponding vertex index
			for (int j=0; j<Sc2.Nrows();j++){
				if (Sc2.element(j,0)==Sc1.element(i,0)){ 
					vert1.element(i,0)=vert2.element(j,0); 
					vert1.element(i,1)=vert2.element(j,1); 
					vert1.element(i,2)=vert2.element(j,2); 
					break; 
				}
			}
			
		}
	}
	
	mesh1->setPoints(vert1);
	mesh1->save(mout1name);
	
}
*/


void meshUtils::do_work_uncentreMesh(const string & inim, const string & inmesh, const string & outmesh){
	//this function is working directly on volumes
  	volume<float> im;
	read_volume(im,inim);
	float cx=(im.xsize()-1)/2.0*abs(im.xdim());
	float cy=(im.ysize()-1)/2.0*abs(im.ydim());
	float cz=(im.zsize()-1)/2.0*abs(im.zdim());
	
	
  	fslvtkIO* mesh = new fslvtkIO;
	mesh->setDataType(static_cast<fslvtkIO::DataType>(0));
	mesh->readPolyData(inmesh);
	
	Matrix MMesh = mesh->getPointsAsMatrix();
	shift3DVertexMatrix(MMesh,cx,cy,cz);
	
	mesh->setPoints(MMesh);
	mesh->save(outmesh);
	
	
	delete mesh;
	
	
	
}

void meshUtils::do_work_MeshReg(const string & inim, const string & inmesh, const string & outmesh){
	//this function is working directly on volumes
  	volume<float> im;
	read_volume(im,inim);
	float cx=(im.xsize()-1)/2.0*abs(im.xdim());
	float cy=(im.ysize()-1)/2.0*abs(im.ydim());
	float cz=(im.zsize()-1)/2.0*abs(im.zdim());
	
	
  	fslvtkIO* mesh = new fslvtkIO;
	mesh->setDataType(static_cast<fslvtkIO::DataType>(0));
	mesh->readPolyData(inmesh);
	
	Matrix MMesh = mesh->getPointsAsMatrix();
	shift3DVertexMatrix(MMesh,cx,cy,cz);
	
	mesh->setPoints(MMesh);
	mesh->save(outmesh);
	
	
	delete mesh;
	
	
	
}

 void meshUtils::subtractMeshes(const string & mesh1name, const string & mesh2name, const string & out){
	
  	fslvtkIO* mesh1 = new fslvtkIO;
	mesh1->setDataType(static_cast<fslvtkIO::DataType>(0));
	mesh1->readPolyData(mesh1name);
	
	fslvtkIO* mesh2 = new fslvtkIO;
	mesh2->setDataType(static_cast<fslvtkIO::DataType>(0));
	mesh2->readPolyData(mesh2name);
	
	Matrix Diff=mesh2->getPointsAsMatrix() -  mesh1->getPointsAsMatrix();
	mesh1->setVectors(Diff);
	mesh1->save(out);
	
}

void meshUtils::warpMeshWithDefField(const string & fieldname, const string & meshname, const string & meshoutname, const float & dx, const float & dy, const float & dz){
	
	fslvtkIO* ugrid = new fslvtkIO;
	ugrid->setDataType(static_cast<fslvtkIO::DataType>(0));
	ugrid->readPolyData(meshname);
	
	Matrix Points=ugrid->getPointsAsMatrix();
	//propogate volumetric mesh through deformation field //and sample intensities
	
	
	//cout<<"done flirt"<<endl;
	volume4D<float> defField;
	read_volume4D(defField, fieldname);
	
	float resx=defField[0].xdim();
	float resy=defField[0].ydim();
	float resz=defField[0].zdim();
	
	for (int i=0;i<Points.Nrows();i++){
		//shifts are for changes between roi and full image
		float px=Points.element(i,0);
		float py=Points.element(i,1);
		float pz=Points.element(i,2);
		Points.element(i,0)+=defField[0].interpolate(px/resx,py/resy,pz/resz) - dx;
		Points.element(i,1)+=defField[1].interpolate(px/resx,py/resy,pz/resz) - dy;
		Points.element(i,2)+=defField[2].interpolate(px/resx,py/resy,pz/resz) - dz;
		
	}
	
	ugrid->setPoints(Points);
	ugrid->save(meshoutname);
	
}

template< class T >
void meshUtils::warpGridWithDefField(const volume4D<T> & defField, const float & dx, const float & dy, const float & dz){
	
	float resx=defField[0].xdim();
	float resy=defField[0].ydim();
	float resz=defField[0].zdim();
	
	for (int i=0;i<Points.Nrows();i++){
		//shifts are for changes between roi and full image
		float px=Points.element(i,0);
		float py=Points.element(i,1);
		float pz=Points.element(i,2);
		
		//cout<<"points "<<px<<" "<<py<<" "<<pz<<endl;
		Points.element(i,0)=px+defField[0].interpolate(px/resx,py/resy,pz/resz) + dx;
		Points.element(i,1)=py+defField[1].interpolate(px/resx,py/resy,pz/resz) + dy;
		Points.element(i,2)=pz+defField[2].interpolate(px/resx,py/resy,pz/resz) + dz;
	}
		
}

template void meshUtils::warpGridWithDefField<float>(const volume4D<float> & fieldname, const float & dx, const float & dy, const float & dz);
template void meshUtils::warpGridWithDefField<double>(const volume4D<double> & fieldname, const float & dx, const float & dy, const float & dz);

	
	template< class T >
	void meshUtils::warpGridWithDefField(const volume4D<T> & defField, vector<float> & points_in, float warpSc, const float & dx, const float & dy, const float & dz){
		float resx=defField[0].xdim();
		float resy=defField[0].ydim();
		float resz=defField[0].zdim();
		
		
		for (vector<float>::iterator i=points_in.begin();i!=points_in.end();i+=3)
		{
			//shifts are for changes between roi and full image
			float px=*i;
			float py=*(i+1);
			float pz=*(i+2);
			
			*i		= px + warpSc * defField[0].interpolate(px/resx,py/resy,pz/resz) + dx; 
			*(i+1)	= py + warpSc * defField[1].interpolate(px/resx,py/resy,pz/resz) + dy;
			*(i+2)	= pz + warpSc * defField[2].interpolate(px/resx,py/resy,pz/resz) + dz;
			
		}
	}
	
	template void meshUtils::warpGridWithDefField<float>(const volume4D<float> & fieldname, vector<float> & points_warped,float warpSc,const float & dx, const float & dy, const float & dz);
	template void meshUtils::warpGridWithDefField<double>(const volume4D<double> & fieldname, vector<float> & points_warped,float warpSc,const float & dx, const float & dy, const float & dz);
	
	
//template<class T>
/*void meshUtils::warpMeshWithDefField( Mesh & m, const volume4D<float> & defField, const Matrix & mat){
	//mat should be 4by 4 lfirt matrix
	//cout<<"done flirt"<<endl;
	
	ColumnVector pdef(4);
	pdef.element(4)=1;
	
	float resx=defField[0].xdim();
	float resy=defField[0].ydim();
	float resz=defField[0].zdim();
	
	for (vector<mesh::Mpoint*>::iterator i=m._points.begin();i!=m._points.end();i++){
		//shifts are for changes between roi and full image
		
		
		float px=(*i)->get_coord().X;
		float py=(*i)->get_coord().Y;
		float pz=(*i)->get_coord().Z;
		pdef.element(0)=px+defField[0].interpolate(px/resx,py/resy,pz/resz);
		pdef.element(1)=py+defField[1].interpolate(px/resx,py/resy,pz/resz) ;
		pdef.element(2)=pz+defField[2].interpolate(px/resx,py/resy,pz/resz) ;
		pdef=mat*pdef;
		
		(*i)->_update_coord.X=pdef.element(0);
		(*i)->_update_coord.Y=pdef.element(1);
		(*i)->_update_coord.Z=pdef.element(2);
	}
}*/



template<class Tdist,class Tim>
void meshUtils::SurfDistToLabels( vector<Tdist> & vdist, const volume<Tim> & image){
cout<<"enter surf dist "<<endl;

		vdist.clear();
cout<<"enter surf dist "<<endl;

	int bounds[6]={0,0,0,0,0,0};
	float xdim=image.xdim();
	float ydim=image.ydim();
	float zdim=image.zdim();
	getBounds(bounds,xdim,ydim,zdim);
	//switch to using vectors later
	for (int i=0; i<Points.Nrows();i++)
	{
		cout<<"point "<<i<<endl;
		float px=Points.element(i,0);
		float py=Points.element(i,1);
		float pz=Points.element(i,2);
		int centre[3]={static_cast<short>(px), static_cast<short>(py),static_cast<short>(pz)};
		Tdist dist=10000000;
		int layer=0;
		bool notFoundD=false;
		while (notFoundD)
		{ 
			for (int x=centre[0]-layer; x<=centre[0]+layer; x++)
				for (int y=centre[1]-layer; y<=centre[1]+layer;y++)
					for (int z=centre[2]-layer; z<=centre[2]+layer;z++)
					{	
						if (image.value(x,y,z)>0)
						{
							float temp=x*x*xdim*xdim+y*y*ydim*ydim+z*z*zdim*zdim;
							if (temp<dist)
							{
								dist=temp;
								notFoundD=false;
							}
						}
					}	
						layer++;
		}
		vdist.push_back(static_cast<Tdist>(sqrt(dist)));	
	}
	
	//setScalars(vdist);
}
template void meshUtils::SurfDistToLabels<float,unsigned int>( vector<float> & dist,const volume<unsigned int>& image);
template void meshUtils::SurfDistToLabels<float,short>( vector<float> & dist,const volume<short>& image);
template void meshUtils::SurfDistToLabels<float,int>( vector<float> & dist,const volume<int>& image);
template void meshUtils::SurfDistToLabels<float,float>( vector<float> & dist,const volume<float>& image);
template void meshUtils::SurfDistToLabels<float,double>( vector<float> & dist,const volume<double>& image);

template<class Tdist,class Tim>
void meshUtils::SurfDistToLabels( vector<Tdist> & vdist, const volume<Tim> & image, const Tim & label){
		vdist.clear();


	float xdim=image.xdim();
	float ydim=image.ydim();
	float zdim=image.zdim();
	for (int i=0; i<Points.Nrows();i++)
	{
			float dist=10000000;
		cout<<"point "<<i<<endl;

		float px=Points.element(i,0);
		float py=Points.element(i,1);
		float pz=Points.element(i,2);
		for (int i=0;i<image.xsize();i++)
			for (int j=0;j<image.ysize();j++)
				for (int k=0;k<image.zsize();k++)
						{
					//		cout<<image.value(i,j,k)<<" "<<label<<endl;
							if (image.value(i,j,k)==label)
							{
							float temp=(i*xdim-px)*(i*xdim-px) + (j*ydim-py)*(j*ydim-py)+ (k*zdim-pz)*(k*zdim-pz);
							//cout<<"temp "<<temp<<endl;
							if (temp<dist)
								dist=temp;
							}
						}
		vdist.push_back(static_cast<Tdist>(sqrt(dist)));
		
	}
	
	//setScalars(vdist);
}
/*
template<class Tdist,class Tim>
void meshUtils::SurfDistToLabels( vector<Tdist> & vdist, const volume<Tim> & image, const Tim & label){
		vdist.clear();


	int bounds[6]={0,0,0,0,0,0};
	float xdim=image.xdim();
	float ydim=image.ydim();
	float zdim=image.zdim();
	getBounds(bounds,xdim,ydim,zdim);
	//switch to using vectors later
	for (int i=0; i<Points.Nrows();i++)
	{
		float px=Points.element(i,0);
		float py=Points.element(i,1);
		float pz=Points.element(i,2);
		int centre[3]={static_cast<short>(px/xdim), static_cast<short>(py/ydim),static_cast<short>(pz/zdim)};
				cout<<"centre "<<centre[0]<<" "<<centre[1]<<" "<<centre[2]<<endl;

		Tdist dist=10000000;
		int layer=0;
		bool notFoundD=true;
		cout<<"test"<<endl;
		while (notFoundD)
		{ 
			cout<<"layer "<<layer<<endl;
			for (int x=centre[0]-layer; x<=centre[0]+layer; x++)
				for (int y=centre[1]-layer; y<=centre[1]+layer;y++)
					for (int z=centre[2]-layer; z<=centre[2]+layer;z++)
					{	
					//	cout<<x<<" "<<y<<" "<<z<<" "<<label<<" "<<image.value(x,y,z)<<endl;
						if (image.value(x,y,z)==label)
						{
							float xdif=px-x*xdim;
							float ydif=py-y*ydim;
							float zdif=pz-z*zdim;
							float temp=xdif*xdif+ydif*ydif+zdif*zdif;
							cout<<"temp "<<temp<<endl;
							if (temp<dist)
							{
								dist=temp;
								notFoundD=false;
							}
						}
					}	
						layer++;
		}
		vdist.push_back(static_cast<Tdist>(sqrt(dist)));	
	}
	
	//setScalars(vdist);
}
*/
template void meshUtils::SurfDistToLabels<float,unsigned int>( vector<float> & dist,const volume<unsigned int>& image, const unsigned & label);
template void meshUtils::SurfDistToLabels<float,short>( vector<float> & dist,const volume<short>& image, const short & label);
template void meshUtils::SurfDistToLabels<float,int>( vector<float> & dist,const volume<int>& image, const int & label);
template void meshUtils::SurfDistToLabels<float,float>( vector<float> & dist,const volume<float>& image, const float & label);
template void meshUtils::SurfDistToLabels<float,double>( vector<float> & dist,const volume<double>& image, const double & label);




void meshUtils::SurfScalarsMeanAndStdev(vector<string> meshList, Matrix & MeanPoints, Matrix & MeanScalars, Matrix & StDevScalars )
{
	Matrix tempScalars;
	Matrix tempScalarsSq;

	for (unsigned int i=0;i<meshList.size();i++)
	{
		loadMesh(meshList.at(i));
		if (i==0)
		{
			MeanPoints=Points;
			tempScalars=Scalars;
			tempScalarsSq=SP(Scalars,Scalars);
		}else
		{
			MeanPoints+=Points;
			tempScalars+=Scalars;
			tempScalarsSq+=SP(Scalars,Scalars);
		}
	}
	float N=static_cast<float>(meshList.size());
	MeanPoints/=N;
	MeanScalars=tempScalars/N;
	StDevScalars=(SP(tempScalars,tempScalars)-tempScalarsSq/N)/(N-1);
	for (int i=0; i<StDevScalars.Nrows();i++)
		StDevScalars.element(i,0)=sqrt(StDevScalars.element(i,0));
	
}


template<class T>
void meshUtils::warpMeshWithDefField(mesh::Mesh & m, const volume4D<T> & defField, const Matrix & mat){
			//mat should be 4by 4 lfirt matrix
	//cout<<"done flirt"<<endl;
	
	ColumnVector pdef(4);
	pdef.element(3)=1;
	
	float resx=defField[0].xdim();
	float resy=defField[0].ydim();
	float resz=defField[0].zdim();
	
	for (vector<mesh::Mpoint*>::iterator i=m._points.begin();i!=m._points.end();i++){
		//shifts are for changes between roi and full image
		
		
		float px=(*i)->get_coord().X;
		float py=(*i)->get_coord().Y;
		float pz=(*i)->get_coord().Z;
		pdef.element(0)=px+defField[0].interpolate(px/resx,py/resy,pz/resz);
		pdef.element(1)=py+defField[1].interpolate(px/resx,py/resy,pz/resz) ;
		pdef.element(2)=pz+defField[2].interpolate(px/resx,py/resy,pz/resz) ;
		pdef=mat*pdef;
		
		(*i)->_update_coord.X=pdef.element(0);
		(*i)->_update_coord.Y=pdef.element(1);
		(*i)->_update_coord.Z=pdef.element(2);
	}
	m.update();
}
template void meshUtils::warpMeshWithDefField<float>( Mesh & m, const volume4D<float> & defField, const Matrix & mat);
template void meshUtils::warpMeshWithDefField<int>( Mesh & m, const volume4D<int> & defField, const Matrix & mat);
template void meshUtils::warpMeshWithDefField<double>( Mesh & m, const volume4D<double> & defField, const Matrix & mat);
template void meshUtils::warpMeshWithDefField<short>( Mesh & m, const volume4D<short> & defField, const Matrix & mat);


//template void meshUtils::warpMeshWithDefField<float>( Mesh & m, const volume4D<float> & defField, const Matrix & mat);




template<class T>
void meshUtils::ugridToImage(volume<T> & im)
{
	//for now use nearest neighbour
	float xres=im.xdim();
	float yres=im.ydim();
	float zres=im.zdim();
	cout<<"scalars size "<<Scalars.Nrows()<<" "<<xres<<" "<<yres<<" "<<zres<<endl;
	//float dist_th=sqrt(xres+xres+yres*yres+zres*zres);
	im=0;
	for (int i=0;i<Points.Nrows();i++)
	{
       	//cout<<"i "<<i<<endl;
	  im.value(static_cast<int>(Points.element(i,0)/xres),static_cast<int>(Points.element(i,1)/yres),static_cast<int>(Points.element(i,2)/zres))=static_cast<T>(Scalars.element(i,0));
		//cout<<Scalars.element(i,0)<<" "<<Scalars.Nrows()<<endl;
	}
	
	
	/*
	for (int i=0;i<im.xsize();i++)
	{
		cout<<"iter "<<i<<" of "<<im.xsize()<<endl;
		for (int j=0;j<im.ysize();j++)
		{
			for (int k=0;k<im.zsize();k++)
			{
				//search for closest point 
				float mindist=(Points.element(0,0)-i*xres)*(Points.element(0,0)-i*xres) \
				+(Points.element(0,1)-j*yres)*(Points.element(0,1)-j*yres) \
				+(Points.element(0,2)-k*zres)*(Points.element(0,2)-k*zres) ;
				int indmin=0;
				for (int x=1;x<Points.Nrows();x++)
				{
					float dx=(Points.element(x,0)-i*xres);
					float dy=(Points.element(x,1)-j*yres);
					float dz=(Points.element(x,2)-k*zres);
					
					if ((dx*dx+dy*dy+dz*dz) < mindist){
						//				cout<<i<<" "<<j<<" "<<k<<" search "<<indmin<<" "<<Points.element(x,0)<<" "<<Points.element(x,1)<<" "<<Points.element(x,2)<<" "<<dx<<" "<<dy<<" "<<dz<<endl;
						indmin=x;
						mindist=(dx*dx+dy*dy+dz*dz); 
						if (mindist<dist_th) break;
					}
				}	
				cout<<"indmind "<<indmin<<" "<<mindist<<" "<<dist_th<<endl;
				im.value(i,j,k)=static_cast<T>(Scalars.element(indmin,0));
		
			}
			//	cout<<i<<" "<<j<<" "<<k<<" INDMIN "<<indmin<<" "<<Points.element(indmin,0)<<" "<<Points.element(indmin,1)<<" "<<Points.element(indmin,2)<<endl;
				
		}
	}
	*/
	
}	
template void meshUtils::ugridToImage<short>(volume<short> & im);
template void meshUtils::ugridToImage<float>(volume<float> & im);
template void meshUtils::ugridToImage<double>(volume<double> & im);


ReturnMatrix meshUtils::getDeformedVector(const ColumnVector & mean, const Matrix & modes, const ColumnVector & eigs, const vector<float> & vars )
{
	ColumnVector def(mean.Nrows());
	def=mean;
	for (unsigned int i=0;i<vars.size();i++)
	{
		def+=modes.Column(i+1)*sqrt(eigs.element(i))*vars.at(i);
		//cout<<"var1 "<<vars.at(i)<<" "<<sqrt(eigs.element(i))<<endl;
	}	
		
	def.Release();
	return def;

}
template<class T,class T2>
	void meshUtils::deformSurface(const volume<T> & image, const float & maxit, const float & wIm, const float & wTang_in, const float & w_area, const float & w_norm, const T & max_thresh, const unsigned int & interRate, const bool & enableInteraction, const string & name)
	{
		const T  max_thresh_sq=max_thresh*max_thresh;
		const T xdim = static_cast<T>(image.xdim()); 
		const T ydim = static_cast<T>(image.ydim()); 
		const T zdim = static_cast<T>(image.zdim()); 
		float wTang=wTang_in;
		
#ifdef USE_VTK
		
		vtkPolyData *meshData = vtkPolyData::New();
		vtkPolyDataReader *vtkmesh = vtkPolyDataReader::New();
		vtkRenderWindow *renWin = vtkRenderWindow::New();
		vtkRenderer *ren1 =vtkRenderer::New();
		vtkPolyDataMapper *meshmap = vtkPolyDataMapper::New();
		vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
		vtkImageData *imStruct =vtkImageData::New();
		vtkImagePlaneWidget *planeWidgetX = vtkImagePlaneWidget::New();
		
		vtkmesh->SetFileName(name.c_str());
		meshData=(vtkmesh->GetOutput());
		meshmap->SetInput(meshData);
		
		//actor coordinates, geome..etc
		vtkActor *aSurf = vtkActor::New();
		aSurf->GetProperty()->SetColor(0.7,0,0);
		aSurf->GetProperty()->SetOpacity(1); 
		aSurf->GetProperty()->SetEdgeColor(0.0,0,8);
		aSurf->GetProperty()->SetEdgeVisibility(1);
		aSurf->GetProperty()->SetPointSize(2);
		aSurf->GetProperty()->SetRepresentationToWireframe();
		aSurf->SetMapper(meshmap);
		
		//interactorv
		iren->SetRenderWindow(renWin);
		
		//add actor
		ren1->AddActor(aSurf);
		//ren1->AddActor(aSphere2);
		ren1->SetBackground(0.4,0.4,0.4);
		
		//Render
		renWin->AddRenderer(ren1);
		
		
		imStruct=newimageTovtkImageData(image);	
		
		
		planeWidgetX->DisplayTextOn();
		planeWidgetX->SetInput(imStruct);
		planeWidgetX->SetPlaneOrientationToYAxes();
		planeWidgetX->SetSliceIndex(120);
		planeWidgetX->SetKeyPressActivationValue('x');
		planeWidgetX->SetInteractor(iren);
		planeWidgetX->TextureVisibilityOn();
		planeWidgetX->GetTexturePlaneProperty()->SetOpacity(0.99);
		planeWidgetX->SetEnabled(1);
		renWin->Render();
		
		if (enableInteraction)
			iren->Start();
		
		
#endif
		
		vector<T> pts_org=getPointsAsVector<T>();
		vector<T> pts=pts_org;
		
		
		vector<  vector<T2> > cells=this->getPolygonsAsVectorOfVectors<T2>();
		vector<  vector<T2> > neighbours=first_mesh::findNeighbours<T2>(cells, pts.size()/3);
		vector<  vector<T2> > neighbour_tris=first_mesh::findNeighbourTriangles<T2>(cells, pts.size()/3);
		
		
		for (int n_iter=1; n_iter<=maxit; n_iter++) 
		{  //iterations
			
			//calculate normals
			vector<T> vnx,vny,vnz;
			first_mesh::normal<T,T2>(pts,neighbour_tris,cells,vnx,vny,vnz);
			
			vector<T> vtri_x,vtri_y,vtri_z;
			first_mesh::maxTriangle<T,T2>(pts,neighbour_tris,cells,vtri_x,vtri_y,vtri_z);
			
			//clauclate medium of neigbours than teh diference
			vector<T> vmx,vmy,vmz;
			first_mesh::medium_neighbours<T,T2>(pts,neighbours,cells,vmx,vmy,vmz);
			
			typename vector<T>::iterator j=vmy.begin();
			typename vector<T>::iterator k=vmz.begin();
			typename vector<T>::iterator l=pts.begin();
			
			typename vector<T>::iterator m=vnx.begin();
			typename vector<T>::iterator n=vny.begin();		
			typename vector<T>::iterator o=vnz.begin();
			//find diffenrece vector
			
			vector<T> v_im_nx,v_im_ny,v_im_nz;
			for (typename vector<T>::iterator i=vmx.begin();i!=vmx.end();i++,j++,k++,m++,n++,o++,l+=3)
			{
				//calculate difference vector
				*i-=*l;
				*j-=*(l+1);
				*k-=*(l+2);
				
				T im_int = image.interpolate( (*l)/xdim, (*(l+1))/ydim,  (*(l+2))/zdim );
				v_im_nx.push_back(im_int * (*m) );
				v_im_ny.push_back(im_int * (*n) );
				v_im_nz.push_back(im_int * (*o) );	
				
				T norm= (*i) * (*m)+ (*j) * (*n) + (*k) * (*o);
				
				(*m) *= norm;
				(*n) *= norm;
				(*o) *= norm;
				
				(*i)-=(*m);
				(*j)-=(*n);
				(*k)-=(*o);
				
			}
			
			
			vector<T> v_dx, v_dy, v_dz;
			typename vector<T>::iterator tang_xi = vmx.begin();
			typename vector<T>::iterator tang_yi = vmy.begin();
			typename vector<T>::iterator tang_zi = vmz.begin();
			typename vector<T>::iterator norm_xi = vnx.begin();
			typename vector<T>::iterator norm_yi = vny.begin();
			typename vector<T>::iterator norm_zi = vnz.begin();
			typename vector<T>::iterator norm_im_xi = v_im_nx.begin();
			typename vector<T>::iterator norm_im_yi = v_im_ny.begin();
			typename vector<T>::iterator norm_im_zi = v_im_nz.begin();
			typename vector<T>::iterator vtri_xi = vtri_x.begin();
			typename vector<T>::iterator vtri_yi = vtri_y.begin();
			typename vector<T>::iterator vtri_zi = vtri_z.begin();
			
			float max_step=0;
			for (typename vector<T>::iterator i=pts.begin();i!=pts.end();i+=3, tang_xi++,tang_yi++,tang_zi++,norm_xi++,norm_yi++,norm_zi++, norm_im_xi++ , norm_im_yi++ , norm_im_zi++, vtri_xi++, vtri_yi++, vtri_zi++ )
			{
				
				float dx = wTang * (*tang_xi) + w_norm * (*norm_xi) + wIm * (*norm_im_xi) + w_area * (*vtri_xi);
				float dy = wTang * (*tang_yi) + w_norm * (*norm_yi) + wIm * (*norm_im_yi) + w_area * (*vtri_yi);
				float dz = wTang * (*tang_zi) + w_norm * (*norm_zi) + wIm * (*norm_im_zi) + w_area * (*vtri_zi);
				v_dx.push_back( dx );
				v_dy.push_back( dy );
				v_dz.push_back( dz );
				
				if ( (dx*dx +dy*dy +dz*dz) > max_step )
					max_step=(dx*dx +dy*dy +dz*dz);
			}
			
			float scale=1;
			if ( max_step > max_thresh_sq)
				scale=max_thresh/sqrt(max_step);
			
			//update points
			typename vector<T>::iterator dx_i = v_dx.begin();
			typename vector<T>::iterator dy_i = v_dy.begin();
			typename vector<T>::iterator dz_i = v_dz.begin();
			int count=0;
			for (typename vector<T>::iterator i=pts.begin();i!=pts.end();i+=3,dx_i++,dy_i++,dz_i++,count++)
			{
				*i+=(*dx_i)*scale;
				*(i+1)+=(*dy_i)*scale;
				*(i+2)+=(*dz_i)*scale;
			}
			
			if ( n_iter % interRate == 0)
				if (first_mesh::self_intersection_test<T,T2>(cells,pts))
				{
					pts=pts_org;
					wTang+=0.1;
					cout<<"THE MESH SELF-INTERSECTS "<<n_iter<<" "<<wTang<<endl;
					n_iter=0;
				}
			
#ifdef USE_VTK
			fslvtkconvname::update_vtkPolyData(pts,meshData);
			meshData->GetPoints()->Modified();
			renWin->Render();
#endif
		}
#ifdef USE_VTK

			
		if (enableInteraction)
			iren->Start();
#endif
		
		
		setPoints(pts);
	}

template void meshUtils::deformSurface<float,short>(const volume<float> & image, const float & maxit, const float & wIm, const float & wTang, const float & maxTri, const float & w_norm, const float & max_thresh, const unsigned int & interRate, const bool & enableInteraction,const string & name);
template void meshUtils::deformSurface<float,unsigned int>(const volume<float> & image, const float & maxit, const float & wIm, const float & wTang, const float & maxTri, const float & w_norm, const float & max_thresh, const unsigned int & interRate,const bool & enableInteraction, const string & name);



void meshUtils::applyReg(Matrix & pts, const Matrix & fmat)
{
	ColumnVector ones(pts.Nrows());
	ones=1;
	Matrix temp= (pts | ones);
//	cout<<"dikms "<<Pts_src_dm.Nrows()<<" "<<Pts_src_dm.Ncols()<<" "<<ones.Nrows()<<" "<<ones.Ncols()<<endl;//temp.Ncols()<<endl;
	temp=(fmat*(pts | ones).t());
	temp=temp.t();
	pts=temp.SubMatrix(1,temp.Nrows(),1,3);

}

ReturnMatrix meshUtils::calculateRotation(const Matrix & Pts_src_dm, const Matrix & Pts_targ_dm)
{	
	Matrix M_Rot(4,4);
	M_Rot=0;
	M_Rot.element(3,3)=1;


		Matrix M=Pts_targ_dm.t()*(Pts_src_dm);
		Matrix U;
		DiagonalMatrix D;
		SVD(M.t()*M,D,U);
		//M should always be a 3x3 matrix
		for (int i=0;i<D.Nrows();i++){
			D.element(i)=1/sqrt(D.element(i));
		}
		
		Matrix R(3,3);
		R=M*(U*D*U.t());
		
		
		SVD(R,D,U);
//cout<<"D "<<D<<endl;
		
		M_Rot.SubMatrix(1,3,1,3)=R;
	
		M_Rot.Release();
		return M_Rot;
}
ReturnMatrix meshUtils::calculateScale(const Matrix & Pts_src_dm, const Matrix & Pts_targ_dm, const bool & global)
{	
	Matrix M_Sc(4,4);
	M_Sc=0;
	M_Sc.element(0,0)=1;M_Sc.element(1,1)=1;M_Sc.element(2,2)=1;M_Sc.element(3,3)=1;
		float sumrx=0, sumry=0,sumrz=0, sumlx=0,sumly=0, sumlz=0;
			for (unsigned int i=0;i<static_cast<unsigned int>(Pts_src_dm.Nrows()) ;i++){
				sumlx += Pts_src_dm.element(i,0)*Pts_src_dm.element(i,0);
				sumly += Pts_src_dm.element(i,1)*Pts_src_dm.element(i,1);
			    sumlz += Pts_src_dm.element(i,2)*Pts_src_dm.element(i,2);
			  
				sumrx += Pts_targ_dm.element(i,0)*Pts_targ_dm.element(i,0);
				sumry += Pts_targ_dm.element(i,1)*Pts_targ_dm.element(i,1);
				sumrz += Pts_targ_dm.element(i,2)*Pts_targ_dm.element(i,2);
			
			}
			if (global)
			{
				float scale=sqrt((sumrx+sumry+sumrz)/(sumlx+sumly+sumlz));
				M_Sc.element(0,0)=scale;M_Sc.element(1,1)=scale;M_Sc.element(2,2)=scale;
			}else {
								M_Sc.element(0,0)=sqrt(sumrx/sumlx);M_Sc.element(1,1)=sqrt(sumry/sumly);M_Sc.element(2,2)=sqrt(sumrz/sumlz);
			}
			M_Sc.Release();
			return M_Sc;
}

Matrix meshUtils::reg_leastsq(const Matrix & SourcePoints, const short & dof)
{
	
	unsigned int Npoints=Points.Nrows(); 
	
		Matrix M_trans_src(4,4);
	M_trans_src=0;
	M_trans_src.element(0,0)=1;	M_trans_src.element(1,1)=1;	M_trans_src.element(2,2)=1;M_trans_src.element(3,3)=1;

	
	Matrix M_trans_targ(4,4);
	M_trans_targ=0;
	M_trans_targ.element(0,0)=1;	M_trans_targ.element(1,1)=1;	M_trans_targ.element(2,2)=1;M_trans_targ.element(3,3)=1;
	
	
	Matrix M_Sc(4,4);
	M_Sc=0;
	M_Sc.element(0,0)=1;M_Sc.element(1,1)=1;M_Sc.element(2,2)=1;M_Sc.element(3,3)=1;
	
	Matrix M_Rot(4,4);
	M_Rot=0;
	M_Rot.element(0,0)=1;M_Rot.element(1,1)=1;M_Rot.element(2,2)=1;M_Rot.element(3,3)=1;
	//---------------------CALCULATE CENTROIDS----------------------//
	double mean_src_x=0, mean_src_y=0, mean_src_z=0;
	double mean_targ_x=0, mean_targ_y=0,mean_targ_z=0;
	
	for (unsigned int i=0; i<Npoints;i++)
	{
		mean_src_x+=SourcePoints.element(i,0);
		mean_src_y+=SourcePoints.element(i,1);
		mean_src_z+=SourcePoints.element(i,2);

		mean_targ_x+=Points.element(i,0);
		mean_targ_y+=Points.element(i,1);
		mean_targ_z+=Points.element(i,2);
	} 
		
		mean_src_x/=Npoints;
		mean_src_y/=Npoints;
		mean_src_z/=Npoints;

		mean_targ_x/=Npoints;
		mean_targ_y/=Npoints;
		mean_targ_z/=Npoints;

	//---------------------TRANSLATION COMPONENTS----------------------//

	M_trans_src.element(0,3)=-mean_src_x;	M_trans_src.element(1,3)=-mean_src_y;	M_trans_src.element(2,3)=-mean_src_z;
	
	M_trans_targ.element(0,3)=mean_targ_x;	M_trans_targ.element(1,3)=mean_targ_y;	M_trans_targ.element(2,3)=mean_targ_z;
	
	//---------------------DEMEAN POINTS ---------------------------//

	Matrix Pts_src_dm=SourcePoints;
	Matrix Pts_targ_dm=Points;
	
	for (unsigned int i=0;i<Npoints;i++)
	{
		 Pts_src_dm.element(i,0)-=mean_src_x;
		 Pts_src_dm.element(i,1)-=mean_src_y;
		 Pts_src_dm.element(i,2)-=mean_src_z;
		 
		 Pts_targ_dm.element(i,0)-=mean_targ_x;
		 Pts_targ_dm.element(i,1)-=mean_targ_y;
		 Pts_targ_dm.element(i,2)-=mean_targ_z;
	}

	
	
		//---------------------ROTATION COMPONENT (OPTIONAL)----------------------//
Matrix M_Sc_final=M_Sc;
Matrix M_Rot_final=M_Rot;

for (unsigned int i=0;i<100;i++)
{ 
	if (dof>=6)
	{
		M_Rot=calculateRotation(Pts_src_dm,Pts_targ_dm);
		M_Rot_final=M_Rot*M_Rot_final;
		if (dof==6)
			break;
	}//----------------------ADD NINE DEGREES OF FREEDON---------------------//
		//---------------------SCALE COMPONENT (OPTIONAL)----------------------//
//	cout<<"do Scale"<<endl;
	if (dof==7)
	{
		M_Sc_final=calculateScale(Pts_src_dm,Pts_targ_dm,true);
		break;
	}else if (dof==9)
	{
		applyReg(Pts_src_dm,M_Rot);
		M_Sc=calculateScale(Pts_src_dm,Pts_targ_dm,false);
		M_Sc_final=M_Sc*M_Sc_final;
		applyReg(Pts_src_dm,M_Sc);

	}
//	cout<<"rot "<<M_Rot<<endl;
//	cout<<"sc "<<M_Sc<<endl;
	}
///	cout<<"trans src "<<M_trans_src<<endl;
//		cout<<"trans targ "<<M_trans_targ<<endl;

	cout<<"rot "<<M_Rot_final<<endl;
	cout<<"sc "<<M_Sc_final<<endl;
//	cout<<M_trans_targ  *  M_Sc_final * M_Rot_final * M_trans_src<<endl;
	
	return M_trans_targ * M_Sc_final * M_Rot_final * M_trans_src;
	
}

Matrix meshUtils::alignSurfaces(const string & src_list, const short & dof, const string & outname )
{

	Matrix allPoints(1,1);
	ifstream innames;
	innames.open(src_list.c_str());
	string fname;

	fslvtkIO* fout = new fslvtkIO();
	while (innames>>fname)
	{
		meshUtils* meshIn=new meshUtils(fname,static_cast<meshUtils::DataType>(0));
		//Matrix Pts_src=;
		Matrix fmat=this->reg_leastsq(meshIn->getPointsAsMatrix(),dof);
		//cout<<"points orh "<<meshIn->getPointsAsMatrix()<<endl;
				meshIn->meshReg(fmat);
		if ((allPoints.Nrows()==1)&&(allPoints.Ncols()==1))
		{
			if (strcmp(outname.c_str(),"nosave"))
			{
				cout<<"adding DDDDDDDpoints"<<endl;
			fout->setPolygons(meshIn->getPolygons());
			fout->setPoints(meshIn->getPointsAsMatrix());						
			}
			allPoints = first_newmat_vector::unwrapMatrix(meshIn->getPointsAsMatrix());
		}else 
		{
			cout<<fout->getPointsAsMatrix().Nrows()<<" "<<fout->getPointsAsMatrix().Ncols()<<" "<<meshIn->getPointsAsMatrix().Nrows()<<" "<<meshIn->getPointsAsMatrix().Ncols()<<" "<<endl;
			if (strcmp(outname.c_str(),"nosave"))
					fout->appendPointsAndPolygons(meshIn->getPointsAsMatrix(), meshIn->getPolygons());						
			allPoints = allPoints |   first_newmat_vector::unwrapMatrix(meshIn->getPointsAsMatrix());
		}
			//	cout<<"points reg "<<meshIn->getPointsAsMatrix()<<endl;

		//meshIn->save(outname.value());
		cout<<"allPoints "<<allPoints.Nrows()<<" "<<allPoints.Ncols()<<endl;
		delete meshIn;
	}
	innames.close();

			if (strcmp(outname.c_str(),"nosave"))
			fout->save(outname);
	
	delete fout;
	
	return allPoints;
}



}	
