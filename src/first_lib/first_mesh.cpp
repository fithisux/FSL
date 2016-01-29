/*
 *  first_mesh.h
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

#include "first_mesh.h"
using namespace std;
using namespace NEWIMAGE;

namespace FIRST_LIB{

	first_mesh::first_mesh(){
	
	}

	first_mesh::~first_mesh(){
	
	}


void first_mesh::getBounds(const vector<float> & m, int *bounds, const float & xdim, const float & ydim, const float & zdim)
{	
	float xmin=1000,xmax=-1000,ymin=1000,ymax=-1000,zmin=1000,zmax=-1000;
	for (vector<float>::const_iterator i=m.begin();i!=m.end();i+=3)
	{
		if ( *i < xmin)		xmin=*i;
		if ( *i > xmax)		xmax=*i;
		if ( *(i+1) < ymin)	ymin=*(i+1);
		if ( *(i+1) > ymax)	ymax=*(i+1);
		if ( *(i+2) < zmin)	zmin=*(i+2);
		if ( *(i+2) > zmax)	zmax=*(i+2);
	}
	
	*bounds=static_cast<int>(floor(xmin/xdim)-1);
	*(bounds+1)=static_cast<int>(ceil(xmax/xdim)+1);
	*(bounds+2)=static_cast<int>(floor(ymin/ydim)-1);
	*(bounds+3)=static_cast<int>(ceil(ymax/ydim)+1);
	*(bounds+4)=static_cast<int>(floor(zmin/zdim)-1);
	*(bounds+5)=static_cast<int>(ceil(zmax/zdim)+1);	
}


void  first_mesh::draw_segment(::volume<short>& image, const float & p1x,const float & p1y,const float & p1z,  const float & p2x,const float & p2y,const float & p2z, int label)
{
	double xdim = (double) image.xdim();
	double ydim = (double) image.ydim();
	double zdim = (double) image.zdim();

	//in new version of bet2
	double mininc = min(xdim,min(ydim,zdim)) * .25;

	float nx = p1x - p2x;
	float ny = p1y - p2y;
	float nz = p1z - p2z;

	double d = sqrt(nx*nx+ny*ny+nz*nz);
	nx/=d;
	ny/=d;
	nz/=d;
	
	for (double i=0; i<=d; i+=mininc)
		image(static_cast<int>(floor((p2x+i*nx)/xdim +.5)),static_cast<int>(floor((p2y+i*ny)/ydim +.5)),static_cast<int>(floor((p2z+i*nz)/zdim +.5))) = label;
}


volume<short>  first_mesh::draw_mesh(const volume<short> & image, const vector<float> & pts, const vector< vector<unsigned int> > & triangles, int label)
{
	volume<short> imout;
	copyconvert(image,imout);
	imout=0;
  	double xdim = (double) image.xdim();
	double ydim = (double) image.ydim();
	double zdim = (double) image.zdim();
	
	//in new version of bet2
	double mininc = min(xdim,min(ydim,zdim)) * 0.5;
	
	for (vector< vector<unsigned int> >::const_iterator i = triangles.begin(); i!=triangles.end(); i++)
	{	
		float px1= pts.at(i->at(1)*3);
		float py1= pts.at(i->at(1)*3+1);
		float pz1= pts.at(i->at(1)*3+2);
						  
		float vx=pts.at(i->at(0)*3) - px1;
		float vy=pts.at(i->at(0)*3+1) - py1;
		float vz=pts.at(i->at(0)*3+2) - pz1;
						  
		float px2=pts.at(i->at(2)*3); 
		float py2=pts.at(i->at(2)*3+1); 
		float pz2=pts.at(i->at(2)*3+2);
						  
		float d=sqrt(vx*vx+vy*vy+vz*vz);
						  
		vx/=d;
		vy/=d;
		vz/=d;
						  
		for (float j=0; j<=d ;  j+=mininc)
			draw_segment(imout, px1+j*vx,py1+j*vy,pz1+j*vz,px2,py2,pz2,label+100);
						  
    }
	return imout;
}


volume<short> first_mesh::make_mask_from_mesh(const volume<float> & image, const vector<float> & m, const vector< vector<unsigned int> > & triangles, const int & label, const int * bounds)
{
	volume<short> mask;
	copyconvert(image,mask);
	
	int boundsNew[6];// = *bounds;
	for (int i=0;i<6;i++)
		boundsNew[i]=bounds[i];
	mask = 0;
	mask = draw_mesh(mask, m, triangles, label);
	
	volume<short> otl=mask;
	
	vector<float> currentX,currentY,currentZ;
	if (boundsNew[0]<0) boundsNew[0]=0;
	if (boundsNew[2]<0) boundsNew[2]=0;
	if (boundsNew[4]<0) boundsNew[4]=0;
	
	if (boundsNew[1]>=image.xsize()-1) boundsNew[1]=image.xsize()-1; 
	if (boundsNew[3]>=image.ysize()-1) boundsNew[3]=image.ysize()-1;
	if (boundsNew[5]>=image.zsize()-1) boundsNew[5]=image.zsize()-1; 
	
	float start_x=boundsNew[0]-2;
	float start_y=boundsNew[2]-2;
	float start_z=boundsNew[4]-2;
	
	if (boundsNew[0]<2) start_x=0;
	if (boundsNew[2]<2) start_y=0;
	if (boundsNew[4]<2) start_z=0;
	
	
	mask.value(start_x,start_y,start_z) = label;
	
	currentX.push_back(start_x);
	currentY.push_back(start_y);
	currentZ.push_back(start_z);
	
	while (!currentX.empty())
	{
		int x, y, z;
		x=static_cast<int>(currentX.back()); 
		y=static_cast<int>(currentY.back()); 
		z=static_cast<int>(currentZ.back()); 

		currentX.pop_back();
		currentY.pop_back();
		currentZ.pop_back();

		if (boundsNew[0]<=x-1 && mask.value(x-1, y, z)==0) 
		{
			mask.value(x-1, y, z) = label;
			currentX.push_back(x-1);
			currentY.push_back(y);
			currentZ.push_back(z);

		}
		if (boundsNew[2]<=y-1 && mask.value(x, y-1, z)==0) 
		{
			mask.value(x, y-1, z) = label;
			currentX.push_back(x);
			currentY.push_back(y-1);
			currentZ.push_back(z);

		}
		if (boundsNew[4]<=z-1 && mask.value(x, y, z-1)==0) 
		{
			mask.value(x, y, z-1) = label;
			currentX.push_back(x);
			currentY.push_back(y);
			currentZ.push_back(z-1);

		}
		if (boundsNew[1]>=x+1 && mask.value(x+1, y, z)==0)
		{
			mask.value(x+1, y, z) = label;
			currentX.push_back(x+1);
			currentY.push_back(y);
			currentZ.push_back(z);

		}
		if (boundsNew[3]>=y+1 && mask.value(x, y+1, z)==0)
		{
			mask.value(x, y+1, z) = label;
			currentX.push_back(x);
			currentY.push_back(y+1);
			currentZ.push_back(z);

		}
		if (boundsNew[5]>=z+1 && mask.value(x, y, z+1)==0)
		{
			mask.value(x, y, z+1) = label;
			currentX.push_back(x); 
			currentY.push_back(y); 
			currentZ.push_back(z+1); 

		}
		
	}

	for (int i=boundsNew[0];i<boundsNew[1];i++)
		for (int j=boundsNew[2];j<boundsNew[3];j++)
			for (int k=boundsNew[4];k<boundsNew[5];k++)
				if (mask.value(i,j,k)==0)
					otl.value(i,j,k)=label;

	return otl;
}

template<class T,class T2>
void first_mesh::normal(const vector<T> & pts, const vector< vector<T2> > & localTri, const vector< vector<T2> > & cells,vector<T> & vx, vector<T> & vy, vector<T> & vz )
{
	
	vx.clear();
	vy.clear();
	vz.clear();
	//calculate the normals for each triangle then references
		
	for (typename vector< vector<T2> >::const_iterator i=localTri.begin();i!=localTri.end();i++)
	{
		float nx=0,ny=0,nz=0;
		float vx1,vy1,vz1,vx2,vy2,vz2;
	
		for (typename vector<T2>::const_iterator j=i->begin();j!=i->end();j++)
		{
			vx1=pts.at(cells.at(*j).at(2)*3)-pts.at(cells.at(*j).at(0)*3);
			vx2=pts.at(cells.at(*j).at(1)*3)-pts.at(cells.at(*j).at(0)*3);
			vy1=pts.at(cells.at(*j).at(2)*3+1)-pts.at(cells.at(*j).at(0)*3+1);
			vy2=pts.at(cells.at(*j).at(1)*3+1)-pts.at(cells.at(*j).at(0)*3+1);
			vz1=pts.at(cells.at(*j).at(2)*3+2)-pts.at(cells.at(*j).at(0)*3+2);
			vz2=pts.at(cells.at(*j).at(1)*3+2)-pts.at(cells.at(*j).at(0)*3+2);
			
		float nxt=(vy1*vz2-vz1*vy2);
			float nyt=(vz1*vx2-vx1*vz2);
			float nzt=(vx1*vy2-vy1*vx2);

			nx+=nxt;
			ny+=nyt;
			nz+=nzt;
		}
		
		float d=sqrt(nx*nx+ny*ny+nz*nz);
		vx.push_back(nx/d);
		vy.push_back(ny/d);
		vz.push_back(nz/d);
	}
}

template<class T>
vector< vector<T> >  first_mesh::findNeighbours(const vector< vector<T> > & triangles, const unsigned int & N)
{
	vector< vector<T> > neighbours;
	for (unsigned int i=0; i<N ; i++)
	{
		vector<T> pt_neigh;
		for (typename vector< vector<T> >::const_iterator tri_i=triangles.begin(); tri_i!=triangles.end(); tri_i++)
			for ( typename vector<T> ::const_iterator vert_i= tri_i->begin(); vert_i!=tri_i->end(); vert_i++)
				if ( static_cast<unsigned int>(*vert_i) == i )//if point belongs to triangle add neigbours
				{	
					for ( typename vector<T> ::const_iterator vert_ii= tri_i->begin(); vert_ii!=tri_i->end(); vert_ii++)
						if  ( static_cast<unsigned int>(*vert_ii) != i )
						{
							bool found=false;
							for ( typename vector<T>::iterator temp_i=pt_neigh.begin(); temp_i!=pt_neigh.end();temp_i++)
									if ( (*vert_ii) == (*temp_i) ) found=true;
							if (!found)
								pt_neigh.push_back(*vert_ii);
						}
					break;
				}
	
		neighbours.push_back(pt_neigh);
		}
	return neighbours;
}

template<class T>
vector< vector<T> >  first_mesh::findNeighbourTriangles(const vector< vector<T> > & cells, const unsigned int & N)
{
	vector< vector<T> > neighbours;
	for (unsigned int i=0; i<N ; i++)
	{
		vector<T> tri_neigh;
		int count=0;
		for (typename vector< vector<T> >::const_iterator tri_i=cells.begin(); tri_i!=cells.end(); tri_i++,count++)
			for ( typename vector<T> ::const_iterator tri_j= tri_i->begin(); tri_j!=tri_i->end(); tri_j++)
				if ( static_cast<unsigned int>(*tri_j) == i )//if point belongs to triangle add neigbours
				{	
					tri_neigh.push_back(count);
					break;
				}
	
		neighbours.push_back(tri_neigh);
		}
	return neighbours;
}



template<class T,class T2>
void  first_mesh::medium_neighbours(const vector<T> & pts, const vector< vector<T2> > & neighbours, const vector< vector<T2> > & cells,\
									vector<T> & vx, vector<T> & vy, vector<T> & vz )
{
	
	vx.clear();
	vy.clear();
	vz.clear();
	//calculate the normals for each triangle then references
		
	for ( typename vector< vector<T2> >::const_iterator i=neighbours.begin();i!=neighbours.end();i++)
	{
		float nx=0,ny=0,nz=0;
			
		for (typename vector<T2>::const_iterator j=i->begin();j!=i->end();j++)
		{
			nx+=pts.at( (*j)*3 );
			ny+=pts.at( (*j)*3 + 1);
			nz+=pts.at( (*j)*3 + 2);
		}
	
		vx.push_back(nx/i->size());
		vy.push_back(ny/i->size());
		vz.push_back(nz/i->size());
	}
}


template<class T, class T2>
void first_mesh::maxTriangle(const vector<T> & pts, const vector< vector<T2> > & localTri, const vector< vector<T2> > & cells,vector<T> & vx, vector<T> & vy, vector<T> & vz )
{
  vx.clear();
	vy.clear();
	vz.clear();

	typename vector<T>::const_iterator pts_i = pts.begin(); 
	
	for (typename vector< vector<T2> >::const_iterator i=localTri.begin();i!=localTri.end();i++ , pts_i+=3)
	{
		float vx1,vy1,vz1,vx2,vy2,vz2;
		float max_area=0, va_x=0, va_y=0, va_z=0;

		for (typename vector<T2>::const_iterator j=i->begin();j!=i->end();j++)
		{
			vx1=pts.at(cells.at(*j).at(2)*3)-pts.at(cells.at(*j).at(0)*3);
			vx2=pts.at(cells.at(*j).at(1)*3)-pts.at(cells.at(*j).at(0)*3);
			vy1=pts.at(cells.at(*j).at(2)*3+1)-pts.at(cells.at(*j).at(0)*3+1);
			vy2=pts.at(cells.at(*j).at(1)*3+1)-pts.at(cells.at(*j).at(0)*3+1);
			vz1=pts.at(cells.at(*j).at(2)*3+2)-pts.at(cells.at(*j).at(0)*3+2);
			vz2=pts.at(cells.at(*j).at(1)*3+2)-pts.at(cells.at(*j).at(0)*3+2);
			
			float nxt=(vy1*vz2-vz1*vy2);
			float nyt=(vz1*vx2-vx1*vz2);
			float nzt=(vx1*vy2-vy1*vx2);

			float area=0.5* sqrt(nxt*nxt+nyt*nyt+nzt*nzt); 
			
			//find bisecting vector
			//cout<<"area "<<area<<endl;
			
			if ( area > max_area )
			{	
				va_x=(      pts.at(cells.at(*j).at(0)*3) + pts.at(cells.at(*j).at(1)*3)      + pts.at(cells.at(*j).at(2)*3 )     ) / 3 - *pts_i ;
				va_y=( pts.at(cells.at(*j).at(0)*3 + 1 ) + pts.at(cells.at(*j).at(1)*3 + 1 ) + pts.at(cells.at(*j).at(2)*3 + 1 ) ) / 3 - *(pts_i+1);
				va_z=( pts.at(cells.at(*j).at(0)*3 + 2 ) + pts.at(cells.at(*j).at(1)*3 + 2 ) + pts.at(cells.at(*j).at(2)*3 + 2 ) ) / 3 - *(pts_i+2);
				max_area=area;
			//	cout<<"MAX AREA "<<max_area<<endl;
			}
			
		}
		float norm=sqrt(va_x*va_x + va_y*va_y + va_z*va_z );
		vx.push_back( va_x/norm * max_area );
		vy.push_back( va_y/norm * max_area );
		vz.push_back( va_z/norm * max_area );
	}

}

template<class T>
bool first_mesh::triangle_intersection(const T* V10, const T* V11, const T* V12, const T* V20, const T* V21, const T* V22 )
{//this is the Moller algorithm for fast traingle intersection
	//calculate triangle plane 2
	//cout<<"test triangles"<<endl;
	short eq=0;
	if ( pt_equal(V10,V20) || pt_equal(V10,V21) ||  pt_equal(V10,V22) ) eq++;
	if ( pt_equal(V11,V20) || pt_equal(V11,V21) ||  pt_equal(V11,V22) ) eq++;
	if ( pt_equal(V12,V20) || pt_equal(V12,V21) ||  pt_equal(V12,V22) ) eq++;

	if (eq==1) return false;
	 
	T* N2 = cross_prod<T>( vec_sub(V21,V20),  vec_sub(V22,V20) );
	T  d2 = -dot_prod(N2,V20);
	
//	T testV[]={1,2,3};
//	T testV2[]= {4,5,6};
//	T* test=cross_prod<T>(testV,testV2);
	//T* N2=cross_prod<T>(testV,testV2);

	//cout<<"test "<<test[0]<<" "<<test[1]<<" "<<test[2]<<" "<<dot_prod(test,testV2)<<endl;
	//calculate distance of vertices from plane 1
	T dist10=dot_prod(N2,V10)+d2;
	T dist11=dot_prod(N2,V11)+d2;
	T dist12=dot_prod(N2,V12)+d2;
	//cout<<"DIST AGAIN "<<dist10<<" "<<dot_prod(N2,V22)+d2<<" "<<dot_prod(N2,V20)+d2<<" "<<dot_prod(N2,V21)+d2<<endl;
	if ( ((dist10>=0) && (dist11>=0) && (dist12>=0)) || ((dist10<=0) && (dist11<=0) && (dist12<=0)) )
	{
		delete N2;
		return false;//not all points are on same side of plan
}
//	}else if (( dist10 == dist11 ) && (dist11 == dist12) && ( dist12 == 0 )) //test for coplanar
//	{
//		delete N2;
//		return false;
//	}else
	if ( (dist10==0) || (dist11==0) || (dist12==0) )
	{
		delete N2;
		return false;
	}else{
		//calculate intersection line
		T* N1 = cross_prod<T>( vec_sub(V11,V10),  vec_sub(V12,V10) );
		T  d1 = -dot_prod(N1,V10);

		T dist20=dot_prod(N1,V20)+d1;
		T dist21=dot_prod(N1,V21)+d1;
		T dist22=dot_prod(N1,V22)+d1;
	//	if ( dist20 == dist21 == dist22 == 0 )
	//	{
	//	cout<<"dist1 "<<dist10<<" "<<dist11<<" "<<dist12<<endl;
	//		cout<<"not intersect 1"<<endl;

	//		delete N1;
	//		delete N2;
	//		return false;
	//	}
			if ( ((dist20>=0) && (dist21>=0) && (dist22>=0)) || ((dist20<=0) && (dist21<=0) && (dist22<=0)) )
	{
		delete N2;
		delete N1;
		return false;//not all points are on same side of plan
}
		//claulcate line in=tersection (without intercept)
		T* D = cross_prod<T>(N1,N2);
		//calculate max axis to project onto
/*		short max_axis = 0;
		if (abs(D[1])>abs(D[0]))
			max_axis=1;
		else if (abs(D[2])>abs(D[max_axis]))
			max_axis=2;
		
		delete N1;
		delete N2;
		delete D;
		
		T p10 = V10[max_axis];
		T p11 = V11[max_axis];
		T p12 = V12[max_axis];

		T p20 = V20[max_axis];
		T p21 = V21[max_axis];
		T p22 = V22[max_axis];
*/
		T p10 = dot_prod(D,V10);//[max_axis];
		T p11 = dot_prod(D,V11);
		T p12 =dot_prod(D,V12);

		T p20 =dot_prod(D,V20);
		T p21 =dot_prod(D,V21);
		T p22 = dot_prod(D,V22);
		
if (abs(dist10)<1e-4) dist10=0;
if (abs(dist11)<1e-4) dist11=0;
if (abs(dist12)<1e-4) dist12=0;
if (abs(dist20)<1e-4) dist20=0;
if (abs(dist21)<1e-4) dist21=0;
if (abs(dist22)<1e-4) dist22=0;

	//	cout<<"dist1 "<<dist10<<" "<<dist11<<" "<<dist12<<" "<<dist20<<" "<<dist21<<" "<<dist22<<endl;

if ( (dist20==0) || (dist21==0) || (dist22==0) )
	{
	//cout<<"not intersect 2"<<endl;
		delete N2;
		delete N1;
		return false;
	}

		T t11, t12, t21, t22; //intersection parameters

		//get triangle 1 interval
		if ( ((dist10>=0) && (dist11>=0)) || ((dist10<=0) && (dist11<=0)) )
		{
			t11= p10 - (p10 - p12)*(dist10/(dist10 - dist12));
			t12= p11 - (p11 - p12)*(dist11/(dist11 - dist12));
		}else if ( ((dist10>=0) && (dist12>=0)) || ((dist10<=0) && (dist12<=0)) )
		{
			t11= p10 - (p10 - p11)*(dist10/(dist10 - dist11));
			t12= p12 - (p12 - p11)*(dist12/(dist12 - dist11));
		}else
		{
			t11= p11 - (p11 - p10)*(dist11/(dist11 - dist10));
			t12= p12 - (p12 - p10)*(dist12/(dist12 - dist10));
		}
		
		//get triangle 2 interval
		if ( ((dist20>=0) && (dist21>=0)) || ((dist20<=0) && (dist21<=0)) )
		{
			t21= p20 - (p20 - p22)*(dist20/(dist20 - dist22));
			t22= p21 - (p21 - p22)*(dist21/(dist21 - dist22));
		}else if ( ((dist20>=0) && (dist22>=0)) || ((dist20<=0) && (dist22<=0)) )
		{
			t21= p20 - (p20 - p21)*(dist20/(dist20 - dist21));
			t22= p22 - (p22 - p21)*(dist22/(dist22 - dist21));
		}else
		{
			t21= p21 - (p21 - p20)*(dist21/(dist21 - dist20));
			t22= p22 - (p22 - p20)*(dist22/(dist22 - dist20));
		}
		
//		cout<<V10[0]<<" "<<V10[1]<<" "<<V10[2]<<endl;
//		cout<<V11[0]<<" "<<V11[1]<<" "<<V11[2]<<endl;
//		cout<<V12[0]<<" "<<V12[1]<<" "<<V12[2]<<endl;
//		cout<<V20[0]<<" "<<V20[1]<<" "<<V20[2]<<endl;
//		cout<<V21[0]<<" "<<V21[1]<<" "<<V21[2]<<endl;
//		cout<<V22[0]<<" "<<V22[1]<<" "<<V22[2]<<endl;
//		cout<<"dist "<<dist10<<" "<<dist11<<" "<<dist12<<" "<<dist20<<" "<<dist21<<" "<<dist22<<endl;
//		cout<<"intervals "<<t11<<" "<<t12<<" "<<t21<<" "<<t22<<endl;
		//check interval overlap 
		if ( ((t21<=t11) && (t21<=t12) &&  (t22<=t11) && (t22<=t12)) || \
			 ((t21>=t11) && (t21>=t12) &&  (t22>=t11) && (t22>=t12)) )
			return false;


	}
	return true;
}


template<class T,class T2>
bool first_mesh::self_intersection_test(const vector< vector<T2> > & tri_cells, const vector<T> & pts  )
{
	vector<T> pts_x, pts_y, pts_z;
	for (typename vector<T>::const_iterator i=pts.begin();i!=pts.end();i+=3)
	{
		pts_x.push_back(*i);
		pts_y.push_back(*(i+1));
		pts_z.push_back(*(i+2));

	}
	
	typename vector< vector<T2> >::const_iterator start_j=tri_cells.begin()+1;

	T* V10 = new T[3];
	T* V11 = new T[3];
	T* V12 = new T[3];
	T* V20 = new T[3];
	T* V21 = new T[3];
	T* V22 = new T[3];
	for (typename vector< vector<T2> >::const_iterator i=tri_cells.begin(); i!=tri_cells.end() ;i++)
	{	//for each triangle, compare against all others
		V10[0]=pts_x.at(i->at(0));
		V10[1]=pts_y.at(i->at(0));
		V10[2]=pts_z.at(i->at(0));
		V11[0]=pts_x.at(i->at(1));
		V11[1]=pts_y.at(i->at(1));
		V11[2]=pts_z.at(i->at(1));
		V12[0]=pts_x.at(i->at(2));
		V12[1]=pts_y.at(i->at(2));
		V12[2]=pts_z.at(i->at(2));
		
		for (typename vector< vector<T2> >::const_iterator j=start_j; j!=tri_cells.end() ;j++)
		{
			V20[0]=pts_x.at(j->at(0));
			V20[1]=pts_y.at(j->at(0));
			V20[2]=pts_z.at(j->at(0));
			V21[0]=pts_x.at(j->at(1));
			V21[1]=pts_y.at(j->at(1));
			V21[2]=pts_z.at(j->at(1));
			V22[0]=pts_x.at(j->at(2));
			V22[1]=pts_y.at(j->at(2));
			V22[2]=pts_z.at(j->at(2));
			//cout<<"tri "<<V20[0]<<" "<<V21[1]<<" "<<V22[0]<<endl;
			if ( triangle_intersection(V10, V11, V12, V20, V21, V22) )
			{
				delete V10;
				delete V11;
				delete V12;
				delete V20;
				delete V21;
				delete V22;
				return true;
			}
		}
		start_j++;
	}
	delete V10;
				delete V11;
				delete V12;
				delete V20;
				delete V21;
				delete V22;
				return false;
				
}



template void first_mesh::normal<float, short>(const vector<float> & pts, const vector< vector<short> > & localTri, \
										const vector< vector<short> > & cells,vector<float> & vx, vector<float> & vy, vector<float> & vz );
template void first_mesh::normal<double, short>(const vector<double> & pts, const vector< vector<short> > & localTri, \
										const vector< vector<short> > & cells,vector<double> & vx, vector<double> & vy, vector<double> & vz );
template void first_mesh::normal<float, unsigned int>(const vector<float> & pts, const vector< vector<unsigned int> > & localTri, \
										const vector< vector<unsigned int> > & cells,vector<float> & vx, vector<float> & vy, vector<float> & vz );
template void first_mesh::normal<double, unsigned int>(const vector<double> & pts, const vector< vector<unsigned int> > & localTri, \
										const vector< vector<unsigned int> > & cells,vector<double> & vx, vector<double> & vy, vector<double> & vz );



//medium neighbours
template void  first_mesh::medium_neighbours<float,short>(const vector<float> & pts, const vector< vector<short> > & localTri, const vector< vector<short> > & cells,\
																  vector<float> & vx, vector<float> & vy, vector<float> & vz );
template void  first_mesh::medium_neighbours<double,short>(const vector<double> & pts, const vector< vector<short> > & localTri, const vector< vector<short> > & cells,\
																  vector<double> & vx, vector<double> & vy, vector<double> & vz );
template void  first_mesh::medium_neighbours<float,unsigned int>(const vector<float> & pts, const vector< vector<unsigned int> > & localTri, const vector< vector<unsigned int> > & cells,\
																  vector<float> & vx, vector<float> & vy, vector<float> & vz );
template void  first_mesh::medium_neighbours<double,unsigned int>(const vector<double> & pts, const vector< vector<unsigned int> > & localTri, const vector< vector<unsigned int> > & cells,\
																  vector<double> & vx, vector<double> & vy, vector<double> & vz );
																  

//find neighbours
template vector< vector<short> >  first_mesh::findNeighbours<short>(const vector< vector<short> > & triangles, const unsigned int & N);
template vector< vector<unsigned int> >  first_mesh::findNeighbours<unsigned int>(const vector< vector<unsigned int> > & triangles, const unsigned int & N);
template vector< vector<int> >  first_mesh::findNeighbours<int>(const vector< vector<int> > & triangles, const unsigned int & N);
																  
template vector< vector<short> >  first_mesh::findNeighbourTriangles<short>(const vector< vector<short> > & triangles, const unsigned int & N);
template vector< vector<unsigned int> >  first_mesh::findNeighbourTriangles<unsigned int>(const vector< vector<unsigned int> > & triangles, const unsigned int & N);
template vector< vector<int> >  first_mesh::findNeighbourTriangles<int>(const vector< vector<int> > & triangles, const unsigned int & N);			

template void first_mesh::maxTriangle<float, short>(const vector<float> & pts, const vector< vector<short> > & localTri, const vector< vector<short> > & cells,vector<float> & vx, vector<float> & vy, vector<float> & vz );
template void first_mesh::maxTriangle<float, unsigned int>(const vector<float> & pts, const vector< vector<unsigned int> > & localTri, const vector< vector<unsigned int> > & cells,vector<float> & vx, vector<float> & vy, vector<float> & vz );
template void first_mesh::maxTriangle<double, short>(const vector<double> & pts, const vector< vector<short> > & localTri, const vector< vector<short> > & cells,vector<double> & vx, vector<double> & vy, vector<double> & vz );
template void first_mesh::maxTriangle<double, unsigned int>(const vector<double> & pts, const vector< vector<unsigned int> > & localTri, const vector< vector<unsigned int> > & cells,vector<double> & vx, vector<double> & vy, vector<double> & vz );


template bool first_mesh::self_intersection_test<float, short>(const vector< vector<short> > & tri_cells, const vector<float> & pts  );
template bool first_mesh::self_intersection_test<float, unsigned int>(const vector< vector<unsigned int> > & tri_cells, const vector<float> & pts  );
template bool first_mesh::self_intersection_test<double, short>(const vector< vector<short> > & tri_cells, const vector<double> & pts  );
template bool first_mesh::self_intersection_test<double, unsigned int>(const vector< vector<unsigned int> > & tri_cells, const vector<double> & pts  );


/*
template short pointTest<float>(const float & px_0, const float & py_0, const float & pz_0, const float & px_1, const float & py_1, const float & pz_1, \
				const float & px_2, const float & py_2, const float & pz_2, const float & nx, const float & ny, const float & nz, \
				float & p0_x, float & p0_y, float & p0_z, float & dx, float & dy, float & dz);
template short pointTest<double>(const double & px_0, const double & py_0, const double & pz_0, const double & px_1, const double & py_1, const double & pz_1, \
				const double & px_2, const double & py_2, const double & pz_2, const double & nx, const double & ny, const double & nz,\
				double & p0_x, double & p0_y, double & p0_z, double & dx, double & dy, double & dz);
*/


}

