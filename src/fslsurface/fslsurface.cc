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
#include "fslsurface.h"
#include <sstream>
#include <cmath>
#include <set>

#define GL_GLEXT_PROTOTYPES

//#include "fslvtkio/fslvtkio.h"

//extern "C" {	
//#include <gifti/gifti_io.h>
//}
using namespace std;
//using namespace fslvtkio;
//using namespace NEWIMAGE;
//using namespace mesh;

//string encoding_str[5] = { "GIFTI_ENCODING_UNDEF", "GIFTI_ENCODING_ASCII", "GIFTI_ENCODING_B64BIN", "GIFTI_ENCODING_B64GZ","GIFTI_ENCODING_EXTBIN"};
//string ind_ord_str[5] = { "GIFTI_IND_ORD_UNDEF", "GIFTI_IND_ORD_ROW_MAJOR", "GIFTI_IND_ORD_COL_MAJOR"};
//string endian_str[3] = { "GIFTI_ENDIAN_UNDEF", "GIFTI_ENDIAN_BIG", "GIFTI_ENDIAN_LITTLE"};
//string data_type_str[] = { }


namespace fslsurface_name {

  string encoding_str[5] = { "GIFTI_ENCODING_UNDEF", "GIFTI_ENCODING_ASCII", "GIFTI_ENCODING_B64BIN", "GIFTI_ENCODING_B64GZ","GIFTI_ENCODING_EXTBIN"};
  string ind_ord_str[5] = { "GIFTI_IND_ORD_UNDEF", "GIFTI_IND_ORD_ROW_MAJOR", "GIFTI_IND_ORD_COL_MAJOR"};
  string endian_str[3] = { "GIFTI_ENDIAN_UNDEF", "GIFTI_ENDIAN_BIG", "GIFTI_ENDIAN_LITTLE"};




	//	void giftWrapper_read_Image(const std::string & filename, gifti_image* gii_surf)
	//	{
	//		gii_surf = gifti_read_image(filename.c_str(), 1);
	//	}
	//--------------------end operators----------------------------------//
    template<class T, class T2>
    void fslSurface<T,T2>::operator*=(const T & mul ) {
        
        for ( typename vector< vertex<T> >::iterator i= vertices.begin(); i!=vertices.end();++i)
        {
            i->x*=mul;//+=i2->x;
            i->y*=mul;//+=i2->y;
            i->z*=mul;//+=i2->z;
            
        }
        calculateNormals();
        //   return *this;
    }
    
    template<class T, class T2>
    void fslSurface<T,T2>::operator/=(const T & div ) {

        for ( typename vector< vertex<T> >::iterator i= vertices.begin(); i!=vertices.end();++i)
        {
            i->x/=div;//+=i2->x;
            i->y/=div;//+=i2->y;
            i->z/=div;//+=i2->z;
            
        }
        calculateNormals();
        //   return *this;
    }

    template<class T, class T2>
    void fslSurface<T,T2>::operator+=(const  fslSurface<T,T2> &surf) {
        
        if (getNumberOfVertices() != surf.getNumberOfVertices())
            throw fslSurfaceException("");
        typename vector< vertex<T> >::const_iterator i2= surf.const_vbegin();
        for ( typename vector< vertex<T> >::iterator i= vertices.begin(); i!=vertices.end();++i,++i2)
        {
            i->x+=i2->x;
            i->y+=i2->y;
            i->z+=i2->z;
            
        }
        calculateNormals();
     //   return *this;
    }
    template<class T, class T2>
    void fslSurface<T,T2>::operator-=(const  fslSurface<T,T2> &surf) {
        
        if (getNumberOfVertices() != surf.getNumberOfVertices())
            throw fslSurfaceException("");
        typename vector< vertex<T> >::const_iterator i2= surf.const_vbegin();
        for ( typename vector< vertex<T> >::iterator i= vertices.begin(); i!=vertices.end();++i,++i2)
        {
            i->x-=i2->x;
            i->y-=i2->y;
            i->z-=i2->z;

        }
        calculateNormals();
        //   return *this;
    }
    
    template<class T, class T2>
    fslSurface<T,T2> fslSurface<T,T2>::operator+(const  fslSurface<T,T2>& surf ) {
        
        fslSurface<T,T2> surf_result(surf);
        if (getNumberOfVertices() != surf.getNumberOfVertices())
            throw fslSurfaceException("");
        typename vector< vertex<T> >::const_iterator i2= surf.const_vbegin();
        for ( typename vector< vertex<T> >::iterator i=  surf_result.vbegin(); i!= surf_result.vend();++i,++i2)
        {
            i->x+=i2->x;
            i->y+=i2->y;
            i->z+=i2->z;
            
        }
        calculateNormals();
        //   return *this;
        return surf_result;
    }
    template<class T, class T2>
    fslSurface<T,T2> fslSurface<T,T2>::operator-(const  fslSurface<T,T2>& surf ) {
        
        fslSurface<T,T2> surf_result(surf);
        if (getNumberOfVertices() != surf.getNumberOfVertices())
            throw fslSurfaceException("");
        typename vector< vertex<T> >::const_iterator i2= surf.const_vbegin();
        for ( typename vector< vertex<T> >::iterator i=  surf_result.vbegin(); i!= surf_result.vend();++i,++i2)
        {
            i->x-=i2->x;
            i->y-=i2->y;
            i->z-=i2->z;
            
        }
        calculateNormals();
        //   return *this;
        return surf_result;
    }
    
    template<class T, class T2>
    fslSurface<T,T2> fslSurface<T,T2>::operator*(const T& mul ) {
        
        fslSurface<T,T2> surf_result(*this);
       
        for ( typename vector< vertex<T> >::iterator i=  surf_result.vbegin(); i!= surf_result.vend();++i)
        {
            i->x*=mul;
            i->y*=mul;
            i->z*=mul;
            
        }
        calculateNormals();
        //   return *this;
        return surf_result;
    }
    
    template<class T, class T2>
    fslSurface<T,T2> fslSurface<T,T2>::operator/(const  T& div ) {
        
        fslSurface<T,T2> surf_result(*this);
       
        for ( typename vector< vertex<T> >::iterator i=  surf_result.vbegin(); i!= surf_result.vend();++i)
        {
            i->x/=div;
            i->y/=div;
            i->z/=div;
            
        }
        calculateNormals();
        //   return *this;
        return surf_result;
    }
    //--------------------operators----------------------------------//

	//float3 funcs
	template<class T>
	vec3<T> subtract( const vec3<T> & a, const vec3<T> & b )
	{
		return vec3<T>( a.x-b.x , a.y-b.y , a.z-b.z );
	}
	template<class T>
	vec3<T> normal(const vec3<T> & v1, const vec3<T> & v2 )
    {
        vec3<T> normal;
        normal.x=v1.y*v2.z - v1.z*v2.y;
        normal.y=v1.z*v2.x - v1.x*v2.z;
        normal.z=v1.x*v2.y - v1.y*v2.x;
        T l=sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
        normal.x/=l;
        normal.y/=l;
        normal.z/=l;
        return normal;
    }
	
    template<class T, class T2>
	fslSurface<T,T2>::fslSurface( const fslSurface<T,T2>& surf ){
        vertices = surf.vertices;
        faces = surf.faces;
    
		 tangents = surf.tangents;
        
		adj_tris = surf.adj_tris;
		adj_verts = surf.adj_verts;
        
		v_coord_sys = surf.v_coord_sys; 
		csys_dspace = surf.csys_dspace;
		v_csys_xfmspace = surf.v_csys_xfmspace; 
        
		topology_name = surf.topology_name;  
		geometry_name = surf.geometry_name;  
		anatomical_name = surf.anatomical_name;  
		anatomical_name2 = surf.anatomical_name2;  

		//if intending to use the openGL use unsigned int (suggested anyways)
		vector_data = surf.vector_data;
		scalar_data = surf.scalar_data;
		scalar_names = surf.scalar_names;
		scalar_indices = surf.scalar_indices;
        
		dataTable = surf.dataTable;
        
		element_sizes = surf.element_sizes;
		
		num_format = surf.num_format;
		format = surf.format;
		//FileType file_type;
		ply_version = surf.ply_version;
		
		N_vertices = surf.N_vertices;
        N_triangles = surf.N_triangles;
		cog = surf.cog;
		datatype_str = surf.datatype_str;
    }
    
    
	template<class T, class T2>
	fslSurface<T,T2>::fslSurface(){
		//endian_format = fslSurface::machineEndianness();
		//cout<<"constructor"<<endl;
		N_vertices=0;
		N_triangles=0;
		csys_dspace="UNKNOWN";
//		vector<double> csys(16,0);
//		csys[0]=1;
//		csys[4]=1;
//		csys[7]=1;
//		csys[11]=1;
//		v_coord_sys.push_back(csys);
//		v_csys_xfmspace.push_back("UNKNOWN");
		//file_type=Unknown;
		//		datatype_str[NIFTI_TYPE_UINT8]= "NIFTI_TYPE_UINT8";
		//datatype_str[NIFTI_TYPE_INT32]= "NIFTI_TYPE_INT32";
		//datatype_str[NIFTI_TYPE_FLOAT32]= "NIFTI_TYPE_FLOAT32";
		
		
	}
	
//	template<class T, class T2>
//	fslSurface<T,T2>::fslSurface( const string & filename){
//		//endian_format = fslSurface::machineEndianness();
//		//cout<<"constructor"<<endl;
//		N_vertices=0;
//		N_triangles=0;
//		//file_type=Unknown;
//		datatype_str[NIFTI_TYPE_UINT8]= "NIFTI_TYPE_UINT8";
//		datatype_str[NIFTI_TYPE_INT32]= "NIFTI_TYPE_INT32";
//		datatype_str[NIFTI_TYPE_FLOAT32]= "NIFTI_TYPE_FLOAT32";
//		read_surface<T,T2>(this, filename);
//	//	v_coord_sys.push_back();
		
//	}
	
	template<class T, class T2>
	fslSurface<T,T2>::~fslSurface(){
		//cout<<"fslsurface detsructor "<<endl;
	}
    
    template<class T, class T2>
    vector<T> fslSurface<T,T2>::getBounds() const
    {
        vector<T> bounds(6,0);
        bounds[0] = bounds[2] = bounds[4] = +1e+12;
        bounds[1] = bounds[3] = bounds[5] = -1e+12;

        for ( typename vector< vertex<T> >::const_iterator i_v = vertices.begin() ; i_v != vertices.end();++i_v)
        {
            if (bounds[0]> i_v->x) bounds[0]=i_v->x;
            if (bounds[2]> i_v->y) bounds[2]=i_v->y;
            if (bounds[4]> i_v->z) bounds[4]=i_v->z;
            if (bounds[1] < i_v->x) bounds[1]=i_v->x;
            if (bounds[3] < i_v->y) bounds[3]=i_v->y;
            if (bounds[5] < i_v->z) bounds[5]=i_v->z;

        }
            
        return bounds;
        
    }

    template<class T, class T2>
    void fslSurface<T,T2>::copyNormalsToVectors(const unsigned int & index)
    {
        vector<T> new_vecs(N_vertices*3);
        
        typename vector<T>::iterator i_new_v = new_vecs.begin();
        for ( typename vector< vertex<T> >::iterator i_v = vertices.begin(); i_v != vertices.end(); ++i_v, i_new_v+=3)
        {
            *i_new_v     = i_v->nx;
            *(i_new_v+1) = i_v->ny;
            *(i_new_v+2) = i_v->nz;
        }
        
        vector_data.insert(vector_data.begin()+index, new_vecs);
        vector_names.insert(vector_names.begin()+index, "normals");
        
    }

    
    template<class T, class T2>
    void fslSurface<T,T2>::copyVerticesToVectors(){
        copyVerticesToVectors("vectors");
    }

    template<class T, class T2>
    void fslSurface<T,T2>::copyVerticesToVectors( const string & name){//copie into first vector location
        vector<T> vec( vertices.size() * 3 );
        vector_names.push_back(name);
        typename vector<T>::iterator i_vec = vec.begin();
        for ( typename vector< vertex<T> >::iterator i_v = vertices.begin(); i_v != vertices.end(); ++i_v, i_vec+=3 )
        {
            *i_vec = i_v->x;
            *(i_vec+1) = i_v->y;
            *(i_vec+2) = i_v->z;

        }
        vector_data.insert(vector_data.begin(),vec);
//        for (typename vector<T>::iterator i = vector_data[0].begin(); i!= vector_data//[0].end(); ++i)
      //  {
        //    cout<<"vec "<<*i<<endl;
            
       // }
            
    }
	template<class T, class T2>
    float fslSurface<T,T2>::L2norm( const unsigned int & index0, const unsigned int & index1)
    {
       // vertex<T> v0=vertices[index0];
       // vertex<T> v1=vertices[index1];
     //   cout<<"l2norm "<<v0.x<<" "<<v1.x<<endl;
    //    cout<<"l2norm "<<v0.y<<" "<<v1.y<<endl;
     //   cout<<"l2norm "<<v0.z<<" "<<v1.z<<endl;
        T difx = vertices[index0].x-vertices[index1].x ;
        T dify = vertices[index0].y-vertices[index1].y ;
        T difz = vertices[index0].z-vertices[index1].z ;

        return sqrtf( difx*difx + dify*dify + difz*difz );

    }
    
    template<class T, class T2>
    fslsurface_name::vec3<T> fslSurface<T,T2>::getVertexCoord(const unsigned int & index)
    {
        return vec3<T>(vertices[index].x,vertices[index].y,vertices[index].z);
    }
    
    template<class T, class T2>
    void fslSurface<T,T2>::setVert(const unsigned int & index, const fslsurface_name::vec3<float> & vert_coord)
    {
        vertices[index].x = vert_coord.x;
        vertices[index].y = vert_coord.y;
        vertices[index].z = vert_coord.z;

    }

    
    template<class T, class T2>
    fslsurface_name::vec3<float> fslSurface<T,T2>::subtractVerts(const unsigned int & index0, const unsigned int & index1)
    {
        return vec3<float>(vertices[index0].x-vertices[index1].x, \
                      vertices[index0].y-vertices[index1].y, \
                      vertices[index0].z-vertices[index1].z) ;
        
    }

    
    
	template<class T, class T2>
	void fslSurface<T,T2>::copyVertices( const fslSurface<T,T2> & surf ){
        
        vertices=surf.vertices;

    }

	template<class T, class T2>
	void fslSurface<T,T2>::setVertices( const vector<T> &  verts )
	{
		vertices.clear();
		for ( typename vector<T>::const_iterator i = verts.begin() ; i != verts.end() ; ++i)
		{
			T x = *(i++);
			T y = *(i++);
			vertices.push_back(vertex<T>(x,y,*i));
			
		}
		//cout<<"nverts "<<N_vertices<<endl;
		N_vertices = vertices.size();
	}
	template<class T, class T2>
	typename vector< vertex<T> >::iterator fslSurface<T,T2>::vbegin()
	{		
		return vertices.begin();
	}
    template<class T, class T2>
	typename vector< vertex<T> >::iterator fslSurface<T,T2>::vend()
	{
		return vertices.end();
	}
    template<class T, class T2>
	typename vector< vertex<T> >::const_iterator fslSurface<T,T2>::const_vbegin() const
	{		
		return vertices.begin();
	}

    template<class T, class T2>
	typename vector< vertex<T> >::const_iterator fslSurface<T,T2>::const_vend() const
	{
		return vertices.end();
	}
    
    
    //face access 
    template<class T, class T2>
    typename vector< T2 >::const_iterator fslSurface<T,T2>::const_facebegin() const
    {
        return faces.begin();
    }
    
    template<class T, class T2>
    typename vector< T2 >::const_iterator fslSurface<T,T2>::const_faceend() const
    {
        return faces.end();
    }

    template<class T, class T2>
    vec3<T2> fslSurface<T,T2>::getFace( const unsigned int & index, int face_size ) const
    {
        return vec3<T2>(faces[index*face_size], faces[index*face_size+1], faces[index*face_size+2]);
        
    }

    template<class T, class T2>
    vector<T2> fslSurface<T,T2>::getFaces() const
    {
        return faces;
    }
    template<class T, class T2>
    vector< vec3<T> > fslSurface<T,T2>::getFaceVertices( const unsigned int & index, int face_size ) const
    {
        vector< vec3<T> > verts(face_size);
        typename vector< vec3<T> >::iterator i_v = verts.begin();
        for (int i=0; i<face_size; ++i,++i_v){
            unsigned int v_index = faces[index+i];
            *i_v = vec3<T>(vertices[v_index].x,vertices[v_index].y,vertices[v_index].z); 
        }
   
        return verts;
    }
    
    template<class T, class T2>
    std::vector< T > fslSurface<T,T2>::getFaceVerticesUnwrapped( const unsigned int & index, int face_size ) const
    {
        vector< T > verts(3*face_size);
        typename vector<T>::iterator i_v = verts.begin();
        for (int i=0; i<face_size; ++i,i_v+=3){
            
            unsigned int v_index = faces[index+i];
            *i_v    = vertices[v_index].x;
            *(i_v +1) = vertices[v_index].y;
            *(i_v +2) = vertices[v_index].z;
            
        }
        
        return verts;
    }

    
    
    template<class T, class T2>
    typename std::vector< T >::const_iterator  fslSurface<T,T2>::const_scbegin(const unsigned int & index ) const
    {
        return scalar_data[index].begin();
    }
    
    template<class T, class T2>
    typename std::vector< T >::const_iterator  fslSurface<T,T2>::const_scend(const unsigned int & index) const
    {
        return scalar_data[index].end();
    }

    
    
    template<class T, class T2>
	typename vector< T >::const_iterator fslSurface<T,T2>::const_vecbegin(const unsigned int & index ) const
	{		
		return vector_data[index].begin();
	}
    
    template<class T, class T2>
	typename vector< T >::const_iterator fslSurface<T,T2>::const_vecend(const unsigned int & index) const
	{
		return vector_data[index].end();
	}
	
	
	template<class T, class T2>
	void fslSurface<T,T2>::setFaces( const vector<T2> &  faces_in )
	{
		faces.clear();
		faces.insert(faces.end(), faces_in.begin(), faces_in.end() );
		//assume all faces are trinagles
		N_triangles=faces.size()/3;
	}
	
	template<class T, class T2>
	void fslSurface<T,T2>::setFaces( const std::vector< std::vector<T2> > &  faces_in )
	{
		faces.clear();
		int count=0;
		for (typename vector< vector<T2> >::const_iterator i = faces_in.begin() ; i != faces_in.end(); ++i) {
			faces.insert(faces.end(),i->begin(), i->end());	
			count++;
		}
	//	cout<<"cout set faces "<<count<<" "<<faces.size()<<endl;
		//assume all faces are trinagles
		N_triangles=faces.size()/3;
	}
	
    template<class T, class T2>
    unsigned int fslSurface<T,T2>::getNumberOfCoordinateSystems()
    {
        return v_coord_sys.size();
    }
    
	template<class T, class T2>
	void fslSurface<T,T2>::copy_coordsystem(fslSurface<T,T2> &  surf_dest, const fslSurface<T,T2>& surf_src)
	{
		surf_dest.v_coord_sys.clear();
		surf_dest.v_csys_xfmspace.clear();
		
		surf_dest.csys_dspace = surf_src.csys_dspace;
		surf_dest.v_coord_sys = surf_src.v_coord_sys;
		surf_dest.v_csys_xfmspace = surf_src.v_csys_xfmspace;

	}

	
	template<class T, class T2>
	void fslSurface<T,T2>::addCoordSystem(const vector<double> & xfm, const string & xfm_space )
	{
		vector<double> csys(16);
        vector<double>::iterator i_csys = csys.begin();
		for ( vector<double>::const_iterator i = xfm.begin(); i!= xfm.end(); ++i,++i_csys)
            *i_csys = (*i);

		v_coord_sys.push_back(csys);
		v_csys_xfmspace.push_back(xfm_space);
	}
	template<class T, class T2>
	void fslSurface<T,T2>::addCoordSystem(const vector<float> & xfm, const string & xfm_space )
	{
		vector<double> csys(16);
		vector<double>::iterator i_csys = csys.begin();
		for ( vector<float>::const_iterator i = xfm.begin(); i!= xfm.end(); ++i,++i_csys)
		{
			//cout<<"coord "<<*i<<" "<<xfm.size()<<endl;
//			csys.push_back(static_cast<double> (*i));
			*i_csys = static_cast<double> (*i);
			//cout<<"coord2 "<<csys.back()<<" "<<endl;

		}
		v_coord_sys.push_back(csys);
		//cout<<"coord sys size "<<v_coord_sys.back().size()<<endl;
		v_csys_xfmspace.push_back(xfm_space);
	}
    
    template<class T, class T2>
    void fslSurface<T,T2>::copyCoordSystems(const fslSurface<T,T2> & surf)
    {
        v_coord_sys.clear(); 
        v_coord_sys.insert(v_coord_sys.end(),surf.v_coord_sys.begin(),surf.v_coord_sys.end());
        v_csys_xfmspace.clear(); 
        v_csys_xfmspace.insert(v_csys_xfmspace.end(),surf.v_csys_xfmspace.begin(),surf.v_csys_xfmspace.end());

        csys_dspace = surf.csys_dspace;

    }
    
    template<class T, class T2>
    void fslSurface<T,T2>::clearCoordSystems()
    {
        v_coord_sys.clear(); 
        v_csys_xfmspace.clear(); 

    }

    
    
     template<class T, class T2> 
    vector<double> fslSurface<T,T2>::getCoordinateSystem( const std::string & name)
    {
//        vector<string>::iterator i_name = v_csys_xfmspace
        vector< vector<double> >::iterator i_c = v_coord_sys.begin(); 
        for (vector<string>::iterator i = v_csys_xfmspace.begin(); i!= v_csys_xfmspace.end();++i,++i_c)
        {
            if ( (*i) == name)
                return *i_c;
        }
        return vector<double>();
        
    }

	
	template<class T, class T2>
	void fslSurface<T,T2>::computeAdjacentTriangles()
	{
		vector<T2> dummy;
		//adj_tris.resize(N_vertices,dummy);
        
	//	unsigned int tri_index=0;
//		for (typename vector< vector<T2> >::iterator i_faces = faces.begin(); i_faces != faces.end(); ++i_faces, ++tri_index) {
//			for (typename vector<T2>::iterator i_faces2 = i_faces->begin();  i_faces2 != i_faces->end(); ++i_faces2)
//				adj_tris[*i_faces].push_back(tri_index);
//		}
		
	}
	
	template<class T, class T2>
	void fslSurface<T,T2>::computeAdjacentVertices( const bool & bidirectional)
	{
        adj_verts_bi = bidirectional;
		//assuming trinagles
		list<T2> dummy;
		adj_verts.resize(N_vertices,dummy);
	//	unsigned int tri_index=0;
		//for (typename vector< vector<T2> >::iterator i_faces = faces.begin(); i_faces != faces.end(); ++i_faces) 
            for (typename vector<T2>::iterator i_faces = faces.begin(); i_faces != faces.end(); i_faces+=3)
        {
            
            //only add connection from lower to higher
            //triangles so only 3 connections
            
            //conn1
          //  T2 ind0 = adj_verts[*i_faces];
          //  T2 ind1 = adj_verts[*(i_faces+1)];
            //T2 ind0 = ( *i_faces < *(i_faces+1) ) *i_faces ? *(i_faces+1);
            // T2 ind1 = ( *i_faces > *(i_faces+1) ) *(i_faces+1) ? *i_faces;
          /*  if ( *i_faces > *(i_faces+1) )
                adj_verts[*(i_faces+1)].push_back(*i_faces);
            else
                adj_verts[*(i_faces)].push_back(*(i_faces+1));

            if ( *(i_faces+1) > *(i_faces+2) )
                adj_verts[*(i_faces+2)].push_back(*(i_faces+1));
            else
                adj_verts[*(i_faces+1)].push_back(*(i_faces+2));

            if ( *(i_faces+2) > *(i_faces) )
                adj_verts[*(i_faces)].push_back(*(i_faces+2));
            else
                adj_verts[*(i_faces+2)].push_back(*(i_faces));
            */
           // cout<<"push nack "<<*i_faces<<" "<<*(i_faces+1)<<endl;
            adj_verts[*i_faces].push_back(*(i_faces+1));
            adj_verts[*(i_faces+1)].push_back(*(i_faces+2));
            adj_verts[*(i_faces+2)].push_back(*i_faces);
        }
        //add other direction plus check to see if there
        
        int count=0;
        for (typename vector< list<T2> >::iterator i_faces = adj_verts.begin(); i_faces != adj_verts.end(); ++i_faces,++count) 
        {
            (i_faces)->sort();
         //   cout<<count<<" 1dir: ";
          //  for (typename list<T2>::iterator ii_faces = i_faces->begin(); ii_faces != i_faces->end(); ++ii_faces)
           // {
            //    cout<<*ii_faces<<" ";                
            //}
            //cout<<endl;
		}
        
        
        if (bidirectional)
        for (typename vector<T2>::iterator i_faces = faces.begin(); i_faces != faces.end(); i_faces+=3)
        {
            {
                //check to see if already been added
                bool found=0;
             //   cout<<"new "<<endl;
                for ( typename list<T2>::iterator i = adj_verts[*(i_faces+1)].begin(); i != adj_verts[*(i_faces+1)].end();++i)
                {
             //       cout<<"tri : "<<adj_verts[*(i_faces+1)].size()<<" "<<*(i_faces)<<" "<<*(i_faces+1)<<" "<<*(i_faces+2)<<" -- "<<*i<<" "<<*(i_faces+1)<<endl;
                    if ( *i == *i_faces )
                    {
                 //       cout<<"found"<<endl;
                        found=true;
                        break;
                    }
                }
                if (!found)
                    adj_verts[*(i_faces+1)].push_back(*i_faces);

                //----------------------------------------------------
                found=0;
                for ( typename list<T2>::iterator i = adj_verts[*(i_faces+2)].begin(); i != adj_verts[*(i_faces+2)].end();++i)
                    if ( *i == *(i_faces+1) )
                    {
                        found=true;
                        break;
                    }
                if (!found)
                    adj_verts[*(i_faces+2)].push_back(*(i_faces+1));
                //----------------------------------------------------
                found=0;
                for ( typename list<T2>::iterator i = adj_verts[*(i_faces)].begin(); i != adj_verts[*(i_faces)].end();++i)
                    if ( *i == *(i_faces+2) )
                    {
                        found=true;
                        break;
                    }
                if (!found)
                    adj_verts[*(i_faces)].push_back(*(i_faces+2));
              //  adj_verts[].push_back(*(i_faces+1));
                //adj_verts[*(i_faces)].push_back(*(i_faces+2));
            }
            
        }
        
         count=0;
        for (typename vector< list<T2> >::iterator i_faces = adj_verts.begin(); i_faces != adj_verts.end(); ++i_faces,++count) 
        {
            (i_faces)->sort();
            //cout<<count<<": ";
           // for (typename list<T2>::iterator ii_faces = i_faces->begin(); ii_faces != i_faces->end(); ++ii_faces)
            //{
             //   cout<<*ii_faces<<" ";                
            //}
            //cout<<endl;
		}
	}
	
	
	template<class T, class T2>
	float3 fslSurface<T,T2>::getCOG()
	{
		return cog;
	}
	
	template<class T,class T2>
	void fslSurface<T,T2>::calculateNormals( bool normalize, bool reverse_winding){
	//	cout<<"normals calculate normals"<<endl;
		
		for (typename vector<T2>::iterator i_faces = faces.begin() ; i_faces != faces.end() ; i_faces+=3)
		{ 
			
			vec3<float> p0 = vec3<float>( vertices[*i_faces].x	,	vertices[*i_faces].y	,	vertices[*i_faces].z	);
            vec3<float> p2 = vec3<float>(	vertices[*(i_faces+1)].x	,	vertices[*(i_faces+1)].y	,	vertices[*(i_faces+1)].z);
            vec3<float> p1 = vec3<float>(	vertices[*(i_faces+2)].x	,	vertices[*(i_faces+2)].y	,	vertices[*(i_faces+2)].z);
            vec3<float> n = normal<float>(subtract<float>(p0,p2),subtract<float>(p0,p1));
			vertices[	*(i_faces)	].nx  +=   n.x;
            vertices[	*(i_faces)	].ny  +=   n.y;
            vertices[	*(i_faces)	].nz  +=   n.z;
			
			vertices[*(i_faces+1)].nx  +=   n.x;
            vertices[*(i_faces+1)].ny  +=   n.y;
            vertices[*(i_faces+1)].nz  +=   n.z;
			
			vertices[*(i_faces+2)].nx  +=   n.x;
            vertices[*(i_faces+2)].ny  +=   n.y;
            vertices[*(i_faces+2)].nz  +=   n.z;
			//	cout<<"normals "<<n.x<<" "<<n.y<<" "<<n.z<<endl;
			
		}
		//cout<<"norm 000 "<<vertices[0].nx<<" "<<vertices[0].ny<<" "<<vertices[0].nz<<endl;
		if (normalize)
		for ( typename vector< vertex<T> >::iterator i_vert = vertices.begin(); i_vert != vertices.end() ; i_vert++)
		{
			float norm = sqrt( i_vert->nx* i_vert->nx +  i_vert->ny* i_vert->ny + i_vert->nz* i_vert->nz);
			i_vert->nx /= norm;
            i_vert->ny /= norm;
			i_vert->nz /= norm;
		}
        
        if (reverse_winding)
        {
            for ( typename vector< vertex<T> >::iterator i_vert = vertices.begin(); i_vert != vertices.end() ; i_vert++)
            {
                i_vert->nx *= -1;
                i_vert->ny *= -1;
                i_vert->nz *= -1;
            }
        }
		
    }
	
	template<class T, class T2>
	unsigned int fslSurface<T,T2>::getNumberOfDataTableEntries()
	{
		return dataTable.size();
	}
	
	template<class T, class T2>
	map<int,float4> fslSurface<T,T2>::getDataTable()
	{
		return dataTable;
	}
	
	
	
	template<class T, class T2>
	std::vector<int> fslSurface<T,T2>::getScalarIndices( const int & index )
	{
		//scalar_indices
		
		//cout<<"get scalar inds "<<index<<" "<<scalar_indices.size()<<endl;
		//	if (!(index < scalar_indices.size()))
		//	{return ;
		
		return scalar_indices[index];
	}
    
    template<class T, class T2>
    void fslSurface<T,T2>::subtractScalars(const unsigned int & index, const T & value )
    {
        vector<T> new_sc( N_vertices , 0 );
        typename vector<T>::iterator i_new = new_sc.begin();
        for ( typename vector<T>::iterator i_sc = scalar_data[index].begin(); i_sc != scalar_data[index].end(); ++i_sc,++i_new)
        {
            *i_new = (*i_sc - value);
            
        }
        
        stringstream ss;
        ss<<value;
        string stemp;
        ss>>stemp;
        insertScalars(new_sc, index, scalar_names[index]+"_sub"+stemp);
        
    }
    
    template<class T, class T2>
    void fslSurface<T,T2>::multiplyScalars(const unsigned int & index, const T & value )
    {
        vector<T> new_sc( N_vertices , 0 );
        typename vector<T>::iterator i_new = new_sc.begin();
        for ( typename vector<T>::iterator i_sc = scalar_data[index].begin(); i_sc != scalar_data[index].end(); ++i_sc,++i_new)
        {
            *i_new = (*i_sc * value);
            
        }
        
        stringstream ss;
        ss<<value;
        string stemp;
        ss>>stemp;
        insertScalars(new_sc, index, scalar_names[index]+"_mul"+stemp);
        
    }
    
    template<class T, class T2>
    void fslSurface<T,T2>::binariseScalars(const unsigned int & index, const T & threshold )
    {
        vector<T> new_sc( N_vertices , 0 );
        typename vector<T>::iterator i_new = new_sc.begin();
        for ( typename vector<T>::iterator i_sc = scalar_data[index].begin(); i_sc != scalar_data[index].end(); ++i_sc,++i_new)
        {
            *i_new = (*i_sc > threshold) ? 1 : 0 ; 
            
        }
        insertScalars(new_sc, index, scalar_names[index]+"_bin");
        
    }

    
    template<class T, class T2>
    void fslSurface<T,T2>::thresholdScalars(const unsigned int & index, const T & threshold )
    {
        vector<T> new_sc( N_vertices , 0 );
        typename vector<T>::iterator i_new = new_sc.begin();
        for ( typename vector<T>::iterator i_sc = scalar_data[index].begin(); i_sc != scalar_data[index].end(); ++i_sc,++i_new)
        {
            if (*i_sc > threshold) *i_new = *i_sc;
        }
        insertScalars(new_sc, index, scalar_names[index]+"_thr");
    }
    
    template<class T, class T2>
    void fslSurface<T,T2>::upperThresholdScalars(const unsigned int & index, const T & threshold )
    {
        vector<T> new_sc( N_vertices , 0 );
        typename vector<T>::iterator i_new = new_sc.begin();
        for ( typename vector<T>::iterator i_sc = scalar_data[index].begin(); i_sc != scalar_data[index].end(); ++i_sc,++i_new)
        {
            if (*i_sc < threshold) *i_new = *i_sc;
        }
        insertScalars(new_sc, index, scalar_names[index]+"_uthr");
    }

    
    template<class T, class T2>
    void fslSurface<T,T2>::insertScalars(const std::vector<T> & scalars,  const unsigned int & index, const std::string & sname )
    {
        scalar_data.insert(scalar_data.begin()+index, scalars);
		scalar_names.insert(scalar_names.begin() + index, sname);
        typename vector<T>::const_iterator i_sc = scalars.begin();
		for ( typename vector< vertex<T> >::iterator i_v = vertices.begin(); i_v != vertices.end(); ++i_sc, ++i_v)
		{
			i_v->sc = *i_sc;
		}
    }
    
    template<class T, class T2>
    void fslSurface<T,T2>::insertNonVertScalars( const std::vector<int> & scalars, const unsigned int & index, const std::string & sname )
    {
        nonvert_int_sc_data.insert(nonvert_int_sc_data.begin()+index, scalars);
		nonvert_int_sc_data_names.insert(nonvert_int_sc_data_names.begin() + index, sname);
 
    }
    
    template<class T, class T2>
    void fslSurface<T,T2>::insertNonVertScalars( const std::vector<float> & scalars, const unsigned int & index, const std::string & sname )
    {
        nonvert_float_sc_data.insert(nonvert_float_sc_data.begin()+index, scalars);
		nonvert_float_sc_data_names.insert(nonvert_float_sc_data_names.begin() + index, sname);

    }

    template<class T, class T2>
    void fslSurface<T,T2>::printScalars(const unsigned int & index )
    {
        
		for ( typename vector<T>::const_iterator i_sc = scalar_data[index].begin(); i_sc != scalar_data[index].end(); ++i_sc, ++i_sc)
		{
			cout<<*i_sc<<" "<<endl;
		}
    }


	template<class T, class T2>
	void fslSurface<T,T2>::addScalars(const vector<T> & scalars, const std::string & sname )
	{
		scalar_data.push_back(scalars);
		scalar_names.push_back(sname);
		
		 typename vector<T>::const_iterator i_sc = scalars.begin();
		for ( typename vector< vertex<T> >::iterator i_v = vertices.begin(); i_v != vertices.end(); ++i_sc, ++i_v)
		{
			i_v->sc = *i_sc;
		}
		
	}
    
    template<class T, class T2>
    void fslSurface<T,T2>::addVectors(const vector<T> & vec, const string & sname, const unsigned int & index ){
        
        vector< vector<float> >::iterator i = vector_data.begin() +index;
        vector<string>::iterator i_n = vector_names.begin() + index;
        vector_data.insert(i, vec );
        vector_names.insert(i_n, sname );
    }
    


    template<class T, class T2>
    std::string fslSurface<T,T2>::getScalarName( const unsigned int & sc_ind )
    {
        return scalar_names[sc_ind];
    }

    
    template<class T, class T2>
    T fslSurface<T,T2>::getScalar( const unsigned int & sc_ind, const unsigned int & vert_ind ) const
    { 
        return scalar_data[sc_ind][vert_ind]; 
    }

	template<class T, class T2>
    std::vector<T> fslSurface<T,T2>::getScalars(const unsigned int & index ) const
    {
        return scalar_data[index];
    }
    
	template<class T, class T2>
    std::vector<T> fslSurface<T,T2>::getScalars(const string & name ) const
    {
        unsigned int index=0;
        for ( vector<string>::const_iterator i_name = scalar_names.begin(); i_name != scalar_names.end(); ++i_name,++index)
        {
            if ( *i_name == name )
                return scalar_data[index];
        }
        vector<T> empty;
        return empty;
            
    }
    template<class T, class T2>
    std::vector<int> fslSurface<T,T2>::getNonVertIntScalars(const unsigned int & index ) const
    {
        return nonvert_int_sc_data[index]; 
    }
    template<class T, class T2>
    std::vector<int> fslSurface<T,T2>::getNonVertIntScalars(const std::string & name ) const
    {
        cout<<"get non vert scalar "<<nonvert_int_sc_data.size()<<" "<<nonvert_int_sc_data_names.size()<<endl;
        unsigned int index=0;
        for ( vector<string>::const_iterator i_name = nonvert_int_sc_data_names.begin(); i_name != nonvert_int_sc_data_names.end(); ++i_name,++index)
        {
            cout<<"name "<<*i_name<<" "<<name<<endl;
            if ( *i_name == name )
                return nonvert_int_sc_data[index];
        }
        vector<int> empty;
        return empty;
    }
    template<class T, class T2>
    std::vector<float> fslSurface<T,T2>::getNonVertFloatScalars(const unsigned int & index ) const
    {
        return nonvert_float_sc_data[index]; 
    }
    template<class T, class T2>
    std::vector<float> fslSurface<T,T2>::getNonVertFloatScalars(const std::string & name ) const
    {
        unsigned int index=0;
        for ( vector<string>::const_iterator i_name = nonvert_int_sc_data_names.begin(); i_name != nonvert_int_sc_data_names.end(); ++i_name,++index)
        {
            if ( *i_name == name )
                return nonvert_float_sc_data[index];
        }
        vector<float> empty;
        return empty;
    }
    
    template<class T, class T2>
    void fslSurface<T,T2>::replaceScalars(const std::vector< std::vector<T> > & scalars_in, const vector<string> & sc_names_in )
    {
//        scalar_data.clear();
        scalar_data = scalars_in;
        scalar_names = sc_names_in;
    }
    
    template<class T, class T2>
    void fslSurface<T,T2>::clearScalars()
    {
        scalar_data.clear();
        scalar_names.clear();
    }

	template<class T, class T2>
	void fslSurface<T,T2>::setScalar(const unsigned int & sc_ind, const unsigned int & vert_ind, const T& value )
	{ // Note this will generate exceptions if scalar_data is accessed outside of bounds
	  scalar_data[sc_ind][vert_ind]=value; 
	}

	template<class T, class T2>
	void fslSurface<T,T2>::reinitialiseScalars(const unsigned int fields)
	{ //assumes data has same number of cols as surface has vertices - use at your peril!!
	  clearScalars();
	  scalar_data.resize(fields); //number of scalar fields
	  scalar_names.resize(fields);
	  for(unsigned int field=0; field < getNumberOfScalarData(); field++)
	    scalar_data[field].resize( getNumberOfVertices() );
	}


	template<class T, class T2>
	void fslSurface<T,T2>::setScalars(const int & index )
	{
		//cout<<"set scalars1 "<<scalar_data.size()<<" "<<scalar_indices.size()<<" "<<index<<endl;
		//	cout<<scalar_indices[index].size()<<endl;
		//	cout<<"set Scalars "<<scalar_indices[index][0]<<" "<<index<<" "<<scalar_data.size()<<" "<<scalar_indices.size()<<" "<<scalar_indices[index].size()<<" "<<scalar_data[scalar_indices[index][0]].size()<<endl;
		//if (index >= static_cast<int>(scalar_data.size()))
		//	return;
		//	if (scalar_indices[index].empty()) {
		//		cout<<"empty"<<endl;
		//			return;
		//		}
		if (scalar_data[index].size() == vertices.size())
		{	
			typename vector< vertex<T> >::iterator i_v = vertices.begin();
			typename vector<T>::iterator i_end = scalar_data[index].end();
			//		cout<<"start scalars "<<endl;
			int count=0;
			for ( typename vector<T>::iterator i = scalar_data[index].begin(); i != i_end; ++i, ++i_v,++count)
			{
				//if (count==0)
					//cout<<"scalars "<<*i<<endl;
				i_v->sc = *i;
			}
		}
		//cout<<"done scalars "<<endl;
	}
		
	
	template class fslSurface<float,int>;
	template class fslSurface<float,unsigned int>;
	
}


