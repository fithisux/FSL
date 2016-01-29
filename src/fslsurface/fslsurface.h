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
#ifndef _fslSurface
#define _fslSurface

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include <map>
#include <list>
#include "fslsurface_structs.h"
//extern "C" {	
//#include <gifti/gifti_io.h>
//}

namespace fslsurface_name{
	//using namespace mesh;
	
	class fslSurfaceException : public std::exception{
		
	public:
		const char* errmesg;
		fslSurfaceException(const char* msg)
		{
		  errmesg=msg;
		}
		virtual const char* what() const throw()
		{
			std::cout<<errmesg<<std::endl;
			return errmesg;
		}
//	private:
		
	};
	
//---------------------------float3-----------------------------//
	template<class T>
	inline vec3<T> subtract( const vec3<T> & a, const vec3<T> & b );
	template<class T>
	inline vec3<T> normal( const vec3<T> & a, const vec3<T> & b );
	//---------------------------end float3-----------------------------//
	
	
//	extern "C" {
//
//		void giftWrapper_read_Image(const std::string & filename, gifti_image* gii_surf);
//	
//	}
//	class fslSurfaceIO;
//	class fslSurfaceFNS;
//	enum MarchingCubesMode {EQUAL_TO, GREATER_THAN, GREATER_THAN_EQUAL_TO};

	
	template<class T,class T2>
	class fslSurface {
		
		
	public:
	
		
		enum DataType { ascii, littleEndian, bigEndian };
		enum NumFormat { float32, float64 };
	
        //-------operators------//
        void operator+=(const fslSurface<T,T2> &surf);
        void operator-=(const fslSurface &surf);
        void operator*=(const T &mul);
        void operator/=(const T &div);
        
        fslSurface<T,T2> operator +(const  fslSurface<T,T2> & surf );
        fslSurface<T,T2> operator -(const  fslSurface<T,T2> & surf );
        fslSurface<T,T2> operator *(const  T & surf );
        fslSurface<T,T2> operator /(const  T & surf );

      //  void operator+=(const fslSurface &surf);
      //  void operator+=(const fslSurface &surf);

        //-----------------------//
        
        
		fslSurface();
        fslSurface( const fslSurface<T,T2> & surf );

		//fslSurface(const std::string & filename);
		virtual ~fslSurface();
		//vertex<T> vbegin3();
        //vertex access
        std::vector<T> getBounds() const;
        //-----------FUNCTIONS TO ACCESS USING ITERATORS, will be faster for multiple access----------------//

        typename std::vector< fslsurface_name::vertex<T> >::iterator vbegin();
        typename std::vector< fslsurface_name::vertex<T> >::const_iterator const_vbegin() const;
	  typename std::vector< fslsurface_name::vertex<T> >::iterator vend() ;
        typename std::vector< fslsurface_name::vertex<T> >::const_iterator const_vend() const;
        
        //face access 
        typename std::vector< T2 >::const_iterator const_facebegin() const;
        typename std::vector< T2 >::const_iterator const_faceend() const;
        
        //assumes triangles 
        //get specific face as a vec3 (x,y,z are vertex indices)
        typename fslsurface_name::vec3<T2> getFace( const unsigned int & index, int face_size = 3 ) const;
        typename std::vector<T2> getFaces() const;
        std::vector< fslsurface_name::vec3<T> > getFaceVertices( const unsigned int & index, int face_size =3 ) const;

        std::vector< T > getFaceVerticesUnwrapped( const unsigned int & index, int face_size =3 ) const;

        
        
        typename std::vector< T >::const_iterator const_scbegin(const unsigned int & index ) const;
        typename std::vector< T >::const_iterator const_scend(const unsigned int & index) const;

        typename std::vector< T >::const_iterator const_vecbegin(const unsigned int & index ) const;
        typename std::vector< T >::const_iterator const_vecend(const unsigned int & index) const;
        //-----------------------------------------------------------------------------------------------//
        //-----------FUNCTIONS FOR Random access to data ACCESS, will be faster for multiple access----------------//

        //various attributes
		float3 getCOG();
		unsigned int getNumberOfFaces() const  { return faces.size()/3; }
		unsigned int getNumberOfVertices() const { return vertices.size(); }
		unsigned int getNumberOfScalarData() const { return scalar_data.size(); }
	  unsigned int getNumberOfVectorData() const { return vector_data.size(); }
        

	  std::string getAnatomicalName() const { return anatomical_name; }
	  std::string getAnatomicalName2() const { return anatomical_name2; }
	  std::string getGeometryName() const { return geometry_name; }
	  std::string getTopologyName() const { return topology_name; }
	  void setAnatomicalName(const std::string& aname) { anatomical_name = aname; }
	  void setAnatomicalName2(const std::string& aname) { anatomical_name2 = aname; }
	  void setGeometryName(const std::string& gname) { geometry_name = gname; }
	  void setTopologyName(const std::string& tname) { topology_name = tname; }


        //access to data (e.g. vertex, scalars, triangles, etc..)
        fslsurface_name::vec3<T> getVertexCoord(const unsigned int & index);
        vertex<T> getVertex( const unsigned int & vert_ind ) { return vertices[ vert_ind ]; }
        
        T getScalar( const unsigned int & sc_ind, const unsigned int & vert_ind ) const;
        std::string getScalarName( const unsigned int & sc_ind );

        std::vector<T> getScalars(const unsigned int & index ) const;
        std::vector<T> getScalars(const std::string & name ) const;
        std::vector<int> getNonVertIntScalars(const unsigned int & index ) const;
        std::vector<int> getNonVertIntScalars(const std::string & name ) const;
        std::vector<float> getNonVertFloatScalars(const unsigned int & index ) const;
        std::vector<float> getNonVertFloatScalars(const std::string & name ) const;

        
		std::vector<int> getScalarIndices( const int & index );

        
        //-----------------------------------------------------------------------------------------------//
        //---------------------------------FUNCTIONS TO SET, ADD OR INSERT DATA-------------------------------------------------------------//
        void setVertices( const std::vector<T> &  verts );
        void setVert(const unsigned int & index, const fslsurface_name::vec3<float> & vert_coord);

		void setFaces( const std::vector<T2> &  faces_in );
		void setFaces( const std::vector< std::vector<T2> > &  faces_in );
        
        
		void setScalars(const int & index );
		void reinitialiseScalars(const unsigned int fields);
		void setScalar( const unsigned int & sc_ind, const unsigned int & vert_ind, const T& value ); 

        void replaceScalars(const std::vector< std::vector<T> > & scalars , const std::vector< std::string > & sc_names);
        void clearScalars();
		void insertScalars(const std::vector<T> & scalars, const unsigned int & index, const std::string & sname);
		void addScalars(const std::vector<T> & scalars, const std::string & sname );
        void insertNonVertScalars( const std::vector<int> & scalars, const unsigned int & index, const std::string & sname );
        void insertNonVertScalars( const std::vector<float> & scalars, const unsigned int & index, const std::string & sname );
        void printScalars(const unsigned int & index );
        
		void addVectors(const std::vector<T> & vec, const std::string & sname, const unsigned int & index );
        
        

        //-----------------------------------------------------------------------------------------------//

        //---------------------------------SCALAR MANIPULATOR FUNCTIONS--------------------------------------------------------------//
        
        void subtractScalars(const unsigned int & index, const T & value );
        void multiplyScalars(const unsigned int & index, const T & value );
        
        void binariseScalars(const unsigned int & index, const T & threshold );
        void thresholdScalars(const unsigned int & index, const T & threshold );
        void upperThresholdScalars(const unsigned int & index, const T & threshold );
        
        //-----------------------------------------------------------------------------------------------//


        void copyNormalsToVectors(const unsigned int & index);
        void copyVerticesToVectors();
        void copyVerticesToVectors( const std::string & name );
        void copyVertices( const fslSurface<T,T2> & surf );

        
        //-------------------vertex operators---------------------------------------//
        float L2norm( const unsigned int & index0, const unsigned int & index1);
        fslsurface_name::vec3<float> subtractVerts(const unsigned int & index0, const unsigned int & index1);

        
        //------------------------------------------coordinate system functions ------------------------------------------//
        unsigned int getNumberOfCoordinateSystems();
		void addCoordSystem(const std::vector<double> & xfm, const std::string & xfm_space );
		void addCoordSystem(const std::vector<float> & xfm, const std::string & xfm_space );
        void copyCoordSystems(const fslSurface<T,T2> & surf);
        void clearCoordSystems();
        std::vector<double> getCoordinateSystem( const std::string & name);

		
        
        
		unsigned int getNumberOfDataTableEntries();
		std::map<int,float4> getDataTable();
		
		void computeAdjacentTriangles();
		void computeAdjacentVertices(const bool & bidirectional);

		void calculateNormals(bool normalize = true, bool reverse_winding=false);

   
		//------------------------------------------------------------//
		static void copy_coordsystem(fslSurface & surf_dest, const fslSurface & surf_src);
		
		
		////////-------------FRIEND FUNCTIONS FSLSURFACE_IO-------------/////////
		template<class U,class U2>
		friend int read_surface( fslSurface<U,U2>& surf, const std::string & filename);
		template<class U,class U2>
		friend int readGIFTI( fslSurface<U,U2>& surf, const std::string & filename );
		template<class U,class U2>
		friend int readVTK( fslSurface<U,U2>& surf, const std::string & filename );
		template<class U,class U2>
		friend int readPLY( fslSurface<U,U2>& surf, const std::string & filename );
		template<class U,class U2>
		friend int writeGIFTI(  const fslSurface<U,U2>& surf, const std::string & filename, int enc);
		////////-------------FRIEND FUNCTIONS FSLSURFACE_FNS-------------/////////
		template<class U,class U2, class U3>
		 friend void runMarchingCubesOnAllLabels( fslSurface<U,U2>& surf , const U3* image,  const fslsurface_name::image_dims & dims, const U3 & thresh_in);
		template<class U,class U2, class U3>
		 friend void marchingCubes( fslSurface<U,U2>& surf , const U3* images,  const fslsurface_name::image_dims & dims, const U3 & thresh, const U & label, const MarchingCubesMode & mode);
		template<class U,class U2>
		friend void maskScalars( fslSurface<U,U2>& surf , const unsigned int & sc_index, const short* mask,fslsurface_name::image_dims dims, bool  setToCurrentSc );
		template<class U, class U2>
		friend void maskSurface( fslSurface<U,U2>& surf , const short* mask, const fslsurface_name::image_dims  dims, bool  setToCurrentSc );
		template<class U,class U2, class U3>
		friend void mapScalars( fslSurface<U,U2>& surf , const unsigned int & sc_index, const std::map<U3,U3> & scMap, bool setToCurrentSc  ); 		
		template<class U,class U2,class U3>
		friend void apply_xfm( fslSurface<U,U2>& surf, const std::vector<U3> & xfm );
    
        template <class U,class U2>
        friend void cluster( fslSurface<U,U2>& surf , const unsigned int & index, const U & threshold );
        
        template <class U,class U2>
        friend void sc_smooth_mean_neighbour( fslSurface<U,U2>& surf , const unsigned int & index );

        template <class U,class U2>
        friend void sc_smooth_gaussian_geodesic( fslSurface<U,U2>& surf ,const unsigned int & index, const U & st_dev , const U & extent , bool run4D);
        
        template<class U, class U2>
        friend void multiplyScalarsByScalars(fslSurface<U,U2>& surf , const unsigned int & index, const fslSurface<U,U2>& surf2 , const unsigned int & index1);
        
        template<class U, class U2>
        friend void subtractScalarsFromScalars(fslSurface<U,U2>& surf , const unsigned int & index, const fslSurface<U,U2>& surf2 , const unsigned int & index1, std::string name );

        
        template<class U,class U2>
        friend void projectVectorsOntoNormals(fslSurface<U,U2>& surf ,const unsigned int & index);
        
		////////-------------FRIEND FUNCTIONS FSLSURFACE_GL-------------/////////
		
		template<class U,class U2>
		friend void glBufferData_Vertices( const fslSurface<U, U2>& surf, const GLuint & vbo );  
		template<class U,class U2>
		friend void glBufferSubData_Vertices( const fslSurface<U, U2>& surf, const GLuint & vbo )  ;
		template<class U,class U2>
		friend void glBufferData_Faces( const fslSurface<U, U2>& surf, const GLuint & vbo )  ;
		template<class U,class U2>
		friend void glBufferSubData_Faces( const fslSurface<U, U2>& surf, const GLuint & vbo )  ;
		template<class U,class U2>
		friend void render( const fslSurface<U, U2>& surf, const GLint & vertexLoc , const GLint &  normalLoc, const GLint & scalarLoc, const GLuint & vbo_verts, const GLuint & vbo_ele_array ) ;
		
        template<class U,class U2>
        friend void depthSortTriangles( const fslSurface<U,U2>& surf, const GLuint & vbo_verts, const GLuint & vbo_tris);

		//friend class fslSurfaceIO;
		//friend class fslSurfaceFNS;x

	protected:
		
		
	private:


		//DataType machineEndianness();
		
	  std::vector< vertex<T> > vertices;
	  std::vector<T2> faces;//assumes triangles
	  std::vector< vec3<T> > tangents;
	  
	  std::vector< std::list<T2> > adj_tris;
	  std::vector< std::list<T2> > adj_verts;
	  bool adj_verts_bi;
	  std::vector< std::vector<double> > v_coord_sys; 
	  std::string csys_dspace;
	  std::vector< std::string > v_csys_xfmspace; 
	  
	  std::string topology_name;
	  std::string geometry_name;
	  std::string anatomical_name;
	  std::string anatomical_name2;


	  //if intending to use the openGL use unsigned int (suggested anyways)
	  std::vector< std::vector<T> > vector_data;
	  std::vector< std::string > vector_names;

	  std::vector< std::vector<T> > scalar_data;
	  std::vector< std::string > scalar_names;
	  std::vector< std::vector<int> > scalar_indices;
	  
	  std::map<int,float4> dataTable;
		
	  std::vector< std::vector<float> > nonvert_float_sc_data;
	  std::vector< std::string > nonvert_float_sc_data_names;
	  
	  std::vector< std::vector<int> > nonvert_int_sc_data;
	  std::vector< std::string > nonvert_int_sc_data_names;

//		std::vector<std::string> element_names;	
//		std::vector< std::vector<std::string> > element_types;//assuming same type
//		std::vector< std::vector<std::string> > element_prop_names;	
//		std::vector< int > element_num_props;	
//		std::vector< std::pair<std::string,std::string> > list_types;


		std::vector<unsigned int> element_sizes;
		
		NumFormat num_format;
		DataType format;
		//FileType file_type;
		std::string ply_version;
		
		unsigned int N_vertices,N_triangles;
		float3 cog;
		std::map<int,std::string> datatype_str;

		
		//		vector<string> celement_names;
		//	DataType format;
//		int format;
		//EndianType endian_format;

	};
}
#endif
