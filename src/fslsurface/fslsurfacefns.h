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
#ifndef FSLSURFACEFNS_H
#define FSLSURFACEFNS_H

#include <map>
#include <list>
#include <vector>
#include <fslsurface_structs.h>
#include <fslsurface.h>
#include <newimage/newimageall.h>
//template<class T, class T2>
//class fslSurface;

namespace fslsurface_name {
	
    //-----------------------float3 functions------------------------//

    //some function for float3 
   // void normalize( float3 & v1);
    float dot(const vec3<float> & v1, const vec3<float> & v2 );
    vec3<float> cross(const vec3<float> & v1, const vec3<float> & v2 );

    vec3<float> subtract(const vec3<float> & v1, const vec3<float> & v2);
    //------------------------SURFACE XFMS----------------------------//
		template<class T,class T2,class T3>
	void apply_xfm( fslSurface<T,T2>& surf, const std::vector<T3> & xfm );
	
    template<class T,class T2>
	void apply_xfm( fslSurface<T,T2>& surf, const std::string & xfm );
    
    template<class T,class T2>
	void apply_warp( fslSurface<T,T2>& surf, const NEWIMAGE::volume4D<float> & warp );
    
	template <class T, class T2>
	void deformSurface( fslSurface<T,T2>& surf , const unsigned int & N, const T & w_norm, const T & w_tan  );
	
    //--------------------SCALAR PROCESSING----------------//
    template <class T2>
    void remove_connection( std::vector< std::list<T2> > & connections, const unsigned int & index );
    template <class T2>
    void remove_item( std::list<T2> & connections, const unsigned int & index );
    
    //computes connections matrix given 2 connection matrices, con1 +con2 gives con 3
    template<class T2>
    std::vector< std::list<T2> > computeNextLevelConnections(  std::vector< std::list<T2> > & connections0, std::vector< std::list<T2> > & connections1 );
    template<class T2>
     std::vector< std::list< std::pair<T2,T2> > > computeNextLevelConnectionsWithMultiWayPoints(  std::vector< std::list<T2> > & connections0, std::vector< std::list<T2> > & connections1 );
    
    template <class T, class T2>
   // std::list< std::pair<T2,T> > 
    void conn2dist(fslSurface<T,T2>& surf, const T2 & vert_index, std::list< std::pair<T2,T2> > & conns, std::vector< std::list< std::pair<T2,T> > > & index_distances );
    
    template <class T,class T2>
    void cluster( fslSurface<T,T2>& surf , const unsigned int & index, const T & threshold );
    
    template <class T, class T2>
    void sc_smooth_mean_neighbour( fslSurface<T,T2>& surf , const unsigned int & index );
    
    template <class T, class T2>
    void sc_smooth_gaussian_geodesic( fslSurface<T,T2>& surf , const unsigned int & index, const T & sigma, const T & extent,bool run4D=false );
    
    template <class T, class T2>
    void multiplyScalarsByScalars(fslSurface<T,T2>& surf , const unsigned int & index, const fslSurface<T,T2>& surf2 , const unsigned int & index1);

    template <class T, class T2>
    void subtractScalarsFromScalars(fslSurface<T,T2>& surf , const unsigned int & index, const fslSurface<T,T2>& surf2 , const unsigned int & index1, std::string name="sc_sub_sc");//surf-surf2
    
    //----------------VECTOR PROCESSING------------------//
    template <class T, class T2>
    void projectVectorsOntoNormals(fslSurface<T,T2>& surf ,const unsigned int & index);

    //-----------------MARCHING CUBES---------------------//
		template<class T,class T2, class T3>
		void runMarchingCubesOnAllLabels( fslSurface<T,T2>& surf , const T3* image,  const fslsurface_name::image_dims & dims, const T3 & thresh_in);

		template<class T,class T2, class T3>
		void marchingCubes( fslSurface<T,T2>& surf , const T3* images,  const fslsurface_name::image_dims & dims, const T3 & thresh, const T & label, const fslsurface_name::MarchingCubesMode & mode);
		
    
    //--------------------DATA MANIPULATIONS---------------//
		template<class T,class T2>
		void maskScalars( fslSurface<T,T2>& surf , const unsigned int & sc_index, const short* mask,fslsurface_name::image_dims dims, bool  setToCurrentSc = true );
		template<class T,class T2>
		void maskSurface( fslSurface<T,T2>& surf , const short* mask, const fslsurface_name::image_dims dims, bool  setToCurrentSc = true );
		
		template<class T,class T2, class T3>
		void mapScalars( fslSurface<T,T2>& surf , const unsigned int & sc_index, const std::map<T3,T3> & scMap, bool setToCurrentSc =true ); 		
		//void maskScalars(const unsigned int & sc_index, const short* mask ,const fslsurface_name::image_dims & dims, bool setToCurrentSc =true); 		
		
		template< class T >
		void vertexInterp( const float & thresh, const vec3<T> & v1 , const float  & val1, const vec3<T> & v2, const float & val2 , vertex<T> & vout);
    
    
        template<class T,class T2, class T3>
    void insertScalars(fslsurface_name::fslSurface<T,T2>& surf , const NEWIMAGE::volume4D<T> & time_series, const std::string & time_name );
    
		template<class T,class T2>
    void sampleTimeSeries(fslsurface_name::fslSurface<T,T2>& surf , const NEWIMAGE::volume4D<T> & time_series, const std::string & time_name );
    
    template<class T,class T2>
    void writeTimeSeries( const fslsurface_name::fslSurface<T,T2>& surf , NEWIMAGE::volume4D<T> & time_series , const std::string & outname);
    
}

#endif
