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
#include <fslsurfacefns.h>
#include <fslsurface.h>
#include <sstream>
//STL includes 
#include <set>
#include <queue>
#include<cmath>
using namespace std;

#define PI 3.141592653589793238462643383

//using namespace fslsurface_name;
namespace fslsurface_name {
	/*
    void normalize( float3 & v1)
    {
        float mag=sqrtf(v1.x*v1.x + v1.y*v1.y + v1.z*v1.z);
        v1.x/=mag;
        v1.y/=mag;
        v1.z/=mag;
    }
    */
    float dot(const vec3<float> & v1, const vec3<float> & v2 )
    {
        return (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
        
    }
    
    vec3<float> cross(const vec3<float> & v1, const vec3<float> & v2 )
    {
        return vec3<float>( v1.x*v2.y - v1.y*v2.x , \
                       v1.y*v2.z - v1.z*v2.y , \
                       v1.z*v2.x - v1.x*v2.z);
    }

    
    vec3<float> subtract(const vec3<float> & v1, const vec3<float> & v2)
    {
        return vec3<float>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z );
    }

    template<class T,class T2>
    void projectVectorsOntoNormals(fslSurface<T,T2>& surf ,const unsigned int & index)
    {
        //caculate normals and normalize
        surf.calculateNormals(true,false);
        vector<T> new_vecs(surf.N_vertices*3);
        vector<T> new_scalars(surf.N_vertices);

        string name= surf.vector_names[index];
        surf.vector_names.push_back(name + "_proj_onto_normals");
        typename vector<T>::iterator i_vec = surf.vector_data[index].begin();
        typename vector<T>::iterator i_vec_new = new_vecs.begin();
        typename vector<T>::iterator i_sc_new = new_scalars.begin();

        for ( typename vector< vertex<T> >::iterator i_v = surf.vbegin(); i_v != surf.vend(); ++i_v, i_vec+=3, i_vec_new+=3,++i_sc_new)
        {
            T mag = i_v->nx * (*i_vec) + i_v->ny * (*(i_vec+1)) + i_v->nz * (*(i_vec+2));
            *i_sc_new = mag;
            *i_vec_new = i_v->nx * mag;
            *(i_vec_new+1) = i_v->ny * mag;
            *(i_vec_new+2) = i_v->nz * mag;

        }
        
        surf.scalar_data.insert(surf.scalar_data.begin(),new_scalars);
        surf.scalar_names.insert(surf.scalar_names.begin(), name + "_proj_onto_normals");


        
        surf.vector_data.insert(surf.vector_data.begin(),new_vecs);
        surf.vector_names.insert(surf.vector_names.begin(), name + "_proj_onto_normals");

    }
    template void projectVectorsOntoNormals<float,unsigned int>(fslSurface<float,unsigned int>& surf ,const unsigned int & index);

    
	template<class T,class T2>
	void apply_xfm( fslSurface<T,T2>& surf, const string & xfm ){
        vector<T> vxfm(16,0);
        ifstream fmat(xfm.c_str());
        T temp;
        for (typename vector<T>::iterator i_v = vxfm.begin(); i_v != vxfm.end(); ++i_v)
        {
            fmat>>temp;
           *i_v = temp;
        }
        fmat.close();
        apply_xfm<T,T2,T>(surf,vxfm );
    }
    
    
	template void apply_xfm<float, unsigned int>( fslSurface<float, unsigned int>& surf, const string & xfm );
	template void apply_xfm<float, int>( fslSurface<float, int>& surf, const string & xfm );
	template void apply_xfm<double, unsigned int>( fslSurface<double, unsigned int>& surf, const string & xfm );
	template void apply_xfm<double, int>( fslSurface<double, int>& surf, const string & xfm );
    

    
    
	template<class T,class T2,class T3>
	void apply_xfm( fslSurface<T,T2>& surf, const vector<T3> & xfm )
	{
		//assume bottom row of xfm is 0 0 0 1
		for ( typename vector< vertex<T> >::iterator i_vert = surf.vertices.begin(); i_vert != surf.vertices.end(); ++i_vert )
		{
			T x_in = i_vert->x;
			T y_in = i_vert->y;
			T z_in = i_vert->z;
			
			i_vert->x = xfm[0] * x_in + xfm[1] * y_in + xfm[2] * z_in + xfm[3] ;
			i_vert->y = xfm[4] * x_in + xfm[5] * y_in + xfm[6] * z_in + xfm[7] ;
			i_vert->z = xfm[8] * x_in + xfm[9] * y_in + xfm[10] * z_in + xfm[11] ;

		}
	}
	template void apply_xfm<float, unsigned int, float>( fslSurface<float, unsigned int>& surf, const vector<float> & xfm );
	template void apply_xfm<float, unsigned int, double>( fslSurface<float, unsigned int>& surf, const vector<double> & xfm );
	template void apply_xfm<float, int, double>( fslSurface<float, int>& surf, const vector<double> & xfm );
	template void apply_xfm<double, unsigned int, float>( fslSurface<double, unsigned int>& surf, const vector<float> & xfm );
	template void apply_xfm<double, int, double>( fslSurface<double, int>& surf, const vector<double> & xfm );

    
    template<class T,class T2>
	void apply_warp( fslSurface<T,T2>& surf, const NEWIMAGE::volume4D<float> & warp )
    {
        vec3<float> dims= vec3<float>(warp.xdim(),warp.ydim(),warp.zdim()); 
        for ( typename vector< vertex<T> >::iterator i_v=surf.vbegin(); i_v != surf.vend(); ++i_v)
        {
//         cout<<i_v->x<<" "<<i_v->y<<" "<<i_v->z<<" "<<
            vec3<float> vox_coord = vec3<float>(i_v->x / dims.x,i_v->y / dims.y, i_v->z / dims.z );
            i_v->x += warp[0].interpolate(vox_coord.x,vox_coord.y,vox_coord.z);
            i_v->y += warp[1].interpolate(vox_coord.x,vox_coord.y,vox_coord.z);
            i_v->z += warp[2].interpolate(vox_coord.x,vox_coord.y,vox_coord.z);
        }
        
    }
    template void apply_warp<float, unsigned int>( fslSurface<float, unsigned int>& surf, const NEWIMAGE::volume4D<float> & warp );

	
	
	template <class T, class T2>
	void deformSurface( fslSurface<T,T2>& surf , const unsigned int & N, const T & w_norm, const T & w_tan  )
	{
		surf.calculateNormals();
		if (surf.adj_tris.empty())
			surf.computeAdjacentTriangles();
		if (surf.adj_verts.empty())
			surf.computeAdjacentVertices(true);
		for ( unsigned int iter = 0 ; iter < N ; ++iter )
		{		
			for ( typename vector< vertex<T> >::iterator i_vert =  surf.vertices.begin() ; i_vert != surf.vertices.end(); ++i_vert )
			{
				i_vert->x += i_vert->nx * w_norm;
				i_vert->y += i_vert->ny * w_norm;
				i_vert->z += i_vert->nz * w_norm;

			}

			
			
		}
	}
  
    template <class T>
    void remove_connection( vector< list<T> > & connections, const unsigned int & index )
    {
        //remove backward connections
        for (typename list<T>::iterator i = connections[index].begin() ; i != connections[index].end(); ++i)
            for (typename list<T>::iterator ii = connections[*i].begin() ; ii != connections[*i].end(); ++ii)
                if (*i == index)
                {
                    connections[index].erase(ii);
                    break;
                }
        //clear row (forward conections)
        connections[index].clear();
    }
    
    template <class T>
    void remove_item( list<T> & connections, const unsigned int & index )
    {
        //remove backward connections
        for (typename list<T>::iterator i = connections.begin() ; i != connections.end(); ++i)
                if (*i == index)
                {
                    connections.erase(i);
                    break;
                }
    }
    
    //this computes the furhter down connections
    template<class T2>
    vector< list<T2> > computeNextLevelConnections( vector< std::list<T2> > & connections0, vector< std::list<T2> > & connections1 )
    {
        vector< list<T2> > new_cons(connections0.size());
        unsigned int curr_ind = 0 ;
        typename vector< list<T2> >::iterator i_new = new_cons.begin();
        for ( typename vector< list<T2> >::iterator i0 = connections0.begin(); i0 != connections0.end(); ++i0,++i_new,++curr_ind)
            for ( typename list<T2>::iterator ii0 = i0->begin(); ii0 != i0->end(); ++ii0)
            {
                for ( typename list<T2>::iterator ii1 = connections1[*ii0].begin(); ii1 != connections1[*ii0].end(); ++ii1)
                {
                    if (curr_ind != *ii1)
                    {
                        //check too make sure it does not map back onto a connection in the first matrix
                        bool exists=false;
                        for ( typename list<T2>::iterator iii0 = i0->begin(); iii0 != i0->end(); ++iii0)
                            if(*iii0 == *ii1)
                                exists=true;
                        //make sure it does already exist in two double connection
                        for ( typename list<T2>::iterator ii_new = i_new->begin(); ii_new != i_new->end(); ++ii_new)
                            if(*ii_new == *ii1)
                                exists=true;
                        
                        
                        if (!exists)
                            i_new->push_back(*ii1);
                    }
                    }
            }

 //       cout<<"new cons"<<endl;
   //     for ( typename vector< list<T2> >::iterator i0 = new_cons.begin(); i0 != new_cons.end(); ++i0)
     //   {
            
       //         for ( typename list<T2>::iterator ii0 = i0->begin(); ii0 != i0->end(); ++ii0)
         //           cout<<*ii0<<" ";
           // cout<<endl;
       //     }
        //cout<<"end new cons"<<endl;
        
        return new_cons;
    }
    template   vector< list<unsigned int> > computeNextLevelConnections( vector< std::list<unsigned int> > & connections0, vector< std::list<unsigned int> > & connections1 );
    
    //this computes the furhter down connections
    template<class T2>
     vector< list< pair<T2,T2> > > computeNextLevelConnectionsWithMultiWayPoints( vector< std::list<T2> > & connections0, vector< std::list<T2> > & connections1 )
    {
        //second in pair 
        vector< list< pair<T2,T2> > > new_cons(connections0.size());
       // vector< list< pair<T2,T2> > > waypoints();
        unsigned int curr_ind = 0 ;
        typename  vector< list< pair<T2,T2> > >::iterator i_new = new_cons.begin();
        for ( typename vector< list<T2> >::iterator i0 = connections0.begin(); i0 != connections0.end(); ++i0,++i_new,++curr_ind)
            for ( typename list<T2>::iterator ii0 = i0->begin(); ii0 != i0->end(); ++ii0)
            {
                for ( typename list<T2>::iterator ii1 = connections1[*ii0].begin(); ii1 != connections1[*ii0].end(); ++ii1)
                {
                    if (curr_ind != *ii1)
                    {
                        //check too make sure it does not map back onto a connection in the first matrix
                        bool exists=false;
                        for ( typename list<T2>::iterator iii0 = i0->begin(); iii0 != i0->end(); ++iii0)
                            if(*iii0 == *ii1)
                                exists=true;
                        //make sure it does already exist in two double connection
                     //   for ( typename list< pair<T2,T2> >::iterator ii_new = i_new->begin(); ii_new != i_new->end(); ++ii_new)
                       //     if(ii_new->first == *ii1)
                         //       exists=true;
                        
                        
                        if (!exists)
                            i_new->push_back( pair<T2,T2>(*ii1,*ii0));
                    }
                }
            }
        
//        cout<<"new cons"<<endl;
  //      int c=0;
    //    for ( typename vector< list< pair<T2,T2> > >::iterator i0 = new_cons.begin(); i0 != new_cons.end(); ++i0, ++c)
      //  {
            
        //    for ( typename list< pair<T2,T2> >::iterator ii0 = i0->begin(); ii0 != i0->end(); ++ii0)
        //        cout<<"conn2 "<<c<<" : "<<ii0->first<<" waypoint: "<<ii0->second<<endl;
         //   cout<<endl;
       // }
        //cout<<"end new cons"<<endl;
        
        return new_cons;
    }
    template vector< list< pair<unsigned int,unsigned int> > > computeNextLevelConnectionsWithMultiWayPoints( vector< std::list<unsigned int> > & connections0, vector< std::list<unsigned int> > & connections1 );
    

    template <class T,class T2>
    void cluster( fslSurface<T,T2>& surf , const unsigned int & index, const T & threshold )
    {
        if (surf.adj_verts_bi)
            surf.adj_verts.clear();
        if (surf.adj_verts.empty()) 
			surf.computeAdjacentVertices(false);
        
        //copy over adjacent vertices
        //going to delete as search through scalars
        vector<T> clusters(surf.N_vertices,0);
        vector<bool> used(surf.N_vertices,0);
        //keep track of which cluster were on
        T cluster_count=1;
        //typename vector<T>::iterator i_cl = clusters.begin();
        vector< list<T2> > connections = surf.adj_verts;
        
        unsigned int Nconns = connections.size();
        if (surf.N_vertices != Nconns)
            throw fslSurfaceException( "cluster : Number of connection does not match number of vertices" );
        
        typename vector<T>::const_iterator i_sc = surf.const_scbegin(index);
        //        for ( typename vector<T>::const_iterator i_sc = surf.const_scbegin(index); i_sc != surf.const_scbegin(index); ++i_sc, ++i_cl, ++i_con, ++v_index)
        for ( T2 v_index = 0 ; v_index <  Nconns; ++v_index)
        {
         //   cout<<"v_index "<<v_index<<endl;
            if (!used[v_index])
            {
                if (*(i_sc+v_index) >= threshold) //will enter this once per cluster
                {
                    //   clusters[v_index]=cluster_count;
                    
                    list<T2> stack;
                    stack.push_back(v_index);
                    while (!stack.empty()){
                        
                        
                        T2 current_vert = stack.back();
                        //cout<<"stack "<<current_vert<<" "<<stack.size()<<" "<<cluster_count<<endl;

                        stack.pop_back();
                        
                        for ( typename list<T2>::const_iterator i_con = connections[current_vert].begin(); i_con !=  connections[current_vert].end(); ++i_con)
                        {
                            //if scalar is above threshold
                            //cout<<"val "<<*i_con<<" "<<*(i_sc + *i_con)<<" "<<threshold<<endl;
                            used[*i_con]=1;

                            if ( (*(i_sc + *i_con)) >= threshold) 
                            {
                                stack.push_back(*i_con);

                                remove_item<T2>(connections[*i_con], current_vert);
                            }    
                            else//other remove connection with the cluster
                            {
                                fslsurface_name::remove_connection<T2>(connections, *i_con);
                               // used[*i_con]=1;
                            }
                        }
                        clusters[current_vert] = cluster_count;
                        //vert has beeen set and connections pushed to stack
                        remove_connection<T2>(connections, current_vert);
                       // used[current_vert];
                        
                    }
                    
                    ++cluster_count;
                }else { //delete all connections
                    fslsurface_name::remove_connection<T2>(connections,v_index);
                    used[v_index];
                }
            }
            
        }
        
     //   surf.scalar_data.insert(surf.scalar_data.begin(),clusters);
        surf.scalar_data.push_back(clusters);

        stringstream ss;
        ss<<threshold;
        string thresh;
        ss>>thresh;

        surf.scalar_names.push_back("clusters_"+thresh);

        
    }

    template void cluster<float, unsigned int>( fslSurface<float, unsigned int>& surf , const unsigned int & index, const float & threshold );

    
    template <class T, class T2>
    void sc_smooth_mean_neighbour( fslSurface<T,T2>& surf , const unsigned int & index )
    {
        cout<<"average using mean neighbours"<<endl;
            surf.adj_verts.clear();
        if (surf.adj_verts.empty()) 
			surf.computeAdjacentVertices(true);
        vector<T> sc_mean = surf.scalar_data[index] ; //(surf.getNumberOfVertices(),0);
        
        vector< list<T2> > connections = surf.adj_verts;
        typename vector< list<T2> >::iterator i_adj = connections.begin();
        typename std::vector< T >::const_iterator i_sc_orig = surf.const_scbegin(index);
//        vector<unsigned int >::iterator i_scN = sc_mean_N.begin();
        for (  vector<float>::iterator i_sc = sc_mean.begin(); i_sc != sc_mean.end(); ++i_sc , ++i_adj)
        {
            for (typename list<T2>::iterator ii_adj = i_adj->begin(); ii_adj != i_adj->end(); ++ii_adj)
            {
                *i_sc += *(i_sc_orig+*ii_adj);
            //    cout<<*(i_sc_orig+*ii_adj)<<endl;
            }
            *i_sc /= i_adj->size()+1;
                
        }
        
        surf.scalar_data.insert( surf.scalar_data.begin(),sc_mean);
        surf.scalar_names.insert( surf.scalar_names.begin(),surf.scalar_names[index] + "_mean_neighbour");

        
    }
    
    template void sc_smooth_mean_neighbour<float, unsigned int>( fslSurface<float, unsigned int>& surf , const unsigned int & index );
    
    
    //return the connection plus the distance
    template <class T, class T2>
    //std::list< pair<T2,T> > 
    void conn2dist(fslSurface<T,T2>& surf,const T2 & vert_index, std::list< std::pair<T2,T2> > & conns, vector< std::list< pair<T2,T> > > & index_distances  )
    {
        //first is connection
        //second is waypoint
        //assumes that we traverse the surface from 0 -> N-1
        //typename vector< vertex<T> >::const_iterator v_iter = surf.const_vbegin();
        
        
       // cout<<"COMPUTE DIST "<<vert_index<<" "<<conns.size()<<endl;
        // list< pair<T2,T> > distances;
        //look for duplicate conns
        while (! conns.empty()) {
            //cout<<"loop"<<endl;
            pair<T2,T> point = conns.front();
            
         //   cout<<point.first<<endl;
            conns.pop_front();
            
            //Only compute distances it hasn't already been computed 
            //takes adantage of the bidirectionalness of the minimal distances
            //        bool not_already_computed=true;
            //        for ( typename std::list< pair<T2,T> >::iterator i = index_distances//[vert_index].begin();  i != index_distances[vert_index].end(); ++i)
            //            if ( i->first == point.first )
            //            {//check to see if exists
            //               not_already_computed = false;
            //               break;
            //           }
            
            // if (not_already_computed)
           // cout<<"pointnew "<<point.first<<" "<<point.second<<" "<<vert_index<<endl;
            //because we traverse the surface from 0->N-i1
            if (point.first > vert_index)
            {
                bool notFoundDup=true;
                float distance = surf.L2norm(vert_index,point.second ) +  surf.L2norm(point.second,point.first );
                for (typename list< std::pair<T2, T2> >::iterator i = conns.begin(); i != conns.end();++i)
                {
                  //  cout<<"found duplicate???? "<<i->first<<" "<<i->second<<endl;

                    if ( point.first == i->first )
                    {
                    //    cout<<"found duplicate "<<i->first<<" "<<i->second<<endl;
                        notFoundDup=false;
                        vec3<float> rot_axis = surf.subtractVerts(i->second, point.second);
//                        normalize(rot_axis);
                        rot_axis.normalize();
                        //subtraction of both way points is
                        //    cout<<"rot axis "<<rot_axis.x<<" "<<rot_axis.y<<" "<<rot_axis.z<<endl;
                        ///find inplane normal by projeting P01 onto rotationa xis
                        //find base axis component by project P12 onto rotational axis
                        ///find in plane perpendicualr component by projecting P12 ontp perpendicular form traingle 1
                        
                        //                    float3 B = surf.getVert(point.second);
                        vec3<float> P0 = surf.getVertexCoord(vert_index);
                        vec3<float> P1 = surf.getVertexCoord(point.second);
                        vec3<float> P2 = surf.getVertexCoord(point.first);
                        vec3<float> P3 = surf.getVertexCoord(i->second);
                        
                        float p10_dot_rot = ((P0.x-P1.x)*rot_axis.x + (P0.y-P1.y)*rot_axis.y + (P0.z-P1.z)*rot_axis.z); 
                        //projection 
                        //      cout<<"P01do Rot "<<p10_dot_rot<<endl;
                        vec3<float> base_perp = vec3<float>(P1.x + rot_axis.x * p10_dot_rot - P0.x , \
                                                  P1.y + rot_axis.y * p10_dot_rot - P0.y , \
                                                  P1.z + rot_axis.z * p10_dot_rot - P0.z  );
                        //        cout<<"base_perp "<<P0.x<<" "<<P0.y<<" "<<P0.z<<endl;
                        //      cout<<"base_perp "<<P1.x<<" "<<P1.y<<" "<<P1.z<<endl;
                        
                        //    cout<<"base_perp "<<P1.x + rot_axis.x * p10_dot_rot<<" "<<P1.y + rot_axis.y * p10_dot_rot<<" "<<P1.z + rot_axis.z * p10_dot_rot<<endl;
                        
                        //     cout<<"base_perp "<<base_perp.x<<" "<<base_perp.y<<" "<<base_perp.z<<endl;
                       // normalize(base_perp);
                        base_perp.normalize();
                        vec3<float> p12 = subtract(P2, P1);
                        float p12_dot_rot = dot(p12,rot_axis);
                        //                    float p12_dot_rot = ((P2.x-P1.x)*rot_axis.x + (P2.y-P1.y)*rot_axis.y + (P2.z-P1.z)*rot_axis.z); 
                        float p12_dot_perp = dot(p12,base_perp);
                        
                        vec3<float> P2_rotated = vec3<float>(P1.x + rot_axis.x * p12_dot_rot + base_perp.x * p12_dot_perp, \
                                                   P1.y + rot_axis.y * p12_dot_rot + base_perp.y * p12_dot_perp, \
                                                   P1.z + rot_axis.z * p12_dot_rot + base_perp.z * p12_dot_perp);
                        
                        
                        ///------------------TEST FOR INTERSECTION WITH ROTATION EDGES---------------//
                        //do a copplanar check all points should be in same plane, don't use rotated point in case 
                        //numerical error. Use a way point and the original, one from each line
                        float t1 =0.0f;
                        
                        //make sure its not in the Z plane (i.e. x or y are all the same)
                        //simplify, dont need to test P2, because will by definition be inplabe with 0,1,2 
                        //usinf roateted could have introduced numerical errors
                        if ((P0.x==P1.x)&&(P1.x==P3.x))
                            t1 = ((P1.z - P0.z)*(P2_rotated.y-P0.y) - (P2_rotated.z-P0.z)*(P1.y-P0.y))/( (P2_rotated.z-P0.z)*(P3.y-P1.y) - (P3.z-P1.z)*(P2_rotated.y-P0.y)  );
                        else if ((P0.y==P1.y)&&(P1.y==P3.y))
                            t1 = ((P1.x - P0.x)*(P2_rotated.z-P0.z) - (P2_rotated.x-P0.x)*(P1.z-P0.z))/( (P2_rotated.x-P0.x)*(P3.z-P1.z) - (P3.x-P1.x)*(P2_rotated.z-P0.z)  );
                        else 
                            t1 = ((P1.x - P0.x)*(P2_rotated.y-P0.y) - (P2_rotated.x-P0.x)*(P1.y-P0.y))/( (P2_rotated.x-P0.x)*(P3.y-P1.y) - (P3.x-P1.x)*(P2_rotated.y-P0.y)  );
                           
                       

                        //*************DEBUGGING*******************//
                    //    cout<<"set roatated point "<<vert_index<<" "<<point.second<<" "<<i->second<<" "<<point.first<<endl;
                    //    surf.setVert(point.first,P2_rotated);
                       // cout<<"t "<<t1<<endl;
                      //  surf.setVert(point.second,float3(P1.x+t1*(P3.x-P1.x),P1.y+t1*(P3.y-P1.y),P1.z+t1*(P3.z-P1.z)) );
                        //*************DEBUGGING*******************//

                        //       cout<<"float intersection1 "<<t1<<" "<<(P1.x)<<" "<<(P1.y)<<" "<<(P1.z)<<" "<<endl;
                        //       cout<<"float intersection12 "<<t1<<" "<<(P3.x)<<" "<<(P3.y)<<" "<<(P3.z)<<" "<<endl;
                        //       cout<<"float intersection13 "<<t1<<" "<<(P3.x-P1.x)<<" "<<(P3.y-P1.y)<<" "<<(P3.z-P1.z)<<" "<<endl;
                        
                        //     cout<<"float intersection2 "<<t1<<" "<<P1.x+t1*(P3.x-P1.x)<<" "<<P1.y+t1*(P3.y-P1.y)<<" "<<P1.z+t1*(P3.z-P1.z)<<" "<<endl;
                        
                        
                        //----------------------------------------------------------------------------//
                        
                        if ((t1>=0)&& (t1<=1))
                        {
                            
                            //   cout<<"Protated "<<P2_rotated.x<<" "<<P2_rotated.y<<" "<<P2_rotated.z<<endl;
                            float dist = sqrtf( (P2_rotated.x-P0.x)*(P2_rotated.x-P0.x) + (P2_rotated.y-P0.y)*(P2_rotated.y-P0.y) + (P2_rotated.z-P0.z)*(P2_rotated.z-P0.z) );
                            //first minus second
                            // cout<<"DISTDISTIDISTI "<<vert_index<<" "<<point.first<<" "<<point.second<<" "<<dist<<endl;
                            //letter refrence to seb point on axis to point rotating
                            //                  float p_min_b_dot_n = dot(subtract(P,B),rot_axis);
                            //                float3 C(B.x + p_min_b_dot_n*rot_axis.x,  B.y + p_min_b_dot_n*rot_axis.y, B.z + p_min_b_dot_n*rot_axis.z );
                            //              float3 u(P.x - C.x, P.y - C.y, P.z - C.z);
                            //            float3 v = cross(n,u);
                            
                            if (dist < distance) distance = dist;
                        }
                        i=conns.erase(i);
                        i--;
                    }
                }
                    //            if (notFoundDup)
                    //          {   //if not more than one way point then cannot get there via two traingle use sum of edges.
                    
                    index_distances[point.first].push_back( pair<T2,T>(vert_index,distance));
                    //                distances.
                    index_distances[vert_index].push_back( pair<T2,T>(point.first,  distance ) );
               // }
            }
        }
       // return distances;
    }
    
    template <class T, class T2>
    void sc_smooth_gaussian_geodesic( fslSurface<T,T2>& surf ,const unsigned int & index, const T & sigma , const T & extent , bool run4D )
    {
        cout<<"do smooth"<<endl;
        //scaling factor for gaussian
      //  float sigma2 = sigma*sigma;
        float invsigma2 = 1.0f/(sigma*sigma);

       // float a = 1.0f / ( sigma * sqrtf(2.0*PI));
        
        vector< list<float> > distance_1;
        surf.adj_verts.clear();
        if (surf.adj_verts.empty()) 
			surf.computeAdjacentVertices(true);//only use double sided, calculate distance in each direction
        //for testing start with onnection is dist 1
       //start with just single neighbours
       // vector< list<T2> > connections = surf.adj_verts;
     //   vector< list<T2> > connections_2 = computeNextLevelConnectionsWithMultiWayPoints<T2>(surf.adj_verts,surf.adj_verts);
        vector< list< pair<T2,T2> > > connections_2 = computeNextLevelConnectionsWithMultiWayPoints<T2>(surf.adj_verts,surf.adj_verts);

        
        vector< list< pair<T2,T> > > index_distance(surf.N_vertices);
        
     //   bool not_already_computed=true;
      //  for ( typename std::list< pair<T2,T> >::iterator i = index_distances[vert_index].begin();  i != index_distances[vert_index].end(); ++i)
//            if ( i->first == point.first )
  //          {//check to see if exists
    //            not_already_computed = false;
      //          break;
       //     }
        
        //if (not_already_computed)
        //{

        
        {//do 1-neighbour ring
            
            unsigned int curr_index=0;
            typename vector< list< pair<T2,T> > >::iterator i_id = index_distance.begin();
            for ( typename vector< list<T2> >::iterator i_adj = surf.adj_verts.begin(); i_adj != surf.adj_verts.end();++i_adj,++curr_index,++i_id)
            {
                
          //      cout<<"current index "<<curr_index<<endl;
                for (typename list<T2>::iterator ii_adj = i_adj->begin(); ii_adj != i_adj->end(); ++ii_adj)
                {
                    
                    if (curr_index<*ii_adj)
                    {
                    float dist = surf.L2norm(curr_index,*ii_adj);
                 //    cout<<"hmm "<<i_adj->size()<<" "<<curr_index<<" "<<*ii_adj<<endl;
                  //   cout<<*ii_adj<<" "<<surf.L2norm(curr_index,*ii_adj) <<endl;
              
                        //Didn't add in both directions, wasn't sure if it would be worht
                        //the search throught the existance
                        i_id->push_back( pair<T2,T>(*ii_adj, dist ) );
                        index_distance[*ii_adj].push_back( pair<T2,T>(curr_index,dist) );
                    }
                }
                // cout<<endl;
            }
        }
        
   //     cout<<"print distances ring1"<<endl;
     //   for ( typename vector< list< pair<T2,T> > >::iterator i_adj = index_distance.begin(); i_adj != index_distance.end();++i_adj )
  //      {
    //        for (typename list< pair<T2,T> >::iterator ii_adj = i_adj->begin(); ii_adj != i_adj->end(); ++ii_adj)
      //      {
        //        cout<<ii_adj->first<<","<<ii_adj->second<<" ";
          //  }
   //         cout<<endl;
     //   }
        {//2-neighbour ring
            unsigned int curr_index=0;
            typename vector< list< pair<T2,T> > >::iterator i_id = index_distance.begin();
            for ( typename vector< list< pair<T2,T2> > >::iterator i_adj = connections_2.begin(); i_adj != connections_2.end();++i_adj,++curr_index,++i_id)
            {
              //  cout<<"current index21 "<<curr_index<<" "<<i_adj->size()<<endl;
               // i_id->push_back( pair<T2, T>(1,1));
                //list< pair<unsigned int,float> > con2s = 
                conn2dist<float,unsigned int>(surf, curr_index, *i_adj, index_distance);
                //i_id->merge(con2s);
                //break;
              //  cout<<"current index2 "<<curr_index<<endl;
         //       for (typename list< pair<T2,T2> >::iterator ii_adj = i_adj->begin(); ii_adj != i_adj->end(); ++ii_adj)
           //     {
         //            cout<<"hmm2 "<<i_adj->size()<<endl;
         //            cout<<ii_adj->first<<" "<<ii_adj->second<<" "<<surf.L2norm(curr_index,ii_adj->first) <<endl;
           //         i_id->push_back( pair<T2,T>(ii_adj->first, surf.conn2dist(curr_index,ii_adj->first) ) );
                    
             //   }
                // cout<<endl;
            }
        }
        cout<<"done distances"<<endl;
 //       cout<<"print distances"<<endl;
 //       for ( typename vector< list< pair<T2,T> > >::iterator i_adj = index_distance.begin(); i_adj != index_distance.end();++i_adj )
//        {
  //          for (typename list< pair<T2,T> >::iterator ii_adj = i_adj->begin(); ii_adj != i_adj->end(); ++ii_adj)
 //           {
  //              cout<<ii_adj->first<<","<<ii_adj->second<<" ";
  ///          }
  //          cout<<endl;
  //          }
        
        ///DO diksjtra algorithm 
        
        
        //now lets go through the graph and find all within a given distance
        ////initializes graph
        vector< T > graph;//(index_distance.size(),1e16);//the distance 
        vector< bool > graph_visited;//(index_distance.size(),false);//the distance 
        vector< bool > in_queue;//(index_distance.size(),false);//the distance 
        
        
        //do firts vertex then will follow with the rest
//        unsigned int curr_index=0;
        //for a given node visit all negihbours
//this will be in large loop around all vertices
      //  T2 vertex_index = 0;
        T2 Nverts = index_distance.size();
    //    Nverts=1;
        vector< vector<float> > smoothed_scalars(1,vector<float>(Nverts,0));

        vector<string> smoothed_names(surf.getNumberOfScalarData());
      //  vector< vector<float> > smoothed_scalars4D;
        if (run4D)
        {
            smoothed_scalars.assign(surf.getNumberOfScalarData(),vector<float>(Nverts,0));
        }else{

        }
        //first smoothed vectors
        typename vector<T>::iterator i_sc_smoothed = smoothed_scalars[0].begin();
        for ( T2 vertex_index =0 ; vertex_index < Nverts; ++vertex_index, ++i_sc_smoothed)
        {
            if ( vertex_index%100 == 0 )
                cout<<"Smoothing at vertex "<<vertex_index<<" ..."<<endl;
            graph.assign(Nverts,1e16);
            graph_visited.assign(Nverts, false);
            in_queue.assign(Nverts, false);
            in_queue[vertex_index]=true;
            graph[vertex_index]=0;
         //   cout<<"graph graph "<<endl;
            
         //   int count2=0;
           // for (typename vector<T>::iterator i = graph.begin(); i != graph.end(); ++i,++count2)
             //   cout<<count2<<" "<<*i<<endl;
                 
          //  typename  vector<T>::iterator i_graph = graph.begin()+vertex_index;
            
            list< pair<T,T2> > points_to_visit;//with distance
            //typename list< pair<T2,T> >::iterator i_neigh_end = index_distance[vertex_index].end();
            //this is actually just copying info from center (directtly from index_distance)
            //this is for immediate neighbours
            for (typename list< pair<T2,T> >::iterator i_neigh = index_distance[vertex_index].begin(); i_neigh != index_distance[vertex_index].end(); ++i_neigh)
            {
                graph[i_neigh->first] = i_neigh->second;//dist; 
                
                points_to_visit.push_back( pair<T,T2>(i_neigh->second,i_neigh->first) );
            }
            graph_visited[vertex_index] = true;
            
    //        cout<<"graph graph graph "<<endl;
      //      count2=0;
        //    for (typename vector<T>::iterator i = graph.begin(); i != graph.end(); ++i,++count2)
          //      cout<<count2<<" "<<*i<<endl;
            
            //as long as there are point to visit keep on going
            while (!points_to_visit.empty())
            {
                //order the queue
                points_to_visit.sort();
                //take the closest point
                pair<T,T2> p_curr = points_to_visit.front();
          //      cout<<"new point of queue "<<p_curr.second<<" "<<p_curr.first<<endl;
                T2 curr_index= p_curr.second;
                T curr_dist = p_curr.first;
                points_to_visit.pop_front();
                in_queue[curr_index] = false;
                
                typename list< pair<T2,T> >::iterator i_neigh_end = index_distance[curr_index].end();
                //this is actually just copying info from center (directtly from index_distance)
                for (typename list< pair<T2,T> >::iterator i_neigh = index_distance[curr_index].begin(); i_neigh != i_neigh_end; ++i_neigh)
                {
                    //set current distance of node to new distance plus distance at current node
                    //don't need to visit if its all ready been totally added
                    if (!graph_visited[i_neigh->first])
                    {
                      //  cout<<"enter if "<<i_neigh->first<<endl;
                        T dist = curr_dist + i_neigh->second;
                        
                        if (graph[i_neigh->first] > dist)
                        {
                       //     cout<<"set distance "<<curr_dist<<" "<<i_neigh->second<<" "<<i_neigh->first<<" "<<graph[i_neigh->first]<<" "<<dist<<endl;
                            graph[i_neigh->first] = dist;
                        }
                        
                        if (dist<= extent)
                        {
                            //if the neighbour has not been marked as visited 
                            if ((!in_queue[i_neigh->first]))
                            {
                                points_to_visit.push_back( pair<T,T2> (dist,i_neigh->first) );
                                in_queue[i_neigh->first]=true;
                            }else
                            {
                                //update distance of point in queue 
                                for (typename list< pair<T,T2> >::iterator i = points_to_visit.begin(); i!= points_to_visit.end(); ++i)
                                {
                                    if ( i->second == curr_index)
                                    {
                                        i->first = dist;
                                        break;
                                    }
                                    
                                }
                            }
                        }
                    }
                }
                graph_visited[curr_index] = true;
                
                
            }
            
        //cout<<"graph results "<<endl;
            
       //     int count=0;
         //   for (typename vector<T>::iterator i = graph.begin(); i != graph.end(); ++i,++count)
           //     cout<<count<<" "<<*i<<endl;
           
           // surf.addScalars(graph,"distances_from_vert_"+stemp);
            if (run4D)
            {
                unsigned int count=0;

                for ( typename vector< vector<T> >::iterator ii_sc = surf.scalar_data.begin(); ii_sc != surf.scalar_data.end(); ++ii_sc,++count)
                {
                    if (vertex_index==0)
                    {
                        
                        stringstream ss;
                        string ssigma,sextent;
                        ss<<sigma;
                        ss>>ssigma;
                        ss<<extent;
                        ss>>sextent;
               //         cout<<"setname "<<(surf.getScalarName(count)+"_smoothed_sigma"+ssigma+"_extent"+sextent)<<endl;
                        smoothed_names[count]=(surf.getScalarName(count)+"_smoothed_sigma"+ssigma+"_extent"+sextent);
                    }


                    //in 4D case need to assign propper smoothed
                    i_sc_smoothed = smoothed_scalars[count].begin() + vertex_index;
                    
                    float norm=0;
                    typename vector<T>::iterator i_sc = ii_sc->begin();
                    for ( typename vector<T>::iterator i_g = graph.begin(); i_g != graph.end(); ++i_g,++i_sc )
                    {
                        //make sure were within the extent
                        if ( *i_g <= extent )
                        {
                            float factor = expf( -0.5*(*i_g)*(*i_g)*invsigma2 );
                            *i_sc_smoothed += (*i_sc)*factor; 
                            norm+=factor;
                        }
                    }
                    *i_sc_smoothed/=norm; 
              
                }
                
            }else{
                
                float norm=0;
                typename vector<T>::iterator i_sc = surf.scalar_data[index].begin();
                for ( typename vector<T>::iterator i_g = graph.begin(); i_g != graph.end(); ++i_g,++i_sc )
                {
                    //make sure were within the extent
                    if ( *i_g <= extent )
                    {
                        float factor = expf( -0.5*(*i_g)*(*i_g)*invsigma2 );
                        *i_sc_smoothed += (*i_sc)*factor; 
                        norm+=factor;
                    }
                }
                *i_sc_smoothed/=norm; 
                
            }
         
           // cout<<"norm "<<norm<<endl;
        }
        
        stringstream ss;
        string ssigma,sextent;
        ss<<sigma;
        ss>>ssigma;
        ss<<extent;
        ss>>sextent;
        if (run4D)
        {
       //     cout<<"replace scalars "<<smoothed_scalars.size()<<" "<<smoothed_names.size()<<endl;
            surf.replaceScalars(smoothed_scalars,smoothed_names);

        
        }else{
            surf.insertScalars(smoothed_scalars[0],0,surf.scalar_names[index] + "_smoothed_sigma"+ssigma+"_extent"+sextent);

        }
            
      // 
      //  typename vector< list<T2> >::iterator i_adj = connections.begin();
//while ( !i_adj->empty() ) {
            
     //   }
        //        for (typename list<T2>::iterator ii_adj = i_adj->begin(); ii_adj != i_adj->end(); ++ii_adj)
//       {
//           *i_sc += *(i_sc_orig+*ii_adj);
//           //    cout<<*(i_sc_orig+*ii_adj)<<endl;
//       }

    }

    template void sc_smooth_gaussian_geodesic( fslSurface<float, unsigned int>& surf ,const unsigned int & index, const float & variance , const float & extent , bool run4D);
    
	template<class T, class T2>
      void multiplyScalarsByScalars(fslSurface<T,T2>& surf , const unsigned int & index, const fslSurface<T,T2>& surf2 , const unsigned int & index1)
    {
        
        vector<T> new_sc( surf.N_vertices , 0 );
        typename vector<T>::iterator i_sc_new= new_sc.begin();

        typename vector<T>::const_iterator i_sc1= surf2.const_scbegin(index);
        for ( typename vector<T>::const_iterator i_sc0= surf.const_scbegin(index); i_sc0 != surf.const_scend(index); ++i_sc0, ++i_sc1, ++i_sc_new)
        {
            *i_sc_new = (*i_sc0) * (*i_sc1);
        }
        surf.insertScalars(new_sc, 0, "sc_mul_sc");
        
    }
    
    template void multiplyScalarsByScalars( fslSurface<float, unsigned int>& surf ,const unsigned int & index, const fslSurface<float, unsigned int>& surf2 ,const unsigned int & index1 );

    
    template <class T, class T2>
    void subtractScalarsFromScalars(fslSurface<T,T2>& surf , const unsigned int & index, const fslSurface<T,T2>& surf2 , const unsigned int & index1, string name ){
        vector<T> new_sc( surf.N_vertices , 0 );
        typename vector<T>::iterator i_sc_new= new_sc.begin();
        
        typename vector<T>::const_iterator i_sc1= surf2.const_scbegin(index);
        for ( typename vector<T>::const_iterator i_sc0= surf.const_scbegin(index); i_sc0 != surf.const_scend(index); ++i_sc0, ++i_sc1, ++i_sc_new)
        {
            cout<<"sc "<<*i_sc0<<" "<<*i_sc1<<endl;
            *i_sc_new = (*i_sc0) - (*i_sc1);
        }
        surf.insertScalars(new_sc, 0, name);
    }
    template void subtractScalarsFromScalars( fslSurface<float, unsigned int>& surf ,const unsigned int & index, const fslSurface<float, unsigned int> & surf2 ,const unsigned int & index1,string name );

    
    
	template<class T, class T2,class T3>
	void runMarchingCubesOnAllLabels( fslSurface<T,T2>& surf, const T3* image,  const fslsurface_name::image_dims & dims, const T3 & thresh_in)
	{
		set<T3> all_labels;
		
		int nvoxels = static_cast<int>(dims.xsize * dims.ysize * dims.zsize);
		
		for (int i = 0 ; i < nvoxels; ++i)
			all_labels.insert(image[i]);
		
		for ( typename set<T3>::iterator i_lb = all_labels.begin(); i_lb != all_labels.end(); ++i_lb)
			if (*i_lb)
				marchingCubes<T,T2,T3>(surf, image, dims, static_cast<T> (*i_lb),*i_lb,EQUAL_TO);
	}
	
	template void runMarchingCubesOnAllLabels<float,unsigned int, float>( fslSurface<float,unsigned int>& surf, const float* image,  const fslsurface_name::image_dims & dims, const float & thresh_in);
	template void runMarchingCubesOnAllLabels<float,unsigned int, int>( fslSurface<float,unsigned int>& surf, const int* image,  const fslsurface_name::image_dims & dims, const int & thresh_in);
	
	
	template<class T, class T2,class T3>
	void marchingCubes( fslSurface<T,T2>& surf, const T3* image,  const fslsurface_name::image_dims & dims, const T3 & thresh_in, const T & label, const MarchingCubesMode & mode)
	{
		//	cout<<"Marching cubes input "<<thresh_in<<" "<<label<<" "<<mode<<endl;
		//	cout<<"EQUAL_TO "<<EQUAL_TO<<endl;
		//	cout<<"GREATER_THAN "<<GREATER_THAN<<endl;
		//	cout<<"GREATER_THAN_EQUAL_TO "<<GREATER_THAN_EQUAL_TO<<endl;
		
		//Tables are taken from http://local.wasp.uwa.edu.au/~pbourke/geometry/polygonise/marchingsource.cpp"
		//Via this website along with algorithm which has been adopted into fslsurface http://local.wasp.uwa.edu.au/~pbourke/geometry/polygonise/
		//--------test data-----// 
		
		
		int edgeTable[256]={
			0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
			0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
			0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
			0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
			0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
			0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
			0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
			0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
			0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
			0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
			0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
			0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
			0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
			0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
			0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
			0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
			0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
			0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
			0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
			0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
			0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
			0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
			0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
			0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
			0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
			0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
			0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
			0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
			0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
			0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
			0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
			0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };
		int triTable[256][16] =
		{{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
			{3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
			{3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
			{3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
			{9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
			{1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
			{9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
			{2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
			{8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
			{9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
			{4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
			{3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
			{1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
			{4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
			{4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
			{9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
			{1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
			{5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
			{2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
			{9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
			{0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
			{2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
			{10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
			{4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
			{5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
			{5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
			{9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
			{0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
			{1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
			{10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
			{8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
			{2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
			{7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
			{9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
			{2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
			{11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
			{9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
			{5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
			{11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
			{11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
			{1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
			{9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
			{5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
			{2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
			{0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
			{5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
			{6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
			{0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
			{3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
			{6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
			{5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
			{1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
			{10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
			{6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
			{1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
			{8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
			{7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
			{3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
			{5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
			{0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
			{9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
			{8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
			{5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
			{0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
			{6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
			{10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
			{10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
			{8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
			{1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
			{3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
			{0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
			{10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
			{0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
			{3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
			{6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
			{9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
			{8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
			{3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
			{6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
			{0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
			{10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
			{10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
			{1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
			{2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
			{7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
			{7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
			{2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
			{1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
			{11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
			{8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
			{0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
			{7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
			{10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
			{2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
			{6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
			{7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
			{2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
			{1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
			{10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
			{10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
			{0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
			{7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
			{6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
			{8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
			{9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
			{6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
			{1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
			{4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
			{10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
			{8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
			{0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
			{1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
			{8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
			{10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
			{4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
			{10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
			{5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
			{11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
			{9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
			{6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
			{7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
			{3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
			{7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
			{9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
			{3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
			{6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
			{9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
			{1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
			{4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
			{7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
			{6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
			{3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
			{0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
			{6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
			{1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
			{0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
			{11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
			{6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
			{5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
			{9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
			{1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
			{1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
			{10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
			{0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
			{5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
			{10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
			{11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
			{0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
			{9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
			{7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
			{2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
			{8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
			{9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
			{9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
			{1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
			{9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
			{9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
			{5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
			{0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
			{10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
			{2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
			{0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
			{0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
			{9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
			{5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
			{3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
			{5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
			{8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
			{0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
			{9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
			{0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
			{1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
			{3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
			{4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
			{9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
			{11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
			{11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
			{2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
			{9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
			{3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
			{1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
			{4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
			{4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
			{0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
			{3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
			{3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
			{0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
			{9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
			{1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
			{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};
		
		//find all the vertices
		//	int xmin = 10e10;
		//	int xmax = -10e10;
		//	int ymin = 10e10;
		//	int ymax = -10e10;
		//	int zmin = 10e10;
		//	int zmax = -10e10;
		//for ( int x = 0 ; x < xsize; ++x)
		//-------------TEST CODE--------------//
		/*	volume<float> im(3,3,3);
		 int tc=0;
		 for (int z=0;z<3;++z)
		 for (int y=0;y<3;++y)
		 for (int x=0;x<3;++x,++tc)
		 im.value(x,y,z)=tc;
		 dims.xsize=3;
		 dims.ysize=3;
		 dims.zsize=3;
		 
		 */
		//-------------TEST CODE--------------//
		//	vertices.clear();
		//		faces.clear();
		//		scalar_data.clear();
		//		vector_data.clear();
		
		if (surf.scalar_data.empty())
		{
			surf.scalar_names.push_back("mcubes_labels");
			
			vector<T> sc_temp;
			surf.scalar_data.push_back(sc_temp);
		}
		
		
		int ystride = static_cast<int>(dims.xsize);
		int zstride = static_cast<int>(dims.xsize * dims.ysize);
		
		int nvoxels = static_cast<int>(dims.xsize * dims.ysize * dims.zsize);
		
		T3 thresh = thresh_in;
		//cout<<"apply threshold "<<nvoxels<<" "<<mode<<endl;
		
		T3* imth = new T3[nvoxels];
		if ( mode == EQUAL_TO ){
			//	cout<<"Run marchig cubes for discrete labels"<<endl;
			for (int i = 0 ; i < nvoxels; ++i){
			//	cout<<"image "<<i<<endl;
	//			cout<<" "<<image[i]<<" "<<thresh<<" "<<(( static_cast<int>(image[i]) ==  static_cast<int>(thresh) ) ? 1 : 0)<<endl;
//if ( (( static_cast<int>(image[i]) == static_cast<int>(thresh) ) ? 1 : 0) == thresh )
//    cout<<"image match "<<i<<" "<<thresh<<endl;
				imth[i] = ( static_cast<int>(image[i]+0.5) == static_cast<int>(thresh+0.5) ) ? 1 : 0 ;
			}
			thresh=(T3)0.5;
			//	imth[i] = ( image[i] == thresh ) ? image[i] : 0 ;
		}else if (mode == GREATER_THAN) {
			for (int i = 0 ; i < nvoxels; ++i)
				imth[i] = ( image[i] > thresh ) ? image[i] : 0 ;
		}else if (mode == GREATER_THAN_EQUAL_TO) {
			for (int i = 0 ; i < nvoxels; ++i)
				imth[i] = ( image[i] >= thresh ) ? image[i] : 0 ;
		}
		//cout<<"done apply threshold "<<endl;

		
		//-------------TEST CODE--------------//
		/*
		 const float* im_ptr = im.fbegin(); 
		 for (int i = 0 ; i < nvoxels; ++i)
		 imth[i] = ( im_ptr[i] > thresh ) ? im_ptr[i] : 0 ;
		 */
		//-------------TEST CODE--------------//
		
		//	//cout<<"Done thresholding the image "<<dims.xsize<<" "<<dims.ysize<<" "<<dims.zsize<<endl;
		/*
		 Determine the index into the edge table which
		 tells us which vertices are inside of the surface
		 */
		
		//	int ElCount =0;
		int cubeindex = 0;
		float* vertVals = new float[8];
		int index  = 0 ; 
		vector< map<int, unsigned int> > index_2_index_all;
		
		for ( int z = 0 ; z < dims.zsize-1; ++z)
		{
		//cout<<"xyz "<<z<<endl;

			for ( int y = 0 ; y < dims.ysize-1; ++y)
			{
				////cout<<"xy "<<y<<endl;

				for ( int x = 0 ; x < dims.xsize-1; ++x, ++index)
				{
				//	if ((z==158)&&(y==153))
				//		//cout<<"vert "<<x<<" "<<y<<" "<<z<<endl;
					//vertVal
					//	//cout<<"new vert "<<x<<" "<<y<<" "<<z<<" "<<index<<" "<<" ."<<static_cast<int>(x + y*ystride	+ z*zstride)<<". "<<imth[index]<<" "<<image[index]<<" "<<image[static_cast<int>(x + y*ystride	+ z*zstride)]<<endl;
					////cout<<"new vert test "<<x<<" "<<y<<" "<<z<<" "<<imth[index]<<" "<<im_ptr[index]<<" "<<im_ptr[static_cast<int>(x + y*ystride	+ z*zstride)]<<endl;
					//	image[static_cast<int>(x + y*ystride	+ z*zstride)]
					cubeindex = 0;
					/*	vertVals[0] = imth[x + (y+1)*ystride	+ z*zstride   ];
					 vertVals[1] = imth[(x+1) + (y+1)*ystride	+ z*zstride   ];
					 vertVals[2] = imth[(x+1) + y*ystride	+ z*zstride   ];
					 vertVals[3] = imth[x + (y+1)*ystride	+ z*zstride   ];
					 vertVals[4] = imth[x + (y+1)*ystride	+ z*zstride   ];
					 vertVals[5] = imth[x + (y+1)*ystride	+ z*zstride   ];
					 vertVals[6] = imth[x + (y+1)*ystride	+ z*zstride   ];
					 vertVals[7] = imth[x + (y+1)*ystride	+ z*zstride   ];
					 
					 */
					
					vertVals[0] = imth[index + ystride   ];
					vertVals[1] = imth[index + ystride  +1];
					vertVals[2] = imth[index + 1];
					vertVals[3] = imth[index ];
					vertVals[4] = imth[index +zstride + ystride    ];
					vertVals[5] = imth[index +zstride + ystride  +1];
					vertVals[6] = imth[index +zstride + 1       ];
					vertVals[7] = imth[index +zstride           ];
				//	if ((z==158)&&(y==153))
				//		//cout<<"cube vals:"<<endl;
					//		for (int i =0 ; i < 8 ; ++i)
					//			//cout<<vertVals[i]<<" ";
					//		//cout<<endl;
					if ( vertVals[0] < thresh) cubeindex |= 1; 
					if ( vertVals[1] < thresh) cubeindex |= 2; 
					if ( vertVals[2] < thresh) cubeindex |= 4; 
					if ( vertVals[3] < thresh) cubeindex |= 8; 
					
					if ( vertVals[4] < thresh) cubeindex |= 16; 
					if ( vertVals[5] < thresh) cubeindex |= 32; 
					if ( vertVals[6] < thresh) cubeindex |= 64; 
					if ( vertVals[7] < thresh) cubeindex |= 128; 
					
					
					map<int, unsigned int> index_2_index;
				//	if ((z==158)&&(y==153))
				//		//cout<<"check edge table "<<index_2_index_all.max_size()<<endl;
					if (edgeTable[cubeindex] == 0)
					{
				//		if ((z==158)&&(y==153))
				//		//cout<<"enter edge table check "<<index_2_index_all.size()<<" "<<index_2_index.size()<<" "<<cubeindex<<endl;
						index_2_index_all.push_back(index_2_index);
				//		if ((z==158)&&(y==153))
				//		//cout<<"go to continue "<<endl;
						continue;
					}
				//	//cout<<"check edge table end"<<endl;
					////cout<<"new vert "<<x<<" "<<y<<" "<<z<<" "<<index<<" "<<" ."<<static_cast<int>(x + y*ystride	+ z*zstride)<<". "<<imth[index]<<" "<<image[index]<<" "<<image[static_cast<int>(x + y*ystride	+ z*zstride)]<<endl;

					T2 vert_count=0;
					T2 vert_index = surf.vertices.size();
					/// Find the vertices where the surface intersects the cube 
					if (edgeTable[cubeindex] & 1)
					{	//edge p0-p1
						if (z==0)
						{
							vertex<T> vnew(label);
							vertexInterp(thresh, \
										 vec3<float>(x*dims.xdim,(y+1)*dims.ydim, z*dims.zdim),\
										 vertVals[0], \
										 vec3<float>((x+1)*dims.xdim,(y+1)*dims.ydim, z*dims.zdim), \
										 vertVals[1], vnew ); 
							surf.vertices.push_back(vnew);
							surf.scalar_data[0].push_back(label);
							index_2_index.insert( pair<int, unsigned int> (0,vert_index+vert_count) );
							++vert_count;
						}else {
							index_2_index.insert( pair<int,unsigned int>(0, index_2_index_all[x +y*ystride + (z-1)*zstride].find(4)->second) );
							
						}
						
					}
					if (edgeTable[cubeindex] & 2)
					{	//edge p1-p2
						if (z==0)
						{
							vertex<T> vnew(label);
							
							vertexInterp(thresh, \
										 vec3<float>((x+1)*dims.xdim,(y+1)*dims.ydim, z*dims.zdim), \
										 vertVals[1], \
										 vec3<float>((x+1)*dims.xdim,y*dims.ydim, z*dims.zdim), \
										 vertVals[2], vnew ); 
							surf.vertices.push_back(vnew);
							surf.scalar_data[0].push_back(label);
							
							index_2_index.insert( pair<int, unsigned int> (1,vert_index+vert_count) );
							++vert_count;
							
						}else {
							index_2_index.insert( pair<int,unsigned int>(1, index_2_index_all[x +y*ystride + (z-1)*zstride].find(5)->second) );
							
						}
						
						
						
					}if (edgeTable[cubeindex] & 4)
					{	//edge p2-p3
						if ((z==0)&&(y==0))
						{
							vertex<T> vnew(label);
							
							vertexInterp(thresh, \
										 vec3<float>((x+1)*dims.xdim,y*dims.ydim, z*dims.zdim), \
										 vertVals[2], \
										 vec3<float>(x*dims.xdim,y*dims.ydim, z*dims.zdim), \
										 vertVals[3], vnew ); 
							surf.vertices.push_back(vnew);
							surf.scalar_data[0].push_back(label);
							
							index_2_index.insert( pair<int, unsigned int> (2,vert_index+vert_count) );
							++vert_count;
							
						}else if (z==0){
							index_2_index.insert( pair<int,unsigned int>(2, index_2_index_all[x + (y-1)*ystride ].find(0)->second) );
						}else{	
							index_2_index.insert( pair<int,unsigned int>(2, index_2_index_all[x +y*ystride + (z-1)*zstride].find(6)->second) );
						}
						
					}if (edgeTable[cubeindex] & 8)
					{	//edge p3-p0
						if ((x==0)&&(z==0))
						{
							vertex<T> vnew(label);
							
							vertexInterp(thresh, \
										 vec3<float>(x*dims.xdim,y*dims.ydim, z*dims.zdim), \
										 vertVals[3], \
										 vec3<float>(x*dims.xdim,(y+1)*dims.ydim, z*dims.zdim), \
										 vertVals[0], vnew ); 
							surf.vertices.push_back(vnew);
							surf.scalar_data[0].push_back(label);
							index_2_index.insert( pair<int, unsigned int> (3,vert_index+vert_count) );
							++vert_count;
							
						}else if (z==0) {
							index_2_index.insert( pair<int,unsigned int>(3, index_2_index_all[x-1 +y*ystride].find(1)->second) );
						}else{
							index_2_index.insert( pair<int,unsigned int>(3, index_2_index_all[x +y*ystride + (z-1)*zstride].find(7)->second) );
						}
						
						
					}	if (edgeTable[cubeindex] & 16)
					{	//edge p4-p5
						vertex<T> vnew(label);
						
						vertexInterp(thresh, \
									 vec3<float>(x*dims.xdim,(y+1)*dims.ydim, (z+1)*dims.zdim), \
									 vertVals[4], \
									 vec3<float>((x+1)*dims.xdim,(y+1)*dims.ydim, (z+1)*dims.zdim), \
									 vertVals[5], vnew ); 
						surf.vertices.push_back(vnew);
						surf.scalar_data[0].push_back(label);
						index_2_index.insert( pair<int, unsigned int> (4,vert_index+vert_count) );
						++vert_count;
						
						
						
					}
					if (edgeTable[cubeindex] & 32)
					{	//edge p1-p2
						//	//cout<<"edge 5 "<<endl;
						
						vertex<T> vnew(label);
						
						vertexInterp(thresh, \
									 vec3<float>((x+1)*dims.xdim,(y+1)*dims.ydim, (z+1)*dims.zdim), \
									 vertVals[5], \
									 vec3<float>((x+1)*dims.xdim,y*dims.ydim, (z+1)*dims.zdim), \
									 vertVals[6], vnew ); 
						surf.vertices.push_back(vnew);
						surf.scalar_data[0].push_back(label);
						index_2_index.insert( pair<int, unsigned int> (5,vert_index+vert_count) );
						++vert_count;
						
						
						
					}if (edgeTable[cubeindex] & 64)
					{	//edge p2-p3
						if (y==0)
						{
							vertex<T> vnew(label);
							
							vertexInterp(thresh, \
										 vec3<float>((x+1)*dims.xdim,y*dims.ydim, (z+1)*dims.zdim), \
										 vertVals[6], \
										 vec3<float>(x*dims.xdim,y*dims.ydim, (z+1)*dims.zdim), \
										 vertVals[7], vnew ); 
							surf.vertices.push_back(vnew);
							surf.scalar_data[0].push_back(label);
							index_2_index.insert( pair<int, unsigned int> (6,vert_index+vert_count) );
							++vert_count;
							
						}else {
							index_2_index.insert( pair<int,unsigned int>(6, index_2_index_all[x +(y-1)*ystride + z*zstride].find(4)->second) );
							
						}
						
						
					}if (edgeTable[cubeindex] & 128)
					{	//edge p3-p0
						////cout<<"edge 7 "<<endl;
						if (x==0)
						{
							vertex<T> vnew(label);
							
							vertexInterp(thresh, \
										 vec3<float>(x*dims.xdim,y*dims.ydim, (z+1)*dims.zdim), \
										 vertVals[7], \
										 vec3<float>(x*dims.xdim,(y+1)*dims.ydim, (z+1)*dims.zdim), \
										 vertVals[4], vnew );
							surf.vertices.push_back(vnew);
							surf.scalar_data[0].push_back(label);
							index_2_index.insert( pair<int, unsigned int> (7,vert_index+vert_count) );
							++vert_count;
							
						}	else {
							////cout<<"find previous index "<<index_2_index_all.size()<<" "<<y*ystride<<" "<<z*zstride<<" "<<(x-1) +y*ystride + z*zstride<<" "<<index<<endl;
							index_2_index.insert( pair<int,unsigned int>(7, index_2_index_all[(x-1) +y*ystride + z*zstride].find(5)->second) );
							
						}
						////cout<<"done edge 7 "<<endl;
						
					}
					if (edgeTable[cubeindex] & 256){	//edge p0-p4
						
						if (x==0)
						{
							vertex<T> vnew(label);
							
							vertexInterp(thresh, \
										 vec3<float>(x*dims.xdim,(y+1)*dims.ydim, z*dims.zdim), \
										 vertVals[0], \
										 vec3<float>(x*dims.xdim,(y+1)*dims.ydim, (z+1)*dims.zdim), \
										 vertVals[4], vnew ); 
							surf.vertices.push_back(vnew);
							surf.scalar_data[0].push_back(label);
							index_2_index.insert( pair<int, unsigned int> (8,vert_index+vert_count) );
							++vert_count;
							
						}else {
							index_2_index.insert( pair<int,unsigned int>(8, index_2_index_all[(x-1) +y*ystride + z*zstride].find(9)->second) );
							
						}
						
						
					}
					if (edgeTable[cubeindex] & 512){	//edge p1-p5
						
						vertex<T> vnew(label);
						
						vertexInterp(thresh, \
									 vec3<float>((x+1)*dims.xdim,(y+1)*dims.ydim, z*dims.zdim), \
									 vertVals[1], \
									 vec3<float>((x+1)*dims.xdim,(y+1)*dims.ydim, (z+1)*dims.zdim), \
									 vertVals[5], vnew ); 
						surf.vertices.push_back(vnew);
						surf.scalar_data[0].push_back(label);
						index_2_index.insert( pair<int, unsigned int> (9,vert_index+vert_count) );
						++vert_count;
						
						
						
					}
					if (edgeTable[cubeindex] & 1024){	//edge p2-p6
						if (y==0) {
							
							vertex<T> vnew(label);
							
							vertexInterp(thresh, \
										 vec3<float>((x+1)*dims.xdim,y*dims.ydim, z*dims.zdim), \
										 vertVals[2], \
										 vec3<float>((x+1)*dims.xdim,y*dims.ydim, (z+1)*dims.zdim), \
										 vertVals[6], vnew );
							surf.vertices.push_back(vnew);
							surf.scalar_data[0].push_back(label);
							index_2_index.insert( pair<int, unsigned int> (10,vert_index+vert_count) );
							++vert_count;
							
						}else {
							index_2_index.insert( pair<int,unsigned int>(10, index_2_index_all[x +(y-1)*ystride + z*zstride].find(9)->second) );
							
						}
						
					}
					if (edgeTable[cubeindex] & 2048){	//edge p3-p7
						if ((y==0)&&(x==0)) {
							
							vertex<T> vnew(label);
							
							vertexInterp(thresh, \
										 vec3<float>(x*dims.xdim,y*dims.ydim, z*dims.zdim), \
										 vertVals[3], \
										 vec3<float>(x*dims.xdim,y*dims.ydim, (z+1)*dims.zdim), \
										 vertVals[7], vnew ); 
							surf.vertices.push_back(vnew);
							surf.scalar_data[0].push_back(label);
							index_2_index.insert( pair<int, unsigned int> (11,vert_index+vert_count) );
							++vert_count;
							
						}else if (x==0) {
							index_2_index.insert( pair<int,unsigned int>(11, index_2_index_all[x +(y-1)*ystride + z*zstride].find(8)->second) );
							
						}else {
							index_2_index.insert( pair<int,unsigned int>(11, index_2_index_all[(x-1) +y*ystride + z*zstride].find(10)->second) );
							
						}
						
						
						
					}
					index_2_index_all.push_back(index_2_index);
					//	//cout<<x<<" "<<y<<" "<<z<<" index2index "<<endl;
					//		for (map<int,unsigned int>::iterator i=index_2_index.begin(); i != index_2_index.end(); ++i)
					//			//cout<<i->first<<" "<<i->second<<endl;
			//	//cout<<"do trinagles "<<x<<" "<<y<<" "<<z<<" "<<cubeindex<<" "<<triTable[cubeindex][0]<<endl;
					//do trinagles now that edges/verts are done
					for (int i_tri = 0 ; triTable[cubeindex][i_tri] != -1 ; i_tri+=3)
					{
						
						//						//cout<<"Push back triangle i_tri "<<i_tri<<endl;
						//	//cout<<index_2_index.find(triTable[cubeindex][i_tri  ])->second<<" "<<index_2_index.find(triTable[cubeindex][i_tri +1 ])->second<<" "<<index_2_index.find(triTable[cubeindex][i_tri +2 ])->second<<endl;
						//						//cout<<triTable[cubeindex][i_tri]<<endl;
						//						//cout<<triTable[cubeindex][i_tri+1]<<endl;
						//						//cout<<triTable[cubeindex][i_tri+2]<<endl;
						//						//cout<<triTable[cubeindex][i_tri]<<endl;
						//						//cout<<triTable[cubeindex][i_tri+1]<<endl;
						//						//cout<<index_2_index.find(triTable[cubeindex][i_tri+0])->second<<endl;
						//						//cout<<index_2_index.find(triTable[cubeindex][i_tri+1])->second<<endl;
						
						//						//cout<<index_2_index.find(triTable[cubeindex][i_tri+2])->second<<endl;
						surf.faces.push_back( index_2_index.find(triTable[cubeindex][i_tri  ])->second);
						surf.faces.push_back( index_2_index.find(triTable[cubeindex][i_tri+1])->second);
						surf.faces.push_back( index_2_index.find(triTable[cubeindex][i_tri+2])->second);
						
					}
				//	//cout<<"done triangles"<<endl;
					/*	ElCount++;
					 if (ElCount==12)
					 {
					 x=1e6;
					 y=1e6;
					 z=1e6;
					 }
					 */
					
				}
				map<int, unsigned int> index_2_index;
				index_2_index_all.push_back(index_2_index);
				
				++index;//because only going to -1
			}
			map<int, unsigned int> index_2_index;
			for (int i=0;i<ystride;++i)
				index_2_index_all.push_back(index_2_index);
			index+=ystride;
		}
		//cout<<"Number of Vertices "<<surf.vertices.size()<<endl;	
		//	cout<<"faces "<<endl;
		//	for (vector<unsigned int>::iterator i = faces.begin(); i!=faces.end(); ++i) {
		//		cout<<*i<<" ";
		//	}
		//cout<<"Number of Triangles "<<faces.size()<<endl;
		//	scalar_data.push_back(sc_label);
		surf.N_vertices = surf.vertices.size();
		surf.N_triangles = surf.faces.size()/3;
		delete[] vertVals;
		
		delete[] imth;
		/* Create the triangle
		 ntriang = 0;
		 for (i=0;triTable[cubeindex][i]!=-1;i+=3) {
		 triangles[ntriang].p[0] = vertlist[triTable[cubeindex][i  ]];
		 triangles[ntriang].p[1] = vertlist[triTable[cubeindex][i+1]];
		 triangles[ntriang].p[2] = vertlist[triTable[cubeindex][i+2]];
		 ntriang++;
		 }
		 
		 return(ntriang); */
		
	}
	template 	void marchingCubes<float,unsigned int, float>( fslSurface<float,unsigned int>& surf, const float* image,  const fslsurface_name::image_dims & dims, const float & thresh, const float & label, const MarchingCubesMode & mode);
	
	
	
	
	
	
	template<class T>
	void vertexInterp( const float & thresh, const vec3<T> & v1 , const float  & val1, const vec3<T> & v2, const float & val2, vertex<T> & vout )
	{
		//vertex vout;
		//cout<<"interp in ("<<v1.x<<","<<v1.y<<","<<v1.z<<") -> ("<<v2.x<<","<<v2.y<<","<<v2.z<<")"<<endl;
		//cout<<"vals "<<val1<<" "<<val2<<endl;
		T t = (thresh - val1)/(val2 - val1);
		vout.x = v1.x + t * (v2.x - v1.x);
		vout.y = v1.y + t * (v2.y - v1.y);
		vout.z = v1.z + t * (v2.z - v1.z);
		//cout<<"interp out ("<<vout.x<<","<<vout.y<<","<<vout.z<<")"<<endl;
		
		//return vout;
		
	}
	
	
	template<class T, class T2,class T3>
	void mapScalars(fslSurface<T,T2>& surf, const unsigned int & sc_index, const std::map<T3,T3> & sc_map , bool  setToCurrentSc )
	{
		vector<T> new_scalars( surf.scalar_data[sc_index].size() );
		typename vector<T>::iterator i_new = new_scalars.begin();
		typename vector<T>::iterator i_end = surf.scalar_data[sc_index].end();
		
		for ( typename vector<T>::iterator i = surf.scalar_data[sc_index].begin(); i != i_end; ++i, ++i_new)
		{
			*i_new = sc_map.find(*i)->second; 
		}
		surf.scalar_data.push_back(new_scalars);
		if (setToCurrentSc)
		{
			surf.setScalars(surf.scalar_data.size()-1);
		}
	}
	template void mapScalars<float, unsigned int, float>(fslSurface<float, unsigned int>& surf, const unsigned int & sc_index, const std::map<float,float> & sc_map , bool  setToCurrentSc );
	template void mapScalars<float, int, float>(fslSurface<float, int>& surf, const unsigned int & sc_index, const std::map<float,float> & sc_map , bool  setToCurrentSc );

	template<class T, class T2>
	void maskScalars(fslSurface<T,T2>& surf, const unsigned int & sc_index, const short* mask, const fslsurface_name::image_dims  dims, bool  setToCurrentSc )
	{
		//Uses nearest neighbour interpolation
		vector<T> new_scalars( surf.scalar_data[sc_index].size() );
		typename vector<T>::iterator i_new = new_scalars.begin();
		typename vector<T>::iterator i_sc = surf.scalar_data[sc_index].begin();
		//	typename vector<T>::iterator i_end = scalar_data[sc_index].end();
		int ystride = dims.xsize;
		int zstride = dims.xsize*dims.ysize;
		//typename vector< vertex >::iterator i_v = vertices.begin();
		for ( typename vector< vertex<T> >::iterator i_v = surf.vertices.begin(); i_v != surf.vertices.end(); ++i_v, ++i_new,++i_sc)
		{
			//cout<<"xyx "<<i_v->x<<" "<<
			int x = static_cast<int>(i_v->x /dims.xdim + 0.5);
			int y = static_cast<int>(i_v->y /dims.ydim + 0.5);
			int z = static_cast<int>(i_v->z /dims.zdim + 0.5);
			
			if ( mask[ x + y*ystride + z * zstride] == 0 )
				*i_new = 0;
			else 
				*i_new = *i_sc;
		}		
		surf.scalar_data.push_back(new_scalars);
		surf.scalar_names.push_back("masked scalars");
	}
	template void maskScalars<float, unsigned int>(fslSurface<float, unsigned int>& surf, const unsigned int & sc_index, const short* mask, const fslsurface_name::image_dims  dims, bool  setToCurrentSc );
	template void maskScalars<float, int>(fslSurface<float, int>& surf, const unsigned int & sc_index, const short* mask, const fslsurface_name::image_dims  dims, bool  setToCurrentSc );

	
	
	template<class T, class T2>
	void maskSurface( fslSurface<T,T2> & surf , const short* mask, const fslsurface_name::image_dims  dims, bool  setToCurrentSc )
	{
		int vert_ind = 0;
		for ( typename vector< vertex<T> >::iterator i_v = surf.vertices.begin(); i_v != surf.vertices.end(); ++i_v, ++vert_ind)
		{
			surf.vertices.erase(i_v);
			--i_v;
		}
	}
	template void maskSurface<float, unsigned int>( fslSurface<float, unsigned int>& surf , const short* mask, const fslsurface_name::image_dims  dims, bool  setToCurrentSc );
	template void maskSurface<float, int>( fslSurface<float, int>& surf , const short* mask, const fslsurface_name::image_dims  dims, bool  setToCurrentSc );

	
    template<class T,class T2, class T3>
    void insertScalars(fslsurface_name::fslSurface<T,T2>& surf , const NEWIMAGE::volume4D<T> & time_series, const string & time_name )
    {
        int nTimesPoints = time_series.tsize();
        int xsize = time_series.xsize();
        int ysize = time_series.ysize();
        int nvals = xsize*ysize;
        if (static_cast<unsigned int>(xsize*ysize) != surf.getNumberOfVertices())
            throw fslSurfaceException( ("Tried to append 4D scalars to surfaces, however unequal number of points. x : " + num2str(xsize) +", y : " + num2str(ysize) + ", n_verts : " + num2str(surf.getNumberOfVertices()) +".").c_str()  );
        
        for ( int t =0 ; t<nTimesPoints; ++t)
        {
            vector<T3> sc(nvals);
            const T* i_ptr = time_series[t].fbegin();
            typename vector<T3>::iterator i_sc = sc.begin();
            for ( int y=0; y<ysize; ++y)
                for ( int x =0; x< xsize ; ++x,++i_ptr,++i_sc)
                {
                    *i_sc = *i_ptr;
                }//cout<<"txy "<<t<<" "<<x<<" "<<y<<" "<<*i_ptr<<" "<<time_series[t].value(x,y,0)<<endl;

            surf.insertScalars(sc,t,time_name+"_tp_"+num2str(t) );
        }
        
    }
    template void insertScalars<float,unsigned int,float>(fslsurface_name::fslSurface<float,unsigned int>& surf , const NEWIMAGE::volume4D<float> & time_series , const string & time_name);
    
    
    template<class T,class T2>
    void sampleTimeSeries(fslsurface_name::fslSurface<T,T2>& surf , const NEWIMAGE::volume4D<T> & time_series, const std::string & time_name )
    {
        //NEWIMAGE::volume4D<T> imtest;
        //imtest = time_series;
       // vector< vertex<T> >::const_iterator i;
        float3 dims=float3(time_series.xdim(), time_series.ydim(), time_series.zdim());
        unsigned int Nt = static_cast<unsigned int>(time_series.tsize());
        unsigned int Nverts = surf.getNumberOfVertices();
       // vector< vector<T> > all_scalars(Nt,vector<T>(Nverts,0));
      //  typename vector< vector<T> >::iterator i_all =  all_scalars.begin();
        for ( unsigned int t = 0 ; t < Nt ; ++t)//, ++i_all)
        {cout<<"time "<<t<<endl;
            vector<T> sc(Nverts);
            typename vector<T>::iterator i_sc = sc.begin();
            for ( typename vector< vertex<T> >::const_iterator i_v = surf.const_vbegin(); i_v!=surf.const_vend(); ++i_v)
            {
                *i_sc = time_series[t].interpolate(i_v->x/dims.x, i_v->y/dims.y, i_v->z/dims.z);
       //         imtest[t].value(i_v->x/dims.x, i_v->y/dims.y, i_v->z/dims.z)=t;

            }
            surf.insertScalars(sc,t,time_name+"_tp_"+num2str(t) );
        }
     //   save_volume4D(imtest,time_name+"_test");
    }
    template void sampleTimeSeries<float,unsigned int>(fslsurface_name::fslSurface<float,unsigned int>& surf , const NEWIMAGE::volume4D<float> & time_series , const string & time_name);
    
    template<class T,class T2>
    void writeTimeSeries( const fslsurface_name::fslSurface<T,T2>& surf , NEWIMAGE::volume4D<T> & time_series, const string & outname )
    {
        
        float3 dims=float3(time_series.xdim(), time_series.ydim(), time_series.zdim());

        unsigned int Nt = static_cast<unsigned int>(time_series.tsize());
        
        NEWIMAGE::volume<float> w_sum(time_series.xsize(),time_series.ysize(),time_series.zsize());
        w_sum=0;
       
       // NEWIMAGE::volume4D<T> time_series_new(time_series.xsize(),time_series.ysize(),time_series.zsize(),time_series.tsize());
        //time_series_new=0;
        unsigned int vert_ind=0;
        for ( typename vector< vertex<T> >::const_iterator i_v = surf.const_vbegin(); i_v!=surf.const_vend(); ++i_v,++vert_ind)
        {
            unsigned int x = static_cast<unsigned int>(i_v->x /dims.x);
            unsigned int y = static_cast<unsigned int>(i_v->y /dims.y);
            unsigned int z = static_cast<unsigned int>(i_v->z /dims.z);
            
          //  vec3<float> centre( (x+0.5)*dims.x, (y+0.5)*dims.y, (z+0.5)*dims.z );
            
            
            vec3<float> diff = i_v->coord() - vec3<float>((x+0.5)*dims.x, (y+0.5)*dims.y, (z+0.5)*dims.z ) ;
            float dist = diff.norm();
            w_sum.value(i_v->x/dims.x, i_v->y/dims.y, i_v->z/dims.z)+=dist;
            
            for (unsigned int t = 0 ; t < Nt ; ++t)
            {
                time_series[t].value(i_v->x/dims.x, i_v->y/dims.y, i_v->z/dims.z) += surf.getScalar(t,vert_ind);
                    
            }
//  *i_sc = time_series[t].interpolate(i_v->x/dims.x, i_v->y/dims.y, i_v->z/dims.z);
            //         imtest[t].value(i_v->x/dims.x, i_v->y/dims.y, i_v->z/dims.z)=t;
            
        }
        for (unsigned int t = 0 ; t < Nt ; ++t)
        time_series[t]/=w_sum;

       save_volume4D(time_series,outname);
    }

       template void writeTimeSeries<float,unsigned int>( const fslsurface_name::fslSurface<float,unsigned int>& surf , NEWIMAGE::volume4D<float> & time_series, const string & outname );
}





