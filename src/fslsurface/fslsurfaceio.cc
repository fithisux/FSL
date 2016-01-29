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

#include <fslsurface.h>
#include <fslsurface_structs.h>
#include <fslsurfaceio.h>
#include <fslvtkio/fslvtkio.h>

extern "C" {
#include <giftiio/gifti_io.h>
}

using namespace fslvtkio;
using namespace std;
using namespace fslsurface_name;

namespace fslsurface_name {
	
  // SUPPORT FUNCTION
  std::string safe_string(char *cptr) {
    if (cptr==NULL) return "";
    return string(cptr);
  }

	template<class T, class T2>
	int read_surface( fslSurface<T,T2>& surf, const string & filename)
	{
		surf.scalar_indices.push_back(vector<int>());
		//cog.x = cog.y = cog.z = 0.0;
		//	endian_format = fslSurface::machineEndianness();
		//cout<<"read surface "<<endl;
		//ifstream fin(filename.c_str());
		string last_3 = filename.substr(filename.size()-3, 3);
		////cout<<"extension "<<last_3<<endl;
		if ( last_3 == "ply" )
			readPLY(surf,filename);
		else if ( last_3 == "vtk" )
			readVTK(surf,filename);
		else if ( last_3 == "gii" )
			readGIFTI(surf,filename);
		else if ( last_3 == ".gz" )
		{
			last_3 = filename.substr(filename.size()-6, 3);
			if ( last_3 == "gii" )
				readGIFTI(surf,filename);
			else {
				throw fslSurfaceException("Unrecognized file format.");
			}
			
		}else {
			throw fslSurfaceException("Unrecognized file format.");
		}
		//cout<<"done read surface "<<endl;
		
		return 0;
	}
	template int read_surface<float,unsigned int>(fslSurface<float,unsigned int>& surf, const string & filename);	
	template int read_surface<float,int>(fslSurface<float,int>& surf, const string & filename);	
	
	
	template<class T, class T2>
	int readPLY( fslSurface<T,T2>& surf, const string & filename){
		/*
		 ifstream fin(filename.c_str());
		 readPLYHeader( fin );
		 unsigned int list_index=0;
		 for (unsigned int i = 0 ; i < element_names.size(); i++)
		 {
		 if (element_names[i]=="vertex")
		 {
		 N_vertices = readPLYVertices(fin,i) ;
		 }else if (element_names[i]=="face")
		 {
		 //assuming all triangles
		 faces = readPLYListElement<T2>(fin,i,list_index) ;
		 N_triangles = faces.size()/3;
		 }else {
		 if (element_types[i][0]=="list"){
		 
		 
		 }else {
		 
		 }
		 }
		 }
		 */
		return 0;
	}

	template<class T, class T2>
	int readVTK( fslSurface<T,T2>& surf,  const string & filename){
		fslvtkIO* f_vtk = new fslvtkIO(filename,fslvtkIO::POLYDATA);
		surf.cog.x = surf.cog.y = surf.cog.z = 0.0;
		////cout<<"read vtk "<<endl;
		vector<T> verts = f_vtk->getPointsAsVector<T>();
		surf.N_vertices = verts.size()/3;
		for ( typename vector<T>::iterator i = verts.begin(); i != verts.end(); i+=3)
		{
			surf.cog.x += *i;
			surf.cog.y += *(i+1);
			surf.cog.z += *(i+2);
			surf.vertices.push_back(vertex<T>(*i,*(i+1),*(i+2)));
		}
		surf.cog.x /= surf.N_vertices;
		surf.cog.y /= surf.N_vertices;
		surf.cog.z /= surf.N_vertices;
		
		vector< vector<T2> > polys = f_vtk->getPolygonsAsVectorOfVectors<T2>();
		for ( typename vector< vector<T2> >::iterator i = polys.begin(); i != polys.end(); ++i )
			surf.faces.insert(surf.faces.end(), i->begin(), i->end());
		surf.N_triangles = surf.faces.size()/3;
        
		// cout<<"add back scalars "<<f_vtk->getScalars<T>().size()<<endl;
		if ( f_vtk->getScalars<T>().size() > 0 )
		  {
		    vector<T> vsc = f_vtk->getScalars<T>();
		    //   for (typename vector<T>::iterator i_sc = vsc.begin(); i_sc != vsc.end(); ++i_sc)
		    //       cout<<*i_sc<<" ";
		    //   cout<<endl;
		    surf.scalar_data.push_back(f_vtk->getScalars<T>());
		    surf.scalar_names.push_back("vtk_scalars");
		    surf.scalar_indices.back().push_back( surf.scalar_data.size()-1 );
		    
		  }
		//	else {
		//		surf->setScalars(0);
		//	}
		if (surf.scalar_data.size()>0)
		  surf.setScalars(0);
		
		delete f_vtk;
		return 0;
	}


  // GIFTI support functions

  int fsl_gifti_add_to_meta( giiMetaData * md, const std::string& name, const std::string& value,
			     int replace ) {
    if ((name!="") && (value!="")) {
      return gifti_add_to_meta(md,name.c_str(),value.c_str(),replace);
    }
    return 1;
  }

	
	template<class T, class T2>
	int readGIFTI( fslSurface<T,T2>& surf, const string & filename){
        cout<<"read gifti "<<endl;
		gifti_image* gii_surf = gifti_read_image(filename.c_str(), 1);
	
		//if (	gii_surf->labeltable == NULL)
		////cout<<"label table "<<gii_surf->labeltable.length<<endl;
	
        int Nentries = gii_surf->labeltable.length;
		int* key = gii_surf->labeltable.key;
		float* rgba = gii_surf->labeltable.rgba;
		//************Read Top-Level Meta Data ***************//
        
      //  giiMetaData surf_meta = gii_surf->meta;
      //  int length= surf_meta.length;//number of meta data fields
      //  for ( int i = 0 ; i < length ; ++i )
       // {
            
            
       // }
        //***********done meta data ************************//
		
		
		for ( int i = 0 ; i < Nentries ; ++i,++key,++rgba)
		{
			float r = *(rgba++);
			float g = *(rgba++);
			float b = *(rgba++);
			////cout<<"rgba "<<*key<<" "<<r<<" "<<g<<" "<<b<<" "<<*rgba<<endl;
			
			surf.dataTable[*key] = float4(r,g,b,*rgba);
		}
		////cout<<"done reading"<<endl;
		if (gii_surf == NULL)
			throw("Error in Reading GIFTI image");
		int NumberOfDataArrays = gii_surf->numDA;
		////cout<<"Numver of Data arrays "<<NumberOfDataArrays<<endl;
		//	////cout<<"GIFTI meta data : length "<<gii_surf->meta.length<<" "<<string(*(gii_surf->meta.name))<<" "<<string(*(gii_surf->meta.value))<<endl;
		int sc_type=0;
		for ( int i_da = 0 ; i_da < NumberOfDataArrays ; ++i_da )
		{
			////cout<<"Read data arrray : "<<i_da<<endl;
			int Nvals = gii_surf->darray[i_da]->nvals;
			
			////cout<<"GIFTI data array ("<<i_da<<") : data array: intent : "<< gifti_intent_to_string(gii_surf->darray[i_da]->intent)<<endl;
			int intent = gii_surf->darray[i_da]->intent;
			
			int datatype = gii_surf->darray[i_da]->datatype;
			int num_dim = gii_surf->darray[i_da]->num_dim;
			
			int* dims = gii_surf->darray[i_da]->dims; 
			
			////cout<<"GIFTI data array ("<<i_da<<") : datatype : "<<surf->datatype_str[datatype]<<endl;

//---------------THESE ARE COMMENTED OUT BECAUSE NO USED --------------------//            
//            int ind_ord = gii_surf->darray[i_da]->ind_ord;
//			int encoding = gii_surf->darray[i_da]->encoding;
//			int endian = gii_surf->darray[i_da]->endian;
//-------------------------------------------------------------------------//
			//	////cout<<"GIFTI data array ("<<i_da<<") : encoding : "<<encoding_str[encoding]<<endl;
			if ( intent == NIFTI_INTENT_POINTSET )
			{	
				//////cout<<"Number of Coordinate Sysstems "<<gii_surf->darray[i_da]->numCS<<endl;
				giiCoordSystem coordsys;
				for (int i_c = 0 ; i_c < gii_surf->darray[i_da]->numCS ; ++i_c)
				{
					vector<double> xfm(16,0);
					int count=0;
					for ( int i = 0 ; i<4 ; ++i)
						for ( int j = 0 ; j<4 ; ++j, ++count)
						{
							xfm[count] = gii_surf->darray[i_da]->coordsys[i_c]->xform[i][j];
							//////cout<<xfm[count]<<endl;
						}
					surf.v_coord_sys.push_back(xfm);
					char *cptr;
					string sname; 
					cptr = gii_surf->darray[i_da]->coordsys[i_c]->dataspace ;
					sname = safe_string(cptr);
					surf.csys_dspace = sname ;
					//	surf->v_csys_dspace.push_back(string( gii_surf->darray[i_da]->coordsys[i_c]->dataspace ));
					cptr = gii_surf->darray[i_da]->coordsys[i_c]->xformspace ;
					sname = safe_string(cptr);
					surf.v_csys_xfmspace.push_back( sname );
					
					
				}
				if (gii_surf->darray[i_da]->numCS >0 ) {
					coordsys =  * (*gii_surf->darray[i_da]->coordsys);
				}

				// READ META DATA
				std::string valstr;
				valstr=safe_string(gifti_get_meta_value(&(gii_surf->darray[i_da]->meta),"AnatomicalStructurePrimary"));
				if (valstr!="") surf.setAnatomicalName(valstr);
				valstr=safe_string(gifti_get_meta_value(&(gii_surf->darray[i_da]->meta),"AnatomicalStructureSecondary"));
				if (valstr!="") surf.setAnatomicalName2(valstr);
				valstr=safe_string(gifti_get_meta_value(&(gii_surf->darray[i_da]->meta),"GeometricType"));
				if (valstr!="") surf.setGeometryName(valstr);
				
				//////cout<<"read coord system  "<<coordsys.dataspace<<" "<<coordsys.xformspace<<endl;
				
				
				//vertices.clear();
				bool valid_dim=false;
				for (int i = 0 ; i<num_dim;i++)
				{
					if ( dims[i] == 3)
						valid_dim = true;
				}
				if (!valid_dim)
					throw	fslSurfaceException("Intent is PointSet, but not 3 dimensions");
				
				float *data = (float*)(gii_surf->darray[i_da]->data);

				surf.cog.x = surf.cog.y = surf.cog.z = 0.0;
				
				//write over vertices
				if ( static_cast<int>(surf.N_vertices) == Nvals/3 )
				{	
				        int noff=surf.N_vertices;
					typename vector< vertex<T> >::iterator i_v = surf.vertices.begin();
					for (unsigned int val = 0 ; val < surf.N_vertices ; ++val, ++data,++i_v)
					{
					  if (gii_surf->darray[i_da]->ind_ord == GIFTI_IND_ORD_COL_MAJOR) {
						i_v->x = *data;
						surf.cog.x+=*data;
						
						i_v->y = *(data + noff);
						surf.cog.y+=*(data + noff);

						i_v->z = *(data + 2*noff);
						surf.cog.z+=*(data + 2*noff);
					  } else {
						i_v->x = *data;
						surf.cog.x+=*data;
						
						++data;
						i_v->y = *data;
						surf.cog.y+=*data;
						
						++data;
						i_v->z = *data;
						surf.cog.z+=*data;
					  }
						
					}
					
				}else {//if different append vertices
		
                    unsigned int Nprev= surf.N_vertices; 

					surf.N_vertices += (dims[0]==3) ? dims[1] : dims[0];
        
                    surf.vertices.resize(surf.N_vertices);
		    int noff=surf.N_vertices;
                    for (typename vector< vertex<T> >::iterator i_v =  surf.vertices.begin()+Nprev ; i_v != surf.vertices.end(); ++i_v,++data)
                    {
					  if (gii_surf->darray[i_da]->ind_ord == GIFTI_IND_ORD_COL_MAJOR) {
						float x = *data;
						surf.cog.x+=x;
						
						float y = *(data + noff);
						surf.cog.y+=y;
						
						float z = *(data + 2*noff);
						surf.cog.z+=z;
						*i_v = vertex<T>(x,y,z);


					  } else {
						float x = *data;
						surf.cog.x+=*data;
						
						++data;
						float y = *data;
						surf.cog.y+=*data;
						
						++data;
						surf.cog.z+=*data;
						*i_v = vertex<T>(x,y,*data);
					  }

		    }
					if (surf.scalar_data.size()>0)
						surf.setScalars(0);
				}
				surf.cog.x /= surf.N_vertices;
				surf.cog.y /= surf.N_vertices;
				surf.cog.z /= surf.N_vertices;
				
			}else if ( intent == NIFTI_INTENT_TRIANGLE )
			{
				bool valid_dim=false;
				for (int i = 0 ; i<num_dim;i++)
				{
					//	////cout<<"dims[i] "<<dims[i]<<endl;
					if ( dims[i] == 3)
						valid_dim = true;
				}
				if (!valid_dim)
					throw	fslSurfaceException("Intent is Triangle, but not 3 dimensions");
				
				int* data = (int*)(gii_surf->darray[i_da]->data);
				surf.N_triangles = Nvals / 3;
				surf.faces.resize(Nvals);
                
				// READ META DATA
				std::string valstr;
				valstr=safe_string(gifti_get_meta_value(&(gii_surf->darray[i_da]->meta),"TopologicalType"));
				if (valstr!="") surf.setTopologyName(valstr);

				if (gii_surf->darray[i_da]->ind_ord == GIFTI_IND_ORD_COL_MAJOR) {
				  int noff=surf.N_triangles;
				  for (typename vector<T2>::iterator i_f = surf.faces.begin(); i_f != surf.faces.end(); ++i_f,++data)
				    { 
				      // awkward read for a transpose matrix
				      *i_f = static_cast<unsigned int> (*data);
				      ++i_f;
				      *i_f = static_cast<unsigned int> (*(data+noff));
				      ++i_f;
				      *i_f = static_cast<unsigned int> (*(data+2*noff));
				    }
				} else {
				  for (typename vector<T2>::iterator i_f = surf.faces.begin(); i_f != surf.faces.end(); ++i_f,++data)
				    {
				      //	////cout<<"loading faces "<<*data<<endl;
				      *i_f = static_cast<unsigned int> (*data);
				    }
				}
		//	}else if ( ( intent == NIFTI_INTENT_SHAPE ) || (intent == NIFTI_INTENT_TIME_SERIES ) || ((intent == NIFTI_INTENT_NONE)&&(datatype == NIFTI_TYPE_FLOAT32)) ){
             //   cout<<"done tris"<<endl;
			}else if ( ( intent == NIFTI_INTENT_SHAPE ) || (intent == NIFTI_INTENT_TIME_SERIES ) || (intent == NIFTI_INTENT_NORMAL) || ((intent == NIFTI_INTENT_NONE)&&(datatype == NIFTI_TYPE_FLOAT32)) ){
			  // TODO - FIGURE OUT WHAT IS GOING ON (!!) AND SEE IF COL MAJOR ORDER FIX IS REQUIRED
				//decide whether vertices have been read already, if not create, if so set to scalar
				cout<<"read scalar data that is float32"<<endl;
				if (datatype != NIFTI_TYPE_FLOAT32 ) 
					throw fslSurfaceException("Invalid Data Type ");
				
				float* data = (float*)(gii_surf->darray[i_da]->data);
				vector<T> scalars(Nvals);
             
                    
                    
				if ( surf.N_vertices == 0 ){
				  surf.vertices.resize(Nvals);
				  
				  typename vector<T>::iterator i_sc = scalars.begin();
				  for (typename vector< vertex<T> >::iterator i_v = surf.vertices.begin(); i_v != surf.vertices.end();++i_v,++i_sc,++data)
				    {
				      *i_v = vertex<T>(*data);
				      *i_sc = static_cast<T>(*data);
				    }
				  
				}else if ( static_cast<unsigned int>(Nvals) == surf.N_vertices ) {
				  if (surf.scalar_data.empty())
				    {
				      typename vector< vertex<T> >::iterator i_v = surf.vertices.begin();
				      for (typename vector<T>::iterator i_sc = scalars.begin(); i_sc != scalars.end(); ++data, ++i_sc,++i_v)
					{
					  i_v->sc = static_cast<T> (*data);
					  *i_sc = static_cast<T>(*data);
					}
				    }else {//non vertex scalara
				    for (typename vector<T>::iterator i_sc = scalars.begin(); i_sc != scalars.end(); ++i_sc,++data)
				      {
					*i_sc = static_cast<T>(*data);
				      }
				  }
				  
				}else {//non vertex scalars for scalar data that does not belong to vertices
				  //ive separated out float and int
				  //throw fslSurfaceException("Tried to add scalars to mismatched number of vertices that is not 0//");
				  vector<float> sc_float(scalars.size());//convert to float
				  vector<float>::iterator i_sc2 = sc_float.begin();
				  for (typename vector<T>::iterator i_sc = scalars.begin(); i_sc != scalars.end(); ++i_sc,++data,++i_sc2)
				    {
				      // sc_float.push_back(static_cast<float>(*data));
				      *i_sc2 = static_cast<float>(*data);
				    }
				  sc_type=1;//used at the moment to know where to put the names
				  surf.nonvert_float_sc_data.push_back( sc_float );
				  //	for (int val = 0 ; val < Nvals ; ++val, ++data )
				  //		scalars.push_back(*data);
				}
				//		////cout<<"meta data"<<endl;
				
            string name;
            giiMetaData sc_meta = gii_surf->darray[i_da]->meta;
            int index=-1;
            for (int i_meta = 0 ; i_meta < sc_meta.length; ++i_meta)
            {
                if (!strcmp(sc_meta.name[i_meta],"name")) {
                    index=i_meta;
                }
            }

                ///index define for meta data above    
                if (index == -1)
				{
					if ( intent == NIFTI_INTENT_SHAPE ) 
						name="NIFTI_INTENT_SHAPE";
					else if (intent == NIFTI_INTENT_TIME_SERIES )
						name="NIFTI_INTENT_TIME_SERIES";
					else if (intent == NIFTI_INTENT_NORMAL )
						name="NIFTI_INTENT_NORMAL";
					else if (intent == NIFTI_INTENT_NONE)
						name="NIFTI_INTENT_NONE";
					else 
						name="name";
					
				}else{
                    name=sc_meta.value[index];
              
				}
                if (sc_type==0)
                {
                    surf.scalar_data.push_back( scalars );
                    surf.scalar_indices.back().push_back( surf.scalar_data.size()-1 );
                    surf.scalar_names.push_back(name);
                
                }else if (sc_type==1)
                {
                    //scalars were added above to float 
                    surf.nonvert_float_sc_data_names.push_back( name );

                }
                     
		
                    
			}else if ( ( intent == NIFTI_INTENT_LABEL ) || ((intent == NIFTI_INTENT_NONE)&&(datatype == NIFTI_TYPE_INT32)) ){
				//////cout<<"intent label found "<<endl;
				if ( num_dim != 1 )
					throw fslSurfaceException("Invalid Number of Dimensions for Label Intent. Should be 1");
				
                int* data = (int*)(gii_surf->darray[i_da]->data);
				vector<T> scalars(Nvals);
                
				if ( ( dims[0] != static_cast<int>(surf.N_vertices) ) && ( surf.N_vertices != 0 ) )
				{
                   // throw fslSurfaceException("Invalid Number for Vertices. Should be equal to number of vertices, unless there are 0 vertices so far.");

                    vector<int> sc_int(scalars.size());//convert to float
                    vector<int>::iterator i_sc2 = sc_int.begin();
                    for (typename vector<T>::iterator i_sc = scalars.begin(); i_sc != scalars.end(); ++i_sc,++data,++i_sc2)
                    {
                        // sc_float.push_back(static_cast<float>(*data));
                        *i_sc2 = static_cast<float>(*data);
                    }
                    sc_type=1;//used at the moment to know where to put the names
                    surf.nonvert_int_sc_data.push_back( sc_int );
                    //	for (int val = 0 ; val < Nvals ; ++val, ++data )
                    
                    
                }else{
                    //	for (int i = 0 ; i<num_dim;i++)
                    //	{
                    int* data = (int*)(gii_surf->darray[i_da]->data);
                    for (typename vector<T>::iterator i_sc = scalars.begin(); i_sc != scalars.end(); ++data, ++i_sc)
                    {
                        *i_sc = static_cast<int>(*data);
                    }
                }
				
                //currently can only get in here in thes case
			//	if ( datatype != NIFTI_TYPE_INT32 )
			//		throw fslSurfaceException("Invalid Data Type. GIFTI standard require labels set to NIFTI_TYPE_INT32");
				
				
		
				
				//		//cout<<"meta data"<<endl;
				giiMetaData sc_meta = gii_surf->darray[i_da]->meta;
				int index = -1;
				for (int i_meta = 0 ; i_meta < sc_meta.length; ++i_meta)
				{
					if (!strcmp(sc_meta.name[i_meta],"name")) {
						//cout<<"found name "<<sc_meta.name[i_meta]<<endl;
						index=i_meta;
						
					}
				//	cout<<"name "<<sc_meta.name[i_meta]<<endl;
				//	cout<<"name "<<sc_meta.value[i_meta]<<endl;
					
				}
                string name="none";
				if (index == -1)
				{
					if (intent == NIFTI_INTENT_LABEL )
						name=("NIFTI_INTENT_LABEL");
					else if (intent == NIFTI_INTENT_NONE)
						name=("NIFTI_INTENT_NONE");
					else 
						name=("name");
					
				}else{
					name=(sc_meta.value[index]);
				}
                
                
                if (sc_type==0)
                {
                    surf.scalar_data.push_back( scalars );
                    surf.scalar_indices.back().push_back( surf.scalar_data.size()-1 );
                    surf.scalar_names.push_back(name);
                    
                }else if (sc_type==1)
                {
                    //scalars were added above to float 
                    surf.nonvert_int_sc_data_names.push_back( name );
                    
                }
                
                
                
				
			}else if ( intent == NIFTI_INTENT_VECTOR ){
			
				if ( num_dim != 2 )
					throw fslSurfaceException("Invalid Number of Dimensions for Vector Intent. Should be 2");
				
				if ( ( dims[0] != static_cast<int>(surf.N_vertices) ) && ( surf.N_vertices != 0 ) )
					throw fslSurfaceException("Invalid Number for Vectors. Should be equal to number of vertices, unless there are 0 vertices so far.");
				
				if ( ( dims[1] != 3 ) )
					throw fslSurfaceException("Invalid Number for Dimension of Vectors. Should be 3.");
				
				
				if ( datatype != NIFTI_TYPE_FLOAT32 )
					throw fslSurfaceException("Invalid Data Type. Briview only supports vectors of type NIFTI_TYPE_FLOAT32 at the moment");
				
				// TODO: NEED TO IMPLEMENT THE FIX FOR COL MAJOR ORDERING HERE!
				float* data = (float*)(gii_surf->darray[i_da]->data);
				
				vector<T> vecs(Nvals);
                for (typename vector<T>::iterator i = vecs.begin(); i != vecs.end(); ++i,++data)
                {
                    *i=*data;

                }
                surf.vector_data.push_back(vecs);
                
                //get names
                giiMetaData sc_meta = gii_surf->darray[i_da]->meta;
				int index = -1;
				for (int i_meta = 0 ; i_meta < sc_meta.length; ++i_meta)
				{
             //       cout<<"meta i "<<i_meta<<" "<<sc_meta.value[i_meta]<<endl;
					if (!strcmp(sc_meta.name[i_meta],"name")) {
						index=i_meta;
					}
				}
				if (index == -1)
				{
                    surf.vector_names.push_back("vectors");
				}else{
               		surf.vector_names.push_back(sc_meta.value[index]);
				}

			}
			
			//				int               datatype;   /* numerical type of Data values         */
			//				int               ind_ord;    /* lowest Dim to highest, or reverse     */
			//				int               num_dim;    /* level of DimX applied                 */
			//				int               dims[6];    /* dimension lengths (first num_dim set) */
			//				int               encoding;   /* format of Data on disk                */
			//				int               endian;     /* endian, if binary Encoding            */
			//				char            * ext_fname;  /* external filename, in cur directory   */
			//				long long         ext_offset; /* offset of data within external file   */
			
			/* elements */
			//				giiMetaData       meta;
			//				giiCoordSystem ** coordsys;   /* array of pointers to giiCoordSystem   */
			//				void            * data;       /* unencoded, uncompressed, swapped      */
			
			/* extras */
			//				long long         nvals;      /* number of values (product of Dims)    */
			///				int               nbyper;     /* number of bytes per value             */
			//				int               numCS;      /* number of giiCoordSystem structs      */
			//				nvpairs           ex_atrs;    /* extra attributes                      */
			
			//				gifti_intent_is_valid
			//				NIFTI_INTENT_POINTSET
			//				NIFTI_INTENT_TRIANGLE
		}
		//cout<<"Sc Indices Size "<<surf->scalar_indices.size()<<" "<<surf->scalar_indices.back().size()<<endl;
		gifti_free_image(gii_surf);
		
        return 1;
		//cout<<"done GIFTI read"<<endl;
	}
	
	template<class T, class T2>
	unsigned int readPLYVertices( fslSurface<T,T2> & surf, ifstream & fin , const unsigned int & index)
	{
		surf.vertices.clear();
		surf.cog.x = surf.cog.y = surf.cog.z = 0.0;
		unsigned int N = surf.element_sizes[index];
		unsigned int N_prop = surf.element_num_props[index];
		for ( unsigned int i=0 ; i < N ;i++)
		{
			vertex<T> vert;
			for ( unsigned int j=0 ; j< N_prop; j++)
			{
				T temp;
				
				fin>>temp;
				//	//cout<<"temp: "<<temp<<endl;
				switch (j) {
					case 0:
						vert.x=temp;
						break;
					case 1:
						vert.y=temp;
						break;
					case 2:
						vert.z=temp;
						break;
					case 3:
						if (N_prop==4)
							vert.sc=temp;
						else {
							vert.nx=temp;
						}
						break;
					case 4:
						vert.ny=temp;
						break;
					case 5:
						vert.nz=temp;
						break;
					case 6:
						vert.sc=temp;
						break;
					default:
						break;
				}
			}
			surf.cog.x+=static_cast<float>(vert.x);
			surf.cog.y+=static_cast<float>(vert.y);
			surf.cog.z+=static_cast<float>(vert.z);
			surf.vertices.push_back(vert);
		}
		surf.cog.x /= N;
		surf.cog.y /= N;
		surf.cog.z /= N;
		
		
		return N;
	}
	
	template<class T, class T2>
	int writePLY( const fslSurface<T,T2> & surf, const string & filename)
	{
		return 1; 
	}
	
	template<class T, class T2>
	int writePLYHeader( const fslSurface<T,T2>& surf, ofstream & fout,const string & version, const string & comment )
	{
		fout<<"ply"<<endl;
		
		//---------write numeric format---------------//
		fout<<"format ";
		switch (fslSurface<T,T2>::format) {
			case fslSurface<T,T2>::ascii:
				fout<<"ascii "<<version<<endl;
				break;
			case fslSurface<T,T2>::bigEndian:
				fout<<"binary_big_endian "<<version<<endl;
				break;
			case fslSurface<T,T2>::littleEndian:
				fout<<"binary_little_endian "<<version<<endl;
				break;				
			default:
				throw fslSurfaceException("Unsupported Data Format");
				break;
		} 
		//-------------------------------------------//
		fout<<"comment "<<comment<<endl;
		//writeElements();
		//----------------write vertex info------------//
		fout<<"element vertex "<<surf.vertices.size()/3<<endl;
		fout<<"property float x"<<endl;
		fout<<"property float y"<<endl;
		fout<<"property float z"<<endl;
		fout<<"property float nx"<<endl;
		fout<<"property float ny"<<endl;
		fout<<"property float nz"<<endl;
		fout<<"property float sc"<<endl;
		
		//assuming only triangles
		fout<<"element face "<<surf.faces.size()/3<<endl; 
		fout<<"property list uchar int vertex_index"<<endl;
		fout<<"end_header"<<endl;
	}
	
	
	template<class T, class T2>
	int writeGIFTI( const fslSurface<T,T2> & surf, const std::string & filename, int enc)
	{
	  //, const int & intent, const int & dtype
	  
	  int numDA = surf.scalar_data.size()   + surf.nonvert_int_sc_data.size() + surf.nonvert_float_sc_data.size()+ surf.vector_data.size() + 2 ; //2 is for vertices and faces
	  //        
	  int Nvertices=surf.N_vertices, intentcode=NIFTI_INTENT_POINTSET;
	  unsigned int da_index = -1;
	  if (Nvertices<=0) { // does not contain geometry info
	    numDA-=2; // no geometry
	    if (surf.vector_data.size()>0) { intentcode=NIFTI_INTENT_VECTOR; Nvertices=surf.vector_data[0].size(); }
	    if (surf.scalar_data.size()>0) { intentcode=NIFTI_INTENT_SHAPE; Nvertices=surf.scalar_data[0].size(); }
	    if (Nvertices<=0) { throw fslSurfaceException("Could not find valid geometry, scalar or vector data to write to file."); }
	  }
	  int dims[] = { Nvertices ,3,0,0,0,0 };
	  cout<<"numda2 "<<numDA<<endl;
	  gifti_image *gim = gifti_create_image( numDA, intentcode, NIFTI_TYPE_FLOAT32, 2, dims, 2 );
	
	  if (surf.N_vertices>0) {
	    da_index++;
	    //    cout<<"write gifti"<<endl;
	    //---------------------POINTSET----------------------//
	    //double xfm[16]={ 1.0,0.0,0.0,0.0, 0.0,1.0,0.0,0.0, 0.0,0.0,1.0,0.0, 0.0,0.0,0.0,1.0 };
	    gim->darray[da_index]->encoding = enc;
	    //giiCoordSystem* csys= new giiCoordSystem();
		  
	    //string temp("NIFTI_XFORM_SCANNER_ANAT");
	    //	(csys->dataspace)="NIFTI_XFORM_SCANNER_ANAT";//const_cast<char*>(temp.c_str());
	    //	(csys->xformspace)="NIFTI_XFORM_SCANNER_ANAT";//const_cast<char*>(temp.c_str());
	    unsigned int count=0;
	    vector< vector<double> >::const_iterator i_csys = surf.v_coord_sys.begin();
	    for (vector<string>::const_iterator i_x = surf.v_csys_xfmspace.begin(); i_x != surf.v_csys_xfmspace.end(); ++i_x,++i_csys, ++count)
	      {
		char* cs=const_cast<char *>(i_x->c_str());
		gifti_add_empty_CS(gim->darray[da_index]);
		gim->darray[da_index]->coordsys[count]->dataspace = cs;
		gim->darray[da_index]->coordsys[count]->xformspace = cs;

		//	cout<<"sys "<<i_csys->at(0)<<" "<<i_csys->at(1)<<endl;
		memcpy(gim->darray[da_index]->coordsys[count]->xform, &(i_csys->at(0)), 16*sizeof(double));
		      
		//cout<<"done memcpy"<<endl;
		      
		//		string type= "NIFTI_XFORM_SCANNER_ANAT";
		//		char* cs=const_cast<char *>(type.c_str());
		//gim->darray[da_index]->coordsys=&(csys);
		//		gifti_add_empty_CS(gim->darray[da_index]);
		//		gim->darray[da_index]->coordsys[da_index]->dataspace = cs;
		//	gim->darray[da_index]->coordsys[da_index]->xformspace = cs;
		//memcpy(gim->darray[da_index]->coordsys[da_index]->xform, xfm, 16*sizeof(double));
	      }
		  
	    cout<<"write vertcies "<<surf.N_vertices<<endl;
		  
	    gim->darray[da_index]->nvals = dims[0]*dims[1];
	    //gim->darray[da_index]->numCS = 0;
	    gim->darray[da_index]->nbyper = sizeof(float);

	    fsl_gifti_add_to_meta(&(gim->darray[da_index]->meta), "AnatomicalStructurePrimary",surf.anatomical_name,1);
	    fsl_gifti_add_to_meta(&(gim->darray[da_index]->meta), "AnatomicalStructureSecondary",surf.anatomical_name2,1);
	    fsl_gifti_add_to_meta(&(gim->darray[da_index]->meta), "GeometricType",surf.geometry_name,1);

	    for (unsigned int i =0; i < surf.N_vertices ; ++i)
	      {
		//memcpy(static_cast<float*>(gim->darray[da_index]->data)+i*3, &(vertices[0]) + i, sizeof(T)*3);
		*(static_cast<float*>(gim->darray[da_index]->data) + 3*i) = static_cast<float> (surf.vertices[i].x);
		*(static_cast<float*>(gim->darray[da_index]->data) + 3*i+1) = static_cast<float> (surf.vertices[i].y);
		*(static_cast<float*>(gim->darray[da_index]->data) + 3*i+2) = static_cast<float> (surf.vertices[i].z);
		//	gim->darray[da_index][1]=i_v.x;
	      }
	    cout<<"write traingles "<<surf.N_triangles<<endl;
	    //---------------------Triangles----------------------//
	    da_index++;
	    gim->darray[da_index]->intent = NIFTI_INTENT_TRIANGLE;	
	    gim->darray[da_index]->encoding = enc;	
	    gim->darray[da_index]->dims[0] = surf.N_triangles;
	    gim->darray[da_index]->dims[1] = 3;
	    gim->darray[da_index]->nvals = surf.N_triangles*3;
	    gim->darray[da_index]->datatype = NIFTI_TYPE_INT32;
	    gim->darray[da_index]->nbyper = sizeof(int);
	    gim->darray[da_index]->numCS=0;
		  
	    fsl_gifti_add_to_meta(&(gim->darray[da_index]->meta), "TopologicalType",surf.topology_name,1);

	    cout<<surf.N_triangles<<" "<<surf.N_vertices<<endl;
	    if (surf.N_triangles != surf.N_vertices)
	      {
		cout<<"resize trinagle array"<<endl;
		delete[] static_cast<int*>(gim->darray[da_index]->data);
		gim->darray[da_index]->data = new int[surf.N_triangles*3];
		//	int dimstr[] = { N_triangles ,3,0,0,0,0 };
	      }
		  
	    //		if (numeric_l) {
	    cout<<"write faces "<<surf.faces.size()<<endl;
		  
	    //		}
	    //		for (unsigned int i = 0 ; i< N_triangles*3 ; ++i)
	    unsigned int i = 0;
	    for (typename vector<T2>::const_iterator i_f = surf.faces.begin(); i_f != surf.faces.end(); ++i_f, ++i)
	      {
		//cout<<"face "<<i<<endl;
		// cout<<"face "<<i<<" "<<surf.faces.size()<<" "<<*i_f<<" "<<static_cast<int>(*i_f)<<endl;
		*(static_cast<int*>(gim->darray[da_index]->data) + i ) = static_cast<int>(*i_f);
		//  cout<<"face2"<<endl;
	      }
	    //cout<<"done face "<< surf.scalar_data.size()<<" "<<surf.scalar_data[0].size()<<endl;
	    //	memcpy(gim->darray[da_index]->data, &(faces[0]), sizeof(T2)*N_triangles);
	    //	gim->darray[0][1]=i_v.x;
	  }
	  //---------------------write scalar data---------------------//
	  da_index++;
	  unsigned int i_sc = 0;
	  for (  ; i_sc < surf.scalar_data.size(); ++i_sc)
	    {
		    
	      //unsigned int da_index=i_sc+2;
	      cout<<"scalar "<<i_sc<<" "<<da_index<<endl;
		    
	      gim->darray[da_index]->intent = NIFTI_INTENT_SHAPE;	
	      gim->darray[da_index]->encoding = enc;	
	      gim->darray[da_index]->num_dim = 1;
	      gim->darray[da_index]->dims[0] = Nvertices;
	      gim->darray[da_index]->dims[1] = 1;
	      gim->darray[da_index]->nvals = Nvertices;
	      gim->darray[da_index]->datatype = NIFTI_TYPE_FLOAT32;
	      gim->darray[da_index]->nbyper = sizeof(float);
	      //string
	      gim->darray[da_index]->numCS=0;
	      //giiMetaData* meta_d = new giiMetaData();
	      cout<<"set scalar names "<<surf.scalar_names.size()<<" "<<surf.scalar_names[i_sc]<<endl;
	      string s2=surf.scalar_names[i_sc];
	      fsl_gifti_add_to_meta(&(gim->darray[da_index]->meta), "name",surf.scalar_names[i_sc],1);
	      cout<<"do memcpy "<<endl;
	      memcpy(static_cast<float*>(gim->darray[da_index]->data), &(surf.scalar_data[i_sc][0]), sizeof(T)*Nvertices);
	      ++da_index;
	    }
	  //  unsigned int i_sc_init=i_sc+1;
	  i_sc=0;
		
	  for (  ; i_sc < ( surf.nonvert_float_sc_data.size()); ++i_sc)
	    {
		    
	      //unsigned int da_index=i_sc+i_sc_init+ 2;
	      cout<<"scalar float "<<i_sc<<" "<<da_index<<endl;
		    
	      gim->darray[da_index]->intent = NIFTI_INTENT_NONE;	
	      gim->darray[da_index]->encoding = enc;	
	      gim->darray[da_index]->num_dim = 1;
	      gim->darray[da_index]->dims[0] = surf.nonvert_float_sc_data[i_sc].size();
	      gim->darray[da_index]->dims[1] = 1;
	      gim->darray[da_index]->nvals = surf.nonvert_float_sc_data[i_sc].size();
	      gim->darray[da_index]->datatype = NIFTI_TYPE_FLOAT32;
	      gim->darray[da_index]->nbyper = sizeof(float);
	      //string
	      gim->darray[da_index]->numCS=0;
	      //giiMetaData* meta_d = new giiMetaData();
	      cout<<"set scalar names "<<surf.nonvert_float_sc_data_names.size()<<endl;
	      string s2=surf.nonvert_float_sc_data_names[i_sc];
	      fsl_gifti_add_to_meta(&(gim->darray[da_index]->meta), "name",surf.nonvert_float_sc_data_names[i_sc],1);
	      cout<<"do memcpy"<<endl;
	      memcpy(static_cast<float*>(gim->darray[da_index]->data), &(surf.nonvert_float_sc_data[i_sc][0]), sizeof(float)*surf.nonvert_float_sc_data_names[i_sc].size());
		    
	      ++da_index;
		    
	    }
	  cout<<"scalar int "<<surf.nonvert_int_sc_data.size()<<endl;
	  //i_sc_init+=i_sc+1;
	  i_sc=0;
	  for (  ; i_sc < (surf.nonvert_int_sc_data.size()); ++i_sc)
	    {
	      cout<<"scalar "<<i_sc<<" "<<surf.nonvert_int_sc_data.size()<<" "<<surf.nonvert_int_sc_data_names.size()<<endl;
	      //unsigned int da_index=i_sc+i_sc_init+2;
	      cout<<"da_index "<<da_index<<endl;
		    
	      gim->darray[da_index]->intent = NIFTI_INTENT_NONE;	
	      gim->darray[da_index]->encoding = enc;	
	      gim->darray[da_index]->num_dim = 1;
	      gim->darray[da_index]->dims[0] = surf.nonvert_int_sc_data[i_sc].size();
	      gim->darray[da_index]->dims[1] = 1;
	      gim->darray[da_index]->nvals = surf.nonvert_int_sc_data[i_sc].size();
	      gim->darray[da_index]->datatype = NIFTI_TYPE_INT32;
	      gim->darray[da_index]->nbyper = sizeof(int);
	      //string
	      gim->darray[da_index]->numCS=0;
	      //giiMetaData* meta_d = new giiMetaData();
	      cout<<"set scalar int names "<<surf.nonvert_int_sc_data_names.size()<<" "<<i_sc<<endl;
	      string s2=surf.nonvert_int_sc_data_names[i_sc];
	      fsl_gifti_add_to_meta(&(gim->darray[da_index]->meta), "name",surf.nonvert_int_sc_data_names[i_sc],1);
	      cout<<"do memcpy "<<surf.nonvert_int_sc_data[i_sc][0]<<" "<<surf.nonvert_int_sc_data[i_sc][1]<<endl;
	      //int count=0;
	      //for ( vector<int>::iterator ii_sc = surf.nonvert_int_sc_data[i_sc].begin(); ii_sc != surf.nonvert_int_sc_data[i_sc].end//();++ii_sc,++count)
	      //        {
	      //              *(gim->darray[da_index]->data + count) = *ii_sc;
		    
	      //            }
	      memcpy(static_cast<int*>(gim->darray[da_index]->data), &(surf.nonvert_int_sc_data[i_sc][0]), sizeof(int)*surf.nonvert_int_sc_data[i_sc].size());
	      ++da_index;
		    
	    }
	  cout<<"done scalar"<<endl;
		
	  //---------------------write vector data---------------------//
	  unsigned int offset = surf.scalar_data.size()   + surf.nonvert_int_sc_data.size();
		
	  for ( unsigned int i_vec = 0 ; i_vec < surf.vector_data.size(); ++i_vec)
	    {
	      cout<<"vector "<<i_vec<<endl;
		    
	      unsigned int da_index=offset + i_vec +2;
	      gim->darray[da_index]->intent = NIFTI_INTENT_VECTOR;	
	      gim->darray[da_index]->encoding = enc;	
	      gim->darray[da_index]->num_dim = 2;
	      gim->darray[da_index]->dims[0] = Nvertices;
	      gim->darray[da_index]->dims[1] = 3;
	      gim->darray[da_index]->nvals = 3*Nvertices;
	      gim->darray[da_index]->datatype = NIFTI_TYPE_FLOAT32;
	      gim->darray[da_index]->nbyper = sizeof(float);
	      //string
	      gim->darray[da_index]->numCS=0;
	      //giiMetaData* meta_d = new giiMetaData();
	      cout<<"access vector names "<<surf.vector_names[i_vec]<<endl;
	      string s2=surf.vector_names[i_vec];
	      cout<<"done access vector names"<<endl;
		    
	      fsl_gifti_add_to_meta(&(gim->darray[da_index]->meta), "name",surf.vector_names[i_vec],1);
	      cout<<"do memcpy"<<endl;
	      memcpy(static_cast<float*>(gim->darray[da_index]->data), &(surf.vector_data[i_vec][0]), sizeof(T)*3*Nvertices);
		    
	    }
		
	  // cout<<"write gifti4"<<endl;
		
	  //cout<<"done scalars "<<endl;
	  gifti_write_image( gim , filename.c_str(), 1) ;
	  //ive added this to manually free the coordinatre system, seems ot give error when using fgifti_free_image
	  //cout<<"free stuf "<<endl;
	  //    cout<<"write gifti4.5"<<endl;
	  if ((gim->darray[0]->coordsys==NULL) || (gim->darray[0]->coordsys[0]==NULL)) { cout<<"null"<<endl; }
	  else { free(gim->darray[0]->coordsys[0]); } //->dataspace = cs;
	  //   cout<<"write gifti4.75"<<endl;
		
	  gim->darray[0]->numCS=0;
	  ////cout<<"free stuf 2"<<endl;
	  //   cout<<"write gifti5"<<endl;
		
	  gifti_free_image(gim);
		
	  return 0;
	}

  template int writeGIFTI<float, unsigned int>( const fslSurface<float, unsigned int> & surf, const std::string & filename, int enc);
	
	
}
