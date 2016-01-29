//
//  fslsurface_dataconv.cpp
//  fslsurface
//
//  Created by Brian Patenaude on 5/30/11.
//  Copyright 2011 University of Oxford. All rights reserved.
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


#include "fslsurface_dataconv.h"

//3rd party includes
#include "newmatap.h"


//fslsurface include
#include "fslsurface.h"
#include "fslsurfacefns.h"

#include "fslsurface_structs.h"

//STL includes
#include <vector>

using namespace std;
using namespace NEWMAT;
using namespace mesh;
//MAY WANT TO REMOVE DEPENDENCY ON THIS LIBRARY IN THE FUTURE, we'll see

//STL
using namespace std;


namespace fslsurface_name {
    
    template<class T, class T2>
    ReturnMatrix getScalars( fslSurface<T,T2> & surf, unsigned int index )
    {
        Matrix Scalars(surf->getNumberOfVertices(),1); 
        unsigned int count = 0;
        for (typename vector<T>::iterator i= surf->const_scbegin(index); i!= surf->const_scend();++i,++count)
        {
            Scalars.element(count, 0) = *i;
        }

        Scalars.Release();
        return Scalars;
    }
    
    template<class T,class T2>
    void convertTo_Mesh( Mesh * mesh, fslSurface<T,T2>* surf )
    {
        //create mesh if it doesn't exists
        if (mesh == NULL)
            mesh = new Mesh();
        
        //clear anything that might be there
        mesh->_points.clear();
        mesh->_triangles.clear();
        unsigned int count=0;
        for ( typename std::vector< fslsurface_name::vertex<T> >::const_iterator  i= surf->const_vbegin(); i!= surf->const_vbegin();++i,++count)
        {
            mesh->_points.push_back( new Mpoint(i->x, i->y, i->z, count));
        }
        
        
        for ( typename vector<T2>::iterator i = surf->const_facesbegin(); i != surf->const_facesend(); ++i)
        {
            Triangle * t = new Triangle(get_point(*i), get_point(*(i+1)), get_point(*(i+2)));
            mesh->_triangles.push_back(t);
            
            //            mesh->_triangles.push_back(new Triangle());
        }
        
    }
    
    template<class T,class T2>
    void convertTo_fslSurface( fslSurface<T,T2>& surf, const Mesh & mesh )
    {
        
        //create mesh if it doesn't exists
      //  if (surf == NULL)
       //     surf = new fslSurface<T,T2>();
        vector<T> vertices(mesh.nvertices()*3);
        vector<T2> triangles(mesh._triangles.size()*3);
        
        typename vector<T>::iterator i_v = vertices.begin();
        for (vector<Mpoint *>::const_iterator i = mesh._points.begin(); i != mesh._points.end();++i,++i_v)
        {
            *i_v = static_cast<T>((*i)->get_coord().X);
            ++i_v;
            *i_v = static_cast<T>((*i)->get_coord().Y);
            ++i_v;
            *i_v = static_cast<T>((*i)->get_coord().Z);
            //cout<<"add vertex "<<
        }
        surf.setVertices(vertices);
        
        typename vector<T2>::iterator i_t = triangles.begin();
        for (list<Triangle* >::const_iterator i = mesh._triangles.begin(); i != mesh._triangles.end();++i,++i_t)
        {
            *i_t = static_cast<T2>(  (*i)->get_vertice(0)->get_no() );
            ++i_t;
            *i_t = static_cast<T2>(  (*i)->get_vertice(1)->get_no() );
            ++i_t;
            *i_t = static_cast<T2>(  (*i)->get_vertice(2)->get_no() );
        }
        surf.setFaces(triangles);
            
        //addin coordinate system
        vector<double> xfm(16,0);
        xfm[0]=xfm[5]=xfm[10]=xfm[15]=1.0;
        
        surf.addCoordSystem(xfm,"UNKNOWN" );

        
    }
    template void convertTo_fslSurface<float,unsigned int>( fslSurface<float,unsigned int>& surf, const Mesh & mesh );

    
    template<class T, class T2,class T3>
    void convertTo_fslSurface( fslsurface_name::fslSurface<T,T2>& dest_surf, const mesh::Mesh & mesh , const NEWIMAGE::volume<T3> & image_ref)
    {
        convertTo_fslSurface( dest_surf, mesh );
        dest_surf.clearCoordSystems();
        //assign appropriate coordinate system
        //        read_volume_hdr_only(*im, *i_imnames );
        Matrix m_vox2mm = image_ref.newimagevox2mm_mat();
        
        Matrix m_mm2vox(4,4);
        m_mm2vox=0;
        m_mm2vox.element(0, 0) = 1.0/image_ref.xdim();
        m_mm2vox.element(1, 1) = 1.0/image_ref.ydim();
        m_mm2vox.element(2, 2) = 1.0/image_ref.zdim();
        m_mm2vox.element(3, 3) = 1.0;
        
        
        Matrix xfm_nii = m_vox2mm * m_mm2vox;
        vector<float> fmat;
        for (int i = 0 ; i < 4; ++i)
            for (int j = 0 ; j < 4; ++j)
                fmat.push_back(xfm_nii.element(i, j));
        apply_xfm(dest_surf, fmat);

        vector< float > csys_eye(16,0);
        csys_eye[0] = csys_eye[5]= csys_eye[10]= csys_eye[15]=1;

        //	//cout<<"add coord system "<<i_xfms->at(0)<<endl;
        dest_surf.addCoordSystem( csys_eye, "NIFTI_XFORM_SCANNER_ANAT");
        //--------------------------back to mni------------------------------
        //cout<<"setmni "<<endl;
        xfm_nii = xfm_nii.i();
        fmat.clear();
        for (int i = 0 ; i < 4; ++i)
            for (int j = 0 ; j < 4; ++j)
                fmat.push_back(xfm_nii.element(i, j));
        
        dest_surf.addCoordSystem( fmat, "FSL_XFORM_SCALEDMM");
        
        
    }

    template void convertTo_fslSurface<float,unsigned int,float>( fslsurface_name::fslSurface<float,unsigned int>& dest_surf, const mesh::Mesh & mesh , const NEWIMAGE::volume<float> & image_ref);


}
