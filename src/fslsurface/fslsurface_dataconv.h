//
//  fslsurface_dataconv.h
//  fslsurface
//
//  Created by Brian Patenaude on 5/30/11.
//  Copyright 2011 University of Oxford. All rights reserved.
//
#ifndef FSLSURFACE_DATACONV_H
#define FSLSURFACE_DATACONV_H

//FSL includes
#include "meshclass/meshclass.h"

//STL INCLUDES
#include <string>

//3rd party includes
#include "newmatap.h"

//fslsurface includes
//#include "fslsurface.h"

//using namespace fslsurface_name;
namespace fslsurface_name {
    template<class T, class T2>
	class fslSurface;
    
    template<class T, class T2>
    NEWMAT::ReturnMatrix getScalars( fslSurface<T,T2>* surf, unsigned int index=0 );
    
    template<class T, class T2>
    void convertTo_Mesh( mesh::Mesh * dest_mesh, fslsurface_name::fslSurface<T,T2>* surf );
    
    template<class T, class T2>
    void convertTo_fslSurface( fslsurface_name::fslSurface<T,T2>* dest_surf, mesh::Mesh * mesh );

    template<class T, class T2,class T3>
    void convertTo_fslSurface( fslsurface_name::fslSurface<T,T2>* dest_surf, mesh::Mesh * mesh , const NEWIMAGE::volume<T3> &image_ref);

}
#endif
