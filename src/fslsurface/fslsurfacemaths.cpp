//
//  main.cpp
//  fslsurfacemaths
//
//  Created by Brian Patenaude on 3/27/11.
//  Copyright 2011 University of Oxford. All rights reserved.
//
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
//FSL includes
#include <newimage/newimageall.h>
//#include <fslsurface/fslsurface.h>
//#include <fslsurface/fslsurfacefns.h>
//#include <fslsurface/fslsurfaceio.h>
//#include <fslsurface/fslsurface_first.h>
#include <fslsurface.h>
#include <fslsurfacefns.h>
#include <fslsurfaceio.h>
#include <fslsurface_first.h>
#ifdef __linux
#else
#include <AGL/agl.h>
#endif
#include <iostream>
#include <sstream>
using namespace fslsurface_name;
using namespace NEWIMAGE;
using namespace std;

float string2float(const string & snum ){
    stringstream ss;
    ss<<snum;
    float num;
    ss>>num;
    return num;
}

void Usage(){
    cout<<endl<<endl;
    cout<<"--------------------------- Usage ----------------------------------------------"<<endl;
    cout<<endl;
    cout<<" fslsurfacemaths <first_input_surface> [operations and inputs] <output> \n or \n fslsurfacemaths -reconFromBvars fike.bvars [operations and inputs] <output> \n\n"<<endl;
    cout<<"Operations: (some inputs can be either surfaces or numbers) \n\n "<<endl;
    cout<<"    -readGIFTI <surface.gii> read GIFTI files, will obey GIFTI IO rules for append and overwriting."<<endl;
    cout<<"    -add <surface> : Add surface vertices to existing vertices."<<endl;
    cout<<"    -sub <surface> : Subtract surface vertices to existing vertices."<<endl;
    cout<<"    -mul <scalar>  : Multiply surface vertices by a scalar."<<endl;
    cout<<"    -div <surface> : Add surface vertices to existing vertices."<<endl;
    cout<<"    -bin <scalars index>  : use (scalars>0) to binarise"<<endl;
    cout<<"    -sc_sub <index> <value>: subtract value from scalars."<<endl;
    cout<<"    -sc_mul <index> <value>: multiply scalar by value."<<endl;
    cout<<"    -sc_mul_sc <index> <surf2> <index2>: multiple scalars of 2 surfaces."<<endl;
    cout<<"    -sc_sub_sc <index> <surf2> <index2> <name_of_new>: subtract scalars of 2 surfaces."<<endl;

    cout<<"    -sc_bin <index> : Add surface vertices to existing vertices."<<endl;

    cout<<"    -sc_thr <thresh> : Thresh scalars at index zero."<<endl;
    cout<<"    -sc_uthr <thresh> : Add surface vertices to existing vertices."<<endl;

    cout<<"    -cluster <index> <threshold> : Cluster scalars at index with threshold."<<endl;
    cout<<"     -sc_smooth <index> : Smooth scalars at index."<<endl;
    cout<<"     -sc_smooth_gauss_geodesic <index> <variance(mm)> : Smoot scalars at index."<<endl;
    cout<<"     -sc_smooth_gauss_geodesic_4D <variance(mm)> : Smoot scalars at index."<<endl;

    cout<<"     --projVecOntoNormals <index> : project vector onto normals"<<endl;
    cout<<"    -applyxfm <xfm.mat> : Apply a linear XFM matrix to the vertices and replace existing vertices (same format as FLIRT)."<<endl;
    cout<<"    -reg_lq <surface> <dof> <output_name_xfm.mat> : Calculate a linear XFM matrix using Least-Squares (vertex differences). XFM'd vertices replace original. Output matrix stored in 'output_name_xfm' (same format as FLIRT)."<<endl;
    cout<<"     -copy_CoordSystems <surf> : Copy coordinate systems for a reference surface. Typically used along side -applyxfm"<<endl;
    cout<<"     -normals : Calculate normals and store in vector field."<<endl;
    cout<<"     -vmag <vector_index> <scalar_name>: Compute vector magnitude and store in scalars"<<endl; 
    cout<<"    -copy_verts <surface.gii> : Copy the vertices from the following surface."<<endl;
    cout<<"    -copy_verts2vec : inserts the vertices as vector into the vector field (no input)."<<endl;
    cout<<"    -copy_scalar2image3D <scalar_index> <output_image> : copies vector into a Nx1x1 NIFTI image."<<endl;
    cout<<"    -copy_scalars2image4D <output_image> : copies vector into a Nx1x1xNt NIFTI image."<<endl;

    cout<<"    -copy_image3D2scalar <image_name> <scalar_name> : insert a Nx1x1 NIFTI into the vector field. Also specify the names of the scalars"<<endl;
    cout<<"    -copy_image4D2scalar <image_name> <scalar_name> : insert a Nx1xM NIFTI into the vector field. Also specify the names of the scalars"<<endl;

    cout<<"    -copy_vec2image3D <vec_index> <output_image> : copies vector into a Nx3x1 NIFTI image."<<endl;
    cout<<"    -copy_image3D2vec <image_name> <vector_name> <vec_index> <xyz>: insert a Nx3x1 NIFTI into the vector field. The last arguments specifies which coomponents of the vector to copy over"<<endl;
    cout<<"      --sampleTimeSeries <4D image> <scalar_name>"<<endl;
    cout<<"     -vertexMVGLM : <list_of_surfaces> <design_matrix>"<<endl;
    cout<<"     -FtoP : convert F values to scalar values"<<endl;
    cout<<"     -fillMesh <label> <ref_im> <save_name> : Draw surface."<<endl;

    cout<<"    -reconFromBvars <file.bvars> : Reconstructs surfaces from .bvars files output (only first surface)."<<endl;
    cout<<"    -reconAllFromBvarsAndSave <file.bvars> <output_basename> : Reconstructs all surfaces from .bvars and save"<<endl;
    cout<<"    -marchingCubes <image> <threshold> <label> <mode> : run marching cubes."<<endl;

    cout<<"     -xfmToFSLScaledMM : Transform to scaled mm coordinates."<<endl;
    cout<<endl;
    cout<<"------------------------------------------------------------------------------"<<endl;

}

int main (int argc, char * argv[])
{
    //Assume we're working with vertices will allow to change mode to scalars
    
    // insert code here...
    std::cout << "Welcome to FSL surface maths !\n";
    
    int i_arg=1;
    if (argc<2)
    {
        cout<<"Too few arguments."<<endl;
        Usage();
        return 1;
    }
    
    fslSurface<float,unsigned int> surf;// = new fslSurface<float,unsigned int>();
    if ( ( string(argv[i_arg]) != "-reconFromBvars") && (string(argv[i_arg]) != "-reconAllFromBvarsAndSave") && ( string(argv[i_arg]) != "-vertexMVGLM" )  && ( string(argv[i_arg]) != "-marchingCubes" ) )
    {
        read_surface(surf,argv[i_arg]);
    
    
        cout<<"done reading surface "<<endl;
        ++i_arg;
    }
    cout<<"iarg "<<i_arg<<" "<<argc<<endl;
    while (i_arg < (argc-1)) {
        cout<<"i_arg "<<i_arg<<" "<<argc<<" "<<argv[i_arg]<<endl;
        string command = string(argv[i_arg]);
        if (command == "-readGIFTI"){
            readGIFTI(surf, argv[i_arg+1]);
            i_arg += 2;
        }else if (command == "-sub"){
            fslSurface<float,unsigned int> surf2;
            cout<<"read surface "<<endl;
            read_surface(surf2, argv[i_arg+1]);
            cout<<"read surface2 "<<endl;

            surf-=surf2;
            cout<<"done sub surface "<<endl;

            i_arg += 2;
        }else if (command == "-sc_sub"){

            surf.subtractScalars(atoi(argv[i_arg+1]),string2float(argv[i_arg+2]));
            cout<<"done sub surface "<<endl;
            
            i_arg += 3;
        }else if (command == "-sc_mul"){
            
            surf.multiplyScalars(atoi(argv[i_arg+1]),string2float(argv[i_arg+2]));
            cout<<"done mul surface "<<endl;
            
            i_arg += 3;
        }else if (command == "-sc_mul_sc"){
            fslSurface<float,unsigned int> surf2;
            cout<<"read surface "<<endl;
            read_surface(surf2, argv[i_arg+2]);
            cout<<"multiply sclars"<<endl;
            multiplyScalarsByScalars(surf,atoi(argv[i_arg+1]),surf2,atoi(argv[i_arg+3]));
            cout<<"done sc mul surface "<<argv[i_arg+3]<<endl;
            
            i_arg += 4;
            cout<<"done sc mul surface2 "<<argv[i_arg]<<endl;

        }else if (command == "-sc_sub_sc"){
            fslSurface<float,unsigned int> surf2;
            cout<<"read surface "<<endl;
            read_surface(surf2, argv[i_arg+2]);
            cout<<"subtract sclars"<<endl;
            string sc_name=argv[i_arg+4];

            subtractScalarsFromScalars(surf,atoi(argv[i_arg+1]),surf2,atoi(argv[i_arg+3]),sc_name);
            cout<<"done sc sub surface "<<argv[i_arg+3]<<endl;
            i_arg += 5;
            cout<<"done sc sub surface2 "<<argv[i_arg]<<endl;
            
        }else if (command == "-add"){
            fslSurface<float,unsigned int> surf2;// = new fslSurface<float,unsigned int>();
            read_surface(surf2, argv[i_arg+1]);
            
            surf+=surf2;
            i_arg += 2;

        }else if (command == "-mul"){
            surf *= string2float(argv[i_arg+1]);
            i_arg += 2;
            cout<<"done multiply"<<endl;
            
        }else if (command == "-div"){
            surf /= string2float(argv[i_arg+1]);
            i_arg += 2;

        }else if (command == "-sc_bin"){
///            surf /= string2float(argv[i_arg+1]);
            surf.binariseScalars(atoi(argv[i_arg+1]),0.0);
            i_arg += 2;
            
        }else if (command == "-sc_thr"){
            float thresh = string2float(argv[i_arg+1]);
            surf.thresholdScalars(0,thresh);
            i_arg += 2;
            
        }else if (command == "-sc_uthr"){
            float thresh = string2float(argv[i_arg+1]);
            surf.upperThresholdScalars(0,thresh);
            i_arg += 2;
            
        }else if (command == "-cluster"){
            cout<<"foudn cluster "<<atoi(argv[i_arg+1])<<" "<<atof(argv[i_arg+2])<<endl;
            fslsurface_name::cluster<float,unsigned int>(surf,atoi(argv[i_arg+1]),atof(argv[i_arg+2]));
            i_arg += 3;
            
        }else if ( command == "-sc_smooth" ){
                
            sc_smooth_mean_neighbour<float,unsigned int>(surf,atoi(argv[i_arg+1]));
            i_arg += 2;

        }else if ( command == "-sc_smooth_gauss_geodesic_4D" ){
            
            sc_smooth_gaussian_geodesic<float,unsigned int>(surf,0,string2float(argv[i_arg+1]),string2float(argv[i_arg+1])*4,true);
            i_arg += 2;
            
        }else if (command == "-sc_smooth_gauss_geodesic" ){
            sc_smooth_gaussian_geodesic<float,unsigned int>(surf,atoi(argv[i_arg+1]),string2float(argv[i_arg+2]),string2float(argv[i_arg+2])*4);
                             i_arg+=3;
            cout<<"done smooth"<<endl;
        }else if (command == "-normals" ){
            surf.calculateNormals(true,false);
            surf.copyNormalsToVectors(0);
            ++i_arg;
        }else if (command == "-projVecOntoNormals"){
            
            projectVectorsOntoNormals(surf,atoi(argv[i_arg+1]));
            i_arg += 2;
            
        }else if (command == "-applyxfm"){
       
            apply_xfm<float, unsigned int>(surf,string(argv[i_arg+1])); 
            i_arg += 2;
            
        }else if (command == "-xfmToFSLScaledMM"){
            // surf /= string2float(argv[i_arg+1]);
            vector<double> fmat = surf.getCoordinateSystem("FSL_XFORM_SCALEDMM");
            cout<<"got coordinate systrem "<<endl;
            if (fmat.size() != 16)
            {
                cerr<<"Invalid XFM size.Most likely FSL_XFORM_SCALEDMM does not exist."<<endl;
                exit (EXIT_FAILURE);
            }
           cout<<"xfmToFSLScaledMM Matrix"<<endl;  
           // vector<double> fmati = invertMatrix(fmat);
          //  for (vector<double>::iterator i = fmati.begin(); i!= fmati.end();++i)
          //      cout<<*i<<" ";
         //   cout<<endl;
            // print<double>(fmati,"inversematrix");
            apply_xfm<float, unsigned int>(surf,fmat); 
           // apply_xfm<float, unsigned int>(&surf,string(argv[i_arg+1])); 
            //apply_xfm<float, unsigned int>(&surf,fmati); 
            
            cout<<"done apply"<<endl;
            ++i_arg;// += 2;
            
        }else if (command == "-applyxfmFlirt"){
           // surf /= string2float(argv[i_arg+1]);
            vector<double> fmat = surf.getCoordinateSystem("FSL_XFORM_SCALEDMM");
            cout<<"got coordinate systrem "<<endl;
            if (fmat.size() != 16)
            {
                cerr<<"Invalid XFM size.Most likely FSL_XFORM_SCALEDMM does not exist."<<endl;
                exit (EXIT_FAILURE);
            }
            cout<<"invert Matrix"<<endl;  
            volume<float> imref;
            read_volume_hdr_only(imref, argv[i_arg+2]);
            
            vector<float> xfm = getFSLtoNIFTIxfm(imref);

            
           // vector<double> fmati = invertMatrix(fmat);
            //for (vector<double>::iterator i = fmati.begin(); i!= fmati.end();++i)
             //   cout<<*i<<" ";
            //cout<<endl;
            // print<double>(fmati,"inversematrix");
            apply_xfm<float, unsigned int>(surf,fmat); 
            apply_xfm<float, unsigned int>(surf,string(argv[i_arg+1])); 
            apply_xfm<float, unsigned int>(surf,xfm); 

            
            i_arg += 3;
            
        }else if (command == "-reg_lq"){
            cout<<"reglq"<<endl;
            // surf /= string2float(argv[i_arg+1]);
            fslSurface<float,unsigned int> surf2;// = new fslSurface<float,unsigned int>();
            cout<<"read surface "<<argv[i_arg+1]<<endl;
            read_surface(surf2, argv[i_arg+1]);
            unsigned int dof = atoi(argv[i_arg+2]);

            cout<<"do mesh reg"<<endl;
            //do registration
            vector<float> xfm = meshRegLeastSq( surf, surf2, dof );
            cout<<"xfm4x4 "<<xfm.size()<<endl;
            //save linear xfm
            ofstream fmat_out(argv[i_arg+3]);
            int col=1;
            for (vector<float>::iterator i = xfm.begin(); i != xfm.end(); ++i,++col)
            {
                if ( (col % 4) ==0)
                {
                    fmat_out<<*i<<endl;
                }else {
                    fmat_out<<*i<<" ";
                }
            }
            //applyxfm to vertices
            apply_xfm<float, unsigned int,float>(surf,xfm); 
            surf.copyCoordSystems(surf2);
            i_arg += 4;
            
        }else if ( command == "-copy_CoordSystems" )
        {
            fslSurface<float,unsigned int> surf2;
            cout<<"read surface "<<argv[i_arg+1]<<endl;
            read_surface(surf2, argv[i_arg+1]);
            surf.copyCoordSystems(surf2);

            i_arg += 2;
        }else if (command == "-copy_verts"){
            fslSurface<float,unsigned int> surf2;// = new fslSurface<float,unsigned int>();
            read_surface(surf2, argv[i_arg+1]);
            surf.copyVertices(surf2);
            //surf+=surf2;
            i_arg += 2;
            
        }else if (command == "-copy_verts2vec"){
            surf.copyVerticesToVectors();
            ++i_arg;
        }else if (command == "-copy_scalar2image3D"){
            volume<float> image(surf.getNumberOfVertices(),1,1);
            unsigned int scind = atoi(argv[i_arg+1]);

            unsigned int index=0;
            for (std::vector< float >::const_iterator i_sc = surf.const_scbegin(scind); i_sc != surf.const_scend(scind); ++i_sc,++index)
            {
                image.value(index,0,0) = *i_sc;
            }

            save_volume(image, argv[i_arg+2]);
            i_arg += 3;
            if (i_arg == argc)
                return 0;
            
        }else if (command == "-copy_scalars2image4D"){
            unsigned int Nt=surf.getNumberOfScalarData();
            volume4D<float> image(surf.getNumberOfVertices(),1,1,Nt);
//            unsigned int scind = atoi(argv[i_arg+1]);
  
            for ( unsigned int scind = 0 ; scind < Nt;++scind)
            {
                unsigned int index=0;
                for (std::vector< float >::const_iterator i_sc = surf.const_scbegin(scind); i_sc != surf.const_scend(scind); ++i_sc,++index)
                {
                    image[scind].value(index,0,0) = *i_sc;
                }
            }
            save_volume4D(image, argv[i_arg+1]);
            i_arg += 2;
            if (i_arg == argc)
                return 0;
            
        }else if (command == "-copy_vec2image3D"){
            cout<<"FOUND copy_vec2image3D "<<endl;
            unsigned int N = surf.getNumberOfVertices();
            volume<float> image(N,3,1);
            unsigned int vecind = atoi(argv[i_arg+1]);
            unsigned int index=0;
            for (vector<float>::const_iterator i_vec = surf.const_vecbegin(vecind); i_vec != surf.const_vecend(vecind); ++i_vec,++index)
            {
                image.value(index,0,0)=*i_vec;
                ++i_vec;
                image.value(index,1,0)=*i_vec;
                ++i_vec;
                image.value(index,2,0)=*i_vec;

                
            }
            cout<<"save image "<<argv[i_arg+2]<<endl;
            save_volume(image, argv[i_arg+2]);
            
            i_arg += 3;
            if (i_arg == argc)
                return 0;
            
        }else if (command == "-vmag"){
            cout<<"FOUND vmag "<<endl;
            unsigned int N = surf.getNumberOfVertices();
            unsigned int vecind = atoi(argv[i_arg+1]);
            vector<float> sc(N);
            vector<float>::iterator i_sc = sc.begin();
            for (vector<float>::const_iterator i_vec = surf.const_vecbegin(vecind); i_vec != surf.const_vecend(vecind); ++i_vec,++i_sc)
            {
                float x = *i_vec;
                ++i_vec;
                float y = *i_vec;
                ++i_vec;
                float z = *i_vec;
                *i_sc = sqrtf( x*x + y*y + z*z );
                
            }
            surf.insertScalars(sc, 0, argv[i_arg+2]);
         //   surf.addScalars(sc, argv[i_arg+2]);
            i_arg += 3;
            
        }else if (command == "-copy_image3D2scalar"){
            unsigned int N = surf.getNumberOfVertices();
            volume<float> image;
            read_volume(image,argv[i_arg+1]);
            if ((static_cast<unsigned int>(image.xsize())!=N) || (static_cast<unsigned int>(image.ysize())!=1) || (static_cast<unsigned int>(image.zsize())!=1))
            {
                cerr<<"Error : copy_image3D2scalar : Incompatible image dimensions "<<endl;
                exit (EXIT_FAILURE);
            }
            

            vector<float> sc(N);
            vector<float>::iterator i_sc = sc.begin();
            for ( unsigned int x = 0 ; x < N; ++x,++i_sc)
            {
                *i_sc = image.value(x,0,0);
                cout<<"sc "<<*i_sc<<endl;

            }
            cout<<"insert the scalars"<<endl;
            surf.insertScalars(sc,0,argv[i_arg+2]);
            i_arg += 3;

        }else if (command == "-copy_image4D2scalar"){
            unsigned int N = surf.getNumberOfVertices();
            volume4D<float> image;
            read_volume4D(image,argv[i_arg+1]);
            if ((static_cast<unsigned int>(image.xsize())!=N) || (static_cast<unsigned int>(image.ysize())!=1) || (static_cast<unsigned int>(image.zsize())!=1))
            {
                cerr<<"Error : copy_image4D2scalar : Incompatible image dimensions "<<endl;
                exit (EXIT_FAILURE);
            }
            insertScalars<float,unsigned int,float>(surf,image,argv[i_arg+2]);
            
          /*
           vector<float> sc(N);
            vector<float>::iterator i_sc = sc.begin();
            for ( unsigned int x = 0 ; x < N; ++x,++i_sc)
            {
                *i_sc = image.value(x,0,0);
                cout<<"sc "<<*i_sc<<endl;
                
            }
            cout<<"insert the scalars"<<endl;
            surf.insertScalars(sc,0,argv[i_arg+2]);
           */
            i_arg += 3;
            
        }else if (command == "-copy_image3D2vec"){
            cout<<"FOUND copy_image3D2vec "<<endl;
            unsigned int N = surf.getNumberOfVertices();
            volume<float> image;
            read_volume(image,argv[i_arg+1]);
            cout<<"N "<<N<<" "<<image.xsize()<<" "<<image.ysize()<<" "<<image.zsize()<<endl;
            
            if ((static_cast<unsigned int>(image.xsize())!=N) || (static_cast<unsigned int>(image.ysize()) !=3) || (static_cast<unsigned int>(image.zsize())!=1))
                        {
                            cerr<<"Error : copy_image3D2vec : Incompatible image dimensions "<<endl;
                            exit (EXIT_FAILURE);
                        }
            
                        //(N,3,1);
            vector<float> vec(N*3);
            
            
            string comp = string(argv[i_arg+4]);
            bool cpx = false;
            bool cpy = false;
            bool cpz = false;
            
            for ( unsigned int i =0 ; i<comp.length();++i)
            {
                string c = comp.substr(i,1);
                cout<<"comp "<<c<<endl;
                if (c=="x")
                    cpx=true;
                else if (c=="y")
                    cpy=true;
                else if (c=="z")
                    cpz=true;
            }
            
            
            vector<float>::iterator i=vec.begin();
            for ( unsigned int x = 0 ; x < N; ++x)
            {
                if (cpx)
                    *i=image.value(x,0,0);
                else
                    *i=0;
                ++i;
                if (cpy)                    
                    *i=image.value(x,1,0);
                else
                    *i=0;
                ++i;
                if (cpz)
                    *i=image.value(x,2,0);
                else
                    *i=0;
                ++i;
                
            }
            unsigned int vecind = atoi(argv[i_arg+3]);

            surf.addVectors(vec, argv[i_arg+2], vecind);       
            i_arg += 5;
            
        }else if ( command == "-sampleTimeSeries" ) {
            volume4D<float> image;
            read_volume4D(image,argv[i_arg+1]);
            
            surf.clearScalars();
            sampleTimeSeries(surf , image, argv[i_arg+2] );
               i_arg += 3;

        }else if ( command == "-writeTimeSeries" ){
            volume4D<float> image;
            read_volume4D(image,argv[i_arg+1]);
            
            writeTimeSeries( surf , image , argv[i_arg+2]);
            i_arg += 3;

            
        }else if ( command == "-marchingCubes"  ){
            volume<float> image;
            read_volume(image,argv[i_arg+1]);
            image_dims dims(image.xsize(), image.ysize(), image.zsize(), image.xdim(), image.ydim(), image.zdim());
            marchingCubes<float,unsigned int, float>( surf, image.fbegin(),  dims, atof(argv[i_arg+2]), atof(argv[i_arg+3]), static_cast<MarchingCubesMode>(atof(argv[i_arg+4])));
            i_arg += 5;

        }else if (command == "-reconFromBvars"){
            cout<<"do recon "<<endl;
            char* fsldir = getenv("FSLDIR");
            if (fsldir == NULL)
            {
                cerr<<"FSLDIR has not been set. "<<endl;
                exit(EXIT_FAILURE); 
            }
            //file.bvars,mni_template.nii.gz
           // string mni = string(fsldir)+"/data/standard/MNI152_T1_1mm";
            reconSurface_from_bvars( surf, string(argv[i_arg+1]));
            i_arg += 2;
            cout<<"done recon"<<endl;

        }else if (command == "-reconAllFromBvarsAndSave"){
            cout<<"do recon+save "<<argc<<" "<<i_arg<<endl;
            char* fsldir = getenv("FSLDIR");
            if (fsldir == NULL)
            {
                cerr<<"FSLDIR has not been set. "<<endl;
                exit(EXIT_FAILURE); 
            }
            //file.bvars,mni_template.nii.gz
            // string mni = string(fsldir)+"/data/standard/MNI152_T1_1mm";
            cout<<"recon "<< string(argv[i_arg+1])<<endl;
            cout<<" "<<string(argv[i_arg+2])<<endl;
            reconAllSurfacesAndSave( string(argv[i_arg+1]), string(argv[i_arg+2]));
            i_arg += 3;
            cout<<"done recon"<<endl;
            return 0;
        }else if ( command == "-FtoP" ){
            cout<<"do FtoP"<<endl;
            vector<int> fdofs=surf.getNonVertIntScalars("FStatDOFs");
            
            cout<<"got dofs "<<fdofs.size()<<endl;
            int df1=fdofs[0];
            int df2=fdofs[1];
            cout << "FStat DOFs = " << df1 << " and " << df2 << endl;
            
            surf.insertScalars(FtoP(surf.getScalars("F-stats_EV_1"),df1,df2),0,"FtoP");
            
            ++i_arg;

            			
		}else if ( command == "-vertexMVGLM" ){
            fslSurface<float,unsigned int> *surf2 = new fslSurface<float,unsigned int>();
            //list of surfaces, deign matrx name, save vertices, norm name
            cout<<"dovertex glm"<<endl;
            vertexMVglm( *surf2, string(argv[i_arg+1]),  string(argv[i_arg+2]), "",string(argv[i_arg+3]));
            cout<<"done glm"<<endl;
            surf=*surf2;
            delete surf2;
            i_arg += 4;

        }else if ( command == "-warp" ){            
            volume4D<float> im_warp;
            cout<<"warp file "<<argv[i_arg+1]<<endl;
            read_volume4D(im_warp,argv[i_arg+1]);
            apply_warp(surf, im_warp);
            i_arg += 2;
        
        }else if ( command == "-deform" ){
            cout<<"deform"<<endl;
        //    int win = glutInitAndCreateWindow(&argc,argv,argv[i_arg+1]);
      //      glutDisplayFunc(glutDisp);
          //  glutIdleFunc(glutIdle);
           // glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            //glutKeyboardFunc(keyb);
            //glutRender(surf);
          //  glBegin(GL_TRIANGLES);
            //glVertex3f(-0.5,-0.5,0.0);
            //glVertex3f(0.5,0.0,0.0);
            //glVertex3f(0.0,0.5,0.0);
            //glEnd();
//            glutMainLoop();
          //  glutDestroyWindow(win);
          //  glutSwapBuffers();
            cout<<"done deform "<<endl;
            i_arg+=2;
        }else if ( command == "-fillMesh" ){

            
            int lb = atoi(argv[i_arg+2]);
            
            
            vector< vec3<float> > test ;
            test.push_back(vec3<float>(4,5,6));
            test.push_back(vec3<float>(1,3,2));
       //     float* fs=static_cast<float*>(&(test[0].x));
        //    for (int i =0 ; i< 6; ++i)
          //      cout<<"fs "<<fs[i]<<" "<<lb<<endl;
            
            volume<short> mask;
            
            read_volume(mask,argv[i_arg+1]);
            mask=0;
            mask=fillMesh( surf, mask, lb);
            save_volume(mask,argv[i_arg+3]);
            
            i_arg += 4;
            
        }else{
            
            cerr<<"Unrecognized command "<<command<<endl;
            Usage();
            exit(EXIT_FAILURE);
        }
    }
    
    
   cout<<"do write "<<i_arg<<" "<<argc<<" "<<argv[i_arg]<<" "<<surf.getNumberOfCoordinateSystems()<<endl;
    if (surf.getNumberOfCoordinateSystems() == 0 )
    {
        vector<double> identity(16,0);
        identity[0] = identity[5] = identity[10] = identity[15] = 1.0; 
        surf.addCoordSystem(identity, "NIFTI_XFORM_UNKNOWN");
    }
    
    if (i_arg > argc)
    {
        cerr<<"Did not specify output name"<<endl;
        exit (EXIT_FAILURE);
    }
    cout<<"write gifti"<<endl;
    writeGIFTI(surf, string(argv[i_arg]));
    cout<<"done writing gifti"<<endl;
    //  delete surf;
    return 0;
}

