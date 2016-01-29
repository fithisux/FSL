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
#include <fslsurface_first.h>

//FSL INCLUDES 
#include <first_lib/first_newmat_vec.h>
#include <fslvtkio/fslvtkio.h>
#include <fslsurface.h>
#include <fslsurfacefns.h>
#include <fslsurfaceio.h>
#include "newmatap.h"
#include "libprob.h"

//FSL namespaces 
using namespace SHAPE_MODEL_NAME;
using namespace fslvtkio;
using namespace NEWIMAGE;

//MAY WANT TO REMOVE DEPENDENCY ON THIS LIBRARY IN THE FUTURE, we'll see
using namespace FIRST_LIB;
//STL
using namespace std;


namespace fslsurface_name {
	
    
    
    
    template<class T>
    void draw_segment(volume<short>& image, const vec3<T> & p1, const vec3<T> &  p2, const int & label)
    {
        vec3<float> dims = vec3<float>(image.xdim(),image.ydim(), image.zdim());
        float mininc = (dims.y < dims.z) ? dims.y : dims.z;
        mininc = ( ( dims.x < mininc ) ? dims.x : mininc ) * 0.5; 
      //  cout<<"dims "<<dims.x<<" "<<dims.y<<" "<<dims.z<<" "<<mininc<<endl;
        vec3<float> n = p1 - p2;
        double d = n.norm();
        n.normalize();
        //	double l = d*4;
        for (double i=0; i<=d; i+=mininc)
        {
        //    cout<<"draw "<<d<<" "<<i<<endl;
            vec3<float> p = p2 + n * i;
          //  cout<<(int) floor((p.x)/dims.x +.5)<<" "<<(int) floor((p.y)/dims.y +.5)<<" "<<(int) floor((p.z)/dims.z +.5)<<" "<<label<<endl;
            image((int) floor((p.x)/dims.x +.5),(int) floor((p.y)/dims.y +.5),(int) floor((p.z)/dims.z +.5)) = label;
            //cout<<image((int) floor((p.x)/dims.x +.5),(int) floor((p.y)/dims.y +.5),(int) floor((p.z)/dims.z +.5))<<endl;
        }
    }
    
    template void draw_segment<float>(volume<short>& image, const vec3<float> & p1, const vec3<float> &  p2, const int & label);

    
    template<class T, class T2>
    void draw_mesh(const fslSurface<T,T2> & surf, volume<short> & image, const int & lb)
    {
        vec3<float> dims = vec3<float>(image.xdim(),image.ydim(), image.zdim());
        
        //in new version of bet2
//        double mininc = min(xdim,min(ydim,zdim)) * .5;
        float mininc = (dims.y < dims.z) ? dims.y : dims.z;
        mininc = ( ( dims.x < mininc ) ? dims.x : mininc ) * 0.5;   
        typename vector< vertex<T> >::const_iterator i_v = surf.const_vbegin();
        for ( typename vector<T2>::const_iterator i_f = surf.const_facebegin(); i_f != surf.const_faceend(); i_f= i_f+3)
        {
            vec3<T> n = vec3<T>( (i_v+ *(i_f+1) )->x - (i_v+ *i_f)->x, \
                          (i_v+ *(i_f+1) )->y - (i_v+ *i_f)->y, \
                          (i_v+ *(i_f+1) )->z - (i_v+ *i_f)->z );
            float d = n.norm();
            n.normalize();
            for (float j=0; j<=d; j+=mininc)
            {
                
                vec3<T> p= n*j;
                p += vec3<T>( (i_v+ *i_f)->x , (i_v+ *i_f)->y , (i_v+ *i_f)->z );
                vec3<T> p2= vec3<T>( (i_v+ *(i_f+2))->x,(i_v+ *(i_f+2))->y,(i_v+ *(i_f+2))->z );
                
                //coord();
                draw_segment(image, p  ,p2,lb);

            }
           // cout<<"f "<<*i_f<<endl;
        }
      //  save_volume(image,"mesh2");

      //  volume<short> res = image;

       /* 
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
        */
        //return res;
    }
    
    template void draw_mesh<float,unsigned int>(const fslSurface<float,unsigned int> & surf, volume<short>& image, const int & label);

    
    template<class T, class T2>
    volume<short> fillMesh( const fslSurface<T,T2> & surf, const volume<short> & image, const int & label )  
    {
        vec3<float> dims = vec3<float>(image.xdim(),image.ydim(), image.zdim());
                
        volume<short> mask;
        copyconvert(image,mask);
        
        mask = 0;
        draw_mesh(surf, mask, label+100);
        volume<short> otl = mask;
        vector<T> bounds = surf.getBounds();
     //   cout<<"got bounds "<<bounds[0]<<" "<<bounds[1]<<" "<<bounds[2]<<" "<<bounds[3]<<" "<<bounds[4]<<" "<<bounds[5]<<endl;
        bounds[0]/=dims.x;
        bounds[1]/=dims.x;
        bounds[2]/=dims.y;
        bounds[3]/=dims.y;
        bounds[4]/=dims.z;
        bounds[5]/=dims.z;
        
        
        vector< vec3<int> > current;

        vec3<int> c(bounds[0]-2, bounds[2]-2, bounds[4]-2);
		
        current.push_back(c);
        
        mask.value(static_cast<int>(c.x),static_cast<int>(c.y),static_cast<int>(c.z)) = label;

        int fillCount=0;
        while (!current.empty())
        {
            vec3<int> pc = current.back();
            int x, y, z;
            x=(int) pc.x; y=(int) pc.y; z=(int) pc.z;
			
            current.pop_back();
            fillCount++;
			
            if (bounds[0]<=x-1 && mask.value(x-1, y, z)==0) {
                mask.value(x-1, y, z) = label;
                current.push_back(vec3<int>(x-1, y, z));
            }
            if (bounds[2]<=y-1 && mask.value(x, y-1, z)==0) {
                mask.value(x, y-1, z) = label;
                current.push_back(vec3<int>(x, y-1, z));
            }
            if (bounds[4]<=z-1 && mask.value(x, y, z-1)==0) {
                mask.value(x, y, z-1) = label;
                current.push_back(vec3<int>(x, y, z-1));
            }
            if (bounds[1]>=x+1 && mask.value(x+1, y, z)==0){
                mask.value(x+1, y, z) = label;
                current.push_back(vec3<int>(x+1, y, z));
            }
            if (bounds[3]>=y+1 && mask.value(x, y+1, z)==0){
                mask.value(x, y+1, z) = label;
                current.push_back(vec3<int>(x, y+1, z));
            }
            if (bounds[5]>=z+1 && mask.value(x, y, z+1)==0){
                mask.value(x, y, z+1) = label;
                current.push_back(vec3<int>(x, y, z+1)); 
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
        
       //save_volume(otl,"mesh");


        return otl;
    }
    
    template NEWIMAGE::volume<short> fillMesh<float,unsigned int>( const fslSurface<float,unsigned int> & surf, const NEWIMAGE::volume<short> & image, const int & label ) ;

    //template volume<short> fillMesh<float,unsigned int>( const fslSurface<float,unsigned int> & surf, volume<short> & image, const int & label );

    
    
	float MVGLM_fit(Matrix G, Matrix D, Matrix contrast, int& df1, int& df2){
		cout<<"enter MVGLMFIT "<<G.Nrows()<<" "<<G.Ncols()<<" "<<D.Nrows()<<" "<<D.Ncols()<<endl;
		//cout<<"contrast "<<contrast<<endl;
		// conversion to "normal" notation is:
		//  G -> X = design matrix (N_subj by N_evs)
		//  D -> Y = data matrix (N_subj by 3)
		
		//Calculate estimated values
		Matrix Yhat=G*(G.t()*G).i()*G.t()*D;
		//caluclate E covariance matrix
        cout<<"enter MVGLMFIT2 "<<endl;

		Matrix E=D-Yhat;
		E=E.t()*E;
        cout<<"enter MVGLMFIT3 "<<endl;

		//calculate H, the sum-square /cross square porduct for hypthosis test
		Matrix YhatH= G*contrast.t()*(contrast*G.t()*G*contrast.t()).i()*contrast*G.t()*D;
		Matrix H=D-YhatH;
		//nto effecient but easy to convert to other statistics
		H=H.t()*H-E;
	
		// Calculate Pillai's Trace
		
		int N=D.Nrows();//number of samples 
		int p=D.Ncols();//number of dimensions (3 for vertex coordinates)
		int N_R=G.Ncols();//number of regressors
		int N_C=contrast.Nrows(); //number of rows in contrast matrix
		
		int v_h = N_R - N_C;
		int v_e = N-N_R;	
        cout<<"mode"<<endl;
		int s=Min(p,v_h);
		float t = (abs(p-v_h)-1)/2.0;
		float u = (v_e-p-1)/2.0;
		df1=MISCMATHS::round(s*(2*t+s+1));
		df2=MISCMATHS::round(s*(2*u+s+1));
		
		float pillai=(H*(H+E).i()).Trace();
		float F=0;
		F=(pillai/(s-pillai))*(df2/df1);
		cout<<"leave MVGLM"<<endl;
      
		//if (verbose.value()){
			//cout<<"Pillai F "<<pillai<<" "<<F<<" "<<df1<<" "<<df2<<endl;
	//	}
		
		return F;
	}
	
	
	float MVGLM_fit(Matrix G, Matrix D, Matrix contrast)
	{
		//cout<<"enter MVGLMFIT2"<<endl;

		int df1, df2;
		float retval=MVGLM_fit(G,D,contrast,df1,df2);
		return retval;
	}
    
    vector<float> FtoP( const std::vector<float> & F, const int & dof1, const int & dof2)
    {
        cout<<"first FtoP "<<endl;
        vector<float> P(F.size());
  
        vector<float>::iterator i_P = P.begin();
        for (vector<float>::const_iterator i_F = F.begin(); i_F != F.end(); ++i_F,++i_P)
            *i_P = 1-MISCMATHS::fdtr(dof1,dof2,*i_F);
        
        return P;
	
    }

	
	shapeModel* loadAndCreateShapeModel( const string & modelname)
	{
		//read in model
		fslvtkIO* fmodel = new fslvtkIO(modelname,static_cast<fslvtkIO::DataType>(0));
		
		//get number of points
		const int Npts=fmodel->getPointsAsMatrix().Nrows();
		
		unsigned int M = static_cast<unsigned int>(fmodel->getField("numSubjects").element(0,0));
		int MaxModes=M;
		
		//Setting uo shape model and appearance model
		//load mean shape into vector
		vector<float> Smean;
		Matrix* Pts = new Matrix;
		*Pts=fmodel->getPointsAsMatrix();
		
		for (int i=0;i< Npts ; i++){
			Smean.push_back(Pts->element(i,0));
			Smean.push_back(Pts->element(i,1));
			Smean.push_back(Pts->element(i,2));
		}   
		Pts->ReleaseAndDelete();
		
		//read polygon data
		vector< vector<unsigned int > > polygons = first_newmat_vector::matrixToVector<unsigned int>(fmodel->getPolygons().t());
		
		//process shape modes and conditional intensity mean modes
		Matrix SmodesM;
		Matrix ImodesM;
		
		
		SmodesM=first_newmat_vector::unwrapMatrix(fmodel->getField("mode0"));
		ImodesM=first_newmat_vector::unwrapMatrix(fmodel->getField("Imode0"));
		
		//only keep max modes of variation 
		for (int i =1; i<MaxModes;i++)
		{
			stringstream ss;
			ss<<i;
			string mode;
			ss>>mode;
			SmodesM=SmodesM | first_newmat_vector::unwrapMatrix(fmodel->getField("mode"+mode));
			ImodesM=ImodesM | first_newmat_vector::unwrapMatrix(fmodel->getField("Imode"+mode));
			
		}
		
		
		
		
		vector< vector<float > > Smodes = first_newmat_vector::matrixToVector<float>(SmodesM);
		vector< vector<float > > Imodes = first_newmat_vector::matrixToVector<float>(ImodesM);
		ImodesM.Release();
		SmodesM.Release();
		
		
		//process rest of information, including intensity variance
		vector< vector<float > > Iprec = first_newmat_vector::matrixToVector<float>(fmodel->getField("iCondPrec0").t());
		vector<float > Errs =  first_newmat_vector::vectorToVector<float>(fmodel->getField("ErrPriors0"));
		vector<float > se =  first_newmat_vector::vectorToVector<float>(fmodel->getField("eigenValues"), MaxModes);
		vector<float > ie =  first_newmat_vector::vectorToVector<float>(fmodel->getField("iCondEigs0"));
		vector<float > Imean =  first_newmat_vector::vectorToVector<float>(first_newmat_vector::unwrapMatrix(fmodel->getField("Imean")));
		vector<int> labels =  first_newmat_vector::vectorToVector<int>(fmodel->getField("labels"));
		
		
		//have read in all data and store in local structures, now delete the reader.
		delete fmodel;
		
		//create shape model
		shapeModel* model1 = new shapeModel(Smean, Smodes, se, Imean, Imodes,Iprec, ie, M,Errs,polygons,labels);
		
		return model1;
		
	}
	void readBvars(const string & filename, const string & image_path, string & modelname, unsigned int & Nsubjects, vector<string> & v_imagenames, vector< vector<float> > & v_bvars, vector< vector<float> > & v_xfms)
	{
		v_bvars.clear();
		v_xfms.clear();
		v_imagenames.clear();
		//for now image name is not used
		string stemp;
		string imagename;
		unsigned int Nbvars;
		
		//lets start reading the mode parameters from the file
		
		ifstream fin;
		fin.open(filename.c_str());
		//cout<<"open file "<<filename<<endl;

		getline(fin,stemp);
		//cout<<"stemp "<<stemp<<endl;
		fin>>modelname;
		fin>>stemp>>Nsubjects;
		//cout<<"Nsubjects "<<Nsubjects<<endl;
		while (Nsubjects!=0)
		{
			vector<float> one_sub_bvars;
			vector<float> one_sub_xfm;
			fin>>imagename>>Nbvars;
			
			//cout<<"imagepath "<<image_path<<endl;
			if (image_path != "")
			{//add slash if it doesn't end with /
				////cout<<"last bit "<<*(image_path.rbegin())<<endl;
				if ( *(image_path.rbegin()) == '/' )
					imagename=image_path +imagename;
				else
					imagename=image_path +"/"+imagename;

			
			}	
			v_imagenames.push_back(imagename);
			char blank;
			fin.read(reinterpret_cast<char*>(&blank),sizeof(char));
			float bvar;
//			//cout<<"BVARS"<<endl;
			//cout<<imagename<<endl;
			for (unsigned int i=0;i<Nbvars;i++)
			{
				fin.read(reinterpret_cast<char*>(&bvar),sizeof(float));
				one_sub_bvars.push_back(bvar);
		//		//cout<<bvar<<" ";
			}
		//	//cout<<endl;
	//		//cout<<"TRANSFORMATION MATRIX"<<endl;
			
			for (int i=1;i<=16;i++)
			{
				float bvar;
				fin.read(reinterpret_cast<char*>(&bvar),sizeof(float));
				one_sub_xfm.push_back(bvar);
	//			//cout<<bvar<<" ";
	//			if (i%4 == 0)
	//				//cout<<endl;
			}
			
			v_bvars.push_back(one_sub_bvars);
			v_xfms.push_back(one_sub_xfm);
			
			
			Nsubjects--;
		}
	}
    ///NEED TO CHECK CORRESPONDANCE BETWEEN THIS AND PREVIOUS METHODS
    //this reconstructs a surface form bvars, will only use first surface
    void reconAllSurfacesAndSave( const string & filename, const string & out_base )
    {
        cout<<"recon and save"<<endl;
        ofstream fsurf_log((out_base+"_list.txt").c_str());

        string modelname;
		unsigned int Nsubjects;
		vector<string> v_imagenames;
		vector< vector<float> > v_bvars;
		vector< vector<float> > v_xfms;
		Matrix xfm_mni(4,4);
		xfm_mni=0;
		//load bvars
        
        
        //"" is in placve of image path
		readBvars(filename, "",  modelname, Nsubjects, v_imagenames, v_bvars, v_xfms);
        shapeModel* model1=loadAndCreateShapeModel(modelname);
        cout<<"shape model loaded"<<endl;
        //this is reconstructing the models mean surface
		//-------------save mean mesh--------------//
		{
            // volume<float> mni;
            
			volume<float>* immni = new volume<float>();
            
            char* fsldir = getenv("FSLDIR");
            read_volume_hdr_only(*immni, string(fsldir) + "/data/standard/MNI152_T1_1mm");
            
			//read_volume_hdr_only(*immni, template_name);
			
            //	vector< float > csys_eye(16,0);
            //	csys_eye[0]=1;
            //	csys_eye[5]=1;
            //	csys_eye[10]=1;
            //	csys_eye[15]=1;
			//	//cout<<"add coord system "<<i_xfms->at(0)<<endl;
            
			Matrix xfm_mni = immni->newimagevox2mm_mat();
			Matrix m_mm2vox_mni(4,4);
			m_mm2vox_mni=0;
			m_mm2vox_mni.element(0, 0) = 1.0/immni->xdim();
			m_mm2vox_mni.element(1, 1) = 1.0/immni->ydim();
			m_mm2vox_mni.element(2, 2) = 1.0/immni->zdim();
			m_mm2vox_mni.element(3, 3) = 1.0;
			
			xfm_mni = xfm_mni * m_mm2vox_mni;
			////cout<<"xfmmni2 "<<xfm_mni<<endl;
            
			vector<float> fmat;
			for (int i = 0 ; i < 4; ++i)
				for (int j = 0 ; j < 4; ++j)
					fmat.push_back(xfm_mni.element(i, j));
			
			delete immni;
            
            
            //-----out mean surface
            //save mean surface
            fslSurface<float, unsigned int> *surf = new fslSurface<float,unsigned int>(); 
            vector<float> bvars(1,0);
			surf->setVertices( model1->getDeformedGrid(bvars) );			
         	surf->setFaces( model1->cells );
            apply_xfm(*surf, fmat);
            //	//cout<<"add coord system "<<i_xfms->at(0)<<endl;
            
            //-----------------------------------------------
            vector< float > csys_eye(16,0);
			csys_eye[0]=1;
			csys_eye[5]=1;
			csys_eye[10]=1;
			csys_eye[15]=1;
            
        
			surf->addCoordSystem( csys_eye, "NIFTI_XFORM_SCANNER_ANAT");
            //--------------------------------------------------------
            
			Matrix flirtmat = xfm_mni.i();
			fmat.clear();
			for (int i = 0 ; i < 4; ++i)
				for (int j = 0 ; j < 4; ++j)
					fmat.push_back(flirtmat.element(i, j));
			
			surf->addCoordSystem( fmat, "FSL_XFORM_SCALEDMM");
            //--------------------------------------------------------
            
            
            writeGIFTI(*surf,out_base+"_modelmean.gii");
            delete surf;
            
            
            
		}
		//------------------------------------------//
        
        
        cout<<"done mean"<<endl;
        
        vector< vector<float> >::iterator i_xfms = v_xfms.begin();
		vector<string>::iterator i_imnames = v_imagenames.begin();
        unsigned int count = 0;
        //only for one set of bvars
       
		for ( vector< vector<float> >::iterator i_bvars = v_bvars.begin(); i_bvars != v_bvars.end(); ++i_bvars, ++i_xfms, ++i_imnames, ++count )
		{
            fslSurface<float, unsigned int> *surf = new fslSurface<float,unsigned int>(); 
            volume<float>* im = new volume<float>();
			read_volume_hdr_only(*im, *i_imnames );
			Matrix m_vox2mm = im->newimagevox2mm_mat();
			Matrix m_mm2vox(4,4);
			m_mm2vox=0;
			m_mm2vox.element(0, 0) = 1.0/im->xdim();
			m_mm2vox.element(1, 1) = 1.0/im->ydim();
			m_mm2vox.element(2, 2) = 1.0/im->zdim();
			m_mm2vox.element(3, 3) = 1.0;
			delete im;
			
			
			Matrix flirtmat(4,4);
			vector<float>::iterator i_one_xfm = i_xfms->begin();
			for (int i=0; i<4;i++){
				for (int j=0; j<4;j++, ++i_one_xfm){
					flirtmat.element(i,j) = *i_one_xfm;
					//cout<<flirtmat.element(i,j)<<" ";
				}
				//cout<<endl;
			}
			
            
			flirtmat=  m_vox2mm * m_mm2vox * flirtmat.i();
			//flirtmat= flirtmat.i();
            
			surf->setVertices( model1->getDeformedGrid(*i_bvars) );			
         	surf->setFaces( model1->cells );
            
            //----------------------------in native space----------------------------
            
			vector< float > csys_eye(16,0);
			csys_eye[0]=1;
			csys_eye[5]=1;
			csys_eye[10]=1;
			csys_eye[15]=1;
            
            // flirtmat = m_mm2vox * m_vox2mm;
            
            vector<float> fmat;
			for (int i = 0 ; i < 4; ++i)
			{
                
				for (int j = 0 ; j < 4; ++j)
				{
					//cout<<flirtmat.element(i, j)<<" ";
					fmat.push_back(flirtmat.element(i, j));
				}
				//cout<<endl;
			}
            apply_xfm(*surf, fmat);
            //	//cout<<"add coord system "<<i_xfms->at(0)<<endl;
			surf->addCoordSystem( csys_eye, "NIFTI_XFORM_SCANNER_ANAT");
            //--------------------------back to mni------------------------------
            //need to transform to scaled voxels, then applyxfm, the transform to NIFTI
            //cout<<"setmni "<<endl;
            //            flirtmat = flirtmat.i();
            
			flirtmat = xfm_mni *flirtmat.i() ;
            
       		fmat.clear();
			for (int i = 0 ; i < 4; ++i)
				for (int j = 0 ; j < 4; ++j)
					fmat.push_back(flirtmat.element(i, j));
			
			//cout<<"mniXFM "<<fmat.size()<<endl;
            //	apply_xfm<float, unsigned int,float>( surf, fmat );
			surf->addCoordSystem( fmat, "NIFTI_XFORM_MNI_152");
            //--------------------------------------------------------
            
			flirtmat = m_vox2mm * m_mm2vox ;
			flirtmat = flirtmat.i();
			fmat.clear();
			for (int i = 0 ; i < 4; ++i)
				for (int j = 0 ; j < 4; ++j)
					fmat.push_back(flirtmat.element(i, j));
			
			surf->addCoordSystem( fmat, "FSL_XFORM_SCALEDMM");
            //--------------------------------------------------------
            stringstream ss;
            ss<<count;
            string snum;
            ss>>snum;
            writeGIFTI(*surf,out_base+snum+".gii");
            fsurf_log<<out_base+snum+".gii"<<endl;
            delete surf;
        }
        
		fsurf_log.close();
		delete model1;
        
        
    }

    
	void reconSurface_from_bvars( fslSurface<float,unsigned int> & surf, const string & filename)//, const string & template_name)
	{
        
        string modelname;
		unsigned int Nsubjects;
		vector<string> v_imagenames;
		vector< vector<float> > v_bvars;
		vector< vector<float> > v_xfms;
		Matrix xfm_mni(4,4);
		xfm_mni=0;
		//load bvars
        
           
        
        //"" is in placve of image path
		readBvars(filename, "",  modelname, Nsubjects, v_imagenames, v_bvars, v_xfms);
        shapeModel* model1=loadAndCreateShapeModel(modelname);
        //this is reconstructing the models mean surface
		//-------------save mean mesh--------------//
		{
           // volume<float> mni;
          
			volume<float>* immni = new volume<float>();
            
            char* fsldir = getenv("FSLDIR");
            read_volume_hdr_only(*immni, string(fsldir) + "/data/standard/MNI152_T1_1mm");
            
			//read_volume_hdr_only(*immni, template_name);
			
		//	vector< float > csys_eye(16,0);
		//	csys_eye[0]=1;
		//	csys_eye[5]=1;
		//	csys_eye[10]=1;
		//	csys_eye[15]=1;
			//	//cout<<"add coord system "<<i_xfms->at(0)<<endl;

			xfm_mni = immni->newimagevox2mm_mat();
			Matrix m_mm2vox_mni(4,4);
			m_mm2vox_mni=0;
			m_mm2vox_mni.element(0, 0) = 1.0/immni->xdim();
			m_mm2vox_mni.element(1, 1) = 1.0/immni->ydim();
			m_mm2vox_mni.element(2, 2) = 1.0/immni->zdim();
			m_mm2vox_mni.element(3, 3) = 1.0;
			
			xfm_mni = xfm_mni * m_mm2vox_mni;
			////cout<<"xfmmni2 "<<xfm_mni<<endl;
            
			vector<float> fmat;
			for (int i = 0 ; i < 4; ++i)
				for (int j = 0 ; j < 4; ++j)
					fmat.push_back(xfm_mni.element(i, j));
			
			delete immni;
            
		}
		//------------------------------------------//

        
        
        
        vector< vector<float> >::iterator i_xfms = v_xfms.begin();
		vector<string>::iterator i_imnames = v_imagenames.begin();
        unsigned int count = 0;
        //only for one set of bvars
		for ( vector< vector<float> >::iterator i_bvars = v_bvars.begin(); i_bvars != v_bvars.begin()+1; ++i_bvars, ++i_xfms, ++i_imnames, ++count )
		{
            volume<float>* im = new volume<float>();
			read_volume_hdr_only(*im, *i_imnames );
			Matrix m_vox2mm = im->newimagevox2mm_mat();
			Matrix m_mm2vox(4,4);
			m_mm2vox=0;
			m_mm2vox.element(0, 0) = 1.0/im->xdim();
			m_mm2vox.element(1, 1) = 1.0/im->ydim();
			m_mm2vox.element(2, 2) = 1.0/im->zdim();
			m_mm2vox.element(3, 3) = 1.0;
			delete im;
			
			
			Matrix flirtmat(4,4);
			vector<float>::iterator i_one_xfm = i_xfms->begin();
			for (int i=0; i<4;i++){
				for (int j=0; j<4;j++, ++i_one_xfm){
					flirtmat.element(i,j) = *i_one_xfm;
					//cout<<flirtmat.element(i,j)<<" ";
				}
				//cout<<endl;
			}
			

			flirtmat=  m_vox2mm * m_mm2vox * flirtmat.i();
			//flirtmat= flirtmat.i();
                            
			surf.setVertices( model1->getDeformedGrid(*i_bvars) );			
         	surf.setFaces( model1->cells );
            
            //----------------------------in native space----------------------------
            
			vector< float > csys_eye(16,0);
			csys_eye[0]=1;
			csys_eye[5]=1;
			csys_eye[10]=1;
			csys_eye[15]=1;
            
           // flirtmat = m_mm2vox * m_vox2mm;
            
            vector<float> fmat;
			for (int i = 0 ; i < 4; ++i)
			{
                
				for (int j = 0 ; j < 4; ++j)
				{
					//cout<<flirtmat.element(i, j)<<" ";
					fmat.push_back(flirtmat.element(i, j));
				}
				//cout<<endl;
			}
            apply_xfm(surf, fmat);
            //	//cout<<"add coord system "<<i_xfms->at(0)<<endl;
			surf.addCoordSystem( csys_eye, "NIFTI_XFORM_SCANNER_ANAT");
            //--------------------------back to mni------------------------------
            //need to transform to scaled voxels, then applyxfm, the transform to NIFTI
            //cout<<"setmni "<<endl;
//            flirtmat = flirtmat.i();
            
			flirtmat = xfm_mni *flirtmat.i() ;
            
       		fmat.clear();
			for (int i = 0 ; i < 4; ++i)
				for (int j = 0 ; j < 4; ++j)
					fmat.push_back(flirtmat.element(i, j));
			
			//cout<<"mniXFM "<<fmat.size()<<endl;
            //	apply_xfm<float, unsigned int,float>( surf, fmat );
			surf.addCoordSystem( fmat, "NIFTI_XFORM_MNI_152");
            //--------------------------------------------------------
            
			flirtmat = m_vox2mm * m_mm2vox ;
			flirtmat = flirtmat.i();
			fmat.clear();
			for (int i = 0 ; i < 4; ++i)
				for (int j = 0 ; j < 4; ++j)
					fmat.push_back(flirtmat.element(i, j));
			
			surf.addCoordSystem( fmat, "FSL_XFORM_SCALEDMM");
            //--------------------------------------------------------
        
          
					}
				
		
		delete model1;


    }
    
    vector<float> getFSLtoNIFTIxfm( const volume<float> & imref)
    {
        Matrix m_vox2mm = imref.newimagevox2mm_mat();
        Matrix m_mm2vox(4,4);
        m_mm2vox=0;
        m_mm2vox.element(0, 0) = 1.0/imref.xdim();
        m_mm2vox.element(1, 1) = 1.0/imref.ydim();
        m_mm2vox.element(2, 2) = 1.0/imref.zdim();
        m_mm2vox.element(3, 3) = 1.0;
        
        Matrix M_fmat=m_vox2mm * m_mm2vox;
        vector<float> fmat;
        for (int i = 0 ; i < 4; ++i)
            for (int j = 0 ; j < 4; ++j)
                fmat.push_back(M_fmat.element(i, j));
        return fmat;
    }
    
	void reconSurfaces_from_bvars( vector< fslSurface<float,unsigned int> > & v_surf, const string & filename, const string & imagepath, const string & template_name, string out_base, bool reconMNI )
	{
//		for (vector<fslSurface*>::iterator i_surf = v_surf.begin(); v_surf.end();
		//will alter object or push back a new one
		ofstream f_surf_list( (out_base + "_surface_list.txt").c_str() );
		
		//-------------------mni tmeplate xfm info-------------//
		
		Matrix xfm_mni(4,4);
		xfm_mni=0;
		//-------------------mni tmeplate xfm info-------------//

		//	//cout<<"mni "<<immni.newimagevox2mm_mat()<<endl;	
		
		v_surf.clear();
		
		string modelname;
		unsigned int Nsubjects;
		vector<string> v_imagenames;
		vector< vector<float> > v_bvars;
		vector< vector<float> > v_xfms;
		
		//load bvars
		readBvars(filename, imagepath,  modelname, Nsubjects, v_imagenames, v_bvars, v_xfms);
		
		//get model name from file
		shapeModel* model1=loadAndCreateShapeModel(modelname);
		//this is reconstructing the models mean surface
		//-------------save mean mesh--------------//
		{
			
			volume<float>* immni = new volume<float>();
			read_volume_hdr_only(*immni, template_name);
			
			fslSurface<float,unsigned int>* surf_mean = new fslSurface<float,unsigned int>();
			surf_mean->setVertices( model1->smean );		
			surf_mean->setFaces( model1->cells );
			
			vector< float > csys_eye(16,0);
			csys_eye[0]=1;
			csys_eye[5]=1;
			csys_eye[10]=1;
			csys_eye[15]=1;
			//	//cout<<"add coord system "<<i_xfms->at(0)<<endl;
			surf_mean->addCoordSystem( csys_eye, "NIFTI_XFORM_MNI_152");
			xfm_mni = immni->newimagevox2mm_mat();
			Matrix m_mm2vox_mni(4,4);
			m_mm2vox_mni=0;
			m_mm2vox_mni.element(0, 0) = 1.0/immni->xdim();
			m_mm2vox_mni.element(1, 1) = 1.0/immni->ydim();
			m_mm2vox_mni.element(2, 2) = 1.0/immni->zdim();
			m_mm2vox_mni.element(3, 3) = 1.0;
			
			xfm_mni = xfm_mni * m_mm2vox_mni;
			////cout<<"xfmmni2 "<<xfm_mni<<endl;

			vector<float> fmat;
			for (int i = 0 ; i < 4; ++i)
				for (int j = 0 ; j < 4; ++j)
					fmat.push_back(xfm_mni.element(i, j));
			apply_xfm<float, unsigned int,float>( *surf_mean, fmat );
			writeGIFTI(*surf_mean, out_base + "_mean_mni.gii");
			delete immni;
			delete surf_mean; 

		}
		//------------------------------------------//
		
		////cout<<"xfmmni "<<xfm_mni<<endl;
		
		vector< vector<float> >::iterator i_xfms = v_xfms.begin();
		vector<string>::iterator i_imnames = v_imagenames.begin();
		unsigned int count = 0;
		for ( vector< vector<float> >::iterator i_bvars = v_bvars.begin(); i_bvars != v_bvars.end(); ++i_bvars, ++i_xfms, ++i_imnames, ++count )
		{
			
		//	////cout<<"imname "<<*i_imnames<<endl;
			volume<float>* im = new volume<float>();
			read_volume_hdr_only(*im, *i_imnames );
			Matrix m_vox2mm = im->newimagevox2mm_mat();
			Matrix m_mm2vox(4,4);
			m_mm2vox=0;
			m_mm2vox.element(0, 0) = 1.0/im->xdim();
			m_mm2vox.element(1, 1) = 1.0/im->ydim();
			m_mm2vox.element(2, 2) = 1.0/im->zdim();
			m_mm2vox.element(3, 3) = 1.0;
			delete im;
			
			
			Matrix flirtmat(4,4);
			vector<float>::iterator i_one_xfm = i_xfms->begin();
			for (int i=0; i<4;i++){
				for (int j=0; j<4;j++, ++i_one_xfm){
					flirtmat.element(i,j) = *i_one_xfm;
					//cout<<flirtmat.element(i,j)<<" ";
				}
				//cout<<endl;
			}
			
			//invert the matrix to go back to native space
			//flirtmat=  flirtmat.i() ;
		//	//cout<<"VOX2MM "<<m_vox2mm<<endl;
		//	//cout<<"MM2VOX "<<m_mm2vox<<endl;
			if (reconMNI){
				flirtmat=  m_vox2mm * m_mm2vox;
			}else {
				flirtmat=  m_vox2mm * m_mm2vox * flirtmat.i();
			}

			//* m_vox2mm * m_mm2vox;
			//reconned surfaces * xfmmed to native
			//flirtmat= m_vox2mm * m_mm2vox * flirtmat.i();
		//	Matrix flirtmat2 = m_vox2mm * m_mm2vox;
		//	cout<<"flirtmat 2 "<<flirtmat2<<endl;
			//flirtmat=  m_vox2mm * m_mm2vox * flirtmat.i();
	//		flirtmat=  flirtmat.i();

			
			vector<float> fmat;
			for (int i = 0 ; i < 4; ++i)
			for (int j = 0 ; j < 4; ++j)
					fmat.push_back(flirtmat.element(i, j));
			//vector<float> v_b(i_bvars->size(), 0);
			fslSurface<float,unsigned int>* surf = new fslSurface<float,unsigned int>();
		//	surf->setVertices( model1->getDeformedGrid(v_b) );//model1->getDeformedGrid(*i_bvars) );			
			surf->setVertices( model1->getDeformedGrid(*i_bvars) );			
	//		apply_xfm<float, unsigned int,float>( surf, fmat );
//			apply_xfm( surf, *i_xfms );
			//for (int j =0 ; j<16 ; ++j)
				//cout<<fmat[j]<<" ";
			//cout<<endl;
			//print<float>(fmat,"inverse fmat");
			surf->setFaces( model1->cells );

//----------------------------in native space----------------------------

			vector< float > csys_eye(16,0);
			csys_eye[0]=1;
			csys_eye[5]=1;
			csys_eye[10]=1;
			csys_eye[15]=1;
		//	//cout<<"add coord system "<<i_xfms->at(0)<<endl;
			surf->addCoordSystem( csys_eye, "NIFTI_XFORM_SCANNER_ANAT");

//--------------------------back to mni------------------------------
			//cout<<"setmni "<<endl;
			flirtmat = flirtmat.i();
		
			flirtmat = xfm_mni *flirtmat ;

			fmat.clear();
			for (int i = 0 ; i < 4; ++i)
			{
				 
				for (int j = 0 ; j < 4; ++j)
				{
					//cout<<flirtmat.element(i, j)<<" ";
					fmat.push_back(flirtmat.element(i, j));
				}
				//cout<<endl;
			}
			
			//cout<<"mniXFM "<<fmat.size()<<endl;
			apply_xfm<float, unsigned int,float>( *surf, fmat );
			surf->addCoordSystem( fmat, "NIFTI_XFORM_MNI_152");
//--------------------------------------------------------
			flirtmat = m_vox2mm * m_mm2vox ;
			flirtmat = flirtmat.i();
			fmat.clear();
			for (int i = 0 ; i < 4; ++i)
				for (int j = 0 ; j < 4; ++j)
					fmat.push_back(flirtmat.element(i, j));
			
			surf->addCoordSystem( fmat, "FSL_XFORM_SCALEDMM");
//--------------------------------------------------------

			
			string outname;
			
			if ( out_base == "NULL" )
			{
	//			//cout<<"imname end "<<i_imnames->substr(i_imnames->length()-3,3)<<endl;
				if ( i_imnames->substr(i_imnames->length()-3,3) == ".gz" )
				{
					i_imnames->erase(i_imnames->length()-3,3);
				}
				
				if  ( i_imnames->substr(i_imnames->length()-4,4) == ".nii" )
				{
					i_imnames->erase(i_imnames->length()-4,4);
				}
				if  ( i_imnames->substr(i_imnames->length()-4,4) == ".hdr" )
				{
					i_imnames->erase(i_imnames->length()-4,4);
				}
				if  ( i_imnames->substr(i_imnames->length()-4,4) == ".img" )
				{
					i_imnames->erase(i_imnames->length()-4,4);
				}

				
			}else {
				stringstream ss;
				ss<<count;
				string temp;
				ss>>temp;
				if (count< 10)
					temp="000"+temp;
				else if (count < 100)
					temp="00"+temp;
				else if (count < 1000)
					temp="0"+temp;
				
				outname=out_base+"_"+temp+".gii";
				
			}
			f_surf_list<<outname<<endl;
			writeGIFTI(*surf, outname);
			v_surf.push_back(*surf);
		}
		
		f_surf_list.close();
		
		
		delete model1;
	/*	vector<string> subjectnames;
		Matrix bvars;
		vector<int> vN;	
		vector< vector< vector<float> > > fmatv;
		read_bvars(inname.value(),bvars,subjectnames, vN, pathname.value(),fmatv);
		vector<float> vars;
		for (int i=0; i<bvars.Nrows();i++){
			vars.push_back(bvars.element(i,0));
		}
		
		
		for (int i=0; i<4;i++){
			for (int j=0; j<4;j++){
				flirtmat.element(i,j)=fmatv.at(0).at(i).at(j);
				//	//cout<<flirtmat.element(i,j)<<" ";
			}
			////cout<<endl;
		}
		//flirtmat=flirtmat.i();
		for (int i=0; i<4;i++){
			for (int j=0; j<4;j++){
				fmatv.at(0).at(i).at(j)=flirtmat.element(i,j);
				//	//cout<<flirtmat.element(i,j)<<" ";
			}
			////cout<<endl;
		}
		
		//	model1->registerModel(fmatv.at(0));
		//cout<<"vars "<<vars.size()<<endl;
		//can do this because FIRST saves bavrs/surfaces relative to original space.
		Mesh m=convertToMesh(model1->getDeformedGrid(vars),model1->cells);	
		meshReg(m,fmatv.at(0));
		
		m.save(outname.value(),3);
	 */
	}	
	
	vector<float> meshRegLeastSq( const fslSurface<float,unsigned int> & m_src, const  fslSurface<float,unsigned int> & m_target, const unsigned int & dof )
	{
		//demenad vertices
		//not most efficient but copied from meshutils
		Matrix Points(m_src.getNumberOfVertices(),3);
		Matrix refPoints(m_target.getNumberOfVertices(),3);
		//cout<<"point "<<m_src->getNumberOfVertices()<<" "<<Points.Nrows()<<endl;
		unsigned int count=0;
		for ( vector< vertex<float> >::const_iterator i_v = m_src.const_vbegin(); i_v != m_src.const_vend(); ++i_v, ++count )
		{
			Points.element(count,0) = i_v->x;
			Points.element(count,1) = i_v->y;
			Points.element(count,2) = i_v->z;

		}
		
		count=0;
		for ( vector< vertex<float> >::const_iterator i_v = m_target.const_vbegin(); i_v != m_target.const_vend(); ++i_v, ++count )
		{
			refPoints.element(count,0) = i_v->x;
			refPoints.element(count,1) = i_v->y;
			refPoints.element(count,2) = i_v->z;
			
		}
		
			//cout<<Points<<endl;
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
		//cout<<"Points2"<<endl;

		//demean data 
		//both meshes must have equal number odf vertices
		for (int i =0; i<Points.Nrows();i++){
			//cout<<"i "<<i<<" "<<Points.Nrows()<<endl;
			PointsDM.element(i,0)=Points.element(i,0)-mean1x;
			PointsDM.element(i,1)=Points.element(i,1)-mean1y;
			PointsDM.element(i,2)=Points.element(i,2)-mean1z;
			
			refPointsDM.element(i,0)= refPoints.element(i,0) -meanRefx;
			refPointsDM.element(i,1)= refPoints.element(i,1) -meanRefy;
			refPointsDM.element(i,2)= refPoints.element(i,2) -meanRefz;
			
		}
		//cout<<"Points3"<<endl;

		//calculate translation component
		float tx=meanRefx - mean1x ;
		float ty=meanRefy - mean1y ;
		float tz=meanRefz - mean1z ;
		
		//calculate scale
		float ssq1[4]={0,0,0,0};
		float ssqRef[4]={0,0,0,0};
		float scale[4]={1,1,1,1};
		//cout<<"Points4"<<endl;

		if (dof==7)
		{
			

			for (int i =0; i<Points.Nrows();i++)
				ssq1[0]+=PointsDM.element(i,0)*PointsDM.element(i,0) + PointsDM.element(i,1)*PointsDM.element(i,1) + PointsDM.element(i,2)*PointsDM.element(i,2);
			for (int i =0; i<refPoints.Nrows();i++)
				ssqRef[0]+=refPointsDM.element(i,0)*refPointsDM.element(i,0) + refPointsDM.element(i,1)*refPointsDM.element(i,1) + refPointsDM.element(i,2)*refPointsDM.element(i,2) ;
			scale[0]=sqrt(ssqRef[0]/ssq1[0]);
		//	cout<<"scale "<<ssqRef[0] << " "<<ssq1[0]<< " "<<scale[0]<<endl;
		//	cout<<"points1 "<<PointsDM<<endl;
		//	cout<<"pointsREf "<<refPointsDM<<endl;
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
		//cout<<"scale calculated "<<scale[0]<<" "<<scale[1]<<" "<<scale[2]<<" "<<scale[3]<<endl;
		//calculate rotation
		Matrix M=refPointsDM.t()*PointsDM;
		//cout<<"translation "<<mean1x<<" "<<mean1y<<" "<<mean1z<<" "<<tx<<" "<<ty<<" "<<tz<<endl;
		//cout<<"translation "<<meanRefx<<" "<<meanRefy<<" "<<meanRefz<<" "<<tx<<" "<<ty<<" "<<tz<<endl;
		
		Matrix U;
		DiagonalMatrix D;
		SVD(M.t()*M,D,U);
		//M should always be a 3x3 matrix
		for (int i=0;i<D.Nrows();i++){
			D.element(i)=1/sqrt(D.element(i));
		}
//		//cout<<"Rotation Matrix ..."<<endl;
		Matrix R(3,3);
		//cout<<D.Nrows()<<" "<<U.Nrows()<<" "<<U.Ncols()<<" "<<M.Nrows()<<" "<<refPointsDM.Nrows()<<" " <<refPointsDM.Ncols()<<endl; 
		R=M*(U*D*U.t());
		//R=R.t();
//		//cout<<"Rotation Matrix "<<R<<endl;
		
		
		Matrix fmat=getIdentityMatrix(4);
//		//cout<<"fmat1 "<<endl;
//		//cout<<fmat<<endl;
		preMultiplyTranslation(fmat,-mean1x,-mean1y,-mean1z);
//		//cout<<"fmat2"<<endl;
//		//cout<<fmat<<endl;
		if (dof==7)
			preMultiplyGlobalScale(fmat,scale[0]);
		else if (dof==9)
			preMultiplyGlobalScale(fmat,scale[1],scale[2],scale[3]);
		
//		//cout<<"fmat3"<<endl;
//		//cout<<fmat<<endl;
		if (dof>3)
			preMultiplyGlobalRotation(fmat,R);	
//		//cout<<"fmat4 "<<mean1x<<" "<<tx<<endl;
//		//cout<<fmat<<endl;
		preMultiplyTranslation(fmat,mean1x+tx,mean1y+ty,mean1z+tz);
//		//cout<<"fmat4"<<endl;
//		//cout<<fmat<<endl;
		
		vector<float> vfmat(16,0);
		for (int i = 0 ; i < 4 ; ++i)
			for (int j = 0; j < 4; ++j) {
				vfmat[4*i +j ] = fmat.element(i,j);
			}
		return vfmat;
	}
	
	//----------------------------This runs the vertex analysis--------------------//
	void vertexMVglm( fslSurface<float,unsigned int> & surf_stats, const string & inname,  const string & designname, const string & saveVertices , string normname ){
		//**********read in bvars and models and lfirt matrices***************//
		Matrix MeshVerts;//this is used when placing t-stats on a mesh
		Matrix target;	//target is only used for glm

		target=read_vest(designname);
		bool useNorm=false;
		if ( normname != "none" )
			useNorm=true;
		
		
		ifstream f_mesh(inname.c_str());
		string meshname;
		unsigned int Nverts = 0;
		unsigned int Nsubjects = 0; 
		//read in surfaces concatenate into one large martrix
		while ( f_mesh >> meshname ) {
			//cout<<"read surface "<<meshname<<endl;
			fslSurface<float, unsigned int> *mesh = new fslSurface<float, unsigned int>();

			read_surface(*mesh,meshname);
			if (Nsubjects==0)
				Nverts = mesh->getNumberOfVertices();
			unsigned int vert_index=0;
			ColumnVector one_mesh_verts(3*Nverts);

			for ( vector< vertex<float> >::iterator i_m= mesh->vbegin(); i_m != mesh->vend(); ++i_m,vert_index+=3) {
			//	//cout<<vert_index<<" "<<3*Nverts<<endl;
				one_mesh_verts.element(vert_index) = i_m->x;
				one_mesh_verts.element(vert_index+1) = i_m->y;
				one_mesh_verts.element(vert_index+2) = i_m->z;

			}

			if ( Nsubjects == 0 )
				MeshVerts = one_mesh_verts;
			else {
				MeshVerts = MeshVerts | one_mesh_verts;
			}

			
			++Nsubjects;
			delete mesh;
		}
		f_mesh.close();

      
        //----------------------Do Design matrix stuff --------------------//
		
		bool isOne=true;
		//checks for mean Column as first column, if it decides there is a column of ones (first column) then 
		//it will not demean the design matrix
		for (int i=0;i<target.Nrows();i++){
			if (target.element(i,0)!=1){
				isOne=false;
			}
		}
		
		int EVmin=1;//ignore mean column first EV to examine
		if(!isOne){
	
			//create demena deisgn, assume a nomean column at start
			Matrix targTemp(target.Nrows(),target.Ncols()+1);
			for (int i=0;i<targTemp.Ncols();i++){
				if (i==0){
					for (int j=0;j<target.Nrows();j++){
						targTemp.element(j,i)=1;
					}
				}else{
					float mean=0;
					for (int j=0;j<target.Nrows();j++){
						mean+=target.element(j,i-1);
					}	
					mean/=target.Nrows();
					for (int j=0;j<target.Nrows();j++){
						targTemp.element(j,i)=target.element(j,i-1)-mean;
					}
				}
			}
			
			target=targTemp;
			EVmin=1;
		}
		//this displays the demean design matrix
        cout<<"Design  Matrix "<<endl;
		for (int j=0;j<target.Nrows();j++){
			for (int i=0;i<target.Ncols();i++){	
				cout<<target.element(j,i)<<" ";
			}	
			cout<<endl;
		}
		
		//-------------------Do Vertex Analysis Stuff --------------------------//
		{
    
			ColumnVector CVnorm(target.Nrows());
			// if chosen can include a normalization EV (i.e. control for size) ...useful for vertex shape statistics
			if (useNorm){
				ifstream normIn;
				float normtemp;
				normIn.open(normname.c_str());
				for (int subject=0;subject<target.Nrows();subject++){
					normIn>>normtemp;
					CVnorm.element(subject)=normtemp;
				}
				target=target | CVnorm ;
			}
			
			
			
			//****************END DEMEAN DESIGN MATRIX AND ADD ONE COLUMS*********************/
			int Nrows = target.Ncols()-1;
		//	if (isOne)
		//		++Nrows;
	
			Matrix contrast(Nrows,target.Ncols());
			
			//int EVmax=1;//EV to include
			int	EVmax=target.Ncols();
	
			for (int EV=EVmin;EV<EVmax;EV++){
                cout<<"EV "<<EV<<endl;
				
				//append the EV number to the output
				stringstream ev2st;
				ev2st<<EV;
				string evnum;
				ev2st>>evnum;
					
				vector<float> mean_verts;
				for ( unsigned int v_ind = 0 ; v_ind < Nverts*3; v_ind += 3 )
				{ 
                   // cout<<"vind "<<v_ind<<endl;
					float meanx1=0,meany1=0,meanz1=0,meanx2=0,meany2=0, meanz2=0;
					// for (int i=0;i<MeshVerts.Nrows();i=i+3){
					int n1=0,n2=0;
					for (int j=0;j<MeshVerts.Ncols();j++){
                     //   cout<<"j "<<j<<endl;
						if (target.element(j,EVmin) <= 0 ){ //handles demeaned data
							meanx1+=MeshVerts.element(v_ind,j);
							meany1+=MeshVerts.element(v_ind+1,j);
							meanz1+=MeshVerts.element(v_ind+2,j);
							n1++;
						}else{
							meanx2+=MeshVerts.element(v_ind,j);
							meany2+=MeshVerts.element(v_ind+1,j);
							meanz2+=MeshVerts.element(v_ind+2,j);
							n2++;
						}
					}
					// }
					meanx1/=n1;
					meany1/=n1;
					meanz1/=n1;
					
					meanx2/=n2;
					meany2/=n2;
					meanz2/=n2;
					
					
					mean_verts.push_back(meanx1);
					mean_verts.push_back(meany1);
					mean_verts.push_back(meanz1);
                }
				
            cout<<"set hsape Mesh"<<endl;
				vector<float> scalarsT;
				
				//this provides the option to do stats on the vertices 
				//create F test contrast matrix (testing each EV separately
				
				unsigned int count=0;
				for (int j=0;j<contrast.Ncols();j++){
					if (j!=EV){
						for (int i=0;i<contrast.Ncols();i++){
							if(i==j){
								//if(i==EV){
								contrast.element(count,i)=1;
							}else{
								contrast.element(count,i)=0;
							}
							
							////cout<<contrast.element(count,i)<<" ";
						}
						////cout<<endl;
						count++;
					}
				}
				cout<<"contrats"<<endl<<contrast<<endl;
				// save vertex Matrix if requested
				if (saveVertices != "") {
					Matrix ReshapedVerts(MeshVerts.Nrows(),MeshVerts.Ncols());
					int nverts=MeshVerts.Nrows()/3;
					int newm=1;
					for (int m=1; m<=ReshapedVerts.Nrows(); m+=3) {
						for (int n=1; n<=ReshapedVerts.Ncols(); n++) {
							ReshapedVerts(newm,n)=MeshVerts(m,n);
							ReshapedVerts(newm+nverts,n)=MeshVerts(m+1,n);
							ReshapedVerts(newm+2*nverts,n)=MeshVerts(m+2,n);
						}
						newm++;
					}
					write_ascii_matrix(fslbasename(saveVertices)+"_mat.txt",ReshapedVerts);
					volume4D<float> vertices(ReshapedVerts.Nrows()/3,3,1,ReshapedVerts.Ncols());
					vertices.setmatrix(ReshapedVerts.t());
					save_volume4D(vertices,saveVertices);
				}
                cout<<"done sabe vertices"<<endl;
				//use multivariate test on each vertex
				int dof1=0, dof2=0;
				for (int i=0;i<MeshVerts.Nrows();i=i+3){
					//use multivariate multiple regression on each veretx 
					float F=MVGLM_fit(target, MeshVerts.SubMatrix(i+1,i+3,1,MeshVerts.Ncols()).t(),contrast,dof1,dof2);
						//	tstatsx.at(i/3)*=wilkL;
					//	tstatsy.at(i/3)*=wilkL;
					//	tstatsz.at(i/3)*=wilkL;	
					//this may 
					scalarsT.push_back(F);
				}
                cout<<"done mvglm"<<endl;
				read_surface(surf_stats,meshname);
				surf_stats.setVertices(mean_verts);//model1->smean);
                vector<int> vdofs(2);
                vdofs[0]=dof1;
                vdofs[1]=dof2;
                cout<<"dof12 "<<dof1<<" "<<dof2<<endl;
                surf_stats.insertNonVertScalars(vdofs,0, "FStatDOFs");
                stringstream ss;
                ss<<EV;
                string snum;
                ss>>snum;
                surf_stats.insertScalars(scalarsT, 0, "F-stats_EV_"+snum);
  

			}
		}
		
	}		
	
	//-----------------------------Transformation Matrix Utilities Begin---------------------------//
	void preMultiplyGlobalRotation(Matrix & fmat, const Matrix & R)
	{
		//assume a 3by 3 rotation matrix
		//SubMatrix command has indices starting from 1
		fmat.SubMatrix(1,3,1,3)=R*fmat.SubMatrix(1,3,1,3);
		fmat.SubMatrix(1,3,4,4)=R*fmat.SubMatrix(1,3,4,4);
	}
	void  preMultiplyGlobalScale(Matrix & fmat, const float & s)
	{
		for (int col=0;col<4;col++)
			for (int row=0;row<3;row++)
				fmat.element(row,col)*=s;
	}
	void  preMultiplyGlobalScale(Matrix & fmat, const float & sx,const float & sy, const float & sz)
	{
		for (int col=0;col<4;col++)
			fmat.element(0,col)*=sx;
		for (int col=0;col<4;col++)
			fmat.element(1,col)*=sy;
		for (int col=0;col<4;col++)
			fmat.element(2,col)*=sz;
	}
	void preMultiplyTranslation(Matrix & fmat, const  float & tx, const float & ty, const float & tz )
	{
		//must be a 4 by 4 matrix 
		fmat.element(0,3)+=tx;
		fmat.element(1,3)+=ty;
		fmat.element(2,3)+=tz;
	}
	ReturnMatrix getIdentityMatrix(const short N)
	{
		Matrix fmat(N,N);
		fmat=0;
		for (int i=0;i<N;i++)
			fmat.element(i,i)=1;
		fmat.Release();
		return fmat;
	}
	
	//-----------------------------Tranformation Matrix Utilities End---------------------------//
	
    vector<double> invertMatrix( const vector<double> & fmat)
    {
        vector<double> fmati(16);
        Matrix Mfmat(4,4);
   // print<doubl)
        vector<double>::const_iterator i_f=fmat.begin();
		for (int i=0;i<4;i++)
        {
            for (int j=0;j<4;j++,++i_f)
            {
                cout<<*i_f<<" ";
                Mfmat.element(i,j)=*i_f;
            }
            cout<<endl;
        }
        Mfmat=Mfmat.i();
        cout<<"inverted"<<endl;
         vector<double>::iterator i_fi=fmati.begin();
		for (int i=0;i<4;i++)
            for (int j=0;j<4;j++,++i_fi)
               *i_fi =  Mfmat.element(i,j);
		
        Mfmat.Release();
		return fmati;
    }


	
	
}
