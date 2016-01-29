/*
 *  meshUtils.h
 *  
 *
 *  Created by Brian Patenaude on 04/04/2008.
 *  Copyright 2008 University of Oxford All rights reserved.
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
#ifndef MESHUTILS_H
#define MESHUTILS_H

#include "newimage/newimageall.h"
#include "fslvtkio/fslvtkio.h"
//#include "shapeModel/shapeModel.h"

//#include "meshclass/meshclass.h"
//using namespace std;
//using namespace fslvtkio;
namespace meshutils {
	
	class meshUtils: public fslvtkio::fslvtkIO {
		
public: 
		meshUtils();
		meshUtils(const string & fname, const fslvtkIO::DataType i);

		//convience method to read polydata
		void loadMesh(const std::string & meshname);

		float getMinX() { return Points.Column(1).Minimum(); }
		float getMaxX() { return Points.Column(1).Maximum(); }
		float getMinY() { return Points.Column(2).Minimum(); }
		float getMaxY() { return Points.Column(2).Maximum(); }	
		float getMinZ() { return Points.Column(3).Minimum(); }
		float getMaxZ() { return Points.Column(3).Maximum(); }
		
		
		static void generateRandomMeshUsingScalar(const mesh::Mesh & m, const string & outname, const vector<bool> & scal, const int & N);
	//	static void addModesToModelUsingMask(shapemodel::shapeModel * min, const vector<bool> & scal);
	//	static void getConditionalMeanAndVariance(shapemodel::shapeModel * min, volume4D<float> & iCondMean, volume4D<float> & iCondVar , const volume<float> & im, const int & mode, const float & bmin, const float & bmax, const float & res, const float & mean_offset);
		static void generateRandom6DOFMatrices( const string & outname, const int & N);
	
		int getNumberOfPolygons(){ return Polygons.Nrows(); }
		
		float interpolateScalar(const unsigned int & tri_index, const float & x, const float & y, const float& z );

		static ReturnMatrix vectorOfVectorsToMatrix(const vector< vector<float> > & vin);
		static void fileToVector(const string & fname, vector<string> list);
		static vector<string> fileToVector(const string & fname );

		static ReturnMatrix readFlirtMat(const string & fname);
		static void writeFlirtMatrix(const Matrix & fmat, const string & fname);

			
		void getBounds(int *bounds, const float & xdim, const float & ydim,const float & zdim);
		static void getBounds(const mesh::Mesh & m, int *bounds,const float & xdim, const float & ydim,const float & zdim) ;
	
	
		template<class Tdist,class Tim> //instantiations are in .cc file
		void SurfDistToLabels(vector<Tdist> & dist, const volume<Tim> & image);
		template<class Tdist,class Tim> //instantiations are in .cc file
		void SurfDistToLabels(vector<Tdist> & dist, const volume<Tim> & image, const Tim & label);
		
		void SurfScalarsMeanAndStdev(vector<string> meshList, Matrix & MeanPoints, Matrix & MeanScalars, Matrix & StDevScalars );
		
		static void meshReg(mesh::Mesh & m, const Matrix & fmat);
		void meshReg(const Matrix & fmat);
		
		static void shift3DVertexMatrix(Matrix & mat, const float & tx, const float & ty, const float & tz );
		static void shift3DVertexColumnVector(ColumnVector & mat, const float & tx, const float & ty, const float & tz );
		static void shift3DMesh(mesh::Mesh & m, const float & tx, const float & ty, const float & tz );
		
		static ReturnMatrix subSampleMatrix(const Matrix & m, const vector<bool> & vmask );
		static ReturnMatrix subSample_Nby1_3D_Matrix(const Matrix & m, const vector<bool> & vmask );


		void shiftPoints(const float & tx, const float & ty, const float & tz );
		void scalePoints( const  float & sx, const float & sy, const float & sz );

		static ReturnMatrix shiftPolygonMatrix(const Matrix & mat, const int & shift );
		void shiftPolygonMatrix( const int & shift );
		static ReturnMatrix meshPointsToMatrix(const mesh::Mesh & m1);
		static bool checkLine(const float & p1, const float & p2, const float & test);
		static bool checkTriangleNeighbour(const short & tri0, const short & tri1, const short & tri2 , const short & ind0, const short & ind1, short & ind0new , short & ind1new);
		static void intersectionPoint(const float & ycut, const float & px0, const float & py0, const float & pz0, const  float & dx, const float & dy, const float & dz, vector<float> & px, vector<float> & py, vector<float> & pz);

		template<class T,class T2>
		void deformSurface(const volume<T> & im, const float & maxit, const float & w_im, const float & wTang, const float & maxTri, const float & w_norm, const T & max_thresh,const unsigned int & interRate, const bool & enableInteraction, const string & name);


		//transformation matrix utilities
		static void preMultiplyGlobalRotation(Matrix & fmat, const Matrix & R);
		static void preMultiplyGlobalScale(Matrix & fmat, const float & s);
		static void preMultiplyGlobalScale(Matrix & fmat, const float & sx,const float & sy, const float & sz);
		static void preMultiplyTranslation(Matrix & fmat, const  float & tx, const float & ty, const float & tz );
		static ReturnMatrix getIdentityMatrix(const short N);
		//end of transofrmation matrix utilities
		
		//this should output the mask values as well as the truncated mesh
		template<class T>
		void sampleImageAtPoints(const volume<T> & immask, vector<T> & vsamples);

		void LQSurfaceReg(const Matrix & refPoints, Matrix & fmat, const int & dof);

		void combineMeshesWithVectorsAndScalars(const vector<string> & meshlist);

		//template<class T>
		void findMidPointOfMidSlice(const volume<char> & im, const Matrix & fmat, float & cx, float & cy, float & cz);
		vector<float> sliceMesh(const float & ycut);
		
		void sampleMeshProfilesFromImage(const volume<float> & image, const float & sample_interval, const unsigned int & ipp);

		
		static void warpMeshWithDefField(const string & fieldname, const string & meshname, const string & meshoutname, const float & dx, const float & dy, const float & dz);
	
		template< class T >
		void warpGridWithDefField(const volume4D<T> & fieldname, const float & dx, const float & dy, const float & dz);

		template< class T >
		static void warpGridWithDefField(const volume4D<T> & defField, vector<float> & points_in, float warpSc,const float & dx, const float & dy, const float & dz);
		

		float drawTriangleScalars(volume<float>& image, volume<int> &count, const unsigned int & tri_index);

		
		
//-----------------------VERTEX ANALYSIS STUFF-----------------------//
//return linear transformation matrix 
Matrix reg_leastsq(const Matrix & TargetPoints,  const short & dof);


static void applyReg(Matrix & pts, const Matrix & fmat);
static ReturnMatrix calculateRotation(const Matrix & Pts_src_dm, const Matrix & Pts_targ_dm);
static ReturnMatrix calculateScale(const Matrix & Pts_src_dm, const Matrix & Pts_targ_dm, const bool & global);

Matrix alignSurfaces(const string & src_list, const short & dof, const string & outname ); 
double maxScalar();
double meanScalar();


//-----------------------VERTEX ANALYSIS STUFF-----------------------//


		
		static float myatan2(const float & y, const float & x);
		
		template<class T>
		static void cartesianToSphericalCoord(vector<T> & verts);
		
		template<class T>
		static void sphericalToCartesianCoord(vector<T> & verts);
			
		static ReturnMatrix addSphericalCorrdinates( const Matrix & m1, const Matrix  & m2 );
		static ReturnMatrix subtractSphericalCoordinates( const Matrix & m1, const Matrix  &  m2 );
static ReturnMatrix  averageSphericalCorrdinates( const Matrix & m1, const Matrix & m2 , int & N1, const int & N2);
 static void SVDarcSpherical( Matrix & m1, DiagonalMatrix & D, Matrix & U, Matrix & V);

		
		
		void getSphericalCoordFromCart(NEWMAT::Matrix & r, NEWMAT::Matrix & theta, NEWMAT::Matrix & phi);
		static void cartesianToSphericalCoord(NEWMAT::Matrix & verts);
		static void sphericalToCartesianCoord(NEWMAT::Matrix & verts);
		
		template<class T>
		static void cartesianToSphericalCoord(mesh::Mesh & m);

		template<class T>
		static void sphericalToCartesianCoord(mesh::Mesh & m);
		
		static void combinedSharedBoundaries(const string & mname1, const string & mname2 );
		static void labelAndCombineSharedBoundaries(const string & mname1, const string & mname2, const string & mout1name );
		//static void appendSharedBoundaryMask(const string & mname1, const string & mname2,const string & mbase, const string & mout1name, const bool & indexed, const bool & useSc2 );
		ReturnMatrix appendSharedBoundaryMask(const Matrix & Points2 );
		
		//this will sample the new set of points based on indices from mask then replace the original
		void sampleSharedBoundaryByMask(const Matrix & Points2);

		static void removeRow(Matrix & mat, int ind );
		static ColumnVector sample_image_at_vertices(string meshname, string imagename);
		static void sample_image_at_verticesNN(const string & meshname, const string & imagename, const string & outname);
		static void sampleSumAndWrite(const string & meshname, const string & imagename, const string & outname);
		static void meshReg(const string & meshname, const string & fname, const string & outname, bool noinv );
		static void draw_segment(volume<short>& image, const mesh::Pt& p1, const mesh::Pt& p2, int label);
		static volume<short> draw_mesh(const volume<short>& image, const mesh::Mesh &m, int label);
		static volume<short> make_mask_from_mesh(const volume<float> & image, const mesh::Mesh& m, int label, int* bounds, const bool & sep_boundary);
		
		static void fillMesh(const string & imname, const string & meshname, const string & outname, const int & label, const bool & sep_boundary  );
	
		static bool findVertex(const Matrix & vert1,const Matrix & vert2, int ind1 );
	
		static void applyFlirtThenSBmask(const string & mname1, const string & mname2,const string & mflirtname, const string & mout1name);
		static void do_work_uncentreMesh(const string & inim, const string & inmesh, const string & outmesh);
		static void do_work_MeshReg(const string & inim, const string & inmesh, const string & outmesh);
		static void subtractMeshes(const string & mesh1name, const string & mesh2name, const string & out);
		
		template<class T>
		static void warpMeshWithDefField(mesh::Mesh & m, const volume4D<T> & defField, const Matrix & mat);
		
		static ReturnMatrix getDeformedVector(const ColumnVector & mean, const Matrix & modes, const ColumnVector & eigs, const vector<float> & vars  );


		template<class T>
		 void ugridToImage(NEWIMAGE::volume<T> & im);
		 
		 
private:
			volume<short> image;
			

	};
}
#endif
