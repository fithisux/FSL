/*  first_utils.cc
    Brian Patenaude, Mark Jenkinson and Matthew Webster
    Copyright (C) 2006-2010 University of Oxford  */
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

#include <stdio.h>
#include <unistd.h>


#include "utils/options.h"
#include "newimage/newimageall.h"
#include "meshclass/meshclass.h"
#include "shapeModel/shapeModel.h"
#include "miscmaths/miscmaths.h"
#include "fslvtkio/fslvtkio.h"

using namespace std;
using namespace NEWIMAGE;
using namespace Utilities;
using namespace mesh;
using namespace SHAPE_MODEL_NAME;
using namespace MISCMATHS;
using namespace fslvtkio;

string title="firt_utils University of Oxford (Brian Patenaude)";
string examples="first_utils [options] -i input -o output ";


Option<bool> verbose(string("-v,--verbose"), false, 
		     string("output F-stats to standard out"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<bool> debug(string("--debug"), false,
		  string("\tturn on debugging mode"),
		  false, no_argument);
Option<bool> overlap(string("--overlap"), false,
		     string("Calculates Dice overlap."),
		     false, no_argument);

Option<bool> useScale(string("--useScale"), false,
		      string("do stats"),
		      false, no_argument);
Option<bool> vertexAnalysis(string("--vertexAnalysis"), false,
			    string("Perform vertex-wise stats from bvars."),
			    false, no_argument);
Option<bool> singleBoundaryCorr(string("--singleBoundaryCorr"), false,
				string("Correct boundary voxels of a single structure."),
				false, no_argument);

Option<bool> usePCAfilter(string("--usePCAfilter"), false,
			  string("Smooths the surface my truncating the mode parameters."),
			  false, no_argument);
Option<bool> usebvars(string("--usebvars"), false,
		      string("Operate using the mode parameters output from FIRST."),
		      false, no_argument);
Option<bool> doMVGLM(string("--doMVGLM"), false,
		      string("doMVGLM."),
		      false, no_argument);
Option<bool> useReconMNI(string("--useReconMNI"), false,
			 string("Reconstruct meshes in MNI space."),
			 false, no_argument);

Option<bool> useReconNative(string("--useReconNative"), false,
			    string("Reconstruct meshes in native space."),
			    false, no_argument);
Option<bool> useRigidAlign(string("--useRigidAlign"), false,
			   string("Register meshes using 6 degree of freedom (7 if useScale is used)."),
			   false, no_argument);

Option<bool> useNorm(string("--useNorm"), false,
		     string("Normalize volumes measurements."),
		     false, no_argument);
Option<bool> surfaceVAout(string("--surfaceout"), false,
		     string("Output vertex analysis on the surface (old method)"),
		     false, no_argument);


Option<bool> reconMeshFromBvars(string("--reconMeshFromBvars"), false,
		       string("Convert bvars to mesh."),
		       false, no_argument);
Option<bool> readBvars(string("--readBvars"), false,
								string("Read bvars from binary format"),
								false, no_argument);
Option<bool> concatBvars(string("--concatBvars"), false,
					   string("Concat bvars from binary format"),
					   false, no_argument);
Option<bool> meshToVol(string("--meshToVol"), false,
					   string("Convert mesh to an image."),
					   false, no_argument);
Option<bool> centreOrigin(string("--centreOrigin"), false,
			  string("Places origin of mesh at the centre of the image"),
			  false, no_argument);
Option<string> saveVertices(string("--saveVertices"), string(""),
			  string("filename for saving matrix of vertex coords: (all x, then all y, then all z) by Nsubjects"),
			  false, requires_argument);
Option<string> inname(string("-i,--in"), string(""),
		      string("filename of input image/mesh/bvars"),
		      true, requires_argument);
Option<string> pathname(string("-a,--in"), string(""),
			string("Specifies extra path to image in .bvars file"),
			false, requires_argument);
Option<string> flirtmatsname(string("-f,--in"), string(""),
			     string("Text file containing filenames of flirt matrices (filenames, not numbers)."),
			     false, requires_argument);
Option<string> meshname(string("-m,--in"), string(""),
			string("Filename of input mesh"),
			false, requires_argument);
Option<string> normname(string("-g,--in"), string(""),
			string("Filename of normalization factors."),
			false, requires_argument);
Option<string> designname(string("-d,--in"), string(""),
			  string("Filename of fsl design matrix"),
			  false, requires_argument);

Option<string> refname(string("-r,--in"), string(""),
		       string("Filename of reference image "),
		       false, requires_argument);

Option<int> meshLabel(string("-l,--meshlabel"), 1,
		      string("Specifies the label used to fill the mesh."),
		      false, requires_argument);

Option<float> thresh(string("-p,--thresh"), 4,
		     string("Threshhold for clean up."),
		     false, requires_argument);
Option<int> numModes(string("-n,--numModes"), 0,
		     string("Number of modes to retain per structure."),
		     false, requires_argument);
Option<string> outname(string("-o,--out"), string(""),
		       string("Output name"),
		       true, requires_argument);

int nonoptarg;

////////////////////////////////////////////////////////////////////////////
//global variables

int refXsize=182;
int refYsize=218;
int refZsize=182;
float refXdim=1.0;
float refYdim=1.0;
float refZdim=1.0;


void setShapeMesh(shapeModel * model1, const Mesh & m){	
  vector<float>::iterator j=model1->smean.begin();
  for (vector<Mpoint*>::const_iterator i = m._points.begin(); i!=m._points.end(); i++ , j+=3){
    *j=(*i)->get_coord().X;
    *(j+1)=(*i)->get_coord().Y;
    *(j+2)=(*i)->get_coord().Z;
  }	
}

void meshReg(Mesh & m, const string & flirtmatname){

  //refsize is actually target image
  int numPoints=m.nvertices();
  Matrix MeshPts(4,numPoints);
  Matrix NewMeshPts(4,numPoints);
  Matrix flirtmat(4,4);
	
  //read in flirt matrix, uses ascii
	
  ifstream fmat;
  fmat.open(flirtmatname.c_str());
  float tmpfloat=0;
  for (int i=0; i<4;i++){
    for (int j=0; j<4;j++){
      fmat>>tmpfloat;
      flirtmat.element(i,j)=tmpfloat;
      //	cout<<flirtmat.element(i,j)<<" ";
    }
    //cout<<endl;
  }
  flirtmat=flirtmat.i();
	
	
  //	cout<<"transform mesh points..."<<endl;
  int count=0;
  for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
    //	cout<<"count "<<count<<endl;	
		
    MeshPts.element(0,count)=(*i)->get_coord().X;
    MeshPts.element(1,count)=(*i)->get_coord().Y;
    MeshPts.element(2,count)=(*i)->get_coord().Z;
    MeshPts.element(3,count)=1;
		
    count++;
  }
  //		cout<<"mesh points loaded into matrix..."<<endl;
  NewMeshPts=flirtmat*MeshPts;
  count=0;
  for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
    Pt newPt(NewMeshPts.element(0,count),NewMeshPts.element(1,count),NewMeshPts.element(2,count));
    (*i)->_update_coord = newPt;
    count++;
  }
  m.update();
	
	
}

void meshReg(Mesh & m, const vector< vector<float> > & flirtmatv){
	
	//refsize is actually target image
	int numPoints=m.nvertices();
	Matrix MeshPts(4,numPoints);
	Matrix NewMeshPts(4,numPoints);
	Matrix flirtmat(4,4);
	
	//read in flirt matrix, uses ascii
	for (int i=0; i<4;i++){
		for (int j=0; j<4;j++){
			
			flirtmat.element(i,j)=flirtmatv.at(i).at(j);
				cout<<flirtmat.element(i,j)<<" ";
		}
		cout<<endl;
	}
	flirtmat=flirtmat.i();
	
	
	//	cout<<"transform mesh points..."<<endl;
	int count=0;
	for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
		//	cout<<"count "<<count<<endl;	
		
		MeshPts.element(0,count)=(*i)->get_coord().X;
		MeshPts.element(1,count)=(*i)->get_coord().Y;
		MeshPts.element(2,count)=(*i)->get_coord().Z;
		MeshPts.element(3,count)=1;
		
		count++;
	}
	//		cout<<"mesh points loaded into matrix..."<<endl;
	NewMeshPts=flirtmat*MeshPts;
	count=0;
	for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
		Pt newPt(NewMeshPts.element(0,count),NewMeshPts.element(1,count),NewMeshPts.element(2,count));
		(*i)->_update_coord = newPt;
		count++;
	}
	m.update();
	
	
}

inline
ReturnMatrix unwrapMatrix(const Matrix & m) 
{
  ColumnVector munwrap(m.Nrows()*m.Ncols());
  unsigned int count=0;
  for (int i =0; i<m.Nrows() ; i++)
    for (int j =0; j<m.Ncols() ; j++,count++)
      munwrap.element(count)=m.element(i,j);
  return munwrap;
}

template<class T>
vector<T> vectorToVector( const Matrix & sm, const int & MaxModes){
  vector<T> vecM;	
  if (sm.Nrows()==1){
    for (int i=0;i<  MaxModes ; i++){
      vecM.push_back(static_cast<T>(sm.element(0,i)));
    }
  }else{
    for (int i=0;i<  MaxModes ; i++){
      vecM.push_back(static_cast<T>(sm.element(i,0)));
    }
		
		
  }
	
  return vecM;
}

template<class T>
vector<T> vectorToVector( const Matrix & sm){
  vector<T> vecM;	
  if (sm.Nrows()==1){
    for (int i=0;i< sm.Ncols() ; i++){
      vecM.push_back(static_cast<T>(sm.element(0,i)));
    }
  }else{
    for (int i=0;i< sm.Nrows() ; i++){
      vecM.push_back(static_cast<T>(sm.element(i,0)));
    }
		
		
  }
	
  return vecM;
}

template<class T>
vector< vector<T> > matrixToVector( const Matrix & sm, const int & MaxModes){
  vector< vector<T> > vecM;	
  for (int j=0;j< MaxModes ; j++){
    vector<T> mode;
    for (int i=0;i< sm.Nrows() ; i++){
      mode.push_back(sm.element(i,j));
    }
    vecM.push_back(mode);
  }   
	
  return vecM;
}

template<class T>
vector< vector<T> > matrixToVector( const Matrix & sm){
  vector< vector<T> > vecM;	
  for (int j=0;j< sm.Ncols() ; j++)
    {
      vector<T> mode;
      for (int i=0;i< sm.Nrows() ; i++)
	mode.push_back(static_cast<T>(sm.element(i,j)));
      vecM.push_back(mode);
    }   

  return vecM;
}



shapeModel* loadAndCreateShapeModel( const string & modelname)
{

  if (verbose.value()) cout<<"read model"<<endl;
  fslvtkIO* fmodel = new fslvtkIO(modelname,static_cast<fslvtkIO::DataType>(0));
  if (verbose.value()) cout<<"done reading model"<<endl;
  const int Npts=fmodel->getPointsAsMatrix().Nrows();
	
  unsigned int M = static_cast<unsigned int>(fmodel->getField("numSubjects").element(0,0));
  int MaxModes=M;
  if (verbose.value()) cout<<"setting up shape/appearance model"<<endl;

  //load mean shape into vector
  vector<float> Smean;
  Matrix* Pts = new Matrix;
  *Pts=fmodel->getPointsAsMatrix();

  if (verbose.value()) cout<<"The shape has "<<Npts<<" vertices."<<endl;
  for (int i=0;i< Npts ; i++){
    Smean.push_back(Pts->element(i,0));
    Smean.push_back(Pts->element(i,1));
    Smean.push_back(Pts->element(i,2));
  }   
  Pts->ReleaseAndDelete();

  //read polygon data
  vector< vector<unsigned int > > polygons = matrixToVector<unsigned int>(fmodel->getPolygons().t());

  //process shape modes and conditional intensity mean modes
  Matrix SmodesM;
  Matrix ImodesM;

	
  SmodesM=unwrapMatrix(fmodel->getField("mode0"));
  ImodesM=unwrapMatrix(fmodel->getField("Imode0"));
	
	
  for (int i =1; i<MaxModes;i++)
    {
      stringstream ss;
      ss<<i;
      string mode;
      ss>>mode;
      SmodesM=SmodesM | unwrapMatrix(fmodel->getField("mode"+mode));
      ImodesM=ImodesM | unwrapMatrix(fmodel->getField("Imode"+mode));
		
    }
  if (verbose.value()) cout<<MaxModes<<" modes of variation are retained."<<endl;

	
  vector< vector<float > > Smodes = matrixToVector<float>(SmodesM);
  vector< vector<float > > Imodes = matrixToVector<float>(ImodesM);
  ImodesM.Release();
  SmodesM.Release();
	
	
  //process rest of information, including intensity variance
  vector< vector<float > > Iprec = matrixToVector<float>(fmodel->getField("iCondPrec0").t());
  vector<float > Errs =  vectorToVector<float>(fmodel->getField("ErrPriors0"));
  vector<float > se =  vectorToVector<float>(fmodel->getField("eigenValues"), MaxModes);
  vector<float > ie =  vectorToVector<float>(fmodel->getField("iCondEigs0"));
  vector<float > Imean =  vectorToVector<float>(unwrapMatrix(fmodel->getField("Imean")));
  vector<int> labels =  vectorToVector<int>(fmodel->getField("labels"));


  if (verbose.value()) cout<<"The model was constructed from "<<M<<" training subjects."<<endl;

  //have read in all data and store in local structures, now delete the reader.
  delete fmodel;
	
  //create shape model
  shapeModel* model1 = new shapeModel(Smean, Smodes, se, Imean, Imodes,Iprec, ie, M,Errs,polygons,labels);

  return model1;

}



Mesh convertToMesh( const vector<float> & pts, const vector< vector<unsigned int> > & polys )
{
  Mesh mesh;

  mesh._points.clear();
  mesh._triangles.clear();
  
  int i=0;
  for (vector<float>::const_iterator p= pts.begin(); p!=pts.end(); p+=3,i++)
    {
      Mpoint * pt = new Mpoint(*p, *(p+1), *(p+2), i);
      mesh._points.push_back(pt);
    }
	
	
  for (vector< vector<unsigned int> >::const_iterator p=polys.begin(); p!=polys.end(); p++)
    {
      Triangle * tr = new Triangle( mesh._points.at( p->at(0) ), mesh._points.at( p->at(1) ), mesh._points.at( p->at(2) ));
      mesh._triangles.push_back(tr);
    }

  return mesh;

}

//%%%%%%%%%%%%%%mesh fill
void getBounds(Mesh m, int *bounds, float xdim, float ydim, float zdim){
	
  float xmin=1000,xmax=-1000,ymin=1000,ymax=-1000,zmin=1000,zmax=-1000;
  for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
    float tempx=(*i)->get_coord().X;
    float tempy=(*i)->get_coord().Y;
    float tempz=(*i)->get_coord().Z;
    if (tempx<xmin){
      xmin=tempx;
    }
    if (tempx>xmax){
      xmax=tempx;
    }
    if (tempy<ymin){
      ymin=tempy;
    }
    if (tempy>ymax){
      ymax=tempy;
    }
    if (tempz<zmin){
      zmin=tempz;
    }
    if (tempz>zmax){
      zmax=tempz;
    }
  }
  *bounds=static_cast<int>(floor(xmin/xdim)-1);
  *(bounds+1)=static_cast<int>(ceil(xmax/xdim)+1);
  *(bounds+2)=static_cast<int>(floor(ymin/ydim)-1);
  *(bounds+3)=static_cast<int>(ceil(ymax/ydim)+1);
  *(bounds+4)=static_cast<int>(floor(zmin/zdim)-1);
  *(bounds+5)=static_cast<int>(ceil(zmax/zdim)+1);
	
}




void draw_segment(volume<short>& image, const Pt& p1, const Pt& p2, int label)
{
  double xdim = (double) image.xdim();
  double ydim = (double) image.ydim();
  double zdim = (double) image.zdim();

  //in new version of bet2
  double mininc = min(xdim,min(ydim,zdim)) * .5;


	
  Vec n = (p1 - p2);
  double d = n.norm();
  n.normalize();
  //	double l = d*4;
  for (double i=0; i<=d; i+=mininc)
    {
      Pt p = p2 + i* n;
      image((int) floor((p.X)/xdim +.5),(int) floor((p.Y)/ydim +.5),(int) floor((p.Z)/zdim +.5)) = label;
    }
}



void draw_scalar_segment(volume<float>& image, volume<float>& nvol, const Pt& p1, const Pt& p2, float s1, float s2)
{
  double xdim = (double) image.xdim();
  double ydim = (double) image.ydim();
  double zdim = (double) image.zdim();

  //in new version of bet2
  double mininc = min(xdim,min(ydim,zdim)) * .25;


	
  Vec n = (p1 - p2);
  double d = n.norm();
  n.normalize();
  //	double l = d*4;
  for (double i=0; i<=d; i+=mininc)
    {
      Pt p = p2 + i* n;
      float interp_s = s2 + i/d*(s1-s2);
      image(MISCMATHS::round(p.X/xdim),MISCMATHS::round(p.Y/ydim),MISCMATHS::round(p.Z/zdim)) += interp_s;
      nvol(MISCMATHS::round(p.X/xdim),MISCMATHS::round(p.Y/ydim),MISCMATHS::round(p.Z/zdim)) += 1.0f;
    }
}


void set_scalars(Mesh& m, const ColumnVector& sc)
{
  if (sc.Nrows()!=m.nvertices()) {
    cerr << "ERROR in set_scalars:: mismatch in number of vertices ("<<m.nvertices()<<") and size of scalar vector ("<<sc.Nrows()<<")"<<endl;
    exit(EXIT_FAILURE);
  }
  int count=1;
  for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){  // loop over vertices
    (*i)->set_value(sc(count++));
  }	
  return;
}


volume<float> draw_scalar_mesh(const volume<float>& image, const Mesh &m, const ColumnVector& sc)
{

  Mesh localMesh;
  localMesh = m;
  if (sc.Nrows()>0) set_scalars(localMesh,sc);  // allow empty vector (zero size) to be passed in and values already in mesh to be used instead of scalars

  double xdim = (double) image.xdim();
  double ydim = (double) image.ydim();
  double zdim = (double) image.zdim();

  double mininc = min(xdim,min(ydim,zdim)) * .25;

  volume<float> res = image;
  res=0.0f;
  volume<float> nvol = res;

  for (list<Triangle*>::const_iterator i = localMesh._triangles.begin(); i!=localMesh._triangles.end(); i++)
    {
      Pt p1 = (*i)->get_vertice(1)->get_coord();
      Pt p2 = (*i)->get_vertice(2)->get_coord();
      Vec n = (*(*i)->get_vertice(0) - *(*i)->get_vertice(1));
      double d = n.norm();
      n.normalize();
	       
      float s0,s1,s2;
      s0=(*i)->get_vertice(0)->get_value();
      s1=(*i)->get_vertice(1)->get_value();
      s2=(*i)->get_vertice(2)->get_value();
		
      for (double j=0; j<=d ;  j+=mininc)
	{
	  Pt p = p1  + (double)j* n;
	  float s = s1 + j/d*(s0-s1);  
	  draw_scalar_segment(res, nvol, p, p2, s, s2);
	} 
    }

  res=divide(res,nvol,nvol);
  return res;
}


volume<float> draw_scalar_mesh(const volume<float>& image, const Mesh &m)
{
  ColumnVector dummy;
  return draw_scalar_mesh(image,m,dummy);
}


volume<short> draw_mesh(const volume<short>& image, const Mesh &m, int label)
{

  double xdim = (double) image.xdim();
  double ydim = (double) image.ydim();
  double zdim = (double) image.zdim();

  //in new version of bet2
  double mininc = min(xdim,min(ydim,zdim)) * .5;


  volume<short> res = image;
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
  return res;
}

volume<short> make_mask_from_meshInOut(const volume<float> & image, const Mesh& m, int label, int* bounds)
{
	
  float xdim = (float) image.xdim();
  float ydim = (float) image.ydim();
  float zdim = (float) image.zdim();
	
  volume<short> mask;
  copyconvert(image,mask);
	
	
  mask = 0;
  mask = draw_mesh(mask, m,label+100);
	
	
	
  // THIS EXCLUDES THE ACTUAL MESH
  volume<short> otl=mask;
  getBounds(m,bounds,xdim,ydim,zdim);
  vector<Pt> current;
  current.clear();
  Pt c(bounds[0]-2, bounds[2]-2, bounds[4]-2);
		
  mask.value(static_cast<int>(c.X),static_cast<int>(c.Y),static_cast<int>(c.Z)) = label;
  current.push_back(c);
  int fillCount=0;
  while (!current.empty())
    {
      Pt pc = current.back();
      int x, y, z;
      x=(int) pc.X; y=(int) pc.Y; z=(int) pc.Z;
			
      current.pop_back();
      fillCount++;
			
			
      if (bounds[0]<=x-1 && mask.value(x-1, y, z)==0) {
	mask.value(x-1, y, z) = label;
	current.push_back(Pt(x-1, y, z));
      }
      if (bounds[2]<=y-1 && mask.value(x, y-1, z)==0) {
	mask.value(x, y-1, z) = label;
	current.push_back(Pt(x, y-1, z));
      }
      if (bounds[4]<=z-1 && mask.value(x, y, z-1)==0) {
	mask.value(x, y, z-1) = label;
	current.push_back(Pt(x, y, z-1));
      }
      if (bounds[1]>=x+1 && mask.value(x+1, y, z)==0){
	mask.value(x+1, y, z) = label;
	current.push_back(Pt(x+1, y, z));
      }
      if (bounds[3]>=y+1 && mask.value(x, y+1, z)==0){
	mask.value(x, y+1, z) = label;
	current.push_back(Pt(x, y+1, z));
      }
      if (bounds[5]>=z+1 && mask.value(x, y, z+1)==0){
	mask.value(x, y, z+1) = label;
	current.push_back(Pt(x, y, z+1)); 
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
  return otl;
}

//*************These are the overlap measures function**************************
int findLabel(int label1,vector<int>* vlabels){
  for (unsigned int i=0;i<vlabels->size();i++){
    if (vlabels->at(i)==label1){
      return i;
    }
  }
  return -1;
}

bool findAddLabel(int label1,int label2,int* indseg, vector<int>* vlabels, vector<int>* vTP, vector<int>* vFN, vector<int>* vFP, vector<int>* segImLabels,vector<int>* minInterX,vector<int>* maxInterX,vector<int>* minInterY,vector<int>* maxInterY,vector<int>* minInterZ,vector<int>* maxInterZ){
  int ind2=-1;//ind1=-1, ind2=-1;
  int ind1=-1;
  //return 1 signifies intersection
	
  *indseg=-1;
  for (unsigned int i=0;i<segImLabels->size();i++){
    if ((segImLabels->at(i)==label1)&&(label1!=0)){
      *indseg=i;
      i=segImLabels->size()+1;
    }
  } 
  if ((*indseg==-1)&&(label1!=0)){
		
    segImLabels->push_back(label1);
    minInterX->push_back(10000);
    maxInterX->push_back(0);
    minInterY->push_back(10000);
    maxInterY->push_back(0);
    minInterZ->push_back(10000);
    maxInterZ->push_back(0);
    *indseg=segImLabels->size()-1;
  }
  for (unsigned int i=0;i<vlabels->size();i++){	
    // cout<<"labels "<<label1<<" "<<label2<<endl;
    if ((vlabels->at(i)==label1)|(vlabels->at(i)==label2)){ 
      if (label1==label2){
	//	cout<<"found ind1=ind2"<<endl;
	ind1=ind2=i;
	i=vlabels->size()+1;
      }else{
	if (label1==vlabels->at(i)){
	  ind1=i;
	}else{
	  ind2=i;
	}
      }
    }
    if(i==vlabels->size()-1){
      // cout<<"do i enter?"<<endl;
      if (ind1==-1){
	//  cout<<"end and ind1=-1 "<<label1<<" "<<label2<<endl;
	ind1=vlabels->size();
	vlabels->push_back(label1);
	vTP->push_back(0);
	vFN->push_back(0);
	vFP->push_back(0);
      }
      if (ind2==-1){
	//  cout<<"end and ind2=-1 "<<label1<<" "<<label2<<endl;
	if (label1==label2){
	  //	cout<<"end and lab1=lab2"<<endl;
	  ind2=ind1;
	}else{
	  //		cout<<"end and lab1!=lab2"<<endl;
	  ind2=vlabels->size();
	  vlabels->push_back(label2);
	  vTP->push_back(0);
	  vFN->push_back(0);
	  vFP->push_back(0);
	}
				
      }
    }
  }
	
  //	cout<<"out of loop"<<endl;
	
  if (label1==label2){
    vTP->at(ind1)+=1;
    return true;
  }else{
    vFP->at(ind1)+=1;
    vFN->at(ind2)+=1;
    return false;
  }
	
  //vCount->at(*ind1)+=1;
	
	
}

Matrix overlaps(const volume<short> segIm, const volume<short> gold){
	
  int sizex= segIm.xsize();
  int sizey=segIm.ysize();
  int sizez=segIm.zsize();
  bool inter;
  int indseg;
  //find union and intersection
  vector<int> vlabels,vTP, vFN, vFP,segLabels;
  //these are used to speed up the distance calculation
  vector<int> minInterX,maxInterX,minInterY,maxInterY,minInterZ,maxInterZ;
  vlabels.push_back(0);
  vTP.push_back(0);
  vFN.push_back(0);
  vFP.push_back(0);
  for (int k=0;k<sizez;k++){
    for (int j=0;j<sizey;j++){ 
      for (int i= 0; i<sizex;i++){
	//    cout<<"prefindlabel "<<i<<" "<<j<<" "<<k<<endl;
	inter=findAddLabel(segIm.value(i,j,k), gold.value(i,j,k),&indseg,&vlabels,&vTP, &vFN, &vFP,&segLabels,&minInterX,&maxInterX,&minInterY,&maxInterY,&minInterZ,&maxInterZ);
	//      cout<<"postfindlabel "<<ind1<<endl;
	if ((inter)&&(segIm.value(i,j,k)!=0)){
	  //	cout<<minInterX.size()<<" "<<segIm.value(i,j,k)<<" "<<gold.value(i,j,k)<<endl;
	  if (i<minInterX.at(indseg)){
	    minInterX.at(indseg)=i;
	  }
	  if (j<minInterY.at(indseg)){
	    minInterY.at(indseg)=j;
	  }
	  if (k<minInterZ.at(indseg)){
	    minInterZ.at(indseg)=k;
	  }
	  if (i>maxInterX.at(indseg)){
	    maxInterX.at(indseg)=i;
	  }
	  if (j>maxInterY.at(indseg)){
	    maxInterY.at(indseg)=j;
	  }
	  if (k>maxInterZ.at(indseg)){
	    maxInterZ.at(indseg)=k;
	  }
	  //	cout<<"leave min max update "<<endl;
	}
      }
    }
  }
	
	
  //this overlap does not give distance weighted stuff
  Matrix simMeasures(static_cast<int>(segLabels.size()),9);
	
  for (unsigned int i=0; i<segLabels.size();i++){
    //cout<<"add measures "<<endl;
    //for each label in segmentation
    //sum intersections
    int ind=findLabel(segLabels.at(i),&vlabels);

    simMeasures.element(i,0)=vlabels.at(ind);
    simMeasures.element(i,1)=vTP.at(ind);
    simMeasures.element(i,2)=vFP.at(ind);
    simMeasures.element(i,3)=vFN.at(ind);
    simMeasures.element(i,4)=static_cast<float>(vTP.at(ind))/(vTP.at(ind)+vFP.at(ind)+vFN.at(ind));
    simMeasures.element(i,5)=2.0*vTP.at(ind)/(2*vTP.at(ind)+vFP.at(ind)+vFN.at(ind));
    simMeasures.element(i,6)=static_cast<float>(vFN.at(ind))/(vFN.at(ind)+vTP.at(ind));
    simMeasures.element(i,7)=static_cast<float>(vFP.at(ind))/(vFN.at(ind)+vTP.at(ind));
    simMeasures.element(i,8)=static_cast<float>(vTP.at(ind))/(vFN.at(ind)+vTP.at(ind));
    cout<<"Dice "<<vlabels.at(ind)<<" "<<2.0*vTP.at(ind)/(2*vTP.at(ind)+vFP.at(ind)+vFN.at(ind))<<endl;
  }
  return simMeasures;
}


//****************************BVARS I/O**************************************

string read_bvars(string fname,Matrix & bvars, vector<string> & vnames, vector<int>& vnvars,string impath, vector< vector< vector<float> > > & fmatv){
	string stemp;
	string modelNames;
	int N;//number of subjects
	fmatv.clear();
	ifstream fin;
	fin.open(fname.c_str());
	//throw away first three lines 
	getline(fin,stemp);//this is bvars file
	getline(fin,modelNames);//modelnames;
	fin>>stemp>>N;
	vnvars.clear();
	//int NmaxBvars=0;
	//vector< vector<float> > all_bvars;
	
	for (int i=0; i<N;i++){
		vector<float> vbvars;

		fin>>stemp;//read in subject id
		vnames.push_back(impath+stemp);
		int nvars;//how many vars written for the subject
		fin>>nvars;

		vnvars.push_back(nvars);
		char blank;
		fin.read(reinterpret_cast<char*>(&blank),sizeof(char));

		if (i==0){
			bvars.ReSize(nvars,N);
		}
		vnvars.push_back(nvars);
		for (int j=0;j<nvars;j++){
			if (j<nvars){
				float ftemp;
			
			//	fin>>ftemp;
				fin.read(reinterpret_cast<char*>(&ftemp),sizeof(float));

				bvars.element(j,i)=ftemp;
			}else{
				bvars.element(j,i)=0;
			}
		}
	//	all_bvars.push_back(vbvars);
		vector< vector<float> > mat;
		for (int k=0;k<4;k++)
		{
			vector<float> row;
			for (int l=0;l<4;l++)
			{
				float ftemp;
			//	fin>>ftemp;
				fin.read(reinterpret_cast<char*>(&ftemp),sizeof(float));

				row.push_back(ftemp);
			}
			mat.push_back(row);
		}
		fmatv.push_back(mat);
	}
	

	return modelNames;
}

void write_bvars(string fname,string modelname,Matrix bvars, int numModes,vector<string> vnames){
  ofstream fout;
	
  fout.open(fname.c_str());
  fout<<"this is a bvars file"<<endl; 
  fout<<modelname<<endl;
  fout<<"NumberOfSubjects "<<bvars.Nrows()<<endl;
	
  for (int i=0;i<bvars.Nrows();i++){
    fout<<vnames.at(i)<<" ";
    fout<<numModes<<" ";
#ifdef PPC64
    int n=0;
#endif
    for (int j=0;j<bvars.Ncols();j++){
      fout<<bvars.element(i,j)<<" ";
#ifdef PPC64
      if ((n++ % 50) == 0) fout.flush();
#endif
    }
    fout<<endl;
  }
	
  fout<<endl;
  fout.close();
}

string read_bvars_ModelName(string fname){
  string stemp;
  string modelNames;
  ifstream fin;
  fin.open(fname.c_str());
  //throw away first three lines 
  getline(fin,stemp);//this is bvars file
  getline(fin,modelNames);//modelnames

  if ( access(modelNames.c_str(),F_OK)!=0 ) {
    string FslDir(getenv("FSLDIR"));
    string modelNameReplaced(modelNames);
    modelNameReplaced.replace(0,14,FslDir);
    cerr << "Warning: " << modelNames << " does not exist. Attempting to switch to: " << modelNameReplaced << endl;
    modelNames=modelNameReplaced;
  }
			
  return modelNames;
}

//****************************END OF BVARS I/O**************************************


//****************************BOUNDARY CORRECTION FUNCTIONS**************************************
float mode(vector<float> vdists, float min, float max){
  int N=static_cast<int>(vdists.size());
  float bins=128;
	
  float binwidth=(max-min)/bins;
  vector<int> bincounts;
  //initialize bincounts to zero
  for (int b=0;b<bins;b++){
    bincounts.push_back(0);
  }
	
  for (int i=0;i<N;i++){
    //search through each bin
    for (int b=0;b<bins;b++){
      if (vdists.at(i)<(min+(b+1)*binwidth)){
	bincounts.at(b)++;
	break;
      }
    }
  }
	
  //search for max bin count
  int maxcount=0;
  int maxind=0;
  for (int b=0;b<bins;b++){
    //cout<<bincounts.at(b)<<endl;
    if (bincounts.at(b)>maxcount){
      maxcount=bincounts.at(b);
      maxind=b;
    }
		
  }
	
  return (min+maxind*binwidth+binwidth/2.0);
}
float mode(vector<float> vdists,int *maxcount){
  //	cout<<"calc mode"<<endl;
  int N=static_cast<int>(vdists.size());
  *maxcount=0;
  float maxlowint=0;
  float bins=256;
  bins=128;
  float binwidth=(vdists.at(N-1)-0)/bins;
  int count=0;	
  int bincount=1;
  float lowint=0;
  for (int i=0;i<N;i++){
		
    //cout<<"imode "<<i<<" "<<N<<endl;
    if (vdists.at(i)<lowint+binwidth){
      count++;
      if (count>(*maxcount)){
	*maxcount=count;
	maxlowint=lowint;
      }
    }else{
      count=0;
      i--;
      bincount++;
      lowint=bincount*binwidth;
    }
		
  }
	
  return (maxlowint+binwidth/2.0);
}	
float fullwidthhalfmax(vector<float> vdists,float halfmaxval, float *halfmin,float *halfmax){
  int N=static_cast<int>(vdists.size());
  int maxcount=0;
  float maxlowint=0;
  float bins=256;
  bins=128;
  float binwidth=(vdists.at(N-1)-0)/bins;
  int countprev=0;
	
  bool foundmin=false;
	
  int count=0;	
  int bincount=1;
  float lowint=0;
  for (int i=0;i<N;i++){
		
		
    if (vdists.at(i)<lowint+binwidth){
      count++;
      if (count>(maxcount)){
	maxcount=count;
	maxlowint=lowint;
      }
    }else{
      countprev=count;
      if ((count>halfmaxval)&&(!foundmin)){
	*halfmin=lowint+binwidth/2.0;
	foundmin=true;
      }
      if ((count>halfmaxval)){
	*halfmax=lowint+binwidth/2.0;
      }
      count=0;
      i--;
      bincount++;
      lowint=bincount*binwidth;
    }
		
  }
	
  return (maxlowint+binwidth/2.0);
}

float boundaryCorr(volume<short>* mask, volume<float>* ref, int label, float zthresh, int* bounds){
  //returns volume
  //build intensity vector 
  vector<float> vgraylevels;
  vector <float>::iterator Iter;
  float dist=10000;
	
  for (int i=bounds[0];i<bounds[1];i++){
		
    for (int j=bounds[2];j<bounds[3];j++){
      for (int k=bounds[4];k<bounds[5];k++){		
				
	if (mask->value(i,j,k)==label){
	  dist=ref->value(i,j,k);
					
	  if (vgraylevels.empty()){
	    vgraylevels.push_back(dist);
	  }else if (dist>=vgraylevels.back()){
	    vgraylevels.push_back(dist);
	  }else {
	    for ( Iter = vgraylevels.begin( ) ; Iter !=vgraylevels.end( ) ; Iter++ ){
	      if (dist<*Iter){
								
		vgraylevels.insert(Iter,dist);
		break;
	      }
							
	    }
						
						
	  }
					
	}
      }
    }
  }
  int maxcount;
  //don't end up using the mode...re-calculate centre from fullwidth half-maximum
  mode(vgraylevels,&maxcount);
  //now find full width half maximum
  //float halfmaxval=maxcount/2.0;
  float halfmin,halfmax;
  fullwidthhalfmax(vgraylevels,maxcount/2.0,&halfmin, &halfmax);
  float mean=(halfmin+halfmax)/2;
  float sdev=abs(halfmax-halfmin)/2.35;
	
  //test for thalamus
  //	volume<float> zvol;
  //	copyconvert(*mask, zvol);
  //	zvol=0;

	
  float vol=0;

  float min=0, max=0;
  //calculates z-value for all lgive a nice place to initialize EM or could act directly on intensity and use range to initialize
  for (int i=bounds[0];i<bounds[1];i++){
    for (int j=bounds[2];j<bounds[3];j++){
      for (int k=bounds[4];k<bounds[5];k++){
	float z=0.0;
	if (mask->value(i,j,k)==(label+100  )){
	  z=(ref->value(i,j,k)-mean)/sdev;
	  //					zvol.value(i,j,k)=(z);
	  if (zthresh>=0){
	    if (abs(z)>zthresh){
	      mask->value(i,j,k)=0;
	    }else{
	      mask->value(i,j,k)=label;
	      vol++;			
	    }
	  }
	}else if (mask->value(i,j,k)==label){
	  vol++;
	  z=(ref->value(i,j,k)-mean)/sdev;
	  //					zvol.value(i,j,k)=(z);
	  if (z>max){ max=z; }
	  if (z<min){ min=z; }
	}
      }
    }
  }
	
  //	save_volume(zvol,outname.value());

  return vol;
}

int findStructLabel(volume<short>* mask, int* bounds){
  //used to find label and set bounds
  int xmin=10000,ymin=10000,zmin=10000;
  int xmax=-1, ymax=-1, zmax=-1;
  int label=999;
  bool found=false;
  for (int i=bounds[0];i<bounds[1];i++){
    for (int j=bounds[2];j<bounds[3];j++){
      for (int k=bounds[4];k<bounds[5];k++){		
	if (mask->value(i,j,k)>0){
	  if (xmin>i){ xmin=i; }
	  if (ymin>j){ ymin=j; }
	  if (zmin>k){ zmin=k; }
	  if (xmax<i){ xmax=i; }
	  if (ymax<j){ ymax=j; }
	  if (zmax<k){ zmax=k; }
					
	}
	if ((mask->value(i,j,k)<100)&&(mask->value(i,j,k)!=0)&&(!found)){
	  label=mask->value(i,j,k);
	  found=true;
	  //	break;
	}
      }
      //	if (found){break;}
    }
    //if (found){break;}
  }
  bounds[0]=xmin-1;
  bounds[1]=xmax+1;
  bounds[2]=ymin-1;
  bounds[3]=ymax+1;
  bounds[4]=zmin-1;
  bounds[5]=zmax+1;
  return label;
}


//****************************END OF BOUNDARY CORRECTION FUNCTIONS**************************************

//********************************GLM for stats******************************//
ColumnVector GLM_fit(Matrix G, Matrix D, ColumnVector contrast){
	
  //start for only well conditioned design matrices
  Matrix A=G.t()*G;
  Matrix Betas(D.Nrows(),D.Ncols());
  Betas=A.i()*G.t()*D;
  Matrix Mres;
  Mres=D-G*(Betas);
  //calculate residual variance
  ColumnVector avgRes(D.Ncols());
  for (int i=0; i<D.Ncols();i++){
    avgRes.element(i)=((Mres.SubMatrix(1,Mres.Nrows(),i+1,i+1)).t()*(Mres.SubMatrix(1,Mres.Nrows(),i+1,i+1))).AsScalar()/(G.Nrows()-G.Ncols());
  } 
  //convert to standard error
  avgRes=avgRes*(contrast.t()*A.i()*contrast).AsScalar();	
  Matrix test;
  test=contrast.t()*Betas.SubMatrix(1,Betas.Nrows(),1,1)*Betas.SubMatrix(1,Betas.Nrows(),1,1).t()*contrast;
  ColumnVector tstats(avgRes.Nrows());
  for (int i=0; i<avgRes.Nrows();i++){
    tstats.element(i)=(contrast.t()*Betas.SubMatrix(1,Betas.Nrows(),i+1,i+1)/sqrt(avgRes.element(i))).AsScalar();
  }
  return tstats;
}


float MVGLM_fit(Matrix G, Matrix D, Matrix contrast, int& df1, int& df2){
	
  // conversion to "normal" notation is:
  //  G -> X = design matrix (N_subj by N_evs)
  //  D -> Y = data matrix (N_subj by 3)
	
  //Calculate estimated values
  Matrix Yhat=G*(G.t()*G).i()*G.t()*D;
  //calculate E covariance matrix
  Matrix E=D-Yhat;
  E=E.t()*E;
//	cout<<"hmm"<<G.t()*G<<endl<<contrast*G.t()*G*contrast.t()<<endl;
  //calculate H, the sum-square /cross square product for hypothesis test
  Matrix YhatH= G*contrast.t()*(contrast*G.t()*G*contrast.t()).i()*contrast*G.t()*D;
  Matrix H=D-YhatH;
  //not efficient but easy to convert to other statistics
  H=H.t()*H-E;

  // Calculate Pillai's Trace
  
  int N=D.Nrows();//number of samples 
  int p=D.Ncols();//number of dimensions (3 for vertex coordinates)
  int N_R=G.Ncols();//number of regressors
  int N_C=contrast.Nrows(); //number of rows in contrast matrix
  
  int v_h = N_R - N_C;
  int v_e = N-N_R;	
  
  int s=Min(p,v_h);
  float t = (abs(p-v_h)-1)/2.0;
  float u = (v_e-p-1)/2.0;
  df1=MISCMATHS::round(s*(2*t+s+1));
  df2=MISCMATHS::round(s*(2*u+s+1));
  
  float pillai=(H*(H+E).i()).Trace();
  float F=0;
  F=(pillai/(s-pillai))*(df2/df1);
  
  if (verbose.value()){
    cout<<"Pillai F "<<pillai<<" "<<F<<" "<<df1<<" "<<df2<<endl;
  }
  
  return F;
}


float MVGLM_fit(Matrix G, Matrix D, Matrix contrast)
{
  int df1, df2;
  float retval=MVGLM_fit(G,D,contrast,df1,df2);
  return retval;
}


//******************************EXECUTION FUNCTIONS******************************


void do_work_SingleClean(){
  //this function is working directly on volumes
  volume<float> t1im;
  volume<short> segim;
  int bounds[6]={0,0,0,0,0,0};
  read_volume(t1im,refname.value());
  read_volume(segim,inname.value());
  //FIND LABEL AND BOUNDS FOR EACH IMAGE
  //need to reset lower bounds as well
  bounds[0]=0;
  bounds[2]=0;
  bounds[4]=0;
  bounds[1]=segim.xsize();
  bounds[3]=segim.ysize();
  bounds[5]=segim.zsize();
  int label=findStructLabel(&segim, bounds);
  volume<float> vz1;

  boundaryCorr(&segim, &t1im,label, thresh.value(), bounds);
  segim.setAuxFile("MGH-Subcortical");  // for automatic colormap
  save_volume(segim,outname.value());
	
}

//*****************************LINEA TRANSFORM********************************************************//
Matrix rigid_linear_xfm(Matrix Data,ColumnVector meanm, Mesh mesh, bool writeToFile){
  //ColumnVector avgM(sub.Ncols());
  //determine translations
  int Nsub=Data.Ncols();
  int Npoints=Data.Nrows()/3;
	
	
  //***********CALCULATE CENTROIDS*************//
  //calculate centroid of mean mesh
  float Mxr=0,Myr=0,Mzr=0;
  for (int i=0;i<meanm.Nrows();i=i+3){
    Mxr+=meanm.element(i);
    Myr+=meanm.element(i+1);
    Mzr+=meanm.element(i+2);
  }
  Mxr/=(meanm.Nrows()/3);
  Myr/=(meanm.Nrows()/3);
  Mzr/=(meanm.Nrows()/3);
	
  //calculate centroid 
  vector<float> vMx,vMy,vMz;
  for (int i=0;i<Nsub;i++){
    float sx=0,sy=0,sz=0;
    //cout<<" i "<<i<<endl;
    for (int j=0;j<Data.Nrows();j=j+3){
      sx+=Data.element(j,i);
      sy+=Data.element(j+1,i);
      sz+=Data.element(j+2,i);
    }
    vMx.push_back(sx/Npoints);
    vMy.push_back(sy/Npoints);
    vMz.push_back(sz/Npoints);
  }
  vector< Matrix > vR;
  vector< float > vscale;
	
  for (int subject=0;subject<Data.Ncols();subject++){
    //***********Demean data*************//
    //Demean the Data and reformat
    //Demean Mean mesh and reformat
    Matrix DataDM(3,Data.Nrows()/3);
    Matrix RefDM(3,Data.Nrows()/3);
    for (int i=0;i<Data.Nrows();i=i+3){
      DataDM.element(0,i/3)=Data.element(i,subject)-vMx.at(subject);
      DataDM.element(1,i/3)=Data.element(i+1,subject)-vMy.at(subject);
      DataDM.element(2,i/3)=Data.element(i+2,subject)-vMz.at(subject);
      RefDM.element(0,i/3)=meanm.element(i)-Mxr;
      RefDM.element(1,i/3)=meanm.element(i+1)-Myr;
      RefDM.element(2,i/3)=meanm.element(i+2)-Mzr;
    }
		
    //*************This includes scale calculation ***********
    float scale=1.0;
    if (useScale.value()){
      //***********calculate scale*************//
      //Data is right (in Horn), Ref is left
      //calculate length vectors and sum
      float sumr=0, suml=0;
      for (int i=0;i<DataDM.Ncols();i++){
	suml+=(DataDM.element(0,i)*DataDM.element(0,i) + DataDM.element(1,i)*DataDM.element(1,i) +DataDM.element(2,i)*DataDM.element(2,i) );
	sumr+=(RefDM.element(0,i)*RefDM.element(0,i) + RefDM.element(1,i)*RefDM.element(1,i) +RefDM.element(2,i)*RefDM.element(2,i) );
      }
      scale=sqrt(sumr/suml);
    }
    vscale.push_back(scale);
    //**********end of scale calculation**************
    //***********calculate rotattions*************//
    //cout<<"reshaped matrices"<<endl;
    Matrix M=RefDM*(DataDM.t());
    Matrix U;
    DiagonalMatrix D;
    SVD(M.t()*M,D,U);
    //M should always be a 3x3 matrix
    for (int i=0;i<D.Nrows();i++){
      D.element(i)=1/sqrt(D.element(i));
    }
		
    Matrix R(3,3);
    R=M*(U*D*U.t());
    vR.push_back(R);
  }
	
  //*****************APPLY TRANSFORMATION TO MESHES**********************//
  //NO SCALE IS CALCULATED
	
	
  Matrix DataNew(Data.Nrows(),Data.Ncols());
  for (int subject=0;subject<Data.Ncols();subject++){
    //	cout<<"subject "<<subject<<endl;
		
    //reshape data
    Matrix DataRS(3,Data.Nrows()/3);
    Matrix RefMean(3,Data.Nrows()/3);
    Matrix DataMean(3,Data.Nrows()/3);
    for (int i=0;i<Data.Nrows();i=i+3){
      //cout<<i/3<<" "<<DataRS.element(i,subject)<<endl;
      DataRS.element(0,i/3)=Data.element(i,subject);
      DataRS.element(1,i/3)=Data.element(i+1,subject);
      DataRS.element(2,i/3)=Data.element(i+2,subject);
      //		if (subject==4){
      //			cout<<i/3<<" "<<DataRS.element(0,i/3)<<endl;
      //		}
      RefMean.element(0,i/3)=Mxr;
      RefMean.element(1,i/3)=Myr;
      RefMean.element(2,i/3)=Mzr;
      DataMean.element(0,i/3)=vMx.at(subject);
      DataMean.element(1,i/3)=vMy.at(subject);
      DataMean.element(2,i/3)=vMz.at(subject);
    }

    // TODO TODO process the normal vectors through this too and formulate a new mean
    ///APPLY TRANSFORMATION
    Matrix Reg(3,Data.Nrows()/3);//=R*DataRS+(RefMean-R*DataMean);
    //difference between 6dof and 7dof
    if (useScale.value()){
      //cout<<"Scale "<<vscale.at(subject)<<endl;
      Reg=vscale.at(subject)*vR.at(subject)*DataRS+(RefMean-vscale.at(subject)*vR.at(subject)*DataMean);
    }else{
      Reg=vR.at(subject)*DataRS+(RefMean-vR.at(subject)*DataMean);
      //Reg=R*DataRS+(RefMean-R*DataMean);
				
    }
			
    for (int i=0;i<Reg.Ncols();i++){
      DataNew.element(3*i,subject)=Reg.element(0,i);
      DataNew.element(3*i+1,subject)=Reg.element(1,i);
      DataNew.element(3*i+2,subject)=Reg.element(2,i);
    }
		
			
    Mesh m=mesh;
    int count=0;
    for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
      (*i)->_update_coord.X=DataNew.element(count,subject);
      (*i)->_update_coord.Y=DataNew.element(count+1,subject);
      (*i)->_update_coord.Z=DataNew.element(count+2,subject);
      count+=3;
    }	
    m.update();
				
//    	string snum;
 //   stringstream ssnum;
   // ssnum<<subject;
  //  ssnum>>snum;
  //  m.save(snum+"reg.vtk",3);
  }
	
  return DataNew;			

}


Mesh setup_avMesh(const Mesh& m, const Matrix& MeshVerts) {
  // input m must have the right topology (i.e. triangles) and MeshVerts must contain the coordinates in the first column
  //   with order of iterated points in m and the rows in MeshVerts corresponding (by a factor of 3 as X,Y,Z are stored sequentially for each point)
  Mesh avMesh;
  Matrix avVerts;
  avVerts=mean(MeshVerts,2);  // form the average coordinate location across the whole group
  avMesh=m;
  int count=0;
  for (vector<Mpoint*>::iterator i = avMesh._points.begin(); i!=avMesh._points.end(); i++ ){  // loop over vertices
    Pt newpt;
    newpt.X = avVerts.element(count,0);
    newpt.Y = avVerts.element(count+1,0);
    newpt.Z = avVerts.element(count+2,0);
    (*i)->_update_coord = newpt;  // <MJrant>  clunky way of doing what should be: set_coord(const Pt &newpt)   </MJrant>
    (*i)->update();
    count+=3;
  }	
  return avMesh;
}


Matrix calc_avNormals(Mesh& avMesh) {  
  // TODO - would like this to be const but can't get it to compile like that!
  // TODO - i.e. Matrix calc_avNormals(Mesh& const avMesh) {  

  // input m must have the right topology (i.e. triangles) and MeshVerts must contain the coordinates in the first column
  //   with order of iterated points in m and the rows in MeshVerts corresponding (by a factor of 3 as X,Y,Z are stored sequentially for each point)
  Matrix avNormals(avMesh.nvertices()*3,1);
  int count=0;
  for (vector<Mpoint*>::iterator i = avMesh._points.begin(); i!=avMesh._points.end(); i++ ){  // loop over vertices
    Vec norm;
    norm=(*i)->local_normal();
    avNormals.element(count,0) = norm.X;
    avNormals.element(count+1,0) = norm.Y;
    avNormals.element(count+2,0) = norm.Z;
    count+=3;
  }	
  return avNormals;
}


Matrix recon_meshesMNI( shapeModel* model1, Matrix bvars, ColumnVector* meanm, Mesh * meshout,vector<Mesh>* vMeshes, Mesh& avMesh){
  // note that meshout may be the same as avMesh  (MJ:not yet verified)

  vMeshes->clear();
  //want to return mean mesh in vector form
  meanm->ReSize(model1->smean.size());
  {
		
    int count=0;
    int cumnum=0; 
    for (int sh=0; sh<1;sh++)
      {
	Mesh m= convertToMesh(model1->smean, model1->cells);//model1->getTranslatedMesh(sh);
	*meshout=m;
	for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
	  meanm->element(3*cumnum+count)=(*i)->get_coord().X;
	  meanm->element(3*cumnum+count+1)=(*i)->get_coord().Y;
	  meanm->element(3*cumnum+count+2)=(*i)->get_coord().Z;
	  count+=3;
	}
	cumnum+=model1->smean.size();	
      }
			
  }


  //need number of subjects (to determine number of modes)
  //int M=model1->getNumberOfSubjects();
  int Tpts=model1->smean.size()/3;
  Matrix MeshVerts(3*Tpts,bvars.Ncols());
  Mesh m;
  //this is different than mnumber of subjects which were used to create model
  //int numSubs=bvars.Ncols();
  //this loads all mesh point into a matrix to do stats on....
  //cout<<"generate and load vertices into matrix "<<endl;
  int NumSubjects=bvars.Ncols();
  for (int j=0;j<NumSubjects;j++){  // loop over subjects
    //for each subject
    vector<float> vars;
    for (int i=0; i<bvars.Nrows();i++){
      vars.push_back(bvars.element(i,j));
    }
    //keep track of number of points preceding
    int cumnum=0;
    for (int sh=0; sh<1;sh++){  // an ex-loop?  Pining for the fjords?
      //cout<<"cumnum "<<cumnum<<endl;
      m=convertToMesh(model1->getDeformedGrid(vars),model1->cells);	
      vMeshes->push_back(m);
      int count=0;
      for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){  // loop over vertices
	MeshVerts.element(3*cumnum+count,j)=(*i)->get_coord().X;
	MeshVerts.element(3*cumnum+count+1,j)=(*i)->get_coord().Y;
	MeshVerts.element(3*cumnum+count+2,j)=(*i)->get_coord().Z;
	count+=3;
      }	
      cumnum+=model1->smean.size();	
    }	
  }
  avMesh = setup_avMesh(m,MeshVerts);
  return MeshVerts;
}


Matrix recon_meshesNative( const string & modelname, const Matrix & bvars, ColumnVector & meanm, Mesh & meshout, const vector<string> & subjectnames, const vector< vector< vector<float> > > & flirtmats,vector<Mesh> & vMeshes, Mesh& avMesh){
  vMeshes.clear();
  //fslvtkIO* fin = new fslvtkIO(modelname,static_cast<fslvtkIO::DataType>(0));
  shapeModel* model1=loadAndCreateShapeModel(modelname);

  //	shapeModel* model1= new shapeModel();
  //	model1->setImageParameters(refXsize,refYsize, refZsize,refXdim, refYdim, refZdim);
  //	cout<<"load model "<<modelname<<endl;
  //	model1->load_bmv_binaryInfo(modelname,1);
  //	model1->load_bmv_binary(modelname,1);
	
  //want to return mean mesh in vector form
  meanm.ReSize(model1->smean.size());
  {
		
    int count=0;
    int cumnum=0; 
    for (int sh=0; sh<1;sh++){
      //difference between origin specification
      Mesh m=convertToMesh( model1->smean , model1->cells );
      //			m.save("targetMesh.vtk",3);
      //Mesh m=model1->getShapeMesh(sh);
      meshout=m;
      for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
	meanm.element(3*cumnum+count)=(*i)->get_coord().X;
	meanm.element(3*cumnum+count+1)=(*i)->get_coord().Y;
	meanm.element(3*cumnum+count+2)=(*i)->get_coord().Z;
	//				cout<<"coord "<<meanm->element(3*cumnum+count)<<" "<<meanm->element(3*cumnum+count+1)<<" "<<meanm->element(3*cumnum+count+2)<<endl;

	count+=3;
      }
      cumnum+=model1->smean.size()/3;	
    }
		
  }
	
  //need number of subjects (to determine number of modes)
  int Tpts=model1->smean.size()/3;
  Matrix MeshVerts(3*Tpts,bvars.Ncols());
  //this is different than mnumber of subjects which were used to create model
  //this loads all mesh point into a matrix to do stats on....
  //	cout<<"generate and load vertices into matrix "<<endl;
  int NumSubjects=bvars.Ncols();
  Mesh m;
  for (int j=0;j<NumSubjects;j++){  // loop over subjects
    volume<float> t1im;
    //for each subject
    vector<float> vars;
    for (int i=0; i<bvars.Nrows();i++){
      vars.push_back(bvars.element(i,j));
    }
    //keep track of number of points preceding
    int cumnum=0;
    for (int sh=0; sh<1;sh++){  // an ex-loop?  Pining for the fjords?
      m=convertToMesh(model1->getDeformedGrid(vars),model1->cells);	
      meshReg(m, flirtmats.at(j));
		
      vMeshes.push_back(m);
      int count=0;
      for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
	MeshVerts.element(3*cumnum+count,j)=(*i)->get_coord().X;
	MeshVerts.element(3*cumnum+count+1,j)=(*i)->get_coord().Y;
	MeshVerts.element(3*cumnum+count+2,j)=(*i)->get_coord().Z;
	count+=3;
      }	
      cumnum+=model1->smean.size()/3;	
    }	
  }
  avMesh = setup_avMesh(m,MeshVerts);
  return MeshVerts;
	
}

Matrix deMeanMatrix(Matrix M){
  //demean rows
  Matrix Mnew(M.Nrows(),M.Ncols());
  for (int i=0; i<M.Nrows();i++){
    float sum=0;
    for (int j=0; j<M.Ncols();j++){
      sum+=M.element(i,j); 
    }
    //sum becomes mean
    sum/=M.Ncols();
    for (int j=0; j<M.Ncols();j++){
      Mnew.element(i,j)=M.element(i,j)-sum; 
    }
  }
	
  return Mnew;
}


void do_work_bvars(){
  //**********read in bvars and models and flirt matrices***************//
  string mname;
  mname=read_bvars_ModelName(inname.value() );	
	
  //load model 
  shapeModel* model1=loadAndCreateShapeModel(mname);
  cout<<"model loaded"<<endl;
  //need number of subjects (to determine number of modes)
  vector<string> subjectnames;
  Matrix bvars;
  vector<int> vN;
  Matrix target;	//target is only used for glm
  //must include a design matrix with bvars to use glm
  //modelname is used to set a path
	vector< vector< vector<float> > > vec_fmats; 
  read_bvars(inname.value(),bvars,subjectnames, vN, pathname.value(),vec_fmats);
  target=read_ascii_matrix(designname.value());
  if (target.Nrows()!=bvars.Ncols()) {
    try {
      target=read_vest(designname.value());
    } catch (...)
    { 
      cerr << "ERROR:: Design matrix incorrect (wrong size or cannot be read)" << endl;
     exit(EXIT_FAILURE);
    }
    if (target.Nrows()!=bvars.Ncols()) {
      cerr << "ERROR:: Design matrix incorrect (wrong size or cannot be read)" << endl;
      exit(EXIT_FAILURE);
    }
  }
  //can filter meshes
  if (usePCAfilter.value()){
    //truncate number of modes to recon mesh (smoothing)
    bvars=bvars.SubMatrix(1,numModes.value(),1,bvars.Ncols());
  }
				
  volume<float> t1im;
  volume<short> segim;
				
  vector<Mesh> vMeshes;
  Matrix MeshVerts;//this is used when placing t-stats on a mesh
					
					
  //need flirt matrices
  //load their names into a vector
  ifstream flirtmats;
  vector<string> flirtmatnames;
					
  //**********done reading in bvars and models and flirt matrices ***************//
  //****************RECONSTRUCTION AND ALIGNMENT*********************/
					
  //Choose the space in which to reconstruct the meshes
  Mesh modelMeanMesh;
  ColumnVector CVmodelMeanMesh;
  Mesh avMesh;
  if (useReconNative.value()){
    flirtmats.open(flirtmatsname.value().c_str());
						
    for (unsigned int i =0; i<subjectnames.size();i++){
      string stemp;
      flirtmats>>stemp;
      flirtmatnames.push_back(stemp); ///inversion of flirtmatrix is handled in shape model function when a string is input
    }
						
    //Reconstruct in native space of the image (it recons the mni then applies flirt matrix)
    MeshVerts=recon_meshesNative( mname, bvars, CVmodelMeanMesh, modelMeanMesh,subjectnames, vec_fmats, vMeshes, avMesh);
	  cout<<"done recon"<<endl;
  }
  else if(useReconMNI.value()){
    //flirt matrix is not applied
    MeshVerts=recon_meshesMNI( model1, bvars, &CVmodelMeanMesh, &modelMeanMesh,&vMeshes, avMesh);
  }else{
    //you must choose a method of reconstruction
    cerr<<"choose a mesh reconstruction method"<<endl;
    return;
  }
  //Added in any realignment of meshes here
  if (useRigidAlign.value()){
    //use a least-squares alignment of the meshes (in this case to the mean as defined by the model)
    MeshVerts=rigid_linear_xfm(MeshVerts,CVmodelMeanMesh, modelMeanMesh, false);
  }
  avMesh=setup_avMesh(avMesh,MeshVerts);  // update for new coords (as rigid_linear_xfm may have been run)
  //****************END RECONSTRUCTION AND ALIGNMENT*********************/
  //****************DEMEAN DESIGN MATRIX AND ADD ONE COLUMNS*********************/
					
					
  cout<<"done recon and reg"<<endl;
  //if the design has not been demeaned and you are doing discriminant analysis
  //then perform demean of matrix
					
					
  bool isOne=true;
  //checks for mean Column as first column, if it decides there is a column of ones (first column) then 
  //it will not demean the design matrix
  for (int i=0;i<target.Nrows();i++){
    if (target.element(i,0)!=1){
      isOne=false;
    }
  }
						
						
  if(!isOne){
    //create demeaned design, assume a non-mean column at start
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
  }
  //this displays the demean design matrix
  cout<<"new design matrix"<<endl;
  for (int j=0;j<target.Nrows();j++){
    for (int i=0;i<target.Ncols();i++){	
      cout<<target.element(j,i)<<" ";
    }	
    cout<<endl;
  }
					
  //cout<<"done processing design"<<endl;

  //**************** NON-VERTEX ANALYSIS HERE *********************/

  if (!vertexAnalysis.value()){
						
    //calculate volumes
    vector<float> vnorm;
						
    if (useNorm.value()){
      ifstream fnorm;
      fnorm.open(normname.value().c_str());
      for (int subject=0;subject<bvars.Ncols();subject++){
	float temp;
	fnorm>>temp;
	vnorm.push_back(temp);							
      }
    }
						
    Matrix Volumes(target.Nrows(),1); //this is only used for GLM
    ofstream fvol_all;
    fvol_all.open((outname.value()+".vols").c_str());
    for (int subject=0;subject<bvars.Ncols();subject++){
      //this conditions allows you to loads volumes
      //if volumes are load the loop is broken
      read_volume(t1im,subjectnames.at(subject));
      Mesh m = vMeshes.at(subject);
							
      //fill mesh
      int bounds[6]={0,0,0,0,0,0};
      segim=make_mask_from_meshInOut(t1im,m,model1->getLabel(0),bounds);
      stringstream sstemp;
      string outnamess;
      sstemp<<subject;
      sstemp>>outnamess;
							
      float voltemp;
      volume<short> segimB;
      segimB=segim;
      voltemp=boundaryCorr(&segim, &t1im,model1->getLabel(0), thresh.value(),bounds);
      stringstream sstemp2;
      sstemp2<<model1->getLabel(0);
      string lbst;
      sstemp2>>lbst;
							
      segim.setAuxFile("MGH-Subcortical");  // for automatic colormap
      save_volume(segim,subjectnames.at(subject)+"FIRSTbcorr_lb"+lbst);
      fvol_all<<subjectnames.at(subject)<<" "<<voltemp<<" "<<voltemp*t1im.xdim()*t1im.ydim()*t1im.zdim()<<endl;
      //save volume to a matrix if you wish to use in GLM
      if (useNorm.value()){
	Volumes.element(subject,0)=voltemp*t1im.xdim()*t1im.ydim()*t1im.zdim()*vnorm.at(subject);
      }else{
	Volumes.element(subject,0)=voltemp*t1im.xdim()*t1im.ydim()*t1im.zdim();
      }
							
							
    }
    //now that volumes are load perform GLM
    cout<<"volumes loaded"<<endl;
    ColumnVector tstats;
    ColumnVector contrast(target.Ncols());
						
    for (int EV=1;EV<target.Ncols();EV++){
							
      stringstream ev2st;
      ev2st<<EV;
      string evnum;
      ev2st>>evnum;
							
      for (int i=0;i<contrast.Nrows();i++){
	if(i==EV){
	  contrast.element(i)=1;
	}else{
	  contrast.element(i)=0;
	}
	cout<<contrast.element(i)<<" ";
      }
							
      tstats=GLM_fit(target, Volumes,contrast);
      ofstream fTs;
      string ftname=outname.value()+evnum+".tstat";
      fTs.open(ftname.c_str());
      for (int i =0 ; i <tstats.Nrows();i++){
	fTs<<tstats.element(i)<<endl;
      }
    }
						
  }else{  

  //**************** VERTEX ANALYSIS HERE *********************/

    // At this point the mesh is appropriately transformed and coordinates are in MeshVerts and corresponding average mesh is avMesh
    
    // Code for the case where we want the scalar volumetric output (dot product with normal vector)
    Matrix avNormals;
    avNormals = calc_avNormals(avMesh);
    int NumSubjects=MeshVerts.Ncols();
    Matrix ScalarVals(MeshVerts.Nrows()/3,MeshVerts.Ncols());
    ScalarVals=0;
    volume4D<float> volscalar;
    int nverts=MeshVerts.Nrows()/3;
    if (!surfaceVAout.value()) {
      // do not output on the surface, instead do the new default of outputting a volume with the scalar normal dot product values (for use with randomise)
      volume<float> refim;
      if (useReconMNI.value()) { read_volume(refim,string(getenv("FSLDIR")) + "/data/standard/MNI152_T1_1mm"); } 
      else { read_volume(refim,string(getenv("FSLDIR")) + "/data/standard/MNI152_T1_1mm"); } 
      volume<float> maskvol(refim);
      maskvol=0.0f;
      volume4D<float> volnormals;
      if (debug.value()) {
	volnormals.addvolume(maskvol);
	volnormals.addvolume(maskvol);
	volnormals.addvolume(maskvol);
      }
      Matrix avCoords;
      avCoords=mean(MeshVerts,2);  // form the average coordinate location across the whole group
      for (int m=1; m<=MeshVerts.Nrows(); m+=3) { // vertex number (times 3)
	float xav=avCoords(m,1);
	float yav=avCoords(m+1,1);
	float zav=avCoords(m+2,1);
	// convert to voxel coordinates
	int xvav=MISCMATHS::round(xav/refim.xdim()), yvav=MISCMATHS::round(yav/refim.ydim()), zvav=MISCMATHS::round(zav/refim.zdim());
	for (int n=1; n<=NumSubjects; n++) {  // subject number
	  int vnum=(m+2)/3;  // careful to make this one for first vertex
	  float dotprod=0.0;
	  // convert to voxel coords (from FSL-style mm coords, which is what I assume everything is in internally - MJ)
	  float xcoord=MeshVerts(m,n);
	  float ycoord=MeshVerts(m+1,n);
	  float zcoord=MeshVerts(m+2,n);
	  // project coordinate difference onto normal
	  // NB: normals are unit length here (based on local_normal() returning this in recon_* functions)
	  dotprod += avNormals(m,1)*(xcoord-xav);
	  dotprod += avNormals(m+1,1)*(ycoord-yav);
	  dotprod += avNormals(m+2,1)*(zcoord-zav);
	  ScalarVals(vnum,n)=dotprod;
	  if (debug.value()) { 
	    volnormals(xvav,yvav,zvav,0)=avNormals(m,1);
	    volnormals(xvav,yvav,zvav,1)=avNormals(m+1,1);
	    volnormals(xvav,yvav,zvav,2)=avNormals(m+2,1);
	  }
	}
      }
      for (int n=1; n<=NumSubjects; n++) {
	volscalar.addvolume(draw_scalar_mesh(refim,avMesh,ScalarVals.Column(n)));
      }
      save_volume4D(volscalar,outname.value());
      maskvol=variancevol(volscalar);
      maskvol.binarise(1e-8);
      save_volume(maskvol,fslbasename(outname.value())+"_mask");
      if (debug.value()) { save_volume4D(volnormals,fslbasename(outname.value())+"_normals"); }
    }

    // save vertex Matrix if requested (format is all x-coords, then all y-coords, then all z-coords in each column - one subject per row)
    if (saveVertices.set()) {
      Matrix ReshapedVerts(MeshVerts.Nrows(),MeshVerts.Ncols());
      int newm=1;
      for (int m=1; m<=ReshapedVerts.Nrows(); m+=3) { // vertex number (times 3)
	for (int n=1; n<=ReshapedVerts.Ncols(); n++) {  // subject number
	  ReshapedVerts(newm,n)=MeshVerts(m,n);
	  ReshapedVerts(newm+nverts,n)=MeshVerts(m+1,n);
	  ReshapedVerts(newm+2*nverts,n)=MeshVerts(m+2,n);
	}
	newm++;
      }
      write_ascii_matrix(fslbasename(saveVertices.value())+"_mat.txt",ReshapedVerts);
      volume4D<float> vertices(ReshapedVerts.Nrows()/3,3,1,ReshapedVerts.Ncols());
      vertices.setmatrix(ReshapedVerts.t());
      save_volume4D(vertices,saveVertices.value());
    }
    
    if (!surfaceVAout.value()) {
      return;  // finish here if not outputting things to the surface (with the MVGLM tests)
    }
    
    
          // Note: "target" is the design matrix (!)
	  ColumnVector CVnorm(target.Nrows());
	  // if chosen can include a normalization EV (i.e. control for size) ...useful for vertex shape statistics
	  if (useNorm.value()){
		  ifstream normIn;
		  float normtemp;
		  normIn.open(normname.value().c_str());
		  for (int subject=0;subject<target.Nrows();subject++){
			  normIn>>normtemp;
			  CVnorm.element(subject)=normtemp;
		  }
		  target=target | CVnorm ;
	  }
	  
	  
	  
	  //****************END DEMEAN DESIGN MATRIX AND ADD ONE COLUMNS*********************/
	  
	  Matrix contrast(target.Ncols()-1,target.Ncols());
	  int EVmin=1;//ignore mean column first EV to examine
	  //int EVmax=1;//EV to include
	  int	EVmax=target.Ncols();
	  
	  for (int EV=EVmin;EV<EVmax;EV++){
		  cout<<"evs "<<EV<<endl;

		  //append the EV number to the output
		  stringstream ev2st;
		  ev2st<<EV;
		  string evnum;
		  ev2st>>evnum;
		  
		  //update average mesh
		  //these vector store tstats to plot on mesh
		  vector<float> meandifx;
		  vector<float> meandify;
		  vector<float> meandifz;
		  
		  Mesh m=vMeshes.at(0);//mesh.at(0) is used to determine topology/number of vertices
		  int count=0;
		  for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ )
		  { 
			  float meanx1=0,meany1=0,meanz1=0,meanx2=0,meany2=0, meanz2=0;
			  // for (int i=0;i<MeshVerts.Nrows();i=i+3){
			  int n1=0,n2=0;
			  for (int j=0;j<MeshVerts.Ncols();j++){
				  if (target.element(j,EVmin) <= 0 ){ //handles demeaned data
					  meanx1+=MeshVerts.element(count,j);
					  meany1+=MeshVerts.element(count+1,j);
					  meanz1+=MeshVerts.element(count+2,j);
					  n1++;
				  }else{
					  meanx2+=MeshVerts.element(count,j);
					  meany2+=MeshVerts.element(count+1,j);
					  meanz2+=MeshVerts.element(count+2,j);
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
			  
			  //as binary classification was assumed we just finish calculating the mean for each group
			  
			  //save vectors point from group1 mean to group2 mean, size 1 (they will be scaled to reflect
			  //the test statistic)
			  float vecx=meanx2-meanx1;
			  float vecy=meany2-meany1;
			  float vecz=meanz2-meanz1;
			  //	float norm=sqrt(vecx*vecx+vecy*vecy+vecz*vecz);
			  
			  //scaling occurs later
			  meandifx.push_back(vecx);///norm);
			  meandify.push_back(vecy);///norm);
			  meandifz.push_back(vecz);///norm);
			  
			  //for each point calculate mean vertex for each group
			  (*i)->_update_coord = Pt(meanx1,meany1,meanz1);
			  count+=3;
		  }
		  m.update();
		  setShapeMesh(model1,m);
		  cout<<"set shape Mesh"<<endl;
		  vector<float> scalarsT;
		  
		  //this provides the option to do stats on the vertices 
		  //create F test contrast matrix (testing each EV separately)
		  //contrast matrix is the identity matrix minus the EV'th row
		  count=0;
		  for (int j=0;j<contrast.Ncols();j++){
			  if (j!=EV){
				  for (int i=0;i<contrast.Ncols();i++){
					  if(i==j){
						  //if(i==EV){
						  contrast.element(count,i)=1;
					  }else{
						  contrast.element(count,i)=0;
					  }
					  
					  cout<<contrast.element(count,i)<<" ";
				  }
				  cout<<endl;
				  count++;
			  }
		  }
		  
		  //use multivariate test on each vertex
		  int dof1=0, dof2=0;
		  for (int i=0;i<MeshVerts.Nrows();i=i+3){
			  //use multivariate multiple regression on each vertex 
			  cout<<"do mvglm"<<endl;
			  float F=MVGLM_fit(target, MeshVerts.SubMatrix(i+1,i+3,1,MeshVerts.Ncols()).t(),contrast,dof1,dof2);
			  cout<<"done mvglm"<<endl;

			  //	tstatsx.at(i/3)*=wilkL;
			  //	tstatsy.at(i/3)*=wilkL;
			  //	tstatsz.at(i/3)*=wilkL;	
			  //this may 
			  scalarsT.push_back(F);
		  }
		  
		  fslvtkIO * fout = new fslvtkIO();
		  fout->setPoints(model1->smean);
		  fout->setPolygons(model1->cells);
		  
		  fout->setScalars(scalarsT); 
		  Matrix Mmeandif(meandifx.size(),3);
		  for (unsigned int mj=0; mj<meandifx.size(); mj++) { Mmeandif.element(mj,0)=meandifx.at(mj);  Mmeandif.element(mj,1)=meandify.at(mj); Mmeandif.element(mj,2)=meandifz.at(mj); }
		  fout->setVectors(Mmeandif);
		  
		  // also save the DOFs
		  Matrix dofvals(2,1);
		  dofvals << dof1 << dof2;
		  fout->addFieldData(dofvals,"FStatDOFs","float");
		  
		  //save the mesg with vectors
		  fout->save(outname.value()+evnum+".vtk");//,5,0);
	  }
  }
}		


void do_work_meshToVol(){
  cout<<"mesh read0"<<endl;

  volume<float> ref;
  volume<short> segim;
  read_volume(ref,inname.value());

  //	shapeModel* mesh1 = new shapeModel;
  //mesh1->setImageParameters(ref.xsize(), ref.ysize(), ref.zsize(), ref.xdim(), ref.ydim(), ref.zdim());
  //	mesh1->load_vtk(meshname.value(),1);
  //	cout<<"mesh read"<<endl;
  Mesh m1;
  m1.load(meshname.value());
  //	if (centreOrigin.value()){
  //	  m1=mesh1->getTranslatedMesh(0);
  //	}else{
  //	  m1=mesh1->getShapeMesh(0);
  //	}
  int bounds[6]={0,0,0,0,0,0};
  int label=meshLabel.value();	
	
  segim=make_mask_from_meshInOut(ref,m1,label,bounds);

  segim.setAuxFile("MGH-Subcortical");  // for automatic colormap
  save_volume(segim,outname.value());

}
void do_work_overlap(){
  //this function is working directly on volumes
  

  volume<short> gold;
  volume<short> segim;

  read_volume(gold,refname.value());
  read_volume(segim,inname.value());
  overlaps(segim,gold);
  
}

void do_work_MVGLM()
{

//	float F=MVGLM_fit(target, MeshVerts.SubMatrix(i+1,i+3,1,MeshVerts.Ncols()).t(),contrast);
	ifstream targ_in;
	targ_in.open(designname.value().c_str());
	int Nevs=numModes.value();
	int k=meshLabel.value();
	float ftemp;
	int count=0;
	while (targ_in>>ftemp)
	{
		count++;
	}
	targ_in.close();
	int N=count/Nevs ;
	
	Matrix target(N,Nevs);
	Matrix Y(N,k);
	Matrix contrast(1,Nevs);
	cout<<"There are "<<N<<" subjects"<<endl;
	targ_in.open(designname.value().c_str());
	int row=0,col=0;
	while (targ_in>>ftemp)
	{
		target.element(row,col)=ftemp;
		col++;
		if (col==Nevs)
		{
		 col=0;
			row++;
		}

	}

cout<<"target "<<target<<endl;
	ifstream y_in;
	y_in.open(inname.value().c_str());
	col=0;
	row=0;
while (y_in>>ftemp)
	{
		Y.element(row,col)=ftemp;
		col++;
		if (col==k)
		{
		 col=0;
			row++;
		}

	}
	cout<<"Y "<<Y<<endl;
	
	
	ifstream c_in;
	c_in.open(refname.value().c_str());
	col=0;
	while ( c_in>>ftemp)
	{
		contrast.element(0,col)=ftemp;
		col++;
		}
		cout<<"contrast "<<contrast<<endl;	
	
		float F=MVGLM_fit(target,Y,contrast);
	cout<<"F "<<F<<endl;
}

void do_work_reconMesh()
{
	string mname;
	mname=read_bvars_ModelName(inname.value() );	
	shapeModel* model1=loadAndCreateShapeModel(mname);
	vector<string> subjectnames;
	Matrix bvars;
	vector<int> vN;	
	Matrix flirtmat(4,4);
	vector< vector< vector<float> > > fmatv;
	read_bvars(inname.value(),bvars,subjectnames, vN, pathname.value(),fmatv);
	vector<float> vars;
	for (int i=0; i<bvars.Nrows();i++){
		vars.push_back(bvars.element(i,0));
	}
	
	
	for (int i=0; i<4;i++){
		for (int j=0; j<4;j++){
			flirtmat.element(i,j)=fmatv.at(0).at(i).at(j);
			//	cout<<flirtmat.element(i,j)<<" ";
		}
		//cout<<endl;
	}
	//flirtmat=flirtmat.i();
	for (int i=0; i<4;i++){
		for (int j=0; j<4;j++){
			fmatv.at(0).at(i).at(j)=flirtmat.element(i,j);
			//	cout<<flirtmat.element(i,j)<<" ";
		}
		//cout<<endl;
	}
	
//	model1->registerModel(fmatv.at(0));
	cout<<"vars "<<vars.size()<<endl;
	Mesh m=convertToMesh(model1->getDeformedGrid(vars),model1->cells);	
	meshReg(m,fmatv.at(0));

	m.save(outname.value(),3);
}	
	
void do_work_readBvars(const string & filename)
{
	string stemp,modelname,imagename;
	int Nsubs,Nbvars;
	ifstream fin;
	fin.open(filename.c_str());
	getline(fin,stemp);
	fin>>modelname;
	fin>>stemp>>Nsubs;
	while (Nsubs!=0)
	{
		fin>>imagename>>Nbvars;
		char blank;
		fin.read(reinterpret_cast<char*>(&blank),sizeof(char));
		float bvar;
		cout<<"BVARS"<<endl;
		cout<<imagename<<" ";
		for (int i=0;i<Nbvars;i++)
		{
			fin.read(reinterpret_cast<char*>(&bvar),sizeof(float));
			cout<<bvar<<" ";
		}
		cout<<endl;
		cout<<"TRANSFORMATION MATRIX"<<endl;
		for (int i=1;i<=16;i++)
		{
			float bvar;
			fin.read(reinterpret_cast<char*>(&bvar),sizeof(float));
			cout<<bvar<<" ";
			if (i%4 == 0)
				cout<<endl;
		}
		
		
		Nsubs--;
	}
}

void do_work_concatBvars(const string & filename, const string & fileout)
{
	
	ifstream flist;
	flist.open(filename.c_str());
	string fname,modelname,imagename;
	vector< vector<float> > all_bvars;
	vector< vector<float> > all_fmats;

	vector<string> vimages;
	int count=0;
	while (flist>>fname)
	{
		count++;
	
		string stemp,modelname_temp;
		int Nsubs,Nbvars;
		ifstream fin;
		fin.open(fname.c_str());
		getline(fin,stemp);
		fin>>modelname_temp;
	
		if (count==1)
			modelname=modelname_temp;
		if (strcmp(modelname.c_str(),modelname_temp.c_str()))
		{
			cerr<<"Bvars were generated with different models"<<endl;
			exit (EXIT_FAILURE);
		}
		fin>>stemp>>Nsubs;
		fin>>imagename>>Nbvars;
	
		vimages.push_back(imagename);
		char blank;
		fin.read(reinterpret_cast<char*>(&blank),sizeof(char));
		float bvar;
		vector<float> vbvars;
		for (int i=0;i<Nbvars;i++)
		{
			fin.read(reinterpret_cast<char*>(&bvar),sizeof(float));
			vbvars.push_back(bvar);
		}
		all_bvars.push_back(vbvars);
		vector<float> fmat;
		for (int i=1;i<=16;i++)
		{
			fin.read(reinterpret_cast<char*>(&bvar),sizeof(float));
			fmat.push_back(bvar);
			if (verbose.value()){
				cout<<bvar<<" ";
				if (i%4 == 0)
					cout<<endl;
			}
		}

		all_fmats.push_back(fmat);
		fin.close();
	}
	
	ofstream fout;
	fout.open(fileout.c_str());
	fout<<"This is a bvars file"<<endl;
	fout<<modelname<<endl;
	fout<<"NumberOfSubjects "<<count<<endl;
	vector< vector<float> >::iterator i_bvars=all_bvars.begin();
	vector< vector<float> >::iterator i_fmat=all_fmats.begin();
	for (vector<string>::iterator i_f=vimages.begin();i_f!=vimages.end(); i_f++,i_bvars++,i_fmat++)
	{
		fout<<*i_f<<" "<<i_bvars->size()<<" ";
		for (vector<float>::iterator i_bvars2=i_bvars->begin();i_bvars2!=i_bvars->end(); i_bvars2++ )
				fout.write(reinterpret_cast<char *>(&(*i_bvars2)),sizeof(*i_bvars2));

		for (vector<float>::iterator i_fmat2=i_fmat->begin(); i_fmat2!=i_fmat->end();i_fmat2++)
				fout.write(reinterpret_cast<char *>(&(*i_fmat2)),sizeof(*i_fmat2));

	}
	

}


int main(int argc,char *argv[])
{
	
  Tracer tr("main");
  OptionParser options(title, examples);
	
  try {
    // must include all wanted options here (the order determines how
    //  the help message is printed)
    options.add(inname);
    options.add(normname);
    options.add(refname);
    options.add(pathname);
    options.add(flirtmatsname);
    options.add(useScale);
    options.add(overlap);
    options.add(meshname);
    options.add(useNorm);
    options.add(surfaceVAout);
    options.add(outname);
    options.add(thresh);
    options.add(meshLabel);
    options.add(usebvars);
    options.add(useReconMNI);
    options.add(vertexAnalysis);
    options.add(useReconNative);
    options.add(useRigidAlign);
    options.add(designname);
    options.add(reconMeshFromBvars);
    options.add(readBvars);
    options.add(meshToVol);
    options.add(centreOrigin);
    options.add(saveVertices);
    options.add(verbose);
    options.add(usePCAfilter);
    options.add(numModes);
    options.add(singleBoundaryCorr);
    options.add(doMVGLM);
    options.add(concatBvars);
    options.add(debug);
    options.add(help);

    nonoptarg = options.parse_command_line(argc, argv);
		
    // line below stops the program if the help was requested or 
    //  a compulsory option was not set
    if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
      {
	options.usage();
	exit(EXIT_FAILURE);
      }
		
    // Call the local functions
    if  (usebvars.value()){
      do_work_bvars();
			
      //	}else if (imfrombvars.value()){
      //	do_work_bvarsToVols();
    }else if (meshToVol.value()){
      do_work_meshToVol();
    }else if(singleBoundaryCorr.value()){
      do_work_SingleClean();
    }else if (overlap.value()){
      do_work_overlap();
    }else if (doMVGLM.value()){
		do_work_MVGLM();
	}else if (reconMeshFromBvars.value())
	{
		do_work_reconMesh();
	}else if (readBvars.value())
	{
		do_work_readBvars(inname.value());
	}else if (concatBvars.value())
	{
		do_work_concatBvars(inname.value(), outname.value());
	}
		
		
  }  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } 
	
  return 0;// do_work(argc,argv);
}

