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
#include "csv_mesh.h"

using namespace fslsurface_name;
using namespace mesh;


//////// CSVMESH
// 1: ascii, 2:vtk, 3: gii, -1: unknown
int  meshFileType(const string& filename){
  if(filename.size()<=5){return -1;}
  string last_3 = filename.substr(filename.size()-3, 3);
  if( last_3 == ".gz" ){
    last_3 = filename.substr(filename.size()-6, 3);
    if( last_3 == "gii" ) {return CSV_GIFTI;}
  }

  if( last_3 == "gii" ) {return CSV_GIFTI;}

  ifstream f(filename.c_str());
  //reading the header
  string header;
  getline(f, header);
  {
    string::size_type pos = header.find("# vtk DataFile Version 3.0");
    if (pos!= string::npos) {
      f.close();
      return CSV_VTK;
    }
  }
  {
    string::size_type pos = header.find("#!ascii");
    if (pos != string::npos) {
      f.close();
      return CSV_ASCII;
    }
  }
  return -1;
}
bool meshExists(const string& filename){
  int type = meshFileType(filename);
  if(type>0){return true;}
  return false;
}

void CsvMesh::load(const string& filename){
  int type=meshFileType(filename);
  if(type==CSV_ASCII){
    load_ascii(filename);
  }
  else if(type==CSV_VTK){
    load_vtk(filename);
  }
  else if(type==CSV_GIFTI){
    // cerr<<" GIFTI format not yet supported"<<endl;exit(1);
   load_gifti(filename);
  }
  else{
    cerr<<"CsvMesh::load:error reading file: "<<filename<<"  ... Unknown format"<<endl;
    exit(1);
  }
  
}

void CsvMesh::load_gifti(const string& filename) {

  fslSurface<float,unsigned int> surf;
  read_surface(surf,filename);        
  
  _points.clear();
  _pvalues.clear();
  unsigned int count=0;
  for ( vector< vertex<float> >::iterator  i= surf.vbegin(); i!= surf.vend();++i,++count){
    CsvMpoint m(i->x, i->y, i->z, count);
    _points.push_back(m);
    _pvalues.push_back(0);
  }            
  count=0;
  _triangles.clear();
  for ( vector<unsigned int>::const_iterator i = surf.const_facebegin(); i != surf.const_faceend(); i+=3,++count){
    CsvMpoint m1 = get_point(*i);    
    CsvMpoint m2 = get_point(*(i+1));
    CsvMpoint m3 = get_point(*(i+2));
    
    CsvTriangle t(m1,m2,m3,count);
    //_triangles.push_back(t);    
    push_triangle(t);
  }        

  reset_pvalues();
  count = 0;
  if (surf.getNumberOfScalarData()>0) {
    for ( vector<float>::const_iterator i= surf.const_scbegin(0); i!= surf.const_scend(0);++i,++count){
      set_pvalue(count,*i);
    }
  }

}

void CsvMesh::load_vtk(const string& filename) {
  clear();
  ifstream f(filename.c_str());
  if (f.is_open())
    {	
      //reading the header
      string header;
      getline(f, header);
      string::size_type pos = header.find("# vtk DataFile Version 3.0");
      if (pos == string::npos) {
	cerr<<"CsvMesh::load_vtk:error in the header"<<endl;exit(1);
      }
      getline(f,header);
      getline(f,header);
      getline(f,header);
      int NVertices, NFaces;
      f>>header>>NVertices>>header;	  
      //reading the points
      for (int i=0; i<NVertices; i++)
	{
	  double x, y, z;
	  f>>x>>y>>z;
	  CsvMpoint m(x,y,z,i);
	  _points.push_back(m);
	  _pvalues.push_back(0);
	}
      f>>header>>NFaces>>header;
      //reading the triangles
      for (int i=0; i<NFaces; i++)
	{
	  int p0, p1, p2;
	  int j;
	  f>>j>>p0>>p1>>p2;
	  CsvTriangle t(get_point(p0), get_point(p1), get_point(p2),i);
	  push_triangle(t);
	  _tvalues.push_back(0);
	}
      f>>header>>header;
      f>>header>>header>>header;
      f>>header>>header;
      //reading the values
      for (int i=0; i<NVertices; i++)
	{
	  int val;
	  f>>val;	      
	  set_pvalue(i,val);
	}      	  
      f.close();
    }
  else {cout<<"CsvMesh::error opening file: "<<filename<<endl; exit(1);}
}

void CsvMesh::load_ascii(const string& filename) { //load a freesurfer ascii mesh
  clear();

  ifstream f(filename.c_str());
  if (f.is_open())
    {	
      //reading the header
      string header;
      getline(f, header);
      string::size_type pos = header.find("#!ascii");
      if (pos == string::npos) {
	cerr<<"CsvMesh::load_ascii:error in the header"<<endl;exit(1);
      }

      //reading the size of the mesh
      int NVertices, NFaces;
      f>>NVertices>>NFaces;
      //reading the points
      for (int i=0; i<NVertices; i++)
	{
	  double x, y, z;
	  float val;
	  f>>x>>y>>z>>val;
	  CsvMpoint m(x, y, z, i);
	  _points.push_back(m);
	  _pvalues.push_back(val);
	}      
      //reading the triangles
      for (int i=0; i<NFaces; i++)
	{
	  int p0, p1, p2;
	  float val;
	  f>>p0>>p1>>p2>>val;
	  CsvTriangle t(get_point(p0), get_point(p1), get_point(p2),i);
	  push_triangle(t);
	  _tvalues.push_back(val);
	}
      f.close();
    }
  else {cout<<"CsvMesh::load_ascii:error opening file: "<<filename<<endl; exit(1);}
 
}

void CsvMesh::save(const string& filename, const int& type)const{
  switch(type)
    {
    case CSV_ASCII:
      save_ascii(filename);break;
    case CSV_VTK:
      save_vtk(filename);break;
    default:
      save_gifti(filename,type);break;
    }
}
void CsvMesh::save_gifti(const string& s,const int& type)const{
  string filename(s);
  string last_3 = filename.substr(filename.size()-3, 3);
  if( last_3 != "gii" ){
    if(last_3 != ".gz"){
      filename=filename+".gii";
    }
  }
  
  fslSurface<float,unsigned int> surf;
  string subtype = filename.substr(filename.size()-8, 4);

  if(subtype == "surf"){
    surf.setVertices(getPointsAsVectors());
    surf.setFaces(getTrianglesAsVectors());
  }
  else if(subtype == "func"){
    surf.insertScalars(getValuesAsVectors(),0,"MyScalars");
  }
  else{
    surf.setVertices(getPointsAsVectors());
    surf.setFaces(getTrianglesAsVectors());
    surf.insertScalars(getValuesAsVectors(),0,"MyScalars");
  }

  writeGIFTI(surf,filename,type);
}
void CsvMesh::save_vtk(const string& s)const{
  string filename(s);
  string last_3 = filename.substr(filename.size()-3, 3);
  if( last_3 != "vtk" ){
    filename=filename+".vtk";
  }
  
  ofstream flot(filename.c_str());
  if(flot.is_open()){
    flot<<"# vtk DataFile Version 3.0"<<endl
	<<"surface file"<<endl
	<<"ASCII"<<endl
	<<"DATASET POLYDATA"<<endl
	<<"POINTS ";
    flot<<_points.size()<<"  float"<<endl;
    
    for (unsigned int i =0; i<_points.size();i++)  { 
      //	flot.precision(6);
      flot<<_points[i].get_coord().X<<" "
	  <<_points[i].get_coord().Y<<" "
	  <<_points[i].get_coord().Z<<endl;
#ifdef PPC64
      if ((n++ % 20) == 0) flot.flush();
#endif
    }
    flot<<"POLYGONS "<<_triangles.size()<<" "<<_triangles.size()*4<<endl;
    for ( unsigned int i=0; i<_triangles.size(); i++) 
      flot<<"3 "
	  <<_triangles[i].get_vertice(0).get_no()<<" "
	  <<_triangles[i].get_vertice(1).get_no()<<" "
	  <<_triangles[i].get_vertice(2).get_no()<<" "<<endl;
#ifdef PPC64
    if ((n++ % 20) == 0) flot.flush();
#endif
  }
  else{
    cerr<<"Csv::save_vtk:error opening file "<<filename<<" for writing"<<endl;
    exit(1);
  }


}
void CsvMesh::save_ascii(const string& s)const{
  string filename(s);
  string last_3 = filename.substr(filename.size()-3, 3);
  if( last_3 != "asc" ){filename=filename+".asc";}
  
  ofstream f(filename.c_str());
  stringstream flot;
  if (f.is_open())
    {
      int ptcount(0), tricount(0);
      for(unsigned int i=0;i<_points.size();i++){
	flot<<_points[i].get_coord().X<<" "
 	    <<_points[i].get_coord().Y<<" "
 	    <<_points[i].get_coord().Z<<" "
 	    <<_pvalues[i]<<endl; 	  
	ptcount++;
      }
      for(unsigned int i=0;i<_triangles.size();i++){
	flot<<_triangles[i].get_vertice(0).get_no()<<" "
	    <<_triangles[i].get_vertice(1).get_no()<<" "
	    <<_triangles[i].get_vertice(2).get_no()<<" "<<0<<endl;
	tricount++;
      }
      f<<"#!ascii from CsvMesh"<<endl;
      f<<ptcount<<" "<<tricount<<endl<<flot.str();
      f.close();
    }
  else cerr<<"CsvMesh::save_ascii:error opening file for writing: "<<s<<endl;

}



// ostream& operator <<(ostream& flot,const Mesh & m){
//   m.stream_mesh(flot,1);
//   return flot;
// }


void CsvMesh::push_triangle(const CsvTriangle& t){
  _triangles.push_back(t);
  for(int i=0;i<3;i++){
    _points[t.get_vertice(i).get_no()].push_triangle(t.get_no());
  }
}


const bool operator ==(const CsvMpoint &p2, const CsvMpoint &p1){
  return (fabs(p1.get_coord().X- p2.get_coord().X)<1e-8 && fabs(p1.get_coord().Y - p2.get_coord().Y)<1e-8 && fabs(p1.get_coord().Z - p2.get_coord().Z)<1e-8);
}
const bool operator ==(const CsvMpoint &p2, const Pt &p1){
  return (fabs(p1.X- p2.get_coord().X)<1e-3 && fabs(p1.Y - p2.get_coord().Y)<1e-3 && fabs(p1.Z - p2.get_coord().Z)<1e-3);
}

const Vec operator -(const CsvMpoint&p1, const CsvMpoint &p2){
  return Vec (p1.get_coord().X - p2.get_coord().X,p1.get_coord().Y - p2.get_coord().Y,p1.get_coord().Z - p2.get_coord().Z );
}

const Vec operator -(const Pt&p1, const CsvMpoint &p2){
  return Vec (p1.X - p2.get_coord().X,p1.Y - p2.get_coord().Y,p1.Z - p2.get_coord().Z );
}

const Vec operator -(const CsvMpoint&p1, const Pt &p2){
  return Vec (p1.get_coord().X - p2.X,p1.get_coord().Y - p2.Y,p1.get_coord().Z - p2.Z );
}

const Vec operator -(const ColumnVector& p1, const CsvMpoint& p2){
  return Vec (p1(1)-p2.get_coord().X,p1(2)-p2.get_coord().Y,p1(3)-p2.get_coord().Z);
}
const Vec operator -(const ColumnVector& p1, const Vec& p2){
  return Vec (p1(1)-p2.X,p1(2)-p2.Y,p1(3)-p2.Z);
}
const Vec operator -(const Vec& p1, const CsvMpoint &p2){
  return Vec(p1.X-p2.get_coord().X,p1.Y-p2.get_coord().Y,p1.Z-p2.get_coord().Z);
}
// calculate on what side of a surface a step goes to
// a step here always starts at a vertex (vertind)
// the sign corresponds to the sign of the dot-product with the 
// normal to te closest tile
int CsvMesh::step_sign(const int& vertind,const Vec& step)const{
  int trid;
  float d=0,dmin=0;
  for(int i=0;i<_points[vertind].ntriangles();i++){
    trid=_points[vertind].get_trID(i);    
    d=(_triangles[trid].normal()|step);
    if(i==0||(fabs(d)<fabs(dmin))){dmin=d;}
  }
  return (dmin>0?1:-1);
}


bool CsvTriangle::isinside(const Vec& x)const{
  Vec v0,v1,v2;
  v0=_vertice[2]-_vertice[0];
  v1=_vertice[1]-_vertice[0];
  v2=x-_vertice[0];
  double dot00=v0|v0;
  double dot01=v0|v1;
  double dot02=v0|v2;
  double dot11=v1|v1;
  double dot12=v1|v2;
  double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
  double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
  double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

  // Check if point is in triangle
  return (u > 0) && (v > 0) && (u + v < 1);

}

double CsvTriangle::dist_to_point(const Vec& x0)const{
  double d;
  Vec x1(_vertice[0].get_coord().X,_vertice[0].get_coord().Y,_vertice[0].get_coord().Z);
  Vec x2(_vertice[1].get_coord().X,_vertice[1].get_coord().Y,_vertice[1].get_coord().Z);
  Vec x3(_vertice[2].get_coord().X,_vertice[2].get_coord().Y,_vertice[2].get_coord().Z);
  Vec u;
  double dmin=1000000;
  // test edges
  u=x2-x1;
  if( ((x0-x1)|u)>0 && ((x0-x2)|u)<0 ){
    d=(((x0-x1)*(x0-x2)).norm()/(x2-x1).norm());
    if(d<dmin)dmin=d;
  }
  u=x3-x1;
  if( ((x0-x1)|u)>0 && ((x0-x3)|u)<0 ){
    d=(((x0-x1)*(x0-x3)).norm()/(x3-x1).norm());  
    if(d<dmin)dmin=d;
  }
  u=x3-x2;
  if( ((x0-x2)|u)>0 && ((x0-x3)|u)<0 ){
   d=(((x0-x2)*(x0-x3)).norm()/(x3-x2).norm());
   if(d<dmin)dmin=d;
  }
  
  d=(x0-x1).norm();if(d<dmin)dmin=d;
  d=(x0-x2).norm();if(d<dmin)dmin=d;
  d=(x0-x3).norm();if(d<dmin)dmin=d;
  return dmin;
}

  // Saad
  // algorithm from:
  // http://softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm#intersect_RayTriangle()

  const bool CsvTriangle::intersect(const vector<Pt> & p) const {
    Vec    u,v,n;   // triangle vectors
    Vec    dir,w0,w; // ray vectors
    double r, a, b;              // params to calc ray-plane intersect

    // check if point is one the vertices
    for(int ii=0;ii<=2;ii++){
      if((_vertice[ii])==p[0])return true;
      if((_vertice[ii])==p[1])return true;
    }

    // get triangle edge vectors and plane normal
    u = _vertice[1]-_vertice[0];
    v = _vertice[2]-_vertice[0];
    n = u*v;             // cross product
    if (n.norm()==0) // triangle is degenerate
      return false;                 
    

    dir = p[1]-p[0];             // ray direction vector
    w0 = p[0]-_vertice[0];
    a = -(n|w0);
    b = (n|dir);
    if (fabs(b) < 0.001) { // ray is parallel to triangle plane
      if (fabs(a) < 0.001)                 // ray lies in triangle plane
	return true;
      else return false;             // ray disjoint from plane
    }
    
    // get intersect point of ray with triangle plane
    r = a / b;
    if (r < 0.0)                   // ray goes away from triangle
      return false;                  // => no intersect

    // for a segment, also test if (r > 1.0) => no intersect
    if (r > 1.0)    
      return false;
    
    Pt I;
    I = p[0] + r * dir;           // intersect point of ray and plane
    
    // is I inside T?
    double    uu, uv, vv, wu, wv, D;
    uu = (u|u);
    uv = (u|v);
    vv = (v|v);
    w = I - _vertice[0];
    wu = (w|u);
    wv = (w|v);
    D = uv * uv - uu * vv;
    
    // get and test parametric coords
    double s, t;
    s = (uv * wv - vv * wu) / D;
    if (s < 0.0 || s > 1.0)        // I is outside T
      return false;
    t = (uv * wu - uu * wv) / D;
    if (t < 0.0 || (s + t) > 1.0)  // I is outside T
      return false;
    
    return true;                      // I is in T
    
  }
  

  const bool CsvTriangle::intersect(const vector<Pt> & p,int& ind) const {
    Tracer_Plus tr("CsvTriangle::intersect");
    Vec    u,v,n;   // triangle vectors
    Vec    dir,w0,w; // ray vectors
    double r, a, b;              // params to calc ray-plane intersect

    // check if point is one the vertices
    for(int ii=0;ii<=2;ii++){
      if((_vertice[ii])==p[0]){ind=ii;return true;}
      if((_vertice[ii])==p[1]){ind=ii;return true;}
    }

    // get triangle edge vectors and plane normal
    u = _vertice[1]-_vertice[0];
    v = _vertice[2]-_vertice[0];
    n = u*v;             // cross product
    if (n.norm()==0) // triangle is degenerate
      return false;                 
    

    dir = p[1]-p[0];             // ray direction vector
    w0 = p[0]-_vertice[0];
    a = -(n|w0);
    b = (n|dir);
    if (fabs(b) < 0.001) { // ray is parallel to triangle plane
      if (fabs(a) < 0.001)                 // ray lies in triangle plane
	{ind=0;return true;}
      else return false;             // ray disjoint from plane
    }
    
    // get intersect point of ray with triangle plane
    r = a / b;
    if (r < 0.0)                   // ray goes away from triangle
      return false;                  // => no intersect
    
    // for a segment, also test if (r > 1.0) => no intersect
    if (r > 1.0)    
      return false;

    Pt I;
    I = p[0] + r * dir;           // intersect point of ray and plane
    
    // is I inside T?
    double    uu, uv, vv, wu, wv, D;
    uu = (u|u);
    uv = (u|v);
    vv = (v|v);
    w = I - _vertice[0];
    wu = (w|u);
    wv = (w|v);
    D = uv * uv - uu * vv;
    
    // get and test parametric coords
    double s, t;
    s = (uv * wv - vv * wu) / D;
    if (s < 0.0 || s > 1.0)        // I is outside T
      return false;
    t = (uv * wu - uu * wv) / D;
    if (t < 0.0 || (s + t) > 1.0)  // I is outside T
      return false;

    // which vertex is closest to where the segment intersects?
    float x=uu-2*wu,y=vv-2*wv;
    if( x<0 ){
      if( x<y ) ind=1;
      else ind=2;
    }
    else{
      if( y<0 ) ind=2;
      else ind=0;
    }
    
    return true;                      // I is in T
    
  }



