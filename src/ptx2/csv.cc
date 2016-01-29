/*  Combined Surfaces and Volumes Class

    Saad Jbabdi  - FMRIB Image Analysis Group
 
    Copyright (C) 2010 University of Oxford  */

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

#include "csv.h"

using namespace NEWIMAGE;
using namespace MISCMATHS;
using namespace mesh;
using namespace fslvtkio;


//   tests if a step from x1 to x2 has crossed the CSV ROI
//   starts by testing if x2 is in the volume ROI
//   then tests whether the segment x1-x2 passes through
//   any of the triangles in the CSV ROI
//   it does not test all the triangles in the ROI, but only those
//   contained within the intersection set of the voxels that contain
//   x1 or x2
//   
bool CSV::has_crossed(const ColumnVector& x1,const ColumnVector& x2,
		      const vector<ColumnVector>& crossed,bool docount,bool docontinue,const float& val){
  //Tracer_Plus tr("has_crossed1");
  // x1 and x2 in voxel space
  int ix,iy,iz;
  if(nvols>0){
    ix=(int)MISCMATHS::round((float)x2(1));
    iy=(int)MISCMATHS::round((float)x2(2));
    iz=(int)MISCMATHS::round((float)x2(3));

    if(roivol(ix,iy,iz)!=0){
      for(int i=1;i<=nvols;i++){
	if(isinroi(vol2mat(ix,iy,iz),i)!=0){
	  if(docount)
	    hitvol(vol2mat(ix,iy,iz),i) += val;
	  if(!docontinue)return true;
	}
      }
      return true;
    }    
  }

  
  if(nsurfs==0){return false;}
  
  //vector<ColumnVector> crossed; // voxels crossed when going from x1 to x2
  //float line[2][3]={{x1(1),x1(2),x1(3)},{x2(1),x2(2),x2(3)}};
  //line_crossed_voxels(line,crossed);//crossed in voxels

  
  vector< pair<int,int> > localtriangles;
  for(unsigned int i=0;i<crossed.size();i++){
    int val=surfvol((int)MISCMATHS::round((float)crossed[i](1)),
		    (int)MISCMATHS::round((float)crossed[i](2)),
		    (int)MISCMATHS::round((float)crossed[i](3)))-1;
    
    if(val<0)continue;
    localtriangles.insert(localtriangles.end(),triangles[val].begin(),triangles[val].end());
  }
  
  if(localtriangles.size()>0){
    ColumnVector x1mm(4),x2mm(4);
    bool hascrossed=false;
    
    // transform into mm space
    x1mm<<x1(1)<<x1(2)<<x1(3)<<1;
    x1mm=vox2mm*x1mm;
    x2mm<<x2(1)<<x2(2)<<x2(3)<<1;
    x2mm=vox2mm*x2mm;
    
    
    // for each triangle, check for intersection and stop if any
    CsvTriangle t;
    vector<Pt> segment(2);
    Pt p1(x1mm(1),x1mm(2),x1mm(3));segment[0]=p1;
    Pt p2(x2mm(1),x2mm(2),x2mm(3));segment[1]=p2;
    
    for(unsigned int i=0;i<localtriangles.size();i++){      
      int indmesh=localtriangles[i].first;
      int indtri=localtriangles[i].second-1;

      // get the triangle
      t = roimesh[indmesh].get_triangle(indtri);
      
      if(t.intersect(segment)){//      
	// 1: update surface hit counts
	if(docount){
	  float v=roimesh[indmesh].get_tvalue(t.get_no());
	  roimesh[indmesh].set_tvalue(t.get_no(),v+val);
	}
	// 2: continue or not depending on flag 
	//    (so it is valid for stop/exclusion/inclusion/classification)
	if(docontinue) hascrossed=true;
	else
	  return true;
      }



    }
    return hascrossed;
  }
  
  return false;
}

// this one fills a vector of which rois are crossed
bool CSV::has_crossed_roi(const ColumnVector& x1,const ColumnVector& x2,
			  const vector<ColumnVector>& crossed,vector<int>& crossedrois)const{
  //Tracer_Plus tr("has_crossed2");
  // x1 and x2 in voxel space
  int ix,iy,iz;
  bool ret=false;

  // start by looking at volume roi
  if(nvols>0){

    ix=(int)MISCMATHS::round((float)x2(1));
    iy=(int)MISCMATHS::round((float)x2(2));
    iz=(int)MISCMATHS::round((float)x2(3));
    
    if(roivol(ix,iy,iz)!=0){
      for(int i=1;i<=nvols;i++){
	if(isinroi(vol2mat(ix,iy,iz),i)!=0){
	  ret=true;
	  crossedrois.push_back(volind[i-1]);	  
	}
      }
    }

  }


  if(nsurfs==0){return ret;}


  vector< pair<int,int> > localtriangles;
  for(unsigned int i=0;i<crossed.size();i++){
    //OUT(crossed[i].t());
    int val=surfvol((int)MISCMATHS::round((float)crossed[i](1)),
		    (int)MISCMATHS::round((float)crossed[i](2)),
		    (int)MISCMATHS::round((float)crossed[i](3)))-1;
    if(val<0)continue;
    localtriangles.insert(localtriangles.end(),triangles[val].begin(),triangles[val].end());
  }

  if(localtriangles.size()>0){
    ColumnVector x1mm(4),x2mm(4);
    // transform into mm space
    x1mm<<x1(1)<<x1(2)<<x1(3)<<1;
    x1mm=vox2mm*x1mm;
    x2mm<<x2(1)<<x2(2)<<x2(3)<<1;
    x2mm=vox2mm*x2mm;


    // for each triangle, check for intersection and stop if any
    CsvTriangle t;
    vector<Pt> segment;
    Pt p1(x1mm(1),x1mm(2),x1mm(3));segment.push_back(p1);
    Pt p2(x2mm(1),x2mm(2),x2mm(3));segment.push_back(p2);

    for(unsigned int i=0;i<localtriangles.size();i++){
      int indmesh=localtriangles[i].first;
      int indtri=localtriangles[i].second-1;
      t = roimesh[indmesh].get_triangle(indtri);

      if(t.intersect(segment)){
	ret=true;
	crossedrois.push_back( surfind[indmesh] );
	
	// return if one surface
	if(nsurfs==1)
	  return ret;

      }
    }
  }
  return ret;
}

// this one fills a vector of which rois are crossed, but only checking volumes (ignoring surfaces)
bool CSV::has_crossed_roi_vols(const ColumnVector& x1,
			  vector<int>& crossedrois)const{
  // x1 in voxel space
  int ix,iy,iz;
  bool ret=false;

  // start by looking at volume roi
  if(nvols>0){

    ix=(int)MISCMATHS::round((float)x1(1));
    iy=(int)MISCMATHS::round((float)x1(2));
    iz=(int)MISCMATHS::round((float)x1(3));
    
    if(roivol(ix,iy,iz)!=0){
      for(int i=1;i<=nvols;i++){
	if(isinroi(vol2mat(ix,iy,iz),i)!=0){
	  ret=true;
	  crossedrois.push_back(volind[i-1]);	  
	}
      }
    }
  }
  return ret;
}

// this one fills a vector of which rois are crossed
// and also fills in which locations are crossed (e.g. for matrix{1,3})
bool CSV::has_crossed_roi(const ColumnVector& x1,const ColumnVector& x2,
			  const vector<ColumnVector>& crossed,
			  vector<int>& crossedrois,
			  vector<int>& crossedlocs,
			  vector< pair<int,int> >& surf_Triangle,
			  bool closestvertex)const{
  //Tracer_Plus tr("CSV::hascrossed3");
  // x1 and x2 in voxel space
  int ix,iy,iz;
  bool ret=false;

  // start by looking at volume roi
  if(nvols>0){
    ix=(int)MISCMATHS::round((float)x2(1));
    iy=(int)MISCMATHS::round((float)x2(2));
    iz=(int)MISCMATHS::round((float)x2(3));
    
    if(roivol(ix,iy,iz)!=0){
      for(int i=1;i<=nvols;i++){
	if(isinroi(vol2mat(ix,iy,iz),i)!=0){
	  ret=true;
	  crossedrois.push_back(volind[i-1]);	  	  
	  crossedlocs.push_back(vol2loc(ix,iy,iz,i-1)); 
	  pair<int,int> mypair;			//only needed for surfaces,
	  mypair.first=-1;			//but vectors surf_Triangle and crossedlocs need to have same size
	  mypair.second=-1;			//it will check later if value is -1, and then it will assume that it is from a volume
	  surf_Triangle.push_back(mypair);
	}
      }
    }

  }
  if(nsurfs==0){return ret;}
  

  vector< pair<int,int> > localtriangles;
    
  for(unsigned int i=0;i<crossed.size();i++){
    //OUT(crossed[i].t());
    int val=surfvol((int)MISCMATHS::round((float)crossed[i](1)),
		    (int)MISCMATHS::round((float)crossed[i](2)),
		    (int)MISCMATHS::round((float)crossed[i](3)))-1;
    if(val<0)continue;
    localtriangles.insert(localtriangles.end(),triangles[val].begin(),triangles[val].end());
  }
  
  //  OUT(localtriangles.size());
  if(localtriangles.size()==0){return ret;}


  ColumnVector x1mm(4),x2mm(4);
  // transform into mm space
  x1mm<<x1(1)<<x1(2)<<x1(3)<<1;
  x1mm=vox2mm*x1mm;
  x2mm<<x2(1)<<x2(2)<<x2(3)<<1;
  x2mm=vox2mm*x2mm;
  
  // for each triangle, check for intersection and stop if any  
  vector<Pt> segment(2);
  Pt p1(x1mm(1),x1mm(2),x1mm(3));segment[0]=p1;
  Pt p2(x2mm(1),x2mm(2),x2mm(3));segment[1]=p2;
  
  
  int indmesh,indtri,ind;
  for(unsigned int i=0;i<localtriangles.size();i++){
    indmesh=localtriangles[i].first;
    indtri=localtriangles[i].second-1;
 
    // by default consider all active vertices as being crossed
    if(!closestvertex){
      if(roimesh[indmesh].triangle_intersect(segment,indtri)){
        ret=true;
        crossedrois.push_back( surfind[indmesh] );
      
        ind = roimesh[indmesh].get_triangle(indtri).get_vertice(0).get_no();	  
        ind = mesh2loc[indmesh][ind];

        if(ind>=0){		
	  crossedlocs.push_back(ind);
	  pair<int,int> mypair;
	  mypair.first=indmesh;
	  mypair.second=indtri;
	  surf_Triangle.push_back(mypair);	//to avoid connections between vertices of the same triangle
	}
        ind = roimesh[indmesh].get_triangle(indtri).get_vertice(1).get_no();	  
        ind = mesh2loc[indmesh][ind];

        if(ind>=0){		
	  crossedlocs.push_back(ind);
	  pair<int,int> mypair;
	  mypair.first=indmesh;
	  mypair.second=indtri;
	  surf_Triangle.push_back(mypair);	//to avoid connections between vertices of the same triangle
	}
        ind = roimesh[indmesh].get_triangle(indtri).get_vertice(2).get_no();	  
        ind = mesh2loc[indmesh][ind];

        if(ind>=0){
	  crossedlocs.push_back(ind);
	  pair<int,int> mypair;
	  mypair.first=indmesh;
	  mypair.second=indtri;
	  surf_Triangle.push_back(mypair);	//to avoid connections between vertices of the same triangle
	}
      }       
    // if closestvertex option, consider only the closest vertex crossed
    }else{
      int closest;
      if(roimesh[indmesh].triangle_intersect(segment,indtri,closest)){
        ret=true;
        crossedrois.push_back( surfind[indmesh] );
      
        ind = roimesh[indmesh].get_triangle(indtri).get_vertice(closest).get_no();	  
        ind = mesh2loc[indmesh][ind];

        if(ind>=0){
	  crossedlocs.push_back(ind);
	  pair<int,int> mypair;
	  mypair.first=indmesh;
	  mypair.second=indtri;
	  surf_Triangle.push_back(mypair);	//to avoid connections between vertices of the same triangle
	}
      }		
    }
  }
  return ret;
}

int CSV::step_sign(const int& loc,const Vec& step)const{
  if(loctype[loc]==VOLUME)return(1);
  int meshind=locsubroi[loc];
  ColumnVector stepmm(4);
  stepmm<<step.X<<step.Y<<step.Z<<1;
  stepmm=vox2mm*stepmm;
  Vec v(stepmm(1),stepmm(2),stepmm(3));
  return roimesh[meshind].step_sign(locsubloc[loc],v);
}
int CSV::coord_sign(const int& loc,const ColumnVector& x2)const{
  if(loctype[loc]==VOLUME)return(1);
  ColumnVector x1mm(4),x2mm(4),x1(4);
  x1mm<<loccoords[loc](1)<<loccoords[loc](2)<<loccoords[loc](3)<<1;
  x1=mm2vox*x1mm;
  x2mm<<x2(1)<<x2(2)<<x2(3)<<1;
  x2mm=vox2mm*x2mm;
//   OUT(x1mm.t());
//   OUT(x2mm.t());

  CsvMpoint p = roimesh[locsubroi[loc]].get_point(locsubloc[loc]);

  CsvTriangle t;
  float d=0,dmin=0;
  float dist=0,distmin=0;
  Vec xx2(x2mm(1),x2mm(2),x2mm(3));
  Vec proj;
  for(int i=0;i<p.ntriangles();i++){
    int trID=p.get_trID(i);
    t = roimesh[locsubroi[loc]].get_triangle(trID);
    d = t.normal()|(x2mm-t.get_vertice(0));
    proj=xx2-((double)d*t.normal());
    if(t.isinside(proj)){
      if(d<0)return(-1);
      else{return(1);}
    }
    //    dist=(proj-t.centroid()).norm();
    dist=t.dist_to_point(xx2);
//     OUT(dist);
//     OUT(d);
    if(i==0||(dist<=distmin)){dmin=d;distmin=dist;}
    if(dist==distmin){if(d<dmin)dmin=d;}
  }
  //OUT(dmin);
  //loc_info(loc);
  //exit(1);
   if(dmin<0)return(-1);

//   // add more triangles 
   //cout<<"MORE!!"<<endl;
  vector<ColumnVector> crossed; // voxels crossed when going from x1 to x2
  float line[2][3]={{x1(1),x1(2),x1(3)},{x2(1),x2(2),x2(3)}};
  line_crossed_voxels(line,crossed);
  vector< pair<int,int> > localtriangles;
  for(unsigned int i=0;i<crossed.size();i++){
    int val=surfvol((int)MISCMATHS::round((float)crossed[i](1)),
		    (int)MISCMATHS::round((float)crossed[i](2)),
		    (int)MISCMATHS::round((float)crossed[i](3)))-1;
    if(val<0)continue;
    localtriangles.insert(localtriangles.end(),triangles[val].begin(),triangles[val].end());
  }
  //  OUT(localtriangles.size());
  bool self;
  for(unsigned int i=0;i<localtriangles.size();i++){    
    int indmesh=localtriangles[i].first;
    int indtri=localtriangles[i].second-1;
    self=false;
    for(int j=0;j<p.ntriangles();j++)
      if(indtri==p.get_trID(j))self=true;break;//ignore self-triangles
    if(self)continue;
    t = roimesh[indmesh].get_triangle(indtri);
    if(locsubroi[loc]!=indmesh)continue; // only look at same mesh

    vector<Pt> segment;
    Pt p1(x1mm(1),x1mm(2),x1mm(3));segment.push_back(p1);
    Pt p2(x2mm(1),x2mm(2),x2mm(3));segment.push_back(p2);
    if(!t.intersect(segment))continue;//ignore if does not intersect?
    return(-1);
    //d = t.normal()|(x2mm-t.get_vertice(0));
    //if((fabs(d)<fabs(dmin)))dmin=d;
  }

  return(1); // (dmin<0?-1:1);
}


// initialise surfvol
// go through all triangles in ROI meshes
// for each triangle, determine which voxels are crossed by each of the three segments
void CSV::init_surfvol(){
  surfvol=0;
  vector<ColumnVector> crossed;

  ColumnVector x1(4),x2(4),x3(4),xx1(3),xx2(3),xx3(3);
  for(unsigned int i=0;i<roimesh.size();i++){
    int tid=0;
    for(int j=0;j<roimesh[i].ntriangles();j++){
      tid++;
      // This bit looks at intersection between triangle and voxels

      // check if triangle is active
      int n1=roimesh[i].get_triangle(j).get_vertice(0).get_no();
      int n2=roimesh[i].get_triangle(j).get_vertice(1).get_no();
      int n3=roimesh[i].get_triangle(j).get_vertice(2).get_no();

      if( roimesh[i].get_pvalue(n1) == 0 &&
	  roimesh[i].get_pvalue(n2) == 0 &&
	  roimesh[i].get_pvalue(n3) == 0 ){ continue; 
      }            
      
      // mm->vox
      x1 << roimesh[i].get_triangle(j).get_vertice(0).get_coord().X 
	 << roimesh[i].get_triangle(j).get_vertice(0).get_coord().Y 
	 << roimesh[i].get_triangle(j).get_vertice(0).get_coord().Z 
	 << 1;
      x2 << roimesh[i].get_triangle(j).get_vertice(1).get_coord().X 
	 << roimesh[i].get_triangle(j).get_vertice(1).get_coord().Y 
	 << roimesh[i].get_triangle(j).get_vertice(1).get_coord().Z 
	 << 1;
      x3 << roimesh[i].get_triangle(j).get_vertice(2).get_coord().X 
	 << roimesh[i].get_triangle(j).get_vertice(2).get_coord().Y 
	 << roimesh[i].get_triangle(j).get_vertice(2).get_coord().Z 
	 << 1;


      x1 = mm2vox*x1;
      x2 = mm2vox*x2;
      x3 = mm2vox*x3;

      xx1=x1.SubMatrix(1,3,1,1);
      xx2=x2.SubMatrix(1,3,1,1);
      xx3=x3.SubMatrix(1,3,1,1);

      float tri[3][3]={{xx1(1),xx1(2),xx1(3)},
		       {xx2(1),xx2(2),xx2(3)},
		       {xx3(1),xx3(2),xx3(3)}};

      tri_crossed_voxels(tri,crossed);
 

      update_surfvol(crossed,tid,i);


    }
  }

}
void CSV::update_surfvol(const vector<ColumnVector>& v,const int& tid,const int& meshid){
  for(unsigned int i=0;i<v.size();i++){// loop over voxels crossed by triangle tid 
    int val = surfvol((int)MISCMATHS::round((float)v[i](1)),
		      (int)MISCMATHS::round((float)v[i](2)),
		      (int)MISCMATHS::round((float)v[i](3)));
    if(val==0){//this voxel hasn't been labeled yet
      vector< pair<int,int> > t(1);
      t[0].first=meshid;      // remember which mesh this triangle is in
      t[0].second=tid;        // triangle id within that mesh
      triangles.push_back(t); // add to set of triangles that cross voxels
      surfvol((int)MISCMATHS::round((float)v[i](1)),
	      (int)MISCMATHS::round((float)v[i](2)),
	      (int)MISCMATHS::round((float)v[i](3)))=(int)triangles.size();
    }
    else{// voxel already labeled as "crossed"
      pair<int,int> mypair(meshid,tid);
      triangles[val-1].push_back(mypair); // add this triangle to the set that cross this voxel
    }
  }
}
void CSV::save_surfvol(const string& filename,const bool& binarise)const{
  if(!binarise){
    volume<float> tmpvol(surfvol.xsize(),surfvol.ysize(),surfvol.zsize());
    volume<int> cnt(surfvol.xsize(),surfvol.ysize(),surfvol.zsize());
    cnt=0;tmpvol=0;
    for(int z=0;z<surfvol.zsize();z++){
      for(int y=0;y<surfvol.ysize();y++){
	for(int x=0;x<surfvol.xsize();x++){
	  int val=surfvol(x,y,z)-1;
	  if(val<0)continue;
	  CsvTriangle t;ColumnVector vox(4),mm(4);
	  for(unsigned int i=0;i<triangles[val].size();i++){
	    t=roimesh[0].get_triangle(triangles[val][i].second-1);

	    tmpvol(x,y,z) += roimesh[0].get_pvalue(t.get_vertice(0).get_no());
	    cnt(x,y,z)    += 1;
	    tmpvol(x,y,z) += roimesh[0].get_pvalue(t.get_vertice(1).get_no());
	    cnt(x,y,z)    += 1;
	    tmpvol(x,y,z) += roimesh[0].get_pvalue(t.get_vertice(2).get_no());
	    cnt(x,y,z)    += 1;


	  }
	}
      }
    }
    for(int z=0;z<surfvol.zsize();z++){
      for(int y=0;y<surfvol.ysize();y++){
	for(int x=0;x<surfvol.xsize();x++){
	  if(cnt(x,y,z)==0){continue;}
	  tmpvol(x,y,z) /= cnt(x,y,z);
	}
      }
    }

    //copyconvert(surfvol,tmpvol);
    copybasicproperties(refvol,tmpvol);
    save_volume(tmpvol,filename);
  }
  else{
    volume<int> tmpvol(surfvol.xsize(),surfvol.ysize(),surfvol.zsize());
    copyconvert(surfvol,tmpvol);
    tmpvol=tmpvol-1;
    tmpvol.binarise(0);
    copybasicproperties(refvol,tmpvol);
    tmpvol.setDisplayMaximumMinimum(1,0);
    save_volume(tmpvol,filename);
  }
    
}
void CSV::save_normalsAsVol(const int& meshind,const string& filename)const{
  volume4D<float> normals(surfvol.xsize(),surfvol.ysize(),surfvol.zsize(),3);
  normals=0;
  ColumnVector n(3);Vec nvec;
  DiagonalMatrix d(3);Matrix v(3,3);SymmetricMatrix dd(3);
  for(int z=0;z<surfvol.zsize();z++){
    for(int y=0;y<surfvol.ysize();y++){
      for(int x=0;x<surfvol.xsize();x++){
	int val=surfvol(x,y,z)-1;
	if(val<0)continue;
	CsvTriangle t;
	for(unsigned int i=0;i<triangles[val].size();i++){
	  t=roimesh[meshind].get_triangle(triangles[val][i].second-1);
	  nvec=t.normal();
	  n<<nvec.X<<nvec.Y<<nvec.Z;
	  dd<<n*n.t();	  
	}
	dd/=float(triangles[val].size());
	EigenValues(dd,d,v);
	normals(x,y,z,0)=v(1,3);
	normals(x,y,z,1)=v(2,3);
	normals(x,y,z,2)=v(3,3);
      }
    }
  }
  copybasicproperties(refvol,normals);
  normals.setDisplayMaximumMinimum(1,-1);
  save_volume4D(normals,filename);
}

void CSV::init_hitvol(const vector<string>& fnames){
  volume<float> tmpvol;
  isinroi.ReSize(nvoxels,nvols);
  hitvol.ReSize(nvoxels,nvols);hitvol=0;
  mat2vol.ReSize(nvoxels,3);
  int indvol=0;
  for(unsigned int i=0;i<fnames.size();i++){
    if(!fsl_imageexists(fnames[i]))continue;    
    read_volume(tmpvol,fnames[i]);
    int cnt=0;indvol++;
    for(int z=0;z<tmpvol.zsize();z++)
      for(int y=0;y<tmpvol.ysize();y++)
	for(int x=0;x<tmpvol.xsize();x++){
	  if(roivol(x,y,z)==0)continue;
	  cnt++;
	  isinroi(cnt,indvol) = tmpvol(x,y,z)!=0?1:0;	  	  
	}
  }
  int cnt=0;
  for(int z=0;z<roivol.zsize();z++)
    for(int y=0;y<roivol.ysize();y++)
      for(int x=0;x<roivol.xsize();x++){
	if(roivol(x,y,z)==0)continue;
	cnt++;
	mat2vol.Row(cnt) << x << y << z;
	vol2mat(x,y,z)=cnt;
      }
  
}
void CSV::reset_values(){
  if(nvols>0){
    hitvol=0;
  }
  if(nsurfs>0){
    for(int i=0;i<nsurfs;i++){
      roimesh[i].reset_pvalues();
      roimesh[i].reset_tvalues();
    }
  }
}
void CSV::set_vol_values(const float& val){
  if(nvols>0){
    hitvol=val;
  }
}

void CSV::reset_values(const vector<int>& locs){
  for(unsigned int i=0;i<locs.size();i++)
    set_value(locs[i],0);
}
void CSV::add_value(const int& loc,const float& val){
  if(loctype[loc]==VOLUME){
    // is locind the right thing to do here?
    int x=(int)MISCMATHS::round((float)loccoords[loc](1));
    int y=(int)MISCMATHS::round((float)loccoords[loc](2));
    int z=(int)MISCMATHS::round((float)loccoords[loc](3));
    hitvol(vol2mat(x,y,z),locsubroi[loc]+1) += val;
  }
  else if(loctype[loc]==SURFACE){
    int   vertind=locsubloc[loc];
    int   meshind=locsubroi[loc];
    float curval=roimesh[meshind].get_pvalue(vertind);
    roimesh[meshind].set_pvalue(vertind,curval+val);
  }
  else{
    cerr<<"CSV::add_value:unrecognised type "<<loctype[loc]<<endl;
    exit(1);
  }  
}
void CSV::set_value(const int& loc,const float& val){
  if(loctype[loc]==VOLUME){
    // is locind the right thing to do here?
    int x=(int)MISCMATHS::round((float)loccoords[loc](1));
    int y=(int)MISCMATHS::round((float)loccoords[loc](2));
    int z=(int)MISCMATHS::round((float)loccoords[loc](3));
    hitvol(vol2mat(x,y,z),locsubroi[loc]+1) = val;
  }
  else if(loctype[loc]==SURFACE){
    int   vertind=locsubloc[loc];
    int   meshind=locsubroi[loc];
    roimesh[meshind].set_pvalue(vertind,val);
  }
  else{
    cerr<<"CSV::set_value:unrecognised type "<<loctype[loc]<<endl;
    exit(1);
  }  
}
float CSV::get_value(const int& loc)const{
  if(loctype[loc]==VOLUME){
    // is locind the right thing to do here?
    int x=(int)MISCMATHS::round((float)loccoords[loc](1));
    int y=(int)MISCMATHS::round((float)loccoords[loc](2));
    int z=(int)MISCMATHS::round((float)loccoords[loc](3));
    return(hitvol(vol2mat(x,y,z),locsubroi[loc]+1));
  }
  else if(loctype[loc]==SURFACE){
    int   vertind=locsubloc[loc];
    int   meshind=locsubroi[loc];
    return(roimesh[meshind].get_pvalue(vertind));
  }
  else{
    OUT(loc);    
    cerr<<"CSV::get_value:unrecognised type "<<loctype[loc]<<endl;
    exit(1);
  }  
}
ReturnMatrix  CSV::get_all_values()const{
  ColumnVector ret(nlocs);
  for (int i=1;i<=nlocs;i++){
    ret(i)=get_value(i-1);
  }
  ret.Release();
  return ret;
}
void CSV::set_all_values(const ColumnVector& vals){
  if(vals.Nrows()!=nlocs){
    cerr<<"CSV::set_all_values:number of locs does not match target"<<endl;
    exit(1);
  }
  for (int loc=0;loc<nlocs;loc++){
    set_value(loc,vals(loc+1));
  }
}
void CSV::set_all_values(const float& val){
  for (int loc=0;loc<nlocs;loc++){
    set_value(loc,val);
  }
}

void CSV::add_map_value(const int& loc,const float& val,const int& map){
  if(map>=(int)maps.size()){
    cerr<<"CSV::add_map_value:trying to access uninitialised map"<<endl;
    exit(1);
  }
  maps[map](loc+1) += val;  
}
void CSV::set_map_value(const int& loc,const float& val,const int& map){
  if(map>=(int)maps.size()){
    cerr<<"CSV::add_map_value:trying to access uninitialised map"<<endl;
    exit(1);
  }
  maps[map](loc+1) = val;  
}

void CSV::loc_info(const int& loc)const{
  if(loctype[loc]==VOLUME){
    cout<<"Volume location"<<endl;
    OUT(loccoords[loc].t());
    OUT(locroi[loc]);
    OUT(locsubroi[loc]);
    OUT(locsubloc[loc]);
  }
  else{
    cout<<"Surface location"<<endl;
    OUT(loccoords[loc].t());
    OUT(locroi[loc]);
    OUT(locsubroi[loc]);
    OUT(locsubloc[loc]);
    cout<<"Triangles"<<endl;
    CsvMpoint p = roimesh[locsubroi[loc]].get_point(locsubloc[loc]);
    CsvTriangle t;
    for(int i=0;i<p.ntriangles();i++){
      int trID=p.get_trID(i);
      t = roimesh[locsubroi[loc]].get_triangle(trID);
      for(int j=0;j<3;j++)
	cout<<t.get_vertice(j).get_coord().X<<" "<<t.get_vertice(j).get_coord().Y<<" "<<t.get_vertice(j).get_coord().Z<<endl;     
    }    
  }
}

void CSV::fill_volume(volume<float>& tmpvol,const int& roiind){
  for(int i=1;i<=nvoxels;i++)
    tmpvol((int)mat2vol(i,1),(int)mat2vol(i,2),(int)mat2vol(i,3)) = hitvol(i,roiind);
}
void CSV::load_volume(const string& filename){
  volume<float> tmpvol;
  read_volume(tmpvol,filename);
  ColumnVector c(3);pair<int,int> mypair;
  for(int z=0;z<tmpvol.zsize();z++)
    for(int y=0;y<tmpvol.ysize();y++)
      for(int x=0;x<tmpvol.xsize();x++){
	if(tmpvol(x,y,z)==0)continue;
	if(roivol(x,y,z)==0)nvoxels++;
	roivol(x,y,z) += 1;
	nlocs++;
	vol2loc(x,y,z,nvols) = nlocs-1;

	loctype.push_back(VOLUME);
	locroi.push_back(nrois);
	locsubroi.push_back(nvols);
	locsubloc.push_back(nvoxels);

	// update loc info
	c<<x<<y<<z;
	loccoords.push_back(c);
      }
  volind.push_back(nrois);
  roitype.push_back(VOLUME);
  roisubind.push_back(nvols);

  nvols++;
  nrois++;
}

bool lookAtAllMesh(CsvMesh& m){
  int nz=0,nnz=0;
  for(int i=0;i<m.nvertices();i++){
    if(m.get_pvalue(i)!=0){nnz++;}
    else{nz++;}
  }
  //OUT(nz);
  //OUT(nnz);
  return (nz==0 || nnz==0);
}

void CSV::load_surface(const string& filename){
  CsvMesh m;
  m.load(filename);

  // update mesh lookup table (in case only part of the surface is used
  vector<int> lu1;
  pair<int,int> mypair;
  ColumnVector c(3);int npts=-1;
  bool allMesh=lookAtAllMesh(m);

  for(int i=0;i<m.nvertices();i++){
    npts++;
    if(allMesh || m.get_pvalue(i)!=0){
      if(allMesh)
	m.set_pvalue(i,1.0);

      nlocs++;
      lu1.push_back(nlocs-1);

      loctype.push_back(SURFACE);
      locroi.push_back(nrois);
      locsubroi.push_back(nsurfs);
      locsubloc.push_back(npts);

      // update locs
      c <<m.get_point(i).get_coord().X
        <<m.get_point(i).get_coord().Y
	<<m.get_point(i).get_coord().Z;
      loccoords.push_back(c);

    }
    else{
      lu1.push_back(-1);
    }
  }
  mesh2loc.push_back(lu1);
  //OUT(nlocs);
  roimesh.push_back(m);
  surfind.push_back(nrois);
  roitype.push_back(SURFACE);
  roisubind.push_back(nsurfs);

  nsurfs++;
  nrois++;

}
void CSV::load_rois(const string& filename,bool do_change_refvol){
  vector<string> fnames;
  roinames.clear();

  // filename is volume
  if(fsl_imageexists(filename)){
    fnames.push_back(filename);
    change_refvol(fnames);
    vol2loc.reinitialize(vol2mat.xsize(),vol2mat.ysize(),vol2mat.zsize(),1);
    vol2loc=0;
    load_volume(fnames[0]);
    roinames.push_back(fnames[0]);
  }
  else if(meshExists(filename)){  
    fnames.push_back(filename);
    vol2loc.reinitialize(vol2mat.xsize(),vol2mat.ysize(),vol2mat.zsize(),1);
    vol2loc=0;
    load_surface(fnames[0]);
    roinames.push_back(fnames[0]);
  }
  else{
    // file name is ascii text file
    ifstream fs(filename.c_str());
    string tmp;
    if(fs){
      fs>>tmp;
      do{
	fnames.push_back(tmp);
	fs>>tmp;
      }while(!fs.eof());
    }
    else{
      cerr<<filename<<" does not exist"<<endl;
      exit(0);
    }
    if(do_change_refvol)
      change_refvol(fnames);
    init_vol2loc(fnames);

    for(unsigned int i=0;i<fnames.size();i++){      
      if(fsl_imageexists(fnames[i])){
	load_volume(fnames[i]);
	roinames.push_back(fnames[i]);	
      }
      else if(meshExists(fnames[i])){
	load_surface(fnames[i]);
	roinames.push_back(fnames[i]);
      }
      else{
	cerr<<"CSV::load_rois:Unknown file type: "<<fnames[i]<<endl;exit(1);
      }
    }
  }
  if(refvolset)
    init_surfvol();

  // second read pass on volumes to update lookup table
  if(refvolset)
    init_hitvol(fnames);


}
void CSV::reload_rois(const string& filename){
  cleanup();
  load_rois(filename,false);
}
void CSV::save_roi(const int& roiind,const string& fname){
  if(roitype[roiind]==VOLUME){//save volume
    volume<float> tmpvol;
    tmpvol.reinitialize(refvol.xsize(),refvol.ysize(),refvol.zsize());
    copybasicproperties(refvol,tmpvol);
    tmpvol=0;
    fill_volume(tmpvol,roisubind[roiind]+1);
    tmpvol.setDisplayMaximumMinimum(tmpvol.max(),tmpvol.min());
    save_volume(tmpvol,fname);
  }
  else{
    roimesh[roisubind[roiind]].save_gifti(fname);
    //roimesh[roisubind[roiind]].save_ascii(fname);
  }
}

void CSV::save_map(const int& roiind,const int& mapind,const string& fname){
  ColumnVector tmpvals(nlocs);
  tmpvals=get_all_values();
  set_all_values(maps[mapind]);
  
  if(roitype[roiind]==VOLUME){//save volume
    volume<float> tmpvol;
    tmpvol.reinitialize(refvol.xsize(),refvol.ysize(),refvol.zsize());
    copybasicproperties(refvol,tmpvol);
    tmpvol=0;
    fill_volume(tmpvol,roisubind[roiind]+1);
    tmpvol.setDisplayMaximumMinimum(tmpvol.max(),tmpvol.min());
    save_volume(tmpvol,fname);
  }
  else{
    roimesh[roisubind[roiind]].save_ascii(fname);
  }
  
  set_all_values(tmpvals);
}

// void CSV::save_as_volume(const string& fname){
//   volume<float> tmpvol;
//   tmpvol.reinitialize(refvol.xsize(),refvol.ysize(),refvol.zsize());
//   copybasicproperties(refvol,tmpvol);
//   tmpvol=0;
//   for(int j=1;j<nvolumes;j++){
//     for(int i=1;i<=nvoxels;i++)
//       tmpvol((int)mat2vol(i,1),(int)mat2vol(i,2),(int)mat2vol(i,3)) += hitvol(i,roiind);
//   }

//   fill_volume(tmpvol,roisubind[roiind]+1);
//   tmpvol.setDisplayMaximumMinimum(tmpvol.max(),tmpvol.min());
//   save_volume(tmpvol,fname);
  
  
// }

void CSV::save_rois(const string& fname){
  if(nrois==1)
    save_roi(0,fname);
  else{
    for(int i=0;i<nrois;i++){
      save_roi(i,fname+num2str(i+1));
    }
  }
}

void CSV::divide_rois(CSV csv2){
  bool divVOL=false;
  for(int i=0;i<nrois;i++){
    if(roitype[i]==VOLUME){//divide volumes
      if(!divVOL){
        for (int x=1; x<=hitvol.Nrows(); x++) {
          for (int y=1; y<=hitvol.Ncols(); y++) {
	    if(csv2.hitvol(x,y)){
	        hitvol(x,y)=hitvol(x,y)/csv2.hitvol(x,y);
            }else{
	      hitvol(x,y)=0;
	    }
	  }
        }
        divVOL=true;
      }
    }else{ 
      for(unsigned int j=0;j<roimesh[roisubind[i]]._pvalues.size();j++){
	if(csv2.roimesh[roisubind[i]]._pvalues[j]){
	  roimesh[roisubind[i]]._pvalues[j]=roimesh[roisubind[i]]._pvalues[j]/csv2.roimesh[roisubind[i]]._pvalues[j];
	}else{
	  roimesh[roisubind[i]]._pvalues[j]=0;
	}
      }
    }
  }
}
void CSV::divide_maps(CSV csv2){
  for(int i=0;i<maps.size();i++){
    for (int j=1; j<=maps[i].Nrows(); j++){
      if(csv2.maps[i](j)){
        maps[i](j)=maps[i](j)/csv2.maps[i](j);
      }else{
	maps[i](j)=0;
      }
    }
  }
}

void CSV::cleanup(){
  nvols=0;nsurfs=0;nvoxels=0;nlocs=0;nrois=0;
  roivol=0;
  vol2mat=0;vol2loc=0;
  surfvol=0;
  clear_all();
}
void CSV::set_convention(const string& conv){
  convention=conv;
  mm2vox.ReSize(4,4);
  vox2mm.ReSize(4,4);
  if(conv=="freesurfer"){
     mm2vox << -1/refvol.xdim() << 0               << 0                << (refvol.xsize())/2
 	   <<  0                << 0               << -1/refvol.ydim() << (refvol.zsize())/2
 	   <<  0                << 1/refvol.zdim() << 0                << (refvol.ysize())/2
 	   <<  0                << 0               << 0                << 1;
//     mm2vox << -1 << 0 << 0  << 128
// 	   <<  0 << 0 << -1 << 128
// 	   <<  0 << 1 << 0  << 128
// 	   <<  0 << 0 << 0  << 1;
    vox2mm=mm2vox.i();
    
  }
  else if(conv=="caret"){
    vox2mm = refvol.sform_mat();
    mm2vox=vox2mm.i();
  }
  else if(conv=="first"){
    vox2mm << _xdim << 0 << 0 << 0
	   << 0 << _ydim << 0 << 0
	   << 0 << 0 << _zdim << 0
	   << 0 << 0 << 0 << 1;

    Matrix Q=refvol.qform_mat();
    if(Q.Determinant()>0){//neuro
      Matrix mat(4,4);
      mat<<-1<<0<<0<<refvol.xsize()-1
	 <<0<<1<<0<<0
	 <<0<<0<<1<<0
	 <<0<<0<<0<<1;
      mat=IdentityMatrix(4); // keep this here till I figure out what to do with Neuro
      mm2vox=(mat*vox2mm).i();
    }
    else{//radio
      mm2vox=vox2mm.i();
    }
  }
  else if(conv=="vox"){
    vox2mm=IdentityMatrix(4);
    mm2vox=IdentityMatrix(4);
  }
  else{
    cerr<<"CSV::set_convention: unknown convention "<<conv<<endl;
    exit(1);
  }
}

void CSV::switch_convention(const string& new_convention,const Matrix& vox2vox,const ColumnVector& old_dims){
  Matrix old_mm2vox,old_vox2mm;
  old_vox2mm=vox2mm;
  old_mm2vox=mm2vox;
  set_convention(new_convention);

  // now transform surfaces coordinates  
  ColumnVector x4(4),x3(3);
  for(int i=0;i<nsurfs;i++){
    for(vector<CsvMpoint>::iterator it = roimesh[i]._points.begin();it!=roimesh[i]._points.end();it++){
      x4 <<(*it).get_coord().X
        <<(*it).get_coord().Y
	<<(*it).get_coord().Z
	<< 1.0;
      x3=vox_to_vox(old_mm2vox*x4,old_dims,_dims,vox2vox);
      x4<<x3(1)<<x3(2)<<x3(3)<<1.0;
      x4=vox2mm*x4;
      Pt coord(x4(1),x4(2),x4(3));
      (*it).set_coord(coord);
    }
  }
}
void CSV::switch_convention(const string& new_convention,const volume4D<float>& new2old_warp,
			    const volume<short int>& oldref,const volume<short int>& newref){
  Matrix old_mm2vox,old_vox2mm;
  old_vox2mm=vox2mm;
  old_mm2vox=mm2vox;
  set_convention(new_convention);

  // now transform surfaces coordinates  
  ColumnVector x4(4),x3(3);
  for(int i=0;i<nsurfs;i++){
    for(vector<CsvMpoint>::iterator it = roimesh[i]._points.begin();it!=roimesh[i]._points.end();it++){
      x4 <<(*it).get_coord().X
        <<(*it).get_coord().Y
	<<(*it).get_coord().Z
	<< 1.0;
      x3 = NewimageCoord2NewimageCoord(new2old_warp,false,oldref,newref,old_mm2vox*x4);

      x4<<x3(1)<<x3(2)<<x3(3)<<1.0;
      x4=vox2mm*x4;
      Pt coord(x4(1),x4(2),x4(3));
      (*it).set_coord(coord);
    }
  }
}

bool CSV::is_near_surface(const ColumnVector& pos,
			  const float& dist,
			  ColumnVector& dir,
			  float& distmesh){
  float xstep=dist/_dims(1),ystep=dist/_dims(2),zstep=dist/_dims(3);
  bool ret=false;
  float d,d2=dist*dist;
  dir.ReSize(3);CsvTriangle t;
  ColumnVector curdir(3),posmm(4),pos1(4);
  pos1<<pos(1)<<pos(2)<<pos(3)<<1;
  posmm = vox2mm*pos1;
  float mindist=-1;
  for(float x=pos(1)-xstep;x<=pos(1)+xstep;x++){
    for(float y=pos(2)-ystep;y<=pos(2)+ystep;y++){
      for(float z=pos(3)-zstep;z<=pos(3)+zstep;z++){
	d=(pos(1)-x)*(pos(1)-x)+(pos(2)-y)*(pos(2)-y)+(pos(3)-z)*(pos(3)-z);
	if(d>d2)continue;
	int ind=surfvol((int)MISCMATHS::round(x),(int)MISCMATHS::round(y),(int)MISCMATHS::round(z))-1;
	if(ind<0)continue;
	// potential surface detected
	float curdist;
	for(unsigned int i=0;i<triangles[ind].size();i++){
	  t = roimesh[ triangles[ind][i].first ].get_triangle( triangles[ind][i].second-1 );
	  Vec cog=t.centroid();
	  curdir<<cog.X-posmm(1)<<cog.Y-posmm(2)<<cog.Z-posmm(3);	  
	  curdist=curdir.SumSquare();
	  if(curdist<=d2){
	    ret=true;
	    if(mindist==-1 || curdist<=mindist){
	      mindist=curdist;
	      dir=mm2vox.SubMatrix(1,3,1,3)*curdir;
	      distmesh=sqrt(mindist);
	      if(dir.MaximumAbsoluteValue()>0)
		dir/=sqrt(dir.SumSquare());
	    }
	  }
	}
	
      }
    }
  }

  return ret;
}


// finds voxels crossed when tracing a line from P1 to P2
// P1 and P2 must be in voxels
// uses Bresenham's line algorithm
void CSV::find_crossed_voxels(const ColumnVector& P1,const ColumnVector& P2,vector<ColumnVector>& crossed){
  crossed.clear();

  float x1 = MISCMATHS::round((float)P1(1));
  float y1 = MISCMATHS::round((float)P1(2));
  float z1 = MISCMATHS::round((float)P1(3));

  float x2 = MISCMATHS::round((float)P2(1));
  float y2 = MISCMATHS::round((float)P2(2));
  float z2 = MISCMATHS::round((float)P2(3));

  float dx = (x2 - x1)*_xdim;
  float dy = (y2 - y1)*_ydim;
  float dz = (z2 - z1)*_zdim;

  float ax = fabs(dx)*2;
  float ay = fabs(dy)*2;
  float az = fabs(dz)*2;

  int sx = (int)sign(dx);
  int sy = (int)sign(dy);
  int sz = (int)sign(dz);

  int x = (int)MISCMATHS::round((float)x1);
  int y = (int)MISCMATHS::round((float)y1);
  int z = (int)MISCMATHS::round((float)z1);

  ColumnVector curVox(3);
  float xd,yd,zd;
  if(ax>=MAX(ay,az)){ // x dominant
    yd = ay - ax/2;
    zd = az - ax/2;

    while(1){
      curVox<<x<<y<<z;
      crossed.push_back(curVox);
      
      if(fabs(x-x2)<1)     // end
	break;
      
      if(yd >= 0){		// move along y
	y = y + sy;
	yd = yd - ax;
      }
      if(zd >= 0){		// move along z
	z = z + sz;
	zd = zd - ax;
      }
      x  = x  + sx;		// move along x
      yd = yd + ay;
      zd = zd + az;  
    }
  }
  else if(ay>=MAX(ax,az)){	// y dominant
    xd = ax - ay/2;
    zd = az - ay/2;

    while(1){
      curVox<<x<<y<<z;
      crossed.push_back(curVox);

      if(fabs(y-y2)<1)
	break;

      if(xd >= 0){
	x = x + sx;
	xd = xd - ay;
      }
      if(zd >= 0){
	z = z + sz;
	zd = zd - ay;
      }
      y  = y  + sy;		
      xd = xd + ax;
      zd = zd + az;
    }
  }     
  else if(az>=MAX(ax,ay)){	// z dominant
    xd = ax - az/2;
    yd = ay - az/2;

    while(1){
      curVox<<x<<y<<z;
      crossed.push_back(curVox);

      if(fabs(z-z2)<1)
	break;	    

      if(xd >= 0){
	x = x + sx;
	xd = xd - az;
      }   
      if(yd >= 0){
	y = y + sy;
	yd = yd - az;
      }
      z  = z  + sz;		
      xd = xd + ax;
      yd = yd + ay;
    
    }
  }
  return;					
}


// finds voxels crossed when tracing a line from P1 to P2
// P1 and P2 must be in voxels
// looks at ALL intersected voxel sides
void CSV::line_crossed_voxels(float line[2][3],vector<ColumnVector>& crossed)const{
  //Tracer_Plus tr("CSV::line_crossed_voxels");
  crossed.clear();
  int minx=(int)MISCMATHS::round(line[0][0]);
  int miny=(int)MISCMATHS::round(line[0][1]);
  int minz=(int)MISCMATHS::round(line[0][2]);
  int maxx=minx,maxy=miny,maxz=minz;
  int i=0;int tmpi;
  do{
    tmpi=(int)MISCMATHS::round(line[i][0]);
    minx=tmpi<minx?tmpi:minx;
    maxx=tmpi>maxx?tmpi:maxx;
    tmpi=(int)MISCMATHS::round(line[i][1]);
    miny=tmpi<miny?tmpi:miny;
    maxy=tmpi>maxy?tmpi:maxy;
    tmpi=(int)MISCMATHS::round(line[i][2]);
    minz=tmpi<minz?tmpi:minz;
    maxz=tmpi>maxz?tmpi:maxz;
    i++;
  }while(i<=1);

  ColumnVector v(3);
  if( minx==maxx && miny==maxy && minz==maxz ){
     v<<minx<<miny<<minz;
     crossed.push_back(v);
     return;
   }

  float invdir[3]={1/(line[1][0]-line[0][0]),
		   1/(line[1][1]-line[0][1]),
		   1/(line[1][2]-line[0][2])};
  float vminmax[2][3];float halfsize=.5;//_sdim/2.0;

  for(int x=minx;x<=maxx;x+=1)
    for(int y=miny;y<=maxy;y+=1)
      for(int z=minz;z<=maxz;z+=1){
	vminmax[0][0]=(float)x-halfsize;
	vminmax[1][0]=(float)x+halfsize;
	vminmax[0][1]=(float)y-halfsize;
	vminmax[1][1]=(float)y+halfsize;
	vminmax[0][2]=(float)z-halfsize;
	vminmax[1][2]=(float)z+halfsize;

	if(rayBoxIntersection(line[0],invdir,vminmax)){
	  v<<x<<y<<z;
	  crossed.push_back(v);
	}
      }

}

// this one searches within the bounding box of the triangle for voxels crossed by the whole triangle
void CSV::tri_crossed_voxels(float tri[3][3],vector<ColumnVector>& crossed){
  int minx=(int)MISCMATHS::round(tri[0][0]);
  int miny=(int)MISCMATHS::round(tri[0][1]);
  int minz=(int)MISCMATHS::round(tri[0][2]);
  int maxx=minx,maxy=miny,maxz=minz;
  crossed.clear();
  int i=0;int tmpi;
  do{
    tmpi=(int)MISCMATHS::round(tri[i][0]);
    minx=tmpi<minx?tmpi:minx;
    maxx=tmpi>maxx?tmpi:maxx;
    tmpi=(int)MISCMATHS::round(tri[i][1]);
    miny=tmpi<miny?tmpi:miny;
    maxy=tmpi>maxy?tmpi:maxy;
    tmpi=(int)MISCMATHS::round(tri[i][2]);
    minz=tmpi<minz?tmpi:minz;
    maxz=tmpi>maxz?tmpi:maxz;
    i++;
  }while(i<=2);
  //OUT(maxx-minx);OUT(maxy-miny);OUT(maxz-minz);

  float boxcentre[3],boxhalfsize[3]={.5,.5,.5};
  ColumnVector v(3);int s=1;
  for(int x=minx-s;x<=maxx+s;x+=1)
    for(int y=miny-s;y<=maxy+s;y+=1)
      for(int z=minz-s;z<=maxz+s;z+=1){
	boxcentre[0]=(float)x;
	boxcentre[1]=(float)y;
	boxcentre[2]=(float)z;
	if(triBoxOverlap(boxcentre,boxhalfsize,tri)){
	  v<<x<<y<<z;
	  crossed.push_back(v);
	}
      }
  //OUT(crossed.size());
}


////////////////////////////////////////////////////////////////////////////////////
//////   UTILS ...
////////////////////////////////////////////////////////////////////////////////////
//  by Tomas Akenine-Möller 
// 
#define X 0
#define Y 1
#define Z 2

#define CROSS(dest,v1,v2) \
  dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
  dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
  dest[2]=v1[0]*v2[1]-v1[1]*v2[0]; 

#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#define SUB(dest,v1,v2) \
  dest[0]=v1[0]-v2[0]; \
  dest[1]=v1[1]-v2[1]; \
  dest[2]=v1[2]-v2[2]; 

#define FINDMINMAX(x0,x1,x2,min,max) \
  min = max = x0;   \
  if(x1<min) min=x1;\
  if(x1>max) max=x1;\
  if(x2<min) min=x2;\
  if(x2>max) max=x2;

int planeBoxOverlap(float normal[3],float d, float maxbox[3])
{
  int q;
  float vmin[3],vmax[3];
  for(q=X;q<=Z;q++)
    {
      if(normal[q]>0.0f)
	{
	  vmin[q]=-maxbox[q];
	  vmax[q]=maxbox[q];
	}
      else
	{
	  vmin[q]=maxbox[q];
	  vmax[q]=-maxbox[q];
	}
    }
  if(DOT(normal,vmin)+d>0.0f) return 0;
  if(DOT(normal,vmax)+d>=0.0f) return 1;
  
  return 0;
}


/*======================== X-tests ========================*/
#define AXISTEST_X01(a, b, fa, fb)   \
  p0 = a*v0[Y] - b*v0[Z];          \
  p2 = a*v2[Y] - b*v2[Z];          \
  if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
  rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
  if(min>rad || max<-rad) return 0;

#define AXISTEST_X2(a, b, fa, fb)   \
  p0 = a*v0[Y] - b*v0[Z];           \
  p1 = a*v1[Y] - b*v1[Z];          \
  if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
  rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
  if(min>rad || max<-rad) return 0;

/*======================== Y-tests ========================*/
#define AXISTEST_Y02(a, b, fa, fb)   \
  p0 = -a*v0[X] + b*v0[Z];         \
  p2 = -a*v2[X] + b*v2[Z];                 \
  if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
  rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
  if(min>rad || max<-rad) return 0;

#define AXISTEST_Y1(a, b, fa, fb)   \
  p0 = -a*v0[X] + b*v0[Z];         \
  p1 = -a*v1[X] + b*v1[Z];               \
  if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
  rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
  if(min>rad || max<-rad) return 0;

/*======================== Z-tests ========================*/

#define AXISTEST_Z12(a, b, fa, fb)   \
  p1 = a*v1[X] - b*v1[Y];           \
  p2 = a*v2[X] - b*v2[Y];          \
  if(p2<p1) {min=p2; max=p1;} else {min=p1; max=p2;} \
  rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
  if(min>rad || max<-rad) return 0;

#define AXISTEST_Z0(a, b, fa, fb)   \
  p0 = a*v0[X] - b*v0[Y];   \
  p1 = a*v1[X] - b*v1[Y];           \
  if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
  rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
  if(min>rad || max<-rad) return 0;

bool triBoxOverlap(float boxcenter[3],float boxhalfsize[3],float triverts[3][3])
{

  /*    use separating axis theorem to test overlap between triangle and box */
  /*    need to test for overlap in these directions: */
  /*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
  /*       we do not even need to test these) */
  /*    2) normal of the triangle */
  /*    3) crossproduct(edge from tri, {x,y,z}-directin) */
  /*       this gives 3x3=9 more tests */
  float v0[3],v1[3],v2[3];
  //float axis[3];
  float min,max,d,p0,p1,p2,rad,fex,fey,fez;  
  float normal[3],e0[3],e1[3],e2[3];

  /* This is the fastest branch on Sun */
  /* move everything so that the boxcenter is in (0,0,0) */
  SUB(v0,triverts[0],boxcenter);
  SUB(v1,triverts[1],boxcenter);
  SUB(v2,triverts[2],boxcenter);

  /* compute triangle edges */
  SUB(e0,v1,v0);      /* tri edge 0 */
  SUB(e1,v2,v1);      /* tri edge 1 */
  SUB(e2,v0,v2);      /* tri edge 2 */

  /* Bullet 3:  */
  /*  test the 9 tests first (this was faster) */
  fex = fabs(e0[X]);
  fey = fabs(e0[Y]);
  fez = fabs(e0[Z]);
  AXISTEST_X01(e0[Z], e0[Y], fez, fey);
  AXISTEST_Y02(e0[Z], e0[X], fez, fex);
  AXISTEST_Z12(e0[Y], e0[X], fey, fex);

  fex = fabs(e1[X]);
  fey = fabs(e1[Y]);
  fez = fabs(e1[Z]);
  AXISTEST_X01(e1[Z], e1[Y], fez, fey);
  AXISTEST_Y02(e1[Z], e1[X], fez, fex);
  AXISTEST_Z0(e1[Y], e1[X], fey, fex);

  fex = fabs(e2[X]);
  fey = fabs(e2[Y]);
  fez = fabs(e2[Z]);
  AXISTEST_X2(e2[Z], e2[Y], fez, fey);
  AXISTEST_Y1(e2[Z], e2[X], fez, fex);
  AXISTEST_Z12(e2[Y], e2[X], fey, fex);

  /* Bullet 1: */
  /*  first test overlap in the {x,y,z}-directions */
  /*  find min, max of the triangle each direction, and test for overlap in */
  /*  that direction -- this is equivalent to testing a minimal AABB around */
  /*  the triangle against the AABB */

  /* test in X-direction */
  FINDMINMAX(v0[X],v1[X],v2[X],min,max);
  if(min>boxhalfsize[X] || max<-boxhalfsize[X]) return false;

  /* test in Y-direction */
  FINDMINMAX(v0[Y],v1[Y],v2[Y],min,max);
  if(min>boxhalfsize[Y] || max<-boxhalfsize[Y]) return false;

  /* test in Z-direction */
  FINDMINMAX(v0[Z],v1[Z],v2[Z],min,max);
  if(min>boxhalfsize[Z] || max<-boxhalfsize[Z]) return false;

  /* Bullet 2: */
  /*  test if the box intersects the plane of the triangle */
  /*  compute plane equation of triangle: normal*x+d=0 */
  CROSS(normal,e0,e1);
  d=-DOT(normal,v0);  /* plane eq: normal.x+d=0 */
  if(!planeBoxOverlap(normal,d,boxhalfsize)) return false;

  return true;   /* box and triangle overlaps */

}

bool rayBoxIntersection(float origin[3],float invdirection[3],float vminmax[2][3])
{
  float tmin,tmax;
  float l1   = (vminmax[0][0] - origin[0]) * invdirection[0];  
  float l2   = (vminmax[1][0] - origin[0]) * invdirection[0];    
  tmin = MIN(l1, l2);  
  tmax = MAX(l1, l2);  

  l1   = (vminmax[0][1] - origin[1]) * invdirection[1];  
  l2   = (vminmax[1][1] - origin[1]) * invdirection[1];    
  tmin = MAX(MIN(l1, l2), tmin);  
  tmax = MIN(MAX(l1, l2), tmax);  
  
  l1   = (vminmax[0][2] - origin[2]) * invdirection[2];  
  l2   = (vminmax[1][2] - origin[2]) * invdirection[2];    
  tmin = MAX(MIN(l1, l2), tmin);  
  tmax = MIN(MAX(l1, l2), tmax);  
  
  return ((tmax >= tmin) && (tmax >= 0.0f));  

}

bool segTriangleIntersection(float seg[2][3],float tri[3][3]){
  //Tracer_Plus tr("segTriangleIntersection");
  float n[3],u[3],v[3],dir[3],w0[3],w[3],I[3];
  float r,a,b,uu,uv,vv,wu,wv,D,s,t;

  SUB(u,tri[1],tri[0]);
  SUB(v,tri[2],tri[0]);
 
  CROSS(n,u,v);
  if(MAX(MAX(n[0],n[1]),n[2])==0){return false;}
  
  SUB(dir,seg[1],seg[0]);
  SUB(w0,seg[0],tri[0]);
  a=-DOT(n,w0);
  b=DOT(n,dir);
  if(fabs(b)<0.0000001){
    if(fabs(a)<0.0000001){return true;}
    else{return false;}
  }
  
  r=a/b;
  if(r<0.0){return false;}
  

  I[0]=seg[0][0]+r*dir[0];
  I[1]=seg[0][1]+r*dir[1];
  I[2]=seg[0][2]+r*dir[2];
  
  uu=DOT(u,u);
  uv=DOT(u,v);
  vv=DOT(v,v);
  SUB(w,I,tri[0]);
  wu=DOT(w,u);
  wv=DOT(w,v);
  D=uv * uv - uu * vv;
  s = (uv * wv - vv * wu) / D;
      
  if (s < 0.0 || s > 1.0)        // I is outside T
    {return false;}
  t = (uv * wu - uu * wv) / D;
  if (t < 0.0 || (s + t) > 1.0)  // I is outside T
    {return false;}
  
  return true;                      // I is in T
}


