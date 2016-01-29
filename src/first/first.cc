/*  first.cc
    Brian Patenaude
    Copyright (C) 2006-2007 University of Oxford  */

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

//#include <math.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <algorithm>
#include "math.h"

#include "utils/options.h"
#include "newimage/newimageall.h"
#include "meshclass/meshclass.h"
#include "shapeModel/shapeModel.h"
#include "fslvtkio/fslvtkio.h"

using namespace std;
using namespace NEWIMAGE;
using namespace Utilities;
using namespace mesh;
using namespace SHAPE_MODEL_NAME;
using namespace fslvtkio;



string title="first University of Oxford (Brian Patenaude)";
string examples="first -i <input image> -l <flirt matrix> -m <model>";


Option<bool> verbose(string("-v,--verbose"), false, 
				 string("switch on diagnostic messages"), 
				 false, no_argument);
Option<bool> help(string("-h,--help"), false,
			  string("display this message"),
			  false, no_argument);
			  
//Standard Inputs
Option<string> inname(string("-i,--in"), string(""),
					  string("Filename of input image to be segmented."),
					  true, requires_argument);
Option<string> outname(string("-k,--outputName"), string(""),
					   string("Output name"),
					   true, requires_argument);
Option<string> flirtmatname(string("-l,--flirtMatrix"), string(""),
						 string("Filename of flirt matrix that transform input image to MNI space (output of first_flirt)"),
						 true, requires_argument);
Option<string> modelname(string("-m,--inputModel"), string(""),
						 string("Filename of input model (the structure to be segmented)."),
						 true, requires_argument);
Option<string> modelname2(string("-p,--inputModel2"), string(""),
						 string("Filename of second input model (the structure to be segmented)."),
						 false, requires_argument);			
Option<string> bmapname(string("-b,--bmapname"), string(""),
						 string("Filename of conditional mapping matrix"),
						 false, requires_argument);
Option<string> bvarsname(string("-o,--bvars"), string(""),
						  string("Initialize using bvars from a previous segmenation. When using with --shcond specifies the \
						  the shape of the structure we are conditioning on."),
						  false, requires_argument);

Option<int> nmodes(string("-n,--nmodes"), 5,
					 string("Robust maximum. Used for global inetensity normalization."),
					 false, requires_argument);

Option<bool> intref(string("--intref,--intref"), false,
			  string("use structure specified by modelname2 as intensity reference"),
			  false, no_argument);
		
Option<bool> multiImageInput(string("--multiImageInput,--multiImageInput"), false,
			  string("use structure specified by modelname2 as intensity reference"),
			  false, no_argument);

Option<bool> binarySurfaceOutput(string("--binarySurfaceOutput,--binarySurfaceOutput"), false,
			  string("use structure specified by modelname2 as intensity reference"),
			  false, no_argument);

Option<bool> loadbvars(string("--loadbvars"), false,
					   string("load intial parameter estimates from a previous segmentation."),
					   false, no_argument);
Option<bool> shcond(string("--shcond"), false,
					   string("Use conditional shape probability"),
					   false, no_argument);
int nonoptarg;

////////////////////////////////////////////////////////////////////////////

//#define NPTS 6166
//#define NMODES 20
#define STDTRUNC 60


	class firstException : public std::exception{
		
public:
		const char* errmesg;
		firstException(const char* msg)
		{
			errmesg=msg;
		}
		
private:
			virtual const char* what() const throw()
		{
				return errmesg;
		}
	};

//global variables
int T=1;


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
ReturnMatrix vectorOfVectorsToMatrix( const vector< vector<T> > & vec){
	Matrix out(vec.size(), vec.at(0).size());
	unsigned int row=0;

	for (typename vector< vector<T> >::const_iterator i=vec.begin() ; i!=vec.end();i++, row++)
	{
		unsigned int col=0;
		for (typename vector<T>::const_iterator j=i->begin() ; j!=i->end();j++,col++)
			out.element(row,col)=*j;
	}		
	out.Release();
	return out;
}


void write_bvars(const string & imagename, const string & modelname, const vector<float> & vars, const string & outname, const vector< vector<float> > & fmatv )
{ 
 	ofstream fout;
	string name=outname+".bvars";
	fout.open(outname.c_str());
	fout.precision(10);
	fout<<"this is a bvars file"<<endl; 
	fout<<modelname<<endl;
	fout<<"NumberOfSubjects "<<1<<endl;
	fout<<imagename<<" ";
	fout<<vars.size()<<" ";

#ifdef PPC64
    int n=0;
#endif
	for (vector<float>::const_iterator i=vars.begin();i!=vars.end();i++){
		fout.write(reinterpret_cast<const char *>(&(*i)),sizeof(*i));

//		fout<<*i<<" ";
#ifdef PPC64
        if ((n++ % 50) == 0) fout.flush();
#endif
	
		
	}
	for (vector< vector<float> >::const_iterator i=fmatv.begin();i!=fmatv.end();i++)
		for (vector<float>::const_iterator i2=i->begin();i2!=i->end();i2++)
			fout.write(reinterpret_cast<const char *>(&(*i2)),sizeof(*i2));
			//fout<<*i2<<" ";
	
	fout<<endl;
	fout.close();
}

void write_vtk(const vector<float> & verts, const vector< vector<unsigned int> >  & cells, const string & name_out, const bool & is_binary)
{
	fslvtkIO *meshout = new fslvtkIO();
	meshout->setPoints(verts);
	meshout->setPolygons(vectorOfVectorsToMatrix<unsigned int>(cells));
	meshout->setBinaryWrite(is_binary);
	meshout->save(name_out);
}

vector<float> bTransform(const vector<float> & vars,const Matrix & Bmat,const unsigned int & M)
{
  	//load vars 2 (predictor) into a columnvector (newmat)
	ColumnVector Cvars2(M);
	for (unsigned int i =0; i<M;i++){
		if (i < vars.size()){
			Cvars2.element(i)=vars.at(i);
		}else{
			Cvars2.element(i)=0;
		}
	}
	
	//transform vars vars to conditional 
	ColumnVector Cvars1(M);
	Cvars1=Bmat*Cvars2;
	vector<float> v_vars1;
	for (unsigned int i =0; i<M;i++){
		v_vars1.push_back(Cvars1.element(i));
	}
	return v_vars1;
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



string read_bvars(const string & fname,vector<float> & bvars,const int & M)
{
	string stemp;
	string modelNames;
	int N;//number of subjects
		ifstream fin;
		fin.open(fname.c_str());
		//throw away first three lines 
		getline(fin,stemp);//this is bvars file
			getline(fin,modelNames);//modelnames
				fin>>stemp>>N;
				bvars.clear();
				
				
				//transform all the bvars
				for (int i=0; i<N;i++)
				{
					fin>>stemp;//read in subject id
					int nvars;//how many vars written for the subject
						fin>>nvars;
						for (int j=0;j<M;j++){
							if (j<nvars){
								float ftemp;
								fin>>ftemp;
								bvars.push_back(ftemp);
							}else
							bvars.push_back(0);
						}
				}
				return modelNames;
}


int readBmap(const string & fname, const vector<float> & b2, Matrix & matBx2, Matrix & matBcx1,unsigned int & kpred )
{
	int M; //number of subjects
	ifstream fin;

	fin.open(fname.c_str());
	string stemp;
	getline(fin,stemp);
	fin>>stemp;
	fin>>M;
	matBx2.ReSize(M,M);
	matBcx1.ReSize(M,M);
	getline(fin,stemp); 
	getline(fin,stemp);
	double ftemp;

	for (int i=0;i<M;i++){
		for (int j=0;j<M;j++){
			fin>>ftemp;
			matBx2.element(i,j)=ftemp;
		}
	}
	
	getline(fin,stemp);
	getline(fin,stemp);

	for (int i=0;i<M;i++){
		for (int j=0;j<M;j++){
			fin>>ftemp;
			matBcx1.element(i,j)=ftemp;

		}
	}
	getline(fin,stemp);
	getline(fin,stemp);
	//float kx2;
	fin>>kpred;//=kx2;
	
	//calculate scale fcator for matBCx1
	float alp=M-1.0/M;
	float gammav=(alp)/(alp-2);
	
	//calculate binner
	float binner=0;
	for (unsigned int i=0;i<b2.size();i++)
		binner+=b2.at(i)*b2.at(i);
	
	
	float scale=sqrt((alp+binner*gammav)/(alp+kpred));
	//this part causes some practical problems
	matBcx1*=scale;

	getline(fin,stemp);
	getline(fin,stemp);
	//float ftemp;
	vector<float> v_ceigs;
	
	for (int i=0;i<M;i++){
	  fin>>ftemp;
	  v_ceigs.push_back(ftemp);
	}
	
	return M;
}



inline
void normal(const vector<float> & pts, const vector< vector<unsigned int> > & localTri, const vector< vector<unsigned int> > & cells,vector<float> & vx, vector<float> & vy, vector<float> & vz )
{
	
	vx.clear();
	vy.clear();
	vz.clear();
	//calculate the normals for each triangle then references
		
	for (vector< vector<unsigned int> >::const_iterator i=localTri.begin();i!=localTri.end();i++)
	{
		float nx=0,ny=0,nz=0;
		float vx1,vy1,vz1,vx2,vy2,vz2;
		
		
//		vector<float>::iterator nx_iter=cell_nx.begin();
//		vector<float>::iterator ny_iter=cell_ny.begin();
//		vector<float>::iterator nz_iter=cell_nz.begin();

		for (vector<unsigned int>::const_iterator j=i->begin();j!=i->end();j++)
		{
			vx1=pts.at(cells.at(*j).at(2)*3)-pts.at(cells.at(*j).at(0)*3);
			vx2=pts.at(cells.at(*j).at(1)*3)-pts.at(cells.at(*j).at(0)*3);
			vy1=pts.at(cells.at(*j).at(2)*3+1)-pts.at(cells.at(*j).at(0)*3+1);
			vy2=pts.at(cells.at(*j).at(1)*3+1)-pts.at(cells.at(*j).at(0)*3+1);
			vz1=pts.at(cells.at(*j).at(2)*3+2)-pts.at(cells.at(*j).at(0)*3+2);
			vz2=pts.at(cells.at(*j).at(1)*3+2)-pts.at(cells.at(*j).at(0)*3+2);
			
		float nxt=(vy1*vz2-vz1*vy2);
			float nyt=(vz1*vx2-vx1*vz2);
			float nzt=(vx1*vy2-vy1*vx2);
			
		//	nx+=*( nx_iter + *j );
		//	ny+=*( ny_iter + *j );
		//	nz+=*( nz_iter + *j );

			nx+=nxt;
			ny+=nyt;
			nz+=nzt;
		}
		
		float d=sqrt(nx*nx+ny*ny+nz*nz);
		vx.push_back(nx/d);
		vy.push_back(ny/d);
		vz.push_back(nz/d);
	}
}


inline
float mode(const vector<float> & v_intens)
{	
	unsigned int bins=128;
	float min=*(v_intens.begin());
	float max=*(v_intens.end()-1);
	float binwidth=( max - min ) /bins;
	vector<unsigned int> bincounts;

	//innitialize bincounts to zero
	for (unsigned int b=0;b<bins;b++)
		bincounts.push_back(0);
	
	for (vector<float>::const_iterator i=v_intens.begin(); i!=v_intens.end();i++)
	{
		float bintest=min+binwidth;
		for (vector<unsigned int>::iterator b=bincounts.begin(); b!=bincounts.end() ; b++, bintest+=binwidth)
			if ((*i)<(bintest)) //find bin in which to place
			{
				(*b)++;
				break;
			}
	}
		
		//search for max bin count
		unsigned int maxcount=0;
		unsigned int maxind=0;
		for (unsigned int b=0;b<bins;b++)
			if (bincounts.at(b)>maxcount)
			{
				maxcount=bincounts.at(b);
				maxind=b;
			}
				
				return (min+maxind*binwidth+binwidth/2.0);//return middle of max bin
}

inline
void getBounds(const vector<float> & m, int *bounds, const float & xdim, const float & ydim, const float & zdim)
{	
	float xmin=1000,xmax=-1000,ymin=1000,ymax=-1000,zmin=1000,zmax=-1000;
	for (vector<float>::const_iterator i=m.begin();i!=m.end();i+=3)
	{
		if ( *i < xmin)		xmin=*i;
		if ( *i > xmax)		xmax=*i;
		if ( *(i+1) < ymin)	ymin=*(i+1);
		if ( *(i+1) > ymax)	ymax=*(i+1);
		if ( *(i+2) < zmin)	zmin=*(i+2);
		if ( *(i+2) > zmax)	zmax=*(i+2);
	}
	
	*bounds=static_cast<int>(floor(xmin/xdim)-1);
	*(bounds+1)=static_cast<int>(ceil(xmax/xdim)+1);
	*(bounds+2)=static_cast<int>(floor(ymin/ydim)-1);
	*(bounds+3)=static_cast<int>(ceil(ymax/ydim)+1);
	*(bounds+4)=static_cast<int>(floor(zmin/zdim)-1);
	*(bounds+5)=static_cast<int>(ceil(zmax/zdim)+1);	
}

inline
void intensity_hist(const volume<float> &image, const volume<short> &  mask, const vector<float> & mesh,int label , vector<float> & vgraylevels, const int* bounds)
{
	//returns a sorted min distance for each vertex 
	vgraylevels.clear();
	float dist=10000;
	
	for (int n=bounds[4];n<=bounds[5];n++)
		for (int m=bounds[2];m<=bounds[3];m++) 
			for (int l=bounds[0]; l<=bounds[1];l++){
				if (mask.value(l,m,n)==label){
					dist=image.value(l,m,n);
					if (vgraylevels.empty())
						vgraylevels.push_back(dist);
					else if (dist>=vgraylevels.back())
						vgraylevels.push_back(dist);
					else 
						for (vector<float>::iterator Iter = vgraylevels.begin( ) ; Iter !=vgraylevels.end( ) ; Iter++ )
							if (dist<*Iter)
							{
								vgraylevels.insert(Iter,dist);
								break;
							}
				}
			}
}

void draw_segment(volume<short>& image, const float & p1x,const float & p1y,const float & p1z,  const float & p2x,const float & p2y,const float & p2z, int label)
{
	double xdim = (double) image.xdim();
	double ydim = (double) image.ydim();
	double zdim = (double) image.zdim();

	//in new version of bet2
	double mininc = min(xdim,min(ydim,zdim)) * .25;

	float nx = p1x - p2x;
	float ny = p1y - p2y;
	float nz = p1z - p2z;

	double d = sqrt(nx*nx+ny*ny+nz*nz);
	nx/=d;
	ny/=d;
	nz/=d;
	
	for (double i=0; i<=d; i+=mininc)
		image(static_cast<int>(floor((p2x+i*nx)/xdim +.5)),static_cast<int>(floor((p2y+i*ny)/ydim +.5)),static_cast<int>(floor((p2z+i*nz)/zdim +.5))) = label;
}


volume<short> draw_mesh(const volume<short> & image, const vector<float> & pts, const vector< vector<unsigned int> > & triangles, int label)
{
	volume<short> imout;
	copyconvert(image,imout);
	imout=0;
  	double xdim = (double) image.xdim();
	double ydim = (double) image.ydim();
	double zdim = (double) image.zdim();
	
	//in new version of bet2
	double mininc = min(xdim,min(ydim,zdim)) * 0.5;
	
	for (vector< vector<unsigned int> >::const_iterator i = triangles.begin(); i!=triangles.end(); i++)
	{	
		float px1= pts.at(i->at(1)*3);
		float py1= pts.at(i->at(1)*3+1);
		float pz1= pts.at(i->at(1)*3+2);
						  
		float vx=pts.at(i->at(0)*3) - px1;
		float vy=pts.at(i->at(0)*3+1) - py1;
		float vz=pts.at(i->at(0)*3+2) - pz1;
						  
		float px2=pts.at(i->at(2)*3); 
		float py2=pts.at(i->at(2)*3+1); 
		float pz2=pts.at(i->at(2)*3+2);
						  
		float d=sqrt(vx*vx+vy*vy+vz*vz);
						  
		vx/=d;
		vy/=d;
		vz/=d;
						  
		for (float j=0; j<=d ;  j+=mininc)
			draw_segment(imout, px1+j*vx,py1+j*vy,pz1+j*vz,px2,py2,pz2,label+100);
						  
    }
	return imout;
}
volume<short> make_mask_from_mesh(const volume<float> & image, const vector<float> & m, const vector< vector<unsigned int> > & triangles, const int & label, const int * bounds)
{
	volume<short> mask;
	copyconvert(image,mask);
	
	mask = 0;
	mask = draw_mesh(mask, m, triangles, label);
	
	volume<short> otl=mask;
	
	vector<Pt> current;
	
	Pt c(bounds[0]-2, bounds[2]-2, bounds[4]-2);
	mask.value(static_cast<int>(c.X),static_cast<int>(c.Y),static_cast<int>(c.Z)) = label;
	current.push_back(c);

	while (!current.empty())
	{
		Pt pc = current.back();
		int x, y, z;
		x=(int) pc.X; y=(int) pc.Y; z=(int) pc.Z;
		current.pop_back();
		
		if (bounds[0]<=x-1 && mask.value(x-1, y, z)==0) 
		{
			mask.value(x-1, y, z) = label;
			current.push_back(Pt(x-1, y, z));
		}
		if (bounds[2]<=y-1 && mask.value(x, y-1, z)==0) 
		{
			mask.value(x, y-1, z) = label;
			current.push_back(Pt(x, y-1, z));
		}
		if (bounds[4]<=z-1 && mask.value(x, y, z-1)==0) 
		{
			mask.value(x, y, z-1) = label;
			current.push_back(Pt(x, y, z-1));
		}
		if (bounds[1]>=x+1 && mask.value(x+1, y, z)==0)
		{
			mask.value(x+1, y, z) = label;
			current.push_back(Pt(x+1, y, z));
		}
		if (bounds[3]>=y+1 && mask.value(x, y+1, z)==0)
		{
			mask.value(x, y+1, z) = label;
			current.push_back(Pt(x, y+1, z));
		}
		if (bounds[5]>=z+1 && mask.value(x, y, z+1)==0)
		{
			mask.value(x, y, z+1) = label;
			current.push_back(Pt(x, y, z+1)); 
		}
		
	}

	for (int i=bounds[0];i<bounds[1];i++)
		for (int j=bounds[2];j<bounds[3];j++)
			for (int k=bounds[4];k<bounds[5];k++)
				if (mask.value(i,j,k)==0)
					otl.value(i,j,k)=label;

	return otl;
}


float costfuncApp(const volume<float> & image, const shapeModel & model1, const vector<float> & vars, const bool & overide_fill)
{
	//do for a single shape
	double cost=0;//
	//  	cout<<"vars "<<endl;
	//for (int i=0;i<10; i++)
	// cout<<vars.at(i)<<" ";
	// cout<<endl;
	vector<float> dif;
	vector<float> shape = model1.getDeformedGrid(vars);
	//	cout<<" got shape"<<endl; 	
	vector<float> igrid = model1.getDeformedIGrid(vars);
	//	cout<<"got deformed grids "<<shape.size()<<" "<<igrid.size()<<endl;
	if ((!model1.getFoundMode()) && (!overide_fill))
	{
		//cout<<" fill mesh"<<endl;
		int bounds[6]={0,0,0,0,0,0};
		getBounds(shape,bounds,image.xdim(),image.ydim(),image.zdim());
		
		volume<short> mask=make_mask_from_mesh(image ,shape, model1.cells,model1.getLabel(0), bounds);
		
		vector<float> v_intens;
		intensity_hist(image,mask,shape,model1.getLabel(0),v_intens, bounds);
		
		if (v_intens.size()<1)
			throw firstException("WARNING: NO INTERIOR VOXELS TO ESTIMATE MODE");
		float mode_val=mode(v_intens);
		//cout<<"modeval "<<mode_val<<endl;
		if (model1.getMode()==mode_val)
		{
			//save_volume(image,"mask");
		  	model1.setFoundMode(true);
			//  model1.setMode(mode_val);
			
			if (verbose.value()) cout<<"found mode "<<mode_val<<endl;
		}else{
			model1.setMode(mode_val);
			cout<<"mode  "<<mode_val<<endl;
		}
	}
	//cout<<"got mode"<<endl;
	//calculate normals for vertices
	vector<float> nx;
	vector<float> ny;
	vector<float> nz;
	normal(shape,model1.localTri, model1.cells,nx,ny,nz);//could possibly optimize further
	//	cout<<"got nrmals"<<endl;
	int ipp=13;//intesnity sample per profile
	float mean=model1.getMode();
	
	const float xdim=image.xdim();
	const float ydim=image.ydim();
	const float zdim=image.zdim();
	
	vector<float>::iterator nx_i=nx.begin();
	vector<float>::iterator ny_i=ny.begin();
	vector<float>::iterator nz_i=nz.begin();
	vector<float>::iterator igrid_i=igrid.begin();
	
	float inc_x= 0.5/xdim;
	float inc_y= 0.5/ydim;
	float inc_z= 0.5/zdim;
	int count=0;
	for (vector<float>::iterator k = shape.begin(); k!= shape.end(); k+=3,nx_i++,ny_i++, nz_i++)
	{
		inc_x = (*nx_i) * 0.5/xdim;
		inc_y = (*ny_i) * 0.5/ydim;
		inc_z = (*nz_i) * 0.5/zdim;
		
		
		(*k)=(*k)/xdim - (ipp-1)*0.5 * inc_x ;
		(*(k+1))=*(k+1)/ydim - (ipp-1)*0.5 * inc_y;
		(*(k+2))= *(k+2)/zdim - (ipp-1)*0.5 * inc_z;
		
		for (int j=0;j<ipp;j++,igrid_i++, (*k)+=inc_x, (*(k+1))+=inc_y, (*(k+2))+=inc_z,count++)
		{
			dif.push_back(image.interpolate(*k,*(k+1),*(k+2))-(*igrid_i) - mean);
			//	cout<<count<<" "<<*k<<" "<<*(k+1)<<" "<<*(k+2)<<" "<<(*igrid_i)<<" "<<igrid.at(count)<<" "<<mean<<" "<<dif.back()<<endl;
		}
	}
	//	cout<<"intesnity samples"<<endl;
	//************Calculate conditional I | s *************************//
	double probIcond=0;
	vector<float> ::const_iterator ieigs_i = model1.ieigs.begin();
	for (vector< vector<float> >::const_iterator col=model1.i_precision.begin(); col != model1.i_precision.end(); col++,ieigs_i++)
	{
		float multemp=0;
		vector<float>::iterator dif_i=dif.begin();
		//		cout<<"ceigs "<<1/(*ieigs_i) - 0.5*(1/model1.Errs.at(1))<<" "<<*ieigs_i<<"  "<<model1.Errs.at(1)<<endl;

		for (vector<float>::const_iterator row=(*col).begin();row!=(*col).end(); row++, dif_i++)
			multemp+=(*dif_i)*(*row);


		probIcond+=multemp*multemp*(1/(*ieigs_i) - 0.5*(1/model1.Errs.at(1)));		
	}
	//cout<<"probicond "<<probIcond<<endl;
	//the multiplication by n-1 or n is left out becomes constant in log cost
	//probI*=M;//this multiplication is performed later
	//calculates inner product of difference between observed intensity and mean 
	//this is todo 
	float sdif=0;
	//for (vector<float>::iterator row=dif.begin();row!=dif.end();row++)
	//		sdif+=(*row)*(*row);
	for (unsigned int row=0;row<dif.size();row++)
	{
	    // cout<<dif.at(row)<<endl;
		sdif+=(dif.at(row))*(dif.at(row));
		///	cout<<"sdif "<<sdif<<endl;
	}
	//work sonly if all are tied together
	float M=model1.NumberOfSubjects;
	//sdif*=M/(model1.getShape(0).getErrs().at(1)*2);
	//the M multiplication is performed later
	sdif*=1/(model1.Errs.at(1)*2);
	probIcond+=sdif;
	
	//P(I) is exlcuded use proportionailty....doesn't work as well in practice
	
	//	cout<<"probicond2 "<<probIcond<<endl;
	/////%%%%%%%%%%%%%%%%%%%%%%%%% Now claculate costfunction
	
	//calculate cumulcative number of points across all meshes
	
	//Define constant
	double alp=M-1.0/M;
	double binner=0;
	double gammav=alp/(alp-2);
	
	float k2=static_cast<float>(model1.smean.size());
	float k1=static_cast<float>(dif.size());
	
	//This is the Mahalnobis distance of the shape
	for (vector<float>::const_iterator i=vars.begin();i!=vars.end();i++)
		binner+=(*i)*(*i);
	
	//this is the intensity portion of the cost function
	cost+= -k1/2.0*log((alp+k2)/(alp+binner*(gammav)))+(alp+k1+k2)/2*log(1+(M-1)*probIcond/(alp+binner*(gammav)));
	//cout<<"cost "<<cost<<endl;
	//this is the shaoe prior ...chooeses between no condition, one conditional, or 2 conditionals
	//add in shape prior term, it can handle 0,1, or 2 conditionals
	if (model1.getCondSet())
	{		
		//cout<<"calc conditional stuff"<<endl;
		ColumnVector mBx2map= vectorOfVectorsToMatrix(model1.getCondMat1());
		ColumnVector bx1temp(vars.size());
		for (unsigned int i=0;i<vars.size();i++)
			bx1temp.element(i)=vars.at(i)-mBx2map.element(i);
		
		//		cout<<"cond2"<<endl;
		ColumnVector Bcx1= vectorOfVectorsToMatrix<float>(model1.getCondMat2()) * bx1temp; //mBx1inv*bx1temp;
		
		float Bcinner=0;
		for (int i=0; i<Bcx1.Nrows();i++)
			Bcinner+=Bcx1.element(i)*Bcx1.element(i);
		//	cout<<"cond3"<<endl;
		cost+=(alp+k2+model1.getKPred())/2*log(1+Bcinner/(alp+model1.getKPred()));
		//and in posteriro bit
		//	cout<<"cond4"<<endl;
	}else
	{
		//p(x) prior--shape prior. No Conditional
		cost+=(alp+k2)/2*log(1+binner*gammav/alp);
	}
	//cout<<"cost "<<cost<<endl;
	return cost;
	
}



bool negGradient(const volume<float> & image, vector<float> & grad, const vector<float> & vars, const shapeModel & model1, const vector<bool> & select, const float & searchRes){
	vector<float> gradtmp;
	
	float sumsq=0;
	float costinit=0;	
	costinit=costfuncApp(image, model1, vars,false);
	
	float opGrad=0;
	//take gradient with respect to each mode
	for (int i=0; i<static_cast<int>(vars.size()); i++){
		//select vector is a bool vector that selecting the mdoes over which to take gradient
		if ( select.at(i)){
		//	cout<<"select "<<i<<endl;
			//save original mode parameters
			gradtmp=vars;	
			//increment the mode parameters
			gradtmp.at(i)=gradtmp.at(i)+searchRes;
			//order is reverse because negative gradient
			//evaluate appropriate cost 
			grad.at(i)=((costinit-costfuncApp(image, model1, gradtmp,true))/(searchRes*sqrt(model1.seigs.at(i))));//(pprev.at(i)))

			
			//this handles hard max on mode parameters, not needed in practice
			if (abs(gradtmp.at(i))>STDTRUNC){
				//ignores gradient if goes beyond truncation
				grad.at(i)=grad.at(i)*1e-11;
			}
			
			//calculate cost difference in opposite direction
			//make sure gradient not positive in both directions
			gradtmp.at(i)=gradtmp.at(i)-2*searchRes;
			
			//evaluate appropriate cost 
			opGrad=((costinit-costfuncApp(image, model1, gradtmp,true))/(-searchRes*sqrt(model1.seigs.at(i))));
			
			
			//impose rules 
			if ((opGrad>0)&&(grad.at(i)<0)){
				//if in a valley
				grad.at(i)=0;
			}else{
				//take central derivative
				grad.at(i)=(opGrad+grad.at(i))/2;
			}
			
			
			sumsq+=grad.at(i)*grad.at(i);
		}else{
			//set gradient to zero when not using that mode
			grad.at(i)=0;
		}
  }
	
	//handles cases of zero gradient
	if (sumsq==0){
		return true;
	}
	//nromalize such that size is one
	for (unsigned int i=0; i<vars.size(); i++){
		grad.at(i)/=sqrt(sumsq);	
	}
	
	return false;
}


void conjGradient(const volume<float> & image, const shapeModel & model1,vector<float> &vars, \
				  const vector<float> & relStd, const  vector<bool> & select, float & searchRes, const float & searchRmax){
	
	//define variable
	vector<float>  svec, res, resPrev;
	double gamma;
	int n=vars.size();
	double dpResPrev=0; 
	vector<float> gradtmp2=res;
	vector<float> svectmp=vars;
	vector<float> gradtmp=vars;
	bool gradzero=false;

	
	//to find cost at truncation
	vector<float> varMax=vars;
	resPrev=vars;
	res=vars;
	
	
	//the statement within test also calculate values
	if (negGradient(image,res,vars,model1,select,searchRes)){
		if (verbose.value())
		cout<<"return, no gradient"<<endl;
			return ;
	}
	resPrev=res;
	svec=svectmp=res;
	vector<float> tmpVar2;
	
	for (int i=0;i<n;i++){
		tmpVar2.push_back(0);
	}
	
	vector<bool> vmax;
	vector<float> costVMAX;///used to see if upper bound hit and that is the max
		for (int i=0;i<n;i++){
			vmax.push_back(false);
			costVMAX.push_back(10e6);
		}
		float maxStd=STDTRUNC;
		for (int iter=0;iter<40;iter++){
			
			//find line mminimizatioon and set svec to vector displavment
			
			//start fitting first mode
			
			float tmpVar=0;
			float tmpCost=0;
			vector<float> varstmp=vars;
			int j=0;
			for (int i=0;i<n;i++){
				tmpVar2.at(i)=0;
			}
			for (int i=0;i<n;i++){
				vmax.at(i)=false;
			}
			int zerocount=0; //to count number of succesive zeros.help speed search
			int searchDist=60;
			
			while ((j<30)){
				
				//enable single select
				for (unsigned int s=0; s<vars.size(); s++){
					if ((abs(varstmp.at(s)+ j*searchRes*svec.at(s))>maxStd)|(abs(gradtmp.at(s)-varstmp.at(s)-j*searchRes*svec.at(s))>relStd.at(s))){
						
						//compare against absoulte allowable maximum std and against relative varaition form start point
						
						//save the value that will render the appropriate max value
						if (!vmax.at(s)){//if not already at its max
							
							if (svec.at(s)==0){//if gradient equals zero, just leave
								tmpVar2.at(s)=0;
							}else{
								if ((abs(varstmp.at(s)+ j*searchRes*svec.at(s))>maxStd)){
									//which max did it hit, handle differently
									if ((varstmp.at(s)+ j*searchRes*svec.at(s))>maxStd){
										tmpVar2.at(s)=(maxStd-varstmp.at(s))/svec.at(s);
									}else{
										tmpVar2.at(s)=(-maxStd-varstmp.at(s))/svec.at(s);
									}
								}else{
									if (gradtmp.at(s)-varstmp.at(s)-j*searchRes*svec.at(s)<relStd.at(s)){
										tmpVar2.at(s)=(gradtmp.at(s)+relStd.at(s)-varstmp.at(s))/svec.at(s);
									}else{
										tmpVar2.at(s)=(gradtmp.at(s)-relStd.at(s)-varstmp.at(s))/svec.at(s);
									}
									
								}
							}
						}
						vmax.at(s)=true;
					}
					
					if (!vmax.at(s)){
						vars.at(s)=varstmp.at(s) + j*searchRes*svec.at(s);
					}else{
						vars.at(s)=varstmp.at(s) + tmpVar2.at(s)*svec.at(s);
					}
				}
				float cost=costfuncApp(image, model1, vars,false);
				
				if ( (j==0) || (cost<tmpCost) )
				{
					tmpVar=j*searchRes;
					tmpCost=cost;
					zerocount=0;//reset count when tmpvar is updated
				}else
					zerocount++;
			
				if (zerocount>3)//causes break in linear search if nothing found in 5
					j=searchDist;//swicth to "break"
				
				j++;
			}
			//assign values of position and displacment
			float delta=0;
			
			//this enables single select
			for (unsigned int s=0; s<svectmp.size(); s++){ 
				if ((!vmax.at(s))||((vmax.at(s))&&(tmpVar<=tmpVar2.at(s)))){
					vars.at(s)=varstmp.at(s)+tmpVar*svec.at(s);
					svectmp.at(s)=tmpVar*svectmp.at(s);
					delta+=tmpVar*svec.at(s)*tmpVar*svec.at(s);
				}else{
					vars.at(s)=varstmp.at(s)+tmpVar2.at(s)*svec.at(s);
					svectmp.at(s)=tmpVar2.at(s)*svectmp.at(s);
					delta+=tmpVar2.at(s)*svec.at(s)*tmpVar2.at(s)*svec.at(s);
				}
				
				//delta is the sum of squares of deviation from current mode parameters, if less than 0.01 leave conjugate gradient
			}
			if (delta<0.01)
				tmpVar=0;
			
			//write new parameters in neggradient call
			
			//incraese search resolution 
			if (verbose.value()){
				cout<<"vars:"<<endl;
				for (unsigned int i=0;i<vars.size();i++)
					cout<<vars.at(i)<<" ";
				cout<<endl;
			}
			
			if (tmpVar==0)
			{
				//play with this maaybe no chane in res
				searchRes=searchRes/2;
				if (searchRes<searchRmax)
					gradzero=true;
			}
			
			if(!gradzero){
				bool breakloop=false;
				
				while(!breakloop){
					negGradient(image, svectmp,vars,model1, select,searchRes);
					if (gradzero==true){
						searchRes=searchRes/2;
						if (searchRes<searchRmax){
							breakloop=true;
						}else{
							gradzero=false;
						}
					}else{
						breakloop=true;
					}
				}	
				
			}
			if (gradzero){
				costfuncApp(image,model1, vars,false);//why do i have this?
				gradzero=false;
				return ;
			}
			
			//calculate new residual dot product
			double dpRes=0;
			//for (int i=0; i<n; i++){
			for (unsigned int i=0; i<select.size(); i++){
				if (select.at(i)){
					dpResPrev+=res.at(i)*res.at(i);
					dpRes+=svectmp.at(i)*svectmp.at(i);
					dpRes+=(-svectmp.at(i)+res.at(i))*-svectmp.at(i);
				}
			}
			gamma= dpRes/dpResPrev;
			
			
			//update update vector
			
			resPrev=res;
			res=svectmp;
			
			//so that single selects may be used			
			for (unsigned int i=0; i<svec.size(); i++)
				svectmp.at(i)=svec.at(i)=res.at(i)+gamma*svectmp.at(i);
			
			}
		}


shapeModel* loadAndCreateShapeModel( const string & modelname, const int & MaxModes)
{

if (verbose.value()) cout<<"read model"<<endl;
	fslvtkIO* fmodel = new fslvtkIO(modelname,static_cast<fslvtkIO::DataType>(0));
	if (verbose.value()) cout<<"done reading model"<<endl;
	const int Npts=fmodel->getPointsAsMatrix().Nrows();
	
	
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
	unsigned int M = static_cast<unsigned int>(fmodel->getField("numSubjects").element(0,0));
	if (verbose.value()) cout<<"The model was constructed from "<<M<<" training subjects."<<endl;

	unsigned int AllModes=M;
	
	for (unsigned int i =1; i<AllModes;i++)
	{
		stringstream ss;
		ss<<i;
		string mode;
		ss>>mode;
		SmodesM=SmodesM | unwrapMatrix(fmodel->getField("mode"+mode));
		ImodesM=ImodesM | unwrapMatrix(fmodel->getField("Imode"+mode));
		
	}
	if (verbose.value()) cout<<AllModes<<" modes of variation are retained."<<endl;

	
	vector< vector<float > > Smodes = matrixToVector<float>(SmodesM);
	vector< vector<float > > Imodes = matrixToVector<float>(ImodesM);
	ImodesM.Release();
	SmodesM.Release();
	
	
	//process rest of information, including intensity variance
	vector< vector<float > > Iprec = matrixToVector<float>(fmodel->getField("iCondPrec0").t());
	vector<float > Errs =  vectorToVector<float>(fmodel->getField("ErrPriors0"));
	vector<float > se =  vectorToVector<float>(fmodel->getField("eigenValues"), AllModes);
	vector<float > ie =  vectorToVector<float>(fmodel->getField("iCondEigs0"));
	vector<float > Imean =  vectorToVector<float>(unwrapMatrix(fmodel->getField("Imean")));
	vector<int> labels =  vectorToVector<int>(fmodel->getField("labels"));



	//have read in all data and store in local structures, now delete the reader.
	delete fmodel;
	cout<<"create shapeModel "<<endl;
	//create shape model
	shapeModel* model1 = new shapeModel(Smean, Smodes, se, Imean, Imodes,Iprec, ie, M,Errs,polygons,labels);
	cout<<"done creating shapeModel "<<endl;

	return model1;

}

inline
void xfm_NEWMAT_To_Vector(const Matrix & fmatM, vector< vector<float> > & fmatv )
{
	fmatv.clear();
	//store in vector of vectors
	for (int i=0; i<4 ; i++){
		vector<float> vf;
		for (int j=0; j<4 ; j++){
			vf.push_back(fmatM.element(i,j));
		}
		fmatv.push_back(vf);
	}
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


void meshReg(Mesh* m, const string & flirtmatname){
	
	//refsize is actually target image
	int numPoints=m->nvertices();
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
	for (vector<Mpoint*>::iterator i = m->_points.begin(); i!=m->_points.end(); i++ ){
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
	for (vector<Mpoint*>::iterator i = m->_points.begin(); i!=m->_points.end(); i++ ){
		Pt newPt(NewMeshPts.element(0,count),NewMeshPts.element(1,count),NewMeshPts.element(2,count));
		(*i)->_update_coord = newPt;
		count++;
	}
	m->update();
	
	
}


int do_work(const string & inname, const string & modelname, const string & modelname2, const string & flirtmatname, const string & outname, const string & bmapname,const string & bvarsname, \
			const int & nmodes, const bool & intref, const bool & multiImageInput, const bool & shcond, const bool & loadbvars, const bool & is_binary, const float & res_in) 
{ 
	//--------------------SET CONSTANTS AND VARIABLES----------------//
	const unsigned int MaxModes=nmodes;
	float searchRes=res_in;	
	unsigned int refModes=10;
	shapeModel* model1;//main model
	shapeModel* modelRef = NULL;//used as a reference for intensity

	//shape model data, used for multiple images
	vector<float> smean1, smeanRef;
	vector< vector<float> > smodes1, smodesRef;

	vector< vector<float> > fmatv, fmatv_org; //transformation matrix
	Matrix fmatM(4,4); //transformation matrix , NEWMAT form
	Matrix fmatM_Prev(4,4); //transformation matrix from previous image , NEWMAT form
	vector<string> image_list;
	vector<string> xfm_list;
	vector<string> out_list;

	//------------------------------READ/set up models ---------------//
	model1=loadAndCreateShapeModel(modelname, MaxModes );
	
	
	if (intref)
	{
		modelRef=loadAndCreateShapeModel(modelname2, refModes );
		if ( multiImageInput )
		{
			smeanRef=modelRef->smean;
			smodesRef=modelRef->smodes;
		}
	}
	
	
	//-------------------READ IN IMAGE, XFM, AND OUTPUT NAMES, and STORE SHAPE INFO------------------//
	if ( multiImageInput)
	{
		//---------------store shape info----------------//
		
		smean1=model1->smean;
		smodes1=model1->smodes;
			
		
	//--------------read in names -----------------//
		ifstream fin;
		fin.open(inname.c_str());
		string stemp;
		while (fin>>stemp)
		{
			image_list.push_back(stemp);
			fin>>stemp;
			xfm_list.push_back(stemp);
			fin>>stemp;
			out_list.push_back(stemp);
		 if (verbose.value()) cout<<"Input list "<<image_list.back()<<" "<<xfm_list.back()<<" "<<out_list.back()<<endl;
		}
	}else
	{
		image_list.push_back(inname);
		xfm_list.push_back(flirtmatname);
		out_list.push_back(outname);
		
	}
	
	//-------------------SET PREVIOUS MATRIX TO IDENTITY--------------------//
	for (unsigned int i=0;i<4;i++)
		for (unsigned int j=0;j<4;j++)
			if (i==j)
				fmatM_Prev.element(i,j)=1;
			else	
				fmatM_Prev.element(i,j)=0;
	
	//--------------------LOOP OVER ALL IMAGES--------------------//
	vector<string>::iterator xfm_iter=xfm_list.begin();
	vector<string>::iterator out_iter=out_list.begin();
	for (vector<string>::iterator im_iter=image_list.begin(); im_iter!=image_list.end();im_iter++, xfm_iter++, out_iter++)
	{
		
		 //--------------IF MULTI IMAGE RESET SHAPE INFO---------------//
		if (multiImageInput)
		{
			searchRes=res_in;	
			model1->smean=smean1;
			model1->smodes=smodes1;
			model1->setFoundMode(false);
			model1->setMode(0);
		}
		//------------------READ IN IMAGE AND NORMALIZE----------------------//
		//load base volume
		if (verbose.value()) cout<<"reading image "<<*im_iter<<endl;

		volume<float> image;
		read_volume(image,*im_iter);
		
		//normalize image intensities
		if (verbose.value()) cout<<"normalize intensity..."<<endl;
	
		image=(image-image.robustmin())*255/(image.robustmax()-image.robustmin());
		
		//------------------READ IN AND APPLY TRANSFORMATION MATRIX-----------------//
		if (verbose.value()) cout<<"Reading transformation matrix "<<*xfm_iter<<endl;
		
		ifstream ifmat;
		ifmat.open(xfm_iter->c_str());
		for (int i=0; i<4 ; i++)
		{
			for (int j=0; j<4 ; j++)
			{
				float ftemp;
				ifmat>>ftemp;
				if (verbose.value()) cout<<ftemp<<" ";
				fmatM.element(i,j)=ftemp;
			}
			if (verbose.value()) cout<<endl;
		}
		
		//invert matrix and assigned previous for purpose of looping
		xfm_NEWMAT_To_Vector(fmatM,fmatv_org);
	fmatM=fmatM.i();
//concert to vector form
		xfm_NEWMAT_To_Vector(fmatM,fmatv);
		
		//------------------FIT REFERENCE MODEL IF NECESSARY-----------------//
		
		if (intref)
		{
			if (multiImageInput)
		{
			searchRes=res_in;	
			modelRef->smean=smeanRef;
			modelRef->smodes=smodesRef;
			modelRef->setFoundMode(false);
			modelRef->setMode(0);
		}
		
			modelRef->registerModel(fmatv);
			cout<<"done registering model"<<endl;
			vector<float> varsRef;
			for (unsigned int i=0; i < refModes;i++)
				varsRef.push_back(0);
			vector<bool> selectRef;
			vector<float> relStdRef;
			for (unsigned int i=0;i<varsRef.size();i++)
			{
				selectRef.push_back(true);
				relStdRef.push_back(STDTRUNC);
				
			}
			conjGradient(image, *modelRef,varsRef, relStdRef, selectRef, searchRes, 0.15);
			//Has fit reference model
			vector<float> shape = modelRef->getDeformedGrid(varsRef);
			
			int bounds[6]={0,0,0,0,0,0};
			getBounds(shape,bounds,image.xdim(),image.ydim(),image.zdim());
			
			volume<short> mask=make_mask_from_mesh(image , shape, modelRef->cells,modelRef->getLabel(0), bounds);
			vector<float> v_intens;
			intensity_hist(image,mask,shape,modelRef->getLabel(0),v_intens, bounds);
			
			if (v_intens.size()<1)
				throw firstException("WARNING: NO INTERIOR VOXELS TO ESTIMATE MODE");
			float mode_val=mode(v_intens);
			
			
			if (verbose.value()) cout<<"found reference mode "<<mode_val<<endl;
			
			model1->setFoundMode(true);
			model1->setMode(mode_val);
		}
		
			//------------------------SET UP NECESSARY VECTORS ----------------------//

		vector<float> vars;
		for (unsigned int i=0; i<MaxModes;i++)
			vars.push_back(0);
		
		vector<bool> select;
		vector<float> relStd;
		for (unsigned int i=0;i<vars.size();i++)
		{
			select.push_back(true);
			relStd.push_back(STDTRUNC);
			
		}
			//-----------------------------LOAD PREVIOUS SEGMENTATION-----------------------//			
		  if (loadbvars)
		  {
			if (verbose.value()) cout<<"Loads previous mode parameters."<<endl;
		     read_bvars(bvarsname,vars,model1->smodes.size());
		  }
				//-------------------------SET UP CONDITIONAL IF NEEDED ----------------------//

		  if (shcond)
		  {
			  	cout<<"shcond"<<endl;
			  int M;
			  Matrix mBx2;
			  Matrix mBcx1;
			  unsigned int kpred;
			  M=readBmap(bmapname,vars,mBx2,mBcx1,kpred);
			  ColumnVector Bx2(M);
			  Bx2=0;
			  //load bx2 vars that were read
			  for (unsigned int i=0;i<vars.size();i++)
				  Bx2.element(i)=vars.at(i);
		
			  Matrix mBx2map=mBx2*Bx2;
		
			  //bmap now directly reflects the transformation t0oo bc
			  //mBx1inv=mBcx1;

			  Matrix mBx1inv=mBcx1.i();
					  
			  model1->setCondMats(matrixToVector<float>(mBx2map), matrixToVector<float>(mBx1inv), kpred);
		
			  vector<float> v_cmean;
			  v_cmean=bTransform(vars,mBx2,M);
			
			  for (unsigned int i=0;i<vars.size();i++)
			    if (i<(unsigned)nmodes){
			      vars.at(i)=v_cmean.at(i);
			    }else{
			      vars.at(i)=0;
			    }
			 
		  }

		
		
		
		//-----------------------------------REGISTER MAIN MODEL AND FIT MAIN MODEL------------------------------//
		//need to implement new reg method, first apply reg can calculate SVD of new modes
		if (verbose.value()) cout<<"Registering model."<<endl;
		
		model1->registerModel(fmatv);
		
		if (verbose.value()) cout<<"Model has been registered to native space now fitting model to image."<<endl;
		
		//		for (int i=0;i<vars.size();i++)
		  //		  {
		    //		    for (int j=0;j<select.size();j++)
		      //		      if (j==i)
			//			select.at(j)=true;
		    //		      else
			//			select.at(j)=false;
			                       	    

		conjGradient(image, *model1,vars, relStd, select, searchRes, 0.15);
		//	  }
		if (verbose.value()) 
		{
			cout<<"Final mode parameters are: "<<endl;
			for (unsigned int i=0; i<MaxModes;i++){
				cout<<vars.at(i)<<" ";
			}
			cout<<endl;
		}
		
		//----------------------------------FILL IMAGE AND WRITE OUTPUT------------------------------//
		
		if (verbose.value()) cout<<"Get deformed surface and fill."<<endl;
		
		//write output image	
		vector<float> shape = model1->getDeformedGrid(vars);
		int bounds[6]={0,0,0,0,0,0};
		getBounds(shape,bounds,image.xdim(),image.ydim(),image.zdim());
		volume<short> mask=make_mask_from_mesh(image, shape, model1->cells,model1->getLabel(0), bounds);
		string name_out=*out_iter;
		if (multiImageInput)
			name_out+=outname;
	
		save_volume(mask, name_out);
	
		vector<float> bvars_old=model1->getOrigSpaceBvars(vars);
				
		write_bvars(*im_iter, modelname, bvars_old, name_out+".bvars",fmatv_org);
		write_vtk(shape, model1->cells, name_out + ".vtk", is_binary);


		if (verbose.value()) cout<<"FIRST has completed subject "<<*im_iter<<endl;
	}
	
		delete model1;
		if (intref) delete modelRef;
	
	if (verbose.value()) cout<<"FIRST has completed for all subjects."<<endl;
	
	return 0;
}




int main(int argc,char *argv[])
{
	
	Tracer tr("main");
	OptionParser options(title, examples);
	
	try {
		// must include all wanted options here (the order determines how
		//  the help message is printed)
		options.add(inname);
		options.add(outname);
		options.add(verbose);
		options.add(help);
		options.add(modelname); 
		options.add(modelname2); 
		options.add(flirtmatname);
		options.add(nmodes); 
		options.add(intref);
		options.add(multiImageInput);
		options.add(binarySurfaceOutput);
		options.add(bmapname);
		options.add(bvarsname);
		options.add(shcond);
		options.add(loadbvars);
		nonoptarg = options.parse_command_line(argc, argv);
		
		// line below stops the program if the help was requested or 
		//  a compulsory option was not set
				if (  (!options.check_compulsory_arguments(true) ))
		{
			options.usage();
			exit(EXIT_FAILURE);
		}
		
		// Call the local functions

		do_work(inname.value(),modelname.value(),modelname2.value(), flirtmatname.value(), outname.value(), bmapname.value(), bvarsname.value(), \
				nmodes.value(), intref.value(), multiImageInput.value(), shcond.value(), loadbvars.value(), binarySurfaceOutput.value(), 0.5);
	
		
	}  catch(X_OptionError& e) {
		options.usage();
		cerr << endl << e.what() << endl;
		exit(EXIT_FAILURE);
	} catch(std::exception &e) {
		cerr << e.what() << endl;
	} 
	
	return 0;// do_work(argc,argv);
}

