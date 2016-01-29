/*
 *  shapeModel.cpp
 *  
 *
 *  Created by Brian Patenaude on 23/06/2008.
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
#include "shapeModel.h"
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdio.h>
//#include <vector>
 #include "newmat.h"
 #include "newmatap.h"
#include "fslvtkio/fslvtkio.h"
#include "first_lib/first_newmat_vec.h"
#include <cmath>
#include <algorithm>
#include "math.h"
using namespace std;
using namespace NEWMAT;
using namespace FIRST_LIB;
using namespace fslvtkio;
namespace SHAPE_MODEL_NAME{
	
	shapeModel::shapeModel()
	{
	
	}
	shapeModel::shapeModel( const vector<float> & mshape, const vector< vector<float> > & modesshape, const vector<float> & se, \
							const vector<float> & ishape, const vector< vector<float> > & modesint, const vector<float> & ie, \
							const int & M, const vector<float> & errs)
{
		smean=mshape;
		smodes=modesshape;
		seigs=se;
		sqrtseigs=se;

		for (vector<float>::iterator i=sqrtseigs.begin();i!=sqrtseigs.end();i++)
			*i = sqrt(*i);
		sqrtseigsi=sqrtseigs;


		imean=ishape;
		imodes=modesint;
		ieigs=ie;
		NumberOfSubjects=M;
		Errs=errs;
		USE_COND=false;
		MODE_FOUND=false;
		mode=0;
	STORE_REG_XFM=true;
		
}

shapeModel::shapeModel( const vector<float> & mshape, const vector< vector<float> > & modesshape, const vector<float> & se, \
						const vector<float> & ishape, const vector< vector<float> > & modesint, const vector< vector<float> > & Iprec, const vector<float> & ie,\
						const int & M, const vector<float> & errs, const vector< vector<unsigned int> > & cellsin, const vector<int> & vlabels )
{
	labels=vlabels;
	smean=mshape;
	smodes=modesshape;
	seigs=se;
	sqrtseigs=se;
	
	for (vector<float>::iterator i=sqrtseigs.begin();i!=sqrtseigs.end();i++)
		*i = sqrt(*i);
	sqrtseigsi=sqrtseigs;
	imean=ishape;
	imodes=modesint;
	i_precision=Iprec;
	ieigs=ie;
	NumberOfSubjects=M;
	Errs=errs;
	cells=cellsin;
	USE_COND=false;
	MODE_FOUND=false;
	mode=0;
	STORE_REG_XFM=true;

	//keeps tracks of neighbourin triangles
	for (unsigned int i=0; i<static_cast<unsigned int>(smean.size()/3); i++)
	{
		vector<unsigned int> inds;
		int count=0;
		for ( vector< vector<unsigned int> >::iterator j=cells.begin(); j!=cells.end();j++,count++)
		{
			for ( vector<unsigned int>::iterator k= j->begin(); k!=j->end();k++)
			{
				if ((*k)==i)
				{
					inds.push_back(count); 
					break;
				}
			}
		}
		localTri.push_back(inds);
	}
	
	
}

shapeModel::shapeModel( const vector<float> & mshape, const vector< vector<float> > & modesshape, const vector<float> & se, \
						const vector<float> & ishape, const vector< vector<float> > & modesint, const vector<float> & ie, \
						const int & M, const vector<float> & errs, const vector<short> & vmaskin)
{
	smean=mshape;
	smodes=modesshape;
	seigs=se;
	imean=ishape;
	imodes=modesint;
	ieigs=ie;
	NumberOfSubjects=M;
	Errs=errs;
	stmask=vmaskin;
	MODE_FOUND=false;
	mode=0;
	STORE_REG_XFM=true;

}
	void shapeModel::printLabel( const unsigned int & i) const 
	{ cout<<"get labvel "<< labels.at(i) <<endl;}
	
	
//	int shapeModel::getLabel( const unsigned int & i) const 
//	{ cout<<"get labvel "<< labels.at(i) <<endl; return labels.at(i); }


std::vector<float> shapeModel::getDeformedGrid( const std::vector<float>  & vars ) const 
{

	std::vector<float> newshape=smean;
	std::vector<float>::const_iterator sqrtseigs_i=sqrtseigs.begin();
	std::vector< std::vector<float> >::const_iterator smodes_i = smodes.begin();
	for (std::vector<float>::const_iterator vars_i=vars.begin(); vars_i!=vars.end(); vars_i++,sqrtseigs_i++, smodes_i++)
	{
		std::vector<float>::iterator new_i=newshape.begin();
		for (std::vector<float>::const_iterator smodes_j= smodes_i->begin(); smodes_j!= smodes_i->end(); smodes_j++,new_i++)
			*new_i += (*vars_i) * (*sqrtseigs_i) * (*smodes_j);//* 
	}
	
	return newshape;
}

vector<float> shapeModel::getDeformedIGrid( const vector<float> & vars) const {

	vector<float> varsnew=vars;//getOrigSpaceBvars(vars);
	vector<float> newigrid=imean;
	for (unsigned int i=0; i< varsnew.size();i++){
		for (unsigned int j=0; j< imean.size(); j++){
		  newigrid.at(j)+=varsnew.at(i)*sqrtseigs.at(i)*imodes.at(i).at(j);
		}
	}
	return newigrid;
}

	vector< vector<float> > shapeModel::registerModeVectors( const vector< vector<float> >& vmodes, const vector< vector<float> >& flirtmat )
	{
		vector< vector<float> > reg_modes = vmodes;
		for (unsigned int i=0; i<vmodes.size();i++){
			for (unsigned int j=0; j<vmodes.at(0).size();j+=3){
				//exclude translation from mode of variation so that it is compatible with implementation	
				float x= flirtmat.at(0).at(0)*vmodes.at(i).at(j)  + flirtmat.at(0).at(1)*vmodes.at(i).at(j+1) + flirtmat.at(0).at(2)*vmodes.at(i).at(j+2);// +  flirtmat.at(0).at(3);
				float y= flirtmat.at(1).at(0)*vmodes.at(i).at(j)  + flirtmat.at(1).at(1)*vmodes.at(i).at(j+1) + flirtmat.at(1).at(2)*vmodes.at(i).at(j+2) ;//+  flirtmat.at(1).at(3);
				float z= flirtmat.at(2).at(0)*vmodes.at(i).at(j)  + flirtmat.at(2).at(1)*vmodes.at(i).at(j+1) + flirtmat.at(2).at(2)*vmodes.at(i).at(j+2) ;//+  flirtmat.at(2).at(3);
				
				 reg_modes.at(i).at(j)=x;
				 reg_modes.at(i).at(j+1)=y;
				 reg_modes.at(i).at(j+2)=z;	
			}
			
		}
		return reg_modes;
	}
	
	void shapeModel::registerModel(const vector< vector<float> > & flirtmat)
	{
		//start by registering mesh 
		float f_11= flirtmat.at(0).at(0);
		float f_12= flirtmat.at(0).at(1);
		float f_13= flirtmat.at(0).at(2);
		float f_14= flirtmat.at(0).at(3);
		float f_21= flirtmat.at(1).at(0);
		float f_22= flirtmat.at(1).at(1);
		float f_23= flirtmat.at(1).at(2);
		float f_24= flirtmat.at(1).at(3);
		float f_31= flirtmat.at(2).at(0);
		float f_32= flirtmat.at(2).at(1);
		float f_33= flirtmat.at(2).at(2);
		float f_34= flirtmat.at(2).at(3);
		//assumes last line is 0 0 0 1
		for (vector<float>::iterator i=smean.begin();i!=smean.end();i+=3)
		{
			float x= f_11 * (*i)  + f_12 * (*(i+1)) + f_13 * (*(i+2)) +  f_14;
			float y= f_21 * (*i)  + f_22 * (*(i+1)) + f_23 * (*(i+2)) +  f_24;
			float z= f_31 * (*i)  + f_32 * (*(i+1)) + f_33 * (*(i+2)) +  f_34;
		
			(*i)=x;
			(*(i+1))=y;
			(*(i+2))=z;				
		}
		
		vector< vector<float> > smodesAinv;
		//Matrix M_smodes_mni= first_newmat_vector::vectorOfVectorsToMatrix<float>(smodes);
		Matrix M_smodesPre= first_newmat_vector::vectorOfVectorsToMatrix<float>(smodes).t();
		smodes=registerModeVectors(smodes,flirtmat);
		
		Matrix M_smodes= first_newmat_vector::vectorOfVectorsToMatrix<float>(smodes).t();

			Matrix U,U2,V;
			DiagonalMatrix D,D2;
		Matrix Mfmat(3,3);
		for (int i=0;i<3;i++)
		{
			for (int j=0;j<3;j++)
			{
				Mfmat.element(i,j)=flirtmat.at(i).at(j);
				cout<<flirtmat.at(i).at(j)<<" ";
			}
			cout<<endl;
		}

		Mfmat=Mfmat.i();

		//invert an transpose
		Matrix A(M_smodes.Nrows(),M_smodes.Nrows());
		A=0;
		for (int i=0;i<A.Nrows()/3;i++)
			A.SubMatrix(3*i+1, 3*(i+1),3*i+1, 3*(i+1))=Mfmat; 
		
		
		if (STORE_REG_XFM)
		{
			u_xfm = first_newmat_vector:: matrixToVectorOfVectors<float>(M_smodesPre.t()*A);
			d_xfm=sqrtseigs;
		}
		
			SVD(M_smodes,D,U,V);
			DiagonalMatrix Eigs= first_newmat_vector::vectorToDiagonalMatrix(seigs);
			
			//DiagonalMatrix SqrtEigs= first_newmat_vector::vectorToDiagonalMatrix(sqrtseigs);
			SVD(D*V.t()*Eigs*V*D,D2,U2);

//		cout<<"st ore reg xfm "<<A.Nrows()<<" "<<A.Ncols()<<" "<<M_smodes.Nrows()<<" "<<M_smodes.Ncols()<<endl;
	
		
		
			smodes= first_newmat_vector::matrixToVector<float>((U*U2).SubMatrix(1,U.Nrows(),1,D2.Nrows()));
			vector<float> veigs, vsqrt_eigs;
			for (unsigned int i=0;i<static_cast<unsigned int>(D2.Nrows());i++)//static_cast<unsigned int>(D2.Nrows());i++)
			{
				veigs.push_back(D2.element(i));
				vsqrt_eigs.push_back(sqrt(D2.element(i)));
			}			
			seigs=veigs;
			sqrtseigs=vsqrt_eigs;

		
	

		
		
		Matrix M_imodes= first_newmat_vector::vectorOfVectorsToMatrix<float>(imodes).t();

//		cout<<"Transform imodes "<<M_imodes.Nrows()<<" "<<M_imodes.Ncols()<<" "<<M_smodes.Nrows()<<" "<<M_smodes.Ncols()<<A.Nrows()<<" "<<A.Ncols()<<endl;
		Matrix C=M_smodesPre.t()*A;
//		cout<<"Transform imodes "<<M_imodes.Nrows()<<" "<<M_imodes.Ncols()<<" "<<M_smodes.Nrows()<<" "<<M_smodes.Ncols()<<A.Nrows()<<" "<<A.Ncols()<<endl;

		C=C*U;
//cout<<"Transform imodes "<<M_imodes.Nrows()<<" "<<M_imodes.Ncols()<<" "<<M_smodes.Nrows()<<" "<<M_smodes.Ncols()<<A.Nrows()<<" "<<A.Ncols()<<endl;
//		cout<<"Transform imodes "<<M_imodes.Nrows()<<" "<<M_imodes.Ncols()<<" "<<M_smodes.Nrows()<<" "<<M_smodes.Ncols()<<A.Nrows()<<" "<<A.Ncols()<<endl;

		imodes= first_newmat_vector::matrixToVector<float>((M_imodes*C*U2));
			
		cout<<"NEw done imodes transform"<<endl;
 

	}
	std::vector<float> shapeModel::getOrigSpaceBvars(const std::vector<float> & bvars ) const
	{
		vector<float> bvars_orig;
		if (STORE_REG_XFM)
		{	
			DiagonalMatrix D2=first_newmat_vector::vectorToDiagonalMatrix(sqrtseigs);
			DiagonalMatrix D=first_newmat_vector::vectorToDiagonalMatrix(d_xfm);
			Matrix UU2 = first_newmat_vector::vectorOfVectorsToMatrix<float>(u_xfm);
			Matrix M_smodes= first_newmat_vector::vectorOfVectorsToMatrix<float>(smodes).t();
			/*
			
			cout<<"Smodes2 "<<UU2.Ncols()<<endl;
			for (int i=0;i<UU2.Ncols();i++)
			{
				cout<<UU2.element(0,i)<<" ";
				cout<<endl;
			}
			
		*/	
			
			ColumnVector Mbvars(D.Nrows()); 
			int count=0;
			Mbvars=0;
			for (vector<float>::const_iterator i=bvars.begin();i!=bvars.end();i++,count++)
				Mbvars.element(count)=*i;
			
			cout<<"D "<<endl;//<<"di"<<endl<<endl;
	//		for (int i=0;i<D.Nrows();i++)
	//		{
			//	if (abs(D.element(i)) < 1e-10)
		//			D.element(i)=0;
		//		else 
	//			cout<<"pre "<<i<<" "<<D.element(i)<<" "<<1/D.element(i)<<endl;
//
//					D.element(i)=1/sqrt(D.element(i));
//				cout<<i<<" "<<D.element(i)<<" "<<1/D.element(i)<<endl;
//			}
			cout<<UU2.Nrows()<<" "<<UU2.Ncols()<<" "<<M_smodes.Nrows()<<" "<<M_smodes.Ncols()<<endl;
			ColumnVector MbvarsOrg = D.i()*UU2*M_smodes*D2*Mbvars;
			
			for (int i=0;i<MbvarsOrg.Nrows();i++)
				bvars_orig.push_back(MbvarsOrg.element(i));
			
		}
		return bvars_orig;
	}
	
	
	shapeModel* shapeModel::loadAndCreateShapeModel( const string & modelname, const bool & verbose)
	{
		
		if (verbose) cout<<"read model"<<endl;
		fslvtkIO* fmodel = new fslvtkIO(modelname,static_cast<fslvtkIO::DataType>(0));
		if (verbose) cout<<"done reading model"<<endl;
		const int Npts=fmodel->getPointsAsMatrix().Nrows();
		
		
		if (verbose) cout<<"setting up shape/appearance model"<<endl;
		
		//load mean shape into vector
		vector<float> Smean;
		Matrix* Pts = new Matrix;
		*Pts=fmodel->getPointsAsMatrix();
		
		if (verbose) cout<<"The shape has "<<Npts<<" vertices."<<endl;
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
		unsigned int M = static_cast<unsigned int>(fmodel->getField("numSubjects").element(0,0));
		if (verbose) cout<<"The model was constructed from "<<M<<" training subjects."<<endl;
		
		unsigned int AllModes=M;
		
		for (unsigned int i =1; i<AllModes;i++)
		{
			stringstream ss;
			ss<<i;
			string mode;
			ss>>mode;
			SmodesM=SmodesM | first_newmat_vector::unwrapMatrix(fmodel->getField("mode"+mode));
			ImodesM=ImodesM | first_newmat_vector::unwrapMatrix(fmodel->getField("Imode"+mode));
			
		}
		if (verbose) cout<<AllModes<<" modes of variation are retained."<<endl;
		
		
		vector< vector<float > > Smodes = first_newmat_vector::matrixToVector<float>(SmodesM);
		vector< vector<float > > Imodes = first_newmat_vector::matrixToVector<float>(ImodesM);
		ImodesM.Release();
		SmodesM.Release();
		
		
		//process rest of information, including intensity variance
		vector< vector<float > > Iprec = first_newmat_vector::matrixToVector<float>(fmodel->getField("iCondPrec0").t());
		vector<float > Errs =  first_newmat_vector::vectorToVector<float>(fmodel->getField("ErrPriors0"));
		vector<float > se =  first_newmat_vector::vectorToVector<float>(fmodel->getField("eigenValues"), AllModes);
		vector<float > ie =  first_newmat_vector::vectorToVector<float>(fmodel->getField("iCondEigs0"));
		vector<float > Imean =  first_newmat_vector::vectorToVector<float>(first_newmat_vector::unwrapMatrix(fmodel->getField("Imean")));
		vector<int> labels =  first_newmat_vector::vectorToVector<int>(fmodel->getField("labels"));
		
		
		
		//have read in all data and store in local structures, now delete the reader.
		delete fmodel;
		cout<<"create shapeModel "<<endl;
		//create shape model
		shapeModel* model1 = new shapeModel(Smean, Smodes, se, Imean, Imodes,Iprec, ie, M,Errs,polygons,labels);
		cout<<"done creating shapeModel "<<endl;
		
		return model1;
		
	}
	
	
}
