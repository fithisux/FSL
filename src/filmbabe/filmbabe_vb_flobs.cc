/*  filmbabe_vb_flobs.cc

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

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

#include "filmbabe_vb_flobs.h"
#include "utils/log.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include "libvis/miscplot.h"
#include "libvis/miscpic.h"
#include "newmat.h"
#include "utils/tracer_plus.h"
#include "filmbabeoptions.h"
#include "newimage/newimagefns.h"
#include <utility>
#include "libprob.h"
#include "miscmaths/sparsefn.h"

using namespace NEWIMAGE;
using namespace Utilities;
using namespace MISCMATHS;
using namespace NEWMAT;
using namespace MISCPLOT;

namespace Filmbabe {

  Filmbabe_Vb_Flobs::Filmbabe_Vb_Flobs(const volume4D<float>& pdata, const volume<int>& pmask, const Matrix& pdesignmatrix, const ColumnVector& pflobsregressors, const volume4D<float>& plocalweights, const vector<Connected_Offset>& pconnected_offsets, int pnum_superthreshold) :
    xsize(pdata.xsize()),
    ysize(pdata.ysize()),
    zsize(pdata.zsize()),
    ntpts(pdata.tsize()),
    data(pdata),
    mask(pmask),
    designmatrix(pdesignmatrix),
    flobsregressors(pflobsregressors),
    QtXXt(),
    localweights(plocalweights),
    connected_offsets(pconnected_offsets),
    indices(pdata.xsize(),pdata.ysize(),pdata.zsize()),
    Y(),
    D(),
    trace_ilambdaDA(FilmbabeOptions::getInstance().ntar.value()),
    trace_ilambdaBeta(pnum_superthreshold),
    trace_ilambdaXXt(pnum_superthreshold),
    diag_ilambdaA(FilmbabeOptions::getInstance().ntar.value()),
    m_Beta(pnum_superthreshold),
    ilambda_Beta(pnum_superthreshold),
    m_A(FilmbabeOptions::getInstance().ntar.value()),
    gam_Beta(pnum_superthreshold),
    gam_A(FilmbabeOptions::getInstance().ntar.value()),
    gam_e(pnum_superthreshold),
    num_superthreshold(pnum_superthreshold),
    maxnumneighs(8),
    ntar(FilmbabeOptions::getInstance().ntar.value()),
    ilambdaDA(FilmbabeOptions::getInstance().ntar.value()),
    ilambdaA(FilmbabeOptions::getInstance().ntar.value()),    
    realevcoord(),
    Q()
  {    
  }

  void Filmbabe_Vb_Flobs::process_flobsregressors()
  {
    Tracer_Plus trace("Filmbabe_Vb_Flobs::process_flobsregressors");

    
    // total number of actual real evs
    nrealevs = flobsregressors.Nrows();

    // assumes real evs belonging to an original flob ev come in contiguously
    // want design matrix reordered so that flob regressors come first 
    // must have same number of basis fns for each original flob ev  
    Matrix designmatrixnew = designmatrix;
    
    nflobsevs=0;// num flobs Original evs
    nnonflobsevs=0;// num of non flobs evs
    nbfs=0;// num basis functions for each flobs Original evs (is the same for all Original evs)
    vector<int> nonflobsindex;
    vector<int> flobsindex;

    // find number of non flobs regressors and put flobs regressors to the start of dm
    int j=1;
    for(int i=1; i<=nrealevs; i++)
      {
	if(flobsregressors(i) == -1)
	  {
	    // non flob ev
	    nonflobsindex.push_back(i);
	  }
	else
	  {
	    designmatrixnew.Column(j) = designmatrix.Column(i);
	    flobsindex.push_back(i);
	    j++;
	  }
      }

    nnonflobsevs = nonflobsindex.size();

    // put non flobs regressors to the end
    for(int i=1; i<=nnonflobsevs; i++)
    {
      designmatrixnew.Column(j) = designmatrix.Column(nonflobsindex[i-1]);
      j++;
    }
    designmatrix = designmatrixnew;

    // output reordered designmatrix
    write_ascii_matrix(designmatrix, LogSingleton::getInstance().appendDir("designmatrixnew"));

    // find nbfs
    int ind = int(flobsregressors(flobsindex[0]));

    for(unsigned int i=1; i<=flobsindex.size(); i++)
      {
	if(flobsregressors(flobsindex[i-1]) == ind)
	  nbfs++;
      }

    nflobsevs = int((nrealevs - nnonflobsevs)/nbfs);

    OUT(nnonflobsevs);
    OUT(nrealevs);
    OUT(nbfs);    
    OUT(nflobsevs);

    if(nrealevs != nbfs*nflobsevs + nnonflobsevs)
      throw Exception("Invalid FLOBS regressors file");

    // setup flobs constraints:    
    if(!FilmbabeOptions::getInstance().flobsprioroff.value())
      {	
	if(FilmbabeOptions::getInstance().flobsdir.value()!=string(""))
	  {
	    m_Beta_0_global = read_vest(FilmbabeOptions::getInstance().flobsdir.value()+"/priormeans.mat");
	    
	  }
	else
	  m_Beta_0_global = read_vest(FilmbabeOptions::getInstance().priormeanfile.value()).AsColumn();
	
	if(FilmbabeOptions::getInstance().verbose.value())
	  OUT(m_Beta_0_global);
	       
	Matrix lambdatmp;
	if(FilmbabeOptions::getInstance().flobsdir.value()!=string(""))
	  lambdatmp = read_vest(FilmbabeOptions::getInstance().flobsdir.value()+"/priorcovars.mat");
	else
	  lambdatmp = read_vest(FilmbabeOptions::getInstance().priorcovarfile.value());

	lambda_Beta_0_global << lambdatmp;

	if(FilmbabeOptions::getInstance().verbose.value())
	  OUT(lambda_Beta_0_global);

	lambda_Beta_0_global = lambda_Beta_0_global.i();

	if(lambda_Beta_0_global.Nrows() != m_Beta_0_global.Nrows() || nbfs != m_Beta_0_global.Nrows())
	   throw Exception("Invalid FLOBS constraints");

      }
    else
      {
	m_Beta_0_global.ReSize(nbfs);
	m_Beta_0_global = 0;
	lambda_Beta_0_global.ReSize(nbfs);
	lambda_Beta_0_global = 0;
	
	for(int b=1; b<=nbfs; b++)
	  {
	    m_Beta_0_global(b) = 0.0001;	
	    lambda_Beta_0_global(b,b) = 1;
	  }
      }

  }

  void Filmbabe_Vb_Flobs::setup()
  {
    Tracer_Plus trace("Filmbabe_Vb_Flobs::setup");   
    
    OUT("Setup");
    process_flobsregressors();

    niters = FilmbabeOptions::getInstance().niters.value();
    OUT(niters);
    OUT(num_superthreshold);

    QemReCquad.resize(nflobsevs);
    Re.resize(nflobsevs);
    Qe.resize(nflobsevs);
    m_Beta_0.resize(nflobsevs);
    lambda_Beta_0.resize(nflobsevs);

    // create squashed voxels vector
    voxels.reserve(num_superthreshold);
    for(int z = 0; z < data.zsize(); z++)    
      for(int y = 0; y < data.ysize(); y++)	
	for(int x = 0; x < data.xsize(); x++)
	  if(mask(x,y,z))
	    {
	      //if(MISCMATHS::var(data.voxelts(x,y,z)).AsScalar()>1e-10)
	      voxels.push_back(Voxel(x,y,z));	      
	    }

    // setup D matrix
    if(FilmbabeOptions::getInstance().verbose.value())
      OUT("Setup D");

    // num neighbours
    ColumnVector num_neigbours(num_superthreshold);
    num_neigbours = 0;
    indices = 0;
    int index=1;
    
    for(int z = 0; z < data.zsize(); z++)    
      for(int y = 0; y < data.ysize(); y++)	
	for(int x = 0; x < data.xsize(); x++)
	  {

	    if(mask(x,y,z))
	      {
		int xi=0,yi=0,zi=0;
		for(unsigned int i = 0; i < connected_offsets.size(); i++) 
		  {
		    xi = x+connected_offsets[i].x;
		    yi = y+connected_offsets[i].y;
		    zi = z+connected_offsets[i].z;
		    
		    if(mask(xi,yi,zi))
		      {
			num_neigbours(index) += localweights(x,y,z,connected_offsets[i].ind);
		      }
		  }
		
		indices(x,y,z) = index;
		index++;
	      }
	  }

    D.ReSize(num_superthreshold,num_superthreshold);
    OUT("Setup D2");

    for(int z = 0; z < data.zsize(); z++)
      for(int y = 0; y < data.ysize(); y++)
	for(int x = 0; x < data.xsize(); x++)
	  if(mask(x,y,z))
	    {
	      int xi=0,yi=0,zi=0;
	      for(unsigned int i = 0; i < connected_offsets.size(); i++) 
		{
		  xi = x+connected_offsets[i].x;
		  yi = y+connected_offsets[i].y;
		  zi = z+connected_offsets[i].z;
		  
		  if(mask(xi,yi,zi))
		    {  
		      D.insert(indices(x,y,z),indices(xi,yi,zi), -1.0/sqrt(num_neigbours(indices(x,y,z))*num_neigbours(indices(xi,yi,zi))));
		      D.insert(indices(xi,yi,zi), indices(x,y,z), -1.0/sqrt(num_neigbours(indices(x,y,z))*num_neigbours(indices(xi,yi,zi))));
		    }
		}

	      D.insert(indices(x,y,z),indices(x,y,z),1);
	    }

    if(FilmbabeOptions::getInstance().verbose.value())
      OUT("Setup Y");

    // Y is T*N, column = spatialmap, row = time-series.
    Y.ReSize(ntpts,num_superthreshold);    
    Y = 0;
    for(int r=0; r < num_superthreshold; r++)
      {	
	const Voxel& vox = voxels[r]; 
	ColumnVector tmp = data.voxelts(vox.x,vox.y,vox.z);
	Y.Column(r+1) = tmp - mean(tmp).AsScalar();
	gam_e(r+1) = 1.0/var(tmp).AsScalar();
      }
   
    if(FilmbabeOptions::getInstance().verbose.value())
      OUT("Initialise params");

    for(int p=1; p<=ntar; p++)
      {
	if(FilmbabeOptions::getInstance().tarmrfprec.value()==-1)
	  gam_A(p) = 1;
	else     
	  gam_A(p) = FilmbabeOptions::getInstance().tarmrfprec.value();
	
	m_A[p-1].ReSize(num_superthreshold);
	m_A[p-1] = 0;
      }

//     gam_e = 0.9;

    // realevcoord gives index of flobs original ev e and basis fn b in beta vector
    realevcoord.ReSize(nflobsevs,nbfs);
    realevcoord = 0;

    if(FilmbabeOptions::getInstance().verbose.value())
      OUT("Initialise FLOBS params");

    int coord = 0;
    for(int e=1; e<=nflobsevs; e++)
      {	
	m_Beta_0[e-1] = m_Beta_0_global;
	lambda_Beta_0[e-1] = lambda_Beta_0_global;	

	for(int b=1; b<=nbfs; b++)
	  {	    
	    coord++;
	    realevcoord(e,b) = coord;
	  }	
      }

    if(FilmbabeOptions::getInstance().verbose.value())
      OUT("Initialise remaining params");

    for(int i=1; i<=num_superthreshold; i++)
      {
	for(int p=1; p<=ntar; p++)
	  m_A[p-1](i) = normrnd().AsScalar()*0.1;

	m_Beta[i-1].ReSize(nrealevs+nflobsevs);
	m_Beta[i-1] = 0;

	gam_Beta[i-1].ReSize(nflobsevs);
	gam_Beta[i-1] = 1;

	ilambda_Beta[i-1].ReSize(nrealevs+nflobsevs);
	ilambda_Beta[i-1] = 0;

	trace_ilambdaBeta[i-1].ReSize(nflobsevs);
	trace_ilambdaBeta[i-1] = 0;
	
//      initialise m_Beta using OLS estimates
//  	ColumnVector betatmp = pinv(designmatrix)*Y.Column(i);
//  	ColumnVector res = Y.Column(i) - designmatrix*betatmp;
//  	gam_e(i) = 1.0/(res.t()*res/(ntpts-nrealevs));
//  	ilambda_Beta[i-1] = pinv(designmatrix)
//  	m_Beta[i-1] = 1;

	for(int e=1; e<=nrealevs; e++)
	  {
	    m_Beta[i-1](e) = normrnd().AsScalar()*0.1;
	    //m_Beta[i-1](e) = betatmp(e);
	  }
	
      }

    for(int p=1; p<=ntar; p++)
      {
	ilambdaDA[p-1].ReSize(num_superthreshold,num_superthreshold);
	ilambdaA[p-1].ReSize(num_superthreshold,num_superthreshold);

	diag_ilambdaA[p-1].ReSize(num_superthreshold);
	diag_ilambdaA[p-1] = 0;
      }

    trace_ilambdaDA = 0;
    trace_ilambdaXXt = 0;

    // beta = {real flobs overall params; real non flobs params; real flobs "size" params}
    // beta has size nrealevs+nflobsevs

    // setup Q stuff
    // Q*beta = {real flobs overall params; real non flobs params} ==> i.e. Y=X*Q*beta+e
    // Q[e]*beta = {real flobs overall param for original flobev e}
    // Re[e]*beta = {real flobs "size" param for original flobev e}
    Q.ReSize(nrealevs,nrealevs+nflobsevs);
    Q = 0;
    for(int e=1; e<=nflobsevs; e++)
      {
	Qe[e-1].ReSize(nbfs,nrealevs+nflobsevs);
	Qe[e-1] = 0;
	Re[e-1].ReSize(nrealevs+nflobsevs);
	Re[e-1] = 0;
	Re[e-1](nrealevs+e) = 1;

	for(int b=1; b<=nbfs; b++)
	  {
	    Q(int(realevcoord(e,b)),int(realevcoord(e,b))) = 1;
	    Qe[e-1](b,int(realevcoord(e,b))) = 1;
	  }
	QemReCquad[e-1] = (Qe[e-1]-m_Beta_0[e-1]*Re[e-1]).t()*lambda_Beta_0[e-1]*(Qe[e-1]-m_Beta_0[e-1]*Re[e-1]);
      }

    for(int e=1; e<=nnonflobsevs; e++)
      {
	Q(nflobsevs*nbfs+e,nflobsevs*nbfs+e) = 1;
      }

    designmatrixQ.ReSize(nrealevs+nflobsevs,ntpts);
    designmatrixQ = 0;
    for(int t=1; t <= ntpts; t++)
      designmatrixQ.Column(t) = (designmatrix.Column(t).AsRow()*Q).AsColumn();

    QtXXt = Q.t()*designmatrix*designmatrix.t();

    if(FilmbabeOptions::getInstance().verbose.value())
      {
	OUT(Q);
	OUT(Qe[0]);
	OUT(Re[0]);
	OUT(QemReCquad[0]);
      }

    OUT("Setup finished");
  }

  void Filmbabe_Vb_Flobs::run()
  {
    Tracer_Plus trace("Filmbabe_Vb_Flobs::run");    

    int i = 1;;
    for(; i<=niters; i++)
      {
	OUT(i);
 	update_Beta();

	update_A();
 	if(FilmbabeOptions::getInstance().tarmrfprec.value()==-1)
 	  update_phiA();

	update_phie();

	for(int i=1; i<=num_superthreshold; i++)
	  for(int e=1; e<=nflobsevs; e++)
	    {		
	      gam_Beta[i-1](e) = 1.0/Sqr((Re[e-1]*m_Beta[i-1]).AsScalar()+ilambda_Beta[i-1](nrealevs+e,nrealevs+e));
	      
	      if(gam_Beta[i-1](e)>1e10)gam_Beta[i-1](e)=1e10;
	      if(gam_Beta[i-1](e)<1e-10)gam_Beta[i-1](e)=1e-10;
	      //OUT(gam_Beta[i-1](e));
	    }
      }

    cout << "Iterations=" << i << endl;
  }

  void Filmbabe_Vb_Flobs::update_Beta()
  {
    Tracer_Plus trace("Filmbabe_Vb_Flobs::update_Beta");

    if(FilmbabeOptions::getInstance().verbose.value())
      OUT("Update Beta");

    Matrix F_Beta(nrealevs+nflobsevs,nrealevs+nflobsevs);
    SymmetricMatrix F_Beta_sym(nrealevs+nflobsevs);
    ColumnVector E_Beta(nrealevs+nflobsevs);

    vector<SymmetricMatrix> etmp(nflobsevs);
    for(int e=1; e<=nflobsevs; e++)
      {
	Matrix tmp = (Qe[e-1]-m_Beta_0[e-1]*Re[e-1]);
	etmp[e-1] << tmp.t()*lambda_Beta_0[e-1]*tmp;
      }

    for(int j = 1; j<=num_superthreshold; j++)
      {
	// cout << j << ",";
	// cout.flush();

	F_Beta = 0;
	E_Beta = 0;	
	
	for(int t=1; t <= ntpts; t++)
	  {
	    ColumnVector sumpx = designmatrixQ.Column(t);

	    float sumpy = Y(t,j);

	    for(int p=1; p<=ntar; p++)
	      if(t>p)
		{
		  //designmatrixQ.Column(t-p) is (designmatrix.Column(t-p).AsRow()*Q).AsColumn();
		  sumpx -= m_A[p-1](j)*designmatrixQ.Column(t-p);
		  sumpy -= m_A[p-1](j)*Y(t-p,j);
		}
	    F_Beta += sumpx*sumpx.t();
	    E_Beta += sumpx*sumpy;
	  }

	F_Beta *= gam_e(j);
	E_Beta *= gam_e(j);

	for(int e=1; e<=nflobsevs; e++)
	  {
//  	    Matrix tmp = (Qe[e-1]-m_Beta_0[e-1]*Re[e-1]);
//  	    F_Beta += tmp.t()*lambda_Beta_0[e-1]*tmp*gam_Beta[j-1](e);
	
	    F_Beta += etmp[e-1]*gam_Beta[j-1](e);
	  }
	
	F_Beta_sym << F_Beta;

	ilambda_Beta[j-1] = F_Beta_sym.i();
	m_Beta[j-1] =  ilambda_Beta[j-1]*E_Beta;
	trace_ilambdaXXt(j) = (Q*ilambda_Beta[j-1]*QtXXt).Trace();

	//	OUT(ilambda_Beta[j-1]);
//  	for(int e=1; e<=neflobsvs; e++)
//  	  {
//  	    trace_ilambdaBeta[j-1](e) = (ilambda_Beta[j-1]*QemReCquad[e-1]).Trace();
//  	  }

      }
  }

  void Filmbabe_Vb_Flobs::update_A()
  {
    Tracer_Plus trace("Filmbabe_Vb_Flobs::update_A");

    SparseMatrix I;
    speye(num_superthreshold,I);

    if(FilmbabeOptions::getInstance().verbose.value())
      OUT("Update A");

    for(int p=1; p<=ntar; p++)
      {
	ColumnVector beta(num_superthreshold);
	beta = 0;
	SparseMatrix lambda_A = D;
	if(FilmbabeOptions::getInstance().verbose.value())
	  OUT("Compute lambda");

	if(FilmbabeOptions::getInstance().verbose.value())
	  {
	    OUT(p);
	    OUT(gam_A(p));
	  }

	lambda_A.multiplyby(gam_A(p));
	for(int r = 1; r<=num_superthreshold; r++)
	  {	    
	    for(int t=p+1; t <= ntpts; t++)
	      {
		float sume = 0;
		for(int e=1; e<=nrealevs; e++)
		  {
		    sume += designmatrix(e,t-p)*m_Beta[r-1](e);
		  }

		float wtp = Y(t-p,r)-sume;
		lambda_A.addto(r,r,gam_e(r)*Sqr(wtp));

		float sumq=0.0;
		for(int q=1; q<=ntar; q++)
		  if(q!=p && t>q)
		    {
		      sume = 0;
		      for(int e=1; e<=nrealevs; e++)
			{		    
			  sume += designmatrix(e,t-q)*m_Beta[r-1](e);
			}
		      sumq += m_A[q-1](r)*(Y(t-q,r)-sume);		      
		    }

		sume = 0;
		for(int e=1; e<=nrealevs; e++)
		  {		    
		    sume += designmatrix(e,t)*m_Beta[r-1](e);
		  }

		beta(r) += gam_e(r)*wtp*(Y(t,r)-sumq-sume);
	      }
	  }

	if(FilmbabeOptions::getInstance().verbose.value())
	  OUT("Compute m_A");
	
	solveforx(lambda_A, beta, m_A[p-1]);

	if(FilmbabeOptions::getInstance().verbose.value())
	  OUT(m_A[p-1](1));

	if(FilmbabeOptions::getInstance().verbose.value())
	  OUT("Compute trace");

	float trace_new = 0;

	if(FilmbabeOptions::getInstance().ntracesamps.value() > 0)
	  {
	    float trace_tol = 0.001;
	    trace_new = solvefortracex(lambda_A, D, ilambdaDA[p-1], FilmbabeOptions::getInstance().ntracesamps.value(), trace_tol);
	  }
	else
	  {
	    for(int r=1; r <= num_superthreshold; r++)
	      {
		trace_new += 1.0/lambda_A(r,r)*D(r,r);		
	      }	    
	  }

	trace_ilambdaDA(p) = trace_new;
	    
	if(FilmbabeOptions::getInstance().verbose.value())
	  OUT(trace_ilambdaDA(p));


      }
  }

  void Filmbabe_Vb_Flobs::update_phiA()
  {
    Tracer_Plus trace("Filmbabe_Vb_Flobs::update_phiA");

    // q(phiA)~Ga(b_A,c_A);

    float c_A = num_superthreshold/2.0;
    if(FilmbabeOptions::getInstance().verbose.value())
      OUT("Update phiA");
    for(int p=1; p<=ntar; p++)
      {
	float b_A = 1.0/(0.5*(quadratic(m_A[p-1],D) + trace_ilambdaDA(p)));
	
	gam_A(p) = exp(log(b_A)+MISCMATHS::lgam(c_A+1)-MISCMATHS::lgam(c_A));
	
	if(gam_A(p) > 1e6) gam_A(p) = 1e6;
	if(gam_A(p) < 1e-6) gam_A(p) = 1e-6;		
      }

    if(ntar>0)
      gamAhist.push_back(gam_A(1));
  }

  void Filmbabe_Vb_Flobs::update_phie()
  {
    Tracer_Plus trace("Filmbabe_Vb_Flobs::update_phie");

    float c_e = (ntpts-1)/2.0;
    if(FilmbabeOptions::getInstance().verbose.value())
      OUT("Update phie");

    for(int r=1; r <= num_superthreshold; r++)
      {
	float sum = 0.0;
	float sum2 = 0.0;

	for(int t=1; t <= ntpts; t++)
	  {
	    float sumq=0.0;
	    for(int q=1; q<=ntar; q++)
	      if(t>q)
		{
		  float sume = 0;
		  for(int e=1; e<=nrealevs; e++)
		    {	
		      sume += designmatrix(e,t-q)*m_Beta[r-1](e);
		    }
		  
		  sumq += m_A[q-1](r)*(Y(t-q,r)-sume);
		  //sum2 += diag_ilambdaA[q-1](r)*Sqr(Y(t-q,r)-sume);
		}

	    float sume = 0;
	    for(int e=1; e<=nrealevs; e++)
	      {		  
		sume += designmatrix(e,t)*m_Beta[r-1](e);
	      }

	    sum += Sqr(Y(t,r)-sumq-sume);
	  }

	sum += trace_ilambdaXXt(r);
	sum += sum2;

// 	if(sum==0)
// 	  {
// 	    OUT(sum);
	
// 	    const Voxel& vox = voxels[r-1]; 
	    
// 	    ColumnVector tmp = data.voxelts(vox.x,vox.y,vox.z);

// 	    write_ascii_matrix(tmp,"tmp");
// 	    OUT(vox.x);
// 	    OUT(vox.y);
// 	    OUT(vox.z);
// 	    OUT(mask(vox.x,vox.y,vox.z));
// 	  }

	float b_e = 1.0/(0.5*(sum));
	
	gam_e(r) = exp(log(b_e)+MISCMATHS::lgam(c_e+1)-MISCMATHS::lgam(c_e));
      }
    if(FilmbabeOptions::getInstance().verbose.value())
      {
	OUT(gam_e(1));
	OUT(gam_e(num_superthreshold));
      }
    
  }
  
  void Filmbabe_Vb_Flobs::save() 
    {
      Tracer_Plus trace("Filmbabe_Vb_Flobs::save");

      for(int p=1; p<=ntar; p++)
	{
          volume<float> A_mean(xsize,ysize,zsize);
	  A_mean = 0.0;
	  
	  for(int z = 0; z < data.zsize(); z++)    
	    for(int y = 0; y < data.ysize(); y++)	
	      for(int x = 0; x < data.xsize(); x++)
		if(mask(x,y,z))
		  {
		    A_mean(x,y,z) = m_A[p-1](indices(x,y,z));
		  }

	  copybasicproperties(data[0],A_mean);
	  save_volume(A_mean, LogSingleton::getInstance().appendDir("A_mean_"+num2str(p)));
	}

      for(int e=1; e<=nrealevs; e++)
	{
          volume<float> Beta_mean(xsize,ysize,zsize);
	  Beta_mean = 0.0;
	  
	  for(int z = 0; z < data.zsize(); z++)    
	    for(int y = 0; y < data.ysize(); y++)	
	      for(int x = 0; x < data.xsize(); x++)
		if(mask(x,y,z))
		  {
		    Beta_mean(x,y,z) = m_Beta[indices(x,y,z)-1](e);		    
		  }

	  copybasicproperties(data[0],Beta_mean);
	  save_volume(Beta_mean, LogSingleton::getInstance().appendDir("pe"+num2str(e)));
	  
	}

      for(int e=1; e<=nflobsevs; e++)
	{
	  volume<float> R_mean(xsize,ysize,zsize);
	  R_mean = 0.0;
	  
	  for(int z = 0; z < data.zsize(); z++)    
	    for(int y = 0; y < data.ysize(); y++)	
	      for(int x = 0; x < data.xsize(); x++)
		if(mask(x,y,z))
		  {
		    R_mean(x,y,z) = m_Beta[indices(x,y,z)-1](nrealevs+e);
		  }

	  copybasicproperties(data[0],R_mean);
	  save_volume(R_mean, LogSingleton::getInstance().appendDir("R_mean_e"+num2str(e)));
      
	}

      volume4D<float> all_Beta_cov(xsize,ysize,zsize,(nrealevs+nflobsevs)*(nrealevs+nflobsevs));
      all_Beta_cov = 0.0;

      volume4D<float> Beta_cov(xsize,ysize,zsize,nrealevs*nrealevs);
      Beta_cov = 0.0;

      volume4D<float> R_cov(xsize,ysize,zsize,nflobsevs*nflobsevs);
      R_cov = 0.0;

      volume<float> sigmasq(xsize,ysize,zsize);
      sigmasq=0.0;

      for(int z = 0; z < data.zsize(); z++)    
	for(int y = 0; y < data.ysize(); y++)	
	  for(int x = 0; x < data.xsize(); x++)
	    if(mask(x,y,z))
	      {
		sigmasq(x,y,z) = 1;

		for(int e=1; e<=nrealevs+nflobsevs; e++)
		  {
		    for(int e2=1; e2<=nrealevs+nflobsevs; e2++)
		      {
			all_Beta_cov(x,y,z,e+(e2-1)*(nrealevs+nflobsevs)-1) = ilambda_Beta[indices(x,y,z)-1](e,e2);
		      }

		  }

		for(int e=1; e<=nrealevs; e++)
		  {
		    for(int e2=1; e2<=nrealevs; e2++)
		      {
			Beta_cov(x,y,z,e+(e2-1)*nrealevs-1) = ilambda_Beta[indices(x,y,z)-1](e,e2);
		      }

		  }

		for(int e=1; e<=nflobsevs; e++)
		  {
		    for(int e2=1; e2<=nflobsevs; e2++)
		      {
			R_cov(x,y,z,e+(e2-1)*nflobsevs-1) = ilambda_Beta[indices(x,y,z)-1](nrealevs+e,nrealevs+e2);
		      }
		  }			
	      }

      copybasicproperties(data[0],all_Beta_cov[0]);
      save_volume4D(all_Beta_cov, LogSingleton::getInstance().appendDir("all_Beta_cov"));
      copybasicproperties(data[0],Beta_cov[0]);
      save_volume4D(Beta_cov, LogSingleton::getInstance().appendDir("corrections"));
      copybasicproperties(data[0],R_cov[0]);
      save_volume4D(R_cov, LogSingleton::getInstance().appendDir("R_cov"));
      copybasicproperties(data[0],sigmasq);
      save_volume(sigmasq, LogSingleton::getInstance().appendDir("sigmasquareds"));

      Matrix a(1,1);
      a =1000;
      write_ascii_matrix(a,LogSingleton::getInstance().appendDir("dof"));

      if(gamAhist.size()>0)
	{
	  miscplot newplot;
	  newplot.add_xlabel("Iterations");    
	  newplot.set_xysize(610,300);
	  newplot.timeseries(vector2ColumnVector(gamAhist).t(), LogSingleton::getInstance().appendDir("gamAhist"), "AR(1) MRF Precision", 0,400,3,0,false);
	}
      
      write_vector(gamAhist, LogSingleton::getInstance().appendDir("gamAhist"));
    }	    

}

