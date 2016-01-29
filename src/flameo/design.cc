/*  design.cc

    Mark Woolrich, Tim Behrens - FMRIB Image Analysis Group

    Copyright (C) 2002 University of Oxford  */

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

#include "design.h"
#include "utils/log.h"
#include "miscmaths/miscmaths.h"
#include "newmat.h"
#include "utils/tracer_plus.h"
#include "gsoptions.h"
#include <string>
#include "newimage/newimageall.h"

using namespace Utilities;
using namespace MISCMATHS;
using namespace NEWMAT;
using namespace NEWIMAGE;

namespace Gs {

  ReturnMatrix Design::getdm(int x, int y, int z, int g) const {

    Tracer_Plus trace("Design::getdm x y z g");  
 
    Matrix pdm = getdm(x,y,z);	 
	
    Matrix zg(getntptsingroup(g),pdm.Ncols());
    zg=0;

    int t2=1;
    for(int t = 1; t <= ntpts; t++)
      {
	if(getgroup(t)==g)
	  {		
	    zg.Row(t2++) = pdm.Row(t);
	  }
      }

    // remove non-zero columns
    Matrix tmp=zg;
    int e2=1;
    for(int e = 1; e<= tmp.Ncols(); e++)
      {	
	if(abs(tmp.Column(e)).Sum()>0)
	  zg.Column(e2++) = tmp.Column(e);
      }

    zg=zg.Columns(1,e2-1);
	       	
    zg.Release();
    return zg;
  }

  // returns design matrix with any voxelwise zero evs removed
  ReturnMatrix Design::getdm(int x, int y, int z) const { 

    Tracer_Plus trace("Design::getdm x y z");  

    Matrix ret=dm;

    if(is_voxelwise_dm())
      {
	// insert voxelwise evs
	for(unsigned int i=0; i<voxelwise_ev_numbers.size(); i++)      
	  {	    
	    ColumnVector ev_i=voxelwise_evs[i].voxelts(x,y,z);
	    ret.Column(voxelwise_ev_numbers[i])=ev_i;
	  }
	
	// remove any zero evs
	ColumnVector zero_ev=zero_evs.voxelts(x,y,z);
	if(zero_ev.Sum()>0)
	  {
	    int index=1;
	    Matrix new_ret=ret;
	    for(int i=1; i<=dm.Ncols(); i++)   
	      {		
		if(zero_ev(i)==0)
		  {
		    new_ret.Column(index)=ret.Column(i);
		    index++;
		  }
	      }   
	    
	    ret=new_ret.Columns(1,index-1);
	  }
      }

    ret.Release(); 
    return ret; 
  }

  ReturnMatrix Design::remove_zeroev_pes(int x, int y, int z, const ColumnVector& petmp)
  {
    ColumnVector ret=petmp;
	
    ColumnVector zero_ev=zero_evs.voxelts(x,y,z);

    if(zero_ev.Sum()>0)
      {
	int index=1;
	for(int i=1; i<=nevs; i++)   
	  {		
	    if(zero_ev(i)==0)
	      {
		ret(index)=petmp(i);
		index++;
	      }
	  }   
    
	ret=ret.Rows(1,index-1);
      }

    ret.Release();
    return ret;
  }

  ReturnMatrix Design::remove_zeroev_covpes(int x, int y, int z, const SymmetricMatrix& covpetmp)
  {
    SymmetricMatrix ret=covpetmp;

    ColumnVector zero_ev=zero_evs.voxelts(x,y,z);

    if(zero_ev.Sum()>0)
      {
	int index=1;
	for(int i=1; i<=nevs; i++) 
	  {
	    int index2=i+1;
	    for(int j=i+1; j<=nevs; j++)   
	      {		
		if(zero_ev(i)==0 && zero_ev(j)==0)
		  {
		    ret(index,index2)=covpetmp(i,j);
		    ret(index2,index)=covpetmp(j,i);
		  }

		if(zero_ev(j)==0)
		  {		
		    index2++;
		  }
	      }

	    if(zero_ev(i)==0)
	      {
		ret(index,index)=covpetmp(i,i);
		index++;
	      }
	    else
	      {
		ret(index,index)=1e32; // set variance very high for this parameter as it is a zero ev
	      }
	  }   

	ret=ret.SymSubMatrix(1,index-1);
      }

    ret.Release();
    return ret;
  }

  ReturnMatrix Design::insert_zeroev_pemcmcsamples(int x, int y, int z, const Matrix& mcmcin)
  {
    Matrix ret=mcmcin;

    ColumnVector zero_ev=zero_evs.voxelts(x,y,z);

    if(zero_ev.Sum()>0)
      {
	ret.ReSize(nevs,mcmcin.Ncols());
	ret=0;
	int index=1;
	for(int i=1; i<=nevs; i++)   
	  {
	    if(zero_ev(i)==0)
	      {
		ret.Row(i)=mcmcin.Row(index);
		index++;
	      }
	  }   
      }
	
    ret.Release();
    return ret;
  }

  ReturnMatrix Design::insert_zeroev_pes(int x, int y, int z, const ColumnVector& petmp)
  {
    ColumnVector ret=petmp;

    ColumnVector zero_ev=zero_evs.voxelts(x,y,z);
	
    if(zero_ev.Sum()>0)
      {
	    
	ret.ReSize(nevs);
	ret=0;
	int index=1;
	for(int i=1; i<=nevs; i++)   
	  {		
	    if(zero_ev(i)==0)
	      {		    
		ret(i)=petmp(index);
		index++;
	      }
	  }   
      }
	
    ret.Release();
    return ret;
  }

  ReturnMatrix Design::insert_zeroev_covpes(int x, int y, int z, const SymmetricMatrix& covpetmp)
  {
    SymmetricMatrix ret=covpetmp;

    ColumnVector zero_ev=zero_evs.voxelts(x,y,z);

    if(zero_ev.Sum()>0) {
      ret.ReSize(nevs);
      ret=0;
      vector<int> realIndex(nevs+1,-1);
      int index(1);
      for(int i=1; i<=nevs; i++) 
	if (zero_ev(i) == 0) 
	  realIndex.at(i)=index++;

      for(int i=1; i<=nevs; i++) {
	for(int j=i; j<=nevs; j++) {
	  if ( realIndex.at(i) > 0 && realIndex.at(j) > 0 ) 
	    ret(i,j)=covpetmp(realIndex.at(i),realIndex.at(j));	  
	}
	if ( realIndex.at(i) ==- 1 )
	  ret(i,i) = 1e32; // set variance very high for this parameter as it is a zero ev
      }
    }    
    ret.Release();
    return ret;
  }

  ReturnMatrix Design::getcopedata(int x, int y, int z, int g){
    ColumnVector Yg(getntptsingroup(g));
    Yg = 0;
    int t2=1;
    for(int t = 1; t <= ntpts; t++)
      {
	if(getgroup(t)==g)
	  {
	    Yg(t2) = copedata(x,y,z,t-1);
	    t2++;
	  }
      }
      
    Yg.Release();
    return Yg;
  }

  ReturnMatrix Design::getvarcopedata(int x, int y, int z, int g){
    ColumnVector Sg(getntptsingroup(g));
    Sg = 0;
    int t2=1;
    for(int t = 1; t <= ntpts; t++)
      {
	if(getgroup(t)==g)
	  {
	    Sg(t2) = varcopedata(x,y,z,t-1);
	    t2++;
	  }
      }
      
    Sg.Release();
    return Sg;
  }

  void Design::setup(bool loadcontrasts)
  {
    Tracer_Plus trace("Design::setup");  
 
    // read data
    read_volume4D(copedata, GsOptions::getInstance().copefile.value());    
    //copedata.read(GsOptions::getInstance().copefile.value());
    
    // mask:
    read_volume(mask, GsOptions::getInstance().maskfile.value());

    if(GsOptions::getInstance().varcopefile.value() != string(""))
      read_volume4D(varcopedata, GsOptions::getInstance().varcopefile.value());
    else
      {
	// set to zero
	varcopedata=copedata;
	varcopedata=0;
      }

    if(GsOptions::getInstance().dofvarcopefile.value() != string(""))
      {
	read_volume4D(dofvarcopedata, GsOptions::getInstance().dofvarcopefile.value());
      }
    else
      {
	dofvarcopedata = varcopedata;
	dofvarcopedata = 0;
      }

    // create sum_dofvarcopedata
    sum_dofvarcopedata=dofvarcopedata[0];
    for(int i=1; i<=dofvarcopedata.maxt(); i++)
      {
	sum_dofvarcopedata+=dofvarcopedata[i];
      }

    dm = read_vest(GsOptions::getInstance().designfile.value());
    //write_ascii_matrix(dm,"dm"); 
    nevs = dm.Ncols();
    ntpts = dm.Nrows();

   // check design and data compatability
    if(getntpts() != copedata.tsize())
      {
	cout << "dm.getntpts()=" << getntpts() << endl;
	cout << "copedata.tsize()=" << copedata.tsize() << endl;
	throw Exception("Cope data and design have different numbers of time points");
      }
    if(getntpts() != varcopedata.tsize())
      {	
	cout << "dm.ntpts()=" << getntpts() << endl;
	cout << "copedata.tsize()=" << varcopedata.tsize() << endl;
	
	throw Exception("Varcope data and design have different numbers of time points");
      }    


    // read in variance groupings
    group_index = read_vest(GsOptions::getInstance().covsplitfile.value());
    nevs = dm.Ncols();
    ntpts = dm.Nrows();
    ngs = int(group_index.Maximum());

    if(nevs == 0 || ntpts == 0)
	throw Exception(string(GsOptions::getInstance().designfile.value() + " is an invalid design matrix file").data());

    if(ngs == 0 || group_index.Nrows() == 0)
      throw Exception(string(GsOptions::getInstance().covsplitfile.value() + " is an invalid covariance split file").data());

    if(ntpts != group_index.Nrows())
      {
	cerr << "design matrix has ntpts=" <<  ntpts << endl;
	cerr << "cov split file has ntpts=" << group_index.Nrows() << endl;
	throw Exception("The ntpts needs to be the same for both");
      }

    if(GsOptions::getInstance().debuglevel.value()==2)
      {
	cout << "******************************************" << endl
	     << "Data Read Complete" << endl << "******************************************"
	     << endl;
      }

    // fill group_index:
    covsplit.ReSize(ntpts,ngs);
    covsplit = 0;
    ntptsing.ReSize(ngs);
    ntptsing = 0;
    nevsing.ReSize(ngs);
    nevsing = 0;
    global_index.resize(ngs);
    index_in_group.resize(ngs);

    for(int t = 1; t <= ntpts; t++)
      {
	covsplit(t,int(group_index(t)))=1;
	ntptsing(int(group_index(t)))++;
	global_index[int(group_index(t))-1].push_back(t);	
      }


    for(int g =1; g<=ngs; g++)
      {

	index_in_group[g-1].resize(ntpts);

	for(int wt = 1; wt <= int(ntptsing(g)); wt++)
	  {	  
	    int gt=global_index[g-1][wt-1];

	    index_in_group[g-1][gt-1]=wt;
	  }
      }

      
    cout << "ntptsing=" << ntptsing << endl;
//     cout << "group_index=" << group_index << endl;

    // assign each EV to a group
    evs_group.ReSize(nevs);
    evs_group=0;
    for(int e=1; e <=nevs; e++)
      {	  
	for(int t = 1; t <= ntpts; t++)
	  {	
	    if(dm(t,e)!=0)
	      {
		if(evs_group(e)==0)
		  {
		    evs_group(e) = group_index(t);
		    nevsing(int(group_index(t)))++;
		  }
		else if(evs_group(e)!=group_index(t))
		  {
		    cerr << "nonseparable design matrix and variance groups" << endl;
		    throw Exception("nonseparable design matrix and variance groups");
		  }
	      }
	  }
      }

    OUT(evs_group);

    // load contrasts:
    if(loadcontrasts)
      {
	try
	  {
	    tcontrasts = read_vest(GsOptions::getInstance().tcontrastsfile.value());

	    numTcontrasts = tcontrasts.Nrows();	
	  }
	catch(Exception exp)
	  {
	    numTcontrasts = 0;
	  }

	// Check contrast matches design matrix
	if(numTcontrasts > 0 && dm.Ncols() != tcontrasts.Ncols())
	  { 
	    cerr << "Num tcontrast cols  = " << tcontrasts.Ncols() << ", design matrix cols = " << dm.Ncols() << endl;
	    throw Exception("size of design matrix does not match t contrasts\n");
	  }

	try
	  {
	    fcontrasts = read_vest(GsOptions::getInstance().fcontrastsfile.value());
	    numFcontrasts = fcontrasts.Nrows();
	  }
	catch(Exception exp)
	  {
	    numFcontrasts = 0;
		cout << "No f contrasts" << endl;
	  }

	if(numFcontrasts > 0)
	  {
	    // Check contrasts match
	    if(tcontrasts.Nrows() != fcontrasts.Ncols())
	      { 
		cerr << "tcontrasts.Nrows()  = " << tcontrasts.Nrows() << endl;
		cerr << "fcontrasts.Ncols()  = " << fcontrasts.Ncols() << endl;
		throw Exception("size of tcontrasts does not match fcontrasts\n");
	      }
	    
	    if(numFcontrasts > 0)
	      setupfcontrasts();
	  }    	
      }  

    // read in any voxelwise EVs
    voxelwise_ev_numbers=GsOptions::getInstance().voxelwise_ev_numbers.value();    
    vector<string> voxelwise_ev_filenames=GsOptions::getInstance().voxelwise_ev_filenames.value();

    if(voxelwise_ev_filenames.size() != voxelwise_ev_numbers.size())
      throw Exception("Number of filenames in voxelwise_ev_filenames command line option needs to be the same as the number of EV number in the voxelwise_ev_numbers command line option");
    
    voxelwise_evs.resize(voxelwise_ev_filenames.size());
    voxelwise_dm=false;
  
    for(unsigned int i=0; i<voxelwise_ev_filenames.size(); i++)      
      {
	if(voxelwise_ev_numbers[i]>nevs)
	  throw Exception("EV number in the voxelwise_ev_numbers command line option is invalid (it's greater than the number of design matrix EVs)");
	
	voxelwise_dm=true;
	//	voxelwise_evs[i].read(voxelwise_ev_filenames[i]);	
	read_volume4D(voxelwise_evs[i], voxelwise_ev_filenames[i]);
      }    
    
    // find if there are any zero voxelwise evs
    // reinforce mask by checking for zeros in lower-level varcopes
    zero_evs.reinitialize(mask.xsize(),mask.ysize(),mask.zsize(),nevs);
    zero_evs=0;
    
    bool global_contains_zero_varcope=false;

    for(int x = 0; x < mask.xsize(); x++)
      for(int y = 0; y < mask.ysize(); y++)
	for(int z = 0; z < mask.zsize(); z++)
	  {
	    if(mask(x,y,z))
	      {
		// reinforce mask by checking for zeros in lower-level varcopes
		if(GsOptions::getInstance().varcopefile.value() != string(""))
		  {
		    bool contains_zero_varcope=false;
		    for(int t = 1; t <= ntpts && !contains_zero_varcope; t++)
		      {	
			if(varcopedata(x,y,z,t-1)<=0)
			  contains_zero_varcope=true;
		      }	
		    if(contains_zero_varcope) 
		      {
			mask(x,y,z)=0; global_contains_zero_varcope=true;
		      }
		  }
		
		// find if there are any zero voxelwise evs
		if(mask(x,y,z))
		  {		
		    for(unsigned int i=0; i<voxelwise_ev_numbers.size(); i++)      
		      {	    
			ColumnVector ev_i=voxelwise_evs[i].voxelts(x,y,z);		     
			
			if(abs(ev_i).Sum()==0)
			{
			  zero_evs(x,y,z,voxelwise_ev_numbers[i]-1)=1;			  
			}
		      }
		  }
	      }
	  }
    
    if(global_contains_zero_varcope)
      {
	cout << endl << "WARNING: The passed in varcope file, "<< GsOptions::getInstance().varcopefile.value()<< ", contains voxels inside the mask with zero (or negative) values. These voxels will be excluded from the analysis." << endl;      
      }

    //save_volume4D(zero_evs, LogSingleton::getInstance().appendDir("zero_evs"));	

  }
  
  void Design::setupfcontrasts()
  {
    Tracer_Plus trace("Design::setupfcontrasts");
    
    fc.resize(numFcontrasts);
    //reduceddms.resize(numFcontrasts);

    for(int f=0; f < numFcontrasts; f++)
    {
      fc[f].ReSize(tcontrasts.Nrows(),nevs);
      fc[f] = 0;

      int count = 0;
    
      for(int c = 1; c <= tcontrasts.Nrows(); c++)
	{	
	  if(fcontrasts(f+1,c) == 1)
	    {	    
	      count++;
	      fc[f].Row(count) = tcontrasts.Row(c);	

	    }
	}
    
      fc[f] = fc[f].Rows(1,count);

//       const Matrix& tfc = fc[f]; 
      
//       reduceddms[f] = getdm()*(IdentityMatrix(nevs)-tfc*pinv(tfc));
      
//       Matrix tmp = zeros(2,2);
//       tmp(1,1) = 1;
//       reduceddms[f] = getdm()*(IdentityMatrix(nevs)-tmp);

//        OUT(f);
//        OUT(fc[f]);
//       OUT(fc[f].t()*fc[f].i()*fc[f].t());	
//       OUT(dm);
//        OUT(reduceddms[f]);


    }
  }

  
}

