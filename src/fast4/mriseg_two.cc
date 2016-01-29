/*  FAST4 - FMRIB's Automated Segmentation Tool v4

    John Vickers, Mark Jenkinson, Steve Smith and Matthew Webster
    FMRIB Image Analysis Group

    Copyright (C) 2005-2007 University of Oxford  */

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

#include "mriseg_two.h" 
#include "newimage/newimageall.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
using namespace NEWIMAGE;
using namespace std;
ZMRISegmentation::ZMRISegmentation(){}

const float PI=3.14159265;
/////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////Tanaka Start/////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

void ZMRISegmentation::TanakaCreate(const NEWIMAGE::volume<float>& image, float fbeta, int nclasses, float nblowpass, bool bbias, int biterationspve, float mixeltypeMRF, int nbiter, int initinitfixed, int winitfixed, int bapused, float Hyp, bool verb,string mansegfle,int typeoffile)
{
  mansegfile=mansegfle; 
  bapusedflag=bapused;
  nClasses=nclasses;
  imagetype=typeoffile;
  m_Mricopy=volume<float>(image);
  m_Mri=volume<float>(image);
  m_nWidth=image.xsize();
  m_nHeight=image.ysize();
  m_nDepth=image.zsize();
  m_nxdim=image.xdim();
  m_nydim=image.ydim();
  m_nzdim=image.zdim();
  m_post=volume4D<float>(m_nWidth, m_nHeight, m_nDepth, nClasses+1);
  m_post.copyproperties(m_Mri);
  m_prob=volume4D<float>(m_nWidth, m_nHeight, m_nDepth, nClasses+1);
  m_prob.copyproperties(m_Mri);
  m_BiasField=volume<float>(m_nWidth, m_nHeight, m_nDepth);
  m_BiasField.copyproperties(m_Mri);
  m_mask=image;
  m_mask.binarise(0,m_mask.max()+1,exclusive);  
  nvoxel=(long)m_mask.sum();
  m_Segment=volume<int>(m_nWidth, m_nHeight, m_nDepth);
  copybasicproperties(m_mask,m_Segment);
  members=volume4D<float>(m_nWidth, m_nHeight, m_nDepth, nClasses);
  iterationspve=biterationspve;
  if(iterationspve>0)
    {
      m_pveSegment=volume<int>(m_nWidth, m_nHeight, m_nDepth);
      m_pveSegment.copyproperties(m_Segment);
    }
  pveBmixeltype=mixeltypeMRF;
  inititerations=winitfixed;
  initfixed=initinitfixed;
  m_nbIter=nbiter;
  biasfieldremoval=bbias;
  beta=fbeta;
  m_nbLowpass=nblowpass;
  if(bapusedflag>0) {
    talpriors=volume4D<float>(m_nWidth, m_nHeight, m_nDepth, nClasses+1);
    talpriors.copyproperties(m_post);
  }
  Hyper=Hyp;
  verboseusage=verb;
}

int ZMRISegmentation::TanakaMain(NEWIMAGE::volume<float>& pcsf, NEWIMAGE::volume<float>& pgm, NEWIMAGE::volume<float>& pwm)
{
  for(int z=0;z<m_nDepth;z++)
    for(int y=0;y<m_nHeight;y++)
      for(int x=0;x<m_nWidth;x++)
	if(m_Mricopy.value(x, y, z)>0.0) 
	  m_Mricopy.value(x, y, z)=log(m_Mricopy.value(x, y, z) + 1.0);
	else 
	  m_Mricopy.value(x, y, z)=0.0f;
 
  m_Mri=volume<float>(m_Mricopy);
  Dimensions();

  long seed=-1;
  srand(seed);

  rhs.resize(100);
  volumequant.resize(nClasses+1);
  m_mean.resize(nClasses+1,0);
  m_variance.resize(nClasses+1,0);
  weight.resize(nClasses,1);
  if(bapusedflag==0)
  {
    try{Initialise();}
    catch (kmeansexception& km){cout << km.what() << endl ; return 1;}
  }
  else
    {
      InitSimple(pcsf, pgm, pwm);
      pcsf.destroy();
      pgm.destroy();
      pwm.destroy();
    }
  if(Hyper<0.0)
    TanakaPriorHyper();
  InitKernel(m_nbLowpass);
  // first loop to remove bias field
  BiasRemoval();

  for (int n=1; n<=nClasses; n++)
    if(isnan(m_mean[n]) || isnan(m_variance[n])  ) cout << "MeaNsK variance nan" << endl;

  for(int iter=0;iter<m_nbIter;iter++)
  {
    if(verboseusage)
      cout<<"Tanaka Iteration "<<iter<<" bias field "<<m_nbIter<<endl;
    TanakaIterations();
    BiasRemoval();
    MeansVariances(nClasses, m_post);
  }
  // second loop to estimate hyperparameter
  for(int iter=0;iter<initfixed;iter++)
  {
    if(verboseusage)
      cout<<"Tanaka Iteration "<<iter<<" hyperparameter "<<initfixed<<endl;
    TanakaIterations();
    if(Hyper>=0.0)
      beta=Hyper;
    else
      TanakaHyper();
    if(verboseusage)
      cout<<" BETA "<<beta<<endl;
    MeansVariances(nClasses, m_post);
  }
  qsort();
  Classification();
  if(iterationspve>0) {
    if(verboseusage)
      cout<< " Starting Partial Volume Estimation \n";
    UpdateMembers(m_post);
    for(int z=0;z<m_nDepth;z++)
      for(int y=0;y<m_nHeight;y++)
	for(int x=0;x<m_nWidth;x++)
	  if(m_mask(x, y, z)==1)
	    m_Mri.value(x, y, z)=exp(m_Mri.value(x, y, z));
    MeansVariances(nClasses, m_post);

    PVClassificationStep();
    PVestimation();
    calculateVolumes(m_pve);
    pveClassification();
  }
  else 
    calculateVolumes(m_post);
  if(verboseusage)
    printVolumeTotals();
  takeexpo();
  return 0;
}


void ZMRISegmentation::BiasRemoval()
{ 
  m_BiasField=0;

  if(biasfieldremoval) {
    volume<float> p_meaninvcov(m_BiasField),p_resmean(m_BiasField);
    for(int z=0;z<m_nDepth;z++)
      for(int y=0;y<m_nHeight;y++)
	for(int x=0;x<m_nWidth;x++)
	  if(m_mask.value(x, y, z)==1)
	    for(int c=1;c<nClasses+1;c++) 
	    {
	      float tempf=m_post(x, y, z, c)/m_variance[c];
	      p_meaninvcov.value(x, y, z)+=tempf;
	      p_resmean.value(x, y, z)+=tempf*(m_Mricopy.value(x, y, z)-m_mean[c]);
	    }
    p_meaninvcov=Convolve(p_meaninvcov);
    p_resmean=Convolve(p_resmean);
    m_BiasField=divide(p_resmean,p_meaninvcov,m_mask);
    m_BiasField-=m_BiasField.mean(m_mask); // subtract off within-mask mean of m_BiasField
    m_BiasField*=m_mask;  //Zero out regions outside of mask ( is this needed?)
    m_Mri=m_Mricopy-m_BiasField;
  }
}


void ZMRISegmentation::Classification()
{
  for(int z=0;z<m_nDepth;z++)
    for(int y=0;y<m_nHeight;y++)
      for(int x=0;x<m_nWidth;x++)
	Classification(x,y,z);
}


void ZMRISegmentation::Classification(int x, int y, int z)
{
  float max=-1e-10;
  m_Segment(x, y, z)=0;
  if(m_mask(x, y, z)==1)
    for(int c=1;c<=nClasses;c++)
       if(m_post(x, y, z, c)>max) {
	 max=m_post(x, y, z, c);
	 m_Segment(x, y, z)=c;
       }
}


NEWIMAGE::volume<float> ZMRISegmentation::Convolve(NEWIMAGE::volume<float>& resfieldimage)
{
  return convolve_separable(resfieldimage, kernelx, kernely, kernelz);
}


void ZMRISegmentation::Dimensions()
{
  amz=1.0/m_nzdim;
  amy=1.0/m_nydim;
  amzy=1.0/sqrt(m_nzdim*m_nzdim+m_nydim*m_nydim);
  amx=1.0/m_nxdim;
  amzx=1.0/sqrt(m_nzdim*m_nzdim+m_nxdim*m_nxdim);
  amxy=1.0/sqrt(m_nxdim*m_nxdim+m_nydim*m_nydim);
}


void ZMRISegmentation::Initialise()
{
  WeightedKMeans();
  Initclass();
  Classification();
}


void ZMRISegmentation::pveClassification(int x, int y, int z)
{
  float max=-1e-10;
  m_pveSegment(x, y, z)=0;
  if(m_mask(x, y, z)==1)
    for(int c=1;c<=nClasses;c++)
      if(m_pve(x, y, z, c)>max) {
	max=m_pve(x, y, z, c);
	m_pveSegment(x, y, z)=c;
      }

}


void ZMRISegmentation::pveClassification()
{
  for(int x=0;x<m_nWidth;x++)
    for(int y=0;y<m_nHeight;y++)
      for(int z=0;z<m_nDepth;z++)
	pveClassification(x, y, z);
}


float ZMRISegmentation::logGaussian(float y, float mu, float sigmasq)
{
  float scaletemp=1.0*log(sqrt(2.0*PI*(sigmasq)));
  float gausstemp=(y-(mu));
  gausstemp=0.5*gausstemp*gausstemp/(sigmasq);
  return scaletemp+gausstemp;
}


float ZMRISegmentation::MRFWeightsTotal()
{
double total=0.0f; //Internally a double to avoid truncation when adding small to large (shouldn't happen anyway with this loop style)
  for(int z=0;z<m_nDepth;z++)
    for(int y=0;y<m_nHeight;y++)
      for(int x=0;x<m_nWidth;x++)
	if(m_mask.value(x, y, z)==1)
	  for(int c=1;c<=nClasses;c++)
	    total+=MRFWeightsInner(x,y,z,c)*m_post(x, y, z, c);
  return (float)total;   
}


double ZMRISegmentation::MRFWeightsInner(const int x, const int y, const int z,const int c)
{
double inner=0.0f;
  for(int n=-1;n<=1;n++)
    for(int m=-1;m<=1;m++)
      for(int l=-1;l<=1;l++)
	if(m_mask(x+l, y+m, z+n)==1)
	  inner+=MRFWeightsAM(l,m,n)*m_post.value(x+l, y+m, z+n, c);
  return inner;
}


double ZMRISegmentation::MRFWeightsAM(const int l, const int m, const int n)
{
double am=0.0f;	
   if((l==0))
   {
      if(((m==0)&&(n!=0)))am=amz;
      if(((m!=0)&&(n==0)))am=amy;
      if(((m!=0)&&(n!=0)))am=amzy;
   }    
   if((m==0))
   {
     if(((l!=0)&&(n==0)))am=amx;
     if(((l!=0)&&(n!=0)))am=amzx;
   }
   if((n==0))
   {
      if(((m!=0)&&(l!=0)))am=amxy;
   }
   return am;
}


void ZMRISegmentation::TanakaHyper()
{
  float total=MRFWeightsTotal();
  float betahtemp=0.0f;
  beta=betahtemp;
  float min=abs(total-rhs[0]);

  cout << " TOTAL " << total << " RHS ";
  for(int betahyp=0;betahyp<15;betahyp++)
    cout << rhs[betahyp] << " ";
  cout << endl;

  for(int betahyp=1;betahyp<15;betahyp++)
    {
      //float betahtemp=(float)(0.05*(float)betahyp);
      betahtemp+=(float)0.05;
      cerr <<  betahtemp << " " << betahtemp << endl;
      if(abs(total-rhs[betahyp])<min)
	{
	  beta=betahtemp;
	  min=abs(total-rhs[betahyp]);
	}
    }
}


void ZMRISegmentation::TanakaPriorHyper()
{
  float betahtemp=0.0f;
  for(int betahyp=0;betahyp<15;betahyp++)
    {
      for(int z=0;z<m_nDepth;z++)
	{
	  for(int y=0;y<m_nHeight;y++)
	    { 
	      for(int x=0;x<m_nWidth;x++)
		{
		  m_prob(x, y, z, 0)=0.0f;
		  float sum=0.0f;
		  for(int c=1;c<=nClasses;c++)
		    {
		      if(m_mask(x, y, z)==1)
			{
			  sum+=m_prob(x, y, z, c)=(float)(rand()/(float)(RAND_MAX));
			}
		      else
			m_prob(x, y, z, c)=0.0f;		    
		    }
		  for(int c=1;c<=nClasses;c++)
		    if(sum>0.0)
		      m_prob(x, y, z, c)/=sum;
		}
	    }
	}
      for(int iteration=0;iteration<5;iteration++)
	{
	  m_post=m_prob;
	  if(verboseusage)
	    cout << "Tanaka-inner-loop-iteration=" << iteration << " MRFWeightsTotal=" << MRFWeightsTotal() << " beta=" << betahtemp << endl;
	  for(int z=0;z<m_nDepth;z++)
	    {
	      for(int y=0;y<m_nHeight;y++)
		{
		  for(int x=0;x<m_nWidth;x++)
		    {
		      float sum=0.0f;
		      if(m_mask.value(x, y, z)==1)
			{
			  for(int c=1;c<=nClasses;c++)
			    {
			      float post=MRFWeightsInner(x,y,z,c);
			      if(bapusedflag<2)
				  sum+=m_prob.value(x, y, z, c)=exp(betahtemp*post);
			      else
				sum+=m_prob.value(x, y, z, c)=talpriors(x, y, z ,c)*exp(betahtemp*post);				
			    }
			  for(int c=1;c<=nClasses;c++)
			    {
			      if(sum>0.0f)
				m_prob.value(x, y, z, c)/=sum;
			      else
				m_prob.value(x, y, z, c)=0.0f;
			    }
			}
		    }
		}
	    }
	}
      m_post=m_prob;
      rhs[betahyp] = MRFWeightsTotal();
      betahtemp+=0.05;
    }
}


void ZMRISegmentation::TanakaIterations()
{  
  m_post=0;
  for(int z=0;z<m_nDepth;z++)
    for(int y=0;y<m_nHeight;y++)
      for(int x=0;x<m_nWidth;x++) 
	if(m_mask(x,y,z)==1) {
	  float sum=0.0f;
	  for(int c=1;c<=nClasses;c++) 
	    sum+=m_post(x, y, z, c)=(float)(rand()/(float)(RAND_MAX));	 		    
	  for(int c=1;c<=nClasses && sum>0;c++)
	    m_post(x, y, z, c)/=sum;
	}
  
  for(int iteration=0;iteration<5;iteration++) {
    for(int z=0;z<m_nDepth;z++) 
      for(int y=0;y<m_nHeight;y++)
	for(int x=0;x<m_nWidth;x++)
	  if(m_mask.value(x, y, z)==1) {
	    double sum=0.0f;       //Very important for sum to be a double for JV's loop ( if sum is a double there as well ) and mine to match results - numerically sensitive
	    for(int c=1;c<=nClasses;c++) {
	      m_post.value(x, y, z, c)=exp(beta*MRFWeightsInner(x,y,z,c)-logGaussian(m_Mri.value(x, y, z), m_mean[c], m_variance[c]));
              if(bapusedflag>=2) 
		m_post.value(x, y, z, c)*=talpriors(x, y, z, c);
	      sum+=m_post.value(x, y, z, c);
	    }
	    for(int c=1;c<=nClasses;c++)
	    {
	      if(sum>0.0f)
		m_post.value(x, y, z, c)/=sum;
	      else
		m_post.value(x, y, z, c)=0.0f;
	    }
	  }

      if(verboseusage)
	cout << "Tanaka-inner-loop-iteration=" << iteration << " MRFWeightsTotal=" << MRFWeightsTotal() << " beta=" << beta << endl;
  }
  if(verboseusage)
    {
      for(int c=1;c<=nClasses;c++)
	cout<<" CLASS "<<c<<" MEAN "<<exp(m_mean[c])<<" STDDEV "<<( exp(m_mean[c] + sqrt(m_variance[c]) ) - exp(m_mean[c]-sqrt(m_variance[c])) )/2;
      cout<<endl;
    }
}


 float ZMRISegmentation::PVEnergy(int x, int y, int z, float mu, float sigmasq)
{
  float temp=m_Mri(x, y, z)-mu;
  float cond=((temp*temp)/sigmasq +log(sigmasq))/2.0;
  return cond;
}


void ZMRISegmentation::PVClassificationStep()
{
  int mixnum=0;
  if(nClasses==3)
    mixnum=6;
  if(nClasses==2 || nClasses>3)
    mixnum=2*nClasses-1;
  PVprob.reinitialize(m_nWidth, m_nHeight, m_nDepth, mixnum);
  PVprob=0.0;
  for(int z=0;z<m_nDepth;z++)
    for(int y=0;y<m_nHeight;y++)
      for(int x=0;x<m_nWidth;x++)
	if(m_mask.value(x, y, z)==1) 
	  for(int type=0;type<mixnum;type++)
	  {
	    if(type<nClasses)
	      PVprob(x, y, z, type)=exp(-1.0*logGaussian(m_Mri(x, y, z), m_mean[type+1], m_variance[type+1]));
	    else if(nClasses==3)
	    {
		float mu=0.0f;
		float sigsq=0.0f;
		for(float delta=0.0;delta<=1.0;delta+=0.01)
		{
		  if(type==3)
		  {
		    mu=delta*m_mean[1]+(1-delta)*m_mean[2];
		    sigsq=delta*delta*m_variance[1]+(1-delta)*(1-delta)*m_variance[2];
		  }
		  if(type==4)
		  {
		    mu=delta*m_mean[1]+(1-delta)*m_mean[3];
		    sigsq=delta*delta*m_variance[1]+(1-delta)*(1-delta)*m_variance[3];
		  }
		  if(type==5)
		  {
		    mu=delta*m_mean[2]+(1-delta)*m_mean[3];
		    sigsq=delta*delta*m_variance[2]+(1-delta)*(1-delta)*m_variance[3];
		  }
		  PVprob(x, y, z, type)+=exp(-1.0*logGaussian(m_Mri(x, y, z), mu, sigsq))*0.01;
		}
	    }
	    else if(nClasses==2 || nClasses>3)
	    {
		float mu=0.0f;
		float sigsq=0.0f;
		for(float delta=0.0;delta<=1.0;delta+=0.01)
		{
		  mu=delta*m_mean[type-nClasses+1]+(1-delta)*m_mean[type-nClasses+2];
		  sigsq=delta*delta*m_variance[type-nClasses+1]+(1-delta)*(1-delta)*m_variance[type-nClasses+2];
		  PVprob(x, y, z, type)+=exp(-1.0*logGaussian(m_Mri(x, y, z), mu, sigsq))*0.01;
		}
	    }
	  } 
  ICMPV();
}


void ZMRISegmentation::ICMPV()
{
  hardPV=volume<int>(m_nWidth, m_nHeight, m_nDepth);
  hardPV.copyproperties(m_Segment);
  int mixnum=0;
  if(nClasses==3)
    mixnum=6;
  if(nClasses==2)
    mixnum=3;
  if(nClasses>3)
    mixnum=2*nClasses-1;
  float* clique=new float[mixnum];
  for(int z=0;z<m_nDepth;z++)
    for(int y=0;y<m_nHeight;y++)
      for(int x=0;x<m_nWidth;x++)
	if(m_mask.value(x, y, z)==1) {
	  float max=-1.0;
	  for(int type=0;type<mixnum;type++)
	    if(max<PVprob(x, y, z, type)) {
	      hardPV(x, y, z)=type;
	      max=PVprob(x, y, z, type);
	    }
	}

  for(int iter=0;iter<1;iter++)
    for(int z=0;z<m_nDepth;z++)
      for(int y=0;y<m_nHeight;y++)
	for(int x=0;x<m_nWidth;x++)
	  if(m_mask.value(x, y, z)==1)
	  {
	    for(int type=0;type<mixnum;type++)
	      clique[type]=0.0f;
	    for(int n=-1;n<=1;n++)
	      for(int m=-1;m<=1;m++)
		for(int l=-1;l<=1;l++)	
		  if(m_mask(x+l, y+m, z+n)==1)
		  {
		    float am=MRFWeightsAM(l, m, n);
		    if(nClasses==3)
		      for(int type=0;type<6;type++)
		      {	
			if(type==hardPV(x+l, y+m, z+n))
			  clique[type]+=am*2;
			else if((type==0)&&((hardPV(x+l, y+m, z+n)==3)||(hardPV(x+l, y+m, z+n)==4)))
			  clique[type]+=am;
			else if((type==1)&&((hardPV(x+l, y+m, z+n)==3)||(hardPV(x+l, y+m, z+n)==5)))
			  clique[type]+=am;
			else if((type==2)&&((hardPV(x+l, y+m, z+n)==4)||(hardPV(x+l, y+m, z+n)==5)))
			  clique[type]+=am;
			else if((type==3)&&((hardPV(x+l, y+m, z+n)==0)||(hardPV(x+l, y+m, z+n)==1)))
			  clique[type]+=am;
			else if((type==4)&&((hardPV(x+l, y+m, z+n)==0)||(hardPV(x+l, y+m, z+n)==2)))
			  clique[type]+=am;
			else if((type==5)&&((hardPV(x+l, y+m, z+n)==1)||(hardPV(x+l, y+m, z+n)==2)))
			  clique[type]+=am;
			else
			  clique[type]-=am;
			}
		    if(nClasses>3)
		      for(int type=0;type<mixnum;type++)
		      {	
			if(type==hardPV(x+l, y+m, z+n))
			  clique[type]+=am*2;
			else if((0<type)&&(type<nClasses-1)&&((hardPV(x+l, y+m, z+n)==(nClasses+type-1))||(hardPV(x+l, y+m, z+n)==(nClasses+type))))
			  clique[type]+=am;
			else if((type==0)&&((hardPV(x+l, y+m, z+n)==nClasses)))
			  clique[type]+=am;
			else if((type==(nClasses-1))&&((hardPV(x+l, y+m, z+n)==mixnum-1)))
			  clique[type]+=am;
			else if((type>(nClasses-1))&&((hardPV(x+l, y+m, z+n)==(type-nClasses))||(hardPV(x+l, y+m, z+n)==(type-nClasses+1))))
			  clique[type]+=am;
		        else 
			  clique[type]-=am;
		      }
		    if(nClasses==2)
		      for(int type=0;type<3;type++)
		      {
			if(type==hardPV(x+l, y+m, z+n))
      			  clique[type]+=am*2;
			else if((type==0)&&((hardPV(x+l, y+m, z+n)==2)))
			  clique[type]+=am;
			else if((type==1)&&((hardPV(x+l, y+m, z+n)==2)))
			  clique[type]+=am;
			else if((type==2)&&((hardPV(x+l, y+m, z+n)==0)||(hardPV(x+l, y+m, z+n)==1)))
			  clique[type]+=am;
			else 
			  clique[type]-=am;
		      }
		  }
	  
	    float max=-1;
	    for(int type=0;type<mixnum;type++)
	    {
	      float prob=PVprob(x, y, z, type)*exp(pveBmixeltype*clique[type]);
	      if(max<prob)
		{
		  hardPV(x, y, z)=type;
		  max=prob;
		}
	    }  
	  }
  delete[] clique;
}


void ZMRISegmentation::PVestimation()
{
  m_pve=volume4D<float>(m_nWidth, m_nHeight, m_nDepth, nClasses+1);
  m_pve.copyproperties(m_Mri);
  m_pve=0.0f;
  float mu=0.0f;
  float sigsq=0.0f;
  double step=(double)(1.0f/(double)(iterationspve));
  for(int z=0;z<m_nDepth;z++)
    {
      for(int y=0;y<m_nHeight;y++)
	{
	  for(int x=0;x<m_nWidth;x++)
	    { 
	     if(m_mask.value(x, y, z)>0)
	       {
		 if(nClasses==3)
		   {
		     float min=1.0e13;
		     if(hardPV(x, y, z)==0)
		       {
			 m_pve(x, y, z, 1)=1.0;
			 m_pve(x, y, z, 2)=0.0;
			 m_pve(x, y, z, 3)=0.0;
			 continue;
		       }
		     if(hardPV(x, y, z)==1)
		       {
			 m_pve(x, y, z, 1)=0.0;
			 m_pve(x, y, z, 2)=1.0;
			 m_pve(x, y, z, 3)=0.0;
			 continue;
		       }
		     if(hardPV(x, y, z)==2)
		       {
			 m_pve(x, y, z, 1)=0.0;
			 m_pve(x, y, z, 2)=0.0;
			 m_pve(x, y, z, 3)=1.0;
			 continue;
		       }
		     for(float delta=0.00;delta<=1.0;delta+=step)
		       {
			 if(hardPV(x, y, z)==3)
			   {
			     mu=delta*m_mean[1]+(1-delta)*m_mean[2];
			     sigsq=delta*delta*m_variance[1]+(1-delta)*(1-delta)*m_variance[2];
			   }
			 if(hardPV(x, y, z)==4)
			   {
			     mu=delta*m_mean[1]+(1-delta)*m_mean[3];
			     sigsq=delta*delta*m_variance[1]+(1-delta)*(1-delta)*m_variance[3];
			   }
			 if(hardPV(x, y, z)==5)
			   {
			     mu=delta*m_mean[2]+(1-delta)*m_mean[3];
			     sigsq=delta*delta*m_variance[2]+(1-delta)*(1-delta)*m_variance[3];
			   }
			 float val=PVEnergy(x, y, z, mu, sigsq);
			 if(min>val)
			   {
			     if(hardPV(x, y, z)==3)
			       {
				 m_pve(x, y, z, 1)=delta;
				 m_pve(x, y, z, 2)=1.0-delta;
				 m_pve(x, y, z, 3)=0.0;
			       }
			     if(hardPV(x, y, z)==4)
			       {
				 m_pve(x, y, z, 1)=delta;
				 m_pve(x, y, z, 2)=0.0;
				 m_pve(x, y, z, 3)=1.0-delta;
			       }
			     if(hardPV(x, y, z)==5)
			       {
				 m_pve(x, y, z, 1)=0.0;
				 m_pve(x, y, z, 2)=delta;
				 m_pve(x, y, z, 3)=1.0-delta;
			       }
			     min=val;
			   }
		       }
		   }
		 if(nClasses>3)
		   {
		     float min=1.0e13;
		     if(hardPV(x, y, z)<nClasses)
		       {
			 for(int c=0;c<nClasses;c++)
			   {
			     if(hardPV(x, y, z)!=c)
			       m_pve(x, y, z, c+1)=0.0;
			   }
			 m_pve(x, y, z, hardPV(x, y, z)+1)=1.0f;
			 continue;
		       }
		     if(hardPV(x, y, z)>=nClasses)
		       {
			 for(float delta=0.00;delta<=1.0;delta+=step)
			   {
			     mu=delta*m_mean[hardPV(x, y, z)-nClasses+1]+(1-delta)*m_mean[hardPV(x, y, z)-nClasses+2];
			     sigsq=delta*delta*m_variance[hardPV(x, y, z)-nClasses+1]+(1-delta)*(1-delta)*m_variance[hardPV(x, y, z)-nClasses+2];
			     float val=PVEnergy(x, y, z, mu, sigsq);
			     if(min>val)
			       {
			     
				 for(int c=0;c<nClasses;c++)
				   {
				     m_pve(x, y, z, c+1)=0.0f;
				   }
				 m_pve(x, y, z, hardPV(x, y, z)-nClasses+1)=delta;
				 m_pve(x, y, z, hardPV(x, y, z)-nClasses+2)=1.0-delta;			     
				 min=val;
			       }
			   }
		       }
		   }

		 if(nClasses==2)
		   {
		     float min=1.0e13;
		     if(hardPV(x, y, z)==0)
		       {
			 m_pve(x, y, z, 1)=1.0;
			 m_pve(x, y, z, 2)=0.0;
			 continue;
		       }
		     if(hardPV(x, y, z)==1)
		       {
			 m_pve(x, y, z, 1)=0.0;
			 m_pve(x, y, z, 2)=1.0;
			 continue;
		       }
		     for(float delta=0.00;delta<=1.0;delta+=step)
		       {
			 if(hardPV(x, y, z)==2)
			   {
			     mu=delta*m_mean[1]+(1-delta)*m_mean[2];
			     sigsq=delta*delta*m_variance[1]+(1-delta)*(1-delta)*m_variance[2];
			   }
		    
			 float val=PVEnergy(x, y, z, mu, sigsq);
			 if(min>val)
			   {
			     if(hardPV(x, y, z)==2)
			       {
				 m_pve(x, y, z, 0)=delta;
				 m_pve(x, y, z, 1)=1.0f-delta;
			       }
			     min=val;
			   }
		       }
		   }
	       }
	    }
	}
    }
}

void ZMRISegmentation::takeexpo()
{
  for(int z=0;z<m_nDepth;z++)
    for(int y=0;y<m_nHeight;y++)
      for(int x=0;x<m_nWidth;x++)
	m_BiasField(x, y, z)=exp(-m_BiasField(x, y, z));
}

void ZMRISegmentation::InitSimple(const NEWIMAGE::volume<float>& pcsf, const NEWIMAGE::volume<float>& pgm, const NEWIMAGE::volume<float>& pwm)
{
  if(verboseusage)
    cout<<"Beginning prior-based initialisation"<<endl;
  m_post=0;
  m_prob=0;
  talpriors=0;
  if(nClasses==3)
  {
    for(int z=0;z<m_nDepth;z++)
      for(int y=0;y<m_nHeight;y++)
	for(int x=0;x<m_nWidth;x++)
	  if(m_mask.value(x, y, z)==1)
	  {
	    float norm2=pcsf.value(x, y, z)+pgm.value(x, y, z)+pwm.value(x, y, z);
	    if(norm2>0.0)
	    {
	      m_post.value(x, y, z, 1)=m_prob.value(x, y, z, 1)=talpriors.value(x, y, z, 1)=pcsf.value(x, y, z)/norm2;
	      m_post.value(x, y, z, 2)=m_prob.value(x, y, z, 2)=talpriors.value(x, y, z, 2)=pgm.value(x, y, z)/norm2;
	      m_post.value(x, y, z, 3)=m_prob.value(x, y, z, 3)=talpriors.value(x, y, z, 3)=pwm.value(x, y, z)/norm2;
	    }
	    else for(int c=1;c<=nClasses;c++)
	      m_post.value(x, y, z, c)=m_prob.value(x, y, z, c)=talpriors.value(x, y, z, c)=1.0/3.0f;	      
	  }
    Classification();
  }
 
  if(nClasses==2)
  {
    for(int z=0;z<m_nDepth;z++)
      for(int y=0;y<m_nHeight;y++)
	for(int x=0;x<m_nWidth;x++)
	  if(m_mask.value(x, y, z)==1)
	    {
	      float norm2=pcsf.value(x, y, z)+pgm.value(x, y, z)+pwm.value(x, y, z);
	      if(norm2>0.0)
	      {
		m_post.value(x, y, z, 1)=m_prob.value(x, y, z, 1)=talpriors.value(x, y, z, 1)=pcsf.value(x, y, z)/norm2;
		m_post.value(x, y, z, 2)=m_prob.value(x, y, z, 2)=talpriors.value(x, y, z, 2)=(pgm.value(x, y, z)+pwm.value(x, y, z))/norm2;
	      }
	      else for(int c=1;c<=nClasses;c++)
		m_post.value(x, y, z, c)=m_prob.value(x, y, z, c)=talpriors.value(x, y, z, c)=1.0/2.0f;
	    }
    Classification();
  }
  MeansVariances(nClasses, m_post);
}


void ZMRISegmentation::calculateVolumes(const NEWIMAGE::volume4D<float>& probs)
{
  for(int c=1;c<nClasses+1;c++)
  {
    volumequant[c]=probs[c].sum(m_mask);
    volumequant[c]*=m_nxdim*m_nydim*m_nzdim;
  }
}

void ZMRISegmentation::printVolumeTotals()
{
  double tot(0);
  for(int c=1;c<nClasses+1;c++)
  {
    tot+=volumequant[c];
    cout<<"\n tissue "<<c<<" " << volumequant[c];
  }
  cout<<"\n total tissue "<<tot<<"\n";
}


void ZMRISegmentation::InitKernel(float kernalsize)
{
  float sigma=0.51*m_nbLowpass/m_nxdim;
  int radius=2*(int)sigma;
  kernelx=gaussian_kernel1D(sigma, radius);
  sigma=0.51*m_nbLowpass/m_nydim;
  radius=2*(int)sigma;
  kernely=gaussian_kernel1D(sigma, radius);
  sigma=0.51*m_nbLowpass/m_nzdim;
  radius=2*(int)sigma;
  kernelz=gaussian_kernel1D(sigma, radius);
}

void ZMRISegmentation::MeansVariances(int numberofclasses, NEWIMAGE::volume4D<float>& probability )
{
  volume4D<float> meanP=m_Mri*probability;
  volume4D<float> varianceP=m_Mri*meanP;
  for(int c=1;c<=numberofclasses;c++)
  {
    m_mean[c] = meanP[c].sum(m_mask)/probability[c].sum(m_mask);
    m_variance[c] = varianceP[c].sum(m_mask)/probability[c].sum(m_mask)-m_mean[c]*m_mean[c];
  }
}

void ZMRISegmentation::Initclass()
{
  m_prob=0;
  for(int z=0;z<m_nDepth;z++)
    for(int y=0;y<m_nHeight;y++)
      for(int x=0;x<m_nWidth;x++)
	if(m_mask.value(x, y, z)==1)
        {
	  float tot=0.0f;
	  for(int c=1;c<nClasses+1;c++)
	  {
	    m_prob.value(x, y, z, c)=logGaussian(m_Mri.value(x, y, z), m_mean[c], m_variance[c]);
	    tot+=m_prob(x, y, z, c)=max(exp(-m_prob(x, y, z, c)),0.0f); // Need max? Can prob be less than zero?
	  }
	  for(int c=1;c<nClasses+1 && tot>0;c++)
	    m_prob.value(x, y, z, c)/=tot; 
	}
  m_post=m_prob;
}

void ZMRISegmentation::UpdateMembers(NEWIMAGE::volume4D<float>& probability)
{
  for(int x=0;x<m_nWidth;x++)
    for(int y=0;y<m_nHeight;y++)
      for(int z=0;z<m_nDepth;z++)
	if(m_mask(x, y, z)==1) {
	  float sum=0.0f;
	  for(int c=0;c<nClasses;c++)
	    sum+=members(x, y, z, c)=max(weight[c]*probability(x, y, z, c+1),0.0f); // Need max? Can weight*prob be less than zero?
	  for(int c=0;c<nClasses && sum>0;c++)
	    members(x, y, z, c)/=sum;
	}
}

void ZMRISegmentation::WeightedKMeans()
{ 
  vector<float> inputMeans;
  if (mansegfile!="") 
  {  
    ifstream inputfile(mansegfile.c_str());
    copy(istream_iterator<float> (inputfile),istream_iterator<float> (),back_inserter(inputMeans));
    inputfile.close();
  }

  m_post=m_prob=0.0f;

  float perc=1.0/((float)(nClasses+1.0));     
  for(int c=1;c<=nClasses;c++)
  {
    if ( (int)inputMeans.size()==nClasses ) m_mean[c]=log(inputMeans[c-1]);
    else m_mean[c]=m_Mricopy.percentile((float)(perc*c), m_mask);
    if (verboseusage) cout << c << " " << m_mean[c] << endl;
  }

  for(int c=1;c<nClasses+1;c++)
    for(int c2=1;c2<c;c2++)
      if(m_mean[c]==m_mean[c2]) throw kmeansexc;
  for(int z=0;z<m_nDepth;z++)
    for(int y=0;y<m_nHeight;y++)
      for(int x=0;x<m_nWidth;x++)
	if(m_mask(x, y, z)==1)
	{
	  int minclass=0;
	  float mindist=1e31;
	  for(int c=1;c<=nClasses;c++)
	  {
	    if((m_Mricopy.value(x, y, z)-m_mean[c])*(m_Mricopy.value(x, y, z)-m_mean[c])<mindist)
	    {
	      mindist=(m_Mricopy.value(x, y, z)-m_mean[c])*(m_Mricopy.value(x, y, z)-m_mean[c]);
	      minclass=c;
	    }
	  }
	  m_post(x, y, z, minclass)=1.0f;
	}

  for(int initfiter=0;initfiter<initfixed+inititerations;initfiter++)
  {
    UpdateMembers(m_post);
    MeansVariances(nClasses, m_post);
    if(verboseusage) 
      cout<<"KMeans Iteration "<<initfiter<<"\n";
    Initclass();
  }
}

void ZMRISegmentation::qsort()
{ // sorts images so that 0 is csf, 1 is grey, 2 is white  
  vector<float> meancopy(m_mean), varcopy(m_variance);
  volume4D<float> postcopy(m_post), probcopy(m_prob);

  sort(m_mean.begin()+1, m_mean.end());
  if(imagetype==2) // for a T2 image reverse the intensity order
    reverse(m_mean.begin()+1, m_mean.end()); 

  for (int n=1; n<=nClasses; n++)
    for (int m=1; m<=nClasses; m++)
      if(m_mean[n]==meancopy[m]) 
      {
	m_variance[n]=varcopy[m];
        m_post[n]=postcopy[m];
        m_prob[n]=probcopy[m];
      }
}
