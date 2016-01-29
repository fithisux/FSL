/*  FAST4 - FMRIB's Automated Segmentation Tool v4

    John Vickers, Mark Jenkinson and Steve Smith
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

#include "multi_mriseg_two.h" 
#include "newimage/newimageall.h"

using namespace NEWIMAGE;

using namespace std;


ZMRIMULTISegmentation::ZMRIMULTISegmentation(){}

const float PI=3.14159265;

void ZMRIMULTISegmentation::TanakaCreate(const NEWIMAGE::volume<float>* images,  int nclasses, bool b2d, int nbiter,float nblowpass, float fbeta, int bapused, bool pveboolean, int nchan, bool bbias, int initinitfixed, bool verb, int biterationspve, int winit, float fpveMRFmixeltyp, float Hyp,string mansegfle, int typeoffile)
{
  mansegfile=mansegfle; 
  verboseusage=verb;
  inititerations=winit;
  m_nbLowpass=nblowpass;
  numberofchannels=nchan;
  m_nWidth=images[0].xsize();
  m_nHeight=images[0].ysize();
  m_mean.ReSize(nclasses,numberofchannels);
  m_nDepth=images[0].zsize();
  m_nxdim=images[0].xdim();
  m_nydim=images[0].ydim();
  m_nzdim=images[0].zdim();
  m_co_variance=new Matrix[nclasses+1];
    for(int i=1;i<=nclasses;i++)
  m_co_variance[i].ReSize(numberofchannels,numberofchannels);
  m_inv_co_variance=new Matrix[nclasses+1];
    for(int i=1;i<=nclasses;i++)
  m_inv_co_variance[i].ReSize(numberofchannels,numberofchannels);

  m_Mricopy=new volume<float>[numberofchannels+1];
  imagetype=typeoffile;
  for(int i=1;i<=numberofchannels;i++)
    m_Mricopy[i]=volume<float>(images[i-1]);
  bapusedflag=bapused;
  maximum=images[0].max();
  minimum=images[0].min();
  iterationspve=biterationspve;
  m_Mri=new volume<float>[numberofchannels+1];
  pveBmixeltype=fpveMRFmixeltyp;
  m_nSlicesize=m_nWidth*m_nHeight;
  noclasses=nclasses;
  m_nbIter=nbiter;
  m_post=volume4D<float>(m_nWidth, m_nHeight, m_nDepth, noclasses+1);
  m_post.copyproperties(images[0]);
  m_prob=volume4D<float>(m_nWidth, m_nHeight, m_nDepth, noclasses+1);
  m_prob.copyproperties(images[0]);
  members=volume4D<float>(m_nWidth, m_nHeight, m_nDepth, noclasses);
  members.copyproperties(images[0]);
  m_mask=volume<int>(m_nWidth, m_nHeight, m_nDepth);
  copybasicproperties(images[0], m_mask);
  m_maskc=volume<float>(m_nWidth, m_nHeight, m_nDepth);
  m_maskc.copyproperties(images[0]);
  m_Segment=volume<int>(m_nWidth, m_nHeight, m_nDepth);
  copybasicproperties(images[0], m_Segment);
  
  biasfieldremoval=bbias; 
  beta=fbeta;
  if(true)
    {
      m_pveSegment=volume<int>(m_nWidth, m_nHeight, m_nDepth);
      m_pveSegment.copyproperties(m_Segment);
    }
  initfixed=initinitfixed;
  m_nbIter=nbiter;
  if(bapusedflag>0)
    {
      talpriors=volume4D<float>(m_nWidth, m_nHeight, m_nDepth, noclasses+1);
      talpriors.copyproperties(m_post);
    }
  Hyper=Hyp;
}

void ZMRIMULTISegmentation::Dimensions()
{
  amz=1.0/m_nzdim;
  amy=1.0/m_nydim;
  amzy=1.0/sqrt(m_nzdim*m_nzdim+m_nydim*m_nydim);
  amx=1.0/m_nxdim;
  amzx=1.0/sqrt(m_nzdim*m_nzdim+m_nxdim*m_nxdim);
  amxy=1.0/sqrt(m_nxdim*m_nxdim+m_nydim*m_nydim);
}

int ZMRIMULTISegmentation::TanakaMain(NEWIMAGE::volume<float>& pcsf, NEWIMAGE::volume<float>& pgm, NEWIMAGE::volume<float>& pwm)
{
  for(int z=0;z<m_nDepth;z++)
    {
      for(int y=0;y<m_nHeight;y++)
	{ 
	  for(int x=0;x<m_nWidth;x++)
	    {
	      m_mask.value(x, y, z)=1;
	      for(int i=1;i<numberofchannels+1;i++)
		{

		  if(m_Mricopy[i].value(x, y, z)>0.0)
		    {
		      m_Mricopy[i].value(x, y, z)=log(m_Mricopy[i].value(x, y, z) + 1.0);
		    }
		  else 
		    {
		      m_Mricopy[i].value(x, y, z)=0.0f;
		      m_mask.value(x, y, z)=0;
		    }
		}
	    }
	}
    }
  Dimensions();
  volumequant=new double[noclasses+1];
  InitKernel();
  m_Finalbias=new volume<float>[numberofchannels+1];
  m_resmean=new volume<float>[numberofchannels+1];
  m_meaninvcov=new volume<float>[numberofchannels+1];
  p_bias=new volume<float>[numberofchannels+1];
  for(int s=1;s<=numberofchannels;s++)
    {
      m_Finalbias[s]=volume<float>(m_nWidth, m_nHeight, m_nDepth);
      m_Finalbias[s].copyproperties(m_Mri[s]);
      m_resmean[s]=volume<float>(m_nWidth, m_nHeight, m_nDepth);
      m_resmean[s].copyproperties(m_Mri[s]);
      m_meaninvcov[s]=volume<float>(m_nWidth, m_nHeight, m_nDepth);
      m_meaninvcov[s].copyproperties(m_Mri[s]);
      p_bias[s]=volume<float>(m_nWidth, m_nHeight, m_nDepth);
      p_bias[s].copyproperties(m_Mri[s]);
    }
  long seed=-1;
  srand(seed);
  rhs=new float[100]; 
  for(int i=1;i<numberofchannels+1;i++)
    m_Mri[i]=volume<float>(m_Mricopy[i]);
  InitWeights();
  if(bapusedflag==0)
  {
    try{Initialise();}
    catch (kmeansexception& km){cout << km.what() << endl ; return -1;}
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

  // first loop to remove bias field
  BiasRemoval();
  for(int iter=0;iter<m_nbIter;iter++)
    {
      if(verboseusage)
	cout<<"Tanaka Iteration "<<iter<<" bias field "<<m_nbIter<<endl;
      TanakaIterations();
      UpdateMembers(m_post);
      BiasRemoval();
      MeansVariances(noclasses);
    }
  // second loop to estimate hyperparameter
  for(int iter=0;iter<initfixed;iter++)
    {
      if(verboseusage)
	cout<<"Tanaka Iteration "<<iter<<" hyperparameter "<<initfixed<<endl;
      TanakaIterations();
      UpdateMembers(m_post);
      if(Hyper>=0.0)
	beta=Hyper;
      else
	TanakaHyper();
     if(verboseusage) cout<<" BETA "<<beta<<endl;
      MeansVariances(noclasses);
    }


  m_Finalbias=p_bias;
  takeexpo();
  qsort();
  Classification();
  UpdateMembers(m_post);
  if(iterationspve>0)
    {
      PVMoreSophisticated();
      Volumesquant(m_pve);
      pveClassification();
    }
  return 0;
}

void ZMRIMULTISegmentation::Initialise()
{
  WeightedKMeans();
  UpdateMembers(m_post);
  MeansVariances(noclasses);
  m_prob=Initclass(noclasses);
  m_post=m_prob;
  UpdateMembers(m_post);
  Classification();
}

void ZMRIMULTISegmentation::Classification(int x, int y, int z)
{
  float max=-1e-10;
  m_Segment(x, y, z)=0;
  if(m_mask(x, y, z)==1)
    {
      for(int c=1;c<noclasses+1;c++)
	{
	  if(weight[c-1]*m_post(x, y, z, c)>max)
	    {
	      max=weight[c-1]*m_post(x, y, z, c);
	      m_Segment(x, y, z)=c;
	    }
	}
    }
}


void ZMRIMULTISegmentation::Classification()
{
  for(int z=0;z<m_nDepth;z++)
      for(int y=0;y<m_nHeight;y++)
	  for(int x=0;x<m_nWidth;x++)
	      Classification(x, y, z);
}


void ZMRIMULTISegmentation::pveClassification(int x, int y, int z)
{
  float max=-1e-10;
  m_pveSegment(x, y, z)=0;
  if(m_mask(x, y, z)==1)
    {
      for(int c=1;c<noclasses+1;c++)
	{
	  if(max<m_pve(x, y, z, c))
	    {
	      max=m_pve(x, y, z, c);
	      m_pveSegment(x, y, z)=c;
	    }
	}
    }
}

void ZMRIMULTISegmentation::pveClassification()
{
   for(int x=0;x<m_nWidth;x++)
     for(int y=0;y<m_nHeight;y++)
       for(int z=0;z<m_nDepth;z++)
	 pveClassification(x, y, z);
}

float ZMRIMULTISegmentation::logGaussian(int classnumber, int x, int y, int z)
{
  if((m_mask.value(x, y, z)>0)&&(classnumber<noclasses+1))
    {
      Matrix submeanrow(1, numberofchannels);
      for(int i=1;i<=numberofchannels;i++)submeanrow(1, i)=m_Mri[i](x, y, z)-m_mean(classnumber, i);
      Matrix submeancol(numberofchannels, 1);
      for(int i=1;i<=numberofchannels;i++)submeancol(i, 1)=m_Mri[i](x, y, z)-m_mean(classnumber, i);
      float sum=(submeanrow*m_inv_co_variance[classnumber]*submeancol).AsScalar();    
      sum=0.5*sum+log(sqrt(abs(m_co_variance[classnumber].Determinant()*M_2PI(numberofchannels))));
      return sum;
    }
  else return 0.0;
}

float ZMRIMULTISegmentation::MRFWeightsTotal()
{
  double total=0.0f; //Internally a double to avoid truncation when adding small to large (shouldn't happen anyway with this loop style)
  for(int z=0;z<m_nDepth;z++)
      for(int y=0;y<m_nHeight;y++)
	  for(int x=0;x<m_nWidth;x++)
	      if(m_mask.value(x, y, z)==1)
		  for(int c=1;c<noclasses+1;c++)
		    total+=MRFWeightsInner(x,y,z,c)*m_post(x, y, z, c);
  return (float)total;   
}

double ZMRIMULTISegmentation::MRFWeightsInner(const int x, const int y, const int z,const int c)
{
double inner=0.0f;
    for(int n=-1;n<=1;n++)
	for(int m=-1;m<=1;m++)
	     for(int l=-1;l<=1;l++)
		if(m_mask(x+l, y+m, z+n)==1)
		  inner+=MRFWeightsAM(l,m,n)*m_post.value(x+l, y+m, z+n, c);
   return inner;
}

double ZMRIMULTISegmentation::MRFWeightsAM(const int l, const int m, const int n)
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

void ZMRIMULTISegmentation::TanakaHyper()
{
  float total=MRFWeightsTotal();
  float betahtemp=0.0f;
  beta=betahtemp;
  float min=abs(total-rhs[0]);
  for(int betahyp=1;betahyp<15;betahyp++)
    {
      betahtemp+=0.05f;
      if(abs(total-rhs[betahyp])<min)
	{
	  beta=betahtemp;
	  min=abs(total-rhs[betahyp]);
	}
    }
}

void ZMRIMULTISegmentation::TanakaPriorHyper()
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
		  for(int c=1;c<noclasses+1;c++)
		    {
		      if(m_mask(x, y, z)==1)
			{
			  sum+=m_prob(x, y, z, c)=(float)(rand()/(float)(RAND_MAX));
			}
		      else
			m_prob(x, y, z, c)=0.0f;		    
		    }
		  for(int c=1;c<noclasses+1;c++)
		    if(sum>0.0)
		      m_prob(x, y, z, c)/=sum;
		}
	    }
	}
      for(int iteration=0;iteration<5;iteration++)
	{
	  m_post=m_prob;
	  for(int z=0;z<m_nDepth;z++)
	    {
	      for(int y=0;y<m_nHeight;y++)
		{
		  for(int x=0;x<m_nWidth;x++)
		    {
		      float sum=0.0f;
		      if(m_mask.value(x, y, z)==1)
			{
			  for(int c=1;c<noclasses+1;c++)
			    {
			      float post=MRFWeightsInner(x,y,z,c);
			      if(bapusedflag<2)
			      sum+=m_prob.value(x, y, z, c)=exp(betahtemp*post);
			      else
				sum+=m_prob.value(x, y, z, c)=talpriors(x, y, z ,c)*exp(betahtemp*post);				
			    }
			  for(int c=1;c<noclasses+1;c++)
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

void ZMRIMULTISegmentation::TanakaIterations()
{
  for(int iteration=0;iteration<5;iteration++)
    {
      for(int z=0;z<m_nDepth;z++)
	{
	  for(int y=0;y<m_nHeight;y++)
	    {
	      for(int x=0;x<m_nWidth;x++)
		{
		  float sum=0.0f;
		  if(m_mask.value(x, y, z)==1)
		    {
		      for(int c=1;c<noclasses+1;c++)
			{
			  float post=MRFWeightsInner(x,y,z,c);
			  if(bapusedflag<2)
			    sum+=m_post.value(x, y, z, c)=exp(beta*post-logGaussian(c, x, y, z));
			  else
			    sum+=m_post.value(x, y, z, c)=talpriors(x, y, z, c)*exp(beta*post-logGaussian(c, x, y, z));
		      	}
		      for(int c=1;c<noclasses+1;c++)
			{
			  if(sum>0.0f)
			    m_post.value(x, y, z, c)/=sum;
			  else
			    m_post.value(x, y, z, c)=0.0f;
			}
		    }
		}
	    }
	}
    }
  if(verboseusage)
    {
      for(int channel=1;channel<=numberofchannels;channel++)
	{
	  cout<<" CHANNEL "<<channel;
	  for(int c=1;c<=noclasses;c++)
	    cout<<" CLASS "<<c<<" MEAN "<<exp(m_mean(c,channel))<<" STDDEV "<<( exp(m_mean(c,channel)+sqrt(m_co_variance[c](channel,channel))) - exp(m_mean(c,channel)-sqrt(m_co_variance[c](channel,channel))) )/2;
	  cout<<endl;
	}
    }
}

float ZMRIMULTISegmentation::PVEnergy(int x, int y, int z, Matrix  mu, Matrix sigma, float detsigma)
{
  if(m_mask(x, y, z)==1)
    {
      Matrix tempval(numberofchannels, 1);
      Matrix rowval(1, numberofchannels);
      for(int i=1;i<=numberofchannels;i++)
	{
	  tempval(i, 1)=rowval(1, i)=m_Mri[i](x, y, z)-mu(i, 1);
	}
      return ((rowval*(sigma.i())*tempval)).AsScalar()+log(detsigma);
    }
  else  
    return 0.0f;
}

void ZMRIMULTISegmentation::PVClassificationStep()
{
  pvvar(0);
  int mixnum=0;
  if(noclasses==3)
    mixnum=6;
  if(noclasses==2)
    mixnum=3;
  if(noclasses>3)
    mixnum=2*noclasses-1;
  float step=(float)(1.0f/(float)(iterationspve));
  PVprob=volume4D<float>(m_nWidth, m_nHeight, m_nDepth, mixnum);
  PVprob=0.0f;
  Matrix mu(numberofchannels, 1);
  Matrix sig(numberofchannels, numberofchannels);
  Matrix invsig(numberofchannels, numberofchannels);
  if(noclasses==3)
    {
      for(int z=0;z<m_nDepth;z++)
	{
	  for(int y=0;y<m_nHeight;y++)
	    {
	      for(int x=0;x<m_nWidth;x++)
		{ 
		  if(m_mask.value(x, y, z)==1)
		    {
		      for(int type=0;type<3;type++)
			{
			  if(noclasses==3)
			    {
			      if(type<3)
				{
				  for(int i=1;i<=numberofchannels;i++)
				    {
				      mu(i, 1)=m_mean(type+1, i);
				      for(int j=1;j<=numberofchannels;j++)
					{
					  sig(i, j)=m_co_variance[type+1](i, j);
					}
				    }
				  //float detsig=sig.Determinant();
				  //invsig=sig.i();
				  PVprob(x, y, z, type)=exp(-1.0*logGaussian(type+1,x, y, z));// mu, invsig, detsig));
				}
			    }
			}
		    }
		}
	    }
	}
      for(int type=3;type<mixnum;type++)
	{
	  for(float delta=0.0;delta<=1.0;delta+=step)
	    {
	      if(type==3)
		{
		  for(int i=1;i<=numberofchannels;i++)
		    {
		      mu(i, 1)=delta*m_mean(1, i)+(1-delta)*m_mean(2, i);
		      {
			for(int j=1;j<=numberofchannels;j++)
			  {
			    sig(i, j)=delta*delta*m_co_variance[1](i, j)+(1-delta)*(1-delta)*m_co_variance[2](i, j);
			  }
		      }
		    }
		}
	      if(type==4)
		{
		  for(int i=1;i<=numberofchannels;i++)
		    {
		      mu(i, 1)=delta*m_mean(1, i)+(1-delta)*m_mean(3, i);
		      {
			for(int j=1;j<=numberofchannels;j++)
			  {
			    sig(i, j)=delta*delta*m_co_variance[1](i, j)+(1-delta)*(1-delta)*m_co_variance[3](i, j);
			  }
		      }
		    }
		}
	      if(type==5)
		{
		  for(int i=1;i<=numberofchannels;i++)
		    {
		      mu(i, 1)=delta*m_mean(2, i)+(1-delta)*m_mean(3, i);
		      {
			for(int j=1;j<=numberofchannels;j++)
			  {
			    sig(i, j)=delta*delta*m_co_variance[2](i, j)+(1-delta)*(1-delta)*m_co_variance[3](i, j);
			  }
		      }
		    }
		}
	      float detsig=sig.Determinant();
	      invsig=sig.i();
	      for(int z=0;z<m_nDepth;z++)
		{
		  for(int y=0;y<m_nHeight;y++)
		    {
		      for(int x=0;x<m_nWidth;x++)
			{ 
			  if(m_mask.value(x, y, z)==1)
			    {
			      PVprob(x, y, z, type)+=exp(-1.0*logpveGaussian(x, y, z, mu, invsig, detsig ))*step;
			    }
			}
		    }
		}
	    }
	}
    }
  if(noclasses>3)
    {
      for(int z=0;z<m_nDepth;z++)
	{
	  for(int y=0;y<m_nHeight;y++)
	    {
	      for(int x=0;x<m_nWidth;x++)
		{ 
		  if(m_mask.value(x, y, z)==1)
		    {
		      for(int type=0;type<noclasses;type++)
			{
			  if(noclasses>3)
			    {
			      if(type<noclasses)
				{
				  for(int i=1;i<=numberofchannels;i++)
				    {
				      mu(i, 1)=m_mean(type+1, i);
				      for(int j=1;j<=numberofchannels;j++)
					{
					  sig(i, j)=m_co_variance[type+1](i, j);
					}
				    }
				  invsig=sig.i();
				  PVprob(x, y, z, type)=exp(-1.0*logpveGaussian(x, y, z, mu, invsig, sig.Determinant()));
				}
			    }
			}
		    }
		}
	    }
	}
      for(int type=noclasses;type<mixnum;type++)
	{
	  for(float delta=0.0;delta<=1.0;delta+=step)
	    {
	      for(int i=1;i<=numberofchannels;i++)
		{
		  mu(i, 1)=delta*m_mean(type-noclasses+1, i)+(1-delta)*m_mean(type-noclasses+2, i);
		  {
		    for(int j=1;j<=numberofchannels;j++)
		      {
			sig(i, j)=delta*delta*m_co_variance[type-noclasses+1](i, j)+(1-delta)*(1-delta)*m_co_variance[type-noclasses+2](i, j);
		      }
		  }
		}
	      float detsig=sig.Determinant();
	      invsig=sig.i();
	      for(int z=0;z<m_nDepth;z++)
		{
		  for(int y=0;y<m_nHeight;y++)
		    {
		      for(int x=0;x<m_nWidth;x++)
			{ 
			  if(m_mask.value(x, y, z)==1)
			    {
			      PVprob(x, y, z, type)+=exp(-1.0*logpveGaussian(x, y, z, mu, invsig, detsig ))*step;
			    }
			}
		    }
		}
	    }
	}
    }
  if(noclasses==2)
    {
      for(int type=0;type<mixnum;type++)    
	{
	  if(type<2)
	    {
	      for(int i=1;i<=numberofchannels;i++)
		mu(i, 1)=m_mean(type+1, i);
	      invsig=m_co_variance[type+1].i();
	      float detsig=m_co_variance[type+1].Determinant();
	      for(int z=0;z<m_nDepth;z++)
		{
		  for(int y=0;y<m_nHeight;y++)
		    {
		      for(int x=0;x<m_nWidth;x++)
			{ 
			  if(m_mask.value(x, y, z)==1)
			    {
			      PVprob(x, y, z, type)=exp(-1.0*logpveGaussian(x, y, z, mu, invsig, detsig));
			    }
			}
		    }
		}
	    }
	  else
	    {
	      for(float delta=0.0;delta<=1.0;delta+=step)
		{
		  {
		    for(int i=1;i<=numberofchannels;i++)
		      {
			mu(i, 1)=delta*m_mean(type-noclasses+1, i)+(1-delta)*m_mean(type-noclasses+2, i);
			{
			  for(int j=1;j<=numberofchannels;j++)
			    {
			      sig(i, j)=delta*delta*m_co_variance[type-noclasses+1](i, j)+(1-delta)*(1-delta)*m_co_variance[type-noclasses+2](i, j);
			    }
			}
		      }
		    float detsig=sig.Determinant();
		    invsig=sig.i();
		    for(int z=0;z<m_nDepth;z++)
		      {
			for(int y=0;y<m_nHeight;y++)
			  {
			    for(int x=0;x<m_nWidth;x++)
			      { 
				if(m_mask.value(x, y, z)==1)
				  {
				    PVprob(x, y, z, type)+=exp(-1.0*logpveGaussian(x, y, z, mu, invsig, detsig ))*step;
				  }
			      }
			  }
		      }
		  }
		}
	    }
	}
    }
  ICMPV();
}

float ZMRIMULTISegmentation::logpveGaussian(int x, int y, int z, Matrix mu, Matrix sig, float detsig)
{
  if((m_mask.value(x, y, z)>0))
    {
      Matrix submeanrow(1, numberofchannels);
      Matrix submeancol(numberofchannels, 1);
      for(int i=1;i<=numberofchannels;i++)submeancol(i, 1)=submeanrow(1, i)=m_Mri[i](x, y, z)-mu(i, 1);
      float sum=(submeanrow*sig*submeancol).AsScalar();    
      sum=0.5*sum+log(sqrt(abs(detsig*M_2PI(numberofchannels))));
      return sum;
    }
  else return 0.0;
}

void ZMRIMULTISegmentation::ICMPV()
{
  hardPV=volume<int>(m_nWidth, m_nHeight, m_nDepth);
  hardPV.copyproperties(m_Segment);
  int mixnum=0;
  if(noclasses==3)
    mixnum=6;
  if(noclasses==2)
    mixnum=3;
  if(noclasses>3)
    mixnum=2*noclasses-1;
  float* clique=new float[mixnum];
  for(int z=0;z<m_nDepth;z++)
    {
      for(int y=0;y<m_nHeight;y++)
	{
	  for(int x=0;x<m_nWidth;x++)
	    { 
	      if(m_mask.value(x, y, z)==1)
		    {
		      float max=-1.0;
		      for(int type=0;type<mixnum;type++)
			{
			  if(max<PVprob(x, y, z, type))
			    {
			      hardPV(x, y, z)=type;
			      max=PVprob(x, y, z, type);
			    }
			}
		    }
	    }
	}
    }
  for(int iter=0;iter<1;iter++)
    {
      for(int z=0;z<m_nDepth;z++)
	{
	  for(int y=0;y<m_nHeight;y++)
	    {
	      for(int x=0;x<m_nWidth;x++)
		{ 
		  if(m_mask.value(x, y, z)==1)
		    {
		      for(int type=0;type<mixnum;type++)
			clique[type]=0.0f;
		      for(int n=-1;n<=1;n++)
			{
			  for(int m=-1;m<=1;m++)
			    {
			      for(int l=-1;l<=1;l++)
				{
				  if(m_mask(x+l, y+m, z+n)==1)
				    {
                                      float am=MRFWeightsAM(l, m, n);
				      if(noclasses==3)
					{
					  for(int type=0;type<6;type++)
					    {	
					      if(type==hardPV(x+l, y+m, z+n))
						{
						  clique[type]+=am*2;
						  continue;
						}
					      if((type==0)&&((hardPV(x+l, y+m, z+n)==3)||(hardPV(x+l, y+m, z+n)==4)))
						{
						  clique[type]+=am;
						  continue;
						}
					      if((type==1)&&((hardPV(x+l, y+m, z+n)==3)||(hardPV(x+l, y+m, z+n)==5)))
						{
						  clique[type]+=am;
						  continue;
						}
					      if((type==2)&&((hardPV(x+l, y+m, z+n)==4)||(hardPV(x+l, y+m, z+n)==5)))
						{
						  clique[type]+=am;
						  continue;
						}
					      if((type==3)&&((hardPV(x+l, y+m, z+n)==0)||(hardPV(x+l, y+m, z+n)==1)))
						{
						  clique[type]+=am;
						  continue;
						}
					      if((type==4)&&((hardPV(x+l, y+m, z+n)==0)||(hardPV(x+l, y+m, z+n)==2)))
						{
						  clique[type]+=am;
						  continue;
						}
					      if((type==5)&&((hardPV(x+l, y+m, z+n)==1)||(hardPV(x+l, y+m, z+n)==2)))
						{
						  clique[type]+=am;
						  continue;
						}
					      clique[type]-=am;
					    }
					}
				      if(noclasses>3)
					{
					  for(int type=0;type<mixnum;type++)
					    {	
					      if(type==hardPV(x+l, y+m, z+n))
						{
						  clique[type]+=am*2;
						  continue;
						}
					      if((0<type)&&(type<noclasses-1)&&((hardPV(x+l, y+m, z+n)==(noclasses+type-1))||(hardPV(x+l, y+m, z+n)==(noclasses+type))))
						{
						  clique[type]+=am;
						  continue;
						}
					      if((type==0)&&((hardPV(x+l, y+m, z+n)==noclasses)))
						{
						  clique[type]+=am;
						  continue;
						}
					      if((type==(noclasses-1))&&((hardPV(x+l, y+m, z+n)==mixnum-1)))
						{
						  clique[type]+=am;
						  continue;
						}
					      if((type>(noclasses-1))&&((hardPV(x+l, y+m, z+n)==(type-noclasses))||(hardPV(x+l, y+m, z+n)==(type-noclasses+1))))
						{
						  clique[type]+=am;
						  continue;
						}
					      clique[type]-=am;
					    }
					}
				      if(noclasses==2)
					{
					  for(int type=0;type<3;type++)
					    {
					      if(type==hardPV(x+l, y+m, z+n))
						{
						  clique[type]+=am*2;
						  continue;
						}
					      if((type==0)&&((hardPV(x+l, y+m, z+n)==2)))
						{
						  clique[type]+=am;
						  continue;
						}
					      if((type==1)&&((hardPV(x+l, y+m, z+n)==2)))
						{
						  clique[type]+=am;
						  continue;
						}
					      if((type==2)&&((hardPV(x+l, y+m, z+n)==0)||(hardPV(x+l, y+m, z+n)==1)))
						{
						  clique[type]+=am;
						  continue;
						}
					      clique[type]-=am;
					    }
					}
				    }
				}
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
		}
	    }
	}
    }
  delete[] clique;
}

void ZMRIMULTISegmentation::PVMoreSophisticated()
{
  if(verboseusage)
    cout<<"Partial Volume Estimation \n";
  for(int iter=0;iter<1;iter++)
    {
      PVClassificationStep();
      PVestimation();
    }
}

void ZMRIMULTISegmentation::PVestimation()
{
  m_pve=volume4D<float>(m_nWidth, m_nHeight, m_nDepth, noclasses+1);
  m_pve.copyproperties(m_Mri[1]);
  m_pve=0.0f;
  float step=(float)(1.0f/(float)(iterationspve));
  Matrix mu(numberofchannels, 1);
  Matrix sig(numberofchannels, numberofchannels);  
  for(int z=0;z<m_nDepth;z++)
    {
      for(int y=0;y<m_nHeight;y++)
	{
	  for(int x=0;x<m_nWidth;x++)
	    { 
	     if(m_mask.value(x, y, z)>0)
	       {
		 float val=0.0;
		 if(noclasses==3)
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
			 continue;;
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
			     for(int j=1;j<=numberofchannels;j++)
			       {
				 mu(j, 1)=delta*m_mean(1, j)+(1-delta)*m_mean(2, j);
				 for(int k=1;k<=numberofchannels;k++)
				   {
				     sig(j, k)=delta*delta*m_co_variance[1](j, k)+(1-delta)*(1-delta)*m_co_variance[2](j, k);
				   }
			       }
			   }
			 if(hardPV(x, y, z)==4)
			   {
			     for(int j=1;j<=numberofchannels;j++)
			       {
				 mu(j, 1)=delta*m_mean(1, j)+(1-delta)*m_mean(3, j);
				 for(int k=1;k<=numberofchannels;k++)
				   {
				     sig(j, k)=delta*delta*m_co_variance[1](j, k)+(1-delta)*(1-delta)*m_co_variance[3](j, k);
				   }
			       }
			   }
			 if(hardPV(x, y, z)==5)
			   {
			     for(int j=1;j<=numberofchannels;j++)
			       {
				 mu(j, 1)=delta*m_mean(2, j)+(1-delta)*m_mean(3, j);
				 for(int k=1;k<=numberofchannels;k++)
				   {
				     sig(j, k)=delta*delta*m_co_variance[2](j, k)+(1-delta)*(1-delta)*m_co_variance[3](j, k);
				   }
			       }
			   }
			 val=PVEnergy(x, y, z, mu, sig, sig.Determinant());
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
		 if(noclasses>3)
		   {
		     float min=1.0e13;
		     if(hardPV(x, y, z)<noclasses)
		       {
			 for(int c=0;c<noclasses;c++)
			   {
			     if(hardPV(x, y, z)!=c)
			       m_pve(x, y, z, c+1)=0.0;
			   }
			 m_pve(x, y, z, hardPV(x, y, z)+1)=1.0f;
			 continue;
		       }
		     if(hardPV(x, y, z)>=noclasses)
		       {
			 for(float delta=0.00;delta<=1.0;delta+=step)
			   {
			     for(int j=1;j<=numberofchannels;j++)
			       {
				 mu(j, 1)=delta*m_mean(hardPV(x, y, z)-noclasses+1, j)+(1-delta)*m_mean(hardPV(x, y, z)-noclasses+2, j);
				 for(int k=1;k<=numberofchannels;k++)
				   {
				     sig(j, k)=delta*delta*m_co_variance[hardPV(x, y, z)-noclasses+1](j, k)+(1-delta)*(1-delta)*m_co_variance[hardPV(x, y, z)-noclasses+2](j, k);

				   }
			       }
			     val=PVEnergy(x, y, z, mu, sig, sig.Determinant());
			     if(min>val)
			       {
				 for(int c=0;c<noclasses;c++)
				   {
				     m_pve(x, y, z, c+1)=0.0f;
				   }
				 m_pve(x, y, z, hardPV(x, y, z)-noclasses+1)=delta;
				 m_pve(x, y, z, hardPV(x, y, z)-noclasses+2)=1.0-delta;			     
				 min=val;
			       }
			   }
		       }
		   }
		 if(noclasses==2)
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
			 continue;;
		       }
		     for(float delta=0.00;delta<=1.0;delta+=step)
		       {
			 if(hardPV(x, y, z)==2)
			   {
			     for(int j=1;j<=numberofchannels;j++)
			       {
				 mu(j, 1)=delta*m_mean(1, j)+(1-delta)*m_mean(2, j);
				 for(int k=1;k<=numberofchannels;k++)
				   {
				     sig(j, k)=delta*delta*m_co_variance[1](j, k)+(1-delta)*(1-delta)*m_co_variance[2](j, k);
				   }
			       }
			   }
			 float val=PVEnergy(x, y, z, mu, sig, sig.Determinant());
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

float ZMRIMULTISegmentation::pvmeans(int clas)
{
return 0.0f;
}

float ZMRIMULTISegmentation::pvvar(int clas)
{
return 0.0f;
}

void ZMRIMULTISegmentation::takeexpo()
{  
  for(int chan=1;chan<=numberofchannels;chan++)
    {
      for(int z=0;z<m_nDepth;z++)
	{
	  for(int y=0;y<m_nHeight;y++)
	    {
	      for(int x=0;x<m_nWidth;x++)
		{
		  m_Finalbias[chan](x, y, z)=exp(-m_Finalbias[chan](x, y, z));
		}
	    }
	}
    }
}

void ZMRIMULTISegmentation::InitWeights()
{
  weight=new float[noclasses];
  for(int c=0;c<noclasses;c++)
    weight[c]=1.0f/(float)(noclasses);
}

void ZMRIMULTISegmentation::InitSimple(const NEWIMAGE::volume<float>& pcsf, const NEWIMAGE::volume<float>& pgm, const NEWIMAGE::volume<float>& pwm)
{
  if(verboseusage)
    cout<<"Beginning prior-based initialisation"<<endl;
  if(noclasses==3)
    {
      for(int z=0;z<m_nDepth;z++)
	{
	  for(int y=0;y<m_nHeight;y++)
	    {
	      for(int x=0;x<m_nWidth;x++)
		{
		  if(m_mask.value(x, y, z)==1)
		    {
		      talpriors.value(x, y, z, 0)=talpriors(x, y, z, 1)=talpriors(x, y, z, 2)=talpriors(x, y, z, 3)=0.0f;
		      float norm2=pcsf.value(x, y, z)+pgm.value(x, y, z)+pwm.value(x, y, z);
		      if(norm2>0.0)
			{
			  m_post.value(x, y, z, 1)=m_prob.value(x, y, z, 1)=talpriors.value(x, y, z, 1)=pcsf.value(x, y, z)/norm2;
			  m_post.value(x, y, z, 2)=m_prob.value(x, y, z, 2)=talpriors.value(x, y, z, 2)=pgm.value(x, y, z)/norm2;
			  m_post.value(x, y, z, 3)=m_prob.value(x, y, z, 3)=talpriors.value(x, y, z, 3)=pwm.value(x, y, z)/norm2;
			}
		      else
			{
			  m_post.value(x, y, z, 1)=m_prob.value(x, y, z, 1)=talpriors.value(x, y, z, 1)=1.0/3.0f;
			  m_post.value(x, y, z, 2)=m_prob.value(x, y, z, 2)=talpriors.value(x, y, z, 2)=1.0/3.0f;
			  m_post.value(x, y, z, 3)=m_prob.value(x, y, z, 3)=talpriors.value(x, y, z, 3)=1.0/3.0f;
			}
		    }
		  else
		    {
		      m_post.value(x, y, z, 1)=m_prob.value(x, y, z, 1)=talpriors.value(x, y, z, 1)=0.0;
		      m_post.value(x, y, z, 2)=m_prob.value(x, y, z, 2)=talpriors.value(x, y, z, 2)=0.0;
		      m_post.value(x, y, z, 3)=m_prob.value(x, y, z, 3)=talpriors.value(x, y, z, 3)=0.0;		  
		    }
		  Classification(x, y, z);
		}
	    }
	}
    }
  if(noclasses==2)
    {
      for(int z=0;z<m_nDepth;z++)
	{
	  for(int y=0;y<m_nHeight;y++)
	    {
	      for(int x=0;x<m_nWidth;x++)
		{
		  if(m_mask.value(x, y, z)==1)
		    {
		      talpriors.value(x, y, z, 0)=talpriors(x, y, z, 1)=talpriors(x, y, z, 2)==0.0f;
		      float norm2=pcsf.value(x, y, z)+pgm.value(x, y, z)+pwm.value(x, y, z);
		      if(norm2>0.0)
			{
			  m_post.value(x, y, z, 1)=m_prob.value(x, y, z, 1)=talpriors.value(x, y, z, 1)=pcsf.value(x, y, z)/norm2;
			  m_post.value(x, y, z, 2)=m_prob.value(x, y, z, 2)=talpriors.value(x, y, z, 2)=(pgm.value(x, y, z)+pwm.value(x, y, z))/norm2;
			}
		      else
			{
			  m_post.value(x, y, z, 1)=m_prob.value(x, y, z, 1)=talpriors.value(x, y, z, 1)=1.0/2.0f;
			  m_post.value(x, y, z, 2)=m_prob.value(x, y, z, 2)=talpriors.value(x, y, z, 2)=1.0/2.0f;
			}
		    }
		  else
		    {
		      m_post.value(x, y, z, 1)=m_prob.value(x, y, z, 1)=talpriors.value(x, y, z, 1)=0.0;
		      m_post.value(x, y, z, 2)=m_prob.value(x, y, z, 2)=talpriors.value(x, y, z, 2)=0.0;		  
		    }
		  Classification(x, y, z);
		}
	    }
	}
    }
  UpdateMembers(m_post);
  MeansVariances(noclasses);
}

void ZMRIMULTISegmentation::Volumesquant(const NEWIMAGE::volume4D<float>& probs)
{
  double tot=0.0;
  for(int c=1;c<noclasses+1;c++)
    {
      volumequant[c]=0.0;
      for(int z=0;z<m_nDepth;z++)
	{
	  for(int y=0;y<m_nHeight;y++)
	    {
	      for(int x=0;x<m_nWidth;x++)
		{
		  if(m_mask.value(x, y, z)==1)
		    {
		      volumequant[c]+=probs.value(x, y, z, c);
		    }
		}
	    }
	}
      volumequant[c]*=m_nxdim;
      volumequant[c]*=m_nydim;
      volumequant[c]*=m_nzdim;
      tot+=volumequant[c];
      if(verboseusage)
	cout<<"\n tissue "<<c<<" " << volumequant[c];
    }
  if(verboseusage)
    cout<<"\n total tissue "<<tot<<"\n";
}

void ZMRIMULTISegmentation::InitKernel()
{
  float sigma=(sqrt(m_nbLowpass/m_nxdim));
  int radius=2*(int)sigma;
  kernelx=gaussian_kernel1D(sigma, radius);
  sigma=(sqrt(m_nbLowpass/m_nydim));
  radius=2*(int)sigma;
  kernely=gaussian_kernel1D(sigma, radius);
  sigma=(sqrt(m_nbLowpass/m_nzdim));
  radius=2*(int)sigma;
  kernelz=gaussian_kernel1D(sigma, radius);
}

NEWIMAGE::volume<float> ZMRIMULTISegmentation::Convolve(NEWIMAGE::volume<float>& resfieldimage)
{
  return convolve_separable(resfieldimage, kernelx, kernely, kernelz);
}

void ZMRIMULTISegmentation::MeansVariances(int numberofclasses)
{
  float* normtemp=new float[noclasses+1];
  for(int c=1;c<noclasses+1;c++)
    normtemp[c]=0.0f;
  for(int c=1;c<noclasses+1;c++)
    {
      m_co_variance[c]=0.0f;
      for(int i=1;i<=numberofchannels;i++)
	{
	  m_mean(c, i)=0.0f;
	}
    }
  for(int z=0;z<m_nDepth;z++)
    {
      for(int y=0;y<m_nHeight;y++)
	{
	  for(int x=0;x<m_nWidth;x++)
	    {
	      if(m_mask(x, y, z)==1)
		{
		  for(int c=1;c<noclasses+1;c++)
		    {
		      float ppp=members(x, y, z, c-1);
		      for(int i=1;i<=numberofchannels;i++)
			{
			  m_mean(c, i)+=ppp*m_Mri[i](x, y, z);
			}
		      normtemp[c]+=ppp;
		    }
		}
	    }
	}
    }
  for(int c=1;c<=noclasses;c++)
    {      
      if(normtemp[c]!=0.0)
	for(int i=1;i<=numberofchannels;i++)
	  m_mean(c, i)=m_mean(c, i)/normtemp[c];
      m_co_variance[c]=covariancematrix(c-1, members);
      m_inv_co_variance[c]=m_co_variance[c].i();
    }
   delete[] normtemp;  
}

void ZMRIMULTISegmentation::BiasRemoval()
{
  if(biasfieldremoval)
    {
      mresmatrix();
      for(int i=1;i<=numberofchannels;i++)
	{
	  m_resmean[i]=Convolve(m_resmean[i]);
	  m_meaninvcov[i]=Convolve(m_meaninvcov[i]);
	}
      for(int ch=1;ch<=numberofchannels;ch++)
	{    
	  for(int z=0;z<m_nDepth;z++)
	    {
	      for(int y=0;y<m_nHeight;y++)
		{
		  for(int x=0;x<m_nWidth;x++)
		    {
		      if(((m_mask(x, y, z)==0)||(m_meaninvcov[ch](x, y, z)==0.0f)))
			{
			  m_resmean[ch](x, y, z)=0.0f;
			  m_meaninvcov[ch](x, y, z)=1.0f;
			}
		      }
		}
	    }
	}
      for(int i=1;i<=numberofchannels;i++)
	{
	  p_bias[i]=m_resmean[i]/m_meaninvcov[i];
	}
      for(int i=1;i<=numberofchannels;i++)
	{
	  for(int z=0;z<m_nDepth;z++)
	    {
	      for(int y=0;y<m_nHeight;y++)
		{
		  for(int x=0;x<m_nWidth;x++)
		    {
		      if(m_mask.value(x, y, z)==0)
			p_bias[i].value(x, y, z)=0.0f;
		    }
		}
	    }
	  m_Mri[i]=m_Mricopy[i]-p_bias[i];
	}
    }
}

NEWIMAGE::volume4D<float> ZMRIMULTISegmentation::InitclassAlt(int numberofsegs)
{
  volume4D<float> probability;
  probability=volume4D<float>(m_nWidth, m_nHeight, m_nDepth, noclasses+1);
  probability=0.0;
  if(numberofsegs==3)
    {
      for(int z=0;z<m_nDepth;z++)
	{
	  for(int y=0;y<m_nHeight;y++)
	    {
	      for(int x=0;x<m_nWidth;x++)
		{
		  if(m_mask.value(x, y, z)==1)
		    {
		      float tot=0.0f;
		      for(int c=1;c<numberofsegs+1;c++)
			{
			  probability.value(x, y, z, c)=logGaussian(c, x, y, z);
			  tot+=probability(x, y, z, c)=exp(-probability(x, y, z, c))*talpriors(x, y, z, c);
			}
		      for(int c=1;c<numberofsegs+1;c++)
			{
			  if((tot>0.0)&&(m_mask.value(x, y, z)==1))
			    {
			      probability.value(x, y, z, c)/=tot;
			    }
			  else 
			    {
			      probability.value(x, y, z, c)=0.0;
			    }
			}
		    }
		}
	    }
	}
    }
  return probability;
}

NEWIMAGE::volume4D<float> ZMRIMULTISegmentation::Initclass(int noclasses)
{
   for(int z=0;z<m_nDepth;z++)
   {
     for(int y=0;y<m_nHeight;y++)
       {
	 for(int x=0;x<m_nWidth;x++)
	   {
	     float tot=0.0f;
	     for(int c=1;c<noclasses+1;c++)
	       {
		 if(m_mask(x, y, z)>0)
		   tot += m_prob(x, y, z, c) = exp(-1.0* logGaussian(c, x, y, z));
	       }
	     for(int c=1;c<noclasses+1;c++)
	       {
		 if((tot>0.0)&&(m_mask(x, y, z)==1))
		   m_prob(x, y, z, c)/=tot;
		 else 
		   m_prob(x, y, z, c)=0.0;
	       }
	     Classification(x, y, z);
	   }
       }
   }
   return m_prob;
}

void ZMRIMULTISegmentation::UpdateMembers(NEWIMAGE::volume4D<float>& probability)
{
  for(int x=0;x<m_nWidth;x++)
    {
      for(int y=0;y<m_nHeight;y++)
	{
	  for(int z=0;z<m_nDepth;z++)
	    { 
	      if(m_mask(x, y, z)>0)
		{
		  float sum=0.0f;
		  for(int c=0;c<noclasses;c++)
		    {
		      sum+=members(x, y, z, c)=weight[c]*probability(x, y, z, c+1);
		    }
		  for(int c=0;c<noclasses;c++)
		    {
		      if(sum>0.0)
			members(x, y, z, c)/=sum;
		      else
			members(x, y, z, c)=0.0f;
		    }
		}
	    }
	}
    }

}

void ZMRIMULTISegmentation::UpdateWeights()
{
  float sum=0.0f;
  for(int c=0;c<noclasses;c++)
    weight[c]=0.0f;
  for(int x=0;x<m_nWidth;x++)
    {
      for(int y=0;y<m_nHeight;y++)
	{
	  for(int z=0;z<m_nDepth;z++)
	    { 
	      if(m_mask(x, y, z)>0)
		{

		  for(int c=0;c<noclasses;c++)
		    {
		      sum+=weight[c]+=members(x, y, z, c);
		    }
		}
	    }
	}
    }
  for(int c=0;c<noclasses;c++)
    {
      weight[c]/=sum;
      if(verboseusage)
	cout<<" weight "<<weight[c];
    }
  if(verboseusage)
    cout<<"\n";
  
}

void ZMRIMULTISegmentation::WeightedKMeans()
{
  vector<float> inputMeans;
  if (mansegfile!="") 
  {  
    ifstream inputfile(mansegfile.c_str());
    copy(istream_iterator<float> (inputfile),istream_iterator<float> (),back_inserter(inputMeans));
    inputfile.close();
  }

  m_post=m_prob=0.0f;
  for(int z=0;z<m_nDepth;z++)
    for(int y=0;y<m_nHeight;y++)
      for(int x=0;x<m_nWidth;x++)
      {
	 if(m_mask(x, y, z)==1)
	   m_maskc(x, y, z)=1.0f;
	 else
	   m_maskc(x, y, z)=0.00f;
      }
	
    
  m_mean.ReSize(noclasses+1, numberofchannels);m_mean=0.0;
  m_co_variance=new Matrix[noclasses+1];
  m_inv_co_variance=new Matrix[noclasses+1];
  for(int i=1;i<=noclasses;i++)
    {
      m_co_variance[i].ReSize(numberofchannels, numberofchannels);
      m_inv_co_variance[i].ReSize(numberofchannels, numberofchannels);
      m_co_variance[i]=0.0f;
      m_inv_co_variance[i]=0.0f;

    }
  for(int n=1;n<=numberofchannels;n++)
   {
      float perc=1.0/((float)(noclasses+1.0));
      for(int c=1;c<=noclasses;c++)
	{
          if ( (int)inputMeans.size() == (noclasses*numberofchannels) ) m_mean(c, n)=log(inputMeans[c+noclasses*(n-1)-1]); 
	  else m_mean(c, n)=m_Mricopy[n].percentile((float)(perc*c), m_maskc);
           if (verboseusage) cout << n << " " << c << " " << m_mean(c, n) << endl ;
	}
    }

  for(int n=1;n<numberofchannels+1;n++)
      for(int c=1;c<noclasses+1;c++)
        for(int c2=1;c2<c;c2++)
	  if(m_mean(c, n)==m_mean(c2, n)) throw kmeansexc;
  
          
       

  for(int z=0;z<m_nDepth;z++)
    {
      for(int y=0;y<m_nHeight;y++)
	{
	  for(int x=0;x<m_nWidth;x++)
	    {
	      if(m_mask(x, y, z)==1)
		{
		  int minclass=0;
		  float mindist=1e31;
		  for(int c=1;c<noclasses+1;c++)
		    {
		      float sum=0.0f;
		      for(int n=1;n<numberofchannels+1;n++)
			for(int m=1;m<numberofchannels+1;m++)
			  sum+=(m_Mricopy[n].value(x, y, z)-m_mean(c, n))*(m_Mricopy[m].value(x, y, z)-m_mean(c, m));
		      if(sum<mindist)
			{
			  mindist=sum;
			  minclass=c;
			}
		    }
		  m_post(x, y, z, minclass)=1.0f;
		}
	    }
	}
    }
 
  for(int initfiter=0;initfiter<inititerations;initfiter++)
    {
      UpdateMembers(m_post);
      MeansVariances(noclasses);
      if(verboseusage)
	cout<<"KMeans Iteration "<<initfiter<<"\n";
      m_post=m_prob=Initclass(noclasses);
    }
}

int ZMRIMULTISegmentation::qsort()
{
    int c;
    Matrix meancopy=m_mean;
    double*  tempmean=new double[noclasses+1];
    Matrix* varcopy=new Matrix[noclasses+1];
    Matrix* tempvar=new Matrix[noclasses+1];
    Matrix* invvarcopy=new Matrix[noclasses+1];
    Matrix* tempinvvar=new Matrix[noclasses+1];
    if(imagetype==2) // for a T2 image reverse the intensity order of _master image_
      for(c=1; c<noclasses+1; c++)
	m_mean(c,1)*=-1.0f;
    for(c=1; c<noclasses+1; c++)
     {	 
        tempmean[c]=m_mean(c,1);
        varcopy[c]=m_co_variance[c];
        tempvar[c]=m_co_variance[c];
 	tempinvvar[c]=m_inv_co_variance[c];
	invvarcopy[c]=m_inv_co_variance[c];
    }   
    volume4D<float> postcopy;
    postcopy=volume4D<float>(m_post);
    volume4D<float> probcopy;
    probcopy=volume4D<float>(m_prob);
    sort(tempmean+1, tempmean+noclasses+1);
    for (c=1; c<noclasses+1; c++)
      for (int m=1; m<noclasses+1; m++)
	if(m_mean(c,1)==tempmean[m])
	{
           for(int n=2;n<numberofchannels+1;n++)
                 m_mean(c,n)=meancopy(m,n);
	   tempvar[c]=varcopy[m];
	   tempinvvar[c]=invvarcopy[m];
	   for(int z=0;z<m_nDepth;z++)
	     for(int y=0;y<m_nHeight;y++)
	       for(int x=0;x<m_nWidth;x++)
		   if(m_mask.value(x, y, z)==1)
		   {
		      m_post.value(x, y, z ,c)=postcopy.value(x, y, z, m);
		      m_prob.value(x, y, z ,c)=probcopy.value(x, y, z, m);
		      Classification(x, y, z);
		    }  
	}
    for(c=1; c<noclasses+1; c++)
      {	
        m_mean(c,1)=tempmean[c];
        m_inv_co_variance[c]=tempinvvar[c];
	m_co_variance[c]=tempvar[c];   
      }
    if(imagetype==2)
      for(c=1; c<noclasses+1; c++)
	  m_mean(c,1)*=-1.0f;
  return 0;
}


//Start of Multichannel only methods

void ZMRIMULTISegmentation:: InitAprioriKMeans()
{
  m_mean.ReSize(noclasses+1, numberofchannels);m_mean=0.0;
  m_co_variance=new Matrix[noclasses+1];
  m_inv_co_variance=new Matrix[noclasses+1];
  for(int i=1;i<noclasses+1;i++)
    {
      m_co_variance[i].ReSize(numberofchannels, numberofchannels);
      m_inv_co_variance[i].ReSize(numberofchannels, numberofchannels);
      m_co_variance[i]=0.0f;
      m_inv_co_variance[i]=0.0f;
    }
  if(bapusedflag==1)
    for(int initfiter=0;initfiter<inititerations;initfiter++)
    {
      UpdateMembers(m_post);
      MeansVariances(noclasses);
      if(verboseusage)
	cout<<"Weighted KMeans Iteration "<<initfiter<<"\n";
      m_post=m_prob=Initclass(noclasses);
    }
  else
    {
      for(int initfiter=0;initfiter<inititerations;initfiter++)
	{
	  UpdateMembers(m_post);
	  MeansVariances(noclasses);
	  if(verboseusage)
	    cout<<"Weighted KMeans Iteration "<<initfiter<<"\n";
	  m_post=m_prob=InitclassAlt(noclasses);
	}
    }
}

void ZMRIMULTISegmentation::mresmatrix()
{
  Matrix submeancol(numberofchannels, 1);
  Matrix onecol(numberofchannels, 1);
  for(int channel=1;channel<=numberofchannels;channel++)
    {
      for(int z=0;z<m_nDepth;z++)
	{
	  for(int y=0;y<m_nHeight;y++)
	    {
	      for(int x=0;x<m_nWidth;x++)
		{ 
		  m_resmean[channel](x, y, z)=0.0f;
		  m_meaninvcov[channel](x, y, z)=0.0f;
		  if(m_mask(x, y, z)==1)
		    {
		      float tempM=0.0f;
		      for(int c=1;c<noclasses+1;c++)
			{
			  tempM=members(x, y, z, c-1)/m_co_variance[c](channel, channel);
			  m_resmean[channel](x, y, z)+=tempM*(m_Mri[channel](x, y, z)-m_mean(c, channel));
			  m_meaninvcov[channel](x, y, z)+=tempM;
			}	  
		    }
		  else
		    {
		      m_meaninvcov[channel](x, y, z)=0.0f;
		      m_resmean[channel](x, y, z)=0.0f;
		    }
		}
	    }
	}
    }
}


void ZMRIMULTISegmentation:: PVMeansVar(const NEWIMAGE::volume4D<float>& m_posttemp)
{
}

Matrix ZMRIMULTISegmentation:: covariancematrix(int classid, volume4D<float> probability)
{
  float tot=0.0f;
  Matrix cov(numberofchannels, numberofchannels);cov=0.0;
  if(classid<noclasses)
    {
      for(int z=0;z<m_nDepth;z++)
	{
	  for(int y=0;y<m_nHeight;y++)
	    {
	      for(int x=0;x<m_nWidth;x++)
		{ 
		  if(m_mask(x, y, z)==1)
		    {
		      for(int i=1;i<=numberofchannels;i++)
			{
			  for(int j=1;j<=numberofchannels;j++)
			    {
			      cov(i,j)+=probability(x, y, z, classid)*(m_Mri[i](x, y, z)-m_mean(classid+1, i))*(m_Mri[j](x, y, z)-m_mean(classid+1, j));
			    }
			}
		      tot+=probability(x, y, z, classid);
		    }
		}
	    }
	}
    }
   
  cov/=tot;
  return cov; 
}

float ZMRIMULTISegmentation::M_2PI(int numberofchan)
{
  float val=1.0;
  for(int i=0;i<numberofchan;i++)
    val=val*2*PI;
  return val;
}
