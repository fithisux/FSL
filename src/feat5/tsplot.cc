/* {{{ Copyright etc. */

/*  tsplot - FMRI time series and model plotting

    Stephen Smith, Mark Woolrich and Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 1999-2008 University of Oxford  */

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

/* }}} */
/* {{{ background theory */

/*

GLM : Y = Xb + e

The "partial model fit" shows, in the case of a contrast which selects
a single EV, the part of the full model fit explained by that EV. In
the case of a more complex contrast, it is basically a plot of the sum
of each EV weighted by the product of that EV's PE and that EV's
contrast weighting.

i.e. the partial model fit is X*diag(c)*b, where X is the design and c
is the contrast vector (renormalised to unit length and then turned
into a diagonal matrix) and b is the parameter estimate vector.

Thus we plot this versus the "reduced data" ie plot

Y - Xb + X*diag(c)*b  vs  X*diag(c)*b
i.e.
residuals + partial fit    vs   partial fit

NOTE: this plot cannot simply be used to generate the t/Z score
associated with this contrast (eg by straight correlation) - this
would not be correct. In order to do that you would need to correlate
the Hansen projection c'*pinv(X) with Y instead.

*/

/* }}} */
/* {{{ defines, includes and typedefs */
#include <iomanip>
#include "featlib.h"
#include "libvis/miscplot.h"
#include <math.h>
#include "utils/fsl_isfinite.h"
#include <vector>
 
using namespace NEWMAT;
using namespace NEWIMAGE;
using namespace MISCPLOT;
using namespace std;
 
/*
void setupq (ColumnVector &x,ColumnVector &dx,ColumnVector &y,int npoint,Matrix &v,Matrix &a )
{
double diff,prev;
   v(1,4) = x(2) - x(1);
   for(int i=2;i<npoint;i++)
   {
      v(i,4) = x(i+1) - x(i);
      v(i,1) = dx(i-1)/v(i-1,4);
      v(i,2) = - dx(i)/v(i,4) - dx(i)/v(i-1,4);
      v(i,3) = dx(i+1)/v(i,4);
   }
   v(npoint,1) = 0;
   for(int i=2;i<npoint;i++) v(i,5) = v(i,1)*v(i,1) + v(i,2)*v(i,2) + v(i,3)*v(i,3);

   if ( npoint >= 4 ) for(int i=3;i<npoint;i++) v(i-1,6) = v(i-1,2)*v(i,1) + v(i-1,3)*v(i,2);

    v(npoint-1,6) = 0;

   if (npoint >= 5) for(int i=4;i<npoint;i++) v(i-2,7) = v(i-2,3)*v(i,1);

   v(npoint-2,7) = 0;
   v(npoint-1,7) = 0;
//construct  q-transp. * y  in  qty.
   prev = (y(2) - y(1))/v(1,4);
   for(int i=2;i<npoint;i++)
   {  
      diff = (y(i+1)-y(i))/v(i,4);
      a(i,4) = diff - prev;
      prev = diff;
   }
}

//  from  * a practical guide to splines *  by c. de boor    
//  from  * a practical guide to splines *  by c. de boor    
// to be called in  s m o o t h constructs the upper three diags. in v(i,j), i=2,,npoint-1, j=1,3, of
//  the matrix  6*(1-p)*q-transp.*(d**2)*q + p*r, then computes its
//  l*l-transp. decomposition and stores it also in v, then applies
//  forward and backsubstitution to the right side q-transp.*y in  qty
//  to obtain the solution in  u .
//  a(1,4)=qty a(1,3)=u a(1,1)=qu
void chol1d (double p, Matrix &v, Matrix &a,int &npoint)
{
double prev,ratio,six1mp,twop;
   six1mp = 6*(1-p);
   twop = 2*p;
   for(int i=2;i<npoint;i++)
   {
      v(i,1) = six1mp*v(i,5) + twop*(v(i-1,4)+v(i,4));
      v(i,2) = six1mp*v(i,6) + p*v(i,4);
      v(i,3) = six1mp*v(i,7);
   }
   if (npoint < 4)
   {     
     a(1,3) = 0;
     a(2,3) = a(2,4)/v(2,1);
     a(3,3) = 0;
     goto fortyone;
   }
//  factorization
   for(int i=2;i<npoint-1;i++)
   {
     ratio = v(i,2)/v(i,1);
     v(i+1,1) = v(i+1,1) - ratio*v(i,2);
     v(i+1,2) = v(i+1,2) - ratio*v(i,3);
     v(i,2) = ratio;
     ratio = v(i,3)/v(i,1);
     v(i+2,1) = v(i+2,1) - ratio*v(i,3);
     v(i,3) = ratio;
   }
//  forward substitution
   a(1,3) = 0;
   v(1,3) = 0;
   a(2,3) = a(2,4);
   for(int i=2;i<npoint-1;i++) a(i+1,3) = a(i+1,4) - v(i,2)*a(i,3) - v(i-1,3)*a(i-1,3);
//  back substitution
   a(npoint,3) = 0;
   a(npoint-1,3) = a(npoint-1,3)/v(npoint-1,1);
   int i = npoint-2;
   do
   {
     a(i,3) = a(i,3)/v(i,1)-a(i+1,3)*v(i,2)-a(i+2,3)*v(i,3);
   } while (--i > 1);
//  construct q*u
   fortyone:
   prev = 0;
   for (int i=2;i<=npoint;i++)
   {
     a(i,1) = (a(i,3) - a(i-1,3))/v(i-1,4);
     a(i-1,1) = a(i,1) - prev;
     prev = a(i,1);
   }
   a(npoint,1) = -a(npoint,1);
}

void smooth(ColumnVector &x,ColumnVector &y,ColumnVector &dy,int npoint,double s)
{
  Matrix a(npoint,4);
  Matrix v(npoint,7);
  double change,ooss,oosf,p,prevsf,prevq,q=0,sfq,sixp,six1mp,utru;
  setupq(x,dy,y,npoint,v,a);

  if ( s > 0 )                    
  {
   p = 0;                     
   chol1d(p,v,a,npoint);
   sfq = 0;
   for (int i=1;i<=npoint;i++) sfq = sfq + pow(a(i,1)*dy(i),2.0);
   sfq*=36;
   if (sfq < s) goto sixty;
   utru = 0;
   for (int i=2;i<=npoint;i++) utru+= v(i-1,4)*(a(i-1,3)*(a(i-1,3)+a(i,3))+pow(a(i,3),2.0));
   ooss = 1./sqrt(s);
   oosf = 1./sqrt(sfq);
   q = -(oosf-ooss)*sfq/(6.*utru*oosf);
   prevq = 0;
   prevsf = oosf;

   thirty:
   chol1d(q/(1.+q),v,a,npoint);
   sfq = 0;
   for(int i=1;i<=npoint;i++) sfq = sfq + pow(a(i,1)*dy(i),2.0);
   sfq*=36.0/pow(1+q,2.0);
   if (abs(sfq-s) < 0.01*s) goto fiftynine;
   oosf = 1.0/sqrt(sfq);
   change = (q-prevq)/(oosf-prevsf)*(oosf-ooss);
   prevq = q;
   q-= change;
   prevsf = oosf;
   goto thirty;
  }
  else 
  {
    p = 1;                     
    chol1d(p,v,a,npoint);
    sfq = 0;
    goto sixty;             
   }


   fiftynine: 
   p = q/(1.0+q);
//correct value of p has been found.
//compute pol.coefficients from  Q*u (in a(.,1)).
   sixty: 
   six1mp = 6./(1.+q);
   for(int i=1;i<=npoint;i++) a(i,1) = y(i) - six1mp*pow(dy(i),2.0)*a(i,1);
   sixp = 6*p;
   for(int i=1;i<=npoint;i++) 
   {
     a(i,3)*=sixp;
     y(i)=a(i,1); 
   }
   for(int i=1;i<npoint;i++)  
   {
     a(i,4) = (a(i+1,3)-a(i,3))/v(i,4);
     a(i,2) = (a(i+1,1)-a(i,1))/v(i,4)- (a(i,3)+a(i,4)/3.*v(i,4))/2.*v(i,4);
   }
}
*/

void usage(const string& message)
{
  if (!message.empty()) cout << message << endl;
  printf("Usage: tsplot <feat_directory.feat> [options]\n");
  printf("[-f <4D_data>] input main filtered data, in case it's not <feat_directory.feat>/filtered_func_data\n");
  printf("[-c <X Y Z>] : use X,Y,Z instead of max Z stat position\n");
  printf("[-C <X Y Z output_file.txt>] : use X,Y,Z to output time series only - no stats or modelling\n");
  printf("[-m <mask>] : use mask image instead of thresholded activation images\n");
  printf("[-o <output_directory>] change output directory from default of input feat directory\n");
  printf("[-n] don't weight cluster averaging with Z stats\n");
  printf("[-p] prewhiten data and model timeseries before plotting\n");
  printf("[-d] don't keep raw data text files\n");
  exit(1);
}

int main(int argc, char **argv)
{
ofstream     outputFile;
int          numEVs, npts, numContrasts=1, nftests=0, GRPHSIZE(600), PSSIZE(600); 
vector<double> normalisedContrasts, model, triggers;
 string       fmriFileName, fslPath, featdir, vType, indexText;
ColumnVector NewimageVoxCoord(4),NiftiVoxCoord(4);
bool outputText(true), useCoordinate(false), prewhiten(false), useTriggers(false), customMask(false), modelFree(false), isHigherLevel(false), outputDataOnly(false);
bool zWeightClusters(true);
volume<float> immask;

  NewimageVoxCoord << 0 << 0 << 0 << 1;
  NiftiVoxCoord << 0 << 0 << 0 << 1;

  /* process arguments */

  if (argc<2) usage("");
  featdir=string(argv[1]);
  fmriFileName=featdir+"/filtered_func_data";
  fslPath=string(getenv("FSLDIR"));

  string outputName(featdir);

  for (int argi=2;argi<argc;argi++)
  {
    if (!strcmp(argv[argi], "-f")) /* alternative fmri data */
    {
      if (argc<argi+2)
	usage("Error: no value given following -f");
      fmriFileName=string(argv[++argi]);
    }
    else if (!strcmp(argv[argi], "-c")) /* alternative voxel position */
    {
      useCoordinate=true;
      if (argc<=(argi+=3)) /* options following c haven't been given */
	usage("Error: incomplete values given following -c");
      NiftiVoxCoord << atoi(argv[argi-2]) << atoi(argv[argi-1]) << atoi(argv[argi]) << 1;
    }
    else if (!strcmp(argv[argi], "-C")) /* output data only */
    {
      outputDataOnly=useCoordinate=true;
      if (argc<=(argi+=4))
	usage("Error: incomplete values given following -C");
      NiftiVoxCoord << atoi(argv[argi-3]) << atoi(argv[argi-2]) << atoi(argv[argi-1]) << 1;
      outputName=string(argv[argi]);
    }
    else if (!strcmp(argv[argi], "-m")) /* alternative mask image */
    {
      customMask=true;
      if (argc<argi+2)
	usage("Error: no mask image given following -m");    
      if ( read_volume(immask,argv[++argi]) )
	usage("Error: mask image chosen doesn't exist");      
    }
    else if (!strcmp(argv[argi], "-o")) /* output dir */
    {
      if (argc<argi+2)
	usage("Error: no value given following -o"); 
      outputName=string(argv[++argi]);
    }
    else if (!strcmp(argv[argi], "-d")) outputText=false;
    else if (!strcmp(argv[argi], "-n")) zWeightClusters=false;   
    else if (!strcmp(argv[argi], "-p")) prewhiten=true; 
  }

  /* read filtered_func_data */

  volume4D<float> im;
  read_volume4D(im, fmriFileName);

  if (useCoordinate) NewimageVoxCoord = im.niftivox2newimagevox_mat()*NiftiVoxCoord;

  if (outputDataOnly && outputText) /* output raw data and exit */
  {
    outputFile.open(outputName.c_str());
    if(!outputFile.is_open())
    {
      cerr << "Can't open output data file " << outputName << endl;
      exit(1);
    }
    for(int t=0; t<im.tsize(); t++) outputFile << scientific << im((int)NewimageVoxCoord(1),(int)NewimageVoxCoord(2),(int)NewimageVoxCoord(3),t) << endl;
    outputFile.close();
    exit(0);
  }

  model=read_model(featdir+"/design.mat",&numEVs,&npts);
  
  if (npts==0)
  {
    modelFree=true;
    nftests=1;
    numContrasts=0;
  }

  npts=im.tsize();

ColumnVector TS_model(npts),TS_copemodel(npts),TS_pemodel(npts*numEVs),TS_data(npts),TS_residuals(npts);

  /* read auto correlation estimates for prewhitening */

vector<double> pwmodel;
volume4D<float> acs;

  if ( prewhiten ) {
    prewhiten=false;
    if(fsl_imageexists(featdir+"/stats/threshac1")) {
      read_volume4D(acs, featdir+"/stats/threshac1");
      if (acs[1].max()!=0) {/* hacky test for whether prewhitening was actually carried out */
	pwmodel.resize(numEVs*npts);
	prewhiten=true;
      }
    }
  }

  /* read design.con and PEs */

  vector< volume<float> > impe(numEVs);
  if (!modelFree)
  {
    Matrix contrasts=read_vest(featdir+"/design.con");
    numEVs=contrasts.Ncols();
    numContrasts=contrasts.Nrows();
    normalisedContrasts.resize( numEVs * numContrasts );
    for(int i=1; i<=numContrasts; i++)
      for(int ev=1; ev<=numEVs; ev++)
	normalisedContrasts[(i-1)*numEVs+(ev-1)] = contrasts(i,ev) / sqrt(contrasts.Row(i).SumSquare());
    
    for(int i=1;i<=numEVs;i++)
      read_volume(impe[i-1],featdir+"/stats/pe"+num2str(i));
  }

  if (!modelFree)
    read_ftests(featdir+"/design.fts",&nftests);
  useTriggers=read_triggers(featdir+"/design.trg",triggers,numEVs,npts);
  /* check analysis level */

  ifstream testFile((featdir+"/design.lev").c_str());
  isHigherLevel=testFile.is_open();
  testFile.close();

  /* create plot(s) for each contrast */
  miscplot newplot;
  for(int type=0;type<2;type++) /* setup stats type */
  {
    int nPlots(nftests);
    string statType("zfstat");
    if (type==0) { statType="zstat";  nPlots=numContrasts; }
    for(int i=1; i<=nPlots; i++)
    {
      volume<float> imcope, imz;
      bool haveclusters=false;
      string graphText("");
      string peristimulusText("");
      /* read COPE and derived stats; test for f-test output */
      /* load zstat or zfstat */
      if (fsl_imageexists(featdir+"/stats/"+statType+num2str(i))) 
	read_volume(imz,featdir+"/stats/"+statType+num2str(i));
      else 
	continue; /* f-test i wasn't valid - no zfstat image */
      /* load cope */
      if ( (type==0) && (!modelFree) )
	read_volume(imcope,featdir+"/stats/cope"+num2str(i));
      
      /* load cluster mask */
      if (!useCoordinate) {
	if (!customMask) {
	  if (fsl_imageexists(featdir+"/cluster_mask_"+statType+num2str(i)))
	    read_volume(immask,featdir+"/cluster_mask_"+statType+num2str(i));
	}
	haveclusters=(immask.max()>0);
      }
      /* find max Z and X,Y,Z */

      double maxz(-1000);
      if (!useCoordinate)
      {
	NewimageVoxCoord << 0 << 0 << 0 << 1;
	for(int z=0; z<im.zsize(); z++)
	  for(int y=0; y<im.ysize(); y++)
	    for(int x=0; x<im.xsize(); x++)
	      if ( (imz(x,y,z)>maxz) && ( (!haveclusters) ||  ( (immask(x,y,z)>0) && (!prewhiten || acs(x,y,z,1)!=0 || acs(x,y,z,2)!=0) ) ) )
	      {
		/* make max Z be inside a cluster if we found a cluster map */
		maxz=imz(x,y,z);
		NewimageVoxCoord << x << y << z << 1;
	      }
      }
      else
	maxz=imz((int)NewimageVoxCoord(1),(int)NewimageVoxCoord(2),(int)NewimageVoxCoord(3));

      /* first do peak voxel plotting then do mask-averaged plotting */
      for(int v=0;v<=1;v++)
      {
	double wtotal=0;
	int maskedVoxels=0;
	if (v==0) vType.clear();
	else vType="c";
	 
	/* {{{ create model and data time series */
	TS_model=0;
	TS_residuals=0;
	TS_copemodel=0;
	TS_data=0;
	TS_pemodel=0;
	ColumnVector prewhitenedTS;
	for(int x=0; x<im.xsize(); x++) 
	  for(int y=0; y<im.ysize(); y++) 
	    for(int z=0; z<im.zsize(); z++)
	      if ( ((v==0 && x==(int)NewimageVoxCoord(1) && y==(int)NewimageVoxCoord(2) && z==(int)NewimageVoxCoord(3)) || (v==1 && immask(x,y,z)>0)) && (!prewhiten || acs(x,y,z,1)!=0 || acs(x,y,z,2)!=0)) {
		maskedVoxels++;
		double weight(1);
		if (v!=0 && zWeightClusters) weight=imz(x,y,z);

		wtotal+=weight;
		if(prewhiten)
		  prewhiten_timeseries(acs.voxelts(x,y,z), im.voxelts(x,y,z), prewhitenedTS, npts);
		else
		   prewhitenedTS = im.voxelts(x,y,z);
		for(int t=1; t<=npts; t++) TS_data(t)+= prewhitenedTS(t)*weight;
		if (!modelFree) {
		  if (prewhiten)
		    prewhiten_model(acs.voxelts(x,y,z), model, pwmodel, numEVs, npts);
		  else
		    pwmodel=model;
		  for(int t=1; t<=npts; t++)
		    for(int ev=0; ev<numEVs; ev++)
		    {
		      double tmpf=pwmodel[(t-1)*numEVs+ev]*impe[ev](x,y,z)*weight;
		      TS_model(t)           += tmpf;
		      TS_copemodel(t)       += tmpf*normalisedContrasts[(i-1)*numEVs+ev];
		      TS_pemodel(ev*npts+t) += tmpf;
		    }
		}
	      }

	TS_data/=wtotal;
	double tsmean(TS_data.Sum()/npts);
	  
	if (isHigherLevel) tsmean=0;
	if (!modelFree)
	  for(int t=1; t<=npts; t++)
	  {
	    TS_model(t) = TS_model(t)/wtotal + tsmean;
	    TS_copemodel(t) = TS_copemodel(t)/wtotal + tsmean;
	    TS_residuals(t)=TS_data(t)-TS_model(t);
	    for(int ev=0; ev<numEVs; ev++)
	      TS_pemodel(ev*npts+t) = TS_pemodel(ev*npts+t)/wtotal + tsmean;
	  }
	/* output data text files */
	if (outputText) 
	  outputFile.open((outputName+"/tsplot"+vType+"_"+statType+num2str(i)+".txt").c_str());
	double ymin(TS_data(1)),ymax(TS_data(1));
	for(int t=1; t<=npts; t++)
	{
	  if (outputText) outputFile << scientific << TS_data(t);
	  ymin=MISCMATHS::Min(TS_data(t),ymin); 
	  ymax=MISCMATHS::Max(TS_data(t),ymax);
	  if (!modelFree)
	  {
	    if (type==0)
	    {
	      if (outputText) outputFile << " " << TS_copemodel(t); 
	      ymin=MISCMATHS::Min(TS_copemodel(t),ymin); 
	      ymax=MISCMATHS::Max(TS_copemodel(t),ymax);
	    }
	    if (outputText) outputFile << " " << TS_model(t); 
	    ymin=MISCMATHS::Min(TS_model(t),ymin); 
	    ymax=MISCMATHS::Max(TS_model(t),ymax);
	    if (type==0) outputFile << " " << TS_residuals(t)+TS_copemodel(t);
	  }
	  if (outputText) outputFile << endl;
	}
	if (outputText) 
	  outputFile.close();
	ymax+=(ymax-ymin)/5;
	ymin-=(ymax-ymin)/20;
	if (ymin==ymax) 
	  ymin-=1; 	
	/* create graphs */
	string graphName="tsplot"+vType+"_"+statType+num2str(i);
	string graphFileName(outputName+"/"+graphName);
	GRPHSIZE= MISCMATHS::Min(MISCMATHS::Max(npts*4,600),3000);
	newplot.set_minmaxscale(1.001);
	newplot.set_xysize(GRPHSIZE,192);
	newplot.set_yrange(ymin,ymax);
	string title=statType+num2str(i);
	if (v==0) 
	{
	  NiftiVoxCoord = im.niftivox2newimagevox_mat().i()*NewimageVoxCoord;
	  if (!useCoordinate) title+= ": max Z stat of "+num2str(maxz)+" at ";
	  else title+= ": Z stat of "+num2str(maxz)+" at selected ";
	  title+="voxel ("+num2str((int)NiftiVoxCoord(1))+" "+num2str((int)NiftiVoxCoord(2))+" "+num2str((int)NiftiVoxCoord(3))+")";
	} 
	else 
	  title+= ": averaged over "+num2str(maskedVoxels)+" voxels";
	Matrix blank=TS_data;
	blank=log(-1.0);
	if (!modelFree)
	{
	  if (type==0)
	  {
	    newplot.add_label("full model fit");
	    newplot.add_label("cope partial model fit");
	    newplot.add_label("data");
	    newplot.timeseries((TS_model | TS_copemodel | TS_data).t(),graphFileName,title,1,GRPHSIZE,4,2,false);
	    newplot.remove_labels(3);
	    newplot.add_label("");
	    newplot.add_label("cope partial model fit");
	    newplot.add_label("reduced data");
	    newplot.timeseries((blank | TS_copemodel | TS_residuals+TS_copemodel ).t(),graphFileName+"p",title,1,GRPHSIZE,4,2,false);
	    newplot.remove_labels(3);
	    graphText+="Full model fit - <a href=\""+graphName+"p.png\">Partial model fit</a> - <a href=\""+graphName+".txt\">Raw data</a><br>\n<IMG BORDER=0 SRC=\""+graphName+".png\"><br><br>\n";
	  }
	  else
	  {
	    newplot.add_label("full model fit");
	    newplot.add_label("");
	    newplot.add_label("data");
	    newplot.timeseries((TS_model|blank|TS_data).t(),graphFileName,title,1,GRPHSIZE,4,2,false);
	    newplot.remove_labels(3);
	    graphText+="Full model fit - <a href=\""+graphName+".txt\">Raw data</a><br>\n<IMG BORDER=0 SRC=\""+graphName+".png\"><br><br>\n";
	  }
	}
	else
	{    
	  newplot.add_label("");
	  newplot.add_label("");
	  newplot.add_label("data");
	  newplot.timeseries((blank | blank | TS_data).t(),graphFileName,title,1,GRPHSIZE,4,2,false);
	  newplot.remove_labels(3);
	  graphText+="Data plot - <a href=\""+graphName+".txt\">Raw data</a>\n<IMG BORDER=0 SRC=\""+graphName+".png\"><br><br>\n";
	}
	/* picture for main web index page */
	if (v==0) 
	  indexText+="<a href=\"" + graphName + ".html\"><IMG BORDER=0 SRC=\"" + graphName + ".png\"></a><br><br>\n";
  	/* peri-stimulus: output text and graphs */
	if (useTriggers)
	{
	  if (!modelFree)
	    peristimulusText+="<table><tr>\n";
	  for(int ev=0; ev<numEVs; ev++) 
	    if (triggers[ev]>0.5) {
	      float ps_period=triggers[((int)triggers[ev]+1)*numEVs+ev];
	      Matrix ps_compact((int)(10*ps_period)+1,3);
	      if (!modelFree) ps_compact.ReSize((int)(10*ps_period)+1,6);
	      Matrix ps_full(0,ps_compact.Ncols()-1);
	      ps_compact=0;
	      for(int which_event=1;which_event<=triggers[ev];which_event++)
	      {
		double min_t=triggers[which_event*numEVs+ev];
		int int_min_t=(int)ceil(min_t-(1e-10*min_t)),max_t=MISCMATHS::Min(npts-1,int_min_t+(int)ps_period);		
		for(int t=int_min_t+1;t<=max_t;t++)
		{
		  RowVector input(ps_compact.Ncols());
		  if (!modelFree) input << (ceil((t-min_t-1)*10))/10.0 << TS_residuals(t)+TS_model(t) << TS_model(t) << TS_pemodel(ev*npts+t) << TS_residuals(t)+TS_pemodel(ev*npts+t) << 1;  //(restricted temporal accuraccy (0.1*TR) must be at least 0.1 ( scatter can not take 0 )
		  else input << t-min_t-1 << TS_residuals(t)+TS_model(t) << 1;
		  ps_compact.Row(((int)((t-min_t-1)*10))+1)+=input; 
		  ps_full &= input.Columns(1,input.Ncols()-1);
		}
	      }
	      graphName="ps_tsplot"+vType+"_"+statType+num2str(i)+"_ev"+num2str(ev+1);
	      graphFileName=outputName+"/"+graphName;
	      if (outputText) {
		outputFile.open((graphFileName+".txt").c_str());
		for(int k=1;k<=ps_full.Nrows();k++)
		{ 
		  outputFile << setprecision(1) << fixed << ps_full(k,1) << setprecision(6) << scientific;
		  for (int j=2;j<=ps_full.Ncols();j++) outputFile << " " << ps_full(k,j);
		  outputFile << endl;
		}
		outputFile.close();
	      }
	      title=statType+num2str(i)+" ev"+num2str(ev+1);
	      for(int j=1;j<=ps_compact.Nrows();j++) 
	      {
		if (isfinite(ps_compact(j,6))) ps_compact.Row(j)/=ps_compact(j,6);
		else  ps_compact.Row(j)=log10(-1.0); //deliberately set to nan
	      }
	      Matrix ps_interp=ps_compact.t();
	      PSSIZE = MISCMATHS::Min(MISCMATHS::Max(ps_period*3,400),3000);
	      newplot.set_minmaxscale(1.001);
	      newplot.add_xlabel("peristimulus time (TRs)");
	      newplot.set_xysize(PSSIZE,192);
	      newplot.set_yrange(ymin,ymax);
	      if (v==0) {
		NiftiVoxCoord = im.niftivox2newimagevox_mat().i()*NewimageVoxCoord;
		if (!useCoordinate) title+= ": max Z stat of "+num2str(maxz)+" at ";
		else title+= ": Z stat of "+num2str(maxz)+" at selected ";
		title+="voxel ("+num2str((int)NiftiVoxCoord(1))+" "+num2str((int)NiftiVoxCoord(2))+" "+num2str((int)NiftiVoxCoord(3))+")"; 
	      } 
	      else title+= ": averaged over "+num2str(maskedVoxels)+" voxels";
	      if (!modelFree)
	      {
		ps_compact=ps_full.SubMatrix(1,ps_full.Nrows(),1,2);
		ps_compact.Column(1)*=10;
		newplot.setscatter(ps_compact,(int)(10*(ps_period+3)));
		newplot.add_label("full model fit");
		newplot.add_label("EV "+num2str(ev+1)+" model fit");
		newplot.add_label("data");
		newplot.timeseries(ps_interp.SubMatrix(3,4,1,ps_interp.Ncols()),graphFileName,title,-0.1,PSSIZE,3,2,false);
		newplot.remove_labels(3);
		ps_compact=ps_full.SubMatrix(1,ps_full.Nrows(),1,1) | ps_full.SubMatrix(1,ps_full.Nrows(),5,5);
		ps_compact.Column(1)*=10;
		newplot.setscatter(ps_compact,(int)(10*(ps_period+3)));
		newplot.add_label("");
		newplot.add_label("EV "+num2str(ev+1)+" model fit");
		newplot.add_label("reduced data");
		ps_interp.Row(3)=log(-1.0);
		newplot.timeseries(ps_interp.SubMatrix(3,4,1,ps_interp.Ncols()),graphFileName+"p",title,-0.1,PSSIZE,3,2,false);
		newplot.deletescatter();
		newplot.remove_labels(3);
		peristimulusText+="<td>Full model fit - <a href=\""+graphName+"p.png\">Partial model fit</a> - <a href=\""+graphName+".txt\">Raw data</a><br>\n<IMG BORDER=0 SRC=\""+graphName+".png\">\n";
	      }
	      else
	      {
		Matrix blank=ps_full.SubMatrix(2,2,1,ps_compact.Ncols());
		blank=log(-1.0);
		newplot.add_label("");
		newplot.add_label("");
		newplot.add_label("data");
		newplot.timeseries(blank & blank & ps_full.SubMatrix(2,2,1,ps_compact.Ncols()),graphFileName,title,-0.1,PSSIZE,3,2,false);
		newplot.remove_labels(3);
		peristimulusText+="%sData plot - <a href=\""+graphName+".txt\">Raw data</a>\n<IMG BORDER=0 SRC=\""+graphName+".png\"><br><br>\n";
	      }
	      newplot.remove_xlabel();
	    }
	  if (!modelFree) peristimulusText+="</tr></table><br><br>\n";
	}
	if (!haveclusters) break;
      }
      /* {{{ web output */
      outputFile.open((outputName+"/tsplot_"+statType+num2str(i)+".html").c_str());
      if(!outputFile.is_open())
      {
	cerr << "Can't open output report file " << outputName << endl;
	exit(1);
      }
      outputFile << "<HTML>\n<TITLE>"<< statType << num2str(i) <<"</TITLE>\n<BODY BACKGROUND=\"file:"<< fslPath <<"/doc/images/fsl-bg.jpg\">\n<hr><CENTER>\n<H1>FEAT Time Series Report - "<< statType << num2str(i) <<"</H1>\n</CENTER>\n<hr><b>Full plots</b><p>\n"<< graphText;
      if (useTriggers) outputFile << "\n<hr><b>Peristimulus plots</b><p>\n"<< peristimulusText <<"\n<HR></BODY></HTML>\n\n";
      else outputFile << "\n</BODY></HTML>\n\n";
      outputFile.close();
    }
  }

  /* main web index page output */
  /* first output full index page (eg for use by featquery) */
  outputFile.open((outputName+"/tsplot_index.html").c_str());
  if(!outputFile.is_open())
  {
      cerr << "Can't open output report file " << outputName << endl;
      exit(1);
  }
  outputFile << "<HTML>\n<TITLE>FEAT Time Series Report</TITLE>\n<BODY BACKGROUND=\"file:" << fslPath << "/doc/images/fsl-bg.jpg\">\n<hr><CENTER>\n<H1>FEAT Time Series Report</H1>\n</CENTER>\n<hr>" << indexText << "<HR></BODY></HTML>" << endl << endl;
  outputFile.close();

  /* now output same thing without start and end, for inclusion in feat report */
  outputFile.open((outputName+"/tsplot_index").c_str());
  if(!outputFile.is_open())
  {
      cerr << "Can't open output report file " << outputName << endl;
      exit(1);
  }
  outputFile << indexText << endl << endl;
  outputFile.close();

  exit(0);
}

