/*  signal2image.cc

    Ivana Drobnjak and Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2003 University of Oxford  */

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

// Application for converting signal output from possum into an image

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

#include "utils/options.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"

using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;

//Tejas-apod
#define pi 3.1416
//Tejas-end

// The two strings below specify the title and example usage that is
//  printed out as the help or usage message

string title="signal2image\nCopyright(c) 2003, University of Oxford (Mark Jenkinson & Ivana Drobnjak)";
string examples="signal2image [options] -i <signal> -p <pulse> -o <image> \n signal2image -p <pulse> -c <kcoord>";


Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<string> opt_pulse(string("-p,--pulse"), string(""),
                  string("8-column pulse_sequence matrix. Expects to find all other pulse sequence files in the same directory."),
		  true, requires_argument);
Option<bool> useabs(string("-a,--abs"), false,
		  string("save absolute magnitude and phase"),
		  false, no_argument);
Option<string> inname(string("-i,--in"), string(""),
		  string("input signal"),
		  false, requires_argument);
Option<string> outname(string("-o,--out"), string(""),
		  string("output image"),
		  false, requires_argument);
Option<string> koutname(string("-k,--kout"), string(""),
		  string("output k-space"),
		  false, requires_argument);
Option<string> opt_kcoord(string("-c, --kcoord"), string(""),
		  string("kspace coordinates"),
		  false, requires_argument);
Option<bool> opt_homo(string("--homo"), false,
		  string("do the homodyne reconstruction"),
		  false, no_argument);

//Tejas-apod
Option<bool> opt_doapod(string("-z,--apodize"),false,
			string("Do apodization"),
			false,no_argument);

Option<int> opt_cutoff(string("-l,--limit"),100,
			string("Apodization with this cutoff; default 100"),
			false,requires_argument);

Option<int> opt_rolloff(string("-r,--roll"),10,
			string("Apodization with this rolloff; default 10"),
			false,requires_argument);

//Debug
Option<string> opt_savewindow(string("-s,--save"),string(""),
				string("(DEBUG!) Save window as ascii matrix"),
				false,requires_argument);
//Tejas-end
int nonoptarg;

//TEJAS-apodization!
int do_apodization(Matrix& signal, const int sizeX, const int sizeY, const int sizeZ)
{
	int kc = opt_cutoff.value();
	int rolloff = opt_rolloff.value();
	int midpointX = sizeX / 2;
	int midpointY = sizeY / 2;
	double width;
	long counter = 1;

	string savewindow = opt_savewindow.value();
	bool dosave = !(savewindow.empty());

	//UPDATE values according to hanning window
	Matrix window(sizeX,sizeY);
	window = 0;

	//Check rolloff value. Should not exceed the gap between kc & boundary
	if ( rolloff > (int)(sizeX-kc)/2 )
		rolloff = (sizeX-kc)/2 - 1;
	if ( rolloff > (int) (sizeY-kc)/2 )
		rolloff = (sizeY-kc)/2 - 1;

	// Check if the cutoff is too large
	if (kc > sizeX || kc > sizeY)
	{
		cout << "ERROR! Cutoff cannot exceed image dimensions!" << endl;
		exit(1);
	}

	if (verbose.value())	cout << "Apodization: kc=" << kc << ";rolloff=" << rolloff << endl;

	//Use hanning window expression to calculate rolloff regions
	//Calculate the rolloff region: Along X
	width = rolloff;

	for (int j=1; j<=sizeY; j++)
	{
		//Values: hanning Left hand side
		for (double i=midpointX-kc/2-width; i<=midpointX-kc/2; i++)
		{
			window((int)i,j) = 0.5*(1-cos(pi*((i-midpointX-kc/2-width-1)/width)));
		}

		//Values: cos Right hand side
		for (double i=midpointX+kc/2+width; i>=midpointX+kc/2; i--)
		{
			window((int)i,j) = 0.5*(1-cos(pi*((i-midpointX+kc/2+width-1)/width)));
		}
	}

	//Calculate the rolloff region: Along Y
	for (int i=1; i<sizeX; i++)
	{
		//Values: Top
		for (double j=midpointY-kc/2-width; j<=midpointY-kc/2; j++)
		{
			if (window(i,(int)j) != 0)
				window((int)i,(int)j) *= 0.5*(1-cos(pi*((j-midpointY-kc/2-width-1)/width)));
			else
				window(i,(int)j) = 0.5*(1-cos(pi*((j-midpointY-kc/2-width-1)/width)));
		}
		//Values: Bottom
		for (double j=midpointY+kc/2+width; j>=midpointY+kc/2; j--)
		{
			if (window(i,(int)j) != 0)
				window(i,(int)j) *= 0.5*(1-cos(pi*((j-midpointY+kc/2+width-1)/width)));
			else
				window(i,(int)j) = 0.5*(1-cos(pi*((j-midpointY+kc/2+width-1)/width)));
		}
	}

	//Rewrite 0's to ensure 0 as outer values of window
	for (int i=1; i<=sizeX; i++)
	{
		for (int j=1; j<=sizeY; j++)
		{
			if ((i < midpointX-kc/2-width) || (i > midpointX+kc/2+width))
				window(i,j) = 0;
			if ((j < midpointY-kc/2-width) || (j > midpointY+kc/2+width))
				window(i,j) = 0;
		}
	}

	//Re-write 1's to ensure 1 from mid-point to kc/2 in both directions
	for (int i=midpointX-kc/2; i<=midpointX+kc/2; i++)
	{
		for (int j=midpointY-kc/2; j<=midpointY+kc/2 ; j++)
		{
			window(i,j) = 1;
		}
	}

	//Apply window to signal === ONLY WORKS FOR GE & EPI
	for (int k=1; k<=sizeZ; k++)
	{
		for (int j=1; j<=sizeY; j++)
		{
			for(int i=1; i<=sizeX; i++)
			{
				signal(1,counter) *= window(i,j);
				signal(2,counter) *= window(i,j);
				counter++;
			}
		}
	}

	//Write window as ascii matrix
	if (dosave)
	{
		write_ascii_matrix(window,savewindow,1);
		if (verbose.value()) cout << "window in ascii file '" << savewindow << "'" << endl;
	}

	return 0;
}
//Tejas-end


////////////////////////////////////////////////////////////////////////////

int ReshapeEpiSignal(const Matrix& signal,
		     const int slcdir,const int nslc,const int phasedir, const int nphase, const int readdir, const int nread, const int startkspace,
             volume4D<double>& kspace_real,volume4D<double>& kspace_imag) {
  cout<<"Reshaping the signal..."<<endl;
  int n=kspace_real.tsize();
  int nsize=signal.Ncols();
  kspace_real=0;
  kspace_imag=0;
  int slchelp=0;
  int simdir=1;
  int zdir1=1;
  int zdir2=1;  
  int ydir1=1;
  int ydir2=1; 
  int xdir1=1; 
  int xdir2=1;
  if (sign(slcdir)<0) {
    simdir=-1;
    slchelp=nslc-1;
  }
  int nphase_new=nphase-startkspace+1;
  for (int nn=1;nn<=n;nn++){
    for (int nzz=1;nzz<=nslc;nzz++){
      for (int k=nphase_new;k>=1;k=k-2){
	     int a=(nphase_new-k)*nread+(nzz-1)*nread*nphase_new+(nn-1)*nread*nphase_new*nslc+1;
	     int c=a+nread;
	     if (verbose.value()) { cout << "a,c,k = " << a << " " << c << " " << k << endl; }
	     for (int m=1;m<=nread;m++){
               if (abs(slcdir)==1) {
                 zdir1=slchelp+simdir*(nzz-1);
		 zdir2=slchelp+simdir*(nzz-1);
	       }
               if (abs(slcdir)==2) {
                 ydir1=slchelp+simdir*(nzz-1);
                 ydir2=slchelp+simdir*(nzz-1);
	       }
               if (abs(slcdir)==3) {
                 xdir1=slchelp+simdir*(nzz-1);
                 xdir2=slchelp+simdir*(nzz-1);
	       }
               if (phasedir==1) {
                 zdir1=(nphase-1)-(k-1);
		 zdir2=(nphase-1)-(k-2);
	       }
               if (phasedir==2) {
                 ydir1=(nphase-1)-(k-1);
		 ydir2=(nphase-1)-(k-2);
	       }
               if (phasedir==3) {
                 xdir1=(nphase-1)-(k-1);
		 xdir2=(nphase-1)-(k-2);
	       }
               if (readdir==1) {
                 zdir1=m-1;
		 zdir2=nread-m;
	       }
               if (readdir==2) {
                 ydir1=m-1;
		 ydir2=nread-m;
	       }
               if (readdir==3) {
                 xdir1=m-1;
		 xdir2=nread-m;
	       }
	       kspace_real(xdir1,ydir1,zdir1,nn-1)=signal(1,a+m-1);
               if (verbose.value()) {
                  cout<<signal(1,a+m-1)<<endl;
                  cout<<"xdir1="<<xdir1<<" ydir1="<<ydir1<<" zdir1="<<zdir1<<endl;
		  cout<<kspace_real(xdir1,ydir1,zdir1,nn-1)<<endl;
		  cout<<""<<endl;
	       }
	       kspace_imag(xdir1,ydir1,zdir1,nn-1)=signal(2,a+m-1);
               if (c < nsize) {
	           kspace_real(xdir2,ydir2,zdir2,nn-1)=signal(1,c+m-1);
	           kspace_imag(xdir2,ydir2,zdir2,nn-1)=signal(2,c+m-1);
	       }
	      }
	  }
      }
  }
  return 0;
}

//--------------------
int ReshapeGradEchoSignal(const Matrix& signal,const int slcdir,const int nslc,const int phasedir, const int nphase, const int readdir, const int nread,  volume4D<double>& kspace_real,
		  volume4D<double>& kspace_imag) 
{
  int n=kspace_real.tsize();
  int slchelp=0;
  int simdir=1;
  int xdir=1;
  int ydir=1;
  int zdir=1;
  if (sign(slcdir)<0) {
    simdir=-1;
    slchelp=nslc-1;
  }
  int counter=1;
  for (int nn=1;nn<=n;nn++){
    for (int nzz=1;nzz<=nslc;nzz++){
      for (int nyy=1;nyy<=nphase;nyy++){
	    for (int nxx=1;nxx<=nread;nxx++){
	      if (abs(slcdir)==1) zdir=slchelp+simdir*(nzz-1);
          if (abs(slcdir)==2) ydir=slchelp+simdir*(nzz-1);
          if (abs(slcdir)==3) xdir=slchelp+simdir*(nzz-1);
          if (phasedir==1) zdir=nyy-1;
          if (phasedir==2) ydir=nyy-1;
          if (phasedir==3) xdir=nyy-1;
          if (readdir==1)  zdir=nxx-1;
          if (readdir==2)  ydir=nxx-1;
          if (readdir==3)  xdir=nxx-1;
	      kspace_real(xdir,ydir,zdir,nn-1)=signal(1,counter);
	      kspace_imag(xdir,ydir,zdir,nn-1)=signal(2,counter);
          if (verbose.value()) { cout << signal(1,counter) << "  "<< signal(2,counter) << endl; }
	      counter=counter+1;
        }
      }
    }
  }
  return 0;
}

int setdir(int& xdir, int& ydir, int& zdir, const int x, const int y, const int z, const int readdir, const int phasedir,const int slcdir){
  if (abs(slcdir)==1){
    zdir=z;
  }
  if (abs(slcdir)==2){
    ydir=z;
  }
  if (abs(slcdir)==3){
    xdir=z;
  }
  if (phasedir==1){
    zdir=y;
  }
  if (phasedir==2){
    ydir=y;
  }
  if (phasedir==3){
    xdir=y;
  }
  if (readdir==1){
    zdir=x;
  }
  if (readdir==2){
    ydir=x;
  }
  if (readdir==3){ 
    xdir=x;
  }
  return 0;
}


int do_work(int argc, char* argv[]) 
{
  RowVector pulseinfo	;
  pulseinfo=read_ascii_matrix(opt_pulse.value()+".info");
  int n=(int) (pulseinfo(12));
  int nslc=(int) (pulseinfo(13));//Nslc
  double dt=pulseinfo(3);//TR
  double dslc=pulseinfo(14)*1e03;//slcthk (mm)
  double dread=pulseinfo(7)*1e03;//read
  double dphase=pulseinfo(8)*1e03;//phase
  int nread=(int) (pulseinfo(5));
  int nphase=(int) (pulseinfo(6));
  int seqnum=(int) (pulseinfo(1));//2 for ge 1 for epi and 0 for none
  int slcdir=(int) pulseinfo(15);//1 for z, 2 for y and 3 for x
  int phasedir=(int) pulseinfo(19);
  int readdir=(int) pulseinfo(20);
  int startkspace=1;
  if (pulseinfo.Ncols() >= 22){
    startkspace=(int) pulseinfo(22);
  } 
  if (slcdir==phasedir || slcdir==readdir || readdir==phasedir){
   cout<<"WARNING: The same gradients used for different directions in the k-space!!"<<endl;
   exit(EXIT_FAILURE);
  }
  int nx=nread;
  int ny=nphase;
  int nz=nslc;
  double dx=dread;
  double dy=dphase;
  double dz=dslc;
  if (abs(slcdir)==1) {
    nz=nslc;
    dz=dslc;
  }
  if (abs(slcdir)==2) {
    ny=nslc;
    dy=dslc;
  }
  if (abs(slcdir)==3) {
    nx=nslc;
    dx=nslc;
  }
  if (phasedir==1) {
    nz=nphase;
    dz=dphase;
  }
  if (phasedir==2) {
    ny=nphase;
    dy=dphase;
  }
  if (phasedir==3) {
    nx=nphase;
    dx=dphase;
  }
  if (readdir==1) {
    nz=nread;
    dz=dread;
  }
  if (readdir==2) {
    ny=nread;
    dy=dread;
  }
  if (readdir==3) {
    nx=nread;
    dx=dread;
  }
  cout<<"nx,ny,nz "<<nx<<" "<<ny<<" "<<nz<<" "<<endl;
  cout<<"dx,dy,dz "<<dx<<" "<<dy<<" "<<dz<<" "<<endl;

  if (opt_kcoord.set()){
    volume4D<double> kcoord_kx(nx,ny,nz,n);
    volume4D<double> kcoord_ky(nx,ny,nz,n);
    kcoord_kx.setxdim(dx);
    kcoord_kx.setydim(dy);
    kcoord_kx.setzdim(dz);
    kcoord_kx.settdim(dt);
    kcoord_ky.setxdim(dx);
    kcoord_ky.setydim(dy);
    kcoord_ky.setzdim(dz);
    kcoord_ky.settdim(dt);
    Matrix kcoord;
    kcoord=read_binary_matrix(opt_kcoord.value());
    if (seqnum==1) ReshapeEpiSignal(kcoord,slcdir,nslc,phasedir,nphase,readdir,nread,startkspace,kcoord_kx,kcoord_ky);
    else if (seqnum==2) ReshapeGradEchoSignal(kcoord,slcdir,nslc,phasedir,nphase,readdir,nread,kcoord_kx,kcoord_ky);
    else cout<<"Do not know the sequence number "<<seqnum<<endl;
    save_volume4D(kcoord_kx,opt_kcoord.value()+"_kx");
    save_volume4D(kcoord_ky,opt_kcoord.value()+"_ky");
  }

  if (inname.set()){
    Matrix signal;
    signal=read_binary_matrix(inname.value());

//Tejas-apod
	//Call apodization here
	if (opt_doapod.value())
	{
		cout << "Doing apodization..." << endl;
		do_apodization(signal,nread,nphase,nslc);
	}
//Tejas-end

    volume4D<double> kspace_real(nx,ny,nz,n);
    volume4D<double> kspace_imag(nx,ny,nz,n);
    kspace_real.setxdim(dx);
    kspace_real.setydim(dy);
    kspace_real.setzdim(dz);
    kspace_real.settdim(dt);
    kspace_imag.setxdim(dx);
    kspace_imag.setydim(dy);
    kspace_imag.setzdim(dz);
    kspace_imag.settdim(dt);
    if (seqnum==1) ReshapeEpiSignal(signal,slcdir,nslc,phasedir,nphase,readdir,nread,startkspace,kspace_real,kspace_imag);
    else if (seqnum==2) ReshapeGradEchoSignal(signal,slcdir,nslc,phasedir,nphase,readdir,nread,kspace_real,kspace_imag);
    else cout<<"Do not know the sequence"<<endl;

    if (koutname.set()) {
      if (useabs.value()) {
        volume4D<double> dummy(kspace_real);
        volume4D<double> dummy_phase(kspace_real);
        for (int nn=0; nn<kspace_real.tsize(); nn++) { 
 	  dummy[nn] = sqrt(kspace_real[nn]*kspace_real[nn] + kspace_imag[nn]*kspace_imag[nn]);
          for (int z=kspace_real.minz(); z<=kspace_real.maxz(); z++) {
	    for (int y=kspace_real.miny(); y<=kspace_real.maxy(); y++) {
	      for (int x=kspace_real.minx(); x<=kspace_real.maxx(); x++) {
	        dummy_phase(x,y,z,nn) = atan2(kspace_real(x,y,z,nn),kspace_imag(x,y,z,nn));
	      }
	    }
	  }
	}
        save_volume4D(dummy,koutname.value()+"_abs");
        save_volume4D(dummy_phase,koutname.value()+"_phase");
      }   else {
        save_volume4D(kspace_real,koutname.value()+"_real");
        save_volume4D(kspace_imag,koutname.value()+"_imag");
        //volume4D<float> kspace_real_float,kspace_imag_float;
        //copyconvert(kspace_real,kspace_real_float);
        //copyconvert(kspace_imag,kspace_imag_float);
        //save_complexvolume4D(kspace_real_float,kspace_imag_float,koutname.value());
      }
    }

    if (outname.set() && startkspace>1 && opt_homo.value()) {
        string aa="x";
        string bb="y";
        string cc="z";
        string aaa="x";
        string bbb="y";
        string ccc="z";
	//in the old version in order to make it be the same orientation as the images from the scanner I had to do swapdimensions("-x","-y","z") after the I did the fft2. The thing is the convention for the scanner is (y,x,z) and for me was (x,y,z) so maybe that had to do.will see...still testing this orientation thing. 
        if (abs(slcdir)==1){
          cc="z";
	  ccc="z";
	}
	if (abs(slcdir)==2){
          cc="y";
	  bbb="z";
	}
        if (abs(slcdir)==3){
	  cc="x";
	  aaa="z";
	}
        if (phasedir==1){
	  bb="z";
	  ccc="y";
	}
        if (phasedir==2){
	  bb="y";
	  bbb="y";
	}
        if (phasedir==3){
	  bb="x";
	  aaa="y";
	}
        if (readdir==1){
	  aa="z";
	  ccc="x";
	}
        if (readdir==2){
	  aa="y";
	  bbb="x";
	}
        if (readdir==3){ 
	  aa="x";
	  aaa="x";
	}
	cout<<"slcdir="<<cc<<"; phasedir="<<bb<<"; readdir="<<aa<<endl;
      //Setting up the the temporary volumes and variables needed for the homodyne recon
	int ck=33;
        if (nphase%2==0) {
	   ck=nphase/2+1;//central line of the k-space
	  } else {
	   ck=(nphase-1)/2+1;
        }
        if (verbose.value()) {
	  cout<<"Central line of the k-space is "<<ck<<endl;
	  cout<<"Starting k-space line is "<<startkspace<<endl;
	}
        int nslope=ck-startkspace+1;
	int npWslope=nslope*2+1;
	volume<float> W(nx,ny,nz);
	volume<float> Ws(nx,ny,nz);
	W=0;Ws=0;
	int np1=startkspace-2;
	for (int s=0;s<nslc;s++){
	  for (int k=0;k<np1;k++){
	    for (int l=0;l<nread;l++){
	      int zdir=s;
	      int ydir=k;
	      int xdir=l;
	      setdir(xdir, ydir, zdir, l,k,s, readdir, phasedir, slcdir);
	      W(xdir,ydir,zdir)=0;Ws(xdir,ydir,zdir)=0;
	    }
	  }
	}
	for (int s=0;s<nslc;s++){
	  for (int k=np1;k<npWslope+np1;k++){
	    for (int l=0;l<nread;l++){
	      int zdir=s;
	      int ydir=k;
	      int xdir=l;
	      setdir(xdir, ydir, zdir, l,k,s, readdir, phasedir, slcdir);
	      W(xdir,ydir,zdir)=(k-np1)*2.0/(nslope*2.0);
	      //W(xdir,ydir,zdir)=Wslope(xdir,ydir-np1,zdir);
	      Ws(xdir,ydir,zdir)=1;
	    }
	  }
	}
	for (int s=0;s<nslc;s++){
	  for (int k=np1+npWslope;k<nphase;k++){
	    for (int l=0;l<nread;l++){
	      int zdir=s;
	      int ydir=k;
	      int xdir=l;
	      setdir(xdir, ydir, zdir, l,k,s, readdir, phasedir, slcdir);
	      W(xdir,ydir,zdir)=2;Ws(xdir,ydir,zdir)=0;
	    }
	  }
	}
	if (verbose.value()){
           save_volume(W,outname.value()+"_W");
           save_volume(Ws,outname.value()+"_Ws");
	   cout<<"Bottom "<<np1<<" lines are zeros. Middle "<<npWslope<<" lines are ones (in case of Ws) or slope increase (in case of W). The top  "<<nphase-np1-npWslope<<" lines are zeros (Ws) or two (W)."<<endl;
	}
	for (int nn=1;nn<=n;nn++){
	  volume<float> kspaceHpc_real(nx,ny,nz);
	  volume<float> kspaceHpc_imag(nx,ny,nz);
	  for (int s=0;s<nslc;s++){
	    for (int k=0;k<nphase;k++){
	      for (int l=0;l<nread;l++){
		int zdir=s;
		int ydir=k;
		int xdir=l;
		setdir(xdir, ydir, zdir, l,k,s, readdir, phasedir, slcdir);
		kspace_real(xdir,ydir,zdir,nn-1)=kspace_real(xdir,ydir,zdir,nn-1)*W(xdir,ydir,zdir);
		kspace_imag(xdir,ydir,zdir,nn-1)=kspace_imag(xdir,ydir,zdir,nn-1)*W(xdir,ydir,zdir);
		kspaceHpc_real(xdir,ydir,zdir)=kspace_real(xdir,ydir,zdir,nn-1)*Ws(xdir,ydir,zdir);
		kspaceHpc_imag(xdir,ydir,zdir)=kspace_imag(xdir,ydir,zdir,nn-1)*Ws(xdir,ydir,zdir);
	      }
	    }
	  }
	  kspace_real[nn-1].swapdimensions(aa,bb,cc); 
	  kspace_imag[nn-1].swapdimensions(aa,bb,cc);
	  fftshift(kspace_real[nn-1]);
	  fftshift(kspace_imag[nn-1]);
	  fft2(kspace_real[nn-1],kspace_imag[nn-1]);
	  fftshift(kspace_real[nn-1]);
	  fftshift(kspace_imag[nn-1]);
	  kspaceHpc_real.swapdimensions(aa,bb,cc); 
	  kspaceHpc_imag.swapdimensions(aa,bb,cc);
	  fftshift(kspaceHpc_real);
	  fftshift(kspaceHpc_imag);
	  fft2(kspaceHpc_real,kspaceHpc_imag);
	  fftshift(kspaceHpc_real);
	  fftshift(kspaceHpc_imag);
	  volume<float> pcA(nx,ny,nz);
	  volume<float> pcB(nx,ny,nz);
	  for (int s=0;s<nslc;s++){
	    for (int k=np1;k<nphase;k++){
	      for (int l=0;l<nread;l++){
		int zdir=s;
		int ydir=k;
		int xdir=l;
		setdir(xdir, ydir, zdir, l,k,s, readdir, phasedir, slcdir);
		float tmp1=1;
		float tmp2=M_PI;
		if (kspaceHpc_real(xdir,ydir,zdir)>0.0000001){
		  tmp1=kspaceHpc_imag(xdir,ydir,zdir)/kspaceHpc_real(xdir,ydir,zdir);
		  tmp2=atan(tmp1);
		} else{
		  if (kspaceHpc_imag(xdir,ydir,zdir)>0){
		    tmp2=M_PI;
		  }else{
		    tmp2=-M_PI;
		  }
		}
		pcA(xdir,ydir,zdir)=cos(tmp2);
		pcB(xdir,ydir,zdir)=-sin(tmp2);
	      }
	    }
	  }
	  for (int s=0;s<nslc;s++){
	    for (int k=np1;k<nphase;k++){
	      for (int l=0;l<nread;l++){
		int zdir=s;
		int ydir=k;
		int xdir=l;
		setdir(xdir, ydir, zdir, l,k,s, readdir, phasedir, slcdir);
		kspace_real(xdir,ydir,zdir,nn-1)=pcA(xdir,ydir,zdir)*kspace_real(xdir,ydir,zdir,nn-1)-pcB(xdir,ydir,zdir)*kspace_imag(xdir,ydir,zdir,nn-1);
	      }
	    }
	  }
	  kspace_real[nn-1].swapdimensions(aaa,bbb,ccc); 
	  //kspace_real[nn-1].swapdimensions("-x","-y","z"); 
	}
	save_volume4D(kspace_real,outname.value()+"_homo");
    } else {
      for (int nn=1;nn<=n;nn++){
        string aa="x";
        string bb="y";
        string cc="z";
        string aaa="x";
        string bbb="y";
        string ccc="z";
	//in the old version in order to make it be the same orientation as the images from the scanner I had to do swapdimensions("-x","-y","z") after the I did the fft2. The thing is the convention for the scanner is (y,x,z) and for me was (x,y,z) so maybe that had to do.will see...still testing this orientation thing. 
        if (abs(slcdir)==1){
          cc="z";
	  ccc="z";
	}
	if (abs(slcdir)==2){
          cc="y";
	  bbb="z";
	}
        if (abs(slcdir)==3){
	  cc="x";
	  aaa="z";
	}
        if (phasedir==1){
	  bb="z";
	  ccc="y";
	}
        if (phasedir==2){
	  bb="y";
	  bbb="y";
	}
        if (phasedir==3){
	  bb="x";
	  aaa="y";
	}
        if (readdir==1){
	  aa="z";
	  ccc="x";
	}
        if (readdir==2){
	  aa="y";
	  bbb="x";
	}
        if (readdir==3){ 
	  aa="x";
	  aaa="x";
	}
	cout<<"slcdir="<<cc<<"; phasedir="<<bb<<"; readdir="<<aa<<endl;
        kspace_real[nn-1].swapdimensions(aa,bb,cc); 
        kspace_imag[nn-1].swapdimensions(aa,bb,cc);
        fftshift(kspace_real[nn-1]);
        fftshift(kspace_imag[nn-1]);
        fft2(kspace_real[nn-1],kspace_imag[nn-1]);
        // WARNING: from now on kspace is actually IMAGE SPACE!
        fftshift(kspace_real[nn-1]);
        fftshift(kspace_imag[nn-1]);
        kspace_real[nn-1].swapdimensions(aaa,bbb,ccc); 
        kspace_imag[nn-1].swapdimensions(aaa,bbb,ccc);
        //kspace_real[nn-1].swapdimensions("-x","-y","z"); 
        //kspace_imag[nn-1].swapdimensions("-x","-y","z");
      }
      if (useabs.value()) {
        volume4D<double> dummy(kspace_real);
        volume4D<double> dummy_phase(kspace_real);
        for (int nn=0; nn<kspace_real.tsize(); nn++) { 
	  dummy[nn] = sqrt(kspace_real[nn]*kspace_real[nn] + kspace_imag[nn]*kspace_imag[nn]);
          for (int z=kspace_real.minz(); z<=kspace_real.maxz(); z++) {
	    for (int y=kspace_real.miny(); y<=kspace_real.maxy(); y++) {
	      for (int x=kspace_real.minx(); x<=kspace_real.maxx(); x++) {
	        dummy_phase(x,y,z,nn) = atan2(kspace_real(x,y,z,nn),kspace_imag(x,y,z,nn));
	      }
	    }
	  }
	}
        save_volume4D(dummy,outname.value()+"_abs");
        save_volume4D(dummy_phase,outname.value()+"_phase");
      } else {
        save_volume4D(kspace_real,outname.value()+"_real");
        save_volume4D(kspace_imag,outname.value()+"_imag");
        //volume4D<float> kspace_real_float,kspace_imag_float;
        //copyconvert(kspace_real,kspace_real_float);
        //copyconvert(kspace_imag,kspace_imag_float);
        //save_complexvolume4D(kspace_real_float,kspace_imag_float,outname.value());
      }
    }
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{
  Tracer tr("main");
  OptionParser options(title, examples);
  try {
    // must include all wanted options here (the order determines how
    //  the help message is printed)
    options.add(inname);
    options.add(outname);
    options.add(opt_kcoord);
    options.add(koutname);
    options.add(useabs);
    options.add(opt_homo);
    options.add(verbose);
    options.add(help);
    options.add(opt_pulse);

//Tejas-apod
	options.add(opt_doapod);
	options.add(opt_cutoff);
	options.add(opt_rolloff);
	//Debug
	options.add(opt_savewindow);
//Tejas-end

    nonoptarg = options.parse_command_line(argc, argv);
    // line below stops the program if the help was requested or 
    //  a compulsory option was not set
    if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
      {
	options.usage();
	exit(EXIT_FAILURE);
      }    
  }  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } 
  // Call the local functions
  return do_work(argc,argv);
}

