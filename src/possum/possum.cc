/*  POSSUM
    Ivana Drobnjak, Mark Jenkinson and Matthew Webster
    Copyright (C) 2005-2010 University of Oxford  */

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

// POSSUM

#include <iostream>
#include <string>
#include <fstream>
#include <unistd.h>
#include <time.h>

#ifdef USE_MPI
#include <mpi.h>
#include <unistd.h>
#endif //USE_MPI

#include "libprob.h"
#include "newmatap.h"
#include "newmatio.h"
#include "newimage/newimageall.h"
#include "possumfns.h"
#include "utils/options.h"
#include "newimage/costfns.h"
#include "miscmaths/miscmaths.h"

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

using namespace NEWIMAGE;
using namespace NEWMAT;
using namespace MISCMATHS;
using namespace Utilities;
//using namespace std;

string title="possum\nCopyright(c) 2007, University of Oxford (Ivana Drobnjak)";
string examples="possum -i <input phantom volume> -x <MR parameters matrix> -p <pulse> -f <RF slice profile> -m <motion file> -o <output signal matrix> [optional arguments]";

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);

//INPUT object and its characteristics (including susc effects on B0 and the RF inhomogeneities)
Option<string> opt_object(string("-i,--inp"), string(""),
		  string("<input4Dvol-filename> (Input object)"),
		  true, requires_argument);
Option<string> opt_tissue(string("-x,--mrpar"),string("") ,
		  string("<inputmatrix-filename> (MR parameters)"),
		  true,requires_argument);
Option<string> opt_b0(string("-b,--b0p"), string(""),
		  string("<input3Dvol-basename> (B0 inhomogeneities due to the susceptibility differences - base name, without extras z_dz, z_dx etc)"),
                  false, requires_argument);
Option<string> opt_b0extra(string("--b0extra"), string(""),
		 string("<input3Dvol-filename> (B0 inhomogeneities due to an extra field - see b0time)"),
		  false, requires_argument);
Option<string> opt_b0timecourse4D(string("--b0time"),string(""),
	          string("<inputmatrix-filename> (B0inhomogeneities_timecourse [time(s) multiply_factor(perc 0 to 1)] - see b0extra) "),
	          false,requires_argument);
Option<string> opt_RFrec(string("-r,--rfr"), string(""),
		  string("<input3Dvol-filename> ( RF inhomogeneity - receive. NOTE: not yet to be used ) "),
		  false, requires_argument);
Option<string> opt_RFtrans(string("-s,--rft"), string(""),
		  string("<input3Dvol-filename> ( RF inhomogeneity - transmit. NOTE: not yet to be used )"),
		  false, requires_argument);

//INPUT motion and activation
Option<string> opt_activation4D(string("-q,--activ4D"),string(""),
				string("<input4Dvol-filename> (Activation volume) "),
	          false,requires_argument);
Option<string> opt_timecourse4D(string("-u,--activt4D"),string(""),
	          string("<inputmatrix-filename> (Activation4D_timecourse [time(s)])"),
	          false,requires_argument);
Option<string> opt_activation(string("-a,--activ"),string(""),
	          string("<input3Dvol-filename> (Activation volume)"),
	          false,requires_argument);
Option<string> opt_timecourse(string("-t,--activt"),string(""),
	          string("<inputmatrix-filename> (Activation_timecourse [time(s) multiply_factor(perc 0 to 1)])"),
	          false,requires_argument);
Option<string> opt_motion(string("-m,--motion"), string(""),
		  string("<inputmatrix-filename> (Motion matrix [time(s) Tx(m) Ty(m) Tz(m) Rx(rad) Ry(rad) Rz(rad)]) "),
		  true, requires_argument);

//INPUT for the pulse sequence
Option<string> opt_pulse(string("-p,--pulse"), string(""),
		  string("<inputmatrix-basename> (Pulse sequence - all additional files .posx,.posy, etc,  expected to be in the same directory)"),
		  true, requires_argument);

Option<string> opt_slcprof(string("-f,--slcprof"), string(""),
		  string("<inputmatrix-filename> (RF slice profile)"),
		  true, requires_argument);

Option<bool> opt_rfavg(string("--rfavg"), false,
           string("If this option is ON it will use RF angle averging"),
                       false, no_argument);

//INPUT for the computational efficiency 
Option<int>    opt_level(string("-l,--lev"), 1,
		  string("{1,2,3,4} (Levels: 1.no motion//basic B0 2.motion//basic B0, 3.motion//full B0, 4.no motion//time changing B0)"),
		  false,requires_argument);
Option<bool> opt_nospeedup(string("--nospeedup"), false, 
		     string("If this option is ON it will NOT do the speedup but will do signal for all the slices for each voxel."), 
		     false, no_argument);

//INPUT for the manual paralelisation -- used with sge_possum
Option<int>    opt_nproc(string("--nproc"), 1,
		  string("<int> (INPUT for the paralelisation -- Number of processors we have available)"),
		  false,requires_argument);
Option<int>    opt_procid(string("--procid"), 0,
		  string("<int> (INPUT for the paralelisation -- ID of the processor we are on)"),
		  false,requires_argument);

//OUTPUT signal
Option<string> opt_signal(string("-o,--out"), string(""),
		  string("<outputmatrix-filename> (Signal - [sreal, simag])"),
		  true, requires_argument);

//INPUT main event matrix
Option<string> opt_mainmatrix(string("-e,--mainmatx"), string(""),
		  string("<inputmatrix-filename> (Main event matrix [t(s),rf_ang(rad),rf_freq_band(Hz),(4)=rf_cent_freq(Hz),read(1/0),Gx,Gy,Gz(T/m),Tx,Ty,Tz(m),angle_of_rot B(rad),rot_axis Bx,By,Bz(m),angle_of_rot A(rad),rot_axis Ax,Ay,Az(m)]) "),
		  true, requires_argument);

//OUTPUT kcoord if needed
Option<bool> opt_kcoord(string("-k,--kcoord"), false,
		  string("If this option is ON it will save the kspace coordinates"),
		  false, no_argument);


int nonoptarg;

/////////////////////////////////////////////////////////////////////////////////////////////////////
int compute_volume(int argc, char *argv[])
{
  cout<<"Starting POSSUM..."<<endl;
  cout<<""<<endl;

//tejas-22.11.12
//start time
time_t theTime = time(NULL);
cout << "possum-start-time\t"<< opt_procid.value() << "\t" << theTime << endl;
//end

  /////////////////////////////////////////////////////////////////////////////
  // SET UP COORDINATE SYSTEM with the CENTER IN THE CENTER OF THE OBJECT    
  /////////////////////////////////////////////////////////////////////////////
  cout<<"Setting up the coordinate system..."<<endl;
  RowVector posx; 
  RowVector posy; 
  RowVector posz;
  posx=read_ascii_matrix(opt_pulse.value()+".posx");//in SI already
  posy=read_ascii_matrix(opt_pulse.value()+".posy");//in SI already
  posz=read_ascii_matrix(opt_pulse.value()+".posz");//in SI already
  int Nxx=posx.Ncols();
  cout<<""<<endl;
  //////////////////////////////////////////////////////////////////////////
  //PARALELL STUFF
  //////////////////////////////////////////////////////////////////////////
  cout<<"Setting up the process ID..."<<endl;
  //int lo=0;
  //int hi=Nxx;
  int myid=0;
  int numprocs=1;
  #ifdef USE_MPI
    MPI_Init(&argc,&argv);//MPI::Init(argc,argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);//int numprocs=MPI::COMM_WORLD.Get_size();
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);//myid=MPI::COMM_WORLD.Get_rank();
    //lo=myid*(hi/numprocs);
    //if (myid<numprocs-1) hi=(myid+1)*(hi/numprocs);
    //cout<<" Process "<<myid<<" does from "<<lo<<" to "<<hi<<".\n";
  #else 
    myid=opt_procid.value();
    numprocs=opt_nproc.value();
  #endif //USE_MPI
  cout<<"Number of processors is "<<numprocs<<".; My ID is "<<myid<<"."<<endl;
  cout<<""<<endl;
  //////////////////////////////////////////////////////////////////////////
  // READ IN THE OBJECT (BRAIN) and TISSUE PROPERTIES                                       
  //////////////////////////////////////////////////////////////////////////
  sleep(myid*10);
  cout<<"Reading the object for ID "<<myid<<"..."<<endl;
  volume4D<double> phantom;//consists of gry,wht,csf,fat,mus,con,gli,skn (in that order)
  read_volume4DROI(phantom,opt_object.value(),myid,0,0,0,Nxx,-1,-1,-1,numprocs,1,1,1);
  //read_volume4D(phantom,opt_object.value());
  int Nx=phantom.xsize();double xdim=phantom.xdim()*0.001; 
  int Ny=phantom.ysize();double ydim=phantom.ydim()*0.001;
  int Nz=phantom.zsize();double zdim=phantom.zdim()*0.001;
  int Nt=phantom.tsize();
  print_volume_info(phantom,"object");
  cout<<""<<endl;
  cout<<"Reading the tissue properties..."<<endl;
  Matrix tissue;  //tissue caracteristics   T1,T2,PD,ChemicalShift=value(ppm)*gammabar*B0
  tissue=read_ascii_matrix(opt_tissue.value());//in SI already
  cout<<"T1,T2,SD,CS: "<<endl;
  cout<<tissue<<endl;
  ///////////////////////////////////////////////////////
  //PULSE & MOTION MATRIX SORT IN MAINMATRIX
  ///////////////////////////////////////////////////////
  RowVector pulseinfo;
  pulseinfo=read_ascii_matrix(opt_pulse.value()+".info");//[SeqType,TE,TR,TRslc,Nx,Ny,dx,dy,maxG,RiseT,BWrec, Nvol,Nslc,SlcThk,SlcDir,Gap,zstart,FlipAngle]
  cout<<"[SeqType,TE,TR,TRslc,Nx,Ny,dx,dy,maxG,RiseT,BW,Nvol,Nslc,SlcThk,SlcDir,Gap,zstart,FA]"<<endl;
  cout<<pulseinfo<<endl;
  cout<<""<<endl;
  cout<<"Reading the motion file..."<<endl;
  Matrix motion;
  motion=read_ascii_matrix(opt_motion.value());
  //cout<<"Motion file is "<<motion<<endl;
  /////////////////////////////////////////////////////////////////////////
  //SLICE PROFILE
  /////////////////////////////////////////////////////////////////////////
  cout<<"Creating table for the slice profile..."<<endl;
  double* table_slcprof;
  //slcprof has two columns: first one is normalised frequency (f-fc)/df -user defined, second is the normalised strength (0 1)- amount of the flip angle 
  Matrix slcprof=read_ascii_matrix(opt_slcprof.value());
  int Nslc=slcprof.Nrows();
  table_slcprof = new double[Nslc];
  for (int n=0;n<=Nslc-1;n++){
    table_slcprof[n]=slcprof(n+1,2);
  }
  double slcpr_range=slcprof(Nslc,1)-slcprof(1,1);
  double dslcp=slcpr_range/(Nslc-1);//step size for table for the slcprof
  double dslcp_first=slcprof(1,1);
  double slcthk=pulseinfo(14);//slcthk(m)
  double slcpr_add=(slcpr_range-1)/2*slcthk;
  cout<<"Stepsize for the table: dslcp= "<<dslcp<<endl;
  cout<<""<<endl;
  //////////////////////////////////////////////////////////////////
  //CALCULATING ZSTART, ZEND AND LEVEL, useful just for the speed up
  //////////////////////////////////////////////////////////////////
  cout<<"Calculating zstart, zend and level..."<<endl;
  int numvol=(int) (pulseinfo(12));
  int numslc=(int) (pulseinfo(13));
  double gap=pulseinfo(16);// gap (m)
  int resX=(int) (pulseinfo(5));
  int resY=(int) (pulseinfo(6));
  int kspacestart=(int) (pulseinfo(22));
  int nrf=numslc*numvol;//number of rf pulses 
  int motionsize=motion.Nrows();
  double tzmax, tzmin, txmax, txmin, tymax,tymin, rxmaxabs, rymaxabs, rzmaxabs; 
  tzmax=0;tzmin=0; txmax=0;txmin=0; tymax=0; tymin=0; rxmaxabs=0;  rymaxabs=0; rzmaxabs=0;
  int level;
  for (int k=1;k<=motionsize;k++){
    if (motion(k,2)>txmax) txmax=motion(k,2);
    if (motion(k,2)<txmin) txmin=motion(k,2);
    if (motion(k,3)>tymax) tymax=motion(k,3);
    if (motion(k,3)<tymin) tymin=motion(k,3);
    if (motion(k,4)>tzmax) tzmax=motion(k,4);
    if (motion(k,4)<tzmin) tzmin=motion(k,4);
    if (fabs(motion(k,5))>rxmaxabs) rxmaxabs=fabs(motion(k,5));
    if (fabs(motion(k,6))>rymaxabs) rymaxabs=fabs(motion(k,6));
    if (fabs(motion(k,7))>rzmaxabs) rzmaxabs=fabs(motion(k,7));
  }
  //slc selection direction stuff
  int slcdir=(int) (pulseinfo(15));
  cout<<"Slice selection direction is "<<slcdir<<". (-+1 is z, -+2 is y and -+3 is x.)"<<endl;
  double ssrxmaxabs,ssrymaxabs, sszdim, sstzmax,sstzmin;
  int ssNz;
  if (fabs((float)slcdir)==1){
    sszdim=zdim;
    ssNz=Nz;
    ssrxmaxabs=rxmaxabs;
    ssrymaxabs=rymaxabs;
    sstzmax=tzmax;
    sstzmin=tzmin;
  }
  else if (fabs((float)slcdir)==2){
    sszdim=ydim;
    ssNz=Ny;
    ssrxmaxabs=rxmaxabs;
    ssrymaxabs=rzmaxabs;
    sstzmax=tymax;
    sstzmin=tymin;
  }
  else {
    sszdim=xdim;
    ssNz=Nx;
    ssrxmaxabs=rymaxabs;
    ssrymaxabs=rzmaxabs;
    sstzmax=txmax;
    sstzmin=txmin;
  }
  //level calculation
  if (opt_level.set()) level=opt_level.value();
  else {
    if (tzmax==0 && tzmin==0 && txmax==0 && txmin==0 && tymax==0 && tymin==0 && rxmaxabs==0 && rymaxabs==0 && rzmaxabs==0) level=1;
    else if ( ssrxmaxabs==0 && ssrymaxabs==0 ) level=2;
    else level=3;
  }
  cout<<"Level is "<<level<<endl;
  cout<<""<<endl;
  cout<<"Extra slc calculation..."<<endl;
    //Basic
    double sszstart_in=pulseinfo(17);// starting point in the volume in the direction of the slice selection (used to be just z)
    double sszend_in=sszstart_in+numslc*slcthk+(numslc-1)*gap;
    double nvox=1/sszdim;//number of voxels per 1m in slice selection direction
    int sszstart_p=(int) (sszstart_in*nvox+0.0001);//zstart in vox, 0.0001 is just a fix so that it rounds it properly 
    int sszend_p=(int)ceil(sszend_in*nvox);
    cout<<"Specified starting point is zstart= "<<sszstart_in<<" (in m) = "<<sszstart_p<<" (in voxels)"<<endl;
    cout<<"Specified end point is zend= "<<sszend_in<<" (in m) = "<<sszend_p<<" (in voxels)"<<endl;
    //Slice profile add ons
    int slcpr_add_vox=(int)ceil(slcpr_add*nvox);
    cout<<"Slice profile add-on up and down is "<<slcpr_add<<" (in m) = "<<slcpr_add_vox<<" (in voxels)"<<endl;
    //Motion add ons
    double rotation_add=max(fabs(posx(1))*sin(rymaxabs),fabs(posy(1))*sin(rxmaxabs));
    double motion_add_down=fabs(sstzmax)+rotation_add;
    double motion_add_up=fabs(sstzmin)+rotation_add;
    int motion_add_down_vox=(int)ceil(motion_add_down*nvox);
    int motion_add_up_vox=(int)ceil(motion_add_up*nvox);
    cout<<"Motion add-on up is "<<motion_add_up<<" (in m) = "<<motion_add_up_vox<<" (in voxels)"<<endl;
    cout<<"Motion add-on down is "<<motion_add_down<<" (in m) = "<<motion_add_down_vox<<" (in voxels)"<<endl;
    //int extra_down_slc=(int)(ceil((motion_add_down+slcpr_add)/slcthk));
    int sszstart=sszstart_p-slcpr_add_vox-motion_add_down_vox;
    if (sszstart<0) sszstart=0; 
    //int extra_up_slc=(int)(ceil((motion_add_up+slcpr_add)/slcthk));
    int sszend=sszend_p+slcpr_add_vox+motion_add_up_vox;  
    if (sszend>(ssNz-1)) sszend=ssNz-1;
    cout<<"FINAL: Begining of the object is "<<sszstart<<" and the end is "<<sszend<<" (in vox, in slice select direction)."<<endl;
  int nospeedup=0;
  if (opt_nospeedup.value())nospeedup=1;//when no speed up used for slices
  //counter endings for the main loop
  int xstart=0;
  int xend=Nx;
  int ystart=0;
  int yend=Ny;
  int zstart=0;
  int zend=Nz;
  int sszz;
  int sszz_slc=1; 
  if (fabs((float)slcdir)==1){
       zstart=sszstart;
       zend=sszend+1;
    }
    else if (fabs((float)slcdir)==2){
       ystart=sszstart;
       yend=sszend+1;
    }
    else {
       xstart=sszstart;
       xend=sszend+1;
    }
  cout<<""<<endl;
  /////////////////////////////////////////////////////////////////////
  // ACTIVATION (T2*)                
  /////////////////////////////////////////////////////////////////////
  cout<<"Creating activation4D volume and timecourse vector..."<<endl;
  volume4D<double> activation4D;
  volume<double> activation(Nx,Ny,Nz);
  double* timecourse;
  double* timecourse_2=0;
  double* activation4D_voxel;
  int Nact;
  if (opt_activation.set()) {
    cout<<"3D activation mode"<<endl;
    read_volumeROI(activation,opt_activation.value(),myid,0,0,Nxx,-1,-1,numprocs,1,1);
    //read_volume(activation,opt_activation.value());
    Matrix timecourse_tmp;
    timecourse_tmp=read_ascii_matrix(opt_timecourse.value());//multiply with the activation file. 
    //cout<<"Timecourse_matrix is "<<timecourse_tmp<<endl;
    Nact=timecourse_tmp.Nrows();
    timecourse=new double[Nact];
    timecourse_2=new double[Nact];
    activation4D_voxel=new double[Nact];
    cout<<"Testing if the values in the activation volumes are smaller than the T2* values"<<endl;
    for (int n=0;n<=Nact-1;n++){
      cout<<"vol= "<<n<<endl;
      timecourse[n]=timecourse_tmp(n+1,1);
      timecourse_2[n]=timecourse_tmp(n+1,2);
      for (int xx=0;xx<Nx;xx++){
        for (int yy=0;yy<Ny;yy++){
          for (int zz=0;zz<Nz;zz++){
            for (int tt=0;tt<Nt;tt++){
              double test=fabs(activation(xx,yy,zz)*timecourse_2[n]);
              if (test>=tissue(tt+1,2)) {
		cout<<"WARNING: Perturbation in the T2* values in the input activation file bigger than the input T2* values."<<endl;
		exit(EXIT_FAILURE);
	      }
	    }
	  }
	}
      }
    }
  }
  else  if (opt_activation4D.set()) {
    cout<<"4D activation mode"<<endl;
    read_volume4DROI(activation4D,opt_activation4D.value(),myid,0,0,0,Nxx,-1,-1,-1,numprocs,1,1,1);
    //read_volume4D(activation4D,opt_activation4D.value());
    Matrix timecourse_tmp;
    timecourse_tmp=read_ascii_matrix(opt_timecourse4D.value());
    Nact=timecourse_tmp.Nrows();
    timecourse=new double[Nact];
    activation4D_voxel=new double[Nact];
    cout<<"Testing if the values in the activation volumes are smaller than the T2* values"<<endl;
    for (int n=0;n<=Nact-1;n++){
      cout<<"vol= "<<n<<endl;
      timecourse[n]=timecourse_tmp(n+1,1);
      for (int xx=0;xx<Nx;xx++){
        for (int yy=0;yy<Ny;yy++){
          for (int zz=0;zz<Nz;zz++){
            for (int tt=0;tt<Nt;tt++){
              double test=activation4D(xx,yy,zz,n);
              if (test>=tissue(tt+1,2)) {
		cout<<"WARNING: Perturbation in the T2* values bigger than the T2* values."<<endl;
		exit(EXIT_FAILURE);
	      }
	    }
	  }
	}
      }
    }             
  }
  else { 
    cout<<"No activation mode"<<endl;
    Nact=1;
    activation4D_voxel=new double[Nact];
    timecourse=new double[Nact];
    timecourse[0]=0.0;
    activation4D_voxel[0]=0.0;
  }
  cout<<""<<endl;
  /////////////////////////////////////////////////////////////////////////
  //RF INHOMOGENEITY
  /////////////////////////////////////////////////////////////////////////
  cout<<"Reading the RF inhomogeneity volumes..."<<endl;
  volume<double> RFrec(Nx,Ny,Nz);//signal=signal*RFrec
  volume<double> RFtrans(Nx,Ny,Nz);//flip_ang=flip_ang*RFtrans
  //if (opt_RFrec.set()) read_volume_new(RFrec,opt_RFrec.value(),myid,numprocs,Nxx);
  if (opt_RFrec.set()) read_volumeROI(RFrec,opt_RFrec.value(),myid,0,0,Nxx,-1,-1,numprocs,1,1);
  else {
    for (int px=0;px<Nx;px++){
      for (int py=0;py<Ny;py++){
        for (int pz=0;pz<Nz;pz++){
          RFrec(px,py,pz)=1;
	}
      }
    }
  }
  //if (opt_RFtrans.set())  read_volume_new(RFtrans,opt_RFtrans.value(),myid,numprocs,Nxx);
  if (opt_RFtrans.set())  read_volumeROI(RFtrans,opt_RFtrans.value(),myid,0,0,Nxx,-1,-1,numprocs,1,1);
  else {
    for (int px=0;px<Nx;px++){
      for (int py=0;py<Ny;py++){
        for (int pz=0;pz<Nz;pz++){
          RFtrans(px,py,pz)=1;
	}
      }
    }
  }
  print_volume_info(RFrec,"RFobject");
  //cout<<Nx<<" "<<Ny<<" "<<Nz<<"; RFrec(2,2,20)= "<<RFrec(2,2,20)<<"; RFrec(2,2,520)= "<<RFrec(2,2,520)<<endl;
  cout<<""<<endl;
  ///////////////////////////////////////////////////////////////////////
  //TESTING
  ///////////////////////////////////////////////////////////////////////
  int opt_test=0;
  if (verbose.value()) {
    opt_test=1;
    cout<<"Verbose is ON"<<endl;
  }
  ///////////////////////////////////////////////////////////////////////
  //K-SPACE COORDINATES
  ///////////////////////////////////////////////////////////////////////
  int save_kcoord=0;
  if (opt_kcoord.value()){ 
    save_kcoord=1;
    cout<<"k-space coordinates will be saved"<<endl;
  } 
  cout<<""<<endl;
  /////////////////////////////////////////////////////////////////////////
  // SIGNAL                                                             
  /////////////////////////////////////////////////////////////////////////
  cout<<"Calculating the signal..."<<endl;
  int nreadp=resX*(resY-kspacestart+1)*numslc*numvol;
  cout<<"Number of read out points is "<<nreadp<<endl;
  Matrix signal(2,nreadp);//two rows, one real and one complex for the signal in sum
  signal=0;
  double* sreal; 
  double* simag;
  sreal= new double[nreadp];
  simag= new double[nreadp];
  for (int i=0;i<=nreadp-1;i++){
    sreal[i]=0.0;
    simag[i]=0.0;
  }
  int voxelcounter=0;
  double cxyz=xdim*ydim*zdim;
  string outputname=opt_signal.value();
  cout<<""<<endl;

	  ///////////////////////////////////////////////////////////
	  //NO MOTION/
	  ///////////////////////////////////////////////////////////
	  if (level==1){
	    cout<<"Reading the pulse sequence..."<<endl;
	    PMatrix pulse;
	    read_binary_matrix(pulse, opt_pulse.value());
	    ////////////////////////
	    // B0 PERTURBATION 
	    ////////////////////////
	    cout<<"LEVEL1"<<endl;
	    cout<<"Creating 1 B0file together with 3 gradient files..."<<endl;
	    volume<double> b0(Nxx,Ny,Nz);
	    volume<double> b0x(Nxx,Ny,Nz);
	    volume<double> b0y(Nxx,Ny,Nz);
	    volume<double> b0z(Nxx,Ny,Nz);
	    if (opt_b0.set()) {
	      read_volume(b0,opt_b0.value()+"z_dz");
	      calc_gradientsROI(b0,b0x,b0y,b0z,myid,Nxx,numprocs);
	    }
	    else {
	      b0=phantom[0]*0;
	      b0x=b0;
	      b0y=b0;
	      b0z=b0;
	    }  
	    print_volume_info(b0,"b0");
	    cout<<""<<endl;
	    //save_volume(b0z,"b0ztest");
	    /////////////
	    //MAIN LOOP
	    /////////////
	    cout<<"Main loop..."<<endl;
	    for (register int tt=0;tt<Nt;tt++){
	      for (register int zz=zstart;zz<zend;zz++){
		cout<<"Tissue type="<<tt<<"; zstart="<<zstart<<"; zz="<<zz<<"; zend="<<zend<<"; Voxelnumber="<<voxelcounter<<endl;
		for (register int yy=ystart;yy<yend;yy++){
		  int xxx=myid+xstart;
		  for (register int xx=xstart;xx<xend;xx++){
		    //slice speed up stuff
		    if (fabs((float)slcdir)==1) sszz=zz;
		    else if (fabs((float)slcdir)==2) sszz=yy;
		    else sszz=xx;
	      	    float slctmp=(sszz-sszstart_p)/(slcthk*nvox);
		    if(ceil(slctmp) == slctmp) {
		      sszz_slc=(int)(ceil(slctmp)+1);
		    } else {
		      sszz_slc=(int) (ceil(slctmp));
		    }
		    if (phantom(xx,yy,zz,tt)!=0){
		      voxelcounter=voxelcounter+1;
		      double den=phantom(xx,yy,zz,tt)*RFrec(xx,yy,zz)*cxyz;
		      if (opt_activation.set()) {
		        for (int n=0;n<=Nact-1;n++){
			  double a=activation(xx,yy,zz)*timecourse_2[n];
			  double b=tissue(tt+1,2);
			  activation4D_voxel[n]=a*b/(a+b);// conversion of beta into beta1 because of the integral (see possum no 7 page 121)
			}
		      }
		     else if (opt_activation4D.set()){
		       for (int n=0;n<=Nact-1;n++){
			 double a=activation4D(xx,yy,zz,n);//zz-zstart_p when having only a pieace of the activation volume so we start from where the phantom starts 
		         double b=tissue(tt+1,2);
		         activation4D_voxel[n]=a*b/(a+b);
		       }              
		     }
		      voxel1(posx(xxx+1),posy(yy+1),posz(zz+1),tissue.Row(tt+1),
		             pulse,nreadp,voxelcounter,xdim,ydim,zdim,
		             b0(xx,yy,zz),b0x(xx,yy,zz),b0y(xx,yy,zz),b0z(xx,yy,zz),
		             timecourse,activation4D_voxel,Nact,outputname,
		             table_slcprof,dslcp,dslcp_first,Nslc,den,RFtrans(xx,yy,zz),
		             opt_test,nospeedup,save_kcoord,sreal,simag);   
		      //srealT+=sreal*coil(xx,yy,zz); //Ivana 01.11.12 -trial		      
		    }
		  xxx=xxx+numprocs;
		  }
		}
	      }
	    }
	  }
	  ////////////////////////////////////////////////////////////
	  //MOTION WHEN ONLY POSSIBLE ROTATION CAN BE IN PLANE
	  ////////////////////////////////////////////////////////////
	  if (level==2){
	    cout<<"LEVEL2"<<endl;
	    ///////////////////////
	    //MAIN MATRIX
	    ///////////////////////
	    cout<<"Sorting pulse sequence and the motion matrix in one large matrix..."<<endl;
	    cout<<"Reading the pulse sequence..."<<endl;
	    PMatrix pulse;
	    read_binary_matrix(pulse,opt_mainmatrix.value());
	    ////////////////////////
	    // B0 PERTURBATION 
	    ////////////////////////
	    cout<<"Creating 1 B0file together with 3 gradient files..."<<endl;
	    volume<double> b0(Nxx,Ny,Nz);
	    volume<double> b0x(Nxx,Ny,Nz);
	    volume<double> b0y(Nxx,Ny,Nz);
	    volume<double> b0z(Nxx,Ny,Nz);
	    if (opt_b0.set()) {
	      read_volume(b0,opt_b0.value()+"z_dz");
	      calc_gradientsROI(b0,b0x,b0y,b0z,myid,Nxx,numprocs);
	    }
	    else {
	      b0=phantom[0]*0;
	      b0x=b0;
	      b0y=b0;
	      b0z=b0;
	    }  
	    print_volume_info(b0,"b0");
	    ////////////////
	    //MAIN LOOP
	    ////////////////
	    cout<<"Main loop..."<<endl;
	    for (register int tt=0;tt<Nt;tt++){
	      for (register int zz=zstart;zz<zend;zz++){
		cout<<"Tissue type="<<tt<<"; zstart="<<zstart<<"; zz="<<zz<<"; zend="<<zend<<"; Voxelnumber="<<voxelcounter<<endl;
		for (register int yy=ystart;yy<yend;yy++){
		  int xxx=myid+xstart;
		  for (register int xx=xstart;xx<xend;xx++){
		      //slice speed up stuff
		    if (fabs((float)slcdir)==1) sszz=zz;
		    else if (fabs((float)slcdir)==2) sszz=yy;
		    else sszz=xx;
		    float slctmp=(sszz-sszstart_p)/(slcthk*nvox);
		    if(ceil(slctmp) == slctmp) {
		      sszz_slc=(int)(ceil(slctmp)+1);
		    } else {
		      sszz_slc=(int) (ceil(slctmp));
		    }
		    if (phantom(xx,yy,zz,tt)!=0){
		      voxelcounter=voxelcounter+1;
		      //cout<<"xx= "<<xx<<"; yy= "<<yy<<"; zz= "<<zz<<"; RFrec(xx,yy,zz)= "<<RFrec(xx,yy,zz)<<"; RFtrans(xx,yy,zz)= "<<RFtrans(xx,yy,zz)<<endl;
		      double den=phantom(xx,yy,zz,tt)*RFrec(xx,yy,zz)*cxyz;
		      if (opt_activation.set()) {
		        for (int n=0;n<=Nact-1;n++){
		          double a=activation(xx,yy,zz)*timecourse_2[n];
			  double b=tissue(tt+1,2);
			  activation4D_voxel[n]=a*b/(a+b);// conversion of beta into beta1 because of the integral (see possum no 7 page 121)
		        }
		      }
		      else if (opt_activation4D.set()){
		        for (int n=0;n<=Nact-1;n++){
		          double a=activation4D(xx,yy,zz,n);//zz-zstart_p when having only a pieace of the activation volume so we start from where the phantom starts 
		          double b=tissue(tt+1,2);
		          activation4D_voxel[n]=a*b/(a+b);
		        }
		      }
		      voxel2(posx(xxx+1),posy(yy+1),posz(zz+1),tissue.Row(tt+1),
		             pulse,nrf,nreadp,voxelcounter,xdim,ydim,zdim,
		             b0(xx,yy,zz),b0x(xx,yy,zz),b0y(xx,yy,zz),b0z(xx,yy,zz),
		             timecourse,activation4D_voxel,Nact,outputname,
		             table_slcprof,dslcp,dslcp_first,Nslc,den,RFtrans(xx,yy,zz),
		             opt_test,nospeedup,save_kcoord,sreal,simag);
		      }
		      xxx=xxx+numprocs;
		  }
		}
	      }
	    }
	  }
	  ////////////////////////////////////////////////////////////
	  //MOTION INVOLVING ROTATION Rx or Ry or both
	  ////////////////////////////////////////////////////////////
	  if (level==3){
	    cout<<"LEVEL3"<<endl;
	    ////////////////////////////////////////////
	    // B0 PERTURBATION 
	    ////////////////////////////////////////////
	    cout<<"Creating 9 B0file together with 27 gradient files..."<<endl;
	    volume<double> b0x_dx, b0x_dy, b0x_dz, b0y_dx, b0y_dy, b0y_dz, b0z_dx, 
	      b0z_dy, b0z_dz;//read in
	    volume<double> b0x_dx_gx, b0x_dx_gy, b0x_dx_gz, b0x_dy_gx, b0x_dy_gy, 
	      b0x_dy_gz, b0x_dz_gx, b0x_dz_gy, b0x_dz_gz;//calculate from, the calc gradients
	    volume<double> b0y_dx_gx, b0y_dx_gy, b0y_dx_gz, b0y_dy_gx, b0y_dy_gy, 
	      b0y_dy_gz, b0y_dz_gx, b0y_dz_gy, b0y_dz_gz;
	    volume<double> b0z_dx_gx, b0z_dx_gy, b0z_dx_gz, b0z_dy_gx, b0z_dy_gy, 
	      b0z_dy_gz, b0z_dz_gx, b0z_dz_gy, b0z_dz_gz;
	    if (opt_b0.set()) {
	      read_volume(b0x_dx,opt_b0.value()+"x_dx");
	      calc_gradientsROI(b0x_dx,b0x_dx_gx,b0x_dx_gy,b0x_dx_gz,myid,Nxx,numprocs);
	      read_volume(b0x_dy,opt_b0.value()+"x_dy");
	      calc_gradientsROI(b0x_dy,b0x_dy_gx,b0x_dy_gy,b0x_dy_gz,myid,Nxx,numprocs);
	      read_volume(b0x_dz,opt_b0.value()+"x_dz");
	      calc_gradientsROI(b0x_dz,b0x_dz_gx,b0x_dz_gy,b0x_dz_gz,myid,Nxx,numprocs);
	      read_volume(b0y_dx,opt_b0.value()+"y_dx");
	      calc_gradientsROI(b0y_dx,b0y_dx_gx,b0y_dx_gy,b0y_dx_gz,myid,Nxx,numprocs);
	      read_volume(b0y_dy,opt_b0.value()+"y_dy");
	      calc_gradientsROI(b0y_dy,b0y_dy_gx,b0y_dy_gy,b0y_dy_gz,myid,Nxx,numprocs);
	      read_volume(b0y_dz,opt_b0.value()+"y_dz");
	      calc_gradientsROI(b0y_dz,b0y_dz_gx,b0y_dz_gy,b0y_dz_gz,myid,Nxx,numprocs);
	      read_volume(b0z_dx,opt_b0.value()+"z_dx");
	      calc_gradientsROI(b0z_dx,b0z_dx_gx,b0z_dx_gy,b0z_dx_gz,myid,Nxx,numprocs);
	      read_volume(b0z_dy,opt_b0.value()+"z_dy");
	      calc_gradientsROI(b0z_dy,b0z_dy_gx,b0z_dy_gy,b0z_dy_gz,myid,Nxx,numprocs);
	      read_volume(b0z_dz,opt_b0.value()+"z_dz");
	      calc_gradientsROI(b0z_dz,b0z_dz_gx,b0z_dz_gy,b0z_dz_gz,myid,Nxx,numprocs);
	    }
	    else {
	      b0x_dx=phantom[0]*0;
	      b0x_dx_gx=b0x_dx; b0x_dx_gy=b0x_dx; b0x_dx_gz=b0x_dx;
	      b0x_dy=b0x_dx; b0x_dy_gx=b0x_dx; b0x_dy_gy=b0x_dx; b0x_dy_gz=b0x_dx;
	      b0x_dz=b0x_dx; b0x_dz_gx=b0x_dx; b0x_dz_gy=b0x_dx; b0x_dz_gz=b0x_dx;
	      b0y_dx=b0x_dx; b0y_dx_gx=b0x_dx; b0y_dx_gy=b0x_dx; b0y_dx_gz=b0x_dx;
	      b0y_dy=b0x_dx; b0y_dy_gx=b0x_dx; b0y_dy_gy=b0x_dx; b0y_dy_gz=b0x_dx;
	      b0y_dz=b0x_dx; b0y_dz_gx=b0x_dx; b0y_dz_gy=b0x_dx; b0y_dz_gz=b0x_dx;
	      b0z_dx=b0x_dx; b0z_dx_gx=b0x_dx; b0z_dx_gy=b0x_dx; b0z_dx_gz=b0x_dx;
	      b0z_dy=b0x_dx; b0z_dy_gx=b0x_dx; b0z_dy_gy=b0x_dx; b0z_dy_gz=b0x_dx;
	      b0z_dz=b0x_dx; b0z_dz_gx=b0x_dx; b0z_dz_gy=b0x_dx; b0z_dz_gz=b0x_dx;
	    }
	    print_volume_info(b0z_dz,"b0z_dz");
	    ///////////////////
	    //MAIN LOOP
	    ///////////////////
	    cout<<"Main loop..."<<endl;
	    int nonzero=0;
	    for (register int tt=0;tt<Nt;tt++){
	      for (register int zz=zstart;zz<zend;zz++){
		for (register int yy=ystart;yy<yend;yy++){
		  for (register int xx=xstart;xx<xend;xx++){
		    if (phantom(xx,yy,zz,tt)>1e-05) nonzero++;
		  }
		}
	      }
	    }
	    cout<<"The number of non-zero voxels is "<<nonzero<<endl;
	    RowVector numpointstmp;
	    numpointstmp=read_ascii_matrix(opt_mainmatrix.value()+".numpoints");
	    int numpoints=(int)numpointstmp(1);
	    int seg=(int)numpointstmp(2);
	    if (seg==0 || seg>numpoints) seg=numpoints;
	    PMatrix pulse(seg,19);
	    pulse=0.0;
	    int segA=1;int segB=segA+seg-1;
	    while (segA<numpoints){
	      voxelcounter=0;
	      cout<<"Reading the main matrix.Portion SegA="<<segA<<" SegB="<<segB<<endl;
	      if ((segB-segA+1)!=pulse.Nrows()) pulse.ReSize(segB-segA+1,19);
	      read_binary_matrix(pulse,opt_mainmatrix.value(),segA,segB,1,19);//add boundaries
	      cout<<"MatrixFile="<<opt_mainmatrix.value()<<endl;
	      nrf=0;
	      for (int tmp=1;tmp<=(segB-segA+1);tmp++){
		if (fabs(pulse(tmp,2))>1e-06) nrf++;
		//cerr<<" Nrf="<<nrf<<"Pulse("<<tmp<<",2)="<<pulse(tmp,2)<<endl;
	      }
	      for (register int tt=0;tt<Nt;tt++){
		for (register int zz=zstart;zz<zend;zz++){
		  cout<<"Tissue type="<<tt<<"; zstart="<<zstart<<"; zz="<<zz<<"; zend=";
		  cout<<zend<<"; Voxelnumber="<<voxelcounter<<endl;
		  for (register int yy=ystart;yy<yend;yy++){
		    int xxx=myid+xstart;
		    for (register int xx=0;xx<xend;xx++){
		      //slice speed up stuff
		      if (fabs((float)slcdir)==1) sszz=zz;
		      else if (fabs((float)slcdir)==2) sszz=yy;
		      else sszz=xx;
		      float slctmp=(sszz-sszstart_p)/(slcthk*nvox);
		      if(ceil(slctmp) == slctmp) {
		        sszz_slc=(int)(ceil(slctmp)+1);
		      } else {
		        sszz_slc=(int) (ceil(slctmp));
		      }
		      if (phantom(xx,yy,zz,tt)>1e-05){
		        voxelcounter=voxelcounter+1;
			double den=phantom(xx,yy,zz,tt)*RFrec(xx,yy,zz)*cxyz;
		        if (opt_activation.set()) {
		          for (int n=0;n<=Nact-1;n++){
		            double a=activation(xx,yy,zz)*timecourse_2[n];
			    double b=tissue(tt+1,2);
			    activation4D_voxel[n]=a*b/(a+b);// conversion of beta 
			    //into beta1 because of the integral (see possum no 7 page 121)
		          }
			}
			else if (opt_activation4D.set()){
		          for (int n=0;n<=Nact-1;n++){
		            double a=activation4D(xx,yy,zz,n);//zz-zstart_p 
			    //when having only a pieace of the activation volume so we start 
			    //from where the phantom starts 
		            double b=tissue(tt+1,2);
		            activation4D_voxel[n]=a*b/(a+b); 
			  }
			}
			voxel3(posx(xxx+1),posy(yy+1),posz(zz+1),tissue.Row(tt+1),
		               pulse ,segA, nrf,nreadp,nonzero,voxelcounter,xdim,ydim,zdim,
		               b0x_dx(xx,yy,zz),b0y_dx(xx,yy,zz),b0z_dx(xx,yy,zz),
		               b0x_dy(xx,yy,zz),b0y_dy(xx,yy,zz),b0z_dy(xx,yy,zz),
			       b0x_dz(xx,yy,zz),b0y_dz(xx,yy,zz),b0z_dz(xx,yy,zz),
		               b0x_dx_gx(xx,yy,zz),b0y_dx_gx(xx,yy,zz),b0z_dx_gx(xx,yy,zz),
			       b0x_dy_gx(xx,yy,zz),b0y_dy_gx(xx,yy,zz),b0z_dy_gx(xx,yy,zz),
			       b0x_dz_gx(xx,yy,zz),b0y_dz_gx(xx,yy,zz),b0z_dz_gx(xx,yy,zz),
			       b0x_dx_gy(xx,yy,zz),b0y_dx_gy(xx,yy,zz),b0z_dx_gy(xx,yy,zz),
			       b0x_dy_gy(xx,yy,zz),b0y_dy_gy(xx,yy,zz),b0z_dy_gy(xx,yy,zz),
			       b0x_dz_gy(xx,yy,zz),b0y_dz_gy(xx,yy,zz),b0z_dz_gy(xx,yy,zz),
		               b0x_dx_gz(xx,yy,zz),b0y_dx_gz(xx,yy,zz),b0z_dx_gz(xx,yy,zz),
			       b0x_dy_gz(xx,yy,zz),b0y_dy_gz(xx,yy,zz),b0z_dy_gz(xx,yy,zz),
			       b0x_dz_gz(xx,yy,zz),b0y_dz_gz(xx,yy,zz),b0z_dz_gz(xx,yy,zz),
		               timecourse,activation4D_voxel,Nact,outputname,table_slcprof,
		               dslcp,dslcp_first,Nslc,den,RFtrans(xx,yy,zz),opt_test,
		               nospeedup,save_kcoord,opt_rfavg.value(),
			       sreal,simag);
		      }
		    xxx=xxx+numprocs;
		    }
		  }
		}
	      }
	      segA=segB;
	      segB=segA+seg-1;
	      if (segB>numpoints) segB=numpoints;
	    }
	  }
	  ///////////////////////////////////////////////////////////
	  //NO MOTION BUT B0 PERTURBATION CHANGES WITH TIME
	  ///////////////////////////////////////////////////////////
	  if (level==4){
	    cout<<"LEVEL4"<<endl;
	    cout<<"Reading the pulse sequence..."<<endl;
	    PMatrix pulse;
	    read_binary_matrix(pulse, opt_pulse.value());  
	    ////////////////////////
	    // B0 PERTURBATION 
	    ////////////////////////
	    cout<<"Creating 1 B0extra file together with 3 gradient files all changing in time..."<<endl;
	    Matrix b0timecourse_tmp;
	    double* b0timecourse;
	    double* b0timecourse_2;
	    double* b0time;
	    double* b0xtime;
	    double* b0ytime; 
	    double* b0ztime;
	    b0timecourse_tmp=read_ascii_matrix(opt_b0timecourse4D.value());
	    int Nb0=b0timecourse_tmp.Nrows();
	    b0timecourse=new double[Nb0];
	    b0timecourse_2=new double[Nb0];
	    b0time=new double[Nb0];
	    b0xtime=new double[Nb0];
	    b0ytime=new double[Nb0];
	    b0ztime=new double[Nb0];
	    for (int n=0;n<=Nb0-1;n++){
	      b0timecourse[n]=b0timecourse_tmp(n+1,1);
	      b0timecourse_2[n]=b0timecourse_tmp(n+1,2);
	    }
	    volume<double> b0extra;
	    volume<double> b0xextra;
	    volume<double> b0yextra;
	    volume<double> b0zextra;
	    read_volume(b0extra,opt_b0extra.value());
	    calc_gradientsROI(b0extra,b0xextra,b0yextra,b0zextra,myid,Nxx,numprocs);
	    print_volume_info(b0extra,"b0extra");
	    cout<<"Creating 1 B0file together with 3 gradient files..."<<endl;
	    volume<double> b0(Nxx,Ny,Nz);
	    volume<double> b0x(Nxx,Ny,Nz);
	    volume<double> b0y(Nxx,Ny,Nz);
	    volume<double> b0z(Nxx,Ny,Nz);
	    if (opt_b0.set()) {
	      read_volume(b0,opt_b0.value()+"z_dz");
	      calc_gradientsROI(b0,b0x,b0y,b0z,myid,Nxx,numprocs);
	    }
	    else {
	      b0=phantom[0]*0;
	      b0x=b0;
	      b0y=b0;
	      b0z=b0;
	    }  
	    print_volume_info(b0,"b0");
	    cout<<""<<endl;
	    /////////////
	    //MAIN LOOP
	    /////////////
	    cout<<"Main loop..."<<endl;
	    for (register int tt=0;tt<Nt;tt++){
	      for (register int zz=zstart;zz<zend;zz++){
		cout<<"Tissue type="<<tt<<"; zstart="<<zstart<<"; zz="<<zz<<"; zend="<<zend<<"; Voxelnumber="<<voxelcounter<<endl;
		for (register int yy=ystart;yy<yend;yy++){
		  int xxx=myid+xstart;
		  for (register int xx=xstart;xx<xend;xx++){
		    //slice speed up stuff
		    if (fabs((float)slcdir)==1) sszz=zz;
		    else if (fabs((float)slcdir)==2) sszz=yy;
		    else sszz=xx;
		    float slctmp=(sszz-sszstart_p)/(slcthk*nvox);
		    if(ceil(slctmp) == slctmp) {
		      sszz_slc=(int)(ceil(slctmp)+1);
		      } else {
		      sszz_slc=(int) (ceil(slctmp));
		      }
		    if (phantom(xx,yy,zz,tt)!=0){
		        voxelcounter=voxelcounter+1;
		        double den=phantom(xx,yy,zz,tt)*RFrec(xx,yy,zz)*cxyz;
		        if (opt_activation.set()) {
		          for (int n=0;n<=Nact-1;n++){
		            double a=activation(xx,yy,zz)*timecourse_2[n];
			    double b=tissue(tt+1,2);
			    activation4D_voxel[n]=a*b/(a+b);// conversion of beta into beta1 because of the integral (see possum no 7 page 121)
			  }
			}
			else if (opt_activation4D.set()){
		          for (int n=0;n<=Nact-1;n++){
		            double a=activation4D(xx,yy,zz,n);//zz-zstart_p when having only a pieace of the activation volume so we start from where the phantom starts 
		            double b=tissue(tt+1,2);
		            activation4D_voxel[n]=a*b/(a+b); 
			  }              
			}
		        for (int n=0;n<=Nb0-1;n++){
		          b0time[n]=b0extra(xx,yy,zz)*b0timecourse_2[n];
		          b0xtime[n]=b0xextra(xx,yy,zz)*b0timecourse_2[n];
		          b0ytime[n]=b0yextra(xx,yy,zz)*b0timecourse_2[n];
		          b0ztime[n]=b0zextra(xx,yy,zz)*b0timecourse_2[n];
			  if (voxelcounter==1){
			    cout<<"b0time[n]="<<b0time[n]<<"b0xtime[n]="<<b0xtime[n]<<"b0ytime[n]="<<b0ytime[n]<<"b0ztime[n]="<<b0ztime[n]<<endl;
			  }
			} 
		    	voxel4(posx(xxx+1),posy(yy+1),posz(zz+1),tissue.Row(tt+1),
		               pulse,nreadp,voxelcounter,xdim,ydim,zdim,
		               b0time,b0xtime,b0ytime,b0ztime, b0timecourse, Nb0,
		               b0(xx,yy,zz),b0x(xx,yy,zz),b0y(xx,yy,zz),b0z(xx,yy,zz),
		               timecourse,activation4D_voxel,Nact,outputname,table_slcprof,
		               dslcp,dslcp_first,Nslc,den,RFtrans(xx,yy,zz),opt_test,
		               nospeedup,save_kcoord,sreal,simag);  
		    }
		    xxx=xxx+numprocs;
		  }
		}
	      }
	    }
	  }

  /////////////////
  //PARALEL STUFF
  #ifdef USE_MPI
    double* sreal_total;
    sreal_total= new double[nreadp];
    double* simag_total;
    simag_total= new double[nreadp];
    MPI_Reduce(sreal, sreal_total, nreadp, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);//MPI::COMM_WORLD.Reduce(sreal,sreal_total,nreadp,MPI::DOUBLE,MPI::SUM,0);
    MPI_Reduce(simag, simag_total, nreadp, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);//MPI::COMM_WORLD.Reduce(simag,simag_total,nreadp,MPI::DOUBLE,MPI::SUM,0);
    cout<<"Element nmb 1091 on  "<<myid<<" is <"<<sreal[1091]<<", "<<simag[1091]<<">.\n"; 
    if (myid==0){
      cout<<"Element nmb 1091 on "<<myid<<" is <"<<sreal_total[1091]<<", "<<simag_total[1091]<<">.\n"; 
      for (int i=1;i<=nreadp;i++) {
        signal(1,i)=sreal_total[i-1]*1e06; //1e06 we need to make signal have larger values so that it more resembles the scanner
        signal(2,i)=simag_total[i-1]*1e06; //1e06 we need to make signal have larger values so that it more resembles the scanner
      }
      write_binary_matrix(signal,opt_signal.value());
    }
    MPI_Finalize();//MPI::Finalize();
  #else
   /////////////////////////
   //OUTPUT
   /////////////////////////
    for (int i=1;i<=nreadp;i++) {
      signal(1,i)=sreal[i-1]*1e06; //1e06 we need to make signal have larger values so that it more resembles the scanner
      signal(2,i)=simag[i-1]*1e06; //1e06 we need to make signal have larger values so that it more resembles the scanner
    }
    write_binary_matrix(signal,opt_signal.value());
    cout<<"Possum finished generating the signal for "<<voxelcounter<<" voxels."<<endl;
  #endif

//tejas-22.11.12
//start time
theTime = time(NULL);
cout << "possum-end-time\t" << opt_procid.value() << "\t"  << theTime << endl;
//end

  return 0;
}


int main (int argc, char *argv[])
{

  Tracer tr("main");
  OptionParser options(title, examples);
  try {
    options.add(verbose);
    options.add(help);
    options.add(opt_kcoord);
    options.add(opt_object);
    options.add(opt_tissue);
    options.add(opt_motion);
    options.add(opt_pulse);
    options.add(opt_b0);
    options.add(opt_b0extra);
    options.add(opt_b0timecourse4D);
    options.add(opt_RFrec);
    options.add(opt_RFtrans);
    options.add(opt_activation);
    options.add(opt_timecourse);
    options.add(opt_activation4D);
    options.add(opt_timecourse4D);
    options.add(opt_level);
    options.add(opt_nproc);
    options.add(opt_procid);
    options.add(opt_signal);
    options.add(opt_slcprof);
    options.add(opt_mainmatrix);
    options.add(opt_nospeedup);
    options.add(opt_rfavg);
    
    nonoptarg = options.parse_command_line(argc, argv);

    // line below stops the program if there are less than 2 non-optional args
    //   or the help was requested or a compulsory option was not set
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
  compute_volume(argc, argv);
  return 0;
}
