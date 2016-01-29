/*  vector_flirt.cc

    Saad Jbabdi, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford */

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

#include "vecreg.h"
#include "utils/options.h"
#include "warpfns/fnirt_file_reader.h"
#include "warpfns/warpfns.h"

using namespace Utilities;

string title="vecreg \nVector Affine/NonLinear Tranformation with Orientation Preservation";
string examples="vecreg -i <input4D> -o <output4D> -r <refvol> -t <transform>";

Option<bool> verbose(string("-v,--verbose"),false,
		       string("switch on diagnostic messages"),
		       false,no_argument);
Option<bool> help(string("-h,--help"),false,
		       string("display this message"),
		       false,no_argument);
Option<string> infilename(string("-i,--input"),string(""),
		       string("filename for input vector or tensor field"),
		       true,requires_argument);
Option<string> outfilename(string("-o,--output"),string(""),
		       string("filename for output registered vector or tensor field"),
		       true,requires_argument);
Option<string> reffname(string("-r,--ref"),string(""),
		       string("filename for reference (target) volume"),
		       true,requires_argument);
Option<string> matrix(string("-t,--affine"),string(""),
		       string("filename for affine transformation matrix"),
		       false,requires_argument);
Option<string> warp(string("-w,--warpfield"),string(""),
		       string("filename for 4D warp field for nonlinear registration"),
		       false,requires_argument);
Option<string> matrix2(string("--rotmat"),string(""),
		       string("filename for secondary affine matrix \n\t\t\tif set, this will be used for the rotation of the vector/tensor field "),
		       false,requires_argument);
Option<string> warp2(string("--rotwarp"),string(""),
		       string("filename for secondary warp field \n\t\t\tif set, this will be used for the rotation of the vector/tensor field "),
		       false,requires_argument);
Option<string> interpmethod(string("--interp"),"",
		       string("interpolation method : nearestneighbour, trilinear (default), sinc or spline"),
		       false,requires_argument);
Option<string> maskfile(string("-m,--mask"),string(""),
		       string("brain mask in input space"),
		       false,requires_argument);
Option<string> omaskfile(string("--refmask"),string(""),
			 string("brain mask in output space (useful for speed up of nonlinear reg)"),
			 false,requires_argument);
////////////////////////////////////////////////////////

ReturnMatrix rodrigues(const float& angle,ColumnVector& w){
  Matrix W(3,3),R(3,3);

  w/=sqrt(w.SumSquare()); // normalise w
  W <<  0.0  << -w(3) << -w(2)
    <<  w(3) << 0.0   << -w(1)
    << -w(2) << w(1)  <<  0.0;
  
  R << 1.0 << 0.0 << 0.0
    << 0.0 << 1.0 << 0.0
    << 0.0 << 0.0 << 1.0;
  R += (sin(angle)*W + (1-cos(angle))*(W*W));

  R.Release();
  return R;
}
ReturnMatrix rodrigues(const float& s,const float& c,ColumnVector& w){
  Matrix W(3,3),R(3,3);

  w /= sqrt(w.SumSquare()); // normalise w
  W <<  0.0  << -w(3) <<  w(2)
    <<  w(3) <<  0.0  << -w(1)
    << -w(2) <<  w(1) <<  0.0;
  
  R << 1.0 << 0.0 << 0.0
    << 0.0 << 1.0 << 0.0
    << 0.0 << 0.0 << 1.0;
  R += (s*W + (1-c)*(W*W));


  R.Release();
  return R;
}
ReturnMatrix rodrigues(const ColumnVector& n1,const ColumnVector& n2){
  ColumnVector w(3);
  
  w=cross(n1,n2);

  if(w.MaximumAbsoluteValue()>0){
    float ca=dot(n1,n2);
    float sa=sqrt(cross(n1,n2).SumSquare());
    
    return rodrigues(sa,ca,w);
  }
  else{
    Matrix R(3,3);
    R << 1.0 << 0.0 << 0.0
      << 0.0 << 1.0 << 0.0
      << 0.0 << 0.0 << 1.0;
    R.Release();
    return R;
  }
}
ReturnMatrix ppd(const Matrix& F,const ColumnVector& e1, const ColumnVector& e2){
  ColumnVector n1(3),n2(3),Pn2(3);
  Matrix R(3,3),R1(3,3),R2(3,3);

  n1=F*e1;
  if(n1.MaximumAbsoluteValue()>0)
    n1/=sqrt(n1.SumSquare());
  n2=F*e2;
  if(n2.MaximumAbsoluteValue()>0)
    n2/=sqrt(n2.SumSquare());

  R1=rodrigues(e1,n1);

  Pn2=cross(n1,n2);
  Pn2=n2-dot(n1,n2)*n1;Pn2=Pn2/sqrt(Pn2.SumSquare());
  R2=rodrigues(R1*e2/sqrt((R1*e2).SumSquare()),Pn2);
  R=R2*R1;

  R.Release();
  return R;
}
ReturnMatrix ppd(const Matrix& F,ColumnVector& e1){
  ColumnVector n1(3);
  Matrix R(3,3);

  e1/=sqrt(e1.SumSquare());

  n1=F*e1;
  n1=n1/sqrt(n1.SumSquare());

  R=rodrigues(e1,n1);

  R.Release();
  return R;
}

void sjgradient(const volume<float>& im,volume4D<float>& grad){
  
  grad.reinitialize(im.xsize(),im.ysize(),im.zsize(),3);
  copybasicproperties(im,grad[0]);

  int fx,fy,fz,bx,by,bz;
  float dx,dy,dz; 
  for (int z=0; z<grad.zsize(); z++){
    fz = z ==(grad.zsize()-1) ? 0 :  1;
    bz = z == 0              ? 0 : -1;
    dz = (fz==0 || bz==0)    ? 1.0 :  2.0;
    for (int y=0; y<grad.ysize(); y++){
      fy = y ==(grad.ysize()-1) ? 0 :  1;
      by = y == 0              ? 0 : -1;
      dy = (fy==0 || by==0)    ? 1.0 :  2.0;
      for (int x=0; x<grad.xsize(); x++){
	fx = x ==(grad.xsize()-1) ? 0 :  1;
	bx = x == 0              ? 0 : -1;
	dx = (fx==0 || bx==0)    ? 1.0 :  2.0;
	grad[0](x,y,z) = (im(x+fx,y,z) - im(x+bx,y,z))/dx;
	grad[1](x,y,z) = (im(x,y+fy,z) - im(x,y+by,z))/dy;
	grad[2](x,y,z) = (im(x,y,z+fz) - im(x,y,z+bz))/dz;
      }
    }
  }

}

void vecreg_aff(const volume4D<float>& tens,
		volume4D<float>& oV1,
		const volume<float>& refvol,
		const Matrix& M,
		const volume<float>& mask){

  Matrix iM(4,4);
  iM=M.i();

  /////////////////////////////////////////////////////////////////////////
  // Where we define a potential affine transfo for the rotation
  Matrix M2;
  if(matrix2.value()!="")
    M2 = read_ascii_matrix(matrix2.value());
  // Where we define a potential warp field transfo for the rotation
  volume4D<float> warpvol2;
  volume4D<float> jx,jy,jz;
  Matrix Jw(3,3),I(3,3);I<<1<<0<<0<<0<<1<<0<<0<<0<<1;
  if(warp2.value()!=""){
    FnirtFileReader ffr2(warp2.value());
    warpvol2=ffr2.FieldAsNewimageVolume4D(true);
    sjgradient(warpvol2[0],jx);
    sjgradient(warpvol2[1],jy);
    sjgradient(warpvol2[2],jz);
  }
  /////////////////////////////////////////////////////////////////////////
  volume<float> omask;
  if(omaskfile.value()!="")
    read_volume(omask,omaskfile.value());


  // extract rotation matrix from M
  Matrix F(3,3),R(3,3),u(3,3),v(3,3);
  DiagonalMatrix d(3);
  if(matrix2.value()=="")
    F=M.SubMatrix(1,3,1,3);
  else
    F=M2.SubMatrix(1,3,1,3);
  SVD(F*F.t(),d,u,v);
  R=(u*sqrt(d)*v.t()).i()*F;

  ColumnVector seeddim(3),targetdim(3);
  seeddim << tens.xdim() << tens.ydim() << tens.zdim();
  targetdim  << refvol.xdim() << refvol.ydim() << refvol.zdim();
  SymmetricMatrix Tens(3);
  Matrix FullTens(3,3);
  ColumnVector X_seed(3),X_target(3); 
  ColumnVector V_seed(3),V_target(3);

  for(int z=0;z<oV1.zsize();z++)
    for(int y=0;y<oV1.ysize();y++)
      for(int x=0;x<oV1.xsize();x++){
	if(omaskfile.value()!="")
	  if(omask(x,y,z)==0)
	    continue;

	// compute seed coordinates
	X_target << x << y << z;
	X_seed=vox_to_vox(X_target,targetdim,seeddim,iM);
	
	if(mask((int)MISCMATHS::round(float(X_seed(1))),(int)MISCMATHS::round(float(X_seed(2))),(int)MISCMATHS::round(float(X_seed(3))))==0)
	  continue;
	

	 // compute interpolated tensor
	if(oV1.tsize()!=9){
	  Tens.Row(1) << tens[0].interpolate(X_seed(1),X_seed(2),X_seed(3));
	  Tens.Row(2) << tens[1].interpolate(X_seed(1),X_seed(2),X_seed(3))
		      << tens[3].interpolate(X_seed(1),X_seed(2),X_seed(3));
	  Tens.Row(3) << tens[2].interpolate(X_seed(1),X_seed(2),X_seed(3))
		      << tens[4].interpolate(X_seed(1),X_seed(2),X_seed(3))
		      << tens[5].interpolate(X_seed(1),X_seed(2),X_seed(3));
	}
	else{
	  FullTens.Row(1) << tens[0].interpolate(X_seed(1),X_seed(2),X_seed(3))
			  << tens[3].interpolate(X_seed(1),X_seed(2),X_seed(3))
			  << tens[6].interpolate(X_seed(1),X_seed(2),X_seed(3));
	  FullTens.Row(2) << tens[1].interpolate(X_seed(1),X_seed(2),X_seed(3))
			  << tens[4].interpolate(X_seed(1),X_seed(2),X_seed(3))
			  << tens[7].interpolate(X_seed(1),X_seed(2),X_seed(3));
	  FullTens.Row(3) << tens[2].interpolate(X_seed(1),X_seed(2),X_seed(3))
			  << tens[5].interpolate(X_seed(1),X_seed(2),X_seed(3))
			  << tens[8].interpolate(X_seed(1),X_seed(2),X_seed(3));
	}

	 if(warp2.value()!=""){
	   // Local Jacobian of the backward warpfield
	   Jw <<   jx(x,y,z,0) <<  jx(x,y,z,1) << jx(x,y,z,2)
	      <<   jy(x,y,z,0) <<  jy(x,y,z,1) << jy(x,y,z,2)
	      <<   jz(x,y,z,0) <<  jz(x,y,z,1) << jz(x,y,z,2); 
	   // compute local forward affine transformation	
	   F = (I + Jw).i();
	 }
	 if(oV1.tsize()==3){ // case where input is a vector
	   // compute first eigenvector
	   EigenValues(Tens,d,v);
	   V_seed = v.Column(3);

	   // rotate vector
	   V_target=F*V_seed;
	   if(V_target.MaximumAbsoluteValue()>0)
	     V_target/=sqrt(V_target.SumSquare());
	   
	   oV1(x,y,z,0)=V_target(1);
	   oV1(x,y,z,1)=V_target(2);
	   oV1(x,y,z,2)=V_target(3);
	 }
	
	 // create Symmetric tensor
	 if(oV1.tsize()==6){
	   if(warp2.value()!=""){
	     EigenValues(Tens,d,v);
	     R=ppd(F,v.Column(3),v.Column(2));
	     Tens << R*Tens*R.t();
	   }
	   else{
	     Tens << R*Tens*R.t();
	   }
	   oV1(x,y,z,0)=Tens(1,1);
	   oV1(x,y,z,1)=Tens(2,1);
	   oV1(x,y,z,2)=Tens(3,1);
	   oV1(x,y,z,3)=Tens(2,2);
	   oV1(x,y,z,4)=Tens(3,2);
	   oV1(x,y,z,5)=Tens(3,3);
	 }

	 // create Non-Symmetric tensor
	 if(oV1.tsize()==9){
	   if(warp2.value()!=""){
	     SVD(FullTens,d,u);
	     R=ppd(F,u.Column(3),u.Column(2));
	     FullTens << R*FullTens*R.t();
	   }
	   else{
	     FullTens << R*FullTens*R.t();
	   }
	   oV1(x,y,z,0)=FullTens(1,1);
	   oV1(x,y,z,1)=FullTens(2,1);
	   oV1(x,y,z,2)=FullTens(3,1);
	   oV1(x,y,z,3)=FullTens(1,2);
	   oV1(x,y,z,4)=FullTens(2,2);
	   oV1(x,y,z,5)=FullTens(3,2);
	   oV1(x,y,z,6)=FullTens(1,3);
	   oV1(x,y,z,7)=FullTens(2,3);
	   oV1(x,y,z,8)=FullTens(3,3);
	 }

      }
}



void vecreg_nonlin(const volume4D<float>& tens,volume4D<float>& oV1,
		   const volume<float>& refvol,volume4D<float>& warpvol,
		   const volume<float>& mask){

  ColumnVector X_seed(3),X_target(3);
  
  
  // read warp field created by Jesper
  FnirtFileReader ffr(warp.value());
  warpvol=ffr.FieldAsNewimageVolume4D(true);
 
  Matrix F(3,3),u(3,3),v(3,3);
  DiagonalMatrix d(3);
    
  /////////////////////////////////////////////////////////////////////////
  // Where we define a potential affine transfo for the rotation
  Matrix M2;
  if(matrix2.value()!=""){
    M2 = read_ascii_matrix(matrix2.value());
    // extract rotation matrix from M
    F=M2.SubMatrix(1,3,1,3);
    SVD(F*F.t(),d,u,v);
    F=(u*sqrt(d)*v.t()).i()*F;
  }
  // Where we define a potential warp field transfo for the rotation
  volume4D<float> warpvol2;
  volume4D<float> jx,jy,jz;
  if(warp2.value()!=""){
    FnirtFileReader ffr2(warp2.value());
    warpvol2=ffr2.FieldAsNewimageVolume4D(true);
    sjgradient(warpvol2[0],jx);
    sjgradient(warpvol2[1],jy);
    sjgradient(warpvol2[2],jz);
  }
  else{
    sjgradient(warpvol[0],jx);
    sjgradient(warpvol[1],jy);
    sjgradient(warpvol[2],jz);
  }
  /////////////////////////////////////////////////////////////////////////
  volume<float> omask;
  if(omaskfile.value()!="")
    read_volume(omask,omaskfile.value());

  ColumnVector V_seed(3),V_target(3);
  ColumnVector V1_seed(3),V2_seed(3);
  ColumnVector V1_target(3),V2_target(3),V3_target(3);
  Matrix R(3,3),I(3,3);I<<1<<0<<0<<0<<1<<0<<0<<0<<1;
  Matrix Jw(3,3);
  SymmetricMatrix Tens(3);
  Matrix FullTens(3,3);
  for(int z=0;z<oV1.zsize();z++)
    for(int y=0;y<oV1.ysize();y++)
       for(int x=0;x<oV1.xsize();x++){
	 if(omaskfile.value()!="")
	   if(omask(x,y,z)==0)
	     continue;

	 X_target << x << y << z;
	 X_seed = NewimageCoord2NewimageCoord(warpvol,false,oV1[0],mask,X_target);


	 if(mask((int)MISCMATHS::round(float(X_seed(1))),(int)MISCMATHS::round(float(X_seed(2))),(int)MISCMATHS::round(float(X_seed(3))))==0)
	  continue;
	
	 // compute interpolated tensor
	if(oV1.tsize()!=9){
	  Tens.Row(1) << tens[0].interpolate(X_seed(1),X_seed(2),X_seed(3));
	  Tens.Row(2) << tens[1].interpolate(X_seed(1),X_seed(2),X_seed(3))
		      << tens[3].interpolate(X_seed(1),X_seed(2),X_seed(3));
	  Tens.Row(3) << tens[2].interpolate(X_seed(1),X_seed(2),X_seed(3))
		      << tens[4].interpolate(X_seed(1),X_seed(2),X_seed(3))
		      << tens[5].interpolate(X_seed(1),X_seed(2),X_seed(3));
	}
	else{
	  FullTens.Row(1) << tens[0].interpolate(X_seed(1),X_seed(2),X_seed(3))
			  << tens[3].interpolate(X_seed(1),X_seed(2),X_seed(3))
			  << tens[6].interpolate(X_seed(1),X_seed(2),X_seed(3));
	  FullTens.Row(2) << tens[1].interpolate(X_seed(1),X_seed(2),X_seed(3))
			  << tens[4].interpolate(X_seed(1),X_seed(2),X_seed(3))
			  << tens[7].interpolate(X_seed(1),X_seed(2),X_seed(3));
	  FullTens.Row(3) << tens[2].interpolate(X_seed(1),X_seed(2),X_seed(3))
			  << tens[5].interpolate(X_seed(1),X_seed(2),X_seed(3))
			  << tens[8].interpolate(X_seed(1),X_seed(2),X_seed(3));
	}

	 // F will not change if matrix2 is set
	 if(matrix2.value()==""){
	   // Local Jacobian of the backward warpfield
	   Jw <<   jx(x,y,z,0) <<  jx(x,y,z,1) << jx(x,y,z,2)
	      <<   jy(x,y,z,0) <<  jy(x,y,z,1) << jy(x,y,z,2)
	      <<   jz(x,y,z,0) <<  jz(x,y,z,1) << jz(x,y,z,2); 
	   // compute local forward affine transformation	
	   F = (I + Jw).i();
	 }

	 if(oV1.tsize()==3){// case where input is a vector
	   // compute first eigenvector
	   EigenValues(Tens,d,v);
	   V_seed = v.Column(3);
	 
	   V_target=F*V_seed;

	   if(V_target.MaximumAbsoluteValue()>0)
	     V_target/=sqrt(V_target.SumSquare());

	   oV1(x,y,z,0)=V_target(1);
	   oV1(x,y,z,1)=V_target(2);
	   oV1(x,y,z,2)=V_target(3);
	 }

	 // create Symmetric tensor
	 if(oV1.tsize()==6){
	   if(matrix2.value()==""){
	     EigenValues(Tens,d,v);
	     R=ppd(F,v.Column(3),v.Column(2));
	     Tens << R*Tens*R.t();
	   }
	   else{
	     Tens << F*Tens*F.t();
	   }
	   oV1(x,y,z,0)=Tens(1,1);
	   oV1(x,y,z,1)=Tens(2,1);
	   oV1(x,y,z,2)=Tens(3,1);
	   oV1(x,y,z,3)=Tens(2,2);
	   oV1(x,y,z,4)=Tens(3,2);
	   oV1(x,y,z,5)=Tens(3,3);
	 }
	
	 // create Non-Symmetric tensor
	 if(oV1.tsize()==9){
	   if(matrix2.value()==""){
	     SVD(FullTens,d,u);
	     R=ppd(F,u.Column(3),u.Column(2));
	     FullTens << R*FullTens*R.t();
	   }
	   else{
	     FullTens << F*FullTens*F.t();
	   }
	   oV1(x,y,z,0)=FullTens(1,1);
	   oV1(x,y,z,1)=FullTens(2,1);
	   oV1(x,y,z,2)=FullTens(3,1);
	   oV1(x,y,z,3)=FullTens(1,2);
	   oV1(x,y,z,4)=FullTens(2,2);
	   oV1(x,y,z,5)=FullTens(3,2);
	   oV1(x,y,z,6)=FullTens(1,3);
	   oV1(x,y,z,7)=FullTens(2,3);
	   oV1(x,y,z,8)=FullTens(3,3);
	 }

       }
}



int do_vecreg(){
  volume4D<float> ivol,warpvol;
  volume<float> refvol,mask;
  Matrix Aff(4,4);

  if((matrix.set())){
    Aff = read_ascii_matrix(matrix.value());
  }
  //if((warp.set())){
  //if(verbose.value()) cerr << "Loading warpfield" << endl;
  //read_volume4D(warpvol,warp.value());
  //}
  if(verbose.value()) cerr << "Loading volumes" << endl;
  read_volume4D(ivol,infilename.value());
  read_volume(refvol,reffname.value());

  volume4D<float> ovol;
  ovol.reinitialize(refvol.xsize(),refvol.ysize(),refvol.zsize(),ivol.tsize());
  copybasicproperties(refvol,ovol);
  ovol=0;

  // set interpolation method
  if(interpmethod.value()=="nearestneighbour")
    ivol.setinterpolationmethod(nearestneighbour);
  else if(interpmethod.value()=="sinc")
    ivol.setinterpolationmethod(sinc);
  else if(interpmethod.value()=="spline")
    ivol.setinterpolationmethod(spline);
  else
    ivol.setinterpolationmethod(trilinear);

  if(maskfile.value()!="")
    read_volume(mask,maskfile.value());
  else{
    mask.reinitialize(ivol[0].xsize(),ivol[0].ysize(),ivol[0].zsize());
    copybasicproperties(ivol,mask);
    for(int z=0;z<mask.zsize();z++)
      for(int y=0;y<mask.ysize();y++)
	for(int x=0;x<mask.xsize();x++){
	  if(abs(ivol(x,y,z,0))==0 && abs(ivol(x,y,z,1))==0 && abs(ivol(x,y,z,2))==0)
	    mask(x,y,z) = 0;
	  else
	    mask(x,y,z) = 1;
	}
  }

  ///////////////////////
  // tensor for interpolation
  volume4D<float> tens(ivol.xsize(),ivol.ysize(),ivol.zsize(),6);
  copybasicproperties(ivol,tens);
  if(ivol.tsize()==3){   //vector
    cout<<"Registering vector..."<<endl;
    for(int z=0;z<ivol.zsize();z++) 
      for(int y=0;y<ivol.ysize();y++)  
	for(int x=0;x<ivol.xsize();x++){
	  tens(x,y,z,0)=ivol(x,y,z,0)*ivol(x,y,z,0);
	  tens(x,y,z,1)=ivol(x,y,z,1)*ivol(x,y,z,0);
	  tens(x,y,z,2)=ivol(x,y,z,2)*ivol(x,y,z,0);
	  tens(x,y,z,3)=ivol(x,y,z,1)*ivol(x,y,z,1);
	  tens(x,y,z,4)=ivol(x,y,z,2)*ivol(x,y,z,1);
	  tens(x,y,z,5)=ivol(x,y,z,2)*ivol(x,y,z,2);
	}
  }
  else if (ivol.tsize()==6){ //symmetric tensor
    cout<<"Registering symmetric tensor..."<<endl;
    tens=ivol;
  }
  else{                     //Nonsymmetric tensor: Input elements are expected to be stored in Column-Order: (D11, D21, D31, D12, D22, D32, D13, D23, D33)
    cout<<"Registering non-symmetric tensor..."<<endl;
    tens.reinitialize(ivol.xsize(),ivol.ysize(),ivol.zsize(),9);
    tens=ivol;
  }

  //time_t _time=time(NULL);
  if(matrix.set()){
    if(verbose.value()) cerr << "Affine registration" << endl;
    vecreg_aff(tens,ovol,refvol,Aff,mask);
  }
  else{
    if(verbose.value()) cerr << "Nonlinear registration" << endl;
    vecreg_nonlin(tens,ovol,refvol,warpvol,mask);
  }
  //cout<<"elapsed time:"<<time(NULL)-_time<<" sec"<<endl;

  ovol.setDisplayMaximumMinimum(ivol.max(),ivol.min());
  save_volume4D(ovol,outfilename.value());

  return 0;

}


int main(int argc,char *argv[]){

  Tracer tr("main");
  OptionParser options(title,examples);

  try{
    options.add(verbose);
    options.add(help);
    options.add(infilename);
    options.add(outfilename);
    options.add(reffname);
    options.add(matrix);
    options.add(warp);
    options.add(matrix2);
    options.add(warp2);
    options.add(interpmethod);
    options.add(maskfile);
    options.add(omaskfile);

    options.parse_command_line(argc,argv);

    
    if ( (help.value()) || (!options.check_compulsory_arguments(true)) ){
      options.usage();
      exit(EXIT_FAILURE);
    }
    if( (matrix.set()) && (warp.set()) ){
      cerr << endl
	   << "Cannot specify both --affine AND --warpfield"
	   << endl << endl;
      exit(EXIT_FAILURE);
    }
    if( (matrix.unset()) && (warp.unset()) ){
      cerr << endl
	   << "Please Specify either --affine OR --warpfield"
	   << endl << endl;
      exit(EXIT_FAILURE);
    }
    if( (warp2.set()) && (matrix2.set()) ){
      cerr << endl
	   << "Cannot Specify both --rotaff AND --rotwarp"
	   << endl << endl;
      exit(EXIT_FAILURE);
    }
    
  }
  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } 
  catch(std::exception &e) {
    cerr << e.what() << endl;
  } 
  
  return do_vecreg();
  
  
}
