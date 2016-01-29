/*    Copyright (C) 2012 University of Oxford  */

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

#include "streamlines.h"
#include "warpfns/fnirt_file_reader.h"
#include "warpfns/warpfns.h"


namespace TRACT{

  ColumnVector mean_sph_pol(ColumnVector& A, ColumnVector& B){
    // A and B contain th, ph f. 
    float th,ph;
    ColumnVector rA(3), rB(3);

    rA << (sin(A(1))*cos(A(2))) << (sin(A(1))*sin(A(2))) << (cos(A(1)));
    rB << (sin(B(1))*cos(B(2))) << (sin(B(1))*sin(B(2))) << (cos(B(1)));
    
    if(sum(SP(rA,rB)).AsScalar()>0)
      cart2sph((rA+rB)/2,th,ph);
    else
      cart2sph((rA-rB)/2,th,ph);
    
    A(1)=th; A(2)=ph;
    return A;
  }
  
  void read_masks(vector<string>& masks, const string& filename){
    ifstream fs(filename.c_str());
    string tmp;
    if(fs){
      fs>>tmp;
      do{
	masks.push_back(tmp);
	fs>>tmp;
      }while(!fs.eof());
    }
    else{
      cerr<<filename<<" does not exist"<<endl;
      exit(0);
    }
  }
  void imgradient(const volume<float>& im,volume4D<float>& grad){
    
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
  void imlookup(const volume<float>& im,volume<int>& vol2mat,Matrix& mat2vol){
    vol2mat.reinitialize(im.xsize(),im.ysize(),im.zsize());
    vol2mat = 0;

    int nrows=0;
    for(int Wz=im.minz();Wz<=im.maxz();Wz++)
      for(int Wy=im.miny();Wy<=im.maxy();Wy++)
	for(int Wx=im.minx();Wx<=im.maxx();Wx++)
	  if(im(Wx,Wy,Wz)!=0){
	    nrows++;
	  }
    mat2vol.ReSize(nrows,3);
    nrows=0;
    for(int Wz=im.minz();Wz<=im.maxz();Wz++)
      for(int Wy=im.miny();Wy<=im.maxy();Wy++)
	for(int Wx=im.minx();Wx<=im.maxx();Wx++)
	  if(im(Wx,Wy,Wz)!=0){
	    nrows++;
	    mat2vol(nrows,1) = Wx;
	    mat2vol(nrows,2) = Wy;
	    mat2vol(nrows,3) = Wz;
	    
	    vol2mat(Wx,Wy,Wz) = nrows;

	  }
  }
  
  void  imfill(volume<float>& im,const ColumnVector& vec,const Matrix& lookup){
    im = 0;
    for(int i=1;i<=lookup.Nrows();i++)
      im((int)lookup(i,1),(int)lookup(i,2),(int)lookup(i,3)) = vec(i);
  }

  
  Streamliner::Streamliner(const volume<float>& seeds):opts(probtrackxOptions::getInstance()),
						       logger(LogSingleton::getInstance()),
						       vols(opts.usef.value()),
						       m_seeds(seeds){
    
    read_volume(m_mask,opts.maskfile.value());
    m_part.initialise(0,0,0,0,0,0,opts.steplength.value(),m_mask.xdim(),m_mask.ydim(),m_mask.zdim(),false);
    if(opts.skipmask.value()!="") read_volume(m_skipmask,opts.skipmask.value());
    m_lcrat=5;
    if(opts.loopcheck.value()){
      m_loopcheck.reinitialize(int(ceil(m_mask.xsize()/m_lcrat)+1),int(ceil(m_mask.ysize()/m_lcrat)+1),int(ceil(m_mask.zsize()/m_lcrat)+1),3);
      m_loopcheck=0;
    }
    if(opts.rubbishfile.value()!=""){
      read_volume(m_rubbish,opts.rubbishfile.value());
    }
    if(opts.stopfile.value()!=""){
      read_volume(m_stop,opts.stopfile.value());
    }
    if(opts.prefdirfile.value()!=""){
      read_volume4D(m_prefdir,opts.prefdirfile.value());
    }
 
    vector<string> masknames;
    if(opts.waypoints.value()!=""){
      if(fsl_imageexists(opts.waypoints.value())){
	masknames.push_back( opts.waypoints.value() );
      }
      else{
	read_masks(masknames,opts.waypoints.value());
      }

      for( unsigned int m = 0; m < masknames.size(); m++ ){
	volume<float>* tmpptr =new volume<float>;
	if(opts.verbose.value()>0)
	  cout<<masknames[m]<<endl;
	read_volume(*tmpptr,masknames[m]);
	m_waymasks.push_back(tmpptr);
	m_passed_flags.push_back(false);
	m_own_waymasks.push_back(true);
      }
    }
    
    
    // Allow for either matrix transform (12dof affine) or nonlinear (warpfield)
    m_Seeds_to_DTI = IdentityMatrix(4);
    m_DTI_to_Seeds = IdentityMatrix(4);
    m_rotdir = IdentityMatrix(3);


    m_IsNonlinXfm = false;
    if(opts.seeds_to_dti.value()!=""){
      if(!fsl_imageexists(opts.seeds_to_dti.value())){
	m_Seeds_to_DTI = read_ascii_matrix(opts.seeds_to_dti.value());
	m_DTI_to_Seeds = m_Seeds_to_DTI.i();

	// set rotation matrix
	Matrix F(3,3),u(3,3),v(3,3);
	DiagonalMatrix d(3);

	if(opts.meshfile.value()!=""){
	  F << -m_Seeds_to_DTI(1,1) << m_Seeds_to_DTI(1,3) << -m_Seeds_to_DTI(1,2)
	    << -m_Seeds_to_DTI(2,1) << m_Seeds_to_DTI(2,3) << -m_Seeds_to_DTI(2,2)
	    << -m_Seeds_to_DTI(3,1) << m_Seeds_to_DTI(3,3) << -m_Seeds_to_DTI(3,2);
	}
	else{
	  F = m_Seeds_to_DTI.SubMatrix(1,3,1,3);
	}
	SVD(F*F.t(),d,u,v);
	m_rotdir.ReSize(3,3);
	m_rotdir = (u*sqrt(d)*v.t()).i()*F;
	
      }
      else{
	m_IsNonlinXfm = true;
	FnirtFileReader ffr(opts.seeds_to_dti.value());
	m_Seeds_to_DTI_warp = ffr.FieldAsNewimageVolume4D(true);
	if(opts.dti_to_seeds.value()==""){
	  cerr << "TRACT::Streamliner:: DTI -> Seeds transform needed" << endl;
	  exit(1);
	}
	FnirtFileReader iffr(opts.dti_to_seeds.value());
	m_DTI_to_Seeds_warp = iffr.FieldAsNewimageVolume4D(true);

	// now calculate the jacobian of this transformation (useful for rotating vectors)
	imgradient(m_Seeds_to_DTI_warp[0],m_jacx);
	imgradient(m_Seeds_to_DTI_warp[1],m_jacy);
	imgradient(m_Seeds_to_DTI_warp[2],m_jacz);

      }
    }
    
    
    vols.initialise(opts.basename.value());
    m_path.reserve(opts.nparticles.value());
    m_x_s_init=0;
    m_y_s_init=0;
    m_z_s_init=0;

    m_inmask3.reserve(opts.nsteps.value());

  }
  
  
  int Streamliner::streamline(const float& x_init,const float& y_init,const float& z_init, const ColumnVector& dim_seeds,const int& fibst,const ColumnVector& dir){ 
    
    //fibst tells tractvolsx which fibre to start with if there are more than one..
    //x_init etc. are in seed space...
    vols.reset(fibst);
    m_x_s_init=x_init; //seed x position in voxels
    m_y_s_init=y_init; // and y
    m_z_s_init=z_init; // and z
    ColumnVector xyz_seeds(3);
    xyz_seeds<<x_init<<y_init<<z_init;
    ColumnVector xyz_dti;
    ColumnVector th_ph_f;
    float xst,yst,zst,x,y,z,tmp2;

    // find xyz in dti space
    if(!m_IsNonlinXfm)
      xyz_dti = vox_to_vox(xyz_seeds,dim_seeds,vols.dimensions(),m_Seeds_to_DTI);    
    else{
      xyz_dti = NewimageCoord2NewimageCoord(m_DTI_to_Seeds_warp,false,m_seeds,m_mask,xyz_seeds);
    }

    xst=xyz_dti(1);yst=xyz_dti(2);zst=xyz_dti(3);
    m_path.clear();
    x=xst;y=yst;z=zst;
    m_part.change_xyz(x,y,z);
    
    if(opts.meshfile.value()!=""){ 
      m_part.set_dir(dir(1),dir(2),dir(3));//Set the start dir so that we track inwards from cortex 
    }

    bool rubbish_passed=false;
    bool stop_flag=false;
    
      //NB - this only goes in one direction!!
    
      
    float pathlength=0;
    for( int it = 1 ; it <= opts.nsteps.value()/2; it++){
      if( (m_mask( MISCMATHS::round(m_part.x()), MISCMATHS::round(m_part.y()), MISCMATHS::round(m_part.z())) > 0) ){

	///////////////////////////////////
	//loopchecking
	///////////////////////////////////
	if(opts.loopcheck.value()){
	  float oldrx=m_loopcheck((int)MISCMATHS::round(m_part.x()/m_lcrat),(int)MISCMATHS::round(m_part.y()/m_lcrat),(int)MISCMATHS::round(m_part.z()/m_lcrat),0);
	  float oldry=m_loopcheck((int)MISCMATHS::round(m_part.x()/m_lcrat),(int)MISCMATHS::round(m_part.y()/m_lcrat),(int)MISCMATHS::round(m_part.z()/m_lcrat),1);
	  float oldrz=m_loopcheck((int)MISCMATHS::round(m_part.x()/m_lcrat),(int)MISCMATHS::round(m_part.y()/m_lcrat),(int)MISCMATHS::round(m_part.z()/m_lcrat),2);
	  if(m_part.rx()*oldrx+m_part.ry()*oldry+m_part.rz()*oldrz<0)
	    {
	      break;
	    }
	    
	  m_loopcheck((int)MISCMATHS::round(m_part.x()/m_lcrat),(int)MISCMATHS::round(m_part.y()/m_lcrat),(int)MISCMATHS::round(m_part.z()/m_lcrat),0)=m_part.rx();
	  m_loopcheck((int)MISCMATHS::round(m_part.x()/m_lcrat),(int)MISCMATHS::round(m_part.y()/m_lcrat),(int)MISCMATHS::round(m_part.z()/m_lcrat),1)=m_part.ry();
	  m_loopcheck((int)MISCMATHS::round(m_part.x()/m_lcrat),(int)MISCMATHS::round(m_part.y()/m_lcrat),(int)MISCMATHS::round(m_part.z()/m_lcrat),2)=m_part.rz();  
	    
	}
	
	if(opts.verbose.value()>1)
	  logger<<m_part;
	
	x=m_part.x();y=m_part.y();z=m_part.z();
	xyz_dti <<x<<y<<z;


	// now find xyz in seeds space
	if(!m_IsNonlinXfm)
	  xyz_seeds = vox_to_vox(xyz_dti,vols.dimensions(),dim_seeds,m_DTI_to_Seeds);    
	else{
	  xyz_seeds = NewimageCoord2NewimageCoord(m_Seeds_to_DTI_warp,false,m_mask,m_seeds,xyz_dti);
	}

	
	int x_s =(int)MISCMATHS::round((float)xyz_seeds(1));
	int y_s =(int)MISCMATHS::round((float)xyz_seeds(2));
	int z_s =(int)MISCMATHS::round((float)xyz_seeds(3));
	
	float pref_x=0,pref_y=0,pref_z=0;
	if(opts.prefdirfile.value()!=""){
	  pref_x = m_prefdir(x_s,y_s,z_s,0);
	  pref_y = m_prefdir(x_s,y_s,z_s,1);
	  pref_z = m_prefdir(x_s,y_s,z_s,2);
	}
	//update every passed_flag
	for( unsigned int wm=0;wm<m_waymasks.size();wm++ ){
	  if( (*m_waymasks[wm])(x_s,y_s,z_s)!=0 ) {
	    m_passed_flags[wm]=true;
	  }
	}
	
	m_path.push_back(xyz_seeds);
	//m_path.push_back(xyz_dti);
	pathlength += opts.steplength.value();
	
	// //////////////////////////////
	// update coordinates for matrix3
	if(opts.matrix3out.value()){
	  if( m_mask3(x_s,y_s,z_s)!=0 ){
	    if( m_beenhere3(x_s,y_s,z_s)==0 ){
	      m_beenhere3(x_s,y_s,z_s)=1;
	      m_inmask3.push_back(xyz_seeds);
	    }
	  }
	}
	// //////////////////////////////
	

	if(opts.rubbishfile.value()!=""){
	  if(m_rubbish(x_s,y_s,z_s)!=0){
	    rubbish_passed=true;
	    break;
	  }
	}
	if(opts.stopfile.value()!=""){
	  if(m_stop(x_s,y_s,z_s)!=0){
	    stop_flag=true;
	  }
	  if(stop_flag)break;
	}	  
	  
	if(opts.skipmask.value() == ""){
	  th_ph_f=vols.sample(m_part.x(),m_part.y(),m_part.z(),
			      m_part.rx(),m_part.ry(),m_part.rz(),
			      pref_x,pref_y,pref_z);
	}
	else{
	  if(m_skipmask(x_s,y_s,z_s)==0)
	    th_ph_f=vols.sample(m_part.x(),m_part.y(),m_part.z(),m_part.rx(),m_part.ry(),m_part.rz(),pref_x,pref_y,pref_z);
	}
	  

	tmp2=rand(); tmp2/=RAND_MAX;
	if(th_ph_f(3)>tmp2){
	  if(!m_part.check_dir(th_ph_f(1),th_ph_f(2),opts.c_thr.value())){
	    
	    break;
	  }
	    
	  if((th_ph_f(1)!=0&&th_ph_f(2)!=0)){
	    if( (m_mask( MISCMATHS::round(m_part.x()), MISCMATHS::round(m_part.y()), MISCMATHS::round(m_part.z())) != 0) ){
	      if(!opts.modeuler.value())
		m_part.jump(th_ph_f(1),th_ph_f(2));
	      else
		{
		  ColumnVector test_th_ph_f;
		  
		  m_part.testjump(th_ph_f(1),th_ph_f(2));
		  test_th_ph_f=vols.sample(m_part.testx(),m_part.testy(),m_part.testz(),m_part.rx(),m_part.ry(),m_part.rz(),pref_x,pref_y,pref_z);
		  test_th_ph_f=mean_sph_pol(th_ph_f,test_th_ph_f);
		  m_part.jump(test_th_ph_f(1),test_th_ph_f(2));
		  
		}
	    }
	    
	    
	  }
	}
	
      }

      
    } // Close Step Number Loop
    if(opts.loopcheck.value()){
      m_loopcheck=0;
    }
    
    int rejflag=0;
    if(m_passed_flags.size()!=0){
      unsigned int numpassed=0;
      for(unsigned int i=0; i<m_passed_flags.size();i++){
	if(m_passed_flags[i])numpassed++;
      }
      if(numpassed==0)rejflag=1;
      else if(numpassed<m_passed_flags.size())rejflag=2;
      else rejflag=0;
    }   
    if(rubbish_passed){
      rejflag=1;
    }
    if(pathlength<opts.distthresh.value()){
      rejflag=1;
    }

    return rejflag;
  }

  void Streamliner::rotdir(const ColumnVector& dir,ColumnVector& rotdir,
			   const float& x,const float& y,const float& z){
    if(!m_IsNonlinXfm){
      rotdir = m_rotdir*dir;       
    }
    else{
      ColumnVector xyz_dti(3),xyz_seeds(3);
      xyz_seeds << x << y << z;
      xyz_dti = NewimageCoord2NewimageCoord(m_DTI_to_Seeds_warp,false,m_seeds,m_mask,xyz_seeds);

      Matrix F(3,3),Jw(3,3);
      Jw << m_jacx((int)xyz_dti(1),(int)xyz_dti(2),(int)xyz_dti(3),0) << m_jacx((int)xyz_dti(1),(int)xyz_dti(2),(int)xyz_dti(3),1) << m_jacx((int)xyz_dti(1),(int)xyz_dti(2),(int)xyz_dti(3),2)
	 << m_jacy((int)xyz_dti(1),(int)xyz_dti(2),(int)xyz_dti(3),0) << m_jacy((int)xyz_dti(1),(int)xyz_dti(2),(int)xyz_dti(3),1) << m_jacy((int)xyz_dti(1),(int)xyz_dti(2),(int)xyz_dti(3),2)
	 << m_jacz((int)xyz_dti(1),(int)xyz_dti(2),(int)xyz_dti(3),0) << m_jacz((int)xyz_dti(1),(int)xyz_dti(2),(int)xyz_dti(3),1) << m_jacz((int)xyz_dti(1),(int)xyz_dti(2),(int)xyz_dti(3),2);
	 	 
      F = (IdentityMatrix(3) + Jw).i();
      rotdir = F*dir;
    }
  }

  
  
  void Counter::initialise(){
    // set lookup table for seed mask

    imlookup(m_seeds,m_seeds_vol2mat,m_seeds_mat2vol);
   
    if(opts.simpleout.value()||opts.matrix1out.value()){
      initialise_path_dist();
    }
    if(opts.s2tout.value()){
      initialise_seedcounts();
    }
    if(opts.matrix1out.value()){
      initialise_matrix1();
    }
    if(opts.matrix2out.value()){
      initialise_matrix2();
    }
    if(opts.matrix3out.value()){
      initialise_matrix3();
    }
    if(opts.maskmatrixout.value()){
      initialise_maskmatrix();
    }
  }

  
  void Counter::initialise_seedcounts(){
    
    // store only nonzero values of the target masks
    volume<float> tmptarget,alltargets;
    volume<int> tmpint;
    
    ColumnVector scounter(m_numseeds);
    scounter=0;

    read_masks(m_targetmasknames,opts.targetfile.value());
    m_targflags.resize(m_targetmasknames.size(),0);
    cout<<"Number of masks "<<m_targetmasknames.size()<<endl;
    cout << "loading target masks - stage 1" << endl;

    alltargets.reinitialize(m_seeds.xsize(),m_seeds.ysize(),m_seeds.zsize());
    alltargets = 0;
    for(unsigned int m=0;m<m_targetmasknames.size();m++){
      cout << m+1 << endl;
      read_volume(tmptarget,m_targetmasknames[m]);
      alltargets += NEWIMAGE::abs(tmptarget);
      m_seedcounts.push_back(scounter);
    }
    imlookup(alltargets,m_targets_vol2mat,m_targets_mat2vol);

    cout << "loading target masks - stage 2" << endl;
    m_targetmasks.ReSize(m_targets_mat2vol.Nrows(),m_targetmasknames.size());
    m_targetmasks = 0;

    for(unsigned int m=0;m<m_targetmasknames.size();m++){
      cout << m+1 << endl;
      read_volume(tmptarget,m_targetmasknames[m]);

      for(int Wz=tmptarget.minz();Wz<=tmptarget.maxz();Wz++)
	for(int Wy=tmptarget.miny();Wy<=tmptarget.maxy();Wy++)
	  for(int Wx=tmptarget.minx();Wx<=tmptarget.maxx();Wx++){
	    if(tmptarget(Wx,Wy,Wz)!=0){
	      m_targetmasks( m_targets_vol2mat(Wx,Wy,Wz), m+1 ) = 1;
	    }
	  }
      
    }

    m_SeedCountMat.ReSize(m_numseeds,m_targetmasknames.size());
    m_SeedCountMat=0;
    m_SeedRow=1;
    

  }
  

  void Counter::initialise_matrix1(){
    m_Conrow=0;
    
    m_ConMat.reinitialize(m_numseeds,m_numseeds,1);
    m_CoordMat.reinitialize(m_numseeds,3,1);
    int myrow=0;
      
    for(int Wz=m_seeds.minz();Wz<=m_seeds.maxz();Wz++){
      for(int Wy=m_seeds.miny();Wy<=m_seeds.maxy();Wy++){
	for(int Wx=m_seeds.minx();Wx<=m_seeds.maxx();Wx++){
	  if(m_seeds(Wx,Wy,Wz)!=0){
	    m_CoordMat(myrow,0,0)=Wx;
	    m_CoordMat(myrow,1,0)=Wy;
	    m_CoordMat(myrow,2,0)=Wz;
	    myrow++;
	  }
	}
      }
    }
    
  }
  
  void Counter::initialise_matrix2(){
    
    m_Conrow2=0;
    read_volume(m_lrmask,opts.lrmask.value());
    m_beenhere2.reinitialize(m_lrmask.xsize(),m_lrmask.ysize(),m_lrmask.zsize());
    m_lrdim.ReSize(3);
    m_lrdim<<m_lrmask.xdim()<<m_lrmask.ydim()<<m_lrmask.zdim();
    int numnz=0;
    
    for(int Wz=m_lrmask.minz();Wz<=m_lrmask.maxz();Wz++)
      for(int Wy=m_lrmask.miny();Wy<=m_lrmask.maxy();Wy++)
	for(int Wx=m_lrmask.minx();Wx<=m_lrmask.maxx();Wx++)
	  if(m_lrmask.value(Wx,Wy,Wz)!=0)
	    numnz++;
    
    
    if(numnz> pow(2,(float)sizeof(short)*8-1)-1){
      cerr<<"Output matrix too big for AVW - stopping."<<endl;
      cerr<<" Remember - you can store your tracts in "<<endl;
      cerr<<" low res even if you want your seeds in high res"<<endl;
      cerr<<" Just subsample the structural space mask"<<endl;
      cerr<<" Although, it must stay in line with the seeds"<<endl;
      exit(-1);
    }

    m_ConMat2.reinitialize(m_numseeds,numnz,1);
    m_ConMat2 = 0;

    m_CoordMat2.reinitialize(m_numseeds,3,1);
    m_CoordMat_tract2.reinitialize(numnz,3,1);
    
    Matrix tempy(numnz,1);
    for(int i=1;i<=numnz;i++){tempy(i,1)=i-1;}
    m_lookup2.addvolume(m_lrmask);
    m_lookup2.setmatrix(tempy.t(),m_lrmask);
      
    int mytrow=0;
    for(int Wz=m_lrmask.minz();Wz<=m_lrmask.maxz();Wz++)
      for(int Wy=m_lrmask.miny();Wy<=m_lrmask.maxy();Wy++)
	for(int Wx=m_lrmask.minx();Wx<=m_lrmask.maxx();Wx++)
	  if(m_lrmask(Wx,Wy,Wz)!=0){
	    m_CoordMat_tract2(mytrow,0,0)=Wx;
	    m_CoordMat_tract2(mytrow,1,0)=Wy;
	    m_CoordMat_tract2(mytrow,2,0)=Wz;
	    mytrow++;
	  }
    
    int myrow=0;
    for(int Wz=m_seeds.minz();Wz<=m_seeds.maxz();Wz++)
      for(int Wy=m_seeds.miny();Wy<=m_seeds.maxy();Wy++)
	for(int Wx=m_seeds.minx();Wx<=m_seeds.maxx();Wx++)
	  if(m_seeds(Wx,Wy,Wz)!=0){
	    m_CoordMat2(myrow,0,0)=Wx;
	    m_CoordMat2(myrow,1,0)=Wy;
	    m_CoordMat2(myrow,2,0)=Wz;
	    myrow++;
	  }
      
      
  }
  
  // the following will use the voxels of a given mask to initialise 
  // an NxN matrix. This matrix will store the number of samples from 
  // each seed voxel that have made it to each pair of voxels in the mask
  void Counter::initialise_matrix3(){
    volume<int>& m_mask3           = m_stline.get_mask3();
    volume<int>& m_beenhere3       = m_stline.get_beenhere3();

    read_volume(m_mask3,opts.maskmatrix3.value());
    m_beenhere3.reinitialize(m_mask3.xsize(),m_mask3.ysize(),m_mask3.zsize());
    m_beenhere3=0;

    int nmask3=0;
    for(int z=0;z<m_mask3.zsize();z++){
      for(int y=0;y<m_mask3.ysize();y++){
	for(int x=0;x<m_mask3.xsize();x++){
	  if(m_mask3(x,y,z)==0)continue;
	  nmask3++;
	}
      }
    }
    m_CoordMat3.ReSize(nmask3,3);
    //m_ConMat3.reinitialize(nmask3,nmask3,1);
    //m_ConMat3=0;
    m_ConMat3 = new SpMat<int> (nmask3,nmask3); // is this how you do it?

    m_Lookup3.reinitialize(m_mask3.xsize(),m_mask3.ysize(),m_mask3.zsize());

    nmask3=1;
    for(int Wz=m_mask3.minz();Wz<=m_mask3.maxz();Wz++)
      for(int Wy=m_mask3.miny();Wy<=m_mask3.maxy();Wy++)
	for(int Wx=m_mask3.minx();Wx<=m_mask3.maxx();Wx++){
	  if(m_mask3(Wx,Wy,Wz)==0)continue;
	  m_CoordMat3(nmask3,1) = Wx;
	  m_CoordMat3(nmask3,2) = Wy;
	  m_CoordMat3(nmask3,3) = Wz;
	  m_Lookup3(Wx,Wy,Wz)  = nmask3;
	  nmask3++;
	}
    //write_ascii_matrix(m_CoordMat3,logger.appendDir("coords_for_fdt_matrix3.txt"));
  }
  void Counter::count_streamline(){
    if(opts.simpleout.value()||opts.matrix1out.value()){
      update_pathdist();
    }
    if(opts.s2tout.value()){
      update_seedcounts();
    }
    if(opts.matrix2out.value()){
      update_matrix2_row();
    }
    if(opts.maskmatrixout.value()){
      update_maskmatrix();
    }
  }
  
  void Counter::count_seed(){
    if(opts.matrix1out.value()){
      update_matrix1();
    }
    if(opts.matrix2out.value()){
      next_matrix2_row();
    }
    if(opts.seedcountastext.value()){
      m_SeedRow++;
    }
  }
  
    
  void Counter::clear_streamline(){
    if(opts.simpleout.value()||opts.matrix1out.value()){
      reset_beenhere();
    }
    if(opts.s2tout.value()){
      reset_targetflags();
    }
    if(opts.matrix2out.value()){
      reset_beenhere2();
    }
    if(opts.maskmatrixout.value()){
      //Do whatever it is you have to do!!
    }
    if(opts.matrix3out.value()){
      reset_beenhere3();
    }
  }
  
  void Counter::update_pathdist(){
    //const vector<ColumnVector>& path=m_stline.get_path_ref();

    if(!opts.pathdist.value()){
      for(unsigned int i=0;i<m_path.size();i++){
	int x_s=int(MISCMATHS::round(float(m_path[i](1))));
	int y_s=int(MISCMATHS::round(float(m_path[i](2))));
	int z_s=int(MISCMATHS::round(float(m_path[i](3))));
	if(m_beenhere(x_s,y_s,z_s)==0){
	  m_prob(x_s,y_s,z_s)+=1; 
	  m_beenhere(x_s,y_s,z_s)=1;
	}
      }
    }
    else{
      int d=1;
      for(unsigned int i=0;i<m_path.size();i++){
	int x_s=int(MISCMATHS::round(float(m_path[i](1))));
	int y_s=int(MISCMATHS::round(float(m_path[i](2))));
	int z_s=int(MISCMATHS::round(float(m_path[i](3))));
	if(m_beenhere(x_s,y_s,z_s)==0){
	  m_prob(x_s,y_s,z_s)+=d;d++;
	  m_beenhere(x_s,y_s,z_s)=1;
	}
      }
    }

    
  }

  void Counter::reset_beenhere(){
    for(unsigned int i=0;i<m_path.size();i++){
      int x_s=int(MISCMATHS::round(float(m_path[i](1)))),y_s=int(MISCMATHS::round(float(m_path[i](2)))),z_s=int(MISCMATHS::round(float(m_path[i](3))));
      m_beenhere(x_s,y_s,z_s)=0;
    }
  }
  
  
  void Counter::update_seedcounts(){

    //const vector<ColumnVector>& path=m_stline.get_path_ref();
    int xseedvox=int(MISCMATHS::round(m_stline.get_x_seed()));
    int yseedvox=int(MISCMATHS::round(m_stline.get_y_seed()));
    int zseedvox=int(MISCMATHS::round(m_stline.get_z_seed()));

    float pathlength;
    if(!opts.pathdist.value()){
      pathlength=0;
      for(unsigned int i=0;i<m_path.size();i++){
	int x_s=int(MISCMATHS::round(float(m_path[i](1)))),y_s=int(MISCMATHS::round(float(m_path[i](2)))),z_s=int(MISCMATHS::round(float(m_path[i](3))));
	for(unsigned int m=0;m<m_targetmasknames.size();m++){
	  if(m_targets_vol2mat(x_s,y_s,z_s)!=0 && m_targflags[m]==0)
	    if(m_targetmasks(m_targets_vol2mat(x_s,y_s,z_s),m+1)!=0){
	      if(pathlength>=opts.distthresh.value()){
		m_seedcounts[m](m_seeds_vol2mat(xseedvox,yseedvox,zseedvox)) += 1;
		if(opts.seedcountastext.value())
		  m_SeedCountMat(m_SeedRow,m+1) += 1;
	      }
	      m_targflags[m]=1;	    
	    }
	}
	pathlength += opts.steplength.value();
      }
    }
    else{
      int x_s,y_s,z_s;
      pathlength=0;
      for(unsigned int i=0;i<m_path.size();i++){
	x_s=int(MISCMATHS::round(float(m_path[i](1))));y_s=int(MISCMATHS::round(float(m_path[i](2))));z_s=int(MISCMATHS::round(float(m_path[i](3))));
	for(unsigned int m=0;m<m_targetmasknames.size();m++){
	  if(m_targets_vol2mat(x_s,y_s,z_s)!=0 && m_targflags[m]==0)
	    if(m_targetmasks(m_targets_vol2mat(x_s,y_s,z_s),m+1)!=0){
	      if(pathlength>=opts.distthresh.value()){
		m_seedcounts[m](m_seeds_vol2mat(xseedvox,yseedvox,zseedvox)) += pathlength;
	      if(opts.seedcountastext.value())
		m_SeedCountMat(m_SeedRow,m+1) += pathlength;
	      }
	      m_targflags[m]=1;	
	    }
	}
	pathlength += opts.steplength.value();
      }
    }
    

  }
  

  
  void Counter::update_matrix1(){
    //after each particle, update_pathdist(), only run this after each voxel
    int Concol=0;
    for(int Wz=m_prob.minz();Wz<=m_prob.maxz();Wz++){
      for(int Wy=m_prob.miny();Wy<=m_prob.maxy();Wy++){
	for(int Wx=m_prob.minx();Wx<=m_prob.maxx();Wx++){
	  if(m_seeds(Wx,Wy,Wz)!=0){
	    if(m_prob(Wx,Wy,Wz)!=0){
	      m_ConMat(m_Conrow,Concol,0)=m_prob(Wx,Wy,Wz);
	    }
	    Concol++;
	  }
	  m_prob(Wx,Wy,Wz)=0;
	  
	}
      }
    }
    
    m_Conrow++;
  }
  
  void Counter::update_matrix2_row(){
    //run this one every streamline - not every voxel..
    //const vector<ColumnVector>& path=m_stline.get_path_ref();

    if(!opts.pathdist.value()){
      for(unsigned int i=0;i<m_path.size();i++){
	ColumnVector xyz_seeds=m_path[i];
	//do something here
	ColumnVector xyz_lr=vox_to_vox(xyz_seeds,m_seedsdim,m_lrdim,m_I);
	
	int x_lr=int(MISCMATHS::round(float(xyz_lr(1)))),y_lr=int(MISCMATHS::round(float(xyz_lr(2)))),z_lr=int(MISCMATHS::round(float(xyz_lr(3))));
	int Concol2=m_lookup2(x_lr,y_lr,z_lr,0);
	if(Concol2!=0){
	  if(m_beenhere2(x_lr,y_lr,z_lr)==0){
	    m_ConMat2(m_Conrow2,Concol2,0)+=1;
	    m_beenhere2(x_lr,y_lr,z_lr)=1;
	  }
	}
	
      }
    }
    else{
      int d=1;
      for(unsigned int i=0;i<m_path.size();i++){
	ColumnVector xyz_seeds=m_path[i];
	ColumnVector xyz_lr=vox_to_vox(xyz_seeds,m_seedsdim,m_lrdim,m_I);
	int x_lr=int(MISCMATHS::round(float(xyz_lr(1)))),y_lr=int(MISCMATHS::round(float(xyz_lr(2)))),z_lr=int(MISCMATHS::round(float(xyz_lr(3))));
	int Concol2=m_lookup2(x_lr,y_lr,z_lr,0);
	if(Concol2!=0){
	  if(m_beenhere2(x_lr,y_lr,z_lr)==0){
	    m_ConMat2(m_Conrow2,Concol2,0)+=d;d++;
	    m_beenhere2(x_lr,y_lr,z_lr)=1;
	  }
	}
      }
    }
    
  }

  void Counter::update_matrix3(){
    vector<ColumnVector>& inmask3   = m_stline.get_inmask3();    
    
    if(inmask3.size()<2)return;
    
    float length;
    for(unsigned int i=0;i<inmask3.size();i++){
      length=0;
      for(unsigned int j=i+1;j<inmask3.size();j++){
	length += 1;//opts.steplength.value();
	int row1 = m_Lookup3((int)MISCMATHS::round(float(inmask3[i](1))),(int)MISCMATHS::round(float(inmask3[i](2))),(int)MISCMATHS::round(float(inmask3[i](3))));
	int row2 = m_Lookup3((int)MISCMATHS::round(float(inmask3[j](1))),(int)MISCMATHS::round(float(inmask3[j](2))),(int)MISCMATHS::round(float(inmask3[j](3))));
	//m_ConMat3(row1,row2,0) += 1;
	//m_ConMat3(row2,row1,0) += 1;

	if(opts.distthresh3.value()>0){
	  if(length<opts.distthresh3.value())
	    continue;
	}
	m_ConMat3->AddTo(row1,row2,1);

      }
    }
  }  
  void Counter::reset_beenhere3(){
    m_stline.clear_beenhere3();
    m_stline.clear_inmask3();
  }  
  
  void Counter::reset_beenhere2(){
    for(unsigned int i=0;i<m_path.size();i++){
      ColumnVector xyz_seeds=m_path[i];
      
      ColumnVector xyz_lr=vox_to_vox(xyz_seeds,m_seedsdim,m_lrdim,m_I);
      
      int x_lr=int(MISCMATHS::round(float(xyz_lr(1)))),y_lr=int(MISCMATHS::round(float(xyz_lr(2)))),z_lr=int(MISCMATHS::round(float(xyz_lr(3))));
      m_beenhere2(x_lr,y_lr,z_lr)=0;
    }    
  }

  void Counter::save_total(const int& keeptotal){
    
    // save total number of particles that made it through the streamlining
    ColumnVector keeptotvec(1);
    keeptotvec(1)=keeptotal;
    write_ascii_matrix(keeptotvec,logger.appendDir("waytotal"));

  }
  void Counter::save_total(const vector<int>& keeptotal){
    
    // save total number of particles that made it through the streamlining
    ColumnVector keeptotvec(keeptotal.size());
    for (int i=1;i<=(int)keeptotal.size();i++)
      keeptotvec(i)=keeptotal[i-1];
    write_ascii_matrix(keeptotvec,logger.appendDir("waytotal"));

  }

  void Counter::save(){
    cout << "now saving various outputs" << endl;
    if(opts.simpleout.value() && opts.mode.value()!="simple"){
      save_pathdist();
    }
    if(opts.s2tout.value()){
      save_seedcounts();
    }
    if(opts.matrix1out.value()){
      save_matrix1();
    }
    if(opts.matrix2out.value()){
      save_matrix2();
    }
    if(opts.matrix3out.value()){
      save_matrix3();
    }
    if(opts.maskmatrixout.value()){
      save_maskmatrix();
    }
    
  }
  
  void Counter::save_pathdist(){  
    m_prob.setDisplayMaximumMinimum(m_prob.max(),m_prob.min());
    save_volume(m_prob,logger.appendDir(opts.outfile.value()));
  }
  
  void Counter::save_pathdist(string add){  //for simple mode
    string thisout=opts.outfile.value();
    make_basename(thisout);
    thisout+=add;
    m_prob.setDisplayMaximumMinimum(m_prob.max(),m_prob.min());
    save_volume(m_prob,logger.appendDir(thisout));
  }

  void Counter::save_seedcounts(){
    volume<float> seedcounts;
    seedcounts.reinitialize(m_seeds.xsize(),m_seeds.ysize(),m_seeds.zsize());
    copybasicproperties(m_seeds,seedcounts);

    for(unsigned int m=0;m<m_targetmasknames.size();m++){
      string tmpname=m_targetmasknames[m];
      
      int pos=tmpname.find("/",0);
      int lastpos=pos;
      
      while(pos>=0){
	lastpos=pos;
	pos=tmpname.find("/",pos);
	// replace / with _
	tmpname[pos]='_';
      }
      
      //only take things after the last pos
      tmpname=tmpname.substr(lastpos+1,tmpname.length()-lastpos-1);
      
      imfill(seedcounts,m_seedcounts[m],m_seeds_mat2vol);

      seedcounts.setDisplayMaximumMinimum(opts.nparticles.value(),0);
      save_volume(seedcounts,logger.appendDir("seeds_to_"+tmpname));
    }

    if(opts.seedcountastext.value()){
      write_ascii_matrix(m_SeedCountMat,logger.appendDir("matrix_seeds_to_all_targets"));
    }

  }
    
  // the following is a helper function for save_matrix*
  //  to convert between nifti coords (external) and newimage coord (internal)
  void applycoordchange(volume<int>& coordvol, const Matrix& old2new_mat)
  {
    for (int n=0; n<=coordvol.maxx(); n++) {
      ColumnVector v(4);
      v << coordvol(n,0,0) << coordvol(n,1,0) << coordvol(n,2,0) << 1.0;
      v = old2new_mat * v;
      coordvol(n,0,0) = MISCMATHS::round(v(1));
      coordvol(n,1,0) = MISCMATHS::round(v(2));
      coordvol(n,2,0) = MISCMATHS::round(v(3));
    }
  }
  void applycoordchange(Matrix& coordvol, const Matrix& old2new_mat)
  {
    for (int n=1; n<=coordvol.Nrows(); n++) {
      ColumnVector v(4);
      v << coordvol(n,1) << coordvol(n,2) << coordvol(n,3) << 1.0;
      v = old2new_mat * v;
      coordvol(n,1) = MISCMATHS::round(v(1));
      coordvol(n,2) = MISCMATHS::round(v(2));
      coordvol(n,3) = MISCMATHS::round(v(3));
    }
  }

  void Counter::save_matrix1(){
    save_volume(m_ConMat,logger.appendDir("fdt_matrix1"));
    applycoordchange(m_CoordMat, m_seeds.niftivox2newimagevox_mat().i());
    save_volume(m_CoordMat,logger.appendDir("coords_for_fdt_matrix1"));
    applycoordchange(m_CoordMat, m_seeds.niftivox2newimagevox_mat());
  }

  void Counter::save_matrix2(){
    if(!opts.splitmatrix2.value()){
      save_volume(m_ConMat2,logger.appendDir("fdt_matrix2"));
      applycoordchange(m_CoordMat2, m_seeds.niftivox2newimagevox_mat().i());
      save_volume(m_CoordMat2,logger.appendDir("coords_for_fdt_matrix2"));
      applycoordchange(m_CoordMat2, m_seeds.niftivox2newimagevox_mat());
      applycoordchange(m_CoordMat_tract2, m_lrmask.niftivox2newimagevox_mat().i());
      save_volume(m_CoordMat_tract2,logger.appendDir("tract_space_coords_for_fdt_matrix2"));
      applycoordchange(m_CoordMat_tract2, m_lrmask.niftivox2newimagevox_mat());
      save_volume4D(m_lookup2,logger.appendDir("lookup_tractspace_fdt_matrix2"));
    }
    else{
      cout << "saving matrix2 into splitted files" << endl;

      int nsplits = 10;
      while( float(m_ConMat2.xsize()/nsplits) >= 32767 ){
	  nsplits++;
      }

      int nrows = std::floor(float(m_ConMat2.xsize()/nsplits))+1;
      volume<int> tmpmat;

      applycoordchange(m_CoordMat2, m_seeds.niftivox2newimagevox_mat().i());

      for(int i=1;i<=nsplits;i++){
	int first_row = (i-1)*nrows+1;
	int last_row  = i*nrows > m_ConMat2.xsize() ? m_ConMat2.xsize() : i*nrows;
	if(first_row > m_ConMat2.xsize()) break;

	// set limits
	m_ConMat2.setROIlimits(first_row-1,m_ConMat2.miny(),m_ConMat2.minz(),last_row-1,m_ConMat2.maxy(),m_ConMat2.maxz());
	m_ConMat2.activateROI();
	tmpmat = m_ConMat2.ROI();
	save_volume(tmpmat,logger.appendDir("fdt_matrix2_"+num2str(i)));

	m_CoordMat2.setROIlimits(first_row-1,m_CoordMat2.miny(),m_CoordMat2.minz(),last_row-1,m_CoordMat2.maxy(),m_CoordMat2.maxz());
	m_CoordMat2.activateROI();
	tmpmat = m_CoordMat2.ROI();
	save_volume(tmpmat,logger.appendDir("coords_for_fdt_matrix2_"+num2str(i)));


      }

      applycoordchange(m_CoordMat_tract2, m_lrmask.niftivox2newimagevox_mat());
      save_volume4D(m_lookup2,logger.appendDir("lookup_tractspace_fdt_matrix2"));

      applycoordchange(m_CoordMat2, m_seeds.niftivox2newimagevox_mat());
      applycoordchange(m_CoordMat_tract2, m_lrmask.niftivox2newimagevox_mat().i());
      save_volume(m_CoordMat_tract2,logger.appendDir("tract_space_coords_for_fdt_matrix2"));

    }
  }
void Counter::save_matrix3(){
  //save_volume(m_ConMat3,logger.appendDir("fdt_matrix3"));
  m_ConMat3->Print(logger.appendDir("fdt_matrix3.dot"));
  applycoordchange(m_CoordMat3, m_seeds.niftivox2newimagevox_mat().i());
  write_ascii_matrix(m_CoordMat3,logger.appendDir("coords_for_fdt_matrix3.txt"));
}
  int Seedmanager::run(const float& x,const float& y,const float& z,bool onewayonly, int fibst){
    ColumnVector dir(3);
    dir=0;
    return run(x,y,z,onewayonly,fibst,dir);
  }
  // this function now returns the total number of pathways that survived a streamlining (SJ)
  int Seedmanager::run(const float& x,const float& y,const float& z,bool onewayonly, int fibst,const ColumnVector& dir){
    //onewayonly for mesh things..
    cout <<x<<" "<<y<<" "<<z<<endl;
    if(opts.fibst.set()){
      fibst=opts.fibst.value()-1;

    }
    else{
      if(fibst == -1){
	fibst=0;//m_seeds(int(MISCMATHS::round(x)),int(MISCMATHS::round(y)),int(MISCMATHS::round(z)))-1;//fibre to start with is taken from seed volume..
      }
      //TB moved randfib option inside tractvols.h 28/10/2009
      // This means that we have access to fsamples when figuring out fibst
      // so we can choose to seed in proportion to f in that voxel. 
    }
    
    // now re-orient dir using xfm transform
    ColumnVector rotdir(3);
    m_stline.rotdir(dir,rotdir,x,y,z);
    //rotdir=dir;

    int nlines=0;
    for(int p=0;p<opts.nparticles.value();p++){
      if(opts.randfib.value()>2){ 
	//This bit of code just does random sampling from all fibre populations - even those whose f value is less than f-thresh. 
	//3 other possibilities - randfib==0 -> use fibst (default first fibre but can be set)
	// randfib==1 - random sampling of fibres bigger than fthresh
	// randfib==2 random sampling of fibres bigger than fthresh in proporthion to their f-values. 
	float tmp=rand()/float(RAND_MAX) * float(m_stline.nfibres()-1);
	fibst = (int)MISCMATHS::round(tmp);
      }
      
      // random sampling within a seed voxel
      float newx=x,newy=y,newz=z;
      if(opts.sampvox.value()){
	newx+=(float)rand()/float(RAND_MAX)-0.5;
	newy+=(float)rand()/float(RAND_MAX)-0.5;
	newz+=(float)rand()/float(RAND_MAX)-0.5;
      }

      if(opts.verbose.value()>1)
	logger.setLogFile("particle"+num2str(p));
   
      m_stline.reset(); //This now includes a vols.reset() in order to get fibst right. 
      bool forwardflag=false,backwardflag=false;

      int  rejflag1=1,rejflag2=1; // 0:accept, 1:reject, 2:wait
      m_counter.clear_path();

      // track in one direction
      if(!onewayonly || opts.matrix3out.value()){//always go both ways in matrix3 mode
	rejflag1 = m_stline.streamline(newx,newy,newz,m_seeddims,fibst,rotdir);
	if(rejflag1==0 || rejflag1==2){ 
	  forwardflag=true;
	  m_counter.append_path();
	}
	m_stline.reverse();
      }

      // track in the other direction
      rejflag2=m_stline.streamline(newx,newy,newz,m_seeddims,fibst,rotdir);
      if(rejflag2==0){	
	backwardflag=true;
      }
      if(rejflag2>0){
	backwardflag=false;
	if(rejflag1>0)
	  forwardflag=false;
      }

      if(!forwardflag)
	m_counter.clear_path();
      if(backwardflag)
	m_counter.append_path();

      if(forwardflag || backwardflag){
	nlines++; 
	m_counter.count_streamline();

	if(opts.matrix3out.value()){
	  m_counter.update_matrix3();
	}
      }

      m_counter.clear_streamline(); 
    }

    m_counter.count_seed();
    

    return nlines;
    
  }


}
  
  

