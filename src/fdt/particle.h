/*  Copyright (C) 2004 University of Oxford  */

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

#ifndef __PARTICLE_H_
#define __PARTICLE_H_



//////////////////////////////////////////////////////////////////
//      class Particle                                          //
//            tract particle..                                  //
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//                                                              //
//           NB - Everything in this Class is in voxels!!       // 
//                                                              //
//////////////////////////////////////////////////////////////////
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

namespace PARTICLE{

  class Particle
    { 
      float m_x;
      float m_y;
      float m_z;
      float m_rx;
      float m_ry;
      float m_rz;
      float m_rx_init;
      float m_ry_init;
      float m_rz_init;
      float m_testx;
      float m_testy;
      float m_testz;
      float m_steplength;
      float m_xdim;
      float m_ydim;
      float m_zdim;
      bool m_has_jumped;
      bool m_simdiff;
      int m_jumpsign;
    public:
      //constructors::
      Particle(const float& xin,const float& yin,
	       const float& zin,const float& rxin,
	       const float& ryin,const float &rzin,
	       const float& steplengthin,
	       const float& xdimin,
	       const float& ydimin,
	       const float& zdimin,
	       const bool& hasjumpedin=false,
	       const bool& simdiffin=false) : 
	m_x(xin), m_y(yin), m_z(zin), m_rx(rxin),m_ry(ryin),m_rz(rzin),m_rx_init(rxin), 
	m_ry_init(ryin),m_rz_init(rzin),m_steplength(steplengthin),
	m_xdim(xdimin),m_ydim(ydimin),m_zdim(zdimin),
	m_has_jumped(hasjumpedin),m_simdiff(false){}
      Particle(){}
      ~Particle(){}
      
      //initialise
      void initialise(const float& xin=0,const float& yin=0,
		      const float& zin=0,const float& rxin=0,
		      const float& ryin=0,const float &rzin=0,
		 const float& steplengthin=0.5,
		 const float& xdimin=2,
		 const float& ydimin=2,
		 const float& zdimin=2,
		 const bool& hasjumpedin=false,
		 const bool& simdiffin=false){
       
	m_x=xin;
	m_y=yin;
	m_z=zin;
	m_rx=rxin; 
	m_ry=ryin;
	m_rz=rzin;
        m_rx_init=rxin;
	m_ry_init=ryin;
	m_rz_init=rzin;
	m_steplength=steplengthin;
	m_xdim=xdimin;
	m_ydim=ydimin;
	m_zdim=zdimin;
	m_has_jumped=hasjumpedin;
	m_simdiff=simdiffin;

      }
      
      
      //return values
      const float& x() const { return m_x; }
      float x() { return m_x; }
      
      const float& y() const { return m_y; }
      float y() { return m_y; }
  
      const float& z() const { return m_z; }
      float z() { return m_z; }
  
      const float& rx() const { return m_rx; }
      float rx() { return m_rx; }
  
      const float& ry() const { return m_ry; }
      float ry() { return m_ry; }
  
      const float& rz() const { return m_rz; }
      float rz() { return m_rz; }

      const float& testx() const { return m_testx; }
      float testx() { return m_testx; }

      const float& testy() const { return m_testy; }
      float testy() { return m_testy; }

      const float& testz() const { return m_testz; }
      float testz() { return m_testz; }
      
      const float& steplength() const { return m_steplength; }
      float steplength() { return m_steplength; }
      
      //change values
      void change_x (float new_x) { m_x=new_x; }
      void change_y (float new_y) { m_y=new_y; }
      void change_z  (float new_z) { m_z=new_z; }
      void change_xyz (float new_x,float new_y,float new_z){
	 m_x=new_x;
	 m_y=new_y;
	 m_z=new_z;
      } 
      void change_steplength (float new_sl) { m_steplength = new_sl; } 
      void reset(){
	m_x=0;m_y=0;m_z=0;m_rx=0;m_ry=0;m_rz=0;m_has_jumped=false;
      }
      //functions
      void jump(const float& theta,const float& phi){
	float rx_new=cos(phi)*sin(theta);
	float ry_new=sin(phi)*sin(theta);
	float rz_new=cos(theta);
	int sign; bool init=false;
	if(!m_simdiff){
	  if(m_has_jumped)
	    {sign=(rx_new*m_rx + ry_new*m_ry + rz_new*m_rz)>0 ? 1:-1;}
	  else{
	    float tmp=rand(); tmp/=RAND_MAX;
	    sign=tmp > 0.5 ? 1:-1;
	    m_jumpsign=sign;
	    m_has_jumped=true;
	    init=true;
	  }
	}
	else{
	  float tmp=rand(); tmp/=RAND_MAX;
	  sign=tmp > 0.5 ? 1:-1;
	}
	m_x += sign*m_steplength/m_xdim*rx_new;
	m_y += sign*m_steplength/m_ydim*ry_new;
	m_z += sign*m_steplength/m_zdim*rz_new;
	m_rx=sign*rx_new; m_ry=sign*ry_new;m_rz=sign*rz_new;
	
	if(init){
	  m_rx_init=m_rx;
	  m_ry_init=m_ry;
	  m_rz_init=m_rz;
	}
	
	
      }
     

      void testjump(const float& theta,const float& phi){
	float rx_new=cos(phi)*sin(theta);
	float ry_new=sin(phi)*sin(theta);
	float rz_new=cos(theta);
	int sign;bool init=false;
	if(!m_simdiff){
	  if(m_has_jumped)
	    {sign=(rx_new*m_rx + ry_new*m_ry + rz_new*m_rz)>0 ? 1:-1;}
	  else{
	    float tmp=rand(); tmp/=RAND_MAX;
	    sign=tmp > 0.5 ? 1:-1;
	    m_jumpsign=sign;
	    // m_has_jumped=true; // causes probtrackx to go in one direction!
	    // init=true;
	  }
	}
	else{
	  float tmp=rand(); tmp/=RAND_MAX;
	  sign=tmp > 0.5 ? 1:-1;
	}
	m_testx = m_x+sign*m_steplength/m_xdim*rx_new;
	m_testy = m_y+sign*m_steplength/m_ydim*ry_new;
	m_testz = m_z+sign*m_steplength/m_zdim*rz_new;
	
	if(init){
	  m_rx_init=m_rx;
	  m_ry_init=m_ry;
	  m_rz_init=m_rz;
	}
	
      }
      
      
      void restart_reverse(){
	if(m_has_jumped){
	  m_rx=-m_rx_init;
	  m_ry=-m_ry_init;
	  m_rz=-m_rz_init;
	}
	
      }

      void set_dir(const float& rx,const float& ry,const float& rz){
	m_rx=rx;m_ry=ry;m_rz=rz;m_has_jumped=true;
      }
      
      
      bool check_dir(const float& theta,const float& phi, const float& thr){
	if(m_has_jumped){
	  float rx_new=cos(phi)*sin(theta);
	  float ry_new=sin(phi)*sin(theta);
	  float rz_new=cos(theta);
	  return fabs(rx_new*m_rx + ry_new*m_ry + rz_new*m_rz)>thr;
	}
	else return true;
      }

      // function added by Saad to choose a direction during deterministic streamlining
      // the choosed direction is the one that is closest to current direction
      unsigned int choose_dir(const vector<float>& th,const vector<float>& ph){
	float ps,tmpps=0;
	unsigned int r=0;

	for(unsigned int i=0;i<th.size();i++){
	  ps=tmpps;
	  tmpps=fabs((sin(th[i])*(cos(ph[i])*m_rx+sin(ph[i])*m_ry)+cos(th[i])*m_rz));
	  r = tmpps > ps ? i : r;
	}
	
	return r;
	
      }


      friend ostream& operator<<(ostream& ostr,const Particle& p);
  

    };

  //overload <<
  inline ostream& operator<<(ostream& ostr,const Particle& p){
    ostr<<p.m_x<<" "<<p.m_y<<" "<<p.m_z<<endl;
    return ostr;
  }

  


}

#endif











