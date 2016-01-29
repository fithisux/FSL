/*  tractvols.h

    Tim Behrens, FMRIB Image Analysis Group

    Copyright (C) 2004 University of Oxford  */

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

#ifndef __TRACTVOLS_H_
#define __TRACTVOLS_H_

/////////////////////////////////////////////////////////
//         Class TractVols                             //
/////////////////////////////////////////////////////////

#include "newimage/newimageall.h"
#include <iostream>
#include "stdlib.h"

using namespace std;
using namespace NEWIMAGE;

namespace TRACTVOLS{
  class TractVols
    {
    private:
      volume4D<float> thsamples;
      volume4D<float> phsamples;
      volume4D<float> fsamples;
      bool usef;
    public:
      //constructors::
      TractVols(const bool& usefin=false):usef(usefin){}
      ~TractVols(){}
      
      //Initialise
      void initialise(const string& basename){
	read_volume4D(thsamples,basename+"_thsamples");
	read_volume4D(phsamples,basename+"_phsamples");
	if(usef)
	  read_volume4D(fsamples,basename+"_fsamples");
      }
      

      ColumnVector sample(const float& x,const float& y,const float&z){
	// 	int r_x=(int) MISCMATHS::round(x);
	// 	int r_y=(int) MISCMATHS::round(y);
	// 	int r_z=(int) MISCMATHS::round(z);
	
	////////Probabilistic interpolation
	int cx =(int) ceil(x),fx=(int) floor(x);
	int cy =(int) ceil(y),fy=(int) floor(y);
	int cz =(int) ceil(z),fz=(int) floor(z);
	
	//cerr<<x<<" "<<y<<" "<<z<<" "<<cx<<" "<<cy<<" "<<cz<<" "<<fx<<" "<<fy<<" "<<fz<<endl;
	float pcx,pcy,pcz;
	if(cx==fx)
	  pcx=1;
	else
	  pcx=(x-fx)/(cx-fx);
	
	if(cy==fy)
	  pcy=1;
	else
	  pcy=(y-fy)/(cy-fy);
	
	if(cz==fz)
	  pcz=1;
	else
	  pcz=(z-fz)/(cz-fz);
	
	///////new xyz values from probabilistic interpolation
	int newx,newy,newz; 
	float tmp=rand(); tmp/=RAND_MAX;
	if(tmp>pcx)
	  newx=fx;
	else
	  newx=cx;
	
	tmp=rand(); tmp/=RAND_MAX;
	if(tmp>pcy)
	  newy=fy;
	else
	  newy=cy;
	
	tmp=rand(); tmp/=RAND_MAX;
	if(tmp>pcz)
	  newz=fz;
	else
	  newz=cz;
 
	ColumnVector th_ph_f(3);	
	
	
	float samp=rand(); samp/=RAND_MAX;
	samp=MISCMATHS::round(samp*(thsamples.tsize()-1));
	//float phi = phsamples(r_x,r_y,r_z,samp);
	//float theta = thsamples(r_x,r_y,r_z,samp);
		
	float phi = phsamples(int(newx),int(newy),int(newz),int(samp));
	float theta = thsamples(int(newx),int(newy),int(newz),int(samp));
	
	float f;
	
	if(usef){
	  f = fsamples(int(newx),int(newy),int(newz),int(samp));
	}
	else
	  f=1;
	
	th_ph_f(1)=theta;
	th_ph_f(2)=phi;
	th_ph_f(3)=f;
	return th_ph_f;
      }

      ColumnVector dimensions(){
	ColumnVector dims(3);
	dims << thsamples.xdim() <<thsamples.ydim() << thsamples.zdim();
	return dims;
      }
    };
}

#endif



