/*  Copyright (C) 2000 University of Oxford  */

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

#include "newimagefns.h"
#include "complexvolume.h"
#include <cstdlib>

namespace NEWIMAGE {

// COMPLEX REF
const complexpoint& complexref::operator=(const complexpoint& val){
  *m_real = val.re();
  *m_imag = val.im();
  return(val);
}
  
// COMPLEX POINT
float complexpoint::operator=(const float val)
{
  m_real=val;
  m_imag=0;
  return(val);
}
complexpoint& complexpoint::operator=(const complexpoint& source)
{
  m_real=source.re();
  m_imag=source.im();
  return *this;
}
complexpoint& complexpoint::operator=(const complexref& source){
  m_real = *(source.re_pointer());
  m_imag = *(source.im_pointer());
  return *this;
}
const complexpoint& complexpoint::operator+=(const complexpoint& val)
{
  m_real+=val.re();
  m_imag+=val.im();
  return *this;
}
const complexpoint& complexpoint::operator-=(const complexpoint& val)
{
  m_real-=val.re();
  m_imag-=val.im();
  return *this;
}
const complexpoint& complexpoint::operator*=(const complexpoint& val)
{
  float r2 = (m_real*val.re()) - (m_imag*val.im());
  float i2 = (m_real*val.im()) + (m_imag*val.re());
  m_real = r2;
  m_imag = i2;
  return *this;
}
const complexpoint& complexpoint::operator/=(const complexpoint& val)
{
  float r2 = ((m_real*val.re())+(m_imag*val.im()))/((val.re()*val.re())+(val.im()*val.im()));
  float i2 = ((m_imag*val.re())-(m_real*val.im()))/((val.re()*val.re())+(val.im()*val.im()));
  m_real = r2;
  m_imag = i2;
  return *this;
}
complexpoint complexpoint::operator+(const complexpoint& val) const
{
  complexpoint tmp = *this;
  tmp += val;
  return(tmp);
}
complexpoint complexpoint::operator-(const complexpoint& val) const
{
  complexpoint tmp = *this;
  tmp -= val;
  return(tmp);
}
complexpoint complexpoint::operator*(const complexpoint& val) const
{
  complexpoint tmp = *this;
  tmp *= val;
  return(tmp);
}
complexpoint complexpoint::operator/(const complexpoint& val) const
{
  complexpoint tmp = *this;
  tmp /= val;
  return(tmp);
}
float complexpoint::abs(void) const
{
  return(sqrt((m_real)*(m_real)+(m_imag)*(m_imag)));
}
float complexpoint::phase(void) const
{
  return(atan2(m_imag,m_real));
}

//ostream& complexpoint::operator<<(ostream& s, const complexpoint& val)
//{
//  if(im()>=0.0){
//    return s << re() << " + " << im() << "i";
//  } else {
//    return s << re() << " - " << fabs(im()) << "i";
//  }
//}

// COMPLEX VOLUME
complexvolume::complexvolume(int xsize, int ysize, int zsize)
{
  volume<float> dummy(xsize,ysize,zsize);
  dummy=0.0;
  real=dummy;
  imag=dummy;
}

complexvolume::complexvolume(const complexvolume& source)
{
  real=source.real;
  imag=source.imag;
}
complexvolume::complexvolume(const volume<float>& r, const volume<float>& i)
{
  real=r;
  imag=i;
  if(!samesize(r,i))
    imthrow("Attempted to create complex volume with non-matching sizes",2);
}
complexvolume::complexvolume(const volume<float>& r)
{
  real=r;
  imag=0;
}
complexvolume::~complexvolume()
{
  this->destroy();
}
void complexvolume::destroy()
{
  real.destroy();
  imag.destroy();
}
float complexvolume::operator=(const float val)
{
  real=val;
  imag=0;
  return(val);
}
complexvolume& complexvolume::operator=(const complexvolume& source)
{
  real=source.real;
  imag=source.imag;
  return *this;
}
int complexvolume::copyproperties(const complexvolume& source)
{
  real.copyproperties(source.real);
  imag.copyproperties(source.real);
  return 0;
}
int complexvolume::copydata(const complexvolume& source)
{
  real.copydata(source.real);
  imag.copydata(source.real);
  return 0;
}
const complexvolume& complexvolume::operator+=(const complexpoint& val)
{
  real+=val.re();
  imag+=val.im();
  return *this;
}
const complexvolume& complexvolume::operator-=(const complexpoint& val)
{
  real-=val.re();
  imag-=val.im();
  return *this;
}
const complexvolume& complexvolume::operator*=(const complexpoint& val)
{
  volume<float> r2 = (real*val.re())-(imag*val.im());
  volume<float> i2 = (real*val.im())+(imag*val.re());
  real=r2;
  imag=i2;
  return *this;
}
const complexvolume& complexvolume::operator/=(const complexpoint& val)
{
  volume<float> r2 = ((real*val.re())+(imag*val.im()))/((val.re()*val.re())+(val.im()*val.im()));
  volume<float> i2 = ((imag*val.re())-(real*val.im()))/((val.re()*val.re())+(val.im()*val.im()));
  real = r2;
  imag = i2;
  return *this;
}
const complexvolume& complexvolume::operator+=(const complexvolume& val)
{
  real+=val.real;
  imag+=val.imag;
  return *this;
}
const complexvolume& complexvolume::operator-=(const complexvolume& val)
{
  real-=val.real;
  imag-=val.imag;
  return *this;
}
const complexvolume& complexvolume::operator*=(const complexvolume& val)
{
  volume<float> r2 = (real*val.real)-(imag*val.imag);
  volume<float> i2 = (real*val.imag)+(imag*val.real);
  real=r2;
  imag=i2;
  return *this;
}
const complexvolume& complexvolume::operator/=(const complexvolume& val)
{
  volume<float> r2 = ((real*val.real)+(imag*val.imag))/((val.real*val.real)+(val.imag*val.imag));
  volume<float> i2 = ((imag*val.real)-(real*val.imag))/((val.real*val.real)+(val.imag*val.imag));
  real = r2;
  imag = i2;
  return *this;
}
complexvolume complexvolume::operator+(const complexpoint& val) const
{
  complexvolume tmp=*this;
  tmp += val;
  return(tmp);
}
complexvolume complexvolume::operator-(const complexpoint& val) const
{
  complexvolume tmp=*this;
  tmp -= val;
  return(tmp);
}
complexvolume complexvolume::operator*(const complexpoint& val) const
{
  complexvolume tmp=*this;
  tmp *= val;
  return(tmp);
}
complexvolume complexvolume::operator/(const complexpoint& val) const
{
  complexvolume tmp=*this;
  tmp /= val;
  return(tmp);
}
complexvolume complexvolume::operator+(const complexvolume& val) const
{
  complexvolume tmp=*this;
  tmp += val;
  return(tmp);
}
complexvolume complexvolume::operator-(const complexvolume& val) const
{
  complexvolume tmp=*this;
  tmp -= val;
  return(tmp);
}
complexvolume complexvolume::operator*(const complexvolume& val) const
{
  complexvolume tmp=*this;
  tmp *= val;
  return(tmp);
}
complexvolume complexvolume::operator/(const complexvolume& val) const
{
  complexvolume tmp=*this;
  tmp /= val;
  return(tmp);
}

volume<float> complexvolume::abs(void) const
{
  return(NEWIMAGE::abs(real,imag));
}

volume<float> complexvolume::phase(void) const
{
  return(NEWIMAGE::phase(real,imag));
}

volume<float>& complexvolume::re(void)
{
  return(real);
}

volume<float>& complexvolume::im(void)
{
  return(imag);
}

const volume<float>& complexvolume::re(void) const
{
  return(real);
}

const volume<float>& complexvolume::im(void) const
{
  return(imag);
}

complexvolume complexvolume::extract_slice(int slice) const
{
  volume<float> tempr(real.xsize(),real.ysize(),1);
  volume<float> tempi(real.xsize(),real.ysize(),1);
  
  for(int x=0;x<real.xsize();x++){
    for(int y=0;y<real.ysize();y++){
      tempr(x,y,0) = real(x,y,slice);
      tempi(x,y,0) = imag(x,y,slice);
    }
  }
  complexvolume out(tempr,tempi);
  return(out);
}
void complexvolume::overwrite_slice(const complexvolume& data,int slice)
{
  for(int x=0;x<real.xsize();x++){
    for(int y=0;y<real.ysize();y++){
      real(x,y,slice) = data.re(x,y,0);
      imag(x,y,slice) = data.im(x,y,0);
    }
  }
}

}
