/*  Copyright (C) 1999-2004 University of Oxford  */

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

#include "profile.h"


Profile::Profile()
{
  v.clear();
  lroi = 0;
  rroi = 1;
  maxdef = false;
  mindef = false;
}


Profile::~Profile()
{
}


int Profile::size() const
{
  int counter = 0;
  for (vector<pro_pair>::const_iterator i = v.begin(); i!=v.end(); i++)
    counter++;
  return counter;
}

void Profile::print() const
{
  vector<pro_pair>::const_iterator i;
  for (i = v.begin(); i!=v.end(); i++)
    cout<<(*i).abs<<" : "<<(*i).val<<endl;

}


void Profile::add (const double d, const double t)
{
  //excuse our french.... leaky-leaky!
  //pro_pair * p = new(pro_pair);
  //p->abs = d;
  //p->val = t;
  //v.push_back(*p);

  //corrected code:
  pro_pair p;
  p.abs=d;
  p.val=t;
  v.push_back(p);

  rroi=v.size();
  maxdef = false;
  mindef = false;
}


void Profile::init_roi()
{
  lroi = 0;
  rroi = v.size();
  maxdef = false;
  mindef = false;
}


void Profile::set_lroi(const double abs)
{
  vector<pro_pair >::const_iterator i = v.begin();
  int counter = 0;
  while ((*i).abs < abs && (i++)!=v.end()) counter ++;
  lroi = counter;
  maxdef = false;
  mindef = false;
  if (rroi < lroi) rroi = lroi;
}


void Profile::set_rroi(const double abs)
{
  vector<pro_pair >::const_iterator i = v.end();
  i--;
  int counter = v.size();
  while ((*i).abs > abs && i!=v.begin()) {counter --; i--;}
  rroi = counter;
  maxdef = false;
  mindef = false;
  if (rroi < lroi) lroi = rroi;
}


const double Profile::value(const double d) const 
{
  vector<pro_pair>::const_iterator i = v.begin();
  while ((*i).abs < d && i!=v.end())
    i++ ;
  if (i == v.end())
    {
      cerr<<"out of range"<<endl;
      exit (-1);
    }
  return (*i).val;
}


const double Profile::min() 
{
  if (mindef) return v[amin].val;
  double result = v[lroi].val;
  int abs = lroi;
  for (int i = lroi; i < rroi; i++)
    {
      if (v[i].val < result) {result = v[i].val; abs = i;}
    }
  mindef = true;
  amin = abs;
  return result;
}


const double Profile::max() 
{
  if (maxdef) {return v[amax - 1].val;};
  double result = v[lroi].val;
  int abs = lroi;
  for (int i = lroi; i < rroi; i++)
    {
      if (v[i].val > result) {result = v[i].val; abs = i;}
    }
  maxdef = true;
  amax = abs + 1;
  return result;
}



const double Profile::minabs()
{
  if (mindef) return v[amin].abs;
  else {
    min();
    return v[amin].abs;
  }
}


const double Profile::maxabs()
{
  if (maxdef) return v[amax - 1].abs;
  max();
  return v[amax - 1].abs;
}


const double Profile::threshold(const double d)
{
  return min() + d* (max() - min()); 
}


const double Profile::begin()
{
  return v[lroi].abs;
}


const double Profile::end()
{
  return v[rroi - 1].abs;
}


const double Profile::next_point_over (const double abs, const double thr)
{
  double t = threshold(thr);
  vector<pro_pair >::const_iterator i = v.begin();
  int counter = 0;
  while ((*i).abs < abs && (i++)!=v.end()) counter ++;

  if (i == v.end()) return -500;

  while ((*i).val < t && counter < rroi) {counter ++; i++; if(i == v.end()) return -500;};
  
  if (counter == rroi) return -500;
  else 
    return (v[counter].abs);
}


const double Profile::next_point_under (const double abs, const double thr)
{
  double t = threshold(thr);

  vector<pro_pair >::const_iterator i = v.begin();
  int counter = 0;
  while ((*i).abs < abs && (i++)!=v.end()) counter ++;

  while ((*i).val > t && counter < rroi) {counter ++; i++; if(i == v.end()) return -500;};
  
  if (counter == rroi) return -500;
  else 
    return (v[counter].abs);
}



const double Profile::last_point_under (const double abs, const double thr)
{
  double t = threshold(thr);

  vector<pro_pair >::const_iterator i = v.end();
  i--;
  int counter = v.size();
  while ((*i).abs > abs && i!=v.begin()) {counter --; i--;}

  while (counter > lroi && (*i).val > t && i!=v.begin()) {counter --; i--;};
  
  if ((counter == lroi) | (i==v.begin())) return -500;
  else 
    return (v[counter - 1].abs);
}



const double Profile::last_point_over (const double abs, const double thr)
{
  double t = threshold(thr);

  vector<pro_pair >::const_iterator i = v.end();
  i--;
  int counter = v.size();
  while ((*i).abs > abs && i!=v.begin()) {counter --; i--;}

  while ((*i).val < t && counter > lroi && i!=v.begin()) {counter --; i--;};
  
  if ((counter == lroi) | (i==v.begin())) return -500;
  else 
    return (v[counter - 1].abs);
}







