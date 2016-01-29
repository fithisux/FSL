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

#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <algorithm>

#include "utils/options.h"
#include "newimage/newimageall.h"
#include "meshclass.h"

using namespace std;
using namespace NEWIMAGE;
using namespace Utilities;
using namespace mesh;

void noMoreMemory()
{
  cerr<<"Unable to satisfy request for memory"<<endl;
  abort();
}


vector<float> empty_vector(0, 0);



void draw_segment(volume<short>& image, const Pt& p1, const Pt& p2)
{
  double xdim = (double) image.xdim();
  double ydim = (double) image.ydim();
  double zdim = (double) image.zdim();
  double mininc = min(xdim,min(ydim,zdim)) * .5;

  Vec n = p1 - p2;
  double d = n.norm();
  n.normalize();

  for (double i=0; i<=d; i+=mininc)
    {
      Pt p = p2 + i* n;
      image((int) floor((p.X)/xdim +.5),(int) floor((p.Y)/ydim +.5),(int) floor((p.Z)/zdim +.5)) = 1;
    }
}


volume<short> draw_mesh(const volume<short>& image, const Mesh &m)
{
  double xdim = (double) image.xdim();
  double ydim = (double) image.ydim();
  double zdim = (double) image.zdim();
  double mininc = min(xdim,min(ydim,zdim)) * .5;
  volume<short> res = image;
  for (list<Triangle*>::const_iterator i = m._triangles.begin(); i!=m._triangles.end(); i++)
    {
      Vec n = (*(*i)->get_vertice(0) - *(*i)->get_vertice(1));
      double d = n.norm();
      n.normalize();

      for (double j=0; j<=d; j+=mininc)
	{
	  Pt p = (*i)->get_vertice(1)->get_coord()  + j* n;
	  draw_segment(res, p, (*i)->get_vertice(2)->get_coord());
	} 
    }
  return res;
}

volume<short> make_mask_from_mesh(const volume<short> & image, const Mesh& m)
{
  double xdim = (double) image.xdim();
  double ydim = (double) image.ydim();
  double zdim = (double) image.zdim();

  volume<short> mask = image;
  
  int xsize = mask.xsize();
  int ysize = mask.ysize();
  int zsize = mask.zsize();
  
  vector<Pt> current;
  current.clear();
  Pt c(0., 0., 0.);
  for (vector<Mpoint *>::const_iterator it=m._points.begin(); it!=m._points.end(); it++)
    c+=(*it)->get_coord();

  c*=(1./m._points.size());
  c.X/=xdim; c.Y/=ydim; c.Z/=zdim;

  current.push_back(c);

  while (!current.empty())
    {
      Pt pc = current.back();
      int x, y, z;
      x=(int) pc.X; y=(int) pc.Y; z=(int) pc.Z;
      current.pop_back();
      mask.value(x, y, z) = 1;
      if (0<=x-1 && mask.value(x-1, y, z)==0) current.push_back(Pt(x-1, y, z));
      if (0<=y-1 && mask.value(x, y-1, z)==0) current.push_back(Pt(x, y-1, z));
      if (0<=z-1 && mask.value(x, y, z-1)==0) current.push_back(Pt(x, y, z-1));
      if (xsize>x+1 && mask.value(x+1, y, z)==0) current.push_back(Pt(x+1, y, z));
      if (ysize>y+1 && mask.value(x, y+1, z)==0) current.push_back(Pt(x, y+1, z));
      if (zsize>z+1 && mask.value(x, y, z+1)==0) current.push_back(Pt(x, y, z+1)); 
    }
  return mask;
}

int main(int argc, char *argv[]) {

  
  if (argc != 3 && argc != 4) {
    cerr<<"Usage : drawmesh  [volume.hdr] [mesh.off] (-m for the mask)"<<endl;
    exit (-1);
  }
  
  bool mesh = false;
  if (argc == 4)
    {
      string ismesh=argv[3];
      if (ismesh.compare("-m") == 0)
	mesh = true;
      else 
	{
	  cerr<<"Usage : drawmesh  [volume.hdr] [mesh.off] (-m for the mask)"<<endl;
	  exit (-1);
	}
    }

  string volumename=argv[1];
  string meshname=argv[2];
  

  string out = meshname;
  if (out.find(".off")!=string::npos) out.erase(out.find(".off"), 4);

  string in = volumename;
  if (in.find(".hdr")!=string::npos) in.erase(in.find(".hdr"), 4);
  if (in.find(".img")!=string::npos) in.erase(in.find(".hdr"), 4);
  if (out == "default__default") {out=in+"_brain";}
  

  //set a memory hanlder that displays an error message
  set_new_handler(noMoreMemory);


  //the real program
  volume<short> testvol;
  volume<double> testvol2;
  
  if (read_volume(testvol2,in.c_str())<0)  return -1;

  copyconvert(testvol2, testvol);

  testvol = 0;

  Mesh m;
  m.load(meshname);  
  
  cout<<"saving volume"<<endl;
  string outlinestr = in+"_outline";
  volume<short> outline = draw_mesh(testvol, m);
  if (save_volume(outline, outlinestr.c_str())<0)  return -1;

  if (mesh) 
    {
      outline = make_mask_from_mesh(outline, m);
      if (save_volume(outline, (in+"_mask"))<0)  return -1;
    }

  return 0;
  
}





