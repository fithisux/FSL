
/*  BET - Brain Extraction Tool

    BETv1 Steve Smith
    BETv2 Mickael Pechaud, Mark Jenkinson, Steve Smith
    FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

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
#include "meshclass/meshclass.h"
#include "betsurf.h"

using namespace std;
using namespace NEWIMAGE;
using namespace Utilities;
using namespace mesh;

int infxm;
int infym;
int infzm;
double l;



void noMoreMemory()
{
  cerr<<"Unable to satisfy request for memory"<<endl;
  abort();
}

string title="BETSURF (BET Surface Finder) v2.1 - FMRIB Analysis Group, Oxford";
string examples=" betsurf          [options] <t1> <t2> <bet_mesh.off> <t1_to_standard.mat> <output>\n betsurf --t1only [options] <t1>      <bet_mesh.off> <t1_to_standard.mat> <output>";

Option<bool> help(string("-h,--help"), false, 
		     string("displays this help, then exits"), 
		     false, no_argument);

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);

Option<bool> t1only(string("-1,--t1only"), false, 
		    string("extraction with t1 only"), false, no_argument);

Option<bool> outline(string("-o,--outline"), false, 
		     string("generates all surface outlines"), 
		     false, no_argument);

Option<bool> mask(string("-m,--mask"), false, 
		  string("generates binary masks from the meshes"), 
		  false, no_argument);

Option<bool> skullmask(string("-s,--skullmask"), false,
		       string("generates skull binary mask"), false, no_argument);

//Option<bool> csf(string("-c,--csf"), false,
//		 string("generates an approximation of csf"), false, no_argument);

Option<int> precision(string("-p,--increased_precision"), 0,
		 string("~<int>\tretessellates the meshes the indicated number of times"), false, requires_argument);


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


void draw_mesh(volume<short>& image, const Mesh &m)
{
  double xdim = (double) image.xdim();
  double ydim = (double) image.ydim();
  double zdim = (double) image.zdim();
  double mininc = min(xdim,min(ydim,zdim)) * .5;

  for (list<Triangle*>::const_iterator i = m._triangles.begin(); i!=m._triangles.end(); i++)
    {
      Vec n = (*(*i)->get_vertice(0) - *(*i)->get_vertice(1));
      double d = n.norm();
      n.normalize();

      for (double j=0; j<=d; j+=mininc)
	{
	  Pt p = (*i)->get_vertice(1)->get_coord()  + j* n;
	  draw_segment(image, p, (*i)->get_vertice(2)->get_coord());
	} 
    }

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

double standard_step_of_computation(const volume<float> & image, Mesh & m, const int iteration_number, const double E,const double F, const float addsmooth, const float speed, const int nb_iter, const int id, const int od, const bool vol, const volume<short> & mask){
  double xdim = image.xdim();
  double ydim = image.ydim();
  double zdim = image.zdim();
  
  if (nb_iter % 50 == 0)
    {
      double l2 = 0;
      int counter = 0;
      for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ )
	{
	  counter++;
	  l2 += (*i)->medium_distance_of_neighbours();
	}
      l = l2/counter;
    }
  if (nb_iter % 100 == 0)
    {
      for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++)
	{
	  Vec n = (*i)->local_normal();
	  Pt point = (*i)->get_coord();
	  Pt ipoint(point.X/xdim, point.Y/ydim, point.Z/zdim);
	  Vec in(n.X/xdim, n.Y/ydim, n.Z/zdim);
	  double max = 0;
	  Pt c_m1 = ipoint + (-1) * in;
	  double current = image.interpolate((c_m1.X),(c_m1.Y),(c_m1.Z));
	  for (double i2 = 1; i2 < 150; i2+=2)
	    {
	      if (max > .1) break;
	      Pt c_p = ipoint + i2 * in;
	      double tmpp = image.interpolate((c_p.X),(c_p.Y),(c_p.Z));
	      double tmp = (tmpp - current) * 100;
	      max = Max(max, tmp);
	      current = tmpp;
	      if (tmpp > .1) {max = 1; break;}
	    }

	  if (max < .1)
	    {   
              //There is a problem here for precision mode, since with the copy, no guarantee that data is non-zero size
              //even if mesh.cpp operator = is modified to copy data, after retesselate "new" points will have zero size data member
	      if ( (*i)->data.size() ) (*i)->data.pop_back();
	      (*i)->data.push_back(1);
	    }
	  else
	    {
	      if ( (*i)->data.size() ) (*i)->data.pop_back();
	      (*i)->data.push_back(0);
	    }
	}
    }

  for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++)
    {
      Vec sn, st, u1, u2, u3, u;
      double f2, f3=0;
      
      Vec n = (*i)->local_normal();
      Vec dv = (*i)->difference_vector();
      
      double tmp = dv|n;
      sn = n * tmp;
      st = dv - sn;
      
      u1 = st*.5;
      
      double rinv = (2 * fabs(sn|n))/(l*l);
      
      f2 = (1+tanh(F*(rinv - E)))*0.5;
      
      u2 = f2 * sn * addsmooth;
      
      if ((*i)->data.back() == 0)
	{
	  //main term of skull_extraction
	  {
	    Pt point = (*i)->get_coord();
	    Pt ipoint(point.X/xdim, point.Y/ydim, point.Z/zdim);
	    Vec in(n.X/xdim, n.Y/ydim, n.Z/zdim);
	    
	    Pt c_m = ipoint + (-1.) * in;
	    Pt c_p = ipoint + 1. * in;
	    
	    double tmp = image.interpolate((c_p.X ),( c_p.Y),(c_p.Z));
	    double gradient = tmp - image.interpolate((c_m.X),(c_m.Y), (c_m.Z));
	    
	    double tmp2 = gradient*100;
	    f3 = max(-1., min(tmp2, 1.));
	    if (tmp2 >= 0 && tmp2 < .1 && tmp < .1 ) f3 = speed;
	    
	    if (vol) 
	      {
		double tmpvol = mask.interpolate((ipoint.X ),(ipoint.Y),(ipoint.Z));
		if (tmpvol > .0)
		  {
		    f3 = Max(Max(tmpvol*.5, .1), f3); 
		    f2 = 0;
		  }
	      }
	    
	    
	  }
	}
      else 
	{
	  f3 = 0;
	  Pt point = (*i)->get_coord();
	  Pt ipoint(point.X/xdim, point.Y/ydim, point.Z/zdim);
	  double tmpvol = mask.interpolate((ipoint.X ),(ipoint.Y),(ipoint.Z));
	  if (tmpvol > .0)
	    {
	      f3 = Max(tmpvol*.5, .1); 
	      f2 = 0;
	    }
	}

      u3 = .05 * f3 * n;
      
      u = u1 + u2 + u3;
      
      (*i)->_update_coord = (*i)->get_coord() + u;
    }

  m.update();
  
  return (0); 
}

//extracts profiles in the area around the eyes
vector<double> t1only_special_extract(const volume<float> & t1, const Pt & point, const Vec & n) {
  vector<double> resul;
  resul.clear();
  bool output = true; 
  const double INNER_DEPTH = 3;
  const double OUTER_DEPTH = 100;
  
  Profile pt1;
  for (double d = -INNER_DEPTH; d < OUTER_DEPTH; d+=.5)
    {
      Pt c = point + d * n;
      double tmp1 = t1.interpolate(c.X, c.Y, c.Z);
      pt1.add (d, tmp1);
    }
  pt1.init_roi();

  //outer skin
  double outskin = pt1.last_point_over(pt1.end(), .2);
  double check = pt1.last_point_over(outskin - 1.5, .2);
  if (outskin - check > 2) output = false;
  
  pt1.set_rroi(outskin);
  
  double inskull = pt1.next_point_under(-INNER_DEPTH, .25);
  if (inskull > 5) inskull = 0;
  
  double outskull = pt1.next_point_over(inskull, .35);
  
  resul.push_back(inskull);
  resul.push_back(outskull);
  resul.push_back(outskin);
  
  return resul;
}

vector<double> t1only_co_ext(const volume<float> & t1, const Pt & point, const Vec & n) {
  vector<double> resul;
  resul.clear();
  bool output = true; 
  bool alloutput = true;
  const double INNER_DEPTH = 3;
  const double OUTER_DEPTH = 60;
  
  Profile pt1;
  for (double d = -INNER_DEPTH; d < OUTER_DEPTH; d+=.5)
    {
      Pt c = point + d * n;
      double tmp1 = t1.interpolate(c.X, c.Y, c.Z);
      pt1.add (d, tmp1);
    }
  pt1.init_roi();

  //outer skin

  double outskin = pt1.last_point_over(pt1.end(), .2);
  double check = pt1.last_point_over(outskin - 1.5, .2);
  if (outskin - check > 2) outskin = check;

  const double OUTER_SKIN = outskin;  

  pt1.set_rroi(OUTER_SKIN);

  double inskull = pt1.next_point_under(pt1.begin(), .25);

  pt1.set_lroi(inskull);

  if (alloutput)
    {
      //outer skull
      
      //starting from the skin
      double outskull2 = 0;
      outskull2 = pt1.last_point_over(outskin, .75); 

      outskull2 = pt1.last_point_under(outskull2, .20);

      //starting from the brain
      double minabs = pt1.next_point_under(pt1.begin(), .30);
      if (minabs == -500) output = false;
      minabs = Max(minabs, -INNER_DEPTH);

      double localminabs = minabs;//0
      if (output)
	{
	  bool stop = false;

	  const double lowthreshold = pt1.threshold(.15);

	  const double upthreshold = pt1.threshold(.60);
	  int test = 0;
	  for (vector<pro_pair>::const_iterator i = pt1.v.begin(); i != pt1.v.end(); i++)
	    {
	      if (!stop && (*i).abs>=minabs && i!=pt1.v.end() && i!=pt1.v.begin())
		{
		  if ((*i).val>upthreshold) stop = true; //avoid climbing skin
		  if ((*i).val>lowthreshold) test++;
		  if ((*i).val<lowthreshold) test--;
		  if (test < 0) test=0;
		  if (test == 12) {stop = true;}
		  if ((*i).val<lowthreshold) 
		    {
		      localminabs = (*i).abs;
		    }
		}
	    }
	}

      double outskull = pt1.next_point_over(localminabs, .15);//.20
      
      if (outskull2 - outskull < -2 || outskull2 - outskull >= 2) {output = false;}
      
      if (outskin - outskull2 < 2) output = false;
      
      if (output) 
	{
	  resul.push_back(inskull);
	  resul.push_back(outskull2);
	  resul.push_back(outskin);
	}
      else resul.push_back(outskin);
    }
  
  return resul;
}


//writes externall skull computed from image on output.
void t1only_write_ext_skull(volume<float> & output_inskull, volume<float> & output_outskull, volume<float> & output_outskin, const volume<float> & t1, const Mesh & m, const trMatrix & M) {
  int glob_counter = 0;
  int rem_counter = 0;

  const double xdim = t1.xdim();
  const double ydim = t1.ydim();
  const double zdim = t1.zdim();

  double imax = t1.max();
  if (imax == 0) imax = 1;

  volume<short> meshimage;
  copyconvert(t1, meshimage);
  meshimage = 0;
  draw_mesh(meshimage, m);

  for (vector<Mpoint*>::const_iterator i = m._points.begin(); i != m._points.end(); i++)
    {
      (*i)->data.clear();

      double max_neighbour = 0;
      const Vec normal = (*i)->local_normal();
      const Vec n = Vec(normal.X/xdim, normal.Y/ydim, normal.Z/zdim);      
      
      for (list<Mpoint*>::const_iterator nei = (*i)->_neighbours.begin(); nei != (*i)->_neighbours.end(); nei++)
	max_neighbour = Max(((**i) - (**nei)).norm(), max_neighbour); 
      
      max_neighbour = ceil((max_neighbour)/2);
      
      const Pt mpoint((*i)->get_coord().X/xdim,(*i)->get_coord().Y/ydim,(*i)->get_coord().Z/zdim);
      for (int ck = (int)floor(mpoint.Z - max_neighbour/zdim); ck <= (int)floor(mpoint.Z + max_neighbour/zdim); ck++)
	for (int cj = (int)floor(mpoint.Y - max_neighbour/ydim); cj <= (int)floor(mpoint.Y + max_neighbour/ydim); cj++)
	  for (int ci = (int)floor(mpoint.X - max_neighbour/xdim); ci <= (int)floor(mpoint.X + max_neighbour/xdim); ci++)
	    {
	      bool compute = false;
	      const Pt point(ci, cj, ck);
	      const Pt realpoint(ci*xdim, cj*ydim, ck*zdim);
	      if (meshimage(ci, cj, ck) == 1) 
		{
		  double mindist = 10000;
		  for (list<Mpoint*>::const_iterator nei = (*i)->_neighbours.begin(); nei != (*i)->_neighbours.end(); nei++)
		    mindist = Min(((realpoint) - (**nei)).norm(), mindist); 
		  if (mindist >= ((realpoint) - (**i)).norm()) compute = true;
		}
	    

	      if (compute)
		{
		  glob_counter ++;
		  vector<double> val;
		  if (!special_case(realpoint, normal, M))
		    val = t1only_co_ext(t1, point, n);
		  else
		    {
		      val = t1only_special_extract(t1, point, n);
		    }

		  if (val.size() == 3)
		    {
		      Pt opoint(point.X, point.Y, point.Z); 
		      Vec on(n.X, n.Y, n.Z); 
		      Pt c0 = opoint + val[0]*on;
		      Pt c1 = opoint + val[1]*on;
		      Pt c2 = opoint + val[2]*on;

		      output_inskull((int)floor(c0.X + .5) + infxm,(int) floor(c0.Y + .5) + infym,(int) floor(c0.Z + .5) + infzm) +=1; 
		      output_outskull((int)floor(c1.X + .5) + infxm,(int) floor(c1.Y + .5) + infym,(int) floor(c1.Z + .5) + infzm)+=1; 
		      output_outskin((int)floor(c2.X + .5) + infxm,(int) floor(c2.Y + .5) + infym,(int) floor(c2.Z + .5) + infzm) +=1; 
		    }
		  else {
		    rem_counter++;

		    if (val.size()==1)
		      {
			Pt opoint(point.X, point.Y, point.Z); 
			Vec on(n.X, n.Y, n.Z); 
			Pt c0 = opoint + val[0]*on;
			
			output_outskin((int)floor(c0.X + .5) + infxm,(int) floor(c0.Y + .5) + infym,(int) floor(c0.Z + .5) + infzm) +=1; 
		      }
		  }
		}      
	    }
    }
  if (verbose.value())
    {
      cout<<" nb of profiles : "<<glob_counter<<endl;
      cout<<" removed profiles : "<<100. * rem_counter/(double) glob_counter<<"%"<<endl;
    }
}

int t1only_main(int argc, char *argv[], int nb_pars, OptionParser & options){
  
  if (argc - nb_pars < 4)
    {
      cerr<<"too few arguments"<<endl;
      options.usage(); return -1;
    }
  
  int count_arg = nb_pars;
  
  const string inputt1(argv[count_arg]);
  count_arg++;
  const string mesh(argv[count_arg]);
  count_arg++;
  const string matrix(argv[count_arg]);
  count_arg++;
  const string outputstr(argv[count_arg]);


  //load the mesh
  Mesh m;
  m.load(mesh);

  //load the matrix
  trMatrix M;
  ifstream f(matrix.c_str());
  if (f.is_open())
    {
      f>>M.m11>>M.m12>>M.m13>>M.m14>>M.m21>>M.m22>>M.m23>>M.m24>>M.m31>>M.m32>>M.m33>>M.m34>>M.m41>>M.m42>>M.m43>>M.m44;
      f.close();
    }
  else {cerr<<"unable to open "<<matrix<<endl;return -1;}

  //load the volume
  volume<float> t1;

  if (read_volume(t1,inputt1.c_str())<0)  return -1;

  t1.setinterpolationmethod(trilinear);

  const double xdim = t1.xdim();
  const double ydim = t1.ydim();
  const double zdim = t1.zdim();
  const int xsize = t1.xsize();
  const int ysize = t1.ysize();
  const int zsize = t1.zsize();

  //founding brain robustmax (useful if skull is too bright)
  double thr;
  
  {
    volume<short> brain_mask;
    copyconvert(t1, brain_mask);
    brain_mask = 0;
    draw_mesh(brain_mask, m);
    brain_mask = make_mask_from_mesh(brain_mask, m);
    vector<double> brain_hist;
    for (int k = 0; k < zsize; k++)
      for(int j = 0; j < ysize; j++)
	for(int i = 0; i < xsize; i++)
	  if (brain_mask.value(i, j, k) == 1) brain_hist.push_back(t1.value(i, j, k));
    
    int size = brain_hist.size();
    int e98 = (int) ceil(.98 * size);
    nth_element(brain_hist.begin(), brain_hist.begin() + e98 - 1, brain_hist.end());
    thr = brain_hist[e98 - 1];
    thr *= 1.1;
  }

  for (int k = 0; k < zsize; k++)
    for(int j = 0; j < ysize; j++)
      for(int i = 0; i < xsize; i++)
	if (t1.value(i, j, k) > thr) t1.value(i, j, k) = thr;




  infxm=0; infym=0; infzm=0;
  int infxp=0, infyp=0, infzp=0;
  //checking if the mesh is inside the volume
  for (vector<Mpoint *>::const_iterator p = m._points.begin(); p != m._points.end(); p++)
    {
      if ((*p)->get_coord().X/xdim < 3) infxm = Max(infxm, (int) ceil(3 - (*p)->get_coord().X/xdim));
      if ((*p)->get_coord().Y/ydim < 3) infym = Max(infym, (int) ceil(3 - (*p)->get_coord().Y/ydim));
      if ((*p)->get_coord().Z/zdim < 3) infzm = Max(infzm, (int) ceil(3 - (*p)->get_coord().Z/zdim));
      if (xsize - (*p)->get_coord().X/xdim< 3) infxp = Max(infxp, (int) ceil(3 - xsize + (*p)->get_coord().X/xdim));
      if (ysize - (*p)->get_coord().Y/ydim< 3) infyp = Max(infyp, (int) ceil(3 - ysize + (*p)->get_coord().Y/ydim));
      if (zsize - (*p)->get_coord().Z/zdim < 3) infzp = Max(infzp, (int) ceil(3 - zsize + (*p)->get_coord().Z/zdim));
    }

  infxp += 1;
  infxm += 1;
  infyp += 1;
  infym += 1;
  infzp += 1;
  infzm += 1;

  volume<float> write0(xsize + infxm + infxp, ysize + infym + infyp, zsize + infzm + infzp);
  write0.setxdim(xdim);
  write0.setydim(ydim);
  write0.setzdim(zdim);
  //write0.setorigin(infxm * xdim, infym * ydim, infzm * zdim);

  write0 = 0;

  volume<float> write2;
  volume<float> write1;
  write1 = write2 = write0;

  if (verbose.value()) cout<<"extracting profiles"<<endl;

  t1only_write_ext_skull(write0, write1, write2, t1, m, M);

  const double rmin=3.33;
  const double rmax=10;
  const double E = (1/rmin + 1/rmax)/2.;
  const double F = 6./(1/rmin - 1/rmax);

  if (verbose.value()) cout<<"blurring"<<endl;

  const double blurring0 = 3;
  const double blurring1 = 3;
  const double blurring2 = 3;

  const double smoothness0 = 1;
  const double smoothness1 = 1;
  const double smoothness2 = 1;

  write0 = smooth(write0, blurring0);

  write2 = smooth(write2, blurring2);

  const int nb_iter0 = 800;
  const int nb_iter1 = 800;
  const int nb_iter2 = 1500;

  if (verbose.value()) cout<<"computation"<<endl;

  for (vector<Mpoint*>::iterator i = m._points.begin(); i!= m._points.end(); i++)
    {
      (*i)->data.push_back(0);
    }
  
  m.translation(xdim * infxm, ydim * infym, zdim * infzm);

  if (verbose.value()) cout<<" inner skull"<<endl;
  for (int c = 0; c < nb_iter0; c++)
    standard_step_of_computation(write0, m, c, E, F, smoothness0, .5, c);

  Mesh mprecise;
  if (precision.value() > 0)
    {
      mprecise = m;
      if (verbose.value()) cout << "  increased precision" << endl;
      for (int i = 0; i < precision.value(); i++)
        mprecise.retessellate();
      for (int c = 0; c < 100; c++)
      standard_step_of_computation(write0, mprecise, c, E, F, smoothness0, .5, c);
      
    }
  
  Mesh realmesh;
    if (precision.value() == 0)
      realmesh = m;
    else realmesh = mprecise;

  realmesh.translation(- xdim * infxm, - ydim * infym, - zdim * infzm);
  realmesh.save((outputstr+"_inskull_mesh.off").c_str());


  volume<short> maskinskull;
  volume<short> maskinskull2;
  maskinskull2 = 0;

  {
    copyconvert(write0, maskinskull);
    draw_mesh(maskinskull, m);
    maskinskull = make_mask_from_mesh(maskinskull, m);
    
    volume<short> output1;
    copyconvert(t1, output1);
    output1 = 0;

    if (outline.value() | mask.value()| skullmask.value())
      {
	draw_mesh(output1, realmesh);
      }
    
    
    if (outline.value())
      {
	if (verbose.value()) cout<<"  outline"<<endl;
	if (save_volume(output1, (outputstr+"_inskull_mesh").c_str())<0)  return -1;
      }
      
    if (mask.value() | skullmask.value())
      {
	volume<short> smask = make_mask_from_mesh(output1, realmesh);
	if (skullmask.value())
	  {
	    maskinskull2 = smask;
	  }
	if (mask.value())
	  {  
	    if (verbose.value()) cout<<"  mask"<<endl;	
	    if (save_volume(smask, (outputstr+"_inskull_mask").c_str())<0)  return -1;
	  }
      }
  }
  
  
  {
    int xsiz = write1.xsize();
    int ysiz = write1.ysize();
    int zsiz = write1.zsize();
    for (int k = 0; k < zsiz; k++)
      for (int j = 0; j < ysiz; j++)
	for (int i = 0; i < xsiz; i++)
	  if (write1.value(i, j, k) > 0 && maskinskull.value(i, j, k) == 1)
	    {
	      write1.value(i, j, k) = 0;
	    }
  }

  write1 = smooth(write1, blurring1);


  if (verbose.value()) cout<<" outer skull"<<endl;
  for (int c = 0; c < nb_iter1; c++)
    standard_step_of_computation(write1, m, c, E, F, smoothness1, .5, c, 5, 15, true, maskinskull);

  maskinskull.destroy();

  if (precision.value() > 0)
    {
      if (verbose.value()) cout<<"  increased precision"<<endl;
      mprecise = m;
      for (int i = 0; i < precision.value(); i++)
	mprecise.retessellate();
      for (int c = 0; c < 100; c++)
	standard_step_of_computation(write1, mprecise, c, E, F, smoothness1, .5, c, 5, 15);
      
    }

  write1.destroy();
  
  Mesh realmesh2;
  if (precision.value() == 0)
    realmesh2 = m;
  else realmesh2 = mprecise;
  
  realmesh2.translation(- xdim * infxm, - ydim * infym, - zdim * infzm);
  realmesh2.save((outputstr+"_outskull_mesh.off").c_str());

  volume<short> maskoutskull;
  maskoutskull = 0;
  volume<short> maskoutskull2;
  maskoutskull2 = 0;
  volume<short> meshoutskull;
  meshoutskull = 0;  

  {
    copyconvert(write0, maskoutskull);
    write0.destroy();
    draw_mesh(maskoutskull, m);
    maskoutskull = make_mask_from_mesh(maskoutskull, m);
    
    volume<short> output1;
    copyconvert(t1, output1);
    output1 = 0;

    if (mask.value() | outline.value() | skullmask.value())
      draw_mesh(output1, realmesh2);
    
    if (skullmask.value())
      meshoutskull = output1;
    
    if (outline.value())
      {
	if (verbose.value()) cout<<"  outline"<<endl;
	if (save_volume(output1, (outputstr+"_outskull_mesh").c_str())<0)  return -1;
      }
    
    if (mask.value()|skullmask.value())
      {      
	volume<short> smask = make_mask_from_mesh(output1, realmesh2);
	
	if (skullmask.value())
	  {
	    maskoutskull2 = smask;
	  }
	if (mask.value())
	  {
	    if (verbose.value()) cout<<"  mask"<<endl;
	    if (save_volume(smask, (outputstr+"_outskull_mask").c_str())<0)  return -1;
	  }
      }
  }
  
  if (skullmask.value())
    {
      volume<short> smask = maskoutskull2;
      smask = 0;
      int xsize = t1.xsize(), ysize = t1.ysize(), zsize = t1.zsize();
      if (verbose.value()) cout<<" skull mask"<<endl;
      for (int k = 0; k < zsize; k++)
	for(int j = 0; j < ysize; j++)
	  for(int i = 0; i < xsize; i++)
	    smask.value(i, j, k) = Min (1, maskoutskull2.value(i, j, k) * (1 - maskinskull2.value(i, j, k)) + meshoutskull.value(i, j, k));
      
      if (save_volume(smask, (outputstr+"_skull_mask").c_str())<0) return -1;
    }

  meshoutskull.destroy();
  maskoutskull2.destroy();
  maskinskull2.destroy();
  
  //computing the starting mesh for outer skull
  if (verbose.value()) cout<<" skin"<<endl;
  
  Mesh m2;
  make_mesh_from_icosa(5, m2);

  for (vector<Mpoint*>::iterator i = m2._points.begin(); i!= m2._points.end(); i++)
    {
      (*i)->data.push_back(0);
    }
  
  Pt p;
  int counter=0;
  for (vector<Mpoint*>::const_iterator i = m._points.begin(); i!= m._points.end(); i++)
    {
      counter++;
      p+=(*i)->get_coord();
    }
  p*=(1./counter);

  double radius = 0;
  for (vector<Mpoint*>::const_iterator i = m._points.begin(); i!= m._points.end(); i++)
    radius+=((*i)->get_coord()-p).norm();

  radius/=(counter);
  radius*=.75;

  m2.rescale(radius);
  m2.translation(p.X, p.Y, p.Z);

  for (int c = 0; c < nb_iter2; c++)
    standard_step_of_computation(write2, m2, c, E, F, smoothness2, 1.5, c, 5, 15, true, maskoutskull);

  maskoutskull.destroy();
  
  if (precision.value() > 0)
    {
      if (verbose.value()) cout<<"  increased precision"<<endl;
      for (int i = 0; i < precision.value(); i++)
	m2.retessellate();
      for (int c = 0; c < 100; c++)
	standard_step_of_computation(write2, m2, c, E, F, smoothness2, 1.5, c, 5, 15);
      
    }

  write2.destroy();

  m2.translation(-xdim * infxm, -ydim * infym, -zdim * infzm);
  m2.save((outputstr+"_outskin_mesh.off").c_str());

  if (outline.value() | mask.value())
    {
      volume<short> output1;
      copyconvert(t1, output1);
      output1 = 0;

      draw_mesh(output1, m2);

      if (outline.value())
	{
	  if (verbose.value()) cout<<"  outline"<<endl;
	  if (save_volume(output1, (outputstr+"_outskin_mesh").c_str())<0)  return -1;  
	}
      if (mask.value())
	{
	  if (verbose.value()) cout<<"  mask"<<endl;
	  volume<short> mask = make_mask_from_mesh(output1, m);
	  if (save_volume(mask, (outputstr+"_outskin_mask").c_str())<0)  return -1;
	}
    }
  

  return 0;
}


//extracts profiles in the area around the eyes
vector<double> special_extract(const volume<float> & t1, const volume<float> & t2, const Pt & point, const Vec & n, volume<short> & csfvolume) {
  vector<double> resul;
  resul.clear();
  bool alloutput = true;
  const double INNER_DEPTH = 3;
  const double OUTER_DEPTH = 100;
  
  Profile pt1;
  Profile pt2;
  for (double d = -30; d < OUTER_DEPTH; d+=.5)
    {
      Pt c = point + d * n;
      double tmp1 = t1.interpolate(c.X, c.Y, c.Z);
      double tmp2 = t2.interpolate(c.X, c.Y, c.Z);
      pt1.add (d, tmp1);
      pt2.add (d, tmp2);
    }
  pt1.init_roi();
  pt2.init_roi();

  pt1.set_lroi(-INNER_DEPTH);
  pt2.set_lroi(-INNER_DEPTH);

  double inskull=pt1.next_point_under(pt1.begin(), .25);
  double inskullt1 = pt2.next_point_under(inskull, .30);
  if (inskullt1 > inskull) inskull = inskullt1;

  if (inskull > 15) inskull = 0;


  //outer skin
  double outskin = pt1.last_point_over(pt1.end(), .2);
  double check = pt1.last_point_over(outskin - 1.5, .2);
  if (outskin - check > 2) alloutput = false;

  pt1.set_rroi(outskin);

  double val = pt1.next_point_under(inskull, .25);

  double outskull = pt1.next_point_over(val, .35);

  if (val - inskull >= 7) outskull = inskull + .5; 
  if (outskull - inskull < 0) outskull = inskull + .5;

  if (alloutput)
    {
      resul.push_back(outskin);
      resul.push_back(outskull);
      resul.push_back(inskull);
    }
  
  //computing csf
  pt1.init_roi();
  pt1.set_rroi(outskull - .5);
  pt2.init_roi();
  pt2.set_rroi(outskin);
  double end_csf = pt2.last_point_over(inskull, .85);
  if (end_csf != -500)
    {
      for (double d = -30; d <= end_csf; d+=.5)
	if (pt2.value(d) > pt2.threshold(.85) & pt1.value(d) < pt1.threshold(.45))
	  {
	    Pt csf = point + d*n;
	    csfvolume((int)floor(csf.X + .5), (int)floor(csf.Y + .5), (int)floor(csf.Z + .5))=5;
	  }
      
    }
  
  return resul;
}

vector<double> co_ext(const volume<float> & t1, const volume<float> & t2, const Pt & point, const Vec & n, volume<short> & csfvolume) {

  vector<double> resul;
  resul.clear();
  bool output = true; 
  bool alloutput = true;
  const double INNER_DEPTH = 3; 
  const double OUTER_DEPTH = 60;
  
  Profile pt1;
  Profile pt2;
  for (double d = -30; d < OUTER_DEPTH; d+=.5)
    {
      Pt c = point + d * n;
      double tmp1 = t1.interpolate(c.X, c.Y, c.Z);
      double tmp2 = t2.interpolate(c.X, c.Y, c.Z);
      pt1.add (d, tmp1);
      pt2.add (d, tmp2);
    }
  pt1.init_roi();
  pt2.init_roi();

  pt1.set_lroi(-INNER_DEPTH);
  pt2.set_lroi(-INNER_DEPTH);

  //outer skin
  double outskin = pt1.last_point_over(pt1.end(), .2);
  double check = pt1.last_point_over(outskin - 1.5, .2);
  if (outskin - check > 2) outskin = check;


  double outskin2 = pt2.last_point_over(pt2.end(), .25);
  if (fabs(outskin - outskin2) > 5) 
    {
      //artefact might be present ...
      bool b = true;
      if (outskin < outskin2) b = false;
      double m = Min(outskin, outskin2) + 15;
      double m2;
      if (b) 
	{
	  m2 = pt1.last_point_over(m, .2);
	  if (fabs(m2 - outskin2) > 10) alloutput = false;
	  else outskin = m2;
	}
      else 
	{
	  m2 = pt2.last_point_over(m, .25);
	  if (fabs(m2 - outskin) > 10) alloutput = false;
	  else outskin = outskin;
          // The line above is very suspicious
	}
    }

  const double OUTER_SKIN = outskin;

  pt1.set_rroi(OUTER_SKIN);
  pt2.set_rroi(OUTER_SKIN);

  if (alloutput)
    {
      //inner skull & beginning of outer skull search
      double OUTER_SKULL_SEARCH_BEGIN = 0;
      double inskull;
      {
	inskull = pt1.next_point_under(pt2.begin(), .25);
	double inskullt1 = pt2.next_point_under(inskull, .30); 
	if (inskullt1 > inskull) inskull = inskullt1;
	pt2.set_lroi(inskull);
	vector<double> gaps;
	double sum = -50;
	double t = inskull;
	int counter = 0;
	bool step = false;
	while (t!=-500 && counter < 6)
	  {	
	    if (step)
	      {
		t = pt2.next_point_under(t, .30);
		counter++;
		gaps.push_back(t);
	      }
	    else
	      {
		t = pt2.next_point_over(t, .40);
		if (sum == -50) sum = t;
	      }
	    step = !step;
	
	  }
	if (counter > 1 && (sum - inskull) < 10)
	  {
	    OUTER_SKULL_SEARCH_BEGIN = gaps[0];
	  }
	else OUTER_SKULL_SEARCH_BEGIN = inskull + .5;
      }

      //outer skull

      //starting from the skin
      pt1.set_lroi(OUTER_SKULL_SEARCH_BEGIN);
      double outskull2 = 0;
      outskull2 = pt1.last_point_over(outskin, .75);
      outskull2 = pt1.last_point_under(outskull2, .30); 
      if (outskull2 == -500) outskull2 = OUTER_SKULL_SEARCH_BEGIN;

      //leftmost min
      double minabs = pt1.next_point_under(pt1.begin(), .30);
      if (minabs == -500) {output = false;exit (-1);}
	
      double localminabs = minabs; 
      if (output)
	{
	  bool stop = false;
	  const double lowthreshold = pt1.threshold(.15);

	  const double upthreshold = pt1.threshold(.25);
	  for (vector<pro_pair>::const_iterator i = pt1.v.begin(); i != pt1.v.end(); i++)
	    {
	      if (!stop && (*i).abs>=minabs && i!=pt1.v.end() && i!=pt1.v.begin())
		{
		  if ((*i).val>upthreshold) stop = true; //avoid climbing skin
		  if ((*i).val<lowthreshold) 
		    {
		      localminabs = (*i).abs;
		    }
		}     
	    }
	}
      double outskull = pt1.next_point_over(localminabs, .15);
      
      if (outskull2 - outskull >= 3 | outskull - outskull2 >=2) {output = false;}

      if (outskin - outskull < 2) output = false; 
      if (outskull - inskull < 0) output = false;

      if (output) 
	{
	  resul.push_back(outskin);
	  resul.push_back(outskull);
	  resul.push_back(inskull);
	}
      else resul.push_back(outskin);

      //computing csf
      pt1.init_roi();
      pt1.set_rroi(outskull - .5);
      pt2.init_roi();
      pt2.set_rroi(outskin);
      double end_csf = pt2.last_point_over(inskull, .85);
      if (end_csf != -500)
	{
	  for (double d = -30; d <= end_csf; d+=.5)
	    if (pt2.value(d) > pt2.threshold(.85) & pt1.value(d) < pt1.threshold(.45))
	      {
		Pt csf = point + d*n;
		csfvolume((int)floor(csf.X + .5), (int)floor(csf.Y + .5), (int)floor(csf.Z + .5))=5;
	      }
	}
    }
  
  return resul;
}



bool special_case(const Pt & point, const Vec & n, const trMatrix & M)
{
  bool result = false;
  bool test2 = false;
  for (int i = 0; i < 11; i+=5)
    {
      Pt p = point + i * n;
      double realx = (M.m11 * p.X) + (M.m12 * p.Y) + (M.m13 * p.Z) + M.m14;
      double realy = (M.m21 * p.X) + (M.m22 * p.Y) + (M.m23 * p.Z) + M.m24;
      double realz = (M.m31 * p.X) + (M.m32 * p.Y) + (M.m33 * p.Z) + M.m34;
      realx -= 88;
      realy -= 128;
      realz -= 74;

      double plan = realy * .0117213 + realz * -.021433355 - .6; 
      if (plan > 0) result = true;
      if (realy > 62 && realx > -16 && realx < 16) test2 = true;
    }
  if (test2) result = false;
  return result;
}

//writes externall skull computed from image on output.
void write_ext_skull(volume<float> & output_inskull, volume<float> & output_outskull, volume<float> & output_outskin, const volume<float> & t1, const volume<float> & t2, const Mesh & m, const trMatrix & M, volume<short> & csfvolume) {
  int glob_counter = 0;
  int rem_counter = 0;

  const double xdim = t1.xdim();
  const double ydim = t1.ydim();
  const double zdim = t1.zdim();

  double imax = t1.max();
  if (imax == 0) imax = 1;

  volume<short> meshimage;
  copyconvert(t1, meshimage);
  meshimage = 0;
  draw_mesh(meshimage, m);

  for (vector<Mpoint*>::const_iterator i = m._points.begin(); i != m._points.end(); i++)
    {
      (*i)->data.clear();

      double max_neighbour = 0;
      const Vec normal = (*i)->local_normal();
      const Vec n = Vec(normal.X/xdim, normal.Y/ydim, normal.Z/zdim);      
      
      for (list<Mpoint*>::const_iterator nei = (*i)->_neighbours.begin(); nei != (*i)->_neighbours.end(); nei++)
	max_neighbour = Max(((**i) - (**nei)).norm(), max_neighbour); 
      
      max_neighbour = ceil((max_neighbour)/2);
      
      const Pt mpoint((*i)->get_coord().X/xdim,(*i)->get_coord().Y/ydim,(*i)->get_coord().Z/zdim);
      for (int ck = (int)floor(mpoint.Z - max_neighbour/zdim); ck <= (int)floor(mpoint.Z + max_neighbour/zdim); ck++)
	for (int cj = (int)floor(mpoint.Y - max_neighbour/ydim); cj <= (int)floor(mpoint.Y + max_neighbour/ydim); cj++)
	  for (int ci = (int)floor(mpoint.X - max_neighbour/xdim); ci <= (int)floor(mpoint.X + max_neighbour/xdim); ci++)
	    {
	      bool compute = false;
	      const Pt point(ci, cj, ck);
	      const Pt realpoint(ci*xdim, cj*ydim, ck*zdim);
	      if (meshimage(ci, cj, ck) == 1) 
		{
		  double mindist = 10000;
		  for (list<Mpoint*>::const_iterator nei = (*i)->_neighbours.begin(); nei != (*i)->_neighbours.end(); nei++)
		    mindist = Min(((realpoint) - (**nei)).norm(), mindist); 
		  if (mindist >= ((realpoint) - (**i)).norm()) compute = true;
		}
	    

	      if (compute)
		{
		  glob_counter ++;
		  vector<double> val;
		  if (!special_case(realpoint, normal, M))
		    val = co_ext(t1, t2, point, n, csfvolume);
		  else
		    val = special_extract(t1, t2, point, n, csfvolume);

		  if (val.size() == 3)
		    {
		      Pt opoint(point.X, point.Y, point.Z); 
		      Vec on(n.X, n.Y, n.Z); 
		      Pt c0 = opoint + val[0]*on;
		      Pt c1 = opoint + val[1]*on;
		      Pt c2 = opoint + val[2]*on;

		      output_outskin((int)floor(c0.X + .5) + infxm,(int) floor(c0.Y + .5) + infym,(int) floor(c0.Z + .5) + infzm) +=1;
		      output_outskull((int)floor(c1.X + .5) + infxm,(int) floor(c1.Y + .5) + infym,(int) floor(c1.Z + .5) + infzm)+=1;
		      output_inskull((int)floor(c2.X + .5) + infxm,(int) floor(c2.Y + .5) + infym,(int) floor(c2.Z + .5) + infzm)+=1;
		    }
		  else {
		    rem_counter++;

		    if (val.size()==1)
		      {
			Pt opoint(point.X, point.Y, point.Z); 
			Vec on(n.X, n.Y, n.Z); 
			Pt c0 = opoint + val[0]*on;
			
			output_outskin((int)floor(c0.X + .5) + infxm ,(int) floor(c0.Y + .5) + infym ,(int) floor(c0.Z + .5) + infzm) +=1; 
		      }
		  }
		}
	    }
    }

  if (verbose.value())
    {
      cout<<" nb of profiles : "<<glob_counter<<endl;
      cout<<" removed profiles : "<<100. * rem_counter/(double) glob_counter<<"%"<<endl;
    }
}

int main(int argc, char *argv[]) {

  //parsing options
  OptionParser options(title, examples);
  options.add(help);
  options.add(verbose);
  options.add(t1only);
  options.add(outline);
  options.add(mask);
  options.add(skullmask);
  //  options.add(csf);
  options.add(precision);

  int nb_pars = 0;

  try {
    nb_pars = options.parse_command_line(argc, argv, 0, true);
  }
  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } 
  catch(std::exception &e) {
    cerr << e.what() << endl;
  } 
  
  if (help.value()) {options.usage(); return 0;};

  for (int i = nb_pars; i < argc; i++)
    {
      string s = argv[i];
    if (s.find("-")==0)
      {
	options.usage(); return -1;
      }
    }

  if (t1only.value()) {
    t1only_main(argc, argv, nb_pars, options);
    exit (0);
  }

  if (argc - nb_pars < 5)
    {
      cerr<<"too few arguments"<<endl;
      options.usage(); return -1;
    }

  int count_arg = nb_pars;

  const string inputt1(argv[count_arg]);
  count_arg++;
  const string inputt2(argv[count_arg]);
  count_arg++;
  const string mesh(argv[count_arg]);
  count_arg++;
  const string matrix(argv[count_arg]);
  count_arg++;
  const string outputstr(argv[count_arg]);

  l = 0;

  //load the mesh
  Mesh m;
  m.load(mesh);

  //load the matrix
  trMatrix M;
  ifstream f(matrix.c_str());
  if (f.is_open())
    {
      f>>M.m11>>M.m12>>M.m13>>M.m14>>M.m21>>M.m22>>M.m23>>M.m24>>M.m31>>M.m32>>M.m33>>M.m34>>M.m41>>M.m42>>M.m43>>M.m44;
      f.close();
    }
  else {cerr<<"unable to open "<<matrix<<endl;return -1;}

  //load the volumes
  volume<float> t1;
  volume<float> t2;

  if (read_volume(t1,inputt1.c_str())<0)  return -1;
  if (read_volume(t2,inputt2.c_str())<0)  return -1;

  t1.setinterpolationmethod(trilinear);
  t2.setinterpolationmethod(trilinear);

  const double xdim = t1.xdim();
  const double ydim = t1.ydim();
  const double zdim = t1.zdim();
  const int xsize = t1.xsize();
  const int ysize = t1.ysize();
  const int zsize = t1.zsize();

  //founding brain robustmax (useful if skull is too bright)
  double thr;
  
  {
    volume<short> brain_mask;
    copyconvert(t1, brain_mask);
    brain_mask = 0;
    draw_mesh(brain_mask, m);
    brain_mask = make_mask_from_mesh(brain_mask, m);
    vector<double> brain_hist;
    for (int k = 0; k < zsize; k++)
      for(int j = 0; j < ysize; j++)
	for(int i = 0; i < xsize; i++)
	  if (brain_mask.value(i, j, k) == 1) brain_hist.push_back(t1.value(i, j, k));
    
    int size = brain_hist.size();
    int e98 = (int) ceil(.98 * size);
    nth_element(brain_hist.begin(), brain_hist.begin() + e98 - 1, brain_hist.end());
    thr = brain_hist[e98 - 1];
    thr *= 1.1;
  }

  for (int k = 0; k < zsize; k++)
    for(int j = 0; j < ysize; j++)
      for(int i = 0; i < xsize; i++)
	if (t1.value(i, j, k) > thr) t1.value(i, j, k) = thr;
  
  infxm=0; infym=0; infzm=0; 
  int infxp=0, infyp=0, infzp=0;
  //checking if the mesh is inside the volume. If not, increasing the size of the volume.
  for (vector<Mpoint *>::const_iterator p = m._points.begin(); p != m._points.end(); p++)
    {
      if ((*p)->get_coord().X/xdim < 3) infxm = Max(infxm, (int) ceil(3 - (*p)->get_coord().X/xdim));
      if ((*p)->get_coord().Y/ydim < 3) infym = Max(infym, (int) ceil(3 - (*p)->get_coord().Y/ydim));
      if ((*p)->get_coord().Z/zdim < 3) infzm = Max(infzm, (int) ceil(3 - (*p)->get_coord().Z/zdim));
      if (xsize - (*p)->get_coord().X/xdim< 3) infxp = Max(infxp, (int) ceil(3 - xsize + (*p)->get_coord().X/xdim));
      if (ysize - (*p)->get_coord().Y/ydim< 3) infyp = Max(infyp, (int) ceil(3 - ysize + (*p)->get_coord().Y/ydim));
      if (zsize - (*p)->get_coord().Z/zdim < 3) infzp = Max(infzp, (int) ceil(3 - zsize + (*p)->get_coord().Z/zdim));
    }

  infxp += 2;
  infxm += 2;
  infyp += 2;
  infym += 2;
  infzp += 2;
  infzm += 2;

  volume<float> write0(xsize + infxm + infxp, ysize + infym + infyp, zsize + infzm + infzp);
  write0.setxdim(xdim);
  write0.setydim(ydim);
  write0.setzdim(zdim);

  write0 = 0;

  volume<float> write2;
  volume<float> write1;
  write1 = write2 = write0;

  volume<short> csfvolume;
  copyconvert(t1, csfvolume);
  csfvolume = 0;
  
  if (verbose.value()) cout<<"extracting profiles"<<endl;
  
  write_ext_skull(write0, write1, write2, t1, t2, m, M, csfvolume);


  const double rmin=3.33;
  const double rmax=10;
  const double E = (1/rmin + 1/rmax)/2.;
  const double F = 6./(1/rmin - 1/rmax);

  if (verbose.value()) cout<<"blurring"<<endl;

  const double blurring0 = 3;
  const double blurring1 = 3;
  const double blurring2 = 3;

  const double smoothness0 = 1;
  const double smoothness1 = 1;
  const double smoothness2 = 1;

  write0 = smooth(write0, blurring0);
  write2 = smooth(write2, blurring2);

  const int nb_iter0 = 800;
  const int nb_iter1 = 800;
  const int nb_iter2 = 1500;

  if (verbose.value()) cout<<"computation"<<endl;
  //N.B. you could remove the pop_back size check in standard_step, by placing this line _after_ mprecise=m
  for (vector<Mpoint*>::iterator i = m._points.begin(); i!= m._points.end(); i++)
    {
      (*i)->data.push_back(0);
    }

  //inner skull

  m.translation(xdim * infxm, ydim * infym, zdim * infzm);

  if (verbose.value()) cout<<" inner skull"<<endl;
  for (int c = 0; c < nb_iter0; c++)
    standard_step_of_computation(write0, m, c, E, F, smoothness0, .5, c);

  Mesh mprecise;
  if (precision.value() > 0)
    {
      mprecise = m;
      if (verbose.value()) cout<<"  increased precision"<<endl;
      for (int i = 0; i < precision.value(); i++)
	mprecise.retessellate(); 
      for (int c = 0; c < 100; c++)
	standard_step_of_computation(write0, mprecise, c, E, F, smoothness0, .5, c);
    }
  
  Mesh realmesh;
    if (precision.value() == 0)
      realmesh = m;
    else realmesh = mprecise;

  realmesh.translation(- xdim * infxm, - ydim * infym, - zdim * infzm);
  realmesh.save((outputstr+"_inskull_mesh.off").c_str());

  volume<short> maskinskull;
  volume<short> maskinskull2;
  volume<short> meshinskull;

  {
    copyconvert(write0, maskinskull);
    draw_mesh(maskinskull, m);
    maskinskull = make_mask_from_mesh(maskinskull, m);

    volume<short> output1;
    copyconvert(t1, output1);
    output1 = 0;

    if (outline.value() | mask.value() | skullmask.value())
      {
	draw_mesh(output1, realmesh);
      }
    
    if (skullmask.value())
      meshinskull = output1;

    if (outline.value())
      {
	if (verbose.value()) cout<<"  outline"<<endl;
	if (save_volume(output1, (outputstr+"_inskull_mesh").c_str())<0)  return -1;
      }

    if (mask.value() | skullmask.value())
      {
	volume<short> smask = make_mask_from_mesh(output1, realmesh);

	if (skullmask.value())
	  {
	    maskinskull2 = smask;
	  }
	if (mask.value())
	  {
	    if (verbose.value()) cout<<"  mask"<<endl;	
	    if (save_volume(smask, (outputstr+"_inskull_mask").c_str())<0)  return -1;
	  }

      }
  }


  {
    int xsiz = write1.xsize();
    int ysiz = write1.ysize();
    int zsiz = write1.zsize();
    for (int k = 0; k < zsiz; k++)
      for (int j = 0; j < ysiz; j++)
	for (int i = 0; i < xsiz; i++)
	  if (write1.value(i, j, k) > 0 && maskinskull.value(i, j, k) == 1)
	    {
	      write1.value(i, j, k) = 0;
	    }
  }

  write1 = smooth(write1, blurring1);

  if (verbose.value()) cout<<" outer skull"<<endl;
  for (int c = 0; c < nb_iter1; c++)
    standard_step_of_computation(write1, m, c, E, F, smoothness1, .5, c, 5, 15, true, maskinskull);

  maskinskull.destroy();

  if (precision.value() > 0)
    {
      if (verbose.value()) cout<<"  increased precision"<<endl;
      mprecise = m;
      for (int i = 0; i < precision.value(); i++)
	mprecise.retessellate();
      for (int c = 0; c < 100; c++)
	standard_step_of_computation(write1, mprecise, c, E, F, smoothness1, .5, c, 5, 15);
      
    }

  write1.destroy();
  
  volume<short> maskoutskull;
  volume<short> maskoutskull2;

  Mesh realmesh2;
  if (precision.value() == 0)
    realmesh2 = m;
  else realmesh2 = mprecise;
  
  realmesh2.translation(- xdim * infxm, - ydim * infym, - zdim * infzm);
  realmesh2.save((outputstr+"_outskull_mesh.off").c_str());


  {
    copyconvert(write0, maskoutskull);
    write0.destroy();
    draw_mesh(maskoutskull, m);
    maskoutskull = make_mask_from_mesh(maskoutskull, m);
    
    volume<short> output1;
    copyconvert(t1, output1);
    output1 = 0;

    if (mask.value() | outline.value() | skullmask.value())
      draw_mesh(output1, realmesh2);
    
    if (outline.value())
      {
	if (verbose.value()) cout<<"  outline"<<endl;
	if (save_volume(output1, (outputstr+"_outskull_mesh").c_str())<0)  return -1;
      }
    
    if (mask.value()|skullmask.value())
      {
	volume<short> smask = make_mask_from_mesh(output1, realmesh2);
	if (skullmask.value())
	  {
	    maskoutskull2 = smask;
	  }
	if (mask.value())
	  {
	    if (verbose.value()) cout<<"  mask"<<endl;
	    if (save_volume(smask, (outputstr+"_outskull_mask").c_str())<0)  return -1;
	  }
      }
    
  }

  
  if (skullmask.value())
    {
      volume<short> smask = maskoutskull2;
      smask = 0;
      int xsize = t1.xsize(), ysize = t1.ysize(), zsize = t1.zsize();
      if (verbose.value()) cout<<" skull mask"<<endl;
      for (int k = 0; k < zsize; k++)
	for(int j = 0; j < ysize; j++)
	  for(int i = 0; i < xsize; i++)
	    smask.value(i, j, k) = Min (1, maskoutskull2.value(i, j, k) * (1 - maskinskull2.value(i, j, k)) + meshinskull.value(i, j, k));
      
      if (save_volume(smask, (outputstr+"_skull_mask").c_str())<0) return -1;
    }

  maskoutskull2.destroy();
  maskinskull2.destroy();
  meshinskull.destroy();
  
  //computing the starting mesh for outer skull
  if (verbose.value()) cout<<" skin"<<endl;
  
  Mesh m2;
  make_mesh_from_icosa(5, m2);

  for (vector<Mpoint*>::iterator i = m2._points.begin(); i!= m2._points.end(); i++)
    {
      (*i)->data.push_back(0);
    }
  
  Pt p;
  int counter=0;
  for (vector<Mpoint*>::const_iterator i = m._points.begin(); i!= m._points.end(); i++)
    {
      counter++;
      p+=(*i)->get_coord();
    }
  p*=(1./counter);

  double radius = 0;
  for (vector<Mpoint*>::const_iterator i = m._points.begin(); i!= m._points.end(); i++)
    radius+=((*i)->get_coord()-p).norm();

  radius/=(counter);
  radius*=.75;

  m2.rescale(radius);
  m2.translation(p.X, p.Y, p.Z);

  for (int c = 0; c < nb_iter2; c++)
    standard_step_of_computation(write2, m2, c, E, F, smoothness2, 1.5, c, 5, 15, true, maskoutskull);

  maskoutskull.destroy();

  if (precision.value() > 0)
    {
      if (verbose.value()) cout<<"  increased precision"<<endl;
      for (int i = 0; i < precision.value(); i++)
	m2.retessellate();
      for (int c = 0; c < 100; c++)
	standard_step_of_computation(write2, m2, c, E, F, smoothness2, 1.5, c, 5, 15);
      
    }

  write2.destroy();
  
  m2.translation(-xdim * infxm, -ydim * infym, -zdim * infzm);
  m2.save((outputstr+"_outskin_mesh.off").c_str());
  
  if (outline.value()| mask.value())
    {
      volume<short> output1;
      copyconvert(t1, output1);
      output1 = 0;

      draw_mesh(output1, m2);
      
      if (outline.value())
	{
	  if (verbose.value()) cout<<"  outline"<<endl;
	  if (save_volume(output1, (outputstr+"_outskin_mesh").c_str())<0)  return -1;  
	}
      if (mask.value())
	{
	  if (verbose.value()) cout<<"  mask"<<endl;
	  volume<short> mask = make_mask_from_mesh(output1, m2);
	  if (save_volume(mask, (outputstr+"_outskin_mask").c_str())<0)  return -1;
	}
    }
  
  return (0);
}

