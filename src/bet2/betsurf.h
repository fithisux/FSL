
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

#ifndef _t1only
#define _t1only

#include "meshclass/meshclass.h"
#include "newimage/newimageall.h"
#include "utils/options.h"

using namespace mesh;
using namespace NEWIMAGE;
using namespace Utilities;
using namespace std;

struct trMatrix
{
  double m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34, m41, m42, m43, m44;
};

volume<short> default_volume;

void draw_segment(volume<short>& image, const Pt& p1, const Pt& p2);

void draw_mesh(volume<short>& image, const Mesh &m);

volume<short> make_mask_from_mesh(const volume<short> & image, const Mesh& m);

double standard_step_of_computation(const volume<float> & image, Mesh & m, const int iteration_number, const double E,const double F, const float addsmooth, const float speed, const int nb_iter, const int id=5, const int od=15, const bool vol=false, const volume<short> & mask=default_volume);

vector<double> t1only_special_extract(const volume<float> & t1, const Pt & point, const Vec & n) ;

vector<double> t1only_co_ext(const volume<float> & t1, const Pt & point, const Vec & n) ;

void t1only_write_ext_skull(volume<float> & output_inskull, volume<float> & output_outskull, volume<float> & output_outskin, const volume<float> & t1, const Mesh & m, const trMatrix & M) ;

int t1only_main(int argc, char *argv[], int nb_pars, OptionParser & options);

vector<double> special_extract(const volume<float> & t1, const volume<float> & t2, const Pt & point, const Vec & n, volume<short> & csfvolume);

vector<double> co_ext(const volume<float> & t1, const volume<float> & t2, const Pt & point, const Vec & n, volume<short> & csfvolume);

bool special_case(const Pt & point, const Vec & n, const trMatrix & M);

void write_ext_skull(volume<float> & output_inskull, volume<float> & output_outskull, volume<float> & output_outskin, const volume<float> & t1, const volume<float> & t2, const Mesh & m, const trMatrix & M, volume<short> & csfvolume);


#endif












