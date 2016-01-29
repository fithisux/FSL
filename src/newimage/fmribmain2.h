/*  General call feature for templated image classes

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2000 University of Oxford  */

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

#if !defined(__fmribmain2_h)
#define __fmribmain2_h

#include <iostream>

template <class T, class S>
int fmrib_main(int argc, char* argv[]);

template <class T>
int call_fmrib_main(short datatype, int argc, char* argv[])
{
  datatype=NEWIMAGE::closestTemplatedType(datatype);
  if ( datatype==DT_UNSIGNED_CHAR ) return fmrib_main<T,char>(argc, argv);
  else if ( datatype==DT_SIGNED_SHORT ) return fmrib_main<T,short>(argc, argv);
  else if ( datatype==DT_SIGNED_INT ) return fmrib_main<T,int>(argc, argv);
  else if ( datatype==DT_FLOAT )  return fmrib_main<T,float>(argc, argv);
  else if ( datatype==DT_DOUBLE ) return fmrib_main<T,double>(argc, argv);
  return -1;
}

int call_fmrib_main(short datatype1, short datatype2, int argc, char* argv[])
{
  datatype1=NEWIMAGE::closestTemplatedType(datatype1);
  if ( datatype1==DT_UNSIGNED_CHAR ) return call_fmrib_main<char>(datatype2, argc, argv);
  else if ( datatype1==DT_SIGNED_SHORT ) return call_fmrib_main<short>(datatype2, argc, argv);
  else if ( datatype1==DT_SIGNED_INT ) return call_fmrib_main<int>(datatype2, argc, argv);
  else if ( datatype1==DT_FLOAT )  return call_fmrib_main<float>(datatype2, argc, argv);
  else if ( datatype1==DT_DOUBLE ) return call_fmrib_main<double>(datatype2, argc, argv);
  return -1;
}

#endif



