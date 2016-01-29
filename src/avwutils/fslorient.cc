/*  fslorient.cc

    Mark Jenkinson and Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 2003-20047 University of Oxford  */

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


#ifndef EXPOSE_TREACHEROUS
#define EXPOSE_TREACHEROUS
#endif

#include "newimage/newimageall.h"
#include "newimage/fmribmain.h"

using namespace NEWIMAGE;

void print_usage() {
  string progname="fslorient";
  cout << "Usage: " << progname << " <main option> <filename>" << endl;
  cout << endl;
  cout << "  where the main option is one of:" << endl;
  cout << "    -getorient             (prints FSL left-right orientation)" << endl;
  cout << "    -getsform              (prints the 16 elements of the sform matrix)" << endl;
  cout << "    -getqform              (prints the 16 elements of the qform matrix)" << endl;
  cout << "    -setsform <m11 m12 ... m44>  (sets the 16 elements of the sform matrix)" << endl;
  cout << "    -setqform <m11 m12 ... m44>  (sets the 16 elements of the qform matrix)" << endl;
  cout << "    -getsformcode          (prints the sform integer code)" << endl;
  cout << "    -getqformcode          (prints the qform integer code)" << endl;
  cout << "    -setsformcode <code>   (sets sform integer code)" << endl;
  cout << "    -setqformcode <code>   (sets qform integer code)" << endl;
  cout << "    -copysform2qform       (sets the qform equal to the sform - code and matrix)" << endl;
  cout << "    -copyqform2sform       (sets the sform equal to the qform - code and matrix)" << endl;
  cout << "    -deleteorient          (removes orient info from header)" << endl;
  cout << "    -forceradiological     (makes FSL radiological header)" << endl;
  cout << "    -forceneurological     (makes FSL neurological header - not Analyze)" << endl;
  cout << "    -swaporient            (swaps FSL radiological and FSL neurological)" << endl;
  cout << endl;
  cout << "       Note: the stored data order is never changed here - only the header info." << endl;
  cout << "       To change the data storage use fslswapdim." << endl;
  cout << endl;
  cout << "  e.g.  " << progname << " -forceradiological myimage" << endl;
  cout << "        " << progname << " -copysform2qform myimage" << endl;
  cout << "        " << progname << " -setsform -2 0 0 90 0 2 0 -126 0 0 2 -72 0 0 0 1 myimage" << endl;
  cout << endl;
}



template <class T>
void swaporient(volume4D<T>& invol) {
  Matrix swapmat(IdentityMatrix(4));
  swapmat(1,1)=-1;
  swapmat(1,4)=invol.xsize()-1;
  invol.set_sform(invol.sform_code(),invol.sform_mat() * swapmat);
  invol.set_qform(invol.qform_code(),invol.qform_mat() * swapmat);
}


template <class T>
int getorient(const volume4D<T>& invol) 
{
  if (invol.left_right_order()==FSL_RADIOLOGICAL) { 
    cout << "RADIOLOGICAL" << endl;
  } else {
    cout << "NEUROLOGICAL" << endl;
  }
  return 0;
}


template <class T>
int getorient(const string& filename) 
{
  volume4D<T> invol;
  read_orig_volume4D(invol,filename);
  return getorient(invol);
}


int showmat(const Matrix& mat) {
  for (int a=1; a<=4; a++) {
    for (int b=1; b<=4; b++) {
      cout << mat(a,b) << " ";
    }
  }
  cout << endl;
  return 0;
}


int testminargs(int argnum, int argc) {
  if (argc < (argnum + 1)) {
    print_usage();
    exit(EXIT_FAILURE);
  }
  return 0;
}


template <class T>
int fmrib_main(int argc,char *argv[])
{
  bool modified=false;
  int retval=0;
  string option=argv[1], filename;

  volume4D<T> invol;

  filename=argv[argc-1];
  read_orig_volume4D(invol,filename);

  if (option=="-getorient") {
    getorient(invol);
  } else if (option=="-getsformcode") {
    cout << invol.sform_code() << endl;
  } else if (option=="-getqformcode") {
    cout << invol.qform_code() << endl;
  } else if (option=="-getqform") {
    showmat(invol.qform_mat());
  } else if (option=="-getsform") {
    showmat(invol.sform_mat());
  } else if (option=="-setsformcode") {
    testminargs(3,argc);
    modified=true;
    int code = (int) atof(argv[2]);
    invol.set_sform(code,invol.sform_mat());
  } else if (option=="-setqformcode") {
    testminargs(3,argc);
    modified=true;
    int code = (int) atof(argv[2]);
    invol.set_qform(code,invol.qform_mat());
  } else if ( (option=="-setqform") || (option=="-setsform") ) {
    testminargs(18,argc);
    modified=true;
    Matrix mat(4,4);
    for (int a=1; a<=4; a++) {
      for (int b=1; b<=4; b++) {
	mat(a,b)=atof(argv[a*4+b-3]);
      }
    }
    // override the last line
    mat(4,1)=0; mat(4,2)=0; mat(4,3)=0; mat(4,4)=1;
    if (option=="-setqform") {
      invol.set_qform(invol.qform_code(),mat);
    }
    if (option=="-setsform") {
      invol.set_sform(invol.sform_code(),mat);
    }
  } else if (option=="-copysform2qform") {
    modified=true;
    invol.set_qform(invol.sform_code(),invol.sform_mat());
  } else if (option=="-copyqform2sform") {
    modified=true;
    invol.set_sform(invol.qform_code(),invol.qform_mat());
  } else if (option=="-swaporient") {
    modified=true;
    swaporient(invol);
  } else if (option=="-deleteorient") {
    modified=true;
    invol.set_sform(NIFTI_XFORM_UNKNOWN,invol.sform_mat());
    invol.set_qform(NIFTI_XFORM_UNKNOWN,invol.qform_mat());
  } else if (option=="-forceradiological") {
    modified=true;
    if (invol.left_right_order()==FSL_NEUROLOGICAL) {
      swaporient(invol);
    }
  } else if (option=="-forceneurological") {
    modified=true;
    if (invol.left_right_order()==FSL_RADIOLOGICAL) {
      swaporient(invol);
    }
  } else{
    // does the first arg exist as an image file?
    if (fsl_imageexists(string(argv[1]))) {
      getorient<T>(string(argv[1]));
    } else {
      cerr << "Unrecognised option: " << option << endl;
      print_usage();
      retval = -1;
    }
  }


  if (modified) {
    if (FslBaseFileType(fslFileType(filename))!=FSL_TYPE_ANALYZE) {
      FslSetOverrideOutputType(fslFileType(filename));
      save_orig_volume4D(invol,filename);
      FslSetOverrideOutputType(-1);  // restore to default
    } else {
      cerr << "Cannot modify orientation for Analyze files" << endl;
      cerr << "  All Analyze files are treated as radiological" << endl;
      cerr << "  To change the data storage use fslswapdim" << endl;
      retval=-1;
    }
  }
  return retval;
}

  
int main(int argc,char *argv[])
{
  if (argc<2) { 
    print_usage();
    return -1; 
  }
  
  string inname;
  inname = argv[argc-1];
  // call the templated main
  return call_fmrib_main(dtype(inname),argc,argv);
}

