/*  applyxfm4D.cc

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

#include "newimage/newimageall.h"
#include "newimage/fmribmain.h"

using namespace NEWIMAGE;

// Globals - needed by fmrib_main

string oname, iname, transname, refname, matprefix="/MAT_0";
bool singlematrix, fourd;


//////////////////////////////////////////////////////////////////////

template <class T>
int fmrib_main(int argc, char* argv[])
{
  if (fourd) {
    // 4D mode
    volume4D<T> invol, outvol;
    volume<T> refvol, dummy;
    read_volume4D(invol,iname);
    for (int t=0; t<invol.tsize(); t++) {
      invol[t].setpadvalue(invol[t].backgroundval());
    }
    invol.setextrapolationmethod(extraslice);
    invol.setinterpolationmethod(sinc);
    invol.definesincinterpolation("b",7);

    // old form used a volume number
    //    refvol = invol[atoi(refname.c_str())];  

    read_volume(refvol,refname);

    Matrix affmat(4,4);
    string matname;
    if (singlematrix) { affmat = read_ascii_matrix(transname); }

    if (invol.maxt() - invol.mint() > 10000) {
      cerr << "WARNING:: More than 10000 volumes - only doing first 10000" << endl;
    }

    for (int m=invol.mint(); m<=Min(invol.maxt(),invol.mint()+10000); m++) {

      if (!singlematrix) {
	matname = transname + matprefix;
	char nc='0';
	int n = m;
	matname += (nc + (char) (n / 1000));
	n -= (n/1000)*1000;
	matname += (nc + (char) (n / 100));
	n -= (n/100)*100;
	matname += (nc + (char) (n / 10));
	n -= (n/10)*10;
	matname += (nc + (char) n);
	cout << matname << endl;
	affmat = read_ascii_matrix(matname);
      }
      
      dummy = refvol;
      affine_transform(invol[m],dummy,affmat);
      outvol.addvolume(dummy);
    }
    outvol.setDisplayMaximumMinimum(0,0);
    save_volume4D(outvol,oname);

  } else {
    // 3D mode
    volume<T> invol, outvol;
    
    read_volume(invol,iname);
    read_volume(outvol,refname);
    invol.setextrapolationmethod(extraslice);
    invol.setinterpolationmethod(sinc);
    invol.definesincinterpolation("r",9);

    Matrix affmat(4,4);
    affmat = read_ascii_matrix(transname);
    
    affine_transform(invol,outvol,affmat);
    outvol.setDisplayMaximumMinimum(0,0);
    save_volume(outvol,oname);
  }

  return 0;
}


int main(int argc,char *argv[])
{

  Tracer tr("main");
  if (argc<5) { 
    cerr << "Usage: " << argv[0] << " <input volume> <ref volume>"
	 << " <output volume> <transformation matrix file/[dir]> [-singlematrix/-fourdigit/-userprefix <prefix>]]\n"; 
    return -1; 
  }
  
  // NB: a hidden option (-3D) exists (must appear after singlematrix)


  // parse the command line
  oname = argv[3];
  iname = argv[1];
  transname = argv[4]; 
  refname = argv[2];
  singlematrix = false;
  if (argc>=6) {
    string option = argv[5];
    if (option == "-singlematrix" )  singlematrix = true;
    if (option == "-fourdigit" )  matprefix = "/MAT_";
    if ( (option == "-userprefix" ) && (argc>=7) ) {
      string uprefix = argv[6];
      matprefix = "/" + uprefix;
    }
  }
  fourd = true;
  if (argc>=7) {
    string option = argv[6];
    if (option == "-3D" )  fourd = false;
  }

  // call the templated main
  return call_fmrib_main(dtype(iname),argc,argv);

}



