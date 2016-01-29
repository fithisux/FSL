//     fslcc.cc cross-correlate two time series timepoint by timepoint
//     Steve Smith and Matthew Webster, FMRIB Image Analysis Group
//     Copyright (C) 2001-2013 University of Oxford  
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
#include "utils/options.h"
#include <iomanip>

using namespace NEWIMAGE;
using namespace Utilities;

Option<string> fnmask(string("-m"), string(""),
		string("mask file name "),
		false, requires_argument);
Option<bool> noabs(string("--noabs"), false, 
		     string("\tDon't return absolute values (keep sign)"), 
		     false, no_argument);
Option<bool> nodemean(string("--nodemean"), false, 
		     string("Don't demean the input files"), 
		     false, no_argument);
Option<float> thresh(string("-t"), 0.1,
		     string("\tThreshhold ( default 0.1 )"),
		     false, requires_argument);
Option<float> precision(string("-p"), 2,
		     string("\tNumber of decimal places to display in output ( default 2 )"),
		     false, requires_argument);

template <class T>
int fmrib_main(int argc, char *argv[])
{
  volume4D<T> input_volume1, input_volume2;
  string input_name1(argv[0]);
  string input_name2(argv[1]);
  read_volume4D(input_volume1,input_name1);
  read_volume4D(input_volume2,input_name2);

  if (input_volume1.maxx() != input_volume2.maxx() ||  input_volume1.maxy() != input_volume2.maxy()  ||  input_volume1.maxz() != input_volume2.maxz())
  {
    cerr << "Error: Mismatch in image dimensions" << endl; 
    return 1;
  }

  volume<T> mask;
  if(fnmask.value().length()>0){
	read_volume(mask,fnmask.value());
	if(!samesize(input_volume1[0],mask)){  	
	  cerr << "Error: Mismatch in mask dimensions" << endl; 
	  return 1;
	}
  } else {
    mask=input_volume1[0];
    mask=1;
  }
 

   if ( !nodemean.value() ) {
    for(int t1=0;t1<=input_volume1.maxt();t1++)
      input_volume1[t1] -= input_volume1[t1].mean(mask);
    for(int t2=0;t2<=input_volume2.maxt();t2++)
      input_volume2[t2] -= input_volume2[t2].mean(mask);
   }

  for(int t1=0;t1<=input_volume1.maxt();t1++)
  {
    double ss1=sqrt(input_volume1[t1].sumsquares(mask));  
    for(int t2=0;t2<=input_volume2.maxt();t2++)
    {
       double ss2=sqrt(input_volume2[t2].sumsquares(mask));  
       double score=0;
       for(int k=0;k<=input_volume1.maxz();k++)
         for(int j=0;j<=input_volume1.maxy();j++)
           for(int i=0;i<=input_volume1.maxx();i++)
	     if ( !fnmask.set() || mask(i,j,k)>0 ) 
	       score+=(double)input_volume1(i,j,k,t1)*(double)input_volume2(i,j,k,t2); 
       if (!noabs.value())
	 score=fabs(score);
       score/=(ss1*ss2);
       if (score>thresh.value())
         cout << setw(3) << t1+1 << " " << setw(3) << t2+1 << " " <<  setiosflags (ios::fixed) << setprecision(precision.value()) << score << endl;
    }
  }

  return 0;
}


int main(int argc,char *argv[])
{
  string title("fslcc: Cross-correlate two time-series, timepoint by timepoint");
  string examples("fslcc [options] <first_input> <second_input> ");
  OptionParser options(title, examples);

  options.add(fnmask);
  options.add(noabs);
  options.add(nodemean);
  options.add(thresh);
  options.add(precision);
  unsigned int done;

  try {
    done=options.parse_command_line(argc, argv,0,true);
  }
  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    return(1);
  } 
  int extraArgs=argc-done;
  argv+=done;

  if ( extraArgs != 2 )
  {
    options.usage();
    return(1);
  }
  return call_fmrib_main(dtype(string(argv[0])),argc,argv); 
}
